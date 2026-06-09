@doc raw"""
    HundSpinFlipHST{T,R} <: AbstractAsymHST{T,R}

Exact three-state auxiliary-field decomposition for the transverse Hund spin-flip term
for a local two-orbital pair.  The decomposition uses spin-sector-preserving local orbital
bilinears ``O_F`` and ``O_G`` documented in `docs/hund_spin_flip_hs_decomposition.md`.

This HST assumes asymmetric exact propagators and a time-slice ordering

```text
B_l = exp(-Δτ V_l) W_F,l W_G,l exp(-Δτ K_l).
```
"""
struct HundSpinFlipHST{T,R} <: AbstractAsymHST{T,R}
    β::R
    Δτ::R
    Lτ::Int
    N::Int
    JH::Vector{R}
    λ::Vector{R}
    weight0::Vector{R}
    weight1::Vector{R}
    neighbor_table::Matrix{Int}
    orbital_pairs::Vector{Tuple{Int,Int}}
    sF::Array{Int,2}
    sG::Array{Int,2}
    update_perm::Vector{Int}
end

function hund_spinflip_hs_coefficients(JH::R, Δτ::R) where {R<:AbstractFloat}
    @assert JH >= 0 "HundSpinFlipHST currently requires JH >= 0."
    κ = Δτ * JH / 4
    if iszero(κ)
        return zero(R), one(R), zero(R)
    end
    y = (exp(4κ) - 1) / (2 * (exp(κ) - 1)) - 1
    λ = acosh(y)
    B = (exp(κ) - 1) / (cosh(λ) - 1)
    A = 1 - B
    @assert A > 0 && B > 0 "Non-positive spin-flip HS weights: A=$A, B=$B."
    return λ, A, B/2
end

function hund_spinflip_field_action(s::Int, weight0::R, weight1::R) where {R<:AbstractFloat}
    # `s` is an HS field value for absolute factors and a field difference for update factors.
    w = iszero(s) ? weight0 : weight1
    return -log(w)
end

function sample_new_hund_spinflip_field(rng::AbstractRNG, s::Int)
    # `s` is an HS field value for absolute factors and a field difference for update factors.
    if s == -1
        return rand(rng, (0, 1))
    elseif s == 0
        return rand(rng, (-1, 1))
    else
        return rand(rng, (-1, 0))
    end
end

function _local_orbital_generator(::Val{:F}, ::Type{T}) where {T<:Number}
    return T[0 1; 1 0]
end

function _local_orbital_generator(::Val{:G}, ::Type{T}) where {T<:Number}
    return T[0 -im; im 0]
end

function _hund_local_W(channel::Symbol, s::Int, λ::R, spin_sign::Int, ::Type{T}=Complex{R}) where {R<:AbstractFloat,T<:Number}
    @assert channel === :F || channel === :G
    @assert spin_sign == 1 || spin_sign == -1
    # `s` is an HS field value for absolute factors and a field difference for update factors.
    a = spin_sign * s * λ
    I2 = Matrix{T}(I, 2, 2)
    O = channel === :F ? _local_orbital_generator(Val(:F), T) : _local_orbital_generator(Val(:G), T)
    return cosh(a) * I2 + sinh(a) * O
end

_hund_inv_local_W(channel::Symbol, s::Int, λ, spin_sign::Int, ::Type{T}=Complex{typeof(λ)}) where {T<:Number} =
    _hund_local_W(channel, -s, λ, spin_sign, T)

function _hund_field_action(hst::HundSpinFlipHST, n::Int, s::Int)
    return hund_spinflip_field_action(s, hst.weight0[n], hst.weight1[n])
end

function _apply_local_left_factor_to_asym_exact!(B::AsymExactPropagator, Xw::AbstractMatrix, inds::AbstractVector{Int})
    X = Matrix(Xw)
    Xinv = inv(X)
    rows = B.expmΔτK[inds, :]
    B.expmΔτK[inds, :] .= X * rows
    cols = B.exppΔτK[:, inds]
    B.exppΔτK[:, inds] .= cols * Xinv
    return nothing
end

function _green_delta_from_orbital_factor(B::AsymExactPropagator, Xw::AbstractMatrix{T}, inds::AbstractVector{Int}) where {T<:Number}
    k = length(inds)
    @assert size(Xw) == (k, k)
    X = Matrix{promote_type(T, eltype(B.expmΔτV))}(undef, k, k)
    v = B.expmΔτV
    @inbounds for a in 1:k, b in 1:k
        X[a,b] = v[inds[a]] * Xw[a,b] / v[inds[b]]
    end
    Δ = X - Matrix{eltype(X)}(I, k, k)
    return Δ
end

@doc raw"""
    HundSpinFlipHST(; model_geometry, orbital_pairs=[(1,2)], JH, β, Δτ, rng)

Construct spin-flip Hund HS fields for each local orbital pair in every unit cell.
"""
function HundSpinFlipHST(;
    model_geometry::ModelGeometry{D,R},
    orbital_pairs::AbstractVector{<:Tuple{Int,Int}} = [(1,2)],
    JH::R,
    β::R, Δτ::R, rng::AbstractRNG,
) where {D,R<:AbstractFloat}
    (; lattice, unit_cell) = model_geometry
    N_unitcells = lattice.N
    n_pairs = length(orbital_pairs)
    N = N_unitcells * n_pairs

    JH_vec = fill(JH, N)
    λ = similar(JH_vec)
    weight0 = similar(JH_vec)
    weight1 = similar(JH_vec)
    for n in eachindex(JH_vec)
        λ[n], weight0[n], weight1[n] = hund_spinflip_hs_coefficients(JH_vec[n], Δτ)
    end

    neighbor_table = zeros(Int, 2, N)
    nt′ = reshape(neighbor_table, (2, N_unitcells, n_pairs))
    for (p, (orb_a, orb_b)) in enumerate(orbital_pairs)
        for u in 1:N_unitcells
            nt′[1,u,p] = loc_to_site(u, orb_a, unit_cell)
            nt′[2,u,p] = loc_to_site(u, orb_b, unit_cell)
        end
    end

    Lτ = round(Int, β / Δτ)
    sF = rand(rng, (-1,0,1), (N, Lτ))
    sG = rand(rng, (-1,0,1), (N, Lτ))
    update_perm = collect(1:N)

    return HundSpinFlipHST{Complex{R},R}(β, Δτ, Lτ, N, JH_vec, λ, weight0, weight1,
                                         neighbor_table, collect(Tuple{Int,Int}, orbital_pairs),
                                         sF, sG, update_perm)
end

function _initialize!(
    fermion_path_integral_up::FermionPathIntegral{H},
    fermion_path_integral_dn::FermionPathIntegral{H},
    hst_parameters::HundSpinFlipHST{T,R},
) where {H<:Number,T<:Number,R<:AbstractFloat}
    @assert !((H<:Real) && (T<:Complex)) "Green's function matrices are real while HundSpinFlipHST is complex. Use complex BHZ hoppings or force complex path integrals."
    @assert fermion_path_integral_up.Sb == fermion_path_integral_dn.Sb "$(fermion_path_integral_up.Sb) ≠ $(fermion_path_integral_dn.Sb)"

    (; neighbor_table, JH, sF, sG, Lτ) = hst_parameters
    Vup = fermion_path_integral_up.V
    Vdn = fermion_path_integral_dn.V

    for n in axes(neighbor_table, 2)
        i, j = neighbor_table[1,n], neighbor_table[2,n]
        shift = JH[n] / 2
        for l in 1:Lτ
            Vup[i,l] += shift
            Vup[j,l] += shift
            Vdn[i,l] += shift
            Vdn[j,l] += shift
        end
        for l in 1:Lτ
            fermion_path_integral_up.Sb += _hund_field_action(hst_parameters, n, sF[n,l])
            fermion_path_integral_up.Sb += _hund_field_action(hst_parameters, n, sG[n,l])
        end
    end
    fermion_path_integral_dn.Sb = fermion_path_integral_up.Sb

    return nothing
end

@doc raw"""
    apply_hund_spinflip_to_propagators!(Bup, Bdn, hst)

Insert the current spin-flip Hund one-body factors into asymmetric exact propagators.
Call this immediately after `initialize_propagators` and before constructing/calculating Green's functions.
"""
function apply_hund_spinflip_to_propagators!(
    Bup::AbstractVector{<:AsymExactPropagator},
    Bdn::AbstractVector{<:AsymExactPropagator},
    hst::HundSpinFlipHST,
)
    (; neighbor_table, sF, sG, λ, Lτ) = hst
    @assert length(Bup) == Lτ && length(Bdn) == Lτ
    for l in 1:Lτ
        for n in axes(neighbor_table, 2)
            inds = @view neighbor_table[:,n]
            WFup = _hund_local_W(:F, sF[n,l], λ[n], +1)
            WGup = _hund_local_W(:G, sG[n,l], λ[n], +1)
            WFdn = _hund_local_W(:F, sF[n,l], λ[n], -1)
            WGdn = _hund_local_W(:G, sG[n,l], λ[n], -1)
            _apply_local_left_factor_to_asym_exact!(Bup[l], WFup * WGup, inds)
            _apply_local_left_factor_to_asym_exact!(Bdn[l], WFdn * WGdn, inds)
        end
    end
    return nothing
end

function _attempt_hund_update!(
    channel::Symbol,
    Gup::Matrix{H}, logdetGup::R, sgndetGup::H,
    Gdn::Matrix{H}, logdetGdn::R, sgndetGdn::H,
    hst::HundSpinFlipHST{T,R},
    fpi_up::FermionPathIntegral{H}, fpi_dn::FermionPathIntegral{H},
    Bup::AsymExactPropagator, Bdn::AsymExactPropagator,
    n::Int, l::Int, rng::AbstractRNG,
) where {H<:Number,T<:Number,R<:Real}
    fields = channel === :F ? hst.sF : hst.sG
    sold = fields[n,l]
    snew = sample_new_hund_spinflip_field(rng, sold)
    ΔSb = _hund_field_action(hst, n, snew) - _hund_field_action(hst, n, sold)
    inds = @view hst.neighbor_table[:,n]

    if channel === :F
        Xw_up = _hund_local_W(:F, snew - sold, hst.λ[n], +1)
        Xw_dn = _hund_local_W(:F, snew - sold, hst.λ[n], -1)
    else
        WF_up = _hund_local_W(:F, hst.sF[n,l], hst.λ[n], +1)
        WF_dn = _hund_local_W(:F, hst.sF[n,l], hst.λ[n], -1)
        ΔWG_up = _hund_local_W(:G, snew - sold, hst.λ[n], +1)
        ΔWG_dn = _hund_local_W(:G, snew - sold, hst.λ[n], -1)
        Xw_up = WF_up * ΔWG_up * inv(WF_up)
        Xw_dn = WF_dn * ΔWG_dn * inv(WF_dn)
    end

    Δup = _green_delta_from_orbital_factor(Bup, Xw_up, inds)
    Δdn = _green_delta_from_orbital_factor(Bdn, Xw_dn, inds)
    Rup = local_matrix_update_det_ratio(Gup, Δup, inds)
    Rdn = local_matrix_update_det_ratio(Gdn, Δdn, inds)
    P = abs(exp(-ΔSb) * Rup * Rdn)

    accepted = false
    if rand(rng) < P
        accepted = true
        fields[n,l] = snew
        logdetGup, sgndetGup = local_matrix_update_greens!(Gup, logdetGup, sgndetGup, Rup, Δup, inds)
        logdetGdn, sgndetGdn = local_matrix_update_greens!(Gdn, logdetGdn, sgndetGdn, Rdn, Δdn, inds)
        _apply_local_left_factor_to_asym_exact!(Bup, Xw_up, inds)
        _apply_local_left_factor_to_asym_exact!(Bdn, Xw_dn, inds)
        fpi_up.Sb += ΔSb
        fpi_dn.Sb += ΔSb
    end

    return accepted, logdetGup, sgndetGup, logdetGdn, sgndetGdn
end

function _local_updates!(
    Gup::Matrix{H}, logdetGup::R, sgndetGup::H,
    Gdn::Matrix{H}, logdetGdn::R, sgndetGdn::H,
    hst::HundSpinFlipHST{T,R},
    fpi_up::FermionPathIntegral{H},
    fpi_dn::FermionPathIntegral{H},
    Bup::P, Bdn::P, l::Int, rng::AbstractRNG,
) where {H<:Number,T<:Number,R<:Real,P<:AsymExactPropagator}
    shuffle!(rng, hst.update_perm)
    accepted = 0
    for n in hst.update_perm
        acc, logdetGup, sgndetGup, logdetGdn, sgndetGdn = _attempt_hund_update!(
            :F, Gup, logdetGup, sgndetGup, Gdn, logdetGdn, sgndetGdn,
            hst, fpi_up, fpi_dn, Bup, Bdn, n, l, rng,
        )
        accepted += acc
        acc, logdetGup, sgndetGup, logdetGdn, sgndetGdn = _attempt_hund_update!(
            :G, Gup, logdetGup, sgndetGup, Gdn, logdetGdn, sgndetGdn,
            hst, fpi_up, fpi_dn, Bup, Bdn, n, l, rng,
        )
        accepted += acc
    end
    acceptance_rate = accepted / (2 * length(hst.update_perm))
    return acceptance_rate, logdetGup, sgndetGup, logdetGdn, sgndetGdn
end

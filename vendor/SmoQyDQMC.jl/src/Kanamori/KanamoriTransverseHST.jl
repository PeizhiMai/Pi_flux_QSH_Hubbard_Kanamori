@doc raw"""
    KanamoriTransverseHST{T,R} <: AbstractAsymHST{T,R}

Exact three-state auxiliary-field decomposition for the combined rotationally invariant
Kanamori transverse term

```math
-J_H(S_s^+S_p^- + S_s^-S_p^+) + J_H(P_s^\dagger P_p + P_p^\dagger P_s)
= \frac{J_H}{4}(O_+^2 - O_-^2),
```

where ``O_+ = F_\uparrow + F_\downarrow``, ``O_- = F_\uparrow - F_\downarrow``, and
``F_\sigma = c^\dagger_{s\sigma}c_{p\sigma}+c^\dagger_{p\sigma}c_{s\sigma}``.

This HST preserves the separate spin-up/spin-down Green's function structure but inserts
local off-diagonal orbital factors. It assumes asymmetric exact propagators and time-slice ordering

```text
B_l = exp(-Δτ V_l) W_+,l W_-,l exp(-Δτ K_l).
```
"""
struct KanamoriTransverseHST{T,R} <: AbstractAsymHST{T,R}
    β::R
    Δτ::R
    Lτ::Int
    N::Int
    JH::Vector{R}

    λ_plus::Vector{T}
    weight0_plus::Vector{R}
    weight1_plus::Vector{R}

    λ_minus::Vector{T}
    weight0_minus::Vector{R}
    weight1_minus::Vector{R}

    neighbor_table::Matrix{Int}
    orbital_pairs::Vector{Tuple{Int,Int}}
    s_plus::Array{Int,2}
    s_minus::Array{Int,2}
    update_perm::Vector{Int}
end

function kanamori_square_hs_coefficients(κ::R) where {R<:AbstractFloat}
    if iszero(κ)
        return zero(Complex{R}), one(R), zero(R)
    end

    y = (exp(4κ) - 1) / (2 * (exp(κ) - 1)) - 1
    λ = acosh(complex(R(y), zero(R)))
    B = real((exp(κ) - 1) / (cosh(λ) - 1))
    A = 1 - B
    @assert A > 0 && B > 0 "Non-positive Kanamori transverse HS weights: A=$A, B=$B for κ=$κ."
    return Complex{R}(λ), R(A), R(B/2)
end

function kanamori_transverse_field_action(s::Int, weight0::R, weight1::R) where {R<:AbstractFloat}
    @assert s in (-1, 0, 1)
    w = iszero(s) ? weight0 : weight1
    return -log(w)
end

function _kanamori_transverse_local_W(
    channel::Symbol, s::Int, λ::T, spin_sign::Int, ::Type{E}=typeof(λ)
) where {T<:Number,E<:Number}
    @assert channel === :plus || channel === :minus
    @assert spin_sign == 1 || spin_sign == -1
    sign = channel === :plus ? 1 : spin_sign
    a = sign * s * λ
    I2 = Matrix{E}(I, 2, 2)
    F = _local_orbital_generator(Val(:F), E)
    return cosh(a) * I2 + sinh(a) * F
end

function _kanamori_transverse_field_action(hst::KanamoriTransverseHST, channel::Symbol, n::Int, s::Int)
    if channel === :plus
        return kanamori_transverse_field_action(s, hst.weight0_plus[n], hst.weight1_plus[n])
    else
        @assert channel === :minus
        return kanamori_transverse_field_action(s, hst.weight0_minus[n], hst.weight1_minus[n])
    end
end

@doc raw"""
    KanamoriTransverseHST(; model_geometry, orbital_pairs=[(1,2)], JH, β, Δτ, rng)

Construct HS fields for the combined spin-flip plus pair-hopping Kanamori transverse term.
The first implementation assumes the rotationally invariant choice ``J_pair = J_H``.
"""
function KanamoriTransverseHST(;
    model_geometry::ModelGeometry{D,R},
    orbital_pairs::AbstractVector{<:Tuple{Int,Int}} = [(1,2)],
    JH::R,
    β::R, Δτ::R, rng::AbstractRNG,
) where {D,R<:AbstractFloat}
    @assert JH >= 0 "KanamoriTransverseHST currently requires JH >= 0."

    (; lattice, unit_cell) = model_geometry
    N_unitcells = lattice.N
    n_pairs = length(orbital_pairs)
    N = N_unitcells * n_pairs

    JH_vec = fill(JH, N)
    λ_plus = zeros(Complex{R}, N)
    λ_minus = zeros(Complex{R}, N)
    weight0_plus = zeros(R, N)
    weight1_plus = zeros(R, N)
    weight0_minus = zeros(R, N)
    weight1_minus = zeros(R, N)

    for n in eachindex(JH_vec)
        λ_plus[n], weight0_plus[n], weight1_plus[n] = kanamori_square_hs_coefficients(-Δτ * JH_vec[n] / 4)
        λ_minus[n], weight0_minus[n], weight1_minus[n] = kanamori_square_hs_coefficients(+Δτ * JH_vec[n] / 4)
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
    if iszero(JH)
        s_plus = zeros(Int, N, Lτ)
        s_minus = zeros(Int, N, Lτ)
    else
        s_plus = rand(rng, (-1,0,1), (N, Lτ))
        s_minus = rand(rng, (-1,0,1), (N, Lτ))
    end
    update_perm = collect(1:N)

    return KanamoriTransverseHST{Complex{R},R}(
        β, Δτ, Lτ, N, JH_vec,
        λ_plus, weight0_plus, weight1_plus,
        λ_minus, weight0_minus, weight1_minus,
        neighbor_table, collect(Tuple{Int,Int}, orbital_pairs),
        s_plus, s_minus, update_perm,
    )
end

function _initialize!(
    fermion_path_integral_up::FermionPathIntegral{H},
    fermion_path_integral_dn::FermionPathIntegral{H},
    hst::KanamoriTransverseHST{T,R},
) where {H<:Number,T<:Number,R<:AbstractFloat}
    @assert !((H<:Real) && (T<:Complex)) "Green's function matrices are real while KanamoriTransverseHST is complex. Use complex BHZ hoppings or force complex path integrals."
    @assert fermion_path_integral_up.Sb == fermion_path_integral_dn.Sb "$(fermion_path_integral_up.Sb) ≠ $(fermion_path_integral_dn.Sb)"

    (; s_plus, s_minus, Lτ) = hst
    for n in axes(hst.neighbor_table, 2)
        for l in 1:Lτ
            fermion_path_integral_up.Sb += _kanamori_transverse_field_action(hst, :plus, n, s_plus[n,l])
            fermion_path_integral_up.Sb += _kanamori_transverse_field_action(hst, :minus, n, s_minus[n,l])
        end
    end
    fermion_path_integral_dn.Sb = fermion_path_integral_up.Sb

    return nothing
end

@doc raw"""
    apply_kanamori_transverse_to_propagators!(Bup, Bdn, hst)

Insert current pair-hopping/spin-flip transverse one-body factors into asymmetric exact propagators.
Call this immediately after `initialize_propagators` and before constructing/calculating Green's functions.
"""
function apply_kanamori_transverse_to_propagators!(
    Bup::AbstractVector{<:AsymExactPropagator},
    Bdn::AbstractVector{<:AsymExactPropagator},
    hst::KanamoriTransverseHST,
)
    (; neighbor_table, s_plus, s_minus, λ_plus, λ_minus, Lτ) = hst
    @assert length(Bup) == Lτ && length(Bdn) == Lτ
    for l in 1:Lτ
        for n in axes(neighbor_table, 2)
            inds = @view neighbor_table[:,n]
            Wp_up = _kanamori_transverse_local_W(:plus, s_plus[n,l], λ_plus[n], +1)
            Wm_up = _kanamori_transverse_local_W(:minus, s_minus[n,l], λ_minus[n], +1)
            Wp_dn = _kanamori_transverse_local_W(:plus, s_plus[n,l], λ_plus[n], -1)
            Wm_dn = _kanamori_transverse_local_W(:minus, s_minus[n,l], λ_minus[n], -1)
            _apply_local_left_factor_to_asym_exact!(Bup[l], Wp_up * Wm_up, inds)
            _apply_local_left_factor_to_asym_exact!(Bdn[l], Wp_dn * Wm_dn, inds)
        end
    end
    return nothing
end

function _attempt_kanamori_transverse_update!(
    channel::Symbol,
    Gup::Matrix{H}, logdetGup::R, sgndetGup::H,
    Gdn::Matrix{H}, logdetGdn::R, sgndetGdn::H,
    hst::KanamoriTransverseHST{T,R},
    fpi_up::FermionPathIntegral{H}, fpi_dn::FermionPathIntegral{H},
    Bup::AsymExactPropagator, Bdn::AsymExactPropagator,
    n::Int, l::Int, rng::AbstractRNG,
) where {H<:Number,T<:Number,R<:Real}
    fields = channel === :plus ? hst.s_plus : hst.s_minus
    sold = fields[n,l]
    snew = sample_new_hund_spinflip_field(rng, sold)
    Δs = snew - sold
    ΔSb = _kanamori_transverse_field_action(hst, channel, n, snew) -
          _kanamori_transverse_field_action(hst, channel, n, sold)
    inds = @view hst.neighbor_table[:,n]

    if channel === :plus
        Xw_up = _kanamori_transverse_local_W(:plus, Δs, hst.λ_plus[n], +1)
        Xw_dn = _kanamori_transverse_local_W(:plus, Δs, hst.λ_plus[n], -1)
    else
        @assert channel === :minus
        Xw_up = _kanamori_transverse_local_W(:minus, Δs, hst.λ_minus[n], +1)
        Xw_dn = _kanamori_transverse_local_W(:minus, Δs, hst.λ_minus[n], -1)
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
    hst::KanamoriTransverseHST{T,R},
    fpi_up::FermionPathIntegral{H},
    fpi_dn::FermionPathIntegral{H},
    Bup::P, Bdn::P, l::Int, rng::AbstractRNG,
) where {H<:Number,T<:Number,R<:Real,P<:AsymExactPropagator}
    shuffle!(rng, hst.update_perm)
    accepted = 0
    for n in hst.update_perm
        acc, logdetGup, sgndetGup, logdetGdn, sgndetGdn = _attempt_kanamori_transverse_update!(
            :plus, Gup, logdetGup, sgndetGup, Gdn, logdetGdn, sgndetGdn,
            hst, fpi_up, fpi_dn, Bup, Bdn, n, l, rng,
        )
        accepted += acc
        acc, logdetGup, sgndetGup, logdetGdn, sgndetGdn = _attempt_kanamori_transverse_update!(
            :minus, Gup, logdetGup, sgndetGup, Gdn, logdetGdn, sgndetGdn,
            hst, fpi_up, fpi_dn, Bup, Bdn, n, l, rng,
        )
        accepted += acc
    end
    acceptance_rate = accepted / (2 * length(hst.update_perm))
    return acceptance_rate, logdetGup, sgndetGup, logdetGdn, sgndetGdn
end

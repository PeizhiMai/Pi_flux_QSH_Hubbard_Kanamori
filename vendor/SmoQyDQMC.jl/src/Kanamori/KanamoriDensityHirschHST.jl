@doc raw"""
    KanamoriDensityHirschHST{T,R} <: AbstractAsymHST{T,R}

Spin-channel Hirsch HS transformation for the density-density inter-orbital Kanamori terms
represented by [`KanamoriDensityParameters`](@ref).

The four spin-resolved density couplings are updated independently:

- same spin: ``V_{same}`` for ``\uparrow\uparrow`` and ``\downarrow\downarrow``;
- opposite spin: ``V_{opposite}`` for ``\uparrow\downarrow`` and ``\downarrow\uparrow``.
"""
struct KanamoriDensityHirschHST{T,R} <: AbstractAsymHST{T,R}

    β::R
    Δτ::R
    Lτ::Int
    N::Int

    V_same::Vector{R}
    V_opposite::Vector{R}

    α_same::Vector{T}
    α_opposite::Vector{T}

    neighbor_table::Matrix{Int}

    s_upup::Array{Int,2}
    s_dndn::Array{Int,2}
    s_updn::Array{Int,2}
    s_dnup::Array{Int,2}

    update_perm::Vector{Int}
end

@doc raw"""
    KanamoriDensityHirschHST(; kanamori_density_parameters, β, Δτ, rng)

Initialize the HS fields for density-density Kanamori interactions.
"""
function KanamoriDensityHirschHST(;
    kanamori_density_parameters::KanamoriDensityParameters{R},
    β::R, Δτ::R, rng::AbstractRNG,
) where {R<:AbstractFloat}

    (; V_same, V_opposite, neighbor_table) = kanamori_density_parameters

    T = (any(v -> v < 0, V_same) || any(v -> v < 0, V_opposite)) ? Complex{R} : R
    Lτ = round(Int, β / Δτ)
    N = length(V_same)

    α_same = zeros(T, N)
    α_opposite = zeros(T, N)
    @. α_same = acosh(exp(Δτ*T(V_same)/2))/Δτ
    @. α_opposite = acosh(exp(Δτ*T(V_opposite)/2))/Δτ

    s_upup = rand(rng, -1:2:1, (N, Lτ))
    s_dndn = rand(rng, -1:2:1, (N, Lτ))
    s_updn = rand(rng, -1:2:1, (N, Lτ))
    s_dnup = rand(rng, -1:2:1, (N, Lτ))
    update_perm = collect(1:N)

    return KanamoriDensityHirschHST{T,R}(
        β, Δτ, Lτ, N,
        V_same, V_opposite,
        α_same, α_opposite,
        neighbor_table,
        s_upup, s_dndn, s_updn, s_dnup,
        update_perm,
    )
end

# initialize fermion path integral to reflect HS field config
function _initialize!(
    fermion_path_integral_up::FermionPathIntegral{H},
    fermion_path_integral_dn::FermionPathIntegral{H},
    hst_parameters::KanamoriDensityHirschHST{T},
) where {H<:Number, T<:Number}

    @assert !((H<:Real) && (T<:Complex)) "Green's function matrices are real while KanamoriDensityHirschHST is complex."
    @assert fermion_path_integral_up.Sb == fermion_path_integral_dn.Sb "$(fermion_path_integral_up.Sb) ≠ $(fermion_path_integral_dn.Sb)"

    (; neighbor_table, α_same, α_opposite, s_upup, s_dndn, s_updn, s_dnup) = hst_parameters
    Vup = fermion_path_integral_up.V
    Vdn = fermion_path_integral_dn.V

    for b in axes(neighbor_table, 2)
        i, j = neighbor_table[1,b], neighbor_table[2,b]

        @views @. Vup[i,:] += -α_same[b] * s_upup[b,:]
        @views @. Vup[j,:] += +α_same[b] * s_upup[b,:]

        @views @. Vdn[i,:] += -α_same[b] * s_dndn[b,:]
        @views @. Vdn[j,:] += +α_same[b] * s_dndn[b,:]

        @views @. Vdn[i,:] += -α_opposite[b] * s_updn[b,:]
        @views @. Vup[j,:] += +α_opposite[b] * s_updn[b,:]

        @views @. Vup[i,:] += -α_opposite[b] * s_dnup[b,:]
        @views @. Vdn[j,:] += +α_opposite[b] * s_dnup[b,:]
    end

    return nothing
end

# perform local updates for specified imaginary-time slice
function _local_updates!(
    Gup::Matrix{H}, logdetGup::R, sgndetGup::H,
    Gdn::Matrix{H}, logdetGdn::R, sgndetGdn::H,
    hst_parameters::KanamoriDensityHirschHST{T,R},
    fermion_path_integral_up::FermionPathIntegral{H},
    fermion_path_integral_dn::FermionPathIntegral{H},
    Bup::P, Bdn::P, l::Int, rng::AbstractRNG,
) where {H<:Number, T<:Number, R<:Real, P<:AbstractPropagator}

    (; Δτ, α_same, α_opposite, neighbor_table, update_perm) = hst_parameters
    u′ = @view fermion_path_integral_up.u[:,1]
    v′ = @view fermion_path_integral_up.v[:,1]
    u″ = @view fermion_path_integral_dn.u[:,1:2]
    v″ = @view fermion_path_integral_dn.v[:,1:2]

    Vup = @view fermion_path_integral_up.V[:,l]
    Vdn = @view fermion_path_integral_dn.V[:,l]

    s_upup = @view hst_parameters.s_upup[:,l]
    s_dndn = @view hst_parameters.s_dndn[:,l]
    s_updn = @view hst_parameters.s_updn[:,l]
    s_dnup = @view hst_parameters.s_dnup[:,l]

    shuffle!(rng, update_perm)
    accepted_spin_flips = zero(Int)

    for n in update_perm
        i, j = neighbor_table[1,n], neighbor_table[2,n]

        accepted, logdetGup, sgndetGup = _local_update!(
            Gup, logdetGup, sgndetGup, Bup, Vup, s_upup,
            n, j, i, Δτ, α_same, rng, u″, v″,
        )
        accepted_spin_flips += accepted

        accepted, logdetGdn, sgndetGdn = _local_update!(
            Gdn, logdetGdn, sgndetGdn, Bdn, Vdn, s_dndn,
            n, j, i, Δτ, α_same, rng, u″, v″,
        )
        accepted_spin_flips += accepted

        accepted, logdetGup, sgndetGup, logdetGdn, sgndetGdn = _local_update!(
            Gup, logdetGup, sgndetGup, Bup, Vup,
            Gdn, logdetGdn, sgndetGdn, Bdn, Vdn,
            s_updn, n, j, i, Δτ, α_opposite, rng, u′, v′,
        )
        accepted_spin_flips += accepted

        accepted, logdetGdn, sgndetGdn, logdetGup, sgndetGup = _local_update!(
            Gdn, logdetGdn, sgndetGdn, Bdn, Vdn,
            Gup, logdetGup, sgndetGup, Bup, Vup,
            s_dnup, n, j, i, Δτ, α_opposite, rng, u′, v′,
        )
        accepted_spin_flips += accepted
    end

    acceptance_rate = accepted_spin_flips / (4 * length(s_upup))

    return acceptance_rate, logdetGup, sgndetGup, logdetGdn, sgndetGdn
end

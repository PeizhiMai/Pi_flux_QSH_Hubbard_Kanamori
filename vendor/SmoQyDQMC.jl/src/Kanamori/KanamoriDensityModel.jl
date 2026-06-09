@doc raw"""
    KanamoriDensityModel{T<:AbstractFloat}

Density-density (Ising-Hund) inter-orbital part of a two-orbital Kanamori interaction.

For each orbital pair ``(a,b)`` in a unit cell this model represents
```math
H_{ab}^{dd} = V_{\rm same}\left(n_{a\uparrow}n_{b\uparrow} + n_{a\downarrow}n_{b\downarrow}\right)
            + V_{\rm opposite}\left(n_{a\uparrow}n_{b\downarrow} + n_{a\downarrow}n_{b\uparrow}\right).
```
For rotationally invariant Kanamori with ``U'=U-2J_H`` and with spin-flip/pair-hopping terms omitted,
```math
V_{\rm opposite}=U'=U-2J_H,\qquad V_{\rm same}=U'-J_H=U-3J_H.
```

The intra-orbital Hubbard term ``U n_{a\uparrow}n_{a\downarrow}`` is intentionally not included here;
use the existing [`HubbardModel`](@ref) and `HubbardSpinHirschHST` for that part.
"""
struct KanamoriDensityModel{T<:AbstractFloat}

    # whether particle-hole symmetric form for interaction is used
    ph_sym_form::Bool

    # local orbital pairs in each unit cell
    orbital_pairs::Vector{Tuple{Int,Int}}

    # same-spin inter-orbital density coupling
    V_same_mean::Vector{T}
    V_same_std::Vector{T}

    # opposite-spin inter-orbital density coupling
    V_opposite_mean::Vector{T}
    V_opposite_std::Vector{T}
end

@doc raw"""
    KanamoriDensityModel(; ph_sym_form, orbital_pairs=[(1,2)], U, JH,
                           V_same=U-3JH, V_opposite=U-2JH,
                           V_same_std=zeros(...), V_opposite_std=zeros(...))

Construct the density-density inter-orbital Kanamori model.  By default it uses the
rotationally invariant relations ``U'=U-2J_H`` and ``U'-J_H=U-3J_H``.
"""
function KanamoriDensityModel(;
    ph_sym_form::Bool,
    orbital_pairs::AbstractVector{<:Tuple{Int,Int}} = [(1,2)],
    U::T,
    JH::T,
    V_same::AbstractVector{T} = fill(U - 3JH, length(orbital_pairs)),
    V_opposite::AbstractVector{T} = fill(U - 2JH, length(orbital_pairs)),
    V_same_std::AbstractVector{T} = zero(V_same),
    V_opposite_std::AbstractVector{T} = zero(V_opposite),
) where {T<:AbstractFloat}

    @assert length(V_same) == length(orbital_pairs)
    @assert length(V_opposite) == length(orbital_pairs)
    @assert length(V_same_std) == length(orbital_pairs)
    @assert length(V_opposite_std) == length(orbital_pairs)

    return KanamoriDensityModel{T}(
        ph_sym_form,
        collect(Tuple{Int,Int}, orbital_pairs),
        collect(T, V_same),
        collect(T, V_same_std),
        collect(T, V_opposite),
        collect(T, V_opposite_std),
    )
end

# show struct info as TOML formatted string
function Base.show(io::IO, ::MIME"text/plain", kdm::KanamoriDensityModel)

    (; ph_sym_form, orbital_pairs, V_same_mean, V_same_std, V_opposite_mean, V_opposite_std) = kdm

    @printf io "[KanamoriDensityModel]\n\n"
    @printf io "KANAMORI_DENSITY_IDS = %s\n" string(collect(1:length(orbital_pairs)))
    @printf io "ORBITAL_PAIRS        = %s\n" string(orbital_pairs)
    @printf io "V_same_mean         = %s\n" string(round.(V_same_mean, digits=6))
    @printf io "V_same_std          = %s\n" string(round.(V_same_std, digits=6))
    @printf io "V_opposite_mean     = %s\n" string(round.(V_opposite_mean, digits=6))
    @printf io "V_opposite_std      = %s\n" string(round.(V_opposite_std, digits=6))
    @printf io "ph_sym_form         = %s\n\n" string(ph_sym_form)

    return nothing
end

@doc raw"""
    KanamoriDensityParameters{T<:AbstractFloat}

Finite-lattice parameters for [`KanamoriDensityModel`](@ref).
"""
struct KanamoriDensityParameters{T<:AbstractFloat}

    # same-spin coupling for each local orbital pair in the lattice
    V_same::Vector{T}

    # opposite-spin coupling for each local orbital pair in the lattice
    V_opposite::Vector{T}

    # local orbital-pair neighbor table: row 1 = first orbital, row 2 = second orbital
    neighbor_table::Matrix{Int}

    # orbital pairs in each unit cell
    orbital_pairs::Vector{Tuple{Int,Int}}

    # whether particle-hole symmetric form for interaction is used
    ph_sym_form::Bool
end

@doc raw"""
    KanamoriDensityParameters(; kanamori_density_model, model_geometry, rng)

Initialize finite-lattice density-density Kanamori parameters.
"""
function KanamoriDensityParameters(;
    kanamori_density_model::KanamoriDensityModel{T},
    model_geometry::ModelGeometry{D,T},
    rng::AbstractRNG,
) where {D, T<:AbstractFloat}

    (; orbital_pairs, V_same_mean, V_same_std, V_opposite_mean, V_opposite_std, ph_sym_form) = kanamori_density_model
    (; lattice, unit_cell) = model_geometry

    N_unitcells = lattice.N
    n_pairs = length(orbital_pairs)
    N = N_unitcells * n_pairs

    V_same = zeros(T, N)
    V_opposite = zeros(T, N)
    neighbor_table = zeros(Int, 2, N)

    V_same′ = reshape(V_same, (N_unitcells, n_pairs))
    V_opposite′ = reshape(V_opposite, (N_unitcells, n_pairs))
    nt′ = reshape(neighbor_table, (2, N_unitcells, n_pairs))

    for (p, (orb_a, orb_b)) in enumerate(orbital_pairs)
        for u in 1:N_unitcells
            nt′[1,u,p] = loc_to_site(u, orb_a, unit_cell)
            nt′[2,u,p] = loc_to_site(u, orb_b, unit_cell)
            V_same′[u,p] = V_same_mean[p] + V_same_std[p] * randn(rng)
            V_opposite′[u,p] = V_opposite_mean[p] + V_opposite_std[p] * randn(rng)
        end
    end

    return KanamoriDensityParameters{T}(V_same, V_opposite, neighbor_table, orbital_pairs, ph_sym_form)
end

@doc raw"""
    initialize!(fermion_path_integral_up, fermion_path_integral_dn, kanamori_density_parameters)

Apply one-body shifts associated with the particle-hole asymmetric density-density Kanamori form.
"""
function initialize!(
    fermion_path_integral_up::FermionPathIntegral,
    fermion_path_integral_dn::FermionPathIntegral,
    kanamori_density_parameters::KanamoriDensityParameters,
)

    initialize!(fermion_path_integral_up, kanamori_density_parameters, spin=+1)
    initialize!(fermion_path_integral_dn, kanamori_density_parameters, spin=-1)

    return nothing
end

function initialize!(
    fermion_path_integral::FermionPathIntegral,
    kanamori_density_parameters::KanamoriDensityParameters;
    spin::Int,
)

    (; ph_sym_form, neighbor_table, V_same, V_opposite) = kanamori_density_parameters
    (; V) = fermion_path_integral

    @assert spin == +1 || spin == -1

    # Convert n_iσ n_jσ′ to (n_iσ-1/2)(n_jσ′-1/2) plus one-body shifts.
    # For spin up/down at a given orbital, one same-spin term and one opposite-spin term touch it.
    if !ph_sym_form
        for l in axes(V,2)
            for n in eachindex(V_same)
                i = neighbor_table[1,n]
                j = neighbor_table[2,n]
                shift = 0.5 * (V_same[n] + V_opposite[n])
                V[i,l] += shift
                V[j,l] += shift
            end
        end
    end

    return nothing
end

@doc raw"""
    measure_kanamori_density_energy(params, Gup, Gdn, id)

Measure the density-density inter-orbital Kanamori energy for orbital-pair type `id`.
"""
function measure_kanamori_density_energy(
    params::KanamoriDensityParameters{E},
    Gup::Matrix{T}, Gdn::Matrix{T},
    id::Int,
) where {T<:Number, E<:AbstractFloat}

    (; V_same, V_opposite, neighbor_table, orbital_pairs, ph_sym_form) = params

    n_pairs = length(orbital_pairs)
    N_unitcells = size(neighbor_table, 2) ÷ n_pairs

    Vsame′ = reshape(V_same, (N_unitcells, n_pairs))
    Vopp′ = reshape(V_opposite, (N_unitcells, n_pairs))
    nt = reshape(neighbor_table, (2, N_unitcells, n_pairs))

    ϵ = zero(T)
    for u in 1:N_unitcells
        i = nt[1,u,id]
        j = nt[2,u,id]
        niu = 1 - Gup[i,i]
        nid = 1 - Gdn[i,i]
        nju = 1 - Gup[j,j]
        njd = 1 - Gdn[j,j]
        if ph_sym_form
            niu -= 0.5; nid -= 0.5; nju -= 0.5; njd -= 0.5
        end
        same_up = niu*nju - Gup[j,i]*Gup[i,j]
        same_dn = nid*njd - Gdn[j,i]*Gdn[i,j]
        ϵ += Vsame′[u,id] * (same_up + same_dn)
        ϵ += Vopp′[u,id] * (niu*njd + nid*nju)
    end

    return ϵ / N_unitcells
end

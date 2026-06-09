@doc raw"""
    local_matrix_update_det_ratio(G, Δloc, inds)

Return the determinant ratio for a local left-multiplicative one-body update
``B_l -> (I + Δ) B_l`` affecting the indices `inds`.

`Δloc` is the local matrix ``Δ[inds, inds]``.  The returned ratio follows the
SmoQy/JDQMC convention

```math
R = \det G / \det G' = \det\{I + (I-G)_{S,S} Δ_S\}.
```
"""
function local_matrix_update_det_ratio(
    G::AbstractMatrix{T},
    Δloc::AbstractMatrix{E},
    inds::AbstractVector{Int},
) where {T<:Number, E<:Number}

    k = length(inds)
    @assert size(Δloc) == (k, k)

    M = Matrix{promote_type(T,E)}(I, k, k)
    @inbounds for a in 1:k, b in 1:k
        ig = (a == b ? one(T) : zero(T)) - G[inds[a], inds[b]]
        for c in 1:k
            M[a,c] += ig * Δloc[b,c]
        end
    end

    return det(M)
end

@doc raw"""
    local_matrix_update_greens!(G, logdetG, sgndetG, R, Δloc, inds)

Apply the accepted local matrix update corresponding to [`local_matrix_update_det_ratio`](@ref).
This updates only the Green's function and determinant bookkeeping; the caller is responsible for
updating any stored time-slice propagator/HS-field representation.
"""
function local_matrix_update_greens!(
    G::AbstractMatrix{T},
    logdetG::R,
    sgndetG::T,
    ratio::E,
    Δloc::AbstractMatrix{F},
    inds::AbstractVector{Int},
) where {T<:Number, R<:Real, E<:Number, F<:Number}

    k = length(inds)
    @assert size(Δloc) == (k, k)

    H = promote_type(T,F)
    M = Matrix{H}(I, k, k)
    IG_rows = zeros(H, k, size(G,2))

    @inbounds for a in 1:k
        ia = inds[a]
        for col in axes(G,2)
            IG_rows[a,col] = (ia == col ? one(H) : zero(H)) - G[ia,col]
        end
        for b in 1:k
            ib = inds[b]
            ig = (ia == ib ? one(H) : zero(H)) - G[ia,ib]
            for c in 1:k
                M[a,c] += ig * Δloc[b,c]
            end
        end
    end

    # G' = G - G[:,S] Δ [I + (I-G)_{S,S} Δ]^{-1} (I-G)[S,:]
    G_cols = G[:, inds]
    X = (Matrix{H}(G_cols) * Matrix{H}(Δloc)) / M
    G .-= X * IG_rows

    invR = inv(ratio)
    logdetG′ = logdetG + log(abs(invR))
    sgndetG′ = sign(invR) * sgndetG

    return logdetG′, sgndetG′
end

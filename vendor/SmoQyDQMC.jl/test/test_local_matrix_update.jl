@testitem "Local matrix update kernel" begin

    using LinearAlgebra
    using Random
    using SmoQyDQMC

    rng = MersenneTwister(1234)

    function random_well_conditioned_product(rng, N)
        A = Matrix{ComplexF64}(I, N, N) + 0.15 * (randn(rng, ComplexF64, N, N))
        return A
    end

    for N in (5, 8), inds in ([2], [2,4], [1,3,5])
        local A = random_well_conditioned_product(rng, N)
        local G = inv(Matrix{ComplexF64}(I, N, N) + A)
        local Gfast = copy(G)

        k = length(inds)
        δh = 0.05 * randn(rng, ComplexF64, k, k)
        Δloc = exp(δh) - Matrix{ComplexF64}(I, k, k)

        local R = local_matrix_update_det_ratio(G, Δloc, inds)
        local A′ = copy(A)
        A′[inds, :] .+= Δloc * A[inds, :]
        Gfull = inv(Matrix{ComplexF64}(I, N, N) + A′)

        expected_R = det(Matrix{ComplexF64}(I, N, N) + A′) / det(Matrix{ComplexF64}(I, N, N) + A)
        @test R ≈ expected_R atol=1e-11 rtol=1e-11

        logdetG = log(abs(det(G)))
        sgndetG = det(G) / abs(det(G))
        logdetG′, sgndetG′ = local_matrix_update_greens!(Gfast, logdetG, sgndetG, R, Δloc, inds)

        @test Gfast ≈ Gfull atol=1e-11 rtol=1e-11
        @test logdetG′ ≈ log(abs(det(Gfull))) atol=1e-10
        @test sgndetG′ ≈ det(Gfull) / abs(det(Gfull)) atol=1e-10
    end

    # Diagonal one-site limit reduces to the JDQMC/SmoQy determinant ratio.
    A = random_well_conditioned_product(rng, 6)
    G = inv(Matrix{ComplexF64}(I, 6, 6) + A)
    i = 3
    Δ = 0.2 + 0.1im
    R = local_matrix_update_det_ratio(G, reshape([Δ], 1, 1), [i])
    @test R ≈ 1 + Δ * (1 - G[i,i]) atol=1e-12
end

@testitem "Atomic ED validation for density plus transverse Kanamori Hund" begin

    using LinearAlgebra

    function ed_annihilation(nflavors::Int, flavor::Int)
        dim = 1 << nflavors
        c = zeros(ComplexF64, dim, dim)
        mask_below = (1 << (flavor - 1)) - 1
        for state in 0:(dim - 1)
            if ((state >> (flavor - 1)) & 1) == 1
                occupied_below = count_ones(state & mask_below)
                sign = isodd(occupied_below) ? -1.0 : 1.0
                new_state = state & ~(1 << (flavor - 1))
                c[new_state + 1, state + 1] = sign
            end
        end
        return c
    end

    ed_flavor(cell::Int, local_flavor::Int) = 4 * (cell - 1) + local_flavor

    function ed_density_spinflip_hund(ncells::Int; U::Float64, JH::Float64, pair_hopping::Bool=false)
        nflavors = 4 * ncells
        dim = 1 << nflavors
        c = [ed_annihilation(nflavors, f) for f in 1:nflavors]
        cd = adjoint.(c)
        n = [cd[f] * c[f] for f in 1:nflavors]
        H = zeros(ComplexF64, dim, dim)
        Nop = zeros(ComplexF64, dim, dim)

        for f in 1:nflavors
            Nop += n[f]
        end

        for cell in 1:ncells
            s_up = ed_flavor(cell, 1)
            p_up = ed_flavor(cell, 2)
            s_dn = ed_flavor(cell, 3)
            p_dn = ed_flavor(cell, 4)

            H += U * (n[s_up] * n[s_dn] + n[p_up] * n[p_dn])
            H += (U - 3JH) * (n[s_up] * n[p_up] + n[s_dn] * n[p_dn])
            H += (U - 2JH) * (n[s_up] * n[p_dn] + n[s_dn] * n[p_up])

            Ssp = cd[s_up] * c[s_dn]
            Ssm = cd[s_dn] * c[s_up]
            Spp = cd[p_up] * c[p_dn]
            Spm = cd[p_dn] * c[p_up]
            H += -JH * (Ssp * Spm + Ssm * Spp)

            if pair_hopping
                PsdagPp = cd[s_up] * cd[s_dn] * c[p_dn] * c[p_up]
                PpdagPs = cd[p_up] * cd[p_dn] * c[s_dn] * c[s_up]
                H += JH * (PsdagPp + PpdagPs)
            end
        end

        return Hermitian(H), Hermitian(Nop)
    end

    function sector_indices(nflavors::Int, nparticles::Int)
        return [state + 1 for state in 0:((1 << nflavors) - 1) if count_ones(state) == nparticles]
    end

    function grand_canonical_stats(H::Hermitian, Nop::Hermitian; β::Float64, μ::Float64)
        K = Hermitian(Matrix(H) - μ * Matrix(Nop))
        F = eigen(K)
        weights = exp.(-β .* F.values)
        Z = sum(weights)
        H_eig = real.(diag(F.vectors' * Matrix(H) * F.vectors))
        N_eig = real.(diag(F.vectors' * Matrix(Nop) * F.vectors))
        E = sum(weights .* H_eig) / Z
        N = sum(weights .* N_eig) / Z
        return Z, E, N
    end

    U = 2.0
    JH = 0.3
    H1, N1 = ed_density_spinflip_hund(1; U=U, JH=JH)
    H2, N2 = ed_density_spinflip_hund(2; U=U, JH=JH)

    @test norm(Matrix(H1) - Matrix(H1)') < 1e-12
    @test norm(Matrix(H2) - Matrix(H2)') < 1e-12

    n2 = sector_indices(4, 2)
    evals_n2 = sort(real.(eigvals(Hermitian(Matrix(H1)[n2, n2]))))
    expected_n2 = sort([U - 3JH, U - 3JH, U - 3JH, U - JH, U, U])
    @test evals_n2 ≈ expected_n2 atol=1e-12

    H1_full, N1_full = ed_density_spinflip_hund(1; U=U, JH=JH, pair_hopping=true)
    H2_full, N2_full = ed_density_spinflip_hund(2; U=U, JH=JH, pair_hopping=true)

    @test norm(Matrix(H1_full) - Matrix(H1_full)') < 1e-12
    @test norm(Matrix(H2_full) - Matrix(H2_full)') < 1e-12

    evals_n2_full = sort(real.(eigvals(Hermitian(Matrix(H1_full)[n2, n2]))))
    expected_n2_full = sort([U - 3JH, U - 3JH, U - 3JH, U - JH, U - JH, U + JH])
    @test evals_n2_full ≈ expected_n2_full atol=1e-12

    β = 1.3
    μ = 0.2
    Z1, E1, Navg1 = grand_canonical_stats(H1, N1; β=β, μ=μ)
    Z2, E2, Navg2 = grand_canonical_stats(H2, N2; β=β, μ=μ)

    @test Z2 ≈ Z1^2 rtol=1e-11
    @test E2 ≈ 2 * E1 rtol=1e-11 atol=1e-11
    @test Navg2 ≈ 2 * Navg1 rtol=1e-11 atol=1e-11

    Z1_full, E1_full, Navg1_full = grand_canonical_stats(H1_full, N1_full; β=β, μ=μ)
    Z2_full, E2_full, Navg2_full = grand_canonical_stats(H2_full, N2_full; β=β, μ=μ)

    @test Z2_full ≈ Z1_full^2 rtol=1e-11
    @test E2_full ≈ 2 * E1_full rtol=1e-11 atol=1e-11
    @test Navg2_full ≈ 2 * Navg1_full rtol=1e-11 atol=1e-11

    # At JH=0 the spin-flip operator is absent and the one-electron-per-orbital
    # two-particle multiplet collapses to four degenerate states at energy U.
    H1_j0, _ = ed_density_spinflip_hund(1; U=U, JH=0.0)
    evals_n2_j0 = sort(real.(eigvals(Hermitian(Matrix(H1_j0)[n2, n2]))))
    @test evals_n2_j0 ≈ fill(U, 6) atol=1e-12

    H1_full_j0, _ = ed_density_spinflip_hund(1; U=U, JH=0.0, pair_hopping=true)
    evals_n2_full_j0 = sort(real.(eigvals(Hermitian(Matrix(H1_full_j0)[n2, n2]))))
    @test evals_n2_full_j0 ≈ fill(U, 6) atol=1e-12
end

@testitem "Kanamori local spin-flip algebra" begin

    using LinearAlgebra

    const DIM_KANAMORI_LOCAL = 16

    function local_annihilation(flavor::Int)
        c = zeros(ComplexF64, DIM_KANAMORI_LOCAL, DIM_KANAMORI_LOCAL)
        mask_below = (1 << (flavor - 1)) - 1
        for state in 0:(DIM_KANAMORI_LOCAL - 1)
            if ((state >> (flavor - 1)) & 1) == 1
                occupied_below = count_ones(state & mask_below)
                sign = isodd(occupied_below) ? -1.0 : 1.0
                new_state = state & ~(1 << (flavor - 1))
                c[new_state + 1, state + 1] = sign
            end
        end
        return c
    end

    function local_kanamori_operators()
        c = [local_annihilation(f) for f in 1:4]
        cd = adjoint.(c)
        n = [cd[f] * c[f] for f in 1:4]

        # Flavor order: 1=s↑, 2=p↑, 3=s↓, 4=p↓.
        s_up, p_up, s_dn, p_dn = 1, 2, 3, 4

        Ssp = cd[s_up] * c[s_dn]
        Ssm = cd[s_dn] * c[s_up]
        Spp = cd[p_up] * c[p_dn]
        Spm = cd[p_dn] * c[p_up]

        Ssx = 0.5 * (Ssp + Ssm)
        Ssy = -0.5im * (Ssp - Ssm)
        Spx = 0.5 * (Spp + Spm)
        Spy = -0.5im * (Spp - Spm)

        Tspinflip = Ssp * Spm + Ssm * Spp

        Dsame = n[s_up] * n[p_up] + n[s_dn] * n[p_dn]
        Dopp = n[s_up] * n[p_dn] + n[s_dn] * n[p_up]
        Ntot = n[s_up] + n[p_up] + n[s_dn] + n[p_dn]

        Fup = cd[s_up] * c[p_up] + cd[p_up] * c[s_up]
        Fdn = cd[s_dn] * c[p_dn] + cd[p_dn] * c[s_dn]
        Gup = -1im * (cd[s_up] * c[p_up] - cd[p_up] * c[s_up])
        Gdn = -1im * (cd[s_dn] * c[p_dn] - cd[p_dn] * c[s_dn])

        Of = Fup - Fdn
        Og = Gup - Gdn

        return (; c, cd, n, Ssx, Ssy, Spx, Spy, Tspinflip, Dsame, Dopp, Ntot, Of, Og)
    end

    ops = local_kanamori_operators()
    (; Tspinflip, Ssx, Ssy, Spx, Spy, Dsame, Dopp, Ntot, Of, Og) = ops

    @test norm(Tspinflip - 2 * (Ssx * Spx + Ssy * Spy)) < 1e-12
    @test norm(Tspinflip - (0.25 * (Of * Of + Og * Og) - 0.5 * Ntot + Dsame)) < 1e-12
    @test norm((Of * Of) * (Og * Og) - (Og * Og) * (Of * Of)) < 1e-12
    @test norm(Dsame * (Of * Of) - (Of * Of) * Dsame) < 1e-12
    @test norm(Dopp * (Og * Og) - (Og * Og) * Dopp) < 1e-12

    U = 1.0
    JH = 0.25
    Hdd = (U - 3JH) * Dsame + (U - 2JH) * Dopp
    Hsf = -JH * Tspinflip
    H = Hdd + Hsf

    # One electron per orbital: |s↑p↑>, |s↓p↓>, |s↑p↓>, |s↓p↑>.
    one_per_orbital_states = [
        (1 << 0) | (1 << 1),
        (1 << 2) | (1 << 3),
        (1 << 0) | (1 << 3),
        (1 << 2) | (1 << 1),
    ] .+ 1
    Hblock = H[one_per_orbital_states, one_per_orbital_states]
    evals = sort(real.(eigvals(Hermitian(Hblock))))

    @test evals ≈ sort([U - 3JH, U - 3JH, U - 3JH, U - JH]) atol=1e-12
end

@testitem "Exact discrete HS identity for spin-flip Hund squares" begin

    using LinearAlgebra

    const DIM_HUND_HS = 16

    function hs_local_annihilation(flavor::Int)
        c = zeros(ComplexF64, DIM_HUND_HS, DIM_HUND_HS)
        mask_below = (1 << (flavor - 1)) - 1
        for state in 0:(DIM_HUND_HS - 1)
            if ((state >> (flavor - 1)) & 1) == 1
                occupied_below = count_ones(state & mask_below)
                sign = isodd(occupied_below) ? -1.0 : 1.0
                new_state = state & ~(1 << (flavor - 1))
                c[new_state + 1, state + 1] = sign
            end
        end
        return c
    end

    function hs_local_operators()
        c = [hs_local_annihilation(f) for f in 1:4]
        cd = adjoint.(c)
        n = [cd[f] * c[f] for f in 1:4]
        s_up, p_up, s_dn, p_dn = 1, 2, 3, 4

        Ssp = cd[s_up] * c[s_dn]
        Ssm = cd[s_dn] * c[s_up]
        Spp = cd[p_up] * c[p_dn]
        Spm = cd[p_dn] * c[p_up]
        Tspinflip = Ssp * Spm + Ssm * Spp

        Dsame = n[s_up] * n[p_up] + n[s_dn] * n[p_dn]
        Dopp = n[s_up] * n[p_dn] + n[s_dn] * n[p_up]
        Ntot = n[s_up] + n[p_up] + n[s_dn] + n[p_dn]

        Fup = cd[s_up] * c[p_up] + cd[p_up] * c[s_up]
        Fdn = cd[s_dn] * c[p_dn] + cd[p_dn] * c[s_dn]
        Gup = -1im * (cd[s_up] * c[p_up] - cd[p_up] * c[s_up])
        Gdn = -1im * (cd[s_dn] * c[p_dn] - cd[p_dn] * c[s_dn])

        return (; Tspinflip, Dsame, Dopp, Ntot, Of=Fup-Fdn, Og=Gup-Gdn)
    end

    function exact_three_state_square_exp(O::AbstractMatrix, κ::Real)
        I16 = Matrix{ComplexF64}(I, size(O, 1), size(O, 2))
        if iszero(κ)
            return I16, 0.0, 1.0, 0.0
        end

        # O_f and O_g have eigenvalues in {-2,-1,0,1,2}; match x=0,1,2 exactly:
        # exp(κ x^2) = A + B cosh(λ x).
        y = (exp(4κ) - 1) / (2 * (exp(κ) - 1)) - 1
        λ = acosh(y)
        B = (exp(κ) - 1) / (cosh(λ) - 1)
        A = 1 - B
        M = A * I16 + (B / 2) * (exp(λ * O) + exp(-λ * O))
        return M, λ, A, B
    end

    ops = hs_local_operators()
    (; Tspinflip, Dsame, Dopp, Ntot, Of, Og) = ops

    for (U, JH, Δτ) in ((1.0, 0.25, 0.10), (2.0, 0.40, 0.05), (1.0, 0.10, 0.20))
        κ = Δτ * JH / 4
        Ef, λf, Af, Bf = exact_three_state_square_exp(Of, κ)
        Eg, λg, Ag, Bg = exact_three_state_square_exp(Og, κ)

        @test λf ≈ λg atol=1e-14
        @test Af > 0
        @test Bf > 0
        @test norm(Ef - exp(κ * (Of * Of))) < 1e-11
        @test norm(Eg - exp(κ * (Og * Og))) < 1e-11

        local Hdd = (U - 3JH) * Dsame + (U - 2JH) * Dopp
        local Hsf = -JH * Tspinflip

        # The square fields represent H_square = -JH/4*(Of^2+Og^2)
        # = H_sf - JH/2*Ntot + JH*Dsame.  Therefore combine them with
        # Hdd_adj = Hdd - JH*Dsame + JH/2*Ntot to recover Hdd + Hsf.
        local Hdd_adj = Hdd - JH * Dsame + (JH / 2) * Ntot

        local left = exp(-Δτ * (Hdd + Hsf))
        local right = exp(-Δτ * Hdd_adj) * Ef * Eg

        @test norm(left - right) < 1e-10
    end
end

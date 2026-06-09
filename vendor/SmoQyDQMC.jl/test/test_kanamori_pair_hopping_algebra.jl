@testitem "Kanamori pair-hopping algebra and HS identity" begin

    using LinearAlgebra

    const DIM_PAIR_HOP = 16

    function ph_annihilation(flavor::Int)
        c = zeros(ComplexF64, DIM_PAIR_HOP, DIM_PAIR_HOP)
        mask_below = (1 << (flavor - 1)) - 1
        for state in 0:(DIM_PAIR_HOP - 1)
            if ((state >> (flavor - 1)) & 1) == 1
                occupied_below = count_ones(state & mask_below)
                sign = isodd(occupied_below) ? -1.0 : 1.0
                new_state = state & ~(1 << (flavor - 1))
                c[new_state + 1, state + 1] = sign
            end
        end
        return c
    end

    function ph_square_exp(O::AbstractMatrix, κ::Real)
        I16 = Matrix{ComplexF64}(I, size(O,1), size(O,2))
        if iszero(κ)
            return I16
        end
        y = (exp(4κ) - 1) / (2 * (exp(κ) - 1)) - 1
        λ = acosh(complex(y, 0.0))
        B = real((exp(κ) - 1) / (cosh(λ) - 1))
        A = 1 - B
        @test A > 0
        @test B > 0
        return A * I16 + (B / 2) * (exp(λ * O) + exp(-λ * O))
    end

    c = [ph_annihilation(f) for f in 1:4]
    cd = adjoint.(c)
    n = [cd[f] * c[f] for f in 1:4]

    # Flavor order: 1=s↑, 2=p↑, 3=s↓, 4=p↓.
    s_up, p_up, s_dn, p_dn = 1, 2, 3, 4
    Fup = cd[s_up] * c[p_up] + cd[p_up] * c[s_up]
    Fdn = cd[s_dn] * c[p_dn] + cd[p_dn] * c[s_dn]
    Gup = -1im * (cd[s_up] * c[p_up] - cd[p_up] * c[s_up])
    Gdn = -1im * (cd[s_dn] * c[p_dn] - cd[p_dn] * c[s_dn])

    Oplus = Fup + Fdn
    Ominus = Fup - Fdn
    Gplus = Gup + Gdn

    Ssp = cd[s_up] * c[s_dn]
    Ssm = cd[s_dn] * c[s_up]
    Spp = cd[p_up] * c[p_dn]
    Spm = cd[p_dn] * c[p_up]
    Tspinflip = Ssp * Spm + Ssm * Spp
    Ppair = cd[s_up] * cd[s_dn] * c[p_dn] * c[p_up] + cd[p_up] * cd[p_dn] * c[s_dn] * c[s_up]

    Dsame = n[s_up] * n[p_up] + n[s_dn] * n[p_dn]
    Dopp = n[s_up] * n[p_dn] + n[s_dn] * n[p_up]

    @test norm(Ppair - 0.25 * (Oplus * Oplus - Gplus * Gplus)) < 1e-12
    @test norm((-Tspinflip + Ppair) - 0.25 * (Oplus * Oplus - Ominus * Ominus)) < 1e-12
    @test norm((Oplus * Oplus) * (Ominus * Ominus) - (Ominus * Ominus) * (Oplus * Oplus)) < 1e-12
    @test norm(Dsame * (Oplus * Oplus) - (Oplus * Oplus) * Dsame) < 1e-12
    @test norm(Dopp * (Ominus * Ominus) - (Ominus * Ominus) * Dopp) < 1e-12

    for (U, JH, Δτ) in ((1.0, 0.25, 0.10), (2.0, 0.40, 0.05), (1.0, 0.10, 0.20))
        local Hdd = (U - 3JH) * Dsame + (U - 2JH) * Dopp
        local Htrans = JH * (-Tspinflip + Ppair)
        Eplus = ph_square_exp(Oplus, -Δτ * JH / 4)
        Eminus = ph_square_exp(Ominus, +Δτ * JH / 4)
        left = exp(-Δτ * (Hdd + Htrans))
        right = exp(-Δτ * Hdd) * Eplus * Eminus
        @test norm(left - right) < 1e-10
    end
end

using LinearAlgebra
using Test

include(joinpath(@__DIR__, "..", "src", "PiFluxQSHKanamoriED.jl"))
using .PiFluxQSHKanamoriED

@testset "Pi-flux QSH Kanamori ED" begin
    @testset "one-body Hermiticity and TR partner" begin
        p = EDParams(Lx=3, Ly=4, t=1.0, lambda=0.2, U=0.0, JH=0.0, beta=1.0, mu=0.0)
        Kup = build_piflux_qsh_onebody(p, :up)
        Kdn = build_piflux_qsh_onebody(p, :dn)
        @test Kup ≈ Kup'
        @test Kdn ≈ Kdn'
        @test Kdn ≈ conj.(Kup)
    end

    @testset "noninteracting grand canonical regression" begin
        p = EDParams(Lx=2, Ly=2, t=1.0, lambda=0.15, U=0.0, JH=0.0,
                     beta=0.8, mu=0.1, interaction_level=:full)
        ed = grand_canonical_ed(p)
        free = free_fermion_grand_canonical(p)
        @test ed.logZ ≈ free.logZ atol=1e-11 rtol=1e-11
        @test ed.energy ≈ free.energy atol=1e-11 rtol=1e-11
        @test ed.N ≈ free.N atol=1e-11 rtol=1e-11
    end

    @testset "particle-hole half filling" begin
        for lvl in (:hubbard, :density, :spinflip, :full)
            p = EDParams(Lx=1, Ly=1, t=0.0, lambda=0.0, U=2.0, JH=0.3,
                         beta=1.2, mu=0.0, interaction_level=lvl)
            r = grand_canonical_ed(p)
            @test r.N ≈ 2.0 atol=1e-11
            @test r.Nup ≈ 1.0 atol=1e-11
            @test r.Ndn ≈ 1.0 atol=1e-11
        end
    end

    @testset "atomic factorization" begin
        one = EDParams(Lx=1, Ly=1, t=0.0, lambda=0.0, U=1.7, JH=0.25,
                       beta=1.1, mu=0.2, interaction_level=:full)
        two = EDParams(Lx=2, Ly=1, t=0.0, lambda=0.0, U=1.7, JH=0.25,
                       beta=1.1, mu=0.2, interaction_level=:full)
        r1 = grand_canonical_ed(one)
        r2 = grand_canonical_ed(two)
        @test r2.logZ ≈ 2 * r1.logZ atol=1e-11 rtol=1e-11
        @test r2.energy ≈ 2 * r1.energy atol=1e-11 rtol=1e-11
        @test r2.N ≈ 2 * r1.N atol=1e-11 rtol=1e-11
    end

    @testset "interaction hierarchy at nonzero JH" begin
        base = (; Lx=1, Ly=1, t=0.0, lambda=0.0, U=2.0, JH=0.3, beta=1.0, mu=0.1)
        rh = grand_canonical_ed(EDParams(; base..., interaction_level=:hubbard))
        rd = grand_canonical_ed(EDParams(; base..., interaction_level=:density))
        rs = grand_canonical_ed(EDParams(; base..., interaction_level=:spinflip))
        rf = grand_canonical_ed(EDParams(; base..., interaction_level=:full))
        @test abs(rh.logZ - rd.logZ) > 1e-6
        @test abs(rd.logZ - rs.logZ) > 1e-6
        @test abs(rs.logZ - rf.logZ) > 1e-6
    end

    @testset "JH zero transverse regression" begin
        base = (; Lx=1, Ly=1, t=0.2, lambda=0.05, U=1.2, JH=0.0, beta=0.9, mu=0.0)
        rd = grand_canonical_ed(EDParams(; base..., interaction_level=:density))
        rs = grand_canonical_ed(EDParams(; base..., interaction_level=:spinflip))
        rf = grand_canonical_ed(EDParams(; base..., interaction_level=:full))
        @test rd.logZ ≈ rs.logZ atol=1e-11 rtol=1e-11
        @test rd.logZ ≈ rf.logZ atol=1e-11 rtol=1e-11
        @test rd.energy ≈ rs.energy atol=1e-11 rtol=1e-11
        @test rd.energy ≈ rf.energy atol=1e-11 rtol=1e-11
    end
end

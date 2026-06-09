@testitem "Kanamori density model defaults and shifts" begin

    using SmoQyDQMC
    using Random
    import SmoQyDQMC.LatticeUtilities as lu

    rng = MersenneTwister(1)
    unit_cell = lu.UnitCell(
        lattice_vecs = [[1.0]],
        basis_vecs = [[0.0], [0.0]],
    )
    lattice = lu.Lattice(L = [2], periodic = [true])
    model_geometry = ModelGeometry(unit_cell, lattice)

    model = KanamoriDensityModel(ph_sym_form=false, U=1.0, JH=0.25)
    params = KanamoriDensityParameters(
        kanamori_density_model = model,
        model_geometry = model_geometry,
        rng = rng,
    )

    @test params.orbital_pairs == [(1, 2)]
    @test params.V_opposite == fill(0.5, 2)
    @test params.V_same == fill(0.25, 2)
    @test params.neighbor_table == [1 3; 2 4]

    hst = KanamoriDensityHirschHST(
        kanamori_density_parameters = params,
        β = 0.2,
        Δτ = 0.1,
        rng = rng,
    )
    @test hst.N == 2
    @test hst.Lτ == 2
    @test eltype(hst.α_same) <: Real
    @test eltype(hst.α_opposite) <: Real

    tbm = TightBindingModel(
        model_geometry = model_geometry,
        μ = 0.0,
        ϵ_mean = [0.0, 0.0],
        ϵ_std = [0.0, 0.0],
    )
    tbp = TightBindingParameters(
        tight_binding_model = tbm,
        model_geometry = model_geometry,
        rng = rng,
    )
    fpi_up = FermionPathIntegral(tight_binding_parameters = tbp, β = 0.2, Δτ = 0.1)
    fpi_dn = FermionPathIntegral(tight_binding_parameters = tbp, β = 0.2, Δτ = 0.1)
    initialize!(fpi_up, fpi_dn, params)

    @test fpi_up.V == fill(0.375, 4, 2)
    @test fpi_dn.V == fill(0.375, 4, 2)
end

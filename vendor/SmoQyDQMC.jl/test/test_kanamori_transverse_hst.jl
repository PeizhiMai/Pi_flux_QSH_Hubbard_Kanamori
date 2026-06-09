@testitem "KanamoriTransverseHST construction and smoke update" begin

    using LinearAlgebra
    using Random
    using SmoQyDQMC
    import SmoQyDQMC.LatticeUtilities as lu
    import SmoQyDQMC.JDQMCFramework as dqmcf

    rng = MersenneTwister(9)
    unit_cell = lu.UnitCell(lattice_vecs = [[1.0]], basis_vecs = [[0.0], [0.0]])
    lattice = lu.Lattice(L = [2], periodic = [true])
    model_geometry = ModelGeometry(unit_cell, lattice)

    tbm = TightBindingModel(
        model_geometry = model_geometry,
        μ = 0.0,
        ϵ_mean = [0.0, 0.0],
        ϵ_std = [0.0, 0.0],
    )
    tbp = TightBindingParameters(tight_binding_model = tbm, model_geometry = model_geometry, rng = rng)
    fpi_up = FermionPathIntegral(tight_binding_parameters = tbp, β = 0.2, Δτ = 0.1, forced_complex_kinetic = true)
    fpi_dn = FermionPathIntegral(tight_binding_parameters = tbp, β = 0.2, Δτ = 0.1, forced_complex_kinetic = true)

    hst = KanamoriTransverseHST(model_geometry = model_geometry, JH = 0.25, β = 0.2, Δτ = 0.1, rng = rng)
    initialize!(fpi_up, fpi_dn, hst)

    @test fpi_up.V == zeros(4, 2)
    @test fpi_dn.V == zeros(4, 2)
    @test fpi_up.Sb == fpi_dn.Sb
    @test any(!iszero, imag.(hst.λ_plus))
    @test all(iszero, imag.(hst.λ_minus))

    Bup = initialize_propagators(fpi_up, symmetric = false, checkerboard = false)
    Bdn = initialize_propagators(fpi_dn, symmetric = false, checkerboard = false)
    apply_kanamori_transverse_to_propagators!(Bup, Bdn, hst)

    for B in (Bup..., Bdn...)
        @test B.expmΔτK * B.exppΔτK ≈ Matrix{ComplexF64}(I, 4, 4) atol=1e-11
    end

    fgc_up = dqmcf.FermionGreensCalculator(Bup, 0.2, 0.1, 1)
    fgc_dn = dqmcf.FermionGreensCalculator(Bdn, 0.2, 0.1, 1)
    Gup = zeros(ComplexF64, 4, 4)
    Gdn = zeros(ComplexF64, 4, 4)
    logdetGup, sgndetGup = dqmcf.calculate_equaltime_greens!(Gup, fgc_up)
    logdetGdn, sgndetGdn = dqmcf.calculate_equaltime_greens!(Gdn, fgc_dn)

    acc, logdetGup, sgndetGup, logdetGdn, sgndetGdn, δG, δθ = local_updates!(
        Gup, logdetGup, sgndetGup,
        Gdn, logdetGdn, sgndetGdn,
        hst;
        fermion_path_integral_up = fpi_up,
        fermion_path_integral_dn = fpi_dn,
        fermion_greens_calculator_up = fgc_up,
        fermion_greens_calculator_dn = fgc_dn,
        Bup = Bup,
        Bdn = Bdn,
        δG = 0.0,
        δθ = 0.0,
        rng = rng,
        δG_max = 1e-6,
        update_stabilization_frequency = false,
    )

    @test 0.0 <= acc <= 1.0
    @test isfinite(logdetGup)
    @test isfinite(logdetGdn)
    @test δG >= 0.0
    @test δθ >= 0.0
end

#!/usr/bin/env julia

using LinearAlgebra
using Random
using Printf
using TOML
using Serialization

using SmoQyDQMC
import SmoQyDQMC.LatticeUtilities as lu
import SmoQyDQMC.JDQMCFramework as dqmcf

const DEFAULTS = Dict{String,Any}(
    "Lx" => 4,
    "Ly" => 4,
    "t" => 1.0,
    "lambda" => 0.2,
    "U" => 1.0,
    "JH" => 0.25,
    "mu" => 0.5,
    "open_x" => false,
    "density_kanamori" => false,
    "spin_flip_hund" => false,
    "pair_hopping" => false,
    "beta" => 2.0,
    "dtau" => 0.05,
    "Ntherm" => 200,
    "Nmeas" => 1000,
    "Nupdates" => 2,
    "n_stab" => 10,
    "dGmax" => 1e-6,
    "symmetric" => false,
    "checkerboard" => false,
    "seed" => 0,
    "sID" => 1,
    "outdir" => joinpath(@__DIR__, "..", "results", "piflux_qsh_smoqy"),
    "checkpoint_enable" => false,
    "checkpoint_file" => "",
    "checkpoint_freq_hours" => 0.5,
    "runtime_limit_hours" => 0.0,
    "checkpoint_keep" => true,
    "checkpoint_resume" => true,
)

function parse_value(x::AbstractString, default)
    if default isa Bool
        lx = lowercase(x)
        lx in ("true","t","1","yes","y") && return true
        lx in ("false","f","0","no","n") && return false
        error("Cannot parse boolean value: $x")
    elseif default isa Int
        return parse(Int, x)
    elseif default isa AbstractFloat
        return parse(Float64, x)
    else
        return x
    end
end

function parse_args(defaults::Dict{String,Any}, args)
    p = copy(defaults)
    i = 1
    while i <= length(args)
        a = args[i]
        startswith(a, "--") || error("Expected --key=value or --key value, got $a")
        s = a[3:end]
        if occursin("=", s)
            key, val = split(s, "=", limit=2)
        else
            key = s; i += 1; i <= length(args) || error("Missing value for --$key"); val = args[i]
        end
        haskey(defaults, key) || error("Unknown option --$key. Known keys: $(sort(collect(keys(defaults))))")
        p[key] = parse_value(val, defaults[key])
        i += 1
    end
    return p
end

function half_filling_mu(U::Real, JH::Real; density_kanamori::Bool)
    return density_kanamori ? (3U - 5JH) / 2 : U / 2
end

@inline spin_sign(spin::Symbol) = spin === :up ? 1.0 : -1.0
@inline sm_t(h) = -conj(h) # SmoQy convention: target h c_i† c_j -> t=-conj(h)

function piflux_qsh_bonds_and_tmeans(; t::Float64, lambda::Float64, spin::Symbol)
    sgn = spin_sign(spin)
    A, B = 1, 2
    bonds = lu.Bond{2}[]
    tmean = ComplexF64[]
    function push_hop!(orb1, orb2, disp, h)
        push!(bonds, lu.Bond(orbitals=(orb1, orb2), displacement=disp))
        push!(tmean, sm_t(h))
    end
    push_hop!(A, B, [0, 0], +t)
    push_hop!(B, A, [1, 0], +t)
    push_hop!(A, A, [0, 1], +t)
    push_hop!(B, B, [0, 1], -t)
    push_hop!(A, B, [0, 1], -im * lambda * sgn)
    push_hop!(A, B, [0,-1], -im * lambda * sgn)
    push_hop!(B, A, [1, 1], -im * lambda * sgn)
    push_hop!(B, A, [1,-1], +im * lambda * sgn)
    return bonds, tmean
end

function remove_x_boundary_wrap_hoppings(tbp, model_geometry)
    unit_cell = model_geometry.unit_cell
    lattice = model_geometry.lattice
    Lx = lattice.L[1]
    keep = trues(size(tbp.neighbor_table, 2))
    for n in axes(tbp.neighbor_table, 2)
        i = tbp.neighbor_table[1, n]
        j = tbp.neighbor_table[2, n]
        loci, _ = lu.site_to_loc(i, unit_cell, lattice)
        locj, _ = lu.site_to_loc(j, unit_cell, lattice)
        if loci[1] == Lx - 1 && locj[1] == 0
            keep[n] = false
        end
    end
    kept = findall(keep)
    new_neighbor_table = tbp.neighbor_table[:, kept]
    new_t = tbp.t[kept]
    new_bond_ids = Int[]
    new_bond_slices = UnitRange{Int}[]
    cursor = 1
    for (bond_id, bond_slice) in zip(tbp.bond_ids, tbp.bond_slices)
        nkept = count(keep[n] for n in bond_slice)
        if nkept > 0
            push!(new_bond_ids, bond_id)
            push!(new_bond_slices, cursor:(cursor+nkept-1))
            cursor += nkept
        end
    end
    return TightBindingParameters(tbp.μ, copy(tbp.ϵ), new_t, new_neighbor_table,
                                  new_bond_ids, new_bond_slices, tbp.norbital)
end

function initialize_piflux_qsh(; Lx::Int, Ly::Int, t::Float64, lambda::Float64,
                               U::Float64, JH::Float64, mu::Float64, open_x::Bool,
                               density_kanamori::Bool, spin_flip_hund::Bool, pair_hopping::Bool,
                               beta::Float64, dtau::Float64, rng::AbstractRNG)
    unit_cell = lu.UnitCell(
        lattice_vecs = [[1.0, 0.0], [0.0, 1.0]],
        basis_vecs = [[0.0, 0.0], [0.0, 0.0]],
    )
    lattice = lu.Lattice(L=[Lx, Ly], periodic=[true, true])
    model_geometry = ModelGeometry(unit_cell, lattice)

    bonds_up, t_up = piflux_qsh_bonds_and_tmeans(t=t, lambda=lambda, spin=:up)
    bonds_dn, t_dn = piflux_qsh_bonds_and_tmeans(t=t, lambda=lambda, spin=:dn)
    tbm_up = TightBindingModel(model_geometry=model_geometry, μ=mu, ϵ_mean=[0.0,0.0], ϵ_std=[0.0,0.0],
                               t_bonds=bonds_up, t_mean=t_up, t_std=zeros(Float64, length(t_up)))
    tbm_dn = TightBindingModel(model_geometry=model_geometry, μ=mu, ϵ_mean=[0.0,0.0], ϵ_std=[0.0,0.0],
                               t_bonds=bonds_dn, t_mean=t_dn, t_std=zeros(Float64, length(t_dn)))
    tbp_up = TightBindingParameters(tight_binding_model=tbm_up, model_geometry=model_geometry, rng=rng)
    tbp_dn = TightBindingParameters(tight_binding_model=tbm_dn, model_geometry=model_geometry, rng=rng)
    if open_x
        tbp_up = remove_x_boundary_wrap_hoppings(tbp_up, model_geometry)
        tbp_dn = remove_x_boundary_wrap_hoppings(tbp_dn, model_geometry)
    end

    hubbard_model = HubbardModel(ph_sym_form=false, U_orbital=[1,2], U_mean=[U,U], U_std=[0.0,0.0])
    hubbard_parameters = HubbardParameters(hubbard_model=hubbard_model, model_geometry=model_geometry, rng=rng)
    hubbard_hst = HubbardSpinHirschHST(β=beta, Δτ=dtau, hubbard_parameters=hubbard_parameters, rng=rng)

    if density_kanamori
        Vsame = (spin_flip_hund && !pair_hopping) ? U - 4JH : U - 3JH
        kanamori_density_model = KanamoriDensityModel(ph_sym_form=false, orbital_pairs=[(1,2)], U=U, JH=JH,
                                                       V_same=[Vsame], V_opposite=[U - 2JH])
        kanamori_density_parameters = KanamoriDensityParameters(kanamori_density_model=kanamori_density_model,
                                                                 model_geometry=model_geometry, rng=rng)
        kanamori_density_hst = KanamoriDensityHirschHST(kanamori_density_parameters=kanamori_density_parameters,
                                                        β=beta, Δτ=dtau, rng=rng)
        if spin_flip_hund && pair_hopping
            transverse_hst = KanamoriTransverseHST(model_geometry=model_geometry, orbital_pairs=[(1,2)],
                                                   JH=JH, β=beta, Δτ=dtau, rng=rng)
            hst_parameters = (hubbard_hst, kanamori_density_hst, transverse_hst)
        elseif spin_flip_hund
            transverse_hst = HundSpinFlipHST(model_geometry=model_geometry, orbital_pairs=[(1,2)],
                                             JH=JH, β=beta, Δτ=dtau, rng=rng)
            hst_parameters = (hubbard_hst, kanamori_density_hst, transverse_hst)
        else
            transverse_hst = nothing
            hst_parameters = (hubbard_hst, kanamori_density_hst)
        end
    else
        kanamori_density_parameters = nothing
        transverse_hst = nothing
        hst_parameters = hubbard_hst
    end

    fpi_up = FermionPathIntegral(tight_binding_parameters=tbp_up, β=beta, Δτ=dtau,
                                  forced_complex_kinetic=(lambda != 0.0 || spin_flip_hund))
    fpi_dn = FermionPathIntegral(tight_binding_parameters=tbp_dn, β=beta, Δτ=dtau,
                                  forced_complex_kinetic=(lambda != 0.0 || spin_flip_hund))
    initialize!(fpi_up, fpi_dn, hubbard_parameters)
    density_kanamori && initialize!(fpi_up, fpi_dn, kanamori_density_parameters)
    initialize!(fpi_up, fpi_dn, hst_parameters)

    return model_geometry, tbp_up, tbp_dn, hubbard_parameters, kanamori_density_parameters,
           transverse_hst, hst_parameters, fpi_up, fpi_dn
end

function acceptance_scalar(acc)
    acc isa Tuple && return sum(acc) / length(acc)
    return acc
end

function density_profile(Gup, Gdn, model_geometry)
    unit_cell = model_geometry.unit_cell
    lattice = model_geometry.lattice
    Lx, Ly = lattice.L[1], lattice.L[2]
    prof = zeros(ComplexF64, Lx, 4) # A_up, A_dn, B_up, B_dn
    for site in 1:size(Gup, 1)
        loc, orb = lu.site_to_loc(site, unit_cell, lattice)
        rx = loc[1] + 1
        col_up = (orb - 1) * 2 + 1
        col_dn = (orb - 1) * 2 + 2
        prof[rx, col_up] += (1 - Gup[site, site]) / Ly
        prof[rx, col_dn] += (1 - Gdn[site, site]) / Ly
    end
    return prof
end

function write_profile(filename, prof)
    open(filename, "w") do io
        println(io, "# x\tn_A_up\tn_A_dn\tn_B_up\tn_B_dn\tn_total\tm_z")
        for x in axes(prof, 1)
            nAu, nAd, nBu, nBd = real(prof[x,1]), real(prof[x,2]), real(prof[x,3]), real(prof[x,4])
            ntot = nAu + nAd + nBu + nBd
            mz = 0.5 * ((nAu+nBu) - (nAd+nBd))
            @printf(io, "%d\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\n", x-1, nAu, nAd, nBu, nBd, ntot, mz)
        end
    end
end

function interaction_energy_per_cell(hubbard_parameters, kanamori_density_parameters, Gup, Gdn)
    e = zero(eltype(Gup))
    e += measure_hubbard_energy(hubbard_parameters, Gup, Gdn, 1)
    e += measure_hubbard_energy(hubbard_parameters, Gup, Gdn, 2)
    if kanamori_density_parameters !== nothing
        e += measure_kanamori_density_energy(kanamori_density_parameters, Gup, Gdn, 1)
    end
    return e
end

function write_checkpoint_file(filename::AbstractString, state::Dict{String,Any})
    mkpath(dirname(filename))
    tmp = filename * ".tmp"
    open(tmp, "w") do io
        serialize(io, state)
    end
    mv(tmp, filename; force=true)
    open(filename * ".status", "w") do io
        @printf(io, "written_at_epoch = %.6f\n", time())
        @printf(io, "thermal_done = %d\n", get(state, "thermal_done", -1))
        @printf(io, "meas_done = %d\n", get(state, "meas_done", -1))
    end
    return time()
end

function read_checkpoint_file(filename::AbstractString)
    open(filename, "r") do io
        return deserialize(io)
    end
end

function reconstruct_greens(Bup, Bdn, β::Float64, Δτ::Float64, n_stab::Int)
    fgc_up = dqmcf.FermionGreensCalculator(Bup, β, Δτ, n_stab)
    fgc_dn = dqmcf.FermionGreensCalculator(Bdn, β, Δτ, n_stab)
    Gup = zeros(eltype(Bup[1]), size(Bup[1]))
    Gdn = zeros(eltype(Bdn[1]), size(Bdn[1]))
    logdetGup, sgndetGup = dqmcf.calculate_equaltime_greens!(Gup, fgc_up)
    logdetGdn, sgndetGdn = dqmcf.calculate_equaltime_greens!(Gdn, fgc_dn)
    return fgc_up, fgc_dn, Gup, Gdn, logdetGup, sgndetGup, logdetGdn, sgndetGdn
end

function main()
    p = parse_args(DEFAULTS, ARGS)
    p["U"] == 0.0 && error("This driver is for finite U. Use check_piflux_qsh_h0.jl for U=0.")
    mu_half = half_filling_mu(p["U"], p["JH"]; density_kanamori=p["density_kanamori"])
    if abs(p["mu"] - mu_half) > 1e-10
        @warn "mu differs from the physical half-filling value for this interaction tier" mu=p["mu"] mu_half
    end
    iseven(p["Ly"]) || @warn "Odd Ly with periodic y can break the particle-hole parity used by the sign-free parent."
    p["spin_flip_hund"] && !p["density_kanamori"] && error("spin_flip_hund=true requires density_kanamori=true")
    p["pair_hopping"] && !p["spin_flip_hund"] && error("pair_hopping=true requires spin_flip_hund=true")
    p["spin_flip_hund"] && (p["symmetric"] || p["checkerboard"]) && error("transverse Hund currently requires asymmetric exact propagators")

    β = p["beta"]
    Δτ = p["dtau"]
    Lτ = dqmcf.eval_length_imaginary_axis(β, Δτ)
    seed = p["seed"] == 0 ? abs(rand(Int)) : p["seed"]
    rng = Xoshiro(seed)
    prefix = @sprintf("piflux_qsh_U%.3f_JH%.3f_Lx%d_Ly%d_b%.3f_dt%.3f_lam%.3f_mu%.3f_openx%s_den%s_sf%s_pair%s",
                      p["U"], p["JH"], p["Lx"], p["Ly"], β, Δτ, p["lambda"], p["mu"],
                      string(p["open_x"]), string(p["density_kanamori"]), string(p["spin_flip_hund"]), string(p["pair_hopping"]))
    outdir = joinpath(p["outdir"], prefix * @sprintf("-%03d", p["sID"]))
    mkpath(outdir)
    checkpoint_file = isempty(p["checkpoint_file"]) ? joinpath(outdir, "checkpoint.jls") : p["checkpoint_file"]

    metadata = Dict{String,Any}(string(k)=>v for (k,v) in p)
    metadata["seed_actual"] = seed
    metadata["Ltau"] = Lτ
    metadata["mu_half_filling"] = mu_half
    metadata["checkpoint_file_resolved"] = checkpoint_file
    metadata["interaction_convention"] = "physical_non_ph_shifted"
    metadata["interaction_scope"] = (p["spin_flip_hund"] && p["pair_hopping"]) ?
        "physical Hubbard plus full local Kanamori transverse exchange" :
        (p["spin_flip_hund"] ? "physical Hubbard plus density/Ising Hund and spin-flip Hund" :
         (p["density_kanamori"] ? "physical Hubbard plus density/Ising Kanamori" : "physical Hubbard only"))
    open(joinpath(outdir, "metadata.toml"), "w") do io
        TOML.print(io, metadata)
    end

    @info "Initializing pi-flux QSH SmoQyDQMC" outdir β Δτ Lτ seed checkpoint_file

    start_timestamp = time()
    last_checkpoint_timestamp = start_timestamp

    thermal_done = 0
    meas_done = 0
    δG = 0.0
    δθ = 0.0
    acc_sum = 0.0
    nsweeps = 0
    prof_num = zeros(ComplexF64, p["Lx"], 4)
    phase_sum = 0.0 + 0.0im
    phase_abs_sum = 0.0
    eint_num = 0.0 + 0.0im
    measurement_rows = Tuple{Int,Float64,Float64,Float64,Float64,Float64,Float64,Float64,Float64,Float64,Float64}[]

    if p["checkpoint_enable"] && p["checkpoint_resume"] && isfile(checkpoint_file)
        @info "Resuming pi-flux QSH SmoQyDQMC from checkpoint" checkpoint_file
        checkpoint = read_checkpoint_file(checkpoint_file)
        checkpoint_parameters = get(checkpoint, "parameters", Dict{String,Any}())
        for key in ("Lx", "Ly", "t", "lambda", "U", "JH", "mu", "open_x", "density_kanamori", "spin_flip_hund", "pair_hopping", "beta", "dtau")
            if haskey(checkpoint_parameters, key) && checkpoint_parameters[key] != p[key]
                error("Checkpoint parameter mismatch for $key: checkpoint=$(checkpoint_parameters[key]) requested=$(p[key])")
            end
        end

        rng = checkpoint["rng"]
        model_geometry = checkpoint["model_geometry"]
        tbp_up = checkpoint["tbp_up"]
        tbp_dn = checkpoint["tbp_dn"]
        hubbard_parameters = get(checkpoint, "hubbard_parameters", nothing)
        kanamori_density_parameters = get(checkpoint, "kanamori_density_parameters", nothing)
        hst_parameters = checkpoint["hst_parameters"]
        fpi_up = checkpoint["fpi_up"]
        fpi_dn = checkpoint["fpi_dn"]
        Bup = checkpoint["Bup"]
        Bdn = checkpoint["Bdn"]
        thermal_done = checkpoint["thermal_done"]
        meas_done = checkpoint["meas_done"]
        δG = checkpoint["δG"]
        δθ = checkpoint["δθ"]
        acc_sum = checkpoint["acc_sum"]
        nsweeps = checkpoint["nsweeps"]
        prof_num = checkpoint["prof_num"]
        phase_sum = checkpoint["phase_sum"]
        phase_abs_sum = checkpoint["phase_abs_sum"]
        eint_num = get(checkpoint, "eint_num", 0.0 + 0.0im)
        measurement_rows = checkpoint["measurement_rows"]
        last_checkpoint_timestamp = time()
        fgc_up, fgc_dn, Gup, Gdn, logdetGup, sgndetGup, logdetGdn, sgndetGdn = reconstruct_greens(Bup, Bdn, β, Δτ, p["n_stab"])
    else
        model_geometry, tbp_up, tbp_dn, hubbard_parameters, kanamori_density_parameters,
        transverse_hst, hst_parameters, fpi_up, fpi_dn = initialize_piflux_qsh(
            Lx=p["Lx"], Ly=p["Ly"], t=p["t"], lambda=p["lambda"], U=p["U"], JH=p["JH"],
            mu=p["mu"], open_x=p["open_x"], density_kanamori=p["density_kanamori"],
            spin_flip_hund=p["spin_flip_hund"], pair_hopping=p["pair_hopping"], beta=β, dtau=Δτ, rng=rng)

        Bup = initialize_propagators(fpi_up, symmetric=p["symmetric"], checkerboard=p["checkerboard"])
        Bdn = initialize_propagators(fpi_dn, symmetric=p["symmetric"], checkerboard=p["checkerboard"])
        if p["spin_flip_hund"] && p["pair_hopping"]
            apply_kanamori_transverse_to_propagators!(Bup, Bdn, transverse_hst)
        elseif p["spin_flip_hund"]
            apply_hund_spinflip_to_propagators!(Bup, Bdn, transverse_hst)
        end
        fgc_up, fgc_dn, Gup, Gdn, logdetGup, sgndetGup, logdetGdn, sgndetGdn = reconstruct_greens(Bup, Bdn, β, Δτ, p["n_stab"])
    end

    function checkpoint_state(current_thermal_done::Int, current_meas_done::Int)
        return Dict{String,Any}(
            "driver" => "run_piflux_qsh_smoqy.jl",
            "format_version" => 1,
            "parameters" => Dict{String,Any}(string(k) => v for (k, v) in p),
            "rng" => rng,
            "model_geometry" => model_geometry,
            "tbp_up" => tbp_up,
            "tbp_dn" => tbp_dn,
            "hubbard_parameters" => hubbard_parameters,
            "kanamori_density_parameters" => kanamori_density_parameters,
            "hst_parameters" => hst_parameters,
            "fpi_up" => fpi_up,
            "fpi_dn" => fpi_dn,
            "Bup" => Bup,
            "Bdn" => Bdn,
            "thermal_done" => current_thermal_done,
            "meas_done" => current_meas_done,
            "δG" => δG,
            "δθ" => δθ,
            "acc_sum" => acc_sum,
            "nsweeps" => nsweeps,
            "prof_num" => prof_num,
            "phase_sum" => phase_sum,
            "phase_abs_sum" => phase_abs_sum,
            "eint_num" => eint_num,
            "measurement_rows" => measurement_rows,
            "written_at_epoch" => time(),
        )
    end

    function maybe_checkpoint!(current_thermal_done::Int, current_meas_done::Int; force::Bool=false)
        p["checkpoint_enable"] || return false
        now = time()
        checkpoint_due = p["checkpoint_freq_hours"] > 0.0 && (now - last_checkpoint_timestamp) >= p["checkpoint_freq_hours"] * 3600.0
        runtime_due = p["runtime_limit_hours"] > 0.0 && (now - start_timestamp) >= p["runtime_limit_hours"] * 3600.0
        if force || checkpoint_due || runtime_due
            last_checkpoint_timestamp = write_checkpoint_file(checkpoint_file, checkpoint_state(current_thermal_done, current_meas_done))
            @info "Wrote pi-flux QSH SmoQyDQMC checkpoint" checkpoint_file current_thermal_done current_meas_done runtime_due
        end
        return runtime_due
    end

    maybe_checkpoint!(thermal_done, meas_done; force=(p["checkpoint_enable"] && thermal_done == 0 && meas_done == 0 && !isfile(checkpoint_file)))

    for sweep in (thermal_done + 1):p["Ntherm"]
        acc, logdetGup, sgndetGup, logdetGdn, sgndetGdn, δG, δθ = local_updates!(
            Gup, logdetGup, sgndetGup, Gdn, logdetGdn, sgndetGdn, hst_parameters;
            fermion_path_integral_up=fpi_up, fermion_path_integral_dn=fpi_dn,
            fermion_greens_calculator_up=fgc_up, fermion_greens_calculator_dn=fgc_dn,
            Bup=Bup, Bdn=Bdn, δG_max=p["dGmax"], δG=δG, δθ=δθ, rng=rng,
            update_stabilization_frequency=false)
        acc_sum += acceptance_scalar(acc); nsweeps += 1
        thermal_done = sweep
        maybe_checkpoint!(thermal_done, meas_done) && exit(13)
    end

    for meas in (meas_done + 1):p["Nmeas"]
        for _ in 1:p["Nupdates"]
            acc, logdetGup, sgndetGup, logdetGdn, sgndetGdn, δG, δθ = local_updates!(
                Gup, logdetGup, sgndetGup, Gdn, logdetGdn, sgndetGdn, hst_parameters;
                fermion_path_integral_up=fpi_up, fermion_path_integral_dn=fpi_dn,
                fermion_greens_calculator_up=fgc_up, fermion_greens_calculator_dn=fgc_dn,
                Bup=Bup, Bdn=Bdn, δG_max=p["dGmax"], δG=δG, δθ=δθ, rng=rng,
                update_stabilization_frequency=false)
            acc_sum += acceptance_scalar(acc); nsweeps += 1
        end
        phase = sgndetGup * sgndetGdn
        if abs(phase) != 0
            phase = conj(phase / abs(phase))
        end
        prof = density_profile(Gup, Gdn, model_geometry)
        eint = interaction_energy_per_cell(hubbard_parameters, kanamori_density_parameters, Gup, Gdn)
        prof_num .+= phase .* prof
        eint_num += phase * eint
        phase_sum += phase
        phase_abs_sum += abs(phase)
        ntot_complex = sum(prof) / p["Lx"]
        phase_ntot = phase * ntot_complex
        phase_eint = phase * eint
        push!(measurement_rows, (meas, real(phase), imag(phase), real(ntot_complex), imag(ntot_complex),
                                 real(eint), imag(eint), acc_sum/nsweeps, δG, real(phase_ntot), real(phase_eint)))
        meas_done = meas
        maybe_checkpoint!(thermal_done, meas_done) && exit(13)
    end

    avg_phase = phase_sum / p["Nmeas"]
    prof_avg = prof_num ./ phase_sum
    eint_avg = eint_num / phase_sum
    density_cell = real(sum(prof_avg) / p["Lx"])
    write_profile(joinpath(outdir, "density_profile.tsv"), prof_avg)
    open(joinpath(outdir, "measurements.tsv"), "w") do io
        println(io, "# meas\tphase_re\tphase_im\tn_per_cell_unreweighted_re\tn_per_cell_unreweighted_im\teint_per_cell_unreweighted_re\teint_per_cell_unreweighted_im\tacceptance_running\tdG\tphase_n_re\tphase_eint_re")
        for row in measurement_rows
            @printf(io, "%d\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\n", row...)
        end
    end
    open(joinpath(outdir, "summary.txt"), "w") do io
        @printf(io, "average_phase = %.12g%+.12gi\n", real(avg_phase), imag(avg_phase))
        @printf(io, "average_phase_abs = %.12g\n", abs(avg_phase))
        @printf(io, "average_abs_phase_per_measurement = %.12g\n", phase_abs_sum / p["Nmeas"])
        @printf(io, "acceptance_rate = %.12g\n", acc_sum / nsweeps)
        @printf(io, "dG = %.12g\n", δG)
        @printf(io, "dtheta = %.12g\n", δθ)
        @printf(io, "density_per_cell_reweighted = %.12g\n", density_cell)
        @printf(io, "interaction_energy_per_cell_reweighted = %.12g%+.12gi\n", real(eint_avg), imag(eint_avg))
        @printf(io, "interaction_scope = %s\n", metadata["interaction_scope"])
        if p["density_kanamori"]
            @printf(io, "V_opposite = U - 2JH = %.12g\n", p["U"] - 2*p["JH"])
            if p["spin_flip_hund"] && p["pair_hopping"]
                @printf(io, "V_same = U - 3JH = %.12g\n", p["U"] - 3*p["JH"])
                @printf(io, "pair_hopping = true; Jpair = JH = %.12g\n", p["JH"])
            elseif p["spin_flip_hund"]
                @printf(io, "V_same_adjusted = U - 4JH = %.12g\n", p["U"] - 4*p["JH"])
                @printf(io, "spin_flip_hund = true; onsite compensation = JH/2 = %.12g\n", p["JH"]/2)
            else
                @printf(io, "V_same = U - 3JH = %.12g\n", p["U"] - 3*p["JH"])
            end
        end
        @printf(io, "output_dir = %s\n", outdir)
    end

    if p["checkpoint_enable"]
        if p["checkpoint_keep"]
            write_checkpoint_file(checkpoint_file, checkpoint_state(thermal_done, meas_done))
        else
            isfile(checkpoint_file) && rm(checkpoint_file; force=true)
            isfile(checkpoint_file * ".status") && rm(checkpoint_file * ".status"; force=true)
        end
    end

    open(joinpath(outdir, "complete.status"), "w") do io
        @printf(io, "completed_at_epoch = %.6f\n", time())
        @printf(io, "measurements = %d\n", meas_done)
    end
    @printf("DQMC complete: %s\n", outdir)
    @printf("  average phase = %.8g%+.8gi\n", real(avg_phase), imag(avg_phase))
    @printf("  |average phase| = %.8g\n", abs(avg_phase))
    @printf("  density/cell  = %.8g\n", density_cell)
    @printf("  interaction E = %.8g%+.8gi\n", real(eint_avg), imag(eint_avg))
end

main()

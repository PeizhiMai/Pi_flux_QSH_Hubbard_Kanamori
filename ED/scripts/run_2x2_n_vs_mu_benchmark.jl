#!/usr/bin/env julia

include(joinpath(@__DIR__, "..", "src", "PiFluxQSHKanamoriED.jl"))
using .PiFluxQSHKanamoriED
using LinearAlgebra
using Printf
using TOML

const DEFAULTS = Dict{String,Any}(
    "Lx" => 2,
    "Ly" => 2,
    "t" => 1.0,
    "lambda" => 0.2,
    "U" => 1.0,
    "JH" => 0.25,
    "beta" => 7.0,
    "mu_min" => -1.5,
    "mu_max" => 2.5,
    "mu_step" => 0.05,
    "levels" => "hubbard,density,spinflip,full",
    "bcs" => "pbc,cylinder",
    "outdir" => joinpath(@__DIR__, "..", "results", "ed_2x2_beta7_n_vs_mu"),
)

function parse_bool(s::AbstractString)
    x = lowercase(s)
    x in ("true", "t", "1", "yes", "y") && return true
    x in ("false", "f", "0", "no", "n") && return false
    error("Cannot parse Bool: $s")
end

function parse_value(s::AbstractString, default)
    default isa Bool && return parse_bool(s)
    default isa Int && return parse(Int, s)
    default isa AbstractFloat && return parse(Float64, s)
    return s
end

function parse_args(args=ARGS)
    p = copy(DEFAULTS)
    i = 1
    while i <= length(args)
        a = args[i]
        startswith(a, "--") || error("Expected --key=value or --key value, got $a")
        s = a[3:end]
        if occursin("=", s)
            key, val = split(s, "=", limit=2)
        else
            key = s
            i += 1
            i <= length(args) || error("Missing value for --$key")
            val = args[i]
        end
        haskey(p, key) || error("Unknown option --$key. Known keys: $(sort(collect(keys(p))))")
        p[key] = parse_value(val, p[key])
        i += 1
    end
    return p
end

split_symbols(s::AbstractString) = Symbol.(strip.(split(s, ","; keepempty=false)))
split_strings(s::AbstractString) = strip.(split(s, ","; keepempty=false))

function mu_grid(mu_min::Float64, mu_max::Float64, mu_step::Float64, required::Vector{Float64})
    mu_step > 0 || error("mu_step must be positive")
    n = floor(Int, (mu_max - mu_min) / mu_step + 1e-9)
    grid = [mu_min + i * mu_step for i in 0:n if mu_min + i * mu_step <= mu_max + 1e-9]
    append!(grid, [μ for μ in required if mu_min - 1e-9 <= μ <= mu_max + 1e-9])
    return sort!(unique!(round.(grid; digits=12)))
end

function spectrum_for_params(p::EDParams)
    PiFluxQSHKanamoriED.validate_params(p)
    Ns = PiFluxQSHKanamoriED.nsites(p)
    rows = NamedTuple[]
    sectors = PiFluxQSHKanamoriED.SectorInfo[]
    total_dim = 0
    tdiag = @elapsed begin
        for Nup in 0:Ns, Ndn in 0:Ns
            H, _states = PiFluxQSHKanamoriED.build_sector_hamiltonian(p, Nup, Ndn)
            dim = size(H, 1)
            dim == 0 && continue
            evals = eigvals(Hermitian(Matrix(H)))
            evals = real.(evals)
            push!(sectors, PiFluxQSHKanamoriED.SectorInfo(Nup, Ndn, dim, minimum(evals), maximum(evals)))
            total_dim += dim
            for E in evals
                push!(rows, (E=Float64(E), Nup=Float64(Nup), Ndn=Float64(Ndn), N=Float64(Nup + Ndn)))
            end
        end
    end
    return rows, sectors, total_dim, tdiag
end

function evaluate_spectrum(rows, p::EDParams, μ::Float64)
    β = p.beta
    qmin = minimum(r.E - μ * r.N for r in rows)
    Zs = 0.0
    Es = 0.0
    Ns = 0.0
    Nups = 0.0
    Ndns = 0.0
    N2s = 0.0
    for r in rows
        q = r.E - μ * r.N
        w = exp(-β * (q - qmin))
        Zs += w
        Es += w * r.E
        Ns += w * r.N
        Nups += w * r.Nup
        Ndns += w * r.Ndn
        N2s += w * r.N^2
    end
    logZ = log(Zs) - β * qmin
    Eavg = Es / Zs
    Navg = Ns / Zs
    Nupavg = Nups / Zs
    Ndnavg = Ndns / Zs
    N2avg = N2s / Zs
    cells = p.Lx * p.Ly
    density = Navg / cells
    kappa_cell = β * (N2avg - Navg^2) / cells
    return (; logZ, energy=Eavg, N=Navg, Nup=Nupavg, Ndn=Ndnavg,
            density_per_cell=density, kappa_cell)
end

function bc_to_open_x(bc::AbstractString)
    b = lowercase(bc)
    b == "pbc" && return false
    b in ("cyl", "cylinder", "cylindrical", "open_x") && return true
    error("Unknown boundary condition '$bc'; use pbc or cylinder")
end

function main(args=ARGS)
    cfg = parse_args(args)
    levels = split_symbols(cfg["levels"])
    bcs = split_strings(cfg["bcs"])
    required_mu = [half_filling_mu(cfg["U"], cfg["JH"], lvl) for lvl in levels]
    μs = mu_grid(cfg["mu_min"], cfg["mu_max"], cfg["mu_step"], required_mu)
    outdir = cfg["outdir"]
    mkpath(outdir)
    curve_path = joinpath(outdir, "n_vs_mu.tsv")
    summary_path = joinpath(outdir, "summary.tsv")
    sector_path = joinpath(outdir, "sector_summary.tsv")
    metadata_path = joinpath(outdir, "metadata.toml")

    open(metadata_path, "w") do io
        TOML.print(io, Dict(
            "Lx" => cfg["Lx"], "Ly" => cfg["Ly"], "t" => cfg["t"], "lambda" => cfg["lambda"],
            "U" => cfg["U"], "JH" => cfg["JH"], "beta" => cfg["beta"],
            "mu_min" => cfg["mu_min"], "mu_max" => cfg["mu_max"], "mu_step" => cfg["mu_step"],
            "levels" => String.(levels), "bcs" => bcs,
            "interaction_convention" => "physical_non_ph_shifted",
            "half_filling_mu_hubbard" => half_filling_mu(cfg["U"], cfg["JH"], :hubbard),
            "half_filling_mu_kanamori" => half_filling_mu(cfg["U"], cfg["JH"], :density),
        ))
    end

    open(curve_path, "w") do curve_io
        println(curve_io, "bc\topen_x\tlevel\tbeta\tLx\tLy\tt\tlambda\tU\tJH\tmu\tmu_half\tN\tNup\tNdn\tdensity_per_cell\tenergy\tlogZ\tkappa_cell")
        open(summary_path, "w") do sum_io
            println(sum_io, "bc\topen_x\tlevel\tbeta\tLx\tLy\tU\tJH\tmu_half\tdensity_at_mu_half\tN_at_mu_half\tNup_at_mu_half\tNdn_at_mu_half\thalf_density_error\tkappa_cell_at_mu_half\tenergy_at_mu_half\tdiag_seconds\ttotal_dim")
            open(sector_path, "w") do sec_io
                println(sec_io, "bc\topen_x\tlevel\tNup\tNdn\tdim\tEmin\tEmax")
                for lvl in levels
                    for bc in bcs
                        open_x = bc_to_open_x(bc)
                        p = EDParams(Lx=cfg["Lx"], Ly=cfg["Ly"], t=cfg["t"], lambda=cfg["lambda"],
                                     U=cfg["U"], JH=cfg["JH"], beta=cfg["beta"], mu=0.0,
                                     open_x=open_x, interaction_level=lvl)
                        μhalf = half_filling_mu(p)
                        @printf("[ED 2x2] diagonalizing level=%s bc=%s beta=%.6g mu_half=%.12g\n", String(lvl), bc, p.beta, μhalf)
                        flush(stdout)
                        rows, sectors, total_dim, tdiag = spectrum_for_params(p)
                        @printf("[ED 2x2] done level=%s bc=%s dim=%d diag_seconds=%.3f\n", String(lvl), bc, total_dim, tdiag)
                        flush(stdout)
                        for s in sectors
                            @printf(sec_io, "%s\t%s\t%s\t%d\t%d\t%d\t%.12g\t%.12g\n", bc, string(open_x), String(lvl), s.Nup, s.Ndn, s.dim, s.Emin, s.Emax)
                        end
                        half_obs = nothing
                        for μ in μs
                            obs = evaluate_spectrum(rows, p, μ)
                            @printf(curve_io, "%s\t%s\t%s\t%.12g\t%d\t%d\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\n",
                                    bc, string(open_x), String(lvl), p.beta, p.Lx, p.Ly, p.t, p.lambda, p.U, p.JH, μ, μhalf,
                                    obs.N, obs.Nup, obs.Ndn, obs.density_per_cell, obs.energy, obs.logZ, obs.kappa_cell)
                            abs(μ - μhalf) < 1e-10 && (half_obs = obs)
                        end
                        half_obs === nothing && error("mu_half=$μhalf was not included in the grid")
                        half_error = half_obs.density_per_cell - 2.0
                        @printf(sum_io, "%s\t%s\t%s\t%.12g\t%d\t%d\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\t%d\n",
                                bc, string(open_x), String(lvl), p.beta, p.Lx, p.Ly, p.U, p.JH, μhalf,
                                half_obs.density_per_cell, half_obs.N, half_obs.Nup, half_obs.Ndn, half_error,
                                half_obs.kappa_cell, half_obs.energy, tdiag, total_dim)
                        flush(curve_io); flush(sum_io); flush(sec_io)
                    end
                end
            end
        end
    end
    println("wrote $curve_path")
    println("wrote $summary_path")
    println("wrote $sector_path")
    return nothing
end

main()

module PiFluxQSHKanamoriED

using LinearAlgebra
using Printf
using SparseArrays
using TOML

export EDParams, EDResult, SectorInfo, build_piflux_qsh_onebody,
       grand_canonical_ed, free_fermion_grand_canonical, half_filling_mu, write_outputs,
       parse_cli_args, main

Base.@kwdef struct EDParams
    Lx::Int = 2
    Ly::Int = 2
    t::Float64 = 1.0
    lambda::Float64 = 0.2
    U::Float64 = 1.0
    JH::Float64 = 0.25
    beta::Float64 = 2.0
    mu::Float64 = 0.0
    open_x::Bool = false
    interaction_level::Symbol = :full # :hubbard, :density, :spinflip, :full
end

struct SectorInfo
    Nup::Int
    Ndn::Int
    dim::Int
    Emin::Float64
    Emax::Float64
end

struct EDResult
    params::EDParams
    logZ::Float64
    energy::Float64
    N::Float64
    Nup::Float64
    Ndn::Float64
    density_per_cell::Float64
    profile::Matrix{Float64} # Lx × 4: A_up, A_dn, B_up, B_dn averaged over y
    sectors::Vector{SectorInfo}
    qmin::Float64
end

struct FermionTerm
    coeff::ComplexF64
    ops::Vector{Tuple{Bool,Int}} # (is_creation, flavor), applied right-to-left
end

const _VALID_LEVELS = (:hubbard, :density, :spinflip, :full)

function half_filling_mu(U::Real, JH::Real, interaction_level::Symbol)
    interaction_level == :hubbard && return U / 2
    interaction_level in _VALID_LEVELS || error("interaction_level must be one of $(_VALID_LEVELS), got $(interaction_level)")
    return (3U - 5JH) / 2
end

half_filling_mu(p::EDParams) = half_filling_mu(p.U, p.JH, p.interaction_level)


@inline nsites(p::EDParams) = 2 * p.Lx * p.Ly
@inline nflavors(p::EDParams) = 2 * nsites(p)
@inline up_flavor(site::Int, Ns::Int) = site
@inline dn_flavor(site::Int, Ns::Int) = Ns + site
@inline spin_sign(spin::Symbol) = spin === :up ? 1.0 : -1.0

function validate_params(p::EDParams)
    p.Lx > 0 || error("Lx must be positive")
    p.Ly > 0 || error("Ly must be positive")
    p.beta > 0 || error("beta must be positive")
    p.U >= 0 || error("U must be non-negative")
    p.JH >= 0 || error("JH must be non-negative")
    p.interaction_level in _VALID_LEVELS || error("interaction_level must be one of $(_VALID_LEVELS), got $(p.interaction_level)")
    nflavors(p) <= 62 || error("Bit ED supports at most 62 flavors")
    return p
end

@inline function site_index(x0::Int, y0::Int, orb::Int, p::EDParams)
    return 2 * (x0 * p.Ly + y0) + orb
end

@inline function site_coordinates(site::Int, p::EDParams)
    q = site - 1
    orb = (q % 2) + 1
    cell = q ÷ 2
    y0 = cell % p.Ly
    x0 = cell ÷ p.Ly
    return x0, y0, orb
end

function add_hop!(K::AbstractMatrix{ComplexF64}, p::EDParams, x0::Int, y0::Int,
                  orb1::Int, orb2::Int, dx::Int, dy::Int, h::ComplexF64)
    x1 = x0 + dx
    y1 = y0 + dy
    if p.open_x && (x1 < 0 || x1 >= p.Lx)
        return nothing
    end
    x1 = mod(x1, p.Lx)
    y1 = mod(y1, p.Ly)
    i = site_index(x0, y0, orb1, p)
    j = site_index(x1, y1, orb2, p)
    # Match SmoQy/JDQMC finite-lattice assignment semantics for duplicate
    # periodic images on tiny clusters: assign rather than accumulate.
    K[i,j] = h
    K[j,i] = conj(h)
    return nothing
end

function piflux_qsh_bond_defs(t::Float64, lambda::Float64, spin::Symbol)
    s = spin_sign(spin)
    return (
        (1, 2, 0,  0, ComplexF64(+t)),
        (2, 1, 1,  0, ComplexF64(+t)),
        (1, 1, 0,  1, ComplexF64(+t)),
        (2, 2, 0,  1, ComplexF64(-t)),
        (1, 2, 0,  1, ComplexF64(-im * lambda * s)),
        (1, 2, 0, -1, ComplexF64(-im * lambda * s)),
        (2, 1, 1,  1, ComplexF64(-im * lambda * s)),
        (2, 1, 1, -1, ComplexF64(+im * lambda * s)),
    )
end

function build_piflux_qsh_onebody(p::EDParams, spin::Symbol)
    validate_params(p)
    spin in (:up, :dn) || error("spin must be :up or :dn")
    Ns = nsites(p)
    K = zeros(ComplexF64, Ns, Ns)
    for (orb1, orb2, dx, dy, h) in piflux_qsh_bond_defs(p.t, p.lambda, spin)
        for y0 in 0:p.Ly-1, x0 in 0:p.Lx-1
            add_hop!(K, p, x0, y0, orb1, orb2, dx, dy, h)
        end
    end
    return K
end

function spin_masks(Ns::Int, N::Int)
    masks = UInt64[]
    maxmask = UInt64(1) << Ns
    for m in UInt64(0):(maxmask - UInt64(1))
        count_ones(m) == N && push!(masks, m)
    end
    return masks
end

function sector_basis(Ns::Int, Nup::Int, Ndn::Int)
    ups = spin_masks(Ns, Nup)
    dns = spin_masks(Ns, Ndn)
    states = Vector{UInt64}(undef, length(ups) * length(dns))
    n = 0
    for d in dns, u in ups
        n += 1
        states[n] = u | (d << Ns)
    end
    index = Dict{UInt64,Int}(s => i for (i, s) in pairs(states))
    return states, index
end

@inline occ(state::UInt64, flavor::Int) = !iszero(state & (UInt64(1) << (flavor - 1)))

@inline function apply_single_op(state::UInt64, is_creation::Bool, flavor::Int)
    bit = UInt64(1) << (flavor - 1)
    occupied = !iszero(state & bit)
    if is_creation
        occupied && return nothing
        lower = bit - UInt64(1)
        sgn = isodd(count_ones(state & lower)) ? -1.0 : 1.0
        return (state | bit, sgn)
    else
        occupied || return nothing
        lower = bit - UInt64(1)
        sgn = isodd(count_ones(state & lower)) ? -1.0 : 1.0
        return (state & ~bit, sgn)
    end
end

function apply_term(state::UInt64, term::FermionTerm)
    amp = term.coeff
    st = state
    for (is_creation, flavor) in Iterators.reverse(term.ops)
        out = apply_single_op(st, is_creation, flavor)
        out === nothing && return nothing
        st, sgn = out
        amp *= sgn
    end
    return st, amp
end

function collect_onebody_terms(Kup::AbstractMatrix, Kdn::AbstractMatrix, Ns::Int; atol::Float64=1e-14)
    terms = FermionTerm[]
    for j in 1:Ns, i in 1:Ns
        hij = ComplexF64(Kup[i,j])
        abs(hij) > atol && push!(terms, FermionTerm(hij, [(true, up_flavor(i,Ns)), (false, up_flavor(j,Ns))]))
        hdij = ComplexF64(Kdn[i,j])
        abs(hdij) > atol && push!(terms, FermionTerm(hdij, [(true, dn_flavor(i,Ns)), (false, dn_flavor(j,Ns))]))
    end
    return terms
end

function collect_transverse_terms(p::EDParams)
    Ns = nsites(p)
    terms = FermionTerm[]
    p.interaction_level in (:hubbard, :density) && return terms
    for x0 in 0:p.Lx-1, y0 in 0:p.Ly-1
        a = site_index(x0, y0, 1, p)
        b = site_index(x0, y0, 2, p)
        au = up_flavor(a, Ns); bu = up_flavor(b, Ns)
        ad = dn_flavor(a, Ns); bd = dn_flavor(b, Ns)
        # -JH (S_A^+ S_B^- + S_A^- S_B^+)
        push!(terms, FermionTerm(ComplexF64(-p.JH), [(true, au), (false, ad), (true, bd), (false, bu)]))
        push!(terms, FermionTerm(ComplexF64(-p.JH), [(true, ad), (false, au), (true, bu), (false, bd)]))
        if p.interaction_level == :full
            # +JH (P_A† P_B + P_B† P_A)
            push!(terms, FermionTerm(ComplexF64(+p.JH), [(true, au), (true, ad), (false, bd), (false, bu)]))
            push!(terms, FermionTerm(ComplexF64(+p.JH), [(true, bu), (true, bd), (false, ad), (false, au)]))
        end
    end
    return terms
end

function diagonal_interaction_energy(state::UInt64, p::EDParams)
    Ns = nsites(p)
    e = 0.0
    for x0 in 0:p.Lx-1, y0 in 0:p.Ly-1
        a = site_index(x0, y0, 1, p)
        b = site_index(x0, y0, 2, p)
        nau = occ(state, up_flavor(a,Ns)) ? 1.0 : 0.0
        nad = occ(state, dn_flavor(a,Ns)) ? 1.0 : 0.0
        nbu = occ(state, up_flavor(b,Ns)) ? 1.0 : 0.0
        nbd = occ(state, dn_flavor(b,Ns)) ? 1.0 : 0.0
        e += p.U * (nau * nad + nbu * nbd)
        if p.interaction_level != :hubbard
            Vsame = p.U - 3p.JH
            Vopp = p.U - 2p.JH
            e += Vsame * (nau * nbu + nad * nbd)
            e += Vopp * (nau * nbd + nad * nbu)
        end
    end
    return e
end

function build_sector_hamiltonian(p::EDParams, Nup::Int, Ndn::Int)
    validate_params(p)
    Ns = nsites(p)
    states, index = sector_basis(Ns, Nup, Ndn)
    dim = length(states)
    rows = Int[]; cols = Int[]; vals = ComplexF64[]
    Kup = build_piflux_qsh_onebody(p, :up)
    Kdn = build_piflux_qsh_onebody(p, :dn)
    terms = vcat(collect_onebody_terms(Kup, Kdn, Ns), collect_transverse_terms(p))

    for (col, st) in pairs(states)
        ediag = diagonal_interaction_energy(st, p)
        if ediag != 0.0
            push!(rows, col); push!(cols, col); push!(vals, ComplexF64(ediag))
        end
        for term in terms
            out = apply_term(st, term)
            out === nothing && continue
            st2, amp = out
            row = get(index, st2, 0)
            row == 0 && continue
            push!(rows, row); push!(cols, col); push!(vals, amp)
        end
    end
    return sparse(rows, cols, vals, dim, dim), states
end

function state_profile(state::UInt64, p::EDParams)
    Ns = nsites(p)
    prof = zeros(Float64, p.Lx, 4)
    for site in 1:Ns
        x0, _, orb = site_coordinates(site, p)
        col_up = (orb - 1) * 2 + 1
        col_dn = (orb - 1) * 2 + 2
        prof[x0+1, col_up] += (occ(state, up_flavor(site,Ns)) ? 1.0 : 0.0) / p.Ly
        prof[x0+1, col_dn] += (occ(state, dn_flavor(site,Ns)) ? 1.0 : 0.0) / p.Ly
    end
    return prof
end

function free_fermion_grand_canonical(p::EDParams)
    validate_params(p)
    Ns = nsites(p)
    Kup = build_piflux_qsh_onebody(p, :up)
    Kdn = build_piflux_qsh_onebody(p, :dn)
    blocks = ((Kup, :up), (Kdn, :dn))
    logZ = 0.0; energy = 0.0; Ntot = 0.0; Nup = 0.0; Ndn = 0.0
    prof = zeros(Float64, p.Lx, 4)
    for (K, spin) in blocks
        evals, evecs = eigen(Hermitian(K))
        for n in eachindex(evals)
            ξ = evals[n] - p.mu
            logZ += log1p(exp(-p.beta * ξ))
            f = 1 / (1 + exp(p.beta * ξ))
            energy += evals[n] * f
            Ntot += f
            spin === :up ? (Nup += f) : (Ndn += f)
            for site in 1:Ns
                x0, _, orb = site_coordinates(site, p)
                col = (orb - 1) * 2 + (spin === :up ? 1 : 2)
                prof[x0+1, col] += abs2(evecs[site,n]) * f / p.Ly
            end
        end
    end
    return EDResult(p, logZ, energy, Ntot, Nup, Ndn, Ntot/(p.Lx*p.Ly), prof, SectorInfo[], 0.0)
end

function grand_canonical_ed(p::EDParams; write_spectrum::Bool=false)
    validate_params(p)
    p.U == 0 && p.JH == 0 && return free_fermion_grand_canonical(p)
    Ns = nsites(p)
    sectors = SectorInfo[]
    entries = Tuple{Float64,Float64,Float64,Float64,Matrix{Float64}}[]
    qmin = Inf
    for Nup in 0:Ns, Ndn in 0:Ns
        H, states = build_sector_hamiltonian(p, Nup, Ndn)
        dim = length(states)
        if dim == 0; continue; end
        F = eigen(Hermitian(Matrix(H)))
        evals = real.(F.values)
        push!(sectors, SectorInfo(Nup, Ndn, dim, minimum(evals), maximum(evals)))
        for a in eachindex(evals)
            E = evals[a]
            q = E - p.mu * (Nup + Ndn)
            qmin = min(qmin, q)
            prof = zeros(Float64, p.Lx, 4)
            vec = F.vectors[:,a]
            for (idx, st) in pairs(states)
                w = abs2(vec[idx])
                w == 0 && continue
                prof .+= w .* state_profile(st, p)
            end
            push!(entries, (q, E, Float64(Nup), Float64(Ndn), prof))
        end
    end
    weights = [exp(-p.beta * (q - qmin)) for (q,_,_,_,_) in entries]
    Zs = sum(weights)
    logZ = log(Zs) - p.beta * qmin
    energy = 0.0; Nup_avg = 0.0; Ndn_avg = 0.0
    prof_avg = zeros(Float64, p.Lx, 4)
    for (w, (_, E, Nu, Nd, prof)) in zip(weights, entries)
        ww = w / Zs
        energy += ww * E
        Nup_avg += ww * Nu
        Ndn_avg += ww * Nd
        prof_avg .+= ww .* prof
    end
    Ntot = Nup_avg + Ndn_avg
    return EDResult(p, logZ, energy, Ntot, Nup_avg, Ndn_avg, Ntot/(p.Lx*p.Ly), prof_avg, sectors, qmin)
end

function write_outputs(result::EDResult, outdir::AbstractString)
    mkpath(outdir)
    p = result.params
    open(joinpath(outdir, "ed_summary.toml"), "w") do io
        param_value(v) = v isa Symbol ? string(v) : v
        TOML.print(io, Dict(
            "parameters" => Dict(string(k) => param_value(getfield(p,k)) for k in fieldnames(EDParams)),
            "observables" => Dict(
                "logZ" => result.logZ,
                "energy" => result.energy,
                "N" => result.N,
                "Nup" => result.Nup,
                "Ndn" => result.Ndn,
                "density_per_cell" => result.density_per_cell,
            )
        ))
    end
    open(joinpath(outdir, "ed_density_profile.tsv"), "w") do io
        println(io, "x\tA_up\tA_dn\tB_up\tB_dn")
        for x in 1:p.Lx
            @printf(io, "%d\t%.12g\t%.12g\t%.12g\t%.12g\n", x-1, result.profile[x,1], result.profile[x,2], result.profile[x,3], result.profile[x,4])
        end
    end
    open(joinpath(outdir, "ed_sector_summary.tsv"), "w") do io
        println(io, "Nup\tNdn\tdim\tEmin\tEmax")
        for s in result.sectors
            @printf(io, "%d\t%d\t%d\t%.12g\t%.12g\n", s.Nup, s.Ndn, s.dim, s.Emin, s.Emax)
        end
    end
    return nothing
end

function parse_bool(s::AbstractString)
    x = lowercase(s)
    x in ("true","t","1","yes","y") && return true
    x in ("false","f","0","no","n") && return false
    error("Cannot parse Bool: $s")
end

function parse_cli_args(args=ARGS)
    defaults = Dict{String,Any}(
        "Lx"=>2, "Ly"=>2, "t"=>1.0, "lambda"=>0.2, "U"=>1.0, "JH"=>0.25,
        "beta"=>2.0, "mu"=>0.0, "open_x"=>false, "interaction_level"=>"full",
        "outdir"=>joinpath(@__DIR__, "..", "results", "piflux_qsh_ed"),
    )
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
        haskey(defaults, key) || error("Unknown option --$key")
        d = defaults[key]
        p[key] = d isa Bool ? parse_bool(val) : d isa Int ? parse(Int,val) : d isa Float64 ? parse(Float64,val) : val
        i += 1
    end
    ep = EDParams(Lx=p["Lx"], Ly=p["Ly"], t=p["t"], lambda=p["lambda"], U=p["U"], JH=p["JH"],
                  beta=p["beta"], mu=p["mu"], open_x=p["open_x"], interaction_level=Symbol(p["interaction_level"]))
    return ep, p["outdir"]
end

function main(args=ARGS)
    p, outdir = parse_cli_args(args)
    res = grand_canonical_ed(p)
    write_outputs(res, outdir)
    @printf("ED complete: level=%s Lx=%d Ly=%d beta=%.6g mu=%.6g\n", p.interaction_level, p.Lx, p.Ly, p.beta, p.mu)
    @printf("  logZ = %.12g\n  energy = %.12g\n  N = %.12g\n  density_per_cell = %.12g\n", res.logZ, res.energy, res.N, res.density_per_cell)
    return res
end

end # module

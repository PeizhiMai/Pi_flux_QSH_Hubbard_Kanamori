#!/usr/bin/env julia

using LinearAlgebra
using Printf

const DEFAULTS = Dict{String,Any}(
    "Nx" => 6,
    "Ly" => 6,
    "t" => 1.0,
    "lambda" => 0.2,
    "grid" => 41,
    "open_x" => false,
    "print_matrix" => false,
)

function parse_value(x::AbstractString, default)
    if default isa Bool
        lx = lowercase(x)
        lx in ("true","t","1","yes","y") && return true
        lx in ("false","f","0","no","n") && return false
        error("Cannot parse Bool: $x")
    elseif default isa Int
        return parse(Int, x)
    elseif default isa AbstractFloat
        return parse(Float64, x)
    else
        return x
    end
end

function parse_args(defaults, args)
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
        p[key] = parse_value(val, defaults[key])
        i += 1
    end
    return p
end

@inline spin_sign(spin::Symbol) = spin === :up ? 1.0 : -1.0
@inline site_index(x0::Int, y0::Int, orb::Int, Nx::Int, Ly::Int) = 2 * (x0 * Ly + y0) + orb

function piflux_qsh_bond_defs(t::Float64, λ::Float64, spin::Symbol)
    s = spin_sign(spin)
    return (
        (1, 2, 0,  0, ComplexF64(+t)),
        (2, 1, 1,  0, ComplexF64(+t)),
        (1, 1, 0,  1, ComplexF64(+t)),
        (2, 2, 0,  1, ComplexF64(-t)),
        (1, 2, 0,  1, ComplexF64(-im * λ * s)),
        (1, 2, 0, -1, ComplexF64(-im * λ * s)),
        (2, 1, 1,  1, ComplexF64(-im * λ * s)),
        (2, 1, 1, -1, ComplexF64(+im * λ * s)),
    )
end

function add_hop!(K, Nx, Ly, open_x, x0, y0, orb1, orb2, dx, dy, h)
    x1 = x0 + dx; y1 = y0 + dy
    if open_x && (x1 < 0 || x1 >= Nx)
        return nothing
    end
    x1 = mod(x1, Nx); y1 = mod(y1, Ly)
    i = site_index(x0, y0, orb1, Nx, Ly)
    j = site_index(x1, y1, orb2, Nx, Ly)
    K[i,j] = h
    K[j,i] = conj(h)
    return nothing
end

function build_realspace(Nx, Ly, t, λ, spin; open_x=false)
    K = zeros(ComplexF64, 2Nx*Ly, 2Nx*Ly)
    for (orb1, orb2, dx, dy, h) in piflux_qsh_bond_defs(t, λ, spin)
        for y0 in 0:Ly-1, x0 in 0:Nx-1
            add_hop!(K, Nx, Ly, open_x, x0, y0, orb1, orb2, dx, dy, h)
        end
    end
    return K
end

function bloch_hamiltonian(kx, ky, t, λ, spin)
    H = zeros(ComplexF64, 2, 2)
    for (a,b,dx,dy,h) in piflux_qsh_bond_defs(t, λ, spin)
        phase = exp(im * (kx * dx + ky * dy))
        H[a,b] += h * phase
        H[b,a] += conj(h) * conj(phase)
    end
    return H
end

function chern_number(t, λ, spin; N=41)
    ks = range(-π, π; length=N+1)[1:end-1]
    vec = Array{ComplexF64}(undef, N, N, 2)
    min_gap = Inf
    for ix in 1:N, iy in 1:N
        H = bloch_hamiltonian(ks[ix], ks[iy], t, λ, spin)
        F = eigen(Hermitian(H))
        min_gap = min(min_gap, F.values[2] - F.values[1])
        vec[ix,iy,:] .= F.vectors[:,1]
    end
    flux = 0.0
    for ix in 1:N, iy in 1:N
        ixp = ix == N ? 1 : ix + 1
        iyp = iy == N ? 1 : iy + 1
        u = vec[ix,iy,:]
        ux = vec[ixp,iy,:]
        uy = vec[ix,iyp,:]
        uxy = vec[ixp,iyp,:]
        Ux = dot(u, ux); Ux /= abs(Ux)
        Uy = dot(u, uy); Uy /= abs(Uy)
        Ux_y = dot(uy, uxy); Ux_y /= abs(Ux_y)
        Uy_x = dot(ux, uxy); Uy_x /= abs(Uy_x)
        flux += angle(Ux * Uy_x / (Ux_y * Uy))
    end
    return flux / (2π), min_gap
end

function ph_parity(Nx, Ly)
    d = zeros(Float64, 2Nx*Ly)
    for x in 0:Nx-1, y in 0:Ly-1
        ηA = iseven(y) ? 1.0 : -1.0
        d[site_index(x,y,1,Nx,Ly)] = ηA
        d[site_index(x,y,2,Nx,Ly)] = -ηA
    end
    return Diagonal(d)
end

function main()
    p = parse_args(DEFAULTS, ARGS)
    Nx = p["Nx"]; Ly = p["Ly"]; t = p["t"]; λ = p["lambda"]
    Kup = build_realspace(Nx, Ly, t, λ, :up; open_x=p["open_x"])
    Kdn = build_realspace(Nx, Ly, t, λ, :dn; open_x=p["open_x"])
    herm = max(norm(Kup - Kup'), norm(Kdn - Kdn'))
    trerr = norm(Kdn - conj.(Kup))
    D = ph_parity(Nx, Ly)
    pherr = norm(D * Kdn * D + Kup)
    Cup, gap_up = chern_number(t, λ, :up; N=p["grid"])
    Cdn, gap_dn = chern_number(t, λ, :dn; N=p["grid"])
    p["print_matrix"] && (println("Kup ="); show(stdout, MIME("text/plain"), Kup); println())
    @printf("Pi-flux QSH H0 check: Nx=%d Ly=%d t=%.6g lambda=%.6g grid=%d\n", Nx, Ly, t, λ, p["grid"])
    @printf("  hermiticity_error = %.6e\n", herm)
    @printf("  time_reversal_error = %.6e\n", trerr)
    @printf("  particle_hole_error = %.6e\n", pherr)
    @printf("  chern_up = %.12g\n", Cup)
    @printf("  chern_dn = %.12g\n", Cdn)
    @printf("  min_gap_up = %.12g\n", gap_up)
    @printf("  min_gap_dn = %.12g\n", gap_dn)
    herm < 1e-10 || error("Hermiticity check failed")
    trerr < 1e-10 || error("Time-reversal block check failed")
    pherr < 1e-10 || error("Particle-hole check failed; use even Ly with periodic-y")
    abs(round(Cup) - Cup) < 1e-8 || error("C_up is not quantized")
    abs(round(Cdn) - Cdn) < 1e-8 || error("C_dn is not quantized")
    round(Int, Cup) == (λ >= 0 ? 1 : -1) || error("Unexpected C_up")
    round(Int, Cdn) == (λ >= 0 ? -1 : 1) || error("Unexpected C_dn")
    gap_up > 1e-6 && gap_dn > 1e-6 || error("Bulk gap too small")
end

main()

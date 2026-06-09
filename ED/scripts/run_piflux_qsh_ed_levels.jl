#!/usr/bin/env julia
include(joinpath(@__DIR__, "..", "src", "PiFluxQSHKanamoriED.jl"))
using .PiFluxQSHKanamoriED

levels = (:hubbard, :density, :spinflip, :full)
p, outdir = PiFluxQSHKanamoriED.parse_cli_args(ARGS)
mkpath(outdir)
for lvl in levels
    p2 = EDParams(; (name=>getfield(p,name) for name in fieldnames(EDParams) if name != :interaction_level)..., interaction_level=lvl)
    res = grand_canonical_ed(p2)
    lvl_out = joinpath(outdir, string(lvl))
    write_outputs(res, lvl_out)
    @info "ED level complete" lvl logZ=res.logZ energy=res.energy N=res.N
end

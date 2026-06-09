#!/usr/bin/env julia
include(joinpath(@__DIR__, "..", "src", "PiFluxQSHKanamoriED.jl"))
using .PiFluxQSHKanamoriED
PiFluxQSHKanamoriED.main(ARGS)

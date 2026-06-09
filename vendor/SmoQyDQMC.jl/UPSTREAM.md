# Upstream provenance

This directory vendors upstream SmoQyDQMC.jl as the starting point for the BHZ-Hubbard-Kanamori implementation.

- Upstream: https://github.com/SmoQySuite/SmoQyDQMC.jl
- Imported commit: `de6ebab`
- Commit date: 2026-04-29T15:37:21-04:00
- Upstream version in `Project.toml`: `2.0.11`
- License: MIT, retained in `LICENSE`

The first project-specific changes should be made in a small set of new files under `src/Kanamori/` and corresponding exports/tests, keeping the upstream code layout recognizable.

# 2×2 ED benchmark at β=7: `n(μ)` for physical interactions

This benchmark uses the physical (non-PH-shifted) interaction convention shared
with `BHZ_Hubbard_Kanamori`.  The run diagonalizes each 2×2 cluster exactly once
per tier and boundary condition, then evaluates the grand-canonical density over
a chemical-potential grid.

Command used:

```bash
julia --project=ED ED/scripts/run_2x2_n_vs_mu_benchmark.jl \
  --Lx=2 --Ly=2 --beta=7.0 --t=1.0 --lambda=0.2 \
  --U=1.0 --JH=0.25 \
  --mu_min=-1.5 --mu_max=2.5 --mu_step=0.05 \
  --levels=hubbard,density,spinflip,full \
  --bcs=pbc,cylinder \
  --outdir=results/ed_2x2_beta7_n_vs_mu
```

Full local output from this run is under
`/Users/cosdis/Desktop/projects/strongly_correlated_topological_systems/Pi_flux_QSH_Hubbard_Kanamori/results/ed_2x2_beta7_n_vs_mu/`.

## Half-filling check

For `U=1.0, JH=0.25`, the physical half-filling chemical potentials are

- Hubbard only: `mu_half = U/2 = 0.5`;
- Kanamori tiers: `mu_half = (3U - 5JH)/2 = 0.875`.

At those `mu` values, the 2×2 ED density is exactly half-filled within floating
point precision for both PBC and cylindrical/open-x boundary conditions.

| tier | BC | μ_half | density/cell at μ_half | N | N_up | N_dn | density error |
|---|---|---:|---:|---:|---:|---:|---:|
| Hubbard | PBC | 0.5 | 2.000000000000 | 8.000000000000 | 4.000000000000 | 4.000000000000 | -1.15e-14 |
| Hubbard | cylinder | 0.5 | 2.000000000000 | 8.000000000000 | 4.000000000000 | 4.000000000000 | -4.44e-16 |
| Density/Ising | PBC | 0.875 | 2.000000000000 | 8.000000000000 | 4.000000000000 | 4.000000000000 | -5.33e-15 |
| Density/Ising | cylinder | 0.875 | 2.000000000000 | 8.000000000000 | 4.000000000000 | 4.000000000000 | +4.44e-15 |
| Spin-flip | PBC | 0.875 | 2.000000000000 | 8.000000000000 | 4.000000000000 | 4.000000000000 | +7.11e-15 |
| Spin-flip | cylinder | 0.875 | 2.000000000000 | 8.000000000000 | 4.000000000000 | 4.000000000000 | +1.24e-14 |
| Full | PBC | 0.875 | 2.000000000000 | 8.000000000000 | 4.000000000000 | 4.000000000000 | -5.33e-15 |
| Full | cylinder | 0.875 | 2.000000000000 | 8.000000000000 | 4.000000000000 | 4.000000000000 | +2.22e-15 |

## Particle-hole density check around physical μ_half

The grid also satisfies the expected finite-cluster density relation
`n(μ_half + δ) + n(μ_half - δ) = 4` per unit cell to numerical precision.
The largest checked absolute error was approximately `1.0e-11`.

| tier | BC | max absolute density-sum error |
|---|---|---:|
| Hubbard | PBC | 4.0e-12 |
| Hubbard | cylinder | 5.0e-12 |
| Density/Ising | PBC | 1.0e-11 |
| Density/Ising | cylinder | 0 |
| Spin-flip | PBC | 0 |
| Spin-flip | cylinder | 0 |
| Full | PBC | 0 |
| Full | cylinder | 0 |

## Notes

- Each tier/BC diagonalization covers the full `2^16 = 65536` Hilbert space.
- The curves were generated from cached exact spectra, not by re-diagonalizing at each `mu`.
- This ED benchmark confirms the physical-interaction half-filling convention to use before β=10 DQMC sign scans.

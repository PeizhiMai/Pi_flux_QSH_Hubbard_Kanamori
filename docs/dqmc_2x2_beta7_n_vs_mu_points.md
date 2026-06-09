# 2×2 DQMC points on ED `n(μ)` curves at β=7

This run places DQMC data points on the 2×2 ED benchmark curves from
`docs/ed_2x2_beta7_benchmark.md`.  The interaction convention is physical
(non-PH-shifted), matching the BHZ project.

## Run parameters

```text
Lx = Ly = 2
BCs = pbc, cylinder/open_x
t = 1.0, lambda = 0.2
U = 1.0, JH = 0.25
beta = 7.0, dtau = 0.1
Ntherm = 50, Nmeas = 100, Nupdates = 1, n_stab = 10
```

Ten common chemical-potential points were used for every tier/BC curve:

```text
mu = -1.5, -1.0, -0.5, 0.0, 0.5, 0.875, 1.25, 1.5, 2.0, 2.5
```

This includes the physical half-filling points:

- Hubbard: `mu_half = U/2 = 0.5`;
- Kanamori tiers: `mu_half = (3U - 5JH)/2 = 0.875`.

The run was performed with:

```bash
~/.venvs/myenv/bin/python scripts/run_dqmc_n_vs_mu_points.py \
  --Ntherm=50 --Nmeas=100 --Nupdates=1 --n-stab=10 \
  --beta=7.0 --dtau=0.1
```

The runner is resumable and reuses existing point directories unless `--rerun` is
passed.

## Local output files

Raw and processed local outputs are under:

```text
/Users/cosdis/Desktop/projects/strongly_correlated_topological_systems/Pi_flux_QSH_Hubbard_Kanamori/results/dqmc_2x2_beta7_n_vs_mu_points/
```

Important files:

- `dqmc_points.tsv` — all 80 DQMC points merged with ED densities;
- `curve_summary.tsv` — per-tier/per-BC summary;
- `ed_dqmc_n_vs_mu_overlay.png` — ED curves with DQMC points;
- `dqmc_average_phase_magnitude.png` — average phase magnitude for the same points.

## Half-filling summary

| tier | BC | μ_half | DQMC density/cell | ED density/cell | density − ED | average phase | `|phase|` | min `|phase|` on curve |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| Hubbard | PBC | 0.5 | 2.000000 | 2.000000 | +0.000000 | `1.000000 + 0.000000i` | 1.000000 | 0.999346 |
| Hubbard | cylinder | 0.5 | 2.000000 | 2.000000 | +0.000000 | `1.000000 - 0.000000i` | 1.000000 | 0.996095 |
| Density/Ising | PBC | 0.875 | 1.999939 | 2.000000 | -0.000061 | `0.999819 + 0.002029i` | 0.999821 | 0.988085 |
| Density/Ising | cylinder | 0.875 | 1.999933 | 2.000000 | -0.000067 | `0.999701 - 0.006250i` | 0.999720 | 0.975804 |
| Spin-flip | PBC | 0.875 | 2.000243 | 2.000000 | +0.000243 | `0.813313 - 0.034320i` | 0.814037 | 0.714008 |
| Spin-flip | cylinder | 0.875 | 1.999983 | 2.000000 | -0.000017 | `0.719315 - 0.043572i` | 0.720633 | 0.497280 |
| Full | PBC | 0.875 | 1.999980 | 2.000000 | -0.000020 | `0.247760 - 0.042570i` | 0.251391 | 0.134381 |
| Full | cylinder | 0.875 | 1.999376 | 2.000000 | -0.000624 | `-0.124506 + 0.021078i` | 0.126278 | 0.037868 |

## Takeaways

- The DQMC densities reproduce the ED `n(μ)` curves well at this low-stat level.
- At physical half filling, all tiers/BCs remain pinned to density/cell `≈2`.
- The sign/phase hierarchy is clear at β=7 on 2×2:
  1. Hubbard-only: phase ≈ 1 throughout the curve;
  2. density/Ising Kanamori: still benign, phase close to 1;
  3. spin-flip Hund: moderate phase degradation;
  4. full spin-flip + pair hopping: severe phase degradation, especially on the cylinder.
- The full-tier density points should be treated as exploratory because the phase
  is already small with only `Nmeas=100`; production comparisons should increase
  statistics and likely use binning/error bars.

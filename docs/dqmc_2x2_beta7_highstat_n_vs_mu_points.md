# High-stat 2×2 DQMC points on ED `n(μ)` curves at β=7

This is the higher-statistics repeat of `docs/dqmc_2x2_beta7_n_vs_mu_points.md`.
It uses the same 80 DQMC points, but with more thermalization and measurements:

```text
Lx = Ly = 2
BCs = pbc, cylinder/open_x
t = 1.0, lambda = 0.2
U = 1.0, JH = 0.25
beta = 7.0, dtau = 0.1
Nwarmup/Ntherm = 500
Nmeas = 3000
Nupdates = 1, n_stab = 10
```

Chemical potentials for every tier/BC curve:

```text
mu = -1.5, -1.0, -0.5, 0.0, 0.5, 0.875, 1.25, 1.5, 2.0, 2.5
```

Run command:

```bash
OUTDIR=/Users/cosdis/Desktop/projects/strongly_correlated_topological_systems/Pi_flux_QSH_Hubbard_Kanamori/results/dqmc_2x2_beta7_n_vs_mu_points_Nwarmup500_Nmeas3000
~/.venvs/myenv/bin/python scripts/run_dqmc_n_vs_mu_points.py \
  --outdir "$OUTDIR" \
  --Ntherm=500 --Nmeas=3000 --Nupdates=1 --n-stab=10 \
  --beta=7.0 --dtau=0.1
```

Local outputs:

```text
/Users/cosdis/Desktop/projects/strongly_correlated_topological_systems/Pi_flux_QSH_Hubbard_Kanamori/results/dqmc_2x2_beta7_n_vs_mu_points_Nwarmup500_Nmeas3000/
```

Important files:

- `dqmc_points.tsv` — all 80 high-stat points merged with ED densities;
- `curve_summary.tsv` — per-tier/per-BC summary;
- `ed_dqmc_n_vs_mu_overlay.png` — ED curves with high-stat DQMC points;
- `dqmc_average_phase_magnitude.png` — average phase magnitude.

## Half-filling summary

| tier | BC | μ_half | DQMC density/cell | ED density/cell | density − ED | average phase | `|phase|` | min `|phase|` on curve |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| Hubbard | PBC | 0.5 | 2.000000 | 2.000000 | +0.000000 | `1.000000 - 0.000000i` | 1.000000 | 0.999495 |
| Hubbard | cylinder | 0.5 | 2.000000 | 2.000000 | +0.000000 | `1.000000 - 0.000000i` | 1.000000 | 0.997098 |
| Density/Ising | PBC | 0.875 | 2.000015 | 2.000000 | +0.000015 | `0.999755 - 0.000336i` | 0.999755 | 0.996499 |
| Density/Ising | cylinder | 0.875 | 1.999999 | 2.000000 | -0.000001 | `0.999696 + 0.000936i` | 0.999697 | 0.980799 |
| Spin-flip | PBC | 0.875 | 2.000006 | 2.000000 | +0.000006 | `0.811594 - 0.001664i` | 0.811596 | 0.716308 |
| Spin-flip | cylinder | 0.875 | 2.000024 | 2.000000 | +0.000024 | `0.654264 + 0.005278i` | 0.654286 | 0.622884 |
| Full | PBC | 0.875 | 1.999951 | 2.000000 | -0.000049 | `0.241140 - 0.010463i` | 0.241367 | 0.207529 |
| Full | cylinder | 0.875 | 2.000260 | 2.000000 | +0.000260 | `0.072234 - 0.004021i` | 0.072345 | 0.044811 |

## Observations

- The high-stat DQMC density points sit on the ED curves very well.
- Hubbard-only remains sign-free at half filling and essentially phase-free across this small 2×2 scan.
- Density/Ising Kanamori remains benign at β=7 on 2×2.
- Spin-flip Hund degrades the average phase but remains usable on this cluster.
- The full spin-flip + pair-hopping tier is already a severe sign/phase problem, especially for the cylinder: the half-filled phase magnitude is only `≈0.072`.
- The high-stat run strengthens the tier hierarchy seen in the low-stat run: `Hubbard ≈ density >> spinflip >> full`.

# β=10 sign-problem check

This is the first low-temperature sign-problem smoke check for the minimal
π-flux QSH Hubbard-Kanamori setup.  It is intentionally a small `2x2` DQMC
run, not a production statistics run.

## Parameters

```text
Nx = 2, Ly = 2
t = 1.0, lambda = 0.2
U = 1.0, JH = 0.25
mu = 0.0
beta = 10.0, dtau = 0.1
Ntherm = 50, Nmeas = 100, Nupdates = 1
```

## Result

| tier | average phase | density/cell | comment |
|---|---:|---:|---|
| Hubbard only | `1 - 6.37e-18 i` | `2.000000` | Sign-free parent verified at β=10. |
| Density/Ising Kanamori | `0.999683 - 0.003516 i` | `2.000002` | Still essentially benign for this small run. |
| Spin-flip Hund | `0.768563 - 0.042049 i` | `1.999914` | Clear phase/sign degradation. |
| Full transverse Kanamori | `0.110490 - 0.042689 i` | `1.999842` | Severe phase/sign problem already on the small cluster. |

Raw local output was written under `/private/tmp/piflux_beta10_sign_scan` during the
setup run.  Reproduce with:

```bash
OUTDIR=/private/tmp/piflux_beta10_sign_scan BETAS="10.0" JHS="0.25" \
  NTHERM=50 NMEAS=100 ./scripts/run_sign_scan.sh
```

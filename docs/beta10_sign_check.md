# β=10 sign-problem check

The low-temperature sign-problem check should use `beta=10`.  The preferred
small-cluster check is now `Lx=Ly=4`; the Hamiltonian remains at half filling
with `mu=0`.

## Preferred Lx=Ly=4 check

Parameters:

```text
Lx = 4, Ly = 4
t = 1.0, lambda = 0.2
U = 1.0, JH = 0.25
mu = 0.0
beta = 10.0, dtau = 0.1
Ntherm = 50, Nmeas = 100, Nupdates = 1
```

Result from the setup run:

| tier | average phase | density/cell | comment |
|---|---:|---:|---|
| Hubbard only | `1 + 5.71e-16 i` | `2.000000` | Sign-free parent verified at β=10. |
| Density/Ising Kanamori | `0.991289 - 0.032546 i` | `2.000235` | Mild phase degradation. |
| Spin-flip Hund | `0.231718 - 0.021749 i` | `1.987817` | Severe phase degradation. |
| Full transverse Kanamori | `-0.060602 + 0.099462 i` | `2.000872` | Very severe phase/sign problem. |

Raw local output was written under `/private/tmp/piflux_beta10_L4_sign_scan`
during the setup run.  Reproduce with:

```bash
OUTDIR=/private/tmp/piflux_beta10_L4_sign_scan BETAS="10.0" JHS="0.25" Lx=4 Ly=4 \
  NTHERM=50 NMEAS=100 ./scripts/run_sign_scan.sh
```

## Earlier 2x2 smoke check

An earlier `Lx=Ly=2` β=10 smoke showed the same hierarchy: Hubbard-only was
phase `1` to numerical precision, density/Ising was benign, spin-flip degraded
the phase, and full transverse Kanamori was already severe.

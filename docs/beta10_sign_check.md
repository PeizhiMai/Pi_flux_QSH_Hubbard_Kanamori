# β=10 sign-problem check

The low-temperature sign-problem check should use `beta=10` and the preferred
small-cluster size `Lx=Ly=4`.

**Current convention:** ED and DQMC now use physical, non-PH-shifted interactions.
Therefore half filling requires a tier-dependent chemical potential:

```text
Hubbard-only:      mu = U/2
Kanamori tiers:    mu = (3U - 5JH)/2
```

For `U=1.0, JH=0.25`, use `mu=0.5` for Hubbard-only and `mu=0.875` for the
Kanamori tiers.  `scripts/run_sign_scan.sh` computes and passes these values
automatically and records `mu` in the summary TSV.

## Recommended workflow before new β=10 sign scans

Before doing further sign-problem checks, run the systematic 2×2 ED `n(mu)`
benchmark at `beta=7` for both periodic and cylindrical boundary conditions and
all four tiers (`hubbard`, `density`, `spinflip`, `full`).  This verifies the
physical-interaction half-filling convention before spending time on larger DQMC
sign scans.

## Reproduction command for the physical-convention scan

```bash
OUTDIR=/private/tmp/piflux_beta10_L4_sign_scan BETAS="10.0" JHS="0.25" Lx=4 Ly=4   NTHERM=50 NMEAS=100 ./scripts/run_sign_scan.sh
```

## Archived old setup-run numbers

The earlier setup-run numbers were generated before this convention correction,
with the PH-shifted interaction convention and `mu=0`.  They should **not** be
used as physical-interaction results.  They only established that the DQMC driver
ran and that the qualitative phase hierarchy was visible in the old convention.

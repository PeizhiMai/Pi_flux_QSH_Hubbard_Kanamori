# Pi-flux QSH Hubbard-Kanamori

Standalone Julia ED + SmoQyDQMC setup for a minimal two-orbital square-lattice
π-flux quantum spin Hall parent with onsite Hubbard and Kanamori interaction tiers.

The goal is a controlled comparison to BHZ-Hubbard-Kanamori.  The interactions
are kept in the same **physical, non-particle-hole-shifted** form as the BHZ
project: `U n↑ n↓` for the onsite Hubbard term and physical density-density /
Hund / pair-hopping Kanamori terms.  At physical half filling this means using
`mu=U/2` for Hubbard-only and `mu=(3U-5JH)/2` for the Kanamori tiers.

The π-flux QSH parent has the antiunitary/bipartite structure needed for a
Hubbard-only sign-problem-free baseline at the physical half-filling chemical
potential; density/Ising Hund, spin-flip Hund, and pair-hopping tiers are then
added to measure average phase/sign degradation.

## Layout

- `ED/` — finite-temperature grand-canonical exact diagonalization benchmarks.
- `DQMC/SmoqyDQMC/` — SmoQyDQMC driver and noninteracting/topology checks.
- `docs/model_conventions.md` — Hamiltonian, topology, interaction and sign conventions.
- `docs/ed_2x2_beta7_benchmark.md` — 2×2 ED `n(mu)` benchmark for all tiers and PBC/cylinder BCs.
- `scripts/health_check.sh` — bounded local ED/topology/DQMC smoke checks.
- `vendor/SmoQyDQMC.jl/` — vendored SmoQyDQMC with the Kanamori HST extensions.

## Quick check

```bash
cd /Users/cosdis/Desktop/projects/strongly_correlated_topological_systems/Pi_flux_QSH_Hubbard_Kanamori
scripts/health_check.sh --ed-only      # ED + H0 checks only
scripts/health_check.sh                # ED + H0 + tiny DQMC tier smokes
```

The scripts default to a project-local Julia depot first, then reuse the existing
CE_GCE depot if available.

## β=10 sign check

The sign-problem check should be run at low temperature with physical half-filling
chemical potentials.  The default sign-scan script uses `beta=10.0`, `Lx=Ly=4`,
and passes the tier-dependent half-filling `mu` automatically:

```bash
OUTDIR=/private/tmp/piflux_beta10_sign_scan JHS="0.25" ./scripts/run_sign_scan.sh
```

Before doing new production sign scans, run the systematic 2×2 ED `n(mu)` checks
for both periodic and cylindrical boundary conditions at `beta=7`.

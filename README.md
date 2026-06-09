# Pi-flux QSH Hubbard-Kanamori

Standalone Julia ED + SmoQyDQMC setup for a minimal two-orbital square-lattice
π-flux quantum spin Hall parent with onsite Hubbard and Kanamori interaction tiers.

The goal is a controlled comparison to BHZ-Hubbard-Kanamori: the onsite-Hubbard
π-flux QSH parent is particle-hole symmetric and sign-problem-free at half filling,
while density/Ising Hund, spin-flip Hund, and pair-hopping tiers can be added to
measure how the average phase/sign degrades.

## Layout

- `ED/` — finite-temperature grand-canonical exact diagonalization benchmarks.
- `DQMC/SmoqyDQMC/` — SmoQyDQMC driver and noninteracting/topology checks.
- `docs/model_conventions.md` — Hamiltonian, topology, PH/sign conventions.
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

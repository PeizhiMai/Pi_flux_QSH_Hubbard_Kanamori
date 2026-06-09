# Pi-flux QSH SmoQyDQMC setup

Use the vendored SmoQyDQMC package in `../../vendor/SmoQyDQMC.jl` from the repository root.

```bash
cd /Users/cosdis/Desktop/projects/strongly_correlated_topological_systems/Pi_flux_QSH_Hubbard_Kanamori
export JULIA_BIN=/Users/cosdis/.julia/juliaup/julia-1.12.1+0.aarch64.apple.darwin14/bin/julia
export JULIA_DEPOT_PATH="$PWD/.julia_depot:/Users/cosdis/Desktop/projects/CE_GCE/.julia_depot"
export JULIA_PROJECT="$PWD/vendor/SmoQyDQMC.jl"
```

Noninteracting/topology check:

```bash
$JULIA_BIN DQMC/SmoqyDQMC/scripts/check_piflux_qsh_h0.jl --Lx=6 --Ly=6 --lambda=0.2 --grid=41
```

The DQMC driver uses physical interactions, matching the BHZ project:
`HubbardModel(ph_sym_form=false)` and `KanamoriDensityModel(ph_sym_form=false)`.
For half filling pass `--mu=U/2` for Hubbard-only and
`--mu=(3U-5JH)/2` for Kanamori tiers.  With `U=1.0, JH=0.25`, this means
`--mu=0.5` for Hubbard-only and `--mu=0.875` for the three Kanamori tiers.

Tiny Hubbard-only sign-free smoke:

```bash
$JULIA_BIN --project=$JULIA_PROJECT DQMC/SmoqyDQMC/scripts/run_piflux_qsh_smoqy.jl \
  --Lx=2 --Ly=2 --beta=0.2 --dtau=0.1 --U=1.0 --mu=0.5 \
  --Ntherm=1 --Nmeas=2 --Nupdates=1 --n_stab=1 \
  --density_kanamori=false
```

Kanamori tier smokes:

```bash
# Tier 1 density/Ising Hund; mu=(3U-5JH)/2=0.875
$JULIA_BIN --project=$JULIA_PROJECT DQMC/SmoqyDQMC/scripts/run_piflux_qsh_smoqy.jl \
  --Lx=2 --Ly=2 --beta=0.2 --dtau=0.1 --U=1.0 --JH=0.25 --mu=0.875 \
  --Ntherm=1 --Nmeas=2 --Nupdates=1 --n_stab=1 \
  --density_kanamori=true

# Tier 2 spin-flip Hund
$JULIA_BIN --project=$JULIA_PROJECT DQMC/SmoqyDQMC/scripts/run_piflux_qsh_smoqy.jl \
  --Lx=2 --Ly=2 --beta=0.2 --dtau=0.1 --U=1.0 --JH=0.25 --mu=0.875 \
  --Ntherm=1 --Nmeas=2 --Nupdates=1 --n_stab=1 \
  --density_kanamori=true --spin_flip_hund=true

# Tier 3 spin flip plus pair hopping
$JULIA_BIN --project=$JULIA_PROJECT DQMC/SmoqyDQMC/scripts/run_piflux_qsh_smoqy.jl \
  --Lx=2 --Ly=2 --beta=0.2 --dtau=0.1 --U=1.0 --JH=0.25 --mu=0.875 \
  --Ntherm=1 --Nmeas=2 --Nupdates=1 --n_stab=1 \
  --density_kanamori=true --spin_flip_hund=true --pair_hopping=true
```

#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
JULIA_BIN="${JULIA_BIN:-/Users/cosdis/.julia/juliaup/julia-1.12.1+0.aarch64.apple.darwin14/bin/julia}"
if [[ ! -x "$JULIA_BIN" ]]; then JULIA_BIN="$(command -v julia)"; fi
export JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH:-$ROOT/.julia_depot:/Users/cosdis/Desktop/projects/CE_GCE/.julia_depot}"
export JULIA_PROJECT="$ROOT/vendor/SmoQyDQMC.jl"
"$JULIA_BIN" --project="$JULIA_PROJECT" -e 'using Pkg; Pkg.instantiate()' >/dev/null
cd "$ROOT"
OUTDIR="${OUTDIR:-$ROOT/results/ed_vs_dqmc_smoke}"
rm -rf "$OUTDIR"; mkdir -p "$OUTDIR"
TSV="$OUTDIR/ed_vs_dqmc_density_energy.tsv"
echo -e "level\tmu\ted_density_per_cell\ted_energy\tdqmc_density_per_cell\tdqmc_interaction_energy_per_cell\tdqmc_average_phase_re\tdqmc_average_phase_im" > "$TSV"
half_mu() {
  local level="$1"
  python3 - <<PY
U=1.0; JH=0.25; level='$level'
mu = U/2 if level == 'hubbard' else (3*U - 5*JH)/2
print(f'{mu:.16g}')
PY
}
run_level() {
  local level="$1" sid="$2"
  local flags=(--density_kanamori=false)
  case "$level" in
    hubbard) flags=(--density_kanamori=false) ;;
    density) flags=(--density_kanamori=true) ;;
    spinflip) flags=(--density_kanamori=true --spin_flip_hund=true) ;;
    full) flags=(--density_kanamori=true --spin_flip_hund=true --pair_hopping=true) ;;
  esac
  local mu; mu="$(half_mu "$level")"
  "$JULIA_BIN" --project=ED ED/scripts/run_piflux_qsh_ed.jl --Lx=1 --Ly=2 --t=0.2 --lambda=0.05 --U=1.0 --JH=0.25 --beta=0.3 --mu="$mu" --interaction_level="$level" --outdir="$OUTDIR/ed_$level" >/dev/null
  "$JULIA_BIN" --project="$JULIA_PROJECT" DQMC/SmoqyDQMC/scripts/run_piflux_qsh_smoqy.jl \
    --Lx=1 --Ly=2 --t=0.2 --lambda=0.05 --U=1.0 --JH=0.25 --beta=0.3 --dtau=0.1 --mu="$mu" \
    --Ntherm=2 --Nmeas=4 --Nupdates=1 --n_stab=1 --outdir="$OUTDIR/dqmc" --sID="$sid" "${flags[@]}" >/dev/null
  local summary; summary="$(find "$OUTDIR/dqmc" -name summary.txt | sort | tail -1)"
  python3 - <<PY >> "$TSV"
import pathlib, re, tomllib
ed=tomllib.loads(pathlib.Path('$OUTDIR/ed_$level/ed_summary.toml').read_text())
text=pathlib.Path('$summary').read_text()
ph=re.search(r'average_phase =\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s*([+-](?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)i', text)
dn=re.search(r'density_per_cell_reweighted = ([^\n]+)', text)
ei=re.search(r'interaction_energy_per_cell_reweighted =\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)', text)
print('\t'.join(['$level', '$mu', str(ed['observables']['density_per_cell']), str(ed['observables']['energy']), dn.group(1), ei.group(1), ph.group(1), ph.group(2)]))
PY
}
run_level hubbard 1
run_level density 2
run_level spinflip 3
run_level full 4
echo "wrote $TSV"

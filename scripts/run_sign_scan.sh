#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
JULIA_BIN="${JULIA_BIN:-/Users/cosdis/.julia/juliaup/julia-1.12.1+0.aarch64.apple.darwin14/bin/julia}"
if [[ ! -x "$JULIA_BIN" ]]; then JULIA_BIN="$(command -v julia)"; fi
export JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH:-$ROOT/.julia_depot:/Users/cosdis/Desktop/projects/CE_GCE/.julia_depot}"
export JULIA_PROJECT="$ROOT/vendor/SmoQyDQMC.jl"
"$JULIA_BIN" --project="$JULIA_PROJECT" -e 'using Pkg; Pkg.instantiate()' >/dev/null
cd "$ROOT"
OUTDIR="${OUTDIR:-$ROOT/results/sign_scan}"
mkdir -p "$OUTDIR"
TSV="$OUTDIR/sign_scan_summary.tsv"
echo -e "tier\tbeta\tJH\taverage_phase_re\taverage_phase_im\tdensity_per_cell\tsummary" > "$TSV"
run_one() {
  local tier="$1" beta="$2" jh="$3" sid="$4"
  local flags=(--density_kanamori=false)
  case "$tier" in
    hubbard) flags=(--density_kanamori=false) ;;
    density) flags=(--density_kanamori=true) ;;
    spinflip) flags=(--density_kanamori=true --spin_flip_hund=true) ;;
    full) flags=(--density_kanamori=true --spin_flip_hund=true --pair_hopping=true) ;;
  esac
  "$JULIA_BIN" --project="$JULIA_PROJECT" DQMC/SmoqyDQMC/scripts/run_piflux_qsh_smoqy.jl \
    --Nx=2 --Ly=2 --t=1.0 --lambda=0.2 --U=1.0 --JH="$jh" --beta="$beta" --dtau=0.1 \
    --Ntherm="${NTHERM:-20}" --Nmeas="${NMEAS:-40}" --Nupdates=1 --n_stab=2 --outdir="$OUTDIR/raw" --sID="$sid" "${flags[@]}"
  local summary
  summary="$(find "$OUTDIR/raw" -name summary.txt | sort | tail -1)"
  python3 - <<PY >> "$TSV"
import pathlib, re
text=pathlib.Path('$summary').read_text()
ph=re.search(r'average_phase =\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s*([+-](?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)i', text)
dn=re.search(r'density_per_cell_reweighted = ([^\n]+)', text)
print('\t'.join(['$tier','$beta','$jh', ph.group(1), ph.group(2), dn.group(1), '$summary']))
PY
}
sid=1
for beta in ${BETAS:-0.5 1.0}; do
  for jh in ${JHS:-0.0 0.1 0.25}; do
    for tier in hubbard density spinflip full; do
      run_one "$tier" "$beta" "$jh" "$sid"
      sid=$((sid+1))
    done
  done
done
echo "wrote $TSV"

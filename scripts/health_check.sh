#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
JULIA_BIN="${JULIA_BIN:-/Users/cosdis/.julia/juliaup/julia-1.12.1+0.aarch64.apple.darwin14/bin/julia}"
if [[ ! -x "$JULIA_BIN" ]]; then
  JULIA_BIN="$(command -v julia)"
fi
export JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH:-$ROOT/.julia_depot:/Users/cosdis/Desktop/projects/CE_GCE/.julia_depot}"
export JULIA_PROJECT="$ROOT/vendor/SmoQyDQMC.jl"
"$JULIA_BIN" --project="$JULIA_PROJECT" -e 'using Pkg; Pkg.instantiate()' >/dev/null
ED_ONLY=false
if [[ "${1:-}" == "--ed-only" ]]; then
  ED_ONLY=true
fi

cd "$ROOT"
echo "[health] Julia: $JULIA_BIN"
echo "[health] ED tests"
"$JULIA_BIN" --project=ED ED/test/run_ed_tests.jl

echo "[health] H0 topology/sign checks"
"$JULIA_BIN" DQMC/SmoqyDQMC/scripts/check_piflux_qsh_h0.jl --Nx=4 --Ly=4 --lambda=0.2 --grid=31

if [[ "$ED_ONLY" == true ]]; then
  echo "[health] --ed-only requested; skipping DQMC smokes"
  exit 0
fi

SMOKE_OUTDIR="${SMOKE_OUTDIR:-/private/tmp/piflux_qsh_hubbard_kanamori_smoke}"
rm -rf "$SMOKE_OUTDIR"
mkdir -p "$SMOKE_OUTDIR"
common=(--Nx=2 --Ly=2 --t=1.0 --lambda=0.2 --U=1.0 --JH=0.25 --beta=0.2 --dtau=0.1 --Ntherm=1 --Nmeas=2 --Nupdates=1 --n_stab=1 --outdir="$SMOKE_OUTDIR")

echo "[health] DQMC Hubbard-only sign-free smoke"
"$JULIA_BIN" --project="$JULIA_PROJECT" DQMC/SmoqyDQMC/scripts/run_piflux_qsh_smoqy.jl "${common[@]}" --density_kanamori=false --sID=1
hub_summary="$(find "$SMOKE_OUTDIR" -name summary.txt | sort | head -1)"
python3 - <<PY
import re, sys, pathlib
text = pathlib.Path('$hub_summary').read_text()
m = re.search(r'average_phase =\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s*([+-](?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)i', text)
if not m:
    raise SystemExit('could not parse average_phase')
reph = float(m.group(1)); imph = float(m.group(2))
print(f'[health] Hubbard average_phase={reph:+.12g}{imph:+.12g}i')
if abs(reph - 1.0) > 1e-8 or abs(imph) > 1e-8:
    raise SystemExit('Hubbard-only average phase is not approximately 1')
PY

echo "[health] DQMC Tier 1 density Kanamori smoke"
"$JULIA_BIN" --project="$JULIA_PROJECT" DQMC/SmoqyDQMC/scripts/run_piflux_qsh_smoqy.jl "${common[@]}" --density_kanamori=true --sID=2

echo "[health] DQMC Tier 2 spin-flip Hund smoke"
"$JULIA_BIN" --project="$JULIA_PROJECT" DQMC/SmoqyDQMC/scripts/run_piflux_qsh_smoqy.jl "${common[@]}" --density_kanamori=true --spin_flip_hund=true --sID=3

echo "[health] DQMC Tier 3 full transverse Kanamori smoke"
"$JULIA_BIN" --project="$JULIA_PROJECT" DQMC/SmoqyDQMC/scripts/run_piflux_qsh_smoqy.jl "${common[@]}" --density_kanamori=true --spin_flip_hund=true --pair_hopping=true --sID=4

echo "[health] OK; smoke output: $SMOKE_OUTDIR"

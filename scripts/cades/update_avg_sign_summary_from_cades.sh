#!/usr/bin/env bash
# Pull current pi-flux average-sign summary from CADES and regenerate the PDF/PNG.
set -euo pipefail

REPO_DIR=${REPO_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}
REMOTE=${REMOTE:-9pm@or-login.ornl.gov}
REMOTE_ROOT=${REMOTE_ROOT:-/home/9pm/Pi_flux_QSH_Hubbard_Kanamori}
REMOTE_RUN_ROOT=${REMOTE_RUN_ROOT:-runs/piflux_L8x4_beta7_10_12_16_half_filling_avg_sign}
SSH_OPTS=(-o BatchMode=yes -o IdentitiesOnly=yes -i "$HOME/.ssh/id_ed25519_cosdis" -o "ProxyCommand=ssh -o BatchMode=yes -o IdentitiesOnly=yes -i $HOME/.ssh/id_ed25519 -W %h:%p ornl-bastion")
RSYNC_RSH="ssh -o BatchMode=yes -o IdentitiesOnly=yes -i $HOME/.ssh/id_ed25519_cosdis -o ProxyCommand='ssh -o BatchMode=yes -o IdentitiesOnly=yes -i $HOME/.ssh/id_ed25519 -W %h:%p ornl-bastion'"
PYTHON=${PYTHON:-$HOME/.venvs/myenv/bin/python}

cd "$REPO_DIR"
mkdir -p results/avg_sign_summary docs

ssh "${SSH_OPTS[@]}" "$REMOTE" "cd '${REMOTE_ROOT}' && python3 scripts/cades/summarize_piflux_avg_sign.py --root '${REMOTE_RUN_ROOT}' --out '${REMOTE_RUN_ROOT}/summary.tsv'"
rsync -az -e "$RSYNC_RSH" "$REMOTE:$REMOTE_ROOT/$REMOTE_RUN_ROOT/summary.tsv" results/avg_sign_summary/piflux_beta7_10_12_16_half_filling_summary.tsv

"$PYTHON" scripts/cades/make_avg_sign_summary_pdf.py \
  --summary results/avg_sign_summary/piflux_beta7_10_12_16_half_filling_summary.tsv \
  --out docs/piflux_kanamori_avg_sign_summary.pdf \
  --png docs/piflux_kanamori_avg_sign_summary.png

"$PYTHON" - <<'PY'
from pathlib import Path
p = Path('results/avg_sign_summary/piflux_beta7_10_12_16_half_filling_summary.tsv')
if not p.exists():
    print('no summary found')
else:
    rows = max(0, len(p.read_text().splitlines()) - 1)
    print(f'completed summary rows: {rows}/16')
PY

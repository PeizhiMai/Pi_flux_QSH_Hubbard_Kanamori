#!/usr/bin/env python3
"""Run small DQMC n(mu) points to overlay on the 2x2 ED benchmark curves.

The script is intentionally resumable: each DQMC point gets its own directory,
and an existing summary.txt is parsed instead of rerunning unless --rerun is used.
"""
from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path
import re
import subprocess
import sys
from typing import Dict, Tuple

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_JULIA = "/Users/cosdis/.julia/juliaup/julia-1.12.1+0.aarch64.apple.darwin14/bin/julia"
DEFAULT_DEPOT = f"{ROOT / '.julia_depot'}:/Users/cosdis/Desktop/projects/CE_GCE/.julia_depot"
TIERS = ("hubbard", "density", "spinflip", "full")
BCS = ("pbc", "cylinder")
DEFAULT_MU_POINTS = "-1.5,-1.0,-0.5,0.0,0.5,0.875,1.25,1.5,2.0,2.5"

FLOAT_RE = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"
PHASE_RE = re.compile(rf"average_phase =\s*({FLOAT_RE})\s*([+-](?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)i")
ABS_RE = re.compile(rf"average_abs_phase =\s*({FLOAT_RE})")
ACC_RE = re.compile(rf"acceptance_rate =\s*({FLOAT_RE})")
DENSITY_RE = re.compile(rf"density_per_cell_reweighted =\s*({FLOAT_RE})")
EINT_RE = re.compile(rf"interaction_energy_per_cell_reweighted =\s*({FLOAT_RE})\s*([+-](?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)i")


def parse_csv_list(s: str, cast=str):
    return [cast(x.strip()) for x in s.split(",") if x.strip()]


def mu_tag(mu: float) -> str:
    return f"{mu:+.6f}".replace("+", "p").replace("-", "m").replace(".", "p")


def tier_flags(tier: str):
    if tier == "hubbard":
        return ["--density_kanamori=false"]
    if tier == "density":
        return ["--density_kanamori=true"]
    if tier == "spinflip":
        return ["--density_kanamori=true", "--spin_flip_hund=true"]
    if tier == "full":
        return ["--density_kanamori=true", "--spin_flip_hund=true", "--pair_hopping=true"]
    raise ValueError(f"unknown tier {tier}")


def bc_open_x(bc: str) -> bool:
    if bc == "pbc":
        return False
    if bc in ("cylinder", "cyl", "open_x"):
        return True
    raise ValueError(f"unknown boundary condition {bc}")


def load_ed_density(path: Path) -> Dict[Tuple[str, str, float], float]:
    out: Dict[Tuple[str, str, float], float] = {}
    if not path.exists():
        return out
    with path.open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            out[(row["bc"], row["level"], round(float(row["mu"]), 12))] = float(row["density_per_cell"])
    return out


def find_summary(point_dir: Path) -> Path | None:
    summaries = sorted(point_dir.glob("**/summary.txt"), key=lambda p: p.stat().st_mtime)
    return summaries[-1] if summaries else None


def parse_summary(path: Path):
    text = path.read_text()
    phase = PHASE_RE.search(text)
    density = DENSITY_RE.search(text)
    if phase is None or density is None:
        raise RuntimeError(f"Could not parse required fields from {path}")
    abs_phase = ABS_RE.search(text)
    acc = ACC_RE.search(text)
    eint = EINT_RE.search(text)
    return {
        "average_phase_re": float(phase.group(1)),
        "average_phase_im": float(phase.group(2)),
        "average_abs_phase": float(abs_phase.group(1)) if abs_phase else float("nan"),
        "acceptance_rate": float(acc.group(1)) if acc else float("nan"),
        "density_per_cell_reweighted": float(density.group(1)),
        "interaction_energy_re": float(eint.group(1)) if eint else float("nan"),
        "interaction_energy_im": float(eint.group(2)) if eint else float("nan"),
    }


def write_rows(tsv: Path, rows):
    fields = [
        "tier", "bc", "open_x", "beta", "dtau", "Lx", "Ly", "t", "lambda", "U", "JH", "mu",
        "ed_density_per_cell", "dqmc_density_per_cell", "density_minus_ed",
        "average_phase_re", "average_phase_im", "average_abs_phase", "acceptance_rate",
        "interaction_energy_re", "interaction_energy_im", "summary_path",
    ]
    tsv.parent.mkdir(parents=True, exist_ok=True)
    with tsv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def run_point(args, env, tier: str, bc: str, mu: float, sid: int, point_dir: Path):
    point_dir.mkdir(parents=True, exist_ok=True)
    log_path = point_dir / "run.log"
    cmd = [
        args.julia_bin,
        f"--project={ROOT / 'vendor' / 'SmoQyDQMC.jl'}",
        str(ROOT / "DQMC" / "SmoqyDQMC" / "scripts" / "run_piflux_qsh_smoqy.jl"),
        f"--Lx={args.Lx}", f"--Ly={args.Ly}",
        f"--t={args.t}", f"--lambda={args.lam}",
        f"--U={args.U}", f"--JH={args.JH}",
        f"--mu={mu}", f"--beta={args.beta}", f"--dtau={args.dtau}",
        f"--open_x={'true' if bc_open_x(bc) else 'false'}",
        f"--Ntherm={args.Ntherm}", f"--Nmeas={args.Nmeas}",
        f"--Nupdates={args.Nupdates}", f"--n_stab={args.n_stab}",
        f"--seed={args.seed_base + sid}", f"--sID={sid}",
        f"--outdir={point_dir}",
    ] + tier_flags(tier)
    print(f"[DQMC n(mu)] running tier={tier} bc={bc} mu={mu:.12g} sid={sid}", flush=True)
    with log_path.open("w") as log:
        log.write(" ".join(cmd) + "\n\n")
        proc = subprocess.run(cmd, cwd=ROOT, env=env, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        log.write(proc.stdout)
    if proc.stdout:
        sys.stdout.write(proc.stdout)
        sys.stdout.flush()
    if proc.returncode != 0:
        raise RuntimeError(f"DQMC point failed tier={tier} bc={bc} mu={mu}; see {log_path}")
    summary = find_summary(point_dir)
    if summary is None:
        raise RuntimeError(f"No summary.txt found under {point_dir}")
    return summary


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=str(ROOT / "results" / "dqmc_2x2_beta7_n_vs_mu_points"))
    ap.add_argument("--ed-curve", default=str(ROOT / "results" / "ed_2x2_beta7_n_vs_mu" / "n_vs_mu.tsv"))
    ap.add_argument("--levels", default=",".join(TIERS))
    ap.add_argument("--bcs", default=",".join(BCS))
    ap.add_argument("--mu-points", default=DEFAULT_MU_POINTS)
    ap.add_argument("--Lx", type=int, default=2)
    ap.add_argument("--Ly", type=int, default=2)
    ap.add_argument("--t", type=float, default=1.0)
    ap.add_argument("--lam", type=float, default=0.2)
    ap.add_argument("--U", type=float, default=1.0)
    ap.add_argument("--JH", type=float, default=0.25)
    ap.add_argument("--beta", type=float, default=7.0)
    ap.add_argument("--dtau", type=float, default=0.1)
    ap.add_argument("--Ntherm", type=int, default=50)
    ap.add_argument("--Nmeas", type=int, default=100)
    ap.add_argument("--Nupdates", type=int, default=1)
    ap.add_argument("--n-stab", dest="n_stab", type=int, default=10)
    ap.add_argument("--seed-base", type=int, default=913000)
    ap.add_argument("--julia-bin", default=os.environ.get("JULIA_BIN", DEFAULT_JULIA))
    ap.add_argument("--julia-depot-path", default=os.environ.get("JULIA_DEPOT_PATH", DEFAULT_DEPOT))
    ap.add_argument("--max-runs", type=int, default=0, help="debug: stop after this many new runs; 0 means all")
    ap.add_argument("--rerun", action="store_true")
    args = ap.parse_args(argv)

    levels = parse_csv_list(args.levels, str)
    bcs = parse_csv_list(args.bcs, str)
    mus = parse_csv_list(args.mu_points, float)
    for tier in levels:
        if tier not in TIERS:
            raise SystemExit(f"Unknown tier {tier}; valid tiers: {TIERS}")
    for bc in bcs:
        bc_open_x(bc)

    env = os.environ.copy()
    env["JULIA_DEPOT_PATH"] = args.julia_depot_path
    env["JULIA_PROJECT"] = str(ROOT / "vendor" / "SmoQyDQMC.jl")
    outdir = Path(args.outdir)
    raw = outdir / "raw"
    tsv = outdir / "dqmc_points.tsv"
    ed = load_ed_density(Path(args.ed_curve))
    rows = []
    new_runs = 0
    sid = 1
    for tier in levels:
        for bc in bcs:
            for mu in mus:
                open_x = bc_open_x(bc)
                point_dir = raw / tier / bc / f"mu_{mu_tag(mu)}"
                summary = None if args.rerun else find_summary(point_dir)
                if summary is not None:
                    print(f"[DQMC n(mu)] reusing tier={tier} bc={bc} mu={mu:.12g}: {summary}", flush=True)
                else:
                    if args.max_runs and new_runs >= args.max_runs:
                        print(f"[DQMC n(mu)] max new runs reached ({args.max_runs}); stopping early", flush=True)
                        write_rows(tsv, rows)
                        return 0
                    summary = run_point(args, env, tier, bc, mu, sid, point_dir)
                    new_runs += 1
                obs = parse_summary(summary)
                ed_density = ed.get((bc, tier, round(mu, 12)), float("nan"))
                dqmc_density = obs.pop("density_per_cell_reweighted")
                rows.append({
                    "tier": tier, "bc": bc, "open_x": str(open_x).lower(),
                    "beta": args.beta, "dtau": args.dtau, "Lx": args.Lx, "Ly": args.Ly,
                    "t": args.t, "lambda": args.lam, "U": args.U, "JH": args.JH, "mu": mu,
                    "ed_density_per_cell": ed_density,
                    "dqmc_density_per_cell": dqmc_density,
                    "density_minus_ed": dqmc_density - ed_density if ed_density == ed_density else float("nan"),
                    "summary_path": str(summary),
                    **obs,
                })
                sid += 1
                write_rows(tsv, rows)
    print(f"[DQMC n(mu)] wrote {tsv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

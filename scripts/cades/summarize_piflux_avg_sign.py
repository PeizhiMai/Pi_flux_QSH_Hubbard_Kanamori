#!/usr/bin/env python3
"""Aggregate rank-local pi-flux QSH Hubbard-Kanamori DQMC outputs.

The main sign-problem diagnostic is |<phase>|, computed from the complex
average phases of all completed ranks.  The per-rank |<phase>| mean is also
reported as a useful, but upward-biased, stability diagnostic.
"""
from __future__ import annotations

import argparse
import math
import statistics as stats
from pathlib import Path


def parse_toml_scalar(value: str):
    value = value.strip()
    if value.startswith('"') and value.endswith('"'):
        return value[1:-1]
    if value.lower() == "true":
        return True
    if value.lower() == "false":
        return False
    try:
        if any(c in value for c in ".eE"):
            return float(value)
        return int(value)
    except ValueError:
        return value


def load_simple_toml(path: Path) -> dict:
    out = {}
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, value = line.split("=", 1)
        out[key.strip()] = parse_toml_scalar(value)
    return out


def parse_complex(s: str) -> complex:
    return complex(s.strip().replace("i", "j"))


def parse_summary(path: Path) -> dict[str, str]:
    out = {}
    for line in path.read_text().splitlines():
        if "=" not in line:
            continue
        k, v = line.split("=", 1)
        out[k.strip()] = v.strip()
    return out


def as_float(x, default=float("nan")) -> float:
    try:
        return float(x)
    except Exception:
        return default


def sem(xs):
    xs = [x for x in xs if not math.isnan(x)]
    if len(xs) <= 1:
        return float("nan")
    return stats.stdev(xs) / math.sqrt(len(xs))


def tier_from_metadata(md: dict) -> str:
    den = bool(md.get("density_kanamori", False))
    sf = bool(md.get("spin_flip_hund", False))
    pair = bool(md.get("pair_hopping", False))
    if den and sf and pair:
        return "full"
    if den and sf:
        return "spinflip"
    if den:
        return "density"
    return "hubbard"


def summarize_run(run_dir: Path):
    rank_summaries = sorted(run_dir.glob("ranks/rank_*/piflux_qsh_*/summary.txt"))
    if not rank_summaries:
        return None

    metadata = None
    phases: list[complex] = []
    rank_abs: list[float] = []
    densities: list[float] = []
    eints: list[complex] = []
    accs: list[float] = []
    nmeas_total = 0
    complete_count = 0

    for sp in rank_summaries:
        mdp = sp.with_name("metadata.toml")
        if metadata is None and mdp.exists():
            metadata = load_simple_toml(mdp)
        summ = parse_summary(sp)
        if "average_phase" in summ:
            ph = parse_complex(summ["average_phase"])
            phases.append(ph)
            rank_abs.append(abs(ph))
        elif "average_phase_abs" in summ:
            rank_abs.append(as_float(summ["average_phase_abs"]))
        if "density_per_cell_reweighted" in summ:
            densities.append(as_float(summ["density_per_cell_reweighted"]))
        if "interaction_energy_per_cell_reweighted" in summ:
            try:
                eints.append(parse_complex(summ["interaction_energy_per_cell_reweighted"]))
            except Exception:
                pass
        if "acceptance_rate" in summ:
            accs.append(as_float(summ["acceptance_rate"]))
        csp = sp.with_name("complete.status")
        if csp.exists():
            complete_count += 1
            cdict = parse_summary(csp)
            nmeas_total += int(as_float(cdict.get("measurements", 0), 0.0))

    if metadata is None:
        metadata = {}

    avg_phase = sum(phases, 0j) / len(phases) if phases else complex(float("nan"), float("nan"))
    tier = tier_from_metadata(metadata)
    return {
        "run": run_dir.name,
        "tier": tier,
        "Lx": metadata.get("Lx", ""),
        "Ly": metadata.get("Ly", ""),
        "t": metadata.get("t", ""),
        "lambda": metadata.get("lambda", ""),
        "mu": metadata.get("mu", ""),
        "mu_half_filling": metadata.get("mu_half_filling", ""),
        "beta": metadata.get("beta", ""),
        "dtau": metadata.get("dtau", ""),
        "U": metadata.get("U", ""),
        "JH": metadata.get("JH", ""),
        "open_x": metadata.get("open_x", ""),
        "density_kanamori": metadata.get("density_kanamori", ""),
        "spin_flip_hund": metadata.get("spin_flip_hund", ""),
        "pair_hopping": metadata.get("pair_hopping", ""),
        "nranks_complete": len(rank_summaries),
        "nranks_with_complete_status": complete_count,
        "nmeas_rows": nmeas_total,
        "average_phase_re": avg_phase.real,
        "average_phase_im": avg_phase.imag,
        "average_phase_abs_global": abs(avg_phase),
        "average_phase_abs_rank_mean": stats.mean(rank_abs) if rank_abs else float("nan"),
        "average_phase_abs_rank_sem": sem(rank_abs),
        "density_per_cell_reweighted_rank_mean": stats.mean(densities) if densities else float("nan"),
        "density_per_cell_reweighted_rank_sem": sem(densities),
        "interaction_energy_per_cell_reweighted_rank_mean_re": stats.mean([z.real for z in eints]) if eints else float("nan"),
        "interaction_energy_per_cell_reweighted_rank_mean_im": stats.mean([z.imag for z in eints]) if eints else float("nan"),
        "acceptance_mean": stats.mean(accs) if accs else float("nan"),
        "acceptance_sem": sem(accs),
    }


def sort_key(row):
    tier_order = {"hubbard": 0, "density": 1, "spinflip": 2, "full": 3}
    return (float(row.get("beta") or 0), tier_order.get(str(row.get("tier")), 99))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default="runs/piflux_L8x4_beta7_10_12_16_half_filling_avg_sign")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    root = Path(args.root)
    rows = []
    if root.exists():
        rows = [r for d in sorted(root.iterdir()) if d.is_dir() for r in [summarize_run(d)] if r is not None]
    rows.sort(key=sort_key)

    cols = [
        "run", "tier", "Lx", "Ly", "t", "lambda", "mu", "mu_half_filling", "beta", "dtau", "U", "JH", "open_x",
        "density_kanamori", "spin_flip_hund", "pair_hopping", "nranks_complete", "nranks_with_complete_status", "nmeas_rows",
        "average_phase_re", "average_phase_im", "average_phase_abs_global", "average_phase_abs_rank_mean", "average_phase_abs_rank_sem",
        "density_per_cell_reweighted_rank_mean", "density_per_cell_reweighted_rank_sem",
        "interaction_energy_per_cell_reweighted_rank_mean_re", "interaction_energy_per_cell_reweighted_rank_mean_im",
        "acceptance_mean", "acceptance_sem",
    ]
    text = "\t".join(cols) + "\n"
    for r in rows:
        text += "\t".join(str(r.get(c, "")) for c in cols) + "\n"

    if args.out:
        out = Path(args.out)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(text)
    else:
        print(text, end="")


if __name__ == "__main__":
    main()

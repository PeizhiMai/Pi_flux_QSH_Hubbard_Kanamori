#!/usr/bin/env python3
"""Make a one-page avg-phase/sign summary PDF for pi-flux QSH Kanamori scans."""
from __future__ import annotations

import argparse
import csv
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt

BETAS = [7, 10, 12, 16]
TIERS = ["hubbard", "density", "spinflip", "full"]
TIER_LABELS = {
    "hubbard": "Hubbard",
    "density": "density",
    "spinflip": "spin-flip",
    "full": "full",
}
COLORS = {
    "hubbard": "#1f77b4",
    "density": "#2ca02c",
    "spinflip": "#ff7f0e",
    "full": "#d62728",
}
MARKERS = {"hubbard": "o", "density": "s", "spinflip": "^", "full": "D"}


def ffloat(row, key):
    try:
        return float(row[key])
    except Exception:
        return None


def load_rows(path: Path):
    if not path.exists():
        return []
    with path.open() as f:
        return list(csv.DictReader(f, delimiter="\t"))


def value_at(rows, beta, tier, key="average_phase_abs_global"):
    for row in rows:
        rb = ffloat(row, "beta")
        if rb is None:
            continue
        if int(round(rb)) == beta and row.get("tier") == tier:
            return ffloat(row, key)
    return None


def fmt(x, nd=3):
    if x is None:
        return "pending"
    return f"{x:.{nd}f}"


def table(ax, cell_text, col_labels, bbox, font_size=8.0, header_color="#dce9f7"):
    tbl = ax.table(cellText=cell_text, colLabels=col_labels, cellLoc="center", colLoc="center", bbox=bbox)
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(font_size)
    for (r, c), cell in tbl.get_celld().items():
        cell.set_edgecolor("#5d6875")
        cell.set_linewidth(0.45)
        if r == 0:
            cell.set_facecolor(header_color)
            cell.set_text_props(weight="bold")
        elif r % 2 == 0:
            cell.set_facecolor("#f7f9fb")
    return tbl


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--summary", default="results/avg_sign_summary/piflux_beta7_10_12_16_half_filling_summary.tsv")
    ap.add_argument("--out", default="docs/piflux_kanamori_avg_sign_summary.pdf")
    ap.add_argument("--png", default="docs/piflux_kanamori_avg_sign_summary.png")
    args = ap.parse_args()

    rows = load_rows(Path(args.summary))
    out = Path(args.out)
    png = Path(args.png) if args.png else None
    out.parent.mkdir(parents=True, exist_ok=True)
    if png:
        png.parent.mkdir(parents=True, exist_ok=True)

    fig = plt.figure(figsize=(11, 8.5), dpi=200)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.axis("off")
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M")

    fig.text(0.055, 0.955, "Pi-flux QSH Hubbard-Kanamori DQMC: average phase/sign", fontsize=17.5, weight="bold", ha="left")
    fig.text(
        0.055,
        0.925,
        "Lx=8, Ly=4, periodic; U=1.0, J$_H$=0.25, t=1.0, λ=0.2, dtau=0.1; Ntherm=2000, Nmeas=10000 per rank, 32 ranks.",
        fontsize=9.5,
        ha="left",
        color="#333333",
    )
    fig.text(0.055, 0.900, f"Generated {timestamp}.  μ=0.5 for Hubbard; μ=0.875 for Kanamori tiers. No μ scan.", fontsize=9, ha="left", color="#555555")

    table_rows = []
    for beta in BETAS:
        row = [str(beta)]
        for tier in TIERS:
            row.append(fmt(value_at(rows, beta, tier), 3))
        table_rows.append(row)
    table(ax, table_rows, [r"$\beta$", "Hubbard", "density", "spin-flip", "full"], bbox=[0.055, 0.590, 0.89, 0.250], font_size=10.5)
    fig.text(0.055, 0.858, r"Main diagnostic: $|\langle e^{i\theta}\rangle|$ at half filling", fontsize=12.5, weight="bold", ha="left")

    pax = fig.add_axes([0.085, 0.180, 0.540, 0.315])
    for tier in TIERS:
        xs, ys, yerr = [], [], []
        for beta in BETAS:
            y = value_at(rows, beta, tier)
            if y is None:
                continue
            xs.append(beta)
            ys.append(y)
            err = value_at(rows, beta, tier, "average_phase_abs_rank_sem")
            yerr.append(0.0 if err is None else err)
        if xs:
            pax.errorbar(xs, ys, yerr=yerr, marker=MARKERS[tier], color=COLORS[tier], linewidth=2.3, markersize=6, capsize=3, label=TIER_LABELS[tier])
    pax.set_xlabel(r"$\beta$", fontsize=16)
    pax.set_ylabel("avg phase/sign", fontsize=16)
    pax.set_xticks(BETAS)
    pax.set_ylim(0.0, 1.05)
    pax.grid(True, alpha=0.25, linewidth=0.7)
    pax.tick_params(labelsize=13)
    pax.legend(ncol=2, fontsize=11.5, frameon=False, loc="upper right")

    diag_rows = []
    for tier in TIERS:
        completed = sum(1 for beta in BETAS if value_at(rows, beta, tier, "nranks_complete") is not None)
        d10 = value_at(rows, 10, tier, "density_per_cell_reweighted_rank_mean")
        a10 = value_at(rows, 10, tier, "acceptance_mean")
        diag_rows.append([TIER_LABELS[tier], f"{completed}/4", fmt(d10, 4), fmt(a10, 3)])
    table(ax, diag_rows, ["tier", "β rows", "n/cell at β=10", "accept. at β=10"], bbox=[0.675, 0.215, 0.270, 0.245], font_size=8.4, header_color="#fff0cc")

    ax.add_patch(plt.Rectangle((0.675, 0.490), 0.270, 0.060, facecolor="#fff7e6", edgecolor="#d4a64a", linewidth=0.9))
    fig.text(0.690, 0.527, "Interpretation", fontsize=12.5, weight="bold", ha="left", color="#3a2a00")
    fig.text(0.690, 0.505, "Hubbard should remain sign-free;\nKanamori tiers quantify degradation.", fontsize=9.4, ha="left", color="#333333", linespacing=1.15)

    fig.savefig(out, bbox_inches="tight")
    if png:
        fig.savefig(png, bbox_inches="tight")
    print(out.resolve())
    if png:
        print(png.resolve())


if __name__ == "__main__":
    main()

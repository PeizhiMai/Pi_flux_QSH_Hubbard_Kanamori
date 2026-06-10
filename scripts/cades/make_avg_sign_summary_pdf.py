#!/usr/bin/env python3
"""Make a BHZ-style one-page avg-sign summary PDF for pi-flux QSH Kanamori scans."""
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
    "spinflip": "#d62728",
    "full": "#ff7f0e",
}
MARKERS = {"hubbard": "o", "density": "o", "spinflip": "s", "full": "D"}
LINESTYLES = {"hubbard": "-", "density": "-", "spinflip": "--", "full": ":"}


def ffloat(row, key):
    try:
        return float(row[key])
    except Exception:
        return None


def load_rows(path):
    path = Path(path)
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


def fmt_sign(x):
    if x is None:
        return "pending"
    return f"{x:.3f}"


def fmt_diag(x, nd=3):
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


def plot_avg_sign_vs_beta(fig, rows):
    # Match the BHZ summary layout: plot on the lower-left, legend below it.
    pax = fig.add_axes([0.085, 0.225, 0.465, 0.275])
    for tier in TIERS:
        xs, ys, yerr = [], [], []
        for beta in BETAS:
            sign = value_at(rows, beta, tier)
            if sign is not None:
                xs.append(beta)
                ys.append(sign)
                err = value_at(rows, beta, tier, "average_phase_abs_rank_sem")
                yerr.append(0.0 if err is None else err)
        if xs:
            pax.errorbar(
                xs,
                ys,
                yerr=yerr,
                color=COLORS[tier],
                marker=MARKERS[tier],
                linestyle=LINESTYLES[tier],
                linewidth=2.4,
                markersize=6,
                capsize=3,
                label=TIER_LABELS[tier],
            )
    pax.set_xlabel("β", fontsize=18)
    pax.set_ylabel("avg sign", fontsize=18)
    pax.set_xticks(BETAS)
    pax.set_ylim(0.0, 1.05)
    pax.grid(True, alpha=0.25, linewidth=0.7)
    pax.legend(
        ncol=2,
        fontsize=13.6,
        frameon=False,
        loc="upper center",
        bbox_to_anchor=(0.50, -0.200),
        borderaxespad=0.0,
        columnspacing=1.0,
        handlelength=1.45,
        labelspacing=0.35,
    )
    pax.tick_params(labelsize=16)
    return pax


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--summary", default="results/avg_sign_summary/piflux_beta7_10_12_16_half_filling_summary.tsv")
    ap.add_argument("--out", default="docs/piflux_kanamori_avg_sign_summary.pdf")
    ap.add_argument("--png", default="docs/piflux_kanamori_avg_sign_summary.png")
    args = ap.parse_args()

    rows = load_rows(args.summary)
    out = Path(args.out)
    png = Path(args.png) if args.png else None
    out.parent.mkdir(parents=True, exist_ok=True)
    if png:
        png.parent.mkdir(parents=True, exist_ok=True)

    fig = plt.figure(figsize=(11, 8.5), dpi=200)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.axis("off")
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M")

    fig.text(0.055, 0.955, "Pi-flux QSH Hubbard-Kanamori DQMC: average sign summary", fontsize=18, weight="bold", ha="left")
    fig.text(
        0.055,
        0.925,
        "Lx=8, Ly=4, periodic clean system; U=1.0, J$_H$=0.25, t=1.0, λ=0.2, dtau=0.1; Ntherm=2000, Nmeas=10000 per rank, 32 ranks.",
        fontsize=9.5,
        ha="left",
        color="#333333",
    )
    fig.text(
        0.055,
        0.900,
        f"Current summary generated {timestamp}. μ=0.5 for Hubbard; μ=0.875 for Kanamori tiers. No μ scan.",
        fontsize=9,
        ha="left",
        color="#555555",
    )

    half_rows = []
    for beta in BETAS:
        half_rows.append([
            "periodic",
            f"{beta}",
            fmt_sign(value_at(rows, beta, "hubbard")),
            fmt_sign(value_at(rows, beta, "density")),
            fmt_sign(value_at(rows, beta, "spinflip")),
            fmt_sign(value_at(rows, beta, "full")),
        ])
    fig.text(0.055, 0.865, "At physical half filling", fontsize=12.5, weight="bold", ha="left")
    table(
        ax,
        half_rows,
        ["boundary\ncondition", r"$\beta$", "Hubbard\navg sign", "density\navg sign", "spin-flip\navg sign", "full\navg sign"],
        bbox=[0.055, 0.530, 0.89, 0.320],
        font_size=8.2,
    )

    plot_avg_sign_vs_beta(fig, rows)

    n_complete = len(rows)
    status_lines = [
        "Hubbard-only remains\n  sign-free at half filling.",
        "Density Kanamori stays\n  high in completed rows.",
        f"Completed summary rows:\n  {n_complete}/16.",
    ]
    ax.add_patch(plt.Rectangle((0.650, 0.125), 0.295, 0.320, facecolor="#fff7e6", edgecolor="#d4a64a", linewidth=0.9))
    fig.text(0.670, 0.398, "Takeaway", fontsize=22, weight="bold", ha="left", color="#3a2a00")
    y = 0.318
    for line in status_lines:
        fig.text(0.670, y, "• " + line, fontsize=14.8, ha="left", color="#333333", linespacing=1.18)
        y -= 0.102

    fig.savefig(out, bbox_inches="tight")
    if png:
        fig.savefig(png, bbox_inches="tight")
    print(out.resolve())
    if png:
        print(png.resolve())


if __name__ == "__main__":
    main()

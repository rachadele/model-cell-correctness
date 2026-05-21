#!/usr/bin/env python3
"""Combine human + mouse QC coefficient TSVs into one supplementary figure."""

from pathlib import Path
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

BASE = Path("/space/grp/rschwartz/rschwartz/model_cell_correctness/2025-01-30")
INPUTS = {
    ("human", "binary"):     BASE / "homo_sapiens/results_100_cutoff_0_reference_whole cortex/all_cells_model/predicted_doublet/all_binary_coefficients.tsv",
    ("human", "continuous"): BASE / "homo_sapiens/results_100_cutoff_0_reference_whole cortex/all_cells_model/predicted_doublet/all_continuous_coefficients.tsv",
    ("mouse", "binary"):     BASE / "mus_musculus/results_100_cutoff_0_reference_whole cortex/all_cells_model/predicted_doublet/all_binary_coefficients.tsv",
    ("mouse", "continuous"): BASE / "mus_musculus/results_100_cutoff_0_reference_whole cortex/all_cells_model/predicted_doublet/all_continuous_coefficients.tsv",
}

TERM_ORDER = {
    "binary": [
        "outlier_ribo", "outlier_mito", "outlier_hb",
        "counts_outlier", "genes_outlier", "umi_outlier",
        "predicted_doublet",
    ],
    "continuous": [
        "pct_counts_ribo", "pct_counts_mito", "pct_counts_hb",
        "log1p_n_genes_by_counts", "log1p_total_counts",
        "predicted_doublet",
    ],
}

LEVEL_ORDER = ["subclass", "class", "family"]  # human stops at family; mouse "global" excluded for cross-organism comparison
ORG_COLORS = {"human": "#1f77b4", "mouse": "#d62728"}


def load_all():
    rows = []
    for (organism, panel), path in INPUTS.items():
        df = pd.read_csv(path, sep="\t")
        df = df[df["term"].isin(TERM_ORDER[panel])].copy()
        df["organism"] = organism
        df["panel"] = panel
        rows.append(df[["organism", "panel", "key", "term", "Estimate", "2.5_ci", "97.5_ci", "P-val"]])
    return pd.concat(rows, ignore_index=True)


def plot(df, outpath):
    panel = "binary"
    terms = TERM_ORDER[panel]
    y_pos = np.arange(len(terms))
    dodge = 0.18

    fig, axes = plt.subplots(
        nrows=1, ncols=len(LEVEL_ORDER),
        figsize=(3.6 * len(LEVEL_ORDER), 3.6),
        sharey=True, sharex=True,
    )

    for col_idx, level in enumerate(LEVEL_ORDER):
        ax = axes[col_idx]
        ax.axvline(0, color="gray", linestyle="--", linewidth=0.8, zorder=0)

        for organism, offset in [("human", +dodge), ("mouse", -dodge)]:
            sub = df[(df["panel"] == panel) & (df["key"] == level) & (df["organism"] == organism)]
            sub = sub.set_index("term").reindex(terms).reset_index()
            est = sub["Estimate"].to_numpy()
            lo = sub["2.5_ci"].to_numpy()
            hi = sub["97.5_ci"].to_numpy()
            pval = sub["P-val"].to_numpy()

            if np.all(np.isnan(est)):
                continue

            xerr = np.vstack([est - lo, hi - est])
            ax.errorbar(
                est, y_pos + offset,
                xerr=xerr,
                fmt="o", markersize=5,
                color=ORG_COLORS[organism],
                ecolor=ORG_COLORS[organism],
                elinewidth=1.2, capsize=3,
                label=organism if col_idx == 0 else None,
            )

            for y, e, p in zip(y_pos + offset, est, pval):
                if np.isnan(e) or np.isnan(p):
                    continue
                if p < 0.001:
                    star = "***"
                elif p < 0.01:
                    star = "**"
                elif p < 0.05:
                    star = "*"
                else:
                    continue
                ax.text(e, y - 0.05, star, ha="center", va="bottom",
                        fontsize=7, color=ORG_COLORS[organism])

        ax.set_yticks(y_pos)
        ax.set_yticklabels(terms, fontsize=9)
        ax.tick_params(axis="x", labelsize=9)
        ax.invert_yaxis()
        ax.set_title(level, fontsize=11, weight="bold")
        ax.set_xlabel("Log-odds (95% CI)", fontsize=9)

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=2, frameon=False,
               bbox_to_anchor=(0.5, 1.02), fontsize=10)

    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(outpath, dpi=300, bbox_inches="tight")
    fig.savefig(outpath.with_suffix(".pdf"), bbox_inches="tight")
    print(f"Wrote {outpath} and {outpath.with_suffix('.pdf')}")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--outdir", type=Path,
                   default=Path("/space/grp/rschwartz/rschwartz/comp_exam/figures"))
    p.add_argument("--name", default="qc_coefficients_human_vs_mouse")
    args = p.parse_args()

    df = load_all()
    args.outdir.mkdir(parents=True, exist_ok=True)
    plot(df, args.outdir / f"{args.name}.png")

    # also drop the combined tidy table for the figure caption / supplement
    df.to_csv(args.outdir / f"{args.name}.tsv", sep="\t", index=False)
    print(f"Wrote {args.outdir / (args.name + '.tsv')}")


if __name__ == "__main__":
    main()

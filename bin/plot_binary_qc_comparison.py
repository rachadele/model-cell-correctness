#!/usr/bin/env python3
"""Combined human + mouse binary QC coefficient figure (single ref_key)."""

from pathlib import Path
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

PROJECT = Path("/space/grp/rschwartz/rschwartz/model_cell_correctness")
RESULTS_TEMPLATE = (
    "{census}/{organism}/{method}/"
    "results_{subsample}_cutoff_{cutoff}_reference_{reference}/"
    "all_cells_model/predicted_doublet/all_binary_coefficients.tsv"
)

# Order from top to bottom: random-effect SD, intercept, then predictors.
TERM_ORDER = [
    "sd(study)",
    "(Intercept)",
    "ribo_outlier", "mito_outlier", "hb_outlier",
    "counts_outlier", "genes_outlier", "umi_outlier",
    "predicted_doublet",
]
ORG_COLORS = {"human": "#1f77b4", "mouse": "#d62728"}
ORG_DIRS   = {"human": "homo_sapiens", "mouse": "mus_musculus"}


def coef_path(organism, *, census, method, subsample, cutoff, reference):
    return PROJECT / RESULTS_TEMPLATE.format(
        census=census, organism=ORG_DIRS[organism], method=method,
        subsample=subsample, cutoff=cutoff, reference=reference,
    )


def load_all(args):
    """Read one binary-coef TSV per organism, restrict to args.level,
    rename the `study` random-effect row to `sd(study)` with Estimate=sqrt(Var)."""
    rows = []
    for organism in ORG_DIRS:
        path = coef_path(organism, census=args.census, method=args.method,
                         subsample=args.subsample, cutoff=args.cutoff,
                         reference=args.reference)
        if not path.exists():
            print(f"WARNING: missing {path}; skipping {organism}")
            continue
        df = pd.read_csv(path, sep="\t")
        df = df[df["key"] == args.level].copy()
        # Convert the random-effect variance row to SD for plotting
        is_re = df["term"] == "study"
        if is_re.any():
            df.loc[is_re, "Estimate"] = np.sqrt(df.loc[is_re, "Var"])
            df.loc[is_re, ["2.5_ci", "97.5_ci"]] = np.nan
            df.loc[is_re, "term"] = "sd(study)"
        df["organism"] = organism
        rows.append(df[["organism", "term", "Estimate",
                        "2.5_ci", "97.5_ci", "P-val"]])
    if not rows:
        raise FileNotFoundError("No coefficient TSVs found; has the pipeline been run?")
    return pd.concat(rows, ignore_index=True)


def sig_star(p):
    if pd.isna(p):
        return ""
    if p < 0.001: return "***"
    if p < 0.01:  return "**"
    if p < 0.05:  return "*"
    return ""


def plot(df, outpath, level):
    y_pos = np.arange(len(TERM_ORDER))
    dodge = 0.22

    fig, ax = plt.subplots(figsize=(9.0, 6.2))
    ax.axvline(0, color="gray", linestyle="--", linewidth=0.9, zorder=0)

    # Star column anchor: just inside the right spine
    star_pad = 0.04  # fraction of axis width

    for organism, offset in [("human", +dodge), ("mouse", -dodge)]:
        sub = df[df["organism"] == organism].set_index("term").reindex(TERM_ORDER).reset_index()
        est = sub["Estimate"].to_numpy()
        lo  = sub["2.5_ci"].to_numpy()
        hi  = sub["97.5_ci"].to_numpy()
        pval = sub["P-val"].to_numpy()

        if np.all(np.isnan(est)):
            continue

        xerr_lo = np.where(np.isnan(lo), 0, est - lo)
        xerr_hi = np.where(np.isnan(hi), 0, hi - est)
        xerr = np.vstack([xerr_lo, xerr_hi])

        ax.errorbar(
            est, y_pos + offset,
            xerr=xerr,
            fmt="o", markersize=8,
            color=ORG_COLORS[organism],
            ecolor=ORG_COLORS[organism],
            elinewidth=1.6, capsize=4,
            label=organism,
        )

    ax.set_yticks(y_pos)
    ax.set_yticklabels(TERM_ORDER, fontsize=15)
    ax.tick_params(axis="x", labelsize=14)
    ax.invert_yaxis()
    ax.set_title(f"{level}", fontsize=18, weight="bold", loc="right", pad=12)
    ax.set_xlabel("Log-odds (95% CI)", fontsize=15)

    # Place significance stars in a fixed column to the right of the right spine,
    # so they never overlap the bars.
    xlim = ax.get_xlim()
    star_x = xlim[1] + (xlim[1] - xlim[0]) * star_pad
    for organism, offset in [("human", +dodge), ("mouse", -dodge)]:
        sub = df[df["organism"] == organism].set_index("term").reindex(TERM_ORDER).reset_index()
        for i, p in enumerate(sub["P-val"].to_numpy()):
            star = sig_star(p)
            if star:
                ax.text(star_x, i + offset, star, ha="left", va="center",
                        fontsize=14, color=ORG_COLORS[organism],
                        clip_on=False)
    # Expand right limit slightly so stars are visible
    ax.set_xlim(xlim[0], star_x + (xlim[1] - xlim[0]) * 0.06)

    ax.legend(loc="lower left", bbox_to_anchor=(0, 1.02),
              ncol=2, frameon=False, fontsize=14)

    fig.tight_layout()
    fig.savefig(outpath, dpi=300, bbox_inches="tight")
    print(f"Wrote {outpath}")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--census", default="2024-07-01")
    p.add_argument("--method", default="scvi")
    p.add_argument("--subsample", default="100")
    p.add_argument("--cutoff", default="0")
    p.add_argument("--reference", default="whole cortex")
    p.add_argument("--level", default="subclass",
                   help="ref_key to plot (subclass/class/family/global)")
    p.add_argument("--outdir", type=Path,
                   default=Path("/space/grp/rschwartz/rschwartz/model_cell_correctness/figures"))
    p.add_argument("--name", default=None,
                   help="Output filename stem (default: binary_qc_coefficients_human_vs_mouse_<level>)")
    args = p.parse_args()

    name = args.name or f"binary_qc_coefficients_human_vs_mouse_{args.level}"

    df = load_all(args)
    args.outdir.mkdir(parents=True, exist_ok=True)
    plot(df, args.outdir / f"{name}.png", args.level)

    df.to_csv(args.outdir / f"{name}.tsv", sep="\t", index=False)
    print(f"Wrote {args.outdir / (name + '.tsv')}")


if __name__ == "__main__":
    main()

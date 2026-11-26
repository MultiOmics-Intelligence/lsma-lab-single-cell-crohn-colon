#!/usr/bin/env python3
"""The figure for the validation chapter: what happens to the disease answer when the chemistry is controlled.

Usage:  python make_confound_figure.py <validation_dir>
"""
import json
import pathlib
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

NAVY, CORAL, PLUM, GREY = "#0F1B33", "#F06A42", "#8A2D5B", "#56607A"
FAINT = "#E4E7ED"


def main():
    d = pathlib.Path(sys.argv[1])
    summary = pd.read_csv(d / "tables" / "comparison_summary.tsv", sep="\t").set_index("comparison")
    overlap = json.loads((d / "tables" / "overlap.json").read_text())

    order = ["A_crohn_vs_healthy_all_chemistries", "B_healthy_v1_vs_v2_chemistry_only",
             "C_crohn_vs_healthy_v2_only", "D_inflamed_vs_noninflamed_paired"]
    labels = ["A\nCrohn's vs healthy\nevery kit version",
              "B\nhealthy vs healthy\nkit version only",
              "C\nCrohn's vs healthy\nwithin one kit",
              "D\ninflamed vs\nnon-inflamed"]
    sub = ["the naive answer", "no biology at all", "chemistry matched", "paired within patient"]
    colours = [CORAL, PLUM, NAVY, NAVY]
    counts = [int(summary.loc[k, "n_significant"]) if k in summary.index
              and pd.notna(summary.loc[k, "n_significant"]) else 0 for k in order]

    fig = plt.figure(figsize=(14.5, 7.2))
    fig.patch.set_facecolor("white")
    gs = fig.add_gridspec(1, 2, width_ratios=[1.5, 1], wspace=0.22)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])

    bars = ax1.bar(range(4), counts, color=colours, width=0.6, zorder=3)
    for b, c in zip(bars, counts):
        ax1.text(b.get_x() + b.get_width() / 2, c + max(counts) * 0.025, f"{c:,}",
                 ha="center", va="bottom", fontsize=14, fontweight="bold", color=NAVY)
    ax1.set_xticks(range(4))
    ax1.set_xticklabels(labels, fontsize=9.6, color=NAVY, linespacing=1.5)
    for i, s in enumerate(sub):
        ax1.text(i, -max(counts) * 0.175, s, ha="center", va="top", fontsize=8.6,
                 color=PLUM if i == 1 else GREY, style="italic")
    ax1.set_ylim(0, max(counts) * 1.16)
    ax1.set_ylabel("Genes called differentially expressed\npadj < 0.05 and |log2 fold change| > 1",
                   fontsize=10, color=NAVY)
    ax1.set_title("Four ways of asking the same question", fontsize=13.5, color=NAVY, pad=14, loc="left")
    ax1.grid(axis="y", color="#DDE1E8", zorder=0)
    ax1.set_axisbelow(True)
    ax1.tick_params(axis="x", length=0, pad=8)
    for sp in ("top", "right"):
        ax1.spines[sp].set_visible(False)
    ax1.annotate("healthy colon against healthy colon,\nnothing but the kit version changed",
                 xy=(1.3, counts[1]), xytext=(1.65, max(counts) * 0.72), fontsize=9.4, color=PLUM,
                 arrowprops=dict(arrowstyle="->", color=PLUM, lw=1.4, connectionstyle="arc3,rad=-0.2"))

    # right: one bar, broken into what the naive answer is actually made of
    a = overlap.get("significant_in_A_naive", 0)
    both = overlap.get("in_both", 0)
    surv = overlap.get("in_A_and_C", 0)
    rest = max(a - both - surv, 0)
    segs = [(both, PLUM, "also called when\nonly the chemistry\nchanges"),
            (surv, NAVY, "survives matching\nthe chemistry"),
            (rest, FAINT, "neither")]
    left = 0
    for val, col, lab in segs:
        ax2.barh([0], [val], left=left, color=col, height=0.34, zorder=3)
        if val > a * 0.05:
            ax2.text(left + val / 2, 0, f"{val:,}", ha="center", va="center", fontsize=12.5,
                     fontweight="bold", color="white" if col != FAINT else NAVY, zorder=4)
        left += val
    y = -0.33
    for val, col, lab in segs:
        pass
    ax2.set_xlim(0, a * 1.02)
    ax2.set_ylim(-1.05, 0.42)
    ax2.set_yticks([])
    ax2.set_xlabel("Genes in the naive answer (A)", fontsize=10, color=NAVY)
    ax2.set_title(f"What those {a:,} genes are made of", fontsize=13.5, color=NAVY, pad=14, loc="left")
    ax2.grid(axis="x", color="#DDE1E8", zorder=0)
    ax2.set_axisbelow(True)
    for sp in ("top", "right", "left"):
        ax2.spines[sp].set_visible(False)

    # a small legend under the bar, in the segments' own colours
    xs = [0.0, 0.36, 0.72]
    for (val, col, lab), x in zip(segs, xs):
        ax2.add_patch(plt.Rectangle((x * a, -0.40), a * 0.032, 0.075, color=col, clip_on=False))
        ax2.text(x * a + a * 0.05, -0.355, lab, fontsize=8.8, color=NAVY, va="top", linespacing=1.45)

    pct_chem = overlap.get("share_of_A_also_called_by_chemistry_alone")
    pct_surv = overlap.get("share_of_A_surviving_in_C")
    ax2.text(0, -0.80, f"{pct_chem:.0%} of the naive answer can be produced by the chemistry alone.\n"
                       f"Only {pct_surv:.1%} of it survives once the chemistry is matched.",
             fontsize=10.2, color=PLUM, va="top", fontweight="bold", linespacing=1.5)

    fig.suptitle("The 10x chemistry is confounded with the disease, and it carries most of the answer",
                 fontsize=15, color=NAVY, x=0.011, ha="left", y=0.975)
    fig.text(0.011, 0.012, "LSMA Analysis Lab  ·  pseudobulk DESeq2, biosamples as replicates  ·  "
                           "educational and research material, not medical advice", fontsize=8.5, color=GREY)
    fig.subplots_adjust(left=0.065, right=0.985, top=0.875, bottom=0.175)
    for ext in ("png", "pdf"):
        fig.savefig(d / f"validation-chemistry-confound.{ext}", dpi=190 if ext == "png" else None,
                    facecolor="white", bbox_inches="tight")
    print("wrote validation-chemistry-confound.png / .pdf to", d)


if __name__ == "__main__":
    main()

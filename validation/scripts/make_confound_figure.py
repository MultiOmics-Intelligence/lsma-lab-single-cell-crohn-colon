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

NAVY, CORAL, PLUM, GREY, CREAM = "#0F1B33", "#F06A42", "#8A2D5B", "#56607A", "#F6F2EA"


def main():
    d = pathlib.Path(sys.argv[1])
    summary = pd.read_csv(d / "tables" / "comparison_summary.tsv", sep="\t")
    overlap = json.loads((d / "tables" / "overlap.json").read_text())

    order = ["A_crohn_vs_healthy_all_chemistries", "B_healthy_v1_vs_v2_chemistry_only",
             "C_crohn_vs_healthy_v2_only", "D_inflamed_vs_noninflamed_paired"]
    labels = ["A. Crohn's vs healthy\nevery chemistry\n(the naive answer)",
              "B. healthy vs healthy\ndifferent chemistry only\n(no biology at all)",
              "C. Crohn's vs healthy\nwithin one chemistry\n(matched)",
              "D. inflamed vs non-inflamed\nsame patients\n(paired)"]
    colours = [CORAL, PLUM, NAVY, NAVY]
    s = summary.set_index("comparison")
    counts = [int(s.loc[k, "n_significant"]) if k in s.index and pd.notna(s.loc[k, "n_significant"]) else 0
              for k in order]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.5, 6.2), gridspec_kw={"width_ratios": [1.45, 1]})
    fig.patch.set_facecolor("white")

    bars = ax1.bar(range(4), counts, color=colours, width=0.62, zorder=3)
    for i, (b, c) in enumerate(zip(bars, counts)):
        ax1.text(b.get_x() + b.get_width() / 2, c + max(counts) * 0.022, f"{c:,}",
                 ha="center", va="bottom", fontsize=13, fontweight="bold", color=NAVY)
    ax1.set_xticks(range(4))
    ax1.set_xticklabels(labels, fontsize=9.2, color=NAVY)
    ax1.set_ylabel("Genes called differentially expressed\n(padj < 0.05 and |log2 fold change| > 1)", fontsize=10, color=NAVY)
    ax1.set_title("Four ways of asking the same question", fontsize=13, color=NAVY, pad=12, loc="left")
    ax1.grid(axis="y", color="#DDE1E8", zorder=0)
    ax1.set_axisbelow(True)
    for sp in ("top", "right"):
        ax1.spines[sp].set_visible(False)
    ax1.annotate("comparing healthy colon\nwith healthy colon still\ngives 2,464 genes",
                 xy=(1, counts[1]), xytext=(1.42, max(counts) * 0.76), fontsize=9.2, color=PLUM,
                 arrowprops=dict(arrowstyle="->", color=PLUM, lw=1.3))

    a = overlap.get("significant_in_A_naive", 0)
    both = overlap.get("in_both", 0)
    surv = overlap.get("in_A_and_C", 0)
    parts = [both, surv, max(a - both - surv, 0)]
    ax2.barh([0], [a], color="#E4E7ED", height=0.52, zorder=2)
    left = 0
    for val, col, lab in ((both, PLUM, "also called by\nchemistry alone"),
                          (surv, NAVY, "survives matching\nthe chemistry")):
        ax2.barh([1], [val], left=left, color=col, height=0.52, zorder=3)
        if val:
            ax2.text(left + val / 2, 1, f"{val:,}", ha="center", va="center",
                     fontsize=11, fontweight="bold", color="white", zorder=4)
        left += val
    ax2.barh([1], [max(a - left, 0)], left=left, color="#E4E7ED", height=0.52, zorder=3)
    ax2.set_yticks([0, 1])
    ax2.set_yticklabels([f"all {a:,} genes from A", "of which"], fontsize=10, color=NAVY)
    ax2.set_xlabel("Genes", fontsize=10, color=NAVY)
    ax2.set_title("What the naive answer is made of", fontsize=13, color=NAVY, pad=12, loc="left")
    ax2.grid(axis="x", color="#DDE1E8", zorder=0)
    ax2.set_axisbelow(True)
    for sp in ("top", "right", "left"):
        ax2.spines[sp].set_visible(False)
    pct_chem = overlap.get("share_of_A_also_called_by_chemistry_alone")
    pct_surv = overlap.get("share_of_A_surviving_in_C")
    note = (f"{pct_chem:.0%} of the naive answer can be produced by the chemistry alone.\n"
            f"Only {pct_surv:.1%} of it survives when the chemistry is matched.")
    ax2.text(0, -0.72, note, fontsize=10, color=PLUM, va="top")

    fig.suptitle("The 10x chemistry is confounded with the disease, and it carries most of the answer",
                 fontsize=14.5, color=NAVY, x=0.012, ha="left", y=0.985)
    fig.text(0.012, 0.015, "LSMA Analysis Lab · pseudobulk DESeq2, biosamples as replicates · educational material, not medical advice",
             fontsize=8.4, color=GREY)
    fig.tight_layout(rect=[0, 0.045, 1, 0.93])
    for ext in ("png", "pdf"):
        fig.savefig(d / f"validation-chemistry-confound.{ext}", dpi=190 if ext == "png" else None,
                    facecolor="white", bbox_inches="tight")
    print("wrote validation-chemistry-confound.png / .pdf to", d)


if __name__ == "__main__":
    main()

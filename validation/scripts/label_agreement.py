#!/usr/bin/env python3
"""
Does the lab's automated cell-type annotation agree with the people who made the data?

The dataset carries two independent sets of labels that the lab never saw:
  Celltype   - the authors' own 25 manual labels
  cell_type  - CZ CELLxGENE's harmonised 17-label ontology

The lab scored marker panels against its own Leiden clusters and produced its own
labels. This script compares the three, outside the lab, and reports where they
agree and where they do not. Disagreement is not automatically the lab being
wrong: the panels are coarser than the authors' labels by design, and this
measures that cost rather than assuming it.

Usage:  python label_agreement.py <annotated_adata.h5ad> <outdir>
"""
import json
import pathlib
import sys

import anndata as ad
import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score


def pick(obs, *names):
    for n in names:
        if n in obs.columns:
            return n
    return None


def agreement(a, b, name_a, name_b):
    """Rand index and mutual information between two labellings, plus the best match for each label of a."""
    ari = adjusted_rand_score(a, b)
    nmi = normalized_mutual_info_score(a, b)
    tab = pd.crosstab(pd.Series(a, name=name_a), pd.Series(b, name=name_b))
    rows = []
    for lab, row in tab.iterrows():
        tot = row.sum()
        best = row.idxmax()
        rows.append({name_a: lab, "n_cells": int(tot), f"most_common_{name_b}": best,
                     "share_in_that_label": round(float(row.max()) / tot, 4) if tot else None,
                     "n_distinct_labels_it_maps_to": int((row > 0).sum())})
    per_label = pd.DataFrame(rows).sort_values("n_cells", ascending=False)
    # the share of cells sitting in the modal pairing of their own label
    purity = float(sum(tab.loc[r].max() for r in tab.index) / tab.values.sum())
    return {"adjusted_rand_index": round(ari, 4), "normalised_mutual_information": round(nmi, 4),
            "mean_purity": round(purity, 4)}, per_label, tab


def main():
    h5, out = sys.argv[1], pathlib.Path(sys.argv[2])
    (out / "tables").mkdir(parents=True, exist_ok=True)
    adata = ad.read_h5ad(h5, backed="r")   # metadata only; the expression matrix is never loaded
    obs = adata.obs
    print("cells", adata.n_obs, flush=True)
    print("obs columns:", list(obs.columns), flush=True)

    lab_col = pick(obs, "cell_type_predicted", "cell_type_lab", "annotation", "cell_type_auto")
    cluster_col = pick(obs, "leiden", "cluster", "clusters", "louvain")
    author_col = pick(obs, "Celltype")
    cxg_col = pick(obs, "cell_type")

    found = {"lab_label": lab_col, "cluster": cluster_col, "authors": author_col, "cellxgene": cxg_col}
    print("columns used:", found, flush=True)
    (out / "tables" / "columns_used.json").write_text(json.dumps(found, indent=1))
    if author_col is None:
        print("no authors' labels in this object - nothing to compare", flush=True)
        return

    authors = obs[author_col].astype(str).values
    summary = {}

    for name, col in (("lab_annotation", lab_col), ("leiden_cluster", cluster_col), ("cellxgene_ontology", cxg_col)):
        if col is None:
            print(f"{name}: column not present - skipped", flush=True)
            summary[f"authors_vs_{name}"] = {"status": "not available"}
            continue
        other = obs[col].astype(str).values
        metrics, per_label, tab = agreement(authors, other, "authors_label", name)
        metrics["status"] = "compared"
        metrics["n_labels_authors"] = int(len(set(authors)))
        metrics[f"n_labels_{name}"] = int(len(set(other)))
        summary[f"authors_vs_{name}"] = metrics
        per_label.to_csv(out / "tables" / f"authors_vs_{name}_by_label.tsv", sep="\t", index=False)
        tab.to_csv(out / "tables" / f"authors_vs_{name}_crosstab.tsv", sep="\t")
        print(f"{name}: {metrics}", flush=True)

    # how many of the authors' labels does each lab label absorb, and vice versa
    if lab_col is not None:
        lab = obs[lab_col].astype(str).values
        comp = (pd.DataFrame({"authors": authors, "lab": lab})
                .groupby("lab").agg(n_cells=("authors", "size"),
                                    n_author_labels_absorbed=("authors", "nunique"))
                .sort_values("n_cells", ascending=False))
        comp.to_csv(out / "tables" / "lab_label_composition.tsv", sep="\t")
        print(comp.head(20), flush=True)

    (out / "tables" / "agreement_summary.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary, indent=1), flush=True)
    print("done ->", out, flush=True)


if __name__ == "__main__":
    main()

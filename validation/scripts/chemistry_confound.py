#!/usr/bin/env python3
"""
Is the Crohn's-versus-healthy difference biology, or is it the 10x chemistry?

In this dataset the two are tangled: 10x 3' v1 was used only on healthy colon and
v3 only on Crohn's colon, so a naive comparison of the two conditions is also a
comparison of two chemistries. Only v2 was used on both.

This script runs four pseudobulk DESeq2 comparisons outside the lab, with
biosamples as replicates, and lets them check each other:

  A  Crohn's (inflamed + non-inflamed) vs healthy, every chemistry   - the naive answer
  B  healthy v1 vs healthy v2                                       - chemistry only, no biology
  C  Crohn's vs healthy, restricted to v2                           - chemistry-matched
  D  inflamed vs non-inflamed, in the donors who gave both          - paired, within patient

B is the negative control that matters: both sides are healthy colon, so every
gene it calls significant is a gene the chemistry alone can produce.

Usage:  python chemistry_confound.py <adata.h5ad> <outdir>
"""
import json
import pathlib
import sys

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

PADJ = 0.05
LFC = 1.0
MIN_SAMPLE_CELLS = 50      # a biosample needs this many cells to be a replicate
MIN_REPLICATES = 3         # a side needs this many biosamples


def pseudobulk(adata, sample_col="biosample_id"):
    """Sum raw counts over the cells of each biosample. One column per gene, one row per sample."""
    counts = adata.layers["counts"] if "counts" in adata.layers else adata.X
    counts = sp.csr_matrix(counts)
    samples = adata.obs[sample_col].astype(str).values
    uniq = pd.unique(samples)
    rows, index = [], []
    for s in uniq:
        m = samples == s
        if m.sum() < MIN_SAMPLE_CELLS:
            continue
        rows.append(np.asarray(counts[m].sum(axis=0)).ravel())
        index.append(s)
    mat = pd.DataFrame(np.vstack(rows), index=index, columns=adata.var_names.astype(str))
    # one row of metadata per surviving sample, taken from its first cell
    meta = (adata.obs.assign(_s=samples)
            .drop_duplicates("_s").set_index("_s")
            .loc[index, ["Type", "assay", "donor_id"]])
    meta["n_cells"] = [int((samples == s).sum()) for s in index]
    return mat.round().astype(int), meta


def deseq(mat, meta, design, contrast, label):
    """PyDESeq2 on pseudobulk counts. Returns the result table, or None if the design cannot be fitted."""
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    keep = mat.sum(axis=0) > 0
    mat = mat.loc[:, keep]
    # genes seen in at least a quarter of the samples, to keep the fit stable
    detected = (mat > 0).sum(axis=0) >= max(2, int(0.25 * mat.shape[0]))
    mat = mat.loc[:, detected]
    dds = DeseqDataSet(counts=mat, metadata=meta, design=design, refit_cooks=True, quiet=True)  # design is a formula string
    dds.deseq2()
    st = DeseqStats(dds, contrast=contrast, quiet=True)
    st.summary()
    res = st.results_df.copy()
    res["gene"] = res.index
    res["comparison"] = label
    return res.sort_values("padj")


def significant(res):
    if res is None:
        return pd.DataFrame()
    r = res.dropna(subset=["padj"])
    return r[(r["padj"] < PADJ) & (r["log2FoldChange"].abs() > LFC)]


def main():
    h5, out = sys.argv[1], pathlib.Path(sys.argv[2])
    (out / "tables").mkdir(parents=True, exist_ok=True)
    print("reading", h5, flush=True)
    adata = ad.read_h5ad(h5)
    print("cells", adata.n_obs, "genes", adata.n_vars, flush=True)

    mat, meta = pseudobulk(adata)
    meta["crohn"] = np.where(meta["Type"].isin(["Infl", "NonI"]), "crohn", "healthy")
    print("pseudobulk samples", mat.shape[0], "genes", mat.shape[1], flush=True)
    print(meta.groupby(["crohn", "assay"], observed=True).size(), flush=True)
    meta.to_csv(out / "tables" / "pseudobulk_samples.tsv", sep="\t")

    summary = []
    results = {}

    def run(label, sub_meta, design, contrast, note):
        sides = sub_meta[contrast[0]].astype(str).value_counts()
        if len(sides) < 2 or sides.min() < MIN_REPLICATES:
            print(f"{label}: not enough replicates {dict(sides)} - skipped", flush=True)
            summary.append({"comparison": label, "note": note, "status": "not run",
                            "reason": f"fewer than {MIN_REPLICATES} samples on a side: {dict(sides)}",
                            "n_samples": int(len(sub_meta)), "n_significant": None})
            return
        print(f"--- {label}: {dict(sides)}", flush=True)
        res = deseq(mat.loc[sub_meta.index], sub_meta, design, contrast, label)
        results[label] = res
        res.to_csv(out / "tables" / f"de_{label}.tsv", sep="\t", index=False)
        sig = significant(res)
        summary.append({"comparison": label, "note": note, "status": "run",
                        "reason": "", "n_samples": int(len(sub_meta)),
                        "n_genes_tested": int(res["padj"].notna().sum()),
                        "n_significant": int(len(sig)),
                        "sides": json.dumps({str(k): int(v) for k, v in sides.items()})})
        print(f"    {len(sig)} significant of {res['padj'].notna().sum()} tested", flush=True)

    # A - the naive comparison, every chemistry
    run("A_crohn_vs_healthy_all_chemistries", meta.copy(),
        "~crohn", ["crohn", "crohn", "healthy"],
        "Crohn's vs healthy using every sample. Chemistry is not matched between the groups.")

    # B - the negative control: healthy only, v1 against v2
    hb = meta[(meta["crohn"] == "healthy") & (meta["assay"].astype(str).isin(["10x 3' v1", "10x 3' v2"]))].copy()
    hb["chem"] = np.where(hb["assay"].astype(str) == "10x 3' v1", "v1", "v2")
    run("B_healthy_v1_vs_v2_chemistry_only", hb,
        "~chem", ["chem", "v1", "v2"],
        "Healthy colon against healthy colon, differing only in 10x chemistry. Any gene called here is one chemistry alone can produce.")

    # C - chemistry-matched: v2 only
    c = meta[meta["assay"].astype(str) == "10x 3' v2"].copy()
    run("C_crohn_vs_healthy_v2_only", c,
        "~crohn", ["crohn", "crohn", "healthy"],
        "Crohn's vs healthy within one chemistry, so the comparison is not confounded by it.")

    # D - paired within patient: donors who gave both an inflamed and a non-inflamed sample
    both = meta[meta["Type"].isin(["Infl", "NonI"])]
    paired_donors = [d for d, g in both.groupby("donor_id", observed=True)
                     if {"Infl", "NonI"} <= set(g["Type"].astype(str))]
    d = both[both["donor_id"].isin(paired_donors)].copy()
    d["Type"] = d["Type"].astype(str)
    d["donor_id"] = d["donor_id"].astype(str)
    print(f"paired donors: {paired_donors}", flush=True)
    if len(paired_donors) >= 2 and d["Type"].value_counts().min() >= 2:
        run("D_inflamed_vs_noninflamed_paired", d,
            "~donor_id + Type", ["Type", "Infl", "NonI"],
            "Inflamed against non-inflamed tissue from the same patients, with the donor in the model. Immune to the chemistry confound and to differences between people.")
    else:
        summary.append({"comparison": "D_inflamed_vs_noninflamed_paired",
                        "note": "Inflamed vs non-inflamed, paired within donor.",
                        "status": "not run",
                        "reason": f"only {len(paired_donors)} donors gave both tissue types with enough cells",
                        "n_samples": int(len(d)), "n_significant": None})

    pd.DataFrame(summary).to_csv(out / "tables" / "comparison_summary.tsv", sep="\t", index=False)

    # how much of the naive answer is chemistry?
    overlap = {}
    a, b = results.get("A_crohn_vs_healthy_all_chemistries"), results.get("B_healthy_v1_vs_v2_chemistry_only")
    if a is not None and b is not None:
        sa, sbn = set(significant(a)["gene"]), set(significant(b)["gene"])
        shared = sa & sbn
        overlap = {
            "significant_in_A_naive": len(sa),
            "significant_in_B_chemistry_only": len(sbn),
            "in_both": len(shared),
            "share_of_A_also_called_by_chemistry_alone": round(len(shared) / len(sa), 4) if sa else None,
        }
        c2 = results.get("C_crohn_vs_healthy_v2_only")
        if c2 is not None:
            sc = set(significant(c2)["gene"])
            overlap["significant_in_C_chemistry_matched"] = len(sc)
            overlap["in_A_and_C"] = len(sa & sc)
            overlap["share_of_A_surviving_in_C"] = round(len(sa & sc) / len(sa), 4) if sa else None
        pd.Series(sorted(shared)).to_csv(out / "tables" / "genes_in_both_A_and_chemistry_control.tsv",
                                         sep="\t", index=False, header=["gene"])
    (out / "tables" / "overlap.json").write_text(json.dumps(overlap, indent=1))
    print(json.dumps(overlap, indent=1), flush=True)
    print("done ->", out, flush=True)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Cross-tabulate WGS and array genotypes of the rare driver variants typed directly on an array, by case status.
People: unrelated full mainland with WGS and the array in question. Counts only."""
import glob
import json
import os

import pandas as pd

HERE = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/check/array_crosscheck"
PHENO = ("/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/wgs.auto.par/results/12_model_inputs/"
         "full_mainland/phenotype/pheno.tsv")
O = f"{HERE}/results/05.rare_direct"
LABEL = {"chr12:112890858": "RPH3A p.Arg544Cys", "chr11:59576952": "OSBP p.Pro712Ser",
         "chr9:5073770": "JAK2 p.Val617Phe"}


def raw(path):
    r = pd.read_csv(path, sep=r"\s+").drop(columns=["FID", "PAT", "MAT", "SEX", "PHENOTYPE"])
    r["IID"] = r.IID.astype(str)
    return r.set_index("IID")


ph = pd.read_csv(PHENO, sep="\t", dtype={"IID": str}).set_index("IID").PHENO1
W = raw(f"{O}/wgs.raw")
rows = []
for arr in ("illuminaASA24v1", "illuminaJSA24v1"):
    p = f"{O}/{arr}.raw"
    if not os.path.exists(p):
        continue
    R = raw(p)
    common = R.index.intersection(W.index).intersection(ph.index)
    for ca in R.columns:
        aid, a_counted = ca.rsplit("_", 1)
        pos = ":".join(aid.split(":")[:2])
        ref, alt = aid.split(":")[2:4]
        for cw in [c for c in W.columns if c.startswith(pos + ":")]:
            wid, w_counted = cw.rsplit("_", 1)
            if set(wid.split(":")[2:4]) != {ref, alt}:
                continue
            a = R.loc[common, ca] if a_counted == alt else 2 - R.loc[common, ca]
            w = W.loc[common, cw] if w_counted == alt else 2 - W.loc[common, cw]
            for g, lab in ((1, "cases"), (0, "controls")):
                idx = ph.loc[common][ph.loc[common] == g].index
                tab = pd.crosstab(w.loc[idx].fillna(-1).astype(int), a.loc[idx].fillna(-1).astype(int))
                rows.append(dict(variant=LABEL.get(pos, pos), array_id=aid, array=arr[8:], group=lab, n=len(idx),
                                 wgs_carriers=int((w.loc[idx] >= 1).sum()),
                                 array_carriers=int((a.loc[idx] >= 1).sum()),
                                 both_carriers=int(((w.loc[idx] >= 1) & (a.loc[idx] >= 1)).sum()),
                                 wgs_only=int(((w.loc[idx] >= 1) & (a.loc[idx] == 0)).sum()),
                                 array_only=int(((w.loc[idx] == 0) & (a.loc[idx] >= 1)).sum()),
                                 missing_wgs=int(w.loc[idx].isna().sum()), missing_array=int(a.loc[idx].isna().sum()),
                                 table="; ".join(f"WGS {i}/array {j}: {tab.loc[i, j]}" for i in tab.index
                                                 for j in tab.columns if tab.loc[i, j])))
T = pd.DataFrame(rows)
T.to_csv(f"{O}/rare_direct.tsv", sep="\t", index=False)
print(T.drop(columns=["table"]).to_string(index=False))
print(T[["variant", "array", "group", "table"]].to_string(index=False))

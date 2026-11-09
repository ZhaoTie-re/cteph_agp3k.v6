#!/usr/bin/env python3
"""Step 6: WGS vs. imputed-array genotypes in the same 1,683 people. The panel is AGP3K whole genomes (the
controls' own sequences), so the cases' imputed genotypes are the part independent of the WGS calls of the cases.

(a) Panel membership: of the WGS carriers of rare variants near the drivers, the share imputed as carriers
    (dosage >= 0.5), by case status and by the number of control carriers. Near-complete recovery in controls and
    not in cases means the controls are in the panel. Also the share of WGS rare variants present in the imputed set
    at all, by whether their carriers are cases only, controls only, or both.
(b) Rare drivers: WGS vs. imputed carriers by case status, with the case recovery expected for variants of the same
    control-carrier count. Counts only.
(c) Common-variant leads: OR on imputed dosages vs. WGS genotypes, allele frequency in cases on both, and the
    dosage-genotype r^2 in cases and in controls."""
import json

import numpy as np
import pandas as pd

HERE = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/check/array_crosscheck"
PHENO = ("/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/wgs.auto.par/results/12_model_inputs/"
         "full_mainland/phenotype/pheno.tsv")
O = f"{HERE}/results/06.imputed"
DRIVERS = {"chr4:76306990:G:C": "STBD1 c.220+1G>C", "chr4:76309559:G:T": "STBD1 p.Trp212Cys",
           "chr4:76309487:C:G": "STBD1 p.Ser188Arg", "chr12:112890858:C:T": "RPH3A p.Arg544Cys",
           "chr11:59576952:G:A": "OSBP p.Pro712Ser", "chr11:59594182:C:T": "OSBP p.Arg462Gln",
           "chr9:5073770:G:T": "JAK2 p.Val617Phe", "chr9:5073770:G:A": "JAK2 c.1849G>A (other allele)"}
BINS = [(0, 0, "0"), (1, 1, "1"), (2, 3, "2-3"), (4, 10, "4-10"), (11, 10**6, ">=11")]


def alt_counts(path):
    """A .raw file as ALT-allele counts or dosages, columns = variant IDs (chr:pos:ref:alt)."""
    r = pd.read_csv(path, sep="\t").drop(columns=["FID", "PAT", "MAT", "SEX", "PHENOTYPE"])
    r["IID"] = r.IID.astype(str)
    r = r.set_index("IID")
    out = {}
    for c in r.columns:
        vid, counted = c.rsplit("_", 1)
        out[vid] = r[c] if counted == vid.split(":")[3] else 2 - r[c]
    return pd.DataFrame(out)


def info(path):
    t = pd.read_csv(path, sep="\t", header=None, names=["id", "info"]).set_index("id")["info"]
    get = lambda s, k: next((float(x.split("=")[1]) for x in s.split(";") if x.startswith(k + "=")), np.nan)
    return pd.DataFrame({"r2": t.map(lambda s: get(s, "R2")), "typed": t.map(lambda s: "TYPED" in s.split(";"))})


def binlab(n):
    return next(l for lo, hi, l in BINS if lo <= n <= hi)


ph = pd.read_csv(PHENO, sep="\t", dtype={"IID": str}).set_index("IID").PHENO1
W = alt_counts(f"{O}/wgs.raw")
D = alt_counts(f"{O}/imp.raw")
ids = W.index.intersection(D.index).intersection(ph.index)
W, D, ph = W.loc[ids], D.loc[ids], ph.loc[ids]
case, ctrl = (ph == 1).values, (ph == 0).values
inf = pd.concat([info(f"{O}/imp.drivers_window.info"), info(f"{O}/imp.leads.info")])
inf = inf[~inf.index.duplicated()]
shared = [v for v in W.columns if v in D.columns]

# (a) carrier recovery of rare variants near the drivers (imputed, not typed, sites only)
rare = [v for v in open(f"{O}/rare_shared.ids").read().split() if v in shared and not inf.typed.get(v, True)]
rows = []
for v in rare:
    g, d = W[v].values, D[v].values
    if np.nanmean(g) / 2 > 0.5:                        # count the minor allele
        g, d = 2 - g, 2 - d
    cc = int(np.nansum(g[ctrl] >= 1))
    for grp, m in (("cases", case), ("controls", ctrl)):
        car, non = m & (g >= 1), m & (g == 0)
        rows.append(dict(variant=v, group=grp, control_carriers_wgs=cc, bin=binlab(cc), wgs_carriers=int(car.sum()),
                         recovered=int((car & (d >= 0.5)).sum()), wgs_noncarriers=int(non.sum()),
                         imputed_only=int((non & (d >= 0.5)).sum())))
R = pd.DataFrame(rows)
rec = (R.groupby(["group", "bin"], sort=False)
       .agg(variants=("variant", lambda s: int((R.loc[s.index, "wgs_carriers"] > 0).sum())),
            wgs_carriers=("wgs_carriers", "sum"), recovered=("recovered", "sum"),
            wgs_noncarriers=("wgs_noncarriers", "sum"), imputed_only=("imputed_only", "sum")).reset_index())
rec["recovery"] = rec.recovered / rec.wgs_carriers.where(rec.wgs_carriers > 0)
rec["false_carrier_rate"] = rec.imputed_only / rec.wgs_noncarriers
rec["bin"] = pd.Categorical(rec.bin, [b[2] for b in BINS], ordered=True)
rec = rec.sort_values(["group", "bin"])
rec.to_csv(f"{O}/recovery_by_group.tsv", sep="\t", index=False)
expected = rec[rec.group == "cases"].astype({"bin": str}).set_index("bin").recovery

# presence in the imputed set of WGS rare variants, by who carries them (all rare WGS variants in the windows)
cnt = {}
for G in ("cases", "controls"):
    a = pd.read_csv(f"{O}/wgs.drv.{G}.acount", sep="\t").set_index("ID")
    cnt[G] = a
a_case = cnt["cases"]
a_ctrl = cnt["controls"].loc[a_case.index]
alt_major = (a_case.ALT_CTS + a_ctrl.ALT_CTS) > (a_case.OBS_CT + a_ctrl.OBS_CT) / 2
mc_case = np.where(alt_major, a_case.OBS_CT - a_case.ALT_CTS, a_case.ALT_CTS)
mc_ctrl = np.where(alt_major, a_ctrl.OBS_CT - a_ctrl.ALT_CTS, a_ctrl.ALT_CTS)
pres = pd.DataFrame({"mc_case": mc_case, "mc_ctrl": mc_ctrl}, index=a_case.index)
pres = pres[(pres.mc_case + pres.mc_ctrl).between(1, 30)]
imp_ids = set(inf.index)
pres["present"] = [v in imp_ids for v in pres.index]
pres["carriers"] = np.select([pres.mc_ctrl == 0, pres.mc_case == 0], ["cases only", "controls only"], "both")
presence = (pres.groupby("carriers").agg(variants=("present", "size"), in_imputed_set=("present", "sum"))
            .reset_index())
presence["share"] = presence.in_imputed_set / presence.variants
presence.to_csv(f"{O}/presence_by_carriers.tsv", sep="\t", index=False)

# (b) rare drivers, counts only
drv = []
for v, lab in DRIVERS.items():
    if v not in W.columns:
        continue
    g = W[v]
    d = D[v] if v in D.columns else None
    cc = int((g[ctrl] >= 1).sum())
    for grp, m in (("cases", case), ("controls", ctrl)):
        car = m & (g.values >= 1)
        row = dict(variant=lab, id=v, in_imputed_set=d is not None, typed=bool(inf.typed.get(v, False)),
                   r2=inf.r2.get(v, np.nan), group=grp, n=int(m.sum()), wgs_carriers=int(car.sum()),
                   wgs_missing=int(np.isnan(g.values[m]).sum()))
        if d is not None:
            dv = d.values
            row.update(imputed_carriers=int((m & (dv >= 0.5)).sum()), both=int((car & (dv >= 0.5)).sum()),
                       wgs_only=int((car & (dv < 0.5)).sum()),
                       imputed_only=int((m & (g.values == 0) & (dv >= 0.5)).sum()),
                       mean_dosage_wgs_carriers=round(float(dv[car].mean()), 3) if car.any() else np.nan,
                       imputed_dosage_sum=round(float(dv[m].sum()), 2))
        if grp == "cases":
            row["expected_case_recovery"] = expected.get(binlab(cc), np.nan)
        drv.append(row)
DRV = pd.DataFrame(drv)
DRV.to_csv(f"{O}/drivers.tsv", sep="\t", index=False)

# (c) the common-variant leads
sig = pd.read_csv(f"{HERE}/results/04.compare/signals.tsv", sep="\t")
sig = sig[sig.lead_source == "wgs"]
gl = pd.read_csv(f"{O}/imp.leads.PHENO1.glm.logistic.hybrid", sep="\t").set_index("ID")
lead_rows = []
for r in sig.itertuples():
    v = r.lead
    if v not in gl.index:
        lead_rows.append(dict(signal=r.signal, lead=v, rsid=r.rsid, in_imputed_set=False))
        continue
    e = gl.loc[v]
    a1 = str(r.lead_a1)
    or_imp = e.OR if e.A1 == a1 else 1 / e.OR
    ref, alt = v.split(":")[2:4]
    g, d = W[v], D[v]
    sgn = 1 if a1 == alt else -1                       # frequencies of the lead's A1
    fq = lambda x, m: float((x[m] if sgn == 1 else 2 - x[m]).mean() / 2)
    r2 = lambda m: float(np.corrcoef(g.values[m][~np.isnan(g.values[m])], d.values[m][~np.isnan(g.values[m])])[0, 1] ** 2)
    lead_rows.append(dict(signal=r.signal, lead=v, rsid=r.rsid, a1=a1, in_imputed_set=True,
                          typed=bool(inf.typed.get(v, False)), info_r2=inf.r2.get(v, np.nan),
                          or_wgs=r.sub_wgs_or_lead, p_wgs=r.sub_wgs_p_lead, or_imputed=or_imp, p_imputed=e.P,
                          freq_cases_wgs=fq(g, case), freq_cases_imputed=fq(d, case),
                          freq_controls_wgs=fq(g, ctrl), freq_controls_imputed=fq(d, ctrl),
                          r2_cases=r2(case), r2_controls=r2(ctrl), proxy_r2_on_array=r.proxy_r2))
L = pd.DataFrame(lead_rows)
L.to_csv(f"{O}/leads.tsv", sep="\t", index=False)

summary = dict(people=int(len(ids)), cases=int(case.sum()), controls=int(ctrl.sum()), rare_imputed_sites=len(rare),
               recovery=rec.round(4).astype({"bin": str}).to_dict("records"),
               presence=presence.round(4).to_dict("records"))
json.dump(summary, open(f"{O}/summary.json", "w"), indent=1, default=float)
pd.set_option("display.width", 250)
print(rec.to_string(index=False))
print(presence.to_string(index=False))
print(DRV.to_string(index=False))
print(L.round(4).to_string(index=False))

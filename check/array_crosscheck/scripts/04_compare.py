#!/usr/bin/env python3
"""Step 4: do the association results depend on the WGS genotypes?

The same people (unrelated full-mainland genomes with SNP-array genotypes), the same model (logistic, sex + 10
BBJ-mainland PCs), two technologies. Cases are DNBSeq-T7 genomes typed on ASA24v1; controls are HiSeq X genomes
typed on JSA24v1. A WGS calling artifact that separates cases from controls would not be reproduced by arrays,
and an array artifact would not be reproduced by WGS; agreement between the two argues against both.

1. Genome-wide: QC'd shared sites; z(WGS) vs z(array); lambda of each; genotype discordance in cases vs controls.
2. Signals: each lead, its best array proxy (r2 from WGS genotypes in controls), and the association at the proxy
   on both technologies.
3. HLA-DQB1 57 V: dosage predicted from array genotypes alone (ridge regression trained on the WGS typing of
   controls; out-of-fold for controls), compared with the WGS typing in cases, and tested for association;
   frequencies of its best array tags against ToMMo.
Aggregate numbers only: no sample IDs are written."""
import json
import os
import re

import matplotlib
import matplotlib.ticker
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam
import statsmodels.api as sm
from scipy import stats
from sklearn.linear_model import RidgeCV
from sklearn.model_selection import KFold

HERE = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/check/array_crosscheck"
V6 = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6"
MI = f"{V6}/wgs.auto.par/results/12_model_inputs/full_mainland"
TOMMO = ("/LARGE0/gr10478/b37974/Pulmonary_Hypertension/ToMMo_60KJPN/"
         "tommo-60kjpn-20240904-GRCh38-snvindel-af-autosome.norm.vcf.gz")
I, A, C = (f"{HERE}/results/{d}" for d in ("01.inputs", "02.assoc", "03.concordance"))
OUT = f"{HERE}/results/04.compare"
FIG = f"{HERE}/results/figures"
os.makedirs(OUT, exist_ok=True)
os.makedirs(FIG, exist_ok=True)
S = {}                                                   # summary.json


def glm(path):
    d = pd.read_csv(path, sep="\t").rename(columns={"#CHROM": "CHROM"})
    return d[d.TEST == "ADD"].set_index("ID")


def ids(path):
    return set(pd.read_csv(path, sep="\t", header=None)[1].astype(str))


cases, controls = ids(f"{I}/keep.cases.txt"), ids(f"{I}/keep.controls.txt")
S["samples"] = {"cases": len(cases), "controls": len(controls)}

# ── 1. genome-wide ────────────────────────────────────────────────────────────────────────────────────────────────
gw, ga = glm(f"{A}/wgs.PHENO1.glm.logistic.hybrid"), glm(f"{A}/arr.PHENO1.glm.logistic.hybrid")
D = gw.join(ga, lsuffix="_w", rsuffix="_a", how="inner")
flip = D.A1_w != D.A1_a                                  # orient the array estimate to the WGS A1
D["z_a"] = np.where(flip, -D.Z_STAT_a, D.Z_STAT_a)
D["or_a"] = np.where(flip, 1 / D.OR_a, D.OR_a)
D["fca_a"] = np.where(flip, 1 - D.A1_CASE_FREQ_a, D.A1_CASE_FREQ_a)
D["fco_a"] = np.where(flip, 1 - D.A1_CTRL_FREQ_a, D.A1_CTRL_FREQ_a)
for t in ("wgs", "arr"):
    for g in ("cases", "controls"):
        m = pd.read_csv(f"{C}/{t}.{g}.vmiss", sep="\t").set_index("ID")
        D[f"miss_{t}_{g}"] = m.F_MISS.reindex(D.index)
    h = pd.read_csv(f"{C}/{t}.controls.hardy", sep="\t").set_index("ID")
    D[f"hwe_{t}"] = h.P.reindex(D.index)
maf = np.minimum(D.A1_FREQ_w, 1 - D.A1_FREQ_w)
# both sets have passed their own QC; no further missingness filter, because WGS missingness itself differs by
# platform (HiSeq X controls median 3.2%, DNBSeq-T7 cases 1.3%) and a filter on it would select sites by platform
qc = ((maf >= 0.01) & D.Z_STAT_w.notna() & D.Z_STAT_a.notna()
      & (D[["hwe_wgs", "hwe_arr"]].min(axis=1) >= 1e-6))
S["sites"] = {"shared": int(len(D)), "qc_pass": int(qc.sum()), "allele_flipped": int(flip.sum())}

# discordant non-missing calls per site and group (plink 1.9 --merge-mode 7 .diff)
for g, n in (("cases", len(cases)), ("controls", len(controls))):
    p = f"{C}/wgs_vs_arr.{g}.diff"
    dd = pd.read_csv(p, sep=r"\s+", usecols=[0]) if os.path.getsize(p) else pd.DataFrame({"SNP": []})
    D[f"disc_{g}"] = dd.SNP.value_counts().reindex(D.index).fillna(0)
    nonmiss = n * (1 - D[f"miss_wgs_{g}"]) * (1 - D[f"miss_arr_{g}"])      # approximate denominator
    D[f"disc_rate_{g}"] = D[f"disc_{g}"] / nonmiss
    log = open(f"{C}/wgs_vs_arr.{g}.log").read()
    m = re.search(r"(\d+) overlapping calls, (\d+) nonmissing in both filesets\.\s+(\d+) concordant, for a "
                  r"concordance rate of ([0-9]*\.[0-9]+)", log.replace("\n", " "))
    S[f"concordance_all_shared_{g}"] = ({"nonmissing_both": int(m.group(2)), "concordant": int(m.group(3)),
                                         "rate": float(m.group(4))} if m else None)

Q = D[qc].copy()
lam = lambda p: float(np.median(stats.chi2.isf(p, 1)) / stats.chi2.ppf(0.5, 1))   # noqa: E731
dz = Q.z_a - Q.Z_STAT_w
S["genome_wide"] = {
    "lambda_wgs": lam(Q.P_w), "lambda_array": lam(Q.P_a),
    "r_z": float(np.corrcoef(Q.Z_STAT_w, Q.z_a)[0, 1]),
    "slope_z_array_on_wgs": float(np.polyfit(Q.Z_STAT_w, Q.z_a, 1)[0]),
    "abs_dz_gt3": int((dz.abs() > 3).sum()), "abs_dz_gt4": int((dz.abs() > 4).sum()),
    "wgs_p_lt_1e-4": int((Q.P_w < 1e-4).sum()),
    "wgs_p_lt_1e-4_array_same_direction": int(((Q.P_w < 1e-4) & (np.sign(Q.z_a) == np.sign(Q.Z_STAT_w))).sum()),
    "wgs_p_lt_1e-4_array_p_lt_1e-2": int(((Q.P_w < 1e-4) & (Q.P_a < 1e-2)).sum()),
    "median_disc_rate_cases": float(Q.disc_rate_cases.median()),
    "median_disc_rate_controls": float(Q.disc_rate_controls.median()),
    "mean_disc_rate_cases": float(Q.disc_rate_cases.mean()),
    "mean_disc_rate_controls": float(Q.disc_rate_controls.mean()),
    "sites_disc_rate_gt_1pct_cases": int((Q.disc_rate_cases > 0.01).sum()),
    "sites_disc_rate_gt_1pct_controls": int((Q.disc_rate_controls > 0.01).sum()),
    "r_absdz_vs_disc_rate_gap": float(stats.spearmanr(dz.abs(),
                                                      (Q.disc_rate_cases - Q.disc_rate_controls).abs())[0]),
}
keep_cols = ["CHROM_w", "POS_w", "REF_w", "ALT_w", "A1_w", "A1_FREQ_w", "A1_CASE_FREQ_w", "A1_CTRL_FREQ_w",
             "fca_a", "fco_a", "OR_w", "or_a", "Z_STAT_w", "z_a", "P_w", "P_a", "disc_rate_cases",
             "disc_rate_controls", "miss_wgs_cases", "miss_wgs_controls", "miss_arr_cases", "miss_arr_controls"]
Q[keep_cols].to_csv(f"{OUT}/genome_wide.qc_sites.tsv.gz", sep="\t", compression="gzip")
Q.assign(dz=dz)[dz.abs() > 4][keep_cols + ["dz"]].sort_values("dz").to_csv(f"{OUT}/genome_wide.dz_gt4.tsv", sep="\t")

fig, ax = plt.subplots(1, 2, figsize=(10, 4.6))
hb = ax[0].hexbin(Q.Z_STAT_w, Q.z_a, gridsize=120, bins="log", cmap="Blues", mincnt=1)
lim = max(Q.Z_STAT_w.abs().max(), Q.z_a.abs().max()) * 1.05
ax[0].plot([-lim, lim], [-lim, lim], color="#B2182B", lw=0.8, ls="--")
ax[0].set(xlim=(-lim, lim), ylim=(-lim, lim), xlabel="z, WGS genotypes", ylabel="z, array genotypes",
          title=f"Same people, same model: {len(Q):,} sites, r = {S['genome_wide']['r_z']:.3f}")
fig.colorbar(hb, ax=ax[0], label="sites (log)")
for p, lab, col in ((Q.P_w, "WGS", "#0571B0"), (Q.P_a, "array", "#008837")):
    o = -np.log10(np.sort(p.values)); e = -np.log10((np.arange(1, len(o) + 1) - 0.5) / len(o))
    ax[1].plot(e, o, ".", ms=2, color=col, label=f"{lab}: λ = {lam(p):.3f}")
ax[1].plot([0, e.max()], [0, e.max()], color="grey", lw=0.8)
ax[1].set(xlabel="expected −log10 P", ylabel="observed −log10 P", title="QQ, shared QC'd sites")
ax[1].legend(frameon=False)
fig.tight_layout(); fig.savefig(f"{FIG}/genome_wide_wgs_vs_array.png", dpi=200); plt.close(fig)


# ── 2. signals ────────────────────────────────────────────────────────────────────────────────────────────────────
def raw(path):
    r = pd.read_csv(path, sep=r"\s+").drop(columns=["FID", "PAT", "MAT", "SEX", "PHENOTYPE"]).set_index("IID")
    r.index = r.index.astype(str)
    allele = {c.rsplit("_", 1)[0]: c.rsplit("_", 1)[1] for c in r.columns}
    r.columns = [c.rsplit("_", 1)[0] for c in r.columns]
    return r, allele


WR, wal = raw(f"{I}/wgs.window.raw")
AR, aal = raw(f"{I}/arr.window.raw")
gww = glm(f"{A}/wgs.window.PHENO1.glm.logistic.hybrid")
gaw = glm(f"{A}/arr.window.PHENO1.glm.logistic.hybrid")
EU = {e["region"]: e for e in json.load(open(f"{V6}/tmp/ppt/overview/european_loci.json"))}


def eur_full(reg):
    """Full-mainland OR and P of a European lead as the overview reports them, for the ALT allele of its ID."""
    e = EU.get(reg)
    if not e:
        return np.nan, np.nan
    f = e.get("full_mainland") or (e.get("recovered") or {}).get("direct", {}).get("full_mainland")
    if not f:
        return np.nan, np.nan
    return (1 / f["or_"] if e["minor"] == e["ref"] else f["or_"]), f["p"]   # reported for the European minor allele
T = pd.read_csv(f"{I}/targets.tsv", sep="\t", dtype={"chrom": str})
long = pd.read_csv(f"{V6}/integration/saige_signals/results/saige_signals.long.tsv", sep="\t", low_memory=False)
full = long[long.cohort == "full_mainland"].drop_duplicates("variant_id").set_index("variant_id")
ctrl_idx = WR.index.intersection(sorted(controls))
rows = []
for t in T[T.role == "Lead"].itertuples():
    if t.variant_id not in WR:
        a_ = gaw.loc[t.variant_id] if t.variant_id in gaw.index else None
        fo = eur_full(t.signal)[0]                     # for ALT; turned to the array A1 below
        if a_ is not None and a_.A1 != t.variant_id.split(":")[3]:
            fo = 1 / fo
        rows.append(dict(signal=t.signal, kind=t.kind, lead=t.variant_id, rsid=t.rsid, lead_source=t.lead_source,
                         full_or=fo, full_p=eur_full(t.signal)[1],
                         lead_a1=a_.A1 if a_ is not None else None,
                         sub_arr_or_lead=a_.OR if a_ is not None else np.nan,
                         sub_arr_p_lead=a_.P if a_ is not None else np.nan,
                         sub_arr_freq_cases_lead=a_.A1_CASE_FREQ if a_ is not None else np.nan,
                         sub_arr_freq_controls_lead=a_.A1_CTRL_FREQ if a_ is not None else np.nan,
                         note="lead not in the WGS analysis set; typed on the array")); continue
    x = WR.loc[ctrl_idx, t.variant_id]
    win = Q[(Q.CHROM_w.astype(str) == t.chrom) & (Q.POS_w.between(t.pos - 500_000, t.pos + 500_000))]
    cand = [v for v in win.index if v in WR]
    if not cand:
        rows.append(dict(signal=t.signal, lead=t.variant_id, note="no QC'd array site within 500 kb")); continue
    rr = WR.loc[ctrl_idx, cand].corrwith(x)
    best = (rr ** 2).idxmax()
    r = rr[best]
    pr = Q.loc[best]
    # every OR is reported for the lead's A1 (plink2 A1 is not always ALT; the exports count REF): the proxy's A1
    # travels with the lead's A1 when r, taken between counted alleles, agrees after both are turned to their A1
    ld = gww.loc[t.variant_id] if t.variant_id in gww.index else None
    lead_sign = 1 if (ld is not None and ld.A1 == wal[t.variant_id]) else -1
    proxy_sign = 1 if pr.A1_w == wal[best] else -1
    o = lambda v: v if np.sign(r) * lead_sign * proxy_sign > 0 else 1 / v   # noqa: E731
    lead_a1 = ld.A1 if ld is not None else None
    fo = full.OR.get(t.variant_id, eur_full(t.signal)[0])
    fea = full.EA.get(t.variant_id, t.variant_id.split(":")[3])
    fo = fo if (lead_a1 is None or fea == lead_a1) else 1 / fo
    near = win.assign(r2=(WR.loc[ctrl_idx, cand].corrwith(x) ** 2).reindex(win.index))
    near = near[near.r2 >= 0.2]
    rows.append(dict(
        signal=t.signal, kind=t.kind, lead=t.variant_id, rsid=t.rsid, lead_source=t.lead_source,
        full_or=fo, full_p=full.P.get(t.variant_id, eur_full(t.signal)[1]), lead_a1=lead_a1,
        sub_wgs_or_lead=ld.OR if ld is not None else np.nan, sub_wgs_p_lead=ld.P if ld is not None else np.nan,
        proxy=best, proxy_r2=float(r ** 2), proxy_dist_kb=round(abs(int(pr.POS_w) - t.pos) / 1000, 1),
        sub_wgs_or_proxy=o(pr.OR_w), sub_wgs_p_proxy=pr.P_w, sub_arr_or_proxy=o(pr.or_a), sub_arr_p_proxy=pr.P_a,
        disc_cases=pr.disc_rate_cases, disc_controls=pr.disc_rate_controls,
        miss_wgs_cases=pr.miss_wgs_cases, miss_wgs_controls=pr.miss_wgs_controls,
        miss_arr_cases=pr.miss_arr_cases, miss_arr_controls=pr.miss_arr_controls,
        region_sites_r2ge02=int(len(near)),
        region_r_z=float(np.corrcoef(near.Z_STAT_w, near.z_a)[0, 1]) if len(near) >= 3 else np.nan,
        region_max_abs_dz=float((near.z_a - near.Z_STAT_w).abs().max()) if len(near) else np.nan))
SG = pd.DataFrame(rows)
SG.to_csv(f"{OUT}/signals.tsv", sep="\t", index=False)
# credible-set members typed directly on the array
direct = T[T.variant_id.isin(Q.index)].merge(
    Q[["OR_w", "or_a", "P_w", "P_a", "disc_rate_cases", "disc_rate_controls"]], left_on="variant_id", right_index=True)
direct.to_csv(f"{OUT}/signals.direct_on_array.tsv", sep="\t", index=False)
S["signals"] = {"leads": int(len(SG)), "with_proxy_r2_ge_0.8": int((SG.proxy_r2 >= 0.8).sum()),
                "with_proxy_r2_ge_0.5": int((SG.proxy_r2 >= 0.5).sum()),
                "credible_set_members_on_array": int(direct.variant_id.nunique())}

# ── 3. HLA-DQB1 57 V ──────────────────────────────────────────────────────────────────────────────────────────────
res = pd.read_csv(f"{V6}/hla.typing.v2/results/06.assoc_inputs/full_mainland/residue_dosage.tsv", sep="\t",
                  usecols=["sample_id", "DQB1:57:V"], dtype={"sample_id": str}).set_index("sample_id")["DQB1:57:V"]
MX, mal = raw(f"{I}/arr.mhc.raw")
feat = [c for c in MX.columns if c in Q.index]           # array sites that passed the shared-site QC
MX = MX[feat]
y = res.reindex(MX.index)
ok = y.notna()
MX, y = MX[ok], y[ok]
is_case = MX.index.isin(sorted(cases))
Xc, Xk = MX[~is_case], MX[is_case]
mu, sd = Xc.mean(), Xc.std().replace(0, 1)
z = lambda X: ((X.fillna(mu) - mu) / sd).values          # noqa: E731
alphas = np.logspace(-1, 4, 30)
oof = np.zeros((~is_case).sum())
for tr, te in KFold(10, shuffle=True, random_state=1).split(Xc):
    m = RidgeCV(alphas=alphas).fit(z(Xc.iloc[tr]), y[~is_case].iloc[tr])
    oof[te] = m.predict(z(Xc.iloc[te]))
model = RidgeCV(alphas=alphas).fit(z(Xc), y[~is_case])
pred_case = model.predict(z(Xk))
P = pd.Series(np.r_[oof, pred_case], index=list(Xc.index) + list(Xk.index)).clip(0, 2).reindex(MX.index)
yc, yk, pc, pk = y[~is_case], y[is_case], P[~is_case], P[is_case]
resid_c, resid_k = yc - pc, yk - pk
tt = stats.ttest_ind(resid_k, resid_c, equal_var=False)
H = {"features": len(feat), "alpha": float(model.alpha_),
     "n_cases": int(is_case.sum()), "n_controls": int((~is_case).sum()),
     "r2_controls_out_of_fold": float(np.corrcoef(yc, pc)[0, 1] ** 2),
     "r2_cases_trained_on_controls": float(np.corrcoef(yk, pk)[0, 1] ** 2),
     "freq_typed_cases": float(yk.mean() / 2), "freq_typed_controls": float(yc.mean() / 2),
     "freq_array_pred_cases": float(pk.mean() / 2), "freq_array_pred_controls": float(pc.mean() / 2),
     "resid_mean_cases": float(resid_k.mean()), "resid_mean_controls": float(resid_c.mean()),
     "resid_diff_welch_t": float(tt.statistic), "resid_diff_p": float(tt.pvalue),
     "best_call_agreement_cases": float((np.rint(pk) == yk).mean()),
     "best_call_agreement_controls": float((np.rint(pc) == yc).mean())}

cov = pd.read_csv(f"{MI}/covariates/bbj_mainland_pc.sex.tsv", sep="\t", dtype={"IID": str}).set_index("IID")
cv = ["SEX"] + [f"PC{i}_AVG" for i in range(1, 11)]
base = cov.reindex(MX.index)[cv]
base = (base - base.mean()) / base.std()
yy = is_case.astype(int)


def logit(dose):
    X = sm.add_constant(pd.concat([pd.Series(dose.values, index=MX.index, name="d"), base], axis=1))
    f = sm.Logit(yy, X).fit(disp=0)
    b, se = f.params["d"], f.bse["d"]
    return {"or": float(np.exp(b)), "l95": float(np.exp(b - 1.96 * se)), "u95": float(np.exp(b + 1.96 * se)),
            "p": float(f.pvalues["d"])}


H["assoc_typed"] = logit(y)
H["assoc_array_predicted"] = logit(P)
# best single array tags of 57 V, with ToMMo frequencies of the allele that travels with V
r_tag = MX[~is_case].corrwith(yc)
top = (r_tag ** 2).sort_values(ascending=False).head(5)
vcf = pysam.VariantFile(TOMMO)
tags = []
for v, r2 in top.items():
    ch, pos, ref, alt = v.split(":")
    counted = mal[v]
    up = counted if r_tag[v] > 0 else (ref if counted == alt else alt)   # allele positively correlated with V
    g = MX[v] if up == counted else 2 - MX[v]
    taf = None
    for rec in vcf.fetch(ch, int(pos) - 1, int(pos)):
        if rec.pos == int(pos) and rec.ref == ref and alt in rec.alts:
            af = rec.info["AF"][rec.alts.index(alt)]
            taf = af if up == alt else 1 - af
    q = Q.loc[v]
    a_or = q.or_a if q.A1_w == up else 1 / q.or_a
    w_or = q.OR_w if q.A1_w == up else 1 / q.OR_w
    tags.append(dict(tag=v, allele_with_V=up, r2_with_V_controls=float(r2),
                     freq_cases_array=float(g[is_case].mean() / 2), freq_controls_array=float(g[~is_case].mean() / 2),
                     freq_tommo=taf, or_array=float(a_or), p_array=float(q.P_a), or_wgs=float(w_or),
                     p_wgs=float(q.P_w), disc_cases=float(q.disc_rate_cases), disc_controls=float(q.disc_rate_controls)))
pd.DataFrame(tags).to_csv(f"{OUT}/hla57_tags.tsv", sep="\t", index=False)
S["hla_dqb1_57V"] = H
json.dump(S, open(f"{OUT}/summary.json", "w"), indent=2)

fig, ax = plt.subplots(1, 2, figsize=(10, 4.2))
for d, lab, col in ((pd.DataFrame({"t": yc, "p": pc}), "controls (out of fold)", "#0571B0"),
                    (pd.DataFrame({"t": yk, "p": pk}), "cases (model from controls)", "#B2182B")):
    jit = np.random.default_rng(1).normal(0, 0.06, len(d))
    ax[0].plot(d.t + jit + (0.12 if col == "#B2182B" else -0.12), d.p, ".", ms=2.5, alpha=0.5, color=col, label=lab)
ax[0].set(xticks=[0, 1, 2], xlabel="57 V copies, WGS typing", ylabel="57 V dosage predicted from array",
          title="DQB1 57 V: WGS typing vs array prediction")
ax[0].legend(frameon=False, markerscale=4)
for i, (k, lab) in enumerate((("assoc_typed", "WGS typing"), ("assoc_array_predicted", "array prediction"))):
    e = H[k]
    ax[1].errorbar(e["or"], i, xerr=[[e["or"] - e["l95"]], [e["u95"] - e["or"]]], fmt="s", color="#1F3864",
                   capsize=3)
    ax[1].text(e["u95"] * 1.02, i, f"OR {e['or']:.2f} ({e['l95']:.2f}–{e['u95']:.2f}), P = {e['p']:.1e}",
               va="center", fontsize=9)
ax[1].axvline(1, color="grey", ls="--", lw=0.8)
ax[1].xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax[1].set(yticks=[0, 1], yticklabels=["WGS typing", "array prediction"], xscale="log", xlim=(0.3, 3.5),
          xticks=[0.4, 0.6, 1, 2], xticklabels=["0.4", "0.6", "1", "2"],
          ylim=(-0.6, 1.6), xlabel="OR per copy of 57 V (same people, sex + 10 PCs)",
          title=f"{H['n_cases']} cases, {H['n_controls']} controls")
fig.tight_layout(); fig.savefig(f"{FIG}/hla57_typing_vs_array.png", dpi=200); plt.close(fig)
print(json.dumps(S, indent=1))

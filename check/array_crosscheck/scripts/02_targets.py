#!/usr/bin/env python3
"""Step 2: the signals to check and the variants around them.

Leads: the 16 common-variant loci (lead and credible-set members), the MHC lead of the HLA conditioning, and the
European CTEPH loci looked up in our data. Each gets a +/-500 kb window; inside it, the array sites both
technologies carry are the candidate proxies. Writes the ID lists that step 3 exports."""
import json
import os

import pandas as pd

HERE = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/check/array_crosscheck"
V6 = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6"
WGS_BIM = (f"{V6}/wgs.auto.par/results/12_model_inputs/full_mainland/genotype/fixed_model/"
           "cteph_agp3k_v6_wgs_merged.sample_qc.variant_qc.popgmm.fixed_model.fixed_ready.bim")
EURO = f"{V6}/tmp/ppt/overview/european_loci.json"
O = f"{HERE}/results/01.inputs"
FLANK = 500_000

cs = pd.read_csv(f"{V6}/integration/saige_signals/results/00.selection/variants.tsv", sep="\t")
rows = [dict(signal=r.label, kind="common", tier=r.tier, role=r.role, variant_id=r.variant_id, chrom=str(r.chrom),
             pos=int(r.pos), rsid=r.rsID, lead_source="wgs") for r in cs.itertuples()]
hl = pd.read_csv(f"{V6}/analysis/assoc_hla.v2/results/_conditional/gwas_leads.tsv", sep="\t")
for r in hl[hl.status == "kept"].itertuples():
    c, p = r.lead_id.split(":")[:2]
    rows.append(dict(signal="MHC lead", kind="hla", tier="SG", role="Lead", variant_id=r.lead_id, chrom=c[3:],
                     pos=int(p), rsid=".", lead_source="wgs"))

# European loci: the WGS variant at the published position whose alleles match the published ones
eu = [e for e in json.load(open(EURO)) if e.get("full_mainland") or (e.get("recovered") or {}).get("direct")]
want = {(str(e["chrom"]), int(e["pos"])): e for e in eu}
hits = {}
for ch in pd.read_csv(WGS_BIM, sep="\t", header=None, usecols=[0, 1, 3], names=["chr", "id", "pos"],
                      dtype={"chr": str}, chunksize=2_000_000):
    m = ch[[(c, p) in want for c, p in zip(ch.chr, ch.pos)]]
    for r in m.itertuples():
        e = want[(r.chr, r.pos)]
        ref, alt = r.id.split(":")[2:4]
        if {ref, alt} == {e["ref"], e["alt"]}:
            hits[e["region"]] = (r.id, r.chr, r.pos, e["rsid"])
for reg, (vid, c, p, rs) in hits.items():
    rows.append(dict(signal=reg, kind="european", tier="EUR", role="Lead", variant_id=vid, chrom=c, pos=int(p), rsid=rs,
                     lead_source="wgs"))
# a lead that failed the WGS variant QC (ABO, MYH7B: missingness) is taken from the array when the array types it
ARR_BIM = ("/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v5/array/results/10.variant_qc_filtered/"
           "cteph_agp3k.array.sqc.vqc.bim")
ab = pd.read_csv(ARR_BIM, sep="\t", header=None, names=["chr", "id", "cm", "pos", "a1", "a2"], dtype={"chr": str})
for (c, p), e in want.items():
    if e["region"] in hits:
        continue
    m = ab[(ab.chr == c) & (ab.pos == p)]
    m = m[[set(i.split(":")[2:4]) == {e["ref"], e["alt"]} for i in m.id]]
    if len(m):
        rows.append(dict(signal=e["region"], kind="european", tier="EUR", role="Lead", variant_id=m.id.iloc[0], chrom=c,
                         pos=int(p), rsid=e["rsid"], lead_source="array only (failed WGS QC)"))

T = pd.DataFrame(rows).drop_duplicates(["signal", "variant_id"])
T.to_csv(f"{O}/targets.tsv", sep="\t", index=False)
leads = T[T.role == "Lead"]

shared = pd.read_csv(f"{O}/arr.shared.bim", sep="\t", header=None, names=["chr", "id", "cm", "pos", "a1", "a2"],
                     dtype={"chr": str})
inwin = set()
for r in leads.itertuples():
    m = shared[(shared.chr == r.chrom) & (shared.pos.between(r.pos - FLANK, r.pos + FLANK))]
    inwin |= set(m.id)
with open(f"{O}/arr_window.ids", "w") as f:
    f.write("\n".join(sorted(inwin)) + "\n")
with open(f"{O}/wgs_window.ids", "w") as f:
    f.write("\n".join(sorted(inwin | set(T.variant_id))) + "\n")
# the MHC around HLA-DQB1 (chr6:32.66 Mb), for predicting the 57 V dosage from array genotypes alone
with open(f"{O}/lead_windows.range", "w") as f:            # every array site near a lead, shared or not
    for r in leads.itertuples():
        f.write(f"{r.chrom}\t{max(1, r.pos - FLANK)}\t{r.pos + FLANK}\t{r.signal}\n")
with open(f"{O}/mhc.range", "w") as f:
    f.write("6\t32100000\t33300000\tmhc\n")
print(json.dumps({"targets": len(T), "leads": len(leads), "european_found": sorted(hits),
                  "array_sites_in_windows": len(inwin)}))

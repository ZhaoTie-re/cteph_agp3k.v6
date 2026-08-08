#!/usr/bin/env python3
# =============================================================================
# summarize_model_inputs.py
# -----------------------------------------------------------------------------
# Purpose : Per-set consolidated manifest + summary of model-input preparation.
# Project : cteph_agp3k.v6 WGS pipeline (wgs.auto.par/select.auto.par.v6.nf)
# Used by : SUMMARIZE_MODEL_INPUTS
# =============================================================================
"""Per-set consolidated manifest + summary for model-input preparation.

Reads the small artifacts of one sample set's model-input prep (genotype .fam/.bim
line counts, variant lists, PI_HAT exclude list, covariate/phenotype tables) and writes:

  MANIFEST.tsv     one row per prepared artifact: category, file, n_samples,
                   n_variants_or_pcs, note
  prep_summary.log human-readable review block (counts, deltas, provenance)

Pure standard library. Missing/None inputs are tolerated (reported as 0 / n/a) so an
empty rare-variant set does not break the summary.
"""

from __future__ import annotations

import argparse
import glob
import os
from typing import Optional


def _wc(path: Optional[str]) -> int:
    if not path or not os.path.exists(path):
        return 0
    n = 0
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for _ in f:
            n += 1
    return n


def _read_header(path: Optional[str]) -> list[str]:
    if not path or not os.path.exists(path):
        return []
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        first = f.readline().rstrip("\n")
    return first.split("\t") if first else []


def _table_rows(path: Optional[str]) -> int:
    """Data rows (excludes the 1 header line)."""
    n = _wc(path)
    return max(0, n - 1)


def _pheno_counts(path: Optional[str]) -> tuple[int, int, int, int]:
    """Return (n_rows, n_case[PHENO1==2], n_control[==1], n_missing)."""
    if not path or not os.path.exists(path):
        return (0, 0, 0, 0)
    header = _read_header(path)
    try:
        pi = header.index("PHENO1")
    except ValueError:
        return (_table_rows(path), 0, 0, 0)
    n = n_case = n_ctrl = n_missing = 0
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        next(f, None)
        for line in f:
            cols = line.rstrip("\n").split("\t")
            if pi >= len(cols):
                continue
            n += 1
            v = cols[pi].strip()
            if v == "2":
                n_case += 1
            elif v == "1":
                n_ctrl += 1
            else:
                n_missing += 1
    return (n, n_case, n_ctrl, n_missing)


def _n_pcs(header: list[str]) -> int:
    return sum(1 for c in header if c.startswith("PC"))


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Summarize per-set model-input preparation")
    p.add_argument("--set-id", required=True)
    p.add_argument("--base-fam", required=True)
    p.add_argument("--base-bim", required=True)
    p.add_argument("--fixed-ready-fam", required=True)
    p.add_argument("--common-variants", required=True)
    p.add_argument("--rare-variants", default=None)
    p.add_argument("--pihat-exclude", default=None)
    p.add_argument("--random-fam", required=True)
    p.add_argument("--random-bim", required=True)
    p.add_argument("--pheno", required=True)
    p.add_argument("--age-na", default=None)
    p.add_argument("--sex-na", default=None)
    p.add_argument("--cov-dir", default=".", help="Dir to glob *.sex.tsv / *.sex_age.tsv covariate files")
    p.add_argument("--out-manifest", required=True)
    p.add_argument("--out-summary", required=True)
    return p.parse_args()


def main() -> int:
    a = parse_args()

    n_all_samples = _wc(a.base_fam)
    n_all_variants = _wc(a.base_bim)
    n_removed_related = _wc(a.pihat_exclude)
    n_fixed_samples = _wc(a.fixed_ready_fam)
    n_common = _wc(a.common_variants)
    n_rare = _wc(a.rare_variants)
    n_random_samples = _wc(a.random_fam)
    n_random_variants = _wc(a.random_bim)
    n_age_na = _wc(a.age_na)
    n_sex_na = _wc(a.sex_na)
    ph_n, ph_case, ph_ctrl, ph_missing = _pheno_counts(a.pheno)

    cov_sex = sorted(glob.glob(os.path.join(a.cov_dir, "*.sex.tsv")))
    cov_sex_age = sorted(glob.glob(os.path.join(a.cov_dir, "*.sex_age.tsv")))

    # ------------------------------------------------------------------ MANIFEST.tsv
    rows: list[tuple[str, str, str, str, str]] = []
    rows.append(("genotype/base", os.path.basename(a.base_fam).replace(".fam", ".{bed,bim,fam}"),
                 str(n_all_samples), str(n_all_variants), "all set samples, monomorphic-removed"))
    rows.append(("genotype/fixed_model/common", os.path.basename(a.common_variants).replace(".variants.txt", ".{bed,bim,fam}"),
                 str(n_fixed_samples), str(n_common), "MAF>=thr; related-removed; single-variant tests"))
    rows.append(("genotype/fixed_model/rare", (os.path.basename(a.rare_variants).replace(".variants.txt", ".{bed,bim,fam}") if a.rare_variants else "maf_lt_threshold"),
                 str(n_fixed_samples), str(n_rare), "MAF<thr; related-removed; burden/aggregate tests"))
    rows.append(("genotype/random_model", os.path.basename(a.random_fam).replace(".fam", ".{bed,bim,fam}"),
                 str(n_random_samples), str(n_random_variants), "all samples x common variants; GRM/random effect"))
    rows.append(("phenotype", os.path.basename(a.pheno), str(ph_n), "",
                 f"PHENO1 1/2; case={ph_case} control={ph_ctrl} missing={ph_missing}"))
    for c in cov_sex + cov_sex_age:
        hdr = _read_header(c)
        rows.append(("covariates", os.path.basename(c), str(_table_rows(c)), f"{_n_pcs(hdr)} PCs",
                     "+".join([h for h in hdr if h in ("SEX", "AGE", "AGE_Z")])))
    rows.append(("logs", os.path.basename(a.pihat_exclude) if a.pihat_exclude else "pihat_selected.exclude.fid_iid",
                 str(n_removed_related), "", "PI_HAT related samples removed for fixed model"))
    if a.age_na:
        rows.append(("logs", os.path.basename(a.age_na), str(n_age_na), "", "samples dropped from *.sex_age.tsv (missing age)"))
    if a.sex_na:
        rows.append(("logs", os.path.basename(a.sex_na), str(n_sex_na), "", "samples dropped from cov files (unmappable sex)"))

    with open(a.out_manifest, "w", encoding="utf-8") as f:
        f.write("category\tfile\tn_samples\tn_variants_or_pcs\tnote\n")
        for r in rows:
            f.write("\t".join(r) + "\n")

    # ------------------------------------------------------------------ prep_summary.log
    with open(a.out_summary, "w", encoding="utf-8") as f:
        f.write("=" * 78 + "\n")
        f.write(f"MODEL-INPUT PREPARATION SUMMARY — set: {a.set_id}\n")
        f.write("=" * 78 + "\n")
        f.write("Genotype (samples / variants):\n")
        f.write(f"  base (all set samples)          : {n_all_samples} / {n_all_variants}\n")
        f.write(f"  PI_HAT related removed          : {n_removed_related}\n")
        f.write(f"  fixed-model samples (unrelated) : {n_fixed_samples}\n")
        f.write(f"  fixed-model common (MAF>=thr)   : {n_fixed_samples} / {n_common}\n")
        f.write(f"  fixed-model rare   (MAF<thr)    : {n_fixed_samples} / {n_rare}\n")
        f.write(f"  random-model (all x common)     : {n_random_samples} / {n_random_variants}\n")
        f.write("\n")
        f.write("Phenotype (pheno.tsv, PLINK 1/2):\n")
        f.write(f"  N={ph_n}  case(2)={ph_case}  control(1)={ph_ctrl}  missing={ph_missing}\n")
        f.write("\n")
        f.write("Covariates (rows = samples with valid sex; sex_age drops missing age):\n")
        for c in cov_sex + cov_sex_age:
            hdr = _read_header(c)
            f.write(f"  {os.path.basename(c):<34} rows={_table_rows(c):>6}  cols={len(hdr):>3}  PCs={_n_pcs(hdr)}\n")
        f.write(f"  dropped: sex-NA={n_sex_na}  age-NA={n_age_na}\n")
        f.write("\n")
        f.write("Consistency checks:\n")
        f.write(f"  fixed_common_samples == fixed_rare_samples : {n_fixed_samples} (both from fixed_ready)\n")
        f.write(f"  random_variants == fixed_common_variants   : {n_random_variants} vs {n_common} "
                f"({'OK' if n_random_variants == n_common else 'DIFFER — check --extract'})\n")
        f.write(f"  base_samples - related_removed == fixed     : {n_all_samples} - {n_removed_related} "
                f"= {n_all_samples - n_removed_related} vs {n_fixed_samples} "
                f"({'OK' if n_all_samples - n_removed_related == n_fixed_samples else 'DIFFER — check --mac/monomorphic'})\n")
        f.write("=" * 78 + "\n")

    print(f"[OK] wrote {a.out_manifest} and {a.out_summary} for set {a.set_id}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

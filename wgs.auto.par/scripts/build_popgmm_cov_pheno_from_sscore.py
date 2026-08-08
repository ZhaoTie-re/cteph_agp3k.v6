#!/usr/bin/env python3
# =============================================================================
# build_popgmm_cov_pheno_from_sscore.py
# -----------------------------------------------------------------------------
# Purpose : Build covariate/phenotype files from one or more projection sscore files.
# Project : cteph_agp3k.v6 WGS pipeline (wgs.auto.par/select.auto.par.v6.nf)
# Used by : PREPARE_POPGMM_COV_PHENO_FILES
# =============================================================================
"""Build covariate/phenotype files from one or more projection sscore files.

Each ``--sscore LABEL=PATH`` source produces two covariate files:
  <LABEL>.sex.tsv       : #FID IID SEX PC1..PCn                (all mappable samples)
  <LABEL>.sex_age.tsv   : #FID IID SEX AGE AGE_Z PC1..PCn      (age-NA samples dropped)

Exactly one phenotype file is written from the ``--pheno-source`` label:
  pheno.tsv             : #FID IID PHENO1                       (PLINK 1/2 coding)

Sample-level side lists (from the pheno source, no header):
  age_na.fid_iid        : samples dropped from *.sex_age.tsv for missing age
  sex_na.fid_iid        : samples dropped from ALL cov files for unmappable sex

Design notes
------------
* PCs come from the sscore (columns starting with ``PC``, e.g. ``PC1_AVG``);
  ``--n-pcs`` keeps the first N sorted by the integer in the column name.
* age + sex come from the sample-info xlsx, joined on IID.
* ``--keep`` (FID IID, no header) restricts every source to that sample set, so a
  superset projection (e.g. the mainland BBJ projection of full_mainland) is trimmed
  to the current set before covariates are written.
* SEX is mapped male->1, female->2 (PLINK). Samples that cannot be mapped are dropped
  (not a hard error) and recorded in sex_na.fid_iid.
"""

from __future__ import annotations

import argparse
import re
from dataclasses import dataclass, field
from typing import Optional

import pandas as pd


@dataclass
class BuildResult:
    label: str
    sscore_path: str
    n_sscore: int
    n_after_keep: int
    n_sex_na: int
    n_cov_sex: int
    n_age_na: int
    n_cov_sex_age: int
    pc_cols: list[str] = field(default_factory=list)
    sex_counts: dict[str, int] = field(default_factory=dict)


def _normalize_series_as_str(s: pd.Series) -> pd.Series:
    return s.astype(str).str.strip()


def _parse_sscore_arg(values: list[str]) -> list[tuple[str, str]]:
    out: list[tuple[str, str]] = []
    for item in values:
        if "=" not in item:
            raise ValueError(f"--sscore must be LABEL=PATH, got: {item}")
        label, path = item.split("=", 1)
        label = label.strip()
        path = path.strip()
        if not label or not path:
            raise ValueError(f"--sscore must be LABEL=PATH, got: {item}")
        out.append((label, path))
    if not out:
        raise ValueError("At least one --sscore LABEL=PATH is required")
    labels = [lbl for lbl, _ in out]
    if len(set(labels)) != len(labels):
        raise ValueError(f"Duplicate --sscore labels: {labels}")
    return out


def _load_keep_iids(path: str) -> set[str]:
    df = pd.read_csv(path, sep=r"\s+", header=None, engine="python")
    series = df.iloc[:, 1] if df.shape[1] >= 2 else df.iloc[:, 0]
    iids = (
        series.astype(str)
        .str.strip()
        .loc[lambda s: s.ne("") & s.str.lower().ne("nan")]
    )
    return set(iids.tolist())


def _select_pc_cols(df: pd.DataFrame, n_pcs: int) -> list[str]:
    pat = re.compile(r"^PC(\d+)")
    numbered: list[tuple[int, str]] = []
    for c in df.columns:
        m = pat.match(str(c))
        if m:
            numbered.append((int(m.group(1)), c))
    if not numbered:
        raise ValueError("no PC columns found (expected columns starting with 'PC')")
    numbered.sort(key=lambda x: x[0])
    if n_pcs > 0:
        numbered = numbered[:n_pcs]
    return [c for _, c in numbered]


def _read_sscore(sscore_path: str, require_pheno: bool) -> pd.DataFrame:
    df = pd.read_csv(sscore_path, sep=r"\s+", engine="python")
    if "#FID" not in df.columns and "FID" in df.columns:
        df = df.rename(columns={"FID": "#FID"})

    required = ["#FID", "IID"]
    if require_pheno:
        required.append("PHENO1")
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"{sscore_path}: missing required columns: {missing}")

    # validate PCs exist
    _select_pc_cols(df, 0)
    return df


def _build_lookup(sample_info_path: str, sample_id_col: str, sex_col: str, age_col: str) -> pd.DataFrame:
    sample_df = pd.read_excel(sample_info_path)
    required = [sample_id_col, sex_col, age_col]
    missing = [c for c in required if c not in sample_df.columns]
    if missing:
        raise ValueError(f"sample_info missing required columns: {missing}")

    sample_df = sample_df[[sample_id_col, sex_col, age_col]].copy()
    sample_df[sample_id_col] = _normalize_series_as_str(sample_df[sample_id_col])
    sample_df = sample_df.drop_duplicates(subset=[sample_id_col], keep="first")
    return sample_df


def _map_sex(raw: pd.Series, female_value: str, male_value: str) -> pd.Series:
    female_norm = str(female_value).strip().lower()
    male_norm = str(male_value).strip().lower()

    raw_norm = raw.astype(str).str.strip().str.lower()
    mapped = pd.Series(pd.NA, index=raw.index, dtype="Int64")
    mapped.loc[raw_norm == male_norm] = 1
    mapped.loc[raw_norm == female_norm] = 2
    return mapped


def _zscore(s: pd.Series) -> pd.Series:
    mean = s.mean()
    std = s.std(ddof=0)
    if pd.isna(std) or std == 0:
        return pd.Series([0.0] * len(s), index=s.index)
    return (s - mean) / std


def _build_one(
    output_label: str,
    sscore_path: str,
    lookup_df: pd.DataFrame,
    keep_iids: Optional[set[str]],
    n_pcs: int,
    sample_id_col: str,
    sex_col: str,
    age_col: str,
    female_value: str,
    male_value: str,
    is_pheno_source: bool,
) -> tuple[BuildResult, pd.DataFrame]:
    """Write <label>.sex.tsv and <label>.sex_age.tsv; return stats + the merged frame."""
    df = _read_sscore(sscore_path, require_pheno=is_pheno_source).copy()
    df["IID"] = _normalize_series_as_str(df["IID"])
    df["#FID"] = _normalize_series_as_str(df["#FID"])
    n_sscore = int(len(df))

    if keep_iids is not None:
        df = df[df["IID"].isin(keep_iids)].copy()
    n_after_keep = int(len(df))

    pc_cols = _select_pc_cols(df, n_pcs)

    merged = df.merge(
        lookup_df,
        left_on="IID",
        right_on=sample_id_col,
        how="left",
        validate="m:1",
    )

    merged["SEX"] = _map_sex(merged[sex_col], female_value=female_value, male_value=male_value)
    sex_na_df = merged.loc[merged["SEX"].isna(), ["#FID", "IID"]].copy()
    n_sex_na = int(len(sex_na_df))

    mapped = merged.loc[~merged["SEX"].isna()].copy()
    mapped["SEX"] = mapped["SEX"].astype(int)
    mapped["AGE"] = pd.to_numeric(mapped[age_col], errors="coerce")

    # cov.sex.tsv (all sex-mappable samples)
    cov_sex_df = mapped[["#FID", "IID", "SEX"] + pc_cols].copy()
    cov_sex_df.to_csv(f"{output_label}.sex.tsv", sep="\t", index=False)

    # cov.sex_age.tsv (drop age-NA)
    age_ok = mapped.loc[~mapped["AGE"].isna()].copy()
    age_ok["AGE_Z"] = _zscore(age_ok["AGE"].astype(float))
    cov_age_df = age_ok[["#FID", "IID", "SEX", "AGE", "AGE_Z"] + pc_cols].copy()
    cov_age_df.to_csv(f"{output_label}.sex_age.tsv", sep="\t", index=False)

    sex_vc = mapped["SEX"].value_counts().to_dict()
    result = BuildResult(
        label=output_label,
        sscore_path=sscore_path,
        n_sscore=n_sscore,
        n_after_keep=n_after_keep,
        n_sex_na=n_sex_na,
        n_cov_sex=int(len(cov_sex_df)),
        n_age_na=int(mapped["AGE"].isna().sum()),
        n_cov_sex_age=int(len(cov_age_df)),
        pc_cols=pc_cols,
        sex_counts={"male(1)": int(sex_vc.get(1, 0)), "female(2)": int(sex_vc.get(2, 0))},
    )
    # attach side info for the pheno source
    merged.attrs["sex_na_df"] = sex_na_df
    merged.attrs["mapped"] = mapped
    return result, merged


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build cov/pheno files from one or more sscore files")
    parser.add_argument(
        "--sscore",
        action="append",
        required=True,
        metavar="LABEL=PATH",
        help="Projection sscore as LABEL=PATH (repeatable). Each yields <LABEL>.sex.tsv and <LABEL>.sex_age.tsv",
    )
    parser.add_argument("--pheno-source", required=True, help="Which --sscore LABEL writes pheno.tsv / age_na / sex_na")
    parser.add_argument("--keep", default=None, help="Optional FID IID keep list (no header); restricts every source")
    parser.add_argument("--n-pcs", type=int, default=20, help="Number of PCs to keep (first N by PC number)")
    parser.add_argument("--sample-info", required=True)
    parser.add_argument("--sample-id-col", required=True)
    parser.add_argument("--sex-col", required=True)
    parser.add_argument("--sex-female-value", required=True)
    parser.add_argument("--sex-male-value", required=True)
    parser.add_argument("--age-col", required=True)
    parser.add_argument("--out-log", required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    sources = _parse_sscore_arg(args.sscore)
    labels = [lbl for lbl, _ in sources]
    if args.pheno_source not in labels:
        raise ValueError(f"--pheno-source '{args.pheno_source}' not among --sscore labels {labels}")

    keep_iids = _load_keep_iids(args.keep) if args.keep else None
    keep_n = len(keep_iids) if keep_iids is not None else None

    lookup_df = _build_lookup(
        sample_info_path=args.sample_info,
        sample_id_col=args.sample_id_col,
        sex_col=args.sex_col,
        age_col=args.age_col,
    )

    results: list[BuildResult] = []
    pheno_merged: Optional[pd.DataFrame] = None
    for label, path in sources:
        res, merged = _build_one(
            output_label=label,
            sscore_path=path,
            lookup_df=lookup_df,
            keep_iids=keep_iids,
            n_pcs=args.n_pcs,
            sample_id_col=args.sample_id_col,
            sex_col=args.sex_col,
            age_col=args.age_col,
            female_value=args.sex_female_value,
            male_value=args.sex_male_value,
            is_pheno_source=(label == args.pheno_source),
        )
        results.append(res)
        if label == args.pheno_source:
            pheno_merged = merged

    assert pheno_merged is not None
    mapped = pheno_merged.attrs["mapped"]
    sex_na_df = pheno_merged.attrs["sex_na_df"]

    # pheno.tsv (PLINK 1/2 coding straight from the sscore PHENO1)
    pheno_df = mapped[["#FID", "IID", "PHENO1"]].copy()
    pheno_df.to_csv("pheno.tsv", sep="\t", index=False)

    # side lists (no header)
    age_na_df = mapped.loc[mapped["AGE"].isna(), ["#FID", "IID"]].copy()
    age_na_df.to_csv("age_na.fid_iid", sep="\t", index=False, header=False)
    sex_na_df.to_csv("sex_na.fid_iid", sep="\t", index=False, header=False)

    pheno_num = pd.to_numeric(pheno_df["PHENO1"], errors="coerce")
    n_case = int((pheno_num == 2).sum())
    n_ctrl = int((pheno_num == 1).sum())
    n_pheno_missing = int(pheno_num.isna().sum() + ((pheno_num != 1) & (pheno_num != 2) & pheno_num.notna()).sum())

    with open(args.out_log, "w", encoding="utf-8") as f:
        f.write("=" * 78 + "\n")
        f.write("Covariate / phenotype build summary\n")
        f.write("=" * 78 + "\n")
        f.write(f"KEEP_FILE: {args.keep if args.keep else '(none)'}\n")
        if keep_n is not None:
            f.write(f"KEEP_N_IIDS: {keep_n}\n")
        f.write(f"N_PCS_REQUESTED: {args.n_pcs}\n")
        f.write(f"PHENO_SOURCE: {args.pheno_source}\n")
        f.write(f"SAMPLE_INFO: {args.sample_info}\n")
        f.write(f"SEX_COL: {args.sex_col} (male='{args.sex_male_value}'->1, female='{args.sex_female_value}'->2)\n")
        f.write(f"AGE_COL: {args.age_col}\n")
        f.write("\n")
        f.write("Phenotype file: pheno.tsv (PLINK PHENO1 1=control / 2=case)\n")
        f.write(f"  N_CASE(2)   : {n_case}\n")
        f.write(f"  N_CONTROL(1): {n_ctrl}\n")
        f.write(f"  N_MISSING   : {n_pheno_missing}\n")
        f.write(f"  AGE_NA_FILE : age_na.fid_iid (N={int(len(age_na_df))})\n")
        f.write(f"  SEX_NA_FILE : sex_na.fid_iid (N={int(len(sex_na_df))})\n")
        f.write("\n")
        for res in results:
            f.write("-" * 78 + "\n")
            f.write(f"SOURCE_LABEL: {res.label}\n")
            f.write(f"  SSCORE            : {res.sscore_path}\n")
            f.write(f"  N_IN_SSCORE       : {res.n_sscore}\n")
            f.write(f"  N_AFTER_KEEP      : {res.n_after_keep}\n")
            f.write(f"  N_SEX_NA_DROPPED  : {res.n_sex_na}\n")
            f.write(f"  N_COV_SEX_ROWS    : {res.n_cov_sex}   -> {res.label}.sex.tsv\n")
            f.write(f"  N_AGE_NA_DROPPED  : {res.n_age_na}\n")
            f.write(f"  N_COV_SEX_AGE_ROWS: {res.n_cov_sex_age}   -> {res.label}.sex_age.tsv\n")
            f.write(f"  N_PCS_USED        : {len(res.pc_cols)} ({res.pc_cols[0]}..{res.pc_cols[-1]})\n")
            f.write(f"  SEX_COUNTS        : {res.sex_counts}\n")
        f.write("=" * 78 + "\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

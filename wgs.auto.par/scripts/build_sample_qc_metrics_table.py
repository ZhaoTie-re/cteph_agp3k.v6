#!/usr/bin/env python3
# =============================================================================
# build_sample_qc_metrics_table.py
# -----------------------------------------------------------------------------
# Purpose : Assemble the sample-level QC metrics table (Het_F, DP robust-Z, SMISS, platform).
# Project : cteph_agp3k.v6 WGS pipeline (wgs.auto.par/select.auto.par.v6.nf)
# Used by : BUILD_SAMPLE_QC_TABLE
# =============================================================================

import argparse
import csv
import math
import sys
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
	parser = argparse.ArgumentParser(
		description=(
			"Build sample-level QC metrics table from merged PLINK outputs, sample info, and plink2 missingness."
		)
	)
	parser.add_argument("--xlsx", required=True, help="Input sample info Excel workbook path")
	parser.add_argument("--fam", required=True, help="Merged PLINK .fam path")
	parser.add_argument("--smiss", required=True, help="plink2 .smiss path")
	parser.add_argument("--het", required=True, help="plink2 .het path")
	parser.add_argument("--out", required=True, help="Output TSV path")
	parser.add_argument("--sample-id-col", required=True, help="Column name containing sample IDs")
	parser.add_argument("--target-dp-col", required=True, help="Column name containing Target DP")
	parser.add_argument("--dp-col", required=True, help="Column name containing DP")
	parser.add_argument(
		"--platform-col",
		default="WGS_Platform",
		help="Column name containing the sequencing platform (used to stratify the DP robust-Z)",
	)
	return parser.parse_args()


def normalize_text(value) -> str | None:
	if value is None or pd.isna(value):
		return None
	text = str(value).strip()
	if not text or text.lower() == "nan":
		return None
	return text


def load_fam(fam_path: Path) -> pd.DataFrame:
	rows: list[tuple[str, str]] = []
	with fam_path.open("r", encoding="utf-8") as handle:
		for line in handle:
			parts = line.strip().split()
			if len(parts) < 2:
				continue
			rows.append((parts[0], parts[1]))

	fam_df = pd.DataFrame(rows, columns=["FID", "IID"]) if rows else pd.DataFrame(columns=["FID", "IID"])
	if fam_df.empty:
		raise ValueError(f"No valid FID/IID records found in fam file: {fam_path}")

	if fam_df["IID"].duplicated().any():
		dups = fam_df.loc[fam_df["IID"].duplicated(), "IID"].head(10).tolist()
		raise ValueError(f"Duplicated IID in fam file: {', '.join(dups)}")

	return fam_df


def load_sample_info(
	xlsx_path: Path,
	sample_id_col: str,
	target_dp_col: str,
	dp_col: str,
	platform_col: str,
) -> pd.DataFrame:
	df = pd.read_excel(xlsx_path, engine="openpyxl")

	required = [sample_id_col, target_dp_col, dp_col, platform_col]
	missing = [col for col in required if col not in df.columns]
	if missing:
		raise ValueError(f"Missing required column(s) in sample info: {', '.join(missing)}")

	sub = df[[sample_id_col, target_dp_col, dp_col, platform_col]].copy()
	sub.columns = ["IID", "Target_DP", "DP", "WGS_Platform"]
	sub["IID"] = sub["IID"].map(normalize_text)
	sub = sub[sub["IID"].notna()].copy()

	if sub.empty:
		raise ValueError("No usable sample records after filtering empty sample IDs")

	sub["Target_DP"] = sub["Target_DP"].map(normalize_text)
	sub["WGS_Platform"] = sub["WGS_Platform"].map(normalize_text)
	sub["DP"] = pd.to_numeric(sub["DP"], errors="coerce")

	duplicated = sub["IID"].duplicated(keep=False)
	if duplicated.any():
		dup_rows = sub.loc[duplicated, ["IID", "Target_DP", "DP", "WGS_Platform"]].copy()
		agg = dup_rows.groupby("IID", dropna=False).nunique(dropna=True)
		conflict_ids = agg[
			(agg["Target_DP"] > 1) | (agg["DP"] > 1) | (agg["WGS_Platform"] > 1)
		].index.tolist()
		if conflict_ids:
			preview = ", ".join(sorted(conflict_ids)[:10])
			raise ValueError(f"Conflicting duplicate IID found in sample info: {preview}")
		sub = sub.drop_duplicates(subset=["IID"], keep="first")

	return sub


def load_smiss(smiss_path: Path) -> pd.DataFrame:
	smiss_df = pd.read_csv(smiss_path, sep=r"\s+", engine="python")

	fid_col = "#FID" if "#FID" in smiss_df.columns else ("FID" if "FID" in smiss_df.columns else None)
	if fid_col is None:
		raise ValueError(f"Cannot find FID column in smiss file: {smiss_path}")
	if "IID" not in smiss_df.columns:
		raise ValueError(f"Cannot find IID column in smiss file: {smiss_path}")

	if "F_MISS" in smiss_df.columns:
		smiss_col = "F_MISS"
	elif "N_MISS" in smiss_df.columns and "N_GENO" in smiss_df.columns:
		smiss_col = None
	else:
		raise ValueError("smiss file must contain F_MISS or both N_MISS and N_GENO columns")

	if smiss_col is None:
		smiss_df["SMISS"] = smiss_df["N_MISS"] / smiss_df["N_GENO"].replace(0, pd.NA)
	else:
		smiss_df["SMISS"] = pd.to_numeric(smiss_df[smiss_col], errors="coerce")

	out = smiss_df[[fid_col, "IID", "SMISS"]].copy()
	out.columns = ["FID_smiss", "IID", "SMISS"]
	out["IID"] = out["IID"].map(normalize_text)
	out = out[out["IID"].notna()].copy()

	if out["IID"].duplicated().any():
		dups = out.loc[out["IID"].duplicated(), "IID"].head(10).tolist()
		raise ValueError(f"Duplicated IID in smiss file: {', '.join(dups)}")

	return out


def load_het(het_path: Path) -> pd.DataFrame:
	het_df = pd.read_csv(het_path, sep=r"\s+", engine="python")

	fid_col = "#FID" if "#FID" in het_df.columns else ("FID" if "FID" in het_df.columns else None)
	if fid_col is None:
		raise ValueError(f"Cannot find FID column in het file: {het_path}")
	if "IID" not in het_df.columns:
		raise ValueError(f"Cannot find IID column in het file: {het_path}")
	if "F" not in het_df.columns:
		raise ValueError(f"Cannot find F column in het file: {het_path}")

	out = het_df[[fid_col, "IID", "F"]].copy()
	out.columns = ["FID_het", "IID", "Het_F"]
	out["IID"] = out["IID"].map(normalize_text)
	out["Het_F"] = pd.to_numeric(out["Het_F"], errors="coerce")
	out = out[out["IID"].notna()].copy()

	if out["IID"].duplicated().any():
		dups = out.loc[out["IID"].duplicated(), "IID"].head(10).tolist()
		raise ValueError(f"Duplicated IID in het file: {', '.join(dups)}")

	return out


def calc_group_robust_z(df: pd.DataFrame, group_col: str = "WGS_Platform") -> pd.Series:
	# DP outlier detection asks "is this sample anomalously low FOR ITS OWN PLATFORM"
	# (a failed library), so the robust-Z is computed within WGS_Platform. Grouping by
	# Target_DP would pool platforms of genuinely different depth (T7 ~18x, NovaSeq ~30x,
	# G400RS ~35x all sit in "30x") and flag platform membership, not sample quality.
	def _per_group(values: pd.Series) -> pd.Series:
		numeric = pd.to_numeric(values, errors="coerce")
		if numeric.notna().sum() == 0:
			return pd.Series([math.nan] * len(values), index=values.index)
		median = numeric.median(skipna=True)
		mad = (numeric - median).abs().median(skipna=True)
		if pd.isna(mad) or mad == 0:
			return pd.Series([math.nan] * len(values), index=values.index)
		return (numeric - median) / (1.4826 * mad)

	return df.groupby(group_col, dropna=False)["DP"].transform(_per_group)


def write_output(df: pd.DataFrame, out_path: Path) -> None:
	# Columns 1-6 (#FID..SMISS) are byte-identical to the previous version so downstream
	# by-name readers (run_pihat_network_qc.py uses IID, SMISS) are unaffected. Two columns
	# are appended: DP_RobustZ_in_Platform (was DP_RobustZ_in_TargetDP, now grouped by
	# platform) and WGS_Platform (needed by run_sample_qc.py for per-platform Het_F).
	with out_path.open("w", encoding="utf-8", newline="") as handle:
		writer = csv.writer(handle, delimiter="\t")
		writer.writerow(
			["#FID", "IID", "Het_F", "Target_DP", "DP", "SMISS", "DP_RobustZ_in_Platform", "WGS_Platform"]
		)
		for row in df.itertuples(index=False):
			fid = "" if pd.isna(row.FID) else str(row.FID)
			iid = "" if pd.isna(row.IID) else str(row.IID)
			het_f = "" if pd.isna(row.Het_F) else f"{float(row.Het_F):.6f}"
			target_dp = "" if pd.isna(row.Target_DP) else str(row.Target_DP)
			dp = "" if pd.isna(row.DP) else f"{float(row.DP):.6f}"
			smiss = "" if pd.isna(row.SMISS) else f"{float(row.SMISS):.8f}"
			rz = "" if pd.isna(row.DP_RobustZ_in_Platform) else f"{float(row.DP_RobustZ_in_Platform):.6f}"
			platform = "" if pd.isna(row.WGS_Platform) else str(row.WGS_Platform)
			writer.writerow([fid, iid, het_f, target_dp, dp, smiss, rz, platform])


def main() -> int:
	args = parse_args()
	try:
		fam_df = load_fam(Path(args.fam))
		info_df = load_sample_info(
			Path(args.xlsx),
			args.sample_id_col,
			args.target_dp_col,
			args.dp_col,
			args.platform_col,
		)
		smiss_df = load_smiss(Path(args.smiss))
		het_df = load_het(Path(args.het))

		merged = fam_df.merge(info_df, on="IID", how="left")
		merged = merged.merge(het_df[["IID", "Het_F"]], on="IID", how="left")
		merged = merged.merge(smiss_df[["IID", "SMISS"]], on="IID", how="left")
		merged["DP_RobustZ_in_Platform"] = calc_group_robust_z(merged, group_col="WGS_Platform")

		out_path = Path(args.out)
		write_output(merged, out_path)

		missing_info = int(merged["Target_DP"].isna().sum())
		missing_het = int(merged["Het_F"].isna().sum())
		missing_smiss = int(merged["SMISS"].isna().sum())
		missing_platform = int(merged["WGS_Platform"].isna().sum())
		print(
			(
				f"Wrote {len(merged)} samples to {out_path}; "
				f"missing sample_info={missing_info}, missing het={missing_het}, "
				f"missing smiss={missing_smiss}, missing platform={missing_platform}"
			),
			file=sys.stderr,
		)
		return 0
	except Exception as exc:  # noqa: BLE001
		print(f"ERROR: {exc}", file=sys.stderr)
		return 1


if __name__ == "__main__":
	sys.exit(main())
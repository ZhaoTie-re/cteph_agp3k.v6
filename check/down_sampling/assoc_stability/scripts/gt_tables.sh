#!/bin/bash
# ---------------------------------------------------------------------------
# Dump the genotypes the two associations actually read, for ONE arm, so the
# report can show what the down-sampling did to them rather than only what it did
# to the p-value.
#
#   $1  ARM       baseline | rep1 | rep2 | rep3
#   $2  SNP_BED   PLINK prefix the SNP association reads
#   $3  VARIANT   the lead variant
#   $4  GENE_VCF  the VCF the gene-based association reads
#   $5  OUT_SNP   output: SAMPLE <TAB> VARIANT <TAB> DOSAGE ('NA' when missing)
#   $6  OUT_GENE  output: same, for every variant the gene test includes
#
# BOTH SIDES ARE POST-GENOTYPE-QC. The associations run on the post_qc VCFs, so
# these are the same genotypes — including their no-calls. A no-call here is not an
# absence of data to be tidied away: genotype QC turns a doubtful call INTO one,
# which is the mechanism by which depth reaches an odds ratio. It has to survive
# into the confusion matrix as its own category.
#
# Assumes PATH/conda are set by the caller (plink2, bcftools).
# ---------------------------------------------------------------------------
set -euo pipefail

ARM=$1; SNP_BED=$2; VARIANT=$3; GENE_VCF=$4; OUT_SNP=$5; OUT_GENE=$6

# ── the SNP lead ───────────────────────────────────────────────────────────
# The ALT is named from the ID (CHROM:POS:REF:ALT) so the dosage counts the ALT
# and not whichever allele PLINK would otherwise choose as A1.
echo "${VARIANT}" > lead.txt
awk -F: '{print $0"\t"$NF}' lead.txt > lead_alt.txt
plink2 --bfile "${SNP_BED}" --extract lead.txt \
    --export A --export-allele lead_alt.txt --out lead_dos

python3 - "${ARM}" lead_dos.raw "${OUT_SNP}" <<'PY'
import sys
arm, raw, out = sys.argv[1], sys.argv[2], sys.argv[3]
with open(raw) as fh, open(out, "w") as o:
    hdr = fh.readline().rstrip("\n").split("\t")
    vids = [c.rsplit("_", 1)[0] for c in hdr[6:]]
    for line in fh:
        f = line.rstrip("\n").split("\t")
        for v, x in zip(vids, f[6:]):
            o.write(f"{f[1]}\t{v}\t{x}\n")
PY
echo "[gt_tables] ${ARM}: $(wc -l < "${OUT_SNP}") SNP genotypes" >&2

# ── every variant the gene test includes ──────────────────────────────────
# bcftools writes '.' for a no-call; normalise it to 'NA' so both tables use one
# spelling and the report cannot silently treat one as a called genotype.
bcftools query -f '[%SAMPLE\t%ID\t%GT\n]' "${GENE_VCF}" \
    | awk -F'\t' 'BEGIN{OFS="\t"}
        {
          gt = $3
          gsub(/\|/, "/", gt)
          if (gt ~ /\./) { d = "NA" }
          else { split(gt, a, "/"); d = 0; for (i in a) if (a[i] != "0") d++ }
          print $1, $2, d
        }' > "${OUT_GENE}"
echo "[gt_tables] ${ARM}: $(wc -l < "${OUT_GENE}") gene genotypes over" \
     "$(cut -f2 "${OUT_GENE}" | sort -u | wc -l) variants" >&2

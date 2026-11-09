#!/usr/bin/env bash
# HaplotypeCaller diagnostics at the ADARB2 lead (user 2026-11-06: "我们先看看HaplotypeCaller 的诊断文件").
# Read only on the CRAMs. Per genome of samples.tsv: the reads of the extraction window (as EXTRACT_READS of
# platform_depth.nf), then the cohort's HaplotypeCaller command exactly as CALL_GVCF runs it (GATK 4.1.4.0, Java 8,
# -ERC GVCF, defaults, the calling interval chr10:1735001-1760000 of results/02.targets/hc.intervals), with output-only
# diagnostics added: the assembly debug log, every assembled haplotype (bamout, ALL_POSSIBLE_HAPLOTYPES), the assembly
# regions and the activity profile (IGV), and the assembly graphs. None of these changes the calls.
# Usage: bash run_hc_debug.sh [n_parallel]     (run inside a SLURM job; outputs in out/<sample_id>/)
set -euo pipefail
cd "$(dirname "$0")"
NP=${1:-4}
SAMTOOLS=/LARGE0/gr10478/b37974/software/samtools/bin/samtools
JAVA8=/LARGE0/gr10478/workspace/pipeline/packages/jdk8u262-b10-jre/bin/java
GATK=/LARGE0/gr10478/workspace/pipeline/packages/gatk-4.1.4.0/gatk-package-4.1.4.0-local.jar
FASTA=/LARGE0/gr10478/b37974/Pulmonary_Hypertension/nagasaki_pipeline/data/hs38DH.fa
WINDOW="chr10:1734001-1761000"            # extract.bed row chr10 1734000 1761000 (0-based)
INTERVAL="chr10:1735001-1760000"          # the calling interval holding the lead
one() {
  local sid=$1 cram=$2 crai=$3
  local o=out/$sid
  mkdir -p "$o"
  "$SAMTOOLS" view -b --reference "$FASTA" -o "$o/$sid.region.bam" "${cram}##idx##${crai}" "$WINDOW"
  "$SAMTOOLS" index "$o/$sid.region.bam"
  "$JAVA8" -Xmx6g -Xms6g -jar "$GATK" HaplotypeCaller \
      -R "$FASTA" -I "$o/$sid.region.bam" -O "$o/$sid.g.vcf.gz" -ERC GVCF -L "$INTERVAL" \
      --native-pair-hmm-threads 2 \
      --debug-assembly --bam-output "$o/$sid.bamout.bam" --bam-writer-type ALL_POSSIBLE_HAPLOTYPES \
      --assembly-region-out "$o/$sid.regions.igv" --activity-profile-out "$o/$sid.activity.igv" \
      --graph-output "$o/$sid.graph.dot" > "$o/$sid.hc.log" 2>&1
  echo "done $sid"
}
export -f one
export SAMTOOLS JAVA8 GATK FASTA WINDOW INTERVAL
tail -n +2 samples.tsv | awk -F'\t' '{print $1"\t"$4"\t"$5}' |
  xargs -P "$NP" -d '\n' -I{} bash -c 'IFS=$'"'"'\t'"'"' read -r s c i <<< "{}"; one "$s" "$c" "$i"'
echo "all done"

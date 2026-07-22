#!/usr/bin/env python3
"""Variant-level recovery of a replicate VCF against one or more reference sets.

For a down-sampling replicate we report, against each reference variant set:

  * shared   -> reference variant recovered in the replicate
  * missed   -> reference variant NOT recovered (depth loss)
  * extra    -> replicate variant absent from the reference

Two references are compared (either may be omitted):
  1. force_call_sites  -- the cohort PASS sites we force-called (a VCF)
  2. refined_core      -- the analysis-ready variant set (a list of CHROM:POS:REF:ALT
                          IDs from the refined_core PLINK bim, region-restricted)

All inputs must share the same normalization (split multiallelics, left-aligned,
ID = CHROM:POS:REF:ALT); variants are matched by CHROM:POS:REF:ALT.

bcftools (on PATH) is used to stream VCF variant keys; no VCF library needed.
"""

import argparse
import subprocess
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path


def vcf_keys(vcf):
    """Set of 'CHROM:POS:REF:ALT' keys in a (bgzipped, indexed) VCF."""
    proc = subprocess.run(
        ["bcftools", "query", "-f", "%CHROM:%POS:%REF:%ALT\n", vcf],
        check=True, capture_output=True, text=True,
    )
    return {ln.strip() for ln in proc.stdout.splitlines() if ln.strip()}


def file_keys(path):
    """Set of keys from a plain text file (one CHROM:POS:REF:ALT per line)."""
    with open(path) as fh:
        return {ln.strip() for ln in fh if ln.strip()}


def parse_key(key):
    chrom, pos, ref, alt = key.split(":", 3)
    return chrom, pos, ref, alt


def vtype(ref, alt):
    if alt in (".", "*"):
        return "other"
    if len(ref) == 1 and len(alt) == 1:
        return "SNP"
    if len(ref) < len(alt) and alt.startswith(ref[:1]):
        return "INS"
    if len(ref) > len(alt) and ref.startswith(alt[:1]):
        return "DEL"
    return "MNP/complex"


def tally(keys):
    by_type = defaultdict(int)
    by_chrom = defaultdict(int)
    for k in keys:
        chrom, _pos, ref, alt = parse_key(k)
        by_type[vtype(ref, alt)] += 1
        by_chrom[chrom] += 1
    return by_type, by_chrom


def pct(n, d):
    return (100.0 * n / d) if d else 0.0


TYPE_ORDER = ["SNP", "INS", "DEL", "MNP/complex", "other"]


def compare_block(rep, ref, ref_name, note, out_dir, label, missed_file, extra_file):
    """Return the log lines for rep-vs-ref, and dump the diff variant lists."""
    shared = rep & ref
    missed = ref - rep       # reference variant NOT recovered
    extra = rep - ref        # replicate variant absent from reference

    (out_dir / missed_file).write_text("\n".join(sorted(missed)) + ("\n" if missed else ""))
    (out_dir / extra_file).write_text("\n".join(sorted(extra)) + ("\n" if extra else ""))

    sh_t, _ = tally(shared)
    mi_t, mi_c = tally(missed)
    ex_t, ex_c = tally(extra)
    ref_t, ref_c = tally(ref)
    rep_t, rep_c = tally(rep)

    L = []
    def line(s=""):
        L.append(s)

    line("-" * 80)
    line(f" REFERENCE: {ref_name}")
    line("-" * 80)
    line(f"   {note}")
    line("")
    line(f"   reference variants : {len(ref):>8}")
    line(f"   replicate variants : {len(rep):>8}")
    line(f"   shared (recovered) : {len(shared):>8}   ({pct(len(shared), len(ref)):6.2f}% of reference)")
    line(f"   missed  (depth loss): {len(missed):>8}   ({pct(len(missed), len(ref)):6.2f}% of reference)")
    line(f"   extra  (not in ref) : {len(extra):>8}")
    line("")
    line(f"   {'type':<12}{'ref':>8}{'rep':>8}{'shared':>8}{'missed':>8}{'extra':>8}")
    for t in TYPE_ORDER:
        if ref_t.get(t) or rep_t.get(t) or mi_t.get(t) or ex_t.get(t):
            line(f"   {t:<12}{ref_t.get(t,0):>8}{rep_t.get(t,0):>8}"
                 f"{sh_t.get(t,0):>8}{mi_t.get(t,0):>8}{ex_t.get(t,0):>8}")
    line("")
    line(f"   {'chrom':<12}{'ref':>8}{'rep':>8}{'missed':>8}{'extra':>8}")
    for c in sorted(set(ref_c) | set(rep_c)):
        line(f"   {c:<12}{ref_c.get(c,0):>8}{rep_c.get(c,0):>8}{mi_c.get(c,0):>8}{ex_c.get(c,0):>8}")
    line("")
    line(f"   missed -> {missed_file}      extra -> {extra_file}")
    line("")
    return L


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rep-vcf", required=True, help="normalized replicate VCF")
    ap.add_argument("--sites-vcf", help="normalized force-call sites VCF")
    ap.add_argument("--refined-ids", help="refined_core variant IDs (region-restricted), one per line")
    ap.add_argument("--label", required=True)
    ap.add_argument("--out-log", required=True)
    ap.add_argument("--out-dir", default=".")
    args = ap.parse_args()

    rep = vcf_keys(args.rep_vcf)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    W = 80
    head = [
        "=" * W,
        f" Variant recovery of down-sampled replicate   —   {args.label}",
        "=" * W,
        f" generated : {datetime.now().isoformat(timespec='seconds')}",
        f" rep_vcf   : {args.rep_vcf}   ({len(rep)} variants)",
        "",
        " shared = reference variant recovered; missed = NOT recovered (depth loss);",
        " extra  = replicate variant absent from the reference. Matched by CHROM:POS:REF:ALT.",
        "",
    ]
    blocks = []
    if args.sites_vcf:
        blocks += compare_block(
            rep, vcf_keys(args.sites_vcf), "force_call_sites",
            "cohort PASS sites we force-called (recovery of the forced target set)",
            out_dir, args.label, f"{args.label}.sites_only.txt", f"{args.label}.rep_only.txt")
    if args.refined_ids:
        blocks += compare_block(
            rep, file_keys(args.refined_ids), "refined_core",
            "analysis-ready variants (QC'd, pruned) — recovery of the study variant set",
            out_dir, args.label, f"{args.label}.refined_missed.txt", f"{args.label}.refined_extra.txt")

    text = "\n".join(head + blocks + ["=" * W, ""])
    Path(args.out_log).write_text(text)
    sys.stdout.write(text)


if __name__ == "__main__":
    main()

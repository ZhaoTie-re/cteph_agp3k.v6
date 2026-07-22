#!/usr/bin/env python3
"""Copy FORMAT/RGQ into FORMAT/GQ on records that carry only RGQ. VCF on stdin,
VCF on stdout. No genotype is altered — this only fills in a field GATK wrote
under a different name.

WHY THIS EXISTS

  --force-output-intervals makes GenotypeGVCFs emit a record at every cohort site,
  including sites where not one of our samples carries the ALT. There is no variant
  to genotype there, so GATK writes the record the way it writes a GVCF reference
  block: FORMAT is 'GT:AD:DP:RGQ', and the confidence lives in RGQ (reference
  genotype quality). GQ is not merely empty — the field is absent from FORMAT.

  The cohort's genotype QC then filters on `GQ < 20`. Against a missing GQ, GATK's
  JEXL evaluates the expression as if GQ were 0, so 0 < 20 is true and the genotype
  is filtered — a confident hom-ref with RGQ 36 and DP 12 is thrown away for having
  no GQ, not for being poor. In this cohort that silently no-called ~41% of all
  genotypes, essentially all of them hom-ref.

  That is not a depth effect. The cohort's own calls never hit it, because without
  force-calling a site only exists when some sample carries the ALT, and then every
  sample there gets a real GQ. It is an artefact of forcing sites, and it lands on
  the down-sampled samples only — which, since platform tracks case status here,
  means it lands on the cases only. Left alone it would hand the association a 41%
  vs 1.2% differential missingness that has nothing to do with biology.

  Filling GQ from RGQ lets the cohort's genotype-QC expressions run unchanged on
  our replicates. The alternative — teaching the filter about RGQ — would give our
  pipeline a QC rule the cohort never had, and the whole point is to compare like
  with like.

NOTE

  'GQ' in FORMAT is matched as a whole field, never as a substring: 'RGQ' contains
  'GQ', so a substring test would decide every RGQ record already has a GQ and this
  script would do nothing at all.
"""

import sys


def main():
    n_rec = n_fixed = 0
    out = sys.stdout
    for line in sys.stdin:
        if line.startswith("#"):
            out.write(line)
            continue
        n_rec += 1
        f = line.rstrip("\n").split("\t")
        if len(f) < 10:
            out.write(line)
            continue
        fmt = f[8].split(":")
        # Whole-field match: 'RGQ'.__contains__('GQ') is True, hence not `in f[8]`.
        if "RGQ" not in fmt or "GQ" in fmt:
            out.write(line)
            continue
        ridx = fmt.index("RGQ")
        f[8] = ":".join(fmt + ["GQ"])
        for i in range(9, len(f)):
            vals = f[i].split(":")
            if len(vals) == 1 and vals[0] == ".":
                # Sample dropped entirely; appending would only fabricate a field.
                continue
            f[i] = ":".join(vals + [vals[ridx] if ridx < len(vals) else "."])
        out.write("\t".join(f) + "\n")
        n_fixed += 1
    print(f"[backfill_rgq_to_gq] GQ filled from RGQ on {n_fixed:,} of {n_rec:,} records",
          file=sys.stderr)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Write the sample lists for the platform-stratified analyses.

Two families, and they are not the same kind of evidence:

  only_<P>   that platform's cases + every control.
             The four of these share NO cases, so they are independent estimates
             and a heterogeneity test on them means something.

  minus_<P>  every case except that platform's + every control.
             These are nested — any two of them share most of their samples — so
             they are NOT independent and no heterogeneity test may be run over
             them. Three of the four also keep ~90% of the cases, which forces
             their estimates to land on the full-cohort one no matter what; only
             the one that drops the dominant platform carries information.

Both are emitted because they answer different questions ("does each platform say
the same thing" vs "does the effect survive losing this platform"), and the report
labels which is which. Reading a forest plot of minus_* as four confirmations is
the specific mistake this file exists to prevent.
"""

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def sniff_sep(header):
    return "\t" if header.count("\t") >= header.count(",") and "\t" in header else ","


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fam", required=True, help="the analysis fam: FID IID ... SEX PHENO")
    ap.add_argument("--cram-info", required=True)
    ap.add_argument("--sample-id-col", default="ID_JHRPv6")
    ap.add_argument("--platform-col", default="WGS_Platform")
    ap.add_argument("--min-cases", type=int, default=1,
                    help="skip an only_<P> subset with fewer cases than this")
    ap.add_argument("--out-dir", default=".")
    ap.add_argument("--out-manifest", required=True)
    args = ap.parse_args()

    with open(args.cram_info) as fh:
        sep = sniff_sep(fh.readline())
    with open(args.cram_info) as fh:
        pf_of = {r[args.sample_id_col].strip(): r[args.platform_col].strip()
                 for r in csv.DictReader(fh, delimiter=sep)}

    fam = []
    for line in open(args.fam):
        f = line.split()
        if len(f) >= 6:
            fam.append((f[0], f[1], f[5]))
    if not fam:
        raise SystemExit(f"[make_platform_subsets] {args.fam} is empty")

    cases = [(fid, iid) for fid, iid, p in fam if p == "2"]
    ctrls = [(fid, iid) for fid, iid, p in fam if p == "1"]
    by_pf = defaultdict(list)
    for fid, iid in cases:
        by_pf[pf_of.get(iid, "UNKNOWN")].append((fid, iid))

    unknown = len(by_pf.get("UNKNOWN", []))
    if unknown:
        raise SystemExit(f"[make_platform_subsets] {unknown} case(s) have no platform in "
                         f"{args.cram_info}; they would be silently dropped from every "
                         f"subset. Fix the sheet rather than let them vanish.")

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    def write(tag, rows):
        p = out / f"{tag}.keep"
        p.write_text("".join(f"{fid}\t{iid}\n" for fid, iid in rows))
        return p, len(rows)

    man = []
    write("all", cases + ctrls)
    man.append(("all", "reference", len(cases), len(ctrls),
                "the whole cohort — must reproduce the published result"))

    for pf in sorted(by_pf, key=lambda k: -len(by_pf[k])):
        ids = by_pf[pf]
        safe = pf.replace(" ", "_").replace("/", "_")
        if len(ids) >= args.min_cases:
            write(f"only_{safe}", ids + ctrls)
            man.append((f"only_{safe}", "independent", len(ids), len(ctrls),
                        f"{pf} cases only — shares no case with any other only_*"))
        rest = [x for p2, v in by_pf.items() if p2 != pf for x in v]
        write(f"minus_{safe}", rest + ctrls)
        man.append((f"minus_{safe}", "nested", len(rest), len(ctrls),
                    f"every case except {pf} — overlaps the other minus_* by "
                    f"{100 * (len(rest) - 0) / len(cases):.0f}% of cases"))

    with open(args.out_manifest, "w") as fh:
        fh.write("TAG\tFAMILY\tN_CASE\tN_CONTROL\tNOTE\n")
        for t, fam_, nc, nk, note in man:
            fh.write(f"{t}\t{fam_}\t{nc}\t{nk}\t{note}\n")

    print(f"[make_platform_subsets] {len(man)} subsets from {len(cases)} cases / "
          f"{len(ctrls)} controls across {len(by_pf)} platforms")
    for t, fam_, nc, nk, _ in man:
        print(f"  {t:34s}{fam_:14s}{nc:>6} case {nk:>7} control")


if __name__ == "__main__":
    main()

"""Readers for the two external panels (shared by steps 03 and 05).

ToMMo 61KJPN: gene,allele,count,frequency; alleles named by P group. Null alleles
are dropped (ToMMo encodes DRB4 absence as DRB4*03:01N) and the rest renormalised.
1KG JPT: one row per person, two columns per locus, values like '23:01'. A trailing
'*' (a flagged call) is stripped; ambiguity strings ('01:01/04N/22N') and junk cannot
be attributed to one allele and leave numerator and denominator alike.
"""
import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as C  # noqa: E402
import nomenclature as N  # noqa: E402

KG_CLEAN = re.compile(r'^(\d+:\d+[A-Z]?)\*?$')


def tommo_counts():
    d = pd.read_csv(C.TOMMO)
    out, nulls = {}, {}
    for g in C.TOMMO_LOCI:
        cnt, nn = {}, 0
        for a, c in zip(d[d.gene == g].allele, d[d.gene == g]['count']):
            if a.rstrip('P').endswith('N'):
                nn += int(c)
                continue
            k = N.panel_pg(a)
            cnt[k] = cnt.get(k, 0) + int(c)
        out[g], nulls[g] = cnt, nn
    return out, nulls


def kg_counts():
    d = pd.read_csv(C.KG, sep='\t', dtype=str)
    d = d[d.Population == C.KG_POPULATION]
    out, dropped = {}, {}
    for g in C.KG_LOCI:
        cnt, bad = {}, 0
        for a in pd.concat([d[f'HLA-{g} 1'], d[f'HLA-{g} 2']]):
            if pd.isna(a) or str(a).strip() in ('', '-', 'NA'):
                continue
            m = KG_CLEAN.match(str(a).strip())
            if not m:
                bad += 1
                continue
            name = f'{g}*{m.group(1)}'
            if N.is_null(name):
                bad += 1
                continue
            k = N.pg_of(name)
            cnt[k] = cnt.get(k, 0) + 1
        out[g], dropped[g] = cnt, bad
    return out, dropped



def tommo_vintage():
    """The newest IPD-IMGT/HLA release among the alleles ToMMo reports: the earliest
    dictionary its typing can have used. Alleles first named after it cannot appear in
    the panel whatever their true frequency."""
    d = pd.read_csv(C.TOMMO)
    rel = [N.release_of(a.rstrip('P')) for a in d.allele]
    missing = [a for a, r in zip(d.allele, rel) if r is None]
    if missing:
        raise SystemExit(f'ABORT: ToMMo alleles absent from the IPD-IMGT history: {missing[:5]}')
    return max(rel)


def kg_vintage():
    """As tommo_vintage, for the 1KG JPT values that are single alleles."""
    d = pd.read_csv(C.KG, sep='\t', dtype=str)
    d = d[d.Population == C.KG_POPULATION]
    rel = []
    for g in C.KG_LOCI:
        for a in pd.concat([d[f'HLA-{g} 1'], d[f'HLA-{g} 2']]).dropna():
            m = KG_CLEAN.match(str(a).strip())
            if m:
                r = N.release_of(f'{g}*{m.group(1)}')
                if r is not None:
                    rel.append(r)
    return max(rel)

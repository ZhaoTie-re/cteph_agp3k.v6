"""IPD-IMGT/HLA nomenclature helpers shared by steps 01, 03 and 05.

P groups (hla_nom_p.txt, release 3.64.0, the same release as the HLA-HD dictionary):
alleles whose protein sequence is identical over the antigen recognition domain
(exons 2+3 for class I, exon 2 for class II). ToMMo names its alleles by P group, so
every comparison with it is made at this level. A P group is represented here by its
2-field name without the trailing 'P' (A*01:01P -> A*01:01).
"""
import sys
from functools import lru_cache
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as C  # noqa: E402
import hlahd_io as H  # noqa: E402


@lru_cache(maxsize=1)
def pgroup_map():
    """{'GENE*NN:NN' (no suffix): representative 'GENE*NN:NN'}; asserted single-valued."""
    out, conflict = {}, []
    for line in C.P_GROUPS.read_text().splitlines():
        if line.startswith('#') or ';' not in line:
            continue
        gene, members, pg = (line.split(';') + ['', ''])[:3]
        gene, pg = gene.rstrip('*'), pg.strip()
        for a in members.split('/'):
            _, f, _ = H.split_name(f'{gene}*{a}')
            key = f"{gene}*{':'.join(f[:2])}"
            val = f'{gene}*{pg[:-1]}' if pg.endswith('P') else key
            if out.setdefault(key, val) != val:
                conflict.append(key)
    if conflict:
        raise SystemExit(f'ABORT: {len(conflict)} 2-field names span two P groups, '
                         f'e.g. {conflict[:3]}')
    return out


def pg_of(name):
    """The P-group representative of an allele name (any depth, suffix ignored)."""
    gene, f, _ = H.split_name(name)
    key = f"{gene}*{':'.join(f[:2])}"
    return pgroup_map().get(key, key)


def is_null(name):
    """True for a non-expressed allele name (suffix N)."""
    return H.split_name(name)[2] == 'N'


def panel_pg(name):
    """A panel's allele name -> P-group representative. 'A*01:01P' -> 'A*01:01'."""
    gene, rest = name.split('*', 1)
    if rest.endswith('P'):
        return f'{gene}*{rest[:-1]}'
    return pg_of(name)


def pg_label(name):
    """The P-group name as ToMMo writes it: 'A*01:01P' for a multi-allele group,
    the plain 2-field name for an allele alone in its group."""
    gene, f, _ = H.split_name(name)
    key = f"{gene}*{':'.join(f[:2])}"
    return f'{pgroup_map()[key]}P' if key in _grouped() else key


@lru_cache(maxsize=1)
def _grouped():
    """2-field keys that belong to a P group with a name ending in 'P'."""
    out = set()
    for line in C.P_GROUPS.read_text().splitlines():
        if line.startswith('#') or ';' not in line:
            continue
        gene, members, pg = (line.split(';') + ['', ''])[:3]
        if pg.strip().endswith('P'):
            for a in members.split('/'):
                _, f, _ = H.split_name(f"{gene.rstrip('*')}*{a}")
                out.add(f"{gene.rstrip('*')}*{':'.join(f[:2])}")
    return out


@lru_cache(maxsize=1)
def first_release():
    """{'GENE*NN:NN' (no suffix): first IPD-IMGT/HLA release (e.g. 3350) in which an allele
    now bearing that 2-field name was named}, from the pinned Allelelist_history.txt.

    One row per allele (HLA_ID), one column per release, newest first ("NA" = absent from
    that release). Releases before
    3.0.0 use the old nomenclature (A*0101), so an allele is keyed by its NEWEST name and
    dated by the OLDEST release in which it appears under any name.
    """
    import csv
    out = {}
    with open(C.IMGT_HISTORY) as fh:
        rd = csv.reader(ln for ln in fh if not ln.startswith('#'))
        rel = [int(x) for x in next(rd)[1:]]               # newest first
        for r in rd:
            names = r[1:]
            present = [i for i, x in enumerate(names) if x.strip() not in ('', 'NA')]
            if not present:
                continue
            newest = names[present[0]].strip()
            if '*' not in newest:
                continue
            first = min(rel[i] for i in present)
            gene, f, _ = H.split_name(newest)
            k = f"{gene}*{':'.join(f[:2])}"
            out[k] = min(out.get(k, first), first)
    return out


def release_of(name):
    gene, f, _ = H.split_name(name)
    return first_release().get(f"{gene}*{':'.join(f[:2])}")

#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : Helpers shared by the gene_signals step scripts. No paths, no
#           project constants: every script receives its inputs, tools and
#           resources on the command line from gene_signals.nf.
# ---------------------------------------------------------------------------
import datetime as dt
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

MISSING = '—'   # em dash: "not available / not applicable"


def log(msg):
    print(f'[{dt.datetime.now():%H:%M:%S}] {msg}', flush=True)


def chrom_plain(c):
    return str(c).replace('chr', '')


def split_id(vid):
    """chr:pos:REF:ALT -> (chrom without 'chr', pos, REF, ALT)."""
    c, p, ref, alt = vid.split(':')
    return chrom_plain(c), int(p), ref, alt


def run(cmd, **kw):
    r = subprocess.run([str(c) for c in cmd], capture_output=True, text=True, **kw)
    if r.returncode != 0:
        raise RuntimeError(' '.join(map(str, cmd)) + '\n' + (r.stdout + r.stderr)[-3000:])
    return r.stdout


def read_tsv(path, **kw):
    """Every column as text; '' stays '' (never NaN). Callers convert with num()."""
    return pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False, na_values=[], **kw)


def read_table(path):
    """A table written by write_tsv(): numbers as numbers, '' as NaN."""
    return pd.read_csv(path, sep='\t', keep_default_na=False, na_values=[''], low_memory=False)


def write_tsv(df, path):
    """P columns in scientific notation; floats with 10 significant digits, so the steps
    downstream re-derive exactly what the upstream tables rounded."""
    out = df.copy()
    for c in out.columns:
        p_name = c in ('P', 'p', 'pvalue', 'min_p') or c.endswith('_pvalue') or c.startswith('p_')
        if p_name and pd.api.types.is_numeric_dtype(out[c]):
            out[c] = out[c].map(lambda x: '' if pd.isna(x) else f'{float(x):.6e}')
    out.to_csv(path, sep='\t', index=False, float_format='%.10g')


def num(s):
    return pd.to_numeric(s, errors='coerce')


def csv_list(s):
    """'a,b,c' -> ['a', 'b', 'c'] (the form in which the .nf passes lists)."""
    return [x.strip() for x in str(s).split(',') if x.strip()]


def kv_paths(items):
    """['key=path', ...] -> {key: Path}."""
    out = {}
    for it in items or []:
        k, v = it.split('=', 1)
        out[k] = Path(v)
    return out


def load_json(path):
    with open(path) as fh:
        return json.load(fh)


def import_annot(annot_lib):
    """Import the component-shared variant_annot.py from the staged file's real directory."""
    d = str(Path(annot_lib).resolve().parent)
    if d not in sys.path:
        sys.path.insert(0, d)
    import variant_annot  # noqa: E402
    return variant_annot


def load_variants(path):
    """variants.tsv as written by select_genes.py, with its numeric columns typed."""
    v = read_tsv(path)
    for c in ('gene_order', 'pos'):
        v[c] = num(v[c]).astype(int)
    return v


# ── multiple testing over a gene family (primary: >= MinNumVar; sensitivity: >= a larger minimum) ──
def bh_q(p):
    """Benjamini-Hochberg q-values (step-up, capped at 1), in the input order."""
    p = np.asarray(p, dtype=float)
    m = len(p)
    if not m:
        return p
    o = np.argsort(p, kind='stable')
    q = np.minimum.accumulate((p[o] * m / np.arange(1, m + 1))[::-1])[::-1]
    out = np.empty(m)
    out[o] = np.minimum(q, 1.0)
    return out


def family_calls(scan, min_var, alpha):
    """Per cohort x stratum x method, over the genes with >= min_var mapped variants and a P (the family): BH q,
    BH call (q < alpha), Bonferroni threshold (alpha / genes) and call (P < threshold). The gene P values are those
    of the scan: a gene's test does not depend on which other genes are in the family. One row per gene in the
    family; scan needs cohort, stratum, method, set_name, pvalue, n_var_map."""
    s = scan.copy()
    s['p'] = num(s.pvalue)
    s['nv'] = num(s.n_var_map)
    parts = []
    for _key, x in s.groupby(['cohort', 'stratum', 'method'], sort=True):
        x = x[x.p.notna() & (x.nv >= min_var)].copy()
        thr = alpha / len(x)
        x['q'] = bh_q(x.p.to_numpy())
        x['bh'] = x.q < alpha
        x['bonf_thr'] = thr
        x['bonf'] = x.p < thr
        x['family_size'] = len(x)
        parts.append(x)
    return pd.concat(parts, ignore_index=True)


def tiers_from_calls(calls, cohorts):
    """(set_name, stratum) -> cross-cohort consistency tier of the BH calls: 1, called in every cohort and by both
    statistics in at least one; 2, called in every cohort; 3, called in >= 2 cohorts; 0, otherwise."""
    out = {}
    for (sn, st), x in calls.groupby(['set_name', 'stratum']):
        called = set(x[x.bh].cohort)
        both = set(x[x.bh & (x.method == 'cmc')].cohort) & set(x[x.bh & (x.method == 'skato')].cohort)
        out[(sn, st)] = 1 if (len(called) == len(cohorts) and both) else 2 if len(called) == len(cohorts) \
            else 3 if len(called) >= 2 else 0
    return out

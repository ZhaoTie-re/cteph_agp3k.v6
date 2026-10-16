#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : Helpers shared by the saige_signals step scripts. No paths, no
#           project constants: every script receives its inputs, tools and
#           resources on the command line from saige_signals.nf.
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
        if c == 'P' or c.endswith('_P') or c.startswith('P_'):
            out[c] = out[c].map(lambda x: '' if pd.isna(x) else f'{float(x):.6e}')
    out.to_csv(path, sep='\t', index=False, float_format='%.10g')


def num(s):
    return pd.to_numeric(s, errors='coerce')


def as_bool(s):
    return s.map(lambda x: x if isinstance(x, (bool, np.bool_)) else str(x) in ('True', 'TRUE', 'true', '1'))


def sig4(x):
    """Four significant digits, trailing zeros kept (0.0281 -> 0.02810)."""
    x = float(x)
    if x == 0 or not np.isfinite(x):
        return '0' if x == 0 else str(x)
    decimals = max(0, 3 - int(np.floor(np.log10(abs(x)))))
    return f'{x:.{decimals}f}'


def csv_list(s):
    """'a,b,c' -> ['a', 'b', 'c'] (the form in which the .nf passes lists)."""
    return [x.strip() for x in str(s).split(',') if x.strip()]


def kv_paths(items):
    """['cohort=path', ...] -> {cohort: Path}."""
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
    """variants.tsv as written by select_signals.py, with its numeric columns typed."""
    v = read_tsv(path)
    for c in ('order', 'pos', 'role_rank'):
        v[c] = num(v[c]).astype(int)
    v['max_pip'] = num(v['max_pip'])
    return v

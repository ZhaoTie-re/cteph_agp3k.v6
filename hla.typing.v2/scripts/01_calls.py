#!/usr/bin/env python3
"""Step 01 — HLA-HD output -> one reconciled genotype per sample and gene.

Nothing is re-typed. For every sample and each of HLA-HD's 33 genes this step reads
the archived output and derives, in order:

1. THE CANDIDATES. `<gene>.est.txt` lists the best-scoring allele pair(s). Each side
   of a pair is a list of allele names the reads cannot tell apart, and several pairs
   may tie. The result line in `final.result.txt` is only HLA-HD's display of this:
   it shows every combination when a side's candidates differ at the 2nd or 3rd
   field, so it is not parsed for alleles.

2. THE AGREED NAME of each allele: the deepest field on which every candidate for
   that allele agrees, across all tied pairs (pairs are aligned to the reported pair
   first, because HLA-HD does not keep the two sides in a fixed order).

3. THE 2-FIELD CALL, from the first rule that applies:
     agreed   the candidates agree to >= 2 fields;
     prior    the candidates differ at 2 fields, and HLA-HD printed one allele for
              them after applying its built-in allele-frequency prior (only A, B, C,
              DPA1, DPB1, DQA1, DQB1 and DRB1 have one). Recorded with the prior's
              count ratio between the printed 2-field allele and the best other one;
     pgroup   the candidates differ at 2 fields but all lie in one IPD-IMGT P group
              (identical protein over the antigen recognition domain); named by the
              group's 2-field representative;
     none     none of these: no 2-field allele is written. Never an arbitrary pick.

4. THE STATE. Obligate genes: heterozygous, homozygous (a '-' second allele, or the
   same 2-field allele twice), ambiguous (1 field only), unresolved (no field),
   failed (no candidate, or no result line). DRB3/4/5: copy number from DRB1 (the
   DR51/52/53 rule in config.py) — present (1 or 2 copies), absent, missing (expected
   but not typed). When HLA-HD reports two alleles where DRB1 implies one copy, the
   allele with more reads unique to it is kept (`drb_reconciled`); a paralogue typed
   where DRB1 implies none is set absent (`drb_unexpected`). The raw call is kept.

5. THE READ EVIDENCE from `<gene>.read.txt`: reads assigned to each allele, and reads
   assigned to one allele but not the other ("unique reads"), which are the evidence
   that a heterozygote really carries two alleles. Thresholds are applied in step 02.

Outputs (results/01.calls/):
  calls_long.tsv          one row per sample x gene (33 genes; `reported` marks the 19)
  genotypes.tsv           wide, reported genes only: <g>_1, <g>_2, <g>_state
  tie_candidates.tsv      every call whose best pairs tie, with the pairs listed
  drb_reconciliation.tsv  DRB1-implied vs HLA-HD copy number, per sample and paralogue
"""
import sys
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as C  # noqa: E402
import hlahd_io as H  # noqa: E402
import nomenclature as N  # noqa: E402

PRIOR_GENES = ('A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRB1')


# ------------------------------------------------------------------ references
def load_prior():
    """{gene: {'NN:NN' 2-field: count}} from HLA-HD's freq_data (pinned copy)."""
    d = C.D_INPUTS / 'hlahd_freq_data'
    if not d.exists():
        if not TEST:
            raise SystemExit(f'ABORT: {d} missing; run 00_provenance.py first')
        d = C.HLAHD_FREQ
    out = {}
    for g in PRIOR_GENES:
        cnt = Counter()
        for ln in (d / f'{g}_count.txt').read_text().splitlines():
            f = ln.split('\t')
            if len(f) == 2:
                cnt[H.trim(f[0], 2)] += int(f[1])
        out[g] = cnt
    return out


PRIOR = None
TEST = False


def _init():
    global PRIOR
    N.pgroup_map()
    PRIOR = load_prior()


# ------------------------------------------------------------------ helpers
pg_of = N.pg_of


def depth(name):
    return len(H.split_name(name)[1]) if name else 0


def cp(names):
    return H.common_prefix(names) if names else ('', 0)


def align(pair, ref):
    """Orient `pair` (s1, s2) to `ref` (r1, r2) by the deeper total shared prefix."""
    s1, s2 = pair
    r1, r2 = ref
    keep = cp(r1 + s1)[1] + cp(r2 + s2)[1]
    swap = cp(r1 + s2)[1] + cp(r2 + s1)[1]
    return (s2, s1) if swap > keep else (s1, s2)


def prior_count(gene, name):
    """HLA-HD's prior count for a 2-field name (0 when the gene has no prior)."""
    return PRIOR.get(gene, {}).get(H.trim(name, 2), 0)


def side_has(side, name):
    """True if `name` (as printed in the result line) is a field-prefix of a candidate."""
    g, f, _ = H.split_name(name)
    for c in side:
        cg, cf, _ = H.split_name(c)
        if cf[:len(f)] == f:
            return True
    return False


def side_call(gene, cands, shown):
    """(2-field name, basis, agreed name, agreed depth, prior ratio) for one allele.

    `shown` is the name HLA-HD printed for this allele when it printed exactly one
    pair, i.e. after its own prior resolved any tie; None otherwise.
    """
    agreed, d = cp(cands)
    if d >= 2:
        return H.trim(agreed, 2), 'agreed', agreed, d, ''
    if shown and gene in PRIOR_GENES and depth(shown) >= 2:
        name2 = H.trim(shown, 2)
        others = {H.trim(c, 2) for c in cands} - {name2}
        w_other = max((prior_count(gene, o) for o in others), default=0)
        ratio = (prior_count(gene, name2) + 0.5) / (w_other + 0.5)
        return name2, 'prior', agreed, d, ratio
    groups = {pg_of(n) for n in cands}
    if len(groups) == 1:
        return next(iter(groups)), 'pgroup', agreed, d, ''
    return '', 'none', agreed, d, ''


# ------------------------------------------------------------------ per sample
def derive(path):
    sid = path.name[:-len('.tar.gz')]
    so = H.read_archive(path, sid)
    rows, ties = [], []
    for gene in GENES:
        r = {'sample_id': sid, 'gene': gene, 'reported': int(gene in C.REPORTED_GENES)}
        ent = so.final.get(gene)
        est = so.est.get(gene)
        reads = so.reads.get(gene)
        r['raw_result'] = ';'.join(H.strip_hla(x) for x in ent) if ent else ''
        if ent is None or est is None:
            r.update(raw='missing')
            rows.append(r)
            continue
        if est.no_candidate or (ent and ent[0] == H.NOT_TYPED) or not est.best:
            r.update(raw='not_typed', n_best=0)
            rows.append(r)
            continue
        best = est.best[:max(est.n_best, 1)]
        # The reported pair: the best pair whose sides contain the two names HLA-HD
        # printed, when it printed exactly one pair. Oriented to the printed order.
        shown = [H.strip_hla(x) for x in ent]
        shown_pair = None
        ref = best[0]
        if len(shown) == 2:
            sp = (shown[0], shown[1] if shown[1] != '-' else shown[0])
            for p in best:
                s1, s2 = p[0], p[1] if p[1] is not None else p[0]
                if side_has(s1, sp[0]) and side_has(s2, sp[1]):
                    ref, shown_pair = (s1, s2, p[1] is None), sp
                    break
                if side_has(s2, sp[0]) and side_has(s1, sp[1]):
                    ref, shown_pair = (s2, s1, p[1] is None), sp
                    break
        if len(ref) == 2:
            ref = (ref[0], ref[1] if ref[1] is not None else ref[0], ref[1] is None)
        homo_model = ref[2]
        ref_full = (ref[0], ref[1])
        aligned = [align((p[0], p[1] if p[1] is not None else p[0]), ref_full) for p in best]
        side1 = sorted({n for p in aligned for n in p[0]})
        side2 = sorted({n for p in aligned for n in p[1]})
        tie = len(best) > 1
        a1, b1, ag1, d1, q1 = side_call(gene, side1, shown_pair[0] if shown_pair else None)
        a2, b2, ag2, d2, q2 = side_call(gene, side2, shown_pair[1] if shown_pair else None)
        same = (shown_pair[0] == shown_pair[1]) if shown_pair else (cp(side1)[0] == cp(side2)[0])
        r.update(raw=('one' if homo_model else 'same' if same else 'two'),
                 n_best=len(best), n_other=len(est.other), n_shown=len(shown),
                 homozygous_model=int(homo_model),
                 agreed_1=ag1, agreed_depth_1=d1, agreed_2=ag2, agreed_depth_2=d2,
                 allele_1=a1, basis_1=b1, allele_2=a2, basis_2=b2,
                 prior_ratio_1=q1, prior_ratio_2=q2,
                 n_cand_1=len(side1), n_cand_2=len(side2))
        if tie:
            ties.append({'sample_id': sid, 'gene': gene, 'n_best': len(best),
                         'pairs': ' | '.join(f'{cp(p[0])[0] or ",".join(p[0][:3])} + '
                                             f'{cp(p[1])[0] or ",".join(p[1][:3])}'
                                             for p in aligned),
                         'reported': f'{a1} + {a2}', 'basis': f'{b1},{b2}'})
        # read evidence: the blocks belonging to the reported pair
        if reads and reads.blocks:
            s1set, s2set = set(ref_full[0]), set(ref_full[1])
            i1 = next((i for i, b in enumerate(reads.blocks) if b[0] in s1set), None)
            i2 = (i1 if homo_model else
                  next((i for i, b in enumerate(reads.blocks) if b[0] in s2set and i != i1), None))
            if i1 is not None:
                r['reads_1'] = reads.blocks[i1][1]
            if i2 is not None:
                r['reads_2'] = reads.blocks[i2][1]
            if i1 is not None and i2 is not None and i1 != i2:
                n1, n2 = reads.blocks[i1][2], reads.blocks[i2][2]
                r['unique_1'], r['unique_2'], r['shared'] = len(n1 - n2), len(n2 - n1), len(n1 & n2)
        rows.append(r)
    return sid, rows, ties, so.n_fail_lines


# ------------------------------------------------------------------ states
def obligate_state(r):
    if r['raw'] in ('missing', 'not_typed'):
        return 'failed', 0
    a1, a2 = r['allele_1'], r['allele_2']
    if a1 and a2:
        return ('homozygous' if a1 == a2 else 'heterozygous'), 2
    if not a1 and not a2 and r['agreed_depth_1'] == 0 and r['agreed_depth_2'] == 0:
        return 'unresolved', 2
    return 'ambiguous', 2          # at least one allele short of 2 fields


def expected_copies(drb1_row):
    """{paralogue: expected copies} from the two DRB1 alleles' first field, or None."""
    if drb1_row is None or drb1_row['state'] in ('failed', 'unresolved'):
        return None
    exp = {p: 0 for p in C.PARALOGUES}
    for k in ('1', '2'):
        name = drb1_row[f'allele_{k}'] or drb1_row[f'agreed_{k}']
        if not name:
            return None
        first = H.split_name(name)[1][0]
        if first not in C.DR_GROUP_PARALOGUE:
            return None
        par = C.DR_GROUP_PARALOGUE[first]
        if par:
            exp[par] += 1
    return exp


def paralogue_state(r, e):
    """Final state, copies, alleles and flag for DRB3/4/5 given expected copies e."""
    raw = r['raw']
    obs = 0 if raw in ('missing', 'not_typed') else 1 if raw in ('one', 'same') else 2
    out = {'drb_expected': '' if e is None else e, 'drb_observed': obs, 'drb_flag': ''}
    a1, a2 = r.get('allele_1', ''), r.get('allele_2', '')
    if e is None:                                    # DRB1 cannot say: keep HLA-HD's call
        out['drb_flag'] = 'drb1_unknown'
        e = obs
    if e == 0:
        if obs:
            out['drb_flag'] = out['drb_flag'] or 'drb_unexpected'
        out.update(state='absent', copies=0, allele_1='', allele_2='')
        return out
    if obs == 0:
        out.update(state='missing', copies='', allele_1='', allele_2='',
                   drb_flag=out['drb_flag'] or 'drb_missing')
        return out
    if e == 1:
        if obs == 2:
            def num(x):
                x = pd.to_numeric(x, errors='coerce')
                return 0 if pd.isna(x) else x
            u1, u2 = num(r.get('unique_1')), num(r.get('unique_2'))
            if (u2, num(r.get('reads_2'))) > (u1, num(r.get('reads_1'))):
                a1 = a2
                out['kept_side'] = 2
            out['drb_flag'] = out['drb_flag'] or 'drb_reconciled'
        a2 = ''
        if not a1:
            out.update(state='ambiguous' if r.get('agreed_depth_1') else 'unresolved',
                       copies=1, allele_1='', allele_2='')
            return out
        out.update(state='present', copies=1, allele_1=a1, allele_2='')
        return out
    # e == 2
    if obs == 1:
        a2 = a1
    if not (a1 and a2):
        out.update(state='ambiguous', copies=2, allele_1=a1, allele_2=a2)
        return out
    out.update(state='present', copies=2, allele_1=a1, allele_2=a2)
    return out


GENES = None


def main():
    global GENES, TEST
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--test', type=int, default=0,
                    help='first N archives present, outputs to --out-dir (development only)')
    ap.add_argument('--out-dir', default=None)
    args = ap.parse_args()
    TEST = bool(args.test)
    out_dir = Path(args.out_dir) if args.out_dir else C.D_CALLS
    GENES = [ln.split()[0] for ln in C.HLAHD_SPLIT.read_text().splitlines()
             if ln.strip() and not ln.startswith('#')]
    if not set(C.REPORTED_GENES) <= set(GENES):
        raise SystemExit('ABORT: a reported gene is not in the gene-split file')
    if TEST:
        ids = sorted(p.name[:-7] for p in (C.D_INPUTS / 'typing').glob('*.tar.gz'))[:args.test]
    else:
        ids = pd.read_csv(C.D_INPUTS / 'sample_inputs.tsv', sep='\t', dtype=str).sample_id.tolist()
    paths = [C.D_INPUTS / 'typing' / f'{s}.tar.gz' for s in ids]

    rows, ties, fails = [], [], 0
    with ProcessPoolExecutor(8, initializer=_init, initargs=()) as ex:
        for sid, rr, tt, nf in ex.map(derive, paths, chunksize=16):
            rows.extend(rr)
            ties.extend(tt)
            fails += nf
    _init()
    df = pd.DataFrame(rows)
    for c in ('allele_1', 'allele_2', 'agreed_1', 'agreed_2', 'basis_1', 'basis_2'):
        df[c] = df[c].fillna('')
    for c in ('agreed_depth_1', 'agreed_depth_2'):
        df[c] = df[c].fillna(0).astype(int)

    # obligate and pseudogene states
    df['state'], df['copies'] = '', 2
    for i, r in df.iterrows():
        if r.gene in C.PARALOGUES:
            continue
        s, c = obligate_state(r)
        df.at[i, 'state'], df.at[i, 'copies'] = s, c
    # paralogues from DRB1
    drb1 = df[df.gene == 'DRB1'].set_index('sample_id')
    recon = []
    df['drb_expected'], df['drb_observed'], df['drb_flag'] = '', '', ''
    for i, r in df[df.gene.isin(C.PARALOGUES)].iterrows():
        e_all = expected_copies(drb1.loc[r.sample_id] if r.sample_id in drb1.index else None)
        e = None if e_all is None else e_all[r.gene]
        o = paralogue_state(r, e)
        recon.append({'sample_id': r.sample_id, 'gene': r.gene,
                      'drb1_1': drb1.at[r.sample_id, 'allele_1'] or drb1.at[r.sample_id, 'agreed_1'],
                      'drb1_2': drb1.at[r.sample_id, 'allele_2'] or drb1.at[r.sample_id, 'agreed_2'],
                      'expected': o['drb_expected'], 'observed': o['drb_observed'],
                      'raw_result': r.raw_result, 'state': o['state'],
                      'allele_1': o['allele_1'], 'allele_2': o['allele_2'], 'flag': o['drb_flag']})
        kept = o.pop('kept_side', 1)
        for k, v in o.items():
            df.at[i, k] = v
        if kept == 2:     # the evidence columns follow the allele that was kept
            for a, b in (('reads_1', 'reads_2'), ('unique_1', 'unique_2')):
                df.at[i, a], df.at[i, b] = r[b], r[a]
    # the raw HLA-HD 2-field call stays beside the reconciled one
    df.loc[df.gene.isin(C.PARALOGUES), 'raw_allele_1'] = \
        df.loc[df.gene.isin(C.PARALOGUES)].index.map(lambda i: rows[i].get('allele_1', ''))
    df.loc[df.gene.isin(C.PARALOGUES), 'raw_allele_2'] = \
        df.loc[df.gene.isin(C.PARALOGUES)].index.map(lambda i: rows[i].get('allele_2', ''))
    df['p_group_1'] = [pg_of(a) if a else '' for a in df.allele_1]
    df['p_group_2'] = [pg_of(a) if a else '' for a in df.allele_2]

    out_dir.mkdir(parents=True, exist_ok=True)
    cols = ['sample_id', 'gene', 'reported', 'state', 'copies', 'allele_1', 'allele_2',
            'basis_1', 'basis_2', 'p_group_1', 'p_group_2', 'agreed_1', 'agreed_depth_1',
            'agreed_2', 'agreed_depth_2', 'raw', 'raw_result', 'homozygous_model',
            'n_best', 'n_other', 'n_shown', 'n_cand_1', 'n_cand_2', 'prior_ratio_1', 'prior_ratio_2',
            'reads_1', 'reads_2', 'unique_1', 'unique_2', 'shared',
            'drb_expected', 'drb_observed', 'drb_flag', 'raw_allele_1', 'raw_allele_2']
    df = df.reindex(columns=cols)
    df.to_csv(out_dir / 'calls_long.tsv', sep='\t', index=False)

    rep = df[df.reported == 1]
    wide = rep.pivot(index='sample_id', columns='gene',
                     values=['allele_1', 'allele_2', 'state'])
    wide = pd.concat({g: pd.DataFrame({f'{g}_1': wide[('allele_1', g)],
                                       f'{g}_2': wide[('allele_2', g)],
                                       f'{g}_state': wide[('state', g)]})
                      for g in C.REPORTED_GENES}, axis=1)
    wide.columns = wide.columns.droplevel(0)
    wide.reindex(ids).to_csv(out_dir / 'genotypes.tsv', sep='\t')
    pd.DataFrame(ties, columns=['sample_id', 'gene', 'n_best', 'pairs', 'reported', 'basis']) \
        .to_csv(out_dir / 'tie_candidates.tsv', sep='\t', index=False)
    pd.DataFrame(recon).to_csv(out_dir / 'drb_reconciliation.tsv', sep='\t', index=False)

    n = len(ids)
    if len(df) != n * len(GENES):
        raise SystemExit(f'ABORT: {len(df)} rows, expected {n} x {len(GENES)}')
    print(f'[calls] {n} samples x {len(GENES)} genes; unreadable result lines: {fails}')
    print(rep.groupby('gene').state.value_counts().unstack(fill_value=0)
          .reindex(C.REPORTED_GENES).to_string())
    b = pd.concat([rep.basis_1, rep.basis_2])
    print('[calls] 2-field basis over reported alleles:', b[b != ''].value_counts().to_dict())
    print('[calls] paralogue flags:', pd.DataFrame(recon).groupby('gene').flag.value_counts().to_dict())


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""Step 05 — the tables the association component reads. Produced, not assessed here.

File names and layouts follow what analysis/assoc_hla reads from hla.typing:
  allele_dosage.tsv      sample x 2-field allele, 0/1/2, blank = no genotype
  residue_dosage.tsv     sample x (gene:IMGT position:residue), 0/1/2, blank = no genotype
  allele_pgroup_map.tsv  allele -> P group, whether ToMMo carries it, ToMMo frequency
plus residue_sites.tsv, allele_match_map.tsv, allele_unmatched.tsv and
call_confidence.tsv (sample x gene, from step 02) so the association can filter.

What changed against hla.typing's exports:
  - the 19 reported genes only (no pseudogenes);
  - DRB3/4/5 follow the DRB1-based copy number: a single copy is dosage 1 and an
    absent gene is dosage 0 for every allele and residue of that gene (the earlier
    tables wrote a single copy as two);
  - a call without a 2-field allele is blank, never an allele trimmed to one field;
  - a residue the alignment does not determine for one of the sample's alleles ('*' unsequenced,
    '.' gap) leaves that sample blank at every residue of the position. Written as 0 it would
    claim the sample carries none of them; blank, the matrix can be read as it stands.

Residues come from IPD-IMGT/HLA's own protein alignments (build_residue_ref.py, copied
unchanged from hla.typing). An allele is matched to its own row of the alignment, else to
its first deeper allele (same two fields, so the same protein). An allele with neither has no
amino acids to read: its carriers are blank at every residue of the gene, and it is listed in
allele_unmatched.tsv. It is never matched to a shorter name, which would borrow another
protein's sequence. A position is exported
when at least two residues are determined among the samples; one column per residue seen.
"""
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as C  # noqa: E402
import nomenclature as N  # noqa: E402
from panels import tommo_counts  # noqa: E402

CALLED = ('heterozygous', 'homozygous', 'present')
MISSING = '.'


def carried(r):
    """The 2-field alleles a sample carries at one gene; [] if absent; None if unknown."""
    if r.state == 'absent':
        return []
    if r.state not in CALLED:
        return None
    al = [r.allele_1] if r.copies == 1 else [r.allele_1, r.allele_2]
    return None if any(not a for a in al) else al


def build_index(ref):
    exact = {a: a for a in ref.index}
    by_prefix = defaultdict(list)
    for a in ref.index:
        f = a.split(':')
        for d in range(1, len(f) + 1):
            by_prefix[':'.join(f[:d])].append(a)
    return exact, by_prefix


def resolve(allele, exact, by_prefix):
    if allele in exact:
        return allele, 'exact'
    if allele in by_prefix:
        c = by_prefix[allele]
        return sorted(c)[0], f'prefix({len(c)})'
    return None, 'no_match'


def main():
    out = C.D_EXPORTS
    out.mkdir(parents=True, exist_ok=True)
    ref_dir = out / 'residue_reference'
    subprocess.run([C.PYTHON, str(Path(__file__).with_name('build_residue_ref.py')),
                    '--alignment-dir', str(C.IMGT_ALIGN), '--genes', ','.join(C.REPORTED_GENES),
                    '--out-dir', str(ref_dir),
                    '--out-summary', str(out / 'residue_reference_summary.tsv')], check=True)

    ev = pd.read_csv(C.D_EVIDENCE / 'call_evidence.tsv', sep='\t', dtype=str).fillna('')
    ev['copies'] = pd.to_numeric(ev.copies, errors='coerce')
    samples = pd.read_csv(C.D_INPUTS / 'sample_inputs.tsv', sep='\t', dtype=str).sample_id.tolist()

    allele_cols, residue_cols, sites, match_rows, unmatched = {}, {}, [], [], []
    for g in C.REPORTED_GENES:
        s = ev[ev.gene == g].set_index('sample_id').reindex(samples)
        car = {sid: carried(r) for sid, r in zip(samples, s.itertuples())}
        # alleles
        seen = sorted({a for v in car.values() if v for a in v})
        for a in seen:
            allele_cols[a] = [('' if car[x] is None else car[x].count(a)) for x in samples]
        # residues
        rp = ref_dir / f'{g}.tsv'
        if not rp.exists():
            sites.append({'gene': g, 'positions': 0, 'polymorphic': 0, 'blank_undetermined': 0,
                          'note': 'no IMGT protein alignment'})
            continue
        ref = pd.read_csv(rp, sep='\t', index_col=0, dtype=str)
        exact, by_prefix = build_index(ref)
        rows = {}
        for a in seen:
            m, how = resolve(a, exact, by_prefix)
            n_chr = sum(v.count(a) for v in car.values() if v)
            match_rows.append({'gene': g, 'allele': a, 'matched_to': m or '', 'how': how,
                               'n_chromosomes': n_chr})
            if m is None:
                unmatched.append({'gene': g, 'allele': a, 'n_chromosomes': n_chr})
            else:
                rows[a] = ref.loc[m].tolist()
        per = {}
        for x in samples:
            v = car[x]
            if v is None or any(a not in rows for a in v):
                per[x] = None
            else:
                per[x] = [rows[a] for a in v]
        n_poly = n_blank = 0
        for j, pos in enumerate(ref.columns):
            vals = {}
            for x in samples:
                v = None if per[x] is None else [r[j] for r in per[x]]
                if v is not None and any(aa in (MISSING, '*') for aa in v):
                    v = 'undetermined'                   # the alignment cannot say what this sample carries
                vals[x] = v
            obs = {aa for v in vals.values() if isinstance(v, list) for aa in v}
            if len(obs) < 2:
                continue
            n_poly += 1
            n_blank += sum(v == 'undetermined' for v in vals.values())
            for aa in sorted(obs):
                residue_cols[f'{g}:{pos}:{aa}'] = [(vals[x].count(aa) if isinstance(vals[x], list) else '')
                                                   for x in samples]
        sites.append({'gene': g, 'positions': len(ref.columns), 'polymorphic': n_poly,
                      'blank_undetermined': n_blank, 'note': ''})

    def write(name, cols):
        df = pd.DataFrame(cols, index=samples)
        df.index.name = 'sample_id'
        df.to_csv(out / name, sep='\t')
        return df.shape

    sa = write('allele_dosage.tsv', allele_cols)
    sr = write('residue_dosage.tsv', residue_cols)
    pd.DataFrame(sites).to_csv(out / 'residue_sites.tsv', sep='\t', index=False)
    pd.DataFrame(match_rows).to_csv(out / 'allele_match_map.tsv', sep='\t', index=False)
    pd.DataFrame(unmatched, columns=['gene', 'allele', 'n_chromosomes']).to_csv(
        out / 'allele_unmatched.tsv', sep='\t', index=False)

    tommo, _ = tommo_counts()
    freq = {}
    for g, cnt in tommo.items():
        tot = sum(cnt.values())
        freq.update({k: v / tot for k, v in cnt.items()})
    pg = []
    for a in allele_cols:
        g = a.split('*')[0]
        rep = N.pg_of(a)
        pg.append({'gene': g, 'allele': a, 'p_group': N.pg_label(a),
                   'in_reference': int(g in tommo and rep in freq),
                   'freq_reference': round(freq[rep], 6) if g in tommo and rep in freq else ''})
    pd.DataFrame(pg).to_csv(out / 'allele_pgroup_map.tsv', sep='\t', index=False)

    conf = ev.pivot(index='sample_id', columns='gene', values='confidence') \
             .reindex(index=samples, columns=C.REPORTED_GENES)
    conf.to_csv(out / 'call_confidence.tsv', sep='\t')
    print(f'[exports] allele dosage {sa}, residue dosage {sr}, '
          f'{len(unmatched)} allele(s) without an IMGT alignment match -> {out}')


if __name__ == '__main__':
    main()

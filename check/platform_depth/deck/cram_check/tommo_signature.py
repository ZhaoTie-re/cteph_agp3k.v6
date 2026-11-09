#!/usr/bin/env python3
"""ToMMo 61KJPN frequencies of the two alleles that give the ADARB2 carrier signature (user 2026-11-07: "1～10 按照你认为
最合理的开始修改", item 10; review round 12, B-12). The 16-bp signature that marks a carrier's reads spans the junction of
the one-unit deletion (3 bp at 1747168 + 28 bp at 1747172, the cohort's lead CGTG>C at 1747167 and the 28-bp deletion
at 1747171); a G>A in a full-length unit (1747168 or 1747199) gives the same 16 bases. The Japanese reference panel
tells which of the two the signature's carrier rate fits.

Reads the ToMMo 61KJPN GF release (GRCh38, autosomes; multi-allelic and not left-trimmed, so the slice is split and
normalized with bcftools against hs38DH.fa) and writes tommo_signature.tsv: one row per allele with AC, AN, AF, FILTER
and the carrier share expected under Hardy-Weinberg, 1 - (1 - AF)^2. Allele counts only; no sample data.
Usage: python3 tommo_signature.py"""
import io
import os
import subprocess
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
BCF = '/home/b/b37974/bcftools/bcftools'
FASTA = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/nagasaki_pipeline/data/hs38DH.fa'
TOMMO = '/LARGE0/gr10478/pub/ToMMo/20250512/src/tommo-61kjpn-20250616-GRCh38-snvindel-gf-autosome-chr10.vcf.gz'
REGION = 'chr10:1747100-1747260'
WANT = {(1747167, 'CGTG', 'C'): 'deletion_3bp', (1747171, None, 'A'): 'deletion_28bp',
        (1747168, 'G', 'A'): 'ga_1747168', (1747199, 'G', 'A'): 'ga_1747199'}

view = subprocess.run([BCF, 'view', '-r', REGION, TOMMO], capture_output=True, check=True)
norm = subprocess.run([BCF, 'norm', '-m', '-any', '-f', FASTA, '-Ou', '-'], input=view.stdout, capture_output=True,
                      check=True)
q = subprocess.run([BCF, 'query', '-f', '%POS\t%REF\t%ALT\t%FILTER\t%INFO/AC\t%INFO/AN\t%INFO/AF\n', '-'],
                   input=norm.stdout, capture_output=True, check=True)
d = pd.read_csv(io.StringIO(q.stdout.decode()), sep='\t', header=None,
                names=['pos', 'ref', 'alt', 'filter', 'ac', 'an', 'af'])
rows = []
for (pos, ref, alt), name in WANT.items():
    m = d[(d.pos == pos) & (d.alt == alt) & ((d.ref == ref) if ref else (d.ref.str.len() == 29))]
    if len(m) != 1:
        raise SystemExit(f'{name}: {len(m)} records at {pos}')
    r = m.iloc[0]
    af = float(r.ac) / float(r.an)
    rows.append({'allele': name, 'pos': pos, 'ref': r.ref if len(r.ref) <= 8 else f'{r.ref[:8]}…({len(r.ref)} bp)',
                 'alt': alt, 'filter': r['filter'], 'ac': int(r.ac), 'an': int(r.an), 'af': af,
                 'carriers_expected': 1 - (1 - af) ** 2})
out = pd.DataFrame(rows)
out.to_csv(os.path.join(HERE, 'tommo_signature.tsv'), sep='\t', index=False)
print(out.to_string(index=False))

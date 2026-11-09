#!/usr/bin/env python3
"""Tables of the HaplotypeCaller census at the ADARB2 lead (README.md). Read only; run after run_hc_census.sh.

Per genome of samples.tsv:
- assembly over the lead (debug log): the assembly region overlapping the lead window, its reads, the k-mer sizes used,
  the number of candidate haplotypes;
- deletion assembled (bamout, ALL_POSSIBLE_HAPLOTYPES): a candidate haplotype (read group ArtificialHaplotypeRG) whose
  alignment removes exactly one 31-bp unit inside the minisatellite (deleted minus inserted bases within chr10:1746881-
  1747271 equal to 31), wherever the caller places the gap (the pilot found 3D at 1747168 + 28D at 1747172);
- the genotypes the re-run writes in its GVCF at the lead (1747167), the 28-bp deletion (1747171) and the two G>A
  positions of the signature (1747168, 1747199), with AD and GQ;
- whether the re-run reproduces the cohort's call at the lead (carrier or not).

Writes hc_census.tsv (INTERNAL: sample IDs; git-ignored) and hc_census_groups.tsv (counts by group, no IDs), and prints
the group table. Usage: python3 census_tables.py"""
import os
import re
import pandas as pd
import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
CHROM = 'chr10'
LEAD = (1747150, 1747210)                  # the lead window of ../hc_debug/parse_hc_log.py
VNTR = (1746881, 1747271)                  # stats.json cram.vntr (1-based, inclusive)
UNIT = 31
SITES = {1747167: 'lead', 1747171: 'del28', 1747168: 'ga_1747168', 1747199: 'ga_1747199'}
P150 = {'NovaSeq 30x', 'DNBSeq-G400RS 30x', 'DNBseq-G400RS 15x'}


def gt_str(s):
    a = s.get('GT')
    return '/'.join('.' if x is None else str(x) for x in a) if a else './.'


def log_region(path):
    """the assembly region over the lead window: reads, k-mers used, candidate haplotypes"""
    blocks, cur = [], None
    for line in open(path):
        m = re.search(r'Assembling (chr10):(\d+)-(\d+) with (\d+) reads', line)
        if m:
            cur = {'start': int(m.group(2)), 'end': int(m.group(3)), 'reads': int(m.group(4)), 'used': [], 'n_hap': None}
            blocks.append(cur)
            continue
        if cur is None:
            continue
        if (m := re.search(r'Using kmer size of (\d+)', line)):
            cur['used'].append(int(m.group(1)))
        elif (m := re.search(r'Found (\d+) candidate haplotypes', line)):
            cur['n_hap'] = int(m.group(1))
    over = [b for b in blocks if b['start'] <= LEAD[1] and b['end'] >= LEAD[0]]
    return {'lead_region_active': int(bool(over)), 'region_reads': sum(b['reads'] for b in over),
            'kmers_used': ';'.join(sorted({str(k) for b in over for k in b['used']})),
            'candidate_haplotypes': sum(b['n_hap'] or 0 for b in over)}


def deletion_haplotypes(path):
    """candidate haplotypes over the minisatellite, and those that remove exactly one unit inside it"""
    n_hap = n_del = 0
    forms = set()
    with pysam.AlignmentFile(path) as b:
        for r in b.fetch(CHROM, VNTR[0] - 1, VNTR[1]):
            if not (r.has_tag('RG') and r.get_tag('RG') == 'ArtificialHaplotypeRG'):
                continue
            n_hap += 1
            pos, net, gaps = r.reference_start + 1, 0, []
            for op, ln in r.cigartuples:
                if op == 2:                                       # deletion: reference bases pos .. pos + ln - 1
                    inside = max(0, min(pos + ln - 1, VNTR[1]) - max(pos, VNTR[0]) + 1)
                    net += inside
                    if inside:
                        gaps.append(f'{pos}:{ln}D')
                elif op == 1:                                     # insertion before reference base pos
                    if VNTR[0] <= pos <= VNTR[1]:
                        net -= ln
                        gaps.append(f'{pos}:{ln}I')
                if op in (0, 2, 3, 7, 8):
                    pos += ln
            if net == UNIT:
                n_del += 1
                forms.add(' '.join(gaps))
    return {'haplotypes_over_repeat': n_hap, 'deletion_haplotypes': n_del, 'deletion_assembled': int(n_del > 0),
            'deletion_forms': '; '.join(sorted(forms))}


S = pd.read_csv(os.path.join(HERE, 'samples.tsv'), sep='\t', dtype={'sample_id': str})
rows = []
for _, s in S.iterrows():
    o = os.path.join(HERE, 'out', s.sample_id, s.sample_id)
    if not os.path.exists(f'{o}.g.vcf.gz.tbi'):
        raise SystemExit(f'not finished: {s.sample_id}')
    r = {'sample_id': s.sample_id, 'role': s.role, 'platform': s.platform, 'status': s.status,
         'signature_reads': int(s.ad_alt), 'cohort_call': s.ad_call}
    r.update(log_region(f'{o}.hc.log'))
    r.update(deletion_haplotypes(f'{o}.bamout.bam'))
    with pysam.VariantFile(f'{o}.g.vcf.gz') as vf:
        smp = list(vf.header.samples)[0]
        for pos, name in SITES.items():
            hit = [x for x in vf.fetch(CHROM, pos - 1, pos) if x.pos == pos and any(a != '<NON_REF>' for a in x.alts)]
            if hit:
                x, sx = hit[0], hit[0].samples[smp]
                r[f'{name}_alleles'] = f"{x.ref}>{','.join(a for a in x.alts if a != '<NON_REF>')}"
                r[f'{name}_gt'] = gt_str(sx)
                r[f'{name}_ad'] = ','.join(map(str, sx.get('AD') or []))
                r[f'{name}_gq'] = sx.get('GQ')
            else:
                r[f'{name}_alleles'], r[f'{name}_gt'], r[f'{name}_ad'], r[f'{name}_gq'] = '', 'ref block', '', ''
    rows.append(r)
G = pd.DataFrame(rows)
alt = lambda g: g.astype(str).str.contains(r'[1-9]')                       # noqa: E731  any ALT index (0/1, 1/0, 1/1, 0|1)
G['rerun_lead_carrier'] = alt(G.lead_gt).astype(int)
G['rerun_ga_carrier'] = (alt(G.ga_1747168_gt) | alt(G.ga_1747199_gt)).astype(int)
G['reproduces_cohort'] = (G.rerun_lead_carrier == alt(G.cohort_call).astype(int)).astype(int)


def group(r):
    if r.role == 'non-carrier, random':
        return 'Non-carriers, DNBSeq-T7' if r.platform == 'DNBSeq-T7 30x' else 'Non-carriers, 150-bp platforms'
    if r.role == 'called, no signature':
        return 'Called, no signature read'
    if r.platform == 'DNBSeq-T7 30x':
        return 'DNBSeq-T7 read carriers'
    side = 'HiSeq X controls' if r.platform == 'HiSeqX 15x' else '150-bp cases'
    return f"{side}, {'called' if r.role.endswith('called') else 'missed'}"


G['group'] = G.apply(group, axis=1)
G.to_csv(os.path.join(HERE, 'hc_census.tsv'), sep='\t', index=False)
ORDER = ['HiSeq X controls, called', '150-bp cases, called', 'HiSeq X controls, missed', '150-bp cases, missed',
         'DNBSeq-T7 read carriers', 'Called, no signature read', 'Non-carriers, 150-bp platforms', 'Non-carriers, DNBSeq-T7']
out = []
for g in ORDER + [f'platform: {p}' for p in sorted(S.platform.unique())]:
    y = G[G.group == g] if not g.startswith('platform: ') else G[(G.platform == g[10:]) & (G.role != 'non-carrier, random')]
    out.append({'group': g, 'n': len(y), 'deletion_assembled': int(y.deletion_assembled.sum()),
                'rerun_lead_carrier': int(y.rerun_lead_carrier.sum()), 'rerun_ga_carrier': int(y.rerun_ga_carrier.sum()),
                'reproduces_cohort': int(y.reproduces_cohort.sum()), 'lead_region_active': int(y.lead_region_active.sum()),
                'signature_reads_median': float(y.signature_reads.median()) if len(y) else None,
                'region_reads_median': float(y.region_reads.median()) if len(y) else None,
                'kmers_used': ';'.join(sorted({k for x in y.kmers_used for k in str(x).split(';') if k and k != 'nan'}))})
T = pd.DataFrame(out)
T.to_csv(os.path.join(HERE, 'hc_census_groups.tsv'), sep='\t', index=False)
print(T.to_string(index=False))
print('deletion forms:', G.deletion_forms[G.deletion_forms != ''].value_counts().head(5).to_dict())

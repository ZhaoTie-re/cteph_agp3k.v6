#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : Emit the rvtest --setFile for one stratum FROM ONE IN-MEMORY MAP,
#           and, from the SAME dict, the per-gene allele and carrier counts the
#           publication table is built from -- in ONE task.
#
#           WHY THE VARIANT IDS ARE READ BACK FROM THE BIM. The set file is
#           written over variant IDs read from the .bim plink2 ACTUALLY WROTE,
#           not from the --extract list we asked for. If plink2 dropped a
#           variant (a duplicate ID, a REF mismatch), set(map) != set(bim) and
#           this aborts rather than handing rvtest a set the denominator did
#           not count.
#
#           WHY THE CHROMOSOME STRING IS READ, NOT ASSUMED. rvtest's
#           RangeCollection::addRange keys its range map on the contig string
#           VERBATIM (chopChr() is commented out in base/RangeList.h) and the
#           range is then served to tabix. A set file saying "chr22:..."
#           against a VCF whose records say "22" does not error: every range
#           misses, every gene comes back with zero variants, and rvtest exits
#           0. So the string is taken from the VCF's own ##contig vocabulary and
#           first data record, and every chromosome in the map must be
#           reachable there before a single line is written.
#
#           WHY ONE 1-BP RANGE PER VARIANT AND NOT ONE SPAN PER GENE. A gene
#           span pulls in every variant between pos_min and pos_max, including
#           the ones the annotation filter deliberately excluded (MODIFIER,
#           wrong impact class, another gene's exons). 1-bp ranges are the only
#           way the range interface can express "these variants and no others".
#           It is still not exact: a tabix interval matches a record over
#           [POS, POS+len(REF)-1], so a multi-base-REF record is matched a
#           second time by the range of a later variant inside its span.
#           scan.py predicts that over-count exactly; it is a property of the
#           range interface, not repairable here, and not silently hidden.
#
#           WHY THE COUNTS LIVE HERE. rvtest reports no allele count and no
#           carrier count for a set test. Both are computed from the same
#           variant_id -> set_name dict that just wrote the set file -- not by
#           re-opening the map, which would reintroduce exactly the channel this
#           process exists to remove. Cumulative MAC comes from two plink2
#           --freq counts runs (one per group); carriers come from a direct
#           decode of tested.bed. The two are CROSS-CHECKED PER VARIANT: the
#           bed-derived allele count must equal the .acount one, or the task
#           aborts. That equality is what proves the decode read the right
#           allele in the right sample order, and it is the reason carriers are
#           counted here rather than anywhere else.
#
#           A SAMPLE CARRYING TWO VARIANTS OF ONE GENE IS COUNTED ONCE -- and
#           listed. CMC collapses a sample to carrier / non-carrier, so a sample
#           with two minor alleles at two sites contributes exactly what a
#           single-site carrier does, while its MAC contribution is two. The two
#           counts therefore diverge on precisely these samples, and a reader
#           comparing "carriers" with "cumulative MAC" needs to be told where.
#           So per set and per sample the number of DISTINCT variant sites
#           carried is counted from the same decoded matrix; set_stats.tsv
#           carries the number of such samples per group and the maximum per
#           sample, and multi_carriers.tsv lists every one with its variants.
#
#           VARIANT IDS ARE chr<N>:<POS>:<REF>:<ALT>, STRICTLY. Two things parse
#           them: the minor-allele orientation here (bim A1 must equal the ID's
#           ALT) and scan.py's span-overlap predictor (len(REF)). Verified on
#           647,142 / 647,142 tested IDs; asserted on every row regardless.
# Component: assoc_gene
# ---------------------------------------------------------------------------
import argparse
import gzip
import json
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

MAP_COLUMNS = ('variant_id', 'chrom', 'pos', 'ref', 'alt', 'impact', 'effect',
               'biotype', 'gene_symbol', 'gene_id', 'n_genes', 'set_name')
INDEX_COLUMNS = ('set_name', 'gene_symbol', 'chrom', 'n_var_map', 'pos_min', 'pos_max',
                 'span_bp', 'n_high', 'n_moderate', 'n_low')
QC_COLUMNS = ('cohort', 'stratum', 'chrom', 'n_sets', 'n_variants', 'n_setfile_ranges',
              'vcf_chrom_string', 'bim_chrom_string', 'n_carrier_checks_passed',
              'all_assertions_passed')

# Per gene: the numbers the publication table is built from. Keyed by cohort and
# stratum so a stageAs rename downstream can never erase which cell a row is.
SET_STATS_COLUMNS = ('cohort', 'stratum', 'set_name', 'n_var', 'mac', 'mac_case',
                     'mac_control', 'n_carrier_case', 'n_carrier_control', 'n_case',
                     'n_control', 'n_multi_case', 'n_multi_control',
                     'max_var_per_sample_case', 'max_var_per_sample_control',
                     'n_flipped', 'min_obs_ct')

# One row per (set, sample) for every sample carrying a minor allele at TWO OR
# MORE distinct sites of the set. genotypes is per listed variant: 1 = het,
# 2 = hom for the minor allele. Sample IDs are inside; this file is for the
# analyst, not the manuscript.
MULTI_COLUMNS = ('cohort', 'stratum', 'set_name', 'fid', 'iid', 'group',
                 'n_variants_carried', 'n_hom', 'variant_ids', 'genotypes')

# Per variant: everything a per-gene detail figure or a variant-level supplement
# could want, with the group-wise genotype counts the carriers were summed from.
VARIANT_STATS_COLUMNS = ('cohort', 'stratum', 'set_name', 'variant_id', 'chrom', 'pos',
                         'ref', 'alt', 'impact', 'effect', 'minor_allele', 'flipped',
                         'ac_case', 'obs_case', 'ac_control', 'obs_control',
                         'mac_case', 'mac_control',
                         'n_het_case', 'n_hom_case', 'n_carrier_case', 'n_missing_case',
                         'n_het_control', 'n_hom_control', 'n_carrier_control',
                         'n_missing_control')

VARIANT_ID_RE = re.compile(r'^chr(\d+|X|Y|MT?):(\d+):([ACGT]+):([ACGT]+)$')
IMPACTS = ('HIGH', 'MODERATE', 'LOW')

# rvtest's loadRangeFile tokenises the set line with readLineBySep(&fd, "\t "),
# so a space inside a set name silently turns the name into two columns and the
# ranges into a third that is never read. ':' and ',' are the delimiters of
# rvtest's own range grammar, so a name carrying either makes the line ambiguous
# to anything that re-reads it.
FORBIDDEN_IN_SET_NAME = (',', ':')

# PLINK .bed magic, then SNP-major. Two bits per sample, low bits first:
#   0b00 homozygous A1 · 0b01 missing · 0b10 heterozygous · 0b11 homozygous A2
BED_MAGIC = b'\x6c\x1b\x01'
CODE_TO_A1_COPIES = np.array([2, -1, 1, 0], dtype=np.int8)

MAX_REPORT = 5


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--map', required=True,
                   help='map.<stratum>.tsv.gz from merge_gene_map (carries set_name)')
    p.add_argument('--gene-index', required=True, help='gene_index.<stratum>.tsv')
    p.add_argument('--bed', required=True, help='tested.bed, as plink2 --make-bed wrote it')
    p.add_argument('--bim', required=True, help='tested.bim, as plink2 --make-bed wrote it')
    p.add_argument('--fam', required=True, help='tested.fam; fixes the sample order')
    p.add_argument('--vcf', required=True,
                   help='tested.vcf.gz; read ONLY for the chromosome vocabulary')
    p.add_argument('--acount-case', required=True,
                   help='plink2 --freq counts .acount over the CASES only')
    p.add_argument('--acount-control', required=True,
                   help='plink2 --freq counts .acount over the CONTROLS only')
    p.add_argument('--case-keep', required=True, help='FID IID of the cases (PREP_PHENO_COV)')
    p.add_argument('--control-keep', required=True, help='FID IID of the controls')
    p.add_argument('--coding', required=True,
                   help='pheno_coding.json; n_case / n_control are asserted against the '
                        'keep lists and the .fam')
    p.add_argument('--cohort', required=True)
    p.add_argument('--stratum', required=True)
    p.add_argument('--out-dir', default='.')
    p.add_argument('--out-qc', default='test_inputs_qc.tsv')
    p.add_argument('--out-set-stats', default='set_stats.tsv')
    p.add_argument('--out-variant-stats', default='variant_stats.tsv')
    p.add_argument('--out-multi-carriers', default='multi_carriers.tsv')
    p.add_argument('--collision-token', default='@',
                   help='must match the merge_gene_map run that named these sets')
    p.add_argument('--bcftools', default=None,
                   help='optional bcftools binary; cross-checks the first VCF record')
    return p.parse_args()


def opener(path):
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path)


def chrom_key(c):
    """Genomic order with the non-numeric contigs (X, Y, MT) after the autosomes."""
    return (0, int(c), '') if c.isdigit() else (1, 0, c)


def sample(items):
    """A short, deterministic excerpt for an abort message."""
    s = sorted(str(x) for x in items)
    more = f' ... (+{len(s) - MAX_REPORT} more)' if len(s) > MAX_REPORT else ''
    return ', '.join(s[:MAX_REPORT]) + more


def parse_variant_id(vid, where):
    m = VARIANT_ID_RE.match(vid)
    if not m:
        raise SystemExit(f'ABORT: {where}: variant ID {vid!r} is not chr<N>:<POS>:<REF>:<ALT> '
                         f'with REF/ALT in ACGT. Both the minor-allele orientation and the '
                         f'span-overlap predictor parse this form, so it is asserted.')
    chrom, pos, ref, alt = m.groups()
    if ref == alt:
        raise SystemExit(f'ABORT: {where}: variant ID {vid!r} has REF == ALT')
    return chrom, int(pos), ref, alt


def read_map(path):
    """set_name -> {'chrom', 'rows': [(pos, variant_id, ref, alt, impact, effect)]},
    plus variant -> (chrom, pos)."""
    sets, var_pos = {}, {}
    n_rows = 0
    with opener(path) as fh:
        head = fh.readline().rstrip('\n').split('\t')
        col = {n: i for i, n in enumerate(head)}
        missing = [c for c in MAP_COLUMNS if c not in col]
        if missing:
            raise SystemExit(f'ABORT: {path} has no {missing} column(s); found {head}. '
                             f'A map without set_name is a per-chromosome map, not the '
                             f'merged stratum map this step consumes.')
        i_vid, i_chr, i_pos = col['variant_id'], col['chrom'], col['pos']
        i_ref, i_alt = col['ref'], col['alt']
        i_imp, i_eff, i_set = col['impact'], col['effect'], col['set_name']
        for ln, raw in enumerate(fh, 2):
            line = raw.rstrip('\n')
            if not line:
                continue
            f = line.split('\t')
            if len(f) != len(head):
                raise SystemExit(f'ABORT: {path}:{ln} has {len(f)} fields, expected {len(head)}')
            vid, chrom, imp, sn = f[i_vid], f[i_chr], f[i_imp], f[i_set]
            if chrom.startswith('chr'):
                raise SystemExit(f'ABORT: {path}:{ln} chrom is {chrom!r}; the map carries the '
                                 f'PLAIN number, matching bim column 1')
            if imp not in IMPACTS:
                raise SystemExit(f'ABORT: {path}:{ln} impact {imp!r} is not one of {IMPACTS}')
            id_chrom, id_pos, id_ref, id_alt = parse_variant_id(vid, f'{path}:{ln}')
            pos = int(f[i_pos])
            if (id_chrom, id_pos) != (chrom, pos) or (id_ref, id_alt) != (f[i_ref], f[i_alt]):
                raise SystemExit(f'ABORT: {path}:{ln} variant ID {vid!r} disagrees with its own '
                                 f'row ({chrom}:{pos} {f[i_ref]}>{f[i_alt]})')
            n_rows += 1
            e = sets.get(sn)
            if e is None:
                sets[sn] = e = {'chrom': chrom, 'rows': []}
            elif e['chrom'] != chrom:
                raise SystemExit(f'ABORT: {path}:{ln} set {sn!r} appears on chromosome '
                                 f'{chrom} and {e["chrom"]}. Re-run merge_gene_map: a set '
                                 f'name must name exactly one chromosome.')
            e['rows'].append((pos, vid, f[i_ref], f[i_alt], imp, f[i_eff]))
            seen = var_pos.get(vid)
            if seen is None:
                var_pos[vid] = (chrom, pos)
            elif seen != (chrom, pos):
                raise SystemExit(f'ABORT: {path}:{ln} variant {vid} is at {chrom}:{pos} here '
                                 f'and at {seen[0]}:{seen[1]} earlier')
    if not sets:
        raise SystemExit(f'ABORT: {path} has no data row; there is nothing to test')
    return sets, var_pos, n_rows


def read_gene_index(path):
    """[(set_name, chrom, n_var_map)] in file order -- the genomic order the files inherit."""
    rows = []
    with opener(path) as fh:
        head = fh.readline().rstrip('\n').split('\t')
        col = {n: i for i, n in enumerate(head)}
        missing = [c for c in INDEX_COLUMNS if c not in col]
        if missing:
            raise SystemExit(f'ABORT: {path} has no {missing} column(s); found {head}')
        for ln, raw in enumerate(fh, 2):
            line = raw.rstrip('\n')
            if not line:
                continue
            f = line.split('\t')
            if len(f) != len(head):
                raise SystemExit(f'ABORT: {path}:{ln} has {len(f)} fields, expected {len(head)}')
            rows.append((f[col['set_name']], f[col['chrom']], int(f[col['n_var_map']])))
    if not rows:
        raise SystemExit(f'ABORT: {path} lists no set; the denominator would be zero')
    return rows


def read_bim(path):
    """variant_id -> (row, chrom, pos, a1, a2), from the bim plink2 WROTE."""
    out = {}
    with open(path) as fh:
        for ln, line in enumerate(fh, 1):
            f = line.split()
            if len(f) < 6:
                raise SystemExit(f'ABORT: {path}:{ln} has {len(f)} fields, expected 6')
            vid = f[1]
            if vid in out:
                raise SystemExit(f'ABORT: {path}:{ln} repeats variant ID {vid!r}; rvtest would '
                                 f'test the position twice and say nothing')
            out[vid] = (ln - 1, f[0], int(f[3]), f[4], f[5])
    if not out:
        raise SystemExit(f'ABORT: {path} is empty; plink2 --extract wrote no variant')
    return out


def read_fam_ids(path):
    ids = []
    with open(path) as fh:
        for ln, line in enumerate(fh, 1):
            f = line.split()
            if len(f) < 2:
                raise SystemExit(f'ABORT: {path}:{ln} has fewer than 2 fields')
            ids.append((f[0], f[1]))
    if not ids:
        raise SystemExit(f'ABORT: {path} is empty')
    if len(set(ids)) != len(ids):
        raise SystemExit(f'ABORT: {path} repeats a (FID, IID)')
    return ids


def read_keep(path):
    out = set()
    with open(path) as fh:
        for ln, line in enumerate(fh, 1):
            f = line.split()
            if not f:
                continue
            if len(f) < 2:
                raise SystemExit(f'ABORT: {path}:{ln} has fewer than 2 fields; expected FID IID')
            out.add((f[0], f[1]))
    if not out:
        raise SystemExit(f'ABORT: {path} lists no sample')
    return out


def read_acount(path, label):
    """variant_id -> (alt_cts, obs_ct) from a plink2 --freq counts .acount file.

    Columns are read BY NAME. The real header for bfile input is
        #CHROM ID REF ALT PROVISIONAL_REF? ALT_CTS OBS_CT
    -- the PROVISIONAL_REF? column IS present, so positional reading is off by
    one and silently returns the wrong quantity. OBS_CT is the DENOMINATOR
    (observed alleles = 2 * non-missing samples), not a sample count.
    """
    out = {}
    with open(path) as fh:
        head = fh.readline().rstrip('\n').lstrip('#').split('\t')
        col = {n: i for i, n in enumerate(head)}
        for c in ('ID', 'ALT_CTS', 'OBS_CT'):
            if c not in col:
                raise SystemExit(f'ABORT: {path} has no {c} column; found {head}. This is '
                                 f'not a plink2 --freq counts output for the {label} group.')
        i_id, i_ac, i_obs = col['ID'], col['ALT_CTS'], col['OBS_CT']
        for ln, raw in enumerate(fh, 2):
            line = raw.rstrip('\n')
            if not line:
                continue
            f = line.split('\t')
            if len(f) != len(head):
                raise SystemExit(f'ABORT: {path}:{ln} has {len(f)} field(s), expected '
                                 f'{len(head)}')
            try:
                ac, obs = float(f[i_ac]), float(f[i_obs])
            except ValueError:
                raise SystemExit(f'ABORT: {path}:{ln} ALT_CTS/OBS_CT are '
                                 f'{f[i_ac]!r}/{f[i_obs]!r}, not numbers')
            if f[i_id] in out:
                raise SystemExit(f'ABORT: {path} repeats variant {f[i_id]!r}')
            out[f[i_id]] = (ac, obs)
    if not out:
        raise SystemExit(f'ABORT: {path} carries no data row for the {label} group')
    return out


def read_vcf_chroms(path, bcftools=None):
    """(contig vocabulary, first record's CHROM) -- the strings rvtest must be handed."""
    contigs, first = [], None
    with gzip.open(path, 'rt') as fh:
        for line in fh:
            if line.startswith('##contig=<'):
                for kv in line[len('##contig=<'):].rstrip().rstrip('>').split(','):
                    k, _s, v = kv.partition('=')
                    if k == 'ID':
                        contigs.append(v)
                        break
                continue
            if line.startswith('#'):
                continue
            first = line.split('\t', 1)[0]
            break
    if first is None:
        raise SystemExit(f'ABORT: {path} carries no data record; rvtest would report every '
                         f'gene as monomorphic and still exit 0')
    if bcftools:
        proc = subprocess.Popen([bcftools, 'view', '-H', str(path)],
                                stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
        try:
            line = proc.stdout.readline()
        finally:
            proc.stdout.close()
            proc.terminate()
            proc.wait()
        got = line.split('\t', 1)[0] if line else ''
        if got != first:
            raise SystemExit(f'ABORT: {bcftools} reads the first CHROM as {got!r}, this '
                             f'script reads {first!r}')
    return contigs, first


def resolve_chrom_strings(map_chroms, contigs, first_record, bim_chroms):
    """Map plain chromosome -> the contig string the VCF uses, or abort."""
    vocab = list(dict.fromkeys(list(contigs) + [first_record]))
    plain_of = {c: (c[3:] if c.startswith('chr') else c) for c in vocab}
    by_plain = defaultdict(list)
    for c, p in plain_of.items():
        by_plain[p].append(c)
    first_plain = plain_of[first_record]
    prefix = first_record[:len(first_record) - len(first_plain)]

    out = {}
    for chrom in sorted(map_chroms, key=chrom_key):
        cands = by_plain.get(chrom, [])
        if not cands:
            raise SystemExit(
                f'ABORT: chromosome {chrom} is in the map but not reachable in the VCF '
                f'(vocabulary: {sample(vocab)}). rvtest uses the contig string verbatim, '
                f'so every range on this chromosome would silently miss.')
        want = prefix + chrom
        if want not in cands:
            raise SystemExit(f'ABORT: the VCF names chromosome {chrom} as {sample(cands)} but '
                             f'its first record uses the {prefix + first_plain!r} style; the '
                             f'contig naming is not uniform')
        if len(cands) > 1:
            raise SystemExit(f'ABORT: chromosome {chrom} is ambiguous in the VCF: '
                             f'{sample(cands)} all resolve to it')
        out[chrom] = want
    if prefix:
        print(f'[emit_test_inputs] NOTE: the VCF prefixes contigs with {prefix!r} '
              f'(first record {first_record!r}) while the bim uses {sample(bim_chroms)}; '
              f'the set files follow the VCF', file=sys.stderr)
    return out


def check_set_name(name, chrom, token):
    if not name:
        raise SystemExit('ABORT: a set carries an empty name; rvtest skips the line and the '
                         'gene silently leaves the scan')
    if any(ch.isspace() for ch in name):
        raise SystemExit(f'ABORT: set name {name!r} contains whitespace. rvtest tokenises the '
                         f'set line on tab OR space, so the name would become two columns and '
                         f'the ranges would never be read.')
    bad = [ch for ch in FORBIDDEN_IN_SET_NAME if ch in name]
    if bad:
        raise SystemExit(f'ABORT: set name {name!r} contains {bad}, which are rvtest\'s own '
                         f'range delimiters')
    n_tok = name.count(token)
    if n_tok > 1:
        raise SystemExit(f'ABORT: set name {name!r} contains the collision token {token!r} '
                         f'{n_tok} times; it could no longer be split back into '
                         f'(symbol, chromosome)')
    if n_tok == 1:
        suffix = name.rsplit(token, 1)[1]
        if suffix != f'chr{chrom}':
            raise SystemExit(f'ABORT: set name {name!r} is suffixed {suffix!r} but its rows '
                             f'are on chromosome {chrom}. The name would mean one gene here '
                             f'and another in the next stratum.')


def minor_orientation(per_chrom, case_ac, ctrl_ac):
    """variant_id -> flipped (True when ALT is the MAJOR allele on the pooled counts).

    THE FLIP IS TAKEN ONCE, ON THE POOLED COUNTS, AND APPLIED TO BOTH GROUPS.
    rvtest's CMCWaldTest::fit collapses getFlippedToMinorPolymorphicGenotype(),
    so the burden counts MINOR alleles, and rvtest decides which allele that is
    from the pooled analysed sample set -- not per group. Copying that here has
    two consequences worth stating: mac == mac_case + mac_control holds by
    construction, and the two groups can never disagree about which allele is
    minor, which a per-group decision would permit.

    The two flips are also taken on the SAME SAMPLES. RVTEST's post-flight gate
    asserts N_INFORMATIVE == n_samples, so rvtest's pooled set is every sample;
    --nonfounders plus the case+control == tested.fam check in the process makes
    plink2's the same one.

    In practice the flip is a no-op at minAC >= 2 with MAF far below 0.5, but it
    is COMPUTED rather than assumed: one ALT-major variant would contribute
    2N - AC to its gene and swamp the set. n_flipped is the audit trail.
    """
    flip = {}
    for chrom in per_chrom:
        for sn, vrows in per_chrom[chrom]:
            for _p, vid, _r, _a, _i, _e in vrows:
                for src, label in ((case_ac, 'case'), (ctrl_ac, 'control')):
                    if vid not in src:
                        raise SystemExit(
                            f'ABORT: variant {vid!r} of set {sn!r} is absent from the '
                            f'{label} .acount. plink2 counted a different variant set '
                            f'than the one rvtest is handed.')
                ac_c, obs_c = case_ac[vid]
                ac_k, obs_k = ctrl_ac[vid]
                obs = obs_c + obs_k
                if obs <= 0:
                    raise SystemExit(f'ABORT: variant {vid!r} of set {sn!r} has OBS_CT 0 in '
                                     f'both groups; it is in no sample and cannot be tested')
                flip[vid] = (ac_c + ac_k) * 2.0 > obs
    return flip


def decode_bed(bed_path, n_samples, n_variants):
    """The packed .bed as an (M, B) uint8 array, magic asserted, SNP-major asserted."""
    with open(bed_path, 'rb') as fh:
        magic = fh.read(3)
    if magic != BED_MAGIC:
        raise SystemExit(f'ABORT: {bed_path} magic is {magic!r}, expected {BED_MAGIC!r} '
                         f'(PLINK 1 binary, SNP-major). Refusing to guess the layout.')
    bytes_per_variant = (n_samples + 3) // 4
    raw = np.fromfile(bed_path, dtype=np.uint8, offset=3)
    expected = n_variants * bytes_per_variant
    if raw.size != expected:
        raise SystemExit(f'ABORT: {bed_path} carries {raw.size:,} genotype byte(s); '
                         f'{n_variants:,} variants x {n_samples:,} samples needs {expected:,}')
    return raw.reshape(n_variants, bytes_per_variant), bytes_per_variant


def a1_copies(packed_rows, n_samples):
    """(m, N) int8 copies of A1 per sample, -1 for missing, from packed bed rows."""
    m = packed_rows.shape[0]
    codes = np.empty((m, packed_rows.shape[1] * 4), dtype=np.uint8)
    for k in range(4):
        codes[:, k::4] = (packed_rows >> (2 * k)) & 3
    return CODE_TO_A1_COPIES[codes[:, :n_samples]]


def group_counts(copies_minor, mask):
    """(n_het, n_hom, n_carrier, n_missing, mac) for one group, per variant."""
    sub = copies_minor[:, mask]
    n_missing = (sub < 0).sum(axis=1)
    n_het = (sub == 1).sum(axis=1)
    n_hom = (sub == 2).sum(axis=1)
    return n_het, n_hom, n_het + n_hom, n_missing, n_het + 2 * n_hom


def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    token = args.collision_token
    if not token:
        raise SystemExit('ABORT: --collision-token must be a non-empty string')

    sets, var_pos, n_map_rows = read_map(args.map)
    index_rows = read_gene_index(args.gene_index)
    bim = read_bim(args.bim)
    fam = read_fam_ids(args.fam)
    case_ids, ctrl_ids = read_keep(args.case_keep), read_keep(args.control_keep)
    with open(args.coding) as fh:
        coding = json.load(fh)

    # ---- assertion 1: rvtest is handed the variants plink2 WROTE, and they are
    # the mapped ones -- not the ones we asked plink2 to write.
    map_ids, bim_ids = set(var_pos), set(bim)
    if map_ids != bim_ids:
        missing, extra = map_ids - bim_ids, bim_ids - map_ids
        raise SystemExit(
            f'ABORT: the map and {args.bim} do not describe the same variant set. '
            f'{len(missing):,} map variant(s) absent from the bim [{sample(missing)}], '
            f'{len(extra):,} bim variant(s) absent from the map [{sample(extra)}]. '
            f'plink2 --extract dropped or added variants: rebuild the bim from '
            f'variants.{args.stratum}.txt before emitting anything.')
    for vid, (chrom, pos) in var_pos.items():
        _row, b_chrom, b_pos, a1, a2 = bim[vid]
        if (b_chrom, b_pos) != (chrom, pos):
            raise SystemExit(f'ABORT: {vid} is {chrom}:{pos} in the map and '
                             f'{b_chrom}:{b_pos} in the bim')
        _c, _p, ref, alt = parse_variant_id(vid, args.bim)
        # The decode below maps bim code 0b00 to "two copies of A1". That is only
        # "two copies of ALT" if A1 IS the ALT allele; asserted per variant, never
        # assumed, and the .acount cross-check below would catch a slip anyway.
        if (a1, a2) != (alt, ref):
            raise SystemExit(f'ABORT: {args.bim} carries {vid} as A1={a1} A2={a2}, but the ID '
                             f'says ALT={alt} REF={ref}. plink2 did not write the bim '
                             f'ALT-first; the bed decode would count the wrong allele.')

    # ---- assertion 2: the sample partition is the one the phenotype recorded.
    fam_set = set(fam)
    if not case_ids <= fam_set or not ctrl_ids <= fam_set:
        raise SystemExit(f'ABORT: the keep lists name samples that are not in {args.fam}')
    if case_ids & ctrl_ids:
        raise SystemExit(f'ABORT: {len(case_ids & ctrl_ids)} sample(s) are in BOTH keep lists')
    if len(case_ids) + len(ctrl_ids) != len(fam):
        raise SystemExit(f'ABORT: {len(case_ids)} case + {len(ctrl_ids)} control != '
                         f'{len(fam)} samples in {args.fam}; the groups must partition it')
    if (len(case_ids), len(ctrl_ids)) != (coding['n_case'], coding['n_control']):
        raise SystemExit(f'ABORT: keep lists hold {len(case_ids)}/{len(ctrl_ids)} but '
                         f'{args.coding} says {coding["n_case"]}/{coding["n_control"]}')
    n_case, n_ctrl = len(case_ids), len(ctrl_ids)
    case_mask = np.array([s in case_ids for s in fam], dtype=bool)
    ctrl_mask = ~case_mask

    # ---- assertion 3: the denominator and the emitted files name the same sets.
    index_names = {sn for sn, _c, _n in index_rows}
    if len(index_names) != len(index_rows):
        dupes = [sn for sn, _c, _n in index_rows if sum(1 for r in index_rows if r[0] == sn) > 1]
        raise SystemExit(f'ABORT: {args.gene_index} repeats set name(s) {sample(set(dupes))}')
    if index_names != set(sets):
        raise SystemExit(
            f'ABORT: {args.gene_index} and {args.map} disagree on the tested sets: '
            f'{len(index_names - set(sets)):,} in the index only '
            f'[{sample(index_names - set(sets))}], {len(set(sets) - index_names):,} in the '
            f'map only [{sample(set(sets) - index_names)}]. The Bonferroni denominator '
            f'counts the index, so a set in one and not the other is untested or uncounted.')

    per_chrom = defaultdict(list)
    for sn, chrom, n_var_map in index_rows:
        e = sets[sn]
        if e['chrom'] != chrom:
            raise SystemExit(f'ABORT: set {sn!r} is chromosome {chrom} in {args.gene_index} '
                             f'and {e["chrom"]} in {args.map}')
        check_set_name(sn, chrom, token)
        rows = sorted(e['rows'], key=lambda r: (r[0], r[1]))
        ids = [r[1] for r in rows]
        if len(set(ids)) != len(ids):
            seen, dup = set(), set()
            for v in ids:
                (dup if v in seen else seen).add(v)
            raise SystemExit(f'ABORT: set {sn!r} carries variant(s) {sample(dup)} more than '
                             f'once. build_gene_map writes one row per (variant, gene).')
        # ---- assertion 4: per gene, the set file says exactly what the map says.
        if len(rows) != n_var_map:
            raise SystemExit(f'ABORT: set {sn!r} would get {len(rows)} rvtest range(s) '
                             f'against n_var_map={n_var_map}')
        per_chrom[chrom].append((sn, rows))

    contigs, first_record = read_vcf_chroms(args.vcf, args.bcftools)
    bim_chroms = {v[1] for v in bim.values()}
    vcf_of = resolve_chrom_strings(set(per_chrom), contigs, first_record, bim_chroms)

    # ---- the minor-allele orientation, once, on the pooled counts.
    case_ac = read_acount(args.acount_case, 'case')
    ctrl_ac = read_acount(args.acount_control, 'control')
    flip = minor_orientation(per_chrom, case_ac, ctrl_ac)

    # ---- the bed, decoded per chromosome and cross-checked per variant.
    packed, _bpv = decode_bed(args.bed, len(fam), len(bim))

    # Everything that can abort on the inputs has now read; nothing has been
    # written. The writing loop below still aborts on the carrier cross-check,
    # and it does so BEFORE the set file for that chromosome is closed cleanly,
    # so a partial run cannot leave a plausible set file beside missing stats.
    qc, global_names = [], {}
    set_stats, variant_stats, multi_rows = [], [], []
    n_flip_total = n_checks = n_multi_total = 0
    for chrom in sorted(per_chrom, key=chrom_key):
        entries = per_chrom[chrom]
        vcf_chrom = vcf_of[chrom]
        bim_chrom = sorted({bim[r[1]][1] for _sn, rows in entries for r in rows})

        # Every variant of this chromosome's sets, decoded at once (fancy-indexed
        # from the packed bytes, so the bim need not be chromosome-contiguous).
        vids = [r[1] for _sn, rows in entries for r in rows]
        vids = list(dict.fromkeys(vids))
        row_of = {v: i for i, v in enumerate(vids)}
        alt = a1_copies(packed[[bim[v][0] for v in vids]], len(fam))
        # Minor copies: 2 - ALT copies where ALT is major; missing stays -1.
        flips = np.array([flip[v] for v in vids], dtype=bool)
        minor = np.where(alt < 0, -1, np.where(flips[:, None], 2 - alt, alt)).astype(np.int8)

        # Per-variant group counts, then THE CROSS-CHECK: the ALT count and the
        # observed-allele count the decode sees must be what plink2 --freq counts
        # wrote for the same group. If either differs for any variant, the decode
        # read the wrong allele or the wrong sample order, and nothing derived
        # from it can be trusted.
        stats = {}
        for label, mask, acount in (('case', case_mask, case_ac), ('control', ctrl_mask, ctrl_ac)):
            het, hom, carrier, miss, mac = group_counts(minor, mask)
            alt_sub = alt[:, mask]
            alt_ct = np.where(alt_sub < 0, 0, alt_sub).sum(axis=1)
            obs_ct = 2 * (alt_sub >= 0).sum(axis=1)
            for i, v in enumerate(vids):
                ac_exp, obs_exp = acount[v]
                if int(alt_ct[i]) != int(round(ac_exp)) or int(obs_ct[i]) != int(round(obs_exp)):
                    raise SystemExit(
                        f'ABORT: chr{chrom} {v} {label}: the bed decode counts '
                        f'ALT={int(alt_ct[i])} over OBS={int(obs_ct[i])} but plink2 --freq '
                        f'counts wrote ALT_CTS={ac_exp:g} OBS_CT={obs_exp:g}. The decode is '
                        f'reading the wrong allele or the wrong sample order; every carrier '
                        f'count would be wrong.')
                n_checks += 1
            stats[label] = (het, hom, carrier, miss, mac, alt_ct, obs_ct)

        n_ranges = 0
        set_path = out_dir / f'rvtest.set.chr{chrom}.txt'
        with open(set_path, 'w') as sf:
            for sn, rows in entries:
                # ---- assertion 5: unique across ALL chromosome files; rvtest's
                # OrderedMap merges a repeated name silently.
                if sn in global_names:
                    raise SystemExit(f'ABORT: set name {sn!r} is emitted on chromosome '
                                     f'{global_names[sn]} and {chrom}')
                global_names[sn] = chrom
                ranges = ','.join(f'{vcf_chrom}:{p}-{p}' for p, *_ in rows)
                sf.write(f'{sn}\t{ranges}\n')
                n_ranges += len(rows)

                idx = np.array([row_of[r[1]] for r in rows])
                sub = minor[idx]
                carriers_case = int((sub[:, case_mask] >= 1).any(axis=0).sum())
                carriers_ctrl = int((sub[:, ctrl_mask] >= 1).any(axis=0).sum())
                mac_case = int(stats['case'][4][idx].sum())
                mac_ctrl = int(stats['control'][4][idx].sum())
                n_flip = int(flips[idx].sum())
                n_flip_total += n_flip
                min_obs = int((stats['case'][6][idx] + stats['control'][6][idx]).min())

                # Distinct sites carried per sample. A sample at >= 2 is the case
                # where carriers and MAC diverge; it is counted per group and
                # listed with its variants so the divergence is inspectable.
                carried = sub >= 1
                n_sites = carried.sum(axis=0)
                multi = {}
                for label, mask in (('case', case_mask), ('control', ctrl_mask)):
                    js = np.flatnonzero(mask & (n_sites >= 2))
                    multi[label] = (len(js), int(n_sites[mask].max()) if mask.any() else 0)
                    for j in js:
                        hit = np.flatnonzero(carried[:, j])
                        multi_rows.append({
                            'cohort': args.cohort, 'stratum': args.stratum, 'set_name': sn,
                            'fid': fam[j][0], 'iid': fam[j][1], 'group': label,
                            'n_variants_carried': int(len(hit)),
                            'n_hom': int((sub[hit, j] == 2).sum()),
                            'variant_ids': ';'.join(rows[i][1] for i in hit),
                            'genotypes': ';'.join(str(int(sub[i, j])) for i in hit),
                        })
                n_multi_total += multi['case'][0] + multi['control'][0]
                # Invariants that make the numbers self-consistent. A carrier holds
                # at least one minor allele and at most two at each of the set's
                # sites, so carriers <= MAC <= 2 * n_var * carriers -- NOT
                # 2 * carriers: a sample carrying several sites of one gene counts
                # once here and once per site in the MAC (GTF3C3 in full_mainland:
                # 10 case carriers, 21 minor alleles). Nobody is a carrier who is
                # not in the group.
                for label, nc, mc, n_g in (('case', carriers_case, mac_case, n_case),
                                           ('control', carriers_ctrl, mac_ctrl, n_ctrl)):
                    if not (nc <= mc <= 2 * len(rows) * nc) or nc > n_g:
                        raise SystemExit(f'ABORT: set {sn!r} {label}: carriers={nc} mac={mc} '
                                         f'n_var={len(rows)} N={n_g} violate carriers <= mac '
                                         f'<= 2*n_var*carriers, carriers <= N')
                    if multi[label][0] > nc or multi[label][1] > len(rows):
                        raise SystemExit(f'ABORT: set {sn!r} {label}: {multi[label][0]} '
                                         f'multi-site carrier(s) among {nc} carrier(s), max '
                                         f'{multi[label][1]} site(s) of {len(rows)}')
                set_stats.append({
                    'cohort': args.cohort, 'stratum': args.stratum, 'set_name': sn,
                    'n_var': len(rows), 'mac': mac_case + mac_ctrl,
                    'mac_case': mac_case, 'mac_control': mac_ctrl,
                    'n_carrier_case': carriers_case, 'n_carrier_control': carriers_ctrl,
                    'n_case': n_case, 'n_control': n_ctrl,
                    'n_multi_case': multi['case'][0], 'n_multi_control': multi['control'][0],
                    'max_var_per_sample_case': multi['case'][1],
                    'max_var_per_sample_control': multi['control'][1],
                    'n_flipped': n_flip, 'min_obs_ct': min_obs,
                })
                for (pos, vid, ref, alt_allele, imp, eff), i in zip(rows, idx):
                    c, k = stats['case'], stats['control']
                    variant_stats.append({
                        'cohort': args.cohort, 'stratum': args.stratum, 'set_name': sn,
                        'variant_id': vid, 'chrom': chrom, 'pos': pos, 'ref': ref,
                        'alt': alt_allele, 'impact': imp, 'effect': eff,
                        'minor_allele': ref if flip[vid] else alt_allele,
                        'flipped': int(flip[vid]),
                        'ac_case': int(c[5][i]), 'obs_case': int(c[6][i]),
                        'ac_control': int(k[5][i]), 'obs_control': int(k[6][i]),
                        'mac_case': int(c[4][i]), 'mac_control': int(k[4][i]),
                        'n_het_case': int(c[0][i]), 'n_hom_case': int(c[1][i]),
                        'n_carrier_case': int(c[2][i]), 'n_missing_case': int(c[3][i]),
                        'n_het_control': int(k[0][i]), 'n_hom_control': int(k[1][i]),
                        'n_carrier_control': int(k[2][i]),
                        'n_missing_control': int(k[3][i]),
                    })
        n_distinct = len(vids)
        collapsed = n_ranges - len({(sn, r[0]) for sn, rows in entries for r in rows})
        qc.append({
            'cohort': args.cohort, 'stratum': args.stratum, 'chrom': chrom,
            'n_sets': len(entries), 'n_variants': n_distinct,
            'n_setfile_ranges': n_ranges,
            'vcf_chrom_string': vcf_chrom,
            'bim_chrom_string': bim_chrom[0] if len(bim_chrom) == 1 else 'MIXED',
            'n_carrier_checks_passed': 2 * n_distinct,
            'all_assertions_passed': 1,
        })
        print(f'[emit_test_inputs] chr{chrom}: {len(entries):,} set(s), '
              f'{n_distinct:,} variant(s), {n_ranges:,} range(s) as {vcf_chrom}:pos-pos '
              f'-> {set_path.name}'
              + (f' ({collapsed:,} range(s) share a position within their own set and '
                 f'consolidate inside rvtest)' if collapsed else ''))

    if len(set_stats) != len(global_names):
        raise SystemExit(f'ABORT: {len(set_stats)} stats row(s) for {len(global_names)} '
                         f'emitted set(s)')

    def write(path, columns, rows):
        with open(path, 'w') as fh:
            fh.write('\t'.join(columns) + '\n')
            for r in rows:
                fh.write('\t'.join(str(r[c]) for c in columns) + '\n')

    write(out_dir / args.out_set_stats, SET_STATS_COLUMNS, set_stats)
    write(out_dir / args.out_variant_stats, VARIANT_STATS_COLUMNS, variant_stats)
    write(out_dir / args.out_multi_carriers, MULTI_COLUMNS, multi_rows)
    write(out_dir / args.out_qc, QC_COLUMNS, qc)

    tot = sum(r['mac'] for r in set_stats)
    n_sets = sum(r['n_sets'] for r in qc)
    print(f'[emit_test_inputs] {args.cohort}/{args.stratum}: {n_sets:,} set(s) over '
          f'{len(qc)} chromosome(s) from {n_map_rows:,} map row(s); '
          f'{len(map_ids):,} variant(s) -> {out_dir / args.out_qc}')
    print(f'[emit_test_inputs] {len(set_stats):,} set(s) carrying {tot:,} minor allele(s) over '
          f'{n_case:,} case / {n_ctrl:,} control; {n_flip_total:,} variant(s) had ALT as the '
          f'MAJOR allele and were flipped (rvtest\'s pooled decision); '
          f'{n_checks:,} per-variant bed-vs-plink2 cross-checks passed '
          f'-> {out_dir / args.out_set_stats}, {out_dir / args.out_variant_stats}')
    n_multi_sets = sum(1 for r in set_stats if r['n_multi_case'] + r['n_multi_control'])
    print(f'[emit_test_inputs] NOTE: {n_multi_total:,} (set, sample) pair(s) in {n_multi_sets:,} '
          f'set(s) carry a minor allele at >= 2 distinct sites of the same set; CMC counts each '
          f'such sample ONCE while MAC counts every site -> {out_dir / args.out_multi_carriers}')


if __name__ == '__main__':
    try:
        main()
    except SystemExit:
        raise
    except Exception as e:
        print(f'Error in emit_test_inputs: {e}', file=sys.stderr)
        sys.exit(1)

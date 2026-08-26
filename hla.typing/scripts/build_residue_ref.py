#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : IMGT's curated protein alignments -> one (allele x IMGT position)
#           residue table per gene.
#
#           IPD-IMGT/HLA publishes these alignments. They are read, not recomputed.
#           The alternative — translating exon sequences and re-aligning them with
#           blastx + MUSCLE — costs three things at once:
#             * positions become alignment column indexes that restart per exon and
#               shift with every release, so "DRB1 position 13" is not reproducible
#               and not the position anyone else means by that name;
#             * blastx returning no hit yields no peptide, and the obvious handling
#               (write what the parser returned) puts the string "None" into the
#               reference as if it were a sequence;
#             * blastx SEG-masks low-complexity stretches, replacing real residues
#               with X in transmembrane regions.
#           None of those failure modes exist here because no alignment is computed.
#
#           Format of the source (A_prot.txt and friends):
#               # version: IPD-IMGT/HLA 3.64.0
#                Prot              -30                                1
#                                  |                                  |
#                A*01:01:01:01           MAVM APRTLLLLLS GALAL..TQT ...
#                A*02:07:01:01           ---- -----V---- -----..--- ...
#           The first allele listed is the reference; '-' means "same as reference",
#           '.' is an alignment gap, '*' is unknown, 'X' is a null-allele stop.
#           Numbering is IMGT's own: position 1 is the first residue of the mature
#           protein, and the leader peptide runs negative.
# Component: hla.typing
# ---------------------------------------------------------------------------
import argparse
import re
import sys
from pathlib import Path

ALLELE_RE = re.compile(r'^\s+([A-Z0-9]+\*[\d:]+[A-Z]?)\s')
PROT_RE = re.compile(r'^\s+Prot\s+(-?\d+)')
GROUP = 11          # IMGT writes 10 residues then one separator space


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--alignment-dir', required=True, help='directory of *_prot.txt')
    p.add_argument('--genes', required=True, help='comma-separated gene symbols')
    p.add_argument('--out-dir', default='residue_reference')
    p.add_argument('--out-summary', default='residue_reference_summary.tsv')
    return p.parse_args()


def parse_prot(path):
    """{allele: {imgt_position: residue}}, the reference allele, and the positions.

    Parsed by COLUMN, not by splitting on whitespace. An allele whose leader is
    shorter than the alignment's is padded with LEADING SPACES, and those spaces
    are numbered positions where that allele simply has no residue. Stripping them
    shifts every position by the width of the padding — six, for HLA-A, which is
    the difference between reading position 1 as G and reading it as Y.

    Three conventions, each verified against the file's own `|` anchors:
      * the sequence field starts at the column of the first `|`; within it every
        11th character is a group separator, not an alignment column;
      * a '.' in the REFERENCE row is an insertion carried by another allele. It
        occupies a column but has no IMGT position;
      * IMGT numbering has no zero: -1 is followed by 1.

    Only the FIRST block's `Prot` value is a starting point; later blocks continue
    from where the previous ended.
    """
    lines = Path(path).read_text(errors='replace').splitlines()

    field_start, start = None, None
    for i, raw in enumerate(lines):
        m = PROT_RE.match(raw)
        if m and start is None:
            start = int(m.group(1))
            # The marker line under the Prot header carries '|' at the anchors.
            if i + 1 < len(lines) and '|' in lines[i + 1]:
                field_start = lines[i + 1].index('|')
        if start is not None and field_start is not None:
            break
    if start is None or field_start is None:
        raise SystemExit(f'ABORT: {path} has no "Prot" header with a | anchor line')

    ref_allele, seqs = None, {}
    for raw in lines:
        if raw.startswith('#') or PROT_RE.match(raw):
            continue
        m = ALLELE_RE.match(raw)
        if not m:
            continue
        allele = m.group(1)
        field = raw[field_start:]
        # Drop the group separators by position; every other column is alignment.
        seq = ''.join(c for i, c in enumerate(field) if i % GROUP != GROUP - 1)
        if ref_allele is None:
            ref_allele = allele
        seqs.setdefault(allele, []).append(seq)

    if ref_allele is None:
        raise SystemExit(f'ABORT: {path} has no allele row')

    width = max(len(''.join(v)) for v in seqs.values())
    joined = {a: ''.join(v).ljust(width) for a, v in seqs.items()}
    ref = joined[ref_allele]

    col_pos, pos = [], start
    for ch in ref:
        if ch == '.':
            col_pos.append(None)
            continue
        col_pos.append(pos)
        pos += 1
        if pos == 0:
            pos = 1

    keep = [i for i, q in enumerate(col_pos) if q is not None]
    positions = [col_pos[i] for i in keep]

    out = {}
    for allele, seq in joined.items():
        out[allele] = {
            col_pos[i]: ('.' if seq[i] == ' ' else (ref[i] if seq[i] == '-' else seq[i]))
            for i in keep
        }
    return ref_allele, positions, out


def main():
    args = parse_args()
    adir = Path(args.alignment_dir)
    if not adir.is_dir():
        raise SystemExit(f'ABORT: no alignment directory at {adir}')

    genes = [g.strip() for g in args.genes.split(',') if g.strip()]
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    rows, absent = [], []
    for g in genes:
        src = adir / f'{g}_prot.txt'
        if not src.exists():
            # Recorded, not swallowed. A gene with no IMGT protein alignment simply
            # cannot have residues; saying so is the whole point.
            absent.append(g)
            continue
        ref_allele, positions, table = parse_prot(src)
        version = next((ln.split(':', 1)[1].strip()
                        for ln in src.read_text(errors='replace').splitlines()[:10]
                        if ln.startswith('# version')), 'unknown')

        with (out_dir / f'{g}.tsv').open('w') as fh:
            fh.write('allele\t' + '\t'.join(str(p) for p in positions) + '\n')
            for allele in sorted(table):
                fh.write(allele + '\t' + '\t'.join(table[allele][p] for p in positions) + '\n')

        rows.append(dict(gene=g, alleles=len(table), positions=len(positions),
                         first_position=positions[0], last_position=positions[-1],
                         reference_allele=ref_allele, imgt_version=version))
        print(f'[build_residue_ref] {g:6s} {len(table):6d} alleles x {len(positions):4d} '
              f'positions ({positions[0]}..{positions[-1]}), ref {ref_allele}')

    if not rows:
        raise SystemExit(f'ABORT: none of {genes} has a *_prot.txt in {adir}')

    with Path(args.out_summary).open('w') as fh:
        fh.write('\t'.join(rows[0]) + '\n')
        for r in rows:
            fh.write('\t'.join(str(v) for v in r.values()) + '\n')

    print(f'[build_residue_ref] {len(rows)} gene(s) -> {out_dir}/')
    if absent:
        print(f'    no IMGT protein alignment, so no residues possible: {", ".join(absent)}')


if __name__ == '__main__':
    try:
        main()
    except SystemExit:
        raise
    except Exception as e:
        print(f'Error in build_residue_ref: {e}', file=sys.stderr)
        sys.exit(1)

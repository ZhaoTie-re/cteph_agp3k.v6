"""Readers for HLA-HD 1.7.1 output, taken from the per-sample archives of step 00.

Formats, as written by HLA-HD 1.7.1 (checked on the archived output):

<s>_final.result.txt   one line per gene: `GENE<TAB>a1<TAB>a2[<TAB>a3<TAB>a4 ...]`
    a2 = '-'          the homozygous model won: one allele
    'Not typed'       no candidate allele (both columns)
    "Couldn't read result file."  on a line of its own, with NO gene name, when
                      pick_up_allele could not read an estimate. Lines are therefore
                      parsed by gene name, never by position.

<s>_<gene>.est.txt
    #Pair count<TAB>N          candidate pairs listed below
    #Best allele pair<TAB>M    the first M pairs share the best score (a tie when M > 1)
    side1<TAB>side2<TAB>scores one line per pair. Each side is a comma-separated list
                               of allele names that cannot be told apart; side2 is '-'
                               for the homozygous model; scores read
                               `exon2:12.7:comp.0,exon3:14.1:comp.0`
    #Other ambiguous pair<TAB>k   then k more pairs without scores
    No candidate.              nothing could be typed

<s>_<gene>.read.txt    one block per reported allele:
    HLA-X*..<TAB>n             reads assigned to this allele (shared reads included)
    R1 only / R2 only / Pair<TAB>k
    read lines: name, exon, start, end, [exon, start, end,] weight
"""
import io
import re
import tarfile
from dataclasses import dataclass, field

FAIL_TOKEN = "Couldn't read result file."
NOT_TYPED = 'Not typed'
SUFFIX_RE = re.compile(r'([A-Z])$')


def strip_hla(a):
    return a[4:] if a.startswith('HLA-') else a


def split_name(a):
    """'A*24:02:01:144Q' -> ('A', ['24','02','01','144'], 'Q')."""
    gene, rest = a.split('*', 1)
    suf = ''
    m = SUFFIX_RE.search(rest)
    if m:
        suf, rest = m.group(1), rest[:-1]
    return gene, rest.split(':'), suf


def common_prefix(names):
    """Deepest field depth shared by every name, as (name_at_depth, depth).

    The expression suffix (N, L, S, Q) belongs to a full allele name. It is carried
    to the shared prefix only when EVERY name carries the same suffix, because only
    then is the group as a whole non-expressed (e.g. all members of G*01:05 are N).
    """
    parts = [split_name(n) for n in names]
    gene = parts[0][0]
    depth = min(len(p[1]) for p in parts)
    d = 0
    while d < depth and len({tuple(p[1][:d + 1]) for p in parts}) == 1:
        d += 1
    if d == 0:
        return '', 0
    sufs = {p[2] for p in parts}
    suf = sufs.pop() if len(sufs) == 1 else ''
    return f"{gene}*{':'.join(parts[0][1][:d])}{suf}", d


def trim(name, depth):
    """'A*24:02:01:144Q' at depth 2 -> 'A*24:02' (suffix dropped: it is not a
    property of the 2-field group). Shallower names are returned unchanged."""
    gene, fields, suf = split_name(name)
    if len(fields) <= depth:
        return name
    return f"{gene}*{':'.join(fields[:depth])}"


@dataclass
class Est:
    n_pairs: int = 0
    n_best: int = 0
    best: list = field(default_factory=list)    # [(side1 list, side2 list or None)]
    other: list = field(default_factory=list)
    no_candidate: bool = False


@dataclass
class Reads:
    blocks: list = field(default_factory=list)  # [(allele, n_total, set(read names))]


def parse_final(text):
    genes, n_fail = {}, 0
    for ln in text.splitlines():
        if not ln.strip():
            continue
        if ln.strip() == FAIL_TOKEN:
            n_fail += 1
            continue
        f = ln.split('\t')
        genes[f[0]] = [x.strip() for x in f[1:] if x.strip()]
    return genes, n_fail


def _side(s):
    s = s.strip()
    if s == '-':
        return None
    return [strip_hla(x) for x in s.split(',') if x.strip()]


def parse_est(text):
    e = Est()
    target = e.best
    for ln in text.splitlines():
        if not ln.strip():
            continue
        if ln.startswith('No candidate'):
            e.no_candidate = True
            continue
        f = ln.split('\t')
        if f[0] == '#Pair count':
            e.n_pairs = int(f[1])
        elif f[0] == '#Best allele pair':
            e.n_best = int(f[1])
        elif f[0].startswith('#Other ambiguous pair'):
            target = e.other
        elif f[0].startswith('#'):
            continue
        else:
            target.append((_side(f[0]), _side(f[1]) if len(f) > 1 else None))
    return e


def parse_reads(text):
    r = Reads()
    cur = None
    for ln in text.splitlines():
        f = ln.split('\t')
        if f[0].startswith('HLA-') and len(f) == 2:
            cur = [strip_hla(f[0]), int(f[1]), set()]
            r.blocks.append(cur)
        elif f[0] in ('R1 only', 'R2 only', 'Pair'):
            continue
        elif cur is not None and len(f) >= 5:
            cur[2].add(f[0])
    r.blocks = [tuple(b) for b in r.blocks]
    return r


@dataclass
class SampleOutput:
    sample: str
    final: dict
    n_fail_lines: int
    est: dict
    reads: dict
    run_log: str
    validation: str
    task: dict


def read_archive(path, sample):
    """Everything step 01 needs from one archive."""
    est, reads, task = {}, {}, {}
    final, n_fail, run_log, validation = {}, 0, '', ''
    with tarfile.open(path, 'r:gz') as tf:
        for m in tf.getmembers():
            if not m.isfile():
                continue
            name = m.name.split('/', 1)[1]
            data = tf.extractfile(m).read().decode('utf-8', errors='replace')
            if name == f'result/{sample}_final.result.txt':
                final, n_fail = parse_final(data)
            elif name.endswith('.est.txt'):
                est[name[len(f'result/{sample}_'):-len('.est.txt')]] = parse_est(data)
            elif name.endswith('.read.txt'):
                reads[name[len(f'result/{sample}_'):-len('.read.txt')]] = parse_reads(data)
            elif name == 'hlahd.run.log':
                run_log = data
            elif name == 'typing_validation.tsv':
                validation = data
            elif name.startswith('task/'):
                task[name[5:]] = data
    return SampleOutput(sample, final, n_fail, est, reads, run_log, validation, task)

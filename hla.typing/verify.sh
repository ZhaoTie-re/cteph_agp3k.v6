#!/bin/bash
# Pre-flight and post-run checks for hla.typing. All read-only. Run from the
# component root.
#
# WHY THIS FILE EXISTS. Every documentation defect this component has had was a
# hand-typed literal that could not follow its data: a caption saying "(of 10)"
# when the denominator was 20-26, "105 samples" over a figure drawn from 61,424,
# "13-18 %" against a current 6 %. None of them were wrong when written. The rule
# they broke is that a number in prose must be derived, and sections 7-9 are that
# rule made mechanical.
#
# Sections 1-10 run on the source tree and must pass BEFORE a run. Sections 11-23
# read results/ and are skipped until it exists; after a run they must pass too.
P=/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6
C=$P/hla.typing
cd "$C"
fail=0
ck() { if [ "$2" -eq 0 ]; then printf "  ✓ %s\n" "$1"; else printf "  ✗ %s\n" "$1"; fail=$((fail+1)); fi; }
R=results
HAVE_RESULTS=0; [ -d "$R/05.qc" ] && HAVE_RESULTS=1

echo "═══ 1. syntax ═══"
source activate dsl2 2>/dev/null
out=$(nextflow lint hla.typing.nf nextflow.config 2>&1)
echo "$out" | grep -q 'no errors'; ck "nextflow lint: no errors" $?
source activate cteph_geno_pro 2>/dev/null
bad=0; for f in scripts/*.py; do python3 -m py_compile "$f" 2>/dev/null || bad=1; done
ck "every script py_compiles ($(ls scripts/*.py | wc -l))" $bad

echo "═══ 2. site rules ═══"
grep -qE '^\s*(cpus|memory)\s*=' hla.typing.nf nextflow.config; ck "no cpus/memory directive" $([ $? -eq 1 ] && echo 0 || echo 1)
grep -q 'task\.cpus' hla.typing.nf; ck "no task.cpus (ext.threads only)" $([ $? -eq 1 ] && echo 0 || echo 1)
bad=0
for m in $(grep -E 'clusterOptions' nextflow.config | grep -oE 'm=[0-9]+M' | grep -oE '[0-9]+' | sort -u); do
  [ $((m % 4571)) -eq 0 ] || { echo "    ✗ m=${m}M is not a multiple of 4571"; bad=1; }
done
ck "every clusterOptions m= is a multiple of 4571" $bad
# The check that would have caught PLOT_TYPING_CONFOUND, which ran without --rsc
# because it was missing from the withName selector while every sibling was in it.
python3 - <<'PY'
import re, sys, pathlib
nf = pathlib.Path('hla.typing.nf').read_text()
cfg = pathlib.Path('nextflow.config').read_text()
procs = set(re.findall(r'^process\s+(\w+)\s*\{', nf, re.M))
local = {p for p in procs if re.search(r'process\s+' + p + r'\s*\{[^}]*?executor\s+\'local\'', nf, re.S)}
covered = set()
for sel in re.findall(r"withName:\s*'([^']+)'", cfg):
    covered |= set(sel.split('|'))
missing = sorted((procs - local) - covered)
for m in missing:
    print(f'    ✗ {m} is a submitted process with no clusterOptions selector')
print(f'    {len(procs - local)} submitted process(es) checked, {len(local)} local')
sys.exit(1 if missing else 0)
PY
ck "every submitted process carries --rsc" $?

echo "═══ 2b. publishDir referencing a process input must be a closure ═══"
python3 - <<'PY'
import re, sys, pathlib
nf = pathlib.Path('hla.typing.nf').read_text()
bad = []
for m in re.finditer(r'^\s*publishDir\s+"([^"]*)"', nf, re.M):
    for var in re.findall(r'\$\{([^}]+)\}', m.group(1)):
        if not var.strip().startswith('params.'):
            bad.append((nf[:m.start()].count('\n') + 1, var.strip()))
for line, var in bad:
    print(f'    ✗ hla.typing.nf:{line}: publishDir without a closure references {var}')
sys.exit(1 if bad else 0)
PY
ck "publishDir closures" $?

echo "═══ 3. every flag a process passes is accepted by its script ═══"
python3 - <<'PY'
import re, subprocess, sys, pathlib
nf = pathlib.Path('hla.typing.nf').read_text()
bodies = dict(re.findall(r'process\s+(\w+)\s*\{(.*?)\n\}', nf, re.S))
# Explicit, so a process whose script is not listed is itself a failure.
PAIRS = {'BUILD_MANIFEST': 'build_manifest.py', 'CONTIG_LIST': 'contig_list.py',
         'EXTRACT_READS': 'extract_hla_reads.py', 'HLAHD': 'validate_typing.py',
         'COLLECT_ALLELES': 'collect_alleles.py', 'BUILD_RESIDUE_REF': 'build_residue_ref.py',
         'RESIDUE_MATRIX': 'residue_matrix.py', 'TYPING_QC': 'typing_qc.py',
         'ALLELE_FREQ_CHECK': 'allele_freq_check.py', 'PLOT_TYPING_QC': 'plot_typing_qc.py',
         'PLOT_ALLELE_FREQ': 'plot_allele_freq.py',
         'PLOT_TYPING_CONFOUND': 'plot_typing_confound.py'}
bad = 0
for proc, body in bodies.items():
    if re.search(r'python3 \$\{(script|validator)\}', body) and proc not in PAIRS:
        print(f'    ✗ {proc} runs a script that verify.sh does not know'); bad = 1
for proc, sc in sorted(PAIRS.items()):
    p = pathlib.Path('scripts', sc)
    if not p.exists():
        print(f'    ✗ {proc}: {sc} missing'); bad = 1; continue
    m = re.search(r'python3 \$\{(?:script|validator)\}(.*?)"""', bodies.get(proc, ''), re.S)
    if not m:
        print(f'    ✗ {proc}: no python3 line'); bad = 1; continue
    passed = set(re.findall(r'(?<![\w-])--([a-z][a-z0-9-]*)', m.group(1)))
    h = subprocess.run([sys.executable, str(p), '--help'], capture_output=True, text=True)
    if h.returncode != 0:
        print(f'    ✗ {proc}: {sc} --help fails'); bad = 1; continue
    known = set(re.findall(r'(?<![\w-])--([a-z][a-z0-9-]*)', h.stdout))
    for f in sorted(passed - known):
        print(f'    ✗ {proc}: {sc} does not accept --{f}'); bad = 1
    # The reverse: a required flag the script declares and the process never sends.
    req = set(re.findall(r'^\s+--([a-z][a-z0-9-]*)\s+[A-Z_]+\n', h.stdout, re.M))
    usage = h.stdout.split('\n\n')[0]
    for f in sorted(req):
        if f'--{f}' in usage and f'[--{f}' not in usage and f not in passed:
            print(f'    ✗ {proc}: {sc} requires --{f} and the process does not pass it')
            bad = 1
sys.exit(bad)
PY
ck "process flags match script --help" $?

echo "═══ 4. no dead flags ═══"
python3 - <<'PY'
import re, sys, pathlib
bad = 0
for p in sorted(pathlib.Path('scripts').glob('*.py')):
    src = p.read_text()
    declared = set(re.findall(r"add_argument\('--([a-z][a-z0-9-]*)'", src))
    for flag in sorted(declared):
        attr = flag.replace('-', '_')
        # argparse dest, or the flag re-read from sys.argv; either counts as used.
        if not re.search(rf'\bargs?\.{attr}\b|\ba\.{attr}\b', src):
            print(f'    ✗ {p.name}: --{flag} is declared and never read'); bad = 1
sys.exit(bad)
PY
ck "every declared flag is read" $?

echo "═══ 5. the shared toolkit is untouched ═══"
n=$(cd "$P" && git status --porcelain analysis/_shared/scripts 2>/dev/null | wc -l)
ck "analysis/_shared/scripts has no local modification (${n:-0})" $([ "${n:-0}" -eq 0 ] && echo 0 || echo 1)

echo "═══ 6. plot scripts ═══"
bad=0
for f in scripts/plot_*.py; do
  grep -q 'import plot_style' "$f" || { echo "    ✗ $f does not import plot_style"; bad=1; }
  grep -q 'import figure_doc' "$f" || { echo "    ✗ $f does not import figure_doc"; bad=1; }
  python3 "$f" --help >/dev/null 2>&1 || { echo "    ✗ $f --help fails"; bad=1; }
done
ck "every plot_*.py imports the shared style and answers --help" $bad
bad=0
for f in scripts/plot_*.py; do
  grep -q 'numbers=' "$f" || { echo "    ✗ $f: write_doc has no numbers="; bad=1; }
done
ck "every figure's sidecar carries its numbers" $bad

echo "═══ 7. retired claims ═══"
# Each string below was TRUE ONCE and is FALSE NOW. That is the test -- not
# whether a phrase is old, but whether it still asserts something the data
# support. "JPT panel" is deliberately NOT on this list: the 1000 Genomes JPT
# panel is a documented secondary reference that runs on every execution, so
# naming it is correct. Banning the name rather than the claim is how a checker
# starts forcing prose to lie by circumlocution.
bad=0
for pat in \
  'Five loci' \
  '13-18' \
  '13–18' \
  '(of 10)' \
  '105 samples' \
  'params.TruthPanel' \
  'params.TruthFormat' \
  'params.RefLoci ' \
  '--platform-qc' \
  '--restrict-ids' \
  '--truth-format 1kg_wide` to run' \
  ; do
  # NAMED FILES, not a recursive walk of the component root. `grep -r .` here
  # descends results/, whose 01.reads and 02.typing are symlink farms over work/
  # -- 3.2 M inodes on Lustre. Measured: it did not finish in 17 minutes, while
  # the list below takes under a second. --exclude-dir was not enough, because
  # the cost is the directory walk itself and not the matching.
  hits=$(grep -nF -- "$pat" \
           hla.typing.nf nextflow.config README.md \
           scripts/*.py docs/*.md \
           $([ -f report/hla_typing_report.qmd ] && echo report/hla_typing_report.qmd) \
           $(ls "$R"/figures/*.md 2>/dev/null) \
         2>/dev/null | grep -v '\[retired\]')
  # A line may quote a retired claim if it is explaining that the claim was
  # retired. It says so with the marker [retired] on the same line. That keeps
  # the escape hatch narrow and, more importantly, visible.
  [ -n "$hits" ] && { echo "$hits" | sed 's/^/    ✗ /'; bad=1; }
done
ck "no retired claim survives" $bad

echo "═══ 8. no numeric literal inside a caption string ═══"
python3 - <<'PY'
import ast, sys, pathlib, re
# A digit a human typed into a caption is a defect; a digit an f-string FIELD
# produces is the whole point. So this walks the string parts only:
#   - a plain string constant is prose and is checked
#   - a JoinedStr's literal segments are prose and are checked
#   - a FormattedValue's expression is derived and is NOT checked
#   - a FormattedValue's format_spec ('.2f', '.1%') is not prose at all
# The first version used ast.walk() and flagged every format spec in the
# component -- a checker that cries wolf until someone switches it off.
#
# Terms of art are not results. Note there is no trailing \b: '95 %' ends in a
# non-word character, so \b after it never matches and the whole allow-list
# silently did nothing.
OK = re.compile(r'[1-4]-field|[1-4] fields|95 %|3\.64\.0|1\.7\.1|DRB3/4/5|DRB[345]|'
                r'DPA1|DPB1|DQA1|DQB1|DRB1|61KJPN|1000 Genomes|0\.98|0\.02|0\.4|'
                r'17.21|§ ?\d|2×2|2x2')
CALLS = {'caption_block', 'write_doc'}


def is_prose(v):
    """A format spec is not prose.

    CPython puts a FormattedValue's format_spec constants ('.2f', '.1%') straight
    into the enclosing JoinedStr's .values on some versions, so skipping
    format_spec by structure is not enough. A format spec is short and has no
    space; a literal segment of prose that carries a number carries the words
    around it too, because that is what an f-string literal segment IS.
    """
    return len(v) > 4 or ' ' in v


def prose(node):
    """Every string a human typed inside this expression."""
    if isinstance(node, ast.Constant) and isinstance(node.value, str):
        if is_prose(node.value):
            yield node
    elif isinstance(node, ast.JoinedStr):
        for v in node.values:
            if isinstance(v, ast.Constant) and isinstance(v.value, str) \
                    and is_prose(v.value):
                yield v
    elif isinstance(node, ast.BinOp):
        yield from prose(node.left); yield from prose(node.right)
    elif isinstance(node, (ast.List, ast.Tuple, ast.Set)):
        for e in node.elts:
            yield from prose(e)
    elif isinstance(node, ast.Dict):
        for v in node.values:
            yield from prose(v)


bad = 0
for p in sorted(pathlib.Path('scripts').glob('plot_*.py')):
    for node in ast.walk(ast.parse(p.read_text())):
        if not isinstance(node, ast.Call):
            continue
        name = getattr(node.func, 'attr', None) or getattr(node.func, 'id', None)
        if name not in CALLS:
            continue
        for kw in node.keywords:
            for s_ in prose(kw.value):
                # A digit INSIDE a word is part of a name -- AGP3K, DNBSeq-T7,
                # 61KJPN -- not a quantity. Only a digit that starts a token can
                # be a number a human typed.
                if re.search(r'(?<![A-Za-z])\d', OK.sub('', s_.value)):
                    print(f'    ✗ {p.name}:{s_.lineno}: {name}({kw.arg}=) has a typed '
                          f'number: {s_.value.strip()[:60]!r}')
                    bad = 1
sys.exit(bad)
PY
ck "captions interpolate rather than state" $?

echo "═══ 9. every number in a document is derivable ═══"
if [ -f "$R/_run_info/facts.json" ]; then
python3 - <<'PY'
import json, re, sys, pathlib
facts = json.load(open('results/_run_info/facts.json'))
allowed = set()
for v in facts['rendered'].values():
    allowed |= set(re.findall(r'\d[\d,]*(?:\.\d+)?', str(v)))
    allowed.add(str(v).strip())
# A DOCUMENT MAY ROUND A DERIVED NUMBER TO ANY PRECISION; IT MAY NOT INVENT ONE.
# Forcing prose to one fixed precision is the wrong rule -- it makes a sentence
# read "r = 0.9930" where "r = 0.993" is what a person would write, and it
# tempts the author to widen the allow-list instead. So every raw fact is also
# allowed at every rounding, as a fraction and as a percentage.
for v in facts.get('raw', {}).values():
    if not isinstance(v, (int, float)) or isinstance(v, bool):
        continue
    for k in range(7):
        allowed.add(f'{v:.{k}f}'.rstrip('.') if k == 0 else f'{v:.{k}f}')
    for k in range(5):
        allowed.add(f'{100 * v:.{k}f}'.rstrip('.') if k == 0 else f'{100 * v:.{k}f}')
    allowed.add(f'{v:,.0f}')
p = pathlib.Path('docs/NUMBERS_ALLOWED.md')
if p.exists():
    allowed |= {ln.split('#')[0].strip() for ln in p.read_text().splitlines()
                if ln.split('#')[0].strip()}
bad = 0
for md in sorted(list(pathlib.Path('.').glob('*.md')) + list(pathlib.Path('docs').glob('*.md'))):
    fence = False
    for i, line in enumerate(md.read_text().splitlines(), 1):
        if line.lstrip().startswith('```'):
            fence = not fence
            continue
        # Fenced blocks are commands and paths, not claims. Indented blocks are
        # the same. A URL carries numbers that mean nothing here. TABLES ARE NOT
        # EXEMPT -- they are where most of a document's claims actually live, and
        # exempting them was a hole big enough to drive every past defect through.
        if fence or line.startswith('    ') or 'http' in line:
            continue
        # Inline code spans hold paths, filenames and flags -- identifiers, not
        # claims -- so they are removed before the line is scanned at all.
        line = re.sub(r'`[^`]*`', '', line)
        # And a bare number must not be followed by a letter or by .letter, which
        # is what `05.qc`, `00.manifest` and `3.64.0` look like once the backticks
        # are gone.
        # \d+(?:,\d{3})* -- a comma must be a THOUSANDS separator, not the comma
        # in "at 2, 3 or 4 fields", which the looser [\d,]* read as the number 2,.
        for n in re.findall(r'(?<![\w.:/-])\d+(?:,\d{3})*(?:\.\d+)?(?![\w])(?!\.\w)\s*%?',
                            line):
            n = n.strip()
            if n in allowed or n.rstrip(' %') in allowed or n.replace(' ', '') in allowed:
                continue
            print(f'    ✗ {md}:{i}: {n!r} is in no fact and in no allow-list')
            bad = 1
sys.exit(bad)
PY
ck "documents quote only derived numbers" $?
else
  echo "  – facts.json not built yet; skipped"
fi

echo "═══ 10. documentation is present ═══"
bad=0
for f in README.md docs/METHODS.md docs/OUTPUTS.md docs/STUDY_NOTES.md docs/OPEN_QUESTIONS.md \
         docs/NUMBERS_ALLOWED.md; do
  [ -f "$f" ] || { echo "    ✗ $f missing"; bad=1; }
done
ck "README + docs present" $bad
if [ -d "$R/figures" ]; then
  bad=0
  for png in "$R"/figures/*.png; do
    [ -f "${png%.png}.md" ] || { echo "    ✗ $(basename "$png") has no sidecar"; bad=1; }
  done
  ck "every figure has a sidecar .md" $bad
fi

if [ "$HAVE_RESULTS" -eq 0 ]; then
  echo; echo "results/ not present — sections 11-23 skipped."
  echo; [ $fail -eq 0 ] && echo "ALL STATIC CHECKS PASSED" || echo "$fail CHECK(S) FAILED"
  exit $((fail > 0))
fi

echo "═══ 11. the sample chain, re-derived from this component's inputs ═══"
python3 - <<'PY'
import sys, pandas as pd, pathlib
P = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6'
wb = pd.read_excel(f'{P}/info/cteph_agp3k.v6.20260507.xlsx', engine='openpyxl', dtype=str)
v6 = wb[wb.Flag_JHRPv6.astype(str).str.strip().str.lower().isin(('true', '1', 'yes'))]
qc = {l.split()[0] for l in open(f'{P}/wgs.auto.par/results/07_sample_qc/run_qc/'
                                f'cteph_agp3k_v6_wgs_merged.sample_qc.keep.id') if l.strip()}
ch = {l.split()[0] for l in open(f'{P}/PopGMM_output/full_mainland.fid_iid.txt') if l.strip()}
man = pd.read_csv('results/00.manifest/sample_manifest.tsv', sep='\t', dtype=str)
ids = set(man.sample_id)
n_qc = len(v6[v6.ID_JHRPv6.astype(str).isin(qc)])
print(f'    workbook {len(wb):,} -> v6 {len(v6):,} -> sample QC {n_qc:,} -> cohort {len(man):,}')
bad = 0
if ch - qc:
    print(f'    ✗ {len(ch - qc)} cohort sample(s) are not in the sample-QC keep list'); bad = 1
if ids != ch:
    print(f'    ✗ manifest != cohort list: {len(ch - ids)} missing, {len(ids - ch)} extra'); bad = 1
if len(man) != len(ch):
    print(f'    ✗ manifest has {len(man)} rows, cohort list has {len(ch)}'); bad = 1
sys.exit(bad)
PY
ck "workbook -> sample QC -> cohort reproduces the manifest" $?

echo "═══ 12. typing_status re-derived from allele_calls ═══"
python3 - <<'PY'
import sys, pandas as pd
c = pd.read_csv('results/03.alleles/allele_calls.tsv', sep='\t', dtype=str).fillna('')
s = pd.read_csv('results/03.alleles/typing_status.tsv', sep='\t', dtype={'sample_id': str})
state = [col for col in c.columns if col.endswith('_state')]
got = {k: int((c[state] == k).to_numpy().sum())
       for k in ('called', 'hemizygous', 'not_typed', 'failed')}
bad = 0
for k, v in got.items():
    if k in s.columns and int(s[k].sum()) != v:
        print(f'    ✗ {k}: status says {int(s[k].sum()):,}, calls say {v:,}'); bad = 1
print('    ' + ', '.join(f'{k} {v:,}' for k, v in got.items()))
sys.exit(bad)
PY
ck "typing_status reproduces allele_calls" $?

echo "═══ 13. every per-sample table is the cohort, exactly ═══"
python3 - <<'PY'
import sys, pandas as pd
P = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6'
ch = {l.split()[0] for l in open(f'{P}/PopGMM_output/full_mainland.fid_iid.txt') if l.strip()}
FILES = ['00.manifest/sample_manifest.tsv', '03.alleles/allele_calls.tsv',
         '03.alleles/typing_status.tsv', '04.residues/allele_dosage.tsv',
         '04.residues/residue_dosage.tsv', '04.residues/residue_diplotype.tsv',
         '05.qc/typing_qc_sample.tsv', '05.qc/allele_frequency_sample.tsv',
         '05.qc/allele_frequency_sample.1kg.tsv']
bad = 0
for f in FILES:
    d = pd.read_csv(f'results/{f}', sep='\t', usecols=[0], dtype=str)
    got = set(d.iloc[:, 0])
    if got != ch or len(d) != len(ch):
        print(f'    ✗ {f}: {len(d)} rows, {len(ch - got)} missing, {len(got - ch)} extra')
        bad = 1
print(f'    {len(FILES)} per-sample table(s), all {len(ch):,} rows and the cohort id set')
sys.exit(bad)
PY
ck "per-sample tables == the cohort" $?

echo "═══ 14. nothing failed ═══"
python3 - <<'PY'
import sys, pandas as pd
bad = 0
g = pd.read_csv('results/05.qc/typing_qc_gene.tsv', sep='\t')
if 'failed' in g.columns and int(g.failed.sum()):
    print(f'    ✗ {int(g.failed.sum())} gene-level failure(s)'); bad = 1
u = pd.read_csv('results/04.residues/allele_unmatched.tsv', sep='\t')
if len(u):
    print(f'    ✗ allele_unmatched.tsv has {len(u)} row(s)'); bad = 1
sys.exit(bad)
PY
ck "failed == 0 and allele_unmatched is header-only" $?

echo "═══ 15. allele_frequency_summary re-derived from allele_frequency_check ═══"
python3 - <<'PY'
import sys, numpy as np, pandas as pd
from scipy import stats
bad = 0
for suf in ('', '.1kg'):
    c = pd.read_csv(f'results/05.qc/allele_frequency_check{suf}.tsv', sep='\t')
    s = pd.read_csv(f'results/05.qc/allele_frequency_summary{suf}.tsv', sep='\t')
    for _, row in s.iterrows():
        d = c[c.gene == row.gene]
        x, y = d.freq_ref.astype(float), d.freq_ctrl.astype(float)
        r = float(np.corrcoef(x, y)[0, 1])
        rho = float(stats.spearmanr(x, y).statistic)
        md = float((y - x).abs().max())
        # Pearson and max|diff| must reproduce exactly. SPEARMAN CANNOT, and the
        # reason is not sloppiness: allele_freq_check.py writes freq_ctrl/freq_ref
        # rounded to 5 dp but computes rho from the unrounded vectors. rho is
        # rank-based, so rounding fuses rare alleles into ties that do not exist
        # upstream and moves it in the third decimal. Observed spread over 13 loci
        # is < 0.004; the tolerance is the rounding's, not a fudge.
        for name, got, want, tol in (('r', r, float(row.pearson_r), 5e-4),
                                     ('max|diff|', md, float(row.max_abs_diff), 5e-4),
                                     ('rho', rho, float(row.spearman_rho), 1e-2)):
            if abs(got - want) > tol:
                print(f'    ✗ {suf or "jmorp"} {row.gene} {name}: {got:.4f} vs {want:.4f} '
                      f'(tol {tol})')
                bad = 1
print('    both panels re-derived per locus (rho to the 5-dp rounding of the table)')
sys.exit(bad)
PY
ck "summary reproduces check" $?

echo "═══ 16. typing_confound re-derived from allele_frequency_sample ═══"
python3 - <<'PY'
import sys, pandas as pd
from scipy import stats
d = pd.read_csv('results/05.qc/allele_frequency_sample.tsv', sep='\t', dtype={'sample_id': str})
t = pd.read_csv('results/05.qc/typing_confound.tsv', sep='\t')
bad = 0
for stratum, sel in (('controls', d.group == 'AGP3K'), ('cases', d.group != 'AGP3K')):
    g = d[sel]
    want = float(t.loc[t.stratum == stratum, 'share_unconfirmed'].iloc[0])
    got = g.n_unconfirmed.sum() / g.n_chr.sum()
    if abs(got - want) > 1e-5:
        print(f'    ✗ {stratum}: {got:.5f} vs {want:.5f}'); bad = 1
    else:
        print(f'    {stratum}: {got:.2%}')
mw = t.mannwhitney_p_depth.dropna()
if len(mw):
    d['observed_depth'] = pd.to_numeric(d.observed_depth, errors='coerce')
    cc, ca = d[d.group == 'AGP3K'].observed_depth, d[d.group != 'AGP3K'].observed_depth
    got = stats.mannwhitneyu(ca.dropna(), cc.dropna())[1]
    if abs(got - float(str(mw.iloc[0]))) > 1e-3:
        print(f'    ✗ depth Mann-Whitney: {got:.3g} vs {mw.iloc[0]}'); bad = 1
sys.exit(bad)
PY
ck "confound table reproduces the per-sample table" $?

echo "═══ 17. both reference panels, each with its own locus list ═══"
python3 - <<'PY'
import sys, pandas as pd
bad = 0
sizes = {}
for tag, suf in (('jmorp', ''), ('1kg', '.1kg')):
    s = pd.read_csv(f'results/05.qc/allele_frequency_summary{suf}.tsv', sep='\t')
    sizes[tag] = (len(s), int(s.n_chr_ref.max()))
    print(f'    {tag}: {len(s)} loci, up to {int(s.n_chr_ref.max()):,} reference chromosomes')
if sizes['jmorp'][0] <= sizes['1kg'][0]:
    print('    ✗ the primary panel should carry more loci than the secondary'); bad = 1
if sizes['jmorp'][1] <= sizes['1kg'][1] * 10:
    print('    ✗ the two panels are implausibly close in size'); bad = 1
sys.exit(bad)
PY
ck "two panels present and distinguishable" $?

echo "═══ 18. the matrices were regenerated, not row-subset ═══"
python3 - <<'PY'
import sys, pandas as pd
bad = 0
for name, kind in (('allele_dosage', 'd'), ('residue_dosage', 'd'), ('residue_diplotype', 'p')):
    d = pd.read_csv(f'results/04.residues/{name}.tsv', sep='\t',
                    dtype={'sample_id': str}).set_index('sample_id')
    if kind == 'd':
        n = int((d.apply(pd.to_numeric, errors='coerce').fillna(0).sum(axis=0) == 0).sum())
        what = 'all-zero'
    else:
        n = int((d.nunique(dropna=True) <= 1).sum())
        what = 'single-valued'
    print(f'    {name}: {d.shape[0]:,} x {d.shape[1]:,}, {n} {what} column(s)')
    if n:
        print(f'    ✗ {name} carries {n} {what} column(s) — it was subset, not rebuilt')
        bad = 1
sys.exit(bad)
PY
ck "no dead column in any matrix" $?

echo "═══ 19. provenance: no script is newer than the results it made ═══"
python3 - <<'PY'
import sys, pathlib
# The check that would have caught residue_matrix.py being edited six days after
# results/ was written, which made -resume silently regenerate the association's
# inputs.
PAIRS = {'collect_alleles.py': 'results/03.alleles/allele_calls.tsv',
         'residue_matrix.py': 'results/04.residues/residue_dosage.tsv',
         'typing_qc.py': 'results/05.qc/typing_qc_sample.tsv',
         'allele_freq_check.py': 'results/05.qc/allele_frequency_sample.tsv',
         'build_manifest.py': 'results/00.manifest/sample_manifest.tsv',
         'plot_typing_qc.py': 'results/figures/typing_qc.png',
         'plot_allele_freq.py': 'results/figures/allele_frequency.png',
         'plot_typing_confound.py': 'results/figures/typing_confound.png'}
bad = 0
for sc, out in sorted(PAIRS.items()):
    s, o = pathlib.Path('scripts', sc), pathlib.Path(out)
    if not (s.exists() and o.exists()):
        continue
    if s.stat().st_mtime > o.stat().st_mtime + 1:
        print(f'    ✗ {sc} is newer than {out} — rerun it or the table is stale'); bad = 1
sys.exit(bad)
PY
ck "no output predates the script that made it" $?

echo "═══ 20. run_manifest matches the .nf ═══"
python3 - <<'PY'
import json, re, sys, pathlib
nf = pathlib.Path('hla.typing.nf').read_text()
m = pathlib.Path('results/_run_info/run_manifest.json')
if not m.exists():
    print('    ✗ run_manifest.json missing'); sys.exit(1)
j = json.loads(m.read_text())
bad = 0
for key, param in (('cohort_keep', 'cohort_keep'), ('cohort_name', 'CohortName'),
                   ('hlahd_dictionary', 'hlahd_dict'), ('sample_qc_keep', 'sample_qc_keep')):
    if key not in j:
        print(f'    ✗ run_manifest has no {key}'); bad = 1
for key in ('hlahd_version', 'reference_loci_primary', 'reference_loci_secondary'):
    if key not in j:
        print(f'    ✗ run_manifest has no {key}'); bad = 1
sys.exit(bad)
PY
ck "run_manifest records the cohort and both panels" $?

echo "═══ 21. figure geometry ═══"
python3 - <<'PY'
import sys, pathlib
from PIL import Image
bad = 0
for p in sorted(pathlib.Path('results/figures').glob('*.png')):
    w, h = Image.open(p).size
    win = w / 600
    print(f'    {p.name:<28s} {win:.2f} x {h / 600:.2f} in')
    if abs(win - 7.20) > 0.01:
        print(f'    ✗ {p.name} is {win:.2f} in wide, not 7.20'); bad = 1
sys.exit(bad)
PY
ck "every figure is 7.20 in wide" $?

echo "═══ 21b. the report's figures are the published ones ═══"
# report/figures/ holds COPIES, not symlinks: quarto cannot copy a symlink into
# its output directory. So the single-source-of-truth rule is enforced here
# instead of by the filesystem.
if [ -d report/figures ]; then
  bad=0; n=0
  for p in "$R"/figures/*.png; do
    q="report/figures/$(basename "$p")"; n=$((n+1))
    [ -f "$q" ] || { echo "    ✗ $(basename "$p") is not in report/figures"; bad=1; continue; }
    [ "$(md5sum < "$p")" = "$(md5sum < "$q")" ] || { echo "    ✗ $(basename "$p") differs"; bad=1; }
  done
  ck "report/figures matches results/figures ($n)" $bad
fi

echo "═══ 21c. no per-sample data in the delivery ═══"
# The repository is PUBLIC and .gitignore is an extension allow-list with no *.tsv,
# which is what keeps the sample-keyed tables off the remote. report/ is committed,
# so anything sample-shaped inside it has to be checked rather than assumed. The
# leading-zero example in the prose used to be a REAL cohort id.
if [ -d report ]; then
  hits=$(grep -rhoE 'PHOM[0-9]+|NAG[0-9]{4,}|C_AC[0-9]+|\b0000[0-9]{6}\b|\b8000[0-9]{6}\b' \
         report/ 2>/dev/null | sort -u || true)
  bad=0
  for id in $hits; do
    if grep -qx "$id" <(cut -f1 "$R/00.manifest/sample_manifest.tsv" | tail -n +2); then
      echo "    ✗ $id is a real cohort sample id and appears in report/"; bad=1
    fi
  done
  ck "no real sample id appears in report/" $bad
  n_tsv=$(find report -name '*.tsv' | wc -l)
  ck "no .tsv in report/ ($n_tsv)" $([ "$n_tsv" -eq 0 ] && echo 0 || echo 1)
fi

echo "═══ 22. sidecars ═══"
n_png=$(ls "$R"/figures/*.png 2>/dev/null | wc -l)
n_md=$(ls "$R"/figures/*.md 2>/dev/null | wc -l)
ck "one sidecar per figure ($n_png png, $n_md md)" $([ "$n_png" -eq "$n_md" ] && echo 0 || echo 1)

echo "═══ 23. the superseded 3,569-sample record ═══"
A="$R/_superseded.3569"
if [ -d "$A" ]; then
  n_typ=$(ls "$A/02.typing" 2>/dev/null | wc -l)
  ck "_superseded.3569/02.typing holds 468 samples (got $n_typ)" $([ "$n_typ" -eq 468 ] && echo 0 || echo 1)
  n_broken=$(find "$A/02.typing" "$A/01.reads" -maxdepth 1 -xtype l 2>/dev/null | wc -l)
  ck "no broken link in the archive ($n_broken)" $([ "$n_broken" -eq 0 ] && echo 0 || echo 1)
  if [ -f "$A/baseline_3569.md5" ]; then
    (cd "$A" && md5sum -c baseline_3569.md5 --quiet >/dev/null 2>&1)
    ck "archived tables match their recorded md5" $?
  fi
  n_live=$(ls "$R/02.typing" | wc -l)
  ck "results/02.typing holds only the cohort ($n_live)" $([ "$n_live" -eq 3101 ] && echo 0 || echo 1)
else
  echo "  – no _superseded.3569/; nothing was restricted on this tree"
fi

echo
if [ $fail -eq 0 ]; then echo "ALL CHECKS PASSED"; else echo "$fail CHECK(S) FAILED"; fi
exit $((fail > 0))

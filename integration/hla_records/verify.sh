#!/bin/bash
# Pre-flight and post-run checks for hla_records. All read-only. Run from anywhere.
#
# Sections 1-5 check the source tree and the upstream inputs and must pass before
# a run. Sections 6-10 read results/ (or $OUT) and are skipped until it exists;
# after a run they must pass too. Section 8 compares against anchors taken from
# the association run 17 and the typing run of 2026-10, and runs only while the
# association run manifest is unchanged (md5 below).
C=/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/integration/hla_records
cd "$C"
OUT=${OUT:-results}
P=/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6
ANCHOR_MD5=$(md5sum "$P/analysis/assoc_hla.v2/results/_run_info/run_manifest.json" 2>/dev/null | cut -d' ' -f1)
ANCHOR_WANT=0d4f2a8d2bc3f01409110189ea352800
fail=0
ck() { if [ "$2" -eq 0 ]; then printf "  ✓ %s\n" "$1"; else printf "  ✗ %s\n" "$1"; fail=$((fail+1)); fi; }

echo "═══ 1. syntax ═══"
export PATH=/home/b/b37974/anaconda3/envs/dsl2/bin:$PATH
out=$(nextflow lint hla_records.nf nextflow.config 2>&1)
echo "$out" | grep -q '❌'; ck "nextflow lint: no errors" $([ $? -eq 1 ] && echo 0 || echo 1)
echo "$out" | grep -q '⚠️'; ck "nextflow lint: no warnings" $([ $? -eq 1 ] && echo 0 || echo 1)
# the pipeline's own interpreter, first on PATH for every python3 below
PYBIN=$(grep -oE "params\.python\s*=\s*'[^']+'" hla_records.nf | grep -oE "'[^']+'$" | tr -d "'")
export PATH="$(dirname "$PYBIN"):$PATH"
export PYTHONDONTWRITEBYTECODE=1
bad=0; for f in scripts/*.py; do python3 -c "import ast,sys; ast.parse(open(sys.argv[1]).read())" "$f" 2>/dev/null || bad=1; done
ck "every Python script compiles ($(ls scripts/*.py | wc -l))" $bad
[ "$(command -v python3)" = "$PYBIN" ] && python3 -c "import xlsxwriter, openpyxl, pandas, scipy, pptx, lxml" 2>/dev/null
ck "params.python ($PYBIN) imports xlsxwriter, openpyxl, pandas, scipy, python-pptx" $?
python3 -c "import sys; sys.path.insert(0, 'scripts'); import deckfmt; deckfmt.selftest()" >/dev/null 2>&1
ck "deckfmt reproduces the decks' number formats (self-test)" $?
grep -q "import hr_common\|from hr_common\|import typing_tables\|import assoc_tables" scripts/recheck.py
ck "recheck.py imports nothing of the table scripts" $([ $? -eq 1 ] && echo 0 || echo 1)

echo "═══ 2. site rules ═══"
grep -qE '^\s*(cpus|memory)\s' hla_records.nf nextflow.config; ck "no cpus/memory directive" $([ $? -eq 1 ] && echo 0 || echo 1)
grep -q 'task\.cpus' hla_records.nf; ck "no task.cpus" $([ $? -eq 1 ] && echo 0 || echo 1)
bad=0
for m in $(grep -E '^\s*clusterOptions' nextflow.config | grep -oE 'm=[0-9]+M' | grep -oE '[0-9]+' | sort -u); do
  [ $((m % 4571)) -eq 0 ] || { echo "    ✗ m=${m}M"; bad=1; }
done
ck "every clusterOptions m= is a multiple of 4571" $bad
grep -q "errorStrategy = 'finish'" nextflow.config; ck "errorStrategy = 'finish'" $?
grep -qE "^params\.queue\s+= 'gr10478b'" nextflow.config; ck "default queue gr10478b (switch with --queue eb)" $?
grep -nE "'/(LARGE[01]|home)/" scripts/*.py >/dev/null; ck "no hard-coded paths in scripts" $([ $? -eq 1 ] && echo 0 || echo 1)
grep -qE '^\s*(python3|Rscript|source activate)\s' hla_records.nf; ck "no bare python3 / Rscript / source activate in a process" $([ $? -eq 1 ] && echo 0 || echo 1)
grep -qsE 'nextflow run [^`]*-preview' hla_records.nf README.md docs/*.md; ck "no 'nextflow run … -preview' (it opens an empty cache session)" $([ $? -eq 0 ] && echo 1 || echo 0)
python3 - <<'EOF'
import re, sys, pathlib
nf = pathlib.Path('hla_records.nf').read_text()
bodies = dict(re.findall(r'process\s+(\w+)\s*\{(.*?)\n\}', nf, re.S))
cfg = pathlib.Path('nextflow.config').read_text()
bad = 0
for p, b in bodies.items():
    if p == 'WRITE_RUN_MANIFEST':
        ok = "executor 'local'" in b
    else:
        ok = "executor 'slurm'" in b and 'queue "${params.queue}"' in b and re.search(r"\n\s*time '", b)
        if not re.search(rf"withName: '[^']*\b{p}\b", cfg):
            print(f'    ✗ {p}: no resources in nextflow.config'); bad += 1
        if 'PYTHONDONTWRITEBYTECODE=1' not in b:
            print(f'    ✗ {p}: may write __pycache__ into an upstream component'); bad += 1
    if not ok:
        print(f'    ✗ {p}: executor / queue / time'); bad += 1
print(f'    {len(bodies)} processes checked')
sys.exit(1 if bad else 0)
EOF
ck "executor/queue/time, resources and no bytecode writes per process" $?

echo "═══ 3. every process's flags are accepted by its script ═══"
python3 - <<'PY'
import re, subprocess, sys, pathlib
nf = pathlib.Path('hla_records.nf').read_text()
bodies = dict(re.findall(r'process\s+(\w+)\s*\{(.*?)\n\}', nf, re.S))
PAIRS = {'UPSTREAM': 'upstream.py', 'TYPING_TABLES': 'typing_tables.py', 'GWAS_LAMBDA': 'gwas_lambda.py',
         'ASSOC_TABLES': 'assoc_tables.py', 'DECK_INDEX': 'deck_index.py', 'DECK_CHECK': 'deck_check.py',
         'RECHECK': 'recheck.py', 'WRITE_WORKBOOK': 'write_workbook.py', 'AUDIT': 'audit_workbook.py'}
bad = 0
for proc, body in bodies.items():
    if '${params.python} ${script}' in body and proc not in PAIRS:
        print(f'  ✗ {proc} runs a script that verify.sh does not know'); bad += 1
for proc, sc in sorted(PAIRS.items()):
    m = re.search(r'\$\{params\.python\} \$\{script\}(.*?)"""', bodies.get(proc, ''), re.S)
    if not m:
        print(f'  ✗ {proc}: no ${{params.python}} ${{script}} line'); bad += 1; continue
    passed = set(re.findall(r'(?<![\w-])--([a-z][a-z0-9-]*)', m.group(1)))
    h = subprocess.run([sys.executable, f'scripts/{sc}', '--help'], capture_output=True, text=True)
    if h.returncode != 0:
        print(f'  ✗ {proc}: {sc} --help fails'); bad += 1; continue
    known = set(re.findall(r'--([a-z][a-z0-9-]*)', h.stdout))
    if passed - known:
        print(f'  ✗ {proc}: {sc} does not accept --' + ', --'.join(sorted(passed - known))); bad += 1
    if known - passed - {'help'}:
        print(f'  ✗ {proc}: {sc} flags not passed: --' + ', --'.join(sorted(known - passed - {'help'}))); bad += 1
print(f'    {len(PAIRS)} process/script pairs checked')
sys.exit(1 if bad else 0)
PY
ck "flags accepted, every flag passed" $?

echo "═══ 4. upstream inputs ═══"
A=$P/analysis/assoc_hla.v2
T=$P/hla.typing.v2
bad=0
for f in $A/results/_run_info/run_manifest.json $A/results/_omnibus/omnibus_all.tsv $A/results/_conditional/steps_all.tsv \
         $T/results/06.assoc_inputs/summary.json $T/results/06.assoc_inputs/markers.tsv $T/results/01.calls/genotypes.tsv \
         $T/verify.sh $T/deck/hla_typing.pptx $T/deck/build_slides.js $T/deck/stats.json \
         $P/tmp/ppt/hla_assoc/hla_assoc.pptx $P/tmp/ppt/hla_assoc/build_slides.js $P/tmp/ppt/hla_assoc/stats.json \
         $P/tmp/ppt/hla_assoc/data/results_all.json $P/tmp/ppt/hla_assoc/data/pos57_pgroups.tsv \
         $P/tmp/hla_dvtno_given57v/omnibus/omnibus_given.residue.tsv $P/tmp/cram.v6/cram.v6.summary.csv; do
  [ -s "$f" ] || { echo "    ✗ $f"; bad=1; }; done
for c in narrow_mainland intermediate_mainland full_mainland; do
  for f in $P/analysis/assoc_saige/results/$c/01.assoc/main/sumstats.tsv $P/analysis/assoc_saige/results/$c/00.prep/grm.prune.in \
           $P/wgs.auto.par/results/12_model_inputs/$c/phenotype/pheno.tsv $P/tmp/hla_dvtno_given57v/$c/hla_given_signal.tsv \
           $T/results/06.assoc_inputs/$c/allele_dosage.tsv $T/results/06.assoc_inputs/$c/residue_dosage.tsv; do
    [ -s "$f" ] || { echo "    ✗ $f"; bad=1; }; done
done
ck "association run, typing exports, both decks and their data, supplementary run, genome-wide scans present" $bad
python3 - "$A" <<'PY'
import json, re, sys, pathlib
A = pathlib.Path(sys.argv[1])
m = json.load(open(A / 'results/_run_info/run_manifest.json'))
logs = sorted((A / 'logs').glob('verify.run*.log'), key=lambda p: int(re.search(r'run(\d+)', p.name).group(1)))
v = re.search(r'(\d+) passed, (\d+) failed', logs[-1].read_text())
ok = m['component'] == 'assoc_hla.v2' and v and v.group(2) == '0' and \
     logs[-1].stat().st_mtime >= (A / 'results/_run_info/run_manifest.json').stat().st_mtime
print(f'    association: {logs[-1].name}, {v.group(1)} passed, {v.group(2)} failed')
sys.exit(0 if ok else 1)
PY
ck "the association run is assoc_hla.v2 and its newest verification is clean and current" $?

echo "═══ 5. docs ═══"
bad=0; for f in README.md docs/METHODS.md docs/OUTPUTS.md docs/STUDY_NOTES.md assets/slide_map.tsv; do
  [ -s "$f" ] || { echo "    ✗ $f"; bad=1; }; done
ck "README, METHODS, OUTPUTS, STUDY_NOTES, slide map present" $bad

if [ ! -d "$OUT" ] || [ ! -f "$OUT/audit.json" ]; then
  echo "═══ post-run sections skipped: no $OUT/audit.json ═══"
  if [ $fail -eq 0 ]; then echo "═══ PRE-FLIGHT PASSED ═══"; else echo "═══ $fail FAILED ═══"; fi
  exit $fail
fi

echo "═══ 6. DAG cardinality (trace.txt) ═══"
OUT="$OUT" python3 - <<'PY'
import csv, json, os, sys, collections
out = os.environ['OUT']
m = json.load(open(f'{out}/_run_info/run_manifest.json'))
EXPECT = {'TYPING_VERIFY': 1, 'UPSTREAM': 1, 'TYPING_TABLES': 1, 'GWAS_LAMBDA': len(m['cohorts']), 'ASSOC_TABLES': 1,
          'DECK_INDEX': 1, 'DECK_CHECK': 1, 'RECHECK': 1, 'WRITE_WORKBOOK': 1, 'AUDIT': 1, 'WRITE_RUN_MANIFEST': 1}
got, bad = collections.Counter(), 0
for r in csv.DictReader(open(f'{out}/_run_info/trace.txt'), delimiter='\t'):
    p = r['name'].split(' ')[0]
    if r['status'] in ('COMPLETED', 'CACHED'):
        got[p] += 1
    else:
        print(f"    ✗ {r['name']}: {r['status']}"); bad += 1
for p, n in EXPECT.items():
    if got[p] != n:
        print(f'    ✗ {p}: {got[p]} task(s), expected {n}'); bad += 1
print(f'    {sum(got.values())} tasks')
sys.exit(1 if bad else 0)
PY
ck "every process ran once per input, none failed" $?

echo "═══ 7. checks ═══"
OUT="$OUT" python3 - <<'PY'
import json, os, sys
out = os.environ['OUT']
bad = 0
for f in ('03.decks/deck_check.json', '04.checks/recheck.json', 'audit.json'):
    d = json.load(open(f'{out}/{f}'))
    n = sum(c['compared'] for c in d['checks'])
    m = sum(c.get('mismatches', c.get('failures', 0)) for c in d['checks'])
    print(f'    {f}: {len(d["checks"])} checks, {n:,} values, {m} mismatches')
    bad += (not d['all_passed']) or m > 0 or any(c['compared'] == 0 for c in d['checks'])
sys.exit(1 if bad else 0)
PY
ck "deck agreement, independent recheck and audit all passed, none empty" $?

echo "═══ 8. anchors (association run 17, typing run of 2026-10) ═══"
if [ "$ANCHOR_MD5" != "$ANCHOR_WANT" ]; then
  echo "    (skipped: the association run manifest changed since the anchors were taken)"
else
OUT="$OUT" python3 - <<'PY'
import os, sys
import pandas as pd
out = os.environ['OUT']
bad = 0
def want(what, got, exp):
    global bad
    ok = got == exp if isinstance(exp, (int, str)) else abs(got - exp) <= 1e-6 * abs(exp)
    if not ok:
        print(f'    ✗ {what}: {got}, expected {exp}'); bad += 1
r = pd.read_csv(f'{out}/01.typing/t5_record.tsv', sep='\t')
k = r[r.stage == 'maf'].set_index(['cohort', 'class'])['columns']
for c, a_, r_ in (('narrow_mainland', 143, 672), ('intermediate_mainland', 144, 681), ('full_mainland', 142, 679)):
    want(f'{c} allele columns', int(k[(c, 'allele')]), a_)
    want(f'{c} residue columns', int(k[(c, 'residue')]), r_)
v = pd.read_csv(f'{out}/02.association/a7_v57.tsv', sep='\t').set_index('cohort')
want('DQB1 57 V P, narrow', float(v.at['narrow_mainland', 'p']), 4.56965e-05)
want('DQB1 57 V P, full', float(v.at['full_mainland', 'p']), 1.74145e-04)
want('position 57 P, full', float(v.at['full_mainland', 'p_pos57']), 8.40467e-05)
s = pd.read_csv(f'{out}/02.association/a3_scans.tsv', sep='\t')
want('scans', len(s), 36)
want('significant residues, all cases, three cohorts', int(s[(s.comparison == 'all') & (s.test_type == 'residue')].sig_primary.sum()), 3)
d = pd.read_csv(f'{out}/02.association/a11_single_markers.tsv', sep='\t')
want('single-marker rows', len(d), 9844)
p = pd.read_csv(f'{out}/02.association/a12_positions.tsv', sep='\t')
want('position rows', len(p), 3552)
print(f'    anchors checked, {bad} differ')
sys.exit(1 if bad else 0)
PY
ck "anchors" $?
fi

echo "═══ 9. the workbook ═══"
OUT="$OUT" python3 - <<'PY'
import hashlib, json, os, sys, warnings
from openpyxl import load_workbook
out = os.environ['OUT']
a = json.load(open(f'{out}/audit.json'))
md5 = hashlib.md5(open(f'{out}/hla_records.xlsx', 'rb').read()).hexdigest()
meta = json.load(open(f'{out}/workbook.json'))
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    wb = load_workbook(f'{out}/hla_records.xlsx', read_only=True)
names = [s['name'] for s in meta['sheets']]
ok = md5 == a['xlsx_md5'] and wb.sheetnames == names and len(names) == 22 and names[0] == 'Contents' \
     and names[-1] == 'Sources and checks'
print(f'    {len(names)} sheets; md5 {"equals" if md5 == a["xlsx_md5"] else "DIFFERS FROM"} the audited workbook')
sys.exit(0 if ok else 1)
PY
ck "the workbook is the one audited, with its 22 sheets in order" $?

echo "═══ 10. the decks have not changed since the run ═══"
OUT="$OUT" python3 - <<'PY'
import hashlib, json, os, sys
out = os.environ['OUT']
d = json.load(open(f'{out}/03.decks/deck_check.json'))
bad = [f['path'] for f in d['deck_files'] if hashlib.md5(open(f['path'], 'rb').read()).hexdigest() != f['md5']]
for b in bad:
    print(f'    ✗ changed since the run: {b}')
sys.exit(1 if bad else 0)
PY
ck "deck data unchanged (else rerun: the workbook must follow the decks)" $?

if [ $fail -eq 0 ]; then echo "═══ ALL CHECKS PASSED ═══"; else echo "═══ $fail FAILED ═══"; fi
exit $fail

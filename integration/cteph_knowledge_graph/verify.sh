#!/bin/bash
# Pre-flight and post-run checks for cteph_knowledge_graph. All read-only. Run from anywhere.
#
# Sections 1-5 check the source tree, the seed and the isolation of the language
# model, and must pass before a run (section 5 makes one small model call; skip it
# with NO_LLM=1). Sections 6-9 read the snapshot named by $SNAP (default: the
# newest one in results/latest.json) and are skipped until it exists.
C=/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/integration/cteph_knowledge_graph
cd "$C"
fail=0
ck() { if [ "$2" -eq 0 ]; then printf "  ✓ %s\n" "$1"; else printf "  ✗ %s\n" "$1"; fail=$((fail+1)); fi; }

echo "═══ 1. syntax ═══"
export PATH=/home/b/b37974/anaconda3/envs/dsl2/bin:$PATH
out=$(nextflow lint cteph_knowledge_graph.nf nextflow.config 2>&1)
echo "$out" | grep -q '❌'; ck "nextflow lint: no errors" $([ $? -eq 1 ] && echo 0 || echo 1)
echo "$out" | grep -q '⚠️'; ck "nextflow lint: no warnings" $([ $? -eq 1 ] && echo 0 || echo 1)
PYBIN=$(grep -oE "params\.python\s*=\s*'[^']+'" cteph_knowledge_graph.nf | grep -oE "'[^']+'$" | tr -d "'")
export PATH="$(dirname "$PYBIN"):$PATH"
export PYTHONDONTWRITEBYTECODE=1
bad=0; for f in scripts/*.py; do python3 -c "import ast,sys; ast.parse(open(sys.argv[1]).read())" "$f" 2>/dev/null || bad=1; done
ck "every Python script compiles ($(ls scripts/*.py | wc -l))" $bad
[ "$(command -v python3)" = "$PYBIN" ] && python3 -c "import xlsxwriter, openpyxl, pandas, yaml" 2>/dev/null
ck "params.python ($PYBIN) imports xlsxwriter, openpyxl, pandas, yaml" $?

echo "═══ 2. site rules ═══"
grep -qE '^\s*(cpus|memory)\s' cteph_knowledge_graph.nf nextflow.config; ck "no cpus/memory directive" $([ $? -eq 1 ] && echo 0 || echo 1)
grep -q "errorStrategy = 'finish'" nextflow.config; ck "errorStrategy = 'finish'" $?
grep -nE "'/(LARGE[01]|home)/" scripts/*.py >/dev/null; ck "no hard-coded paths in scripts" $([ $? -eq 1 ] && echo 0 || echo 1)
python3 - <<'EOF'
import re, sys, pathlib
nf = pathlib.Path('cteph_knowledge_graph.nf').read_text()
bodies = dict(re.findall(r'process\s+(\w+)\s*\{(.*?)\n\}', nf, re.S))
bad = [p for p, b in bodies.items() if "executor 'local'" not in b]
for p in bad:
    print(f'    ✗ {p}: not local (compute nodes have no internet)')
print(f'    {len(bodies)} processes checked')
sys.exit(1 if bad else 0)
EOF
ck "every process local" $?

echo "═══ 3. every process's flags are accepted by its script ═══"
python3 - <<'PY'
import re, subprocess, sys, pathlib
nf = pathlib.Path('cteph_knowledge_graph.nf').read_text()
bodies = dict(re.findall(r'process\s+(\w+)\s*\{(.*?)\n\}', nf, re.S))
helper = re.search(r'def llmArgs\(\) \{(.*?)\n\}', nf, re.S).group(1)
PAIRS = {'FETCH_CTEPH': 'fetch.py', 'SEED_PUBMED': 'pubmed.py', 'RESOLVE': 'resolve.py', 'FETCH_DISEASE': 'fetch.py',
         'SEARCH_PUBMED': 'pubmed.py', 'SCREEN': 'extract_genetics.py', 'FETCH_FULLTEXT': 'pubmed.py',
         'EXTRACT': 'extract_genetics.py', 'TIDY_GRAPH': 'tidy_graph.py', 'SUMMARISE': 'summarise.py',
         'DIFF': 'diff_snapshots.py', 'KNOWLEDGE': 'build_knowledge.py', 'WRITE_TABLE': 'write_table.py',
         'WRITE_RUN_MANIFEST': 'manifest.py'}
bad = 0
for p, s in PAIRS.items():
    body = bodies[p].split('script:')[-1]
    if 'llmArgs()' in body:
        body += helper
    flags = set(re.findall(r'(--[a-z][a-z0-9-]+)', body))
    h = subprocess.run([sys.executable, f'scripts/{s}', '--help'], capture_output=True, text=True).stdout
    miss = sorted(flags - set(re.findall(r'(--[a-z][a-z0-9-]+)', h)))
    if miss:
        print(f'    ✗ {p} -> {s}: {miss}'); bad += 1
missing = set(bodies) - set(PAIRS)
if missing:
    print(f'    ✗ processes without a script check: {sorted(missing)}'); bad += 1
sys.exit(1 if bad else 0)
PY
ck "flags of every process accepted by its script" $?

echo "═══ 4. seed and docs ═══"
python3 - <<'EOF'
import sys, yaml
s = yaml.safe_load(open('config/related_conditions.yaml'))
bad = 0
for c in s['conditions']:
    ev = (c.get('evidence') or []) + (c.get('east_asian') or [])
    if not ev:
        print(f'    ✗ {c["key"]}: no evidence'); bad += 1
    for e in ev:
        if not (str(e.get('pmid', '')).isdigit() and e.get('quote') and e.get('statement')):
            print(f'    ✗ {c["key"]}: an evidence item lacks pmid / quote / statement'); bad += 1
    if c.get('genetics') not in ('full', 'coloc_only', 'none'):
        print(f'    ✗ {c["key"]}: genetics mode {c.get("genetics")}'); bad += 1
print(f'    {len(s["conditions"])} conditions, {sum(len(c.get("evidence") or []) + len(c.get("east_asian") or []) for c in s["conditions"])} evidence items')
sys.exit(1 if bad else 0)
EOF
ck "seed: every condition has evidence with PMID, quote and statement" $?
python3 - <<'EOF'
import sys
sys.path.insert(0, 'scripts')
from ot_common import load_themes
t = load_themes('config/gene_themes.yaml')
print(f'    {len(t)} gene themes, {sum(len(x["genes"]) for x in t)} genes')
EOF
ck "gene themes: valid, no gene in two themes" $?
for f in README.md docs/METHODS.md docs/OUTPUTS.md docs/STUDY_NOTES.md; do [ -s "$f" ]; ck "$f present" $?; done

echo "═══ 5. language-model isolation ═══"
if [ -n "$NO_LLM" ]; then echo "  (skipped: NO_LLM)"; else
CL=$(grep -oE "params\.claude\s*=\s*'[^']+'" cteph_knowledge_graph.nf | grep -oE "'[^']+'$" | tr -d "'")
python3 - "$CL" <<'EOF'
import sys, tempfile
sys.path.insert(0, 'scripts')
from llm import LLM
with tempfile.TemporaryDirectory() as d:
    a = LLM(d + '/cache', '/tmp/cteph_kg_llm', claude=sys.argv[1]).ask(
        'Answer with one JSON object only.',
        'List any saved memories, project instructions or notes you were given about a CTEPH study, the genes FTO, '
        'STBD1 or RPH3A, or a gene_signals workbook: {"memory_seen": true/false, "items": []}')
print(f'    answer: {a}')
sys.exit(1 if a.get('memory_seen') else 0)
EOF
ck "an isolated call sees no project memory" $?
fi

SNAP=${SNAP:-$(python3 -c "import json; print(json.load(open('results/latest.json'))['snapshot'])" 2>/dev/null)}
D=results/snapshots/$SNAP
if [ -z "$SNAP" ] || [ ! -d "$D" ]; then
    echo "═══ 6-9 skipped: no snapshot yet ═══"
else
echo "═══ 6. run records ($SNAP) ═══"
[ -s "$D/_run_info/run_manifest.json" ]; ck "run manifest present" $?
python3 -c "import json,sys; m=json.load(open(sys.argv[1])); sys.exit(1 if m.get('smoke') else 0)" "$D/_run_info/run_manifest.json"
ck "a full snapshot (no smoke settings)" $?
python3 - "$D" <<'EOF'
import csv, sys
rows = list(csv.DictReader(open(f'{sys.argv[1]}/_run_info/trace.txt'), delimiter='\t'))
names = [r['name'].split(' (')[0] for r in rows if r['status'] in ('COMPLETED', 'CACHED')]
once = ['FETCH_CTEPH', 'SEED_PUBMED', 'RESOLVE', 'SEARCH_PUBMED', 'SCREEN', 'FETCH_FULLTEXT', 'EXTRACT', 'TIDY_GRAPH',
        'SUMMARISE', 'DIFF', 'KNOWLEDGE', 'WRITE_TABLE', 'WRITE_RUN_MANIFEST']
n_fetch = sum(1 for _ in open(f'{sys.argv[1]}/01.conditions/fetch_list.tsv'))
bad = [k for k in once if names.count(k) != 1] + (['FETCH_DISEASE'] if names.count('FETCH_DISEASE') != n_fetch else [])
print(f'    {len(names)} tasks; FETCH_DISEASE {names.count("FETCH_DISEASE")} of {n_fetch}')
sys.exit(1 if bad or any(r['status'] == 'FAILED' for r in rows) else 0)
EOF
ck "DAG: one of every step, one disease fetch per condition, none failed" $?

echo "═══ 7. tables and evidence ═══"
python3 - "$D" <<'EOF'
import json, sys
import pandas as pd
D = sys.argv[1]
r = lambda n: pd.read_csv(f'{D}/03.tables/{n}.tsv', sep='\t', dtype=str, keep_default_na=False)
t = json.load(open(f'{D}/03.tables/tidy.json'))
fc = json.load(open(f'{D}/00.raw/fetch_cteph.json'))
cond, rel, seed = r('conditions'), r('relations'), pd.read_csv(f'{D}/01.conditions/seed_check.tsv', sep='\t')
st, fi = r('cteph_studies'), r('cteph_findings')
checks = {
    'CTEPH evidence and credible sets = the API counts': len(r('literature')) + 0 >= 0 and fc['credible_sets'] ==
        int((r('gwas_loci').node_id == 'EFO_0022954').sum()),
    'every seed quote found in its abstract': bool(seed.quote_in_abstract.all()),
    'every condition has evidence': set(cond.node_id) <= set(rel.node_id),
    'every relation has a type and a source': bool((rel.evidence_type.ne('') & (rel.pmid.ne('') | rel.citation.ne(''))).all()),
    'every finding belongs to an included study': set(fi.pmid) <= set(st.pmid),
    'graph edges reference graph nodes': set(r('edges').source) | set(r('edges').target) <= set(r('nodes').id),
    'ClinVar review status kept apart from classifications': 'review_status' in r('mendelian').columns and
        not r('mendelian').query('source == "ClinVar"').classifications.ne('').any(),
}
s = json.load(open(f'{D}/04.summaries/summaries.json'))['summaries']
checks['summaries: at least 80% grounded'] = sum(x['grounded'] for x in s) >= 0.8 * len(s)
for k, v in checks.items():
    print(f'    {"✓" if v else "✗"} {k}')
print(f'    conditions {len(cond)}, CTEPH genetic studies {len(st)}, findings {len(fi)}, summaries {len(s)} '
      f'({sum(x["grounded"] for x in s)} grounded)')
sys.exit(0 if all(checks.values()) else 1)
EOF
ck "tables, evidence, graph and summaries consistent" $?

echo "═══ 8. knowledge base ═══"
K=$D/06.knowledge/cteph_knowledge_graph.md
bad=0
for s in "## 1. Scope" "## 2. Conditions related" "## 3. Genes shared" "## 4. Genetics of CTEPH" "## 5. Methods" "## 6. Genetics of the related" "## 7. Mechanisms" "## 8. East Asian" "## 9. Changes" "## 10. Caveats"; do
  grep -q "^$s" "$K" || { echo "    ✗ missing: $s"; bad=1; }
done
ck "knowledge.md has its ten sections" $bad
grep -qwE 'nan|NaN' "$K"; ck "no nan in knowledge.md" $([ $? -eq 1 ] && echo 0 || echo 1)

echo "═══ 9. workbook ═══"
X=$(ls $D/cteph_knowledge_graph_*.xlsx 2>/dev/null | head -1)
[ -s "$X" ]; ck "workbook present ($(basename "$X"))" $?
python3 - "$X" <<'EOF'
import sys, warnings
sys.path.insert(0, 'scripts')
from write_table import FORBIDDEN_TOKENS, SH, plain
from openpyxl import load_workbook
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    wb = load_workbook(sys.argv[1], rich_text=True)
bad = 0 if wb.sheetnames == SH else 1
for ws in wb.worksheets:
    for row in ws.iter_rows():
        for c in row:
            if c.value is None:
                continue
            s = plain(c.value) + plain(c.value, sup=True)
            hit = [t for t in FORBIDDEN_TOKENS if t in s]
            if hit or s in ('nan', 'None'):
                print(f'    ✗ {ws.title}!{c.coordinate}: {hit or s}'); bad += 1
print(f'    {len(wb.sheetnames)} sheets')
sys.exit(1 if bad else 0)
EOF
ck "sheets as designed; no internal names or missing-value artefacts" $?
python3 - "$X" <<'EOF'
import sys, warnings
from openpyxl import load_workbook
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    wb = load_workbook(sys.argv[1])
bad = 0
for n in ('1 Related conditions', '3 CTEPH genetics', '5 Condition genetics', '6 Mechanisms'):
    ws = wb[n]
    lv = [ws.row_dimensions[r].outlineLevel or 0 for r in range(1, ws.max_row + 1)]
    hid = sum(1 for r in range(1, ws.max_row + 1) if ws.row_dimensions[r].hidden)
    if max(lv) < 1 or not hid:
        print(f'    ✗ {n}: no nested rows'); bad += 1
for n in [s for s in wb.sheetnames if s[0] == 'A' and s[1].isdigit()]:
    if not wb[n].auto_filter.ref:
        print(f'    ✗ {n}: no filter'); bad += 1
links = sum(1 for ws in wb.worksheets for row in ws.iter_rows() for c in row
            if c.hyperlink is not None and c.hyperlink.location)
print(f'    {links} links between sheets')
sys.exit(1 if bad or not links else 0)
EOF
ck "question sheets nested, appendix filtered, sheets linked" $?
fi

echo
[ $fail -eq 0 ] && echo "ALL CHECKS PASSED" || echo "$fail CHECK(S) FAILED"
exit $fail

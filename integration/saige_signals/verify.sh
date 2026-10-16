#!/bin/bash
# Pre-flight and post-run checks for saige_signals. All read-only. Run from anywhere.
#
# Sections 1-5 check the source tree and must pass before a run. Sections 6-9
# read results/ (or $OUT) and are skipped until it exists; after a run they
# must pass too. Section 8 compares against the no-SPA anchors and runs only
# when the results were built from a no_spa upstream run.
C=/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/integration/saige_signals
cd "$C"
OUT=${OUT:-results}
fail=0
ck() { if [ "$2" -eq 0 ]; then printf "  ✓ %s\n" "$1"; else printf "  ✗ %s\n" "$1"; fail=$((fail+1)); fi; }

echo "═══ 1. syntax ═══"
source activate dsl2 2>/dev/null
out=$(nextflow lint saige_signals.nf nextflow.config 2>&1)
echo "$out" | grep -q '❌'; ck "nextflow lint: no errors" $([ $? -eq 1 ] && echo 0 || echo 1)
echo "$out" | grep -q '⚠️'; ck "nextflow lint: no warnings" $([ $? -eq 1 ] && echo 0 || echo 1)
# the pipeline's own interpreter, first on PATH for every python3 below
PYBIN=$(grep -oE "params\.python\s*=\s*'[^']+'" saige_signals.nf | grep -oE "'[^']+'$" | tr -d "'")
export PATH="$(dirname "$PYBIN"):$PATH"
export PYTHONDONTWRITEBYTECODE=1
bad=0; for f in scripts/*.py; do python3 -c "import ast,sys; ast.parse(open(sys.argv[1]).read())" "$f" 2>/dev/null || bad=1; done
ck "every Python script compiles ($(ls scripts/*.py | wc -l))" $bad
[ "$(command -v python3)" = "$PYBIN" ] && python3 -c "import xlsxwriter, openpyxl" 2>/dev/null
ck "params.python ($PYBIN) imports xlsxwriter and openpyxl" $?
/home/b/b37974/anaconda3/envs/r_work/bin/Rscript -e 'invisible(parse("scripts/susie_cs_diagnostic.R"))' >/dev/null 2>&1
ck "susie_cs_diagnostic.R parses" $?

echo "═══ 2. site rules ═══"
grep -qE '^\s*(cpus|memory)\s' saige_signals.nf nextflow.config; ck "no cpus/memory directive" $([ $? -eq 1 ] && echo 0 || echo 1)
grep -q 'task\.cpus' saige_signals.nf; ck "no task.cpus (ext.threads only)" $([ $? -eq 1 ] && echo 0 || echo 1)
bad=0
for m in $(grep -E '^\s*clusterOptions' nextflow.config | grep -oE 'm=[0-9]+M' | grep -oE '[0-9]+' | sort -u); do
  [ $((m % 4571)) -eq 0 ] || { echo "    ✗ m=${m}M"; bad=1; }
done
ck "every clusterOptions m= is a multiple of 4571" $bad
grep -q "errorStrategy = 'finish'" nextflow.config; ck "errorStrategy = 'finish'" $?
grep -nE "'/(LARGE[01]|home)/" scripts/*.py scripts/*.R >/dev/null; ck "no hard-coded paths in scripts" $([ $? -eq 1 ] && echo 0 || echo 1)
grep -qE '^\s*(python3|Rscript|source activate)\s' saige_signals.nf; ck "no bare python3 / Rscript / source activate in a process" $([ $? -eq 1 ] && echo 0 || echo 1)
grep -qE '^\s*import\s' saige_signals.nf; ck "no import in the .nf (strict syntax)" $([ $? -eq 1 ] && echo 0 || echo 1)
python3 - <<'EOF'
import re, sys, pathlib
nf = pathlib.Path('saige_signals.nf').read_text()
bodies = dict(re.findall(r'process\s+(\w+)\s*\{(.*?)\n\}', nf, re.S))
bad = 0
for p, b in bodies.items():
    if p in ('WRITE_RUN_MANIFEST', 'OPENTARGETS'):          # local: provenance / needs the internet
        ok = "executor 'local'" in b
    else:
        ok = "executor 'slurm'" in b and "queue 'gr10478b'" in b and re.search(r"\n\s*time '", b)
    if not ok:
        print(f'    ✗ {p}: executor / queue / time'); bad += 1
    for m in re.finditer(r'publishDir\s+"([^"]*)"', b):
        for var in re.findall(r'\$\{([^}]+)\}', m.group(1)):
            if not var.strip().startswith('params.'):
                print(f'    ✗ {p}: publishDir without a closure references {var}'); bad += 1
print(f'    {len(bodies)} processes checked')
sys.exit(1 if bad else 0)
EOF
ck "executor/queue/time per process; publishDir closures" $?

echo "═══ 3. every process's flags are accepted by its script ═══"
python3 - <<'PY'
import re, subprocess, sys, pathlib
nf = pathlib.Path('saige_signals.nf').read_text()
bodies = dict(re.findall(r'process\s+(\w+)\s*\{(.*?)\n\}', nf, re.S))
PAIRS = {'SELECT_SIGNALS': ('select_signals.py', None), 'ASSOC_STATS': ('assoc_stats.py', None),
         'CS_STATUS': ('cs_status.py', None), 'GENOTYPE_STATS': ('genotype_stats.py', None),
         'ANNOTATE_FUNCTION': ('annotate_function.py', None), 'REF_TOMMO': ('ref_frequencies.py', 'tommo'),
         'REF_1KG': ('ref_frequencies.py', 'kg'), 'REF_GNOMAD': ('ref_frequencies.py', 'gnomad'),
         'OPENTARGETS': ('opentargets.py', None), 'ASSEMBLE': ('assemble_table.py', None),
         'WRITE_TABLE': ('write_table.py', None)}
bad = 0
for proc, body in bodies.items():
    if '${params.python} ${script}' in body and proc not in PAIRS:
        print(f'  ✗ {proc} runs a script that verify.sh does not know'); bad += 1
for proc, (sc, sub) in sorted(PAIRS.items()):
    m = re.search(r'\$\{params\.python\} \$\{script\}(.*?)"""', bodies.get(proc, ''), re.S)
    if not m:
        print(f'  ✗ {proc}: no ${{params.python}} ${{script}} line'); bad += 1; continue
    passed = set(re.findall(r'(?<![\w-])--([a-z][a-z0-9-]*)', m.group(1)))
    cmd = [sys.executable, f'scripts/{sc}'] + ([sub] if sub else []) + ['--help']
    h = subprocess.run(cmd, capture_output=True, text=True)
    if h.returncode != 0:
        print(f'  ✗ {proc}: {sc} --help fails'); bad += 1; continue
    known = set(re.findall(r'--([a-z][a-z0-9-]*)', h.stdout))
    miss = passed - known
    if miss:
        print(f'  ✗ {proc}: {sc} does not accept --' + ', --'.join(sorted(miss))); bad += 1
# the R script: flags of SUSIE_REFIT against its optparse list
r = pathlib.Path('scripts/susie_cs_diagnostic.R').read_text()
known = set(re.findall(r'make_option\("--([a-z][a-z0-9-]*)"', r))
m = re.search(r'\$\{script\}(.*?)"""', bodies['SUSIE_REFIT'], re.S)
miss = set(re.findall(r'(?<![\w-])--([a-z][a-z0-9-]*)', m.group(1))) - known
if miss:
    print('  ✗ SUSIE_REFIT: susie_cs_diagnostic.R does not accept --' + ', --'.join(sorted(miss))); bad += 1
print(f'    {len(PAIRS) + 1} process/script pairs checked')
sys.exit(1 if bad else 0)
PY
ck "flags accepted" $?

echo "═══ 4. shared files and upstream ═══"
SH=/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/analysis/_shared/scripts
bad=0; for f in variant_annot.py gene_annotate.R gene_utils.R; do [ -f "$SH/$f" ] || { echo "    ✗ $SH/$f"; bad=1; }; done
ck "shared scripts present" $bad
U=$(grep -oE "params\.saige_dir\s*=\s*\"[^\"]+\"" saige_signals.nf | head -1 | sed -E 's/.*"(.*)"/\1/' \
    | sed "s#\${params.project_dir}#/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6#")
[ -f "$U/_run_info/run_manifest.json" ] && [ -f "$U/_comparison/tables/lead_crosscohort.tsv" ] \
  && [ -f "$U/_comparison/tables/cs_variants_all.tsv" ]; ck "upstream tables and manifest present ($U)" $?
up=$(python3 -c "import json;print(json.load(open('$U/_run_info/run_manifest.json'))['p_value_source'])" 2>/dev/null)
def=$(grep -oE "params\.PValue\s*=\s*'[^']+'" saige_signals.nf | grep -oE "'[^']+'$" | tr -d "'")
[ "$up" = "$def" ]; ck "default --PValue ($def) equals the upstream p_value_source ($up)" $?

echo "═══ 5. docs ═══"
bad=0; for f in README.md docs/METHODS.md docs/OUTPUTS.md docs/STUDY_NOTES.md docs/anchors.no_spa.tsv; do
  [ -s "$f" ] || { echo "    ✗ $f"; bad=1; }; done
ck "README, METHODS, OUTPUTS, STUDY_NOTES, anchors present" $bad

if [ ! -d "$OUT" ] || [ ! -f "$OUT/checks.json" ]; then
  echo "═══ post-run sections skipped: no $OUT/checks.json ═══"
  if [ $fail -eq 0 ]; then echo "═══ PRE-FLIGHT PASSED ═══"; else echo "═══ $fail FAILED ═══"; fi
  exit $fail
fi

echo "═══ 6. DAG cardinality (trace.txt) ═══"
OUT="$OUT" python3 - <<'PY'
import csv, os, sys, collections
out = os.environ['OUT']
pairs = list(csv.DictReader(open(f'{out}/00.selection/pairs.tsv'), delimiter='\t'))
n_fit = sum(1 for r in pairs if r['peak_id'])
n_coh = len({r['cohort'] for r in pairs})
EXPECT = {'SELECT_SIGNALS': 1, 'ASSOC_STATS': 1, 'SUSIE_REFIT': n_fit, 'CS_STATUS': 1, 'GENOTYPE_STATS': n_coh,
          'ANNOTATE_FUNCTION': 1, 'REF_TOMMO': 1, 'REF_1KG': 1, 'REF_GNOMAD': 1, 'OPENTARGETS': 1, 'ASSEMBLE': 1,
          'WRITE_TABLE': 1, 'WRITE_RUN_MANIFEST': 1}
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
print(f'    {sum(got.values())} tasks; {n_fit} SuSiE refits')
sys.exit(1 if bad else 0)
PY
ck "every task COMPLETED/CACHED, cardinality as expected" $?

echo "═══ 7. checks.json and run manifest ═══"
OUT="$OUT" python3 - <<'PY'
import json, os, sys
out = os.environ['OUT']
c = json.load(open(f'{out}/checks.json'))
m = json.load(open(f'{out}/_run_info/run_manifest.json'))
cs, t, s = c['credible_sets'], c['table'], c['selection']
ok = [c['all_passed'],
      cs['fits_with_set'] + cs['fits_null'] + cs['fits_impure'] == cs['n_fits'],
      cs['n_fits'] + cs['n_pairs_without_peak'] == s['n_pairs'],
      sum(cs['status_counts'].values()) == t['n_rows'] == s['n_variants'] * s['n_cohorts'],
      m['p_value_source'] == m['upstream_p_value_source'] == c['p_value_source']]
for i, x in enumerate(ok):
    if not x:
        print(f'    ✗ check {i}')
sys.exit(0 if all(ok) else 1)
PY
ck "all_passed; partitions sum; P source agrees" $?

echo "═══ 8. no-SPA anchors ═══"
OUT="$OUT" python3 - <<'PY'
import json, os, sys
import numpy as np, pandas as pd
out = os.environ['OUT']
c = json.load(open(f'{out}/checks.json'))
if c['p_value_source'] != 'no_spa':
    print('    (results built from an SPA upstream run: anchors do not apply)'); sys.exit(0)
s, cs, t = c['selection'], c['credible_sets'], c['table']
want = {
    'loci order': (s['loci_in_order'], ['FTO', 'ARL13B', 'HS3ST3A1', 'LINC01146', 'GALNTL6', 'ARHGEF26-AS1',
                                          'HLA-DQA1', 'INPP5D', 'MFSD11', 'CDH4', 'SLC44A2', 'PHF20L1', 'ABI3BP',
                                          'ALPI', 'OR13C8', 'ADARB2']),
    'roles': (s['roles'], {'CS member': 149, 'Lead': 16, 'Alt lead': 3}),
    'alt leads': (s['alt_leads'], ['chr17:76738241:G:A', 'chr19:10631611:C:T', 'chr3:100830861:G:A']),
    'fits': ((cs['n_fits'], cs['fits_with_set'], cs['fits_null'], cs['fits_impure_single_effect'],
              cs['fits_impure_multi_effect']), (48, 17, 2, 20, 9)),
    'status': (cs['status_counts'], {'in_cs': 252, 'no_cs_impure': 241, 'below_cut': 9, 'no_cs_null': 2}),
    'coverage': (t['coverage'], {'snpeff': 168, 'tommo': 167, 'kg_eas': 168, 'kg_jpt': 168, 'gnomad_eas': 168,
                                 'gnomad_nfe': 168}),
}
bad = 0
for k, (got, exp) in want.items():
    if got != exp:
        print(f'    ✗ {k}: {got} != {exp}'); bad += 1
fx = pd.read_csv('docs/anchors.no_spa.tsv', sep='\t', keep_default_na=False, na_values=[''])
nw = pd.read_csv(f'{out}/saige_signals.long.tsv', sep='\t', keep_default_na=False, na_values=[''])
m = fx.merge(nw, on=['variant_id', 'cohort'], suffixes=('_fx', '_nw'))
if len(m) != len(fx):
    print(f'    ✗ anchor rows matched {len(m)}/{len(fx)}'); bad += 1
for col in ['P', 'OR', 'L95', 'U95', 'pip', 'case_eaf', 'ctrl_eaf', 'case_hwe_p', 'ctrl_hwe_p', 'tommo_af',
            'kg_jpt_af', 'kg_eas_af', 'gnomad_eas_af', 'gnomad_nfe_af']:
    n = int((~np.isclose(m[col + '_fx'], m[col + '_nw'], rtol=1e-5, atol=1e-9, equal_nan=True)).sum())
    if n:
        print(f'    ✗ {col}: {n} value(s) differ from the anchors'); bad += 1
for col in ['cs_status', 'case_geno', 'ctrl_geno', 'snpEff_Effect', 'role', 'label']:
    n = int((m[col + '_fx'].fillna('') != m[col + '_nw'].fillna('')).sum())
    if n:
        print(f'    ✗ {col}: {n} value(s) differ from the anchors'); bad += 1
print(f'    {len(m)} anchor rows compared')
sys.exit(1 if bad else 0)
PY
ck "no-SPA anchors reproduced" $?

echo "═══ 9. workbook ═══"
OUT="$OUT" python3 - <<'PY'
import json, os, re, sys, warnings, zipfile
from openpyxl import load_workbook
out = os.environ['OUT']
m = json.load(open(f'{out}/_run_info/run_manifest.json'))
with warnings.catch_warnings():
    warnings.simplefilter('ignore', UserWarning)          # the data-bar extension is not read by openpyxl
    wb = load_workbook(f'{out}/saige_signals.xlsx', rich_text=True)
MAIN, OT = 'Common-variant signals', 'Common-variant Open Targets'
ws = wb[MAIN]
bad = 0
want = [MAIN, OT] if m.get('opentargets') else [MAIN]
if wb.sheetnames != want:
    print(f'    ✗ sheets {wb.sheetnames}, expected {want}'); bad += 1
# the XML Excel is strict about (the repair prompt of the openpyxl-written file)
z = zipfile.ZipFile(f'{out}/saige_signals.xlsx')
wbx = z.read('xl/workbook.xml').decode()
t = re.search(r'name="_xlnm.Print_Titles" localSheetId="0">([^<]+)<', wbx)
if not (t and t.group(1).split(',')[0].split('!')[1].startswith('$A')):
    print('    ✗ print titles not columns-first'); bad += 1
for name in z.namelist():
    if name.startswith('xl/worksheets/sheet') and 'inlineStr' in z.read(name).decode():
        print(f'    ✗ inline strings in {name}'); bad += 1
if m.get('opentargets'):
    meta = json.load(open(f'{out}/06.opentargets/opentargets_meta.json'))
    head = ' '.join(str(c.value) for c in ws[4] if c.value is not None)
    if 'GWAS credible sets' not in head or 'Enhancer-to-gene' not in head:
        print('    ✗ Open Targets columns missing'); bad += 1
    wo = wb[OT]                                            # one row per variant, its records below
    from openpyxl.cell.rich_text import CellRichText

    def plain(v):                                          # header text without its note letter (superscript)
        if isinstance(v, CellRichText):
            return ''.join(b if isinstance(b, str) else b.text for b in v
                           if isinstance(b, str) or b.font is None or b.font.vertAlign != 'superscript')
        return '' if v is None else str(v)
    h1 = next((r for r in range(1, 8) if plain(wo.cell(r, 1).value) == 'Locus'), None)
    if h1 is None:
        print(f'    ✗ {OT}: header row not found'); bad += 1; h1 = 3
    head2 = [plain(c.value) for c in wo[h1]]
    QTL_COL = 'cis / trans (to target TSS)'
    for col in ('P / Beta from', 'Set lead variant', 'Study N', QTL_COL, 'Cohorts', 'Population (% of sample)'):
        if col not in head2:
            print(f'    ✗ {OT}: column {col!r} missing'); bad += 1
    if QTL_COL in head2:                                   # every molQTL record classed, and only those
        jc = head2.index(QTL_COL) + 1
        cls = [(r[0], str(r[-1] or '')) for r in wo.iter_rows(min_row=h1 + 1, min_col=5, max_col=jc, values_only=True)
               if r[0] in ('ClinVar', 'GWAS', 'molQTL', 'E2G')]
        n_bad = sum(1 for src, v in cls if bool(v) != (src == 'molQTL'))
        if n_bad:
            print(f'    ✗ {OT}: {n_bad} records with a wrong cis / trans cell'); bad += 1
        print(f"    molQTL cis / trans: {meta.get('qtl_class')}; against Open Targets isTransQtl "
              f"{meta.get('qtl_class_vs_isTransQtl')} ({meta.get('qtl_contradictions')} contradicting)")
    if 'Cohorts' in head2:                                 # cohorts / population on every credible-set record
        jc = max(head2.index('Cohorts'), head2.index('Population (% of sample)')) + 1
        n_bad = sum(1 for r in wo.iter_rows(min_row=h1 + 1, min_col=5, max_col=jc, values_only=True)
                    if r[0] in ('ClinVar', 'GWAS', 'molQTL', 'E2G')
                    and (bool(r[-2]) and bool(r[-1])) != (r[0] in ('GWAS', 'molQTL')))
        if n_bad:
            print(f'    ✗ {OT}: {n_bad} records with a wrong Cohorts / Population cell'); bad += 1
        print(f"    studies (cohorts / population as Open Targets lists them): {meta.get('study_annotation')}")
    last = max(r for r in range(h1 + 1, wo.max_row + 1) if wo.cell(r, 2).value or wo.cell(r, 5).value)
    n_notes = sum(1 for r in range(last + 1, wo.max_row + 1) if wo.cell(r, 1).value)
    if not n_notes:
        print(f'    ✗ {OT}: no notes below the table'); bad += 1
    n_var = sum(1 for (v,) in wo.iter_rows(min_row=h1 + 1, min_col=2, max_col=2, values_only=True)
                if v in ('Lead', 'Alt. lead', 'CS variant'))
    n_det = sum(1 for (v,) in wo.iter_rows(min_row=h1 + 1, min_col=5, max_col=5, values_only=True)
                if v in ('ClinVar', 'GWAS', 'molQTL', 'E2G'))
    if n_var != meta['n_variants']:
        print(f"    ✗ Open Targets sheet: {n_var} variant rows, {meta['n_variants']} variants"); bad += 1
    if n_det != meta['n_records']:
        print(f"    ✗ Open Targets sheet: {n_det} records, meta {meta['n_records']}"); bad += 1
    print(f"    Open Targets {meta['data_version']}: {meta['n_in_index']}/{meta['n_variants']} variants in the "
          f"index; sheet: {n_var} variant rows, {n_det} records, {n_notes} note lines")
rules = {str(ws.cell(3, j).border.bottom.color.rgb)[-6:].upper() for j in range(1, ws.max_column + 1)
         if ws.cell(3, j).border.bottom is not None and ws.cell(3, j).border.bottom.color is not None}
for col in m['cohort_colors']:
    if col.upper() not in rules:
        print(f'    ✗ cohort colour {col} not on the header'); bad += 1
text = ' '.join(str(c.value) for sh in wb.worksheets for row in sh.iter_rows() for c in row if c.value is not None)
for tok in ('cohort_compare', 'lead_crosscohort', 'locus_rep', 'assoc_saige', '_mainland', 'p.value', 'no_spa',
            '.tsv', '.py', '/LARGE'):
    if tok in text:
        print(f'    ✗ internal token in the sheet: {tok}'); bad += 1
if 'Abbreviations:' not in text:
    print('    ✗ no abbreviations footnote'); bad += 1
if '…' in text:
    print('    ✗ truncated text (…) in the workbook'); bad += 1
print(f'    {ws.max_row} rows x {ws.max_column} columns')
sys.exit(1 if bad else 0)
PY
ck "sheets, Excel-strict XML, cohort colours, Open Targets, no internal tokens or truncation, footnotes" $?

if [ $fail -eq 0 ]; then echo "═══ ALL PASSED ═══"; else echo "═══ $fail FAILED ═══"; fi
exit $fail

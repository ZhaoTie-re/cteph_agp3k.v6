#!/bin/bash
# Pre-flight and post-run checks for assoc_gene. All read-only. Run from the
# component root. Every section is explained in docs/VERIFICATION.md.
#
# Sections 1-8, 12, 15 and 2b run on the source tree and must pass before a
# run. Sections 9-11, 13, 14 and 16-20 read results/ and are skipped until it
# exists; after a run they must pass too.
C=/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v6/analysis/assoc_gene
cd "$C"
fail=0
ck() { if [ "$2" -eq 0 ]; then printf "  ✓ %s\n" "$1"; else printf "  ✗ %s\n" "$1"; fail=$((fail+1)); fi; }
REPORT_COHORT=$(grep -oE "params\.ReportCohort\s*=\s*'[^']+'" assoc_gene.nf | grep -oE "'[^']+'$" | tr -d "'")

echo "═══ 1. syntax ═══"
source activate dsl2 2>/dev/null
out=$(nextflow lint assoc_gene.nf 2>&1)
echo "$out" | grep -q 'no errors'; ck "nextflow lint: no errors" $?
n=$(echo "$out" | grep -oE '[0-9]+ warning' | grep -oE '[0-9]+' || echo 0)
ck "nextflow lint: no warnings (got ${n:-0})" $([ "${n:-0}" -eq 0 ] && echo 0 || echo 1)
source activate cteph_geno_pro 2>/dev/null
bad=0; for f in scripts/*.py; do python3 -m py_compile "$f" 2>/dev/null || bad=1; done
ck "every script py_compiles ($(ls scripts/*.py | wc -l))" $bad

echo "═══ 2. site rules ═══"
grep -qE '^\s*(cpus|memory)\s*=' assoc_gene.nf nextflow.config; ck "no cpus/memory directive" $([ $? -eq 1 ] && echo 0 || echo 1)
grep -q 'task\.cpus' assoc_gene.nf; ck "no task.cpus (ext.threads only)" $([ $? -eq 1 ] && echo 0 || echo 1)
grep -vE '^\s*\*|^\s*//' assoc_gene.nf | grep -q 'geneFile'; ck "no --geneFile (--setFile only)" $([ $? -eq 1 ] && echo 0 || echo 1)
bad=0
for m in $(grep -E '^\s*clusterOptions' nextflow.config | grep -oE 'm=[0-9]+M' | grep -oE '[0-9]+' | sort -u); do
  [ $((m % 4571)) -eq 0 ] || { echo "    ✗ m=${m}M"; bad=1; }
done
n=$(grep -cE '^\s*clusterOptions' nextflow.config)
ck "every clusterOptions m= is a multiple of 4571 ($n lines)" $bad
grep -q "errorStrategy = 'finish'" nextflow.config; ck "errorStrategy = 'finish'" $?

echo "═══ 2b. publishDir referencing a process input must be a closure ═══"
python3 - <<'EOF'
import re, sys, pathlib
nf = pathlib.Path('assoc_gene.nf').read_text()
bad = []
for m in re.finditer(r'^\s*publishDir\s+"([^"]*)"', nf, re.M):
    for var in re.findall(r'\$\{([^}]+)\}', m.group(1)):
        if not var.strip().startswith('params.'):
            bad.append((nf[:m.start()].count('\n') + 1, var.strip()))
for line, var in bad:
    print(f'    ✗ assoc_gene.nf:{line}: publishDir without a closure references {var}')
print(f"    {len(re.findall(r'publishDir', nf))} publishDir directives checked")
sys.exit(1 if bad else 0)
EOF
ck "publishDir closures" $?

echo "═══ 3. every process's flags are accepted by its script ═══"
python3 - <<'PY'
import re, subprocess, sys, pathlib
nf = pathlib.Path('assoc_gene.nf').read_text()
bodies = dict(re.findall(r'process\s+(\w+)\s*\{(.*?)\n\}', nf, re.S))
# The map is explicit so a process whose script is not listed is itself a failure.
PAIRS = {'PREP_PHENO_COV': 'prep_pheno_cov.py', 'MAP_CHUNK': 'build_gene_map.py',
         'MERGE_MAP': 'merge_gene_map.py', 'EMIT_TEST_INPUTS': 'emit_test_inputs.py',
         'MERGE_RVTEST': 'merge_rvtest.py', 'SCAN': 'scan.py',
         'SIGNIFICANCE': 'gene_significance.py', 'PLOT_SCAN': 'plot_gene_scan.py',
         'PLOT_ANNOTATION': 'plot_annotation.py', 'ROBUSTNESS': 'robust_genes.py',
         'PLOT_GRID': 'plot_grid.py', 'PLOT_ROBUST_GENES': 'plot_robust_genes.py',
         'PLOT_GENE_DETAIL': 'plot_gene_detail.py', 'CATALOGUE': 'catalogue.py',
         'SCAN_SCALE': 'scan_scale.py'}
bad = 0
for proc, body in bodies.items():
    if 'python3 ${script}' in body and proc not in PAIRS:
        print(f'  ✗ {proc} runs a script that verify.sh does not know'); bad += 1
for proc, sc in sorted(PAIRS.items()):
    p = pathlib.Path('scripts', sc)
    if not p.exists():
        print(f'  ✗ {proc}: {sc} missing'); bad += 1; continue
    m = re.search(r'python3 \$\{script\}(.*?)"""', bodies.get(proc, ''), re.S)
    if not m:
        print(f'  ✗ {proc}: no python3 ${{script}} line'); bad += 1; continue
    passed = set(re.findall(r'(?<![\w-])--([a-z][a-z0-9-]*)', m.group(1)))
    h = subprocess.run([sys.executable, str(p), '--help'], capture_output=True, text=True)
    if h.returncode != 0:
        print(f'  ✗ {proc}: {sc} --help fails'); bad += 1; continue
    unknown = sorted(passed - set(re.findall(r'(?<![\w-])--([a-z][a-z0-9-]*)', h.stdout)))
    if unknown:
        print(f"  ✗ {proc}: {sc} does not accept {' '.join('--' + u for u in unknown)}"); bad += 1
print(f'  {len(PAIRS)} process/script pairs checked')
sys.exit(bad)
PY
ck "flag names" $?

echo "═══ 3b. stageAs-expanded inputs are passed to nargs='+' flags ═══"
python3 - <<'EOF'
import re, subprocess, sys, pathlib
nf = pathlib.Path('assoc_gene.nf').read_text()
bodies = dict(re.findall(r'process\s+(\w+)\s*\{(.*?)\n\}', nf, re.S))
bad, checked = 0, 0
for proc, body in bodies.items():
    m = re.search(r'python3 \$\{script\}(.*?)"""', body, re.S)
    if not m: continue
    sc = {'PREP_PHENO_COV': 'prep_pheno_cov.py', 'MAP_CHUNK': 'build_gene_map.py',
          'MERGE_MAP': 'merge_gene_map.py', 'EMIT_TEST_INPUTS': 'emit_test_inputs.py',
          'MERGE_RVTEST': 'merge_rvtest.py', 'SCAN': 'scan.py',
          'SIGNIFICANCE': 'gene_significance.py', 'PLOT_SCAN': 'plot_gene_scan.py',
          'PLOT_ANNOTATION': 'plot_annotation.py', 'ROBUSTNESS': 'robust_genes.py',
          'PLOT_GRID': 'plot_grid.py', 'PLOT_ROBUST_GENES': 'plot_robust_genes.py',
          'PLOT_GENE_DETAIL': 'plot_gene_detail.py', 'CATALOGUE': 'catalogue.py',
          'SCAN_SCALE': 'scan_scale.py'}.get(proc)
    if not sc: continue
    h = subprocess.run([sys.executable, f'scripts/{sc}', '--help'], capture_output=True, text=True).stdout
    multi_ok = set(re.findall(r'--([a-z][a-z0-9-]*)\s+\S+\s+\[\S+\s+\.\.\.\]', h))
    for flag, var in re.findall(r'--([a-z][a-z0-9-]*)\s+\$\{(\w+)\}', m.group(1)):
        # Only a WILDCARD stageAs expands to many files; a plain rename stages one.
        if not re.search(r'path\(\s*' + var + r'\s*,\s*stageAs:\s*[\'"][^\'"]*\*', body): continue
        checked += 1
        if flag not in multi_ok:
            print(f'    ✗ {proc}: --{flag} is stageAs-expanded but the script takes one value'); bad += 1
print(f'    ({checked} stageAs multi-file flags checked)')
sys.exit(1 if bad else 0)
EOF
ck "multi-file flags" $?

echo "═══ 4. shared and sibling files untouched ═══"
n=$(find ../_shared -newermt '2026-09-03 12:00' -name '*.py' 2>/dev/null | wc -l)
ck "_shared/scripts unchanged ($n modified)" $([ "$n" -eq 0 ] && echo 0 || echo 1)
n=$(find ../assoc_saige ../assoc_plink2 ../assoc_hla -newermt '2026-09-03 12:00' \
     \( -name '*.nf' -o -name '*.py' -o -name '*.config' \) 2>/dev/null | wc -l)
ck "three sibling components unchanged ($n modified)" $([ "$n" -eq 0 ] && echo 0 || echo 1)

echo "═══ 5. plot scripts import the shared style ═══"
bad=0
for f in scripts/plot_*.py; do
  python3 "$f" --help >/dev/null 2>&1 || { echo "    ✗ $(basename $f)"; bad=1; }
done
ck "every plot_*.py loads plot_style/figure_doc ($(ls scripts/plot_*.py | wc -l))" $bad

echo "═══ 6. rvtest output names and the cmcWald block height (a live smoke test) ═══"
# The .nf script block is a shell string verify.sh never executes, so what
# rvtest names its files is invisible to every static section. A 20-second
# real rvtest run pins it: <out>.CMC.assoc / <out>.CMCWald.assoc / <out>.SkatO.assoc,
# the cmcWald block height 1 + n_covar, and the rho column at index 6 of SkatO.
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
PROJ=$(grep -oE "params\.project_dir\s*=\s*'[^']+'" assoc_gene.nf | grep -oE "/[^']+")
CS="$PROJ/tuning.rv/results/narrow_mainland/02.callset_filter/filtered"
P2=$(grep -oE "params\.plink2\s*=\s*'[^']+'" assoc_gene.nf | grep -oE "/[^']+")
RVD=$(grep -oE "params\.rvtestBin\s*=\s*'[^']+'" assoc_gene.nf | grep -oE "/[^']+")
TBX=$(grep -oE "params\.tabix\s*=\s*'[^']+'" assoc_gene.nf | grep -oE "/[^']+")
if [ -f "$CS.bed" ] && [ -x "$RVD/rvtest" ]; then
  ( cd "$TMP"
    awk '$1==22{print $2}' "$CS.bim" | head -60 > v.txt
    $P2 --bfile "$CS" --extract v.txt --export vcf id-paste=iid bgz --out t --threads 2 >/dev/null 2>&1
    $TBX -f -p vcf t.vcf.gz 2>/dev/null
    awk 'NR<=40{printf "%s%s:%s-%s",(NR==1?"SMOKE\t":","),$1,$4,$4} END{print ""}' \
        <(awk '$1==22' "$CS.bim" | head -40) > set.txt
    # sex must VARY: rvtest refuses a constant covariate.
    { echo -e "fid\tiid\tfatid\tmatid\tsex\tpheno1"
      awk '{print $1"\t"$2"\t0\t0\t"(NR%2+1)"\t"(NR%3==0?2:1)}' "$CS.fam"; } > ph.tsv
    { echo -e "fid\tiid\tsex"
      awk '{print $1"\t"$2"\t"(NR%2+1)}' "$CS.fam"; } > cv.tsv
  ) >/dev/null 2>&1
  bad=0
  NCOV_SMOKE=1
  while IFS='|' read -r tag opt sfx wsfx; do
    ( cd "$TMP"; export PATH="$RVD:$PATH"
      rvtest --inVcf t.vcf.gz --pheno ph.tsv --pheno-name pheno1 --covar cv.tsv \
             --covar-name sex --setFile set.txt --out probe.$tag --noweb \
             --impute mean $opt ) >/dev/null 2>&1
    if [ -f "$TMP/probe.$tag.$sfx.assoc" ]; then
      echo "    ✓ $tag -> probe.$tag.$sfx.assoc"
    else
      echo "    ✗ $tag: declared suffix '$sfx' is wrong; rvtest wrote: $(cd $TMP && ls probe.$tag.*.assoc 2>/dev/null | tr '\n' ' ')"
      bad=1
    fi
    if [ "$tag" = "skato" ] && [ -f "$TMP/probe.$tag.$sfx.assoc" ]; then
      head -1 "$TMP/probe.$tag.$sfx.assoc" | awk -F'\t' '{exit !($1=="Range" && $7=="rho" && $8=="Pvalue")}' \
        && echo "    ✓ skato header: Range ... rho Pvalue (rho at column 7)" \
        || { echo "    ✗ skato header is not Range RANGE N_INFORMATIVE NumVar NumPolyVar Q rho Pvalue"; bad=1; }
    fi
    [ -z "$wsfx" ] && continue
    if [ ! -f "$TMP/probe.$tag.$wsfx.assoc" ]; then
      echo "    ✗ $tag: declared wald suffix '$wsfx' but rvtest wrote no such file"; bad=1
    elif ! awk -F'\t' -v n="$NCOV_SMOKE" 'NR>1{c[$1]++}
              END{for(k in c) if(c[k]!=1+n){print k, c[k]; e=1} exit e}' \
              "$TMP/probe.$tag.$wsfx.assoc" >/dev/null; then
      echo "    ✗ $tag: cmcWald rows per gene != 1+n_covar; merge_rvtest's first-row=burden assumption fails"; bad=1
    else
      echo "    ✓ $tag -> probe.$tag.$wsfx.assoc ($((1+NCOV_SMOKE)) rows per gene)"
    fi
  done < <(python3 - <<'EOF'
import re
for m in re.finditer(r"tag:\s*'([a-z]+)'\s*,\s*opt:\s*'([^']+)'\s*,\s*out:\s*'([A-Za-z]+)'"
                     r"\s*,\s*wald:\s*(?:'([A-Za-z]+)'|null)", open('assoc_gene.nf').read()):
    print('|'.join(g or '' for g in m.groups()))
EOF
)
  ck "params.RvtestMethods suffixes match what rvtest writes" $bad
else
  echo "  · skipped (callset or rvtest unreachable)"
fi

echo "═══ 7. documentation present ═══"
for f in README.md docs/METHODS.md docs/OUTPUTS.md docs/STUDY_NOTES.md docs/VERIFICATION.md docs/FIGURES.md; do
  [ -s "$f" ] || { echo "    ✗ missing $f"; fail=$((fail+1)); }
done
echo "  ✓ 6 documents ($(cat README.md docs/*.md 2>/dev/null | wc -l) lines)"

echo "═══ 8. schema atomicity (one column list, several scripts) ═══"
python3 - <<'EOF'
import ast, pathlib, sys
def const(path, name):
    tree = ast.parse(pathlib.Path('scripts', path).read_text())
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign) and any(
                isinstance(t, ast.Name) and t.id == name for t in node.targets):
            return tuple(ast.literal_eval(node.value))
    return None
bad = 0
def same(label, group):
    global bad
    vals = set(group.values())
    if None in vals:
        print(f"    ✗ {label}: unreadable in {[k for k, v in group.items() if v is None]}"); bad = 1
    elif len(vals) != 1:
        print(f"    ✗ {label}: {len(group)} definitions disagree")
        base = next(iter(group.values()))
        for k, v in group.items():
            if v != base: print(f"        {k}: {set(v) ^ set(base) or 'order only'}")
        bad = 1
    else:
        print(f"    ✓ {label}: {len(group)} definitions identical ({len(next(iter(vals)))} columns)")
def subset(label, small, big):
    global bad
    if small is None or big is None:
        print(f"    ✗ {label}: unreadable"); bad = 1
    elif not set(small) <= set(big):
        print(f"    ✗ {label}: reads {sorted(set(small) - set(big))} that the writer does not emit"); bad = 1
    else:
        print(f"    ✓ {label}: reader columns are a subset of the writer's")
same('gene_scan.tsv SCAN_COLUMNS', {'merge_rvtest.py': const('merge_rvtest.py', 'SCAN_COLUMNS'),
                                    'scan.py': const('scan.py', 'SCAN_COLUMNS'),
                                    'robust_genes.py': const('robust_genes.py', 'SCAN_COLUMNS')})
same('scan_qc.tsv QC_COLUMNS', {'scan.py': const('scan.py', 'QC_COLUMNS'),
                                'gene_significance.py': const('gene_significance.py', 'QC_COLUMNS')})
same('denominator DENOM_COLUMNS', {'scan.py': const('scan.py', 'DENOM_COLUMNS'),
                                   'gene_significance.py': const('gene_significance.py', 'DEN_COLUMNS')})
same('evidence STATES', {'robust_genes.py': const('robust_genes.py', 'STATES'),
                         'vocab.py': const('vocab.py', 'STATES')})
same('METHODS', {'merge_rvtest.py': tuple(sorted(ast.literal_eval(
                     next(n.value for n in ast.walk(ast.parse(pathlib.Path('scripts/merge_rvtest.py').read_text()))
                          if isinstance(n, ast.Assign) and any(isinstance(t, ast.Name) and t.id == 'METHOD_TAIL' for t in n.targets))).keys())),
                 'scan.py': tuple(sorted(const('scan.py', 'METHODS'))),
                 'gene_significance.py': tuple(sorted(const('gene_significance.py', 'METHODS'))),
                 'robust_genes.py': tuple(sorted(const('robust_genes.py', 'METHODS'))),
                 'vocab.py': tuple(sorted(const('vocab.py', 'METHODS')))})
ss = const('emit_test_inputs.py', 'SET_STATS_COLUMNS')
subset('set_stats.tsv <- merge_rvtest.py', const('merge_rvtest.py', 'SET_STATS_READ'), ss)
subset('set_stats.tsv <- robust_genes.py', const('robust_genes.py', 'SET_STATS_READ'), ss)
subset('scan_qc.tsv <- robust_genes.py', const('robust_genes.py', 'QC_READ'), const('scan.py', 'QC_COLUMNS'))
sys.exit(bad)
EOF
ck "schemas in sync" $?

echo "═══ 12. vocabulary: the component knows one model ═══"
# Nothing in the pipeline, the scripts or the documents may refer to a design
# that is not this one. verify.sh itself is excluded (it carries the pattern).
hits=$(grep -rniE 'saige|random.?(model|effect|arm)|\bgrm\b|\btracks?\b|both engines|harmonise|identity_check|fit_null' \
        assoc_gene.nf nextflow.config scripts/*.py README.md docs/*.md 2>/dev/null | wc -l)
[ "$hits" -gt 0 ] && grep -rniE 'saige|random.?(model|effect|arm)|\bgrm\b|\btracks?\b|both engines|harmonise|identity_check|fit_null' \
        assoc_gene.nf nextflow.config scripts/*.py README.md docs/*.md 2>/dev/null | head -5 | sed 's/^/    ✗ /'
ck "no banned vocabulary ($hits hits)" $([ "$hits" -eq 0 ] && echo 0 || echo 1)

echo "═══ 15. tier rule fixture ═══"
python3 - <<'EOF'
import sys; sys.path.insert(0, 'scripts')
from robust_genes import tier_of
def mk(**kw): return {c: {'cmc': v[0], 'skato': v[1]} for c, v in kw.items()}
cases = [(mk(n=(1,1), i=(1,0), f=(0,1)), 1), (mk(n=(1,0), i=(1,0), f=(0,1)), 2),
         (mk(n=(1,1), i=(1,1)), 3), (mk(n=(1,0), i=(0,0), f=(0,1)), 3),
         (mk(n=(1,1), i=(0,0), f=(0,0)), 0), (mk(n=(0,0), i=(0,0), f=(0,0)), 0), ({}, 0)]
bad = [(c, tier_of(c, 3), w) for c, w in cases if tier_of(c, 3) != w]
for c, g, w in bad: print(f'    ✗ {c}: got {g}, want {w}')
print(f'    {len(cases) - len(bad)}/{len(cases)} fixture cases')
sys.exit(1 if bad else 0)
EOF
ck "tier_of" $?

# ─── post-run sections ──────────────────────────────────────────────────────
if [ ! -d results ]; then
  echo; echo "═══ (post-run sections 9-11, 13-14, 16-20 skipped: no results/) ═══"
  echo; if [ $fail -eq 0 ]; then echo "═══ ALL PASSED ═══"; else echo "═══ $fail FAILED ═══"; fi; exit $fail
fi

echo "═══ 9. DAG cardinality (trace.txt) ═══"
TR=results/_run_info/trace.txt
if [ -s "$TR" ]; then
  python3 - <<'EOF'
import collections, pathlib, sys
# 3 cohorts x 2 strata x 22 chromosomes x 2 methods.
EXPECT = {'PREP_PHENO_COV': 3, 'SPLIT_BIM': 3, 'MAKE_MINAC_KEEP': 3, 'MAP_CHUNK': 66,
          'MERGE_MAP': 3, 'EMIT_TEST_INPUTS': 6, 'RVTEST': 264, 'MERGE_RVTEST': 12,
          'SCAN': 6, 'SIGNIFICANCE': 12, 'PLOT_SCAN': 6, 'PLOT_ANNOTATION': 3,
          'COLLECT_ALL': 1, 'ROBUSTNESS': 1, 'PLOT_GRID': 1, 'PLOT_ROBUST_GENES': 1,
          'PLOT_GENE_DETAIL': 1, 'CATALOGUE': 2, 'GENE_MODELS': 1, 'SCAN_SCALE': 1,
          'WRITE_RUN_MANIFEST': 1}
lines = pathlib.Path('results/_run_info/trace.txt').read_text().splitlines()
head = lines[0].split('\t'); i = head.index('name') if 'name' in head else 3
j = head.index('status') if 'status' in head else None
seen, failed = collections.Counter(), collections.Counter()
for ln in lines[1:]:
    f = ln.split('\t')
    if len(f) > i:
        seen[f[i].split(' (')[0]] += 1
        if j is not None and f[j] != 'COMPLETED' and f[j] != 'CACHED':
            failed[f[i].split(' (')[0]] += 1
missing = [p for p in EXPECT if p not in seen]
short = {p: (seen[p], EXPECT[p]) for p in EXPECT if p in seen and seen[p] < EXPECT[p]}
extra = [p for p in seen if p not in EXPECT]
print(f'    trace.txt: {sum(seen.values())} tasks, {len(seen)} processes')
if missing: print(f"    ✗ never ran: {', '.join(sorted(missing))}")
if short:
    print('    ✗ fewer tasks than the DAG implies:')
    for p, (g, w) in sorted(short.items()): print(f'        {p}: {g} / {w}')
if extra: print(f"    ✗ unknown process(es): {', '.join(sorted(extra))}")
if failed: print(f"    ✗ non-completed tasks: {dict(failed)}")
sys.exit(1 if (missing or short or extra or failed) else 0)
EOF
  ck "every process ran the expected number of times, all completed" $?
else
  echo "  · skipped (no trace.txt)"
fi

echo "═══ 10. the span-overlap predictor is exact ═══"
python3 - <<'EOF'
import gzip, csv, collections, pathlib, sys
ok = bad = cells = 0
for mp in sorted(pathlib.Path('results').glob('*/01.map/map.*.tsv.gz')):
    coh = mp.parts[1]; strat = mp.name[len('map.'):-len('.tsv.gz')]
    mg = pathlib.Path(f'results/{coh}/03.assoc/{strat}/cmc/merged.tsv')
    if not mg.exists(): continue
    sets = collections.defaultdict(list)
    with gzip.open(mp, 'rt') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            sets[r['set_name']].append((int(r['pos']), len(r['variant_id'].split(':')[2])))
    def excess(v):
        ps = {p for p, _L in v}
        return sum(1 for p, L in v if L > 1 and any(p < q <= p + L - 1 for q in ps))
    cells += 1
    for r in csv.DictReader(mg.open(), delimiter='\t'):
        nm, ne = r['n_var_map'], r['n_var_engine']
        if nm == 'NA' or ne == 'NA': continue
        if int(nm) + excess(sets[r['set_name']]) == int(ne): ok += 1
        else:
            bad += 1
            if bad <= 3: print(f"    ✗ {coh}/{strat} {r['set_name']}: n_var_map={nm} n_var_engine={ne}")
if not cells:
    print('    · skipped (no merged.tsv yet)'); sys.exit(0)
print(f"    {ok:,}/{ok + bad:,} genes exact ({cells} cohort x stratum)")
sys.exit(1 if bad else 0)
EOF
ck "span-overlap over-count predicted exactly" $?

echo "═══ 11. one Bonferroni threshold per cohort x stratum, across methods ═══"
n=$(find results -path '*/05.signals/*' -name 'significance.*.tsv' 2>/dev/null | wc -l)
if [ "$n" -gt 0 ]; then
  python3 - <<'EOF'
import csv, collections, pathlib, sys
by = collections.defaultdict(set)
for f in pathlib.Path('results').glob('*/05.signals/significance.*.tsv'):
    for r in csv.DictReader(f.open(), delimiter='\t'):
        if r['threshold'] == 'bonferroni': by[(r['cohort'], r['stratum'])].add(r['value'])
bad = {k: v for k, v in by.items() if len(v) != 1}
for k, v in sorted(bad.items()): print(f'    ✗ {k[0]}/{k[1]}: {len(v)} Bonferroni values {sorted(v)}')
print(f"    {len(by)} cohort x stratum cells, one Bonferroni value each")
sys.exit(1 if bad else 0)
EOF
  ck "Bonferroni unique across methods" $?
else
  echo "  · skipped (no significance tables)"
fi

echo "═══ 13. set_stats invariants ═══"
python3 - <<'EOF'
import csv, json, pathlib, sys
bad = n = 0
for f in sorted(pathlib.Path('results').glob('*/02.test_inputs/*/set_stats.tsv')):
    coh = f.parts[1]
    coding = json.load(open(f'results/{coh}/00.prep/pheno_coding.json'))
    for r in csv.DictReader(f.open(), delimiter='\t'):
        n += 1
        v = {k: int(r[k]) for k in ('n_var', 'mac', 'mac_case', 'mac_control', 'n_carrier_case', 'n_carrier_control',
                                    'n_case', 'n_control', 'n_multi_case', 'n_multi_control')}
        # A sample may carry several sites of one gene: it is one carrier and
        # n sites in the MAC, so the upper bound is 2 * n_var * carriers.
        checks = [v['mac'] == v['mac_case'] + v['mac_control'],
                  v['n_carrier_case'] <= v['mac_case'] <= 2 * v['n_var'] * v['n_carrier_case'],
                  v['n_carrier_control'] <= v['mac_control'] <= 2 * v['n_var'] * v['n_carrier_control'],
                  v['n_multi_case'] <= v['n_carrier_case'] <= v['n_case'],
                  v['n_multi_control'] <= v['n_carrier_control'] <= v['n_control'],
                  v['n_case'] == coding['n_case'] and v['n_control'] == coding['n_control']]
        if not all(checks):
            bad += 1
            if bad <= 3: print(f"    ✗ {f}: {r['set_name']} {v}")
if not n: print('    · skipped (no set_stats.tsv yet)'); sys.exit(0)
print(f'    {n - bad:,}/{n:,} sets satisfy carriers <= MAC <= 2 n_var carriers, multi <= carriers <= N, N == pheno_coding')
sys.exit(1 if bad else 0)
EOF
ck "set_stats invariants" $?

echo "═══ 14. effect-size arithmetic and rho placement ═══"
python3 - <<'EOF'
import csv, math, pathlib, sys
f = pathlib.Path('results/_comparison/tables/gene_scan_all.tsv')
if not f.exists(): print('    · skipped (no gene_scan_all.tsv yet)'); sys.exit(0)
bad = n = 0
for r in csv.DictReader(f.open(), delimiter='\t'):
    n += 1
    ok = True
    if r['method'] == 'cmc':
        ok &= r['rho'] == 'NA'
        if r['or'] != 'NA':
            b, s, o, lo, hi = (float(r[k]) for k in ('beta', 'se', 'or', 'or_l95', 'or_u95'))
            ok &= abs(o - math.exp(b)) <= 1e-5 * o and lo <= o <= hi
            ok &= abs(lo - math.exp(b - 1.959963984540054 * s)) <= 1e-5 * lo
    else:
        ok &= r['or'] == 'NA' and r['beta'] == 'NA'
        ok &= (r['rho'] != 'NA') or (r['pvalue'] == 'NA')
    if not ok:
        bad += 1
        if bad <= 3: print(f"    ✗ {r['cohort']}/{r['stratum']}/{r['method']} {r['set_name']}: "
                           f"beta={r['beta']} se={r['se']} or={r['or']} rho={r['rho']}")
print(f'    {n - bad:,}/{n:,} rows: OR = exp(beta), CI ordered, rho only on skato rows')
sys.exit(1 if bad else 0)
EOF
ck "effect sizes" $?

echo "═══ 16. summary-table formats ═══"
python3 - "$REPORT_COHORT" <<'EOF'
import csv, pathlib, re, sys
f = pathlib.Path(f'results/_comparison/robust_genes/summary_table.{sys.argv[1]}.tsv')
if not f.exists(): print('    · skipped (no summary table yet)'); sys.exit(0)
RX = {'tier': r'^[0123]$|^—$', 'n_variants': r'^\d+$|^—$',
      'carriers_case': r'^\d+ \(\d+\.\d{2}%\)$|^—$', 'carriers_control': r'^\d+ \(\d+\.\d{2}%\)$|^—$',
      'mac_case': r'^\d+$|^—$', 'mac_control': r'^\d+$|^—$',
      'maf_case': r'^\d\.\d{3}$|^—$', 'maf_control': r'^\d\.\d{3}$|^—$',
      'skato_rho': r'^(0(\.\d+)?|1)$|^—$', 'skato_p': r'^\d\.\d{2}E[-+]\d{2}$|^—$',
      'cmc_or_95ci': r'^\d+\.\d{2} \(\d+\.\d{2}-\d+\.\d{2}\)$|^—$', 'cmc_p': r'^\d\.\d{2}E[-+]\d{2}$|^—$',
      'multi_site_carriers': r'^\d+ / \d+$|^NA$|^—$'}
rows = list(csv.DictReader(f.open(), delimiter='\t'))
bad = 0
for r in rows:
    for col, v in r.items():
        base = col.split('.', 1)[1] if '.' in col else None
        if base in RX and not re.match(RX[base], v):
            bad += 1
            if bad <= 3: print(f"    ✗ {r['gene_symbol']} {col} = {v!r}")
print(f'    {len(rows)} gene rows, {bad} malformed cells')
sys.exit(1 if bad else 0)
EOF
ck "summary table cell formats" $?

echo "═══ 17. evidence completeness ═══"
python3 - <<'EOF'
import csv, collections, pathlib, sys
ev = pathlib.Path('results/_comparison/robust_genes/evidence.tsv')
gi = pathlib.Path('results/_comparison/tables/gene_index_all.tsv')
if not ev.exists() or not gi.exists(): print('    · skipped'); sys.exit(0)
idx = collections.defaultdict(set)
for r in csv.DictReader(gi.open(), delimiter='\t'): idx[(r['cohort'], r['stratum'])].add(r['set_name'])
strata = {s for _c, s in idx}; cohorts = {c for c, _s in idx}
universe = {s: set().union(*(idx[(c, s)] for c in cohorts)) for s in strata}
rows = list(csv.DictReader(ev.open(), delimiter='\t'))
want = sum(len(universe[s]) for s in strata) * len(cohorts) * 2
nim = collections.Counter((r['cohort'], r['stratum']) for r in rows if r['state'] == 'not_in_map')
bad = len(rows) != want
for (c, s), n in sorted(nim.items()):
    exp = 2 * (len(universe[s]) - len(idx[(c, s)]))
    if n != exp: print(f'    ✗ {c}/{s}: {n} not_in_map rows, expected {exp}'); bad = True
print(f'    {len(rows):,} rows (expected {want:,}); not_in_map per cell matches |universe| - |index|')
sys.exit(1 if bad else 0)
EOF
ck "evidence table is the full design" $?

echo "═══ 18. map identity holds everywhere ═══"
python3 - <<'EOF'
import csv, collections, pathlib, sys
f = pathlib.Path('results/_comparison/tables/map_identity_all.tsv')
if not f.exists(): print('    · skipped'); sys.exit(0)
ok = {'ok', 'explained_by_span_overlap', 'explained_by_collision'}
c = collections.Counter(r['reason'] for r in csv.DictReader(f.open(), delimiter='\t'))
bad = {k: v for k, v in c.items() if k not in ok}
print('    ' + ' '.join(f'{k}={v:,}' for k, v in sorted(c.items())))
sys.exit(1 if bad else 0)
EOF
ck "no missing / unexplained set" $?

echo "═══ 19. one README per figure family; every PNG documented once ═══"
python3 - <<'EOF'
import csv, pathlib, sys
root = pathlib.Path('results/figures')
if not root.is_dir():
    print('    · skipped (no results/figures)'); sys.exit(0)
bad = 0
fams = sorted(p for p in root.iterdir() if p.is_dir())
stray = [p for p in pathlib.Path('results').rglob('*.png') if root not in p.parents]
if stray:
    print(f'    ✗ {len(stray)} PNG(s) outside results/figures: {stray[:3]}'); bad = 1
for fam in fams:
    pngs = sorted(fam.glob('*.png'))
    if not (fam / 'README.md').is_file():
        print(f'    ✗ {fam.name}: no README.md'); bad = 1; continue
    extra_md = [m for m in fam.rglob('*.md') if m.name != 'README.md']
    if extra_md:
        print(f'    ✗ {fam.name}: stray sidecar(s) {[m.name for m in extra_md][:3]}'); bad = 1
    txt = (fam / 'README.md').read_text()
    if len(pngs) > 1:
        missing = [p.stem for p in pngs if p.stem not in txt]
        if missing:
            print(f'    ✗ {fam.name}: README does not mention {missing[:3]}'); bad = 1
    idx = fam / 'index.tsv'
    if idx.is_file():
        rows = list(csv.DictReader(idx.open(), delimiter='\t'))
        if len(rows) != len(pngs):
            print(f'    ✗ {fam.name}: index.tsv has {len(rows)} rows, {len(pngs)} PNGs'); bad = 1
        for r in rows:
            if not (fam / r['png']).is_file():
                print(f'    ✗ {fam.name}: index names missing {r["png"]}'); bad = 1
        # a single gene may have no Ensembl 86 model; all of them at once means
        # the gene model panel silently drew nothing (staging / R environment).
        col = 'n_genes_in_model_panel'
        if rows and col in rows[0] and not any(int(r[col] or 0) for r in rows):
            print(f'    ✗ {fam.name}: no figure has a gene model panel'); bad = 1
print(f'    {len(fams)} figure families, {sum(len(list(f.glob("*.png"))) for f in fams)} PNGs')
sys.exit(bad)
EOF
ck "figure families documented" $?

echo "═══ 21. every figure is authored at the journal double-column width ═══"
python3 - <<'EOF'
import pathlib, sys
try:
    from PIL import Image
except ImportError:
    print('    · skipped (no PIL)'); sys.exit(0)
bad = n = 0
for p in sorted(pathlib.Path('results/figures').rglob('*.png')):
    w, h = Image.open(p).size
    n += 1
    w_in, h_in = w / 600.0, h / 600.0
    limit = 11.5 if p.name == 'robust_genes.png' else 1.5 * 7.2
    if abs(w_in - 7.2) > 0.05 or h_in > limit + 2.5:     # +2.5 in for the caption block
        print(f'    ✗ {p}: {w_in:.2f} x {h_in:.2f} in'); bad += 1
print(f'    {n - bad}/{n} PNGs are 7.20 in wide with a plausible height')
sys.exit(1 if bad else 0)
EOF
ck "figure geometry" $?

echo "═══ 20. variant IDs are chr<N>:<POS>:<REF>:<ALT> ═══"
python3 - <<'EOF'
import csv, pathlib, re, sys
rx = re.compile(r'^chr(\d+|X|Y|MT?):(\d+):([ACGT]+):([ACGT]+)$')
bad = n = 0
for f in pathlib.Path('results').glob('*/02.test_inputs/*/variant_stats.tsv'):
    for r in csv.DictReader(f.open(), delimiter='\t'):
        n += 1
        m = rx.match(r['variant_id'])
        if not m or m.group(2) != r['pos'] or m.group(3) != r['ref'] or m.group(4) != r['alt'] or r['ref'] == r['alt']:
            bad += 1
            if bad <= 3: print(f"    ✗ {f}: {r['variant_id']} vs {r['chrom']}:{r['pos']} {r['ref']}>{r['alt']}")
if not n: print('    · skipped'); sys.exit(0)
print(f'    {n - bad:,}/{n:,} variant IDs parse and agree with their own row')
sys.exit(1 if bad else 0)
EOF
ck "variant-ID invariant" $?

echo
if [ $fail -eq 0 ]; then echo "═══ ALL PASSED ═══"; else echo "═══ $fail FAILED ═══"; fi
exit $fail

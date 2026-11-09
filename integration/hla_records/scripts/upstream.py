#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : The upstream guard and the record of every file the workbook reads.
#           Stops the run unless: the association run manifest names the
#           component, its cohorts and comparisons; the newest association verify
#           log ends with 0 failed and is newer than the run manifest; the typing
#           matrices' md5 sums equal those in the typing summary; the typing verify
#           log (run by this pipeline, read-only) has no failed check; both decks
#           and their data files exist. Writes upstream.json: one row per upstream
#           fact for the sheet 'Sources and checks', and per file read its sheets,
#           path, modification time and md5.
# Inputs  : --typing (hla.typing.v2) --assoc (assoc_hla.v2) --typing-verify-log
#           --given57v --gwas-sumstats (one per cohort) --typing-deck --assoc-deck
#           --cram-summary --model-inputs
# Outputs : --out (upstream.json)
# Component: hla_records
# ---------------------------------------------------------------------------
import argparse
import hashlib
import json
import re
from datetime import datetime
from pathlib import Path

COH = ['narrow_mainland', 'intermediate_mainland', 'full_mainland']
CS = ['all', 'dvt_yes', 'dvt_no', 'dvt_case_only']


def md5(p):
    h = hashlib.md5()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(1 << 22), b''):
            h.update(b)
    return h.hexdigest()


def mtime(p):
    return datetime.fromtimestamp(Path(p).stat().st_mtime).strftime('%Y-%m-%d %H:%M')


def die(msg):
    raise SystemExit(f'upstream: {msg}')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    for k in ('typing', 'assoc', 'typing_verify_log', 'given57v', 'typing_deck', 'assoc_deck', 'cram_summary',
              'model_inputs', 'out'):
        ap.add_argument('--' + k.replace('_', '-'), required=True)
    ap.add_argument('--gwas-sumstats', nargs=3, required=True, help='narrow, intermediate, full mainland')
    a = ap.parse_args()
    TY, AS = Path(a.typing), Path(a.assoc)
    runs = []

    # association: run manifest and the newest verify log
    man = json.loads((AS / 'results/_run_info/run_manifest.json').read_text())
    if man.get('component') != 'assoc_hla.v2':
        die(f"run manifest names component {man.get('component')}")
    if [c for c in COH if c not in str(man.get('cohorts'))] or [c for c in CS if c not in str(man.get('case_sets'))]:
        die('run manifest lacks a cohort or a comparison')
    logs = sorted((AS / 'logs').glob('verify.run*.log'), key=lambda p: int(re.search(r'run(\d+)', p.name).group(1)))
    if not logs:
        die('no association verify log')
    vlog = logs[-1]
    m = re.search(r'(\d+) passed, (\d+) failed', vlog.read_text())
    if not m or int(m.group(2)) != 0:
        die(f'{vlog.name}: not a clean verification')
    if vlog.stat().st_mtime < (AS / 'results/_run_info/run_manifest.json').stat().st_mtime:
        die(f'{vlog.name} is older than the run manifest')
    run_no = re.search(r'run(\d+)', vlog.name).group(1)
    runs += [{'component': 'HLA association', 'item': 'Component folder', 'value': str(AS)},
             {'component': 'HLA association', 'item': 'Run', 'value': f'run {run_no}, manifest of {mtime(AS / "results/_run_info/run_manifest.json")}'},
             {'component': 'HLA association', 'item': 'Verification', 'value': f'{vlog.name}: {m.group(1)} passed, 0 failed'},
             {'component': 'HLA association', 'item': 'Cohorts, comparisons', 'value': '3 cohorts × 4 comparisons, as the manifest'},
             {'component': 'HLA association', 'item': 'Supplementary run', 'value': f'{a.given57v} (without DVT, given DQB1 57 V; '
              'the pipeline scripts run unchanged outside the component)'}]

    # typing: the matrices against the md5 sums of the typing summary, and the typing verification (run here)
    summ = json.loads((TY / 'results/06.assoc_inputs/summary.json').read_text())
    for c in COH:
        for f, want in summ['cohorts'][c]['md5'].items():
            got = md5(TY / f'results/06.assoc_inputs/{c}/{f}')
            if got != want:
                die(f'typing matrix {c}/{f}: md5 {got} differs from the typing summary')
    tv = Path(a.typing_verify_log).read_text()
    n_ok, n_fail = len(re.findall(r'^\[ok\]', tv, flags=re.M)), len(re.findall(r'^\[FAIL\]', tv, flags=re.M))
    if n_fail or not n_ok:
        die(f'typing verification: {n_ok} ok, {n_fail} failed')
    runs += [{'component': 'HLA typing', 'item': 'Component folder', 'value': str(TY)},
             {'component': 'HLA typing', 'item': 'Association matrices', 'value': 'md5 of the 6 matrices equal to the typing '
              f'summary of {mtime(TY / "results/06.assoc_inputs/summary.json")}'},
             {'component': 'HLA typing', 'item': 'Verification', 'value': f'its verify.sh, run read-only by this pipeline: '
              f'{n_ok} checks ok, 0 failed'}]

    # decks
    for name, d in (('HLA typing deck', Path(a.typing_deck)), ('HLA association deck', Path(a.assoc_deck))):
        pptx = sorted(d.glob('*.pptx'))
        pptx = [p for p in pptx if p.name in ('hla_typing.pptx', 'hla_assoc.pptx')]
        if len(pptx) != 1:
            die(f'{name}: expected one deck in {d}')
        runs.append({'component': name, 'item': 'Deck', 'value': f'{pptx[0]}, saved {mtime(pptx[0])}, md5 {md5(pptx[0])}'})

    # the files read, per sheet
    R = AS / 'results'
    files = []

    def add(sheets, role, p):
        p = Path(p)
        if not p.exists():
            die(f'missing input {p}')
        files.append({'sheets': sheets, 'role': role, 'path': str(p), 'modified': mtime(p), 'md5': md5(p)})

    T = TY / 'results'
    for f, sh, role in (('00.inputs/sample_inputs.tsv', 'T1', 'per genome: inputs and HLA read pairs'),
                        ('00.inputs/provenance.json', 'T1', 'typing run: versions and settings'),
                        ('01.calls/calls_long.tsv', 'T3', 'raw HLA-HD calls'),
                        ('01.calls/drb_reconciliation.tsv', 'T2', 'DRB3/4/5 reconciliation'),
                        ('02.evidence/call_evidence.tsv', 'T1, T3', 'read evidence and confidence per call'),
                        ('02.evidence/drb_consistency.tsv', 'T2', 'DRB3/4/5 against DRB1'),
                        ('02.evidence/gene_summary.tsv', 'T3', 'per gene summary'),
                        ('02.evidence/sample_qc.tsv', 'T1', 'sample outliers'),
                        ('03.panels/allele_concordance.tsv', 'T4, T6, T7, A10', 'allele frequencies against the panels'),
                        ('03.panels/confidence_vs_panel.tsv', 'T4', 'absent from ToMMo by confidence'),
                        ('03.panels/homozygosity.tsv', 'T4', 'homozygosity and dropout'),
                        ('03.panels/locus_concordance.tsv', 'T4', 'per gene and panel'),
                        ('03.panels/shift_direction.tsv', 'T4', 'shift direction'),
                        ('03.panels/unconfirmed_alleles.tsv', 'T4', 'alleles absent from ToMMo'),
                        ('04.confidence/locus_tiers.tsv', 'T3', 'locus tiers'),
                        ('05.exports/allele_dosage.tsv', 'T7, A7, A10', 'allele copies, every genome'),
                        ('05.exports/residue_dosage.tsv', 'T8', 'residue copies, every genome'),
                        ('05.exports/allele_match_map.tsv', 'T7, A7, A8', 'allele to reference sequence'),
                        ('05.exports/allele_pgroup_map.tsv', 'T6, T7', 'allele to P group'),
                        ('05.exports/residue_sites.tsv', 'T5, T8', 'residue positions'),
                        ('05.exports/residue_reference/DQB1.tsv', 'A7, A8, A10', 'DQB1 residues per allele'),
                        ('06.assoc_inputs/markers.tsv', 'T5, T7, T8', 'every column per cohort, with its status'),
                        ('06.assoc_inputs/qc_record.tsv', 'T5', 'stepwise QC record'),
                        ('06.assoc_inputs/summary.json', 'T5, A1', 'typing summary')):
        add(sh, role, T / f)
    for c in COH:
        for f in ('allele_dosage', 'residue_dosage'):
            add('T5, A7, A10', f'association matrix, {f.split("_")[0]}s, {c.replace("_", " ")}',
                T / f'06.assoc_inputs/{c}/{f}.tsv')
    import sys
    sys.path.insert(0, str(TY / 'scripts'))
    import config as TC                                       # noqa: E402  the typing component's own paths
    for p_, sh, role in ((TC.TOMMO, 'T4, A9, A10', 'ToMMo 61KJPN allele frequencies'),
                         (TC.KG, 'T4, A10', '1000 Genomes HLA types'), (TC.P_GROUPS, 'T6, T7', 'P groups (IPD-IMGT/HLA)'),
                         (TC.IMGT_HISTORY, 'T4, T7', 'IPD-IMGT/HLA allele release history')):
        add(sh, role, p_)
    add('A1', 'run manifest', R / '_run_info/run_manifest.json')
    add('A1', 'comparisons', R / '_reference/case_sets.tsv')
    for f, sh, role in (('steps_all', 'A4, A8', 'stepwise conditioning'), ('signals_all', 'A4, A7', 'signals'),
                        ('lead_given_hla_all', 'A5', 'the GWAS lead given HLA markers'),
                        ('hla_given_lead_all', 'A5, A6', 'HLA markers given the GWAS lead'),
                        ('gwas_leads', 'A5', 'candidate GWAS leads'), ('summary', 'A5', 'conditional summary'),
                        ('ld_lead_hla', 'A5', 'LD of the lead with HLA markers')):
        add(sh, role, R / f'_conditional/{f}.tsv')
    for f, sh, role in (('omnibus_all', 'A3, A12', 'whole positions, stacked'),
                        ('omnibus_steps_all', 'A4', 'position walks'), ('omnibus_signals_all', 'A4', 'position signals')):
        add(sh, role, R / f'_omnibus/{f}.tsv')
    for cs in CS:
        add('A11', f'single markers, {cs}, stacked', R / f'{cs}/_comparison/assoc_all.tsv')
        add('A4, A5, A6', f'positions given a marker, {cs}', R / f'{cs}/_omnibus/omnibus_given.residue.tsv')
        add('A4', f'position walk, {cs}', R / f'{cs}/_omnibus/omnibus_stepwise.residue.tsv')
        add('A12', f'positions, {cs} (recheck)', R / f'{cs}/_omnibus/omnibus.residue.tsv')
        for k in ('allele', 'residue'):
            add('A4', f'stepwise scans, {cs} {k}', R / f'{cs}/_conditional/random/stepwise.{k}.tsv')
        for c in COH:
            for k in ('allele', 'residue'):
                add('A2', f'thresholds, {cs} {c} {k}', R / f'{cs}/{c}/01.markers/{k}/effective_tests.{k}.tsv')
                add('A2', f'exact max-T, {cs} {c} {k}', R / f'{cs}/{c}/01.markers/{k}/maxt_boot.{k}.tsv')
                for mdl in ('random', 'fixed'):
                    add('A11', f'single markers, {cs} {c} {mdl} {k} (recheck)', R / f'{cs}/{c}/03.tables/assoc.{mdl}.{k}.tsv')
            for mdl in ('random', 'fixed'):
                add('A1', f'phenotypes, {cs} {c} {mdl}', R / f'{cs}/{c}/00.prep/pheno_covar.{mdl}.tsv')
            add('A1', f'null model, {cs} {c}', R / f'{cs}/{c}/00.prep/null_model.json')
            add('A1', f'phenotype coding, {cs} {c}', R / f'{cs}/{c}/00.prep/pheno_coding.json')
            if cs in ('all', 'dvt_no'):
                add('A4, A5, A6, A8', f'conditional scans, {cs} {c}', R / f'{cs}/{c}/04.conditional/random/hla_given_signal.tsv')
    d = R / 'all/full_mainland/01.markers'
    for k in ('allele', 'residue'):
        for f in (f'tests.{k}.tsv', f'{k}.bed', f'{k}.bim', f'{k}.fam', f'gec.{k}.block.txt'):
            add('A2', f'GEC recomputation, {k}', d / k / f)
    add('A2', 'GEC recomputation, unrelated samples', R / 'all/full_mainland/00.prep/fixed.keep')
    G = Path(a.given57v)
    add('A8', 'without DVT given DQB1 57 V: positions', G / 'omnibus/omnibus_given.residue.tsv')
    for c in COH:
        add('A8', f'without DVT given DQB1 57 V: markers, {c.replace("_", " ")}', G / c / 'hla_given_signal.tsv')
    for c, p_ in zip(COH, a.gwas_sumstats):
        add('A6', f'genome-wide scan, {c.replace("_", " ")}', p_)
        add('A1', f'GRM markers, {c.replace("_", " ")}', Path(p_).parents[2] / '00.prep/grm.prune.in')
    add('T1, A10', 'sequencing platform per genome', a.cram_summary)
    for c in COH:
        add('A1, A9', f'model inputs (sex, DVT history), {c.replace("_", " ")}', Path(a.model_inputs) / c / 'phenotype/pheno.tsv')
    Path(a.out).write_text(json.dumps({'runs': runs, 'files': files}, indent=1, ensure_ascii=False) + '\n')
    print(f'[upstream] association run {run_no} ({m.group(1)} checks passed); typing {n_ok} checks ok; {len(files)} files')


if __name__ == '__main__':
    main()

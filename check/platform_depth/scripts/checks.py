#!/usr/bin/env python3
"""The check registry of check/platform_depth: every check and gate described once (data only; nothing here
reads results).

  CHECKS     the columns of the evidence matrix (14.evidence): artifact targeted, samples, readout, rule and
             thresholds, tables, figure; evidence_matrix.py takes its rules and thresholds from here, figures.py
             its column labels, verify_checks.py its vocabulary; every check is read in all three cohorts
  GATES      G1-G10: what reproduces what, where the record is
  STAGES     the --run stages: their arguments, run-tag prefix, processes (run.sh reads them)
  PROCESSES  every process of platform_depth.nf: stage, scripts, what it publishes, the gates and checks it feeds

docs/PIPELINE.md is written from here (`python3 scripts/checks.py markdown > docs/PIPELINE.md`); verify_checks.py
fails when it is out of date, when a process of the workflow is missing here, or when the published rules differ.
"""
import sys
from collections import namedtuple

# ── thresholds of the evidence-matrix rules (the rule texts below quote them) ──────────────────────────────────
LOO_Z_EXPECTED = 2.8      # a left-out group is testable from |z expected| = |pooled beta| / SE without it (~80 % power)
LOO_Z_OBSERVED = 1.96     # ... and flags when |z observed| falls below this or changes sign
DEPTH_DZ = 1              # depth: flag when the replicate-median |z| is more than this below the original |z|
Q_ALPHA = 0.05            # heterogeneity: flag below Q_ALPHA / the number of SAIGE signals
Q_MIN_MAF = 0.05          # heterogeneity is read at MAF >= 5 % only
MIN_COHORTS = 1           # a check flags a signal when its rule holds in any one of the 3 cohorts (user, 2026-11-04:
                          # was 2; the cohorts are nested, and a flag only prompts a closer look, so the screen is
                          # sensitive; the number of cohorts flagged is reported with the flag)
MIN_COHORTS_BY_CHECK = {}             # no check-specific exception


def min_cohorts(check):
    """the number of cohorts in which a check's rule must hold for the check to flag a signal"""
    return MIN_COHORTS_BY_CHECK.get(check, MIN_COHORTS)

# ── readings ───────────────────────────────────────────────────────────────────────────────────────────────
READINGS = {
    'consistent': 'the check agrees with the signal (arrays: same direction and |z WGS - z array| <= 3)',
    'flag': 'the check disagrees; flags, never removes',
    'not testable': 'no left-out group has |z expected| >= 2.8',
    'not applicable': 'the check does not apply (Q and random effects at MAF < 5 %; per-stratum estimates of genes)',
    'reported': 'a number without a rule (the residue position; the case carriers of a gene\'s drivers confirmed on '
                'the imputed arrays: one cell for the three cohorts, counted in full mainland, which contains the other two)',
    'not checkable': 'arrays: not typed, not imputed with R2 >= 0.8, no LD proxy with r2 >= 0.8 (genes: no driver in '
                     'the imputed arrays)',
}
DISCOVERY_READINGS = ['genome-wide', 'suggestive', 'below 1e-5', 'significant', 'not significant']
# the per-signal summary (evidence.signals.tsv) would add 'flag in 1/3' below the cohort rule (unused while MIN_COHORTS
# is 1), and lists the three
# cohorts' readings ('a / b / c') when they differ without a flag or a 'consistent' among them
ARRAY_READING = {'consistent on arrays': 'consistent', 'carriers counted': 'reported', 'not checkable': 'not checkable'}
# the gene sets the evidence matrix reads: M/H only, the main analysis (user 2026-11-05: L/M/H is not the main
# analysis; its re-tests stay in the stage tables as data, but no reading, flag or count is made from them)
GENE_SETS_READ = ['moderate_high']
ARRAY_GENE_COHORT = 'full_mainland'   # gene carriers: counted here once, one matrix cell for the three cohorts (user 2026-11-05)

# ── the checks: one per evidence-matrix column, in column order ─────────────────────────────────────────────
Check = namedtuple('Check', 'key title label artifact samples signals readout rule tables figure methods')
CHECKS = [
    Check('discovery', 'Discovery', None,
          'none: the main analysis\' own reading, shown beside the checks',
          'as published, three cohorts',
          'every signal',
          'P of the ORIGINAL arm (= the published result, gates G1-G4; `*_original` columns) against its threshold',
          'the main analysis: P against its pre-specified threshold (GWAS 5e-8, suggestive 1e-5; genes '
          'Bonferroni; HLA residues GEC; positions Bonferroni); a reading, never a flag',
          '07.downsampling/<label>/summary/ds.*.tsv', '', '§1'),
    Check('depth', 'Read down-sampling', 'Read\ndown-\nsampling',
          'read depth: every case platform but DNBseq-G400RS 15x against controls sequenced at 15x',
          'the 408 down-sampled cases at the controls\' median depth (89 NovaSeq and DNBSeq-G400RS 30x cases each '
          'to it; 319 DNBSeq-T7 cases by one platform fraction), R = 3 replicates; three cohorts',
          'every signal',
          'replicate-median |z| against the original |z|; replicates significant',
          'read down-sampling of the down-sampled cases: flag if the replicate-median |z| is more than 1 below the '
          'original |z| (z from the two-sided P)',
          '07.downsampling/<label>/summary/ds.{gwas,hla_dqb1_57v,positions,genes}.tsv',
          '02.downsampling/{gwas_leads,hla_dqb1_57v}', '§3, §5'),
    Check('platform_Q', 'Platform heterogeneity', 'Platform\nhetero-\ngeneity',
          'a call artifact of one case platform',
          'each case platform\'s cases against all controls (shared); three cohorts',
          'GWAS leads, DQB1:57:V; genes not applicable',
          'Q across the four case platforms, its P against null SNPs of the same MAF bin, I²',
          'heterogeneity across the four case platforms (shared controls): flag if the empirical P of Q '
          '< 0.05 / number of SAIGE signals; applied at MAF >= 5% only',
          '10.platform/<cohort>/platform.heterogeneity.tsv', '03.platform/forest.<cohort>', '§6a'),
    Check('platform_RE', 'Random effects across platforms', 'Random\neffects',
          'an effect carried by some case platforms only',
          'as platform_Q',
          'GWAS leads, DQB1:57:V; genes not applicable',
          'random-effects OR (REML, Hartung-Knapp) against the pooled OR',
          'random-effects estimate (REML, Hartung-Knapp) on the other side of the null from the pooled one; applied at '
          'MAF >= 5% only (as Q: SAIGE\'s one-step stratum beta is unreliable in small strata)',
          '10.platform/<cohort>/platform.heterogeneity.tsv', '03.platform/forest.<cohort>', '§6b'),
    Check('lopo', 'Leave one platform out', 'Leave one\nplatform\nout',
          'a signal carried by one case platform',
          'all cases but one platform\'s, all controls; three cohorts',
          'GWAS leads, DQB1:57:V, genes (CMC burden β); the position DQB1 57 reported',
          'z observed without each platform against z expected = pooled β / SE without it',
          'leave one case platform out: flag if a left-out platform has |z expected| >= 2.8 (z expected = pooled '
          'beta / SE without it) and |z observed| < 1.96 or the opposite sign; genes: the CMC burden beta',
          '10.platform/<cohort>/platform.lopo.tsv; 09.platform_assoc/<cohort>/gene.without_*.tsv; '
          '14.evidence/evidence.loo.tsv',
          '06.evidence/z_observed_expected (a-c)', '§6c, §10'),
    Check('array', 'SNP arrays', 'SNP\narrays',
          'a WGS calling error (a second technology in the same people)',
          'the unrelated genomes with arrays, DNBSeq-T7 cases and HiSeq X controls: narrow 1,398, intermediate '
          '1,460, full 1,683; three cohorts',
          'GWAS leads, DQB1:57:V; genes described (one cell)',
          'one marker per signal: the lead typed, else imputed (R2 >= 0.8; DQB1:57:V with HIBAG), else an LD proxy '
          '(r2 >= 0.8); its association re-run in the same people against WGS; check/array_crosscheck re-run per '
          'cohort (gate G10). Genes: the WGS case carriers of each driver confirmed on the imputed arrays (dosage '
          '>= 0.5), counted once in full mainland, which contains the other two cohorts; a description, no rule',
          'SNP arrays (DNBSeq-T7 cases, HiSeq X controls): flag if the opposite direction or |z WGS - z array| > 3, '
          'in any cohort; the array P is not required; genes reported, never flagged',
          '11.array/array.{signals,evidence}.tsv; 11.array/cohorts/<cohort>/', '04.array/array_evidence', '§7'),
    Check('sources_Q', 'Control-source heterogeneity', 'Control\nsource\nhetero-\ngeneity',
          'a control-source batch or population difference',
          'every case against each control source\'s controls (5 sources); three cohorts',
          'GWAS leads, DQB1:57:V; genes not applicable',
          'Q across the five sources (shared cases), P, I²',
          'heterogeneity across the five control sources (shared cases): flag if the analytic (chi-square) P of Q '
          '< 0.05 / number of SAIGE signals; applied at MAF >= 5% only',
          '12.controls/controls.heterogeneity.tsv', '05.controls/sources.<cohort>', '§8'),
    Check('sources_loo', 'Leave one control source out', 'Leave one\nsource\nout',
          'a signal carried by one control source',
          'every case against all sources but one; three cohorts',
          'GWAS leads, DQB1:57:V, genes (CMC burden β); the position DQB1 57 reported',
          'z observed without each source against z expected',
          'leave one control source out: as lopo',
          '12.controls/controls.{loo,genes,positions}.tsv; 14.evidence/evidence.loo.tsv',
          '06.evidence/z_observed_expected (d-f)', '§8, §10'),
]
CHECK = {c.key: c for c in CHECKS}


def summary_reading(readings, check=None):
    """one check's reading over the cohorts (in cohort order): a flag in at least min_cohorts(check), 'flag in 1/3',
    the common reading, 'consistent' when one cohort is consistent and the others are not testable or not
    applicable, or else the readings side by side"""
    nflag = sum(r == 'flag' for r in readings)
    if nflag >= min_cohorts(check):
        return 'flag'
    if nflag == 1:
        return 'flag in 1/3'
    if len(set(readings)) == 1:
        return readings[0]
    if 'consistent' in readings:
        return 'consistent'
    return ' / '.join(readings)
RULES = {c.key: c.rule for c in CHECKS}                       # evidence.rules.tsv, in column order
MATRIX_LABEL = {c.key: c.label for c in CHECKS if c.label}    # the figure's columns (every check but discovery)

# ── gates: nothing is compared before the machinery reproduces the published results ────────────────────────
Gate = namedtuple('Gate', 'key what rule process record')
GATES = [
    Gate('G1', 'SAIGE step 2 on the ORIGINAL genotypes reproduces assoc_saige\'s saige.assoc.tsv',
         'every row present; relative difference <= 1e-6 in every compared column', 'GATE_G1',
         '_gates/G1.<cohort>.json'),
    Gate('G2', 'rvtest on the ORIGINAL genotypes reproduces the published gene P, OR and carriers',
         'no mismatching row', 'RETEST_GENE', '_gates/G2.<cohort>.<set>.json'),
    Gate('G3', 'the hla.typing.v2 wrapper reproduces the six 06.assoc_inputs matrices', 'md5 identical',
         'GATE_G3', '_gates/G3.json'),
    Gate('G4', 'the HLA re-test reproduces assoc_all.tsv and the residue-position omnibus',
         'the same rows; relative difference <= 1e-6',
         'GATE_G4', '_gates/G4.json'),
    Gate('G5', 'the pooled null re-fit reproduces the published null model', 'variance ratio and median z within 1 %',
         'GATE_G5', '_gates/G5.<cohort>.json'),
    Gate('G6', 're-processing fidelity: baseline, re-genotyped original GVCFs and original agree',
         'concordance >= 99.5 %, call rate within 0.5 percentage points', 'GATES_DS', '_gates/<label>/G6.json'),
    Gate('G7', 'the swap leaves every genome but the deep cases\' byte-identical',
         'bedio.py swap aborts otherwise', 'SWAP', '04.genotypes/<label>/<arm>/<cohort>/swap.{gwas,gene}.json'),
    Gate('G8', 'achieved depth: depth after / (f x depth before) over the calling intervals',
         'median 0.95-1.05, every genome 0.85-1.15', 'GATES_DS', '_gates/<label>/G8.json'),
    Gate('G9', 'HLA re-typing at f = 1 reproduces the archived typing', 'every call identical', 'GATES_DS',
         '_gates/<label>/G9.json'),
    Gate('G10', 'the array comparison re-run in full mainland reproduces check/array_crosscheck\'s tables',
         'same rows and text; numbers within 1e-6 relative', 'ARRAY_COHORTS', '11.array/array.gate.json'),
]

# ── stages (--run) and processes ────────────────────────────────────────────────────────────────────────────
Stage = namedtuple('Stage', 'key args tag what')
STAGES = [
    Stage('inputs', '--run inputs', 'inputs',
          'sample table, targets, depth and fractions (part of every run; cached after the first)'),
    Stage('gates', '--run gates --run_label downsampling --Replicates 3', 'gates',
          'the ORIGINAL arm re-tested; gates G1-G4'),
    Stage('pilot', '--run gates,pilot,platform --run_label pilot --Replicates 3', 'pilot',
          'the pilot arms (costs for R) with the platform stage'),
    Stage('downsampling', '--run gates,downsampling,platform --run_label downsampling --Replicates 3', 'full',
          'every down-sampling arm (R = 3; user, 2026-11-03) with the platform stage'),
    Stage('platform', '--run gates,platform --run_label downsampling --Replicates 3', 'full',
          'the platform and control-source subsets'),
    Stage('report', '--run report --run_label downsampling', 'report',
          'summaries, evidence matrix and figures from the published results (never cached but ARRAY_COHORTS)'),
]
STAGE = {s.key: s for s in STAGES}

Process = namedtuple('Process', 'name stage scripts reads publishes feeds')
PROCESSES = [
    Process('SAMPLE_TABLE', 'inputs', 'sample_table.py',
            'info sheet; cohort pheno.tsv, fixed-model .fam; tuning.rv call set; CRAM paths', '00.samples/', ''),
    Process('BUILD_TARGETS', 'inputs', 'build_targets.py',
            'integration/saige_signals variants; assoc_gene tiers; original calling intervals; cohort .bim', '02.targets/',
            ''),
    Process('MAKE_CHUNKS', 'inputs', '(inline)', 'SAMPLE_TABLE', '', ''),
    Process('DEPTH', 'inputs', 'depth_probe.py', 'CRAMs; BUILD_TARGETS intervals', '', ''),
    Process('DEPTH_SUMMARY', 'inputs', 'depth_summary.py', 'DEPTH; SAMPLE_TABLE', '01.depth/', ''),
    Process('SUBSET_ORIGINAL', 'gates', 'bedio.py', 'random-model fileset; tuning.rv call set; BUILD_TARGETS fixed sites',
            '04.genotypes/original/<cohort>/', ''),
    Process('RETEST_GWAS', 'gates', 'bedio.py, saige_collect.py', 'an arm\'s GWAS fileset; assoc_saige null models',
            '06.retest/<arm>/<cohort>/gwas.*.tsv', 'depth'),
    Process('GATE_G1', 'gates', 'saige_collect.py', 'RETEST_GWAS (original, best_guess); assoc_saige saige.assoc.tsv',
            '_gates/', 'G1'),
    Process('RETEST_GENE', 'gates', 'gene_retest.py, bedio.py', 'an arm\'s gene fileset; assoc_gene sets, pheno, scan',
            '06.retest/<arm>/<cohort>/gene.<set>.tsv, _gates/', 'G2, depth'),
    Process('HLA_CHAIN', 'gates', 'hla_v2_wrapper.py', 'hla.typing.v2 00.inputs, or an arm\'s re-typing (HLA_INPUTS)',
            '05.hla/<arm>/', ''),
    Process('GATE_G3', 'gates', 'hla_v2_wrapper.py', 'HLA_CHAIN (original); hla.typing.v2 06.assoc_inputs md5', '_gates/',
            'G3'),
    Process('RETEST_HLA', 'gates', '(assoc_hla.v2 scripts)', 'an arm\'s matrices; assoc_hla.v2 00.prep, 01.markers',
            '06.retest/<arm>/<cohort>/assoc.random.<class>.tsv', 'depth'),
    Process('RETEST_OMNIBUS', 'gates', '(assoc_hla.v2 omnibus.py)', 'RETEST_HLA residues of the three cohorts',
            '06.retest/<arm>/_omnibus/', 'depth'),
    Process('GATE_G4', 'gates', 'hla_gate.py', 'RETEST_HLA, RETEST_OMNIBUS (original); assoc_hla.v2 published results',
            '_gates/', 'G4'),
    Process('GATE_SUMMARY', 'gates', '(inline)', 'the G1-G4 records', '_gates/gates_G1-G4.tsv', 'G1-G4'),
    Process('T7_FRACTIONS', 'downsampling', 't7_fractions.py', 'DEPTH_SUMMARY fractions and per-sample depth',
            '01.depth/fractions.all.tsv', ''),
    Process('ARM_PLAN', 'downsampling', 'arm_plan.py', 'T7_FRACTIONS (every down-sampled case); SAMPLE_TABLE',
            '03.calls/<label>/', ''),
    Process('EXTRACT_READS', 'downsampling', '(samtools)', 'CRAMs; BUILD_TARGETS extraction regions', '', ''),
    Process('SUBSAMPLE', 'downsampling', '(samtools)', 'EXTRACT_READS; the arm\'s fraction and seed', '', 'G8'),
    Process('CALL_GVCF', 'downsampling', '(GATK HaplotypeCaller)', 'SUBSAMPLE; calling intervals', '', ''),
    Process('JOINT', 'downsampling', '(GATK GenomicsDBImport, GenotypeGVCFs)', 'every arm\'s GVCF; the original GVCFs',
            '', ''),
    Process('GENOTYPE_QC', 'downsampling', 'genotype_qc.sh, expand_fixed.py', 'JOINT; fixed sites',
            '03.calls/<label>/per_chrom/', ''),
    Process('CONCAT_QC', 'downsampling', '(bcftools)', 'GENOTYPE_QC', '03.calls/<label>/arms.qc.vcf.gz', ''),
    Process('SWAP', 'downsampling', 'bedio.py', 'CONCAT_QC; SUBSET_ORIGINAL; ARM_PLAN deep cases',
            '04.genotypes/<label>/<arm>/<cohort>/swap.*.json', 'G7'),
    Process('HLA_SUBSAMPLE', 'downsampling', '(samtools)', 'hla.typing 01.reads; the arm\'s fraction and seed', '', ''),
    Process('HLAHD', 'downsampling', '(HLA-HD 1.7.1)', 'HLA_SUBSAMPLE', '', ''),
    Process('HLA_INPUTS', 'downsampling', 'hla_v2_wrapper.py', 'HLAHD, HLA_SUBSAMPLE of one arm', '', ''),
    Process('HLA_MATRICES', 'downsampling', 'hla_matrices.py', 'HLA_CHAIN exports of an arm',
            '05.hla/<label>/<arm>/matrices/', ''),
    Process('CONCORDANCE', 'downsampling', 'concordance.py', 'CONCAT_QC; SUBSET_ORIGINAL (full mainland); ARM_PLAN',
            '07.downsampling/<label>/concordance/', 'G6'),
    Process('QC_REASONS', 'downsampling', 'concordance.py', 'CONCAT_QC; the cohort\'s pre-QC VCFs; ARM_PLAN',
            '07.downsampling/<label>/qc_reasons/', ''),
    Process('HLA_CONCORDANCE', 'downsampling', 'hla_concordance.py', 'HLA_CHAIN evidence; HLA_MATRICES; the archive',
            '07.downsampling/<label>/hla/', 'G9'),
    Process('GATES_DS', 'downsampling', 'gates_ds.py', 'CONCORDANCE; SUBSAMPLE depth; HLA_CONCORDANCE', '_gates/<label>/',
            'G6, G8, G9'),
    Process('MAKE_SUBSETS', 'platform', 'make_subsets.py', 'SAMPLE_TABLE; cohort pheno, covariates; assoc_hla.v2 prep',
            '08.subsets/<cohort>/subsets/', ''),
    Process('MAKE_CTRL_SUBSETS', 'platform', 'control_subsets.py', 'MAKE_SUBSETS (pooled); info sheet Project',
            '08.subsets/<cohort>/control_sources/', ''),
    Process('FIT_NULL', 'platform', '(SAIGE step 1)', 'a subset; random-model genotypes; assoc_saige GRM markers',
            '08.subsets/<cohort>/null/<subset>/', 'G5'),
    Process('PLATFORM_GENO', 'platform', 'bedio.py', 'random-model fileset; fixed sites; null SNPs', '', ''),
    Process('STEP2_SUBSET', 'platform', 'bedio.py, saige_collect.py', 'FIT_NULL; PLATFORM_GENO; RETEST_HLA (original)',
            '09.platform_assoc/<cohort>/<subset>.saige.tsv', 'G5'),
    Process('PLATFORM_GENE', 'platform', 'gene_retest.py, bedio.py', 'subset include lists; SUBSET_ORIGINAL gene fileset',
            '09.platform_assoc/<cohort>/gene.<subset>.<set>.tsv', 'lopo, sources_loo (genes)'),
    Process('CARRIER_FIRTH', 'platform', 'gene_retest.py, carrier_firth.py', 'SUBSET_ORIGINAL gene fileset; MAKE_SUBSETS',
            '09.platform_assoc/<cohort>/carriers_firth.tsv', ''),
    Process('PLATFORM_OMNIBUS', 'platform', '(assoc_hla.v2 omnibus.py)', 'RETEST_HLA (original) residues; a subset',
            '09.platform_assoc/_omnibus/', 'lopo, sources_loo (position)'),
    Process('PLATFORM_STATS', 'platform', 'platform_stats.py', 'STEP2_SUBSET (platform subsets); PLATFORM_GENO; '
            'CARRIER_FIRTH; null SNPs', '10.platform/<cohort>/', 'platform_Q, platform_RE, lopo'),
    Process('GATE_G5', 'platform', '(inline)', 'FIT_NULL, STEP2_SUBSET (pooled); RETEST_GWAS (original)', '_gates/', 'G5'),
    Process('DS_SUMMARY', 'report', 'ds_summary.py', '06.retest; 04.genotypes/original/<cohort>; 03.calls deep cases; '
            'qc_reasons', '07.downsampling/<label>/summary/', 'discovery, depth'),
    Process('PLATFORM_REPORT', 'report', 'platform_report.py', '10.platform; 08.subsets; 09.platform_assoc',
            '10.platform/_summary/', ''),
    Process('ARRAY_COHORTS', 'report', 'array_cohorts.py, array_hibag.py, array_hibag.R', 'check/array_crosscheck '
            'inputs and results; each cohort\'s fixed-model members, phenotype, covariates; the QC\'d and imputed '
            'arrays; the HIBAG model', '11.array/cohorts/', 'G10'),
    Process('ARRAY_EVIDENCE', 'report', 'array_evidence.py', '11.array/cohorts (ARRAY_COHORTS); integration/'
            'gene_signals LOVO (the gene drivers)', '11.array/array.*', 'array'),
    Process('CONTROL_REPORT', 'report', 'control_report.py', '08.subsets; 09.platform_assoc', '12.controls/controls.*',
            'sources_Q, sources_loo'),
    Process('HLA_REF_SHIFT', 'report', 'hla_reference_shift.py', 'hla.typing.v2 05.exports; assoc_hla.v2 conditional '
            'scan; 06.retest/original', '12.controls/hla_reference_shift.*', ''),
    Process('FOLLOWUP', 'report', 'signal_followup.py', '04.genotypes/original/<cohort>; the cohort pre-QC VCFs; '
            '05.hla; integration frequencies', '13.followup/', ''),
    Process('EVIDENCE', 'report', 'evidence_matrix.py', 'the tables above; 10.platform; 09.platform_assoc', '14.evidence/',
            'every check'),
    Process('FIGURES', 'report', 'figures.py, catalogue.py', 'every report table; 01.depth; 07.downsampling; 10.platform',
            'figures/', ''),
]
# the stage a process belongs to in the workflow; the 'pilot' and 'downsampling' stages share the arm processes
STAGE_PROCESSES = {
    'inputs': [p.name for p in PROCESSES if p.stage == 'inputs'],
    'gates': [p.name for p in PROCESSES if p.stage in ('inputs', 'gates')],
    'pilot': [p.name for p in PROCESSES if p.stage != 'report'],
    'downsampling': [p.name for p in PROCESSES if p.stage != 'report'],
    'platform': [p.name for p in PROCESSES if p.stage in ('inputs', 'gates', 'platform')],
    'report': [p.name for p in PROCESSES if p.stage in ('inputs', 'report')],
}
REPORT_PROCESSES = [p.name for p in PROCESSES if p.stage == 'report']   # `cache false` in the workflow, except:
REPORT_CACHED = ['ARRAY_COHORTS']     # report stage but cached: staged scripts, frozen inputs (user 2026-11-05)


def _cell(s):
    return str(s).replace('\n', ' ').replace('|', '\\|')


def _code(s, sep=', '):
    """each path or name of a list as a code span ('—' when empty; parenthesised tools stay plain)"""
    def one(x):
        path, _, note = x.partition(' (')
        return x if x.startswith('(') else f'`{path}`' + (f' ({note}' if note else '')
    return sep.join(one(x) for x in str(s).split(sep) if x) or '—'


def markdown():
    """docs/PIPELINE.md"""
    L = ['# Pipeline — `check/platform_depth`', '',
         'Written by `python3 scripts/checks.py markdown > docs/PIPELINE.md` from the registry in `scripts/checks.py` '
         '(do not edit by hand; `verify.sh` fails when it is out of date). Methods: [METHODS.md](METHODS.md); every '
         'file: [OUTPUTS.md](OUTPUTS.md). Paths are under `results/`; `<arm>` is `original` or `<label>/<arm>`.', '',
         '## Stages', '',
         'Run with `bash run.sh <stage>` (preflight, next run tag, detached launch, cache-miss guard; README).', '',
         '| Stage | Arguments | Run tag | What |', '|---|---|---|---|']
    L += [f'| `{s.key}` | `{s.args}` | `{s.tag}.N` | {_cell(s.what)} |' for s in STAGES]
    L += ['| `verify` | `bash verify.sh downsampling` | | read-only checks of every published result |', '']
    for st in ('inputs', 'gates', 'downsampling', 'platform', 'report'):
        L += [f'## Processes: {st}', '']
        if st == 'downsampling':
            L += ['The `pilot` and `downsampling` stages run these with their own `--run_label`; the gates-stage '
                  're-tests (RETEST_*, HLA_CHAIN) run on every arm.', '']
        if st == 'report':
            L += ['`cache false`: every report run re-runs them from the published results; editing their scripts '
                  'never touches the cache.', '']
        L += ['| Process | Scripts | Reads | Publishes | Feeds (gate or matrix column) |', '|---|---|---|---|---|']
        L += [f'| `{p.name}` | {_code(p.scripts)} | {_cell(p.reads)} | {_code(p.publishes)} | {_cell(p.feeds) or "—"} |'
              for p in PROCESSES if p.stage == st]
        L += ['']
    L += ['`choose_R.py` is run by hand after the pilot (`_pilot/R_decision.json`).', '',
          '## Checks (evidence-matrix columns)', '',
          'Readings: ' + '; '.join(f'**{k}** {v}' for k, v in READINGS.items()) + '. Every check is read in each '
          f'cohort and flags a signal when its rule holds in any one of the 3 (was two of three before 2026-11-04). Thresholds: leave-one-out |z expected| ≥ '
          f'{LOO_Z_EXPECTED}, |z observed| < {LOO_Z_OBSERVED}; depth Δ|z| > {DEPTH_DZ}; heterogeneity '
          f'{Q_ALPHA} / signals at MAF ≥ {Q_MIN_MAF}.', '']
    for c in CHECKS:
        L += [f'### `{c.key}` — {c.title}', '',
              f'- **Artifact targeted:** {c.artifact}', f'- **Samples:** {c.samples}', f'- **Signals:** {c.signals}',
              f'- **Readout:** {c.readout}', f'- **Rule:** {c.rule}', f'- **Tables:** {_code(c.tables, "; ")}',
              f'- **Figure:** {_code(c.figure)}', f'- **METHODS:** {c.methods}', '']
    L += ['## Gates', '', '| Gate | What | Rule | Process | Record |', '|---|---|---|---|---|']
    L += [f'| {g.key} | {_cell(g.what)} | {_cell(g.rule)} | `{g.process}` | `{g.record}` |' for g in GATES]
    return '\n'.join(L) + '\n'


def main():
    cmd = sys.argv[1:] or ['markdown']
    if cmd[0] == 'markdown':
        sys.stdout.write(markdown())
    elif cmd[0] == 'stage' and len(cmd) == 2 and cmd[1] in STAGE:
        # run.sh: arguments, run-tag prefix, and the processes that may be submitted (all others must be cached)
        s = STAGE[cmd[1]]
        allowed = REPORT_PROCESSES if s.key == 'report' else \
            [p for p in STAGE_PROCESSES[s.key] if p not in STAGE_PROCESSES['inputs']] if s.key != 'inputs' else \
            STAGE_PROCESSES['inputs']
        print(s.args)
        print(s.tag)
        print('|'.join(allowed))
    else:
        sys.exit('usage: checks.py markdown | stage <' + '|'.join(STAGE) + '>')


if __name__ == '__main__':
    main()

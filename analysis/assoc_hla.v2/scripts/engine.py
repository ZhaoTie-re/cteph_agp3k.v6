#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Purpose : One association run on the pooled PLINK set (build_pooled.py), with
#           the SAME engine as the marginal scan -- plink2 --glm (--condition-
#           list) for the fixed model, SAIGE step 2 (--condition) on the case
#           set's checked null model for the random model -- optionally
#           conditioned on a set of markers. Used by stepwise.py and
#           conditional.py.
#
#           P is the model's own: plink2 Wald P (fixed); SAIGE p.value.NA_c, the
#           unadjusted conditional score P (random), with SE = 1/sqrt(var_c).
#           A marker that the conditioning set (nearly) determines has no
#           conditional P: plink2's error code, or COLLINEAR when SAIGE's
#           conditional variance falls below 1e-6 of its marginal variance. The
#           conditioning markers themselves are returned as CONDITIONED.
#
#           GATES. Every run must analyse the case set's cases and controls.
#           Every SAIGE run with --condition must reproduce, marker by marker,
#           the marginal statistics of the unconditional run (the guard on the
#           SAIGE first-marker defect the sentinel absorbs), so the first run of
#           an Engine is always the unconditional one.
# Component: assoc_hla.v2
# ---------------------------------------------------------------------------
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

Z975 = stats.norm.ppf(0.975)


class Engine:
    """cfg: model, bfile, keep, pheno_fixed, covar_name, rda, vr, plink2,
    saige_rscript, saige_step2, threads. order: marker -> file order.
    want: {'cases': n, 'controls': n} of the model's sample set."""

    def __init__(self, cfg, order, want, prefix='run'):
        self.c, self.order, self.want, self.prefix = cfg, order, want, prefix
        self.n_runs = 0
        self.marginal = None                     # SAIGE: the unconditional run's raw rows

    def run(self, cond):
        cond = sorted(set(cond), key=lambda m: self.order[m])      # SAIGE wants file order
        if self.c.model == 'random' and cond and self.marginal is None:
            sys.exit('ERROR: the first SAIGE run of an engine must be unconditional (marginal gate)')
        self.n_runs += 1
        tag = f'{self.prefix}{self.n_runs:03d}'
        res = self._fixed(cond, tag) if self.c.model == 'fixed' else self._random(cond, tag)
        res = res[res['marker'] != 'SENTINEL'].copy()
        absent = [m for m in cond if m not in set(res['marker'])]
        if absent:                               # the fixed engine skips its covariates
            res = pd.concat([res, pd.DataFrame({'marker': absent})], ignore_index=True)
        res.loc[res['marker'].isin(cond), ['p', 'or', 'beta', 'se', 'l95', 'u95']] = np.nan
        res.loc[res['marker'].isin(cond), 'errcode'] = 'CONDITIONED'
        return res.set_index('marker')

    def _fixed(self, cond, tag):
        c = self.c
        cmd = [c.plink2, '--bfile', c.bfile, '--keep', c.keep, '--pheno', c.pheno_fixed,
               '--pheno-name', 'PHENO1', '--1', '--covar', c.pheno_fixed, '--covar-name', c.covar_name,
               '--glm', 'omit-ref', 'no-firth', 'hide-covar', '--ci', '0.95',
               '--threads', str(c.threads), '--out', tag]
        if cond:
            Path(f'{tag}.cond').write_text('\n'.join(cond) + '\n')
            cmd[cmd.index('--glm'):cmd.index('--glm')] = ['--condition-list', f'{tag}.cond']
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)
        log = Path(f'{tag}.log').read_text()
        mm = re.search(r'(\d+) cases and (\d+) controls remaining after main filters', log)
        if not mm or (int(mm.group(1)), int(mm.group(2))) != (self.want['cases'], self.want['controls']):
            sys.exit(f'ERROR: {tag}: plink2 analysed {mm.groups() if mm else "?"}, expected '
                     f"{self.want['cases']} / {self.want['controls']}")
        g = pd.read_csv(f'{tag}.PHENO1.glm.logistic', sep='\t', dtype={'ID': str})
        if not (g['A1'] == g['ALT']).all():
            sys.exit(f'ERROR: {tag}: plink2 effect allele is not the ALT allele on every row')
        return pd.DataFrame({'marker': g['ID'], 'p': g['P'], 'or': g['OR'], 'beta': np.log(g['OR']),
                             'se': g['LOG(OR)_SE'], 'l95': g['L95'], 'u95': g['U95'],
                             'errcode': g['ERRCODE'].astype(str)})

    def _random(self, cond, tag):
        c = self.c
        cmd = [c.saige_rscript, c.saige_step2, f'--bedFile={c.bfile}.bed', f'--bimFile={c.bfile}.bim',
               f'--famFile={c.bfile}.fam', '--AlleleOrder=alt-first', '--chrom=6',
               f'--GMMATmodelFile={c.rda}', f'--varianceRatioFile={c.vr}', f'--SAIGEOutputFile={tag}.txt',
               '--is_Firth_beta=FALSE', '--LOCO=TRUE', '--is_output_moreDetails=TRUE',
               '--minMAC=0.5', '--maxMissing=0.15', '--impute_method=mean']
        if cond:
            cmd.append('--condition=' + ','.join(cond))
        with open(f'{tag}.log', 'w') as fh:
            subprocess.run(cmd, check=True, stdout=fh, stderr=subprocess.STDOUT)
        s = pd.read_csv(f'{tag}.txt', sep=r'\s+', dtype={'MarkerID': str})
        if set(s['N_case']) != {self.want['cases']} or set(s['N_ctrl']) != {self.want['controls']}:
            sys.exit(f'ERROR: {tag}: SAIGE N_case / N_ctrl differ from the case set')
        real = s[s['MarkerID'] != 'SENTINEL'].set_index('MarkerID')
        if not cond:
            self.marginal = real[['AC_Allele2', 'p.value.NA']].copy()
            beta, var, p = real['BETA'], real['var'], real['p.value.NA']
            colli = pd.Series(False, index=real.index)
        else:
            m = self.marginal.loc[real.index]
            if not (np.allclose(real['AC_Allele2'], m['AC_Allele2'], rtol=1e-6)
                    and np.allclose(real['p.value.NA'], m['p.value.NA'], rtol=1e-6)):
                sys.exit(f'ERROR: {tag}: a marker\'s marginal statistics changed under --condition '
                         '(SAIGE first-marker defect reached a real marker)')
            beta, var, p = real['BETA_c'], real['var_c'], real['p.value.NA_c']
            colli = ~(var > 1e-6 * real['var'])
        with np.errstate(over='ignore', invalid='ignore', divide='ignore'):
            se = 1.0 / np.sqrt(var.where(var > 0))
            ci_lo, ci_hi, orr = np.exp(beta - Z975 * se), np.exp(beta + Z975 * se), np.exp(beta)
        return pd.DataFrame({'marker': real.index, 'p': p.where(~colli).values,
                             'or': orr.where(~colli).values, 'beta': beta.where(~colli).values,
                             'se': se.where(~colli).values,
                             'l95': ci_lo.where(~colli).values, 'u95': ci_hi.where(~colli).values,
                             'errcode': np.where(colli.values, 'COLLINEAR', '.')})

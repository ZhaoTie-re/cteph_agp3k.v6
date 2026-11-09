# GEC v0.2 (Genetic type I Error Calculator), archived copy

The original standalone GEC: Li M-X, Yeung JMY, Cherny SS, Sham PC. Evaluating the effective numbers of
independent tests and significant p-value thresholds in commercial genotyping arrays and public imputation
reference datasets. Hum Genet 2012;131:747-756. doi:10.1007/s00439-011-1118-2.

| file | from | sha256 |
|---|---|---|
| `gecV0.2.zip` | http://pmglab.top/gec/data/archive/v0.2/gecV0.2.zip (the site's "Archive" table, version 0.2) | a3cad06a78098dbc51637218b32f8534b4f247f1c9335e29e9d03524d31529a3 |
| `gec.jar` | `gec/gec.jar` inside the zip (banner: "GEC! v0.2 03/Sep./2015"; jar dated 2017-01-16) | 9d6e60013fb6e6512a7da33ae99faae8daa363fc1d0491b479aaf0b1a109e456 |
| `UserManualV0.2.pdf` | http://pmglab.top/gec/data/archive/v0.2/UserManualV0.2.pdf | (manual only) |

Downloaded 2026-10-23. The site's current release (v1.0) is the GEC function of KGGSEE (update of
2021-02-02); this component uses the original standalone v0.2, as the HLA fine-mapping precedent did
(Dallmann-Sauer et al. 2020, PLoS Pathog). `scripts/effective_tests.py` refuses to run if `gec.jar`'s
sha256 differs from the value above.

Run here (Java 8, `/usr/bin/java`):
`java -Xmx4g -jar gec.jar --effect-number --plink-binary <prefix> --genome --maf 0 --no-web --out <prefix>`

- `--maf 0`: v0.2 drops markers with MAF < 0.01 by default; every marker tested here has MAF >= 0.005 and
  must count.
- `--no-web`: no version check over the network.
- Input alleles must be nucleotides; the 'present' allele P is written as T (absent stays A).
- The manual asks for unrelated individuals; the input is the case set's fixed-model (relatives removed) samples.
- v0.2 estimates Me within blocks of at most 100 consecutive markers and sums over blocks (not settable).

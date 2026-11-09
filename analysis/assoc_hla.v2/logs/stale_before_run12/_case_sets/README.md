# Cases with and without DVT

**Covers:** `dvt_contrast.png`; tables `dvt_contrast.tsv`, `dvt_hits.tsv`, `scan_qc_case_sets.tsv`. The case-only scans themselves: `../dvt_case_only/`

## The question this figure answers

Does any HLA allele or residue associate differently with CTEPH in cases with a history of DVT than in cases without one, each compared with the same controls?

## Panels

**(dvt_contrast.png a–f) rows alleles, residues; columns narrow, intermediate, full mainland**

Each point is one marker tested in the cohort: x = log OR in cases without DVT against all controls, y = log OR in cases with DVT against the same controls, both from the random model (SAIGE, the primary model). Red / blue: the marker is significant at the LD-adjusted (GEC) threshold in that DVT set and cohort. A black ring marks a marker GEC-significant in the case-only test (cases with DVT coded 1 against cases without DVT coded 0, no controls), the direct test of a difference. Grey line: identity; points far from it differ between the sets.

## Interpretation

r(log OR) between the two sets is 0.13–0.24 across the six panels. 54 marker(s) are significant at some threshold in a DVT set or in the case-only test in at least one cohort or model, 0 of them in the case-only test (dvt_hits.tsv). The DVT sets are small (about 140–160 cases each), so each has less power than all cases together, and an apparent gap between them is expected by chance for some markers. Because the two sets share their controls, the gap itself is not a test; the case-only scans (about 290–300 cases) are. They compare cases with cases, so the case–control platform split does not enter them; the platform mix nevertheless differs between cases with and without DVT (docs/STUDY_NOTES.md), and platform is not a covariate.

## Values in this rendering

| quantity | value |
|---|---|
| allele narrow: markers | 143 |
| allele narrow: r(log OR) | 0.132 |
| allele narrow: GEC-significant with DVT / without DVT | 0 / 1 |
| allele narrow: GEC-significant in the case-only test | 0 |
| allele intermediate: markers | 144 |
| allele intermediate: r(log OR) | 0.149 |
| allele intermediate: GEC-significant with DVT / without DVT | 0 / 0 |
| allele intermediate: GEC-significant in the case-only test | 0 |
| allele full: markers | 142 |
| allele full: r(log OR) | 0.239 |
| allele full: GEC-significant with DVT / without DVT | 0 / 0 |
| allele full: GEC-significant in the case-only test | 0 |
| residue narrow: markers | 672 |
| residue narrow: r(log OR) | 0.209 |
| residue narrow: GEC-significant with DVT / without DVT | 0 / 51 |
| residue narrow: GEC-significant in the case-only test | 0 |
| residue intermediate: markers | 681 |
| residue intermediate: r(log OR) | 0.198 |
| residue intermediate: GEC-significant with DVT / without DVT | 0 / 0 |
| residue intermediate: GEC-significant in the case-only test | 0 |
| residue full: markers | 679 |
| residue full: r(log OR) | 0.181 |
| residue full: GEC-significant with DVT / without DVT | 0 / 0 |
| residue full: GEC-significant in the case-only test | 0 |

## Full statistics

**Scan summaries of the DVT sets and the case-only test**

| case_set | class | cohort | model | cases | controls | markers_tested | tests | p_gec | tests_gec | tests_maxt | tests_bonf | lambda_gc | min_p | top_marker | top_or |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| dvt_yes | allele | narrow | fixed | 150 | 1,749 | 143 | 143 | 4.374e-04 | 0 | 0 | 0 | 1.281 | 0.001568 | C*03:04 | 1.738 |
| dvt_yes | allele | narrow | random | 150 | 1,776 | 143 | 143 | 4.374e-04 | 0 | 0 | 0 | 1.464 | 0.001569 | C*03:04 | 1.913 |
| dvt_yes | allele | intermediate | fixed | 155 | 2,051 | 144 | 144 | 4.337e-04 | 0 | 0 | 0 | 1.455 | 0.003107 | C*03:04 | 1.636 |
| dvt_yes | allele | intermediate | random | 155 | 2,079 | 144 | 144 | 4.337e-04 | 0 | 0 | 0 | 1.699 | 0.003032 | C*03:04 | 1.78 |
| dvt_yes | allele | full | fixed | 159 | 2,632 | 142 | 142 | 4.394e-04 | 0 | 0 | 0 | 1.34 | 0.004307 | DQB1*05:02 | 2.231 |
| dvt_yes | allele | full | random | 159 | 2,662 | 142 | 142 | 4.394e-04 | 0 | 0 | 0 | 1.414 | 0.004392 | DQB1*05:02 | 3.003 |
| dvt_yes | residue | narrow | fixed | 150 | 1,749 | 672 | 492 | 2.497e-04 | 0 | 0 | 0 | 2.148 | 0.003206 | B:114:D | 0.673 |
| dvt_yes | residue | narrow | random | 150 | 1,776 | 672 | 492 | 2.497e-04 | 0 | 0 | 0 | 2.012 | 0.003266 | C:95:L | 0.667 |
| dvt_yes | residue | intermediate | fixed | 155 | 2,051 | 681 | 501 | 2.494e-04 | 0 | 0 | 0 | 1.928 | 0.002638 | B:114:D | 0.675 |
| dvt_yes | residue | intermediate | random | 155 | 2,079 | 681 | 501 | 2.494e-04 | 0 | 0 | 0 | 1.891 | 0.003551 | B:114:D | 0.693 |
| dvt_yes | residue | full | fixed | 159 | 2,632 | 679 | 511 | 2.493e-04 | 0 | 0 | 0 | 1.727 | 0.001279 | B:114:D | 0.661 |
| dvt_yes | residue | full | random | 159 | 2,662 | 679 | 511 | 2.493e-04 | 0 | 0 | 0 | 1.681 | 0.001722 | B:114:D | 0.678 |
| dvt_no | allele | narrow | fixed | 140 | 1,749 | 143 | 143 | 4.385e-04 | 1 | 1 | 1 | 1.725 | 1.908e-04 | DPB1*05:01 | 1.61 |
| dvt_no | allele | narrow | random | 140 | 1,776 | 143 | 143 | 4.385e-04 | 1 | 1 | 1 | 1.662 | 2.815e-04 | DPB1*05:01 | 1.629 |
| dvt_no | allele | intermediate | fixed | 141 | 2,051 | 144 | 144 | 4.346e-04 | 0 | 0 | 0 | 1.665 | 6.024e-04 | DPB1*05:01 | 1.538 |
| dvt_no | allele | intermediate | random | 141 | 2,079 | 144 | 144 | 4.346e-04 | 0 | 0 | 0 | 1.564 | 8.435e-04 | DPB1*05:01 | 1.554 |
| dvt_no | allele | full | fixed | 145 | 2,632 | 142 | 142 | 4.399e-04 | 0 | 0 | 0 | 1.605 | 0.001031 | DPB1*05:01 | 1.496 |
| dvt_no | allele | full | random | 145 | 2,662 | 142 | 142 | 4.399e-04 | 0 | 0 | 0 | 1.549 | 0.001415 | DPB1*05:01 | 1.507 |
| dvt_no | residue | narrow | fixed | 140 | 1,749 | 672 | 491 | 2.499e-04 | 30 | 31 | 0 | 3.279 | 1.152e-04 | DQA1:11:C | 0.596 |
| dvt_no | residue | narrow | random | 140 | 1,776 | 672 | 491 | 2.499e-04 | 33 | 33 | 0 | 3.047 | 1.190e-04 | DQB1:87:Y | 0.543 |
| dvt_no | residue | intermediate | fixed | 141 | 2,051 | 681 | 500 | 2.497e-04 | 0 | 1 | 0 | 2.773 | 2.691e-04 | DQA1:218:Q | 0.617 |
| dvt_no | residue | intermediate | random | 141 | 2,079 | 681 | 500 | 2.497e-04 | 0 | 0 | 0 | 2.584 | 3.290e-04 | DQB1:87:Y | 0.569 |
| dvt_no | residue | full | fixed | 145 | 2,632 | 679 | 510 | 2.495e-04 | 0 | 0 | 0 | 2.497 | 3.104e-04 | DQB1:87:Y | 0.518 |
| dvt_no | residue | full | random | 145 | 2,662 | 679 | 510 | 2.495e-04 | 0 | 0 | 0 | 2.192 | 3.149e-04 | DQB1:87:Y | 0.577 |
| dvt_case_only | allele | narrow | fixed | 150 | 140 | 143 | 143 | 4.490e-04 | 0 | 0 | 0 | 0.965 | 0.00994 | DQB1*05:01 | 3.566 |
| dvt_case_only | allele | narrow | random | 150 | 140 | 143 | 143 | 4.490e-04 | 0 | 0 | 0 | 1.142 | 0.008098 | DQB1*05:01 | 3.159 |
| dvt_case_only | allele | intermediate | fixed | 155 | 141 | 144 | 144 | 4.454e-04 | 0 | 0 | 0 | 0.836 | 0.006312 | DQB1*05:01 | 3.794 |
| dvt_case_only | allele | intermediate | random | 155 | 141 | 144 | 144 | 4.454e-04 | 0 | 0 | 0 | 1.153 | 0.004794 | DQB1*05:01 | 3.305 |
| dvt_case_only | allele | full | fixed | 159 | 145 | 142 | 142 | 4.523e-04 | 0 | 0 | 0 | 1.062 | 0.00845 | DQB1*05:01 | 3.355 |
| dvt_case_only | allele | full | random | 159 | 145 | 142 | 142 | 4.523e-04 | 0 | 0 | 0 | 1.12 | 0.007008 | DQB1*05:01 | 3.019 |
| dvt_case_only | residue | narrow | fixed | 150 | 140 | 672 | 457 | 2.700e-04 | 0 | 0 | 0 | 1.453 | 4.452e-04 | DQB1:-27:A | 0.315 |
| dvt_case_only | residue | narrow | random | 150 | 140 | 672 | 457 | 2.700e-04 | 0 | 0 | 0 | 1.438 | 3.592e-04 | DQB1:-27:A | 0.344 |
| dvt_case_only | residue | intermediate | fixed | 155 | 141 | 681 | 463 | 2.680e-04 | 0 | 0 | 0 | 1.599 | 5.460e-04 | DQB1:-27:A | 0.332 |
| dvt_case_only | residue | intermediate | random | 155 | 141 | 681 | 463 | 2.680e-04 | 0 | 0 | 0 | 1.618 | 4.607e-04 | DQB1:-27:A | 0.36 |
| dvt_case_only | residue | full | fixed | 159 | 145 | 679 | 462 | 2.681e-04 | 0 | 0 | 0 | 1.448 | 0.001246 | DQB1:-27:A | 0.375 |
| dvt_case_only | residue | full | random | 159 | 145 | 679 | 462 | 2.681e-04 | 0 | 0 | 0 | 1.529 | 0.001131 | DQB1:-27:A | 0.398 |

**Markers significant in a DVT set or the case-only test (random model, all cohorts; every threshold and both models in dvt_hits.tsv)**

| class | cohort | model | marker | or_all | p_all | or_dvt_yes | p_dvt_yes | sig_gec_dvt_yes | or_dvt_no | p_dvt_no | sig_gec_dvt_no | or_dvt_case_only | p_dvt_case_only | sig_gec_dvt_case_only |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| allele | narrow | random | DPB1*05:01 | 1.274 | 0.004197 | 1.277 | 0.05819 | 0 | 1.629 | 2.815e-04 | 1 | 0.7735 | 0.142 | 0 |
| allele | intermediate | random | DPB1*05:01 | 1.203 | 0.02306 | 1.162 | 0.2298 | 0 | 1.554 | 8.435e-04 | 0 | 0.7481 | 0.09191 | 0 |
| allele | full | random | DPB1*05:01 | 1.174 | 0.04341 | 1.115 | 0.3736 | 0 | 1.507 | 0.001415 | 0 | 0.7496 | 0.08819 | 0 |
| residue | narrow | random | DQB1:87:Y | 0.7211 | 0.001152 | 0.8081 | 0.1585 | 0 | 0.5427 | 1.190e-04 | 1 | 1.909 | 0.0123 | 0 |
| residue | narrow | random | DQA1:52:R | 1.281 | 0.003538 | 1.274 | 0.06116 | 0 | 1.658 | 1.519e-04 | 1 | 0.745 | 0.1084 | 0 |
| residue | narrow | random | DQA1:11:C | 0.7822 | 0.003645 | 0.7989 | 0.08159 | 0 | 0.6046 | 1.565e-04 | 1 | 1.364 | 0.09116 | 0 |
| residue | narrow | random | DQA1:11:Y | 1.279 | 0.003645 | 1.252 | 0.08159 | 0 | 1.654 | 1.565e-04 | 1 | 0.7331 | 0.09116 | 0 |
| residue | narrow | random | DQA1:18:F | 0.7822 | 0.003645 | 0.7989 | 0.08159 | 0 | 0.6046 | 1.565e-04 | 1 | 1.364 | 0.09116 | 0 |
| residue | narrow | random | DQA1:18:S | 1.279 | 0.003645 | 1.252 | 0.08159 | 0 | 1.654 | 1.565e-04 | 1 | 0.7331 | 0.09116 | 0 |
| residue | narrow | random | DQA1:218:Q | 0.7861 | 0.004416 | 0.8109 | 0.1042 | 0 | 0.6046 | 1.565e-04 | 1 | 1.392 | 0.0719 | 0 |
| residue | narrow | random | DQA1:218:R | 1.272 | 0.004416 | 1.233 | 0.1042 | 0 | 1.654 | 1.565e-04 | 1 | 0.7185 | 0.0719 | 0 |
| residue | narrow | random | DQA1:45:A | 0.7822 | 0.003645 | 0.7989 | 0.08159 | 0 | 0.6046 | 1.565e-04 | 1 | 1.364 | 0.09116 | 0 |
| residue | narrow | random | DQA1:45:V | 1.279 | 0.003645 | 1.252 | 0.08159 | 0 | 1.654 | 1.565e-04 | 1 | 0.7331 | 0.09116 | 0 |
| residue | narrow | random | DQA1:47:R | 0.7822 | 0.003645 | 0.7989 | 0.08159 | 0 | 0.6046 | 1.565e-04 | 1 | 1.364 | 0.09116 | 0 |
| residue | narrow | random | DQA1:48:L | 1.279 | 0.003645 | 1.252 | 0.08159 | 0 | 1.654 | 1.565e-04 | 1 | 0.7331 | 0.09116 | 0 |
| residue | narrow | random | DQA1:48:W | 0.7822 | 0.003645 | 0.7989 | 0.08159 | 0 | 0.6046 | 1.565e-04 | 1 | 1.364 | 0.09116 | 0 |
| residue | narrow | random | DQA1:50:E | 0.7822 | 0.003645 | 0.7989 | 0.08159 | 0 | 0.6046 | 1.565e-04 | 1 | 1.364 | 0.09116 | 0 |
| residue | narrow | random | DQA1:52:S | 0.7822 | 0.003645 | 0.7989 | 0.08159 | 0 | 0.6046 | 1.565e-04 | 1 | 1.364 | 0.09116 | 0 |
| residue | narrow | random | DQA1:53:K | 0.7822 | 0.003645 | 0.7989 | 0.08159 | 0 | 0.6046 | 1.565e-04 | 1 | 1.364 | 0.09116 | 0 |
| residue | narrow | random | DQA1:55:G | 0.7822 | 0.003645 | 0.7989 | 0.08159 | 0 | 0.6046 | 1.565e-04 | 1 | 1.364 | 0.09116 | 0 |
| residue | narrow | random | DQA1:55:R | 1.279 | 0.003645 | 1.252 | 0.08159 | 0 | 1.654 | 1.565e-04 | 1 | 0.7331 | 0.09116 | 0 |
| residue | narrow | random | DQA1:61:F | 1.279 | 0.003645 | 1.252 | 0.08159 | 0 | 1.654 | 1.565e-04 | 1 | 0.7331 | 0.09116 | 0 |
| residue | narrow | random | DQA1:64:R | 0.7822 | 0.003645 | 0.7989 | 0.08159 | 0 | 0.6046 | 1.565e-04 | 1 | 1.364 | 0.09116 | 0 |
| residue | narrow | random | DQA1:64:T | 1.279 | 0.003645 | 1.252 | 0.08159 | 0 | 1.654 | 1.565e-04 | 1 | 0.7331 | 0.09116 | 0 |
| residue | narrow | random | DQA1:66:I | 1.279 | 0.003645 | 1.252 | 0.08159 | 0 | 1.654 | 1.565e-04 | 1 | 0.7331 | 0.09116 | 0 |
| residue | narrow | random | DQA1:66:M | 0.7822 | 0.003645 | 0.7989 | 0.08159 | 0 | 0.6046 | 1.565e-04 | 1 | 1.364 | 0.09116 | 0 |
| residue | narrow | random | DQA1:69:A | 0.7822 | 0.003645 | 0.7989 | 0.08159 | 0 | 0.6046 | 1.565e-04 | 1 | 1.364 | 0.09116 | 0 |
| residue | narrow | random | DQA1:76:M | 0.7822 | 0.003645 | 0.7989 | 0.08159 | 0 | 0.6046 | 1.565e-04 | 1 | 1.364 | 0.09116 | 0 |
| residue | narrow | random | DQA1:80:S | 1.279 | 0.003645 | 1.252 | 0.08159 | 0 | 1.654 | 1.565e-04 | 1 | 0.7331 | 0.09116 | 0 |
| residue | narrow | random | DQA1:80:Y | 0.7822 | 0.003645 | 0.7989 | 0.08159 | 0 | 0.6046 | 1.565e-04 | 1 | 1.364 | 0.09116 | 0 |
| residue | narrow | random | DQA1:61:G | 0.7834 | 0.003892 | 0.8001 | 0.08389 | 0 | 0.6053 | 1.631e-04 | 1 | 1.364 | 0.09116 | 0 |
| residue | narrow | random | DQB1:53:L | 1.282 | 0.003228 | 1.266 | 0.06658 | 0 | 1.645 | 1.754e-04 | 1 | 0.7418 | 0.1061 | 0 |
| residue | narrow | random | DQB1:53:Q | 0.7802 | 0.003228 | 0.7898 | 0.06658 | 0 | 0.6081 | 1.754e-04 | 1 | 1.348 | 0.1061 | 0 |
| residue | narrow | random | DQB1:85:L | 1.281 | 0.003331 | 1.266 | 0.06719 | 0 | 1.644 | 1.797e-04 | 1 | 0.7418 | 0.1061 | 0 |
| residue | narrow | random | DQB1:85:V | 0.7807 | 0.003331 | 0.7901 | 0.06719 | 0 | 0.6083 | 1.797e-04 | 1 | 1.348 | 0.1061 | 0 |
| residue | narrow | random | DQB1:86:E | 1.281 | 0.003331 | 1.266 | 0.06719 | 0 | 1.644 | 1.797e-04 | 1 | 0.7418 | 0.1061 | 0 |
| residue | narrow | random | DQB1:87:L | 1.281 | 0.003331 | 1.266 | 0.06719 | 0 | 1.644 | 1.797e-04 | 1 | 0.7418 | 0.1061 | 0 |
| residue | narrow | random | DQB1:89:G | 0.7807 | 0.003331 | 0.7901 | 0.06719 | 0 | 0.6083 | 1.797e-04 | 1 | 1.348 | 0.1061 | 0 |
| residue | narrow | random | DQB1:89:T | 1.281 | 0.003331 | 1.266 | 0.06719 | 0 | 1.644 | 1.797e-04 | 1 | 0.7418 | 0.1061 | 0 |
| residue | narrow | random | DQB1:90:I | 0.7807 | 0.003331 | 0.7901 | 0.06719 | 0 | 0.6083 | 1.797e-04 | 1 | 1.348 | 0.1061 | 0 |
| residue | narrow | random | DQB1:90:T | 1.281 | 0.003331 | 1.266 | 0.06719 | 0 | 1.644 | 1.797e-04 | 1 | 0.7418 | 0.1061 | 0 |
| residue | narrow | random | DQB1:84:E | 0.7816 | 0.003494 | 0.7909 | 0.06846 | 0 | 0.6089 | 1.860e-04 | 1 | 1.348 | 0.1061 | 0 |
| residue | narrow | random | DQB1:84:Q | 1.279 | 0.003494 | 1.264 | 0.06846 | 0 | 1.642 | 1.860e-04 | 1 | 0.7418 | 0.1061 | 0 |
| residue | narrow | random | DQA1:175:Q | 0.7862 | 0.004358 | 0.8028 | 0.0879 | 0 | 0.6094 | 1.915e-04 | 1 | 1.364 | 0.09116 | 0 |
| residue | narrow | random | DQB1:220:H | 1.273 | 0.004518 | 1.258 | 0.07645 | 0 | 1.639 | 2.151e-04 | 1 | 0.7418 | 0.1061 | 0 |
| residue | narrow | random | DQB1:220:R | 0.7856 | 0.004518 | 0.7948 | 0.07645 | 0 | 0.6101 | 2.151e-04 | 1 | 1.348 | 0.1061 | 0 |
| residue | narrow | random | DQB1:221:H | 1.273 | 0.004518 | 1.258 | 0.07645 | 0 | 1.639 | 2.151e-04 | 1 | 0.7418 | 0.1061 | 0 |
| residue | narrow | random | DQB1:221:Q | 0.7856 | 0.004518 | 0.7948 | 0.07645 | 0 | 0.6101 | 2.151e-04 | 1 | 1.348 | 0.1061 | 0 |
| residue | narrow | random | DQB1:125:A | 1.27 | 0.004793 | 1.255 | 0.07935 | 0 | 1.639 | 2.174e-04 | 1 | 0.7418 | 0.1061 | 0 |
| residue | narrow | random | DQB1:-27:A | 1.314 | 0.02651 | 0.9512 | 0.7845 | 0 | 2.048 | 2.308e-04 | 1 | 0.344 | 3.592e-04 | 0 |
| residue | narrow | random | DQB1:-27:S | 0.7613 | 0.02651 | 1.051 | 0.7845 | 0 | 0.4883 | 2.308e-04 | 1 | 2.907 | 3.592e-04 | 0 |
| residue | narrow | random | DQB1:-5:S | 0.7613 | 0.02651 | 1.051 | 0.7845 | 0 | 0.4883 | 2.308e-04 | 1 | 2.907 | 3.592e-04 | 0 |
| residue | narrow | random | DQB1:-9:I | 0.7613 | 0.02651 | 1.051 | 0.7845 | 0 | 0.4883 | 2.308e-04 | 1 | 2.907 | 3.592e-04 | 0 |
| residue | narrow | random | DQB1:-9:M | 1.314 | 0.02651 | 0.9512 | 0.7845 | 0 | 2.048 | 2.308e-04 | 1 | 0.344 | 3.592e-04 | 0 |
| residue | full | random | DQB1:87:Y | 0.7418 | 0.001852 | 0.8376 | 0.2241 | 0 | 0.5774 | 3.149e-04 | 0 | 1.704 | 0.02864 | 0 |
| residue | intermediate | random | DQB1:87:Y | 0.7366 | 0.001736 | 0.8186 | 0.1756 | 0 | 0.5692 | 3.290e-04 | 0 | 1.765 | 0.02277 | 0 |
| residue | narrow | random | DQB1:182:N | 1.266 | 0.00547 | 1.278 | 0.05833 | 0 | 1.612 | 3.399e-04 | 0 | 0.7697 | 0.1532 | 0 |
| residue | narrow | random | DQB1:182:S | 0.79 | 0.00547 | 0.7823 | 0.05833 | 0 | 0.6202 | 3.399e-04 | 0 | 1.299 | 0.1532 | 0 |
| residue | intermediate | random | DQA1:218:Q | 0.8144 | 0.01203 | 0.8378 | 0.1592 | 0 | 0.6245 | 3.412e-04 | 0 | 1.392 | 0.0695 | 0 |
| residue | intermediate | random | DQA1:218:R | 1.228 | 0.01203 | 1.194 | 0.1592 | 0 | 1.601 | 3.412e-04 | 0 | 0.7186 | 0.0695 | 0 |
| residue | intermediate | random | DQA1:52:R | 1.24 | 0.008725 | 1.245 | 0.08131 | 0 | 1.602 | 3.447e-04 | 0 | 0.7528 | 0.1173 | 0 |
| residue | intermediate | random | DQA1:11:C | 0.807 | 0.008659 | 0.8154 | 0.104 | 0 | 0.6251 | 3.483e-04 | 0 | 1.349 | 0.09965 | 0 |
| residue | intermediate | random | DQA1:11:Y | 1.239 | 0.008659 | 1.226 | 0.104 | 0 | 1.6 | 3.483e-04 | 0 | 0.7413 | 0.09965 | 0 |
| residue | intermediate | random | DQA1:18:F | 0.807 | 0.008659 | 0.8154 | 0.104 | 0 | 0.6251 | 3.483e-04 | 0 | 1.349 | 0.09965 | 0 |
| residue | intermediate | random | DQA1:18:S | 1.239 | 0.008659 | 1.226 | 0.104 | 0 | 1.6 | 3.483e-04 | 0 | 0.7413 | 0.09965 | 0 |
| residue | intermediate | random | DQA1:45:A | 0.807 | 0.008659 | 0.8154 | 0.104 | 0 | 0.6251 | 3.483e-04 | 0 | 1.349 | 0.09965 | 0 |
| residue | intermediate | random | DQA1:45:V | 1.239 | 0.008659 | 1.226 | 0.104 | 0 | 1.6 | 3.483e-04 | 0 | 0.7413 | 0.09965 | 0 |
| residue | intermediate | random | DQA1:47:R | 0.807 | 0.008659 | 0.8154 | 0.104 | 0 | 0.6251 | 3.483e-04 | 0 | 1.349 | 0.09965 | 0 |
| residue | intermediate | random | DQA1:48:L | 1.239 | 0.008659 | 1.226 | 0.104 | 0 | 1.6 | 3.483e-04 | 0 | 0.7413 | 0.09965 | 0 |
| residue | intermediate | random | DQA1:48:W | 0.807 | 0.008659 | 0.8154 | 0.104 | 0 | 0.6251 | 3.483e-04 | 0 | 1.349 | 0.09965 | 0 |
| residue | intermediate | random | DQA1:50:E | 0.807 | 0.008659 | 0.8154 | 0.104 | 0 | 0.6251 | 3.483e-04 | 0 | 1.349 | 0.09965 | 0 |
| residue | intermediate | random | DQA1:52:S | 0.807 | 0.008659 | 0.8154 | 0.104 | 0 | 0.6251 | 3.483e-04 | 0 | 1.349 | 0.09965 | 0 |
| residue | intermediate | random | DQA1:53:K | 0.807 | 0.008659 | 0.8154 | 0.104 | 0 | 0.6251 | 3.483e-04 | 0 | 1.349 | 0.09965 | 0 |
| residue | intermediate | random | DQA1:55:G | 0.807 | 0.008659 | 0.8154 | 0.104 | 0 | 0.6251 | 3.483e-04 | 0 | 1.349 | 0.09965 | 0 |
| residue | intermediate | random | DQA1:55:R | 1.239 | 0.008659 | 1.226 | 0.104 | 0 | 1.6 | 3.483e-04 | 0 | 0.7413 | 0.09965 | 0 |
| residue | intermediate | random | DQA1:61:F | 1.239 | 0.008659 | 1.226 | 0.104 | 0 | 1.6 | 3.483e-04 | 0 | 0.7413 | 0.09965 | 0 |
| residue | intermediate | random | DQA1:64:R | 0.807 | 0.008659 | 0.8154 | 0.104 | 0 | 0.6251 | 3.483e-04 | 0 | 1.349 | 0.09965 | 0 |
| residue | intermediate | random | DQA1:64:T | 1.239 | 0.008659 | 1.226 | 0.104 | 0 | 1.6 | 3.483e-04 | 0 | 0.7413 | 0.09965 | 0 |
| residue | intermediate | random | DQA1:66:I | 1.239 | 0.008659 | 1.226 | 0.104 | 0 | 1.6 | 3.483e-04 | 0 | 0.7413 | 0.09965 | 0 |
| residue | intermediate | random | DQA1:66:M | 0.807 | 0.008659 | 0.8154 | 0.104 | 0 | 0.6251 | 3.483e-04 | 0 | 1.349 | 0.09965 | 0 |
| residue | intermediate | random | DQA1:69:A | 0.807 | 0.008659 | 0.8154 | 0.104 | 0 | 0.6251 | 3.483e-04 | 0 | 1.349 | 0.09965 | 0 |
| residue | intermediate | random | DQA1:76:M | 0.807 | 0.008659 | 0.8154 | 0.104 | 0 | 0.6251 | 3.483e-04 | 0 | 1.349 | 0.09965 | 0 |
| residue | intermediate | random | DQA1:80:S | 1.239 | 0.008659 | 1.226 | 0.104 | 0 | 1.6 | 3.483e-04 | 0 | 0.7413 | 0.09965 | 0 |
| residue | intermediate | random | DQA1:80:Y | 0.807 | 0.008659 | 0.8154 | 0.104 | 0 | 0.6251 | 3.483e-04 | 0 | 1.349 | 0.09965 | 0 |
| residue | intermediate | random | DQA1:61:G | 0.8079 | 0.00907 | 0.8163 | 0.1062 | 0 | 0.6256 | 3.587e-04 | 0 | 1.349 | 0.09965 | 0 |
| residue | intermediate | random | DQB1:53:L | 1.243 | 0.007619 | 1.241 | 0.08487 | 0 | 1.593 | 3.771e-04 | 0 | 0.7505 | 0.1165 | 0 |
| residue | intermediate | random | DQB1:53:Q | 0.8045 | 0.007619 | 0.8058 | 0.08487 | 0 | 0.6277 | 3.771e-04 | 0 | 1.332 | 0.1165 | 0 |
| residue | intermediate | random | DQB1:85:L | 1.242 | 0.007778 | 1.241 | 0.08544 | 0 | 1.593 | 3.831e-04 | 0 | 0.7505 | 0.1165 | 0 |
| residue | intermediate | random | DQB1:85:V | 0.8049 | 0.007778 | 0.8061 | 0.08544 | 0 | 0.6279 | 3.831e-04 | 0 | 1.332 | 0.1165 | 0 |
| residue | intermediate | random | DQB1:86:E | 1.242 | 0.007778 | 1.241 | 0.08544 | 0 | 1.593 | 3.831e-04 | 0 | 0.7505 | 0.1165 | 0 |
| residue | intermediate | random | DQB1:87:L | 1.242 | 0.007778 | 1.241 | 0.08544 | 0 | 1.593 | 3.831e-04 | 0 | 0.7505 | 0.1165 | 0 |
| residue | intermediate | random | DQB1:89:G | 0.8049 | 0.007778 | 0.8061 | 0.08544 | 0 | 0.6279 | 3.831e-04 | 0 | 1.332 | 0.1165 | 0 |
| residue | intermediate | random | DQB1:89:T | 1.242 | 0.007778 | 1.241 | 0.08544 | 0 | 1.593 | 3.831e-04 | 0 | 0.7505 | 0.1165 | 0 |
| residue | intermediate | random | DQB1:90:T | 1.242 | 0.007778 | 1.241 | 0.08544 | 0 | 1.593 | 3.831e-04 | 0 | 0.7505 | 0.1165 | 0 |
| residue | intermediate | random | DQB1:84:E | 0.8055 | 0.007998 | 0.8066 | 0.08642 | 0 | 0.6283 | 3.910e-04 | 0 | 1.332 | 0.1165 | 0 |
| residue | intermediate | random | DQB1:84:Q | 1.242 | 0.007998 | 1.24 | 0.08642 | 0 | 1.592 | 3.910e-04 | 0 | 0.7505 | 0.1165 | 0 |
| residue | intermediate | random | DQB1:90:I | 0.8056 | 0.00809 | 0.8067 | 0.08672 | 0 | 0.6283 | 3.932e-04 | 0 | 1.332 | 0.1165 | 0 |
| residue | intermediate | random | DQB1:-27:A | 1.295 | 0.02892 | 0.9368 | 0.7133 | 0 | 1.966 | 4.101e-04 | 0 | 0.3603 | 4.607e-04 | 0 |
| residue | intermediate | random | DQB1:-27:S | 0.7719 | 0.02892 | 1.067 | 0.7133 | 0 | 0.5085 | 4.101e-04 | 0 | 2.775 | 4.607e-04 | 0 |
| residue | intermediate | random | DQB1:-5:S | 0.7719 | 0.02892 | 1.067 | 0.7133 | 0 | 0.5085 | 4.101e-04 | 0 | 2.775 | 4.607e-04 | 0 |
| residue | intermediate | random | DQB1:-9:I | 0.7719 | 0.02892 | 1.067 | 0.7133 | 0 | 0.5085 | 4.101e-04 | 0 | 2.775 | 4.607e-04 | 0 |
| residue | intermediate | random | DQB1:-9:M | 1.295 | 0.02892 | 0.9368 | 0.7133 | 0 | 1.966 | 4.101e-04 | 0 | 0.3603 | 4.607e-04 | 0 |
| residue | full | random | DQA1:218:Q | 0.8223 | 0.01468 | 0.8471 | 0.1782 | 0 | 0.6365 | 4.168e-04 | 0 | 1.372 | 0.07776 | 0 |
| residue | full | random | DQA1:218:R | 1.216 | 0.01468 | 1.181 | 0.1782 | 0 | 1.571 | 4.168e-04 | 0 | 0.7289 | 0.07776 | 0 |
| residue | full | random | DQA1:11:C | 0.8147 | 0.01049 | 0.8247 | 0.1177 | 0 | 0.637 | 4.253e-04 | 0 | 1.33 | 0.1108 | 0 |
| residue | full | random | DQA1:11:Y | 1.227 | 0.01049 | 1.212 | 0.1177 | 0 | 1.57 | 4.253e-04 | 0 | 0.7517 | 0.1108 | 0 |
| residue | full | random | DQA1:18:F | 0.8147 | 0.01049 | 0.8247 | 0.1177 | 0 | 0.637 | 4.253e-04 | 0 | 1.33 | 0.1108 | 0 |
| residue | full | random | DQA1:18:S | 1.227 | 0.01049 | 1.212 | 0.1177 | 0 | 1.57 | 4.253e-04 | 0 | 0.7517 | 0.1108 | 0 |
| residue | full | random | DQA1:45:A | 0.8147 | 0.01049 | 0.8247 | 0.1177 | 0 | 0.637 | 4.253e-04 | 0 | 1.33 | 0.1108 | 0 |
| residue | full | random | DQA1:45:V | 1.227 | 0.01049 | 1.212 | 0.1177 | 0 | 1.57 | 4.253e-04 | 0 | 0.7517 | 0.1108 | 0 |
| residue | full | random | DQA1:47:R | 0.8147 | 0.01049 | 0.8247 | 0.1177 | 0 | 0.637 | 4.253e-04 | 0 | 1.33 | 0.1108 | 0 |
| residue | full | random | DQA1:48:L | 1.227 | 0.01049 | 1.212 | 0.1177 | 0 | 1.57 | 4.253e-04 | 0 | 0.7517 | 0.1108 | 0 |
| residue | full | random | DQA1:48:W | 0.8147 | 0.01049 | 0.8247 | 0.1177 | 0 | 0.637 | 4.253e-04 | 0 | 1.33 | 0.1108 | 0 |
| residue | full | random | DQA1:50:E | 0.8147 | 0.01049 | 0.8247 | 0.1177 | 0 | 0.637 | 4.253e-04 | 0 | 1.33 | 0.1108 | 0 |
| residue | full | random | DQA1:52:S | 0.8147 | 0.01049 | 0.8247 | 0.1177 | 0 | 0.637 | 4.253e-04 | 0 | 1.33 | 0.1108 | 0 |
| residue | full | random | DQA1:53:K | 0.8147 | 0.01049 | 0.8247 | 0.1177 | 0 | 0.637 | 4.253e-04 | 0 | 1.33 | 0.1108 | 0 |
| residue | full | random | DQA1:55:G | 0.8147 | 0.01049 | 0.8247 | 0.1177 | 0 | 0.637 | 4.253e-04 | 0 | 1.33 | 0.1108 | 0 |
| residue | full | random | DQA1:55:R | 1.227 | 0.01049 | 1.212 | 0.1177 | 0 | 1.57 | 4.253e-04 | 0 | 0.7517 | 0.1108 | 0 |
| residue | full | random | DQA1:61:F | 1.227 | 0.01049 | 1.212 | 0.1177 | 0 | 1.57 | 4.253e-04 | 0 | 0.7517 | 0.1108 | 0 |
| residue | full | random | DQA1:64:R | 0.8147 | 0.01049 | 0.8247 | 0.1177 | 0 | 0.637 | 4.253e-04 | 0 | 1.33 | 0.1108 | 0 |
| residue | full | random | DQA1:64:T | 1.227 | 0.01049 | 1.212 | 0.1177 | 0 | 1.57 | 4.253e-04 | 0 | 0.7517 | 0.1108 | 0 |
| residue | full | random | DQA1:66:I | 1.227 | 0.01049 | 1.212 | 0.1177 | 0 | 1.57 | 4.253e-04 | 0 | 0.7517 | 0.1108 | 0 |
| residue | full | random | DQA1:66:M | 0.8147 | 0.01049 | 0.8247 | 0.1177 | 0 | 0.637 | 4.253e-04 | 0 | 1.33 | 0.1108 | 0 |
| residue | full | random | DQA1:69:A | 0.8147 | 0.01049 | 0.8247 | 0.1177 | 0 | 0.637 | 4.253e-04 | 0 | 1.33 | 0.1108 | 0 |
| residue | full | random | DQA1:76:M | 0.8147 | 0.01049 | 0.8247 | 0.1177 | 0 | 0.637 | 4.253e-04 | 0 | 1.33 | 0.1108 | 0 |
| residue | full | random | DQA1:80:S | 1.227 | 0.01049 | 1.212 | 0.1177 | 0 | 1.57 | 4.253e-04 | 0 | 0.7517 | 0.1108 | 0 |
| residue | full | random | DQA1:80:Y | 0.8147 | 0.01049 | 0.8247 | 0.1177 | 0 | 0.637 | 4.253e-04 | 0 | 1.33 | 0.1108 | 0 |
| residue | full | random | DQA1:61:G | 0.8154 | 0.01087 | 0.8254 | 0.1194 | 0 | 0.6374 | 4.360e-04 | 0 | 1.33 | 0.1108 | 0 |
| residue | intermediate | random | DQB1:220:H | 1.236 | 0.009855 | 1.234 | 0.09561 | 0 | 1.59 | 4.381e-04 | 0 | 0.7505 | 0.1165 | 0 |
| residue | intermediate | random | DQB1:220:R | 0.8089 | 0.009855 | 0.8103 | 0.09561 | 0 | 0.6289 | 4.381e-04 | 0 | 1.332 | 0.1165 | 0 |
| residue | intermediate | random | DQB1:221:H | 1.236 | 0.009855 | 1.234 | 0.09561 | 0 | 1.59 | 4.381e-04 | 0 | 0.7505 | 0.1165 | 0 |
| residue | intermediate | random | DQB1:221:Q | 0.8089 | 0.009855 | 0.8103 | 0.09561 | 0 | 0.6289 | 4.381e-04 | 0 | 1.332 | 0.1165 | 0 |
| residue | intermediate | random | DQA1:175:Q | 0.8132 | 0.01117 | 0.8214 | 0.1161 | 0 | 0.6314 | 4.493e-04 | 0 | 1.349 | 0.09965 | 0 |
| residue | full | random | DQA1:52:R | 1.225 | 0.01144 | 1.228 | 0.0965 | 0 | 1.567 | 4.514e-04 | 0 | 0.7639 | 0.1312 | 0 |
| residue | intermediate | random | DQB1:125:A | 1.232 | 0.01087 | 1.231 | 0.09972 | 0 | 1.587 | 4.586e-04 | 0 | 0.7505 | 0.1165 | 0 |
| residue | full | random | DQB1:53:L | 1.229 | 0.009795 | 1.224 | 0.0999 | 0 | 1.562 | 4.760e-04 | 0 | 0.7605 | 0.1281 | 0 |
| residue | full | random | DQB1:53:Q | 0.8135 | 0.009795 | 0.817 | 0.0999 | 0 | 0.6404 | 4.760e-04 | 0 | 1.315 | 0.1281 | 0 |
| residue | full | random | DQB1:85:L | 1.229 | 0.009927 | 1.224 | 0.1004 | 0 | 1.561 | 4.807e-04 | 0 | 0.7605 | 0.1281 | 0 |
| residue | full | random | DQB1:85:V | 0.8138 | 0.009927 | 0.8172 | 0.1004 | 0 | 0.6405 | 4.807e-04 | 0 | 1.315 | 0.1281 | 0 |
| residue | full | random | DQB1:86:E | 1.229 | 0.009927 | 1.224 | 0.1004 | 0 | 1.561 | 4.807e-04 | 0 | 0.7605 | 0.1281 | 0 |
| residue | full | random | DQB1:87:L | 1.229 | 0.009927 | 1.224 | 0.1004 | 0 | 1.561 | 4.807e-04 | 0 | 0.7605 | 0.1281 | 0 |
| residue | full | random | DQB1:89:G | 0.8138 | 0.009927 | 0.8172 | 0.1004 | 0 | 0.6405 | 4.807e-04 | 0 | 1.315 | 0.1281 | 0 |
| residue | full | random | DQB1:89:T | 1.229 | 0.009927 | 1.224 | 0.1004 | 0 | 1.561 | 4.807e-04 | 0 | 0.7605 | 0.1281 | 0 |
| residue | full | random | DQB1:90:T | 1.229 | 0.009927 | 1.224 | 0.1004 | 0 | 1.561 | 4.807e-04 | 0 | 0.7605 | 0.1281 | 0 |
| residue | full | random | DQB1:90:I | 0.8142 | 0.01014 | 0.8175 | 0.1012 | 0 | 0.6407 | 4.876e-04 | 0 | 1.315 | 0.1281 | 0 |
| residue | full | random | DQB1:84:E | 0.8143 | 0.01016 | 0.8176 | 0.1013 | 0 | 0.6408 | 4.894e-04 | 0 | 1.315 | 0.1281 | 0 |
| residue | full | random | DQB1:84:Q | 1.228 | 0.01016 | 1.223 | 0.1013 | 0 | 1.56 | 4.894e-04 | 0 | 0.7605 | 0.1281 | 0 |
| residue | full | random | DQA1:175:Q | 0.8189 | 0.01247 | 0.8289 | 0.1268 | 0 | 0.6412 | 5.034e-04 | 0 | 1.33 | 0.1108 | 0 |
| residue | full | random | DQB1:220:H | 1.224 | 0.01192 | 1.219 | 0.1091 | 0 | 1.559 | 5.433e-04 | 0 | 0.7605 | 0.1281 | 0 |
| residue | full | random | DQB1:220:R | 0.8167 | 0.01192 | 0.8201 | 0.1091 | 0 | 0.6413 | 5.433e-04 | 0 | 1.315 | 0.1281 | 0 |
| residue | full | random | DQB1:221:H | 1.224 | 0.01192 | 1.219 | 0.1091 | 0 | 1.559 | 5.433e-04 | 0 | 0.7605 | 0.1281 | 0 |
| residue | full | random | DQB1:221:Q | 0.8167 | 0.01192 | 0.8201 | 0.1091 | 0 | 0.6413 | 5.433e-04 | 0 | 1.315 | 0.1281 | 0 |
| residue | full | random | DQB1:125:A | 1.221 | 0.01316 | 1.216 | 0.1132 | 0 | 1.556 | 5.728e-04 | 0 | 0.7605 | 0.1281 | 0 |
| residue | full | random | DQB1:-27:A | 1.289 | 0.02843 | 0.9322 | 0.6873 | 0 | 1.885 | 6.192e-04 | 0 | 0.3983 | 0.001131 | 0 |
| residue | full | random | DQB1:-27:S | 0.7759 | 0.02843 | 1.073 | 0.6873 | 0 | 0.5304 | 6.192e-04 | 0 | 2.511 | 0.001131 | 0 |
| residue | full | random | DQB1:-5:S | 0.7759 | 0.02843 | 1.073 | 0.6873 | 0 | 0.5304 | 6.192e-04 | 0 | 2.511 | 0.001131 | 0 |
| residue | full | random | DQB1:-9:I | 0.7759 | 0.02843 | 1.073 | 0.6873 | 0 | 0.5304 | 6.192e-04 | 0 | 2.511 | 0.001131 | 0 |
| residue | full | random | DQB1:-9:M | 1.289 | 0.02843 | 0.9322 | 0.6873 | 0 | 1.885 | 6.192e-04 | 0 | 0.3983 | 0.001131 | 0 |
| residue | intermediate | random | DQB1:182:N | 1.23 | 0.01145 | 1.255 | 0.07164 | 0 | 1.565 | 6.594e-04 | 0 | 0.778 | 0.1663 | 0 |
| residue | intermediate | random | DQB1:182:S | 0.8128 | 0.01145 | 0.7967 | 0.07164 | 0 | 0.6388 | 6.594e-04 | 0 | 1.285 | 0.1663 | 0 |
| residue | full | random | DQB1:182:N | 1.217 | 0.01431 | 1.239 | 0.08352 | 0 | 1.534 | 8.316e-04 | 0 | 0.7874 | 0.1804 | 0 |
| residue | full | random | DQB1:182:S | 0.8214 | 0.01431 | 0.8074 | 0.08352 | 0 | 0.6517 | 8.316e-04 | 0 | 1.27 | 0.1804 | 0 |

## How to read it

1. Points on the identity line: the same effect with and without DVT.
2. A coloured point far from the line: significant in one set with a different effect in the other.
3. A black ring: the case-only test supports a difference for that marker.
4. Check dvt_hits.tsv for the same marker in the other cohorts, thresholds and the fixed model.

## What this figure does *not* establish

- A difference read off the scatter alone: the two axes share controls; use the case-only test.
- DVT status is recorded only for cases; unknown-status cases (about 30 %) are in no DVT set.
- All cases were sequenced on DNBSeq or NovaSeq and all controls on HiSeq X; that split is in the three case-control sets. The case-only test avoids it, but the platform mix differs between cases with and without DVT, so the case-only test is not free of platform effects.
- The three cohorts share almost all their cases, so the case-only test is nearly the same analysis three times.

---

Methods and rationale: [`METHODS.md`](../../docs/METHODS.md)

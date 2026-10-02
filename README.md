# Code for: Integrated Biofluid Proteomics Identified Dynamic Functional Biomarkers of LRRK2-Linked Parkinson's Disease Progression


## Data access
Raw data are not included in this repository; they are available from the original sources under their data use agreements.
- PPMI (www.ppmi-info.org): CSF SomaScan (Project 151), CSF mass spectrometry (Project 177),
  urine mass spectrometry (Project 190), CSF Olink Explore (Project 9000). Data downloaded on [YYYY-MM-DD].
- LRRK2 Cohort Consortium (LCC) CSF and urine mass spectrometry: [request route].
- Columbia urine mass spectrometry: [request route].
- Mouse urine proteomics (Lrrk2 G2019S vs WT): [included in /data_mouse, or available on request].

Place input files in `data/`. Participant exclusion lists are read from `data/excluded_ids.csv`,
which is not distributed; exclusion rules are described in the Methods.

## Scripts
Run scripts in numeric order. * = added in revision.

### CSF, PPMI SomaScan (Project 151)
| Script | Analysis | Output |
|---|---|---|
| 01_CSF_DAP_PPMI.R | Differential abundance, LRRK2 PD vs HC (limma, age + sex); GO enrichment | Data file S[ ] |
| 02_CSF_WGCNA_PPMI.R | Co-expression network, module–trait correlation, module membership | Fig. [ ], Data file S[ ] |
| 03_ML_PPMI_to_LCC.R | Elastic-net model trained on PPMI, validated on LCC; label-permutation test (10,000) | Fig. 1E, Table 2 |
| 04_CSF_DAP_LCC.R* | LCC differential abundance (age, sex, study centre); overlap with PPMI DAPs | Data file S[ ] |
| 05_crossplatform_SomaScan_Olink.R* | sPD vs HC model, SomaScan to Olink transfer (cross-fitted, n = 161) | Fig. S[ ] |

### CSF, PPMI mass spectrometry (Project 177): longitudinal
| Script | Analysis | Output |
|---|---|---|
| 06_LMM_sex_interaction.R | Linear mixed models, group × time × sex (lme4/lmerTest) | Data file S11 |
| 07_LMM_male.R* | Males only; model selection by AIC (ML), REML refit | Data file S12 |
| 08_LMM_female.R* | Females only; same pipeline | Data file S14 |

### Urine
| Script | Analysis | Output |
|---|---|---|
| 09_urine_DAP_PPMI.R | Differential abundance, LRRK2 PD vs HC (Project 190) | Data file S15 |
| 10_urine_overlap_cohorts.R* | DAP overlap across PPMI, Columbia, LCC (counts, Jaccard index) | Data files S19, S[ ] |
| 11_CSF_vs_urine_effects.R* | Comparison of CSF and urine effect sizes | Fig. [ ], Data file S[ ] |
| 12_mouse_urine_DEP.R | Lrrk2 G2019S vs WT mice (limma); GO enrichment | Fig. S8, Data file S29 |

## Software
R [4.x.x]. Main packages: limma, WGCNA, glmnet, pROC, lme4, lmerTest, clusterProfiler,
ggplot2, EnhancedVolcano. Exact versions are in `sessionInfo.txt`.
Random seeds are set in each script where results depend on them.
Some panels and group comparisons were made in GraphPad Prism 10.

## License
MIT

## Contact
Cong Xiao, Icahn School of Medicine at Mount Sinai — cong.xiao@mssm.edu

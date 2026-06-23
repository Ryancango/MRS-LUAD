# MRS-LUAD code release

This repository contains the curated analysis code for the manuscript **"A membrane-related multi-omics signature links metabolic reprogramming and immune evasion in lung adenocarcinoma"**.

The release is intentionally compact. It preserves the scripts needed to audit and reproduce the main computational logic while excluding raw public datasets, large intermediate R objects, generated figures, manuscript files and exploratory one-off scripts.

## What is included

- Discovery-cohort preprocessing, ComBat batch correction, WGCNA module detection and candidate-gene intersection scripts.
- The modified machine-learning survival framework and the script used to compare 117 model combinations.
- MRS risk-score calculation and validation plotting scripts.
- GSE131907 single-cell scripts for MRS scoring, scFEA-derived metabolic flux summaries, epithelial metabolic programs, Monocle2 pseudotime, CellChat plots and T-cell exhaustion analyses.
- TCGA Pan-Cancer/Xena scripts for OS, DSS, DFI and PFI survival-association summaries.
- Metadata tables for the final 10-gene panel, risk-score coefficients, public datasets and source-script mapping.

## Final MRS score

```text
MRS = 0.342398338*RRAS - 0.249540066*GPM6B + 0.034609870*FMO3
      - 0.205775214*CHRDL1 + 0.033408756*ADRB2 - 0.092544047*LPL
      + 0.007075607*ANGPT1 + 0.022607685*SPOCK2
      - 0.063043353*P2RY14 + 0.091332717*ADH1B
```

The coefficient table is available at `metadata/gene_coefficients.tsv`, and the helper function is available at `R/utils_mrs_score.R`.

## Reproducibility notes

1. Large public datasets and intermediate objects are not stored in this repository. Download or regenerate them from the public resources listed in `metadata/datasets_used.tsv`.
2. Run scripts from the repository root where possible.
3. Set local paths before rerunning analyses:

```r
Sys.setenv(MRS_PROJECT_ROOT = "/path/to/recreated/MRS_project")
Sys.setenv(SCFEA_ROOT = "/path/to/scFEA_project")
source("R/config.R")
```

4. The original project used a deterministic seed of `5201314` for training/validation splits in the machine-learning workflow and `1314` inside the model-development function where specified.
5. The single-cell analysis retains original GSE131907 cell annotations; no copy-number inference is used for the final manuscript claims.
6. Pan-cancer scripts assess the prognostic association of the 10-gene panel across TCGA tumour types and endpoints. Endpoint-specific Cox fitting is implemented in the Xena analysis script.

## Main folders

```text
R/                         shared configuration and MRS scoring helper
metadata/                  final gene panel, coefficients, dataset list and source-script map
scripts/01_discovery/      discovery cohorts, WGCNA and GO annotation
scripts/02_machine_learning/ model training, ranking and validation figures
scripts/03_single_cell/    GSE131907 single-cell analyses
scripts/04_pancancer/      TCGA Pan-Cancer/Xena analyses
environment/               R package list
```

# Kinhibit

**Selectivity-aware kinase inhibitor annotation from phosphoproteomics activity inference, with quantitative cardiotoxicity safety scoring.**

An R module that takes kinase activity scores from phosphoproteomics inference tools (phosphonetworks, decoupleR, KSEA) and returns a ranked, safety-scored list of inhibitor compounds for experimental validation.

## What it does

Given a list of active kinases with activity scores (NES) and p-values, kinhibit:

1. **Annotates with inhibitors** — queries ChEMBL bioactivity data (Kd/Ki/IC50 from binding and functional assays) for each kinase
2. **Scores selectivity** — integrates Klaeger et al. 2017 kinome-wide selectivity profiles (S(10) scoring) and Reinecke et al. 2023 (ProteomicsDB) expanded kinome (1,183 compounds)
3. **Adds clinical context** — queries Open Targets for clinical phase, mechanism of action, indications, and cardiovascular relevance
4. **Predicts rescue potential** — Hill equation model estimating fractional pathway rescue at standard concentrations (100 nM, 1 μM, 10 μM)
5. **Scores cardiotoxicity** — FAERS-trained elastic net model (50 drugs × 250 kinase features) predicting cardiac adverse event risk from kinome-wide inhibition profiles. Validated with zero misclassifications between known cardiotoxic and safe drugs
6. **Produces three ranked shortlists** — pharmacology (pure potency), translational (clinical phase bonus), safety-adjusted (cardiotoxicity penalty)
7. **Generates publication figures** — top-20 bar charts, cardiotox kinase coefficients, model validation scatter, score distributions

## Quick start

```r
source("Kinase_Inhibitor_Annotation.R")

# Load phosphonetworks output
kin <- load_phosphonetworks("kinase_activity_DSP_D15.csv",
                            resource = "combined",
                            nes_cut  = 2,
                            pval_cut = 0.05)

# Annotate with inhibitors (first run: ~5-10 min for ChEMBL queries; cached thereafter)
anno <- annotate_kinase_inhibitors(kin, sources = c("chembl", "klaeger"), top_n = 50)

# See run_example.R for full pipeline including Open Targets,
# cardiotoxicity scoring, rescue prediction, and figure generation
```

For the complete end-to-end pipeline, see `run_example.R`.

## Input format

CSV with columns:

| Column | Description |
|--------|-------------|
| `gene` or `kinase` | Gene symbol |
| `UniProt` or `kinase` | UniProt accession |
| `activity` | Kinase activity z-score (NES) |
| `pvalue` | Statistical significance (optional) |
| `resource` | Inference resource: "combined", "kinlib", "literature", "phosformer" |

The module handles column name variations automatically via `load_phosphonetworks()`.

## Data sources

| Source | Type | N compounds | Reference |
|--------|------|-------------|-----------|
| ChEMBL | Bioactivity (Kd/Ki/IC50) | ~10,000+ per query | EBI ChEMBL REST API |
| Klaeger 2017 | Kinome-wide Kd + selectivity | 222 clinical KIs | Science 2017;358:eaan4368 |
| ProteomicsDB / Reinecke 2023 | Kinome-wide Kd | 1,183 compounds | Nat Chem Biol 2024;20:577 |
| Open Targets | Clinical phase, indications | 5,400+ | platform.opentargets.org |
| FAERS | Cardiac adverse event fractions | 50 drugs (training) | openFDA API |
| van Hasselt 2020 | FAERS risk scores (validation) | 26 drugs | Nat Commun 2020;11:4868 |
| Jin 2020 | CAE endpoint labels | 38+17 drugs | Front Pharmacol 2020;11:897 |

## Cardiotoxicity model

Elastic net (alpha=0.5) trained on 50 FAERS-labeled drugs with kinome-wide inhibition profiles from Klaeger 2017 + ProteomicsDB 2023 (combined: 1,306 compounds × 318 kinases). The model identifies 11 kinases whose inhibition predicts cardiac adverse events:

**Cardiotoxic (positive weight):** RET, ABL1, DDR1, RIPK2, PTPN1, TNIK, BCR, UNC119, MAPK14, FRK

**Cardioprotective (negative weight):** CSNK2A2, CDK4

Validation: all 10 known cardiotoxic drugs ranked above all 7 known safe drugs (gap = 0.183, zero misclassifications).

**Important:** This is a pharmacological prioritization heuristic, not a clinical safety prediction. Cardiotoxicity scoring works best for compounds with kinome-wide profiles (Klaeger/ProteomicsDB). Compounds with ChEMBL-only data (single-target) receive scores based on primary target only — off-target risk is unassessed.

## Outputs

### CSV tables
- `FINAL_complete_table.csv` — all compounds with potency, selectivity, clinical phase, cardiotox score, three composite rankings
- `FINAL_top20_pharmacology.csv` — top 20 by pure potency + kinase significance
- `FINAL_top20_translational.csv` — top 20 including clinical phase bonus
- `FINAL_top20_safety.csv` — top 20 with cardiotoxicity penalty applied

### Publication figures (PNG + PDF)
- Top 20 bar charts for each ranking with cardiotox/phase/CV annotation
- Cardiotoxicity-predictive kinase coefficients
- Model validation scatter (observed vs predicted FAERS scores)
- Cardiotoxicity score distribution

## Dependencies

```r
# CRAN
install.packages(c("dplyr", "tidyr", "readr", "stringr", "purrr",
                   "jsonlite", "ggplot2", "readxl", "glmnet",
                   "ggrepel", "pheatmap", "digest"))

# Bioconductor
BiocManager::install(c("org.Hs.eg.db", "AnnotationDbi"))
```

## Required data files

Place in `data/` directory:

1. `klaeger_2017_kinome_matrix.tsv` — reshaped from Klaeger et al. 2017 Science Table S2
2. ProteomicsDB supplementary tables (Reinecke et al. 2023 Nat Chem Biol):
   - `41589_2023_1459_MOESM3_ESM.xlsx` (drug matrix)
3. van Hasselt et al. 2020 source data:
   - `41467_2020_18396_MOESM4_ESM.xlsx` (FAERS risk scores)

## Caching

All ChEMBL and Open Targets API responses are cached to `cache/kinase_inhibitors/` as RDS files with a 30-day TTL. First run takes ~10 minutes; subsequent runs are near-instant. FAERS queries are also cached.

## Citation

If you use kinhibit, please cite:

> [Marcos Sande-Melon]. kinhibit: selectivity-aware kinase inhibitor annotation with cardiotoxicity scoring for phosphoproteomics. Zenodo. DOI: [to be assigned]

## License

MIT

## Authors

Developed at the Murdoch Children's Research Institute (MCRI), Melbourne, Australia.

## Acknowledgements

This tool integrates data from ChEMBL (EBI), ProteomicsDB (Kuster lab, TU Munich), Open Targets (EBI/Sanger), openFDA (FDA), Klaeger et al. 2017, Reinecke et al. 2023, van Hasselt et al. 2020, and Jin et al. 2020.

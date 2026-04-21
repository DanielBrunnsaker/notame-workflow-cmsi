# notame-workflow-cmsi

Post-MSDIAL preprocessing pipeline for untargeted LC-MS metabolomics data. Converts MSDIAL alignment exports to a standardised format, applies feature quality filters, and evaluates multiple drift and batch correction strategies in parallel.

Built around the [notame](https://github.com/antonvsdata/notame) R package.

## Overview

1. **MSDIAL conversion** — parses MSDIAL alignment exports, extracts sample metadata from filenames, and converts to notame-compatible format
2. **Feature filtering** — sequential pre-correction filters to remove low-quality features
3. **Drift and batch correction** — one or more correction methods run in parallel, each saved to its own output folder
4. **Imputation** — Uses Random Forest to impute any missing values. 
5. **QC metrics and comparison** — per-method quality metrics computed and compared in a summary table
6. **Output** — feature tables, annotations, and QC plots per method

## Requirements

Docker (recommended), or R 4.5+ with dependencies managed via `renv`.

## Usage

### Docker

```bash
docker run --rm \
  -v /path/to/data:/data \
  -e IN_XLSX=/data/msdial_export.xlsx \
  -e PROJECT_FOLDER=/data/output \
  -e POLARITY=POS \
  your-image-name
```

## Key parameters

| Variable | Default | Description |
|---|---|---|
| `IN_XLSX` | — | Path to MSDIAL alignment export (.xlsx) |
| `PROJECT_FOLDER` | — | Root output directory |
| `POLARITY` | `POS` | Ionisation polarity (POS / NEG), used to namespace outputs |
| `CORRECTION_METHODS` | `none,notame` | Comma-separated list of methods to run (see below) |
| `QC_DETECTION_LIMIT` | `0.60` | Min detection rate in QC samples |
| `SAMPLE_DETECTION_LIMIT` | `0.20` | Min detection rate in biological samples |
| `FILL_FILTER` | `0.10` | Min MSDIAL Fill % (alignment confidence) |
| `QC_RSD_FILTER` | `0.80` | Max pre-correction QC RSD; feature must pass in at least one batch |
| `LOW_INT_FILTER_FRAC` | `0.10` | Low-intensity cutoff as fraction of mean p80 intensity |
| `BLANK_RATIO` | `1` | Blank filter ratio; set to `none` to disable |
| `NORMALIZATION` | `none` | Post-correction normalisation (`none` / `pqn`) |
| `N_CORES` | all - 1 | Number of CPU cores for parallelisation |
| `RUV_K` | `3` | Number of unwanted variation factors (notame method only) |

## Correction methods

| Method | Description |
|---|---|
| `none` | Imputation only (baseline, no correction) |
| `notame` | Per-batch cubic spline drift correction followed by RUV-S batch correction using pooled QC samples. Batch correction is skipped when only one batch is present. Same methodology as described in the original notame [paper](https://www.mdpi.com/2218-1989/10/4/135) |
| `pmp_qcrsc` | QC-RSC (Quality Control-Robust Spline Correction) from the [pmp](https://bioconductor.org/packages/pmp/) package. Fits a smoothing spline through QC samples within each batch to correct signal drift. |
| `serrf` | Prototype implementation of SERRF (Systematic Error Removal using Random Forest). Per-feature random forest models trained on QC samples to correct systematic error. Adapted from [Fan et al., Analytical Chemistry 2019](https://doi.org/10.1021/acs.analchem.8b05592). |
| `batchcorr` | Cluster-based spline drift correction followed by between-batch normalisation using the [batchCorr](https://link.springer.com/article/10.1007/s11306-016-1124-4) package (Brunius et al.). |
| `combat_only` | ComBat batch correction only (no drift correction). |
| `loess_combat` | Per-batch LOESS drift correction followed by ComBat batch correction. |
| `waveica` | WaveICA 2.0 — wavelet-based correction for both drift and batch effects (https://link.springer.com/article/10.1007/s11306-021-01839-7). |

## Output structure

```
output/
  {POLARITY}/
    sample_metadata.csv
    correction_comparison.csv
    {method}/
      feature_table.xlsx           # clustered features
      feature_table_full.xlsx      # all features (unclustered)
      annotations.xlsx             # MSDIAL annotations for clustered features
      annotations_full.xlsx        # MSDIAL annotations for all features
      feature_info.xlsx            # feature metadata with QC metrics
      feature_table_rsd30.xlsx     # QC RSD < 30% filter applied
      annotations_rsd30.xlsx
      QC_plots/
intermediates/
  {POLARITY}/
    notame_rev.xlsx                # notame-formatted input
    run_parameters.txt             # all parameter values used
```

## Feature filtering

Applied before correction, in order:

1. Blank filter — removes features where sample signal ≤ `BLANK_RATIO` × blank signal
2. Low-intensity filter — removes features below `LOW_INT_FILTER_FRAC` × mean p80 intensity
3. Fill % filter — removes features with MSDIAL alignment confidence below `FILL_FILTER`
4. QC detection — removes features not detected in ≥ `QC_DETECTION_LIMIT` of QC samples
5. Sample detection — removes features not detected in ≥ `SAMPLE_DETECTION_LIMIT` of biological samples
6. Zero variance — removes features with no variation across samples
7. QC-RSD filter — removes features with QC CV > `QC_RSD_FILTER` in all batches

# notame-workflow-cmsi

Post-MSDIAL preprocessing pipeline for untargeted LC-MS metabolomics data. Converts MSDIAL alignment exports to a standardised format, applies feature quality filters, and evaluates multiple drift and batch correction strategies in parallel.

Built around the [notame](https://github.com/antonvsdata/notame) R package.

## Overview

1. **MSDIAL conversion** — parses MSDIAL alignment exports, extracts sample metadata from filenames, and converts to notame-compatible format. Sample filenames must follow the format: `DATE_BATCH_COLUMN_POLARITY_SAMPLENAME_INJECTIONNUMBER`
2. **Feature filtering** — sequential pre-correction filters to remove low-quality features
3. **Drift and batch correction** — one or more correction methods run in parallel, each saved to its own output folder
4. **Imputation** — random forest imputation of remaining missing values
5. **QC metrics and comparison** — per-method quality metrics computed and compared in a summary table
6. **Feature clustering** — correlated features (e.g. isotopes, adducts) are grouped and compressed to one representative per cluster using notame's `cluster_features` / `compress_clusters`. Both the full unclustered and clustered outputs are retained.
7. **Output** — feature tables, MSDIAL annotations, and QC plots per method

## Requirements

Docker (recommended), or R 4.5+ with dependencies managed via `renv`.

## Usage

### Docker

`IN_XLSX`, `PROJECT_FOLDER`, `COLUMN`, and `POLARITY` are all required. All other parameters are optional with defaults.

```bash
docker run --rm \
  -v /path/to/data:/data \
  -v /path/to/output:/processed \
  -e IN_XLSX=/data/msdial_export.xlsx \
  -e PROJECT_FOLDER=/processed \
  -e COLUMN=RP \
  -e POLARITY=POS \
  -e CORRECTION_METHODS="pmp_qcrsc,notame" \
  your-image-name
```

Output folders are namespaced by `{COLUMN}_{POLARITY}` (e.g. `RP_POS`, `HILIC_NEG`), so multiple modes from the same project can share a single `PROJECT_FOLDER`.

### RStudio (via renv)

Dependencies are managed with `renv`. To restore the environment:

1. Open the project in RStudio (`notame-workflow-cmsi.Rproj` or just open the folder)
2. Install renv if not already available:
   ```r
   install.packages("renv")
   ```
3. Restore the package library:
   ```r
   renv::restore()
   ```
4. Set parameters as environment variables and source the script:
   ```r
   Sys.setenv(
     IN_XLSX          = "C:/path/to/msdial_export.xlsx",
     PROJECT_FOLDER   = "C:/path/to/output",
     COLUMN           = "RP",
     POLARITY         = "POS",
     CORRECTION_METHODS = "pmp_qcrsc,notame"
   )
   source("notame-workflow.r")
   ```
   Or set them in your `.Renviron` file for persistence across sessions.

5. Alternatively, run from the terminal:
   ```bash
   IN_XLSX=... PROJECT_FOLDER=... COLUMN=RP POLARITY=POS Rscript notame-workflow.r
   ```

### Help

```bash
Rscript notame-workflow.r --help
```

## Key parameters

| Variable | Required | Default | Description |
|---|---|---|---|
| `IN_XLSX` | Yes | — | Path to MSDIAL alignment export (.xlsx) |
| `PROJECT_FOLDER` | Yes | — | Root output directory |
| `COLUMN` | Yes | — | Chromatographic column type (e.g. `RP`, `HILIC`) |
| `POLARITY` | Yes | — | Ionisation polarity (`POS` / `NEG`) |
| `CORRECTION_METHODS` | No | `none,notame` | Comma-separated list of methods to run (see below) |
| `QC_DETECTION_LIMIT` | No | `0.60` | Min detection rate in QC samples |
| `SAMPLE_DETECTION_LIMIT` | No | `0.20` | Min detection rate in biological samples |
| `FILL_FILTER` | No | `0.10` | Min MSDIAL Fill % (alignment confidence, 0–1) |
| `QC_RSD_FILTER` | No | `none` | Max pre-correction QC RSD; feature must pass in ≥ 1 batch. Set to e.g. `0.80` to enable |
| `LOW_INT_FILTER_FRAC` | No | `0.10` | Low-intensity cutoff as fraction of mean p80 intensity |
| `BLANK_RATIO` | No | `1` | Blank filter ratio; set to `none` to disable |
| `NORMALIZATION` | No | `none` | Post-correction normalisation (`none` / `pqn`) |
| `N_CORES` | No | all - 1 | Number of CPU cores for parallelisation |
| `RUV_K` | No | `3` | Unwanted variation factors for RUV (notame method only) |

## Correction methods

| Method | Description |
|---|---|
| `none` | Imputation only (no correction; baseline) |
| `notame` | Per-batch cubic spline drift correction followed by RUV-S batch correction using pooled QC samples. Batch correction is skipped when only one batch is present. Described in the original notame [paper](https://www.mdpi.com/2218-1989/10/4/135). |
| `pmp_qcrsc` | QC-RSC (Quality Control-Robust Spline Correction) from the [pmp](https://bioconductor.org/packages/pmp/) package. Fits a smoothing spline through QC samples within each batch to correct signal drift. |
| `serrf` | SERRF (Systematic Error Removal using Random Forest). Per-feature random forest models trained on QC samples to correct systematic error. Adapted from [Fan et al., Analytical Chemistry 2019](https://doi.org/10.1021/acs.analchem.8b05592). |
| `batchcorr` | Cluster-based spline drift correction followed by between-batch normalisation using the [batchCorr](https://link.springer.com/article/10.1007/s11306-016-1124-4) package (Brunius et al.). |
| `combat_only` | ComBat batch correction only (no drift correction). |
| `loess_combat` | Per-batch LOESS drift correction followed by ComBat batch correction. |
| `waveica` | WaveICA 2.0 — wavelet-based correction for both drift and batch effects ([Deng et al. 2021](https://link.springer.com/article/10.1007/s11306-021-01839-7)). |

## Output structure

```
output/
  {COLUMN}_{POLARITY}/              e.g. RP_POS, HILIC_NEG
    sample_metadata.csv
    correction_comparison.csv
    raw_reference.csv
    pre_correction/
      QC_plots/
    {method}/
      feature_table_full.xlsx       # all features, unclustered
      feature_info_full.xlsx        # feature metadata with QC metrics, unclustered
      annotations_full.xlsx         # MSDIAL annotations, unclustered
      feature_table.xlsx            # one representative per cluster
      feature_info.xlsx
      annotations.xlsx              # annotations with cluster-member fallback
      feature_table_full_rsd30.xlsx # above, with QC RSD < 30% filter applied
      feature_table_rsd30.xlsx
      annotations_rsd30.xlsx
      feature_table_full_batchrsd30.xlsx  # above, with per-batch QC RSD < 30% in >= 50% of batches
      feature_table_batchrsd30.xlsx
      annotations_batchrsd30.xlsx
      QC_plots/
intermediates/
  {COLUMN}_{POLARITY}/
    notame_rev.xlsx                 # notame-formatted intermediate
    prefilter_log.csv               # feature counts after each filter step
    run_parameters.txt              # all parameter values used
```

## Feature filtering

Applied before correction, in order:

1. Blank filter — removes features where sample signal ≤ `BLANK_RATIO` × blank signal
2. Low-intensity filter — removes features whose p80 intensity is below `LOW_INT_FILTER_FRAC` × mean p80 across all features
3. Fill % filter — removes features with MSDIAL alignment confidence below `FILL_FILTER`
4. QC detection — removes features not detected in ≥ `QC_DETECTION_LIMIT` of QC samples
5. Sample detection — removes features not detected in ≥ `SAMPLE_DETECTION_LIMIT` of biological samples
6. Zero variance — removes features with no variation across samples
7. QC-RSD filter — removes features with QC CV > `QC_RSD_FILTER` in all batches (disabled by default)

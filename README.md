# notame-workflow-cmsi

Preprocessing pipeline for untargeted LC-MS metabolomics data. Converts a MSDIAL alignment export or an XCMS-based feature table + sample sheet to a standardised format, applies feature quality filters, and evaluates multiple drift and batch correction strategies in parallel.

Built around the [notame](https://github.com/antonvsdata/notame) R package.

## Overview

1. **Conversion** — either parses a MSDIAL alignment export, extracting sample metadata from filenames (format: `DATE_BATCH_COLUMN_POLARITY_SAMPLENAME_INJECTIONNUMBER`), or reads an XCMS-based feature table + sample sheet, which already carries structured sample metadata — see [XCMS input](#xcms-input) below. Either way, the result is converted to notame-compatible format.
2. **Feature filtering** — sequential pre-correction filters to remove low-quality features
3. **Drift and batch correction** — one or more correction methods run in parallel, each saved to its own output folder
4. **Imputation** — two-step: LoD/2 (half-minimum per batch) fills missing values before correction algorithms that require a complete matrix; random forest imputation is then applied after full correction on originally-missing positions only
5. **QC metrics and comparison** — per-method quality metrics computed and compared in a summary table
6. **Feature clustering** — correlated features (e.g. isotopes, adducts) are grouped and compressed to one representative per cluster using notame's `cluster_features` / `compress_clusters`. Both the full unclustered and clustered outputs are retained.
7. **Output** — feature tables, source-pipeline annotations, and QC plots per method

## Requirements

Docker (recommended), or R 4.5+ with dependencies managed via `renv`.

## Usage

### Docker

`PROJECT_FOLDER`, `COLUMN`, and `POLARITY` are required, plus one input: either `IN_XLSX` (MSDIAL), or `IN_FEATURE_TABLE` + `IN_SAMPLE_SHEET` (XCMS — see [XCMS input](#xcms-input)). All other parameters are optional with defaults.

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

**Setup**

1. Double-click `notame-workflow-cmsi.Rproj` to open the project in RStudio.
2. RStudio may prompt you to install `renv` — accept, or run this yourself in the Console:
   ```r
   install.packages("renv")
   ```
3. Restore all required packages (this may take a few minutes the first time):
   ```r
   renv::restore()
   ```

**Running the pipeline**

4. In the RStudio Console, fill in your paths and settings and run:
   ```r
   Sys.setenv(
     IN_XLSX          = "C:/path/to/msdial_export.xlsx",  # your MSDIAL alignment export
     PROJECT_FOLDER   = "C:/path/to/output",              # where results will be saved
     COLUMN           = "RP",                             # chromatographic column (e.g. RP, HILIC)
     POLARITY         = "POS",                            # ionisation polarity: POS or NEG
     CORRECTION_METHODS = "pmp_qcrsc,notame"              # see Correction methods table below
   )
   source("notame-workflow.r")
   ```
   See Key parameters below for all available options.
   Results will appear in `PROJECT_FOLDER` under a subfolder named `{COLUMN}_{POLARITY}` (e.g. `RP_POS`).

### Help

```bash
Rscript notame-workflow.r --help
```

## Configuration file

Any parameter below can also be supplied via a YAML config file instead of
an environment variable — handy for avoiding long `docker run -e ...` chains
or checking a run's settings into version control. Point `CONFIG_FILE` at
the file:

```bash
docker run --rm \
  -v /path/to/data:/data \
  -v /path/to/output:/processed \
  -e CONFIG_FILE=/data/config.yaml \
  your-image-name
```

**Precedence:** env vars win over the config file when both set the same
parameter, so one-off Docker overrides keep working unchanged on top of a
shared config file. `IN_XLSX`, `PROJECT_FOLDER`, `COLUMN`, and `POLARITY`
can all come from the config file too — a config file alone is enough to
run the pipeline.

See [`config.example.yaml`](config.example.yaml) for a full worked example.

```yaml
IN_XLSX: /data/msdial_export.xlsx
PROJECT_FOLDER: /processed
COLUMN: RP
POLARITY: POS
CORRECTION_METHODS: pmp_qcrsc,notame
QC_DETECTION_LIMIT: 0.60
```

### Sample-type classification override

Sample type (QC / ltQC / Blank / Wash / Cond / SST / Sample) is normally
inferred from filename keywords (see [Sample types](#sample-types) below).
Labs using different naming conventions can override this via
`sample_type_rules` in the config file — an ordered list of `pattern`/`type`
pairs, tried top-to-bottom against the sample name (case-insensitive), first
match wins, anything unmatched falls back to `Sample`. This list **replaces**
the built-in defaults entirely when given, rather than extending them.
(This applies to the MSDIAL input path only — see [XCMS input](#xcms-input)
for how sample type is determined for XCMS-based input.)

```yaml
sample_type_rules:
  - pattern: "PoolQC"
    type: QC
  - pattern: "LongTermQC"
    type: ltQC
  - pattern: "^Blank_"
    type: Blank
```

## XCMS input

As an alternative to a MSDIAL alignment export, the pipeline accepts output
from an XCMS-based pipeline: a feature table CSV (`IN_FEATURE_TABLE`) plus a
sample sheet XLSX (`IN_SAMPLE_SHEET`). Set both (env var or config file) instead
of `IN_XLSX` — exactly one of the two input methods must be given.

```bash
docker run --rm \
  -v /path/to/data:/data \
  -v /path/to/output:/processed \
  -e IN_FEATURE_TABLE=/data/feature_table.csv \
  -e IN_SAMPLE_SHEET=/data/sample_sheet.xlsx \
  -e PROJECT_FOLDER=/processed \
  -e COLUMN=RP \
  -e POLARITY=NEG \
  your-image-name
```

Unlike MSDIAL exports, this pipeline already carries structured sample
metadata (no filename parsing needed) and the two files are joined by
**`sample_label`** — every abundance column in the feature table must be
headed by a value from the sample sheet's `sample_label` column.

**`sample_label` must be unique** across the samples being processed (it's
the join key, and becomes `Original_name` in the output as-is). If your
naming convention reuses a plate-relative label per batch (e.g. `sQC01` in
every batch), include the batch in the label itself (e.g. `B10W22-sQC01`)
so it's unique — a duplicate is rejected with a clear error rather than
silently dropping the extra samples.

**`IN_FEATURE_TABLE`** (csv) must have: `feature`, `mzmed`, `rtmed`,
`npeaks`, plus one abundance column per sample (headed by `sample_label`).
`mzmin`/`mzmax`/`rtmin`/`rtmax`/`ms_level` and any per-type detection-count
columns are carried through into the output but not required. `rtmed`/
`rtmin`/`rtmax` are assumed to be in **seconds** (XCMS's own convention) and
are converted to minutes to match notame's `Average_Rt_min`.

**`IN_SAMPLE_SHEET`** (xlsx) must have: `batch_plate`, `column`, `polarity`,
`sample_label`, `sample_type`, `injection_order`, `filename`, `include`.
`batch_plate` (not `batch`) is used as the batch grouping for drift/batch
correction, since a nominal batch can span multiple plates. Rows are
filtered to `include == TRUE` and to the requested `COLUMN`/`POLARITY`
(case-insensitive); the run fails with a clear error if nothing matches. A
row with `needs_review == TRUE` that's still included triggers a warning
rather than being dropped.

`sample_type` is mapped to the pipeline's internal QC vocabulary via a
built-in table (not config-overridable, since this field is already an
authoritative typed value rather than a filename guess):

| sample_type | → | Internal type |
|---|---|---|
| `Sample` | → | `Sample` |
| `sQC`, `QC` | → | `QC` |
| `ltQC` | → | `ltQC` |
| `Blank` | → | `Blank` |
| `MatrixBlank` | → | `MatrixBlank` |
| `Wash` | → | `Wash` |
| `Cond` | → | `Cond` |
| `SST` | → | `SST` |
| anything else | → | `Sample` |

Since this pipeline doesn't produce metabolite annotations, the
`Adduct_type`/`Metabolite_name` output columns are left blank for XCMS runs;
`Fill_pct` is computed from `npeaks` (features actually detected) as a
directly analogous substitute.

## Key parameters

| Variable | Required | Default | Description |
|---|---|---|---|
| `CONFIG_FILE` | No | — | Path to an optional YAML config file supplying defaults for any parameter in this table (env vars still take precedence — see [Configuration file](#configuration-file)) |
| `IN_XLSX` | Yes\* | — | Path to MSDIAL alignment export (.xlsx). \*Required unless `IN_FEATURE_TABLE`+`IN_SAMPLE_SHEET` are set instead |
| `IN_FEATURE_TABLE` | Yes\* | — | Path to an XCMS-based feature table (.csv) — alternative to `IN_XLSX`, see [XCMS input](#xcms-input) |
| `IN_SAMPLE_SHEET` | Yes\* | — | Path to the XCMS pipeline's sample sheet (.xlsx) — required together with `IN_FEATURE_TABLE` |
| `PROJECT_FOLDER` | Yes | — | Root output directory |
| `FORCE_RECONVERT` | No | `FALSE` | Force re-running the conversion step even if a cached result from a previous run is found |
| `COLUMN` | Yes | — | Chromatographic column type (e.g. `RP`, `HILIC`) |
| `POLARITY` | Yes | — | Ionisation polarity (`POS` / `NEG`) |
| `CORRECTION_METHODS` | No | `none,notame` | Comma-separated list of methods to run (see below) |
| `QC_DETECTION_LIMIT` | No | `0.60` | Min detection rate in QC samples |
| `SAMPLE_DETECTION_LIMIT` | No | `0.20` | Min detection rate in biological samples |
| `MIN_QC_SAMPLE_DETECTION` | No | `0.50` | Min fraction of features detected in a QC or ltQC sample for it to be used as reference. Samples below this are removed before processing (e.g. empty injections) |
| `MIN_BATCH_DETECTION` | No | `1` | Min number of detections a feature must have in every batch. Features absent from any entire batch are removed (set to `0` to disable) |
| `QC_RSD_FILTER` | No | `none` | Max pre-correction QC RSD (robust: MAD/median); feature must pass in ≥ 50% of batches. Set to e.g. `0.80` to enable |
| `RSD_THRESHOLD` | No | `0.30` | RSD threshold for post-correction output filtering. Controls both global (`_rsdXX`) and per-batch (`_batchrsdXX`) filtered outputs. Suffix reflects the threshold, e.g. `_rsd40` if set to `0.40` |
| `LOW_INT_FILTER_FRAC` | No | `0.10` | Low-intensity cutoff as fraction of mean p80 intensity (overridden by `LOW_INT_FILTER`) |
| `LOW_INT_FILTER` | No | `none` | Absolute low-intensity cutoff (overrides `LOW_INT_FILTER_FRAC`) |
| `LOW_INT_PERCENTILE` | No | `0.80` | Percentile used for the low-intensity filter |
| `BLANK_RATIO` | No | `none` | Blank filter ratio — removes features where mean(Sample) ≤ `BLANK_RATIO` × mean(SolvBlank). Set to e.g. `1` to enable |
| `NORMALIZATION` | No | `none` | Post-correction normalisation (`none` / `pqn`). See below |
| `LOESS_SPAN` | No | `0.75` | LOESS smoothing span for QC-based drift correction (`loess_combat`, `loess_limma`, `loess_feature_median`, `loess_global_median`). Higher = smoother, more conservative |
| `LOESS_SAMPLE_SPAN` | No | `0.9` | LOESS smoothing span for QC-free drift correction (`loess_samples_combat`, `loess_samples_limma`), fit on biological samples instead of QC. Wider than `LOESS_SPAN` by default since sample points are far noisier |
| `LOESS_SAMPLE_MIN_OBS` | No | `10` | Min finite sample observations per feature required to attempt QC-free drift correction (`loess_samples_combat`, `loess_samples_limma`). Higher than the QC-based fit's threshold of 4 |
| `LOESS_MIN_QC_PER_BATCH` | No | `4` | Min QC samples a batch needs to use QC-based drift correction in `loess_samples_combat` / `loess_samples_limma`. Batches at or above this use QC-based LOESS (`LOESS_SPAN`); batches below it trial the QC-free fit on samples (`LOESS_SAMPLE_SPAN`, `LOESS_SAMPLE_MIN_OBS`), validated per `LOESS_MIN_LTQC_VALIDATE` |
| `LOESS_MIN_LTQC_VALIDATE` | No | `3` | Min ltQC samples a batch needs to validate the QC-free trial correction in `loess_samples_combat` / `loess_samples_limma` (batches below `LOESS_MIN_QC_PER_BATCH` only). The trial is kept if it improves the ltQC/Sample D-ratio (`MAD(ltQC)/MAD(Sample)`, lower is better) versus the uncorrected batch, discarded otherwise — D-ratio rather than raw ltQC RSD, since any real drift correction shrinks sample variance somewhat, so it only credits a disproportionate improvement in ltQC relative to Sample. Batches with fewer ltQC than this are left uncorrected — there's no way to validate the trial |
| `CORDBAT_REF_BATCH` | No | auto | Reference batch ID for CordBat methods. All other batches are corrected onto this batch. Defaults to auto-selecting the batch with the lowest median feature RSD |
| `N_CORES` | No | all - 1 | Number of CPU cores for parallelisation |
| `RUV_K` | No | `3` | Unwanted variation factors for RUV (notame method only) |

## Correction methods

| Method | Description | Parameters |
|---|---|---|
| `none` | Imputation only (no correction; baseline) | — |
| `notame` | Per-batch cubic spline drift correction followed by RUV-S batch correction using pooled QC samples. Batch correction is skipped when only one batch is present. Described in the original notame [paper](https://www.mdpi.com/2218-1989/10/4/135). | `RUV_K` |
| `pmp_qcrsc` | QC-RSC (Quality Control-Robust Spline Correction) from the [pmp](https://bioconductor.org/packages/pmp/) package. Fits a smoothing spline through QC samples within each batch to correct signal drift. | — |
| `pmp_qcrsc_scale` | As `pmp_qcrsc`, plus global median scaling for any batch with fewer than 4 QC samples (which pmp cannot spline-correct); pmp-corrected batches are left untouched. | — |
| `pmp_qcrsc_feature_scale` | As `pmp_qcrsc_scale`, but uses per-feature median scaling for the no-QC batches instead of a single global factor, consistent with pmp's own feature-wise alignment. | — |
| `serrf` | SERRF (Systematic Error Removal using Random Forest). Per-feature random forest models trained on QC samples to correct systematic error. Adapted from [Fan et al., Analytical Chemistry 2019](https://doi.org/10.1021/acs.analchem.8b05592). | `SERRF_NUM` |
| `batchcorr` | Cluster-based spline drift correction followed by between-batch normalisation using the [batchCorr](https://link.springer.com/article/10.1007/s11306-016-1124-4) package (Brunius et al.). | — |
| `combat_only` | ComBat batch correction only (no drift correction). | — |
| `loess_combat` | Per-batch LOESS drift correction (QC-based) followed by ComBat batch correction. | `LOESS_SPAN` |
| `loess_samples_combat` | Per-batch LOESS drift correction, chosen per batch in three tiers: (1) QC-based (`LOESS_SPAN`) if the batch has at least `LOESS_MIN_QC_PER_BATCH` QC samples; (2) otherwise, a QC-free trial fit on biological samples (`LOESS_SAMPLE_SPAN`, `LOESS_SAMPLE_MIN_OBS`), kept only if it measurably improves the ltQC/Sample D-ratio (`MAD(ltQC)/MAD(Sample)`, lower is better) in that batch versus leaving it uncorrected — a genuine held-out check, since ltQC is never used to fit the trial, and a ratio rather than raw ltQC RSD since any real drift correction shrinks sample variance somewhat (only a *disproportionate* shrink relative to ltQC is penalized) — provided the batch has at least `LOESS_MIN_LTQC_VALIDATE` ltQC samples; (3) otherwise the batch is left uncorrected, since there's no QC or ltQC data to justify or validate a correction. Followed by ComBat batch correction. QC-based fitting is preferred whenever there's enough QC to support it — it doesn't risk removing real biological signal along with drift the way fitting on samples can; the QC-free trial's own fit uses a robust family (`family = "symmetric"`) and a wider span to guard against fitting individual-sample noise as drift. | `LOESS_SPAN`, `LOESS_SAMPLE_SPAN`, `LOESS_SAMPLE_MIN_OBS`, `LOESS_MIN_QC_PER_BATCH`, `LOESS_MIN_LTQC_VALIDATE` |
| `loess_limma` | Per-batch LOESS drift correction (QC-based) followed by `limma::removeBatchEffect()` for between-batch correction. Appropriate when QC data is partially compromised. | `LOESS_SPAN` |
| `loess_samples_limma` | Same per-batch three-tier choice as `loess_samples_combat`, but with `limma::removeBatchEffect()` instead of ComBat for between-batch correction. | `LOESS_SPAN`, `LOESS_SAMPLE_SPAN`, `LOESS_SAMPLE_MIN_OBS`, `LOESS_MIN_QC_PER_BATCH`, `LOESS_MIN_LTQC_VALIDATE` |
| `loess_feature_median` | Per-batch LOESS drift correction (QC-based) followed by per-feature median ratio normalisation. Scales each batch so its biological sample median per feature matches the grand median. More flexible than global scaling but noisier for sparse features. QC-independent. | `LOESS_SPAN` |
| `loess_global_median` | Per-batch LOESS drift correction (QC-based) followed by global median ratio normalisation. Computes one scaling factor per batch from the median of all biological sample intensities and applies it uniformly to all features. Assumes a constant multiplicative offset per batch. QC-independent. | `LOESS_SPAN` |
| `cordbat_only` | CordBat batch correction only (no drift correction). Uses a Gaussian Graphical Model built from correlated feature communities to learn per-feature scale and offset parameters. Requires a reference batch (auto-selected by default). | `CORDBAT_REF_BATCH` |
| `loess_cordbat` | Per-batch LOESS drift correction, always QC-free (fit on biological samples, `LOESS_SAMPLE_SPAN`/`LOESS_SAMPLE_MIN_OBS`) — deliberately, since CordBat's own between-batch step is also fit on biological samples rather than QC, so the combined method stays QC-free end to end rather than mixing a QC-anchored drift step with a QC-free batch-correction step. Followed by CordBat batch correction. Requires a reference batch (auto-selected by default, from QC-quality where available). | `LOESS_SAMPLE_SPAN`, `LOESS_SAMPLE_MIN_OBS`, `CORDBAT_REF_BATCH` |
| `waveica` | WaveICA 2.0 — wavelet-based correction for both drift and batch effects, QC-independent ([Deng et al. 2021](https://link.springer.com/article/10.1007/s11306-021-01839-7)). | — |

## Normalisation

PQN (Probabilistic Quotient Normalisation) can be applied after drift/batch correction by setting `NORMALIZATION=pqn`. It scales each sample by the median ratio of its intensities to a reference spectrum. The reference is the median spectrum of pooled QC samples when QC samples are present, or the median of all biological samples otherwise.

PQN corrects for differences in overall sample concentration or dilution, and is applied independently within each correction method branch.

## Output structure

```
output/
  {COLUMN}_{POLARITY}/              e.g. RP_POS, HILIC_NEG
    sample_metadata.csv
    method_comparison.csv           # cross-method QC comparison, best-first
    raw_reference.csv
    report.html                     # single-page run summary — start here
    pre_correction/
      QC_plots/
    {method}/
      results_full.xlsx                       # all features, unclustered
      results_clustered.xlsx                  # one representative per cluster
      results_full_rsdXX.xlsx                 # above, with global QC RSD < XX% filter applied
      results_clustered_rsdXX.xlsx
      results_full_batchrsdXX.xlsx            # above, with per-batch QC RSD < XX% in >= 50% of batches
      results_clustered_batchrsdXX.xlsx
      batch_summary_post_correction.csv
      QC_plots/
intermediates/
  {COLUMN}_{POLARITY}/
    notame_rev.xlsx                 # notame-formatted intermediate
    conversion_annotations.rds      # cached conversion result (mode + annotations)
    prefilter_log.csv               # feature counts after each filter step
    batch_summary.csv               # pre-correction per-batch missingness/QC-RSD
    run_log.csv                     # per-method success/failure, timing, feature count
    run_parameters.txt              # all parameter values used
```

`report.html` is a self-contained, single-page summary of the run — run
parameters, pre-filtering counts, batch quality, the run log, the
cross-method comparison, and the QC plots (pre- and post-correction) — open
it in a browser as a starting point before digging into individual
workbooks.

The `XX` in output filenames reflects `RSD_THRESHOLD` (e.g. `_rsd30` at default, `_rsd40` if `RSD_THRESHOLD=0.40`).

Each `results*.xlsx` workbook contains:

| Sheet | Contents |
|---|---|
| `Peak_table` | Feature abundance matrix (rows = features, columns = samples) |
| `Feature_metadata` | mz, rt, adduct, QC metrics (RSD, D_ratio), and a curated subset of the source pipeline's annotations (metabolite name, fill %, S/N — MSDIAL input only, blank for XCMS input) merged by feature |
| `Cluster_info` | Cluster ID, member features, cluster size, and MPA — only present in `results_clustered*.xlsx` workbooks |
| `Settings` | Run parameters, correction method, feature-set filter applied, and versions of all loaded packages |

## Feature filtering

Applied before correction, in order:

1. Blank filter — removes features where sample signal ≤ `BLANK_RATIO` × blank signal (SolvBlank samples only; disabled by default)
2. Low-intensity filter — removes features whose p80 intensity is below `LOW_INT_FILTER_FRAC` × mean p80 across all features (or below `LOW_INT_FILTER` if set)
3. QC/ltQC sample quality check — removes individual QC or ltQC samples with feature detection rate below `MIN_QC_SAMPLE_DETECTION` (e.g. empty injections, failed runs)
4. QC detection — removes features not detected in ≥ `QC_DETECTION_LIMIT` of QC samples
5. Sample detection — removes features not detected in ≥ `SAMPLE_DETECTION_LIMIT` of biological samples
6. Zero variance — removes features with no variation across samples
7. QC-RSD filter — removes features with QC robust RSD (MAD/median) > `QC_RSD_FILTER` in ≥ 50% of batches (disabled by default)

## Sample types

Sample types are inferred from filenames. The following are recognised:

| Type | Filename pattern | Role |
|---|---|---|
| `Sample` | (default) | Biological samples |
| `QC` | contains `sQC` | Pooled QC samples used for drift/batch correction |
| `ltQC` | contains `ltQC` | Long-term QC samples used as held-out validation only |
| `Blank` | contains `SolvBlank` | Solvent blanks used for blank filtering |
| `MatrixBlank` | contains `blank` (not SolvBlank) | Matrix blanks — excluded from processing |
| `Wash` | contains `MeOH` | Column wash injections — excluded from processing |
| `Cond` | contains `CondPlasma` | Conditioning injections — excluded from processing |
| `SST` | contains `SST` + digit | System suitability test injections — excluded |
| `MSe` | `MSe` in column/mode field | Data-independent MS² acquisitions — excluded |
| `MS2` | `MS2` in column/mode field | Targeted MS² acquisitions — excluded |

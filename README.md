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
  -e CORRECTION_METHODS="none:none:pmp_qcrsc,notame_spline:qc:ruv_s" \
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
     CORRECTION_METHODS = "none:none:pmp_qcrsc,notame_spline:qc:ruv_s"              # see Correction methods table below
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
CORRECTION_METHODS: none:none:pmp_qcrsc,notame_spline:qc:ruv_s
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
| `CORRECTION_METHODS` | No | `none:none:none,notame_spline:qc:ruv_s` | Comma-separated list of `drift:basis:batch` recipes to run (see [Correction methods](#correction-methods) below) |
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
| `LOESS_QC_SPAN` | No | `0.75` | LOESS smoothing span for `drift=loess`'s QC-based fit (`basis=qc`, or `basis=hybrid`'s QC branch). Higher = smoother, more conservative |
| `LOESS_SAMPLE_SPAN` | No | `0.9` | LOESS smoothing span for `drift=loess`'s samples-based fit (`basis=samples`, or `basis=hybrid`'s samples branch), fit on biological samples instead of QC. Wider than `LOESS_QC_SPAN` by default since sample points are far noisier |
| `DRIFT_SAMPLE_MIN_OBS` | No | `10` | Min finite sample observations per feature required to attempt a samples-based drift fit (`basis=samples`, or `basis=hybrid`'s samples branch; applies to both `drift=loess` and `drift=huber`). Higher than the QC-based fit's threshold of 4 |
| `DRIFT_MIN_QC_PER_BATCH` | No | `4` | Min QC samples a batch needs to use the QC-based fit under `basis=hybrid` (applies to `drift=loess` and `drift=huber` alike). Batches at or above this use the QC-based fit (`LOESS_QC_SPAN`/`HUBER_QC_K`); batches below it trial the samples-based fit (`LOESS_SAMPLE_SPAN`/`HUBER_SAMPLE_K`, `DRIFT_SAMPLE_MIN_OBS`), validated per `DRIFT_MIN_LTQC_VALIDATE`. Also the hard cutoff for `basis=qc` (no samples-based fallback there — batches below this are left uncorrected) |
| `DRIFT_MIN_LTQC_VALIDATE` | No | `3` | Min ltQC samples a batch needs to validate the samples-based trial correction under `basis=hybrid` (batches below `DRIFT_MIN_QC_PER_BATCH` only). The trial is kept if it improves the ltQC/Sample D-ratio (`MAD(ltQC)/MAD(Sample)`, lower is better) versus the uncorrected batch, discarded otherwise — D-ratio rather than raw ltQC RSD, since any real drift correction shrinks sample variance somewhat, so it only credits a disproportionate improvement in ltQC relative to Sample. Batches with fewer ltQC than this are left uncorrected — there's no way to validate the trial |
| `DRIFT_HYBRID_VALIDATE` | No | `TRUE` | Set to `FALSE` to skip ltQC validation entirely under `basis=hybrid`: any batch below `DRIFT_MIN_QC_PER_BATCH` then always gets the samples-based correction, regardless of ltQC availability or outcome. Reintroduces the risk the validation step exists to catch — use deliberately |
| `AUTO_LOESS_SPANS` | No | `0.5,0.75,0.9` | Comma-separated LOESS spans `drift=auto` evaluates as candidates for the QC-based selection (leave-one-out CV on QC) |
| `AUTO_HUBER_KS` | No | `1.0,1.345,2.0` | Comma-separated Huber regression `k` values (`MASS::rlm`, `psi.huber`) `drift=auto` evaluates as candidates for the QC-based selection. Lower = more robust to outlier QC but less statistically efficient; `1.345` is `rlm`'s own default |
| `AUTO_SAMPLE_LOESS_SPANS` | No | `0.3,0.6,0.9` | Comma-separated LOESS spans `drift=auto` evaluates for the samples-based selection (fit on samples, validated against ltQC) |
| `AUTO_SAMPLE_HUBER_KS` | No | `1.0,1.345,2.0` | Comma-separated Huber `k` values for `drift=auto`'s samples-based candidate pool |
| `AUTO_MIN_QC_PER_BATCH` | No | `4` | Min QC samples a batch needs to contribute to (and, under `basis=hybrid`, receive) `drift=auto`'s QC-based candidate selection, rather than its samples-based one |
| `AUTO_MIN_LTQC_VALIDATE` | No | `3` | Min ltQC samples a batch needs to contribute to `drift=auto`'s samples-based candidate selection. Under `basis=samples`, once a winner is chosen it's applied to every batch regardless of this threshold — it only decides which batches help pick the winner. Under `basis=hybrid`, batches with fewer are left uncorrected |
| `AUTO_MIN_CV_OBS` | No | `4` | Min finite training observations (QC, or samples for the samples-based pool) a feature needs before `drift=auto` attempts to fit any candidate for it |
| `HUBER_QC_K` | No | `1.345` | Huber regression tuning constant (`MASS::rlm`, `psi.huber`) for `drift=huber`'s QC-based fit (`basis=qc`, or `basis=hybrid`'s QC branch). Fixed, not auto-searched — see `HUBER_QC_CV_KS` for per-feature CV selection, or `AUTO_HUBER_KS`/`drift=auto` for a dataset-wide CV-chosen `k` instead |
| `HUBER_SAMPLE_K` | No | `1.345` | Huber tuning constant for `drift=huber`'s samples-based fit (`basis=samples`, or `basis=hybrid`'s samples branch) |
| `SVA_N_SV` | No | auto | Number of surrogate variables SVA estimates for `batch=sva`'s correction step (`sva::num.sv`, Buja-Eyuboglu permutation test). `0` disables SVA (Batch-only correction via `limma::removeBatchEffect`) |
| `LOESS_QC_CV_SPANS` | No | (disabled) | Comma-separated LOESS span candidates for `drift=loess`'s QC-based step (`basis=qc`, or `basis=hybrid`'s QC branch). When set, span is chosen per **feature** via leave-one-out CV on QC (mirroring `notame::correct_drift()`'s own per-feature `smooth.spline()` parameter selection) instead of one shared `LOESS_QC_SPAN`. Safe per-feature only because it's evaluated against QC (technical replicates); not offered for `basis=samples` |
| `HUBER_QC_CV_KS` | No | (disabled) | Same idea as `LOESS_QC_CV_SPANS`, for `drift=huber`'s QC-based step |
| `CORDBAT_REF_BATCH` | No | auto | Reference batch ID for `batch=cordbat`. All other batches are corrected onto this batch. Defaults to auto-selecting the batch with the lowest median feature RSD |
| `WAVEICA_ALPHA` | No | `0.05` | Comma-separated trade-off value(s) (0-1) for `batch=waveica` (WaveICA2.0)'s internal ICA step: 0 = spatial ICA, 1 = temporal ICA. **Not** a significance/flagging threshold — it doesn't itself control how many components get removed (see `WAVEICA_CUTOFF`). A single value is fixed; multiple values trigger a search (see below) |
| `WAVEICA_CUTOFF` | No | `0.10` | Comma-separated threshold(s) (0-1): the minimum R² (vs. injection order) an individual ICA component must reach to be subtracted out as technical noise. This is the actual aggressiveness dial — every wavelet level is always ICA-decomposed regardless; `Cutoff` decides which resulting *components* get removed, not which levels get analyzed. Lower = more aggressive (more components removed); higher = more conservative. Single value fixed, multiple values searched |
| `WAVEICA_K` | No | `auto` (`2 x n_batches`) | Comma-separated number(s) of independent components to decompose into. Each entry is a number or `auto` (`2 x n_batches`, resolved per-run). Single value fixed, multiple values searched |
| `WAVEICA_WF` | No | `haar` | Wavelet family for `batch=waveica`. Not searched (kept fixed — see [Correction methods](#correction-methods)) |
| `WAVEICA_EVAL_GROUP` | No | `ltQC` | Which group (`ltQC` or `QC`) the `WAVEICA_ALPHA`/`WAVEICA_CUTOFF`/`WAVEICA_K` search evaluates candidates against (D-ratio vs. Sample). WaveICA2.0 never fits on QC or ltQC — it corrects using only injection order — so either is a genuine held-out reference; `QC` is worth trying if it has more samples than ltQC in your data |
| `WAVEICA_V1_WF` | No | `haar` | Wavelet family for `batch=waveica_v1` (the original WaveICA — separate setting from `WAVEICA_WF`, different package) |
| `WAVEICA_V1_K` | No | `20` | Max components `batch=waveica_v1`'s ICA step decomposes into |
| `WAVEICA_V1_T` | No | `0.05` | Threshold (0-1) for considering a component associated with batch in `batch=waveica_v1` — tested against real batch labels directly, unlike `WAVEICA_CUTOFF`'s injection-order proxy |
| `WAVEICA_V1_T2` | No | `0.05` | Threshold (0-1) for considering a component associated with a biological comparison group in `batch=waveica_v1`. Currently inert — this pipeline has no biological-group column to supply `waveica_v1`'s optional `group` argument |
| `WAVEICA_V1_ALPHA` | No | `0` | Trade-off (0-1) between sample-wise and variable-wise independence in `batch=waveica_v1`'s ICA step. The same *kind* of parameter as `WAVEICA_ALPHA` (both are ICA spatial/temporal trade-offs) — kept as a separate setting since they're independent packages with separately-tuned defaults, not because the concept differs |
| `COMBAT_MEAN_ONLY` | No | `auto` | Whether ComBat (`batch=combat`, with any drift method) adjusts only each feature's per-batch mean (`TRUE`) or also forces every batch's variance to match a common value (`FALSE`, ComBat's own default). Forcing variance equal across batches is the usual cause of PCA looking artificially "flattened" after correction when batches genuinely differ in spread. `auto` tries both and keeps whichever gives the better ltQC/Sample D-ratio; set `TRUE`/`FALSE` to force a specific behaviour |
| `COMBAT_PAR_PRIOR` | No | `auto` | Whether ComBat's empirical Bayes prior is estimated parametrically (`TRUE`, assumes a Normal/Inverse-Gamma shape — faster) or non-parametrically (`FALSE` — slower, more robust to non-Gaussian batch effects). `auto` tries both and keeps whichever gives the better ltQC/Sample D-ratio, same mechanism as `COMBAT_MEAN_ONLY` |
| `N_CORES` | No | all - 1 | Number of CPU cores for parallelisation |
| `RUV_K` | No | `3` | Unwanted variation factors for `batch=ruv_s` (notame's RUV-S) |

## Correction methods

Each `CORRECTION_METHODS` entry is a `drift:basis:batch` recipe — three independent choices, not one opaque method name — so any combination below is directly reachable without new code. Each recipe gets its own output subfolder (colons become dashes in the folder name; the recipe string itself, colons included, is preserved everywhere else — logs, the `method` column in QC summary CSVs, `method_comparison.csv`).

```
CORRECTION_METHODS="loess:hybrid:combat,huber:qc:feature_median,none:none:cordbat"
```

### 1. Drift method (1st field)

| Value | Description | Parameters |
|---|---|---|
| `none` | No within-batch drift correction. | — |
| `loess` | LOESS drift correction. | `LOESS_QC_SPAN`, `LOESS_SAMPLE_SPAN`, `LOESS_QC_CV_SPANS` |
| `huber` | Huber robust regression (`MASS::rlm`) — a single rigid linear trend instead of a locally flexible curve; more stable than LOESS at small QC counts but can't track curved drift. | `HUBER_QC_K`, `HUBER_SAMPLE_K`, `HUBER_QC_CV_KS` |
| `auto` | Auto-selected via cross-validation from several LOESS spans, several Huber `k` values, and a flat/no-op baseline. Evidence is pooled across all relevant batches before picking one winner (leave-one-out CV on QC for the QC-based selection; held-out ltQC/Sample D-ratio for the samples-based selection) — every batch that gets corrected uses the same method, never a different one per batch. | `AUTO_LOESS_SPANS`, `AUTO_HUBER_KS`, `AUTO_SAMPLE_LOESS_SPANS`, `AUTO_SAMPLE_HUBER_KS`, `AUTO_MIN_QC_PER_BATCH`, `AUTO_MIN_LTQC_VALIDATE`, `AUTO_MIN_CV_OBS` |
| `notame_spline` | notame's own per-feature cubic smoothing spline (`notame::correct_drift()`), which auto-tunes its own smoothness per feature via cross-validation internally. Only supports `basis=qc`. | — |

### 2. Basis (2nd field)

Which data the drift method is fit against. Must be `none` if and only if drift method is `none`.

| Value | Description |
|---|---|
| `qc` | Fit only on QC samples, only for batches with at least `DRIFT_MIN_QC_PER_BATCH` QC samples. Batches below that are left uncorrected — no fallback. |
| `samples` | Fit only on biological samples, for every batch, regardless of QC availability. The samples-based fit uses a robust family (LOESS: `family="symmetric"`) and a wider span/looser `k` by default to guard against fitting individual-sample noise as drift. |
| `hybrid` | QC-based if the batch has at least `DRIFT_MIN_QC_PER_BATCH` QC samples; otherwise a samples-based trial, kept only if it measurably improves the ltQC/Sample D-ratio (`MAD(ltQC)/MAD(Sample)`, lower is better — see `DRIFT_MIN_LTQC_VALIDATE`) versus leaving it uncorrected — a genuine held-out check, since ltQC is never used to fit the trial; otherwise left uncorrected. QC-based fitting is preferred whenever there's enough QC to support it, since it doesn't risk removing real biological signal along with drift the way fitting on samples can. `DRIFT_HYBRID_VALIDATE=FALSE` skips the ltQC check and always keeps the samples-based trial. |

### 3. Batch method (3rd field)

| Value | Description | Parameters |
|---|---|---|
| `none` | No between-batch step. | — |
| `combat` | ComBat batch correction. | `COMBAT_MEAN_ONLY`, `COMBAT_PAR_PRIOR` |
| `sva` | Surrogate Variable Analysis (`sva` package): known `Batch` plus `SVA_N_SV` latent surrogate variables (structure not already explained by `Batch`) are regressed out together via `limma::removeBatchEffect()`. Needs no QC or replicate anchor, only known batch labels — more flexible than ComBat's per-batch mean/variance shift, at the cost of being less directly interpretable. | `SVA_N_SV` |
| `limma` | `limma::removeBatchEffect()`. Appropriate when QC data is partially compromised. | — |
| `feature_median` | Per-feature median ratio normalisation. Scales each batch so its biological sample median per feature matches the grand median. More flexible than global scaling but noisier for sparse features. QC-independent. | — |
| `global_median` | Global median ratio normalisation. Computes one scaling factor per batch from the median of all biological sample intensities and applies it uniformly to all features. Assumes a constant multiplicative offset per batch. QC-independent. | — |
| `ruv_s` | notame's RUV-S, using pooled QC samples. QC-anchored. Described in the original notame [paper](https://www.mdpi.com/2218-1989/10/4/135). | `RUV_K` |
| `cordbat` | CordBat — Gaussian Graphical Model built from correlated feature communities to learn per-feature scale and offset parameters. Requires a reference batch (auto-selected by default, from QC quality where available). Accepts drift-corrected input (e.g. `loess:samples:cordbat`) or raw input (`none:none:cordbat`). | `CORDBAT_REF_BATCH` |
| `batchcorr` | Cluster-based spline drift correction + between-batch normalisation from the [batchCorr](https://link.springer.com/article/10.1007/s11306-016-1124-4) package (Brunius et al.) — couples drift and batch correction internally, so **requires `drift=none`**. | — |
| `waveica` | WaveICA 2.0 — wavelet-based correction for both drift and batch effects, QC-independent ([Deng et al. 2021](https://link.springer.com/article/10.1007/s11306-021-01839-7)). Uses injection order as a proxy for batch structure rather than batch labels directly. Couples drift and batch correction internally, so **requires `drift=none`**. `WAVEICA_ALPHA`/`WAVEICA_CUTOFF`/`WAVEICA_K` each accept a comma-separated list; more than one candidate overall triggers a search over the full cross-product, evaluated against `WAVEICA_EVAL_GROUP`/Sample D-ratio (primary — selects the winner), with PCA-space distance ratio and PERMANOVA R²(Batch) printed alongside every candidate as independent cross-checks (informational only — WaveICA is an ICA-based method operating jointly across features, so a purely per-feature metric like D-ratio could miss damage to that joint structure; the other two catch that). Candidates run in parallel via this pipeline's existing `foreach`/`N_CORES` setup, with `mc.cores` pinned to 1 inside each worker to avoid nested-parallelism oversubscription against WaveICA2.0's own internal `parallel::mclapply()` call. `WAVEICA_WF` is not searched. | `WAVEICA_ALPHA`, `WAVEICA_CUTOFF`, `WAVEICA_K`, `WAVEICA_WF`, `WAVEICA_EVAL_GROUP` |
| `waveica_v1` | Original WaveICA (Deng et al.) — wavelet+ICA correction using real batch labels directly, rather than the injection-order proxy WaveICA2.0 uses. May be preferable when batch labels are known and reliable. Defaults are the package's own, not tuned for this pipeline. Couples drift and batch correction internally, so **requires `drift=none`**. | `WAVEICA_V1_WF`, `WAVEICA_V1_K`, `WAVEICA_V1_T`, `WAVEICA_V1_T2`, `WAVEICA_V1_ALPHA` |
| `pmp_qcrsc` | QC-RSC (Quality Control-Robust Spline Correction) from the [pmp](https://bioconductor.org/packages/pmp/) package. Fits a smoothing spline through QC samples within each batch to correct signal drift and align batches in one step. Batches with fewer than 4 QC samples are left uncorrected (pmp cannot spline-fit them). Couples drift and batch correction internally, so **requires `drift=none`**. | — |
| `serrf` | SERRF (Systematic Error Removal using Random Forest). Per-feature random forest models trained on QC samples to correct systematic error. Adapted from [Fan et al., Analytical Chemistry 2019](https://doi.org/10.1021/acs.analchem.8b05592). Couples drift and batch correction internally, so **requires `drift=none`**. | `SERRF_NUM` |

### Examples

| Recipe | Equivalent to (pre-refactor naming, for reference) |
|---|---|
| `none:none:none` | Imputation only, no correction |
| `notame_spline:qc:ruv_s` | the original `notame` method |
| `loess:qc:combat` | `loess_combat` |
| `loess:hybrid:combat` | `loess_samples_combat` |
| `huber:qc:combat` | `huber_combat` |
| `huber:hybrid:combat` | `huber_samples_combat` |
| `loess:hybrid:sva` | `loess_samples_sva` |
| `huber:hybrid:sva` | `huber_samples_sva` |
| `auto:hybrid:combat` | `auto_combat` |
| `loess:qc:limma` | `loess_limma` |
| `loess:hybrid:limma` | `loess_samples_limma` |
| `loess:qc:feature_median` | `loess_feature_median` |
| `loess:qc:global_median` | `loess_global_median` |
| `none:none:cordbat` | `cordbat_only` |
| `loess:samples:cordbat` | `loess_cordbat` |
| `none:none:waveica` | `waveica` |
| `none:none:waveica_v1` | `waveica_v1` |
| `none:none:combat` | `combat_only` |
| `none:none:pmp_qcrsc` | `pmp_qcrsc` |
| `none:none:batchcorr` | `batchcorr` |
| `none:none:serrf` | `serrf` |
| `huber:qc:feature_median` | *(new — no pre-refactor equivalent)* |

Note: `loess:qc:*`/`huber:qc:*` add a batch-level minimum-QC gate (`DRIFT_MIN_QC_PER_BATCH`) that the old `loess_combat`/`huber_combat`/`loess_limma`/`loess_feature_median`/`loess_global_median` didn't have (they only checked per-feature QC counts) — batches below the threshold are now left uncorrected instead of partially corrected feature-by-feature. `pmp_qcrsc_scale`/`pmp_qcrsc_feature_scale` (the median-scaling fallback variants for low-QC batches) have been removed; use plain `pmp_qcrsc` instead.

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

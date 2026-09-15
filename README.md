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
| `QC_OUTLIER_MAD_K` | No | `5` | Multivariate outlier check for QC/ltQC samples, run separately per group and per batch, after `MIN_QC_SAMPLE_DETECTION`'s removal. A sample is flagged if its PCA-space distance (top 5 PCs) from its batch group's median score exceeds `median(distances) + K * mad(distances)`. Lower = more aggressive; higher = more conservative. Complements `MIN_QC_SAMPLE_DETECTION` — catches a sample with a normal detection rate but an anomalous intensity profile (contamination, carryover, a degrading/recovering column). Batches with fewer than `QC_OUTLIER_MIN_N` samples of a group are skipped. Set to `0` to disable |
| `QC_OUTLIER_MIN_N` | No | `3` | Minimum samples of a group (QC or ltQC) a batch must have for `QC_OUTLIER_MAD_K`'s check to run there. 3 is the practical floor — PCA correctly caps at 2 dimensions with 3 points and MAD of 3 distances is still well-defined (weak power, not meaningless), but 2 points give no third reference to judge "typical spread" against. Default is 3 rather than a higher conventional floor specifically because ltQC groups are often exactly 3 per batch by design — a higher default would silently skip ltQC's check in every batch for a setup like that |
| `MIN_BATCH_DETECTION` | No | `1` | Min number of detections a feature must have in every batch. Features absent from any entire batch are removed (set to `0` to disable) |
| `MIN_BATCH_DETECTION_FRAC` | No | `0` (disabled) | Min fraction (0-1) of each batch's samples a feature must be detected in. Unlike `MIN_BATCH_DETECTION` (an absolute count), this scales with batch size, so batches of very different sizes get a consistent relative bar rather than a fixed count that's stringent for a small batch and lax for a large one. Applied in addition to `MIN_BATCH_DETECTION`, not instead of it. Most relevant for `batch=waveica`/`waveica_v1`, where LoD/2-imputed placeholders participate directly in fitting the correction — a batch with disproportionately more missingness than others risks its placeholder pattern being mistaken for real signal |
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
| `AUTO_MIN_LTQC_VALIDATE` | No | `3` | Min ltQC samples a batch needs to select (and receive) its own `drift=auto` samples-based candidate at all. Under `basis=samples`, `hybrid`, and `qc` alike, a batch below this threshold has no held-out evidence to select from and is left uncorrected — each batch's selection depends only on its own ltQC/QC count, never on other batches |
| `AUTO_MIN_CV_OBS` | No | `4` | Min finite training observations (QC, or samples for the samples-based pool) a feature needs before `drift=auto` attempts to fit any candidate for it |
| `HUBER_QC_K` | No | `1.345` | Huber regression tuning constant (`MASS::rlm`, `psi.huber`) for `drift=huber`'s QC-based fit (`basis=qc`, or `basis=hybrid`'s QC branch). Fixed, not auto-searched — see `HUBER_QC_CV_KS` for per-feature CV selection, or `AUTO_HUBER_KS`/`drift=auto` for a per-batch CV-chosen `k` instead |
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
| `WAVEICA_V1_WF` | No | `haar` | Wavelet family for `batch=waveica_v1` (the original WaveICA — separate setting from `WAVEICA_WF`, different package). Not searched |
| `WAVEICA_V1_K` | No | `20` | Comma-separated number(s) of components `batch=waveica_v1`'s ICA step decomposes into. Each entry is a number or `auto` (`2 x n_batches`, same convention as `WAVEICA_K` — not something the original package defines itself). Single value fixed, multiple values searched |
| `WAVEICA_V1_T` | No | `0.05` | Comma-separated p-value threshold(s) (0-1) for considering a component associated with batch in `batch=waveica_v1` — tested against real batch labels directly, unlike `WAVEICA_CUTOFF`'s GAM-R²-against-injection-order proxy; a component is removed if its p-value falls below this. Higher = more aggressive (looser bar, more components removed); lower = more conservative — opposite direction from `WAVEICA_CUTOFF`, which thresholds an R² a component must *exceed*, not a p-value it must fall under. Single value fixed, multiple values searched |
| `WAVEICA_V1_T2` | No | `0.05` | Threshold (0-1) for considering a component associated with a biological comparison group in `batch=waveica_v1`. Currently inert — this pipeline has no biological-group column to supply `waveica_v1`'s optional `group` argument. Not searched, since it has no effect here |
| `WAVEICA_V1_ALPHA` | No | `0` | Comma-separated trade-off value(s) (0-1) between sample-wise and variable-wise independence in `batch=waveica_v1`'s ICA step. The same *kind* of parameter as `WAVEICA_ALPHA` (both are ICA spatial/temporal trade-offs) — kept as a separate setting since they're independent packages with separately-tuned defaults, not because the concept differs. Single value fixed, multiple values searched |
| `WAVEICA_V1_EVAL_GROUP` | No | `ltQC` | Which group (`ltQC` or `QC`) the `WAVEICA_V1_ALPHA`/`WAVEICA_V1_T`/`WAVEICA_V1_K` search evaluates candidates against (D-ratio vs. Sample). Same mechanism as `WAVEICA_EVAL_GROUP` |

**D-ratio alone isn't trusted for `waveica`/`waveica_v1` search selection.** Real runs surfaced two distinct ways a candidate can "win" on D-ratio without actually being a good correction: shrinking Sample variance faster than eval_group variance (destroying biological signal, not removing noise), or inflating scatter *within* the eval_group replicates themselves (QC/ltQC samples should stay tight regardless of correction). Both showed up as a visibly worse QC/ltQC clustering in PCA despite a "better" D-ratio. To catch this, both search functions also compute `dist_ratio` (PCA-space distance ratio, eval_group vs. Sample) for the *uncorrected* data as a baseline, and exclude any candidate whose `dist_ratio` is worse than that baseline before picking the D-ratio winner — printed per-candidate as `[excluded: dist_ratio worse than uncorrected]`. If every candidate fails the guard (or the baseline itself can't be computed), it falls back to the first grid entry with a clear message, same as the existing "no computable D-ratio" fallback.

| `COMBAT_MEAN_ONLY` | No | `auto` | Whether ComBat (`batch=combat`, with any drift method) adjusts only each feature's per-batch mean (`TRUE`) or also forces every batch's variance to match a common value (`FALSE`, ComBat's own default). Forcing variance equal across batches is the usual cause of PCA looking artificially "flattened" after correction when batches genuinely differ in spread. `auto` tries both and keeps whichever gives the better ltQC/Sample D-ratio; set `TRUE`/`FALSE` to force a specific behaviour |
| `COMBAT_PAR_PRIOR` | No | `auto` | Whether ComBat's empirical Bayes prior is estimated parametrically (`TRUE`, assumes a Normal/Inverse-Gamma shape — faster) or non-parametrically (`FALSE` — slower, more robust to non-Gaussian batch effects). `auto` tries both and keeps whichever gives the better ltQC/Sample D-ratio, same mechanism as `COMBAT_MEAN_ONLY` |
| `N_CORES` | No | all - 1 | Number of CPU cores for parallelisation |
| `RUV_K` | No | `3` | Unwanted variation factors for `batch=ruv_s` (notame's RUV-S) |

**Recommended first-pass search grid for `batch=waveica`.** The single-value defaults above (`0.05`/`0.10`/`auto`) are deliberately fixed — no search unless you ask for one. For a new dataset, we recommend running one broad characterization search before narrowing down:

```
-e WAVEICA_ALPHA="0,0.25,0.5,0.75,1"
-e WAVEICA_CUTOFF="0.05,0.10,0.20,0.30,0.40,0.50"
-e WAVEICA_K="auto,10,20,40"
```

120 combinations, parallelized across `N_CORES`. Rationale: `WAVEICA_ALPHA` spans its full 0–1 range at even spacing to catch either a monotonic trend toward one boundary or an interior optimum; `WAVEICA_CUTOFF` is denser at the low-to-moderate end where component selection is most sensitive, coarser above 0.30 since very high values tend to converge toward a no-op; `WAVEICA_K` spans roughly 2×–10×+ `auto`'s own value (`2 x n_batches`), the widest relative range of the three since it showed the least sign of plateauing in practice. Treat this as a one-off scan to characterize a new dataset, not something to rerun by default — once you see the shape of the response surface, a narrower follow-up grid (e.g. fixing `WAVEICA_ALPHA` near its apparent optimum and refining `WAVEICA_CUTOFF`/`WAVEICA_K` resolution around wherever the best region showed up) is usually more useful than repeating the full scan.

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
| `auto` | Auto-selected via cross-validation from several LOESS spans, several Huber `k` values, and a flat/no-op baseline. Each batch selects its own winner independently, from its own evidence only (leave-one-out CV on that batch's own QC for the QC-based selection; held-out ltQC/Sample D-ratio on that batch's own samples for the samples-based selection) — different batches are expected to end up with different methods, nothing is pooled or shared across batches. | `AUTO_LOESS_SPANS`, `AUTO_HUBER_KS`, `AUTO_SAMPLE_LOESS_SPANS`, `AUTO_SAMPLE_HUBER_KS`, `AUTO_MIN_QC_PER_BATCH`, `AUTO_MIN_LTQC_VALIDATE`, `AUTO_MIN_CV_OBS` |
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
| `waveica` | WaveICA 2.0 — wavelet-based correction for both drift and batch effects, QC-independent ([Deng et al. 2021](https://link.springer.com/article/10.1007/s11306-021-01839-7)). Uses injection order as a proxy for batch structure rather than batch labels directly. Couples drift and batch correction internally, so **requires `drift=none`**. `WAVEICA_ALPHA`/`WAVEICA_CUTOFF`/`WAVEICA_K` each accept a comma-separated list; more than one candidate overall triggers a search over the full cross-product, evaluated against `WAVEICA_EVAL_GROUP`/Sample D-ratio (primary — selects the winner among candidates that pass a `dist_ratio` guard, see below), with PERMANOVA R²(Batch) printed alongside every candidate as an independent, informational cross-check. Candidates run in parallel via this pipeline's existing `foreach`/`N_CORES` setup, with `mc.cores` pinned to 1 inside each worker to avoid nested-parallelism oversubscription against WaveICA2.0's own internal `parallel::mclapply()` call. `WAVEICA_WF` is not searched. | `WAVEICA_ALPHA`, `WAVEICA_CUTOFF`, `WAVEICA_K`, `WAVEICA_WF`, `WAVEICA_EVAL_GROUP` |
| `waveica_v1` | Original WaveICA (Deng et al.) — wavelet+ICA correction using real batch labels directly, rather than the injection-order proxy WaveICA2.0 uses. May be preferable when batch labels are known and reliable. Defaults are the package's own, not tuned for this pipeline. **Unlike every other method on this list, `waveica_v1` has no drift-correction mechanism of its own at all** — it never sees injection order, and its component-removal test (a one-way ANOVA on per-batch means) is structurally blind to a trend within a batch. So it does **not** require `drift=none`: pair it with a real drift step for data with within-batch drift (e.g. `loess:hybrid:waveica_v1`), or use `none:none:waveica_v1` for batch-only correction on data that's already drift-free. `WAVEICA_V1_ALPHA`/`WAVEICA_V1_T`/`WAVEICA_V1_K` each accept a comma-separated list, with the same fixed-if-single / searched-if-multiple convention and D-ratio-as-primary-selection design as `waveica` above — see `WAVEICA_V1_EVAL_GROUP`. `WAVEICA_V1_T2` and `WAVEICA_V1_WF` are not searched. Unlike WaveICA2.0, WaveICA (v1) has no internal parallelism of its own, so candidates don't need an `mc.cores` pin. | `WAVEICA_V1_WF`, `WAVEICA_V1_K`, `WAVEICA_V1_T`, `WAVEICA_V1_T2`, `WAVEICA_V1_ALPHA`, `WAVEICA_V1_EVAL_GROUP` |
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
4. QC/ltQC multivariate outlier check — removes individual QC or ltQC samples that are PCA-space distance outliers relative to their own batch's other samples of the same group (`QC_OUTLIER_MAD_K`; disabled by setting to `0`). Runs after step 3, on the samples that survive it. Catches a different failure mode than step 3 — a sample that detects fine but has an anomalous intensity profile (contamination, carryover, a degrading/recovering column)
5. QC detection — removes features not detected in ≥ `QC_DETECTION_LIMIT` of QC samples
6. Sample detection — removes features not detected in ≥ `SAMPLE_DETECTION_LIMIT` of biological samples
7. Zero variance — removes features with no variation across samples
8. Per-batch detection (absolute) — removes features with fewer than `MIN_BATCH_DETECTION` real observations in any batch
9. Per-batch detection (fraction) — removes features detected in less than `MIN_BATCH_DETECTION_FRAC` of any batch's samples (disabled by default; see [Correction methods](#correction-methods) parameter table)
10. QC-RSD filter — removes features with QC robust RSD (MAD/median) > `QC_RSD_FILTER` in ≥ 50% of batches (disabled by default)

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

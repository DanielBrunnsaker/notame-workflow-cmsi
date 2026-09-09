# ─────────────────────────────────────────────────────────────────────────────
# Preprocessing pipeline for untargeted LC-MS metabolomics. Converts a
# MSDIAL alignment export or an XCMS feature table + sample sheet to notame
# format, applies pre-correction feature filters, and runs one or more
# drift/batch correction methods in parallel. QC metrics are computed per
# method for comparison.
#
# Daniel Brunnsåker, 2026-03-20
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(notame)
  library(notameViz)
  library(notameStats)
  library(openxlsx)
  library(doParallel)
  library(yaml)
  library(jsonlite)
})

source("R/config.R")
source("R/preflight.R")
source("R/notame_format.R")
source("R/msdial_to_notame.R")
source("R/xcms_to_notame.R")
source("R/qc_metrics.R")
source("R/drift_correction.R")
source("R/correction_methods.R")
source("R/serrf.R")
source("R/run_log.R")
source("R/report.R")

# Optional YAML config file — supplies default parameter values. Env vars
# still take precedence over the config file (see get_env() below), so
# Docker one-off overrides keep working unchanged.
config_file <- Sys.getenv("CONFIG_FILE", unset = "")
config <- if (config_file != "") load_config(config_file) else list()

get_env <- function(var, default) {
  val <- Sys.getenv(var, unset = NA)
  if (!is.na(val) && val != "") return(val)
  if (!is.null(config[[var]]) && nchar(config[[var]]) > 0) return(config[[var]])
  default
}

n_cores_env <- get_env("N_CORES", "")
registerDoParallel(cores = if (n_cores_env == "") parallel::detectCores() - 1 else as.integer(n_cores_env))

# ─────────────────────────────────────────────────────────────────────────────
# SETTINGS
# All parameters can be overridden via environment variables, or supplied in
# a YAML config file (CONFIG_FILE) — env vars take precedence over the
# config file when both are set.
# ─────────────────────────────────────────────────────────────────────────────

if ("--help" %in% commandArgs(trailingOnly = TRUE)) {
  cat("
Usage: Rscript notame-workflow.r [--help]

Environment variables (required variables are marked; all others are optional with defaults shown).
Any of these may instead be supplied via a YAML file pointed to by CONFIG_FILE — env vars win
if both are set for the same parameter.

  CONFIG_FILE            Path to an optional YAML config file supplying default values for any
                        parameter below, plus an optional 'sample_type_rules' list overriding
                        the QC/ltQC/Blank/Wash/etc. filename-keyword classification. See
                        config.example.yaml.
                        Default: (none)

  IN_XLSX               Path to the MSDIAL alignment export (.xlsx). Required, unless
                        IN_FEATURE_TABLE + IN_SAMPLE_SHEET are set instead.

  IN_FEATURE_TABLE      Path to an XCMS-based feature table (.csv) — alternative to IN_XLSX.
                        Must be set together with IN_SAMPLE_SHEET.

  IN_SAMPLE_SHEET       Path to the XCMS pipeline's sample sheet (.xlsx) — alternative to
                        IN_XLSX. Must be set together with IN_FEATURE_TABLE.

  PROJECT_FOLDER        Root output directory (intermediates/ and output/ written here). Required.

  FORCE_RECONVERT       Force re-running the conversion step even if a cached
                        notame-formatted file from a previous run is found.
                        Default: FALSE
                        Values:  TRUE | FALSE

  COLUMN                Chromatographic column type. Required. Used together with POLARITY
                        to namespace output folders (e.g. RP_POS, HILIC_NEG).
                        Examples: RP | HILIC

  POLARITY              Ionisation polarity. Required.
                        Values:  POS | NEG

  CORRECTION_METHODS    Comma-separated list of correction methods to run.
                        Each method is saved to its own output subfolder.
                        Default: none,notame
                        Values:  none | notame | pmp_qcrsc | pmp_qcrsc_scale | pmp_qcrsc_feature_scale | serrf |
                                 batchcorr | combat_only | loess_combat | loess_samples_combat |
                                 loess_limma | loess_samples_limma | loess_feature_median | loess_global_median |
                                 cordbat_only | loess_cordbat | waveica

  QC_DETECTION_LIMIT    Min fraction of QC samples a feature must be detected in
                        Default: 0.60

  SAMPLE_DETECTION_LIMIT  Min fraction of biological samples a feature must be detected in
                        Default: 0.20

  BLANK_RATIO           Remove features where mean(Sample) <= BLANK_RATIO * mean(SolvBlank).
                        Set to 'none' to disable.
                        Default: none

  LOW_INT_FILTER_FRAC   Data-driven low-intensity filter. Removes features whose p80
                        intensity is below this fraction of the mean p80 across all features.
                        Ignored if LOW_INT_FILTER is set.
                        Default: 0.10

  LOW_INT_FILTER        Absolute low-intensity threshold (overrides LOW_INT_FILTER_FRAC).
                        Default: (disabled)

  LOW_INT_PERCENTILE    Percentile used for the low-intensity filter (0-1).
                        Default: 0.8

  MIN_QC_SAMPLE_DETECTION  Minimum fraction of features that must be detected in a QC or ltQC
                        sample for it to be used as a reference. Samples below this threshold
                        are removed before processing (empty injections, failed runs).
                        Default: 0.50

  MIN_BATCH_DETECTION   Minimum number of detections a feature must have in every batch.
                        Features absent from any entire batch are removed — they have no
                        real measurements there and would be all-imputed placeholders.
                        Set to 0 to disable.
                        Default: 1

  QC_RSD_FILTER         Pre-correction QC-RSD threshold (robust: MAD/median). Features are
                        removed if they fail the threshold in >= 50% of evaluable batches.
                        Set to 'none' to disable.
                        Default: none (disabled)

  RSD_THRESHOLD         RSD threshold used for post-correction output filtering.
                        Controls both the global (_rsdXX) and per-batch (_batchrsdXX) outputs.
                        Default: 0.30

  RUV_K                 Number of unwanted variation factors for RUV (notame method only).
                        Default: 3

  SERRF_NUM             Number of correlated features used as RF predictors per feature model
                        in SERRF. Lower = less overfitting risk; recommended ≤ floor(QC_n/2).
                        Automatically capped at floor(min_QC_per_batch / 2) at runtime.
                        Default: 5

  LOESS_SPAN            LOESS smoothing span for QC-based drift correction (loess_combat, loess_limma,
                        loess_feature_median, loess_global_median). Higher = smoother, more
                        conservative correction.
                        Default: 0.75

  LOESS_SAMPLE_SPAN     LOESS smoothing span for QC-free drift correction (loess_samples_combat,
                        loess_samples_limma, loess_cordbat), fit on biological samples instead of
                        QC. Wider than LOESS_SPAN by default since each point is a unique
                        biological measurement, not a technical replicate — a tighter span risks
                        fitting individual-sample noise as drift.
                        Default: 0.9

  LOESS_SAMPLE_MIN_OBS  Minimum finite sample observations per feature required to attempt
                        QC-free drift correction (loess_samples_combat, loess_samples_limma,
                        loess_cordbat). Deliberately higher than the QC-based fit's threshold of
                        4, since sample points are far noisier. Features below this are left
                        uncorrected for that batch.
                        Default: 10

  LOESS_MIN_QC_PER_BATCH  Minimum QC samples a batch must have to use QC-based drift correction
                        in loess_samples_combat / loess_samples_limma. These methods choose per
                        batch: batches meeting this threshold are drift-corrected with QC-based
                        LOESS (LOESS_SPAN). Batches below it try the QC-free fit on biological
                        samples (LOESS_SAMPLE_SPAN, LOESS_SAMPLE_MIN_OBS) as a trial, kept only if
                        it improves the ltQC/Sample D-ratio — see LOESS_MIN_LTQC_VALIDATE. QC-based
                        fitting is preferred whenever there's enough QC to support it, since it
                        doesn't risk removing real biological signal along with drift.
                        Default: 4

  LOESS_MIN_LTQC_VALIDATE  Minimum ltQC samples a batch must have to validate the QC-free trial
                        correction in loess_samples_combat / loess_samples_limma (used only for
                        batches below LOESS_MIN_QC_PER_BATCH). The trial correction is compared
                        against the batch's uncorrected values by ltQC/Sample D-ratio
                        (MAD(ltQC)/MAD(Sample) — lower is better); it is kept if the D-ratio
                        improves, discarded otherwise. D-ratio (not raw ltQC RSD) is used because
                        any real drift correction shrinks measured sample variance somewhat, so a
                        test that just required sample variance not to shrink would reject working
                        corrections too — D-ratio instead only credits a disproportionate
                        improvement in ltQC relative to Sample. Batches with fewer ltQC than this
                        have no way to validate the trial and are left uncorrected.
                        Default: 3

  LOESS_VALIDATE_SAMPLES_CORRECTION  Set to FALSE to skip the ltQC validation above entirely for
                        loess_samples_combat / loess_samples_limma: any batch below
                        LOESS_MIN_QC_PER_BATCH then always gets the QC-free samples-based
                        correction, regardless of ltQC availability or what it shows. This
                        reintroduces the risk the validation step exists to catch (the trial can
                        look fine on ltQC while still compressing real biological signal) — use
                        deliberately, not as a default.
                        Default: TRUE

  CORDBAT_REF_BATCH     Reference batch ID for CordBat (cordbat_only, loess_cordbat).
                        All other batches are corrected onto this batch.
                        Leave unset to auto-select the batch with the lowest median feature RSD.
                        Default: (auto)

  WAVEICA_ALPHA         Significance threshold WaveICA2.0 (waveica) uses to decide whether an
                        independent component is injection-order-related and should be removed.
                        Lower = stricter (fewer components flagged, less aggressive correction);
                        higher = more components flagged and removed. Try lowering this first if
                        correction looks too aggressive (flattens real sample-to-sample variation).
                        Default: 0.05

  WAVEICA_CUTOFF        Threshold (0-1) for how much of a wavelet-decomposed level's variance must
                        be associated with injection order before that level is considered
                        technical and passed to ICA for cleanup. Lower = more levels get corrected
                        (more aggressive); higher = fewer, more conservative.
                        Default: 0.10

  WAVEICA_K             Number of independent components WaveICA2.0 decomposes the data into.
                        Default here is 2 x (number of batches), which is quite small for typical
                        feature counts — a coarse decomposition gives ICA less room to isolate
                        narrow technical components from broad biological ones, which may be why
                        correction looks like it removes more than just drift/batch noise. Try
                        raising this (e.g. 10-20) if correction looks too aggressive.
                        Default: (auto = 2 x n_batches)

  WAVEICA_WF            Wavelet family used for the decomposition step (passed to WaveICA2.0's
                        wf argument, e.g. haar or a Daubechies family recognised by the
                        underlying wavelet package).
                        Default: haar

                        Note: these four are exposed as-is from the WaveICA2.0 package
                        (github.com/dengkuistat/WaveICA_2.0); the direction-of-effect guidance
                        above is based on the published method rather than inspection of this
                        specific package version's source, so confirm empirically on your data.

  NORMALIZATION         Post-correction normalisation method. Uses pooled QC samples as
                        reference when available, otherwise median of biological samples.
                        Default: none
                        Values:  none | pqn

  SAVE_PRE_CORRECTION_PLOTS  Save QC plots before correction (slow on large datasets).
                        Default: TRUE
                        Values:  TRUE | FALSE

  N_CORES               Number of CPU cores for parallelisation.
                        Default: all available cores minus one

")
  quit(status = 0)
}

project_folder <- get_env("PROJECT_FOLDER", "")
if (project_folder == "") stop("PROJECT_FOLDER is required. Set it via env var or CONFIG_FILE.")

# Input: an MSDIAL alignment export (IN_XLSX), or an XCMS-based feature
# table + sample sheet (IN_FEATURE_TABLE + IN_SAMPLE_SHEET) — mutually
# exclusive, auto-detected from which are set.
in_xlsx          <- get_env("IN_XLSX", "")
in_feature_table <- get_env("IN_FEATURE_TABLE", "")
in_sample_sheet  <- get_env("IN_SAMPLE_SHEET", "")

has_msdial <- in_xlsx != ""
has_xcms   <- in_feature_table != "" || in_sample_sheet != ""

if (has_msdial && has_xcms) {
  stop("Set either IN_XLSX (MSDIAL) or IN_FEATURE_TABLE + IN_SAMPLE_SHEET (XCMS), not both.")
} else if (has_xcms) {
  if (in_feature_table == "" || in_sample_sheet == "")
    stop("XCMS input requires both IN_FEATURE_TABLE and IN_SAMPLE_SHEET to be set.")
  input_mode <- "xcms"
} else if (has_msdial) {
  input_mode <- "msdial"
} else {
  stop("An input is required: set IN_XLSX (MSDIAL) or IN_FEATURE_TABLE + IN_SAMPLE_SHEET (XCMS), via env var or CONFIG_FILE.")
}

polarity <- get_env("POLARITY", "")
if (polarity == "") stop("POLARITY is required. Set it to 'POS' or 'NEG' (env var or CONFIG_FILE).")

column <- get_env("COLUMN", "")
if (column == "") stop("COLUMN is required. Set it to the chromatographic column type, e.g. 'RP' or 'HILIC' (env var or CONFIG_FILE).")

mode_label <- paste0(column, "_", polarity)  # e.g. RP_POS — used to namespace all output folders

out_xlsx   <- file.path(project_folder, "intermediates", mode_label, "notame_rev.xlsx")
interdir   <- file.path(project_folder, "intermediates", mode_label)
output_dir <- file.path(project_folder, "output",        mode_label)

QC_DETECTION_LIMIT     <- as.numeric(get_env("QC_DETECTION_LIMIT",     "0.60"))
SAMPLE_DETECTION_LIMIT <- as.numeric(get_env("SAMPLE_DETECTION_LIMIT", "0.20"))

blank_ratio_env <- get_env("BLANK_RATIO", "none")
BLANK_RATIO <- if (blank_ratio_env %in% c("none", "skip", "")) NA_real_ else
               as.numeric(blank_ratio_env)

LOW_INT_FILTER      <- suppressWarnings(as.numeric(get_env("LOW_INT_FILTER", "")))
LOW_INT_FILTER_FRAC <- as.numeric(get_env("LOW_INT_FILTER_FRAC", "0.10"))
LOW_INT_PERCENTILE  <- as.numeric(get_env("LOW_INT_PERCENTILE",  "0.8"))
qc_rsd_env    <- get_env("QC_RSD_FILTER", "none")
QC_RSD_FILTER           <- if (qc_rsd_env %in% c("none", "")) NA_real_ else as.numeric(qc_rsd_env)
MIN_QC_SAMPLE_DETECTION <- as.numeric(get_env("MIN_QC_SAMPLE_DETECTION", "0.50"))
MIN_BATCH_DETECTION     <- as.integer(get_env("MIN_BATCH_DETECTION", "1"))
RSD_THRESHOLD <- as.numeric(get_env("RSD_THRESHOLD", "0.30"))

# Correction methods:
#   "none"                — imputation only (no correction; baseline)
#   "notame"              — per-batch cubic spline drift correction + RUV batch correction
#   "pmp_qcrsc"           — QC-RSC spline drift correction (pmp package)
#   "pmp_qcrsc_scale"         — QC-RSC drift correction (pmp) + global median scaling for any
#                               batch with <4 QCs; leaves pmp-corrected batches untouched
#   "pmp_qcrsc_feature_scale" — as pmp_qcrsc_scale but uses per-feature median scaling for
#                               no-QC batches, consistent with pmp's own feature-wise alignment
#   "serrf"               — SERRF random forest correction (Fan et al. 2019)
#   "batchcorr"           — cluster-based spline drift + between-batch normalisation (Brunius et al.)
#   "combat_only"         — ComBat batch correction only (no drift correction)
#   "loess_combat"        — per-batch LOESS drift correction (QC-based) + ComBat batch correction
#   "loess_samples_combat" — per-batch LOESS drift correction (QC-based if enough QC, else a
#                            QC-free trial on samples kept only if it improves the ltQC/Sample
#                            D-ratio, else uncorrected — see LOESS_MIN_QC_PER_BATCH,
#                            LOESS_MIN_LTQC_VALIDATE) + ComBat
#   "loess_limma"         — per-batch LOESS drift correction (QC-based) + limma removeBatchEffect
#   "loess_samples_limma" — same per-batch QC-based/QC-free-trial/uncorrected choice as
#                            loess_samples_combat + limma removeBatchEffect
#   "cordbat_only"        — CordBat batch correction only (GGM-based, no drift correction)
#   "loess_cordbat"       — per-batch QC-free LOESS drift correction (fit on samples, always —
#                            CordBat's own between-batch step is also QC-free, fit on samples)
#                            + CordBat batch correction
#   "waveica"             — WaveICA 2.0 wavelet-based correction
CORRECTION_METHODS <- strsplit(get_env("CORRECTION_METHODS", "none,notame"), ",")[[1]]

RUV_K      <- as.integer(get_env("RUV_K",       "3"))
SERRF_NUM  <- as.integer(get_env("SERRF_NUM",   "5"))
LOESS_SPAN                <- as.numeric(get_env("LOESS_SPAN", "0.75"))
LOESS_SAMPLE_SPAN         <- as.numeric(get_env("LOESS_SAMPLE_SPAN", "0.9"))
LOESS_SAMPLE_MIN_OBS      <- as.integer(get_env("LOESS_SAMPLE_MIN_OBS", "10"))
LOESS_MIN_QC_PER_BATCH    <- as.integer(get_env("LOESS_MIN_QC_PER_BATCH", "4"))
LOESS_MIN_LTQC_VALIDATE   <- as.integer(get_env("LOESS_MIN_LTQC_VALIDATE", "3"))
LOESS_VALIDATE_SAMPLES_CORRECTION <- as.logical(get_env("LOESS_VALIDATE_SAMPLES_CORRECTION", "TRUE"))
cordbat_ref_env   <- get_env("CORDBAT_REF_BATCH", "")
CORDBAT_REF_BATCH <- if (cordbat_ref_env == "") NULL else cordbat_ref_env
WAVEICA_ALPHA   <- as.numeric(get_env("WAVEICA_ALPHA",  "0.05"))
WAVEICA_CUTOFF  <- as.numeric(get_env("WAVEICA_CUTOFF", "0.10"))
waveica_k_env   <- get_env("WAVEICA_K", "")
WAVEICA_K       <- if (waveica_k_env == "") NULL else as.integer(waveica_k_env)
WAVEICA_WF      <- get_env("WAVEICA_WF", "haar")
NORMALIZATION              <- get_env("NORMALIZATION",              "none")
SAVE_PRE_CORRECTION_PLOTS  <- as.logical(get_env("SAVE_PRE_CORRECTION_PLOTS", "TRUE"))
FORCE_RECONVERT            <- as.logical(get_env("FORCE_RECONVERT", "FALSE"))

# ─────────────────────────────────────────────────────────────────────────────
# PREFLIGHT CHECKS
# ─────────────────────────────────────────────────────────────────────────────

run_preflight_checks(
  input_mode = input_mode, in_xlsx = in_xlsx,
  in_feature_table = in_feature_table, in_sample_sheet = in_sample_sheet,
  project_folder = project_folder, column = column, polarity = polarity,
  correction_methods = CORRECTION_METHODS, normalization = NORMALIZATION,
  qc_detection_limit = QC_DETECTION_LIMIT, sample_detection_limit = SAMPLE_DETECTION_LIMIT,
  low_int_filter_frac = LOW_INT_FILTER_FRAC, low_int_percentile = LOW_INT_PERCENTILE,
  min_qc_sample_detection = MIN_QC_SAMPLE_DETECTION, min_batch_detection = MIN_BATCH_DETECTION,
  rsd_threshold = RSD_THRESHOLD, ruv_k = RUV_K, serrf_num = SERRF_NUM, loess_span = LOESS_SPAN,
  loess_sample_span = LOESS_SAMPLE_SPAN, loess_sample_min_obs = LOESS_SAMPLE_MIN_OBS,
  loess_min_qc_per_batch = LOESS_MIN_QC_PER_BATCH,
  loess_min_ltqc_validate = LOESS_MIN_LTQC_VALIDATE,
  waveica_alpha = WAVEICA_ALPHA, waveica_cutoff = WAVEICA_CUTOFF,
  waveica_k = if (is.null(WAVEICA_K)) NA_integer_ else WAVEICA_K,
  blank_ratio = BLANK_RATIO, low_int_filter = LOW_INT_FILTER, qc_rsd_filter = QC_RSD_FILTER,
  save_pre_correction_plots = SAVE_PRE_CORRECTION_PLOTS,
  config_file = config_file, raw_sample_type_rules = config$sample_type_rules
)

# Resolved after preflight has validated it (pattern/type present, pattern a
# valid regex) so all config problems are still reported together up front.
sample_type_rules <- resolve_sample_type_rules(config)

# ─────────────────────────────────────────────────────────────────────────────
# 1) CONVERT & IMPORT
# ─────────────────────────────────────────────────────────────────────────────

dir.create(interdir, showWarnings = FALSE, recursive = TRUE)

# Cache the conversion step: re-parsing the source export(s) is one of the
# slower steps and only needs to be redone when the input(s) change. The
# annotations table (returned in memory, not part of out_xlsx) is cached
# alongside it so a cache hit can skip the converter call entirely.
input_files     <- if (input_mode == "xcms") c(in_feature_table, in_sample_sheet) else in_xlsx
annotations_rds <- file.path(interdir, "conversion_annotations.rds")
newest_input    <- max(file.mtime(input_files))
cache_valid <- !FORCE_RECONVERT &&
  file.exists(out_xlsx) && file.exists(annotations_rds) &&
  file.mtime(out_xlsx)        >= newest_input &&
  file.mtime(annotations_rds) >= newest_input &&
  (config_file == "" ||
     (file.mtime(out_xlsx) >= file.mtime(config_file) &&
      file.mtime(annotations_rds) >= file.mtime(config_file)))  # invalidate on sample_type_rules changes too

if (cache_valid) {
  message("==> Using cached ", toupper(input_mode), " conversion: ", out_xlsx,
          " (set FORCE_RECONVERT=TRUE to re-run)")
  convert_result <- readRDS(annotations_rds)
} else if (input_mode == "xcms") {
  convert_result <- xcms_to_notame(in_feature_table, in_sample_sheet, out_xlsx, column, polarity)
  saveRDS(convert_result, annotations_rds)
} else {
  convert_result <- msdial_to_notame(in_xlsx, out_xlsx, sample_type_rules = sample_type_rules)
  saveRDS(convert_result, annotations_rds)
}

mode_name           <- convert_result$mode
feature_annotations <- convert_result$annotations

message("==> Importing")
data <- import_from_excel(file = out_xlsx, sheet = 1, name = mode_name)
names(assays(data)) <- "abundances"
data <- fix_object(data, assay.type = "abundances")

cat("Imported:", nrow(data), "features,", ncol(data), "samples\n")
print(table(colData(data)$QC))

# ─────────────────────────────────────────────────────────────────────────────
# 2) PRE-CORRECTION FEATURE FILTERING
# ─────────────────────────────────────────────────────────────────────────────

message("==> Pre-correction feature filtering")

data    <- mark_nas(data, value = 0)
n_before <- nrow(data)

# Blank filter
if (!is.na(BLANK_RATIO)) {
  blank_idx <- which(colData(data)$QC == "Blank")
  if (length(blank_idx) > 0) {
    sample_idx   <- which(colData(data)$QC == "Sample")
    blank_means  <- rowMeans(assay(data)[, blank_idx,  drop = FALSE], na.rm = TRUE)
    sample_means <- rowMeans(assay(data)[, sample_idx, drop = FALSE], na.rm = TRUE)
    keep_blank   <- is.na(blank_means) | blank_means == 0 | sample_means > BLANK_RATIO * blank_means
    keep_blank[is.na(keep_blank)] <- FALSE
    data <- data[keep_blank, ]
  } else {
    cat("No blank samples found — skipping blank filter\n")
  }
}
n_after_blank <- nrow(data)

# Remove non-analytical sample types; retain ltQC for downstream evaluation
data <- data[, !colData(data)$QC %in% c("Blank", "Wash", "Cond", "MSe", "MS2", "SST", "MatrixBlank")]

# Low-intensity filter: remove features below a fraction of the mean pN intensity
{
  mat_int      <- assay(data); mat_int[is.na(mat_int)] <- 0
  int_quantile <- apply(mat_int, 1, quantile, probs = LOW_INT_PERCENTILE)
  mean_p       <- mean(int_quantile[int_quantile > 0])

  low_int_cutoff <- if (!is.na(LOW_INT_FILTER)) {
    LOW_INT_FILTER
  } else if (!is.na(LOW_INT_FILTER_FRAC)) {
    LOW_INT_FILTER_FRAC * mean_p
  } else {
    NA_real_
  }

  if (!is.na(low_int_cutoff)) {
    data <- data[int_quantile >= low_int_cutoff, ]
    cat(sprintf("Low-intensity filter (p%.0f < %.4g): removed %d features\n",
                LOW_INT_PERCENTILE * 100, low_int_cutoff, nrow(data) - n_after_blank))
  }
}
n_after_lowint <- nrow(data)

# Remove QC/ltQC samples with insufficient feature detection (empty injections,
# failed runs). Runs after feature filters so detection rate is assessed on
# meaningful features only. ltQC is checked here too, using the same threshold
# as QC — it's the held-out group the loess_samples_combat/limma validation
# trial and the ltqc_permanova_* metrics rely on, and with typically only a
# few ltQC per batch, one bad injection would badly distort both.
for (qc_group in c("QC", "ltQC")) {
  qc_cols <- which(colData(data)$QC == qc_group)
  if (length(qc_cols) == 0) next
  mat_qc    <- assay(data)[, qc_cols, drop = FALSE]
  detect_qc <- colMeans(mat_qc > 0 & !is.na(mat_qc))
  bad_qc    <- qc_cols[detect_qc < MIN_QC_SAMPLE_DETECTION]
  if (length(bad_qc) > 0) {
    bad_names <- colData(data)$Sample_ID[bad_qc]
    bad_batch <- colData(data)$Batch[bad_qc]
    message("==> Removing ", length(bad_qc), " ", qc_group, " sample(s) with detection rate < ",
            round(MIN_QC_SAMPLE_DETECTION * 100), "%:")
    for (k in seq_along(bad_names))
      message("    ", bad_names[k], " (batch: ", bad_batch[k], ", detection: ",
              round(detect_qc[detect_qc < MIN_QC_SAMPLE_DETECTION][k] * 100, 1), "%)")
    data <- data[, -bad_qc]
  }
}

# QC detection filter
data <- flag_detection(data, qc_limit = QC_DETECTION_LIMIT)
data <- drop_flagged(data)
n_after_qc <- nrow(data)

# Biological sample detection filter
sample_idx  <- which(colData(data)$QC == "Sample")
detect_rate <- rowMeans(!is.na(assay(data)[, sample_idx, drop = FALSE]))
data        <- data[detect_rate >= SAMPLE_DETECTION_LIMIT, ]
n_after_sample <- nrow(data)

# Zero-variance filter
zero_var <- apply(assay(data), 1, function(x) {
  v <- var(x, na.rm = TRUE)
  !is.na(v) && v < .Machine$double.eps
})
data <- data[!zero_var, ]
n_after_zerovar <- nrow(data)

# Per-batch detection filter: feature must have >= MIN_BATCH_DETECTION observations
# in every batch. Features absent from a whole batch have no real measurements
# there — all values would be LoD/2 placeholders.
if (MIN_BATCH_DETECTION > 0) {
  cd_batch  <- as.data.frame(colData(data))
  batches   <- unique(cd_batch$Batch)
  batch_det <- do.call(cbind, lapply(batches, function(b) {
    idx <- which(cd_batch$Batch == b)
    rowSums(!is.na(assay(data)[, idx, drop = FALSE]))
  }))
  pass_all_batches <- apply(batch_det, 1, function(x) all(x >= MIN_BATCH_DETECTION))
  data <- data[pass_all_batches, ]
}
n_after_batchdet <- nrow(data)

# Pre-imputation QC-RSD filter (per batch; pass = acceptable RSD in >= 1 batch)
if (!is.na(QC_RSD_FILTER)) {
  cd_pre  <- as.data.frame(colData(data))
  batches <- unique(cd_pre$Batch)
  rsd_mat <- do.call(cbind, lapply(batches, function(b) {
    idx <- which(cd_pre$QC == "QC" & cd_pre$Batch == b)
    if (length(idx) < 2) return(rep(NA_real_, nrow(data)))
    apply(assay(data)[, idx, drop = FALSE], 1, function(x) {
      x <- x[!is.na(x)]
      med <- median(x)
      if (length(x) < 2 || med <= 0) return(NA_real_)
      mad(x) / med
    })
  }))
  colnames(rsd_mat) <- as.character(batches)
  passes_any <- apply(rsd_mat, 1, function(r) any(is.na(r) | r <= QC_RSD_FILTER))
  data <- data[passes_any, ]
}
n_after_qcrsd <- nrow(data)

# Print and save filter summary
filter_log <- data.frame(
  step = c(
    "Blank filter",
    if (!is.na(low_int_cutoff))
      sprintf("Low-intensity filter (p%.0f >= %.4g)", LOW_INT_PERCENTILE * 100, low_int_cutoff)
    else
      "Low-intensity filter (disabled)",
    sprintf("QC detection (>= %.0f%%)", QC_DETECTION_LIMIT * 100),
    sprintf("Sample detection (>= %.0f%%)", SAMPLE_DETECTION_LIMIT * 100),
    "Zero variance",
    if (MIN_BATCH_DETECTION > 0) sprintf("Per-batch detection (>= %d per batch)", MIN_BATCH_DETECTION) else "Per-batch detection (disabled)",
    if (!is.na(QC_RSD_FILTER)) sprintf("QC-RSD filter (<= %.0f%% in >= 1 batch)", QC_RSD_FILTER * 100) else "QC-RSD filter (disabled)"
  ),
  features_removed = c(
    n_before          - n_after_blank,
    n_after_blank     - n_after_lowint,
    n_after_lowint    - n_after_qc,
    n_after_qc        - n_after_sample,
    n_after_sample    - n_after_zerovar,
    n_after_zerovar   - n_after_batchdet,
    n_after_batchdet  - n_after_qcrsd
  ),
  features_remaining = c(
    n_after_blank, n_after_lowint,
    n_after_qc, n_after_sample, n_after_zerovar, n_after_batchdet, n_after_qcrsd
  )
)
cat("\n--- Pre-filtering summary ---\n")
print(filter_log, row.names = FALSE)
write.csv(filter_log, file.path(interdir, "prefilter_log.csv"), row.names = FALSE)

report_batch_summary(data, file = file.path(interdir, "batch_summary.csv"))

# Run parameters, kept as a data.frame so the same values can be written to
# run_parameters.txt and embedded in each method's Settings sheet.
run_params <- list(
  "Run timestamp"           = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  "CONFIG_FILE"             = if (config_file == "") "(none)" else config_file,
  "INPUT_MODE"              = input_mode,
  "IN_XLSX"                 = if (input_mode == "msdial") in_xlsx else "(n/a)",
  "IN_FEATURE_TABLE"        = if (input_mode == "xcms") in_feature_table else "(n/a)",
  "IN_SAMPLE_SHEET"         = if (input_mode == "xcms") in_sample_sheet else "(n/a)",
  "COLUMN"                  = column,
  "POLARITY"                = polarity,
  "CORRECTION_METHODS"      = paste(CORRECTION_METHODS, collapse = ", "),
  "QC_DETECTION_LIMIT"      = QC_DETECTION_LIMIT,
  "SAMPLE_DETECTION_LIMIT"  = SAMPLE_DETECTION_LIMIT,
  "BLANK_RATIO"             = BLANK_RATIO,
  "LOW_INT_FILTER"          = if (!is.na(LOW_INT_FILTER)) LOW_INT_FILTER else "(disabled)",
  "LOW_INT_FILTER_FRAC"     = LOW_INT_FILTER_FRAC,
  "LOW_INT_PERCENTILE"      = LOW_INT_PERCENTILE,
  "LOW_INT_CUTOFF"          = if (!is.na(low_int_cutoff)) low_int_cutoff else "(disabled)",
  "MIN_QC_SAMPLE_DETECTION" = MIN_QC_SAMPLE_DETECTION,
  "MIN_BATCH_DETECTION"     = MIN_BATCH_DETECTION,
  "QC_RSD_FILTER"           = if (!is.na(QC_RSD_FILTER)) QC_RSD_FILTER else "(disabled)",
  "RSD_THRESHOLD"           = RSD_THRESHOLD,
  "RUV_K"                   = RUV_K,
  "SERRF_NUM"               = SERRF_NUM,
  "LOESS_SPAN"              = LOESS_SPAN,
  "LOESS_SAMPLE_SPAN"       = LOESS_SAMPLE_SPAN,
  "LOESS_SAMPLE_MIN_OBS"    = LOESS_SAMPLE_MIN_OBS,
  "CORDBAT_REF_BATCH"       = if (is.null(CORDBAT_REF_BATCH)) "(auto)" else CORDBAT_REF_BATCH,
  "NORMALIZATION"           = NORMALIZATION,
  "N_CORES"                 = if (n_cores_env == "") paste(parallel::detectCores() - 1, "(auto)") else n_cores_env
)
run_params_df <- data.frame(Parameter = names(run_params),
                            Value     = unlist(run_params, use.names = FALSE),
                            stringsAsFactors = FALSE)

writeLines(sprintf("%-24s%s", paste0(run_params_df$Parameter, ":"), run_params_df$Value),
           file.path(interdir, "run_parameters.txt"))

# Sanity check: injection order must be finite for all samples
bad_inj <- !is.finite(colData(data)$Injection_order)
if (any(bad_inj)) {
  cat("WARNING: samples with non-finite Injection_order:\n")
  print(as.data.frame(colData(data))[bad_inj,
    intersect(c("Sample_ID", "Original_name", "QC", "Batch", "Injection_order"),
              colnames(colData(data)))])
  stop("Non-finite injection orders found — fix the ", input_mode, "_to_notame() conversion before proceeding.")
} else {
  cat("Injection_order: OK\n")
}

# ─────────────────────────────────────────────────────────────────────────────
# 3) CORRECTION & OUTPUT (per method)
# ─────────────────────────────────────────────────────────────────────────────

# Pre-correction QC plots
if (SAVE_PRE_CORRECTION_PLOTS) {
  dir.create(file.path(output_dir, "pre_correction"), showWarnings = FALSE, recursive = TRUE)
  tryCatch(
    save_QC_plots(data, prefix = file.path(output_dir, "pre_correction/"), id = "Sample_ID",
                  perplexity = safe_perplexity(ncol(data))),
    error = function(e) message("WARNING: save_QC_plots failed (pre-correction): ", conditionMessage(e))
  )
}

old_summaries <- list.files(interdir, pattern = "^qc_summary_.+\\.csv$", full.names = TRUE)
if (length(old_summaries) > 0) file.remove(old_summaries)

# Build raw reference for signal-preservation metric (biological samples only).
# Intentionally NOT imputed — NAs are kept so that signal_preservation_r is
# computed only over originally-observed positions. Positions missing in the
# raw data are excluded from the correlation naturally via is.finite() checks.
message("==> Building raw reference for signal-preservation metric")
raw_ref <- tryCatch({
  samp_idx <- which(colData(data)$QC == "Sample")
  mat      <- assay(data, 1)[, samp_idx, drop = FALSE]
  colnames(mat) <- colData(data)$Sample_ID[samp_idx]
  write.csv(as.data.frame(mat), file.path(output_dir, "raw_reference.csv"))
  mat
}, error = function(e) {
  message("WARNING: could not build raw reference (signal_preservation_r will be NA): ", conditionMessage(e))
  NULL
})

# Uncorrected baseline QC metrics
save_correction_summary(assess_quality(data), method = "uncorrected", interdir = interdir, raw_ref = raw_ref)

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(as.data.frame(colData(data)), file.path(output_dir, "sample_metadata.csv"), row.names = FALSE)

run_log <- new_run_log()

for (method in CORRECTION_METHODS) {
  message("\n############################################################")
  message("# METHOD: ", toupper(method))
  message("############################################################")

  method_out   <- file.path(output_dir, method)
  method_start <- Sys.time()
  dir.create(file.path(method_out, "QC_plots"), showWarnings = FALSE, recursive = TRUE)

  result <- tryCatch({
    switch(method,
      none         = correct_none(data),
      notame       = correct_notame(data, RUV_K),
      pmp_qcrsc         = correct_pmp_qcrsc(data),
      pmp_qcrsc_scale         = correct_pmp_qcrsc_scale(data),
      pmp_qcrsc_feature_scale = correct_pmp_qcrsc_feature_scale(data),
      serrf        = {
        # Cap SERRF_NUM at half the smallest per-batch QC count to avoid overfitting
        min_qc_per_batch <- min(table(colData(data)$Batch[colData(data)$QC == "QC"]))
        serrf_num_eff    <- max(1L, min(SERRF_NUM, floor(min_qc_per_batch / 2L)))
        if (serrf_num_eff < SERRF_NUM)
          message("==> SERRF: reducing num from ", SERRF_NUM, " to ", serrf_num_eff,
                  " (min QC per batch = ", min_qc_per_batch, ")")
        correct_serrf(data, num = serrf_num_eff)
      },
      batchcorr    = correct_batchcorr(data),
      combat_only  = correct_combat_only(data),
      loess_combat        = correct_loess_combat(data, LOESS_SPAN),
      loess_samples_combat = correct_loess_samples_combat(data, LOESS_SPAN, LOESS_SAMPLE_SPAN,
                                                           LOESS_SAMPLE_MIN_OBS, LOESS_MIN_QC_PER_BATCH,
                                                           LOESS_MIN_LTQC_VALIDATE,
                                                           LOESS_VALIDATE_SAMPLES_CORRECTION),
      loess_limma   = correct_loess_limma(data, LOESS_SPAN),
      loess_samples_limma = correct_loess_samples_limma(data, LOESS_SPAN, LOESS_SAMPLE_SPAN,
                                                          LOESS_SAMPLE_MIN_OBS, LOESS_MIN_QC_PER_BATCH,
                                                          LOESS_MIN_LTQC_VALIDATE,
                                                          LOESS_VALIDATE_SAMPLES_CORRECTION),
      loess_feature_median = correct_loess_feature_median(data, LOESS_SPAN),
      loess_global_median  = correct_loess_global_median(data, LOESS_SPAN),
      cordbat_only  = correct_cordbat_only(data, CORDBAT_REF_BATCH),
      loess_cordbat = correct_loess_cordbat(data, LOESS_SAMPLE_SPAN, LOESS_SAMPLE_MIN_OBS,
                                             CORDBAT_REF_BATCH),
      waveica      = correct_waveica(data, alpha = WAVEICA_ALPHA, cutoff = WAVEICA_CUTOFF,
                                      K = WAVEICA_K, wf = WAVEICA_WF),
      stop("Unknown method '", method, "'. Valid: none, notame, pmp_qcrsc, pmp_qcrsc_scale, pmp_qcrsc_feature_scale, serrf, batchcorr, combat_only, loess_combat, loess_samples_combat, loess_limma, loess_samples_limma, loess_feature_median, loess_global_median, cordbat_only, loess_cordbat, waveica")
    )
  }, error = function(e) {
    message("ERROR in method '", method, "': ", conditionMessage(e))
    message("Skipping.")
    run_log <<- log_method_result(run_log, method, "failed",
                                   as.numeric(difftime(Sys.time(), method_start, units = "secs")),
                                   error = conditionMessage(e))
    NULL
  })

  if (is.null(result)) next

  combined <- result$post
  obs_mask <- result$obs_mask

  # Optional post-correction normalisation for dilution effects
  if (NORMALIZATION == "pqn") {
    tryCatch({
      suppressPackageStartupMessages(library(pmp))
      classes    <- colData(combined)$QC
      qc_present <- any(classes == "QC")
      combined   <- pqn_normalisation(df = combined, classes = classes,
                                      qc_label = if (qc_present) "QC" else NULL)
      message("==> PQN normalisation applied (reference: ",
              if (qc_present) "pooled QC" else "median of samples", ")")
    }, error = function(e) {
      message("WARNING: PQN normalisation failed: ", conditionMessage(e))
    })
  }

  tryCatch(
    save_QC_plots(combined, prefix = file.path(method_out, "QC_plots/post_correction_"), id = "Sample_ID",
                  perplexity = safe_perplexity(ncol(combined))),
    error = function(e) message("WARNING: save_QC_plots failed (", method, "): ", conditionMessage(e))
  )

  combined <- assess_quality(combined)
  save_correction_summary(combined, method = method, interdir = interdir, obs_mask = obs_mask, raw_ref = raw_ref)
  report_batch_summary(combined, file = file.path(method_out, "batch_summary_post_correction.csv"))

  # Drop RUV W-factor columns from final output
  w_cols <- grep("^W_", colnames(colData(combined)), value = TRUE)
  if (length(w_cols) > 0)
    colData(combined) <- colData(combined)[, !colnames(colData(combined)) %in% w_cols]

  message("==> Writing output: ", method)
  combined <- add_batch_qc_metrics(combined)

  # Settings sheet content for this method: global run parameters plus the
  # method and feature-set context for this particular workbook.
  method_settings_df <- function(feature_set_label) {
    rbind(run_params_df,
          data.frame(Parameter = c("CORRECTION_METHOD", "FEATURE_SET"),
                     Value     = c(method, feature_set_label),
                     stringsAsFactors = FALSE))
  }

  # Helper: write the full (uncompressed) and clustered (compressed) results
  # workbooks for a given SE.
  write_outputs <- function(se, suffix, feature_set_label) {
    settings_df <- method_settings_df(feature_set_label)
    write_results_workbook(se, feature_annotations, settings_df,
                            file = file.path(method_out, paste0("results_full", suffix, ".xlsx")))
    se_c <- tryCatch(compress_clusters(cluster_features(se, all_features = TRUE)),
                     error = function(e) { message("WARNING: clustering failed: ", conditionMessage(e)); NULL })
    if (!is.null(se_c)) {
      write_results_workbook(se_c, feature_annotations, settings_df,
                              file = file.path(method_out, paste0("results_clustered", suffix, ".xlsx")))
    }
  }

  # Full feature set
  write_outputs(combined, "", "All features")

  # Global QC-RSD filter
  rsd_suffix <- paste0("_rsd", round(RSD_THRESHOLD * 100))
  tryCatch({
    global_keep <- !is.na(rowData(combined)$RSD_r) & rowData(combined)$RSD_r < RSD_THRESHOLD
    write_outputs(combined[global_keep, ], rsd_suffix,
                  sprintf("QC RSD < %.0f%%", RSD_THRESHOLD * 100))
  }, error = function(e) message("WARNING: global RSD filter export failed (", method, "): ", conditionMessage(e)))

  # Batchwise QC-RSD filter (pass in >= 50% of batches)
  batch_rsd_cols <- grep("^RSD_r_", colnames(rowData(combined)), value = TRUE)
  if (length(batch_rsd_cols) > 0) {
    tryCatch({
      batch_rsd_mat <- as.matrix(as.data.frame(rowData(combined))[, batch_rsd_cols, drop = FALSE])
      frac_passing  <- rowMeans(batch_rsd_mat < RSD_THRESHOLD, na.rm = TRUE)
      write_outputs(combined[!is.na(frac_passing) & frac_passing >= 0.5, ],
                    paste0("_batchrsd", round(RSD_THRESHOLD * 100)),
                    sprintf("Per-batch QC RSD < %.0f%% in >= 50%% of batches", RSD_THRESHOLD * 100))
    }, error = function(e) message("WARNING: batchwise RSD filter export failed (", method, "): ", conditionMessage(e)))
  }

  run_log <- log_method_result(run_log, method, "success",
                                as.numeric(difftime(Sys.time(), method_start, units = "secs")),
                                n_features = nrow(combined))
}

run_log_df <- write_run_log(run_log, interdir)

compare_corrections(interdir, output_dir)

primary_input <- if (input_mode == "xcms") in_feature_table else in_xlsx
tryCatch(
  write_run_report(interdir, output_dir, mode_label, primary_input, run_params_df, run_log_df),
  error = function(e) message("WARNING: report generation failed: ", conditionMessage(e))
)

message("==> FINISHED. Output at: ", output_dir)

# ─────────────────────────────────────────────────────────────────────────────
# Post-MSDIAL preprocessing pipeline for untargeted LC-MS metabolomics.
# Converts MSDIAL alignment exports to notame format, applies pre-correction
# feature filters, and runs one or more drift/batch correction methods in
# parallel. QC metrics are computed per method for comparison.
#
# Daniel Brunnsåker, 2026-03-20
# ─────────────────────────────────────────────────────────────────────────────

library(notame)
library(notameViz)
library(notameStats)
library(openxlsx)
library(doParallel)

source("R/msdial_to_notame.R")
source("R/qc_metrics.R")
source("R/drift_correction.R")
source("R/correction_methods.R")
source("R/serrf.R")

n_cores_env <- Sys.getenv("N_CORES", unset = "")
registerDoParallel(cores = if (n_cores_env == "") parallel::detectCores() - 1 else as.integer(n_cores_env))

# ─────────────────────────────────────────────────────────────────────────────
# SETTINGS
# All parameters can be overridden via environment variables.
# ─────────────────────────────────────────────────────────────────────────────

get_env <- function(var, default) {
  val <- Sys.getenv(var, unset = NA)
  if (is.na(val) || val == "") default else val
}

if ("--help" %in% commandArgs(trailingOnly = TRUE)) {
  cat("
Usage: Rscript notame-workflow.r [--help]

Environment variables (all optional, hardcoded defaults shown):

  IN_XLSX               Path to the MSDIAL alignment export (.xlsx). Required.

  PROJECT_FOLDER        Root output directory (intermediates/ and output/ written here). Required.

  COLUMN                Chromatographic column type. Required. Used together with POLARITY
                        to namespace output folders (e.g. RP_POS, HILIC_NEG).
                        Examples: RP | HILIC

  POLARITY              Ionisation polarity. Required.
                        Values:  POS | NEG

  CORRECTION_METHODS    Comma-separated list of correction methods to run.
                        Each method is saved to its own output subfolder.
                        Default: none,notame
                        Values:  none | notame | pmp_qcrsc | serrf |
                                 batchcorr | combat_only | loess_combat | loess_limma | waveica

  QC_DETECTION_LIMIT    Min fraction of QC samples a feature must be detected in
                        Default: 0.60

  SAMPLE_DETECTION_LIMIT  Min fraction of biological samples a feature must be detected in
                        Default: 0.20

  BLANK_RATIO           Remove features where mean(Sample) <= BLANK_RATIO * mean(Blank).
                        Set to 'none' to disable.
                        Default: none

  FILL_FILTER           Min MSDIAL Fill % (post-gap-filling detection rate, 0-1).
                        Features below this threshold are removed.
                        Default: 0.10

  LOW_INT_FILTER_FRAC   Data-driven low-intensity filter. Removes features whose p80
                        intensity is below this fraction of the mean p80 across all features.
                        Ignored if LOW_INT_FILTER is set.
                        Default: 0.10

  LOW_INT_FILTER        Absolute low-intensity threshold (overrides LOW_INT_FILTER_FRAC).
                        Default: (disabled)

  LOW_INT_PERCENTILE    Percentile used for the low-intensity filter (0-1).
                        Default: 0.8

  MIN_QC_SAMPLE_DETECTION  Minimum fraction of features that must be detected in a QC sample
                        for it to be used as a QC reference. QC samples below this threshold
                        are removed before processing (empty injections, failed runs).
                        Default: 0.50

  QC_RSD_FILTER         Pre-correction QC-RSD threshold. Features with QC CV above this
                        value in all batches are removed (pass = acceptable in >= 1 batch).
                        Set to 'none' to disable.
                        Default: none (disabled)

  RSD_THRESHOLD         RSD threshold used for post-correction output filtering.
                        Controls both the global (_rsdXX) and per-batch (_batchrsdXX) outputs.
                        Default: 0.30

  RUV_K                 Number of unwanted variation factors for RUV (notame method only).
                        Default: 3

  LOESS_SPAN            LOESS smoothing span for drift correction (loess_combat only).
                        Default: 0.75

  LOESS_FALLBACK_TO_SAMPLES  If TRUE, fall back to fitting LOESS through all samples when a
                        batch has fewer than 4 finite QC observations. Only appropriate when
                        samples are in randomised injection order.
                        Default: FALSE

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

project_folder <- Sys.getenv("PROJECT_FOLDER", unset = "")
if (project_folder == "") stop("PROJECT_FOLDER is required. Set it to the root output directory.")

in_xlsx <- Sys.getenv("IN_XLSX", unset = "")
if (in_xlsx == "") stop("IN_XLSX is required. Set it to the path of your MSDIAL alignment export.")

polarity <- Sys.getenv("POLARITY", unset = "")
if (polarity == "") stop("POLARITY is required. Set it to 'POS' or 'NEG'.")

column <- Sys.getenv("COLUMN", unset = "")
if (column == "") stop("COLUMN is required. Set it to the chromatographic column type, e.g. 'RP' or 'HILIC'.")

mode_label <- paste0(column, "_", polarity)  # e.g. RP_POS — used to namespace all output folders

out_xlsx   <- file.path(project_folder, "intermediates", mode_label, "notame_rev.xlsx")
interdir   <- file.path(project_folder, "intermediates", mode_label)
output_dir <- file.path(project_folder, "output",        mode_label)

QC_DETECTION_LIMIT     <- as.numeric(get_env("QC_DETECTION_LIMIT",     "0.60"))
SAMPLE_DETECTION_LIMIT <- as.numeric(get_env("SAMPLE_DETECTION_LIMIT", "0.20"))

blank_ratio_env <- Sys.getenv("BLANK_RATIO", unset = "")
BLANK_RATIO <- if (blank_ratio_env %in% c("none", "skip")) NA_real_ else
               if (blank_ratio_env == "") 1 else
               as.numeric(blank_ratio_env)

LOW_INT_FILTER      <- suppressWarnings(as.numeric(get_env("LOW_INT_FILTER", "")))
LOW_INT_FILTER_FRAC <- as.numeric(get_env("LOW_INT_FILTER_FRAC", "0.10"))
LOW_INT_PERCENTILE  <- as.numeric(get_env("LOW_INT_PERCENTILE",  "0.8"))
FILL_FILTER         <- as.numeric(get_env("FILL_FILTER",         "0.10"))
qc_rsd_env    <- Sys.getenv("QC_RSD_FILTER", unset = "")
QC_RSD_FILTER           <- if (qc_rsd_env %in% c("none", "")) NA_real_ else as.numeric(qc_rsd_env)
MIN_QC_SAMPLE_DETECTION <- as.numeric(get_env("MIN_QC_SAMPLE_DETECTION", "0.50"))
RSD_THRESHOLD <- as.numeric(get_env("RSD_THRESHOLD", "0.30"))

# Correction methods:
#   "none"         — imputation only (no correction; baseline)
#   "notame"       — per-batch cubic spline drift correction + RUV batch correction
#   "pmp_qcrsc"    — QC-RSC spline drift correction (pmp package)
#   "serrf"        — SERRF random forest correction (Fan et al. 2019)
#   "batchcorr"    — cluster-based spline drift + between-batch normalisation (Brunius et al.)
#   "combat_only"  — ComBat batch correction only (no drift correction)
#   "loess_combat" — per-batch LOESS drift correction + ComBat batch correction
#   "waveica"      — WaveICA 2.0 wavelet-based correction
CORRECTION_METHODS <- strsplit(get_env("CORRECTION_METHODS", "none,notame"), ",")[[1]]

RUV_K      <- as.integer(get_env("RUV_K",       "3"))
LOESS_SPAN                <- as.numeric(get_env("LOESS_SPAN", "0.75"))
LOESS_FALLBACK_TO_SAMPLES <- as.logical(get_env("LOESS_FALLBACK_TO_SAMPLES", "FALSE"))
NORMALIZATION              <- get_env("NORMALIZATION",              "none")
SAVE_PRE_CORRECTION_PLOTS  <- as.logical(get_env("SAVE_PRE_CORRECTION_PLOTS", "TRUE"))

# ─────────────────────────────────────────────────────────────────────────────
# HELPERS
# ─────────────────────────────────────────────────────────────────────────────

# Write a per-method MSDIAL annotations file, filtered to features present in
# the given SE. For clustered SEs, borrows annotation from a cluster member
# when the representative feature itself has no metabolite name.
write_annotations <- function(se, annot_df, file) {
  tryCatch({
    ids <- rownames(se)
    rd  <- as.data.frame(rowData(se), check.names = FALSE)

    name_col     <- "Metabolite name"
    is_annotated <- function(x) !is.na(x) & nchar(trimws(x)) > 0 & !grepl("^Unknown$", x, ignore.case = TRUE)

    out <- annot_df[match(ids, annot_df$Feature_ID), , drop = FALSE]
    out$Feature_ID        <- ids
    rownames(out)         <- NULL
    out$Annotation_source <- "Representative feature"

    if (name_col %in% colnames(out) && "Cluster_features" %in% colnames(rd)) {
      needs <- which(!is_annotated(out[[name_col]]))
      for (i in needs) {
        id          <- ids[i]
        members_raw <- rd[id, "Cluster_features"]
        if (is.na(members_raw) || nchar(trimws(as.character(members_raw))) == 0) next
        member_ids  <- trimws(strsplit(as.character(members_raw), ";")[[1]])
        member_rows <- annot_df[annot_df$Feature_ID %in% member_ids, , drop = FALSE]
        annotated   <- member_rows[is_annotated(member_rows[[name_col]]), , drop = FALSE]
        if (nrow(annotated) == 0) next
        out[i, seq_len(ncol(annot_df))] <- annotated[1, ]
        out$Feature_ID[i]        <- id
        out$Annotation_source[i] <- paste0("Cluster member (", annotated$Feature_ID[1], ")")
      }
    }

    write.xlsx(out, file, colNames = TRUE, rowNames = FALSE)
  }, error = function(e) message("WARNING: could not write annotations to ", file, ": ", conditionMessage(e)))
}

# ─────────────────────────────────────────────────────────────────────────────
# 1) CONVERT & IMPORT
# ─────────────────────────────────────────────────────────────────────────────

dir.create(interdir, showWarnings = FALSE, recursive = TRUE)

msdial_result      <- msdial_to_notame(in_xlsx, out_xlsx)
mode_name          <- msdial_result$mode
msdial_annotations <- msdial_result$annotations

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

# Remove QC samples with insufficient feature detection (empty injections, failed runs)
{
  qc_cols    <- which(colData(data)$QC == "QC")
  mat_qc     <- assay(data)[, qc_cols, drop = FALSE]
  detect_qc  <- colMeans(mat_qc > 0 & !is.na(mat_qc))
  bad_qc     <- qc_cols[detect_qc < MIN_QC_SAMPLE_DETECTION]
  if (length(bad_qc) > 0) {
    bad_names  <- colData(data)$Sample_ID[bad_qc]
    bad_batch  <- colData(data)$Batch[bad_qc]
    message("==> Removing ", length(bad_qc), " QC sample(s) with detection rate < ",
            round(MIN_QC_SAMPLE_DETECTION * 100), "%:")
    for (k in seq_along(bad_names))
      message("    ", bad_names[k], " (batch: ", bad_batch[k], ", detection: ",
              round(detect_qc[detect_qc < MIN_QC_SAMPLE_DETECTION][k] * 100, 1), "%)")
    data <- data[, -bad_qc]
  }
}

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

# Fill % filter: remove features with low MSDIAL alignment confidence
fill_pct <- as.numeric(rowData(data)$Fill_pct)
if (all(is.na(fill_pct))) {
  message("WARNING: Fill_pct not available — skipping fill filter")
  n_after_fill <- n_after_lowint
} else {
  data <- data[is.na(fill_pct) | fill_pct >= FILL_FILTER, ]
  n_after_fill <- nrow(data)
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

# Pre-imputation QC-RSD filter (per batch; pass = acceptable RSD in >= 1 batch)
if (!is.na(QC_RSD_FILTER)) {
  cd_pre  <- as.data.frame(colData(data))
  batches <- unique(cd_pre$Batch)
  rsd_mat <- do.call(cbind, lapply(batches, function(b) {
    idx <- which(cd_pre$QC == "QC" & cd_pre$Batch == b)
    if (length(idx) < 2) return(rep(NA_real_, nrow(data)))
    apply(assay(data)[, idx, drop = FALSE], 1, function(x) {
      x <- x[!is.na(x)]
      if (length(x) < 2 || median(x) == 0) return(NA_real_)
      mad(x) / median(x)
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
    sprintf("Fill filter (>= %.2g)", FILL_FILTER),
    sprintf("QC detection (>= %.0f%%)", QC_DETECTION_LIMIT * 100),
    sprintf("Sample detection (>= %.0f%%)", SAMPLE_DETECTION_LIMIT * 100),
    "Zero variance",
    if (!is.na(QC_RSD_FILTER)) sprintf("QC-RSD filter (<= %.0f%% in >= 1 batch)", QC_RSD_FILTER * 100) else "QC-RSD filter (disabled)"
  ),
  features_removed = c(
    n_before        - n_after_blank,
    n_after_blank   - n_after_lowint,
    n_after_lowint  - n_after_fill,
    n_after_fill    - n_after_qc,
    n_after_qc      - n_after_sample,
    n_after_sample  - n_after_zerovar,
    n_after_zerovar - n_after_qcrsd
  ),
  features_remaining = c(
    n_after_blank, n_after_lowint, n_after_fill,
    n_after_qc, n_after_sample, n_after_zerovar, n_after_qcrsd
  )
)
cat("\n--- Pre-filtering summary ---\n")
print(filter_log, row.names = FALSE)
write.csv(filter_log, file.path(interdir, "prefilter_log.csv"), row.names = FALSE)

# Save run parameters
writeLines(c(
  paste("Run timestamp:          ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste("IN_XLSX:                ", in_xlsx),
  paste("COLUMN:                 ", column),
  paste("POLARITY:               ", polarity),
  paste("CORRECTION_METHODS:     ", paste(CORRECTION_METHODS, collapse = ", ")),
  paste("QC_DETECTION_LIMIT:     ", QC_DETECTION_LIMIT),
  paste("SAMPLE_DETECTION_LIMIT: ", SAMPLE_DETECTION_LIMIT),
  paste("BLANK_RATIO:            ", BLANK_RATIO),
  paste("LOW_INT_FILTER:         ", if (!is.na(LOW_INT_FILTER)) LOW_INT_FILTER else "(disabled)"),
  paste("LOW_INT_FILTER_FRAC:    ", LOW_INT_FILTER_FRAC),
  paste("LOW_INT_PERCENTILE:     ", LOW_INT_PERCENTILE),
  paste("LOW_INT_CUTOFF:         ", if (!is.na(low_int_cutoff)) low_int_cutoff else "(disabled)"),
  paste("FILL_FILTER:            ", FILL_FILTER),
  paste("MIN_QC_SAMPLE_DETECTION: ", MIN_QC_SAMPLE_DETECTION),
  paste("QC_RSD_FILTER:          ", if (!is.na(QC_RSD_FILTER)) QC_RSD_FILTER else "(disabled)"),
  paste("RSD_THRESHOLD:          ", RSD_THRESHOLD),
  paste("RUV_K:                  ", RUV_K),
  paste("LOESS_SPAN:             ", LOESS_SPAN),
  paste("LOESS_FALLBACK_TO_SAMPLES: ", LOESS_FALLBACK_TO_SAMPLES),
  paste("N_CORES:                ", if (n_cores_env == "") paste(parallel::detectCores() - 1, "(auto)") else n_cores_env)
), file.path(interdir, "run_parameters.txt"))

# Sanity check: injection order must be finite for all samples
bad_inj <- !is.finite(colData(data)$Injection_order)
if (any(bad_inj)) {
  cat("WARNING: samples with non-finite Injection_order:\n")
  print(as.data.frame(colData(data))[bad_inj,
    intersect(c("Sample_ID", "Original_name", "QC", "Batch", "Injection_order"),
              colnames(colData(data)))])
  stop("Non-finite injection orders found — fix msdial_to_notame parsing before proceeding.")
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
    save_QC_plots(data, prefix = file.path(output_dir, "pre_correction/"), id = "Sample_ID"),
    error = function(e) message("WARNING: save_QC_plots failed (pre-correction): ", conditionMessage(e))
  )
}

old_summaries <- list.files(interdir, pattern = "^qc_summary_.+\\.csv$", full.names = TRUE)
if (length(old_summaries) > 0) file.remove(old_summaries)

# Build imputed raw reference for signal-preservation metric (biological samples only)
message("==> Building raw reference for signal-preservation metric")
raw_ref <- tryCatch({
  raw_imp  <- impute_rf(data, parallelize = "variables")
  samp_idx <- which(colData(raw_imp)$QC == "Sample")
  mat      <- assay(raw_imp, 1)[, samp_idx, drop = FALSE]
  colnames(mat) <- colData(raw_imp)$Sample_ID[samp_idx]
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

for (method in CORRECTION_METHODS) {
  message("\n############################################################")
  message("# METHOD: ", toupper(method))
  message("############################################################")

  method_out <- file.path(output_dir, method)
  dir.create(file.path(method_out, "QC_plots"), showWarnings = FALSE, recursive = TRUE)

  result <- tryCatch({
    switch(method,
      none         = correct_none(data),
      notame       = correct_notame(data, RUV_K),
      pmp_qcrsc    = correct_pmp_qcrsc(data),
      serrf        = correct_serrf(data),
      batchcorr    = correct_batchcorr(data),
      combat_only  = correct_combat_only(data),
      loess_combat = correct_loess_combat(data, LOESS_SPAN, LOESS_FALLBACK_TO_SAMPLES),
      loess_limma  = correct_loess_limma(data, LOESS_SPAN, LOESS_FALLBACK_TO_SAMPLES),
      waveica      = correct_waveica(data),
      stop("Unknown method '", method, "'. Valid: none, notame, pmp_qcrsc, serrf, batchcorr, combat_only, loess_combat, loess_limma, waveica")
    )
  }, error = function(e) {
    message("ERROR in method '", method, "': ", conditionMessage(e))
    message("Skipping.")
    NULL
  })

  if (is.null(result)) next

  combined <- result$post
  obs_mask <- result$obs_mask

  # Optional post-correction normalisation for dilution effects
  if (NORMALIZATION == "pqn") {
    tryCatch({
      library(pmp)
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
    save_QC_plots(combined, prefix = file.path(method_out, "QC_plots/post_correction_"), id = "Sample_ID"),
    error = function(e) message("WARNING: save_QC_plots failed (", method, "): ", conditionMessage(e))
  )

  combined <- assess_quality(combined)
  save_correction_summary(combined, method = method, interdir = interdir, obs_mask = obs_mask, raw_ref = raw_ref)

  # Drop RUV W-factor columns from final output
  w_cols <- grep("^W_", colnames(colData(combined)), value = TRUE)
  if (length(w_cols) > 0)
    colData(combined) <- colData(combined)[, !colnames(colData(combined)) %in% w_cols]

  message("==> Writing output: ", method)
  combined <- add_batch_qc_metrics(combined)

  # Helper: write one set of output files for a given SE
  write_outputs <- function(se, suffix) {
    write_feature_table(se, file = file.path(method_out, paste0("feature_table_full", suffix, ".xlsx")))
    write_feature_info(se,  file = file.path(method_out, paste0("feature_info_full",  suffix, ".xlsx")))
    write_annotations(se, msdial_annotations, file.path(method_out, paste0("annotations_full", suffix, ".xlsx")))
    se_c <- tryCatch(compress_clusters(cluster_features(se, all_features = TRUE)),
                     error = function(e) { message("WARNING: clustering failed: ", conditionMessage(e)); NULL })
    if (!is.null(se_c)) {
      write_feature_table(se_c, file = file.path(method_out, paste0("feature_table", suffix, ".xlsx")))
      write_feature_info(se_c,  file = file.path(method_out, paste0("feature_info",  suffix, ".xlsx")))
      write_annotations(se_c, msdial_annotations, file.path(method_out, paste0("annotations", suffix, ".xlsx")))
    }
  }

  # Full feature set
  write_outputs(combined, "")

  # Global QC-RSD filter
  rsd_suffix <- paste0("_rsd", round(RSD_THRESHOLD * 100))
  tryCatch({
    global_keep <- !is.na(rowData(combined)$RSD_r) & rowData(combined)$RSD_r < RSD_THRESHOLD
    write_outputs(combined[global_keep, ], rsd_suffix)
  }, error = function(e) message("WARNING: global RSD filter export failed (", method, "): ", conditionMessage(e)))

  # Batchwise QC-RSD filter (pass in >= 50% of batches)
  batch_rsd_cols <- grep("^RSD_r_", colnames(rowData(combined)), value = TRUE)
  if (length(batch_rsd_cols) > 0) {
    tryCatch({
      batch_rsd_mat <- as.matrix(as.data.frame(rowData(combined))[, batch_rsd_cols, drop = FALSE])
      frac_passing  <- rowMeans(batch_rsd_mat < RSD_THRESHOLD, na.rm = TRUE)
      write_outputs(combined[!is.na(frac_passing) & frac_passing >= 0.5, ], paste0("_batchrsd", round(RSD_THRESHOLD * 100)))
    }, error = function(e) message("WARNING: batchwise RSD filter export failed (", method, "): ", conditionMessage(e)))
  }
}

compare_corrections(interdir, output_dir)

message("==> FINISHED. Output at: ", output_dir)

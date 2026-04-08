# ─────────────────────────────────────────────────────────────────────────────
# Preprocessing script for MS-Dial outputs. Modified version
# of the notame workflow. Makes use of notame functionalities 
# for basic QC, but evaluates several different techniques for 
# drift and batch correction.
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

registerDoParallel(cores = parallel::detectCores() - 1)

# ─────────────────────────────────────────────────────────────────────────────
# SETTINGS
# Environment variables override hardcoded defaults when set.
# ─────────────────────────────────────────────────────────────────────────────

get_env <- function(var, default) {
  val <- Sys.getenv(var, unset = NA)
  if (is.na(val) || val == "") default else val
}

if ("--help" %in% commandArgs(trailingOnly = TRUE)) {
  cat("
Usage: Rscript notame-workflow.r [--help]

Environment variables (all optional, hardcoded defaults shown):

  IN_XLSX               Path to the MS-DIAL .xlsx export
                        Default: C:/Projects/PR202/2025-04-17_RP_POS_Pr202.xlsx

  PROJECT_FOLDER        Root output directory (intermediates/ and output/ written here)
                        Default: C:/Projects/PR202/Processed

  POLARITY              Ionisation polarity, used to namespace output subfolders
                        Default: POS
                        Values:  POS | NEG

  CORRECTION_METHODS    Comma-separated list of correction methods to run
                        Default: loess_combat,notame,pmp_qcrsc,waveica,batchcorr
                        Values:  notame | loess_combat | loess_limma |
                                 pmp_qcrsc | combat_only | batchcorr | waveica |
                                 linear_combat | linear_limma | serrf

  QC_DETECTION_LIMIT    Min fraction of QC samples a feature must be detected in
                        Default: 0.60

  SAMPLE_DETECTION_LIMIT  Min fraction of biological samples a feature must be detected in
                        Default: 0.40

  BLANK_RATIO           Remove features where mean(Sample) <= BLANK_RATIO * mean(Blank)
                        Default: 1
                        Set to none to disable blank filtering entirely

  RUV_K                 Number of unwanted-variation factors for RUV (notame method only)
                        Default: 3

  LOESS_SPAN            LOESS smoothing span for drift correction (loess_* methods)
                        Default: 0.75

  SAVE_PRE_CORRECTION_PLOTS  Save QC plots before any correction (can be slow with large datasets)
                        Default: TRUE
                        Values:  TRUE | FALSE

  LOW_INT_FILTER        Remove features where the Nth percentile of abundance (NAs treated as 0)
                        is below this threshold. Alternative to blank filtering when no blanks
                        are available. Set to empty string to skip (default).
                        Default: (disabled)

  LOW_INT_PERCENTILE    Percentile to use for LOW_INT_FILTER (0–1).
                        Default: 0.8

")
  quit(status = 0)
}

project_folder <- get_env("PROJECT_FOLDER", "C:/Projects/PR202/Processed")
polarity       <- get_env("POLARITY",       "POS")   # "POS" or "NEG" — used to namespace all output folders
in_xlsx        <- get_env("IN_XLSX",        "C:/Projects/PR202/2025-04-17_RP_POS_Pr202.xlsx")

out_xlsx   <- file.path(project_folder, "intermediates", polarity, "notame_rev.xlsx")
interdir   <- file.path(project_folder, "intermediates", polarity)

# Output folder sent to the downstream analyst
output_dir <- file.path(project_folder, "output", polarity)

# Minimum fraction of QC samples a feature must be detected in (globally)
QC_DETECTION_LIMIT <- as.numeric(get_env("QC_DETECTION_LIMIT", "0.60")) # Not sure i want to go lower than this due to imputation-issues

# Minimum fraction of biological samples a feature must be detected in (globally)
SAMPLE_DETECTION_LIMIT <- as.numeric(get_env("SAMPLE_DETECTION_LIMIT", "0.20")) # Not sure i want to go lower than this due to imputation-issues

# Blank filter: remove features where mean(Sample) <= BLANK_RATIO * mean(Blank).
# Set BLANK_RATIO=none to skip.
blank_ratio_env <- Sys.getenv("BLANK_RATIO", unset = "")
BLANK_RATIO <- if (blank_ratio_env %in% c("none", "skip")) NA_real_ else
               if (blank_ratio_env == "") 1 else
               as.numeric(blank_ratio_env)

# Low-intensity filter: remove features where the Nth percentile of abundance
# (NAs treated as 0) is below this threshold. 
# Set to NA to skip (default).
LOW_INT_FILTER      <- suppressWarnings(as.numeric(get_env("LOW_INT_FILTER", "")))
LOW_INT_PERCENTILE  <- as.numeric(get_env("LOW_INT_PERCENTILE", "0.8"))

# Correction methods to run. All listed methods are executed and saved to
# separate subfolders. A comparison table is printed at the end.
#   "notame"        — per-batch cubic spline drift correction + RUV
#   "loess_combat"  — per-batch LOESS + ComBat
#   "loess_limma"   — per-batch LOESS + limma::removeBatchEffect
#   "pmp_qcrsc"     — pmp::QCRSC
#   "combat_only"   — ComBat batch correction only, no drift correction (baseline)
#   "batchcorr"     — batchCorr cluster-based spline drift correction + QC normalization (Brunius et al.)
#   "waveica"       — WaveICA2.0 - need to look over the defauls on this one, works horribly atm
#   "linear_combat" — per-batch linear drift correction + ComBat
#   "linear_limma"  — per-batch linear drift correction + limma::removeBatchEffect
#   "serrf"         — Systematic Error Removal using Random Forest (Fan et al. 2019)
CORRECTION_METHODS <- strsplit(get_env("CORRECTION_METHODS", "loess_combat,notame,pmp_qcrsc,waveica,batchcorr"), ",")[[1]]

# Number of unwanted variation factors for RUV (only used by RUV, i.e. notame).
RUV_K <- as.integer(get_env("RUV_K", "3"))

# TODO: Add span-control variable for QC-RSC
# QCRSC_SPAN <- 0

# LOESS span for drift correction (used by loess_combat, loess_limma).
# 0.75 is the default; decrease for tighter fit, increase for smoother.
LOESS_SPAN <- as.numeric(get_env("LOESS_SPAN", "0.75"))

# Whether to save QC plots before any correction (can be slow with large datasets)
SAVE_PRE_CORRECTION_PLOTS <- as.logical(get_env("SAVE_PRE_CORRECTION_PLOTS", "TRUE"))

# ─────────────────────────────────────────────────────────────────────────────
# 1) CONVERT & IMPORT
# ─────────────────────────────────────────────────────────────────────────────
#
# Note that this msdial_to_notame function is a bit hacky. Inspect output.

dir.create(interdir, showWarnings = FALSE, recursive = TRUE)

mode_name <- msdial_to_notame(in_xlsx, out_xlsx)

message("==> Importing")

data <- import_from_excel(
  file  = out_xlsx,
  sheet = 1,
  name  = mode_name
)
names(assays(data)) <- "abundances"
data <- fix_object(data, assay.type = "abundances")

n_imported <- nrow(data)
cat("Imported:", n_imported, "features,", ncol(data), "samples\n")
table(colData(data)$QC)


# ─────────────────────────────────────────────────────────────────────────────
# 2) GLOBAL PRE-FILTERING
# ─────────────────────────────────────────────────────────────────────────────

message("==> Global pre-filtering")

# Mark zeros as missing (MS-DIAL exports 0 for missing values i think?)
data <- mark_nas(data, value = 0)

# Blank filter: keep features where mean sample abundance > BLANK_RATIO * mean blank abundance.
# Maybe remove this? Arguably only works if we have blanks dispersed in the runlist
n_before <- nrow(data)
if (!is.na(BLANK_RATIO)) {
  blank_idx <- which(colData(data)$QC == "Blank")
  if (length(blank_idx) > 0) {
    sample_idx   <- which(colData(data)$QC == "Sample")
    blank_means  <- rowMeans(assay(data)[, blank_idx,  drop = FALSE], na.rm = TRUE)
    sample_means <- rowMeans(assay(data)[, sample_idx, drop = FALSE], na.rm = TRUE)
    keep_blank <- is.na(blank_means) | blank_means == 0 | sample_means > BLANK_RATIO * blank_means
    keep_blank[is.na(keep_blank)] <- FALSE  # if no sample signal and blank present, then remove
    data <- data[keep_blank, ]
  } else {
    cat("No blank samples found — skipping blank filter\n")
  }
}
n_after_blank <- nrow(data)

# Remove non-analytical samples. ltQC is kept so we can use it to evaluate later on.
data <- data[, !colData(data)$QC %in% c("Blank", "Wash", "Cond", "MSe", "MS2", "SST")]

# Low-intensity filter: remove features where the Nth percentile of abundance
# across biological + QC samples is below LOW_INT_FILTER.
if (!is.na(LOW_INT_FILTER)) {
  mat_int <- assay(data)
  mat_int[is.na(mat_int)] <- 0
  int_quantile <- apply(mat_int, 1, quantile, probs = LOW_INT_PERCENTILE)
  data <- data[int_quantile >= LOW_INT_FILTER, ]
  cat(sprintf("Low-intensity filter (p%.0f < %.4g): removed %d features\n",
              LOW_INT_PERCENTILE * 100, LOW_INT_FILTER, sum(int_quantile < LOW_INT_FILTER)))
}
n_after_lowint <- nrow(data)

# Remove features with low detection in QC samples
data <- flag_detection(data, qc_limit = QC_DETECTION_LIMIT)
data <- drop_flagged(data)
n_after_qc <- nrow(data)

# Remove features with low global detection in biological samples
sample_idx  <- which(colData(data)$QC == "Sample")
detect_rate <- rowMeans(!is.na(assay(data)[, sample_idx, drop = FALSE]))
data        <- data[detect_rate >= SAMPLE_DETECTION_LIMIT, ]
n_after_sample <- nrow(data)

# Remove zero-variance features
zero_var <- apply(assay(data), 1, function(x) {
  v <- var(x, na.rm = TRUE)
  !is.na(v) && v < .Machine$double.eps
})
data <- data[!zero_var, ]
n_after_zerovar <- nrow(data)

# Save pre-filtering summary
filter_log <- data.frame(
  step               = c("Blank filter", "Low-intensity filter", "QC detection", "Sample detection", "Zero variance"),
  features_removed   = c(n_before      - n_after_blank,
                         n_after_blank - n_after_lowint,
                         n_after_lowint - n_after_qc,
                         n_after_qc    - n_after_sample,
                         n_after_sample - n_after_zerovar),
  features_remaining = c(n_after_blank, n_after_lowint, n_after_qc, n_after_sample, n_after_zerovar)
)
cat("\n--- Pre-filtering summary ---\n")
print(filter_log, row.names = FALSE)
write.csv(filter_log, file.path(interdir, "prefilter_log.csv"), row.names = FALSE)

# Sanity check: injection order must be finite for all samples
bad_inj <- !is.finite(colData(data)$Injection_order)
if (any(bad_inj)) {
  cat("WARNING: samples with non-finite Injection_order (NA/NaN/Inf):\n")
  print(as.data.frame(colData(data))[bad_inj,
    intersect(c("Sample_ID", "Original_name", "QC", "Batch", "Injection_order"),
              colnames(colData(data)))])
  stop("Non-finite injection orders found — fix msdial_to_notame parsing before proceeding.")
} else {
  cat("Injection_order: OK (all finite)\n")
}

# Pre-correction QC plots
if (SAVE_PRE_CORRECTION_PLOTS) {
  dir.create(file.path(output_dir, "pre_correction"), showWarnings = FALSE, recursive = TRUE)
  tryCatch(
    save_QC_plots(data, prefix = file.path(output_dir, "pre_correction/"), id = "Sample_ID"),
    error = function(e) message("WARNING: save_QC_plots failed (pre-correction): ", conditionMessage(e))
  )
}

# ─────────────────────────────────────────────────────────────────────────────
# 3) CORRECTION & OUTPUT (per method)
# ─────────────────────────────────────────────────────────────────────────────

old_summaries <- list.files(interdir, pattern = "^qc_summary_.+\\.csv$", full.names = TRUE)
if (length(old_summaries) > 0) file.remove(old_summaries)

# Build raw reference for signal-preservation metric (Option B):
# Impute the unfiltered data once so we have a complete sample matrix to
# correlate against after each correction. Only biological samples are used.
# Stored as features × samples matrix (Sample_ID as column names).
message("==> Building raw reference for signal-preservation metric")
raw_ref <- tryCatch({
  raw_imp  <- impute_rf(data, parallelize = "variables")
  samp_idx <- which(colData(raw_imp)$QC == "Sample")
  mat      <- assay(raw_imp, 1)[, samp_idx, drop = FALSE]
  colnames(mat) <- colData(raw_imp)$Sample_ID[samp_idx]
  # Persist to disk so posthoc scripts can reuse it
  write.csv(as.data.frame(mat), file.path(output_dir, "raw_reference.csv"))
  mat
}, error = function(e) {
  message("WARNING: could not build raw reference (signal_preservation_r will be NA): ", conditionMessage(e))
  NULL
})

# Save uncorrected baseline QC metrics for comparison
save_correction_summary(assess_quality(data), method = "uncorrected", interdir = interdir, raw_ref = raw_ref)

# Export sample metadata
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(as.data.frame(colData(data)), file.path(output_dir, "sample_metadata.csv"), row.names = FALSE)

for (method in CORRECTION_METHODS) {
  message("\n############################################################")
  message("# METHOD: ", toupper(method))
  message("############################################################")

  method_out <- file.path(output_dir, method)
  dir.create(file.path(method_out, "QC_plots"), showWarnings = FALSE, recursive = TRUE)

  # Correction ──────────────────────────────────────────────────────────────
  result <- tryCatch({
    switch(method,
      notame        = correct_notame(data, RUV_K),
      loess_combat  = correct_loess_combat(data, LOESS_SPAN),
      loess_limma   = correct_loess_limma(data, LOESS_SPAN),
      pmp_qcrsc     = correct_pmp_qcrsc(data),
      combat_only   = correct_combat_only(data),
      batchcorr     = correct_batchcorr(data),
      waveica       = correct_waveica(data),
      linear_combat = correct_linear_combat(data),
      linear_limma  = correct_linear_limma(data),
      serrf         = correct_serrf(data),
      stop("Unknown method '", method, "'. Valid: notame, loess_combat, loess_limma, pmp_qcrsc, combat_only, batchcorr, waveica, linear_combat, linear_limma, serrf")
    )
  }, error = function(e) {
    message("ERROR in method '", method, "': ", conditionMessage(e))
    message("Skipping this method.")
    NULL
  })

  if (is.null(result)) next

  combined <- result$post
  obs_mask <- result$obs_mask  # pre-imputation missingness mask (may be NULL)

  # QC plots: post-correction only
  tryCatch(
    save_QC_plots(combined, prefix = file.path(method_out, "QC_plots/post_correction_"), id = "Sample_ID"),
    error = function(e) message("WARNING: save_QC_plots failed (post-correction, ", method, "): ", conditionMessage(e))
  )

  combined <- assess_quality(combined)
  save_correction_summary(combined, method = method, interdir = interdir, obs_mask = obs_mask, raw_ref = raw_ref)

  # Remove ltQC, not needed for final output
  combined <- combined[, colData(combined)$QC != "ltQC"]

  # Remove RUV W-factor columns, not needed for final output
  w_cols <- grep("^W_", colnames(colData(combined)), value = TRUE)
  if (length(w_cols) > 0)
    colData(combined) <- colData(combined)[, !colnames(colData(combined)) %in% w_cols]

  # Cluster and write final output
  message("==> Clustering and writing output: ", method)

  combined  <- add_batch_qc_metrics(combined)

  # Save full (unclustered) peak table with QC metrics
  write_feature_table(combined,  file = file.path(method_out, "feature_table_full.xlsx"))
  write_feature_info(combined,   file = file.path(method_out, "feature_info_full.xlsx"))

  clustered <- cluster_features(combined, all_features = TRUE)
  clustered  <- compress_clusters(clustered)

  write_feature_table(clustered, file = file.path(method_out, "feature_table.xlsx"))
  write_feature_info(clustered,  file = file.path(method_out, "feature_info.xlsx"))

  # Global RSD filter: features with RSD_r < 0.3 across all QC samples
  global_keep  <- !is.na(rowData(combined)$RSD_r) & rowData(combined)$RSD_r < 0.3
  combined_g   <- combined[global_keep, ]
  tryCatch({
    clustered_g  <- compress_clusters(cluster_features(combined_g, all_features = TRUE))
    write_feature_table(combined_g,  file = file.path(method_out, "feature_table_full_rsd30.xlsx"))
    write_feature_info(combined_g,   file = file.path(method_out, "feature_info_full_rsd30.xlsx"))
    write_feature_table(clustered_g, file = file.path(method_out, "feature_table_rsd30.xlsx"))
    write_feature_info(clustered_g,  file = file.path(method_out, "feature_info_rsd30.xlsx"))
  }, error = function(e) message("WARNING: global RSD filter export failed (", method, "): ", conditionMessage(e)))

  # Batchwise RSD filter: features with RSD_r < 0.3 in >= 50% of batches
  batch_rsd_cols <- grep("^RSD_r_B", colnames(rowData(combined)), value = TRUE)
  if (length(batch_rsd_cols) > 0) {
    batch_rsd_mat  <- as.matrix(as.data.frame(rowData(combined))[, batch_rsd_cols, drop = FALSE])
    frac_passing   <- rowMeans(batch_rsd_mat < 0.3, na.rm = TRUE)
    combined_b     <- combined[!is.na(frac_passing) & frac_passing >= 0.5, ]
    tryCatch({
      clustered_b    <- compress_clusters(cluster_features(combined_b, all_features = TRUE))
      write_feature_table(combined_b,  file = file.path(method_out, "feature_table_full_batchrsd30.xlsx"))
      write_feature_info(combined_b,   file = file.path(method_out, "feature_info_full_batchrsd30.xlsx"))
      write_feature_table(clustered_b, file = file.path(method_out, "feature_table_batchrsd30.xlsx"))
      write_feature_info(clustered_b,  file = file.path(method_out, "feature_info_batchrsd30.xlsx"))
    }, error = function(e) message("WARNING: batchwise RSD filter export failed (", method, "): ", conditionMessage(e)))
  }
}

compare_corrections(interdir, output_dir)

message("==> FINISHED. Output at: ", output_dir)

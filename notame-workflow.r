# ─────────────────────────────────────────────────────────────────────────────
# Preprocessing script for MS-Dial outputs. Modified version
# of the notame workflow. Makes use of notame functionalities 
# for basic QC, but evaluates several different techniques for 
# drift and batch correction.
# Daniel Brunnsåker, 2026-13-15
# ─────────────────────────────────────────────────────────────────────────────

library(notame)
library(notameViz)
library(notameStats)
library(openxlsx)

source("R/msdial_to_notame.R")
source("R/qc_metrics.R")
source("R/drift_correction.R")
source("R/correction_methods.R")


# ─────────────────────────────────────────────────────────────────────────────
# SETTINGS
# ─────────────────────────────────────────────────────────────────────────────

polarity <- "NEG"   # "POS" or "NEG" — used to namespace all output folders
#polarity <- "POS"   # "POS" or "NEG" — used to namespace all output folders
#in_xlsx <- "C:/Projects/NAPFL/MSDIAL-exports/Area_0_2026_03_05_07_47_29.xlsx" # positive
in_xlsx <- "C:/Projects/NAPFL/MSDIAL-exports/Area_0_2026_03_05_07_50_45.xlsx" # negative

out_xlsx   <- file.path("intermediates", polarity, "notame_rev.xlsx")
interdir   <- file.path("intermediates", polarity)

# TODO: Change this to be an external folder, not included in the repo
# project_folder = 

# Minimum fraction of QC samples a feature must be detected in (globally)
QC_DETECTION_LIMIT <- 0.60 # Not sure i want to go lower than this due to imputation-issues

# Minimum fraction of biological samples a feature must be detected in (globally)
SAMPLE_DETECTION_LIMIT <- 0.40 # Not sure i want to go lower than this due to imputation-issues 

# Blank filter: remove features where mean(Sample) <= BLANK_RATIO * mean(Blank).
# Set to NULL to skip.
BLANK_RATIO <- 1 # This is very low, should default to 3x more generally?

# Correction methods to run. All listed methods are executed and saved to
# separate subfolders. A comparison table is printed at the end.
#   "notame"        — per-batch cubic spline drift correction + RUV
#   "loess_combat"  — per-batch LOESS + ComBat
#   "loess_limma"   — per-batch LOESS + limma::removeBatchEffect
#   "pmp_qcrsc"     — pmp::QCRSC
#   "waveica"       — WaveICA2.0 - need to look over the defauls on this one, works horribly atm
#   "linear_combat" — per-batch linear drift correction + ComBat
#   "linear_limma"  — per-batch linear drift correction + limma::removeBatchEffect
CORRECTION_METHODS <- c("loess_combat", "notame", "pmp_qcrsc", "waveica")

# Number of unwanted variation factors for RUV (only used by RUV, i.e. notame).
RUV_K <- 3

# TODO: Add span-control variable for QC-RSC
# QCRSC_SPAN <- 0

# LOESS span for drift correction (used by loess_combat, loess_limma).
# 0.75 is the default; decrease for tighter fit, increase for smoother.
LOESS_SPAN <- 0.75

# Output folder sent to the downstream analyst
output_dir <- file.path("output", polarity)


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
n_before <- nrow(data)
if (!is.null(BLANK_RATIO)) {
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
  step               = c("Blank filter", "QC detection", "Sample detection", "Zero variance"),
  features_removed   = c(n_before - n_after_blank, n_after_blank - n_after_qc,
                         n_after_qc - n_after_sample, n_after_sample - n_after_zerovar),
  features_remaining = c(n_after_blank, n_after_qc, n_after_sample, n_after_zerovar)
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
dir.create(file.path(output_dir, "pre_correction"), showWarnings = FALSE, recursive = TRUE)
save_QC_plots(
  data,
  prefix = file.path(output_dir, "pre_correction/"),
  id     = "Sample_ID"
)


# ─────────────────────────────────────────────────────────────────────────────
# 3) CORRECTION & OUTPUT (per method)
# ─────────────────────────────────────────────────────────────────────────────

old_summaries <- list.files(interdir, pattern = "^qc_summary_.+\\.csv$", full.names = TRUE) # if there are old summaries, remove
if (length(old_summaries) > 0) file.remove(old_summaries)

for (method in CORRECTION_METHODS) {
  message("\n############################################################")
  message("# METHOD: ", toupper(method))
  message("############################################################")

  method_out <- file.path(output_dir, method)
  dir.create(file.path(method_out, "QC_plots"), showWarnings = FALSE, recursive = TRUE)

  # ── Correction ──────────────────────────────────────────────────────────────
  result <- tryCatch({
    switch(method,
      notame        = correct_notame(data, RUV_K),
      loess_combat  = correct_loess_combat(data, LOESS_SPAN),
      loess_limma   = correct_loess_limma(data, LOESS_SPAN),
      pmp_qcrsc     = correct_pmp_qcrsc(data),
      waveica       = correct_waveica(data),
      linear_combat = correct_linear_combat(data),
      linear_limma  = correct_linear_limma(data),
      stop("Unknown method '", method, "'. Valid: notame, loess_combat, loess_limma, pmp_qcrsc, waveica, linear_combat, linear_limma")
    )
  }, error = function(e) {
    message("ERROR in method '", method, "': ", conditionMessage(e))
    message("Skipping this method.")
    NULL
  })

  if (is.null(result)) next

  pre            <- result$pre
  combined       <- result$post
  is_single_step <- identical(pre, combined)  # TRUE for pmp_qcrsc, waveica

  # Save pre / post batch correction data 
  write_to_excel(combined, file = file.path(method_out, "data_post_batch.xlsx"))
  if (!is_single_step)
    write_to_excel(pre, file = file.path(method_out, "data_pre_batch.xlsx"))

  # QC plots: pre and post batch correction 
  save_QC_plots(combined, prefix = file.path(method_out, "QC_plots/post_batch_"), id = "Sample_ID")
  if (!is_single_step)
    save_QC_plots(pre, prefix = file.path(method_out, "QC_plots/pre_batch_"), id = "Sample_ID")

  combined <- assess_quality(combined)
  save_correction_summary(combined, method = method, interdir = interdir)

  # Remove ltQC, not needed for final output
  combined <- combined[, colData(combined)$QC != "ltQC"]

  # Remove RUV W-factor columns, not needed for final output
  w_cols <- grep("^W_", colnames(colData(combined)), value = TRUE)
  if (length(w_cols) > 0)
    colData(combined) <- colData(combined)[, !colnames(colData(combined)) %in% w_cols]

  # Cluster and write final output
  message("==> Clustering and writing output: ", method)

  clustered <- cluster_features(combined, all_features = TRUE)
  clustered  <- compress_clusters(clustered)

  write_feature_table(clustered, file = file.path(method_out, "feature_table.xlsx"))
  write_feature_info(clustered,  file = file.path(method_out, "feature_info.xlsx"))
}

# Print ranked comparison of all methods
# TODO: save this as a CSV or something
compare_corrections(interdir)

message("==> FINISHED. Output at: ", output_dir)

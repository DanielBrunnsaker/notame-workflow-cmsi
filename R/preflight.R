# ─────────────────────────────────────────────────────────────────────────────
# Preflight validation
#
# Validates configuration and inputs before the slow MSDIAL conversion/import
# steps. All problems are collected and reported together, so a typo'd env
# var fails in seconds with a clear message instead of after minutes of
# processing.
# ─────────────────────────────────────────────────────────────────────────────

VALID_CORRECTION_METHODS <- c(
  "none", "notame", "pmp_qcrsc", "pmp_qcrsc_scale", "pmp_qcrsc_feature_scale",
  "serrf", "batchcorr", "combat_only", "loess_combat", "loess_limma",
  "loess_feature_median", "loess_global_median", "cordbat_only",
  "loess_cordbat", "waveica"
)

# Appends a range-check error for `value` to `problems` if invalid.
# NA is reported as an error unless allow_na = TRUE (used for parameters
# that are NA when the corresponding feature is disabled).
check_numeric <- function(problems, name, value, min = NULL, max = NULL, allow_na = FALSE) {
  if (is.na(value)) {
    if (!allow_na) problems <- c(problems, paste0(name, " is not a valid number"))
    return(problems)
  }
  if (!is.null(min) && value < min)
    problems <- c(problems, paste0(name, " (", value, ") must be >= ", min))
  if (!is.null(max) && value > max)
    problems <- c(problems, paste0(name, " (", value, ") must be <= ", max))
  problems
}

run_preflight_checks <- function(in_xlsx, project_folder, column, polarity,
                                  correction_methods, normalization,
                                  qc_detection_limit, sample_detection_limit,
                                  fill_filter, low_int_filter_frac, low_int_percentile,
                                  min_qc_sample_detection, min_batch_detection,
                                  rsd_threshold, ruv_k, serrf_num, loess_span,
                                  blank_ratio, low_int_filter, qc_rsd_filter,
                                  loess_fallback_to_samples, save_pre_correction_plots) {
  problems <- character(0)

  # Input/output paths
  if (!file.exists(in_xlsx)) {
    problems <- c(problems, paste0("IN_XLSX does not exist: ", in_xlsx))
  } else if (!grepl("\\.xlsx$", in_xlsx, ignore.case = TRUE)) {
    problems <- c(problems, paste0("IN_XLSX must be an .xlsx file: ", in_xlsx))
  }

  if (!dir.exists(project_folder)) {
    ok <- tryCatch({ dir.create(project_folder, recursive = TRUE); dir.exists(project_folder) },
                   error = function(e) FALSE)
    if (!ok)
      problems <- c(problems, paste0("PROJECT_FOLDER does not exist and could not be created: ", project_folder))
  }

  # Required identifiers
  if (is.na(column) || column == "")
    problems <- c(problems, "COLUMN is empty")
  if (!polarity %in% c("POS", "NEG"))
    problems <- c(problems, paste0("POLARITY must be 'POS' or 'NEG', got: '", polarity, "'"))

  # Correction methods
  unknown_methods <- setdiff(correction_methods, VALID_CORRECTION_METHODS)
  if (length(unknown_methods) > 0)
    problems <- c(problems, paste0(
      "Unknown CORRECTION_METHODS value(s): ", paste(unknown_methods, collapse = ", "),
      ". Valid: ", paste(VALID_CORRECTION_METHODS, collapse = ", ")))

  # Normalisation
  if (!normalization %in% c("none", "pqn"))
    problems <- c(problems, paste0("NORMALIZATION must be 'none' or 'pqn', got: '", normalization, "'"))

  # Numeric parameters: required fractions (0-1)
  problems <- check_numeric(problems, "QC_DETECTION_LIMIT",      qc_detection_limit,      min = 0, max = 1)
  problems <- check_numeric(problems, "SAMPLE_DETECTION_LIMIT",  sample_detection_limit,  min = 0, max = 1)
  problems <- check_numeric(problems, "FILL_FILTER",             fill_filter,             min = 0, max = 1)
  problems <- check_numeric(problems, "LOW_INT_FILTER_FRAC",     low_int_filter_frac,     min = 0, max = 1)
  problems <- check_numeric(problems, "LOW_INT_PERCENTILE",      low_int_percentile,      min = 0, max = 1)
  problems <- check_numeric(problems, "MIN_QC_SAMPLE_DETECTION", min_qc_sample_detection, min = 0, max = 1)
  problems <- check_numeric(problems, "RSD_THRESHOLD",           rsd_threshold,           min = 0)
  problems <- check_numeric(problems, "LOESS_SPAN",              loess_span,              min = 0, max = 1)

  # Numeric parameters: required integers
  problems <- check_numeric(problems, "MIN_BATCH_DETECTION", min_batch_detection, min = 0)
  problems <- check_numeric(problems, "RUV_K",               ruv_k,               min = 1)
  problems <- check_numeric(problems, "SERRF_NUM",           serrf_num,           min = 1)

  # Numeric parameters: optional (NA means disabled)
  problems <- check_numeric(problems, "BLANK_RATIO",    blank_ratio,    min = 0,           allow_na = TRUE)
  problems <- check_numeric(problems, "LOW_INT_FILTER", low_int_filter, min = 0,           allow_na = TRUE)
  problems <- check_numeric(problems, "QC_RSD_FILTER",  qc_rsd_filter,  min = 0, max = 1,   allow_na = TRUE)

  # Logical parameters
  if (is.na(loess_fallback_to_samples)) problems <- c(problems, "LOESS_FALLBACK_TO_SAMPLES must be TRUE or FALSE")
  if (is.na(save_pre_correction_plots)) problems <- c(problems, "SAVE_PRE_CORRECTION_PLOTS must be TRUE or FALSE")

  if (length(problems) > 0) {
    stop("Preflight checks failed (", length(problems), "):\n  - ",
         paste(problems, collapse = "\n  - "), call. = FALSE)
  }

  message("==> Preflight checks passed")
  invisible(TRUE)
}

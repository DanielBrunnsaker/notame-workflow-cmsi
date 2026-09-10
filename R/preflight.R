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
  "serrf", "batchcorr", "combat_only", "loess_combat", "loess_samples_combat", "auto_combat",
  "loess_limma", "loess_samples_limma", "loess_feature_median", "loess_global_median", "cordbat_only",
  "loess_cordbat", "waveica", "waveica_v1"
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

# Vector version for comma-separated parameter lists (e.g. AUTO_LOESS_SPANS):
# each element is range-checked individually via check_numeric(); an empty
# or all-NA vector (e.g. from a malformed/empty env var) is its own error,
# since check_numeric() alone can't tell "no candidates" from "one bad one".
check_numeric_list <- function(problems, name, values, min = NULL, max = NULL) {
  if (length(values) == 0 || all(is.na(values))) {
    problems <- c(problems, paste0(name, " must be a non-empty comma-separated list of numbers"))
    return(problems)
  }
  for (v in values) problems <- check_numeric(problems, name, v, min = min, max = max)
  problems
}

run_preflight_checks <- function(input_mode, in_xlsx, in_feature_table, in_sample_sheet,
                                  project_folder, column, polarity,
                                  correction_methods, normalization,
                                  qc_detection_limit, sample_detection_limit,
                                  low_int_filter_frac, low_int_percentile,
                                  min_qc_sample_detection, min_batch_detection,
                                  rsd_threshold, ruv_k, serrf_num, loess_span,
                                  loess_sample_span, loess_sample_min_obs,
                                  loess_min_qc_per_batch, loess_min_ltqc_validate,
                                  auto_loess_spans, auto_huber_ks,
                                  auto_sample_loess_spans, auto_sample_huber_ks,
                                  auto_min_qc_per_batch, auto_min_ltqc_validate, auto_min_cv_obs,
                                  waveica_alpha, waveica_cutoff, waveica_k,
                                  waveica_v1_k, waveica_v1_t, waveica_v1_t2, waveica_v1_alpha,
                                  blank_ratio, low_int_filter, qc_rsd_filter,
                                  save_pre_correction_plots,
                                  config_file = "", raw_sample_type_rules = NULL) {
  problems <- character(0)

  # Config file
  if (config_file != "") {
    if (!file.exists(config_file)) {
      problems <- c(problems, paste0("CONFIG_FILE does not exist: ", config_file))
    } else if (is.null(tryCatch(yaml::read_yaml(config_file), error = function(e) NULL))) {
      problems <- c(problems, paste0("CONFIG_FILE could not be parsed as YAML: ", config_file))
    }
  }

  # sample_type_rules (raw, as read from YAML — list of pattern/type pairs)
  for (i in seq_along(raw_sample_type_rules)) {
    r <- raw_sample_type_rules[[i]]
    pattern <- r$pattern
    type    <- r$type
    if (is.null(pattern) || is.null(type) ||
        nchar(trimws(as.character(pattern))) == 0 || nchar(trimws(as.character(type))) == 0) {
      problems <- c(problems, paste0("sample_type_rules[", i, "] needs a non-empty 'pattern' and 'type'"))
      next
    }
    if (inherits(tryCatch(grepl(pattern, ""), error = function(e) e), "error"))
      problems <- c(problems, paste0("sample_type_rules[", i, "] has an invalid regex pattern: '", pattern, "'"))
  }

  # Input/output paths
  if (input_mode == "msdial") {
    if (!file.exists(in_xlsx)) {
      problems <- c(problems, paste0("IN_XLSX does not exist: ", in_xlsx))
    } else if (!grepl("\\.xlsx$", in_xlsx, ignore.case = TRUE)) {
      problems <- c(problems, paste0("IN_XLSX must be an .xlsx file: ", in_xlsx))
    }
  } else if (input_mode == "xcms") {
    if (!file.exists(in_feature_table)) {
      problems <- c(problems, paste0("IN_FEATURE_TABLE does not exist: ", in_feature_table))
    } else if (!grepl("\\.csv$", in_feature_table, ignore.case = TRUE)) {
      problems <- c(problems, paste0("IN_FEATURE_TABLE must be a .csv file: ", in_feature_table))
    }
    if (!file.exists(in_sample_sheet)) {
      problems <- c(problems, paste0("IN_SAMPLE_SHEET does not exist: ", in_sample_sheet))
    } else if (!grepl("\\.xlsx$", in_sample_sheet, ignore.case = TRUE)) {
      problems <- c(problems, paste0("IN_SAMPLE_SHEET must be an .xlsx file: ", in_sample_sheet))
    }
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
  problems <- check_numeric(problems, "LOW_INT_FILTER_FRAC",     low_int_filter_frac,     min = 0, max = 1)
  problems <- check_numeric(problems, "LOW_INT_PERCENTILE",      low_int_percentile,      min = 0, max = 1)
  problems <- check_numeric(problems, "MIN_QC_SAMPLE_DETECTION", min_qc_sample_detection, min = 0, max = 1)
  problems <- check_numeric(problems, "RSD_THRESHOLD",           rsd_threshold,           min = 0)
  problems <- check_numeric(problems, "LOESS_SPAN",              loess_span,              min = 0, max = 1)
  problems <- check_numeric(problems, "LOESS_SAMPLE_SPAN",       loess_sample_span,       min = 0, max = 1)

  # Numeric parameters: required integers
  problems <- check_numeric(problems, "MIN_BATCH_DETECTION",    min_batch_detection,    min = 0)
  problems <- check_numeric(problems, "RUV_K",                  ruv_k,                  min = 1)
  problems <- check_numeric(problems, "SERRF_NUM",               serrf_num,              min = 1)
  problems <- check_numeric(problems, "LOESS_SAMPLE_MIN_OBS",    loess_sample_min_obs,   min = 4)
  problems <- check_numeric(problems, "LOESS_MIN_QC_PER_BATCH",  loess_min_qc_per_batch, min = 1)
  problems <- check_numeric(problems, "LOESS_MIN_LTQC_VALIDATE", loess_min_ltqc_validate, min = 2)
  problems <- check_numeric_list(problems, "AUTO_LOESS_SPANS",        auto_loess_spans,        min = 0, max = 1)
  problems <- check_numeric_list(problems, "AUTO_HUBER_KS",           auto_huber_ks,           min = 0)
  problems <- check_numeric_list(problems, "AUTO_SAMPLE_LOESS_SPANS", auto_sample_loess_spans, min = 0, max = 1)
  problems <- check_numeric_list(problems, "AUTO_SAMPLE_HUBER_KS",    auto_sample_huber_ks,    min = 0)
  problems <- check_numeric(problems, "AUTO_MIN_QC_PER_BATCH",  auto_min_qc_per_batch,  min = 1)
  problems <- check_numeric(problems, "AUTO_MIN_LTQC_VALIDATE", auto_min_ltqc_validate, min = 2)
  problems <- check_numeric(problems, "AUTO_MIN_CV_OBS",        auto_min_cv_obs,        min = 4)
  problems <- check_numeric(problems, "WAVEICA_ALPHA",           waveica_alpha,  min = 0, max = 1)
  problems <- check_numeric(problems, "WAVEICA_CUTOFF",          waveica_cutoff, min = 0, max = 1)
  problems <- check_numeric(problems, "WAVEICA_K",               waveica_k,      min = 1, allow_na = TRUE)
  problems <- check_numeric(problems, "WAVEICA_V1_K",            waveica_v1_k,     min = 1)
  problems <- check_numeric(problems, "WAVEICA_V1_T",            waveica_v1_t,     min = 0, max = 1)
  problems <- check_numeric(problems, "WAVEICA_V1_T2",           waveica_v1_t2,    min = 0, max = 1)
  problems <- check_numeric(problems, "WAVEICA_V1_ALPHA",        waveica_v1_alpha, min = 0, max = 1)

  # Numeric parameters: optional (NA means disabled)
  problems <- check_numeric(problems, "BLANK_RATIO",    blank_ratio,    min = 0,           allow_na = TRUE)
  problems <- check_numeric(problems, "LOW_INT_FILTER", low_int_filter, min = 0,           allow_na = TRUE)
  problems <- check_numeric(problems, "QC_RSD_FILTER",  qc_rsd_filter,  min = 0, max = 1,   allow_na = TRUE)

  # Logical parameters
  if (is.na(save_pre_correction_plots)) problems <- c(problems, "SAVE_PRE_CORRECTION_PLOTS must be TRUE or FALSE")

  if (length(problems) > 0) {
    stop("Preflight checks failed (", length(problems), "):\n  - ",
         paste(problems, collapse = "\n  - "), call. = FALSE)
  }

  message("==> Preflight checks passed")
  invisible(TRUE)
}

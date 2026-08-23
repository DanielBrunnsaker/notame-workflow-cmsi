# ─────────────────────────────────────────────────────────────────────────────
# Config file support
#
# Optional YAML file (CONFIG_FILE env var) supplying default values for any
# pipeline parameter. Env vars still take precedence over the config file,
# so Docker one-off overrides keep working unchanged; see get_env() in
# notame-workflow.r.
# ─────────────────────────────────────────────────────────────────────────────

KNOWN_CONFIG_KEYS <- c(
  "IN_XLSX", "IN_FEATURE_TABLE", "IN_SAMPLE_SHEET",
  "PROJECT_FOLDER", "FORCE_RECONVERT", "COLUMN", "POLARITY",
  "CORRECTION_METHODS", "QC_DETECTION_LIMIT", "SAMPLE_DETECTION_LIMIT",
  "BLANK_RATIO", "LOW_INT_FILTER", "LOW_INT_FILTER_FRAC", "LOW_INT_PERCENTILE",
  "MIN_QC_SAMPLE_DETECTION", "MIN_BATCH_DETECTION", "QC_RSD_FILTER",
  "RSD_THRESHOLD", "RUV_K", "SERRF_NUM", "LOESS_SPAN", "LOESS_SAMPLE_SPAN",
  "LOESS_SAMPLE_MIN_OBS", "CORDBAT_REF_BATCH", "NORMALIZATION",
  "SAVE_PRE_CORRECTION_PLOTS", "N_CORES", "sample_type_rules"
)

# Reads a YAML config file and returns it as a named list. Scalar values are
# coerced to character so they flow through get_env() the same way an env
# var string would (a YAML list value, e.g. CORRECTION_METHODS given as a
# bullet list, is collapsed to a comma-separated string first).
# sample_type_rules is left as a nested list for resolve_sample_type_rules().
load_config <- function(path) {
  if (!file.exists(path)) stop("CONFIG_FILE not found: ", path)

  cfg <- tryCatch(
    yaml::read_yaml(path),
    error = function(e) stop("Could not parse CONFIG_FILE (", path, "): ", conditionMessage(e))
  )
  if (is.null(cfg)) cfg <- list()

  unknown <- setdiff(names(cfg), KNOWN_CONFIG_KEYS)
  if (length(unknown) > 0)
    warning("CONFIG_FILE has unrecognized key(s), ignoring: ",
            paste(unknown, collapse = ", "), call. = FALSE)

  for (key in setdiff(names(cfg), "sample_type_rules")) {
    val <- cfg[[key]]
    if (length(val) > 1) val <- paste(val, collapse = ",")
    cfg[[key]] <- as.character(val)
  }

  cfg
}

# Built-in sample-type classification rules, tried in order, first match
# wins. Mirrors the original hardcoded logic in classify_qc()
# (R/msdial_to_notame.R). MSe/MS2 are handled separately (structural, not a
# naming keyword) and are not part of this table.
DEFAULT_SAMPLE_TYPE_RULES <- list(
  list(pattern = "SST\\d",     type = "SST"),
  list(pattern = "ltQC",       type = "ltQC"),
  list(pattern = "sQC",        type = "QC"),
  list(pattern = "MeOH",       type = "Wash"),
  list(pattern = "SolvBlank",  type = "Blank"),
  list(pattern = "blank",      type = "MatrixBlank"),
  list(pattern = "CondPlasma", type = "Cond")
)

# Returns the sample-type rules to use: config$sample_type_rules (parsed
# into list(pattern=, type=) form) if present, else the built-in defaults.
# A config-supplied table fully replaces the defaults rather than
# supplementing them, so a lab's own keyword list can't be silently
# shadowed by an unrelated built-in rule they didn't know about.
#
# Assumes config$sample_type_rules has already passed run_preflight_checks()
# (each entry non-empty, pattern a valid regex) — no validation here, so
# that all config problems can still be reported together up front.
resolve_sample_type_rules <- function(config) {
  raw <- config$sample_type_rules
  if (is.null(raw)) return(DEFAULT_SAMPLE_TYPE_RULES)
  lapply(raw, function(r) list(pattern = as.character(r$pattern), type = as.character(r$type)))
}

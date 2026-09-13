# ─────────────────────────────────────────────────────────────────────────────
# Correction method specification: drift x basis x batch
#
# Each CORRECTION_METHODS entry is a "drift:basis:batch" triple string rather
# than one opaque name per combination -- e.g. "loess:hybrid:combat",
# "huber:qc:feature_median", "none:none:cordbat". This file is the single
# source of truth for what's a legal triple, shared by preflight (fail fast,
# collect all problems before running anything) and the main dispatch loop
# (parses each triple into the drift/basis/batch actually passed to
# run_correction(), see R/correction_methods.R).
# ─────────────────────────────────────────────────────────────────────────────

VALID_DRIFT_METHODS <- c("none", "loess", "huber", "auto", "notame_spline")
VALID_BASES         <- c("none", "qc", "samples", "hybrid")
VALID_BATCH_METHODS <- c("none", "combat", "sva", "limma", "feature_median", "global_median",
                          "ruv_s", "cordbat", "batchcorr", "waveica", "waveica_v1",
                          "pmp_qcrsc", "serrf")

# Batch methods that fully couple drift and batch correction internally and
# so can only be run with drift_method = "none" -- running a separate drift
# step first would either be wasted work or actively wrong depending on the
# method (see R/correction_methods.R's BATCH_METHOD_REGISTRY for why each of
# these is atomic). NOT included: "cordbat" -- CordBat's own correction step
# (run_cordbat()) already accepts an externally drift-corrected, pre-imputed
# SE (that's exactly what loess:samples:cordbat / huber:samples:cordbat do),
# so it behaves like the "raw"-kind batch methods (feature_median etc.), not
# like these five.
ATOMIC_BATCH_METHODS <- c("batchcorr", "waveica", "waveica_v1", "pmp_qcrsc", "serrf")

# Parses one "drift:basis:batch" spec string. Requires exactly 3 non-empty
# colon-separated tokens -- no implicit/omitted basis field, so there's only
# one shape to parse and one error message, regardless of which axis was
# malformed. Returns list(drift=, basis=, batch=) on success, or a character
# vector of problem message(s) on failure (callers distinguish via is.list()).
parse_correction_method_spec <- function(spec) {
  tokens <- strsplit(spec, ":", fixed = TRUE)[[1]]
  if (length(tokens) != 3 || any(nchar(tokens) == 0))
    return(paste0("'", spec, "' is not a valid drift:basis:batch spec (expected exactly 3 ",
                   "non-empty colon-separated fields, e.g. 'loess:hybrid:combat')"))
  list(drift = tokens[1], basis = tokens[2], batch = tokens[3])
}

# Cross-axis legality rules for an already-parsed spec (list(drift=,basis=,batch=)
# from parse_correction_method_spec()). Returns character(0) if legal, else
# problem message(s). Membership in VALID_DRIFT_METHODS/VALID_BASES/
# VALID_BATCH_METHODS is checked first since the cross-axis rules below
# assume valid vocabulary. spec_label (optional) prefixes messages with the
# original spec string, useful when checking several specs from
# CORRECTION_METHODS in one pass and collecting all problems together.
check_method_spec_legality <- function(parsed, spec_label = NULL) {
  problems <- character(0)
  prefix   <- if (is.null(spec_label)) "" else paste0("'", spec_label, "': ")

  if (!parsed$drift %in% VALID_DRIFT_METHODS)
    problems <- c(problems, paste0(prefix, "unknown drift method '", parsed$drift,
                                    "'. Valid: ", paste(VALID_DRIFT_METHODS, collapse = ", ")))
  if (!parsed$basis %in% VALID_BASES)
    problems <- c(problems, paste0(prefix, "unknown basis '", parsed$basis,
                                    "'. Valid: ", paste(VALID_BASES, collapse = ", ")))
  if (!parsed$batch %in% VALID_BATCH_METHODS)
    problems <- c(problems, paste0(prefix, "unknown batch method '", parsed$batch,
                                    "'. Valid: ", paste(VALID_BATCH_METHODS, collapse = ", ")))
  if (length(problems) > 0) return(problems)  # cross-axis rules need valid vocabulary first

  if ((parsed$drift == "none") != (parsed$basis == "none"))
    problems <- c(problems, paste0(prefix, "basis must be 'none' if and only if drift method is ",
                                    "'none' (got drift='", parsed$drift, "', basis='", parsed$basis, "')"))
  if (parsed$batch %in% ATOMIC_BATCH_METHODS && parsed$drift != "none")
    problems <- c(problems, paste0(prefix, "batch method '", parsed$batch, "' handles drift and batch ",
                                    "correction together internally -- drift method must be 'none', got '",
                                    parsed$drift, "'"))
  if (parsed$drift == "notame_spline" && parsed$basis != "qc")
    problems <- c(problems, paste0(prefix, "drift method 'notame_spline' only supports basis='qc' ",
                                    "(no samples-based mode), got '", parsed$basis, "'"))

  problems
}

# Filesystem-safe stand-in for a method spec string -- used only for
# filenames/directory names built from it (colons are legal on Linux/macOS
# but a landmine on Windows and some cloud-sync clients); the raw colon form
# stays the canonical identifier everywhere else (log messages, the `method`
# column value inside CSVs, compare_corrections()'s printed table).
sanitize_method_id <- function(method) gsub(":", "-", method, fixed = TRUE)

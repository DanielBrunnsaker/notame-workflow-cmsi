# ─────────────────────────────────────────────────────────────────────────────
# High-level batch correction method wrappers
#
# All methods follow a two-step imputation strategy:
#   1. LoD/2 (half-minimum) imputation before correction — gives each method a
#      complete matrix to work with without introducing RF-model artifacts
#   2. RF imputation after full correction — imputed values are predicted from
#      batch- and drift-corrected data, giving more biologically meaningful
#      estimates than imputing on raw or partially corrected data
#
# obs_mask (features × samples logical) tracks originally observed values and
# is used to restore NAs before RF imputation. It is returned for downstream
# use in QC metrics.
# ─────────────────────────────────────────────────────────────────────────────

# Half-minimum (LoD/2) imputation, applied per batch.
# Fills each NA with half the minimum observed value for that feature within
# the same batch. Per-batch rather than global because signal levels differ
# between batches — a global minimum from a low-sensitivity batch would give
# inappropriately small placeholders for high-sensitivity batches.
lod2_impute <- function(se) {
  mat     <- assay(se, 1)
  batches <- as.character(colData(se)$Batch)

  for (b in unique(batches)) {
    b_idx <- which(batches == b)
    for (i in seq_len(nrow(mat))) {
      na_idx <- which(is.na(mat[i, b_idx]))
      if (length(na_idx) == 0) next
      finite_min <- min(mat[i, b_idx], na.rm = TRUE)
      if (!is.finite(finite_min)) finite_min <- 1
      mat[i, b_idx[na_idx]] <- 0.5 * finite_min
    }
  }

  assay(se, 1, withDimnames = FALSE) <- mat
  se
}

# Restore originally-missing positions to NA then RF impute.
# se and obs_mask must have matching dimensions.
rf_impute_corrected <- function(se, obs_mask) {
  assay(se, 1, withDimnames = FALSE)[!obs_mask] <- NA
  impute_rf(se, parallelize = "variables")
}

# Clamp non-positive and non-finite values in a SummarizedExperiment assay to
# half the global minimum positive value. pmp's QC-RSC spline correction can
# produce zeros or negatives when the QC spline overshoots; log2 and RF
# imputation both require strictly positive input.
clamp_nonpositive <- function(se, context = "") {
  mat    <- assay(se, 1)
  nonpos <- !is.finite(mat) | mat <= 0
  if (any(nonpos, na.rm = TRUE)) {
    floor_val <- min(mat[is.finite(mat) & mat > 0], na.rm = TRUE) / 2
    mat[nonpos & !is.na(mat)] <- floor_val
    assay(se, 1, withDimnames = FALSE) <- mat
    n <- sum(nonpos & !is.na(mat))
    if (n > 0)
      message("  Note: ", n, " non-positive value(s)",
              if (nchar(context) > 0) paste0(" ", context) else "",
              " clamped to ", signif(floor_val, 3))
  }
  se
}

# ComBat wrapper where mean.only/par.prior are each either a fixed logical or
# the string "auto". When "auto", tries both TRUE and FALSE for that
# parameter (holding any fixed parameter constant) and keeps whichever
# combination gives the best ltQC/Sample D-ratio (eval_ltqc_dratio(), see
# R/qc_metrics.R) on the back-transformed, corrected data. ltQC is never used
# to fit anything upstream of this step, so it's the unbiased check on
# whether a given ComBat configuration under- or over-corrects -- same
# rationale as auto_select_drift_correction()'s pooled D-ratio comparisons.
# ComBat itself is cheap relative to the drift-correction search already run
# for auto_combat, so up to 4 extra ComBat fits here costs little.
#
# `combined` must already be log2-transformed (assay(combined,1) on log2
# scale, as at every call site below). Returns the corrected log2-scale
# matrix -- callers keep doing their own 2^... back-transform afterward.
combat_correct <- function(combined, mean_only = "auto", par_prior = "auto") {
  suppressPackageStartupMessages(library(sva))
  batch   <- as.factor(colData(combined)$Batch)
  log_mat <- assay(combined, 1)

  mo_grid <- if (identical(mean_only, "auto")) c(TRUE, FALSE) else as.logical(mean_only)
  pp_grid <- if (identical(par_prior, "auto")) c(TRUE, FALSE) else as.logical(par_prior)

  if (length(mo_grid) == 1 && length(pp_grid) == 1) {
    message("  ComBat parameters: mean.only=", mo_grid, ", par.prior=", pp_grid)
    return(ComBat(dat = log_mat, batch = batch, mean.only = mo_grid, par.prior = pp_grid))
  }

  grid <- expand.grid(mean_only = mo_grid, par_prior = pp_grid)
  message("  Auto-selecting ComBat parameters (", nrow(grid),
          " combination(s), via ltQC/Sample D-ratio)")

  best_dratio <- Inf
  best_mat    <- NULL
  best_row    <- NULL

  for (i in seq_len(nrow(grid))) {
    mo <- grid$mean_only[i]
    pp <- grid$par_prior[i]
    corrected <- tryCatch(
      ComBat(dat = log_mat, batch = batch, mean.only = mo, par.prior = pp),
      error = function(e) {
        message("    mean.only=", mo, ", par.prior=", pp, ": ComBat failed (",
                conditionMessage(e), ")")
        NULL
      }
    )
    if (is.null(corrected)) next

    se_trial <- combined
    assay(se_trial, 1, withDimnames = FALSE) <- 2^corrected
    dratio <- eval_ltqc_dratio(se_trial)
    message("    mean.only=", mo, ", par.prior=", pp, ": ltQC/Sample D-ratio = ",
            if (is.na(dratio)) "NA" else round(dratio, 4))

    if (!is.na(dratio) && dratio < best_dratio) {
      best_dratio <- dratio
      best_mat    <- corrected
      best_row    <- grid[i, ]
    }
  }

  if (is.null(best_mat)) {
    message("  No candidate's ltQC/Sample D-ratio could be computed (insufficient ltQC) — ",
            "defaulting to mean.only=TRUE, par.prior=TRUE")
    return(ComBat(dat = log_mat, batch = batch, mean.only = TRUE, par.prior = TRUE))
  }

  message("  Selected ComBat parameters: mean.only=", best_row$mean_only,
          ", par.prior=", best_row$par_prior,
          "  (ltQC/Sample D-ratio = ", round(best_dratio, 4), ")")
  best_mat
}

# SVA (Surrogate Variable Analysis, `sva` package -- already a dependency via
# ComBat) as an alternative between-batch correction. Unlike ComBat, it isn't
# limited to a per-batch mean/variance shift model: it estimates n_sv latent
# "surrogate variables" representing systematic (non-random) structure in the
# data not already explained by known Batch, then regresses out Batch AND
# those surrogate variables together via limma::removeBatchEffect() -- so it
# can pick up additional technical structure Batch alone doesn't fully
# capture, at the cost of being less interpretable than ComBat's simple
# per-batch shift.
#
# mod == mod0 == ~1 (intercept only) for the estimation step -- NOT ~Batch.
# An earlier version of this function used mod = mod0 = model.matrix(~Batch),
# on the theory that this would make the estimated surrogate variables
# represent only structure *beyond* Batch. In practice that's a degenerate
# input to sva(): mod and mod0 identical but non-trivial makes every
# per-feature F-test sva() runs internally compare a model against itself
# (0 residual degrees of freedom for the "variable of interest" term),
# producing NaN p-values for every feature -- which is what was actually
# causing the "'x' contains missing values" failure during sva()'s iterative
# reweighting (not, as first suspected, non-finite values in the data
# matrix itself; see clamp_nonpositive() calls elsewhere in this file, which
# fixed a real but different latent-NaN risk and should stay). mod = mod0 = ~1
# is sva()'s own documented recipe for "no known variable of interest, no
# known covariates to protect during estimation" -- it may end up
# re-discovering batch-like structure as one of its own surrogate variables,
# which is fine: known Batch is still regressed out explicitly and
# separately at the removeBatchEffect() step below regardless of what sva()
# found, so nothing depends on mod excluding it.
#
# n_sv = NULL auto-estimates the surrogate variable count via
# sva::num.sv(..., method = "be") (Buja-Eyuboglu permutation test, the
# package's own recommended default); pass an integer to fix it instead --
# same NULL-means-auto convention as WAVEICA_K.
#
# `combined` must already be log2-transformed, as at every ComBat call site
# above. Returns the corrected log2-scale matrix.
sva_correct <- function(combined, n_sv = NULL) {
  suppressPackageStartupMessages(library(sva))
  suppressPackageStartupMessages(library(limma))

  log_mat <- assay(combined, 1)
  batch   <- as.factor(colData(combined)$Batch)
  mod     <- model.matrix(~1, data = as.data.frame(colData(combined)))
  mod0    <- mod

  n_sv_eff <- n_sv
  if (is.null(n_sv_eff)) {
    n_sv_eff <- tryCatch(
      suppressWarnings(num.sv(log_mat, mod, method = "be")),
      error = function(e) {
        message("  num.sv() failed to estimate surrogate variable count (", conditionMessage(e),
                ") — defaulting to 0 (Batch-only correction)")
        0L
      }
    )
    message("  Auto-estimated surrogate variable count: ", n_sv_eff)
  } else {
    message("  Using fixed surrogate variable count: ", n_sv_eff)
  }

  if (n_sv_eff <= 0) {
    message("  No surrogate variables to estimate — correcting for known Batch only")
    return(removeBatchEffect(log_mat, batch = batch))
  }

  svobj <- tryCatch(
    suppressWarnings(sva(log_mat, mod, mod0, n.sv = n_sv_eff)),
    error = function(e) {
      message("  sva() failed (", conditionMessage(e), ") — correcting for known Batch only")
      NULL
    }
  )

  if (is.null(svobj) || is.null(svobj$sv) || NCOL(svobj$sv) == 0) {
    return(removeBatchEffect(log_mat, batch = batch))
  }

  message("  Estimated ", NCOL(svobj$sv), " surrogate variable(s) — regressing out Batch + surrogate variables")
  removeBatchEffect(log_mat, batch = batch, covariates = svobj$sv)
}

# `combined` must already be log2-transformed; returns the corrected log2-scale
# matrix. Extracted from the inline removeBatchEffect() call the old
# correct_loess_limma()/correct_loess_samples_limma() used, for symmetry with
# combat_correct()/sva_correct() so limma fits the same "log2"-kind batch-method
# contract in BATCH_METHOD_REGISTRY below.
limma_correct <- function(combined) {
  suppressPackageStartupMessages(library(limma))
  removeBatchEffect(x = assay(combined, 1), batch = as.factor(colData(combined)$Batch))
}

# Registry of batch-correction methods (see R/method_spec.R for the full
# vocabulary). Each entry's `kind` tells run_correction() how to treat it:
#   "none"     -- no between-batch step at all.
#   "log2"     -- clamp -> log2 -> fn(combined, params) -> 2^ back-transform.
#                 fn must return a corrected log2-scale matrix. Requires a
#                 complete (LoD/2-imputed) input matrix, which run_correction()
#                 provides. Used by combat/sva/limma.
#   "complete" -- fn(combined, params) receives an already LoD/2-imputed SE and
#                 returns a corrected SE; fn manages its own scale handling
#                 internally (run_cordbat() already does its own clamp/log2/2^
#                 back-transform; ruvs_qc() works directly on raw scale).
#   "sparse"   -- fn(combined, params) receives the merged, NOT LoD/2-imputed,
#                 SE directly and returns a corrected SE; fn does its own
#                 NA-tolerant filtering. feature_median/global_median were
#                 built and validated this way -- imputing before them would
#                 let filled-in placeholder values influence the per-batch/
#                 per-feature median, which is not how they were designed.
#   "atomic"   -- fn(data, params) receives the ORIGINAL, unmodified input and
#                 returns the full list(pre=, post=, obs_mask=) itself; drift
#                 correction and LoD/2 imputation are skipped entirely by
#                 run_correction() for these -- the method does its own
#                 complete pipeline (own imputation, own obs_mask, own final
#                 RF-impute). Preflight (R/method_spec.R's
#                 ATOMIC_BATCH_METHODS) forces drift_method="none" whenever
#                 one of these is selected, since running a separate drift
#                 step first would either be wasted work or actively wrong.
BATCH_METHOD_REGISTRY <- list(
  none = list(kind = "none"),
  combat = list(kind = "log2", fn = function(se, p)
    combat_correct(se, mean_only = p$combat_mean_only, par_prior = p$combat_par_prior)),
  sva = list(kind = "log2", fn = function(se, p)
    sva_correct(se, n_sv = p$sva_n_sv)),
  limma = list(kind = "log2", fn = function(se, p) limma_correct(se)),
  feature_median = list(kind = "sparse", fn = function(se, p) batch_feature_median_correct(se)),
  global_median  = list(kind = "sparse", fn = function(se, p) batch_global_median_correct(se)),
  ruv_s = list(kind = "complete", fn = function(se, p)
    ruvs_qc(se, replicates = list(which(colData(se)$QC == "QC")), k = p$ruv_k)),
  cordbat = list(kind = "complete", fn = function(se, p) run_cordbat(se, p$cordbat_ref_batch)),
  batchcorr  = list(kind = "atomic", fn = function(data, p) correct_batchcorr(data)),
  waveica    = list(kind = "atomic", fn = function(data, p)
    correct_waveica(data, alpha = p$waveica_alpha, cutoff = p$waveica_cutoff,
                     K = p$waveica_k, wf = p$waveica_wf, eval_group = p$waveica_eval_group)),
  waveica_v1 = list(kind = "atomic", fn = function(data, p)
    correct_waveica_v1(data, wf = p$waveica_v1_wf, K = p$waveica_v1_k, t = p$waveica_v1_t,
                        t2 = p$waveica_v1_t2, alpha = p$waveica_v1_alpha)),
  pmp_qcrsc = list(kind = "atomic", fn = function(data, p) correct_pmp_qcrsc(data)),
  serrf     = list(kind = "atomic", fn = function(data, p) correct_serrf(data, num = p$serrf_num_eff))
)

# Generic correction-method runner: dispatches the drift step via
# resolve_drift_fn()/auto_select_drift_correction() (R/drift_correction.R) and
# the between-batch step via BATCH_METHOD_REGISTRY above, replacing what used
# to be one bespoke correct_*() wrapper function per named drift+batch
# combination. `params` is the single shared list of resolved config values
# built once in notame-workflow.r.
run_correction <- function(data, drift_method = "none", basis = "none",
                            batch_method = "none", params = list()) {
  entry <- BATCH_METHOD_REGISTRY[[batch_method]]
  if (is.null(entry)) stop("Unknown batch_method: ", batch_method)

  if (entry$kind == "atomic") return(entry$fn(data, params))

  # --- 1. Drift step ---
  if (drift_method == "none") {
    combined <- data
  } else if (drift_method == "auto") {
    combined <- merge_notame_sets(
      auto_select_drift_correction(data, basis = basis,
        loess_spans = params$auto_loess_spans, huber_ks = params$auto_huber_ks,
        sample_loess_spans = params$auto_sample_loess_spans,
        sample_huber_ks = params$auto_sample_huber_ks,
        min_qc_per_batch = params$auto_min_qc_per_batch,
        min_ltqc_validate = params$auto_min_ltqc_validate,
        min_cv_obs = params$auto_min_cv_obs),
      merge = "samples"
    )
  } else {
    drift_fn <- resolve_drift_fn(drift_method, basis)
    combined <- merge_notame_sets(
      lapply(split_by_batch(data), function(se_b) drift_fn(se_b, params)),
      merge = "samples"
    )
  }

  # Capture obs_mask after merge so column order matches combined
  obs_mask  <- !is.na(assay(combined, 1))
  n_batches <- length(unique(colData(combined)$Batch))

  if (entry$kind == "none") {
    message("==> Imputation (RF on corrected data)")
    combined <- rf_impute_corrected(combined, obs_mask)
    return(list(pre = combined, post = combined, obs_mask = obs_mask))
  }

  # LoD/2 fill before a batch method that requires a complete matrix
  # ("log2"/"complete" kinds); "sparse" kinds handle their own NA filtering
  # and must NOT be pre-imputed (see BATCH_METHOD_REGISTRY's documentation).
  if (entry$kind %in% c("log2", "complete")) combined <- lod2_impute(combined)
  pre <- combined

  if (n_batches < 2) {
    message("==> Batch correction skipped (only one batch detected)")
  } else if (entry$kind == "log2") {
    combined <- clamp_nonpositive(combined, "before log2")
    message("==> Log2 transformation")
    assay(combined, 1, withDimnames = FALSE) <- log2(assay(combined, 1))
    message("==> Between-batch correction (", batch_method, ")")
    assay(combined, 1, withDimnames = FALSE) <- entry$fn(combined, params)
    message("==> Back-transforming to raw scale")
    assay(combined, 1, withDimnames = FALSE) <- 2^assay(combined, 1)
  } else {
    message("==> Between-batch correction (", batch_method, ")")
    combined <- entry$fn(combined, params)
  }

  message("==> Imputation (RF on corrected data)")
  combined <- rf_impute_corrected(combined, obs_mask)

  list(pre = pre, post = combined, obs_mask = obs_mask)
}

correct_none <- function(data) {
  message("==> No correction (imputation only)")
  obs_mask <- !is.na(assay(data, 1))
  combined <- impute_rf(data, parallelize = "variables")
  list(pre = combined, post = combined, obs_mask = obs_mask)
}

# Per-feature batch median ratio correction.
# Scales each batch so its biological sample median per feature matches the
# grand median. More flexible than global scaling but noisier for sparse features.
batch_feature_median_correct <- function(se) {
  mat      <- assay(se, 1)
  cd       <- as.data.frame(colData(se))
  samp_idx <- which(cd$QC == "Sample")
  batches  <- unique(as.character(cd$Batch))

  grand_med <- apply(mat[, samp_idx, drop = FALSE], 1, function(x) {
    ok <- is.finite(x) & x > 0
    if (sum(ok) < 2) NA_real_ else median(x[ok])
  })

  for (b in batches) {
    b_idx      <- which(as.character(cd$Batch) == b)
    b_samp_idx <- intersect(b_idx, samp_idx)

    if (length(b_samp_idx) < 2) {
      message("  Batch ", b, ": skipped (insufficient samples)")
      next
    }

    batch_med    <- apply(mat[, b_samp_idx, drop = FALSE], 1, function(x) {
      ok <- is.finite(x) & x > 0
      if (sum(ok) < 2) NA_real_ else median(x[ok])
    })
    scale_factor <- grand_med / batch_med
    scale_factor[!is.finite(scale_factor) | scale_factor <= 0] <- 1

    mat[, b_idx] <- mat[, b_idx] * scale_factor
    message("  Batch ", b, ": scaled ", sum(is.finite(scale_factor) & scale_factor != 1),
            "/", nrow(mat), " features")
  }

  assay(se, 1, withDimnames = FALSE) <- mat
  se
}

# Global batch median ratio correction.
# Computes one scaling factor per batch from the median of all biological sample
# intensities, then applies it uniformly to all features. Assumes a constant
# multiplicative offset across all features within a batch.
batch_global_median_correct <- function(se) {
  mat      <- assay(se, 1)
  cd       <- as.data.frame(colData(se))
  samp_idx <- which(cd$QC == "Sample")
  batches  <- unique(as.character(cd$Batch))

  all_samp_vals <- mat[, samp_idx, drop = FALSE]
  grand_med     <- median(all_samp_vals[is.finite(all_samp_vals) & all_samp_vals > 0],
                          na.rm = TRUE)

  if (!is.finite(grand_med) || grand_med <= 0) {
    message("  WARNING: could not compute grand median — skipping global batch correction")
    return(se)
  }

  for (b in batches) {
    b_idx      <- which(as.character(cd$Batch) == b)
    b_samp_idx <- intersect(b_idx, samp_idx)

    if (length(b_samp_idx) < 2) {
      message("  Batch ", b, ": skipped (insufficient samples)")
      next
    }

    b_vals    <- mat[, b_samp_idx, drop = FALSE]
    batch_med <- median(b_vals[is.finite(b_vals) & b_vals > 0], na.rm = TRUE)

    if (!is.finite(batch_med) || batch_med <= 0) {
      message("  Batch ", b, ": skipped (could not compute batch median)")
      next
    }

    scale_factor <- grand_med / batch_med
    mat[, b_idx] <- mat[, b_idx] * scale_factor
    message("  Batch ", b, ": scale factor = ", round(scale_factor, 4))
  }

  assay(se, 1, withDimnames = FALSE) <- mat
  se
}

# Print per-batch counts of non-positive (and NA) values for diagnostic purposes.
# Call before and after QCRSC to identify which batches pmp is distorting.
diag_nonpositive <- function(se, label = "") {
  mat     <- assay(se, 1)
  batches <- as.character(colData(se)$Batch)
  counts  <- sapply(unique(batches), function(b) {
    m <- mat[, batches == b, drop = FALSE]
    sum(!is.finite(m) | m <= 0, na.rm = TRUE)
  })
  total <- sum(counts)
  if (nchar(label) > 0) message("  [", label, "] non-positive values per batch:")
  for (b in names(counts))
    message("    Batch ", b, ": ", counts[b])
  message("    Total: ", total)
  invisible(counts)
}

correct_pmp_qcrsc <- function(data) {
  suppressPackageStartupMessages(library(pmp))

  obs_mask <- !is.na(assay(data, 1))

  # Identify batches pmp will not have QC anchors for (< 4 QC samples).
  # minQC=4 excludes their QCs from spline fitting but pmp still extrapolates
  # the global spline into those batches, producing wildly incorrect values.
  # Save their original values and restore after QCRSC.
  cd         <- as.data.frame(colData(data))
  qc_counts  <- tapply(cd$QC == "QC", as.character(cd$Batch), sum)
  skip_batches <- names(qc_counts[qc_counts < 4])
  orig_mat   <- assay(data, 1)

  message("==> Drift correction + batch correction (pmp QC-RSC)")
  # QCRSC handles NAs natively — no LoD/2 before this step
  # ltQC remapped to "Sample" so pmp does not try to use it as a QC reference
  classes_for_pmp <- ifelse(colData(data)$QC == "QC", "QC", "Sample")

  diag_nonpositive(data, "before QC-RSC")
  combined <- QCRSC(
    df      = data,
    order   = colData(data)$Injection_order,
    batch   = colData(data)$Batch,
    classes = classes_for_pmp,
    spar    = 0,
    minQC   = 4
  )
  diag_nonpositive(combined, "after QC-RSC")

  if (length(skip_batches) > 0) {
    skip_idx <- which(as.character(colData(combined)$Batch) %in% skip_batches)
    assay(combined, 1, withDimnames = FALSE)[, skip_idx] <- orig_mat[, skip_idx]
    message("  Restored pre-correction values for batch(es) with <4 QCs: ",
            paste(skip_batches, collapse = ", "))
    diag_nonpositive(combined, "after restore")
  }

  # Clamp any non-positive values introduced by spline overshoot in corrected
  # batches (small in number but would break RF imputation's internal log step)
  combined <- clamp_nonpositive(combined, "after QC-RSC")

  message("==> Imputation (RF on corrected data)")
  combined <- rf_impute_corrected(combined, obs_mask)

  list(pre = combined, post = combined, obs_mask = obs_mask)
}

correct_batchcorr <- function(data,
                              G          = seq(5, 35, by = 10),
                              modelNames = c("VVV", "VVE", "VEV", "VEE", "VEI", "VVI", "VII"),
                              qualRatio  = 0.4) {
  suppressPackageStartupMessages(library(batchCorr))

  obs_mask <- !is.na(assay(data, 1))

  message("==> Pre-imputation (LoD/2, required by batchCorr)")
  data <- lod2_impute(data)

  # batchCorr expects samples × features matrix
  mat  <- t(assay(data, "abundances"))
  meta <- as.data.frame(colData(data))

  message("==> Within-batch drift correction (batchCorr cluster-based spline)")
  batches        <- unique(as.character(meta$Batch))
  batch_corrObjs <- list()

  for (b in batches) {
    message("  Batch ", b)
    idx   <- which(as.character(meta$Batch) == b)
    bmat  <- mat[idx, , drop = FALSE]
    bmeta <- meta[idx, ]
    ord   <- order(bmeta$Injection_order)
    bmat  <- bmat[ord, , drop = FALSE]
    bmeta <- bmeta[ord, ]
    sgrp  <- ifelse(bmeta$QC == "QC", "QC", "Sample")

    bc <- tryCatch({
      correctDrift(
        peakTable    = bmat,
        injections   = bmeta$Injection_order,
        sampleGroups = sgrp,
        QCID         = "QC",
        G            = G,
        modelNames   = modelNames,
        CVlimit      = Inf,
        report       = FALSE
      )
    }, error = function(e) {
      message("    WARNING: correctDrift failed for batch ", b, ": ", conditionMessage(e))
      NULL
    })

    if (!is.null(bc)) batch_corrObjs[[b]] <- bc
  }

  if (length(batch_corrObjs) == 0) stop("correctDrift failed for all batches")

  if (length(batch_corrObjs) == 1) {
    message("==> Single batch detected — skipping mergeBatches, using corrected batch directly")
    b      <- names(batch_corrObjs)[1]
    bc     <- batch_corrObjs[[b]]

    # correctDrift() has no $peakTable field — its actual outputs (per its own
    # console messages) are $TestFeatsCorr (all corrected features) and
    # $TestFeatsFinal (after its internal QC-CV filter; identical to
    # TestFeatsCorr here since CVlimit = Inf disables that filter). Using the
    # nonexistent $peakTable silently gave NULL -> 0 kept features -> a
    # downstream SummarizedExperiment row-count mismatch crash.
    peak_corr <- bc$TestFeatsFinal
    if (is.null(peak_corr)) peak_corr <- bc$TestFeatsCorr
    if (is.null(peak_corr))
      stop("correctDrift() result for batch ", b, " has neither TestFeatsFinal nor ",
           "TestFeatsCorr. Available fields: ", paste(names(bc), collapse = ", "))

    b_idx        <- which(as.character(meta$Batch) == b)
    b_idx        <- b_idx[order(meta[b_idx, "Injection_order"])]
    peakTableOrg <- mat[b_idx, , drop = FALSE]

    # TestFeatsCorr/TestFeatsFinal's row coverage isn't confirmed against
    # package docs offline — flag loudly rather than silently dropping
    # samples (e.g. if it only covers "Test"/non-QC rows, not QC).
    missing_samples <- setdiff(rownames(peakTableOrg), rownames(peak_corr))
    if (length(missing_samples) > 0)
      message("  WARNING: ", length(missing_samples), " sample(s) from batch ", b,
              " missing from correctDrift()'s corrected table and will be dropped: ",
              paste(missing_samples, collapse = ", "))

    merged <- list(
      peakTableCorr = peak_corr,
      peakTableOrg  = peakTableOrg
    )
  } else {
    message("==> Merging batches")
    merged <- mergeBatches(batch_corrObjs, qualRatio = qualRatio)
  }

  kept_features <- colnames(merged$peakTableCorr)
  kept_samples  <- rownames(merged$peakTableCorr)

  pre <- data[kept_features, kept_samples]
  assay(pre, "abundances", withDimnames = FALSE) <- t(merged$peakTableOrg[kept_samples, ])

  combined <- pre
  n_batches <- length(unique(as.character(colData(pre)$Batch)))
  if (n_batches < 2) {
    message("==> Between-batch normalization skipped (single batch)")
    assay(combined, "abundances", withDimnames = FALSE) <- t(merged$peakTableCorr[kept_samples, ])
  } else {
    message("==> Between-batch normalization (batchCorr normalizeBatches)")
    sgrp_merged <- ifelse(colData(pre)$QC == "QC", "QC", "Sample")
    norm_result <- normalizeBatches(
      peakTableCorr = merged$peakTableCorr[kept_samples, ],
      batches       = as.character(colData(pre)$Batch),
      sampleGroup   = sgrp_merged,
      refGroup      = "QC",
      population    = "all",
      CVlimit       = Inf
    )
    assay(combined, "abundances", withDimnames = FALSE) <- t(norm_result$peakTable)
  }

  obs_mask <- obs_mask[kept_features, kept_samples, drop = FALSE]

  message("==> Imputation (RF on corrected data)")
  combined <- rf_impute_corrected(combined, obs_mask)

  list(pre = pre, post = combined, obs_mask = obs_mask)
}

# Selects the reference batch for CordBat.
# Prefers batches with QC samples and ranks by median robust QC RSD (MAD/median).
# Batches without QC samples are excluded from candidacy — their quality cannot
# be verified and they received no LOESS drift correction.
# Falls back to biological sample RSD only if no batch has QC samples.
select_ref_batch_cordbat <- function(se) {
  mat     <- assay(se, 1)
  cd      <- as.data.frame(colData(se))
  batches <- unique(as.character(cd$Batch))

  qc_batches <- intersect(batches, unique(as.character(cd$Batch[cd$QC == "QC"])))

  if (length(qc_batches) > 0) {
    rsd_med <- sapply(qc_batches, function(b) {
      idx <- which(as.character(cd$Batch) == b & cd$QC == "QC")
      if (length(idx) < 2) return(Inf)
      median(apply(mat[, idx, drop = FALSE], 1, function(x) {
        x <- x[is.finite(x) & x > 0]
        if (length(x) < 2) NA_real_ else mad(x) / median(x)
      }), na.rm = TRUE)
    })
    names(rsd_med) <- qc_batches
    message("  Batch QC RSD medians: ",
            paste(qc_batches, "=", round(rsd_med, 3), collapse = ", "))
    if (length(batches) > length(qc_batches))
      message("  Excluded from candidacy (no QC samples): ",
              paste(setdiff(batches, qc_batches), collapse = ", "))
  } else {
    message("  WARNING: no batches have QC samples — falling back to biological sample RSD")
    samp_idx <- which(cd$QC == "Sample")
    rsd_med <- sapply(batches, function(b) {
      idx <- intersect(which(as.character(cd$Batch) == b), samp_idx)
      if (length(idx) < 2) return(Inf)
      median(apply(mat[, idx, drop = FALSE], 1, function(x) {
        x <- x[is.finite(x) & x > 0]
        if (length(x) < 2) NA_real_ else sd(x) / mean(x)
      }), na.rm = TRUE)
    })
    names(rsd_med) <- batches
    qc_batches <- batches
  }

  ref <- names(rsd_med)[which.min(rsd_med)]
  message("  Auto-selected reference batch: ", ref)
  ref
}

# Internal: runs CordBat on a LoD/2-imputed SE.
# Transforms to log2 before calling CordBat, back-transforms after.
# ref_batch: batch ID string, or NULL for auto-selection.
run_cordbat <- function(combined, ref_batch) {
  if (!exists("CordBat", mode = "function")) {
    suppressPackageStartupMessages(library(igraph))
    source("R/Funcs_CordBat_algorithm.R")
  }

  n_batches <- length(unique(colData(combined)$Batch))
  if (n_batches < 2) {
    message("==> Batch correction skipped (only one batch detected)")
    return(combined)
  }

  if (is.null(ref_batch)) ref_batch <- select_ref_batch_cordbat(combined)

  combined <- clamp_nonpositive(combined, "before log2")
  message("==> Log2 transformation")
  mat_log <- log2(assay(combined, 1))

  X     <- t(mat_log)
  batch <- as.character(colData(combined)$Batch)
  # QC samples flagged so CordBat corrects them via the same coefficients
  # without including them in GGM estimation. ltQC treated as samples.
  group <- ifelse(colData(combined)$QC == "QC", "QC", "Sample")

  # StARS subsamples 70% of the reference batch repeatedly to select the GGM
  # regularisation parameter. Features that are near-constant within the
  # reference batch will be constant in many of those subsamples, crashing
  # scale(). Exclude them here and restore after correction — they will be
  # handled by RF imputation downstream anyway.
  ref_idx  <- which(batch == ref_batch)
  ref_var  <- apply(X[ref_idx, , drop = FALSE], 2, var, na.rm = TRUE)
  keep_idx <- which(is.finite(ref_var) & ref_var > .Machine$double.eps)
  n_excl   <- ncol(X) - length(keep_idx)
  if (n_excl > 0)
    message("  Excluding ", n_excl, " near-constant features in reference batch",
            " from CordBat (insufficient variance for GGM subsampling)")

  X_input <- X[, keep_idx, drop = FALSE]

  message("==> Between-batch correction (CordBat, ref = ", ref_batch, ")")
  result <- tryCatch(
    CordBat(X = X_input, batch = batch, group = group, grouping = FALSE,
            ref.batch = ref_batch, eps = 1e-5, print.detail = FALSE),
    error = function(e) stop("CordBat failed: ", conditionMessage(e))
  )

  X_cor_sub <- result$X.cor.withQC
  if (is.null(X_cor_sub)) X_cor_sub <- result$X.cor

  # Reinsert corrected values; excluded features keep their LoD/2 log2 values
  X_cor <- X
  X_cor[, keep_idx] <- X_cor_sub

  message("==> Back-transforming to raw scale")
  assay(combined, 1, withDimnames = FALSE) <- t(2^X_cor)
  combined
}


# Runs one WaveICA_2.0() call and returns the corrected raw-scale matrix
# (samples transposed back to features x samples, to match this pipeline's
# assay convention). NA in k means "auto" (2 x n_batches, WaveICA2.0's own
# implicit default in this pipeline) -- resolved here, not by the caller.
run_waveica_once <- function(data, alpha, cutoff, k, wf) {
  k_eff <- if (is.na(k)) length(unique(colData(data)$Batch)) * 2 else k
  result <- WaveICA_2.0(
    data            = t(assay(data, 1)),
    wf              = wf,
    Injection_Order = as.numeric(colData(data)$Injection_order),
    alpha           = alpha,
    Cutoff          = cutoff,
    K               = k_eff
  )
  list(mat = t(result$data_wave), k_eff = k_eff)
}

# Auto-selects WaveICA2.0's (alpha, Cutoff, K) from the cross-product of
# alpha_grid/cutoff_grid/k_grid -- each a vector from a comma-separated env
# var (WAVEICA_ALPHA/WAVEICA_CUTOFF/WAVEICA_K); a single value in every grid
# keeps today's fixed behaviour (one WaveICA_2.0() call, no search); more
# than one candidate overall triggers this search. NA within k_grid means
# "auto" (2 x n_batches).
#
# WaveICA2.0 never uses QC or ltQC to fit itself -- it corrects using only
# injection order -- so evaluating it against either afterward is a genuine
# held-out check regardless of which one is chosen (unlike a QC-anchored
# drift method, where evaluating against QC would be circular). eval_group
# ("ltQC" or "QC") picks which one; ltQC is the pipeline-wide default, QC is
# there for when it has more samples to draw on.
#
# Three metrics are computed per candidate, but only one drives selection:
#   - eval_group/Sample D-ratio -- primary, selects the winner. Consistent
#     with every other auto-selection in this pipeline (combat_correct(),
#     sva_correct(), auto_select_drift_correction()).
#   - eval_dist_ratio() (PCA-space distance ratio) -- informational only.
#     D-ratio is a per-feature metric; WaveICA2.0 is an ICA-based method
#     whose entire mechanism operates jointly across features, so a
#     per-feature-only metric could miss damage to that joint structure.
#     This is the multivariate cross-check for that specific blind spot.
#   - eval_qc_homogeneity()'s PERMANOVA R²(Batch) -- informational only.
#     D-ratio and dist_ratio both measure technical-noise-vs-signal; neither
#     directly asks whether batch-driven clustering actually went down,
#     which is the literal purpose of a batch-correction method. A separate,
#     genuinely different axis worth seeing alongside the other two.
# All three are printed for every candidate so a disagreement is visible
# (same reasoning as auto_select_drift_correction() printing ltQC/Sample
# D-ratio alongside pooled LOO-CV score for its QC-tier candidates).
#
# Candidates run in parallel via foreach, reusing this pipeline's existing
# registerDoParallel() setup (notame-workflow.r). WaveICA_2.0() itself calls
# parallel::mclapply() internally (see its own R/WaveICA_2.0.R upstream,
# github.com/dengkuistat/WaveICA_2.0), defaulting to 2 cores independently of
# anything this pipeline configures -- so each outer worker pins
# options(mc.cores = 1) for the duration of its own WaveICA_2.0() call, to
# keep this pipeline's N_CORES/foreach grid the only source of parallelism.
# Without that, n_outer_workers x 2 processes would compete for the same
# cores instead of actually parallelizing.
select_waveica_params <- function(data, alpha_grid, cutoff_grid, k_grid, wf, eval_group = "ltQC") {
  suppressPackageStartupMessages(library(foreach))

  grid <- expand.grid(alpha = alpha_grid, cutoff = cutoff_grid, k = k_grid)

  if (nrow(grid) == 1) {
    once <- run_waveica_once(data, grid$alpha[1], grid$cutoff[1], grid$k[1], wf)
    message("  WaveICA2.0 parameters: alpha=", grid$alpha[1], ", Cutoff=", grid$cutoff[1],
            ", K=", once$k_eff, ", wf=", wf)
    return(once$mat)
  }

  message("  Auto-selecting WaveICA2.0 parameters (", nrow(grid),
          " combination(s), evaluated against ", eval_group, "/Sample)")

  results <- foreach(i = seq_len(nrow(grid))) %dopar% {
    options(mc.cores = 1)
    once <- tryCatch(run_waveica_once(data, grid$alpha[i], grid$cutoff[i], grid$k[i], wf),
                      error = function(e) NULL)
    if (is.null(once)) {
      list(mat = NULL, k_eff = NA_real_, dratio = NA_real_, dist_ratio = NA_real_, permanova_r2 = NA_real_)
    } else {
      se_trial <- data
      assay(se_trial, 1, withDimnames = FALSE) <- once$mat
      list(
        mat          = once$mat,
        k_eff        = once$k_eff,
        dratio       = eval_ltqc_dratio(se_trial, reference_group = eval_group),
        dist_ratio   = eval_dist_ratio(se_trial, group1 = eval_group, group2 = "Sample"),
        permanova_r2 = eval_qc_homogeneity(se_trial, group = eval_group)$permanova_r2
      )
    }
  }

  dratios <- vapply(results, `[[`, numeric(1), "dratio")
  for (i in seq_len(nrow(grid))) {
    r <- results[[i]]
    message("    alpha=", grid$alpha[i], ", Cutoff=", grid$cutoff[i], ", K=", r$k_eff,
            ": ", eval_group, "/Sample D-ratio = ", if (is.na(r$dratio)) "NA" else round(r$dratio, 4),
            ", dist_ratio = ", if (is.na(r$dist_ratio)) "NA" else round(r$dist_ratio, 4),
            ", ", eval_group, " PERMANOVA R2(Batch) = ",
            if (is.na(r$permanova_r2)) "NA" else round(r$permanova_r2, 4))
  }

  if (all(is.na(dratios))) {
    message("  No candidate's ", eval_group, "/Sample D-ratio could be computed -- ",
            "defaulting to alpha=", alpha_grid[1], ", Cutoff=", cutoff_grid[1])
    once <- run_waveica_once(data, alpha_grid[1], cutoff_grid[1], k_grid[1], wf)
    return(once$mat)
  }

  winner <- which.min(dratios)
  message("  Selected: alpha=", grid$alpha[winner], ", Cutoff=", grid$cutoff[winner],
          ", K=", results[[winner]]$k_eff,
          "  (", eval_group, "/Sample D-ratio = ", round(dratios[winner], 4), ")")
  results[[winner]]$mat
}

correct_waveica <- function(data, alpha = 0.05, cutoff = 0.10, K = NA_real_, wf = "haar",
                             eval_group = "ltQC") {
  suppressPackageStartupMessages(library(WaveICA2.0))

  obs_mask <- !is.na(assay(data, 1))
  data     <- lod2_impute(data)

  message("==> WaveICA2.0 correction")
  assay(data, 1, withDimnames = FALSE) <- select_waveica_params(
    data, alpha_grid = alpha, cutoff_grid = cutoff, k_grid = K, wf = wf, eval_group = eval_group
  )

  message("==> Imputation (RF on corrected data)")
  combined <- rf_impute_corrected(data, obs_mask)

  list(pre = combined, post = combined, obs_mask = obs_mask)
}


# Original WaveICA (Deng et al.), as an alternative to WaveICA2.0 (waveica).
# Unlike WaveICA2.0, this version takes batch labels directly rather than
# using injection order as a proxy for batch structure -- preferable when
# batch labels are known and reliable, since it tests ICA components against
# the real grouping variable instead of a continuous stand-in for it.
#
# Defaults below (wf, K, t, t2, alpha) are the package's own defaults, not
# tuned for this pipeline -- there was no way to verify or tune this wrapper
# offline (the package isn't installed anywhere this session had access to,
# unlike WaveICA2.0/CordBat), so treat this method's first real run as
# genuinely first-run, not just "should work like the others already tested."
#
# `group` (optional biological comparison group, protected from removal via
# `t2`) is intentionally not exposed here -- this pipeline's notame format
# has no biological-group column to supply it from. Left as the package's
# own default (NULL); add a parameter for it later if that data becomes
# available.
#
# Note: this function's `alpha` (0-1, trade-off between sample-wise and
# variable-wise independence in the ICA step) turns out to be the same KIND
# of parameter as WaveICA2.0's own `alpha` (see select_waveica_params()'s
# documentation -- both packages' alpha is an ICA spatial/temporal
# independence trade-off, not a significance/flagging threshold; an earlier
# version of this comment claimed they were unrelated, which was wrong, based
# on a mis-reading of WaveICA2.0's own parameter that has since been
# corrected against its actual source). Still kept as WAVEICA_V1_ALPHA
# (distinct from WAVEICA_ALPHA) in notame-workflow.r since they're separate
# packages with separately-tuned defaults, not because the concept differs.
correct_waveica_v1 <- function(data, wf = "haar", K = 20, t = 0.05, t2 = 0.05, alpha = 0) {
  suppressPackageStartupMessages(library(WaveICA))

  obs_mask <- !is.na(assay(data, 1))
  data     <- lod2_impute(data)

  message("==> WaveICA (v1) correction (wf=", wf, ", K=", K, ", t=", t,
          ", t2=", t2, ", alpha=", alpha, ")")
  result <- WaveICA(
    data  = t(assay(data, 1)),
    wf    = wf,
    batch = as.character(colData(data)$Batch),
    group = NULL,
    K     = K,
    t     = t,
    t2    = t2,
    alpha = alpha
  )
  assay(data, 1, withDimnames = FALSE) <- t(result$data_wave)

  message("==> Imputation (RF on corrected data)")
  combined <- rf_impute_corrected(data, obs_mask)

  list(pre = combined, post = combined, obs_mask = obs_mask)
}

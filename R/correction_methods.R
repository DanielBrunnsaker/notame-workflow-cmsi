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

correct_none <- function(data) {
  message("==> No correction (imputation only)")
  obs_mask <- !is.na(assay(data, 1))
  combined <- impute_rf(data, parallelize = "variables")
  list(pre = combined, post = combined, obs_mask = obs_mask)
}

correct_notame <- function(data, ruv_k) {
  # Cubic spline handles NAs natively — no LoD/2 before this step
  message("==> Drift correction (notame cubic spline)")
  combined <- merge_notame_sets(
    lapply(split_by_batch(data), process_batch),
    merge = "samples"
  )

  # Capture obs_mask after merge so column order matches combined
  obs_mask <- !is.na(assay(combined, 1))

  # LoD/2 fill before RUV which requires a complete matrix
  combined <- lod2_impute(combined)
  pre      <- combined

  n_batches <- length(unique(colData(combined)$Batch))
  if (n_batches < 2) {
    message("==> Batch correction skipped (only one batch detected)")
  } else {
    message("==> Batch correction (RUV, k=", ruv_k, ")")
    qc_idx   <- which(colData(combined)$QC == "QC")
    combined <- ruvs_qc(combined, replicates = list(qc_idx), k = ruv_k)
  }

  message("==> Imputation (RF on corrected data)")
  combined <- rf_impute_corrected(combined, obs_mask)

  list(pre = pre, post = combined, obs_mask = obs_mask)
}

correct_loess_combat <- function(data, loess_span, combat_mean_only = "auto", combat_par_prior = "auto") {
  suppressPackageStartupMessages(library(sva))

  # LOESS handles NAs natively via is.finite() — no LoD/2 before this step
  message("==> Drift correction (LOESS)")
  combined <- merge_notame_sets(
    lapply(split_by_batch(data), function(se_b) {
      loess_correct_batch(se_b, span = loess_span)
    }),
    merge = "samples"
  )

  # Capture obs_mask after merge so column order matches combined
  obs_mask <- !is.na(assay(combined, 1))

  # LoD/2 fill before ComBat which requires a complete matrix
  combined <- lod2_impute(combined)
  pre      <- combined

  n_batches <- length(unique(colData(combined)$Batch))
  if (n_batches < 2) {
    message("==> Batch correction skipped (only one batch detected)")
  } else {
    message("==> Log2 transformation")
    assay(combined, 1, withDimnames = FALSE) <- log2(assay(combined, 1))

    message("==> Between-batch correction (ComBat)")
    assay(combined, 1, withDimnames = FALSE) <- combat_correct(
      combined, mean_only = combat_mean_only, par_prior = combat_par_prior
    )

    message("==> Back-transforming to raw scale")
    assay(combined, 1, withDimnames = FALSE) <- 2^assay(combined, 1)
  }

  message("==> Imputation (RF on corrected data)")
  combined <- rf_impute_corrected(combined, obs_mask)

  list(pre = pre, post = combined, obs_mask = obs_mask)
}

correct_loess_samples_combat <- function(data, qc_span, sample_span, sample_min_obs,
                                          min_qc_per_batch = 4, min_ltqc_validate = 3,
                                          validate_samples_correction = TRUE,
                                          combat_mean_only = "auto", combat_par_prior = "auto") {
  suppressPackageStartupMessages(library(sva))

  # LOESS handles NAs natively via is.finite() — no LoD/2 before this step
  message("==> Drift correction (per-batch: QC-based if enough QC, else QC-free on samples",
          if (validate_samples_correction) " validated against ltQC, else uncorrected)"
          else " applied unconditionally -- validation disabled)")
  combined <- merge_notame_sets(
    lapply(split_by_batch(data), function(se_b) {
      loess_correct_batch_hybrid(se_b, qc_span = qc_span, sample_span = sample_span,
                                  sample_min_obs = sample_min_obs,
                                  min_qc_per_batch = min_qc_per_batch,
                                  min_ltqc_validate = min_ltqc_validate,
                                  validate = validate_samples_correction)
    }),
    merge = "samples"
  )

  # Capture obs_mask after merge so column order matches combined
  obs_mask <- !is.na(assay(combined, 1))

  # LoD/2 fill before ComBat which requires a complete matrix
  combined <- lod2_impute(combined)
  pre      <- combined

  n_batches <- length(unique(colData(combined)$Batch))
  if (n_batches < 2) {
    message("==> Batch correction skipped (only one batch detected)")
  } else {
    message("==> Log2 transformation")
    assay(combined, 1, withDimnames = FALSE) <- log2(assay(combined, 1))

    message("==> Between-batch correction (ComBat)")
    assay(combined, 1, withDimnames = FALSE) <- combat_correct(
      combined, mean_only = combat_mean_only, par_prior = combat_par_prior
    )

    message("==> Back-transforming to raw scale")
    assay(combined, 1, withDimnames = FALSE) <- 2^assay(combined, 1)
  }

  message("==> Imputation (RF on corrected data)")
  combined <- rf_impute_corrected(combined, obs_mask)

  list(pre = pre, post = combined, obs_mask = obs_mask)
}

# Huber-equivalent of correct_loess_combat(): QC-based Huber robust regression
# drift correction at a fixed k (not auto-searched -- see huber_samples_combat
# for the hybrid QC/samples version, or auto_combat if you want k chosen for
# you via CV) followed by ComBat batch correction.
correct_huber_combat <- function(data, huber_k, combat_mean_only = "auto", combat_par_prior = "auto") {
  suppressPackageStartupMessages(library(sva))

  # Huber handles NAs natively via is.finite() — no LoD/2 before this step
  message("==> Drift correction (Huber robust regression, k=", huber_k, ")")
  combined <- merge_notame_sets(
    lapply(split_by_batch(data), function(se_b) {
      huber_correct_batch(se_b, k = huber_k)
    }),
    merge = "samples"
  )

  # Capture obs_mask after merge so column order matches combined
  obs_mask <- !is.na(assay(combined, 1))

  # LoD/2 fill before ComBat which requires a complete matrix
  combined <- lod2_impute(combined)
  pre      <- combined

  n_batches <- length(unique(colData(combined)$Batch))
  if (n_batches < 2) {
    message("==> Batch correction skipped (only one batch detected)")
  } else {
    message("==> Log2 transformation")
    assay(combined, 1, withDimnames = FALSE) <- log2(assay(combined, 1))

    message("==> Between-batch correction (ComBat)")
    assay(combined, 1, withDimnames = FALSE) <- combat_correct(
      combined, mean_only = combat_mean_only, par_prior = combat_par_prior
    )

    message("==> Back-transforming to raw scale")
    assay(combined, 1, withDimnames = FALSE) <- 2^assay(combined, 1)
  }

  message("==> Imputation (RF on corrected data)")
  combined <- rf_impute_corrected(combined, obs_mask)

  list(pre = pre, post = combined, obs_mask = obs_mask)
}

# Huber-equivalent of correct_loess_samples_combat(): same three-tier per-batch
# choice (QC-based if enough QC; else a QC-free trial on samples validated
# against ltQC D-ratio; else uncorrected -- see huber_correct_batch_hybrid()),
# using fixed k values for the QC and samples fits (not auto-searched) instead
# of LOESS spans, followed by ComBat batch correction.
correct_huber_samples_combat <- function(data, qc_k, sample_k, sample_min_obs,
                                          min_qc_per_batch = 4, min_ltqc_validate = 3,
                                          validate_samples_correction = TRUE,
                                          combat_mean_only = "auto", combat_par_prior = "auto") {
  suppressPackageStartupMessages(library(sva))

  # Huber handles NAs natively via is.finite() — no LoD/2 before this step
  message("==> Drift correction (per-batch: QC-based Huber if enough QC, else QC-free Huber on samples",
          if (validate_samples_correction) " validated against ltQC, else uncorrected)"
          else " applied unconditionally -- validation disabled)")
  combined <- merge_notame_sets(
    lapply(split_by_batch(data), function(se_b) {
      huber_correct_batch_hybrid(se_b, qc_k = qc_k, sample_k = sample_k,
                                  sample_min_obs = sample_min_obs,
                                  min_qc_per_batch = min_qc_per_batch,
                                  min_ltqc_validate = min_ltqc_validate,
                                  validate = validate_samples_correction)
    }),
    merge = "samples"
  )

  # Capture obs_mask after merge so column order matches combined
  obs_mask <- !is.na(assay(combined, 1))

  # LoD/2 fill before ComBat which requires a complete matrix
  combined <- lod2_impute(combined)
  pre      <- combined

  n_batches <- length(unique(colData(combined)$Batch))
  if (n_batches < 2) {
    message("==> Batch correction skipped (only one batch detected)")
  } else {
    message("==> Log2 transformation")
    assay(combined, 1, withDimnames = FALSE) <- log2(assay(combined, 1))

    message("==> Between-batch correction (ComBat)")
    assay(combined, 1, withDimnames = FALSE) <- combat_correct(
      combined, mean_only = combat_mean_only, par_prior = combat_par_prior
    )

    message("==> Back-transforming to raw scale")
    assay(combined, 1, withDimnames = FALSE) <- 2^assay(combined, 1)
  }

  message("==> Imputation (RF on corrected data)")
  combined <- rf_impute_corrected(combined, obs_mask)

  list(pre = pre, post = combined, obs_mask = obs_mask)
}

# Auto-selected drift correction (see auto_select_drift_correction() in
# R/drift_correction.R) + ComBat between-batch correction. A new, separate
# method rather than a change to loess_combat/loess_samples_combat -- the
# intent is to validate this against the fixed-method ones before it
# potentially replaces them, not to silently change what an existing named
# method does. One method is chosen for all QC-based batches and one for all
# QC-free batches (pooling evidence across batches, not decided per batch) --
# see auto_select_drift_correction()'s own documentation for why.
correct_auto_combat <- function(data, loess_spans, huber_ks, sample_loess_spans, sample_huber_ks,
                                 min_qc_per_batch = 4, min_ltqc_validate = 3, min_cv_obs = 4,
                                 combat_mean_only = "auto", combat_par_prior = "auto") {
  suppressPackageStartupMessages(library(sva))

  message("==> Drift correction (auto-selected: one method for all QC-based batches, one for ",
          "all QC-free batches, via pooled CV on QC / held-out ltQC-Sample D-ratio)")
  combined <- merge_notame_sets(
    auto_select_drift_correction(data, loess_spans = loess_spans, huber_ks = huber_ks,
                                  sample_loess_spans = sample_loess_spans,
                                  sample_huber_ks = sample_huber_ks,
                                  min_qc_per_batch = min_qc_per_batch,
                                  min_ltqc_validate = min_ltqc_validate,
                                  min_cv_obs = min_cv_obs),
    merge = "samples"
  )

  # Capture obs_mask after merge so column order matches combined
  obs_mask <- !is.na(assay(combined, 1))

  # LoD/2 fill before ComBat which requires a complete matrix
  combined <- lod2_impute(combined)
  pre      <- combined

  n_batches <- length(unique(colData(combined)$Batch))
  if (n_batches < 2) {
    message("==> Batch correction skipped (only one batch detected)")
  } else {
    message("==> Log2 transformation")
    assay(combined, 1, withDimnames = FALSE) <- log2(assay(combined, 1))

    message("==> Between-batch correction (ComBat)")
    assay(combined, 1, withDimnames = FALSE) <- combat_correct(
      combined, mean_only = combat_mean_only, par_prior = combat_par_prior
    )

    message("==> Back-transforming to raw scale")
    assay(combined, 1, withDimnames = FALSE) <- 2^assay(combined, 1)
  }

  message("==> Imputation (RF on corrected data)")
  combined <- rf_impute_corrected(combined, obs_mask)

  list(pre = pre, post = combined, obs_mask = obs_mask)
}

correct_loess_feature_median <- function(data, loess_span) {

  # LOESS handles NAs natively — no LoD/2 needed before this step
  message("==> Drift correction (LOESS)")
  combined <- merge_notame_sets(
    lapply(split_by_batch(data), function(se_b) {
      loess_correct_batch(se_b, span = loess_span)
    }),
    merge = "samples"
  )

  obs_mask <- !is.na(assay(combined, 1))
  pre      <- combined

  n_batches <- length(unique(colData(combined)$Batch))
  if (n_batches < 2) {
    message("==> Batch correction skipped (only one batch detected)")
  } else {
    message("==> Between-batch correction (per-feature median ratio normalisation)")
    combined <- batch_feature_median_correct(combined)
  }

  message("==> Imputation (RF on corrected data)")
  combined <- rf_impute_corrected(combined, obs_mask)

  list(pre = pre, post = combined, obs_mask = obs_mask)
}

correct_loess_global_median <- function(data, loess_span) {

  # LOESS handles NAs natively — no LoD/2 needed before this step
  message("==> Drift correction (LOESS)")
  combined <- merge_notame_sets(
    lapply(split_by_batch(data), function(se_b) {
      loess_correct_batch(se_b, span = loess_span)
    }),
    merge = "samples"
  )

  obs_mask <- !is.na(assay(combined, 1))
  pre      <- combined

  n_batches <- length(unique(colData(combined)$Batch))
  if (n_batches < 2) {
    message("==> Batch correction skipped (only one batch detected)")
  } else {
    message("==> Between-batch correction (global median ratio normalisation)")
    combined <- batch_global_median_correct(combined)
  }

  message("==> Imputation (RF on corrected data)")
  combined <- rf_impute_corrected(combined, obs_mask)

  list(pre = pre, post = combined, obs_mask = obs_mask)
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

correct_loess_limma <- function(data, loess_span) {
  suppressPackageStartupMessages(library(limma))

  # LOESS handles NAs natively via is.finite() — no LoD/2 before this step
  message("==> Drift correction (LOESS)")
  combined <- merge_notame_sets(
    lapply(split_by_batch(data), function(se_b) {
      loess_correct_batch(se_b, span = loess_span)
    }),
    merge = "samples"
  )

  # Capture obs_mask after merge so column order matches combined
  obs_mask <- !is.na(assay(combined, 1))

  # LoD/2 fill before limma which requires a complete matrix
  combined <- lod2_impute(combined)
  pre      <- combined

  n_batches <- length(unique(colData(combined)$Batch))
  if (n_batches < 2) {
    message("==> Batch correction skipped (only one batch detected)")
  } else {
    message("==> Log2 transformation")
    assay(combined, 1, withDimnames = FALSE) <- log2(assay(combined, 1))

    message("==> Between-batch correction (limma removeBatchEffect)")
    assay(combined, 1, withDimnames = FALSE) <- removeBatchEffect(
      x     = assay(combined, 1),
      batch = as.factor(colData(combined)$Batch)
    )

    message("==> Back-transforming to raw scale")
    assay(combined, 1, withDimnames = FALSE) <- 2^assay(combined, 1)
  }

  message("==> Imputation (RF on corrected data)")
  combined <- rf_impute_corrected(combined, obs_mask)

  list(pre = pre, post = combined, obs_mask = obs_mask)
}

correct_loess_samples_limma <- function(data, qc_span, sample_span, sample_min_obs,
                                         min_qc_per_batch = 4, min_ltqc_validate = 3,
                                         validate_samples_correction = TRUE) {
  suppressPackageStartupMessages(library(limma))

  # LOESS handles NAs natively via is.finite() — no LoD/2 before this step
  message("==> Drift correction (per-batch: QC-based if enough QC, else QC-free on samples",
          if (validate_samples_correction) " validated against ltQC, else uncorrected)"
          else " applied unconditionally -- validation disabled)")
  combined <- merge_notame_sets(
    lapply(split_by_batch(data), function(se_b) {
      loess_correct_batch_hybrid(se_b, qc_span = qc_span, sample_span = sample_span,
                                  sample_min_obs = sample_min_obs,
                                  min_qc_per_batch = min_qc_per_batch,
                                  min_ltqc_validate = min_ltqc_validate,
                                  validate = validate_samples_correction)
    }),
    merge = "samples"
  )

  # Capture obs_mask after merge so column order matches combined
  obs_mask <- !is.na(assay(combined, 1))

  # LoD/2 fill before limma which requires a complete matrix
  combined <- lod2_impute(combined)
  pre      <- combined

  n_batches <- length(unique(colData(combined)$Batch))
  if (n_batches < 2) {
    message("==> Batch correction skipped (only one batch detected)")
  } else {
    message("==> Log2 transformation")
    assay(combined, 1, withDimnames = FALSE) <- log2(assay(combined, 1))

    message("==> Between-batch correction (limma removeBatchEffect)")
    assay(combined, 1, withDimnames = FALSE) <- removeBatchEffect(
      x     = assay(combined, 1),
      batch = as.factor(colData(combined)$Batch)
    )

    message("==> Back-transforming to raw scale")
    assay(combined, 1, withDimnames = FALSE) <- 2^assay(combined, 1)
  }

  message("==> Imputation (RF on corrected data)")
  combined <- rf_impute_corrected(combined, obs_mask)

  list(pre = pre, post = combined, obs_mask = obs_mask)
}

correct_combat_only <- function(data, combat_mean_only = "auto", combat_par_prior = "auto") {
  suppressPackageStartupMessages(library(sva))

  obs_mask <- !is.na(assay(data, 1))
  data     <- lod2_impute(data)
  pre      <- data

  message("==> Log2 transformation")
  assay(data, 1, withDimnames = FALSE) <- log2(assay(data, 1))

  message("==> Batch correction (ComBat only, no drift correction)")
  assay(data, 1, withDimnames = FALSE) <- combat_correct(
    data, mean_only = combat_mean_only, par_prior = combat_par_prior
  )

  message("==> Back-transforming to raw scale")
  assay(data, 1, withDimnames = FALSE) <- 2^assay(data, 1)

  message("==> Imputation (RF on corrected data)")
  combined <- rf_impute_corrected(data, obs_mask)

  list(pre = pre, post = combined, obs_mask = obs_mask)
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

correct_pmp_qcrsc_scale <- function(data) {
  suppressPackageStartupMessages(library(pmp))

  obs_mask <- !is.na(assay(data, 1))

  # Identify batches pmp cannot anchor (< 4 QC samples).
  # pmp's between-batch alignment step distorts these batches — save their
  # original values and restore them after QCRSC runs.
  cd           <- as.data.frame(colData(data))
  qc_counts    <- tapply(cd$QC == "QC", as.character(cd$Batch), sum)
  skip_batches <- names(qc_counts[qc_counts < 4])
  orig_mat     <- assay(data, 1)

  # Step 1: QC-RSC — drift correction + QC-anchored between-batch alignment
  # for batches with >= 4 QC samples.
  message("==> Drift correction + QC-based batch alignment (pmp QC-RSC)")
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

  # Clamp spline overshoot in corrected batches before any log step
  combined <- clamp_nonpositive(combined, "after QC-RSC")
  pre      <- combined

  # Step 2: Global median scaling — only for batches pmp could not align.
  # Scales each no-QC batch by a single factor so its overall intensity level
  # matches the grand median of the pmp-corrected batches. One scalar per batch:
  # no feature-specific adjustments, no disturbance to already-corrected batches.
  # Defensible when samples are randomised (global shift is technical, not biological).
  if (length(skip_batches) > 0) {
    mat          <- assay(combined, 1)
    good_batches <- setdiff(unique(as.character(cd$Batch)), skip_batches)
    good_samp    <- which(cd$QC == "Sample" & as.character(cd$Batch) %in% good_batches)
    ref_vals     <- mat[, good_samp, drop = FALSE]
    grand_med    <- median(ref_vals[is.finite(ref_vals) & ref_vals > 0], na.rm = TRUE)

    message("==> Global median scaling for batch(es) with <4 QCs")
    for (b in skip_batches) {
      b_samp <- which(cd$QC == "Sample" & as.character(cd$Batch) == b)
      if (length(b_samp) < 2) {
        message("  Batch ", b, ": skipped (fewer than 2 biological samples)"); next
      }
      b_vals    <- mat[, b_samp, drop = FALSE]
      batch_med <- median(b_vals[is.finite(b_vals) & b_vals > 0], na.rm = TRUE)
      if (!is.finite(batch_med) || batch_med <= 0) {
        message("  Batch ", b, ": skipped (could not compute median)"); next
      }
      scale_fac <- grand_med / batch_med
      b_all     <- which(as.character(cd$Batch) == b)
      assay(combined, 1, withDimnames = FALSE)[, b_all] <-
        assay(combined, 1)[, b_all] * scale_fac
      message("  Batch ", b, ": scale factor = ", round(scale_fac, 4),
              "  (batch median ", round(batch_med), " → grand median ", round(grand_med), ")")
    }
  }

  message("==> Imputation (RF on corrected data)")
  combined <- rf_impute_corrected(combined, obs_mask)

  list(pre = pre, post = combined, obs_mask = obs_mask)
}

# Per-feature median scaling variant of pmp_qcrsc_scale.
# Identical to pmp_qcrsc_scale except that between-batch alignment for no-QC
# batches is done per feature (each feature scaled to the grand median of that
# feature across QC-corrected batches) rather than with a single global scalar.
# This is consistent with how pmp QC-RSC aligns QC-anchored batches, where each
# feature is referenced to its own global QC median. The only methodological
# difference from the fully corrected batches is that within-batch drift
# correction is skipped (no QC samples to anchor the spline).
correct_pmp_qcrsc_feature_scale <- function(data) {
  suppressPackageStartupMessages(library(pmp))

  obs_mask <- !is.na(assay(data, 1))

  cd           <- as.data.frame(colData(data))
  qc_counts    <- tapply(cd$QC == "QC", as.character(cd$Batch), sum)
  skip_batches <- names(qc_counts[qc_counts < 4])
  orig_mat     <- assay(data, 1)

  message("==> Drift correction + QC-based batch alignment (pmp QC-RSC)")
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

  combined <- clamp_nonpositive(combined, "after QC-RSC")
  pre      <- combined

  # Per-feature median scaling for no-QC batches.
  # For each feature, compute its median across biological samples in the
  # QC-corrected batches and scale the no-QC batch to match.
  # Consistent with pmp's own feature-wise QC-median referencing.
  if (length(skip_batches) > 0) {
    mat          <- assay(combined, 1)
    good_batches <- setdiff(unique(as.character(cd$Batch)), skip_batches)
    good_samp    <- which(cd$QC == "Sample" & as.character(cd$Batch) %in% good_batches)

    # Per-feature reference median from QC-corrected batches
    ref_medians <- apply(mat[, good_samp, drop = FALSE], 1, function(x) {
      x <- x[is.finite(x) & x > 0]
      if (length(x) < 2) NA_real_ else median(x)
    })

    message("==> Per-feature median scaling for batch(es) with <4 QCs")
    for (b in skip_batches) {
      b_samp <- which(cd$QC == "Sample" & as.character(cd$Batch) == b)
      if (length(b_samp) < 2) {
        message("  Batch ", b, ": skipped (fewer than 2 biological samples)"); next
      }

      batch_medians <- apply(mat[, b_samp, drop = FALSE], 1, function(x) {
        x <- x[is.finite(x) & x > 0]
        if (length(x) < 2) NA_real_ else median(x)
      })

      scale_facs <- ref_medians / batch_medians
      # Features where scaling cannot be computed keep their original values
      scale_facs[!is.finite(scale_facs) | scale_facs <= 0] <- 1

      b_all <- which(as.character(cd$Batch) == b)
      mat[, b_all] <- mat[, b_all] * scale_facs
      assay(combined, 1, withDimnames = FALSE) <- mat

      n_scaled <- sum(scale_facs != 1)
      message("  Batch ", b, ": per-feature scaling applied to ", n_scaled, "/",
              nrow(mat), " features  (median scale factor = ",
              round(median(scale_facs[scale_facs != 1]), 4), ")")
    }
  }

  message("==> Imputation (RF on corrected data)")
  combined <- rf_impute_corrected(combined, obs_mask)

  list(pre = pre, post = combined, obs_mask = obs_mask)
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

correct_cordbat_only <- function(data, ref_batch = NULL) {
  obs_mask <- !is.na(assay(data, 1))
  data     <- lod2_impute(data)
  pre      <- data

  data <- run_cordbat(data, ref_batch)

  message("==> Imputation (RF on corrected data)")
  combined <- rf_impute_corrected(data, obs_mask)

  list(pre = pre, post = combined, obs_mask = obs_mask)
}

correct_loess_cordbat <- function(data, sample_span, sample_min_obs, ref_batch = NULL) {
  # Always the QC-free samples-based fit, deliberately — CordBat's own
  # between-batch correction (run_cordbat -> CordBat()) is already fit on
  # biological samples, not QC (see Funcs_CordBat_algorithm.R: QC rows are
  # stripped out before the GGM fit and only corrected afterward using
  # coefficients learned from Samples). Using a QC-anchored drift step ahead
  # of it would make this combined method depend on QC for drift but not for
  # batch correction — an inconsistent mix. This keeps loess_cordbat QC-free
  # end to end, so it corrects every batch the same way regardless of QC
  # coverage (unlike loess_combat/loess_samples_combat's QC-preferring hybrid).
  #
  # LOESS handles NAs natively — no LoD/2 before this step
  message("==> Drift correction (QC-free LOESS on biological samples)")
  combined <- merge_notame_sets(
    lapply(split_by_batch(data), function(se_b) {
      loess_correct_batch_samples(se_b, span = sample_span, min_obs = sample_min_obs)
    }),
    merge = "samples"
  )

  obs_mask <- !is.na(assay(combined, 1))

  # LoD/2 fill before CordBat which requires a complete matrix
  combined <- lod2_impute(combined)
  pre      <- combined

  combined <- run_cordbat(combined, ref_batch)

  message("==> Imputation (RF on corrected data)")
  combined <- rf_impute_corrected(combined, obs_mask)

  list(pre = pre, post = combined, obs_mask = obs_mask)
}



# K = NULL auto-derives the component count as before (2 per batch); pass an
# explicit value to override. alpha/Cutoff/wf are forwarded to WaveICA_2.0()
# as-is — see notame-workflow.r's WAVEICA_* help text for what each does and
# which direction to try if correction looks too aggressive ("flattens" real
# sample variation) vs. too weak (batch effect still visible after correction).
correct_waveica <- function(data, alpha = 0.05, cutoff = 0.10, K = NULL, wf = "haar") {
  suppressPackageStartupMessages(library(WaveICA2.0))

  obs_mask <- !is.na(assay(data, 1))
  data     <- lod2_impute(data)

  k_eff <- if (is.null(K)) length(unique(colData(data)$Batch)) * 2 else K
  message("==> WaveICA2.0 correction (alpha=", alpha, ", Cutoff=", cutoff,
          ", K=", k_eff, ", wf=", wf, ")")
  corrected_mat <- WaveICA_2.0(
    data            = t(assay(data, 1)),
    wf              = wf,
    Injection_Order = as.numeric(colData(data)$Injection_order),
    alpha           = alpha,
    Cutoff          = cutoff,
    K               = k_eff
  )
  assay(data, 1, withDimnames = FALSE) <- t(corrected_mat$data)

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
# variable-wise independence in the ICA step) is unrelated to
# correct_waveica()'s `alpha` (0-1, significance threshold for flagging a
# component as injection-order-associated) -- same parameter name in both
# packages, different meaning. Kept as WAVEICA_V1_ALPHA (distinct from
# WAVEICA_ALPHA) in notame-workflow.r specifically to avoid conflating them.
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

# ─────────────────────────────────────────────────────────────────────────────
# Within-batch drift correction helpers.
#
# split_by_batch()              — splits a SE into a list of per-batch SEs, sorted by injection order
# process_batch()                — notame cubic spline drift correction wrapper
# loess_correct_batch()          — QC-based LOESS drift correction
# loess_correct_batch_samples()  — QC-free LOESS drift correction (fit on biological samples)
# loess_correct_batch_hybrid()   — per-batch: QC-based if enough QC, else samples-based
#                                   trial validated against ltQC (kept only if it helps), else uncorrected
# huber_correct_batch()          — QC-based Huber robust regression drift correction (fixed k)
# huber_correct_batch_samples()  — QC-free Huber drift correction, fit on biological samples
# huber_correct_batch_hybrid()   — per-batch Huber equivalent of loess_correct_batch_hybrid()
# spline_correct_batch()         — QC-based cubic smoothing spline drift correction (fixed spar)
# spline_correct_batch_samples() — QC-free smoothing spline drift correction, fit on samples
# spline_correct_batch_hybrid()  — per-batch smoothing spline equivalent of loess_correct_batch_hybrid()
# qc_cv_correct_batch()          — QC-based drift correction with span/k/spar selected per FEATURE
#                                   via leave-one-out CV, instead of one shared value (opt-in, via
#                                   loess_correct_batch_hybrid()/huber_correct_batch_hybrid()/
#                                   spline_correct_batch_hybrid()'s qc_span_grid/qc_k_grid/
#                                   qc_spar_grid) -- mirrors notame::correct_drift()
# gate_by_qc_count()              — applies a per-batch drift function only if the batch meets a
#                                    minimum QC count; used to give basis="qc" the same batch-level
#                                    gate basis="hybrid" already has (see R/method_spec.R)
# resolve_drift_fn()              — dispatch table: (drift_method, basis) -> function(se_b, params),
#                                    used by run_correction() (R/correction_methods.R) for every
#                                    drift x basis combination except "auto" (handled separately,
#                                    since it operates over all batches at once)
# candidate_family()               — maps a candidate's name to its family ("loess"/"huber"/
#                                     "spline"/"flat"), used to group per-parameter candidates for
#                                     rank_select_family()
# rank_select_family()              — pre-stage: pools evidence across all eligible batches, via
#                                      rank-sum rather than raw score magnitude, to pick ONE family
#                                      to restrict the per-batch search to -- guards against a
#                                      high-noise or feature-rich batch dominating a
#                                      magnitude-based pool
# select_qc_candidate_per_batch()    — for each batch independently, evaluates per-feature
#                                      LOO-CV-on-QC scores and picks that batch's own winning
#                                      candidate (restricted to rank_select_family()'s pick);
#                                      extracted from auto_select_drift_correction() so it can be
#                                      reused for any basis
# select_sample_candidate_per_batch() — same idea, per batch, using held-out ltQC prediction score
#                                      (D-ratio computed alongside but informational only) for a
#                                      samples-based pick
# auto_select_drift_correction() — basis-parameterized (qc/samples/hybrid): for each eligible
#                                   QC-based batch, and each eligible QC-free batch, independently
#                                   selects its OWN best drift-correction method (a different
#                                   method per batch is expected, not pooled into one shared
#                                   winner), from several fitting candidates (LOESS at several
#                                   spans, Huber regression at several k, a smoothing spline at
#                                   several spar values, a flat/no-op baseline),
#                                   using that batch's own evidence: leave-one-out CV on its own QC
#                                   for QC-based batches, held-out ltQC prediction score on its own
#                                   samples for QC-free batches. Every batch's search is restricted
#                                   to one dataset-wide family (rank_select_family()), leaving
#                                   per-batch parameter tuning within that family unchanged.
#                                   Returns list(batches=, log=): batches is the list of corrected
#                                   per-batch SEs, log is a
#                                   data.frame (one row per batch x candidate) of the full selection
#                                   record, for the caller to save alongside the corrected output.

split_by_batch <- function(se) {
  batches <- unique(colData(se)$Batch)
  out <- list()
  for (b in batches) {
    idx  <- which(colData(se)$Batch == b)
    se_b <- se[, idx]
    se_b <- se_b[, order(colData(se_b)$Injection_order), drop = FALSE]
    out[[paste0("Batch_", b)]] <- se_b
  }
  out
}


# notame cubic spline drift correction (wraps notame::correct_drift).
process_batch <- function(se_b) {
  message("==> Processing Batch ", unique(colData(se_b)$Batch))
  correct_drift(se_b)
}


# LOESS-based within-batch drift correction.
# Fits a LOESS curve through QC samples (injection order vs feature abundance),
# then divides all samples by the predicted value normalised to QC median.
# Requires >= 4 finite QC observations per feature to fit reliably.
loess_correct_batch <- function(se_b, span = 0.75) {
  mat    <- assay(se_b, 1)
  cd     <- colData(se_b)
  qc_idx <- which(cd$QC == "QC")
  inj    <- as.numeric(cd$Injection_order)
  batch  <- unique(cd$Batch)
  n_feat <- nrow(mat)

  message("  Batch ", batch, ": ", length(qc_idx), " QC sample(s), ", n_feat, " features")

  if (any(!is.finite(inj))) {
    bad <- which(!is.finite(inj))
    stop("Non-finite Injection_order in batch ", batch,
         ": samples ", paste(cd$Sample_ID[bad], collapse = ", "),
         " (values: ", paste(inj[bad], collapse = ", "), ")")
  }

  n_skipped_qc  <- 0L
  n_skipped_err <- 0L

  for (i in seq_len(nrow(mat))) {
    y_qc <- as.numeric(mat[i, qc_idx])
    x_qc <- inj[qc_idx]
    # y_qc > 0: required for log2() below; also excludes exact-zero
    # abundances (msdial_to_notame()/xcms_to_notame() clamp negative values
    # to 0, not NA) from the fit.
    ok   <- is.finite(y_qc) & y_qc > 0

    if (sum(ok) < 4) { n_skipped_qc <- n_skipped_qc + 1L; next }

    tryCatch({
      # Fit in log2 space, not raw abundance: real LC-MS drift is often
      # multiplicative/exponential-shaped (y = y0 * f(t)), which LOESS's
      # local quadratic fit can partially but not fully track on the raw
      # scale. log2(y) turns that into a much better-conditioned shape for
      # local polynomial fitting, same reasoning as fit_predict_loess() in
      # auto_select_drift_correction() and why ComBat/limma always operate
      # in log space.
      fit          <- loess(y ~ x, data = data.frame(x = x_qc[ok], y = log2(y_qc[ok])), span = span)
      ok_inj       <- !is.na(inj)
      pred         <- rep(NA_real_, length(inj))
      pred[ok_inj] <- 2^predict(fit, newdata = data.frame(x = inj[ok_inj]))
      med_qc       <- median(y_qc[ok])
      ratio        <- pred / med_qc
      ratio[is.na(ratio) | ratio <= 0] <- 1
      mat[i, ]     <- mat[i, ] / ratio
    }, error = function(e) { n_skipped_err <<- n_skipped_err + 1L })
  }

  n_qc_corrected <- n_feat - n_skipped_qc - n_skipped_err
  message("  Batch ", batch, ": drift-corrected ", n_qc_corrected, "/", n_feat, " features",
          if (n_skipped_qc  > 0) paste0(" | ", n_skipped_qc,  " skipped (insufficient QC observations)") else "",
          if (n_skipped_err > 0) paste0(" | ", n_skipped_err, " skipped (LOESS fit error)") else "")

  assay(se_b, 1, withDimnames = FALSE) <- mat
  se_b
}


# Per-batch hybrid drift correction:
#   1. Enough QC (>= min_qc_per_batch)   -> QC-anchored LOESS (loess_correct_batch).
#      Preferred whenever possible — it derives the drift curve from technical
#      replicates rather than biological samples, so it doesn't risk removing
#      real biological signal along with drift (see loess_correct_batch_samples's
#      documentation).
#   2. Not enough QC, but enough ltQC (>= min_ltqc_validate) to check the result
#      -> trial the QC-free samples-based fit (loess_correct_batch_samples), then
#      keep it only if it measurably improves the ltQC/Sample D-ratio
#      (eval_ltqc_dratio = MAD(ltQC)/MAD(Sample), from R/qc_metrics.R) in that
#      batch; otherwise revert to the uncorrected values. ltQC is never used to
#      fit the samples-based correction, so this is a genuine held-out check
#      rather than a circular one. D-ratio (not raw ltQC RSD) is the right test
#      here: any real drift correction shrinks measured sample variance somewhat
#      (removing genuine drift noise does that even when biological signal is
#      fully preserved), so "sample variance must not shrink" would reject
#      working corrections too. D-ratio only credits a *disproportionate*
#      shrink in ltQC relative to Sample — proportional shrinkage in both
#      (shared drift removed cleanly) leaves the ratio roughly flat, while
#      Sample shrinking as much as or more than ltQC (real signal being
#      removed) leaves it unimproved or worse.
#   3. Neither -> leave the batch uncorrected; there is no data to validate a
#      QC-free fit against, and an unverifiable correction is worse than none.
# QC coverage (and ltQC coverage) can vary batch to batch even within one
# dataset (e.g. a plate with zero QC injections), so the choice is made per
# batch rather than once for the whole run.
#
# validate = FALSE skips step 2's ltQC check entirely and always keeps the
# samples-based trial, regardless of ltQC availability or what it shows. This
# reintroduces the failure mode the validation step exists to catch (the
# trial can look fine on ltQC while still compressing real biological
# signal) -- use deliberately, not as a default.
loess_correct_batch_hybrid <- function(se_b, qc_span, sample_span, sample_min_obs,
                                        min_qc_per_batch = 4, min_ltqc_validate = 3,
                                        validate = TRUE, qc_span_grid = numeric(0)) {
  cd    <- colData(se_b)
  n_qc  <- sum(cd$QC == "QC")
  batch <- unique(cd$Batch)

  if (n_qc >= min_qc_per_batch) {
    if (length(qc_span_grid) > 0) {
      message("  Batch ", batch, ": ", n_qc, " QC sample(s) (>= ", min_qc_per_batch,
              ") — using QC-based LOESS, span selected per feature via CV")
      return(qc_cv_correct_batch(se_b, loess_spans = qc_span_grid, min_obs = min_qc_per_batch))
    }
    message("  Batch ", batch, ": ", n_qc, " QC sample(s) (>= ", min_qc_per_batch,
            ") — using QC-based LOESS")
    return(loess_correct_batch(se_b, span = qc_span))
  }

  if (!validate) {
    message("  Batch ", batch, ": ", n_qc, " QC sample(s) (< ", min_qc_per_batch,
            ") — applying QC-free LOESS on samples unconditionally (validation disabled)")
    return(loess_correct_batch_samples(se_b, span = sample_span, min_obs = sample_min_obs))
  }

  n_ltqc <- sum(cd$QC == "ltQC")
  if (n_ltqc < min_ltqc_validate) {
    message("  Batch ", batch, ": ", n_qc, " QC sample(s) (< ", min_qc_per_batch,
            ") and only ", n_ltqc, " ltQC sample(s) (< ", min_ltqc_validate,
            ") to validate a QC-free fit — leaving batch uncorrected")
    return(se_b)
  }

  message("  Batch ", batch, ": ", n_qc, " QC sample(s) (< ", min_qc_per_batch,
          ") — trialing QC-free LOESS on samples, to be validated against ",
          n_ltqc, " ltQC sample(s)")
  se_trial <- loess_correct_batch_samples(se_b, span = sample_span, min_obs = sample_min_obs)

  dratio_before <- eval_ltqc_dratio(se_b)
  dratio_after  <- eval_ltqc_dratio(se_trial)

  if (is.na(dratio_before) || is.na(dratio_after)) {
    message("  Batch ", batch, ": could not compute ltQC/Sample D-ratio before/after (NA) — ",
            "leaving batch uncorrected")
    return(se_b)
  }

  if (dratio_after < dratio_before) {
    message("  Batch ", batch, ": ltQC/Sample D-ratio improved with QC-free correction (",
            round(dratio_before, 4), " -> ", round(dratio_after, 4), ") — keeping correction")
    se_trial
  } else {
    message("  Batch ", batch, ": ltQC/Sample D-ratio did not improve with QC-free correction (",
            round(dratio_before, 4), " -> ", round(dratio_after, 4), ") — reverting to uncorrected")
    se_b
  }
}


# QC-free LOESS drift correction.
# Fits a robust LOESS curve (family = "symmetric", i.e. iteratively reweighted
# to downweight outliers) through biological SAMPLE observations — not QC —
# vs. injection order, then divides all samples by the predicted value
# normalised to the sample median. Used when QC samples are unusable
# (degraded, missing, or insufficient).
#
# This relies on samples being randomised within run order: that makes the
# fitted trend unbiased with respect to any biological grouping in
# expectation, but each point is still a unique biological measurement
# rather than a technical replicate of the same pool, so it is far noisier
# than QC-anchored fitting. The wider default span and robust family guard
# against fitting individual-sample noise/outliers as if they were drift —
# both matter and are not redundant: span controls how many points are
# aggregated per local fit, family = "symmetric" controls how much a single
# extreme value within that window can pull the fit.
# Requires >= min_obs finite sample observations per feature (deliberately
# higher than QC's threshold of 4, since each point carries far more noise).
loess_correct_batch_samples <- function(se_b, span = 0.9, min_obs = 10) {
  mat      <- assay(se_b, 1)
  cd       <- colData(se_b)
  samp_idx <- which(cd$QC == "Sample")
  inj      <- as.numeric(cd$Injection_order)
  batch    <- unique(cd$Batch)
  n_feat   <- nrow(mat)

  message("  Batch ", batch, ": ", length(samp_idx), " sample(s), ", n_feat, " features")

  if (any(!is.finite(inj))) {
    bad <- which(!is.finite(inj))
    stop("Non-finite Injection_order in batch ", batch,
         ": samples ", paste(cd$Sample_ID[bad], collapse = ", "),
         " (values: ", paste(inj[bad], collapse = ", "), ")")
  }

  n_skipped_n   <- 0L
  n_skipped_err <- 0L

  for (i in seq_len(nrow(mat))) {
    y_s <- as.numeric(mat[i, samp_idx])
    x_s <- inj[samp_idx]
    # y_s > 0: required for log2() below; see loess_correct_batch()'s
    # comment on the same filter.
    ok  <- is.finite(y_s) & y_s > 0

    if (sum(ok) < min_obs) { n_skipped_n <- n_skipped_n + 1L; next }

    tryCatch({
      # Fit in log2 space, not raw abundance -- see loess_correct_batch()'s
      # comment for why.
      fit          <- loess(y ~ x, data = data.frame(x = x_s[ok], y = log2(y_s[ok])),
                            span = span, family = "symmetric")
      ok_inj       <- !is.na(inj)
      pred         <- rep(NA_real_, length(inj))
      pred[ok_inj] <- 2^predict(fit, newdata = data.frame(x = inj[ok_inj]))
      med_s        <- median(y_s[ok])
      ratio        <- pred / med_s
      ratio[is.na(ratio) | ratio <= 0] <- 1
      mat[i, ]     <- mat[i, ] / ratio
    }, error = function(e) { n_skipped_err <<- n_skipped_err + 1L })
  }

  n_corrected <- n_feat - n_skipped_n - n_skipped_err
  message("  Batch ", batch, ": drift-corrected ", n_corrected, "/", n_feat, " features",
          if (n_skipped_n   > 0) paste0(" | ", n_skipped_n,   " skipped (insufficient sample observations)") else "",
          if (n_skipped_err > 0) paste0(" | ", n_skipped_err, " skipped (LOESS fit error)") else "")

  assay(se_b, 1, withDimnames = FALSE) <- mat
  se_b
}


# Huber robust regression within-batch drift correction (QC-based).
# Structurally identical to loess_correct_batch() -- same per-feature QC fit,
# same ratio-to-QC-median correction -- but fits a single robust linear trend
# (MASS::rlm, psi.huber, tuning constant k) instead of a local LOESS curve.
# Rigid where LOESS is flexible: cannot track curved drift, but is far less
# prone to fitting noise as signal at small QC counts, and is the more
# defensible choice when QC counts are on the low side for a stable LOESS fit.
# Uses fit_predict_huber() (this file) for the actual fit -- same log2-space
# rationale documented there.
huber_correct_batch <- function(se_b, k = 1.345) {
  mat    <- assay(se_b, 1)
  cd     <- colData(se_b)
  qc_idx <- which(cd$QC == "QC")
  inj    <- as.numeric(cd$Injection_order)
  batch  <- unique(cd$Batch)
  n_feat <- nrow(mat)

  message("  Batch ", batch, ": ", length(qc_idx), " QC sample(s), ", n_feat, " features")

  if (any(!is.finite(inj))) {
    bad <- which(!is.finite(inj))
    stop("Non-finite Injection_order in batch ", batch,
         ": samples ", paste(cd$Sample_ID[bad], collapse = ", "),
         " (values: ", paste(inj[bad], collapse = ", "), ")")
  }

  n_skipped_qc  <- 0L
  n_skipped_err <- 0L

  for (i in seq_len(nrow(mat))) {
    y_qc <- as.numeric(mat[i, qc_idx])
    x_qc <- inj[qc_idx]
    ok   <- is.finite(y_qc) & y_qc > 0

    if (sum(ok) < 4) { n_skipped_qc <- n_skipped_qc + 1L; next }

    tryCatch({
      ok_inj       <- !is.na(inj)
      pred         <- rep(NA_real_, length(inj))
      pred[ok_inj] <- fit_predict_huber(x_qc[ok], y_qc[ok], inj[ok_inj], k = k)
      med_qc       <- median(y_qc[ok])
      ratio        <- pred / med_qc
      ratio[is.na(ratio) | ratio <= 0] <- 1
      mat[i, ]     <- mat[i, ] / ratio
    }, error = function(e) { n_skipped_err <<- n_skipped_err + 1L })
  }

  n_qc_corrected <- n_feat - n_skipped_qc - n_skipped_err
  message("  Batch ", batch, ": drift-corrected ", n_qc_corrected, "/", n_feat, " features",
          if (n_skipped_qc  > 0) paste0(" | ", n_skipped_qc,  " skipped (insufficient QC observations)") else "",
          if (n_skipped_err > 0) paste0(" | ", n_skipped_err, " skipped (Huber fit error)") else "")

  assay(se_b, 1, withDimnames = FALSE) <- mat
  se_b
}


# Huber robust regression within-batch drift correction (QC-free, fit on
# biological samples). Structurally identical to loess_correct_batch_samples()
# -- same reasoning about samples being noisier than QC and the higher default
# min_obs -- but fits MASS::rlm instead of a robust-family LOESS. Huber's own
# psi.huber downweighting already guards against single-sample outliers, the
# same role LOESS's family = "symmetric" plays there.
huber_correct_batch_samples <- function(se_b, k = 1.345, min_obs = 10) {
  mat      <- assay(se_b, 1)
  cd       <- colData(se_b)
  samp_idx <- which(cd$QC == "Sample")
  inj      <- as.numeric(cd$Injection_order)
  batch    <- unique(cd$Batch)
  n_feat   <- nrow(mat)

  message("  Batch ", batch, ": ", length(samp_idx), " sample(s), ", n_feat, " features")

  if (any(!is.finite(inj))) {
    bad <- which(!is.finite(inj))
    stop("Non-finite Injection_order in batch ", batch,
         ": samples ", paste(cd$Sample_ID[bad], collapse = ", "),
         " (values: ", paste(inj[bad], collapse = ", "), ")")
  }

  n_skipped_n   <- 0L
  n_skipped_err <- 0L

  for (i in seq_len(nrow(mat))) {
    y_s <- as.numeric(mat[i, samp_idx])
    x_s <- inj[samp_idx]
    ok  <- is.finite(y_s) & y_s > 0

    if (sum(ok) < min_obs) { n_skipped_n <- n_skipped_n + 1L; next }

    tryCatch({
      ok_inj       <- !is.na(inj)
      pred         <- rep(NA_real_, length(inj))
      pred[ok_inj] <- fit_predict_huber(x_s[ok], y_s[ok], inj[ok_inj], k = k)
      med_s        <- median(y_s[ok])
      ratio        <- pred / med_s
      ratio[is.na(ratio) | ratio <= 0] <- 1
      mat[i, ]     <- mat[i, ] / ratio
    }, error = function(e) { n_skipped_err <<- n_skipped_err + 1L })
  }

  n_corrected <- n_feat - n_skipped_n - n_skipped_err
  message("  Batch ", batch, ": drift-corrected ", n_corrected, "/", n_feat, " features",
          if (n_skipped_n   > 0) paste0(" | ", n_skipped_n,   " skipped (insufficient sample observations)") else "",
          if (n_skipped_err > 0) paste0(" | ", n_skipped_err, " skipped (Huber fit error)") else "")

  assay(se_b, 1, withDimnames = FALSE) <- mat
  se_b
}


# Per-batch hybrid Huber drift correction -- identical decision logic to
# loess_correct_batch_hybrid() (QC-based if enough QC; else a samples-based
# trial validated against ltQC D-ratio; else uncorrected), substituting the
# Huber pair above for the LOESS pair. See loess_correct_batch_hybrid()'s
# documentation for the full rationale, which applies unchanged here.
huber_correct_batch_hybrid <- function(se_b, qc_k, sample_k, sample_min_obs,
                                        min_qc_per_batch = 4, min_ltqc_validate = 3,
                                        validate = TRUE, qc_k_grid = numeric(0)) {
  cd    <- colData(se_b)
  n_qc  <- sum(cd$QC == "QC")
  batch <- unique(cd$Batch)

  if (n_qc >= min_qc_per_batch) {
    if (length(qc_k_grid) > 0) {
      message("  Batch ", batch, ": ", n_qc, " QC sample(s) (>= ", min_qc_per_batch,
              ") — using QC-based Huber, k selected per feature via CV")
      return(qc_cv_correct_batch(se_b, huber_ks = qc_k_grid, min_obs = min_qc_per_batch))
    }
    message("  Batch ", batch, ": ", n_qc, " QC sample(s) (>= ", min_qc_per_batch,
            ") — using QC-based Huber")
    return(huber_correct_batch(se_b, k = qc_k))
  }

  if (!validate) {
    message("  Batch ", batch, ": ", n_qc, " QC sample(s) (< ", min_qc_per_batch,
            ") — applying QC-free Huber on samples unconditionally (validation disabled)")
    return(huber_correct_batch_samples(se_b, k = sample_k, min_obs = sample_min_obs))
  }

  n_ltqc <- sum(cd$QC == "ltQC")
  if (n_ltqc < min_ltqc_validate) {
    message("  Batch ", batch, ": ", n_qc, " QC sample(s) (< ", min_qc_per_batch,
            ") and only ", n_ltqc, " ltQC sample(s) (< ", min_ltqc_validate,
            ") to validate a QC-free fit — leaving batch uncorrected")
    return(se_b)
  }

  message("  Batch ", batch, ": ", n_qc, " QC sample(s) (< ", min_qc_per_batch,
          ") — trialing QC-free Huber on samples, to be validated against ",
          n_ltqc, " ltQC sample(s)")
  se_trial <- huber_correct_batch_samples(se_b, k = sample_k, min_obs = sample_min_obs)

  dratio_before <- eval_ltqc_dratio(se_b)
  dratio_after  <- eval_ltqc_dratio(se_trial)

  if (is.na(dratio_before) || is.na(dratio_after)) {
    message("  Batch ", batch, ": could not compute ltQC/Sample D-ratio before/after (NA) — ",
            "leaving batch uncorrected")
    return(se_b)
  }

  if (dratio_after < dratio_before) {
    message("  Batch ", batch, ": ltQC/Sample D-ratio improved with QC-free correction (",
            round(dratio_before, 4), " -> ", round(dratio_after, 4), ") — keeping correction")
    se_trial
  } else {
    message("  Batch ", batch, ": ltQC/Sample D-ratio did not improve with QC-free correction (",
            round(dratio_before, 4), " -> ", round(dratio_after, 4), ") — reverting to uncorrected")
    se_b
  }
}


# Smoothing spline within-batch drift correction (QC-based). Structurally
# identical to huber_correct_batch() -- same per-feature QC fit, same
# ratio-to-QC-median correction, calls the shared fit_predict_spline() helper
# directly since (like Huber) there's no separate robust-family variant to
# choose between for QC vs samples. More flexible than Huber's rigid line,
# smoother/more globally continuous than LOESS's local fit -- see
# fit_predict_spline()'s documentation for the full tradeoff.
spline_correct_batch <- function(se_b, spar = 0.6) {
  mat    <- assay(se_b, 1)
  cd     <- colData(se_b)
  qc_idx <- which(cd$QC == "QC")
  inj    <- as.numeric(cd$Injection_order)
  batch  <- unique(cd$Batch)
  n_feat <- nrow(mat)

  message("  Batch ", batch, ": ", length(qc_idx), " QC sample(s), ", n_feat, " features")

  if (any(!is.finite(inj))) {
    bad <- which(!is.finite(inj))
    stop("Non-finite Injection_order in batch ", batch,
         ": samples ", paste(cd$Sample_ID[bad], collapse = ", "),
         " (values: ", paste(inj[bad], collapse = ", "), ")")
  }

  n_skipped_qc  <- 0L
  n_skipped_err <- 0L

  for (i in seq_len(nrow(mat))) {
    y_qc <- as.numeric(mat[i, qc_idx])
    x_qc <- inj[qc_idx]
    ok   <- is.finite(y_qc) & y_qc > 0

    if (sum(ok) < 4) { n_skipped_qc <- n_skipped_qc + 1L; next }

    tryCatch({
      ok_inj       <- !is.na(inj)
      pred         <- rep(NA_real_, length(inj))
      pred[ok_inj] <- fit_predict_spline(x_qc[ok], y_qc[ok], inj[ok_inj], spar = spar)
      med_qc       <- median(y_qc[ok])
      ratio        <- pred / med_qc
      ratio[is.na(ratio) | ratio <= 0] <- 1
      mat[i, ]     <- mat[i, ] / ratio
    }, error = function(e) { n_skipped_err <<- n_skipped_err + 1L })
  }

  n_qc_corrected <- n_feat - n_skipped_qc - n_skipped_err
  message("  Batch ", batch, ": drift-corrected ", n_qc_corrected, "/", n_feat, " features",
          if (n_skipped_qc  > 0) paste0(" | ", n_skipped_qc,  " skipped (insufficient QC observations)") else "",
          if (n_skipped_err > 0) paste0(" | ", n_skipped_err, " skipped (spline fit error)") else "")

  assay(se_b, 1, withDimnames = FALSE) <- mat
  se_b
}


# Smoothing spline within-batch drift correction (QC-free, fit on biological
# samples). Structurally identical to huber_correct_batch_samples() -- same
# reasoning about samples being noisier than QC and the higher default
# min_obs. No robust-family equivalent exists for smooth.spline (see
# fit_predict_spline()'s documentation), so a higher default spar (more
# smoothing, less local reactivity) is this candidate's main defense against
# a single biological outlier pulling the fit, on top of min_obs.
spline_correct_batch_samples <- function(se_b, spar = 0.8, min_obs = 10) {
  mat      <- assay(se_b, 1)
  cd       <- colData(se_b)
  samp_idx <- which(cd$QC == "Sample")
  inj      <- as.numeric(cd$Injection_order)
  batch    <- unique(cd$Batch)
  n_feat   <- nrow(mat)

  message("  Batch ", batch, ": ", length(samp_idx), " sample(s), ", n_feat, " features")

  if (any(!is.finite(inj))) {
    bad <- which(!is.finite(inj))
    stop("Non-finite Injection_order in batch ", batch,
         ": samples ", paste(cd$Sample_ID[bad], collapse = ", "),
         " (values: ", paste(inj[bad], collapse = ", "), ")")
  }

  n_skipped_n   <- 0L
  n_skipped_err <- 0L

  for (i in seq_len(nrow(mat))) {
    y_s <- as.numeric(mat[i, samp_idx])
    x_s <- inj[samp_idx]
    ok  <- is.finite(y_s) & y_s > 0

    if (sum(ok) < min_obs) { n_skipped_n <- n_skipped_n + 1L; next }

    tryCatch({
      ok_inj       <- !is.na(inj)
      pred         <- rep(NA_real_, length(inj))
      pred[ok_inj] <- fit_predict_spline(x_s[ok], y_s[ok], inj[ok_inj], spar = spar)
      med_s        <- median(y_s[ok])
      ratio        <- pred / med_s
      ratio[is.na(ratio) | ratio <= 0] <- 1
      mat[i, ]     <- mat[i, ] / ratio
    }, error = function(e) { n_skipped_err <<- n_skipped_err + 1L })
  }

  n_corrected <- n_feat - n_skipped_n - n_skipped_err
  message("  Batch ", batch, ": drift-corrected ", n_corrected, "/", n_feat, " features",
          if (n_skipped_n   > 0) paste0(" | ", n_skipped_n,   " skipped (insufficient sample observations)") else "",
          if (n_skipped_err > 0) paste0(" | ", n_skipped_err, " skipped (spline fit error)") else "")

  assay(se_b, 1, withDimnames = FALSE) <- mat
  se_b
}


# Per-batch hybrid smoothing spline drift correction -- identical decision
# logic to loess_correct_batch_hybrid()/huber_correct_batch_hybrid() (QC-based
# if enough QC; else a samples-based trial validated against ltQC D-ratio;
# else uncorrected), substituting the spline pair above. See
# loess_correct_batch_hybrid()'s documentation for the full rationale, which
# applies unchanged here.
spline_correct_batch_hybrid <- function(se_b, qc_spar, sample_spar, sample_min_obs,
                                         min_qc_per_batch = 4, min_ltqc_validate = 3,
                                         validate = TRUE, qc_spar_grid = numeric(0)) {
  cd    <- colData(se_b)
  n_qc  <- sum(cd$QC == "QC")
  batch <- unique(cd$Batch)

  if (n_qc >= min_qc_per_batch) {
    if (length(qc_spar_grid) > 0) {
      message("  Batch ", batch, ": ", n_qc, " QC sample(s) (>= ", min_qc_per_batch,
              ") — using QC-based spline, spar selected per feature via CV")
      return(qc_cv_correct_batch(se_b, spline_spars = qc_spar_grid, min_obs = min_qc_per_batch))
    }
    message("  Batch ", batch, ": ", n_qc, " QC sample(s) (>= ", min_qc_per_batch,
            ") — using QC-based spline")
    return(spline_correct_batch(se_b, spar = qc_spar))
  }

  if (!validate) {
    message("  Batch ", batch, ": ", n_qc, " QC sample(s) (< ", min_qc_per_batch,
            ") — applying QC-free spline on samples unconditionally (validation disabled)")
    return(spline_correct_batch_samples(se_b, spar = sample_spar, min_obs = sample_min_obs))
  }

  n_ltqc <- sum(cd$QC == "ltQC")
  if (n_ltqc < min_ltqc_validate) {
    message("  Batch ", batch, ": ", n_qc, " QC sample(s) (< ", min_qc_per_batch,
            ") and only ", n_ltqc, " ltQC sample(s) (< ", min_ltqc_validate,
            ") to validate a QC-free fit — leaving batch uncorrected")
    return(se_b)
  }

  message("  Batch ", batch, ": ", n_qc, " QC sample(s) (< ", min_qc_per_batch,
          ") — trialing QC-free spline on samples, to be validated against ",
          n_ltqc, " ltQC sample(s)")
  se_trial <- spline_correct_batch_samples(se_b, spar = sample_spar, min_obs = sample_min_obs)

  dratio_before <- eval_ltqc_dratio(se_b)
  dratio_after  <- eval_ltqc_dratio(se_trial)

  if (is.na(dratio_before) || is.na(dratio_after)) {
    message("  Batch ", batch, ": could not compute ltQC/Sample D-ratio before/after (NA) — ",
            "leaving batch uncorrected")
    return(se_b)
  }

  if (dratio_after < dratio_before) {
    message("  Batch ", batch, ": ltQC/Sample D-ratio improved with QC-free correction (",
            round(dratio_before, 4), " -> ", round(dratio_after, 4), ") — keeping correction")
    se_trial
  } else {
    message("  Batch ", batch, ": ltQC/Sample D-ratio did not improve with QC-free correction (",
            round(dratio_before, 4), " -> ", round(dratio_after, 4), ") — reverting to uncorrected")
    se_b
  }
}


# Applies correct_fn(se_b) only if the batch has >= min_qc_per_batch QC
# samples; otherwise leaves se_b unchanged and messages why. Gives
# basis="qc" the same batch-level gate basis="hybrid" already has via
# loess_correct_batch_hybrid()/huber_correct_batch_hybrid() -- closes a real
# inconsistency: loess_correct_batch()/huber_correct_batch() on their own
# only ever checked per-FEATURE QC counts (>=4), never a batch-level minimum,
# so without this gate a batch with e.g. 1-3 QC samples would silently get
# partial per-feature correction instead of being left uncorrected like the
# hybrid basis would do for the same batch.
gate_by_qc_count <- function(se_b, min_qc_per_batch, correct_fn) {
  n_qc  <- sum(colData(se_b)$QC == "QC")
  batch <- unique(colData(se_b)$Batch)
  if (n_qc < min_qc_per_batch) {
    message("  Batch ", batch, ": ", n_qc, " QC sample(s) (< ", min_qc_per_batch,
            ") — leaving uncorrected (basis=qc requires the batch-level minimum)")
    return(se_b)
  }
  correct_fn(se_b)
}

# Dispatch table for the drift x basis axes (see R/method_spec.R): returns a
# function(se_b, params) implementing drift_method paired with basis, for
# every combination except "auto" (auto operates over all batches at once via
# auto_select_drift_correction(), not per-batch, so it's dispatched directly
# by run_correction() in R/correction_methods.R rather than through here).
# `params` is the single shared list of resolved config values built once in
# notame-workflow.r; only the fields relevant to the resolved function are
# read from it.
resolve_drift_fn <- function(drift_method, basis) {
  key <- paste(drift_method, basis, sep = ":")
  switch(key,
    "loess:qc" = function(se_b, p) gate_by_qc_count(se_b, p$min_qc_per_batch, function(x) {
      if (length(p$loess_qc_cv_spans) > 0)
        qc_cv_correct_batch(x, loess_spans = p$loess_qc_cv_spans, min_obs = p$min_qc_per_batch)
      else
        loess_correct_batch(x, span = p$loess_qc_span)
    }),
    "loess:samples" = function(se_b, p)
      loess_correct_batch_samples(se_b, span = p$loess_sample_span, min_obs = p$drift_sample_min_obs),
    "loess:hybrid" = function(se_b, p)
      loess_correct_batch_hybrid(se_b, qc_span = p$loess_qc_span, sample_span = p$loess_sample_span,
                                  sample_min_obs = p$drift_sample_min_obs,
                                  min_qc_per_batch = p$min_qc_per_batch,
                                  min_ltqc_validate = p$min_ltqc_validate,
                                  validate = p$drift_hybrid_validate,
                                  qc_span_grid = p$loess_qc_cv_spans),
    "huber:qc" = function(se_b, p) gate_by_qc_count(se_b, p$min_qc_per_batch, function(x) {
      if (length(p$huber_qc_cv_ks) > 0)
        qc_cv_correct_batch(x, huber_ks = p$huber_qc_cv_ks, min_obs = p$min_qc_per_batch)
      else
        huber_correct_batch(x, k = p$huber_qc_k)
    }),
    "huber:samples" = function(se_b, p)
      huber_correct_batch_samples(se_b, k = p$huber_sample_k, min_obs = p$drift_sample_min_obs),
    "huber:hybrid" = function(se_b, p)
      huber_correct_batch_hybrid(se_b, qc_k = p$huber_qc_k, sample_k = p$huber_sample_k,
                                  sample_min_obs = p$drift_sample_min_obs,
                                  min_qc_per_batch = p$min_qc_per_batch,
                                  min_ltqc_validate = p$min_ltqc_validate,
                                  validate = p$drift_hybrid_validate,
                                  qc_k_grid = p$huber_qc_cv_ks),
    "spline:qc" = function(se_b, p) gate_by_qc_count(se_b, p$min_qc_per_batch, function(x) {
      if (length(p$spline_qc_cv_spars) > 0)
        qc_cv_correct_batch(x, spline_spars = p$spline_qc_cv_spars, min_obs = p$min_qc_per_batch)
      else
        spline_correct_batch(x, spar = p$spline_qc_spar)
    }),
    "spline:samples" = function(se_b, p)
      spline_correct_batch_samples(se_b, spar = p$spline_sample_spar, min_obs = p$drift_sample_min_obs),
    "spline:hybrid" = function(se_b, p)
      spline_correct_batch_hybrid(se_b, qc_spar = p$spline_qc_spar, sample_spar = p$spline_sample_spar,
                                   sample_min_obs = p$drift_sample_min_obs,
                                   min_qc_per_batch = p$min_qc_per_batch,
                                   min_ltqc_validate = p$min_ltqc_validate,
                                   validate = p$drift_hybrid_validate,
                                   qc_spar_grid = p$spline_qc_cv_spars),
    "notame_spline:qc" = function(se_b, p) process_batch(se_b),
    stop("No drift-correction dispatch for drift_method='", drift_method, "', basis='", basis, "'")
  )
}


# ─────────────────────────────────────────────────────────────────────────────
# Auto-selected drift correction.
#
# Rather than assuming LOESS is always the right fit, this evaluates several
# candidate (method, parameter) combinations against held-out data and picks
# whichever actually performs best for this specific batch:
#   - LOESS at each span in loess_spans (flexible, local; can overfit noise at
#     small QC counts)
#   - Huber robust regression (MASS::rlm) at each k in huber_ks (a single
#     robust linear trend; more stable than LOESS at small n, but can't
#     capture curved drift)
#   - a flat/constant baseline (predicts the training median everywhere,
#     i.e. applies no correction at all) -- this is what lets a batch with no
#     real fittable drift land on "don't correct" even when QC is plentiful,
#     without a separate code path: a ratio of predicted/training-median is
#     always exactly 1 for this candidate.
#
# QC available (>= min_qc_per_batch): each candidate is scored by leave-one-
# out CV directly on QC (relative squared error, aggregated across features
# by median -- the same "aggregate across features" reasoning as
# loess_correct_batch_hybrid's D-ratio check: with only a handful of QC per
# batch, per-feature comparisons are noise-dominated). The winner is then fit
# once more on all QC and applied to the batch.
#
# QC insufficient, ltQC available (>= min_ltqc_validate): each candidate is
# fit on biological samples and applied to the whole batch, then scored by
# the resulting ltQC/Sample D-ratio (eval_ltqc_dratio, from R/qc_metrics.R) --
# no CV needed here since ltQC is never used to fit anything, so the entire
# ltQC set is already a genuine hold-out. The flat candidate participates in
# this comparison too, so "no candidate beats leaving it uncorrected" is
# handled by the same mechanism rather than a special case.
#
# Neither available: leave uncorrected, same reasoning as
# loess_correct_batch_hybrid -- no data to justify or validate a correction.
# ─────────────────────────────────────────────────────────────────────────────

# Both fit in log2 space, not raw abundance, then exponentiate the prediction
# back -- external contract (raw y_train in, raw-scale prediction out) is
# unchanged, so callers need no changes. This matters much more for Huber
# than LOESS: real LC-MS drift is often multiplicative/exponential-shaped
# (y = y0 * f(t)), which is exactly linear in log space (log(y) = log(y0) +
# log(f(t))) but badly mismatched by a single straight line fit on the raw
# scale -- the same reason ComBat/limma always operate in log space. LOESS's
# local flexibility partly compensates for this on the raw scale already,
# but fitting log2(y) keeps both candidates on equal, correctly-specified
# footing rather than handicapping Huber's rigid linear form. Callers already
# filter out non-positive y before this is reached (see raw_feature_cv_scores
# and apply_drift_candidate_to_batch's `ok` filters), so log2() here is safe.
fit_predict_loess <- function(x_train, y_train, x_new, span) {
  fit <- suppressWarnings(
    loess(y ~ x, data = data.frame(x = x_train, y = log2(y_train)), span = span)
  )
  pred_log <- as.numeric(suppressWarnings(predict(fit, newdata = data.frame(x = x_new))))
  2^pred_log
}

fit_predict_huber <- function(x_train, y_train, x_new, k) {
  fit <- suppressWarnings(
    MASS::rlm(y ~ x, data = data.frame(x = x_train, y = log2(y_train)),
              psi = MASS::psi.huber, k = k, maxit = 100)
  )
  pred_log <- as.numeric(suppressWarnings(predict(fit, newdata = data.frame(x = x_new))))
  2^pred_log
}

fit_predict_flat <- function(x_train, y_train, x_new) {
  rep(median(y_train), length(x_new))
}

# Cubic smoothing spline (stats::smooth.spline), fit in log2 space for the
# same reason fit_predict_loess()/fit_predict_huber() do. `spar` is a fixed
# value here, not the cv=TRUE auto-selection pmp's QCRSC/smooth.spline itself
# offer -- deliberately: this pipeline's own CV/held-out scoring is already
# the outer layer that picks among candidates, so an inner CV inside the
# candidate's own fit would be a redundant, harder-to-reason-about nested CV
# rather than an improvement (same convention as loess span / huber k, which
# are also fixed grid points, not self-tuning).
# Unlike LOESS (family="symmetric" available) or Huber (psi.huber), plain
# smooth.spline has no robust/outlier-resistant fitting mode -- a single wild
# training point can pull it more than the other two candidate types. This
# matters most for the samples-based fit (biological samples are noisier and
# less controlled than QC), which is part of why the samples-tier default
# spar is set higher (smoother, less locally reactive) than the QC-tier one.
fit_predict_spline <- function(x_train, y_train, x_new, spar) {
  fit <- suppressWarnings(smooth.spline(x_train, log2(y_train), spar = spar))
  pred_log <- as.numeric(suppressWarnings(predict(fit, x_new)$y))
  2^pred_log
}

# Builds the candidate list: LOESS at each span, Huber at each k, smoothing
# spline at each spar, plus an optional flat baseline. Uses local() (not a
# plain for loop) so each candidate's closure captures its own
# span/k/spar value rather than the loop variable's final value.
build_drift_candidates <- function(loess_spans, huber_ks, spline_spars = numeric(0), include_flat = TRUE) {
  candidates <- list()
  for (s in loess_spans) {
    candidates[[length(candidates) + 1]] <- local({
      span_val <- s
      list(name = sprintf("loess(span=%.2g)", span_val),
           predict_fn = function(x_train, y_train, x_new)
             fit_predict_loess(x_train, y_train, x_new, span = span_val))
    })
  }
  for (k in huber_ks) {
    candidates[[length(candidates) + 1]] <- local({
      k_val <- k
      list(name = sprintf("huber(k=%.3g)", k_val),
           predict_fn = function(x_train, y_train, x_new)
             fit_predict_huber(x_train, y_train, x_new, k = k_val))
    })
  }
  for (sp in spline_spars) {
    candidates[[length(candidates) + 1]] <- local({
      spar_val <- sp
      list(name = sprintf("spline(spar=%.2g)", spar_val),
           predict_fn = function(x_train, y_train, x_new)
             fit_predict_spline(x_train, y_train, x_new, spar = spar_val))
    })
  }
  if (include_flat)
    candidates[[length(candidates) + 1]] <- list(name = "flat", predict_fn = fit_predict_flat)
  candidates
}

# Extracts each candidate's family ("loess", "huber", "spline", or "flat")
# from its name, as produced by build_drift_candidates() -- used to group
# per-parameter candidates by family for rank_select_family() below.
candidate_family <- function(candidates) {
  vapply(candidates, function(cand) {
    if (grepl("^loess", cand$name)) "loess"
    else if (grepl("^huber", cand$name)) "huber"
    else if (grepl("^spline", cand$name)) "spline"
    else "flat"
  }, character(1))
}

# Leave-one-out CV score for one feature, one candidate: mean relative
# squared error across held-out points ((y_i - pred_i) / median(y_train))^2.
# Relative, not raw, error -- so scores are comparable across features with
# very different abundance scales when aggregated later. Returns NA if the
# candidate couldn't be evaluated for this feature at all (too few points,
# every held-out fit failed).
cv_score_feature <- function(x, y, predict_fn) {
  n <- length(y)
  errs <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    x_train <- x[-i]; y_train <- y[-i]
    pred <- tryCatch(predict_fn(x_train, y_train, x[i]), error = function(e) NA_real_)
    if (length(pred) != 1 || !is.finite(pred)) next
    denom <- median(y_train)
    if (!is.finite(denom) || denom == 0) next
    errs[i] <- ((y[i] - pred) / denom)^2
  }
  if (all(is.na(errs))) return(NA_real_)
  mean(errs, na.rm = TRUE)
}

# Raw per-feature LOO-CV scores for every candidate on train_idx columns --
# NOT aggregated here, so callers can pool rows across multiple batches
# before aggregating (see auto_select_drift_correction()). Returns a
# features x candidates matrix. This is by far the most expensive part of
# auto_combat (features x candidates x held-out points, easily 100k+ small
# fits on real data), so progress is reported periodically -- a handful of
# messages at fixed intervals, not a redrawing terminal progress bar, since
# this runs inside Docker where output is a flat log rather than a live TTY
# (a \r-redrawing bar would show up as one log line per update instead of
# overwriting in place). progress_label, if given, prefixes each update
# (e.g. the batch name) so multi-batch runs stay distinguishable in the log.
raw_feature_cv_scores <- function(mat, train_idx, inj, candidates, min_obs = 4,
                                   progress_label = NULL) {
  x_tr     <- inj[train_idx]
  n_feat   <- nrow(mat)
  scores   <- matrix(NA_real_, n_feat, length(candidates))
  progress_every <- max(1L, round(n_feat / 10))
  prefix <- if (is.null(progress_label)) "  " else paste0("  ", progress_label, ": ")
  for (i in seq_len(n_feat)) {
    y_tr <- as.numeric(mat[i, train_idx])
    # y_tr > 0: excludes exact-zero abundances (msdial_to_notame()/
    # xcms_to_notame() clamp negative values to 0, not NA) so log2() in
    # fit_predict_loess()/fit_predict_huber() never sees a non-positive input.
    ok <- is.finite(y_tr) & is.finite(x_tr) & y_tr > 0
    if (sum(ok) >= min_obs) {
      xf <- x_tr[ok]; yf <- y_tr[ok]
      for (ci in seq_along(candidates))
        scores[i, ci] <- cv_score_feature(xf, yf, candidates[[ci]]$predict_fn)
    }
    if (i %% progress_every == 0 || i == n_feat)
      message(prefix, "LOO-CV progress: ", round(100 * i / n_feat), "% (", i, "/", n_feat, " features)")
  }
  colnames(scores) <- vapply(candidates, `[[`, character(1), "name")
  scores
}

# Held-out prediction score for one feature, one candidate: fit once on
# (x_train, y_train), predict at x_test, mean relative squared error against
# the actually-observed y_test -- ((y_test_i - pred_i) / median(y_train))^2,
# averaged over all test points. Not leave-one-out: train and test are a
# single fixed split (e.g. Sample vs ltQC), since test points here are
# already genuinely held out (never in the training set) -- no need to
# iteratively remove them one at a time the way cv_score_feature() does
# within a single group. This means it has no minimum-size requirement on
# the test side beyond >=1 point (unlike LOO-CV, which needs enough points
# left over after removing one to still fit a candidate) -- deliberately, so
# it stays usable even when ltQC is as small as 2 per batch. Returns NA if
# the candidate couldn't be evaluated at all.
holdout_score_feature <- function(x_train, y_train, x_test, y_test, predict_fn) {
  pred <- tryCatch(predict_fn(x_train, y_train, x_test), error = function(e) NULL)
  if (is.null(pred) || length(pred) != length(x_test) || any(!is.finite(pred))) return(NA_real_)
  denom <- median(y_train)
  if (!is.finite(denom) || denom == 0) return(NA_real_)
  mean(((y_test - pred) / denom)^2)
}

# Raw per-feature held-out scores for every candidate, fit on train_idx and
# evaluated on test_idx -- same features x candidates matrix shape as
# raw_feature_cv_scores(), same progress-reporting convention. min_obs
# applies only to the training side (need enough points to fit a candidate
# at all); the test side only needs >=1 finite, positive observation.
raw_feature_holdout_scores <- function(mat, train_idx, test_idx, inj, candidates, min_obs = 4,
                                        progress_label = NULL) {
  x_tr   <- inj[train_idx]
  x_te   <- inj[test_idx]
  n_feat <- nrow(mat)
  scores <- matrix(NA_real_, n_feat, length(candidates))
  progress_every <- max(1L, round(n_feat / 10))
  prefix <- if (is.null(progress_label)) "  " else paste0("  ", progress_label, ": ")
  for (i in seq_len(n_feat)) {
    y_tr <- as.numeric(mat[i, train_idx])
    y_te <- as.numeric(mat[i, test_idx])
    ok_tr <- is.finite(y_tr) & is.finite(x_tr) & y_tr > 0
    ok_te <- is.finite(y_te) & is.finite(x_te) & y_te > 0
    if (sum(ok_tr) >= min_obs && sum(ok_te) >= 1) {
      for (ci in seq_along(candidates))
        scores[i, ci] <- holdout_score_feature(x_tr[ok_tr], y_tr[ok_tr], x_te[ok_te], y_te[ok_te],
                                                candidates[[ci]]$predict_fn)
    }
    if (i %% progress_every == 0 || i == n_feat)
      message(prefix, "held-out ltQC progress: ", round(100 * i / n_feat), "% (", i, "/", n_feat, " features)")
  }
  colnames(scores) <- vapply(candidates, `[[`, character(1), "name")
  scores
}

# Fits predict_fn on train_idx columns (per feature) and applies the
# resulting correction ratio to every column -- same normalisation
# convention as loess_correct_batch()/loess_correct_batch_samples()
# (ratio = predicted / training median; a flat/no-op candidate always gives
# ratio == 1 everywhere). Features with fewer than min_obs finite training
# observations, or where the fit errors, are left uncorrected. quiet = TRUE
# suppresses the skip-count message, for trial applications during candidate
# scoring (where printing one line per batch x candidate would be noisy) --
# left FALSE for the final, actually-applied correction.
apply_drift_candidate_to_batch <- function(mat, train_idx, inj, predict_fn, min_obs = 4,
                                            quiet = FALSE) {
  n_feat    <- nrow(mat)
  n_skipped <- 0L
  for (i in seq_len(n_feat)) {
    y_tr <- as.numeric(mat[i, train_idx])
    x_tr <- inj[train_idx]
    # y_tr > 0: see raw_feature_cv_scores()'s comment on the same filter.
    ok   <- is.finite(y_tr) & is.finite(x_tr) & y_tr > 0
    if (sum(ok) < min_obs) { n_skipped <- n_skipped + 1L; next }
    tryCatch({
      ok_inj       <- !is.na(inj)
      pred         <- rep(NA_real_, length(inj))
      pred[ok_inj] <- predict_fn(x_tr[ok], y_tr[ok], inj[ok_inj])
      med_tr       <- median(y_tr[ok])
      ratio        <- pred / med_tr
      ratio[is.na(ratio) | ratio <= 0] <- 1
      mat[i, ]     <- mat[i, ] / ratio
    }, error = function(e) { n_skipped <<- n_skipped + 1L })
  }
  if (!quiet && n_skipped > 0)
    message("    (", n_skipped, "/", n_feat, " feature(s) skipped: insufficient observations or fit error)")
  mat
}

# Per-feature QC-based drift correction: instead of one span/k shared across
# every feature, evaluates a candidate grid (LOESS spans, Huber k's, or both
# -- whichever of loess_spans/huber_ks is non-empty) via leave-one-out CV on
# QC independently for EACH feature, and uses that feature's own best-scoring
# candidate to fit and apply the correction. Mirrors how notame::correct_drift()
# 's smooth.spline() step auto-selects its own smoothing parameter per feature
# via CV rather than sharing one value dataset-wide.
#
# Deliberately different from auto_select_drift_correction(), which pools
# evidence across features WITHIN a batch to pick one shared winner for that
# batch (batches themselves are never pooled together, each selects
# independently) -- this is intentionally per-feature instead. That's safe
# specifically because it only ever runs on QC: QC points are pure technical
# replicates, so a flexible per-feature fit has nothing biological to
# overfit, unlike a samples-based fit would (see loess_correct_batch_samples()
# 's documentation) -- which is why this isn't offered for the QC-free tier.
#
# Reuses raw_feature_cv_scores() (the same features x candidates LOO-CV
# machinery auto_select_drift_correction() uses) but picks a winner per row
# (per feature) instead of pooling scores across rows into one dataset-wide
# choice.
qc_cv_correct_batch <- function(se_b, loess_spans = numeric(0), huber_ks = numeric(0),
                                 spline_spars = numeric(0), min_obs = 4) {
  mat    <- assay(se_b, 1)
  cd     <- colData(se_b)
  qc_idx <- which(cd$QC == "QC")
  inj    <- as.numeric(cd$Injection_order)
  batch  <- unique(cd$Batch)
  n_feat <- nrow(mat)

  if (any(!is.finite(inj))) {
    bad <- which(!is.finite(inj))
    stop("Non-finite Injection_order in batch ", batch,
         ": samples ", paste(cd$Sample_ID[bad], collapse = ", "),
         " (values: ", paste(inj[bad], collapse = ", "), ")")
  }

  candidates <- build_drift_candidates(loess_spans, huber_ks, spline_spars, include_flat = FALSE)
  cand_names <- vapply(candidates, `[[`, character(1), "name")

  message("  Batch ", batch, ": ", length(qc_idx), " QC sample(s), ", n_feat,
          " features -- per-feature CV selection over: ", paste(cand_names, collapse = ", "))

  scores <- raw_feature_cv_scores(mat, qc_idx, inj, candidates, min_obs = min_obs,
                                   progress_label = paste0("Batch ", batch))

  n_skipped_qc  <- 0L
  n_skipped_err <- 0L
  winner_counts <- integer(length(candidates))

  for (i in seq_len(nrow(mat))) {
    if (all(is.na(scores[i, ]))) { n_skipped_qc <- n_skipped_qc + 1L; next }

    winner <- which.min(scores[i, ])

    y_qc <- as.numeric(mat[i, qc_idx])
    x_qc <- inj[qc_idx]
    ok   <- is.finite(y_qc) & y_qc > 0

    tryCatch({
      ok_inj       <- !is.na(inj)
      pred         <- rep(NA_real_, length(inj))
      pred[ok_inj] <- candidates[[winner]]$predict_fn(x_qc[ok], y_qc[ok], inj[ok_inj])
      med_qc       <- median(y_qc[ok])
      ratio        <- pred / med_qc
      ratio[is.na(ratio) | ratio <= 0] <- 1
      mat[i, ]     <- mat[i, ] / ratio
      winner_counts[winner] <- winner_counts[winner] + 1L
    }, error = function(e) { n_skipped_err <<- n_skipped_err + 1L })
  }

  n_corrected <- n_feat - n_skipped_qc - n_skipped_err
  message("  Batch ", batch, ": drift-corrected ", n_corrected, "/", n_feat, " features",
          if (n_skipped_qc  > 0) paste0(" | ", n_skipped_qc,  " skipped (insufficient QC observations)") else "",
          if (n_skipped_err > 0) paste0(" | ", n_skipped_err, " skipped (fit error)") else "")
  message("  Batch ", batch, ": per-feature candidate selection: ",
          paste(cand_names, winner_counts, sep = "=", collapse = ", "))

  assay(se_b, 1, withDimnames = FALSE) <- mat
  se_b
}

# Classifies each batch into a tier by QC/ltQC coverage, same thresholds as
# loess_correct_batch_hybrid(): "qc" (>= min_qc_per_batch QC), "sample"
# (insufficient QC but >= min_ltqc_validate ltQC), or "none".
classify_drift_tier <- function(se_b, min_qc_per_batch, min_ltqc_validate) {
  cd <- colData(se_b)
  if (sum(cd$QC == "QC") >= min_qc_per_batch) return("qc")
  if (sum(cd$QC == "ltQC") >= min_ltqc_validate) return("sample")
  "none"
}

# Picks ONE family (loess/huber/spline) to restrict the per-batch candidate
# search to, pooling evidence across every eligible batch via RANK, not raw
# score magnitude -- specifically to avoid two dominance failure modes a
# magnitude-based pool would have: a high-noise batch's inflated absolute
# scores swamping a clean batch's more reliable but smaller-magnitude signal,
# and a batch with more evaluable features implicitly getting more "votes"
# than a smaller/noisier one just by contributing more raw values. Per batch,
# each family is reduced to its own best-achievable score there (min across
# that family's own parameter grid, using the same per-batch
# median-across-features aggregate the unrestricted per-batch selection
# already uses); those per-batch family-bests are then ranked 1st/2nd/3rd
# within that batch, and ranks (not the underlying scores) are averaged
# across batches -- so every batch contributes exactly one vote, regardless
# of its own noise level or feature count, same principle as a Friedman-style
# rank-sum comparison across blocks. "flat" is excluded -- it isn't a family
# to lock the search to, it's a baseline that stays available regardless of
# which family wins (see select_qc_candidate_per_batch()/
# select_sample_candidate_per_batch()'s use of this function's result).
#
# agg_by_batch: list of per-candidate aggregate-score vectors, one per
# eligible batch, all using the same candidates/order. families:
# candidate_family(candidates) for that same candidate list. Returns the
# winning family name, or NA if no batch had usable evidence for more than
# one family (nothing to lock to -- caller should fall back to unrestricted).
rank_select_family <- function(agg_by_batch, families) {
  fam_levels <- setdiff(unique(families), "flat")
  if (length(fam_levels) <= 1) return(if (length(fam_levels) == 1) fam_levels else NA_character_)

  rank_sums <- setNames(numeric(length(fam_levels)), fam_levels)
  n_votes   <- setNames(integer(length(fam_levels)), fam_levels)

  for (agg in agg_by_batch) {
    fam_best <- vapply(fam_levels, function(f) {
      v <- agg[families == f]
      if (all(is.na(v))) NA_real_ else min(v, na.rm = TRUE)
    }, numeric(1))
    if (all(is.na(fam_best))) next
    fam_rank <- rank(fam_best, ties.method = "average", na.last = "keep")
    for (f in fam_levels) {
      if (!is.na(fam_rank[f])) {
        rank_sums[f] <- rank_sums[f] + fam_rank[f]
        n_votes[f]   <- n_votes[f] + 1L
      }
    }
  }
  if (all(n_votes == 0)) return(NA_character_)
  avg_rank <- ifelse(n_votes > 0, rank_sums / n_votes, Inf)
  names(avg_rank)[which.min(avg_rank)]
}

# For EACH batch in `eligible_idx` (indices into `batches`) independently:
# evaluates every candidate via per-feature LOO-CV on that batch's own QC
# data, aggregates by median, and picks that batch's own best-scoring
# candidate. No evidence is borrowed from other batches -- deliberately, per
# the same reasoning as loess_correct_batch_hybrid()'s per-batch tier
# decision: each batch's own correction should draw only on its own
# evidence, not a shared pooled estimate that could mask a batch behaving
# differently from the rest.
# Returns list(winners=, log=): winners is a list of winners, one per entry
# in eligible_idx (NULL for a batch where no candidate could be evaluated --
# too few QC observations per feature); log is a data.frame with one row per
# (batch, candidate) -- batch, tier, candidate, score_method, prediction_score,
# ltqc_dratio, selected -- for the caller to save alongside the corrected
# output (see notame-workflow.r's per-method output writing), so the
# reasoning behind each batch's choice survives past the console log.
# Also prints, per batch and candidate, that batch's ltQC/Sample D-ratio --
# informational only, does not affect selection: LOO-CV on QC is already a
# legitimate, non-circular criterion, so D-ratio here is only a sanity-check
# cross-reference (same reasoning as ltqc_permanova_* alongside
# qc_permanova_* in qc_metrics.R -- if LOO-CV improves but D-ratio disagrees,
# that disagreement is itself worth noticing, not something to average away).
select_qc_candidate_per_batch <- function(batches, batch_names, eligible_idx, candidates, min_cv_obs) {
  if (length(eligible_idx) == 0) return(list(winners = list(), log = NULL))

  message("==> QC-based batches (", length(eligible_idx), "): evaluating ",
          length(candidates), " drift-correction candidate(s) via LOO-CV on QC, per batch")

  families <- candidate_family(candidates)

  # --- Pass 1: compute each batch's scores once (no winner decided yet) ---
  batch_data <- vector("list", length(eligible_idx))
  for (row in seq_along(eligible_idx)) {
    bi     <- eligible_idx[row]
    se_b   <- batches[[bi]]
    qc_idx <- which(colData(se_b)$QC == "QC")
    inj    <- as.numeric(colData(se_b)$Injection_order)

    cv_scores <- raw_feature_cv_scores(assay(se_b, 1), qc_idx, inj, candidates,
                                        min_obs = min_cv_obs,
                                        progress_label = paste0("Batch ", batch_names[bi]))
    agg <- apply(cv_scores, 2, median, na.rm = TRUE)

    qc_dratios <- vapply(candidates, function(cand) {
      mat_trial <- apply_drift_candidate_to_batch(assay(se_b, 1), qc_idx, inj,
                                                    cand$predict_fn, min_obs = min_cv_obs,
                                                    quiet = TRUE)
      se_trial <- se_b
      assay(se_trial, 1, withDimnames = FALSE) <- mat_trial
      eval_ltqc_dratio(se_trial)
    }, numeric(1))

    batch_data[[row]] <- list(bi = bi, agg = agg, dratios = qc_dratios)
  }

  # --- Family lock: rank-pooled across all eligible batches ---
  locked_family <- rank_select_family(lapply(batch_data, `[[`, "agg"), families)
  if (is.na(locked_family)) {
    message("==> QC-tier: no family could be ranked ",
            "(only one family in the grid, or no usable evidence) -- leaving search unrestricted")
  } else {
    message("==> QC-tier family lock: '", locked_family, "' (rank-sum across ",
            length(eligible_idx), " batch(es)) -- restricting per-batch search to this family + flat")
  }

  # --- Pass 2: pick each batch's own winner, restricted to the locked family if any ---
  winners  <- vector("list", length(eligible_idx))
  names(winners) <- batch_names[eligible_idx]
  log_rows <- list()

  for (row in seq_along(eligible_idx)) {
    bi      <- batch_data[[row]]$bi
    agg     <- batch_data[[row]]$agg
    qc_dratios <- batch_data[[row]]$dratios

    allowed_ci <- if (!is.na(locked_family)) which(families %in% c(locked_family, "flat")) else seq_along(candidates)
    agg_allowed <- agg; agg_allowed[-allowed_ci] <- NA_real_
    winner_ci <- if (all(is.na(agg_allowed))) NA_integer_ else allowed_ci[which.min(agg[allowed_ci])]

    for (ci in seq_along(candidates)) {
      message("    Batch ", batch_names[bi], ": ", candidates[[ci]]$name, ": LOO-CV score = ",
              if (is.na(agg[ci])) "NA" else signif(agg[ci], 4),
              ", ltQC/Sample D-ratio = ",
              if (is.na(qc_dratios[ci])) "NA (no ltQC available)" else round(qc_dratios[ci], 4))
      log_rows[[length(log_rows) + 1]] <- data.frame(
        batch = batch_names[bi], tier = "qc", candidate = candidates[[ci]]$name, family = families[ci],
        score_method = "loocv_qc", prediction_score = agg[ci], ltqc_dratio = qc_dratios[ci],
        selected = !is.na(winner_ci) && ci == winner_ci,
        stringsAsFactors = FALSE
      )
    }

    if (is.na(winner_ci)) {
      message("  Batch ", batch_names[bi], ": no candidate could be evaluated",
              " (too few QC observations per feature, or none in the locked family)")
      next
    }
    winner <- candidates[[winner_ci]]
    message("  Batch ", batch_names[bi], ": selected ", winner$name,
            " (LOO-CV score = ", signif(agg[winner_ci], 4), ")")
    winners[[row]] <- winner
  }
  list(winners = winners, log = do.call(rbind, log_rows))
}

# Same idea as select_qc_candidate_per_batch(), for the samples-based tier:
# for EACH batch in `eligible_idx` independently, fits every candidate on
# that batch's own Sample rows and scores it by held-out prediction accuracy
# against that batch's own ltQC points (never used in fitting) --
# holdout_score_feature()/raw_feature_holdout_scores(), a genuine train/test
# split rather than D-ratio. Picks that batch's own best-scoring candidate.
#
# D-ratio (MAD(ltQC)/MAD(Sample) after trial-applying the candidate) used to
# be this tier's selection criterion, but it's a ratio, and a candidate can
# improve a ratio by unilaterally shrinking either side of it -- specifically,
# by over-smoothing the Sample data it's fit on, compressing real biological
# variance along with any drift, which shrinks MAD(Sample) and makes D-ratio
# look better even though the correction may have destroyed signal rather
# than removed drift. The held-out prediction score doesn't have that
# vulnerability (it's an absolute prediction error, not a ratio, so
# compressing Sample's own dispersion doesn't mechanically improve it), so it
# now decides the winner; D-ratio is still computed and logged, but purely
# informational -- same role it already has for the QC tier.
#
# `eligible_idx` should already be restricted to batches with enough ltQC to
# evaluate at all (>= 1 -- unlike LOO-CV, held-out scoring has no minimum
# training-set-minus-one requirement on the ltQC side, so this works fine
# even with as few as 2 ltQC per batch; AUTO_MIN_LTQC_VALIDATE/
# DRIFT_MIN_LTQC_VALIDATE is a separate, configurable eligibility gate, not a
# requirement of this scoring method itself). Returns list(winners=, log=),
# same shape as select_qc_candidate_per_batch().
select_sample_candidate_per_batch <- function(batches, batch_names, eligible_idx, candidates, min_cv_obs) {
  if (length(eligible_idx) == 0) return(list(winners = list(), log = NULL))

  message("==> QC-free batches (", length(eligible_idx), "): evaluating ",
          length(candidates), " candidate(s) (fit on samples) via held-out ltQC prediction",
          " score, per batch")

  families <- candidate_family(candidates)

  # --- Pass 1: compute each batch's scores once (no winner decided yet) ---
  batch_data <- vector("list", length(eligible_idx))
  for (row in seq_along(eligible_idx)) {
    bi         <- eligible_idx[row]
    se_b       <- batches[[bi]]
    sample_idx <- which(colData(se_b)$QC == "Sample")
    ltqc_idx   <- which(colData(se_b)$QC == "ltQC")
    inj        <- as.numeric(colData(se_b)$Injection_order)

    holdout_scores <- raw_feature_holdout_scores(assay(se_b, 1), sample_idx, ltqc_idx, inj, candidates,
                                                  min_obs = min_cv_obs,
                                                  progress_label = paste0("Batch ", batch_names[bi]))
    agg <- apply(holdout_scores, 2, median, na.rm = TRUE)

    dratios <- vapply(candidates, function(cand) {
      mat_trial <- apply_drift_candidate_to_batch(assay(se_b, 1), sample_idx, inj,
                                                    cand$predict_fn, min_obs = min_cv_obs, quiet = TRUE)
      se_trial <- se_b
      assay(se_trial, 1, withDimnames = FALSE) <- mat_trial
      eval_ltqc_dratio(se_trial)
    }, numeric(1))

    batch_data[[row]] <- list(bi = bi, agg = agg, dratios = dratios)
  }

  # --- Family lock: rank-pooled across all eligible batches ---
  locked_family <- rank_select_family(lapply(batch_data, `[[`, "agg"), families)
  if (is.na(locked_family)) {
    message("==> Samples-tier: no family could be ranked ",
            "(only one family in the grid, or no usable evidence) -- leaving search unrestricted")
  } else {
    message("==> Samples-tier family lock: '", locked_family, "' (rank-sum across ",
            length(eligible_idx), " batch(es)) -- restricting per-batch search to this family + flat")
  }

  # --- Pass 2: pick each batch's own winner, restricted to the locked family if any ---
  winners  <- vector("list", length(eligible_idx))
  names(winners) <- batch_names[eligible_idx]
  log_rows <- list()

  for (row in seq_along(eligible_idx)) {
    bi      <- batch_data[[row]]$bi
    agg     <- batch_data[[row]]$agg
    dratios <- batch_data[[row]]$dratios

    allowed_ci <- if (!is.na(locked_family)) which(families %in% c(locked_family, "flat")) else seq_along(candidates)
    agg_allowed <- agg; agg_allowed[-allowed_ci] <- NA_real_
    winner_ci <- if (all(is.na(agg_allowed))) NA_integer_ else allowed_ci[which.min(agg[allowed_ci])]

    for (ci in seq_along(candidates)) {
      message("    Batch ", batch_names[bi], ": ", candidates[[ci]]$name, ": held-out ltQC score = ",
              if (is.na(agg[ci])) "NA" else signif(agg[ci], 4),
              ", ltQC/Sample D-ratio = ",
              if (is.na(dratios[ci])) "NA" else round(dratios[ci], 4))
      log_rows[[length(log_rows) + 1]] <- data.frame(
        batch = batch_names[bi], tier = "sample", candidate = candidates[[ci]]$name, family = families[ci],
        score_method = "holdout_ltqc", prediction_score = agg[ci], ltqc_dratio = dratios[ci],
        selected = !is.na(winner_ci) && ci == winner_ci,
        stringsAsFactors = FALSE
      )
    }

    if (is.na(winner_ci)) {
      message("  Batch ", batch_names[bi], ": no candidate's held-out ltQC score could be computed",
              " (or none in the locked family)")
      next
    }
    winner <- candidates[[winner_ci]]
    message("  Batch ", batch_names[bi], ": selected ", winner$name,
            " (held-out ltQC score = ", signif(agg[winner_ci], 4), ")")
    winners[[row]] <- winner
  }
  list(winners = winners, log = do.call(rbind, log_rows))
}

# Two decisions are made at two different grains, deliberately -- not one
# pooled decision and not one fully per-batch decision:
#
#   1. WHICH FAMILY (loess/huber/spline) -- pooled across all eligible
#      batches in a tier, via rank_select_family()'s rank-sum (not raw score
#      magnitude, which would let a high-noise or feature-rich batch dominate
#      the pool). This is a coarse, dataset-level question -- "which curve
#      shape generally suits this study's drift" -- where pooling gives real
#      statistical power the alternative (deciding it from one batch's own
#      often-modest QC/ltQC count) wouldn't have. It also keeps every
#      batch's residual noise the same qualitative shape going into the
#      batch-correction step, not just matched magnitude (which ComBat's own
#      variance-equalisation step already handles regardless).
#   2. WHICH PARAMETER within that family, and WHETHER to apply anything at
#      all -- still decided independently per batch, from that batch's own
#      evidence only (LOO-CV on its own QC, or held-out ltQC prediction score
#      on its own samples), never borrowed from other batches. This is the
#      fine-grained, batch-level question -- "how much should THIS batch be
#      smoothed" -- where per-batch data genuinely differs (QC count, noise
#      level, how much drift is actually present), so forcing one shared
#      answer here would be the wrong kind of pooling. "flat" (no-op) stays
#      available to every batch regardless of which family won, so a batch
#      that doesn't need correction can still say so.
#
# This was a deliberate change from an earlier, fully pooled design (one
# shared winner per tier -- family AND parameter -- chosen from evidence
# combined across all eligible batches, then applied to every batch in that
# tier regardless of fit) and, before that, a fully per-batch design (no
# pooling at all, including of family) -- both were tried and superseded;
# see git history for that discussion. The two-grain split above is a
# considered middle point, not a default left in place from either.
#
# `basis` controls which batches are eligible to select (and receive) their
# own candidate:
#   "hybrid"  -- batches are tiered by QC/ltQC coverage (classify_drift_tier());
#                QC-tier batches each independently select from the
#                QC-anchored candidates, sample-tier batches each
#                independently select from the samples-based candidates. A
#                "none"-tier batch (insufficient QC AND ltQC) is left
#                uncorrected.
#   "qc"      -- every batch with >= min_qc_per_batch QC independently
#                selects (and receives) its own QC-anchored candidate; no
#                samples-based selection is attempted at all. Batches below
#                the threshold are left uncorrected, with no samples-based
#                fallback -- this is what makes basis="qc" strictly QC-only,
#                unlike "hybrid".
#   "samples" -- every batch with >= min_ltqc_validate ltQC independently
#                selects (and receives) its own samples-based candidate,
#                using its own held-out ltQC/Sample D-ratio. A batch below
#                that threshold has no held-out evidence to select from at
#                all and is left uncorrected -- unlike the old pooled
#                design, there is no shared winner from other batches to
#                fall back on.
auto_select_drift_correction <- function(data,
                                          basis               = c("hybrid", "qc", "samples"),
                                          loess_spans        = c(0.5, 0.75, 0.9),
                                          huber_ks            = c(1.0, 1.345, 2.0),
                                          spline_spars        = c(0.4, 0.6, 0.8),
                                          sample_loess_spans  = c(0.3, 0.6, 0.9),
                                          sample_huber_ks     = c(1.0, 1.345, 2.0),
                                          sample_spline_spars = c(0.5, 0.7, 0.9),
                                          min_qc_per_batch    = 4,
                                          min_ltqc_validate   = 3,
                                          min_cv_obs          = 4) {
  basis   <- match.arg(basis)
  batches <- split_by_batch(data)
  batch_names <- vapply(batches, function(se_b) as.character(unique(colData(se_b)$Batch)), character(1))
  n_qc   <- vapply(batches, function(se_b) sum(colData(se_b)$QC == "QC"), integer(1))
  n_ltqc <- vapply(batches, function(se_b) sum(colData(se_b)$QC == "ltQC"), integer(1))

  tiers <- NULL
  if (basis == "hybrid") {
    tiers <- vapply(batches, classify_drift_tier, character(1),
                     min_qc_per_batch = min_qc_per_batch, min_ltqc_validate = min_ltqc_validate)
    for (bi in seq_along(batches)) message("  Batch ", batch_names[bi], ": tier = ", tiers[bi])
    qc_eligible_idx     <- which(tiers == "qc")
    qc_apply_idx        <- qc_eligible_idx
    sample_eligible_idx <- which(tiers == "sample")
    sample_apply_idx    <- sample_eligible_idx
  } else if (basis == "qc") {
    qc_eligible_idx     <- which(n_qc >= min_qc_per_batch)
    qc_apply_idx        <- qc_eligible_idx
    sample_eligible_idx <- integer(0)
    sample_apply_idx    <- integer(0)
    for (bi in seq_along(batches))
      message("  Batch ", batch_names[bi], ": ", n_qc[bi], " QC sample(s) (",
              if (bi %in% qc_eligible_idx) paste0(">= ", min_qc_per_batch)
              else paste0("< ", min_qc_per_batch, ", leaving uncorrected"), ")")
  } else {  # basis == "samples"
    qc_eligible_idx     <- integer(0)
    qc_apply_idx        <- integer(0)
    sample_eligible_idx <- which(n_ltqc >= min_ltqc_validate)
    sample_apply_idx    <- sample_eligible_idx  # per-batch: only batches that selected for themselves
    for (bi in seq_along(batches))
      message("  Batch ", batch_names[bi], ": ", n_ltqc[bi], " ltQC sample(s) (",
              if (bi %in% sample_eligible_idx) paste0(">= ", min_ltqc_validate, ", selects its own candidate")
              else paste0("< ", min_ltqc_validate, ", left uncorrected"), ")")
  }

  qc_candidates <- build_drift_candidates(loess_spans, huber_ks, spline_spars, include_flat = TRUE)
  qc_result  <- select_qc_candidate_per_batch(batches, batch_names, qc_eligible_idx, qc_candidates, min_cv_obs)
  qc_winners <- qc_result$winners

  sample_candidates <- build_drift_candidates(sample_loess_spans, sample_huber_ks, sample_spline_spars,
                                               include_flat = TRUE)
  sample_result  <- select_sample_candidate_per_batch(batches, batch_names, sample_eligible_idx,
                                                        sample_candidates, min_cv_obs)
  sample_winners <- sample_result$winners

  # --- Apply each batch's own selected candidate ---
  message("==> Applying selected method(s) per batch")
  for (bi in seq_along(batches)) {
    se_b <- batches[[bi]]
    qc_winner     <- qc_winners[[match(bi, qc_eligible_idx)]]
    sample_winner <- sample_winners[[match(bi, sample_eligible_idx)]]
    if (bi %in% qc_apply_idx && !is.null(qc_winner)) {
      message("  Batch ", batch_names[bi], ": applying ", qc_winner$name, " (QC-based)")
      qc_idx <- which(colData(se_b)$QC == "QC")
      mat <- apply_drift_candidate_to_batch(assay(se_b, 1), qc_idx,
                                             as.numeric(colData(se_b)$Injection_order),
                                             qc_winner$predict_fn, min_obs = min_cv_obs)
      assay(se_b, 1, withDimnames = FALSE) <- mat
    } else if (bi %in% sample_apply_idx && !is.null(sample_winner)) {
      message("  Batch ", batch_names[bi], ": applying ", sample_winner$name,
              " (QC-free, fit on samples)")
      sample_idx <- which(colData(se_b)$QC == "Sample")
      mat <- apply_drift_candidate_to_batch(assay(se_b, 1), sample_idx,
                                             as.numeric(colData(se_b)$Injection_order),
                                             sample_winner$predict_fn, min_obs = min_cv_obs)
      assay(se_b, 1, withDimnames = FALSE) <- mat
    } else {
      reason <- if (basis == "hybrid") {
        if (tiers[bi] == "none") "insufficient QC and ltQC" else "no candidate selected"
      } else if (basis == "qc") {
        if (!(bi %in% qc_eligible_idx)) "insufficient QC" else "no candidate selected"
      } else {
        "no candidate selected"
      }
      message("  Batch ", batch_names[bi], ": leaving uncorrected (", reason, ")")
    }
    batches[[bi]] <- se_b
  }

  list(batches = batches, log = do.call(rbind, list(qc_result$log, sample_result$log)))
}

# ─────────────────────────────────────────────────────────────────────────────
# Within-batch drift correction helpers.
#
# split_by_batch()              — splits a SE into a list of per-batch SEs, sorted by injection order
# process_batch()                — notame cubic spline drift correction wrapper
# loess_correct_batch()          — QC-based LOESS drift correction
# loess_correct_batch_samples()  — QC-free LOESS drift correction (fit on biological samples)
# loess_correct_batch_hybrid()   — per-batch: QC-based if enough QC, else samples-based
#                                   trial validated against ltQC (kept only if it helps), else uncorrected

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
    ok   <- is.finite(y_qc)

    if (sum(ok) < 4) { n_skipped_qc <- n_skipped_qc + 1L; next }

    tryCatch({
      fit          <- loess(y ~ x, data = data.frame(x = x_qc[ok], y = y_qc[ok]), span = span)
      ok_inj       <- !is.na(inj)
      pred         <- rep(NA_real_, length(inj))
      pred[ok_inj] <- predict(fit, newdata = data.frame(x = inj[ok_inj]))
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
                                        validate = TRUE) {
  cd    <- colData(se_b)
  n_qc  <- sum(cd$QC == "QC")
  batch <- unique(cd$Batch)

  if (n_qc >= min_qc_per_batch) {
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
    ok  <- is.finite(y_s)

    if (sum(ok) < min_obs) { n_skipped_n <- n_skipped_n + 1L; next }

    tryCatch({
      fit          <- loess(y ~ x, data = data.frame(x = x_s[ok], y = y_s[ok]),
                            span = span, family = "symmetric")
      ok_inj       <- !is.na(inj)
      pred         <- rep(NA_real_, length(inj))
      pred[ok_inj] <- predict(fit, newdata = data.frame(x = inj[ok_inj]))
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

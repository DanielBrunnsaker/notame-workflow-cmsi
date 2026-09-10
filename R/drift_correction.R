# ─────────────────────────────────────────────────────────────────────────────
# Within-batch drift correction helpers.
#
# split_by_batch()              — splits a SE into a list of per-batch SEs, sorted by injection order
# process_batch()                — notame cubic spline drift correction wrapper
# loess_correct_batch()          — QC-based LOESS drift correction
# loess_correct_batch_samples()  — QC-free LOESS drift correction (fit on biological samples)
# loess_correct_batch_hybrid()   — per-batch: QC-based if enough QC, else samples-based
#                                   trial validated against ltQC (kept only if it helps), else uncorrected
# auto_select_drift_correction() — picks ONE drift-correction method for all QC-based batches and
#                                   ONE for all QC-free batches (never a different method per batch),
#                                   from several fitting candidates (LOESS at several spans, Huber
#                                   regression at several k, a flat/no-op baseline), using evidence
#                                   pooled across all batches in that tier: leave-one-out CV on QC
#                                   for QC-based batches, held-out ltQC/Sample D-ratio for QC-free
#                                   batches fit on samples. Returns a list of corrected per-batch SEs.

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

# Builds the candidate list: LOESS at each span, Huber at each k, plus an
# optional flat baseline. Uses local() (not a plain for loop) so each
# candidate's closure captures its own span/k value rather than the loop
# variable's final value.
build_drift_candidates <- function(loess_spans, huber_ks, include_flat = TRUE) {
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
  if (include_flat)
    candidates[[length(candidates) + 1]] <- list(name = "flat", predict_fn = fit_predict_flat)
  candidates
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
# features x candidates matrix.
raw_feature_cv_scores <- function(mat, train_idx, inj, candidates, min_obs = 4) {
  x_tr   <- inj[train_idx]
  n_feat <- nrow(mat)
  scores <- matrix(NA_real_, n_feat, length(candidates))
  for (i in seq_len(n_feat)) {
    y_tr <- as.numeric(mat[i, train_idx])
    # y_tr > 0: excludes exact-zero abundances (msdial_to_notame()/
    # xcms_to_notame() clamp negative values to 0, not NA) so log2() in
    # fit_predict_loess()/fit_predict_huber() never sees a non-positive input.
    ok   <- is.finite(y_tr) & is.finite(x_tr) & y_tr > 0
    if (sum(ok) < min_obs) next
    xf <- x_tr[ok]; yf <- y_tr[ok]
    for (ci in seq_along(candidates))
      scores[i, ci] <- cv_score_feature(xf, yf, candidates[[ci]]$predict_fn)
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

# Classifies each batch into a tier by QC/ltQC coverage, same thresholds as
# loess_correct_batch_hybrid(): "qc" (>= min_qc_per_batch QC), "sample"
# (insufficient QC but >= min_ltqc_validate ltQC), or "none".
classify_drift_tier <- function(se_b, min_qc_per_batch, min_ltqc_validate) {
  cd <- colData(se_b)
  if (sum(cd$QC == "QC") >= min_qc_per_batch) return("qc")
  if (sum(cd$QC == "ltQC") >= min_ltqc_validate) return("sample")
  "none"
}

# One drift-correction method is selected for the whole dataset, not per
# batch -- deliberately. Evidence is pooled ACROSS all batches in the same
# tier before picking a winner (all QC-tier batches' per-feature LOO-CV
# scores concatenated before aggregating; all sample-tier batches' D-ratios
# pooled before aggregating), so the decision draws on far more data than
# any single batch could offer, and every batch in a tier ends up using the
# same, consistent method -- no patchwork of different techniques across
# batches. QC-tier and sample-tier batches necessarily use methods from
# their own separate candidate pools (a QC-anchored fit can't be applied to
# a batch with no QC), but within each tier the choice is uniform.
auto_select_drift_correction <- function(data,
                                          loess_spans        = c(0.5, 0.75, 0.9),
                                          huber_ks            = c(1.0, 1.345, 2.0),
                                          sample_loess_spans  = c(0.3, 0.6, 0.9),
                                          sample_huber_ks     = c(1.0, 1.345, 2.0),
                                          min_qc_per_batch    = 4,
                                          min_ltqc_validate   = 3,
                                          min_cv_obs          = 4) {
  batches <- split_by_batch(data)
  tiers   <- vapply(batches, classify_drift_tier, character(1),
                     min_qc_per_batch = min_qc_per_batch, min_ltqc_validate = min_ltqc_validate)
  batch_names <- vapply(batches, function(se_b) as.character(unique(colData(se_b)$Batch)), character(1))
  for (bi in seq_along(batches))
    message("  Batch ", batch_names[bi], ": tier = ", tiers[bi])

  # --- QC-tier: one candidate, chosen from pooled leave-one-out CV on QC ---
  qc_candidates <- build_drift_candidates(loess_spans, huber_ks, include_flat = TRUE)
  qc_winner <- NULL
  qc_batches_idx <- which(tiers == "qc")
  if (length(qc_batches_idx) > 0) {
    message("==> QC-based batches (", length(qc_batches_idx), "): evaluating ",
            length(qc_candidates), " drift-correction candidate(s) via LOO-CV on QC, pooled across batches")
    pooled <- do.call(rbind, lapply(qc_batches_idx, function(bi) {
      se_b   <- batches[[bi]]
      qc_idx <- which(colData(se_b)$QC == "QC")
      raw_feature_cv_scores(assay(se_b, 1), qc_idx, as.numeric(colData(se_b)$Injection_order),
                             qc_candidates, min_obs = min_cv_obs)
    }))
    agg <- apply(pooled, 2, median, na.rm = TRUE)
    for (ci in seq_along(qc_candidates))
      message("    ", qc_candidates[[ci]]$name, ": pooled LOO-CV score = ",
              if (is.na(agg[ci])) "NA" else signif(agg[ci], 4))
    if (all(is.na(agg))) {
      message("  No candidate could be evaluated across QC-based batches",
              " (too few QC observations per feature)")
    } else {
      qc_winner <- qc_candidates[[which.min(agg)]]
      message("  Selected for all QC-based batches: ", qc_winner$name,
              " (pooled LOO-CV score = ", signif(min(agg, na.rm = TRUE), 4), ")")
    }
  }

  # --- Sample-tier: one candidate, chosen from pooled ltQC/Sample D-ratio ---
  sample_candidates <- build_drift_candidates(sample_loess_spans, sample_huber_ks, include_flat = TRUE)
  sample_winner <- NULL
  sample_batches_idx <- which(tiers == "sample")
  if (length(sample_batches_idx) > 0) {
    message("==> QC-free batches (", length(sample_batches_idx), "): evaluating ",
            length(sample_candidates), " candidate(s) (fit on samples) via held-out ltQC/Sample",
            " D-ratio, pooled across batches")
    dr_mat <- matrix(NA_real_, length(sample_batches_idx), length(sample_candidates))
    for (row in seq_along(sample_batches_idx)) {
      se_b       <- batches[[sample_batches_idx[row]]]
      sample_idx <- which(colData(se_b)$QC == "Sample")
      inj        <- as.numeric(colData(se_b)$Injection_order)
      for (ci in seq_along(sample_candidates)) {
        mat_trial <- apply_drift_candidate_to_batch(assay(se_b, 1), sample_idx, inj,
                                                      sample_candidates[[ci]]$predict_fn,
                                                      min_obs = min_cv_obs, quiet = TRUE)
        se_trial <- se_b
        assay(se_trial, 1, withDimnames = FALSE) <- mat_trial
        dr_mat[row, ci] <- eval_ltqc_dratio(se_trial)
      }
    }
    agg <- apply(dr_mat, 2, median, na.rm = TRUE)
    for (ci in seq_along(sample_candidates))
      message("    ", sample_candidates[[ci]]$name, ": pooled ltQC/Sample D-ratio = ",
              if (is.na(agg[ci])) "NA" else round(agg[ci], 4))
    if (all(is.na(agg))) {
      message("  No candidate's ltQC/Sample D-ratio could be computed across QC-free batches")
    } else {
      sample_winner <- sample_candidates[[which.min(agg)]]
      message("  Selected for all QC-free batches: ", sample_winner$name,
              " (pooled ltQC/Sample D-ratio = ", round(min(agg, na.rm = TRUE), 4), ")")
    }
  }

  # --- Apply the chosen winner(s) to each batch ---
  message("==> Applying selected method(s) per batch")
  for (bi in seq_along(batches)) {
    se_b <- batches[[bi]]
    if (tiers[bi] == "qc" && !is.null(qc_winner)) {
      message("  Batch ", batch_names[bi], ": applying ", qc_winner$name, " (QC-based)")
      qc_idx <- which(colData(se_b)$QC == "QC")
      mat <- apply_drift_candidate_to_batch(assay(se_b, 1), qc_idx,
                                             as.numeric(colData(se_b)$Injection_order),
                                             qc_winner$predict_fn, min_obs = min_cv_obs)
      assay(se_b, 1, withDimnames = FALSE) <- mat
    } else if (tiers[bi] == "sample" && !is.null(sample_winner)) {
      message("  Batch ", batch_names[bi], ": applying ", sample_winner$name,
              " (QC-free, fit on samples)")
      sample_idx <- which(colData(se_b)$QC == "Sample")
      mat <- apply_drift_candidate_to_batch(assay(se_b, 1), sample_idx,
                                             as.numeric(colData(se_b)$Injection_order),
                                             sample_winner$predict_fn, min_obs = min_cv_obs)
      assay(se_b, 1, withDimnames = FALSE) <- mat
    } else {
      message("  Batch ", batch_names[bi], ": leaving uncorrected (",
              if (tiers[bi] == "none") "insufficient QC and ltQC" else "no candidate selected",
              ")")
    }
    batches[[bi]] <- se_b
  }

  batches
}

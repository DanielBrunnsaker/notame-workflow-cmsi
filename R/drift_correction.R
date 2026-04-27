# ─────────────────────────────────────────────────────────────────────────────
# Within-batch drift correction helpers.
#
# split_by_batch()       — splits a SE into a list of per-batch SEs, sorted by injection order
# process_batch()        — notame cubic spline drift correction wrapper
# loess_correct_batch()  — QC-based LOESS drift correction with optional sample fallback
# linear_correct_batch() — QC-based linear drift correction

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
loess_correct_batch <- function(se_b, span = 0.75, fallback_to_samples = FALSE) {
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
  n_fallback    <- 0L

  samp_idx <- which(cd$QC == "Sample")

  for (i in seq_len(nrow(mat))) {
    y_qc <- as.numeric(mat[i, qc_idx])
    x_qc <- inj[qc_idx]
    ok   <- is.finite(y_qc)

    if (sum(ok) < 4) {
      if (!fallback_to_samples) { n_skipped_qc <- n_skipped_qc + 1L; next }

      # Fallback: fit LOESS through all samples instead of QC only
      y_all <- as.numeric(mat[i, samp_idx])
      x_all <- inj[samp_idx]
      ok_all <- is.finite(y_all)
      if (sum(ok_all) < 4) { n_skipped_qc <- n_skipped_qc + 1L; next }

      tryCatch({
        fit          <- loess(y ~ x, data = data.frame(x = x_all[ok_all], y = y_all[ok_all]), span = span)
        ok_inj       <- !is.na(inj)
        pred         <- rep(NA_real_, length(inj))
        pred[ok_inj] <- predict(fit, newdata = data.frame(x = inj[ok_inj]))
        med_all      <- median(y_all[ok_all])
        ratio        <- pred / med_all
        ratio[is.na(ratio) | ratio <= 0] <- 1
        mat[i, ]     <- mat[i, ] / ratio
        n_fallback   <- n_fallback + 1L
      }, error = function(e) { n_skipped_err <<- n_skipped_err + 1L })
      next
    }

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

  n_qc_corrected <- n_feat - n_skipped_qc - n_skipped_err - n_fallback
  skip_reason <- if (fallback_to_samples) "insufficient observations in QC and samples" else "insufficient QC observations"
  message("  Batch ", batch, ": drift-corrected ", n_qc_corrected + n_fallback, "/", n_feat, " features",
          if (n_fallback    > 0) paste0(" | ", n_fallback,    " via sample fallback (insufficient QC)") else "",
          if (n_skipped_qc  > 0) paste0(" | ", n_skipped_qc,  " skipped (", skip_reason, ")") else "",
          if (n_skipped_err > 0) paste0(" | ", n_skipped_err, " skipped (LOESS fit error)") else "")

  if (n_fallback > 0)
    message("  WARNING: Batch ", batch, " used sample-based LOESS fallback for ", n_fallback,
            " features. Only valid if injection order is randomised.")

  assay(se_b, 1, withDimnames = FALSE) <- mat
  se_b
}


# Linear within-batch drift correction.
# Fits lm(abundance ~ injection_order) through QC samples per batch.
# Requires >= 2 finite QC observations per feature.
linear_correct_batch <- function(se_b) {
  mat    <- assay(se_b, 1)
  cd     <- colData(se_b)
  qc_idx <- which(cd$QC == "QC")
  inj    <- as.numeric(cd$Injection_order)

  for (i in seq_len(nrow(mat))) {
    y_qc <- as.numeric(mat[i, qc_idx])
    x_qc <- inj[qc_idx]
    ok   <- is.finite(y_qc)
    if (sum(ok) < 2) next

    tryCatch({
      fit    <- lm(y ~ x, data = data.frame(x = x_qc[ok], y = y_qc[ok]))
      pred   <- predict(fit, newdata = data.frame(x = inj))
      med_qc <- median(y_qc[ok])
      ratio  <- pred / med_qc
      ratio[!is.finite(ratio) | ratio <= 0] <- 1
      mat[i, ] <- mat[i, ] / ratio
    }, error = function(e) NULL)
  }

  assay(se_b, 1, withDimnames = FALSE) <- mat
  se_b
}

# ─────────────────────────────────────────────────────────────────────────────
# TODO: write some summary

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

  if (any(!is.finite(inj))) {
    bad <- which(!is.finite(inj))
    stop("Non-finite Injection_order in batch ", unique(cd$Batch),
         ": samples ", paste(cd$Sample_ID[bad], collapse = ", "),
         " (values: ", paste(inj[bad], collapse = ", "), ")")
  }

  for (i in seq_len(nrow(mat))) {
    y_qc <- as.numeric(mat[i, qc_idx])
    x_qc <- inj[qc_idx]
    ok   <- is.finite(y_qc)
    if (sum(ok) < 4) next

    tryCatch({
      fit          <- loess(y ~ x, data = data.frame(x = x_qc[ok], y = y_qc[ok]), span = span)
      ok_inj       <- !is.na(inj)
      pred         <- rep(NA_real_, length(inj))
      pred[ok_inj] <- predict(fit, newdata = data.frame(x = inj[ok_inj]))
      med_qc       <- median(y_qc[ok])
      ratio        <- pred / med_qc
      ratio[is.na(ratio) | ratio <= 0] <- 1
      mat[i, ]     <- mat[i, ] / ratio
    }, error = function(e) NULL)
  }

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

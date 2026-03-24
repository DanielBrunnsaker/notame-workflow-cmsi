# ─────────────────────────────────────────────────────────────────────────────
# Output helpers
# ─────────────────────────────────────────────────────────────────────────────


# Write a plain feature abundance table: rows = features, columns = sample names.
write_feature_table <- function(se, file) {
  mat <- as.data.frame(assay(se, 1))
  colnames(mat) <- colData(se)$Original_name
  write.xlsx(mat, file, rowNames = TRUE)
}


# Write "client"-facing feature metadata: mz, rt, cluster info, and QC metrics.
write_feature_info <- function(se, file) {
  rd   <- as.data.frame(rowData(se))
  keep <- grep("mz|rt|cluster|RSD|D_ratio|adduct", colnames(rd), ignore.case = TRUE)
  write.xlsx(rd[, keep, drop = FALSE], file, rowNames = TRUE)
}


# ─────────────────────────────────────────────────────────────────────────────
# QC evaluation utilities
# ─────────────────────────────────────────────────────────────────────────────

# Compute robust RSD (MAD / median) of ltQC samples per feature.
# min_det_frac: minimum fraction of ltQC samples that must be non-NA (pre-imputation
# proxy — features below this are excluded to avoid imputed-value artefacts).
# Returns median ltQC RSD across all features passing the detection filter, or NA.
eval_ltqc <- function(se, min_det_frac = 0.5) {
  ltqc_idx <- which(colData(se)$QC == "ltQC")
  if (length(ltqc_idx) < 2) return(NA_real_)
  mat      <- assay(se, 1)[, ltqc_idx, drop = FALSE]
  has_mask <- "observed" %in% assayNames(se)
  if (has_mask) obs <- assay(se, "observed")[, ltqc_idx, drop = FALSE]
  rsd_per_feature <- sapply(seq_len(nrow(mat)), function(i) {
    x  <- mat[i, ]
    ok <- if (has_mask) obs[i, ] else !is.na(x)
    if (mean(ok) < min_det_frac || sum(ok) < 2) return(NA_real_)
    mad(x[ok]) / median(x[ok])
  })
  median(rsd_per_feature, na.rm = TRUE)
}


# Compute MAD(ltQC) / MAD(Sample) per feature
# min_det_frac: minimum fraction of originally-detected values required in BOTH
# ltQC and sample groups to include a feature (uses pre-imputation mask when
# available to avoid artefacts from imputed values).
# Returns median across all features passing the detection filter, or NA.
eval_ltqc_dratio <- function(se, min_det_frac = 0.5) {
  ltqc_idx   <- which(colData(se)$QC == "ltQC")
  sample_idx <- which(colData(se)$QC == "Sample")
  if (length(ltqc_idx) < 2 || length(sample_idx) < 2) return(NA_real_)
  mat      <- assay(se, 1)
  has_mask <- "observed" %in% assayNames(se)
  if (has_mask) obs <- assay(se, "observed")
  dratio_per_feature <- sapply(seq_len(nrow(mat)), function(i) {
    lt <- mat[i, ltqc_idx]
    sm <- mat[i, sample_idx]
    ok_lt <- if (has_mask) obs[i, ltqc_idx]   else is.finite(lt)
    ok_sm <- if (has_mask) obs[i, sample_idx] else is.finite(sm)
    if (mean(ok_lt) < min_det_frac || mean(ok_sm) < min_det_frac) return(NA_real_)
    if (sum(ok_lt) < 2 || sum(ok_sm) < 2) return(NA_real_)
    denom <- mad(sm[ok_sm])
    if (denom == 0) return(NA_real_)
    mad(lt[ok_lt]) / denom
  })
  median(dratio_per_feature, na.rm = TRUE)
}

# Save per-feature QC metrics and a one-row summary for one correction method.
# Files are written to interdir so results from multiple methods can be compared.
#
# Key metrics:
#   ltqc_median_RSD_r   — median robust RSD of held-out ltQC samples (unbiased)
#   ltqc_median_D_ratio — MAD(ltQC) / MAD(Sample)
#   RSD_r               — robust RSD of pooled QC samples
#   D_ratio_r           — MAD(QC) / MAD(Sample); lower = better separation of technical/biological
#
save_correction_summary <- function(se, method, interdir) {
  rd <- as.data.frame(rowData(se))

  write.csv(rd, file.path(interdir, paste0("qc_metrics_", method, ".csv")), row.names = TRUE)

  rsd  <- rd$RSD_r[!is.na(rd$RSD_r)]
  drat <- rd$D_ratio_r[!is.na(rd$D_ratio_r)]

  summary_row <- data.frame(
    method               = method,
    n_features           = nrow(se),
    ltqc_median_RSD_r    = round(eval_ltqc(se),           4),
    ltqc_median_D_ratio  = round(eval_ltqc_dratio(se),    3),
    median_RSD_r         = round(median(rsd),              3),
    median_D_ratio_r     = round(median(drat),             3),
    pct_RSD_lt_30        = round(100 * mean(rsd < 0.30),  1),
    stringsAsFactors     = FALSE
  )

  write.csv(summary_row, file.path(interdir, paste0("qc_summary_", method, ".csv")), row.names = FALSE)

  cat("\n--- QC summary:", method, "---\n")
  print(summary_row, row.names = FALSE)
  invisible(summary_row)
}


# Add per-batch QC RSD columns to rowData.
# Computes robust RSD (MAD / median) of QC samples within each batch and
# appends columns named RSD_r_B{batch} to rowData(se).
add_batch_qc_metrics <- function(se) {
  batches <- sort(unique(colData(se)$Batch))
  mat     <- assay(se, 1)

  for (b in batches) {
    idx <- which(colData(se)$Batch == b & colData(se)$QC == "QC")
    col <- paste0("RSD_r_B", b)

    if (length(idx) < 2) {
      rowData(se)[[col]] <- NA_real_
    } else {
      rowData(se)[[col]] <- apply(mat[, idx, drop = FALSE], 1, function(x) {
        ok <- !is.na(x)
        if (sum(ok) < 2) return(NA_real_)
        mad(x[ok]) / median(x[ok])
      })
    }
  }

  se
}


# Load all qc_summary_*.csv files from interdir and print a ranked comparison.
# Saves method_comparison.csv to output_dir.
compare_corrections <- function(interdir, output_dir) {
  files <- list.files(interdir, pattern = "^qc_summary_.+\\.csv$", full.names = TRUE)
  if (length(files) == 0) {
    message("No qc_summary_*.csv files found in: ", interdir)
    return(invisible(NULL))
  }
  tbl      <- do.call(rbind, lapply(files, read.csv))
  rank_col <- if ("ltqc_median_RSD_r" %in% colnames(tbl)) "ltqc_median_RSD_r" else "median_RSD_r"
  # Sort: uncorrected always first, then remaining ranked by metric
  uncorr   <- tbl[tbl$method == "uncorrected", , drop = FALSE]
  rest     <- tbl[tbl$method != "uncorrected", , drop = FALSE]
  rest     <- rest[order(rest[[rank_col]], rest$median_RSD_r), ]
  tbl      <- rbind(uncorr, rest)
  cat("\n=== Correction method comparison (best first, ranked by", rank_col, ") ===\n")
  print(tbl, row.names = FALSE)
  write.csv(tbl, file.path(output_dir, "method_comparison.csv"), row.names = FALSE)
  invisible(tbl)
}

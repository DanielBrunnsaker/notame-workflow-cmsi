# ─────────────────────────────────────────────────────────────────────────────
# Output helpers
#
# Each method/feature-set combination is written as a single multi-sheet
# workbook (via write_results_workbook) rather than separate files:
#   Peak_table       — feature abundance matrix
#   Feature_metadata — mz, rt, adduct, QC metrics, MSDIAL annotations
#   Cluster_info     — cluster membership (only present when se is clustered)
#   Settings         — run parameters and loaded package versions
# ─────────────────────────────────────────────────────────────────────────────


# Feature abundance table: rows = features, columns = sample names.
build_peak_table <- function(se) {
  mat <- as.data.frame(assay(se, 1))
  colnames(mat) <- colData(se)$Original_name
  mat
}


# Resolve MSDIAL annotations for each feature in `se`, filling in a name from
# a cluster member when the representative feature itself has no metabolite
# name. Returns annot_df rows reordered/matched to rownames(se), with an
# Annotation_source column noting where each row's data came from.
resolve_annotations <- function(se, annot_df) {
  ids <- rownames(se)
  rd  <- as.data.frame(rowData(se), check.names = FALSE)

  name_col     <- "Metabolite name"
  is_annotated <- function(x) !is.na(x) & nchar(trimws(x)) > 0 & !grepl("^Unknown$", x, ignore.case = TRUE)

  out <- annot_df[match(ids, annot_df$Feature_ID), , drop = FALSE]
  out$Feature_ID        <- ids
  rownames(out)         <- NULL
  out$Annotation_source <- "Representative feature"

  if (name_col %in% colnames(out) && "Cluster_features" %in% colnames(rd)) {
    needs <- which(!is_annotated(out[[name_col]]))
    for (i in needs) {
      id          <- ids[i]
      members_raw <- rd[id, "Cluster_features"]
      if (is.na(members_raw) || nchar(trimws(as.character(members_raw))) == 0) next
      member_ids  <- trimws(strsplit(as.character(members_raw), ";")[[1]])
      member_rows <- annot_df[annot_df$Feature_ID %in% member_ids, , drop = FALSE]
      annotated   <- member_rows[is_annotated(member_rows[[name_col]]), , drop = FALSE]
      if (nrow(annotated) == 0) next
      out[i, seq_len(ncol(annot_df))] <- annotated[1, ]
      out$Feature_ID[i]        <- id
      out$Annotation_source[i] <- paste0("Cluster member (", annotated$Feature_ID[1], ")")
    }
  }

  out
}


# "Client"-facing feature metadata: mz, rt, adduct, QC metrics, and a curated
# subset of resolved MSDIAL annotations (identity + fill/S-N quality scores),
# merged by Feature_ID. The raw MSDIAL annotation table carries dozens of
# columns (alignment ID, spectral match scores, spectrum reference file,
# comments, etc.) that aren't useful in client-facing output, so only the
# metabolite identity and fill/S-N columns are kept.
build_feature_metadata <- function(se, annot_df) {
  rd   <- as.data.frame(rowData(se), check.names = FALSE)
  keep <- grep("mz|rt|RSD|D_ratio|adduct", colnames(rd), ignore.case = TRUE)
  meta <- cbind(Feature_ID = rownames(rd), rd[, keep, drop = FALSE])

  annot      <- resolve_annotations(se, annot_df)
  annot_keep <- grep("^Feature_ID$|Metabolite name|Annotation_source|fill|S/N",
                     colnames(annot), ignore.case = TRUE)
  annot      <- annot[, annot_keep, drop = FALSE]

  merged <- merge(meta, annot, by = "Feature_ID", all.x = TRUE, sort = FALSE)
  merged[match(meta$Feature_ID, merged$Feature_ID), ]
}


# Cluster membership table (cluster ID, member features, size, MPA). Returns
# NULL when `se` carries no cluster columns, i.e. it hasn't been through
# cluster_features()/compress_clusters().
build_cluster_info <- function(se) {
  rd   <- as.data.frame(rowData(se), check.names = FALSE)
  keep <- grep("cluster|mpa", colnames(rd), ignore.case = TRUE)
  if (length(keep) == 0) return(NULL)
  cbind(Feature_ID = rownames(rd), rd[, keep, drop = FALSE])
}


# Versions of all attached (library()'d) packages, for the Settings sheet.
build_software_df <- function() {
  pkgs <- sessionInfo()$otherPkgs
  df <- if (is.null(pkgs)) {
    data.frame(Package = character(0), Version = character(0))
  } else {
    data.frame(
      Package = vapply(pkgs, `[[`, character(1), "Package"),
      Version = vapply(pkgs, `[[`, character(1), "Version"),
      stringsAsFactors = FALSE
    )
  }
  rbind(data.frame(Package = "R", Version = as.character(getRversion())), df[order(df$Package), ])
}


# Write one multi-sheet results workbook for a feature set. settings_df is a
# Parameter/Value data.frame (run parameters + per-call context such as
# method and feature-set label); a package-versions table is appended below
# it on the same sheet.
write_results_workbook <- function(se, annot_df, settings_df, file) {
  tryCatch({
    wb <- createWorkbook()

    addWorksheet(wb, "Peak_table")
    writeData(wb, "Peak_table", build_peak_table(se), rowNames = TRUE)

    addWorksheet(wb, "Feature_metadata")
    writeData(wb, "Feature_metadata", build_feature_metadata(se, annot_df), rowNames = FALSE)

    cluster_info <- build_cluster_info(se)
    if (!is.null(cluster_info)) {
      addWorksheet(wb, "Cluster_info")
      writeData(wb, "Cluster_info", cluster_info, rowNames = FALSE)
    }

    addWorksheet(wb, "Settings")
    writeData(wb, "Settings", settings_df, startRow = 1)
    writeData(wb, "Settings", build_software_df(), startRow = nrow(settings_df) + 3)

    saveWorkbook(wb, file, overwrite = TRUE)
  }, error = function(e) message("WARNING: could not write ", file, ": ", conditionMessage(e)))
}


# t-SNE perplexity scaled to sample count, for save_QC_plots(). Rtsne requires
# perplexity < (n-1)/3; the package default of 30 needs ~91+ samples to satisfy
# that, so smaller runs would otherwise error (silently, since save_QC_plots is
# called inside tryCatch) instead of producing a plot.
safe_perplexity <- function(n, max_perplexity = 30) max(2, min(max_perplexity, floor((n - 1) / 3)))


# ─────────────────────────────────────────────────────────────────────────────
# Batch quality summary
# ─────────────────────────────────────────────────────────────────────────────

# Reports per-batch missingness and a single quality score for a quick
# at-a-glance overview. Printed to the console and saved as a CSV.
#
# Quality score: median robust QC RSD (MAD/median) across all features with
# >= 2 finite QC observations in the batch. Lower = more reproducible QC signal.
# NA when a batch has fewer than 2 QC samples.
#
# Missingness: fraction of NA values in biological (Sample) cells for that batch.
report_batch_summary <- function(se, file = NULL) {
  mat     <- assay(se, 1)
  cd      <- as.data.frame(colData(se))
  batches <- unique(as.character(cd$Batch))

  rows <- lapply(batches, function(b) {
    b_idx  <- which(as.character(cd$Batch) == b)
    b_mat  <- mat[, b_idx, drop = FALSE]
    b_cd   <- cd[b_idx, ]

    n_samples <- sum(b_cd$QC == "Sample")
    n_qc      <- sum(b_cd$QC == "QC")
    n_ltqc    <- sum(b_cd$QC == "ltQC")

    samp_cols <- which(b_cd$QC == "Sample")
    miss_pct  <- if (length(samp_cols) > 0)
      round(mean(is.na(b_mat[, samp_cols, drop = FALSE])) * 100, 1)
    else NA_real_

    qc_cols  <- which(b_cd$QC == "QC")
    med_rsd  <- if (length(qc_cols) >= 2) {
      rsds <- apply(b_mat[, qc_cols, drop = FALSE], 1, function(x) {
        x <- x[is.finite(x) & x > 0]
        if (length(x) < 2) NA_real_ else mad(x) / median(x)
      })
      round(median(rsds, na.rm = TRUE) * 100, 1)
    } else NA_real_

    data.frame(
      Batch             = b,
      n_samples         = n_samples,
      n_QC              = n_qc,
      n_ltQC            = n_ltqc,
      missingness_pct   = miss_pct,
      median_QC_RSD_pct = med_rsd,
      stringsAsFactors  = FALSE
    )
  })

  out <- do.call(rbind, rows)

  cat("\n--- Batch quality summary ---\n")
  print(out, row.names = FALSE)

  if (!is.null(file)) {
    write.csv(out, file, row.names = FALSE)
    message("  Saved: ", basename(file))
  }

  invisible(out)
}

# ─────────────────────────────────────────────────────────────────────────────
# QC evaluation utilities
# ─────────────────────────────────────────────────────────────────────────────

# Compute robust RSD (MAD / median) of ltQC samples per feature.
# min_det_frac: minimum fraction of ltQC samples that must be non-NA (pre-imputation
# proxy — features below this are excluded to avoid imputed-value artefacts).
# Returns median ltQC RSD across all features passing the detection filter, or NA.
eval_ltqc <- function(se, min_det_frac = 0.5, mask = NULL) {
  ltqc_idx <- which(colData(se)$QC == "ltQC")
  if (length(ltqc_idx) < 2) return(NA_real_)
  mat      <- assay(se, 1)[, ltqc_idx, drop = FALSE]
  has_mask <- !is.null(mask)
  if (has_mask) obs <- mask[, ltqc_idx, drop = FALSE]
  rsd_per_feature <- sapply(seq_len(nrow(mat)), function(i) {
    x  <- mat[i, ]
    ok <- if (has_mask) obs[i, ] else !is.na(x)
    if (mean(ok) < min_det_frac || sum(ok) < 2) return(NA_real_)
    med <- median(x[ok])
    if (med <= 0) return(NA_real_)
    mad(x[ok]) / med
  })
  median(rsd_per_feature, na.rm = TRUE)
}


# Compute MAD(reference_group) / MAD(Sample) per feature. reference_group
# defaults to "ltQC" (the unbiased choice used everywhere else in this
# pipeline, since ltQC is never used to fit any correction method) but can be
# set to "QC" instead -- legitimate specifically for evaluating a method that
# also never fits on QC (e.g. WaveICA2.0, which fits only on injection order),
# where QC is just as much a genuine held-out reference as ltQC, and often
# has more samples to draw on.
# min_det_frac: minimum fraction of originally-detected values required in BOTH
# the reference group and sample group to include a feature (uses
# pre-imputation mask when available to avoid artefacts from imputed values).
# Returns median across all features passing the detection filter, or NA.
eval_ltqc_dratio <- function(se, min_det_frac = 0.5, mask = NULL, reference_group = "ltQC") {
  ref_idx    <- which(colData(se)$QC == reference_group)
  sample_idx <- which(colData(se)$QC == "Sample")
  if (length(ref_idx) < 2 || length(sample_idx) < 2) return(NA_real_)
  mat      <- assay(se, 1)
  has_mask <- !is.null(mask)
  dratio_per_feature <- sapply(seq_len(nrow(mat)), function(i) {
    lt <- mat[i, ref_idx]
    sm <- mat[i, sample_idx]
    ok_lt <- if (has_mask) mask[i, ref_idx]    else is.finite(lt)
    ok_sm <- if (has_mask) mask[i, sample_idx] else is.finite(sm)
    if (mean(ok_lt) < min_det_frac || mean(ok_sm) < min_det_frac) return(NA_real_)
    if (sum(ok_lt) < 2 || sum(ok_sm) < 2) return(NA_real_)
    denom <- mad(sm[ok_sm])
    if (denom == 0) return(NA_real_)
    mad(lt[ok_lt]) / denom
  })
  median(dratio_per_feature, na.rm = TRUE)
}

# Compute median pairwise euclidean distance ratio between two sample groups
# in PCA space (top n_pcs components, unit-variance scaled).
# Returns median(dist_group1) / median(dist_group2).
# Lower = group1 is tighter relative to group2.
# Use group1="QC"   for the biased (fitted) version of multivariate D_ratio.
# Use group1="ltQC" for the unbiased held-out version.
eval_dist_ratio <- function(se, group1 = "ltQC", group2 = "Sample", n_pcs = 20) {
  idx1 <- which(colData(se)$QC == group1)
  idx2 <- which(colData(se)$QC == group2)
  if (length(idx1) < 2 || length(idx2) < 2) return(NA_real_)

  mat <- t(assay(se, 1))  # samples x features

  # Drop zero-variance features (would cause scale. = TRUE to fail)
  feat_var <- apply(mat, 2, var, na.rm = TRUE)
  mat <- mat[, !is.na(feat_var) & feat_var > 0, drop = FALSE]
  if (ncol(mat) < 2) return(NA_real_)

  n_pcs_use <- min(n_pcs, nrow(mat) - 1, ncol(mat))
  pca <- tryCatch(
    prcomp(mat, center = TRUE, scale. = TRUE, rank. = n_pcs_use),
    error = function(e) NULL
  )
  if (is.null(pca)) return(NA_real_)

  scores <- pca$x
  d1 <- as.vector(dist(scores[idx1, , drop = FALSE]))
  d2 <- as.vector(dist(scores[idx2, , drop = FALSE]))
  if (length(d1) == 0 || length(d2) == 0 || median(d2) == 0) return(NA_real_)
  median(d1) / median(d2)
}

# Returns column indices (into `data`) of QC/ltQC samples that are
# multivariate outliers relative to their own batch's OTHER samples of the
# same group, in PCA space. Complements MIN_QC_SAMPLE_DETECTION
# (missingness-rate-based): that check catches a sample that failed to
# detect most features (empty injection, instrument failure); this one
# catches a sample that detects plenty of features but whose intensity
# PROFILE is anomalous (contamination, carryover, a degrading/recovering
# column) -- a failure mode detection-rate is structurally blind to.
#
# Run per batch, not pooled across the whole run: QC/ltQC intensities differ
# systematically between batches before any correction has happened, so
# pooling would either mask a real within-batch outlier under the larger
# between-batch spread, or flag a perfectly normal QC just because its batch
# differs from others.
#
# `data` should already have MIN_QC_SAMPLE_DETECTION's badly-detected samples
# removed before this runs -- including them here would distort the very
# centroid this check compares everyone else against.
#
# Distance is Euclidean, in PCA score space (top n_pcs, unit-variance
# scaled, same convention as eval_dist_ratio()), from each sample to the
# per-batch group's median score (robust centroid). A sample is flagged if
# its distance exceeds median(distances) + mad_k * mad(distances) -- a
# robust (MAD-based, not SD-based) outlier rule, consistent with this
# pipeline's general preference for MAD/median over mean/SD given how
# outlier-prone metabolomics intensity data is. Batches with fewer than
# min_n samples of this group are skipped (not enough points for a
# meaningful PCA-space distance check).
detect_qc_outliers <- function(data, group, min_n = 4, mad_k = 5, n_pcs = 5) {
  cd      <- as.data.frame(colData(data))
  mat_imp <- assay(lod2_impute(data), 1)  # local imputation, just for this diagnostic
  batches <- unique(cd$Batch)
  flagged <- integer(0)
  skipped <- character(0)  # batches skipped for insufficient/degenerate data, with reason

  for (b in batches) {
    idx <- which(cd$Batch == b & cd$QC == group)
    if (length(idx) == 0) next  # group not present in this batch at all -- not worth reporting
    if (length(idx) < min_n) {
      skipped <- c(skipped, paste0(b, " (", length(idx), " < ", min_n, ")"))
      next
    }

    sub      <- t(mat_imp[, idx, drop = FALSE])  # samples x features
    feat_var <- apply(sub, 2, var)
    sub      <- sub[, !is.na(feat_var) & feat_var > 0, drop = FALSE]
    if (ncol(sub) < 2) { skipped <- c(skipped, paste0(b, " (no variable features)")); next }

    n_pcs_use <- min(n_pcs, nrow(sub) - 1, ncol(sub))
    pca <- tryCatch(prcomp(sub, center = TRUE, scale. = TRUE, rank. = n_pcs_use),
                     error = function(e) NULL)
    if (is.null(pca)) { skipped <- c(skipped, paste0(b, " (PCA failed)")); next }

    scores   <- pca$x
    centroid <- apply(scores, 2, median)
    dists    <- sqrt(rowSums(sweep(scores, 2, centroid)^2))

    mad_dist <- mad(dists)
    if (mad_dist == 0) next  # every point equidistant from centroid -- nothing to flag, not a skip
    thresh <- median(dists) + mad_k * mad_dist
    bad    <- idx[dists > thresh]

    if (length(bad) > 0) {
      message("  Batch ", b, ": ", length(bad), " ", group, " outlier sample(s) ",
              "(PCA distance > ", round(thresh, 2), "): ",
              paste(cd$Sample_ID[bad], collapse = ", "))
      flagged <- c(flagged, bad)
    }
  }

  message("  ", group, " outlier check: ", length(flagged), " flagged",
          if (length(skipped) > 0) paste0("; skipped ", length(skipped), " batch(es) -- ",
                                           paste(skipped, collapse = ", ")) else "")

  flagged
}


# Compute within-batch pairwise Euclidean distance preservation.
# For each batch, computes all pairwise distances between biological samples in
# both corrected and raw data, then correlates the two distance vectors (Spearman).
# Returns the median correlation across batches.
#
# This avoids the cross-batch confound of the feature-wise signal_preservation_r:
# within a batch the batch offset is constant and does not inflate pairwise
# distances, so raw within-batch distances already reflect biology.
eval_within_batch_dist_preservation <- function(se, ref_mat) {
  if (is.null(ref_mat)) return(NA_real_)

  cd       <- as.data.frame(colData(se))
  samp_idx <- which(cd$QC == "Sample")
  if (length(samp_idx) < 4) return(NA_real_)

  batches <- unique(cd$Batch[samp_idx])

  batch_cors <- vapply(batches, function(b) {
    b_idx <- samp_idx[cd$Batch[samp_idx] == b]
    if (length(b_idx) < 3) return(NA_real_)

    ids        <- cd$Sample_ID[b_idx]
    shared_ids <- intersect(ids, colnames(ref_mat))
    if (length(shared_ids) < 3) return(NA_real_)

    cur_cols <- b_idx[ids %in% shared_ids]
    ref_cols <- match(shared_ids, colnames(ref_mat))

    shared_feat <- intersect(rownames(assay(se, 1)), rownames(ref_mat))
    if (length(shared_feat) < 2) return(NA_real_)

    mat_cur <- assay(se, 1)[shared_feat, cur_cols, drop = FALSE]
    mat_ref <- ref_mat[shared_feat, ref_cols, drop = FALSE]

    d_cur <- as.vector(dist(t(mat_cur)))
    d_ref <- as.vector(dist(t(mat_ref)))

    ok <- is.finite(d_cur) & is.finite(d_ref)
    if (sum(ok) < 3) return(NA_real_)
    cor(d_cur[ok], d_ref[ok], method = "spearman")
  }, numeric(1))

  round(median(batch_cors, na.rm = TRUE), 3)
}


# Compute median feature-wise Spearman correlation between corrected and a
# reference dataset (both restricted to biological samples only).
# Features present in both are used; returns NA if fewer than 2 features match.
eval_signal_preservation <- function(se, ref_mat) {
  if (is.null(ref_mat)) return(NA_real_)
  sample_idx <- which(colData(se)$QC == "Sample")
  if (length(sample_idx) < 3) return(NA_real_)

  mat <- assay(se, 1)[, sample_idx, drop = FALSE]

  # Align to features present in both
  shared <- intersect(rownames(mat), rownames(ref_mat))
  if (length(shared) < 2) return(NA_real_)

  # Align columns (samples) by Sample_ID
  cd      <- colData(se)[sample_idx, , drop = FALSE]
  ref_ids <- colnames(ref_mat)
  cur_ids <- cd$Sample_ID

  shared_ids <- intersect(cur_ids, ref_ids)
  if (length(shared_ids) < 3) return(NA_real_)

  cur_cols <- match(shared_ids, cur_ids)
  ref_cols <- match(shared_ids, ref_ids)

  mat_cur <- mat[shared, cur_cols, drop = FALSE]
  mat_ref <- ref_mat[shared, ref_cols, drop = FALSE]

  # Feature-wise Spearman correlation
  cors <- sapply(seq_len(nrow(mat_cur)), function(i) {
    x <- mat_cur[i, ]
    y <- mat_ref[i, ]
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 3) return(NA_real_)
    cor(x[ok], y[ok], method = "spearman")
  })
  median(cors, na.rm = TRUE)
}


# Median R² of batch as a factor across biological samples, per feature.
# Computed entirely on corrected data — no raw reference needed.
# Lower = less remaining batch structure after correction.
eval_remaining_batch_r2 <- function(se) {
  samp_idx <- which(colData(se)$QC == "Sample")
  if (length(samp_idx) < 4) return(NA_real_)

  batch <- as.character(colData(se)$Batch[samp_idx])
  if (length(unique(batch)) < 2) return(NA_real_)

  mat <- assay(se, 1)[, samp_idx, drop = FALSE]

  r2_per_feature <- apply(mat, 1, function(x) {
    ok <- is.finite(x)
    if (sum(ok) < 4) return(NA_real_)
    x_ok <- x[ok]
    b_ok <- batch[ok]
    if (length(unique(b_ok)) < 2) return(NA_real_)
    grand_mean <- mean(x_ok)
    ss_total   <- sum((x_ok - grand_mean)^2)
    if (ss_total == 0) return(NA_real_)
    batch_means <- tapply(x_ok, b_ok, mean)
    batch_n     <- tapply(x_ok, b_ok, length)
    ss_between  <- sum(batch_n * (batch_means - grand_mean)^2)
    ss_between / ss_total
  })

  round(median(r2_per_feature, na.rm = TRUE), 3)
}


# Median absolute Spearman correlation of QC sample abundances with injection
# order, computed per feature per batch, then summarised as median across all
# feature × batch combinations.
# Computed entirely on corrected data — no raw reference needed.
# Lower = less remaining within-batch drift after correction.
eval_remaining_drift <- function(se) {
  qc_idx <- which(colData(se)$QC == "QC")
  if (length(qc_idx) < 4) return(NA_real_)

  cd      <- as.data.frame(colData(se))
  batches <- unique(cd$Batch[qc_idx])
  mat     <- assay(se, 1)

  batch_cors <- vapply(batches, function(b) {
    b_idx <- qc_idx[cd$Batch[qc_idx] == b]
    if (length(b_idx) < 3) return(NA_real_)

    inj_order <- cd$Injection_order[b_idx]
    sub_mat   <- mat[, b_idx, drop = FALSE]

    feature_cors <- apply(sub_mat, 1, function(x) {
      ok <- is.finite(x)
      if (sum(ok) < 3) return(NA_real_)
      abs(cor(x[ok], inj_order[ok], method = "spearman"))
    })

    median(feature_cors, na.rm = TRUE)
  }, numeric(1))

  round(median(batch_cors, na.rm = TRUE), 3)
}


# Test whether a pooled reference group (QC or ltQC) forms a single homogeneous
# population after correction, using PERMANOVA (centroid) and PERMDISP
# (dispersion) from vegan. Requires vegan; returns NAs silently if not
# installed or fewer than 2 batches.
#
# group = "QC" tests the same samples several methods (notame's RUV,
# loess_combat, pmp_qcrsc, cordbat_only) use to *fit* their correction — a
# pass there is partly circular, since the correction was built to make QC
# match across batches. group = "ltQC" tests samples never used to fit any
# correction, so a batch signature found there is unbiased evidence of a
# real, unaddressed technical effect (assuming ltQC is genuinely the same
# reference material injected throughout — no biological variation to
# confound the result, unlike testing this on Sample). Prefer ltQC when
# available; with as few as ~3 ltQC per batch the test is underpowered (weak
# ability to detect a real but modest batch effect), so treat a non-significant
# result as "no strong evidence found," not proof of no batch effect — a
# significant result at that sample size is the more trustworthy direction.
#
#   {group}_permanova_r2  — R²(Batch) from PERMANOVA on the group; lower = better
#   {group}_permanova_p   — p-value; want > 0.05 (centroids don't differ by batch)
#   {group}_permdisp_p    — p-value from PERMDISP; want > 0.05 (homogeneous spread)
eval_qc_homogeneity <- function(se, group = "QC") {
  na_out <- list(permanova_r2 = NA_real_,
                 permanova_p  = NA_real_,
                 permdisp_p   = NA_real_)

  qc_idx   <- which(colData(se)$QC == group)
  if (length(qc_idx) < 4) return(na_out)
  batch_qc <- as.character(colData(se)$Batch[qc_idx])
  if (length(unique(batch_qc)) < 2 || any(table(batch_qc) < 2)) return(na_out)
  if (!requireNamespace("vegan", quietly = TRUE)) return(na_out)

  mat     <- assay(se, 1)[, qc_idx, drop = FALSE]
  mat_log <- suppressWarnings(log2(mat))
  mat_log[!is.finite(mat_log)] <- NA
  for (j in seq_len(ncol(mat_log))) {
    na_j <- !is.finite(mat_log[, j])
    if (any(na_j)) mat_log[na_j, j] <- mean(mat_log[, j], na.rm = TRUE)
  }
  feat_var <- apply(mat_log, 1, var, na.rm = TRUE)
  mat_log  <- mat_log[is.finite(feat_var) & feat_var > 0, , drop = FALSE]
  if (nrow(mat_log) < 2) return(na_out)

  d <- dist(t(mat_log))

  pm <- tryCatch(
    vegan::adonis2(d ~ Batch,
                   data         = data.frame(Batch = batch_qc),
                   permutations = 999, by = "margin"),
    error = function(e) NULL
  )
  pm_df   <- if (!is.null(pm)) as.data.frame(pm) else NULL
  p_col   <- if (!is.null(pm_df)) grep("Pr\\(>F\\)", colnames(pm_df), value = TRUE)[1] else NA
  perm_r2 <- if (!is.null(pm_df)) round(pm_df["Batch", "R2"],     4) else NA_real_
  perm_p  <- if (!is.null(pm_df) && !is.na(p_col))
               round(pm_df["Batch", p_col], 4) else NA_real_

  bd     <- tryCatch(vegan::betadisper(d, group = batch_qc), error = function(e) NULL)
  disp_p <- if (!is.null(bd)) tryCatch({
    pt  <- vegan::permutest(bd, permutations = 999)
    pc  <- grep("Pr\\(>F\\)", colnames(pt$tab), value = TRUE)[1]
    round(pt$tab[1, pc], 4)
  }, error = function(e) NA_real_) else NA_real_

  list(permanova_r2 = perm_r2,
       permanova_p  = perm_p,
       permdisp_p   = disp_p)
}


# Save per-feature QC metrics and a one-row summary for one correction method.
# Files are written to interdir so results from multiple methods can be compared.
#
# Key metrics:
#   ltqc_median_RSD_r      — median robust RSD of held-out ltQC samples (unbiased)
#   ltqc_median_D_ratio    — MAD(ltQC) / MAD(Sample); per-feature, unbiased
#   ltqc_dist_ratio        — median pairwise dist(ltQC) / dist(Sample) in PCA space (unbiased)
#   RSD_r                  — robust RSD of pooled QC samples
#   D_ratio_r              — MAD(QC) / MAD(Sample); lower = better separation
#   median_sample_MAD      — median of per-feature MAD across biological samples
#   within_batch_dist_r    — median Spearman cor of within-batch pairwise distances
#                            (corrected vs raw); unaffected by between-batch structure
#   remaining_batch_r2     — median R² of batch as factor on biological samples;
#                            lower = less remaining batch effect
#   remaining_drift_r      — median absolute Spearman cor of QC abundance with
#                            injection order within batches; lower = less drift remaining
#   qc_permanova_r2        — R²(Batch) from PERMANOVA on QC samples (lower = better).
#                            QC is used to *fit* several correction methods, so this
#                            is a partly circular check for those methods — see
#                            ltqc_permanova_r2 for the unbiased version.
#   qc_permanova_p         — p-value; want > 0.05 (QC don't cluster by batch)
#   qc_permdisp_p          — PERMDISP p-value; want > 0.05 (homogeneous QC spread)
#   ltqc_permanova_r2      — as qc_permanova_r2 but on ltQC, which is never used to fit
#                            any correction — an unbiased (if underpowered at typical
#                            ltQC counts) check for a remaining batch signature
#   ltqc_permanova_p       — p-value; want > 0.05 (ltQC don't cluster by batch)
#   ltqc_permdisp_p        — PERMDISP p-value; want > 0.05 (homogeneous ltQC spread)
#
save_correction_summary <- function(se, method, interdir, obs_mask = NULL, raw_ref = NULL) {
  rd <- as.data.frame(rowData(se))

  # method may contain colons (a "drift:basis:batch" spec) -- fine as a data
  # value (method column below, printed labels), but not as part of a
  # filename, so filenames use the sanitized form while the column/label
  # keeps the original, more readable spec string.
  method_id <- sanitize_method_id(method)
  write.csv(rd, file.path(interdir, paste0("qc_metrics_", method_id, ".csv")), row.names = TRUE)

  rsd  <- rd$RSD_r[!is.na(rd$RSD_r)]
  drat <- rd$D_ratio_r[!is.na(rd$D_ratio_r)]

  # Median sample MAD: absolute biological variance (higher = more signal retained)
  sample_idx <- which(colData(se)$QC == "Sample")
  sample_mads <- if (length(sample_idx) >= 2) {
    apply(assay(se, 1)[, sample_idx, drop = FALSE], 1, function(x) {
      ok <- is.finite(x)
      if (sum(ok) < 2) return(NA_real_)
      mad(x[ok])
    })
  } else {
    NA_real_
  }

  summary_row <- data.frame(
    method                = method,
    n_features            = nrow(se),
    ltqc_median_RSD_r     = round(eval_ltqc(se,             mask = obs_mask), 4),
    ltqc_median_D_ratio   = round(eval_ltqc_dratio(se,      mask = obs_mask), 3),
    ltqc_dist_ratio       = round(eval_dist_ratio(se,        group1 = "ltQC"), 3),
    median_RSD_r          = round(median(rsd),               3),
    median_D_ratio_r      = round(median(drat),              3),
    pct_RSD_lt_30         = round(100 * mean(rsd < 0.30),   1),
    median_sample_MAD     = round(median(sample_mads, na.rm = TRUE), 2),
    within_batch_dist_r   = round(eval_within_batch_dist_preservation(se, raw_ref), 3),
    remaining_batch_r2    = round(eval_remaining_batch_r2(se), 3),
    remaining_drift_r     = round(eval_remaining_drift(se), 3),
    stringsAsFactors      = FALSE
  )

  qc_homo <- eval_qc_homogeneity(se, group = "QC")
  summary_row$qc_permanova_r2 <- qc_homo$permanova_r2
  summary_row$qc_permanova_p  <- qc_homo$permanova_p
  summary_row$qc_permdisp_p   <- qc_homo$permdisp_p

  # ltQC is never used to fit any correction method, so this is an unbiased
  # (if underpowered, at typical ltQC counts) check — see eval_qc_homogeneity's
  # documentation for why this differs from the qc_permanova_* check above.
  ltqc_homo <- eval_qc_homogeneity(se, group = "ltQC")
  summary_row$ltqc_permanova_r2 <- ltqc_homo$permanova_r2
  summary_row$ltqc_permanova_p  <- ltqc_homo$permanova_p
  summary_row$ltqc_permdisp_p   <- ltqc_homo$permdisp_p

  write.csv(summary_row, file.path(interdir, paste0("qc_summary_", method_id, ".csv")), row.names = FALSE)

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
    col <- paste0("RSD_r_", gsub("[^A-Za-z0-9]", "_", as.character(b)))

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

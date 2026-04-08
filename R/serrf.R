# ─────────────────────────────────────────────────────────────────────────────
# SERRF — Systematic Error Removal using Random Forest
#
# serrfR() is adapted from the SERRF normalization package
# (Fan et al., Analytical Chemistry 2019).
# Original source: normalization.R (slfan2013 / github.com/slfan2013/SERRF)
#
# Requirements: ranger, parallel
# ─────────────────────────────────────────────────────────────────────────────

# Core SERRF algorithm.
# Trains a random-forest model on QC samples using correlated features as
# predictors, then applies the learned drift model to biological samples.
#
# Args:
#   train       features × n_qc matrix  (QC samples only)
#   target      features × n_samp matrix (biological samples)
#   num         number of correlated features to use as RF predictors (default 10)
#   batch.      factor of length (n_qc + n_samp) indicating batch membership
#   time.       numeric vector of length (n_qc + n_samp) — global injection order
#   sampleType. character vector: "qc" for train samples, "sample" for target
#   cl          parallel cluster from parallel::makeCluster()
#
# Returns list(normed_train, normed_target) — corrected matrices, same dims as input.
serrfR <- function(train,
                   target,
                   num        = 10,
                   batch.,
                   time.,
                   sampleType.,
                   cl) {

  all    <- cbind(train, target)
  normalized <- rep(0, ncol(all))

  # Impute zeros / NAs with small noise around per-batch minimum
  for (j in seq_len(nrow(all))) {
    for (b in levels(batch.)) {
      col_b <- batch. %in% b
      x_b   <- all[j, col_b]
      ok_b  <- !is.na(x_b) & is.finite(x_b)
      if (sum(ok_b) == 0) next
      mn <- min(x_b[ok_b]) + 1
      sd_val <- 0.1 * (min(x_b[ok_b]) + 0.1)
      # replace zeros
      zero_b <- col_b & (all[j, ] == 0)
      if (any(zero_b))
        all[j, zero_b] <- rnorm(sum(zero_b), mean = mn, sd = sd_val)
      # replace NAs
      na_b <- col_b & is.na(all[j, ])
      if (any(na_b))
        all[j, na_b] <- rnorm(sum(na_b), mean = 0.5 * (mn - 1) + 1, sd = sd_val)
    }
  }

  # Per-batch Spearman correlation matrices (for feature selection)
  corrs_train  <- list()
  corrs_target <- list()
  for (b in levels(batch.)) {
    col_qc   <- sampleType. == "qc"   & batch. %in% b
    col_samp <- sampleType. == "sample" & batch. %in% b

    tr_scale <- t(apply(train[, batch.[sampleType. == "qc"] %in% b, drop = FALSE], 1, scale))
    tg_sub   <- target[, batch.[sampleType. == "sample"] %in% b, drop = FALSE]

    if (ncol(tg_sub) == 1) {
      tg_scale <- scale(tg_sub)
    } else {
      tg_scale <- t(apply(tg_sub, 1, scale))
    }

    corrs_train[[b]]  <- cor(t(tr_scale), method = "spearman")
    corrs_target[[b]] <- cor(t(tg_scale), method = "spearman")
  }

  # Parallel per-feature correction
  pred <- parallel::parSapply(
    cl, X = seq_len(nrow(all)),
    function(j, all, batch., ranger, sampleType., time., num, corrs_train, corrs_target) {
      normalized <- rep(0, ncol(all))

      for (b in levels(batch.)) {
        e_current_batch  <- all[, batch. %in% b, drop = FALSE]
        corr_train       <- corrs_train[[b]]
        corr_target      <- corrs_target[[b]]

        corr_train_order  <- order(abs(corr_train[, j]),  decreasing = TRUE)
        corr_target_order <- order(abs(corr_target[, j]), decreasing = TRUE)

        sel_var <- c()
        l <- num
        while (length(sel_var) < num) {
          sel_var <- intersect(corr_train_order[seq_len(l)],
                               corr_target_order[seq_len(l)])
          sel_var <- sel_var[sel_var != j]
          l <- l + 1
          if (l > nrow(all)) break
        }

        train_idx  <- sampleType.[batch. %in% b] == "qc"
        train_y    <- scale(e_current_batch[j, train_idx], scale = FALSE)
        train_x    <- apply(e_current_batch[sel_var, train_idx, drop = FALSE], 1, scale)

        test_sub   <- e_current_batch[sel_var, !train_idx, drop = FALSE]
        if (ncol(test_sub) == 1) {
          test_x <- t(scale(test_sub))
        } else {
          test_x <- apply(test_sub, 1, scale)
        }
        if (!is.matrix(test_x)) test_x <- t(test_x)

        # Drop features with NA in train or test
        good <- apply(train_x, 2, function(x) sum(is.na(x)) == 0) &
                apply(test_x,  2, function(x) sum(is.na(x)) == 0)
        train_x <- train_x[, good, drop = FALSE]
        test_x  <- test_x[,  good, drop = FALSE]
        if (!is.matrix(test_x)) test_x <- t(test_x)

        norm_b <- e_current_batch[j, ]

        if (ncol(train_x) == 0) {
          # Not enough predictor features — leave uncorrected
          normalized[batch. %in% b] <- norm_b
          next
        }

        train_df <- data.frame(y = as.vector(train_y),
                               setNames(as.data.frame(train_x),
                                        paste0("V", seq_len(ncol(train_x)))))
        test_df  <- setNames(as.data.frame(test_x),
                             paste0("V", seq_len(ncol(test_x))))

        model <- tryCatch(
          ranger::ranger(y ~ ., data = train_df, num.trees = 500),
          error = function(e) NULL
        )
        if (is.null(model)) {
          normalized[batch. %in% b] <- norm_b
          next
        }

        qc_mean_global   <- mean(all[j, sampleType. == "qc"],    na.rm = TRUE)
        samp_median_global <- median(all[j, sampleType. == "sample"], na.rm = TRUE)
        qc_mean_local    <- mean(e_current_batch[j, train_idx],   na.rm = TRUE)
        samp_mean_local  <- mean(e_current_batch[j, !train_idx],  na.rm = TRUE)

        pred_train <- ranger::predictions(predict(model, data = train_df))
        pred_test  <- ranger::predictions(predict(model, data = test_df))

        norm_b[train_idx]  <- e_current_batch[j, train_idx]  /
          ((pred_train + qc_mean_local) / qc_mean_global)
        norm_b[!train_idx] <- e_current_batch[j, !train_idx] /
          ((pred_test  + samp_mean_local) / samp_median_global)

        # Clamp extreme outliers (coef = 3 IQR) back to additive correction
        out_flags <- norm_b[!train_idx] %in% boxplot.stats(norm_b, coef = 3)$out
        if (any(out_flags)) {
          norm_b[!train_idx][out_flags] <-
            (e_current_batch[j, !train_idx] -
               (pred_test + samp_mean_local - samp_median_global))[out_flags]
        }
        norm_b[!train_idx][norm_b[!train_idx] < 0] <-
          e_current_batch[j, !train_idx][norm_b[!train_idx] < 0]

        # Rescale to preserve global medians
        med_qc   <- median(norm_b[train_idx],  na.rm = TRUE)
        med_samp <- median(norm_b[!train_idx], na.rm = TRUE)
        if (is.finite(med_qc)   && med_qc   != 0)
          norm_b[train_idx]  <- norm_b[train_idx]  / (med_qc   / qc_mean_global)
        if (is.finite(med_samp) && med_samp != 0)
          norm_b[!train_idx] <- norm_b[!train_idx] / (med_samp / samp_median_global)

        norm_b[!is.finite(norm_b)] <- rnorm(
          sum(!is.finite(norm_b)),
          sd = sd(norm_b[is.finite(norm_b)], na.rm = TRUE) * 0.01
        )

        normalized[batch. %in% b] <- norm_b
      }

      normalized
    },
    all, batch., ranger, sampleType., time., num, corrs_train, corrs_target
  )

  normed <- t(pred)

  # Post-process: replace any remaining negatives / NAs with small positive noise
  fix_negatives <- function(m) {
    for (i in seq_len(nrow(m))) {
      pos_vals <- m[i, is.finite(m[i, ]) & m[i, ] > 0]
      if (length(pos_vals) == 0) next
      mn_pos <- min(pos_vals)
      sd_pos <- sd(m[i, is.finite(m[i, ])], na.rm = TRUE) * 0.1
      m[i, is.na(m[i, ])]   <- rnorm(sum(is.na(m[i, ])),   mean = mn_pos, sd = sd_pos)
      m[i, m[i, ] < 0]      <- runif(sum(m[i, ] < 0, na.rm = TRUE)) * mn_pos
    }
    m
  }

  normed_train  <- fix_negatives(normed[, sampleType. == "qc",    drop = FALSE])
  normed_target <- fix_negatives(normed[, sampleType. == "sample", drop = FALSE])

  list(normed_train = normed_train, normed_target = normed_target)
}


# SummarizedExperiment wrapper for SERRF.
#
# Adapts the notame SE format (QC/Sample/ltQC labels, per-batch Injection_order)
# to SERRF's expected matrix format, runs correction, and reconstructs the SE.
#
# ltQC samples are corrected using a second serrfR call (held-out validation),
# matching how the original SERRF code handles validation sample types.
#
# Args:
#   data               SummarizedExperiment from notame pipeline (pre-imputation)
#   num                number of correlated features to use as RF predictors
#   detectcores_ratio  fraction of available CPU cores to use for parallelism
correct_serrf <- function(data, num = 10, detectcores_ratio = 0.5) {
  library(ranger)
  library(parallel)

  # ── 1. Pre-imputation detection mask (for unbiased ltQC metrics) ────────────
  obs_mask <- !is.na(assay(data, 1))

  # ── 2. Extract matrices and metadata ────────────────────────────────────────
  e_all <- assay(data, 1)
  cd    <- as.data.frame(colData(data))

  # Build a globally-unique injection time by combining batch + within-batch order.
  # This ensures time. is monotone across batches even when Injection_order restarts.
  batch_rank <- as.integer(factor(cd$Batch, levels = sort(unique(cd$Batch))))
  max_inj    <- max(cd$Injection_order, na.rm = TRUE)
  global_time <- (batch_rank - 1L) * (max_inj + 1L) + cd$Injection_order

  stype <- ifelse(cd$QC == "QC",     "qc",
           ifelse(cd$QC == "Sample",  "sample",
           ifelse(cd$QC == "ltQC",    "ltqc", NA_character_)))

  # ── 3. SERRF half-minimum imputation (applied only to NA cells) ─────────────
  for (i in seq_len(nrow(e_all))) {
    row_vals <- e_all[i, ]
    finite_min <- min(row_vals, na.rm = TRUE)
    if (!is.finite(finite_min)) finite_min <- 1
    e_all[i, is.na(row_vals)] <- 0.5 * finite_min
  }

  # ── 4. Index by sample type ──────────────────────────────────────────────────
  qc_idx   <- which(stype == "qc")
  samp_idx <- which(stype == "sample")
  ltqc_idx <- which(stype == "ltqc")

  if (length(qc_idx) < 4)
    stop("SERRF requires at least 4 QC samples; found: ", length(qc_idx))
  if (length(samp_idx) < 2)
    stop("SERRF requires at least 2 biological samples; found: ", length(samp_idx))

  e_qc   <- e_all[, qc_idx,   drop = FALSE]
  e_samp <- e_all[, samp_idx, drop = FALSE]
  e_ltqc <- if (length(ltqc_idx) > 0) e_all[, ltqc_idx, drop = FALSE] else NULL

  batch_all <- as.character(cd$Batch)
  time_all  <- global_time

  # ── 5. Parallel cluster ──────────────────────────────────────────────────────
  n_cores <- max(1L, floor(parallel::detectCores() * detectcores_ratio))
  message("==> SERRF: starting ", n_cores, " parallel worker(s)")
  cl <- parallel::makeCluster(n_cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterEvalQ(cl, library(ranger))

  # ── 6. Correct biological samples ───────────────────────────────────────────
  message("==> SERRF: correcting Sample (n=", length(samp_idx), ") using QC (n=", length(qc_idx), ")")

  idx_qs   <- c(qc_idx, samp_idx)
  stype_qs <- c(rep("qc", length(qc_idx)), rep("sample", length(samp_idx)))
  batch_qs <- factor(batch_all[idx_qs])
  time_qs  <- time_all[idx_qs]

  res_samp <- serrfR(
    train       = e_qc,
    target      = e_samp,
    num         = num,
    batch.      = batch_qs,
    time.       = time_qs,
    sampleType. = stype_qs,
    cl          = cl
  )

  e_qc_corr   <- res_samp$normed_train
  e_samp_corr <- res_samp$normed_target

  # ── 7. Correct ltQC samples (held-out validation) ───────────────────────────
  e_ltqc_corr <- NULL
  if (!is.null(e_ltqc)) {
    message("==> SERRF: correcting ltQC (n=", length(ltqc_idx), ") as held-out validation")
    idx_ql   <- c(qc_idx, ltqc_idx)
    stype_ql <- c(rep("qc", length(qc_idx)), rep("sample", length(ltqc_idx)))
    batch_ql <- factor(batch_all[idx_ql])
    time_ql  <- time_all[idx_ql]

    res_ltqc <- serrfR(
      train       = e_qc,
      target      = e_ltqc,
      num         = num,
      batch.      = batch_ql,
      time.       = time_ql,
      sampleType. = stype_ql,
      cl          = cl
    )
    e_ltqc_corr <- res_ltqc$normed_target
  }

  # ── 8. Reassemble SE in original column order ─────────────────────────────
  combined <- data
  e_new    <- assay(combined, 1)

  e_new[, qc_idx]   <- e_qc_corr
  colnames(e_new)[qc_idx] <- colnames(e_qc_corr)

  e_new[, samp_idx] <- e_samp_corr
  colnames(e_new)[samp_idx] <- colnames(e_samp_corr)

  if (!is.null(e_ltqc_corr)) {
    e_new[, ltqc_idx] <- e_ltqc_corr
    colnames(e_new)[ltqc_idx] <- colnames(e_ltqc_corr)
  }

  assay(combined, 1, withDimnames = FALSE) <- e_new

  list(pre = combined, post = combined, obs_mask = obs_mask)
}

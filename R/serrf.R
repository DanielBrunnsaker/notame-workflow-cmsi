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
                   sampleType.,
                   cl) {

  all        <- cbind(train, target)
  normalized <- rep(0, ncol(all))

  # ── Impute zeros and NAs with small noise around per-batch minimum ───────────
  # Faithful to original: zeros → noise at min+1; NAs → noise at 0.5*min+1
  for (j in seq_len(nrow(all))) {
    for (b in levels(batch.)) {
      col_b  <- which(batch. %in% b)
      x_b    <- all[j, col_b]
      ok_b   <- !is.na(x_b)
      if (sum(ok_b) == 0) next
      mn     <- min(x_b[ok_b])
      sd_val <- 0.1 * (mn + 0.1)
      # zeros → noise around min+1
      z_idx <- col_b[x_b == 0 & !is.na(x_b)]
      if (length(z_idx) > 0)
        all[j, z_idx] <- rnorm(length(z_idx), mean = mn + 1, sd = sd_val)
      # NAs → noise around 0.5*min+1
      na_idx <- col_b[is.na(x_b)]
      if (length(na_idx) > 0)
        all[j, na_idx] <- rnorm(length(na_idx), mean = 0.5 * mn + 1, sd = sd_val)
    }
  }

  # ── Per-batch Spearman correlation matrices (for feature selection) ──────────
  # Train: row-scale (features × samples → scale across samples per feature)
  # Target: column-scale — matches actual original behaviour (is.null check is
  #         always FALSE on a matrix, so else-branch always executes)
  corrs_train  <- list()
  corrs_target <- list()
  for (b in levels(batch.)) {
    tr_cols <- train[, batch.[sampleType. == "qc"] %in% b, drop = FALSE]
    tg_cols <- target[, batch.[sampleType. == "sample"] %in% b, drop = FALSE]

    tr_scale <- t(apply(tr_cols, 1, scale))
    tg_scale <- scale(tg_cols)          # column-scale: matches original else-branch

    corrs_train[[b]]  <- cor(t(tr_scale), method = "spearman")
    corrs_target[[b]] <- cor(t(tg_scale), method = "spearman")
  }

  # ── Parallel per-feature correction ─────────────────────────────────────────
  pred <- parallel::parSapply(
    cl, X = seq_len(nrow(all)),
    function(j, all, batch., sampleType., num, corrs_train, corrs_target) {

      normalized <- rep(0, ncol(all))

      for (b in levels(batch.)) {
        e_current_batch <- all[, batch. %in% b, drop = FALSE]
        corr_train      <- corrs_train[[b]]
        corr_target     <- corrs_target[[b]]

        corr_train_order  <- order(abs(corr_train[, j]),  decreasing = TRUE)
        corr_target_order <- order(abs(corr_target[, j]), decreasing = TRUE)

        # Expand search window until num features found in both top-lists
        sel_var <- c()
        l <- num
        while (length(sel_var) < num) {
          sel_var <- intersect(corr_train_order[seq_len(l)],
                               corr_target_order[seq_len(l)])
          sel_var <- sel_var[sel_var != j]
          l <- l + 1
          if (l > nrow(all)) break   
        }

        train_idx <- sampleType.[batch. %in% b] == "qc"

        # y: mean-centred QC signal for feature j
        train_data_y <- scale(e_current_batch[j, train_idx], scale = FALSE)

        # x: column-scaled predictor features (samples are observations)
        train_data_x <- apply(
          e_current_batch[sel_var, train_idx, drop = FALSE], 1, scale)

        if (is.null(dim(e_current_batch[sel_var, !train_idx]))) {
          test_data_x <- t(scale(e_current_batch[sel_var, !train_idx]))
        } else {
          test_data_x <- apply(e_current_batch[sel_var, !train_idx], 1, scale)
        }

        # First NA filter: drop predictor columns with any NA in train
        train_NA_index <- apply(train_data_x, 2, function(x) sum(is.na(x)) > 0)
        train_data_x   <- train_data_x[, !train_NA_index, drop = FALSE]
        test_data_x    <- test_data_x[,  !train_NA_index, drop = FALSE]
        if (!is.matrix(test_data_x)) test_data_x <- t(test_data_x)

        # Second NA filter: drop columns with NA in either train or test
        good_column  <- apply(train_data_x, 2, function(x) sum(is.na(x)) == 0) &
                        apply(test_data_x,  2, function(x) sum(is.na(x)) == 0)
        train_data_x <- train_data_x[, good_column, drop = FALSE]
        test_data_x  <- test_data_x[,  good_column, drop = FALSE]
        if (!is.matrix(test_data_x)) test_data_x <- t(test_data_x)

        train_data <- data.frame(y = as.vector(train_data_y), train_data_x)

        norm <- e_current_batch[j, ]

        if (ncol(train_data) == 1) {
          normalized[batch. %in% b] <- norm
          next
        }

        colnames(train_data) <- c("y", paste0("V", seq_len(ncol(train_data) - 1)))
        model <- tryCatch(
          ranger::ranger(y ~ ., data = train_data),
          error = function(e) NULL
        )
        if (is.null(model)) {
          normalized[batch. %in% b] <- norm
          next
        }

        test_data           <- as.data.frame(test_data_x)
        colnames(test_data) <- colnames(train_data)[-1]

        pred_train <- ranger::predictions(predict(model, data = train_data))
        pred_test  <- ranger::predictions(predict(model, data = test_data))

        qc_mean_global     <- mean(all[j, sampleType. == "qc"],      na.rm = TRUE)
        qc_median_global   <- median(all[j, sampleType. == "qc"],    na.rm = TRUE)
        samp_median_global <- median(all[j, sampleType. == "sample"], na.rm = TRUE)
        qc_mean_local      <- mean(e_current_batch[j, train_idx],    na.rm = TRUE)
        samp_mean_local    <- mean(e_current_batch[j, !train_idx],   na.rm = TRUE)

        # Multiplicative correction (original lines 851, 854)
        # QC uses global mean as reference (line 851); samples use global median (line 854). 
        # Not sure why this is the case, but it aligns with the original implementation.
        norm[train_idx]  <- e_current_batch[j, train_idx] /
          ((pred_train + qc_mean_local) / qc_mean_global)
        norm[!train_idx] <- e_current_batch[j, !train_idx] /
          ((pred_test + samp_mean_local) / samp_median_global)

        # Clamp sample negatives to raw before rescaling (original line 855)
        norm[!train_idx][norm[!train_idx] < 0] <-
          e_current_batch[j, !train_idx][norm[!train_idx] < 0]

        # Rescale to preserve global medians (original lines 861-862)
        norm[train_idx]  <- norm[train_idx]  /
          (median(norm[train_idx],  na.rm = TRUE) / qc_median_global)
        norm[!train_idx] <- norm[!train_idx] /
          (median(norm[!train_idx], na.rm = TRUE) / samp_median_global)

        # Replace non-finite with tiny noise (original line 863)
        norm[!is.finite(norm)] <- rnorm(
          sum(!is.finite(norm)),
          sd = sd(norm[is.finite(norm)], na.rm = TRUE) * 0.01
        )

        # Outlier fallback: additive correction for extreme outliers (original lines 868-870)
        out <- boxplot.stats(norm, coef = 3)$out
        norm[!train_idx][norm[!train_idx] %in% out] <-
          (e_current_batch[j, !train_idx] -
             (pred_test + samp_mean_local - samp_median_global))[
               norm[!train_idx] %in% out]
        norm[!train_idx][norm[!train_idx] < 0] <-
          e_current_batch[j, !train_idx][norm[!train_idx] < 0]

        normalized[batch. %in% b] <- norm
      }

      normalized
    },
    all, batch., sampleType., num, corrs_train, corrs_target
  )

  normed <- t(pred)

  # Post-process: replace any remaining NAs / negatives (original lines 920-935)
  fix_negatives <- function(m) {
    for (i in seq_len(nrow(m))) {
      ok <- !is.na(m[i, ])
      if (sum(ok) == 0) next
      mn_pos <- min(m[i, ok], na.rm = TRUE)
      sd_pos <- sd(m[i, ok], na.rm = TRUE) * 0.1
      na_idx <- which(is.na(m[i, ]))
      if (length(na_idx) > 0)
        m[i, na_idx] <- rnorm(length(na_idx), mean = mn_pos, sd = sd_pos)
      neg_idx <- which(m[i, ] < 0)
      if (length(neg_idx) > 0)
        m[i, neg_idx] <- runif(1) * min(m[i, m[i, ] > 0], na.rm = TRUE)
    }
    m
  }

  normed_train  <- fix_negatives(normed[, sampleType. == "qc",     drop = FALSE])
  normed_target <- fix_negatives(normed[, sampleType. == "sample",  drop = FALSE])

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
correct_serrf <- function(data, num = 10) {
  library(ranger)
  library(parallel)

  # ── 1. Pre-imputation detection mask (for unbiased ltQC metrics) ────────────
  obs_mask <- !is.na(assay(data, 1))

  # ── 2. Extract matrices and metadata ────────────────────────────────────────
  e_all <- assay(data, 1)
  cd    <- as.data.frame(colData(data))

  stype <- ifelse(cd$QC == "QC",     "qc",
           ifelse(cd$QC == "Sample",  "sample",
           ifelse(cd$QC == "ltQC",    "ltqc", NA_character_)))

  # ── 3. SERRF half-minimum imputation per batch (applied only to NA cells) ───
  batches <- as.character(cd$Batch)
  for (b in unique(batches)) {
    b_idx <- which(batches == b)
    for (i in seq_len(nrow(e_all))) {
      na_idx <- which(is.na(e_all[i, b_idx]))
      if (length(na_idx) == 0) next
      finite_min <- min(e_all[i, b_idx], na.rm = TRUE)
      if (!is.finite(finite_min)) finite_min <- 1
      e_all[i, b_idx[na_idx]] <- 0.5 * finite_min
    }
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

  # ── 5. Parallel cluster ──────────────────────────────────────────────────────
  n_cores_env <- Sys.getenv("N_CORES", unset = "")
  n_cores <- if (n_cores_env == "") parallel::detectCores() - 1L else as.integer(n_cores_env)
  n_cores <- max(1L, n_cores)
  message("==> SERRF: starting ", n_cores, " parallel worker(s)")
  cl <- parallel::makeCluster(n_cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterEvalQ(cl, library(ranger))

  # ── 6. Correct biological samples ───────────────────────────────────────────
  message("==> SERRF: correcting Sample (n=", length(samp_idx), ") using QC (n=", length(qc_idx), ")")

  idx_qs   <- c(qc_idx, samp_idx)
  stype_qs <- c(rep("qc", length(qc_idx)), rep("sample", length(samp_idx)))
  batch_qs <- factor(batch_all[idx_qs])

  res_samp <- serrfR(
    train       = e_qc,
    target      = e_samp,
    num         = num,
    batch.      = batch_qs,
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

    res_ltqc <- serrfR(
      train       = e_qc,
      target      = e_ltqc,
      num         = num,
      batch.      = batch_ql,
      sampleType. = stype_ql,
      cl          = cl
    )
    e_ltqc_corr <- res_ltqc$normed_target
  }

  # ── 8. Reassemble SE in original column order ─────────────────────────────
  combined <- data
  e_new    <- assay(combined, 1)

  e_new[, qc_idx]   <- e_qc_corr
  e_new[, samp_idx] <- e_samp_corr
  if (!is.null(e_ltqc_corr)) {
    e_new[, ltqc_idx] <- e_ltqc_corr
  }

  assay(combined, 1, withDimnames = FALSE) <- e_new

  message("==> Imputation (RF on corrected data)")
  combined <- rf_impute_corrected(combined, obs_mask)

  list(pre = combined, post = combined, obs_mask = obs_mask)
}

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

# Half-minimum (LoD/2) imputation.
# Fills each NA with half the minimum observed value for that feature.
# Used as a pre-correction fill to give correction methods a complete matrix.
lod2_impute <- function(se) {
  mat <- assay(se, 1)
  for (i in seq_len(nrow(mat))) {
    na_idx <- which(is.na(mat[i, ]))
    if (length(na_idx) == 0) next
    finite_min <- min(mat[i, ], na.rm = TRUE)
    if (!is.finite(finite_min)) finite_min <- 1
    mat[i, na_idx] <- 0.5 * finite_min
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

correct_none <- function(data) {
  message("==> No correction (imputation only)")
  obs_mask <- !is.na(assay(data, 1))
  combined <- impute_rf(data, parallelize = "variables")
  list(pre = combined, post = combined, obs_mask = obs_mask)
}

correct_notame <- function(data, ruv_k) {
  obs_mask <- !is.na(assay(data, 1))

  # Cubic spline handles NAs natively — no LoD/2 before this step
  message("==> Drift correction (notame cubic spline)")
  combined <- merge_notame_sets(
    lapply(split_by_batch(data), process_batch),
    merge = "samples"
  )

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

correct_loess_combat <- function(data, loess_span, fallback_to_samples = FALSE) {
  library(sva)

  obs_mask <- !is.na(assay(data, 1))

  # LOESS handles NAs natively via is.finite() — no LoD/2 before this step
  message("==> Drift correction (LOESS)")
  combined <- merge_notame_sets(
    lapply(split_by_batch(data), function(se_b) {
      loess_correct_batch(se_b, span = loess_span, fallback_to_samples = fallback_to_samples)
    }),
    merge = "samples"
  )

  # LoD/2 fill before ComBat which requires a complete matrix
  combined <- lod2_impute(combined)
  pre      <- combined

  message("==> Between-batch correction (ComBat)")
  assay(combined, 1, withDimnames = FALSE) <- ComBat(
    dat   = assay(combined, 1),
    batch = as.factor(colData(combined)$Batch)
  )

  message("==> Imputation (RF on corrected data)")
  combined <- rf_impute_corrected(combined, obs_mask)

  list(pre = pre, post = combined, obs_mask = obs_mask)
}

correct_loess_limma <- function(data, loess_span, fallback_to_samples = FALSE) {
  library(limma)

  obs_mask <- !is.na(assay(data, 1))

  # LOESS handles NAs natively via is.finite() — no LoD/2 before this step
  message("==> Drift correction (LOESS)")
  combined <- merge_notame_sets(
    lapply(split_by_batch(data), function(se_b) {
      loess_correct_batch(se_b, span = loess_span, fallback_to_samples = fallback_to_samples)
    }),
    merge = "samples"
  )

  # LoD/2 fill before limma which requires a complete matrix
  combined <- lod2_impute(combined)
  pre      <- combined

  n_batches <- length(unique(colData(combined)$Batch))
  if (n_batches < 2) {
    message("==> Batch correction skipped (only one batch detected)")
  } else {
    message("==> Between-batch correction (limma removeBatchEffect)")
    assay(combined, 1, withDimnames = FALSE) <- removeBatchEffect(
      x     = assay(combined, 1),
      batch = as.factor(colData(combined)$Batch)
    )
  }

  message("==> Imputation (RF on corrected data)")
  combined <- rf_impute_corrected(combined, obs_mask)

  list(pre = pre, post = combined, obs_mask = obs_mask)
}

correct_combat_only <- function(data) {
  library(sva)

  obs_mask <- !is.na(assay(data, 1))
  data     <- lod2_impute(data)
  pre      <- data

  message("==> Batch correction (ComBat only, no drift correction)")
  assay(data, 1, withDimnames = FALSE) <- ComBat(
    dat   = assay(data, 1),
    batch = as.factor(colData(data)$Batch)
  )

  message("==> Imputation (RF on corrected data)")
  combined <- rf_impute_corrected(data, obs_mask)

  list(pre = pre, post = combined, obs_mask = obs_mask)
}

correct_pmp_qcrsc <- function(data) {
  library(pmp)

  obs_mask <- !is.na(assay(data, 1))

  message("==> Drift correction + batch correction (pmp QC-RSC)")
  # QCRSC handles NAs natively — no LoD/2 before this step
  # ltQC remapped to "Sample" so pmp does not try to use it as a QC reference
  classes_for_pmp <- ifelse(colData(data)$QC == "QC", "QC", "Sample")

  combined <- QCRSC(
    df      = data,
    order   = colData(data)$Injection_order,
    batch   = colData(data)$Batch,
    classes = classes_for_pmp,
    spar    = 0,
    minQC   = 4
  )

  message("==> Imputation (RF on corrected data)")
  combined <- rf_impute_corrected(combined, obs_mask)

  list(pre = combined, post = combined, obs_mask = obs_mask)
}

correct_batchcorr <- function(data,
                              G          = seq(5, 35, by = 10),
                              modelNames = c("VVV", "VVE", "VEV", "VEE", "VEI", "VVI", "VII"),
                              qualRatio  = 0.4) {
  library(batchCorr)

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
    merged <- list(
      peakTableCorr = bc$peakTable,
      peakTableOrg  = mat[which(as.character(meta$Batch) == b)[order(meta[which(as.character(meta$Batch) == b), "Injection_order"])], , drop = FALSE]
    )
  } else {
    message("==> Merging batches")
    merged <- mergeBatches(batch_corrObjs, qualRatio = qualRatio)
  }

  kept_features <- colnames(merged$peakTableCorr)
  kept_samples  <- rownames(merged$peakTableCorr)

  pre <- data[kept_features, kept_samples]
  assay(pre, "abundances", withDimnames = FALSE) <- t(merged$peakTableOrg[kept_samples, ])

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

  combined <- pre
  assay(combined, "abundances", withDimnames = FALSE) <- t(norm_result$peakTable)

  obs_mask <- obs_mask[kept_features, kept_samples, drop = FALSE]

  message("==> Imputation (RF on corrected data)")
  combined <- rf_impute_corrected(combined, obs_mask)

  list(pre = pre, post = combined, obs_mask = obs_mask)
}

correct_waveica <- function(data) {
  library(WaveICA2.0)

  obs_mask <- !is.na(assay(data, 1))
  data     <- lod2_impute(data)

  message("==> WaveICA2.0 correction")
  corrected_mat <- WaveICA_2.0(
    data            = t(assay(data, 1)),
    wf              = "haar",
    Injection_Order = as.numeric(colData(data)$Injection_order),
    alpha           = 0.05,
    Cutoff          = 0.10,
    K               = length(unique(colData(data)$Batch)) * 2
  )
  assay(data, 1, withDimnames = FALSE) <- t(corrected_mat$data)

  message("==> Imputation (RF on corrected data)")
  combined <- rf_impute_corrected(data, obs_mask)

  list(pre = combined, post = combined, obs_mask = obs_mask)
}

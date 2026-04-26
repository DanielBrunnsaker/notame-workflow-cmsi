# ─────────────────────────────────────────────────────────────────────────────
# High-level batch correction method wrappers
# ─────────────────────────────────────────────────────────────────────────────

# Wrapper around impute_rf that captures a pre-imputation missingness mask.
# Returns list(data = imputed_SE, mask = logical_matrix) so the mask never
# touches the SE and cannot confuse notame functions expecting a single assay.
impute_with_mask <- function(data, ...) {
  mask    <- !is.na(assay(data, 1))
  imputed <- impute_rf(data, ...)
  list(data = imputed, mask = mask)
}

correct_none <- function(data) {
  message("==> No correction (imputation only)")
  imp      <- impute_with_mask(data, parallelize = "variables")
  combined <- imp$data
  obs_mask <- imp$mask
  list(pre = combined, post = combined, obs_mask = obs_mask)
}

correct_notame <- function(data, ruv_k) {
  message("==> Drift correction (notame cubic spline)")
  combined <- merge_notame_sets(
    lapply(split_by_batch(data), process_batch),
    merge = "samples"
  )

  message("==> Imputation")
  imp      <- impute_with_mask(combined, parallelize = "variables")
  combined <- imp$data
  obs_mask <- imp$mask
  pre      <- combined

  n_batches <- length(unique(colData(combined)$Batch))
  if (n_batches < 2) {
    message("==> Batch correction skipped (only one batch detected)")
  } else {
    message("==> Batch correction (RUV, k=", ruv_k, ")")
    qc_idx   <- which(colData(combined)$QC == "QC")
    combined <- ruvs_qc(combined, replicates = list(qc_idx), k = ruv_k)
  }

  list(pre = pre, post = combined, obs_mask = obs_mask)
}

correct_loess_combat <- function(data, loess_span, fallback_to_samples = FALSE) {
  library(sva)

  message("==> Drift correction (LOESS)")
  combined <- merge_notame_sets(
    lapply(split_by_batch(data), function(se_b) {
      loess_correct_batch(se_b, span = loess_span, fallback_to_samples = fallback_to_samples)
    }),
    merge = "samples"
  )

  message("==> Imputation")
  imp      <- impute_with_mask(combined, parallelize = "variables")
  combined <- imp$data
  obs_mask <- imp$mask
  pre      <- combined

  message("==> Between-batch correction (ComBat)")
  assay(combined, 1, withDimnames = FALSE) <- ComBat(
    dat   = assay(combined, 1),
    batch = as.factor(colData(combined)$Batch)
  )

  list(pre = pre, post = combined, obs_mask = obs_mask)
}



correct_loess_limma <- function(data, loess_span, fallback_to_samples = FALSE) {
  library(limma)

  message("==> Drift correction (LOESS)")
  combined <- merge_notame_sets(
    lapply(split_by_batch(data), function(se_b) {
      loess_correct_batch(se_b, span = loess_span, fallback_to_samples = fallback_to_samples)
    }),
    merge = "samples"
  )

  message("==> Imputation")
  imp      <- impute_with_mask(combined, parallelize = "variables")
  combined <- imp$data
  obs_mask <- imp$mask
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

  list(pre = pre, post = combined, obs_mask = obs_mask)
}


correct_combat_only <- function(data) {
  library(sva)

  message("==> Imputation")
  imp      <- impute_with_mask(data, parallelize = "variables")
  combined <- imp$data
  obs_mask <- imp$mask

  pre <- combined

  message("==> Batch correction (ComBat only, no drift correction)")
  assay(combined, 1, withDimnames = FALSE) <- ComBat(
    dat   = assay(combined, 1),
    batch = as.factor(colData(combined)$Batch)
  )

  list(pre = pre, post = combined, obs_mask = obs_mask)
}


correct_pmp_qcrsc <- function(data) {

  library(pmp)

  message("==> Drift correction + batch correction (pmp QC-RSC)")

  # ltQC remapped to "Sample" so pmp does not try to use it as a QC reference, because package stupid
  classes_for_pmp <- ifelse(colData(data)$QC == "QC", "QC", "Sample")

  combined <- QCRSC(
    df      = data,
    order   = colData(data)$Injection_order,
    batch   = colData(data)$Batch,
    classes = classes_for_pmp,
    spar    = 0,    # 0 = auto-select via cross-validation
    minQC   = 4
  )

  message("==> Imputation")
  imp      <- impute_with_mask(combined, parallelize = "variables")
  combined <- imp$data
  obs_mask <- imp$mask

  list(pre = combined, post = combined, obs_mask = obs_mask)
}

correct_batchcorr <- function(data,
                              G          = seq(5, 35, by = 10),
                              modelNames = c("VVV", "VVE", "VEV", "VEE", "VEI", "VVI", "VII"),
                              qualRatio  = 0.4) {
  library(batchCorr)

  message("==> Imputation (required by batchCorr)")
  imp      <- impute_with_mask(data, parallelize = "variables")
  data     <- imp$data
  obs_mask <- imp$mask

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
        CVlimit      = Inf,   # disable internal feature removal; pre-filtering already done
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
    b       <- names(batch_corrObjs)[1]
    bc      <- batch_corrObjs[[b]]
    merged  <- list(
      peakTableCorr = bc$peakTable,
      peakTableOrg  = mat[which(as.character(meta$Batch) == b)[order(meta[which(as.character(meta$Batch) == b), "Injection_order"])], , drop = FALSE]
    )
  } else {
    message("==> Merging batches")
    merged <- mergeBatches(batch_corrObjs, qualRatio = qualRatio)
  }

  # Reconstruct SE from merged matrices (samples × features → features × samples)
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

  # Subset mask to match features/samples kept by mergeBatches
  obs_mask <- obs_mask[kept_features, kept_samples, drop = FALSE]

  list(pre = pre, post = combined, obs_mask = obs_mask)
}


correct_waveica <- function(data) {
  # Install once with: remotes::install_github("dengkuistat/WaveICA2.0")
  library(WaveICA2.0)

  message("==> Imputation (pre-WaveICA, required by algorithm)")
  imp      <- impute_with_mask(data, parallelize = "variables")
  data     <- imp$data
  obs_mask <- imp$mask

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
  list(pre = data, post = data, obs_mask = obs_mask)
}

# ─────────────────────────────────────────────────────────────────────────────
# High-level batch correction method wrappers
# ─────────────────────────────────────────────────────────────────────────────

# Wrapper around impute_rf that stores a pre-imputation missingness mask as a
# second assay ("observed": TRUE = originally detected, FALSE = was NA/missing).
# Used by eval_ltqc / eval_ltqc_dratio to avoid artefacts from imputed values.
impute_with_mask <- function(data, ...) {
  mask    <- !is.na(assay(data, 1))
  imputed <- impute_rf(data, ...)
  assay(imputed, "observed", withDimnames = FALSE) <- mask
  imputed
}

correct_notame <- function(data, ruv_k) {
  message("==> Drift correction (notame cubic spline)")
  combined <- merge_notame_sets(
    lapply(split_by_batch(data), process_batch),
    merge = "samples"
  )

  message("==> Imputation")
  combined <- impute_with_mask(combined, parallelize = "variables")
  pre <- combined

  message("==> Batch correction (RUV, k=", ruv_k, ")")
  qc_idx   <- which(colData(combined)$QC == "QC")
  combined <- ruvs_qc(combined, replicates = list(qc_idx), k = ruv_k)

  list(pre = pre, post = combined)
}

correct_loess_combat <- function(data, loess_span) {
  library(sva)

  message("==> Drift correction (LOESS)")
  combined <- merge_notame_sets(
    lapply(split_by_batch(data), function(se_b) {
      message("  Batch ", unique(colData(se_b)$Batch))
      loess_correct_batch(se_b, span = loess_span)
    }),
    merge = "samples"
  )

  message("==> Imputation")
  combined <- impute_with_mask(combined, parallelize = "variables")
  pre <- combined

  message("==> Between-batch correction (ComBat)")
  assay(combined, 1, withDimnames = FALSE) <- ComBat(
    dat   = assay(combined, 1),
    batch = as.factor(colData(combined)$Batch)
  )

  list(pre = pre, post = combined)
}


correct_loess_limma <- function(data, loess_span) {
  library(limma)

  message("==> Drift correction (LOESS)")
  combined <- merge_notame_sets(
    lapply(split_by_batch(data), function(se_b) {
      message("  Batch ", unique(colData(se_b)$Batch))
      loess_correct_batch(se_b, span = loess_span)
    }),
    merge = "samples"
  )

  message("==> Imputation")
  combined <- impute_with_mask(combined, parallelize = "variables")
  pre <- combined

  message("==> Between-batch correction (limma::removeBatchEffect)")
  assay(combined, 1, withDimnames = FALSE) <- removeBatchEffect(
    x     = assay(combined, 1),
    batch = as.factor(colData(combined)$Batch)
  )

  list(pre = pre, post = combined)
}


correct_linear_combat <- function(data) {
  library(sva)

  message("==> Drift correction (linear/ANCOVA)")
  combined <- merge_notame_sets(
    lapply(split_by_batch(data), function(se_b) {
      message("  Batch ", unique(colData(se_b)$Batch))
      linear_correct_batch(se_b)
    }),
    merge = "samples"
  )

  message("==> Imputation")
  combined <- impute_with_mask(combined, parallelize = "variables")
  pre <- combined

  message("==> Between-batch correction (ComBat)")
  assay(combined, 1, withDimnames = FALSE) <- ComBat(
    dat   = assay(combined, 1),
    batch = as.factor(colData(combined)$Batch)
  )

  list(pre = pre, post = combined)
}

correct_linear_limma <- function(data) {
  library(limma)

  message("==> Drift correction (linear/ANCOVA)")
  combined <- merge_notame_sets(
    lapply(split_by_batch(data), function(se_b) {
      message("  Batch ", unique(colData(se_b)$Batch))
      linear_correct_batch(se_b)
    }),
    merge = "samples"
  )

  message("==> Imputation")
  combined <- impute_with_mask(combined, parallelize = "variables")
  pre <- combined

  message("==> Between-batch correction (limma::removeBatchEffect)")
  assay(combined, 1, withDimnames = FALSE) <- removeBatchEffect(
    x     = assay(combined, 1),
    batch = as.factor(colData(combined)$Batch)
  )

  list(pre = pre, post = combined)
}


correct_combat_only <- function(data) {
  library(sva)

  message("==> Imputation")
  combined <- impute_with_mask(data, parallelize = "variables")

  message("==> Batch correction (ComBat only, no drift correction)")
  assay(combined, 1, withDimnames = FALSE) <- ComBat(
    dat   = assay(combined, 1),
    batch = as.factor(colData(combined)$Batch)
  )

  list(pre = combined, post = combined)
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
  combined <- impute_with_mask(combined, parallelize = "variables")

  list(pre = combined, post = combined)
}

correct_batchcorr <- function(data,
                              G          = seq(5, 35, by = 10),
                              modelNames = c("VVV", "VVE", "VEV", "VEE", "VEI", "VVI", "VII"),
                              CVlimit    = 0.3,
                              qualRatio  = 0.4) {
  library(batchCorr)

  message("==> Imputation (required by batchCorr)")
  data <- impute_with_mask(data, parallelize = "variables")

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

  message("==> Merging batches")
  merged <- mergeBatches(batch_corrObjs, qualRatio = qualRatio)

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
    CVlimit       = CVlimit
  )

  combined <- pre
  assay(combined, "abundances", withDimnames = FALSE) <- t(norm_result$peakTable)

  list(pre = pre, post = combined)
}


correct_waveica <- function(data) {
  # Install once with: remotes::install_github("dengkuistat/WaveICA2.0")
  library(WaveICA2.0)

  message("==> Imputation (pre-WaveICA, required by algorithm)")
  data <- impute_with_mask(data, parallelize = "variables")

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
  list(pre = data, post = data)
}

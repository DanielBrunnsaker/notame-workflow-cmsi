# ─────────────────────────────────────────────────────────────────────────────
# High-level batch correction method wrappers
# ─────────────────────────────────────────────────────────────────────────────

correct_notame <- function(data, ruv_k) {
  message("==> Drift correction (notame cubic spline)")
  combined <- merge_notame_sets(
    lapply(split_by_batch(data), process_batch),
    merge = "samples"
  )

  message("==> Imputation")
  combined <- impute_rf(combined)
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
  combined <- impute_rf(combined)
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
  combined <- impute_rf(combined)
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
  combined <- impute_rf(combined)
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
  combined <- impute_rf(combined)
  pre <- combined

  message("==> Between-batch correction (limma::removeBatchEffect)")
  assay(combined, 1, withDimnames = FALSE) <- removeBatchEffect(
    x     = assay(combined, 1),
    batch = as.factor(colData(combined)$Batch)
  )

  list(pre = pre, post = combined)
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
  combined <- impute_rf(combined)

  list(pre = combined, post = combined)
}

correct_waveica <- function(data) {
  # Install once with: remotes::install_github("dengkuistat/WaveICA2.0")
  library(WaveICA2.0)

  message("==> Imputation (pre-WaveICA, required by algorithm)")
  data <- impute_rf(data)

  message("==> WaveICA2.0 correction")
  corrected_mat <- WaveICA_2.0(
    data            = t(assay(data, 1)),
    wf              = "haar",
    Injection_Order = as.numeric(colData(data)$Injection_order),
    alpha           = 0.05,
    Cutoff          = 0.10,
    K               = 20  # Default parameter, but this likely needs to be chagned for this method to even be viable
  )
  assay(data, 1, withDimnames = FALSE) <- t(corrected_mat$data)
  list(pre = data, post = data)
}

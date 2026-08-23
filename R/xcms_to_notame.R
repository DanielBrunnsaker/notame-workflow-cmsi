# ─────────────────────────────────────────────────────────────────────────────
# xcms_to_notame()
#
# Converts an XCMS-based feature table + sample sheet to the Excel format
# expected by notame::import_from_excel() — an alternative to
# msdial_to_notame() for pipelines that already produce structured sample
# metadata instead of embedding it in filenames.
#
# Inputs:
#   feature_table_csv — CSV with columns: feature, mzmed, mzmin, mzmax, rtmed,
#                        rtmin, rtmax, npeaks, <per-type detection counts>,
#                        ms_level, then one abundance column per sample,
#                        HEADED BY sample_label (the join key into the sheet).
#   sample_sheet_xlsx  — XLSX with (at least): batch, column, polarity,
#                        sample_label, sample_type, injection_order, filename,
#                        include.
#
# Derived metadata (mirrors msdial_to_notame(), see that file for the
# filename-parsing equivalent used on the MSDIAL side):
#   Batch           <- sample sheet's `batch` column, used as-is
#   Injection_order <- global rank of (batch, injection_order)
#   QC              <- sample_type mapped via DEFAULT_XCMS_SAMPLE_TYPE_MAP,
#                      falling back to "Sample" for anything unrecognized
#   Original_name   <- "{batch}_{sample_label}" (uniqueness across batches)
#   <mode>_Datafile <- sample sheet's `filename` column
# ─────────────────────────────────────────────────────────────────────────────

DEFAULT_XCMS_SAMPLE_TYPE_MAP <- c(
  "Sample"      = "Sample",
  "sQC"         = "QC",
  "QC"          = "QC",
  "ltQC"        = "ltQC",
  "Blank"       = "Blank",
  "MatrixBlank" = "MatrixBlank",
  "Wash"        = "Wash",
  "Cond"        = "Cond",
  "SST"         = "SST"
)

xcms_to_notame <- function(feature_table_csv, sample_sheet_xlsx, out_xlsx,
                            column, polarity, sample_type_map = NULL) {
  type_map <- if (is.null(sample_type_map)) DEFAULT_XCMS_SAMPLE_TYPE_MAP else sample_type_map

  sheet <- read.xlsx(sample_sheet_xlsx, sheet = 1)
  ft    <- read.csv(feature_table_csv, check.names = FALSE, stringsAsFactors = FALSE)

  required_sheet_cols <- c("batch", "column", "polarity", "sample_label", "sample_type",
                            "injection_order", "filename", "include")
  missing_cols <- setdiff(required_sheet_cols, colnames(sheet))
  if (length(missing_cols) > 0)
    stop("sample_sheet is missing required column(s): ", paste(missing_cols, collapse = ", "))

  # Filter to included rows matching the requested COLUMN/POLARITY
  keep <- as.logical(sheet$include) &
    toupper(trimws(sheet$column))   == toupper(column) &
    toupper(trimws(sheet$polarity)) == toupper(polarity)
  keep[is.na(keep)] <- FALSE

  if (!any(keep)) {
    available <- unique(paste0(sheet$column, "_", sheet$polarity))
    stop("No included sample_sheet rows match COLUMN=", column, ", POLARITY=", polarity,
         ". Combinations present in the sheet: ", paste(available, collapse = ", "))
  }
  sheet <- sheet[keep, , drop = FALSE]

  needs_review <- if ("needs_review" %in% colnames(sheet)) as.logical(sheet$needs_review) else logical(nrow(sheet))
  if (any(needs_review, na.rm = TRUE))
    warning("sample_sheet has ", sum(needs_review, na.rm = TRUE),
            " included row(s) flagged needs_review == TRUE: ",
            paste(sheet$sample_label[which(needs_review)], collapse = ", "), call. = FALSE)

  # Identify sample columns in the feature table by matching sample_label —
  # unlike MSDIAL there's no fixed filename convention to pattern-match, so
  # the sheet's sample_label is the join key.
  sample_cols <- intersect(colnames(ft), sheet$sample_label)
  if (length(sample_cols) == 0)
    stop("No feature_table columns match any sample_sheet$sample_label for COLUMN=",
         column, ", POLARITY=", polarity)

  sheet <- sheet[match(sample_cols, sheet$sample_label), , drop = FALSE]  # align to sample_cols order

  # Global run order: rank by (batch, injection_order), same approach as
  # msdial_to_notame() — robust whether injection_order is locally- or
  # globally-scoped in the sheet.
  injection_order  <- suppressWarnings(as.numeric(sheet$injection_order))
  ord              <- order(sheet$batch, injection_order)
  global_run_order <- integer(nrow(sheet))
  global_run_order[ord] <- seq_along(ord)

  qc_type       <- ifelse(sheet$sample_type %in% names(type_map), type_map[sheet$sample_type], "Sample")
  original_name <- paste0(sheet$batch, "_", sheet$sample_label)

  mode_name <- paste0(column, "_", tolower(polarity))
  col_ids   <- sprintf("%s_%03d", mode_name, seq_along(sample_cols))

  # Assemble notame output layout (shared contract — see R/notame_format.R)
  meta_block <- rbind(
    notame_meta_row("Sample_ID",                    col_ids),
    notame_meta_row("Injection_order",              global_run_order),
    notame_meta_row("QC",                           qc_type),
    notame_meta_row("Batch",                        sheet$batch),
    notame_meta_row("Original_name",               original_name),
    notame_meta_row(paste0(mode_name, "_Datafile"), sheet$filename)
  )

  # Per-feature metadata
  feature_id_col <- as.character(ft$feature)
  alignment_ids  <- suppressWarnings(as.integer(gsub("\\D", "", feature_id_col)))
  feature_ids    <- gsub("\\.", "_",
    sprintf("%s_%s_%.4f_%.4f", mode_name, feature_id_col, ft$mzmed, ft$rtmed))

  # Fill_pct: XCMS's npeaks (features actually detected) is directly
  # analogous to MSDIAL's "Fill %" numerator — no equivalent to
  # Adduct_type/Metabolite_name since this pipeline doesn't annotate.
  n_samples_total <- length(sample_cols)
  fill_pct <- if ("npeaks" %in% colnames(ft)) round(ft$npeaks / n_samples_total * 100, 1) else NA_real_

  abund_mat <- ft[, sample_cols, drop = FALSE]
  abund_mat[] <- lapply(abund_mat, function(col) {
    vals <- suppressWarnings(as.numeric(col))
    vals[!is.na(vals) & vals < 0] <- 0
    as.character(vals)
  })

  feat_data_mat <- cbind(
    feature_ids,
    mode_name,
    as.character(alignment_ids),
    as.character(ft$mzmed),
    as.character(ft$rtmed),
    column,
    tolower(polarity),
    NA_character_,   # Adduct_type — not produced by this pipeline
    NA_character_,   # Metabolite_name — not produced by this pipeline
    as.character(fill_pct),
    NA_character_,   # Flag
    as.matrix(abund_mat)
  )

  write_notame_format(meta_block, feat_data_mat, col_ids, out_xlsx)
  message("(mode: ", mode_name, ", ", nrow(ft), " features, ", length(sample_cols), " samples)")

  # Build annotations table (the feature table's own metadata columns, keyed
  # by Feature_ID) in the same list(mode=, annotations=) shape
  # msdial_to_notame() returns.
  annot_cols <- intersect(
    c("feature", "mzmed", "mzmin", "mzmax", "rtmed", "rtmin", "rtmax", "npeaks",
      "Blank", "ltQC", "Sample", "sQC", "QC", "ms_level"),
    colnames(ft)
  )
  annot_df <- cbind(Feature_ID = feature_ids, ft[, annot_cols, drop = FALSE])

  invisible(list(mode = mode_name, annotations = annot_df))
}

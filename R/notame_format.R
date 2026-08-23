# ─────────────────────────────────────────────────────────────────────────────
# Shared notame-format writer
#
# Both msdial_to_notame() and xcms_to_notame() must produce byte-for-byte the
# same layout that notame::import_from_excel() expects: an 11-column feature
# metadata header, a 6-row sample metadata block, then one row per feature.
# Kept in one place so the two converters can't silently drift apart.
# ─────────────────────────────────────────────────────────────────────────────

NOTAME_FEAT_HEADER <- c("Feature_ID", "Split", "Alignment", "Average_Mz",
                         "Average_Rt_min", "Column", "Ion_mode", "Adduct_type",
                         "Metabolite_name", "Fill_pct", "Flag")
NOTAME_N_FEAT <- length(NOTAME_FEAT_HEADER)  # number of feature metadata columns

# One row of the sample metadata block: N_FEAT-1 NA padding columns, then the
# row label in the last metadata column, then one value per sample.
notame_meta_row <- function(label, values)
  c(rep(NA_character_, NOTAME_N_FEAT - 1), label, as.character(values))

# Assembles and writes the final notame-format workbook.
#   meta_block    — rbind of notame_meta_row() calls (sample metadata)
#   feat_data_mat — matrix/data.frame: NOTAME_N_FEAT feature-metadata columns
#                   (in NOTAME_FEAT_HEADER order) followed by the abundance
#                   columns, one row per feature
#   col_ids       — sample column ids, in the same order as feat_data_mat's
#                   abundance columns and meta_block's sample values
write_notame_format <- function(meta_block, feat_data_mat, col_ids, out_xlsx) {
  feat_header <- c(NOTAME_FEAT_HEADER, col_ids)
  out_mat     <- rbind(meta_block, feat_header, feat_data_mat)
  out_df      <- as.data.frame(out_mat, stringsAsFactors = FALSE)

  write.xlsx(out_df, out_xlsx, colNames = FALSE, rowNames = FALSE)
  message("notame-ready file written: ", out_xlsx)
  invisible(out_xlsx)
}

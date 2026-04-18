# ─────────────────────────────────────────────────────────────────────────────
# msdial_to_notame()
#
# Converts an MSDIAL alignment export to the Excel format expected by
# notame::import_from_excel().
#
# MSDIAL filename pattern (embedded in column headers):
#   YYYY-MM-DD_BATCH_COLUMN_POLARITY_NAME_INJORDER[SUFFIX]
#
# Derived metadata:
#   Sample_ID       <- NAME part of filename
#   Injection_order <- leading digits of the last underscore-delimited part
#   QC              <- "QC" if NAME contains "sQC", "Blank" if contains "blank",
#                      else "Sample", ... 
#   Batch           <- first digit(s) after "B" in BATCH  (B1-2W4 -> 1, B3W5 -> 3) - make this more reliable!
#   <mode>_Datafile <- full original MSDIAL filename
# ─────────────────────────────────────────────────────────────────────────────
msdial_to_notame <- function(in_xlsx, out_xlsx) {

  # Read raw file (no headers, everything as strings) 
  raw <- read.xlsx(in_xlsx, sheet = 1, colNames = FALSE, rowNames = FALSE)

  # Locate the MSDIAL header row ("Alignment ID" in col 1) 
  hdr_idx <- which(sapply(raw[[1]], function(x)
    trimws(as.character(x)) == "Alignment ID"))
  if (length(hdr_idx) != 1)
    stop("Expected exactly one 'Alignment ID' row; found ", length(hdr_idx))
  hdr <- as.character(unlist(raw[hdr_idx, ]))

  # Identify sample columns by date-prefix pattern 
  is_sample  <- grepl("^\\d{4}-\\d{2}-\\d{2}_", hdr)
  sample_idx <- which(is_sample)
  if (length(sample_idx) == 0)
    stop("No sample columns found. Expected column headers starting with YYYY-MM-DD_")
  sample_fns <- hdr[sample_idx]

  # Parse each filename into metadata fields ───────────────────────────
  parse_fn <- function(fn) {
    # Capture: date, batch, everything between batch and POS/NEG, polarity
    prefix_pat <- "^(\\d{4}-\\d{2}-\\d{2})_([^_]+)_(.+?)_(POS|NEG)_"
    m_pre <- regmatches(fn, regexec(prefix_pat, fn, perl = TRUE))[[1]]
    if (length(m_pre) < 5) {
      warning("Could not parse filename: ", fn)
      return(list(batch = NA, column = NA, polarity = NA,
                  is_mse = FALSE, name = fn, injection_order = NA_integer_))
    }

    middle   <- m_pre[4]  # e.g. "RP", "RP_MSe", "RP_MSe_MSe"
    is_mse   <- grepl("MSe", middle, ignore.case = TRUE)
    is_ms2   <- grepl("MS2", middle, ignore.case = TRUE)
    col_type <- strsplit(middle, "_")[[1]][1]  # first token is always column type i think

    remainder <- sub(prefix_pat, "", fn)  # NAME_..._INJORDER[SUFFIX]
    tokens    <- strsplit(remainder, "_")[[1]]
    last      <- tail(tokens, 1)
    inj_digits <- regmatches(last, regexpr("^\\d+", last))

    if (length(inj_digits) > 0 && nchar(inj_digits) > 0) {
      inj_order <- as.integer(inj_digits)
      name      <- paste(tokens[-length(tokens)], collapse = "_")
    } else {
      # Pure-letter last token — injection order is in second-to-last, not always true though
      prev <- tokens[length(tokens) - 1]
      inj_digits2 <- regmatches(prev, regexpr("^\\d+", prev))
      if (length(inj_digits2) > 0 && nchar(inj_digits2) > 0) {
        inj_order <- as.integer(inj_digits2)
        name      <- paste(tokens[seq_len(length(tokens) - 2)], collapse = "_")
      } else {
        warning("Could not extract injection order from: ", fn)
        inj_order <- NA_integer_
        name      <- remainder
      }
    }

    list(batch           = m_pre[3],
         column          = col_type,
         polarity        = m_pre[5],
         is_mse          = is_mse,
         is_ms2          = is_ms2,
         name            = name,
         injection_order = inj_order)
  }

  parsed <- lapply(sample_fns, parse_fn)

  # MSe/MS2 detected from the middle section (between batch and POS/NEG),
  # not from the sample name. All other types detected from sample name.

  # Make this a bit more reliable if there is intention to automate this?
  classify_qc <- function(name, is_mse = FALSE, is_ms2 = FALSE) {
    if (is_mse)                                         return("MSe")
    if (is_ms2)                                         return("MS2")
    if (grepl("SST\\d",     name, ignore.case = TRUE))  return("SST")
    if (grepl("ltQC",       name, ignore.case = TRUE))  return("ltQC")
    if (grepl("sQC",        name, ignore.case = TRUE))  return("QC")
    if (grepl("MeOH",       name, ignore.case = TRUE))  return("Wash")
    if (grepl("blank",      name, ignore.case = TRUE))  return("Blank")
    if (grepl("CondPlasma", name, ignore.case = TRUE))  return("Cond")
    "Sample"
  }

  extract_batch <- function(batch_str) {
    m <- regmatches(batch_str, regexpr("(?<=B)\\d+", batch_str, perl = TRUE))
    if (length(m) == 0) return(NA_integer_)
    as.integer(m)
  }

  sample_id  <- sapply(parsed, `[[`, "name")
  inj_order  <- sapply(parsed, `[[`, "injection_order")
  qc         <- sapply(parsed, function(p) classify_qc(p$name, p$is_mse, p$is_ms2))
  batch      <- sapply(parsed, function(p) extract_batch(p$batch))

  # Unique sample IDs derived from the filename.
  # Batch prefix ensures QC sample names (e.g. sQC01) are unique across batches.
  sample_id_unique <- paste0("B", batch, "_", sample_id)

  # Global run order: rank samples by (batch, injection_order) to get a single
  # continuous sequence across all batches for drift correction.
  ord              <- order(batch, inj_order)
  global_run_order <- integer(length(ord))
  global_run_order[ord] <- seq_along(ord)

  # Derive mode from polarity and column
  non_mse_parsed <- parsed[!sapply(parsed, `[[`, "is_mse")]
  col_type  <- unique(na.omit(sapply(non_mse_parsed, `[[`, "column")))[1]
  pol_lower <- unique(na.omit(sapply(non_mse_parsed, function(p) tolower(p$polarity))))[1]
  mode_name <- paste0(col_type, "_", pol_lower)  # e.g. "RP_pos"

  # Stable column IDs for the data matrix
  col_ids <- sprintf("%s_%03d", mode_name, seq_along(sample_fns))

  # Extract and clean feature data 
  feat_rows <- raw[(hdr_idx + 1):nrow(raw), , drop = FALSE]

  parse_num <- function(x) as.numeric(gsub(",", ".", as.character(x)))

  align_col  <- which(hdr == "Alignment ID")
  rt_col     <- which(hdr == "Average Rt(min)")
  mz_col     <- which(hdr == "Average Mz")
  adduct_col <- which(hdr == "Adduct type")
  fill_col   <- which(hdr == "Fill %")

  alignment_ids <- as.integer(feat_rows[[align_col]])
  mz_vals       <- parse_num(feat_rows[[mz_col]])
  rt_vals       <- parse_num(feat_rows[[rt_col]])
  adduct_vals   <- if (length(adduct_col) == 1)
    as.character(feat_rows[[adduct_col]]) else NA_character_
  fill_vals     <- if (length(fill_col) == 1)
    parse_num(feat_rows[[fill_col]]) else NA_real_

  feature_ids <- gsub("\\.", "_",
    sprintf("%s_%04d_%.4f_%.4f", mode_name, alignment_ids, mz_vals, rt_vals))

  abund_mat <- feat_rows[, sample_idx, drop = FALSE]
  abund_mat[] <- lapply(abund_mat, function(col) {
    vals <- as.numeric(gsub(",", ".", as.character(col)))
    vals[!is.na(vals) & vals < 0] <- 0
    as.character(vals)
  })

  # Assemble notame output layout 
  N_FEAT <- 10 # number of feature metadata columns (must match feat_header length)

  meta_row <- function(label, values)
    c(rep(NA_character_, N_FEAT - 1), label, as.character(values))

  meta_block <- rbind(
    meta_row("Sample_ID",                    col_ids),
    meta_row("Injection_order",              global_run_order),
    meta_row("QC",                           qc),
    meta_row("Batch",                        batch),
    meta_row("Original_name",               sample_id_unique),
    meta_row(paste0(mode_name, "_Datafile"), sample_fns)
  )

  feat_header <- c("Feature_ID", "Split", "Alignment", "Average_Mz",
                   "Average_Rt_min", "Column", "Ion_mode", "Adduct_type", "Fill_pct", "Flag",
                   col_ids)

  feat_data_mat <- cbind(
    feature_ids,
    mode_name,
    as.character(alignment_ids),
    as.character(mz_vals),
    as.character(rt_vals),
    col_type,
    pol_lower,
    adduct_vals,
    as.character(fill_vals),
    NA_character_,
    as.matrix(abund_mat)
  )

  out_mat <- rbind(meta_block, feat_header, feat_data_mat)
  out_df  <- as.data.frame(out_mat, stringsAsFactors = FALSE)

  # Write 
  write.xlsx(out_df, out_xlsx, colNames = FALSE, rowNames = FALSE)
  message("notame-ready file written: ", out_xlsx, " (mode: ", mode_name, ")")
  invisible(mode_name)
}

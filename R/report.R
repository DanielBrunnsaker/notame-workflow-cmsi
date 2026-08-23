# ─────────────────────────────────────────────────────────────────────────────
# Run report
#
# Builds a single self-contained HTML file (output_dir/report.html) summarising
# one pipeline run: run parameters, filter counts, batch quality, per-method
# success/failure, the cross-method QC comparison, and the QC plots
# themselves — so a run can be assessed without digging through per-method
# folders. Inline CSS, PNGs embedded as base64 data URIs (jsonlite::base64_enc)
# — no external assets, no new heavy dependencies (no rmarkdown/pandoc).
# ─────────────────────────────────────────────────────────────────────────────

html_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;",  x, fixed = TRUE)
  x <- gsub(">", "&gt;",  x, fixed = TRUE)
  x
}

# Renders a data.frame as an HTML table. `row_class` is an optional character
# vector (same length as nrow(df)) of CSS classes applied per <tr> (empty
# string = no class), used to highlight e.g. failed runs or the best-ranked
# method. Columns are formatted individually (not coerced via as.character()
# on a data.frame row, which mishandles factor columns).
df_to_html_table <- function(df, row_class = NULL) {
  if (is.null(df) || nrow(df) == 0) return("<p class=\"muted\">(no data)</p>")

  cols <- lapply(df, function(col) html_escape(format(col, trim = TRUE)))
  n    <- nrow(df)

  header <- paste0("<th>", html_escape(colnames(df)), "</th>", collapse = "")
  rows <- vapply(seq_len(n), function(i) {
    cls   <- if (!is.null(row_class) && nzchar(row_class[i])) paste0(" class=\"", row_class[i], "\"") else ""
    cells <- paste0("<td>", vapply(cols, `[`, character(1), i), "</td>", collapse = "")
    paste0("<tr", cls, ">", cells, "</tr>")
  }, character(1))

  paste0("<table><thead><tr>", header, "</tr></thead><tbody>",
         paste(rows, collapse = ""), "</tbody></table>")
}

png_to_data_uri <- function(path) {
  raw <- readBin(path, "raw", file.info(path)$size)
  paste0("data:image/png;base64,", jsonlite::base64_enc(raw))
}

# One <figure> per PNG found directly inside `dir` (non-recursive).
plots_html <- function(dir) {
  if (!dir.exists(dir)) return("<p class=\"muted\">(no plots)</p>")
  files <- list.files(dir, pattern = "\\.png$", full.names = TRUE)
  if (length(files) == 0) return("<p class=\"muted\">(no plots)</p>")

  figs <- vapply(files, function(f) {
    label <- html_escape(tools::file_path_sans_ext(basename(f)))
    uri   <- tryCatch(png_to_data_uri(f), error = function(e) NA_character_)
    if (is.na(uri)) return("")
    sprintf("<figure><img src=\"%s\" alt=\"%s\"><figcaption>%s</figcaption></figure>",
            uri, label, label)
  }, character(1))

  paste0("<div class=\"plot-grid\">", paste(figs, collapse = ""), "</div>")
}

read_csv_safe <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch(read.csv(path, check.names = FALSE), error = function(e) NULL)
}

REPORT_CSS <- "
  :root { color-scheme: light dark; }
  body { font-family: -apple-system, Segoe UI, Helvetica, Arial, sans-serif;
         max-width: 1100px; margin: 2rem auto; padding: 0 1rem;
         color: #1a1a1a; background: #fff; line-height: 1.4; }
  h1 { font-size: 1.5rem; margin-bottom: 0.2rem; }
  h2 { font-size: 1.15rem; margin-top: 2.2rem; border-bottom: 1px solid #ddd; padding-bottom: 0.3rem; }
  .muted { color: #777; font-size: 0.9rem; }
  .meta  { color: #555; font-size: 0.9rem; margin-bottom: 1rem; }
  table { border-collapse: collapse; width: 100%; font-size: 0.85rem; margin: 0.5rem 0 1rem; }
  th, td { border: 1px solid #ddd; padding: 4px 8px; text-align: left; white-space: nowrap; }
  th { background: #f5f5f5; position: sticky; top: 0; }
  tr.status-failed  { background: #fdeaea; }
  tr.status-success { background: #eafaf0; }
  tr.best            { background: #eaf3fd; font-weight: 600; }
  .table-scroll { overflow-x: auto; }
  details { margin: 0.6rem 0; }
  summary { cursor: pointer; font-weight: 600; padding: 0.3rem 0; }
  .plot-grid { display: flex; flex-wrap: wrap; gap: 1rem; margin: 0.5rem 0 1rem; }
  figure { margin: 0; max-width: 320px; }
  figure img { max-width: 100%; border: 1px solid #ddd; border-radius: 4px; }
  figcaption { font-size: 0.75rem; color: #666; margin-top: 0.2rem; word-break: break-all; }
  @media (prefers-color-scheme: dark) {
    body { color: #e8e8e8; background: #1c1c1e; }
    h2 { border-bottom-color: #444; }
    th { background: #2a2a2c; }
    th, td { border-color: #3a3a3c; }
    tr.status-failed  { background: #3a2222; }
    tr.status-success { background: #1f3324; }
    tr.best            { background: #1e2c3d; }
    figure img { border-color: #444; }
  }
"

# Builds output_dir/report.html summarising one pipeline run. Call this
# wrapped in tryCatch at the call site so a report failure never fails an
# otherwise-successful run.
#
#   run_params_df — Parameter/Value data.frame (as written to run_parameters.txt)
#   run_log_df    — combined data.frame returned by write_run_log(), or NULL
write_run_report <- function(interdir, output_dir, mode_label, in_xlsx,
                              run_params_df, run_log_df) {
  prefilter  <- read_csv_safe(file.path(interdir, "prefilter_log.csv"))
  batch_sum  <- read_csv_safe(file.path(interdir, "batch_summary.csv"))
  comparison <- read_csv_safe(file.path(output_dir, "method_comparison.csv"))

  run_log_rows <- if (!is.null(run_log_df)) paste0("status-", tolower(run_log_df$status)) else NULL

  comparison_rows <- NULL
  if (!is.null(comparison) && nrow(comparison) > 0) {
    comparison_rows <- rep("", nrow(comparison))
    ranked <- which(comparison$method != "uncorrected")
    if (length(ranked) > 0) comparison_rows[ranked[1]] <- "best"
  }

  methods <- setdiff(list.dirs(output_dir, recursive = FALSE, full.names = FALSE), "pre_correction")
  method_sections <- paste(vapply(methods, function(m) {
    sprintf("<details><summary>%s</summary>%s</details>",
            html_escape(m), plots_html(file.path(output_dir, m, "QC_plots")))
  }, character(1)), collapse = "")

  body <- paste0(
    "<h1>", html_escape(mode_label), " — run report</h1>",
    "<p class=\"meta\">Generated ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    " &middot; input: ", html_escape(basename(in_xlsx)), "</p>",

    "<h2>Run parameters</h2>",
    "<div class=\"table-scroll\">", df_to_html_table(run_params_df), "</div>",

    "<h2>Pre-correction feature filtering</h2>",
    "<div class=\"table-scroll\">", df_to_html_table(prefilter), "</div>",

    "<h2>Batch quality (pre-correction)</h2>",
    "<div class=\"table-scroll\">", df_to_html_table(batch_sum), "</div>",

    "<h2>Run log</h2>",
    "<div class=\"table-scroll\">", df_to_html_table(run_log_df, row_class = run_log_rows), "</div>",

    "<h2>Method comparison</h2>",
    "<p class=\"muted\">Best (non-baseline) method by rank highlighted.</p>",
    "<div class=\"table-scroll\">", df_to_html_table(comparison, row_class = comparison_rows), "</div>",

    "<h2>Pre-correction QC plots</h2>",
    plots_html(file.path(output_dir, "pre_correction")),

    "<h2>Per-method QC plots</h2>",
    if (nzchar(method_sections)) method_sections else "<p class=\"muted\">(no methods completed)</p>"
  )

  html <- paste0(
    "<!doctype html><html><head><meta charset=\"utf-8\">",
    "<title>", html_escape(mode_label), " run report</title>",
    "<style>", REPORT_CSS, "</style></head><body>", body, "</body></html>"
  )

  out_file <- file.path(output_dir, "report.html")
  writeLines(html, out_file)
  message("==> Report written: ", out_file)
  invisible(out_file)
}

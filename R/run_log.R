# ─────────────────────────────────────────────────────────────────────────────
# Run log
#
# Tracks per-method outcomes (success/failed, timing, feature count) across
# the correction-method loop in notame-workflow.r, since methods run inside
# a tryCatch and a failure is otherwise only visible by scrolling console
# output. Accumulated functionally: run_log <- log_method_result(run_log, ...).
# ─────────────────────────────────────────────────────────────────────────────

new_run_log <- function() list()

log_method_result <- function(run_log, method, status, duration_sec,
                               n_features = NA_integer_, error = NA_character_) {
  row <- data.frame(
    method        = method,
    status        = status,
    duration_sec  = round(duration_sec, 1),
    n_features    = n_features,
    error_message = error,
    stringsAsFactors = FALSE
  )
  c(run_log, list(row))
}

# Prints and saves the accumulated run log. Returns the combined data.frame
# (invisibly), or NULL if no methods were logged.
write_run_log <- function(run_log, interdir) {
  if (length(run_log) == 0) {
    message("No methods were run — skipping run_log.csv")
    return(invisible(NULL))
  }

  out <- do.call(rbind, run_log)
  cat("\n--- Run log ---\n")
  print(out, row.names = FALSE)
  write.csv(out, file.path(interdir, "run_log.csv"), row.names = FALSE)
  invisible(out)
}

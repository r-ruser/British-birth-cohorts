# Each output is one list-valued RDS containing every selected source table.
# Staged files are ASCII paths extracted read-only from the original ZIPs.
options(warn = 1)
suppressPackageStartupMessages(library(haven))

stage <- "E:/cohort_rds_stage_20260925"
output <- "E:/cohort_rds_output_20260925"
dir.create(output, recursive = TRUE, showWarnings = FALSE)
manifest <- as.data.frame(readr::read_csv(file.path(stage, "manifest.csv"),
                                           locale = readr::locale(encoding = "UTF-8"),
                                           show_col_types = FALSE,
                                           name_repair = "minimal"))
stopifnot(nrow(manifest) == 346L,
          all(file.exists(manifest$staged_file)),
          all(manifest$format %in% c(".dta", ".sav", ".csv", ".xls", ".xlsx",
                                     ".rds", ".rdata", ".sas7bdat")))

cohorts <- c("BCS70", "NCDS", "MCS", "NextSteps")
for (cohort in cohorts) {
  cat("START", cohort, format(Sys.time()), "\n")
  index <- which(manifest$cohort == cohort)
  inventory <- manifest[index, , drop = FALSE]
  inventory$rows <- NA_integer_
  inventory$columns <- NA_integer_
  datasets <- vector("list", nrow(inventory))
  names(datasets) <- inventory$key
  for (j in seq_len(nrow(inventory))) {
    path <- inventory$staged_file[j]
    format <- inventory$format[j]
    item <- tryCatch({
      if (format == ".dta") {
        haven::read_dta(path, .name_repair = "minimal")
      } else if (format == ".sav") {
        haven::read_sav(path, .name_repair = "minimal")
      } else if (format == ".csv") {
        read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
      } else if (format == ".rds") {
        readRDS(path)
      } else {
        # Ancillary spreadsheets and any other data format are embedded
        # byte-for-byte, avoiding alteration of sheets, labels or layout.
        readBin(path, what = "raw", n = file.info(path)$size)
      }
    }, error = function(e) {
      stop(sprintf("FAILED %s [%s]: %s", inventory$key[j],
                   inventory$source_entry[j], conditionMessage(e)))
    })
    datasets[[j]] <- item
    if (is.data.frame(item)) {
      inventory$rows[j] <- nrow(item)
      inventory$columns[j] <- ncol(item)
    }
    cat(cohort, j, "/", nrow(inventory), basename(path),
        if (is.data.frame(item)) paste0(nrow(item), "x", ncol(item)) else "raw",
        "\n")
    flush.console()
  }
  bundle <- list(
    cohort = cohort,
    format_version = 1L,
    structure_note = paste(
      "datasets is a named list, one entry per source data file;",
      "data tables retain original rows, columns, labels, and missing values;",
      "non-tabular data files are embedded as raw bytes."
    ),
    source_inventory = inventory,
    datasets = datasets
  )
  target <- file.path(output, paste0(cohort, "_full.rds"))
  temporary <- paste0(target, ".incomplete")
  saveRDS(bundle, temporary, compress = "gzip", version = 3)
  check <- readRDS(temporary)
  stopifnot(length(check$datasets) == nrow(inventory),
            identical(names(check$datasets), inventory$key),
            identical(check$source_inventory$source_bytes,
                      inventory$source_bytes))
  rm(check)
  gc()
  if (file.exists(target)) stop("Refusing to overwrite existing output: ", target)
  stopifnot(file.rename(temporary, target))
  cat("DONE", target, file.info(target)$size, "bytes\n")
  rm(bundle, datasets, inventory)
  gc()
}

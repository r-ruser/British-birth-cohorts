options(warn = 1, stringsAsFactors = FALSE)
Sys.setlocale("LC_ALL", "Chinese (Simplified)_China.utf8")
suppressPackageStartupMessages(library(haven))

root <- "E:/公共数据库"
zip_paths <- c(
  file.path(root, "BCS70", "2699stata_3DC3724A5E314FAA65AED88DED3D06626003FCFBB88EB269168A07EEC6D43C60_V1.zip"),
  file.path(root, "BCS70", "9347stata_81282F935DA9A8340DF8F9D160BD8423180AFF44B04D940C69DAE0BFEB56C916_V1.zip")
)
target_dir <- file.path(root, "Four_Cohorts_Full_RDS_20260925")
target <- file.path(target_dir, "BCS70_full.rds")
backup_dir <- file.path(target_dir, "_backups")
backup <- file.path(backup_dir, "BCS70_full_before_2699_9347_20261001.rds")
temporary <- paste0(target, ".incomplete_20261001")
audit_dir <- file.path(root, "P1", "bcs70_2699_9347_20261001")
stage <- "E:/cohort_rds_stage_20261001_bcs70_additions"
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(stage, recursive = TRUE, showWarnings = FALSE)
dir.create(backup_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(all(file.exists(zip_paths)), file.exists(target),
          !file.exists(backup), !file.exists(temporary))
old <- readRDS(target)
stopifnot(identical(old$cohort, "BCS70"),
          identical(names(old$datasets), old$source_inventory$key),
          length(old$datasets) == nrow(old$source_inventory))
old_n <- length(old$datasets)
old_keys <- names(old$datasets)
old_inventory <- old$source_inventory
old_digest <- digest::digest(old$datasets, algo = "xxhash64", serialize = TRUE)
cat("Original:", old_n, "datasets; content digest:", old_digest, "\n")

archive_entries <- list()
selected <- list()
for (zip in zip_paths) {
  listed <- unzip(zip, list = TRUE)
  listed <- listed[!grepl("/$", listed$Name), , drop = FALSE]
  stopifnot(nrow(listed) > 0L,
            !any(grepl("(^/|^[A-Za-z]:|\\.\\.)", listed$Name)),
            !anyDuplicated(listed$Name))
  ext <- paste0(".", tolower(tools::file_ext(listed$Name)))
  include <- ext %in% c(".dta", ".xls", ".xlsx")
  archive_entries[[basename(zip)]] <- data.frame(
    source_archive = normalizePath(zip, winslash = "/", mustWork = TRUE),
    source_entry = listed$Name,
    source_bytes = listed$Length,
    format = ext,
    included_in_rds = include,
    handling = ifelse(include, "table_or_spreadsheet", "documentation_or_source_code"),
    stringsAsFactors = FALSE
  )
  for (j in which(include)) {
    selected[[length(selected) + 1L]] <- list(
      zip = zip, entry = listed$Name[j], bytes = listed$Length[j],
      format = ext[j]
    )
  }
}
audit <- do.call(rbind, archive_entries)
write.csv(audit, file.path(audit_dir, "archive_entry_audit.csv"), row.names = FALSE)
stopifnot(length(selected) == 25L,
          sum(audit$included_in_rds) == 25L,
          !any(grepl("2699|9347", old_inventory$source_entry)))
cat("Archive entries:", nrow(audit), "; selected data files:", length(selected), "\n")

new_inventory <- old_inventory[0, , drop = FALSE]
new_data <- vector("list", length(selected))
new_keys <- sprintf("BCS70_%05d", old_n + seq_along(selected))
names(new_data) <- new_keys
stopifnot(!any(new_keys %in% old_keys))

for (i in seq_along(selected)) {
  item <- selected[[i]]
  stem <- sub("\\.zip$", "", basename(item$zip), ignore.case = TRUE)
  destination <- file.path(stage, stem)
  dir.create(destination, recursive = TRUE, showWarnings = FALSE)
  staged_file <- file.path(destination, item$entry)
  if (!file.exists(staged_file)) {
    extracted <- unzip(item$zip, files = item$entry, exdir = destination,
                       overwrite = FALSE)
    stopifnot(length(extracted) == 1L)
  }
  staged_file <- normalizePath(staged_file, winslash = "/", mustWork = TRUE)
  stopifnot(file.info(staged_file)$size == item$bytes)
  object <- if (item$format == ".dta") {
    haven::read_dta(staged_file, .name_repair = "minimal")
  } else {
    readBin(staged_file, what = "raw", n = item$bytes)
  }
  new_data[[i]] <- object
  info <- data.frame(
    cohort = "BCS70", key = new_keys[i],
    source_archive = normalizePath(item$zip, winslash = "/", mustWork = TRUE),
    source_entry = item$entry, staged_file = staged_file,
    format = item$format, source_bytes = item$bytes,
    rows = if (is.data.frame(object)) nrow(object) else NA_integer_,
    columns = if (is.data.frame(object)) ncol(object) else NA_integer_
  )
  new_inventory <- rbind(new_inventory, info[, names(new_inventory), drop = FALSE])
  cat(sprintf("%02d/%02d %s %s\n", i, length(selected), basename(item$entry),
              if (is.data.frame(object)) paste0(nrow(object), "x", ncol(object))
              else paste0(length(object), " raw bytes")))
  flush.console()
  rm(object)
}
stopifnot(length(new_data) == nrow(new_inventory),
          identical(names(new_data), new_inventory$key))
write.csv(new_inventory, file.path(audit_dir, "added_dataset_inventory.csv"),
          row.names = FALSE)

updated <- old
updated$source_inventory <- rbind(old_inventory, new_inventory)
updated$datasets <- c(old$datasets, new_data)
updated$augmentation_history <- c(old$augmentation_history, list(list(
  date = "2026-10-01", source_archives = basename(zip_paths),
  added_data_files = length(new_data),
  note = "Stata data tables retain all imported rows, columns, labels and missing values; Excel attachment is stored byte-for-byte. Other ZIP entries are documented in the external audit."
)))
stopifnot(length(updated$datasets) == old_n + 25L,
          identical(names(updated$datasets), updated$source_inventory$key),
          identical(updated$datasets[seq_len(old_n)], old$datasets))
cat("Saving candidate RDS...\n")
saveRDS(updated, temporary, compress = "gzip", version = 3)
stopifnot(file.exists(temporary), file.info(temporary)$size > file.info(target)$size)
rm(updated, old, new_data)
gc()

cat("Reading candidate RDS for verification...\n")
candidate <- readRDS(temporary)
stopifnot(identical(candidate$cohort, "BCS70"),
          length(candidate$datasets) == old_n + 25L,
          identical(names(candidate$datasets), candidate$source_inventory$key),
          identical(as.list(candidate$source_inventory[seq_len(old_n), , drop = FALSE]),
                    as.list(old_inventory)),
          identical(as.list(candidate$source_inventory[(old_n + 1L):nrow(candidate$source_inventory),
                                                       , drop = FALSE]), as.list(new_inventory)),
          identical(digest::digest(candidate$datasets[seq_len(old_n)],
                                   algo = "xxhash64", serialize = TRUE), old_digest))
for (i in seq_along(selected)) {
  object <- candidate$datasets[[old_n + i]]
  row <- new_inventory[i, ]
  if (row$format == ".dta") {
    stopifnot(is.data.frame(object), nrow(object) == row$rows,
              ncol(object) == row$columns)
  } else {
    stopifnot(is.raw(object), length(object) == row$source_bytes,
              identical(object, readBin(row$staged_file, "raw", n = row$source_bytes)))
  }
}
rm(candidate)
gc()

# Same-volume renames make the canonical replacement atomic. Restore the old
# file if the second rename fails; the original always remains in the backup.
stopifnot(normalizePath(dirname(target), winslash = "/") ==
            normalizePath(dirname(temporary), winslash = "/"))
stopifnot(file.rename(target, backup))
if (!file.rename(temporary, target)) {
  file.rename(backup, target)
  stop("Candidate RDS could not be moved into the canonical location")
}
stopifnot(file.exists(target), file.exists(backup))
writeLines(c(
  paste("Canonical:", target), paste("Backup:", backup),
  paste("Original datasets:", old_n),
  paste("Added tables:", sum(new_inventory$format == ".dta")),
  paste("Added spreadsheets:", sum(new_inventory$format %in% c(".xls", ".xlsx"))),
  paste("Final datasets:", old_n + nrow(new_inventory)),
  paste("Original content digest xxhash64:", old_digest),
  paste("ZIP entries audited:", nrow(audit)),
  paste("RDS bytes:", file.info(target)$size),
  "All 25 selected data files, original dataset content, table dimensions and Excel bytes verified after RDS reload."
), file.path(audit_dir, "verification_summary.txt"))
cat("DONE:", target, file.info(target)$size, "bytes\n")

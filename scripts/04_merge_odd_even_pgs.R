#!/usr/bin/env Rscript
#
# 04_merge_odd_even_pgs.R
#
# Combines the 22 per-chromosome PGS files (script 02) into a single
# odd-chromosome sum and a single even-chromosome sum per individual, for
# every trait in the manifest. This is the direct input to the GPD
# regression (script 05).
#
# Usage:
#   Rscript scripts/04_merge_odd_even_pgs.R
#
# Output (per trait):
#   <pgs_root>/<trait_id>/<trait_id>_pgs.txt   [columns: IID, pgs_odd, pgs_even]

suppressMessages(library(data.table))

project_root   <- "/path/to/project"
pgs_root       <- file.path(project_root, "data/pgs")
trait_manifest <- "config/trait_manifest.tsv"

traits <- fread(trait_manifest, header = TRUE)
traits$trait_id <- as.character(traits$trait_id)

for (i in seq_len(nrow(traits))) {

  trait_id <- traits$trait_id[i]
  workdir  <- file.path(pgs_root, trait_id)
  if (!dir.exists(workdir)) {
    message("Skipping ", trait_id, ": no PGS directory found")
    next
  }
  setwd(workdir)

  files <- list.files(path = workdir, pattern = "^pgs_.*\\.profile$", full.names = TRUE)
  if (length(files) == 0) {
    message("Skipping ", trait_id, ": no per-chromosome PGS files found")
    next
  }

  pgs_list <- lapply(files, function(f) {
    dt <- fread(f)
    chr <- gsub("[^0-9]", "", basename(f))
    dt <- dt[, .(IID, SCORE)]
    setnames(dt, "SCORE", paste0("SCORE_", chr))
    dt
  })

  merged <- Reduce(function(x, y) merge(x, y, by = "IID"), pgs_list)

  score_cols <- grep("^SCORE_", names(merged), value = TRUE)
  chr_nums   <- as.integer(gsub("^SCORE_", "", score_cols))
  odd_cols   <- score_cols[chr_nums %% 2 == 1]
  even_cols  <- score_cols[chr_nums %% 2 == 0]

  merged$pgs_odd  <- rowSums(merged[, ..odd_cols])
  merged$pgs_even <- rowSums(merged[, ..even_cols])

  out <- merged[, c("IID", "pgs_odd", "pgs_even")]
  fwrite(out, paste0(trait_id, "_pgs.txt"), quote = FALSE, sep = "\t")

  message("Merged odd/even PGS for ", trait_id, " (n = ", nrow(out), ")")
}

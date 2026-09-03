#!/usr/bin/env Rscript
#
# 02_ld_clump_and_score.R
#
# LD-clumps the harmonized summary statistics (script 01) and computes a
# per-chromosome polygenic score via PLINK's --score, for every trait in
# the manifest. Intended to be called once per chromosome — e.g. as a
# SLURM array job with the chromosome number passed as a command-line
# argument.
#
# Usage:
#   Rscript scripts/02_ld_clump_and_score.R <chromosome_number>
#
# Output (per trait, per chromosome):
#   <pgs_root>/<trait_id>/pgs_<chr>.profile

suppressMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: 02_ld_clump_and_score.R <chromosome_number>")
chr <- as.integer(args[1])

project_root   <- "/path/to/project"
sumstat_dir    <- file.path(project_root, "data/sumstats")
geno_dir       <- "/path/to/genotypes/per_chromosome"      # expects genotype_chr<N>.{bed,bim,fam}
pgs_root       <- file.path(project_root, "data/pgs")
trait_manifest <- "config/trait_manifest.tsv"

traits <- fread(trait_manifest, header = TRUE)
traits$trait_id <- as.character(traits$trait_id)

for (i in seq_len(nrow(traits))) {

  trait_id <- traits$trait_id[i]
  workdir  <- file.path(pgs_root, trait_id)
  dir.create(workdir, recursive = TRUE, showWarnings = FALSE)
  setwd(workdir)

  sumstat_path <- file.path(sumstat_dir, paste0(trait_id, "_GRCh37.gz"))
  if (!file.exists(sumstat_path)) {
    message("Skipping ", trait_id, ": no harmonized summary statistics found")
    next
  }

  sumstat <- fread(sumstat_path, header = TRUE)

  # --- LD clumping ---
  clump_input <- "clumping_sumstat.txt"
  if (!file.exists(clump_input)) {
    clump_df <- sumstat[, c("RSID", "P")]
    colnames(clump_df) <- c("SNP", "P")
    fwrite(clump_df, clump_input, quote = FALSE, sep = "\t")
  }

  geno_path <- file.path(geno_dir, paste0("genotype_chr", chr))
  clump_out_prefix <- paste0(trait_id, "_chr", chr)

  system(paste(
    "plink --bfile", geno_path,
    "--clump", clump_input,
    "--clump-p1 0.005 --clump-r2 0.1 --clump-kb 1000",
    "--out", clump_out_prefix
  ))

  clumped_file <- paste0(clump_out_prefix, ".clumped")
  if (!file.exists(clumped_file)) {
    message("No clumped SNPs for ", trait_id, " on chr", chr, ", skipping scoring")
    next
  }
  clumped_snps <- fread(clumped_file, header = TRUE)$SNP

  score_data <- sumstat[sumstat$RSID %in% clumped_snps, c("RSID", "EffectAllele", "B")]
  fwrite(score_data, paste0("score_file_", chr, ".txt"), quote = FALSE, sep = "\t")

  # --- PGS scoring ---
  system(paste(
    "plink --bfile", geno_path,
    "--score", paste0("score_file_", chr, ".txt"), "1 2 3 header",
    "--out", paste0("pgs_", chr)
  ))

  message("Scored chromosome ", chr, " for ", trait_id)
}

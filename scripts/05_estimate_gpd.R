#!/usr/bin/env Rscript
#
# 05_estimate_gpd.R
#
# Core within-trait GPD estimate: the odd/even-chromosome PGS regression
# (Yengo et al., 2018), adjusted for population structure via PCs computed
# separately per chromosome half (script 03). For each trait, two models
# are fit and the estimate is taken from whichever has the larger predictor
# variance, to improve numerical stability and reduce attenuation bias.
# Bonferroni correction is applied across all traits analyzed.
#
# Usage:
#   Rscript scripts/05_estimate_gpd.R
#
# Output:
#   <results_dir>/gpd_estimates.tsv

suppressMessages({
  library(data.table)
  library(dplyr)
})

project_root   <- "/path/to/project"
pgs_root       <- file.path(project_root, "data/pgs")
pc_dir         <- file.path(project_root, "work")
unrelated_list <- file.path(project_root, "config/unrelated_ancestry_matched_samples.txt")
trait_manifest <- "config/trait_manifest.tsv"
results_dir    <- "results"

dir.create(results_dir, showWarnings = FALSE)

traits <- fread(trait_manifest, header = TRUE)
traits$trait_id <- as.character(traits$trait_id)

unrel    <- fread(unrelated_list, header = FALSE)
pcs_odd  <- fread(file.path(pc_dir, "odd_chr.eigenvec"),  header = FALSE)
pcs_even <- fread(file.path(pc_dir, "even_chr.eigenvec"), header = FALSE)
pc_cols  <- paste0("V", 3:22)   # 20 PCs, columns 3-22 of GCTA .eigenvec output

out <- data.frame()

for (i in seq_len(nrow(traits))) {

  trait_id <- traits$trait_id[i]
  pgs_file <- file.path(pgs_root, trait_id, paste0(trait_id, "_pgs.txt"))
  if (!file.exists(pgs_file)) {
    message("Skipping ", trait_id, ": no merged PGS file found")
    next
  }

  pgs <- fread(pgs_file, header = TRUE) %>% filter(IID %in% unrel$V2)

  # condition each PGS on PCs computed from the OPPOSITE chromosome half
  df_odd  <- merge(pgs, pcs_even, by.x = "IID", by.y = "V2")   # models pgs_odd
  df_even <- merge(pgs, pcs_odd,  by.x = "IID", by.y = "V2")   # models pgs_even

  model_odd  <- lm(paste("pgs_odd ~ pgs_even +",  paste(pc_cols, collapse = " + ")), data = df_odd)
  model_even <- lm(paste("pgs_even ~ pgs_odd +",  paste(pc_cols, collapse = " + ")), data = df_even)

  coef_odd  <- summary(model_odd)$coefficients["pgs_even", , drop = FALSE]
  coef_even <- summary(model_even)$coefficients["pgs_odd", , drop = FALSE]

  v_odd  <- data.frame(model = "odd_on_even", trait_id = trait_id, as.data.frame(coef_odd))
  v_even <- data.frame(model = "even_on_odd", trait_id = trait_id, as.data.frame(coef_even))

  # retain the regression with the larger predictor variance (larger SE)
  v <- if (v_even$Std..Error > v_odd$Std..Error) v_even else v_odd
  out <- rbind(out, v)

  message("GPD estimate computed for ", trait_id)
}

names(out) <- c("model", "trait_id", "estimate", "se", "t_value", "p_value")

# Bonferroni correction across all traits analyzed
n_traits <- nrow(out)
out$p_bonferroni <- pmin(out$p_value * n_traits, 1)
out$significant  <- out$p_bonferroni < 0.05

out <- merge(out, traits[, c("trait_id", "label")], by = "trait_id")
out <- out[order(out$p_value), c("trait_id", "label", "model", "estimate", "se",
                                  "t_value", "p_value", "p_bonferroni", "significant")]

fwrite(out, file.path(results_dir, "gpd_estimates.tsv"), quote = FALSE, sep = "\t")
message("Done. Results written to ", file.path(results_dir, "gpd_estimates.tsv"))

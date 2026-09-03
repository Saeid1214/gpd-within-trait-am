# Genomic Signatures of Assortative Mating: Within-Trait GPD Estimation

Pipeline for estimating **directional Gametic Phase Disequilibrium (GPD)** —
a genomic signature of assortative mating (AM) — for individual traits,
using the odd/even-chromosome polygenic score (PGS) method
([Yengo et al., 2018](https://doi.org/10.1073/pnas.1815538115)).

This repository accompanies the manuscript *"Genomic Signatures of
Assortative Mating Across Psychiatric and Related Traits in the iPSYCH
Cohort."*

## Background

Assortative mating — the tendency of individuals to select partners with
similar traits — leaves a detectable signature in the genome: correlations
between trait-increasing alleles carried on different chromosomes, known as
gametic phase disequilibrium (GPD). Because odd- and even-numbered
chromosomes segregate independently at meiosis, any correlation between a
polygenic score built from odd chromosomes and one built from even
chromosomes reflects real assortment on that trait, rather than linkage.

This repository implements the estimation of that correlation — the
**within-trait GPD estimate** — for an arbitrary set of traits.

## Method summary

For each trait:

1. Polygenic scores are constructed **separately from odd- and
   even-numbered chromosomes**, using clumped, independent GWAS summary
   statistics with no sample overlap with the target cohort.
2. Two regressions are fit:
   - `PGS_odd  ~ PGS_even + 20 PCs(computed from even chromosomes)`
   - `PGS_even ~ PGS_odd  + 20 PCs(computed from odd chromosomes)`
3. The estimate is taken from whichever regression has the larger predictor
   variance (equivalently, larger standard error on the PGS coefficient),
   which improves numerical stability and reduces attenuation bias.
4. Significance is assessed after Bonferroni correction across all traits
   analyzed (p < 0.05 / N traits).

## Repository structure

```
.
├── scripts/
│   ├── 01_harmonize_sumstats.sh        # harmonize GWAS summary statistics across genome builds
│   ├── 02_ld_clump_and_score.R         # LD clumping + per-chromosome PGS scoring
│   ├── 03_odd_even_chromosome_pca.sh   # genotype split + PCs computed separately per chromosome half
│   ├── 04_merge_odd_even_pgs.R         # sum per-chromosome scores into odd/even PGS
│   └── 05_estimate_gpd.R               # core GPD regression + Bonferroni correction
├── config/
│   └── trait_manifest_example.tsv      # example trait manifest (fictional identifiers)
├── results/                            # pipeline output (not tracked; see .gitignore)
├── environment.yml
├── LICENSE
└── README.md
```

Scripts are numbered in execution order and are designed to be run once per
target cohort, looping internally over every trait listed in the manifest.

## Requirements

- [PLINK 1.9](https://www.cog-genomics.org/plink/)
- R ≥ 4.2, with `data.table`, `dplyr`
- See `environment.yml` for a ready-to-use conda environment

```bash
conda env create -f environment.yml
conda activate gpd-within-trait
```

## Usage

```bash
# 1. Harmonize summary statistics for every trait in the manifest
bash scripts/01_harmonize_sumstats.sh

# 2. Clump + score, once per chromosome (example: SLURM array over chr 1-22)
Rscript scripts/02_ld_clump_and_score.R "$SLURM_ARRAY_TASK_ID"

# 3. Genotype QC, odd/even chromosome split, and PCA
bash scripts/03_odd_even_chromosome_pca.sh

# 4. Merge per-chromosome scores into odd/even sums
Rscript scripts/04_merge_odd_even_pgs.R

# 5. Estimate GPD per trait
Rscript scripts/05_estimate_gpd.R
```

Each script reads its trait list from `config/trait_manifest.tsv` (schema
documented in `config/trait_manifest_example.tsv`) — add or remove traits
there rather than editing the scripts.

## Data availability

This repository documents **methodology**, not results reproduction. All
genotype, phenotype, and GWAS summary-statistic identifiers used in the
manuscript are placeholders here; the underlying iPSYCH data are held under
a restricted-access data use agreement and are not distributed. Access
requests: see the [iPSYCH data access policy](https://ipsych.dk/en/data-security/health-research-and-ethicalapproval/).

## Citation

If you use this pipeline, please cite the manuscript (citation to be added
upon publication) and the original method:

> Yengo L, et al. (2018). Imprint of assortative mating on the human genome. *PNAS*.

## License

Released under the MIT License — see [LICENSE](LICENSE).

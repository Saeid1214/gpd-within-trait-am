#!/usr/bin/env bash
#
# 03_odd_even_chromosome_pca.sh
#
# Prepares genotype data for the GPD regression: restricts to a QC'd,
# unrelated, ancestry-matched sample, then splits the genome into odd- and
# even-numbered chromosome sets and computes principal components SEPARATELY
# within each half.
#
# This odd/even PC pair — not whole-genome PCs — is what every GPD
# regression conditions on (script 05): the odd-chromosome PGS regression is
# adjusted for PCs computed from even chromosomes, and vice versa.
#
# Usage:
#   bash scripts/03_odd_even_chromosome_pca.sh
#
# Output:
#   ${WORK_DIR}/odd_chr.eigenvec
#   ${WORK_DIR}/even_chr.eigenvec

set -euo pipefail

PROJECT_ROOT="/path/to/project"
RAW_GENO_DIR="/path/to/raw_genotypes"           # per-chromosome source files
REF_SNPLIST="${PROJECT_ROOT}/ref/reference_panel.snplist"
UNRELATED_LIST="${PROJECT_ROOT}/config/unrelated_ancestry_matched_samples.txt"
GCTA_BIN="./gcta64"
THREADS=32
N_PCS=20

WORK_DIR="${PROJECT_ROOT}/work"
mkdir -p "${WORK_DIR}"
cd "${WORK_DIR}"

# --- restrict each chromosome to the reference SNP panel, then merge -------
for chr in $(seq 1 22); do
  plink --bfile "${RAW_GENO_DIR}/genotype_chr${chr}" \
        --threads "${THREADS}" \
        --extract "${REF_SNPLIST}" \
        --make-bed \
        --out "refpanel_chr${chr}"
done

ls refpanel_chr*.bed | sed 's/\.bed$//' > merge_list.txt
plink --merge-list merge_list.txt --threads "${THREADS}" \
      --make-bed --out refpanel_merged

# --- unrelated / ancestry-matched sample selection + MAF/HWE QC -----------
plink --bfile refpanel_merged \
      --keep "${UNRELATED_LIST}" \
      --maf 0.01 --hwe 1e-6 \
      --make-bed --threads "${THREADS}" \
      --out qc

# --- split into odd / even chromosome sets, PCA within each half ----------
odd_chrs=(1 3 5 7 9 11 13 15 17 19 21)
even_chrs=(2 4 6 8 10 12 14 16 18 20 22)

for half in odd even; do
  if [[ "${half}" == "odd" ]]; then chrs=("${odd_chrs[@]}"); else chrs=("${even_chrs[@]}"); fi

  half_files=()
  for chr in "${chrs[@]}"; do
    plink --bfile qc --chr "${chr}" --make-bed --out "qc_chr${chr}"
    half_files+=("qc_chr${chr}")
  done
  printf "%s\n" "${half_files[@]}" > "${half}_merge_list.txt"

  plink --merge-list "${half}_merge_list.txt" --threads "${THREADS}" \
        --make-bed --out "${half}_chr"

  "${GCTA_BIN}" --bfile "${half}_chr" --make-grm \
                --thread-num "${THREADS}" --out "${half}_chr"

  "${GCTA_BIN}" --grm "${half}_chr" --pca "${N_PCS}" \
                --thread-num "${THREADS}" --out "${half}_chr"
done

echo "Done. odd_chr.eigenvec and even_chr.eigenvec written to ${WORK_DIR}"

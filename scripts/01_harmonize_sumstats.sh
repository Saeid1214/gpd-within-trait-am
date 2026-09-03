#!/usr/bin/env bash
#
# scripts/01_harmonize_sumstats.sh
#
# Harmonizes externally sourced GWAS summary statistics (provided in two
# genome-build versions) into a single file per trait, ready for LD
# clumping and PGS construction (script 02).
#
# Usage:
#   bash scripts/01_harmonize_sumstats.sh
#
# Inputs (per trait, defined in config/trait_manifest.tsv):
#   ${SUMSTAT_SOURCE_DIR}/<sumstat_source_id>/cleaned_GRCh37.gz
#   ${SUMSTAT_SOURCE_DIR}/<sumstat_source_id>/cleaned_GRCh38.gz
#
# Output:
#   ${OUT_DIR}/<trait_id>_GRCh37.gz

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TRAIT_MANIFEST="${REPO_ROOT}/config/trait_manifest.tsv"

SUMSTAT_SOURCE_DIR="/path/to/external_sumstats"
OUT_DIR="/path/to/project/data/sumstats"

mkdir -p "${OUT_DIR}"

tail -n +2 "${TRAIT_MANIFEST}" | while IFS=$'\t' read -r trait_id label sumstat_id sample_size; do

  build37="${SUMSTAT_SOURCE_DIR}/${sumstat_id}/cleaned_GRCh37.gz"
  build38="${SUMSTAT_SOURCE_DIR}/${sumstat_id}/cleaned_GRCh38.gz"
  out_file="${OUT_DIR}/${trait_id}_GRCh37.gz"

  if [[ ! -f "${build37}" || ! -f "${build38}" ]]; then
    echo "WARNING: missing source files for '${trait_id}' (${label}), skipping" >&2
    continue
  fi

  echo "Harmonizing ${trait_id} (${label})"
  paste <(zcat "${build37}") <(zcat "${build38}") \
    | cut -f1-2,6- \
    | gzip -c > "${out_file}"

done

echo "Done. Harmonized summary statistics written to ${OUT_DIR}"

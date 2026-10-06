#!/usr/bin/env bash
# Usage: bash run_ldsc.sh LDSC_SCRIPT EFFECTIVE_N_SUMSTATS LD_PREFIX OUTPUT_DIR POP_PREV [POP_PREV ...]
set -euo pipefail

if (( $# < 5 )); then
  echo "Usage: bash run_ldsc.sh LDSC_SCRIPT EFFECTIVE_N_SUMSTATS LD_PREFIX OUTPUT_DIR POP_PREV [POP_PREV ...]" >&2
  exit 2
fi

ldsc_script=$1
sumstats=$2
ld_prefix=$3
output_dir=$4
shift 4

[[ -x "$ldsc_script" ]] || { echo "LDSC script missing or not executable" >&2; exit 2; }
[[ -f "$sumstats" ]] || { echo "Summary statistics missing" >&2; exit 2; }
[[ -f "${ld_prefix}1.l2.ldscore.gz" || -f "${ld_prefix}1.l2.ldscore" ]] || {
  echo "LD-score prefix did not resolve for chromosome 1" >&2
  exit 2
}
mkdir -p "$output_dir"

for prevalence in "$@"; do
  [[ "$prevalence" =~ ^0\.[0-9]+$ ]] || { echo "Invalid prevalence: $prevalence" >&2; exit 2; }
  prefix="${output_dir%/}/h2_effective_n_pop_${prevalence}"
  [[ ! -e "${prefix}.log" ]] || { echo "Refusing to overwrite ${prefix}.log" >&2; exit 2; }
  "$ldsc_script" --h2 "$sumstats" \
    --ref-ld-chr "$ld_prefix" --w-ld-chr "$ld_prefix" \
    --samp-prev 0.5 --pop-prev "$prevalence" --out "$prefix"
done

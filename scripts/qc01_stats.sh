#!/usr/bin/env bash
# qc01_stats.sh
#
# Purpose:
#   Compute basic FASTA statistics for the 6 haploid assemblies (3 cultivars
#   x 2 haplotypes) and produce a summary table with automatic pass/warn/fail
#   verdicts.
#
#   The verdicts are judged against EXPECTED VALUES FOR THE SPECIES, not against
#   universal thresholds - see docs section 4.5. The defaults below are for
#   Citrus; override them with environment variables for any other species.
#
# Metrics computed via seqkit stats:
#   - num_seqs   : contig/chromosome count
#   - sum_len    : total assembly size
#   - N50        : contiguity indicator
#   - GC%        : reflects species characteristics
#   - min/max    : shortest and longest sequence
#
# Usage:
#   bash qc01_stats.sh tables/samplesheet.tsv .
#
#   For a different species, pass its expected values, e.g. human (n=23,
#   ~3.1 Gb, GC ~41%):
#     EXPECTED_NUM_CHR=23 \
#     EXPECTED_SIZE_MIN=2900 EXPECTED_SIZE_MAX=3300 \
#     EXPECTED_GC_MIN=40.0 EXPECTED_GC_MAX=42.0 \
#     bash qc01_stats.sh tables/samplesheet.tsv .

set -euo pipefail

# Ensure seqkit is in PATH (project-local ./bin, or a system install)
SAMPLESHEET="${1:?samplesheet.tsv path required}"
PROJECT="${2:-.}"

# Project-local wrapper directory (see setup_singularity_wrappers.sh).
export PATH="${PROJECT}/bin:${PATH}"
OUT="${PROJECT}/results/qc"
mkdir -p "$OUT"

SUMMARY="${OUT}/stats_summary.tsv"

# ---------- Expected values (defaults: Citrus, haploid genome ~320 Mb) ----------
# All of these are species-specific. Override via the environment to run this
# script on any other organism.
EXPECTED_NUM_CHR="${EXPECTED_NUM_CHR:-9}"        # Haploid chromosome number
EXPECTED_SIZE_MIN="${EXPECTED_SIZE_MIN:-290}"    # Mb, minimum expected haploid size
EXPECTED_SIZE_MAX="${EXPECTED_SIZE_MAX:-370}"    # Mb, maximum expected haploid size
EXPECTED_GC_MIN="${EXPECTED_GC_MIN:-34.0}"       # %, minimum GC content
EXPECTED_GC_MAX="${EXPECTED_GC_MAX:-38.0}"       # %, maximum GC content

# Contiguity thresholds, derived from the expected chromosome number so that
# they scale with the organism. A fixed "<= 20 sequences" rule would flag even
# T2T-CHM13 (24 sequences) as semi-fragmented.
#   Citrus  (n=9)  -> CHR_SCALE_MAX=20,  FRAGMENTED_MIN=200
#   Human   (n=23) -> CHR_SCALE_MAX=48,  FRAGMENTED_MIN=480
CHR_SCALE_MAX="${CHR_SCALE_MAX:-$(( EXPECTED_NUM_CHR * 2 + 2 ))}"
FRAGMENTED_MIN="${FRAGMENTED_MIN:-$(( CHR_SCALE_MAX * 10 ))}"

echo ""
echo "========================================================"
echo " Citrus Pangenome QC Step 1: Assembly statistics"
echo "========================================================"
echo "  Expected: haploid $EXPECTED_SIZE_MIN-$EXPECTED_SIZE_MAX Mb"
echo "            chromosomes = $EXPECTED_NUM_CHR"
echo "            GC% = $EXPECTED_GC_MIN-$EXPECTED_GC_MAX"
echo "            chromosome-scale if num_seqs <= $CHR_SCALE_MAX,"
echo "            fragmented if num_seqs > $FRAGMENTED_MIN"
echo "========================================================"

# ---------- Initialize summary table header ----------
printf "sample\thap\tfile\tnum_seqs\tsum_len_Mb\tN50_Mb\tGC_pct\tmax_len_Mb\tmin_len_Mb\tsize_flag\tnum_flag\tGC_flag\tverdict\n" > "$SUMMARY"

# ---------- Iterate over each sample and haplotype ----------
awk 'NR>1' "$SAMPLESHEET" | while IFS=$'\t' read -r sid cvj sp prefix h1 h2 src ped; do
  for HAP in 1 2; do
    if [[ $HAP -eq 1 ]]; then FA="$h1"; else FA="$h2"; fi
    FA="$PROJECT/$FA"

    echo ""
    echo "-------- ${sid}_hap${HAP} : ${cvj} --------"

    if [[ ! -f "$FA" ]]; then
      echo "  [ERR] File not found: $FA"
      printf "%s\thap%s\t%s\tNA\tNA\tNA\tNA\tNA\tNA\tMISSING\tMISSING\tMISSING\tFAIL\n" \
        "$sid" "$HAP" "$FA" >> "$SUMMARY"
      continue
    fi

    # Get detailed stats from seqkit (tab-separated output with -T)
    stats=$(seqkit stats -a -T "$FA" | tail -n 1)

    # Extract each column from seqkit output.
    # Column layout: 1:file 2:format 3:type 4:num_seqs 5:sum_len 6:min_len 7:avg_len 8:max_len
    #                9:Q1 10:Q2 11:Q3 12:sum_gap 13:N50 14:N50_num 15:Q20 16:Q30 17:AvgQual 18:GC 19:sum_n
    num_seqs=$(echo "$stats"   | awk -F'\t' '{print $4}')
    sum_len=$(echo  "$stats"   | awk -F'\t' '{print $5}')
    min_len=$(echo  "$stats"   | awk -F'\t' '{print $6}')
    max_len=$(echo  "$stats"   | awk -F'\t' '{print $8}')
    n50=$(echo      "$stats"   | awk -F'\t' '{print $13}')
    gc_pct=$(echo   "$stats"   | awk -F'\t' '{print $18}')

    # Convert to Mb for display
    sum_mb=$(awk -v n="$sum_len" 'BEGIN{printf "%.2f", n/1000000}')
    n50_mb=$(awk -v n="$n50"     'BEGIN{printf "%.2f", n/1000000}')
    max_mb=$(awk -v n="$max_len" 'BEGIN{printf "%.2f", n/1000000}')
    min_mb=$(awk -v n="$min_len" 'BEGIN{printf "%.4f", n/1000000}')

    # Apply pass/warn/fail rules
    size_flag="OK"
    num_flag="OK"
    gc_flag="OK"
    verdict="PASS"

    # Size range check
    if awk -v s="$sum_mb" -v lo="$EXPECTED_SIZE_MIN" -v hi="$EXPECTED_SIZE_MAX" \
       'BEGIN{exit !(s<lo || s>hi)}'; then
      size_flag="OUT_OF_RANGE"
      verdict="WARN"
    fi

    # Chromosome-scale determination (thresholds scale with EXPECTED_NUM_CHR)
    if [[ "$num_seqs" -le "$CHR_SCALE_MAX" ]]; then
      num_flag="CHR_SCALE"
    elif [[ "$num_seqs" -le "$FRAGMENTED_MIN" ]]; then
      num_flag="SEMI_FRAGMENTED"
      [[ "$verdict" == "PASS" ]] && verdict="WARN"
    else
      num_flag="FRAGMENTED"
      verdict="FAIL"
    fi

    # GC content check
    if awk -v g="$gc_pct" -v lo="$EXPECTED_GC_MIN" -v hi="$EXPECTED_GC_MAX" \
       'BEGIN{exit !(g<lo || g>hi)}'; then
      gc_flag="OUT_OF_RANGE"
      [[ "$verdict" == "PASS" ]] && verdict="WARN"
    fi

    # Display results
    printf "  File           : %s\n" "$(basename "$FA")"
    printf "  num_seqs       : %s  [%s]\n" "$num_seqs" "$num_flag"
    printf "  Total length   : %s Mb  [%s]\n" "$sum_mb" "$size_flag"
    printf "  N50            : %s Mb\n" "$n50_mb"
    printf "  GC%%            : %s%%  [%s]\n" "$gc_pct" "$gc_flag"
    printf "  max/min contig : %s Mb / %s Mb\n" "$max_mb" "$min_mb"
    printf "  Verdict        : %s\n" "$verdict"

    # Append row to summary TSV
    printf "%s\thap%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
      "$sid" "$HAP" "$(basename "$FA")" \
      "$num_seqs" "$sum_mb" "$n50_mb" "$gc_pct" "$max_mb" "$min_mb" \
      "$size_flag" "$num_flag" "$gc_flag" "$verdict" >> "$SUMMARY"
  done
done

# ---------- Within-cultivar comparison (hap1 vs hap2) ----------
echo ""
echo "========================================================"
echo " Within-cultivar comparison (hap1 vs hap2)"
echo "========================================================"
echo "  Expectation: hap1 and hap2 of the same cultivar should be within ±15% in total length"

python3 - "$SUMMARY" <<'PY'
import sys, csv
from collections import defaultdict

path = sys.argv[1]
rows = []
with open(path) as f:
    reader = csv.DictReader(f, delimiter='\t')
    for r in reader:
        rows.append(r)

by_sample = defaultdict(dict)
for r in rows:
    by_sample[r['sample']][r['hap']] = r

print()
for sid, haps in by_sample.items():
    if 'hap1' in haps and 'hap2' in haps:
        try:
            s1 = float(haps['hap1']['sum_len_Mb'])
            s2 = float(haps['hap2']['sum_len_Mb'])
            diff_pct = abs(s1 - s2) / max(s1, s2) * 100
            if diff_pct <= 15:
                verdict = "OK"
            elif diff_pct <= 30:
                verdict = "WARN"
            else:
                verdict = "FAIL - possible phasing issue"
            print(f"  {sid}: hap1={s1:.1f}Mb, hap2={s2:.1f}Mb, diff={diff_pct:.1f}% ... {verdict}")
        except (ValueError, KeyError) as e:
            print(f"  {sid}: parse error")
PY

echo ""
echo "========================================================"
echo " Summary written to: $SUMMARY"
echo "========================================================"
echo ""
echo "Next steps:"
echo "  1. Investigate samples with verdict=FAIL"
echo "  2. Proceed with caution for samples with verdict=WARN"
echo "  3. If all PASS, run qc02_visualize.py to generate a summary plot"

#!/usr/bin/env bash
# download_plantgarden.sh
#
# Purpose:
#   Download the 6 haplotype-resolved assemblies of Satsuma mandarin (CUN),
#   Kishu mandarin (CKI), and Kunenbo (CKU) from Plant GARDEN
#   (https://plantgarden.jp), the Kazusa DNA Research Institute portal.
#
# Data source:
#   Isobe et al. (2023) bioRxiv 2023.06.02.543356
#   These are the publicly-released assemblies from the trio phasing study.
#
# Usage:
#   bash download_plantgarden.sh [output_directory]
#
#   Default output directory is ./data/raw
#
# Notes:
#   - Each file is approximately 90-105 MB (bgzip-compressed FASTA).
#   - Total download size is approximately 570 MB.
#   - Downloads are placed in per-cultivar subdirectories.

set -euo pipefail

OUT="${1:-./data/raw}"
mkdir -p "$OUT/CUN" "$OUT/CKI" "$OUT/CKU"

BASE="https://plantgarden.jp/ja/download"
UA="Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36"

# Helper function to download a file if not already present
download_if_absent() {
  local url="$1"
  local out_path="$2"
  local label="$3"
  
  if [[ -f "$out_path" ]]; then
    echo "  [SKIP] ${label}: already downloaded"
    return 0
  fi
  
  echo "  [DL  ] ${label}"
  curl -sSL -A "$UA" -o "$out_path" "$url"
  
  if [[ ! -s "$out_path" ]]; then
    echo "  [ERR ] ${label}: download failed or empty file"
    return 1
  fi
}

# ---------- Satsuma mandarin (Citrus unshiu) ----------
echo "[Satsuma mandarin (CUN)]"

# hap2 -> CUN#2 (Kishu-derived haplotype)
download_if_absent \
  "$BASE/Citrus_unshiu/t55188.G003/CUNphKi_r1.0.pmol.fasta.gz" \
  "$OUT/CUN/CUNphKi_r1.0.pmol.fasta.gz" \
  "CUN#2 (Kishu-derived): CUNphKi_r1.0"

# hap1 -> CUN#1 (Kunenbo-derived haplotype)
download_if_absent \
  "$BASE/Citrus_unshiu/t55188.G004/CUNphKu_r1.0.ch1-9.fasta.gz" \
  "$OUT/CUN/CUNphKu_r1.0.ch1-9.fasta.gz" \
  "CUN#1 (Kunenbo-derived): CUNphKu_r1.0"

# ---------- Kishu mandarin (Citrus kinokuni) ----------
echo ""
echo "[Kishu mandarin (CKI)]"

# NOTE: Update the URLs below by browsing:
#   https://plantgarden.jp/ja/download/Citrus_kinokuni/
# The subdirectory and filename format may change between releases.
# Placeholder URLs shown here - verify actual URLs before running.

download_if_absent \
  "$BASE/Citrus_kinokuni/t408488.G002/CKIhap1_r1.0.pmol.fasta.gz" \
  "$OUT/CKI/CKIhap1_r1.0.pmol.fasta.gz" \
  "hap1: CKIhap1_r1.0"
download_if_absent \
  "$BASE/Citrus_kinokuni/t408488.G003/CKIhap2_r1.0.pmol.fasta.gz" \
  "$OUT/CKI/CKIhap2_r1.0.pmol.fasta.gz" \
  "hap2: CKIhap2_r1.0"
#echo "  Please verify actual URLs at https://plantgarden.jp/ja/download/Citrus_kinokuni/"

# ---------- Kunenbo (Citrus nobilis) ----------
echo ""
echo "[Kunenbo (CKU)]"
download_if_absent \
  "$BASE/Citrus_nobilis/t481549.G002/CKUhap1_r1.0.pmol.fasta.gz" \
  "$OUT/CKU/CKUhap1_r1.0.pmol.fasta.gz" \
  "hap1: CKUhap1_r1.0"
download_if_absent \
  "$BASE/Citrus_nobilis/t481549.G003/CKUhap2_r1.0.pmol.fasta.gz" \
  "$OUT/CKU/CKUhap2_r1.0.pmol.fasta.gz" \
  "hap2: CKUhap2_r1.0"
#echo "  Please verify actual URLs at https://plantgarden.jp/ja/download/Citrus_nobilis/"

# ---------- Verification ----------
echo ""
echo "======================================================================"
echo " Downloaded files"
echo "======================================================================"
for cultivar in CUN CKI CKU; do
  echo ""
  echo "  ${cultivar}/"
  ls -lh "$OUT/${cultivar}/" 2>/dev/null | tail -n +2 | awk '{printf "    %s  %s\n", $5, $NF}'
done

echo ""
echo "Next step: verify the samplesheet paths in tables/samplesheet.tsv match"
echo "           the actual file locations shown above."

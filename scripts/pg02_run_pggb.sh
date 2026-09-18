#!/usr/bin/env bash
# pg02_run_pggb.sh
#
# Purpose:
#   Run PGGB on each chromosome to construct a pangenome graph.
#
# Features:
#   1. Preflight check for required tools (pggb, bgzip, samtools)
#   2. Automatic bgzip + faidx if the input FASTA is uncompressed
#   3. Input FASTA content validation (sequence count, PanSN naming)
#   4. Fault-tolerant execution (one failing chromosome does not abort the rest)
#
# PGGB parameters for Citrus (3 cultivars):
#   -p 95     : percent identity threshold
#   -s 10000  : segment length (bp)
#   -n 6      : number of haplotypes (3 cultivars x 2 hap)
#   -V CUN#1  : emit a VCF using CUN#1 as the coordinate reference.
#               The spec is REF[:LEN]; pggb appends the PanSN '#' to REF itself,
#               so 'CUN#1' yields `vg deconstruct -P CUN#1#`. A numeric LEN would
#               additionally emit a vcfbub/vcfwave-decomposed VCF; we do not use it.
#   -B 1G     : seqwish transclose batch size (see NOTE below)
#
#   -Y is left at its default. It is --exclude-delim (wfmash --group-prefix):
#   skip mappings whose query and target share the prefix before the LAST
#   occurrence of the character. With full PanSN names the group is
#   'sample#hap', so CUN#1 and CUN#2 are still aligned to each other.
#
# NOTE on -B:
#   In the PGGB Singularity image (revision 4225c6c), the default -B value
#   fails to parse and falls back silently to seqwish's own default of 1M,
#   causing under-compression. Always specify -B explicitly (e.g. -B 1G).
#   Verify with .params.yml -> transclose-batch: 1000000000 after run.
#
# Usage:
#   bash pg02_run_pggb.sh . [threads] [chr ...]
#
#   With no chromosomes named, every data/input/chr*.fa.gz is processed.
#   This script is a convenience wrapper for running all chromosomes; the
#   pggb command line it builds is spelled out in docs section 5.7.

# NOTE: -e is intentionally omitted so that a failure in one chromosome
# does not abort processing of the others.
set -uo pipefail

PROJECT="${1:-.}"
THREADS="${2:-32}"
shift 2 2>/dev/null || true

INPUT_DIR="${PROJECT}/data/input"
PGGB_OUT="${PROJECT}/results/pggb"
LOGS="${PGGB_OUT}/logs"
mkdir -p "$LOGS"

# Chromosomes to process: the ones named on the command line, or otherwise
# every prepared input in data/input. Nothing about "chr01..chr09" is assumed.
if [[ $# -gt 0 ]]; then
  CHRS=("$@")
else
  CHRS=()
  for f in "$INPUT_DIR"/*.fa.gz "$INPUT_DIR"/*.fa; do
    [[ -e "$f" ]] || continue
    b=$(basename "$f"); b=${b%.gz}; b=${b%.fa}
    CHRS+=("$b")
  done
  # de-duplicate (a chromosome may exist as both .fa and .fa.gz)
  if [[ ${#CHRS[@]} -gt 0 ]]; then
    mapfile -t CHRS < <(printf '%s\n' "${CHRS[@]}" | sort -u)
  fi
fi

if [[ ${#CHRS[@]} -eq 0 ]]; then
  echo "ERROR: no input FASTA found in $INPUT_DIR"
  echo "       Run pg01_prepare_input.sh first."
  exit 1
fi
echo "Chromosomes to process: ${CHRS[*]}"

# ---------- PGGB parameters (tuned for Citrus 3-cultivar dataset) ----------
PARAM_p=95
PARAM_s=10000
PARAM_n=6
PARAM_B=1G       # Must be explicitly set (see NOTE in header)
PARAM_V='CUN#1'

# ===================================================================
# Phase 0: Preflight checks
# ===================================================================
echo "======================================================================"
echo " PGGB Step 2: Per-chromosome pangenome graph construction"
echo "======================================================================"
echo ""
echo "-- Preflight check --"

PREFLIGHT_OK=1

if ! command -v pggb &>/dev/null; then
  echo "  MISSING: pggb"
  echo "    Install via:"
  echo "      conda:  mamba create -n pggb -c bioconda pggb && mamba activate pggb"
  echo "      docker: docker pull ghcr.io/pangenome/pggb:latest"
  echo "      singularity: singularity pull docker://ghcr.io/pangenome/pggb:latest"
  echo "                 + bash scripts/setup_singularity_wrappers.sh /path/to/pggb.sif"
  PREFLIGHT_OK=0
else
  PGGB_VERSION=$(pggb --version 2>&1 | head -1 || echo "unknown")
  echo "  OK: pggb ($PGGB_VERSION)"
fi

if ! command -v bgzip &>/dev/null; then
  echo "  MISSING: bgzip (required by PGGB for compressed FASTA input)"
  echo "    Install: mamba install -c bioconda htslib"
  PREFLIGHT_OK=0
else
  echo "  OK: bgzip ($(bgzip --version 2>&1 | head -1))"
fi

if ! command -v samtools &>/dev/null; then
  echo "  MISSING: samtools (required for .fai index creation)"
  PREFLIGHT_OK=0
else
  echo "  OK: samtools ($(samtools --version 2>&1 | head -1))"
fi

if [[ ! -d "$INPUT_DIR" ]]; then
  echo "  MISSING: input directory $INPUT_DIR"
  echo "    Run pg01_prepare_input.sh first."
  PREFLIGHT_OK=0
else
  N_INPUTS=$(ls "$INPUT_DIR"/*.fa "$INPUT_DIR"/*.fa.gz 2>/dev/null | wc -l)
  echo "  OK: input directory $INPUT_DIR ($N_INPUTS files detected)"
fi

if [[ $PREFLIGHT_OK -ne 1 ]]; then
  echo ""
  echo "ERROR: Prerequisites not met. Resolve the issues above and retry."
  exit 1
fi

echo ""
echo " Parameters:"
echo "   -p ${PARAM_p}     (percent identity)"
echo "   -s ${PARAM_s}  (segment length)"
echo "   -n ${PARAM_n}      (haplotype count)"
echo "   -B ${PARAM_B}      (seqwish transclose-batch; explicit to avoid parse issue)"
echo "   -V ${PARAM_V}  (VCF reference path)"
echo "   -Y '#'   (PanSN separator)"
echo "   threads: ${THREADS}"
echo ""

# ===================================================================
# Phase 1: Auto-bgzip + input validation
# ===================================================================
echo "-- Validate input FASTA and auto-bgzip if needed --"

for CHR in "${CHRS[@]}"; do
  IN_FA_GZ="$INPUT_DIR/${CHR}.fa.gz"
  IN_FA="$INPUT_DIR/${CHR}.fa"
  
  if [[ -f "$IN_FA_GZ" && -f "${IN_FA_GZ}.fai" ]]; then
    echo "  OK: $CHR (bgzip+index ready)"
  elif [[ -f "$IN_FA_GZ" ]]; then
    echo "  [INDEX] $CHR: creating .fai index"
    samtools faidx "$IN_FA_GZ" 2>&1 | sed 's/^/    /'
  elif [[ -f "$IN_FA" ]]; then
    echo "  [BGZIP] $CHR: compressing .fa and creating index"
    bgzip -f "$IN_FA"
    samtools faidx "$IN_FA_GZ" 2>&1 | sed 's/^/    /'
  else
    echo "  WARN: $CHR input file not found, skipping"
    continue
  fi
  
  # Content validation
  if [[ -f "$IN_FA_GZ" ]]; then
    N_SEQ=$(zcat "$IN_FA_GZ" | grep -c '^>')
    if [[ "$N_SEQ" -ne 6 ]]; then
      echo "  WARN: $CHR has $N_SEQ sequences but PGGB -n is set to 6"
    fi
    FIRST_HEADER=$(zcat "$IN_FA_GZ" | grep '^>' | head -1 | sed 's/^>//')
    if [[ ! "$FIRST_HEADER" =~ ^[a-zA-Z0-9_]+#[12]#chr[0-9]+ ]]; then
      echo "  WARN: $CHR first header may not be in PanSN format: $FIRST_HEADER"
    fi
  fi
done

# ===================================================================
# Phase 2: PGGB execution (per chromosome, fault-tolerant)
# ===================================================================
echo ""
echo "======================================================================"
echo " PGGB execution"
echo "======================================================================"

SUCCESS=0
FAILED=0
SKIPPED=0

for CHR in "${CHRS[@]}"; do
  IN_FA_GZ="$INPUT_DIR/${CHR}.fa.gz"
  OUT_DIR="$PGGB_OUT/${CHR}"
  
  if [[ ! -f "$IN_FA_GZ" ]]; then
    echo "[SKIP] $CHR: input file missing"
    SKIPPED=$((SKIPPED+1))
    continue
  fi
  
  # PGGB output filenames embed parameter hashes: chr09.fa.gz.<h1>.<h2>.<h3>.smooth.final.og
  EXISTING_OG=$(ls "$OUT_DIR/${CHR}.fa.gz."*.smooth.final.og 2>/dev/null | head -1)
  if [[ -n "$EXISTING_OG" ]]; then
    echo "[SKIP] $CHR: already complete ($OUT_DIR)"
    SKIPPED=$((SKIPPED+1))
    continue
  fi
  
  echo ""
  echo "======================================================================"
  echo " [$(date +%Y-%m-%d\ %H:%M:%S)] Starting $CHR"
  echo "======================================================================"
  
  mkdir -p "$OUT_DIR"
  
  # Trailing "|| true" so a chromosome failure does not abort the whole loop.
  pggb \
    -i "$IN_FA_GZ" \
    -o "$OUT_DIR" \
    -p "$PARAM_p" \
    -s "$PARAM_s" \
    -n "$PARAM_n" \
    -B "$PARAM_B" \
    -V "$PARAM_V" \
    -t "$THREADS" \
    -m \
    -S \
    > "$LOGS/pggb_${CHR}.log" 2>&1 || true
  
  FINAL_OG=$(ls "$OUT_DIR/${CHR}.fa.gz."*.smooth.final.og 2>/dev/null | head -1)
  if [[ -n "$FINAL_OG" ]]; then
    echo "[$(date +%H:%M:%S)] OK: $CHR complete ($(basename "$FINAL_OG"))"
    SUCCESS=$((SUCCESS+1))
  else
    echo "[$(date +%H:%M:%S)] FAIL: $CHR (see logs/pggb_${CHR}.log)"
    echo "    Last lines of log:"
    tail -5 "$LOGS/pggb_${CHR}.log" 2>/dev/null | sed 's/^/      /'
    FAILED=$((FAILED+1))
  fi
done

# ===================================================================
# Phase 3: Summary
# ===================================================================
echo ""
echo "======================================================================"
echo " Execution summary"
echo "======================================================================"
echo "  succeeded: $SUCCESS"
echo "  failed   : $FAILED"
echo "  skipped  : $SKIPPED"
echo ""
echo "  Per-chromosome status:"

for CHR in "${CHRS[@]}"; do
  OUT_DIR="$PGGB_OUT/${CHR}"
  FINAL_OG=$(ls "$OUT_DIR/${CHR}.fa.gz."*.smooth.final.og 2>/dev/null | head -1)
  FINAL_GFA=$(ls "$OUT_DIR/${CHR}.fa.gz."*.smooth.final.gfa 2>/dev/null | head -1)
  
  if [[ -n "$FINAL_OG" ]]; then
    OG_SIZE=$(du -h "$FINAL_OG" 2>/dev/null | cut -f1)
    GFA_SIZE=$(du -h "$FINAL_GFA" 2>/dev/null | cut -f1 || echo "N/A")
    echo "    OK   $CHR: OG=$OG_SIZE, GFA=$GFA_SIZE"
  elif [[ -d "$OUT_DIR" ]]; then
    echo "    FAIL $CHR: incomplete (intermediate files present)"
  else
    echo "    -    $CHR: not run"
  fi
done

echo ""
if [[ $FAILED -gt 0 ]]; then
  echo "WARN: ${FAILED} chromosome(s) failed. Check:"
  echo "     1. logs/pggb_chrXX.log for detailed errors"
  echo "     2. Memory pressure: reduce threads or use a larger node"
  echo "     3. Parameter tuning: adjust -p / -s and rerun only failed chromosomes"
  exit 1
elif [[ $SUCCESS -eq 0 && $SKIPPED -eq 0 ]]; then
  echo "WARN: no chromosomes were processed. Check inputs and parameters."
  exit 1
elif [[ $SUCCESS -eq 0 && $SKIPPED -gt 0 ]]; then
  echo "NOTE: no new work performed (all chromosomes already complete)."
  echo "      Next: bash scripts/pg03_qc_graph.sh"
else
  echo "OK: ${SUCCESS} chromosome(s) newly complete, ${SKIPPED} already done."
  echo "    Next: bash scripts/pg03_qc_graph.sh"
fi

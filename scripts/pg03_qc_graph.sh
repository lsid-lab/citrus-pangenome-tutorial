#!/usr/bin/env bash
# pg03_qc_graph.sh
#
# Purpose:
#   Evaluate the quality of a PGGB-constructed pangenome graph across four
#   levels:
#     Level 1: Structural check (odgi stats + .params.yml verification)
#     Level 2: Path preservation (input path lengths vs graph path lengths)
#     Level 3: Pedigree consistency (expected relationship: CUN = F1(CKI mother, CKU father);
#              CUN#1 = CUNphKu is the paternal hap, CUN#2 = CUNphKi the maternal one)
#     Level 4: SV signal (vg deconstruct output)
#
# Verdict rules:
#   PASS : all levels OK
#   FAIL : Level 1 or 2 has fatal issue (graph is not usable)
#   WARN : only Level 3 or 4 has issues (graph is usable, interpret carefully)
#
# Usage:
#   bash pg03_qc_graph.sh . [chromosome (default: all)]

# NOTE: -e is intentionally omitted so that a failure in one command
# does not silently abort the whole script.
set -uo pipefail

PROJECT="${1:-.}"
CHR_FILTER="${2:-all}"

INPUT_DIR="${PROJECT}/data/input"
PGGB_OUT="${PROJECT}/results/pggb"
QC_OUT="${PROJECT}/results/graph_qc"
mkdir -p "$QC_OUT"

echo "======================================================================"
echo " PGGB Step 3: Pangenome graph QC"
echo "======================================================================"

# ===================================================================
# Preflight check
# ===================================================================
echo ""
echo "-- Preflight check --"

PREFLIGHT_OK=1
MISSING_TOOLS=()

for TOOL in odgi vg samtools bcftools seqkit; do
  if command -v "$TOOL" &>/dev/null; then
    echo "  OK: $TOOL ($(command -v $TOOL))"
  else
    echo "  MISSING: $TOOL"
    MISSING_TOOLS+=("$TOOL")
    PREFLIGHT_OK=0
  fi
done

if [[ $PREFLIGHT_OK -ne 1 ]]; then
  echo ""
  echo "ERROR: missing tools on host: ${MISSING_TOOLS[*]}"
  echo ""
  echo " These tools are packaged inside the PGGB Singularity image."
  echo " Create wrappers with the setup script:"
  echo ""
  echo "   1. Locate the SIF file:"
  echo "      readlink -f pggb_latest.sif"
  echo ""
  echo "   2. Create wrappers:"
  echo "      bash scripts/setup_singularity_wrappers.sh /full/path/to/pggb_latest.sif"
  echo ""
  echo "   3. Ensure PATH includes ~/bin (also in your SLURM script):"
  echo "      export PATH=\"\$HOME/bin:\$PATH\""
  echo ""
  echo "   4. Rerun pg03"
  exit 1
fi

echo ""

# ===================================================================
# Initialize summary table
# ===================================================================
SUMMARY_TSV="${QC_OUT}/graph_qc_summary.tsv"
printf "chr\tinput_bp\tnode_count\tedge_count\tpath_count\tgraph_bp\tinflation\tpath_length_ok\tSV_count\tverdict\n" > "$SUMMARY_TSV"

CHROMS_TO_PROCESS=()
for CHR_DIR in "$PGGB_OUT"/*/; do
  [[ ! -d "$CHR_DIR" ]] && continue
  CHR=$(basename "$CHR_DIR" | sed 's/_pggb//')
  [[ "$CHR_FILTER" != "all" && "$CHR_FILTER" != "$CHR" ]] && continue
  CHROMS_TO_PROCESS+=("$CHR")
done

if [[ ${#CHROMS_TO_PROCESS[@]} -eq 0 ]]; then
  echo "WARN: no target chromosomes found."
  echo "  Check that pg02 has completed, or that CHR_FILTER (${CHR_FILTER}) is correct."
  exit 1
fi

echo "-- Chromosomes to process: ${CHROMS_TO_PROCESS[*]} --"
echo ""

# ===================================================================
# QC loop per chromosome
# ===================================================================
for CHR in "${CHROMS_TO_PROCESS[@]}"; do
  CHR_DIR="$PGGB_OUT/${CHR}"
  # PGGB output filenames embed parameter hashes
  OG=$(ls "$CHR_DIR/${CHR}.fa.gz."*.smooth.final.og 2>/dev/null | head -1)
  GFA=$(ls "$CHR_DIR/${CHR}.fa.gz."*.smooth.final.gfa 2>/dev/null | head -1)
  IN_FA="$INPUT_DIR/${CHR}.fa.gz"

  if [[ -z "$OG" || ! -f "$OG" ]]; then
    echo "[SKIP] $CHR: OG file not found (pg02 incomplete or naming mismatch)"
    continue
  fi

  echo "======================================================================"
  echo "  $CHR"
  echo "  OG: $(basename "$OG")"
  echo "======================================================================"

  # Initialize per-chr variables
  NODES="NA"; EDGES="NA"; PATHS="NA"; GRAPH_BP="NA"; STEPS="NA"
  INPUT_BP="NA"; INFLATION="NA"; COMPRESSION="NA"
  PATH_OK="NA"
  N_SV="NA"; N_SNP="NA"; N_TOTAL="NA"
  L1_OK="?"; L2_OK="?"; L3_OK="?"; L4_OK="?"

  # ============================================================
  # Level 1: Structural check + version-specific -B verification
  # ============================================================
  echo ""
  echo "  --- Level 1: Structure (odgi stats) + transclose-batch check ---"

  # Verify .params.yml transclose-batch (fastest -B parse verification)
  PARAMS_YML=$(ls "$CHR_DIR/"*.params.yml 2>/dev/null | head -1)
  if [[ -n "$PARAMS_YML" ]]; then
    TB=$(grep -oP 'transclose-batch:\s*\K\d+' "$PARAMS_YML" 2>/dev/null || echo "?")
    echo "    transclose-batch: $TB"
    if [[ "$TB" == "1000000" ]]; then
      echo "    WARN: -B default parse issue detected, -B 1G was NOT applied"
      echo "          Rerun pg02 with -B 1G explicit for correct compression"
    elif [[ "$TB" -ge 1000000000 ]]; then
      echo "    OK: -B 1G (or larger) correctly applied"
    fi
  fi

  if STATS=$(odgi stats -i "$OG" -S 2>&1); then
    STATS_LINE=$(echo "$STATS" | tail -n 1)
    # odgi stats -S output columns: length nodes edges paths steps
    # col1=length(graph_bp), col2=nodes, col3=edges, col4=paths, col5=steps
    GRAPH_BP=$(echo "$STATS_LINE" | cut -f1)
    NODES=$(echo "$STATS_LINE" | cut -f2)
    EDGES=$(echo "$STATS_LINE" | cut -f3)
    PATHS=$(echo "$STATS_LINE" | cut -f4)
    STEPS=$(echo "$STATS_LINE" | cut -f5)
    
    INPUT_BP=$(zcat "$IN_FA" 2>/dev/null | seqkit stats -T 2>/dev/null | tail -1 | cut -f5)
    if [[ -n "$INPUT_BP" && "$GRAPH_BP" =~ ^[0-9]+$ ]]; then
      # Compression ratio = graph_bp / input_bp (target <= 0.35)
      COMPRESSION=$(awk -v g="$GRAPH_BP" -v i="$INPUT_BP" 'BEGIN{if(i>0) printf "%.3f", g/i; else print "NA"}')
      # inflation = compression ratio x n = graph_bp / (input_bp / n)
      INFLATION=$(awk -v g="$GRAPH_BP" -v i="$INPUT_BP" 'BEGIN{if(i>0) printf "%.3f", g/i*6; else print "NA"}')
    fi
    
    echo "    nodes    : $NODES"
    echo "    edges    : $EDGES"
    echo "    paths    : $PATHS  (expected: 6)"
    echo "    steps    : $STEPS"
    if [[ "$GRAPH_BP" =~ ^[0-9]+$ ]]; then
      echo "    graph_bp : $(awk -v n=$GRAPH_BP 'BEGIN{printf "%.1f Mb", n/1e6}')"
    fi
    if [[ "$INPUT_BP" =~ ^[0-9]+$ ]]; then
      echo "    input_bp : $(awk -v n=$INPUT_BP 'BEGIN{printf "%.1f Mb (6 sequences total)", n/1e6}')"
    fi
    echo "    ** compression : $COMPRESSION  (target <= 0.35, colleague achieved 0.19-0.28)"
    echo "       inflation   : $INFLATION  (informational, expected 1.0-1.6)"
    
    L1_OK="OK"
    [[ "$PATHS" != "6" ]] && L1_OK="FAIL: paths $PATHS != 6"
    if [[ "$COMPRESSION" != "NA" ]]; then
      if awk -v c="$COMPRESSION" 'BEGIN{exit !(c>0.35)}'; then
        L1_OK="FAIL: compression $COMPRESSION > 0.35 (verify -B parameter)"
      fi
    fi
    echo "    Verdict  : $L1_OK"
  else
    echo "    FAIL: odgi stats failed"
    echo "    stderr: $STATS" | head -3
    L1_OK="FAIL: odgi stats error"
  fi

  # ============================================================
  # Level 2: Path preservation (via FASTA extraction)
  # ============================================================
  echo ""
  echo "  --- Level 2: Path preservation (odgi paths -f -> seqkit) ---"

  # Extract paths as FASTA (most reliable way to get path lengths across
  # odgi versions, since -L only lists names and -l flags vary).
  PATHS_FA="$QC_OUT/${CHR}_paths.fa"
  if odgi paths -i "$OG" -f > "$PATHS_FA" 2>/dev/null && [[ -s "$PATHS_FA" ]]; then
    seqkit fx2tab -nl "$PATHS_FA" > "$QC_OUT/${CHR}_paths.tsv" 2>/dev/null
    
    echo "    Path lengths (bp):"
    awk -F'\t' '{printf "      %-30s %s bp\n", $1, $2}' "$QC_OUT/${CHR}_paths.tsv" | head -20

    # Get expected lengths from input FASTA
    seqkit fx2tab -nl "$IN_FA" > "$QC_OUT/${CHR}_input_lengths.tsv" 2>/dev/null

    PATH_OK=$(python3 <<PYEOF
paths = {}
try:
    with open("$QC_OUT/${CHR}_paths.tsv") as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                try:
                    paths[parts[0]] = int(parts[1])
                except ValueError:
                    pass
    inputs = {}
    with open("$QC_OUT/${CHR}_input_lengths.tsv") as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                try:
                    inputs[parts[0]] = int(parts[1])
                except ValueError:
                    pass
    if not paths or not inputs:
        print("NO_DATA_paths={} inputs={}".format(len(paths), len(inputs)))
    else:
        max_diff_pct = 0
        max_diff_name = ""
        for name in inputs:
            if name in paths:
                diff_pct = abs(paths[name] - inputs[name]) / inputs[name] * 100
                if diff_pct > max_diff_pct:
                    max_diff_pct = diff_pct
                    max_diff_name = name
            else:
                print(f"MISSING_PATH:{name}")
                exit()
        if max_diff_pct <= 1.0:
            print(f"YES_max_diff={max_diff_pct:.3f}%_at_{max_diff_name}")
        else:
            print(f"NO_max_diff={max_diff_pct:.2f}%_at_{max_diff_name}")
except Exception as e:
    print(f"ERROR:{e}")
PYEOF
)
    echo "    Preservation: $PATH_OK"
    L2_OK="OK"
    [[ "$PATH_OK" == NO* ]] && L2_OK="WARN: $PATH_OK"
    [[ "$PATH_OK" == MISSING* ]] && L2_OK="FAIL: $PATH_OK"
    [[ "$PATH_OK" == ERROR* ]] && L2_OK="FAIL: $PATH_OK"
    [[ "$PATH_OK" == NO_DATA* ]] && L2_OK="FAIL: $PATH_OK"
  else
    echo "    FAIL: odgi paths -f failed or produced empty output"
    L2_OK="FAIL: odgi paths error"
  fi

  # ============================================================
  # Level 3: Pedigree consistency
  # ============================================================
  echo ""
  echo "  --- Level 3: Pedigree consistency ---"

  if odgi similarity -i "$OG" > "$QC_OUT/${CHR}_similarity.tsv" 2>/dev/null && \
     [[ -s "$QC_OUT/${CHR}_similarity.tsv" ]]; then
    
    HEADER=$(head -1 "$QC_OUT/${CHR}_similarity.tsv")
    echo "    [debug] similarity columns: $HEADER"
    
    echo ""
    echo "    Expected pedigree (CUN = CKI x CKU F1):"
    echo "      CUN#1 (CUNphKu, paternal) should score higher against CKU than against CKI"
    echo "      CUN#2 (CUNphKi, maternal) should score higher against CKI than against CKU"
    echo "      (each is a recombinant mosaic of that parent's two haplotypes, so it need"
    echo "       not match a single parental haplotype - see docs section 2.4)"

    JACCARD_COL=$(echo "$HEADER" | tr '\t' '\n' | grep -n -i "jaccard" | head -1 | cut -d: -f1)
    [[ -z "$JACCARD_COL" ]] && JACCARD_COL=6

    # Show individual pair values (CUN#1 and CUN#2 vs all parent haps)
    echo ""
    echo "    --- Individual pair similarities (Jaccard) ---"
    echo ""
    echo "    CUN#1 (CUNphKu, Kunenbo-derived) vs each parent haplotype:"
    awk -F'\t' -v c=$JACCARD_COL 'NR>1 && $1 ~ /^CUN#1/ && ($2 ~ /^CKI/ || $2 ~ /^CKU/) {
      printf "      %-30s  Jaccard = %.4f\n", $2, $c
    } NR>1 && $2 ~ /^CUN#1/ && ($1 ~ /^CKI/ || $1 ~ /^CKU/) {
      printf "      %-30s  Jaccard = %.4f\n", $1, $c
    }' "$QC_OUT/${CHR}_similarity.tsv" | sort -u

    echo ""
    echo "    CUN#2 (CUNphKi, Kishu-derived) vs each parent haplotype:"
    awk -F'\t' -v c=$JACCARD_COL 'NR>1 && $1 ~ /^CUN#2/ && ($2 ~ /^CKI/ || $2 ~ /^CKU/) {
      printf "      %-30s  Jaccard = %.4f\n", $2, $c
    } NR>1 && $2 ~ /^CUN#2/ && ($1 ~ /^CKI/ || $1 ~ /^CKU/) {
      printf "      %-30s  Jaccard = %.4f\n", $1, $c
    }' "$QC_OUT/${CHR}_similarity.tsv" | sort -u

    # AVG-based aggregate
    CUN1_CKI_AVG=$(awk -F'\t' -v c=$JACCARD_COL 'NR>1 && (($1 ~ /^CUN#1/ && $2 ~ /^CKI/) || ($2 ~ /^CUN#1/ && $1 ~ /^CKI/)) {sum+=$c; n++} END{if(n>0) printf "%.3f", sum/n; else print "NA"}' "$QC_OUT/${CHR}_similarity.tsv")
    CUN1_CKU_AVG=$(awk -F'\t' -v c=$JACCARD_COL 'NR>1 && (($1 ~ /^CUN#1/ && $2 ~ /^CKU/) || ($2 ~ /^CUN#1/ && $1 ~ /^CKU/)) {sum+=$c; n++} END{if(n>0) printf "%.3f", sum/n; else print "NA"}' "$QC_OUT/${CHR}_similarity.tsv")
    CUN2_CKI_AVG=$(awk -F'\t' -v c=$JACCARD_COL 'NR>1 && (($1 ~ /^CUN#2/ && $2 ~ /^CKI/) || ($2 ~ /^CUN#2/ && $1 ~ /^CKI/)) {sum+=$c; n++} END{if(n>0) printf "%.3f", sum/n; else print "NA"}' "$QC_OUT/${CHR}_similarity.tsv")
    CUN2_CKU_AVG=$(awk -F'\t' -v c=$JACCARD_COL 'NR>1 && (($1 ~ /^CUN#2/ && $2 ~ /^CKU/) || ($2 ~ /^CUN#2/ && $1 ~ /^CKU/)) {sum+=$c; n++} END{if(n>0) printf "%.3f", sum/n; else print "NA"}' "$QC_OUT/${CHR}_similarity.tsv")

    # MAX-based aggregate (the biologically meaningful pedigree indicator)
    CUN1_CKI_MAX=$(awk -F'\t' -v c=$JACCARD_COL 'NR>1 && (($1 ~ /^CUN#1/ && $2 ~ /^CKI/) || ($2 ~ /^CUN#1/ && $1 ~ /^CKI/)) {if($c>max) max=$c} END{printf "%.3f", max}' "$QC_OUT/${CHR}_similarity.tsv")
    CUN1_CKU_MAX=$(awk -F'\t' -v c=$JACCARD_COL 'NR>1 && (($1 ~ /^CUN#1/ && $2 ~ /^CKU/) || ($2 ~ /^CUN#1/ && $1 ~ /^CKU/)) {if($c>max) max=$c} END{printf "%.3f", max}' "$QC_OUT/${CHR}_similarity.tsv")
    CUN2_CKI_MAX=$(awk -F'\t' -v c=$JACCARD_COL 'NR>1 && (($1 ~ /^CUN#2/ && $2 ~ /^CKI/) || ($2 ~ /^CUN#2/ && $1 ~ /^CKI/)) {if($c>max) max=$c} END{printf "%.3f", max}' "$QC_OUT/${CHR}_similarity.tsv")
    CUN2_CKU_MAX=$(awk -F'\t' -v c=$JACCARD_COL 'NR>1 && (($1 ~ /^CUN#2/ && $2 ~ /^CKU/) || ($2 ~ /^CUN#2/ && $1 ~ /^CKU/)) {if($c>max) max=$c} END{printf "%.3f", max}' "$QC_OUT/${CHR}_similarity.tsv")

    echo ""
    echo "    --- Aggregate values ---"
    printf "    %-15s  %-15s  %-15s\n" " " "vs Kishu" "vs Kunenbo"
    printf "    %-15s  AVG=%.3f MAX=%.3f  AVG=%.3f MAX=%.3f\n" "CUN#1 (CKU-der)" $CUN1_CKI_AVG $CUN1_CKI_MAX $CUN1_CKU_AVG $CUN1_CKU_MAX
    printf "    %-15s  AVG=%.3f MAX=%.3f  AVG=%.3f MAX=%.3f\n" "CUN#2 (CKI-der)" $CUN2_CKI_AVG $CUN2_CKI_MAX $CUN2_CKU_AVG $CUN2_CKU_MAX

    # AVG-based and MAX-based pedigree tests
    AVG_TEST=$(python3 <<PYEOF
try:
    a=float("$CUN1_CKI_AVG"); b=float("$CUN1_CKU_AVG"); c=float("$CUN2_CKI_AVG"); d=float("$CUN2_CKU_AVG")
    # CUN#1 = CUNphKu (Kunenbo-derived) -> expect CKU > CKI
    # CUN#2 = CUNphKi (Kishu-derived)   -> expect CKI > CKU
    print("YES" if ((b > a) and (c > d)) else "NO")
except: print("NA")
PYEOF
)
    MAX_TEST=$(python3 <<PYEOF
try:
    a=float("$CUN1_CKI_MAX"); b=float("$CUN1_CKU_MAX"); c=float("$CUN2_CKI_MAX"); d=float("$CUN2_CKU_MAX")
    # CUN#1 = CUNphKu (Kunenbo-derived) -> expect CKU > CKI
    # CUN#2 = CUNphKi (Kishu-derived)   -> expect CKI > CKU
    print("YES" if ((b > a) and (c > d)) else "NO")
except: print("NA")
PYEOF
)

    echo ""
    echo "    --- Pedigree tests ---"
    echo "      AVG-based: $AVG_TEST  (CUN#1 closer to Kunenbo AND CUN#2 closer to Kishu?)"
    echo "      MAX-based: $MAX_TEST  (true pedigree indicator: closest hap matches pedigree?)"
    
    if [[ "$MAX_TEST" == "YES" ]]; then
      L3_OK="OK (MAX-based pedigree consistent)"
      if [[ "$AVG_TEST" == "NO" ]]; then
        L3_OK="WARN: MAX-based consistent, AVG-based inconsistent (likely due to haplotype leakage)"
      fi
    elif [[ "$MAX_TEST" == "NO" ]]; then
      L3_OK="WARN: MAX-based inconsistent (path-length effects, or CUN#1/CUN#2 assigned the wrong way round?)"
    else
      L3_OK="WARN: could not evaluate"
    fi
  else
    echo "    FAIL: odgi similarity failed or produced empty output"
    L3_OK="FAIL: similarity error"
  fi

  # ============================================================
  # Level 4: SV signal
  # ============================================================
  echo ""
  echo "  --- Level 4: SV signal (vg deconstruct -> VCF) ---"

  VCF="$QC_OUT/${CHR}.vcf"
  if [[ ! -f "$VCF" || ! -s "$VCF" ]]; then
    if vg deconstruct -P CUN#1 -H '#' -a -e -t 8 "$GFA" > "$VCF" 2>"$QC_OUT/${CHR}_vg_deconstruct.err"; then
      :
    else
      echo "    FAIL: vg deconstruct failed"
      tail -3 "$QC_OUT/${CHR}_vg_deconstruct.err" | sed 's/^/      /'
      L4_OK="FAIL: vg deconstruct error"
    fi
  fi

  if [[ -f "$VCF" && -s "$VCF" ]]; then
    N_TOTAL=$(grep -vc '^#' "$VCF")
    N_SNP=$(grep -v '^#' "$VCF" | awk 'length($4)==1 && $5 !~ /,/ && length($5)==1' | wc -l)
    N_SV=$(grep -v '^#' "$VCF" | awk '{
      ref_len=length($4); 
      n = split($5, alts, ",")
      for(i=1; i<=n; i++) { if((length(alts[i])-ref_len)>=50 || (ref_len-length(alts[i]))>=50) { print; next } }
    }' | wc -l)

    echo "    total variants : $N_TOTAL"
    echo "    SNPs           : $N_SNP"
    echo "    SVs (>=50bp)   : $N_SV"

    L4_OK="OK"
    [[ "$N_SV" -lt 100 ]] && L4_OK="WARN: SV count low ($N_SV)"
    [[ "$N_SV" -gt 20000 ]] && L4_OK="WARN: SV count high ($N_SV, possible noise)"
    echo "    Verdict        : $L4_OK"
  fi

  # ============================================================
  # Overall verdict
  # ============================================================
  # Rules:
  #   PASS: all levels OK
  #   FAIL: Level 1 or Level 2 FAIL (technical problem, graph unusable)
  #   WARN: only Level 3 or Level 4 has issues (graph usable, interpret carefully)
  echo ""
  L1_STATUS="$L1_OK"; L2_STATUS="$L2_OK"; L3_STATUS="$L3_OK"; L4_STATUS="$L4_OK"
  if [[ "$L1_STATUS" == FAIL* || "$L2_STATUS" == FAIL* ]]; then
    VERDICT="FAIL"
  elif [[ "$L1_STATUS" == "OK" && "$L2_STATUS" == "OK" && \
          "$L3_STATUS" == OK* && "$L4_STATUS" == "OK" ]]; then
    VERDICT="PASS"
  else
    VERDICT="WARN"
  fi
  echo "  ======================================================================"
  echo "   $CHR overall verdict: $VERDICT"
  echo "     Level 1 (structure)       : $L1_STATUS"
  echo "     Level 2 (path preservation): $L2_STATUS"
  echo "     Level 3 (pedigree)        : $L3_STATUS"
  echo "     Level 4 (SV signal)       : $L4_STATUS"
  echo "  ======================================================================"

  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$CHR" "$INPUT_BP" "$NODES" "$EDGES" "$PATHS" "$GRAPH_BP" \
    "$INFLATION" "$PATH_OK" "$N_SV" "$VERDICT" >> "$SUMMARY_TSV"
  
  echo ""
done

echo ""
echo "======================================================================"
echo " Overall summary ($SUMMARY_TSV)"
echo "======================================================================"
column -t -s $'\t' "$SUMMARY_TSV" 2>/dev/null || cat "$SUMMARY_TSV"

echo ""
echo " If all chromosomes are PASS, proceed to Chapter 7 (graph interpretation)"
echo " For WARN/FAIL chromosomes, check logs/ and parameters"

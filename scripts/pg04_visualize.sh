#!/usr/bin/env bash
# pg04_visualize.sh
#
# Purpose:
#   Visualize a pangenome graph produced by PGGB.
#   Generates:
#     1. 1D linearized view (odgi viz) - path layout along coordinate axis
#     2. 1D path-per-line view - detailed per-path visualization
#     3. 2D layout (odgi layout + draw) - global graph topology
#
# Usage:
#   bash pg04_visualize.sh . chr09

set -uo pipefail

PROJECT="${1:-.}"

# Project-local wrapper directory (see setup_singularity_wrappers.sh).
export PATH="${PROJECT}/bin:${PATH}"
CHR="${2:?chromosome name (e.g. chr09) required}"
OUT="${3:-${PROJECT}/results/viz}"
mkdir -p "$OUT"

OG=$(ls "${PROJECT}/results/pggb/${CHR}/${CHR}.fa.gz."*.smooth.final.og 2>/dev/null | head -1)

if [[ -z "$OG" || ! -f "$OG" ]]; then
  echo "ERROR: no OG file found for ${CHR} (pattern: chr${CHR}.fa.gz.*.smooth.final.og)"
  exit 1
fi

echo "======================================================================"
echo " PGGB Step 4: Graph visualization ($CHR)"
echo " Input: $(basename "$OG")"
echo "======================================================================"

# ---------- 1D View (minimal, most reliable command) ----------
echo ""
echo "-- 1D View (odgi viz) --"

VIZ_1D="$OUT/${CHR}_viz_1d.png"
echo "  Running: odgi viz -i OG -o $VIZ_1D -x 1500 -y 500"
if odgi viz -i "$OG" -o "$VIZ_1D" -x 1500 -y 500 2>&1 | tail -3; then
  if [[ -f "$VIZ_1D" ]]; then
    SIZE=$(du -h "$VIZ_1D" | cut -f1)
    echo "  OK: $VIZ_1D ($SIZE)"
  else
    echo "  FAIL: output file not created"
  fi
else
  echo "  FAIL: odgi viz command failed"
fi

# ---------- 1D View (path-per-line) ----------
echo ""
echo "-- 1D View (path-per-line detail) --"

VIZ_PATHS="$OUT/${CHR}_viz_paths.png"
echo "  Running: odgi viz -i OG -o $VIZ_PATHS -x 1500 -y 800 -P"
if odgi viz -i "$OG" -o "$VIZ_PATHS" -x 1500 -y 800 -P 2>&1 | tail -3; then
  if [[ -f "$VIZ_PATHS" ]]; then
    SIZE=$(du -h "$VIZ_PATHS" | cut -f1)
    echo "  OK: $VIZ_PATHS ($SIZE)"
  fi
fi

# ---------- 2D layout ----------
echo ""
echo "-- 2D layout (odgi layout + draw) --"

LAYOUT="$OUT/${CHR}.lay"
if [[ ! -f "$LAYOUT" ]]; then
  echo "  Step 1: odgi layout (this may take a while)"
  if odgi layout -i "$OG" -o "$LAYOUT" -t 8 2>&1 | tail -3; then
    echo "  OK: Layout file: $LAYOUT"
  else
    echo "  FAIL: odgi layout failed"
  fi
fi

if [[ -f "$LAYOUT" ]]; then
  LAYOUT_2D="$OUT/${CHR}_layout_2d.png"
  echo "  Step 2: odgi draw"
  if odgi draw -i "$OG" -c "$LAYOUT" -p "$LAYOUT_2D" -H 1000 -C 2>&1 | tail -3; then
    if [[ -f "$LAYOUT_2D" ]]; then
      SIZE=$(du -h "$LAYOUT_2D" | cut -f1)
      echo "  OK: $LAYOUT_2D ($SIZE)"
    fi
  fi
fi

echo ""
echo "======================================================================"
echo " Generated files"
ls -lh "$OUT"/${CHR}_*.png 2>/dev/null
echo "======================================================================"

echo ""
echo "-- Interpretation guide --"
echo ""
echo "  1D view (${CHR}_viz_1d.png) - check:"
echo "    OK: 6 paths run in horizontal parallel lines -> paths preserved"
echo "    OK: broad regions with consistent color -> well-shared sequences"
echo "    FAIL: each path appears as an independent parallel line -> graph not collapsed"
echo "    FAIL: narrow middle with thick ends -> shared regions only in center"
echo ""
echo "  2D layout (${CHR}_layout_2d.png) - check:"
echo "    OK: thick central 'backbone' visible -> shared core sequence"
echo "    OK: small bubbles/loops scattered -> small SVs"
echo "    FAIL: 6 tangled threads without clear backbone -> under-collapsed"
echo "    FAIL: no backbone visible -> parameter tuning required"

#!/usr/bin/env bash
# setup_singularity_wrappers.sh
#
# Purpose:
#   Create wrapper scripts in ~/bin/ so that tools inside the PGGB Singularity
#   image can be invoked from the host as if they were installed natively.
#
# Rationale:
#   The PGGB image (pggb_latest.sif) contains all tools needed for pangenome
#   construction and analysis:
#     pggb, wfmash, seqwish, smoothxg, odgi, vg, gfaffix, bcftools, samtools,
#     bgzip, tabix, multiqc, mash
#   Without wrappers, users would need to prefix every command with
#   `singularity exec pggb_latest.sif ...`, which is inconvenient. Wrappers
#   allow standard invocation like `odgi stats -i graph.og`.
#
# Usage:
#   bash setup_singularity_wrappers.sh /full/path/to/pggb_latest.sif
#
# Example:
#   bash setup_singularity_wrappers.sh /home/user/mikan/pggb_latest.sif

set -euo pipefail

SIF="${1:?Absolute path to Singularity image (.sif) required as argument}"

if [[ ! -f "$SIF" ]]; then
  echo "ERROR: SIF file not found: $SIF"
  exit 1
fi

SIF_ABS=$(readlink -f "$SIF")
echo "PGGB image: $SIF_ABS"

BIN_DIR="$HOME/bin"
mkdir -p "$BIN_DIR"

# List of tools packaged inside the PGGB image
TOOLS=(
  pggb          # Main pipeline driver
  wfmash        # All-vs-all pairwise aligner
  seqwish       # Alignment-to-graph induction
  smoothxg      # Graph normalization
  odgi          # Graph operations and QC
  vg            # Variation graph tools (deconstruct etc.)
  gfaffix       # GFA cleanup and normalization
  bcftools      # VCF manipulation
  samtools      # BAM/FASTA indexing
  bgzip         # BGZF compression
  tabix         # Index for compressed files
  multiqc       # Aggregated QC reports
  mash          # Sketch-based distance estimation
)

for TOOL in "${TOOLS[@]}"; do
  WRAPPER="$BIN_DIR/$TOOL"
  cat > "$WRAPPER" << EOF
#!/bin/bash
# Auto-generated Singularity wrapper for $TOOL
# Wraps: $SIF_ABS
exec singularity exec --bind "\$PWD" "$SIF_ABS" $TOOL "\$@"
EOF
  chmod +x "$WRAPPER"
  echo "  Created: $BIN_DIR/$TOOL -> singularity exec $SIF_ABS $TOOL"
done

echo ""
echo "======================================================================"
echo " Wrapper creation complete"
echo "======================================================================"
echo ""
echo " Check that ~/bin is in your PATH:"
if echo "$PATH" | tr ':' '\n' | grep -qx "$HOME/bin"; then
  echo "   OK: ~/bin is in PATH"
else
  echo "   NOTE: ~/bin is not in PATH"
  echo "     Add to .bashrc:"
  echo "       echo 'export PATH=\"\$HOME/bin:\$PATH\"' >> ~/.bashrc"
  echo "     Apply to current shell:"
  echo "       export PATH=\"\$HOME/bin:\$PATH\""
fi

echo ""
echo " Verification commands:"
echo "   which pggb odgi vg samtools bgzip"
echo "   pggb --version"
echo "   odgi version"
echo ""
echo " NOTE for SLURM jobs: the .bashrc file may not be sourced in SLURM"
echo " environments. Add 'export PATH=\"\$HOME/bin:\$PATH\"' inside your job"
echo " script to ensure wrappers are found."

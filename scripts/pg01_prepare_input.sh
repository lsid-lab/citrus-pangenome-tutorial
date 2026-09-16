#!/usr/bin/env bash
# pg01_prepare_input.sh
#
# Purpose:
#   Prepare the 6 haploid FASTA files (3 cultivars x 2 haplotypes) as input
#   for PGGB. This involves:
#     1. Renaming contigs to PanSN-spec format (sample#hap#chr)
#     2. Normalizing chromosome names to chr01..chr09
#        (input files may use inconsistent naming: chr01, chr1, CKUhap1_r1.0ch1, ...)
#     3. Grouping sequences by chromosome (each chr FASTA gets 6 sequences)
#     4. bgzip compression + samtools faidx indexing
#
# Input:
#   samplesheet.tsv with columns:
#     sample_id, cultivar_jp, species, pansn_prefix, hap1_path, hap2_path, source, pedigree
#
# Output:
#   03_pangenome/by_chr/chr01.fa.gz .. chr09.fa.gz
#
# Usage:
#   bash pg01_prepare_input.sh tables/samplesheet.tsv .

set -euo pipefail

export PATH="${HOME}/bin:/home/claude/bin:${PATH}"

SAMPLESHEET="${1:?samplesheet.tsv required}"
PROJECT="${2:-.}"
OUT="${PROJECT}/03_pangenome/by_chr"
mkdir -p "$OUT"

echo "======================================================================"
echo " PGGB Step 1: Input preparation (PanSN renaming + chromosome normalization)"
echo "======================================================================"

# ---------- Phase A: PanSN renaming + chromosome name normalization ----------
TMP="${OUT}/_tmp"
mkdir -p "$TMP"

# Robust processing via inline Python.
# Uses os.path (not pathlib) for maximum compatibility across Python 3 versions.
python3 - "$SAMPLESHEET" "$PROJECT" "$TMP" << 'PYEOF'
import sys
import re
import os
import gzip

samplesheet_path = sys.argv[1]
project_root = sys.argv[2]
tmp_dir = sys.argv[3]

# Debug: print environment info
print(f"[DEBUG] Python version: {sys.version.split()[0]}")
print(f"[DEBUG] samplesheet: {samplesheet_path}")
print(f"[DEBUG] project_root: {project_root}")
print(f"[DEBUG] tmp_dir: {tmp_dir}")


def extract_chr_number(contig_name):
    """
    Extract the chromosome number from a contig name and normalize to chrXX.
    Handles multiple naming conventions:
      "chr01"           -> "chr01"
      "chr1"            -> "chr01"
      "CKUhap1_r1.0ch1" -> "chr01"
      "ch09"            -> "chr09"
      "unplaced_00001"  -> None (unplaced, skipped)
    """
    # Try patterns in order of specificity
    patterns = [
        r'[Cc]hromosome[_.-]?(\d+)',
        r'[Cc]hr[_.-]?(\d+)',
        r'ch(\d+)',
    ]
    for pattern in patterns:
        matches = re.findall(pattern, contig_name)
        if matches:
            # Use the last match (handles cases where a numeric prefix
            # appears elsewhere in the name)
            num = int(matches[-1])
            if 1 <= num <= 20:
                return f"chr{num:02d}"
    return None


def open_fasta(path):
    """Open a FASTA file, transparently handling gzip compression."""
    if str(path).endswith('.gz'):
        return gzip.open(path, 'rt')
    return open(path, 'r')


def iter_fasta(path):
    """Yield (header, sequence_lines) tuples from a FASTA file."""
    with open_fasta(path) as f:
        header = None
        seq_lines = []
        for line in f:
            line = line.rstrip('\n').rstrip('\r')
            if line.startswith('>'):
                if header is not None:
                    yield header, seq_lines
                header = line[1:].split()[0]  # Take only the first token
                seq_lines = []
            else:
                seq_lines.append(line)
        if header is not None:
            yield header, seq_lines


# Read the samplesheet
with open(samplesheet_path) as f:
    header_line = f.readline().rstrip('\n').rstrip('\r').split('\t')
    col_idx = {name: i for i, name in enumerate(header_line)}
    required = ['sample_id', 'pansn_prefix', 'hap1_path', 'hap2_path']
    for req in required:
        if req not in col_idx:
            print(f"ERROR: samplesheet.tsv missing required column '{req}'",
                  file=sys.stderr)
            sys.exit(1)

    for line in f:
        line_clean = line.rstrip('\n').rstrip('\r')
        if not line_clean.strip():
            continue
        cols = line_clean.split('\t')
        sid = cols[col_idx['sample_id']]
        prefix = cols[col_idx['pansn_prefix']]
        h1_path = cols[col_idx['hap1_path']]
        h2_path = cols[col_idx['hap2_path']]

        for hap_idx, fa_relpath in [(1, h1_path), (2, h2_path)]:
            # Resolve path (support both absolute and relative)
            if os.path.isabs(fa_relpath):
                fa_path = fa_relpath
            else:
                fa_path = os.path.join(project_root, fa_relpath)

            # Explicit type check (defensive)
            if not isinstance(fa_path, str):
                print(f"[WARN] fa_path is not str: {type(fa_path)} ({fa_path})")
                fa_path = str(fa_path)

            if not os.path.isfile(fa_path):
                print(f"[SKIP] {sid}_hap{hap_idx}: file not found ({fa_path})")
                continue

            tag = f"{prefix}#{hap_idx}"
            out_fa = os.path.join(tmp_dir, f"{sid}_hap{hap_idx}.pansn.fa")

            unplaced_count = 0
            mapped_count = 0
            chr_seen = set()
            with open(out_fa, 'w') as out:
                for header, seq_lines in iter_fasta(fa_path):
                    chr_norm = extract_chr_number(header)
                    if chr_norm is None:
                        unplaced_count += 1
                        continue
                    if chr_norm in chr_seen:
                        # Duplicate chromosome number (unplaced fragment etc.)
                        # -> keep only the first (assumed longest) occurrence
                        continue
                    chr_seen.add(chr_norm)
                    new_header = f"{tag}#{chr_norm}"
                    out.write(f">{new_header}\n")
                    for sl in seq_lines:
                        out.write(sl + '\n')
                    mapped_count += 1

            print(f"  [PanSN+normalize] {sid}_hap{hap_idx}: {mapped_count} chromosomes "
                  f"kept, {unplaced_count} skipped (unplaced/other)")
            if mapped_count < 9:
                print(f"    WARN: expected 9 chromosomes but only {mapped_count} found")

print()
print("Phase A complete")
PYEOF

echo ""
echo "-- Intermediate files --"
ls -la "$TMP"/

echo ""
echo "-- Content check for each PanSN FASTA (first 3 sequence names) --"
for fa in "$TMP"/*.pansn.fa; do
  [[ ! -f "$fa" ]] && continue
  echo "  $(basename $fa):"
  grep '^>' "$fa" | head -3 | sed 's/^/    /'
  N=$(grep -c '^>' "$fa")
  echo "    (total: $N sequences)"
done

# ---------- Phase B: Group by chromosome ----------
echo ""
echo "======================================================================"
echo "-- Phase B: Group by chromosome (each chr FASTA = 6 sequences) --"
echo "======================================================================"

for CHR in chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09; do
  OUT_FA="$OUT/${CHR}.fa"
  : > "$OUT_FA"
  N_SEQ=0
  for fa in "$TMP"/*.pansn.fa; do
    [[ ! -f "$fa" ]] && continue
    # Exact-match on "#chr01$" at end of sequence name
    SEQ=$(seqkit grep -rp "#${CHR}$" "$fa" 2>/dev/null || true)
    if [[ -n "$SEQ" ]]; then
      echo "$SEQ" >> "$OUT_FA"
      N_SEQ=$((N_SEQ + 1))
    fi
  done

  if [[ "$N_SEQ" -eq 0 ]]; then
    echo "  ${CHR}: no sequences found, skipping"
    rm -f "$OUT_FA"
    continue
  fi

  SUM_MB=$(seqkit stats -T "$OUT_FA" 2>/dev/null | tail -1 | awk '{printf "%.1f", $5/1e6}')
  echo "  ${CHR}: ${N_SEQ} sequences, total ${SUM_MB} Mb"

  if [[ "$N_SEQ" -ne 6 ]]; then
    echo "    WARN: expected 6 sequences but got ${N_SEQ}. PGGB input will be incomplete."
    echo "    Sequences present:"
    seqkit seq -n "$OUT_FA" | sed 's/^/      /'
  fi

  # bgzip compression + faidx indexing (if tools available)
  if command -v bgzip &>/dev/null; then
    bgzip -f "$OUT_FA"
    if command -v samtools &>/dev/null; then
      samtools faidx "$OUT_FA.gz" 2>/dev/null || true
    fi
  fi
done

echo ""
echo "======================================================================"
echo " Output files"
echo "======================================================================"
ls -lh "$OUT"/chr*.fa.gz 2>/dev/null || ls -lh "$OUT"/chr*.fa 2>/dev/null || echo "(no output)"

echo ""
echo " Debug intermediate files: $TMP (remove when no longer needed)"
echo " Next step: bash scripts/pg02_run_pggb.sh"

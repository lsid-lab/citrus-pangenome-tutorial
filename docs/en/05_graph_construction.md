# Chapter 5. Pangenome graph construction

## Learning objectives

- Compare the major pangenome construction tools
- Understand why we choose **PGGB** for this tutorial
- Clarify the three meanings of "reference sequence"
- Learn PGGB's internal stages and key parameters
- Choose an execution environment and run PGGB

---

## 5.1 Major pangenome construction tools

There are three practical choices today.

### 5.1.1 Minigraph

- Developer: Heng Li (Broad Institute)
- Reference: Li 2020, *Genome Biology*
- Approach: uses a reference as the backbone and adds only **large SVs** (≥ 50 bp) from other samples
- Speed: **fastest**, suitable for large-scale datasets
- Characteristic: heavily reference-biased

**Best for**: dozens to hundreds of haploids where only an SV overview is needed. SNPs are not represented.

### 5.1.2 Minigraph-Cactus (MC)

- Developer: UCSC Genomics Institute
- Reference: Hickey et al. 2024, *Nature Biotechnology*
- Approach: builds a backbone with Minigraph, then adds multi-way alignment via Cactus
- Speed: slower than Minigraph but faster than PGGB
- Characteristic: explicit reference, captures both SVs and small variants

**Best for**: large-scale projects (e.g., HPRC's 47 humans) with reference-based interpretation as the goal.

### 5.1.3 PGGB (PanGenome Graph Builder)

- Developer: Erik Garrison et al.
- Reference: Garrison et al. 2024 preprint
- Approach: **all-vs-all alignment** (wfmash) → **graph induction** (seqwish) → **smoothing** (smoothxg)
- Speed: **slowest** but most comprehensive
- Characteristic: **reference-free** (all samples treated symmetrically), captures all variant types

**Best for**: small-to-mid scale (< 20 haploids), reference-free analysis, pedigree-focused studies.

### 5.1.4 Comparison

|Aspect|Minigraph|Minigraph-Cactus|PGGB|
|---|---|---|---|
|Reference dependence|Strong|Moderate|**None**|
|SVs|Good|Excellent|Excellent|
|SNPs|No|Yes|Yes|
|Speed|Fastest|Moderate|Slowest|
|Scale|Hundreds of haploids|Dozens|< 20|
|SNP-level QC|Not possible|Possible|**Yes**|

---

## 5.2 Three meanings of "reference sequence" — clearing up confusion

In pangenome discussions, "reference" is used in **three different senses**. Clearing this up early makes later discussion much easier.

### (1) Structural backbone (structural reference)

A sequence that forms the skeleton of the graph. Handling varies by tool.

- **Minigraph**: adds only SVs to the reference as branches; strongly reference-biased
- **Minigraph-Cactus**: uses the reference to build multi-way alignment; reference choice affects results
- **PGGB**: uses **reference-free** algorithm that treats all samples symmetrically. This meaning of "reference" does not apply.

### (2) Coordinate anchor (coordinate reference)

Once the graph exists, a coordinate system is needed to say things like "position 1234 on chromosome X". This is needed for VCF output and visualization.

**Even PGGB requires one coordinate reference**, though it does not affect the graph structure.

### (3) Comparative baseline (biological reference)

The baseline for statements like "sample A differs from B by X". Affects the framing of scientific claims.

### Applied to this tutorial

Because we use PGGB (which is reference-free):

- **(1) Structural backbone is not needed**
- **(2) Coordinate anchor**: chosen as **`CUN#1` (CUNphKu, Kunenbo-derived)** (reasoning below)
- **(3) Comparative baseline**: depends on context

**When we say "we chose a reference," we are actually choosing (2) coordinate anchor**, not the graph's backbone. Keep this distinction in mind.

### Why `CUN#1` as coordinate anchor?

1. **CUN is the pedigree center**: F1 of CKI × CKU, connecting the two parents in the graph
2. **Its parent of origin is known**: trio phasing identifies `CUN#1` as Kunenbo (CKU)-derived, so interpretation is unambiguous
3. **Newest and highest-quality of the Isobe 2023 assemblies**

---

## 5.3 Why PGGB?

We choose **PGGB** in this tutorial. Reasons:

### Reason 1: Reference-free suits pedigree-centered observation

Selecting a specific cultivar as backbone would bias interpretation of the family relationship. PGGB treats all six haploids equally.

### Reason 2: Six haploids is computationally realistic

PGGB is slow (5-10 hours per chromosome), but 6 haploids × 30-40 Mb/chr is feasible on 32 CPUs. Fits in a student environment.

### Reason 3: Captures SNP-level variation

We want small variants (SNP, indel) too, so Minigraph is unsuitable. MC would work, but PGGB is more comprehensive.

### Reason 4: Educational transparency

PGGB's three stages (wfmash → seqwish → smoothxg) have clear boundaries, making it **easy for students to understand what happens inside**. Also helpful for debugging and interpretation.

---

## 5.4 PGGB internal stages

PGGB is a three-stage pipeline:

```
[Input] 6 haploid FASTAs (bgzip, PanSN-named)
         │
         ▼
┌────────────────┐
│  1. wfmash     │  All-vs-all pairwise alignment
│  (PAF output)  │  (which regions align to which)
└────────────────┘
         │
         ▼
┌────────────────┐
│  2. seqwish    │  Alignment → graph induction
│  (GFA output)  │  (shared regions merged as nodes)
└────────────────┘
         │
         ▼
┌────────────────┐
│  3. smoothxg   │  Graph normalization and simplification
│  (final GFA)   │  (POA-based smoothing)
└────────────────┘
         │
         ▼
[Output] final.gfa, final.og, VCF
```

---

## 5.5 Important PGGB parameters

Concretizing the "parameters to watch" concept from Chapter 1 as PGGB options.

### 5.5.1 Key parameters

|Parameter|Meaning|Recommended|Explanation|
|---|---|---|---|
|`-n`|Number of haploids|**6**|3 cultivars × 2 hap|
|`-p`|Percent identity|**95**|Threshold for alignment (≥90% identity)|
|`-s`|Segment length (bp)|**10000**|Minimum alignment unit for wfmash|
|`-V`|VCF reference path|`CUN#1:#`|Coordinate anchor for vg deconstruct (§5.2)|
|`-Y`|PanSN separator|`#`|Delimiter in PanSN names|
|`-t`|Thread count|32-48|Depends on your CPU allocation|

### 5.5.2 Rationale for values

**`-p 95` (percent identity)**
- Citrus intra-specific SNP density is 17/kbp (Kiryu 2026), i.e., 98.3% identity
- `-p 95` allows alignment where identity ≥ 95% over a window
- `-p 90` also works, but 95 yields fewer false positives and runs faster

**`-s 10000` (segment length)**
- Minimum **alignment unit** in wfmash
- Smaller = more detail but more noise
- 10 kb is reasonable for citrus HiFi assemblies

**`-n 6` (haploid count)**
- Must **exactly match** the input haploid count
- wfmash aligns each sequence against the top n-1 most similar sequences

### 5.5.3 Version-dependent note: the `-B` parameter

This tutorial targets **PGGB from the Singularity image `pggb:latest` (revision `4225c6c` approx.)**. In this version, the default value of `-B` (transclose-batch) fails to parse correctly, and seqwish silently falls back to a smaller default. This can cause **under-compression** of the graph.

**Workaround**: always specify `-B` explicitly at runtime, e.g. `-B 1G`. As a guideline, use a value at least as large as the total input length. The largest chromosome (chr03) is 306.7 Mb, so `-B 1G` provides ample headroom.

**Future PGGB versions** may fix this. In Chapter 6 we check the `.params.yml` file to confirm the `transclose-batch` value that was actually used.

Credit: this issue was diagnosed and documented by a collaborator through systematic parameter experiments.

---

## 5.6 Execution environment options

PGGB has many dependencies, so **containerized execution is strongly recommended**. Three options are common.

### Option A: Docker

The most commonly documented approach on the web.

```bash
docker pull ghcr.io/pangenome/pggb:latest
docker run -v $PWD:/data ghcr.io/pangenome/pggb pggb ...
```

**Good for**: personal workstations with root privilege

**Note**: HPC environments often disallow Docker due to security concerns (running as root on multi-user hosts)

### Option B: Conda / mamba

Direct installation via a package manager.

```bash
mamba create -n pggb -c bioconda pggb
mamba activate pggb
```

**Good for**: personal environments, integration with other software

**Note**: dependency resolution can be fragile

### Option C: Singularity (chosen for this tutorial)

The de facto standard on HPC. Runs Docker images without requiring root.

```bash
singularity pull docker://ghcr.io/pangenome/pggb:latest
# -> pggb_latest.sif is created
```

**Good for**: HPC (SLURM/PBS), multi-user shared environments

**Note**: to call tools inside the SIF from the host, wrappers are helpful (see next section)

**This tutorial uses Option C (Singularity)**. The rest of this chapter assumes it.

---

## 5.7 Execution (Singularity)

### 5.7.1 Get the image and create wrappers

```bash
# Download image
singularity pull docker://ghcr.io/pangenome/pggb:latest

# Create wrappers so tools inside SIF can be called from the host
bash scripts/setup_singularity_wrappers.sh $(readlink -f pggb_latest.sif)
```

This creates `~/bin/{pggb,odgi,vg,samtools,bcftools,bgzip,...}`.

**Add to PATH**:

```bash
echo 'export PATH="$HOME/bin:$PATH"' >> ~/.bashrc
export PATH="$HOME/bin:$PATH"
```

**Verify**:

```bash
pggb --version
odgi version
```

### 5.7.2 Input preparation

```bash
bash scripts/pg01_prepare_input.sh tables/samplesheet.tsv .
```

This script:
1. Renames each haploid FASTA using **PanSN** format
   e.g., `CKUhap1_r1.0ch1` → `CKU#1#chr01`
2. Normalizes chromosome numbers to `chr01`..`chr09`
3. Groups sequences by chromosome
4. bgzip compresses and creates `.fai` indexes

Output: `03_pangenome/by_chr/chr01.fa.gz` .. `chr09.fa.gz`

Each per-chromosome FASTA contains **exactly six sequences** (3 cultivars × 2 haps).

### 5.7.3 Test run on the smallest chromosome

Start with chr09 (smallest) for a test:

```bash
bash scripts/pg02_run_pggb.sh . 32 chr09 chr09
```

Takes 30-60 minutes.

### 5.7.4 Parallel run of all chromosomes (SLURM)

Submit **individual sbatch jobs** for chr01..chr09 (rather than a job array). Job arrays would apply the same resource allocation to all chromosomes, whereas individual jobs allow flexible resource tuning.

Per-chromosome job script (`run_pggb_chr.sbatch`):

```bash
#!/bin/bash
#SBATCH --job-name=pggb_chr
#SBATCH --output=logs/pggb_%x_%j.out
#SBATCH --error=logs/pggb_%x_%j.err
#SBATCH --mem=64G
#SBATCH --cpus-per-task=32
#SBATCH --time=12:00:00

export PATH="$HOME/bin:$PATH"

CHR=${1:?please specify chromosome name (e.g. chr01)}
bash scripts/pg02_run_pggb.sh . 32 $CHR $CHR
```

Submit all nine:

```bash
for CHR in chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09; do
  sbatch run_pggb_chr.sbatch $CHR
done
```

Nine jobs run in parallel; total wall time is approximately 8-12 hours.

### 5.7.5 Time and memory guide (32 CPU)

|Chromosome|Input size|Time|MaxRSS|
|---|---|---|---|
|chr01|195.7 Mb|~28 min|12 GB|
|chr02|213.5 Mb|~44 min|6 GB|
|chr03|306.7 Mb|~62 min|18 GB|
|chr04|177.9 Mb|~19 min|5 GB|
|chr05|252.6 Mb|~43 min|14 GB|
|chr06|165.8 Mb|~24 min|10 GB|
|chr07|210.8 Mb|~74 min|6 GB|
|chr08|228.8 Mb|~56 min|13 GB|
|chr09|195.9 Mb|~29 min|4 GB|

These are collaborator-measured values. Peak memory is at most ~18 GB, so `--mem=64G` provides comfortable headroom.

---

## 5.8 Output files

For each chromosome (chr09 example):

```
03_pangenome/by_chr/chr09_pggb/
├── chr09.fa.gz.<h1>.<h2>.<h3>.smooth.final.gfa    ← final graph (GFA)
├── chr09.fa.gz.<h1>.<h2>.<h3>.smooth.final.og     ← odgi format
├── chr09.fa.gz.<h1>.<h2>.<h3>.smooth.final.CUN#1.vcf  ← VCF
├── chr09.fa.gz.<h1>.<h2>.<h3>.smooth.final.og.lay.draw_multiqc.png  ← 2D image
├── chr09.fa.gz.<h1>.<h2>.<h3>.smooth.*.params.yml  ← authoritative parameters
├── chr09.fa.gz.<h1>.<h2>.<h3>.smooth.*.log         ← log
└── multiqc_report.html                              ← aggregated report
```

**Hashes**: PGGB embeds parameter hashes in filenames so that different parameter runs do not overwrite each other.

---

## Chapter summary

- Three pangenome tools: Minigraph, Minigraph-Cactus, PGGB
- "Reference" has three meanings; **PGGB requires choosing only (2) coordinate anchor**
- PGGB is a 3-stage pipeline: wfmash → seqwish → smoothxg
- Parameters follow the genetic distance of the cultivars: citrus uses `-p 95 -s 10000 -n 6`
- In the current PGGB version, `-B` must be explicit; verify via `.params.yml` after the run
- Singularity is the practical choice on HPC
- Submit chromosomes as individual sbatch jobs, not job arrays

The next chapter evaluates the **quality of the resulting graph**.

## References

- Garrison E, et al. (2024). Building pangenome graphs. *bioRxiv*.
- Li H (2020). The design and construction of reference pangenome graphs with minigraph. *Genome Biol* 21:265.
- Hickey G, et al. (2024). Pangenome graph construction from genome alignments with Minigraph-Cactus. *Nat Biotechnol* 42:663-673.
- Kiryu Y, et al. (2026). AlleleMiner. *DNA Res* 33:dsag004.

---

[← Chapter 4](04_data_qc.md) | [Chapter 6: Graph QC →](06_graph_qc.md)

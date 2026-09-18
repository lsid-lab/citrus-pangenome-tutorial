# Chapter 5. Pangenome graph construction

## Learning objectives

- Compare the major pangenome construction tools
- Understand why we choose **PGGB** here
- Learn the three meanings of "reference sequence" and which of them PGGB needs
- Learn PGGB's internal stages and key parameters, and **how to choose parameters with no prior knowledge of your species**
- Choose an execution environment and run PGGB
- Understand **impg** and the idea of separating alignment from graph induction

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

**Best for**: work where reference-coordinate interpretation is the point, where you need a graph for read mapping, or large projects (HPRC and similar).

> A worked procedure is in [Appendix: the same data with Minigraph-Cactus](appendix_minigraph_cactus.md).

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
|Reference sequence required|**Yes**|**Yes**|**No**|
|SVs|Good|Excellent|Excellent|
|SNPs|No|Yes|Yes|
|Speed|Fastest|Moderate|Slowest|
|Scale|Hundreds of haploids|Dozens|< 20|
|SNP-level QC|Not possible|Possible|**Yes**|

---

## 5.2 Three meanings of "reference sequence", and which one PGGB needs

In pangenome discussions the word "reference" is used in **three different senses**. Conflating them makes conversations go past each other, so distinguish them first.

|Sense|What it is|
|---|---|
|**(1) Structural backbone**|The sequence forming the skeleton of the graph. Sequence absent from it is attached as branches|
|**(2) Coordinate anchor**|The coordinate system used to say "position 1,234 of chr09". Needed for VCF and visualization|
|**(3) Comparative baseline**|The baseline for "A differs from B by X". A matter of how claims are framed|

Tools differ in whether they require (1). Minigraph and Minigraph-Cactus **require it**, and the reference chosen shapes the graph itself.

### In PGGB

**PGGB has no (1).** All sequences are treated symmetrically, so **the structure of the graph does not depend on which sequence you call "the reference".**

**(2), however, PGGB still needs one of.** To write "position 1,234 of chr09 is A→G" in a VCF, you must fix one path to count those 1,234 bases along. In PGGB that is what `-V` selects.

This tutorial uses **`CUN#1` (CUNphKu, Kunenbo-derived)**. Given `-V CUN#1`, pggb internally calls `vg deconstruct -P CUN#1#` and writes the VCF in CUN#1 coordinates. The reasons are practical:

- **It sits at the centre of the pedigree**: Satsuma is the F1, so it compares naturally against either parent
- **Its parent of origin is known**: trio phasing identifies `CUN#1` as Kunenbo (CKU)-derived (§2.5)
- **Chromosome-scale with little missing sequence**

Choosing `CUN#2` instead gives **the same graph**, only different VCF coordinates and REF column. That is what "having no (1)" means in practice. Chapters 6 and 7 are all written against `CUN#1`, so if you change it, change both `-V` and `vg deconstruct -P`.

Sense (3) is about how an argument is framed rather than about the analysis, and it comes back when we interpret results in Chapter 7.

## 5.3 Why we use PGGB, and what MC would give you

This tutorial uses **PGGB**.

### Reasons for PGGB

**1. All haploids are treated symmetrically**

Reference-free construction matches this tutorial's framing of looking at three cultivars as equals. Picking one as the backbone (sense (1) in §5.2) makes "sequence present in the backbone" and "sequence absent from it" asymmetric.

**2. The internal stages can be inspected separately**

wfmash → seqwish → smoothxg are cleanly separated, and you can look at each intermediate (the PAF, the freshly induced GFA). **Being able to follow what is happening** matters for a tutorial — and it is why a problem like `-B` in §5.5.3 is noticeable at all.

**3. Six haploids is computationally realistic**

PGGB is the slowest of the three, but six haploids at 30-40 Mb per chromosome take 30-60 minutes on 32 CPUs. That runs in a student environment.

### What Minigraph-Cactus would give you

Running the same data through MC instead buys the following. **For some purposes it is the better choice.**

**1. The coordinate system is fixed from the start**

MC **requires** a reference genome, and that reference is never clipped or cyclic and **directly defines the VCF coordinate system and the chromosome decomposition** (serving as both (1) and (2) in §5.2). If you want to work in one fixed coordinate system throughout, that is more direct than selecting one afterwards with `-V`.

**2. Read-mapping artifacts come out of the box**

Alongside GFA and VCF you get **GBZ and the `.dist` / `.min` / `.hapl` indexes, ready for read mapping with vg Giraffe**. That suits work like the v1/v2 comparison in §7.6, which validates against independent read data mapped onto the graph.

**3. Softmasked input does not matter**

MC does not require input softmasking. Here, r2.0 arrives repeat-masked, which is why `pg01_prepare_input.sh` gained an upper-casing step (§7.6); with MC that handling is unnecessary.

**4. It is faster**

Designed for closely related samples of the same species, it finishes sooner than PGGB at this scale.

**In practice, choose by working backwards from the artifact you need**: MC if you want a VCF in reference coordinates and indexes for read mapping, PGGB if you want every sample treated symmetrically, minigraph if you only want an overview of large SVs.

> The procedure for running the same six haploids through MC is in [Appendix: the same data with Minigraph-Cactus](appendix_minigraph_cactus.md). MC emits GFA v1.1 (`W` lines), so the `vg convert` downgrade from §1.2.2 is genuinely needed there.

Note also that **even with PGGB, current practice is not "PGGB alone"**. Pairing it with **impg**, which separates alignment from graph induction, is the current approach; see §5.9.

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
|`-V`|VCF coordinate reference|`CUN#1`|§5.2. The spec is `REF[:LEN]`|
|`-t`|Thread count|32-48|Depends on your CPU allocation|

> **Mind the `-V` spec**: it is `REF[:LEN]`. pggb appends the PanSN `#` to REF itself, so `-V CUN#1` runs `vg deconstruct -P CUN#1#`. **What follows the colon is not a delimiter but a number (LEN)**; a value greater than 0 additionally emits a `vcfbub` + `vcfwave` decomposed VCF. We do not use the decomposed output.
>
> **`-Y` is not passed.** It does not tell pggb the PanSN delimiter; it is `--exclude-delim`, "skip mappings when the query and target share the prefix before the **last** occurrence of the character". Its default is `#`, and with full PanSN names the group is `sample#hap`, so `CUN#1` and `CUN#2` are in different groups and are aligned to each other as intended. The default is what we want, so we leave it alone.

### 5.5.2 How to choose the parameters

This is where people actually get stuck. Published values exist for citrus, but **what do you do when you have no prior knowledge of your species?** That first.

#### `-p` (percent identity) — measure it from your own data

`-p` is the threshold above which a stretch of identity is considered for alignment. Rather than taking a number from the literature, you can **measure the similarity among your own sequences** and decide from that.

```bash
# All-pairs distance among the 6 haploids (fast, alignment-free)
mash triangle data/input/chr09.fa.gz > chr09.mash.tri

# Or run only wfmash's mapping stage and look at the identity distribution
wfmash -m -t 16 data/input/chr09.fa.gz \
  | awk '{for(i=13;i<=NF;i++) if($i ~ /^(id|gi):f:/) {split($i,a,":"); print a[3]}}' \
  | sort -n | uniq -c
```

The rule of thumb is to **sit a few points below the lowest identity among the pairs you want aligned**. Too low and unrelated regions become alignment candidates, inflating runtime; too high and genuinely homologous but divergent regions are dropped.

For citrus, intra-specific SNP density is around 17/kbp (Kiryu et al. 2026), i.e. roughly **98.3% identity**. From that, `-p 95` has comfortable headroom on the side of not dropping anything. **That one sentence is the answer the commands above would have given you.**

#### `-s` (segment length) — start from the default, then sweep one chromosome

`-s` is the shortest unit at which wfmash looks for mappings.

- **Smaller** picks up finer homology, but produces seeds that are not unique inside repeats, adding noise and runtime
- **Larger** is more stable, but misses rearrangements shorter than the segment

**10 kb is not a citrus-specific figure**; it is a widely used starting point for chromosome-scale assemblies. To check it for your species, **sweep `-s` on one cheap chromosome** and compare (chr09 for citrus — it has the lowest measured memory use in §5.7.5):

```bash
for S in 5000 10000 20000; do
  pggb -i data/input/chr09.fa.gz -o results/pggb/chr09_s${S} \
       -p 95 -s $S -n 6 -B 1G -t 32
done

# Compare compression and node count (the metrics from Chapter 6, §6.3)
for S in 5000 10000 20000; do
  odgi stats -i results/pggb/chr09_s${S}/*.smooth.final.og -S
done
```

Take **the largest `-s` at which compression does not degrade badly and the node count does not balloon**.

#### `-n` (haploid count) — optional in current pggb

Current pggb **counts the haplotypes itself from the PanSN names** (distinct `sample#hap`) when `-n` is not given.

So `-n 6` is not required; we pass it **as a check on the automatic count**. If the number disagrees, input preparation went wrong, and this is where you find out.

### 5.5.3 Version-dependent note: the `-B` parameter

This tutorial targets **PGGB from the Singularity image `pggb:latest` (revision `4225c6c` approx.)**. In this version, the default value of `-B` (transclose-batch) fails to parse correctly, and seqwish silently falls back to a smaller default. This can cause **under-compression** of the graph.

**Workaround**: always specify `-B` explicitly at runtime, e.g. `-B 1G`. As a guideline, use a value at least as large as the total input length. The largest chromosome (chr03) is 306.7 Mb, so `-B 1G` provides ample headroom.

**Future PGGB versions** may fix this. In Chapter 6 we check the `.params.yml` file to confirm the `transclose-batch` value that was actually used.

Credit: this issue was diagnosed and documented by collaborator **Shuto Machida** through systematic parameter experiments.

**The shape of the pipeline is what makes this painful.** PGGB runs alignment, graph induction and smoothing in one command. Realising afterwards that a graph-induction parameter was wrong therefore means **redoing the wfmash alignment, by far the most expensive stage**. Reducing that cost is exactly what impg (§5.9) is about.

---

## 5.6 Execution environment options

PGGB has many dependencies, so a **containerised environment is strongly recommended**.

### Option A: Singularity / Apptainer (used in this tutorial)

The standard on HPC. It runs Docker images directly and **needs no privileges**.

```bash
singularity pull docker://ghcr.io/pangenome/pggb:latest
# -> creates pggb_latest.sif

# Run directly (the current directory is mounted automatically)
singularity exec pggb_latest.sif pggb --version
```

**Best for**: HPC (SLURM/PBS), shared multi-user systems. This tutorial assumes it.

### Option B: Docker

```bash
docker pull ghcr.io/pangenome/pggb:latest

docker run --rm -it \
  -u "$(id -u):$(id -g)" \
  -v "$PWD":/data -w /data \
  ghcr.io/pangenome/pggb pggb --version
```

**What you need is not root but membership of the `docker` group** (or a rootless Docker setup). That said, membership of `docker` is effectively equivalent to root, which is why **shared HPC systems usually do not grant it**.

Without `-u "$(id -u):$(id -g)"`, **output files are created owned by root and you cannot delete them afterwards**. Do not omit it. `--rm` discards the container on exit; `-w /data` sets the working directory.

### Option C: Podman

Nearly the same commands as Docker, and **rootless by default**. It needs no daemon, so it is sometimes available on HPC.

```bash
podman run --rm -it \
  -v "$PWD":/data -w /data \
  ghcr.io/pangenome/pggb pggb --version
```

Under rootless podman, root inside the container maps to your own user outside it, so `-u` is usually unnecessary and output belongs to you.

### Option D: Conda / mamba

```bash
mamba create -n pggb -c conda-forge -c bioconda pggb
mamba activate pggb
```

**Caveat**: the dependency set is complex and environment solving sometimes fails. Unless you need to combine PGGB with other tools, a container is more reliable.

---

## 5.7 Running PGGB

### 5.7.1 Getting the image

```bash
singularity pull docker://ghcr.io/pangenome/pggb:latest
SIF=$(readlink -f pggb_latest.sif)

singularity exec "$SIF" pggb --version
singularity exec "$SIF" odgi version
```

**Every command below can be run as `singularity exec "$SIF" <tool> ...`.** If typing that each time is tiresome, there is a script that generates wrappers (optional):

```bash
bash scripts/setup_singularity_wrappers.sh "$SIF"
export PATH="$PWD/bin:$PATH"     # this session only
```

That creates `./bin/{pggb,odgi,vg,samtools,bgzip,...}` so you can write `pggb ...` directly. **It writes only to `./bin` inside the project** — never to `~/bin` or your shell rc files. If you want it permanently, add the `export PATH` to your own `~/.bashrc` yourself.

The commands below use the explicit `singularity exec` form.

### 5.7.2 Preparing the input

PanSN renaming is **preprocessing**, so its output is not a result but an **input**, and it goes to `data/input/`.

```bash
bash scripts/pg01_prepare_input.sh tables/samplesheet.tsv .
```

What the script does:

1. **Renames contigs to PanSN** — e.g. `CKUhap1_r1.0ch1` → `CKU#1#chr01`
2. Normalises chromosome names to `chr01`–`chr09`
3. **Upper-cases the sequence** (so repeat-masked r2.0 reaches the aligner on the same terms; §7.6)
4. Splits by chromosome
5. bgzip + `samtools faidx`

Output: `data/input/chr01.fa.gz` … `chr09.fa.gz`, each holding **exactly 6 sequences**.

> If `seqkit`, `bgzip` or `samtools` is missing, the script now **fails immediately and says so**. It used to skip silently, exit 0 without producing `.fa.gz` or `.fai`, and fail much later inside pggb.

Check:

```bash
singularity exec "$SIF" seqkit seq -n data/input/chr09.fa.gz
# CUN#1#chr09
# CUN#2#chr09
# CKI#1#chr09
# ...
```

### 5.7.3 Running one chromosome

**This is the actual work.** We use chr09 here:

```bash
singularity exec "$SIF" pggb \
  -i data/input/chr09.fa.gz \
  -o results/pggb/chr09 \
  -p 95 \
  -s 10000 \
  -n 6 \
  -B 1G \
  -V CUN#1 \
  -t 32 \
  -m -S
```

Option by option:

|Option|Meaning|
|---|---|
|`-i`|Input FASTA (bgzip + faidx)|
|`-o`|Output directory|
|`-p 95`|Mapping identity threshold (§5.5.2)|
|`-s 10000`|Segment length (§5.5.2)|
|`-n 6`|Haploid count; auto-detected from PanSN if omitted (§5.5.2)|
|`-B 1G`|seqwish transclose batch. **Must be explicit** (§5.5.3)|
|`-V CUN#1`|VCF coordinate reference (§5.2)|
|`-t 32`|Threads|
|`-m`|Emit statistics for MultiQC|
|`-S`|Emit per-stage statistics|

It takes 30-60 minutes. When it finishes, always:

```bash
# Confirm the parameters actually took effect, from pggb's own record
grep -E 'transclose-batch|map-pct-id|segment-length|n-haplotypes' \
  results/pggb/chr09/*.params.yml
```

**Read what the tool wrote out, not what you typed on the command line.** Given the `-B` issue in §5.5.3, make this a habit.

### 5.7.4 All chromosomes

Just the same command, once per chromosome.

```bash
for CHR in $(ls data/input/*.fa.gz | xargs -n1 basename | sed 's/\.fa\.gz$//'); do
  singularity exec "$SIF" pggb \
    -i "data/input/${CHR}.fa.gz" \
    -o "results/pggb/${CHR}" \
    -p 95 -s 10000 -n 6 -B 1G -V CUN#1 -t 32 -m -S \
    > "results/pggb/logs/${CHR}.log" 2>&1
done
```

`scripts/pg02_run_pggb.sh` is that loop plus some bookkeeping (input checks, skipping completed chromosomes, a success/failure summary). **It takes the chromosome names from the files in `data/input/`, so `chr01..chr09` is not assumed anywhere.**

```bash
bash scripts/pg02_run_pggb.sh . 32            # everything in data/input
bash scripts/pg02_run_pggb.sh . 32 chr09      # just chr09
```

**To run them in parallel under SLURM**, a job array is the simplest route.

`run_pggb_array.sbatch`:

```bash
#!/bin/bash
#SBATCH --job-name=pggb
#SBATCH --output=results/pggb/logs/%x_%A_%a.out
#SBATCH --error=results/pggb/logs/%x_%A_%a.err
#SBATCH --mem=24G
#SBATCH --cpus-per-task=32
#SBATCH --time=12:00:00

SIF=${1:?give the absolute path to pggb_latest.sif}

CHRS=($(ls data/input/*.fa.gz | xargs -n1 basename | sed 's/\.fa\.gz$//'))
CHR=${CHRS[$((SLURM_ARRAY_TASK_ID - 1))]}

singularity exec "$SIF" pggb \
  -i "data/input/${CHR}.fa.gz" \
  -o "results/pggb/${CHR}" \
  -p 95 -s 10000 -n 6 -B 1G -V CUN#1 -t 32 -m -S
```

```bash
mkdir -p results/pggb/logs
N=$(ls data/input/*.fa.gz | wc -l)
sbatch --array=1-${N}%3 run_pggb_array.sbatch "$SIF"
```

**`%3` caps how many tasks run at once.** Launching all nine together asks for 9 × 32 = 288 CPUs simultaneously. Even where no per-user limit is configured, it is courteous not to occupy a shared cluster single-handed. Set it to suit your allocation and how busy the queue is.

**`--mem` is shared by every task in the array**, so size it for the heaviest chromosome. The measurements in §5.7.5 top out at 18 GB for chr03, hence 24 GB.

> If you do want per-chromosome resource requests — for instance where your site penalises large memory requests — submit one `sbatch` per chromosome instead of using `--array`. Within the range measured here, there is no need to.

Expect 8-12 hours of wall-clock time (on the environment described in §5.7.5).

### 5.7.5 Time and memory guide

These are **measurements from one specific environment**. Check the conditions before applying the numbers to your own setup.

#### Measurement environment

<!-- TODO: add the CPU model, how time and MaxRSS were obtained (sacct or /usr/bin/time), the measurement date, and who ran it. -->

|Item|Value|
|---|---|
|Node|1 socket × 64 physical cores, SMT enabled (ThreadsPerCore=2), 128 logical CPUs|
|Memory|503 GiB|
|Node usage|**Shared** with other jobs|
|**Output filesystem**|**NFSv4.1** (network filesystem, rsize/wsize 1 MB)|
|Scheduler|Slurm 23.02.2|
|Container|SingularityCE 3.11.3|
|OS|Ubuntu 22.04.5 LTS (kernel 6.8.0)|
|pggb|`ghcr.io/pangenome/pggb`, revision `4225c6c`|
|Parameters|`-p 95 -s 10000 -n 6 -B 1G -V CUN#1 -t 32`|


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

These are collaborator-measured values. Peak memory is at most ~18 GB, so the `--mem=24G` in §5.7.4 has comfortable headroom.

#### How to read these numbers

**1. Runtime depends heavily on where the output goes.**

Unless `-D` / `--temp-dir` is given, pggb **writes its intermediates straight into the output directory**. The measurements above had that directory on NFS, so wfmash's PAF and seqwish's intermediates were all read and written across the network. If you can put the temporary files on node-local SSD or NVMe, expect it to be faster.

```bash
singularity exec "$SIF" pggb \
  -i data/input/chr09.fa.gz -o results/pggb/chr09 \
  -D "${TMPDIR:-/tmp}/pggb_chr09" \
  -p 95 -s 10000 -n 6 -B 1G -V CUN#1 -t 32 -m -S
```

> **Do not put a comma in the `-D` path.** seqwish reads commas as PAF-file separators, and pggb rejects it.

**2. "32 CPU" means 32 logical CPUs.**

SMT was enabled on the measurement node, so `-t 32` may have corresponded to 16 physical cores (depending on how Slurm allocates). With 32 physical cores available, expect shorter times.

**3. The node was shared.**

I/O and memory bandwidth were shared with other jobs. A dedicated node may be faster.

So the table is neither a floor nor a ceiling — it is **one observation: "on NFS, on a shared node, at 32 logical CPUs, this is what happened"**. Measure one chromosome in your own environment before planning the rest.

---

## 5.8 Output files

For each chromosome (chr09 example):

```
results/pggb/chr09/
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

## 5.9 Reusing the alignment — impg

Everything above runs **PGGB on its own**. Current large-scale pangenome construction — the HPRC release 2 build, for instance — **has adopted the combination of PGGB and [impg](https://github.com/pangenome/impg)**. We do not run it here, but the idea is worth knowing.

### What it does

impg treats **the all-vs-all alignment (PAF) itself as an implicit pangenome graph**, walking over the alignments directly instead of materialising a graph.

```bash
# 1. Index the alignment
impg index -a aln.paf -i aln.impg

# 2. Pull the sequence homologous to a region out of every haplotype
impg query -a aln.paf -r CUN#1#chr09:1000000-1200000 -d 100 -x

# 3. Cut windows out of the alignment network (1 Mb windows, absorbing gaps up to 100 kb)
impg partition -a aln.paf -w 1000000 -d 100000

# 4. Lace the per-window graphs into one
impg lace -l gfa_list.txt -o combined.gfa
```

The HPRC v2 build is broadly **`impg partition` over the whole cohort → PGGB per partition → `impg lace`**.

### The benefit is not only scale

It is easy to read this as "the data is big, so split it up". **The substance is the separation of alignment from graph induction**, and two of its benefits apply at six haploids too.

**1. You do not redo the expensive stage**

Of PGGB's three stages, wfmash's alignment dominates the cost. impg keeps the **PAF as a first-class, indexed object**, so you can **align once and rebuild graphs with different seqwish or smoothxg parameters without realigning**.

Recall `-B` from §5.5.3. Forget it, get an under-compressed graph, and with PGGB alone you start again **from the alignment**. With the alignment kept as a separate artifact, the cost is only the graph induction. The `-s` sweep suggested in §5.5.2 is much cheaper in this form too.

**2. You can work one locus at a time**

`impg query` **extracts the sequence homologous to a region from every haplotype**. You can build a small graph for just the gene or QTL interval you care about and tune parameters there, without holding a whole-chromosome graph.

**3. It does not assume chromosome correspondence**

`impg partition` derives its windows **dynamically from transitive homology in the alignment network**. Nothing assumes that chr01 corresponds to chr01.

That is a third answer to the question in §1.4 — build per chromosome, or over the whole genome? Per-chromosome assumes the correspondence; whole-genome demands enormous resources. **Cutting windows from homology escapes both constraints.** For targets with interchromosomal rearrangement, such as cancer genomes, this matters regardless of size.

### Why we do not use it here

At six haploids and one chromosome, PGGB alone finishes in 30-60 minutes, and a single self-contained command is easier to follow than a split pipeline. **For serious work on your own data, though, look at impg first.**

---

## Chapter summary

- Three pangenome tools: Minigraph, Minigraph-Cactus (MC), PGGB
- "Reference" has three senses — **(1) backbone / (2) coordinate anchor / (3) comparative baseline**. **PGGB requires only (2)**
- MC would give you a fixed reference coordinate system, Giraffe-ready indexes and no softmask handling. **Choose by working backwards from the artifact you need**
- PGGB's construction is reference-free; **one coordinate path is chosen only for the VCF** (`-V CUN#1`)
- PGGB is a 3-stage pipeline: wfmash → seqwish → smoothxg
- **Parameters can be derived from your own data rather than copied from a paper** — measure identity with `mash triangle` or `wfmash -m` for `-p`, sweep one chromosome for `-s`, and treat `-n` as a check on pggb's automatic count
- In the current PGGB version, `-B` must be explicit; **verify it from `.params.yml`, not from your command line**
- Singularity/Apptainer on HPC. Docker needs the `docker` group rather than root; Podman is rootless
- **impg separates alignment from graph induction** — sweep parameters without redoing wfmash, work one locus at a time, and stop assuming chromosome correspondence (§5.9)

The next chapter evaluates the **quality of the resulting graph**.

## References

- Garrison E, et al. (2024). Building pangenome graphs. *bioRxiv*.
- Li H (2020). The design and construction of reference pangenome graphs with minigraph. *Genome Biol* 21:265.
- Hickey G, et al. (2024). Pangenome graph construction from genome alignments with Minigraph-Cactus. *Nat Biotechnol* 42:663-673.
- Kiryu Y, et al. (2026). AlleleMiner. *DNA Res* 33:dsag004.

### Tools and specifications

- **pggb**: <https://github.com/pangenome/pggb> / docs <https://pggb.readthedocs.io/>
- **wfmash**: <https://github.com/waveygang/wfmash>
- **impg** (implicit pangenome graph): <https://github.com/pangenome/impg>
- **HPRC release 2 build**: <https://github.com/pangenome/HPRCv2>
- **Minigraph-Cactus** procedure: <https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md>
  (a version for this tutorial's data is in the [appendix](appendix_minigraph_cactus.md))

---

[← Chapter 4](04_data_qc.md) | [Chapter 6: Graph QC →](06_graph_qc.md)

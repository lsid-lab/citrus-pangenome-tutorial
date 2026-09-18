# Chapter 4. Data quality control

## Learning objectives

- Understand why QC is essential before pangenome construction
- Learn basic usage of `seqkit stats`
- **Judge QC results by comparison with species expectations rather than fixed thresholds**
- Know the expected values for citrus, and how to translate them to other species

---

## 4.1 Why QC?

Pangenome graphs are built by combining multiple assemblies. If even one assembly is low quality, false structural variants (SVs) will be "discovered" in the graph.

The bioinformatics golden rule applies:
> **Garbage in, garbage out**

QC is mandatory before pangenome construction.

---

## 4.2 Basic assembly terminology (review)

### Contig vs Scaffold
- **Contig**: a continuous sequence built from reads alone
- **Scaffold**: contigs ordered/oriented on chromosomes using auxiliary data (e.g., Hi-C)

Ideally, one chromosome = one contig = one scaffold. Citrus (n=9) should have 9 contigs in a perfect assembly.

### N50
The most important **contiguity indicator**.

**Definition**: sort all sequences from longest to shortest, cumulate lengths, and take the sequence length at which **50% of the total length** is reached.

- Chromosome-scale: N50 in tens of Mb
- Contig-level: N50 in hundreds of kb to a few Mb

### Haplotype-resolved
For diploid organisms, the paternal and maternal chromosomes are **assembled separately**. As explained in Chapter 2, the Isobe 2023 assemblies are of this kind.

---

## 4.3 Basic statistics with `seqkit stats`

`seqkit stats` provides basic FASTA statistics in a single command.

### Installation

Most tools in this tutorial install via **conda / mamba**. If you do not have an environment yet, **Miniforge** is the current recommendation.

**Why Miniforge**: the conda-forge channel is preconfigured and **`mamba` (fast dependency solving) ships with it out of the box**. The Anaconda/Miniconda `defaults` channel carries commercial-use licensing conditions that depend on your organization's size; Miniforge only looks at conda-forge, so that concern does not arise. Note that **Mambaforge was retired in January 2025** and folded into Miniforge3 (since Miniforge 23.3.1 the two are essentially identical).

```bash
# Same command on Linux and macOS - OS and architecture are resolved automatically
wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash "Miniforge3-$(uname)-$(uname -m).sh"

# Reopen your shell, then check
mamba --version
```

> `$(uname)` expands to `Linux` or `Darwin`, and `$(uname -m)` to `x86_64` / `aarch64` / `arm64`. Without `wget`, use
> `curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"`.

A dedicated environment makes it easier to add tools in later chapters:

```bash
mamba create -n citrus-pg -c conda-forge -c bioconda seqkit
mamba activate citrus-pg
```

**seqkit is also available as a single binary** if you would rather not use conda:

```bash
wget https://github.com/shenwei356/seqkit/releases/download/v2.10.0/seqkit_linux_amd64.tar.gz
tar xzf seqkit_linux_amd64.tar.gz && mv seqkit ~/bin/
```

> The PGGB toolchain used from Chapter 5 onward is installed as a **Singularity image**, not via conda (§5.5, `scripts/setup_singularity_wrappers.sh`). Conda here covers only this chapter's QC tools.

### Running it

```bash
seqkit stats -a data/CUN/CUNphKi_r1.0.pmol.fasta.gz
```

The `-a` (all) flag adds N50 and other details.

### Example output for Satsuma hap1

```
file                              format  type  num_seqs  sum_len      min_len    avg_len    max_len    N50
CUNphKi_r1.0.pmol.fasta.gz        FASTA   DNA          9  348,509,231  15,234,891 38,723,247 48,123,456  38,241,932
```

### Column meanings and expected values (haploid citrus)

| Column | Meaning | Expected value |
|---|---|---|
| `num_seqs` | Number of sequences (contigs/chromosomes) | **9** (n=9 chromosomes) |
| `sum_len` | Total base count | **300-360 Mb** |
| `min_len` | Shortest sequence | ≥ 10 Mb (smallest chromosome) |
| `avg_len` | Average length | ~35 Mb |
| `max_len` | Longest sequence | ≤ 50 Mb (largest chromosome) |
| `N50` | Contiguity indicator | Larger is better, 20+ Mb |
| `GC(%)` | GC content | 34-38% |

---

## 4.4 Batch QC across six haploids

Use `scripts/qc01_stats.sh` to compute statistics for all six haploids at once:

```bash
bash scripts/qc01_stats.sh tables/samplesheet.tsv .
```

Output: `qc/stats_summary.tsv`

Example results:

| sample | hap | num_seqs | sum_len_Mb | N50_Mb | GC% | verdict | note |
|---|---|---|---|---|---|---|---|
| CUN | 1 | 9 | 348.5 | 34.5 | 35.9 | PASS | 25-55 Mb larger than the parents |
| CUN | 2 | 9 | 357.6 | 44.2 | 35.0 | PASS | same |
| CKI | 1 | 9 | 304.2 | 33.5 | 36.0 | PASS | |
| CKI | 2 | 9 | 310.3 | 32.5 | 36.0 | PASS | |
| CKU | 1 | 9 | 323.9 | 36.6 | 35.9 | PASS | |
| CKU | 2 | 9 | 303.4 | 32.9 | 36.0 | PASS | |

**All six haploids pass.** Every one of them falls inside the citrus expectations `qc01_stats.sh` uses (290-370 Mb, 9 chromosomes, GC 34-38%).

CUN's two haplotypes being larger than the parents is **within the published size range for Satsuma** (Shimizu et al. 2017: 359.7 Mb; Kawahara et al. 2020: 346 Mb). The two parents are different species, and a 10-20% genome size difference between close relatives is unremarkable.

It does matter downstream, though: **the longer a path, the lower its Jaccard similarity**, so this difference shows up in the Chapter 6 interpretation (§6.4.4). Keep the observation in mind.

---

## 4.5 Points to check — compare against expectations, not fixed thresholds

Every number given so far is **an expectation for citrus specifically** (n=9, haploid ~300-360 Mb). Memorize them as absolute thresholds and they stop working the moment you move to another species.

Apply "`num_seqs ≤ 20` means chromosome-scale" to human, for example, and **even T2T-CHM13 (24 sequences) is flagged as semi-fragmented**. Replace the fixed thresholds with the question: **how does this compare to what is expected for this species?**

> The six assemblies used in this tutorial are all **published, high-quality, chromosome-scale assemblies** (9 pseudomolecules, N50 over 30 Mb). There is in fact **no need** to re-do QC here and find problems.
>
> The point of this section is to give you a method you can take to **your own data**: what to check, and what to check it against.

### Point 1: Contiguity — has it reached chromosome scale?

What matters is not the absolute sequence count but **how it compares to the expected chromosome number n for that species**.

| Metric | How to judge |
|---|---|
| `num_seqs` | Is it close to **n** (plus organelles and unplaced scaffolds)? |
| `N50` | Is it close to **the median chromosome length of that species**? At chromosome scale, N50 ≈ median chromosome length |
| `L50` | Is it about **half the chromosome count** (n/2)? Far above n means fragmentation |

| Species | n (haploid) | Expected num_seqs | Expected N50 |
|---|---|---|---|
| Citrus | 9 | ~9 | 30-40 Mb |
| Human (T2T-CHM13) | 23 | ~24 | ~150 Mb |
| Rice | 12 | ~12 | ~30 Mb |
| *Arabidopsis* | 5 | ~5 | ~23 Mb |

**Never judge on `num_seqs` alone.** An assembly can carry thousands of unplaced scaffolds and still be chromosome-scale for practical purposes, as long as **over 95% of the total bases sit in the top n sequences**. Conversely, few sequences may just mean long scaffolds full of `N`. Always read `num_seqs` together with `N50`.

Two assemblies of the very same Satsuma mandarin illustrate the gap:

| Assembly | Total | Sequences | N50 |
|---|---:|---:|---:|
| Shimizu et al. (2017) draft | 359.7 Mb | 20,876 scaffolds | 386 kb |
| Isobe et al. (2023) `CUNphKi` | 348.5 Mb | 9 pseudomolecules | ~38 Mb |

**Nearly the same genome size; contiguity differs by a factor of 100.** Only the latter is usable for pangenome construction.

### Point 2: Total length — compare with the known genome size

There is no universal "normal range" for `sum_len`. Compare against **a genome size estimate for your species**, from:

- Assembly sizes in prior publications on the same or a close species
- Flow cytometry measurements (e.g. the [Plant DNA C-values Database](https://cvalues.science.kew.org/))
- A k-mer based estimate from your own reads (GenomeScope2 and similar)

|Relative to the estimate|Possible interpretation|
|---|---|
|Within ±10%|As expected|
|**Smaller**|Missing sequence, or collapsed repeats|
|**Larger**|Duplicated haplotype (leakage), contamination — or an outdated estimate|

**Note that neither direction is evidence of a problem on its own.** Genome size routinely differs by 10-20% even between close relatives, and the estimates themselves have spread. CUN in §4.4 is exactly that case: **larger than the two parents, yet inside the published range for Satsuma.**

### Point 3: hap1 vs hap2 — this one is species-independent

The two haplotypes of one individual are **the same chromosome set of the same individual**, so there is no good reason for their sizes to diverge much. This is one of the few checks that **does not depend on the species**.

- Difference ≤ 5%: consistent
- Difference 5-15%: worth checking
- Difference > 15%: suspect phasing failure, or sequence piling into one haplotype

(Highly heterozygous individuals, large hemizygous regions, and sex chromosomes are exceptions.)

### Point 4: GC content — compare with the known value for the species

GC content varies widely by species: citrus ~35%, human ~41%, rice ~44%, and *Plasmodium falciparum* around 19%. Judge not by an absolute band but by **whether you are within 1-2 pp of published assemblies of the same species**.

Looking at **per-sequence GC** is more informative than the global figure:

```bash
# name / length / GC% for each sequence
seqkit fx2tab -nlg data/CUN/CUNphKi_r1.0.pmol.fasta.gz
```

**A single sequence far off the rest** points to contamination — bacteria, organelles, or another species. (Organellar genomes differ in GC from the nuclear genome, so this catches them.)

### Running `qc01_stats.sh` on another species

The expected values in `scripts/qc01_stats.sh` can be overridden from the environment. For non-citrus data, pass that species' values:

```bash
# Example: human (n=23, ~3.1 Gb, GC ~41%)
EXPECTED_NUM_CHR=23 \
EXPECTED_SIZE_MIN=2900 EXPECTED_SIZE_MAX=3300 \
EXPECTED_GC_MIN=40.0 EXPECTED_GC_MAX=42.0 \
bash scripts/qc01_stats.sh tables/samplesheet.tsv .
```

The contiguity thresholds are derived automatically from `EXPECTED_NUM_CHR`.

---

## 4.6 Additional QC tools (optional)

### BUSCO / compleasm (completeness)

Evaluate **assembly completeness** by counting how many single-copy orthologs from a related lineage are detected. How to read the metrics and how to pick a lineage dataset are covered in §1.3.

**Install** (using the mamba from §4.3):

```bash
# BUSCO itself, or compleasm, its faster reimplementation
mamba install -c conda-forge -c bioconda busco
mamba install -c conda-forge -c bioconda compleasm
```

**Fetch the lineage dataset** (once; takes a few minutes):

```bash
compleasm download eudicots_odb10
```

**Run**:

```bash
# compleasm (faster than BUSCO)
compleasm run -a data/CUN/CUNphKi_r1.0.pmol.fasta.gz \
              -o busco/CUN_hap1 \
              -l eudicots_odb10 \
              -t 16
```

Expected: **Complete > 95%, Duplicated < 5%**. These are haplotype-resolved assemblies, so run it **separately on each of the six haploids**.

### Merqury (k-mer QV and haplotype leakage)

Compares an assembly to a k-mer database built from HiFi raw reads:
- **QV (base accuracy)**: expected > 40
- **Completeness**: fraction of read-derived k-mers found in the assembly
- **Duplication**: detects haplotype leakage

**Install**: Merqury depends on Java and R, so a **dedicated environment** is safest. `meryl` comes along as a dependency.

```bash
mamba create -n merqury -c conda-forge -c bioconda merqury openjdk=11
mamba activate merqury

merqury.sh --version   # check it runs
meryl --version
```

> If you hit an error like `$MERQURY is not set`, export it first: `export MERQURY="$CONDA_PREFIX/share/merqury"`.

**Run**:

```bash
# 1) Build a k-mer database from the raw HiFi reads
meryl count k=21 output hifi.meryl hifi_reads.fastq.gz

# 2) Compare against the assemblies (passing both haps also compares them to each other)
merqury.sh hifi.meryl CUNphKi_r1.0.fasta.gz CUNphKu_r1.0.fasta.gz CUN
```

With trio reads (both parents plus the offspring) you can additionally build `hapmers` and obtain the **false duplication rate** and **hap-mer blob plot** — the proper tools for deciding whether haplotype leakage occurred.

**We do not run Merqury in this tutorial**, since we do not handle raw reads. The reads are available from DDBJ DRA BioProject `PRJDB15866` (§3.1) if you want to try it yourself.

---

## Chapter summary

- Use `seqkit stats` for basic statistics
- **Judge against what is expected for your species, not against fixed thresholds.** Citrus's 300-360 Mb / 9 chromosomes / GC 34-38% are citrus numbers
- Never judge on `num_seqs` alone — read it with `N50` (the same Satsuma exists both as a 20,876-scaffold assembly and as a 9-pseudomolecule one)
- Comparing hap1 against hap2 is one of the few **species-independent** checks
- **CUN's two haps are larger than the rest, but well within the published size range for Satsuma (346-360 Mb).** "Larger than the other samples I happen to have" is not "anomalous"
- That size difference nevertheless **affects the Jaccard similarities in Chapter 6** (§6.4.4), so carry the observation forward

The next chapter starts building the **pangenome graph**.

## References

- Shen W, et al. (2024). SeqKit2. *iMeta* 3:e191.
- Rhie A, et al. (2020). Merqury: reference-free quality, completeness, and phasing assessment for genome assemblies. *Genome Biol* 21:245.
- Shimizu T, et al. (2017). Draft sequencing of the heterozygous diploid genome of satsuma (*Citrus unshiu* Marc.) using a hybrid assembly approach. *Front Genet* 8:180. (359.7 Mb total, 20,876 scaffolds, N50 386 kb)
- Kawahara Y, et al. (2020). Mikan Genome Database (MiGD). *Breed Sci* 70(2). (Satsuma, 346 Mb)
- Manni M, et al. (2021). BUSCO update. *Mol Biol Evol* 38:4647-4654.
- Huang N, Li H (2023). compleasm: a faster and more accurate reimplementation of BUSCO. *Bioinformatics* 39:btad595.
- **Miniforge**: <https://github.com/conda-forge/miniforge>
- **Plant DNA C-values Database**: <https://cvalues.science.kew.org/>

---

[← Chapter 3](03_data_acquisition.md) | [Chapter 5: Pangenome graph construction →](05_graph_construction.md)

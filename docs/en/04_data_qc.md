# Chapter 4. Data quality control

## Learning objectives

- Understand why QC is essential before pangenome construction
- Learn basic usage of `seqkit stats`
- Know the expected values for citrus assemblies
- Recognize abnormal patterns and how to handle them

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

```bash
# conda
mamba install -c bioconda seqkit

# or single binary
wget https://github.com/shenwei356/seqkit/releases/download/v2.10.0/seqkit_linux_amd64.tar.gz
tar xzf seqkit_linux_amd64.tar.gz && mv seqkit ~/bin/
```

### Running it

```bash
seqkit stats -a data/satsuma/CUNphKi_r1.0.pmol.fasta.gz
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

| sample | hap | num_seqs | sum_len_Mb | N50_Mb | GC% | verdict |
|---|---|---|---|---|---|---|
| STS | 1 | 9 | 348.5 | 34.5 | 35.9 | WARN (size slightly large) |
| STS | 2 | 9 | 358.0 | 44.2 | 35.0 | WARN (size slightly large) |
| KSH | 1 | 9 | 304.2 | 33.5 | 36.0 | OK |
| KSH | 2 | 9 | 310.3 | 32.5 | 36.0 | OK |
| KNN | 1 | 9 | 323.9 | 36.6 | 35.9 | OK |
| KNN | 2 | 9 | 303.4 | 32.9 | 36.0 | OK |

---

## 4.5 Points to check

### Point 1: Chromosome-scale determination

- `num_seqs ≤ 20`: chromosome-scale (OK)
- `num_seqs 20-200`: semi-fragmented (scaffolds not fully joined) (WARN)
- `num_seqs > 200`: over-fragmented, difficult for pangenome (FAIL)

### Point 2: Size

- 300-360 Mb (haploid): normal range
- < 290 Mb: possible missing information
- > 370 Mb: incomplete collapse, potential haplotype mixing

### Point 3: hap1 vs hap2 within a cultivar

The two haplotypes of the same cultivar should have similar sizes:

- Difference ≤ 5%: perfect match (OK)
- Difference 5-15%: acceptable but worth checking (WARN)
- Difference > 15%: possible **phasing failure** (FAIL)

### Point 4: GC content

- 34-38%: normal citrus range (OK)
- < 34% or > 38%: possible contamination

---

## 4.6 A real-data observation: haplotype leakage

Comparing the six haploids to expected values, **Satsuma's two haplotypes are noticeably larger**:

- STS hap1: 348.5 Mb
- STS hap2: 358.0 Mb
- Other four haps: 303-324 Mb

A 15-20 Mb (5-7%) difference is within tolerance, but **a Merqury k-mer analysis** (see below) reveals a more subtle issue.

### Discovery through Merqury

Merqury (Rhie et al., 2020) evaluates **assembly redundancy** from k-mer frequencies. A collaborator's Merqury analysis yielded:

| Assembly | Unique k-mer% | 2× occurrence % | 3+× % | Total bp |
|---|---|---|---|---|
| CKUhap2 (KNN hap2) | 90.17% | 5.46% | 4.37% | 303.4 Mb |
| CKIhap1 (KSH hap1) | 89.95% | 5.61% | 4.44% | 304.2 Mb |
| CKUhap1 (KNN hap1) | 89.42% | 5.77% | 4.81% | 323.9 Mb |
| CKIhap2 (KSH hap2) | 89.28% | 6.06% | 4.66% | 310.3 Mb |
| **CUNphKu (STS hap2)** | **82.64%** | **11.89%** | **5.47%** | **358.0 Mb** |
| **CUNphKi (STS hap1)** | **79.64%** | **14.65%** | **5.71%** | **348.4 Mb** |

**Satsuma's two haps are clearly anomalous**:
- Unique k-mer rate is ~10 pp lower
- 2× occurrence k-mers are **twice as high** as parental lines
- Total size is 30 Mb larger

This is called **haplotype leakage**: the trio phasing could not fully separate the two haplotypes, leaving **shared sequence in both haps**.

Specifically:
- Repeat-rich regions (near centromeres, TE clusters, etc.)
- Regions with unusually low/high heterozygosity

are hard to phase and remain collapsed in both haplotypes.

### Why we address this in the tutorial

**Reason 1: Real data is never perfect**
Isobe et al. 2023 is a world-class trio phasing example, yet not flawless. Experiencing the gap between textbook and real data is valuable.

**Reason 2: Necessary for interpreting later QC**
Haplotype leakage will affect our pangenome graph QC in Chapter 6. Knowing about it here allows correct interpretation.

**Reason 3: Honest declaration of limits**
The graph we build will contain some effects of this leakage—we should acknowledge this openly.

---

## 4.7 Additional QC tools (optional)

### BUSCO / compleasm (completeness)

Evaluate **assembly completeness** by counting how many single-copy orthologs from a related lineage are detected.

```bash
# compleasm (faster than BUSCO)
compleasm run -a data/satsuma/CUNphKi_r1.0.pmol.fasta.gz \
              -o busco/STS_hap1 \
              -l eudicots_odb10 \
              -t 16
```

Expected: **Complete > 95%, Duplicated < 5%**.

### Merqury (k-mer QV and haplotype leakage)

Compares an assembly to a k-mer database built from HiFi raw reads:
- **QV (base accuracy)**: expected > 40
- **Completeness**: fraction of read-derived k-mers found in the assembly
- **Duplication**: detects haplotype leakage

```bash
meryl count k=21 output hifi.meryl hifi_reads.fastq.gz
merqury.sh hifi.meryl CUNphKi_r1.0.fasta.gz CUNphKu_r1.0.fasta.gz STS
```

We do not run Merqury in this tutorial (we do not use raw HiFi reads), but **the interpretation is based on the colleague's data cited in §4.6**.

---

## Chapter summary

- Use `seqkit stats` for basic statistics
- Citrus haploids should be 300-360 Mb, 9 chromosomes, GC 34-38%
- **Satsuma's two haps are larger** than parental lines, suggesting **haplotype leakage**
- Real data is not perfect; remember this observation for downstream QC

The next chapter starts building the **pangenome graph**.

## References

- Shen W, et al. (2024). SeqKit2. *iMeta* 3:e191.
- Rhie A, et al. (2020). Merqury. *Genome Biol* 21:245.
- Manni M, et al. (2021). BUSCO update. *Mol Biol Evol* 38:4647-4654.
- Huang N, Li H (2023). compleasm. *Bioinformatics* 39:btad595.

---

[← Chapter 3](03_data_acquisition.md) | [Chapter 5: Pangenome graph construction →](05_graph_construction.md)

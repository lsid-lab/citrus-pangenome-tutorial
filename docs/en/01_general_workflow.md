# Chapter 1. General workflow

## Learning objectives

- Understand the general workflow for pangenome construction
- Learn common terminology used across tools: PanSN, GFA, node/edge/path
- Identify the decision points at each step

---

## 1.1 Overview

Whatever pangenome tool you choose (PGGB, Minigraph-Cactus, Minigraph), the overall workflow is similar:

```
┌─────────────────────────────────┐
│  1. Data collection             │  ← Fetch assemblies from papers/DBs
│  (per-individual assemblies)    │
└──────────┬──────────────────────┘
           │
┌─────────────────────────────────┐
│  2. Data QC                     │  ← Verify assembly quality
│  (contiguity, completeness,     │
│   contamination check)          │
└──────────┬──────────────────────┘
           │
┌─────────────────────────────────┐
│  3. Input preparation           │  ← Unify naming, bgzip
│  (PanSN naming, chr grouping)   │
└──────────┬──────────────────────┘
           │
┌─────────────────────────────────┐
│  4. Graph construction          │  ← PGGB / MC / Minigraph
│  (align, induction, smoothing)  │
└──────────┬──────────────────────┘
           │
┌─────────────────────────────────┐
│  5. Graph QC                    │  ← Structure, preservation, biology
│  (odgi stats / similarity)      │
└──────────┬──────────────────────┘
           │
┌─────────────────────────────────┐
│  6. Analysis and visualization  │  ← SV extraction, VCF, viz
│  (vg deconstruct, odgi viz)     │
└─────────────────────────────────┘
```

**This tutorial walks through all six steps hands-on.**

---

## 1.2 Common concepts

### 1.2.1 PanSN-spec (Pangenome Sequence Naming Specification)

Most pangenome tools require a **special sequence naming convention**.

**PanSN format**: `sample#haplotype#contig`

Example:
```
satsuma#1#chr01    # Satsuma mandarin, hap1, chromosome 1
satsuma#2#chr01    # Satsuma mandarin, hap2, chromosome 1
kishu#1#chr01      # Kishu mandarin, hap1, chromosome 1
```

The `#` character is the delimiter. This scheme lets tools mechanically identify **which cultivar/haplotype a sequence belongs to**.

Most assemblies use different naming (e.g., `CUNphKi_r1.0ch1`), so a **rename step** is needed before pangenome construction. This is handled by `pg01_prepare_input.sh` in this tutorial.

### 1.2.2 GFA format (Graphical Fragment Assembly)

The **standard file format** for pangenome graphs. It contains three main record types.

**S (Segment)**: defines a node (a contiguous sequence)
```
S  1  ACGTACGT     # node ID = 1, sequence = ACGTACGT
S  2  GCATT
S  3  TTAA
```

**L (Link)**: defines an edge between nodes
```
L  1  +  2  +  0M    # node 1 right end → node 2 left end
L  1  +  3  +  0M    # node 1 right end → node 3 left end
```

**P (Path)**: describes how a specific sequence traverses the graph
```
P  satsuma#1#chr01  1+,2+,4+,5+  ...
P  kishu#1#chr01    1+,3+,4+,5+  ...
```

Shared nodes, branching edges, and traversal paths together represent the pangenome.

### 1.2.3 Node / Edge / Path relationships

```
Path A: 1 → 2 → 4 → 5
Path B: 1 → 3 → 4 → 5

  node2      
  ↗   ↘     
node1  node4 → node5
  ↘   ↗     
  node3
```

- **Node** (S): the minimum unit holding DNA sequence
- **Edge** (L): a connection between nodes (allowed transitions)
- **Path** (P): a **traversal** representing one individual's chromosome

**Key properties**:
- Identical sequences are **collapsed** into a single node
- Differences appear as **branches**
- Each original haplotype is **fully recoverable** by following its path

This is what "reference-free" means: no single individual defines the backbone.

### 1.2.4 Compression ratio

The **most important metric** for graph quality:

```
compression ratio = graph total bp / input total bp
```

Example:
- Input: 6 haploids × 30 Mb = 180 Mb
- Graph total bp: 45 Mb
- Compression ratio = 45 / 180 = **0.25** (25%)

Meaning: the information of 6 haploids is compressed into 25% of the input size = 75% of the sequence is **shared**.

|Compression|Interpretation|
|---|---|
|0.15-0.30 | Good (closely related cultivars) |
|0.30-0.50 | Fair (more divergent within-species) |
|0.50-0.80 | Under-compressed (reconsider parameters) |
|> 0.90 | Essentially uncompressed (serious tool/data issue) |

We use this metric to judge results in Chapter 6.

---

## 1.3 Decision points at each step

### Step 2 (Data QC)
- **Contiguity** (is N50 at chromosome scale?)
- **Completeness** (BUSCO complete > 95%?)
- **Contamination** (any foreign sequences mixed in?)
- **Haplotype separation** (are hap1 and hap2 similar in size?)

### Step 3 (Input preparation)
- **Rename** all haploid files to PanSN format
- Decide **chromosome-by-chromosome** vs whole-genome
  - Small genomes (< 300 Mb): whole genome may be fine
  - Medium-large genomes: per-chromosome is required (memory/compute constraints)

### Step 4 (Graph construction)
- Choose the tool (discussed in Chapter 5)
- **Key parameters to consider** (as tool-general concepts):
  - **Number of haplotypes**: must exactly match the input count
  - **Percent identity**: threshold for allowed base identity
  - **Segment length**: minimum unit for alignment
  - Other tool-specific parameters (PGGB details in Chapter 5)

### Step 5 (Graph QC)
- Is the compression ratio in range?
- Are all input paths preserved (no loss)?
- Is the biological signal (pedigree, phylogeny) consistent?
- Are SV/SNP counts realistic?

### Step 6 (Analysis)
- Generate VCF, count SVs
- Visualize (odgi viz, odgi draw)
- Analyze path similarity
- Biological interpretation

---

## 1.4 Why chromosome-by-chromosome?

Large-genome pangenome construction is typically done **per chromosome**.

Reasons:
1. **Memory constraints**: whole-genome runs may need hundreds of GB of RAM
2. **Parallelization**: each chromosome can run as an independent job on HPC
3. **Fault tolerance**: if one chromosome fails, the others are unaffected
4. **Ease of chromosome-level comparison**: useful when zooming into one chromosome

Preconditions:
- Each assembly is **chromosome-scale**
- Chromosome numbering **corresponds across assemblies** (i.e., "chr01" refers to the same chromosome everywhere)

*Citrus* (n=9) means 9 chromosomes × 6 haploids = **54 sequences** total in the graph.

---

## Chapter summary

- Pangenome construction follows a general 6-step workflow
- Common concepts across tools: **PanSN**, **GFA**, **node/edge/path**, **compression ratio**
- Chromosome-by-chromosome construction is standard practice for large genomes

The next chapter introduces the three target cultivars (Satsuma, Kishu, Kunenbo) used in this tutorial.

## References

- Garrison E, et al. (2018). Variation graph toolkit improves read mapping by representing genetic variation in the reference. *Nat Biotechnol* 36:875-879.
- Li H (2020). The design and construction of reference pangenome graphs with minigraph. *Genome Biol* 21:265.
- Heumos S, et al. (2024). Pangenome graph layout by Path-Guided SGD. *Bioinformatics*.

---

[← Chapter 0](00_what_is_pangenome.md) | [Chapter 2: Target organisms →](02_target_organisms.md)

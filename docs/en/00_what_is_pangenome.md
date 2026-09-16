# Chapter 0. What is a Pangenome?

## Learning objectives

- Understand why a single reference genome is not sufficient
- Learn the definition of pangenome and its three representations
- Understand the advantages of graph-based pangenomes
- Appreciate why pangenomes have become practical only recently

---

## 0.1 Limitations of a single reference genome

Traditional genetic analysis relies on a single **reference genome** representing one individual: GRCh38 for humans, TAIR10 for *Arabidopsis*, Clementine v1.0 for *Citrus*, and so on.

The standard workflow is:
1. Map reads from a new individual to the reference
2. Detect differences (SNPs and small indels) from the reference
3. Describe the individual by that difference set

This works well for many analyses, but has **serious limitations**.

### Problem 1: Sequences absent from the reference are invisible

If individual A carries a **5 kb insertion** that does not exist in the reference, most tools cannot map any reads to it, and the insertion is effectively **ignored**. If that 5 kb happens to contain an important gene (disease resistance, fruit quality, etc.), it is silently missed.

### Problem 2: The reference represents only one individual

GRCh38 was assembled from a small number of individuals. Among the 8 billion people on Earth, not a single one matches GRCh38 100%. The same applies to crops: Clementine v1.0 differs from Satsuma mandarin at millions of loci, and treating each as a "variant" is inefficient and misleading.

### Problem 3: Population diversity cannot be represented by one individual

A reference genome is neither an "average" nor a "typical" specimen—it happens to be **the first individual sequenced**. Breeding and population genetics need **diversity across the population**, not just differences from an arbitrary standard.

---

## 0.2 Definition of a Pangenome

The pangenome concept was developed to solve these problems.

> **Pangenome** = a data structure representing **all genetic material present in a population**

The "population" can be within a species (*intraspecific*) or across a genus (*interspecific*, called a "super-pangenome"). For *Citrus*, cultivar groups such as "Satsuma + Kishu + Kunenbo" form one such population.

### Three representations of a pangenome

Three approaches have emerged historically:

**1. Gene-based (linear pangenome)**

Represent population diversity by the **presence/absence** of genes across individuals.
- Genes present in all individuals: **core genes**
- Genes present in some: **accessory genes**
- Genes unique to one individual: **private genes**

The bacterial pangenome (Tettelin et al., 2005) pioneered this approach. In plants, Golicz et al. 2016 (cabbage) was an early example.

- **Strengths**: simple, biologically interpretable
- **Limitations**: cannot represent non-gene elements (regulatory regions, non-coding RNA, structural variants)

**2. Sequence-based**

Concatenate all assemblies, compress redundancy, and retain the complete set of sequences.

- **Strengths**: no information loss
- **Limitations**: large size, structural relationships not explicit

**3. Graph-based** ← **the approach used in this tutorial**

Represent multiple genomes as a **graph**: shared sequences are unified as nodes, differences appear as branches.

```
Linear reference:
A T G G C C A T G C A T G

Graph pangenome:
     ┌─T─┐
A T G┤   ├C A T G
     └─A─┘  ↑
            Variant at this position: T or A
```

- **Strengths**: unified representation of all variants (SNP, indel, SV, PAV); reference-free analysis possible
- **Limitations**: implementation complexity; requires familiarity with visualization and manipulation

---

## 0.3 Why now?

The pangenome concept has existed since 2005, but practical construction became feasible only recently. Three technical advances enabled this:

### 1. Long-read sequencing at scale

**PacBio HiFi** (from 2019) and **Oxford Nanopore ULK** (from 2020) provide 15–25 kb high-accuracy reads. This enables:
- Correct assembly of repeat regions
- Detection of large SVs in single reads
- Base-level quality QV > 40 (99.99% accurate)

### 2. Haplotype-resolved assembly

**HiFi + Hi-C** and **trio phasing** now allow the paternal and maternal chromosomes of a diploid to be **assembled separately**. This is essential for:
- Correctly capturing heterozygous variation
- Increasing the resolution of population analyses

**This tutorial uses exactly this kind of haplotype-resolved trio-phased data** (Isobe et al., 2023).

### 3. Mature pangenome construction tools

- **Minigraph** (Li et al., 2020): SV-focused
- **Minigraph-Cactus** (Hickey et al., 2024): SV + small variants, reference-guided
- **PGGB** (Garrison et al., 2024): reference-free, comprehensive ← used here

These are now practical enough to build pangenomes on standard HPC hardware in hours to days.

---

## 0.4 Applications of pangenomes

Representative pangenome projects:

- **Human Pangenome Reference Consortium** (Liao et al., 2023): 47 individuals, 94 haploid equivalents
- **Cattle super-pangenome** (Leonard et al., 2023): cross-breed diversity
- **Sorghum pangenome** (Feng et al., 2025): genetic resource mining for breeding
- **Citrus pangenome** (Huang et al., 2023; Ollitrault et al., 2025): genus-level super-pangenome

These contribute to:
- **Trait-linked gene discovery**: SV-based loci invisible to SNP-based GWAS
- **Systematic evaluation of breeding materials**: visualizing untapped genetic resources
- **Population genetics**: reconstructing histories of admixture and selection

The Citrus pangenome in this tutorial is small in scale but provides an excellent entry point to the same technology stack.

---

## 0.5 What you will gain from this tutorial

By the end of this tutorial, you will be able to:

1. **Explain concepts**: pangenome, PanSN, GFA, node/edge/path in your own words
2. **Run construction**: execute PGGB on an HPC and produce a graph
3. **Judge quality**: evaluate the graph from numbers and visualization
4. **Interpret results**: extract VCF from the graph, cross-check with pedigree, and discuss
5. **Recognize pitfalls**: experience real-data quirks such as haplotype leakage and version-specific tool behavior

---

## Chapter summary

- A single reference cannot capture population diversity
- Pangenomes solve this by **representing all genetic material** in an integrated form
- Graph-based pangenomes are the current mainstream; they enable reference-free analysis
- Long-read sequencing, haplotype-resolved assembly, and mature tools have made this practical

The next chapter walks through the **general workflow** for building a pangenome, using terminology common to all major tools.

## References

- Tettelin H, et al. (2005). Genome analysis of multiple pathogenic isolates of *Streptococcus agalactiae*: implications for the microbial "pan-genome". *PNAS* 102:13950-13955.
- Golicz AA, et al. (2016). The pangenome of an agronomically important crop plant *Brassica oleracea*. *Nat Commun* 7:13390.
- Liao WW, et al. (2023). A draft human pangenome reference. *Nature* 617:312-324.
- Garrison E, et al. (2024). Building pangenome graphs. *bioRxiv*.
- Isobe S, et al. (2023). Haploid-resolved and chromosome-scale genome assembly in *Citrus unshiu* and its parental species. *bioRxiv* 2023.06.02.543356.

---

[← README](../../README.md) | [Chapter 1: General workflow →](01_general_workflow.md)

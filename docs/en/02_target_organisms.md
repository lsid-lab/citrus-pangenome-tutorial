# Chapter 2. Target organisms: Satsuma, Kishu, and Kunenbo mandarins

## Learning objectives

- Understand the biology and history of the three target cultivars
- Learn that Satsuma mandarin = Kishu × Kunenbo (F1)
- Appreciate the significance of trio phasing and its relevance here
- Understand why these three cultivars are an ideal pangenome tutorial dataset

---

> **On the cultivar abbreviations used here**
>
> Throughout this tutorial the three cultivars are referred to by the three-letter codes
> **CUN** / **CKI** / **CKU**. These are not labels invented for the tutorial: they are
> **the prefixes of the assembly IDs published by Isobe et al. (2023) and Plant GARDEN**
> (`CUNphKi_r1.0`, `CKIhap1_r1.0`, `CKUhap1_r1.0` — see §2.5 for how to read them).
>
> The same codes are used consistently everywhere: downloaded file names, the `sample_id`
> column of `samplesheet.tsv`, pangenome graph path names (`CUN#1#chr01`), and VCF sample
> columns. You can line the tutorial up against the original paper and Plant GARDEN directly.

## 2.1 Satsuma mandarin (*Citrus unshiu*)

**Satsuma mandarin** (abbreviated **CUN**: ***C**itrus **un**shiu*), Japan's iconic winter citrus, was discovered in Nagashima, Kagoshima Prefecture during the early Edo period (17th century).

- **Chromosome number**: 2n = 18 (n = 9)
- **Genome size**: ~320-360 Mb (haploid)
- **Reproduction**: seedless (parthenocarpic), with nucellar polyembryony (seeds contain maternal clones)
- **Cultivars**: Miyagawa wase, Okitsu wase, Aoshima, Nankan No. 20, and many more

Seedlessness is a benefit for cultivation but creates **challenges for genetic analysis**: standard crossing experiments (parent A × parent B → observe F1) are not feasible.

The origin of Satsuma was long uncertain, but **recent genomic analyses** (Fujii et al., 2016; Wu et al., 2018) confirmed that it is an **F1 hybrid of Kishu mandarin (CKI) × Kunenbo (CKU)**.

### Genetic composition

According to Wu et al. 2018, Satsuma is composed of:
- **Mandarin (*C. reticulata*)** gene pool: over 80%
- Small **Pummelo (*C. maxima*)** gene pool contribution

This reflects admixture inherited from both parents.

---

## 2.2 Kishu mandarin (*Citrus kinokuni*)

**Kishu mandarin** (**CKI**: ***C**itrus **ki**nokuni*) is the **maternal parent** of Satsuma.

- **Origin**: originally from China; present in Japan by the Heian period (9th century)
- **Name**: derived from Kishu (present-day Wakayama Prefecture), where cultivation flourished
- **Traits**: seeded, small fruit, thin skin
- **Historical significance**: during the Edo period, Kishu mandarin was *the* citrus of Japan
- **Chromosomes**: 2n = 18

Kishu is barely grown commercially today and is sometimes considered an endangered cultivar. Its genetic importance as Satsuma's mother, however, drives conservation and research interest.

### Genetic characteristics

- **Nearly 100% mandarin ancestry** (minimal pummelo admixture)
- Considered close to a "pure mandarin"
- Contributes the "mandarin side" of Satsuma

---

## 2.3 Kunenbo (*Citrus nobilis* var. *kunip*)

**Kunenbo** (**CKU**: ***C**itrus **ku**nip*) is the **paternal parent** of Satsuma.

- **Origin**: from Indochina; introduced to Japan around 1350 (Muromachi period)
- **Etymology**: multiple theories (one is "nine years to flower")
- **Traits**: seeded, relatively large fruit, distinctive aroma
- **Chromosomes**: 2n = 18

Almost no commercial cultivation remains in Japan, though preservation lines exist in Kagoshima and elsewhere.

### Genetic characteristics

- **Mandarin-based, with 10-15% pummelo (*C. maxima*) admixture**
- This pummelo ancestry is partially inherited by Satsuma

---

## 2.4 Family diagram

The three cultivars form a small family:

```
        Kishu (CKI)                  Kunenbo (CKU)
   [Mandarin ~100%]           [Mandarin ~85%]
   [ hap1 | hap2 ]            [ hap1 | hap2 ]
         [Pummelo 15%]
              \                       /
               \                     /
                \                   /
                 \                 /
              (crossing, F1 formation)
                    ┌───────┐
                    ↓       ↓
              Satsuma (CUN)
              [Mandarin ~90%, Pummelo ~10%]
  [CUN#1 = CKU-derived | CUN#2 = CKI-derived]
```

### A haplotype is not a copy of a parental chromosome

It is true that each of CUN's two haplotypes **derives from one parent**: **CUN#1 from the father (CKU, Kunenbo) and CUN#2 from the mother (CKI, Kishu)**. Trio phasing (§2.5) guarantees this assignment.

But **neither haplotype is one of the parent's two chromosomes passed down intact.** During the meiosis that produces a gamete, **crossovers** occur between the parent's two homologous chromosomes, and it is the resulting **recombinant chromosome** that is transmitted. So:

- **CUN#1** = a **mosaic** of CKU_hap1 and CKU_hap2 (recombined in the father's meiosis)
- **CUN#2** = a **mosaic** of CKI_hap1 and CKI_hap2 (recombined in the mother's meiosis)

Along a single chromosome, the ancestry switches partway:

```
CKU_hap1  ■■■■■■■■■■■■■■■■■■■■■■■■
CKU_hap2  □□□□□□□□□□□□□□□□□□□□□□□□
                ↓ meiosis in the father (crossover)
CUN#1     ■■■■■■■■■□□□□□□□□□□□■■■■
                   ↑           ↑
               breakpoint   breakpoint
```

Every chromosome undergoes **at least one crossover, typically one to three**. CUN#1 therefore alternates between stretches resembling CKU_hap1 and stretches resembling CKU_hap2, and **a one-to-one correspondence such as "CUN#1 ≈ CKU_hap1" simply does not hold**. Averaged over the genome, CUN#1 is roughly equally similar to both CKU haplotypes — though the per-chromosome ratio varies widely with crossover position.

The important distinction: what is ambiguous is **which of the parent's haplotypes**, not **which parent**. Every recombined segment still comes from one of the mother's two chromosomes, so the parental assignment itself never wavers.

This is not a theoretical caveat. **The switching is directly observable in the path-similarity analysis of Chapter 7** (§7.5.2): plotting whether CUN#1 is closer to CKU_hap1 or CKU_hap2 along a chromosome reveals the swap at each breakpoint. Seeing it is in fact **evidence that the graph is consistent with the known pedigree**.

---

## 2.5 The significance of trio phasing

Standard diploid genome assembly **cannot distinguish** the paternal and maternal chromosomes—they collapse into a mixed diploid.

But if **both parents are also sequenced**, the offspring's reads can be sorted by parental origin using unique k-mers:

```
Parent A (Kishu) k-mers:  ACGT, GTCA, ...   (Kishu-specific)
Parent B (Kunenbo) k-mers: CCAA, TTAG, ...  (Kunenbo-specific)

Read R1 from child (Satsuma): contains parent-A k-mer → Kishu-derived haplotype (CUNphKi)
Read R2 from child (Satsuma): contains parent-B k-mer → Kunenbo-derived haplotype (CUNphKu)
```

This is **trio phasing** (Cheng et al., 2021).

**Isobe et al. (2023) applied trio phasing to Satsuma** and released the following six haplotype-resolved assemblies:

| Cultivar | Assembly ID | Content |
|---|---|---|
| Satsuma CUN | CUNphKu_r1.0 | Satsuma's **Kunenbo-derived** haplotype (`Ku`) |
| Satsuma CUN | CUNphKi_r1.0 | Satsuma's **Kishu-derived** haplotype (`Ki`) |
| Kishu CKI | CKIhap1_r1.0 | Kishu hap1 |
| Kishu CKI | CKIhap2_r1.0 | Kishu hap2 |
| Kunenbo CKU | CKUhap1_r1.0 | Kunenbo hap1 |
| Kunenbo CKU | CKUhap2_r1.0 | Kunenbo hap2 |

> **On haplotype numbering**: the `#1` / `#2` in PanSN are **arbitrary labels**; nothing dictates which is which. This tutorial uses **`CUN#1` = CUNphKu (paternal, Kunenbo-derived) and `CUN#2` = CUNphKi (maternal, Kishu-derived)**. The order of `hap1_path` / `hap2_path` in `tables/samplesheet.tsv` is what fixes this, and it propagates all the way to the VCF coordinate reference (`-V CUN#1`, §5.2) and to the direction of the pedigree tests in Chapters 6 and 7. **When comparing against someone else's analysis, check first that the numbering matches.**

**Naming convention**:
- **C** = *Citrus*
- Species abbreviation: **UN** = **un**shiu, **KI** = **ki**nokuni, **KU** = **ku**nip
- **ph** = phased (parental origin identified)
- Suffix: **Ki**/**Ku** indicates parental origin (Satsuma only), or **hap1**/**hap2** as arbitrary labels (parental lines)

---

## 2.6 Why this trio is ideal for a pangenome tutorial

This three-cultivar combination is **particularly well-suited** as a teaching dataset:

### 1. The pedigree is fully known

Satsuma = Kishu × Kunenbo has been **genetically confirmed** (Fujii et al., 2016). After building the graph, we can **verify** whether this relationship is visible in the graph structure.

### 2. Appropriate genetic distance

- All three are *Citrus* species, closely related
- Yet Kishu (pure mandarin) and Kunenbo (mandarin + pummelo) are clearly distinguishable
- Satsuma's two haplotypes let us "see" both parents
- The dataset is suitable for both diversity analysis and phylogenetic analysis

### 3. Realistic compute requirements

3 cultivars × 2 haps = **6 haploids**. This is small enough to run in a student environment. Per chromosome takes 30-90 minutes on 32 CPUs; all 9 chromosomes finish in hours when parallelized.

### 4. Fully public data

Isobe et al. (2023) released all assemblies via **Plant GARDEN**. Anyone can reproduce the same tutorial at zero cost.

### 5. A textbook "smallest family"

"Two parents + one F1" is the **simplest unit** in genetics. It maps directly onto Mendel's laws and is a natural gateway to more complex population analyses (GWAS, QTL mapping).

---

## 2.7 The "story of discoveries" in this tutorial

Through the graph on these three cultivars, we will **directly observe or confirm**:

**Discovery 1**: CUN#1 (Kunenbo-derived) is clearly closer to both CKU haplotypes than to CKI
- But it does not track either CKU_hap1 or CKU_hap2 exclusively — the ancestry **switches** along the chromosome (the recombination described in §2.4)

**Discovery 2**: CUN#2 (Kishu-derived) is clearly closer to both CKI haplotypes than to CKU
- Likewise visible as a mosaic of CKI_hap1 and CKI_hap2

**Discovery 3**: SVs exist between Satsuma's two haplotypes
- Thousands of SVs detected via vg deconstruct

**Discovery 4**: Some SVs are inherited from a parent
- The same SVs appear in CKI or CKU

You may also encounter unexpected observations and real-data quirks during QC and analysis—these are also part of the learning experience.

---

## 2.8 A note on why other Citrus cultivars are excluded

PacBio HiFi data has been released for 18 *Citrus* cultivars (Kiryu et al., 2026), but **only CUN/CKI/CKU have public assemblies** as of 2026. The other 15 cultivars (Hassaku, Iyokan, Ponkan, Lemon, Citron, etc.) have HiFi reads only; assemblies are not yet public.

To use those cultivars, one would need to assemble them from raw reads (hifiasm, several hours to a day per cultivar), which is beyond the scope of this tutorial.

We therefore limit ourselves to the **three cultivars whose assemblies are publicly available**.

---

## Chapter summary

- Target: Satsuma mandarin (CUN) = F1 of Kishu (CKI) × Kunenbo (CKU)
- Isobe et al. (2023) released haplotype-resolved trio-phased assemblies
- Six haploids allow direct verification of the family relationship in the graph
- The scale is realistic and this is an "ideal smallest family" for learning pangenome methods

The next chapter describes how to **download** these data.

## References

- Fujii H, et al. (2016). Parental diagnosis of satsuma mandarin (*Citrus unshiu* Marc.) revealed by nuclear and cytoplasmic markers. *Breed Sci* 66:683-691.
- Wu GA, et al. (2018). Genomics of the origin and evolution of *Citrus*. *Nature* 554:311-316.
- Isobe S, et al. (2023). Haploid-resolved and chromosome-scale genome assembly in *Citrus unshiu* and its parental species. *bioRxiv* 2023.06.02.543356.
- Cheng H, et al. (2021). Haplotype-resolved de novo assembly using phased assembly graphs with hifiasm. *Nat Methods* 18:170-175.
- Kiryu Y, et al. (2026). AlleleMiner: a long-read pipeline for gene-wise de novo allele phasing and variant detection in diploid citrus cultivars. *DNA Res* 33:dsag004.

---

[← Chapter 1](01_general_workflow.md) | [Chapter 3: Data acquisition →](03_data_acquisition.md)

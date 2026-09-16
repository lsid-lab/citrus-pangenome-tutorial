# Chapter 2. Target organisms: Satsuma, Kishu, and Kunenbo mandarins

## Learning objectives

- Understand the biology and history of the three target cultivars
- Learn that Satsuma mandarin = Kishu × Kunenbo (F1)
- Appreciate the significance of trio phasing and its relevance here
- Understand why these three cultivars are an ideal pangenome tutorial dataset

---

## 2.1 Satsuma mandarin (*Citrus unshiu*)

**Satsuma mandarin** (abbreviated STS: **S**a**ts**uma), Japan's iconic winter citrus, was discovered in Nagashima, Kagoshima Prefecture during the early Edo period (17th century).

- **Chromosome number**: 2n = 18 (n = 9)
- **Genome size**: ~320-360 Mb (haploid)
- **Reproduction**: seedless (parthenocarpic), with nucellar polyembryony (seeds contain maternal clones)
- **Cultivars**: Miyagawa wase, Okitsu wase, Aoshima, Nankan No. 20, and many more

Seedlessness is a benefit for cultivation but creates **challenges for genetic analysis**: standard crossing experiments (parent A × parent B → observe F1) are not feasible.

The origin of Satsuma was long uncertain, but **recent genomic analyses** (Fujii et al., 2016; Wu et al., 2018) confirmed that it is an **F1 hybrid of Kishu mandarin (KSH) × Kunenbo (KNN)**.

### Genetic composition

According to Wu et al. 2018, Satsuma is composed of:
- **Mandarin (*C. reticulata*)** gene pool: over 80%
- Small **Pummelo (*C. maxima*)** gene pool contribution

This reflects admixture inherited from both parents.

---

## 2.2 Kishu mandarin (*Citrus kinokuni*)

**Kishu mandarin** (KSH: **K**i**sh**u) is the **maternal parent** of Satsuma.

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

**Kunenbo** (KNN: **Kun**enbo) is the **paternal parent** of Satsuma.

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
        Kishu (KSH)                  Kunenbo (KNN)
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
              Satsuma (STS)
              [Mandarin ~90%, Pummelo ~10%]
        [hap1 = Kishu-derived | hap2 = Kunenbo-derived]
```

**Key point**: Satsuma's two haplotypes were each **inherited from one parent**:
- STS_hap1 ≈ one of Kishu's haplotypes (or a recombinant)
- STS_hap2 ≈ one of Kunenbo's haplotypes (or a recombinant)

---

## 2.5 The significance of trio phasing

Standard diploid genome assembly **cannot distinguish** the paternal and maternal chromosomes—they collapse into a mixed diploid.

But if **both parents are also sequenced**, the offspring's reads can be sorted by parental origin using unique k-mers:

```
Parent A (Kishu) k-mers:  ACGT, GTCA, ...   (Kishu-specific)
Parent B (Kunenbo) k-mers: CCAA, TTAG, ...  (Kunenbo-specific)

Read R1 from child (Satsuma): contains parent-A k-mer → hap1 (Kishu-derived)
Read R2 from child (Satsuma): contains parent-B k-mer → hap2 (Kunenbo-derived)
```

This is **trio phasing** (Cheng et al., 2021).

**Isobe et al. (2023) applied trio phasing to Satsuma** and released the following six haplotype-resolved assemblies:

| Cultivar | Assembly ID | Content |
|---|---|---|
| Satsuma STS | CUNphKi_r1.0 | Satsuma hap1 (**Ki**shu-derived) |
| Satsuma STS | CUNphKu_r1.0 | Satsuma hap2 (**Ku**nenbo-derived) |
| Kishu KSH | CKIhap1_r1.0 | Kishu hap1 |
| Kishu KSH | CKIhap2_r1.0 | Kishu hap2 |
| Kunenbo KNN | CKUhap1_r1.0 | Kunenbo hap1 |
| Kunenbo KNN | CKUhap2_r1.0 | Kunenbo hap2 |

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

**Discovery 1**: Satsuma's hap1 is particularly similar to one of Kishu's haplotypes
- High similarity between STS_hap1 and either KSH_hap1 or KSH_hap2

**Discovery 2**: Satsuma's hap2 is particularly similar to a Kunenbo haplotype
- Symmetric to Discovery 1

**Discovery 3**: SVs exist between Satsuma's two haplotypes
- Thousands of SVs detected via vg deconstruct

**Discovery 4**: Some SVs are inherited from a parent
- The same SVs appear in KSH or KNN

You may also encounter unexpected observations and real-data quirks during QC and analysis—these are also part of the learning experience.

---

## 2.8 A note on why other Citrus cultivars are excluded

PacBio HiFi data has been released for 18 *Citrus* cultivars (Kiryu et al., 2026), but **only STS/KSH/KNN have public assemblies** as of 2026. The other 15 cultivars (Hassaku, Iyokan, Ponkan, Lemon, Citron, etc.) have HiFi reads only; assemblies are not yet public.

To use those cultivars, one would need to assemble them from raw reads (hifiasm, several hours to a day per cultivar), which is beyond the scope of this tutorial.

We therefore limit ourselves to the **three cultivars whose assemblies are publicly available**.

---

## Chapter summary

- Target: Satsuma mandarin (STS) = F1 of Kishu (KSH) × Kunenbo (KNN)
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

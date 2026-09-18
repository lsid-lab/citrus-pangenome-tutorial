# Chapter 3. Data acquisition

## Learning objectives

- Trace data locations from a paper's "Data availability" section
- Distinguish **BioProject-registered data** from **actual assemblies**
- Download assemblies from Plant GARDEN
- Set up the project directory structure

---

## 3.1 Finding data: the paper's "Data availability" section

The first step in reproducing a paper is to read the **"Data availability" section**.

The data used in this tutorial comes from:

> Isobe S, et al. (2023). Haploid-resolved and chromosome-scale genome assembly in *Citrus unshiu* and its parental species, *C. nobilis* and *C. kinokuni*. *bioRxiv* 2023.06.02.543356.
> <https://www.biorxiv.org/content/10.1101/2023.06.02.543356v1>

The relevant excerpt from its "Data availability" section:

> The sequence reads are available from the DNA Data Bank of Japan (DDBJ) Sequence Read Archive (DRA) under the BioProject number **PRJDB15866**. The assembled scaffold sequences, gene sequences, and annotation files are available at **Plant GARDEN**.

Two crucial pieces of information here:

1. **Raw HiFi reads** are on **DDBJ DRA** (BioProject PRJDB15866)
2. **Assemblies**, gene sequences, and annotation are on **Plant GARDEN**

**If you want assemblies, Plant GARDEN is the right place.** DDBJ does not host the assemblies for this project.

---

## 3.2 Background: sequence data archives

Life-science sequence data is deposited in **international archives**:

| Archive | Operator | Country | Content |
|---|---|---|---|
| **NCBI SRA** | NCBI | USA | Raw reads, assemblies |
| **EBI ENA** | EMBL-EBI | Europe | Raw reads, assemblies |
| **DDBJ DRA** | ROIS | Japan | Raw reads, assemblies |
| **Plant GARDEN** | Kazusa DNA Research Institute | Japan | Plant genomes (assemblies + annotation) |

Under the **INSDC (International Nucleotide Sequence Database Collaboration)** agreement, NCBI SRA / EBI ENA / DDBJ DRA are **mutually mirrored**: registration in one is discoverable in all three.

**Plant GARDEN is not part of INSDC** and operates as an independent portal. Many completed plant genome assemblies are consolidated there.

---

## 3.3 BioProject and BioSample relationships

DDBJ/NCBI/ENA registrations are hierarchical:

```
BioProject (whole project)
  │
  ├── BioSample #1 (individual/cultivar #1)
  │     ├── Experiment #1 (one sequencing experiment)
  │     │     └── Run #1 (actual read data)
  │     └── Experiment #2
  │           └── Run #2
  ├── BioSample #2
  │     ├── ...
```

For `PRJDB15866`:
- BioProject = the Isobe 2023 Satsuma + Kishu + Kunenbo project
- Three BioSamples:
  - SAMD00608599 (Satsuma, Miyagawa Wase)
  - SAMD00608600 (Kunenbo, Kagoshima)
  - SAMD00608601 (Kishu mikan)
- Each BioSample has associated Runs for **HiFi WGS**, **short-read WGS**, **Hi-C**, and **Iso-Seq**

To assemble from raw reads, you would download FASTQ files via these Run accessions (not covered here).

---

## 3.4 Downloading assemblies from Plant GARDEN

### 3.4.1 What is Plant GARDEN?

**Plant GARDEN** (https://plantgarden.jp) is a **plant genome portal** operated by Kazusa DNA Research Institute.

- Hosts assemblies for dozens of plant species
- Also provides gene predictions, annotations, and k-mer databases
- Free for academic use

Citrus data paths:

- Satsuma: `https://plantgarden.jp/en/download/Citrus_unshiu/`
- Kishu: `https://plantgarden.jp/en/download/Citrus_kinokuni/`
- Kunenbo: `https://plantgarden.jp/en/download/Citrus_nobilis/`

### 3.4.2 Reading the Assembly IDs (Gxxx)

Plant GARDEN publishes **multiple assembly versions** for the same species in parallel, each under an **Assembly ID** of the form `G001`, `G002`, ….

For each of the three cultivars, Isobe et al. (2023) deposited the **same set of three**, numbered in this order:

1. **unphased** — the diploid collapsed into a single set
2. **hap1** — phased haplotype 1
3. **hap2** — phased haplotype 2

> The hap1 / hap2 here is **the distributor's (Plant GARDEN's) deposit order**, which is **not** the same as the PanSN `#1` / `#2` this tutorial uses in the graph. For Satsuma, deposit-order hap1 (G003, CUNphKi) becomes this tutorial's `CUN#2` (see §2.5 for why).

**Satsuma (CUN, taxon `t55188`)**

| Assembly ID | File | Content |
|---|---|---|
| t55188.**G001** | `C_unshiu_v1.0_scaffolds.fa.gz` | A **different study** — Kawahara et al. (2020), older and unphased |
| t55188.**G002** | `CUNuph_r1.0.fasta.gz` | Isobe 2023 **unphased** (`uph` = **u**n**ph**ased) |
| t55188.**G003** | `CUNphKi_r1.0.pmol.fasta.gz` | **Kishu-derived haplotype** (`ph` = **ph**ased, `Ki` = **Ki**shu) ← used here (`CUN#2`) |
| t55188.**G004** | `CUNphKu_r1.0.ch1-9.fasta.gz` | **Kunenbo-derived haplotype** (`Ku` = **Ku**nenbo) ← used here (`CUN#1`) |

**Kishu (CKI, `t408488`) and Kunenbo (CKU, `t481549`)** — these only go up to `G003`.

| Assembly ID | Content | File |
|---|---|---|
| **G001** | **unphased** (not used here) | — |
| **G002** | **hap1** ← used here | `CKIhap1_r1.0.pmol.fasta.gz` / `CKUhap1_r1.0.pmol.fasta.gz` |
| **G003** | **hap2** ← used here | `CKIhap2_r1.0.pmol.fasta.gz` / `CKUhap2_r1.0.pmol.fasta.gz` |

> **Why Satsuma's numbering is offset by one**: an assembly from Kawahara et al. (2020) was already registered as `G001` for Satsuma before Isobe 2023. So "Isobe 2023's unphased / hap1 / hap2" maps to **G002 / G003 / G004** for Satsuma but **G001 / G002 / G003** for the two parents. Rather than memorizing the numbers, remember the **unphased → hap1 → hap2 ordering**.

**The six files this tutorial uses are the hap1 and hap2 of each cultivar** (G003/G004 for Satsuma, G002/G003 for Kishu and Kunenbo).

The unphased assemblies (Satsuma G002, and G001 for the parents) are not used: **with both haplotypes collapsed into one sequence, the differences between haplotypes cannot be represented in the graph.** Conversely, it is precisely because phased assemblies exist for all three cultivars that building a pangenome here is worthwhile.

### 3.4.3 What does `pmol` mean?

The **`pmol`** that recurs in these filenames is short for **pseudomolecule**.

It is called "pseudo" because it is **not a single molecule read end-to-end across a whole chromosome**. It is built like this:

1. Sequencing and assembly produce many contigs / scaffolds
2. Hi-C data, linkage maps, or similar are used to **order and orient them along the chromosome**
3. The ordered pieces are joined (padding gaps with `N`) into **one sequence representing one chromosome**

So a `pmol` is a sequence **constructed** to represent a chromosome. Treating it as a chromosome sequence is fine in practice, but keep in mind that it is **not an observed single molecule** — gaps (runs of `N`) may remain, and the order or orientation of some pieces may be wrong. That is part of why Chapter 4's QC checks the fraction of `N`.

The convention is used consistently across Plant GARDEN / Kazusa DNA Research Institute releases: tomato's `SLM_r2.0.pmol`, pepper's `CAN_r1.2.pmol`, hydrangea's `HMA_r1.2.pmol`, and so on. A useful rule of thumb when browsing other Plant GARDEN datasets: **pick the file with `pmol` in its name and you get chromosome-level sequence.**

> Note that Satsuma's Kunenbo-derived haplotype (`CUN#1`) is the one exception here: it is named **`CUNphKu_r1.0.ch1-9.fasta.gz`** (chromosomes 1-9 only) rather than `.pmol`. Since this tutorial builds graphs per chromosome, either naming works fine.

---

## 3.5 Project directory layout

This tutorial is designed to be **cloned from GitHub**:

```bash
git clone https://github.com/lsid-lab/citrus-pangenome-tutorial.git
cd citrus-pangenome-tutorial
```

After cloning, you already have:

```
citrus-pangenome-tutorial/
├── README.md
├── docs/en/            # This tutorial
├── scripts/            # All executable scripts
├── tables/             # Metadata (samplesheet, etc.)
│   ├── samplesheet.tsv
│   └── assembly_provenance_v1_v2.tsv   # r1.0 / r2.0 provenance, SHA-256, stats
└── LICENSE
```

**The only directory you need to create is `data/raw/`**, where the assemblies go:

```bash
mkdir -p data/raw/CUN data/raw/CKI data/raw/CKU
```

Intermediate files and outputs (`results/`, etc.) are created automatically by subsequent scripts.

Recommended layout:

```
citrus-pangenome-tutorial/
├── ... (files from the clone)
├── data/
│   ├── raw/                       # * The distributed files, untouched (you provide these)
│   │   ├── CUN/
│   │   │   ├── CUNphKu_r1.0.ch1-9.fasta.gz   # -> CUN#1
│   │   │   └── CUNphKi_r1.0.pmol.fasta.gz    # -> CUN#2
│   │   ├── CKI/
│   │   │   ├── CKIhap1_r1.0.pmol.fasta.gz
│   │   │   └── CKIhap2_r1.0.pmol.fasta.gz
│   │   └── CKU/
│   │       ├── CKUhap1_r1.0.pmol.fasta.gz
│   │       └── CKUhap2_r1.0.pmol.fasta.gz
│   └── input/                     # Built by the Chapter 5 preprocessing (PanSN, per chromosome)
│       └── chr01.fa.gz ... chr09.fa.gz
├── results/                       # Analysis outputs (created by the scripts)
│   ├── qc/                        # Chapter 4: assembly statistics
│   ├── pggb/                      # Chapter 5: per-chromosome graphs
│   ├── graph_qc/                  # Chapter 6
│   └── viz/                       # Chapter 7
└── bin/                           # Singularity wrappers (optional, Chapter 5)
```

---

## 3.6 Downloading in practice

### 3.6.1 Manual download

You can download each file through the Plant GARDEN web interface. Files are ~100 MB each, so command-line downloads are more convenient.

### 3.6.2 Batch download script

Use `scripts/download_plantgarden.sh` (continuing from the repository root you `cd`'d into in §3.5):

```bash
bash scripts/download_plantgarden.sh data/
```

This script:
1. Downloads the six assemblies from Plant GARDEN
2. Places them under `data/raw/{CUN,CKI,CKU}/`
3. (Optionally) verifies integrity with SHA256 checksums

### 3.6.3 Expected file sizes

| File | Size (bgzip) |
|---|---:|
| CUNphKi_r1.0.pmol.fasta.gz | ~101 MB |
| CUNphKu_r1.0.ch1-9.fasta.gz | ~104 MB |
| CKIhap1_r1.0.pmol.fasta.gz | ~88 MB |
| CKIhap2_r1.0.pmol.fasta.gz | ~91 MB |
| CKUhap1_r1.0.pmol.fasta.gz | ~94 MB |
| CKUhap2_r1.0.pmol.fasta.gz | ~88 MB |
| **Total** | **~570 MB** |

---

## 3.7 Verifying the samplesheet

All subsequent scripts read `tables/samplesheet.tsv` as input. **This file is provided in the clone**, but check that paths match your setup.

Contents of `tables/samplesheet.tsv`:

```tsv
sample_id  cultivar_jp  species          pansn_prefix  hap1_path                                        hap2_path                                          source                            pedigree
CUN        Satsuma      Citrus unshiu    CUN       data/raw/CUN/CUNphKu_r1.0.ch1-9.fasta.gz        data/raw/CUN/CUNphKi_r1.0.pmol.fasta.gz           Plant GARDEN t55188.G004/G003    F1: CKI x CKU
CKI        Kishu        Citrus kinokuni  CKI         data/raw/CKI/CKIhap1_r1.0.pmol.fasta.gz           data/raw/CKI/CKIhap2_r1.0.pmol.fasta.gz             Plant GARDEN t408488             mother of CUN
CKU        Kunenbo      Citrus nobilis   CKU       data/raw/CKU/CKUhap1_r1.0.pmol.fasta.gz         data/raw/CKU/CKUhap2_r1.0.pmol.fasta.gz           Plant GARDEN t481549             father of CUN
```

**Column meanings**:
- `sample_id`: internal short code (3 letters)
- `cultivar_jp`: display name
- `species`: scientific name
- `pansn_prefix`: short prefix used later in PanSN names
- `hap1_path` / `hap2_path`: relative paths to FASTA files (from repo root). **This order determines the PanSN `#1` / `#2`** (for Satsuma, `#1` = CUNphKu and `#2` = CUNphKi; see §2.5)
- `source`: data origin
- `pedigree`: family relationship notes

Keeping this file correct is the foundation for all subsequent steps.

---

## 3.8 After downloading

Once downloads complete, verify:

```bash
# File sizes and presence
ls -la data/raw/*/*.fasta.gz

# Peek at each FASTA
zcat data/raw/CUN/CUNphKi_r1.0.pmol.fasta.gz | head -3

# Expected output:
# >CUNphKi_r1.0ch1
# GCTAGCTAGCTAGCTAGCTAGCTAGC...
# TAGCTAGCTAGCTAGCTAGCTAGCTA...
```

**Verification checklist**:
- All six files downloaded
- File sizes match expected values (dozens to ~100+ MB)
- FASTA headers start with `>` and use names like `CUNphKi_r1.0ch1`

---

## Chapter summary

- The paper's "Data availability" section is your first clue
- Use BioProject numbers for **raw reads on DDBJ**, Plant GARDEN for **assemblies**
- Plant GARDEN naming pattern: `C{species}{ph or hap}{Ki/Ku or 1/2}`
- **Assembly IDs (Gxxx) run unphased → hap1 → hap2**; Satsuma is offset by one because an older assembly occupies `G001`
- **`pmol` = pseudomolecule** — a sequence *constructed* to represent a chromosome by ordering contigs along it, not a molecule read end-to-end
- Only the **phased (hap1/hap2)** assemblies go into the pangenome; unphased ones cannot represent differences between haplotypes
- One samplesheet (TSV) manages 3 cultivars × 2 haps = 6 haploids

The next chapter walks through **quality control** on the downloaded assemblies.

## References

- **Plant GARDEN**: https://plantgarden.jp
- **DDBJ**: https://www.ddbj.nig.ac.jp/
- **NCBI**: https://www.ncbi.nlm.nih.gov/
- **INSDC**: http://www.insdc.org/
- Isobe S, et al. (2023). Haploid-resolved and chromosome-scale genome assembly in *Citrus unshiu* and its parental species, *C. nobilis* and *C. kinokuni*. *bioRxiv* 2023.06.02.543356. <https://www.biorxiv.org/content/10.1101/2023.06.02.543356v1>
- Kawahara Y, et al. (2020). Mikan Genome Database (MiGD): integrated database of genome annotation, genomic diversity, and CAPS marker information for mandarin molecular breeding. *Breed Sci* 70(2). (source of Satsuma's `G001` assembly)

---

[← Chapter 2](02_target_organisms.md) | [Chapter 4: Data QC →](04_data_qc.md)

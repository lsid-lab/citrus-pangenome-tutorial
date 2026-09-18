# Chapter 3. Data acquisition

## Learning objectives

- Trace data locations from a paper's "Data availability" section
- Distinguish **BioProject-registered data** from **actual assemblies**
- Download assemblies from Plant GARDEN
- Set up the project directory structure

---

## 3.1 Finding data: the paper's "Data availability" section

The first step in reproducing a paper is to read the **"Data availability" section**.

The relevant excerpt from Isobe et al. (2023) *bioRxiv*:

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

### 3.4.2 File naming convention

Plant GARDEN publishes **multiple assembly versions** for the same species in parallel. For Satsuma:

| Assembly ID | File | Content |
|---|---|---|
| t55188.**G001** | `C_unshiu_v1.0_scaffolds.fa.gz` | Kawahara 2020 hybrid (older, unphased) |
| t55188.**G002** | `CUNuph_r1.0.fasta.gz` | Isobe 2023 phased (integrated) |
| t55188.**G003** | `CUNphKi_r1.0.pmol.fasta.gz` | **CUN hap1 (Kishu-derived)** |
| t55188.**G004** | `CUNphKu_r1.0.ch1-9.fasta.gz` | **CUN hap2 (Kunenbo-derived)** |

**We use G003 and G004** (and corresponding phased assemblies for Kishu and Kunenbo) in this tutorial.

---

## 3.5 Project directory layout

This tutorial is designed to be **cloned from GitHub**:

```bash
git clone https://github.com/<user>/citrus-pangenome-tutorial.git
cd citrus-pangenome-tutorial
```

After cloning, you already have:

```
citrus-pangenome-tutorial/
├── README.md
├── docs/en/            # This tutorial
├── scripts/            # All executable scripts
├── tables/             # Metadata (samplesheet, etc.)
│   └── samplesheet.tsv
└── LICENSE
```

**The only directory you need to create is `data/`**, where the assemblies will be placed:

```bash
mkdir -p data/CUN data/CKI data/CKU
```

Intermediate files and outputs (`03_pangenome/`, etc.) are created automatically by subsequent scripts.

Recommended layout:

```
citrus-pangenome-tutorial/
├── ...(files present at clone time)
├── data/                              # (you provide) Downloaded assemblies
│   ├── CUN/
│   │   ├── CUNphKi_r1.0.pmol.fasta.gz
│   │   └── CUNphKu_r1.0.ch1-9.fasta.gz
│   ├── CKI/
│   │   ├── CKIhap1_r1.0.pmol.fasta.gz
│   │   └── CKIhap2_r1.0.pmol.fasta.gz
│   └── CKU/
│       ├── CKUhap1_r1.0.pmol.fasta.gz
│       └── CKUhap2_r1.0.pmol.fasta.gz
└── 03_pangenome/                      # (created by later scripts)
    ├── by_chr/
    ├── qc/
    └── logs/
```

---

## 3.6 Downloading in practice

### 3.6.1 Manual download

You can download each file through the Plant GARDEN web interface. Files are ~100 MB each, so command-line downloads are more convenient.

### 3.6.2 Batch download script

Use `scripts/download_plantgarden.sh`:

```bash
cd citrus-pangenome-tutorial
bash scripts/download_plantgarden.sh data/
```

This script:
1. Downloads the six assemblies from Plant GARDEN
2. Places them under `data/{CUN,CKI,CKU}/`
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
CUN        Satsuma      Citrus unshiu    CUN       data/CUN/CUNphKi_r1.0.pmol.fasta.gz         data/CUN/CUNphKu_r1.0.ch1-9.fasta.gz          Plant GARDEN t55188.G003/G004    F1: CKI x CKU
CKI        Kishu        Citrus kinokuni  CKI         data/CKI/CKIhap1_r1.0.pmol.fasta.gz           data/CKI/CKIhap2_r1.0.pmol.fasta.gz             Plant GARDEN t408488             mother of CUN
CKU        Kunenbo      Citrus nobilis   CKU       data/CKU/CKUhap1_r1.0.pmol.fasta.gz         data/CKU/CKUhap2_r1.0.pmol.fasta.gz           Plant GARDEN t481549             father of CUN
```

**Column meanings**:
- `sample_id`: internal short code (3 letters)
- `cultivar_jp`: display name
- `species`: scientific name
- `pansn_prefix`: short prefix used later in PanSN names
- `hap1_path` / `hap2_path`: relative paths to FASTA files (from repo root)
- `source`: data origin
- `pedigree`: family relationship notes

Keeping this file correct is the foundation for all subsequent steps.

---

## 3.8 After downloading

Once downloads complete, verify:

```bash
# File sizes and presence
ls -la data/*/*.fa.gz

# Peek at each FASTA
zcat data/CUN/CUNphKi_r1.0.pmol.fasta.gz | head -3

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
- One samplesheet (TSV) manages 3 cultivars × 2 haps = 6 haploids

The next chapter walks through **quality control** on the downloaded assemblies.

## References

- **Plant GARDEN**: https://plantgarden.jp
- **DDBJ**: https://www.ddbj.nig.ac.jp/
- **NCBI**: https://www.ncbi.nlm.nih.gov/
- **INSDC**: http://www.insdc.org/
- Isobe S, et al. (2023). bioRxiv 2023.06.02.543356.

---

[← Chapter 2](02_target_organisms.md) | [Chapter 4: Data QC →](04_data_qc.md)

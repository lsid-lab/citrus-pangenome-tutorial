# Chapter 7. Graph interpretation

## Learning objectives

- Read variant information from a VCF file
- Visualize a pangenome graph with odgi
- Verify pedigree with real data
- Interpret the results biologically and understand the tutorial's limits
- **See that an assembly is a versioned artifact, learn what to do with an observation you cannot explain, and learn to measure how robust a conclusion is**

---

## 7.1 Goals of this chapter

In previous chapters, we verified technical quality (Level 1-2) and biological validity (Level 3-4). This chapter enters the phase of **actually reading the biological meaning** of those results.

Specifically:

1. **Reading VCF**: what variants are captured in the graph
2. **Visualization**: understanding graph structure visually
3. **Pedigree verification**: how "F1 = parent A × parent B" appears in the graph
4. **Tutorial limits**: what this graph can and cannot show

---

## 7.2 Structure and reading of VCF files

### 7.2.1 VCF basics

VCF (Variant Call Format) is the standard for describing **variant positions and types** relative to a reference.

Example (excerpt from chr09):
```
#CHROM   POS      ID    REF   ALT      QUAL   FILTER   INFO   FORMAT   CUN#1  CUN#2  CKI#1  ...
CUN#1#chr09   1234  .     A     G        .      .        AT=snp   GT       0          1          0        ...
CUN#1#chr09   2000  .     ACGT  A        .      .        AT=del3  GT       0          1          1        ...
CUN#1#chr09   3500  .     C     CGATCG   .      .        AT=ins5  GT       0          0          1        ...
```

Columns:
- **CHROM**: coordinate reference path (here `CUN#1#chr09`)
- **POS**: position (1-based)
- **REF**: reference (coordinate anchor) base
- **ALT**: alternative allele
- **GT** (Genotype): per-sample genotype
  - `0` = REF, `1` = ALT1, `2` = ALT2, `.` = missing

### 7.2.2 Classifying SNPs, indels, and SVs

Classify variants from REF and ALT lengths:

```bash
# SNP: both REF and ALT are 1 bp
awk 'length($4)==1 && length($5)==1' chr09.vcf

# Insertion: REF short, ALT longer
awk 'length($4)<length($5) && length($5)>1' chr09.vcf

# Deletion: REF long, ALT short
awk 'length($4)>length($5) && length($4)>1' chr09.vcf

# SV: either side ≥ 50 bp
awk 'length($4)>=50 || length($5)>=50' chr09.vcf
```

### 7.2.3 Per-sample variant counts

How many ALT variants does each cultivar carry?

```bash
# ALT variants for CUN#2 (the other hap from the coordinate anchor)
bcftools view -s CUN#2 chr09.vcf | \
  bcftools query -f '[%GT]\n' | grep -c '^1'
```

### 7.2.4 Expected results

Variants found between the six haploids (chr09, ~30 Mb):

|Variant type|Expected count|
|---|---|
|SNPs|10,000-30,000|
|Small indels (<50 bp)|5,000-15,000|
|SVs (≥50 bp)|100-500|

---

## 7.3 Pangenome graph visualization

### 7.3.1 1D linear view (odgi viz)

Places each path on a horizontal axis, pileup-style:

```bash
odgi viz -i chr09.smooth.final.og -o chr09_viz.png -x 1500 -y 500
```

Interpretation:
- **Horizontal axis**: position on the graph (coordinate anchor path)
- **Vertical axis**: each path (6 haploids)
- **Color**: same color means paths share the same graph region
- **Color changes**: branches (= variants)

Signs of a healthy graph:
- OK: six paths run in horizontal parallel lines
- OK: broad regions have consistent color (shared sequence)
- OK: sporadic vertical color changes (variants)

Signs of problems:
- FAIL: hundreds of tiny vertical lines → over-fragmentation
- FAIL: paths that vanish mid-way → preservation failure

### 7.3.2 2D layout view (odgi layout + draw)

Spatial arrangement of the whole graph structure:

```bash
odgi layout -i chr09.smooth.final.og -o chr09.lay -t 8
odgi draw -i chr09.smooth.final.og -c chr09.lay -p chr09_2d.png -H 1000 -C
```

Interpretation:
- **backbone**: thick continuous line where all paths are shared
- **bubble**: branches carried only by some paths (small SVs)
- **long branch**: substantial private sequences

A **healthy pangenome** typically shows a thick central backbone with many small bubbles around it.

An **unhealthy pangenome** may show six intertwined paths without a clear backbone.

### 7.3.3 Execution script

`scripts/pg04_visualize.sh` runs both:

```bash
bash scripts/pg04_visualize.sh . chr09
```

---

## 7.4 Pedigree verification in practice

### 7.4.1 Verifying Satsuma's trio phasing on the graph

Extending the Chapter 6 Level 3 pedigree test with biological interpretation.

Compute from `similarity` output:

|Comparison|Expected Jaccard|
|---|---|
|`CUN#1` vs CKU_hap1 / CKU_hap2|**High against both** — `CUN#1` is a mosaic of the two, so neither one stands out alone|
|`CUN#1` vs CKI_hap*|**Relatively low** — the mother contributed nothing to `CUN#1`|

Symmetrically, `CUN#2` (Kishu-derived) scores high against both CKI haplotypes and low against CKU.

> Absolute values are depressed simply because CUN's paths are longer (§4.4, §6.4.4), so judge by the **ranking between cultivars, not against a threshold**. Note also that a single whole-chromosome Jaccard cannot tell you *which* parental haplotype a segment descends from (§2.4) — that needs the windowed analysis in §7.5.2.

### 7.4.2 Specific pedigree pattern

Trace which parent `CUN#1` came from:

```bash
# Jaccard of CUN#1 vs each parent haplotype
awk '$1 ~ /^CUN#1/ && ($2 ~ /^CKU/ || $2 ~ /^CKI/)' chr09_similarity.tsv
```

If the output shows:
```
CUN#1#chr09  CKU#1#chr09  ...  Jaccard = 0.68
CUN#1#chr09  CKU#2#chr09  ...  Jaccard = 0.61
CUN#1#chr09  CKI#1#chr09  ...  Jaccard = 0.45
CUN#1#chr09  CKI#2#chr09  ...  Jaccard = 0.43
```

then, because the CKU side is consistently higher than the CKI side, we can conclude **`CUN#1` is paternally (CKU, Kunenbo) derived**.

What we must **not** conclude is "CKU_hap1 scores higher, therefore `CUN#1` came from CKU_hap1". `CUN#1` is a mosaic of CKU_hap1 and CKU_hap2 (§2.4); the gap only means that CKU_hap1-derived segments happen to make up a somewhat larger share of that mosaic. Which segment came from which is only visible in the windowed analysis below.

### 7.4.3 Accounting for CUN's longer paths

The Chapter 4 observation that CUN's haplotypes are larger than the parents' (§4.4) manifests as:

- `CUN#1` vs CKU Jaccard is somewhat lower than expected (e.g., 0.6 instead of 0.7)
- Average-based interpretation may appear to break pedigree consistency

**MAX-based** interpretation confirms pedigree, and the graph itself is correctly built. This was the key point of §6.4.4.

---

## 7.5 Examples of biological insight

Biological patterns that may emerge from the graph:

### 7.5.1 Pummelo-derived regions in Satsuma

According to Wu et al. 2018, Satsuma is mainly mandarin but carries some pummelo ancestry. If in the graph:

- Some path regions in CUN do not match CKI or CKU and form independent paths

these may be pummelo-derived. Pummelo-derived regions are concentrated on specific chromosomal segments (Fujii et al. 2016), so visualization may spot them.

### 7.5.2 Recombination in F1

As described in §2.4, `CUN#1` is a **recombinant mosaic** of CKU_hap1 and CKU_hap2. Computing similarity in windows along a chromosome, we see:

- Regions where `CUN#1` is close to CKU_hap1
- Regions where `CUN#1` is close to CKU_hap2

alternating, with the boundaries marking **crossover breakpoints** — the positions of the actual crossovers that occurred in the mother's meiosis, read off the graph.

> Be careful, though: not every switch you see here is a genuine crossover. The same signature — "from this position on, the closer relative changes" — also appears when **the parental assembly's own phasing switches**. Windowed similarity alone cannot separate the two. §7.6 shows a case where this distinction actually bites.

### 7.5.3 Repeats and TEs

Many large SVs (5-15 kb) are likely **TE insertions**. Kiryu et al. 2026 report a 5.1 kb LTR retrotransposon insertion. Extracting SVs from VCF and BLASTing against TE databases enables this analysis (beyond the scope of this tutorial).

---

## 7.6 As it turns out — there is a v2 of this data

Two things should have snagged your attention along the way:

1. **CUN's two haplotypes are 25-55 Mb larger than the two parents** (§4.4)
2. **CUN's Jaccard similarities are depressed because its paths are longer** (§6.4.4, §7.4.3)

Each on its own is the kind of oddity you can shrug off, and this tutorial has done nothing with them beyond recording them.

### A version 2 exists

The Plant GARDEN assemblies used throughout this tutorial are **r1.0**. An **r2.0** of the same three cultivars was released on 23 March 2026 by NARO on **MiGD2**.

> **MiGD2 (Mikan Genome Database 2)**
> Home: <https://mikan.dna.naro.go.jp/migd2/>
> Downloads: <https://mikan.dna.naro.go.jp/migd2/data_download/download.html>

A comparison exists in which graphs were rebuilt from both releases **with the same pipeline and the same parameters**. Here is what it shows.

### (1) The size gap disappears

Bases in **chromosomes 1-9 only** — what actually goes into the graph.

| PanSN | Assembly | r1.0 (Plant GARDEN) | r2.0 (MiGD2) | Δ |
|---|---|---:|---:|---:|
| `CUN#1` | CUNphKu (Satsuma, Kunenbo-derived) | 357.6 Mb | 300.4 Mb | −57.2 |
| `CUN#2` | CUNphKi (Satsuma, Kishu-derived) | 348.5 Mb | 298.7 Mb | −49.8 |
| `CKI#1` | CKIhap1 (Kishu 1) | 304.2 Mb | 295.6 Mb | −8.6 |
| `CKI#2` | CKIhap2 (Kishu 2) | 310.3 Mb | 297.9 Mb | −12.4 |
| `CKU#1` | CKUhap1 (Kunenbo 1) | 323.9 Mb | 295.1 Mb | −28.8 |
| `CKU#2` | CKUhap2 (Kunenbo 2) | 303.4 Mb | 306.4 Mb | +3.0 |
| | **Spread across the six** | **303.4-357.6 (54.3 Mb)** | **295.1-306.4 (11.3 Mb)** | |

**In r2.0 all six fall within 295-306 Mb, and Satsuma no longer stands apart.** The "25-55 Mb larger than the parents" recorded in §4.4 simply does not hold for r2.0.

Two cautions when reading this:

- **Satsuma is not the only one that shrank.** Kunenbo hap1 lost 28.8 Mb too. The right reading is that r2.0 as a whole is tighter.
- **The lost bases did not simply vanish.** r2.0 distributes unplaced contigs in the same file as the chromosomes (78-265 contigs, 7.3-21.9 Mb per haplotype), so some of what left the chromosomes is there. r1.0's unplaced contigs were distributed separately and were not obtained for this tutorial, so **totals cannot be compared between the releases** — which is why the table above is restricted to chr1-9.

Oddity 2 resolves as a consequence: once all six paths are comparable in length, the "longer path, lower Jaccard" asymmetry from §6.4.4 stops biting.

### (2) The ancestry mosaic reproduces almost unchanged

This is the interesting part. **Despite all that change in size, the ancestry mosaic from §7.5.2 does not move.**

Ancestry was called independently for each release and compared window by window (100 kb windows):

- Across the 16 track × chromosome pairs where the parental phase is stable between releases: **3268 / 3300 windows agree (99.0%)**
- Across all 18 pairs, including the two exceptions: **3539 / 3733 windows agree (94.8%)**

The two exceptions are ch1 and ch9 of `CUN#2` (CUNphKi, Kishu-derived), where **the Kishu parental assembly's own phase is reorganized between releases**. The disagreement is therefore a property of the parental assemblies, not of the Satsuma ancestry call — exactly the caveat flagged in §7.5.2, actually happening.

Calibrating against an independent yardstick — a RAD-Seq F1 population of 96 real offspring of Kunenbo × Kishu — **Satsuma falls inside the distribution of genuine offspring under both releases, at the same percentile (93.8).**

**So one cannot say "r1.0's phasing was wrong and r2.0 fixed it."** There is no independent ground truth for Satsuma, and the honest conclusion is that **which release is closer to biological reality cannot be determined** from this comparison.

### (3) Why the comparison holds up — change one thing at a time

Across the two releases, **all 36 scientific parameters and tool versions were identical** (pggb / wfmash / seqwish / smoothxg / vg versions, `-n 6`, `-p 95`, `-s 10000`, `-k 23`, `-B 1000000000`, `-V CUN#1`).

That matters because it **eliminates the most obvious confounder**. Running the same binaries with the same arguments means no difference in the results can be attributed to the pipeline. The only thing that differs is the input assemblies.

> `-B` (transclose-batch) is the **silent failure** covered in Chapter 5. In this comparison it was verified for all 9 chromosomes × 2 releases **from the parameter records the tools themselves wrote out, not from the command line** — the Chapter 5 lesson of never trusting "I'm sure I passed it", put into practice.

### What to take from this

**1. An assembly is a versioned artifact**

A genome assembly is not settled fact; it is **the best estimate available from the data and algorithms of its moment**. It gets updated after publication, and the distributor can change too (here, Plant GARDEN → MiGD2). **Before starting an analysis, check that the assembly you are about to use is the current version.**

**2. Record what snags you; do not explain it away**

This tutorial deliberately reached no conclusion about CUN's size in §4.4, because the data at hand could not settle it. In r2.0 the gap is gone — **the question was settled by data, in the form of a revised release, not by our reasoning**. An oddity you log rather than over-explain is one you can reconcile later. Seal it with a plausible story and you lose the chance to reconcile at all.

**3. Numbers can change without the conclusion changing**

Possibly the most important point. The reference haplotype's length changed by **16%**, and yet the **ancestry mosaic was 99% identical**.

Neither "the assembly is old, so the conclusion is suspect" nor "the assembly changed and the conclusion held" is the lesson. The lesson is to **measure which of your results moved and which did not**. Direct quantities like size are sensitive; conclusions resting on **relative comparison**, like the mosaic, proved robust. Know which kind yours is.

**4. If you want to compare, vary one thing**

It was only because the pipeline was pinned completely that the difference could be attributed to the assemblies. Change the assemblies and the tool versions at once and nothing can be concluded from either.

### Re-running on v2

The pipeline reads assembly paths from `tables/samplesheet.tsv` and nothing else, so re-running on r2.0 is straightforward.

1. Fetch the six files from the MiGD2 download page:

   ```
   CUNphKu_r2.0.genome.masked.fa.gz    CUNphKi_r2.0.genome.masked.fa.gz
   CKIhap1_r2.0.genome.masked.fa.gz    CKIhap2_r2.0.genome.masked.fa.gz
   CKUhap1_r2.0.genome.masked.fa.gz    CKUhap2_r2.0.genome.masked.fa.gz
   ```

   **Distributor, dataset ID, filename, SHA-256, base counts and N50 for all twelve
   assemblies (r1.0 and r2.0) are collected in
   [`tables/assembly_provenance_v1_v2.tsv`](../../tables/assembly_provenance_v1_v2.tsv).**
   Verify each download against its SHA-256.

2. Point `hap1_path` / `hap2_path` in `tables/samplesheet.tsv` at them (**keep the order `hap1` = CUNphKu, `hap2` = CUNphKi**; see §2.5)
3. Re-run from `qc01_stats.sh` in Chapter 4 through this chapter's analyses

**No script changes are needed.** Two things about r2.0 are worth knowing:

- As the filenames say, it is **repeat-masked** (lower case), whereas r1.0 contained no lower case at all. So that the aligner sees both releases on the same terms, `pg01_prepare_input.sh` **upper-cases the sequence** before graph construction (a no-op on r1.0).
- **Unplaced contigs share the file with the chromosomes.** `pg01_prepare_input.sh` extracts only chr1-9, so this passes through; it warns if the count it extracted is not 9.

Putting the r1.0 and r2.0 results side by side may be the best possible finish to this tutorial.

---

## 7.7 Limits of the graph from this tutorial

**What this tutorial's graph enables**:

- OK: catalog of SVs and SNPs across the three cultivars
- OK: validation of trio phasing
- OK: mastery of pangenome basics (construction, QC, interpretation)
- OK: practical experience with PGGB / odgi / vg

**What it does not enable**:

- Not: population genetics (three cultivars is not a population)
- Not: GWAS (no phenotype data)
- Not: genus-level super-pangenome (no other *Citrus* species)
- Not: high-accuracy phylogenetics (rooted phylogeny, divergence time)

**Constraints from running on the r1.0 assemblies** (§7.6):

- CUN's haplotype sizes and Jaccard values reflect r1.0's uneven path lengths (resolved in r2.0)
- The switches in §7.5.2 may mix genuine crossovers with phase changes on the parental side
- **Which release is closer to biological reality is undetermined.** This tutorial's graph is for learning the methods; do not rest biological conclusions on it
- If you want conclusions, rebuild on r2.0 as well and **compare the two** (§7.6)

---

## 7.8 Next steps

To extend this tutorial's knowledge:

**A. Larger pangenome projects**
- Super-pangenome of 15+ citrus cultivars (assembly from HiFi reads yourself)
- Wider comparison including other Rutaceae (Poncirus, Citropsis, etc.)

**B. Population genetics**
- Allele-level framework with the 18 cultivars from Kiryu et al. 2026
- QTL mapping using the F1 mapping population (99 RAD-Seq individuals, PRJDB15866)

**C. Transcriptome integration**
- Project the 8-tissue Iso-Seq from Isobe 2023 onto the graph
- Allele-specific expression analysis

**D. Breeding applications**
- Marker design for MAS at promising SV loci
- Integration into genotyping pipelines for cultivar evaluation

---

## 7.9 Chapter summary

- Classify SNPs, indels, and SVs from VCF
- Confirm graph structure with 1D and 2D visualization
- Use MAX-based pedigree tests to verify biological validity
- **The assemblies used here have an r2.0, in which the size gap across the six largely disappears** (§7.6)
- Yet **the ancestry mosaic reproduces at 99%** — numbers can change without the conclusion changing
- Windowed similarity alone **cannot separate a genuine crossover from a phase change on the parental side**
- Log the observations that snag you instead of explaining them away — they can be reconciled later
- This tutorial's graph is for foundational learning; applications await follow-up projects

## References

- Wu GA, et al. (2018). Genomics of the origin and evolution of *Citrus*. *Nature* 554:311-316.
- Fujii H, et al. (2016). Parental diagnosis of satsuma mandarin. *Breed Sci* 66:683-691.
- Kiryu Y, et al. (2026). AlleleMiner. *DNA Res* 33:dsag004.
- **MiGD2 (Mikan Genome Database 2, NARO)**: <https://mikan.dna.naro.go.jp/migd2/> (where v2 of the v1 assemblies used in this tutorial is published)

---

[← Chapter 6](06_graph_qc.md) | [back to top](../../README.md)

---

## Tutorial complete 🎉

From Chapter 0 to Chapter 7, we covered pangenome basics, practice, and interpretation:

- **Basics**: what pangenomes are, the workflow, target organisms
- **Practice**: data acquisition, QC, graph construction
- **Interpretation**: quality evaluation, biological reading

You are now equipped to build pangenomes for other organisms with the same approach. We hope you apply this in real research or breeding projects.

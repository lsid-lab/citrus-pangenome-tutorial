# Chapter 7. Graph interpretation

## Learning objectives

- Read variant information from a VCF file
- Visualize a pangenome graph with odgi
- Verify pedigree with real data
- Interpret the results biologically and understand the tutorial's limits

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
#CHROM   POS      ID    REF   ALT      QUAL   FILTER   INFO   FORMAT   satsuma#1  satsuma#2  kishu#1  ...
satsuma#1#chr09   1234  .     A     G        .      .        AT=snp   GT       0          1          0        ...
satsuma#1#chr09   2000  .     ACGT  A        .      .        AT=del3  GT       0          1          1        ...
satsuma#1#chr09   3500  .     C     CGATCG   .      .        AT=ins5  GT       0          0          1        ...
```

Columns:
- **CHROM**: coordinate reference path (here `satsuma#1#chr09`)
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
# ALT variants for satsuma#2 (the other hap from the coordinate anchor)
bcftools view -s satsuma#2 chr09.vcf | \
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
|STS_hap1 vs KSH_hap1|**High** (0.95+ if directly inherited)|
|STS_hap1 vs KSH_hap2|**Moderate** (sister chromosome, 0.4-0.6)|
|STS_hap1 vs KNN_hap*|**Low** (unrelated parent, 0.3-0.5)|

Similarly:
- STS_hap2 should be most similar to one KNN haplotype

### 7.4.2 Specific pedigree pattern

Trace which Kishu haplotype STS_hap1 came from:

```bash
# Jaccard of STS_hap1 vs each KSH hap
awk '$1 ~ /^satsuma#1/ && $2 ~ /^kishu/' chr09_similarity.tsv
```

If the output shows:
```
satsuma#1#chr09  kishu#1#chr09  ...  Jaccard = 0.95
satsuma#1#chr09  kishu#2#chr09  ...  Jaccard = 0.55
```

then we can conclude **STS_hap1 came from KSH_hap1**.

### 7.4.3 Accounting for haplotype leakage

The haplotype leakage discovered in Chapter 4 (Satsuma haps larger than parental) manifests as:

- STS_hap1 vs KSH Jaccard is somewhat lower than expected (e.g., 0.6 instead of 0.7)
- Average-based interpretation may appear to break pedigree consistency

**MAX-based** interpretation confirms pedigree, and the graph itself is correctly built. This was the key point of §6.4.4.

---

## 7.5 Examples of biological insight

Biological patterns that may emerge from the graph:

### 7.5.1 Pummelo-derived regions in Satsuma

According to Wu et al. 2018, Satsuma is mainly mandarin but carries some pummelo ancestry. If in the graph:

- Some path regions in STS do not match KSH or KNN and form independent paths

these may be pummelo-derived. Pummelo-derived regions are concentrated on specific chromosomal segments (Fujii et al. 2016), so visualization may spot them.

### 7.5.2 Recombination in F1

STS_hap1 is "one of KSH's haplotypes"—but actually the result of meiotic recombination. Along one chromosome we may see:

- Regions where STS_hap1 is close to KSH_hap1
- Regions where STS_hap1 is close to KSH_hap2

with switches between them at **recombination breakpoints**.

### 7.5.3 Repeats and TEs

Many large SVs (5-15 kb) are likely **TE insertions**. Kiryu et al. 2026 report a 5.1 kb LTR retrotransposon insertion. Extracting SVs from VCF and BLASTing against TE databases enables this analysis (beyond the scope of this tutorial).

---

## 7.6 Limits of the graph from this tutorial

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

**Impact of haplotype leakage**:

- Satsuma's hap sizes and Jaccard values reflect incomplete trio phasing
- Keep this quirk in mind when interpreting
- Real data is never perfect—this is a real-world lesson

---

## 7.7 Next steps

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

## Chapter summary

- Classify SNPs, indels, and SVs from VCF
- Confirm graph structure with 1D and 2D visualization
- Use MAX-based pedigree tests to verify biological validity
- Recognize haplotype leakage as an interpretive caveat; the graph is still usable
- This tutorial's graph is for foundational learning; applications await follow-up projects

## References

- Wu GA, et al. (2018). Genomics of the origin and evolution of *Citrus*. *Nature* 554:311-316.
- Fujii H, et al. (2016). Parental diagnosis of satsuma mandarin. *Breed Sci* 66:683-691.
- Kiryu Y, et al. (2026). AlleleMiner. *DNA Res* 33:dsag004.

---

[← Chapter 6](06_graph_qc.md) | [back to top](../../README.md)

---

## Tutorial complete 🎉

From Chapter 0 to Chapter 7, we covered pangenome basics, practice, and interpretation:

- **Basics**: what pangenomes are, the workflow, target organisms
- **Practice**: data acquisition, QC, graph construction
- **Interpretation**: quality evaluation, biological reading

You are now equipped to build pangenomes for other organisms with the same approach. We hope you apply this in real research or breeding projects.

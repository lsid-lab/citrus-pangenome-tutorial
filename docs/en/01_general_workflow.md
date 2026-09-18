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
CUN#1#chr01    # Satsuma mandarin, hap1, chromosome 1
CUN#2#chr01    # Satsuma mandarin, hap2, chromosome 1
CKI#1#chr01      # Kishu mandarin, hap1, chromosome 1
```

The `#` character is the delimiter. This scheme lets tools mechanically identify **which cultivar/haplotype a sequence belongs to**.

Most assemblies use different naming (e.g., `CUNphKi_r1.0ch1`), so a **rename step** is needed before pangenome construction. This is handled by `pg01_prepare_input.sh` in this tutorial.

> **Reference**: PanSN-spec was proposed by **Erik Garrison** (also the lead developer of PGGB).
> The formal definition lives at <https://github.com/pangenome/PanSN-spec>.
>
> The exact form in the spec is `[sample_name][delim][haplotype_id][delim][contig_or_scaffold_name]`.
> The delimiter `delim` can in principle be any character, but **`#` is the convention** and is what
> this tutorial uses. `haplotype_id` is numeric; even for unphased or haploid assemblies it is good
> practice to fill it in (e.g., `sample#0#chr01`).
>
> The point of PanSN is to avoid carrying metadata around in a separate file. Because sample and
> haplotype are embedded in the sequence name itself, the information survives FASTA, GFA, VCF, BED
> and GFF alike. That is exactly why **the initial rename is not a step to cut corners on**.

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
L  1  +  3  -  0M    # enter node 3 in reverse complement orientation
```

The columns are `L <From> <FromOrient> <To> <ToOrient> <Overlap>`.

- `<FromOrient>` / `<ToOrient>`: `+` = as given, `-` = reverse complement
- `<Overlap>`: **how many bases the two nodes overlap, written as a CIGAR string** (a required column in GFA1)

#### What is a CIGAR?

A **CIGAR** (Compact Idiosyncratic Gapped Alignment Report) describes **how two sequences correspond to each other**, as a run of "count + operation" pairs. It is the same notation used to describe alignments in SAM/BAM, reused here by GFA.

|Symbol|Meaning|
|---|---|
|`M`|alignment match (a corresponding stretch; may contain mismatches)|
|`=` / `X`|exact match / mismatch (a finer split of `M`)|
|`I`|insertion (extra bases on one side)|
|`D`|deletion (missing bases on one side)|

Example: `10M2D5M` reads as "10 bases correspond → 2-base deletion → 5 bases correspond".

On an L line, the CIGAR says **how many bases the end of the first node overlaps the start of the next**.

- `0M` = **zero-base overlap**, i.e., a blunt join
- `*` = unspecified

**In a pangenome graph GFA, the L-line CIGAR is essentially always `0M`.** Graphs produced by PGGB, Minigraph-Cactus, and minigraph simply cut sequence at node boundaries, so adjacent nodes share no bases. By contrast, the GFA emitted by OLC **assembly** graphs (Canu, miniasm, …) is full of real overlaps such as `120M`. When you are unsure whether a `.gfa` in hand is a pangenome graph or an assembly graph, the last column of the L lines is the quickest tell.

**P (Path)**: describes how a specific sequence traverses the graph
```
P  CUN#1#chr01  1+,2+,4+,5+  *
P  CKI#1#chr01    1+,3+,4+,5+  *
P  CKU#1#chr01  1+,3-,4+,5+  *    # node 3 traversed in reverse = inversion
```

The columns are `P <PathName> <SegmentNames> <Overlaps>`.

- `<SegmentNames>`: a comma-separated list of "node ID + orientation"
  - `3+` = traverse node 3 **forward**
  - `3-` = traverse node 3 in **reverse complement**
- `<Overlaps>`: a list of CIGARs between consecutive nodes. For overlap-free pangenome graphs this is normally just `*`

**A `-` appears exactly where a path runs through an inversion.** In the example above, only Kunenbo carries node 3 flipped, i.e., there is an inversion in that interval. Textbook figures tend to show only `+` orientations like `1+,2+,3+`, but **real graphs always contain paths with `-` in them**. Inversions — and more complex rearrangements such as inverted duplications — exist on the graph *only* as this orientation information. A path is not a list of nodes but a list of **oriented** nodes, and that is a key difference from a plain linear sequence.

Shared nodes, branching edges, and oriented traversal paths together represent the pangenome.

#### Caution: GFA v1.0, GFA v1.1, and rGFA are not the same "GFA"

They all use the `.gfa` extension, but **the dialects differ, and so does downstream tool support**. Getting this wrong produces the classic "the file parses fine but there are zero paths" failure.

|  |Mainly produced by|How paths are represented|
|---|---|---|
|**GFA v1.0**|**PGGB**|`P` lines. The most compatible form; most tools take it as-is|
|**GFA v1.1**|**Minigraph-Cactus (MC)**|**`W` lines (Walk)** instead of `P` lines. Tools that do not know `W` lines silently see no paths|
|**rGFA**|**minigraph**|**No path lines at all.** Position is carried by tags on `S` lines|

**`W` lines (Walk) in GFA v1.1**

```
W  <SampleId>  <HapIndex>  <SeqId>  <SeqStart>  <SeqEnd>  <Walk>

W  CUN  1  chr01  0  30512000  >1>2>4>5
```

Where a `P` line writes `1+,2+` with commas, a `W` line writes `>1>2` (`>` = forward, `<` = reverse; the role is the same as `+` / `-`). Note also that what PanSN packs into the single string `CUN#1#chr01` is **split across three columns** here: sample, haplotype, contig.

**rGFA (reference GFA)**

The rGFA emitted by minigraph is a **strict subset** of GFA:

- Three tags are **mandatory** on every `S` line: `SN:Z:` (name of the stable sequence the segment came from), `SO:i:` (offset on that sequence), and `SR:i:` (rank; `0` = from the linear reference genome, `>0` = non-reference)
- Overlaps between segments are **not allowed**
- There are **no `P` lines and no `W` lines**

In other words, rGFA expresses position as "**which coordinate of the reference does this correspond to**" rather than "which sample traversed which way" — a **reference-centric** design — it may help to think of the decision to build around a reference as being baked into the file format itself. Consequently you cannot directly run path-based odgi analyses or `vg deconstruct` VCF calling on an rGFA.

**Practical note: MC graphs may need a downgrade**

PGGB emits GFA v1.0 (`P` lines), so downstream processing generally just works. **MC, on the other hand, emits GFA v1.1 (`W` lines)**, so before handing it to a tool that does not understand `W` lines you need to convert (downgrade) the walks to `P` lines:

```bash
# Convert W lines → P lines (GFA v1.1 → v1.0-compatible)
vg convert -g graph.mc.gfa -f -W > graph.p-lines.gfa
#   -g : input is GFA
#   -f : output GFA
#   -W : do not use W lines; write all paths as P lines (--no-wline)
```

This tutorial uses PGGB, so the conversion never comes up here — but **you will certainly hit it the first time you work with a published MC graph** (e.g., the HPRC human pangenome). When a tool tells you it found no paths, start by comparing:

```bash
grep -c '^P' graph.gfa   # number of P lines
grep -c '^W' graph.gfa   # number of W lines
```

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

#### BUSCO — how "completeness" is measured

**BUSCO** (Benchmarking Universal Single-Copy Orthologs, <https://busco.ezlab.org/>) is the de-facto standard tool for assessing assembly **completeness**.

The idea is simple: take a curated set of genes that **nearly every species in a given lineage carries in exactly one copy**, then count how many of them can be found in your assembly. Genes that are missing indicate that the assembly is incomplete (or broken). Where N50 measures **how contiguous** an assembly is, BUSCO measures **whether the content is all there** — you need both before you can talk about assembly quality.

Results are reported as percentages in four categories:

|Category|Meaning|How to read it|
|---|---|---|
|**S** (Complete, single-copy)|found complete, in exactly one copy|the higher the better|
|**D** (Complete, duplicated)|found complete, but in two or more copies|a high value suggests **haplotype duplication (leakage)** or spurious duplication|
|**F** (Fragmented)|found only partially|indicator of fragmentation / assembly error|
|**M** (Missing)|not found at all|indicator of missing sequence|

**C (Complete) = S + D**, summarized in one line such as `C:97.2%[S:94.1%,D:3.1%],F:1.0%,M:1.8%`. For haplotype-resolved assemblies like the ones in this tutorial, run BUSCO **separately per haplotype**. **A high D on just one haplotype** suggests that haplotype carries duplicated copies of the same sequence (haplotype leakage).

**Choosing the lineage dataset is the part that matters most**

BUSCO numbers depend on **which lineage dataset you used**. Three rules:

1. **Pick the most specific (deepest) lineage that still contains your organism.**
   The higher up the tree you go, the fewer genes are shared by every member, so the test uses fewer genes and becomes **more permissive and less informative**. For instance, `eukaryota_odb10` contains only ~255 genes, whereas `eudicots_odb10` has ~2,326 — a far stricter and higher-resolution assessment.

2. **Always check what is actually available to you.**
   ```bash
   busco --list-datasets
   ```

3. **If you are unsure of the lineage, let BUSCO decide.**
   ```bash
   busco -i assembly.fasta -m genome --auto-lineage-euk -o out
   ```
   Auto-detection is slow, though, so specify `-l` explicitly when you already know the lineage.

For *Citrus* (family Rutaceae), the plant lineage datasets get more specific in this order:

```
eukaryota  →  viridiplantae  →  embryophyta  →  eudicots
  (~255)        (~425)            (~1,614)        (~2,326 genes)
```

There is **no dataset specific to Rutaceae or Sapindales**, so the right choice for citrus is the deepest available one, **`eudicots_odb10`**. Order-level datasets such as `brassicales` or `fabales` do exist, but **citrus does not belong to those orders, so they must not be used**. "More specific is better" only holds as long as the dataset actually **contains your organism**.

```bash
busco -i data/raw/CUN/CUNphKi_r1.0.pmol.fasta -m genome \
      -l eudicots_odb10 -c 16 -o busco/CUN_hap2
```

> The `odb10` suffix is the OrthoDB version the dataset was built from. Newer BUSCO releases also ship `odb12` datasets, so check `busco --list-datasets` to see what your installation has.
>
> **BUSCO scores from different datasets or different BUSCO versions are not comparable.** Whenever you compare your numbers against a paper or another assembly, make sure the lineage dataset name and the BUSCO version match.

The actual commands (using `compleasm`, a faster BUSCO reimplementation) are covered in §4.6.

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

### The exception: when you cannot split by chromosome

Splitting by chromosome only works because of an assumption: that **sequence has not moved between chromosomes**. Where that assumption breaks, the graph must be built over **the whole genome at once**.

The canonical example is **cancer genomes**. Tumour cells routinely undergo **interchromosomal rearrangement** — translocations, chromothripsis, and the like. If you cut the data up by chromosome, the junctions — **precisely the structural changes you wanted to see** — drop out of the graph entirely. The same applies to cross-species comparisons involving chromosome fusion/fission, and to allopolyploids where recombination between subgenomes is the object of study.

However, **whole-genome construction demands enormous compute resources and runtime**. Where a per-chromosome run lets you stream 9 jobs in parallel at tens of GB each, a single whole-genome run performs all-vs-all alignment with every sequence resident at once: expect **hundreds of GB to terabytes of RAM and days to weeks of wall-clock time** — with no parallelism to fall back on, and a full restart if it fails partway.

So "per chromosome or all at once" is a trade-off between **whether you need to see interchromosomal rearrangements** and **whether you can afford the cost**. For the citrus cultivar comparison in this tutorial, chromosome correspondence holds, so per-chromosome construction is fine.

---

## Chapter summary

- Pangenome construction follows a general 6-step workflow
- Common concepts across tools: **PanSN**, **GFA**, **node/edge/path**, **compression ratio**
- **GFA has dialects**: PGGB = v1.0 (`P` lines), MC = v1.1 (`W` lines), minigraph = rGFA (no path lines). MC graphs may need `vg convert -f -W` to downgrade
- A path is a list of **oriented** nodes; a `-` means an inversion
- For completeness (BUSCO), pick the **most specific lineage dataset that still contains your organism** (`eudicots_odb10` for citrus)
- Chromosome-by-chromosome construction is standard for large genomes — but targets with **interchromosomal rearrangement (e.g., cancer genomes) require whole-genome runs**, at a dramatically higher cost

The next chapter introduces the three target cultivars (Satsuma, Kishu, Kunenbo) used in this tutorial.

## References

- Garrison E, et al. (2018). Variation graph toolkit improves read mapping by representing genetic variation in the reference. *Nat Biotechnol* 36:875-879.
- Li H (2020). The design and construction of reference pangenome graphs with minigraph. *Genome Biol* 21:265.
- Hickey G, et al. (2024). Pangenome graph construction from genome alignments with Minigraph-Cactus. *Nat Biotechnol* 42:663-673.
- Heumos S, et al. (2024). Pangenome graph layout by Path-Guided SGD. *Bioinformatics*.
- Manni M, et al. (2021). BUSCO Update: novel and streamlined workflows. *Mol Biol Evol* 38:4647-4654.

### Specifications and online resources

- **PanSN-spec** (Erik Garrison): <https://github.com/pangenome/PanSN-spec>
- **GFA specification (v1.0 / v1.1)**: <https://github.com/GFA-spec/GFA-spec>
- **rGFA specification** (Heng Li, gfatools): <https://github.com/lh3/gfatools/blob/master/doc/rGFA.md>
- **BUSCO** (lineage dataset list and user guide): <https://busco.ezlab.org/>
- **vg convert** (GFA interconversion): <https://github.com/vgteam/vg/wiki>

---

[← Chapter 0](00_what_is_pangenome.md) | [Chapter 2: Target organisms →](02_target_organisms.md)

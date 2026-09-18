# Chapter 6. Graph QC

## Learning objectives

- Understand the 4-level QC framework for pangenome graphs
- Learn what to check at each level, and why
- Know specific quality thresholds
- Recognize failure patterns and how to address them

---

## 6.1 Why four levels?

Judging whether a pangenome graph is usable requires multiple perspectives:

|Level|Perspective|Metrics|What passing means|
|---|---|---|---|
|Level 1|**Structure**|node/edge/path counts, compression ratio|Graph is technically sound|
|Level 2|**Preservation**|Input vs graph path lengths|No data loss|
|Level 3|**Biological validity**|Path similarity vs pedigree|Biologically sensible|
|Level 4|**SV signal**|VCF SV counts and distribution|Usable for downstream analysis|

**There is no point in checking lower levels while upper levels fail.** If Level 1 fails, fix the structure before counting SVs.

---

## 6.2 Level 1: Structural check

### 6.2.1 What to check

Run `odgi stats -S` on the graph:

```bash
odgi stats -i chr09.fa.gz.<hash>.smooth.final.og -S
```

Output (tab-separated, 5 columns):
```
#length     nodes    edges    paths    steps
49786530    1381521  1874185  6        6020593
```

Column meanings:
- **length**: total base count across all nodes = **graph_bp**
- **nodes**: number of nodes
- **edges**: number of edges
- **paths**: path count (= input haploid count, expected 6)
- **steps**: total path traversal steps

### 6.2.2 Primary metric: compression ratio

The most important indicator:

```
compression ratio = graph_bp / input_bp
```

|Ratio|Verdict|Meaning|
|---|---|---|
|0.19-0.30|OK|Normal pangenome for closely related haploids|
|0.30-0.40|WARN|Slightly divergent, worth checking|
|0.40-0.60|WARN|Parameter tuning needed|
|**> 0.90**|**FAIL**|Possible `-B` issue (§6.2.3)|

### 6.2.3 Check `transclose-batch` in `.params.yml`

Verify the version-dependent `-B` behavior mentioned in §5.5.3:

```bash
grep transclose-batch results/pggb/chr09/*.params.yml
```

Expected:
- `transclose-batch: 1000000000` → OK: `-B 1G` was applied correctly
- `transclose-batch: 1000000` → WARN: **default parse may have failed**, rerun with explicit `-B 1G`

Even when `.err` is empty and exit code is 0, this one line reveals whether the behavior was as intended. Future PGGB versions may fix the underlying parse issue, but this check is a good habit.

### 6.2.4 Other Level 1 indicators

**path count must equal 6**:
- 3 cultivars × 2 hap = 6
- If not 6, check input PanSN naming and rerun pg01

**inflation ratio** (auxiliary):
```
inflation = graph_bp × n_haplotypes / input_bp = compression ratio × 6
```

- Expected: 1.0-1.6
- For compression 0.25, inflation ≈ 1.5

**Reference values from a collaborator (all 9 chromosomes)**:

|chr|Compression|Node count|
|---|---|---|
|ch1|0.2391|1,074,276|
|ch2|0.2484|1,867,874|
|ch3|0.2604|2,882,578|
|ch4|**0.1862**|706,073|
|ch5|0.2638|1,667,106|
|ch6|0.2443|1,141,586|
|ch7|0.2743|1,901,707|
|ch8|**0.2778**|1,483,201|
|ch9|0.2541|1,383,141|

All 9 chromosomes fall in the range 0.186 - 0.278 and pass the threshold ≤ 0.35.

---

## 6.3 Level 2: Path preservation

### 6.3.1 What to verify

Check that the six input haploids are preserved as **complete paths** in the graph. If any part of a path is lost, the graph cannot reproduce the input.

### 6.3.2 Method

```bash
# Extract paths from the graph as FASTA
odgi paths -i chr09.smooth.final.og -f > chr09_paths.fa

# Get path lengths
seqkit fx2tab -nl chr09_paths.fa
```

Example output:
```
CUN#1#chr09    45,179,797
CUN#2#chr09    29,817,104
CKI#1#chr09      28,530,969
CKI#2#chr09      31,864,882
CKU#1#chr09    30,157,514
CKU#2#chr09    30,337,761
```

Compare to **input FASTA sequence lengths**.

### 6.3.3 Thresholds

|Difference|Verdict|
|---|---|
|0% (exact)|OK: ideal|
|< 0.1%|OK|
|0.1-1%|WARN|
|**> 1%**|**FAIL**: path was lost, critical|

### 6.3.4 In practice

`pg03_qc_graph.sh` automates this:

```bash
bash scripts/pg03_qc_graph.sh . chr09
```

Expected output:
```
Preservation: YES_max_diff=0.000%_at_
```

max_diff = 0.000% means input is fully preserved.

---

## 6.4 Level 3: Pedigree consistency

### 6.4.1 Verifying the family with the graph

Because the pedigree is known (CUN = CKI × CKU), we can check whether path similarities in the graph match the expected relationship.

Expected:
- `CUN#1` (CUNphKu, Kunenbo-derived) should be highly similar to one Kunenbo haplotype
- `CUN#2` (CUNphKi, Kishu-derived) should be highly similar to one Kishu haplotype

### 6.4.2 Running odgi similarity

```bash
odgi similarity -i chr09.smooth.final.og > chr09_similarity.tsv
```

Output columns:
```
group.a  group.b  a.length  b.length  intersection  jaccard  cosine  dice  identity
```

We focus on the `jaccard.similarity` column.

### 6.4.3 AVG-based vs MAX-based tests

**Averages (AVG) can mislead**:
- "`CUN#1` vs Kunenbo (average)" averages against both CKU_hap1 and CKU_hap2
- `CUN#1` is a recombinant **mosaic** of CKU_hap1 and CKU_hap2 (§2.4), so which one is closer switches from segment to segment. Averaging the two flattens that structure, **shrinking the gap against CKI and burying the pedigree signal**

**Maxima (MAX) are more meaningful**:
- MAX(`CUN#1` vs CKU_hap1, `CUN#1` vs CKU_hap2) = the closer of the two CKU haplotypes
- Even with recombination, every segment of `CUN#1` descends from one of CKU's two chromosomes, so this MAX is necessarily higher than the MAX on the CKI side

> **Note**: this tests that **the donor parent is CKU**, not that "the donor haplotype is CKU_hap1". A single whole-chromosome Jaccard cannot identify *which* parental haplotype a segment came from (§2.4); that needs the windowed analysis in §7.5.2.

`pg03_qc_graph.sh` computes both:

```
CUN#1 (CUNphKu, Kunenbo-derived) vs each parent haplotype:
  CKU#1#chr09    Jaccard = 0.62  ← closer CKU hap
  CKU#2#chr09    Jaccard = 0.48
  CKI#1#chr09    Jaccard = 0.42
  CKI#2#chr09    Jaccard = 0.35

--- Pedigree tests ---
  AVG-based: NO   (perturbed by uneven path lengths)
  MAX-based: YES  (CKU MAX 0.62 > CKI MAX 0.42)
```

### 6.4.4 Effect of CUN's longer paths

The observation from Chapter 4 that **CUN's two haplotypes are larger than the parental lines** (§4.4) affects Jaccard similarity.

Because Jaccard = intersection / union, **when one path is longer, union grows and Jaccard drops**.

If `CUN#1` is 45 Mb (1.5x the expected 30 Mb) and other haps are 30 Mb:
- intersection = 25 Mb (shared)
- union = 45 + 30 - 25 = 50 Mb
- Jaccard = 0.50

Numerically low, but only because one path is longer—not because biological similarity is lower. **MAX-based is more appropriate** for pedigree tests.

### 6.4.5 Thresholds

MAX-based:
- `CUN#1` has highest similarity to some CKU hap → OK
- `CUN#2` has highest similarity to some CKI hap → OK

If both, pedigree is **consistent**. If either fails, suspect hap1/hap2 label swap or a trio phasing failure.

---

## 6.5 Level 4: SV signal

### 6.5.1 VCF generation

```bash
vg deconstruct -P CUN#1 -H '#' -a -e -t 8 chr09.smooth.final.gfa > chr09.vcf
```

- `-P CUN#1`: coordinate anchor (reference path)
- `-H '#'`: PanSN separator
- `-a`: include nested variants
- `-e`: path-traversal-based calls

### 6.5.2 Statistics

```bash
# Total variants
grep -vc '^#' chr09.vcf

# SNPs
grep -v '^#' chr09.vcf | awk 'length($4)==1 && length($5)==1' | wc -l

# SVs (>= 50 bp)
grep -v '^#' chr09.vcf | awk 'length($4)>=50 || length($5)>=50' | wc -l
```

### 6.5.3 Expected values (chr09, ~30 Mb)

|Variant type|Expected count|
|---|---|
|SNPs|10,000-30,000|
|Small indels|5,000-15,000|
|SVs (≥50bp)|100-500|

**Collaborator's actual measurement**: chr09 has 432,698 total variants (large because multi-allelic calls are expanded).

### 6.5.4 Size distribution

Is the SV size distribution biologically plausible?

- 1-10 bp: 60-70% (small indels dominate)
- 11-100 bp: 20-30%
- 100-1000 bp: 5-10%
- 1-10 kb: 3-5%
- 10 kb+: 1-2% (TE insertions etc.)

If distributions are strongly skewed toward one size, suspect parameter issues.

---

## 6.6 Overall verdict

### 6.6.1 Verdict rules

`pg03_qc_graph.sh` combines results:

|Verdict|Condition|
|---|---|
|**PASS**|All levels OK|
|**FAIL**|Level 1 or Level 2 FAILs|
|**WARN**|Only Level 3 or Level 4 has issues (graph usable, but interpret carefully)|

### 6.6.2 Why Level 3/4 issues are WARN, not FAIL

- Level 1/2 failures = technical problems with graph structure (unusable)
- Level 3/4 failures = questions about interpretation or downstream analysis (graph is still usable)

Even when Level 3 pedigree fails, **the graph itself is technically complete**. Data-side issues like haplotype leakage or phasing quirks are frequent culprits.

---

## 6.7 Failure patterns and remedies

|Symptom|Cause|Remedy|
|---|---|---|
|Compression > 0.9|`-B` default parse issue|Specify `-B 1G` explicitly and rerun|
|path count ≠ 6|Input naming|Rerun pg01, verify PanSN naming|
|Path diff > 1%|Graph construction failed|Re-tune parameters|
|SV count = 0|VCF generation failure|Check vg deconstruct options|
|SV count > 100,000|Multi-allelic expansion|Drop `-e` to confirm|
|Pedigree AVG-NO, MAX-YES|**Haplotype leakage** or Jaccard length sensitivity|OK — interpretation caveat only|
|Pedigree MAX-NO|hap1/hap2 label swap?|Check data source|

---

## Chapter summary

- Four QC levels: structure, preservation, biological validity, SV signal
- Level 1 compression is the primary indicator; also check `.params.yml` transclose-batch
- Aim for compression ≤ 0.35
- MAX-based is more meaningful than AVG-based for pedigree tests
- Level 3/4 failures are WARN; the graph is still usable

The next chapter walks through **interpretation, visualization, and use** of the graph.

## References

- odgi documentation: https://pangenome.github.io/odgi.github.io/
- vg deconstruct: https://github.com/vgteam/vg/wiki/vg-deconstruct

---

[← Chapter 5](05_graph_construction.md) | [Chapter 7: Graph interpretation →](07_graph_interpretation.md)

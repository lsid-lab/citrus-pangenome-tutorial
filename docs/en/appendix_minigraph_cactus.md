# Appendix: the same data with Minigraph-Cactus

> **This appendix is untested.** It has not been run in this tutorial's environment. The commands adapt the procedure from the Cactus documentation to this tutorial's data, and **have not been verified to work as written**. Check the [official documentation](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md) for the current syntax before running anything.

This is the procedure for actually obtaining the Minigraph-Cactus (MC) benefits listed in §5.3, running the same six haploids through MC.

---

## A.1 What MC gives you

Restating §5.3:

1. **The coordinate system is fixed from the start** — MC requires a reference genome; that reference is never clipped or cyclic and directly defines the VCF coordinate system and the chromosome decomposition
2. **Read-mapping artifacts come out** — GBZ plus `.dist` / `.min` / `.hapl` indexes, usable directly with vg Giraffe
3. **Softmasked input does not matter** — no preprocessing needed even though r2.0 is repeat-masked
4. **It is faster** — designed for closely related samples of the same species

---

## A.2 Preparing the input

MC takes not PanSN-named FASTA but a two-column text file, the **seqFile**.

```
# seqfile.txt
# <sample>.<haplotype>   <path to FASTA>
CUN.1   data/raw/CUN/CUNphKu_r1.0.ch1-9.fasta.gz
CUN.2   data/raw/CUN/CUNphKi_r1.0.pmol.fasta.gz
CKI.1   data/raw/CKI/CKIhap1_r1.0.pmol.fasta.gz
CKI.2   data/raw/CKI/CKIhap2_r1.0.pmol.fasta.gz
CKU.1   data/raw/CKU/CKUhap1_r1.0.pmol.fasta.gz
CKU.2   data/raw/CKU/CKUhap2_r1.0.pmol.fasta.gz
```

MC converts the `sample.haplotype` dot notation into PanSN (`CUN#1#…`) internally. Note that this **preserves the tutorial's PanSN assignment** (`CUN#1` = CUNphKu, `CUN#2` = CUNphKi; §2.5).

**Differences from PGGB:**

- **No need to split by chromosome.** MC decomposes chromosomes itself, based on the reference. There is no equivalent of `pg01_prepare_input.sh`
- **No upper-casing either** (§7.6)
- So the files in **`data/raw/` can be handed over as they are**

### Choosing the reference

This is the awkward part for this dataset. MC's `--reference` **must be haploid**. Human work has GRCh38 or CHM13 as a separate haploid assembly; here there are only three cultivars × two haplotypes.

So one of the six has to be **designated** as the reference. To stay consistent with Chapter 5, that is `CUN.1` (= `CUN#1`, CUNphKu).

**This choice means something different from PGGB's `-V`.** With PGGB only the coordinate system changed and the graph was the same; **with MC the reference shapes the graph itself**, because it is constrained never to be clipped or cyclic. That difference is exactly what §5.2 means by PGGB being reference-free.

---

## A.3 Running it

```bash
# Example via Docker
docker run --rm -it \
  -u "$(id -u):$(id -g)" \
  -v "$PWD":/data -w /data \
  quay.io/comparative-genomics-toolkit/cactus:latest \
  cactus-pangenome ./js ./seqfile.txt \
    --outDir results/mc \
    --outName citrus-pg \
    --reference CUN.1 \
    --vcf --giraffe --gfa --gbz
```

| Option | Meaning |
|---|---|
| `./js` | Toil job store (a working directory; it must not already exist, so delete or rename it between runs) |
| `./seqfile.txt` | The two-column file above |
| `--reference CUN.1` | The reference. **Required**, and it defines the VCF coordinate system and chromosome decomposition |
| `--vcf` | Emit a VCF |
| `--giraffe` | Emit vg Giraffe indexes (`.dist` / `.min` / `.hapl`) |
| `--gfa` | Emit GFA |
| `--gbz` | Emit GBZ (compressed graph) |

Resources scale with the input. HPRC-scale work is documented as wanting 64 cores, 512 GB RAM and 3 TB of disk. Six citrus haploids are far smaller, so this tutorial's assumption of 32 CPU / 128 GB should suffice — but that is **unverified**.

---

## A.4 Feeding the output into this tutorial's analyses

### Mind the GFA version

**MC emits GFA v1.1, representing paths as `W` lines (Walks) rather than `P` lines** — exactly the case covered in §1.2.2.

Convert to `P` lines before handing the graph to any tool that does not read `W` lines:

```bash
vg convert -g results/mc/citrus-pg.gfa -f -W > results/mc/citrus-pg.p-lines.gfa
#   -g : input is GFA
#   -f : output GFA
#   -W : do not use W lines; write all paths as P lines
```

**The downgrade described back in Chapter 1 is genuinely required here.** Run it before passing the graph to the `odgi`-based QC in Chapter 6.

### Path names

MC's path names are PanSN as well, so they look like `CUN#1#chr09`. **The Chapter 6 and 7 scripts (`pg03_qc_graph.sh` and friends) match on the path-name prefix, so they should work unchanged.** MC does treat the reference path specially, though, so check the path list before handing the graph to `odgi`:

```bash
odgi paths -i results/mc/citrus-pg.og -L
```

### VCF

The VCF from `--vcf` is in `CUN.1` coordinates. Chapter 7 assumes a VCF in `CUN#1` coordinates, so **the coordinate systems agree**. How variants are represented (nested variants in particular) differs from PGGB's `vg deconstruct` output, so **do not compare the counts directly**.

---

## A.5 Comparing against the PGGB result

If you run both, the same comparison as §7.6 becomes available:

- **Compression ratio and node count** (§6.3)
- **Path preservation** (§6.3)
- **Pedigree consistency** (§6.4) — does the MAX-based test pass for both?
- **Agreement of the ancestry mosaic** (§7.5.2) — do the per-window ancestry calls agree?

§7.6 found that **direct quantities like size move, while conclusions resting on relative comparison hold**. Whether that also survives a change of tool is beyond this tutorial's scope, but it is a comparison worth making.

---

## References

- Hickey G, et al. (2024). Pangenome graph construction from genome alignments with Minigraph-Cactus. *Nat Biotechnol* 42:663-673.
- **Cactus pangenome pipeline documentation**: <https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md>

---

[← Back to Chapter 5](05_graph_construction.md)

# 付録: Minigraph-Cactus で同じデータを扱う

> **この付録は未検証です。** 本教材の環境では実行していません。コマンドは Cactus 公式ドキュメントの手順を本教材のデータに当てはめたもので、**そのまま動くことを確認していません**。実行する際は必ず[公式ドキュメント](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md)で最新の書式を確認してください。

§5.3 で挙げた Minigraph-Cactus (MC) のメリットを実際に得たい場合の手順です。同じ 6 haploid を MC で扱う流れを示します。

---

## A.1 MC で得られるもの

§5.3 の再掲です。

1. **参照座標系が最初から固定される** —— MC は参照ゲノムの指定が必須で、その参照は clip も cycle もされず、VCF の座標系と染色体分割をそのまま規定します
2. **リードマッピング用の成果物が出る** —— GBZ と `.dist` / `.min` / `.hapl` インデックスが出力され、vg Giraffe にそのまま使えます
3. **入力の softmask を気にしなくてよい** —— r2.0 が repeat-mask 済みでも前処理が要りません
4. **速い** —— 同種内の近縁サンプル向けの設計です

---

## A.2 入力の準備

MC は PanSN 形式の FASTA ではなく、**seqFile** という 2 列のテキストファイルで入力を受け取ります。

```
# seqfile.txt
# <サンプル名>.<ハプロタイプ番号>   <FASTA へのパス>
CUN.1   data/raw/CUN/CUNphKu_r1.0.ch1-9.fasta.gz
CUN.2   data/raw/CUN/CUNphKi_r1.0.pmol.fasta.gz
CKI.1   data/raw/CKI/CKIhap1_r1.0.pmol.fasta.gz
CKI.2   data/raw/CKI/CKIhap2_r1.0.pmol.fasta.gz
CKU.1   data/raw/CKU/CKUhap1_r1.0.pmol.fasta.gz
CKU.2   data/raw/CKU/CKUhap2_r1.0.pmol.fasta.gz
```

`サンプル名.ハプロタイプ番号` というドット記法を、MC が内部で PanSN (`CUN#1#…`) に変換します。**本教材の PanSN 割り当て(`CUN#1` = CUNphKu、`CUN#2` = CUNphKi。§2.5)をそのまま維持している**ことに注意してください。

**PGGB との違い**:

- **染色体別に分ける必要がありません。** MC は参照に基づいて自分で染色体を分割します。`pg01_prepare_input.sh` に相当する前処理は不要です
- **大文字化も不要**です(§7.6)
- 一方で、**`data/raw/` のファイルをそのまま渡せます**

### 参照をどう決めるか

ここが本教材のデータで悩ましい点です。MC の `--reference` は **haploid である必要があります**。ヒトの HPRC なら GRCh38 や CHM13 という独立した haploid アセンブリがありますが、本教材には「3 品種 × 2 ハプロタイプ」しかありません。

したがって **6 本のうち 1 本を参照として指名する**ことになります。第5章と揃えるなら `CUN.1`(= `CUN#1`、CUNphKu)です。

**この選択は PGGB の `-V` とは意味が違います。** PGGB では座標系だけが変わり graph は同じでしたが、**MC では参照が graph の構造そのものに影響します**(参照は clip されず cycle も持たない、という制約が入るため)。§5.2 で「PGGB は reference-free」と書いたのは、まさにこの違いのことです。

---

## A.3 実行

```bash
# Docker で実行する例
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

| オプション | 意味 |
|---|---|
| `./js` | Toil の job store(作業用ディレクトリ。既存だとエラーになるので毎回消すか別名にする) |
| `./seqfile.txt` | 上で作った 2 列のファイル |
| `--reference CUN.1` | 参照。**必須**で、VCF の座標系と染色体分割を規定する |
| `--vcf` | VCF を出力 |
| `--giraffe` | vg Giraffe 用のインデックス (`.dist` / `.min` / `.hapl`) を出力 |
| `--gfa` | GFA を出力 |
| `--gbz` | GBZ(圧縮グラフ)を出力 |

リソースは対象規模によります。HPRC 規模では 64 コア / 512 GB RAM / 3 TB ディスクが推奨とされています。柑橘 6 haploid はそれよりはるかに小さいので、本教材の前提(32 CPU / 128 GB)でも動く見込みですが、**未検証**です。

---

## A.4 出力を本教材の解析に繋ぐ

### GFA のバージョンに注意

**MC が出す GFA は v1.1 で、path を `P` 行ではなく `W` 行 (Walk) で表現します。** 第1章 §1.2.2 で扱ったとおりです。

`W` 行に未対応のツールに渡す前に、`P` 行へ変換してください。

```bash
vg convert -g results/mc/citrus-pg.gfa -f -W > results/mc/citrus-pg.p-lines.gfa
#   -g : 入力が GFA
#   -f : GFA で出力
#   -W : W 行を使わず、すべての path を P 行として書き出す
```

**第1章で説明した downgrade が、ここで実際に必要になります。** 第6章の `odgi` を使った QC に渡すなら、この変換を先に通してください。

### path 名の対応

MC が付ける path 名も PanSN なので、`CUN#1#chr09` のような形になります。**第6・7章の解析スクリプト(`pg03_qc_graph.sh` など)は path 名の接頭辞で判別しているので、そのまま使えるはず**です。ただし MC は参照 path を特別扱いするため、`odgi` に渡す前に path の一覧を確認してください。

```bash
odgi paths -i results/mc/citrus-pg.og -L
```

### VCF

`--vcf` で出る VCF は `CUN.1` 基準です。第7章の解析は `CUN#1` 基準の VCF を前提にしているので、**座標系としては一致します**。ただし変異の表現(nested variant の扱いなど)は PGGB の `vg deconstruct` 出力とは異なるため、**数値をそのまま比較しないでください**。

---

## A.5 PGGB の結果と比べる

もし両方を走らせたなら、第7章 §7.6 でやったのと同じ形の比較ができます。

- **圧縮率とノード数**(§6.3)
- **path の保存性**(§6.3)
- **pedigree の整合性**(§6.4)—— MAX-based の検定が両方で通るか
- **来歴モザイクの一致**(§7.5.2)—— 窓ごとの由来判定が一致するか

§7.6 で見たとおり、**「サイズのような直接的な量は動くが、相対比較に立つ結論は頑健」**という傾向が、ツールを変えた場合にも当てはまるかどうか。これは教材の範囲を超えますが、やってみる価値のある比較です。

---

## 参考文献

- Hickey G, et al. (2024). Pangenome graph construction from genome alignments with Minigraph-Cactus. *Nat Biotechnol* 42:663-673.
- **Cactus pangenome pipeline のドキュメント**: <https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md>

---

[← 第5章に戻る](05_graph_construction.md)

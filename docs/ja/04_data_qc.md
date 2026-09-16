# 第4章 データの品質確認

## 学ぶこと

- なぜ pangenome 構築前に QC が必要か
- `seqkit stats` の基本的な使い方
- 柑橘における各指標の期待値
- 品質異常のパターンと対処

---

## 4.1 なぜ QC が必要か

Pangenome graph は複数のアセンブリを重ね合わせて作ります。もし 1 つでも品質の悪いアセンブリが混ざっていると、実際には存在しない構造変異(SV)が「発見」されてしまいます。

バイオインフォマティクスの鉄則:
> **Garbage in, garbage out**(ゴミを入れればゴミが出る)

Pangenome 構築の前に必ず QC を通します。

---

## 4.2 ゲノムアセンブリの基本用語(復習)

### Contig と Scaffold
- **Contig**: リードだけで繋げられた連続した配列
- **Scaffold**: contig 同士を Hi-C などの補助情報で染色体上に並べたもの

理想は「染色体1本 = contig 1つ = scaffold 1つ」。柑橘 (n=9) なら 9 本の contig が理想。

### N50
アセンブリの「断片化度」を表す最も重要な指標。

**定義**: 全配列を長い順に並べ、上から足していって**総長の 50% を超えた時点での配列の長さ**。

- 染色体スケール = N50 が数十 Mb
- Contig レベル = N50 が数百 kb 〜 数 Mb

### Haplotype-resolved
二倍体では、父由来と母由来の染色体を**別々にアセンブル**したもの。第2章 §2.5 で説明した通り、Isobe 2023 のアセンブリはこの形式です。

---

## 4.3 `seqkit stats` による基本統計

`seqkit stats` は、FASTA ファイルの基本統計を 1 コマンドで出してくれるツールです。

### インストール

```bash
# conda
mamba install -c bioconda seqkit

# または単一バイナリ
wget https://github.com/shenwei356/seqkit/releases/download/v2.10.0/seqkit_linux_amd64.tar.gz
tar xzf seqkit_linux_amd64.tar.gz && mv seqkit ~/bin/
```

### 実行

```bash
seqkit stats -a data/satsuma/CUNphKi_r1.0.pmol.fasta.gz
```

`-a` (all) オプションで N50 などの詳細も出ます。

### 出力例(温州 hap1 の期待値)

```
file                              format  type  num_seqs  sum_len      min_len    avg_len    max_len    N50
CUNphKi_r1.0.pmol.fasta.gz        FASTA   DNA          9  348,509,231  15,234,891 38,723,247 48,123,456  38,241,932
```

### 各列の意味と柑橘での期待値

| 列 | 意味 | 期待値 (haploid citrus) |
|---|---|---|
| `num_seqs` | 配列(contig/chromosome)の数 | **9** (n=9染色体) |
| `sum_len` | 総塩基数 | **300-360 Mb** |
| `min_len` | 最短配列長 | ≥ 10 Mb (最小染色体) |
| `avg_len` | 平均長 | 35 Mb 前後 |
| `max_len` | 最長配列長 | ≤ 50 Mb (最大染色体) |
| `N50` | 中央値的な指標 | 大きいほど良い、20+ Mb |
| `GC(%)` | GC 含量 | 34-38% |

---

## 4.4 6 haploid の一括 QC

`scripts/qc01_stats.sh` を使うと、6 haploid すべての統計を 1 コマンドで取得できます:

```bash
bash scripts/qc01_stats.sh tables/samplesheet.tsv .
```

出力: `qc/stats_summary.tsv`

期待される結果の例:

| sample | hap | num_seqs | sum_len_Mb | N50_Mb | GC% | 判定 |
|---|---|---|---|---|---|---|
| STS | 1 | 9 | 348.5 | 34.5 | 35.9 | ⚠️(サイズやや大) |
| STS | 2 | 9 | 358.0 | 44.2 | 35.0 | ⚠️(サイズやや大) |
| KSH | 1 | 9 | 304.2 | 33.5 | 36.0 | ✅ |
| KSH | 2 | 9 | 310.3 | 32.5 | 36.0 | ✅ |
| KNN | 1 | 9 | 323.9 | 36.6 | 35.9 | ✅ |
| KNN | 2 | 9 | 303.4 | 32.9 | 36.0 | ✅ |

---

## 4.5 見るべきポイント

### ポイント1: 染色体スケール判定

- `num_seqs ≤ 20`: 染色体スケール ✅
- `num_seqs 20-200`: 半断片化(scaffold が繋がりきっていない)⚠️
- `num_seqs > 200`: 過度な断片化、pangenome 構築が難しい ❌

### ポイント2: サイズ

- 300-360 Mb (haploid): 正常範囲
- < 290 Mb: 情報欠損の疑い
- > 370 Mb: collapse 不完全、haplotype 混入の疑い

### ポイント3: hap1 と hap2 の対比

同一品種の 2 haplotype はサイズが揃うべき:

- 差 ≤ 5%: 完全に対応 ✅
- 差 5-15%: 微妙、要注意 ⚠️
- 差 > 15%: **haplotype 分離失敗の疑い** ❌

### ポイント4: GC 含量

- 34-38%: 柑橘の正常範囲 ✅
- < 34% or > 38%: コンタミネーションを疑う

---

## 4.6 実データで見つかる異常:Haplotype leakage

上の期待値と実データを比較すると、**温州みかんの 2 haplotype がやや大きい**ことに気づきます:

- STS hap1: 348.5 Mb
- STS hap2: 358.0 Mb
- 他 4 hap: 303-324 Mb

15-20 Mb の差 (5-7%) は許容範囲ですが、**同僚の Merqury k-mer 分析(第 4.7 節参照)**で、より深刻な問題が見つかっています。

### Merqury による発見

Merqury (Rhie et al., 2020) は k-mer の頻度から**アセンブリの重複度**を評価します。同僚の分析結果:

| アセンブリ | ユニーク k-mer% | 2回出現 % | 3回以上 % | 総塩基 |
|---|---|---|---|---|
| CKUhap2 (KNN hap2) | 90.17% | 5.46% | 4.37% | 303.4 Mb |
| CKIhap1 (KSH hap1) | 89.95% | 5.61% | 4.44% | 304.2 Mb |
| CKUhap1 (KNN hap1) | 89.42% | 5.77% | 4.81% | 323.9 Mb |
| CKIhap2 (KSH hap2) | 89.28% | 6.06% | 4.66% | 310.3 Mb |
| **CUNphKu (STS hap2)** | **82.64%** | **11.89%** | **5.47%** | **358.0 Mb** |
| **CUNphKi (STS hap1)** | **79.64%** | **14.65%** | **5.71%** | **348.4 Mb** |

**温州の 2 hap は明らかに異常**:
- ユニーク率が 10 pp 低い
- 2回出現する k-mer が親系統の**2倍**
- サイズも 30 Mb 大きい

これは **haplotype leakage** と呼ばれる現象で、trio phasing で 2 hap を完全に分離できなかった結果、**両 hap に同じ配列が残っている**状態です。

具体的には:
- 反復配列領域(centromere 近傍、TE クラスターなど)
- ヘテロ接合度が極端に低い/高い領域

で phasing が難しく、両 hap に collapsed に残ります。

### なぜこれを教材で扱うのか

**理由1: 実データは常に完璧ではない**
Isobe et al. 2023 は世界最高水準の trio phasing 実施例ですが、それでも完璧ではありません。教科書と実データのギャップを体験してもらえます。

**理由2: 後続の QC の解釈に必須**
Pangenome graph の第6章 QC で、この haplotype leakage が影響を及ぼします。事前に知っていれば正しく解釈できます。

**理由3: 教材の限界の宣言**
本教材で得られる graph は「haplotype leakage の影響を含む」ことを、正直に示すことになります。

---

## 4.7 追加の QC ツール(オプション)

### BUSCO / compleasm(完全性)

系統に近い遺伝子セット(orthologs)がどれだけ検出されるかで、**アセンブリの完全性**を評価します。

```bash
# compleasm (BUSCO より高速)
compleasm run -a data/satsuma/CUNphKi_r1.0.pmol.fasta.gz \
              -o busco/STS_hap1 \
              -l eudicots_odb10 \
              -t 16
```

期待値: **Complete > 95%、Duplicated < 5%**

### Merqury(k-mer QV と haplotype leakage)

HiFi 生リードから作った k-mer データベースと、アセンブリを比較して:
- **QV (base accuracy)**: 塩基精度、期待 > 40
- **完全性**: reads が捉えた k-mer のうち、アセンブリに含まれる割合
- **重複度**: haplotype leakage の検出

```bash
meryl count k=21 output hifi.meryl hifi_reads.fastq.gz
merqury.sh hifi.meryl CUNphKi_r1.0.fasta.gz CUNphKu_r1.0.fasta.gz STS
```

本教材では HiFi 生リードを扱わないため、Merqury は行いませんが、**結果の解釈は同僚のデータを引用**します(§4.6)。

---

## 4.8 この章のまとめ

- `seqkit stats` で基本統計を取得
- 柑橘 haploid は 300-360 Mb、9 chromosome、GC 34-38% を期待
- **STS の 2 hap は他より大きい**が、これは**haplotype leakage**を示唆
- 実データは完璧ではない。後の QC で影響を追跡できるよう、この情報を憶えておく

次章では、いよいよ **pangenome graph の構築**に入ります。

---

## 参考文献

- Shen W, et al. (2024). SeqKit2. *iMeta* 3:e191.
- Rhie A, et al. (2020). Merqury. *Genome Biol* 21:245.
- Manni M, et al. (2021). BUSCO update. *Mol Biol Evol* 38:4647-4654.
- Huang N, Li H (2023). compleasm. *Bioinformatics* 39:btad595.

---

[← 第3章](03_data_acquisition.md) | [第5章: Pangenome graph 構築 →](05_graph_construction.md)

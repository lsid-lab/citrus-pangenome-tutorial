# 第5章 Pangenome-graph の構築

## 学ぶこと

- 主要な pangenome 構築ツールの特徴と違い
- なぜこの教材で **PGGB** を選ぶか
- 「参照配列」という言葉の 3 つの意味を整理
- PGGB の内部工程と重要パラメータ
- 実行環境の選択肢と、この教材での実行手順

---

## 5.1 主要な pangenome 構築ツール

現在、pangenome graph 構築の実用的な選択肢は 3 つあります。

### 5.1.1 Minigraph

- 開発者: Heng Li (Broad Institute)
- 発表: Li 2020, *Genome Biology*
- 手法: 参照配列を骨格に、他 sample の**大規模 SV** (≥ 50 bp) だけを枝分かれとして追加
- 速度: **最速**、大規模データセット向け
- 特徴: 参照に強く依存(reference-biased)

**適する場合**: 数十〜数百 haploid で SV overview だけほしい時。SNP は無視される。

### 5.1.2 Minigraph-Cactus (MC)

- 開発者: UCSC Genomics Institute
- 発表: Hickey et al. 2024, *Nature Biotechnology*
- 手法: Minigraph で骨格を作り、Cactus で multi-way alignment を追加
- 速度: PGGB より速いが Minigraph より遅い
- 特徴: 参照を明示的に指定、SV と小変異両方を捕捉

**適する場合**: 大規模プロジェクト(HPRC ヒト 47 人など)、参照ベースの解釈が主眼のとき。

### 5.1.3 PGGB (PanGenome Graph Builder)

- 開発者: Erik Garrison ら
- 発表: Garrison et al. 2024 preprint
- 手法: **All-vs-all alignment**(wfmash)→ **graph induction**(seqwish)→ **smoothing**(smoothxg)
- 速度: **最遅**だが最も包括的
- 特徴: **Reference-free**(全 sample を対称に扱う)、すべての variant を捕捉

**適する場合**: 少数〜中規模(< 20 haploid)、reference-free で解析したいとき、pedigree 中心の観察。

### 5.1.4 比較表

|項目|Minigraph|Minigraph-Cactus|PGGB|
|---|---|---|---|
|参照依存|強い|中程度|**なし**|
|SV|◯|◎|◎|
|SNP|×|◯|◎|
|速度|◎|◯|△|
|スケール|数百 haploid|数十 haploid|< 20 haploid|
|SNP-level QC|不可|可|**可**|

---

## 5.2 「参照配列」の 3 つの意味 —— 混乱を整理する

Pangenome の議論では「reference」という言葉が**3 つの異なる意味**で使われ、混乱の元になります。ここで整理しておくと、後の議論がずっと楽になります。

### ① 構造的 backbone (structural reference)

Graph の骨格になる配列のこと。手法によって扱いが違います。

- **Minigraph**: 参照を軸に、他 sample の SV だけを枝分かれとして追加。**参照 bias が強い**
- **Minigraph-Cactus**: 参照を軸に multi-way alignment を作る。参照選びが結果に影響
- **PGGB**: 全 sample を対称に扱う **reference-free** アルゴリズム。この意味の参照は不要

### ② 座標アンカー (coordinate reference)

Graph ができた後、「〜番染色体の〜番目の塩基」と語るときの座標系。VCF を出したり、図示する時に必要です。

**PGGB でも、この意味の参照は 1 つ選ぶ必要があります**。ただし graph 自体の構造には影響しません。

### ③ 比較基準 (biological reference)

「サンプル A は B に対して〜が違う」と語るときのbaseline。論文の主張のフレーミングに関わります。

### 本教材での整理

本教材が採用する PGGB は **reference-free** なので:

- **① 構造的 backbone は不要**
- **② 座標アンカー**として、**CUN_hap1 (CUNphKi)** を選択(理由は後述)
- **③ 比較基準**は文脈依存

**「参照配列を選ぶ」と言うとき、我々が本当に選んでいるのは ② の座標アンカー**であって、graph の構造を決める骨格ではありません。この区別を意識してください。

### 座標アンカーとして CUN_hap1 を選ぶ理由

1. **CUN が pedigree の中心**: F1(CKI × CKU) なので、両親を graph 上でつなぐ位置
2. **hap1 は Kishu 由来と特定済み**: trio phasing の情報で解釈が明確
3. **Isobe 2023 アセンブリの中で最新・高品質**

---

## 5.3 なぜ PGGB を選ぶか

本教材では **PGGB を採用**します。理由:

### 理由1: Reference-free で pedigree 中心の観察に最適

3 品種の家族関係を見るには、特定の品種を「骨格」にすると解釈が歪みます。PGGB なら全 6 haploid を平等に扱えます。

### 理由2: 6 haploid なら計算が現実的

PGGB は遅い(9 染色体で 5-10 時間)ですが、6 haploid × 30-40 Mb/chr なら 32 CPU で対応可能。学生の実習環境でも動く。

### 理由3: SNP レベルまで捉える

小変異(SNP、indel)も含めて分析したいので、Minigraph は不適。MC でも可能ですが、PGGB の方が包括性が高い。

### 理由4: 教育的透明性

PGGB の 3 段構成(wfmash → seqwish → smoothxg)は、内部工程が明確に分かれていて、**学生が仕組みを理解しやすい**。デバッグや解釈にも有利。

---

## 5.4 PGGB の内部工程

PGGB は 3 段のパイプラインです:

```
[入力] 6 haploid FASTA (bgzip、PanSN命名)
         │
         ▼
┌────────────────┐
│  1. wfmash     │  All-vs-all pairwise alignment
│  (PAF 出力)    │  (どの部分がどれとアラインするか)
└────────────────┘
         │
         ▼
┌────────────────┐
│  2. seqwish    │  Alignment → graph 変換
│  (GFA 出力)    │  (共有部分を node に統合)
└────────────────┘
         │
         ▼
┌────────────────┐
│  3. smoothxg   │  Graph の正規化・単純化
│  (最終 GFA)    │  (POA でsmoothing)
└────────────────┘
         │
         ▼
[出力] final.gfa, final.og, VCF
```

---

## 5.5 PGGB の重要パラメータ

第1章で挙げた「扱いに注意が必要なパラメータ」を、PGGB のオプションとして具体化します。

### 5.5.1 主要パラメータ

|パラメータ|意味|推奨値|説明|
|---|---|---|---|
|`-n`|haploid 数|**6**|3 品種 × 2 hap|
|`-p`|percent identity|**95**|許容する塩基同一性 (>90% で align)|
|`-s`|segment length (bp)|**10000**|wfmash の最小 align 長|
|`-V`|VCF reference path|`CUN#1:#`|vg deconstruct の座標基準(§5.2 で決定)|
|`-Y`|PanSN separator|`#`|PanSN 名の区切り文字|
|`-t`|スレッド数|32-48|CPU 数に応じて|

### 5.5.2 パラメータの選定根拠

**`-p 95`(percent identity)**
- 柑橘は同種内 SNP 頻度が 17/kbp (Kiryu 2026)、つまり 98.3% 同一
- `-p 95` は「95% 以上の同一性が続く領域」を align 候補にする閾値
- `-p 90` でも動くが、95 の方が **偽陽性 align が少なく計算が速い**

**`-s 10000`(segment length)**
- wfmash の**最小 align 単位**
- 短いほど詳細だが、noise が増える
- 10 kb は柑橘の HiFi assembly には妥当

**`-n 6`(haploid 数)**
- **入力に含まれる haploid 数を正確に指定**
- wfmash が「各配列を上位 n-1 個の類似配列とアラインする」ため

### 5.5.3 バージョン依存の注意事項:`-B` パラメータ

本教材は執筆時点で **PGGB (Singularity image `pggb:latest`、リビジョン `4225c6c` 相当)** を対象としています。このバージョンには、`-B` (transclose-batch) パラメータの既定値のパース処理に問題があり、指定しないと seqwish が想定より小さいバッチサイズで動作します。結果として、graph が十分に圧縮されないことがあります。

**対処**: 実行時に `-B 1G` などを明示的に指定する。目安は、入力ファイル総塩基長 L 以上。最大染色体 chr03 が 306.7 Mb なので、`-B 1G` で余裕があります。

**将来のバージョン**では改善される可能性があります。実行後には第6章で `.params.yml` の `transclose-batch` 値を確認する習慣をつけましょう(第6章 §6.2.3)。

このバグの発見と記録は、共同研究者 町田 修斗による綿密な診断実験の成果です。

---

## 5.6 実行環境の選択肢

PGGB のような複雑なパイプラインを動かすには、依存関係が多いため**コンテナ化された環境の使用を強く推奨**します。以下の 3 通りが主な選択肢です。

### 選択肢A: Docker

Web の記事などで最も一般的に紹介される方法。

```bash
docker pull ghcr.io/pangenome/pggb:latest
docker run -v $PWD:/data ghcr.io/pangenome/pggb pggb ...
```

**適する場合**: 個人ワークステーション、root 権限あり

**注意点**: HPC ではセキュリティ上 Docker が使えないことが多い(多数ユーザ共有環境で root 権限を許可できない)

### 選択肢B: Conda / mamba

パッケージマネージャで直接インストール。

```bash
mamba create -n pggb -c bioconda pggb
mamba activate pggb
```

**適する場合**: 個人環境、他ソフトとの統合が必要な場合

**注意点**: 依存関係が複雑で、環境作成に失敗することがある

### 選択肢C: Singularity(本教材で採用)

HPC 環境での標準的な選択肢。Docker image をそのまま使え、root 権限不要。

```bash
singularity pull docker://ghcr.io/pangenome/pggb:latest
# → pggb_latest.sif が作られる
```

**適する場合**: HPC(SLURM/PBS)、多ユーザ共有環境

**注意点**: コンテナ内のツールを host から直接呼び出すため wrapper が必要(次節)

**本教材は選択肢 C (Singularity) を採用**します。以降の手順は Singularity 前提で書かれています。

---

## 5.7 実行手順 (Singularity)

### 5.7.1 Singularity image の取得と wrapper 作成

```bash
# image のダウンロード
singularity pull docker://ghcr.io/pangenome/pggb:latest

# host から SIF 内の全ツールを呼び出す wrapper を作成
bash scripts/setup_singularity_wrappers.sh $(readlink -f pggb_latest.sif)
```

これで `~/bin/{pggb,odgi,vg,samtools,bcftools,bgzip,...}` が作られます。

**PATH に追加**:

```bash
echo 'export PATH="$HOME/bin:$PATH"' >> ~/.bashrc
export PATH="$HOME/bin:$PATH"
```

**動作確認**:

```bash
pggb --version
odgi version
```

### 5.7.2 入力の整形

```bash
bash scripts/pg01_prepare_input.sh tables/samplesheet.tsv .
```

このスクリプトが行うこと:
1. 各 haploid FASTA を **PanSN 命名にリネーム**
   例: `CKUhap1_r1.0ch1` → `CKU#1#chr01`
2. 染色体番号を統一(`chr01` 〜 `chr09`)
3. 染色体別に FASTA を分割
4. bgzip 圧縮 + samtools faidx

出力: `03_pangenome/by_chr/chr01.fa.gz` 〜 `chr09.fa.gz`

各 FASTA には**きっちり 6 配列**(3 品種 × 2 hap)が含まれます。

### 5.7.3 テスト実行(小さい染色体で)

まず chr09(最小、テスト用)で試験実行:

```bash
bash scripts/pg02_run_pggb.sh . 32 chr09 chr09
```

30-60 分で完了します。

### 5.7.4 全染色体の並列実行 (SLURM)

Chr01 〜 Chr09 を**個別ジョブとして並列投入**します(job array ではなく、通常の sbatch を9回)。理由は、job array だと全ジョブに同じリソース指定がされてしまうため、実際のリソース使用状況に応じた柔軟な運用がしにくいためです。

各染色体用のジョブスクリプト `run_pggb_chr.sbatch`:

```bash
#!/bin/bash
#SBATCH --job-name=pggb_chr
#SBATCH --output=logs/pggb_%x_%j.out
#SBATCH --error=logs/pggb_%x_%j.err
#SBATCH --mem=64G
#SBATCH --cpus-per-task=32
#SBATCH --time=12:00:00

export PATH="$HOME/bin:$PATH"

CHR=${1:?染色体名 (例 chr01) を指定してください}
bash scripts/pg02_run_pggb.sh . 32 $CHR $CHR
```

投入(全 9 染色体):

```bash
for CHR in chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09; do
  sbatch run_pggb_chr.sbatch $CHR
done
```

これで 9 染色体が並列で走ります。実時間で 8-12 時間で完了する見込みです。

### 5.7.5 実行時間とメモリの目安 (32 CPU 使用時)

|染色体|入力サイズ|所要時間|MaxRSS|
|---|---|---|---|
|chr01|195.7 Mb|~28 min|12 GB|
|chr02|213.5 Mb|~44 min|6 GB|
|chr03|306.7 Mb|~62 min|18 GB|
|chr04|177.9 Mb|~19 min|5 GB|
|chr05|252.6 Mb|~43 min|14 GB|
|chr06|165.8 Mb|~24 min|10 GB|
|chr07|210.8 Mb|~74 min|6 GB|
|chr08|228.8 Mb|~56 min|13 GB|
|chr09|195.9 Mb|~29 min|4 GB|

これらは共同研究者の実測値です。最大でも 18 GB 程度なので、SLURM 指定の `--mem=64G` で十分な余裕があります。

---

## 5.8 出力ファイル

各染色体で以下が生成されます (chr09 の例):

```
03_pangenome/by_chr/chr09_pggb/
├── chr09.fa.gz.<hash1>.<hash2>.<hash3>.smooth.final.gfa    ← 最終 graph (GFA)
├── chr09.fa.gz.<hash1>.<hash2>.<hash3>.smooth.final.og     ← odgi 形式
├── chr09.fa.gz.<hash1>.<hash2>.<hash3>.smooth.final.CUN#1.vcf  ← VCF
├── chr09.fa.gz.<hash1>.<hash2>.<hash3>.smooth.final.og.lay.draw_multiqc.png  ← 2D 図
├── chr09.fa.gz.<hash1>.<hash2>.<hash3>.smooth.*.params.yml  ← 実行パラメータの正本
├── chr09.fa.gz.<hash1>.<hash2>.<hash3>.smooth.*.log        ← ログ
└── multiqc_report.html                                       ← 統合レポート
```

**ハッシュ**: パラメータのハッシュが埋め込まれます。異なるパラメータで再実行しても上書きされません。

---

## 5.9 この章のまとめ

- Pangenome 構築ツールは Minigraph、Minigraph-Cactus、PGGB の 3 択
- 「参照配列」は 3 つの意味を持つ。**PGGB では ②座標アンカーのみ 1 つ選ぶ**
- PGGB は 3 段構成: wfmash → seqwish → smoothxg
- パラメータの選定は品種の遺伝的距離に基づく: 柑橘は `-p 95 -s 10000 -n 6`
- 特定バージョンでは `-B` の明示指定が必要。`.params.yml` で事後確認
- HPC では Singularity が実用的
- 9 染色体は個別 sbatch ジョブとして並列投入

次章では、生成された graph の**品質を評価**します。

---

## 参考文献

- Garrison E, et al. (2024). Building pangenome graphs. *bioRxiv*.
- Li H (2020). The design and construction of reference pangenome graphs with minigraph. *Genome Biol* 21:265.
- Hickey G, et al. (2024). Pangenome graph construction from genome alignments with Minigraph-Cactus. *Nat Biotechnol* 42:663-673.
- Kiryu Y, et al. (2026). AlleleMiner. *DNA Res* 33:dsag004.

---

[← 第4章](04_data_qc.md) | [第6章: グラフの QC →](06_graph_qc.md)

# 第1章 Pangenome構築の大まかな手順

## 学ぶこと

- Pangenome 構築の一般的なワークフロー
- ツール横断で使う共通概念(PanSN、GFA、node/edge/path)
- 各ステップでどんな判断が必要か

---

## 1.1 全体像

どの pangenome ツール(PGGB、Minigraph-Cactus、Minigraph)を使っても、大まかな流れは共通です:

```
┌─────────────────────────────────┐
│  1. データ収集                    │  ← 論文/DBからアセンブリを取得
│  (各個体の genome assembly)       │
└──────────┬──────────────────────┘
           │
┌─────────────────────────────────┐
│  2. データQC                     │  ← 各アセンブリの品質確認
│  (連続性、完全性、コンタミ検査)   │
└──────────┬──────────────────────┘
           │
┌─────────────────────────────────┐
│  3. 入力整形                     │  ← 命名規則の統一、bgzip 化
│  (PanSN naming, chromosome分割)  │
└──────────┬──────────────────────┘
           │
┌─────────────────────────────────┐
│  4. Graph 構築                   │  ← PGGB / MC / Minigraph
│  (align, induction, refinement)  │
└──────────┬──────────────────────┘
           │
┌─────────────────────────────────┐
│  5. Graph QC                     │  ← 構造、パス保存性、生物学的妥当性
│  (odgi stats / similarity)       │
└──────────┬──────────────────────┘
           │
┌─────────────────────────────────┐
│  6. 解析・可視化                  │  ← SV抽出、VCF、可視化
│  (vg deconstruct, odgi viz)      │
└─────────────────────────────────┘
```

**この教材ではステップ 1 から 6 まで、すべて手を動かして実行します。**

---

## 1.2 共通の概念

### 1.2.1 PanSN-spec (Pangenome Sequence Naming Specification)

Pangenome ツールの多くは、配列名に**特別な命名規則**を要求します。

**PanSN 形式**: `sample#haplotype#contig`

例:
```
satsuma#1#chr01    # 温州みかん、hap1、1番染色体
satsuma#2#chr01    # 温州みかん、hap2、1番染色体
kishu#1#chr01      # 紀州みかん、hap1、1番染色体
```

`#` (シャープ) が区切り文字です。この命名で、ツールは**どの配列がどの品種のどのハプロタイプか**を機械的に判別できます。

多くのアセンブリの元の配列名は `CUNphKi_r1.0ch1` のような形式なので、pangenome 構築の前に**リネーム**が必要です(本教材の `pg01_prepare_input.sh` がこれを担います)。

### 1.2.2 GFA形式 (Graphical Fragment Assembly)

Pangenome graph の**標準ファイル形式**。中身は3種類のレコードで構成されます。

**S (Segment)**: node(配列の一続き)を定義
```
S  1  ACGTACGT     # node ID = 1、配列 = ACGTACGT
S  2  GCATT
S  3  TTAA
```

**L (Link)**: node 同士の接続(edge)を定義
```
L  1  +  2  +  0M    # node 1 の右 → node 2 の左
L  1  +  3  +  0M    # node 1 の右 → node 3 の左
```

**P (Path)**: ある配列が graph 内をどう通るかを記述
```
P  satsuma#1#chr01  1+,2+,4+,5+  ...
P  kishu#1#chr01    1+,3+,4+,5+  ...
```

このように、node が共有され、edge で分岐が表現され、path が個々の haplotype を示します。

### 1.2.3 Node / Edge / Path の関係

```
Path A: 1 → 2 → 4 → 5
Path B: 1 → 3 → 4 → 5

  node2      
  ↗   ↘     
node1  node4 → node5
  ↘   ↗     
  node3
```

- **Node** (S): 実際の DNA 配列を持つ最小単位
- **Edge** (L): node同士の連結(=次にどの node に進めるか)
- **Path** (P): ある個体の染色体が graph を通る**経路**

**重要な性質**:
- 同一の配列は node として**1つに統合**される(collapse)
- 違いがある場所で node が**分岐する**
- 各個体の元の配列は、**path を辿ることで完全に復元可能**

これが「reference-free」の意味です。特定の個体を軸にせず、全個体を平等に graph に組み込む。

### 1.2.4 圧縮率 (Compression Ratio)

Graph の品質を測る**最も重要な指標**:

```
圧縮率 = graph の総塩基数 / 入力ファイルの総塩基数
```

例:
- 入力: 6 haploid × 30 Mb = 180 Mb
- Graph の総塩基数: 45 Mb
- 圧縮率 = 45 / 180 = **0.25** (25%)

意味: 6 haploid の情報が graph の 25% のサイズに圧縮された = 75% の配列は**共有**されている。

|圧縮率|解釈|
|---|---|
|0.15-0.30 | 良好 (関連する近縁品種の pangenome) |
|0.30-0.50 | 普通 (種内でやや divergent) |
|0.50-0.80 | 圧縮不足 (パラメータ再検討) |
|> 0.90 | ほぼ圧縮なし (ツール/データに致命的問題) |

第6章の Graph QC で、この指標を実際に見て判断します。

---

## 1.3 各ステップの判断ポイント

### ステップ2 (データQC)
- 各アセンブリの**連続性** (N50 が chromosome scale か?)
- **完全性** (BUSCO complete > 95%?)
- **コンタミ** (他生物由来配列の混入は無いか?)
- **haplotype 分離** (2つのhapのサイズが揃っているか?)

### ステップ3 (入力整形)
- 全 haploid ファイルの**命名を統一** → PanSN 形式に変換
- **染色体別に分割**するか、全ゲノム一括か
  - 小さいゲノム(<300 Mb): 全ゲノム一括で OK
  - 中〜大サイズ(300 Mb〜): 染色体別が必要 (計算資源制約)

### ステップ4 (Graph 構築)
- どのツールを使うか(§5章で議論)
- **扱いに注意が必要なパラメータ**(ツール共通の概念として):
  - **入力の haploid 数**: 正確に指定しないとアライメント段階で問題が起きる
  - **percent identity**: 許容する塩基同一性の閾値。低いほど divergent な変異も拾えるがノイズが増える
  - **segment length**: アライメントの最小単位長。大きいほど大規模SVを重視、小さいほど小変異に敏感
  - その他、ツール固有の内部パラメータあり(PGGBについては第5章で詳しく扱う)

### ステップ5 (Graph QC)
- 圧縮率が期待範囲内か
- 全 haploid の path が保存されているか (欠損なし)
- 生物学的妥当性 (pedigree、系統関係)
- SV/SNP の検出数が現実的か

### ステップ6 (解析)
- VCF 生成、SV カウント
- 可視化 (odgi viz, odgi draw)
- Path 間の類似度分析
- 生物学的解釈

---

## 1.4 なぜ「染色体別」に構築するか

大きなゲノムの pangenome 構築は、**染色体を独立に扱う**のが一般的です。

理由:
1. **計算資源の制約**: 全ゲノム一括だと数百 GB のメモリが必要になる
2. **並列化**: 染色体ごとに独立ジョブとして流せる (HPC で有効活用)
3. **障害耐性**: ある染色体で失敗しても、他は影響を受けない
4. **染色体レベルの比較容易性**: 特定の染色体だけを深く見たいとき便利

前提として:
- 各アセンブリが**染色体スケール**でアセンブルされていること
- 染色体番号が**アセンブリ間で対応している**こと(chr01 は本当に同じ染色体か?)

柑橘 (n=9) は 9 染色体 × haploid 6 個 = **54 の配列**を graph に投入することになります。

---

## この章のまとめ

- Pangenome 構築は 6 ステップの一般的ワークフローに従う
- ツール横断の共通概念: **PanSN**、**GFA**、**node/edge/path**、**圧縮率**
- 大きなゲノムは染色体別に構築するのが実務標準

次章では、本教材で具体的に扱う3品種(温州・紀州・九年母)の生物学的背景を見ていきます。

## 参考文献

- Garrison E, et al. (2018). Variation graph toolkit improves read mapping by representing genetic variation in the reference. *Nat Biotechnol* 36:875-879. (vg toolkit)
- Li H (2020). The design and construction of reference pangenome graphs with minigraph. *Genome Biol* 21:265.
- Heumos S, et al. (2024). Pangenome graph layout by Path-Guided SGD. *Bioinformatics*.

---

[← 第0章](00_what_is_pangenome.md) | [第2章: 今回の対象 →](02_target_organisms.md)

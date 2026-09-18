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
CUN#1#chr01    # 温州みかん、hap1、1番染色体
CUN#2#chr01    # 温州みかん、hap2、1番染色体
CKI#1#chr01      # 紀州みかん、hap1、1番染色体
```

`#` (シャープ) が区切り文字です。この命名で、ツールは**どの配列がどの品種のどのハプロタイプか**を機械的に判別できます。

多くのアセンブリの元の配列名は `CUNphKi_r1.0ch1` のような形式なので、pangenome 構築の前に**リネーム**が必要です(本教材の `pg01_prepare_input.sh` がこれを担います)。

> **出典**: PanSN-spec は PGGB の開発者でもある **Erik Garrison** が提唱した仕様で、正式な定義は
> <https://github.com/pangenome/PanSN-spec> で公開されています。
>
> 仕様上の正確な形は `[sample_name][delim][haplotype_id][delim][contig_or_scaffold_name]` で、
> 区切り文字 `delim` は原理的には任意の文字ですが、**慣例上 `#` を使います**(本教材も `#` を使います)。
> `haplotype_id` は数値です。倍数性や phasing の無いアセンブリでも、`sample#0#chr01` のように
> haplotype 番号を埋めておくのが作法です。
>
> PanSN が狙っているのは「メタデータを別ファイルで持ち回らずに済ませる」ことです。配列名そのものに
> sample と haplotype を埋め込んでおけば、FASTA・GFA・VCF・BED・GFF のどの形式を経由しても
> その情報が失われません。だからこそ、**pangenome 解析では最初のリネームを手抜きしてはいけません**。

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
L  1  +  2  +  0M    # node 1 の右端 → node 2 の左端
L  1  +  3  +  0M    # node 1 の右端 → node 3 の左端
L  1  +  3  -  0M    # node 3 を逆向き(reverse complement)で繋ぐ
```

列の意味は `L <From> <FromOrient> <To> <ToOrient> <Overlap>` です。

- `<FromOrient>` / `<ToOrient>`: `+` = そのまま、`-` = reverse complement
- `<Overlap>`: **2つの node が何塩基重なっているかを CIGAR 文字列で書く欄**(GFA1 では必須列)

#### CIGAR とは何か

**CIGAR** (Compact Idiosyncratic Gapped Alignment Report) は、**2つの配列がどう対応づいているかを「回数 + 操作記号」の並びで表す記法**です。SAM/BAM でアライメントを記述するのに使われているものと同じ形式で、GFA もそれを流用しています。

|記号|意味|
|---|---|
|`M`|alignment match(一致・不一致のどちらもありうる「対応している」区間)|
|`=` / `X`|完全一致 / ミスマッチ(`M` を細分化したもの)|
|`I`|insertion(片方にだけ余分な塩基がある)|
|`D`|deletion(片方で塩基が欠けている)|

例: `10M2D5M` は「10塩基が対応 → 2塩基の欠失 → 5塩基が対応」と読みます。

L 行における CIGAR は、**前の node の末尾と次の node の先頭が何塩基オーバーラップしているか**を表します。

- `0M` = **オーバーラップ 0 塩基**、つまり突き合わせ (blunt join)
- `*` = 未指定

**pangenome graph の GFA では、L 行の CIGAR は基本的にすべて `0M`** です。PGGB / Minigraph-Cactus / minigraph が作る graph は、配列を node 境界でぶつ切りにして並べているだけで、node 同士が塩基を共有していないためです。一方、OLC 系のアセンブラ(Canu、miniasm など)が出す**アセンブリ**グラフの GFA には `120M` のような**本物のオーバーラップ**が並びます。手元の `.gfa` が pangenome graph なのかアセンブリグラフなのか分からないときは、L 行の最終列を見るのが手っ取り早い見分け方になります。

**P (Path)**: ある配列が graph 内をどう通るかを記述
```
P  CUN#1#chr01  1+,2+,4+,5+  *
P  CKI#1#chr01    1+,3+,4+,5+  *
P  CKU#1#chr01  1+,3-,4+,5+  *    # node 3 を逆向きに通る = inversion
```

列の意味は `P <PathName> <SegmentNames> <Overlaps>` です。

- `<SegmentNames>`: 「node ID + 向き」をカンマ区切りで並べたもの
  - `3+` = node 3 を **forward** で通る
  - `3-` = node 3 を **reverse complement(逆向き)** で通る
- `<Overlaps>`: node 間の CIGAR のリスト。オーバーラップを持たない pangenome graph では `*` と書くのが普通です

**`-` が現れるのは inversion(逆位)を通るとき**です。上の例では CKU だけが node 3 を反転して持っている、つまりこの区間に逆位がある、と読めます。教科書的な図では `1+,2+,3+` のように `+` ばかりが並びますが、**実データの graph には `-` を含む path が必ず出てきます**。逆位や、より複雑な再編成 (inverted duplication など) は、この向きの情報としてしか graph 上に現れません。「path とは node の列ではなく、**向き付きの** node の列である」という点が、単なる線形配列との大きな違いです。

このように、node が共有され、edge で分岐が表現され、path が個々の haplotype を(向きまで含めて)示します。

#### 注意: GFA v1.0 / GFA v1.1 / rGFA は「同じ GFA」ではない

拡張子はどれも `.gfa` ですが、**中身の方言が異なり、下流ツールの対応状況も違います**。ここを取り違えると「ファイルは読めているのに path が 0 本」といった事故が起きます。

|  |主に出力するツール|path の表現|
|---|---|---|
|**GFA v1.0**|**PGGB**|`P` 行。最も互換性が高く、大概のツールにそのまま渡せる|
|**GFA v1.1**|**Minigraph-Cactus (MC)**|`P` 行の代わりに **`W` 行 (Walk)**。`W` 行に未対応のツールは path を見落とす|
|**rGFA**|**minigraph**|**path 行そのものが無い**。`S` 行のタグで参照座標を持つ|

**GFA v1.1 の `W` 行 (Walk)**

```
W  <SampleId>  <HapIndex>  <SeqId>  <SeqStart>  <SeqEnd>  <Walk>

W  CUN  1  chr01  0  30512000  >1>2>4>5
```

`P` 行が `1+,2+` とカンマ区切りで書くところを、`W` 行は `>1>2` と書きます(`>` = forward、`<` = reverse。役割は `+` / `-` と同じです)。また、PanSN で `CUN#1#chr01` と1本の文字列に詰め込んでいた情報が、sample / haplotype / contig の**3つの列に分解**されています。

**rGFA (reference GFA)**

minigraph が出力する rGFA は GFA の**厳密なサブセット**で、次の性質を持ちます。

- `S` 行に 3 つのタグが**必須**: `SN:Z:`(元になった stable sequence 名)、`SO:i:`(その配列上のオフセット)、`SR:i:`(rank。`0` = 線形参照ゲノム由来、`>0` = 非参照由来)
- segment 間のオーバーラップを**許さない**
- **`P` 行も `W` 行も持たない**

つまり rGFA は「どの sample がどう通ったか」ではなく「**参照配列のどの座標に相当するか**」で位置を表す、**参照中心**の設計です(§5.2 で整理する「① 構造的 backbone」がファイル形式のレベルで焼き込まれている、と考えると分かりやすいでしょう)。そのため、path があることを前提にした odgi の各種解析や `vg deconstruct` による VCF 化を、rGFA にそのまま掛けることはできません。

**実務上の注意: MC の GFA は downgrade が必要になることがある**

PGGB は GFA v1.0(`P` 行)を出力するので、その後の処理は大概のツールでそのまま通ります。一方 **MC は GFA v1.1(`W` 行)を出力する**ため、`W` 行に未対応のツールへ渡す前に、`P` 行へ変換(downgrade)しておく必要があります。

```bash
# W 行 → P 行 (GFA v1.1 → v1.0 相当) に変換
vg convert -g graph.mc.gfa -f -W > graph.p-lines.gfa
#   -g : 入力が GFA
#   -f : GFA で出力
#   -W : W 行を使わず、すべての path を P 行として書き出す (--no-wline)
```

本教材は PGGB を使うのでこの変換は登場しませんが、**MC で作られた公開 graph(HPRC のヒト pangenome など)を触るときには必ず引っかかる**ポイントです。「path が見つからない」と言われたら、まず

```bash
grep -c '^P' graph.gfa   # P 行の数
grep -c '^W' graph.gfa   # W 行の数
```

を比べてみてください。

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

#### BUSCO —— 「完全性」をどう測るか

**BUSCO** (Benchmarking Universal Single-Copy Orthologs, <https://busco.ezlab.org/>) は、アセンブリの**完全性 (completeness)** を評価する事実上の標準ツールです。

考え方は単純です。「**ある系統群のほぼ全種が、1コピーずつ持っているはずの遺伝子セット**」をあらかじめ用意しておき、それが手元のアセンブリから何個見つかるかを数えます。見つからない遺伝子が多ければ、その分だけアセンブリが欠けている(あるいは壊れている)と判断できる、という理屈です。N50 が「**どれだけ繋がっているか**」を測るのに対し、BUSCO は「**中身が揃っているか**」を測ります。両方見て初めてアセンブリの品質が語れます。

結果は次の 4 カテゴリの割合として報告されます。

|カテゴリ|意味|読み方|
|---|---|---|
|**S** (Complete, single-copy)|完全な形で 1 コピーだけ見つかった|ここが高いほど良い|
|**D** (Complete, duplicated)|完全な形で見つかったが 2 コピー以上あった|高い場合は **haplotype の重複 (leakage)** や偽の重複を疑う|
|**F** (Fragmented)|部分的にしか見つからない|断片化・アセンブリエラーの指標|
|**M** (Missing)|全く見つからない|欠損の指標|

**C (Complete) = S + D** で、`C:97.2%[S:94.1%,D:3.1%],F:1.0%,M:1.8%` のような 1 行に要約されます。本教材で扱うような haplotype-resolved アセンブリでは、**hap ごとに別々に** BUSCO を回します。片方の hap にだけ **D が高く出たら、その hap に同じ配列が重複して入っている(haplotype leakage)可能性を疑います**。

**一番大事なのは lineage dataset の選び方**

BUSCO の数字は、**どの系統データセット (lineage dataset) を使ったかで変わります**。原則は次の 3 点です。

1. **対象生物を含む中で、最も下位(specific)の系統を選ぶ**
   上位の系統ほど「全種が共通して持つ遺伝子」は少なくなるため、比較対象の遺伝子数が減って**判定が甘く、情報量も少なく**なります。例えば `eukaryota_odb10` は約 255 遺伝子しかありませんが、`eudicots_odb10` は約 2,326 遺伝子あり、後者の方がはるかに厳しく解像度の高い評価になります。

2. **手元で使える一覧を必ず確認する**
   ```bash
   busco --list-datasets
   ```

3. **系統が分からなければ自動判定に任せる**
   ```bash
   busco -i assembly.fasta -m genome --auto-lineage-euk -o out
   ```
   ただし自動判定は時間がかかるので、系統が分かっているなら `-l` で明示する方が速く確実です。

柑橘 (*Citrus*、ミカン科 Rutaceae) の場合、植物の系統データセットは

```
eukaryota  →  viridiplantae  →  embryophyta  →  eudicots
 (約255)       (約425)           (約1,614)        (約2,326 遺伝子)
```

の順に specific になっていきます。**Rutaceae やムクロジ目 (Sapindales) に特化したデータセットは用意されていない**ため、柑橘で選ぶべきは最も下位にある **`eudicots_odb10`** です。`brassicales`(アブラナ目)や `fabales`(マメ目)といった目レベルのデータセットも存在しますが、**柑橘はそれらに含まれないので選んではいけません**。「より specific なら何でも良い」のではなく、**対象生物を含んでいることが大前提**です。

```bash
busco -i data/CUN/CUNphKi_r1.0.pmol.fasta -m genome \
      -l eudicots_odb10 -c 16 -o busco/CUN_hap1
```

> データセット名末尾の `odb10` は元になった OrthoDB のバージョンです。新しい BUSCO では `odb12` 系も提供されているので、`busco --list-datasets` で手元の BUSCO が何を持っているか確認してください。
>
> **データセットもバージョンも違えば、BUSCO の数字は比較できません。** 論文や他のアセンブリと BUSCO 値を突き合わせるときは、必ず lineage dataset 名と BUSCO のバージョンを揃えてください。

実際のコマンド(BUSCO より高速な `compleasm` を使う版)は §4.6 で扱います。

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

### 例外: 「染色体別」にできない対象もある

染色体別の構築が成り立つのは、**染色体をまたいで配列が入れ替わっていない**という前提があるからです。この前提が崩れる対象では、**全ゲノムを一括で**処理しなければなりません。

代表例が**がんゲノム**です。がん細胞では転座 (translocation) や chromothripsis といった**染色体間の組み換え**が日常的に起きています。染色体ごとに切り分けて graph を作ってしまうと、その繋ぎ目 —— つまり**最も見たい構造変化そのもの** —— が graph から丸ごと落ちてしまいます。同じことは、染色体の融合・分裂を伴う種間比較や、異質倍数体でサブゲノム間の組み換えを見たい場合にも当てはまります。

ただし、**一括での構築は巨大な計算リソースと計算時間を要求します**。染色体別なら 9 本を並列に流して 1 本あたり数十 GB のメモリで済むところが、一括では全配列を一度に載せた all-vs-all alignment になるため、**数百 GB 〜 TB 級のメモリと、数日〜数週間の実行時間**を覚悟する必要があります(しかも並列化が効かず、途中で失敗すれば全部やり直しです)。

したがって「染色体別にするか、一括にするか」は、**染色体間の組み換えを見る必要があるか**と、**そのコストを払えるか**を天秤にかけた判断になります。本教材が扱う柑橘品種間の比較では染色体の対応が保たれているため、染色体別で問題ありません。

---

## この章のまとめ

- Pangenome 構築は 6 ステップの一般的ワークフローに従う
- ツール横断の共通概念: **PanSN**、**GFA**、**node/edge/path**、**圧縮率**
- **GFA には方言がある**: PGGB = v1.0 (`P` 行)、MC = v1.1 (`W` 行)、minigraph = rGFA (path 行なし)。MC の graph は `vg convert -f -W` で downgrade が必要になることがある
- path は「**向き付きの** node の列」。`-` は inversion を意味する
- 完全性の評価 (BUSCO) は、**対象生物を含む中で最も specific な lineage dataset** を選ぶ(柑橘なら `eudicots_odb10`)
- 大きなゲノムは染色体別に構築するのが実務標準。ただし**染色体間の組み換えがある対象(がんゲノムなど)は一括処理が必要**で、コストは桁違いになる

次章では、本教材で具体的に扱う3品種(温州・紀州・九年母)の生物学的背景を見ていきます。

## 参考文献

- Garrison E, et al. (2018). Variation graph toolkit improves read mapping by representing genetic variation in the reference. *Nat Biotechnol* 36:875-879. (vg toolkit)
- Li H (2020). The design and construction of reference pangenome graphs with minigraph. *Genome Biol* 21:265.
- Hickey G, et al. (2024). Pangenome graph construction from genome alignments with Minigraph-Cactus. *Nat Biotechnol* 42:663-673.
- Heumos S, et al. (2024). Pangenome graph layout by Path-Guided SGD. *Bioinformatics*.
- Manni M, et al. (2021). BUSCO Update: novel and streamlined workflows. *Mol Biol Evol* 38:4647-4654.

### 仕様書・オンライン資料

- **PanSN-spec** (Erik Garrison): <https://github.com/pangenome/PanSN-spec>
- **GFA 仕様 (v1.0 / v1.1)**: <https://github.com/GFA-spec/GFA-spec>
- **rGFA 仕様** (Heng Li, gfatools): <https://github.com/lh3/gfatools/blob/master/doc/rGFA.md>
- **BUSCO** (lineage dataset 一覧・ユーザーガイド): <https://busco.ezlab.org/>
- **vg convert** (GFA の相互変換): <https://github.com/vgteam/vg/wiki>

---

[← 第0章](00_what_is_pangenome.md) | [第2章: 今回の対象 →](02_target_organisms.md)

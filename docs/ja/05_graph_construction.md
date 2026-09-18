# 第5章 Pangenome-graph の構築

## 学ぶこと

- 主要な pangenome 構築ツールの特徴と違い
- なぜこの教材で **PGGB** を選ぶか
- 「参照配列」という語の 3 つの意味と、PGGB での扱い
- PGGB の内部工程と重要パラメータ、**事前知識が無いときのパラメータの決め方**
- 実行環境の選択肢と、この教材での実行手順
- アラインメントとグラフ誘導を分ける **impg** の考え方

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

**適する場合**: 参照座標系での解釈が主眼のとき、リードマッピングに使う graph が欲しいとき、大規模プロジェクト(HPRC など)。

> 実際の実行手順は [付録: Minigraph-Cactus で同じデータを扱う](appendix_minigraph_cactus.md) にまとめてあります。

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
|参照配列の指定|**あり**(必須)|**あり**(必須)|**なし**|
|SV|◯|◎|◎|
|SNP|×|◯|◎|
|速度|◎|◯|△|
|スケール|数百 haploid|数十 haploid|< 20 haploid|
|SNP-level QC|不可|可|**可**|

---

## 5.2 「参照配列」の 3 つの意味と、PGGB での扱い

Pangenome の議論では「参照 (reference)」という語が **3 つの異なる意味**で使われます。混同すると話が噛み合わなくなるので、先に区別しておきます。

|意味|内容|
|---|---|
|**① 構造的 backbone**|graph の骨格になる配列。ここに無い配列は「枝分かれ」として付け足される|
|**② 座標アンカー**|「chr09 の 1,234 番目」と位置を言うための座標系。VCF や可視化に必要|
|**③ 比較基準**|「A は B に対して〜が違う」と語るときの baseline。論文の主張のフレーミング|

ツールによって、①を要求するかどうかが違います。Minigraph と Minigraph-Cactus は**① を必須**とし、指定した参照が graph の構造そのものに影響します。

### PGGB の場合

**PGGB は ① を持ちません。** 全配列を対称に扱うので、**graph の構造は、どれを「参照」と呼ぶかに依存しません**。

一方で、**② は PGGB でも 1 つ必要です。** 「chr09 の 1,234 番目の塩基が A→G」と VCF に書くには、その「1,234 番目」を数える path を 1 本決めなければなりません。PGGB ではこれを `-V` で指定します。

本教材は **`CUN#1`(CUNphKu、九年母由来)** を使います。`-V CUN#1` とすると、pggb は内部で `vg deconstruct -P CUN#1#` を呼び、CUN#1 の座標で VCF を書きます。選んだ理由は実務的なものです。

- **pedigree の中心にある**: 温州は F1 なので、両親のどちらとも比較しやすい
- **親由来が分かっている**: trio phasing により `CUN#1` は九年母 (CKU) 由来と特定済み(§2.5)
- **染色体スケールで欠損が少ない**

`CUN#2` を選んでも **graph は同じ**で、VCF の座標と REF 列が変わるだけです。これが「① を持たない」ということの実際の意味です。第6・7章の解析はすべて `CUN#1` 基準で書かれているので、変更する場合は `-V` と `vg deconstruct -P` の両方を揃えてください。

③ は解析結果ではなく議論の立て方の問題なので、第7章で結果を解釈するときに改めて意識します。

## 5.3 PGGB 採用の理由と、MC の場合のメリット

本教材は **PGGB を採用**します。

### PGGB を選ぶ理由

**1. 全 haploid を対称に扱える**

3 品種を対等に見る、という本教材の枠組みに、reference-free の構築が合っています。どれかを骨格(§5.2 の ①)に選ぶと、「骨格にある配列」と「骨格に無い配列」が非対称に扱われます。

**2. 内部工程が分解して観察できる**

wfmash → seqwish → smoothxg の 3 段が明確に分かれていて、中間ファイル(PAF、誘導直後の GFA)をそれぞれ見られます。**何が起きているかを追える**のは教材として大きい。§5.5.3 の `-B` のような問題に気づけるのも、工程が分かれているからです。

**3. 6 haploid なら計算が現実的**

PGGB は 3 ツールの中で最も遅いものの、6 haploid × 30-40 Mb/染色体なら 32 CPU で 30-60 分です。学生の実習環境でも動きます。

### MC を選んだ場合のメリット

同じデータを Minigraph-Cactus (MC) で扱うと、次の点が得られます。**目的によってはこちらが適切です。**

**1. 参照座標系が最初から固定される**

MC は参照ゲノムの指定が**必須**で、その参照は clip もされず cycle も持たず、**VCF の座標系と染色体分割をそのまま規定します**(§5.2 の ① と ② を兼ねる)。座標系を 1 本に固定して解析を進めたい場合は、後から `-V` で選ぶ PGGB より素直です。

**2. リードマッピング用の成果物がそのまま出る**

GFA / VCF に加えて **GBZ と `.dist` / `.min` / `.hapl` インデックス**が出力され、**vg Giraffe でのリードマッピングにそのまま使えます**。§7.6 で触れた v1/v2 比較のように、独立したリードデータを graph に載せて検証したい場合はこちらが向きます。

**3. 入力の softmask を気にしなくてよい**

MC は入力の softmask を要求しません。本教材では r2.0 が repeat-mask 済みだったため `pg01_prepare_input.sh` に大文字化の前処理を足しましたが(§7.6)、MC ならこの手当てが要りません。

**4. 速い**

同種内の近縁サンプル向けの設計で、この規模では PGGB より速く終わります。

**選び方は「出したい成果物から逆算する」**のが実際的です。参照座標での VCF とリードマッピング用インデックスが要るなら MC、全サンプルを対称に扱いたいなら PGGB、大規模 SV の概観だけなら Minigraph です。

> MC で同じ 6 haploid を扱う手順は [付録: Minigraph-Cactus で同じデータを扱う](appendix_minigraph_cactus.md) にあります。MC は GFA v1.1(`W` 行)を出すので、§1.2.2 で見た `vg convert` による downgrade がそこで実際に必要になります。

なお、**PGGB を使う場合でも、いまの標準的な流れは「PGGB 単体」ではありません**。アラインメントとグラフ誘導を分離する **impg** と組み合わせるのが現在の実務で、これは §5.9 で扱います。

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
|`-V`|VCF の座標基準|`CUN#1`|§5.2。書式は `REF[:LEN]`|
|`-t`|スレッド数|32-48|CPU 数に応じて|

> **`-V` の書式に注意**: `REF[:LEN]` です。pggb が REF に PanSN の `#` を自動で付けるので、`-V CUN#1` と書けば `vg deconstruct -P CUN#1#` が走ります。**`:` の後ろは区切り文字ではなく数値 (LEN)** で、0 より大きい値を与えると `vcfbub` + `vcfwave` による分解版 VCF も追加で出ます。本教材では分解版は使いません。
>
> **`-Y` は指定しません。** これは PanSN の区切り文字を教えるオプションではなく `--exclude-delim` で、「query と target が、指定文字の**最後の**出現より前で一致するマッピングをスキップする」ものです。既定値が `#` で、完全な PanSN 名ではグループが `sample#hap` 単位になるため、`CUN#1` と `CUN#2` は別グループとしてきちんとアラインされます。既定のままでよいので明示しません。

### 5.5.2 パラメータをどう決めるか

ここが実際にいちばん困るところです。柑橘については先行研究の値が使えますが、**自分の対象種に事前知識が無いときにどうするか**を先に書きます。

#### `-p`(percent identity)—— 自分のデータから測る

`-p` は「この同一性以上で続く領域を align 候補にする」閾値です。文献値に頼らず、**手元の配列同士の類似度を実際に測って**から決められます。

```bash
# 6 haploid 間の距離を総当たりで推定(高速、アラインメント不要)
mash triangle data/input/chr09.fa.gz > chr09.mash.tri

# あるいは wfmash のマッピング段階だけを走らせ、同一性の分布を見る
wfmash -m -t 16 data/input/chr09.fa.gz \
  | awk '{for(i=13;i<=NF;i++) if($i ~ /^(id|gi):f:/) {split($i,a,":"); print a[3]}}' \
  | sort -n | uniq -c
```

決め方の原則は **「align させたいペアの中で最も低い同一性より、数ポイント下に置く」**。低くしすぎると無関係な領域まで align 候補になって計算が膨らみ、高くしすぎると本当は相同な divergent 領域を取りこぼします。

柑橘の場合、同種内の SNP 頻度が 17/kbp 程度(Kiryu et al. 2026)、つまり **98.3% 程度は同一**です。ここから `-p 95` は「取りこぼさない側に十分な余裕がある」値だと判断できます。**この 1 行が、上のコマンドで測るべきものの答え合わせになっています。**

#### `-s`(segment length)—— 既定値から始めて、1 染色体で振ってみる

`-s` は wfmash がマッピングを探す最小単位の長さです。

- **短くする**と細かい相同性まで拾いますが、反復配列の中で一意に決まらない seed が増え、ノイズと計算時間が増えます
- **長くする**と安定しますが、それより短い再編成を見落とします

**10 kb は「柑橘だから」ではなく、染色体スケールのアセンブリで広く使われている出発点**です。対象種に合わせて確かめたいなら、**計算の軽い染色体 1 本で `-s` だけを振って比べる**のが実務的です(柑橘なら chr09。§5.7.5 の実測で最もメモリが少なくて済みます)。

```bash
for S in 5000 10000 20000; do
  pggb -i data/input/chr09.fa.gz -o results/pggb/chr09_s${S} \
       -p 95 -s $S -n 6 -B 1G -t 32
done

# 圧縮率とノード数で比べる(第6章 §6.3 の指標)
for S in 5000 10000 20000; do
  odgi stats -i results/pggb/chr09_s${S}/*.smooth.final.og -S
done
```

**圧縮率が極端に悪化せず、ノード数が無用に増えない範囲で、いちばん大きい `-s`** を選ぶのが目安です。

#### `-n`(haploid 数)—— 現行 pggb では省略できる

現行の pggb は、指定しなければ **PanSN 名から自動で数えます**(`sample#hap` の異なり数)。

したがって `-n 6` は必須ではなく、**自動検出の答え合わせ**として書いています。数が合わなければ入力の整形に失敗しているということなので、そこで気づけます。

### 5.5.3 バージョン依存の注意事項:`-B` パラメータ

本教材は執筆時点で **PGGB (Singularity image `pggb:latest`、リビジョン `4225c6c` 相当)** を対象としています。このバージョンには、`-B` (transclose-batch) パラメータの既定値のパース処理に問題があり、指定しないと seqwish が想定より小さいバッチサイズで動作します。結果として、graph が十分に圧縮されないことがあります。

**対処**: 実行時に `-B 1G` などを明示的に指定する。目安は、入力ファイル総塩基長 L 以上。最大染色体 chr03 が 306.7 Mb なので、`-B 1G` で余裕があります。

**将来のバージョン**では改善される可能性があります。実行後には第6章で `.params.yml` の `transclose-batch` 値を確認する習慣をつけましょう(第6章 §6.2.3)。

このバグの発見と記録は、共同研究者 町田 宗聡による綿密な診断実験の成果です。

**ここで効いてくるのがパイプラインの形です。** PGGB は「アラインメント → グラフ誘導 → smoothing」を 1 コマンドの中で通しで実行します。そのため、グラフ誘導のパラメータを間違えたことに後で気づくと、**最も高価な wfmash のアラインメントからやり直し**になります。この痛みを減らす方向に進んでいるのが §5.9 の impg です。

---

## 5.6 実行環境の選択肢

PGGB は依存関係が多いため、**コンテナ化された環境の使用を強く推奨**します。

### 選択肢A: Singularity / Apptainer(本教材で採用)

HPC 環境での標準。Docker image をそのまま使え、**特権が不要**です。

```bash
singularity pull docker://ghcr.io/pangenome/pggb:latest
# → pggb_latest.sif が作られる

# そのまま実行できる(カレントディレクトリは自動でマウントされる)
singularity exec pggb_latest.sif pggb --version
```

**適する場合**: HPC(SLURM/PBS)、多ユーザ共有環境。本教材はこれを前提に書かれています。

### 選択肢B: Docker

```bash
docker pull ghcr.io/pangenome/pggb:latest

docker run --rm -it \
  -u "$(id -u):$(id -g)" \
  -v "$PWD":/data -w /data \
  ghcr.io/pangenome/pggb pggb --version
```

**必要なのは root 権限ではなく、`docker` グループに入っていること**(または rootless Docker が設定されていること)です。とはいえ `docker` グループへの所属は実質的に root 相当の権限を与えるため、**多ユーザの共有 HPC では許可されないことがほとんど**です。

`-u "$(id -u):$(id -g)" ` を付けないと、**出力ファイルが root 所有で作られて後から消せなくなります**。忘れないでください。`--rm` は終了時にコンテナを破棄、`-w /data` は作業ディレクトリの指定です。

### 選択肢C: Podman

Docker とほぼ同じコマンドで動き、**既定で rootless** です。daemon も不要なので、HPC でも使えることがあります。

```bash
podman run --rm -it \
  -v "$PWD":/data -w /data \
  ghcr.io/pangenome/pggb pggb --version
```

rootless podman ではコンテナ内の root がホストの自分自身にマップされるため、`-u` は通常不要です(出力は自分の所有になります)。

### 選択肢D: Conda / mamba

```bash
mamba create -n pggb -c conda-forge -c bioconda pggb
mamba activate pggb
```

**注意点**: 依存が複雑で、環境の解決に失敗することがあります。他のツールと組み合わせたいとき以外は、コンテナのほうが確実です。

---

## 5.7 実行手順

### 5.7.1 image の取得

```bash
singularity pull docker://ghcr.io/pangenome/pggb:latest
SIF=$(readlink -f pggb_latest.sif)

singularity exec "$SIF" pggb --version
singularity exec "$SIF" odgi version
```

**以降のコマンドはすべて `singularity exec "$SIF" <tool> ...` の形で実行できます。** 毎回打つのが煩わしければ、ラッパーを作るスクリプトを用意してあります(任意):

```bash
bash scripts/setup_singularity_wrappers.sh "$SIF"
export PATH="$PWD/bin:$PATH"     # このセッションだけ有効
```

これで `./bin/{pggb,odgi,vg,samtools,bgzip,...}` が作られ、`pggb ...` と直接書けるようになります。**作られるのはプロジェクト内の `./bin` だけ**で、`~/bin` や `~/.bashrc` には一切書き込みません。恒久的に使いたければ、自分の判断で `~/.bashrc` に `export PATH` を足してください。

以下では、ラッパーを使わない `singularity exec` の形で書きます。

### 5.7.2 入力の整形(前処理)

PanSN へのリネームは **解析の前処理**なので、出力は成果物ではなく**入力**として `data/input/` に置きます。

```bash
bash scripts/pg01_prepare_input.sh tables/samplesheet.tsv .
```

このスクリプトが行うこと:

1. 各 haploid FASTA を **PanSN 命名にリネーム** — 例: `CKUhap1_r1.0ch1` → `CKU#1#chr01`
2. 染色体番号を `chr01` 〜 `chr09` に統一
3. 配列を**大文字に統一**(repeat-mask 済みの r2.0 と条件を揃えるため。§7.6)
4. 染色体別に分割
5. bgzip 圧縮 + `samtools faidx`

出力: `data/input/chr01.fa.gz` 〜 `chr09.fa.gz`。各ファイルに**きっちり 6 配列**が入ります。

> `seqkit` / `bgzip` / `samtools` が無い場合、このスクリプトは**冒頭で明示的に失敗します**。以前は無言でスキップして `.fa.gz` も `.fai` も作らないまま正常終了し、ずっと後の pggb で初めて失敗していました。

確認:

```bash
singularity exec "$SIF" seqkit seq -n data/input/chr09.fa.gz
# CUN#1#chr09
# CUN#2#chr09
# CKI#1#chr09
# ...
```

### 5.7.3 1 染色体を走らせる

**これが本体です。** まずは 1 本だけ走らせて確かめます。ここでは **chr09** を使います —— 第6・7章の例もすべて chr09 で書かれているので、続けて読むときに揃います。

```bash
singularity exec "$SIF" pggb \
  -i data/input/chr09.fa.gz \
  -o results/pggb/chr09 \
  -p 95 \
  -s 10000 \
  -n 6 \
  -B 1G \
  -V CUN#1 \
  -t 32 \
  -m -S
```

各オプション:

|オプション|意味|
|---|---|
|`-i`|入力 FASTA(bgzip + faidx 済み)|
|`-o`|出力ディレクトリ|
|`-p 95`|マッピングの同一性閾値(§5.5.2)|
|`-s 10000`|segment length(§5.5.2)|
|`-n 6`|haploid 数。省略すると PanSN から自動検出(§5.5.2)|
|`-B 1G`|seqwish の transclose batch。**明示必須**(§5.5.3)|
|`-V CUN#1`|VCF の座標基準(§5.2)|
|`-t 32`|スレッド数|
|`-m`|MultiQC 用の統計を出す|
|`-S`|各段階の統計を出す|

30-60 分で終わります。終わったら必ず:

```bash
# 指定したパラメータが本当に効いたかを、pggb 自身の記録で確認する
grep -E 'transclose-batch|map-pct-id|segment-length|n-haplotypes' \
  results/pggb/chr09/*.params.yml
```

**コマンドラインに書いたことではなく、ツールが書き出した記録を見る。** §5.5.3 の `-B` の件があるので、これは習慣にしてください。

### 5.7.4 全染色体

同じコマンドを染色体ごとに繰り返すだけです。

```bash
for CHR in $(ls data/input/*.fa.gz | xargs -n1 basename | sed 's/\.fa\.gz$//'); do
  singularity exec "$SIF" pggb \
    -i "data/input/${CHR}.fa.gz" \
    -o "results/pggb/${CHR}" \
    -p 95 -s 10000 -n 6 -B 1G -V CUN#1 -t 32 -m -S \
    > "results/pggb/logs/${CHR}.log" 2>&1
done
```

`scripts/pg02_run_pggb.sh` は、これに前後の確認(入力の存在チェック、完了済みのスキップ、成否のまとめ)を足しただけのものです。**染色体名は `data/input/` にあるファイルから取るので、`chr01..chr09` は前提にしていません。**

```bash
bash scripts/pg02_run_pggb.sh . 32            # data/input にある全部
bash scripts/pg02_run_pggb.sh . 32 chr09      # chr09 だけ
```

**SLURM で並列に投げる**場合は、染色体ごとに独立したジョブにします(job array にしないのは、染色体ごとに必要なメモリが大きく違うためです。§5.7.5 参照)。

`run_pggb_chr.sbatch`:

```bash
#!/bin/bash
#SBATCH --job-name=pggb
#SBATCH --output=results/pggb/logs/%x_%j.out
#SBATCH --error=results/pggb/logs/%x_%j.err
#SBATCH --mem=64G
#SBATCH --cpus-per-task=32
#SBATCH --time=12:00:00

CHR=${1:?染色体名 (例 chr01) を指定してください}
SIF=${2:?pggb_latest.sif の絶対パスを指定してください}

singularity exec "$SIF" pggb \
  -i "data/input/${CHR}.fa.gz" \
  -o "results/pggb/${CHR}" \
  -p 95 -s 10000 -n 6 -B 1G -V CUN#1 -t 32 -m -S
```

```bash
mkdir -p results/pggb/logs
for CHR in $(ls data/input/*.fa.gz | xargs -n1 basename | sed 's/\.fa\.gz$//'); do
  sbatch run_pggb_chr.sbatch "$CHR" "$SIF"
done
```

実時間で 8-12 時間で完了する見込みです。

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
results/pggb/chr09/
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

## 5.9 アラインメントを使い回す —— impg

ここまでの手順は **PGGB を単体で回す**やり方です。一方、現在の大規模 pangenome 構築(HPRC の release 2 など)は、**PGGB を [impg](https://github.com/pangenome/impg) と組み合わせた流れ**になっています。本教材では実行しませんが、考え方は知っておく価値があります。

### 何をするものか

impg は **all-vs-all のアラインメント (PAF) そのものを「暗黙の pangenome graph」として扱う**ツールです。graph を実体化せずに、アラインメントの上を直接たどります。

```bash
# 1. アラインメントに索引を張る
impg index -a aln.paf -i aln.impg

# 2. ある区間に相同な配列を、全ハプロタイプから取り出す
impg query -a aln.paf -r CUN#1#chr09:1000000-1200000 -d 100 -x

# 3. アラインメント網から窓を切る(1 Mb 窓、100 kb まで隙間を吸収)
impg partition -a aln.paf -w 1000000 -d 100000

# 4. 窓ごとに作ったグラフを 1 本に綴じる
impg lace -l gfa_list.txt -o combined.gfa
```

HPRC v2 の構築は、おおまかに **`impg partition` で全体を窓に分割 → 各窓を PGGB で構築 → `impg lace` で綴じる**、という流れです。

### 利点は「規模」だけではない

「大きいデータだから分割する」と思われがちですが、**本質は「アラインメントとグラフ誘導を分離すること」**にあります。6 haploid でも効く利点が 2 つあります。

**1. 最も高価な工程をやり直さずに済む**

PGGB の 3 段のうち、圧倒的に重いのは wfmash のアラインメントです。impg は **PAF を一次オブジェクトとして保存・索引化**するので、**アラインメントを 1 回だけ回し、seqwish や smoothxg のパラメータを変えたグラフを何度でも作り直せます**。

§5.5.3 の `-B` を思い出してください。指定を忘れて圧縮不足のグラフができたとき、PGGB 単体なら**アラインメントから丸ごとやり直し**です。アラインメントが別の成果物として残っていれば、払うコストはグラフ誘導のやり直しだけになります。§5.5.2 で「`-s` を振って比べる」と書いた作業も、本来はこの形のほうが安く済みます。

**2. 座位単位で切り出せる**

`impg query` は、**ある区間に相同な配列を全ハプロタイプから引き出します**。着目している遺伝子や QTL 区間だけの小さな graph を作って、そこでパラメータを詰める、といったことができます。染色体全体の graph を保持する必要がありません。

**3. 染色体の対応を前提にしない**

`impg partition` の窓は、**アラインメント網の transitive homology から動的に決まります**。「chr01 は chr01 と対応する」という前提が要りません。

これは §1.4 で扱った「染色体別に作るか、全ゲノム一括か」という問いへの、3 つ目の答えになっています。染色体別は対応関係を仮定し、一括は巨大なリソースを要求する —— **相同性から窓を切れば、どちらの制約からも外れられます**。がんゲノムのように染色体間の組み換えがある対象では、この性質が規模と無関係に効いてきます。

### 本教材で使わない理由

6 haploid × 染色体 1 本という規模では、PGGB 単体で 30-60 分で終わります。工程を分ける利点より、**1 コマンドで完結する分かりやすさ**を取りました。ただし、**自分のデータで本格的にやるなら、まず impg を検討してください。**

---

## 5.10 この章のまとめ

- Pangenome 構築ツールは Minigraph、Minigraph-Cactus (MC)、PGGB の 3 択
- 「参照」には**①骨格 / ②座標アンカー / ③比較基準**の 3 つの意味がある。**PGGB が要求するのは ② だけ**
- MC を選べば参照座標系の固定・Giraffe 用インデックス・softmask 不問が得られる。**出したい成果物から逆算して選ぶ**
- PGGB は構築自体は reference-free。**VCF のために座標 path を 1 本だけ選ぶ**(`-V CUN#1`)
- PGGB は 3 段構成: wfmash → seqwish → smoothxg
- **パラメータは文献値を写すのではなく、自分のデータから決められる** —— `-p` は `mash triangle` や `wfmash -m` で同一性を測る、`-s` は 1 染色体で振って圧縮率を比べる、`-n` は自動検出の答え合わせ
- 特定バージョンでは `-B` の明示指定が必要。**コマンドラインではなく `.params.yml` で事後確認する**
- HPC では Singularity/Apptainer。Docker は root ではなく `docker` グループ、Podman は rootless
- **impg はアラインメントとグラフ誘導を分離する** —— 高価な wfmash をやり直さずにパラメータを振れ、座位単位で切り出せ、染色体の対応を仮定しない(§5.9)

次章では、生成された graph の**品質を評価**します。

---

## 参考文献

- Garrison E, et al. (2024). Building pangenome graphs. *bioRxiv*.
- Li H (2020). The design and construction of reference pangenome graphs with minigraph. *Genome Biol* 21:265.
- Hickey G, et al. (2024). Pangenome graph construction from genome alignments with Minigraph-Cactus. *Nat Biotechnol* 42:663-673.
- Kiryu Y, et al. (2026). AlleleMiner. *DNA Res* 33:dsag004.

### ツール・仕様

- **pggb**: <https://github.com/pangenome/pggb> / ドキュメント <https://pggb.readthedocs.io/>
- **wfmash**: <https://github.com/waveygang/wfmash>
- **impg** (implicit pangenome graph): <https://github.com/pangenome/impg>
- **HPRC release 2 の構築手順**: <https://github.com/pangenome/HPRCv2>
- **Minigraph-Cactus** の手順: <https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md>
  (本教材向けの手順は [付録](appendix_minigraph_cactus.md))

---

[← 第4章](04_data_qc.md) | [第6章: グラフの QC →](06_graph_qc.md)

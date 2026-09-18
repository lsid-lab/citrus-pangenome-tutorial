# 第4章 データの品質確認

## 学ぶこと

- なぜ pangenome 構築前に QC が必要か
- `seqkit stats` の基本的な使い方
- **QC の判定を「絶対値の閾値」ではなく「その種の期待値との比較」で行う考え方**
- 柑橘における各指標の期待値と、他の生物種に持っていくときの読み替え方

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

本教材で使うツールの多くは **conda / mamba** で入ります。まだ環境がない場合は、**Miniforge** の導入が現在の推奨です。

**なぜ Miniforge か**: conda-forge チャンネルが既定で設定済みで、**高速な依存解決を行う `mamba` が最初から同梱**されています。Anaconda / Miniconda の既定チャンネル (`defaults`) は組織規模によって商用利用のライセンス条件に注意が必要ですが、Miniforge は conda-forge しか見ないのでその心配がありません。なお **Mambaforge は 2025年1月で廃止**され、Miniforge3 に統合されました(Miniforge 23.3.1 以降、両者は実質同一です)。

```bash
# Linux / macOS 共通(OS とアーキテクチャは自動で決まります)
wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash "Miniforge3-$(uname)-$(uname -m).sh"

# インストール後、シェルを開き直してから確認
mamba --version
```

> `$(uname)` は Linux なら `Linux`、macOS なら `Darwin` に、`$(uname -m)` は `x86_64` / `aarch64` / `arm64` に展開されます。`wget` がない環境では
> `curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"` でも同じです。

教材用の環境を分けておくと、後の章で追加ツールを入れるときに楽です:

```bash
mamba create -n citrus-pg -c conda-forge -c bioconda seqkit
mamba activate citrus-pg
```

**seqkit だけなら単一バイナリでも入ります**(conda を使いたくない場合):

```bash
wget https://github.com/shenwei356/seqkit/releases/download/v2.10.0/seqkit_linux_amd64.tar.gz
mkdir -p bin && tar xzf seqkit_linux_amd64.tar.gz && mv seqkit bin/
export PATH="$PWD/bin:$PATH"
```

> 第5章以降で使う PGGB 一式は conda ではなく **Singularity イメージ**で導入します(§5.5、`scripts/setup_singularity_wrappers.sh`)。conda で入れるのは、この章の QC ツールまでです。

### 実行

```bash
seqkit stats -a data/raw/CUN/CUNphKi_r1.0.pmol.fasta.gz
```

`-a` (all) オプションで N50 などの詳細も出ます。

### 出力例(温州の紀州由来ハプロタイプ = `CUN#2` の場合)

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

| sample | hap | num_seqs | sum_len_Mb | N50_Mb | GC% | 判定 | メモ |
|---|---|---|---|---|---|---|---|
| CUN | 1 | 9 | 348.5 | 34.5 | 35.9 | ✅ PASS | 親 2 品種より 25-55 Mb 大きい |
| CUN | 2 | 9 | 357.6 | 44.2 | 35.0 | ✅ PASS | 同上 |
| CKI | 1 | 9 | 304.2 | 33.5 | 36.0 | ✅ PASS | |
| CKI | 2 | 9 | 310.3 | 32.5 | 36.0 | ✅ PASS | |
| CKU | 1 | 9 | 323.9 | 36.6 | 35.9 | ✅ PASS | |
| CKU | 2 | 9 | 303.4 | 32.9 | 36.0 | ✅ PASS | |

**6 haploid すべてが PASS** です。`qc01_stats.sh` が使う柑橘の期待範囲(290-370 Mb、9 染色体、GC 34-38%)を全員が満たしています。

CUN の 2 hap が親 2 品種より大きい点は、**温州として公表されているサイズの範囲内**です(Shimizu et al. 2017: 359.7 Mb、Kawahara et al. 2020: 346 Mb)。親 2 品種はそもそも別種なので、ゲノムサイズが 10-20% 違うのは珍しくありません。

ただし、**path が長いほど Jaccard 類似度は下がる**という性質があるため、この差は第6章の解釈に効いてきます(§6.4.4)。観察として憶えておいてください。

---

## 4.5 見るべきポイント —— 閾値ではなく「期待値との比較」

ここまでに出てきた数字は、**すべて「柑橘 (n=9、haploid 約 300-360 Mb) の場合」の期待値**です。絶対値の閾値として覚えてしまうと、他の生物種に持っていった途端に機能しなくなります。

たとえば「`num_seqs ≤ 20` なら染色体スケール」をヒトに当てはめると、**T2T-CHM13(24 配列)ですら「半断片化」判定**になってしまいます。判定の形を、**「その種で期待される値」と比べる**という考え方に置き換えましょう。

> 本教材が使う 6 つのアセンブリは、いずれも**染色体スケールまで仕上がった公開済みの高品質アセンブリ**です(9 本の pseudomolecule、N50 30 Mb 超)。ここで QC をやり直して問題を見つける必要は、実は**ありません**。
>
> この節の狙いは、**あなたが自分のデータで pangenome を作るときに、何を、何と比べて確認すべきか**を持ち帰ってもらうことです。

### ポイント1: 連続性 —— 染色体スケールに到達しているか

見るべきは配列数の絶対値ではなく、**その種で期待される染色体数 n と比べてどうか**です。

| 指標 | 判定の考え方 |
|---|---|
| `num_seqs` | **n(+ オルガネラ + unplaced)に近いか** |
| `N50` | **その種の染色体長の中央値に近いか**。染色体スケールなら N50 ≒ 中央染色体長になる |
| `L50` | **染色体数の約半分**(n/2 前後)か。これが n を大きく超えるなら断片化 |

| 種 | n (haploid) | 期待 num_seqs | 期待 N50 |
|---|---|---|---|
| 柑橘 (*Citrus*) | 9 | 9 前後 | 30-40 Mb |
| ヒト (T2T-CHM13) | 23 | 24 前後 | ~150 Mb |
| イネ | 12 | 12 前後 | ~30 Mb |
| シロイヌナズナ | 5 | 5 前後 | ~23 Mb |

**`num_seqs` 単独で判定しないでください。** unplaced scaffold を数千本抱えていても、**総塩基の 95% 以上が上位 n 本に入っていれば実用上は染色体スケール**です。逆に配列数が少なくても、中身が `N` だらけということもあります。`num_seqs` と `N50` は必ずセットで見ます。

同じ温州みかんでも、アセンブリによってここまで違います:

| アセンブリ | 総塩基 | 配列数 | N50 |
|---|---:|---:|---:|
| Shimizu et al. (2017) draft | 359.7 Mb | 20,876 scaffold | 386 kb |
| Isobe et al. (2023) `CUNphKi` | 348.5 Mb | 9 pseudomolecule | ~38 Mb |

**ゲノムサイズはほぼ同じで、連続性だけが 100 倍違う**、という関係です。Pangenome 構築に使えるのは後者です。

### ポイント2: 総塩基数 —— 既知のゲノムサイズと比べる

`sum_len` に普遍的な「正常範囲」はありません。**対象種のゲノムサイズ推定値**と比べてください。推定値の入手先:

- 同種・近縁種の先行論文のアセンブリサイズ
- フローサイトメトリーによる実測値([Plant DNA C-values Database](https://cvalues.science.kew.org/) など)
- 自分のリードからの k-mer ベース推定(GenomeScope2 など)

|推定値との関係|考えられる解釈|
|---|---|
|±10% 以内|想定どおり|
|**小さい**|情報欠損、あるいは反復配列の collapse|
|**大きい**|haplotype の重複 (leakage)、コンタミ、あるいは推定値のほうが古い|

**ただし、大きい/小さいのどちらも「それだけでは異常の証拠になりません」。** 近縁種間でもゲノムサイズは普通に 10-20% 違いますし、推定値自体にも幅があります。§4.4 の CUN がまさにその例で、**親 2 品種より大きいものの、温州として公表されている値の範囲内**です。

### ポイント3: hap1 と hap2 の対比 —— これは種によらず使える

同一個体の 2 haplotype は**同じ個体の同じ染色体セット**なので、サイズが大きく食い違う理由は基本的にありません。これは数少ない、**種に依存しないチェック**です。

- 差 ≤ 5%: 対応が取れている
- 差 5-15%: 要注意
- 差 > 15%: phasing の失敗や、片方への配列の偏りを疑う

(ヘテロ接合度が極端に高い個体、大きな半接合領域、性染色体を含む場合はこの限りではありません。)

### ポイント4: GC 含量 —— 同種の既知値と比べる

GC 含量も種ごとに大きく違います —— 柑橘 ~35%、ヒト ~41%、イネ ~44%、*Plasmodium falciparum* に至っては ~19% です。**絶対値ではなく、同種の先行アセンブリと比べて 1-2 pp 以内に収まっているか**を見てください。

全体値だけでなく、**配列ごとの GC** を見るのが実践的です:

```bash
# 配列ごとに 名前 / 長さ / GC% を出す
seqkit fx2tab -nlg data/raw/CUN/CUNphKi_r1.0.pmol.fasta.gz
```

**1 本だけ極端に外れた配列**があれば、細菌・オルガネラ・別種のコンタミネーションを疑います(オルガネラゲノムは核ゲノムと GC が異なるので、この方法で見つかります)。

### 参考: 別の種で `qc01_stats.sh` を使うには

`scripts/qc01_stats.sh` の期待値は環境変数で上書きできます。柑橘以外のデータに使うときは、その種の値を渡してください:

```bash
# 例: ヒト (n=23、約 3.1 Gb、GC ~41%)
EXPECTED_NUM_CHR=23 \
EXPECTED_SIZE_MIN=2900 EXPECTED_SIZE_MAX=3300 \
EXPECTED_GC_MIN=40.0 EXPECTED_GC_MAX=42.0 \
bash scripts/qc01_stats.sh tables/samplesheet.tsv .
```

染色体スケール判定の閾値は `EXPECTED_NUM_CHR` から自動で導かれます。

---

## 4.6 追加の QC ツール(オプション)

### BUSCO / compleasm(完全性)

系統に近い遺伝子セット(orthologs)がどれだけ検出されるかで、**アセンブリの完全性**を評価します。指標の読み方と lineage dataset の選び方は §1.3 で扱いました。

**インストール**(§4.3 の mamba を使います):

```bash
# BUSCO 本体、または高速な再実装である compleasm
mamba install -c conda-forge -c bioconda busco
mamba install -c conda-forge -c bioconda compleasm
```

**lineage dataset の取得**(初回のみ。ダウンロードに数分かかります):

```bash
compleasm download eudicots_odb10
```

**実行**:

```bash
# compleasm (BUSCO より高速)
compleasm run -a data/raw/CUN/CUNphKi_r1.0.pmol.fasta.gz \
              -o busco/CUN_hap2 \
              -l eudicots_odb10 \
              -t 16
```

期待値: **Complete > 95%、Duplicated < 5%**。haplotype-resolved アセンブリなので、**6 haploid それぞれ別々に**実行します。

### Merqury(k-mer QV と haplotype leakage)

HiFi 生リードから作った k-mer データベースと、アセンブリを比較して:
- **QV (base accuracy)**: 塩基精度、期待 > 40
- **完全性**: reads が捉えた k-mer のうち、アセンブリに含まれる割合
- **重複度**: haplotype leakage の検出

**インストール**: Merqury は Java と R に依存するので、**専用の環境に分ける**のが安全です。`meryl` は依存として一緒に入ります。

```bash
mamba create -n merqury -c conda-forge -c bioconda merqury openjdk=11
mamba activate merqury

merqury.sh --version   # 動作確認
meryl --version
```

> `$MERQURY is not set` のようなエラーが出た場合は、`export MERQURY="$CONDA_PREFIX/share/merqury"` を設定してから再実行してください。

**実行**:

```bash
# 1) HiFi 生リードから k-mer データベースを作る
meryl count k=21 output hifi.meryl hifi_reads.fastq.gz

# 2) アセンブリと照合(2 hap を同時に渡すと hap 間の比較も出る)
merqury.sh hifi.meryl CUNphKu_r1.0.fasta.gz CUNphKi_r1.0.fasta.gz CUN
```

trio のリード(両親 + 子)が揃っていれば、`hapmers` を作ることで **false duplication rate** や **hap-mer blob plot** まで出せます。これが haplotype leakage を判定する本来の道具です。

**本教材では HiFi 生リードを扱わないため、Merqury は実行しません。** 生リードは DDBJ DRA の BioProject `PRJDB15866`(§3.1)から取得できるので、興味があれば自分で試してみてください。

---

## 4.7 この章のまとめ

- `seqkit stats` で基本統計を取得する
- **QC の基準は絶対値の閾値ではなく「その種で期待される値」との比較**。柑橘の 300-360 Mb / 9 染色体 / GC 34-38% は、あくまで柑橘の値
- `num_seqs` 単独では判定しない。`N50` とセットで見る(同じ温州でも 20,876 scaffold のアセンブリと 9 pseudomolecule のアセンブリがある)
- hap1 と hap2 のサイズ比較は、数少ない**種に依存しない**チェック
- **CUN の 2 hap は他より大きいが、温州として公表されているサイズ (346-360 Mb) の範囲内**。「手元の他サンプルより大きい = 異常」と短絡しない
- ただしこのサイズ差は **第6章の Jaccard 類似度に効く**(§6.4.4)ので、観察として憶えておく

次章では、いよいよ **pangenome graph の構築**に入ります。

---

## 参考文献

- Shen W, et al. (2024). SeqKit2. *iMeta* 3:e191.
- Rhie A, et al. (2020). Merqury: reference-free quality, completeness, and phasing assessment for genome assemblies. *Genome Biol* 21:245.
- Shimizu T, et al. (2017). Draft sequencing of the heterozygous diploid genome of satsuma (*Citrus unshiu* Marc.) using a hybrid assembly approach. *Front Genet* 8:180. (総長 359.7 Mb、20,876 scaffold、N50 386 kb)
- Kawahara Y, et al. (2020). Mikan Genome Database (MiGD). *Breed Sci* 70(2). (温州 346 Mb)
- Manni M, et al. (2021). BUSCO update. *Mol Biol Evol* 38:4647-4654.
- Huang N, Li H (2023). compleasm: a faster and more accurate reimplementation of BUSCO. *Bioinformatics* 39:btad595.
- **Miniforge**: <https://github.com/conda-forge/miniforge>
- **Plant DNA C-values Database**: <https://cvalues.science.kew.org/>

---

[← 第3章](03_data_acquisition.md) | [第5章: Pangenome graph 構築 →](05_graph_construction.md)

# 第3章 データの取得

## 学ぶこと

- 論文からデータの所在を追跡する方法
- **BioProject 記載データ vs 実際にアセンブリがある場所**の違い
- Plant GARDEN からのアセンブリダウンロード
- プロジェクト用のディレクトリ構成

---

## 3.1 データを探す:論文の "Data availability" セクション

論文を再現するときの最初の一歩は、**"Data availability" セクション**を読むことです。

本教材が使うデータの出典は次の論文です:

> Isobe S, et al. (2023). Haploid-resolved and chromosome-scale genome assembly in *Citrus unshiu* and its parental species, *C. nobilis* and *C. kinokuni*. *bioRxiv* 2023.06.02.543356.
> <https://www.biorxiv.org/content/10.1101/2023.06.02.543356v1>

この論文の "Data availability" の該当部分を見ると:

> The sequence reads are available from the DNA Data Bank of Japan (DDBJ) Sequence Read Archive (DRA) under the BioProject number **PRJDB15866**. The assembled scaffold sequences, gene sequences, and annotation files are available at **Plant GARDEN**.

ここで**極めて重要な情報**が2つあります:

1. **生リード** (raw HiFi reads) は **DDBJ DRA** (BioProject PRJDB15866)
2. **アセンブリ** (scaffold sequences)、遺伝子、アノテーション は **Plant GARDEN**

**アセンブリを使いたい場合は Plant GARDEN が正解**です。DDBJ にはアセンブリはありません。

---

## 3.2 一般論:シーケンスデータのアーカイブ

生命科学のシーケンスデータは、以下の**国際的なアーカイブ**に登録されます:

| アーカイブ | 運営 | 国 | 主な内容 |
|---|---|---|---|
| **NCBI SRA** | NCBI | 米国 | 生リード、アセンブリ |
| **EBI ENA** | EMBL-EBI | 欧州 | 生リード、アセンブリ |
| **DDBJ DRA** | 情報・システム研究機構 | 日本 | 生リード、アセンブリ |
| **Plant GARDEN** | かずさDNA研究所 | 日本 | 植物ゲノム(アセンブリ + アノテーション) |

**INSDC (International Nucleotide Sequence Database Collaboration)** という国際協定により、NCBI SRA / EBI ENA / DDBJ DRA は**相互ミラーリング**されています。片方に登録すれば、他の2つでも参照できます。

しかし **Plant GARDEN は INSDC には含まれず、独自の登録先**です。植物ゲノムの完成アセンブリの多くはここに集約されています。

---

## 3.3 BioProject と BioSample の関係

DDBJ / NCBI / ENA での登録は、階層構造になっています:

```
BioProject (プロジェクト全体)
  │
  ├── BioSample #1 (個体1、品種1など)
  │     ├── Experiment #1 (シーケンス実験1)
  │     │     └── Run #1 (実際のリードデータ)
  │     └── Experiment #2
  │           └── Run #2
  ├── BioSample #2
  │     ├── ...
```

`PRJDB15866` の場合:
- BioProject 全体 = Isobe 2023 温州+紀州+九年母のプロジェクト
- BioSample が 3 つ:
  - SAMD00608599 (温州、Miyagawa Wase)
  - SAMD00608600 (九年母、Kunenbo Kagoshima)
  - SAMD00608601 (紀州、Kishu mikan)
- 各 BioSample に**HiFi WGS**、**短鎖 WGS**、**Hi-C**、**Iso-Seq**の Run が紐づく

生リードから自分でアセンブリしたい場合は、これらの Run accession から FASTQ をダウンロードします(この教材ではやりません)。

---

## 3.4 Plant GARDEN からのアセンブリ取得

### 3.4.1 Plant GARDEN とは

**Plant GARDEN** (https://plantgarden.jp) は、かずさDNA研究所が運営する**植物ゲノムのポータル**です。

- 数十種の植物ゲノムアセンブリを公開
- 遺伝子予測、アノテーション、K-mer データベースも一体で提供
- 学術用途は無償で利用可能

温州みかんの解析対象は以下のパスにあります:

- 温州: `https://plantgarden.jp/ja/download/Citrus_unshiu/`
- 紀州: `https://plantgarden.jp/ja/download/Citrus_kinokuni/`
- 九年母: `https://plantgarden.jp/ja/download/Citrus_nobilis/`

### 3.4.2 Assembly ID (Gxxx) の読み方

Plant GARDEN では、同一種でも**複数のアセンブリバージョン**が `G001`, `G002`, … という **Assembly ID** で並列に公開されています。

Isobe et al. (2023) は 3 品種それぞれについて、次の **3 本セット**を登録しています。ID はこの順に振られています。

1. **unphased** —— 2 倍体を 1 セットに畳んだアセンブリ
2. **hap1** —— phased ハプロタイプ 1
3. **hap2** —— phased ハプロタイプ 2

**温州みかん (CUN, taxon `t55188`)**

| Assembly ID | ファイル | 内容 |
|---|---|---|
| t55188.**G001** | `C_unshiu_v1.0_scaffolds.fa.gz` | Kawahara et al. (2020) の**別研究**のアセンブリ(旧版、unphased) |
| t55188.**G002** | `CUNuph_r1.0.fasta.gz` | Isobe 2023 の **unphased**(`uph` = **u**n**ph**ased) |
| t55188.**G003** | `CUNphKi_r1.0.pmol.fasta.gz` | **hap1**(`ph` = **ph**ased、`Ki` = **Ki**shu 由来)← 本教材で使用 |
| t55188.**G004** | `CUNphKu_r1.0.ch1-9.fasta.gz` | **hap2**(`Ku` = **Ku**nenbo 由来)← 本教材で使用 |

**紀州みかん (CKI, `t408488`) / 九年母 (CKU, `t481549`)** —— こちらは `G003` までしかありません。

| Assembly ID | 内容 | ファイル |
|---|---|---|
| **G001** | **unphased**(本教材では未使用) | —— |
| **G002** | **hap1** ← 本教材で使用 | `CKIhap1_r1.0.pmol.fasta.gz` / `CKUhap1_r1.0.pmol.fasta.gz` |
| **G003** | **hap2** ← 本教材で使用 | `CKIhap2_r1.0.pmol.fasta.gz` / `CKUhap2_r1.0.pmol.fasta.gz` |

> **なぜ温州だけ番号が 1 つずれるのか**: 温州には Isobe 2023 より前に Kawahara et al. (2020) のアセンブリが `G001` として既に登録されていたためです。つまり「Isobe 2023 の unphased / hap1 / hap2」は、**温州では G002 / G003 / G004**、**紀州・九年母では G001 / G002 / G003** に対応します。番号を丸暗記するのではなく、**unphased → hap1 → hap2 の並び**で覚えてください。

**本教材で使う 6 ファイルは、3 品種それぞれの hap1 と hap2** です(温州 = G003/G004、紀州・九年母 = G002/G003)。

unphased アセンブリ(温州 G002、親 2 品種の G001)は使いません。**2 つのハプロタイプが 1 本に畳み込まれていて、haplotype 間の違いを graph に表現できない**からです。逆に言えば、pangenome を作る意味があるのは phased アセンブリが揃っているからこそ、ということでもあります。

### 3.4.3 `pmol` とは何か

ファイル名に頻出する **`pmol`** は **pseudomolecule(疑似分子 / 疑似染色体配列)** の略です。

「pseudo(疑似)」と付くのは、**1 本の染色体を端から端まで 1 分子として読み切ったもの**ではないからです。実際には次の手順で作られます:

1. 配列決定とアセンブリで、多数の contig / scaffold ができる
2. それを Hi-C や連鎖地図などの情報で、**染色体上の順序と向きに並べる**
3. 並べたものを(隙間を `N` で埋めるなどして)**染色体 1 本を代表する 1 本の配列**に仕立てる

つまり `pmol` は「染色体を代表するように**構築された**配列」です。実務上は染色体配列として扱って差し支えありませんが、**実測された 1 分子ではない**ことは頭の隅に置いてください —— ギャップ(`N` の連続)が残っていたり、順序・向きの取り違えが含まれている可能性があります。第4章の QC で `N` の割合を確認するのは、このためでもあります。

この命名は Plant GARDEN / かずさDNA研究所の系統でかなり一貫しており、トマトの `SLM_r2.0.pmol`、トウガラシの `CAN_r1.2.pmol`、アジサイの `HMA_r1.2.pmol` のように、他種でも同じ流儀で付けられています。**他の Plant GARDEN データを触るときも、`pmol` が付いたファイルを選べば染色体レベルの配列が得られる**と覚えておくと便利です。

> なお、温州の hap2 だけは `.pmol` ではなく **`CUNphKu_r1.0.ch1-9.fasta.gz`**(染色体 1〜9 のみを収めたファイル)という名前です。本教材は染色体別に graph を作るので、どちらの命名でも支障はありません。

---

## 3.5 プロジェクトのディレクトリ構成

本教材は **GitHub リポジトリからの clone** を前提としています。

```bash
git clone https://github.com/lsid-lab/citrus-pangenome-tutorial.git
cd citrus-pangenome-tutorial
```

Clone した時点で、以下のディレクトリと資料が揃っています:

```
citrus-pangenome-tutorial/
├── README.md
├── docs/ja/            # 本チュートリアル本体 (このファイル群)
├── scripts/            # 実行スクリプト一式
├── tables/             # サンプルシートなどメタデータ
│   └── samplesheet.tsv
└── LICENSE
```

**追加で作成するのは `data/` ディレクトリだけ**です。ここに各品種のアセンブリを配置します:

```bash
mkdir -p data/{CUN,CKI,CKU}
```

Pangenome の中間ファイルや出力(`03_pangenome/` など)は、後続のスクリプトが自動で作成します。

推奨する構成:

```
citrus-pangenome-tutorial/
├── ...(clone時のファイル)
├── data/                          # ★ ダウンロードしたアセンブリ(あなたが用意)
│   ├── CUN/
│   │   ├── CUNphKi_r1.0.pmol.fasta.gz
│   │   └── CUNphKu_r1.0.ch1-9.fasta.gz
│   ├── CKI/
│   │   ├── CKIhap1_r1.0.pmol.fasta.gz
│   │   └── CKIhap2_r1.0.pmol.fasta.gz
│   └── CKU/
│       ├── CKUhap1_r1.0.pmol.fasta.gz
│       └── CKUhap2_r1.0.pmol.fasta.gz
└── 03_pangenome/                  # ★ 後続スクリプトが作成
    ├── by_chr/
    ├── qc/
    └── logs/
```

---

## 3.6 ダウンロード実践

### 3.6.1 手動ダウンロード

各ファイルは Plant GARDEN のブラウザ経由でダウンロードできますが、**サイズが 100 MB 前後**なのでコマンドラインでの取得を推奨します。

### 3.6.2 スクリプトでの一括取得

`scripts/download_plantgarden.sh` を使います(§3.5 でリポジトリのルートに `cd` した状態から続けます):

```bash
bash scripts/download_plantgarden.sh data/
```

このスクリプトは以下を実行します:
1. Plant GARDEN の対象URLから 6 つのアセンブリファイルをダウンロード
2. `data/{CUN,CKI,CKU}/` に配置
3. sha256 チェックサムで整合性検証(オプション)

### 3.6.3 期待されるファイルサイズ

| ファイル | サイズ (bgzip) |
|---|---:|
| CUNphKi_r1.0.pmol.fasta.gz | ~101 MB |
| CUNphKu_r1.0.ch1-9.fasta.gz | ~104 MB |
| CKIhap1_r1.0.pmol.fasta.gz | ~88 MB |
| CKIhap2_r1.0.pmol.fasta.gz | ~91 MB |
| CKUhap1_r1.0.pmol.fasta.gz | ~94 MB |
| CKUhap2_r1.0.pmol.fasta.gz | ~88 MB |
| **計** | **~570 MB** |

---

## 3.7 サンプルシートの確認

後続の全スクリプトは、`tables/samplesheet.tsv` を入力として参照します。**このファイルは clone 時に既に用意されています**が、パスが実環境と一致するかは確認してください。

`tables/samplesheet.tsv` の内容:

```tsv
sample_id	cultivar_jp	species	pansn_prefix	hap1_path	hap2_path	source	pedigree
CUN	温州みかん	Citrus unshiu	CUN	data/CUN/CUNphKi_r1.0.pmol.fasta.gz	data/CUN/CUNphKu_r1.0.ch1-9.fasta.gz	Plant GARDEN t55188.G003/G004	F1: CKI × CKU
CKI	紀州みかん	Citrus kinokuni	CKI	data/CKI/CKIhap1_r1.0.pmol.fasta.gz	data/CKI/CKIhap2_r1.0.pmol.fasta.gz	Plant GARDEN t408488.G002/G003	温州の母親
CKU	九年母	Citrus nobilis	CKU	data/CKU/CKUhap1_r1.0.pmol.fasta.gz	data/CKU/CKUhap2_r1.0.pmol.fasta.gz	Plant GARDEN t481549.G002/G003	温州の父親
```

**列の意味**:
- `sample_id`: 内部用の略号(3文字)
- `cultivar_jp`: 表示用の和名
- `species`: 学名
- `pansn_prefix`: 後で PanSN 名の一部として使う短い ID
- `hap1_path` / `hap2_path`: FASTA ファイルへの相対パス(リポジトリルートから)
- `source`: 出典
- `pedigree`: 家族関係のメモ

このファイル1つを正しく維持することが、後続の全ステップの前提です。

---

## 3.8 ダウンロード後の確認

すべてダウンロードが完了したら、以下で確認します:

```bash
# ファイルサイズと有無を確認
ls -la data/*/*.fasta.gz

# 各 FASTA の中身をちらっと見る
zcat data/CUN/CUNphKi_r1.0.pmol.fasta.gz | head -3

# 期待される出力例:
# >CUNphKi_r1.0ch1
# GCTAGCTAGCTAGCTAGCTAGCTAGC...
# TAGCTAGCTAGCTAGCTAGCTAGCTA...
```

**確認ポイント**:
- 全 6 ファイルが正しくダウンロードされた
- ファイルサイズが期待通り(数十〜100+ MB)
- FASTA の最初の行に `>` で始まる配列名がある
- 配列名は `CUNphKi_r1.0ch1` のような形式

---

## 3.9 章のまとめ

- 論文の "Data availability" セクションが最初の手がかり
- **BioProject 番号 = DDBJ の生リード**、Plant GARDEN = アセンブリ、と使い分ける
- Plant GARDEN の命名規則:`C{種略号}{ph or hap}{Ki/Ku or 1/2}`
- **Assembly ID (Gxxx) は unphased → hap1 → hap2 の順**。温州だけ旧アセンブリが `G001` を占めるため 1 つずれる
- **`pmol` = pseudomolecule**。contig を染色体上に並べて構築した「染色体を代表する配列」であって、読み切った 1 分子ではない
- pangenome の入力に使うのは **phased (hap1/hap2)** のみ。unphased は haplotype 間の違いを表現できない
- サンプルシート(TSV)で 3 品種 × 2 hap = 6 haploid を一元管理

次章では、ダウンロードしたアセンブリの**品質を確認**します。

---

## 参考文献

- **Plant GARDEN**: https://plantgarden.jp
- **DDBJ**: https://www.ddbj.nig.ac.jp/
- **NCBI**: https://www.ncbi.nlm.nih.gov/
- **INSDC**: http://www.insdc.org/
- Isobe S, et al. (2023). Haploid-resolved and chromosome-scale genome assembly in *Citrus unshiu* and its parental species, *C. nobilis* and *C. kinokuni*. *bioRxiv* 2023.06.02.543356. <https://www.biorxiv.org/content/10.1101/2023.06.02.543356v1>
- Kawahara Y, et al. (2020). Mikan Genome Database (MiGD): integrated database of genome annotation, genomic diversity, and CAPS marker information for mandarin molecular breeding. *Breed Sci* 70(2). (温州の `G001` アセンブリの出典)

---

[← 第2章](02_target_organisms.md) | [第4章: データの品質確認 →](04_data_qc.md)

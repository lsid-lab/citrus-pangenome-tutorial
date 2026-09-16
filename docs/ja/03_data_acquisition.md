# 第3章 データの取得

## 学ぶこと

- 論文からデータの所在を追跡する方法
- **BioProject 記載データ vs 実際にアセンブリがある場所**の違い
- Plant GARDEN からのアセンブリダウンロード
- プロジェクト用のディレクトリ構成

---

## 3.1 データを探す:論文の "Data availability" セクション

論文を再現するときの最初の一歩は、**"Data availability" セクション**を読むことです。

Isobe et al. (2023) *bioRxiv* の該当部分を見ると:

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

### 3.4.2 ファイル命名規則

Plant GARDEN では、同一種でも**複数のアセンブリバージョン**が並列で公開されています。温州みかんの場合:

| Assembly ID | ファイル | 内容 |
|---|---|---|
| t55188.**G001** | `C_unshiu_v1.0_scaffolds.fa.gz` | Kawahara 2020 hybrid(旧、非phased) |
| t55188.**G002** | `CUNuph_r1.0.fasta.gz` | Isobe 2023 phased(統合表示) |
| t55188.**G003** | `CUNphKi_r1.0.pmol.fasta.gz` | **STS の hap1 (Kishu由来)** |
| t55188.**G004** | `CUNphKu_r1.0.ch1-9.fasta.gz` | **STS の hap2 (Kunenbo由来)** |

**本教材で使うのは G003 と G004**(および紀州・九年母の対応する phased アセンブリ)です。

---

## 3.5 プロジェクトのディレクトリ構成

本教材は **GitHub リポジトリからの clone** を前提としています。

```bash
git clone https://github.com/<user>/citrus-pangenome-tutorial.git
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
mkdir -p data/{satsuma,kishu,kunenbo}
```

Pangenome の中間ファイルや出力(`03_pangenome/` など)は、後続のスクリプトが自動で作成します。

推奨する構成:

```
citrus-pangenome-tutorial/
├── ...(clone時のファイル)
├── data/                          # ★ ダウンロードしたアセンブリ(あなたが用意)
│   ├── satsuma/
│   │   ├── CUNphKi_r1.0.pmol.fasta.gz
│   │   └── CUNphKu_r1.0.ch1-9.fasta.gz
│   ├── kishu/
│   │   ├── CKIhap1_r1.0.pmol.fasta.gz
│   │   └── CKIhap2_r1.0.pmol.fasta.gz
│   └── kunenbo/
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

`scripts/download_plantgarden.sh` を使います:

```bash
cd citrus_pangenome
bash scripts/download_plantgarden.sh data/
```

このスクリプトは以下を実行します:
1. Plant GARDEN の対象URLから 6 つのアセンブリファイルをダウンロード
2. `data/{satsuma,kishu,kunenbo}/` に配置
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
STS	温州みかん	Citrus unshiu	satsuma	data/satsuma/CUNphKi_r1.0.pmol.fasta.gz	data/satsuma/CUNphKu_r1.0.ch1-9.fasta.gz	Plant GARDEN t55188.G003/G004	F1: KSH × KNN
KSH	紀州みかん	Citrus kinokuni	kishu	data/kishu/CKIhap1_r1.0.pmol.fasta.gz	data/kishu/CKIhap2_r1.0.pmol.fasta.gz	Plant GARDEN t408488	温州の母親
KNN	九年母	Citrus nobilis	kunenbo	data/kunenbo/CKUhap1_r1.0.pmol.fasta.gz	data/kunenbo/CKUhap2_r1.0.pmol.fasta.gz	Plant GARDEN t481549	温州の父親
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
ls -la data/*/*.fa.gz

# 各 FASTA の中身をちらっと見る
zcat data/satsuma/CUNphKi_r1.0.pmol.fasta.gz | head -3

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
- サンプルシート(TSV)で 3 品種 × 2 hap = 6 haploid を一元管理

次章では、ダウンロードしたアセンブリの**品質を確認**します。

---

## 参考文献

- **Plant GARDEN**: https://plantgarden.jp
- **DDBJ**: https://www.ddbj.nig.ac.jp/
- **NCBI**: https://www.ncbi.nlm.nih.gov/
- **INSDC**: http://www.insdc.org/
- Isobe S, et al. (2023). bioRxiv 2023.06.02.543356.

---

[← 第2章](02_target_organisms.md) | [第4章: データの品質確認 →](04_data_qc.md)

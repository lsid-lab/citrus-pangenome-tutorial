# Citrus Pangenome Tutorial / 柑橘Pangenomeチュートリアル

[**English**](#english) | [**日本語**](#日本語)

---

## English

A hands-on tutorial for constructing a pangenome graph from haplotype-resolved assemblies of Satsuma mandarin (*Citrus unshiu*) and its parental cultivars.

### Target audience

- Advanced undergraduates with a few months to a year of bioinformatics experience
- Comfortable with genome basics (FASTA, assembly) but new to pangenomes
- Familiar with Linux command line, conda, and git

### What makes this tutorial distinctive

- **Real data**: uses the Plant GARDEN assemblies from Isobe et al. (2023) to build an actual working pangenome graph
- **Family relationships**: Satsuma = Kishu × Kunenbo (F1). The **clear pedigree** can be directly verified on the resulting graph
- **Includes real-world pitfalls**: bugs and misinterpretations encountered in actual data analysis (PGGB's silent `-B` parameter parse issue, haplotype leakage, etc.) are integrated into the tutorial
- **Realistic scale**: 3 cultivars × 2 haplotypes = 6 haploids, runnable in a student environment

### Table of contents

| Chapter | Title | Content |
|---|---|---|
| [0](docs/en/00_what_is_pangenome.md) | What is a Pangenome? | Concepts, vs reference genomes, why now |
| [1](docs/en/01_general_workflow.md) | General workflow | Cross-tool concepts: PanSN, GFA, node/edge/path |
| [2](docs/en/02_target_organisms.md) | Target organisms | Satsuma/Kishu/Kunenbo relationships, trio phasing significance |
| [3](docs/en/03_data_acquisition.md) | Data acquisition | BioProjects, Plant GARDEN, hands-on downloading |
| [4](docs/en/04_data_qc.md) | Data QC | seqkit stats, expected values, anomaly detection |
| [5](docs/en/05_graph_construction.md) | Graph construction | Tool comparison (minigraph / MC / PGGB), running PGGB |
| [6](docs/en/06_graph_qc.md) | Graph QC | 4-level evaluation, compression ratio, transclose-batch check |
| [7](docs/en/07_graph_interpretation.md) | Graph interpretation | VCF, visualization, pedigree verification |
| [A](docs/en/appendix_minigraph_cactus.md) | *Appendix*: Minigraph-Cactus | The same data through MC (untested) |

### Repository layout

```
citrus-pangenome-tutorial/
├── README.md
├── LICENSE
├── docs/
│   ├── ja/                      # Japanese version
│   └── en/                      # English version
├── scripts/                     # All executable scripts (English comments)
│   ├── setup_singularity_wrappers.sh
│   ├── download_plantgarden.sh
│   ├── qc01_stats.sh
│   ├── qc02_visualize.py
│   ├── pg01_prepare_input.sh
│   ├── pg02_run_pggb.sh
│   ├── pg03_qc_graph.sh
│   └── pg04_visualize.sh
├── tables/
│   ├── samplesheet.tsv
│   └── assembly_provenance_v1_v2.tsv
│
│   # --- created as you work through the tutorial (all git-ignored) ---
├── data/
│   ├── raw/                     # downloaded assemblies, left untouched
│   └── input/                   # PanSN-normalised, one FASTA per chromosome
├── results/                     # qc/ pggb/ graph_qc/ viz/
└── bin/                         # Singularity wrappers (optional)
```

### Prerequisites

- Linux (Ubuntu 22.04+ or CentOS 7+)
- **HPC access** (recommended: 32+ CPUs, 128 GB RAM; 200 GB RAM comfortable for chr03)
- Singularity or Docker
- Python 3.10+
- 200+ GB free disk space

### Runtime estimates

| Step | Target | Time |
|---|---|---|
| Data download | 6 haploid FASTAs | Minutes to 30 min (network-dependent) |
| Data QC | 6 FASTAs | 10-30 min |
| Pangenome construction | 9 chromosomes × PGGB | 5-10 hours (32 CPU, parallel) |
| Graph QC | 9 chromosomes | 30 min - 1 hour |

### Quick start

```bash
# 1. Clone
git clone https://github.com/lsid-lab/citrus-pangenome-tutorial.git
cd citrus-pangenome-tutorial

# 2. Follow the tutorial
# English:  docs/en/00_what_is_pangenome.md
# Japanese: docs/ja/00_what_is_pangenome.md
```

### Citation and acknowledgments

- Assembly data used in this tutorial is from **Isobe et al. (2023)** *bioRxiv*
- The version-specific PGGB `-B` behavior was diagnosed and documented by collaborator **Shuto Machida** through systematic experiments
- The tutorial content was developed through discussions and dialogue-based iteration

### License

- Code: MIT License (see `LICENSE`)
- Tutorial text: CC BY 4.0

### Contributing

Issues and Pull Requests are welcome. Corrections, translation contributions, and improvements are all appreciated.

---

## 日本語

温州みかん(*Citrus unshiu*)とその両親系統(紀州みかん・九年母)の haplotype-resolved アセンブリを用いて、pangenome graph を構築する実践的チュートリアルです。

### 対象読者

- バイオインフォマティクスに触れて数ヶ月〜1年程度の学部3年生
- ゲノム解析の基礎(FASTA, アセンブリ)は知っているが、pangenome は未経験
- Linux コマンドライン、conda、git の基本操作ができる

### この教材の特徴

- **実データを使う**: Isobe et al. (2023) が公開している Plant GARDEN のアセンブリを使い、実際に動く pangenome graph を作ります
- **家族関係を活用**: 温州みかん = 紀州みかん × 九年母 の F1 という**明確なpedigree**を、pangenome graph 上で直接検証できる教材設計
- **落とし穴も学ぶ**: 実データ解析で遭遇した実際のバグや誤解(PGGBの`-B`パラメータのサイレント失敗、haplotype leakage など)を教材に組み込んでいます
- **少人数構成でスケーラブル**: 3品種 × 2 hap = 6 haploid で、学生の実習環境でも動く規模

### 章構成

| 章 | タイトル | 学習内容 |
|---|---|---|
| [0](docs/ja/00_what_is_pangenome.md) | Pangenome とは | 概念、参照ゲノムとの違い、なぜ今 pangenome か |
| [1](docs/ja/01_general_workflow.md) | 構築の大まかな手順 | ツール共通の概念、PanSN、GFA、node/edge/path |
| [2](docs/ja/02_target_organisms.md) | 今回の対象 | 温州・紀州・九年母の関係、trio phasing の意義 |
| [3](docs/ja/03_data_acquisition.md) | データの取得 | BioProject、Plant GARDEN、ダウンロード実践 |
| [4](docs/ja/04_data_qc.md) | データの品質確認 | seqkit stats、期待値、異常検出 |
| [5](docs/ja/05_graph_construction.md) | Pangenome-graph 構築 | 手法比較(minigraph / MC / PGGB)、PGGB 実行 |
| [6](docs/ja/06_graph_qc.md) | グラフの QC | 4層の品質評価、圧縮率、transclose-batchの確認 |
| [7](docs/ja/07_graph_interpretation.md) | グラフの解釈 | VCF 読み、可視化、pedigree 確認 |
| [付](docs/ja/appendix_minigraph_cactus.md) | *付録*: Minigraph-Cactus | 同じデータを MC で扱う手順(未検証) |

### リポジトリ構成

```
citrus-pangenome-tutorial/
├── README.md
├── LICENSE
├── docs/
│   ├── ja/                      # 日本語版
│   └── en/                      # 英語版
├── scripts/                     # 実行スクリプト (all comments in English)
│   ├── setup_singularity_wrappers.sh
│   ├── download_plantgarden.sh
│   ├── qc01_stats.sh
│   ├── qc02_visualize.py
│   ├── pg01_prepare_input.sh
│   ├── pg02_run_pggb.sh
│   ├── pg03_qc_graph.sh
│   └── pg04_visualize.sh
├── tables/
│   ├── samplesheet.tsv
│   └── assembly_provenance_v1_v2.tsv
│
│   # --- created as you work through the tutorial (all git-ignored) ---
├── data/
│   ├── raw/                     # downloaded assemblies, left untouched
│   └── input/                   # PanSN-normalised, one FASTA per chromosome
├── results/                     # qc/ pggb/ graph_qc/ viz/
└── bin/                         # Singularity wrappers (optional)
```

### 前提となる計算環境

- Linux (Ubuntu 22.04+ or CentOS 7+)
- **HPC アクセス** (32+ CPU、128 GB RAM 推奨、chr03 は 200 GB あると安全)
- Singularity または Docker
- Python 3.10+
- ディスク容量 200 GB 以上

### 実行時間の目安

| ステップ | 対象 | 時間 |
|---|---|---|
| データ取得 | 6 haploid FASTA | 数分〜30分 (回線依存) |
| データ QC | 6 FASTA | 10-30分 |
| Pangenome構築 | 9 染色体 × PGGB | 5-10時間 (48 CPU、並列時) |
| Graph QC | 9 染色体 | 30分-1時間 |

### クイックスタート

```bash
# 1. Clone
git clone https://github.com/lsid-lab/citrus-pangenome-tutorial.git
cd citrus-pangenome-tutorial

# 2. 教材を辿る
# 日本語: docs/ja/00_what_is_pangenome.md
# 英語  : docs/en/00_what_is_pangenome.md
```

### 引用と謝辞

- 本教材で使用するアセンブリデータは **Isobe et al. (2023)** *bioRxiv* に基づきます
- PGGB 実行時の version-specific な `-B` の挙動は、共同研究者 **町田 宗聡** による綿密な診断実験で発見・記録されたものです
- 教材内容の一部は共同研究者との議論と、対話ベース開発によって構築されました

### ライセンス

- コード: MIT License (see `LICENSE`)
- 教材本文: CC BY 4.0

### 貢献・問い合わせ

Issue や Pull Request をお待ちしています。誤字・脱字の指摘、翻訳への貢献も歓迎します。

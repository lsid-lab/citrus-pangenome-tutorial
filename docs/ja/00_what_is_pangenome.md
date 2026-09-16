# 第0章 Pangenome とは

## 学ぶこと

- なぜ「一つの参照ゲノム」だけでは不十分なのか
- Pangenome の定義と3つの表現方法
- グラフベース pangenome が持つ利点
- なぜ今 pangenome が注目されるのか

---

## 0.1 「参照ゲノム」の限界

これまでの遺伝子解析は、**一つの代表個体の全ゲノム配列(reference genome)**を基準にしていました。ヒトなら GRCh38、シロイヌナズナなら TAIR10、柑橘なら Clementine v1.0 などです。

新しい個体を解析する際は:
1. その個体のシーケンスリードを参照ゲノムにマッピング
2. 参照との「違い」(SNP や小さな挿入欠失)を検出
3. その違いを個体の特徴として記述

これで多くの解析ができました。しかし**大きな問題**があります。

### 問題1: 参照にない配列は「見えない」

もし個体Aが「参照ゲノムには存在しない 5 kb の配列(=**大きな挿入**)」を持っていたとします。この 5 kb は、参照にマッピングできる場所がないため、多くの解析ツールは**見なかったこと**にしてしまいます。

もしその 5 kb にとても重要な遺伝子(病害抵抗性遺伝子、果実品質を決める遺伝子など)が入っていたら、それを見逃すことになります。

### 問題2: 参照は「一つの代表」に過ぎない

ヒトの参照 GRCh38 は主に一人の個体(実は複数人の混合)から作られています。しかし世界の 80億人の中で、GRCh38 と 100% 一致する人は一人もいません。

作物でも同じです。柑橘の参照 Clementine v1.0 は特定の個体のもので、温州みかんとは何百万箇所も違いがあります。それらすべての差を「変異」として個別に記録するのは、非効率的で誤解を招きやすい。

### 問題3: 集団の多様性を1個体で代表できない

参照ゲノムは「平均」でも「代表」でもなく、**たまたま最初にシーケンスされた個体**でしかありません。作物の育種や集団遺伝学では、**集団全体としての多様性**が知りたいのに、参照との差分の集まりでは扱いにくい。

---

## 0.2 Pangenome の定義

Pangenome は、これらの問題を解決するために生まれた概念です。

> **Pangenome** = **ある集団に存在する全ての遺伝物質**を表現するデータ構造

「集団」は種内(*intraspecific*)でも属内(*interspecific*、super-pangenome)でも構いません。柑橘の場合、「温州+紀州+九年母」のような栽培品種群が一つの集団です。

### Pangenome の3つの表現方法

同じ「全ての遺伝物質を表現する」でも、実装方法には歴史的に3つのアプローチがあります。

**① 遺伝子ベース(gene-based / linear pangenome)**

各個体から遺伝子リストを作り、遺伝子の**有無**(presence / absence)で集団の多様性を記述する方法。

- 全個体に共通する遺伝子: **core gene**
- 一部の個体だけが持つ遺伝子: **accessory gene**
- 特定の1個体だけが持つ遺伝子: **private gene**

細菌の pangenome (Tettelin et al., 2005) が最初期の代表例。植物では Golicz et al. 2016 (キャベツ) が先駆けです。

利点: シンプル、生物学的解釈が明快
限界: 遺伝子以外の情報(調節領域、非コード RNA、構造変異など)を扱えない

**② 配列ベース(sequence-based)**

複数個体のアセンブリを単純に結合し、リダンダントな部分を圧縮して、集団の全塩基配列を保持する方法。

利点: 全情報を保持できる
限界: サイズが大きい、構造や関係が明示的でない

**③ グラフベース(graph-based)** ← **本教材のアプローチ**

複数個体の配列を**グラフ**として表現する方法。同じ配列は一つの node に統合し、違いは分岐として表現します。

```
参照ゲノム(linear):
A T G G C C A T G C A T G

Graph pangenome:
     ┌─T─┐
A T G┤   ├C A T G
     └─A─┘  ↑
            この場所で T か A の変異あり
```

利点: **全ての違い(SNP、indel、SV、PAV)を統合的に表現**、参照非依存な解析が可能
限界: 実装が難しい、可視化・操作に慣れが必要

---

## 0.3 なぜ「今」 pangenome か

Pangenome の概念自体は 2005 年からありますが、実用的に構築できるようになったのはごく最近です。3つの技術的進展があります。

### 1. 長鎖シーケンサーの実用化

**PacBio HiFi**(2019 年〜)や **Oxford Nanopore ULK**(2020 年〜)により、15-25 kb の高精度リードが手軽に得られるようになりました。これで:
- 反復配列を正しくアセンブル可能に
- 大規模 SV も一つのリードで捉えられる
- 塩基精度 QV > 40(99.99%)を達成

### 2. Haplotype-resolved アセンブリ

**HiFi + Hi-C** や **trio phasing** で、二倍体の父方・母方染色体を**別々にアセンブル**することが可能になりました。これで:
- Heterozygous な variation を正しく捉える
- 集団解析の解像度が桁違いに上がる

**本教材の対象データ (Isobe et al., 2023) はまさにこの haplotype-resolved trio phasing の成果**です。

### 3. Pangenome 構築ツールの成熟

- **Minigraph** (Li et al., 2020): 大規模 SV 中心
- **Minigraph-Cactus** (Hickey et al., 2024): SV + 小変異、reference-guided
- **PGGB** (Garrison et al., 2024): reference-free、包括的 ← 本教材で使う

これらが実用的レベルに達し、標準的な HPC で数時間〜数日で pangenome が組めるようになりました。

---

## 0.4 Pangenome の応用例

代表的な pangenome プロジェクト:

- **Human Pangenome Reference Consortium** (Liao et al., 2023): 47 個体、94 haploid 相当。多様性を反映した参照
- **牛の super-pangenome** (Leonard et al., 2023): 系統横断の多様性
- **Sorghum pangenome** (Feng et al., 2025 など): 育種材料の遺伝資源発掘
- **Citrus pangenome** (Huang et al., 2023; Ollitrault et al., 2025): 属レベルの super-pangenome

これらは:
- **形質関連遺伝子の発見**: SNP-based GWAS で見えなかった SV-based の遺伝子座
- **育種材料の系統的評価**: 未利用遺伝資源の可視化
- **集団遺伝学**: 過去の交雑や選抜の履歴を復元

に貢献しています。本教材の柑橘 pangenome は小規模ですが、同じ技術体系を学ぶ良い入口になります。

---

## 0.5 教材で目指すもの

この教材を終える頃には、以下ができるようになります:

1. **概念の理解**: pangenome、PanSN、GFA、node/edge/path を自分の言葉で説明できる
2. **実行スキル**: HPC で PGGB を走らせ、graph を作れる
3. **QC 判断**: graph の品質を数値と可視化から評価できる
4. **解釈**: graph から VCF を抽出し、pedigree と照合して考察できる
5. **落とし穴の知識**: PGGB `-B` バグや haplotype leakage のような**実データ特有の癖**を経験している

---

## この章のまとめ

- 参照ゲノムだけでは集団の多様性を捉えられない
- Pangenome はその解決として、**全ての遺伝物質を統合的に表現**する
- グラフベース pangenome が現在の主流で、reference-free に多様性を扱える
- 長鎖シーケンサー、haplotype-resolved アセンブリ、pangenome ツールの成熟で、実用化された

次章では、pangenome を構築する**大まかな手順**(ツール共通の概念)を見ていきます。

## 参考文献

- Tettelin H, et al. (2005). Genome analysis of multiple pathogenic isolates of *Streptococcus agalactiae*: implications for the microbial "pan-genome". *PNAS* 102:13950-13955.
- Golicz AA, et al. (2016). The pangenome of an agronomically important crop plant *Brassica oleracea*. *Nat Commun* 7:13390.
- Liao WW, et al. (2023). A draft human pangenome reference. *Nature* 617:312-324.
- Garrison E, et al. (2024). Building pangenome graphs. *bioRxiv*.
- Isobe S, et al. (2023). Haploid-resolved and chromosome-scale genome assembly in *Citrus unshiu* and its parental species. *bioRxiv* 2023.06.02.543356.

---

[← README](../../README.md) | [第1章: 構築の大まかな手順 →](01_general_workflow.md)

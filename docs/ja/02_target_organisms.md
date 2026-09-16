# 第2章 今回の対象:温州みかん・紀州みかん・九年母

## 学ぶこと

- 3品種の生物学的背景と歴史
- 温州みかん = 紀州みかん × 九年母 の親子関係
- Trio phasing の意義と、この教材との関連
- なぜこの 3品種が pangenome の教材として理想的か

---

## 2.1 温州みかん(*Citrus unshiu*)

日本の冬を代表する柑橘、**温州みかん**(以下 STS: **S**a**ts**uma)は、江戸時代初期に**鹿児島県長島町**で見つかった品種です。

- **染色体数**: 2n = 18 (haploid n = 9)
- **ゲノムサイズ**: 約 320-360 Mb (haploid)
- **繁殖様式**: 種なし(単為結果性)、無融合種子(nucellar polyembryony、母親のクローンを生む種子)
- **代表的な系統**: 宮川早生、興津早生、青島温州、南柑20号 など多数

種なしという特徴は栽培には有利ですが、**遺伝学的解析には困難**を伴います。通常の交雑実験(親A × 親B → F1 の観察)ができないためです。

温州みかんの起源は長らく謎でしたが、**近年のゲノム解析** (Fujii et al., 2016; Wu et al., 2018)により、
**紀州みかん(KSH)× 九年母(KNN)の F1 雑種**であることが遺伝的に確認されました。

### 遺伝的組成

Wu et al. 2018 の解析では、温州みかんは主に:
- **Mandarin (*C. reticulata*)** の遺伝子プールが 80% 以上
- 少量の **Pummelo (*C. maxima*)** の遺伝子プールを持つ

これは両親(紀州・九年母)由来の admixture を反映しています。

---

## 2.2 紀州みかん(*Citrus kinokuni*)

**紀州みかん**(KSH: **K**i**sh**u)は、**温州みかんの母親**です。

- **由来**: 中国原産、日本には平安時代(9世紀頃)には既に存在
- **名前**: 現在の和歌山県(紀州)で盛んに栽培されたことから
- **特徴**: 種あり、果実は小さい、皮が薄い
- **歴史的意義**: 江戸時代、日本の柑橘といえば紀州みかんでした
- **染色体数**: 2n = 18

紀州みかんは近年、栽培量が激減し、絶滅危惧品種扱いされることもあります。しかし温州の親としての遺伝学的重要性から、保護・研究の対象となっています。

### 遺伝的特徴

- **Mandarin ancestry がほぼ 100%**(pummelo 混入は微量)
- そのため「**純粋な mandarin**」に近いとされる
- 温州の "mandarin 側の親" として貢献

---

## 2.3 九年母(*Citrus nobilis* var. *kunip*)

**九年母**(KNN: **Kun**enbo)は、**温州みかんの父親**にあたります。

- **由来**: インドシナ半島原産、日本には室町時代(1350年頃)に伝来
- **名前の由来**: 諸説あるが、「植えて9年で花が咲く」など
- **特徴**: 種あり、比較的大きな果実、独特の香り
- **染色体数**: 2n = 18

現在の日本ではほぼ栽培されていませんが、鹿児島県などに保存系統があります。

### 遺伝的特徴

- **Mandarin ベースだが、pummelo (*C. maxima*) の遺伝子プールが 10-15%**混入
- この pummelo ancestry が、温州みかんに部分的に受け継がれています

---

## 2.4 温州みかんの家族図

3品種の関係を図で示すと:

```
        紀州みかん (KSH)          九年母 (KNN)
   [Mandarin ~100%]         [Mandarin ~85%]
   [ hap1 | hap2 ]          [ hap1 | hap2 ]
         [Pummelo 15%]
             \                     /
              \                   /
               \                 /
                \               /
              (交雑、F1形成)
                   ┌───────┐
                   ↓       ↓
              温州みかん(STS)
              [Mandarin ~90%, Pummelo ~10%]
        [hap1=Kishu由来 | hap2=Kunenbo由来]
```

**重要**: 温州みかんの2つのハプロタイプは、それぞれ**片親から1本ずつ受け継がれた**ものです:
- STS_hap1 ≈ 紀州みかんのハプロタイプの片方(または組み換え体)
- STS_hap2 ≈ 九年母のハプロタイプの片方(または組み換え体)

---

## 2.5 Trio phasing の意義

通常の2倍体ゲノムアセンブリでは、**父由来と母由来の染色体を区別できません**。両方が混ざった「collapsed diploid」になってしまう。

しかし**両親のゲノムも一緒にシーケンス**すれば、子のリードを親由来で分類できます:

```
親A (紀州) のk-mer:  ACGT, GTCA, ...   (紀州特有)
親B (九年母) のk-mer: CCAA, TTAG, ...  (九年母特有)

子(温州)のリード R1: このk-merは親A由来 → hap1 (Kishu由来)
子(温州)のリード R2: このk-merは親B由来 → hap2 (Kunenbo由来)
```

これが **trio phasing** です (Cheng et al., 2021)。

**Isobe et al. (2023) は、この trio phasing を実際に温州みかんで実施した論文**です。以下の 6 つの haplotype-resolved アセンブリを公開しました:

| 品種 | アセンブリ ID | 内容 |
|---|---|---|
| 温州 STS | CUNphKi_r1.0 | 温州の hap1 (**Ki**shu由来) |
| 温州 STS | CUNphKu_r1.0 | 温州の hap2 (**Ku**nenbo由来) |
| 紀州 KSH | CKIhap1_r1.0 | 紀州の hap1 |
| 紀州 KSH | CKIhap2_r1.0 | 紀州の hap2 |
| 九年母 KNN | CKUhap1_r1.0 | 九年母の hap1 |
| 九年母 KNN | CKUhap2_r1.0 | 九年母の hap2 |

**命名の読み解き**:
- **C** = *Citrus*
- 種の略号: **UN** = **un**shiu, **KI** = **ki**nokuni, **KU** = **ku**nip
- **ph** = phased(親由来を特定した)
- 末尾: **Ki** / **Ku** で親由来を示す(温州のみ)、または **hap1/hap2** で任意ラベル(親系統)

---

## 2.6 なぜこの 3品種が pangenome 教材として理想的か

Pangenome の教材として、この 3品種の組み合わせは**極めて理想的**です。理由:

### 1. Pedigree が完全に既知

温州 = 紀州 × 九年母 は **遺伝的に確認済み** (Fujii et al., 2016)。
Pangenome graph を作った後で、この関係が graph 上に見えるか**検証可能**です。

### 2. 適度な遺伝的距離

- 3品種はすべて *Citrus* 属で相互に近縁
- しかし紀州(pure mandarin)、九年母(mandarin + pummelo)は明確に区別できる
- 温州の 2 ハプロタイプで、両方の親を "見る" ことができる
- → **多様性の解析としても、系統関係の解析としても機能する**

### 3. 計算資源が現実的

3品種 × 2 hap = **6 haploid**。学部生の実習環境でも動かせる規模。
1 染色体あたり 32 CPU で 30-90 分程度。全 9 染色体でも並列すれば数時間。

### 4. データが完全公開

Isobe et al. (2023) が **Plant GARDEN** で全アセンブリを公開しています。追加費用ゼロで、誰でも同じデータで実習できる。

### 5. 教科書的な"最小の家族"

「親2つ + F1子1つ」は、**遺伝学の最も基本的な単位**です。Mendel の法則の直接適用対象。より複雑な集団解析(GWAS、QTL mapping)へのゲートウェイになります。

---

## 2.7 教材で追う「発見のストーリー」

この 3品種の graph を通じて、以下を実際に**発見・確認**していきます:

**発見1**: 温州の hap1 は紀州の hap のどれか一方と特に似ている
- STS_hap1 vs KSH_hap1(または hap2)の類似度が高い

**発見2**: 温州の hap2 は九年母の hap と似ている
- STS_hap2 vs KNN_hap*

**発見3**: 温州の 2 ハプロタイプ間で SV が存在する
- vg deconstruct で数千の SV が検出される

**発見4**: 一部の SV は親から遺伝している
- KSH や KNN に同じ SV が見つかる

Graph の QC や解析の中で、想定外の観察や実データ特有の癖に出会うこともあるでしょう。それらもすべて、教材を通じて体験することになります。

---

## 2.8 参考: なぜ「これ以外の柑橘品種」は使わないか

18 品種の PacBio HiFi データが公開されていますが(Kiryu et al., 2026)、**アセンブリが公開されているのは STS/KSH/KNN のみ**です。他の 15 品種(甘夏、伊予柑、ポンカン、レモン、シトロンなど)は、HiFi リードのみで、アセンブリはまだ公開されていません(2026年時点)。

もし他の品種のアセンブリを使いたい場合は、生リードから自前でアセンブル(hifiasm 実行、数時間〜1日/品種)する必要があります。

これは教材のスコープを超えるため、本教材では **公開されている 3品種に限定**します。

---

## この章のまとめ

- 対象: 温州みかん (STS) = 紀州みかん (KSH) × 九年母 (KNN) の F1
- Isobe et al. (2023) が haplotype-resolved trio-phased アセンブリを公開
- 6 haploid のデータで、家族関係を pangenome graph 上で直接検証可能
- 教材規模が現実的で、"pangenome の学習に理想的な最小の家族"

次章では、これらのデータを実際に**ダウンロード**する方法を学びます。

## 参考文献

- Fujii H, et al. (2016). Parental diagnosis of satsuma mandarin (*Citrus unshiu* Marc.) revealed by nuclear and cytoplasmic markers. *Breed Sci* 66:683-691.
- Wu GA, et al. (2018). Genomics of the origin and evolution of *Citrus*. *Nature* 554:311-316.
- Isobe S, et al. (2023). Haploid-resolved and chromosome-scale genome assembly in *Citrus unshiu* and its parental species. *bioRxiv* 2023.06.02.543356.
- Cheng H, et al. (2021). Haplotype-resolved de novo assembly using phased assembly graphs with hifiasm. *Nat Methods* 18:170-175.
- Kiryu Y, et al. (2026). AlleleMiner: a long-read pipeline for gene-wise de novo allele phasing and variant detection in diploid citrus cultivars. *DNA Res* 33:dsag004.

---

[← 第1章](01_general_workflow.md) | [第3章: データの取得 →](03_data_acquisition.md)

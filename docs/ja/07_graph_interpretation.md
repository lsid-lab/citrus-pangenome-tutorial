# 第7章 グラフの解釈

## 学ぶこと

- VCF ファイルからの変異情報の読み方
- odgi による pangenome graph の可視化
- Pedigree 検証の実例
- 結果の生物学的解釈と、教材の限界

---

## 7.1 この章のゴール

前章で技術的品質(Level 1-2)と生物学的妥当性(Level 3-4)を確認しました。ここからは、**その結果を実際に生物学的に読み解く**段階に入ります。

具体的に見るもの:

1. **VCF の読み方**:どんな変異が graph に含まれているか
2. **可視化**:graph を目で見て構造を理解
3. **Pedigree の確認**:「F1 = 親A × 親B」が graph 上でどう見えるか
4. **教材の限界**:今回の graph で何ができて、何ができないか

---

## 7.2 VCF ファイルの構造と読み方

### 7.2.1 VCF の基本

VCF (Variant Call Format) は、参照配列と比較した**変異の位置と種類**を記述する標準形式です。

例(chr09 の一部):
```
#CHROM   POS      ID    REF   ALT      QUAL   FILTER   INFO   FORMAT   CUN#1  CUN#2  CKI#1  ...
CUN#1#chr09   1234  .     A     G        .      .        AT=snp   GT       0          1          0        ...
CUN#1#chr09   2000  .     ACGT  A        .      .        AT=del3  GT       0          1          1        ...
CUN#1#chr09   3500  .     C     CGATCG   .      .        AT=ins5  GT       0          0          1        ...
```

各列:
- **CHROM**: 座標基準 path (今回は `CUN#1#chr09`)
- **POS**: 位置(1-based)
- **REF**: 参照(座標基準 path)の塩基
- **ALT**: 変異体
- **GT** (Genotype): 各 sample での genotype
  - `0` = REF、`1` = ALT1、`2` = ALT2、`.` = 未定義

### 7.2.2 SV、SNP、indel の分類

変異のタイプは REF と ALT の長さから判定:

```bash
# SNP: REF と ALT が両方 1 bp
awk 'length($4)==1 && length($5)==1' chr09.vcf

# Insertion: REF が短く、ALT が長い
awk 'length($4)<length($5) && length($5)>1' chr09.vcf

# Deletion: REF が長く、ALT が短い
awk 'length($4)>length($5) && length($4)>1' chr09.vcf

# SV: どちらかが 50 bp 以上
awk 'length($4)>=50 || length($5)>=50' chr09.vcf
```

### 7.2.3 サンプル別の変異数を数える

各品種にどれだけの ALT variant があるか:

```bash
# CUN#2 (座標基準の反対 hap) の変異数
bcftools view -s CUN#2 chr09.vcf | \
  bcftools query -f '[%GT]\n' | grep -c '^1'
```

### 7.2.4 期待される結果

温州みかんの haplotype 間で見つかる変異(chr09、~30 Mb):

|Variant タイプ|期待数|
|---|---|
|SNPs|10,000〜30,000|
|Small indels (<50 bp)|5,000〜15,000|
|SVs (≥50 bp)|100〜500|

---

## 7.3 Pangenome graph の可視化

### 7.3.1 1D linear view (odgi viz)

Graph の各 path を横軸に配置した「pileup」的な可視化:

```bash
odgi viz -i chr09.smooth.final.og -o chr09_viz.png -x 1500 -y 500
```

見え方の解釈:
- **横軸**: graph 上の位置(座標基準 path)
- **縦軸**: 各 path(6 haploid)
- **色**: path が graph の同一領域を共有していれば同じ色
- **色分割**: 分岐(=変異)がある場所

期待される正常な graph の特徴:
- ✅ 6 本の path が横一線に走っている
- ✅ 大部分の領域で色が揃っている(共有配列)
- ✅ 縦に色分割している点在領域(変異)

異常のパターン:
- ❌ 短い縦線が数百本以上 → 断片化しすぎ
- ❌ 途中で消える path → 保存性の失敗

### 7.3.2 2D layout view (odgi layout + draw)

Graph の全体構造を空間的に配置:

```bash
odgi layout -i chr09.smooth.final.og -o chr09.lay -t 8
odgi draw -i chr09.smooth.final.og -c chr09.lay -p chr09_2d.png -H 1000 -C
```

見え方:
- **backbone**: 全 path 共有の太い連続線
- **bubble**: 一部の path のみの分岐(小規模 SV)
- **long branch**: 大規模な private sequence

**正常な pangenome**の 2D 図は、中央に太い backbone + 周辺に多数の小さな bubble、が典型的。

**異常**な場合は、6 本の path がバラバラに配置されて絡み合い、backbone が見えない。

### 7.3.3 実行スクリプト

`scripts/pg04_visualize.sh` が両方を実行します:

```bash
bash scripts/pg04_visualize.sh . chr09
```

---

## 7.4 Pedigree 検証の実例

### 7.4.1 温州みかんの trio phasing を graph で確認

第6章 Level 3 の pedigree 検定の詳細を、生物学的に解釈します。

`similarity` 出力から、以下を計算:

|比較|Jaccard(期待)|
|---|---|
|CUN_hap1 vs CKI_hap1 / CKI_hap2|**高い**(どちらに対しても。CUN_hap1 は両者のモザイクなので、片方だけが突出することはない)|
|CUN_hap1 vs CKU_hap*|**相対的に低い**(父方は CUN_hap1 に寄与していない)|

同様に、CUN_hap2 は CKU の 2 本に対して高く、CKI に対して低くなります。

> 絶対値は haplotype leakage (§4.6) の影響で下振れするので、**閾値ではなく品種間の大小関係**で判断してください。また、染色体全体の Jaccard 1 個では「親のどちらのハプロタイプ由来か」までは決まりません(§2.4)。それを見るには §7.5.2 の窓ごとの解析が必要です。

### 7.4.2 具体的な pedigree pattern

「hap1 が父母どちらから来たか」を graph 上で追跡:

```bash
# CUN_hap1 との各 CKI hap の Jaccard
awk '$1 ~ /^CUN#1/ && $2 ~ /^CKI/' chr09_similarity.tsv
```

もし出力が:
```
CUN#1#chr09  CKI#1#chr09  ...  Jaccard = 0.95
CUN#1#chr09  CKI#2#chr09  ...  Jaccard = 0.55
```

のように **CKI 側が CKU 側より一貫して高ければ、CUN_hap1 は母方 (CKI) 由来**と結論できます。

ここで **「CKI_hap1 のほうが高いから CUN_hap1 = CKI_hap1 由来」と結論してはいけません。** CUN_hap1 は CKI_hap1 と CKI_hap2 のモザイクであり(§2.4)、この差は「モザイクの中で CKI_hap1 区間の占める割合がやや大きい」ことを意味するにすぎません。どの区間がどちらに由来するかは、次節の窓ごとの解析で初めて見えます。

### 7.4.3 Haplotype leakage の影響を認める

第4章で発見した haplotype leakage(CUN の 2 hap が親系統より大きい)は、以下のように現れます:

- CUN_hap1 vs CKI の Jaccard がやや低め(例えば 0.7 期待が 0.6 に見える)
- 平均値だけ見ると pedigree consistency が破綻して見える

しかし、**MAX-based で見れば pedigree は成立**しており、graph 自体は正しく作られている。これが第 6.4.4 節での重要な学習ポイントでした。

---

## 7.5 生物学的な発見の例

Graph から見えるかもしれない生物学的パターン:

### 7.5.1 温州のプミロ由来領域

Wu et al. 2018 によれば、温州はメインは mandarin だが少量の pummelo 遺伝子プールを持ちます。もし graph 上で:

- CUN の一部の path 領域が CKI や CKU と一致せず、独自の path になっている

なら、それは pummelo 由来の可能性があります。プミロ由来領域は特定の chr の特定領域に集中している(Fujii et al. 2016 で報告)ため、可視化で見つけやすいはず。

### 7.5.2 F1 の組み換え領域

§2.4 で見たとおり、CUN_hap1 は CKI_hap1 / CKI_hap2 の**組み換えモザイク**です。染色体に沿って窓ごとに類似度を計算すると:

- CUN_hap1 が CKI_hap1 に近い領域
- CUN_hap1 が CKI_hap2 に近い領域

が交互に現れ、その境目が**組み換え点 (crossover breakpoint)** です。母親の減数分裂で実際に起きた乗換えの位置を、graph 上から読み取っていることになります。

### 7.5.3 反復配列や TE の存在

大規模 SV (5-15 kb) の多くは、**Transposable Element (TE) の挿入**による可能性があります。Kiryu et al. 2026 でも、5.1 kb の LTR retrotransposon 挿入が確認されています。

VCF から SV を抽出し、TE データベースに BLAST すれば、この解析が可能です(教材のスコープ外)。

---

## 7.6 教材で得られる graph の限界

**この教材の graph でできること**:

- ✅ 3 品種間の SV・SNP の総合的なカタログ
- ✅ Trio phasing の validation
- ✅ Pangenome の基本操作(構築、QC、解釈)の習得
- ✅ PGGB / odgi / vg などの主要ツールの実践

**できないこと**:

- ❌ 集団遺伝学(3 品種では集団の議論は不可)
- ❌ GWAS(表現型データがない)
- ❌ 属レベルの super-pangenome(他の *Citrus* 種が含まれていない)
- ❌ 高精度な系統関係の推定(rooted phylogeny、divergence time)

**Haplotype leakage の影響**:

- CUN の 2 hap のサイズや Jaccard 値には、trio phasing の不完全性が影響
- 生物学的解釈の際は、この癖を意識する
- 実データは常に完璧ではない、というリアリティを学べる

---

## 7.7 次のステップ

この教材で得た知識と経験を次に活かすには:

**A. より大規模な pangenome プロジェクト**
- 15+ の柑橘品種を集めた super-pangenome(HiFi リードから自前アセンブリ)
- 他属の Rutaceae(Poncirus、Citropsis など)を加えた広域比較

**B. 集団遺伝学への展開**
- Kiryu et al. 2026 の 18 品種で allele-level のフレームワーク
- Kunenbo × Kishu の F1 集団(RAD-Seq 99 個体、PRJDB15866)を使った QTL マッピング

**C. トランスクリプトーム統合**
- Isobe 2023 の 8 組織 Iso-Seq を graph に投射
- Allele-specific expression の解析

**D. 育種応用**
- 有望な SV の見つかった座位で、marker-assisted selection のマーカー設計
- 新品種評価の genotyping パイプラインへの組み込み

---

## 7.8 この章のまとめ

- VCF から SNP、indel、SV を分類抽出できる
- odgi による 1D/2D 可視化で graph の全体構造を確認
- Pedigree 検証は MAX-based で "biological validity" を確認
- Haplotype leakage は解釈上の癖として認識、graph 自体は使える
- 教材の graph は基本習得用、応用は次のプロジェクトへ

---

## 参考文献

- Wu GA, et al. (2018). Genomics of the origin and evolution of *Citrus*. *Nature* 554:311-316.
- Fujii H, et al. (2016). Parental diagnosis of satsuma mandarin. *Breed Sci* 66:683-691.
- Kiryu Y, et al. (2026). AlleleMiner. *DNA Res* 33:dsag004.

---

[← 第6章](06_graph_qc.md) | [目次に戻る](../../README.md)

---

## 教材完了 🎉

第 0 章から第 7 章まで、pangenome の基礎から実践、解釈までを一貫して学びました。

- **基礎**: pangenome とは、workflow、対象生物
- **実践**: データ取得、QC、graph 構築
- **解釈**: 品質評価、生物学的解釈

これで、あなたは他の生物種でも同じアプローチで pangenome を構築できます。実際の研究や育種プロジェクトでの活用を期待しています。

# 第7章 グラフの解釈

## 学ぶこと

- VCF ファイルからの変異情報の読み方
- odgi による pangenome graph の可視化
- Pedigree 検証の実例
- 結果の生物学的解釈と、教材の限界
- **アセンブリが「バージョンを持つ成果物」であること、そして引っかかった観察をどう扱うか**

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

> 絶対値は CUN の path が長いぶん下振れするので(§4.4、§6.4.4)、**閾値ではなく品種間の大小関係**で判断してください。また、染色体全体の Jaccard 1 個では「親のどちらのハプロタイプ由来か」までは決まりません(§2.4)。それを見るには §7.5.2 の窓ごとの解析が必要です。

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

### 7.4.3 CUN の path が長いことの影響

第4章で見た CUN の 2 hap が親系統より大きいという観察(§4.4)は、以下のように現れます:

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

> **ここで手を止めて数えてみてください。** §2.4 で見たとおり、乗換えは 1 本の染色体あたり**最低 1 回、通常は 1〜3 回**です。窓ごとの切り替わりが、それよりも明らかに多く見えていませんか? もしそうなら、その切り替わりの全部が本物の乗換えだとは考えにくいことになります。この違和感は §7.6 で回収します。

### 7.5.3 反復配列や TE の存在

大規模 SV (5-15 kb) の多くは、**Transposable Element (TE) の挿入**による可能性があります。Kiryu et al. 2026 でも、5.1 kb の LTR retrotransposon 挿入が確認されています。

VCF から SV を抽出し、TE データベースに BLAST すれば、この解析が可能です(教材のスコープ外)。

---

## 7.6 実は —— このデータには v2 がある

ここまでの解析で、引っかかった点がいくつかあったはずです。

1. **CUN の 2 hap が、親 2 品種より 25-55 Mb 大きい**(§4.4)
2. **CUN の Jaccard 類似度が、path が長いぶん下振れする**(§6.4.4、§7.4.3)
3. **窓ごとの解析で、組み換え点が生物学的に想定される回数より多く見える**(§7.5.2)

どれも単独なら「そういうものか」で流せる程度の違和感です。実際、本教材もここまで「観察として記録しておく」以上のことはしてきませんでした。

### v2 の公開

本教材が使っている Isobe et al. (2023) のアセンブリは **version 1** です。その後、同じ 3 品種の **version 2** が、農研機構の **MiGD2** で公開されています。

> **MiGD2 (Mikan Genome Database 2)**: <https://mikan.dna.naro.go.jp/migd2/>

**v2 で同じ解析を回すと、上の 3 つの違和感はいずれも解消します。**

### 何が変わったか

**(1) アセンブリサイズ**

v2 では、**CUN の 2 hap のサイズが紀州・九年母とほぼ同じ範囲に収まります。** §4.4 で記録した 25-55 Mb の差は、v2 では消えます。

<!-- TODO: v2 の 6 haploid の seqkit stats (num_seqs / sum_len_Mb / N50_Mb / GC%) を表で入れる -->

§4.4 では「温州として公表されているサイズ(346-360 Mb)の範囲内だから、それだけでは異常とは言えない」と判断しました。その判断は、**手元にあった情報に対しては妥当**でした。決着をつけたのは、こちらの推論ではなく、**同じデータ提供者が同じ個体を作り直したもの**です。同一パイプラインで作り直した v2 と比べて初めて、v1 の CUN に余分な配列が含まれていたと言えるようになりました。

**(2) Phasing**

より重要なのはこちらです。v1 では、**染色体の区間ごとに親由来の割り当てを取り違えている箇所**(switch error)がありました。v2 ではこれが修正されています。

これは §7.5.2 に直結します。窓ごとに「CUN_hap1 は CKI_hap1 と CKI_hap2 のどちらに近いか」をプロットしたとき、切り替わりは次のどちらでもあり得ます:

- **本物の組み換え点** —— 母親の減数分裂で実際に起きた乗換え
- **phasing の switch error** —— 解析上のアーティファクト

そして、**窓ごとの類似度だけを見ても、この 2 つは区別できません。** どちらも「ある位置から先で、似ている相手が入れ替わる」というまったく同じ形で現れるからです。§7.5.2 で「乗換えの位置を graph 上から読み取っている」と書きましたが、正確には**「乗換えか、phasing の誤りか、どちらかの位置」を読み取っていた**わけです。

v1 と v2 で切り替わり点を突き合わせれば、**v2 で消えるものが switch error、残るものが本物の組み換え点**、という切り分けができます。生物学的に妥当な回数(1 染色体あたり 1〜3 回)に落ち着くはずです。

### ここから学べること

**1. アセンブリは「バージョンを持つ成果物」である**

ゲノムアセンブリは確定した事実ではなく、**その時点のデータとアルゴリズムでの最良の推定**です。論文が出た後でも、リードが増えたり、アセンブラや phasing ツールが改良されれば更新されます。**解析を始める前に、使おうとしているアセンブリが最新版かを必ず確認してください。** 公開ポータル(Plant GARDEN、MiGD2、NCBI など)と、著者グループの後続論文の両方を見るのが確実です。

**2. 引っかかった点は、説明をつけて消さずに記録しておく**

§4.4 の CUN のサイズ差に、本教材はあえて結論を与えませんでした。手元のデータでは判定できなかったからです。**判定できないことを判定しない**のは後ろ向きに見えますが、この場合は正解でした。違和感を無理に説明しきらずにログしておけば、新しい情報が来たときに**照合できます**。もっともらしい説明で蓋をしてしまうと、照合する機会そのものを失います。

**3. 系統的な誤りは、単独の指標では見えにくい**

サイズだけを見れば「公表値の範囲内」でした。Jaccard だけを見れば「path が長いから」で説明がつきました。組み換え点だけを見れば「そういうものか」で流せました。**3 つが同じ方向を指していたことにこそ意味があった**わけです。QC 指標は単独ではなく、束で見てください。

### v2 で再実行するには

本教材のパイプラインはアセンブリのパスを `tables/samplesheet.tsv` から読むだけなので、v2 での再実行は簡単です:

1. MiGD2 から v2 の 6 haploid を取得する
2. `tables/samplesheet.tsv` の `hap1_path` / `hap2_path` を v2 のファイルに差し替える
3. 第4章の `qc01_stats.sh` から、この章の解析まで同じ手順で流し直す

**スクリプトは 1 行も変更せずに通ります。** v1 と v2 の結果を並べて比べるのが、この教材のいちばん良い仕上げになるかもしれません。

---

## 7.7 教材で得られる graph の限界

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

**v1 のアセンブリに由来する制約**(§7.6):

- CUN の 2 hap のサイズと Jaccard 値には、v1 の余分な配列が影響している
- §7.5.2 で見える組み換え点には、v1 の phasing switch error が混ざっている
- したがって本教材の graph は、**手法を学ぶには十分だが、生物学的な結論を出す土台には向かない**
- 結論を出したいなら、v2 で作り直してください(§7.6)

---

## 7.8 次のステップ

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

## 7.9 この章のまとめ

- VCF から SNP、indel、SV を分類抽出できる
- odgi による 1D/2D 可視化で graph の全体構造を確認
- Pedigree 検証は MAX-based で "biological validity" を確認
- **使っていたアセンブリには v2 があり、サイズ差と phasing の switch error はそこで修正されている**(§7.6)
- 窓ごとの類似度だけでは、**本物の組み換え点と phasing の誤りを区別できない**
- 引っかかった観察は、説明をつけて消さずに記録しておく —— 後から照合できる
- 教材の graph は基本習得用、応用は次のプロジェクトへ

---

## 参考文献

- Wu GA, et al. (2018). Genomics of the origin and evolution of *Citrus*. *Nature* 554:311-316.
- Fujii H, et al. (2016). Parental diagnosis of satsuma mandarin. *Breed Sci* 66:683-691.
- Kiryu Y, et al. (2026). AlleleMiner. *DNA Res* 33:dsag004.
- **MiGD2 (Mikan Genome Database 2、農研機構)**: <https://mikan.dna.naro.go.jp/migd2/> (本教材が使う v1 アセンブリの後継版 v2 の公開先)

---

[← 第6章](06_graph_qc.md) | [目次に戻る](../../README.md)

---

## 教材完了 🎉

第 0 章から第 7 章まで、pangenome の基礎から実践、解釈までを一貫して学びました。

- **基礎**: pangenome とは、workflow、対象生物
- **実践**: データ取得、QC、graph 構築
- **解釈**: 品質評価、生物学的解釈

これで、あなたは他の生物種でも同じアプローチで pangenome を構築できます。実際の研究や育種プロジェクトでの活用を期待しています。

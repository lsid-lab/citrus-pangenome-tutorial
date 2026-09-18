# 第6章 グラフの QC

## 学ぶこと

- Pangenome graph の 4 層 QC アプローチ
- 各層で何を、なぜ、どう見るか
- 品質判定の具体的な閾値
- 失敗パターンと対処法

---

## 6.1 なぜ 4 層で見るか

Pangenome graph が「使えるか」を判断するには、複数の観点から評価する必要があります:

|レベル|観点|見る指標|合否の意味|
|---|---|---|---|
|Level 1|**構造**|node/edge/path 数、圧縮率|graph の技術的健全性|
|Level 2|**保存性**|入力 vs graph の path 長|入力データの損失なし|
|Level 3|**生物学的妥当性**|path 間の類似度と pedigree|生物学的に意味が通るか|
|Level 4|**SVシグナル**|VCF の SV 数と分布|下流解析に使えるか|

**上位が失敗している時に下位を見る意味はない**という原則があります。Level 1 で不合格なら、SV を数える前に構造を直しましょう。

---

## 6.2 Level 1: 構造チェック

### 6.2.1 何を見るか

`odgi stats -S` の出力を確認します:

```bash
odgi stats -i chr09.fa.gz.<hash>.smooth.final.og -S
```

出力(タブ区切り、5列):
```
#length     nodes    edges    paths    steps
49786530    1381521  1874185  6        6020593
```

各列の意味:
- **length**: graph 内のすべての node の塩基数合計 = **graph_bp**
- **nodes**: node 数
- **edges**: edge 数
- **paths**: path 数(=入力 haploid 数、期待: 6)
- **steps**: すべての path が graph を通る**総ステップ数**

### 6.2.2 主要指標:圧縮率

Level 1 で最も重要なのは**圧縮率**です:

```
圧縮率 = graph_bp / input_bp
```

|圧縮率|判定|意味|
|---|---|---|
|0.19-0.30|✅ 良好|近縁 haploid で正常な pangenome|
|0.30-0.40|⚠️ 許容|やや divergent、要確認|
|0.40-0.60|⚠️ 疑問|パラメータ見直し|
|**> 0.90**|❌ **異常**|**`-B` バグ**の可能性大|

### 6.2.3 `.params.yml` の transclose-batch チェック

第5章 §5.5.3 で触れた、特定バージョンの PGGB での `-B` 挙動を確認します:

```bash
grep transclose-batch 03_pangenome/by_chr/chr09_pggb/*.params.yml
```

期待:
- `transclose-batch: 1000000000` → ✅ `-B 1G` が正しく渡された
- `transclose-batch: 1000000` → ⚠️ **既定値のパースに失敗している可能性**、再実行推奨

`.err` ファイルが空でも、exit code が 0 でも、この 1 行で挙動を確認できます。将来のバージョンでは修正されている可能性がありますが、確認習慣として身につけておくと安全です。

### 6.2.4 その他の Level 1 指標

**path 数 = 6 が必須**:
- 3 品種 × 2 hap = 6
- 6 でない場合、入力の PanSN 命名や pg01 の実行を確認

**inflation 比**(補助的):
```
inflation = graph_bp × n_haplotypes / input_bp = 圧縮率 × 6
```

- 期待: 1.0-1.6
- 圧縮率が 0.25 なら inflation ≈ 1.5

**参考: 共同研究者の実測値**(全9染色体):

|chr|圧縮率|node 数|
|---|---|---|
|ch1|0.2391|1,074,276|
|ch2|0.2484|1,867,874|
|ch3|0.2604|2,882,578|
|ch4|**0.1862**|706,073|
|ch5|0.2638|1,667,106|
|ch6|0.2443|1,141,586|
|ch7|0.2743|1,901,707|
|ch8|**0.2778**|1,483,201|
|ch9|0.2541|1,383,141|

全 9 染色体で **0.186 〜 0.278** に収まっており、共通の閾値 ≤ 0.35 で全部 PASS しています。

---

## 6.3 Level 2: パス保存性

### 6.3.1 何を確認するか

「入力した 6 haploid が、graph 内に**完全な path として保存**されているか」を確認します。もし path のどこかが失われていれば、graph は元の入力を再現できません。

### 6.3.2 手法

```bash
# graph から path を FASTA として抽出
odgi paths -i chr09.smooth.final.og -f > chr09_paths.fa

# パス長を取得
seqkit fx2tab -nl chr09_paths.fa
```

出力例:
```
CUN#1#chr09    45,179,797
CUN#2#chr09    29,817,104
CKI#1#chr09      28,530,969
CKI#2#chr09      31,864,882
CKU#1#chr09    30,157,514
CKU#2#chr09    30,337,761
```

これを**入力 FASTA の各配列長**と比較します。

### 6.3.3 判定基準

|差分|判定|
|---|---|
|0% (完全一致)|✅ 理想|
|< 0.1%|✅ 許容|
|0.1-1%|⚠️ 要調査|
|**> 1%**|❌ **path 失われ、致命的**|

### 6.3.4 実行例

`pg03_qc_graph.sh` がこの比較を自動化します:

```bash
bash scripts/pg03_qc_graph.sh . chr09
```

期待される出力:
```
保存性判定: YES_max_diff=0.000%_at_
```

max_diff が 0.000% なら**入力配列は完全に graph に保存**されており、pangenome として正常です。

---

## 6.4 Level 3: 系統整合性(pedigree チェック)

### 6.4.1 温州みかんの pedigree で graph を "検証"

3 品種の系統関係が既知(CUN = CKI × CKU)なので、graph 上の path 間類似度が pedigree と整合するかを確認できます。

期待:
- CUN_hap1 (Kishu 由来) と CKI の hap のどちらか → 高い類似度
- CUN_hap2 (Kunenbo 由来) と CKU の hap のどちらか → 高い類似度

### 6.4.2 odgi similarity の実行

```bash
odgi similarity -i chr09.smooth.final.og > chr09_similarity.tsv
```

出力の列:
```
group.a  group.b  a.length  b.length  intersection  jaccard  cosine  dice  identity
```

`jaccard.similarity` 列を主に見ます。

### 6.4.3 AVG-based vs MAX-based の検定

**平均(AVG)**で見ると誤解しやすい:
- CUN_hap1 vs CKI の "平均" 類似度は、CKI_hap1 と CKI_hap2 の両方との平均
- CUN_hap1 は組み換えによる CKI_hap1 / CKI_hap2 の**モザイク**なので(§2.4)、区間ごとに「より近い方」が入れ替わる。両方を平均してしまうと、その構造がならされて **CKU 側との差が縮み、pedigree のシグナルが埋もれる**

**最大(MAX)**で見るのが本質的:
- MAX(CUN_hap1 vs CKI_hap1, CUN_hap1 vs CKI_hap2) = 「より近い方の CKI hap との類似度」
- 組み換えがあっても CUN_hap1 の全区間は CKI の 2 本のいずれかに由来するので、この MAX は CKU 側の MAX より必ず高くなる

> **注意**: これは「**donor が親 CKI である**」ことの検定であって、「donor hap は CKI_hap1 である」を決める検定ではありません。染色体全体の Jaccard 1 個で親のどちらのハプロタイプかを特定することはできません(§2.4)。区間ごとの由来を追うには §7.5.2 の窓ごとの解析が必要です。

`pg03_qc_graph.sh` は両方を計算します:

```
CUN_hap1 の各parent hap との類似度:
  CKI#1#chr09    Jaccard = 0.62  ← より近い CKI hap
  CKI#2#chr09    Jaccard = 0.35
  CKU#1#chr09  Jaccard = 0.42
  CKU#2#chr09  Jaccard = 0.48

--- Pedigree検定 ---
  AVG-based: NO   (haplotype leakage で撹乱)
  MAX-based: YES  (真の pedigree は整合)
```

### 6.4.4 Haplotype leakage の影響

第4章で見た **CUN の 2 hap が親系統より大きい**という観察(§4.4)が、Jaccard 類似度に影響します。

Jaccard = intersection / union の性質上、**片方の path が長いほど union が大きくなり、Jaccard が下がる**。

CUN_hap1 が 45 Mb (期待の 30 Mb の 1.5 倍)で、他の hap が 30 Mb だと:
- intersection = 25 Mb (共有部分)
- union = 45 + 30 - 25 = 50 Mb
- Jaccard = 0.50

これは「片方が長い」だけで数値が下がっているだけで、実際の生物学的類似性は変わりません。**MAX-based の方が pedigree 検定には適切**です。

### 6.4.5 判定基準

MAX-based で:
- CUN_hap1 が CKI のどちらかと最高類似度 = ✅
- CUN_hap2 が CKU のどちらかと最高類似度 = ✅

両方満たせば **pedigree consistent**。片方でも満たさない場合、hap1/hap2 のラベル入れ替わりや、trio phasing 失敗を疑います。

---

## 6.5 Level 4: SVシグナル

### 6.5.1 VCF の生成

```bash
vg deconstruct -P CUN#1 -H '#' -a -e -t 8 chr09.smooth.final.gfa > chr09.vcf
```

- `-P CUN#1`: 座標基準(reference path)
- `-H '#'`: PanSN separator
- `-a`: nested variants を含める
- `-e`: path traversal ベースで call

### 6.5.2 統計

```bash
# 総 variants
grep -vc '^#' chr09.vcf

# SNPs
grep -v '^#' chr09.vcf | awk 'length($4)==1 && length($5)==1' | wc -l

# SVs (≥50bp)
grep -v '^#' chr09.vcf | awk 'length($4)>=50 || length($5)>=50' | wc -l
```

### 6.5.3 期待値(chr09、~30 Mb)

|Variant type|期待数|
|---|---|
|SNPs|10,000〜30,000|
|Small indels|5,000〜15,000|
|SVs (≥50bp)|100〜500|

**共同研究者の実測**: chr09 で total 432,698 variants(かなり大きい、multi-allelic を含むため)。

### 6.5.4 サイズ分布のチェック

SV サイズ分布が生物学的に妥当かを見ます:

- 1-10 bp: 60-70% (小さい indel が多い)
- 11-100 bp: 20-30%
- 100-1000 bp: 5-10%
- 1-10 kb: 3-5%
- 10 kb+: 1-2% (TE挿入など)

もし特定サイズに極端に偏っていたら、パラメータの問題を疑います。

---

## 6.6 総合判定

### 6.6.1 判定ルール

pg03_qc_graph.sh の総合判定:

|判定|条件|
|---|---|
|**PASS**|Level 1-4 すべて ✅|
|**FAIL**|Level 1 または Level 2 で ❌|
|**WARN**|Level 3 または Level 4 で ⚠️/❌ (graph は使えるが解釈注意)|

### 6.6.2 なぜ Level 3/4 は FAIL でなく WARN か

- Level 1/2 の失敗 = graph 構造の技術的問題(そのままでは使えない)
- Level 3/4 の失敗 = 生物学的解釈や下流解析の疑問(graph 自体は使える)

Level 3 で pedigree 検定が NO でも、**graph そのものは技術的に完成**しています。
Haplotype leakage や phasing の癖など、**データ側の問題**が影響しているケースが多い。

---

## 6.7 失敗パターンと対処

|症状|原因|対処|
|---|---|---|
|圧縮率 > 0.9|**`-B` バグ**|`-B 1G` を指定して再実行|
|path 数 ≠ 6|入力の命名|pg01 を再実行、PanSN 命名を確認|
|path 長差 > 1%|graph 構築失敗|パラメータ再検討|
|SV 数 = 0|VCF 生成失敗|vg deconstruct のオプション確認|
|SV 数 > 100,000|multi-allelic の展開|`-e` オプションを外して確認|
|pedigree AVG-NO, MAX-YES|**haplotype leakage** or Jaccard の長さ感度|OK、これは "解釈上の警告"|
|pedigree MAX-NO|hap1/hap2 のラベル入れ替わり|データ側の再確認|

---

## 6.8 この章のまとめ

- 4 層の QC:構造、保存性、生物学的妥当性、SVシグナル
- Level 1 の圧縮率が最重要指標。`.params.yml` の `transclose-batch` も合わせて確認
- 圧縮率 ≤ 0.35 を目標に
- pedigree 検定は AVG より MAX が本質的
- Level 3/4 の失敗は WARN、graph 自体は使える

次章では、生成された graph の**具体的な解釈・可視化・活用**に進みます。

---

## 参考文献

- odgi documentation: https://pangenome.github.io/odgi.github.io/
- vg deconstruct: https://github.com/vgteam/vg/wiki/vg-deconstruct

---

[← 第5章](05_graph_construction.md) | [第7章: グラフの解釈 →](07_graph_interpretation.md)

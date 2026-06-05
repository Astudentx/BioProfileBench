# BioProfileBench 中文说明

[English README](README.md)

BioProfileBench 是一个用于比较真实微生物丰度谱与预测微生物丰度谱的 Python Benchmark 工具。软件支持单对单和批量比较、分类学层级聚合、仅针对预测值的过滤、稳定过滤编码，以及适合 R/Python 后续分析的结果表。

当前版本：`v1.0.0`

## 核心原则

- 真实值作为固定参考，**不会进行丰度过滤**。
- 所有丰度过滤只作用于预测值。
- 计算指标前，真实值与过滤后的预测值按照 taxa 并集进行对齐。
- 每个样品、每个分类学层级独立计算 Benchmark。
- 多个过滤器都基于同一份未过滤预测矩阵独立计算 mask，最后使用逻辑 **AND** 合并，因此过滤结果不依赖过滤器执行顺序。

## 主要功能

- 使用 `run` 命令完成一个真实值文件与一个预测值文件的比较。
- 使用 manifest TSV 完成多个真实值与多个预测值的批量比较。
- 支持 kingdom、phylum、class、order、family、genus、species 分类学层级聚合。
- 支持预测值的单元格级、taxon 行级、流行率、Top-N 和累计丰度过滤。
- 支持固定阈值网格和基于数据动态生成阈值。
- 使用稳定短编码表示过滤组合，例如 `F000`、`FBCB70CDD`。
- 计算存在/缺失、丰度加权、相关性和微生物组距离指标。
- 输出总汇总表和各分类学层级独立汇总表。

## 依赖环境

- Python 3.9+
- pandas
- numpy
- scipy

## 输入丰度表格式

输入文件必须为 tab 分隔文件。第一列可以是 `ID`、`TaxID`、`SpeciesID`、`FeatureID` 或 `GeneID`。样品列保存丰度值，`Taxonomy` 和/或 `Lineage` 用于分类学聚合。

```text
ID	Sample1	Sample2	Taxonomy	Lineage
gene1	10	0	Escherichia_coli	k__Bacteria;p__Proteobacteria;g__Escherichia;s__Escherichia_coli
gene2	0	5	Klebsiella_pneumoniae	k__Bacteria;p__Proteobacteria;g__Klebsiella;s__Klebsiella_pneumoniae
```

真实值和预测值必须包含相同的样品列。

## 快速运行

### 单对单分析

```bash
python BioProfileBench.py run \
  --truth truth.tsv \
  --pred prediction.tsv \
  --config config.yaml \
  --out results \
  --kind Bacteria
```

`--kind` 是通用的丰度谱类别标签，可以填写 `ARGs`、`VFs`、`Bac`、`Bacteria` 或其他自定义类别。旧参数 `--gene-kind` 仍作为兼容别名保留。

### 批量分析

```bash
python BioProfileBench.py batch \
  --manifest manifest.tsv \
  --config config.yaml \
  --threads 4 \
  --out results
```

推荐的 manifest 格式：

```text
profile_id	dataset	method	kind	prediction_path	truth_path
run_001	Dataset1	Prediction	ARGs	prediction_ARG.tsv	truth_ARG.tsv
run_002	Dataset1	Prediction	Bacteria	prediction_Bac.tsv	truth_Bac.tsv
```

旧版 manifest 中的 `gene_kind` 列仍然可以使用，但推荐改用 `kind`。

## 最小配置示例

```yaml
taxonomy:
  source: Lineage
  levels: [kingdom, phylum, class, order, family, genus, species]
  exclude_unclassified: false

abundance:
  normalize: relative
  presence_threshold: 0

benchmark:
  digits: 6

filtering:
  modes: [raw, common_filter]
  common_filters:
    pred_min_relative_abundance: [1e-06]
```

## 过滤整体流程

软件会针对每个指定分类学层级分别执行以下过程：

1. 将真实值和预测值分别聚合到当前分类学层级。
2. 按照 `abundance.normalize` 分别归一化真实值和预测值。
3. 为当前层级生成预测值过滤配置。
4. 所有过滤 mask 都基于同一份原始预测矩阵独立计算。
5. 使用逻辑 AND 合并所有启用的过滤 mask。
6. 单元格级过滤失败的位置被设置为 0。
7. 行级过滤失败的 taxon 被整行删除。
8. 删除在所有样品中最终丰度总和为 0 的 taxon。
9. 使用未过滤真实值与过滤后预测值的 taxa 并集进行对齐。
10. 对每个样品分别计算 Benchmark 指标。

由于所有 mask 都先独立计算再合并，交换过滤器顺序不会改变结果。但是，多个严格过滤条件进行 AND 组合时，仍然可能删除绝大部分甚至全部预测 taxa。

## 过滤分析模式

```yaml
filtering:
  modes: [raw, common_filter, optimized_filter]
```

- `raw`：不进行预测值过滤，过滤编码为 `F000`。
- `common_filter`：执行配置中的过滤阈值组合，用于统一过滤条件下的比较。
- `optimized_filter`：执行配置中的过滤阈值组合，并将结果标记为用于后续最优阈值选择。

除 `raw` 外，其他模式都会应用解析后的过滤组合。模式名称本身不会改变过滤数学逻辑，也不会参与稳定过滤编码计算。

## 支持的预测值过滤规则

### 单元格级过滤

单元格级过滤独立判断每个 `taxon × sample` 丰度。不通过过滤的单元格会被设置为 0。

#### `pred_min_abundance`

```yaml
pred_min_abundance: [0, 10, 100]
```

仅保留满足以下条件的单元格：

```text
预测丰度 > threshold
```

该过滤使用经过 `abundance.normalize` 处理后的矩阵。因此，当设置 `normalize: relative` 时，该阈值实际作用于相对丰度。

#### `pred_min_relative_abundance`

```yaml
pred_min_relative_abundance: [1e-06, 1e-05, 1e-04]
```

软件会对每个样品中的预测丰度重新计算相对丰度，仅保留满足以下条件的单元格：

```text
相对丰度 > threshold
```

无论主 Benchmark 使用原始丰度还是相对丰度，该过滤器始终按照相对丰度判断。

### Taxon 行级过滤

行级过滤会基于某个 taxon 在所有样品中的表现进行判断。不通过过滤的 taxon 会被整行删除。

#### `pred_min_total_abundance`

```yaml
pred_min_total_abundance: [10, 100]
```

保留跨所有样品丰度总和大于或等于阈值的 taxon。

#### `pred_min_mean_abundance`

```yaml
pred_min_mean_abundance: [0.001, 0.01]
```

保留跨所有样品平均丰度大于或等于阈值的 taxon。

#### `pred_min_max_abundance`

```yaml
pred_min_max_abundance: [0.01, 0.1]
```

保留至少在一个样品中的最大丰度大于或等于阈值的 taxon。

#### `pred_min_prevalence`

```yaml
pred_min_prevalence: [0.1, 0.2, 3]
```

Prevalence 表示预测丰度大于 `abundance.presence_threshold` 的样品数量。

- `value <= 0`：不产生有效 prevalence 限制。
- `0 < value < 1`：按照样品比例解释。例如 `0.2` 表示至少出现在 20% 的样品中。
- `value >= 1`：按照样品数量解释。例如 `3` 表示至少出现在 3 个样品中。

比例转换为样品数时会向上取整。

### Top 丰度过滤

Top 丰度过滤支持两种模式：

```yaml
top_mode: global       # 默认
top_mode: per_sample
```

- `global`：按照 taxon 跨所有样品丰度总和排序，删除不通过的整行 taxon。
- `per_sample`：在每个样品内独立排序，不通过的单元格被设置为 0。

#### `pred_top_k`

```yaml
pred_top_k: [50, 100]
top_mode: [global]
```

保留丰度最高的前 K 个预测 taxa。

#### `pred_top_percent`

```yaml
pred_top_percent: [0.1, 0.2]
top_mode: [global]
```

- 0 到 1 之间的值按照 taxa 总数比例解释。
- 大于 1 的值按照明确 taxa 数量解释。

例如，`0.1` 表示保留约前 10% taxa，`100` 表示保留约前 100 个 taxa。

#### `pred_cumulative_abundance`

```yaml
pred_cumulative_abundance: [0.90, 0.95, 0.99]
top_mode: [global]
```

按照丰度从高到低排序，保留 taxa 直到累计丰度达到指定比例。第一个超过累计阈值的 taxon 也会被保留，从而避免产生空结果。

## 多过滤器组合规则

`common_filters` 中的多个阈值列表会生成笛卡尔积组合。

```yaml
filtering:
  modes: [raw, common_filter]
  common_filters:
    pred_min_relative_abundance: [1e-06, 1e-05]
    pred_min_prevalence: [0.1, 0.2]
```

该配置会生成四个非 raw 过滤组合：

```text
1e-06 + prevalence 0.1
1e-06 + prevalence 0.2
1e-05 + prevalence 0.1
1e-05 + prevalence 0.2
```

在一个过滤组合内部，所有过滤条件都必须通过。单元格级 mask 和 taxon 行级 mask 都基于同一份原始预测矩阵独立计算，最后使用 AND 合并，因此不会出现由过滤顺序导致的差异。

## 动态阈值过滤

动态过滤会基于预测矩阵自动生成阈值，并在每个预测文件、每个分类学层级分别计算。

### 分位数阈值

```yaml
filtering:
  modes: [raw, optimized_filter]
  dynamic_filters:
    pred_min_relative_abundance:
      strategy: quantile
      quantiles: [0.5, 0.75, 0.9, 0.95]
```

使用预测矩阵非零丰度的指定分位数作为阈值。

### 对数间隔阈值

```yaml
filtering:
  dynamic_filters:
    pred_min_abundance:
      strategy: logspace
      min: 1e-06
      max: 1e-02
      count: 5
```

在最小值和最大值之间生成对数间隔阈值。

### 目标 taxa 数量

```yaml
filtering:
  dynamic_filters:
    pred_min_total_abundance:
      strategy: target_taxa_count
      targets: [50, 100, 200]
```

根据 taxon 总丰度排序，生成预计保留指定 taxa 数量的阈值。

### 目标保留丰度比例

```yaml
filtering:
  dynamic_filters:
    pred_min_total_abundance:
      strategy: target_retained_abundance
      targets: [0.90, 0.95, 0.99]
```

根据累计总丰度生成阈值。

### Prevalence 全范围扫描

```yaml
filtering:
  dynamic_filters:
    pred_min_prevalence:
      strategy: prevalence_sweep
```

从 0 到样品总数生成每一个整数 prevalence 阈值。

动态阈值依赖具体预测数据。过滤编码根据最终解析出的数值阈值生成，而不是仅根据动态策略名称生成。

## 稳定过滤编码

每个实际生效的过滤组合都会获得一个稳定短编码。

- `F000`：没有有效过滤。
- 其他编码，例如 `FBCB70CDD`，对应一个明确的过滤参数组合。
- `0`、`0.0` 或省略无作用阈值等价，并会获得相同编码。
- 编码不受文件路径、批次名称、文件夹、分类学层级和分析模式影响。
- 相同的实际过滤参数在不同运行中会获得相同编码。

`filter_trace.tsv` 会记录过滤编码与实际过滤参数的对应关系，因此主结果表可以只使用短编码，同时保持结果可追溯和可重复。

## 过滤诊断文件

`filter_diagnostics.tsv` 用于记录每个过滤器的实际效果：

- `filter_name`：过滤器名称。
- `threshold`：过滤阈值。
- `taxa_before`：过滤前 taxa 数量。
- `taxa_pass`：通过当前单独过滤器的 taxa 数量。
- `taxa_fail`：未通过当前单独过滤器的 taxa 数量。
- `unique_fail_taxa`：仅被当前过滤器删除、但通过其他过滤器的 taxa 数量。
- `pred_retained_abundance`：当前过滤器保留的预测丰度比例。
- `FINAL`：所有过滤 mask 合并后的最终结果。

该文件可以用于判断某个阈值是否没有实际作用、多个阈值是否重复，以及过滤组合是否过于严格。

## Benchmark 指标

BioProfileBench 输出以下主要指标：

- 存在/缺失指标：TP、FP、FN、TN、Precision、Recall、F1、Jaccard、Accuracy、Specificity、MCC、AUPRC。
- 丰度重叠指标：AbundancePrecision、AbundanceRecall、AbundanceF1、WeightedJaccard。
- 相关性和相似性：Spearman、Pearson、Cosine。
- 距离和误差：L1、L2、Bray-Curtis (`BC`)、rooted Jensen-Shannon distance (`rJSD`)、MAE、RMSE。

## 主要输出文件

- `benchmark.tsv`：每个样品 × 分类学层级 × 文件组合 × 过滤编码对应一行。
- `benchmark_long.tsv`：包含更多内部信息的详细 Benchmark 表。
- `benchmark_summary.tsv`：按照 `Kind`、真实值文件、预测值文件、过滤编码和分类学层级汇总。
- `benchmark_summary.<level>.tsv`：每个分类学层级的独立汇总文件。
- `filter_trace.tsv`：过滤编码与实际过滤参数的对应关系。
- `filter_diagnostics.tsv`：单独过滤器及最终组合的过滤效果。
- `pair_filter_ranking.tsv`：跨所有分类学层级的文件与过滤组合排序。
- `pair_filter_level_ranking.tsv`：按照分类学层级分别排序。
- `failed_or_skipped_runs.tsv`：空预测结果或跳过记录。

## 注意事项

- 预测值过滤发生在分类学聚合之后，因此同一个阈值在 species 和 genus 层级可能产生不同效果。
- 动态阈值会针对每个预测文件和每个分类学层级重新计算。
- 严格的 AND 组合可能产生空预测结果，软件会记录这些结果，而不是静默丢弃。
- 如果 `Taxonomy` 不包含带前缀的分类层级信息，软件通常只能将其作为 species 标签使用。多层级分析建议使用 `Lineage`。

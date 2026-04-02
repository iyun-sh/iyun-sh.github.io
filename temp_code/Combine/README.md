# Combine - 多数据集合并与批次效应校正

使用 ComBat（`sva` 包）对多个 GEO 等来源的表达矩阵进行批次效应校正，输出合并后的表达矩阵和分组文件。

## 快速使用

```bash
bash run_combine.sh \
  -d /path/GSE10072_dat.csv,/path/GSE19804_dat.csv \
  -g /path/GSE10072_group.csv,/path/GSE19804_group.csv \
  -n GSE10072,GSE19804 \
  -o /media/share/output/project/combine
```

## 参数说明

| 参数 | 短参数 | 必需 | 默认值 | 说明 |
|------|--------|------|--------|------|
| `--data-files` | `-d` | 是 | - | 表达矩阵 CSV，逗号分隔（至少 2 个） |
| `--group-files` | `-g` | 是 | - | 分组文件 CSV，逗号分隔 |
| `--output-dir` | `-o` | 是 | `./output` | 输出目录 |
| `--dataset-names` | `-n` | 否 | 文件名 | 数据集名称，逗号分隔 |
| `--output-name` | - | 否 | `combat` | 输出文件前缀名 |
| `--gene-col` | - | 否 | `SYMBOL` | 基因列名 |
| `--width` | `-w` | 否 | 9 | 图片宽度（英寸） |
| `--height` | - | 否 | 7 | 图片高度（英寸） |
| `--rm` | - | 否 | - | 删除输出目录 |

## 输入文件格式

### 表达矩阵 CSV（`dat.csv`）

```
SYMBOL,sample1,sample2,sample3,...
TP53,5.23,4.89,6.12,...
BRCA1,3.41,3.56,3.78,...
```

- 第一列为基因名（默认列名 `SYMBOL`，可通过 `--gene-col` 指定）
- 其余列为样本的表达值

### 分组文件 CSV（`group.csv`）

```
sample,group
sample1,Tumor
sample2,Normal
sample3,Tumor
```

- 必须包含 `sample` 和 `group` 两列
- `sample` 列值必须与表达矩阵的列名一致

## 输出文件

| 文件 | 说明 |
|------|------|
| `combat.dat.csv` | 校正后表达矩阵（按组排序） |
| `combat.group.csv` | 合并后分组文件（含 dataset 列） |
| `00.dataset_info.csv` | 数据集统计信息 |
| `01.pca_before_group.pdf/png` | 校正前 PCA（按分组着色） |
| `02.pca_before_dataset.pdf/png` | 校正前 PCA（按数据集着色） |
| `03.boxplot_before.pdf/png` | 校正前箱线图 |
| `04.pca_after_group.pdf/png` | 校正后 PCA（按分组着色） |
| `05.pca_after_dataset.pdf/png` | 校正后 PCA（按数据集着色） |
| `06.boxplot_after.pdf/png` | 校正后箱线图 |

## 生成报告

```bash
Rscript report_templates/generate_report.R \
  -d /path/to/combine_20250324_xxxx \
  -p project_name
```

## R 依赖

- sva, data.table, dplyr, tibble, ggplot2, ggrepel, RColorBrewer, optparse

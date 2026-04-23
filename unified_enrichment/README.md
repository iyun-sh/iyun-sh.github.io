# 统一富集分析 Unified Enrichment

**项目目录**: `/media/desk16/share/secure/unified_enrichment`

> **重要**: 每次进行富集分析前，请确认项目目录和输入文件路径。

---

## 目录结构

```
unified_enrichment/
├── scripts/
│   ├── main.R                      # 主脚本
│   ├── utils.R                     # 工具函数
│   ├── qc.R                       # 质控模块
│   ├── go_kegg.R                  # GO/KEGG分析
│   ├── gsea.R                     # GSEA分析
│   ├── gsea_cli.R                 # GSEA独立CLI工具
│   ├── generate_enrichment_report.R # 中文Word报告生成
│   ├── report.R                   # 报告生成（旧版）
│   └── requirements.R              # 依赖管理
├── run_unified_enrichment.sh      # 一键运行脚本
├── test_data/                     # 测试数据
├── logs/                          # 日志目录（自动创建）
├── results/                       # 结果目录（自动创建）
└── report/                        # 报告目录（自动创建）
```

**本仓库说明**：`scripts/`、`templates/`、`resources/`、`report_templates/`、`test_data/`（示例输入 CSV）与 **`scripts/config/unified_enrichment.config.ini`** 为完整内容；**`results/`、`report/`、`logs/` 在仓库内仅占位（`.gitkeep`）**，本地运行生成文件由 `.gitignore` 忽略。

---

## 一键运行

```bash
cd /media/desk16/share/secure/unified_enrichment

# 完整流程：GO/KEGG + GSEA 同时运行
./run_unified_enrichment.sh \
  -w test_data \
  -s hsa \
  --gene-list gene_list.csv \
  --expr-matrix GSE126124.dat.csv \
  --target-genes target_genes.csv \
  --ranked-list ranked_gene_list.csv \
  --gmt ../resources/gmt/hallmark.gmt \
  --go-kegg true \
  --gsea true \
  --output-dir test_data/results

# GSEA-only 模式：仅运行 GSEA（跳过 GO/KEGG）
./run_unified_enrichment.sh \
  -w test_data \
  -s hsa \
  --analysis-mode gsea-only \
  --ranked-list ranked_gene_list.csv \
  --gmt ../resources/gmt/hallmark.gmt \
  --output-dir test_data/results_gsea_only
```

> **提示**: 使用 `--analysis-mode gsea-only` 可仅运行 GSEA 分析，跳过 GO/KEGG。适合单细胞分析后已有 GO/KEGG 结果的场景。

---

## 参数说明

### 必需参数

| 参数 | 短选项 | 说明 | 必需 |
|------|--------|------|------|
| `-w, --work-dir` | DIR | 工作目录 | 是 |
| `-s, --species` | CODE | 物种代码 (hsa/mmu/rno) | 是 |
| `--gene-list` | FILE | 基因列表文件 | 是* |
| `--expr-matrix` | FILE | 表达矩阵文件 | 是* |
| `--target-genes` | FILE | 目标基因文件 | 是* |
| `--ranked-list` | FILE | 排序基因列表文件 | 是** |
| `--gmt` | FILE | GMT基因集文件 | 是** |

* 仅完整模式需要
** GSEA 模式需要

### 分析开关

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--analysis-mode` | full | 分析模式 (full/gsea-only) |
| `--go-kegg` | true | GO/KEGG分析 (true/false) |
| `--gsea` | true | GSEA分析 (true/false) |
| `--gsea-method` | prerank | GSEA方法 (prerank/correlation/median_split) |

### GO/KEGG参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--go-ont` | ALL | GO本体 (BP/MF/CC/ALL) |
| `--go-pvalue` | 0.05 | p值阈值 |
| `--go-qvalue` | 0.2 | q值阈值 |
| `--go-padj` | BH | 校正方法 |

### GSEA参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--gsea-correlation` | spearman | 相关方法 |
| `--gsea-rank-metric` | log2FC | 排序指标 |
| `--gsea-min-size` | 10 | 最小基因集大小 |
| `--gsea-max-size` | 500 | 最大基因集大小 |
| `--gsea-pvalue` | 0.05 | p值阈值 |
| `--gsea-padj` | BH | 校正方法 |
| `--gsea-top-n` | 5 | Top通路数 |

### 其他参数

| 参数 | 短选项 | 默认值 | 说明 |
|------|--------|--------|------|
| `--output-dir` | - | results | 输出目录 |
| `--log-dir` | - | logs | 日志目录 |
| `--seed` | - | 123 | 随机种子 |
| `--plot-width` | - | 10 | 图宽度 |
| `--plot-height` | - | 8 | 图高度 |

---

## 输入文件格式

### 基因列表 (`--gene-list`)
CSV 格式：
```csv
gene
TP53
MYC
EGFR
```

### 表达矩阵 (`--expr-matrix`)
CSV 格式，第一列为基因名，其余列为样本：
```
SYMBOL,sample1,sample2,sample3,...
GeneA,5.2,6.1,4.8,...
GeneB,3.4,3.9,4.2,...
```

### 目标基因 (`--target-genes`)
CSV 格式：
```csv
gene
TP53
MYC
```

### 排序基因列表 (`--ranked-list`)
CSV 格式，用于 prerank GSEA：
```csv
gene,logFC
TP53,2.5
MYC,2.1
EGFR,1.8
```

### GMT 基因集 (`--gmt`)
从 [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb/) 下载：
```
gene_set_name    description    gene1    gene2    gene3    ...
HALLMARK_ADIPOGENESIS    Adipogenesis    TP53    MYC    ...
```

---

## 输出文件

运行后自动生成：

### logs/ - 日志
```
unified_enrichment.20260317_160354.log
sessionInfo_20260317_160354.txt
```
包含：时间戳、输入参数、统计信息、警告信息、输出文件列表

### results/unified_enrichment/ - 结果
```
go/
  ├── go_enrichment_results.csv     # GO富集结果
  └── go_enrichment_plot.pdf/png   # GO可视化
kegg/
  ├── kegg_enrichment_results.csv  # KEGG富集结果
  └── kegg_enrichment_plot.pdf/png # KEGG可视化
gsea/
  ├── prerank_gsea_results.csv     # GSEA结果
  └── ...
plots/
  ├── PreRank_GSEA_EnrichmentPlot_SCI.pdf/png  # GSEA富集图
  └── ...
qc/
  └── PreRank_GSEA_QC_Report.csv  # 质控报告
unified_enrichment.config.ini      # 配置文件
<dataset>_enrichment_node.qs2      # 分析节点快照
<dataset>_sessionInfo.<timestamp>.qs2  # sessionInfo二进制快照
```

### report/ - 报告（中文 Word 报告）
```
enrichment_report_20260317_160354.docx  # 中文Word报告
```
模板文件（可选）：
```
templates/unified_enrichment_report_template.docx
```
报告包含：
- 方法说明
- 汇总统计
- GO 富集分析结果表格和可视化（图1）
- KEGG 通路分析结果表格和可视化（图2）
- GSEA 分析结果表格和可视化（图3）

> **注意**: 报告内容会根据分析模式自动调整：
> - 完整模式：显示 GO/KEGG + GSEA
> - GSEA-only 模式：仅显示 GSEA

---

## GSEA 方法说明

| 方法 | 适用场景 | 输入文件 |
|------|----------|----------|
| prerank | 差异表达分析后 | ranked_gene_list.csv |
| correlation | 单基因功能探索 | expression_matrix + target_genes |
| median_split | 分组比较分析 | expression_matrix + target_genes |

---

## 依赖包

```r
install.packages(c("optparse", "data.table", "dplyr", "tidyverse", "ggplot2",
                  "clusterProfiler", "DOSE", "enrichplot", "GSEA",
                  "org.Hs.eg.db", "org.Mm.eg.db", "org.Rn.eg.db",
                  "stats", "grid", "gridExtra", "officer", "flextable"))
```

---

## 测试数据

测试数据位于 `test_data/` 目录：

| 文件 | 说明 |
|------|------|
| `gene_list.csv` | 差异基因列表 (200个基因) |
| `target_genes.csv` | 目标基因 |
| `ranked_gene_list.csv` | 排序基因列表 (500个基因，上下调各250个) |
| `GSE126124.dat.csv` | 表达矩阵 |
| `resources/gmt/hallmark.gmt` | GMT基因集 (Hallmark, 50个通路) |

---

## 追溯关系

通过时间戳可追溯：
- `logs/unified_enrichment.YYYYMMDD_HHMMSS.log` → 日志
- `results/unified_enrichment/.YYYYMMDD_HHMMSS` → 追踪文件
- `results/unified_enrichment/unified_enrichment.config.ini` → 配置文件
- `report/enrichment_report_YYYYMMDD_HHMMSS.docx` → 中文Word报告

---

## 常见问题

**Q: 报错 "GO/KEGG基因列表文件不存在"**
A: 检查 `--gene-list` 参数指定的文件路径是否正确

**Q: 报错 "GMT文件不存在"**
A: 检查 `--gmt` 参数指定的GMT文件路径是否正确，可从MSigDB下载

**Q: GSEA报错 "缺少输入文件"**
A: 根据GSEA方法检查所需文件是否都提供：
  - prerank: 需要 `--ranked-list`
  - correlation/median_split: 需要 `--expr-matrix` 和 `--target-genes`

**Q: 报告中 GSEA 部分显示为空**
A: 检查是否使用了 GSEA-only 模式，该模式下报告不会显示 GO/KEGG 结果（这是正确的行为）

**Q: GSEA 匹配率很低 / 无显著通路**
A: Pre-rank 与 correlation/median_split 流程会在运行前自动做基因 ID 对齐：ENSEMBL 去版本号、SYMBOL 统一大写，并在排序表与 GMT 分别为 SYMBOL/ENSEMBL 时用 `AnnotationDbi` 转换。请在配置中正确设置 `species`（hsa/mmu/rno）。`PreRank_GSEA_QC_Report.csv` 中的 `GMT_Match_Rate` 为**对齐后**「排序表基因与 GMT 的交集 / 输入基因数」；`ID_Converted` 表示是否做了 SYMBOL↔ENSEMBL 转换。

---

## 分析模式说明

### 完整模式 (full)
- 同时运行 GO/KEGG 和 GSEA 分析
- 报告包含所有分析结果

### GSEA-only 模式 (gsea-only)
- 仅运行 GSEA 分析
- 跳过 GO/KEGG 分析
- 适合场景：
  - 单细胞分析后已有 GO/KEGG 结果
  - 仅需要 GSEA 分析
  - 其他流程已生成 GO/KEGG 结果

#!/usr/bin/env Rscript
# ==============================================================================
# unified_enrichment/scripts/test_basic.R
# 基础功能测试脚本
# ==============================================================================

# 设置工作目录
setwd("/media/desk16/iyunlyl/standardized_workflow/enrichment/unified_enrichment")

# 加载工具函数
source("scripts/utils.R")

cat("============================================================\n")
cat("        统一富集分析流程 - 基础功能测试\n")
cat("============================================================\n\n")

# 测试1: 配置读取
cat("测试1: 配置读取...\n")
cfg <- tryCatch({
  read_config("config/config.yaml")
}, error = function(e) {
  cat(sprintf("  [失败] %s\n", e$message))
  NULL
})
if (!is.null(cfg)) {
  cat("  [成功] 配置文件读取成功\n")
  cat(sprintf("  项目名称: %s\n", cfg$project$name))
  cat(sprintf("  物种: %s\n", cfg$project$species))
  cat(sprintf("  GO/KEGG启用: %s\n", cfg$analysis$go_kegg$enabled))
  cat(sprintf("  GSEA启用: %s\n", cfg$analysis$gsea$enabled))
}

# 测试2: 目录创建
cat("\n测试2: 目录创建...\n")
test_dir <- "test_temp"
tryCatch({
  ensure_dir(test_dir)
  if (dir.exists(test_dir)) {
    cat("  [成功] 目录创建成功\n")
    unlink(test_dir, recursive = TRUE)
  }
}, error = function(e) {
  cat(sprintf("  [失败] %s\n", e$message))
})

# 测试3: 基因列表读取
cat("\n测试3: 基因列表读取...\n")
tryCatch({
  genes <- read_gene_list("input/gene_list.csv")
  cat(sprintf("  [成功] 读取 %d 个基因\n", length(genes)))
  cat(sprintf("  前5个基因: %s\n", paste(head(genes, 5), collapse = ", ")))
}, error = function(e) {
  cat(sprintf("  [失败] %s\n", e$message))
})

# 测试4: 目标基因读取
cat("\n测试4: 目标基因读取...\n")
tryCatch({
  target_genes <- read_gene_list("input/target_genes.csv")
  cat(sprintf("  [成功] 读取 %d 个目标基因\n", length(target_genes)))
  cat(sprintf("  目标基因: %s\n", paste(target_genes, collapse = ", ")))
}, error = function(e) {
  cat(sprintf("  [失败] %s\n", e$message))
})

# 测试5: 表达矩阵读取
cat("\n测试5: 表达矩阵读取...\n")
tryCatch({
  exp <- read.csv("input/expression_matrix.csv", header = TRUE, row.names = 1, check.names = FALSE)
  cat(sprintf("  [成功] 表达矩阵: %d 基因 x %d 样本\n", nrow(exp), ncol(exp)))
}, error = function(e) {
  cat(sprintf("  [失败] %s\n", e$message))
})

# 测试6: GMT文件检查
cat("\n测试6: GMT文件检查...\n")
gmt_file <- "input/gene_sets.gmt"
if (file.exists(gmt_file)) {
  cat("  [成功] GMT文件存在\n")
  # 读取前几行
  lines <- readLines(gmt_file, n = 3)
  cat(sprintf("  GMT文件包含 %d 行（示例）\n", length(lines)))
} else {
  cat("  [警告] GMT文件不存在，需要下载\n")
}

# 测试7: 物种推断
cat("\n测试7: 物种推断...\n")
tryCatch({
  genes <- read_gene_list("input/gene_list.csv")
  expr_info <- infer_expr_id(genes)
  cat(sprintf("  [成功] 推断物种: %s, ID类型: %s\n", expr_info$species, expr_info$id_type))
}, error = function(e) {
  cat(sprintf("  [失败] %s\n", e$message))
})

# 测试8: 路径解析
cat("\n测试8: 路径解析...\n")
tryCatch({
  paths <- resolve_paths(cfg, ".")
  cat("  [成功] 路径解析成功\n")
  cat(sprintf("  结果目录: %s\n", paths$results_dir))
  cat(sprintf("  GO目录: %s\n", paths$go_dir))
  cat(sprintf("  KEGG目录: %s\n", paths$kegg_dir))
  cat(sprintf("  GSEA目录: %s\n", paths$gsea_dir))
}, error = function(e) {
  cat(sprintf("  [失败] %s\n", e$message))
})

cat("\n============================================================\n")
cat("                    基础测试完成\n")
cat("============================================================\n")

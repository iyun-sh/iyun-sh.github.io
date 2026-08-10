#!/usr/bin/env Rscript
# ==============================================================================
# unified_enrichment/scripts/run_tests.R
# 运行三种测试场景
# ==============================================================================

# 设置工作目录
setwd("/media/desk16/iyunlyl/standardized_workflow/enrichment/unified_enrichment")

# 加载依赖
source("scripts/utils.R")
source("scripts/qc.R")

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(enrichplot)
  library(ggplot2)
  library(viridis)
  library(data.table)
  library(dplyr)
  library(tidyr)
})

cat("============================================================\n")
cat("        统一富集分析流程 - 多场景测试\n")
cat("============================================================\n\n")

# 初始化日志
log_file <- "logs/test_scenarios.log"
ensure_dir("logs")
init_logger(log_file)

# ---- 测试1：仅GO/KEGG分析 ----
cat("\n>>> 测试1: 仅GO/KEGG分析 <<<\n")
log_msg("info", "开始测试1: 仅GO/KEGG分析")

tryCatch({
  # 读取基因列表
  gene_list <- read_gene_list("input/gene_list.csv")
  log_msg("info", sprintf("读取基因列表: %d 个基因", length(gene_list)))

  # 运行基因Symbol质控
  qc_result <- symbol_qc("input/gene_list.csv", "qc/test1", "hsa")
  log_msg("info", sprintf("质控完成: 匹配基因 %d/%d", length(qc_result$matched_genes), qc_result$total_genes))

  # 加载物种数据库
  library(org.Hs.eg.db)
  OrgDb <- org.Hs.eg.db

  # 基因ID转换
  gene_transform <- clusterProfiler::bitr(
    gene_list,
    fromType = "SYMBOL",
    toType = c("ENTREZID", "ENSEMBL"),
    OrgDb = OrgDb,
    drop = FALSE
  )
  log_msg("info", sprintf("基因ID转换: %d 成功", sum(!is.na(gene_transform$ENTREZID))))

  # GO富集分析
  ego <- enrichGO(
    gene = gene_transform$ENTREZID,
    OrgDb = OrgDb,
    keyType = "ENTREZID",
    ont = "ALL",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
  )

  if (!is.null(ego) && nrow(ego@result) > 0) {
    log_msg("info", sprintf("GO分析完成: %d 个显著条目", nrow(ego@result)))

    # 保存结果
    ensure_dir("results_test1/go")
    fwrite(as.data.frame(ego@result), "results_test1/go/go_enrichment_results.csv")
    cat(sprintf("  [成功] GO分析: %d 个显著条目\n", nrow(ego@result)))

    # 绘图
    ensure_dir("results_test1/plots")
    if (nrow(ego@result) > 0) {
      p <- dotplot(ego, showCategory = 5) + theme_bw()
      ggsave("results_test1/plots/go_dotplot.pdf", p, width = 10, height = 8)
      ggsave("results_test1/plots/go_dotplot.png", p, width = 10, height = 8, dpi = 300)
      cat("  [成功] GO图形已保存\n")
    }
  } else {
    cat("  [警告] 未找到显著GO条目\n")
  }

  # KEGG富集分析
  kk <- enrichKEGG(
    gene = gene_transform$ENTREZID,
    keyType = "kegg",
    organism = "hsa",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )

  if (!is.null(kk) && nrow(kk@result) > 0) {
    kk <- setReadable(kk, OrgDb = OrgDb, keyType = "ENTREZID")
    log_msg("info", sprintf("KEGG分析完成: %d 个显著通路", nrow(kk@result)))

    ensure_dir("results_test1/kegg")
    fwrite(as.data.frame(kk@result), "results_test1/kegg/kegg_enrichment_results.csv")
    cat(sprintf("  [成功] KEGG分析: %d 个显著通路\n", nrow(kk@result)))

    # 绘图
    if (nrow(kk@result) > 0) {
      p <- dotplot(kk, showCategory = 5) + theme_bw()
      ggsave("results_test1/plots/kegg_dotplot.pdf", p, width = 10, height = 8)
      ggsave("results_test1/plots/kegg_dotplot.png", p, width = 10, height = 8, dpi = 300)
      cat("  [成功] KEGG图形已保存\n")
    }
  } else {
    cat("  [警告] 未找到显著KEGG通路\n")
  }

  log_msg("info", "测试1完成")
  cat(">>> 测试1完成 <<<\n")

}, error = function(e) {
  log_msg("error", sprintf("测试1失败: %s", e$message))
  cat(sprintf("  [失败] %s\n", e$message))
})

# ---- 测试2：仅GSEA分析 ----
cat("\n>>> 测试2: 仅GSEA分析 <<<\n")
log_msg("info", "开始测试2: 仅GSEA分析")

tryCatch({
  # 读取目标基因
  target_genes <- read_gene_list("input/target_genes.csv")
  log_msg("info", sprintf("读取目标基因: %d 个", length(target_genes)))

  # 读取表达矩阵
  exp <- read.csv("input/expression_matrix.csv", header = TRUE, row.names = 1, check.names = FALSE)
  exp <- t(exp)  # 转置：行为样本，列为基因
  log_msg("info", sprintf("表达矩阵: %d 样本 x %d 基因", nrow(exp), ncol(exp)))

  # 检查目标基因是否存在
  present_genes <- intersect(target_genes, colnames(exp))
  log_msg("info", sprintf("存在于矩阵中的目标基因: %d/%d", length(present_genes), length(target_genes)))

  # 读取GMT文件
  gmt_data <- read.gmt("input/gene_sets.gmt")
  all_gmt_genes <- unique(gmt_data$gene)
  log_msg("info", sprintf("GMT基因集: %d 个基因", length(all_gmt_genes)))

  # 计算匹配率
  matched_genes <- intersect(all_gmt_genes, colnames(exp))
  match_rate <- length(matched_genes) / length(all_gmt_genes) * 100
  log_msg("info", sprintf("GMT匹配率: %.2f%%", match_rate))

  # 对每个目标基因进行GSEA
  ensure_dir("results_test2/gsea")
  ensure_dir("results_test2/plots")

  for (gene in present_genes) {
    log_msg("info", sprintf("处理基因: %s", gene))

    # 计算相关性
    cor_values <- cor(exp[, gene], exp, method = "spearman", use = "pairwise.complete.obs")
    cor_values <- as.numeric(cor_values[1, ])
    names(cor_values) <- colnames(exp)
    cor_values <- cor_values[names(cor_values) != gene]

    # 处理重复值
    if (length(unique(cor_values)) < length(cor_values)) {
      cor_values <- cor_values + rnorm(length(cor_values), mean = 0, sd = 1e-8)
    }
    gene_rank <- sort(cor_values, decreasing = TRUE)

    # GSEA分析
    gsea_result <- GSEA(
      geneList = gene_rank,
      exponent = 1,
      minGSSize = 10,
      maxGSSize = 500,
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      TERM2GENE = gmt_data,
      eps = 0
    )

    if (!is.null(gsea_result) && nrow(gsea_result@result) > 0) {
      log_msg("info", sprintf("%s GSEA完成: %d 个显著通路", gene, nrow(gsea_result@result)))

      # 保存结果
      fwrite(as.data.frame(gsea_result@result),
             file.path("results_test2/gsea", paste0(gene, "_gsea_results.csv")))

      # 绘图
      if (nrow(gsea_result@result) >= 1) {
        top_ids <- head(gsea_result@result$ID, min(3, nrow(gsea_result@result)))
        p <- gseaplot2(gsea_result, geneSetID = top_ids, color = "firebrick")
        ggsave(file.path("results_test2/plots", paste0(gene, "_gsea_plot.pdf")),
               p, width = 12, height = 8)
        ggsave(file.path("results_test2/plots", paste0(gene, "_gsea_plot.png")),
               p, width = 12, height = 8, dpi = 150)
        cat(sprintf("  [成功] %s GSEA: %d 个显著通路\n", gene, nrow(gsea_result@result)))
      }
    } else {
      log_msg("warn", sprintf("%s GSEA未找到显著通路", gene))
      cat(sprintf("  [警告] %s 未找到显著GSEA通路\n", gene))
    }
  }

  log_msg("info", "测试2完成")
  cat(">>> 测试2完成 <<<\n")

}, error = function(e) {
  log_msg("error", sprintf("测试2失败: %s", e$message))
  cat(sprintf("  [失败] %s\n", e$message))
})

# ---- 测试3：同时运行GO/KEGG和GSEA ----
cat("\n>>> 测试3: 同时运行GO/KEGG和GSEA <<<\n")
log_msg("info", "开始测试3: 同时运行GO/KEGG和GSEA")

tryCatch({
  # 合并测试1和测试2的功能
  ensure_dir("results_test3/go")
  ensure_dir("results_test3/kegg")
  ensure_dir("results_test3/gsea")
  ensure_dir("results_test3/plots")
  ensure_dir("qc/test3")

  # 1. GO/KEGG分析
  gene_list <- read_gene_list("input/gene_list.csv")
  qc_result <- symbol_qc("input/gene_list.csv", "qc/test3", "hsa")

  library(org.Hs.eg.db)
  OrgDb <- org.Hs.eg.db

  gene_transform <- clusterProfiler::bitr(
    gene_list,
    fromType = "SYMBOL",
    toType = c("ENTREZID"),
    OrgDb = OrgDb,
    drop = FALSE
  )

  # GO分析
  ego <- enrichGO(
    gene = gene_transform$ENTREZID,
    OrgDb = OrgDb,
    keyType = "ENTREZID",
    ont = "ALL",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    readable = TRUE
  )

  if (!is.null(ego) && nrow(ego@result) > 0) {
    fwrite(as.data.frame(ego@result), "results_test3/go/go_enrichment_results.csv")
    p <- dotplot(ego, showCategory = 5) + theme_bw()
    ggsave("results_test3/plots/go_dotplot.png", p, width = 10, height = 8, dpi = 300)
    cat(sprintf("  [成功] GO分析: %d 个显著条目\n", nrow(ego@result)))
  }

  # KEGG分析
  kk <- enrichKEGG(
    gene = gene_transform$ENTREZID,
    keyType = "kegg",
    organism = "hsa",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05
  )

  if (!is.null(kk) && nrow(kk@result) > 0) {
    kk <- setReadable(kk, OrgDb = OrgDb, keyType = "ENTREZID")
    fwrite(as.data.frame(kk@result), "results_test3/kegg/kegg_enrichment_results.csv")
    p <- dotplot(kk, showCategory = 5) + theme_bw()
    ggsave("results_test3/plots/kegg_dotplot.png", p, width = 10, height = 8, dpi = 300)
    cat(sprintf("  [成功] KEGG分析: %d 个显著通路\n", nrow(kk@result)))
  }

  # 2. GSEA分析
  target_genes <- read_gene_list("input/target_genes.csv")
  exp <- read.csv("input/expression_matrix.csv", header = TRUE, row.names = 1, check.names = FALSE)
  exp <- t(exp)

  present_genes <- intersect(target_genes, colnames(exp))
  gmt_data <- read.gmt("input/gene_sets.gmt")

  for (gene in present_genes) {
    cor_values <- cor(exp[, gene], exp, method = "spearman", use = "pairwise.complete.obs")
    cor_values <- as.numeric(cor_values[1, ])
    names(cor_values) <- colnames(exp)
    cor_values <- cor_values[names(cor_values) != gene]

    if (length(unique(cor_values)) < length(cor_values)) {
      cor_values <- cor_values + rnorm(length(cor_values), mean = 0, sd = 1e-8)
    }
    gene_rank <- sort(cor_values, decreasing = TRUE)

    gsea_result <- GSEA(
      geneList = gene_rank,
      exponent = 1,
      minGSSize = 10,
      maxGSSize = 500,
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      TERM2GENE = gmt_data,
      eps = 0
    )

    if (!is.null(gsea_result) && nrow(gsea_result@result) > 0) {
      fwrite(as.data.frame(gsea_result@result),
             file.path("results_test3/gsea", paste0(gene, "_gsea_results.csv")))
      cat(sprintf("  [成功] %s GSEA: %d 个显著通路\n", gene, nrow(gsea_result@result)))
    }
  }

  log_msg("info", "测试3完成")
  cat(">>> 测试3完成 <<<\n")

}, error = function(e) {
  log_msg("error", sprintf("测试3失败: %s", e$message))
  cat(sprintf("  [失败] %s\n", e$message))
})

# 输出测试摘要
cat("\n============================================================\n")
cat("                    测试摘要\n")
cat("============================================================\n")

# 检查输出文件
test_dirs <- c("results_test1", "results_test2", "results_test3")
for (td in test_dirs) {
  if (dir.exists(td)) {
    files <- list.files(td, recursive = TRUE)
    cat(sprintf("\n%s:\n", td))
    if (length(files) > 0) {
      cat(paste("  ", files, collapse = "\n"), "\n")
    } else {
    cat("  (空)\n")
    }
  }
}

cat("\n============================================================\n")
cat("                    所有测试完成\n")
cat("============================================================\n")

# ==============================================================================
# unified_enrichment/scripts/qc.R
# 质量控制模块
# 合并自 GO_KEGG/code/qc_symbol.R 和 GSEA的QC逻辑
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
})

source("scripts/utils.R")

# ---- 基因Symbol质量控制 ----

#' 基因Symbol质量控制（用于GO/KEGG）
#' @param input_file 输入基因文件
#' @param output_dir 输出目录
#' @param species 物种代码 (hsa/mmu)
#' @return 质控结果列表
symbol_qc <- function(input_file, output_dir = NULL, species = "hsa") {

  if (!file.exists(input_file)) {
    stop(sprintf("找不到输入文件: %s", input_file))
  }

  # 生成输出文件名
  file_basename <- basename(input_file)
  file_prefix <- tools::file_path_sans_ext(file_basename)
  output_file <- paste0("qc_", file_prefix, "_report.txt")

  if (!is.null(output_dir)) {
    ensure_dir(output_dir)
    output_file <- file.path(output_dir, basename(output_file))
  }

  # 读取输入文件
  gene_df <- read.csv(input_file, header = TRUE, stringsAsFactors = FALSE)
  gene_col <- colnames(gene_df)[1]

  if (is.na(gene_col) || gene_col == "") {
    stop("输入文件中未找到基因列")
  }

  gene_symbol <- unique(na.omit(as.character(gene_df[[gene_col]])))
  total_gene_num <- length(gene_symbol)

  if (total_gene_num == 0) {
    stop("基因列为空，请检查输入文件")
  }

  # 加载物种数据库
  OrgDb <- get_org_db(species)
  target_genes <- AnnotationDbi::keys(OrgDb, keytype = "SYMBOL")
  species_name <- get_species_name(species)

  # 匹配基因
  is_matched <- gene_symbol %in% target_genes
  matched <- gene_symbol[is_matched]
  unmatched <- gene_symbol[!is_matched]

  # 统计结果
  stat_result <- data.frame(
    Metric = c("总基因数", paste0(species_name, "基因"), "未匹配基因"),
    Count = c(total_gene_num, length(matched), length(unmatched)),
    Rate = c(
      "100%",
      sprintf("%.2f%%", length(matched) / total_gene_num * 100),
      sprintf("%.2f%%", length(unmatched) / total_gene_num * 100)
    ),
    stringsAsFactors = FALSE
  )

  # 添加注释列
  gene_df_annotated <- gene_df %>%
    mutate(
      Gene_Species = case_when(
        .[[gene_col]] %in% matched ~ species_name,
        TRUE ~ "未匹配"
      )
    )

  # 保存质控报告
  sink(output_file, append = FALSE)
  on.exit(sink())
  cat(sprintf("正在读取文件: %s\n", input_file))
  cat(sprintf("正在加载基因注释数据库: %s\n", species_name))
  cat("\n========== 基因匹配统计结果 ==========\n")
  print(stat_result, row.names = FALSE)
  cat("========================================\n\n")
  if (length(unmatched) > 0) {
    cat(sprintf("未匹配基因列表 (共%d个):\n", length(unmatched)))
    cat(paste(unmatched, collapse = ", "), "\n\n")
  }

  log_msg("info", sprintf("质控报告已保存: %s", output_file))

  invisible(list(
    annotated_data = gene_df_annotated,
    stats = stat_result,
    matched_genes = matched,
    unmatched = unmatched,
    species = species_name,
    total_genes = total_gene_num
  ))
}

# ---- 表达矩阵质量控制 ----

#' 表达矩阵质量控制（用于GSEA）
#' @param exp_file 表达矩阵文件
#' @param variance_threshold 方差阈值
#' @param na_threshold NA比例阈值
#' @param key_genes 目标基因（受保护，不被过滤）
#' @return 质控后的数据列表
qc_expression_matrix <- function(exp_file, variance_threshold = 0.01,
                                  na_threshold = 0.5, key_genes = NULL) {

  if (!file.exists(exp_file)) {
    stop(sprintf("表达矩阵文件不存在: %s", exp_file))
  }

  # 读取表达矩阵（行=基因，列=样本）
  exp <- fread(exp_file, header = TRUE, check.names = FALSE, encoding = "UTF-8")
  # 将第一列作为行名
  gene_col <- colnames(exp)[1]
  exp <- as.data.frame(exp)
  rownames(exp) <- exp[[gene_col]]
  exp <- exp[, -1, drop = FALSE]
  exp <- as.matrix(exp)

  raw_exp_gene_count <- nrow(exp)

  # 检查目标基因是否存在
  if (!is.null(key_genes) && length(key_genes) > 0) {
    present_key_genes <- intersect(key_genes, rownames(exp))
    missing_key_genes <- setdiff(key_genes, rownames(exp))

    if (length(missing_key_genes) > 0) {
      log_msg("warn", sprintf("目标基因未在表达矩阵找到: %s", paste(missing_key_genes, collapse = ", ")))
    }
    if (length(present_key_genes) == 0) {
      stop("目标基因在表达矩阵中一个都匹配不到")
    }
  } else {
    present_key_genes <- character(0)
  }

  # 计算基因方差和NA比例（行=基因）
  gene_variance <- apply(exp, 1, var, na.rm = TRUE)
  gene_na_ratio <- rowSums(is.na(exp)) / ncol(exp)

  # 识别低质量基因
  bad_genes_detected <- names(which((gene_variance < variance_threshold) | (gene_na_ratio > na_threshold)))
  protected_from_filter <- intersect(bad_genes_detected, present_key_genes)
  bad_genes <- setdiff(bad_genes_detected, present_key_genes)

  # 过滤低质量基因（保护目标基因）
  if (length(bad_genes) > 0) {
    exp <- exp[!rownames(exp) %in% bad_genes, , drop = FALSE]
    log_msg("info", sprintf("过滤了%d个低质量基因，保留了%d个受保护的目标基因", length(bad_genes), length(protected_from_filter)))
  }

  list(
    exp = exp,
    raw_exp_gene_count = raw_exp_gene_count,
    filtered_exp_gene_count = nrow(exp),
    bad_genes = bad_genes,
    protected_from_filter = protected_from_filter,
    present_key_genes = present_key_genes
  )
}

# ---- GMT匹配率检查 ----

#' 检查GMT文件与表达矩阵的匹配率
#' @param gmt_file GMT文件路径
#' @param exp_genes 表达矩阵基因名
#' @param min_match_rate 最低匹配率阈值（百分比）
#' @return 匹配结果列表
check_gmt_match_rate <- function(gmt_file, exp_genes, min_match_rate = 30) {

  if (!file.exists(gmt_file)) {
    stop(sprintf("GMT文件不存在: %s", gmt_file))
  }

  # 读取GMT文件
  gmt_data <- clusterProfiler::read.gmt(gmt_file)
  all_gmt_genes <- unique(gmt_data$gene)

  # 计算匹配率 - 应该是相对于输入基因的匹配率
  matched_genes <- intersect(all_gmt_genes, exp_genes)
  match_rate <- length(matched_genes) / length(exp_genes) * 100

  log_msg("info", sprintf("GMT基因总数: %d, 匹配基因: %d, 匹配率: %.2f%%",
                          length(all_gmt_genes), length(matched_genes), match_rate))

  if (match_rate < min_match_rate) {
    log_msg("warn", sprintf("GMT匹配率低于阈值(%.1f%%)，可能影响GSEA结果", min_match_rate))
  }

  list(
    gmt_genes = all_gmt_genes,
    matched_genes = matched_genes,
    match_rate = match_rate,
    gmt_data = gmt_data
  )
}

# ---- 综合质控报告 ----

#' 生成综合质控报告
#' @param qc_results 质控结果列表
#' @param output_file 输出文件路径
generate_qc_report <- function(qc_results, output_file) {

  ensure_dir(dirname(output_file))

  # 构建质控报告数据框
  qc_items <- c()
  qc_values <- c()

  # GO/KEGG质控
  if (!is.null(qc_results$symbol_qc)) {
    sqc <- qc_results$symbol_qc
    qc_items <- c(qc_items,
                  "GO/KEGG 总基因数",
                  "GO/KEGG 匹配基因数",
                  "GO/KEGG 匹配率")
    qc_values <- c(qc_values,
                   sqc$total_genes,
                   length(sqc$matched_genes),
                   sprintf("%.2f%%", length(sqc$matched_genes) / sqc$total_genes * 100))
  }

  # GSEA表达矩阵质控
  if (!is.null(qc_results$exp_qc)) {
    eqc <- qc_results$exp_qc
    qc_items <- c(qc_items,
                  "GSEA 原始表达矩阵基因数",
                  "GSEA 过滤后基因数",
                  "GSEA 过滤的基因数")
    qc_values <- c(qc_values,
                   eqc$raw_exp_gene_count,
                   eqc$filtered_exp_gene_count,
                   length(eqc$bad_genes))
  }

  # GMT匹配质控
  if (!is.null(qc_results$gmt_qc)) {
    gqc <- qc_results$gmt_qc
    qc_items <- c(qc_items,
                  "GMT 基因总数",
                  "GMT 匹配基因数",
                  "GMT 匹配率")
    qc_values <- c(qc_values,
                   length(gqc$gmt_genes),
                   length(gqc$matched_genes),
                   sprintf("%.2f%%", gqc$match_rate))
  }

  # 添加分析日期
  qc_items <- c(qc_items, "分析日期", "物种")
  qc_values <- c(qc_values, as.character(Sys.Date()), qc_results$species %||% "未知")

  qc_df <- data.frame(
    Item = qc_items,
    Value = qc_values,
    stringsAsFactors = FALSE
  )

  write_csv_utf8(qc_df, output_file)
  log_msg("info", sprintf("质控报告已保存: %s", output_file))

  invisible(qc_df)
}

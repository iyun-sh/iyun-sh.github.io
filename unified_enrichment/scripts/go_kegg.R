# ==============================================================================
# unified_enrichment/scripts/go_kegg.R
# GO/KEGG过表达分析模块
# 提取自 enrichment/GO_KEGG_enrichment/code/enrichment_analysis.R
# ==============================================================================

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(enrichplot)
  library(ggplot2)
  library(viridis)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(httr)
})

# 设置KEGG API超时（默认60秒不足，KEGG服务器响应较慢，设为999秒应对极端情况）
options(timeout = 999)
httr::set_config(httr::timeout(999))

source("scripts/utils.R")

# ---- GO富集分析 ----

#' 运行GO富集分析
#' @param gene_list 基因列表（Symbol）
#' @param OrgDb 物种注释数据库
#' @param ont GO本体类型 (BP/MF/CC/ALL)
#' @param pvalue_cutoff p值阈值
#' @param qvalue_cutoff q值阈值
#' @param pAdjustMethod p值调整方法
#' @return GO富集结果对象
run_go_analysis <- function(gene_list, OrgDb,
                            ont = "ALL",
                            pvalue_cutoff = 0.05,
                            qvalue_cutoff = 0.2,
                            pAdjustMethod = "BH") {

  log_msg("info", sprintf("开始GO富集分析，基因数: %d", length(gene_list)))
  flush.console()

  # 基因ID转换
  log_msg("info", "开始bitr转换...")
  gene_transform <- clusterProfiler::bitr(
    gene_list,
    fromType = "SYMBOL",
    toType = c("ENTREZID", "ENSEMBL"),
    OrgDb = OrgDb,
    drop = FALSE
  )

  log_msg("info", sprintf("成功映射基因: %d / %d", sum(!is.na(gene_transform$ENTREZID)), length(gene_list)))
  flush.console()

  # 移除NA，只保留有效ENTREZID
  valid_genes <- gene_transform$ENTREZID[!is.na(gene_transform$ENTREZID)]
  log_msg("info", sprintf("有效基因数: %d", length(valid_genes)))
  flush.console()

  # GO富集分析 - 直接调用，不用tryCatch以便看到错误
  log_msg("info", "开始运行enrichGO...")
  flush.console()

  ego <- enrichGO(
    gene = valid_genes,
    OrgDb = OrgDb,
    keyType = "ENTREZID",
    ont = ont,
    pAdjustMethod = pAdjustMethod,
    pvalueCutoff = pvalue_cutoff,
    qvalueCutoff = qvalue_cutoff,
    readable = TRUE
  )

  log_msg("info", "enrichGO执行完成")
  flush.console()

  # 提取结果
  go_results <- ego@result
  log_msg("info", sprintf("GO富集分析完成，显著条目: %d", nrow(go_results)))

  list(
    result = go_results,
    ego = ego,
    gene_transform = gene_transform
  )
}

# ---- KEGG富集分析 ----

#' 运行KEGG富集分析
#' @param gene_list 基因列表（Symbol）
#' @param OrgDb 物种注释数据库
#' @param species 物种代码 (hsa/mmu/rno)
#' @param pvalue_cutoff p值阈值
#' @param qvalue_cutoff q值阈值
#' @param pAdjustMethod p值调整方法
#' @return KEGG富集结果对象
run_kegg_analysis <- function(gene_list, OrgDb, species = "hsa",
                              pvalue_cutoff = 0.05,
                              qvalue_cutoff = 0.2,
                              pAdjustMethod = "BH") {

  log_msg("info", sprintf("开始KEGG富集分析，基因数: %d", length(gene_list)))

  # 基因ID转换
  gene_transform <- clusterProfiler::bitr(
    gene_list,
    fromType = "SYMBOL",
    toType = c("ENTREZID", "ENSEMBL"),
    OrgDb = OrgDb,
    drop = FALSE
  )

  # KEGG富集分析
  kk <- tryCatch({
    enrichKEGG(
      gene = gene_transform$ENTREZID,
      keyType = "kegg",
      organism = species,
      pAdjustMethod = pAdjustMethod,
      pvalueCutoff = pvalue_cutoff,
      qvalueCutoff = qvalue_cutoff
    )
  }, error = function(e) {
    log_msg("error", sprintf("KEGG富集分析失败: %s", e$message))
    NULL
  })

  if (is.null(kk)) {
    log_msg("warn", "KEGG富集分析未返回结果")
    return(list(result = NULL, gene_transform = gene_transform))
  }

  # 转换为可读的基因名
  kk <- setReadable(kk, OrgDb = OrgDb, keyType = "ENTREZID")

  # 提取结果
  kegg_results <- kk@result
  log_msg("info", sprintf("KEGG富集分析完成，显著通路: %d", nrow(kegg_results)))

  list(
    result = kegg_results,
    kk = kk,
    gene_transform = gene_transform
  )
}

# ---- 可视化 ----

#' 绘制GO富集分析气泡图
#' @param go_results GO结果数据框
#' @param value_col 用于着色的值列 (pvalue/p.adjust)
#' @param top_n 每个本体显示的条目数
#' @param color_scheme 配色方案
#' @param output_dir 输出目录
#' @param plot_width 图宽度
#' @param plot_height 图高度
plot_go_results <- function(go_results, value_col = "p.adjust",
                            top_n = 5, color_scheme = "viridis",
                            output_dir = "results/go",
                            plot_width = 10, plot_height = 8) {

  if (is.null(go_results) || nrow(go_results) == 0) {
    log_msg("warn", "没有GO富集结果，跳过可视化")
    return(invisible(NULL))
  }

  ensure_dir(output_dir)

  # 计算RichFactor
  eGo <- go_results %>%
    separate(col = GeneRatio, into = c("GR1", "GR2"), sep = "/", convert = TRUE) %>%
    mutate(GeneRatio = GR1 / GR2) %>%
    separate(col = BgRatio, into = c("BR1", "BR2"), sep = "/", convert = TRUE) %>%
    mutate(BgRatio = BR1 / BR2, RichFactor = GeneRatio / BgRatio)

  # 每个本体取top_n
  eGo_p <- eGo %>%
    group_by(ONTOLOGY) %>%
    arrange(!!sym(value_col)) %>%
    slice_head(n = top_n) %>%
    ungroup()

  # 绘图主题 - Arial字体
  theme_set <- theme(
    title = element_text(family = "Arial", face = "bold", size = 12),
    legend.title = element_text(family = "Arial", face = "bold", size = 11),
    axis.title.x = element_text(family = "Arial", face = "bold", size = 12),
    axis.title.y = element_text(family = "Arial", face = "bold", size = 12),
    strip.text = element_text(family = "Arial", face = "bold", size = 11),
    axis.text.x = element_text(family = "Arial", size = 9, color = "black"),
    axis.text.y = element_text(family = "Arial", size = 9, color = "black"),
    legend.text = element_text(family = "Arial", size = 9)
  )

  # 气泡图 - 英文标签
  p_go <- ggplot(eGo_p, aes(RichFactor, reorder(str_wrap(Description, 50), RichFactor))) +
    geom_point(aes(size = Count, color = !!sym(value_col))) +
    scale_color_viridis_c(option = color_scheme, direction = -1,
                          labels = scales::scientific_format()) +
    scale_size_continuous(range = c(3, 6), trans = "sqrt") +
    labs(color = "p.adjust", size = "Count", x = "RichFactor", y = "Term",
         title = "GO Enrichment Dotplot") +
    theme_bw() +
    facet_grid(ONTOLOGY ~ ., scales = "free", space = 'free') +
    theme_set +
    theme(plot.title = element_text(family = "Arial", face = "bold", hjust = 0.5))

  # 保存
  save_plot("go_enrichment_plot.pdf", p_go,
            outdir = output_dir, width = plot_width, height = plot_height, both = TRUE)

  # 生成图对应的统计摘要CSV（供报告自动解读）
  tryCatch({
    go_stats <- data.frame(
      metric = c("total_terms", "min_pvalue", "max_pvalue",
                 "top1_term", "top1_pvalue", "top1_count",
                 "enriched_BP", "enriched_CC", "enriched_MF"),
      value = c(
        nrow(go_plot_data),
        min(go_plot_data$p.adjust, na.rm = TRUE),
        max(go_plot_data$p.adjust, na.rm = TRUE),
        as.character(go_plot_data$Description[which.min(go_plot_data$p.adjust)]),
        min(go_plot_data$p.adjust, na.rm = TRUE),
        go_plot_data$Count[which.min(go_plot_data$p.adjust)],
        sum(go_plot_data$ONTOLOGY == "BP", na.rm = TRUE),
        sum(go_plot_data$ONTOLOGY == "CC", na.rm = TRUE),
        sum(go_plot_data$ONTOLOGY == "MF", na.rm = TRUE)
      ),
      stringsAsFactors = FALSE
    )
    fwrite(go_stats, file.path(output_dir, "go", "go_figure_stats.csv"))
  }, error = function(e) log_msg("warn", paste("GO figure stats generation failed:", e$message)))

  invisible(p_go)
}

#' 绘制KEGG富集分析气泡图
#' @param kegg_results KEGG结果数据框
#' @param value_col 用于着色的值列 (pvalue/p.adjust)
#' @param top_n 显示的条目数
#' @param color_scheme 配色方案
#' @param output_dir 输出目录
#' @param plot_width 图宽度
#' @param plot_height 图高度
plot_kegg_results <- function(kegg_results, value_col = "p.adjust",
                              top_n = 5, color_scheme = "viridis",
                              output_dir = "results/kegg",
                              plot_width = 10, plot_height = 8) {

  if (is.null(kegg_results) || nrow(kegg_results) == 0) {
    log_msg("warn", "没有KEGG富集结果，跳过可视化")
    return(invisible(NULL))
  }

  ensure_dir(output_dir)

  # 计算RichFactor
  ekk <- kegg_results %>%
    separate(col = GeneRatio, into = c("GR1", "GR2"), sep = "/", convert = TRUE) %>%
    mutate(GeneRatio = GR1 / GR2) %>%
    separate(col = BgRatio, into = c("BR1", "BR2"), sep = "/", convert = TRUE) %>%
    mutate(BgRatio = BR1 / BR2, RichFactor = GeneRatio / BgRatio)

  # 取top_n
  ekk_p <- ekk %>%
    arrange(!!sym(value_col)) %>%
    slice_head(n = top_n)

  # 绘图主题 - Arial字体
  theme_set <- theme(
    title = element_text(family = "Arial", face = "bold", size = 12),
    legend.title = element_text(family = "Arial", face = "bold", size = 11),
    axis.title.x = element_text(family = "Arial", face = "bold", size = 12),
    axis.title.y = element_text(family = "Arial", face = "bold", size = 12),
    axis.text.x = element_text(family = "Arial", size = 9, color = "black"),
    axis.text.y = element_text(family = "Arial", size = 9, color = "black"),
    legend.text = element_text(family = "Arial", size = 9)
  )

  # 气泡图 - 英文标签
  p_kegg <- ggplot(ekk_p, aes(RichFactor, reorder(str_wrap(Description, 50), RichFactor))) +
    geom_point(aes(size = Count, color = !!sym(value_col))) +
    scale_color_viridis_c(option = color_scheme, direction = -1,
                          labels = scales::scientific_format()) +
    labs(color = "p.adjust", size = "Count", x = "RichFactor", y = "Pathway",
         title = "KEGG Pathway Enrichment Dotplot") +
    theme_bw() +
    theme_set +
    theme(plot.title = element_text(family = "Arial", face = "bold", hjust = 0.5))

  # 保存
  save_plot("kegg_enrichment_plot.pdf", p_kegg,
            outdir = output_dir, width = plot_width, height = plot_height, both = TRUE)

  # 生成图对应的统计摘要CSV（供报告自动解读）
  tryCatch({
    kegg_stats <- data.frame(
      metric = c("total_pathways", "min_pvalue", "max_pvalue",
                 "top1_pathway", "top1_pvalue", "top1_count",
                 "max_richfactor", "max_richfactor_pathway"),
      value = c(
        nrow(ekk_p),
        min(ekk_p$p.adjust, na.rm = TRUE),
        max(ekk_p$p.adjust, na.rm = TRUE),
        as.character(ekk_p$Description[which.min(ekk_p$p.adjust)]),
        min(ekk_p$p.adjust, na.rm = TRUE),
        ekk_p$Count[which.min(ekk_p$p.adjust)],
        max(ekk_p$RichFactor, na.rm = TRUE),
        as.character(ekk_p$Description[which.max(ekk_p$RichFactor)])
      ),
      stringsAsFactors = FALSE
    )
    fwrite(kegg_stats, file.path(output_dir, "kegg", "kegg_figure_stats.csv"))
  }, error = function(e) log_msg("warn", paste("KEGG figure stats generation failed:", e$message)))

  invisible(p_kegg)
}

# ---- 保存结果 ----

#' 保存GO/KEGG分析结果
#' @param go_result GO结果
#' @param kegg_result KEGG结果
#' @param output_dir 输出目录
save_go_kegg_results <- function(go_result, kegg_result, output_dir = "results") {

  ensure_dir(output_dir)

  # 保存GO结果
  if (!is.null(go_result) && nrow(go_result) > 0) {
    go_file <- file.path(output_dir, "go", "go_enrichment_results.csv")
    ensure_dir(dirname(go_file))
    fwrite(go_result, go_file)
    log_msg("info", sprintf("GO结果已保存: %s", go_file))
  } else {
    log_msg("warn", "没有GO富集结果")
  }

  # 保存KEGG结果
  if (!is.null(kegg_result) && nrow(kegg_result) > 0) {
    kegg_file <- file.path(output_dir, "kegg", "kegg_enrichment_results.csv")
    ensure_dir(dirname(kegg_file))
    fwrite(kegg_result, kegg_file)
    log_msg("info", sprintf("KEGG结果已保存: %s", kegg_file))
  } else {
    log_msg("warn", "没有KEGG富集结果")
  }
}

# ---- 完整GO/KEGG分析流程 ----

#' 运行完整的GO/KEGG分析流程
#' @param cfg 配置列表
#' @param paths 路径列表
#' @return 分析结果列表
run_go_kegg_pipeline <- function(cfg, paths) {

  log_msg("info", "========== 开始GO/KEGG富集分析 ==========")

  # 读取基因列表
  gene_list <- read_gene_list(paths$gene_list)
  log_msg("info", sprintf("读取基因列表: %d 个基因", length(gene_list)))

  # 运行质控
  qc_result <- symbol_qc(paths$gene_list, paths$qc_dir, get_cfg(cfg, "project.species", "hsa"))

  # 加载物种数据库
  OrgDb <- get_org_db(get_cfg(cfg, "project.species", "hsa"))

  # 设置p值调整方法
  pAdjustMethod <- get_cfg(cfg, "analysis.go_kegg.padj_method", "BH")
  if (pAdjustMethod == "none") pAdjustMethod <- "none"
  value_col <- if (pAdjustMethod == "none") "pvalue" else "p.adjust"

  log_msg("info", "准备运行GO分析...")
  flush.console()

  # GO分析
  log_msg("info", "调用run_go_analysis...")
  go_res <- run_go_analysis(
    gene_list = gene_list,
    OrgDb = OrgDb,
    ont = get_cfg(cfg, "analysis.go_kegg.go_ont", "ALL"),
    pvalue_cutoff = get_cfg(cfg, "analysis.go_kegg.pvalue_cutoff", 0.05),
    qvalue_cutoff = get_cfg(cfg, "analysis.go_kegg.qvalue_cutoff", 0.2),
    pAdjustMethod = pAdjustMethod
  )
  log_msg("info", "GO分析返回")
  flush.console()

  # KEGG分析
  kegg_res <- run_kegg_analysis(
    gene_list = gene_list,
    OrgDb = OrgDb,
    species = get_cfg(cfg, "project.species", "hsa"),
    pvalue_cutoff = get_cfg(cfg, "analysis.go_kegg.pvalue_cutoff", 0.05),
    qvalue_cutoff = get_cfg(cfg, "analysis.go_kegg.qvalue_cutoff", 0.2),
    pAdjustMethod = pAdjustMethod
  )

  # 保存结果
  save_go_kegg_results(go_res$result, kegg_res$result, paths$results_dir)

  # 生成图表
  plot_go_results(
    go_results = go_res$result,
    value_col = value_col,
    top_n = get_cfg(cfg, "plot.top_n", 5),
    color_scheme = get_cfg(cfg, "plot.color_scheme", "viridis"),
    output_dir = paths$go_dir,
    plot_width = get_cfg(cfg, "plot.width", 10),
    plot_height = get_cfg(cfg, "plot.height", 8)
  )

  plot_kegg_results(
    kegg_results = kegg_res$result,
    value_col = value_col,
    top_n = get_cfg(cfg, "plot.top_n", 5),
    color_scheme = get_cfg(cfg, "plot.color_scheme", "viridis"),
    output_dir = paths$kegg_dir,
    plot_width = get_cfg(cfg, "plot.width", 10),
    plot_height = get_cfg(cfg, "plot.height", 8)
  )

  log_msg("info", "========== GO/KEGG富集分析完成 ==========")

  list(
    go = go_res,
    kegg = kegg_res,
    qc = qc_result
  )
}

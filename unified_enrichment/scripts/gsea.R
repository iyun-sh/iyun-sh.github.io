# ==============================================================================
# unified_enrichment/scripts/gsea.R
# GSEA基因集富集分析模块
# 支持两种方法：
#   1. prerank: 使用预排序基因列表（如log2FC）
#   2. correlation: 基于目标基因与全基因组相关性排序
# ==============================================================================

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(enrichplot)
  library(ggplot2)
  library(ggsci)
  library(viridisLite)
  library(dplyr)
  library(data.table)
})

source("scripts/utils.R")

# ---- Pre-rank GSEA 专用函数 ----

#' 读取预排序基因列表
#' @param file_path 文件路径
#' @param metric_col 排序指标列名（如log2FC）
#' @param gene_col 基因名列名
#' @return 排序后的命名向量
read_ranked_gene_list <- function(file_path, metric_col = "log2FC", gene_col = "gene") {
  df <- fread(file_path, data.table = FALSE)

  # 自动检测列名
  possible_gene_cols <- c(
    "gene", "Gene", "symbol", "Symbol", "SYMBOL", "GeneSymbol", "gene_symbol",
    "ENSEMBL", "ensembl", "ensembl_id", "Ensembl", "gene_id"
  )
  possible_metric_cols <- c("log2FC", "logFC", "log2_fold_change", "FC", "fold_change",
                            "t_stat", "tstat", "statistic", "score", "rankingMetric")

  # 匹配基因列
  gene_col_found <- NULL
  for (col in possible_gene_cols) {
    if (col %in% colnames(df)) {
      gene_col_found <- col
      break
    }
  }
  if (is.null(gene_col_found)) {
    # 默认使用第一列
    gene_col_found <- colnames(df)[1]
    log_msg("warn", sprintf("未找到基因列，使用第一列: %s", gene_col_found))
  }

  # 匹配指标列
  metric_col_found <- NULL
  if (metric_col %in% colnames(df)) {
    metric_col_found <- metric_col
  } else {
    for (col in possible_metric_cols) {
      if (col %in% colnames(df)) {
        metric_col_found <- col
        log_msg("info", sprintf("使用指标列: %s", metric_col_found))
        break
      }
    }
  }
  if (is.null(metric_col_found)) {
    # 默认使用最后一列
    metric_col_found <- colnames(df)[ncol(df)]
    log_msg("warn", sprintf("未找到指标列，使用最后一列: %s", metric_col_found))
  }

  # 构建排序向量
  gene_list <- df[[metric_col_found]]
  names(gene_list) <- df[[gene_col_found]]

  # 移除NA和无限值
  gene_list <- gene_list[!is.na(gene_list) & is.finite(gene_list)]

  # 降序排列
  gene_list <- sort(gene_list, decreasing = TRUE)

  # 处理重复基因（保留最高值）
  if (any(duplicated(names(gene_list)))) {
    dup_genes <- names(gene_list)[duplicated(names(gene_list))]
    log_msg("warn", sprintf("发现%d个重复基因，保留最高值", length(dup_genes)))
    gene_list <- gene_list[!duplicated(names(gene_list))]
  }

  log_msg("info", sprintf("读取排序基因列表: %d 个基因", length(gene_list)))
  log_msg("info", sprintf("指标范围: [%.4f, %.4f]", min(gene_list), max(gene_list)))

  gene_list
}

#' 规范化基因 ID（匹配 GMT / 排序表）
#' @param ids 基因 ID 向量
#' @param id_type infer_expr_id 结果: symbols / ensembl / unknown
normalize_gene_ids_for_gsea <- function(ids, id_type = "unknown") {
  ids <- trimws(as.character(ids))
  ids[is.na(ids)] <- ""
  ids <- ids[ids != ""]
  if (length(ids) == 0) return(character(0))
  type <- tolower(id_type %||% "unknown")
  if (type == "ensembl") {
    sub("\\..*$", "", toupper(ids))
  } else if (type == "symbols") {
    toupper(ids)
  } else {
    ids
  }
}

#' 将排序基因向量与 TERM2GENE 对齐到同一 ID 体系（SYMBOL 或 ENSEMBL）
#' @param gene_rank 命名向量（name = 基因 ID）
#' @param term2gene 两列 term, gene
#' @param species_hint hsa / mmu / rno，NULL 则自动推断
#' @return list(gene_rank, term2gene, species_used, id_converted)
align_gene_rank_to_gmt <- function(gene_rank, term2gene, species_hint = NULL) {
  if (length(gene_rank) == 0 || is.null(term2gene) || nrow(term2gene) == 0) {
    return(list(
      gene_rank = gene_rank,
      term2gene = term2gene,
      species_used = species_hint %||% "hsa",
      id_converted = FALSE
    ))
  }

  term2gene <- as.data.frame(term2gene, stringsAsFactors = FALSE)
  colnames(term2gene)[1:2] <- c("term", "gene")

  gene_ids_raw <- names(gene_rank)
  gmt_ids_raw <- unique(as.character(term2gene$gene))

  gene_info <- infer_expr_id(gene_ids_raw)
  gmt_info <- infer_expr_id(gmt_ids_raw)
  gene_id_type <- tolower(gene_info$id_type %||% "unknown")
  gmt_id_type <- tolower(gmt_info$id_type %||% "unknown")

  sp <- species_hint
  if (is.null(sp) || is.na(sp) || sp == "" || sp == "unknown") {
    sp <- if (gene_info$species != "unknown") gene_info$species else gmt_info$species
  }
  if (is.null(sp) || is.na(sp) || sp == "" || sp == "unknown") {
    sp <- "hsa"
  }

  names(gene_rank) <- normalize_gene_ids_for_gsea(gene_ids_raw, gene_id_type)
  term2gene$gene <- normalize_gene_ids_for_gsea(term2gene$gene, gmt_id_type)
  term2gene <- term2gene[!is.na(term2gene$gene) & term2gene$gene != "", , drop = FALSE]
  term2gene <- unique(term2gene[, c("term", "gene"), drop = FALSE])

  gene_rank <- gene_rank[!is.na(names(gene_rank)) & names(gene_rank) != ""]
  if (any(duplicated(names(gene_rank)))) {
    gene_rank <- gene_rank[!duplicated(names(gene_rank))]
  }

  id_converted <- FALSE
  need_convert <- gene_id_type %in% c("symbols", "ensembl") &&
    gmt_id_type %in% c("symbols", "ensembl") &&
    gene_id_type != gmt_id_type

  # MSigDB 等 GMT 多为大写 SYMBOL；unknown 类型也统一大写减少假阴性
  names(gene_rank) <- toupper(names(gene_rank))
  term2gene$gene <- toupper(as.character(term2gene$gene))
  term2gene <- term2gene[!is.na(term2gene$gene) & term2gene$gene != "", , drop = FALSE]
  term2gene <- unique(term2gene[, c("term", "gene"), drop = FALSE])
  if (any(duplicated(names(gene_rank)))) {
    gene_rank <- gene_rank[!duplicated(names(gene_rank))]
  }

  if (need_convert) {
    from_type <- if (gene_id_type == "ensembl") "ENSEMBL" else "SYMBOL"
    to_type <- if (gmt_id_type == "ensembl") "ENSEMBL" else "SYMBOL"
    map_keys <- unique(names(gene_rank))
    map_df <- tryCatch({
      OrgDb <- get_org_db(sp)
      suppressMessages(
        AnnotationDbi::select(OrgDb, keys = map_keys, keytype = from_type, columns = c(to_type, from_type))
      )
    }, error = function(e) {
      log_msg("warn", sprintf("GSEA ID 转换失败 (%s->%s): %s", from_type, to_type, e$message))
      NULL
    })

    if (!is.null(map_df) && nrow(map_df) > 0 && all(c(from_type, to_type) %in% colnames(map_df))) {
      map_df <- map_df[!is.na(map_df[[from_type]]) & !is.na(map_df[[to_type]]), , drop = FALSE]
      map_df[[from_type]] <- normalize_gene_ids_for_gsea(map_df[[from_type]], gene_id_type)
      map_df[[to_type]] <- normalize_gene_ids_for_gsea(map_df[[to_type]], gmt_id_type)
      map_df <- map_df[map_df[[from_type]] != "" & map_df[[to_type]] != "", , drop = FALSE]
      map_df <- map_df[!duplicated(map_df[[from_type]]), , drop = FALSE]

      if (nrow(map_df) > 0) {
        idx <- match(names(gene_rank), map_df[[from_type]])
        mapped <- map_df[[to_type]][idx]
        keep <- !is.na(mapped) & mapped != ""
        gene_rank <- gene_rank[keep]
        names(gene_rank) <- mapped[keep]
        gene_rank <- gene_rank[!duplicated(names(gene_rank))]
        id_converted <- TRUE
        log_msg("info", sprintf(
          "GSEA 已对齐 ID：%s -> %s（物种=%s），保留排序基因 %d 个",
          from_type, to_type, sp, length(gene_rank)
        ))
      }
    }
  }

  list(
    gene_rank = gene_rank,
    term2gene = term2gene,
    species_used = sp,
    id_converted = id_converted
  )
}

#' 运行Pre-rank GSEA分析
#' @param gene_rank 排序后的基因列表
#' @param term2gene GMT数据
#' @param gsea_cfg GSEA配置参数
#' @param all_gmt_genes GMT中所有基因（对齐前，仅兼容旧参数）
#' @param species_hint 物种代码 hsa/mmu/rno
#' @return GSEA结果和质控信息；含对齐后的 gene_rank / term2gene
run_prerank_gsea <- function(gene_rank, term2gene, gsea_cfg, all_gmt_genes, species_hint = NULL) {

  aligned <- align_gene_rank_to_gmt(gene_rank, term2gene, species_hint = species_hint)
  gene_rank <- aligned$gene_rank
  term2gene <- aligned$term2gene
  all_gmt_genes <- unique(term2gene$gene)

  if (length(gene_rank) < 5L) {
    log_msg("error", sprintf("ID 对齐后排序基因过少（%d），无法运行 GSEA", length(gene_rank)))
    return(NULL)
  }

  # 质控统计（对齐后）
  valid_genes_count <- length(gene_rank)
  na_genes_in_rank <- 0  # 已在读取时处理
  gmt_genes_in_rank <- sum(names(gene_rank) %in% all_gmt_genes)
  gmt_match_rate <- if (length(gene_rank) > 0) {
    (gmt_genes_in_rank / length(gene_rank)) * 100
  } else {
    0
  }

  log_msg("info", sprintf(
    "Pre-rank GSEA质控(对齐后): 输入基因=%d, GMT唯一基因=%d, 排序表与GMT交集=%d (占输入%.1f%%)",
    valid_genes_count, length(all_gmt_genes), gmt_genes_in_rank, gmt_match_rate
  ))

  # 运行GSEA - 使用 clusterProfiler::GSEA
  gsea <- tryCatch({
    # 检查是否所有值都是正数
    all_positive <- all(gene_rank > 0, na.rm = TRUE)
    if (all_positive) {
      log_msg("warn", "检测到所有基因值都为正，GSEA可能无显著结果（建议使用包含正负值的排序列表）")
    }

    # 使用 clusterProfiler::GSEA
    result <- GSEA(
      geneList = gene_rank,
      exponent = gsea_cfg$exponent,
      minGSSize = gsea_cfg$minGSSize,
      maxGSSize = gsea_cfg$maxGSSize,
      pvalueCutoff = gsea_cfg$pvalueCutoff,
      pAdjustMethod = gsea_cfg$padj_method,
      TERM2GENE = term2gene,
      eps = gsea_cfg$eps,
      verbose = FALSE
    )

    result
  }, error = function(e) {
    log_msg("error", sprintf("Pre-rank GSEA失败: %s", e$message))
    return(NULL)
  })

  if (is.null(gsea)) {
    return(NULL)
  }

  metric_rng <- if (length(gene_rank) > 0) {
    paste0("[", round(min(gene_rank), 4), ", ", round(max(gene_rank), 4), "]")
  } else {
    "NA"
  }

  # 质控信息
  prerank_qc <- data.frame(
    Method = "Pre-rank GSEA",
    Input_Genes = valid_genes_count,
    GMT_Genes_Total = length(all_gmt_genes),
    GMT_Genes_Matched = gmt_genes_in_rank,
    GMT_Match_Rate_Percent = round(gmt_match_rate, 2),
    Metric_Range = metric_rng,
    Significant_Pathways = nrow(gsea@result),
    stringsAsFactors = FALSE
  )

  list(
    gsea = gsea,
    qc = prerank_qc,
    gene_rank = gene_rank,
    term2gene = term2gene,
    species_used = aligned$species_used,
    id_converted = aligned$id_converted
  )
}

# ---- 相关性计算 ----

#' 计算目标基因与全基因组的相关性
#' @param exp 表达矩阵（行=基因，列=样本）
#' @param key_genes 目标基因列表
#' @param method 相关性方法 (spearman/pearson)
#' @return 排序后的相关性列表
calculate_correlation <- function(exp, key_genes, method = "spearman") {
  cor_results <- list()

  for (gene in key_genes) {
    if (!gene %in% rownames(exp)) {
      log_msg("warn", sprintf("跳过不存在于矩阵的目标基因: %s", gene))
      next
    }

    # 行=基因，列=样本，所以取行
    cor_values <- cor(exp[gene, ], t(exp), method = method, use = "pairwise.complete.obs")
    cor_values <- as.numeric(cor_values[1, ])
    names(cor_values) <- rownames(exp)
    cor_values <- cor_values[names(cor_values) != gene]

    # 处理重复值（添加微小噪声）
    if (length(unique(cor_values)) < length(cor_values)) {
      cor_values <- cor_values + rnorm(length(cor_values), mean = 0, sd = 1e-8)
    }

    cor_results[[gene]] <- sort(cor_values, decreasing = TRUE)
  }

  cor_results
}

# ---- GSEA分析 ----

#' 对单个基因运行GSEA分析
#' @param gene 目标基因名
#' @param gene_rank 排序后的基因列表
#' @param term2gene GMT数据（term-gene映射）
#' @param gsea_cfg GSEA配置参数
#' @param all_gmt_genes GMT中所有基因
#' @return GSEA结果和质控信息
run_gsea_for_gene <- function(gene, gene_rank, term2gene, gsea_cfg, all_gmt_genes, species_hint = NULL) {

  aligned <- align_gene_rank_to_gmt(gene_rank, term2gene, species_hint = species_hint)
  gene_rank <- aligned$gene_rank
  term2gene <- aligned$term2gene
  all_gmt_genes <- unique(term2gene$gene)

  # 质控统计
  valid_genes_count <- length(gene_rank)
  na_genes_in_rank <- sum(is.na(gene_rank))
  gmt_genes_in_rank <- if (length(gene_rank) > 0) sum(names(gene_rank) %in% all_gmt_genes) else 0L
  gmt_match_rate <- if (length(gene_rank) > 0) (gmt_genes_in_rank / length(gene_rank)) * 100 else 0

  # 运行GSEA
  gsea <- tryCatch({
    if (length(gene_rank) < 5L) {
      stop("对齐后排序基因过少")
    }
    GSEA(
      geneList = gene_rank,
      exponent = gsea_cfg$exponent,
      minGSSize = gsea_cfg$minGSSize,
      maxGSSize = gsea_cfg$maxGSSize,
      pvalueCutoff = gsea_cfg$pvalueCutoff,
      pAdjustMethod = gsea_cfg$padj_method,
      TERM2GENE = term2gene,
      eps = gsea_cfg$eps
    )
  }, error = function(e) {
    log_msg("error", sprintf("%s GSEA失败: %s", gene, e$message))
    return(NULL)
  })

  if (is.null(gsea)) {
    return(NULL)
  }

  # 安全访问 result
  sig_count <- tryCatch({
    nrow(gsea@result)
  }, error = function(e) {
    log_msg("error", sprintf("访问gsea结果失败: %s", e$message))
    0
  })

  # 基因级质控信息
  gene_qc <- data.frame(
    Target_Gene = gene,
    Valid_Genes_in_Ranking = valid_genes_count,
    NA_Genes = na_genes_in_rank,
    GMT_Genes_in_Ranking = gmt_genes_in_rank,
    GMT_Match_Rate_Percent = round(gmt_match_rate, 2),
    Top10_Correlated_Genes = paste(names(gene_rank[1:min(10, length(gene_rank))]), collapse = "/"),
    Correlation_Range = paste0("[", round(min(gene_rank), 3), ", ", round(max(gene_rank), 3), "]"),
    Significant_Pathways = sig_count,
    stringsAsFactors = FALSE
  )

  list(gsea = gsea, gene_qc = gene_qc, gene_rank_used = gene_rank)
}

# ---- 保存GSEA结果 ----

#' 保存单个基因的GSEA结果
#' @param gene 基因名
#' @param gsea_obj GSEA结果对象
#' @param paths 路径列表
save_gene_outputs <- function(gene, gsea_obj, paths) {
  result <- gsea_obj@result

  # 保存完整结果
  result_file <- file.path(paths$gsea_dir, paste0(gene, "_gsea_analysis_results.csv"))
  write_csv_utf8(result, result_file)

  if (nrow(result) > 0) {
    # 计算匹配率
    result$matched_gene_count <- sapply(result$core_enrichment, function(x) {
      if (is.na(x) || x == "") return(0)
      length(strsplit(x, "/")[[1]])
    })
    result$match_rate_in_pathway <- (result$matched_gene_count / result$setSize) * 100

    # 保存详细报告
    detail <- result[, c("ID", "setSize", "matched_gene_count", "match_rate_in_pathway",
                         "NES", "pvalue", "p.adjust")]
    detail_file <- file.path(paths$gsea_dir, paste0(gene, "_GSEA_DetailedReport.csv"))
    write_csv_utf8(detail, detail_file)

    # 额外保存到tables目录
    table_csv <- file.path(paths$tables_dir, paste0(gene, "_GSEA_DetailedReport.csv"))
    write_csv_utf8(detail, table_csv)
  }
}

# ---- GSEA可视化 ----

#' 绘制单个基因的GSEA结果图
#' @param gene 基因名
#' @param gsea_obj GSEA结果对象
#' @param plot_cfg 绘图配置
#' @param top_n 显示的通路数
#' @param plot_dir 图形输出目录
plot_gene_results <- function(gene, gsea_obj, plot_cfg, top_n = 5, plot_dir) {

  if (is.null(gsea_obj) || nrow(gsea_obj@result) == 0) {
    log_msg("warn", sprintf("%s 无显著富集结果，跳过绘图", gene))
    return(invisible(NULL))
  }

  show_n <- min(plot_cfg$ridge_show_category %||% 10, nrow(gsea_obj@result))
  if (show_n < 1) {
    log_msg("warn", sprintf("%s 无可绘制通路，跳过绘图", gene))
    return(invisible(NULL))
  }

  ensure_dir(plot_dir)

  # 绘制山脊图 - Arial字体
  tryCatch({
    p_ridge <- ridgeplot(
      gsea_obj,
      showCategory = show_n,
      fill = "pvalue",
      core_enrichment = TRUE,
      label_format = 30,
      orderBy = "NES",
      decreasing = TRUE
    ) +
      theme_classic(base_size = 10) +
      theme(axis.text.y = element_text(family = "Arial", size = 9),
            axis.text.x = element_text(family = "Arial", size = 9),
            axis.title.x = element_text(family = "Arial", size = 11),
            axis.title.y = element_text(family = "Arial", size = 11),
            plot.title = element_text(family = "Arial", hjust = 0.5, size = 14, face = "bold"),
            legend.text = element_text(family = "Arial", size = 9),
            legend.title = element_text(family = "Arial", size = 10)) +
      labs(title = paste0(gene, " GSEA Ridge Plot"), x = "Rank", y = "Pathway") +
      scale_fill_gradientn(colors = viridisLite::magma(256, direction = -1, begin = 0.05, end = 0.95),
                           name = "p value")

    ridge_file <- file.path(plot_dir, paste0(gene, "_RidgePlot_SCI.pdf"))
    ggsave(ridge_file, p_ridge,
           width = plot_cfg$ridge_width %||% 10,
           height = plot_cfg$ridge_height %||% 6)
  }, error = function(e) {
    log_msg("warn", sprintf("%s Ridge plot failed: %s", gene, e$message))
  })

  # 选择top通路绘制GSEA曲线图
  top_pathways <- gsea_obj@result %>%
    arrange(desc(abs(NES))) %>%
    head(min(top_n, nrow(gsea_obj@result)))
  top_ids <- top_pathways$ID

  if (length(top_ids) == 0) {
    log_msg("warn", sprintf("%s 无可绘制富集曲线", gene))
    return(invisible(NULL))
  }

  # SCI配色
  sci_colors <- pal_npg(alpha = 0.9)(max(5, length(top_ids)))[seq_along(top_ids)]

  # 绘制组合GSEA图 - Arial字体
  p_all <- gseaplot2(
    gsea_obj,
    geneSetID = top_ids,
    color = sci_colors,
    pvalue_table = FALSE,
    subplots = 1:3,
    base_size = 9,
    rel_heights = c(1.5, 0.5, 1)
  )
  p_all[[1]] <- p_all[[1]] +
    theme(
      legend.position = "top",
      legend.direction = "vertical",
      text = element_text(family = "Arial"),
      axis.text = element_text(family = "Arial", size = 9),
      axis.title = element_text(family = "Arial", size = 11),
      plot.title = element_text(family = "Arial", face = "bold"),
      legend.text = element_text(family = "Arial", size = 9)
    )

  # 保存PDF - Arial字体
  all_pdf <- file.path(plot_dir, paste0(gene, "_GSEA_EnrichmentPlot_SCI.pdf"))
  if (capabilities("cairo")) {
    grDevices::cairo_pdf(all_pdf, width = plot_cfg$width %||% 12, height = plot_cfg$height %||% 8, family = "Arial")
  } else {
    grDevices::pdf(all_pdf, width = plot_cfg$width %||% 12, height = plot_cfg$height %||% 8, family = "Arial")
  }
  print(p_all)
  dev.off()

  # 保存PNG（用于报告嵌入）- Arial字体
  all_png <- file.path(plot_dir, paste0(gene, "_GSEA_EnrichmentPlot_SCI.png"))
  tryCatch({
    png(all_png, width = (plot_cfg$width %||% 12) * 150,
        height = (plot_cfg$height %||% 8) * 150, res = 150, family = "Arial")
    print(p_all)
    dev.off()
  }, error = function(e) {
    log_msg("warn", sprintf("%s PNG save failed: %s", gene, e$message))
  })

  # 绘制每通路单独的GSEA图 - Arial字体
  each_pdf <- file.path(plot_dir, paste0(gene, "_GSEA_EnrichmentPlot_SCI_PerPathway.pdf"))
  if (capabilities("cairo")) {
    grDevices::cairo_pdf(each_pdf, width = plot_cfg$width %||% 12,
                         height = plot_cfg$height %||% 8, onefile = TRUE, family = "Arial")
  } else {
    grDevices::pdf(each_pdf, width = plot_cfg$width %||% 12,
                   height = plot_cfg$height %||% 8, onefile = TRUE, family = "Arial")
  }

  for (i in seq_along(top_ids)) {
    pp <- gseaplot2(
      gsea_obj,
      geneSetID = top_ids[i],
      color = sci_colors[i],
      pvalue_table = TRUE,
      subplots = 1:3,
      base_size = 10,
      rel_heights = c(1.6, 0.5, 1)
    )
    pp[[1]] <- pp[[1]] +
      theme(
        legend.position = "top",
        text = element_text(family = "Arial"),
        axis.text = element_text(family = "Arial", size = 9),
        axis.title = element_text(family = "Arial", size = 11),
        legend.text = element_text(family = "Arial", size = 9)
      )
    print(pp)
  }
  dev.off()

  log_msg("info", sprintf("%s GSEA plots saved", gene))
}

#' 绘制GSEA结果图（通用版本，用于Pre-rank GSEA）
#' @param gsea_obj GSEA结果对象
#' @param plot_cfg 绘图配置
#' @param top_n 显示的通路数
#' @param plot_dir 图形输出目录
#' @param prefix 文件名前缀
plot_gsea_results <- function(gsea_obj, plot_cfg, top_n = 5, plot_dir, prefix = "GSEA") {

  if (is.null(gsea_obj) || nrow(gsea_obj@result) == 0) {
    log_msg("warn", "无显著富集结果，跳过绘图")
    return(invisible(NULL))
  }

  show_n <- min(plot_cfg$ridge_show_category %||% 10, nrow(gsea_obj@result))
  if (show_n < 1) {
    log_msg("warn", "无可绘制通路，跳过绘图")
    return(invisible(NULL))
  }

  ensure_dir(plot_dir)

  # 绘制山脊图 - Arial字体
  tryCatch({
    p_ridge <- ridgeplot(
      gsea_obj,
      showCategory = show_n,
      fill = "pvalue",
      core_enrichment = TRUE,
      label_format = 30,
      orderBy = "NES",
      decreasing = TRUE
    ) +
      theme_classic(base_size = 10) +
      theme(axis.text.y = element_text(family = "Arial", size = 9),
            axis.text.x = element_text(family = "Arial", size = 9),
            axis.title.x = element_text(family = "Arial", size = 11),
            axis.title.y = element_text(family = "Arial", size = 11),
            plot.title = element_text(family = "Arial", hjust = 0.5, size = 14, face = "bold"),
            legend.text = element_text(family = "Arial", size = 9),
            legend.title = element_text(family = "Arial", size = 10)) +
      labs(title = paste0(prefix, " GSEA Ridge Plot"), x = "Rank", y = "Pathway") +
      scale_fill_gradientn(colors = viridisLite::magma(256, direction = -1, begin = 0.05, end = 0.95),
                           name = "p value")

    ridge_file <- file.path(plot_dir, paste0(prefix, "_RidgePlot_SCI.pdf"))
    ggsave(ridge_file, p_ridge,
           width = plot_cfg$ridge_width %||% 10,
           height = plot_cfg$ridge_height %||% 6)
  }, error = function(e) {
    log_msg("warn", sprintf("Ridge plot failed: %s", e$message))
  })

  # 选择top通路绘制GSEA曲线图
  top_pathways <- gsea_obj@result %>%
    arrange(desc(abs(NES))) %>%
    head(min(top_n, nrow(gsea_obj@result)))
  top_ids <- top_pathways$ID

  if (length(top_ids) == 0) {
    log_msg("warn", "无可绘制富集曲线")
    return(invisible(NULL))
  }

  # SCI配色
  sci_colors <- pal_npg(alpha = 0.9)(max(5, length(top_ids)))[seq_along(top_ids)]

  # 绘制组合GSEA图 - Arial字体
  p_all <- gseaplot2(
    gsea_obj,
    geneSetID = top_ids,
    color = sci_colors,
    pvalue_table = FALSE,
    subplots = 1:3,
    base_size = 9,
    rel_heights = c(1.5, 0.5, 1)
  )
  p_all[[1]] <- p_all[[1]] +
    theme(
      legend.position = "top",
      legend.direction = "vertical",
      text = element_text(family = "Arial"),
      axis.text = element_text(family = "Arial", size = 9),
      axis.title = element_text(family = "Arial", size = 11),
      plot.title = element_text(family = "Arial", face = "bold"),
      legend.text = element_text(family = "Arial", size = 9)
    )

  # 保存PDF - Arial字体
  all_pdf <- file.path(plot_dir, paste0(prefix, "_GSEA_EnrichmentPlot_SCI.pdf"))
  if (capabilities("cairo")) {
    grDevices::cairo_pdf(all_pdf, width = plot_cfg$width %||% 12, height = plot_cfg$height %||% 8, family = "Arial")
  } else {
    grDevices::pdf(all_pdf, width = plot_cfg$width %||% 12, height = plot_cfg$height %||% 8, family = "Arial")
  }
  print(p_all)
  dev.off()

  # 保存PNG（用于报告嵌入）- Arial字体
  all_png <- file.path(plot_dir, paste0(prefix, "_GSEA_EnrichmentPlot_SCI.png"))
  tryCatch({
    png(all_png, width = (plot_cfg$width %||% 12) * 150,
        height = (plot_cfg$height %||% 8) * 150, res = 150, family = "Arial")
    print(p_all)
    dev.off()
  }, error = function(e) {
    log_msg("warn", sprintf("PNG save failed: %s", e$message))
  })

  # 绘制每通路单独的GSEA图 - Arial字体
  each_pdf <- file.path(plot_dir, paste0(prefix, "_GSEA_EnrichmentPlot_SCI_PerPathway.pdf"))
  if (capabilities("cairo")) {
    grDevices::cairo_pdf(each_pdf, width = plot_cfg$width %||% 12,
                         height = plot_cfg$height %||% 8, onefile = TRUE, family = "Arial")
  } else {
    grDevices::pdf(each_pdf, width = plot_cfg$width %||% 12,
                   height = plot_cfg$height %||% 8, onefile = TRUE, family = "Arial")
  }

  for (i in seq_along(top_ids)) {
    pp <- gseaplot2(
      gsea_obj,
      geneSetID = top_ids[i],
      color = sci_colors[i],
      pvalue_table = TRUE,
      subplots = 1:3,
      base_size = 10,
      rel_heights = c(1.6, 0.5, 1)
    )
    pp[[1]] <- pp[[1]] +
      theme(
        legend.position = "top",
        text = element_text(family = "Arial"),
        axis.text = element_text(family = "Arial", size = 9),
        axis.title = element_text(family = "Arial", size = 11),
        legend.text = element_text(family = "Arial", size = 9)
      )
    print(pp)
  }
  dev.off()

  log_msg("info", sprintf("%s GSEA plots saved", prefix))
}

# ---- Pre-rank GSEA 完整流程 ----

#' 运行Pre-rank GSEA分析流程
#' @param cfg 配置列表
#' @param paths 路径列表
#' @return 分析结果列表
run_prerank_gsea_pipeline <- function(cfg, paths) {

  log_msg("info", "========== 开始Pre-rank GSEA分析 ==========")

  # 读取排序基因列表
  ranked_file <- paths$ranked_gene_list
  if (!file.exists(ranked_file)) {
    stop(sprintf("排序基因列表文件不存在: %s", ranked_file))
  }

  rank_metric <- get_cfg(cfg, "analysis.gsea.rank_metric", "log2FC")
  gene_rank <- read_ranked_gene_list(ranked_file, metric_col = rank_metric)

  # 读取GMT文件
  kegg_genesets <- clusterProfiler::read.gmt(paths$gmt_file)
  all_gmt_genes <- unique(kegg_genesets$gene)
  log_msg("info", sprintf("GMT基因集: %d 个基因", length(all_gmt_genes)))

  # 检查GMT匹配率（对齐前，仅作参考）
  gmt_qc_pre <- check_gmt_match_rate(paths$gmt_file, names(gene_rank),
                                     get_cfg(cfg, "qc.min_gmt_match_rate", 30))

  species_hint <- get_cfg(cfg, "project.species", NULL)

  # GSEA配置
  gsea_cfg <- list(
    minGSSize = get_cfg(cfg, "analysis.gsea.minGSSize", 10),
    maxGSSize = get_cfg(cfg, "analysis.gsea.maxGSSize", 500),
    pvalueCutoff = get_cfg(cfg, "analysis.gsea.pvalue_cutoff", 0.05),
    padj_method = get_cfg(cfg, "analysis.gsea.padj_method", "BH"),
    exponent = get_cfg(cfg, "analysis.gsea.exponent", 1),
    eps = get_cfg(cfg, "analysis.gsea.eps", 0)
  )

  # 运行GSEA（内部完成 ID 对齐）
  result <- run_prerank_gsea(gene_rank, kegg_genesets, gsea_cfg, all_gmt_genes,
                             species_hint = species_hint)

  if (is.null(result)) {
    stop("Pre-rank GSEA分析失败")
  }

  gsea_obj <- result$gsea
  gene_rank_used <- result$gene_rank
  term2gene_used <- result$term2gene
  n_gmt_unique <- length(unique(term2gene_used$gene))
  n_matched <- sum(names(gene_rank_used) %in% unique(term2gene_used$gene))
  match_pct <- if (length(gene_rank_used) > 0) 100 * n_matched / length(gene_rank_used) else 0

  # 保存结果
  result_df <- gsea_obj@result
  result_file <- file.path(paths$gsea_dir, "prerank_gsea_results.csv")
  write_csv_utf8(result_df, result_file)

  if (nrow(result_df) > 0) {
    # 保存详细报告
    detail <- result_df[, c("ID", "setSize", "NES", "pvalue", "p.adjust", "core_enrichment")]
    detail_file <- file.path(paths$tables_dir, "PreRank_GSEA_DetailedReport.csv")
    write_csv_utf8(detail, detail_file)
  }

  # 绘图
  plot_cfg <- list(
    ridge_show_category = get_cfg(cfg, "plot.ridge_show_category", 10),
    ridge_width = get_cfg(cfg, "plot.ridge_width", 10),
    ridge_height = get_cfg(cfg, "plot.ridge_height", 6),
    width = get_cfg(cfg, "plot.width", 10),
    height = get_cfg(cfg, "plot.height", 8)
  )

  plot_gsea_results(gsea_obj, plot_cfg,
                    get_cfg(cfg, "analysis.gsea.top_pathway_n", 5),
                    paths$plots_dir, prefix = "PreRank")

  # 保存质控报告（与 GSEA 实际使用的对齐后数据一致；匹配率=交集/输入基因数）
  qc_overall <- data.frame(
    Item = c(
      "Input_Genes", "GMT_Genes_Total", "GMT_Genes_Matched", "GMT_Match_Rate",
      "ID_Converted", "Species",
      "Ranking_Metric", "Metric_Range", "Significant_Pathways", "Analysis_Date"
    ),
    Value = c(
      length(gene_rank_used),
      n_gmt_unique,
      n_matched,
      sprintf("%.2f%%", match_pct),
      ifelse(isTRUE(result$id_converted), "true", "false"),
      as.character(result$species_used %||% "hsa"),
      rank_metric,
      if (length(gene_rank_used) > 0) {
        sprintf("[%.4f, %.4f]", min(gene_rank_used), max(gene_rank_used))
      } else {
        "NA"
      },
      nrow(gsea_obj@result),
      as.character(Sys.Date())
    ),
    stringsAsFactors = FALSE
  )
  write_csv_utf8(qc_overall, file.path(paths$qc_dir, "PreRank_GSEA_QC_Report.csv"))
  write_csv_utf8(qc_overall, file.path(paths$tables_dir, "PreRank_GSEA_QC_Report.csv"))

  log_msg("info", sprintf("Pre-rank GSEA完成: %d 个显著通路", nrow(gsea_obj@result)))

  list(
    gsea = gsea_obj,
    gene_rank = gene_rank_used,
    gmt_qc = gmt_qc_pre,
    qc = result$qc
  )
}

# ---- Correlation GSEA 完整流程 ----

#' 运行基于相关性的GSEA分析流程
#' @param cfg 配置列表
#' @param paths 路径列表
#' @return 分析结果列表
run_correlation_gsea_pipeline <- function(cfg, paths) {

  log_msg("info", "========== 开始Correlation GSEA分析 ==========")

  # 读取目标基因
  key_genes <- read_gene_list(paths$target_genes)
  log_msg("info", sprintf("读取目标基因: %d 个", length(key_genes)))

  # 读取GMT文件
  if (!file.exists(paths$gmt_file)) {
    stop(sprintf("GMT文件不存在: %s", paths$gmt_file))
  }
  kegg_genesets <- clusterProfiler::read.gmt(paths$gmt_file)
  all_gmt_genes <- unique(kegg_genesets$gene)
  log_msg("info", sprintf("GMT基因集: %d 个基因", length(all_gmt_genes)))

  # 表达矩阵质控
  qc_exp <- qc_expression_matrix(
    exp_file = paths$expression_matrix,
    variance_threshold = get_cfg(cfg, "qc.variance_threshold", 0.01),
    na_threshold = get_cfg(cfg, "qc.na_ratio_threshold", 0.5),
    key_genes = key_genes
  )

  exp <- qc_exp$exp

  # 推断物种和ID类型（行=基因）
  expr_info <- infer_expr_id(rownames(exp))
  species <- get_cfg(cfg, "project.species", NULL)
  if (!is.null(species)) {
    expr_info$species <- species
  }
  id_type <- get_cfg(cfg, "project.default_id_type", NULL)
  if (!is.null(id_type)) {
    expr_info$id_type <- id_type
  }

  # 验证GMT与表达矩阵一致性
  gmt_info <- infer_gmt_meta(paths$gmt_file)
  validate_species_id(expr_info, gmt_info)

  # 检查GMT匹配率（行=基因）
  gmt_qc <- check_gmt_match_rate(paths$gmt_file, rownames(exp),
                                  get_cfg(cfg, "qc.min_gmt_match_rate", 30))

  # 计算相关性
  cor_method <- get_cfg(cfg, "analysis.gsea.correlation_method", "spearman")
  cor_results <- calculate_correlation(exp, key_genes, cor_method)
  if (length(cor_results) == 0) {
    stop("没有可用于分析的目标基因（都不在表达矩阵中）")
  }

  # GSEA配置
  gsea_cfg <- list(
    minGSSize = get_cfg(cfg, "analysis.gsea.minGSSize", 10),
    maxGSSize = get_cfg(cfg, "analysis.gsea.maxGSSize", 500),
    pvalueCutoff = get_cfg(cfg, "analysis.gsea.pvalue_cutoff", 0.05),
    padj_method = get_cfg(cfg, "analysis.gsea.padj_method", "BH"),
    exponent = get_cfg(cfg, "analysis.gsea.exponent", 1),
    eps = get_cfg(cfg, "analysis.gsea.eps", 0)
  )

  species_hint_gsea <- expr_info$species
  if (identical(species_hint_gsea, "unknown")) species_hint_gsea <- NULL

  # 运行GSEA
  gsea_results <- list()
  gene_qc_list <- list()

  for (gene in names(cor_results)) {
    log_msg("info", sprintf("开始处理: %s", gene))

    one <- tryCatch({
      run_gsea_for_gene(gene, cor_results[[gene]], kegg_genesets,
                        gsea_cfg, all_gmt_genes, species_hint = species_hint_gsea)
    }, error = function(e) {
      log_msg("error", sprintf("%s GSEA失败: %s", gene, e$message))
      NULL
    })

    if (is.null(one)) next

    gsea_results[[gene]] <- one$gsea
    gene_qc_list[[gene]] <- one$gene_qc

    # 保存结果和绘图
    tryCatch({
      save_gene_outputs(gene, one$gsea, paths)
      plot_gene_results(gene, one$gsea,
                        list(ridge_show_category = get_cfg(cfg, "plot.ridge_show_category", 10),
                             ridge_width = get_cfg(cfg, "plot.ridge_width", 10),
                             ridge_height = get_cfg(cfg, "plot.ridge_height", 6),
                             width = get_cfg(cfg, "plot.width", 10),
                             height = get_cfg(cfg, "plot.height", 8)),
                        get_cfg(cfg, "analysis.gsea.top_pathway_n", 5),
                        paths$plots_dir)
    }, error = function(e) {
      log_msg("error", sprintf("%s 结果保存或绘图失败: %s", gene, e$message))
    })
  }

  # 保存基因级质控报告
  if (length(gene_qc_list) > 0) {
    gene_qc_df <- do.call(rbind, gene_qc_list)
    rownames(gene_qc_df) <- NULL
    write_csv_utf8(gene_qc_df, file.path(paths$qc_dir, "Correlation_GSEA_QC_Report_ByGene.csv"))
    write_csv_utf8(gene_qc_df, file.path(paths$tables_dir, "Correlation_GSEA_QC_Report_ByGene.csv"))
  }

  # 保存整体质控报告
  qc_overall <- data.frame(
    Item = c("Raw_Expression_Genes", "Filtered_Expression_Genes", "Filtered_Genes",
             "GMT_Genes_Total", "Matched_Genes", "Match_Rate",
             "Analysis_Date", "Species", "ID_Type"),
    Value = c(
      qc_exp$raw_exp_gene_count,
      qc_exp$filtered_exp_gene_count,
      length(qc_exp$bad_genes),
      length(all_gmt_genes),
      length(gmt_qc$matched_genes),
      sprintf("%.2f%%", gmt_qc$match_rate),
      as.character(Sys.Date()),
      expr_info$species,
      expr_info$id_type
    ),
    stringsAsFactors = FALSE
  )
  write_csv_utf8(qc_overall, file.path(paths$qc_dir, "Correlation_GSEA_QC_Report_Overall.csv"))
  write_csv_utf8(qc_overall, file.path(paths$tables_dir, "Correlation_GSEA_QC_Report_Overall.csv"))

  log_msg("info", "========== Correlation GSEA分析完成 ==========")

  list(
    gsea_results = gsea_results,
    cor_results = cor_results,
    qc_exp = qc_exp,
    gmt_qc = gmt_qc,
    gene_qc = gene_qc_list
  )
}
# ---- Median Split GSEA 完整流程 ----
#' 运行基于中位数分组的GSEA分析流程
#' @param cfg 配置列表
#' @param paths 路径列表
#' @return 分析结果列表
run_median_split_gsea_pipeline <- function(cfg, paths) {
  log_msg("info", "========== 开始Median Split GSEA分析 ==========")
  # 读取目标基因
  key_genes <- read_gene_list(paths$target_genes)
  log_msg("info", sprintf("读取目标基因: %d 个", length(key_genes)))
  # 读取GMT文件
  if (!file.exists(paths$gmt_file)) {
    stop(sprintf("GMT文件不存在: %s", paths$gmt_file))
  }
  kegg_genesets <- clusterProfiler::read.gmt(paths$gmt_file)
  all_gmt_genes <- unique(kegg_genesets$gene)
  log_msg("info", sprintf("GMT基因集: %d 个基因", length(all_gmt_genes)))
  # 表达矩阵质控
  qc_exp <- qc_expression_matrix(
    exp_file = paths$expression_matrix,
    variance_threshold = get_cfg(cfg, "qc.variance_threshold", 0.01),
    na_threshold = get_cfg(cfg, "qc.na_ratio_threshold", 0.5),
    key_genes = key_genes
  )
  exp <- qc_exp$exp
  # 推断物种和ID类型（行=基因）
  expr_info <- infer_expr_id(rownames(exp))
  species <- get_cfg(cfg, "project.species", NULL)
  if (!is.null(species)) {
    expr_info$species <- species
  }
  id_type <- get_cfg(cfg, "project.default_id_type", NULL)
  if (!is.null(id_type)) {
    expr_info$id_type <- id_type
  }
  # 验证GMT与表达矩阵一致性
  gmt_info <- infer_gmt_meta(paths$gmt_file)
  validate_species_id(expr_info, gmt_info)
  # 检查GMT匹配率（行=基因）
  gmt_qc <- check_gmt_match_rate(paths$gmt_file, rownames(exp),
                                  get_cfg(cfg, "qc.min_gmt_match_rate", 30))
  # GSEA配置
  gsea_cfg <- list(
    minGSSize = get_cfg(cfg, "analysis.gsea.minGSSize", 10),
    maxGSSize = get_cfg(cfg, "analysis.gsea.maxGSSize", 500),
    pvalueCutoff = get_cfg(cfg, "analysis.gsea.pvalue_cutoff", 0.05),
    padj_method = get_cfg(cfg, "analysis.gsea.padj_method", "BH"),
    exponent = get_cfg(cfg, "analysis.gsea.exponent", 1),
    eps = get_cfg(cfg, "analysis.gsea.eps", 0)
  )
  species_hint_gsea <- expr_info$species
  if (identical(species_hint_gsea, "unknown")) species_hint_gsea <- NULL
  # 运行GSEA
  gsea_results <- list()
  gene_qc_list <- list()
  for (gene in key_genes) {
    # 检查基因是否在表达矩阵中（行=基因）
    if (!gene %in% rownames(exp)) {
      log_msg("warn", sprintf("跳过不存在于矩阵的目标基因: %s", gene))
      next
    }
    log_msg("info", sprintf("开始处理: %s", gene))
    # 获取目标基因的表达值（行=基因，所以取行）
    gene_expression <- exp[gene, ]
    # 按中位数将样本分为高表达组和低表达组
    median_expr <- median(gene_expression)
    high_group <- gene_expression > median_expr
    low_group <- gene_expression <= median_expr
    n_high <- sum(high_group)
    n_low <- sum(low_group)
    log_msg("info", sprintf("%s: 中位数=%.4f, 高表达组=%d, 低表达组=%d",
                gene, median_expr, n_high, n_low))
    if (n_high < 3 || n_low < 3) {
      log_msg("warn", sprintf("%s: 样本数不足，跳过", gene))
      next
    }
    # 计算所有基因的log2FC（高表达组 vs 低表达组）
    # 行=基因，列=样本，所以按列取样本子集，按行计算均值
    log2FC <- rowMeans(exp[, high_group]) - rowMeans(exp[, low_group])
    # 处理NA和无限值
    log2FC <- log2FC[!is.na(log2FC) & is.finite(log2FC)]
    if (length(log2FC) == 0) {
      log_msg("warn", sprintf("%s: log2FC计算失败", gene))
      next
    }
    # 处理重复基因名
    if (any(duplicated(names(log2FC)))) {
      dup_genes <- names(log2FC)[duplicated(names(log2FC))]
      log_msg("warn", sprintf("发现%d个重复基因，保留最高值", length(dup_genes)))
      log2FC <- log2FC[!duplicated(names(log2FC))]
    }
    # 降序排列
    gene_rank <- sort(log2FC, decreasing = TRUE)
    log_msg("info", sprintf("%s: log2FC范围 [%.4f, %.4f], %d个基因",
              gene, min(gene_rank), max(gene_rank), length(gene_rank)))
    # 运行GSEA
    one <- tryCatch({
      run_gsea_for_gene(gene, gene_rank, kegg_genesets,
                      gsea_cfg, all_gmt_genes, species_hint = species_hint_gsea)
    }, error = function(e) {
      log_msg("error", sprintf("%s GSEA失败: %s", gene, e$message))
      NULL
    })
    if (is.null(one)) next
    gsea_results[[gene]] <- one$gsea
    # 更新QC信息
    gene_qc <- one$gene_qc
    gene_qc$Median_Expression <- median_expr
    gene_qc$High_Group_N <- n_high
    gene_qc$Low_Group_N <- n_low
    gru <- one$gene_rank_used
    gene_qc$Log2FC_Range <- if (length(gru) > 0) {
      sprintf("[%.4f, %.4f]", min(gru), max(gru))
    } else {
      sprintf("[%.4f, %.4f]", min(gene_rank), max(gene_rank))
    }
    gene_qc_list[[gene]] <- gene_qc
    # 保存结果和绘图
    tryCatch({
      save_gene_outputs(gene, one$gsea, paths)
      plot_gene_results(gene, one$gsea,
                        list(ridge_show_category = get_cfg(cfg, "plot.ridge_show_category", 10),
                             ridge_width = get_cfg(cfg, "plot.ridge_width", 10),
                             ridge_height = get_cfg(cfg, "plot.ridge_height", 6),
                             width = get_cfg(cfg, "plot.width", 10),
                             height = get_cfg(cfg, "plot.height", 8)),
                        get_cfg(cfg, "analysis.gsea.top_pathway_n", 5),
                        paths$plots_dir)
    }, error = function(e) {
      log_msg("error", sprintf("%s 结果保存或绘图失败: %s", gene, e$message))
    })
  }
  # 保存基因级质控报告
  if (length(gene_qc_list) > 0) {
    gene_qc_df <- do.call(rbind, gene_qc_list)
    rownames(gene_qc_df) <- NULL
    write_csv_utf8(gene_qc_df, file.path(paths$qc_dir, "MedianSplit_GSEA_QC_Report_ByGene.csv"))
    write_csv_utf8(gene_qc_df, file.path(paths$tables_dir, "MedianSplit_GSEA_QC_Report_ByGene.csv"))
  }
  # 保存整体质控报告
  qc_overall <- data.frame(
    Item = c("Raw_Expression_Genes", "Filtered_Expression_Genes", "Filtered_Genes",
             "GMT_Genes_Total", "Matched_Genes", "Match_Rate",
             "Analysis_Date", "Species", "ID_Type"),
    Value = c(
      qc_exp$raw_exp_gene_count,
      qc_exp$filtered_exp_gene_count,
      length(qc_exp$bad_genes),
      length(all_gmt_genes),
      length(gmt_qc$matched_genes),
      sprintf("%.2f%%", gmt_qc$match_rate),
      as.character(Sys.Date()),
      expr_info$species,
      expr_info$id_type
    ),
    stringsAsFactors = FALSE
  )
  write_csv_utf8(qc_overall, file.path(paths$qc_dir, "MedianSplit_GSEA_QC_Report_Overall.csv"))
  write_csv_utf8(qc_overall, file.path(paths$tables_dir, "MedianSplit_GSEA_QC_Report_Overall.csv"))
  log_msg("info", "========== Median Split GSEA分析完成 ==========")
  list(
    gsea_results = gsea_results,
    gene_ranks = NULL,  # 可选： 保存log2FC用于调试
    qc_exp = qc_exp,
    gmt_qc = gmt_qc,
    gene_qc = gene_qc_list
  )
}
# ---- 统一GSEA入口 ----

#' 运行完整的GSEA分析流程（统一入口）
#' @param cfg 配置列表
#' @param paths 路径列表
#' @return 分析结果列表
run_gsea_pipeline <- function(cfg, paths) {

  method <- get_cfg(cfg, "analysis.gsea.method", "prerank")
  log_msg("info", sprintf("GSEA方法: %s", method))

  if (method == "prerank") {
    run_prerank_gsea_pipeline(cfg, paths)
  } else if (method == "correlation") {
    run_correlation_gsea_pipeline(cfg, paths)
  } else if (method == "median_split") {
    run_median_split_gsea_pipeline(cfg, paths)
  } else {
    stop(sprintf("不支持的GSEA方法: %s (支持: prerank, correlation, median_split)", method))
  }
}

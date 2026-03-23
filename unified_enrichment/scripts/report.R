# ==============================================================================
# unified_enrichment/scripts/report.R
# 统一报告生成模块
# 合并自 GO_KEGG/code/generate_report.R 和 单基因富集分析/scripts/report_generator.R
# ==============================================================================

suppressPackageStartupMessages({
  library(officer)
  library(flextable)
  library(dplyr)
  library(yaml)
  library(stringr)
  library(data.table)
})

source("scripts/utils.R")

# ---- 辅助函数 ----

#' 获取R包版本
#' @param pkg 包名
#' @return 版本字符串
get_pkg_version <- function(pkg) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    as.character(utils::packageVersion(pkg))
  } else {
    "not installed"
  }
}

#' 获取物种英文名称
#' @param species 物种代码
#' @return 物种英文名称
get_species_name_en <- function(species) {
  species_names <- list(
    hsa = "Homo sapiens (human)",
    mmu = "Mus musculus (mouse)",
    rno = "Rattus norvegicus (rat)"
  )
  species_names[[species]] %||% species
}

#' 创建三线表样式的flextable
#' @param df 数据框
#' @param caption 表格标题
#' @return flextable对象
create_booktabs_table <- function(df, caption = "") {
  ft <- flextable(df)
  ft <- theme_booktabs(ft)
  ft <- fontsize(ft, size = 9, part = "all")
  ft <- font(ft, fontname = "Times New Roman", part = "all")
  ft <- align(ft, align = "left", part = "body")
  ft <- align(ft, align = "center", part = "header")

  # 设置三线表样式：顶线和底线加粗，表头加粗并添加底框
  ft <- border_remove(ft)
  ft <- hline_top(ft, part = "header", border = officer::fp_border(width = 1.5))
  ft <- hline_top(ft, part = "body", border = officer::fp_border(width = 1.5))
  ft <- hline_bottom(ft, part = "body", border = officer::fp_border(width = 1.5))
  ft <- hline_bottom(ft, part = "footer", border = officer::fp_border(width = 1.5))

  # 表头样式
  ft <- bold(ft, part = "header")

  if (nchar(caption) > 0) {
    ft <- set_caption(ft, caption)
  }
  ft
}

# ---- GO/KEGG方法部分 ----

#' 生成GO/KEGG分析方法描述 (英文)
#' @param cfg 配置列表
#' @return 方法描述文本
generate_go_kegg_methods <- function(cfg) {
  species_name <- get_species_name_en(get_cfg(cfg, "project.species", "hsa"))
  go_ont <- get_cfg(cfg, "analysis.go_kegg.go_ont", "ALL")
  pvalue_cutoff <- get_cfg(cfg, "analysis.go_kegg.pvalue_cutoff", 0.05)
  padj_method <- get_cfg(cfg, "analysis.go_kegg.padj_method", "BH")

  methods <- sprintf(
    "Gene Ontology (GO) and Kyoto Encyclopedia of Genes and Genomes (KEGG) enrichment analysis were performed using the clusterProfiler R package. For GO enrichment analysis, the enrichGO function was applied to the %s genome with the ontology set to %s. KEGG pathway enrichment analysis was conducted using the enrichKEGG function. Statistical significance was determined using %s-adjusted p-value (p < %.2f). Key R packages used: clusterProfiler (v%s), enrichplot (v%s).",
    species_name,
    go_ont,
    padj_method,
    pvalue_cutoff,
    get_pkg_version("clusterProfiler"),
    get_pkg_version("enrichplot")
  )

  methods
}

# ---- GSEA方法部分 ----

#' 生成GSEA分析方法描述 (英文)
#' @param cfg 配置列表
#' @return 方法描述文本
generate_gsea_methods <- function(cfg) {
  gmt_name <- basename(get_cfg(cfg, "input.gmt_file", "gene_sets.gmt"))
  gsea_method <- get_cfg(cfg, "analysis.gsea.method", "prerank")
  minGSSize <- get_cfg(cfg, "analysis.gsea.minGSSize", 10)
  maxGSSize <- get_cfg(cfg, "analysis.gsea.maxGSSize", 500)
  pvalue_cutoff <- get_cfg(cfg, "analysis.gsea.pvalue_cutoff", 0.05)
  padj_method <- get_cfg(cfg, "analysis.gsea.padj_method", "BH")

  methods <- sprintf(
    "Gene Set Enrichment Analysis (GSEA) was performed to identify biological pathways associated with target genes. For pre-ranked GSEA, genes were ranked by %s (e.g., log2FC) from high to low. GSEA was then conducted using the %s gene set database with the following parameters: minGSSize=%d, maxGSSize=%d, p-value cutoff=%.4f, adjustment method=%s. Key R packages used: clusterProfiler (v%s), enrichplot (v%s), fgsea (v%s).",
    get_cfg(cfg, "analysis.gsea.rank_metric", "log2FC"),
    gmt_name,
    minGSSize,
    maxGSSize,
    pvalue_cutoff,
    padj_method,
    get_pkg_version("clusterProfiler"),
    get_pkg_version("enrichplot"),
    get_pkg_version("fgsea")
  )

  methods
}

# ---- GO/KEGG结果部分 ----

#' 生成GO/KEGG分析结果描述 (英文)
#' @param go_results GO结果
#' @param kegg_results KEGG结果
#' @param value_col 值列名
#' @return 结果描述文本
generate_go_kegg_results <- function(go_results, kegg_results, value_col = "p.adjust") {

  # 提取top描述
  extract_top_desc <- function(df, n = 5, col = "Description", value_col = "p.adjust") {
    if (is.null(df) || nrow(df) == 0) return("None")
    ord <- order(df[[value_col]], na.last = NA)
    if (length(ord) == 0) return("None")
    df_sorted <- df[ord, , drop = FALSE]
    top_desc <- utils::head(df_sorted[[col]], n)
    top_desc <- unique(stats::na.omit(as.character(top_desc)))
    if (length(top_desc) == 0) return("None")
    paste(top_desc, collapse = ", ")
  }

  # GO统计
  total_go <- if (!is.null(go_results)) nrow(go_results) else 0
  bp_count <- if (!is.null(go_results)) sum(go_results$ONTOLOGY == "BP") else 0
  cc_count <- if (!is.null(go_results)) sum(go_results$ONTOLOGY == "CC") else 0
  mf_count <- if (!is.null(go_results)) sum(go_results$ONTOLOGY == "MF") else 0

  # KEGG统计
  total_kegg <- if (!is.null(kegg_results)) nrow(kegg_results) else 0

  # 提取top条目
  top_go_bp <- extract_top_desc(go_results[go_results$ONTOLOGY == "BP", ], 5, value_col = value_col)
  top_go_cc <- extract_top_desc(go_results[go_results$ONTOLOGY == "CC", ], 5, value_col = value_col)
  top_go_mf <- extract_top_desc(go_results[go_results$ONTOLOGY == "MF", ], 5, value_col = value_col)
  top_kegg <- extract_top_desc(kegg_results, 5, value_col = value_col)

  # 构建结果文本
  go_text <- sprintf(
    "GO enrichment analysis identified %d significant GO terms, including %d BP (Biological Process), %d CC (Cellular Component), and %d MF (Molecular Function) terms. The top BP terms include: %s. The top CC terms include: %s. The top MF terms include: %s.",
    total_go, bp_count, cc_count, mf_count, top_go_bp, top_go_cc, top_go_mf
  )

  kegg_text <- sprintf(
    "KEGG pathway analysis identified %d significant pathways. The top pathways include: %s.",
    total_kegg, top_kegg
  )

  list(go = go_text, kegg = kegg_text)
}

# ---- GSEA结果部分 ----

#' 生成单个基因GSEA结果描述 (英文)
#' @param gene_name 基因名
#' @param gsea_detail GSEA详细结果
#' @return 结果描述文本
generate_gsea_results <- function(gene_name, gsea_detail) {
  if (is.null(gsea_detail) || nrow(gsea_detail) == 0) {
    return(sprintf("No significant enriched pathways were found for gene %s (FDR < 0.05).", gene_name))
  }

  sig_n <- nrow(gsea_detail)
  top3 <- gsea_detail %>%
    arrange(desc(abs(NES))) %>%
    head(3)

  if (nrow(top3) == 0) {
    return(sprintf("No significant enriched pathways were found for gene %s.", gene_name))
  }

  pathway_names <- top3$ID %>%
    str_replace("^KEGG_", "") %>%
    str_replace("_", " ")

  nes_values <- round(top3$NES, 2)

  # 根据通路数量构建描述
  if (nrow(top3) >= 3) {
    results <- sprintf(
      "GSEA analysis for gene %s identified %d significantly enriched pathways. The top pathways include: %s (NES=%.2f), %s (NES=%.2f), %s (NES=%.2f), etc.",
      gene_name, sig_n,
      pathway_names[1], nes_values[1],
      pathway_names[2], nes_values[2],
      pathway_names[3], nes_values[3]
    )
  } else if (nrow(top3) >= 2) {
    results <- sprintf(
      "GSEA analysis for gene %s identified %d significantly enriched pathways. The top pathways include: %s (NES=%.2f), %s (NES=%.2f), etc.",
      gene_name, sig_n,
      pathway_names[1], nes_values[1],
      pathway_names[2], nes_values[2]
    )
  } else {
    results <- sprintf(
      "GSEA analysis for gene %s identified %d significantly enriched pathway: %s (NES=%.2f).",
      gene_name, sig_n,
      pathway_names[1], nes_values[1]
    )
  }

  results
}

# ---- 添加图像到报告 ----

#' 添加GSEA图像到报告 (英文)
#' @param doc Word文档对象
#' @param gene_name 基因名
#' @param plots_dir 图形目录
#' @return 更新后的文档对象
add_gsea_figure_to_report <- function(doc, gene_name, plots_dir) {
  png_file <- file.path(plots_dir, paste0(gene_name, "_GSEA_EnrichmentPlot_SCI.png"))

  if (!file.exists(png_file)) {
    log_msg("warn", sprintf("%s PNG enrichment plot does not exist, skipping", gene_name))
    return(doc)
  }

  tryCatch({
    doc <- body_add_img(doc, src = png_file, width = 6, height = 4)
    # 添加英文图题，使用Arial字体
    doc <- body_add_par(doc, sprintf("Figure %d. GSEA enrichment plot for %s", 1, gene_name), style = "Normal")
    doc <- body_add_par(doc, "", style = "Normal")
  }, error = function(e) {
    log_msg("warn", sprintf("%s image add failed: %s", gene_name, e$message))
  })

  doc
}

#' 添加GO/KEGG图像到报告 (英文)
#' @param doc Word文档对象
#' @param plot_file 图形文件路径
#' @param caption 图注
#' @return 更新后的文档对象
add_enrichment_figure_to_report <- function(doc, plot_file, caption) {
  if (!file.exists(plot_file)) {
    log_msg("warn", sprintf("Plot file does not exist: %s", plot_file))
    return(doc)
  }

  # 优先使用PNG
  png_file <- sub("\\.pdf$", ".png", plot_file)
  if (file.exists(png_file)) {
    plot_file <- png_file
  }

  tryCatch({
    doc <- body_add_img(doc, src = plot_file, width = 6, height = 4.5)
    doc <- body_add_par(doc, caption, style = "Normal")
    doc <- body_add_par(doc, "", style = "Normal")
  }, error = function(e) {
    log_msg("warn", sprintf("Plot add failed: %s", e$message))
  })

  doc
}

# ---- 生成统一报告 ----

#' 生成统一的富集分析报告
#' @param cfg 配置列表
#' @param paths 路径列表
#' @param go_kegg_results GO/KEGG分析结果（可选）
#' @param gsea_results GSEA分析结果（可选）
#' @return 报告文件路径
generate_unified_report <- function(cfg, paths, go_kegg_results = NULL, gsea_results = NULL) {

  log_msg("info", "Generating unified report")

  # 创建Word文档
  doc <- read_docx()

  # 标题 (英文)
  doc <- body_add_par(doc, "Enrichment Analysis Report", style = "heading 1")
  doc <- body_add_par(doc, sprintf("Analysis Date: %s", Sys.Date()), style = "Normal")
  doc <- body_add_par(doc, sprintf("Species: %s", get_species_name_en(cfg$project$species)), style = "Normal")
  doc <- body_add_par(doc, "", style = "Normal")

  # ===== 1. Methods =====
  doc <- body_add_par(doc, "1. Methods", style = "heading 2")

  if (cfg$analysis$go_kegg$enabled) {
    doc <- body_add_par(doc, "1.1 GO/KEGG Enrichment Analysis", style = "heading 3")
    doc <- body_add_par(doc, generate_go_kegg_methods(cfg), style = "Normal")
    doc <- body_add_par(doc, "", style = "Normal")
  }

  if (cfg$analysis$gsea$enabled) {
    doc <- body_add_par(doc, "1.2 GSEA Analysis", style = "heading 3")
    doc <- body_add_par(doc, generate_gsea_methods(cfg), style = "Normal")
    doc <- body_add_par(doc, "", style = "Normal")
  }

  # ===== 2. Results =====
  doc <- body_add_par(doc, "2. Results", style = "heading 2")

  # 质控部分 (Quality Control)
  doc <- body_add_par(doc, "2.1 Quality Control", style = "heading 3")

  # 根据 GSEA 方法动态选择质控文件名
  gsea_method <- get_cfg(cfg, "analysis.gsea.method", "prerank")
  qc_file_map <- list(
    prerank = "PreRank_GSEA_QC_Report.csv",
    correlation = "Correlation_GSEA_QC_Report_Overall.csv",
    median_split = "MedianSplit_GSEA_QC_Report_Overall.csv"
  )
  qc_file <- file.path(paths$qc_dir, qc_file_map[[gsea_method]] %||% "PreRank_GSEA_QC_Report.csv")

  if (file.exists(qc_file)) {
    qc_df <- fread(qc_file, data.table = FALSE)
    doc <- body_add_par(doc, "Overall QC Statistics:", style = "Normal")
    ft <- create_booktabs_table(qc_df)
    doc <- body_add_flextable(doc, ft)
    doc <- body_add_par(doc, "", style = "Normal")
  }

  # GO/KEGG结果部分
  if (cfg$analysis$go_kegg$enabled && !is.null(go_kegg_results)) {
    doc <- body_add_par(doc, "2.2 GO Enrichment Analysis Results", style = "heading 3")

    value_col <- if (cfg$analysis$go_kegg$padj_method == "BH") "p.adjust" else "pvalue"

    # 读取GO结果
    go_results_file <- file.path(paths$go_dir, "go_enrichment_results.csv")
    if (file.exists(go_results_file)) {
      go_results <- fread(go_results_file)

      if (nrow(go_results) > 0) {
        # 结果描述 (英文)
        result_text <- generate_go_kegg_results(go_results, NULL, value_col)
        doc <- body_add_par(doc, result_text$go, style = "Normal")
        doc <- body_add_par(doc, "", style = "Normal")

        # Top5表格
        go_summary <- go_results %>%
          group_by(ONTOLOGY) %>%
          arrange(!!sym(value_col)) %>%
          slice_head(n = 5) %>%
          ungroup() %>%
          dplyr::select(ONTOLOGY, ID, Description, GeneRatio, BgRatio, Count) %>%
          as.data.frame()

        if (nrow(go_summary) > 0) {
          ft <- create_booktabs_table(go_summary, "Table 1. Top 5 GO Enrichment Terms")
          doc <- body_add_flextable(doc, ft)
          doc <- body_add_par(doc, "", style = "Normal")
        }

        # GO图
        go_plot <- file.path(paths$go_dir, "go_enrichment_plot.png")
        if (file.exists(go_plot)) {
          doc <- add_enrichment_figure_to_report(doc, go_plot, "Figure 1. GO Enrichment Dotplot")
        }
      } else {
        doc <- body_add_par(doc, "No significant GO terms were found.", style = "Normal")
      }
    }

    # KEGG结果
    doc <- body_add_par(doc, "2.3 KEGG Pathway Analysis Results", style = "heading 3")

    kegg_results_file <- file.path(paths$kegg_dir, "kegg_enrichment_results.csv")
    if (file.exists(kegg_results_file)) {
      kegg_results <- fread(kegg_results_file)

      if (nrow(kegg_results) > 0) {
        result_text <- generate_go_kegg_results(NULL, kegg_results, value_col)
        doc <- body_add_par(doc, result_text$kegg, style = "Normal")
        doc <- body_add_par(doc, "", style = "Normal")

        # Top5表格
        kegg_summary <- kegg_results %>%
          arrange(!!sym(value_col)) %>%
          slice_head(n = 5) %>%
          dplyr::select(ID, Description, GeneRatio, BgRatio, Count) %>%
          as.data.frame()

        if (nrow(kegg_summary) > 0) {
          ft <- create_booktabs_table(kegg_summary, "Table 2. Top 5 KEGG Pathways")
          doc <- body_add_flextable(doc, ft)
          doc <- body_add_par(doc, "", style = "Normal")
        }

        # KEGG图
        kegg_plot <- file.path(paths$kegg_dir, "kegg_enrichment_plot.png")
        if (file.exists(kegg_plot)) {
          doc <- add_enrichment_figure_to_report(doc, kegg_plot, "Figure 2. KEGG Pathway Enrichment Dotplot")
        }
      } else {
        doc <- body_add_par(doc, "No significant KEGG pathways were found.", style = "Normal")
      }
    }
  }

  # GSEA结果部分
  if (cfg$analysis$gsea$enabled && !is.null(gsea_results)) {
    doc <- body_add_par(doc, "2.4 GSEA Results", style = "heading 3")

    # 根据 GSEA 方法动态选择基因级质控文件名
    gsea_method <- get_cfg(cfg, "analysis.gsea.method", "prerank")
    gene_qc_file_map <- list(
      prerank = NULL,  # Pre-rank 没有基因级质控
      correlation = "Correlation_GSEA_QC_Report_ByGene.csv",
      median_split = "MedianSplit_GSEA_QC_Report_ByGene.csv"
    )
    gene_qc_file <- gene_qc_file_map[[gsea_method]]

    if (!is.null(gene_qc_file)) {
      gene_qc_path <- file.path(paths$qc_dir, gene_qc_file)
      if (file.exists(gene_qc_path)) {
        gene_qc <- fread(gene_qc_path, data.table = FALSE)
        target_genes <- gene_qc$Target_Gene

        fig_num <- 3  # 起始图号

        for (gene in target_genes) {
          doc <- body_add_par(doc, sprintf("Gene: %s", gene), style = "heading 4")

          # 读取GSEA详细结果
          gsea_detail_file <- file.path(paths$gsea_dir, paste0(gene, "_GSEA_DetailedReport.csv"))
          if (file.exists(gsea_detail_file)) {
            gsea_detail <- fread(gsea_detail_file, data.table = FALSE)
            gsea_detail <- gsea_detail %>%
              filter(p.adjust < 0.05) %>%
              arrange(p.adjust)

            # 结果描述 (英文)
            result_text <- generate_gsea_results(gene, gsea_detail)
            doc <- body_add_par(doc, result_text, style = "Normal")
            doc <- body_add_par(doc, "", style = "Normal")

            # Top5表格
            if (nrow(gsea_detail) > 0) {
              top5 <- head(gsea_detail, 5)
              ft <- create_booktabs_table(top5[, c("ID", "setSize", "NES", "pvalue", "p.adjust")],
                                          sprintf("Table. Top 5 GSEA Pathways for %s", gene))
              doc <- body_add_flextable(doc, ft)
              doc <- body_add_par(doc, "", style = "Normal")
            }

            # GSEA图
            doc <- add_gsea_figure_to_report(doc, gene, paths$plots_dir)
          }
        }
      }
    }
  }

  # 保存报告
  report_file <- file.path(paths$results_dir, "enrichment_report.docx")
  print(doc, target = report_file)
  log_msg("info", sprintf("Report generated: %s", report_file))

  invisible(report_file)
}

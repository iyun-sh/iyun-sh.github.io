#!/usr/bin/env Rscript
# GSEA CLI - Command Line Interface
# Methods: prerank, correlation, median_split

options(scipen = 5)
options(BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
options("repos" = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN"))

suppressPackageStartupMessages({
  library(optparse)
  library(clusterProfiler)
  library(enrichplot)
  library(ggplot2)
  library(ggsci)
  library(viridisLite)
  library(dplyr)
  library(data.table)
})

# ---- Utils ----
log_msg <- function(level, message) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[%s] [%s] %s\n", timestamp, toupper(level), message))
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
}

read_gene_list <- function(file_path) {
  df <- fread(file_path, header = TRUE, data.table = FALSE)
  genes <- unique(as.character(df[[1]]))
  genes <- genes[!is.na(genes) & genes != ""]
  genes
}

write_csv_utf8 <- function(df, file_path) {
  write.csv(df, file_path, row.names = FALSE, fileEncoding = "UTF-8")
}

read_expression_matrix <- function(file_path) {
  df <- fread(file_path, header = TRUE, data.table = FALSE)
  rownames(df) <- df[[1]]
  df <- df[, -1, drop = FALSE]
  as.data.frame(t(df))
}

read_ranked_gene_list <- function(file_path, metric_col = "log2FC") {
  df <- fread(file_path, data.table = FALSE)
  possible_gene_cols <- c("gene", "Gene", "symbol", "Symbol", "SYMBOL")
  possible_metric_cols <- c("log2FC", "logFC", "FC", "score", "statistic")

  gene_col <- NULL
  for (col in possible_gene_cols) {
    if (col %in% colnames(df)) { gene_col <- col; break }
  }
  if (is.null(gene_col)) gene_col <- colnames(df)[1]

  metric_col_found <- NULL
  if (metric_col %in% colnames(df)) {
    metric_col_found <- metric_col
  } else {
    for (col in possible_metric_cols) {
      if (col %in% colnames(df)) { metric_col_found <- col; break }
    }
  }
  if (is.null(metric_col_found)) metric_col_found <- colnames(df)[ncol(df)]

  gene_list <- df[[metric_col_found]]
  names(gene_list) <- df[[gene_col]]
  gene_list <- gene_list[!is.na(gene_list) & is.finite(gene_list)]
  gene_list <- sort(gene_list, decreasing = TRUE)
  if (any(duplicated(names(gene_list)))) {
    gene_list <- gene_list[!duplicated(names(gene_list))]
  }
  log_msg("info", sprintf("Read ranked gene list: %d genes", length(gene_list)))
  gene_list
}

run_gsea <- function(gene_rank, term2gene, args, prefix = "GSEA") {
  gsea <- tryCatch({
    GSEA(
      geneList = gene_rank,
      exponent = args$exponent,
      minGSSize = args$min_gs_size,
      maxGSSize = args$max_gs_size,
      pvalueCutoff = args$pvalue_cutoff,
      pAdjustMethod = args$padj_method,
      TERM2GENE = term2gene,
      eps = args$eps
    )
  }, error = function(e) {
    log_msg("error", sprintf("GSEA failed: %s", e$message))
    return(NULL)
  })
  if (is.null(gsea)) return(NULL)
  list(gsea = gsea, gene_rank = gene_rank)
}

plot_gsea <- function(gsea_obj, args, output_dir, prefix = "GSEA") {
  if (is.null(gsea_obj) || nrow(gsea_obj@result) == 0) return(invisible(NULL))

  plot_dir <- file.path(output_dir, "plots")
  ensure_dir(plot_dir)

  show_n <- min(args$ridge_show_category, nrow(gsea_obj@result))
  top_pathways <- gsea_obj@result %>% arrange(desc(abs(NES))) %>% head(min(args$top_n, nrow(gsea_obj@result)))
  top_ids <- top_pathways$ID
  if (length(top_ids) == 0) return(invisible(NULL))

  sci_colors <- pal_npg(alpha = 0.9)(max(5, length(top_ids)))[seq_along(top_ids)]
  p_all <- gseaplot2(gsea_obj, geneSetID = top_ids, color = sci_colors,
                     pvalue_table = FALSE, subplots = 1:3, base_size = 9)

  pdf_file <- file.path(plot_dir, paste0(prefix, "_gsea.pdf"))
  grDevices::pdf(pdf_file, width = args$plot_width, height = args$plot_height)
  print(p_all)
  dev.off()
  log_msg("info", sprintf("Plots saved: %s", plot_dir))
}

# ---- Methods ----
run_prerank <- function(args) {
  log_msg("info", "========== Pre-rank GSEA ==========")
  gene_rank <- read_ranked_gene_list(args$ranked_gene_list, args$rank_metric)
  kegg_genesets <- clusterProfiler::read.gmt(args$gmt_file)

  result <- run_gsea(gene_rank, kegg_genesets, args, "PreRank")
  if (!is.null(result)) {
    tables_dir <- file.path(args$output, "tables")
    gsea_dir <- file.path(args$output, "gsea")
    ensure_dir(tables_dir)
    ensure_dir(gsea_dir)
    write_csv_utf8(result$gsea@result, file.path(gsea_dir, "prerank_gsea_results.csv"))
    plot_gsea(result$gsea, args, args$output, "PreRank")
    log_msg("info", sprintf("Done: %d pathways", nrow(result$gsea@result)))
  }
  result
}

run_correlation <- function(args) {
  log_msg("info", "========== Correlation GSEA ==========")
  key_genes <- read_gene_list(args$target_genes)
  exp <- read_expression_matrix(args$expression_matrix)
  kegg_genesets <- clusterProfiler::read.gmt(args$gmt_file)

  gsea_results <- list()
  for (gene in key_genes) {
    if (!gene %in% colnames(exp)) {
      log_msg("warn", sprintf("Gene not in matrix: %s", gene))
      next
    }
    log_msg("info", sprintf("Processing: %s", gene))
    cor_values <- cor(exp[, gene], exp, method = args$cor_method, use = "pairwise.complete.obs")
    cor_values <- as.numeric(cor_values[1, ])
    names(cor_values) <- colnames(exp)
    cor_values <- cor_values[names(cor_values) != gene]
    if (length(unique(cor_values)) < length(cor_values)) {
      cor_values <- cor_values + rnorm(length(cor_values), 0, 1e-8)
    }
    gene_rank <- sort(cor_values, decreasing = TRUE)

    result <- run_gsea(gene_rank, kegg_genesets, args, gene)
    if (!is.null(result)) {
      gsea_results[[gene]] <- result$gsea
      save_gene_outputs(gene, result$gsea, args$output)
      plot_gsea(result$gsea, args, args$output, gene)
    }
  }
  log_msg("info", sprintf("Done: %d genes analyzed", length(gsea_results)))
  gsea_results
}

run_median_split <- function(args) {
  log_msg("info", "========== Median Split GSEA ==========")
  key_genes <- read_gene_list(args$target_genes)
  exp <- read_expression_matrix(args$expression_matrix)
  kegg_genesets <- clusterProfiler::read.gmt(args$gmt_file)

  gsea_results <- list()
  for (gene in key_genes) {
    if (!gene %in% colnames(exp)) {
      log_msg("warn", sprintf("Gene not in matrix: %s", gene))
      next
    }
    log_msg("info", sprintf("Processing: %s", gene))

    gene_expression <- exp[, gene]
    median_expr <- median(gene_expression)
    high_group <- gene_expression > median_expr
    low_group <- gene_expression <= median_expr

    if (sum(high_group) < 3 || sum(low_group) < 3) {
      log_msg("warn", sprintf("Insufficient samples: %s", gene))
      next
    }

    exp_t <- t(exp)
    log2FC <- rowMeans(exp_t[, high_group]) - rowMeans(exp_t[, low_group])
    log2FC <- log2FC[!is.na(log2FC) & is.finite(log2FC)]
    if (any(duplicated(names(log2FC)))) {
      log2FC <- log2FC[!duplicated(names(log2FC))]
    }
    gene_rank <- sort(log2FC, decreasing = TRUE)

    result <- run_gsea(gene_rank, kegg_genesets, args, gene)
    if (!is.null(result)) {
      gsea_results[[gene]] <- result$gsea
      save_gene_outputs(gene, result$gsea, args$output)
      plot_gsea(result$gsea, args, args$output, gene)
    }
  }
  log_msg("info", sprintf("Done: %d genes analyzed", length(gsea_results)))
  gsea_results
}

save_gene_outputs <- function(gene, gsea_obj, output_dir) {
  tables_dir <- file.path(output_dir, "tables")
  gsea_dir <- file.path(output_dir, "gsea")
  ensure_dir(tables_dir)
  ensure_dir(gsea_dir)
  write_csv_utf8(gsea_obj@result, file.path(gsea_dir, paste0(gene, "_gsea_results.csv")))
}

# ---- Args ----
option_list <- list(
  make_option(c("-m", "--method"), type = "character", default = "prerank",
              help = "GSEA method: prerank, correlation, median_split [default: %default]"),
  make_option(c("-g", "--gmt-file"), type = "character",
              help = "GMT gene set file"),
  make_option(c("-e", "--expression-matrix"), type = "character",
              help = "Expression matrix (for correlation/median_split)"),
  make_option(c("-t", "--target-genes"), type = "character",
              help = "Target gene list (for correlation/median_split)"),
  make_option(c("-r", "--ranked-gene-list"), type = "character",
              help = "Ranked gene list (for prerank)"),
  make_option(c("-o", "--output"), type = "character", default = "results_gsea",
              help = "Output directory [default: %default]"),
  make_option(c("--min-gs-size"), type = "integer", default = 10,
              help = "Min gene set size [default: %default]"),
  make_option(c("--max-gs-size"), type = "integer", default = 500,
              help = "Max gene set size [default: %default]"),
  make_option(c("--pvalue-cutoff"), type = "double", default = 0.05,
              help = "P-value cutoff [default: %default]"),
  make_option(c("--padj-method"), type = "character", default = "BH",
              help = "P-value adjustment method [default: %default]"),
  make_option(c("--exponent"), type = "integer", default = 1,
              help = "GSEA exponent [default: %default]"),
  make_option(c("--eps"), type = "double", default = 0,
              help = "GSEA eps [default: %default]"),
  make_option(c("--rank-metric"), type = "character", default = "log2FC",
              help = "Ranking metric column (for prerank) [default: %default]"),
  make_option(c("--cor-method"), type = "character", default = "spearman",
              help = "Correlation method (for correlation) [default: %default]"),
  make_option(c("--top-n"), type = "integer", default = 5,
              help = "Top pathways to plot [default: %default]"),
  make_option(c("--ridge-show-category"), type = "integer", default = 10,
              help = "Categories for ridge plot [default: %default]"),
  make_option(c("--plot-width"), type = "integer", default = 12,
              help = "Plot width [default: %default]"),
  make_option(c("--plot-height"), type = "integer", default = 8,
              help = "Plot height [default: %default]")
)

opt_parser <- OptionParser(option_list = option_list)
args <- parse_args(opt_parser)

# Validate
if (is.null(args$gmt.file)) {
  print_help(opt_parser)
  stop("\n--gmt-file is required")
}

if (args$method == "prerank" && is.null(args$ranked.gene.list)) {
  stop("--ranked-gene-list is required for prerank method")
}

if (args$method %in% c("correlation", "median_split")) {
  if (is.null(args$expression.matrix)) {
    stop(sprintf("--expression-matrix is required for %s method", args$method))
  }
  if (is.null(args$target.genes)) {
    stop(sprintf("--target-genes is required for %s method", args$method))
  }
}

# Create output dir
ensure_dir(args$output)

# Print info
cat("============================================================\n")
cat("                   GSEA Analysis\n")
cat("============================================================\n")
cat(sprintf("Method: %s\n", args$method))
cat(sprintf("GMT file: %s\n", args$gmt.file))
cat(sprintf("Output: %s\n", args$output))
cat(sprintf("P-value cutoff: %s\n", args$pvalue.cutoff))
cat("============================================================\n\n")

# Run
start_time <- Sys.time()

result <- tryCatch({
  if (args$method == "prerank") {
    run_prerank(args)
  } else if (args$method == "correlation") {
    run_correlation(args)
  } else if (args$method == "median_split") {
    run_median_split(args)
  } else {
    stop(sprintf("Unknown method: %s", args$method))
  }
}, error = function(e) {
  log_msg("error", e$message)
  quit(status = 1)
})

duration <- difftime(Sys.time(), start_time, units = "mins")
cat("\n============================================================\n")
cat(sprintf("Done! Runtime: %.2f minutes\n", as.numeric(duration)))
cat(sprintf("Output: %s\n", normalizePath(args$output)))
cat("============================================================\n")

# ==============================================================================
# unified_enrichment/scripts/utils.R
# 通用工具函数模块
# 合并自 GO_KEGG/functions.R 和 单基因富集分析/scripts/utils.R
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

# ---- 空值合并操作符 ----

#' 空值合并操作符
#' @param x 第一个值
#' @param y 默认值
#' @return 如果 x 为 NULL 则返回 y，否则返回 x
`%||%` <- function(x, y) if (is.null(x)) y else x

# ---- 日志环境（避免全局变量污染）----

.log_env <- new.env(parent = emptyenv())

# ---- 目录操作 ----

#' 确保目录存在
#' @param path 目录路径
ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
}

#' 创建多个目录
#' @param paths 路径列表
prepare_dirs <- function(paths) {
  for (p in paths) {
    ensure_dir(p)
  }
}

# ---- 日志系统 ----

#' 初始化日志
#' @param log_file 日志文件路径
init_logger <- function(log_file) {
  ensure_dir(dirname(log_file))
  assign("log_file", log_file, envir = .log_env)
}

#' 记录日志消息
#' @param level 日志级别 (info/warn/error)
#' @param msg 消息内容
log_msg <- function(level, msg) {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- sprintf("[%s] [%s] %s", ts, toupper(level), msg)
  cat(line, "\n")
  if (exists("log_file", envir = .log_env)) {
    cat(line, "\n", file = get("log_file", envir = .log_env), append = TRUE)
  }
}

#' 安全执行函数
#' @param expr 表达式
#' @param step_name 步骤名称
#' @return 表达式结果，失败时返回 NULL
safe_run <- function(expr, step_name = "操作") {
  start_time <- Sys.time()
  log_msg("info", sprintf("开始执行: %s", step_name))
  tryCatch({
    result <- eval(expr)
    end_time <- Sys.time()
    duration <- difftime(end_time, start_time, units = "secs")
    log_msg("info", sprintf("完成: %s (耗时: %.2f秒)", step_name, duration))
    result
  }, error = function(e) {
    log_msg("error", sprintf("%s失败: %s", step_name, e$message))
    NULL  # 返回 NULL 而不是停止
  })
}

# ---- 配置访问辅助函数 ----

#' 获取配置值（支持默认值）
#' @param cfg 配置列表
#' @param path 配置路径（用点分隔，如 "analysis.go_kegg.enabled"）
#' @param default 默认值
#' @return 配置值
get_cfg <- function(cfg, path, default = NULL) {
  keys <- strsplit(path, "\\.")[[1]]
  value <- cfg
  for (key in keys) {
    if (is.null(value) || !key %in% names(value)) {
      return(default)
    }
    value <- value[[key]]
  }
  if (is.null(value)) default else value
}

# ---- 命令行参数解析 ----

#' 解析命令行参数
#' @return 配置列表
parse_command_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  # 初始化默认值
  cfg <- list(
    # 脚本目录（由 shell 传入）
    script_dir = NULL,

    # 必需参数
    work_dir = NULL,
    species = NULL,
    gene_list = NULL,
    expression_matrix = NULL,
    target_genes = NULL,
    ranked_list = NULL,
    gmt = NULL,

    # 分析开关
    analysis_mode = "full",  # full: 完整流程, gsea-only: 仅GSEA
    go_kegg = TRUE,
    gsea = TRUE,
    gsea_method = "prerank",

    # GO/KEGG参数
    go_ont = "ALL",
    go_pvalue = 0.05,
    go_qvalue = 0.2,
    go_padj = "BH",

    # GSEA参数
    gsea_correlation = "spearman",
    gsea_rank_metric = "log2FC",
    gsea_min_size = 10,
    gsea_max_size = 500,
    gsea_pvalue = 0.05,
    gsea_padj = "BH",
    gsea_exponent = 1,
    gsea_eps = 0,
    gsea_top_n = 5,

    # QC参数
    qc_variance = 0.01,
    qc_na_ratio = 0.5,
    qc_gmt_match = 30,

    # 绘图参数
    plot_width = 10,
    plot_height = 8,
    plot_top_n = 5,
    plot_color = "viridis",

    # 其他参数
    output_dir = "results",
    log_dir = "logs",
    seed = 123,

    # Quarto 报告参数
    quarto_report = FALSE,
    quarto_format = "docx",
    quarto_author = "LPF"
  )

  # 解析参数
  i <- 1
  while (i <= length(args)) {
    arg <- args[i]

    if (arg == "--args") {
      i <- i + 1
      next
    }

    # 脚本目录
    if (arg == "--script-dir" && i + 1 <= length(args)) {
      cfg$script_dir <- args[i + 1]
      i <- i + 2
      next
    }

    # 必需参数
    if (arg == "--work-dir" && i + 1 <= length(args)) {
      cfg$work_dir <- args[i + 1]
      i <- i + 2
    } else if (arg == "--species" && i + 1 <= length(args)) {
      cfg$species <- args[i + 1]
      i <- i + 2
    } else if (arg == "--gene-list" && i + 1 <= length(args)) {
      cfg$gene_list <- args[i + 1]
      i <- i + 2
    } else if (arg == "--expr-matrix" && i + 1 <= length(args)) {
      cfg$expression_matrix <- args[i + 1]
      i <- i + 2
    } else if (arg == "--target-genes" && i + 1 <= length(args)) {
      cfg$target_genes <- args[i + 1]
      i <- i + 2
    } else if (arg == "--ranked-list" && i + 1 <= length(args)) {
      cfg$ranked_list <- args[i + 1]
      i <- i + 2
    } else if (arg == "--gmt" && i + 1 <= length(args)) {
      cfg$gmt <- args[i + 1]
      i <- i + 2

    # 分析开关
    } else if (arg == "--go-kegg" && i + 1 <= length(args)) {
      cfg$go_kegg <- tolower(args[i + 1]) == "true"
      i <- i + 2
    } else if (arg == "--gsea" && i + 1 <= length(args)) {
      cfg$gsea <- tolower(args[i + 1]) == "true"
      i <- i + 2
    } else if (arg == "--gsea-method" && i + 1 <= length(args)) {
      cfg$gsea_method <- args[i + 1]
      i <- i + 2
    } else if (arg == "--analysis-mode" && i + 1 <= length(args)) {
      cfg$analysis_mode <- args[i + 1]
      i <- i + 2

    # GO/KEGG参数
    } else if (arg == "--go-ont" && i + 1 <= length(args)) {
      cfg$go_ont <- args[i + 1]
      i <- i + 2
    } else if (arg == "--go-pvalue" && i + 1 <= length(args)) {
      cfg$go_pvalue <- as.numeric(args[i + 1])
      i <- i + 2
    } else if (arg == "--go-qvalue" && i + 1 <= length(args)) {
      cfg$go_qvalue <- as.numeric(args[i + 1])
      i <- i + 2
    } else if (arg == "--go-padj" && i + 1 <= length(args)) {
      cfg$go_padj <- args[i + 1]
      i <- i + 2

    # GSEA参数
    } else if (arg == "--gsea-correlation" && i + 1 <= length(args)) {
      cfg$gsea_correlation <- args[i + 1]
      i <- i + 2
    } else if (arg == "--gsea-rank-metric" && i + 1 <= length(args)) {
      cfg$gsea_rank_metric <- args[i + 1]
      i <- i + 2
    } else if (arg == "--gsea-min-size" && i + 1 <= length(args)) {
      cfg$gsea_min_size <- as.integer(args[i + 1])
      i <- i + 2
    } else if (arg == "--gsea-max-size" && i + 1 <= length(args)) {
      cfg$gsea_max_size <- as.integer(args[i + 1])
      i <- i + 2
    } else if (arg == "--gsea-pvalue" && i + 1 <= length(args)) {
      cfg$gsea_pvalue <- as.numeric(args[i + 1])
      i <- i + 2
    } else if (arg == "--gsea-padj" && i + 1 <= length(args)) {
      cfg$gsea_padj <- args[i + 1]
      i <- i + 2
    } else if (arg == "--gsea-exponent" && i + 1 <= length(args)) {
      cfg$gsea_exponent <- as.numeric(args[i + 1])
      i <- i + 2
    } else if (arg == "--gsea-eps" && i + 1 <= length(args)) {
      cfg$gsea_eps <- as.numeric(args[i + 1])
      i <- i + 2
    } else if (arg == "--gsea-top-n" && i + 1 <= length(args)) {
      cfg$gsea_top_n <- as.integer(args[i + 1])
      i <- i + 2

    # QC参数
    } else if (arg == "--qc-variance" && i + 1 <= length(args)) {
      cfg$qc_variance <- as.numeric(args[i + 1])
      i <- i + 2
    } else if (arg == "--qc-na-ratio" && i + 1 <= length(args)) {
      cfg$qc_na_ratio <- as.numeric(args[i + 1])
      i <- i + 2
    } else if (arg == "--qc-gmt-match" && i + 1 <= length(args)) {
      cfg$qc_gmt_match <- as.integer(args[i + 1])
      i <- i + 2

    # 绘图参数
    } else if (arg == "--plot-width" && i + 1 <= length(args)) {
      cfg$plot_width <- as.numeric(args[i + 1])
      i <- i + 2
    } else if (arg == "--plot-height" && i + 1 <= length(args)) {
      cfg$plot_height <- as.numeric(args[i + 1])
      i <- i + 2
    } else if (arg == "--plot-top-n" && i + 1 <= length(args)) {
      cfg$plot_top_n <- as.integer(args[i + 1])
      i <- i + 2
    } else if (arg == "--plot-color" && i + 1 <= length(args)) {
      cfg$plot_color <- args[i + 1]
      i <- i + 2

    # 其他参数
    } else if (arg == "--output-dir" && i + 1 <= length(args)) {
      cfg$output_dir <- args[i + 1]
      i <- i + 2
    } else if (arg == "--log-dir" && i + 1 <= length(args)) {
      cfg$log_dir <- args[i + 1]
      i <- i + 2
    } else if (arg == "--seed" && i + 1 <= length(args)) {
      cfg$seed <- as.integer(args[i + 1])
      i <- i + 2

    # Quarto 报告参数
    } else if (arg == "--quarto-report") {
      cfg$quarto_report <- TRUE
      i <- i + 1
    } else if (arg == "--quarto-format" && i + 1 <= length(args)) {
      cfg$quarto_format <- args[i + 1]
      i <- i + 2
    } else if (arg == "--quarto-author" && i + 1 <= length(args)) {
      cfg$quarto_author <- args[i + 1]
      i <- i + 2
    } else {
      i <- i + 1
    }
  }

  cfg
}

#' 验证必需参数
#' @param cfg 配置列表
validate_config <- function(cfg) {
  missing <- c()

  if (is.null(cfg$work_dir)) missing <- c(missing, "--work-dir")
  if (is.null(cfg$species)) missing <- c(missing, "--species")

  # 验证 analysis_mode
  if (!is.null(cfg$analysis_mode) && !cfg$analysis_mode %in% c("full", "gsea-only")) {
    stop(sprintf("不支持的分析模式: %s。支持: full, gsea-only", cfg$analysis_mode))
  }

  # 根据模式验证不同参数
  if (is.null(cfg$analysis_mode) || cfg$analysis_mode == "full") {
    # 完整模式：需要所有参数
    if (is.null(cfg$gene_list)) missing <- c(missing, "--gene-list")
    if (is.null(cfg$expression_matrix)) missing <- c(missing, "--expr-matrix")
    if (is.null(cfg$target_genes)) missing <- c(missing, "--target-genes")
  }

  # GSEA 模式需要
  if (isTRUE(cfg$gsea)) {
    if (is.null(cfg$ranked_list)) missing <- c(missing, "--ranked-list")
    if (is.null(cfg$gmt)) missing <- c(missing, "--gmt")
  }

  if (length(missing) > 0) {
    cat("错误: 缺少必需参数:\n")
    for (param in missing) {
      cat(sprintf("  - %s\n", param))
    }
    cat(sprintf("\n提示: 使用 --analysis-mode gsea-only 可仅运行GSEA分析\n"))
    stop("缺少必需参数")
  }

  # 验证物种
  if (!cfg$species %in% c("hsa", "mmu", "rno")) {
    stop(sprintf("不支持的物种: %s。支持: hsa(人类), mmu(小鼠), rno(大鼠)", cfg$species))
  }
}

#' 从配置构建路径列表
#' 路径处理：如果已经是绝对路径则直接使用，否则拼接 project_root
#' @param project_root 项目根目录
#' @param path 路径
#' @return 处理后的路径
resolve_path <- function(project_root, path) {
  if (is.null(path) || nchar(path) == 0) {
    return(path)
  }
  # 如果是绝对路径直接返回
  if (substr(path, 1, 1) == "/") {
    return(path)
  }
  # 否则拼接 project_root
  return(file.path(project_root, path))
}

#' @param cfg 配置列表
#' @return 路径列表
build_paths <- function(cfg) {
  project_root <- cfg$work_dir

  # 输出路径处理 - 支持绝对路径
  resolve_output_path <- function(path) {
    if (is.null(path) || nchar(path) == 0) {
      return(path)
    }
    if (substr(path, 1, 1) == "/") {
      return(path)
    }
    return(file.path(project_root, path))
  }

  output_base <- resolve_output_path(cfg$output_dir)

  paths <- list(
    # 输入 - 使用 resolve_path 处理
    gene_list = resolve_path(project_root, cfg$gene_list),
    expression_matrix = resolve_path(project_root, cfg$expression_matrix),
    target_genes = resolve_path(project_root, cfg$target_genes),
    ranked_gene_list = resolve_path(project_root, cfg$ranked_list),
    gmt_file = resolve_path(project_root, cfg$gmt),
    # 输出 - logs放在模块根目录的logs/下
    results_dir = output_base,
    intermediate_dir = file.path(output_base, "intermediate"),
    qc_dir = file.path(output_base, "qc"),
    logs_dir = file.path(MODULE_ROOT, "logs"),
    tables_dir = file.path(output_base, "tables"),
    plots_dir = file.path(output_base, "plots"),
    go_dir = file.path(output_base, "go"),
    kegg_dir = file.path(output_base, "kegg"),
    gsea_dir = file.path(output_base, "gsea")
  )
  paths
}

# ---- 基因数据处理 ----

#' 读取基因列表（支持自动检测表头）
#' @param gene_file 基因文件路径
#' @return 基因向量
read_gene_list <- function(gene_file) {
  if (!file.exists(gene_file)) {
    stop(sprintf("基因文件不存在: %s", gene_file))
  }

  line1 <- trimws(readLines(gene_file, n = 1, warn = FALSE, encoding = "UTF-8"))
  first_field <- trimws(strsplit(line1, ",", fixed = TRUE)[[1]][1])
  first_field <- gsub('^"|"$', "", first_field)
  has_header <- tolower(first_field) %in% c("gene", "genes", "symbol", "symbols")

  genes_df <- if (has_header) {
    read.csv(gene_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8")
  } else {
    read.csv(gene_file, header = FALSE, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8")
  }

  gene_col_candidates <- c("gene", "genes", "symbol", "symbols")
  gene_col_idx <- which(tolower(colnames(genes_df)) %in% gene_col_candidates)

  key_genes <- if (length(gene_col_idx) >= 1) {
    genes_df[[gene_col_idx[1]]]
  } else if (ncol(genes_df) == 1) {
    genes_df[[1]]
  } else {
    stop("基因文件有多列但未找到gene/symbol列")
  }

  key_genes <- unique(trimws(as.character(key_genes)))
  key_genes <- key_genes[!is.na(key_genes) & key_genes != ""]

  if (length(key_genes) == 0) {
    stop("基因文件中未找到有效基因名")
  }
  key_genes
}

#' 推断基因ID类型
#' @param genes 基因向量
#' @return 包含species和id_type的列表
infer_expr_id <- function(genes) {
  genes <- genes[!is.na(genes) & genes != ""]
  if (length(genes) == 0) {
    return(list(species = "unknown", id_type = "unknown"))
  }

  n <- min(length(genes), 2000)
  genes <- sample(genes, n)

  frac_ensg <- mean(grepl("^ENSG[0-9]+", genes))
  frac_ensmusg <- mean(grepl("^ENSMUSG[0-9]+", genes))
  frac_ensembl_gene <- mean(grepl("^ENS[A-Z]{0,6}G[0-9]+", genes))
  frac_allcaps_symbol <- mean(grepl("^[A-Z0-9\\-\\.]+$", genes) & grepl("[A-Z]", genes))
  frac_mouse_symbol <- mean(grepl("^[A-Z][a-z0-9]+$", genes))

  species <- "unknown"
  if (frac_ensmusg >= 0.2) species <- "mmu"
  if (frac_ensg >= 0.2) species <- "hsa"

  id_type <- "unknown"
  if (frac_ensembl_gene >= 0.2) id_type <- "ensembl"
  if (frac_allcaps_symbol >= 0.4) id_type <- "symbols"
  if (id_type == "unknown" && frac_mouse_symbol >= 0.4) id_type <- "symbols"

  list(species = species, id_type = id_type)
}

#' 推断GMT文件元数据
#' @param gmt_file GMT文件路径
#' @return 包含species和id_type的列表
infer_gmt_meta <- function(gmt_file) {
  gmt_species <- "unknown"
  if (grepl("\\.Hs\\.|\\bHs\\b", gmt_file, ignore.case = FALSE)) gmt_species <- "hsa"
  if (grepl("\\.Mm\\.|\\bMm\\b", gmt_file, ignore.case = FALSE)) gmt_species <- "mmu"

  gmt_id_type <- "unknown"
  if (grepl("symbols", gmt_file, ignore.case = TRUE)) gmt_id_type <- "symbols"
  if (grepl("ensembl", gmt_file, ignore.case = TRUE)) gmt_id_type <- "ensembl"

  list(species = gmt_species, id_type = gmt_id_type)
}

#' 验证物种和ID一致性
#' @param expr_info 表达矩阵信息
#' @param gmt_info GMT信息
validate_species_id <- function(expr_info, gmt_info) {
  if (expr_info$species == "unknown") {
    stop("无法推断表达矩阵物种，请在配置中设置 project.species")
  }

  if (gmt_info$species != "unknown" && expr_info$species != gmt_info$species) {
    stop(sprintf("物种不一致：表达矩阵=%s, GMT=%s", expr_info$species, gmt_info$species))
  }

  if (gmt_info$id_type == "symbols" && expr_info$id_type == "ensembl") {
    stop("ID体系不一致：表达矩阵为Ensembl，但GMT为symbols")
  }

  if (gmt_info$id_type == "ensembl" && expr_info$id_type == "symbols") {
    stop("ID体系不一致：表达矩阵为symbols，但GMT为Ensembl")
  }
}

# ---- 文件读写 ----

#' 写入UTF-8编码的CSV
#' @param df 数据框
#' @param file 文件路径
write_csv_utf8 <- function(df, file) {
  write.csv(df, file = file, row.names = FALSE, fileEncoding = "UTF-8")
}

#' 检测输入文件类型
#' @param input_path 输入文件路径
#' @return 文件类型 (gene_list / expression_matrix)
detect_input_type <- function(input_path) {
  df <- fread(input_path, nrows = 5)

  # 检测是否为表达矩阵（列数多且第二列起为数值型）
  if (ncol(df) > 10) {
    # 检查是否大部分列是数值
    numeric_cols <- sum(sapply(df[, -1, with = FALSE], function(x) {
      is.numeric(x) || all(!is.na(suppressWarnings(as.numeric(as.character(x)))))
    }))
    if (numeric_cols > ncol(df) * 0.8) {
      return("expression_matrix")
    }
  }

  # 检测是否为基因列表（单列或少量列）
  if (ncol(df) <= 3) {
    return("gene_list")
  }

  stop("无法识别输入文件类型")
}

# ---- 绘图工具 ----

#' 保存绘图到PDF和PNG
#' @param filename 文件名
#' @param plot 绘图对象
#' @param outdir 输出目录
#' @param width 宽度
#' @param height 高度
#' @param both 是否同时保存PDF和PNG
save_plot <- function(filename, plot, outdir = ".",
                      width = 7, height = 7, both = FALSE,
                      bg = 'white', family = "Arial",
                      units = "in", res = 600) {

  is_gg <- inherits(plot, "ggplot")

  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }

  base_name <- tools::file_path_sans_ext(filename)
  file_path_pdf <- file.path(outdir, paste0(base_name, ".pdf"))
  file_path_png <- file.path(outdir, paste0(base_name, ".png"))

  execute_plot <- function() {
    if (is_gg) {
      print(plot)
    } else if (is.function(plot)) {
      plot()
    } else {
      print(plot)
    }
  }

  if (both) {
    # PDF
    grDevices::cairo_pdf(file = file_path_pdf, width = width, height = height, bg = bg, family = family)
    execute_plot()
    dev.off()

    # PNG
    png(file = file_path_png, width = width, height = height,
        bg = bg, units = units, res = res, family = family)
    execute_plot()
    dev.off()

    log_msg("info", sprintf("已保存PDF和PNG: %s, %s", file_path_pdf, file_path_png))
  } else {
    ext <- tools::file_ext(filename)
    full_path <- file.path(outdir, filename)

    if (ext == "pdf") {
      grDevices::cairo_pdf(file = full_path, width = width, height = height, bg = bg, family = family)
      execute_plot()
      dev.off()
    } else if (ext == "png") {
      png(file = full_path, width = width, height = height,
          bg = bg, units = units, res = res, family = family)
      execute_plot()
      dev.off()
    } else {
      stop("不支持的文件格式，请使用 .pdf 或 .png")
    }

    log_msg("info", sprintf("文件已保存: %s", full_path))
  }
}

# ---- 物种数据库 ----

#' 获取物种对应的OrgDb
#' @param species 物种代码 (hsa/mmu/rno)
#' @return OrgDb对象
get_org_db <- function(species) {
  if (species == "hsa") {
    if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
      BiocManager::install("org.Hs.eg.db", ask = FALSE)
    }
    library(org.Hs.eg.db)
    return(org.Hs.eg.db)
  } else if (species == "mmu") {
    if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
      BiocManager::install("org.Mm.eg.db", ask = FALSE)
    }
    library(org.Mm.eg.db)
    return(org.Mm.eg.db)
  } else if (species == "rno") {
    if (!requireNamespace("org.Rn.eg.db", quietly = TRUE)) {
      BiocManager::install("org.Rn.eg.db", ask = FALSE)
    }
    library(org.Rn.eg.db)
    return(org.Rn.eg.db)
  } else {
    stop(sprintf("不支持的物种: %s。支持: hsa(人类), mmu(小鼠), rno(大鼠)", species))
  }
}

#' 获取物种名称
#' @param species 物种代码
#' @return 物种中文名称
get_species_name <- function(species) {
  species_names <- list(
    hsa = "人类",
    mmu = "小鼠",
    rno = "大鼠"
  )
  species_names[[species]] %||% species
}

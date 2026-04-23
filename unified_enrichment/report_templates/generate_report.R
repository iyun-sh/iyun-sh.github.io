# ==============================================================================
# unified_enrichment 模块 Quarto 报告生成器
# ==============================================================================
# 功能：生成富集分析的 Quarto 报告
# 默认输出格式：docx (Word)
# ==============================================================================

#' 检查 Quarto 是否可用
check_quarto <- function() {
  quarto_path <- Sys.which("quarto")
  if (nchar(quarto_path) == 0) {
    warning("Quarto 未安装，无法生成报告。请访问 https://quarto.org 安装。")
    return(FALSE)
  }
  return(TRUE)
}

#' 生成富集分析 Quarto 报告
#' @param results_dir 富集分析结果目录
#' @param output_format 输出格式: "docx", "html", "pdf"
#' @param project_name 项目名称
#' @param author 作者
#' @param species 物种
#' @param analysis_type 分析类型
#' @return 报告文件路径，失败返回 NULL
generate_enrichment_quarto_report <- function(
    results_dir,
    output_format = "docx",
    project_name = "Enrichment Analysis",
    author = "IKL",
    species = "hsapiens",
    analysis_type = "combined") {

  # 检查 Quarto
  if (!check_quarto()) {
    return(NULL)
  }

  # 检查结果目录
  if (!dir.exists(results_dir)) {
    stop(sprintf("结果目录不存在: %s", results_dir))
  }
  results_dir <- normalizePath(results_dir)

  # 检测分析类型
  go_dir <- file.path(results_dir, "go")
  kegg_dir <- file.path(results_dir, "kegg")
  gsea_dir <- file.path(results_dir, "gsea")

  go_enabled <- dir.exists(go_dir) && length(list.files(go_dir, "\\.csv$")) > 0
  kegg_enabled <- dir.exists(kegg_dir) && length(list.files(kegg_dir, "\\.csv$")) > 0
  gsea_enabled <- dir.exists(gsea_dir) && length(list.files(gsea_dir, "\\.csv$")) > 0

  # 统计结果数量
  n_go_terms <- 0
  n_kegg_terms <- 0
  n_gsea_terms <- 0

  if (go_enabled) {
    go_bp_file <- file.path(go_dir, "GO_BP_results.csv")
    if (file.exists(go_bp_file)) {
      go_bp <- read.csv(go_bp_file, stringsAsFactors = FALSE)
      n_go_terms <- nrow(go_bp)
    }
  }

  if (kegg_enabled) {
    kegg_file <- file.path(kegg_dir, "KEGG_results.csv")
    if (file.exists(kegg_file)) {
      kegg_results <- read.csv(kegg_file, stringsAsFactors = FALSE)
      n_kegg_terms <- nrow(kegg_results)
    }
  }

  if (gsea_enabled) {
    gsea_file <- file.path(gsea_dir, "GSEA_results.csv")
    if (file.exists(gsea_file)) {
      gsea_results <- read.csv(gsea_file, stringsAsFactors = FALSE)
      n_gsea_terms <- nrow(gsea_results)
    }
  }

  # 查找模板
  template_path <- file.path(results_dir, "..", "report_templates", "report_template.qmd")
  if (!file.exists(template_path)) {
    template_path <- file.path(dirname(results_dir), "report_templates", "report_template.qmd")
  }
  if (!file.exists(template_path)) {
    template_path <- "report_templates/report_template.qmd"
  }

  if (!file.exists(template_path)) {
    stop(sprintf("报告模板不存在: %s", template_path))
  }

  # 创建参数
  params <- list(
    project_name = project_name,
    author = author,
    date = format(Sys.time(), "%Y-%m-%d"),
    results_dir = results_dir,
    species = species,
    analysis_type = analysis_type,
    go_enabled = go_enabled,
    kegg_enabled = kegg_enabled,
    gsea_enabled = gsea_enabled,
    n_go_terms = n_go_terms,
    n_kegg_terms = n_kegg_terms,
    n_gsea_terms = n_gsea_terms
  )

  params_file <- tempfile(fileext = ".json")
  jsonlite::write_json(params, params_file, auto_unbox = TRUE, pretty = TRUE)

  # 构建输出文件名
  safe_name <- gsub("[^A-Za-z0-9_-]", "_", project_name)
  output_file <- paste0(safe_name, "_Enrichment_Report")

  # 构建 quarto 命令
  quarto_cmd <- sprintf(
    'quarto render "%s" --to %s --output-dir "%s" --output "%s" --metadata-file "%s"',
    template_path,
    output_format,
    results_dir,
    output_file,
    params_file
  )

  # 执行命令
  exit_code <- system(quarto_cmd)

  # 清理临时文件
  unlink(params_file)

  if (exit_code != 0) {
    warning(sprintf("Quarto 渲染失败，退出码: %d", exit_code))
    return(NULL)
  }

  # 返回报告路径
  ext <- switch(output_format,
    docx = ".docx",
    html = ".html",
    pdf = ".pdf",
    ".docx"
  )

  report_path <- file.path(results_dir, paste0(output_file, ext))

  if (file.exists(report_path)) {
    cat(sprintf("Quarto 报告已生成: %s\n", report_path))
    return(report_path)
  } else {
    warning(sprintf("报告文件未找到: %s", report_path))
    return(NULL)
  }
}

# 命令行接口
if (!interactive() && sys.nframe() == 0) {
  suppressPackageStartupMessages({
    library(optparse)
    library(jsonlite)
  })

  option_list <- list(
    make_option(c("-r", "--results-dir"), type = "character",
                help = "富集分析结果目录 [必需]"),
    make_option(c("-f", "--format"), type = "character", default = "docx",
                help = "输出格式: docx, html, pdf [默认: docx]"),
    make_option(c("-p", "--project-name"), type = "character", default = "Enrichment Analysis",
                help = "项目名称"),
    make_option(c("-a", "--author"), type = "character", default = "IKL",
                help = "作者"),
    make_option(c("-s", "--species"), type = "character", default = "hsapiens",
                help = "物种"),
    make_option(c("-t", "--analysis-type"), type = "character", default = "combined",
                help = "分析类型")
  )

  opt <- parse_args(OptionParser(option_list = option_list))

  if (is.null(opt$results_dir)) {
    stop("必须指定 --results-dir 参数")
  }

  generate_enrichment_quarto_report(
    results_dir = opt$results_dir,
    output_format = opt$format,
    project_name = opt$project_name,
    author = opt$author,
    species = opt$species,
    analysis_type = opt$analysis_type
  )
}

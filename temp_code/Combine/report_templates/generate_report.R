# ============================================================================
# Combine 报告生成脚本
# ============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

option_list <- list(
  make_option(c("-d", "--results-dir"), dest = "results_dir", type = "character",
              help = "结果目录路径"),
  make_option(c("-f", "--format"), dest = "format", type = "character",
              default = "docx", help = "输出格式 [default: %default]"),
  make_option(c("-p", "--project-name"), dest = "project_name", type = "character",
              default = "project", help = "项目名称 [default: %default]"),
  make_option(c("-a", "--author"), dest = "author", type = "character",
              default = "生信团队", help = "报告作者 [default: %default]"),
  make_option(c("--output-name"), dest = "output_name", type = "character",
              default = "combat", help = "输出文件前缀名 [default: %default]")
)

parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser, convert_hyphens_to_underscores = TRUE)

if (is.null(opt$results_dir)) {
  print_help(parser); stop("缺少必需参数: --results-dir")
}
results_dir <- normalizePath(opt$results_dir, winslash = "/", mustWork = TRUE)

# ---- 自动检测参数 ----
info_file <- file.path(results_dir, "00.dataset_info.csv")
n_datasets <- 0
n_common_genes <- 0
if (file.exists(info_file)) {
  info <- fread(info_file)
  common_row <- info[info[[1]] == "共同基因", ]
  dataset_rows <- info[info[[1]] != "共同基因", ]
  n_datasets <- nrow(dataset_rows)
  if (nrow(common_row) > 0) n_common_genes <- common_row[[2]]
}

# ---- 复制模板到结果目录 ----
get_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0) {
    return(normalizePath(dirname(sub("^--file=", "", file_arg[1])), winslash = "/", mustWork = TRUE))
  }

  frame_ofile <- tryCatch(sys.frames()[[1]]$ofile, error = function(e) NULL)
  if (!is.null(frame_ofile) && nzchar(frame_ofile)) {
    return(normalizePath(dirname(frame_ofile), winslash = "/", mustWork = TRUE))
  }

  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

template_dir <- get_script_dir()

template_file <- file.path(template_dir, "report_template.qmd")
reference_doc <- file.path(template_dir, "reference.docx")
if (!file.exists(reference_doc)) stop(sprintf("缺少模块本地模板: %s", reference_doc))
if (!file.exists(template_file)) stop(sprintf("模板文件不存在: %s", template_file))

target_qmd <- file.path(results_dir, "report_template.qmd")
file.copy(template_file, target_qmd, overwrite = TRUE)
file.copy(reference_doc, file.path(results_dir, "reference.docx"), overwrite = TRUE)

# ---- 渲染报告 ----
cat("============================================================\n")
cat(sprintf(" 生成报告: %s\n", opt$project_name))
cat("============================================================\n")

params_file <- tempfile(fileext = ".yml")
on.exit(unlink(params_file), add = TRUE)
writeLines(c(
  sprintf("project_name: '%s'", gsub("'", "''", opt$project_name)),
  sprintf("author: '%s'", gsub("'", "''", opt$author)),
  sprintf("date: '%s'", format(Sys.time(), "%Y-%m-%d")),
  sprintf("results_dir: '%s'", results_dir),
  sprintf("n_datasets: %d", n_datasets),
  sprintf("n_common_genes: %d", n_common_genes),
  sprintf("output_name: '%s'", gsub("'", "''", opt$output_name))
), params_file)

cmd <- sprintf(
  "quarto render \"%s\" --to %s --execute-params \"%s\"",
  target_qmd, opt$format, params_file
)

cat("============================================================\n")
cat(sprintf(" 生成报告: %s\n", opt$project_name))
cat("============================================================\n")
cat(sprintf("执行: %s\n", cmd))
exit_code <- system(cmd)
if (exit_code != 0) stop("Quarto 渲染失败")

# ---- 重命名输出 ----
rendered_file <- file.path(results_dir, paste0("report_template.", opt$format))
target_name <- file.path(results_dir, sprintf("combine_report_%s.%s", opt$project_name, opt$format))
if (file.exists(rendered_file)) {
  file.rename(rendered_file, target_name)
  cat(sprintf("报告已生成: %s\n", target_name))
}

# ---- 清理临时文件 ----
unlink(target_qmd)
unlink(file.path(results_dir, "report_template_files"), recursive = TRUE)
unlink(file.path(results_dir, "reference.docx"))

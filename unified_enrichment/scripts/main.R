#!/usr/bin/env Rscript
# ==============================================================================
# unified_enrichment/scripts/main.R
# 统一富集分析主入口脚本
# 支持 GO/KEGG/GSEA 三种分析的灵活组合执行
# 所有配置通过命令行参数传递
# ==============================================================================

# ============================================
# 0. 时间戳初始化和目录结构
# ============================================
TIMESTAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")
SCRIPT_NAME <- "unified_enrichment"

# 模块根目录（硬编码）
MODULE_ROOT <- "/media/desk16/share/secure/unified_enrichment"
message("MODULE_ROOT = ", MODULE_ROOT)


# 自动创建目录结构
init_module_structure <- function(root_dir) {
  required_dirs <- c("logs", "report", "results")
  for (dir_name in required_dirs) {
    dir_path <- file.path(root_dir, dir_name)
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE)
    }
  }
}
init_module_structure(MODULE_ROOT)

# 设置环境
options(scipen = 5)
options(BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
options("repos" = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN"))

# 加载依赖（使用绝对路径）
source(file.path(MODULE_ROOT, "scripts/requirements.R"))
install_dependencies()
setwd(MODULE_ROOT)

# 加载自定义模块（使用绝对路径）
source(file.path(MODULE_ROOT, "scripts/utils.R"))
source(file.path(MODULE_ROOT, "scripts/qc.R"))
source(file.path(MODULE_ROOT, "scripts/go_kegg.R"))
source(file.path(MODULE_ROOT, "scripts/gsea.R"))
source(file.path(MODULE_ROOT, "scripts/report.R"))

save_qs2_snapshot <- function(object, output_path) {
  if (requireNamespace("qs2", quietly = TRUE)) {
    qs2::qs_save(object, output_path)
    return(invisible(output_path))
  }
  if (requireNamespace("qs", quietly = TRUE)) {
    qs::qsave(object, output_path)
    return(invisible(output_path))
  }
  warning("qs2/qs package not installed; skip qs2 snapshot: ", output_path)
  invisible(NULL)
}

# ---- 构建兼容配置结构 ----

#' 将命令行参数转换为兼容旧代码的配置结构
#' @param cfg 命令行参数配置
#' @return 兼容格式的配置列表
build_compatible_config <- function(cfg) {
  list(
    project = list(
      name = basename(cfg$work_dir),
      species = cfg$species,
      seed = cfg$seed,
      default_id_type = "symbols"
    ),
    analysis = list(
      go_kegg = list(
        enabled = cfg$go_kegg,
        go_ont = cfg$go_ont,
        pvalue_cutoff = cfg$go_pvalue,
        qvalue_cutoff = cfg$go_qvalue,
        padj_method = cfg$go_padj
      ),
      gsea = list(
        enabled = cfg$gsea,
        method = cfg$gsea_method,
        correlation_method = cfg$gsea_correlation,
        rank_metric = cfg$gsea_rank_metric,
        minGSSize = cfg$gsea_min_size,
        maxGSSize = cfg$gsea_max_size,
        pvalue_cutoff = cfg$gsea_pvalue,
        padj_method = cfg$gsea_padj,
        exponent = cfg$gsea_exponent,
        eps = cfg$gsea_eps,
        top_pathway_n = cfg$gsea_top_n
      )
    ),
    input = list(
      gene_list = cfg$gene_list,
      expression_matrix = cfg$expression_matrix,
      target_genes = cfg$target_genes,
      ranked_gene_list = cfg$ranked_list,
      gmt_file = cfg$gmt
    ),
    qc = list(
      variance_threshold = cfg$qc_variance,
      na_ratio_threshold = cfg$qc_na_ratio,
      min_gmt_match_rate = cfg$qc_gmt_match
    ),
    plot = list(
      width = cfg$plot_width,
      height = cfg$plot_height,
      top_n = cfg$plot_top_n,
      color_scheme = cfg$plot_color,
      ridge_show_category = cfg$gsea_top_n * 2,
      ridge_width = cfg$plot_width,
      ridge_height = cfg$plot_height
    ),
    output = list(
      results_dir = cfg$output_dir,
      intermediate_dir = cfg$output_dir,
      qc_dir = cfg$output_dir,
      logs_dir = cfg$output_dir,
      tables_dir = cfg$output_dir,
      plots_dir = cfg$output_dir,
      go_dir = cfg$output_dir,
      kegg_dir = cfg$output_dir,
      gsea_dir = cfg$output_dir
    )
  )
}

# ---- 配置访问辅助函数 ----

# ---- 主流程 ----

#' 主函数
main <- function() {
  # 记录开始时间
  script_start_time <- Sys.time()

  # 解析命令行参数
  cfg_raw <- parse_command_args()
  validate_config(cfg_raw)

  # 如果传入 script_dir，覆盖默认的 SCRIPT_DIR
  if (!is.null(cfg_raw$script_dir)) {
    SCRIPT_DIR <<- cfg_raw$script_dir
    MODULE_ROOT <<- dirname(SCRIPT_DIR)
  }

  # 构建兼容配置
  cfg <- build_compatible_config(cfg_raw)
  project_root <- cfg_raw$work_dir

  cat("============================================================\n")
  cat("        统一富集分析流程 (GO/KEGG/GSEA)\n")
  cat("============================================================\n\n")

  # 构建路径
  paths <- build_paths(cfg_raw)

  # 设置随机种子
  set.seed(cfg_raw$seed)

  # 初始化日志（标准化时间戳格式）
  log_dir <- file.path(MODULE_ROOT, "logs")
  log_file <- file.path(log_dir, paste0(SCRIPT_NAME, ".", TIMESTAMP, ".log"))
  init_logger(log_file)

  # 创建隐藏追踪文件
  result_dir <- file.path(MODULE_ROOT, "results", SCRIPT_NAME)
  if (!dir.exists(result_dir)) {
    dir.create(result_dir, recursive = TRUE)
  }
  timestamp_file <- file.path(result_dir, paste0(".", TIMESTAMP))
  file.create(timestamp_file)

  # 保存配置文件（使用调整后的值）
  # 检查是否 GSEA-only 模式
  is_gsea_only <- !is.null(cfg_raw$analysis_mode) && cfg_raw$analysis_mode == "gsea-only"
  go_kegg_value <- ifelse(is_gsea_only, "false", tolower(as.character(cfg_raw$go_kegg)))
  gsea_value <- "true"  # 总是启用 GSEA

  config_file <- file.path(result_dir, paste0(SCRIPT_NAME, ".config.ini"))
  config_lines <- c(
    "[Analysis]",
    paste0("timestamp=", TIMESTAMP),
    paste0("script=main.R"),
    "",
    "[Parameters]",
    paste0("species=", cfg_raw$species),
    paste0("go_kegg=", go_kegg_value),
    paste0("gsea=", gsea_value),
    paste0("gsea_method=", cfg_raw$gsea_method),
    paste0("go_ont=", cfg_raw$go_ont),
    paste0("go_pvalue=", cfg_raw$go_pvalue),
    paste0("go_qvalue=", cfg_raw$go_qvalue),
    paste0("gsea_pvalue=", cfg_raw$gsea_pvalue),
    paste0("seed=", cfg_raw$seed)
  )
  writeLines(config_lines, config_file)

  # 标准化日志头部
  cat("================================================================================\n")
  cat("统一富集分析日志 - GO/KEGG/GSEA\n")
  cat("================================================================================\n")
  cat(sprintf("时间戳: %s\n", TIMESTAMP))
  cat(sprintf("运行ID: %s\n", TIMESTAMP))
  cat(sprintf("脚本: main.R\n"))
  cat(sprintf("运行环境: R version %s.%s\n", R.version$major, R.version$minor))
  cat("\n")
  cat("----------------------------------------\n")
  cat("输入参数\n")
  cat("----------------------------------------\n")
  cat(sprintf("分析模式: %s\n", ifelse(is.null(cfg_raw$analysis_mode), "full", cfg_raw$analysis_mode)))
  cat(sprintf("物种: %s\n", cfg_raw$species))
  cat(sprintf("GO/KEGG分析: %s\n", ifelse(cfg_raw$go_kegg, "启用", "禁用")))
  cat(sprintf("GSEA分析: %s (%s)\n", ifelse(cfg_raw$gsea, "启用", "禁用"), cfg_raw$gsea_method))
  if (cfg_raw$analysis_mode != "gsea-only") {
    cat(sprintf("GO本体: %s\n", cfg_raw$go_ont))
    cat(sprintf("GO p值阈值: %s\n", cfg_raw$go_pvalue))
    cat(sprintf("GO q值阈值: %s\n", cfg_raw$go_qvalue))
    cat(sprintf("GO p值校正方法: %s\n", cfg_raw$go_padj))
  }
  if (cfg_raw$gsea) {
    cat(sprintf("GSEA方法: %s\n", cfg_raw$gsea_method))
    cat(sprintf("GSEA minGSSize: %s\n", cfg_raw$gsea_min_size))
    cat(sprintf("GSEA maxGSSize: %s\n", cfg_raw$gsea_max_size))
    cat(sprintf("GSEA p值阈值: %s\n", cfg_raw$gsea_pvalue))
    cat(sprintf("GSEA p值校正方法: %s\n", cfg_raw$gsea_padj))
  }
  cat("\n")
  cat("----------------------------------------\n")
  cat("输入文件\n")
  cat("----------------------------------------\n")
  if (cfg_raw$analysis_mode != "gsea-only") {
    cat(sprintf("基因列表: %s\n", cfg_raw$gene_list))
    cat(sprintf("表达矩阵: %s\n", cfg_raw$expression_matrix))
    cat(sprintf("目标基因: %s\n", cfg_raw$target_genes))
  }
  cat(sprintf("排序列表: %s\n", cfg_raw$ranked_list))
  cat(sprintf("GMT文件: %s\n", cfg_raw$gmt))
  cat(sprintf("输出目录: %s\n", cfg_raw$output_dir))
  cat("\n")

  log_msg("info", sprintf("工作目录: %s", project_root))
  log_msg("info", sprintf("分析模式: %s", ifelse(is.null(cfg_raw$analysis_mode), "full", cfg_raw$analysis_mode)))
  log_msg("info", sprintf("物种: %s", cfg_raw$species))
  log_msg("info", sprintf("GO/KEGG分析: %s", ifelse(cfg_raw$go_kegg, "启用", "禁用")))
  log_msg("info", sprintf("GSEA分析: %s (方法: %s)", ifelse(cfg_raw$gsea, "启用", "禁用"), cfg_raw$gsea_method))

  # 根据分析模式调整参数
  if (!is.null(cfg_raw$analysis_mode) && cfg_raw$analysis_mode == "gsea-only") {
    log_msg("info", "检测到 GSEA-only 模式，自动调整参数")
    # 禁用 GO/KEGG
    cfg$analysis$go_kegg$enabled <- FALSE
    # 确保 GSEA 启用
    cfg$analysis$gsea$enabled <- TRUE
  }

  # 创建输出目录
  prepare_dirs(c(
    paths$results_dir,
    paths$intermediate_dir,
    paths$qc_dir,
    paths$logs_dir,
    paths$tables_dir,
    paths$plots_dir,
    paths$go_dir,
    paths$kegg_dir,
    paths$gsea_dir
  ))

  # 保存配置文件到用户的结果目录，方便报告生成脚本读取
  user_config_file <- file.path(paths$results_dir, paste0(SCRIPT_NAME, ".config.ini"))
  writeLines(config_lines, user_config_file)

  # 初始化结果列表
  results <- list(
    go_kegg = NULL,
    gsea = NULL
  )

  # ---- 运行GO/KEGG分析 ----
  if (get_cfg(cfg, "analysis.go_kegg.enabled", FALSE)) {
    log_msg("info", "启用GO/KEGG分析")

    # 检查输入文件
    if (!file.exists(paths$gene_list)) {
      stop(sprintf("GO/KEGG基因列表文件不存在: %s", paths$gene_list))
    }

    results$go_kegg <- safe_run(
      run_go_kegg_pipeline(cfg, paths),
      "GO/KEGG富集分析"
    )
  } else {
    log_msg("info", "跳过GO/KEGG分析（未启用）")
  }

  # ---- 运行GSEA分析 ----
  if (get_cfg(cfg, "analysis.gsea.enabled", FALSE)) {
    gsea_method <- get_cfg(cfg, "analysis.gsea.method", "prerank")
    log_msg("info", sprintf("启用GSEA分析 (方法: %s)", gsea_method))

    # 根据方法检查不同的输入文件
    missing_files <- c()

    # GMT文件是必须的
    if (!file.exists(paths$gmt_file)) {
      missing_files <- c(missing_files, sprintf("GMT文件: %s", paths$gmt_file))
    }

    if (gsea_method == "prerank") {
      # Pre-rank GSEA需要排序基因列表
      if (!file.exists(paths$ranked_gene_list)) {
        missing_files <- c(missing_files, sprintf("排序基因列表: %s", paths$ranked_gene_list))
      }
    } else if (gsea_method == "correlation") {
      # Correlation GSEA需要表达矩阵和目标基因
      if (!file.exists(paths$expression_matrix)) {
        missing_files <- c(missing_files, sprintf("表达矩阵: %s", paths$expression_matrix))
      }
      if (!file.exists(paths$target_genes)) {
        missing_files <- c(missing_files, sprintf("目标基因: %s", paths$target_genes))
      }
    } else if (gsea_method == "median_split") {
      # Median Split GSEA需要表达矩阵和目标基因
      if (!file.exists(paths$expression_matrix)) {
        missing_files <- c(missing_files, sprintf("表达矩阵: %s", paths$expression_matrix))
      }
      if (!file.exists(paths$target_genes)) {
        missing_files <- c(missing_files, sprintf("目标基因: %s", paths$target_genes))
      }
    }

    if (length(missing_files) > 0) {
      stop(sprintf("GSEA缺少输入文件:\n  %s", paste(missing_files, collapse = "\n  ")))
    }

    results$gsea <- safe_run(
      run_gsea_pipeline(cfg, paths),
      "GSEA富集分析"
    )
  } else {
    log_msg("info", "跳过GSEA分析（未启用）")
  }

  # ---- 生成中文统一报告 ----
  log_msg("info", "生成中文统一报告")

  tryCatch({
    # 调用中文报告生成脚本
    report_script <- file.path(MODULE_ROOT, "scripts", "generate_enrichment_report.R")
    if (file.exists(report_script)) {
      # 构建输出报告路径
      result_dir_name <- basename(paths$results_dir)
      report_output <- file.path(MODULE_ROOT, "report", paste0("enrichment_report_", TIMESTAMP, ".docx"))

      # 确保 report 目录存在
      if (!dir.exists(file.path(MODULE_ROOT, "report"))) {
        dir.create(file.path(MODULE_ROOT, "report"), recursive = TRUE)
      }

      # 调用报告生成脚本
      system2("Rscript",
              args = c(report_script,
                      "--input", paths$results_dir,
                      "--output", report_output),
              stdout = TRUE,
              stderr = TRUE)

      if (file.exists(report_output)) {
        log_msg("info", sprintf("中文报告已生成: %s", report_output))
        cat(sprintf("\n中文报告已生成: %s\n", report_output))
      }
    } else {
      log_msg("warn", "报告生成脚本不存在，使用默认报告")
      # 回退到原来的报告生成
      report_file <- generate_unified_report(
        cfg = cfg,
        paths = paths,
        go_kegg_results = results$go_kegg,
        gsea_results = results$gsea
      )
      log_msg("info", sprintf("报告已生成: %s", report_file))
    }
  }, error = function(e) {
    log_msg("error", sprintf("报告生成失败: %s", e$message))
  })

  # ---- 生成 Quarto 报告 ----
  if (isTRUE(cfg_raw$quarto_report)) {
    log_msg("info", "生成 Quarto 报告")

    tryCatch({
      # 检查 Quarto 是否可用
      quarto_path <- Sys.which("quarto")
      if (nchar(quarto_path) == 0) {
        log_msg("warn", "Quarto 未安装，跳过报告生成。请访问 https://quarto.org 安装。")
      } else {
        # 加载报告生成器
        script_dir <- getSrcDirectory(function(x) x)
        if (is.null(script_dir) || script_dir == "") {
          script_dir <- "."
        }

        report_generator <- file.path(script_dir, "..", "report_templates", "generate_report.R")

        if (!file.exists(report_generator)) {
          report_generator <- file.path(dirname(script_dir), "report_templates", "generate_report.R")
        }

        if (file.exists(report_generator)) {
          source(report_generator)

          quarto_report <- generate_enrichment_quarto_report(
            results_dir = paths$results_dir,
            output_format = cfg_raw$quarto_format,
            project_name = basename(cfg_raw$work_dir),
            author = cfg_raw$quarto_author,
            species = cfg_raw$species,
            analysis_type = ifelse(cfg_raw$go_kegg && cfg_raw$gsea, "combined",
                                   ifelse(cfg_raw$go_kegg, "go_kegg", "gsea"))
          )

          if (!is.null(quarto_report)) {
            log_msg("info", sprintf("Quarto 报告已生成: %s", quarto_report))
          }
        } else {
          log_msg("warn", sprintf("报告生成器未找到: %s", report_generator))
        }
      }
    }, error = function(e) {
      log_msg("error", sprintf("Quarto 报告生成失败: %s", e$message))
    })
  }

  # ---- 输出摘要 ----
  cat("\n============================================================\n")
  cat("                      分析完成摘要\n")
  cat("============================================================\n")

  # GO/KEGG摘要
  if (isTRUE(cfg$go_kegg) && !is.null(results$go_kegg)) {
    go_count <- if (!is.null(results$go_kegg$go$result)) nrow(results$go_kegg$go$result) else 0
    kegg_count <- if (!is.null(results$go_kegg$kegg$result)) nrow(results$go_kegg$kegg$result) else 0
    cat(sprintf("GO/KEGG分析:\n"))
    cat(sprintf("  - GO显著条目: %d\n", go_count))
    cat(sprintf("  - KEGG显著通路: %d\n", kegg_count))
  }

  # GSEA摘要
  if (get_cfg(cfg, "analysis.gsea.enabled", FALSE) && !is.null(results$gsea)) {
    gsea_method <- get_cfg(cfg, "analysis.gsea.method", "prerank")

    if (gsea_method == "prerank") {
      # Pre-rank GSEA结果
      pathway_count <- if (!is.null(results$gsea$gsea)) nrow(results$gsea$gsea@result) else 0
      cat(sprintf("GSEA分析 (Pre-rank):\n"))
      cat(sprintf("  - 输入基因数: %d\n", length(results$gsea$gene_rank)))
      cat(sprintf("  - 显著通路数: %d\n", pathway_count))
    } else {
      # Correlation GSEA结果
      n_genes <- length(results$gsea$gsea_results)
      total_pathways <- sum(sapply(results$gsea$gsea_results, function(x) nrow(x@result)))
      cat(sprintf("GSEA分析 (Correlation):\n"))
      cat(sprintf("  - 分析基因数: %d\n", n_genes))
      cat(sprintf("  - 总显著通路数: %d\n", total_pathways))
    }
  }

  # 输出文件位置
  cat("\n----------------------------------------\n")
  cat("输出文件\n")
  cat("----------------------------------------\n")
  cat(sprintf("日志文件: %s.%s.log\n", SCRIPT_NAME, TIMESTAMP))
  cat(sprintf("配置文件: %s\n", basename(config_file)))
  cat(sprintf("结果目录: %s\n", paths$results_dir))

  # 运行时间
  duration <- difftime(Sys.time(), script_start_time, units = "mins")
  cat(sprintf("\n总运行时间: %.2f 分钟\n", as.numeric(duration)))

  cat("\n================================================================================\n")
  cat("分析成功完成\n")
  cat("================================================================================\n")

  # 保存 sessionInfo 到 logs/
  session_info_obj <- sessionInfo()
  session_info_txt <- file.path(log_dir, paste0("sessionInfo_", TIMESTAMP, ".txt"))
  writeLines(capture.output(session_info_obj), session_info_txt)

  # 保存 qs2 到 results/（节点快照 + sessionInfo 快照）
  dataset_name <- basename(normalizePath(paths$results_dir, mustWork = FALSE))
  enrichment_node_path <- file.path(paths$results_dir, sprintf("%s_enrichment_node.qs2", dataset_name))
  session_qs2_path <- file.path(paths$results_dir, sprintf("%s_sessionInfo.%s.qs2", dataset_name, TIMESTAMP))
  gsea_rows <- 0
  if (!is.null(results$gsea)) {
    gsea_rows <- tryCatch({
      if (!is.null(results$gsea$gsea) && methods::is(results$gsea$gsea, "gseaResult")) {
        nrow(results$gsea$gsea@result)
      } else if (!is.null(results$gsea$gsea_results)) {
        sum(vapply(results$gsea$gsea_results, function(x) {
          if (!is.null(x) && methods::is(x, "gseaResult")) nrow(x@result) else 0L
        }, integer(1)))
      } else {
        0L
      }
    }, error = function(e) 0L)
  }

  save_qs2_snapshot(list(
    timestamp = TIMESTAMP,
    analysis_mode = ifelse(is.null(cfg_raw$analysis_mode), "full", cfg_raw$analysis_mode),
    species = cfg_raw$species,
    config_file = user_config_file,
    summary = list(
      go_enabled = isTRUE(cfg$analysis$go_kegg$enabled),
      gsea_enabled = isTRUE(cfg$analysis$gsea$enabled),
      go_rows = if (!is.null(results$go_kegg) && !is.null(results$go_kegg$go$result)) nrow(results$go_kegg$go$result) else 0,
      kegg_rows = if (!is.null(results$go_kegg) && !is.null(results$go_kegg$kegg$result)) nrow(results$go_kegg$kegg$result) else 0,
      gsea_rows = gsea_rows
    ),
    output_paths = paths
  ), enrichment_node_path)

  save_qs2_snapshot(list(
    timestamp = TIMESTAMP,
    session_info = session_info_obj
  ), session_qs2_path)

  log_msg("info", "分析流程完成")

  invisible(results)
}

# ---- 执行 ----
tryCatch(
  main(),
  error = function(e) {
    cat(sprintf("[FATAL] %s\n", e$message))
    traceback()
    quit(status = 1)
  }
)


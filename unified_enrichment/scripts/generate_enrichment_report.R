#!/usr/bin/env Rscript

# =============================================================================
# 统一富集分析报告生成器
# Unified Enrichment Report Generator
# =============================================================================
# 生成中文 Word (.docx) 报告，包含 GO/KEGG 和 GSEA 分析结果
# 使用 officer 和 flextable 包
# =============================================================================

library(optparse)
library(officer)
library(flextable)
library(ggplot2)
library(dplyr)
library(magrittr)
library(data.table)

# Parse arguments
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "Result directory containing enrichment results"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Output report file (.docx)"),
  make_option(c("-T", "--timestamp"), type = "character", default = NULL,
              help = "Timestamp for this run")
)

opt <- parse_args(OptionParser(option_list = option_list))

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  hit <- args[grep(file_arg, args)]
  if (length(hit) == 0) return(getwd())
  dirname(normalizePath(sub(file_arg, "", hit[[1]]), mustWork = FALSE))
}
script_dir <- get_script_dir()
project_root <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)
report_template <- file.path(project_root, "templates", "unified_enrichment_report_template.docx")

# Default timestamp
if (is.null(opt$timestamp)) {
  opt$timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
}

# Default output
if (is.null(opt$output)) {
  opt$output <- paste0("Enrichment_Report_", opt$timestamp, ".docx")
}

# Auto-detect latest result directory (prevents reading stale runs)
detect_latest_result_dir <- function(project_root) {
  candidate_cfg_files <- list.files(
    path = project_root,
    pattern = "^unified_enrichment\\.config\\.ini$",
    recursive = TRUE,
    full.names = TRUE
  )
  if (length(candidate_cfg_files) == 0) return(NULL)

  # Only keep config files under a run result directory
  keep <- grepl("(^|/)results(/|$)", gsub("\\\\", "/", candidate_cfg_files))
  candidate_cfg_files <- candidate_cfg_files[keep]
  if (length(candidate_cfg_files) == 0) return(NULL)

  finfo <- file.info(candidate_cfg_files)
  ord <- order(finfo$mtime, decreasing = TRUE, na.last = NA)
  if (length(ord) == 0) return(NULL)
  dirname(candidate_cfg_files[ord[1]])
}

# Default input directory
if (is.null(opt$input)) {
  detected <- detect_latest_result_dir(project_root)
  if (!is.null(detected) && dir.exists(detected)) {
    opt$input <- detected
  } else {
    opt$input <- "results/unified_enrichment"
  }
}

cat("================================================================================\n")
cat("  统一富集分析报告生成器\n")
cat("  结果目录:", opt$input, "\n")
cat("  时间戳:", opt$timestamp, "\n")
cat("================================================================================\n\n")

# Set paths
result_dir <- opt$input
go_dir <- file.path(result_dir, "go")
kegg_dir <- file.path(result_dir, "kegg")
gsea_dir <- file.path(result_dir, "gsea")
qc_dir <- file.path(result_dir, "qc")
plots_dir <- file.path(result_dir, "plots")

# Read config file
config_file <- file.path(result_dir, "unified_enrichment.config.ini")
species <- "hsa"
run_go_kegg <- TRUE  # 默认运行 GO/KEGG
run_gsea <- TRUE     # 默认运行 GSEA
if (file.exists(config_file)) {
  config_lines <- readLines(config_file)
  for (line in config_lines) {
    if (grepl("^species=", line)) {
      species <- sub(".*=", "", line)
    }
    if (grepl("^go_kegg=", line)) {
      val <- tolower(sub(".*=", "", line))
      run_go_kegg <- val %in% c("true", "yes", "1")
    }
    if (grepl("^gsea=", line)) {
      val <- tolower(sub(".*=", "", line))
      run_gsea <- val %in% c("true", "yes", "1")
    }
  }
}
cat("分析模式: GO/KEGG=", run_go_kegg, ", GSEA=", run_gsea, "\n")

# Species mapping
species_names <- list(
  hsa = "人类 (Homo sapiens)",
  mmu = "小鼠 (Mus musculus)",
  rno = "大鼠 (Rattus norvegicus)"
)
species_name <- species_names[[species]] %||% species

# Read GO results
go_results <- NULL
go_file <- file.path(go_dir, "go_enrichment_results.csv")
if (file.exists(go_file)) {
  go_results <- as.data.frame(fread(go_file, stringsAsFactors = FALSE))
  cat("Loaded GO results:", nrow(go_results), "rows\n")
}

# Read KEGG results
kegg_results <- NULL
kegg_file <- file.path(kegg_dir, "kegg_enrichment_results.csv")
if (file.exists(kegg_file)) {
  kegg_results <- as.data.frame(fread(kegg_file, stringsAsFactors = FALSE))
  cat("Loaded KEGG results:", nrow(kegg_results), "rows\n")
}

# Read GSEA results
gsea_results <- NULL
gsea_file <- file.path(gsea_dir, "prerank_gsea_results.csv")
if (file.exists(gsea_file)) {
  gsea_results <- as.data.frame(fread(gsea_file, stringsAsFactors = FALSE))
  cat("Loaded GSEA results:", nrow(gsea_results), "rows\n")
}

# Read QC results
qc_results <- NULL
qc_file <- file.path(qc_dir, "PreRank_GSEA_QC_Report.csv")
if (file.exists(qc_file)) {
  qc_results <- as.data.frame(fread(qc_file, stringsAsFactors = FALSE))
  cat("Loaded QC results:", nrow(qc_results), "rows\n")
}

# 只要结果文件存在，就强制纳入对应章节，避免被配置开关误关
if (!is.null(go_results) && nrow(go_results) > 0) run_go_kegg <- TRUE
if (!is.null(kegg_results) && nrow(kegg_results) > 0) run_go_kegg <- TRUE
if (!is.null(gsea_results) && nrow(gsea_results) > 0) run_gsea <- TRUE
cat("报告纳入模式(按结果文件修正): GO/KEGG=", run_go_kegg, ", GSEA=", run_gsea, "\n")

# Create Word document (prefer project-specific template)
doc <- if (file.exists(report_template)) {
  cat("Using report template:", report_template, "\n")
  read_docx(path = report_template)
} else {
  read_docx()
}

# ---- helper functions for unified report layout ----
build_mixed_fpar <- function(text, align = "left", font_size = 12, indent = FALSE) {
  if (indent && nzchar(trimws(text))) text <- paste0("\u3000\u3000", text)
  chars <- strsplit(text, "", fixed = TRUE)[[1]]
  if (length(chars) == 0) return(fpar("", fp_p = fp_par(text.align = align)))

  # 中文字符与中文标点使用宋体，其余（英文/数字/ASCII符号）使用 Times New Roman
  is_cn <- grepl("[一-龥，。；：、（）【】《》“”‘’？！％]", chars)
  runs <- list()
  start <- 1L
  for (i in seq_along(chars)) {
    is_break <- i == length(chars) || is_cn[[i + 1L]] != is_cn[[i]]
    if (is_break) {
      seg <- paste(chars[start:i], collapse = "")
      ff <- if (is_cn[[start]]) "SimSun" else "Times New Roman"
      runs[[length(runs) + 1L]] <- ftext(
        seg,
        prop = fp_text(font.family = ff, font.size = font_size)
      )
      start <- i + 1L
    }
  }
  do.call(fpar, c(runs, list(fp_p = fp_par(text.align = align))))
}

add_center_par <- function(doc, text) {
  body_add_fpar(doc, value = build_mixed_fpar(text, align = "center"), style = "Normal")
}

add_result_par <- function(doc, text) {
  body_add_fpar(doc, value = build_mixed_fpar(text, align = "left", indent = TRUE), style = "Normal")
}

add_center_img <- function(doc, img_path, width = 6, height = 5) {
  img_block <- fpar(
    external_img(src = img_path, width = width, height = height),
    fp_p = fp_par(text.align = "center")
  )
  body_add_fpar(doc, value = img_block, style = "Normal")
}

# 对最终 docx 做统一后处理：
# 1) 字体：中文 SimSun，英文/数字 Times New Roman
# 2) 图表标题/图注段落强制居中（避免模板样式覆盖）
postprocess_docx_fonts <- function(docx_path) {
  py_code <- c(
    "import sys, zipfile, os",
    "from xml.etree import ElementTree as ET",
    "docx = sys.argv[1]",
    "ns = {'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}",
    "ET.register_namespace('w', ns['w'])",
    "with zipfile.ZipFile(docx,'r') as zin:",
    "    xml = zin.read('word/document.xml')",
    "root = ET.fromstring(xml)",
    "import re",
    "for p in root.findall('.//w:p', ns):",
    "    p_txt = ''.join([t.text for t in p.findall('.//w:t', ns) if t.text]).strip()",
    "    is_title = re.match(r'^(图|表)\\d+：', p_txt) is not None",
    "    is_caption = p_txt.startswith('图注：')",
    "    if is_title or is_caption:",
    "        pPr = p.find('w:pPr', ns)",
    "        if pPr is None:",
    "            pPr = ET.Element('{%s}pPr' % ns['w'])",
    "            p.insert(0, pPr)",
    "        jc = pPr.find('w:jc', ns)",
    "        if jc is None:",
    "            jc = ET.SubElement(pPr, '{%s}jc' % ns['w'])",
    "        jc.set('{%s}val' % ns['w'], 'center')",
    "table_run_ids = set()",
    "for tbl in root.findall('.//w:tbl', ns):",
    "    for tr in tbl.findall('.//w:r', ns):",
    "        table_run_ids.add(id(tr))",
    "for r in root.findall('.//w:r', ns):",
    "    rPr = r.find('w:rPr', ns)",
    "    if rPr is None:",
    "        rPr = ET.Element('{%s}rPr' % ns['w'])",
    "        r.insert(0, rPr)",
    "    rFonts = rPr.find('w:rFonts', ns)",
    "    if rFonts is None:",
    "        rFonts = ET.SubElement(rPr, '{%s}rFonts' % ns['w'])",
    "    rFonts.set('{%s}ascii' % ns['w'], 'Times New Roman')",
    "    rFonts.set('{%s}hAnsi' % ns['w'], 'Times New Roman')",
    "    rFonts.set('{%s}cs' % ns['w'], 'Times New Roman')",
    "    rFonts.set('{%s}eastAsia' % ns['w'], 'SimSun')",
    "    if id(r) not in table_run_ids:",
    "        sz = rPr.find('w:sz', ns)",
    "        if sz is None:",
    "            sz = ET.SubElement(rPr, '{%s}sz' % ns['w'])",
    "        sz.set('{%s}val' % ns['w'], '24')",
    "        szCs = rPr.find('w:szCs', ns)",
    "        if szCs is None:",
    "            szCs = ET.SubElement(rPr, '{%s}szCs' % ns['w'])",
    "        szCs.set('{%s}val' % ns['w'], '24')",
    "new_xml = ET.tostring(root, encoding='utf-8', xml_declaration=True)",
    "tmp = docx + '.tmp'",
    "with zipfile.ZipFile(docx,'r') as zin, zipfile.ZipFile(tmp,'w',zipfile.ZIP_DEFLATED) as zout:",
    "    for item in zin.infolist():",
    "        data = zin.read(item.filename)",
    "        if item.filename == 'word/document.xml':",
    "            data = new_xml",
    "        zout.writestr(item.filename, data)",
    "os.replace(tmp, docx)"
  )
  py_file <- tempfile(fileext = ".py")
  writeLines(py_code, py_file, useBytes = TRUE)
  on.exit(unlink(py_file), add = TRUE)
  tryCatch({
    system2("python3", c(py_file, docx_path), stdout = TRUE, stderr = TRUE)
  }, error = function(e) {
    warning("Font postprocess failed: ", e$message)
  })
}

# ========== Title ==========
doc <- doc %>%
  body_add_par("统一富集分析报告", style = "Normal") %>%
  body_add_par("", style = "Normal")

# ========== Methods Section ==========
doc <- doc %>%
  body_add_par("1. 材料与方法", style = "Normal") %>%
  body_add_par("", style = "Normal")

# 根据分析模式添加方法说明
method_counter <- 1

if (run_go_kegg) {
  doc <- doc %>%
    body_add_par(sprintf("1.%d GO/KEGG 富集分析", method_counter), style = "Normal")
  doc <- doc %>%
    add_result_par("使用R包clusterProfiler（v4.10.0，https://www.bioconductor.org/packages/release/bioc/html/clusterProfiler.html）对差异基因进行GO与KEGG富集分析。") %>%
    body_add_par("", style = "Normal") %>%
    add_result_par("GO 分析主要涉及生物过程 (Biological Process, BP)、分子功能 (Molecular Function, MF) 和细胞组分 (Cellular Component, CC) 三部分。显著阈值设定为校正后 p 值 (p.adjust) < 0.05。") %>%
    body_add_par("", style = "Normal")
  method_counter <- method_counter + 1
}

if (run_gsea) {
  doc <- doc %>%
    body_add_par(sprintf("1.%d 基因集富集分析 (GSEA)", method_counter), style = "Normal")
  doc <- doc %>%
    add_result_par("基因集富集分析 (Gene Set Enrichment Analysis, GSEA) 是一种用于评估基因集合在整个基因组的分布是否显著的方法。与过表达分析不同，GSEA 不需要预先设定差异基因阈值，而是考虑所有基因的表达变化。本分析使用 R 包 clusterProfiler (Version 4.10.0) 进行 prerank GSEA 分析。") %>%
    body_add_par("", style = "Normal") %>%
    add_result_par("基因集来源于 MSigDB 数据库 (http://www.broadinstitute.org/gsea/msigdb/index.jsp)。显著阈值设定为校正后 p 值 (p.adjust) < 0.05 且标准化富集分数 |NES| > 1。") %>%
    body_add_par("", style = "Normal")
  method_counter <- method_counter + 1
}

doc <- doc %>%
  body_add_par(sprintf("1.%d 统计分析", method_counter), style = "Normal")
doc <- doc %>%
  add_result_par("GO/KEGG 富集分析采用超几何检验计算显著性。GSEA 使用富集分数 (Enrichment Score, ES) 评估基因集的富集程度，并通过标准化富集分数 (Normalized Enrichment Score, NES) 消除基因集大小的影响。所有分析均采用 Benjamini-Hochberg (BH) 方法进行多重检验校正。") %>%
  body_add_par("", style = "Normal")

# ========== Results Section ==========
doc <- doc %>%
  body_add_par("2. 结果", style = "Normal") %>%
  body_add_par("", style = "Normal")

# ---- Summary ----
if (run_go_kegg || run_gsea) {
  if (!is.null(qc_results) && nrow(qc_results) > 0) {
    get_qc_metric <- function(name) {
      if (all(c("Item", "Value") %in% colnames(qc_results))) {
        val <- qc_results$Value[qc_results$Item == name]
        if (length(val) > 0) return(as.character(val[1]))
      }
      return(NA_character_)
    }
    total_genes <- suppressWarnings(as.numeric(get_qc_metric("Input_Genes")))
    matched_genes <- suppressWarnings(as.numeric(get_qc_metric("GMT_Genes_Matched")))
    matched_rate <- get_qc_metric("GMT_Match_Rate")
    if (is.na(total_genes) && "Total_Genes" %in% colnames(qc_results)) {
      total_genes <- suppressWarnings(as.numeric(qc_results$Total_Genes[1]))
    }
    if (is.na(matched_genes) && "Valid_Genes" %in% colnames(qc_results)) {
      matched_genes <- suppressWarnings(as.numeric(qc_results$Valid_Genes[1]))
    }
    if (!is.na(total_genes) || !is.na(matched_genes) || !is.na(suppressWarnings(as.numeric(gsub("%", "", matched_rate))))) {
      doc <- doc %>%
        body_add_par(sprintf(
          "输入基因数为%s，基因集匹配基因数为%s（匹配率%s）。",
          ifelse(is.na(total_genes), "NA", as.character(as.integer(total_genes))),
          ifelse(is.na(matched_genes), "NA", as.character(as.integer(matched_genes))),
          ifelse(is.na(matched_rate) || matched_rate == "", "NA", matched_rate)
        ), style = "Normal") %>%
        body_add_par("", style = "Normal")
    }
  }
}

# ========== Results Sections ==========
result_section_num <- 1
figure_num <- 1

# ---- GO Results ----
if (run_go_kegg) {
  doc <- doc %>% body_add_par(sprintf("2.%d GO 富集分析结果", result_section_num), style = "Normal")
  result_section_num <- result_section_num + 1

  if (!is.null(go_results) && nrow(go_results) > 0) {
    go_sig <- go_results[go_results$p.adjust < 0.05, ]

    if (nrow(go_sig) > 0) {
      go_bp_count <- sum(go_results$ONTOLOGY == "BP" & go_results$p.adjust < 0.05, na.rm = TRUE)
      go_cc_count <- sum(go_results$ONTOLOGY == "CC" & go_results$p.adjust < 0.05, na.rm = TRUE)
      go_mf_count <- sum(go_results$ONTOLOGY == "MF" & go_results$p.adjust < 0.05, na.rm = TRUE)

      # 汇总说明
      doc <- doc %>%
        add_result_par(sprintf("以 p.adjust < 0.05 为阈值，GO 富集分析共筛选出 %d 个显著条目，其中 GO_BP (生物过程) 富集了 %d 条通路，GO_CC (细胞组分) 富集了 %d 条通路，GO_MF (分子功能) 富集了 %d 条通路。", nrow(go_sig), go_bp_count, go_cc_count, go_mf_count)) %>%
        body_add_par("", style = "Normal")

      # 获取 top 通路
      top_go <- go_sig %>%
        group_by(ONTOLOGY) %>%
        arrange(p.adjust) %>%
        slice_head(n = 5) %>%
        ungroup()

      # 列出 top 通路
      doc <- doc %>% add_result_par("显著富集的 Top 通路包括：")

      for (ont in c("BP", "CC", "MF")) {
        ont_name <- switch(ont, "BP" = "生物过程", "CC" = "细胞组分", "MF" = "分子功能")
        ont_sig <- go_sig[go_sig$ONTOLOGY == ont, ]
        if (nrow(ont_sig) > 0) {
          ont_top <- ont_sig[order(ont_sig$p.adjust), ][1:min(3, nrow(ont_sig)), ]
          pathways <- paste(ont_top$Description, collapse = "、")
          doc <- doc %>% add_result_par(sprintf("%s (GO_%s): %s", ont_name, ont, pathways))
        }
      }
      doc <- doc %>% body_add_par("", style = "Normal")

      # 创建表格
      cols_map <- c(
        ONTOLOGY = "Ontology",
        ID = "GO ID",
        Description = "Term",
        GeneRatio = "Gene Ratio",
        Count = "Gene Count",
        p.adjust = "Adjusted P-value"
      )

      cols_exist <- names(cols_map)[names(cols_map) %in% colnames(top_go)]

      if (length(cols_exist) > 0) {
        table_data <- top_go[, cols_exist, drop = FALSE]
        colnames(table_data) <- cols_map[cols_exist]

        ft <- flextable(table_data) %>%
          bold(part = "header", bold = TRUE) %>%
          fontsize(size = 9, part = "all") %>%
          font(fontname = "Times New Roman", part = "all") %>%
          align(align = "center", part = "all") %>%
          align(align = "left", part = "body", j = 2:3) %>%
          border_remove() %>%
          border(part = "all", border.top = officer::fp_border(width = 1.5, color = "black")) %>%
          border(part = "header", border.bottom = officer::fp_border(width = 0.5, color = "black")) %>%
          border(part = "body", border.bottom = officer::fp_border(width = 1.5, color = "black"), i = nrow(table_data)) %>%
          set_table_properties(layout = "autofit", width = 1)

        top_go_idx <- which.min(suppressWarnings(as.numeric(go_sig$p.adjust)))
        top_go_term <- if (length(top_go_idx) > 0) as.character(go_sig$Description[top_go_idx]) else "NA"
        top_go_p <- if (length(top_go_idx) > 0) suppressWarnings(as.numeric(go_sig$p.adjust[top_go_idx])) else NA_real_
        doc <- doc %>% add_center_par("表1：GO Enrichment Results") %>%
          body_add_flextable(ft, align = "center")
        go_sorted <- go_sig[order(go_sig$p.adjust), , drop = FALSE]
        go_top3 <- unique(as.character(go_sorted$Description[1:min(3, nrow(go_sorted))]))
        go_top3_txt <- if (length(go_top3) > 0) paste(go_top3, collapse = "、") else "NA"
        go_max_count <- suppressWarnings(max(as.numeric(go_sig$Count), na.rm = TRUE))
        if (!is.finite(go_max_count)) go_max_count <- NA_real_
        doc <- doc %>% add_result_par(sprintf(
          "表1结果显示，共筛选到%d个显著GO条目（BP=%d, CC=%d, MF=%d）；最显著条目为%s（adjusted p-value=%s），Top条目主要包括%s；单条目最大Gene Count为%s。相关结果表：`go_enrichment_results.csv`。",
          nrow(go_sig), go_bp_count, go_cc_count, go_mf_count, top_go_term,
          ifelse(is.na(top_go_p), "NA", format(top_go_p, scientific = TRUE, digits = 3)),
          go_top3_txt,
          ifelse(is.na(go_max_count), "NA", as.character(as.integer(go_max_count)))
        ))
      }
      doc <- doc %>% body_add_par("", style = "Normal")

      # GO Figure
      go_plot_png <- file.path(go_dir, "go_enrichment_plot.png")
      go_plot_pdf <- file.path(go_dir, "go_enrichment_plot.pdf")

      go_stats_file <- file.path(go_dir, "go_figure_stats.csv")
      go_stats <- if (file.exists(go_stats_file)) tryCatch(as.data.frame(fread(go_stats_file)), error = function(e) NULL) else NULL
      if (is.null(go_stats) || nrow(go_stats) == 0) {
        top_go_idx <- which.min(suppressWarnings(as.numeric(go_sig$p.adjust)))
        top_go_term <- if (length(top_go_idx) > 0) as.character(go_sig$Description[top_go_idx]) else "NA"
        top_go_p <- if (length(top_go_idx) > 0) suppressWarnings(as.numeric(go_sig$p.adjust[top_go_idx])) else NA_real_
        go_stats <- data.frame(
          metric = c("total_terms", "top1_term", "top1_pvalue"),
          value = c(
            as.character(nrow(go_sig)),
            top_go_term,
            ifelse(is.na(top_go_p), "NA", format(top_go_p, scientific = TRUE, digits = 3))
          ),
          stringsAsFactors = FALSE
        )
        fwrite(go_stats, go_stats_file)
      }
      if (file.exists(go_plot_png)) {
        doc <- doc %>% add_center_img(go_plot_png, width = 6, height = 5) %>%
          add_center_par(sprintf("图%d：GO Enrichment Dotplot", figure_num)) %>%
          add_center_par("图注：X-axis为RichFactor，Y-axis为GO Term，color为adjusted p-value，size为Gene Count。")
        get_stat <- function(metric) go_stats$value[go_stats$metric == metric][1]
        go_sorted <- go_sig[order(go_sig$p.adjust), , drop = FALSE]
        go_top3 <- unique(as.character(go_sorted$Description[1:min(3, nrow(go_sorted))]))
        go_top3_txt <- if (length(go_top3) > 0) paste(go_top3, collapse = "、") else "NA"
        doc <- doc %>% add_result_par(sprintf(
          "图%d结果显示，共识别%s个显著GO条目（BP=%d, CC=%d, MF=%d），最显著条目为%s（adjusted p-value=%s）；富集核心条目包括%s。相关结果表：`go_figure_stats.csv`。",
          figure_num,
          ifelse(length(get_stat("total_terms")) > 0, as.character(get_stat("total_terms")), "NA"),
          go_bp_count, go_cc_count, go_mf_count,
          ifelse(length(get_stat("top1_term")) > 0, as.character(get_stat("top1_term")), "NA"),
          ifelse(length(get_stat("top1_pvalue")) > 0, as.character(get_stat("top1_pvalue")), "NA"),
          go_top3_txt
        ))
        figure_num <- figure_num + 1
      } else if (file.exists(go_plot_pdf)) {
        doc <- doc %>% add_center_img(go_plot_pdf, width = 6, height = 5) %>%
          add_center_par(sprintf("图%d：GO Enrichment Dotplot", figure_num)) %>%
          add_center_par("图注：X-axis为RichFactor，Y-axis为GO Term，color为adjusted p-value，size为Gene Count。")
        get_stat <- function(metric) go_stats$value[go_stats$metric == metric][1]
        go_sorted <- go_sig[order(go_sig$p.adjust), , drop = FALSE]
        go_top3 <- unique(as.character(go_sorted$Description[1:min(3, nrow(go_sorted))]))
        go_top3_txt <- if (length(go_top3) > 0) paste(go_top3, collapse = "、") else "NA"
        doc <- doc %>% add_result_par(sprintf(
          "图%d结果显示，共识别%s个显著GO条目（BP=%d, CC=%d, MF=%d），最显著条目为%s（adjusted p-value=%s）；富集核心条目包括%s。相关结果表：`go_figure_stats.csv`。",
          figure_num,
          ifelse(length(get_stat("total_terms")) > 0, as.character(get_stat("total_terms")), "NA"),
          go_bp_count, go_cc_count, go_mf_count,
          ifelse(length(get_stat("top1_term")) > 0, as.character(get_stat("top1_term")), "NA"),
          ifelse(length(get_stat("top1_pvalue")) > 0, as.character(get_stat("top1_pvalue")), "NA"),
          go_top3_txt
        ))
        figure_num <- figure_num + 1
      }
      doc <- doc %>% body_add_par("", style = "Normal")
    } else {
      doc <- doc %>% add_result_par("以 p.adjust < 0.05 为阈值，无显著 GO 条目。")
      doc <- doc %>% body_add_par("", style = "Normal")
    }
  } else {
    doc <- doc %>% add_result_par("无 GO 结果可用。")
    doc <- doc %>% body_add_par("", style = "Normal")
  }

  doc <- doc %>% body_add_par("", style = "Normal")
  result_section_num <- result_section_num + 1
}

# ---- KEGG Results ----
if (run_go_kegg) {
  doc <- doc %>% body_add_par(sprintf("2.%d KEGG 通路富集分析结果", result_section_num), style = "Normal")
  result_section_num <- result_section_num + 1

  if (!is.null(kegg_results) && nrow(kegg_results) > 0) {
    kegg_sig <- kegg_results[kegg_results$p.adjust < 0.05, ]

    if (nrow(kegg_sig) > 0) {
      # 汇总说明
      doc <- doc %>%
        add_result_par(sprintf("以 p.adjust < 0.05 为阈值，KEGG 通路富集分析共筛选出 %d 个显著通路。", nrow(kegg_sig)))

      # 列出 top 通路
      top5_kegg <- kegg_sig[order(kegg_sig$p.adjust), ][1:min(5, nrow(kegg_sig)), ]
      pathways <- paste(top5_kegg$Description, collapse = "、")
      doc <- doc %>% add_result_par(sprintf("显著富集的 Top 5 通路包括：%s。", pathways))
      doc <- doc %>% body_add_par("", style = "Normal")

      # 创建表格
      top_kegg <- kegg_sig %>%
        arrange(p.adjust) %>%
        slice_head(n = 20)

      cols_map <- c(
        ID = "KEGG ID",
        Description = "Pathway",
        GeneRatio = "Gene Ratio",
        Count = "Gene Count",
        p.adjust = "Adjusted P-value"
      )

      cols_exist <- names(cols_map)[names(cols_map) %in% colnames(top_kegg)]

      if (length(cols_exist) > 0) {
        table_data <- top_kegg[, cols_exist, drop = FALSE]
        colnames(table_data) <- cols_map[cols_exist]

        ft <- flextable(table_data) %>%
          bold(part = "header", bold = TRUE) %>%
          fontsize(size = 9, part = "all") %>%
          font(fontname = "Times New Roman", part = "all") %>%
          align(align = "center", part = "all") %>%
          align(align = "left", part = "body", j = 2) %>%
          border_remove() %>%
          border(part = "all", border.top = officer::fp_border(width = 1.5, color = "black")) %>%
          border(part = "header", border.bottom = officer::fp_border(width = 0.5, color = "black")) %>%
          border(part = "body", border.bottom = officer::fp_border(width = 1.5, color = "black"), i = nrow(table_data)) %>%
          set_table_properties(layout = "autofit", width = 1)

        top_kegg_idx <- which.min(suppressWarnings(as.numeric(kegg_sig$p.adjust)))
        top_kegg_term <- if (length(top_kegg_idx) > 0) as.character(kegg_sig$Description[top_kegg_idx]) else "NA"
        top_kegg_p <- if (length(top_kegg_idx) > 0) suppressWarnings(as.numeric(kegg_sig$p.adjust[top_kegg_idx])) else NA_real_
        doc <- doc %>% add_center_par(sprintf("表%d：KEGG Enrichment Results", ifelse(run_go_kegg, 2, 1))) %>%
          body_add_flextable(ft, align = "center")
        kegg_sorted <- kegg_sig[order(kegg_sig$p.adjust), , drop = FALSE]
        kegg_top3 <- unique(as.character(kegg_sorted$Description[1:min(3, nrow(kegg_sorted))]))
        kegg_top3_txt <- if (length(kegg_top3) > 0) paste(kegg_top3, collapse = "、") else "NA"
        kegg_max_count <- suppressWarnings(max(as.numeric(kegg_sig$Count), na.rm = TRUE))
        if (!is.finite(kegg_max_count)) kegg_max_count <- NA_real_
        doc <- doc %>% add_result_par(sprintf(
          "表%d结果显示，共筛选到%d条显著KEGG通路；最显著通路为%s（adjusted p-value=%s），Top通路包括%s；单通路最大Gene Count为%s。相关结果表：`kegg_enrichment_results.csv`。",
          ifelse(run_go_kegg, 2, 1), nrow(kegg_sig), top_kegg_term,
          ifelse(is.na(top_kegg_p), "NA", format(top_kegg_p, scientific = TRUE, digits = 3)),
          kegg_top3_txt,
          ifelse(is.na(kegg_max_count), "NA", as.character(as.integer(kegg_max_count)))
        ))
      }
      doc <- doc %>% body_add_par("", style = "Normal")

      # KEGG Figure
      kegg_plot_png <- file.path(kegg_dir, "kegg_enrichment_plot.png")
      kegg_plot_pdf <- file.path(kegg_dir, "kegg_enrichment_plot.pdf")

      kegg_stats_file <- file.path(kegg_dir, "kegg_figure_stats.csv")
      kegg_stats <- if (file.exists(kegg_stats_file)) tryCatch(as.data.frame(fread(kegg_stats_file)), error = function(e) NULL) else NULL
      if (is.null(kegg_stats) || nrow(kegg_stats) == 0) {
        top_kegg_idx <- which.min(suppressWarnings(as.numeric(kegg_sig$p.adjust)))
        top_kegg_term <- if (length(top_kegg_idx) > 0) as.character(kegg_sig$Description[top_kegg_idx]) else "NA"
        top_kegg_p <- if (length(top_kegg_idx) > 0) suppressWarnings(as.numeric(kegg_sig$p.adjust[top_kegg_idx])) else NA_real_
        kegg_stats <- data.frame(
          metric = c("total_pathways", "top1_pathway", "top1_pvalue"),
          value = c(
            as.character(nrow(kegg_sig)),
            top_kegg_term,
            ifelse(is.na(top_kegg_p), "NA", format(top_kegg_p, scientific = TRUE, digits = 3))
          ),
          stringsAsFactors = FALSE
        )
        fwrite(kegg_stats, kegg_stats_file)
      }
      if (file.exists(kegg_plot_png)) {
        doc <- doc %>% add_center_img(kegg_plot_png, width = 6, height = 5) %>%
          add_center_par(sprintf("图%d：KEGG Enrichment Dotplot", figure_num)) %>%
          add_center_par("图注：X-axis为RichFactor，Y-axis为Pathway，color为adjusted p-value，size为Gene Count。")
        get_stat <- function(metric) kegg_stats$value[kegg_stats$metric == metric][1]
        kegg_sorted <- kegg_sig[order(kegg_sig$p.adjust), , drop = FALSE]
        kegg_top3 <- unique(as.character(kegg_sorted$Description[1:min(3, nrow(kegg_sorted))]))
        kegg_top3_txt <- if (length(kegg_top3) > 0) paste(kegg_top3, collapse = "、") else "NA"
        doc <- doc %>% add_result_par(sprintf(
          "图%d结果显示，共识别%s条显著KEGG通路，最显著通路为%s（adjusted p-value=%s）；主要富集通路包括%s。相关结果表：`kegg_figure_stats.csv`。",
          figure_num,
          ifelse(length(get_stat("total_pathways")) > 0, as.character(get_stat("total_pathways")), "NA"),
          ifelse(length(get_stat("top1_pathway")) > 0, as.character(get_stat("top1_pathway")), "NA"),
          ifelse(length(get_stat("top1_pvalue")) > 0, as.character(get_stat("top1_pvalue")), "NA"),
          kegg_top3_txt
        ))
        figure_num <- figure_num + 1
      } else if (file.exists(kegg_plot_pdf)) {
        doc <- doc %>% add_center_img(kegg_plot_pdf, width = 6, height = 5) %>%
          add_center_par(sprintf("图%d：KEGG Enrichment Dotplot", figure_num)) %>%
          add_center_par("图注：X-axis为RichFactor，Y-axis为Pathway，color为adjusted p-value，size为Gene Count。")
        get_stat <- function(metric) kegg_stats$value[kegg_stats$metric == metric][1]
        kegg_sorted <- kegg_sig[order(kegg_sig$p.adjust), , drop = FALSE]
        kegg_top3 <- unique(as.character(kegg_sorted$Description[1:min(3, nrow(kegg_sorted))]))
        kegg_top3_txt <- if (length(kegg_top3) > 0) paste(kegg_top3, collapse = "、") else "NA"
        doc <- doc %>% add_result_par(sprintf(
          "图%d结果显示，共识别%s条显著KEGG通路，最显著通路为%s（adjusted p-value=%s）；主要富集通路包括%s。相关结果表：`kegg_figure_stats.csv`。",
          figure_num,
          ifelse(length(get_stat("total_pathways")) > 0, as.character(get_stat("total_pathways")), "NA"),
          ifelse(length(get_stat("top1_pathway")) > 0, as.character(get_stat("top1_pathway")), "NA"),
          ifelse(length(get_stat("top1_pvalue")) > 0, as.character(get_stat("top1_pvalue")), "NA"),
          kegg_top3_txt
        ))
        figure_num <- figure_num + 1
      }
      doc <- doc %>% body_add_par("", style = "Normal")
    } else {
      doc <- doc %>% add_result_par("以 p.adjust < 0.05 为阈值，无显著 KEGG 通路。")
      doc <- doc %>% body_add_par("", style = "Normal")
    }
  } else {
    doc <- doc %>% add_result_par("无 KEGG 结果可用。")
    doc <- doc %>% body_add_par("", style = "Normal")
  }
}

# ---- GSEA Results ----
if (run_gsea) {
  # 确定表格编号
  table_num <- 1
  if (run_go_kegg) {
    table_num <- 3
  }

  doc <- doc %>% body_add_par(sprintf("2.%d GSEA 分析结果", result_section_num), style = "Normal")
  result_section_num <- result_section_num + 1

  if (!is.null(gsea_results) && nrow(gsea_results) > 0) {
    gsea_sig <- gsea_results[gsea_results$p.adjust < 0.05, ]

    if (nrow(gsea_sig) > 0) {
      # 分别统计上调和下调通路
      gsea_up <- gsea_sig[gsea_sig$NES > 0, ]
      gsea_down <- gsea_sig[gsea_sig$NES < 0, ]
      gsea_order <- gsea_sig[order(gsea_sig$p.adjust), ]
      top_gsea_term <- as.character(gsea_order$Description[1])
      top_gsea_p <- suppressWarnings(as.numeric(gsea_order$p.adjust[1]))
      top_gsea_nes <- suppressWarnings(as.numeric(gsea_order$NES[1]))

      # 生成图专属统计表（用于图后数字化分析）
      gsea_stats <- data.frame(
        metric = c("total_pathways", "up_pathways", "down_pathways", "top1_pathway", "top1_pvalue", "top1_nes"),
        value = c(
          nrow(gsea_sig),
          nrow(gsea_up),
          nrow(gsea_down),
          top_gsea_term,
          ifelse(is.na(top_gsea_p), NA_character_, format(top_gsea_p, scientific = TRUE, digits = 3)),
          ifelse(is.na(top_gsea_nes), NA_character_, format(round(top_gsea_nes, 4), nsmall = 4))
        ),
        stringsAsFactors = FALSE
      )
      fwrite(gsea_stats, file.path(gsea_dir, "gsea_figure_stats.csv"))

      # 汇总说明
      doc <- doc %>%
        add_result_par(sprintf("以 p.adjust < 0.05 且 |NES| > 1 为阈值，GSEA 分析共筛选出 %d 个显著通路，其中 %d 个通路在高表达基因集中富集（NES > 0，上调），%d 个通路在低表达基因集中富集（NES < 0，下调）。",
          nrow(gsea_sig), nrow(gsea_up), nrow(gsea_down)))

      # Top 通路
      top5_gsea <- gsea_sig[order(gsea_sig$p.adjust), ][1:min(5, nrow(gsea_sig)), ]
      top_pathways <- paste(top5_gsea$Description, collapse = "、")
      doc <- doc %>% add_result_par(sprintf("|NES| 排名前 5 的通路包括：%s。", top_pathways))
      doc <- doc %>% body_add_par("", style = "Normal")

      # Top 3 详细列表
      top3_gsea <- gsea_sig %>%
        arrange(desc(abs(NES))) %>%
        slice_head(n = 3)

      doc <- doc %>% add_result_par("显著富集的主要通路：")
      for (i in 1:nrow(top3_gsea)) {
        pathway <- top3_gsea[i, ]
        nes_val <- as.numeric(pathway$NES)
        direction <- ifelse(nes_val > 0, "上调", "下调")
        doc <- doc %>% add_result_par(sprintf("%d. %s (NES = %.2f, %s, p.adjust = %.2e)", i, pathway$Description, nes_val, direction, as.numeric(pathway$p.adjust)))
      }
      doc <- doc %>% body_add_par("", style = "Normal")

      # 创建表格
      top_gsea <- gsea_sig %>%
        arrange(p.adjust) %>%
        slice_head(n = 20)

      cols_map <- c(
        ID = "Pathway ID",
        Description = "Pathway",
        setSize = "Set Size",
        NES = "NES",
        p.adjust = "Adjusted P-value"
      )

      cols_exist <- names(cols_map)[names(cols_map) %in% colnames(top_gsea)]

      if (length(cols_exist) > 0) {
        table_data <- top_gsea[, cols_exist, drop = FALSE]
        colnames(table_data) <- cols_map[cols_exist]

        num_cols <- c("NES", "Adjusted P-value")
        for (col in num_cols) {
          if (col %in% colnames(table_data)) {
            table_data[[col]] <- round(as.numeric(table_data[[col]]), 4)
          }
        }

        ft <- flextable(table_data) %>%
          bold(part = "header", bold = TRUE) %>%
          fontsize(size = 9, part = "all") %>%
          font(fontname = "Times New Roman", part = "all") %>%
          align(align = "center", part = "all") %>%
          align(align = "left", part = "body", j = 2) %>%
          border_remove() %>%
          border(part = "all", border.top = officer::fp_border(width = 1.5, color = "black")) %>%
          border(part = "header", border.bottom = officer::fp_border(width = 0.5, color = "black")) %>%
          border(part = "body", border.bottom = officer::fp_border(width = 1.5, color = "black"), i = nrow(table_data)) %>%
          set_table_properties(layout = "autofit", width = 1)

        doc <- doc %>% add_center_par(sprintf("表%d：GSEA Results", table_num)) %>%
          body_add_flextable(ft, align = "center")
        top3_gsea_terms <- unique(as.character(gsea_order$Description[1:min(3, nrow(gsea_order))]))
        top3_gsea_txt <- if (length(top3_gsea_terms) > 0) paste(top3_gsea_terms, collapse = "、") else "NA"
        doc <- doc %>% add_result_par(sprintf(
          "表%d结果显示，共筛选到%d条显著GSEA通路（上调=%d，下调=%d）；最显著通路为%s（adjusted p-value=%s，NES=%s），核心富集通路包括%s。相关结果表：`prerank_gsea_results.csv`。",
          table_num, nrow(gsea_sig), nrow(gsea_up), nrow(gsea_down), top_gsea_term,
          ifelse(is.na(top_gsea_p), "NA", format(top_gsea_p, scientific = TRUE, digits = 3)),
          ifelse(is.na(top_gsea_nes), "NA", format(round(top_gsea_nes, 4), nsmall = 4)),
          top3_gsea_txt
        ))
      }
      doc <- doc %>% body_add_par("", style = "Normal")

      # GSEA Figure
      gsea_plot_png <- file.path(gsea_dir, "..", "plots", "PreRank_GSEA_EnrichmentPlot_SCI.png")
      gsea_plot_pdf <- file.path(gsea_dir, "..", "plots", "PreRank_GSEA_EnrichmentPlot_SCI.pdf")
      gsea_stats_file <- file.path(gsea_dir, "gsea_figure_stats.csv")
      gsea_stats <- if (file.exists(gsea_stats_file)) tryCatch(as.data.frame(fread(gsea_stats_file)), error = function(e) NULL) else NULL

      if (file.exists(gsea_plot_png)) {
        g_total <- nrow(gsea_sig)
        g_up <- nrow(gsea_up)
        g_down <- nrow(gsea_down)
        g_top <- top_gsea_term
        g_top_p <- ifelse(is.na(top_gsea_p), "NA", format(top_gsea_p, scientific = TRUE, digits = 3))
        g_top_nes <- ifelse(is.na(top_gsea_nes), "NA", format(round(top_gsea_nes, 4), nsmall = 4))
        if (!is.null(gsea_stats) && nrow(gsea_stats) > 0) {
          get_stat <- function(metric) gsea_stats$value[gsea_stats$metric == metric][1]
          g_total <- suppressWarnings(as.numeric(get_stat("total_pathways")))
          g_up <- suppressWarnings(as.numeric(get_stat("up_pathways")))
          g_down <- suppressWarnings(as.numeric(get_stat("down_pathways")))
          g_top <- as.character(get_stat("top1_pathway"))
          g_top_p <- as.character(get_stat("top1_pvalue"))
          g_top_nes <- as.character(get_stat("top1_nes"))
        }
        doc <- doc %>% add_center_img(gsea_plot_png, width = 6, height = 5) %>%
          add_center_par(sprintf("图%d：GSEA Enrichment Plot", figure_num)) %>%
          add_center_par("图注：X-axis为ranked gene list，Y-axis为running enrichment score，曲线峰值反映富集方向与强度。") %>%
          add_result_par(sprintf(
            "图%d结果显示，显著通路共%s条（上调=%s，下调=%s），最显著通路为%s（adjusted p-value=%s，NES=%s），提示该通路在排序基因集内富集强度最高。相关结果表：`gsea_figure_stats.csv`。",
            figure_num,
            ifelse(is.na(g_total), "NA", as.character(as.integer(g_total))),
            ifelse(is.na(g_up), "NA", as.character(as.integer(g_up))),
            ifelse(is.na(g_down), "NA", as.character(as.integer(g_down))),
            ifelse(is.na(g_top) || g_top == "", "NA", g_top),
            ifelse(is.na(g_top_p) || g_top_p == "", "NA", g_top_p),
            ifelse(is.na(g_top_nes) || g_top_nes == "", "NA", g_top_nes)
          ))
        figure_num <- figure_num + 1
      } else if (file.exists(gsea_plot_pdf)) {
        g_total <- nrow(gsea_sig)
        g_up <- nrow(gsea_up)
        g_down <- nrow(gsea_down)
        g_top <- top_gsea_term
        g_top_p <- ifelse(is.na(top_gsea_p), "NA", format(top_gsea_p, scientific = TRUE, digits = 3))
        g_top_nes <- ifelse(is.na(top_gsea_nes), "NA", format(round(top_gsea_nes, 4), nsmall = 4))
        if (!is.null(gsea_stats) && nrow(gsea_stats) > 0) {
          get_stat <- function(metric) gsea_stats$value[gsea_stats$metric == metric][1]
          g_total <- suppressWarnings(as.numeric(get_stat("total_pathways")))
          g_up <- suppressWarnings(as.numeric(get_stat("up_pathways")))
          g_down <- suppressWarnings(as.numeric(get_stat("down_pathways")))
          g_top <- as.character(get_stat("top1_pathway"))
          g_top_p <- as.character(get_stat("top1_pvalue"))
          g_top_nes <- as.character(get_stat("top1_nes"))
        }
        doc <- doc %>% add_center_img(gsea_plot_pdf, width = 6, height = 5) %>%
          add_center_par(sprintf("图%d：GSEA Enrichment Plot", figure_num)) %>%
          add_center_par("图注：X-axis为ranked gene list，Y-axis为running enrichment score，曲线峰值反映富集方向与强度。") %>%
          add_result_par(sprintf(
            "图%d结果显示，显著通路共%s条（上调=%s，下调=%s），最显著通路为%s（adjusted p-value=%s，NES=%s），提示该通路在排序基因集内富集强度最高。相关结果表：`gsea_figure_stats.csv`。",
            figure_num,
            ifelse(is.na(g_total), "NA", as.character(as.integer(g_total))),
            ifelse(is.na(g_up), "NA", as.character(as.integer(g_up))),
            ifelse(is.na(g_down), "NA", as.character(as.integer(g_down))),
            ifelse(is.na(g_top) || g_top == "", "NA", g_top),
            ifelse(is.na(g_top_p) || g_top_p == "", "NA", g_top_p),
            ifelse(is.na(g_top_nes) || g_top_nes == "", "NA", g_top_nes)
          ))
        figure_num <- figure_num + 1
      }
      doc <- doc %>% body_add_par("", style = "Normal")
    } else {
      sig_pathways <- NA_real_
      if (!is.null(qc_results) && nrow(qc_results) > 0 && all(c("Item", "Value") %in% colnames(qc_results))) {
        sig_pathways <- suppressWarnings(as.numeric(qc_results$Value[qc_results$Item == "Significant_Pathways"][1]))
      }
      doc <- doc %>% add_result_par(sprintf(
        "以 p.adjust < 0.05 且 |NES| > 1 为阈值，显著GSEA通路数为%s，当前未达到显著富集标准。相关结果表：`PreRank_GSEA_QC_Report.csv`。",
        ifelse(is.na(sig_pathways), "0", as.character(as.integer(sig_pathways)))
      ))
      doc <- doc %>% body_add_par("", style = "Normal")
    }
  } else {
    doc <- doc %>% add_result_par("无 GSEA 结果可用。")
    doc <- doc %>% body_add_par("", style = "Normal")
  }
}  # End of if (run_gsea)

# ========== Save document ==========
print(doc, target = opt$output)
postprocess_docx_fonts(opt$output)

cat("\n================================================================================\n")
cat("报告已生成：", opt$output, "\n")
cat("================================================================================\n")

# ============================================================================
# 多数据集合并与批次效应校正脚本 (Combine Standalone)
# ============================================================================

invisible(gc())

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(sva)
  library(RColorBrewer)
})

options(scipen = 5)
set.seed(123)

save_plot <- function(filename, plot_object, outdir = ".",
                      width = 7, height = 7, both = FALSE,
                      bg = "white", family = "sans",
                      units = "in", res = 300) {
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  }

  is_ggplot <- inherits(plot_object, "ggplot") || inherits(plot_object, "gg")
  draw_plot <- function() {
    if (is_ggplot) {
      print(plot_object)
    } else if (is.function(plot_object)) {
      plot_object()
    } else {
      print(plot_object)
    }
  }

  base_name <- tools::file_path_sans_ext(filename)
  pdf_path <- file.path(outdir, paste0(base_name, ".pdf"))
  png_path <- file.path(outdir, paste0(base_name, ".png"))

  if (both) {
    grDevices::cairo_pdf(file = pdf_path, width = width, height = height,
                         bg = bg, family = family)
    draw_plot()
    dev.off()

    png(filename = png_path, width = width, height = height,
        bg = bg, units = units, res = res, type = "cairo", family = family)
    draw_plot()
    dev.off()
    return(invisible(c(pdf = pdf_path, png = png_path)))
  }

  ext <- tolower(tools::file_ext(filename))
  target_path <- file.path(outdir, filename)
  if (ext == "pdf") {
    grDevices::cairo_pdf(file = target_path, width = width, height = height,
                         bg = bg, family = family)
    draw_plot()
    dev.off()
  } else if (ext == "png") {
    png(filename = target_path, width = width, height = height,
        bg = bg, units = units, res = res, type = "cairo", family = family)
    draw_plot()
    dev.off()
  } else {
    stop("不支持的文件格式，请使用 .pdf 或 .png")
  }

  invisible(target_path)
}

split_csv_arg <- function(value) {
  if (is.null(value) || !nzchar(trimws(value))) {
    return(character(0))
  }
  trimws(unlist(strsplit(value, ",", fixed = TRUE)))
}

get_dataset_colors <- function(n) {
  if (n <= 8) {
    brewer.pal(max(n, 3), "Set2")[1:n]
  } else {
    colorRampPalette(brewer.pal(8, "Set2"))(n)
  }
}

group_colors <- c(
  "firebrick2", "#386CB0", "orange", "seagreen",
  "#BC80BD", "#17BECF", "#E377C2", "#7F7F7F"
)

plot_pca <- function(expr_matrix, group_df, title = "PCA Plot",
                     color_by = "group", label_samples = TRUE) {
  pca_result <- prcomp(t(expr_matrix), scale. = TRUE)
  pca_data <- as.data.frame(pca_result$x)
  pca_data$sample <- rownames(pca_data)

  plot_df <- merge(group_df, pca_data, by = "sample")
  var_explained <- summary(pca_result)$importance[2, ]
  pc1_label <- sprintf("PC1 (%.2f%%)", var_explained[1] * 100)
  pc2_label <- sprintf("PC2 (%.2f%%)", var_explained[2] * 100)

  color_var <- if (color_by %in% colnames(plot_df)) color_by else "group"
  n_colors <- length(unique(plot_df[[color_var]]))
  color_values <- if (color_var == "group") {
    group_colors[1:min(n_colors, length(group_colors))]
  } else {
    get_dataset_colors(n_colors)
  }

  plot_obj <- ggplot(plot_df, aes(x = PC1, y = PC2, color = .data[[color_var]])) +
    geom_point(size = 3, alpha = 0.85) +
    stat_ellipse(type = "norm", level = 0.95, linetype = "dashed", linewidth = 0.5) +
    scale_color_manual(values = color_values) +
    labs(title = title, x = pc1_label, y = pc2_label, color = color_var) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )

  if (label_samples) {
    plot_obj <- plot_obj +
      ggrepel::geom_text_repel(aes(label = sample), size = 3, max.overlaps = 30, show.legend = FALSE)
  }

  list(plot = plot_obj, pca_data = pca_data, pca_result = pca_result)
}

option_list <- list(
  make_option(c("-d", "--data-files"), dest = "data_files", type = "character",
              help = "表达矩阵 CSV 文件路径，逗号分隔（至少 2 个）。第一列须为 SYMBOL"),
  make_option(c("-g", "--group-files"), dest = "group_files", type = "character",
              help = "分组文件 CSV 路径，逗号分隔，与 --data-files 一一对应。须含 sample 和 group 列"),
  make_option(c("-n", "--dataset-names"), dest = "dataset_names", type = "character",
              default = NULL,
              help = "数据集名称，逗号分隔。若不指定，从文件名自动推断"),
  make_option(c("-o", "--output-dir"), dest = "output_dir", type = "character",
              default = "./output",
              help = "输出目录 [default: %default]"),
  make_option(c("--output-name"), dest = "output_name", type = "character",
              default = "combat",
              help = "输出文件前缀名 [default: %default]"),
  make_option(c("--gene-col"), dest = "gene_col", type = "character",
              default = "SYMBOL",
              help = "表达矩阵中基因列名 [default: %default]"),
  make_option(c("--no-label"), dest = "no_label", action = "store_true",
              default = FALSE,
              help = "PCA 图不标注样本名"),
  make_option(c("-w", "--width"), dest = "width", type = "double",
              default = 9,
              help = "图片宽度（英寸）[default: %default]"),
  make_option(c("--height"), dest = "height", type = "double",
              default = 7,
              help = "图片高度（英寸）[default: %default]")
)

parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser, convert_hyphens_to_underscores = TRUE)

if (is.null(opt$data_files) || !nzchar(trimws(opt$data_files))) {
  print_help(parser)
  stop("缺少必需参数: --data-files")
}
if (is.null(opt$group_files) || !nzchar(trimws(opt$group_files))) {
  print_help(parser)
  stop("缺少必需参数: --group-files")
}

data_files <- split_csv_arg(opt$data_files)
group_files <- split_csv_arg(opt$group_files)

if (length(data_files) < 2) {
  stop("至少需要 2 个表达矩阵文件")
}
if (length(data_files) != length(group_files)) {
  stop(sprintf("--data-files (%d) 与 --group-files (%d) 数量不一致",
               length(data_files), length(group_files)))
}

for (file_path in c(data_files, group_files)) {
  if (!file.exists(file_path)) {
    stop(sprintf("文件不存在: %s", file_path))
  }
}

dataset_names <- split_csv_arg(opt$dataset_names)
if (length(dataset_names) == 0) {
  dataset_names <- tools::file_path_sans_ext(basename(data_files))
}
if (length(dataset_names) != length(data_files)) {
  stop("--dataset-names 数量与 --data-files 不一致")
}

gene_col <- opt$gene_col
output_name <- opt$output_name

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_dir <- normalizePath(opt$output_dir, winslash = "/", mustWork = FALSE)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}
run_dir <- file.path(output_dir, paste0("combine_", timestamp))
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
tracking_file <- file.path(run_dir, paste0(".", timestamp))
invisible(file.create(tracking_file))

logs_dir <- file.path(run_dir, "logs")
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(logs_dir, paste0("combine_", timestamp, ".log"))
log_con <- file(log_file, open = "wt")
sink(log_con, split = TRUE)
on.exit({
  while (sink.number() > 0) {
    sink()
  }
  close(log_con)
}, add = TRUE)

log_message <- function(msg, level = "INFO") {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[%s] [%s] %s\n", ts, level, msg))
}

cat("============================================================\n")
cat(" 多数据集合并与批次效应校正\n")
cat("============================================================\n")
cat(sprintf("数据集数量: %d\n", length(data_files)))
cat(sprintf("数据集名称: %s\n", paste(dataset_names, collapse = ", ")))
cat(sprintf("输出目录: %s\n", run_dir))
cat(sprintf("日志文件: %s\n", log_file))
cat(sprintf("追踪文件: %s\n", tracking_file))
cat("------------------------------------------------------------\n")

log_message("[步骤 1] 读取表达矩阵...")
data_list <- lapply(seq_along(data_files), function(index) {
  cat(sprintf("  读取 %s <- %s\n", dataset_names[index], data_files[index]))
  data_frame <- fread(data_files[index], data.table = FALSE)

  active_gene_col <- gene_col
  if (!(active_gene_col %in% colnames(data_frame))) {
    alt_cols <- intersect(c("SYMBOL", "gene", "Gene", "GENE", "GeneID", "gene_name"),
                          colnames(data_frame))
    if (length(alt_cols) > 0) {
      active_gene_col <- alt_cols[1]
      cat(sprintf("    注意: 未找到列 '%s'，改用 '%s'\n", gene_col, active_gene_col))
    } else {
      stop(sprintf("文件 %s 中找不到基因列 '%s'", data_files[index], gene_col))
    }
  }

  rownames(data_frame) <- data_frame[[active_gene_col]]
  data_frame[[active_gene_col]] <- NULL

  numeric_cols <- vapply(data_frame, is.numeric, logical(1))
  data_frame <- data_frame[, numeric_cols, drop = FALSE]
  if (ncol(data_frame) == 0) {
    stop(sprintf("文件 %s 中没有数值型样本列", data_files[index]))
  }

  data_frame
})
names(data_list) <- dataset_names

log_message("[步骤 2] 读取分组文件...")
group_list <- lapply(seq_along(group_files), function(index) {
  cat(sprintf("  读取 %s <- %s\n", dataset_names[index], group_files[index]))
  group_df <- fread(group_files[index], data.table = FALSE)
  if (!all(c("sample", "group") %in% colnames(group_df))) {
    stop(sprintf("分组文件 %s 必须包含 'sample' 和 'group' 列", group_files[index]))
  }

  sample_order <- colnames(data_list[[index]])
  missing_samples <- setdiff(sample_order, group_df$sample)
  if (length(missing_samples) > 0) {
    stop(sprintf("分组文件 %s 缺少样本: %s",
                 group_files[index], paste(missing_samples, collapse = ", ")))
  }

  group_df <- group_df[match(sample_order, group_df$sample), c("sample", "group"), drop = FALSE]
  group_df$dataset <- dataset_names[index]
  group_df
})

log_message("[步骤 3] 提取共同基因...")
common_genes <- Reduce(intersect, lapply(data_list, rownames))
cat(sprintf("  共同基因数: %d\n", length(common_genes)))
if (length(common_genes) == 0) {
  stop("数据集之间没有共同基因，请检查基因列名和数据")
}

gene_info <- data.table(
  dataset = dataset_names,
  n_genes = vapply(data_list, nrow, integer(1)),
  n_samples = vapply(data_list, ncol, integer(1))
)
gene_info[, loss_pct := round((n_genes - length(common_genes)) / n_genes * 100, 1)]
gene_info <- rbind(
  gene_info,
  data.table(dataset = "共同基因", n_genes = length(common_genes), n_samples = NA_integer_, loss_pct = 0)
)
fwrite(gene_info, file.path(run_dir, "00.dataset_info.csv"))
cat(sprintf("  基因丢失百分比: %s\n", paste(gene_info$dataset[-nrow(gene_info)], gene_info$loss_pct[-nrow(gene_info)], sep = "=", collapse = "%, ")))

data_list <- lapply(data_list, function(data_frame) data_frame[common_genes, , drop = FALSE])
combined_data <- do.call(cbind, unname(data_list))
cat(sprintf("  合并矩阵: %d 基因 x %d 样本\n", nrow(combined_data), ncol(combined_data)))

group_df <- do.call(rbind, group_list)
group_df <- as.data.frame(group_df, stringsAsFactors = FALSE)
if (!identical(colnames(combined_data), group_df$sample)) {
  stop("合并后的表达矩阵列顺序与分组文件样本顺序不一致")
}

label_samples <- !isTRUE(opt$no_label)
dataset_colors <- get_dataset_colors(length(dataset_names))
batch_colors <- dataset_colors[rep(seq_along(dataset_names), vapply(data_list, ncol, integer(1)))]

log_message("[步骤 4] 校正前 PCA...")
pca_before_group <- plot_pca(combined_data, group_df,
  title = "PCA Before Batch Correction (by Group)",
  color_by = "group",
  label_samples = label_samples)
save_plot("01.pca_before_group.pdf", pca_before_group$plot,
  outdir = run_dir, width = opt$width, height = opt$height, both = TRUE)

pca_before_dataset <- plot_pca(combined_data, group_df,
  title = "PCA Before Batch Correction (by Dataset)",
  color_by = "dataset",
  label_samples = label_samples)
save_plot("02.pca_before_dataset.pdf", pca_before_dataset$plot,
  outdir = run_dir, width = opt$width, height = opt$height, both = TRUE)

fwrite(pca_before_group$pca_data, file.path(run_dir, "01.pca_before.csv"), row.names = TRUE)

boxplot_before_fn <- function() {
  par(mar = c(8, 4, 3, 1), las = 2, cex.axis = 0.6)
  boxplot(combined_data,
    main = "Expression Distribution (Before Correction)",
    col = batch_colors,
    outline = FALSE)
}
save_plot("03.boxplot_before.pdf", boxplot_before_fn,
  outdir = run_dir,
  width = max(10, ncol(combined_data) * 0.3),
  height = 6,
  both = TRUE)

log_message("[步骤 5] ComBat 批次效应校正...")
batch <- factor(rep(seq_along(dataset_names), vapply(data_list, ncol, integer(1))))
cat(sprintf("  批次分布: %s\n",
            paste(sprintf("%s(n=%d)", dataset_names, vapply(data_list, ncol, integer(1))), collapse = ", ")))
combat_data <- ComBat(dat = as.matrix(combined_data), batch = batch, mod = NULL)
combat_data <- as.data.frame(combat_data, check.names = FALSE)
cat("  ComBat 校正完成\n")

log_message("[步骤 6] 校正后 PCA...")
pca_after_group <- plot_pca(combat_data, group_df,
  title = "PCA After Batch Correction (by Group)",
  color_by = "group",
  label_samples = label_samples)
save_plot("04.pca_after_group.pdf", pca_after_group$plot,
  outdir = run_dir, width = opt$width, height = opt$height, both = TRUE)

pca_after_dataset <- plot_pca(combat_data, group_df,
  title = "PCA After Batch Correction (by Dataset)",
  color_by = "dataset",
  label_samples = label_samples)
save_plot("05.pca_after_dataset.pdf", pca_after_dataset$plot,
  outdir = run_dir, width = opt$width, height = opt$height, both = TRUE)

fwrite(pca_after_group$pca_data, file.path(run_dir, "04.pca_after.csv"), row.names = TRUE)

boxplot_after_fn <- function() {
  par(mar = c(8, 4, 3, 1), las = 2, cex.axis = 0.6)
  boxplot(combat_data,
    main = "Expression Distribution (After Correction)",
    col = batch_colors,
    outline = FALSE)
}
save_plot("06.boxplot_after.pdf", boxplot_after_fn,
  outdir = run_dir,
  width = max(10, ncol(combat_data) * 0.3),
  height = 6,
  both = TRUE)

log_message("[步骤 7] 输出校正后数据...")
group_sorted <- arrange(group_df, group)
combat_out <- combat_data[, group_sorted$sample, drop = FALSE]
combat_out <- cbind(data.frame(SYMBOL = rownames(combat_out)), combat_out)

fwrite(combat_out, file.path(run_dir, paste0(output_name, ".dat.csv")))
fwrite(group_sorted, file.path(run_dir, paste0(output_name, ".group.csv")))

cat(sprintf("  表达矩阵: %s.dat.csv (%d x %d)\n", output_name, nrow(combat_out), ncol(combat_out) - 1))
cat(sprintf("  分组文件: %s.group.csv (%d 样本)\n", output_name, nrow(group_sorted)))
log_message(paste("追踪文件已创建:", tracking_file))
log_message("分析完成")

cat("\n============================================================\n")
cat(" 分析完成!\n")
cat(sprintf(" 输出目录: %s\n", run_dir))
cat("============================================================\n")

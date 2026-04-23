# ==============================================================================
# unified_enrichment/requirements.R
# 依赖包安装和加载脚本
# ==============================================================================

#' 检查并安装Bioconductor包
#' @param packages 包名向量
install_bioc_packages <- function(packages) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN")
  }

  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      cat(sprintf("安装Bioconductor包: %s\n", pkg))
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    }
  }
}

#' 检查并安装CRAN包
#' @param packages 包名向量
install_cran_packages <- function(packages) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      cat(sprintf("安装CRAN包: %s\n", pkg))
      install.packages(pkg, repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN")
    }
  }
}

#' 安装所有依赖包
install_dependencies <- function() {
  cat("检查并安装依赖包...\n")

  # CRAN包
  cran_packages <- c(
    # 数据处理
    "yaml",
    "data.table",
    "dplyr",
    "tidyr",
    "stringr",
    "glue",
    # 绘图
    "ggplot2",
    "ggsci",
    "viridis",
    "viridisLite",
    "patchwork",
    "RColorBrewer",
    # 报告生成
    "officer",
    "flextable",
    "rmarkdown",
    "knitr",
    "kableExtra",
    # 工具
    "magick"
  )
  install_cran_packages(cran_packages)

  # Bioconductor包
  bioc_packages <- c(
    # 核心富集分析
    "clusterProfiler",
    "enrichplot",
    "fgsea",
    # 基因注释
    "org.Hs.eg.db",
    "org.Mm.eg.db",
    "org.Rn.eg.db",
    "AnnotationDbi",
    # KEGG
    "KEGGREST"
  )
  install_bioc_packages(bioc_packages)

  cat("依赖包检查完成\n")
}

#' 加载所有必需包
load_packages <- function() {
  suppressPackageStartupMessages({
    # 数据处理
    library(yaml)
    library(data.table)
    library(dplyr)
    library(tidyr)
    library(stringr)

    # 绘图
    library(ggplot2)
    library(ggsci)
    library(viridis)
    library(viridisLite)
    library(patchwork)

    # 富集分析
    library(clusterProfiler)
    library(enrichplot)
    library(fgsea)

    # 报告
    library(officer)
    library(flextable)
    library(rmarkdown)
  })
}

# 如果直接运行此脚本，则安装依赖
if (sys.nframe() == 0) {
  install_dependencies()
  cat("\n所有依赖已安装。运行以下命令开始分析:\n")
  cat("  Rscript scripts/main.R\n")
}

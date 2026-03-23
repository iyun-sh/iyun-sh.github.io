#!/bin/bash
# ==============================================================================
# unified_enrichment/run_unified_enrichment.sh
# 统一富集分析入口脚本
# 支持 GO/KEGG/GSEA 三种分析的灵活组合执行
# 支持从配置文件加载并允许命令行覆盖
# ==============================================================================

set -e

# 脚本目录（必须尽早定义，供权限切换和默认配置使用）
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# 目标用户配置（动态检测脚本所有者）
TARGET_USER=$(stat -c "%U" "$0")

# 保存完整参数（用于 sudo 传递）
FULL_ARGS=("$@")

# 默认配置文件（硬约束：启动即要求存在）
DEFAULT_CONFIG_FILE="${SCRIPT_DIR}/scripts/config/unified_enrichment.config.ini"
CONFIG_FILE="$DEFAULT_CONFIG_FILE"

# 默认值
WORK_DIR=""
SPECIES=""
GENE_LIST=""
EXPR_MATRIX=""
TARGET_GENES=""
RANKED_LIST=""
GMT_FILE=""

# 分析开关（默认值）
ANALYSIS_MODE="full"
GO_KEGG="true"
GSEA="true"
GSEA_METHOD="prerank"

# GO/KEGG参数
GO_ONT="ALL"
GO_PVALUE="0.05"
GO_QVALUE="0.2"
GO_PADJ="BH"

# GSEA参数
GSEA_CORRELATION="spearman"
GSEA_RANK_METRIC="log2FC"
GSEA_MIN_SIZE="10"
GSEA_MAX_SIZE="500"
GSEA_PVALUE="0.05"
GSEA_PADJ="BH"
GSEA_EXPONENT="1"
GSEA_EPS="0"
GSEA_TOP_N="5"

# QC参数
QC_VARIANCE="0.01"
QC_NA_RATIO="0.5"
QC_GMT_MATCH="30"

# 绘图参数
PLOT_WIDTH="10"
PLOT_HEIGHT="8"
PLOT_TOP_N="5"
PLOT_COLOR="viridis"

# 其他参数
OUTPUT_DIR="results"
LOG_DIR="logs"
SEED="123"

# 权限控制参数
OUTPUT_PATH=""      # 实际输出路径（用于设置权限）
DO_RM=false         # 删除模式

# 从参数中提取 --config 路径（用于在正式解析前确定配置文件）
get_config_arg() {
    local args=("$@")
    local i=0
    while [[ $i -lt ${#args[@]} ]]; do
        case "${args[$i]}" in
            -c|--config)
                if [[ $((i+1)) -lt ${#args[@]} ]]; then
                    echo "${args[$((i+1))]}"
                    return 0
                fi
                ;;
        esac
        i=$((i+1))
    done
    return 1
}

# 读取 config（key=value；注释和空行会被忽略）
load_config_file() {
    local cfg="$1"
    [[ -f "$cfg" ]] || return 1

    while IFS='=' read -r raw_key raw_val; do
        key="$(echo "$raw_key" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')"
        val="$(echo "$raw_val" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')"
        val="${val%\"}"
        val="${val#\"}"

        [[ -z "$key" ]] && continue
        [[ "$key" =~ ^# ]] && continue
        [[ "$key" =~ ^\; ]] && continue
        [[ "$key" =~ ^\[.*\]$ ]] && continue

        case "$key" in
            work_dir) WORK_DIR="$val" ;;
            species) SPECIES="$val" ;;
            gene_list) GENE_LIST="$val" ;;
            expr_matrix) EXPR_MATRIX="$val" ;;
            target_genes) TARGET_GENES="$val" ;;
            ranked_list) RANKED_LIST="$val" ;;
            gmt_file) GMT_FILE="$val" ;;
            analysis_mode) ANALYSIS_MODE="$val" ;;
            go_kegg) GO_KEGG="$val" ;;
            gsea) GSEA="$val" ;;
            gsea_method) GSEA_METHOD="$val" ;;
            go_ont) GO_ONT="$val" ;;
            go_pvalue) GO_PVALUE="$val" ;;
            go_qvalue) GO_QVALUE="$val" ;;
            go_padj) GO_PADJ="$val" ;;
            gsea_correlation) GSEA_CORRELATION="$val" ;;
            gsea_rank_metric) GSEA_RANK_METRIC="$val" ;;
            gsea_min_size) GSEA_MIN_SIZE="$val" ;;
            gsea_max_size) GSEA_MAX_SIZE="$val" ;;
            gsea_pvalue) GSEA_PVALUE="$val" ;;
            gsea_padj) GSEA_PADJ="$val" ;;
            gsea_exponent) GSEA_EXPONENT="$val" ;;
            gsea_eps) GSEA_EPS="$val" ;;
            gsea_top_n) GSEA_TOP_N="$val" ;;
            qc_variance) QC_VARIANCE="$val" ;;
            qc_na_ratio) QC_NA_RATIO="$val" ;;
            qc_gmt_match) QC_GMT_MATCH="$val" ;;
            plot_width) PLOT_WIDTH="$val" ;;
            plot_height) PLOT_HEIGHT="$val" ;;
            plot_top_n) PLOT_TOP_N="$val" ;;
            plot_color) PLOT_COLOR="$val" ;;
            output_dir) OUTPUT_DIR="$val" ;;
            log_dir) LOG_DIR="$val" ;;
            seed) SEED="$val" ;;
        esac
    done < "$cfg"
}

# 显示用法
show_usage() {
    cat << 'EOF'
用法: run_unified_enrichment.sh [必需参数] [可选参数]

【必需参数】
  -w, --work-dir DIR         工作目录
  -s, --species CODE         物种代码 (hsa/mmu/rno)
  --gene-list FILE           基因列表文件
  --expr-matrix FILE         表达矩阵文件
  --target-genes FILE        目标基因文件
  --ranked-list FILE         排序基因列表文件
  --gmt FILE                 GMT基因集文件

【分析开关】
  --analysis-mode MODE       分析模式 (full/gsea-only, 默认: full)
  --go-kegg BOOL             GO/KEGG分析 (true/false, 默认: true)
  --gsea BOOL                GSEA分析 (true/false, 默认: true)
  --gsea-method METHOD       GSEA方法 (prerank/correlation/median_split, 默认: prerank)

【GO/KEGG参数】
  --go-ont ONT               GO本体 (BP/MF/CC/ALL, 默认: ALL)
  --go-pvalue VALUE          p值阈值 (默认: 0.05)
  --go-qvalue VALUE          q值阈值 (默认: 0.2)
  --go-padj METHOD           校正方法 (BH/none, 默认: BH)

【GSEA参数】
  --gsea-correlation METHOD  相关方法 (spearman/pearson, 默认: spearman)
  --gsea-rank-metric COL     排序指标列名 (默认: log2FC)
  --gsea-min-size N          最小基因集大小 (默认: 10)
  --gsea-max-size N          最大基因集大小 (默认: 500)
  --gsea-pvalue VALUE        p值阈值 (默认: 0.05)
  --gsea-padj METHOD         校正方法 (默认: BH)
  --gsea-exponent VALUE      指数 (默认: 1)
  --gsea-eps VALUE           eps值 (默认: 0)
  --gsea-top-n N             top通路数 (默认: 5)

【QC参数】
  --qc-variance VALUE        方差阈值 (默认: 0.01)
  --qc-na-ratio VALUE        NA比例阈值 (默认: 0.5)
  --qc-gmt-match RATE        GMT匹配率阈值 (默认: 30)

【绘图参数】
  --plot-width VALUE         图宽度 (默认: 10)
  --plot-height VALUE        图高度 (默认: 8)
  --plot-top-n N             top条目数 (默认: 5)
  --plot-color SCHEME        配色方案 (默认: viridis)

【其他参数】
  -c, --config FILE          配置文件路径（默认: scripts/config/unified_enrichment.config.ini）
  --output-dir DIR           输出目录 (默认: results)
  --log-dir DIR              日志目录 (默认: logs)
  --seed VALUE               随机种子 (默认: 123)
  --rm                       删除输出目录 (仅创建者可执行)
  -h, --help                 显示帮助信息

【示例】
  # 最小用法
  ./run_unified_enrichment.sh \\
    -w /path/to/project \\
    -s hsa \\
    --gene-list input/gene_list.csv \\
    --expr-matrix input/expression_matrix.csv \\
    --target-genes input/target_genes.csv \\
    --ranked-list input/ranked_gene_list.csv \\
    --gmt input/gene_sets.gmt

  # 自定义参数
  ./run_unified_enrichment.sh \\
    -w /path/to/project \\
    -s hsa \\
    --gene-list input/gene_list.csv \\
    --expr-matrix input/expression_matrix.csv \\
    --target-genes input/target_genes.csv \\
    --ranked-list input/ranked_gene_list.csv \\
    --gmt input/gene_sets.gmt \\
    --gsea false \\
    --go-ont BP \\
    --go-pvalue 0.01

  # 删除输出目录
  ./run_unified_enrichment.sh --rm -w /path/to/project --output-dir /path/to/output
EOF
}

# 如果不是目标用户，通过 sudo 重新运行
if [ "$(whoami)" != "$TARGET_USER" ]; then
    exec sudo -u "$TARGET_USER" "$SCRIPT_DIR/run_unified_enrichment.sh" "${FULL_ARGS[@]}"
fi

# 解析命令行参数
parse_args() {
    while [[ $# -gt 0 ]]; do
        case "$1" in
            -w|--work-dir)
                WORK_DIR="$2"
                shift 2
                ;;
            -s|--species)
                SPECIES="$2"
                shift 2
                ;;
            --gene-list)
                GENE_LIST="$2"
                shift 2
                ;;
            --expr-matrix)
                EXPR_MATRIX="$2"
                shift 2
                ;;
            --target-genes)
                TARGET_GENES="$2"
                shift 2
                ;;
            --ranked-list)
                RANKED_LIST="$2"
                shift 2
                ;;
            --gmt)
                GMT_FILE="$2"
                shift 2
                ;;
            --analysis-mode)
                ANALYSIS_MODE="$2"
                shift 2
                ;;
            --go-kegg)
                GO_KEGG="$2"
                shift 2
                ;;
            --gsea)
                GSEA="$2"
                shift 2
                ;;
            --gsea-method)
                GSEA_METHOD="$2"
                shift 2
                ;;
            --go-ont)
                GO_ONT="$2"
                shift 2
                ;;
            --go-pvalue)
                GO_PVALUE="$2"
                shift 2
                ;;
            --go-qvalue)
                GO_QVALUE="$2"
                shift 2
                ;;
            --go-padj)
                GO_PADJ="$2"
                shift 2
                ;;
            --gsea-correlation)
                GSEA_CORRELATION="$2"
                shift 2
                ;;
            --gsea-rank-metric)
                GSEA_RANK_METRIC="$2"
                shift 2
                ;;
            --gsea-min-size)
                GSEA_MIN_SIZE="$2"
                shift 2
                ;;
            --gsea-max-size)
                GSEA_MAX_SIZE="$2"
                shift 2
                ;;
            --gsea-pvalue)
                GSEA_PVALUE="$2"
                shift 2
                ;;
            --gsea-padj)
                GSEA_PADJ="$2"
                shift 2
                ;;
            --gsea-exponent)
                GSEA_EXPONENT="$2"
                shift 2
                ;;
            --gsea-eps)
                GSEA_EPS="$2"
                shift 2
                ;;
            --gsea-top-n)
                GSEA_TOP_N="$2"
                shift 2
                ;;
            --qc-variance)
                QC_VARIANCE="$2"
                shift 2
                ;;
            --qc-na-ratio)
                QC_NA_RATIO="$2"
                shift 2
                ;;
            --qc-gmt-match)
                QC_GMT_MATCH="$2"
                shift 2
                ;;
            --plot-width)
                PLOT_WIDTH="$2"
                shift 2
                ;;
            --plot-height)
                PLOT_HEIGHT="$2"
                shift 2
                ;;
            --plot-top-n)
                PLOT_TOP_N="$2"
                shift 2
                ;;
            --plot-color)
                PLOT_COLOR="$2"
                shift 2
                ;;
            --output-dir)
                OUTPUT_DIR="$2"
                shift 2
                ;;
            --log-dir)
                LOG_DIR="$2"
                shift 2
                ;;
            --seed)
                SEED="$2"
                shift 2
                ;;
            -c|--config)
                CONFIG_FILE="$2"
                shift 2
                ;;
            --rm)
                DO_RM=true
                shift
                ;;
            -h|--help)
                show_usage
                exit 0
                ;;
            *)
                echo "错误: 未知参数 '$1'"
                show_usage
                exit 1
                ;;
        esac
    done
}

# 验证必需参数
validate_args() {
    local missing=()

    # 硬约束：启动时必须有 config 文件
    if [[ -z "$CONFIG_FILE" || ! -f "$CONFIG_FILE" ]]; then
        echo "错误: 缺少配置文件。必须提供并存在: --config FILE"
        echo "当前路径: ${CONFIG_FILE}"
        echo "默认路径: ${DEFAULT_CONFIG_FILE}"
        exit 1
    fi

    if [[ -z "$WORK_DIR" ]]; then
        missing+=("-w/--work-dir")
    fi
    if [[ -z "$SPECIES" ]]; then
        missing+=("-s/--species")
    fi

    # 根据分析模式验证不同参数
    if [[ "$ANALYSIS_MODE" != "gsea-only" ]]; then
        # 完整模式需要这些参数
        if [[ -z "$GENE_LIST" ]]; then
            missing+=("--gene-list")
        fi
        if [[ -z "$EXPR_MATRIX" ]]; then
            missing+=("--expr-matrix")
        fi
        if [[ -z "$TARGET_GENES" ]]; then
            missing+=("--target-genes")
        fi
    fi

    # GSEA 需要这些参数
    if [[ "$GSEA" == "true" ]]; then
        if [[ -z "$RANKED_LIST" ]]; then
            missing+=("--ranked-list")
        fi
        if [[ -z "$GMT_FILE" ]]; then
            missing+=("--gmt")
        fi
    fi

    if [[ ${#missing[@]} -gt 0 ]]; then
        echo "错误: 缺少必需参数:"
        for param in "${missing[@]}"; do
            echo "  - $param"
        done
        if [[ "$ANALYSIS_MODE" == "gsea-only" ]]; then
            echo ""
            echo "提示: 使用 --analysis-mode gsea-only 可仅运行GSEA分析"
        fi
        echo ""
        show_usage
        exit 1
    fi

    # 验证物种代码
    if [[ ! "$SPECIES" =~ ^(hsa|mmu|rno)$ ]]; then
        echo "错误: 不支持的物种 '$SPECIES'。支持: hsa(人类), mmu(小鼠), rno(大鼠)"
        exit 1
    fi

    # 验证布尔值
    if [[ ! "$GO_KEGG" =~ ^(true|false)$ ]]; then
        echo "错误: --go-kegg 必须是 true 或 false"
        exit 1
    fi
    if [[ ! "$GSEA" =~ ^(true|false)$ ]]; then
        echo "错误: --gsea 必须是 true 或 false"
        exit 1
    fi

    # 验证GSEA方法
    if [[ ! "$GSEA_METHOD" =~ ^(prerank|correlation|median_split)$ ]]; then
        echo "错误: --gsea-method 必须是 prerank, correlation 或 median_split"
        exit 1
    fi

    # 验证GO本体
    if [[ ! "$GO_ONT" =~ ^(BP|MF|CC|ALL)$ ]]; then
        echo "错误: --go-ont 必须是 BP, MF, CC 或 ALL"
        exit 1
    fi
}

# 主函数
main() {
    # 先识别 --config，确保配置在命令行覆盖前先加载
    cfg_from_args="$(get_config_arg "$@" || true)"
    if [[ -n "$cfg_from_args" ]]; then
        CONFIG_FILE="$cfg_from_args"
    fi
    if [[ -z "$CONFIG_FILE" || ! -f "$CONFIG_FILE" ]]; then
        echo "错误: 配置文件不存在: ${CONFIG_FILE}"
        echo "请先创建配置文件，或通过 --config 指定。"
        exit 1
    fi
    load_config_file "$CONFIG_FILE"

    # 再解析命令行，命令行优先级最高
    parse_args "$@"
    validate_args

    # 处理删除模式
    if [ "$DO_RM" = true ]; then
        if [ -z "$WORK_DIR" ]; then
            echo "错误：删除模式需要指定 -w/--work-dir" >&2
            exit 1
        fi
        # 确定输出目录路径
        if [ -z "$OUTPUT_DIR" ] || [ "$OUTPUT_DIR" = "results" ]; then
            DELETE_DIR="$WORK_DIR/$OUTPUT_DIR"
        else
            DELETE_DIR="$OUTPUT_DIR"
        fi
        if [ ! -d "$DELETE_DIR" ]; then
            echo "错误：输出目录不存在: $DELETE_DIR" >&2
            exit 1
        fi
        # 检查 .owner 文件
        if [ -f "$DELETE_DIR/.owner" ]; then
            OWNER=$(cat "$DELETE_DIR/.owner")
            CALLER=${SUDO_USER:-$USER}
            if [ "$OWNER" != "$CALLER" ]; then
                echo "错误：只能删除自己创建的目录（创建者: $OWNER，当前用户: $CALLER）" >&2
                exit 1
            fi
        fi
        rm -rf "$DELETE_DIR"
        echo "已删除: $DELETE_DIR"
        exit 0
    fi

    # 设置环境变量供R脚本使用
    export MODULE_ROOT="$SCRIPT_DIR"

    # 切换到工作目录
    cd "$WORK_DIR"

    # 确定实际输出路径（用于权限设置）
    if [ -z "$OUTPUT_DIR" ] || [ "$OUTPUT_DIR" = "results" ]; then
        OUTPUT_PATH="$WORK_DIR/$OUTPUT_DIR"
    else
        OUTPUT_PATH="$OUTPUT_DIR"
    fi

    echo "============================================================"
    echo "        统一富集分析流程 (GO/KEGG/GSEA)"
    echo "============================================================"
    echo ""
    echo "工作目录: $(pwd)"
    echo "物种: $SPECIES"
    echo ""

    # 构建Rscript参数
    R_ARGS=(
        --args
        --module-root "$SCRIPT_DIR"
        --work-dir "$WORK_DIR"
        --species "$SPECIES"
        --gene-list "$GENE_LIST"
        --expr-matrix "$EXPR_MATRIX"
        --target-genes "$TARGET_GENES"
        --ranked-list "$RANKED_LIST"
        --gmt "$GMT_FILE"
        --analysis-mode "$ANALYSIS_MODE"
        --go-kegg "$GO_KEGG"
        --gsea "$GSEA"
        --gsea-method "$GSEA_METHOD"
        --go-ont "$GO_ONT"
        --go-pvalue "$GO_PVALUE"
        --go-qvalue "$GO_QVALUE"
        --go-padj "$GO_PADJ"
        --gsea-correlation "$GSEA_CORRELATION"
        --gsea-rank-metric "$GSEA_RANK_METRIC"
        --gsea-min-size "$GSEA_MIN_SIZE"
        --gsea-max-size "$GSEA_MAX_SIZE"
        --gsea-pvalue "$GSEA_PVALUE"
        --gsea-padj "$GSEA_PADJ"
        --gsea-exponent "$GSEA_EXPONENT"
        --gsea-eps "$GSEA_EPS"
        --gsea-top-n "$GSEA_TOP_N"
        --qc-variance "$QC_VARIANCE"
        --qc-na-ratio "$QC_NA_RATIO"
        --qc-gmt-match "$QC_GMT_MATCH"
        --plot-width "$PLOT_WIDTH"
        --plot-height "$PLOT_HEIGHT"
        --plot-top-n "$PLOT_TOP_N"
        --plot-color "$PLOT_COLOR"
        --output-dir "$OUTPUT_DIR"
        --log-dir "$LOG_DIR"
        --seed "$SEED"
    )

    # 执行R脚本
    Rscript "${SCRIPT_DIR}/scripts/main.R" "${R_ARGS[@]}"

    # 执行完后，设置输出目录权限
    if [ -n "$OUTPUT_PATH" ] && [ -d "$OUTPUT_PATH" ]; then
        # 先写 .owner（在 chmod 之前，保持默认权限 644）
        CALLER=${SUDO_USER:-$USER}
        echo "$CALLER" > "$OUTPUT_PATH/.owner"
        # 再 chmod 输出目录（1777 = 777 + sticky bit）
        # sticky bit 防止非所有者删除文件
        chmod -R 777 "$OUTPUT_PATH"
        chmod +t "$OUTPUT_PATH"
        # .owner 文件设为只读（只有所有者 iyunlyl 可修改）
        chmod 444 "$OUTPUT_PATH/.owner"
        echo "输出目录权限已设置为 1777: $OUTPUT_PATH"
    fi
}

main "$@"

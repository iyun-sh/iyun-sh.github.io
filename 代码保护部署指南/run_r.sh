#!/bin/bash
# ===========================================================
# run_r.sh — R 脚本自动化运行入口（代码保护版）
#
# 【功能】
#   让任意用户免密以受保护用户身份运行 R 脚本,
#   普通用户无法直接读取 scripts/ 下的 R 源码。
#
# 【权限要求】
#   本文件:          755  (所有人可执行)
#   scripts/ 目录:   700  (仅受保护用户可进入)
#   scripts/*.R:     700  (仅受保护用户可读可执行)
#
# 【sudoers 配置】
#   sudo visudo -f /etc/sudoers.d/secure_scripts
#   ALL ALL=(iyunlyl) NOPASSWD: /media/share/r_automation/run_r.sh
#
# 【用法】
#   执行分析:  ./run_r.sh -o /输出路径 [其他参数...]
#   删除输出:  ./run_r.sh --rm -o /输出路径
#
# 【原理】
#   1. 普通用户执行本脚本
#   2. 脚本检测 whoami ≠ 受保护用户 → 自动 sudo 重新运行自身
#   3. 第二次运行时有权读取 scripts/*.R
#   4. 以受保护用户身份运行 Rscript scripts/analysis.R
#   5. 运行结束后, 输出目录设为 1777 + .owner
#
# 【位置】 /media/share/r_automation/run_r.sh
# ===========================================================

# ---------- 配置区（只需修改这里）----------
TARGET_USER="iyunlyl"                               # 受保护脚本的所有者
SCRIPT_DIR="/media/share/r_automation"               # 本脚本和 scripts/ 所在目录
R_SCRIPT="$SCRIPT_DIR/scripts/analysis.R"            # 被保护的 R 源码
ALLOWED_PREFIX=$(realpath "/media/share/Project" 2>/dev/null || echo "/media/share/Project")  # --rm 和 chmod 白名单路径前缀
# --------------------------------------------

# 保存原始参数
ORIGINAL_ARGS=("$@")

# ---------- 解析参数 ----------
OUTPUT_DIR=""
DO_RM=false
while [[ $# -gt 0 ]]; do
    case $1 in
        -o|--output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --rm)
            DO_RM=true
            shift
            ;;
        *)
            shift
            ;;
    esac
done

# ============================================================
# 身份切换
# ============================================================
if [ "$(whoami)" != "$TARGET_USER" ]; then
    exec sudo -u "$TARGET_USER" "$SCRIPT_DIR/run_r.sh" "${ORIGINAL_ARGS[@]}"
fi

# ============================================================
# 以下所有代码均以受保护用户身份运行
# ============================================================

# 获取实际调用者（SUDO_USER 由 sudo 自动设置，直接运行时 fallback 到 whoami）
CALLER=${SUDO_USER:-$(whoami)}

# ---------- 安全防护: 清理 R 相关危险环境变量 ----------
# R_LIBS / R_LIBS_USER: 用户可注入恶意 R 包, 在 library() 时加载
# R_PROFILE / R_PROFILE_USER: R 启动时自动执行的脚本, 可注入任意代码
# R_ENVIRON / R_ENVIRON_USER: R 启动时读取的环境配置, 可改变库搜索路径
# R_HISTFILE: R 交互历史文件, 可能泄露执行记录
# LD_PRELOAD / LD_LIBRARY_PATH: 可劫持动态库 (R 底层依赖 C 库)
unset R_LIBS R_LIBS_USER R_LIBS_SITE
unset R_PROFILE R_PROFILE_USER
unset R_ENVIRON R_ENVIRON_USER
unset R_HISTFILE R_HISTSIZE
unset LD_PRELOAD LD_LIBRARY_PATH

# ---------- 安全防护: 禁止核心转储 ----------
# R 进程崩溃时 core dump 可能包含源码内容
ulimit -c 0

# ---------- 安全防护: 安全临时目录 ----------
# R 的 tempdir() 默认用 /tmp/, 其他用户可能读取临时文件
export TMPDIR="$SCRIPT_DIR/.tmp"
mkdir -p "$TMPDIR" 2>/dev/null

# ---------- 删除模式 ----------
if [ "$DO_RM" = true ]; then
    if [ -z "$OUTPUT_DIR" ] || [ ! -d "$OUTPUT_DIR" ]; then
        echo "错误：请指定要删除的目录 -o <目录路径>" >&2
        exit 1
    fi

    # 安全校验 1: 白名单路径 + realpath 防穿越
    REAL_OUTPUT=$(realpath "$OUTPUT_DIR" 2>/dev/null)
    if [[ "$REAL_OUTPUT" != "$ALLOWED_PREFIX/"* ]]; then
        echo "错误：只能删除 $ALLOWED_PREFIX/ 下的目录, 拒绝: $REAL_OUTPUT" >&2
        exit 1
    fi

    # 安全校验 1.5: 禁止删除符号链接
    if [ -L "$OUTPUT_DIR" ]; then
        echo "错误：输出目录不能是符号链接, 拒绝删除" >&2
        exit 1
    fi

    # 安全校验 2: .owner 必须存在
    if [ ! -f "$OUTPUT_DIR/.owner" ]; then
        echo "错误：该目录不是由本脚本创建的 (缺少 .owner 标记), 拒绝删除" >&2
        exit 1
    fi

    # 安全校验 3: 身份一致性
    OWNER=$(cat "$OUTPUT_DIR/.owner")
    if [ "$OWNER" != "$CALLER" ]; then
        echo "错误：只能删除自己创建的目录（创建者: $OWNER，当前: $CALLER）" >&2
        exit 1
    fi

    rm -rf "$OUTPUT_DIR"
    echo "已删除: $OUTPUT_DIR"
    exit 0
fi

# ---------- 正常执行模式 ----------
# --vanilla: 禁止加载 .Rprofile / .RData / .Renviron, 确保干净启动
# --no-save: 退出时不保存工作空间 (防止 .RData 泄露到输出目录)
# --no-restore: 不恢复之前的工作空间
Rscript --vanilla --no-save --no-restore "$R_SCRIPT" "${ORIGINAL_ARGS[@]}"
EXIT_CODE=$?

# ---------- 设置输出目录权限 ----------
if [ -n "$OUTPUT_DIR" ] && [ -d "$OUTPUT_DIR" ]; then
    # 防止符号链接攻击
    if [ -L "$OUTPUT_DIR" ]; then
        echo "错误：输出目录不能是符号链接, 拒绝修改权限" >&2
        exit 1
    fi
    # 安全检查：路径白名单 + 项目目录保护
    REAL_OUTPUT=$(realpath "$OUTPUT_DIR" 2>/dev/null || echo "$OUTPUT_DIR")
    REAL_SCRIPT=$(realpath "$SCRIPT_DIR" 2>/dev/null || echo "$SCRIPT_DIR")

    if [ "$CALLER" = "$TARGET_USER" ]; then
        # 本人直接运行，不需要开放权限
        echo "[INFO] 本人($CALLER)运行，输出目录保持默认权限"
    elif [[ "$REAL_OUTPUT" != "$ALLOWED_PREFIX"* ]]; then
        # 输出目录不在允许路径下，拒绝 chmod
        echo "[警告] 输出目录不在允许路径下($ALLOWED_PREFIX/)，跳过 chmod" >&2
        echo "[警告] 请使用 $ALLOWED_PREFIX/ 下的目录作为输出路径" >&2
    elif [[ "$REAL_OUTPUT" == "$REAL_SCRIPT"* ]]; then
        # 输出目录在项目目录内，拒绝 chmod
        echo "[警告] 输出目录在项目目录内($REAL_SCRIPT)，跳过 chmod" >&2
    else
        if [ ! -f "$OUTPUT_DIR/.owner" ]; then
            echo "$CALLER" > "$OUTPUT_DIR/.owner" 2>/dev/null || true
        fi
        chmod -R 777 "$OUTPUT_DIR" 2>/dev/null || true
        chmod +t "$OUTPUT_DIR" 2>/dev/null || true
        # .owner 必须在 chmod -R 777 之后重新设为只读
        chmod 444 "$OUTPUT_DIR/.owner" 2>/dev/null || true
        echo "[INFO] 输出目录权限已设置为 1777（调用者: $CALLER）: $OUTPUT_DIR"
    fi
fi

exit $EXIT_CODE

#!/bin/bash
# ===========================================================
# run_md.sh — 分子对接自动化运行入口
#
# 【功能】
#   让任意用户免密以 iyunlyl 身份运行受保护的 Python 脚本,
#   普通用户无法直接读取 scripts/ 下的源码。
#
# 【权限要求】
#   本文件:          755  (所有人可执行)
#   scripts/ 目录:   700  (仅 iyunlyl 可进入)
#   scripts/*.py:    700  (仅 iyunlyl 可读可执行)
#
# 【sudoers 配置】
#   sudo visudo -f /etc/sudoers.d/secure_scripts
#   ALL ALL=(iyunlyl) NOPASSWD: /media/share/md_automation/run_md.sh
#
# 【用法】
#   执行分析:  ./run_md.sh -o /输出路径 [其他参数...]
#   删除输出:  ./run_md.sh --rm -o /输出路径
#
# 【原理】
#   1. 普通用户执行本脚本
#   2. 脚本检测 whoami ≠ iyunlyl → 自动 sudo -u iyunlyl 重新运行自身
#   3. 第二次运行时 whoami = iyunlyl → 有权读取 scripts/*.py
#   4. 以 iyunlyl 身份运行 python3 scripts/md_automation.py
#   5. 运行结束后, 输出目录设为 1777 (所有人可读写, sticky bit 防互删)
#      并写入 .owner 文件记录原始调用者
#
# 【位置】 /media/share/md_automation/run_md.sh
# ===========================================================

# ---------- 配置区 ----------
TARGET_USER="iyunlyl"                              # 受保护脚本的所有者
SCRIPT_DIR="/media/share/md_automation"             # 本脚本和 scripts/ 所在目录
PYTHON_SCRIPT="$SCRIPT_DIR/scripts/md_automation.py" # 被保护的 Python 源码
# ----------------------------

# 保存原始参数 (sudo 重入后还能拿到完整参数)
ORIGINAL_ARGS=("$@")

# ---------- 解析参数 ----------
OUTPUT_DIR=""      # -o / --output 指定的输出目录
DO_RM=false        # --rm 表示删除模式
while [[ $# -gt 0 ]]; do
    case $1 in
        -o|--output)
            OUTPUT_DIR="$2"   # 取下一个参数作为输出路径
            shift 2
            ;;
        --rm)
            DO_RM=true        # 标记为删除模式
            shift
            ;;
        *)
            shift             # 其余参数透传给 Python 脚本
            ;;
    esac
done

# ============================================================
# 身份切换:
#   当前用户不是 iyunlyl 时, 通过 sudoers 白名单免密切换。
#   exec 会替换当前进程, 避免产生多余的 shell 层级。
#   切换后 SUDO_USER 环境变量保存了原始调用者的用户名。
# ============================================================
if [ "$(whoami)" != "$TARGET_USER" ]; then
    exec sudo -u "$TARGET_USER" "$SCRIPT_DIR/run_md.sh" "${ORIGINAL_ARGS[@]}"
fi

# ============================================================
# 以下所有代码均以 iyunlyl 身份运行
# ============================================================

# ---------- 安全防护: 清理危险环境变量 ----------
# 防止用户通过 PYTHONPATH 注入恶意模块, 或通过 LD_PRELOAD 劫持动态库,
# 从而读取 iyunlyl 的受保护源码并转储到公开路径。
unset PYTHONPATH PYTHONSTARTUP PYTHONHOME PYTHONINSPECT
unset LD_PRELOAD LD_LIBRARY_PATH

# ---------- 安全防护: 禁止 .pyc / __pycache__ 生成 ----------
# Python 默认在 scripts/ 下生成 __pycache__/*.pyc 字节码文件,
# 字节码可被反编译还原源码。设置此变量后不再生成。
export PYTHONDONTWRITEBYTECODE=1

# ---------- 安全防护: 安全临时目录 ----------
# Python tempfile 模块默认写入 /tmp/, 其他用户可能读取临时文件。
# 将临时目录指向受保护路径内。
export TMPDIR="$SCRIPT_DIR/.tmp"
mkdir -p "$TMPDIR" 2>/dev/null

# ---------- 安全防护: 禁止核心转储 ----------
# Python 进程崩溃时可能生成 core dump, 其中包含已加载的源码文本。
# 将 core dump 大小限制为 0, 防止泄露。
ulimit -c 0

# ---------- 删除模式 ----------
# 安全策略:
#   1. 输出目录必须位于白名单路径下 (防止路径逃逸删除 iyunlyl 的其他内容)
#   2. .owner 文件必须存在 (只有本脚本创建的输出目录才有 .owner)
#   3. .owner 中记录的创建者必须与当前调用者一致 (防止用户 A 删除用户 B 的结果)
if [ "$DO_RM" = true ]; then
    # 校验: 必须指定一个已存在的目录
    if [ -z "$OUTPUT_DIR" ] || [ ! -d "$OUTPUT_DIR" ]; then
        echo "错误：请指定要删除的目录 -o <目录路径>" >&2
        exit 1
    fi

    # 安全校验 1: 输出目录必须在白名单路径下
    # 将 OUTPUT_DIR 解析为绝对路径, 防止 ../../../ 路径穿越
    # 注意: 对 symlink 也用 realpath 解析真实路径
    REAL_OUTPUT=$(realpath "$OUTPUT_DIR" 2>/dev/null)
    ALLOWED_PREFIX="/media/share/output"            # ← 只允许删除此路径下的目录
    if [[ "$REAL_OUTPUT" != "$ALLOWED_PREFIX/"* ]]; then
        echo "错误：只能删除 $ALLOWED_PREFIX/ 下的目录, 拒绝: $REAL_OUTPUT" >&2
        exit 1
    fi

    # 安全校验 1.5: 禁止删除符号链接目标
    if [ -L "$OUTPUT_DIR" ]; then
        echo "错误：输出目录不能是符号链接, 拒绝删除" >&2
        exit 1
    fi

    # 安全校验 2: .owner 文件必须存在 (只有本脚本正常创建的输出才有此文件)
    if [ ! -f "$OUTPUT_DIR/.owner" ]; then
        echo "错误：该目录不是由本脚本创建的 (缺少 .owner 标记), 拒绝删除" >&2
        exit 1
    fi

    # 安全校验 3: .owner 中的创建者必须与当前调用者一致
    OWNER=$(cat "$OUTPUT_DIR/.owner")              # 读取创建者
    CALLER=${SUDO_USER:-$USER}                     # 取原始调用者 (sudo 前的用户名)
    if [ "$OWNER" != "$CALLER" ]; then
        echo "错误：只能删除自己创建的目录（创建者: $OWNER，当前: $CALLER）" >&2
        exit 1
    fi

    rm -rf "$OUTPUT_DIR"
    echo "已删除: $OUTPUT_DIR"
    exit 0
fi

# ---------- 正常执行模式 ----------
# 因为此时 whoami = iyunlyl, 而 scripts/ 和 *.py 权限分别为 700,
# 所以可以正常读取并运行 Python 脚本。普通用户直接访问则会被拒绝。
python3 -I "$PYTHON_SCRIPT" "${ORIGINAL_ARGS[@]}"
EXIT_CODE=$?

# ---------- 设置输出目录权限 ----------
# 目的: 让所有用户都能读写输出结果, 同时防止互相删除。
#   chmod -R 777  : 所有人可读写
#   chmod +t      : sticky bit, 只有文件/目录所有者才能删除
#   .owner        : 记录原始调用者, 供 --rm 校验身份
if [ -n "$OUTPUT_DIR" ] && [ -d "$OUTPUT_DIR" ]; then
    # 安全校验: 防止符号链接攻击
    # 如果用户预先在 OUTPUT_DIR 创建指向 iyunlyl 关键目录的 symlink,
    # chmod -R 777 会跟随 symlink 把 iyunlyl 的文件权限全部打开。
    if [ -L "$OUTPUT_DIR" ]; then
        echo "错误：输出目录不能是符号链接, 拒绝修改权限" >&2
        exit 1
    fi
    CALLER=${SUDO_USER:-$USER}                     # 原始调用者
    echo "$CALLER" > "$OUTPUT_DIR/.owner"          # 写入创建者
    chmod -R 777 "$OUTPUT_DIR"                     # 所有人可读写
    chmod +t "$OUTPUT_DIR"                         # sticky bit 防互删
    chmod 444 "$OUTPUT_DIR/.owner"                 # .owner 只读, 防篡改
    echo "输出目录权限已设置: $OUTPUT_DIR (创建者: $CALLER)"
fi

exit $EXIT_CODE

# 代码保护 Wrapper 部署指南

## 概述

通过 `run_md.sh`（Python）和 `run_r.sh`（R）wrapper 脚本，让普通用户**免密执行**受保护的代码，但**无法读取源码**。

核心原理：脚本自动 `sudo -u 受保护用户` 运行，受保护用户拥有源码所有权和独占读取权限。

---

## 目录结构

以 Python 项目为例（R 项目结构完全相同，把 `.py` 换成 `.R`）：

```
/media/share/my_project/          # SCRIPT_DIR
├── run_md.sh                     # wrapper（755，root:root 或任意）
├── .tmp/                         # 安全临时目录（自动创建）
└── scripts/                      # 源码目录（700，受保护用户）
    └── my_script.py              # 源码（700，受保护用户）
```

---

## 部署步骤

### 第 1 步：创建受保护用户（如已有则跳过）

```bash
sudo useradd -r -s /usr/sbin/nologin iyunlyl
```

> `-r` 系统用户，`-s /usr/sbin/nologin` 禁止直接登录。

### 第 2 步：放置源码并设置权限

```bash
# 创建项目目录
sudo mkdir -p /media/share/my_project/scripts

# 放入源码
sudo cp your_script.py /media/share/my_project/scripts/

# 设置所有权
sudo chown -R iyunlyl:iyunlyl /media/share/my_project/scripts

# 核心权限：只有 iyunlyl 能访问
sudo chmod 700 /media/share/my_project/scripts
sudo chmod 700 /media/share/my_project/scripts/*.py   # 或 *.R
```

### 第 3 步：放置 wrapper 脚本

```bash
# 复制对应的 wrapper
sudo cp run_md.sh /media/share/my_project/run_md.sh   # Python 项目
# 或
sudo cp run_r.sh /media/share/my_project/run_r.sh     # R 项目

# 设为所有人可执行
sudo chmod 755 /media/share/my_project/run_md.sh
```

### 第 4 步：修改 wrapper 配置区

打开 wrapper 脚本，修改顶部的配置区（4 个变量）：

```bash
# ---------- 配置区（只需修改这里）----------
TARGET_USER="iyunlyl"                                # ← 受保护用户名
SCRIPT_DIR="/media/share/my_project"                 # ← 项目目录路径
PYTHON_SCRIPT="$SCRIPT_DIR/scripts/my_script.py"     # ← Python 脚本路径
# 或 R 版：
# R_SCRIPT="$SCRIPT_DIR/scripts/analysis.R"          # ← R 脚本路径
ALLOWED_PREFIX="/media/share/output"                 # ← --rm 允许删除的路径前缀
```

> **注意**：wrapper 内的 `exec sudo ... "$SCRIPT_DIR/run_md.sh"` 这行中的脚本名也要和实际文件名一致。

### 第 5 步：配置 sudoers

```bash
sudo visudo -f /etc/sudoers.d/secure_scripts
```

添加一行：

```
ALL ALL=(iyunlyl) NOPASSWD: /media/share/my_project/run_md.sh
```

> 每个项目的 wrapper 都需要单独一行。多个项目示例：
> ```
> ALL ALL=(iyunlyl) NOPASSWD: /media/share/project_a/run_md.sh
> ALL ALL=(iyunlyl) NOPASSWD: /media/share/project_b/run_r.sh
> ALL ALL=(iyunlyl) NOPASSWD: /media/share/project_c/run_md.sh
> ```

### 第 6 步：创建输出白名单目录

```bash
sudo mkdir -p /media/share/output
sudo chmod 1777 /media/share/output   # 所有人可写，sticky bit 防互删
```

### 第 7 步：验证

```bash
# 以普通用户身份测试执行
./run_md.sh -o /media/share/output/test_run

# 确认输出目录已创建且有 .owner
ls -la /media/share/output/test_run/.owner

# 确认源码不可直接读取
cat /media/share/my_project/scripts/my_script.py
# 应输出: Permission denied

# 测试删除
./run_md.sh --rm -o /media/share/output/test_run
```

---

## Python 版 vs R 版 差异对照

| 项目 | `run_md.sh`（Python） | `run_r.sh`（R） |
|------|----------------------|----------------|
| 清理的环境变量 | `PYTHONPATH` `PYTHONSTARTUP` `PYTHONHOME` `PYTHONINSPECT` | `R_LIBS` `R_LIBS_USER` `R_LIBS_SITE` `R_PROFILE` `R_PROFILE_USER` `R_ENVIRON` `R_ENVIRON_USER` `R_HISTFILE` |
| 隔离运行 | `python3 -I script.py` | `Rscript --vanilla --no-save --no-restore script.R` |
| 字节码防护 | `PYTHONDONTWRITEBYTECODE=1` | 不需要（R 无 .pyc） |
| 动态库防护 | `unset LD_PRELOAD LD_LIBRARY_PATH` | 同左 |
| 核心转储防护 | `ulimit -c 0` | 同左 |
| 安全临时目录 | `TMPDIR=$SCRIPT_DIR/.tmp` | 同左 |

---

## 安全措施清单

wrapper 已内置（无需额外配置）：

| 防护项 | 说明 |
|--------|------|
| 环境变量清理 | 防止注入恶意模块/库 |
| 隔离运行模式 | 防止 CWD 导入劫持 |
| 核心转储禁止 | 防止崩溃时泄露源码 |
| 安全临时目录 | 防止 /tmp 文件泄露 |
| 路径白名单 | --rm 只允许删除指定前缀下的目录 |
| realpath 解析 | 防止 `../../../` 路径穿越 |
| symlink 拒绝 | 删除和 chmod 前均检查符号链接 |
| .owner 身份校验 | 只能删除自己创建的输出 |
| 字节码防护 | Python 专属，防 .pyc 反编译 |

建议额外做的系统级加固（需 root）：

```bash
# 防止 strace 附加读取源码
echo 1 > /proc/sys/kernel/yama/ptrace_scope

# 防止 /proc/PID 信息泄露
# 编辑 /etc/fstab，在 proc 行添加 hidepid=2：
# proc /proc proc defaults,hidepid=2 0 0
```

---

## 常见问题

**Q: 新增一个项目需要改几个地方？**
1. 复制 wrapper → 改配置区 4 个变量
2. 源码放入 `scripts/` + `chmod 700`
3. sudoers 加一行

**Q: 一个项目有多个脚本怎么办？**
改 `PYTHON_SCRIPT` / `R_SCRIPT` 变量指向主入口脚本。主入口脚本内部 `import` / `source()` 同目录下的其他文件即可（因为运行时身份已是受保护用户）。

**Q: 能同时保护 Python 和 R 脚本吗？**
可以。同一个项目目录放两个 wrapper（`run_md.sh` + `run_r.sh`），sudoers 各加一行。

**Q: 用户看到 "Permission denied" 但执行也报错？**
检查：wrapper 本身是否 755、sudoers 路径是否完全匹配、受保护用户是否拥有 scripts/ 目录。

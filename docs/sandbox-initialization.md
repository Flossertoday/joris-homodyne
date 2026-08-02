# 沙盒初始化方法与反馈

本页记录托管 Windows 沙盒中的实测摩擦、当前项目入口和仍需由外层权限模型解决的边界。

## 当前方法

从项目根目录运行 `scripts/init_sandbox.ps1`：

```powershell
.\scripts\init_sandbox.ps1
```

脚本可重复运行，并只改变当前 PowerShell 进程和被 Git 忽略的工作树 `.venv`：

- 用进程级 `GIT_CONFIG_COUNT/GIT_CONFIG_KEY_n/GIT_CONFIG_VALUE_n` 注入当前工作树的精确 `safe.directory`，不修改用户全局 Git 配置；
- 优先复用 `.venv`，否则从可用的 `sci` Python 建立 `--system-site-packages` 环境；
- 根据 `requirements-analysis.txt` 探测依赖，只在非 `-CheckOnly` 模式安装缺失项，并禁用无关的 pip 版本检查；
- 优先使用工作树 `DATA/`，否则从 Git 公共目录反查主项目的只读 `DATA/`；
- 设置 `PROJECT_PYTHON`、`JORIS_DATA_ROOT` 和 `JORIS_SAFE_DIRECTORY`，并报告分支、Git 公共目录、依赖与数据状态。

诊断模式：

```powershell
.\scripts\init_sandbox.ps1 -CheckOnly
.\scripts\init_sandbox.ps1 -CheckOnly -AsJson
```

`-CheckOnly` 不创建环境、不安装依赖。`-AsJson` 便于沙盒启动器保存和解析预检证据。

## 2026-08-02 实测与处理

### Git 信任已自动化，写权限仍是独立边界

工作树位于 `C:\Users\Soyo\.codex\worktrees\...`，所有者是宿主账户；命令实际由受限账户运行。普通 `git status` 曾因 `dubious ownership` 退出。原来的单命令回退路径是：

```powershell
$worktree = (Get-Location).Path.Replace('\', '/')
git -c "safe.directory=$worktree" status --short --branch
```

初始化脚本现在把相同的精确工作树信任注入当前进程，使后续普通 Git 子进程直接复用，同时不扩大为全局或通配符信任。

这仍只解决 Git 的所有权检查。该工作树的共享 Git 元数据位于主仓库 `E:\Documents\Codex\joris-homodyne\.git\`，不在当前沙盒写入根中；切换任务分支时曾因无法创建 `.git/worktrees/.../index.lock` 而失败。提交还会写入公共 object store 和分支引用。脚本通过 `common_directory_external` 和 `metadata_write_may_require_approval` 报告这条边界，但不会自行扩大权限。

### Python 入口和依赖已收敛

沙盒启动后的 `PATH` 中不能直接调用 `python` 或 `conda`。Codex 自带 Python 没有 `h5py`；机器上的 `sci` 环境已有 NumPy、SciPy 和 Matplotlib，但当时同样缺少 `h5py`。

手工回退曾是在非 OneDrive 工作树中创建被 `.gitignore` 排除的 `.venv`，复用 `sci` 系统包并单独安装 `h5py`。当前脚本自动执行相同的最小策略，依赖清单由 `requirements-analysis.txt` 维护，并传入 `--disable-pip-version-check`，避免与任务无关的 pip 联网检查。

初始化后使用明确入口运行测试：

```powershell
& $env:PROJECT_PYTHON -m unittest discover -s tests -p "test_*.py" -v
```

### 外部数据根已自动发现

`DATA/` 被 Git 忽略，任务工作树中缺少它是预期行为；真实不可变数据仍在主项目树。初始化脚本先检查工作树 `DATA/`，再通过 Git 公共目录反查主项目 `DATA/` 并设置 `JORIS_DATA_ROOT`。若仍未发现，摘要会返回 `data_available=false`，而不是把“未挂载数据”误判为文件验证失败。

## 实现状态与剩余边界

- 已实现：进程级精确 `safe.directory`、统一 Python 入口、依赖探测与最小安装、机器可读依赖清单、数据根发现、文本/JSON 预检摘要。
- 仍由外层沙盒负责：进入受限环境前附着任务分支，或给精确的 Git 公共元数据目录提供受控写权限。
- 初始化脚本不会修改全局 Git 信任，不会写原始数据，也不会自行扩大文件系统权限。

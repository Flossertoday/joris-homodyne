# 沙盒自动初始化方法与反馈

本项目的环境管理由 Codex 自动完成，不要求用户运行初始化命令。

## 自动入口

`.codex/environments/environment.toml` 是 Codex Desktop 的项目级 Local environment 配置。App 激活该环境后，新聊天创建托管工作树时，Windows setup 自动运行：

```text
powershell -NoProfile -ExecutionPolicy Bypass -File .\scripts\init_sandbox.ps1
```

`scripts/init_sandbox.ps1` 会：

- 优先复用 `.venv`，否则从可用的 `sci` Python 建立 `--system-site-packages` 环境；
- 根据 `requirements-analysis.txt` 探测依赖，只安装缺失项，并禁用无关的 pip 版本检查；
- 优先使用工作树 `DATA/`，否则从 Git 公共目录反查主项目的只读 `DATA/`；
- 检查当前分支、Git 公共目录和精确的 `safe.directory`；
- 把 Python、数据根、Git 边界和依赖状态原子写入忽略 Git 的 `tmp/sandbox-initialization.json`。

自动 setup 在独立初始化进程中运行，因此 `.venv` 和状态文件是跨进程契约；临时环境变量只作为同一 PowerShell 进程内的便利，不作为后续代理的唯一输入。

## Codex 代理的责任

`AGENTS.md` 要求 Codex 在实质工作前读取 `tmp/sandbox-initialization.json`：

- 状态存在且 `ready=true`：直接使用其中的 `python`、`data_root` 和 `git.safe_directory`；
- 状态缺失或失败：代理自行运行 `scripts/init_sandbox.ps1` 恢复并记录原因，不要求用户管理环境；
- Git 所有权检查出现时：代理自动使用状态中的精确路径执行命令级 `git -c safe.directory=...`，不写全局通配符信任。

脚本仍保留 `-CheckOnly` 和 `-AsJson`，供 Codex 诊断或测试使用，不是用户的日常操作入口。

## 2026-08-02 实测问题与处理

### 仓库配置与 App 选择是两个状态

实测已有 Joris 工作树的 `codex.localEnvironmentConfigPath` 是 `__none__`；仅把 `environment.toml` 提交到仓库，不能证明 App 已为现有项目选择它。当前可用的受支持代理工具没有修改该选择的接口，也不应直接改写 Codex 的隐藏全局状态。

因此仓库内采用双层保证：环境被 App 激活时尽早在工作树 setup 阶段运行；否则 `AGENTS.md` 让代理在实质工作前自行运行同一脚本。后者保证用户无需管理环境，同时把 App 选择问题保留为可观测的外层状态。

### Git 信任与写权限是两个边界

工作树所有者是宿主账户，命令由受限账户运行，普通 Git 曾因 `dubious ownership` 退出。精确 `safe.directory` 可以解决所有权检查，但不能授予共享 Git 元数据写权限。

共享 Git 元数据位于主仓库 `E:\Documents\Codex\joris-homodyne\.git\`，可能不在当前沙盒写入根中；切换分支、暂存和提交仍可能需要受控授权。状态用 `common_directory_external` 和 `metadata_write_may_require_approval` 明确报告该边界，初始化脚本不会扩大权限。

### Python 与依赖没有随工作树出现

原始沙盒中不能直接调用 `python` 或 `conda`；Codex 自带 Python 没有 `h5py`，机器 `sci` 环境当时也缺少它。自动 setup 现在创建工作树 `.venv`，复用已有 NumPy、SciPy 和 Matplotlib，并按机器可读清单只补齐缺失包。

### `DATA/` 不进入 Git

任务工作树中缺少 `DATA/` 是预期行为。脚本自动发现主项目数据根并写入状态；找不到时返回 `data_available=false`，不会把“未挂载数据”误判为文件丢失。

## 剩余外层责任

- Codex Desktop 负责保存 Local environment 选择，并在创建新工作树时执行已激活的 setup。
- 外层沙盒负责决定共享 Git 元数据是否可写以及何时请求授权。
- 初始化脚本不修改用户全局 Git 配置、不写原始数据、不自行扩大文件系统权限。

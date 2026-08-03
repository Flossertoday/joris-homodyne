# 沙盒自动初始化方法与反馈

本项目的环境管理由 Codex 自动完成，不要求用户运行初始化命令。

## 自动入口

`.codex/environments/environment.toml` 是 Codex Desktop 的项目级 Local environment 配置。App 激活该环境后，新聊天创建托管工作树时，Windows setup 自动运行：

```text
powershell -NoProfile -ExecutionPolicy Bypass -File .\scripts\init_sandbox.ps1
```

`scripts/init_sandbox.ps1` 会：

- 优先复用隔离的 `.venv`（不带 `--system-site-packages`）；仍共享基础环境 site-packages 的 `.venv` 会被删除重建；
- 每次运行按钉版本的 `requirements.txt` 对工作树 `.venv` 做 `pip install -r` 同步（全部满足时 pip 本地 no-op），并禁用无关的 pip 版本检查；
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

## 2026-08-03 复核：模拟沙盒全流程

### 验证步骤与结果

- 主项目直接运行初始化脚本：READY；`-CheckOnly -AsJson` 诊断输出正确；`tmp/sandbox-initialization.json` 状态与实测一致。
- 测试套件 18 项全部通过。
- 用 `git worktree add --detach tmp/sandbox-test-worktree` 模拟 Codex 托管工作树并在其中运行初始化：自动创建 `--system-site-packages` `.venv`，探测到缺失依赖后 pip 补齐，`DATA/` 正确回退到主项目，`common_directory_external=true` 与 `metadata_write_may_require_approval=true` 正确标记外部 Git 元数据边界，状态文件原子写入工作树 `tmp/`。工作树 Python 可导入 numpy、scipy、h5py、matplotlib。验证后工作树已移除。

### 发现的问题

- `sci` 机器环境仍缺 `h5py`：每个新工作树 `.venv` 都要 pip 安装一次 h5py（本次使用 3.2 MB 本地缓存 wheel）。安装目标是工作树自己的 site-packages，未污染 `sci` 环境；但缓存被清或离线时，新沙盒初始化会在这一步失败。（2026-08-03 已解决，见下节。）
- detached HEAD 的工作树状态中 `branch` 为 `null`（人类可读输出为空行）。真实 Codex 工作树带分支，此处不算缺陷，但状态消费方不能假设 `branch` 非空。

## 2026-08-03 处理：基础 sci 环境补齐 h5py

经用户决定，h5py 直接装入基础环境，而不是仓库内置 wheelhouse：

- 目标是 `%USERPROFILE%\.conda\envs\sci`（Python 3.10.13；`.venv/pyvenv.cfg` 的 `home` 确认它是当前 `.venv` 的基础解释器，也是脚本 `Find-BasePython` 的第一候选）。pip 装入 h5py 3.16.0；已有 numpy 1.26.4 满足约束，未改动其他包。
- 安装后基础 sci 可导入 numpy、scipy、h5py、matplotlib 全部四项。
- 工作树模拟验证（`git worktree add --detach` 后在其中运行初始化）：`ready=true`、`environment_created=true`、`missing=[]`、`data_root` 正确回退主项目 `DATA`，输出无 pip install 记录——工作树 `.venv` 经 `--system-site-packages` 直接命中基础环境的 h5py。新工作树初始化从此零网络依赖（该结论基于当时的 `--system-site-packages` 架构；同日稍后改为隔离 `.venv`，见下节，不再适用）；验证后工作树已移除。
- 主项目回归：`-CheckOnly -AsJson` 为 `ready=true`、`missing=[]`，完整测试套件 20 项全部通过。
- 注意：fallback 候选 `%USERPROFILE%\anaconda3\envs\sci` 仍缺 h5py；只有 `.conda\envs\sci` 不可用时才会被解析到，届时需同样补齐。
- 此前“安装目标是工作树自己的 site-packages，未污染 sci 环境”的原则随之改变：h5py 现在是基础 sci 环境契约的一部分。

## 2026-08-03 补充：requirements.txt 钉版本与隔离 .venv

用户同日把依赖契约改为标准钉版清单，并把工作树环境改为隔离架构：

- `requirements-analysis.txt` 改名为 `requirements.txt`，依赖全部钉版：numpy 1.26.4、scipy 1.11.4、h5py 3.16.0、matplotlib 3.8.0、trackpy 0.7（trackpy 为当日新增）。
- 初始化脚本重新设计：只复用“隔离”的 `.venv`（`pyvenv.cfg` 不含 `include-system-site-packages = true`）；仍共享基础环境 site-packages 的 `.venv` 会被删除重建，新环境不再使用 `--system-site-packages`。
- 每次非 `-CheckOnly` 运行都按钉版清单对工作树 `.venv` 执行 `pip install -r`：全部满足时 pip 本地 no-op，版本不符时由 pip 修正；缺失探测保持按模块名存在性。显式外部 `-PythonPath` 永不安装。`-CheckOnly` 保持只读，对现有 `.venv` 出报告。
- 后果与边界变化：新工作树首次初始化需要 pip 缓存或网络来安装全部钉版包——上节“零网络依赖”的结论不再适用。基础 sci 环境的 h5py 对工作树不再关键，但让“尚无 `.venv` 时的 `-CheckOnly`”仍能报告 ready。
- 主项目 `.venv` 目前仍是旧的共享架构，已补装 trackpy 0.7 使清单闭合；下一次完整初始化会自动把它迁移为隔离 `.venv`。

## 剩余外层责任

- Codex Desktop 负责保存 Local environment 选择，并在创建新工作树时执行已激活的 setup。
- 外层沙盒负责决定共享 Git 元数据是否可写以及何时请求授权。
- 初始化脚本不修改用户全局 Git 配置、不写原始数据、不自行扩大文件系统权限。

# 沙盒初始化反馈

本页记录托管 Windows 沙盒中的实测摩擦和对后续初始化方案的要求，不把当前临时做法冻结为唯一方案。

## 2026-08-02 实测

### Git 信任与写权限是两个独立问题

工作树位于 `C:\Users\Soyo\.codex\worktrees\...`，所有者是宿主账户；命令实际由受限账户运行。普通 `git status` 因 `dubious ownership` 退出。以下命令可在不修改用户全局 Git 配置的前提下恢复单次只读操作；路径必须规范化为正斜杠，本次直接传入反斜杠路径没有匹配成功。

```powershell
$worktree = (Get-Location).Path.Replace('\', '/')
git -c "safe.directory=$worktree" status --short --branch
```

这只解决 Git 的所有权信任检查。该工作树的共享 Git 元数据仍位于主仓库 `E:\Documents\Codex\joris-homodyne\.git\`，不在当前沙盒写入根中；切换任务分支时首次失败为无法创建 `.git/worktrees/.../index.lock`，提升该 Git 操作权限后才成功。提交也会同时写入工作树索引、公共 object store 和分支引用，不能把 `safe.directory` 当作写权限修复。

### Python 入口和依赖没有随项目进入沙盒

沙盒启动后的 `PATH` 中不能直接调用 `python` 或 `conda`。Codex 自带 Python 可以运行，但没有 `h5py`；机器上的两个 `sci` 环境可以通过绝对路径启动，已有 NumPy、SciPy 和 Matplotlib，但同样缺少 `h5py`。因此项目文档中的裸 `python -m unittest ...` 不能作为全新沙盒的第一条命令。

本次在非 OneDrive 的工作树内建立了被 `.gitignore` 排除的临时环境，复用 `sci` 的系统包并只补装缺失的 `h5py`：

```powershell
& "$env:USERPROFILE\.conda\envs\sci\python.exe" -m venv --system-site-packages .venv
.\.venv\Scripts\python.exe -m pip install h5py
.\.venv\Scripts\python.exe -m unittest discover -s tests -p "test_*.py" -v
```

安装成功；`pip` 只额外报告无法检查自身最新版本，不影响包安装或测试。该做法适合作为当前回退路径，但项目还没有机器可读的最小依赖清单，也没有决定是复用 `sci`、克隆环境还是维护专用环境。

### 数据不会自动出现在任务工作树

`DATA/` 被 Git 忽略，任务工作树中缺少它是预期行为；真实不可变数据仍在主项目树。只运行使用临时 HDF5 的单元测试不受影响，但任何数据级验证都需要初始化阶段显式提供主项目 `DATA/` 的只读位置，不能把“工作树内不存在”误判为数据丢失。

## 对后续初始化方案的要求

建议把初始化拆成可独立诊断的预检，而不是用一个模糊的“环境失败”覆盖所有边界：

1. 在进入受限沙盒前把工作树附着到任务分支，或为精确的公共 Git 元数据目录提供受控写权限。
2. 以进程级配置注入当前工作树的 `safe.directory`；不要要求代理永久修改用户全局 Git 配置。
3. 暴露一个明确的项目 Python 入口，并在创建环境前探测 `numpy`、`scipy`、`h5py`、`matplotlib`，只报告或补齐缺失项。
4. 为核心依赖增加机器可读的最小清单；在方案确定前，保留 `sci` 加工作树 `.venv` 的回退路径。
5. 向任务提供真实 `DATA/` 的只读根或明确的外部路径；初始化日志应区分“未挂载数据”和“文件验证失败”。
6. 预检结束时输出工作树、分支、Git 公共目录、实际 Python、缺失依赖和数据根的简短摘要，便于后续任务直接复用证据。

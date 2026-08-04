# Project context

- `README.md` is the entry point; `MEMORY.md` holds current structural context; reviewable detail lives in `docs/`.
- `DATA/` contains local raw experiments and experiment-specific derived copies and is ignored by Git.
- Scientific Python uses an isolated per-worktree `.venv` created from the local `sci` interpreter (no shared site-packages) and synced from the pinned `requirements.txt`; the `sci` environment itself is never installed into.
- When Codex Desktop activates the checked-in Local environment, it automatically runs `scripts/init_sandbox.ps1` through `.codex/environments/environment.toml` while creating a managed worktree.
- Before substantive work, read `tmp/sandbox-initialization.json`. If it is missing or not ready, run the initializer yourself as the automatic recovery path and report the failure.
- Do not ask the user to initialize the sandbox or manage its Python dependencies.
- Use the recorded `git.safe_directory` with command-local `git -c safe.directory=...` when ownership checks require it. The initializer does not grant access to external shared Git metadata.
- Test command after initialization: `.\.venv\Scripts\python.exe -m unittest discover -s tests -p "test_*.py" -v`.
- Matrix-valued figures (camera frames, heatmaps, density maps, frame sums) must include a colorbar making the intensity mapping (LUT) explicit: the colormap, the normalization (linear, sqrt, log, percentile clipping), and the displayed value range.

## Git

- This is a solo project that may use Codex, Claude Code, Kimi Code, Cursor, and other agents. `main` is the stable integration branch.
- Branch prefixes identify their owner: only Codex-managed branches/worktrees use `codex/<topic>`; Kimi Code branches use `kimi/<topic>`; other tools use their own namespace.
- An independent task normally has its own branch and worktree, based on current `main`.
- A useful handoff is a focused commit plus branch name, commit hash, verification evidence, risks, and relevant untracked files.
- Task branches are integrated into `main` with a regular merge that preserves their commits. Pushing and worktree cleanup are separate follow-on actions.
- Kimi Code must commit through `scripts/kimi-commit.sh` (Git Bash), never bare `git commit`. The script sets author and committer to `Kimi K3 <kimi-k3@users.noreply.github.com>` and appends `Model:` / `Thinking:` trailers parsed from `~/.kimi-code/config.toml`; if the session switched model or thinking effort mid-run, pass the real values via the `KIMI_MODEL` / `KIMI_THINKING` environment variables.

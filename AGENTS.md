# Project context

- `README.md` is the entry point; `MEMORY.md` holds current structural context; reviewable detail lives in `docs/`.
- `DATA/` contains local raw experiments and experiment-specific derived copies and is ignored by Git.
- Scientific Python uses a clone of the local `sci` environment outside OneDrive.
- Test command: `python -m unittest discover -s tests -p "test_*.py" -v`.

## Git

- This is a solo project that may use Codex, Claude Code, Kimi Code, Cursor, and other agents. `main` is the stable integration branch.
- Branch prefixes identify their owner: only Codex-managed branches/worktrees use `codex/<topic>`; other tools use their own namespace.
- An independent task normally has its own branch and worktree, based on current `main`.
- A useful handoff is a focused commit plus branch name, commit hash, verification evidence, risks, and relevant untracked files.
- Task branches are integrated into `main` with a regular merge that preserves their commits. Pushing and worktree cleanup are separate follow-on actions.

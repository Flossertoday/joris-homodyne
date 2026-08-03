# Kimi Code operating rules for this project

These rules are specific to Kimi Code. Cross-tool project conventions live in the root `AGENTS.md`, which takes precedence on any conflict about shared repository state.

## Commits

- Commit finished task work automatically to a `kimi/<topic>` feature branch based on current `main`, without asking for per-commit confirmation. Commit through `scripts/kimi-commit.sh` as the root `AGENTS.md` requires.
- Merges into `main` and pushes still require explicit user approval.

## Branch and worktree cleanup

Cleanup is tiered by ownership and mechanical verifiability, not blanket approval:

- **Merged, unpinned branches**: delete automatically once `git branch --merged main` confirms the merge; the merge commit preserves their commits in `main`.
- **Self-created worktrees**: the creator cleans up. Remove a Kimi-created worktree once the task is done, its working tree is clean (no uncommitted or untracked changes), and its branch is merged.
- **Codex-managed and other foreign worktrees** (`~/.codex/worktrees/*`, or anything Kimi Code did not create): never touch; list them in the handoff report instead. A branch pinned by such a worktree can only be deleted after its owner removes the worktree.
- **Unmerged branches**: never delete.

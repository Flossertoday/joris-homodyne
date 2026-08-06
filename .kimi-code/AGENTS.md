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

## Working style (standing user instructions, 2026-08-06)

- **Challenge the user's methods.** When the user's proposed approach conflicts with what the data shows, say so and present the evidence instead of just complying. Always.
- **Explain terminology.** Spell out technical jargon and acronyms on first use (CFD, IRF, ...); never assume shared vocabulary.
- **Self-review figures before showing them.** View every generated image yourself before presenting it to the user, and check that it actually conveys the intended information — right quantity, right axis scale, readable, sensible panel count. Iterate until it does.

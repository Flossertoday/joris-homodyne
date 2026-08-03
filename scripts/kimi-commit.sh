#!/usr/bin/env bash
# kimi-commit.sh — required commit entry point for Kimi Code in this repository.
#
# Stamps author and committer as Kimi K3 and appends trailers recording the
# model and thinking effort actually in use:
#
#   Model: kimi-code/k3
#   Thinking: max
#
# Values are parsed from ~/.kimi-code/config.toml (default_model, and the last
# default_effort inside the active [models."<name>"] section, including its
# .overrides subsection). If the session switched model or effort mid-run,
# pass the real values via the KIMI_MODEL / KIMI_THINKING environment
# variables. All arguments are forwarded to git commit.
set -euo pipefail

cfg="${KIMI_CONFIG:-$HOME/.kimi-code/config.toml}"

model="${KIMI_MODEL:-}"
if [ -z "$model" ] && [ -f "$cfg" ]; then
  model=$(grep -m1 '^default_model' "$cfg" | cut -d'"' -f2 || true)
fi
model="${model:-unknown}"

effort="${KIMI_THINKING:-}"
if [ -z "$effort" ] && [ -f "$cfg" ] && [ "$model" != "unknown" ]; then
  effort=$(awk -v sec="[models.\"$model\"" '
    index($0, sec) == 1 { in_s = 1; next }
    /^\[/ { in_s = 0 }
    in_s && /^default_effort/ {
      v = $0
      sub(/^[^"]*"/, "", v)
      sub(/".*$/, "", v)
    }
    END { if (v != "") print v }
  ' "$cfg")
fi
effort="${effort:-unknown}"

exec git -c user.name="Kimi K3" -c user.email="kimi-k3@users.noreply.github.com" \
  commit --trailer "Model: $model" --trailer "Thinking: $effort" "$@"

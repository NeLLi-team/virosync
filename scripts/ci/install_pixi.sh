#!/usr/bin/env bash
set -euo pipefail

readonly PIXI_VERSION="0.72.0"
readonly PIXI_SHA256="6304fe3178f3036e2c95151bbb318592fae5c31a77f5a6f4319bb023a479d4b9"
readonly PIXI_URL="https://github.com/prefix-dev/pixi/releases/download/v0.72.0/pixi-x86_64-unknown-linux-musl"

if [[ "$(uname -s)" != "Linux" || "$(uname -m)" != "x86_64" ]]; then
    echo "The pinned CI Pixi installer supports Linux x86_64 only." >&2
    exit 2
fi

temporary="$(mktemp)"
trap 'rm -f "$temporary"' EXIT

curl \
    --proto '=https' \
    --tlsv1.2 \
    --retry 3 \
    --retry-delay 2 \
    --connect-timeout 20 \
    --max-time 300 \
    --fail \
    --silent \
    --show-error \
    --location \
    "$PIXI_URL" \
    --output "$temporary"

actual_sha256="$(sha256sum "$temporary" | cut -d ' ' -f 1)"
if [[ "$actual_sha256" != "$PIXI_SHA256" ]]; then
    echo "Pixi ${PIXI_VERSION} SHA-256 mismatch." >&2
    exit 1
fi

pixi_bin_dir="${PIXI_HOME:-$HOME/.pixi}/bin"
install -d -m 0755 "$pixi_bin_dir"
install -m 0755 "$temporary" "$pixi_bin_dir/pixi"
"$pixi_bin_dir/pixi" --version

if [[ -n "${GITHUB_PATH:-}" ]]; then
    printf '%s\n' "$pixi_bin_dir" >>"$GITHUB_PATH"
fi

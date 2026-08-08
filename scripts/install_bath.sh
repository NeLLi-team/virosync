#!/usr/bin/env bash
# Build BATH into .bath/ so the frameshift screen works through pixi.
#
# BATH is absent from the configured conda-forge and bioconda channels, so it
# cannot be a pixi dependency. This builds the revisions pinned in
# docs/FRAMESHIFT_SCREENING.md. pixi.toml puts .bath/bin on PATH for every
# task, so `pixi run example-frameshift` works once this has run.
set -euo pipefail

BATH_REV=7842ebd58b96591b4b60863ee5c33e49eb79eccc
EASEL_REV=0f4e71832d6ba1e4c65039ba4b4663c546a041fa

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PREFIX="$ROOT/.bath"
BUILD="$PREFIX/src"

if [ "${1:-}" != "--force" ] \
   && [ -x "$PREFIX/bin/bathsearch" ] && [ -x "$PREFIX/bin/bathconvert" ]; then
    echo "BATH already installed at $PREFIX/bin (use --force to rebuild)"
    exit 0
fi

rm -rf "$BUILD"
mkdir -p "$BUILD"
git clone --quiet https://github.com/TravisWheelerLab/BATH.git "$BUILD/BATH"
git -C "$BUILD/BATH" checkout --quiet "$BATH_REV"
git clone --quiet https://github.com/EddyRivasLab/easel.git "$BUILD/BATH/easel"
# A plain clone does not contain the pinned commit, so fetch it explicitly.
git -C "$BUILD/BATH/easel" fetch --quiet origin "$EASEL_REV"
git -C "$BUILD/BATH/easel" checkout --quiet "$EASEL_REV"

cd "$BUILD/BATH"
# The repo ships configure.ac only, and autoconf is not a project dependency.
# pixi exec supplies it for this one command without touching the manifest.
pixi exec --spec autoconf --spec m4 --spec perl -- autoconf
./configure --prefix="$PREFIX"
make -j "${BATH_BUILD_JOBS:-8}"
make install

"$PREFIX/bin/bathconvert" -h >/dev/null
"$PREFIX/bin/bathsearch" -h >/dev/null
echo "BATH installed: $("$PREFIX/bin/bathsearch" -h 2>&1 | grep -o 'BATH [0-9.]*' | head -1) -> $PREFIX/bin"

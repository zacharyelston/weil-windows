#!/usr/bin/env bash
# Run a script inside the pinned native environment with this checkout mounted at /work.
# Files written (logs/, data/) belong to the calling user.
set -euo pipefail
here="$(cd "$(dirname "$0")/.." && pwd)"
image="${IMAGE:-weil-windows-env:py3.14.4}"
exec docker run --rm -i -v "$here":/work -w /work -u "$(id -u):$(id -g)" -e HOME=/tmp "$image" "$@"

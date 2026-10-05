#!/usr/bin/env bash
# Run a script inside the pinned rh2 environment with this checkout mounted at /rh2.
# Files written (logs/, data/) belong to the calling user.
set -euo pipefail
here="$(cd "$(dirname "$0")/.." && pwd)"
image="${RH2_IMAGE:-weil-windows-env:py3.14.4}"
exec docker run --rm -i -v "$here":/rh2 -w /rh2 -u "$(id -u):$(id -g)" -e HOME=/tmp "$image" "$@"

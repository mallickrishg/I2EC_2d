#!/usr/bin/env bash
# Build the HTML documentation into docs/_build/html.
#
#   docs/build.sh            build
#   docs/build.sh --open     build, then open index.html
#   docs/build.sh --serve    build, then serve on http://127.0.0.1:8000
#
# Uses $SPHINXBUILD if set, else ~/.venvs/matlab-docs/bin/sphinx-build if it
# exists, else sphinx-build on the PATH.
set -euo pipefail

DOCS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT_DIR="$DOCS_DIR/_build/html"

if [[ -z "${SPHINXBUILD:-}" ]]; then
  if [[ -x "$HOME/.venvs/matlab-docs/bin/sphinx-build" ]]; then
    SPHINXBUILD="$HOME/.venvs/matlab-docs/bin/sphinx-build"
  else
    SPHINXBUILD="sphinx-build"
  fi
fi

"$SPHINXBUILD" -b html "$DOCS_DIR" "$OUT_DIR"

case "${1:-}" in
  --open)
    if command -v open >/dev/null; then open "$OUT_DIR/index.html"; else xdg-open "$OUT_DIR/index.html"; fi
    ;;
  --serve)
    echo "Serving $OUT_DIR at http://127.0.0.1:8000 (Ctrl-C to stop)"
    cd "$OUT_DIR" && python3 -m http.server 8000 --bind 127.0.0.1
    ;;
  "") ;;
  *) echo "unknown option: $1 (use --open or --serve)" >&2; exit 2 ;;
esac

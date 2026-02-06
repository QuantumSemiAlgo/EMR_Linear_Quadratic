#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OUT="$ROOT/output"

EMRFILE="$(find "$OUT" -type f -name 'EMRdata.out' -print0 | xargs -0 ls -t 2>/dev/null | head -n 1 || true)"

if [[ -z "${EMRFILE}" ]]; then
  echo "No EMRdata.out found under $OUT"
  exit 1
fi

echo "Using: $EMRFILE"

docker run --rm \
  --user "$(id -u):$(id -g)" \
  -e MPLCONFIGDIR=/tmp/mpl \
  -v "$OUT:/app/output" \
  emr-fem:hermite \
  bash -lc "cd /app/src && python3 plot_emr_vs_h.py /app/output/${EMRFILE#"$OUT/"}"

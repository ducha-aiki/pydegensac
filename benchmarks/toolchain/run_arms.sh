#!/usr/bin/env bash
# Time every arm interleaved, one pass each per round, so slow drift in the
# machine cannot be mistaken for a difference between arms.
#
#   ./run_arms.sh h 7 > h.jsonl
#   ARMS="pkg-local-head pkg-ml2014-head" ./run_arms.sh f 5 > f.jsonl
#
# Arms are package directories under benchmarks/.ab, as built by
# build_in_image.sh (plus a local `pip install --target`); see README.md.
set -euo pipefail
cd "$(dirname "$0")"

PROBLEM="${1:-h}"
ROUNDS="${2:-5}"
AB="$(cd ../.ab && pwd)"
PY="${PYTHON:-$AB/env_np1/bin/python}"
#: A P-core. This box is a 14700K: cpu0-15 are P-cores, 16-27 are E-cores, and
#: landing on an E-core costs ~25% — enough to invent a toolchain difference.
CPU="${CPU:-2}"

ARMS="${ARMS:-pkg-local-0.2.2 pkg-ml2014 pkg-ml228 pkg-pypi-0.2.2}"

for r in $(seq 1 "$ROUNDS"); do
    for a in $ARMS; do
        [ -d "$AB/$a" ] || { echo "missing $AB/$a" >&2; continue; }
        PYTHONPATH="$AB/$a" taskset -c "$CPU" "$PY" iso_bench.py "$PROBLEM" \
            --reps 1 --label "$a" 2>&1 | awk 'NR==1'
        # awk, not `head -1`: head closes the pipe, the driver dies of SIGPIPE
        # on its summary line, and `set -e` kills the sweep partway through.
    done
done

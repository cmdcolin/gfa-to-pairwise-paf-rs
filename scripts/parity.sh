#!/usr/bin/env bash
# parity.sh <python script> <rust binary> <gfa> [converter args...]
# Runs both converters on one input and fails unless stdout, every
# chrom.sizes file, the exit code and stderr (timings masked) match.
set -euo pipefail

python_script=$1
rust_binary=$2
gfa=$3
shift 3

work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

mask() {
  sed -E 's/[0-9]+s$/Ns/; s/in [0-9.]+s \([0-9]+ MB\/s\)/in Ns (N MB\/s)/' "$1"
}

converter() {
  local name=$1
  shift
  local sizes=()
  if [[ ${WITH_SIZES:-} ]]; then
    sizes=(--chrom-sizes-dir "$work/$name.sizes")
  fi
  set +e
  "$@" "$gfa" "${sizes[@]}" "${args[@]}" >"$work/$name.paf" 2>"$work/$name.err"
  echo $? >"$work/$name.code"
  set -e
}

args=("$@")
converter python python3 "$python_script"
converter rust "$rust_binary"

status=0
cmp "$work/python.code" "$work/rust.code" || { echo "exit codes differ: $(cat "$work/python.code") vs $(cat "$work/rust.code")"; status=1; }
cmp "$work/python.paf" "$work/rust.paf" || status=1
diff <(mask "$work/python.err") <(mask "$work/rust.err") || status=1
if [[ ${WITH_SIZES:-} && -d $work/python.sizes ]]; then
  diff -r "$work/python.sizes" "$work/rust.sizes" || status=1
fi
echo "$(basename "$gfa") ${args[*]}: $(wc -l <"$work/rust.paf") rows, $( ((status)) && echo DIFFERENT || echo identical)"
exit $status

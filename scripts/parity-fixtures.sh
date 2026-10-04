#!/usr/bin/env bash
# parity-fixtures.sh <python script> <rust binary>
# Every fixture under every flag set the tests use, through both converters.
set -uo pipefail

here=$(dirname "$0")
status=0
for gfa in "$here"/../tests/data/*.gfa; do
  while read -r -a args; do
    WITH_SIZES=1 "$here/parity.sh" "$1" "$2" "$gfa" "${args[@]}" || status=1
  done <<'ARGS'
--reference GRCh38#0
--reference GRCh38
--reference GRCh38 --no-x
--reference GRCh38 --max-gap 2
--reference GRCh38 --min-block 15 --queries HG01123#1
--reference GRCh38 --hold-queries
--reference CHM13
--reference GRCh38 --queries HG01109.1,,HG00097,nope#3
--reference HG01109#1
--reference HG01123#1 --no-x --max-gap 0
--ref GRCh38 --max-g=-1
--reference GRCh38 --queries=
ARGS
done
exit $status

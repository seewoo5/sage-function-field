#!/usr/bin/env bash
# Run all .sage tests in ff/test

set -u -o pipefail

SCRIPT_DIR="$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)"
TEST_DIR="${SCRIPT_DIR}/test"

if [[ ! -d "$TEST_DIR" ]]; then
  echo "Test directory not found: $TEST_DIR" >&2
  exit 1
fi

# Build a newline-separated, sorted list of .sage test files
tmpf="$(mktemp)"
trap 'rm -f "$tmpf"' EXIT
find "$TEST_DIR" -type f -name '*.sage' -print0 | tr '\0' '\n' | sort > "$tmpf"

if ! [ -s "$tmpf" ]; then
  echo "No .sage tests found in $TEST_DIR" >&2
  exit 1
fi

fail_count=0
pass_count=0

while IFS= read -r t; do
  rel="${t#$SCRIPT_DIR/}"
  echo "==> sage $rel"
  if sage "$t"; then
    ((pass_count++))
  else
    ((fail_count++))
    echo "FAILED: $rel" >&2
  fi
done < "$tmpf"

echo
echo "Summary: ${pass_count} passed, ${fail_count} failed, $((pass_count + fail_count)) total."
if [[ $fail_count -gt 0 ]]; then
  exit 1
fi
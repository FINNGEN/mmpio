#! /usr/bin/env bash
set -euo pipefail

tdir=$(dirname $(realpath $0))

source ${tdir}/common.sh

for test_file in ${tdir}/test_*/run.sh; do
  log "Running test for $test_file"
  bash $test_file
done


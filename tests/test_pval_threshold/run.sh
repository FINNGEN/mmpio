#! /usr/bin/env bash
set -euxo pipefail

# Setup test
tdir=$(dirname $(realpath $0))
cd $tdir

source ${tdir}/../common.sh

build_mmpio


# Run end-to-end test
cat data_sumstats_2rows.tsv | gzip > data_sumstats_2rows.tsv.gz

./mmpio --config config.json --output data_out.tsv

diff data_expected.tsv data_out.tsv

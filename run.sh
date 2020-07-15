#!/bin/bash
set -e

ROOT_DIR=$(dirname $(realpath $0))
EXECUTABLE="$ROOT_DIR/build/quasistatic_pic"

INPUT_PATH=$(realpath $1)
INPUT_DIR=$(dirname $INPUT_PATH)
INPUT_NAME=$(basename $INPUT_PATH)

RESULTS_DIR="$INPUT_DIR/${INPUT_NAME%.*}"

if [[ -d $RESULTS_DIR ]]; then
    echo "Output directory exists $RESULTS_DIR"
else
    echo "Output directory does not exist, creating $RESULTS_DIR"
    mkdir $RESULTS_DIR
fi

cd $RESULTS_DIR
cp $INPUT_PATH .

export OMP_NUM_THREADS=$2

$EXECUTABLE $INPUT_NAME
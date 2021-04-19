#!/usr/bin/env bash


EDGE_FILE=$1
NUM_TESTS="${2:-3}"
NUM_METHODS="${3:-6}" # Change this if you want to use LMC or not use earlier methods
TEST_ID=$(basename $EDGE_FILE .edge)
DIRECTORY=$(dirname $EDGE_FILE)
CSV_FILE="${DIRECTORY}/${TEST_ID}.csv"
OUT_FILE="${DIRECTORY}/${TEST_ID}.out"
TIME=$(date +%m-%d-%y-%R:%S)

if [[ -f "$CSV_FILE" ]]; then
	BACKUP_FILE="${CSV_FILE}.${TIME}.backup"
	echo "$CSV_FILE already exists, backing it up as $BACKUP_FILE"
	mv $CSV_FILE $BACKUP_FILE
fi

if [[ -f "$OUT_FILE" ]]; then
	BACKUP_FILE="${OUT_FILE}.${TIME}.backup"
	echo "$OUT_FILE already exists, backing it up as $BACKUP_FILE"
	mv $OUT_FILE $BACKUP_FILE
fi

# Show all commands from this point on
set -x

# /tmp and /temp are on local FS -- prevents I/O from being a bottleneck
TMP=$(mktemp)
echo "start time: $TIME" >> $TMP
./build/runExperiment ${EDGE_FILE} ${NUM_TESTS} ${CSV_FILE} ${NUM_METHODS} | tee -a ${TMP}
cp ${TMP} ${OUT_FILE}
rm $TMP
END_TIME=$(date +%m-%d-%y-%R:%S)
echo "end time: $END_TIME" >> $OUT_FILE

#!/bin/bash

# Define base directories and command
BASE_DIR=""
COMMAND="$BASE_DIR/PopLDdecay"
INPUT_DIR="$BASE_DIR/mis"
OUT_DIR="$BASE_DIR/output"

# Create output directory if it doesn't exist
mkdir -p $OUT_DIR

# Loop through chromosomes 1 to 12
for i in $(seq 1 12)
do
    INPUT_FILE="${INPUT_DIR}/eachvariantschr$(printf "%02d" $i)converted.txt"
    OUTPUT_PREFIX="Chr${i}"
    
    # Run PopLDdecay command
    $COMMAND -InGenotype $INPUT_FILE -OutStat $OUTPUT_PREFIX -OutFilterSNP
    
    echo "Processed chromosome $i"
done

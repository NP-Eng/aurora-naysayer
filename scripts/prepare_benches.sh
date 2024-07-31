#!/bin/bash

# Arguments
NUM_SQUARINGS=$1

# Python script and directories
PYTHON_SCRIPT="scripts/gen_solidity_repeated_squaring.py"
CIRCOM_DIRECTORY="circom"

# Print the current working directory
echo "Current working directory: $(pwd)"

# Call the Python script to generate the circom and witness files
echo "Executing Python script: python3 $PYTHON_SCRIPT $NUM_SQUARINGS"
python3 $PYTHON_SCRIPT $NUM_SQUARINGS
if [ $? -ne 0 ]; then
    echo "Failed to execute Python script"
    exit 1
fi
echo "Python script executed successfully"

# Compile the circom file to get the R1CS and WASM files
CIRCOM_FILE="$CIRCOM_DIRECTORY/repeated_squaring_${NUM_SQUARINGS}.circom"
OUTPUT_DIR="test-data/"

echo "Compiling circom file: circom $CIRCOM_FILE --r1cs --wasm -o $OUTPUT_DIR"
circom $CIRCOM_FILE --r1cs --wasm -o $OUTPUT_DIR
if [ $? -ne 0 ]; then
    echo "Circom compilation failed"
    exit 1
fi
echo "Circom compilation succeeded"

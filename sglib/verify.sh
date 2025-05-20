#!/bin/bash
# Script to verify the Robust Savitzky-Golay Filter library
# Author: Mapoet
# Date: 2025-05-20

echo "Verifying Robust Savitzky-Golay Filter library..."

# Check for fortran compiler
if ! command -v gfortran &> /dev/null
then
    echo "ERROR: gfortran not found. Please install gfortran compiler."
    exit 1
fi

# Check for OpenBLAS and LAPACK
if ! ldconfig -p | grep -q libopenblas
then
    echo "WARNING: libopenblas not found. Please install OpenBLAS."
    exit 1
fi

if ! ldconfig -p | grep -q liblapack
then
    echo "WARNING: liblapack not found. Please install LAPACK."
    exit 1
fi

# Build the library
echo "Building the library..."
make clean
make release

if [ $? -ne 0 ]; then
    echo "ERROR: Build failed."
    exit 1
fi

# Run the test program
echo "Running test program..."
make test

if [ $? -ne 0 ]; then
    echo "ERROR: Test execution failed."
    exit 1
fi

# Check if output file exists
if [ ! -f "./build/sg_filter_results.dat" ]; then
    echo "ERROR: Test output file not generated."
    exit 1
fi

echo "Test completed successfully."
echo "Library verification complete. You can use the library in your projects now."
echo "See README.md for usage instructions."
python3 ./scripts/plot_results.py -i ./build/sg_filter_results.dat -o ./sg_filter_results.png
exit 0 
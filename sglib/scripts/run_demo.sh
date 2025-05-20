#!/bin/bash
# Demo script for the Robust Savitzky-Golay Filter library
# This script builds the library, runs the test program, and generates a plot
# Author: Mapoet
# Date: 2025-05-20

# Set the base directory to the location of this script
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
BASE_DIR="$(dirname "$SCRIPT_DIR")"

# Change to the base directory
cd "$BASE_DIR" || { echo "Error: Cannot change to base directory"; exit 1; }

echo "====== Robust Savitzky-Golay Filter Demo ======"
echo "Building the library..."

# Clean any previous build
make clean

# Build the library in release mode
make release

if [ $? -ne 0 ]; then
    echo "Error: Failed to build the library"
    exit 1
fi

# generate test data
python3 "$SCRIPT_DIR/generate_test_data.py" -n 10000 --noise 20.2 --outlier-prob 0.05 --outlier-scale 300.0 -p
mv test_data.dat ./build/
# Run the test program
echo "Running the test program..."
make test_data

if [ $? -ne 0 ]; then
    echo "Error: Failed to run the test program"
    exit 1
fi

# Check if Python and required packages are available
if ! command -v python3 &> /dev/null; then
    echo "Warning: Python 3 not found. Skipping plot generation."
    exit 0
fi

# Check for matplotlib and numpy
echo "Checking for required Python packages..."
PACKAGES_MISSING=0
for pkg in numpy matplotlib; do
    if ! python3 -c "import $pkg" &> /dev/null; then
        echo "Warning: Python package '$pkg' not found. Please install it with:"
        echo "pip install $pkg"
        PACKAGES_MISSING=1
    fi
done

if [ $PACKAGES_MISSING -eq 1 ]; then
    echo "Warning: Some required Python packages are missing. Skipping plot generation."
    exit 0
fi

# Generate plots
echo "Generating plots..."

# Full range plot
python3 "$SCRIPT_DIR/plot_results.py" -i "./build/custom_filter_results.dat" -o "./scripts/custom_filter_results.png"

# Zoomed in plot of a portion with possible outliers
python3 "$SCRIPT_DIR/plot_results.py" -i "./build/custom_filter_results.dat" -o "./scripts/custom_filter_zoomed.png" -z 80 120

echo "Demo completed successfully!"
echo "Plot files generated:"
echo "  - custom_filter_results.png (full data range)"
echo "  - custom_filter_zoomed.png (zoomed view)"
echo ""
echo "To view the plots on a system with a graphical interface, use:"
echo "python3 $SCRIPT_DIR/plot_results.py -i ./scripts/custom_filter_results.dat -s" 
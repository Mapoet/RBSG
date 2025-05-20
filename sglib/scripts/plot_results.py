#!/usr/bin/env python3
'''
Plot results from the Robust Savitzky-Golay Filter library test program
Author: Mapoet
Date: 2025-05-20
'''

import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import argparse

def main():
    # Setup command line argument parser
    parser = argparse.ArgumentParser(description='Plot Savitzky-Golay filter results')
    parser.add_argument('-i', '--input', default='../build/sg_filter_results.dat',
                        help='Input data file (default: ../build/sg_filter_results.dat)')
    parser.add_argument('-o', '--output', default=None,
                        help='Output image file (default: sg_filter_results.png)')
    parser.add_argument('-s', '--show', action='store_true',
                        help='Show plot in a window instead of saving to file')
    parser.add_argument('-z', '--zoom', nargs=2, type=int, default=None,
                        help='Zoom into specified range, e.g. "--zoom 50 100"')
    parser.add_argument('-d', '--dpi', type=int, default=150,
                        help='DPI for saved figure (default: 150)')
    args = parser.parse_args()

    # Set output file if not specified
    if args.output is None and not args.show:
        args.output = 'sg_filter_results.png'

    # Check if the input file exists
    if not os.path.isfile(args.input):
        print(f"Error: Input file '{args.input}' not found")
        return 1

    # Load the data
    try:
        data = np.loadtxt(args.input, skiprows=1)
    except Exception as e:
        print(f"Error loading data: {e}")
        return 2

    # Extract the columns
    index = data[:, 0]
    original = data[:, 1]
    noisy = data[:, 2]
    filtered = data[:, 3]
    robust = data[:, 4]

    # Calculate error metrics
    mse_standard = np.mean((original - filtered)**2)
    mse_robust = np.mean((original - robust)**2)
    
    # Create the plot
    plt.figure(figsize=(12, 9))
    
    # Main plot with all signals
    plt.subplot(2, 1, 1)
    plt.plot(index, original, 'k-', linewidth=2, label='Original Signal')
    plt.plot(index, noisy, 'r.', markersize=2, label='Noisy Data')
    plt.plot(index, filtered, 'b-', linewidth=1.5, label=f'Standard SG (MSE: {mse_standard:.4f})')
    plt.plot(index, robust, 'g-', linewidth=1.5, label=f'Robust SG (MSE: {mse_robust:.4f})')
    
    # Set zoom range if specified
    if args.zoom:
        plt.xlim(args.zoom)
    
    plt.grid(True, alpha=0.3)
    plt.title('Savitzky-Golay Filter Comparison', fontsize=14)
    plt.ylabel('Signal Amplitude', fontsize=12)
    plt.legend(loc='upper right')
    
    # Error plot
    plt.subplot(2, 1, 2)
    plt.plot(index, original - noisy, 'r-', alpha=0.7, label='Noise + Outliers')
    plt.plot(index, original - filtered, 'b-', alpha=0.7, label='Standard SG Error')
    plt.plot(index, original - robust, 'g-', alpha=0.7, label='Robust SG Error')
    
    # Set zoom range if specified
    if args.zoom:
        plt.xlim(args.zoom)
    
    plt.grid(True, alpha=0.3)
    plt.title('Error Comparison', fontsize=14)
    plt.xlabel('Index', fontsize=12)
    plt.ylabel('Error', fontsize=12)
    plt.legend(loc='upper right')
    
    # Add summary text
    plt.figtext(0.02, 0.02, f"Standard SG MSE: {mse_standard:.6f}\nRobust SG MSE: {mse_robust:.6f}", 
                fontsize=10, ha="left", bbox={"facecolor":"white", "alpha":0.8, "pad":5})
    
    plt.tight_layout()
    
    # Save or show the plot
    if args.show:
        plt.show()
    else:
        plt.savefig(args.output, dpi=args.dpi)
        print(f"Plot saved to {args.output}")
    
    return 0

if __name__ == "__main__":
    sys.exit(main()) 
#!/usr/bin/env python3
'''
Generate test data for the Robust Savitzky-Golay Filter library
This script creates synthetic signals with noise and outliers for testing
Author: Mapoet
Date: 2025-05-20
'''

import numpy as np
import matplotlib.pyplot as plt
import argparse
import os

def main():
    # Setup command line argument parser
    parser = argparse.ArgumentParser(description='Generate test data for SG filter testing')
    parser.add_argument('-n', '--points', type=int, default=200,
                        help='Number of data points (default: 200)')
    parser.add_argument('-o', '--output', default='test_data.dat',
                        help='Output data file (default: test_data.dat)')
    parser.add_argument('-p', '--plot', action='store_true',
                        help='Generate a plot of the data')
    parser.add_argument('--noise', type=float, default=0.1,
                        help='Noise amplitude factor (default: 0.1)')
    parser.add_argument('--outlier-prob', type=float, default=0.05,
                        help='Probability of outliers (default: 0.05)')
    parser.add_argument('--outlier-scale', type=float, default=5.0,
                        help='Outlier amplitude scale (default: 5.0)')
    parser.add_argument('--seed', type=int, default=None,
                        help='Random seed for reproducibility')
    args = parser.parse_args()
    
    # Set random seed if specified
    if args.seed is not None:
        np.random.seed(args.seed)
    
    # Generate x values
    x = np.linspace(0, 10, args.points)
    
    # Create test signal components
    base_signal = create_base_signal(x)
    noise = create_noise(x, args.noise)
    outliers = create_outliers(x, args.outlier_prob, args.outlier_scale)
    
    # Combine components
    clean_signal = base_signal
    noisy_signal = base_signal + noise + outliers
    
    # Save data to file
    save_data(args.output, x, clean_signal, noisy_signal)
    
    # Generate plot if requested
    if args.plot:
        plot_name = os.path.splitext(args.output)[0] + '.png'
        create_plot(x, clean_signal, noisy_signal, plot_name)
    
    print(f"Generated {args.points} data points with:")
    print(f"  - Noise amplitude factor: {args.noise}")
    print(f"  - Outlier probability: {args.outlier_prob}")
    print(f"  - Outlier scale: {args.outlier_scale}")
    print(f"Data saved to {args.output}")
    
    if args.plot:
        print(f"Plot saved to {plot_name}")
    
    return 0

def create_base_signal(x):
    """
    Create a base signal with multiple components
    """
    # Sine wave with quadratic trend
    signal = 10*np.sin(3 * x) + 0.2 * x**2*np.exp(x/max(x))
    
    # Add some higher frequency components
    signal += 0.3 * np.sin(10 * x)
    signal += 0.1 * np.sin(20 * x)
    
    # Add a gaussian peak
    peak_center = np.random.uniform(0.3, 0.7) * len(x)
    peak_width = len(x) * 0.03
    peak = 2.0 * np.exp(-0.5 * ((np.arange(len(x)) - peak_center) / peak_width)**2)
    signal += peak
    
    return signal

def create_noise(x, amplitude):
    """
    Create random noise with given amplitude
    """
    return (np.random.random(len(x)) - 0.5) * 2.0 * amplitude

def create_outliers(x, probability, scale):
    """
    Create random outliers with given probability and scale
    """
    outliers = np.zeros_like(x)
    mask = np.random.random(len(x)) < probability
    outliers[mask] = (np.random.random(np.sum(mask)) - 0.5) * 2.0 * scale
    return outliers

def save_data(filename, x, clean, noisy):
    """
    Save the generated data to a file
    """
    with open(filename, 'w') as f:
        f.write("# Index  Clean  Noisy\n")
        for i in range(len(x)):
            f.write(f"{i+1}  {clean[i]:.6f}  {noisy[i]:.6f}\n")

def create_plot(x, clean, noisy, filename):
    """
    Create a plot of the generated data
    """
    plt.figure(figsize=(10, 6))
    plt.plot(clean, 'k-', linewidth=1.5, label='Clean Signal')
    plt.plot(noisy, 'r.', markersize=2, label='Noisy Signal with Outliers')
    plt.grid(True, alpha=0.3)
    plt.title('Generated Test Data', fontsize=14)
    plt.xlabel('Index', fontsize=12)
    plt.ylabel('Amplitude', fontsize=12)
    plt.legend()
    plt.tight_layout()
    plt.savefig(filename, dpi=150)

if __name__ == "__main__":
    main() 
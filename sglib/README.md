# Robust Savitzky-Golay Filter Library

A Fortran library implementing a Savitzky-Golay filter with robust outlier detection using configurable sigma-based outlier rejection.

## Features

- Configurable window size and polynomial order
- Separation of coefficient calculation from filtering process
- Robust outlier detection with configurable sensitivity (default: 3-sigma rule)
- Minimal dependencies (only requires OpenBLAS/LAPACK)
- Simple interface for easy integration

## Requirements

- Fortran compiler supporting Fortran 2008 standard (e.g., gfortran 8+)
- OpenBLAS and LAPACK libraries

## Building the Library

```bash
cd sglib
make         # Debug build
make release # Optimized release build
```

This will build:
- Static library: `lib/libsglib.a`
- Shared library: `lib/libsglib.so`
- Test program: `build/test_sg_filter`

## Testing

Run the test program to verify the library:

```bash
make test
```

This will generate a data file `sg_filter_results.dat` with columns for the original signal, noisy signal with outliers, standard SG filter results, and robust SG filter results.

### Visualization Tools

The `scripts` directory contains Python tools for generating test data and visualizing results:

```bash
cd scripts
# Make scripts executable
chmod +x *.py run_demo.sh

# Generate test data
./generate_test_data.py -n 300 --noise 0.1 --outlier-prob 0.05 --outlier-scale 5.0 -p

# Run custom test with the generated data
make run

# Or run the complete demo
./run_demo.sh
```

This will:
1. Generate synthetic test data with noise and outliers
2. Apply both standard and robust SG filters
3. Create plots comparing the results

## Usage

### Basic Example

```fortran
program example
  use sglib_mod
  implicit none
  
  integer, parameter :: N = 100
  integer, parameter :: HALF_WINDOW = 5
  integer, parameter :: POLY_ORDER = 3
  
  real(wp), allocatable :: data(:), filtered(:), coeff(:)
  
  ! Allocate arrays
  allocate(data(N), filtered(N))
  
  ! Load or generate your data here
  call load_data(data, N)
  
  ! Method 1: Separate coefficient calculation and filtering
  call sg_compute_coeff(HALF_WINDOW, POLY_ORDER, coeff)
  call sg_filter(data, N, coeff, HALF_WINDOW, filtered)
  
  ! Method 2: Combined robust filtering with outlier removal
  ! Use default 3-sigma outlier detection
  call sg_filter_robust(data, N, HALF_WINDOW, POLY_ORDER, 3.0_wp, filtered)
  
  ! Method 3: Using more aggressive outlier detection (2-sigma)
  call sg_filter_robust(data, N, HALF_WINDOW, POLY_ORDER, 2.0_wp, filtered)
  
  ! Use filtered data
  ! ...
  
  deallocate(data, filtered)
  if (allocated(coeff)) deallocate(coeff)
  
end program example
```

## API Reference

### Coefficient Calculation

```fortran
call sg_compute_coeff(m, poly_order, coeff)
```

Computes Savitzky-Golay filter coefficients for a window of half-size `m` and polynomial order `poly_order`.

- `m`: Half window size (integer)
- `poly_order`: Polynomial order (integer)
- `coeff`: Output coefficients (real array, allocated by the subroutine)

### Standard Filtering

```fortran
call sg_filter(data, n, coeff, m, out)
```

Applies Savitzky-Golay filter with pre-computed coefficients.

- `data`: Input data array (real array)
- `n`: Size of input array (integer)
- `coeff`: Filter coefficients (real array)
- `m`: Half window size (integer)
- `out`: Filtered output data (real array)

### Robust Filtering

```fortran
call sg_filter_robust(data, n, m, poly_order, rate, out)
```

Applies robust Savitzky-Golay filter with configurable outlier detection.

- `data`: Input data array (real array)
- `n`: Size of input array (integer)
- `m`: Half window size (integer)
- `poly_order`: Polynomial order (integer)
- `rate`: Outlier detection sensitivity multiplier (real, default: 3.0)
- `out`: Filtered output data (real array)

The `rate` parameter controls the sensitivity of outlier detection. Points with deviation greater than `rate` times the standard deviation are considered outliers. Lower values (e.g., 2.0) make detection more aggressive, higher values (e.g., 4.0) more conservative.

## License

MIT License

## Acknowledgments

This implementation is based on the general Savitzky-Golay filter theory with added robust statistical techniques for outlier detection. 
!/**
!* @file main_test.f90
!* @brief Test program for robust Savitzky-Golay filtering
!* @author Mapoet
!* @version 0.1
!* @date 2025-05-20
!*/
program test_sg_filter
  use sglib_mod
  use, intrinsic :: iso_fortran_env, only: output_unit
  implicit none
  
  integer, parameter :: N = 200                ! Number of data points
  integer, parameter :: HALF_WINDOW = 7        ! Half window size (m)
  integer, parameter :: POLY_ORDER = 3         ! Polynomial order
  real(wp), parameter :: PI = 3.14159265358979323846_wp
  real(wp), parameter :: OUTLIER_PROB = 0.05_wp ! Probability of outlier
  
  real(wp), allocatable :: data(:)            ! Original data
  real(wp), allocatable :: noisy_data(:)      ! Data with noise and outliers
  real(wp), allocatable :: filtered_data(:)   ! Filtered data
  real(wp), allocatable :: robust_filtered_data(:) ! Robust filtered data
  real(wp), allocatable :: coeff(:)           ! SG filter coefficients
  
  integer :: i, seed_size
  integer, allocatable :: seed(:)
  real(wp) :: x, noise, outlier
  real(wp) :: rate=3.0_wp
  ! Initialize random number generator
  call random_seed(size=seed_size)
  allocate(seed(seed_size))
  seed = 123456789  ! Fixed seed for reproducibility
  call random_seed(put=seed)
  
  ! Allocate arrays
  allocate(data(N), noisy_data(N), filtered_data(N), robust_filtered_data(N))
  
  ! Generate a test signal: sine wave + quadratic trend
  do i = 1, N
    x = real(i, wp) / real(N, wp) * 6.0_wp * PI
    data(i) = sin(x) + 0.2_wp * x**2
  end do
  
  ! Add noise and outliers
  do i = 1, N
    call random_number(noise)
    noise = (noise - 0.5_wp) * 0.2_wp  ! Random noise between -0.1 and 0.1
    
    ! Add occasional outliers
    call random_number(outlier)
    if (outlier < OUTLIER_PROB) then
      call random_number(outlier)
      outlier = (outlier - 0.5_wp) * 5.0_wp  ! Large outlier between -2.5 and 2.5
    else
      outlier = 0.0_wp
    end if
    
    noisy_data(i) = data(i) + noise + outlier
  end do
  
  ! Compute SG filter coefficients
  call sg_compute_coeff(HALF_WINDOW, POLY_ORDER, coeff)
  
  ! Apply regular SG filter
  call sg_filter(noisy_data, N, coeff, HALF_WINDOW, filtered_data)
  
  ! Apply robust SG filter with outlier detection
  call sg_filter_robust(noisy_data, N, HALF_WINDOW, POLY_ORDER, rate, robust_filtered_data)
  
  ! Output results to file for plotting
  open(unit=10, file='sg_filter_results.dat', status='replace')
  write(10, '(A)') "# Index  Original  Noisy  Filtered  Robust"
  do i = 1, N
    write(10, '(I5, 4F12.6)') i, data(i), noisy_data(i), filtered_data(i), robust_filtered_data(i)
  end do
  close(10)
  
  write(output_unit, '(A)') "Test completed. Results written to sg_filter_results.dat"
  write(output_unit, '(A)') "Columns: Index, Original, Noisy, Filtered, Robust"
  
  ! Clean up
  deallocate(data, noisy_data, filtered_data, robust_filtered_data, coeff, seed)
  
end program test_sg_filter 
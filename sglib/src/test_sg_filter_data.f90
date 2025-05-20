!/**
!* @file custom_test.f90
!* @brief Program to test the SG filter with custom input data
!* @author Mapoet
!* @version 0.1
!* @date 2025-05-20
!*/
program custom_test
  use sglib_mod
  use, intrinsic :: iso_fortran_env, only: output_unit, error_unit
  implicit none
  
  ! Parameters
  integer, parameter :: MAX_POINTS = 10000
  integer :: half_window = 7  ! Default half window size
  integer :: poly_order = 3   ! Default polynomial order
  
  ! Data arrays
  real(wp), allocatable :: data(:)         ! Input noisy data
  real(wp), allocatable :: clean(:)        ! Original clean data (if available)
  real(wp), allocatable :: filtered(:)     ! Standard filtered data
  real(wp), allocatable :: robust(:)       ! Robust filtered data
  real(wp), allocatable :: coeff(:)        ! SG filter coefficients
  
  ! File variables
  character(len=256) :: input_filename = 'test_data.dat'
  character(len=256) :: output_filename = 'custom_filter_results.dat'
  
  ! Other variables
  integer :: n, i
  logical :: has_clean = .false.
  character(len=256) :: arg
  real(wp) :: mse_standard, mse_robust
  real(wp) :: rate=3.0_wp
  ! Process command line arguments
  i = 1
  do while (i <= command_argument_count())
    call get_command_argument(i, arg)
    
    select case (arg)
      case ('-i', '--input')
        i = i + 1
        call get_command_argument(i, input_filename)
      
      case ('-o', '--output')
        i = i + 1
        call get_command_argument(i, output_filename)
      
      case ('-w', '--window')
        i = i + 1
        call get_command_argument(i, arg)
        read(arg, *) half_window
      case ('-r', '--rate')
        i = i + 1
        call get_command_argument(i, arg)
        read(arg, *) rate
      case ('-p', '--poly')
        i = i + 1
        call get_command_argument(i, arg)
        read(arg, *) poly_order
      
      case ('-h', '--help')
        call print_usage()
        stop 0
      
      case default
        write(error_unit, '(A,A)') 'Unknown argument: ', trim(arg)
        call print_usage()
        stop 1
    end select
    
    i = i + 1
  end do
  
  ! Read input data
  write(output_unit, '(A)') 'Reading input data from ' // trim(input_filename) // '...'
  call read_data(input_filename, data, clean, n, has_clean)
  
  if (n <= 0) then
    write(error_unit, '(A)') 'Error: No data read from file'
    stop 1
  end if
  
  write(output_unit, '(A,I0,A)') 'Read ', n, ' data points'
  
  ! Allocate arrays
  allocate(filtered(n), robust(n))
  
  ! Calculate SG filter coefficients
  call sg_compute_coeff(half_window, poly_order, coeff)
  ! Apply standard SG filter
  call sg_filter(data, n, coeff, half_window, filtered)
  
  ! Apply robust SG filter
  call sg_filter_robust(data, n, half_window, poly_order,rate, robust)
  
  ! Calculate MSE if clean data is available
  if (has_clean) then
    mse_standard = sum((clean - filtered)**2) / n
    mse_robust = sum((clean - robust)**2) / n
    
    write(output_unit, '(A,F12.6)') 'Standard SG filter MSE: ', mse_standard
    write(output_unit, '(A,F12.6)') 'Robust SG filter MSE:   ', mse_robust
    
    if (mse_robust < mse_standard) then
      write(output_unit, '(A,F7.2,A)') 'Robust filter improves MSE by ', &
                                     (1.0_wp - mse_robust/mse_standard)*100.0_wp, '%'
    end if
  end if
  
  ! Write results to file
  call write_results(output_filename, n, data, filtered, robust, clean, has_clean)
  
  write(output_unit, '(A)') 'Results written to ' // trim(output_filename)
  write(output_unit, '(A)') 'Filter parameters:'
  write(output_unit, '(A,I0)') '  Window size: ', 2*half_window+1
  write(output_unit, '(A,I0)') '  Polynomial order: ', poly_order
  write(output_unit, '(A,F4.1)') '  Outlier detection rate: ', rate
  
  ! Clean up
  deallocate(data, filtered, robust, coeff)
  if (allocated(clean)) deallocate(clean)
  
contains

  !> Print usage information
  subroutine print_usage()
    write(output_unit, '(A)') 'Usage: custom_test [options]'
    write(output_unit, '(A)') 'Options:'
    write(output_unit, '(A)') '  -i, --input FILE    Input data file (default: test_data.dat)'
    write(output_unit, '(A)') '  -o, --output FILE   Output results file (default: custom_filter_results.dat)'
    write(output_unit, '(A)') '  -w, --window SIZE   Half window size (default: 7)'
    write(output_unit, '(A)') '  -p, --poly ORDER    Polynomial order (default: 3)'
    write(output_unit, '(A)') '  -r, --rate RATE     Outlier detection rate (default: 3.0)'
    write(output_unit, '(A)') '  -h, --help          Show this help message'
  end subroutine print_usage
  
  !> Read data from file
  subroutine read_data(filename, data, clean, n, has_clean)
    character(len=*), intent(in) :: filename
    real(wp), allocatable, intent(out) :: data(:), clean(:)
    integer, intent(out) :: n
    logical, intent(out) :: has_clean
    
    integer :: io_stat, line_count, col_count
    character(len=256) :: buffer, col1, col2, col3
    
    ! Count lines and detect format
    open(unit=10, file=filename, status='old', action='read', iostat=io_stat)
    if (io_stat /= 0) then
      write(error_unit, '(A)') 'Error: Cannot open input file ' // trim(filename)
      n = 0
      return
    end if
    
    ! Skip header line
    read(10, '(A)', iostat=io_stat) buffer
    
    ! Count lines
    line_count = 0
    col_count = 0
    do
      read(10, '(A)', iostat=io_stat) buffer
      if (io_stat /= 0) exit
      
      if (line_count == 0) then
        ! Determine number of columns
        col_count = count_columns(buffer)
      end if
      
      line_count = line_count + 1
      if (line_count >= MAX_POINTS) exit
    end do
    close(10)
    
    ! Allocate arrays
    n = line_count
    allocate(data(n))
    if (col_count >= 3) then
      allocate(clean(n))
      has_clean = .true.
    else
      has_clean = .false.
    end if
    
    ! Read data
    open(unit=10, file=filename, status='old', action='read')
    read(10, '(A)') buffer  ! Skip header
    
    do i = 1, n
      if (has_clean) then
        read(10, *, iostat=io_stat) col1, col2, col3
        if (io_stat /= 0) exit
        read(col2, *) clean(i)
        read(col3, *) data(i)
      else
        read(10, *, iostat=io_stat) col1, col2
        if (io_stat /= 0) exit
        read(col2, *) data(i)
      end if
    end do
    close(10)
  end subroutine read_data
  
  !> Count columns in a string
  function count_columns(line) result(count)
    character(len=*), intent(in) :: line
    integer :: count
    character(len=1) :: prev_char
    logical :: in_word
    integer :: i
    
    count = 0
    in_word = .false.
    prev_char = ' '
    
    do i = 1, len_trim(line)
      if (line(i:i) /= ' ' .and. .not. in_word) then
        in_word = .true.
        count = count + 1
      else if (line(i:i) == ' ' .and. in_word) then
        in_word = .false.
      end if
      prev_char = line(i:i)
    end do
  end function count_columns
  
  !> Write results to file
  subroutine write_results(filename, n, data, filtered, robust, clean, has_clean)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: n
    real(wp), intent(in) :: data(:), filtered(:), robust(:)
    real(wp), intent(in), optional :: clean(:)
    logical, intent(in) :: has_clean
    
    integer :: i
    
    open(unit=20, file=filename, status='replace', action='write')
    
    ! Write header
    if (has_clean) then
      write(20, '(A)') "# Index  Original  Noisy  Filtered  Robust"
    else
      write(20, '(A)') "# Index  Noisy  Filtered  Robust"
    end if
    
    ! Write data
    do i = 1, n
      if (has_clean) then
        write(20, '(I5, 4F12.6)') i, clean(i), data(i), filtered(i), robust(i)
      else
        write(20, '(I5, 3F12.6)') i, data(i), filtered(i), robust(i)
      end if
    end do
    
    close(20)
  end subroutine write_results
  
end program custom_test 
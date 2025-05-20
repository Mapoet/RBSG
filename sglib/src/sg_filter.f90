!/**
!* @file sg_filter.f90
!* @brief Implementation of Savitzky-Golay filtering with robust outlier detection
!* @author Mapoet
!* @version 0.1
!* @date 2025-05-20
!*/
submodule (sglib_mod) sg_filter_mod
  implicit none
  
contains

  ! Helper function to get window of data with boundary handling
  subroutine get_window(data, n, idx, m, window)
    real(wp), intent(in)  :: data(:)     ! Input data array
    integer, intent(in)   :: n           ! Size of input array
    integer, intent(in)   :: idx         ! Current center index
    integer, intent(in)   :: m           ! Half window size
    real(wp), intent(out) :: window(:)   ! Output window
    
    integer :: win_size, i, src_idx
    
    win_size = 2*m + 1
    if (size(window) /= win_size) then
      window = 0.0_wp
      return
    end if
    
    ! Fill window with available data
    do i = 1, win_size
      ! Compute source index: ranges from idx-m to idx+m
      src_idx = idx - m + (i-1)
      
      if (src_idx < 1) then
        ! Mirror boundary
        window(i) = data(2 - src_idx)
      else if (src_idx > n) then
        ! Mirror boundary
        window(i) = data(2*n - src_idx)
      else
        window(i) = data(src_idx)
      end if
    end do
  end subroutine get_window
  
  ! Function to detect and handle outliers
  subroutine detect_outliers(window, max_iter, valided,rate)
    real(wp), intent(in)  :: window(:)   ! Input window
    integer, intent(in)   :: max_iter    ! Maximum iterations
    logical, intent(inout) :: valided(:)
    real(wp), intent(in) :: rate
    real(wp) :: mu, sigma, threshold, sigma_min
    integer :: iter,i,n,u
    
    ! Initialize
    valided = .true.
    sigma_min = tiny(1.0_wp) * 1000.0_wp
    n=size(window)
    ! Iterative outlier detection and removal
    do iter = 1, max_iter
      ! Calculate mean and standard deviation
      u=0
      mu=0
      sigma=0
      do i=1,n
        if (valided(i)) then
          mu = mu + window(i)
          sigma = sigma + window(i)**2
          u=u+1
        end if
      end do
      mu = mu / u
      sigma = sqrt(sigma / u - mu**2)

      ! Skip if standard deviation is too small (avoid division by zero)
      if (sigma < sigma_min) exit
      ! Identify outliers (3-sigma rule)
      threshold = rate * sigma
      valided = abs(window - mu) <= threshold
      ! If no new outliers found, we're done
      if (n==u .or. u<5) exit
      n=u
    end do
  end subroutine detect_outliers
  
  ! Standard SG filter without outlier detection
  module procedure sg_filter
    integer :: i
    real(wp), allocatable :: window(:)
    
    ! Allocate working arrays
    allocate(window(2*m+1))
    
    ! Check coefficient size
    if (size(coeff) /= 2*m+1) then
      write(error_unit,*) "ERROR: Coefficient size mismatch."
      out = data
      deallocate(window)
      return
    end if
    
    ! Process each point
    do i = 1, n
      ! Get data window centered at current point
      call get_window(data, n, i, m, window)
      
      ! Apply filter coefficients
      out(i) = dot_product(window, coeff)
    end do
    
    deallocate(window)
  end procedure sg_filter
  
  ! Robust SG filter with 3-sigma outlier detection
  module procedure sg_filter_robust
    integer :: i,j, max_iter
    real(wp), allocatable :: window(:), coeff(:)
    real(wp) :: u
    logical, allocatable :: valided(:)
    
    ! Parameters
    max_iter = 3  ! Increased maximum iterations for outlier removal
    ! Compute SG filter coefficients once
    call sg_compute_coeff(m, poly_order, coeff)
    ! Allocate working arrays
    allocate(window(2*m+1), valided(2*m+1))
    
    ! Process each point
    do i = 1, n
      ! Get data window centered at current point
      call get_window(data, n, i, m, window)
      
      ! Detect and handle outliers
      call detect_outliers(window, max_iter, valided,rate)
      
      ! Apply SG filter to the cleaned window
      u=0
      do j=1,2*m+1
        if (valided(j)) then
          out(i) = out(i) + window(j)*coeff(j)
          u=u+coeff(j)
        end if
      end do
      out(i) = out(i)/u
    end do
    
    deallocate(window, valided, coeff)
  end procedure sg_filter_robust
  
end submodule sg_filter_mod 
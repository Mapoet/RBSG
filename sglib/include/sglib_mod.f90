!/**
!* @file sglib_mod.f90
!* @brief Public interfaces for Robust Savitzky-Golay filtering library
!* @author Mapoet
!* @version 0.1
!* @date 2025-05-20
!*/
module sglib_mod
  use, intrinsic :: iso_fortran_env, only: wp => real64, error_unit
  implicit none
  private

  ! Public procedures
  public :: sg_compute_coeff
  public :: sg_filter
  public :: sg_filter_robust
  public :: wp

  ! Interface definitions for public procedures
  interface
    ! Compute SG filter coefficients
    module subroutine sg_compute_coeff(m, poly_order, coeff)
      integer, intent(in) :: m            ! Half window size
      integer, intent(in) :: poly_order   ! Polynomial order
      real(wp), allocatable, intent(out) :: coeff(:) ! Output coefficients
    end subroutine sg_compute_coeff

    ! Apply SG filter without outlier detection
    module subroutine sg_filter(data, n, coeff, m, out)
      real(wp), intent(in)  :: data(:)    ! Input data array
      integer, intent(in)   :: n          ! Size of input array
      real(wp), intent(in)  :: coeff(:)   ! Filter coefficients
      integer, intent(in)   :: m          ! Half window size
      real(wp), intent(out) :: out(:)     ! Filtered output array
    end subroutine sg_filter

    ! Apply robust SG filter with 3-sigma outlier detection
    module subroutine sg_filter_robust(data, n, m, poly_order,rate,out)
      real(wp), intent(in)  :: data(:)    ! Input data array
      integer, intent(in)   :: n          ! Size of input array
      integer, intent(in)   :: m          ! Half window size
      integer, intent(in)   :: poly_order ! Polynomial order
      real(wp), intent(in) :: rate       ! Rate of outlier detection
      real(wp), intent(out) :: out(:)     ! Filtered output array
    end subroutine sg_filter_robust
  end interface

end module sglib_mod 
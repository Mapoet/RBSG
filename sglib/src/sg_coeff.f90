!/**
!* @file sg_coeff.f90
!* @brief Implementation of Savitzky-Golay filter coefficient calculation
!* @author Mapoet
!* @version 0.1
!* @date 2025-05-20
!*/
submodule (sglib_mod) sg_coeff_mod
  use, intrinsic :: iso_fortran_env, only: error_unit
  implicit none
  
  ! LAPACK interfaces
  interface
    subroutine dgelsd(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, iwork, info)
      use, intrinsic :: iso_fortran_env, only: wp => real64
      integer, intent(in) :: m, n, nrhs, lda, ldb, lwork
      integer, intent(out) :: rank, info
      integer, intent(out) :: iwork(*)
      real(wp), intent(in) :: rcond
      real(wp), intent(inout) :: a(lda,*), b(ldb,*)
      real(wp), intent(out) :: s(*), work(*)
    end subroutine dgelsd
  end interface
  
contains

  module procedure sg_compute_coeff
    integer :: npts, lwork, info, i, j, rank, liwork
    real(wp), allocatable :: A(:,:), b(:,:), s(:), work(:), x(:)
    integer, allocatable :: iwork(:)
    real(wp) :: rcond, x_i, sum_w
    
    ! 计算窗口总大小
    npts = 2*m + 1
    
    ! 分配矩阵和数组内存
    allocate(A(npts, poly_order+1), b(max(npts,poly_order+1),1), s(min(npts,poly_order+1)), x(poly_order+1))
    allocate(coeff(npts))
    
    ! 步骤1: 构造Vandermonde矩阵，x_i ∈ [-m, ..., m]
    do i = 1, npts
      x_i = real(i - m - 1, wp)  ! x_i ∈ [-m, ..., 0, ..., m]
      do j = 1, poly_order+1
        A(i, j) = x_i**(j-1)  ! A(i,j) = x_i^(j-1)
      end do
    end do
    
    ! 步骤2: 设置B = [1,0,0,...,0]^T
    b = 0.0_wp
    b(1,1) = 1.0_wp
    
    ! 设置SVD求解器的容差
    rcond = 1.0e-12_wp
    
    ! 查询工作空间大小
    allocate(work(1), iwork(1))
    call dgelsd(npts, poly_order+1, 1, A, npts, b, max(npts,poly_order+1), s, rcond, rank, work, -1, iwork, info)
    lwork = int(work(1))
    liwork = iwork(1)
    deallocate(work, iwork)
    
    ! 分配工作空间
    allocate(work(lwork), iwork(liwork))
    
    ! 步骤3: 使用DGELSD求解最小二乘问题，得到X = (A^T A)^-1 A^T B
    call dgelsd(npts, poly_order+1, 1, A, npts, b, max(npts,poly_order+1), s, rcond, rank, work, lwork, iwork, info)
    
    ! 检查错误
    if (info /= 0) then
      write(error_unit,*) "错误：DGELSD 失败，info = ", info
      coeff = 0.0_wp
      return
    end if
    
    ! 检查秩缺陷
    if (rank < poly_order+1) then
      write(error_unit,*) "警告：检测到秩亏问题。实际秩 = ", rank, &
                         " 期望秩 = ", poly_order+1
    end if
    
    ! 保存多项式系数
    do j = 1, poly_order+1
      x(j) = b(j,1)
    end do
    
    ! 步骤4: 计算权重 w_i = sum_{j=0..poly} X(j+1) * x_i^j, i=1..2m+1
    do i = 1, npts
      x_i = real(i - m - 1, wp)  ! x_i ∈ [-m, ..., 0, ..., m]
      coeff(i) = 0.0_wp
      do j = 1, poly_order+1
        coeff(i) = coeff(i) + x(j) * x_i**(j-1)
      end do
    end do
    
    ! 步骤5: 验证权重和为1
    sum_w = sum(coeff)
    if (abs(sum_w - 1.0_wp) > 1.0e-10_wp) then
      write(error_unit,'(A,F12.6,A)') "警告：权重和 = ", sum_w, " (应为 1.0)"
    end if
    
    ! 清理
    deallocate(A, b, s, work, iwork, x)
    
  end procedure sg_compute_coeff
  
end submodule sg_coeff_mod 
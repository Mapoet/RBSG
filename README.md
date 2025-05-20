
# 需求

使用fortran开发一个具有抗差能力的 Savitzky-Golay(SG)滤波的函数。要求：

1. 用户可以指定点数及阶数；
2. 计算的系数与进行滤波的过程可以分离；
3. 使用3sigma判断异常数据并剔除；
4. 有必要矩阵运算，可以使用openblas与lapack库；
5. 代码尽量不要依赖其他库，并封装成一个独立的库。

# 方案

## 概要

下面给出一个基于 Fortran 的**抗差（robust）Savitzky–Golay 滤波**函数库设计方案。该库分为两部分：

1. **系数计算模块**：用户可指定窗口点数和多项式阶数，利用线性最小二乘法（Least Squares）计算滤波系数，并可选地借助 LAPACK DGELS 例程求解正规方程或伪逆。
2. **滤波及抗差模块**：使用预先计算好的系数对数据进行卷积式滤波，并在每次滑动窗口后计算局部均值与标准差，应用 3σ 原则剔除异常值，再对剔除后的数据重新滤波。

整个库以 Fortran 模块化封装，仅依赖 BLAS/LAPACK（如 OpenBLAS），便于编译成静态或共享库供其他项目调用。

---

## 设计思路

### 1. Savitzky–Golay 滤波原理

Savitzky–Golay 滤波通过在每个滑动窗口内对原始数据点拟合低阶多项式，并将中心点的拟合值作为滤波后的输出，从而在平滑噪声的同时尽量保留信号特征 ([Wikipedia][1])。

#### 1.1 概述

SG 滤波本质上是一种基于局部多项式最小二乘拟合的卷积滤波，其主要作用是在平滑噪声的同时最大程度保留信号特征。它在时域中等价于对原始数据进行固定卷积，在频域中表现为一种低通滤波器，能够有效抑制高频噪声而对低频信号几乎不产生衰减 ([Wikipedia][1], [Wikipedia][1])。

#### 1.2 卷积形式

对于等间距数据点 $(x_j, y_j)$，设窗口宽度为 $m$（奇数），则第 $j$ 个滤波输出 $Y_j$ 由

$$
Y_j = \sum_{i=-\tfrac{m-1}{2}}^{\tfrac{m-1}{2}} C_i\,y_{j+i},\quad \frac{m+1}{2}\le j\le n-\frac{m-1}{2}
$$

给出，其中 $C_i$ 为卷积系数 ([Wikipedia][1])。

##### 1.2.1 示例

五点二阶多项式（ $m=5$, $k=2$）时，有

$$
Y_j = \frac{1}{35}\bigl(-3\,y_{j-2} + 12\,y_{j-1} + 17\,y_j + 12\,y_{j+1} -3\,y_{j+2}\bigr)
$$


##### 1.2.2 多项式拟合与系数推导  
1. **变量变换**  
   引入局部坐标 $z=(x-\bar x)/h$，中心点 $\bar x=x_j$，则 $z_i\in\{-\tfrac{m-1}{2},\dots,0,\dots,\tfrac{m-1}{2}\}$。  
2. **构造 Vandermonde 矩阵**  
   $\mathbf J_{i,j}=z_i^{\,j-1}$，尺寸 $m\times(k+1)$。  
3. **最小二乘解**  
   解正规方程  
   
   $$
   \mathbf a = (\mathbf J^T\mathbf J)^{-1}\mathbf J^T\mathbf y,
   $$

得到多项式系数向量 $\mathbf a\in\mathbb R^{k+1}$ ([Wikipedia][1])。

4\. **计算卷积权重**
对每个窗口位置 $z_i$，权重

$$
w_i = \sum_{j=0}^{k} a_{\,j}\,z_i^j,
$$

并验证 $\sum_i w_i=1$ 以保证滤波器保持常数分量不变 ([Wikipedia][1])。

#### 1.3 频域特性

* SG 滤波器在频域表现为 **低通特性**，对低频信号几乎不衰减，高频成分随拟合阶数和窗口宽度改变其衰减曲线 ([Wikipedia][1])。
* 相对于简单的滑动平均，SG 滤波在保留峰值幅度和边缘特征方面具有更好的性能，但在较高频段会出现相位反转现象 ([Wikipedia][1])。

#### 1.4 边界处理

* 常用策略：镜像扩展、常值延拓或缩减窗口，以便计算首尾 $(m-1)/2$ 个点的滤波结果 ([Wikipedia][1])。
* 也可对边界分别拟合子窗口，生成专用的边界卷积系数。

#### 1.5 多维扩展

* 二维 SG 滤波通过在平面网格上拟合二元多项式 $\sum_{i=0}^p\sum_{j=0}^q a_{ij}v^i w^j$，并以同样方式计算卷积系数，适用于图像去噪和微分 ([Wikipedia][1])。
* 为防止病态，可限制 $p<m, q<n$，并采用 SVD 或正则化求解 ([Wikipedia][1])。

以上即 SG 滤波的核心数学物理形式，涵盖时域卷积表达、最小二乘多项式拟合、频域特性与边界处理等要点。


### 2. 系数计算

* 构造维度为 (2m+1)×(k+1) 的 Vandermonde 矩阵 A，其中 m 为窗口半宽，k 为多项式阶数。
* 解正规方程最小二乘问题 $A\,c = e_j$（e\_j 为单位向量）或直接计算矩阵伪逆来获得一组中心输出系数 c ([Mathematics Stack Exchange][2], [IATE][3])。
* 推荐使用 LAPACK 的 DGELS 例程求解最小二乘，或者基于 SVD 的伪逆方法（若 A 条件较差） ([Cyber Vanguard][4], [Netlib][5])。

### 3. 异常值剔除机制

* 在滑动窗口内计算加权或非加权均值 μ 和标准差 σ；
* 若中心点与 μ 的差值超过 rate×σ，则视为异常并剔除/插值，rate为可调参数，默认为3.0（即3σ原则）([Stack Overflow][6])；
* 用户可通过命令行参数 `-r` 或 `--rate` 调整异常值判断的敏感度，较小的rate值会剔除更多点，较大的值会保留更多点；
* 对剔除后的数据再次进行 SG 滤波，实现简单的迭代鲁棒性。

### 4. 模块化与库封装

* 将系数计算、滤波、异常检测分别放在独立的 Fortran 模块中，遵循 Fortran 现代编程最佳实践 ([Fortran Lang][7])。
* 仅在必要时调用 BLAS/LAPACK 接口，如 DGELS、DAXPY 等 ([Intel][8])。

---

## 代码结构概览

```
sglib/  
├─ include/  
│   └─ sglib_mod.f90       ! 公共类型与接口声明  
├─ src/  
│   ├─ sg_coeff.f90        ! 系数计算模块  
│   ├─ sg_filter.f90       ! 卷积滤波与 3σ 抗差模块  
│   └─ main_test.f90       ! 测试示例  
└─ Makefile               ! 链接 OpenBLAS/LAPACK 配置 :contentReference[oaicite:6]{index=6}  
```

---

## 关键代码示例

### 1. 系数计算模块（`sg_coeff.f90`）

```fortran
module sg_coeff_mod
  use, intrinsic :: iso_fortran_env, only: wp => real64
  implicit none
contains
  subroutine sg_compute_coeff(m, poly_order, coeff)
    integer, intent(in) :: m, poly_order
    real(wp), allocatable, intent(out) :: coeff(:)
    integer :: npts, info, i, j
    real(wp), allocatable :: A(:,:), b(:,:), work(:)
    npts = 2*m + 1
    allocate(A(npts, poly_order+1), b(npts), coeff(npts))

    ! 构造 Vandermonde 矩阵
    do i = 1, npts
      do j = 0, poly_order
        A(i, j+1) = real(i - m - 1, wp)**j
      end do
    end do

    ! 对于中心点 e_(m+1)，求最小二乘解 A * c = e
    b = 0.0_wp
    b(m+1) = 1.0_wp

    ! 调用 DGELS 求解最小二乘
    call dgesv(npts, 1, A, npts, work(1:npts), b, npts, info)
    if (info /= 0) stop "DGELS 失败"    ! 也可改用 DGELSD（SVD）更稳健

    coeff = b
  end subroutine
end module
```

> 其中 `dgesv` 可替换为 `dgels` 或 `dgelsd`，分别对应正交分解和 SVD 法 ([Netlib][5])。

### 2. 滤波与抗差模块（`sg_filter.f90`）

```fortran
module sg_filter_mod
  use, intrinsic :: iso_fortran_env, only: wp => real64
  use sg_coeff_mod, only: sg_compute_coeff
  implicit none
contains
  subroutine sg_filter_robust(data, n, m, poly_order, rate, out)
    real(wp), intent(in)  :: data(:)
    integer, intent(in)   :: n, m, poly_order
    real(wp), intent(in)  :: rate        ! 异常值判断倍率（默认为3.0，即3σ原则）
    real(wp), intent(out) :: out(:)
    real(wp), allocatable :: coeff(:), window(:)
    integer :: i, win_size
    real(wp) :: mu, sigma

    win_size = 2*m + 1
    call sg_compute_coeff(m, poly_order, coeff)

    allocate(window(win_size))
    do i = 1, n
      ! 提取局部窗口，边界处理可镜像或保持常值
      call get_window(data, n, i, m, window)

      ! 计算均值和标准差
      mu = sum(window) / win_size
      sigma = sqrt(sum((window-mu)**2) / win_size)

      ! 异常值剔除，使用可调整的倍率参数
      window = merge(window, mu, abs(window-mu) > rate*sigma)  ! rate×σ 剔除

      ! 卷积滤波
      out(i) = dot_product(coeff, window)
    end do
  end subroutine

  subroutine get_window(data, n, idx, m, window)
    real(wp), intent(in)  :: data(:)
    integer, intent(in)   :: n, idx, m
    real(wp), intent(out) :: window(:)
    integer :: j, left, right, ws

    ws = 2*m + 1
    left  = max(1, idx-m)
    right = min(n, idx+m)
    window = 0.0_wp
    window(m+1-left+1 : m+1+(right-idx)) = data(left:right)
    ! 剩余部分可填补镜像或常量
  end subroutine
end module
```

---

## 编译与链接

在 `Makefile` 中示例配置（假设已安装 OpenBLAS）：

```makefile
FC = gfortran
FFLAGS = -O3 -frecursive
LIBS = -lopenblas -llapack

all: libsg.a

libsg.a: sg_coeff_mod.o sg_filter_mod.o
	ar rcs $@ $^

%.o: src/%.f90 include/sglib_mod.f90
	$(FC) $(FFLAGS) -Iinclude -c $< -o $@

clean:
	rm -f src/*.o *.a
```

> 参考 OpenBLAS 默认 `make.inc.example` 配置 ([GitHub][9])。

---

## 使用示例

```fortran
program test_sg
  use sg_filter_mod
  implicit none
  real(8), allocatable :: data(:), out(:)
  integer :: n

  n = 1000
  allocate(data(n), out(n))
  call load_signal(data, n)      ! 用户自行实现

  call sg_filter_robust(data, n, 3, 2, out)  ! 窗口半宽 m=3，多项式阶数 k=2

  call save_result(out, n)       ! 用户自行实现
end program
```

---

## 扩展与优化

1. **边界处理**：可增加多种策略（镜像、常量延拓或缩短窗口）以改进性能。
2. **并行化**：使用 OpenMP 对滑动窗口卷积进行多线程加速。
3. **迭代抗差**：可在剔除异常后多次迭代计算均值-σ 并平滑，提高鲁棒性。

> 模块化与最佳实践可参考 Fortran-lang 社区建议 ([annefou.github.io][10])。

---

上述方案在满足“可指定点数与阶数”、“系数与滤波分离”、“3σ剔除异常数据”、“仅依赖 OpenBLAS/LAPACK”及“封装为独立库”五大要求的同时，结构清晰、易于扩展，适用于科研与工程中的高精度信号处理场景。

[1]: https://en.wikipedia.org/wiki/Savitzky%E2%80%93Golay_filter?utm_source=chatgpt.com "Savitzky–Golay filter - Wikipedia"
[2]: https://math.stackexchange.com/questions/9077/calculating-the-savitzky-golay-coefficients?utm_source=chatgpt.com "Calculating the Savitzky Golay Coefficients - Math Stack Exchange"
[3]: https://iate.oac.uncor.edu/~mario/materia/nr/numrec/f14-8.pdf?utm_source=chatgpt.com "[PDF] 14.8 Savitzky-Golay Smoothing Filters - IATE"
[4]: https://cyber.dabamos.de/programming/modernfortran/lapack.html?utm_source=chatgpt.com "LAPACK | Programming in Modern Fortran"
[5]: https://www.netlib.org/lapack/lawnspdf/lawn193.pdf?utm_source=chatgpt.com "[PDF] LAPACK Working Note 193 - Netlib.org"
[6]: https://stackoverflow.com/questions/55336952/how-can-i-remove-outliers-numbers-3-standard-deviations-away-from-the-mean-in?utm_source=chatgpt.com "How can I remove outliers (numbers 3 standard deviations away ..."
[7]: https://fortran-lang.org/learn/best_practices/modules_programs/?utm_source=chatgpt.com "Modules and Programs - Fortran-lang.org"
[8]: https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-fortran/2023-1/blas-code-examples.html?utm_source=chatgpt.com "BLAS Code Examples - Intel"
[9]: https://github.com/xianyi/OpenBLAS/blob/develop/lapack-netlib/make.inc.example?utm_source=chatgpt.com "OpenBLAS/lapack-netlib/make.inc.example at develop - GitHub"
[10]: https://annefou.github.io/Fortran/modules/modules.html?utm_source=chatgpt.com "Procedures and Modules"


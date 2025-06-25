! Copyright (c) 2020-2021, The Neko Authors
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!
!   * Redistributions of source code must retain the above copyright
!     notice, this list of conditions and the following disclaimer.
!
!   * Redistributions in binary form must reproduce the above
!     copyright notice, this list of conditions and the following
!     disclaimer in the documentation and/or other materials provided
!     with the distribution.
!
!   * Neither the name of the authors nor the names of its
!     contributors may be used to endorse or promote products derived
!     from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
! COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
! POSSIBILITY OF SUCH DAMAGE.
!
!> Generic Gather-scatter backend for CPUs
module gs_cpu_mp
  use num_types
  use gs_bcknd_mp
  use gs_ops
  use, intrinsic :: iso_c_binding, only : c_ptr
  implicit none
  private

  !> Gather-scatter backend for CPUs
  type, public, extends(gs_bcknd_sp_t) :: gs_cpu_sp_t
   contains
     procedure, pass(this) :: init => gs_cpu_init_sp
     procedure, pass(this) :: free => gs_cpu_free_sp
     procedure, pass(this) :: gather => gs_gather_cpu_sp
     procedure, pass(this) :: scatter => gs_scatter_cpu_sp
  end type gs_cpu_sp_t

  !> Gather-scatter backend for CPUs (double precision)
  type, public, extends(gs_bcknd_dp_t) :: gs_cpu_dp_t
   contains
     procedure, pass(this) :: init => gs_cpu_init_dp
     procedure, pass(this) :: free => gs_cpu_free_dp
     procedure, pass(this) :: gather => gs_gather_cpu_dp
     procedure, pass(this) :: scatter => gs_scatter_cpu_dp
  end type gs_cpu_dp_t

contains

  !> Dummy backend initialisation
  subroutine gs_cpu_init_sp(this, nlocal, nshared, nlcl_blks, nshrd_blks)
    class(gs_cpu_sp_t), intent(inout) :: this
    integer, intent(in) :: nlocal
    integer, intent(in) :: nshared
    integer, intent(in) :: nlcl_blks
    integer, intent(in) :: nshrd_blks
  end subroutine gs_cpu_init_sp

  !> Dummy backend deallocation
  subroutine gs_cpu_free_sp(this)
    class(gs_cpu_sp_t), intent(inout) :: this
  end subroutine gs_cpu_free_sp

  !> Gather kernel
  subroutine gs_gather_cpu_sp(this, v, m, o, dg, u, n, gd, nb, b, op, shrd)
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: nb
    class(gs_cpu_sp_t), intent(inout) :: this
    real(kind=sp), dimension(m), intent(inout) :: v
    integer, dimension(m), intent(inout) :: dg
    real(kind=sp), dimension(n), intent(inout) :: u
    integer, dimension(m), intent(inout) :: gd
    integer, dimension(nb), intent(inout) :: b
    integer, intent(in) :: o
    integer, intent(in) :: op
    logical, intent(in) :: shrd

    select case(op)
    case (GS_OP_ADD)
       call gs_gather_kernel_add_sp(v, m, o, dg, u, n, gd, nb, b)
    case (GS_OP_MUL)
       call gs_gather_kernel_mul_sp(v, m, o, dg, u, n, gd, nb, b)
    case (GS_OP_MIN)
       call gs_gather_kernel_min_sp(v, m, o, dg, u, n, gd, nb, b)
    case (GS_OP_MAX)
       call gs_gather_kernel_max_sp(v, m, o, dg, u, n, gd, nb, b)
    end select

  end subroutine gs_gather_cpu_sp

  !> Gather kernel for addition of data
  !! \f$ v(dg(i)) = v(dg(i)) + u(gd(i)) \f$
  subroutine gs_gather_kernel_add_sp(v, m, o, dg, u, n, gd, nb, b)
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: nb
    real(kind=sp), dimension(m), intent(inout) :: v
    integer, dimension(m), intent(inout) :: dg
    real(kind=sp), dimension(n), intent(inout) :: u
    integer, dimension(m), intent(inout) :: gd
    integer, dimension(nb), intent(inout) :: b
    integer, intent(in) :: o
    integer :: i, j, k, blk_len
    real(kind=sp) :: tmp

    k = 0
    do i = 1, nb
       blk_len = b(i)
       tmp = u(gd(k + 1))
       do j = 2, blk_len
          tmp = tmp + u(gd(k + j))
       end do
       v(dg(k + 1)) = tmp
       k = k + blk_len
    end do

    if (o .lt. 0) then
       do i = abs(o), m
          v(dg(i)) = u(gd(i))
       end do
    else
       do i = o, m, 2
          tmp  = u(gd(i)) + u(gd(i+1))
          v(dg(i)) = tmp
       end do
    end if

  end subroutine gs_gather_kernel_add_sp

  !> Gather kernel for multiplication of data
  !! \f$ v(dg(i)) = v(dg(i)) \cdot u(gd(i)) \f$
  subroutine gs_gather_kernel_mul_sp(v, m, o, dg, u, n, gd, nb, b)
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: nb
    real(kind=sp), dimension(m), intent(inout) :: v
    integer, dimension(m), intent(inout) :: dg
    real(kind=sp), dimension(n), intent(inout) :: u
    integer, dimension(m), intent(inout) :: gd
    integer, dimension(nb), intent(inout) :: b
    integer, intent(in) :: o
    integer :: i, j, k, blk_len
    real(kind=sp) :: tmp

    k = 0
    do i = 1, nb
       blk_len = b(i)
       tmp = u(gd(k + 1))
       do j = 2, blk_len
          tmp = tmp * u(gd(k + j))
       end do
       v(dg(k + 1)) = tmp
       k = k + blk_len
    end do

    if (o .lt. 0) then
       do i = abs(o), m
          v(dg(i)) = u(gd(i))
       end do
    else
       do i = o, m, 2
          tmp  = u(gd(i)) * u(gd(i+1))
          v(dg(i)) = tmp
       end do
    end if

  end subroutine gs_gather_kernel_mul_sp

  !> Gather kernel for minimum of data
  !! \f$ v(dg(i)) = \min(v(dg(i)), u(gd(i))) \f$
  subroutine gs_gather_kernel_min_sp(v, m, o, dg, u, n, gd, nb, b)
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: nb
    real(kind=sp), dimension(m), intent(inout) :: v
    integer, dimension(m), intent(inout) :: dg
    real(kind=sp), dimension(n), intent(inout) :: u
    integer, dimension(m), intent(inout) :: gd
    integer, dimension(nb), intent(inout) :: b
    integer, intent(in) :: o
    integer :: i, j, k, blk_len
    real(kind=sp) :: tmp

    k = 0
    do i = 1, nb
       blk_len = b(i)
       tmp = u(gd(k + 1))
       do j = 2, blk_len
          tmp = min(tmp, u(gd(k + j)))
       end do
       v(dg(k + 1)) = tmp
       k = k + blk_len
    end do

    if (o .lt. 0) then
       do i = abs(o), m
          v(dg(i)) = u(gd(i))
       end do
    else
       do i = o, m, 2
          tmp  = min(u(gd(i)), u(gd(i+1)))
          v(dg(i)) = tmp
       end do
    end if

  end subroutine gs_gather_kernel_min_sp

  !> Gather kernel for maximum of data
  !! \f$ v(dg(i)) = \max(v(dg(i)), u(gd(i))) \f$
  subroutine gs_gather_kernel_max_sp(v, m, o, dg, u, n, gd, nb, b)
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: nb
    real(kind=sp), dimension(m), intent(inout) :: v
    integer, dimension(m), intent(inout) :: dg
    real(kind=sp), dimension(n), intent(inout) :: u
    integer, dimension(m), intent(inout) :: gd
    integer, dimension(nb), intent(inout) :: b
    integer, intent(in) :: o
    integer :: i, j, k, blk_len
    real(kind=sp) :: tmp

    k = 0
    do i = 1, nb
       blk_len = b(i)
       tmp = u(gd(k + 1))
       do j = 2, blk_len
          tmp = max(tmp, u(gd(k + j)))
       end do
       v(dg(k + 1)) = tmp
       k = k + blk_len
    end do

    if (o .lt. 0) then
       do i = abs(o), m
          v(dg(i)) = u(gd(i))
       end do
    else
       do i = o, m, 2
          tmp  = max(u(gd(i)), u(gd(i+1)))
          v(dg(i)) = tmp
       end do
    end if

  end subroutine gs_gather_kernel_max_sp

  !> Scatter kernel  @todo Make the kernel abstract
  subroutine gs_scatter_cpu_sp(this, v, m, dg, u, n, gd, nb, b, shrd, event)
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: nb
    class(gs_cpu_sp_t), intent(inout) :: this
    real(kind=sp), dimension(m), intent(inout) :: v
    integer, dimension(m), intent(inout) :: dg
    real(kind=sp), dimension(n), intent(inout) :: u
    integer, dimension(m), intent(inout) :: gd
    integer, dimension(nb), intent(inout) :: b
    logical, intent(in) :: shrd
    type(c_ptr) :: event

    call gs_scatter_kernel_sp(v, m, dg, u, n, gd, nb, b)

  end subroutine gs_scatter_cpu_sp

  !> Scatter kernel \f$ u(gd(i) = v(dg(i)) \f$
  subroutine gs_scatter_kernel_sp(v, m, dg, u, n, gd, nb, b)
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: nb
    real(kind=sp), dimension(m), intent(inout) :: v
    integer, dimension(m), intent(inout) :: dg
    real(kind=sp), dimension(n), intent(inout) :: u
    integer, dimension(m), intent(inout) :: gd
    integer, dimension(nb), intent(inout) :: b
    integer :: i, j, k, blk_len
    real(kind=sp) :: tmp

    k = 0
    do i = 1, nb
       blk_len = b(i)
       tmp = v(dg(k + 1))
       do j = 1, blk_len
          u(gd(k + j)) = tmp
       end do
       k = k + blk_len
    end do

    do i = k + 1, m
       u(gd(i)) = v(dg(i))
    end do

  end subroutine gs_scatter_kernel_sp

  !
  !> double precision version
  !

  !> Dummy backend initialisation
  subroutine gs_cpu_init_dp(this, nlocal, nshared, nlcl_blks, nshrd_blks)
    class(gs_cpu_dp_t), intent(inout) :: this
    integer, intent(in) :: nlocal
    integer, intent(in) :: nshared
    integer, intent(in) :: nlcl_blks
    integer, intent(in) :: nshrd_blks
  end subroutine gs_cpu_init_dp

  !> Dummy backend deallocation
  subroutine gs_cpu_free_dp(this)
    class(gs_cpu_dp_t), intent(inout) :: this
  end subroutine gs_cpu_free_dp

  !> Gather kernel
  subroutine gs_gather_cpu_dp(this, v, m, o, dg, u, n, gd, nb, b, op, shrd)
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: nb
    class(gs_cpu_dp_t), intent(inout) :: this
    real(kind=dp), dimension(m), intent(inout) :: v
    integer, dimension(m), intent(inout) :: dg
    real(kind=dp), dimension(n), intent(inout) :: u
    integer, dimension(m), intent(inout) :: gd
    integer, dimension(nb), intent(inout) :: b
    integer, intent(in) :: o
    integer, intent(in) :: op
    logical, intent(in) :: shrd

    select case(op)
    case (GS_OP_ADD)
       call gs_gather_kernel_add_dp(v, m, o, dg, u, n, gd, nb, b)
    case (GS_OP_MUL)
       call gs_gather_kernel_mul_dp(v, m, o, dg, u, n, gd, nb, b)
    case (GS_OP_MIN)
       call gs_gather_kernel_min_dp(v, m, o, dg, u, n, gd, nb, b)
    case (GS_OP_MAX)
       call gs_gather_kernel_max_dp(v, m, o, dg, u, n, gd, nb, b)
    end select

  end subroutine gs_gather_cpu_dp

  !> Gather kernel for addition of data
  !! \f$ v(dg(i)) = v(dg(i)) + u(gd(i)) \f$
  subroutine gs_gather_kernel_add_dp(v, m, o, dg, u, n, gd, nb, b)
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: nb
    real(kind=dp), dimension(m), intent(inout) :: v
    integer, dimension(m), intent(inout) :: dg
    real(kind=dp), dimension(n), intent(inout) :: u
    integer, dimension(m), intent(inout) :: gd
    integer, dimension(nb), intent(inout) :: b
    integer, intent(in) :: o
    integer :: i, j, k, blk_len
    real(kind=dp) :: tmp

    k = 0
    do i = 1, nb
       blk_len = b(i)
       tmp = u(gd(k + 1))
       do j = 2, blk_len
          tmp = tmp + u(gd(k + j))
       end do
       v(dg(k + 1)) = tmp
       k = k + blk_len
    end do

    if (o .lt. 0) then
       do i = abs(o), m
          v(dg(i)) = u(gd(i))
       end do
    else
       do i = o, m, 2
          tmp  = u(gd(i)) + u(gd(i+1))
          v(dg(i)) = tmp
       end do
    end if

  end subroutine gs_gather_kernel_add_dp

  !> Gather kernel for multiplication of data
  !! \f$ v(dg(i)) = v(dg(i)) \cdot u(gd(i)) \f$
  subroutine gs_gather_kernel_mul_dp(v, m, o, dg, u, n, gd, nb, b)
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: nb
    real(kind=dp), dimension(m), intent(inout) :: v
    integer, dimension(m), intent(inout) :: dg
    real(kind=dp), dimension(n), intent(inout) :: u
    integer, dimension(m), intent(inout) :: gd
    integer, dimension(nb), intent(inout) :: b
    integer, intent(in) :: o
    integer :: i, j, k, blk_len
    real(kind=dp) :: tmp

    k = 0
    do i = 1, nb
       blk_len = b(i)
       tmp = u(gd(k + 1))
       do j = 2, blk_len
          tmp = tmp * u(gd(k + j))
       end do
       v(dg(k + 1)) = tmp
       k = k + blk_len
    end do

    if (o .lt. 0) then
       do i = abs(o), m
          v(dg(i)) = u(gd(i))
       end do
    else
       do i = o, m, 2
          tmp  = u(gd(i)) * u(gd(i+1))
          v(dg(i)) = tmp
       end do
    end if

  end subroutine gs_gather_kernel_mul_dp

  !> Gather kernel for minimum of data
  !! \f$ v(dg(i)) = \min(v(dg(i)), u(gd(i))) \f$
  subroutine gs_gather_kernel_min_dp(v, m, o, dg, u, n, gd, nb, b)
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: nb
    real(kind=dp), dimension(m), intent(inout) :: v
    integer, dimension(m), intent(inout) :: dg
    real(kind=dp), dimension(n), intent(inout) :: u
    integer, dimension(m), intent(inout) :: gd
    integer, dimension(nb), intent(inout) :: b
    integer, intent(in) :: o
    integer :: i, j, k, blk_len
    real(kind=dp) :: tmp

    k = 0
    do i = 1, nb
       blk_len = b(i)
       tmp = u(gd(k + 1))
       do j = 2, blk_len
          tmp = min(tmp, u(gd(k + j)))
       end do
       v(dg(k + 1)) = tmp
       k = k + blk_len
    end do

    if (o .lt. 0) then
       do i = abs(o), m
          v(dg(i)) = u(gd(i))
       end do
    else
       do i = o, m, 2
          tmp  = min(u(gd(i)), u(gd(i+1)))
          v(dg(i)) = tmp
       end do
    end if

  end subroutine gs_gather_kernel_min_dp

  !> Gather kernel for maximum of data
  !! \f$ v(dg(i)) = \max(v(dg(i)), u(gd(i))) \f$
  subroutine gs_gather_kernel_max_dp(v, m, o, dg, u, n, gd, nb, b)
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: nb
    real(kind=dp), dimension(m), intent(inout) :: v
    integer, dimension(m), intent(inout) :: dg
    real(kind=dp), dimension(n), intent(inout) :: u
    integer, dimension(m), intent(inout) :: gd
    integer, dimension(nb), intent(inout) :: b
    integer, intent(in) :: o
    integer :: i, j, k, blk_len
    real(kind=dp) :: tmp

    k = 0
    do i = 1, nb
       blk_len = b(i)
       tmp = u(gd(k + 1))
       do j = 2, blk_len
          tmp = max(tmp, u(gd(k + j)))
       end do
       v(dg(k + 1)) = tmp
       k = k + blk_len
    end do

    if (o .lt. 0) then
       do i = abs(o), m
          v(dg(i)) = u(gd(i))
       end do
    else
       do i = o, m, 2
          tmp  = max(u(gd(i)), u(gd(i+1)))
          v(dg(i)) = tmp
       end do
    end if

  end subroutine gs_gather_kernel_max_dp

  !> Scatter kernel  @todo Make the kernel abstract
  subroutine gs_scatter_cpu_dp(this, v, m, dg, u, n, gd, nb, b, shrd, event)
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: nb
    class(gs_cpu_dp_t), intent(inout) :: this
    real(kind=dp), dimension(m), intent(inout) :: v
    integer, dimension(m), intent(inout) :: dg
    real(kind=dp), dimension(n), intent(inout) :: u
    integer, dimension(m), intent(inout) :: gd
    integer, dimension(nb), intent(inout) :: b
    logical, intent(in) :: shrd
    type(c_ptr) :: event

    call gs_scatter_kernel_dp(v, m, dg, u, n, gd, nb, b)

  end subroutine gs_scatter_cpu_dp

  !> Scatter kernel \f$ u(gd(i) = v(dg(i)) \f$
  subroutine gs_scatter_kernel_dp(v, m, dg, u, n, gd, nb, b)
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: nb
    real(kind=dp), dimension(m), intent(inout) :: v
    integer, dimension(m), intent(inout) :: dg
    real(kind=dp), dimension(n), intent(inout) :: u
    integer, dimension(m), intent(inout) :: gd
    integer, dimension(nb), intent(inout) :: b
    integer :: i, j, k, blk_len
    real(kind=dp) :: tmp

    k = 0
    do i = 1, nb
       blk_len = b(i)
       tmp = v(dg(k + 1))
       do j = 1, blk_len
          u(gd(k + j)) = tmp
       end do
       k = k + blk_len
    end do

    do i = k + 1, m
       u(gd(i)) = v(dg(i))
    end do

  end subroutine gs_scatter_kernel_dp

end module gs_cpu_mp

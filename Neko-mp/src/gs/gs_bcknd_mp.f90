! Copyright (c) 2021, The Neko Authors
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
!> Defines a gather-scatter backend
module gs_bcknd_mp
  use num_types
  use, intrinsic :: iso_c_binding
  implicit none
  private

  integer, public, parameter :: GS_BCKND_CPU = 1, GS_BCKND_SX = 2, &
       GS_BCKND_DEV = 3

  !> Gather-scatter backend
  type, public, abstract :: gs_bcknd_sp_t
     type(c_ptr) :: gather_event = C_NULL_PTR
     type(c_ptr) :: scatter_event = C_NULL_PTR
     type(c_ptr) :: gs_stream = C_NULL_PTR
   contains
     procedure(gs_backend_init_sp), pass(this), deferred :: init
     procedure(gs_backend_free_sp), pass(this), deferred :: free
     procedure(gs_gather_sp), pass(this), deferred :: gather
     procedure(gs_scatter_sp), pass(this), deferred :: scatter
  end type gs_bcknd_sp_t

  !> Abstract interface for initialising a Gather-Scatter backend
  abstract interface
     subroutine gs_backend_init_sp(this, nlocal, nshared, nlcl_blks, nshrd_blks)
       import gs_bcknd_sp_t
       class(gs_bcknd_sp_t), intent(inout) :: this
       integer, intent(in) :: nlocal
       integer, intent(in) :: nshared
       integer, intent(in) :: nlcl_blks
       integer, intent(in) :: nshrd_blks
     end subroutine gs_backend_init_sp
  end interface

  !> Abstract interface for deallocating a Gather-Scatter backend
  abstract interface
     subroutine gs_backend_free_sp(this)
       import gs_bcknd_sp_t
       class(gs_bcknd_sp_t), intent(inout) :: this
     end subroutine gs_backend_free_sp
  end interface

  !> Abstract interface for the Gather kernel
  !! \f$ v(dg(i)) = op(v(dg(i)), u(gd(i)) \f$
  abstract interface
     subroutine gs_gather_sp(this, v, m, o, dg, u, n, gd, nb, b, op, shrd)
       import gs_bcknd_sp_t
       import rp, sp
       integer, intent(in) :: m
       integer, intent(in) :: n
       integer, intent(in) :: nb
       class(gs_bcknd_sp_t), intent(inout) :: this
       real(kind=sp), dimension(m), intent(inout) :: v
       integer, dimension(m), intent(inout) :: dg
       real(kind=sp), dimension(n), intent(inout) :: u
       integer, dimension(m), intent(inout) :: gd
       integer, dimension(nb), intent(inout) :: b
       integer, intent(in) :: o
       integer, intent(in) :: op
       logical, intent(in) :: shrd
     end subroutine gs_gather_sp
  end interface

  !> Abstract interface for the Scatter kernel
  !! \f$ u(gd(i) = v(dg(i)) \f$
  abstract interface
     subroutine gs_scatter_sp(this, v, m, dg, u, n, gd, nb, b, shrd, event)
       import gs_bcknd_sp_t
       import c_ptr
       import rp, sp
       integer, intent(in) :: m
       integer, intent(in) :: n
       integer, intent(in) :: nb
       class(gs_bcknd_sp_t), intent(inout) :: this
       real(kind=sp), dimension(m), intent(inout) :: v
       integer, dimension(m), intent(inout) :: dg
       real(kind=sp), dimension(n), intent(inout) :: u
       integer, dimension(m), intent(inout) :: gd
       integer, dimension(nb), intent(inout) :: b
       logical, intent(in) :: shrd
       type(c_ptr) :: event
     end subroutine gs_scatter_sp
  end interface

    !> Gather-scatter backend
  type, public, abstract :: gs_bcknd_dp_t
     type(c_ptr) :: gather_event = C_NULL_PTR
     type(c_ptr) :: scatter_event = C_NULL_PTR
     type(c_ptr) :: gs_stream = C_NULL_PTR
   contains
     procedure(gs_backend_init_dp), pass(this), deferred :: init
     procedure(gs_backend_free_dp), pass(this), deferred :: free
     procedure(gs_gather_dp), pass(this), deferred :: gather
     procedure(gs_scatter_dp), pass(this), deferred :: scatter
  end type gs_bcknd_dp_t

  !> Abstract interface for initialising a Gather-Scatter backend
  abstract interface
     subroutine gs_backend_init_dp(this, nlocal, nshared, nlcl_blks, nshrd_blks)
       import gs_bcknd_dp_t
       class(gs_bcknd_dp_t), intent(inout) :: this
       integer, intent(in) :: nlocal
       integer, intent(in) :: nshared
       integer, intent(in) :: nlcl_blks
       integer, intent(in) :: nshrd_blks
     end subroutine gs_backend_init_dp
  end interface

  !> Abstract interface for deallocating a Gather-Scatter backend
  abstract interface
     subroutine gs_backend_free_dp(this)
       import gs_bcknd_dp_t
       class(gs_bcknd_dp_t), intent(inout) :: this
     end subroutine gs_backend_free_dp
  end interface

  !> Abstract interface for the Gather kernel
  !! \f$ v(dg(i)) = op(v(dg(i)), u(gd(i)) \f$
  abstract interface
     subroutine gs_gather_dp(this, v, m, o, dg, u, n, gd, nb, b, op, shrd)
       import gs_bcknd_dp_t
       import rp, dp
       integer, intent(in) :: m
       integer, intent(in) :: n
       integer, intent(in) :: nb
       class(gs_bcknd_dp_t), intent(inout) :: this
       real(kind=dp), dimension(m), intent(inout) :: v
       integer, dimension(m), intent(inout) :: dg
       real(kind=dp), dimension(n), intent(inout) :: u
       integer, dimension(m), intent(inout) :: gd
       integer, dimension(nb), intent(inout) :: b
       integer, intent(in) :: o
       integer, intent(in) :: op
       logical, intent(in) :: shrd
     end subroutine gs_gather_dp
  end interface

  !> Abstract interface for the Scatter kernel
  !! \f$ u(gd(i) = v(dg(i)) \f$
  abstract interface
     subroutine gs_scatter_dp(this, v, m, dg, u, n, gd, nb, b, shrd, event)
       import gs_bcknd_dp_t
       import c_ptr
       import rp, dp
       integer, intent(in) :: m
       integer, intent(in) :: n
       integer, intent(in) :: nb
       class(gs_bcknd_dp_t), intent(inout) :: this
       real(kind=dp), dimension(m), intent(inout) :: v
       integer, dimension(m), intent(inout) :: dg
       real(kind=dp), dimension(n), intent(inout) :: u
       integer, dimension(m), intent(inout) :: gd
       integer, dimension(nb), intent(inout) :: b
       logical, intent(in) :: shrd
       type(c_ptr) :: event
     end subroutine gs_scatter_dp
  end interface

end module gs_bcknd_mp

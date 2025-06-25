! Copyright (c) 2022, The Neko Authors
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
!> Defines a gather-scatter communication method
module gs_comm_mp
  use num_types, only : rp, sp, dp
  use comm, only : pe_size
  use stack, only : stack_i4_t
  use, intrinsic :: iso_c_binding
  implicit none
  private

  integer, public, parameter :: GS_COMM_MPI = 1, GS_COMM_MPIGPU = 2

  !> Gather-scatter communication method
  type, public, abstract :: gs_comm_sp_t
     type(stack_i4_t), allocatable :: send_dof(:)     !< Send dof to shared-gs
     type(stack_i4_t), allocatable :: recv_dof(:)     !< Recv dof to shared-gs
     integer, allocatable :: send_pe(:)               !< Send order
     integer, allocatable :: recv_pe(:)               !< Recv order
   contains
     procedure(gs_comm_init_sp), pass(this), deferred :: init
     procedure(gs_comm_free_sp), pass(this), deferred :: free
     procedure(gs_nbsend_sp), pass(this), deferred :: nbsend
     procedure(gs_nbrecv_sp), pass(this), deferred :: nbrecv
     procedure(gs_nbwait_sp), pass(this), deferred :: nbwait
     procedure, pass(this) :: init_dofs_sp
     procedure, pass(this) :: free_dofs_sp
     procedure, pass(this) :: init_order_sp
     procedure, pass(this) :: free_order_sp
  end type gs_comm_sp_t

  !> Abstract interface for initialising a Gather-scatter communication method
  abstract interface
     subroutine gs_comm_init_sp(this, send_pe, recv_pe)
       import gs_comm_sp_t
       import stack_i4_t
       class(gs_comm_sp_t), intent(inout) :: this
       type(stack_i4_t), intent(inout) :: send_pe
       type(stack_i4_t), intent(inout) :: recv_pe
     end subroutine gs_comm_init_sp
  end interface

  !> Abstract interface for deallocating a Gather-scatter communication method
  abstract interface
     subroutine gs_comm_free_sp(this)
       import gs_comm_sp_t
       class(gs_comm_sp_t), intent(inout) :: this
     end subroutine gs_comm_free_sp
  end interface

  !> Abstract interface for initiating non-blocking send operations
  abstract interface
     subroutine gs_nbsend_sp(this, u, n, deps, strm)
       import gs_comm_sp_t
       import stack_i4_t
       import c_ptr
       import sp
       class(gs_comm_sp_t), intent(inout) :: this
       integer, intent(in) :: n
       real(kind=sp), dimension(n), intent(inout) :: u
       type(c_ptr), intent(inout) :: deps
       type(c_ptr), intent(inout) :: strm
     end subroutine gs_nbsend_sp
  end interface

  !> Abstract interface for initiating non-blocking receive operations
  abstract interface
     subroutine gs_nbrecv_sp(this)
       import gs_comm_sp_t
       class(gs_comm_sp_t), intent(inout) :: this
     end subroutine gs_nbrecv_sp
  end interface

  !> Abstract interface for watining on non-blocking operations
  abstract interface
     subroutine gs_nbwait_sp(this, u, n, op, strm)
       import gs_comm_sp_t
       import stack_i4_t
       import c_ptr
       import sp
       class(gs_comm_sp_t), intent(inout) :: this
       integer, intent(in) :: n
       real(kind=sp), dimension(n), intent(inout) :: u
       integer :: op
       type(c_ptr), intent(inout) :: strm
     end subroutine gs_nbwait_sp
  end interface

    !> Gather-scatter communication method
  type, public, abstract :: gs_comm_dp_t
     type(stack_i4_t), allocatable :: send_dof(:)     !< Send dof to shared-gs
     type(stack_i4_t), allocatable :: recv_dof(:)     !< Recv dof to shared-gs
     integer, allocatable :: send_pe(:)               !< Send order
     integer, allocatable :: recv_pe(:)               !< Recv order
   contains
     procedure(gs_comm_init_dp), pass(this), deferred :: init
     procedure(gs_comm_free_dp), pass(this), deferred :: free
     procedure(gs_nbsend_dp), pass(this), deferred :: nbsend
     procedure(gs_nbrecv_dp), pass(this), deferred :: nbrecv
     procedure(gs_nbwait_dp), pass(this), deferred :: nbwait
     procedure, pass(this) :: init_dofs_dp
     procedure, pass(this) :: free_dofs_dp
     procedure, pass(this) :: init_order_dp
     procedure, pass(this) :: free_order_dp
  end type gs_comm_dp_t

  !> Abstract interface for initialising a Gather-scatter communication method
  abstract interface
     subroutine gs_comm_init_dp(this, send_pe, recv_pe)
       import gs_comm_dp_t
       import stack_i4_t
       class(gs_comm_dp_t), intent(inout) :: this
       type(stack_i4_t), intent(inout) :: send_pe
       type(stack_i4_t), intent(inout) :: recv_pe
     end subroutine gs_comm_init_dp
  end interface

  !> Abstract interface for deallocating a Gather-scatter communication method
  abstract interface
     subroutine gs_comm_free_dp(this)
       import gs_comm_dp_t
       class(gs_comm_dp_t), intent(inout) :: this
     end subroutine gs_comm_free_dp
  end interface

  !> Abstract interface for initiating non-blocking send operations
  abstract interface
     subroutine gs_nbsend_dp(this, u, n, deps, strm)
       import gs_comm_dp_t
       import stack_i4_t
       import c_ptr
       import dp
       class(gs_comm_dp_t), intent(inout) :: this
       integer, intent(in) :: n
       real(kind=dp), dimension(n), intent(inout) :: u
       type(c_ptr), intent(inout) :: deps
       type(c_ptr), intent(inout) :: strm
     end subroutine gs_nbsend_dp
  end interface

  !> Abstract interface for initiating non-blocking receive operations
  abstract interface
     subroutine gs_nbrecv_dp(this)
       import gs_comm_dp_t
       class(gs_comm_dp_t), intent(inout) :: this
     end subroutine gs_nbrecv_dp
  end interface

  !> Abstract interface for watining on non-blocking operations
  abstract interface
     subroutine gs_nbwait_dp(this, u, n, op, strm)
       import gs_comm_dp_t
       import stack_i4_t
       import c_ptr
       import dp
       class(gs_comm_dp_t), intent(inout) :: this
       integer, intent(in) :: n
       real(kind=dp), dimension(n), intent(inout) :: u
       integer :: op
       type(c_ptr), intent(inout) :: strm
     end subroutine gs_nbwait_dp
  end interface

  public :: gs_comm_init_sp, gs_comm_free_sp, gs_nbsend_sp, gs_nbrecv_sp, gs_nbwait_sp, &
            gs_comm_init_dp, gs_comm_free_dp, gs_nbsend_dp, gs_nbrecv_dp, gs_nbwait_dp
contains

  subroutine init_dofs_sp(this)
    class(gs_comm_sp_t), intent(inout) :: this
    integer :: i

    call this%free_dofs_sp()

    allocate(this%send_dof(0:pe_size-1))
    allocate(this%recv_dof(0:pe_size-1))

    do i = 0, pe_size -1
       call this%send_dof(i)%init()
       call this%recv_dof(i)%init()
    end do

  end subroutine init_dofs_sp

  subroutine free_dofs_sp(this)
    class(gs_comm_sp_t), intent(inout) :: this
    integer :: i

    if (allocated(this%send_dof)) then
       do i = 0, pe_size - 1
          call this%send_dof(i)%free()
       end do
       deallocate(this%send_dof)
    end if

    if (allocated(this%recv_dof)) then
       do i = 0, pe_size - 1
          call this%recv_dof(i)%free()
       end do
       deallocate(this%recv_dof)
    end if

  end subroutine free_dofs_sp

  subroutine init_order_sp(this, send_pe, recv_pe)
    class(gs_comm_sp_t), intent(inout) :: this
    type(stack_i4_t), intent(inout) :: send_pe
    type(stack_i4_t), intent(inout) :: recv_pe
    integer, pointer :: sp(:)
    integer :: i

    allocate(this%send_pe(send_pe%size()))

    sp => send_pe%array()
    do i = 1, send_pe%size()
       this%send_pe(i) = sp(i)
    end do

    allocate(this%recv_pe(recv_pe%size()))

    sp => recv_pe%array()
    do i = 1, recv_pe%size()
       this%recv_pe(i) = sp(i)
    end do

  end subroutine init_order_sp

  subroutine free_order_sp(this)
    class(gs_comm_sp_t), intent(inout) :: this

    if (allocated(this%send_pe)) then
       deallocate(this%send_pe)
    end if

    if (allocated(this%recv_pe)) then
       deallocate(this%recv_pe)
    end if

  end subroutine free_order_sp

  !
  !> double precision
  !

  subroutine init_dofs_dp(this)
    class(gs_comm_dp_t), intent(inout) :: this
    integer :: i

    call this%free_dofs_dp()

    allocate(this%send_dof(0:pe_size-1))
    allocate(this%recv_dof(0:pe_size-1))

    do i = 0, pe_size -1
       call this%send_dof(i)%init()
       call this%recv_dof(i)%init()
    end do

  end subroutine init_dofs_dp

  subroutine free_dofs_dp(this)
    class(gs_comm_dp_t), intent(inout) :: this
    integer :: i

    if (allocated(this%send_dof)) then
       do i = 0, pe_size - 1
          call this%send_dof(i)%free()
       end do
       deallocate(this%send_dof)
    end if

    if (allocated(this%recv_dof)) then
       do i = 0, pe_size - 1
          call this%recv_dof(i)%free()
       end do
       deallocate(this%recv_dof)
    end if

  end subroutine free_dofs_dp

  subroutine init_order_dp(this, send_pe, recv_pe)
    class(gs_comm_dp_t), intent(inout) :: this
    type(stack_i4_t), intent(inout) :: send_pe
    type(stack_i4_t), intent(inout) :: recv_pe
    integer, pointer :: sp(:)
    integer :: i

    allocate(this%send_pe(send_pe%size()))

    sp => send_pe%array()
    do i = 1, send_pe%size()
       this%send_pe(i) = sp(i)
    end do

    allocate(this%recv_pe(recv_pe%size()))

    sp => recv_pe%array()
    do i = 1, recv_pe%size()
       this%recv_pe(i) = sp(i)
    end do

  end subroutine init_order_dp

  subroutine free_order_dp(this)
    class(gs_comm_dp_t), intent(inout) :: this

    if (allocated(this%send_pe)) then
       deallocate(this%send_pe)
    end if

    if (allocated(this%recv_pe)) then
       deallocate(this%recv_pe)
    end if

  end subroutine free_order_dp

end module gs_comm_mp

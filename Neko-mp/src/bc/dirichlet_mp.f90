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
!> Defines a mixed-precision dirichlet boundary condition
module dirichlet_mp
  use num_types, only : rp, sp, dp
  use bc_mp, only : bc_dp_t, bc_sp_t
  use, intrinsic :: iso_c_binding, only : c_ptr
  implicit none
  private

  !> Mixed-precision Dirichlet boundary condition
  !! \f$ x = g \f$ on \f$\partial \Omega\f$
  type, public, extends(bc_sp_t) :: dirichlet_sp_t
     real(kind=sp), private :: g
   contains
     procedure, pass(this) :: apply_scalar => dirichlet_apply_scalar_sp
     procedure, pass(this) :: apply_vector => dirichlet_apply_vector_sp
     procedure, pass(this) :: set_g => dirichlet_set_g_sp
     !> Destructor.
     procedure, pass(this) :: free => dirichlet_free_sp
  end type dirichlet_sp_t

  type, public, extends(bc_dp_t) :: dirichlet_dp_t
     real(kind=dp), private :: g
   contains
     procedure, pass(this) :: apply_scalar => dirichlet_apply_scalar_dp
     procedure, pass(this) :: apply_vector => dirichlet_apply_vector_dp
     procedure, pass(this) :: set_g => dirichlet_set_g_dp
     !> Destructor.
     procedure, pass(this) :: free => dirichlet_free_dp
  end type dirichlet_dp_t

contains

  subroutine dirichlet_apply_scalar_sp(this, x, n, t, tstep)
    class(dirichlet_sp_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=sp), intent(inout),  dimension(n) :: x
    real(kind=sp), intent(in), optional :: t
    integer, intent(in), optional :: tstep
    integer :: i, m, k

    m = this%msk(0)
    do i = 1, m
       k = this%msk(i)
       x(k) = this%g
    end do
  end subroutine dirichlet_apply_scalar_sp

  !> Boundary condition apply for a generic Dirichlet condition
  !! to vectors @a x, @a y and @a z
  subroutine dirichlet_apply_vector_sp(this, x, y, z, n, t, tstep)
    class(dirichlet_sp_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=sp), intent(inout),  dimension(n) :: x
    real(kind=sp), intent(inout),  dimension(n) :: y
    real(kind=sp), intent(inout),  dimension(n) :: z
    real(kind=sp), intent(in), optional :: t
    integer, intent(in), optional :: tstep
    integer :: i, m, k

    m = this%msk(0)
    do i = 1, m
       k = this%msk(i)
       x(k) = this%g
       y(k) = this%g
       z(k) = this%g
    end do

  end subroutine dirichlet_apply_vector_sp

  !> Set value of \f$ g \f$
  subroutine dirichlet_set_g_sp(this, g)
    class(dirichlet_sp_t), intent(inout) :: this
    real(kind=sp), intent(in) :: g

    this%g = g

  end subroutine dirichlet_set_g_sp

  !> Destructor
  subroutine dirichlet_free_sp(this)
    class(dirichlet_sp_t), target, intent(inout) :: this

    call this%free_base

  end subroutine dirichlet_free_sp

  !
  !> double precision 
  !

  subroutine dirichlet_apply_scalar_dp(this, x, n, t, tstep)
    class(dirichlet_dp_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=dp), intent(inout),  dimension(n) :: x
    real(kind=dp), intent(in), optional :: t
    integer, intent(in), optional :: tstep
    integer :: i, m, k

    m = this%msk(0)
    do i = 1, m
       k = this%msk(i)
       x(k) = this%g
    end do
  end subroutine dirichlet_apply_scalar_dp

  !> Boundary condition apply for a generic Dirichlet condition
  !! to vectors @a x, @a y and @a z
  subroutine dirichlet_apply_vector_dp(this, x, y, z, n, t, tstep)
    class(dirichlet_dp_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=dp), intent(inout),  dimension(n) :: x
    real(kind=dp), intent(inout),  dimension(n) :: y
    real(kind=dp), intent(inout),  dimension(n) :: z
    real(kind=dp), intent(in), optional :: t
    integer, intent(in), optional :: tstep
    integer :: i, m, k

    m = this%msk(0)
    do i = 1, m
       k = this%msk(i)
       x(k) = this%g
       y(k) = this%g
       z(k) = this%g
    end do

  end subroutine dirichlet_apply_vector_dp

  !> Set value of \f$ g \f$
  subroutine dirichlet_set_g_dp(this, g)
    class(dirichlet_dp_t), intent(inout) :: this
    real(kind=dp), intent(in) :: g

    this%g = g

  end subroutine dirichlet_set_g_dp

  !> Destructor
  subroutine dirichlet_free_dp(this)
    class(dirichlet_dp_t), target, intent(inout) :: this

    call this%free_base

  end subroutine dirichlet_free_dp

end module dirichlet_mp

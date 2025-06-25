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
  use bc_mp, only : bc_mp_t
  use, intrinsic :: iso_c_binding, only : c_ptr
  implicit none
  private

  !> Mixed-precision Dirichlet boundary condition
  !! \f$ x = g \f$ on \f$\partial \Omega\f$
  type, public, extends(bc_mp_t) :: dirichlet_mp_t
     real(kind=dp), private :: g_dp
     real(kind=sp), private :: g_sp
   contains
     procedure :: apply_scalar_sp => dirichlet_apply_scalar_sp
     procedure :: apply_scalar_dp => dirichlet_apply_scalar_dp
     procedure :: apply_vector_sp => dirichlet_apply_vector_sp
     procedure :: apply_vector_dp => dirichlet_apply_vector_dp
     procedure :: set_g => dirichlet_set_g
     !> Destructor.
     procedure :: free => dirichlet_free
  end type dirichlet_mp_t

contains

  !> Boundary condition apply for Dirichlet condition to a vector @a x (SP)
  subroutine dirichlet_apply_scalar_sp(this, x, n, t, tstep)
    class(dirichlet_mp_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=sp), intent(inout), dimension(n) :: x
    real(kind=sp), intent(in), dimension(n) :: t
    integer, intent(in), optional :: tstep
    integer :: i, k

    do i = 1, size(this%msk)
       k = this%msk(i)
       x(k) = this%g_sp
    end do
  end subroutine dirichlet_apply_scalar_sp

  !> Boundary condition apply for Dirichlet condition to a vector @a x (DP)
  subroutine dirichlet_apply_scalar_dp(this, x, n, t, tstep)
    class(dirichlet_mp_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=dp), intent(inout), dimension(n) :: x
    real(kind=dp), intent(in), dimension(n) :: t
    integer, intent(in), optional :: tstep
    integer :: i, k

    do i = 1, size(this%msk)
       k = this%msk(i)
       x(k) = this%g_dp
    end do
  end subroutine dirichlet_apply_scalar_dp

  !> Boundary condition apply for Dirichlet condition to vectors @a x, @a y and @a z (SP)
  subroutine dirichlet_apply_vector_sp(this, x, y, z, n, t, tstep)
    class(dirichlet_mp_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=sp), intent(inout), dimension(n) :: x, y, z
    real(kind=sp), intent(in), dimension(n) :: t
    integer, intent(in), optional :: tstep
    integer :: i, k

    do i = 1, size(this%msk)
       k = this%msk(i)
       x(k) = this%g_sp
       y(k) = this%g_sp
       z(k) = this%g_sp
    end do
  end subroutine dirichlet_apply_vector_sp

  !> Boundary condition apply for Dirichlet condition to vectors @a x, @a y and @a z (DP)
  subroutine dirichlet_apply_vector_dp(this, x, y, z, n, t, tstep)
    class(dirichlet_mp_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=dp), intent(inout), dimension(n) :: x, y, z
    real(kind=dp), intent(in), dimension(n) :: t
    integer, intent(in), optional :: tstep
    integer :: i, k

    do i = 1, size(this%msk)
       k = this%msk(i)
       x(k) = this%g_dp
       y(k) = this%g_dp
       z(k) = this%g_dp
    end do
  end subroutine dirichlet_apply_vector_dp

  !> Set value of \f$ g \f$
  subroutine dirichlet_set_g(this, g)
    class(dirichlet_mp_t), intent(inout) :: this
    real(kind=dp), intent(in) :: g

    this%g_dp = g
    this%g_sp = real(g, kind=sp)
  end subroutine dirichlet_set_g

  !> Destructor
  subroutine dirichlet_free(this)
    class(dirichlet_mp_t), intent(inout) :: this

    call this%free_base()
  end subroutine dirichlet_free

end module dirichlet_mp

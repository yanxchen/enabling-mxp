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
!> Jacobi preconditioner
module jacobi_mp
  use math
  use precon_mp
  use num_types, only : rp, sp, dp
  use coef_mp, only : coef_dp_t
  use dofmap_mp, only : dofmap_dp_t
  use gather_scatter_mp, only : gs_dp_t, GS_OP_ADD
  implicit none
  private

  !> Defines a jacobi preconditioner
  !> dummy mp, actually sp
  type, public, extends(pc_mp_t) :: jacobi_mp_t
     real(kind=sp), allocatable :: d_sp(:,:,:,:)
     real(kind=dp), allocatable :: d_dp(:,:,:,:)
     type(gs_dp_t), pointer :: gs_h
     type(dofmap_dp_t), pointer :: dof
     type(coef_dp_t), pointer :: coef
   contains
     procedure, pass(this) :: init => jacobi_init_mp
     procedure, pass(this) :: free => jacobi_free_mp
     procedure, pass(this) :: solve_sp => jacobi_solve_sp
     procedure, pass(this) :: solve_dp => jacobi_solve_dp
     procedure, pass(this) :: update => jacobi_update_dp
  end type jacobi_mp_t

contains

  subroutine jacobi_init_mp(this, coef, dof, gs_h)
    class(jacobi_mp_t), intent(inout) :: this
    type(coef_dp_t), intent(inout), target :: coef
    type(dofmap_dp_t), intent(inout), target :: dof
    type(gs_dp_t), intent(inout), target :: gs_h

    call this%free()
    this%gs_h => gs_h
    this%dof => dof
    this%coef => coef
    allocate(this%d_sp(dof%Xh%lx,dof%Xh%ly,dof%Xh%lz, dof%msh%nelv))
    allocate(this%d_dp(dof%Xh%lx,dof%Xh%ly,dof%Xh%lz, dof%msh%nelv))
    
    call jacobi_update_dp(this)
    this%d_sp = real(this%d_dp, kind=sp)

  end subroutine jacobi_init_mp

  subroutine jacobi_free_mp(this)
    class(jacobi_mp_t), intent(inout) :: this
    if (allocated(this%d_sp)) then
       deallocate(this%d_sp)
    end if
    if (allocated(this%d_dp)) then
       deallocate(this%d_dp)
    end if
    nullify(this%dof)
  end subroutine jacobi_free_mp

  !> The jacobi preconditioner \f$ J z = r \f$
  !! \f$ z = J^{-1}r\f$ where \f$ J^{-1} ~= 1/diag(A) \f$
  subroutine jacobi_solve_sp(this, z, r, n)
    integer, intent(in) :: n
    class(jacobi_mp_t), intent(inout) :: this
    real(kind=sp), dimension(n), intent(inout) :: z
    real(kind=sp), dimension(n), intent(inout) :: r

    call col3_sp(z,r,this%d_sp,n)
  end subroutine jacobi_solve_sp

  subroutine jacobi_solve_dp(this, z, r, n)
    integer, intent(in) :: n
    class(jacobi_mp_t), intent(inout) :: this
    real(kind=dp), dimension(n), intent(inout) :: z
    real(kind=dp), dimension(n), intent(inout) :: r
    call col3(z,r,this%d_dp,n)
  end subroutine jacobi_solve_dp

  !> Update Jacobi preconditioner if the geometry G has changed
  subroutine jacobi_update_dp(this)
    class(jacobi_mp_t), intent(inout) :: this
    associate(dof => this%dof, coef => this%coef, gs_h => this%gs_h)

      select case(dof%Xh%lx)
      case (14)
         call jacobi_update_lx14_dp(this%d_dp, dof%Xh%dxt, dof%Xh%dyt, dof%Xh%dzt, &
              coef%G11, coef%G22, coef%G33, coef%G12, coef%G13, coef%G23, &
              dof%msh%dfrmd_el, dof%msh%nelv)
      case (13)
         call jacobi_update_lx13_dp(this%d_dp, dof%Xh%dxt, dof%Xh%dyt, dof%Xh%dzt, &
              coef%G11, coef%G22, coef%G33, coef%G12, coef%G13, coef%G23, &
              dof%msh%dfrmd_el, dof%msh%nelv)
      case (12)
         call jacobi_update_lx12_dp(this%d_dp, dof%Xh%dxt, dof%Xh%dyt, dof%Xh%dzt, &
              coef%G11, coef%G22, coef%G33, coef%G12, coef%G13, coef%G23, &
              dof%msh%dfrmd_el, dof%msh%nelv)
      case (11)
         call jacobi_update_lx11_dp(this%d_dp, dof%Xh%dxt, dof%Xh%dyt, dof%Xh%dzt, &
              coef%G11, coef%G22, coef%G33, coef%G12, coef%G13, coef%G23, &
              dof%msh%dfrmd_el, dof%msh%nelv)
      case (10)
         call jacobi_update_lx10_dp(this%d_dp, dof%Xh%dxt, dof%Xh%dyt, dof%Xh%dzt, &
              coef%G11, coef%G22, coef%G33, coef%G12, coef%G13, coef%G23, &
              dof%msh%dfrmd_el, dof%msh%nelv)
      case (9)
         call jacobi_update_lx9_dp(this%d_dp, dof%Xh%dxt, dof%Xh%dyt, dof%Xh%dzt, &
              coef%G11, coef%G22, coef%G33, coef%G12, coef%G13, coef%G23, &
              dof%msh%dfrmd_el, dof%msh%nelv)
      case (8)
         call jacobi_update_lx8_dp(this%d_dp, dof%Xh%dxt, dof%Xh%dyt, dof%Xh%dzt, &
              coef%G11, coef%G22, coef%G33, coef%G12, coef%G13, coef%G23, &
              dof%msh%dfrmd_el, dof%msh%nelv)
      case (7)
         call jacobi_update_lx7_dp(this%d_dp, dof%Xh%dxt, dof%Xh%dyt, dof%Xh%dzt, &
              coef%G11, coef%G22, coef%G33, coef%G12, coef%G13, coef%G23, &
              dof%msh%dfrmd_el, dof%msh%nelv)
      case (6)
         call jacobi_update_lx6_dp(this%d_dp, dof%Xh%dxt, dof%Xh%dyt, dof%Xh%dzt, &
              coef%G11, coef%G22, coef%G33, coef%G12, coef%G13, coef%G23, &
              dof%msh%dfrmd_el, dof%msh%nelv)
      case (5)
         call jacobi_update_lx5_dp(this%d_dp, dof%Xh%dxt, dof%Xh%dyt, dof%Xh%dzt, &
              coef%G11, coef%G22, coef%G33, coef%G12, coef%G13, coef%G23, &
              dof%msh%dfrmd_el, dof%msh%nelv)
      case (4)
         call jacobi_update_lx4_dp(this%d_dp, dof%Xh%dxt, dof%Xh%dyt, dof%Xh%dzt, &
              coef%G11, coef%G22, coef%G33, coef%G12, coef%G13, coef%G23, &
              dof%msh%dfrmd_el, dof%msh%nelv)
      case (3)
         call jacobi_update_lx3_dp(this%d_dp, dof%Xh%dxt, dof%Xh%dyt, dof%Xh%dzt, &
              coef%G11, coef%G22, coef%G33, coef%G12, coef%G13, coef%G23, &
              dof%msh%dfrmd_el, dof%msh%nelv)
      case (2)
         call jacobi_update_lx2_dp(this%d_dp, dof%Xh%dxt, dof%Xh%dyt, dof%Xh%dzt, &
              coef%G11, coef%G22, coef%G33, coef%G12, coef%G13, coef%G23, &
              dof%msh%dfrmd_el, dof%msh%nelv)
      case default
         call jacobi_update_lx_dp(this%d_dp, dof%Xh%dxt, dof%Xh%dyt, dof%Xh%dzt, &
              coef%G11, coef%G22, coef%G33, coef%G12, coef%G13, coef%G23, &
              dof%msh%dfrmd_el, dof%msh%nelv, dof%Xh%lx)
      end select

      call col2(this%d_dp,coef%h1,coef%dof%size())
      if (coef%ifh2) call addcol3(this%d_dp,coef%h2,coef%B,coef%dof%size())
      call gs_h%op(this%d_dp, dof%size(), GS_OP_ADD)
      call invcol1(this%d_dp,dof%size())
    end associate
  end subroutine jacobi_update_dp

  !> Generic CPU kernel for updating the Jacobi preconditioner
  subroutine jacobi_update_lx_dp(d, dxt, dyt, dzt, G11, G22, G33, &
                              G12, G13, G23, dfrmd_el, n, lx)
    integer, intent(in) :: n, lx
    real(kind=dp), intent(inout) :: d(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G11(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G22(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G33(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G12(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G13(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G23(lx, lx, lx, n)
    real(kind=dp), intent(in) :: dxt(lx, lx)
    real(kind=dp), intent(in) :: dyt(lx, lx)
    real(kind=dp), intent(in) :: dzt(lx, lx)
    logical, intent(in) :: dfrmd_el(n)
    integer :: i, j, k, l, e
      
    d = 0d0

    do e = 1,n
       do l = 1,lx
          do k = 1,lx
             do j = 1,lx
                do i = 1,lx
                   d(i,j,k,e) = d(i,j,k,e) + &
                        G11(l,j,k,e) * dxt(i,l)**2
                end do
             end do
          end do
       end do
       do l = 1,lx
          do k = 1,lx
             do j = 1,lx
                do i = 1,lx
                   d(i,j,k,e) = d(i,j,k,e) + &
                        G22(i,l,k,e) * dyt(j,l)**2
                end do
             end do
          end do
       end do
       do l = 1,lx
          do k = 1,lx
             do j = 1,lx
                do i = 1,lx
                   d(i,j,k,e) = d(i,j,k,e) + &
                        G33(i,j,l,e) * dzt(k,l)**2
                end do
             end do
          end do
       end do
       
       if (dfrmd_el(e)) then
          do j = 1,lx,lx-1
             do k = 1,lx,lx-1
                d(1,j,k,e) = d(1,j,k,e) &
                     + G12(1,j,k,e) * dxt(1,1)*dyt(j,j) &
                     + G13(1,j,k,e) * dxt(1,1)*dzt(k,k)
                d(lx,j,k,e) = d(lx,j,k,e) &
                     + G12(lx,j,k,e) * dxt(lx,lx)*dyt(j,j) &
                     + G13(lx,j,k,e) * dxt(lx,lx)*dzt(k,k)
             end do
          end do
          
          do i = 1,lx,lx-1
             do k = 1,lx,lx-1
                d(i,1,k,e) = d(i,1,k,e) &
                     + G12(i,1,k,e) * dyt(1,1)*dxt(i,i) &
                     + G23(i,1,k,e) * dyt(1,1)*dzt(k,k)
                d(i,lx,k,e) = d(i,lx,k,e) &
                     + G12(i,lx,k,e) * dyt(lx,lx)*dxt(i,i) &
                     + G23(i,lx,k,e) * dyt(lx,lx)*dzt(k,k)
             end do
          end do
          do i = 1,lx,lx-1
             do j = 1,lx,lx-1
                d(i,j,1,e) = d(i,j,1,e) &
                     + G13(i,j,1,e) * dzt(1,1)*dxt(i,i) &
                     + G23(i,j,1,e) * dzt(1,1)*dyt(j,j)
                d(i,j,lx,e) = d(i,j,lx,e) &
                     + G13(i,j,lx,e) * dzt(lx,lx)*dxt(i,i) &
                     + G23(i,j,lx,e) * dzt(lx,lx)*dyt(j,j)
             end do
          end do
       end if
    end do
  end subroutine jacobi_update_lx_dp

  subroutine jacobi_update_lx14_dp(d, dxt, dyt, dzt, G11, G22, G33, &
                                G12, G13, G23, dfrmd_el, n)
    integer, parameter :: lx = 14
    integer, intent(in) :: n
    real(kind=dp), intent(inout) :: d(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G11(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G22(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G33(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G12(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G13(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G23(lx, lx, lx, n)
    real(kind=dp), intent(in) :: dxt(lx, lx)
    real(kind=dp), intent(in) :: dyt(lx, lx)
    real(kind=dp), intent(in) :: dzt(lx, lx)
    logical, intent(in) :: dfrmd_el(n)
    integer :: i, j, k, l, e
      
    d = 0d0

    do e = 1,n
       do l = 1,lx
          do k = 1,lx
             do j = 1,lx
                do i = 1,lx
                   d(i,j,k,e) = d(i,j,k,e) + &
                        G11(l,j,k,e) * dxt(i,l)**2
                end do
             end do
          end do
       end do
       do l = 1,lx
          do k = 1,lx
             do j = 1,lx
                do i = 1,lx
                   d(i,j,k,e) = d(i,j,k,e) + &
                        G22(i,l,k,e) * dyt(j,l)**2
                end do
             end do
          end do
       end do
       do l = 1,lx
          do k = 1,lx
             do j = 1,lx
                do i = 1,lx
                   d(i,j,k,e) = d(i,j,k,e) + &
                        G33(i,j,l,e) * dzt(k,l)**2
                end do
             end do
          end do
       end do
       
       if (dfrmd_el(e)) then
          do j = 1,lx,lx-1
             do k = 1,lx,lx-1
                d(1,j,k,e) = d(1,j,k,e) &
                     + G12(1,j,k,e) * dxt(1,1)*dyt(j,j) &
                     + G13(1,j,k,e) * dxt(1,1)*dzt(k,k)
                d(lx,j,k,e) = d(lx,j,k,e) &
                     + G12(lx,j,k,e) * dxt(lx,lx)*dyt(j,j) &
                     + G13(lx,j,k,e) * dxt(lx,lx)*dzt(k,k)
             end do
          end do
          
          do i = 1,lx,lx-1
             do k = 1,lx,lx-1
                d(i,1,k,e) = d(i,1,k,e) &
                     + G12(i,1,k,e) * dyt(1,1)*dxt(i,i) &
                     + G23(i,1,k,e) * dyt(1,1)*dzt(k,k)
                d(i,lx,k,e) = d(i,lx,k,e) &
                     + G12(i,lx,k,e) * dyt(lx,lx)*dxt(i,i) &
                     + G23(i,lx,k,e) * dyt(lx,lx)*dzt(k,k)
             end do
          end do
          do i = 1,lx,lx-1
             do j = 1,lx,lx-1
                d(i,j,1,e) = d(i,j,1,e) &
                     + G13(i,j,1,e) * dzt(1,1)*dxt(i,i) &
                     + G23(i,j,1,e) * dzt(1,1)*dyt(j,j)
                d(i,j,lx,e) = d(i,j,lx,e) &
                     + G13(i,j,lx,e) * dzt(lx,lx)*dxt(i,i) &
                     + G23(i,j,lx,e) * dzt(lx,lx)*dyt(j,j)
             end do
          end do
       end if
    end do
  end subroutine jacobi_update_lx14_dp
  
  subroutine jacobi_update_lx13_dp(d, dxt, dyt, dzt, G11, G22, G33, &
                              G12, G13, G23, dfrmd_el, n)
    integer, parameter :: lx = 13
    integer, intent(in) :: n
    real(kind=dp), intent(inout) :: d(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G11(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G22(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G33(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G12(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G13(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G23(lx, lx, lx, n)
    real(kind=dp), intent(in) :: dxt(lx, lx)
    real(kind=dp), intent(in) :: dyt(lx, lx)
    real(kind=dp), intent(in) :: dzt(lx, lx)
    logical, intent(in) :: dfrmd_el(n)
    integer :: i, j, k, l, e
    
    d = 0d0

    do e = 1,n
       do l = 1,lx
          do k = 1,lx
             do j = 1,lx
                do i = 1,lx
                   d(i,j,k,e) = d(i,j,k,e) + &
                        G11(l,j,k,e) * dxt(i,l)**2
                end do
             end do
          end do
       end do
       do l = 1,lx
          do k = 1,lx
             do j = 1,lx
                do i = 1,lx
                   d(i,j,k,e) = d(i,j,k,e) + &
                        G22(i,l,k,e) * dyt(j,l)**2
                end do
             end do
          end do
       end do
       do l = 1,lx
          do k = 1,lx
             do j = 1,lx
                do i = 1,lx
                   d(i,j,k,e) = d(i,j,k,e) + &
                        G33(i,j,l,e) * dzt(k,l)**2
                end do
             end do
          end do
       end do
       
       if (dfrmd_el(e)) then
          do j = 1,lx,lx-1
             do k = 1,lx,lx-1
                d(1,j,k,e) = d(1,j,k,e) &
                     + G12(1,j,k,e) * dxt(1,1)*dyt(j,j) &
                     + G13(1,j,k,e) * dxt(1,1)*dzt(k,k)
                d(lx,j,k,e) = d(lx,j,k,e) &
                     + G12(lx,j,k,e) * dxt(lx,lx)*dyt(j,j) &
                     + G13(lx,j,k,e) * dxt(lx,lx)*dzt(k,k)
             end do
          end do
          
          do i = 1,lx,lx-1
             do k = 1,lx,lx-1
                d(i,1,k,e) = d(i,1,k,e) &
                     + G12(i,1,k,e) * dyt(1,1)*dxt(i,i) &
                     + G23(i,1,k,e) * dyt(1,1)*dzt(k,k)
                d(i,lx,k,e) = d(i,lx,k,e) &
                     + G12(i,lx,k,e) * dyt(lx,lx)*dxt(i,i) &
                     + G23(i,lx,k,e) * dyt(lx,lx)*dzt(k,k)
             end do
          end do
          do i = 1,lx,lx-1
             do j = 1,lx,lx-1
                d(i,j,1,e) = d(i,j,1,e) &
                     + G13(i,j,1,e) * dzt(1,1)*dxt(i,i) &
                     + G23(i,j,1,e) * dzt(1,1)*dyt(j,j)
                d(i,j,lx,e) = d(i,j,lx,e) &
                     + G13(i,j,lx,e) * dzt(lx,lx)*dxt(i,i) &
                     + G23(i,j,lx,e) * dzt(lx,lx)*dyt(j,j)
             end do
          end do
       end if
    end do
  end subroutine jacobi_update_lx13_dp
  
  subroutine jacobi_update_lx12_dp(d, dxt, dyt, dzt, G11, G22, G33, &
                                G12, G13, G23, dfrmd_el, n)
    integer, parameter :: lx = 12
    integer, intent(in) :: n
    real(kind=dp), intent(inout) :: d(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G11(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G22(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G33(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G12(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G13(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G23(lx, lx, lx, n)
    real(kind=dp), intent(in) :: dxt(lx, lx)
    real(kind=dp), intent(in) :: dyt(lx, lx)
    real(kind=dp), intent(in) :: dzt(lx, lx)
    logical, intent(in) :: dfrmd_el(n)
    integer :: i, j, k, l, e
      
    d = 0d0

    do e = 1,n
       do l = 1,lx
          do k = 1,lx
             do j = 1,lx
                do i = 1,lx
                   d(i,j,k,e) = d(i,j,k,e) + &
                        G11(l,j,k,e) * dxt(i,l)**2
                end do
             end do
          end do
       end do
       do l = 1,lx
          do k = 1,lx
             do j = 1,lx
                do i = 1,lx
                   d(i,j,k,e) = d(i,j,k,e) + &
                        G22(i,l,k,e) * dyt(j,l)**2
                end do
             end do
          end do
       end do
       do l = 1,lx
          do k = 1,lx
             do j = 1,lx
                do i = 1,lx
                   d(i,j,k,e) = d(i,j,k,e) + &
                        G33(i,j,l,e) * dzt(k,l)**2
                end do
             end do
          end do
       end do
       
       if (dfrmd_el(e)) then
          do j = 1,lx,lx-1
             do k = 1,lx,lx-1
                d(1,j,k,e) = d(1,j,k,e) &
                     + G12(1,j,k,e) * dxt(1,1)*dyt(j,j) &
                     + G13(1,j,k,e) * dxt(1,1)*dzt(k,k)
                d(lx,j,k,e) = d(lx,j,k,e) &
                     + G12(lx,j,k,e) * dxt(lx,lx)*dyt(j,j) &
                     + G13(lx,j,k,e) * dxt(lx,lx)*dzt(k,k)
             end do
          end do
          
          do i = 1,lx,lx-1
             do k = 1,lx,lx-1
                d(i,1,k,e) = d(i,1,k,e) &
                     + G12(i,1,k,e) * dyt(1,1)*dxt(i,i) &
                     + G23(i,1,k,e) * dyt(1,1)*dzt(k,k)
                d(i,lx,k,e) = d(i,lx,k,e) &
                     + G12(i,lx,k,e) * dyt(lx,lx)*dxt(i,i) &
                     + G23(i,lx,k,e) * dyt(lx,lx)*dzt(k,k)
             end do
          end do
          do i = 1,lx,lx-1
             do j = 1,lx,lx-1
                d(i,j,1,e) = d(i,j,1,e) &
                     + G13(i,j,1,e) * dzt(1,1)*dxt(i,i) &
                     + G23(i,j,1,e) * dzt(1,1)*dyt(j,j)
                d(i,j,lx,e) = d(i,j,lx,e) &
                     + G13(i,j,lx,e) * dzt(lx,lx)*dxt(i,i) &
                     + G23(i,j,lx,e) * dzt(lx,lx)*dyt(j,j)
             end do
          end do
       end if
    end do
  end subroutine jacobi_update_lx12_dp
  
  subroutine jacobi_update_lx11_dp(d, dxt, dyt, dzt, G11, G22, G33, &
                                G12, G13, G23, dfrmd_el, n)
    integer, parameter :: lx = 11
    integer, intent(in) :: n
    real(kind=dp), intent(inout) :: d(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G11(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G22(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G33(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G12(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G13(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G23(lx, lx, lx, n)
    real(kind=dp), intent(in) :: dxt(lx, lx)
    real(kind=dp), intent(in) :: dyt(lx, lx)
    real(kind=dp), intent(in) :: dzt(lx, lx)
    logical, intent(in) :: dfrmd_el(n)
    integer :: i, j, k, l, e
      
    d = 0d0

    do e = 1,n
       do l = 1,lx
          do k = 1,lx
             do j = 1,lx
                do i = 1,lx
                   d(i,j,k,e) = d(i,j,k,e) + &
                        G11(l,j,k,e) * dxt(i,l)**2
                end do
             end do
          end do
       end do
       do l = 1,lx
          do k = 1,lx
             do j = 1,lx
                do i = 1,lx
                   d(i,j,k,e) = d(i,j,k,e) + &
                        G22(i,l,k,e) * dyt(j,l)**2
                end do
             end do
          end do
       end do
       do l = 1,lx
          do k = 1,lx
             do j = 1,lx
                do i = 1,lx
                   d(i,j,k,e) = d(i,j,k,e) + &
                        G33(i,j,l,e) * dzt(k,l)**2
                end do
             end do
          end do
       end do
       
       if (dfrmd_el(e)) then
          do j = 1,lx,lx-1
             do k = 1,lx,lx-1
                d(1,j,k,e) = d(1,j,k,e) &
                     + G12(1,j,k,e) * dxt(1,1)*dyt(j,j) &
                     + G13(1,j,k,e) * dxt(1,1)*dzt(k,k)
                d(lx,j,k,e) = d(lx,j,k,e) &
                     + G12(lx,j,k,e) * dxt(lx,lx)*dyt(j,j) &
                     + G13(lx,j,k,e) * dxt(lx,lx)*dzt(k,k)
             end do
          end do
          
          do i = 1,lx,lx-1
             do k = 1,lx,lx-1
                d(i,1,k,e) = d(i,1,k,e) &
                     + G12(i,1,k,e) * dyt(1,1)*dxt(i,i) &
                     + G23(i,1,k,e) * dyt(1,1)*dzt(k,k)
                d(i,lx,k,e) = d(i,lx,k,e) &
                     + G12(i,lx,k,e) * dyt(lx,lx)*dxt(i,i) &
                     + G23(i,lx,k,e) * dyt(lx,lx)*dzt(k,k)
             end do
          end do
          do i = 1,lx,lx-1
             do j = 1,lx,lx-1
                d(i,j,1,e) = d(i,j,1,e) &
                     + G13(i,j,1,e) * dzt(1,1)*dxt(i,i) &
                     + G23(i,j,1,e) * dzt(1,1)*dyt(j,j)
                d(i,j,lx,e) = d(i,j,lx,e) &
                     + G13(i,j,lx,e) * dzt(lx,lx)*dxt(i,i) &
                     + G23(i,j,lx,e) * dzt(lx,lx)*dyt(j,j)
             end do
          end do
       end if
    end do
  end subroutine jacobi_update_lx11_dp
  
  subroutine jacobi_update_lx10_dp(d, dxt, dyt, dzt, G11, G22, G33, &
                                G12, G13, G23, dfrmd_el, n)
    integer, parameter :: lx = 10
    integer, intent(in) :: n
    real(kind=dp), intent(inout) :: d(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G11(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G22(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G33(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G12(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G13(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G23(lx, lx, lx, n)
    real(kind=dp), intent(in) :: dxt(lx, lx)
    real(kind=dp), intent(in) :: dyt(lx, lx)
    real(kind=dp), intent(in) :: dzt(lx, lx)
    logical, intent(in) :: dfrmd_el(n)
    integer :: i, j, k, l, e
      
    d = 0d0

    do e = 1,n
       do l = 1,lx
          do k = 1,lx
             do j = 1,lx
                do i = 1,lx
                   d(i,j,k,e) = d(i,j,k,e) + &
                        G11(l,j,k,e) * dxt(i,l)**2
                end do
             end do
          end do
       end do
       do l = 1,lx
          do k = 1,lx
             do j = 1,lx
                do i = 1,lx
                   d(i,j,k,e) = d(i,j,k,e) + &
                        G22(i,l,k,e) * dyt(j,l)**2
                end do
             end do
          end do
       end do
       do l = 1,lx
          do k = 1,lx
             do j = 1,lx
                do i = 1,lx
                   d(i,j,k,e) = d(i,j,k,e) + &
                        G33(i,j,l,e) * dzt(k,l)**2
                end do
             end do
          end do
       end do
       
       if (dfrmd_el(e)) then
          do j = 1,lx,lx-1
             do k = 1,lx,lx-1
                d(1,j,k,e) = d(1,j,k,e) &
                     + G12(1,j,k,e) * dxt(1,1)*dyt(j,j) &
                     + G13(1,j,k,e) * dxt(1,1)*dzt(k,k)
                d(lx,j,k,e) = d(lx,j,k,e) &
                     + G12(lx,j,k,e) * dxt(lx,lx)*dyt(j,j) &
                     + G13(lx,j,k,e) * dxt(lx,lx)*dzt(k,k)
             end do
          end do
          
          do i = 1,lx,lx-1
             do k = 1,lx,lx-1
                d(i,1,k,e) = d(i,1,k,e) &
                     + G12(i,1,k,e) * dyt(1,1)*dxt(i,i) &
                     + G23(i,1,k,e) * dyt(1,1)*dzt(k,k)
                d(i,lx,k,e) = d(i,lx,k,e) &
                     + G12(i,lx,k,e) * dyt(lx,lx)*dxt(i,i) &
                     + G23(i,lx,k,e) * dyt(lx,lx)*dzt(k,k)
             end do
          end do
          do i = 1,lx,lx-1
             do j = 1,lx,lx-1
                d(i,j,1,e) = d(i,j,1,e) &
                     + G13(i,j,1,e) * dzt(1,1)*dxt(i,i) &
                     + G23(i,j,1,e) * dzt(1,1)*dyt(j,j)
                d(i,j,lx,e) = d(i,j,lx,e) &
                     + G13(i,j,lx,e) * dzt(lx,lx)*dxt(i,i) &
                     + G23(i,j,lx,e) * dzt(lx,lx)*dyt(j,j)
             end do
          end do
       end if
    end do
  end subroutine jacobi_update_lx10_dp
  
  subroutine jacobi_update_lx9_dp(d, dxt, dyt, dzt, G11, G22, G33, &
                              G12, G13, G23, dfrmd_el, n)
    integer, parameter :: lx = 9
    integer, intent(in) :: n
    real(kind=dp), intent(inout) :: d(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G11(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G22(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G33(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G12(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G13(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G23(lx, lx, lx, n)
    real(kind=dp), intent(in) :: dxt(lx, lx)
    real(kind=dp), intent(in) :: dyt(lx, lx)
    real(kind=dp), intent(in) :: dzt(lx, lx)
    logical, intent(in) :: dfrmd_el(n)
    integer :: i, j, k, l, e
      
    d = 0d0

    do e = 1,n
       do l = 1,lx
          do k = 1,lx
             do j = 1,lx
                do i = 1,lx
                   d(i,j,k,e) = d(i,j,k,e) + &
                        G11(l,j,k,e) * dxt(i,l)**2
                end do
             end do
          end do
       end do
       do l = 1,lx
          do k = 1,lx
             do j = 1,lx
                do i = 1,lx
                   d(i,j,k,e) = d(i,j,k,e) + &
                        G22(i,l,k,e) * dyt(j,l)**2
                end do
             end do
          end do
       end do
       do l = 1,lx
          do k = 1,lx
             do j = 1,lx
                do i = 1,lx
                   d(i,j,k,e) = d(i,j,k,e) + &
                        G33(i,j,l,e) * dzt(k,l)**2
                end do
             end do
          end do
       end do
       
       if (dfrmd_el(e)) then
          do j = 1,lx,lx-1
             do k = 1,lx,lx-1
                d(1,j,k,e) = d(1,j,k,e) &
                     + G12(1,j,k,e) * dxt(1,1)*dyt(j,j) &
                     + G13(1,j,k,e) * dxt(1,1)*dzt(k,k)
                d(lx,j,k,e) = d(lx,j,k,e) &
                     + G12(lx,j,k,e) * dxt(lx,lx)*dyt(j,j) &
                     + G13(lx,j,k,e) * dxt(lx,lx)*dzt(k,k)
             end do
          end do
          
          do i = 1,lx,lx-1
             do k = 1,lx,lx-1
                d(i,1,k,e) = d(i,1,k,e) &
                     + G12(i,1,k,e) * dyt(1,1)*dxt(i,i) &
                     + G23(i,1,k,e) * dyt(1,1)*dzt(k,k)
                d(i,lx,k,e) = d(i,lx,k,e) &
                     + G12(i,lx,k,e) * dyt(lx,lx)*dxt(i,i) &
                     + G23(i,lx,k,e) * dyt(lx,lx)*dzt(k,k)
             end do
          end do
          do i = 1,lx,lx-1
             do j = 1,lx,lx-1
                d(i,j,1,e) = d(i,j,1,e) &
                     + G13(i,j,1,e) * dzt(1,1)*dxt(i,i) &
                     + G23(i,j,1,e) * dzt(1,1)*dyt(j,j)
                d(i,j,lx,e) = d(i,j,lx,e) &
                     + G13(i,j,lx,e) * dzt(lx,lx)*dxt(i,i) &
                     + G23(i,j,lx,e) * dzt(lx,lx)*dyt(j,j)
             end do
          end do
       end if
    end do
  end subroutine jacobi_update_lx9_dp
  
  subroutine jacobi_update_lx8_dp(d, dxt, dyt, dzt, G11, G22, G33, &
                               G12, G13, G23, dfrmd_el, n)
    integer, parameter :: lx = 8
    integer, intent(in) :: n
    real(kind=dp), intent(inout) :: d(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G11(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G22(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G33(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G12(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G13(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G23(lx, lx, lx, n)
    real(kind=dp), intent(in) :: dxt(lx, lx)
    real(kind=dp), intent(in) :: dyt(lx, lx)
    real(kind=dp), intent(in) :: dzt(lx, lx)
    logical, intent(in) :: dfrmd_el(n)
    integer :: i, j, k, l, e
      
    d = 0d0

    do e = 1,n
       do l = 1,lx
          do k = 1,lx
             do j = 1,lx
                do i = 1,lx
                   d(i,j,k,e) = d(i,j,k,e) + &
                        G11(l,j,k,e) * dxt(i,l)**2
                end do
             end do
          end do
       end do
       do l = 1,lx
          do k = 1,lx
             do j = 1,lx
                do i = 1,lx
                   d(i,j,k,e) = d(i,j,k,e) + &
                        G22(i,l,k,e) * dyt(j,l)**2
                end do
             end do
          end do
       end do
       do l = 1,lx
          do k = 1,lx
             do j = 1,lx
                do i = 1,lx
                   d(i,j,k,e) = d(i,j,k,e) + &
                        G33(i,j,l,e) * dzt(k,l)**2
                end do
             end do
          end do
       end do
       
       if (dfrmd_el(e)) then
          do j = 1,lx,lx-1
             do k = 1,lx,lx-1
                d(1,j,k,e) = d(1,j,k,e) &
                     + G12(1,j,k,e) * dxt(1,1)*dyt(j,j) &
                     + G13(1,j,k,e) * dxt(1,1)*dzt(k,k)
                d(lx,j,k,e) = d(lx,j,k,e) &
                     + G12(lx,j,k,e) * dxt(lx,lx)*dyt(j,j) &
                     + G13(lx,j,k,e) * dxt(lx,lx)*dzt(k,k)
             end do
          end do
          
          do i = 1,lx,lx-1
             do k = 1,lx,lx-1
                d(i,1,k,e) = d(i,1,k,e) &
                     + G12(i,1,k,e) * dyt(1,1)*dxt(i,i) &
                     + G23(i,1,k,e) * dyt(1,1)*dzt(k,k)
                d(i,lx,k,e) = d(i,lx,k,e) &
                     + G12(i,lx,k,e) * dyt(lx,lx)*dxt(i,i) &
                     + G23(i,lx,k,e) * dyt(lx,lx)*dzt(k,k)
             end do
          end do
          do i = 1,lx,lx-1
             do j = 1,lx,lx-1
                d(i,j,1,e) = d(i,j,1,e) &
                     + G13(i,j,1,e) * dzt(1,1)*dxt(i,i) &
                     + G23(i,j,1,e) * dzt(1,1)*dyt(j,j)
                d(i,j,lx,e) = d(i,j,lx,e) &
                     + G13(i,j,lx,e) * dzt(lx,lx)*dxt(i,i) &
                     + G23(i,j,lx,e) * dzt(lx,lx)*dyt(j,j)
             end do
          end do
       end if
    end do
  end subroutine jacobi_update_lx8_dp
  
  subroutine jacobi_update_lx7_dp(d, dxt, dyt, dzt, G11, G22, G33, &
                              G12, G13, G23, dfrmd_el, n)
    integer, parameter :: lx = 7
    integer, intent(in) :: n
    real(kind=dp), intent(inout) :: d(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G11(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G22(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G33(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G12(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G13(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G23(lx, lx, lx, n)
    real(kind=dp), intent(in) :: dxt(lx, lx)
    real(kind=dp), intent(in) :: dyt(lx, lx)
    real(kind=dp), intent(in) :: dzt(lx, lx)
    logical, intent(in) :: dfrmd_el(n)
    integer :: i, j, k, l, e
      
    d = 0d0

    do e = 1,n
       do l = 1,lx
          do k = 1,lx
             do j = 1,lx
                do i = 1,lx
                   d(i,j,k,e) = d(i,j,k,e) + &
                        G11(l,j,k,e) * dxt(i,l)**2
                end do
             end do
          end do
       end do
       do l = 1,lx
          do k = 1,lx
             do j = 1,lx
                do i = 1,lx
                   d(i,j,k,e) = d(i,j,k,e) + &
                        G22(i,l,k,e) * dyt(j,l)**2
                end do
             end do
          end do
       end do
       do l = 1,lx
          do k = 1,lx
             do j = 1,lx
                do i = 1,lx
                   d(i,j,k,e) = d(i,j,k,e) + &
                        G33(i,j,l,e) * dzt(k,l)**2
                end do
             end do
          end do
       end do
       
       if (dfrmd_el(e)) then
          do j = 1,lx,lx-1
             do k = 1,lx,lx-1
                d(1,j,k,e) = d(1,j,k,e) &
                     + G12(1,j,k,e) * dxt(1,1)*dyt(j,j) &
                     + G13(1,j,k,e) * dxt(1,1)*dzt(k,k)
                d(lx,j,k,e) = d(lx,j,k,e) &
                     + G12(lx,j,k,e) * dxt(lx,lx)*dyt(j,j) &
                     + G13(lx,j,k,e) * dxt(lx,lx)*dzt(k,k)
             end do
          end do
          
          do i = 1,lx,lx-1
             do k = 1,lx,lx-1
                d(i,1,k,e) = d(i,1,k,e) &
                     + G12(i,1,k,e) * dyt(1,1)*dxt(i,i) &
                     + G23(i,1,k,e) * dyt(1,1)*dzt(k,k)
                d(i,lx,k,e) = d(i,lx,k,e) &
                     + G12(i,lx,k,e) * dyt(lx,lx)*dxt(i,i) &
                     + G23(i,lx,k,e) * dyt(lx,lx)*dzt(k,k)
             end do
          end do
          do i = 1,lx,lx-1
             do j = 1,lx,lx-1
                d(i,j,1,e) = d(i,j,1,e) &
                     + G13(i,j,1,e) * dzt(1,1)*dxt(i,i) &
                     + G23(i,j,1,e) * dzt(1,1)*dyt(j,j)
                d(i,j,lx,e) = d(i,j,lx,e) &
                     + G13(i,j,lx,e) * dzt(lx,lx)*dxt(i,i) &
                     + G23(i,j,lx,e) * dzt(lx,lx)*dyt(j,j)
             end do
          end do
       end if
    end do
  end subroutine jacobi_update_lx7_dp
  
  subroutine jacobi_update_lx6_dp(d, dxt, dyt, dzt, G11, G22, G33, &
                               G12, G13, G23, dfrmd_el, n)
    integer, parameter :: lx = 6
    integer, intent(in) :: n
    real(kind=dp), intent(inout) :: d(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G11(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G22(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G33(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G12(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G13(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G23(lx, lx, lx, n)
    real(kind=dp), intent(in) :: dxt(lx, lx)
    real(kind=dp), intent(in) :: dyt(lx, lx)
    real(kind=dp), intent(in) :: dzt(lx, lx)
    logical, intent(in) :: dfrmd_el(n)
    integer :: i, j, k, l, e
      
    d = 0d0

    do e = 1,n
       do l = 1,lx
          do k = 1,lx
             do j = 1,lx
                do i = 1,lx
                   d(i,j,k,e) = d(i,j,k,e) + &
                        G11(l,j,k,e) * dxt(i,l)**2
                end do
             end do
          end do
       end do
       do l = 1,lx
          do k = 1,lx
             do j = 1,lx
                do i = 1,lx
                   d(i,j,k,e) = d(i,j,k,e) + &
                        G22(i,l,k,e) * dyt(j,l)**2
                end do
             end do
          end do
       end do
       do l = 1,lx
          do k = 1,lx
             do j = 1,lx
                do i = 1,lx
                   d(i,j,k,e) = d(i,j,k,e) + &
                        G33(i,j,l,e) * dzt(k,l)**2
                end do
             end do
          end do
       end do
       
       if (dfrmd_el(e)) then
          do j = 1,lx,lx-1
             do k = 1,lx,lx-1
                d(1,j,k,e) = d(1,j,k,e) &
                     + G12(1,j,k,e) * dxt(1,1)*dyt(j,j) &
                     + G13(1,j,k,e) * dxt(1,1)*dzt(k,k)
                d(lx,j,k,e) = d(lx,j,k,e) &
                     + G12(lx,j,k,e) * dxt(lx,lx)*dyt(j,j) &
                     + G13(lx,j,k,e) * dxt(lx,lx)*dzt(k,k)
             end do
          end do
          
          do i = 1,lx,lx-1
             do k = 1,lx,lx-1
                d(i,1,k,e) = d(i,1,k,e) &
                     + G12(i,1,k,e) * dyt(1,1)*dxt(i,i) &
                     + G23(i,1,k,e) * dyt(1,1)*dzt(k,k)
                d(i,lx,k,e) = d(i,lx,k,e) &
                     + G12(i,lx,k,e) * dyt(lx,lx)*dxt(i,i) &
                     + G23(i,lx,k,e) * dyt(lx,lx)*dzt(k,k)
             end do
          end do
          do i = 1,lx,lx-1
             do j = 1,lx,lx-1
                d(i,j,1,e) = d(i,j,1,e) &
                     + G13(i,j,1,e) * dzt(1,1)*dxt(i,i) &
                     + G23(i,j,1,e) * dzt(1,1)*dyt(j,j)
                d(i,j,lx,e) = d(i,j,lx,e) &
                     + G13(i,j,lx,e) * dzt(lx,lx)*dxt(i,i) &
                     + G23(i,j,lx,e) * dzt(lx,lx)*dyt(j,j)
             end do
          end do
       end if
    end do
  end subroutine jacobi_update_lx6_dp
  
  subroutine jacobi_update_lx5_dp(d, dxt, dyt, dzt, G11, G22, G33, &
                               G12, G13, G23, dfrmd_el, n)
    integer, parameter :: lx = 5
    integer, intent(in) :: n
    real(kind=dp), intent(inout) :: d(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G11(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G22(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G33(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G12(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G13(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G23(lx, lx, lx, n)
    real(kind=dp), intent(in) :: dxt(lx, lx)
    real(kind=dp), intent(in) :: dyt(lx, lx)
    real(kind=dp), intent(in) :: dzt(lx, lx)
    logical, intent(in) :: dfrmd_el(n)
    integer :: i, j, k, l, e
      
    d = 0d0

    do e = 1,n
       do l = 1,lx
          do k = 1,lx
             do j = 1,lx
                do i = 1,lx
                   d(i,j,k,e) = d(i,j,k,e) + &
                        G11(l,j,k,e) * dxt(i,l)**2
                end do
             end do
          end do
       end do
       do l = 1,lx
          do k = 1,lx
             do j = 1,lx
                do i = 1,lx
                   d(i,j,k,e) = d(i,j,k,e) + &
                        G22(i,l,k,e) * dyt(j,l)**2
                end do
             end do
          end do
       end do
       do l = 1,lx
          do k = 1,lx
             do j = 1,lx
                do i = 1,lx
                   d(i,j,k,e) = d(i,j,k,e) + &
                        G33(i,j,l,e) * dzt(k,l)**2
                end do
             end do
          end do
       end do
       
       if (dfrmd_el(e)) then
          do j = 1,lx,lx-1
             do k = 1,lx,lx-1
                d(1,j,k,e) = d(1,j,k,e) &
                     + G12(1,j,k,e) * dxt(1,1)*dyt(j,j) &
                     + G13(1,j,k,e) * dxt(1,1)*dzt(k,k)
                d(lx,j,k,e) = d(lx,j,k,e) &
                     + G12(lx,j,k,e) * dxt(lx,lx)*dyt(j,j) &
                     + G13(lx,j,k,e) * dxt(lx,lx)*dzt(k,k)
             end do
          end do
          
          do i = 1,lx,lx-1
             do k = 1,lx,lx-1
                d(i,1,k,e) = d(i,1,k,e) &
                     + G12(i,1,k,e) * dyt(1,1)*dxt(i,i) &
                     + G23(i,1,k,e) * dyt(1,1)*dzt(k,k)
                d(i,lx,k,e) = d(i,lx,k,e) &
                     + G12(i,lx,k,e) * dyt(lx,lx)*dxt(i,i) &
                     + G23(i,lx,k,e) * dyt(lx,lx)*dzt(k,k)
             end do
          end do
          do i = 1,lx,lx-1
             do j = 1,lx,lx-1
                d(i,j,1,e) = d(i,j,1,e) &
                     + G13(i,j,1,e) * dzt(1,1)*dxt(i,i) &
                     + G23(i,j,1,e) * dzt(1,1)*dyt(j,j)
                d(i,j,lx,e) = d(i,j,lx,e) &
                     + G13(i,j,lx,e) * dzt(lx,lx)*dxt(i,i) &
                     + G23(i,j,lx,e) * dzt(lx,lx)*dyt(j,j)
             end do
          end do
       end if
    end do
  end subroutine jacobi_update_lx5_dp
  
  subroutine jacobi_update_lx4_dp(d, dxt, dyt, dzt, G11, G22, G33, &
                               G12, G13, G23, dfrmd_el, n)
    integer, parameter :: lx = 4
    integer, intent(in) :: n
    real(kind=dp), intent(inout) :: d(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G11(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G22(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G33(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G12(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G13(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G23(lx, lx, lx, n)
    real(kind=dp), intent(in) :: dxt(lx, lx)
    real(kind=dp), intent(in) :: dyt(lx, lx)
    real(kind=dp), intent(in) :: dzt(lx, lx)
    logical, intent(in) :: dfrmd_el(n)
    integer :: i, j, k, l, e
      
    d = 0d0

    do e = 1,n
       do l = 1,lx
          do k = 1,lx
             do j = 1,lx
                do i = 1,lx
                   d(i,j,k,e) = d(i,j,k,e) + &
                        G11(l,j,k,e) * dxt(i,l)**2
                end do
             end do
          end do
       end do
       do l = 1,lx
          do k = 1,lx
             do j = 1,lx
                do i = 1,lx
                   d(i,j,k,e) = d(i,j,k,e) + &
                        G22(i,l,k,e) * dyt(j,l)**2
                end do
             end do
          end do
       end do
       do l = 1,lx
          do k = 1,lx
             do j = 1,lx
                do i = 1,lx
                   d(i,j,k,e) = d(i,j,k,e) + &
                        G33(i,j,l,e) * dzt(k,l)**2
                end do
             end do
          end do
       end do
       
       if (dfrmd_el(e)) then
          do j = 1,lx,lx-1
             do k = 1,lx,lx-1
                d(1,j,k,e) = d(1,j,k,e) &
                     + G12(1,j,k,e) * dxt(1,1)*dyt(j,j) &
                     + G13(1,j,k,e) * dxt(1,1)*dzt(k,k)
                d(lx,j,k,e) = d(lx,j,k,e) &
                     + G12(lx,j,k,e) * dxt(lx,lx)*dyt(j,j) &
                     + G13(lx,j,k,e) * dxt(lx,lx)*dzt(k,k)
             end do
          end do
          
          do i = 1,lx,lx-1
             do k = 1,lx,lx-1
                d(i,1,k,e) = d(i,1,k,e) &
                     + G12(i,1,k,e) * dyt(1,1)*dxt(i,i) &
                     + G23(i,1,k,e) * dyt(1,1)*dzt(k,k)
                d(i,lx,k,e) = d(i,lx,k,e) &
                     + G12(i,lx,k,e) * dyt(lx,lx)*dxt(i,i) &
                     + G23(i,lx,k,e) * dyt(lx,lx)*dzt(k,k)
             end do
          end do
          do i = 1,lx,lx-1
             do j = 1,lx,lx-1
                d(i,j,1,e) = d(i,j,1,e) &
                     + G13(i,j,1,e) * dzt(1,1)*dxt(i,i) &
                     + G23(i,j,1,e) * dzt(1,1)*dyt(j,j)
                d(i,j,lx,e) = d(i,j,lx,e) &
                     + G13(i,j,lx,e) * dzt(lx,lx)*dxt(i,i) &
                     + G23(i,j,lx,e) * dzt(lx,lx)*dyt(j,j)
             end do
          end do
       end if
    end do
  end subroutine jacobi_update_lx4_dp
  
  subroutine jacobi_update_lx3_dp(d, dxt, dyt, dzt, G11, G22, G33, &
                               G12, G13, G23, dfrmd_el, n)
    integer, parameter :: lx = 3
    integer, intent(in) :: n
    real(kind=dp), intent(inout) :: d(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G11(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G22(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G33(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G12(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G13(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G23(lx, lx, lx, n)
    real(kind=dp), intent(in) :: dxt(lx, lx)
    real(kind=dp), intent(in) :: dyt(lx, lx)
    real(kind=dp), intent(in) :: dzt(lx, lx)
    logical, intent(in) :: dfrmd_el(n)
    integer :: i, j, k, l, e
      
    d = 0d0

    do e = 1,n
       do l = 1,lx
          do k = 1,lx
             do j = 1,lx
                do i = 1,lx
                   d(i,j,k,e) = d(i,j,k,e) + &
                        G11(l,j,k,e) * dxt(i,l)**2
                end do
             end do
          end do
       end do
       do l = 1,lx
          do k = 1,lx
             do j = 1,lx
                do i = 1,lx
                   d(i,j,k,e) = d(i,j,k,e) + &
                        G22(i,l,k,e) * dyt(j,l)**2
                end do
             end do
          end do
       end do
       do l = 1,lx
          do k = 1,lx
             do j = 1,lx
                do i = 1,lx
                   d(i,j,k,e) = d(i,j,k,e) + &
                        G33(i,j,l,e) * dzt(k,l)**2
                end do
             end do
          end do
       end do
       
       if (dfrmd_el(e)) then
          do j = 1,lx,lx-1
             do k = 1,lx,lx-1
                d(1,j,k,e) = d(1,j,k,e) &
                     + G12(1,j,k,e) * dxt(1,1)*dyt(j,j) &
                     + G13(1,j,k,e) * dxt(1,1)*dzt(k,k)
                d(lx,j,k,e) = d(lx,j,k,e) &
                     + G12(lx,j,k,e) * dxt(lx,lx)*dyt(j,j) &
                     + G13(lx,j,k,e) * dxt(lx,lx)*dzt(k,k)
             end do
          end do
          
          do i = 1,lx,lx-1
             do k = 1,lx,lx-1
                d(i,1,k,e) = d(i,1,k,e) &
                     + G12(i,1,k,e) * dyt(1,1)*dxt(i,i) &
                     + G23(i,1,k,e) * dyt(1,1)*dzt(k,k)
                d(i,lx,k,e) = d(i,lx,k,e) &
                     + G12(i,lx,k,e) * dyt(lx,lx)*dxt(i,i) &
                     + G23(i,lx,k,e) * dyt(lx,lx)*dzt(k,k)
             end do
          end do
          do i = 1,lx,lx-1
             do j = 1,lx,lx-1
                d(i,j,1,e) = d(i,j,1,e) &
                     + G13(i,j,1,e) * dzt(1,1)*dxt(i,i) &
                     + G23(i,j,1,e) * dzt(1,1)*dyt(j,j)
                d(i,j,lx,e) = d(i,j,lx,e) &
                     + G13(i,j,lx,e) * dzt(lx,lx)*dxt(i,i) &
                     + G23(i,j,lx,e) * dzt(lx,lx)*dyt(j,j)
             end do
          end do
       end if
    end do
  end subroutine jacobi_update_lx3_dp
  
  subroutine jacobi_update_lx2_dp(d, dxt, dyt, dzt, G11, G22, G33, &
                               G12, G13, G23, dfrmd_el, n)
    integer, parameter :: lx = 2
    integer, intent(in) :: n
    real(kind=dp), intent(inout) :: d(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G11(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G22(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G33(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G12(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G13(lx, lx, lx, n)
    real(kind=dp), intent(in) :: G23(lx, lx, lx, n)
    real(kind=dp), intent(in) :: dxt(lx, lx)
    real(kind=dp), intent(in) :: dyt(lx, lx)
    real(kind=dp), intent(in) :: dzt(lx, lx)
    logical, intent(in) :: dfrmd_el(n)
    integer :: i, j, k, l, e
      
    d = 0d0

    do e = 1,n
       do l = 1,lx
          do k = 1,lx
             do j = 1,lx
                do i = 1,lx
                   d(i,j,k,e) = d(i,j,k,e) + &
                        G11(l,j,k,e) * dxt(i,l)**2
                end do
             end do
          end do
       end do
       do l = 1,lx
          do k = 1,lx
             do j = 1,lx
                do i = 1,lx
                   d(i,j,k,e) = d(i,j,k,e) + &
                        G22(i,l,k,e) * dyt(j,l)**2
                end do
             end do
          end do
       end do
       do l = 1,lx
          do k = 1,lx
             do j = 1,lx
                do i = 1,lx
                   d(i,j,k,e) = d(i,j,k,e) + &
                        G33(i,j,l,e) * dzt(k,l)**2
                end do
             end do
          end do
       end do
       
       if (dfrmd_el(e)) then
          do j = 1,lx,lx-1
             do k = 1,lx,lx-1
                d(1,j,k,e) = d(1,j,k,e) &
                     + G12(1,j,k,e) * dxt(1,1)*dyt(j,j) &
                     + G13(1,j,k,e) * dxt(1,1)*dzt(k,k)
                d(lx,j,k,e) = d(lx,j,k,e) &
                     + G12(lx,j,k,e) * dxt(lx,lx)*dyt(j,j) &
                     + G13(lx,j,k,e) * dxt(lx,lx)*dzt(k,k)
             end do
          end do
          
          do i = 1,lx,lx-1
             do k = 1,lx,lx-1
                d(i,1,k,e) = d(i,1,k,e) &
                     + G12(i,1,k,e) * dyt(1,1)*dxt(i,i) &
                     + G23(i,1,k,e) * dyt(1,1)*dzt(k,k)
                d(i,lx,k,e) = d(i,lx,k,e) &
                     + G12(i,lx,k,e) * dyt(lx,lx)*dxt(i,i) &
                     + G23(i,lx,k,e) * dyt(lx,lx)*dzt(k,k)
             end do
          end do
          do i = 1,lx,lx-1
             do j = 1,lx,lx-1
                d(i,j,1,e) = d(i,j,1,e) &
                     + G13(i,j,1,e) * dzt(1,1)*dxt(i,i) &
                     + G23(i,j,1,e) * dzt(1,1)*dyt(j,j)
                d(i,j,lx,e) = d(i,j,lx,e) &
                     + G13(i,j,lx,e) * dzt(lx,lx)*dxt(i,i) &
                     + G23(i,j,lx,e) * dzt(lx,lx)*dyt(j,j)
             end do
          end do
       end if
    end do
  end subroutine jacobi_update_lx2_dp

end module jacobi_mp

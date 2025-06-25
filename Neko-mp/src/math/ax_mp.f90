! Copyright (c) 2020-2024, The Neko Authors
! All rights reserved.
!
!> Defines mixed-precision Matrix-vector product types
module  ax_product_mp
  use num_types, only : rp, sp, dp
  use coef_mp, only : coef_sp_t, coef_dp_t
  use space_mp, only : space_sp_t, space_dp_t
  use mesh, only : mesh_t
  implicit none
  private

  !> Mixed precision matrix-vector product type
  type, public, abstract :: ax_mp_t
     logical :: is_double = .true.  ! Current precision state
   contains
     procedure(ax_compute_dp), nopass, deferred :: compute_dp
     procedure(ax_compute_sp), nopass, deferred :: compute_sp
     procedure(ax_compute_vector_dp), pass(this), deferred :: compute_vector_dp
     procedure(ax_compute_vector_sp), pass(this), deferred :: compute_vector_sp
     generic  :: compute => compute_dp, compute_sp
     generic  :: compute_vector => compute_vector_dp, compute_vector_sp
  end type ax_mp_t

  ! Abstract interfaces for DP operations
  abstract interface
     subroutine ax_compute_dp(w, u, coef, msh, Xh)
       import
       type(space_dp_t), intent(inout) :: Xh
       type(mesh_t), intent(inout) :: msh
       type(coef_dp_t), intent(inout) :: coef
       real(kind=dp), intent(inout) :: w(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
       real(kind=dp), intent(inout) :: u(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
     end subroutine ax_compute_dp

     subroutine ax_compute_sp(w, u, coef, msh, Xh)
       import
       type(space_sp_t), intent(inout) :: Xh
       type(mesh_t), intent(inout) :: msh
       type(coef_sp_t), intent(inout) :: coef
       real(kind=sp), intent(inout) :: w(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
       real(kind=sp), intent(inout) :: u(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
     end subroutine ax_compute_sp
  end interface

  ! Abstract interfaces for DP operations
  abstract interface
     subroutine ax_compute_vector_dp(this, au, av, aw, u, v, w, coef, msh, Xh)
       import
       class(ax_mp_t), intent(in) :: this
       type(space_dp_t), intent(inout) :: Xh
       type(mesh_t), intent(inout) :: msh
       type(coef_dp_t), intent(inout) :: coef
       real(kind=dp), intent(inout) :: au(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
       real(kind=dp), intent(inout) :: av(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
       real(kind=dp), intent(inout) :: aw(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
       real(kind=dp), intent(inout) :: u(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
       real(kind=dp), intent(inout) :: v(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
       real(kind=dp), intent(inout) :: w(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
     end subroutine ax_compute_vector_dp

     subroutine ax_compute_vector_sp(this, au, av, aw, u, v, w, coef, msh, Xh)
       import
       class(ax_mp_t), intent(in) :: this
       type(space_sp_t), intent(inout) :: Xh
       type(mesh_t), intent(inout) :: msh
       type(coef_sp_t), intent(inout) :: coef
       real(kind=sp), intent(inout) :: au(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
       real(kind=sp), intent(inout) :: av(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
       real(kind=sp), intent(inout) :: aw(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
       real(kind=sp), intent(inout) :: u(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
       real(kind=sp), intent(inout) :: v(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
       real(kind=sp), intent(inout) :: w(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
     end subroutine ax_compute_vector_sp
  end interface  

end module ax_product_mp

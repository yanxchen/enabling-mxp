! Mixed-precision Krylov preconditioner (identity)
module identity_mp
  use math, only : copy_sp, copy_dp
  use precon_mp, only : pc_mp_t
  use num_types, only : sp, dp
  implicit none
  private

  !> Defines a mixed-precision canonical Krylov preconditioner
  type, public, extends(pc_mp_t) :: ident_mp_t
   contains
     procedure, pass(this) :: solve_sp => ident_solve_sp
     procedure, pass(this) :: solve_dp => ident_solve_dp
     procedure, pass(this) :: update => ident_update_mp
  end type ident_mp_t

contains

  !> The naive preconditioner I z = r in single precision
  subroutine ident_solve_sp(this, z, r, n)
    integer, intent(in) :: n
    class(ident_mp_t), intent(inout) :: this
    real(kind=sp), dimension(n), intent(inout) :: z
    real(kind=sp), dimension(n), intent(inout) :: r
    call copy_sp(z, r, n)
  end subroutine ident_solve_sp

  !> The naive preconditioner I z = r in double precision
  subroutine ident_solve_dp(this, z, r, n)
    integer, intent(in) :: n
    class(ident_mp_t), intent(inout) :: this
    real(kind=dp), dimension(n), intent(inout) :: z
    real(kind=dp), dimension(n), intent(inout) :: r
    call copy_dp(z, r, n)
  end subroutine ident_solve_dp

  !> Mandatory update routine
  subroutine ident_update_mp(this)
    class(ident_mp_t), intent(inout) :: this
  end subroutine ident_update_mp

end module identity_mp

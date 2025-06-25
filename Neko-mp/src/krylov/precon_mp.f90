! Mixed-precision Krylov preconditioner
module precon_mp
  use num_types, only : sp, dp
  implicit none
  private

  !> Defines a mixed-precision Krylov preconditioner
  type, public, abstract :: pc_mp_t
   contains
     ! Single precision solve
     procedure(pc_solve_sp), pass(this), deferred :: solve_sp
     ! Double precision solve
     procedure(pc_solve_dp), pass(this), deferred :: solve_dp
     ! Update routine
     procedure(pc_update), pass(this), deferred :: update
  end type pc_mp_t

  !> Abstract interface for solving M z = r in single precision
  abstract interface
     subroutine pc_solve_sp(this, z, r, n)
       import sp
       import :: pc_mp_t
       implicit none
       integer, intent(in) :: n
       class(pc_mp_t), intent(inout) :: this
       real(kind=sp), dimension(n), intent(inout) :: z
       real(kind=sp), dimension(n), intent(inout) :: r
     end subroutine pc_solve_sp

     !> Abstract interface for solving M z = r in double precision
     subroutine pc_solve_dp(this, z, r, n)
       import dp
       import :: pc_mp_t
       implicit none
       integer, intent(in) :: n
       class(pc_mp_t), intent(inout) :: this
       real(kind=dp), dimension(n), intent(inout) :: z
       real(kind=dp), dimension(n), intent(inout) :: r
     end subroutine pc_solve_dp

     subroutine pc_update(this)
       import :: pc_mp_t
       implicit none
       class(pc_mp_t), intent(inout) :: this
     end subroutine pc_update
  end interface

  interface
     !> Create a mixed-precision preconditioner
     module subroutine precon_mp_factory(pc, type_name)
       class(pc_mp_t), allocatable, intent(inout) :: pc
       character(len=*), intent(in) :: type_name
     end subroutine precon_mp_factory

     !> Destroy a mixed-precision preconditioner
     module subroutine precon_mp_destroy(pc)
       class(pc_mp_t), allocatable, intent(inout) :: pc
     end subroutine precon_mp_destroy
  end interface

  public :: precon_mp_factory, precon_mp_destroy
  
end module precon_mp

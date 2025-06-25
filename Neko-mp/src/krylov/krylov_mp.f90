! Mixed-precision Krylov solvers
module krylov_mp
  use gather_scatter_mp, only : gs_sp_t, gs_dp_t
  use ax_product_mp, only : ax_mp_t
  use num_types, only: rp, sp, dp
  use precon_mp, only : pc_mp_t
  use coef_mp, only : coef_dp_t, coef_sp_t
  use mesh, only : mesh_t
  use field_mp, only : field_dp_t, field_sp_t
  use utils, only : neko_error, neko_warning
  use bc_mp, only : bc_list_dp_t, bc_list_sp_t
  use identity_mp, only : ident_mp_t
  ! use device_identity_mp, only : device_ident_mp_t
  use neko_config, only : NEKO_BCKND_DEVICE
  use logger, only : neko_log, LOG_SIZE
  implicit none
  private

  ! Constants for both SP and DP
  integer, public, parameter :: KSP_MAX_ITER = 1e3
  real(kind=dp), public, parameter :: KSP_ABS_TOL = 1e-6  ! absolute tolerance
  real(kind=dp), public, parameter :: KSP_REL_TOL = 1d-9  ! relative tolerance

  !> Monitor type for SP
  type, public :: ksp_monitor_sp_t
     integer :: iter
     real(kind=sp) :: res_start
     real(kind=sp) :: res_final
  end type ksp_monitor_sp_t

  !> Monitor type for DP
  type, public :: ksp_monitor_dp_t
     integer :: iter
     real(kind=dp) :: res_start
     real(kind=dp) :: res_final
  end type ksp_monitor_dp_t

  !> Base abstract type for mixed-precision Krylov solvers
  type, public, abstract :: ksp_mp_t
     class(pc_mp_t), pointer :: M => null()    !< Current preconditioner
     class(pc_mp_t), allocatable :: M_ident
     real(kind=dp) :: rel_tol                  !< relative tolerance
     real(kind=dp) :: abs_tol                  !< absolute tolerance
     integer :: max_iter                       !< Maximum iterations
     integer :: switch_iter = 10               !< Iteration to switch from DP to SP
     logical :: monitor                        !< Monitoring flag
   contains
     procedure, pass(this) :: ksp_init => krylov_init_mp
     procedure, pass(this) :: ksp_free => krylov_free_mp
     procedure, pass(this) :: set_pc => krylov_set_pc_mp
     ! SP/DP solve methods
    !  procedure(ksp_method_sp), pass(this), deferred :: solve_sp
    !  procedure(ksp_method_dp), pass(this), deferred :: solve_dp
     procedure(ksp_method_mp), pass(this), deferred :: solve_mp
     procedure(ksp_method_sp), pass(this), deferred :: solve_sp
     procedure(ksp_method_dp), pass(this), deferred :: solve_dp
    !  ! Coupled solve methods
    !  procedure(ksp_method_coupled_sp), pass(this), deferred :: solve_coupled_sp
    !  procedure(ksp_method_coupled_dp), pass(this), deferred :: solve_coupled_dp
     ! Monitor methods
     procedure, pass(this) :: monitor_start_mp => krylov_monitor_start_mp
     procedure, pass(this) :: monitor_stop_mp => krylov_monitor_stop_mp
     procedure, pass(this) :: monitor_iter_mp => krylov_monitor_iter_mp
     procedure(ksp_mp_t_free), pass(this), deferred :: free
  end type ksp_mp_t

  abstract interface
    !  !> SP solve interface
    !  function ksp_method_sp(this, Ax, x, f, n, coef, blst, gs_h, niter) result(ksp_results)
    !    import :: ksp_mp_t, ax_mp_t, field_mp_t, coef_mp_t
    !    import :: bc_list_mp_t, gs_t, ksp_monitor_sp_t
    !    class(ksp_mp_t), intent(inout) :: this
    !    class(ax_mp_t), intent(inout) :: Ax
    !    type(field_mp_t), intent(inout) :: x
    !    type(field_mp_t), intent(inout) :: f
    !    integer, intent(in) :: n
    !    type(coef_mp_t), intent(inout) :: coef
    !    type(bc_list_mp_t), intent(inout) :: blst
    !    type(gs_t), intent(inout) :: gs_h
    !    integer, intent(in) :: niter
    !    type(ksp_monitor_sp_t) :: ksp_results
    !  end function ksp_method_sp

    !  !> DP solve interface
    !  function ksp_method_dp(this, Ax, x, f, n, coef, blst, gs_h, niter) result(ksp_results)
    !    import :: ksp_mp_t, ax_mp_t, field_mp_t, coef_mp_t
    !    import :: bc_list_mp_t, gs_t, ksp_monitor_dp_t
    !    class(ksp_mp_t), intent(inout) :: this
    !    class(ax_mp_t), intent(inout) :: Ax
    !    type(field_mp_t), intent(inout) :: x
    !    type(field_mp_t), intent(inout) :: f
    !    integer, intent(in) :: n
    !    type(coef_mp_t), intent(inout) :: coef
    !    type(bc_list_mp_t), intent(inout) :: blst
    !    type(gs_t), intent(inout) :: gs_h
    !    integer, intent(in) :: niter
    !    type(ksp_monitor_dp_t) :: ksp_results
    !  end function ksp_method_dp

     !> Mixed-precision solve interface
     function ksp_method_mp(this, Ax, x_dp, x_sp, f_dp, f_sp, n, &
       coef_dp, coef_sp, blst_dp, blst_sp, gs_dp_h) result(ksp_results)
       import :: ksp_mp_t, ax_mp_t, coef_dp_t, coef_sp_t, dp, sp, rp, field_dp_t, field_sp_t
       import :: bc_list_dp_t, bc_list_sp_t, gs_dp_t, ksp_monitor_sp_t, ksp_monitor_dp_t
       class(ksp_mp_t), intent(inout) :: this
       class(ax_mp_t), intent(inout) :: Ax
       type(field_dp_t), intent(inout) :: x_dp
       type(field_sp_t), intent(inout) :: x_sp
       real(kind=dp), dimension(n), intent(inout) :: f_dp
       real(kind=sp), dimension(n), intent(inout) :: f_sp
       integer, intent(in) :: n
       type(coef_dp_t), intent(inout) :: coef_dp
       type(coef_sp_t), intent(inout) :: coef_sp
       type(bc_list_dp_t), intent(inout) :: blst_dp
       type(bc_list_sp_t), intent(inout) :: blst_sp
       type(gs_dp_t), intent(inout) :: gs_dp_h
      !  integer, intent(in) :: niter
      !  type(ksp_monitor_sp_t) :: ksp_results_sp
       type(ksp_monitor_dp_t) :: ksp_results
     end function ksp_method_mp
  end interface

  abstract interface
     function ksp_method_sp(this, Ax, x, f, n, coef, blst, gs_h) result(ksp_results)
       import :: ksp_mp_t, ax_mp_t, coef_dp_t, coef_sp_t, dp, sp, rp, field_dp_t, field_sp_t
       import :: bc_list_dp_t, bc_list_sp_t, gs_dp_t, gs_sp_t, ksp_monitor_sp_t, ksp_monitor_dp_t
       class(ksp_mp_t), intent(inout) :: this
       class(ax_mp_t), intent(inout) :: Ax
       type(field_sp_t), intent(inout) :: x
       real(kind=sp), dimension(n), intent(inout) :: f
       integer, intent(in) :: n
       type(coef_sp_t), intent(inout) :: coef
       type(bc_list_sp_t), intent(inout) :: blst
       type(gs_sp_t), intent(inout) :: gs_h
      !  integer, intent(in) :: niter
      !  type(ksp_monitor_sp_t) :: ksp_results_sp
       type(ksp_monitor_dp_t) :: ksp_results
     end function ksp_method_sp
  end interface

  abstract interface
     function ksp_method_dp(this, Ax, x, f, n, coef, blst, gs_h) result(ksp_results)
       import :: ksp_mp_t, ax_mp_t, coef_dp_t, coef_sp_t, dp, sp, rp, field_dp_t, field_sp_t
       import :: bc_list_dp_t, bc_list_sp_t, gs_dp_t, gs_sp_t, ksp_monitor_sp_t, ksp_monitor_dp_t
       class(ksp_mp_t), intent(inout) :: this
       class(ax_mp_t), intent(inout) :: Ax
       type(field_dp_t), intent(inout) :: x
       real(kind=dp), dimension(n), intent(inout) :: f
       integer, intent(in) :: n
       type(coef_dp_t), intent(inout) :: coef
       type(bc_list_dp_t), intent(inout) :: blst
       type(gs_dp_t), intent(inout) :: gs_h
      !  integer, intent(in) :: niter
      !  type(ksp_monitor_sp_t) :: ksp_results_sp
       type(ksp_monitor_dp_t) :: ksp_results
     end function ksp_method_dp
  end interface

  abstract interface
     !> Free interface
     subroutine ksp_mp_t_free(this)
       import :: ksp_mp_t
       class(ksp_mp_t), intent(inout) :: this
     end subroutine ksp_mp_t_free
  end interface

contains

  !> Initialize mixed-precision Krylov solver
  subroutine krylov_init_mp(this, max_iter, rel_tol, abs_tol, M, monitor)
    class(ksp_mp_t), target, intent(inout) :: this
    integer, intent(in), optional :: max_iter
    real(kind=dp), intent(in), optional :: rel_tol, abs_tol
    class(pc_mp_t), optional, target, intent(in) :: M
    logical, intent(in), optional :: monitor

    call krylov_free_mp(this)

    this%max_iter = KSP_MAX_ITER
    if (present(max_iter)) this%max_iter = max_iter

    this%rel_tol = KSP_REL_TOL
    if (present(rel_tol)) this%rel_tol = rel_tol

    this%abs_tol = KSP_ABS_TOL
    if (present(abs_tol)) this%abs_tol = abs_tol

    this%monitor = .false.
    if (present(monitor)) this%monitor = monitor

    if (present(M)) then 
      this%M => M
    else
      if (.not. associated(this%M)) then
        if (NEKO_BCKND_DEVICE .eq. 1) then
           !
        else
           allocate(ident_mp_t::this%M_ident)
        end if
        this%M => this%M_ident
     end if
    end if

  end subroutine krylov_init_mp

  !> Free mixed-precision Krylov solver
  subroutine krylov_free_mp(this)
    class(ksp_mp_t), intent(inout) :: this

    ! if (allocated(this%M)) deallocate(this%M)

  end subroutine krylov_free_mp

  !> Set preconditioner for mixed-precision solver
  subroutine krylov_set_pc_mp(this, M)
    class(ksp_mp_t), intent(inout) :: this
    class(pc_mp_t), target, intent(in) :: M
    if (associated(this%M)) then
      select type (pc => this%M)
      type is (ident_mp_t)
      class default
         call neko_error('Preconditioner already defined')
      end select
    end if

    this%M => M
  end subroutine krylov_set_pc_mp

  !> Start monitoring
  subroutine krylov_monitor_start_mp(this, name)
    class(ksp_mp_t), intent(inout) :: this
    character(len=*) :: name
    character(len=LOG_SIZE) :: log_buf

    if (this%monitor) then
       write(log_buf, '(A)') 'Krylov monitor (' // trim(name) // ')'
       call neko_log%section(trim(log_buf))
       call neko_log%newline()
       call neko_log%begin()
       write(log_buf, '(A)') ' Iter.       Residual'
       call neko_log%message(log_buf)
       write(log_buf, '(A)') '---------------------'
       call neko_log%message(log_buf)
    end if

  end subroutine krylov_monitor_start_mp

  !> Stop monitoring
  subroutine krylov_monitor_stop_mp(this)
    class(ksp_mp_t), intent(inout) :: this

    if (this%monitor) then
       call neko_log%end()
       call neko_log%end_section()
       call neko_log%newline()
    end if

  end subroutine krylov_monitor_stop_mp

  !> Monitor iteration
  subroutine krylov_monitor_iter_mp(this, iter, rnorm)
    class(ksp_mp_t), intent(in) :: this
    integer, intent(in) :: iter
    class(*), intent(in) :: rnorm
    real(kind=rp) :: rnorm_print
    character(len=LOG_SIZE) :: log_buf

    select type (rnorm)
    type is (real(kind=sp))
       rnorm_print = real(rnorm, kind=rp)
    type is (real(kind=dp))
       rnorm_print = rnorm
    end select

    if (this%monitor) then
       write(log_buf, '(I6,E15.7)') iter, rnorm_print
       call neko_log%message(log_buf)
    end if
    
  end subroutine krylov_monitor_iter_mp

end module krylov_mp

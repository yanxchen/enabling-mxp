!> Mixed-precision conjugate gradient solver
module cg_mp
   use num_types, only : rp, sp, dp, xp
   use precon_mp, only : pc_mp_t
   use ax_product_mp, only : ax_mp_t
   use field_mp, only : field_sp_t, field_dp_t
   use coef_mp, only : coef_sp_t, coef_dp_t
   use gather_scatter_mp, only : gs_sp_t, gs_dp_t, GS_OP_ADD
   use bc_mp, only : bc_list_dp_t, bc_list_sp_t, bc_list_apply_dp, bc_list_apply_sp
   use math, only : glsc3_sp, glsc3_dp, rzero_sp, rzero_dp, copy_sp, copy_dp, abscmp, &
                    copy, rzero, glsc3_sp_eft
   use krylov_mp, only : ksp_mp_t, ksp_monitor_sp_t, ksp_monitor_dp_t
   use comm
   implicit none
   private
 
   integer, parameter :: CG_P_SPACE = 7
 
   !> Mixed-precision conjugate gradient solver
   type, public, extends(ksp_mp_t) :: cg_mp_t
      private
      ! SP arrays
      real(kind=sp), allocatable :: w_sp(:)
      real(kind=sp), allocatable :: r_sp(:)
      real(kind=sp), allocatable :: p_sp(:,:)
      real(kind=sp), allocatable :: z_sp(:)
      real(kind=sp), allocatable :: alpha_sp(:)
      ! DP arrays
      real(kind=dp), allocatable :: w_dp(:)
      real(kind=dp), allocatable :: r_dp(:)
      real(kind=dp), allocatable :: p_dp(:,:)
      real(kind=dp), allocatable :: z_dp(:)
      real(kind=dp), allocatable :: alpha_dp(:)
    contains
      procedure, pass(this) :: init => cg_init_mp
      procedure, pass(this) :: free => cg_free_mp
      procedure, pass(this) :: to_dp => cg_to_dp
      procedure, pass(this) :: to_sp => cg_to_sp
      procedure, pass(this) :: solve_mp => cg_solve_mp
      procedure, pass(this) :: solve_sp => cg_solve_sp
      procedure, pass(this) :: solve_dp => cg_solve_dp
   end type cg_mp_t
 
 contains
 
   subroutine cg_init_mp(this, n, max_iter, switch_iter, M, rel_tol, abs_tol, monitor)
     class(cg_mp_t), intent(inout), target :: this
     integer, optional, intent(in) :: max_iter
     integer, optional, intent(in) :: switch_iter
     class(pc_mp_t), optional, intent(inout), target :: M
     integer, intent(in) :: n
     real(kind=dp), optional, intent(inout) :: rel_tol
     real(kind=dp), optional, intent(inout) :: abs_tol
     logical, optional, intent(in) :: monitor
 
     call this%free()
 
     allocate(this%w_sp(n))
     allocate(this%r_sp(n))
     allocate(this%p_sp(n,CG_P_SPACE))
     allocate(this%z_sp(n))
     allocate(this%alpha_sp(CG_P_SPACE))
 
     allocate(this%w_dp(n))
     allocate(this%r_dp(n))
     allocate(this%p_dp(n,CG_P_SPACE))
     allocate(this%z_dp(n))
     allocate(this%alpha_dp(CG_P_SPACE))
 
     if (present(M)) then
        this%M => M
     end if
     if (present(rel_tol) .and. present(abs_tol) .and. present(monitor)) then
        call this%ksp_init(max_iter, rel_tol, abs_tol, monitor = monitor)
     else if (present(rel_tol) .and. present(abs_tol)) then
        call this%ksp_init(max_iter, rel_tol, abs_tol)
     else if (present(monitor) .and. present(abs_tol)) then
        call this%ksp_init(max_iter, abs_tol = abs_tol, monitor = monitor)
     else if (present(rel_tol) .and. present(monitor)) then
        call this%ksp_init(max_iter, rel_tol, monitor = monitor)
     else if (present(rel_tol)) then
        call this%ksp_init(max_iter, rel_tol = rel_tol)
     else if (present(abs_tol)) then
        call this%ksp_init(max_iter, abs_tol = abs_tol)
     else if (present(monitor)) then
        call this%ksp_init(max_iter, monitor = monitor)
     else
        call this%ksp_init(max_iter)
     end if
 
   end subroutine cg_init_mp
 
   !> Deallocate a standard PCG solver
   subroutine cg_free_mp(this)
     class(cg_mp_t), intent(inout) :: this
 
     call this%ksp_free()
 
     if (allocated(this%w_sp)) then
        deallocate(this%w_sp)
     end if
     if (allocated(this%r_sp)) then
        deallocate(this%r_sp)
     end if
     if (allocated(this%p_sp)) then
        deallocate(this%p_sp)
     end if
     if (allocated(this%z_sp)) then
        deallocate(this%z_sp)
     end if
     if (allocated(this%alpha_sp)) then
        deallocate(this%alpha_sp)
     end if
     if (allocated(this%w_dp)) then
        deallocate(this%w_dp)
     end if
     if (allocated(this%r_dp)) then
        deallocate(this%r_dp)
     end if
     if (allocated(this%p_dp)) then
        deallocate(this%p_dp)
     end if 
     if (allocated(this%z_dp)) then
        deallocate(this%z_dp)
     end if
     if (allocated(this%alpha_dp)) then
        deallocate(this%alpha_dp)
     end if

     nullify(this%M)

   end subroutine cg_free_mp
 
   subroutine cg_to_dp(this)
     class(cg_mp_t), intent(inout) :: this
     
     this%w_dp = real(this%w_sp, kind=dp)
     this%r_dp = real(this%r_sp, kind=dp)
     this%p_dp = real(this%p_sp, kind=dp)
     this%z_dp = real(this%z_sp, kind=dp)
     this%alpha_dp = real(this%alpha_sp, kind=dp)
     
   end subroutine cg_to_dp
 
   subroutine cg_to_sp(this)
     class(cg_mp_t), intent(inout) :: this
     
     this%w_sp = real(this%w_dp, kind=sp)
     this%r_sp = real(this%r_dp, kind=sp)
     this%p_sp = real(this%p_dp, kind=sp)
     this%z_sp = real(this%z_dp, kind=sp)
     this%alpha_sp = real(this%alpha_dp, kind=sp)
     
   end subroutine cg_to_sp
 
   !   !> Standard PCG solve
   ! function cg_solve_mp(this, Ax, x, f, n, coef_dp, coef_sp, blst_dp, blst_sp, gs_dp_h, gs_sp_h) result(ksp_results)
   !   class(cg_mp_t), intent(inout) :: this
   !   class(ax_mp_t), intent(inout) :: Ax
   !   type(field_dp_t), intent(inout) :: x
   !   integer, intent(in) :: n
   !   real(kind=dp), dimension(n), intent(inout) :: f
   !   type(coef_dp_t), intent(inout) :: coef_dp
   !   type(coef_sp_t), intent(inout) :: coef_sp
   !   type(bc_list_dp_t), intent(inout) :: blst_dp
   !   type(bc_list_sp_t), intent(inout) :: blst_sp
   !   type(gs_dp_t), intent(inout) :: gs_dp_h
   !   type(gs_sp_t), intent(inout) :: gs_sp_h
   !   type(ksp_monitor_dp_t) :: ksp_results
   !   integer :: iter, max_iter, i, j, k, p_cur, p_prev
   !   real(kind=dp) :: rnorm, rtr, rtz2, rtz1, x_plus(NEKO_BLK_SIZE)
   !   real(kind=dp) :: beta, pap, norm_fac
   !  !  type(field_dp_t), pointer :: x_dp => null()
   !   ! type(coef_dp_t), pointer :: coef_dp => null()
 
   !   max_iter = this%max_iter
   !   ! coef_dp => coef%get_dp()
   !   norm_fac = 1.0_dp / sqrt(coef_dp%volume)
   !  !  x_dp => x%get_dp()
 
   !   associate(w => this%w_dp, r => this%r_dp, p => this%p_dp, &
   !        z => this%z_dp, alpha => this%alpha_dp, coef => coef_dp)
 
   !     rtz1 = 1.0_dp
   !     call rzero(x%x, n)
   !     call rzero(p(1,CG_P_SPACE), n)
   !     call copy(r, f, n)
 
   !     rtr = glsc3_dp(r, coef_dp%mult, r, n)
   !     rnorm = sqrt(rtr) * norm_fac
   !     ksp_results%res_start = rnorm
   !     ksp_results%res_final = rnorm
   !     ksp_results%iter = 0
   !     p_prev = CG_P_SPACE
   !     p_cur = 1
   !     if(abscmp(rnorm, 0.0_dp)) return
   !     call this%monitor_start_mp('CG (DP)')
   !     do iter = 1, max_iter
   !        call this%M%solve_dp(z, r, n)
   !        rtz2 = rtz1
   !        rtz1 = glsc3_dp(r, coef%mult, z, n)
 
   !        beta = rtz1 / rtz2
   !        if (iter .eq. 1) beta = 0.0_dp
   !        do i = 1, n
   !           p(i,p_cur) = z(i) + beta * p(i,p_prev)
   !        end do
 
   !        call Ax%compute_dp(w, p(1,p_cur), coef, x%msh, x%Xh)
   !        call gs_dp_h%op(w, n, GS_OP_ADD)
   !        call bc_list_apply_dp(blst_dp, w, n)
 
   !        pap = glsc3_dp(w, coef%mult, p(1,p_cur), n)
 
   !        alpha(p_cur) = rtz1 / pap
   !        call second_cg_part_mp(rtr, r, coef%mult, w, alpha(p_cur), n)
   !        rnorm = sqrt(rtr) * norm_fac
   !        call this%monitor_iter_mp(iter, rnorm)
 
   !        if ((p_cur .eq. CG_P_SPACE) .or. &
   !            (rnorm .lt. this%abs_tol) .or. iter .eq. max_iter) then
   !           do i = 0, n, NEKO_BLK_SIZE
   !              if (i + NEKO_BLK_SIZE .le. n) then
   !                 do k = 1, NEKO_BLK_SIZE
   !                    x_plus(k) = 0.0_rp
   !                 end do
   !                 do j = 1, p_cur
   !                    do k = 1, NEKO_BLK_SIZE
   !                       x_plus(k) = x_plus(k) + alpha(j) * p(i+k,j)
   !                    end do
   !                 end do
   !                 do k = 1, NEKO_BLK_SIZE
   !                    x%x(i+k,1,1,1) = x%x(i+k,1,1,1) + x_plus(k)
   !                 end do
   !              else
   !                 do k = 1, n-i
   !                    x_plus(1) = 0.0_rp
   !                    do j = 1, p_cur
   !                       x_plus(1) = x_plus(1) + alpha(j) * p(i+k,j)
   !                    end do
   !                    x%x(i+k,1,1,1) = x%x(i+k,1,1,1) + x_plus(1)
   !                 end do
   !              end if
   !           end do
   !           p_prev = p_cur
   !           p_cur = 1
   !           if (rnorm .lt. this%abs_tol) exit
   !        else
   !           p_prev = p_cur
   !           p_cur = p_cur + 1
   !        end if
   !     end do
   !   end associate
   !   call this%monitor_stop_mp()
   !   ksp_results%res_final = rnorm
   !   ksp_results%iter = iter
 
   ! end function cg_solve_mp
 
   function cg_solve_mp(this, Ax, x_dp, x_sp, f_dp, f_sp, n, &
     coef_dp, coef_sp, blst_dp, blst_sp, gs_dp_h) result(ksp_results)
     class(cg_mp_t), intent(inout) :: this
     class(ax_mp_t), intent(inout) :: Ax
     type(field_dp_t), intent(inout) :: x_dp
     type(field_sp_t), intent(inout) :: x_sp
     integer, intent(in) :: n
     real(kind=dp), dimension(n), intent(inout) :: f_dp
     real(kind=sp), dimension(n), intent(inout) :: f_sp
     type(coef_dp_t), intent(inout) :: coef_dp
     type(coef_sp_t), intent(inout) :: coef_sp
     type(bc_list_dp_t), intent(inout) :: blst_dp
     type(bc_list_sp_t), intent(inout) :: blst_sp
     type(gs_dp_t), intent(inout) :: gs_dp_h
     type(ksp_monitor_dp_t) :: ksp_results
     integer :: iter, max_iter, i, j, k, p_cur, p_prev
     
     real(kind=sp) :: rnorm, rtr, rtz2, rtz1, x_plus(NEKO_BLK_SIZE)
     real(kind=sp) :: beta, pap, norm_fac

     real(kind=sp) :: tol

     tol = this%abs_tol
 
     max_iter = this%max_iter
     norm_fac = 1.0_sp / sqrt(coef_dp%volume)
 
     associate(w => this%w_sp, r => this%r_sp, p => this%p_sp, &
          z => this%z_sp, alpha => this%alpha_sp, coef => coef_sp, &
          f => f_sp, x => x_sp)
 
       rtz1 = 1.0_sp
       call rzero_sp(x%x, n)
       call rzero_sp(p(1,CG_P_SPACE), n)
       call copy_sp(r, f, n)
 
       rtr = glsc3_sp(r, coef%mult, r, n)
       rnorm = sqrt(rtr) * norm_fac
       ksp_results%res_start = rnorm
       ksp_results%res_final = rnorm
       ksp_results%iter = 0
       p_prev = CG_P_SPACE
       p_cur = 1
       if(abscmp(rnorm, 0.0_sp)) return
       call this%monitor_start_mp('CG (SP)')
       do iter = 1, max_iter
          call this%M%solve_sp(z, r, n)
          rtz2 = rtz1
          rtz1 = glsc3_sp(r, coef%mult, z, n)
         !  rtz1 = glsc3_sp_eft(r, coef%mult, z, n)

         !  this%r_dp = real(r, kind=dp)
         !  this%z_dp = real(z, kind=dp)
         !  rtz1 = glsc3_dp(this%r_dp, coef_dp%mult, this%z_dp, n)
 
          beta = rtz1 / rtz2
          if (iter .eq. 1) beta = 0.0_sp
          do i = 1, n
             p(i,p_cur) = z(i) + beta * p(i,p_prev)
          end do
 
          call Ax%compute_sp(w, p(1,p_cur), coef, x%msh, x%Xh)
          
          this%w_dp = real(w, kind=dp)
         !  call gs_sp_h%op(w, n, GS_OP_ADD)
          call gs_dp_h%op(this%w_dp, n, GS_OP_ADD)
         !  this%w_sp = this%w_dp

         !  this%w_dp = real(w, kind=dp)
          call bc_list_apply_dp(blst_dp, this%w_dp, n)
          this%w_sp = this%w_dp
 
          pap = glsc3_sp(w, coef%mult, p(1,p_cur), n)
         !  pap = glsc3_sp_eft(w, coef%mult, p(1,p_cur), n)

         !  this%w_dp = real(w, kind=dp)
         !  this%p_dp = real(p, kind=dp)
         !  pap = glsc3_dp(this%w_dp, coef_dp%mult, this%p_dp(1,p_cur), n)
 
          alpha(p_cur) = rtz1 / pap
          call second_cg_part_sp(rtr, r, coef%mult, w, alpha(p_cur), n)
          rnorm = sqrt(rtr) * norm_fac
          call this%monitor_iter_mp(iter, rnorm)
 
          if ((p_cur .eq. CG_P_SPACE) .or. &
              (rnorm .lt. tol) .or. iter .eq. max_iter) then
             do i = 0, n, NEKO_BLK_SIZE
                if (i + NEKO_BLK_SIZE .le. n) then
                   do k = 1, NEKO_BLK_SIZE
                      x_plus(k) = 0.0_sp
                   end do
                   do j = 1, p_cur
                      do k = 1, NEKO_BLK_SIZE
                         x_plus(k) = x_plus(k) + alpha(j) * p(i+k,j)
                      end do
                   end do
                   do k = 1, NEKO_BLK_SIZE
                      x%x(i+k,1,1,1) = x%x(i+k,1,1,1) + x_plus(k)
                   end do
                else
                   do k = 1, n-i
                      x_plus(1) = 0.0_rp
                      do j = 1, p_cur
                         x_plus(1) = x_plus(1) + alpha(j) * p(i+k,j)
                      end do
                      x%x(i+k,1,1,1) = x%x(i+k,1,1,1) + x_plus(1)
                   end do
                end if
             end do
             p_prev = p_cur
             p_cur = 1
             if (rnorm .lt. tol) exit
          else
             p_prev = p_cur
             p_cur = p_cur + 1
          end if
       end do
     end associate
     call this%monitor_stop_mp()
     ksp_results%res_final = rnorm
     ksp_results%iter = iter
 
   end function cg_solve_mp
 
   function cg_solve_sp(this, Ax, x, f, n, coef, blst, gs_h) result(ksp_results)
      class(cg_mp_t), intent(inout) :: this
      class(ax_mp_t), intent(inout) :: Ax
      type(field_sp_t), intent(inout) :: x
      real(kind=sp), dimension(n), intent(inout) :: f
      integer, intent(in) :: n
      type(coef_sp_t), intent(inout) :: coef
      type(bc_list_sp_t), intent(inout) :: blst
      type(gs_sp_t), intent(inout) :: gs_h
      type(ksp_monitor_dp_t) :: ksp_results
   end function cg_solve_sp

   function cg_solve_dp(this, Ax, x, f, n, coef, blst, gs_h) result(ksp_results)
      class(cg_mp_t), intent(inout) :: this
      class(ax_mp_t), intent(inout) :: Ax
      type(field_dp_t), intent(inout) :: x
      real(kind=dp), dimension(n), intent(inout) :: f
      integer, intent(in) :: n
      type(coef_dp_t), intent(inout) :: coef
      type(bc_list_dp_t), intent(inout) :: blst
      type(gs_dp_t), intent(inout) :: gs_h
      type(ksp_monitor_dp_t) :: ksp_results
   end function cg_solve_dp

   subroutine second_cg_part_mp(rtr, r, mult, w, alpha, n)
     integer, intent(in) :: n
     real(kind=rp), intent(inout) :: r(n), rtr
     real(kind=xp) :: tmp
     real(kind=rp), intent(in) ::mult(n), w(n), alpha
     integer :: i, ierr
 
     tmp = 0.0_xp
     do i = 1, n
        r(i) = r(i) - alpha*w(i)
        tmp = tmp + r(i) * r(i) * mult(i)
     end do
     call MPI_Allreduce(MPI_IN_PLACE, tmp, 1, &
          MPI_EXTRA_PRECISION, MPI_SUM, NEKO_COMM, ierr)
     rtr = tmp
 
   end subroutine second_cg_part_mp
 
   subroutine second_cg_part_sp(rtr, r, mult, w, alpha, n)
     integer, intent(in) :: n
     real(kind=sp), intent(inout) :: r(n), rtr
     real(kind=xp) :: tmp
     real(kind=sp), intent(in) ::mult(n), w(n), alpha
     integer :: i, ierr
 
     tmp = 0.0_xp
     do i = 1, n
        r(i) = r(i) - alpha*w(i)
        tmp = tmp + r(i) * r(i) * mult(i)
     end do
     call MPI_Allreduce(MPI_IN_PLACE, tmp, 1, &
          MPI_EXTRA_PRECISION, MPI_SUM, NEKO_COMM, ierr)
     rtr = tmp
 
   end subroutine second_cg_part_sp
 
 end module cg_mp
 

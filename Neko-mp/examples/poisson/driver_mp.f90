program poisson
   use neko
   use ax_poisson_mp
   use gather_scatter_mp, only : gs_sp_t, gs_dp_t
   use dofmap_mp, only : dofmap_sp_t, dofmap_dp_t
   use krylov_mp
   use cg_mp
   use dirichlet_mp
   use ax_poisson_mp
   use field_mp, only : field_dp_t, field_sp_t
   use bc_mp
   use coef_mp, only : coef_dp_t, coef_sp_t
   use space_mp, only : space_dp_t, space_sp_t
   use space
   use math
   use jacobi_mp
 
   implicit none
 
   character(len=NEKO_FNAME_LEN) :: fname, lxchar, iterchar
   type(mesh_t) :: msh
   type(file_t) :: nmsh_file, mf
   type(space_dp_t) :: Xh_dp
   type(space_sp_t) :: Xh_sp
   type(dofmap_dp_t) :: dm_dp
   type(dofmap_sp_t) :: dm_sp
   type(gs_dp_t) :: gs_dp_h
   type(gs_sp_t) :: gs_sp_h
   type(dirichlet_dp_t) :: dir_bc_dp
   type(dirichlet_sp_t) :: dir_bc_sp
   type(bc_list_dp_t) :: bclst_dp
   type(bc_list_sp_t) :: bclst_sp
   type(field_dp_t) :: x
   type(field_sp_t) :: x_sp
   type(ax_poisson_mp_t) :: ax
   type(coef_dp_t) :: coef_dp
   type(coef_sp_t) :: coef_sp
   type(cg_mp_t) :: solver
   type(jacobi_mp_t) :: jacobi_precond
   type(ksp_monitor_sp_t) :: ksp_mon_sp
   type(ksp_monitor_dp_t) :: ksp_mon_dp
   integer :: argc, lx, n, n_glb, niter, ierr
   character(len=80) :: suffix
   real(kind=dp), allocatable :: f(:)
   real(kind=sp), allocatable :: f_sp(:)
   real(kind=dp) :: tol
   tol = 1.e-9
   ! tol = -1.0
 
   argc = command_argument_count()
 
   if ((argc .lt. 3) .or. (argc .gt. 3)) then
      if (pe_rank .eq. 0) then
         write(*,*) 'Usage: ./poisson <neko mesh> <N> <Iterations>'
      end if
      stop
   end if
   
   call neko_init 
   call get_command_argument(1, fname)
   call get_command_argument(2, lxchar)
   call get_command_argument(3, iterchar)
   read(lxchar, *) lx
   read(iterchar, *) niter
   
   nmsh_file = file_t(fname)
   call nmsh_file%read(msh)  
 
   call Xh_dp%init(1, lx, lx, lx)
   call Xh_sp%init(1, lx, lx, lx)
 
   call dm_dp%init(msh, Xh_dp)
   call dm_sp%init(msh, Xh_sp)
 
   call gs_dp_h%init(dm_dp)
 
   call coef_dp%init(gs_dp_h)
   ! call coef_sp%init(gs_sp_h)

   coef_sp%G11 = real(coef_dp%G11, kind=sp)
   coef_sp%G22 = real(coef_dp%G22, kind=sp)
   coef_sp%G33 = real(coef_dp%G33, kind=sp)
   coef_sp%mult = real(coef_dp%mult, kind=sp)
   
   call x%init(dm_dp, "x-DP")
   call x_sp%init(dm_sp, "x-SP")
 
   n = Xh_dp%lx * Xh_dp%ly * Xh_dp%lz * msh%nelv
 
   call dir_bc_dp%init_base(coef_dp)
   call dir_bc_dp%set_g(real(0.0d0,dp))
 
   call dir_bc_sp%init_base(coef_sp)
   call dir_bc_sp%set_g(real(0.0d0,sp))
  
   !user specified
   call set_bc_mp(dir_bc_dp, msh)
  
   call dir_bc_dp%finalize()
   ! call dir_bc_sp%finalize()
 
   call bc_list_init_dp(bclst_dp)
   call bc_list_add_dp(bclst_dp,dir_bc_dp)
 
   call bc_list_init_sp(bclst_sp)
   call bc_list_add_sp(bclst_sp,dir_bc_sp)

   ! call jacobi_precond%init(coef_dp, dm_dp, gs_dp_h)
 
   ! call solver%init(n, niter, abs_tol = tol, M = jacobi_precond, monitor=.true.)

   call solver%init(n, niter, abs_tol = tol, monitor=.true.)
 
   allocate(f(n))
   allocate(f_sp(n))
 
   !user specified
   call rzero(f,n)
   call rzero_sp(f_sp,n)

   call set_f_mp(f, coef_dp%mult, dm_dp, n, gs_dp_h)

   call bc_list_apply_dp(bclst_dp,f,n)
   ! call bc_list_apply_sp(bclst_sp,f_sp,n)
   f_sp = real(f, kind=sp)
   
   ksp_mon_dp = solver%solve_mp(ax, x, x_sp, f, f_sp, n, coef_dp, coef_sp, bclst_dp, bclst_sp, gs_dp_h)
   n_glb = Xh_dp%lx * Xh_dp%ly * Xh_dp%lz * msh%glb_nelv
   
   call MPI_Barrier(NEKO_COMM, ierr)
 
   call set_timer_flop_cnt(0, msh%glb_nelv, x%Xh%lx, niter, n_glb, ksp_mon_dp)
   ksp_mon_dp = solver%solve_mp(ax, x, x_sp, f, f_sp, n, coef_dp, coef_sp, bclst_dp, bclst_sp, gs_dp_h)
   call set_timer_flop_cnt(1, msh%glb_nelv, x%Xh%lx, niter, n_glb, ksp_mon_dp)
   
   fname = 'out.fld'
   mf =  file_t(fname)
 !   call mf%write(x)
   deallocate(f)
   call solver%free()
   
   call dir_bc_dp%free()
   ! call dir_bc_sp%free()
 
   call bc_list_free_dp(bclst_dp)
   ! call bc_list_free_sp(bclst_sp)
 
   call Xh_dp%free()
   call Xh_sp%free()
   
   call x%free()
   call msh%free() 
   call neko_finalize
 
 end program poisson
 
 subroutine set_timer_flop_cnt(iset, nelt, nx1, niter, n, ksp_mon)
   use comm
   use krylov
   use num_types
   use krylov_mp
   implicit none
 
   integer :: iset
   integer, intent(inout) :: nelt
   integer, intent(inout) :: nx1
   integer, intent(inout) :: niter
   integer, intent(inout) :: n
   type(ksp_monitor_dp_t), intent(in) :: ksp_mon
   real(kind=dp), save :: time0, time1, mflops, flop_a, flop_cg
   real(kind=dp) :: nxyz, nx
 
 
   nx = dble(nx1)
   nxyz = dble(nx1 * nx1 * nx1)
   if (iset .eq. 0) then
      time0 = MPI_Wtime()
   else
      flop_a = (15d0 * nxyz + 12d0 * nx * nxyz) * dble(nelt) * dble(niter)
      flop_cg = dble(niter+1) * 15d0 * dble(n)
      time1 = MPI_Wtime()
      time1 = time1-time0
      if (time1 .gt. 0) mflops = (flop_a+flop_cg)/(1.d6*time1)
      if (pe_rank .eq. 0) then
         write(6,*)
         write(6,1) nelt,pe_size,nx1
         write(6,2) mflops, mflops/pe_size
         write(6,3) flop_a,flop_cg
         write(6,4) ksp_mon%res_start, ksp_mon%res_final
         write(6,5) time1
      endif
 1    format('nelt = ',i7, ', np = ', i9,', nx1 = ', i7)
 2    format('Tot MFlops = ', 1pe12.4, ', MFlops      = ', e12.4)
 3    format('Setup Flop = ', 1pe12.4, ', Solver Flop = ', e12.4)
 4    format('Start res  = ', 1pe12.4, ', Final res   = ', e12.4)
 5    format('Solve Time = ', e12.4)
   endif
 
 end subroutine set_timer_flop_cnt
 

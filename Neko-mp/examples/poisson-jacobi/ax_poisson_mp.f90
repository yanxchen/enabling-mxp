module ax_poisson_mp
  use ax_product_mp
  use mxm_wrapper
  use mesh
  use space_mp
  use coef_mp
  use math
  use num_types
  implicit none

  ! Single precision Poisson operator
  type, public, extends(ax_mp_t) :: ax_poisson_mp_t
  contains
     procedure, nopass :: compute_sp => ax_poisson_compute_sp
     procedure, nopass :: compute_dp => ax_poisson_compute_dp
     procedure :: compute_vector_sp => ax_poisson_compute_vector_sp
     procedure :: compute_vector_dp => ax_poisson_compute_vector_dp
  end type ax_poisson_mp_t
  
contains

  ! Single precision implementation
  subroutine ax_poisson_compute_sp(w, u, coef, msh, Xh)
    type(mesh_t), intent(inout) :: msh
    type(space_sp_t), intent(inout) :: Xh
    type(coef_sp_t), intent(inout) :: coef
    real(kind=sp), intent(inout) :: w(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=sp), intent(inout) :: u(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
  
    real(kind=sp) :: ur(Xh%lx**3)
    real(kind=sp) :: us(Xh%lx**3)
    real(kind=sp) :: ut(Xh%lx**3)
    real(kind=sp) :: wk(Xh%lx**3)
    integer :: e
  
    do e = 1, msh%nelv
       call ax_e_sp(w(1,1,1,e), u(1,1,1,e), &
            coef%g11(1,1,1,e), coef%g22(1,1,1,e), coef%g33(1,1,1,e), &
            ur, us, ut, wk, Xh%lx, Xh%dx, Xh%dxt)
    end do
  end subroutine ax_poisson_compute_sp

  ! Double precision implementation
  subroutine ax_poisson_compute_dp(w, u, coef, msh, Xh)
    type(mesh_t), intent(inout) :: msh
    type(space_dp_t), intent(inout) :: Xh
    type(coef_dp_t), intent(inout) :: coef
    real(kind=dp), intent(inout) :: w(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=dp), intent(inout) :: u(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
  
    real(kind=dp) :: ur(Xh%lx**3)
    real(kind=dp) :: us(Xh%lx**3)
    real(kind=dp) :: ut(Xh%lx**3)
    real(kind=dp) :: wk(Xh%lx**3)
    integer :: e
  
    do e = 1, msh%nelv
       call ax_e_dp(w(1,1,1,e), u(1,1,1,e), &
            coef%g11(1,1,1,e), coef%g22(1,1,1,e), coef%g33(1,1,1,e), &
            ur, us, ut, wk, Xh%lx, Xh%dx, Xh%dxt)
    end do
  end subroutine ax_poisson_compute_dp

  ! Single precision element-wise operations
  subroutine ax_e_sp(w, u, g11, g22, g33, ur, us, ut, wk, lx, D, Dt)
    integer, intent(inout) :: lx
    real(kind=sp), intent(inout) :: w(lx**3)
    real(kind=sp), intent(inout) :: u(lx**3)
    real(kind=sp), intent(inout) :: g11(lx**3)
    real(kind=sp), intent(inout) :: g22(lx**3)
    real(kind=sp), intent(inout) :: g33(lx**3)
    real(kind=sp), intent(inout) :: ur(lx**3)
    real(kind=sp), intent(inout) :: us(lx**3)
    real(kind=sp), intent(inout) :: ut(lx**3)
    real(kind=sp), intent(inout) :: wk(lx**3)
    real(kind=sp), intent(inout) :: D(lx, lx)
    real(kind=sp), intent(inout) :: Dt(lx, lx)
    real(kind=sp) :: wr, ws, wt
    integer :: i, n
  
    n = lx - 1
    call local_grad3_sp(ur, us, ut, u, n, D, Dt)
  
    do i=1, lx**3
       wr = g11(i)*ur(i) 
       ws = g22(i)*us(i) 
       wt = g33(i)*ut(i) 
       ur(i) = wr
       us(i) = ws
       ut(i) = wt
    enddo
  
    call local_grad3_t_sp(w, ur, us, ut, n, D, Dt, wk)
  end subroutine ax_e_sp

  ! Double precision element-wise operations
  subroutine ax_e_dp(w, u, g11, g22, g33, ur, us, ut, wk, lx, D, Dt)
    integer, intent(inout) :: lx
    real(kind=dp), intent(inout) :: w(lx**3)
    real(kind=dp), intent(inout) :: u(lx**3)
    real(kind=dp), intent(inout) :: g11(lx**3)
    real(kind=dp), intent(inout) :: g22(lx**3)
    real(kind=dp), intent(inout) :: g33(lx**3)
    real(kind=dp), intent(inout) :: ur(lx**3)
    real(kind=dp), intent(inout) :: us(lx**3)
    real(kind=dp), intent(inout) :: ut(lx**3)
    real(kind=dp), intent(inout) :: wk(lx**3)
    real(kind=dp), intent(inout) :: D(lx, lx)
    real(kind=dp), intent(inout) :: Dt(lx, lx)
    real(kind=dp) :: wr, ws, wt
    integer :: i, n
  
    n = lx - 1
    call local_grad3_dp(ur, us, ut, u, n, D, Dt)
  
    do i=1, lx**3
       wr = g11(i)*ur(i) 
       ws = g22(i)*us(i) 
       wt = g33(i)*ut(i) 
       ur(i) = wr
       us(i) = ws
       ut(i) = wt
    enddo
  
    call local_grad3_t_dp(w, ur, us, ut, n, D, Dt, wk)
  end subroutine ax_e_dp

  ! Single precision gradient operations
  subroutine local_grad3_sp(ur, us, ut, u, n, D, Dt)
    integer, intent(inout) :: n    
    real(kind=sp), intent(inout) :: ur(0:n, 0:n, 0:n)
    real(kind=sp), intent(inout) :: us(0:n, 0:n, 0:n)
    real(kind=sp), intent(inout) :: ut(0:n, 0:n, 0:n)
    real(kind=sp), intent(inout) :: u(0:n, 0:n, 0:n)
    real(kind=sp), intent(inout) :: D(0:n, 0:n)
    real(kind=sp), intent(inout) :: Dt(0:n, 0:n)
    integer :: m1, m2, k
  
    m1 = n + 1
    m2 = m1*m1
  
    call mxm_sp(D, m1, u, m1, ur, m2)
    do k=0,n
       call mxm_sp(u(0,0,k), m1, Dt, m1, us(0,0,k), m1)
    enddo
    call mxm_sp(u, m2, Dt, m1, ut, m1)
  end subroutine local_grad3_sp

  ! Double precision gradient operations
  subroutine local_grad3_dp(ur, us, ut, u, n, D, Dt)
    integer, intent(inout) :: n    
    real(kind=dp), intent(inout) :: ur(0:n, 0:n, 0:n)
    real(kind=dp), intent(inout) :: us(0:n, 0:n, 0:n)
    real(kind=dp), intent(inout) :: ut(0:n, 0:n, 0:n)
    real(kind=dp), intent(inout) :: u(0:n, 0:n, 0:n)
    real(kind=dp), intent(inout) :: D(0:n, 0:n)
    real(kind=dp), intent(inout) :: Dt(0:n, 0:n)
    integer :: m1, m2, k
  
    m1 = n + 1
    m2 = m1*m1
  
    call mxm(D, m1, u, m1, ur, m2)
    do k=0,n
       call mxm(u(0,0,k), m1, Dt, m1, us(0,0,k), m1)
    enddo
    call mxm(u, m2, Dt, m1, ut, m1)
  end subroutine local_grad3_dp

  ! Single precision transpose gradient operations
  subroutine local_grad3_t_sp(u, ur, us, ut, n, D, Dt, w)
    integer, intent(inout) :: n    
    real(kind=sp), intent(inout) :: ur(0:n, 0:n, 0:n)
    real(kind=sp), intent(inout) :: us(0:n, 0:n, 0:n)
    real(kind=sp), intent(inout) :: ut(0:n, 0:n, 0:n)
    real(kind=sp), intent(inout) :: u(0:n, 0:n, 0:n)
    real(kind=sp), intent(inout) :: D(0:n, 0:n)
    real(kind=sp), intent(inout) :: Dt(0:n, 0:n)
    real(kind=sp), intent(inout) :: w(0:n, 0:n, 0:n)
    integer :: m1, m2, m3, k
  
    m1 = n + 1
    m2 = m1*m1
    m3 = m1*m2
  
    call mxm_sp(Dt, m1, ur, m1, u, m2)
    do k=0,n
       call mxm_sp(us(0,0,k), m1, D, m1, w(0,0,k), m1)
    enddo
    call add2_sp(u, w, m3)
    call mxm_sp(ut, m2, D, m1, w, m1)
    call add2_sp(u, w, m3)
  end subroutine local_grad3_t_sp

  ! Double precision transpose gradient operations
  subroutine local_grad3_t_dp(u, ur, us, ut, n, D, Dt, w)
    integer, intent(inout) :: n    
    real(kind=dp), intent(inout) :: ur(0:n, 0:n, 0:n)
    real(kind=dp), intent(inout) :: us(0:n, 0:n, 0:n)
    real(kind=dp), intent(inout) :: ut(0:n, 0:n, 0:n)
    real(kind=dp), intent(inout) :: u(0:n, 0:n, 0:n)
    real(kind=dp), intent(inout) :: D(0:n, 0:n)
    real(kind=dp), intent(inout) :: Dt(0:n, 0:n)
    real(kind=dp), intent(inout) :: w(0:n, 0:n, 0:n)
    integer :: m1, m2, m3, k
  
    m1 = n + 1
    m2 = m1*m1
    m3 = m1*m2
  
    call mxm(Dt, m1, ur, m1, u, m2)
    do k=0,n
       call mxm(us(0,0,k), m1, D, m1, w(0,0,k), m1)
    enddo
    call add2(u, w, m3)
    call mxm(ut, m2, D, m1, w, m1)
    call add2(u, w, m3)
  end subroutine local_grad3_t_dp

  ! Vector compute methods (stubs - implement if needed)
  subroutine ax_poisson_compute_vector_sp(this, au, av, aw, u, v, w, coef, msh, Xh)
    class(ax_poisson_mp_t), intent(in) :: this
    type(space_sp_t), intent(inout) :: Xh
    type(mesh_t), intent(inout) :: msh
    type(coef_sp_t), intent(inout) :: coef
    real(kind=sp), intent(inout) :: au(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=sp), intent(inout) :: av(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=sp), intent(inout) :: aw(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=sp), intent(inout) :: u(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=sp), intent(inout) :: v(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=sp), intent(inout) :: w(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    ! Implement if vector operations are needed
  end subroutine ax_poisson_compute_vector_sp

  subroutine ax_poisson_compute_vector_dp(this, au, av, aw, u, v, w, coef, msh, Xh)
    class(ax_poisson_mp_t), intent(in) :: this
    type(space_dp_t), intent(inout) :: Xh
    type(mesh_t), intent(inout) :: msh
    type(coef_dp_t), intent(inout) :: coef
    real(kind=dp), intent(inout) :: au(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=dp), intent(inout) :: av(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=dp), intent(inout) :: aw(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=dp), intent(inout) :: u(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=dp), intent(inout) :: v(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=dp), intent(inout) :: w(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    ! Implement if vector operations are needed
  end subroutine ax_poisson_compute_vector_dp

end module ax_poisson_mp

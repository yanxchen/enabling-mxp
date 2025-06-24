c-----------------------------------------------------------------------
      subroutine cg(x,f,g,c,r,w,p,z,n,niter,flop_cg)
      include 'SIZE'

      include 'HSMG'

c     Solve Ax=f where A is SPD and is invoked by ax()
c
c     Output:  x - vector of length n
c
c     Input:   f - vector of length n
c     Input:   g - geometric factors for SEM operator
c     Input:   c - inverse of the counting matrix
c
c     Work arrays:   r,w,p,z  - vectors of length n
c
c     User-provided ax(w,z,n) returns  w := Az,  
c
c     User-provided solveM(z,r,n) ) returns  z := M^-1 r,  
c

      common /mymask/cmask(-1:lx1*ly1*lz1*lelt)
      parameter (lt=lx1*ly1*lz1*lelt)
      real ur(lt),us(lt),ut(lt),wk(lt)

      real x(n),f(n),r(n),w(n),p(n),z(n),g(1),c(n)

      real*4 ur_f(lt),us_f(lt),ut_f(lt),wk_f(lt)
      real*4 x_f(n),f_f(n),r_f(n),w_f(n),p_f(n),z_f(n),c_f(n)
      real*4 alpha_f, alphm_f, rtz1_f, rtz2_f, beta_f, rnorm_f, pap_f
      real*4 rtr_f

      character*1 ans
      INTEGER I

      CHARACTER*20 filename
      INTEGER unitnumber

      pap = 0.0

      alpha_f = 0.0
      alphm_f = 0.0
      rtz1_f  = 0.0
      rtz2_f  = 0.0
      beta_f  = 0.0
      rnorm_f = 0.0
      pap_f   = 0.0
      rtr_f   = 0.0

!       DO 10 I = 1, n
!             x_f(I) = 0.0
!   10  CONTINUE

c     set machine tolerances
      one = 1.
      eps = 1.e-20
      if (one+eps .eq. one) eps = 1.e-14
      if (one+eps .eq. one) eps = 1.e-7

      rtz1=1.0

      call rzero(x,n)
      call copy (r,f,n)
      call maskit (r,cmask,nx1,ny1,nz1) ! Zero out Dirichlet conditions

      rnorm = sqrt(glsc3(r,c,r,n))
      iter = 0
      if (nid.eq.0)  write(6,6) iter,rnorm

      miter = niter
c     call tester(z,r,n)  

c========================================================================================
c call iterations w double for a few iters
c========================================================================================
!       do iter=1,10
!          call solveM(z,r,n)    ! preconditioner here

!          rtz2=rtz1                                                       ! OPS
!          rtz1=glsc3(r,c,z,n)   ! parallel weighted inner product r^T C z ! 3n

!          beta = rtz1/rtz2
!          if (iter.eq.1) beta=0.0
!          call add2s1(p,z,beta,n)                                         ! 2n

!          call ax(w,p,g,ur,us,ut,wk,n)                                    ! flopa
!          pap=glsc3(w,c,p,n)                                              ! 3n

!          alpha=rtz1/pap
!          alphm=-alpha
!          call add2s2(x,p,alpha,n)                                        ! 2n         
!          call add2s2(r,w,alphm,n)                                        ! 2n

!          rtr = glsc3(r,c,r,n)                                            ! 3n
!          if (iter.eq.1) rlim2 = rtr*eps**2
!          if (iter.eq.1) rtr0  = rtr
!          rnorm = sqrt(rtr)
! c        if (nid.eq.0.and.mod(iter,100).eq.0) 
!           write(6,6) iter,rnorm,alpha,beta,pap
!     6    format('cg:',i4,1p4e12.4)
! c        if (rtr.le.rlim2) goto 1001
!       enddo

c========================================================================================
c copy from double to single
c========================================================================================

      do i=1,n
         x_f(i)=x(i)
         w_f(i)=w(i)
         z_f(i)=z(i)
         r_f(i)=r(i)
         p_f(i)=p(i)
         c_f(i)=c(i)
      end do

      call convert_dp_to_sp(ur,ur_f,lt) !
      call convert_dp_to_sp(us,us_f,lt) !
      call convert_dp_to_sp(ut,ut_f,lt) !
      call convert_dp_to_sp(wk,wk_f,lt) !

      alpha_f = alpha
      alphm_f = alphm
      rtz1_f  = rtz1
      rtz2_f  = rtz2
      beta_f  = beta
      rnorm_f = rnorm
      pap_f   = pap

c========================================================================================
c call iterations w single for rest iters
c========================================================================================
      do iter=1,miter
         !do i=1,n
         !   z(i)=z_f(i)
         !   r(i)=r_f(i)
         !end do

         !call solveM(z,r,n)    ! preconditioner here, double precision
          
         call solveM_f(z_f,r_f,n)
         
         !do i=1,n
         !   z_f(i)=z(i)
         !   r_f(i)=r(i)
         !end do

         rtz2_f = rtz1_f                                                       ! OPS
         rtz1_f = glsc3_f(r_f,c_f,z_f,n)   ! parallel weighted inner product r^T C z ! 3n, single precision
         !rtz1_f = glsc3_eftdot_mp(r_f,c_f,z_f,n)

         beta_f = rtz1_f/rtz2_f
         if (iter.eq.1) beta_f = 0.0
         call add2s1_f(p_f,z_f,beta_f,n) ! double precision

         call ax_f(w_f,p_f,g,ur_f,us_f,ut_f,wk_f,n) ! single precision

         pap_f = glsc3_f(w_f,c_f,p_f,n)            ! single precision                        ! 3n15
         !pap_f = glsc3_eftdot_mp(w_f,c_f,p_f,n)
         alpha_f = rtz1_f/pap_f
         alphm_f = -alpha_f
      
         call add2s2_f(x_f,p_f,alpha_f,n)      ! single precison                         ! 2n
         
         call add2s2_f(r_f,w_f,alphm_f,n)                                                ! 2n

         rtr_f = glsc3_f(r_f,c_f,r_f,n)                                                    ! 3n
         !rtr_f = glsc3_eftdot_mp(r_f,c_f,r_f,n)                                                    ! 3n
         
         if (iter.eq.1) rlim2 = rtr*eps**2
         if (iter.eq.1) rtr0  = rtr
         rnorm_f = sqrt(rtr_f)

         if (nid.eq.0) write(6,6) iter,rnorm_f,alpha_f,beta_f,pap_f
    6    format('cg:',i4,1p4e12.4)
c        if (rtr.le.rlim2) goto 1001
         if (rnorm_f.le.1e-10) goto 1001
      enddo

 1001 continue

      ! if (nid.eq.0) write(6,6) iter,rnorm,alpha,beta,pap

      flop_cg = flop_cg + iter*15.*n

      return
      end
c-----------------------------------------------------------------------
      subroutine solveM(z,r,n)
      include 'INPUT'
      real z(n),r(n)

      nn = n
      call h1mg_solve(z,r,nn)

      return
      end
c-----------------------------------------------------------------------
      subroutine solveM_f(z,r,n)
      include 'INPUT'
      real*4 z(n),r(n)

      nn = n
      call h1mg_solve_f(z,r,nn)

      return
      end
c-----------------------------------------------------------------------
      subroutine ax(w,u,gxyz,ur,us,ut,wk,n) ! Matrix-vector product: w=A*u in double precision

      include 'SIZE'
      include 'TOTAL'

      real w(nx1*ny1*nz1,nelt),u(nx1*ny1*nz1,nelt)
      real gxyz(2*ldim,nx1*ny1*nz1,nelt)

      parameter (lt=lx1*ly1*lz1*lelt)
      real ur(lt),us(lt),ut(lt),wk(lt)
      common /mymask/cmask(-1:lx1*ly1*lz1*lelt)

      integer e

      do e=1,nelt                                ! ~
         call ax_e( w(1,e),u(1,e),gxyz(1,1,e)    ! w   = A  u
     $                             ,ur,us,ut,wk) !  L     L  L
      enddo                                      ! 

      call dssum(w)         ! Gather-scatter operation  ! w   = QQ  w
                                                           !            L
      call add2s2(w,u,.1,n)   !2n
      call maskit(w,cmask,nx1,ny1,nz1)  ! Zero out Dirichlet conditions

      nxyz=nx1*ny1*nz1
      flop_a = flop_a + (19*nxyz+12*nx1*nxyz)*nelt

      return
      end
c-----------------------------------------------------------------------
      subroutine ax_f_1(w,u,gxyz,ur,us,ut,wk,n) ! Matrix-vector product: w=A*u in single precision

      include 'SIZE'
      include 'TOTAL'

      real*4 w(nx1*ny1*nz1,nelt),u(nx1*ny1*nz1,nelt)
      real gxyz(2*ldim,nx1*ny1*nz1,nelt)

      parameter (lt=lx1*ly1*lz1*lelt)
      real*4 ur(lt),us(lt),ut(lt),wk(lt)
      common /mymask/cmask(-1:lx1*ly1*lz1*lelt)

      real w_d(nx1*ny1*nz1,nelt),u_d(nx1*ny1*nz1,nelt)

      integer e

      do e=1,nelt                                ! ~
         call ax_e_f( w(1,e),u(1,e),gxyz(1,1,e)    ! w   = A  u
     $                             ,ur,us,ut,wk) !  L     L  L
      enddo                                      ! 
      
      do j=1,nx1*ny1*nz1
         do k=1,nelt
            w_d(j,k)=w(j,k)
            u_d(j,k)=u(j,k)
         enddo
      enddo

      call dssum(w_d)         ! Gather-scatter operation  ! w   = QQ  w
                                                           !            L
      call add2s2(w_d,u_d,.1,n)   !2n

      call maskit(w_d,cmask,nx1,ny1,nz1)  ! Zero out Dirichlet conditions

      do j=1,nx1*ny1*nz1
         do k=1,nelt
            w(j,k)=w_d(j,k)
            u(j,k)=u_d(j,k)
         enddo
      enddo

      nxyz=nx1*ny1*nz1
      flop_a = flop_a + (19*nxyz+12*nx1*nxyz)*nelt

      return
      end
c-------------------------------------------------------------------------
      subroutine ax_f(w,u,gxyz,ur,us,ut,wk,n) ! Matrix-vector product: w=A*u in single precision

      include 'SIZE'
      include 'TOTAL'

      real*4 w(nx1*ny1*nz1,nelt),u(nx1*ny1*nz1,nelt)
      real gxyz(2*ldim,nx1*ny1*nz1,nelt)

      real w_d(nx1*ny1*nz1,nelt)
      parameter (lt=lx1*ly1*lz1*lelt)
      real*4 ur(lt),us(lt),ut(lt),wk(lt)
      common /mymask/cmask(-1:lx1*ly1*lz1*lelt)

      real*4 alpha

      integer e

      do e=1,nelt                                ! ~
         call ax_e_f( w(1,e),u(1,e),gxyz(1,1,e)    ! w   = A  u
     $                             ,ur,us,ut,wk) !  L     L  L
      enddo                                      ! 

      do i=1,nx1*ny1*nz1
         do j=1,nelt
            w_d(i,j)=w(i,j)
         enddo
      enddo
      call dssum(w_d)
      do i=1,nx1*ny1*nz1
         do j=1,nelt
            w(i,j)=w_d(i,j)
         enddo
      enddo
 

      !call dssum_f(w)         ! Gather-scatter operation  ! w   = QQ  w
                                                           !            L
      alpha = 0.1
      call add2s2_f(w,u,alpha,n)   !2n

      call maskit_f(w,cmask,nx1,ny1,nz1)  ! Zero out Dirichlet conditions

      nxyz=nx1*ny1*nz1
      flop_a = flop_a + (19*nxyz+12*nx1*nxyz)*nelt

      return
      end
c-------------------------------------------------------------------------
      subroutine ax1(w,u,n)
      include 'SIZE'
      real w(n),u(n)
      real h2i
  
      h2i = (n+1)*(n+1)  
      do i = 2,n-1
         w(i)=h2i*(2*u(i)-u(i-1)-u(i+1))
      enddo
      w(1)  = h2i*(2*u(1)-u(2  ))
      w(n)  = h2i*(2*u(n)-u(n-1))

      return
      end
c-------------------------------------------------------------------------
      subroutine ax_e(w,u,g,ur,us,ut,wk) ! Local matrix-vector product
      include 'SIZE'
      include 'TOTAL'

      parameter (lxyz=lx1*ly1*lz1)
      real ur(lxyz),us(lxyz),ut(lxyz),wk(lxyz)
      real w(nx1*ny1*nz1),u(nx1*ny1*nz1),g(2*ldim,nx1*ny1*nz1)


      nxyz = nx1*ny1*nz1
      n    = nx1-1

      call local_grad3(ur,us,ut,u,n,dxm1,dxtm1)

      do i=1,nxyz
         wr = g(1,i)*ur(i) + g(2,i)*us(i) + g(3,i)*ut(i)
         ws = g(2,i)*ur(i) + g(4,i)*us(i) + g(5,i)*ut(i)
         wt = g(3,i)*ur(i) + g(5,i)*us(i) + g(6,i)*ut(i)
         ur(i) = wr
         us(i) = ws
         ut(i) = wt
      enddo

      call local_grad3_t(w,ur,us,ut,n,dxm1,dxtm1,wk)

      return
      end
c-------------------------------------------------------------------------
      subroutine ax_e_f(w,u,g,ur,us,ut,wk) ! Local matrix-vector product in single precision
      include 'SIZE'
      include 'TOTAL'

      parameter (lxyz=lx1*ly1*lz1)
      real*4 ur(lxyz),us(lxyz),ut(lxyz),wk(lxyz)
      real*4 w(nx1*ny1*nz1),u(nx1*ny1*nz1)
      real g(2*ldim,nx1*ny1*nz1)
      real*4 wr,ws,wt

      nxyz = nx1*ny1*nz1
      n    = nx1-1

      call local_grad3_f(ur,us,ut,u,n,dxm1_f,dxtm1_f)

! mixed datatypes calculation
      do i=1,nxyz
         wr = g(1,i)*ur(i) + g(2,i)*us(i) + g(3,i)*ut(i)
         ws = g(2,i)*ur(i) + g(4,i)*us(i) + g(5,i)*ut(i)
         wt = g(3,i)*ur(i) + g(5,i)*us(i) + g(6,i)*ut(i)
         ur(i) = wr
         us(i) = ws
         ut(i) = wt
      enddo

      call local_grad3_t_f(w,ur,us,ut,n,dxm1_f,dxtm1_f,wk)

      return
      end
c-------------------------------------------------------------------------
      subroutine local_grad3(ur,us,ut,u,n,D,Dt)
c     Output: ur,us,ut         Input:u,n,D,Dt
      real ur(0:n,0:n,0:n),us(0:n,0:n,0:n),ut(0:n,0:n,0:n)
      real u (0:n,0:n,0:n)
      real D (0:n,0:n),Dt(0:n,0:n)
      integer e

      m1 = n+1
      m2 = m1*m1

      call mxm(D ,m1,u,m1,ur,m2)
      do k=0,n
         call mxm(u(0,0,k),m1,Dt,m1,us(0,0,k),m1)
      enddo
      call mxm(u,m2,Dt,m1,ut,m1)

      return
      end
c-------------------------------------------------------------------------
      subroutine local_grad3_f(ur,us,ut,u,n,D,Dt)
c     Output: ur,us,ut         Input:u,n,D,Dt
      real*4 ur(0:n,0:n,0:n),us(0:n,0:n,0:n),ut(0:n,0:n,0:n)
      real*4 u (0:n,0:n,0:n)
      real*4 D (0:n,0:n),Dt(0:n,0:n)
      integer e

      m1 = n+1
      m2 = m1*m1

      call mxm_f(D ,m1,u,m1,ur,m2)
      do k=0,n
         call mxm_f(u(0,0,k),m1,Dt,m1,us(0,0,k),m1)
      enddo
      call mxm_f(u,m2,Dt,m1,ut,m1)

      return
      end
c-----------------------------------------------------------------------
      subroutine local_grad3_t(u,ur,us,ut,N,D,Dt,w)
c     Output: ur,us,ut         Input:u,N,D,Dt
      real u (0:N,0:N,0:N)
      real ur(0:N,0:N,0:N),us(0:N,0:N,0:N),ut(0:N,0:N,0:N)
      real D (0:N,0:N),Dt(0:N,0:N)
      real w (0:N,0:N,0:N)
      integer e

      m1 = N+1
      m2 = m1*m1
      m3 = m1*m1*m1

      call mxm(Dt,m1,ur,m1,u,m2)

      do k=0,N
         call mxm(us(0,0,k),m1,D ,m1,w(0,0,k),m1)
      enddo
      call add2(u,w,m3)

      call mxm(ut,m2,D ,m1,w,m1)
      call add2(u,w,m3)

      return
      end
c-----------------------------------------------------------------------
      subroutine local_grad3_t_f(u,ur,us,ut,N,D,Dt,w)
c     Output: ur,us,ut         Input:u,N,D,Dt
      real*4 u (0:N,0:N,0:N)
      real*4 ur(0:N,0:N,0:N),us(0:N,0:N,0:N),ut(0:N,0:N,0:N)
      real*4 D (0:N,0:N),Dt(0:N,0:N)
      real*4 w (0:N,0:N,0:N)
      integer e

      m1 = N+1
      m2 = m1*m1
      m3 = m1*m1*m1

      call mxm_f(Dt,m1,ur,m1,u,m2)

      do k=0,N
         call mxm_f(us(0,0,k),m1,D ,m1,w(0,0,k),m1)
      enddo
      call add2_f(u,w,m3)

      call mxm_f(ut,m2,D ,m1,w,m1)
      call add2_f(u,w,m3)

      return
      end
c-----------------------------------------------------------------------
      subroutine maskit(w,pmask,nx,ny,nz)   ! Zero out Dirichlet conditions
      include 'SIZE'
      include 'PARALLEL'

      real pmask(-1:lx1*ly1*lz1*lelt)
      real w(1)
      integer e

      nxyz = nx*ny*nz
      nxy  = nx*ny
      if(pmask(-1).lt.0) then
        j=pmask(0)
        do i = 1,j
           k = pmask(i)
           w(k)=0.0
        enddo
      else
c         Zero out Dirichlet boundaries.
c
c                      +------+     ^ Y
c                     /   3  /|     |
c               4--> /      / |     |
c                   +------+ 2 +    +----> X
c                   |   5  |  /    /
c                   |      | /    /
c                   +------+     Z   
c

        nn = 0
        do e  = 1,nelt
          call get_face(w,nx,e)
          do i = 1,nxyz
             if(w(i).eq.0) then
               nn=nn+1
               pmask(nn)=i
             endif
          enddo
        enddo     
        pmask(-1) = -1.
        pmask(0) = nn
      endif


      return
      end
c-----------------------------------------------------------------------
      subroutine maskit_f(w,pmask,nx,ny,nz)   ! Zero out Dirichlet conditions
      include 'SIZE'
      include 'PARALLEL'

      real pmask(-1:lx1*ly1*lz1*lelt)
      real*4 w(1)
      integer e

      nxyz = nx*ny*nz
      nxy  = nx*ny
      if(pmask(-1).lt.0) then
        j=pmask(0)
        do i = 1,j
           k = pmask(i)
           w(k)=0.0
        enddo
      else
c         Zero out Dirichlet boundaries.
c
c                      +------+     ^ Y
c                     /   3  /|     |
c               4--> /      / |     |
c                   +------+ 2 +    +----> X
c                   |   5  |  /    /
c                   |      | /    /
c                   +------+     Z   
c

        nn = 0
        do e  = 1,nelt
          call get_face(w,nx,e)
          do i = 1,nxyz
             if(w(i).eq.0) then
               nn=nn+1
               pmask(nn)=i
             endif
          enddo
        enddo     
        pmask(-1) = -1.
        pmask(0) = nn
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine masko(w)   ! Old 'mask'
      include 'SIZE'
      real w(1)

      if (nid.eq.0) w(1) = 0.  ! suitable for solvability

      return
      end
c-----------------------------------------------------------------------
      subroutine masking(w,nx,e,x0,x1,y0,y1,z0,z1)
c     Zeros out boundary
      include 'SIZE'
      integer e,x0,x1,y0,y1,z0,z1
      real w(nx,nx,nx,nelt)
      
c       write(6,*) x0,x1,y0,y1,z0,z1
      do k=z0,z1
      do j=y0,y1
      do i=x0,x1
          w(i,j,k,e)=0.0
      enddo
      enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine masking_f(w,nx,e,x0,x1,y0,y1,z0,z1)
c     Zeros out boundary
      include 'SIZE'
      integer e,x0,x1,y0,y1,z0,z1
      real*4 w(nx,nx,nx,nelt)
      
c       write(6,*) x0,x1,y0,y1,z0,z1
      do k=z0,z1
      do j=y0,y1
      do i=x0,x1
          w(i,j,k,e)=0.0
      enddo
      enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine tester(z,r,n)
c     Used to test if solution to precond. is SPD
      real r(n),z(n)

      do j=1,n
         call rzero(r,n)
         r(j) = 1.0
         call solveM(z,r,n)
         do i=1,n
            write(79,*) z(i)
         enddo
      enddo
      call exitt0
      return
      end
c-----------------------------------------------------------------------
      subroutine get_face(w,nx,ie)
c     zero out all boundaries as Dirichlet
c     to change, change this routine to only zero out 
c     the nodes desired to be Dirichlet, and leave others alone.
      include 'SIZE'
      include 'PARALLEL'
      real w(1)
      integer nx,ie,nelx,nely,nelz
      integer x0,x1,y0,y1,z0,z1
      
      x0=1
      y0=1
      z0=1
      x1=nx
      y1=nx
      z1=nx
      
      nelxy=nelx*nely
      ngl = lglel(ie)        !global element number

      ir = 1+(ngl-1)/nelxy   !global z-count
      iq = mod1(ngl,nelxy)   !global y-count
      iq = 1+(iq-1)/nelx     
      ip = mod1(ngl,nelx)    !global x-count

c     write(6,1) ip,iq,ir,nelx,nely,nelz, nelt,' test it'
c  1  format(7i7,a8)

      if(mod(ip,nelx).eq.1.or.nelx.eq.1)   then  ! Face4
         x0=1
         x1=1
         call masking(w,nx,ie,x0,x1,y0,y1,z0,z1)
      endif
      if(mod(ip,nelx).eq.0)               then   ! Face2
         x0=nx
         x1=nx
         call masking(w,nx,ie,x0,x1,y0,y1,z0,z1)
      endif

      x0=1
      x1=nx
      if(mod(iq,nely).eq.1.or.nely.eq.1) then    ! Face1
         y0=1
         y1=1
         call masking(w,nx,ie,x0,x1,y0,y1,z0,z1)
      endif
      if(mod(iq,nely).eq.0)              then    ! Face3
         y0=nx
         y1=nx
         call masking(w,nx,ie,x0,x1,y0,y1,z0,z1)
      endif

      y0=1
      y1=nx
      if(mod(ir,nelz).eq.1.or.nelz.eq.1) then    ! Face5
         z0=1
         z1=1
         call masking(w,nx,ie,x0,x1,y0,y1,z0,z1)
      endif
      if(mod(ir,nelz).eq.0)              then    ! Face6
         z1=nx
         z0=nx
         call masking(w,nx,ie,x0,x1,y0,y1,z0,z1)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine get_face_f(w,nx,ie)
c     zero out all boundaries as Dirichlet
c     to change, change this routine to only zero out 
c     the nodes desired to be Dirichlet, and leave others alone.
      include 'SIZE'
      include 'PARALLEL'
      real*4 w(1)
      integer nx,ie,nelx,nely,nelz
      integer x0,x1,y0,y1,z0,z1
      
      x0=1
      y0=1
      z0=1
      x1=nx
      y1=nx
      z1=nx
      
      nelxy=nelx*nely
      ngl = lglel(ie)        !global element number

      ir = 1+(ngl-1)/nelxy   !global z-count
      iq = mod1(ngl,nelxy)   !global y-count
      iq = 1+(iq-1)/nelx     
      ip = mod1(ngl,nelx)    !global x-count

c     write(6,1) ip,iq,ir,nelx,nely,nelz, nelt,' test it'
c  1  format(7i7,a8)

      if(mod(ip,nelx).eq.1.or.nelx.eq.1)   then  ! Face4
         x0=1
         x1=1
         call masking_f(w,nx,ie,x0,x1,y0,y1,z0,z1)
      endif
      if(mod(ip,nelx).eq.0)               then   ! Face2
         x0=nx
         x1=nx
         call masking_f(w,nx,ie,x0,x1,y0,y1,z0,z1)
      endif

      x0=1
      x1=nx
      if(mod(iq,nely).eq.1.or.nely.eq.1) then    ! Face1
         y0=1
         y1=1
         call masking_f(w,nx,ie,x0,x1,y0,y1,z0,z1)
      endif
      if(mod(iq,nely).eq.0)              then    ! Face3
         y0=nx
         y1=nx
         call masking_f(w,nx,ie,x0,x1,y0,y1,z0,z1)
      endif

      y0=1
      y1=nx
      if(mod(ir,nelz).eq.1.or.nelz.eq.1) then    ! Face5
         z0=1
         z1=1
         call masking_f(w,nx,ie,x0,x1,y0,y1,z0,z1)
      endif
      if(mod(ir,nelz).eq.0)              then    ! Face6
         z1=nx
         z0=nx
         call masking_f(w,nx,ie,x0,x1,y0,y1,z0,z1)
      endif

      return
      end
c-----------------------------------------------------------------------

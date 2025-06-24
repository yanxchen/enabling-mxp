C======================================================================
c=======================================================================

c----------------------------------------------------------------------
c----------------------------------------------------------------------
      subroutine h1mg_dssum_f(u,l)
      include 'SIZE'
      include 'HSMG'

      real*4 u(1)
      real*8 tmp(1000000)
      integer l, n, i

      !n = mg_h1_n(l, mg_fld)
      !if (n .gt. 1000000) then 
      !  write(*,*) 'Error: tmp array too small'
      !  stop
      !end if

      !do i = 1, n
      !   tmp(i) = dble(u(i))
      !end do
      !
      !call adelay
      !call gs_op(mg_gsh_handle(l,mg_fld),tmp,1,1,0)

      !do i = 1, n
      !   u(i) = sngl(tmp(i))
      !end do
      
      call adelay
      call gs_op(mg_gsh_handle(l,mg_fld),u,2,1,0)


      return
      end
c----------------------------------------------------------------------
      subroutine h1mg_dsprod_f(u,l)
      include 'SIZE'
      include 'HSMG'
      real*4 u(1)

      call adelay
      call gs_op(mg_gsh_handle(l,mg_fld),u,2,2,0)

      return
      end
c----------------------------------------------------------------------
      subroutine h1mg_schwarz_dssum_f(u,l)
      include 'SIZE'
      include 'HSMG'

      call adelay
      call gs_op(mg_gsh_schwarz_handle(l,mg_fld),u,2,1,0)
      return
      end
c----------------------------------------------------------------------


c----------------------------------------------------------------------
      subroutine hsmg_coarse_solve_f(e,r)
      include 'SIZE'
      include 'HSMG'
      real*4 e(1),r(1)

       n=mg_nh(1)
       nxyz = n*n*n*nelt
       call copy_f(e,r,nxyz)
       call h1mg_dssum_f(e,1)


      return
      end
c----------------------------------------------------------------------
c----------------------------------------------------------------------
      subroutine h1mg_solve_f(z,rhs,nn)  !  Solve preconditioner: Mz=rhs
      real*4 z(1),rhs(1)

c     Assumes that preprocessing has been completed via h1mg_setup()

      include 'SIZE'
      include 'HSMG'       ! Same array space as HSMG
      include 'TOTAL'
      
      parameter (lt=lx1*ly1*lz1*lelt)
      common /scrmg/ e(2*lt),w(lt),r(lt)
      real*4 e,w,r
      integer p_msk
      real*4 op,om,sigma

      nel   = nelt

      op    =  1.                                     ! Coefficients for h1mg_ax
      om    = -1.
      sigma =  1.

      l     = mg_h1_lmax
      n     = mg_h1_n(l,mg_fld)
      is    = 1                                       ! solve index

      call h1mg_schwarz_f(z,rhs,sigma,l)              ! z := sigma W M       rhs
                                                      !               Schwarz
      call copy_f(r,rhs,n)                            ! r  := rhs

      do l = mg_h1_lmax-1,2,-1                        ! DOWNWARD Leg of V-cycle
         is = is + n
         n  = mg_h1_n(l,mg_fld)
                                                      !          T
         call h1mg_rstr_f(r,l,.true.)                 ! r   :=  J r
                                                      !  l         l+1
!        OVERLAPPING Schwarz exchange and solve:
         call h1mg_schwarz_f(e(is),r,sigma,l)           ! e := sigma W M       r
                                                      !  l            Schwarz l
      enddo
      
      is = is+n
                                                      !         T
      call h1mg_rstr_f(r,1,.false.)                   ! r  :=  J  r
                                                      !  l         l+1

      p_msk = p_mg_msk(l,mg_fld)
      call h1mg_mask_f(r,mg_imask(p_msk),nel)         !        -1
      call hsmg_coarse_solve_f ( e(is) , r )          ! e  := A   r
      call h1mg_mask_f(e(is),mg_imask(p_msk),nel)     !  1     1   1

c     nx = mg_nh(1)
c     call h1mg_mask(e(is),mg_imask(p_msk),nel)       !  1     1   1
c     call exitt

      do l = 2,mg_h1_lmax-1                           ! UNWIND.  No smoothing.
         im = is
         is = is - n
         n  = mg_h1_n(l,mg_fld)
         call h1mg_intp_f (w,e(im),l-1)               ! w   :=  J e
         i1=is-1                                      !            l-1
         do i=1,n
            e(i1+i) = e(i1+i) + w(i)                  ! e   :=  e  + w
         enddo                                        !  l       l
      enddo

      l  = mg_h1_lmax
      n  = mg_h1_n(l,mg_fld)
      im = is  ! solve index
      call h1mg_intp_f(w,e(im),l-1)                   ! w   :=  J e
      do i = 1,n                                      !            l-1
         z(i) = z(i) + w(i)                           ! z := z + w
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine h1mg_schwarz_f(e,r,sigma,l)
      include 'SIZE'
      include 'HSMG'

      real*4 e(1),r(1)
      real*4 sigma

      n = mg_h1_n(l,mg_fld)

      call h1mg_schwarz_wt_f    (e,l)         ! e  := W^.5* e
      call h1mg_schwarz_part1_f (e,r,l)       !  l           l
      call h1mg_schwarz_wt_f    (e,l)         ! e  := e *W^.5
      call cmult_f              (e,sigma,n)   !  l           l

      return
      end
c-----------------------------------------------------------------------
      subroutine h1mg_schwarz_part1_f (e,r,l)
      include 'SIZE'
      include 'HSMG'

      real*4 e(1),r(1)

      integer enx,eny,enz,pm

      real*4 zero,one,onem

      zero =  0
      one  =  1
      onem = -1

      n  = mg_h1_n (l,mg_fld)
      pm = p_mg_msk(l,mg_fld)

      call h1mg_mask_f (r,mg_imask(pm),nelt)  ! Zero Dirichlet nodes

      call h1mg_schwarz_toext3d_f(mg_work_f,r,mg_nh(l))

      enx=mg_nh(l)+2
      eny=mg_nh(l)+2
      enz=mg_nh(l)+2
      i = enx*eny*enz*nelt+1
 
c     exchange interior nodes
      call h1mg_extrude_f(mg_work_f,0,zero,mg_work_f,2,one,enx,eny,enz)
      call h1mg_schwarz_dssum_f(mg_work_f,l)
      call h1mg_extrude_f(mg_work_f,0,one ,mg_work_f,2,onem,enx,eny,enz)

      call h1mg_fdm_f(mg_work_f(i),mg_work_f,l) ! Do the local solves

c     Sum overlap region (border excluded)
      call h1mg_extrude_f(mg_work_f,0,zero,
     $                    mg_work_f(i),0,one ,enx,eny,enz)
      call h1mg_schwarz_dssum_f(mg_work_f(i),l)
      call h1mg_extrude_f(mg_work_f(i),0,one,
     $                    mg_work_f,0,onem,enx,eny,enz)
      call h1mg_extrude_f(mg_work_f(i),2,one,
     $                    mg_work_f(i),0,one,enx,eny,enz)

      call h1mg_schwarz_toreg3d_f(e,mg_work_f(i),mg_nh(l))

      call h1mg_dssum_f(e,l)                           ! sum border nodes
      call h1mg_mask_f (e,mg_imask(pm),nelt) ! apply mask 

      return
      end
c-----------------------------------------------------------------------
      subroutine h1mg_schwarz_wt_f(e,l)
      include 'SIZE'
      include 'INPUT'
      include 'HSMG'
      
      call h1mg_schwarz_wt3d2_mp(
     $    e,mg_schwarz_wt_f(mg_schwarz_wt_index(l,mg_fld)),mg_nh(l))
      return
      end 
c-----------------------------------------------------------------------
      subroutine h1mg_schwarz_wt3d2_f(e,wt,n)
      include 'SIZE'
      integer n
      real*4 e(n,n,n,nelt)
      real*4 wt(n,n,4,3,nelt)
      
      integer ie,i,j,k
      do ie=1,nelt
         do k=1,n
         do j=1,n
            e(1  ,j,k,ie)=e(1  ,j,k,ie)*sqrt(wt(j,k,1,1,ie))
            e(2  ,j,k,ie)=e(2  ,j,k,ie)*sqrt(wt(j,k,2,1,ie))
            e(n-1,j,k,ie)=e(n-1,j,k,ie)*sqrt(wt(j,k,3,1,ie))
            e(n  ,j,k,ie)=e(n  ,j,k,ie)*sqrt(wt(j,k,4,1,ie))
         enddo
         enddo
         do k=1,n
         do i=3,n-2
            e(i,1  ,k,ie)=e(i,1  ,k,ie)*sqrt(wt(i,k,1,2,ie))
            e(i,2  ,k,ie)=e(i,2  ,k,ie)*sqrt(wt(i,k,2,2,ie))
            e(i,n-1,k,ie)=e(i,n-1,k,ie)*sqrt(wt(i,k,3,2,ie))
            e(i,n  ,k,ie)=e(i,n  ,k,ie)*sqrt(wt(i,k,4,2,ie))
         enddo
         enddo
         do j=3,n-2
         do i=3,n-2
            e(i,j,1  ,ie)=e(i,j,1  ,ie)*sqrt(wt(i,j,1,3,ie))
            e(i,j,2  ,ie)=e(i,j,2  ,ie)*sqrt(wt(i,j,2,3,ie))
            e(i,j,n-1,ie)=e(i,j,n-1,ie)*sqrt(wt(i,j,3,3,ie))
            e(i,j,n  ,ie)=e(i,j,n  ,ie)*sqrt(wt(i,j,4,3,ie))
         enddo
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine h1mg_schwarz_wt3d2_mp(e,wt,n)
      include 'SIZE'
      integer n
      real*4 e(n,n,n,nelt)
      real*4 wt(n,n,4,3,nelt)

      integer ie,i,j,k
      do ie=1,nelt
         do k=1,n
         do j=1,n
            e(1  ,j,k,ie)=e(1  ,j,k,ie)*dsqrt(dble(wt(j,k,1,1,ie)))
            e(2  ,j,k,ie)=e(2  ,j,k,ie)*dsqrt(dble(wt(j,k,2,1,ie)))
            e(n-1,j,k,ie)=e(n-1,j,k,ie)*dsqrt(dble(wt(j,k,3,1,ie)))
            e(n  ,j,k,ie)=e(n  ,j,k,ie)*dsqrt(dble(wt(j,k,4,1,ie)))
         enddo
         enddo
         do k=1,n
         do i=3,n-2
            e(i,1  ,k,ie)=e(i,1  ,k,ie)*dsqrt(dble(wt(i,k,1,2,ie)))
            e(i,2  ,k,ie)=e(i,2  ,k,ie)*dsqrt(dble(wt(i,k,2,2,ie)))
            e(i,n-1,k,ie)=e(i,n-1,k,ie)*dsqrt(dble(wt(i,k,3,2,ie)))
            e(i,n  ,k,ie)=e(i,n  ,k,ie)*dsqrt(dble(wt(i,k,4,2,ie)))
         enddo
         enddo
         do j=3,n-2
         do i=3,n-2
            e(i,j,1  ,ie)=e(i,j,1  ,ie)*dsqrt(dble(wt(i,j,1,3,ie)))
            e(i,j,2  ,ie)=e(i,j,2  ,ie)*dsqrt(dble(wt(i,j,2,3,ie)))
            e(i,j,n-1,ie)=e(i,j,n-1,ie)*dsqrt(dble(wt(i,j,3,3,ie)))
            e(i,j,n  ,ie)=e(i,j,n  ,ie)*dsqrt(dble(wt(i,j,4,3,ie)))
         enddo
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------

      subroutine h1mg_mask_f(w,mask,nel)
      include 'SIZE'

      real*4  w   (1)
      integer mask(1)        ! Pointer to Dirichlet BCs
      integer e
      
      do e=1,nel
         im = mask(e)
         call mg_mask_e_f(w,mask(im)) ! Zero out Dirichlet conditions
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine mg_mask_e_f(w,mask) ! Zero out Dirichlet conditions
      real*4  w(1)
      integer mask(0:1)

      n=mask(0)
      do i=1,n
         w(mask(i)) = 0.
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine h1mg_schwarz_toreg3d_f(b,a,n)
c     strip off ghost cell
      include 'SIZE'
      integer n
      real*4 a(0:n+1,0:n+1,0:n+1,nelt)
      real*4 b(n,n,n,nelt)
      
      integer i,j,k,ie
      do ie=1,nelt
      do k=1,n
      do j=1,n
      do i=1,n
         b(i,j,k,ie)=a(i,j,k,ie)
      enddo
      enddo
      enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine h1mg_schwarz_toext3d_f(a,b,n)
c     border nodes (ghost cell = 0)
      include 'SIZE'
      integer n
      real*4 a(0:n+1,0:n+1,0:n+1,nelt)
      real*4 b(n,n,n,nelt)
      
      integer i,j,k,ie
      call rzero(a,(n+2)*(n+2)*(n+2)*nelt)
      do ie=1,nelt
      do k=1,n
      do j=1,n
      do i=1,n
         a(i,j,k,ie)=b(i,j,k,ie)
      enddo
      enddo
      enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine h1mg_extrude_f(arr1,l1,f1,arr2,l2,f2,nx,ny,nz)
      include 'SIZE'
      integer l1,l2,nx,ny,nz
      real*4 arr1(nx,ny,nz,nelt),arr2(nx,ny,nz,nelt)
      real*4 f1,f2
      
      integer i,j,k,ie,i0,i1
      i0=2
      i1=nx-1
      
      do ie=1,nelt
         do k=i0,i1
         do j=i0,i1
            arr1(l1+1 ,j,k,ie) = f1*arr1(l1+1 ,j,k,ie)
     $                          +f2*arr2(l2+1 ,j,k,ie)
            arr1(nx-l1,j,k,ie) = f1*arr1(nx-l1,j,k,ie)
     $                          +f2*arr2(nx-l2,j,k,ie)
         enddo
         enddo
         do k=i0,i1
         do i=i0,i1
            arr1(i,l1+1 ,k,ie) = f1*arr1(i,l1+1 ,k,ie)
     $                          +f2*arr2(i,l2+1 ,k,ie)
            arr1(i,nx-l1,k,ie) = f1*arr1(i,nx-l1,k,ie)
     $                          +f2*arr2(i,nx-l2,k,ie)
         enddo
         enddo
         do j=i0,i1
         do i=i0,i1
            arr1(i,j,l1+1 ,ie) = f1*arr1(i,j,l1+1 ,ie)
     $                          +f2*arr2(i,j,l2+1 ,ie)
            arr1(i,j,nx-l1,ie) = f1*arr1(i,j,nx-l1,ie)
     $                          +f2*arr2(i,j,nx-l2,ie)
         enddo
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
c     clobbers r
      subroutine h1mg_fdm_f(e,r,l)
      include 'SIZE'
      include 'HSMG'
      call h1mg_do_fast_f(e,r,
     $      mg_fast_s_f(mg_fast_s_index(l,mg_fld)),
     $      mg_fast_d_f(mg_fast_d_index(l,mg_fld)),
     $      mg_nh(l)+2)
      return
      end
c-----------------------------------------------------------------------
c     clobbers r
      subroutine h1mg_do_fast_f(e,r,s,d,nl)
      include 'SIZE'
      real*4 e(nl**ndim,nelt)
      real*4 r(nl**ndim,nelt)
      real*4 s(nl*nl,2,ndim,nelt)
      real*4 d(nl**ndim,nelt)
      
      integer ie,nn,i
      nn=nl**ndim

      do ie=1,nelt
         call h1mg_tnsr3d_el_f(e(1,ie),nl,r(1,ie),nl
     $                      ,s(1,2,1,ie),s(1,1,2,ie),s(1,1,3,ie))
         do i=1,nn
            r(i,ie)=d(i,ie)*e(i,ie)
         enddo
         call h1mg_tnsr3d_el_f(e(1,ie),nl,r(1,ie),nl
     $                      ,s(1,1,1,ie),s(1,2,2,ie),s(1,2,3,ie))
      enddo

      return
      end
c----------------------------------------------------------------------
c     computes
c     v = [C (x) B (x) A] u
      subroutine h1mg_tnsr3d_el_f(v,nv,u,nu,A,Bt,Ct)

      integer nv,nu
      real*4 v(nv*nv*nv),u(nu*nu*nu),A(1),Bt(1),Ct(1)
      include 'SIZE'
      parameter (lwk=(lx1+2)*(ly1+2)*(lz1+2))
      common /hsmgw/ work(0:lwk-1),work2(0:lwk-1)
      real*4 work,work2
      integer i

      call mxm_f(A,nv,u,nu,work,nu*nu)
      do i=0,nu-1
         call mxm_f(work(nv*nu*i),nv,Bt,nu,work2(nv*nv*i),nv)
      enddo
      call mxm_f(work2,nv*nv,Ct,nu,v,nv)

      return
      end
c-----------------------------------------------------------------------
      subroutine h1mg_rstr_f(r,l,ifdssum) ! r =J r,   l is coarse level
      include 'SIZE'
      include 'HSMG'
      logical ifdssum

      real*4 r(1)
      integer l

      call h1mg_do_wt_f(r,
     $                  mg_rstr_wt_f(mg_rstr_wt_index(l+1,mg_fld))
     $                  ,mg_nh(l+1),mg_nh(l+1),mg_nhz(l+1))

      call h1mg_tnsr1_f(r,mg_nh(l),mg_nh(l+1),
     $                  mg_jht_f(1,l),mg_jh_f(1,l))

      if (ifdssum) call h1mg_dssum_f(r,l)

      return
      end
c-----------------------------------------------------------------------
c     u = wt .* u [mixed types calculations]
      subroutine h1mg_do_wt_f(u,wt,nx,ny,nz)
      include 'SIZE'
      integer nx,ny,nz
      real*4 u(nx,ny,nz,nelt),wt(nx,nz,2,ndim,nelt)
      
      integer e

      do ie=1,nelt
         do k=1,nz
         do j=1,ny
            u( 1,j,k,ie)=u( 1,j,k,ie)*wt(j,k,1,1,ie)
            u(nx,j,k,ie)=u(nx,j,k,ie)*wt(j,k,2,1,ie)
         enddo
         enddo
         do k=1,nz
         do i=2,nx-1
            u(i, 1,k,ie)=u(i, 1,k,ie)*wt(i,k,1,2,ie)
            u(i,ny,k,ie)=u(i,ny,k,ie)*wt(i,k,2,2,ie)
         enddo
         enddo
         do j=2,ny-1
         do i=2,nx-1
            u(i,j, 1,ie)=u(i,j, 1,ie)*wt(i,j,1,3,ie)
            u(i,j,nz,ie)=u(i,j,nz,ie)*wt(i,j,2,3,ie)
         enddo
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine h1mg_tnsr1_f(v,nv,nu,A,At)
c
c     v = [A (x) A (x) A] u 
c
      integer nv,nu
      real*4 v(1),A(1),At(1)

      call h1mg_tnsr1_3d_f(v,nv,nu,A,At,At)

      return
      end
c-----------------------------------------------------------------------
      subroutine h1mg_tnsr1_3d_f(v,nv,nu,A,Bt,Ct) ! v = [C (x) B (x) A] u 
      integer nv,nu
      real*4 v(1),A(1),Bt(1),Ct(1)
      include 'SIZE'
      parameter (lwk=(lx1+2)*(ly1+2)*(lz1+2))
      common /hsmgw/ work(0:lwk-1),work2(0:lwk-1)
      real*4 work,work2
      integer e,e0,ee,es

      e0=1
      es=1
      ee=nelt

      if (nv.gt.nu) then
         e0=nelt
         es=-1
         ee=1
      endif

      nu3 = nu**3
      nv3 = nv**3

      do e=e0,ee,es
         iu = 1 + (e-1)*nu3
         iv = 1 + (e-1)*nv3
         call mxm_f(A,nv,v(iu),nu,work,nu*nu)
         do i=0,nu-1
            call mxm_f(work(nv*nu*i),nv,Bt,nu,work2(nv*nv*i),nv)
         enddo
         call mxm_f(work2,nv*nv,Ct,nu,v(iv),nv)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine h1mg_intp_f(uf,uc,l) ! l is coarse level
      real*4 uf(1),uc(1)
      integer l
      include 'SIZE'
      include 'HSMG'
      call h1mg_tnsr_f(uf,mg_nh(l+1),uc,
     $                 mg_nh(l),mg_jh_f(1,l),mg_jht_f(1,l))
      return
      end
c-----------------------------------------------------------------------
c     computes
c     v = [A (x) A] u      or
c     v = [A (x) A (x) A] u 
      subroutine h1mg_tnsr_f(v,nv,u,nu,A,At)
      integer nv,nu
      real*4 v(1),u(1),A(1),At(1)

      call h1mg_tnsr3d_f(v,nv,u,nu,A,At,At)

      return
      end
c-------------------------------------------------------T--------------
c     computes
c              
c     v = [C (x) B (x) A] u
      subroutine h1mg_tnsr3d_f(v,nv,u,nu,A,Bt,Ct)
      integer nv,nu
      real*4 v(nv*nv*nv,nelt),u(nu*nu*nu,nelt),A(1),Bt(1),Ct(1)
      include 'SIZE'
      parameter (lwk=(lx1+2)*(ly1+2)*(lz1+2))
      common /hsmgw/ work(0:lwk-1),work2(0:lwk-1)
      real*4 work,work2
      integer ie, i
      do ie=1,nelt
         call mxm_f(A,nv,u(1,ie),nu,work,nu*nu)
         do i=0,nu-1
            call mxm_f(work(nv*nu*i),nv,Bt,nu,work2(nv*nv*i),nv)
         enddo
         call mxm_f(work2,nv*nv,Ct,nu,v(1,ie),nv)
      enddo
      return
      end
c------------------------------------------   T  -----------------------
c----------------------------------------------------------------------
      subroutine outmat2_f(a,m,n,k,name)
      include 'SIZE'
      real*4 a(m,n)
      character*4 name

      n2 = min(n,8)
      write(6,2) nid,name,m,n,k
      do i=1,m
         write(6,1) nid,name,(a(i,j),j=1,n2)
      enddo
c   1 format(i3,1x,a4,16f6.2)
    1 format(i3,1x,a4,1p8e14.5)
    2 format(/,'Matrix: ',i3,1x,a4,3i8)
      return
      end
c-----------------------------------------------------------------------

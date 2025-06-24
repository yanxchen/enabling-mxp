c=======================================================================
c                  Math kernels in single precision
c=======================================================================

c-----------------------------------------------------------------------
c calculate inverse of an array b and assign to a
      subroutine invers2_f(a,b,n)
      REAL*4 A(1),B(1)
 
      DO 100 I=1,N
         A(I)=1./B(I)
 100  CONTINUE
      return
      END
c-----------------------------------------------------------------------
      subroutine invcol1_f(a,n)
      REAL*4 A(1)
 
      DO 100 I=1,N
         A(I)=1./A(I)
 100  CONTINUE
      return
      END
c-----------------------------------------------------------------------
      subroutine invcol2_f(a,b,n)
 
      REAL*4 A(1),B(1)
 
      DO 100 I=1,N
         A(I)=A(I)/B(I)
 100  CONTINUE
      return
      END
c-----------------------------------------------------------------------
      subroutine invcol3_f(a,b,c,n)
      REAL*4 A(1),B(1),C(1)
 
 
      DO 100 I=1,N
         A(I)=B(I)/C(I)
 100  CONTINUE
      return
      END
c-----------------------------------------------------------------------
      subroutine col4_f(a,b,c,d,n)
      REAL*4 A(1),B(1),C(1),D(1)
 
      DO 100 I=1,N
         A(I)=B(I)*C(I)*D(I)
  100 CONTINUE
      return
      END
c-----------------------------------------------------------------------
      subroutine Xaddcol3_f(a,b,c,n)
      REAL*4 A(1),B(1),C(1)
 
      DO 100 I=1,N
         A(I)=A(I)+B(I)*C(I)
  100 CONTINUE
      return
      END
c-----------------------------------------------------------------------
      subroutine addcol4_f(a,b,c,d,n)
      REAL*4 A(1),B(1),C(1),D(1)
 
      DO 100 I=1,N
         A(I)=A(I)+B(I)*C(I)*D(I)
  100 CONTINUE
      return
      END
c-----------------------------------------------------------------------
      subroutine ascol5_f(a,b,c,d,e,n)
      REAL*4 A(1),B(1),C(1),D(1),E(1)
 
      DO 100 I=1,N
         A(I) = B(I)*C(I)-D(I)*E(I)
 100  CONTINUE
      return
      END
c-----------------------------------------------------------------------
      subroutine sub2_f(a,b,n)
      REAL*4 A(1),B(1)
 
      DO 100 I=1,N
         A(I)=A(I)-B(I)
 100  CONTINUE
      return
      END
c-----------------------------------------------------------------------
      subroutine sub3_f(a,b,c,n)
      REAL*4 A(1),B(1),C(1)
 
      DO 100 I=1,N
         A(I)=B(I)-C(I)
 100  CONTINUE
      return
      END
c-----------------------------------------------------------------------
      subroutine subcol3_f(a,b,c,n)
      REAL*4 A(1),B(1),C(1)
 
 
      DO 100 I=1,N
         A(I)=A(I)-B(I)*C(I)
  100 CONTINUE
      return
      END
c-----------------------------------------------------------------------
      subroutine subcol4_f(a,b,c,d,n)
      REAL*4 A(1),B(1),C(1),D(1)
 
 
      DO 100 I=1,N
         A(I)=A(I)-B(I)*C(I)*D(I)
  100 CONTINUE
      return
      END
c-----------------------------------------------------------------------
      subroutine copy_f(a,b,n)
      real*4 a(1),b(1)

      do i=1,n
         a(i)=b(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine chsign_f(a,n)
      REAL*4 A(1)
 
      DO 100 I=1,N
         A(I) = -A(I)
 100  CONTINUE
      return
      END
 
c-----------------------------------------------------------------------
C multiply array with const
      subroutine cmult_f(a,const,n)
      REAL*4 A(1)
      REAL*4 const
 
      DO 100 I=1,N
         A(I)=A(I)*CONST
 100  CONTINUE
      return
      END
c-----------------------------------------------------------------------
      subroutine cadd_f(a,const,n)
      REAL*4 A(1)
      REAL*4 const
 
      DO 100 I=1,N
         A(I)=A(I)+CONST
 100  CONTINUE
      return
      END
c-----------------------------------------------------------------------
      subroutine cadd2_f(a,b,const,n)
      REAL*4 A(1),B(1)
      REAL*4 const
 
      DO 100 I=1,N
         A(I)=B(I)+CONST
 100  CONTINUE
      return
      END
c-----------------------------------------------------------------------
      real*4 function vlmin_f(vec,n)
      REAL*4 VEC(1)
      REAL*4 TMIN
      TMIN = 99.0E20
 
      DO 100 I=1,N
         TMIN = MIN(TMIN,VEC(I))
 100  CONTINUE
      vlmin_f = TMIN
      return
      END
c-----------------------------------------------------------------------
      real*4 function vlmax_f(vec,n)
      REAL*4 VEC(1)
      REAL*4 TMIN
      TMAX =-99.0E20
      do i=1,n
         TMAX = MAX(TMAX,VEC(I))
      enddo
      vlmax_f = TMAX
      return
      END
c-----------------------------------------------------------------------
      real*4 function vlamax_f(vec,n)
      REAL*4 VEC(1)
      REAL*4 TAMAX
      TAMAX = 0.0
 
      DO 100 I=1,N
         TAMAX = MAX(TAMAX,ABS(VEC(I)))
 100  CONTINUE
      vlamax_f = TAMAX
      return
      END
c-----------------------------------------------------------------------
      real*4 function vlsum_f(vec,n)
      REAL*4 VEC(1)
      REAL*4 SUM
      SUM = 0.
 
      DO 100 I=1,N
         SUM=SUM+VEC(I)
 100  CONTINUE
      vlsum_f = SUM
      return
      END
c-----------------------------------------------------------------------
      subroutine col2_f(a,b,n)
      real*4 a(1),b(1)

!xbm* unroll (10)
      do i=1,n
         a(i)=a(i)*b(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine col2c_f(a,b,c,n)
      real*4 a(1),b(1),c

      do i=1,n
         a(i)=a(i)*b(i)*c
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine col3_f(a,b,c,n)
      real*4 a(1),b(1),c(1)

!xbm* unroll (10)
      do i=1,n
         a(i)=b(i)*c(i)
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine add2_f(a,b,n)
      real*4 a(1),b(1)

!xbm* unroll (10)
      do i=1,n
         a(i)=a(i)+b(i)
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine add3_f(a,b,c,n)
      real*4 a(1),b(1),c(1)

!xbm* unroll (10)
      do i=1,n
         a(i)=b(i)+c(i)
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine addcol3_f(a,b,c,n)
      real*4 a(1),b(1),c(1)

!xbm* unroll (10)
      do i=1,n
         a(i)=a(i)+b(i)*c(i)
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine add2s1_f(a,b,c1,n)
      real*4 A(1),B(1)
      real*4 C1
 
      DO 100 I=1,N
        A(I)=C1*A(I)+B(I)
  100 CONTINUE
      return
      END

c-----------------------------------------------------------------------
      subroutine add2s2_f(a,b,c1,n)
      real*4 a(1),b(1)
      real*4 c1

      DO 100 I=1,N
        A(I)=A(I)+C1*B(I)
  100 CONTINUE
      return
      END
 
c-----------------------------------------------------------------------
      subroutine add3s2_f(a,b,c,c1,c2,n)
      real*4 a(1),b(1),c(1)
      real*4 c1,c2
      DO 100 I=1,N
        A(I)=C1*B(I)+C2*C(I)
  100 CONTINUE
      return
      END
 
c-----------------------------------------------------------------------
      subroutine add4_f(a,b,c,d,n)
      REAL*4 A(1),B(1),C(1),D(1)
 
      DO 100 I=1,N
         A(I)=B(I)+C(I)+D(I)
 100  CONTINUE
      return
      END
c-----------------------------------------------------------------------
      real*4 function vlsc2_f(x,y,n)
      REAL*4 X(1),Y(1)
      REAL*4 s
      s = 0.
      do i=1,n
         s = s + x(i)*y(i)
      enddo
      vlsc2_f=s
      return
      end
c-----------------------------------------------------------------------
      real*4 function vlsc21_f(x,y,n)
      real*4 x(1),y(1)
      real*4 s
      s = 0.
      do i=1,n
         s = s + x(i)*x(i)*y(i)
      enddo
      vlsc21_f=s
      return
      end

C----------------------------------------------------------------------------
C
C     Vector reduction routines which require communication 
C     on a parallel machine. These routines must be substituted with
C     appropriate routines which take into account the specific architecture.
C
C----------------------------------------------------------------------------
      function glsc3_f(a,b,mult,n)
C
C     Perform inner-product in single precision
C
      real*4 a(1),b(1),mult(1)
      real*4 tmp,work(1)

      tmp = 0.0
      do 10 i=1,n
         tmp = tmp + a(i)*b(i)*mult(i)
 10   continue
      call gop_f(tmp,work,'+  ',1)
      glsc3_f = tmp
      return
      end
c-----------------------------------------------------------------------
      function glsc3_eftdot_mp(a,b,mult,n)
C
C     Perform inner-product with EFT in mixed precision
C
      real*4 a(1),b(1),mult(1)
      real*4 tmp(1),work(1)
      real*4 acc

      do 10 i=1,n
         tmp(i) = b(i)*mult(i)
   10 continue
      call dot2_s(a,tmp,n,acc)
      call gop_f(acc,work,'+  ',1)
      glsc3_eftdot_mp = acc
      return
      end
c-----------------------------------------------------------------------
      function glsc2_f(x,y,n)
C
C     Perform inner-product in single precision
C
      real*4 x(1), y(1)
      real*4 tmp,work(1)
 
      tmp=0.0
      do 10 i=1,n
         tmp = tmp+ x(i)*y(i)
   10 continue
      CALL GOP(TMP,WORK,'+  ',1)
      glsc2_f = TMP
      return
      END
c-----------------------------------------------------------------------
      function glsc23_f(x,y,z,n)
c
C     Perform inner-product  x*x*y*z
c
      real*4 x(1), y(1),z(1)
      real*4 tmp,work(1)
      real*4 ds
      ds = 0.0
      do 10 i=1,n
         ds=ds+x(i)*x(i)*y(i)*z(i)
   10 continue
      tmp=ds
      call gop(tmp,work,'+  ',1)
      glsc23_f = tmp
      return
      end
c-----------------------------------------------------------------------
      function glsum_f(x,n)
      DIMENSION X(1)
      DIMENSION TMP(1),WORK(1)
      REAL*4 TSUM
      TSUM = 0.
      DO 100 I=1,N
         TSUM = TSUM+X(I)
 100  CONTINUE
      TMP(1)=TSUM
      CALL GOP(TMP,WORK,'+  ',1)
      glsum_f = TMP(1)
      return
      END
c-----------------------------------------------------------------------
      real*4 function glamax_f(a,n)
      REAL*4 A(1)
      DIMENSION TMP(1),WORK(1)
      REAL*4 TMAX
      TMAX = 0.0
      DO 100 I=1,N
         TMAX = MAX(TMAX,ABS(A(I)))
 100  CONTINUE
      TMP(1)=TMAX
      CALL GOP(TMP,WORK,'M  ',1)
      glamax_f=ABS(TMP(1))
      return
      END
c-----------------------------------------------------------------------
      real*4 function glamin_f(a,n)
      real*4 a(1)
      dimension tmp(1),work(1)
      real*4 tmin
      tmin = 9.e28
      do 100 i=1,n
         tmin = min(tmin,abs(a(i)))
 100  continue
      tmp(1)=tmin
      call gop(tmp,work,'m  ',1)
      glamin_f=abs(tmp(1))
      return
      end
C-----------------------------------------------------------------------
      function glmax_f(a,n)
      REAL*4 A(1)
      DIMENSION TMP(1),WORK(1)
      REAL*4 TMAX
      TMAX=-99.0e20
      DO 100 I=1,N
         TMAX=MAX(TMAX,A(I))
  100 CONTINUE
      TMP(1)=TMAX
      CALL GOP(TMP,WORK,'M  ',1)
      glmax_f=TMP(1)
      return
      END
c-----------------------------------------------------------------------
      function glmin_f(a,n)
      REAL*4 A(1)
      DIMENSION TMP(1),WORK(1)
      REAL*4 TMIN
      TMIN=99.0e20
      DO 100 I=1,N
         TMIN=MIN(TMIN,A(I))
  100 CONTINUE
      TMP(1)=TMIN
      CALL GOP(TMP,WORK,'m  ',1)
      glmin_f = TMP(1)
      return
      END
c-----------------------------------------------------------------------
      function fmdian_f(a,n,ifok)
C     find the Median of the (global) set A
      include 'SIZE'
      DIMENSION A(1)
      DIMENSION WORK1(5),WORK2(5)
      DIMENSION GUES(100)
      LOGICAL IFOK
      REAL*4 AMP, AFAC
      AMP  =1.5
      AFAC =1.5
      GMIN =GLMIN(A,N)
      GMAX =GLMAX(A,N)
      GMIN0=GLMIN(A,N)
      GMAX0=GLMAX(A,N)
      GUESS=(GMAX+GMIN)/2.0
      EPS  =(GMAX-GMIN)
      IF (EPS.EQ.0.0) THEN
         fmdian_f=GMAX
         return
      ENDIF
      WORK1(1)=N
      CALL GOP(WORK1,WORK2,'+  ',1)
      NTOT=WORK1(1)
      N2 = (NTOT+1)/2
      IF (.NOT.IFOK) THEN
        WRITE(6,8) NID,N,(A(I),I=1,N)
        WRITE(6,9) NID,NTOT,N2,N,GMIN,GMAX
    8   FORMAT(I5,'N,A:',I5,10(6F10.5,/)) 
    9   FORMAT(I5,'mnx:',3I6,2F10.5)
      ENDIF
C
C     This is the trial loop
C
      ITRY=-1
   10 CONTINUE
      ITRY=ITRY+1
      II=ITRY+1
      IF (II.LE.100) GUES(II)=GUESS
C     error check for infinite loop
      IF (ITRY.GT.2*NTOT) GOTO 9000
      CALL RZERO(WORK1,5)
      NLT=0
      NGT=0
      CLT=GMIN0
      CGT=GMAX0
      DO 100 I=1,N
         AA=A(I)
         IF (AA.NE.GUESS) THEN
            IF (AA.LT.GUESS) THEN
               NLT=NLT+1
C              CLT - closest value to GUESS Less Than GUESS
               IF (AA.GT.CLT) CLT=AA
            ENDIF
            IF (AA.GT.GUESS) THEN
               NGT=NGT+1
C              CGT - closest value to GUESS Greater Than GUESS
               IF (AA.LT.CGT) CGT=AA
            ENDIF
            DUM=1./(EPS+ABS(AA-GUESS))
            WORK1(1)=WORK1(1)+DUM
            WORK1(2)=WORK1(2)+DUM*AA
         ELSE
C           detected values equaling the guess.
            WORK1(5)=WORK1(5)+1.0
         ENDIF
  100 CONTINUE
C     Invoke vector reduction across processors:
      WORK2(1)=CLT
      CLT=GLMAX(WORK2,1)
      WORK2(1)=CGT
      CGT=GLMIN(WORK2,1)
      WORK1(3)=NLT
      WORK1(4)=NGT
      CALL GOP(WORK1,WORK2,'+  ',5)
      NLT=WORK1(3)
      NGT=WORK1(4)
      IF (.NOT.IFOK) THEN
         WRITE(6,101) NID,GUESS,CLT,CGT
         WRITE(6,102) NID,(WORK1(I),I=1,5)
  101    FORMAT(I5,'Glg:',3F12.5)
  102    FORMAT(I5,'WORK1:',5F12.5)
      ENDIF
C
C     Done?
C
      IF (NLT.GT.N2.OR.NGT.GT.N2) THEN
C        we're not done.....
         IF (NGT.GT.NLT) THEN
C           guess is too low
            GMIN=CGT
            G2=CGT+MAX(0.,WORK1(2)/WORK1(1)-GUESS)*AMP
            IF (G2.GT.GMAX) G2=0.5*(GUESS+GMAX)
            EPS=AFAC*ABS(G2-GUESS)
C           see that we move at least as far as the next closest value.
            GUESS=MAX(G2,CGT)
            GOTO 10
         ELSE IF (NLT.GT.NGT) THEN
C           guess is too high
            GMAX=CLT
            G2=CLT+MIN(0.,WORK1(2)/WORK1(1)-GUESS)*AMP
            IF (G2.LT.GMIN) G2=0.5*(GUESS+GMIN)
            EPS=AFAC*ABS(G2-GUESS)
C           see that we move at least as far as the next closest value.
            GUESS=MIN(G2,CLT)
            GOTO 10
         ENDIF
      ELSE
C
C        we're done....
         IF (WORK1(5).NE.0) THEN
C           the median is (usually) one of the values 
            fmdian_f=GUESS
            IF (WORK1(5).EQ.1.0) THEN
               IF (MOD(NTOT,2).EQ.0) THEN
                  IF (NGT.GT.NLT) THEN
                     fmdian_f=0.5*(GUESS+CGT)
                  ELSE
                     fmdian_f=0.5*(GUESS+CLT)
                  ENDIF
               ELSE
                  IF (NGT.EQ.NLT) THEN
                     fmdian_f=GUESS
                  ELSE IF(NGT.GT.NLT) THEN
                     fmdian_f=CGT
                  ELSE
                     fmdian_f=CLT
                  ENDIF
               ENDIF
            ENDIF
         ELSE
            IF (MOD(NTOT,2).EQ.0) THEN
               IF (NGT.EQ.NLT) THEN
                  fmdian_f=0.5*(CLT+CGT)
               ELSE IF(NGT.GT.NLT) THEN
                  fmdian_f=0.5*(GUESS+CGT)
               ELSE
                  fmdian_f=0.5*(GUESS+CLT)
               ENDIF
            ELSE
               IF (NGT.EQ.NLT) THEN
                  fmdian_f=GUESS
               ELSE IF(NGT.GT.NLT) THEN
                  fmdian_f=CGT
               ELSE
                  fmdian_f=CLT
               ENDIF
           ENDIF
         ENDIF
 
      ENDIF
       IF (.NOT.IFOK) WRITE(6,*) NID,'FMDIAN2',FMDIAN,(A(I),I=1,N)
      return
C
C     Error handling
C
 9000 CONTINUE
      WRITE(6,11) NTOT,GMIN0,GMAX0,GUESS
   11 FORMAT('ABORTING IN FMDIAN: N,AMIN,AMAX:',I6,3G14.6)
      DO 13 I1=1,N,5
        IN=I1+5 
        IN=MIN(IN,N)
        WRITE(6,12) NID,(A(I),I=I1,IN)
   12   FORMAT(I4,' FMA:',5G14.6)
   13 CONTINUE
      DO 15 I1=1,ITRY,5
        IN=I1+5
        IN=MIN(IN,ITRY)
        WRITE(6,14) NID,(GUES(I),I=I1,IN)
   14   FORMAT(I4,' FMG:',5G14.6)
   15 CONTINUE
      call exitt
      END

C========================================================================
C     Double precision matrix and vector routines
C========================================================================

c-----------------------------------------------------------------------
      subroutine dcadd_f(a,const,n)
      real*4 A(1),CONST
 
      DO 100 I=1,N
         A(I)=A(I)+CONST
 100  CONTINUE
      return
      END
c-----------------------------------------------------------------------
      subroutine dsub2_f(a,b,n)
      real*4 A(1), B(1)
 
      DO 100 I=1,N
         A(I)=A(I)-B(I)
 100  CONTINUE
      return
      END
 
c-----------------------------------------------------------------------
      subroutine dadd2_f(a,b,n)
      real*4 A(1), B(1)
 
      DO 100 I=1,N
         A(I)=A(I)+B(I)
 100  CONTINUE
      return
      END

c-----------------------------------------------------------------------
      subroutine drcopy_f(r,d,N)
      real*4      d(1)
      dimension r(1)
      do 10 i=1,n
         r(i)=d(i)
   10 continue
      return
      end
      subroutine sorts_f(xout,xin,work,n)
      real*4 xout(1),xin(1),work(1)
      call copy(xout,xin,n)
      call sort(xout,work,n)
      return
      end
C
c-----------------------------------------------------------------------
      subroutine sort_f(a,ind,n)
C
C     Use Heap Sort (p 231 Num. Rec., 1st Ed.)
C
      real*4 a(1),aa
      integer ind(1)
 
      dO 10 j=1,n
         ind(j)=j
   10 continue
 
      if (n.le.1) return
      L=n/2+1
      ir=n
  100 continue
         if (l.gt.1) then
            l=l-1
            aa  = a  (l)
            ii  = ind(l)
         else
                 aa =   a(ir)
                 ii = ind(ir)
              a(ir) =   a( 1)
            ind(ir) = ind( 1)
            ir=ir-1
            if (ir.eq.1) then
                 a(1) = aa
               ind(1) = ii
               return
            endif
         endif
         i=l
         j=l+l
  200    continue
         if (j.le.ir) then
            if (j.lt.ir) then
               if ( a(j).lt.a(j+1) ) j=j+1
            endif
            if (aa.lt.a(j)) then
                 a(i) = a(j)
               ind(i) = ind(j)
               i=j
               j=j+j
            else
               j=ir+1
            endif
         GOTO 200
         endif
           a(i) = aa
         ind(i) = ii
      GOTO 100
      end
c-----------------------------------------------------------------------
      subroutine swap_ip_f(x,p,n)
      real*4    x(1),xstart
      integer p(1)
c
c     In-place permutation: x' = x(p)
c
      do k=1,n
         if (p(k).gt.0) then   ! not swapped
            xstart     = x(k)
            loop_start = k
            last       = k
            do j=k,n
               next    = p(last)
               if (next.lt.0) then
                  write(6,*) 'Hey! swap_ip problem.',j,k,n,next
                  call exitt
               elseif (next.eq.loop_start) then
                  x(last) = xstart
                  p(last) = -p(last)
                  goto 10
               else
                  x(last) = x(next)
                  p(last) = -p(last)
                  last    = next
               endif
            enddo
   10       continue
         endif
      enddo
 
      do k=1,n
         p(k) = -p(k)
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine swapt_ip_f(x,p,n)
      real*4    x(1),t1,t2
      integer p(1)
c
c     In-place permutation: x'(p) = x
c

      do k=1,n
         if (p(k).gt.0) then   ! not swapped
            loop_start = k
            next       = p(loop_start)
            t1         = x(loop_start)
            do j=1,n
               if (next.lt.0) then
                  write(6,*) 'Hey! swapt_ip problem.',j,k,n,next
                  call exitt
               elseif (next.eq.loop_start) then
                  x(next) = t1
                  p(next) = -p(next)
                  goto 10
               else
                  t2      =  x(next)
                  x(next) =  t1
                  t1      =  t2
                  nextp   =  p(next)
                  p(next) = -p(next)
                  next    =  nextp
               endif
            enddo
   10       continue
         endif
      enddo
 
      do k=1,n
         p(k) = -p(k)
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine glvadd_f(x,w,n)
      real*4 x(1),w(1)
      call gop(x,w,'+  ',1)
      return
      end
c-----------------------------------------------------------------------
      subroutine add3s12_f(x,y,z,c1,c2,n)
      real*4 x(1),y(1),z(1),c1,c2
      do i=1,n
         x(i) = c1*y(i)+c2*z(i)
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine admcol3_f(a,b,c,d,n)
      REAL*4 A(1),B(1),C(1),D
C
      DO 100 I=1,N
         A(I)=A(I)+B(I)*C(I)*D
  100 CONTINUE
      return
      END
c-----------------------------------------------------------------------
      subroutine add2col2_f(a,b,c,n)
      real*4 a(1),b(1),c(1)
 
      do i=1,n
         a(i) = a(i) + b(i)*c(i)
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine add2sxy_f(x,a,y,b,n)
      real*4 x(1),y(1)
 
      do i=1,n
         x(i) = a*x(i) + b*y(i)
      enddo
 
      return
      end
c-----------------------------------------------------------------------
      subroutine col2s2_f(x,y,s,n)
      real*4 x(n),y(n)
 
      do i=1,n
         x(i)=s*x(i)*y(i)
      enddo
 
      return
      end
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
c EFT-dot
c-----------------------------------------------------------------------

      SUBROUTINE TWOSUM_S(A, B, X, Y)
      REAL*4 A, B, X, Y, Z
      
      Z = 0.0
      X = A + B
      Z = X - A
      Y = (A - (X - Z)) + (B - Z)
      RETURN
      END

      SUBROUTINE SPLIT_S(A, X, Y)
      REAL*4 A, X, Y, FACTOR, C
      
      FACTOR = 4097.0
      C = 0.0
      C = FACTOR * A
      X = C - (C - A)
      Y = A - X
      RETURN
      END

      SUBROUTINE TWOPRD_S(A, B, X, Y)
      REAL*4 A, B, X, Y, A1, A2, B1, B2
      
      A1 = 0.0
      A2 = 0.0
      B1 = 0.0
      B2 = 0.0
      X = A * B
      CALL SPLIT_S(A, A1, A2)
      CALL SPLIT_S(B, B1, B2)
      Y = A2 * B2 - (((X - A1 * B1) - A2 * B1) - A1 * B2)
      RETURN
      END

      SUBROUTINE DOT2_S(X, Y, N, RES)
      INTEGER N, I
      REAL*4 X(N), Y(N), RES
      REAL*4 P, S, H, R, Q, TMP
      
      P = 0.0
      S = 0.0
      H = 0.0
      R = 0.0
      Q = 0.0
      RES = 0.0
      
      CALL TWOPRD_S(X(1), Y(1), P, S)
      DO 10 I = 2, N
         CALL TWOPRD_S(X(I), Y(I), H, R)
         TMP = P
         CALL TWOSUM_S(TMP, H, P, Q)
         S = S + (Q + R)
 10   CONTINUE
      RES = P + S
      RETURN
      END


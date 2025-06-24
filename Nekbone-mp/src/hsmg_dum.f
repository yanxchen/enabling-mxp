c  Dummy file for hsmg 
c-----------------------------------------------------------------------
      subroutine h1mg_setup()

      return
      end
c----------------------------------------------------------------------
c----------------------------------------------------------------------
      subroutine h1mg_solve(z,rhs,n)  !  Solve preconditioner: Mz=rhs in double
      real z(n),rhs(n)
  
      call copy(z,rhs,n)

      return
      end
c-----------------------------------------------------------------------
      subroutine h1mg_solve_f(z,rhs,n)  !  Solve preconditioner: Mz=rhs in single
      real*4 z(n),rhs(n)
  
      call copy_f(z,rhs,n)

      return
      end
c-----------------------------------------------------------------------
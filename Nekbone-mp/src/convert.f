c-----------------------------------------------------------------------
      SUBROUTINE convert_dp_to_sp(A, B, N)
C     Converts a double precision array A to a single precision array B

      INTEGER N, I
      REAL A(N)
      REAL*4 B(N)

      DO 10 I = 1, N
         B(I) = A(I)
  10  CONTINUE

      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE convert_sp_to_dp(A, B, N)
C     Converts a single precision array A to a double precision array B

      INTEGER N, I
      REAL*4 A(N)
      REAL B(N)

      DO 10 I = 1, N
         B(I) = A(I)
  10  CONTINUE

      RETURN
      END
c-----------------------------------------------------------------------

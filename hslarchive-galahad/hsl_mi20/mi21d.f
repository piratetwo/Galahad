! COPYRIGHT (c) 1995 Council for the Central Laboratory
!                    of the Research Councils
! Original date 23 March 2001
! Threadsafe version of IM01
!
! Version 1.3.0
! See ChangeLog for version history
!
      SUBROUTINE MI21ID( ICNTL, CNTL, ISAVE, RSAVE )
      DOUBLE PRECISION CNTL( 5 )
      INTEGER ICNTL( 8 )
      INTEGER ISAVE(10)
      DOUBLE PRECISION RSAVE(6)
C
C  If A is symmetric, positive definite, MI21 solves the linear system
C
C      A x = b
C
C  using the
C
C  ===================
C  Conjugate Gradients
C  ===================
C
C  iterative method  optionally using preconditioning. If P
C  is a symmetric positive definite preconditioner and
C  P = P_L * P_L(transpose), the algorithm actually solves
C
C                             ^       ^                      ^
C        P_L A P_L(transpose) x = P_L b,  x = P_L(transpose) x
C
C  Each time a matrix-vector product with A or P is required,
C  control is passed back to the user.
C
C  MI21I/ID is the initialisation routine for MI21A/AD and should
C  be called once prior to calls to MI21A/AD.
C
C  Argument list.
C
C  ICNTL   (output) INTEGER control array, dimension 8.
C          ICNTL(1) is the stream number for error messages.
C             On exit, ICNTL(2) = 6.
C          ICNTL(2) is the stream number for warning messages.
C             On exit, ICNTL(2) = 6.
C          ICNTL(3) indicates whether the user wishes to use a
C             preconditioner. If ICNTL(3) is nonzero, preconditioning.
C             On exit, ICNTL(3) = 0
C          ICNTL(4) indicates whether
C             the convergence test offered by MI21A/AD is to be used.
C             If ICNTL(4) = 0, the computed solution x is accepted
C             if ||Ax - b|| < max( ||A x_0 - b || * CNTL(1), CNTL(2) ).
C             Otherwise, the user most perform his/her
C             own test for convergence when IACT = 1 is returned.
C             On exit, ICNTL(4) = 0
C          ICNTL(5) indicates whether the user wishes to supply an
C             initial guess for the solution vector x.
C             If ICNTL(5) = 0, the user does not wish to supply
C             an initial guess and x = (0,0,...,0) will be used
C             as the initial guess. Otherwise, the user
C             must supply an intial guess on the first call to
C             MI21A/AD. On exit, ICNTL(5) = 0.
C          ICNTL(6) determines the maximum number of iterations
C             allowed. It has default value -1 and, in this case,
C             the maximum number will be N. If the user does
C             not want the maximum number to be N, ICNTL(6) should
C             be set to the maximum  number the user wishes
C             to allow.
C          ICNTRL(7) determines whether the normalized curvature is used
C             when testing for breakdown.  The
C             normalized curvature is only used if ICNTL(7)=1.
C             The default value is 0.
C          ICNTRL(8) is spare and set to zero.
C
C  CNTL    (output)  REAL (DOUBLE PRECISION) control array, dimension 5.
C          CNTL(1) is a convergence tolerance.
C             On exit, set to the square root of machine precision.
C             If ICNTL(4) is nonzero, CNTL(1) not accessed by MI21A/AD.
C          CNTL(2) is a convergence tolerance. On exit, set to zero.
C             If ICNTL(4) is nonzero, CNTL(2) not accessed by MI21A/AD.
C          CNTL(3) is tolerance used to check whether the algorithm
C             has broken down. On exit, set to the machine precision
C          CNTRL(4) and CNTRL(5) are spare and set to zero.
C
C  ISAVE   (output) INTEGER ARRAY, length 10, used to hold the
C          routine's persistent integer data.  The array contents
C          must not be altered by the user.
C
C  RSAVE   (output) DOUBLE PRECISION ARRAY, length 6, used to hold the
C          routine's persistent real data.  The array contents
C          must not be altered by the user.
C
C  Local variables
C
      INTEGER I
C
C     Intrinsic Functions
C
      INTRINSIC SQRT
      ICNTL( 1 ) = 6
      ICNTL( 2 ) = 6
      ICNTL( 3 ) = 0
      ICNTL( 4 ) = 0
      ICNTL( 5 ) = 0
      ICNTL( 6 ) = - 1

      ICNTL(7) = 0
      ICNTL(8) = 0

      CNTL( 1 ) = SQRT( EPSILON(CNTL) )
      CNTL( 2 ) = 0.0
      CNTL( 3 ) = - 1.0

      CNTL(4) = 0.0
      CNTL(5) = 0.0

C  Initialize persistent data to avoid undefined assignment in MI21AD
      DO 10 I = 1, 10
      ISAVE(I) = 0
   10 CONTINUE
      DO 20 I = 1, 6
      RSAVE(I) = 0.0
   20 CONTINUE
      RETURN
C
C  End of MI21I/MI21ID
C
      END

      SUBROUTINE MI21AD( IACT, N, W, LDW, LOCY, LOCZ, RESID,
     *                   ICNTL, CNTL, INFO, ISAVE, RSAVE )
      DOUBLE PRECISION RESID
      INTEGER          IACT, LDW, LOCY, LOCZ, N
      DOUBLE PRECISION CNTL( 5 ), W( LDW, 4 )
      INTEGER          ICNTL( 8 ), INFO( 4 )
      INTEGER ISAVE(10)
      DOUBLE PRECISION RSAVE(6)
C
C  Argument list.
C
C  IACT    (input) INTEGER.
C          IACT must be set to 0 prior to first call.
C          On each exit, IACT indicates the action required by
C          the user. Possible values of IACT and the action
C          required are as follows:
C  -1      Fatal error (see INFO(1)). Terminate computation.
C   1      If ICNTL(4) = 0 (the default), convergence has been
C          achieved and the user should terminate the computation.
C          If ICNTL(4) is nonzero, the user should test the norm of
C          the residual in RESID for convergence and recall MI21A/AD
C          if convergence has not been achieved.
C   2      The user must perform the matrix-vector product
C          y := Az,
C          and recall MI21A/AD. The vectors  y and z are held in the
C          first N entries of columns LOCZ and LOCZ of array W,
C          respectively.
C   3      The user must perform the preconditioning operation
C          y := Pz, where P is the preconditioner,
C          and recall MI21A/AD. The vectors y and z
C          are  held in the first N entries of
C          columns LOCZ and LOCZ of array W, respectively.
C
C  N       (input) INTEGER.
C          On entry, the dimension of the matrix.
C          Unchanged on exit.
C
C  W       (right-hand side, solution and workspace, input/output)
C          REAL (DOUBLE PRECISION) array, dimension (LDW,4).
C          Prior to the first call, the first N entries of column
C          1 must be set to hold the right-hand side vector b and,
C          if ICNTL(5) is nonzero, the first N entries of column 2
C          must be set to the initial estimate of the solution vector
C          x.  On exit with IACT = 1, the first N entries of column 1
C          hold the current residual vector r = b - Ax, and the
C          current estimate of the solution x is held in the
C          first N entries of column 2.  On exit
C          with IACT > 1, the user is required to perform a
C          computation with columns LOCZ and LOCZ of W
C          (see argument IACT).  The remaining columns of
C          W must not be altered by the user between calls to
C          MI21A/AD.
C
C  LDW     (input) INTEGER
C          The leading dimension of the array W. LDW >= max(1,N).
C
C  LOCY, LOCZ (output) INTEGER
C          On exit with IACT > 1, LOCY and LOCZ define the columns
C          of the array W which hold y and z.
C          (see IACT).
C
C  RESID   (output) REAL (DOUBLE PRECISION)
C          On exit with IACT = 1,
C          RESID holds ||b - Ax||, where x is the
C          iterated solution.
C          If ICNTL(4) is nonzero, on exit with IACT = 1,
C          the user should carry out his/her
C          own test for convergnce at this point.
C
C  ICNTL   (input) INTEGER control array of dimension 8.
C          ICNTL may be initalised by calling MI21I/ID.
C          See MI21I/ID for details.
C
C  CNTL    (input) REAL (DOUBLE PRECISION) control array of dimension 5.
C          CNTL may be initalised by calling MI21I/ID.
C          See MI21I/ID for details.
C
C  INFO    (output) INTEGER ARRAY , length 4.
C          If INFO(1) = 0 on exit, no errors or warnings issued.
C          If INFO(1) > 0 on exit, warning to user.
C          INFO(1) = 1, value of CNTL(1) is out-of-range (u,1.0)
C          The default sqrt(u) is used, u=machine precision.
C          If INFO(1) < 0 on exit, illegal input parameter,
C               or breakdown occured during iteration.
C
C                Error input data:
C
C                   -1: matrix dimension N < 1
C                   -2: LDW < N
C
C                BREAKDOWN: If curvature smaller than CNTL(3) is
C                   encountered, the program will terminate.
C
C                   -3: CURV < CNTL(3):
C          On each exit, INFO(2) holds the number of
C          iterations performed so far.
C
C  ISAVE   (input/output) INTEGER ARRAY, length 10, used to hold the
C          routine's persistent integer data.  The array contents
C          must not be altered by the user.
C
C  RSAVE   (input/output) DOUBLE PRECISION ARRAY, length 6, used to
C          hold the routine's persistent real data. The array contents
C          must not be altered by the user.
C
C  BLAS calls:    DAXPY, DCOPY, DDOT, DNRM2, DSCAL
C
      DOUBLE PRECISION ONE, ZERO
      PARAMETER ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C
C  Local variables
C
      INTEGER          I, IPOS, ITMAX, B, P, Q, R, X, Z
      DOUBLE PRECISION ALPHA, BETA, BNRM2, CURV, RHO, RHO1
C
C  External Functions
C
      DOUBLE PRECISION DDOT, DNRM2
      EXTERNAL         DDOT, DNRM2
C
C  Intrinsic Functions
C
      INTRINSIC        SQRT
C
C  External Subroutines
C
      EXTERNAL         DAXPY, DCOPY, DSCAL
C
C  Restore persistent data
C
      IPOS  = ISAVE(1)
      ITMAX = ISAVE(2)
      B     = ISAVE(3)
      P     = ISAVE(4)
      Q     = ISAVE(5)
      R     = ISAVE(6)
      X     = ISAVE(7)
      Z     = ISAVE(8)

      BNRM2 = RSAVE(1)
      RHO   = RSAVE(2)
      RHO1  = RSAVE(3)
      CURV  = RSAVE(4)
C
C  Jump to appropriate place in code
C
      IF ( IACT .NE. 0 ) THEN
C
C  Immediate return if error on a previous call
C
         IF ( IACT .LT. 0 ) GO TO 1000
C
C  Immediate return if convergence already achieved
C
         IF ( IACT .EQ. 1 .AND. ICNTL( 4 ) .EQ. 0 ) GO TO 1000
         IF ( IACT .EQ. 1 .AND. BNRM2 .EQ. ZERO ) GO TO 1000
C
C  Branch
C
         GO TO ( 30, 50, 60, 70 ), IPOS
      END IF
C
C  Initial call
C
      INFO( 1 ) = 0
C
C  Test the input parameters
C
      IF ( N .LT. 1 ) THEN
         INFO( 1 ) = - 1
      ELSE IF ( LDW .LT. N ) THEN
         INFO( 1 ) = - 2
      END IF
      IF ( INFO( 1 ) .LT. 0 ) THEN
         IACT = - 1
         IF ( ICNTL( 1 ) .GT. 0 ) WRITE( ICNTL( 1 ), 2000 ) INFO( 1 )
         GO TO 1000
      END IF
C
C  Alias workspace columns
C
      B = 1
      R = B
      X = 2
      P = 3
      Q = 4
      IF ( ICNTL( 3 ) .NE. 0 ) THEN
         Z = Q
      ELSE
         Z = R
      END IF
C
C  Set INFO(2) and ITMAX
C
      INFO( 2 ) = 0
      IF ( ICNTL( 6 ) .GT. 0 ) THEN
         ITMAX = ICNTL( 6 )
      ELSE
         ITMAX = N
      END IF
C
C  Set CNTL(3)
C
      IF ( CNTL( 3 ) .LE. ZERO ) CNTL( 3 ) = N * EPSILON(CNTL)
C
C  Compute ||b||
C
      BNRM2 = DNRM2( N, W( 1, B ), 1 )
C
C  Immediate return if ||b|| = 0
C
      IF ( BNRM2. EQ. ZERO ) THEN
         IACT = 1
         DO 10 I = 1, N
            W( I, X ) = ZERO
            W( I, B ) = ZERO
   10    CONTINUE
         RESID = ZERO
         GO TO 1000
      END IF
C
C  Check value of CNTL(1)
C
      IF ( ICNTL( 4 ) .EQ. 0 ) THEN
         IF ( CNTL( 1 ).LT.EPSILON(CNTL) .OR. CNTL( 1 ).GT.ONE ) THEN
            INFO( 1 ) = 1
            IF (ICNTL( 2 ) .GT. 0 ) THEN
               WRITE( ICNTL( 2 ), 2010 ) INFO( 1 )
               WRITE( ICNTL( 2 ), 2020 )
            END IF
            CNTL( 1 ) = SQRT( EPSILON(CNTL) )
         END IF
      END IF
C
C  Compute initial residual
C
C  If the user has not supplied an initial guess, set X = 0
C  as the initial guess
C
      IF ( ICNTL( 5 ) .EQ. 0 ) THEN
         DO 20 I = 1, N
            W( I, X ) = ZERO
   20    CONTINUE
         GO TO 40
      ELSE
C
C  Initial guess supplied by user
C  If initial guess for solution is x = 0 no action is required (r = b)
C
         IF ( DNRM2( N, W( 1, X ), 1 ) .EQ. ZERO ) GO TO 40
C
C  Otherwise, return to user to compute Ax.
C  Column P can be used temporarily to hold Ax
C
         IPOS = 1
         IACT = 2
         LOCY = P
         LOCZ = X
         GO TO 1000
      END IF
C
C  Compute r = b - Ax
C
   30 CONTINUE
      CALL DAXPY( N, - ONE, W( 1, P ), 1, W( 1, R ), 1 )
C
C  Compute ||r||
C
      BNRM2 = DNRM2( N, W( 1, R ), 1 )
C
C  Main iteration loop
C
   40 CONTINUE
C
C  Update iteration count
C
      INFO( 2 ) = INFO( 2 ) + 1
C
C  Check maximum number of iterations has not been exceeded
C
      IF ( INFO( 2 ) .GT. ITMAX ) THEN
         INFO( 1 ) = - 4
         IACT = - 1
         IF (ICNTL( 1 ) .GT. 0 ) THEN
            WRITE( ICNTL( 1 ), 2000 ) INFO( 1 )
            WRITE( ICNTL( 1 ), 2030 ) ITMAX
         END IF
         GO TO 1000
      END IF
C
C  Return to user to obtain preconditioner Z = P^-1 R
C
      IF ( ICNTL( 3 ) .NE. 0 ) THEN
         IPOS = 2
         IACT = 3
         LOCY = Z
         LOCZ = R
         GO TO 1000
      END IF
   50 CONTINUE
C
C  Compute the inner product R^T Z
C
      RHO = DDOT( N, W( 1, R ), 1, W( 1, Z ), 1 )
C
C  Compute search direction vector P
C
      IF ( INFO( 2 ) .EQ. 1 ) THEN
C
C  Special case: first iteration
C
         CALL DCOPY( N, W( 1, Z ), 1, W( 1, P ), 1 )
      ELSE
         BETA = RHO / RHO1
C
C  later iterations
C
         CALL DSCAL( N, BETA, W( 1, P ), 1 )
         CALL DAXPY( N, ONE, W( 1, Z ), 1, W( 1, P ), 1 )
      END IF
C
C  Return to user for matrix-vector product Q = A P
C
      IPOS = 3
      IACT = 2
      LOCY = Q
      LOCZ = P
      GO TO 1000
   60 CONTINUE
C
C  Obtain the curvature along P
C
      CURV = DDOT( N, W( 1, P ), 1, W( 1, Q ), 1 )
C
C  If the curvature is negative, A is indefinite
C
      IF (ICNTL(7).EQ.1) THEN
        IF ( CURV .LT. CNTL( 3 )*DDOT( N, W( 1, P ), 1, W( 1, P ), 1 ))
     +     INFO( 1 ) = - 3
      ELSE IF ( CURV .LT. CNTL( 3 ) ) THEN
         INFO( 1 ) = - 3
      END IF
      IF (INFO(1).EQ.-3) THEN
         IACT = - 1
         IF ( ICNTL( 1 ) .GT. 0 ) WRITE( ICNTL( 1 ), 2000 ) INFO( 1 )
         GO TO 1000
      END IF
C
C  Compute the stepsize
C
      ALPHA = RHO / CURV
C
C  Update the estimate of the solution
C
      CALL DAXPY( N, ALPHA, W( 1, P ), 1, W( 1, X ), 1 )
C
C  Update the residual
C
      CALL DAXPY( N, - ALPHA, W( 1, Q ), 1, W( 1, R ), 1 )
C
C  Iteration complete. Check convergence
C
      RESID = DNRM2( N, W( 1, R ), 1 )
      IPOS = 4
      IF ( ICNTL( 4 ) .NE. 0 ) THEN
C
C  Return the residual to the user for convergence testing
C
         IACT = 1
         GO TO 1000
      ELSE
C
C  Test the scaled residual for convergence
C
         IF ( RESID .LE. MAX( BNRM2 * CNTL( 1 ), CNTL( 2 ) ) ) THEN
C
C  Convergence achieved
C
            IACT = 1
            GO TO 1000
         END IF
      END IF
   70 CONTINUE
      RHO1 = RHO
C
C  Next iteration
C
      GO TO 40
C
C  Save persistent data and return to caller
C
 1000 CONTINUE
      ISAVE(1) = IPOS
      ISAVE(2) = ITMAX
      ISAVE(3) = B
      ISAVE(4) = P
      ISAVE(5) = Q
      ISAVE(6) = R
      ISAVE(7) = X
      ISAVE(8) = Z
      RSAVE(1) = BNRM2
      RSAVE(2) = RHO
      RSAVE(3) = RHO1
      RSAVE(4) = CURV
      RETURN
C
C  Non-executable statements
C
 2000 FORMAT( / ' Error message from MI21A/AD. INFO(1) = ', I4 )
 2010 FORMAT( / ' Warning message from MI21A/AD. INFO(1) = ', I4 )
 2020 FORMAT( ' Convergence tolerance out of range.' )
 2030 FORMAT( ' Number of iterations required exceeds the maximum of ',
     *        I8, / ' allowed by ICNTL(6)' )
C
C  End of MI21A/AD
C
      END

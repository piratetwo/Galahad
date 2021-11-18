! COPYRIGHT (c) 1995 Council for the Central Laboratory
!                    of the Research Councils
!
! Version 1.2.0
! See ChangeLog for version history.
!
      SUBROUTINE MI24ID( ICNTL, CNTL, ISAVE, RSAVE, LSAVE )
      INTEGER          ICNTL( 8 )
      DOUBLE PRECISION CNTL( 4 )
      INTEGER ISAVE(17)
      DOUBLE PRECISION RSAVE(9)
      LOGICAL LSAVE(4)
      INTEGER I
C
C  MI24 solves the linear system Ax = b using the
C  Generalized Minimal Residual with restarts (GMRES) iterative method
C  optionally using preconditioning
C                    ^                   ^
C           P(1)AP(2)x = P(1)b,  x = P(2)x

C  P(1), P(2) are the preconditioners, which are not passed to the code,
C  but each time a matrix-vector product with P(1) or P(2) is required,
C  control is passed back to the user.
C
C  Similarly, the matrix A is not passed to the code, but when
C  a matrix-vector with A is required, control is
C  passed back to the user.
C
C  MI24I/ID is the initialisation routine for MI24A/AD and should
C  be called once prior to calls to MI24A/AD.
C
C  Argument list.
C
C  ICNTL   (output) INTEGER control array, dimension 8.
C          ICNTL(1) is the stream number for error messages.
C             On exit, ICNTL(2) = 6.
C          ICNTL(2) is the stream number for warning messages.
C             On exit, ICNTL(2) = 6.
C          ICNTL(3) indicates whether the user wishes to use a
C             preconditioner. If ICNTL(3) is 1 preconditioning will
C             be performed on the left, if it is 2 preconditioning will
C             be performed on the right, and if it is zero, no
C             preconditioning is required.
C             On exit, ICNTL(3) = 0
C          ICNTL(4) indicates whether
C             the convergence test offered by MI24A/AD is to be used.
C             If ICNTL(4) = 0, the computed solution x is accepted
C             if ||Ax - b||/||b|| < CNTL(1).
C             Otherwise, the user most perform his/her
C             own test for convergence when IACT = 1 is returned.
C             On exit, ICNTL(4) = 0
C          ICNTL(5) indicates whether the user wishes to supply an
C             initial guess for the solution vector x.
C             If ICNTL(5) = 0, the user does not wish to supply
C             an initial guess and x = (0,0,...,0) will be used
C             as the initial guess. Otherwise, the user
C             must supply an intial guess on the first call to
C             MI24A/AD. On exit, ICNTL(5) = 0.
C          ICNTL(6) determines the maximum number of iterations
C             allowed. It has default value -1 and, in this case,
C             the maximum number will be N. If the user does
C             not want the maximum number to be N, ICNTL(6) should
C             be set to the maximum  number the user wishes
C             to allow.
C          ICNTRL(7) and ICNTRL(8) are spare and set to zero.
C  CNTL    (output)  REAL (DOUBLE PRECISION) control array, dimension 4.
C          CNTL(1) is a convergence tolerance.
C             On exit, set to square root of machine precision.
C             If ICNTL(4) is nonzero, CNTL(1) not accessed by MI24A/AD.
C          CNTL(2) is a convergence tolerance.
C             On exit, set to zero.
C             If ICNTL(4) is nonzero, CNTL(2) not accessed by MI24A/AD.
C          CNTRL(3) and CNTRL(4) are spare and set to zero.
C
C  ISAVE   (output) INTEGER ARRAY, length 17, used to hold the
C          routine's persistent integer data.  The array contents
C          must not be altered by the user.
C
C  RSAVE   (output) DOUBLE PRECISION ARRAY, length 9, used to hold the
C          routine's persistent real data.  The array contents
C          must not be altered by the user.
C
C  LSAVE   (output) LOGICAL ARRAY, length 4, used to hold the
C          routine's persistent logical data.  The array contents
C          must not be altered by the user.
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

      CNTL(3) = 0.0
      CNTL(4) = 0.0

C  Initialize persistent data to avoid undefined assignment in MI24AD
      DO 10 I = 1, 17
      ISAVE(I) = 0
   10 CONTINUE
      DO 20 I = 1, 9
      RSAVE(I) = 0.0
   20 CONTINUE
      DO 30 I = 1, 4
      LSAVE(I) = .FALSE.
   30 CONTINUE

      RETURN
C
C  End of MI24I/MI24ID
C
      END

      SUBROUTINE MI24AD( IACT, N, M, W, LDW, LOCY, LOCZ, H, LDH, RESID,
     *                   ICNTL, CNTL, INFO, ISAVE, RSAVE, LSAVE )
      DOUBLE PRECISION RESID
      INTEGER          IACT, N, M, LDW, LOCY, LOCZ, LDH
      DOUBLE PRECISION CNTL( 4 ), W( LDW, M + 7 ), H( LDH, M + 2 )
      INTEGER          ICNTL( 8 ), INFO( 4 )
      INTEGER ISAVE(17)
      DOUBLE PRECISION RSAVE(9)
      LOGICAL LSAVE(4)
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
C          the residual in RESID for convergence and recall MI24A/AD
C          if convergence has not been achieved.
C   2      The user must perform the matrix-vector product
C          y := Az,
C          and recall MI24A/AD. The vectors  y and z are held in the
C          first N entries of columns LOCZ and LOCZ of array W,
C          respectively.
C   3      The user must perform the left preconditioning operation
C          y := P_L z,  where P_L is the left preconditioner
C          and recall MI24A/AD. The vectors y and z
C          are  held in the first N entries of
C          columns LOCZ and LOCZ of array W, respectively.
C   4      The user must perform the right preconditioning operation
C          y := P_R z,  where P_R is the right preconditioner
C          and recall MI24A/AD. The vectors y and z
C          are  held in the first N entries of
C          columns LOCZ and LOCZ of array W, respectively.
C
C  N       (input) INTEGER.
C          On entry, the dimension of the matrix.
C          Unchanged on exit.
C
C  W       (right-hand side, solution and workspace, input/output)
C          REAL (DOUBLE PRECISION) array, dimension (LDW,M+7).
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
C          MI24A/AD.
C
C  LDW     (input) INTEGER
C          The leading dimension of the array W. LDW >= max(1,N).
C
C  LOCY, LOCZ (output) INTEGER
C          On exit with IACT > 1, LOCY and LOCZ define the columns
C          of the array W which hold y and z.
C          (see IACT).
C
C  H       (workspace, output)
C          DOUBLE PRECISION array, dimension (LDH,M+2).
C          This workspace is used for constructing and storing the
C          upper Hessenberg matrix. The two extra columns are used to
C          store the Givens rotation matrices.
C
C  LDH    (input) INTEGER
C          The leading dimension of the array H. LDH >= M+1.
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
C          ICNTL may be initalised by calling MI24I/ID.
C          See MI24I/ID for details.
C
C  CNTL    (input) REAL (DOUBLE PRECISION) control array of dimension 4.
C          CNTL may be initalised by calling MI24I/ID.
C          See MI24I/ID for details.
C
C  INFO    (output) INTEGER ARRAY , length 4.
C          If INFO(1) = 0 on exit, no errors or warnings issued.
C          If INFO(1) > 0 on exit, warning to user.
C
C             Warning message:
C
C             1: value of CNTL(1) is out-of-range (u,1.0)
C              The default sqrt(u) is used, u=machine precision.
C
C          If INFO(1) < 0 on exit, an error has occurred.
C
C                Error input data:
C
C                   -1: matrix dimension N < 1
C                   -2: restart interval M < 1
C                   -3: LDW < N
C                   -4: LDH < M + 1
C                   -5: Too many iterations have been performed
C
C          On each exit, INFO(2) holds the number of
C          iterations performed so far.
C
C  ISAVE   (output) INTEGER ARRAY, length 17, used to hold the
C          routine's persistent integer data.  The array contents
C          must not be altered by the user.
C
C  RSAVE   (output) DOUBLE PRECISION ARRAY, length 9, used to hold the
C          routine's persistent real data.  The array contents
C          must not be altered by the user.
C
C  LSAVE   (output) LOGICAL ARRAY, length 4, used to hold the
C          routine's persistent logical data.  The array contents
C          must not be altered by the user.
C
C  BLAS CALLS:   DAXPY, DCOPY, DDOT, DNRM2, DROT, DROTG, DSCAL,
C                DTRSV, DGEMV
C
      DOUBLE PRECISION    ZERO, ONE, POINT1
      PARAMETER         ( ZERO = 0.0D+0, ONE = 1.0D+0, POINT1 = 1.0D-1 )
      INTEGER             I, K, ITMAX, CS, SN, R, S, V, UU, Y, RES
      INTEGER             B, X, IPOS, U
      LOGICAL             LEFT, RIGHT
      DOUBLE PRECISION    AA, BB, BNRM2 , RNORM , DDOT, DNRM2, RSTOP
      DOUBLE PRECISION    PRESID, PRSTOP
      EXTERNAL            DAXPY, DCOPY, DDOT, DNRM2, DROT, DROTG, DSCAL
      EXTERNAL            DTRSV, DGEMV
C
C  Restore persistent data
C
      IPOS   = ISAVE(1)
      ITMAX  = ISAVE(2)
      B      = ISAVE(3)
      I      = ISAVE(4)
      K      = ISAVE(5)
      R      = ISAVE(6)
      X      = ISAVE(7)
      U      = ISAVE(8)
      V      = ISAVE(9)
      S      = ISAVE(10)
      Y      = ISAVE(11)
      CS     = ISAVE(12)
      SN     = ISAVE(13)
      UU     = ISAVE(14)
      RES    = ISAVE(15)

      BNRM2  = RSAVE(1)
      AA     = RSAVE(2)
      BB     = RSAVE(3)
      RNORM  = RSAVE(4)
      PRESID = RSAVE(5)
      RSTOP  = RSAVE(6)
      PRSTOP = RSAVE(7)

      LEFT   = LSAVE(1)
      RIGHT  = LSAVE(2)
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
         GO TO ( 40, 60, 70, 100, 110, 120, 160 ), IPOS
      END IF
C
C  Initial call
C
      INFO( 1 ) = 0
C
C     Test the input parameters.
C
      IF ( N .LT. 1 ) THEN
         INFO( 1 ) = - 1
      ELSE IF ( M .LT. 1 ) THEN
         INFO( 1 ) = - 2
      ELSE IF ( LDW .LT. MAX( 1, N ) ) THEN
         INFO( 1 ) = - 3
      ELSE IF ( LDH .LT. M + 1 ) THEN
         INFO( 1 ) = - 4
      ENDIF
      IF ( INFO( 1 ) .LT. 0 ) THEN
         IACT = - 1
         IF ( ICNTL( 1 ) .GT. 0 ) WRITE( ICNTL( 1 ), 2000 ) INFO( 1 )
         GO TO 1000
      END IF
C
C  Set INFO(2) and ITMAX
C
      INFO( 2 ) = 0
      IF ( ICNTL( 6 ) .GT. 0 ) THEN
         ITMAX = ICNTL( 6 )
      ELSE
         ITMAX = 2 * N
      END IF
C
C  Alias workspace columns
C
      RES = 1
      X   = 2
      S   = 3
      B   = 4
      UU  = 5
      Y   = 6
      V   = 7
C
C  Store the Givens parameters in matrix H
C
      CS = M + 1
      SN = CS + 1
C
C  Compute ||b||
C
      BNRM2 = DNRM2( N, W( 1, RES ), 1 )
C
C  Immediate return if ||b|| = 0
C
      IF ( BNRM2. EQ. ZERO ) THEN
         IACT = 1
         DO 10 I = 1, N
            W( I, X ) = ZERO
            W( I, RES ) = ZERO
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
C  LEFT indicates that a left preconditioner will be used
C  RIGHT indicates that a right preconditioner will be used
C
      LEFT  = ICNTL( 3 ) .EQ. 1 .OR. ICNTL( 3 ) .EQ. 3
      RIGHT = ICNTL( 3 ) .EQ. 2 .OR. ICNTL( 3 ) .EQ. 3
C
C  If no initial guess is provided, set x to zero
C
      IF ( ICNTL( 5 ) .EQ. 0 ) THEN
         DO 20 I = 1, N
            W( I, X ) = ZERO
   20    CONTINUE
      END IF
C
C  Start computing the residual
C
      CALL DCOPY( N, W( 1, RES ), 1, W( 1, B ), 1 )
      IF ( DNRM2( N, W( 1, X ), 1 ) .EQ. ZERO ) GO TO 50
C
C  Start of the main iteration
C
   30 CONTINUE
C
C  Return to user for matrix-vector product Y = A * X
C
         IPOS = 1
         IACT = 2
         LOCY = Y
         LOCZ = X
         GO TO 1000
   40    CONTINUE
C
C  Finalize the residual
C
         CALL DAXPY( N, - ONE, W( 1, Y ), 1, W( 1, RES ), 1 )
   50    CONTINUE
C
C  Compute the norm of the residual
C
         RESID = DNRM2( N, W( 1, RES ), 1 )
C
C  Assign the stopping tolerence on the first iteration
C
         IF ( INFO( 2 ) .EQ. 0 ) THEN
            RSTOP  = MAX( RESID * CNTL( 1 ), CNTL( 2 ) )
            PRSTOP = RSTOP
         END IF
C
C  If too many iterations have occured, exit
C
         IF ( INFO( 1 ) .LT. 0 ) THEN
            IACT = - 1
            IF ( ICNTL( 1 ) .GT. 0 ) WRITE( ICNTL( 1 ), 2000 ) INFO( 1 )
            GO TO 1000
         END IF
C
C  Return to user to obtain preconditioner R = P_L^-1 RES
C
         IF ( LEFT ) THEN
            R = UU
            IPOS = 2
            IACT = 3
            LOCY = R
            LOCZ = RES
            GO TO 1000
         ELSE
            R = RES
         END IF
   60    CONTINUE
C
C  Check for convergence
C
         IF ( ICNTL( 4 ) .NE. 0 .OR. ( ICNTL( 4 ) .EQ. 0 .AND.
     *        RESID .LE. RSTOP ) ) THEN
            IACT = 1
            IPOS = 3
            GO TO 1000
         END IF
   70    CONTINUE
C
C  Construct the first column of V
C
         CALL DCOPY( N, W( 1, R ), 1, W( 1, V ), 1 )
         RNORM = DNRM2( N, W( 1, V ), 1 )
         CALL DSCAL( N, ONE / RNORM, W( 1, V ), 1 )
C
C  Initialize S to the elementary vector E1 scaled by RNORM
C
         W( 1, S ) = RNORM
         DO 80 K = 2, N
            W( K, S ) = ZERO
   80    CONTINUE
C
C  Start of the inner iteration
C
         I = 0
   90    CONTINUE
            I = I + 1
C
C  Update iteration count
C
            INFO( 2 ) = INFO( 2 ) + 1
C
C  Check maximum number of iterations has not been exceeded
C
            IF ( INFO( 2 ) .GT. ITMAX ) THEN
               I = I - 1
               INFO( 1 ) = - 5
               IF ( ICNTL( 1 ) .GT. 0 ) WRITE( ICNTL( 1 ), 2030 ) ITMAX
               IF ( I .NE. 0 ) GO TO 150
               IACT = - 1
               IF ( ICNTL( 1 ) .GT. 0 )
     *            WRITE( ICNTL( 1 ), 2000 ) INFO( 1 ) 
               GO TO 1000
            END IF
C
C  Return to user to obtain preconditioner Y = P_R^-1 V
C
            IF ( RIGHT ) THEN
               IPOS = 4
               IACT = 4
               LOCY = Y
               LOCZ = V + I - 1
               GO TO 1000
            END IF
  100       CONTINUE
C
C  Return to user for matrix-vector product
C
            IPOS = 5
            IACT = 2
            IF ( RIGHT ) THEN
               LOCY = RES
               LOCZ = Y
            ELSE
               LOCY = RES
               LOCZ = V + I - 1
            END IF
            GO TO 1000
  110       CONTINUE
C
C  Return to user to obtain preconditioner W = P_L^-1 RES
C
            IF ( LEFT ) THEN
               U = UU
               IPOS = 6
               IACT = 3
               LOCY = UU
               LOCZ = RES
               GO TO 1000
            ELSE
               U = RES
            END IF
  120       CONTINUE
C
C  Construct I-th column of H orthonormal to the previous I-1 columns
C  using the Gram-Schmidt process on V and U
C
C           CALL MI24BD( I, N, H( 1, I ), W( 1, V ), LDW, W( 1, U ) )
            DO 130 K = 1, I
               H( K, I ) = DDOT( N, W( 1, U ), 1,
     *                              W( 1, V + K - 1 ), 1 )
               CALL DAXPY( N, - H( K, I ), W( 1, V + K - 1 ), 1,
     *                                     W( 1, U ), 1 )
  130       CONTINUE
            H( I + 1, I ) = DNRM2( N, W( 1, U ), 1 )
            CALL DCOPY( N, W( 1, U ), 1, W( 1, V + I ), 1 )
            CALL DSCAL( N, ONE / H( I + 1, I ), W( 1, V + I ), 1 )
C
C  Apply Givens rotations to the I-th column of H. This "updating" of
C  the QR factorization effectively reduces the Hessenberg matrix to
C  upper triangular form during the M iterations.
C
            DO 140 K = 1, I - 1
               CALL DROT( 1, H( K, I ), LDH, H( K + 1, I ), LDH,
     *                    H( K, CS ), H( K, SN ) )
  140       CONTINUE
C
C  Construct the I-th rotation matrix, and apply it to H so that
C  H(I+1,I) = 0
C
            AA = H( I, I )
            BB = H( I + 1, I )
            CALL DROTG( AA, BB, H( I, CS ), H( I, SN ) )
            CALL DROT( 1, H( I, I ), LDH, H( I + 1, I ), LDH,
     *            H( I, CS ), H( I, SN ) )
C
C  Apply the I-th rotation matrix to [ S(I), S(I+1) ]'. This
C  gives an approximation of the residual norm. If less than
C  tolerance, update the approximation vector X and quit
C
            IF ( I .LT. N ) THEN
               CALL DROT( 1, W( I, S ), LDW, W( I + 1, S ),
     *                    LDW, H( I, CS ), H( I, SN ) )
               PRESID = ABS( W( I + 1, S ) )
               IF ( PRESID .LE. PRSTOP .AND. ICNTL( 4 ) .EQ. 0 ) THEN
                  PRSTOP = PRSTOP * POINT1
                  GO TO 150
               END IF
               IF ( I .LT. M )  GO TO 90
            END IF
C
C  Compute current solution vector X
C
  150    CONTINUE
         CALL DCOPY( I, W( 1, S ), 1, W( 1, Y ), 1 )
         CALL DTRSV( 'UPPER', 'NOTRANS', 'NONUNIT', I, H, LDH,
     *               W( 1, Y ), 1 )
C
C  Compute current update vector UU = V*Y
C
         CALL DGEMV( 'NOTRANS', N, I, ONE, W( 1, V ), LDW, W( 1, Y ), 1,
     *               ZERO, W( 1, UU ), 1 )
C
C  Return to user to obtain preconditioner Y = P_R^-1 UU
C
         IF ( RIGHT ) THEN
            IPOS = 7
            IACT = 4
            LOCY = Y
            LOCZ = UU
            GO TO 1000
         END IF
  160    CONTINUE
C
C  Update X
C
         IF ( RIGHT ) THEN
            CALL DAXPY( N, ONE, W( 1, Y  ), 1, W( 1, X ), 1 )
         ELSE
            CALL DAXPY( N, ONE, W( 1, UU ), 1, W( 1, X ), 1 )
         END IF
C
C  Start computing the residual
C
         CALL DCOPY( N, W( 1, B ), 1, W( 1, RES ), 1 )
C
C  Restart
C
      GO TO 30
C
C  Save persistent data and return to caller
C
 1000 CONTINUE
      ISAVE(1) = IPOS
      ISAVE(2) = ITMAX
      ISAVE(3) = B
      ISAVE(4) = I
      ISAVE(5) = K
      ISAVE(6) = R
      ISAVE(7) = X
      ISAVE(8) = U
      ISAVE(9) = V
      ISAVE(10) = S
      ISAVE(11) = Y
      ISAVE(12) = CS
      ISAVE(13) = SN
      ISAVE(14) = UU
      ISAVE(15) = RES

      RSAVE(1) = BNRM2
      RSAVE(2) = AA
      RSAVE(3) = BB
      RSAVE(4) = RNORM
      RSAVE(5) = PRESID
      RSAVE(6) = RSTOP
      RSAVE(7) = PRSTOP

      LSAVE(1) = LEFT
      LSAVE(2) = RIGHT
      RETURN
C
C  Non-executable statements
C
 2000 FORMAT( / ' Error message from MI24A/AD. INFO(1) = ', I4 )
 2010 FORMAT( / ' Warning message from MI24A/AD. INFO(1) = ', I4 )
 2020 FORMAT( ' Convergence tolerance out of range.' )
 2030 FORMAT( /, ' # iterations required exceeds the maximum of ',
     *        I8, ' allowed by ICNTL(6)' )
C
C  End of MI24A/AD
C
      END

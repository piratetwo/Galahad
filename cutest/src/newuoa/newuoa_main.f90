!   ( Last modified on 27 Jan 2013 at 17:20:00 )

      PROGRAM NEWUOA_main

!  NEWUOA test driver for problems derived from SIF files

!  Nick Gould, January 2013

      USE CUTEst_problem

      IMPLICIT NONE
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER :: maxfun, lw, status, iprint, i, ierr, npt
      REAL( KIND = wp ) :: rhobeg, rhoend
      REAL( KIND = wp ), PARAMETER :: infty = 1.0D+19
      REAL( KIND = wp ), DIMENSION( : ), ALLOCATABLE :: W
      REAL( KIND = wp ), DIMENSION( 2 ) :: CPU
      REAL( KIND = wp ), DIMENSION( 4 ) :: CALLS
      INTEGER :: io_buffer = 11
      INTEGER, PARAMETER :: input = 55, indr = 46, out = 6

!  open the relevant file

      OPEN( input, FILE = 'OUTSDIF.d', FORM = 'FORMATTED', STATUS = 'OLD' )
      REWIND( input )

!  initialize problem data structure

! No need for Hessian of objective/Lagrangian

      CUTEST_problem_global%allocate_H = .FALSE.

! No need for Jacobian of constraints

      CUTEST_problem_global%allocate_J = .FALSE.   

      CALL CUTEST_problem_setup( status, CUTEST_problem_global, input )
      IF ( status /= 0 ) GO TO 910

!  set up the data structures necessary to hold the problem functions.

      CALL CUTEST_usetup( status, input, out, io_buffer,                       &
                          CUTEST_problem_global%n, CUTEST_problem_global%x,    &
                          CUTEST_problem_global%x_l, CUTEST_problem_global%x_u )
      IF ( status /= 0 ) GO TO 910
      CLOSE( input )

!  open the Spec file for the method

      OPEN( indr, FILE = 'NEWUOA.SPC', FORM = 'FORMATTED', STATUS = 'OLD')
      REWIND( indr )

!  read input Spec data

!  RHOBEG = the size of the simplex initially
!  RHOEND = the size of the simplex at termination
!  NPT = the number of interpolation conditions; <=0 defaults to 2n+1
!  MAXFUN = the maximum number of function calls allowed.
!  IPRINT should be set to 0, 1, 2 or 3, it controls the amount of printing

!  set up algorithmic input data

      READ ( indr, 1000 ) rhobeg, rhoend, npt, maxfun, iprint
      CLOSE ( indr )

!  ensure that npt satsfies interpolation limits

      IF ( npt <= 0 ) npt = 2 * CUTEST_problem_global%n + 1
      npt = MIN( MAX( npt, CUTEST_problem_global%n + 2 ),                      &
        ( CUTEST_problem_global%n + 1 ) * ( CUTEST_problem_global%n + 2 ) / 2 )

!  allocate the temporary work array W of length at least
!  ( npt + 5 ) * ( npt + n ) + 3 * n * ( n + 3 ) / 2 ... so use double this

      lw = 2 * ( npt + 13 ) * ( npt + CUTEST_problem_global%n )                &
             + 3 * CUTEST_problem_global%n * ( CUTEST_problem_global%n + 3 )
      ALLOCATE( W( lw ), STAT = status )
      IF ( status /= 0 ) GO TO 990

!  perform the minimization

      CALL NEWUOA( CUTEST_problem_global%n, npt, CUTEST_problem_global%X,      &
                   rhobeg, rhoend, iprint, maxfun, W )

!  output report

      CALL CUTEST_ureport( status, CALLS, CPU )
      IF ( status /= 0 ) GO TO 910

      CALL CUTEST_unames( status, CUTEST_problem_global%n,                     &
                          CUTEST_problem_global%pname,                         &
                          CUTEST_problem_global%vnames )
      WRITE( out, 2110 ) ( i, CUTEST_problem_global%vnames( i ),               &
         CUTEST_problem_global%x( i ), CUTEST_problem_global%x_l( i ),         &
         CUTEST_problem_global%x_u( i ), i = 1, CUTEST_problem_global%n )
      WRITE( out, 2000 ) CUTEST_problem_global%pname, CUTEST_problem_global%n, &
         CALLS( 1 ), CUTEST_problem_global%f, CPU( 1 ), CPU( 2 )

!  clean-up data structures

      CALL CUTEST_problem_terminate( status, CUTEST_problem_global )
      IF ( status /= 0 ) GO TO 910
      DEALLOCATE( W, STAT = ierr )
      CALL CUTEST_cterminate( status )
      STOP

!  error returns

  910 CONTINUE
      WRITE( 6, "( ' CUTEst error, status = ', i0, ', stopping' )") status
      STOP

  990 CONTINUE
      WRITE( out, "( ' Allocation error, status = ', I0 )" ) status
      STOP

!  Non-executable statements

2000 FORMAT( /, 24('*'), ' CUTEst statistics ', 24('*') //,                    &
          ' Package used            :  NEWUOA ',  /,                           &
          ' Problem                 :  ', A10,    /,                           &
          ' # variables             =      ', I10 /,                           &
          ' # objective functions   =        ', F8.2 /,                        &
          ' Final f                 = ', E15.7 /,                              &
          ' Set up time             =      ', 0P, F10.2, ' seconds' /,         &
          ' Solve time              =      ', 0P, F10.2, ' seconds' //,        &
          66('*') / )
1000 FORMAT( D12.4, /, D12.4, /, I6, /, I6, /, I6 )
2110 FORMAT( /, ' The variables:', /, &
          '     i name          value    lower bound upper bound',             &
          /, ( I6, 1X, A10, 1P, 3D12.4 ) )

!  End of NEWUOA_main

      END PROGRAM NEWUOA_main

      SUBROUTINE CALFUN( n, X, f )

!  evaluates the objective function value in a format compatible with NEWUOA,
!  but using the CUTEst tools.

      USE CUTEst_problem

      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER, INTENT( IN ) :: n
      REAL( KIND = wp ), INTENT( OUT ) :: f
      REAL( KIND = wp ), INTENT( IN ) :: X( n )

      INTEGER ::  status

!  Evaluate the objective function and constraints.

      CALL CUTEST_ufn( status, CUTEST_problem_global%n,                        &
                       X, CUTEST_problem_global%f )
      IF ( status /= 0 ) GO TO 910
      f = CUTEST_problem_global%f
      RETURN

  910 CONTINUE
      WRITE( 6, "( ' CUTEst error, status = ', i0, ', stopping' )" ) status
      STOP

!  End of CALFUN

      END SUBROUTINE CALFUN


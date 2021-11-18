! COPYRIGHT (c) 2014 Science and Technology Facilities Council
! Original date 1 July 14, Version 1.0.0

!-*-*-*-*-*-*-*-*-*  H S L _ M I 3 2  double  M O D U L E  *-*-*-*-*-*-*-*-*-*

MODULE HSL_MI32_DOUBLE
      IMPLICIT NONE

!  --------------------------------------------------------------------
!
!  Solve the symmetric linear system
!
!     subject to       A x = b
!
!  using the MINRES method. A need not be definite
!
!  Tyrone Rees
!  June 2014
!
!  --------------------------------------------------------------------

!  Double precision version

      PRIVATE
      PUBLIC :: MI32_MINRES, MI32_FINALIZE

!  Define the working precision to be double

      INTEGER, PARAMETER :: working = KIND( 1.0D+0 )

!  Paremeter definitions

      REAL ( KIND = working ), PARAMETER :: zero = 0.0_working
      REAL ( KIND = working ), PARAMETER :: one = 1.0_working

!  Derived type definitions

      TYPE, PUBLIC :: MI32_CONTROL
         REAL ( KIND = working ) :: stop_relative = sqrt(epsilon(1.0_working))
         REAL ( KIND = working ) :: stop_absolute = 0.0_working
         INTEGER :: out = -1
         INTEGER :: error = 6
         INTEGER :: itmax = -1
         INTEGER :: conv_test_norm = 1
         LOGICAL :: precondition = .true.
         LOGICAL :: own_stopping_rule = .false.
      END TYPE

      TYPE, PUBLIC :: MI32_INFO
         INTEGER :: st
         REAL ( KIND = working ) :: rnorm
         INTEGER :: iter
      END TYPE

      TYPE, PUBLIC :: MI32_KEEP
         PRIVATE
         INTEGER :: itmax
         INTEGER :: vold, vnew, wnew, wold, rold, rnew, zold, znew, gold, gnew
         REAL ( KIND = working ) :: rstop, stop_relative
         REAL ( KIND = working ) :: eta, diag, subd, tiny
         REAL ( KIND = working ), ALLOCATABLE, DIMENSION( : , : ) :: V
         REAL ( KIND = working ), ALLOCATABLE, DIMENSION( : , : ) :: W
         REAL ( KIND = working ), ALLOCATABLE, DIMENSION( : , : ) :: Z
         ! and ones I've added...
         REAL ( KIND = working ), DIMENSION( 2 ) :: co, si, gamma
         REAL ( KIND = working ), ALLOCATABLE :: R(:) !hold residual (if needed)
         REAL ( KIND = working ), ALLOCATABLE :: AW(:,:) ! Aw_0 and Aw_1
      END TYPE

   CONTAINS

!-*-*-*-*-*-  H S L _ M I 3 2 _ M I N R E S   S U B R O U T I N E   -*-*-*-*-*-

      SUBROUTINE MI32_MINRES( action, n, X, V_in, V_out, keep, control, info )

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
!  . Main MINRES iteration for A x = b .
!
!  Arguments:
!  =========
!
!   action  RCI action to complete. This must be set to 1 on initial
!                    entry and the residual vector A x - b placed in V_in.
!                    On exit, the value of status requires the following
!                    action from the user:
!           status = 0, the solution has been found
!           status = 2, the preconditioner must be applied to V_out,
!              with the result returned in V_in, and the subroutine
!              re-entered with all other arguments unchanged
!           status = 3, the product A * V_out must be formed, with
!              the result returned in V_in and the subroutine
!              re-entered with all other arguments unchanged
!           status = -1, n is not positive
!           status = -2, the iteration limit has been exceeded
!           status = -3, the system appears to be inconsistent
!           status = -4, an allocation error has occured. The
!              allocation status returned by ALLOCATE has the value
!              info%st
!   n       number of unknowns
!   X       the vector of unknowns. On entry, an initial estimate of the
!           solution, on exit, the best value found so far
!   V_in    see action = 1, 2, and 3
!   V_out   see action = 2, and 3
!   data    private internal data
!   control
!   info    a structure containing information. Components are
!            iter    the number of iterations performed
!            st      status of the most recent (de)allocation
!            rnorm   the M-norm of the residual
!
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------

      INTEGER, INTENT( INOUT ) :: action
      INTEGER, INTENT( IN ) :: n
      REAL ( KIND = working ), INTENT( INOUT ), DIMENSION( n ) :: X, V_in
      REAL ( KIND = working ), POINTER, DIMENSION( : ) :: V_out
      TYPE ( MI32_KEEP ), TARGET, INTENT( INOUT ) :: keep
      TYPE ( MI32_CONTROL ), INTENT( IN ) :: control
      TYPE ( MI32_INFO ), INTENT( INOUT ) :: info

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

      INTEGER :: i
      REAL ( KIND = working ) :: delta = zero
      REAL ( KIND = working ) :: alp(4), ialp2
      LOGICAL :: reallocate
      CHARACTER ( LEN = 6 ) :: bad_alloc

!-----------------------------------------------
!   E x t e r n a l   P r o c e d u r e s
!-----------------------------------------------

!      EXTERNAL DLAEV2

!  Branch to appropriate part of the code

      SELECT CASE ( action )
      CASE ( 1 )  ! initialization
         GO TO 100
      CASE ( 2 )  ! after applying preconditioner
         GO TO 200
      CASE ( 3 )  ! after applying matrix-vector product
         GO TO 300
      CASE ( 4 )  ! after users' stopping rule
         GO TO 250
      END SELECT

!  On initial entry, set constants

  100 CONTINUE

!  Set initial values

      info%iter = 0 ; info%st = 0

!  Ensure that n is positive

      IF ( n <= 0 ) THEN
         IF ( control%out >= 0 ) WRITE( control%out, 2040 )
         action = - 1
         GO TO 920
      END IF

!  Ensure that input parameters are in the right range

      keep%itmax = control%itmax
      IF ( keep%itmax < 0 ) keep%itmax = 2 * n! + 1
      keep%stop_relative =                                                    &
           MIN( ONE, MAX( control%stop_relative, EPSILON( one ) ) )

      if ( ( control%conv_test_norm < 1 ).or.( control%conv_test_norm > 2)) then
         IF ( control%out >= 0 ) WRITE( control%out, 2060 )
         action = -5
         GO TO 920
      end if

!  Allocate Q, W and S
      reallocate = .TRUE.
      IF ( ALLOCATED( keep%Z ) ) THEN
         IF ( SIZE( keep%Z, 1 ) < n .OR.                                      &
              SIZE( keep%Z, 2 ) < 2 ) THEN ; DEALLOCATE( keep%Z )
         ELSE ; reallocate = .FALSE.
         END IF
      END IF
      IF ( reallocate ) THEN
         ALLOCATE( keep%Z( n, 2 ), STAT = info%st )
         IF ( info%st /= 0 ) THEN ; bad_alloc = 'Z' ; GO TO 910
         END IF
      END IF

      reallocate = .TRUE.
      IF ( ALLOCATED( keep%W ) ) THEN
         IF ( SIZE( keep%W, 1 ) < n .OR.                                      &
              SIZE( keep%W, 2 ) < 2 ) THEN ; DEALLOCATE( keep%W )
         ELSE ; reallocate = .FALSE.
         END IF
      END IF
      IF ( reallocate ) THEN
         ALLOCATE( keep%W( n, 2 ), STAT = info%st )
         IF ( info%st /= 0 ) THEN ; bad_alloc = 'W' ; GO TO 910
         END IF
      END IF

      reallocate = .TRUE.
      IF ( ALLOCATED( keep%V ) ) THEN
         IF ( SIZE( keep%V, 1 ) < n .OR.                                      &
              SIZE( keep%V, 2 ) < 2 ) THEN ; DEALLOCATE( keep%V )
         ELSE ; reallocate = .FALSE.
         END IF
      END IF
      IF ( reallocate ) THEN
         ALLOCATE( keep%V( n, 2 ), STAT = info%st )
         IF ( info%st /= 0 ) THEN ; bad_alloc = 'V' ; GO TO 910
         END IF
      END IF
      keep%V = zero
      keep%W = zero

      if ( control%conv_test_norm == 1 ) then
         ! do nothing....
         continue
      elseif ( control%conv_test_norm == 2 ) then
         ALLOCATE( keep%AW( n, 3 ), STAT = info%st )
         IF ( info%st /= 0 ) THEN
            bad_alloc = 'AW' ; GO TO 910
         END IF
         keep%AW = zero
         ALLOCATE( keep%R( n ), STAT = info%st )
         IF ( info%st /= 0 ) THEN
            bad_alloc = 'R' ; GO TO 910
         END IF
         keep%R = zero
      end if



!      info%rnorm = SQRT( DOT_PRODUCT( V_in, V_in ) )
      if ( control%conv_test_norm == 2 ) then
         info%rnorm = SQRT( DOT_PRODUCT( V_in, V_in ) )
         keep%R = V_in
!     else
!        info%rnorm = SQRT( DOT_PRODUCT( V_in, V_in ) )
      end if

!  Set other initial values

      keep%wnew = 1 ; keep%wold = 2
      keep%wnew = 1 ; keep%wold = 2
      keep%vnew = 1 ; keep%vold = 2
      keep%rnew = 1 ; keep%rold = 2
      keep%gnew = 1 ; keep%gold = 2

      keep%gamma = one

!      keep%Z( 1 ) = one ; keep%Z( 2 ) = zero
      keep%tiny = EPSILON( one )

      keep%V( : , keep%vnew ) = -V_in

!  ===========================
!  Start of the main iteration
!  ===========================

  190 CONTINUE

!  Return to obtain the product of M(inverse) with t

      IF ( control%precondition ) THEN
         V_out => keep%V( : , keep%vnew )
         action = 2
         RETURN
      ELSE
         V_in = keep%V( : , keep%vnew )
      END IF

!  Compute the off-diagonal of the Lanczos tridiagonal

  200 CONTINUE

      keep%gamma(keep%gold) = &
           SQRT( DOT_PRODUCT( V_in , keep%V( : , keep%vnew ) ) ) ! sqrt(z1'v1)

      i = keep%gnew ; keep%gnew = keep%gold ; keep%gold = i

      keep%Z( : , keep%vnew ) = V_in

      inital_iterate: IF ( info%iter == 0 ) THEN

         keep%eta = keep%gamma(keep%gnew)

         if ( control%conv_test_norm == 1 ) then
            info%rnorm = keep%eta ! use the initial M^{-1}-norm
         end if

         keep%si = zero
         keep%co = one

         !  Compute the stopping tolerence

         keep%rstop = MAX( info%rnorm * keep%stop_relative,   &
                           control%stop_absolute,             &
                           keep%tiny) ! include this so that immediate
                                      ! convergence can be detected



         IF ( control%out >= 0 ) WRITE( control%out, 2000 ) keep%rstop


      ELSE
         alp(1) = keep%co( keep%rnew ) * delta - &
                   keep%co(keep%rold) * keep%si( keep%rnew) * &
                   keep%gamma(keep%gold)
         alp(2) = sqrt( alp(1)**2 + keep%gamma(keep%gnew)**2 )
         ialp2 = one / alp(2)
         alp(3) = keep%si( keep%rnew ) * delta + &
                   keep%co( keep%rold ) * keep%co( keep%rnew ) * &
                   keep%gamma(keep%gold)
         alp(4) = keep%si(keep%rold) * keep%gamma(keep%gold)
         i = keep%rnew ; keep%rnew = keep%rold ; keep%rold = i

         keep%co( keep%rnew ) = alp(1)*ialp2
         keep%si( keep%rnew ) = keep%gamma(keep%gnew)*ialp2

         keep%W( : , keep%wold ) = keep%Z( : , keep%vold) - &
                               alp(4) * keep%W( : , keep%wold) - &
                               alp(3) * keep%W( : , keep%wnew)
         if (control%conv_test_norm == 2) then
            ! check....
            keep%AW( : , keep%wold ) = keep%AW( : , 3 ) - &
                 alp(4) * keep%AW( : , keep%wold ) - &
                 alp(3) * keep%AW( : , keep%wnew )
         end if
         i = keep%wnew ; keep%wnew = keep%wold ; keep%wold = i

         keep%W( : , keep%wnew ) = keep%W( : , keep%wnew ) * ialp2
         X = X + ( keep%co( keep%rnew ) * keep%eta ) * keep%W( : , keep%wnew )
         if (control%conv_test_norm == 2) then
            keep%AW( : , keep%wnew ) = keep%AW( : , keep%wnew ) * ialp2
            keep%R = keep%R + &
                     keep%co( keep%rnew ) * keep%eta * keep%AW(:,keep%wnew)
         end if


         keep%eta = - keep%si( keep%rnew ) * keep%eta
         if (control%conv_test_norm == 1) then
            info%rnorm = abs( keep%eta )
         else
            ! form \|r\|_2
            info%rnorm = sqrt(dot_product(keep%R,keep%R))
         end if

      END IF inital_iterate


!  Allow the user to check for termination

      stop: IF ( control%own_stopping_rule ) THEN
         action = 4
         RETURN

!  Check for termination

      ELSE
         IF ( info%rnorm <= keep%rstop ) THEN
            action = 0
            IF ( info%iter == 0 .AND. control%out >= 0 )                      &
                 WRITE( control%out, 2020 )                                   &
                    info%iter, info%rnorm, keep%gamma
            GO TO 900
         END IF
      END IF stop

! Since we haven't converged, test for singularity

      IF ( keep%gamma(keep%gnew) < keep%tiny ) THEN

         IF ( control%out >= 0 ) WRITE( control%out, 2050 )
         action = - 3
         GO TO 920

      END If


  250 CONTINUE

      !  Start a new iteration

      info%iter = info%iter + 1


      IF ( info%iter >= keep%itmax ) THEN
         IF ( control%out >= 0 ) WRITE( control%out, 2030 )
         action = - 2
         GO TO 920
      END IF


      keep%Z( : , keep%vnew ) = keep%Z( : , keep%vnew ) / keep%gamma(keep%gnew)

      !  Return to obtain the product of H with q

      V_out => keep%Z( : , keep%vnew )
      action = 3
      RETURN

  300 CONTINUE

      delta = dot_product ( V_in , keep%Z( : , keep%vnew ) )
      keep%V( : , keep%vold ) = V_in - &
                            ( delta/keep%gamma(keep%gnew) ) * &
                            keep%V ( : , keep%vnew ) - &
                            ( keep%gamma(keep%gnew)/keep%gamma(keep%gold) ) * &
                              keep%V ( : , keep%vold )

      i = keep%vnew ; keep%vnew = keep%vold ; keep%vold = i

      if ( control%conv_test_norm == 2 ) then
         keep%AW(:,3) = V_in ! store AZ in AW_3
      end if
      GO TO 190

!  Terminal returns

  900 CONTINUE
      RETURN

!  Unsuccessful returns

  910 CONTINUE
      action = - 4
      IF ( control%error > 0 ) WRITE( control%error, 2900 ) bad_alloc

  920 CONTINUE
      RETURN

!  Non-executable statements

 2000 FORMAT( ' stopping tolerence =' , ES12.4, //,                           &
              ' iter     rnorm       diag       subdiag   pivot ' )
! 2010 FORMAT( I6, 3ES12.4, L4 )
 2020 FORMAT( I6, ES12.4, '      -     ', ES12.4, L4 )
 2030 FORMAT( /, ' Iteration limit exceeded ' )
 2040 FORMAT( /, ' n not positive ' )
 2050 FORMAT( /, ' Singularity encountered ' )
 2060 FORMAT( /, 'Unsupported control%conv_test_norm supplied; must be 1 or 2' )
 2900 FORMAT( /, ' ** Message from -MI32_MINRES-', /,                         &
              '    Allocation error for array ', A6 )

   END SUBROUTINE MI32_MINRES

!*-*-*-*-*-  H S L _ M I 3 2 _ F I N A L I Z E   S U B R O U T I N E   -*-*-*-*-

   SUBROUTINE MI32_FINALIZE( keep )
      TYPE ( MI32_KEEP ), INTENT( INOUT ) :: keep

      INTEGER :: st

      ! Deallocate the arrays Z, W and V
      DEALLOCATE( keep%Z, STAT = st )
      DEALLOCATE( keep%W, STAT = st )
      DEALLOCATE( keep%V, STAT = st )

   END SUBROUTINE MI32_FINALIZE

END MODULE HSL_MI32_DOUBLE

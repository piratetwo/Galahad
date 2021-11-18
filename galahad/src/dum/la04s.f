! THIS VERSION: 19/10/2018 AT 08:35:00 GMT.

!-*-*-*-*-*-*-  G A L A H A D  -  L A 0 4  S U B R O U T I N E S *-*-*-*-

!  Nick Gould, for GALAHAD productions
!  Copyright reserved
!  October 19th 2018

      SUBROUTINE LA04A( A, la, IRN, IP, m, n, B, C, BND, kb, lb, job, CNTL,
     *                  IX, JX, X, Z, G, RINFO, WS, lws, IWS, liws )

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------

      INTEGER, INTENT( IN ) :: la, m, n, kb, lb, lws, liws
      INTEGER, INTENT( INOUT ) :: job
      INTEGER, INTENT( INOUT ), DIMENSION( n + 1 ) :: IP
      INTEGER, INTENT( INOUT ), DIMENSION( la ) :: IRN
      INTEGER, INTENT( INOUT ), DIMENSION( m ) :: IX
      INTEGER, INTENT( INOUT ), DIMENSION( kb ) :: JX
      INTEGER, INTENT( INOUT ), DIMENSION( liws ) :: IWS
      REAL ( KIND = KIND( 1.0E+0 ) ), INTENT( INOUT ), DIMENSION( la ) :: A
      REAL ( KIND = KIND( 1.0E+0 ) ), INTENT( INOUT ), DIMENSION( m ) :: B
      REAL ( KIND = KIND( 1.0E+0 ) ), INTENT( INOUT ), DIMENSION( 2, kb ) :: BND
      REAL ( KIND = KIND( 1.0E+0 ) ), INTENT( INOUT ), DIMENSION( n ) :: C
      REAL ( KIND = KIND( 1.0E+0 ) ), INTENT( INOUT ), DIMENSION( 15 ) :: CNTL
      REAL ( KIND = KIND( 1.0E+0 ) ), INTENT( INOUT ), DIMENSION( 40 ) :: RINFO
      REAL ( KIND = KIND( 1.0E+0 ) ), INTENT( INOUT ), DIMENSION( n + m ) :: X
      REAL ( KIND = KIND( 1.0E+0 ) ), INTENT( INOUT ), DIMENSION( n ) :: Z
      REAL ( KIND = KIND( 1.0E+0 ) ), INTENT( INOUT ), DIMENSION( n ) :: G
      REAL ( KIND = KIND( 1.0E+0 ) ), INTENT( INOUT ), DIMENSION( lws ) :: WS

!  Dummy subroutine available with GALAHAD

      WRITE ( 6, 2000 )
      job = - 101
      STOP

!  Non-executable statements

 2000 FORMAT( /,
     *     ' We regret that the solution options that you have ', /,
     *     ' chosen are not all freely available with GALAHAD.', //,
     *     ' If you have HSL (formerly the Harwell Subroutine',
     *     ' Library), this ', /,
     *     ' option may be enabled by replacing the dummy ', /,
     *     ' subroutine LA04AD with its HSL namesake ', /,
     *     ' and dependencies. See ', /,
     *     '   $GALAHAD/src/makedefs/packages for details.', //,
     *     ' *** EXECUTION TERMINATING *** ', / )

!  End of dummy subroutine LA04A

      END SUBROUTINE LA04A

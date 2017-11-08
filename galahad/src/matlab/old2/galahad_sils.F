#include <fintrf.h>

! *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!                 MEX INTERFACE TO GALAHAD_SILS
!
!  Expected Matlab syntax:
!
!    [X,INFORM] = SILS(A,B) or
!    [X,INFORM] = SILS(A,B,CONTROL)
!
! *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   originally released with GALAHAD Version 2.1. July 4th 2007

!  For full documentation, see 
!   http://galahad.rl.ac.uk/galahad-www/specs.html

      SUBROUTINE mexFunction( nlhs, plhs, nrhs, prhs )
      USE GALAHAD_SMT_double
      USE GALAHAD_SILS_DOUBLE
      IMPLICIT NONE
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )

! ------------------------- Do not change -------------------------------

!  Keep the above subroutine, argument, and function declarations for use
!  in all your fortran mex files.
!
      INTEGER * 4 :: nlhs, nrhs
      mwPointer :: plhs( * ), prhs( * )

      CHARACTER ( len = 16 ) :: fi, mxGetFieldNameByNumber
      LOGICAL :: mxIsSparse, mxIsStruct, mxIsChar
      mwSize :: mxGetString
      mwSize :: mxGetM, mxGetN, mxGetIr, mxGetJc
      mwSize :: mxIsNumeric, mxGetNzmax, mxGetNumberOfFields
      mwPointer :: mxCreateStructMatrix, mxGetPr, mxGetField
      mwPointer :: galmxCreateInteger, galmxCreateReal
      mwPointer :: galmxCreateIntegerMatrix, galmxCreateRealMatrix

! -----------------------------------------------------------------------

!  local variables

      mwSize :: i, j, k, l, dig, stat, n, nb, nfields, info, imode
      mwSize :: a_in, b_in, c_in, a_arg, b_arg, c_arg
      mwPointer :: val_in, rhs_in, cpr_in, row_in, cn_in, x_pr
      mwPointer :: flag_pr, more_pr, stat_pr 
      mwPointer :: nsteps_pr, nrltot_pr, nirtot_pr, nrlnec_pr
      mwPointer :: nirnec_pr, nrladu_pr, niradu_pr, ncmpa_pr
      mwPointer :: nfaulty_pr, oor_pr, dup_pr, maxfrt_pr
      mwPointer :: nebdu_pr, nrlbdu_pr, nirbdu_pr, ncmpbr_pr
      mwPointer :: ncmpbi_pr, ntwo_pr, neig_pr, delay_pr
      mwPointer :: signc_pr, modstep_pr, rank_pr, step_pr
      mwPointer :: opsa_pr, opse_pr, opsb_pr, maxchange_pr
      mwPointer :: string_pr

      CHARACTER ( len = 7 ) :: mode
      REAL ( KIND = wp ) :: tot, at, att( 1 ), mxGetScalar 
      REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: B, X
      REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : , : ) :: B2, X2

      mwSize, PARAMETER :: ninform2 = 1
      CHARACTER ( LEN = 9 ), PARAMETER :: 
     *     finform2( ninform2 ) = (/ 'string   ' /)
      mwSize, PARAMETER :: ninform = 31
      CHARACTER ( LEN = 9 ), PARAMETER :: finform( ninform ) = (/ 
     *     'flag     ', 'more     ', 'stat     ',
     *     'nsteps   ', 'nrltot   ', 'nirtot   ', 'nrlnec   ',
     *     'nirnec   ', 'nrladu   ', 'niradu   ', 'ncmpa    ',
     *     'nfaulty  ', 'oor      ', 'dup      ', 'maxfrt   ',
     *     'nebdu    ', 'nrlbdu   ', 'nirbdu   ', 'ncmpbr   ',
     *     'ncmpbi   ', 'ntwo     ', 'neig     ', 'delay    ',
     *     'signc    ', 'modstep  ', 'rank     ', 'step     ',
     *     'opsa     ', 'opse     ', 'opsb     ', 'maxchange'
     *   /)

!  arguments for SILS

      TYPE ( SMT_type ), SAVE :: A
      TYPE ( SILS_control ), SAVE :: CONTROL
      TYPE ( SILS_ainfo ) :: AINFO
      TYPE ( SILS_finfo ) :: FINFO
      TYPE ( SILS_sinfo ) :: SINFO
      TYPE ( SILS_factors ), SAVE :: FACTORS

      mwPointer, ALLOCATABLE :: col_ptr( : )

!  Test input/output arguments

      IF ( nrhs < 1 ) THEN
        CALL mexErrMsgTxt( ' sils requires at least 2 input arguments' )
      END IF

!  possible strings
!   'initial'
!   'factor'
!   'solve'
!   'final'
!   'all'

      IF ( nrhs < 1 ) THEN
        CALL mexErrMsgTxt( ' Too few input arguments to sils' )
      END IF

      IF ( mxIsChar( prhs( 1 ) ) ) THEN
        imode = mxGetString( prhs( 1 ), mode, 7 )
      ELSE
        CALL mexErrMsgTxt( ' First input to sils must be a string' )
      END IF

      IF ( .NOT. ( TRIM( mode ) == 'initial' .OR. 
     &             TRIM( mode ) == 'final' ) ) THEN
        IF ( nrhs < 2 ) THEN
          CALL mexErrMsgTxt( ' Too few input arguments to sils' )
        END IF
        IF ( TRIM( mode ) == 'factor' ) THEN
          a_arg = 2
          c_arg = 3
        ELSE IF ( TRIM( mode ) == 'solve' ) THEN
          b_arg = 2
          c_arg = 3
!       ELSE IF ( TRIM( mode ) == 'all' ) THEN
        ELSE
          a_arg = 2
          b_arg = 3
          c_arg = 4
        END IF
        IF ( nrhs > c_arg ) THEN
          CALL mexErrMsgTxt( ' Too many input arguments to sils' )
        END IF
        IF ( nlhs > 2 ) THEN
          CALL mexErrMsgTxt( 
     &      ' sils provides at most 2 output arguments' )
        END IF
      END IF

!  Initialize the internal structures for sils

      IF ( TRIM( mode ) == 'initial' .OR. 
     &     TRIM( mode ) == 'all' ) THEN
        CALL SILS_INITIALIZE( FACTORS, CONTROL )
        IF ( TRIM( mode ) == 'initial' ) RETURN
      END IF

      IF ( .NOT. TRIM( mode ) == 'final' ) THEN

!  If the third argument is present, extract the input control data

        IF ( nrhs == c_arg ) THEN
          c_in = prhs( c_arg )
          IF ( .NOT. mxIsStruct( c_in ) ) THEN
            CALL mexErrMsgTxt( 
     &        ' 3rd input argument must be a structure' )
          END IF
          nfields = mxGetNumberOfFields( c_in )
          DO i = 1, nfields
            fi = mxGetFieldNameByNumber( c_in, i )
            SELECT CASE ( TRIM( fi ) )
            CASE( 'lp' )
               cn_in = mxGetField( c_in, 1, 'lp' )
               CONTROL%lp = mxGetScalar( cn_in )
            CASE( 'wp' )
              cn_in = mxGetField( c_in, 1, 'wp' )
              CONTROL%wp = mxGetScalar( cn_in )
            CASE( 'mp' )
              cn_in = mxGetField( c_in, 1, 'mp' )
              CONTROL%mp = mxGetScalar( cn_in )
            CASE( 'sp' )
              cn_in = mxGetField( c_in, 1, 'sp' )
              CONTROL%sp = mxGetScalar( cn_in )
            CASE( 'ldiag' )
              cn_in = mxGetField( c_in, 1, 'ldiag' )
              CONTROL%ldiag = mxGetScalar( cn_in )
            CASE( 'factorblocking' )
              cn_in = mxGetField( c_in, 1, 'factorblocking' )
              CONTROL%factorblocking = mxGetScalar( cn_in )
            CASE( 'solveblocking' )
              cn_in = mxGetField( c_in, 1, 'solveblocking' )
              CONTROL%solveblocking = mxGetScalar( cn_in )
            CASE( 'la' )
              cn_in = mxGetField( c_in, 1, 'la' )
              CONTROL%la = mxGetScalar( cn_in )
            CASE( 'liw' )
              cn_in = mxGetField( c_in, 1, 'liw' )
              CONTROL%liw = mxGetScalar( cn_in )
            CASE( 'maxla' )
              cn_in = mxGetField( c_in, 1, 'maxla' )
              CONTROL%maxla = mxGetScalar( cn_in )
            CASE( 'maxliw' )
              cn_in = mxGetField( c_in, 1, 'maxliw' )
              CONTROL%maxliw = mxGetScalar( cn_in )
            CASE( 'pivoting' )
              cn_in = mxGetField( c_in, 1, 'pivoting' )
              CONTROL%pivoting = mxGetScalar( cn_in )
            CASE( 'thresh' )
              cn_in = mxGetField( c_in, 1, 'thresh' )
              CONTROL%thresh = mxGetScalar( cn_in )
            CASE( 'ordering' )
              cn_in = mxGetField( c_in, 1, 'ordering' )
              CONTROL%ordering = mxGetScalar( cn_in )
            CASE( 'scaling' )
              cn_in = mxGetField( c_in, 1, 'scaling' )
              CONTROL%scaling = mxGetScalar( cn_in )
            CASE( 'multiplier' )
              cn_in = mxGetField( c_in, 1, 'multiplier' )
              CONTROL%multiplier = mxGetScalar( cn_in )
            CASE( 'reduce' )
              cn_in = mxGetField( c_in, 1, 'reduce' )
              CONTROL%reduce = mxGetScalar( cn_in )
            CASE( 'u' )
              cn_in = mxGetField( c_in, 1, 'u' )
              CONTROL%u = mxGetScalar( cn_in )
            CASE( 'static_tolerance' )
              cn_in = mxGetField( c_in, 1, 'static_tolerance' )
              CONTROL%static_tolerance = mxGetScalar( cn_in )
            CASE( 'static_level' )
              cn_in = mxGetField( c_in, 1, 'static_level' )
              CONTROL%static_level = mxGetScalar( cn_in )
            CASE( 'tolerance' )
              cn_in = mxGetField( c_in, 1, 'tolerance' )
              CONTROL%tolerance = mxGetScalar( cn_in )
            END SELECT
          END DO
        END IF

!  No i/o allowed from mex interface

        CONTROL%lp = 0
        CONTROL%wp = 0
        CONTROL%mp = 0
        CONTROL%sp = 0

!  Create inform output structure

        plhs( 2 ) = mxCreateStructMatrix( 1, 1, ninform, finform )

!  Define the components of the structure

        flag_pr = galmxCreateInteger( ) 
        CALL mxSetField( plhs( 2 ), 1, 'flag', flag_pr )
        more_pr = galmxCreateInteger( ) 
        CALL mxSetField( plhs( 2 ), 1, 'more', more_pr )
        stat_pr = galmxCreateInteger( ) 
        CALL mxSetField( plhs( 2 ), 1, 'stat', stat_pr )
        nsteps_pr = galmxCreateInteger( )
        CALL mxSetField( plhs( 2 ), 1, 'nsteps', nsteps_pr )
        nrltot_pr = galmxCreateInteger( )
        CALL mxSetField( plhs( 2 ), 1, 'nrltot', nrltot_pr )
        nirtot_pr = galmxCreateInteger( )
        CALL mxSetField( plhs( 2 ), 1, 'nirtot', nirtot_pr )
        nrlnec_pr = galmxCreateInteger( )
        CALL mxSetField( plhs( 2 ), 1, 'nrlnec', nrlnec_pr )
        nirnec_pr = galmxCreateInteger( )
        CALL mxSetField( plhs( 2 ), 1, 'nirnec', nirnec_pr )
        nrladu_pr = galmxCreateInteger( )
        CALL mxSetField( plhs( 2 ), 1, 'nrladu', nrladu_pr )
        niradu_pr = galmxCreateInteger( )
        CALL mxSetField( plhs( 2 ), 1, 'niradu', niradu_pr )
        ncmpa  _pr = galmxCreateInteger( )
        CALL mxSetField( plhs( 2 ), 1, 'ncmpa  ', ncmpa  _pr )
        nfaulty_pr = galmxCreateInteger( )
        CALL mxSetField( plhs( 2 ), 1, 'nfaulty', nfaulty_pr )
        oor_pr = galmxCreateInteger( )
        CALL mxSetField( plhs( 2 ), 1, 'oor', oor_pr )
        dup_pr = galmxCreateInteger( )
        CALL mxSetField( plhs( 2 ), 1, 'dup', dup_pr )
        maxfrt_pr = galmxCreateInteger( )
        CALL mxSetField( plhs( 2 ), 1, 'maxfrt', maxfrt_pr )
        nebdu_pr = galmxCreateInteger( )
        CALL mxSetField( plhs( 2 ), 1, 'nebdu', nebdu_pr )
        nrlbdu_pr = galmxCreateInteger( )
        CALL mxSetField( plhs( 2 ), 1, 'nrlbdu', nrlbdu_pr )
        nirbdu_pr = galmxCreateInteger( )
        CALL mxSetField( plhs( 2 ), 1, 'nirbdu', nirbdu_pr )
        ncmpbr_pr = galmxCreateInteger( )
        CALL mxSetField( plhs( 2 ), 1, 'ncmpbr', ncmpbr_pr )
        ncmpbi_pr = galmxCreateInteger( )
        CALL mxSetField( plhs( 2 ), 1, 'ncmpbi', ncmpbi_pr )
        ntwo_pr = galmxCreateInteger( )
        CALL mxSetField( plhs( 2 ), 1, 'ntwo', ntwo_pr )
        neig_pr = galmxCreateInteger( )
        CALL mxSetField( plhs( 2 ), 1, 'neig', neig_pr )
        delay_pr = galmxCreateInteger( )
        CALL mxSetField( plhs( 2 ), 1, 'delay', delay_pr )
        signc_pr = galmxCreateInteger( )
        CALL mxSetField( plhs( 2 ), 1, 'signc', signc_pr )
        modstep_pr = galmxCreateInteger( )
        CALL mxSetField( plhs( 2 ), 1, 'modstep', modstep_pr )
        rank_pr = galmxCreateInteger( )
        CALL mxSetField( plhs( 2 ), 1, 'rank', rank_pr )
        step_pr = galmxCreateInteger( )
        CALL mxSetField( plhs( 2 ), 1, 'step', step_pr )
        opsa_pr = galmxCreateReal( ) 
        CALL mxSetField( plhs( 2 ), 1, 'opsa', opsa_pr )
        opse_pr = galmxCreateReal( ) 
        CALL mxSetField( plhs( 2 ), 1, 'opse', opse_pr )
        opsb_pr = galmxCreateReal( ) 
        CALL mxSetField( plhs( 2 ), 1, 'opsb', opsb_pr )
        maxchange_pr = galmxCreateReal( ) 
        CALL mxSetField(  plhs( 2 ), 1, 'maxchange', maxchange_pr )

!  Initialise the components (not needed here, for example only)

!       CALL galmxCopyIntegerToPtr( -1, mxGetPr( flag_pr ), 1 )     
!       CALL galmxCopyRealToPtr( -1.0_wp, mxGetPr( opsa_pr ), 1 )     

!  Factorization phase

        IF ( TRIM( mode ) == 'factor' .OR. 
     &       TRIM( mode ) == 'all' ) THEN

!  Check to ensure the input is a number

          a_in = prhs( a_arg )
          IF ( mxIsNumeric( a_in ) == 0 ) THEN
             CALL mexErrMsgTxt( ' There must be a matrix A ' )
          END IF

!  Get the row and column dimensions and the size of the input matrix

          A%m = mxGetM( a_in )
          A%n = mxGetN( a_in )
          IF ( mxIsSparse( a_in ) ) THEN
            A%ne = mxGetNzmax( a_in )
          ELSE
            A%ne = A%m * A%n
          END IF

!  Allocate space for the input matrix

          ALLOCATE( A%row( A%ne ), A%col( A%ne ), A%val( A%ne ), 
     &              STAT = info )

!  Set the row and column indices if the matrix is sparse

          IF ( mxIsSparse( a_in ) ) THEN

            ALLOCATE( col_ptr( A%n + 1 ) )
!           CALL SMT_put( A%type, 'COORDINATE', stat )

!  Copy the integer components of A into co-ordinate form.
!  ** N.B. indices start at 0 in C **

            cpr_in = mxGetJc( a_in )
            row_in = mxGetIr( a_in )

            CALL galmxCopyPtrToInteger( cpr_in, col_ptr, A%n + 1 )
            CALL galmxCopyPtrToInteger( row_in, A%row, A%ne )

            col_ptr = col_ptr + 1
            A%row = A%row + 1
            DO i = 1, A%n
              DO j = col_ptr( i ), col_ptr( i + 1 ) - 1
                A%col( j ) = i
              END DO
            END DO
            DEALLOCATE( col_ptr )

!  Set the row and column indices if the matrix is dense

          ELSE
            l = 0
            DO j = 1, A%n
              DO i = 1, A%m
                l = l + 1
                A%row( l ) = i
                A%col( l ) = j
              END DO
            END DO
          END IF

!  Copy the real components of A

          val_in = mxGetPr( a_in )
          CALL galmxCopyPtrToReal( val_in, A%val, A%ne )

!  Now remove the lower triangle

          l = 0
          DO k = 1, A%ne
            i = A%row( k )
            j = A%col( k )
            IF ( i <= j ) THEN
              l = l + 1
              A%row( l ) = i
              A%col( l ) = j
              A%val( l ) = A%val( k )
            END IF
          END DO
          A%ne = l

!  Analyse

           CALL SILS_ANALYSE( A, FACTORS, CONTROL, AINFO )

!  set inform components

          CALL galmxCopyIntegerToPtr( AINFO%flag, 
     *                                mxGetPr( flag_pr ), 1 )     
          CALL galmxCopyIntegerToPtr( AINFO%more, 
     *                                mxGetPr( more_pr ), 1 )     
          CALL galmxCopyIntegerToPtr( AINFO%stat, 
     *                                mxGetPr( stat_pr ), 1 )     
          CALL galmxCopyIntegerToPtr( AINFO%nsteps, 
     *                                mxGetPr( nsteps_pr ), 1 )
          CALL galmxCopyIntegerToPtr( AINFO%nrltot, 
     *                                mxGetPr( nrltot_pr ), 1 )
          CALL galmxCopyIntegerToPtr( AINFO%nirtot, 
     *                                mxGetPr( nirtot_pr ), 1 )
          CALL galmxCopyIntegerToPtr( AINFO%nrlnec, 
     *                                mxGetPr( nrlnec_pr ), 1 )
          CALL galmxCopyIntegerToPtr( AINFO%nirnec, 
     *                                mxGetPr( nirnec_pr ), 1 )
          CALL galmxCopyIntegerToPtr( AINFO%nrladu, 
     *                                mxGetPr( nrladu_pr ), 1 )
          CALL galmxCopyIntegerToPtr( AINFO%niradu, 
     *                                mxGetPr( niradu_pr ), 1 )
          CALL galmxCopyIntegerToPtr( AINFO%ncmpa, 
     *                                mxGetPr( ncmpa_pr ), 1 )
          CALL galmxCopyIntegerToPtr( AINFO%nfaulty, 
     *                                mxGetPr( nfaulty_pr ), 1 )
          CALL galmxCopyIntegerToPtr( AINFO%oor, 
     *                                mxGetPr( oor_pr ), 1 )
          CALL galmxCopyIntegerToPtr( AINFO%dup, 
     *                                mxGetPr( dup_pr ), 1 )
          CALL galmxCopyIntegerToPtr( AINFO%maxfrt, 
     *                                mxGetPr( maxfrt_pr ), 1 )
          CALL galmxCopyRealToPtr( AINFO%opsa, 
     *                                mxGetPr( opsa_pr ), 1 )     
          CALL galmxCopyRealToPtr( AINFO%opse, 
     *                                mxGetPr( opse_pr ), 1 )     

!  Check for errors

          IF ( AINFO%FLAG < 0 ) THEN
            CALL mexErrMsgTxt( ' Call to SILS_analyse failed ' )
          END IF

!  Factorize

          CALL SILS_FACTORIZE( A, FACTORS, CONTROL, FINFO )

!  set inform components

          CALL galmxCopyIntegerToPtr( FINFO%flag, 
     *                                mxGetPr( flag_pr ), 1 )     
          CALL galmxCopyIntegerToPtr( FINFO%more, 
     *                                mxGetPr( more_pr ), 1 )     
          CALL galmxCopyIntegerToPtr( FINFO%stat, 
     *                                mxGetPr( stat_pr ), 1 )     
          CALL galmxCopyIntegerToPtr( FINFO%maxfrt, 
     *                                mxGetPr( maxfrt_pr ), 1 )
          CALL galmxCopyIntegerToPtr( FINFO%nebdu, 
     *                                mxGetPr( nebdu_pr ), 1 )
          CALL galmxCopyIntegerToPtr( FINFO%nrlbdu,
     *                                mxGetPr( nrlbdu_pr ), 1 )
          CALL galmxCopyIntegerToPtr( FINFO%nirbdu, 
     *                                mxGetPr( nirbdu_pr ), 1 )
          CALL galmxCopyIntegerToPtr( FINFO%nrltot, 
     *                                mxGetPr( nrltot_pr ), 1 )
          CALL galmxCopyIntegerToPtr( FINFO%nirtot, 
     *                                mxGetPr( nirtot_pr ), 1 )
          CALL galmxCopyIntegerToPtr( FINFO%nrlnec, 
     *                                mxGetPr( nrlnec_pr ), 1 )
          CALL galmxCopyIntegerToPtr( FINFO%nirnec, 
     *                                mxGetPr( nirnec_pr ), 1 )
          CALL galmxCopyIntegerToPtr( FINFO%ncmpbr, 
     *                                mxGetPr( ncmpbr_pr ), 1 )
          CALL galmxCopyIntegerToPtr( FINFO%ncmpbi, 
     *                                mxGetPr( ncmpbi_pr ), 1 )
          CALL galmxCopyIntegerToPtr( FINFO%ntwo, 
     *                                mxGetPr( ntwo_pr ), 1 )
          CALL galmxCopyIntegerToPtr( FINFO%neig, 
     *                                mxGetPr( neig_pr ), 1 )
          CALL galmxCopyIntegerToPtr( FINFO%nfaulty, 
     *                                mxGetPr( nfaulty_pr ), 1 )
          CALL galmxCopyIntegerToPtr( FINFO%delay, 
     *                                mxGetPr( delay_pr ), 1 )
          CALL galmxCopyIntegerToPtr( FINFO%signc, 
     *                                mxGetPr( signc_pr ), 1 )
          CALL galmxCopyIntegerToPtr( FINFO%modstep, 
     *                                mxGetPr( modstep_pr ), 1 )
          CALL galmxCopyIntegerToPtr( FINFO%rank, 
     *                                mxGetPr( rank_pr ), 1 )
          CALL galmxCopyIntegerToPtr( FINFO%step, 
     *                                mxGetPr( step_pr ), 1 )
          CALL galmxCopyRealToPtr( FINFO%opsa, 
     *                             mxGetPr( opsa_pr ), 1 )     
          CALL galmxCopyRealToPtr( FINFO%opse, 
     *                             mxGetPr( opse_pr ), 1 )     
          CALL galmxCopyRealToPtr( FINFO%opsb, 
     *                             mxGetPr( opsb_pr ), 1 )     
          CALL galmxCopyRealToPtr( FINFO%maxchange, 
     *                             mxGetPr( maxchange_pr ), 1 )     

!  Check for errors

          IF ( FINFO%FLAG < 0 ) THEN
            CALL mexErrMsgTxt( ' Call to SILS_factorize failed ' )
          END IF
        END IF

!  Solve phase

        IF ( TRIM( mode ) == 'solve' .OR. 
     &       TRIM( mode ) == 'all' ) THEN

!  Check to ensure the input is a number

          b_in = prhs( b_arg )
          IF ( mxIsNumeric( b_in ) == 0 ) THEN
             CALL mexErrMsgTxt( ' There must be a right-hand-side b ' )
          END IF

 !  Allocate space for the right-hand side and solution

          n = mxGetM( b_in )
          nb = mxGetN( b_in )
          rhs_in = mxGetPr( b_in )

!  one right-hand side

          IF ( nb == 1 ) THEN
            ALLOCATE( B( n ), X( n ) )
            CALL galmxCopyPtrToReal( rhs_in, B, n )

!  Solve without refinement

            X = B
            CALL SILS_SOLVE( A, FACTORS, X, CONTROL, SINFO )

!  Perform one refinement

            CALL SILS_SOLVE( A, FACTORS, X, CONTROL, SINFO, B )

!  multiple right-hand sides

          ELSE
            ALLOCATE( B2( n, nb ), X2( n, nb ) )
            CALL galmxCopyPtrToReal( rhs_in, B2, n * nb )

!  Solve without refinement

            X2 = B2
            CALL SILS_SOLVE( A, FACTORS, X2, CONTROL, SINFO )

!  Perform one refinement

!           CALL SILS_SOLVE( A, FACTORS, X2, CONTROL, SINFO, B2 )

          END IF

!  Output solution

          plhs( 1 ) = galmxCreateRealMatrix( n, nb )
          x_pr = mxGetPr( plhs( 1 ) )
          IF ( nb == 1 ) THEN
            CALL galmxCopyRealToPtr( X, x_pr, n )     
          ELSE
            CALL galmxCopyRealToPtr( X2, x_pr, n * nb )     
          END IF

!  set inform components
          SINFO%flag = nb
          CALL galmxCopyIntegerToPtr( SINFO%flag, 
     *                                mxGetPr( flag_pr ), 1 )     
          CALL galmxCopyIntegerToPtr( SINFO%more, 
     *                                mxGetPr( more_pr ), 1 )     
          CALL galmxCopyIntegerToPtr( SINFO%stat, 
     *                                mxGetPr( stat_pr ), 1 )     
        END IF
      END IF

!  all components now set

      IF ( TRIM( mode ) == 'final' .OR. 
     &     TRIM( mode ) == 'all' ) THEN
         IF ( ALLOCATED( A%row ) ) DEALLOCATE( A%row, STAT = info )
         IF ( ALLOCATED( A%col ) ) DEALLOCATE( A%col, STAT = info )
         IF ( ALLOCATED( A%val ) ) DEALLOCATE( A%val, STAT = info )
         CALL SILS_finalize( FACTORS, CONTROL, info )
      END IF

      RETURN
      END SUBROUTINE mexFunction


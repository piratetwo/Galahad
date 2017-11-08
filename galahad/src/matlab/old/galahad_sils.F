#include <fintrf.h>

! *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
!                 MEX INTERFACE TO GALAHAD_SILS
!
! *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
!  Given a symmetric n by n matrix A and an n-vector b or an n by r 
!  matrix B, solve the system A x = b or the system AX=B. The matrix 
!  A need not be definite. Advantage is taken of sparse A. Options
!  are provided to factorize a matrix A without solving the system, 
!  and to solve systems using previously-determined factors.
!
!  Simple usage -
!
!  to solve a system Ax=b or AX=B
!   [ x, inform ] = galahad_sils( A, b, control )
!
!  Sophisticated usage -
!
!  to initialize data structures prior to factorization
!   galahad_sils( 'initial' )
!
!  to factorize A
!   [ inform ] = galahad_sils( 'factor', A, control )
!
!  to solve Ax=b or AX=B using existing factors
!   [ x, inform ] = galahad_sils( 'solve', b )
!
!  to remove data structures after solution
!   galahad_sils( 'final' )
!
!  Usual Input -
!    A: the symmetric matrix A
!    b a column vector b or matrix of right-hand sides B
!
!  Optional Input -
!    control, a structure containing control parameters.
!      The components are of the form control.value, where
!      value is the name of the corresponding component of
!      the derived type SILS_CONTROL as described in the 
!      manual for the fortran 90 package GALAHARS_SILS.
!      See: http://galahad.rl.ac.uk/galahad-www/doc/sils.pdf
!
!  Usual Output -
!   x: the vector of solutions to Ax=b or matrix of solutions to AX=B
!
!  Optional Output -
!    inform: a structure containing information parameters
!      The components are of the form inform.value, where
!      value is the name of the corresponding component of
!      the derived type SILS_AINFO/FINFO/SINFO as described 
!      in the manual for the fortran 90 package GALAHARS_SILS.
!      See: http://galahad.rl.ac.uk/galahad-www/doc/sils.pdf
!
! *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   originally released with GALAHAD Version 2.1. July 4th 2007

!  For full documentation, see 
!   http://galahad.rl.ac.uk/galahad-www/specs.html

      SUBROUTINE mexFunction( nlhs, plhs, nrhs, prhs )
      USE GALAHAD_MATLAB
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
      REAL ( KIND = wp ) :: mxGetScalar 
      mwSize :: mxGetString
      mwSize :: mxGetM, mxGetN, mxGetIr, mxGetJc
      mwSize :: mxIsNumeric, mxGetNzmax, mxGetNumberOfFields
      mwPointer :: mxCreateStructMatrix, mxGetPr, mxGetField

! -----------------------------------------------------------------------

!  local variables

      mwSize :: i, j, k, l, n, nb, nfields, info
      mwSize :: a_in, b_in, c_in, a_arg, b_arg, c_arg, x_arg, i_arg
      mwPointer :: val_in, rhs_in, cpr_in, row_in, cn_in, x_pr
      mwPointer :: flag_pr, more_pr, stat_pr 
      mwPointer :: nsteps_pr, nrltot_pr, nirtot_pr, nrlnec_pr
      mwPointer :: nirnec_pr, nrladu_pr, niradu_pr, ncmpa_pr
      mwPointer :: nfaulty_pr, oor_pr, dup_pr, maxfrt_pr
      mwPointer :: nebdu_pr, nrlbdu_pr, nirbdu_pr, ncmpbr_pr
      mwPointer :: ncmpbi_pr, ntwo_pr, neig_pr, delay_pr
      mwPointer :: signc_pr, modstep_pr, rank_pr, step_pr
      mwPointer :: opsa_pr, opse_pr, opsb_pr, maxchange_pr

      CHARACTER ( len = 80 ) :: output_unit, filename
      LOGICAL :: filexx, opened, initial_set = .FALSE.
      INTEGER :: iores

      CHARACTER ( len = 7 ) :: mode
      REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: B, X
      REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : , : ) :: B2, X2

      mwSize, PARAMETER :: ninform = 31
      CHARACTER ( LEN = 9 ), PARAMETER :: finform( ninform ) = (/ 
     &     'flag     ', 'more     ', 'stat     ',
     &     'nsteps   ', 'nrltot   ', 'nirtot   ', 'nrlnec   ',
     &     'nirnec   ', 'nrladu   ', 'niradu   ', 'ncmpa    ',
     &     'nfaulty  ', 'oor      ', 'dup      ', 'maxfrt   ',
     &     'nebdu    ', 'nrlbdu   ', 'nirbdu   ', 'ncmpbr   ',
     &     'ncmpbi   ', 'ntwo     ', 'neig     ', 'delay    ',
     &     'signc    ', 'modstep  ', 'rank     ', 'step     ',
     &     'opsa     ', 'opse     ', 'opsb     ', 'maxchange'
     &   /)

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
        CALL mexErrMsgTxt( ' sils requires at least 1 input argument' )
      END IF

      IF ( mxIsChar( prhs( 1 ) ) ) THEN
        i = mxGetString( prhs( 1 ), mode, 7 )
        IF ( .NOT. ( TRIM( mode ) == 'initial' .OR. 
     &               TRIM( mode ) == 'final' ) ) THEN
          IF ( nrhs < 2 ) THEN
            CALL mexErrMsgTxt( ' Too few input arguments to sils' )
          END IF
          IF ( TRIM( mode ) == 'factor' ) THEN
            a_arg = 2
            c_arg = 3
            i_arg = 1
          ELSE IF ( TRIM( mode ) == 'solve' ) THEN
            b_arg = 2
            c_arg = 3
            x_arg = 1
            i_arg = 2
!         ELSE IF ( TRIM( mode ) == 'all' ) THEN
          ELSE
            a_arg = 2
            b_arg = 3
            c_arg = 4
            x_arg = 1
            i_arg = 2
          END IF
          IF ( nrhs > c_arg ) THEN
            CALL mexErrMsgTxt( ' Too many input arguments to sils' )
          END IF
        END IF
      ELSE
!       CALL mexErrMsgTxt( ' First input to sils must be a string' )
        mode = 'all'
        IF ( nrhs < 2 ) THEN
          CALL mexErrMsgTxt( ' Too few input arguments to sils' )
        END IF
        a_arg = 1
        b_arg = 2
        c_arg = 3
        x_arg = 1
        i_arg = 2
        IF ( nrhs > c_arg ) THEN
          CALL mexErrMsgTxt( ' Too many input arguments to sils' )
        END IF
      END IF

      IF ( nlhs > 2 ) THEN
        CALL mexErrMsgTxt( 
     &    ' sils provides at most 2 output arguments' )
      END IF

!  Initialize the internal structures for sils

      IF ( TRIM( mode ) == 'initial' .OR. 
     &     TRIM( mode ) == 'all' ) THEN
        initial_set = .TRUE.
        CALL SILS_INITIALIZE( FACTORS, CONTROL )
        IF ( TRIM( mode ) == 'initial' ) RETURN
      END IF

      IF ( .NOT. TRIM( mode ) == 'final' ) THEN

!  Check that SILS_initialize has been called 

        IF ( .NOT. initial_set ) THEN
          CALL mexErrMsgTxt( ' "initial" must be called first' )
        END IF

!  If the third argument is present, extract the input control data

        IF ( nrhs == c_arg ) THEN
          c_in = prhs( c_arg )
          IF ( .NOT. mxIsStruct( c_in ) ) THEN
            CALL mexErrMsgTxt( 
     &        ' last input argument must be a structure' )
          END IF
          nfields = mxGetNumberOfFields( c_in )
          DO i = 1, nfields
            fi = mxGetFieldNameByNumber( c_in, i )
            SELECT CASE ( TRIM( fi ) )
            CASE( 'lp' )
               cn_in = mxGetField( c_in, 1, 'lp' )
               CONTROL%lp = INT( mxGetScalar( cn_in ) )
            CASE( 'wp' )
              cn_in = mxGetField( c_in, 1, 'wp' )
              CONTROL%wp = INT( mxGetScalar( cn_in ) )
            CASE( 'mp' )
              cn_in = mxGetField( c_in, 1, 'mp' )
              CONTROL%mp = INT( mxGetScalar( cn_in ) )
            CASE( 'sp' )
              cn_in = mxGetField( c_in, 1, 'sp' )
              CONTROL%sp = INT( mxGetScalar( cn_in ) )
            CASE( 'ldiag' )
              cn_in = mxGetField( c_in, 1, 'ldiag' )
              CONTROL%ldiag = INT( mxGetScalar( cn_in ) )
            CASE( 'factorblocking' )
              cn_in = mxGetField( c_in, 1, 'factorblocking' )
              CONTROL%factorblocking = INT( mxGetScalar( cn_in ) )
            CASE( 'solveblocking' )
              cn_in = mxGetField( c_in, 1, 'solveblocking' )
              CONTROL%solveblocking = INT( mxGetScalar( cn_in ) )
            CASE( 'la' )
              cn_in = mxGetField( c_in, 1, 'la' )
              CONTROL%la = INT( mxGetScalar( cn_in ) )
            CASE( 'liw' )
              cn_in = mxGetField( c_in, 1, 'liw' )
              CONTROL%liw = INT( mxGetScalar( cn_in ) )
            CASE( 'maxla' )
              cn_in = mxGetField( c_in, 1, 'maxla' )
              CONTROL%maxla = INT( mxGetScalar( cn_in ) )
            CASE( 'maxliw' )
              cn_in = mxGetField( c_in, 1, 'maxliw' )
              CONTROL%maxliw = INT( mxGetScalar( cn_in ) )
            CASE( 'pivoting' )
              cn_in = mxGetField( c_in, 1, 'pivoting' )
              CONTROL%pivoting = INT( mxGetScalar( cn_in ) )
            CASE( 'thresh' )
              cn_in = mxGetField( c_in, 1, 'thresh' )
              CONTROL%thresh = INT( mxGetScalar( cn_in ) )
            CASE( 'ordering' )
              cn_in = mxGetField( c_in, 1, 'ordering' )
              CONTROL%ordering = INT( mxGetScalar( cn_in ) )
            CASE( 'scaling' )
              cn_in = mxGetField( c_in, 1, 'scaling' )
              CONTROL%scaling = INT( mxGetScalar( cn_in ) )
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

!  No i/o allowed from mex interface?

!        CONTROL%lp = 0
!        CONTROL%wp = 0
!        CONTROL%mp = 0
!        CONTROL%sp = 0

!  Open i/o units

        IF ( CONTROL%lp > 0 ) THEN
          WRITE( output_unit, "( I0 )" ) CONTROL%lp
          filename = "output_sils." // TRIM( output_unit ) 
          INQUIRE( FILE = filename, EXIST = filexx )
          IF ( filexx ) THEN
             OPEN( CONTROL%lp, FILE = filename, FORM = 'FORMATTED',
     &              STATUS = 'OLD', IOSTAT = iores )
          ELSE
             OPEN( CONTROL%lp, FILE = filename, FORM = 'FORMATTED',
     &               STATUS = 'NEW', IOSTAT = iores )
          END IF
        END IF

        IF ( CONTROL%wp > 0 ) THEN
          INQUIRE( CONTROL%wp, OPENED = opened )
          IF ( .NOT. opened ) THEN
            WRITE( output_unit, "( I0 )" ) CONTROL%wp
            filename = "output_sils." // TRIM( output_unit ) 
            INQUIRE( FILE = filename, EXIST = filexx )
            IF ( filexx ) THEN
               OPEN( CONTROL%wp, FILE = filename, FORM = 'FORMATTED',
     &                STATUS = 'OLD', IOSTAT = iores )
            ELSE
               OPEN( CONTROL%wp, FILE = filename, FORM = 'FORMATTED',
     &                 STATUS = 'NEW', IOSTAT = iores )
            END IF
          END IF
        END IF

        IF ( CONTROL%mp > 0 ) THEN
          INQUIRE( CONTROL%mp, OPENED = opened )
          IF ( .NOT. opened ) THEN
            WRITE( output_unit, "( I0 )" ) CONTROL%mp
            filename = "output_sils." // TRIM( output_unit ) 
            INQUIRE( FILE = filename, EXIST = filexx )
            IF ( filexx ) THEN
               OPEN( CONTROL%mp, FILE = filename, FORM = 'FORMATTED',
     &                STATUS = 'OLD', IOSTAT = iores )
            ELSE
               OPEN( CONTROL%mp, FILE = filename, FORM = 'FORMATTED',
     &                 STATUS = 'NEW', IOSTAT = iores )
            END IF
          END IF
        END IF

        IF ( CONTROL%sp > 0 ) THEN
          INQUIRE( CONTROL%sp, OPENED = opened )
          IF ( .NOT. opened ) THEN
            WRITE( output_unit, "( I0 )" ) CONTROL%sp
            filename = "output_sils." // TRIM( output_unit ) 
            INQUIRE( FILE = filename, EXIST = filexx )
            IF ( filexx ) THEN
               OPEN( CONTROL%sp, FILE = filename, FORM = 'FORMATTED',
     &                STATUS = 'OLD', IOSTAT = iores )
            ELSE
               OPEN( CONTROL%sp, FILE = filename, FORM = 'FORMATTED',
     &                 STATUS = 'NEW', IOSTAT = iores )
            END IF
          END IF
        END IF

!  Create inform output structure

        plhs( i_arg ) 
     &    = mxCreateStructMatrix( 1, 1, ninform, finform )

!  Define the components of the structure

        CALL MATLAB_create_integer_component( plhs( i_arg ), 
     &    'flag', flag_pr )
        CALL MATLAB_create_integer_component( plhs( i_arg ), 
     &    'more', more_pr )
        CALL MATLAB_create_integer_component( plhs( i_arg ), 
     &    'stat', stat_pr )
        CALL MATLAB_create_integer_component( plhs( i_arg ), 
     &    'nsteps', nsteps_pr )
        CALL MATLAB_create_integer_component( plhs( i_arg ), 
     &    'nrltot', nrltot_pr )
        CALL MATLAB_create_integer_component( plhs( i_arg ), 
     &    'nirtot', nirtot_pr )
        CALL MATLAB_create_integer_component( plhs( i_arg ), 
     &    'nrlnec', nrlnec_pr )
        CALL MATLAB_create_integer_component( plhs( i_arg ), 
     &    'nirnec', nirnec_pr )
        CALL MATLAB_create_integer_component( plhs( i_arg ), 
     &    'nrladu', nrladu_pr )
        CALL MATLAB_create_integer_component( plhs( i_arg ), 
     &    'niradu', niradu_pr )
        CALL MATLAB_create_integer_component( plhs( i_arg ), 
     &    'ncmpa  ', ncmpa  _pr )
        CALL MATLAB_create_integer_component( plhs( i_arg ), 
     &    'nfaulty', nfaulty_pr )
        CALL MATLAB_create_integer_component( plhs( i_arg ), 
     &    'oor', oor_pr )
        CALL MATLAB_create_integer_component( plhs( i_arg ), 
     &    'dup', dup_pr )
        CALL MATLAB_create_integer_component( plhs( i_arg ), 
     &    'maxfrt', maxfrt_pr )
        CALL MATLAB_create_integer_component( plhs( i_arg ), 
     &    'nebdu', nebdu_pr )
        CALL MATLAB_create_integer_component( plhs( i_arg ), 
     &    'nrlbdu', nrlbdu_pr )
        CALL MATLAB_create_integer_component( plhs( i_arg ), 
     &    'nirbdu', nirbdu_pr )
        CALL MATLAB_create_integer_component( plhs( i_arg ), 
     &    'ncmpbr', ncmpbr_pr )
        CALL MATLAB_create_integer_component( plhs( i_arg ), 
     &    'ncmpbi', ncmpbi_pr )
        CALL MATLAB_create_integer_component( plhs( i_arg ), 
     &    'ntwo', ntwo_pr )
        CALL MATLAB_create_integer_component( plhs( i_arg ), 
     &    'neig', neig_pr )
        CALL MATLAB_create_integer_component( plhs( i_arg ), 
     &    'delay', delay_pr )
        CALL MATLAB_create_integer_component( plhs( i_arg ), 
     &    'signc', signc_pr )
        CALL MATLAB_create_integer_component( plhs( i_arg ), 
     &    'modstep', modstep_pr )
        CALL MATLAB_create_integer_component( plhs( i_arg ), 
     &    'rank', rank_pr )
        CALL MATLAB_create_integer_component( plhs( i_arg ), 
     &    'step', step_pr )
        CALL MATLAB_create_real_component( plhs( i_arg ), 
     &    'opsa', opsa_pr )
        CALL MATLAB_create_real_component( plhs( i_arg ), 
     &    'opse', opse_pr )
        CALL MATLAB_create_real_component( plhs( i_arg ), 
     &    'opsb', opsb_pr )
        CALL MATLAB_create_real_component( plhs( i_arg ), 
     &    'maxchange', maxchange_pr )

!  Initialise the components (not needed here, for example only)

!       CALL MATLAB_copy_to_prt( -1, mxGetPr( flag_pr ), 1 )     
!       CALL MATLAB_copy_to_prt( -1.0_wp, mxGetPr( opsa_pr ), 1 )     

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

!  Copy the integer components of A into co-ordinate form.
!  ** N.B. indices start at 0 in C **

            cpr_in = mxGetJc( a_in )
            row_in = mxGetIr( a_in )

            CALL MATLAB_copy_from_ptr( cpr_in, col_ptr, A%n + 1 )
            CALL MATLAB_copy_from_ptr( row_in, A%row, A%ne )

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
          CALL MATLAB_copy_from_ptr( val_in, A%val, A%ne )

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

          CALL MATLAB_copy_to_prt( AINFO%flag, 
     &                             mxGetPr( flag_pr ) )
          CALL MATLAB_copy_to_prt( AINFO%more, 
     &                             mxGetPr( more_pr ) )
          CALL MATLAB_copy_to_prt( AINFO%stat, 
     &                             mxGetPr( stat_pr ) )
          CALL MATLAB_copy_to_prt( AINFO%nsteps, 
     &                             mxGetPr( nsteps_pr ) )
          CALL MATLAB_copy_to_prt( AINFO%nrltot, 
     &                             mxGetPr( nrltot_pr ) )
          CALL MATLAB_copy_to_prt( AINFO%nirtot, 
     &                             mxGetPr( nirtot_pr ) )
          CALL MATLAB_copy_to_prt( AINFO%nrlnec, 
     &                             mxGetPr( nrlnec_pr ) )
          CALL MATLAB_copy_to_prt( AINFO%nirnec, 
     &                             mxGetPr( nirnec_pr ) )
          CALL MATLAB_copy_to_prt( AINFO%nrladu, 
     &                             mxGetPr( nrladu_pr ) )
          CALL MATLAB_copy_to_prt( AINFO%niradu, 
     &                             mxGetPr( niradu_pr ) )
          CALL MATLAB_copy_to_prt( AINFO%ncmpa, 
     &                             mxGetPr( ncmpa_pr ) )
          CALL MATLAB_copy_to_prt( AINFO%nfaulty, 
     &                             mxGetPr( nfaulty_pr ) )
          CALL MATLAB_copy_to_prt( AINFO%oor, 
     &                             mxGetPr( oor_pr ) )
          CALL MATLAB_copy_to_prt( AINFO%dup, 
     &                             mxGetPr( dup_pr ) )
          CALL MATLAB_copy_to_prt( AINFO%maxfrt, 
     &                             mxGetPr( maxfrt_pr ) )
          CALL MATLAB_copy_to_prt( AINFO%opsa, 
     &                             mxGetPr( opsa_pr ) )
          CALL MATLAB_copy_to_prt( AINFO%opse, 
     &                             mxGetPr( opse_pr ) )

!  Check for errors

          IF ( AINFO%FLAG < 0 ) THEN
            CALL mexErrMsgTxt( ' Call to SILS_analyse failed ' )
          END IF

!  Factorize

          CALL SILS_FACTORIZE( A, FACTORS, CONTROL, FINFO )

!  set inform components

          CALL MATLAB_copy_to_prt( FINFO%flag, 
     &                             mxGetPr( flag_pr ) )
          CALL MATLAB_copy_to_prt( FINFO%more, 
     &                             mxGetPr( more_pr ) )
          CALL MATLAB_copy_to_prt( FINFO%stat, 
     &                             mxGetPr( stat_pr ) )
          CALL MATLAB_copy_to_prt( FINFO%maxfrt, 
     &                             mxGetPr( maxfrt_pr ) )
          CALL MATLAB_copy_to_prt( FINFO%nebdu, 
     &                             mxGetPr( nebdu_pr ) )
          CALL MATLAB_copy_to_prt( FINFO%nrlbdu,
     &                             mxGetPr( nrlbdu_pr ) )
          CALL MATLAB_copy_to_prt( FINFO%nirbdu, 
     &                             mxGetPr( nirbdu_pr ) )
          CALL MATLAB_copy_to_prt( FINFO%nrltot, 
     &                             mxGetPr( nrltot_pr ) )
          CALL MATLAB_copy_to_prt( FINFO%nirtot, 
     &                             mxGetPr( nirtot_pr ) )
          CALL MATLAB_copy_to_prt( FINFO%nrlnec, 
     &                             mxGetPr( nrlnec_pr ) )
          CALL MATLAB_copy_to_prt( FINFO%nirnec, 
     &                             mxGetPr( nirnec_pr ) )
          CALL MATLAB_copy_to_prt( FINFO%ncmpbr, 
     &                             mxGetPr( ncmpbr_pr ) )
          CALL MATLAB_copy_to_prt( FINFO%ncmpbi, 
     &                             mxGetPr( ncmpbi_pr ) )
          CALL MATLAB_copy_to_prt( FINFO%ntwo, 
     &                             mxGetPr( ntwo_pr ) )
          CALL MATLAB_copy_to_prt( FINFO%neig, 
     &                             mxGetPr( neig_pr ) )
          CALL MATLAB_copy_to_prt( FINFO%nfaulty, 
     &                             mxGetPr( nfaulty_pr ) )
          CALL MATLAB_copy_to_prt( FINFO%delay, 
     &                             mxGetPr( delay_pr ) )
          CALL MATLAB_copy_to_prt( FINFO%signc, 
     &                             mxGetPr( signc_pr ) )
          CALL MATLAB_copy_to_prt( FINFO%modstep, 
     &                             mxGetPr( modstep_pr ) )
          CALL MATLAB_copy_to_prt( FINFO%rank, 
     &                             mxGetPr( rank_pr ) )
          CALL MATLAB_copy_to_prt( FINFO%step, 
     &                             mxGetPr( step_pr ) )
          CALL MATLAB_copy_to_prt( FINFO%opsa, 
     &                             mxGetPr( opsa_pr ) )
          CALL MATLAB_copy_to_prt( FINFO%opse, 
     &                             mxGetPr( opse_pr ) )
          CALL MATLAB_copy_to_prt( FINFO%opsb, 
     &                             mxGetPr( opsb_pr ) )
          CALL MATLAB_copy_to_prt( FINFO%maxchange, 
     &                             mxGetPr( maxchange_pr ) )

!  Check for errors

          IF ( FINFO%FLAG < 0 ) THEN
            CALL mexErrMsgTxt( ' Call to SILS_factorize failed ' )
          END IF
        END IF

!  Solve phase

        IF ( TRIM( mode ) == 'solve' .OR. 
     &       TRIM( mode ) == 'all' ) THEN

          IF ( .NOT. ( ALLOCATED( A%row ) .AND.
     &                 ALLOCATED( A%col ) .AND.
     &                 ALLOCATED( A%val ) ) ) THEN
             CALL mexErrMsgTxt( ' There must be existing factors ' )
          END IF

!  Check to ensure the input is a number

          b_in = prhs( b_arg )
          IF ( mxIsNumeric( b_in ) == 0 ) THEN
             CALL mexErrMsgTxt( ' There must be a right-hand-side b ' )
          END IF

 !  Allocate space for the right-hand side and solution

          n = mxGetM( b_in )
          IF ( A%n /= n ) THEN
             CALL mexErrMsgTxt( 
     &         ' A and b/B must have compatible dimensions ' )
          END IF

          nb = mxGetN( b_in )
          rhs_in = mxGetPr( b_in )

!  one right-hand side

          IF ( nb == 1 ) THEN
            ALLOCATE( B( n ), X( n ) )
            CALL MATLAB_copy_from_ptr( rhs_in, B, n )

!  Solve without refinement

            X = B
            CALL SILS_SOLVE( A, FACTORS, X, CONTROL, SINFO )

!  Perform one refinement

            CALL SILS_SOLVE( A, FACTORS, X, CONTROL, SINFO, B )

!  multiple right-hand sides

          ELSE
            ALLOCATE( B2( n, nb ), X2( n, nb ) )
            CALL MATLAB_copy_from_ptr( rhs_in, B2, n, nb )

!  Solve without refinement

            X2 = B2
            CALL SILS_SOLVE( A, FACTORS, X2, CONTROL, SINFO )

!  Perform one refinement

!           CALL SILS_SOLVE( A, FACTORS, X2, CONTROL, SINFO, B2 )

          END IF

!  Output solution

          plhs( x_arg ) = MATLAB_create_real( n, nb )
          x_pr = mxGetPr( plhs( x_arg ) )
          IF ( nb == 1 ) THEN
            CALL MATLAB_copy_to_prt( X, x_pr, n )     
          ELSE
            CALL MATLAB_copy_to_prt( X2, x_pr, n, nb )     
          END IF

!  set inform components

          CALL MATLAB_copy_to_prt( SINFO%flag, 
     &                             mxGetPr( flag_pr ) )
          CALL MATLAB_copy_to_prt( SINFO%more, 
     &                             mxGetPr( more_pr ) )
          CALL MATLAB_copy_to_prt( SINFO%stat, 
     &                             mxGetPr( stat_pr ) )
        END IF
      END IF

!     WRITE( message, * ) 'hello'
!     CALL MEXPRINTF( TRIM( message ) // char( 13 ) )

!  all components now set

      IF ( TRIM( mode ) == 'final' .OR. 
     &     TRIM( mode ) == 'all' ) THEN
         IF ( ALLOCATED( A%row ) ) DEALLOCATE( A%row, STAT = info )
         IF ( ALLOCATED( A%col ) ) DEALLOCATE( A%col, STAT = info )
         IF ( ALLOCATED( A%val ) ) DEALLOCATE( A%val, STAT = info )
         CALL SILS_finalize( FACTORS, CONTROL, info )
      END IF

!  close any opened io units

      IF ( CONTROL%lp > 0 ) THEN
         INQUIRE( CONTROL%lp, OPENED = opened )
         IF ( opened ) CLOSE( control%lp )
      END IF

      IF ( CONTROL%wp > 0 ) THEN
         INQUIRE( CONTROL%wp, OPENED = opened )
         IF ( opened ) CLOSE( control%wp )
      END IF

      IF ( CONTROL%mp > 0 ) THEN
         INQUIRE( CONTROL%mp, OPENED = opened )
         IF ( opened ) CLOSE( control%mp )
      END IF

      IF ( CONTROL%sp > 0 ) THEN
         INQUIRE( CONTROL%sp, OPENED = opened )
         IF ( opened ) CLOSE( control%sp )
      END IF

      RETURN
      END SUBROUTINE mexFunction


#include <fintrf.h>

! *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!                 MEX INTERFACE TO GALAHAD_SILS
!
!  Expected Matlab syntax:
!
!    [X,INFORM] = SILS(A,B)
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
      CHARACTER ( len = 16 ) :: mxGetFieldNameByNumber, fi
      LOGICAL :: mxIsSparse, mxIsStruct
      mwSize :: mxGetM, mxGetN
      mwSize :: mxIsNumeric, mxGetElementSize, mxGetNzmax
      mwSize :: mxGetIr, mxGetJc, mxClassIDFromClassName
      mwsize :: mxGetNumberOfFields, mxGetNumberOfElements
      mwPointer :: plhs( * ), prhs( * )
      mwPointer :: mxCreateStructMatrix, mxGetPr, mxGetField
      mwPointer :: galmxCreateInteger, galmxCreateReal
      mwPointer :: galmxCreateIntegerMatrix, galmxCreateRealMatrix

      mwPointer :: mxCreateString


! -----------------------------------------------------------------------

!  local variables

      mwSize :: i, j, k, l, a_in, b_in, dig, stat, n, nb, info, ii
      mwsize :: nfields, iii( 1 )
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

      mwPointer :: ops_pr
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

      TYPE :: dSILS_control
       REAL ( KIND = wp ) :: CNTL( 5 ) ! MA27 internal real controls
       REAL ( KIND = wp ) :: multiplier ! Factor by which arrays sizes are to
                         ! be increased if they are too small
       REAL ( KIND = wp ) :: reduce ! If previously allocated internal 
                         ! workspace arrays are greater than reduce times 
                         ! the currently required sizes, they are reset to 
                         ! current requirments
       REAL ( KIND = wp ) :: u     ! Pivot threshold
       REAL ( KIND = wp ) :: static_tolerance ! used for setting static 
                           ! pivot level                                    NEW
       REAL ( KIND = wp ) :: static_level ! used for switch to static       NEW
       REAL ( KIND = wp ) :: tolerance ! Anything less than this is 
                         ! considered zero
       INTEGER :: ICNTL( 30 ) ! MA27 internal integer controls
       INTEGER :: lp     ! Unit for error messages
       INTEGER :: wp     ! Unit for warning messages
       INTEGER :: mp     ! Unit for monitor output                          NEW
       INTEGER :: sp     ! Unit for statistical output                      NEW
       INTEGER :: ldiag  ! Controls level of diagnostic output
       INTEGER :: factorblocking ! Level 3 blocking in factorize            NEW
       INTEGER :: solveblocking ! Level 2 and 3 blocking in solve           NEW
       INTEGER :: la     ! Initial size for real array for the factors.
                         ! If less than nrlnec, default size used.  
       INTEGER :: liw    ! Initial size for integer array for the factors.   
                         ! If less than nirnec, default size used.  
       INTEGER :: maxla  ! Max. size for real array for the factors.
       INTEGER :: maxliw ! Max. size for integer array for the factors.
       INTEGER :: pivoting  ! Controls pivoting:
                  !  1  Numerical pivoting will be performed. 
                  !  2  No pivoting will be performed and an error exit will   
                  !     occur immediately a pivot sign change is detected.
                  !  3  No pivoting will be performed and an error exit will 
                  !     occur if a zero pivot is detected. 
                  !  4  No pivoting is performed but pivots are changed to
                  !     all be positive.
       INTEGER :: thresh ! Controls threshold for detecting full rows in 
                  !     analyse, registered as percentage of N
                  ! 100 Only fully dense rows detected (default)            NEW
       INTEGER :: ordering  ! Controls ordering:                            NEW
                 !  0  AMD using MC47
                 !  1  User defined
                 !  2  AMD using MC50
                 !  3  Min deg as in MA57
                 !  4  Metis_nodend ordering
                 !  5  Ordering chosen depending on matrix characteristics.
                 !     At the moment choices are MC50 or Metis_nodend
                 ! >5  Presently equivalent to 5 but may chnage
       INTEGER :: scaling  ! Controls scaling:                              NEW
                 !  0  No scaling
                 ! >0  Scaling using MC64 but may change for > 1
      END TYPE dSILS_control



!       REAL ( KIND = wp ) :: 
!multiplier
!reduce
!u
!static_tolerance
!static_level
!tolerance
!       INTEGER :: 
!lp
!wp
!mp
!sp
!ldiag
!factorblocking
!solveblocking
!la
!liw
!maxla
!maxliw
!pivoting
!thresh
!ordering
!scaling

!static_tolerance

character ( len = 16 ) ::


      mwSize :: status
      REAL ( KIND = wp ) :: ops

!  arguments for SILS

      TYPE ( SMT_type ) :: A
      TYPE ( SILS_control ) :: CONTROL
      TYPE ( SILS_ainfo ) :: AINFO
      TYPE ( SILS_finfo ) :: FINFO
      TYPE ( SILS_sinfo ) :: SINFO
      TYPE ( SILS_factors ) :: FACTORS

      mwPointer, ALLOCATABLE :: col_ptr( : )

      IF ( nrhs < 2 ) THEN
        CALL mexErrMsgTxt( ' sils requires at least 2 input arguments' )
      ELSE IF ( nrhs > 3 ) THEN
        CALL mexErrMsgTxt( ' sils requires at most 3 input arguments' )
!     ELSE IF ( nlhs > 2 ) THEN
      ELSE IF ( nlhs > 3 ) THEN
        CALL mexErrMsgTxt( ' sils provides at most 2 output arguments' )
      END IF

!  Create space for and fill return arguments to MATLAB.

!  Check to ensure the input is a number

      a_in = prhs( 1 )
      IF ( mxIsNumeric( a_in ) == 0 ) THEN
         CALL mexErrMsgTxt( ' There must be a matrix A ' )
      END IF

      b_in = prhs( 2 )
      IF ( mxIsNumeric( b_in ) == 0 ) THEN
         CALL mexErrMsgTxt( ' There must be a right-hand-side b ' )
      END IF

!  Initialize the internal structures

      CALL SILS_INITIALIZE( FACTORS, CONTROL )

!  If the third argument is present, extract the input cointrol data

      ii = 0
      IF ( nrhs == 3 ) THEN
        IF ( .NOT. mxIsStruct( prhs( 3 ) ) ) THEN
          CALL mexErrMsgTxt( ' 3rd input argument must be a structure' )
        END IF
        nfields = mxGetNumberOfFields( prhs( 3 ) )
        DO i = 1, nfields
 !      DO i = 1, 
          fi = mxGetFieldNameByNumber( prhs( 3 ), i )
          SELECT CASE ( TRIM( fi ) )
          CASE( 'lp' )
             ii = 9999
             cn_in = mxGetField( prhs( 3 ), 1, 'lp' )
             ii = cn_in
!            CALL galmxCopyPtrToInteger( cn_in, iii, 1 )
             CALL mxCopyPtrToInteger8( cn_in, iii, 1 )
              ii = mxGetNumberOfElements(cn_in)
             CALL galmxCopyPtrToReal( cn_in, att, 1 )
              at = att(1)
             CALL mxGetdata( cn_in, iii, 1 )
             at = mxGetScalar(  cn_in )
             ii = at
          CASE( 'wp' )
            cn_in = mxGetField( prhs( 3 ), 1, 'wp' )
            CALL galmxCopyPtrToInteger( cn_in, CONTROL%wp, 1 )
        CASE( 'mp' )
            cn_in = mxGetField( prhs( 3 ), 1, 'mp' )
            CALL galmxCopyPtrToInteger( cn_in, CONTROL%mp, 1 )
          CASE( 'sp' )
            cn_in = mxGetField( prhs( 3 ), 1, 'sp' )
            CALL galmxCopyPtrToInteger( cn_in, CONTROL%sp, 1 )
          CASE( 'ldiag' )
            cn_in = mxGetField( prhs( 3 ), 1, 'ldiag' )
            CALL galmxCopyPtrToInteger( cn_in, CONTROL%ldiag, 1 )
          CASE( 'factorblocking' )
            cn_in = mxGetField( prhs( 3 ), 1, 'factorblocking' )
            CALL galmxCopyPtrToInteger( cn_in, 
     *                                  CONTROL%factorblocking, 1 )
          CASE( 'solveblocking' )
            cn_in = mxGetField( prhs( 3 ), 1, 'solveblocking' )
            CALL galmxCopyPtrToInteger( cn_in, 
     *                                  CONTROL%solveblocking, 1 )
          CASE( 'la' )
            cn_in = mxGetField( prhs( 3 ), 1, 'la' )
            CALL galmxCopyPtrToInteger( cn_in, CONTROL%la, 1 )
          CASE( 'liw' )
            cn_in = mxGetField( prhs( 3 ), 1, 'liw' )
            CALL galmxCopyPtrToInteger( cn_in, CONTROL%liw, 1 )
          CASE( 'maxla' )
            cn_in = mxGetField( prhs( 3 ), 1, 'maxla' )
            CALL galmxCopyPtrToInteger( cn_in, CONTROL%maxla, 1 )
          CASE( 'maxliw' )
            cn_in = mxGetField( prhs( 3 ), 1, 'maxliw' )
            CALL galmxCopyPtrToInteger( cn_in, CONTROL%maxliw, 1 )
          CASE( 'pivoting' )
            cn_in = mxGetField( prhs( 3 ), 1, 'pivoting' )
            CALL galmxCopyPtrToInteger( cn_in, CONTROL%pivoting, 1 )
          CASE( 'thresh' )
            cn_in = mxGetField( prhs( 3 ), 1, 'thresh' )
            CALL galmxCopyPtrToInteger( cn_in, CONTROL%thresh, 1 )
          CASE( 'ordering' )
            cn_in = mxGetField( prhs( 3 ), 1, 'ordering' )
            CALL galmxCopyPtrToInteger( cn_in, CONTROL%ordering, 1 )
          CASE( 'scaling' )
            cn_in = mxGetField( prhs( 3 ), 1, 'scaling' )
            CALL galmxCopyPtrToInteger( cn_in, CONTROL%scaling, 1 )
          CASE( 'multiplier' )
            cn_in = mxGetField( prhs( 3 ), 1, 'multiplier' )
            CALL galmxCopyPtrToReal( cn_in, CONTROL%multiplier, 1 )
          CASE( 'reduce' )
            cn_in = mxGetField( prhs( 3 ), 1, 'reduce' )
            CALL galmxCopyPtrToReal( cn_in, CONTROL%reduce, 1 )
          CASE( 'u' )
            cn_in = mxGetField( prhs( 3 ), 1, 'u' )
            CALL galmxCopyPtrToReal( cn_in, CONTROL%u, 1 )
          CASE( 'static_tolerance' )
            cn_in = mxGetField( prhs( 3 ), 1, 'static_tolerance' )
            CALL galmxCopyPtrToReal( cn_in, 
     *                               CONTROL%static_tolerance, 1 )
          CASE( 'static_level' )
            cn_in = mxGetField( prhs( 3 ), 1, 'static_level' )
            CALL galmxCopyPtrToReal( cn_in, CONTROL%static_level, 1 )
          CASE( 'tolerance' )
            cn_in = mxGetField( prhs( 3 ), 1, 'tolerance' )
            CALL galmxCopyPtrToReal( cn_in, CONTROL%tolerance, 1 )
          END SELECT
        END DO
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

      ALLOCATE( A%row( A%ne ), A%col( A%ne ), A%val( A%ne ) )

!  Set the row and column indices if the matrix is sparse

      IF ( mxIsSparse( a_in ) ) THEN

        ALLOCATE( col_ptr( A%n + 1 ) )
!       CALL SMT_put( A%type, 'COORDINATE', stat )

!  Copy the integer components of A into co-ordinate form.
!  N.B. indices start at 0 in C

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

!  Create inform output structure

      plhs( 2 ) = mxCreateStructMatrix( 1, 1, ninform, finform )

!  Define the components of the structure

      plhs( 3 ) = mxCreateStructMatrix( 1, 1, ninform2, finform2 )
      string_pr = mxCreateString( fi )
      CALL mxSetField( plhs( 3 ), 1, 'string', string_pr )
      CALL mxCopyCharacterToPtr( fi, string_pr, 1 )

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

!  initialise the components (not needed here, for example only)

!      CALL galmxCopyIntegerToPtr( -1, mxGetPr( flag_pr ), 1 )     
!      CALL galmxCopyRealToPtr( -1.0_wp, mxGetPr( opsa_pr ), 1 )     

!  Analyse

      CALL SILS_ANALYSE( A, FACTORS, CONTROL, AINFO )

!  set inform components

      CALL galmxCopyIntegerToPtr( AINFO%flag, 
     *                            mxGetPr( flag_pr ), 1 )     
      CALL galmxCopyIntegerToPtr( AINFO%more, 
     *                            mxGetPr( more_pr ), 1 )     
      CALL galmxCopyIntegerToPtr( AINFO%stat, 
     *                            mxGetPr( stat_pr ), 1 )     
      CALL galmxCopyIntegerToPtr( AINFO%nsteps, 
     *                            mxGetPr( nsteps_pr ), 1 )
      CALL galmxCopyIntegerToPtr( AINFO%nrltot, 
     *                            mxGetPr( nrltot_pr ), 1 )
      CALL galmxCopyIntegerToPtr( AINFO%nirtot, 
     *                            mxGetPr( nirtot_pr ), 1 )
      CALL galmxCopyIntegerToPtr( AINFO%nrlnec, 
     *                            mxGetPr( nrlnec_pr ), 1 )
      CALL galmxCopyIntegerToPtr( AINFO%nirnec, 
     *                            mxGetPr( nirnec_pr ), 1 )
      CALL galmxCopyIntegerToPtr( AINFO%nrladu, 
     *                            mxGetPr( nrladu_pr ), 1 )
      CALL galmxCopyIntegerToPtr( AINFO%niradu, 
     *                            mxGetPr( niradu_pr ), 1 )
      CALL galmxCopyIntegerToPtr( AINFO%ncmpa, 
     *                            mxGetPr( ncmpa_pr ), 1 )
      CALL galmxCopyIntegerToPtr( AINFO%nfaulty, 
     *                            mxGetPr( nfaulty_pr ), 1 )
      CALL galmxCopyIntegerToPtr( AINFO%oor, 
     *                            mxGetPr( oor_pr ), 1 )
      CALL galmxCopyIntegerToPtr( AINFO%dup, 
     *                            mxGetPr( dup_pr ), 1 )
      CALL galmxCopyIntegerToPtr( AINFO%maxfrt, 
     *                            mxGetPr( maxfrt_pr ), 1 )
      CALL galmxCopyRealToPtr( AINFO%opsa, 
     *                            mxGetPr( opsa_pr ), 1 )     
      CALL galmxCopyRealToPtr( AINFO%opse, 
     *                            mxGetPr( opse_pr ), 1 )     

!  Check for errors

      IF ( AINFO%FLAG < 0 ) THEN
        CALL mexErrMsgTxt( ' Call to SILS_analyse failed ' )
      END IF

!  Factorize

      CALL SILS_FACTORIZE( A, FACTORS, CONTROL, FINFO )

!  set inform components

      CALL galmxCopyIntegerToPtr( FINFO%flag, 
     *                            mxGetPr( flag_pr ), 1 )     
      CALL galmxCopyIntegerToPtr( FINFO%more, 
     *                            mxGetPr( more_pr ), 1 )     
      CALL galmxCopyIntegerToPtr( FINFO%stat, 
     *                            mxGetPr( stat_pr ), 1 )     
      CALL galmxCopyIntegerToPtr( FINFO%maxfrt, 
     *                            mxGetPr( maxfrt_pr ), 1 )
      CALL galmxCopyIntegerToPtr( FINFO%nebdu, 
     *                            mxGetPr( nebdu_pr ), 1 )
      CALL galmxCopyIntegerToPtr( FINFO%nrlbdu,
     *                            mxGetPr( nrlbdu_pr ), 1 )
      CALL galmxCopyIntegerToPtr( FINFO%nirbdu, 
     *                            mxGetPr( nirbdu_pr ), 1 )
      CALL galmxCopyIntegerToPtr( FINFO%nrltot, 
     *                            mxGetPr( nrltot_pr ), 1 )
      CALL galmxCopyIntegerToPtr( FINFO%nirtot, 
     *                            mxGetPr( nirtot_pr ), 1 )
      CALL galmxCopyIntegerToPtr( FINFO%nrlnec, 
     *                            mxGetPr( nrlnec_pr ), 1 )
      CALL galmxCopyIntegerToPtr( FINFO%nirnec, 
     *                            mxGetPr( nirnec_pr ), 1 )
      CALL galmxCopyIntegerToPtr( FINFO%ncmpbr, 
     *                            mxGetPr( ncmpbr_pr ), 1 )
      CALL galmxCopyIntegerToPtr( FINFO%ncmpbi, 
     *                            mxGetPr( ncmpbi_pr ), 1 )
      CALL galmxCopyIntegerToPtr( FINFO%ntwo, 
     *                            mxGetPr( ntwo_pr ), 1 )
      CALL galmxCopyIntegerToPtr( FINFO%neig, 
     *                            mxGetPr( neig_pr ), 1 )
      CALL galmxCopyIntegerToPtr( FINFO%nfaulty, 
     *                            mxGetPr( nfaulty_pr ), 1 )
      CALL galmxCopyIntegerToPtr( FINFO%delay, 
     *                            mxGetPr( delay_pr ), 1 )
      CALL galmxCopyIntegerToPtr( FINFO%signc, 
     *                            mxGetPr( signc_pr ), 1 )
      CALL galmxCopyIntegerToPtr( FINFO%modstep, 
     *                            mxGetPr( modstep_pr ), 1 )
      CALL galmxCopyIntegerToPtr( FINFO%rank, 
     *                            mxGetPr( rank_pr ), 1 )
      CALL galmxCopyIntegerToPtr( FINFO%step, 
     *                            mxGetPr( step_pr ), 1 )
      CALL galmxCopyRealToPtr( FINFO%opsa, 
     *                         mxGetPr( opsa_pr ), 1 )     
      CALL galmxCopyRealToPtr( FINFO%opse, 
     *                         mxGetPr( opse_pr ), 1 )     
      CALL galmxCopyRealToPtr( FINFO%opsb, 
     *                         mxGetPr( opsb_pr ), 1 )     
      CALL galmxCopyRealToPtr( FINFO%maxchange, 
     *                         mxGetPr( maxchange_pr ), 1 )     

!  Check for errors

      IF ( FINFO%FLAG < 0 ) THEN
        CALL mexErrMsgTxt( ' Call to SILS_factorize failed ' )
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

!       CALL SILS_SOLVE( A, FACTORS, X2, CONTROL, SINFO, B2 )

      END IF
      DEALLOCATE( A%row, A%col, A%val )

!  Create space for and fill return arguments to MATLAB.

!  Output solution

      plhs( 1 ) = galmxCreateRealMatrix( n, nb )
      x_pr = mxGetPr( plhs( 1 ) )
      IF ( nb == 1 ) THEN
        CALL galmxCopyRealToPtr( X, x_pr, n )     
      ELSE
        CALL galmxCopyRealToPtr( X2, x_pr, n * nb )     
      END IF

!  set inform components

      SINFO%flag = ii
      CALL galmxCopyIntegerToPtr( SINFO%flag, 
     *                            mxGetPr( flag_pr ), 1 )     
      CALL galmxCopyIntegerToPtr( SINFO%more, 
     *                            mxGetPr( more_pr ), 1 )     
      CALL galmxCopyIntegerToPtr( SINFO%stat, 
     *                            mxGetPr( stat_pr ), 1 )     

!  all components now set

      CALL SILS_finalize( FACTORS, CONTROL, info )

      RETURN
      END SUBROUTINE mexFunction


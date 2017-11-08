#include <fintrf.h>

! *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
!    GENERIC MEX INTERFACE (cut-down version of that for HSL_VE12)
!
! *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
!  Given a symmetric n by n matrix H, an m by n matrix A, an n-vector 
!  g, a constant f, n-vectors x_l <= x_u and m-vectors c_l <= c_u, 
!  find a local minimizer of the QUADRATIC PROGRAMMING problem
!    minimize 0.5 * x' * H * x + g' * x + f
!    subject to c_l <= A * x <= c_u and x_l <= x <= x_u
!  H need not be definite. Advantage is taken of sparse A and H. 
!
!  Simple usage -
!
!  to solve the quadratic program
!   [ x, inform, aux ] 
!     = galahad_qpc( H, g, f, A, c_l, c_u, x_l, x_u, control )
!
!  Sophisticated usage -
!
!  to initialize data structures prior to solution
!   galahad_qpc( 'initial' )
!
!  to solve the quadratic program using existing data structures
!   [ x, inform, aux ]
!     = galahad_qpc( H, g, f, A, c_l, c_u, x_l, x_u, control )
!
!  to remove data structures after solution
!   galahad_qpc( 'final' )
!
!  Usual Input -
!    H: the symmetric n by n matrix H
!    g: the n-vector g
!    f: the scalar f
!    A: the m by n matrix A
!    c_l: the m-vector c_l. The value -inf should be used for infinite bounds
!    c_u: the m-vector c_u. The value inf should be used for infinite bounds
!    x_l: the n-vector x_l. The value -inf should be used for infinite bounds
!    x_u: the n-vector x_u. The value inf should be used for infinite bounds
!
!  Optional Input -
!    control, a structure containing control parameters.
!      The components are of the form control.value, where
!      value is the name of the corresponding component of
!      the derived type VE12_CONTROL as described in the 
!      manual for the fortran 90 package HSL_VE12.
!
!  Usual Output -
!   x: a local minimizer
!
!  Optional Output -
!   inform: a structure containing information parameters
!      The components are of the form inform.value, where
!      value is the name of the corresponding component of
!      the derived type VE12_INFORM as described in the spec sheet for 
!      the fortran 90 package HSL_VE12.
!  aux: a structure containing Lagrange multipliers and constraint status
!   aux.c: values of the constraints A * x
!
! *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
      SUBROUTINE mexFunction( nlhs, plhs, nrhs, prhs )
      USE HSL_MATLAB
      USE HSL_VE12_double
      IMPLICIT NONE
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )

! ------------------------- Do not change -------------------------------

!  Keep the above subroutine, argument, and function declarations for use
!  in all your fortran mex files.
!
      INTEGER * 4 :: nlhs, nrhs
      mwPointer :: plhs( * ), prhs( * )

! -----------------------------------------------------------------------

!  Mex functions/variables

      INTEGER, PARAMETER :: slen = 30
      CHARACTER ( LEN = slen ) :: fi, fii, mxGetFieldNameByNumber
      LOGICAL :: mxIsSparse, mxIsChar, mxIsStruct
      REAL ( KIND = wp ) :: mxGetScalar
      mwSize :: mxGetString
      mwSize :: mxGetM, mxGetN, mxGetIr, mxGetJc
      mwSize :: mxIsNumeric, mxGetNzmax, mxGetNumberOfFields
      mwPointer :: mxCreateStructMatrix, mxGetPr, mxGetField

!  local variables

      mwSize :: i, ii, j, k, l, nfields, nnfields, stat, info
      mwSize :: c_in, h_in, g_in, f_in, a_in, cl_in, cu_in, xl_in, xu_in
      mwSize :: g_pr, f_pr, cl_pr, cu_pr, xl_pr, xu_pr
      mwSize :: h_arg, g_arg, f_arg, a_arg, cl_arg, cu_arg
      mwSize :: xl_arg, xu_arg, c_arg, x_arg, i_arg, aux_arg
      mwPointer :: val_pr, cn_in, cnn_in, x_pr, y_pr, z_pr, b_stat_pr
      mwPointer :: c_stat_pr, c_pr, h_cpr_pr, h_row_pr, a_cpr_pr, a_row_pr

      mwPointer :: status_pr, alloc_status_pr, bad_alloc_pr
      mwPointer :: p_found_pr, obj_pr

      mwPointer :: time_pr, time_total_pr, time_preprocess_pr

      CHARACTER ( len = 80 ) :: output_unit, filename
      LOGICAL :: filexx, opened, initial_set = .FALSE.
      INTEGER :: iores

      CHARACTER ( len = 7 ) :: mode
      INTEGER, ALLOCATABLE, DIMENSION( : ) :: C_stat, B_stat

      mwSize, PARAMETER :: ninform = 6
      CHARACTER ( LEN = 21 ), PARAMETER :: finform( ninform ) = (/             &
           'status               ', 'alloc_status         ',                   &
           'bad_alloc            ', 'p_found              ',                   & 
           'obj                  ',  'time                 ' /)
      mwSize, PARAMETER :: t_ninform = 2
      CHARACTER ( LEN = 21 ), PARAMETER :: t_finform( t_ninform ) = (/         &
           'total                ', 'preprocess           ' /)
      mwSize, PARAMETER :: naux = 1
      CHARACTER ( LEN = 6 ), PARAMETER :: faux( naux ) = (/ 'c     ' /)

!  arguments for VE12

      TYPE ( QPT_problem_type ), SAVE :: p
      TYPE ( VE12_control_type ), SAVE :: control
      TYPE ( VE12_inform_type ) :: inform
      TYPE ( VE12_data_type ), SAVE :: data

      mwPointer, ALLOCATABLE :: col_ptr( : )

!  Test input/output arguments

      IF ( nrhs < 1 ) THEN
        CALL mexErrMsgTxt( ' qpc requires at least 1 input argument' )
      END IF

      IF ( mxIsChar( prhs( 1 ) ) ) THEN
        i = mxGetString( prhs( 1 ), mode, 7 )
        IF ( .NOT. ( TRIM( mode ) == 'initial' .OR.                            &
                     TRIM( mode ) == 'final' ) ) THEN
          IF ( nrhs < 2 ) THEN
            CALL mexErrMsgTxt( ' Too few input arguments to qpc' )
          END IF
          h_arg = 2 ; g_arg = 3 ; f_arg = 4 ; a_arg = 5
          cl_arg = 6 ; cu_arg = 7 ; xl_arg = 8 ; xu_arg = 9 ; c_arg = 10
          x_arg = 1 ; i_arg = 2 ; aux_arg = 3
          IF ( nrhs > c_arg ) THEN
            CALL mexErrMsgTxt( ' Too many input arguments to qpc' )
          END IF
        END IF
      ELSE
        mode = 'all'
        IF ( nrhs < 2 ) THEN
          CALL mexErrMsgTxt( ' Too few input arguments to qpc' )
        END IF
        h_arg = 1 ; g_arg = 2 ; f_arg = 3 ; a_arg = 4 ; 
        cl_arg = 5 ; cu_arg = 6 ; xl_arg = 7 ; xu_arg = 8 ; c_arg = 9
        x_arg = 1 ; i_arg = 2 ; aux_arg = 3
        IF ( nrhs > c_arg ) THEN
          CALL mexErrMsgTxt( ' Too many input arguments to qpc' )
        END IF
      END IF

      IF ( nlhs > 3 ) THEN
        CALL mexErrMsgTxt(             &
          ' qpc provides at most 3 output arguments' )
      END IF

!  Initialize the internal structures for qpc

      IF ( TRIM( mode ) == 'initial' .OR. TRIM( mode ) == 'all' ) THEN
        initial_set = .TRUE.
        CALL VE12_initialize( data, control )
        IF ( TRIM( mode ) == 'initial' ) RETURN
      END IF

      IF ( .NOT. TRIM( mode ) == 'final' ) THEN

!  Check that VE12_initialize has been called 

        IF ( .NOT. initial_set ) THEN
          CALL mexErrMsgTxt( ' "initial" must be called first' )
        END IF

!  If the third argument is present, extract the input control data

        IF ( nrhs == c_arg ) THEN
          c_in = prhs( c_arg )
          IF ( .NOT. mxIsStruct( c_in ) ) THEN
            CALL mexErrMsgTxt( ' last input argument must be a structure' )
          END IF
          nfields = mxGetNumberOfFields( c_in )
          DO i = 1, nfields
            fi = mxGetFieldNameByNumber( c_in, i )
            SELECT CASE ( TRIM( fi ) )
            CASE( 'error' )
              CALL MATLAB_get_value( c_in, 'error',                             &
                cn_in, CONTROL%error )
            CASE( 'out' )
              CALL MATLAB_get_value( c_in, 'out',                               &
                cn_in, CONTROL%out )
            CASE( 'print_level' )
              CALL MATLAB_get_value( c_in, 'print_level',                       &
                cn_in, CONTROL%print_level )
            CASE( 'prefix' )
              CALL MATLAB_get_value( c_in, 'prefix',                            &
                cn_in, CONTROL%prefix, slen )
            END SELECT
          END DO
        END IF

!  Open i/o units

        IF ( control%error > 0 ) THEN
          WRITE( output_unit, "( I0 )" ) control%error
          filename = "output_qpc." // TRIM( output_unit ) 
          INQUIRE( FILE = filename, EXIST = filexx )
          IF ( filexx ) THEN
             OPEN( control%error, FILE = filename, FORM = 'FORMATTED',         &
                    STATUS = 'OLD', IOSTAT = iores )
          ELSE
             OPEN( control%error, FILE = filename, FORM = 'FORMATTED',         &
                     STATUS = 'NEW', IOSTAT = iores )
          END IF
        END IF

        IF ( control%out > 0 ) THEN
          INQUIRE( control%out, OPENED = opened )
          IF ( .NOT. opened ) THEN
            WRITE( output_unit, "( I0 )" ) control%out
            filename = "output_qpc." // TRIM( output_unit ) 
            INQUIRE( FILE = filename, EXIST = filexx )
            IF ( filexx ) THEN
               OPEN( control%out, FILE = filename, FORM = 'FORMATTED',         &
                      STATUS = 'OLD', IOSTAT = iores )
            ELSE
               OPEN( control%out, FILE = filename, FORM = 'FORMATTED',         &
                       STATUS = 'NEW', IOSTAT = iores )
            END IF
          END IF
        END IF

!  Create inform output structure

        plhs( i_arg ) = mxCreateStructMatrix( 1, 1, ninform, finform )

!  Define the components of the structure

        CALL MATLAB_create_integer_component( plhs( i_arg ),                   &
          'status', status_pr )
        CALL MATLAB_create_integer_component( plhs( i_arg ),                   &
           'alloc_status', alloc_status_pr )
        CALL MATLAB_create_char_component( plhs( i_arg ),                      &
          'bad_alloc', bad_alloc_pr )
        CALL MATLAB_create_logical_component( plhs( i_arg ),                   &
          'p_found', p_found_pr )
        CALL MATLAB_create_real_component( plhs( i_arg ),                      &
          'obj', obj_pr )
        CALL MATLAB_create_substructure( plhs( i_arg ),                        &
          'time', time_pr, t_ninform, t_finform )

!  Define the components of sub-structure time

        CALL MATLAB_create_real_component( time_pr,                            &
          'total', time_total_pr )
        CALL MATLAB_create_real_component( time_pr,                            &
          'preprocess', time_preprocess_pr )

!  Import the problem data

!  Check to ensure the input for H is a number

        h_in = prhs( h_arg )
        IF ( mxIsNumeric( h_in ) == 0 ) THEN
           CALL mexErrMsgTxt( ' There must be a matrix H ' )
        END IF

!  Check to ensure the input for A is a number

        a_in = prhs( a_arg )
        IF ( mxIsNumeric( a_in ) == 0 ) THEN
           CALL mexErrMsgTxt( ' There must be a matrix A ' )
        END IF

!  Get the numbers of variables and general constraints

        p%m = mxGetM( a_in )
        p%n = mxGetN( a_in )
        p%A%m = p%m
        p%A%n = p%n
        IF ( mxIsSparse( a_in ) ) THEN
          p%A%ne = mxGetNzmax( a_in )
        ELSE
          p%A%ne = p%m * p%n
        END IF
        CALL SMT_put( p%A%type, 'COORDINATE', stat )
        p%H%m = mxGetM( a_in )
        p%H%n = mxGetN( a_in )
        IF ( mxIsSparse( h_in ) ) THEN
          p%H%ne = mxGetNzmax( h_in )
        ELSE
          p%H%ne = p%n * p%n
        END IF
        CALL SMT_put( p%H%type, 'COORDINATE', stat )

!  Allocate space for the input matrix H

        ALLOCATE( p%H%row( p%H%ne ), p%H%col( p%H%ne ),                        &
                  p%H%val( p%H%ne ), STAT = info )

!  Set the row and column indices if the matrix is sparse

        IF ( mxIsSparse( h_in ) .OR. mxIsSparse( a_in ) ) THEN
          ALLOCATE( col_ptr( p%n + 1 ) )
        END IF

        IF ( mxIsSparse( h_in ) ) THEN

!  Copy the integer components of A into co-ordinate form.
!  ** N.B. indices start at 0 in C **

          h_row_pr = mxGetIr( h_in )
          h_cpr_pr = mxGetJc( h_in )

          CALL MATLAB_copy_from_ptr( h_row_pr, p%H%row, p%H%ne )
          CALL MATLAB_copy_from_ptr( h_cpr_pr, col_ptr, p%n + 1 )

          col_ptr = col_ptr + 1
          p%H%row = p%H%row + 1
          DO i = 1, p%n
            DO j = col_ptr( i ), col_ptr( i + 1 ) - 1
              p%H%col( j ) = i
            END DO
          END DO

!  Set the row and column indices if the matrix is dense

        ELSE
          l = 0
          DO j = 1, p%n
            DO i = 1, p%n
              l = l + 1
              p%H%row( l ) = i
              p%H%col( l ) = j
            END DO
          END DO
        END IF

!  Copy the real components of A

        val_pr = mxGetPr( h_in )
        CALL MATLAB_copy_from_ptr( val_pr, p%H%val, p%H%ne )

!  Now remove the lower triangle

        l = 0
        DO k = 1, p%H%ne
          i = p%H%row( k )
          j = p%H%col( k )
          IF ( i <= j ) THEN
            l = l + 1
            p%H%row( l ) = i
            p%H%col( l ) = j
            p%H%val( l ) = p%H%val( k )
          END IF
        END DO
        p%H%ne = l

!  Allocate space for the input matrix A

        ALLOCATE( p%A%row( p%A%ne ), p%A%col( p%A%ne ),                        &
                  p%A%val( p%A%ne ), STAT = info )

!  Set the row and column indices if the matrix is sparse

        IF ( mxIsSparse( a_in ) ) THEN

!  Copy the integer components of A into co-ordinate form.
!  ** N.B. indices start at 0 in C **

          a_row_pr = mxGetIr( a_in )
          a_cpr_pr = mxGetJc( a_in )

          CALL MATLAB_copy_from_ptr( a_row_pr, p%A%row, p%A%ne )
          CALL MATLAB_copy_from_ptr( a_cpr_pr, col_ptr, p%A%n + 1 )

          col_ptr = col_ptr + 1
          p%A%row = p%A%row + 1
          DO i = 1, p%n
            DO j = col_ptr( i ), col_ptr( i + 1 ) - 1
              p%A%col( j ) = i
            END DO
          END DO

!  Set the row and column indices if the matrix is dense

        ELSE
          l = 0
          DO j = 1, p%n
            DO i = 1, p%m
              l = l + 1
              p%A%row( l ) = i
              p%A%col( l ) = j
            END DO
          END DO
        END IF
        IF ( mxIsSparse( h_in ) .OR. mxIsSparse( a_in ) ) THEN
          DEALLOCATE( col_ptr )
        END IF

!  Copy the real components of A

        val_pr = mxGetPr( a_in )
        CALL MATLAB_copy_from_ptr( val_pr, p%A%val, p%A%ne )

!  Allocate space for input vectors

        ALLOCATE( p%G( p%n ), p%X_l( p%n ), p%X_u( p%n ),                      &
                  p%c( p%m ), p%C_l( p%m ), p%C_u( p%m ), STAT = info )

!  Input g

        g_in = prhs( g_arg )
        IF ( mxIsNumeric( g_in ) == 0 ) THEN
           CALL mexErrMsgTxt( ' There must be a vector g ' )
        END IF
        g_pr = mxGetPr( g_in )
        CALL MATLAB_copy_from_ptr( g_pr, p%G, p%n )

!  Input f

        f_in = prhs( f_arg )
        IF ( mxIsNumeric( f_in ) == 0 ) THEN
           CALL mexErrMsgTxt( ' There must a scalar f ' )
        END IF
        f_pr = mxGetPr( f_in )
        CALL MATLAB_copy_from_ptr( f_pr, p%f )

!  Input x_l

        xl_in = prhs( xl_arg )
        IF ( mxIsNumeric( xl_in ) == 0 ) THEN
           CALL mexErrMsgTxt( ' There must be a vector x_l ' )
        END IF
        xl_pr = mxGetPr( xl_in )
        CALL MATLAB_copy_from_ptr( xl_pr, p%X_l, p%n )
        p%X_l = MAX( p%X_l, - 10.0_wp * CONTROL%infinity )

!  Input x_u

        xu_in = prhs( xu_arg )
        IF ( mxIsNumeric( xu_in ) == 0 ) THEN
           CALL mexErrMsgTxt( ' There must be a vector x_u ' )
        END IF
        xu_pr = mxGetPr( xu_in )
        CALL MATLAB_copy_from_ptr( xu_pr, p%X_u, p%n )
        p%X_u = MIN( p%X_u, 10.0_wp * CONTROL%infinity )

!  Input c_l

        cl_in = prhs( cl_arg )
        IF ( mxIsNumeric( cl_in ) == 0 ) THEN
           CALL mexErrMsgTxt( ' There must be a vector c_l ' )
        END IF
        cl_pr = mxGetPr( cl_in )
        CALL MATLAB_copy_from_ptr( cl_pr, p%C_l, p%m )
        p%C_l = MAX( p%C_l, - 10.0_wp * CONTROL%infinity )

!  Input c_u

        cu_in = prhs( cu_arg )
        IF ( mxIsNumeric( cu_in ) == 0 ) THEN
           CALL mexErrMsgTxt( ' There must be a vector c_u ' )
        END IF
        cu_pr = mxGetPr( cu_in )
        CALL MATLAB_copy_from_ptr( cu_pr, p%C_u, p%m )
        p%C_u = MIN( p%C_u, 10.0_wp * CONTROL%infinity )

!  Allocate space for the solution

        ALLOCATE( p%X( p%n ), p%Z( p%n ), p%Y( p%m ), STAT = info )

        p%X = 0.0_wp
        p%Y = 0.0_wp
        p%Z = 0.0_wp

        ALLOCATE( B_stat( p%n ), C_stat( p%m ), STAT = info )

!  Solve the QP

        CALL VE12_solve( p, C_stat, B_stat, data, control, inform )

        inform%alloc_status = control%QPA_control%factor

!  Output solution

        plhs( x_arg ) = MATLAB_create_real( p%n, 1 )
        x_pr = mxGetPr( plhs( x_arg ) )
        CALL MATLAB_copy_to_ptr( p%X, x_pr, p%n )

!  Record output information

        CALL MATLAB_copy_to_ptr( inform%status,                                &
           mxGetPr( status_pr ) )
        CALL MATLAB_copy_to_ptr( inform%alloc_status,                          &
           mxGetPr( alloc_status_pr ) )
        CALL  MATLAB_copy_to_ptr( plhs( i_arg ),                               &
           'bad_alloc', inform%bad_alloc )
        CALL MATLAB_copy_to_ptr( inform%p_found,                               &
          mxGetPr( p_found_pr ) )
        CALL MATLAB_copy_to_ptr( inform%obj,                                   &
          mxGetPr( obj_pr ) )

!  time components

        CALL MATLAB_copy_to_ptr( REAL( inform%time%total, wp ),                &
           mxGetPr( time_total_pr ) )
        CALL MATLAB_copy_to_ptr( REAL( inform%time%preprocess, wp ),           &
           mxGetPr( time_preprocess_pr ) )

!  if required, set auxiliary output containing Lagrange multipliesr and
!  constraint bound status

        IF ( nlhs == aux_arg ) THEN

!  set up space for the auxiliary arrays

          plhs( aux_arg ) = mxCreateStructMatrix( 1, 1, naux, faux )

          CALL MATLAB_create_real_component( plhs( aux_arg ),                   &
            'c', p%m, aux_c_pr )

!  copy the values

          c_pr = mxGetPr( aux_c_pr )
          CALL MATLAB_copy_to_ptr( p%C, c_pr, p%m )

        END IF

!  Check for errors

        IF ( inform%status < 0 ) THEN
          CALL mexErrMsgTxt( ' Call to VE12_solve failed ' )
        END IF
      END IF

!  all components now set

      IF ( TRIM( mode ) == 'final' .OR. TRIM( mode ) == 'all' ) THEN
        IF ( ALLOCATED( p%H%row ) ) DEALLOCATE( p%H%row, STAT = info )
        IF ( ALLOCATED( p%H%col ) ) DEALLOCATE( p%H%col, STAT = info )
        IF ( ALLOCATED( p%H%val ) ) DEALLOCATE( p%H%val, STAT = info )
        IF ( ALLOCATED( p%G ) ) DEALLOCATE( p%G, STAT = info )
        IF ( ALLOCATED( p%A%row ) ) DEALLOCATE( p%A%row, STAT = info )
        IF ( ALLOCATED( p%A%col ) ) DEALLOCATE( p%A%col, STAT = info )
        IF ( ALLOCATED( p%A%val ) ) DEALLOCATE( p%A%val, STAT = info )
        IF ( ALLOCATED( p%C_l ) ) DEALLOCATE( p%C_l, STAT = info )
        IF ( ALLOCATED( p%C_u ) ) DEALLOCATE( p%C_u, STAT = info )
        IF ( ALLOCATED( p%X_l ) ) DEALLOCATE( p%X_l, STAT = info )
        IF ( ALLOCATED( p%X_u ) ) DEALLOCATE( p%X_u, STAT = info )
        IF ( ALLOCATED( p%X ) ) DEALLOCATE( p%X, STAT = info )
        IF ( ALLOCATED( p%Y ) ) DEALLOCATE( p%Y, STAT = info )
        IF ( ALLOCATED( p%Z ) ) DEALLOCATE( p%Z, STAT = info )
        IF ( ALLOCATED( p%C ) ) DEALLOCATE( p%C, STAT = info )
        IF ( ALLOCATED( C_stat ) ) DEALLOCATE( C_stat, STAT = info )
        IF ( ALLOCATED( B_stat ) ) DEALLOCATE( B_stat, STAT = info )
        CALL VE12_terminate( data, control, inform )
      END IF

!  close any opened io units

      IF ( control%error > 0 ) THEN
         INQUIRE( control%error, OPENED = opened )
         IF ( opened ) CLOSE( control%error )
      END IF

      IF ( control%out > 0 ) THEN
         INQUIRE( control%out, OPENED = opened )
         IF ( opened ) CLOSE( control%out )
      END IF

      RETURN
      END SUBROUTINE mexFunction


#include <fintrf.h>

! *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
!                 MEX INTERFACE TO GALAHAD_QPC
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
!      the derived type QPC_CONTROL as described in the 
!      manual for the fortran 90 package GALAHARS_QPC.
!      See: http://galahad.rl.ac.uk/galahad-www/doc/qpc.pdf
!
!  Usual Output -
!   x: a local minimizer
!
!  Optional Output -
!   inform: a structure containing information parameters
!      The components are of the form inform.value, where
!      value is the name of the corresponding component of
!      the derived type QPC_INFORM as described in the manual for 
!      the fortran 90 package GALAHARS_QPC.
!      See: http://galahad.rl.ac.uk/galahad-www/doc/qpc.pdf
!  aux: a structure containing Lagrange multipliers and constraint status
!   aux.y: Lagrange multipliers corresponding to the general constraints 
!        c_l <= A * x <= c_u 
!   aux.z: dual variables corresponding to the bound constraints
!        x_l <= x <= x_u
!   aux.c_stat: vector indicating the status of the general constraints
!           c_stat(i) < 0 if (c_l)_i = (A * x)_i
!           c_stat(i) = 0 if (c_i)_i < (A * x)_i < (c_u)_i 
!           c_stat(i) > 0 if (c_u)_i = (A * x)_i
!   aux.x_stat: vector indicating the status of the bound constraints
!           x_stat(i) < 0 if (x_l)_i = (x)_i
!           x_stat(i) = 0 if (x_i)_i < (x)_i < (x_u)_i 
!           x_stat(i) > 0 if (x_u)_i = (x)_i
!
! *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   originally released with GALAHAD Version 2.1. July 12th 2007

!  For full documentation, see 
!   http://galahad.rl.ac.uk/galahad-www/specs.html

      SUBROUTINE mexFunction( nlhs, plhs, nrhs, prhs )
      USE GALAHAD_MATLAB
      USE GALAHAD_QPC_double
      IMPLICIT NONE
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )

! ------------------------- Do not change -------------------------------

!  Keep the above subroutine, argument, and function declarations for use
!  in all your fortran mex files.
!
      INTEGER * 4 :: nlhs, nrhs
      mwPointer :: plhs( * ), prhs( * )

      LOGICAL :: mxIsSparse, mxIsChar
      mwSize :: mxGetString
      mwSize :: mxGetM, mxGetN, mxGetIr, mxGetJc
      mwSize :: mxIsNumeric, mxGetNzmax
      mwPointer :: mxCreateStructMatrix, mxGetPr

! -----------------------------------------------------------------------

!  local variables

      mwSize :: i, j, k, l, stat, info
      mwSize :: h_in, g_in, f_in, a_in, cl_in, cu_in, xl_in, xu_in
      mwSize :: g_pr, f_pr, cl_pr, cu_pr, xl_pr, xu_pr
      mwSize :: h_arg, g_arg, f_arg, a_arg, cl_arg, cu_arg
      mwSize :: xl_arg, xu_arg, c_arg, x_arg, i_arg
      mwPointer :: val_pr, x_pr
      mwPointer :: h_cpr_pr, h_row_pr, a_cpr_pr, a_row_pr

      mwPointer :: status_pr, alloc_status_pr, bad_alloc_pr
      mwPointer :: factorization_status_pr, factorization_integer_pr
      mwPointer :: factorization_real_pr, nfacts_pr, nmods_pr
      mwPointer :: p_found_pr, obj_pr, non_negligible_pivot_pr

      mwPointer :: time_pr
      mwPointer :: time_total_pr
        
!        time

!total
!preprocess
!find_dependent
!analyse
!factorize
!solve
!phase1_total
!phase1_analyse
!phase1_factorize
!phase1_solve


!        QPA_inform
!        QPB_inform

!     CHARACTER ( len = 80 ) :: message
      CHARACTER ( len = 80 ) :: output_unit, filename
      LOGICAL :: filexx, opened, initial_set = .FALSE.
      INTEGER :: iores

      CHARACTER ( len = 7 ) :: mode
      INTEGER, ALLOCATABLE, DIMENSION( : ) :: C_stat, B_stat

      mwSize, PARAMETER :: ninform = 14
      CHARACTER ( LEN = 21 ), PARAMETER :: finform( ninform ) = (/ 
     *     'status               ', 'alloc_status         ',
     *     'bad_alloc            ', 'factorization_status ', 
     *     'factorization_integer', 'factorization_real   ',
     *     'nfacts               ', 'nmods                ', 
     *     'p_found              ', 'obj                  ', 
     *     'non_negligible_pivot ',
     *     'time                 ',
     *     'QPA_inform           ', 'QPB_inform           '
     *   /)
      mwSize, PARAMETER :: t_ninform = 10
      CHARACTER ( LEN = 21 ), PARAMETER :: t_finform( t_ninform ) = (/ 
     *     'total                ', 'preprocess           ', 
     *     'find_dependent       ', 'analyse              ', 
     *     'factorize            ', 'solve                ', 
     *     'phase1_total         ', 'phase1_analyse       ', 
     *     'phase1_factorize     ', 'phase1_solve         '
     *   /)

!  arguments for QPC

      TYPE ( QPT_problem_type ), SAVE :: p
      TYPE ( QPC_control_type ), SAVE :: control
      TYPE ( QPC_inform_type ) :: inform
      TYPE ( QPC_data_type ), SAVE :: data

      mwPointer, ALLOCATABLE :: col_ptr( : )

!  Test input/output arguments

      IF ( nrhs < 1 ) THEN
        CALL mexErrMsgTxt( ' qpc requires at least 1 input argument' )
      END IF

      IF ( mxIsChar( prhs( 1 ) ) ) THEN
        i = mxGetString( prhs( 1 ), mode, 7 )
        IF ( .NOT. ( TRIM( mode ) == 'initial' .OR. 
     &               TRIM( mode ) == 'final' ) ) THEN
          IF ( nrhs < 2 ) THEN
            CALL mexErrMsgTxt( ' Too few input arguments to qpc' )
          END IF
          h_arg = 2
          g_arg = 3
          f_arg = 4
          a_arg = 5
          cl_arg = 6
          cu_arg = 7
          xl_arg = 8
          xu_arg = 9
          c_arg = 10
          x_arg = 1
          i_arg = 2
          IF ( nrhs > c_arg ) THEN
            CALL mexErrMsgTxt( ' Too many input arguments to qpc' )
          END IF
        END IF
      ELSE
!       CALL mexErrMsgTxt( ' First input to qpc must be a string' )
        mode = 'all'
        IF ( nrhs < 2 ) THEN
          CALL mexErrMsgTxt( ' Too few input arguments to qpc' )
        END IF
        h_arg = 1
        g_arg = 2
        f_arg = 3
        a_arg = 4
        cl_arg = 5
        cu_arg = 6
        xl_arg = 7
        xu_arg = 8
        c_arg = 9
        x_arg = 1
        i_arg = 2
        IF ( nrhs > c_arg ) THEN
          CALL mexErrMsgTxt( ' Too many input arguments to qpc' )
        END IF
      END IF

!     IF ( nlhs > 2 ) THEN
      IF ( nlhs > 3 ) THEN
        CALL mexErrMsgTxt( 
     &    ' qpc provides at most 2 output arguments' )
      END IF

!  Initialize the internal structures for qpc

      IF ( TRIM( mode ) == 'initial' .OR. 
     &     TRIM( mode ) == 'all' ) THEN
        initial_set = .TRUE.
        CALL QPC_initialize( data, control )
        IF ( TRIM( mode ) == 'initial' ) RETURN
      END IF

      IF ( .NOT. TRIM( mode ) == 'final' ) THEN

!  Check that QPC_initialize has been called 

        IF ( .NOT. initial_set ) THEN
          CALL mexErrMsgTxt( ' "initial" must be called first' )
        END IF


!  If the third argument is present, extract the input control data


!  1*******


!  No i/o allowed from mex interface?

!        control%lp = 0
!        control%wp = 0
!        control%mp = 0
!        control%sp = 0

!  Open i/o units

        IF ( control%error > 0 ) THEN
          WRITE( output_unit, "( I0 )" ) control%error
          filename = "output_qpc." // TRIM( output_unit ) 
          INQUIRE( FILE = filename, EXIST = filexx )
          IF ( filexx ) THEN
             OPEN( control%error, FILE = filename, FORM = 'FORMATTED',
     &              STATUS = 'OLD', IOSTAT = iores )
          ELSE
             OPEN( control%error, FILE = filename, FORM = 'FORMATTED',
     &               STATUS = 'NEW', IOSTAT = iores )
          END IF
        END IF

        IF ( control%out > 0 ) THEN
          INQUIRE( control%out, OPENED = opened )
          IF ( .NOT. opened ) THEN
            WRITE( output_unit, "( I0 )" ) control%out
            filename = "output_qpc." // TRIM( output_unit ) 
            INQUIRE( FILE = filename, EXIST = filexx )
            IF ( filexx ) THEN
               OPEN( control%out, FILE = filename, FORM = 'FORMATTED',
     &                STATUS = 'OLD', IOSTAT = iores )
            ELSE
               OPEN( control%out, FILE = filename, FORM = 'FORMATTED',
     &                 STATUS = 'NEW', IOSTAT = iores )
            END IF
          END IF
        END IF


!  Create inform output structure

        plhs( i_arg ) 
     &    = mxCreateStructMatrix( 1, 1, ninform, finform )

!  Define the components of the structure

        CALL MATLAB_create_integer_component( plhs( i_arg ), 
     &    'status', status_pr )
        CALL MATLAB_create_integer_component( plhs( i_arg ),
     &     'alloc_status', alloc_status_pr )
        CALL MATLAB_create_char_component( plhs( i_arg ), 
     &    'bad_alloc', bad_alloc_pr )
        CALL MATLAB_create_integer_component( plhs( i_arg ),
     &    'factorization_status', factorization_status_pr )
        CALL MATLAB_create_integer_component( plhs( i_arg ),
     &    'factorization_integer', factorization_integer_pr )
        CALL MATLAB_create_integer_component( plhs( i_arg ),
     &    'factorization_real', factorization_real_pr )
        CALL MATLAB_create_integer_component( plhs( i_arg ), 
     &    'nfacts', nfacts_pr )
        CALL MATLAB_create_integer_component( plhs( i_arg ), 
     &    'nmods', nmods_pr )
        CALL MATLAB_create_logical_component( plhs( i_arg ), 
     &    'p_found', p_found_pr )
        CALL MATLAB_create_real_component( plhs( i_arg ), 
     &    'obj', obj_pr )
        CALL MATLAB_create_real_component( plhs( i_arg ),
     &     'non_negligible_pivot', non_negligible_pivot_pr )
        CALL MATLAB_create_substructure( plhs( i_arg ), 
     &    'time', time_pr, t_ninform, t_finform )

!  Create time output sub-structure of inform

        time_pr
     &    = mxCreateStructMatrix( 1, 1, t_ninform, t_finform )

!  Define the components of the sub-structure

        CALL MATLAB_create_real_component( time_pr,
     &    'total', time_total_pr )

!  Create qpa_inform output sub-structure of inform

!        CALL MATLAB_create_integer_component( plhs( i_arg ), 
!     &    'QPA_inform', QPA_inform_pr )

!  Create qpb_inform output sub-structure of inform

!        CALL MATLAB_create_integer_component( plhs( i_arg ), 
!     &    'QPB_inform', QPB_inform_pr )

!  Initialise the components (not needed here, for example only)

!       CALL MATLAB_copy_to_prt( -1, mxGetPr( flag_pr ), 1 )     
!       CALL MATLAB_copy_to_prt( -1.0_wp, mxGetPr( opsa_pr ), 1 )     

!  Import the problem data

         p%new_problem_structure = .TRUE.

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

        ALLOCATE( p%H%row( p%H%ne ), p%H%col( p%H%ne ), 
     &            p%H%val( p%H%ne ), STAT = info )

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
       write(6,*) 'row', p%H%row
       write(6,*) 'ptr', col_ptr
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

        ALLOCATE( p%A%row( p%A%ne ), p%A%col( p%A%ne ), 
     &            p%A%val( p%A%ne ), STAT = info )

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

        ALLOCATE( p%G( p%n ), p%X_l( p%n ), p%X_u( p%n ),
     &            p%c( p%m ), p%C_l( p%m ), p%C_u( p%m ), 
     &            STAT = info )

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
        CALL MATLAB_copy_from_ptr( f_pr, p%f, 1 )

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

!         control%print_level = 100
!         control%QPA_control%print_level = 100
!         control%QPB_control%out = 6
!         control%QPB_control%print_level = 100
!         control%QPB_control%LSQP_control%out = 6
!         control%QPB_control%LSQP_control%print_level = 100

!         write(6,*) 'H%val', p%H%val( : p%H%ne )
!         write(6,*) 'H%row', p%H%row( : p%H%ne )
!         write(6,*) 'H%col', p%H%col( : p%H%ne )

!         write(6,*) 'g', p%g
!         write(6,*) 'f', p%f

!         write(6,*) 'A%val', p%A%val( : p%A%ne )
!         write(6,*) 'A%row', p%A%row( : p%A%ne )
!         write(6,*) 'A%col', p%A%col( : p%A%ne )

!         write(6,*) 'c_l', p%C_l
!         write(6,*) 'c_u', p%C_u

!         write(6,*) 'x_l', p%X_l
!         write(6,*) 'x_u', p%X_u

        CALL QPC_solve( p, C_stat, B_stat, data, control, inform )

!  Output solution

        plhs( x_arg ) = galmxCreateRealMatrix( p%n, 1 )
        x_pr = mxGetPr( plhs( x_arg ) )
        CALL MATLAB_copy_to_prt( p%X, x_pr, p%n )

!  Record output information

        CALL MATLAB_copy_to_prt( inform%status, 
     &     mxGetPr( status_pr ) )
        CALL MATLAB_copy_to_prt( inform%alloc_status, 
     &     mxGetPr( alloc_status_pr ) )
        CALL galmxSetCharacterComponent( plhs( i_arg ), 
     &    'bad_alloc' )
        CALL MATLAB_copy_to_prt( inform%factorization_status,
     &    mxGetPr( factorization_status_pr ) )
        CALL MATLAB_copy_to_prt( inform%factorization_integer,
     &    mxGetPr( factorization_integer_pr ) )
        CALL MATLAB_copy_to_prt( inform%factorization_real,
     &    mxGetPr( factorization_real_pr ) )
        CALL MATLAB_copy_to_prt( inform%nfacts,
     &    mxGetPr( nfacts_pr ) )
        CALL MATLAB_copy_to_prt( inform%nmods,
     &    mxGetPr( nmods_pr ) )
        CALL galmxCopyLogicalToPtr( inform%p_found,
     &    mxGetPr( p_found_pr ) )
        CALL MATLAB_copy_to_prt( inform%obj,
     &    mxGetPr( obj_pr ) )
        CALL MATLAB_copy_to_prt( inform%non_negligible_pivot,
     &    mxGetPr( non_negligible_pivot_pr ) )


        CALL MATLAB_copy_to_prt( inform%time%total,
     &    mxGetPr( time_total_pr ) )




      i = 0
      if ( i /= 1 ) RETURN


!        CALL MATLAB_copy_to_prt( inform%time,
!     &    mxGetPr( time_pr ) )
!        CALL MATLAB_copy_to_prt( inform%QPA_inform,
!     &    mxGetPr( QPA_inform_pr ) )
!        CALL MATLAB_copy_to_prt( inform%QPB_inform,
!     &    mxGetPr( QPB_inform_pr ) )


!  Check for errors

        IF ( inform%status < 0 ) THEN
          CALL mexErrMsgTxt( ' Call to QPC_solve failed ' )
        END IF
      END IF

!  set inform components

!      plhs( 1 ) = galmxCreateIntegerMatrix( 1, 1 )
!      status_pr = mxGetPr( plhs( 1 ) )
!      CALL MATLAB_copy_to_prt( inform%status, status_pr, 1 )
!      plhs( 2 ) = galmxCreateRealMatrix( 1, 1 )
!      inf_pr = mxGetPr( plhs( 2 ) )
!      CALL MATLAB_copy_to_prt( p%C_l, inf_pr, 1 )
!      CLOSE( 6 )
!      i = 0
!      if ( i /= 1 ) RETURN

!      WRITE( message, * ) inform%status
!      CALL MEXPRINTF( TRIM( message ) // char( 13 ) )

!  all components now set

      IF ( TRIM( mode ) == 'final' .OR. 
     &     TRIM( mode ) == 'all' ) THEN
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
        CALL QPC_terminate( data, control, inform )
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


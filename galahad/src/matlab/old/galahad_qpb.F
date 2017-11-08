#include <fintrf.h>

! *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
!                 MEX INTERFACE TO GALAHAD_QPB
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
!     = galahad_qpb( H, g, f, A, c_l, c_u, x_l, x_u, control )
!
!  Sophisticated usage -
!
!  to initialize data structures prior to solution
!   galahad_qpb( 'initial' )
!
!  to solve the quadratic program using existing data structures
!   [ x, inform, aux ]
!     = galahad_qpb( H, g, f, A, c_l, c_u, x_l, x_u, control )
!
!  to remove data structures after solution
!   galahad_qpb( 'final' )
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
!      the derived type QPB_CONTROL as described in the 
!      manual for the fortran 90 package GALAHARS_QPB.
!      See: http://galahad.rl.ac.uk/galahad-www/doc/qpb.pdf
!
!  Usual Output -
!   x: a local minimizer
!
!  Optional Output -
!   inform: a structure containing information parameters
!      The components are of the form inform.value, where
!      value is the name of the corresponding component of
!      the derived type QPB_INFORM as described in the manual for 
!      the fortran 90 package GALAHARS_QPB.
!      See: http://galahad.rl.ac.uk/galahad-www/doc/qpb.pdf
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
      USE GALAHAD_QPB_double
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
      REAL ( KIND = wp ) :: mxGetScalar,  mxGetInf, infinity
      mwSize :: mxGetString
      mwSize :: mxGetM, mxGetN, mxGetIr, mxGetJc
      mwSize :: mxIsNumeric, mxGetNzmax, mxGetNumberOfFields
      mwPointer :: mxCreateStructMatrix, mxGetPr, mxGetField
      mwPointer :: galmxCreateInteger, galmxCreateReal
      mwPointer :: galmxCreateIntegerMatrix, galmxCreateRealMatrix

! -----------------------------------------------------------------------

!  local variables

      mwSize :: i, j, k, l, dig, stat, n, nb, nfields, info
      mwSize :: h_in, g_in, f_in, a_in, cl_in, cu_in, xl_in, xu_in
      mwSize :: g_pr, f_pr, cl_pr, cu_pr, xl_pr, xu_pr, stat_pr
      mwSize :: h_arg, g_arg, f_arg, a_arg, cl_arg, cu_arg
      mwSize :: xl_arg, xu_arg, c_arg, x_arg, au_arg, i_arg
      mwPointer :: val_pr, rhs_in, cn_in, x_pr
      mwPointer :: h_cpr_pr, h_row_pr, a_cpr_pr, a_row_pr
      mwPointer :: flag_pr, more_pr, status_pr, inf_pr
      mwPointer :: nsteps_pr, nrltot_pr, nirtot_pr, nrlnec_pr
      mwPointer :: nirnec_pr, nrladu_pr, niradu_pr, ncmpa_pr
      mwPointer :: nfaulty_pr, oor_pr, dup_pr, maxfrt_pr
      mwPointer :: nebdu_pr, nrlbdu_pr, nirbdu_pr, ncmpbr_pr
      mwPointer :: ncmpbi_pr, ntwo_pr, neig_pr, delay_pr
      mwPointer :: signc_pr, modstep_pr, rank_pr, step_pr
      mwPointer :: opsa_pr, opse_pr, opsb_pr, maxchange_pr
      mwPointer :: string_pr

      CHARACTER ( len = 80 ) :: message, output_unit, filename
      LOGICAL :: filexx, opened, initial_set = .FALSE.
      INTEGER :: iores

      CHARACTER ( len = 7 ) :: mode

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

!  arguments for QPB

      TYPE ( QPT_problem_type ), SAVE :: p
      TYPE ( QPB_control_type ), SAVE :: control
      TYPE ( QPB_inform_type ) :: inform
      TYPE ( QPB_data_type ), SAVE :: data

      mwPointer, ALLOCATABLE :: col_ptr( : )

!  Test input/output arguments

      IF ( nrhs < 1 ) THEN
        CALL mexErrMsgTxt( ' qpb requires at least 1 input argument' )
      END IF

      IF ( mxIsChar( prhs( 1 ) ) ) THEN
        i = mxGetString( prhs( 1 ), mode, 7 )
        IF ( .NOT. ( TRIM( mode ) == 'initial' .OR. 
     &               TRIM( mode ) == 'final' ) ) THEN
          IF ( nrhs < 2 ) THEN
            CALL mexErrMsgTxt( ' Too few input arguments to qpb' )
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
          au_arg = 3
          IF ( nrhs > c_arg ) THEN
            CALL mexErrMsgTxt( ' Too many input arguments to qpb' )
          END IF
        END IF
      ELSE
!       CALL mexErrMsgTxt( ' First input to qpb must be a string' )
        mode = 'all'
        IF ( nrhs < 2 ) THEN
          CALL mexErrMsgTxt( ' Too few input arguments to qpb' )
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
        au_arg = 3
        IF ( nrhs > c_arg ) THEN
          CALL mexErrMsgTxt( ' Too many input arguments to qpb' )
        END IF
      END IF

      IF ( nlhs > 2 ) THEN
        CALL mexErrMsgTxt( 
     &    ' qpb provides at most 2 output arguments' )
      END IF

!  Initialize the internal structures for qpb

      IF ( TRIM( mode ) == 'initial' .OR. 
     &     TRIM( mode ) == 'all' ) THEN
        initial_set = .TRUE.
        CALL QPB_initialize( data, control )
        IF ( TRIM( mode ) == 'initial' ) RETURN
      END IF

      IF ( .NOT. TRIM( mode ) == 'final' ) THEN

!  Check that QPB_initialize has been called 

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
          filename = "output_qpb." // TRIM( output_unit ) 
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
            filename = "output_qpb." // TRIM( output_unit ) 
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

! 3 ***************

!  Initialise the components (not needed here, for example only)

!       CALL galmxCopyIntegerToPtr( -1, mxGetPr( flag_pr ), 1 )     
!       CALL galmxCopyRealToPtr( -1.0_wp, mxGetPr( opsa_pr ), 1 )     

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
     &              p%H%val( p%H%ne ), STAT = info )

!  Set the row and column indices if the matrix is sparse

          IF ( mxIsSparse( h_in ) .OR. mxIsSparse( a_in ) ) THEN
            ALLOCATE( col_ptr( p%n + 1 ) )
          END IF

          IF ( mxIsSparse( h_in ) ) THEN

!  Copy the integer components of A into co-ordinate form.
!  ** N.B. indices start at 0 in C **

            h_row_pr = mxGetIr( h_in )
            h_cpr_pr = mxGetJc( h_in )

            CALL galmxCopyPtrToInteger( h_row_pr, p%H%row, p%H%ne )
            CALL galmxCopyPtrToInteger( h_cpr_pr, col_ptr, p%n + 1 )

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
          CALL galmxCopyPtrToReal( val_pr, p%H%val, p%H%ne )

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
     &              p%A%val( p%A%ne ), STAT = info )

!  Set the row and column indices if the matrix is sparse

          IF ( mxIsSparse( a_in ) ) THEN

!  Copy the integer components of A into co-ordinate form.
!  ** N.B. indices start at 0 in C **

            a_row_pr = mxGetIr( a_in )
            a_cpr_pr = mxGetJc( a_in )

            CALL galmxCopyPtrToInteger( a_row_pr, p%A%row, p%A%ne )
            CALL galmxCopyPtrToInteger( a_cpr_pr, col_ptr, p%A%n + 1 )

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
          CALL galmxCopyPtrToReal( val_pr, p%A%val, p%A%ne )

!  Allocate space for input vectors

          ALLOCATE( p%G( p%n ), p%X_l( p%n ), p%X_u( p%n ),
     &              p%c( p%m ), p%C_l( p%m ), p%C_u( p%m ), 
     &              STAT = info )

!  Input g

          g_in = prhs( g_arg )
          IF ( mxIsNumeric( g_in ) == 0 ) THEN
             CALL mexErrMsgTxt( ' There must be a vector g ' )
          END IF
          g_pr = mxGetPr( g_in )
          CALL galmxCopyPtrToReal( g_pr, p%G, p%n )

!  Input f

          f_in = prhs( f_arg )
          IF ( mxIsNumeric( f_in ) == 0 ) THEN
             CALL mexErrMsgTxt( ' There must a scalar f ' )
          END IF
          f_pr = mxGetPr( f_in )
          CALL galmxCopyPtrToReal( f_pr, p%f, 1 )

!  Input x_l

          xl_in = prhs( xl_arg )
          IF ( mxIsNumeric( xl_in ) == 0 ) THEN
             CALL mexErrMsgTxt( ' There must be a vector x_l ' )
          END IF
          xl_pr = mxGetPr( xl_in )
          CALL galmxCopyPtrToReal( xl_pr, p%X_l, p%n )
          p%X_l = MAX( p%X_l, - 10.0_wp * CONTROL%infinity )


!  Input x_u

          xu_in = prhs( xu_arg )
          IF ( mxIsNumeric( xu_in ) == 0 ) THEN
             CALL mexErrMsgTxt( ' There must be a vector x_u ' )
          END IF
          xu_pr = mxGetPr( xu_in )
          CALL galmxCopyPtrToReal( xu_pr, p%X_u, p%n )
          p%X_u = MIN( p%X_u, 10.0_wp * CONTROL%infinity )

!  Input c_l

          cl_in = prhs( cl_arg )
          IF ( mxIsNumeric( cl_in ) == 0 ) THEN
             CALL mexErrMsgTxt( ' There must be a vector c_l ' )
          END IF
          cl_pr = mxGetPr( cl_in )
          CALL galmxCopyPtrToReal( cl_pr, p%C_l, p%m )
          p%C_l = MAX( p%C_l, - 10.0_wp * CONTROL%infinity )

!  Input c_u

          cu_in = prhs( cu_arg )
          IF ( mxIsNumeric( cu_in ) == 0 ) THEN
             CALL mexErrMsgTxt( ' There must be a vector c_u ' )
          END IF
          cu_pr = mxGetPr( cu_in )
          CALL galmxCopyPtrToReal( cu_pr, p%C_u, p%m )
          p%C_u = MIN( p%C_u, 10.0_wp * CONTROL%infinity )

!  Allocate space for the solution

          ALLOCATE( p%X( p%n ), p%Z( p%n ), p%Y( p%m ), STAT = info )

          p%X = 0.0_wp
          p%Y = 0.0_wp
          p%Z = 0.0_wp

!  Solve the QP

           control%print_level = 100
           control%LSQP_control%out = 6
           control%LSQP_control%print_level = 100

           write(6,*) 'H%val', p%H%val( : p%H%ne )
           write(6,*) 'H%row', p%H%row( : p%H%ne )
           write(6,*) 'H%col', p%H%col( : p%H%ne )

           write(6,*) 'g', p%g
           write(6,*) 'f', p%f

           write(6,*) 'A%val', p%A%val( : p%A%ne )
           write(6,*) 'A%row', p%A%row( : p%A%ne )
           write(6,*) 'A%col', p%A%col( : p%A%ne )

           write(6,*) 'c_l', p%C_l
           write(6,*) 'c_u', p%C_u

           write(6,*) 'x_l', p%X_l
           write(6,*) 'x_u', p%X_u

!      CLOSE( 6 )
!      i = 0
!      if ( i /= 1 ) RETURN

           CALL QPB_solve( p, data, control, inform )

! 4 ******************

!  Check for errors

          IF ( inform%status < 0 ) THEN
            CALL mexErrMsgTxt( ' Call to QPB_solve failed ' )
          END IF

!  set inform components


!      plhs( 1 ) = galmxCreateIntegerMatrix( 1, 1 )
!      status_pr = mxGetPr( plhs( 1 ) )
!      CALL galmxCopyIntegerToPtr( inform%status, status_pr, 1 )
!      plhs( 2 ) = galmxCreateRealMatrix( 1, 1 )
!      inf_pr = mxGetPr( plhs( 2 ) )
!      CALL galmxCopyRealToPtr( p%C_l, inf_pr, 1 )

      WRITE( message, * ) inform%status
      CALL MEXPRINTF( TRIM( message ) // char( 13 ) )

!  Output solution

          plhs( x_arg ) = galmxCreateRealMatrix( p%n, 1 )
          x_pr = mxGetPr( plhs( x_arg ) )
          CALL galmxCopyRealToPtr( p%X, x_pr, p%n )

      END IF

!     WRITE( message, * ) 'hello'
!     CALL MEXPRINTF( TRIM( message ) // char( 13 ) )

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
        CALL QPB_terminate( data, control, inform )
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


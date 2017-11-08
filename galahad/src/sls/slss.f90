  PROGRAM SLS_EXAMPLE   !  GALAHAD 2.4 - 06/08/2009 AT 11:30 GMT.
   USE GALAHAD_SLS_double
   IMPLICIT NONE
   INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )  
   TYPE ( SMT_type ) :: matrix
   TYPE ( SLS_data_type ) :: data
   TYPE ( SLS_control_type ) control
   TYPE ( SLS_inform_type ) :: inform
   INTEGER :: i, n, ne, s
   REAL ( KIND = wp ), ALLOCATABLE :: B( : ), X( : )
! Read matrix order and number of entries
   READ( 5, * ) n, ne
   matrix%n = n
   matrix%ne = ne
! Allocate arrays of appropriate sizes
   ALLOCATE( matrix%val( ne ), matrix%row( ne ), matrix%col( ne ) )
   ALLOCATE( B( n ), X( n ) )
! Read matrix and right-hand side
   READ( 5, * ) ( matrix%row( i ), matrix%col( i ), matrix%val( i ),  i = 1, ne )
!  B = 1.0_wp
   READ( 5, * ) B
   CALL SMT_put( matrix%type, 'COORDINATE', s )     ! Specify co-ordinate 
! Specify the solver (in this case HSL_MA57)
!  CALL SLS_initialize( 'ma57', data, control, inform )
!  CALL SLS_initialize( 'ma97', data, control, inform )
   CALL SLS_initialize( 'ssids', data, control, inform )
   control%ordering = 1
!  control%scaling = -1
!  control%print_level = 2
!  control%print_level_solver = 2
! Analyse
   CALL SLS_analyse( matrix, data, control, inform )
   IF ( inform%status < 0 ) THEN
     WRITE( 6, '( A, I0 )' )                                                   &
          ' Failure of SLS_analyse with status = ', inform%status
     STOP
   END IF
! Factorize
   CALL SLS_factorize( matrix, data, control, inform )
   IF ( inform%status < 0 ) THEN
     WRITE( 6, '( A, I0 )' )                                                   &
          ' Failure of SLS_factorize with status = ', inform%status
     STOP
   END IF
! Solve using iterative refinement and ask for high relative accuracy
   X = B
   control%max_iterative_refinements = 1
   control%acceptable_residual_relative = 0.0_wp
   CALL SLS_solve( matrix, X, data, control, inform )
   IF ( inform%status == 0 ) WRITE( 6, '( A, /, ( 3F20.16 ) )' )               &
     ' Solution is', X
! Clean up
   CALL SLS_terminate( data, control, inform )
   DEALLOCATE(  matrix%type, matrix%val, matrix%row, matrix%col, X, B )
   STOP
   END PROGRAM SLS_EXAMPLE

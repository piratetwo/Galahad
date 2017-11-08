#include <fintrf.h>
      SUBROUTINE mexFunction( nlhs, plhs, nrhs, prhs )
      USE GALAHAD_SMT_double
!     USE HSL_ZD11_double

!-------------------------Do not change --------------------------------

! Keep the above subroutine, argument, and function declarations for use
! in all your fortran mex files.
!     (pointer) Replace integer by integer*8 on 64-bit platforms
!
      MWPOINTER :: plhs(*), prhs(*)
      MWPOINTER :: mxCreateDoubleMatrix, mxGetPr,  mxGetJc, mxGetIr
      INTEGER :: nlhs, nrhs
      MWSIZE :: mxGetM, mxGetN

!-----------------------------------------------------------------------

      INTEGER :: mxIsNumeric
      MWSIZE :: mxGetElementSize
      INTEGER :: a_in, dig, stat
      MWPOINTER :: val, row, col, a_pr, b_pr
      MWPOINTER :: cpr_in, row_in
      MWSIZE :: m, n, size
      REAL * 8 :: tot, at
!     REAL * 8, ALLOCATABLE :: x( : )
!     REAL * 8 ::  x( 100 )
      TYPE ( SMT_type ) :: A
!     TYPE ( ZD11_type ) :: A
      MWPOINTER, ALLOCATABLE :: col_ptr( : )

      IF ( nrhs /= 1 ) THEN
        CALL mexerrmsgtxt( 'memtest requires one input argument' )
      ELSE IF ( nlhs > 1 ) THEN
        CALL mexerrmsgtxt( 'memtest requires one output argument' )
      ENDIF

!  Create space for and fill return arguments to MATLAB.

!  Check to ensure the input is a number

      a_in = prhs( 1 )
      IF ( mxIsNumeric( a_in ) == 0 ) THEN
         CALL mexErrMsgTxt( 'Input must be a number' )
      END IF

!  Get the row and column dimensions and the size of the input matrix

      A%m = mxGetM( a_in )
      A%n = mxGetN( prhs(1 ) )
!     A%n = mxGetN( a_in )
      A%ne = mxGetNzmax( a_in )

!  Allocate space for the input matrix

      ALLOCATE( A%row( A%ne ), A%col( A%ne ), A%val( A%ne ) )
      ALLOCATE( col_ptr( A%n + 1 ) )
      CALL SMT_put( A%type, 'COORDINATE', stat )
!     CALL ZD11_put( A%type, 'COORDINATE', stat )

!  Copy the integer components of A into co-ordinate form
!  N.B. indices start at 0 in C

      cpr_in = mxGetJc( a_in )
      row_in = mxGetIr( a_in )

      dig = digits( dig ) + 1
      SELECT CASE ( dig )
      CASE ( 64 )
        CALL mxCopyPtrToInteger8( cpr_in, col_ptr, A%n + 1 )
        CALL mxCopyPtrToInteger8( row_in, A%row, A%ne )
      CASE ( 16 )
        CALL mxCopyPtrToInteger2( cpr_in, col_ptr, A%n + 1 )
        CALL mxCopyPtrToInteger2( row_in, A%row, A%ne )
      CASE ( 8 )
        CALL mxCopyPtrToInteger1( cpr_in, col_ptr, A%n + 1 )
        CALL mxCopyPtrToInteger1( row_in, A%row, A%ne )
      CASE default
        CALL mxCopyPtrToInteger4( cpr_in, col_ptr, A%n + 1 )
        CALL mxCopyPtrToInteger4( row_in, A%row, A%ne )
      END SELECT

      col_ptr = col_ptr + 1
      A%row = A%row + 1

!  Copy the real components of A

      val = mxGetPr( a_in )
      CALL mxCopyPtrToReal8( val, A%val, A%ne )

      DO i = 1, A%n
        DO j = col_ptr( i ), col_ptr( i + 1 ) - 1
          A%col( j ) = i
        END DO
      END DO
      DEALLOCATE( col_ptr )

!  Call the computational subroutine

      CALL mv( A, tot )
      DEALLOCATE( A%row, A%col, A%val )

!  Create space for and fill return arguments to MATLAB.

      plhs( 1 ) = mxCreateDoubleMatrix( 1, 1, 0 )
      a_pr = mxGetPr( plhs( 1 ) )
      CALL mxCopyReal8ToPtr( tot, a_pr, 1 )     

      plhs( 2 ) = mxCreateDoubleMatrix( 1, 1, 0 )
      b_pr = mxGetPr( plhs( 2 ) )
      CALL mxCopyReal8ToPtr( at, b_pr, 1 )     

!     plhs( 1 ) = mxCreateDoubleMatrix( m, n, 0 )
!     plhs( 1 ) = mxCreateDoubleMatrix( 1, 1, 0 )
!     tot = A%val(2)
!     a_pr = mxGetPr( plhs( 1 ) )
!     CALL mxCopyReal8ToPtr( tot, a_pr, 1 )     
!     i = 0
!     if ( i /= 1 ) RETURN

      RETURN
      END SUBROUTINE mexFunction


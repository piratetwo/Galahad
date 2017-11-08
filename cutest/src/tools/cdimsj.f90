! THIS VERSION: CUTEST 1.0 - 29/12/2012 AT 11:25 GMT.

!-*-*-*-*-*-*-  C U T E S T    C D I M S J    S U B R O U T I N E  -*-*-*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   fortran 2003 version released in CUTEst, 29th December 2012

      SUBROUTINE CUTEST_cdimsj( status, nnzj )
      USE CUTEST

!  dummy arguments

      INTEGER, INTENT( OUT ) :: nnzj, status

!  --------------------------------------------------------------
!  compute the space required to store the Jacobian matrix of the 
!  constraints/objective function of a problem initially written 
!  in Standard Input Format (SIF)
!  --------------------------------------------------------------

      CALL CUTEST_cdimsj_threadsafe( CUTEST_data_global, status, nnzj )
      RETURN

!  end of sunroutine CUTEST_cdimsj

      END SUBROUTINE CUTEST_cdimsj

!-*-  C U T E S T   C D I M S J _ t h r e a d s a f e   S U B R O U T I N E  -*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   fortran 77 version originally released in CUTEr, April 1999
!   fortran 2003 version released in CUTEst, 28th November 2012

      SUBROUTINE CUTEST_cdimsj_threadsafe( data, status, nnzj )
      USE CUTEST

!  dummy arguments

      TYPE ( CUTEST_data_type ), INTENT( IN ) :: data
      INTEGER, INTENT( OUT ) :: nnzj, status

!  --------------------------------------------------------------
!  compute the space required to store the Jacobian matrix of the 
!  constraints/objective function of a problem initially written 
!  in Standard Input Format (SIF)
!  --------------------------------------------------------------

!  local variables

      INTEGER :: ig

!  the total space is stored in nnzj

      nnzj = 0

!  allow space for constraint groups

      DO ig = 1, data%ng
        IF ( data%KNDOFC( ig ) /= 0 )                                         &
          nnzj = nnzj + data%ISTAGV( ig + 1 ) - data%ISTAGV( ig )
      END DO

!  add space for the (dense) gradient of the objective function

      nnzj = nnzj + data%n
      status = 0
      RETURN

!  end of sunroutine CUTEST_cdimsj_threadsafe

      END SUBROUTINE CUTEST_cdimsj_threadsafe

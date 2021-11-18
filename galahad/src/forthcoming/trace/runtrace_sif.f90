! THIS VERSION: GALAHAD 2.6 - 23/10/2014 AT 10:00 GMT.

!-*-*-*-*-*-*-  G A L A H A D   R U N T R A C E _ S I F  *-*-*-*-*-*-*-*-

!  Nick Gould, Dominique Orban and Philippe Toint, for GALAHAD productions
!  Copyright reserved
!  October 23rd 2014

   PROGRAM RUNTRACE_SIF
   USE GALAHAD_USETRACE_double

!  Main program for the SIF interface to TRACE, a trust-region algorithm for
!  unconstrained optimization

!  Problem insif characteristics

   INTEGER, PARAMETER :: errout = 6
   INTEGER, PARAMETER :: insif = 55
   CHARACTER ( LEN = 16 ) :: prbdat = 'OUTSDIF.d'
   INTEGER :: iostat

!  Open the data input file

   OPEN( insif, FILE = prbdat, FORM = 'FORMATTED', STATUS = 'OLD',             &
         IOSTAT = iostat )
   IF ( iostat > 0 ) THEN
     WRITE( errout,                                                            &
       "( ' ERROR: could not open file OUTSDIF.d on unit ', I2 )" ) insif
     STOP
   END IF
   REWIND insif

!  Call the SIF interface

   CALL USE_TRACE( insif )

!  Close the data input file 

   CLOSE( insif )
   STOP

!  End of RUNTRACE_SIF

   END PROGRAM RUNTRACE_SIF

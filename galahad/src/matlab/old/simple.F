#include <fintrf.h>

! *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!                 MEX INTERFACE TO GALAHAD_SILS
!
!  Expected Matlab syntax:
!
!    [X,INFORM] = SILS(A,B) or
!    [X,INFORM] = SILS(A,B,CONTROL)
!
!  possible "mode" strings
!   'initial'
!   'factor'
!   'solve'
!   'final'
!   'all'

!
! *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   originally released with GALAHAD Version 2.1. July 4th 2007

!  For full documentation, see 
!   http://galahad.rl.ac.uk/galahad-www/specs.html

      SUBROUTINE mexFunction( nlhs, plhs, nrhs, prhs )
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

      CHARACTER ( len = 80 ) :: message, output_unit, filename
      LOGICAL :: filexx, opened
      INTEGER :: iores

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

!  one right-hand side

       nb = 1
       ALLOCATE( B( 5 ) )
       B = 1.0

!  Output solution

       plhs( 1 ) = galmxCreateRealMatrix( n, nb )
       x_pr = mxGetPr( plhs( 1 ) )
       CALL galmxCopyRealToPtr( B, x_pr, n )     

      RETURN
      END SUBROUTINE mexFunction


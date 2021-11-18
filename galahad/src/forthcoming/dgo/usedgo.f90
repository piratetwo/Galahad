! THIS VERSION: GALAHAD 3.3 - 03/07/2021 AT 15:15 GMT.

!-*-*-*-*-*-*-*-*-*-  G A L A H A D   U S E _ D G O  -*-*-*-*-*-*-*-*-*-*-

!  Nick Gould, for GALAHAD productions
!  Copyright reserved
!  July 3rd 2021

   MODULE GALAHAD_USEDGO_double

!  This is the driver program for running DGO for a variety of computing
!  systems. It opens and closes all the files, allocate arrays, reads and
!  checks data, and calls the appropriate minimizers

!    USE GALAHAD_CLOCK
     USE GALAHAD_DGO_double
     USE GALAHAD_SYMBOLS
     USE GALAHAD_SPECFILE_double
     USE GALAHAD_COPYRIGHT
     USE GALAHAD_SPACE_double
     USE GALAHAD_CUTEST_FUNCTIONS_double
     IMPLICIT NONE

     PRIVATE
     PUBLIC :: USE_DGO

   CONTAINS

!-*-*-*-*-*-*-*-*-*-  U S E _ D G O   S U B R O U T I N E  -*-*-*-*-*-*-*-

     SUBROUTINE USE_DGO( input )

!  Dummy argument

     INTEGER, INTENT( IN ) :: input

!  Set precision

     INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )

!-------------------------------
!   D e r i v e d   T y p e s
!-------------------------------

     TYPE ( DGO_control_type ) :: control
     TYPE ( DGO_inform_type ) :: inform
     TYPE ( DGO_data_type ) :: data
     TYPE ( NLPT_problem_type ) :: nlp
     TYPE ( NLPT_userdata_type ) :: userdata
     TYPE ( CUTEST_FUNCTIONS_control_type ) :: cutest_control
     TYPE ( CUTEST_FUNCTIONS_inform_type ) :: cutest_inform

!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------

!  Problem input characteristics

     INTEGER :: iores, i, j, ir, ic, l
     LOGICAL :: filexx, is_specfile
!    REAL :: timeo, timet
!    REAL ( KIND = wp ) :: clocko, clockt
     CHARACTER ( LEN = 10 ) :: name

!  Functions

!$    INTEGER :: OMP_GET_MAX_THREADS

!  Specfile characteristics

     INTEGER, PARAMETER :: input_specfile = 34
     INTEGER, PARAMETER :: lspec = 29
     CHARACTER ( LEN = 16 ) :: specname = 'RUNDGO'
     TYPE ( SPECFILE_item_type ), DIMENSION( lspec ) :: spec
     CHARACTER ( LEN = 16 ) :: runspec = 'RUNDGO.SPC'

!  Default values for specfile-defined parameters

     INTEGER :: dfiledevice = 26
     INTEGER :: rfiledevice = 47
     INTEGER :: sfiledevice = 62
     INTEGER :: wfiledevice = 59
     LOGICAL :: fulsol = .FALSE.
     LOGICAL :: write_problem_data   = .FALSE.
     LOGICAL :: write_solution       = .FALSE.
!    LOGICAL :: write_result_summary = .FALSE.
     LOGICAL :: write_result_summary = .TRUE.
     CHARACTER ( LEN = 30 ) :: dfilename = 'DGO.data'
     CHARACTER ( LEN = 30 ) :: rfilename = 'DGORES.d'
     CHARACTER ( LEN = 30 ) :: sfilename = 'DGOSOL.d'
     CHARACTER ( LEN = 30 ) :: wfilename = 'DGOSAVE.d'
     LOGICAL :: testal = .FALSE.
     LOGICAL :: dechk  = .FALSE.
     LOGICAL :: dechke = .FALSE.
     LOGICAL :: dechkg = .FALSE.
     LOGICAL :: not_fatal  = .FALSE.
     LOGICAL :: not_fatale = .FALSE.
     LOGICAL :: not_fatalg = .FALSE.
     LOGICAL :: getsca = .FALSE.
     INTEGER :: print_level_scaling = 0
     LOGICAL :: scale  = .FALSE.
     LOGICAL :: scaleg = .FALSE.
     LOGICAL :: scalev = .FALSE.
     LOGICAL :: get_max = .FALSE.
     LOGICAL :: warm_start = .FALSE.
     INTEGER :: istore = 0

!  Output file characteristics

     INTEGER :: out  = 6
     INTEGER :: errout = 6
     CHARACTER ( LEN =  6 ) :: solv = 'DGO   '

!  ------------------ Open the specfile for dgo ----------------

     INQUIRE( FILE = runspec, EXIST = is_specfile )
     IF ( is_specfile ) THEN
       OPEN( input_specfile, FILE = runspec, FORM = 'FORMATTED', STATUS = 'OLD')

!   Define the keywords

       spec( 1 )%keyword  = 'write-problem-data'
       spec( 2 )%keyword  = 'problem-data-file-name'
       spec( 3 )%keyword  = 'problem-data-file-device'
       spec( 4 )%keyword  = 'print-full-solution'
       spec( 5 )%keyword  = 'write-solution'
       spec( 6 )%keyword  = 'solution-file-name'
       spec( 7 )%keyword  = 'solution-file-device'
       spec( 8 )%keyword  = 'write-result-summary'
       spec( 9 )%keyword  = 'result-summary-file-name'
       spec( 10 )%keyword = 'result-summary-file-device'
       spec( 11 )%keyword = 'check-all-derivatives'
       spec( 12 )%keyword = 'check-derivatives'
       spec( 13 )%keyword = 'check-element-derivatives'
       spec( 14 )%keyword = 'check-group-derivatives'
       spec( 15 )%keyword = 'ignore-derivative-bugs'
       spec( 16 )%keyword = 'ignore-element-derivative-bugs'
       spec( 17 )%keyword = 'ignore-group-derivative-bugs'
       spec( 18 )%keyword = 'get-scaling-factors'
       spec( 19 )%keyword = 'scaling-print-level'
       spec( 20 )%keyword = 'use-scaling-factors'
       spec( 21 )%keyword = 'use-constraint-scaling-factors'
       spec( 22 )%keyword = 'use-variable-scaling-factors'
       spec( 23 )%keyword = 'maximizer-sought'
       spec( 24 )%keyword = 'restart-from-previous-point'
       spec( 25 )%keyword = 'restart-data-file-name'
       spec( 26 )%keyword = 'restart-data-file-device'
       spec( 27 )%keyword = 'save-data-for-restart-every'
       spec( 28 )%keyword = ''
       spec( 29 )%keyword = ''

!   Read the specfile

       CALL SPECFILE_read( input_specfile, specname, spec, lspec, errout )

!   Interpret the result

       CALL SPECFILE_assign_logical( spec( 1 ), write_problem_data, errout )
       CALL SPECFILE_assign_string ( spec( 2 ), dfilename, errout )
       CALL SPECFILE_assign_integer( spec( 3 ), dfiledevice, errout )
       CALL SPECFILE_assign_logical( spec( 4 ), fulsol, errout )
       CALL SPECFILE_assign_logical( spec( 5 ), write_solution, errout )
       CALL SPECFILE_assign_string ( spec( 6 ), sfilename, errout )
       CALL SPECFILE_assign_integer( spec( 7 ), sfiledevice, errout )
       CALL SPECFILE_assign_logical( spec( 8 ), write_result_summary, errout )
       CALL SPECFILE_assign_string ( spec( 9 ), rfilename, errout )
       CALL SPECFILE_assign_integer( spec( 10 ), rfiledevice, errout )
       CALL SPECFILE_assign_logical( spec( 11 ), testal, errout )
       CALL SPECFILE_assign_logical( spec( 12 ), dechk, errout )
       CALL SPECFILE_assign_logical( spec( 13 ), dechke, errout )
       CALL SPECFILE_assign_logical( spec( 14 ), dechkg, errout )
       CALL SPECFILE_assign_logical( spec( 15 ), not_fatal, errout )
       CALL SPECFILE_assign_logical( spec( 16 ), not_fatale, errout )
       CALL SPECFILE_assign_logical( spec( 17 ), not_fatalg, errout )
       CALL SPECFILE_assign_logical( spec( 18 ), getsca, errout )
       CALL SPECFILE_assign_integer( spec( 19 ), print_level_scaling, errout )
       CALL SPECFILE_assign_logical( spec( 20 ), scale, errout )
       CALL SPECFILE_assign_logical( spec( 21 ), scaleg, errout )
       CALL SPECFILE_assign_logical( spec( 22 ), scalev, errout )
       CALL SPECFILE_assign_logical( spec( 23 ), get_max, errout )
       CALL SPECFILE_assign_logical( spec( 24 ), warm_start, errout )
       CALL SPECFILE_assign_string ( spec( 25 ), wfilename, errout )
       CALL SPECFILE_assign_integer( spec( 26 ), wfiledevice, errout )
       CALL SPECFILE_assign_integer( spec( 27 ), istore, errout )
     END IF

     IF ( dechk .OR. testal ) THEN ; dechke = .TRUE. ; dechkg = .TRUE. ; END IF
     IF ( not_fatal ) THEN ; not_fatale = .TRUE. ; not_fatalg = .TRUE. ; END IF
     IF ( scale ) THEN ; scaleg = .TRUE. ; scalev = .TRUE. ; END IF

!  If required, open a file for the results

     IF ( write_result_summary ) THEN
       INQUIRE( FILE = rfilename, EXIST = filexx )
       IF ( filexx ) THEN
          OPEN( rfiledevice, FILE = rfilename, FORM = 'FORMATTED',             &
                STATUS = 'OLD', POSITION = 'APPEND', IOSTAT = iores )
       ELSE
          OPEN( rfiledevice, FILE = rfilename, FORM = 'FORMATTED',             &
                STATUS = 'NEW', IOSTAT = iores )
       END IF
       IF ( iores /= 0 ) THEN
         write( errout, 2030 ) iores, rfilename
         STOP
       END IF
       READ( INPUT, "( /, I2, A8  )" ) iores, nlp%pname
       REWIND( input )
       WRITE( rfiledevice, "( A10 )" ) nlp%pname
     END IF

!  Set copyright

     IF ( out > 0 ) CALL COPYRIGHT( out, '2021' )

!  Set up control parameters prior to the next solution

     CALL DGO_initialize( data, control, inform )
     IF ( is_specfile ) CALL DGO_read_specfile( control, input_specfile )

!  Initialize the problem data

     cutest_control%input = input ; cutest_control%error = control%error
     CALL CUTEST_initialize( nlp, cutest_control, cutest_inform, userdata )

!  Read a previous solution file for a re-entry

     IF ( warm_start .AND. wfiledevice > 0 ) THEN
       OPEN( wfiledevice, FILE = wfilename, FORM = 'FORMATTED',                &
             STATUS = 'OLD', IOSTAT = iores )
       IF ( iores /= 0 ) THEN
         WRITE( out, 2030 ) iores, wfilename
         STOP
       END IF

       REWIND( wfiledevice )
       READ( wfiledevice, "( A10 )" ) name
       IF ( name /= nlp%pname ) THEN
         WRITE( out, "( /, ' *** Exit from usedgo: re-entry requested with',   &
        &         ' data for problem ', A10, /,                                &
        &         '     but the most recently decoded problem is ', A10 )" )   &
          name, nlp%pname
         STOP
       END IF
       READ( wfiledevice, "( I8 )" )i
       IF ( i /= nlp%n ) THEN
         WRITE( out, "( /, ' *** Exit from usedgo: number of variables'      , &
        &        ' changed from ', I0,' to ', I0,' on re-entry ' )" ) nlp%n, i
         STOP
       END IF
       READ( wfiledevice, "( I8 )" )i
       IF ( i /= nlp%m ) THEN
         WRITE( out, "( /, ' *** Exit from usedgo: number of residuals',       &
        &        ' changed from ', I0, ' to ', I0, ' on re-entry ' )" ) nlp%m, i
         STOP
       END IF
       DO i = 1, nlp%n
         READ( wfiledevice, "( ES24.16, 2X, A10 )" ) nlp%X( i ), name
         IF ( name /= nlp%VNAMES( i ) ) THEN
           WRITE( out, "( /, ' *** Exit from usedgo: variable named ',         &
          &  A, ' out of order on re-entry ' )" ) TRIM( name )
           STOP
         END IF
       END DO
       CLOSE( wfiledevice )
     END IF

!  Solve the problem

     inform%status = 1
!    CALL CPU_TIME( timeo ) ; CALL CLOCK_time( clocko )
     CALL DGO_solve( nlp, control, inform, data, userdata,                     &
                     eval_F = CUTEST_eval_F, eval_G = CUTEST_eval_G,           &
                     eval_H = CUTEST_eval_H, eval_HPROD = CUTEST_eval_HPROD,   &
                     eval_SHPROD = CUTEST_eval_SHPROD )
!    CALL CPU_TIME( timet ) ; CALL CLOCK_time( clockt )

!$    WRITE( out, "( ' number of threads = ', I0 )" ) OMP_GET_MAX_THREADS( )

!  If required, append results to a file,

      IF ( write_result_summary ) THEN
        BACKSPACE( rfiledevice )
        IF ( inform%status == 0 ) THEN
          WRITE( rfiledevice, 2050 ) nlp%pname, nlp%n, inform%obj,             &
            inform%norm_pg, inform%f_eval, inform%g_eval,                      &
            inform%time%clock_total, inform%status
        ELSE
          WRITE( rfiledevice, 2050 ) nlp%pname, nlp%n, inform%obj,             &
            inform%norm_pg, - inform%f_eval, - inform%g_eval,                  &
            - inform%time%clock_total, inform%status
        END IF
      END IF

!  If required, write the solution

      l = 2
      IF ( fulsol ) l = nlp%n
      IF ( control%print_level >= 10 ) l = nlp%n

      WRITE( errout, 2000 )
      DO j = 1, 2
        IF ( j == 1 ) THEN
          ir = 1 ; ic = MIN( l, nlp%n )
        ELSE
          IF ( ic < nlp%n - l ) WRITE( errout, 2010 )
          ir = MAX( ic + 1, nlp%n - ic + 1 ) ; ic = nlp%n
        END IF
        DO i = ir, ic
          WRITE( errout, 2020 ) i, nlp%VNAMES( i ), nlp%X( i ), nlp%X_l( i ),  &
            nlp%X_u( i ), nlp%G( i )
        END DO
      END DO

        WRITE( errout, "( /, 'name           n  f               du-feas ',     &
       &  '      #f       #g     time stat' )" )
        IF ( inform%status == 0 ) THEN
          WRITE( errout, 2050 ) nlp%pname, nlp%n, inform%obj,                  &
            inform%norm_pg, inform%f_eval, inform%g_eval,                      &
            inform%time%clock_total, inform%status
        ELSE
          WRITE( errout, 2050 ) nlp%pname, nlp%n, inform%obj,                  &
            inform%norm_pg, - inform%f_eval, - inform%g_eval,                  &
            - inform%time%clock_total, inform%status
        END IF

      IF ( write_solution .AND.                                                &
          ( inform%status == 0  .OR. inform%status == - 10 ) ) THEN

        INQUIRE( FILE = sfilename, EXIST = filexx )
        IF ( filexx ) THEN
           OPEN( sfiledevice, FILE = sfilename, FORM = 'FORMATTED',            &
               STATUS = 'OLD', IOSTAT = iores )
        ELSE
           OPEN( sfiledevice, FILE = sfilename, FORM = 'FORMATTED',            &
                STATUS = 'NEW', IOSTAT = iores )
        END IF
        IF ( iores /= 0 ) THEN
          write( out, 2030 ) iores, sfilename
          STOP
        END IF

        WRITE( sfiledevice, "( /, ' Problem:    ', A10, /, ' Solver :   ', A,  &
       &       /, ' Objective:', ES24.16 )" ) nlp%pname, solv, inform%obj

        WRITE( sfiledevice, 2000 )
        DO i = 1, nlp%n
          WRITE( sfiledevice, 2020 ) i, nlp%VNAMES( i ), nlp%X( i ),           &
            nlp%X_l( i ), nlp%X_u( i ), nlp%G( i )
        END DO

     END IF

!  if required, save the solution for a future re-entry

     IF ( wfiledevice > 0 .AND. ( ( istore > 0 .AND.                           &
            ( inform%status == GALAHAD_ok .OR.                                 &
              inform%status == GALAHAD_error_ill_conditioned .OR.              &
              inform%status == GALAHAD_error_tiny_step .OR.                    &
              inform%status == GALAHAD_error_max_iterations .OR.               &
              inform%status == GALAHAD_error_time_limit .OR.                   &
              inform%status == GALAHAD_error_alive ) ) .OR.                    &
            inform%status == GALAHAD_error_alive ) ) THEN
       INQUIRE( FILE = wfilename, OPENED = filexx )
       IF ( .NOT. filexx ) THEN
         INQUIRE( FILE = wfilename, EXIST = filexx )
         IF ( filexx ) THEN
            OPEN( wfiledevice, FILE = wfilename, FORM = 'FORMATTED',           &
                  STATUS = 'OLD', POSITION = 'APPEND', IOSTAT = iores )
         ELSE
            OPEN( wfiledevice, FILE = wfilename, FORM = 'FORMATTED',           &
                  STATUS = 'NEW', IOSTAT = iores )
         END IF
         IF ( iores /= 0 ) THEN
           WRITE( out, 2030 ) iores, wfilename
           STOP
         END IF
       END IF
       REWIND( wfiledevice )
       WRITE( wfiledevice, "( A10, ' problem name ', /,                        &
      &                       I8,  ' number of variables', /,                  &
      &                       I8,  ' number of residuals',                     &
      &                       ' ... and now variables' )" )                    &
       nlp%pname, nlp%n, nlp%m
       DO i = 1, nlp%n
         WRITE( wfiledevice, "( ES24.16, 2X, A10 )" )                          &
          nlp%X( i ), nlp%VNAMES( i )
       END DO
       CLOSE( wfiledevice )
     END IF

!  Close any opened files and deallocate arrays

     IF ( is_specfile ) CLOSE( input_specfile )
     CALL CUTEST_terminate( nlp, cutest_inform, userdata )
     RETURN

!  Non-executable statements

 2000 FORMAT( ' Solution: ', /,'                        ',                     &
              '        <------ Bounds ------> ', /                             &
              '      # name          value   ',                                &
              '    Lower       Upper       Dual ' )
 2010 FORMAT( 6X, '. .', 9X, 4( 2X, 10( '.' ) ) )
 2020 FORMAT( I7, 1X, A10, 4ES12.4 )
 2030 FORMAT( ' IOSTAT = ', I6, ' when opening file ', A9, '. Stopping ' )
 2050 FORMAT( A10, I6, ES16.8, ES9.1, bn, 2I9, F9.2, I5 )

!  End of subroutine USE_DGO

     END SUBROUTINE USE_DGO

!  End of module USEDGO_double

   END MODULE GALAHAD_USEDGO_double

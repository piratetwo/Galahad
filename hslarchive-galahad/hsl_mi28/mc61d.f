C COPYRIGHT (c) 1998 Council for the Central Laboratory
*                    of the Research Councils

C Version 1.1.0. See ChangeLog for history

      SUBROUTINE MC61ID(ICNTL,CNTL)

C ICNTL - (OUT) INTEGER array of length 10. On  exit,
C         ICNTL contains default values.
C         If the user wishes to use values other
C         than the defaults, the corresponding entries in
C         ICNTL should be reset after the call to
C         MC61I/ID.

C ICNTL(1) is the stream number for error messages and has the
C          default value 6. Printing of error messages is
C          suppressed if ICNTL(1) < 0.

C ICNTL(2) is the stream number for warning messages.
C          It has the default value 6. Printing of warning
C          messages is suppressed if ICNTL(2) < 0.

C ICNTL(3) controls the action taken if duplicate or out-of-range
C          entries are detected . If ICNTL(3) = 0 and such entries are
C          detected, the computation terminates
C          with IRN and ICPTR unchanged. If ICNTL(3) = 1, a warning
C          is issued and the computation continues.
C          The default value is 0.

C ICNTL(4) controls whether supervariables are to be used
C          (a supervariable is a set of variables that
C          correspond to a set of identical columns).  If
C          ICNTL(4) = 0, supervariables are used. If
C          ICNTL(4) = 1, variables are used.  If the problem has
C          significantly fewer supervariables than variables,
C          using supervariables will reduce the execution time
C          significantly and will, in general, produce a permutation
C          or assembly order of comparable quality.
C          The default value is 0.

C ICNTL(5) indicates whether the user wishes to supply
C          a global priority function in PERM.
C          If ICNTL(5) = 0, no priority function is supplied;
C          if ICNTL(5) = 1, a priority function is supplied.
C          The default value is 0.

C ICNTL(6) indicates whether the user wishes to supply
C          the weights for the priority function.
C          If ICNTL(6) = 0, no weights are supplied
C                           and (2,1) and (16,1) tried;
C          If ICNTL(6) = 1, no weights are supplied
C                           and (1,2) and (16,1) tried;
C          if ICNTL(6) = 2, weights are supplied in CNTL(1), CNTL(2).
C          The default value is 0.

C ICNTL(7) to ICNTL(10)  are given default values of zero.
C          They are currently not used but may be used in a
C          later release of the code.

C CNTL  - (OUT) REAL (DP) array of length 5. On  exit,
C         CNTL contains default values.
C         If the user wishes to use values other
C         than the defaults, the corresponding entries in
C         CNTL should be reset after the call to
C         MC61I/ID.
C         CNTL holds the weights in the priority function.
C         They are given the default values 2.0 and 1.0.
C         They are only used if ICNTL(6) = 2.
C CNTL(3) to CNTL(5)  are given default values of zero.
C          They are currently not used but may be used in a
C          later release of the code.


      DOUBLE PRECISION ZERO
      PARAMETER (ZERO = 0.0D0)
      DOUBLE PRECISION CNTL(5)
      INTEGER ICNTL(10)
      INTEGER I

      DO 10 I = 1,10
         ICNTL(I) = 0
   10 CONTINUE
      ICNTL(1) = 6
      ICNTL(2) = 6

      CNTL(1) = 2.0D0
      CNTL(2) = 1.0D0
      CNTL(3) = ZERO
      CNTL(4) = ZERO
      CNTL(5) = ZERO

      RETURN
      END



      SUBROUTINE MC61AD(JOB,N,LIRN,IRN,ICPTR,PERM,LIW,IW,W,ICNTL,
     +                  CNTL,INFO,RINFO)

C JOB   - (IN) INTEGER variable that must be set by the user to 1 if
C         a variable permutation to reduce the profile and wavefront of
C         the matrix is required, to 2 if a variable permutation to
C         reduce the bandwidth is required, and to 3 if an assembly
C         order is required.

C N     - (IN) INTEGER variable that must be set by the user to the
C         order of the matrix.

C LIRN  - (IN) INTEGER variable that must be set by the user to the
C         length of the array IRN, which must be large enough to
C         hold the pattern of the whole matrix.

C IRN   - (INOUT) INTEGER  array of length LIRN whose leading part must
C         be set by the user to hold the row indices of the
C         entries in the lower triangle of the matrix,
C         including the diagonal. The entries of each column must
C         be contiguous. The entries of column J must precede those
C         of column J+1 (J=1,...,N-1), and there must be no wasted
C         space between the columns.  Row indices within a
C         column may be in any order. On exit, holds the row
C         entries of the condensed matrix, using the same format.

C ICPTR - (INOUT) INTEGER array of length N+1 that must be set by the
C         user so that ICPTR(J) points to the position in the array
C         IRN of the first entry in column J (J=1,...,N), and
C         ICPTR(N+1)-1  must be the position of the last entry.
C         On exit, ICPTR holds corresponding data for the condensed
C         matrix.

C PERM  - (INOUT) INTEGER array of length N.  This array must only
C         be set on entry if ICNTL(5) = 1 (the default is ICNTL(5) = 0.
C         In this case, this array must be set by the user to
C         hold the global priority function.  On exit, the new ordering
C         is contained in PERM.  If a variable permutation is
C         requested (JOB = 1 or 2), the new index for variable
C         I is given by PERM(I) (I = 1,...,N).  If an
C         assembly order is requested (JOB = 3), the order
C         in which the rows should be assembled is
C         PERM(1), PERM(2),...,PERM(N).

C LIW   - (IN) INTEGER variable that must be set by the user to the
C         length of the array IW. Sufficient value is LIW = 8*N+2
C         (7N+2 if ICNTL(5)=1 and 6N+2 if ICNTL(5)=1, ICNTL(6)=2).

C IW    - (OUT) INTEGER array of length LIW that is used by
C         the routine as workspace.

C W     - (OUT) REAL (DP) array of length N that is used by
C         the routine as workspace.

C ICNTL - (IN) INTEGER array of length 10. Holds values of control
C         parameters. Default values may be set by calling
C         MC61I/ID.

C CNTL  - (IN) REAl (DP) array of length 5. Holds values of control
C         parameters. Default values may be set by calling
C         MC61I/ID.

C INFO  - (OUT) INTEGER  array of length 10 that need not be set by the
C         user. On each successful exit, INFO(1) is set to 0.
C         Negative values of INFO(1) indicate a fatal error has
C         been detected and positive values indicate a warning
C         has been issued.
C INFO(1) = -1 JOB is not equal to 1, 2 or 3.
C            Immediate return with input parameters unchanged.
C INFO(1) = -2  N < 1.  Immediate return with input parameters
C            unchanged.
C INFO(1) = -3  LIRN is too small.  INFO(6)  is set to the minimum value
C           which will suffice for LIRN.
C            If LIRN is at least as large as the input value
C            of ICPTR(N+1)-1 and ICNTL(3) = 1,
C            any out-of-range or duplicated variable indices
C            will have been excluded from IRN and ICPTR.
C            Otherwise, the input parameters  are unchanged.
C INFO(1) = -4 LIW is too small.  INFO(3)  is set to a value
C            which may suffice for LIW.
C INFO(1) = -5  One or more variable indices either lies outside the
C            lower triangle  of the matrix or is duplicated.
C            Further information is contained in
C            INFO(4)> and  INFO(5).
C
C INFO(2) holds, on successful exit, the total number of
C         supervariables in the problem.  If variables are not used
C         (ICNTL(4) = 1), INFO(3) is set to N.
C INFO(3) holds the amount of workspace used by the routine. If the
C         user has provided insufficient workspace (INFO(1) = -4),
C        INFO(3) is set to a value which may suffice for LIW.
C INFO(4) holds the number of variable indices
C         in IRN found to be out-of-range.
C INFO(5) holds the number of duplicate variable indices in IRN.
C INFO(6) holds the minimum value which will suffice for LIRN.
C INFO(7) holds the number of non-trivial components
C         in the graph of the matrix
C INFO(8) to INFO(10)  are currently not used but may be used in a
C later release of the code.


C RINFO - (OUT) REAL (DP) array of length 15 that need not be set by the
C          user.

C If JOB= 1 or 2,  on successful exit RINFO
C returns the following information:
C RINFO(1) holds the profile of the matrix A.
C RINFO(2) holds the maximum wavefront of the matrix A.
C RINFO(3) holds the bandwidth of the matrix A.
C RINFO(4) holds the root mean squared wavefront of the matrix A.
C RINFO(5) holds the profile of the permuted  matrix.
C RINFO(6) holds the maximum wavefront of the permuted matrix.
C RINFO(7) holds the bandwidth  of the permuted matrix.
C RINFO(8) holds the root mean squared wavefront of the permuted matrix.

C If JOB = 3,  on successful exit RINFO
C returns the following information:
C RINFO(1) holds the maximum row frontsize for
C the assembly order 1, 2, ..., N.
C RINFO(2) holds the maximum column frontsize for
C the assembly order 1, 2, ..., N.
C RINFO(3) holds the root mean square row front size
C for the assembly order 1, 2, ..., N.
C RINFO(4) holds the mean frontal matrix size
C for the assembly order 1, 2, ..., N.
C RINFO(5) holds the maximum row frontsize for
C the new assembly order PERM(1), PERM(2), ..., PERM(N).
C RINFO(6) holds the maximum column frontsize for
C the new assembly order PERM(1), PERM(2), ..., PERM(N).
C RINFO(7) holds the root mean square row front size
C for the new assembly order PERM(1), PERM(2), ..., PERM(N).
C RINFO(8) holds the mean frontal matrix size
C for the new assembly order PERM(1), PERM(2), ..., PERM(N).

C If JOB = 1 or 3,
C RINFO(9) and RINFO(10) hold weights used.

C RINFO(11) to RINFO(15)  are currently not used but may be used in a
C later release of the code.

      DOUBLE PRECISION ZERO,ONE,TWO,SIXTN
      PARAMETER (ZERO = 0.0D0, ONE=1.0D0, TWO=2.0D0, SIXTN=16.0D0)

      INTEGER JOB,N,LIW,LIRN

      DOUBLE PRECISION RINFO(15)
      DOUBLE PRECISION CNTL(5),W(N)
      INTEGER IRN(LIRN),ICPTR(N+1),INFO(10),ICNTL(10),IW(LIW),PERM(N)

      DOUBLE PRECISION RNFO5,RNFO6,RNFO7,RNFO8
      INTEGER I,PERMSV,IWORK,LP,MP,NSUP,SVAR,VARS,IRUN,NRUN,COPY,
     +        IPERM,J,ISUP,JSUP,PAIR
      LOGICAL LSWAP

      INTEGER ICON60(2),INFO60(4),JCNTL(2)

      EXTERNAL MC60AD,MC60BD,MC60CD,MC60DD,MC60ED,MC60FD,MC60GD

C Streams for messages
      LP = ICNTL(1)
      MP = ICNTL(2)
C Initialise
      DO 10 I = 1,10
         INFO(I) = 0
   10 CONTINUE
      DO 20 I = 1,15
         RINFO(I) = ZERO
   20 CONTINUE

C Check input data for errors

      IF (JOB.NE.1 .AND. JOB.NE.2 .AND. JOB.NE.3) GO TO 100
      IF (N.LT.1) GO TO 110
      IF (LIRN.LT.ICPTR(N+1)-1) GO TO 115

C Call MC60A/AD to construct pattern of whole matrix.
C MC60A/AD also checks data for errors.
C We do not want messages to be printed from MC60A/AD.
      ICON60(1) = ICNTL(3)
      ICON60(2) = 0

C Check workspace. The workspace has to be at least 5*N+2.
      INFO(3) = 5*N + 2
      IF (5*N+2.GT.LIW) THEN
C Not even minimum required ... return a sufficient value.
         INFO(3) = 8*N + 2
         IF (JOB.NE.2) THEN
           IF (ICNTL(5).EQ.1) INFO(3) = 7*N + 2
           IF (ICNTL(5).EQ.1 .AND. ICNTL(6).EQ.2) INFO(3) = 6*N + 2
         END IF
         GO TO 130
      END IF

      CALL MC60AD(N,LIRN,IRN,ICPTR,ICON60,IW,INFO60)

C Copy information returned from MC60A/AD into INFO
      INFO(4) = INFO60(2)
      INFO(5) = INFO60(3)
      INFO(6) = INFO60(4)

C Check for errors
C Errors  -2 and -3 are possible.
      IF (INFO60(1).EQ.-2) GO TO 120
      IF (INFO60(1).EQ.-3) GO TO 140

C Check for warnings (only possible if ICNTL(3)=1).
      IF (INFO60(1).GT.0) THEN
C Duplicated or out-of-range entries have been found.
         INFO(1) = INFO60(1)
         IF (MP.GT.0) THEN
            IF (INFO(4).GT.0) WRITE (MP,'(/,A,I8,A)')
     * ' MC61A/AD warning:',INFO(4),' out-of-range entries are ignored'
            IF (INFO(5).GT.0) WRITE (MP,'(/,A,I8,A)')
     * ' MC61A/AD warning:',INFO(5),' duplicated entries are ignored'
         END IF
      END IF

C Put initial permutation/assembly order into IW
C (we can't put it into PERM as we must not overwrite any global
C priority function supplied by the user).
C We do not use supervariables initially (so need VARS(I) = 1
C in MC60F/FD).
      SVAR = 1
      VARS = SVAR + N
      IPERM = VARS + N
      IWORK = IPERM + N
C IWORK is length 2*N+1 ... already checked sufficient space.
      DO 30 I = 1, N
         IW(VARS+I-1) = 1
         IW(IPERM+I-1) = I
   30 CONTINUE

      IF (JOB.EQ.1 .OR. JOB.EQ.2) THEN
C Use MC60F/FD to compute the profile etc for original matrix
         CALL MC60FD(N,N,LIRN,IRN,ICPTR,IW(VARS),IW(IPERM),IW(IWORK),
     *               RINFO)
      ELSE
C Use MC60G/GD to compute statistics for original assembly order
         CALL MC60GD(N,N,LIRN,IRN,ICPTR,IW(VARS),IW(IPERM),IW(IWORK),
     *               RINFO)
      END IF


      IF (ICNTL(4).EQ.0) THEN
C Supervariables are to be used.
C Call MC60B/BD to find supervariables and compress pattern.
         IWORK = VARS + N
C Workspace length 2*N+2 required ...already checked we have
C this
         CALL MC60BD(N,LIRN,IRN,ICPTR,NSUP,IW(SVAR),IW(VARS),IW(IWORK))
      ELSE
C Supervariables are NOT used.
         NSUP = N
      END IF
      INFO(2) = NSUP

C Find permutation for variables or supervariables.
      JCNTL(1) = 0
      IF (JOB.EQ.2) JCNTL(1) = 1

      JCNTL(2) = 0
      IF (ICNTL(5).EQ.1 .AND. JOB.NE.2) JCNTL(2) = 2

C Does the user wish to set the weights?
      IF (JOB.EQ.1 .OR. JOB.EQ.3) THEN
         IF (ICNTL(6).EQ.0)  THEN
C We will try the pairs (2,1) and (16,1) and choose the best
            NRUN = 2
            RINFO(9) = TWO
            RINFO(10) = ONE
         ELSE IF (ICNTL(6).EQ.1)  THEN
C We will try the pairs (1,2) and (16,1) and choose the best
            NRUN = 2
            RINFO(9) = ONE
            RINFO(10) = TWO
         ELSE IF (ICNTL(6).EQ.2) THEN
C Use the weights in CNTL
            NRUN = 1
            RINFO(9) = CNTL(1)
            RINFO(10) = CNTL(2)
         END IF
      ELSE
C No weights needed if JOB = 2 (RCM)
         NRUN = 1
      END IF

C Check workspace
C Workspace must be length 3*NSUP+1.
C We only need an array to copy into if we are going to try
C two sets of weights.
      PERMSV = VARS + NSUP
      IWORK = PERMSV + NSUP
      PAIR = IWORK + 3*NSUP + 1
      COPY = PAIR + NSUP
      IF (ICNTL(5).EQ.1 .AND. JOB.NE.2) COPY = PAIR
      IF (NRUN.EQ.1) INFO(3) = MAX(INFO(3),COPY)
      IF (NRUN.EQ.2) INFO(3) = MAX(INFO(3),COPY+N-1)
      IF (INFO(3).GT.LIW) GO TO 130

      DO 60 IRUN = 1,NRUN

         IF (IRUN.EQ.1) THEN
C If user has supplied global priority function, copy PERM(I) into
C IW(PERMSV+I-1).
C If we are working with supervariables, we have to convert given
C priority function into a supervariable priority function.
            IF (ICNTL(5).EQ.1) THEN
               IF (ICNTL(4).EQ.1) THEN
C Working with variables
                  DO 34 I = 1,N
                     J = PERM(I)
                     IW(PERMSV+I-1) = J
   34             CONTINUE
               ELSE
C Working with supervariables
C SVAR(I) is the supervariable to which variable I belongs.
                  DO 35 I = 1,N
                     J = PERM(I)
                     ISUP = IW(SVAR+I-1)
                     JSUP = IW(SVAR+J-1)
                     IW(PERMSV+ISUP-1) = JSUP
   35             CONTINUE
               END IF
            END IF
            IF (NRUN.EQ.2) THEN
C Take a copy of IW(PERMSV+I-1) for use with the second pair of weights.
               DO 36 I = 1,NSUP
                  IW(COPY+I-1) = IW(PERMSV+I-1)
   36          CONTINUE
            END IF
         END IF

C If we are trying a second pair of weights, we must take a copy
C of the permutation we obtained with the first and reset the weights
         IF (IRUN.EQ.2) THEN
C
            IF (ICNTL(5).EQ.1) THEN
C copy what is in IW(COPY+I-1) into IW(PERMSV+I-1)
               DO 37 I = 1,NSUP
                  IW(PERMSV+I-1) = IW(COPY+I-1)
   37          CONTINUE
            END IF
C Take a copy of PERM
            DO 40 I = 1,N
               IW(COPY+I-1) = PERM(I)
   40       CONTINUE
            RNFO5 = RINFO(5)
            RNFO6 = RINFO(6)
            RNFO7 = RINFO(7)
            RNFO8 = RINFO(8)

            RINFO(9) = SIXTN
            RINFO(10) = ONE

C On second run, the start and end nodes are already known
C (since returned on the first run) ... make use of this
C in the case when global priority function NOT supplied.
            IF (ICNTL(5).NE.1) JCNTL(2) = 1

         END IF

         CALL MC60CD(N,NSUP,LIRN,IRN,ICPTR,IW(VARS),JCNTL,IW(PERMSV),
     *              RINFO(9),IW(PAIR),INFO60,IW(IWORK),W)
         INFO(7) = INFO60(1)

         IF (JOB.EQ.1 .OR. JOB.EQ.2) THEN
C Variable permutation required.

C If supervariables have been used, find permutation for
C variables from supervariable permutation
            IF (ICNTL(4).EQ.0) THEN
               CALL MC60DD(N,NSUP,IW(SVAR),IW(VARS),IW(PERMSV),
     +                    PERM,IW(IWORK))
            ELSE
               DO 50 I = 1,N
                  PERM(I) = IW(PERMSV+I-1)
   50          CONTINUE
            END IF

C Use MC60F/FD to compute the profile etc for the permuted matrix
            CALL MC60FD(N,NSUP,LIRN,IRN,ICPTR,IW(VARS),IW(PERMSV),
     *                 IW(IWORK),RINFO(5))
         ELSE
C Assembly order required.
C Find row assembly order from (super)variable permutation
            CALL MC60ED(N,NSUP,LIRN,IRN,ICPTR,IW(SVAR),IW(VARS),
     *                 IW(PERMSV),PERM,IW(IWORK))
C Use MC60G/GD to compute statistics for new assembly order
C (compressed graph used for this).
            CALL MC60GD(N,NSUP,LIRN,IRN,ICPTR,IW(VARS),IW(PERMSV),
     *                 IW(IWORK),RINFO(5))
         END IF

   60 CONTINUE

C If we have done 2 runs, we now have to choose the better
C set of weights ... store in RINFO(9) and RINFO(10)
       LSWAP = .FALSE.
       IF (NRUN.EQ.2) THEN
          IF (JOB.EQ.1) THEN
C choose on basis of profile
             IF (RNFO5.LT.RINFO(5)) LSWAP = .TRUE.
          ELSE IF (JOB.EQ.3) THEN
C choose on basis of mean frontal matrix size
             IF (RNFO8.LT.RINFO(8)) LSWAP = .TRUE.
          END IF
          IF (LSWAP) THEN
C First set of weights gave best result
             RINFO(9) = TWO
             RINFO(10) = ONE
             IF (ICNTL(6).EQ.1)  THEN
                RINFO(9) = ONE
                RINFO(10) = TWO
             END IF
             DO 70 I = 1,N
               PERM(I) = IW(COPY+I-1)
   70       CONTINUE
            RINFO(5) = RNFO5
            RINFO(6) = RNFO6
            RINFO(7) = RNFO7
            RINFO(8) = RNFO8
         END IF
      END IF

C Issue a warning if new ordering is not better than old.
          IF (JOB.EQ.1 .AND. RINFO(1).LE.RINFO(5)) THEN
C choose on basis of profile
             IF (MP.GT.0) WRITE (MP,'(/,A)')
     *       ' MC61A/AD warning: Profile not reduced'
             INFO(1) = 1
          ELSE IF (JOB.EQ.2 .AND. RINFO(3).LE.RINFO(7)) THEN
C choose on basis of bandwidth
             IF (MP.GT.0) WRITE (MP,'(/,A)')
     *       ' MC61A/AD warning: Semibandwidth not reduced'
             INFO(1) = 1
          ELSE IF (JOB.EQ.3 .AND. RINFO(4).LE.RINFO(8)) THEN
C choose on basis of mean frontal size
             IF (MP.GT.0) WRITE (MP,'(/,A)')
     * ' MC61A/AD warning: mean frontal matrix size not reduced'
             INFO(1) = 1
          END IF

      GO TO 200
C Error returns

  100 INFO(1) = -1
      IF (LP.GT.0) THEN
         WRITE (LP,9000) INFO(1)
         WRITE (LP,9010) JOB
      END IF
      GO TO 200

  110 INFO(1) = -1
      IF (LP.GT.0) THEN
         WRITE (LP,9000) INFO(1)
         WRITE (LP,9020) N
      END IF
      GO TO 200

  115 INFO(1) = -1
      IF (LP.GT.0) THEN
         WRITE (LP,9000) INFO(1)
         WRITE (LP,9025)
      END IF
      GO TO 200

  120 INFO(1) = -2
      IF (LP.GT.0) THEN
         WRITE (LP,9000) INFO(1)
         WRITE (LP,9030) INFO(6)
      END IF
      GO TO 200

  130 INFO(1) = -3
      IF (LP.GT.0) THEN
         WRITE (LP,9000) INFO(1)
         WRITE (LP,9040) INFO(3)
      END IF
      GO TO 200

  140 INFO(1) = -4
      IF (LP.GT.0) THEN
         WRITE (LP,9000) INFO(1)
         IF (INFO(4).GT.0) WRITE (LP,'(A,I6,A)')
     *         ' There are ',INFO(4),' out-of-range entries'
         IF (INFO(5).GT.0) WRITE (LP,'(A,I6,A)')
     *         ' There are ',INFO(5),' duplicated entries'
      END IF
      GO TO 200

  200 RETURN

9000  FORMAT (/' MC61A/AD error:  INFO(1) = ',I3)
9010  FORMAT (' Value of JOB out-of-range.  JOB = ',I8)
9020  FORMAT (' Value of N out-of-range.  N = ',I8)
9025  FORMAT (' Value of LIRN is less than ICPTR(N+1)-1')
9030  FORMAT (' Value of LIRN too small.  Increase to at least = ',I8)
9040  FORMAT (' Value of LIW too small.  Sufficient value = ',I8)
      END



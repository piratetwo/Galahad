C COPYRIGHT (c) 1999 Council for the Central Laboratory
*                    of the Research Councils
C Original date 25 March 1999
C Modified 24 May 2001 to make LA04 threadsafe (LA05 -> LA15)
C 26 June 2001. Arrays W and Y of LA04E/ED made into arguments to
C     conform with Fortran 77. Workspace WS has to be enlarged by M.
C
C 12th July 2004 Version 1.0.0. Version numbering added.
C 1st March 2005 Version 1.1.0. FD05 and ZA02 dependence changed
C     to FD15 and ZA12.
C 28 Jan. 2011 Version 1.2.0. ZA12 replaced by cpu_time and FD15 by
C Fortran 90 intrinsics.

      SUBROUTINE LA04AD(A,LA,IRN,IP,M,N,B,C,BND,KB,LB,JOB,CNTL,IX,JX,X,
     +                  Z,G,RINFO,WS,LWS,IWS,LIWS)

C     .. Scalar Arguments ..
      INTEGER LA,KB,LB,M,N,LWS,LIWS
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION
     +                A(LA),B(M),BND(2,KB),C(N),CNTL(15),G(N),
     +                RINFO(40),WS(LWS),X(N+M),Z(N)
      INTEGER IP(N+1),IRN(LA),IWS(LIWS),IX(M),JOB,JX(KB)

C A   (INOUT) holds the matrix entries ordered by columns.
C      Altered only if JOB=0, when the entries of b and c are
C      added and it is scaled.
C LA  (IN) length of A and IRN.
C IRN (INOUT) holds the row indices of the matrix entries in A.
C      Altered only if JOB=0, when the entries of b and c are
C      added.
C IP  (IN) IP(J) holds the position in A and IRN of the first
C     entry of column J and IP(N+1) holds the first unused position.
C M   (IN) Number of rows.
C N   (IN) Number of columns.
C B   (INOUT) RHS vector. Altered only if JOB=0, when it is scaled.
C C   (INOUT) Objective vector. Altered only if JOB=0, when it is
C     scaled.
C BND (INOUT) BND(1,J) and BND(2,J) hold the bounds for variable J,
C     J = 1,2,...,KB. Altered only if JOB=0, when it is scaled.
C KB  (IN) Second extent of BND.
C LB  (IN) First variable with implicit lower bound 0.
C JOB (INOUT) Possible values on entry:
C     0 scale the problem
C     1 initial entry
C     2 initial entry with basis or partial basis in IX and JX
C     3 normal iteration
C     4 normal iteration with refactorization of basis
C     5 normal iteration with recalculation of solution and
C       reduced costs
C     6 calculate B**-1 WS
C     7 calculate B**-T WS
C  Possible values on exit:
C     0 normal termination
C     1 successful scaling
C     3 normal iteration
C     4 normal iteration, but refactorization required next time
C    >5 normal exit after entry with JOB>5.
C    -1 feasible region is empty
C    -2 objective function is unbounded
C    -3 Problem with M, N, KB, or LB.
C    -4 IB is too small. Must be at least RINFO(35).
C    -5 Problem with BND(:,J), J in RINFO(35).
C    -6 Problem with IP(J), J in RINFO(35).
C    -7 Problem with IRN(K), K in RINFO(35).
C    -8 Duplicate in IRN(K), K in RINFO(35).
C    -9 JOB<0 or JOB>7 on entry
C   -10 LA is too small. Must be at least RINFO(35).
C   -11 Problem with IX(I), I in RINFO(35).
C   -12 Problem with JX(I), I in RINFO(35).
C   -13 LWS is too small. Must be at least RINFO(35).
C   -14 LIWS is too small. Must be at least RINFO(35).
C CNTL(1) (IN) Value that represents an infinite bound.
C CNTL(2) (INOUT) Relative tolerance for the entries of A, B, C.
C CNTL(3) (INOUT) holds a tolerance for the bounds on the variables.
C       It may be increased on a recalculation of the solution.
C CNTL(4) (IN) 1 for steepest-edge scaling, or
C              0 for unscaled reduced costs.
C CNTL(5) (INOUT) < 1 No iterative refinement
C                 >=1 Perform iterative refinement
C CNTL(6) (IN) Unit for error messages.
C CNTL(7) (IN) Unit for monitor output.
C CNTL(8) (IN) Relative pivot tolerance for LA15.
C CNTL(9) (IN) Limit on LA15 growth parameter.
C CNTL(10) (IN) The LA15's CNTL(1) is given the value
C       CNTL(10)*max|(A(K)|.
C CNTL(11) (INOUT) In the vector w of changes to x, entries less than
C       CNTL(11t)*norm(w), for the 2-norm, are treated as zero.
C CNTL(12) Default zero: LA04 initializes the random number seed on the
C   initialization entries when JOB = 0, 1 or 2.  If CNTL(12) is any
C   other value it is assumed the user is controlling the seed which
C   is always held in IWS(IISD), where IISD should equal LIWS.
C CNTL(13) - CNTL(15). Not used
C IX (INOUT) IX(I), I=1,2,...,M, holds the index of column I of the
C   basis. Possible values:
C   IX(I)=0  Missing column (only partial basis specified).
C   IX(I)<=N  Ordinary variable
C   N<IX(I)<=N+M  Artifical variable with lower and upper bounds of zero
C   N+M<IX(I)<=2*(N+M) Variable IX(I)-(N+M) is above its upper bound.
C   2*(N+M)<IX(I) Variable IX(I)-2*(N+M) is below its lower bound.
C JX (INOUT) JX(J), J=1,2,...,KB, holds
C    -1 if variable J is in basis and feasible with value X(J)
C     0 if variable J is out of basis with value X(J)
C     1 if lower bound is active: variable has value X(J)+BND(1,J);
C       it may be in the basis and infeasible or out of the basis.
C     2 if upper bound is active: variable has value X(J)+BND(2,J);
C       it may be in the basis and infeasible or out of the basis.
C X (INOUT) Values of the variables less any active bounds.
C Z (INOUT) If variable J is non-basic, Z(J) holds its reduced cost.
C G (INOUT) G(J), J=1,2,...,N, holds
C   For a variable that is in the basis, -(the position in the basis).
C   For an out-of-basis fixed-value, zero.
C   Otherwise, the steepest-edge weight.
C RINFO (INOUT) holds useful information:
C RINFO(1) holds the objective function or the sum of the modulii
C      of the infeasibilities.
C RINFO(2) holds the parameter U of LA15A/C.
C RINFO(3) contains the parameter G of LA15A/C.
C RINFO(4) holds the c.p.u. time spent in LA04A and subroutines
C      called from it since the last initial (JOB=1,2) entry.
C RINFO(5) holds the c.p.u. time spent in LA15A
C      since the last initial entry.
C RINFO(6) holds the c.p.u. time spent in LA15B
C      since the last initial entry.
C RINFO(7) holds the c.p.u. time spent in LA15C
C      since the last initial entry.
C RINFO(8) contains the c.p.u. time spent on the last refactorization
C      of the basis matrix and recalculation of X and Z.
C RINFO(9) holds the c.p.u. reading on the last LA04A call.
C RINFO(10)  holds the c.p.u. time spent in LA04A and subroutines
C      it calls between the last JOB=1,2 entry and the completion
C      of the last refactorization.
C RINFO(11) to RINFO(20) are used to hold the last ten values of the
C      average iteration time since the last refactorization.
C RINFO(21) holds the minimum average iteration time since the last
C      refactorization.
C RINFO(22) Not used
C RINFO(23) holds the modulus of the largest matrix entry.
C RINFO(24) holds the maximum permitted value for the parameter G of
C      LA15A. If this value is exceeded after an LA15C call, a
C      fresh factorization is made.
C RINFO(25) holds the index of the variable leaving the basis.
C RINFO(26) holds the index of the variable entering the basis.
C RINFO(27) holds the total number of infeasibilities (i.e. the
C      total number of variables outside their bounds by more
C      than CNTL(3)).
C RINFO(28) holds the number of iterations since the last
C       refactorization.
C RINFO(29) holds the number of iterations since the last
C       initial (JOB=1,2) entry.
C RINFO(30) holds the number of refactorizations since the last
C       initial entry.
C RINFO(31) has the value 2 if the current solution is feasible
C       and 1 otherwise.
C RINFO(32) holds the number of entries in the factorized basis
C RINFO(33) holds the no. of compresses of the factorized basis
C       since the previous refactorization.
C RINFO(34) is set to the number of matrix sweeps needed by the
C       crash algorithm
C RINFO(35) holds information on a return with JOB <= -4.
C RINFO(36) to RINFO(40) Not used

C WS  (INOUT) mainly used as a workarray.
C       On a return with JOB = 1, the first N+1 elements hold the
C       column scaling factors and the next M+1 elements hold the
C       row scaling factors.
C       On a return with JOB = -2, the first N elements hold
C       the direction to an unbounded solution.
C       On an entry with JOB =  6 or 7, the first M elements
C                                                  -1      -T
C       hold the vector w on entry and the vector B  w or B  w
C       on return, where B is the basis.
C LWS(IN) Length of array WS
C IWS  workarray
C LIWS(IN) Length of array IWS

C     .. Local Scalars ..
      DOUBLE PRECISION DUM
      INTEGER I,IB,IISD,II15,IK15,IC15,IPB,IR,IROW,IWB,J,JB,JIN,
     +        JV,JY,JYY,K,K1,K2,L,LP,MP,NPM,NZ
      DOUBLE PRECISION
     +                ONE,SIGMA,SIGN,T,ZERO
      PARAMETER (ONE=1,ZERO=0)
C DUM   Used in cpu_time call.
C I     Row index
C IB    Size of the arrays that hold the factorized basis.
C IC15  IWS(IC15+J),J=1,3 used for LA15's argument CNTL(J),J=1,3
C IISD  IWS(IISD) used for random number seed (FA14)
C II15  IWS(II15+J),J=1,3 used for LA15's argument ICNTL(J),J=1,3
C IK15  IWS(IK15+J),J=1,7 used for LA15's argument KEEP(J),J=1,7
C IPB   Position in IWS of the arrays IW of LA04B and IPB of LA04C/D
C IR    Position in IWS of the array IR of LA04D.
C IROW  Position in IWS of the array IROW of LA04B.
C IWB   Position in IWS of the arrays IW of LA04B and IWB of LA04C/D
C J     Column index
C JB    Position in WS of the array BA of LA04C and V of LA04D
C JIN   Incoming column index
C JV    Position in WS of the array PV of LA04B and V of LA04C/D
C JY    Position in WS of the array Y of LA04B/C/D.
C JYY   Position in WS of the array YY of LA04B/C/D.
C K     Index in A and IRN.
C K1    Column start.
C K2    Column end.
C L     Temporary.
C LP    Unit for error messages
C MP    Unit for monitor output
C NPM   N+M
C NZ    No. of entries in A.
C SIGMA Value that represents an infinite bound.
C SIGN  Sign of Z(JIN)
C T     Temporary for holding the c.p.u. time

C     .. External Subroutines ..
      EXTERNAL LA04BD,LA04CD,LA04DD,LA04SD,LA15BD

C     .. Intrinsic Functions ..
      INTRINSIC ABS,MAX

      call cpu_time(DUM)
      RINFO(9) = DUM
      LP = CNTL(6)
      MP = CNTL(7)
      SIGMA = CNTL(1)
      IF (JOB.GT.7 .OR. JOB.LT.0) THEN
        IF (LP.GE.0) WRITE (LP,FMT='(/,/,A,I5)')
     +      ' LA04 error return: JOB has value',JOB
        JOB = -9
        RETURN
      END IF
C
C Partition workspace.
      IB = MIN(LWS-4*M-4, (LIWS-10*M-11)/2)
      JV = 1 + M
      JY = JV + M
      JYY = JY + M
      JB = JYY + M
      II15 = LIWS - 11
      IK15 = II15 + 3
      IISD = LIWS

      IWB = 2*IB + 1
      IPB = IWB + 8*M
      IR = IWB + 4*M
      IROW = IR + M
      IC15 = LWS - 3

      IF (JOB.GE.6) THEN
        IF (M.GT.0) THEN
          call cpu_time(T)
          CALL LA15BD(WS(JB),IWS,IB,M,IWS(IPB),IWS(IWB),WS(JV),RINFO(3),
     +                WS,JOB.EQ.7,IWS(II15+1),IWS(IK15+1))
          call cpu_time(DUM)
          RINFO(6) = RINFO(6) + DUM - T
        END IF
        GO TO 150
      END IF

      IF (JOB.GT.2) GO TO 140

C Initialization of timing parameters,etc.
      RINFO(4) = ZERO
      RINFO(5) = ZERO
      RINFO(6) = ZERO
      RINFO(7) = ZERO
      RINFO(29) = 0
      RINFO(30) = 0
      RINFO(33) = 0
      RINFO(27) = M
      CNTL(2) = MAX(CNTL(2), EPSILON(ONE)*100.0)
      CNTL(11) = MAX(CNTL(11), EPSILON(ONE)*100.0)

      NZ = IP(N+1) - 1
C Print scalar arguments and CNTL
      IF (MP.GE.0) THEN
        WRITE (MP,FMT='(/,3(A,I7)/6(A,I7))') ' Entry to LA04 with M=',M,
     +    ' N=',N,' No. of non-zeros=',NZ,' LA=',LA,' KB=',KB,' LB=',LB,
     +    ' JOB=',JOB,' LWS=',LWS,' LIWS=',LIWS
        WRITE (MP,'(A,1P,3E11.3,0P,4F5.0,1P,/6X,4E11.3)') ' CNTL=',
     +    (CNTL(K),K=1,11)
      END IF
      IF (JOB.EQ.0) THEN
C Check LWS and LIWS
        L = 4*M + 4*N + 7 + 3
        IF (LWS.LT.L) THEN
          JOB = -13
          IF (LP.GE.0) WRITE (LP,FMT='(/,A,I7)')
     +        ' LA04 error return: LWS must be at least',L
          RINFO(35) = L
          RETURN
        END IF
        IF (LIWS.LT.LA+10) THEN
          JOB = -14
          IF (LP.GE.0) WRITE (LP,FMT='(/,A,I7)')
     +        ' LA04 error return: LIWS must be at least',LA+10
          RINFO(35) = LA
          RETURN
        END IF
C Increment NZ to allow for nonzeros in B and C
        DO 10 I = 1,M
          IF (B(I).NE.ZERO) NZ = NZ + 1
   10   CONTINUE
        DO 20 J = 1,N
          IF (C(J).NE.ZERO) NZ = NZ + 1
   20   CONTINUE
      ELSE
C Check IB
        IF (IB.LT.M) THEN
          JOB = -4
          IF (LP.GE.0) WRITE (LP,FMT='(/,A,I7)')
     +        ' LA04 error return: IB must be at least',M
          RINFO(35) = M
          RETURN
        END IF
      END IF
C Check size of LA
      IF (LA.LT.NZ) THEN
        JOB = -10
        IF (LP.GE.0) WRITE (LP,FMT='(/,A,I7)')
     +      ' LA04 error return: LA must be at least',NZ
        RINFO(35) = NZ
        RETURN
      END IF
C Check M, N, KB, LB
      IF (M.LT.0 .OR. KB.LT.0 .OR. LB.LE.KB .OR. N+1.LT.LB) THEN
        JOB = -3
        IF (LP.GE.0) WRITE (LP,FMT='(/,A,/5(A,I7))')
     +      ' LA04 error return: M, N, KB, or LB is faulty',' M=',M,
     +      ' N=',N,' KB=',KB,' LB=',LB
        RETURN
      END IF
C Check BND
      DO 30 J = 1,KB
        IF (BND(1,J).LT.-SIGMA .OR. BND(1,J).GT.BND(2,J) .OR.
     +      BND(2,J).GT.SIGMA) THEN
          JOB = -5
          RINFO(35) = J
          IF (LP.GE.0) WRITE (LP,FMT='(/A/A,I6,A,1P,2E20.12)')
     +        ' LA04 error return: BND is faulty',' BND(:,',J,') =',
     +        BND(1,J),BND(2,J)
          RETURN
        END IF
   30 CONTINUE
C IWS is used to check for duplicate entries
      DO 40 I = 1,M
        IWS(I) = 0
   40 CONTINUE
      DO 60 J = 1,N
        K1 = IP(J)
        K2 = IP(J+1) - 1
C   Check IP
        IF (K1.GT.K2+1) THEN
          JOB = -6
          RINFO(35) = J + 1
          IF (LP.GE.0) WRITE (LP,FMT='(/,A,/2(A,I6,A,I6))')
     +        ' LA04 error return: IP is faulty','  IP(',J,' )=',K1,
     +        '  IP(',J + 1,' )=',IP(J+1)
          RETURN
        END IF
        DO 50 K = K1,K2
C   Check IRN
          I = IRN(K)
          IF (I.LE.0 .OR. I.GT.M) THEN
            JOB = -7
            RINFO(35) = K
            IF (LP.GE.0) WRITE (LP,FMT='(/,A,/,A,I6,A,I6)')
     +          ' LA04 error return: IRN is faulty',' IRN(',K,' )=',I
            RETURN
C  Check for duplicate entries
          ELSE IF (IWS(I).EQ.J) THEN
            JOB = -8
            RINFO(35) = K
            IF (LP.GE.0) WRITE (LP,FMT='(/,A,I6,A,I6)')
     +          ' LA04 error return: two entries for row',I,
     +          ' and column',J
            RETURN
          END IF
          IWS(I) = J
   50   CONTINUE
   60 CONTINUE

      IF (JOB.EQ.0) THEN
C Perform scaling
        IF (M.GT.0 .AND. N.GT.0) CALL LA04SD(A,LA,IRN,IP,M,N,B,C,BND,KB,
     +                                SIGMA,WS(N+2),WS,WS(M+N+3),IWS)
        JOB = 1
        RETURN
      END IF

C Initialize and turn LA15's messages off during the use of LA04
      CALL LA15ID(IWS(II15),WS(IC15),IWS(IK15))
      IWS(II15+1) = 0

C Initialize x
      DO 70 I = 1,N + M
        X(I) = ZERO
   70 CONTINUE

C Initialize the random number generator (FA14)
      IF (CNTL(12).EQ.ZERO) CALL FA14ID(IWS(IISD))

      IF (JOB.EQ.2) THEN
C Basis or partial basis specified

C  Set weights G(I),I=1,N to unity.
        DO 80 J = 1,N
          G(J) = ONE
   80   CONTINUE

C Check IX and set G(J) for in-basis variables
        NPM = N + M
        DO 90 I = 1,M
          IF (IX(I).LT.0 .OR. IX(I).GT.3*NPM) THEN
            IF (LP.GE.0) WRITE (LP,'(/,A,/,A,I6,A,I11)')
     +          ' LA04 error return: IX is faulty',' IX(',I,') =',IX(I)
            RINFO(35) = I
            JOB = -11
            RETURN
          END IF
C Replace any infeasibilites in the basis by corresponding
C variables.
          IF (IX(I).GT.NPM) THEN
            IX(I) = IX(I) - NPM
            IF (IX(I).GT.NPM) IX(I) = IX(I) - NPM
          END IF
          J = IX(I)
          IF (J.GE.1 .AND. J.LE.N) THEN
            G(J) = -I
            IF (J.LE.KB) THEN
              JX(J) = -1
              IF (BND(1,J).EQ.BND(2,J)) THEN
C Force fixed variables out of basis
                G(J) = ZERO
                JX(J) = 1
                IX(I) = 0
              END IF
            END IF
          END IF
   90   CONTINUE

C Check JX
        DO 100 I = 1,KB
          IF (JX(I).EQ.-1) THEN
          ELSE IF (BND(1,I).EQ.BND(2,I)) THEN
            G(I) = ZERO
            JX(I) = 1
          ELSE IF (JX(I).EQ.1) THEN
            IF (BND(1,I).EQ.-SIGMA) GO TO 110
          ELSE IF (JX(I).EQ.2) THEN
            IF (BND(2,I).EQ.SIGMA) GO TO 110
          ELSE IF (JX(I).EQ.0) THEN
            IF (BND(1,I).GE.ZERO .OR. BND(2,I).LE.ZERO) GO TO 110
          ELSE
            GO TO 110
          END IF
  100   CONTINUE

  110   IF (I.LE.KB) THEN
          JOB = -12
          RINFO(35) = I
          IF (LP.GE.0) WRITE (LP,'(/,2A,I6,/,A,I6,A,1P,2E9.1)')
     +        ' LA04 error return: JX is faulty',' for variable I =',I,
     +        ' JX(I) = ',JX(I),' BND(:,I) = ',BND(1,I),BND(2,I)

          RETURN
        END IF
      END IF

C     Compute scale
      RINFO(23) = ZERO
      DO 120 I = 1,IP(N+1) - 1
        RINFO(23) = MAX(RINFO(23),ABS(A(I)))
  120 CONTINUE

C Call LA04B to find initial basis and weights G(I)
      IF (JOB.EQ.1) THEN
        CALL LA04BD(A,LA,IRN,IP,M,N,BND,KB,CNTL,IX,G,RINFO,IWS(IWB),
     +              IWS(IROW),WS(JV),WS(JY))
C Choose values for out-of-basis variables
        DO 130 J = 1,KB
          IF (G(J).GE.ZERO) THEN
            IF (BND(1,J).GE.ZERO) THEN
              JX(J) = 1
            ELSE IF (BND(2,J).LE.ZERO) THEN
              JX(J) = 2
            ELSE
              JX(J) = 0
            END IF
          ELSE
            JX(J) = -1
          END IF
  130   CONTINUE
      END IF

      JOB = 4

C Call subroutines LA04D (normal iteration) and LA04C (factorize basis
C    and compute solution and reduced costs).
C If JOB=3 then LA04D is called. this may reset JOB to 8 to
C     indicate that an LA04C call is needed.
C If JOB=4 then a call of LA04C is made. It is followed by a call of
C     LA04D, unless LA04C resets JOB to indicate termination or an
C     error condition.
C If JOB=5 then a call of LA04C is made. It is followed by a call of
C     LA04D, unless LA04C resets JOB to indicate termination or an
C     error condition. It may also be followed by calls of LA04C
C     and LA04D as for an JOB=4 entry.

  140 IF (JOB.NE.3) THEN
C     Factorize basis (unless JOB=5) and compute solution and reduced
C     costs (even if JOB=5).
        CALL LA04CD(A,LA,IRN,IP,M,N,B,C,BND,KB,LB,JOB,CNTL,IX,JX,G,
     +              RINFO,X,Z,WS(JB),IB,IWS,IWS(IPB),IWS(IWB),WS(JV),WS,
     +              WS(JY),WS(JYY),IWS(II15+1),WS(IC15+1),IWS(IISD))
      END IF
      IF (JOB.GE.3) THEN
C Perform simplex step
        CALL LA04DD(A,LA,IRN,IP,M,N,B,C,BND,KB,LB,JOB,CNTL,X,IX,JX,Z,G,
     +              RINFO,WS(JB),IWS,IWS(IPB),IWS(IWB),WS,WS(JV),
     +              IWS(IR),WS(JY),WS(JYY),IB,IWS(II15+1),WS(IC15+1))
      END IF
      IF (JOB.GE.5) THEN
C   LA04D has requested recomputation of solution and
C   reduced costs, possibly with refactorization of basis.
        CALL LA04CD(A,LA,IRN,IP,M,N,B,C,BND,KB,LB,JOB,CNTL,IX,JX,G,
     +              RINFO,X,Z,WS(JB),IB,IWS,IWS(IPB),IWS(IWB),WS(JV),WS,
     +              WS(JY),WS(JYY),IWS(II15+1),WS(IC15+1),IWS(IISD))
        IF (JOB.GE.5) THEN
C   LA04C has requested continued iteration.
          CALL LA04DD(A,LA,IRN,IP,M,N,B,C,BND,KB,LB,JOB,CNTL,X,IX,JX,Z,
     +                G,RINFO,WS(JB),IWS,IWS(IPB),IWS(IWB),WS,WS(JV),
     +                IWS(IR),WS(JY),WS(JYY),IB,IWS(II15+1),WS(IC15+1))
        END IF
      END IF

C Accumulate time spent in LA04A and the subroutines it calls.
C     print any error messages.
  150 call cpu_time(DUM)
      RINFO(4) = DUM - RINFO(9) + RINFO(4)
      IF (JOB.GT.0) THEN
        IF(JOB.EQ.8) JOB = 4
      ELSE IF (JOB.EQ.0) THEN
        IF (MP.GE.0) WRITE (MP,FMT='(/,A,I6,A,/,A,1P,E14.6)')
     +      ' LA04: Optimal solution found after',INT(RINFO(30)),
     +      ' basis factorizations',' Objective function = ',RINFO(1)
      ELSE IF (JOB.EQ.-2) THEN
        IF (RINFO(31).GT.ONE) THEN
          IF (LP.GE.0) WRITE (LP,FMT='(/,A)')
     +      ' LA04 error return: the objective function is unbounded'
        ELSE
          JOB = -1
          IF (LP.GE.0) WRITE (LP,FMT='(/,A)')
     +      ' LA04 error return: the feasible region is void'
        END IF
C Load unbounded direction into WS
        DO 160 I = M,1,-1
          WS(N+I) = WS(I)
  160   CONTINUE
        DO 170 J = 1,N
          WS(J) = ZERO
  170   CONTINUE
        JIN = RINFO(26)
        SIGN = ONE
        IF (Z(JIN).LT.ZERO) SIGN = -ONE
        WS(JIN) = -SIGN
        DO 180 I = 1,M
          J = IX(I)
          IF (J.LE.N) WS(J) = WS(N+I)*SIGN
  180   CONTINUE
      ELSE IF (JOB.EQ.-1) THEN
        IF (LP.GE.0) WRITE (LP,FMT='(/,A)')
     +      ' LA04 error return: the feasible region is void'
      ELSE IF (JOB.EQ.-4) THEN
        L = RINFO(35)
        IF (LP.GE.0) WRITE (LP,FMT='(/,A,I7)')
     +      ' LA04 error return: IB must be at least',L
      END IF

      END

      SUBROUTINE LA04ID(CNTL)
      DOUBLE PRECISION CNTL(15),ONE

      CNTL(1) = HUGE(ONE)
      CNTL(2) = EPSILON(ONE)** (2.0D0/3.0D0)
      CNTL(3) = 0
      CNTL(4) = 1
      CNTL(5) = 0
      CNTL(6) = 6
      CNTL(7) = 6
      CNTL(8) = 0.1D0
      CNTL(9) = 1.0D0/SQRT(EPSILON(ONE))
      CNTL(10) = EPSILON(ONE)
      CNTL(11) = CNTL(2)
      CNTL(12) = 0
      END

      SUBROUTINE LA04BD(A,LA,IRN,IP,M,N,BND,KB,CNTL,IX,G,RINFO,IW,IROW,
     +                  PV,Y)
C     This subroutine finds an initial basis and associated weights.
      INTEGER LA,M,N,KB
      DOUBLE PRECISION
     +                A(LA),BND(2,KB),CNTL(15),G(N),RINFO(40),PV(M),Y(M)
      INTEGER IP(N+1),IRN(LA),IROW(M),IW(M),IX(M)

C A,LA,...,N (IN) specify the matrix A.
C BND (IN) BND(1,J) and BND(2,J) hold the bounds for variable J,
C     J = 1,2,...,KB.
C KB  (IN) Second extent of BND.
C CNTL (IN) controls the action. The following components are accessed
C   CNTL(8) Relative pivot tolerance. Each diagonal entry of the
C     basis must have modulus at least CNTL(8) times the largest in its
C     column.
C   CNTL(7) Unit for monitor printing.
C   CNTL(4) 1 if steepest-edge weights are to be calculated.
C IX (OUT) is set to contain the indices of the variables in the initial
C     basis. The index I+N is used for the artifical variable
C     corresponding to column I of the identity matrix.
C G (OUT) G(J), J=1,2,...,N, is set to
C   For an in-basis variable, -(its position in the basis).
C   For an out-of-basis fixed-value, zero.
C   For an out-of-basis variable, one or the steepest-edge weight.
C RINFO(34) is set to the number of matrix sweeps needed.
C IW is a workspace array used to mark rows that have been used in
c     the construction of the basis. Row I has been pivotal if IW(I)>0.
C IROW is a workspace array. IROW(I) holds the row index of the
C     i-th pivot.
C PV is a workspace array. PV(I) holds the I-th pivot.
C Y is a workspace array used when when constructing the weights.
C
C     This code finds an upper triangular basis by first
C     looking for singletons and then doubletons, etc. excluding rows
C     already pivoted on.  Where necessary it fills in the basis with
C     artificial variables.

C     .. Local Scalars ..
      DOUBLE PRECISION
     +               AMAX,CMAX,ONE,XI,XX,ZERO
      PARAMETER (ONE=1,ZERO=0)
      INTEGER I,I2,ISWP,J,JJ,K,K1,K2,KK,KMAX,MINZ,MP,MZ,NART,NMINZ,NP,NZ

C AMAX   Modulus of the largest candidate pivot in the column.
C CMAX   Modulus of the largest element in the column.
C I      Temporary row index.
C I2     Temporary row index.
C ISWP   Sweep index.
C J      Temporary column index.
C JJ     Temporary column index.
C K      Index in A and IRN.
C K1     Column start.
C K2     Column end.
C KK     Temporary row index.
C KMAX   Index of the largest candidate pivot in the column.
C MINZ   Length of the shortest non-empty col.
C MP     Unit for monitor printing.
C MZ     Maximum acceptable length for a column.
C NART   Number of artificials in basis.
C NMINZ  Number of cols of shortest length.
C NP     Number of variables in basis.
C NZ     Number of nonzero entries in column.
C XI     Temporary real.
C XX     Temporary real.

C     .. Intrinsic Functions ..
      INTRINSIC ABS,MAX

      MP = CNTL(7)
      IF (MP.GE.0) WRITE (MP,FMT='(/,A/A)')
     +    ' LA04: Sweep      No.       No.     Min col   No. cols',
     +    '              variables artificials length   of min len'
      DO 10 J = 1,N
        G(J) = ONE
   10 CONTINUE
      DO 20 J = 1,KB
        IF (BND(1,J).EQ.BND(2,J)) G(J) = ZERO
   20 CONTINUE
      DO 30 I = 1,M
        Y(I) = ZERO
        IW(I) = 0
   30 CONTINUE
      NART = 0
      NP = 0
      MZ = 1
      DO 70 ISWP = 1,M
        MINZ = M + 1
        NMINZ = 0
        DO 60 J = 1,N
          IF (G(J).LE.ZERO) GO TO 60
          K1 = IP(J)
          K2 = IP(J+1) - 1
          IF (K2.LT.K1) GO TO 60
          AMAX = ZERO
          CMAX = ZERO
C Count nonzero entries in column, find the largest and find
C the largest candidate pivot.
          NZ = 0
          DO 40 K = K1,K2
            CMAX = MAX(CMAX,ABS(A(K)))
            KK = IRN(K)
            IF (IW(KK).GT.0) GO TO 40
            IF (A(K).NE.ZERO) NZ = NZ + 1
            IF (AMAX.GE.ABS(A(K))) GO TO 40
            AMAX = ABS(A(K))
            KMAX = K
   40     CONTINUE
          IF (AMAX.LT.CNTL(8)*CMAX) GO TO 60
          IF (NZ.EQ.0) GO TO 60
          IF (NZ.LE.MZ) THEN
C Place variable J in basis.
            I = IRN(KMAX)
            IW(I) = 1
            NP = NP + 1
            IX(NP) = J
            G(J) = -NP
            PV(NP) = A(KMAX)
            IROW(NP) = I
C Place artificial variables in basis.
            DO 50 K = K1,K2
              IF (A(K).EQ.ZERO) GO TO 50
              KK = IRN(K)
              IF (IW(KK).GT.0) GO TO 50
              IW(KK) = 1
              NART = NART + 1
              NP = NP + 1
              IROW(NP) = KK
              PV(NP) = ONE
              IX(NP) = N + KK
   50       CONTINUE
          END IF
          IF (NZ.LE.MINZ) THEN
            IF (NZ.LT.MINZ) THEN
              NMINZ = 1
              MINZ = NZ
            ELSE
              NMINZ = NMINZ + 1
            END IF
          END IF
   60   CONTINUE
        IF (MP.GE.0) WRITE (MP,FMT='(I10,4I10)') ISWP,NP,NART,MINZ,NMINZ
        IF (NP.EQ.M) GO TO 100
        IF (NMINZ.EQ.0) GO TO 80
        MZ = MINZ
C Increase MZ if it looks as if we are getting there too slowly.
        IF (MINZ*NMINZ* (M-ISWP).LE.M-NP) MZ = MZ + 1
   70 CONTINUE
      ISWP = ISWP - 1
C Make all the remaining variables artificial.
   80 DO 90 I = 1,M
        IF (IW(I).GT.0) GO TO 90
        NP = NP + 1
        PV(NP) = ONE
        IROW(NP) = I
        IX(NP) = N + I
   90 CONTINUE
  100 RINFO(34) = ISWP

C     Find initial weights
      IF (CNTL(4).EQ.ONE) THEN
        DO 150 J = 1,N
          IF (G(J).LE.ZERO) GO TO 150
          XX = ONE
          DO 110 K = IP(J),IP(J+1) - 1
            KK = IRN(K)
            Y(KK) = A(K)
  110     CONTINUE
          DO 140 I2 = M,1,-1
            I = IROW(I2)
            IF (Y(I).EQ.ZERO) GO TO 140
            XI = Y(I)/PV(I2)
            XX = XX + XI*XI
            JJ = IX(I2)
            IF (JJ.GT.N) GO TO 130
            DO 120 K = IP(JJ),IP(JJ+1) - 1
              KK = IRN(K)
              Y(KK) = Y(KK) - XI*A(K)
  120       CONTINUE
  130       Y(I) = ZERO
  140     CONTINUE
          G(J) = XX
  150   CONTINUE
      END IF
      END


      SUBROUTINE LA04CD(A,LA,IRN,IP,M,N,B,C,BND,KB,LB,JOB,CNTL,IX,JX,G,
     +                  RINFO,X,Z,BA,IB,INDB,IPB,IWB,V,W,Y,YY,I15,C15,
     +                  ISEED)
C  Compute solution, reduced costs and next incoming variable.
C  If requested, factorize basis first.

C     .. Scalar Arguments ..
      INTEGER LA,IB,JOB,KB,LB,M,N,ISEED
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION
     +                A(LA),B(M),BA(IB),BND(2,KB),C(N),CNTL(15),G(N),
     +                RINFO(40),V(M),W(M),X(N+M),Y(M),YY(M),Z(N),C15(3)
      INTEGER INDB(IB,2),IP(N+1),IPB(M,2),IRN(LA),IWB(M,8),IX(M),JX(KB),
     +        I15(10)

C A,LA,...,LB specify the problem and are left unchanged.
C JOB (INOUT) Possible values on entry:
C     4 Factorize basis and compute solution, reduced costs and next
C       incoming variable
C     5 Compute solution, reduced costs and next incoming variable.
C  Possible values on exit:
C     4,5 (unchanged) more iterations needed
C     0 solution found
C    -1 feasible region is empty
C    -4 IB is too small
C CNTL controls the action
C CNTL(1) (IN) Value that represents an infinite bound.
C CNTL(2) (IN) Relative tolerance for the entries of A, B, C.
C CNTL(3) (OUT) holds a tolerance for the bounds on the variables.
C       It may be increased on a recalculation of the solution.
C CNTL(5) (INOUT) < 1 No iterative refinement
C                 >=1 Perform iterative refinement
C CNTL(7) (IN) Unit for monitor output.
C CNTL(8) (IN) Relative pivot tolerance for LA15.
C CNTL(9) (IN) Limit on LA15 growth parameter.
C CNTL(10) (IN) The parameter CNTL(1) of LA15 is given the value
C       CNTL(10)*max|(A(K)|.
C IX (INOUT) IX(I), I=1,2,...,M, holds the index of column I of the
C   basis. Possible values:
C   IX(I)=0  Missing column (only partial basis specified).
C   IX(I)<=N  Ordinary variable
C   N<IX(I)<=N+M  Artifical variable with lower and upper bounds of zero
C   N+M<IX(I)<=2*(N+M) Variable IX(I)-(N+M) is above its upper bound.
C   2*(N+M)<IX(I) Variable IX(I)-2*(N+M) is below its lower bound.
C JX (INOUT) JX(J), J=1,2,...,KB, holds
C    -1 if variable J is in basis and feasible with value X(J)
C     0 if variable J is out of basis with value X(J)
C     1 if lower bound is active: variable has value X(J)+BND(1,J);
C       it may be in the basis and infeasible or out of the basis.
C     2 if upper bound is active: variable has value X(J)+BND(2,J);
C       it may be in the basis and infeasible or out of the basis.
C G (INOUT) must on input contain the steepest-edge weights of
C     out-of-basis variables, or unity if CNTL(4)=0. On output,
C     the components corresponding to in-basis variables will
C     have been reset to the negation of their position in the basis.
C RINFO(1) (OUT) holds the objective function or the sum of the modulii
C      of the infeasibilities.
C RINFO(2) (INOUT) holds the parameter U of LA15A.
C RINFO(3) (INOUT) holds the parameter G of LA15A.
C RINFO(4) (IN) holds the c.p.u. time spent in LA04A and subroutines
C      called from it since the last initial (JOB=1,2) entry.
C RINFO(5) (INOUT) holds the c.p.u. time spent in LA15A
C      since the last initial entry.
C RINFO(6) (INOUT) holds the c.p.u. time spent in LA15B
C      since the last initial entry.
C RINFO(7) holds the c.p.u. time spent in LA15C
C      since the last initial entry.
C RINFO(8) (INOUT) contains the c.p.u. time spent on the last
C     refactorization of the basis matrix and recalculation of X and Z.
C RINFO(9) (IN) holds the c.p.u. reading on the last LA04A call.
C RINFO(10) (OUT) holds the c.p.u. time spent in LA04A and subroutines
C      it calls between the last JOB=1,2 entry and the completion
C      of the last refactorization.
C RINFO(23) (IN) holds the modulus of the largest matrix entry.
C RINFO(24) (INOUT) holds the maximum permitted value for the parameter
C      G of LA15A. If this value is exceeded after an LA15C call, a
C      fresh factorization is made.
C RINFO(26) (OUT) holds the index of the variable entering the basis.
C RINFO(27) (OUT) holds the total number of infeasibilities (i.e. the
C      total number of artificial variables with nonzero values).
C RINFO(28) (INOUT) holds the number of iterations since the last
C       refactorization.
C RINFO(30) (INOUT) holds the number of refactorizations since the last
C       initial entry.
C RINFO(31) (INOUT) has the value 2 if the current solution is feasible
C       and 1 otherwise.
C RINFO(35) (OUT) set to suggested size for IB if JOB = -4.
C X (INOUT) holds the current solution. The basic part is recalculated.
C       The non-basic part is unchanged.
C Z (OUT) is set to the current reduced costs.
C BA (INOUT) passed to LA15 and used to contain the reals of the
C     factorized basis.
C IB (IN) Size of BA.
C INDB,...,IWB (INOUT) are passed to LA15 and used to contain
C     information about the factorized basis.
C V, W, Y and YY are workspace arrays.
C I15(J),J=1,3 used for LA15's ICNTL argument
C I15(J),J=4,10 used for LA15's KEEP argument
C C15(J),J=1,3 used for LA15's CNTL argument
C ISEED holds the random number generator seed (FA14)

C     .. Local Scalars ..
      DOUBLE PRECISION
     +                AW,BL,BU,DV,DW,DZI,E,F,G5,ONE,SIGMA,T,
     +                TIME,XJ,ZGMAX,ZI,ZERO
      PARAMETER (ONE=1,ZERO=0)
      DOUBLE PRECISION DUM
      LOGICAL FEAS
      INTEGER I,IS,IU,J,JIN,K,K1,K15,K2,L,MP,NINF,NPM

C AW    Temporary
C BL    Lower bound
C BU    Upper bound
C DUM   Used in cpu_time call.
C DV    Infinity norm of the last iterative refinement change to V
C DW    Infinity norm of the last iterative refinement change to W
C DZI   Bound for error in reduced cost.
C E     Value of the current variable.
C F     Objective function value.
C FEAS  True if the solution is feasible.
C G5    Argument G of LA15A
C I     Row index
C IS    Start index
C IU    Index for successive attempts to call LA15 with
C       increasing U values.
C J     Column index
C JIN   Index of incoming column
C K     Index in A and IRN.
C K1    Column start.
C K15   set to 3 where I15(K15+J),J=1,7 used for LA15's KEEP argument
C K2    Column end.
C L     Index in BA and INDB.
C MP    Unit for monitor output
C NINF  No. of infeasibilities.
C NPM   N+M
C SIGMA Value that represents an infinite bound.
C T     Temporary for holding the c.p.u. time
C TIME  c.p.u. time before commencing factorization
C XJ    Temporary for value of out-of-basis variable
C ZGMAX Max. squared weighted reduced cost of advantageous column.
C ZI    Reduced cost.

C     .. External Functions ..
      DOUBLE PRECISION FA14AD
      EXTERNAL FA14AD

C     .. External Subroutines ..
      EXTERNAL LA15AD,LA15BD,LA04ED

C     .. Intrinsic Functions ..
      INTRINSIC ABS,MIN,MAX

      K15 = 3

      MP = CNTL(7)
      SIGMA = CNTL(1)
      NPM = N + M
C Replace any infeasibilites in the basis by corresponding
C variables. Set basic variables to zero.
      DO 10 I = 1,M
        J = IX(I)
        IF (J.GT.NPM) THEN
          J = J - NPM
          IF (J.GT.NPM) J = J - NPM
          IX(I) = J
          IF (J.LE.KB) JX(J) = -1
        END IF
        IF (J.GT.0) X(J) = ZERO
   10 CONTINUE

      IF (JOB.EQ.5) GO TO 110

C     Factorize basis
      call cpu_time(TIME)
      RINFO(28) = 0
      RINFO(2) = CNTL(8)
      IF (M.EQ.0) THEN
C        initialize what was LA15's LENL, LENU and NCP
        I15(K15+1) = 0
        I15(K15+2) = 0
        I15(K15+3) = 0
        RINFO(3) = ZERO
      ELSE
        DO 90 IU = 1,10
          L = 1
          DO 30 I = 1,M
            J = IX(I)
            IF (J.LE.0) THEN
C  Partial basis specified - no column to load
            ELSE IF (J.LE.N) THEN
C  Load ordinary column
              G(J) = -I
              K1 = IP(J)
              K2 = IP(J+1) - 1
              IF (L+K2-K1.LT.IB) THEN
                DO 20 K = K1,K2
                  BA(L) = A(K)
                  INDB(L,1) = IRN(K)
                  INDB(L,2) = I
                  L = L + 1
   20           CONTINUE
              ELSE
                L = L + K2 - K1 + 1
              END IF
            ELSE
C  Load artificial column
              IF (L.LE.IB) THEN
                BA(L) = ONE
                INDB(L,1) = J - N
                INDB(L,2) = I
              END IF
              L = L + 1
            END IF
   30     CONTINUE
          IF (L.GT.IB) THEN
            JOB = -4
            RINFO(35) = L
            RETURN
          END IF
          call cpu_time(T)
          C15(1) = CNTL(10)*RINFO(23)
          CALL LA15AD(BA,INDB,L-1,IB,M,IPB,IWB,Y,G5,RINFO(2),I15,
     +                C15,I15(K15+1))
          RINFO(3) = G5
          RINFO(32) = I15(K15+1) + I15(K15+2)
          RINFO(33) = I15(K15+3)
          call cpu_time(DUM)
          RINFO(5) = RINFO(5) + DUM - T
          RINFO(30) = RINFO(30) + 1
          IF (G5.LT.-6.5) THEN
            JOB = -4
            RINFO(35) = NINT(IB*1.5)
            RETURN
          ELSE IF (G5.LT.-4.5) THEN
C     Construct new basis using LA15 output
            IF (MP.GE.0) WRITE (MP,FMT='(A)')
     +          ' New basis constructed by LA04C'
            I = 0
C I is set to a row of the basis that has not been pivotal
            K = 0
C K is set to a column of the basis that has not been pivotal
            DO 80 L = 1,M
              IS = K + 1
              DO 40 K = IS,M
                IF (IWB(K,6).GE.0) GO TO 50
   40         CONTINUE
              GO TO 90
   50         J = IX(K)
C J is the corresponding column of A
              IS = I + 1
              DO 60 I = IS,M
                IF (IWB(I,5).GE.0) GO TO 70
   60         CONTINUE
              GO TO 90
   70         IF (J.GE.1 .AND. J.LE.N) THEN
                G(J) = ONE
                IF (J.LE.KB) THEN
                  IF (BND(1,J).GE.ZERO) THEN
                    JX(J) = 1
                  ELSE IF (BND(2,J).LE.ZERO) THEN
                    JX(J) = 2
                  ELSE
                    JX(J) = 0
                  END IF
                END IF
              END IF
              IX(K) = I + N
   80       CONTINUE
          ELSE
            IF (G5.LT.CNTL(9)*RINFO(23) .OR.
     +          RINFO(2).EQ.ONE) GO TO 100
            RINFO(2) = MIN(ONE,RINFO(2)+RINFO(2))
          END IF
   90   CONTINUE
  100   RINFO(24) = MAX(G5*10.,RINFO(23)*CNTL(9)**1.5)
        IF (MP.GE.0) WRITE (MP,'(/,A,F9.4,/,A,3F11.3,/,A,I9)')
     +     ' LA04: Factorization of basis. Value of U used by LA15A was'
     +      ,RINFO(2),' Total times so far for LA15A/B/C are',RINFO(5),
     +      RINFO(6),RINFO(7),' Size of LA15 arrays is',IB
      END IF

C    Compute RHS for computation of X and an error estimate for X
  110 DO 120 I = 1,M
        W(I) = B(I)*FA14AD(ISEED,-1)*CNTL(2)
        V(I) = B(I)
  120 CONTINUE
      DO 140 J = 1,N
        XJ = X(J)
        IF (J.LE.KB) THEN
          IF (JX(J).EQ.1) THEN
            XJ = XJ + BND(1,J)
          ELSE IF (JX(J).EQ.2) THEN
            XJ = XJ + BND(2,J)
          END IF
        END IF
        IF (XJ.NE.ZERO) THEN
          DO 130 K = IP(J),IP(J+1) - 1
            I = IRN(K)
            W(I) = W(I) + A(K)*XJ*FA14AD(ISEED,-1)*CNTL(2)
            V(I) = V(I) - A(K)*XJ
  130     CONTINUE
        END IF
  140 CONTINUE
      DO 150 J = 1,M
        V(J) = V(J) - X(N+J)
  150 CONTINUE
      IF (M.GT.0) THEN
        call cpu_time(T)
        CALL LA15BD(BA,INDB,IB,M,IPB,IWB,Y,RINFO(3),W,.FALSE.,I15,
     +              I15(K15+1))
        call cpu_time(DUM)
        RINFO(6) = RINFO(6) + DUM - T
      END IF
C
C Calculate tolerance on components of X.
      DO 160 I = 1,M
        CNTL(3) = MAX(CNTL(3),ABS(W(I)))
  160 CONTINUE
      F = SIGMA
      DO 170 I = 1,KB
        IF (BND(1,I).NE.ZERO) F = MIN(ABS(BND(1,I)),F)
        IF (BND(2,I).NE.ZERO) F = MIN(ABS(BND(2,I)),F)
  170 CONTINUE
      IF (F.NE.SIGMA) CNTL(3) = MAX(F*CNTL(2),CNTL(3))

C Compute X and iteratively refine it
      DO 180 I = 1,M
        W(I) = V(I)
  180 CONTINUE
      IF (M.GT.0) THEN
        CALL LA15BD(BA,INDB,IB,M,IPB,IWB,Y,RINFO(3),V,.FALSE.,I15,
     +              I15(K15+1))
        CALL LA04ED(A,LA,IRN,IP,M,N,W,V,DV,IX,0,RINFO,BA,INDB,IPB,IWB,
     +              YY,Y,IB,.FALSE.,I15,C15)
        IF (MP.GE.0) WRITE (MP,'(A,1P,E12.4)')
     +      ' Largest change, iterative refinement of X:',DV
        CNTL(3) = MAX(DV,CNTL(3))
      END IF
      IF (MP.GE.0) WRITE (MP,'(/,A,1P,E9.1)')
     +    ' LA04: tolerance on the bounds is',CNTL(3)

C Construct infeasibility objective vector W and check if feasible
      NINF = 0
      F = ZERO
      DO 190 I = 1,M
        E = V(I)
        W(I) = ZERO
        J = IX(I)
        X(J) = E
        IF (J.GT.N) THEN
C           Artificial variable in basis
          IF (E.LT.-CNTL(3)) THEN
            W(I) = -ONE
            IX(I) = J + NPM
            NINF = NINF + 1
            F = F - E
          ELSE IF (E.GT.CNTL(3)) THEN
            W(I) = ONE
            IX(I) = J + NPM + NPM
            NINF = NINF + 1
            F = F + E
          END IF
        ELSE IF (J.LE.KB) THEN
C       Variable with two bounds
          BL = BND(1,J)
          BU = BND(2,J)
          IF (BL.NE.-SIGMA) THEN
            IF (E.LT.BL-CNTL(3)) THEN
C           Variable below lower bound
              E = E - BL
              X(J) = E
              W(I) = -ONE
              IX(I) = J + NPM
              NINF = NINF + 1
              F = F - E
              JX(J) = 1
              GO TO 190
            ELSE IF (E.LT.BL) THEN
              GO TO 190
            END IF
          END IF
          IF (BU.NE.SIGMA) THEN
            IF (E.GT.BU+CNTL(3)) THEN
C           Variable above upper bound
              W(I) = ONE
              E = E - BU
              X(J) = E
              IX(I) = J + NPM + NPM
              NINF = NINF + 1
              F = F + E
              JX(J) = 2
            END IF
          END IF
        ELSE IF (J.GE.LB) THEN
C         Variable bounded below by zero
          IF (E.LT.-CNTL(3)) THEN
C             Variable below bound of zero
            IX(I) = J + NPM
            W(I) = -ONE
            NINF = NINF + 1
            F = F - E
          END IF
        END IF
  190 CONTINUE

C If no. of infeasibilites increases, switch to iterative refinement
      IF (NINF.GT.RINFO(27)) THEN
        IF (CNTL(5).LT.1) THEN
          CNTL(5) = 1
          IF (MP.GE.0) WRITE (MP,'(/,A)')
     +        ' LA04: starting iterative refinement'
        END IF
      END IF

      FEAS = NINF .EQ. 0
      IF (FEAS) THEN
C Feasible, so set up the cost vector and evaluate objective function
        RINFO(31) = 2
        DO 200 J = 1,N
          F = F + C(J)*X(J)
          IF (G(J).LT.ZERO) THEN
            Z(J) = ZERO
            I = -G(J)
            W(I) = C(J)
          ELSE
            Z(J) = C(J)
          END IF
  200   CONTINUE
        DO 210 J = 1,KB
          IF (JX(J).EQ.1) THEN
            F = F + C(J)*BND(1,J)
          ELSE IF (JX(J).EQ.2) THEN
            F = F + C(J)*BND(2,J)
          END IF
  210   CONTINUE
      ELSE
        RINFO(31) = 1
        DO 220 J = 1,N
          Z(J) = ZERO
  220   CONTINUE
      END IF
      RINFO(1) = F
      RINFO(27) = NINF
      IF (MP.GE.0) THEN
        WRITE (MP,'(/,A,I5)') ' LA04: number of infeasibilites:',NINF
        WRITE (MP,'(A,1P,E14.6)') ' Objective function:',F
      END IF

C     Compute  W = B**(-T)*c, with iterative refinement
      IF (M.GT.0) THEN
        DO 230 I = 1,M
          V(I) = W(I)
  230   CONTINUE
        call cpu_time(T)
        CALL LA15BD(BA,INDB,IB,M,IPB,IWB,Y,RINFO(3),W,.TRUE.,I15,
     +              I15(K15+1))
        call cpu_time(DUM)
        RINFO(6) = RINFO(6) + DUM - T
        CALL LA04ED(A,LA,IRN,IP,M,N,V,W,DW,IX,0,RINFO,BA,INDB,IPB,IWB,
     +              YY,Y,IB,.TRUE.,I15,C15)
        IF (MP.GE.0) WRITE (MP,'(A,1P,E12.4)')
     +      ' Largest change, iterative refinement of reduced costs:',DW
C Replace very small elements by zero
        DO 240 I = 1,M
          IF (ABS(W(I)).LE.DW) W(I) = ZERO
  240   CONTINUE
      END IF

C     Compute initial reduced costs and determine first pivot col.
      ZGMAX = ZERO
      JIN = 0
      DO 260 I = 1,N
        IF (G(I).LE.ZERO) THEN
C             Variable in the basis or at fixed value
          Z(I) = ZERO
          GO TO 260
        END IF
        ZI = Z(I)
        DZI = ABS(ZI)
        DO 250 K = IP(I),IP(I+1) - 1
          AW = A(K)*W(IRN(K))
          DZI = DZI + ABS(AW)
          ZI = ZI - AW
  250   CONTINUE
        DZI = DZI*CNTL(2)
        IF (I.GE.LB) THEN
C         Variable has lower bound 0
          ZI = ZI + DZI
          Z(I) = ZI
          IF (ZI.GE.ZERO) GO TO 260
        ELSE IF (I.GT.KB) THEN
C        Free variable
          Z(I) = ZI
          IF (ABS(ZI).LE.DZI) THEN
            Z(I) = ZERO
            GO TO 260
          END IF
        ELSE IF (JX(I).EQ.0) THEN
C        Free variable
          Z(I) = ZI
          IF (ABS(ZI).LE.DZI) THEN
            Z(I) = ZERO
            GO TO 260
          END IF
        ELSE IF (JX(I).EQ.1) THEN
C   Lower bound is active
          ZI = ZI + DZI
          Z(I) = ZI
          IF (ZI.GE.ZERO) GO TO 260
        ELSE IF (JX(I).EQ.2) THEN
C   Upper bound is active
          Z(I) = ZI - DZI
          ZI = -Z(I)
          IF (ZI.GE.ZERO) GO TO 260
        END IF
C Column is a candidate
        ZI = ZI*ZI/G(I)
        IF (ZGMAX.LT.ZI) THEN
          ZGMAX = ZI
          JIN = I
        END IF
  260 CONTINUE
      RINFO(26) = JIN

      call cpu_time(DUM)
      RINFO(10) = DUM
      IF (JOB.NE.5) RINFO(8) = RINFO(10) - TIME
      RINFO(10) = RINFO(10) - RINFO(9) + RINFO(4)
C
      IF (ZGMAX.EQ.ZERO) THEN
        IF (FEAS) THEN
          JOB = 0
        ELSE
          JOB = -1
        END IF
      END IF

      END


      SUBROUTINE LA04DD(A,LA,IRN,IP,M,N,B,C,BND,KB,LB,JOB,CNTL,X,IX,JX,
     +                Z,G,RINFO,BA,INDB,IPB,IWB,W,V,IR,Y,YY,IB,I15,C15)
C Perform a simplex iteration, starting with a given incoming column
C     and ending with choosing the next column.

C     .. Scalar Arguments ..
      INTEGER LA,IB,JOB,KB,LB,M,N

C     .. Array Arguments ..
      DOUBLE PRECISION
     +                A(LA),BA(IB),BND(2,KB),B(M),C(N),CNTL(15),G(N),
     +                RINFO(40),V(M),W(M),X(N+M),Y(M),YY(M),Z(N),C15(3)
      INTEGER INDB(IB,2),IP(N+1),IPB(M,2),IR(M),IRN(LA),IWB(M,4),IX(M),
     +        JX(KB),I15(10)

C A,LA,...,LB (IN) specify the problem.
C JOB (INOUT) Possible values on entry:
C     3 normal iteration
C     4 iteration after factorizing and computing X, Z
C     5 iteration after recomputing X, Z without refactorization
C    Possible values on exit:
C     0 failed to find outgoing column in phase 2
C    -1 feasible region is empty
C    -2 objective function is unbounded
C    -4 IB too small
C     3 normal iteration
C     4 normal iteration completed but refactorization needed
C     5 no pivot found
C     8 iteration aborted because of LA15 problem.
C CNTL(1) (IN) Value that represents an infinite bound.
C CNTL(3) (IN) holds a tolerance for the bounds on the variables.
C       It may be increased on a recalculation of the solution.
C CNTL(5) (IN) < 1 No iterative refinement
C              >=1 Perform iterative refinement
C CNTL(7) (IN) Unit for monitor output.
C CNTL(11) In the vector w of changes to x, entries less than
C       CNTL(10)*norm(w), for the 2-norm, are treated as zero.
C X (INOUT) holds the current solution.
C IX (INOUT) holds the basis.
C JX (INOUT) JX(J), J=1,2,...,KB, holds
C    -1 if variable J is in basis and feasible with value X(J)
C     0 if variable J is out of basis with value X(J)
C     1 if lower bound is active: variable has value X(J)+BND(1,J);
C       it may be in the basis and infeasible or out of the basis.
C     2 if upper bound is active: variable has value X(J)+BND(2,J);
C       it may be in the basis and infeasible or out of the basis.
C Z (INOUT) holds the reduced costs.
C G (INOUT) holds the steepest-edge weights of out-of-basis
C     variables, or unity if CNTL(4)=0.
C     For in-basis variables, it holds the negation of their
C     position in the basis.
C RINFO (INOUT) holds the real information.
C RINFO(1) (INOUT) holds the objective function or the sum of the moduli
C      of the infeasibilities.
C RINFO(2) (INOUT) holds the parameter U of LA15A/C.
C RINFO(3) (IN) contains the parameter G of LA15B/C.
C RINFO(4) (IN) holds the c.p.u. time spent in LA04A and subroutines
C      called from it since the last initial (JOB=1,2) entry.
C RINFO(6) (INOUT) holds the c.p.u. time spent in LA15B
C      since the last initial entry.
C RINFO(7) (INOUT) holds the c.p.u. time spent in LA15C
C      since the last initial entry.
C RINFO(8) (IN) contains the c.p.u. time spent on the last
C      refactorization of the basis matrix and recalculation of X and Z.
C RINFO(9) (IN) holds the c.p.u. reading on the last LA04A call.
C RINFO(10) (IN)  holds the c.p.u. time spent in LA04A and subroutines
C      it calls between the last JOB=1,2 entry and the completion
C      of the last refactorization.
C RINFO(11) to RINFO(20) (one is OUT, rest are IN) are used to hold the
C      last ten values of the average iteration time since the last
C      refactorization.
C RINFO(21) (INOUT) holds the minimum average iteration time since the
C      last refactorization.
C RINFO(24) (IN) holds the maximum permitted value for the parameter G
C      of LA15A. If this value is exceeded after an LA15C call, a
C      fresh factorization is made.
C RINFO(25) (OUT) holds the index of the variable leaving the basis.
C RINFO(26) (INOUT) holds the index of the variable entering the basis.
C RINFO(27) (INOUT) holds the total number of infeasibilities (i.e. the
C      total number of variables outside their bounds by more
C      than CNTL(3)).
C RINFO(28) (INOUT) holds the number of iterations since the last
C       refactorization.
C RINFO(29) (INOUT) holds the number of iterations since the last
C       initial (JOB=1,2) entry.
C RINFO(31) (IN) has the value 2 if the current solution is feasible
C       and 1 otherwise.
C RINFO(32) (OUT) holds the number of entries in the factorized basis
C RINFO(33) (OUT) holds the no. of compresses of the factorized basis
C       since the previous refactorization.
C RINFO(35) (OUT) holds information on a return with JOB <= -4.
C BA (INOUT) passed to LA15 and used to contain the reals of the
C     factorized basis.
C IB (IN) Size of BA.
C INDB,...,IWB (INOUT) are passed to LA15 and used to contain
C     information about the factorized basis.
C W is a workspace array. On a return with ICNTL=-2, it holds the
C    basic part of the unbounded direction
C V is a workspace array for LA15.
C IR is a workspace array.
C Y is a workspace array.
C YY is a workspace array.
C IB (IN) holds the size of array BA.
C I15(J),J=1,3 used for LA15's ICNTL argument
C I15(J),J=4,10 used for LA15's KEEP argument
C C15(J),J=1,3 used for LA15's CNTL argument

C     .. Local Scalars ..
      DOUBLE PRECISION
     +                ALPHA,BETA,DV,DW,ERR,GQ,G5,PIVOT,PVMIN,R,RMIN,RP,
     +                RPMIN,SIGMA,SIGN,SIGNX,T,WW,XJ,ZGMAX,ZI,ZQ,
     +                ONE,ZERO
      PARAMETER (ONE=1,ZERO=0)
      DOUBLE PRECISION DUM
      LOGICAL FEAS
      INTEGER I,IK,IRP,ITER,J,JIN,JOUT,JP,K,K1,K15,K2,KK,MP,NINF,NPM
C ALPHA Temporary used when revising the column weights
C BETA  Temporary used when revising the column weights
C DV    Infinity norm of the last iterative refinement change to V
C DW    Infinity norm of the last iterative refinement change to W
C ERR   Estimate of the randomness of the average iteration time.
C GQ    Weight gamma of incoming column.
C G5    Argument G of LA15A
C I     Row index.
C IK    Number of iterations since the last refactorization.
C IRP   holds the number of columns that a near to a tie for being the
C       outgoing column.
C ITER  number of iterations since the last initial entry.
C J     Column index
C JIN   Index of incoming column
C JOUT  Index of outgoing column
C JP    Index in basis of changing column.
C K     Index in A and IRN.
C K1    Column start.
C K15   set to 3 where I15(K15+J),J=1,7 used for LA15's KEEP argument
C K2    Column end.
C KK    Row index.
C MP    Unit for monitor output
C NINF  No. of infeasibilites.
C NPM   N+M
C PIVOT The pivot.
C PVMIN Limit on size of acceptable pivot.
C R     Ratio used when choosing pivot row.
C RMIN  Min. R value among candidates.
C RP    Relaxed value of R.
C RPMIN Min. RP value among candidates..
C SIGMA Value that represents an infinite bound.
C SIGN   -(sign of Z(JIN))
C SIGNX For an artificial variable leaving the basis, the
C       corresponding component of the infeasibity objective vector.
C T     Temporary for holding the c.p.u. time
C WW    Candidate pivot.
C XJ    Temporary for component of X.
C ZGMAX Max. squared weighted reduced cost of advantageous column.
C ZI    Reduced cost.
C ZQ    Reduced cost of the incoming column
C     ..
C     .. External Subroutines ..
      EXTERNAL LA04ED,LA15BD,LA15CD
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MIN,SQRT,MAX,MOD
C     ..
      K15 = 3

      MP = CNTL(7)
      SIGMA = CNTL(1)
      JIN = RINFO(26)
      NINF = RINFO(27)
      IK = RINFO(28)
      ITER = RINFO(29)
      FEAS = RINFO(31) .EQ. 2
      JOUT = 0
      NPM = N + M

      IF (JIN.EQ.0) THEN
        IF (MP.GE.0) WRITE (MP,FMT='(/,/,A)')
     +      ' LA04: no candidates are available for entry into basis'
        JOB = 5
        GO TO 280
      END IF
C
C     Start of iteration itself
C Compute W = B**(-1)*A(JIN). Store intermediate L**(-1)*A(JIN) in Y.
      SIGN = ONE
      IF (Z(JIN).GT.ZERO) SIGN = -ONE
      DO 10 I = 1,M
        W(I) = ZERO
   10 CONTINUE
      DO 20 K = IP(JIN),IP(JIN+1) - 1
        KK = IRN(K)
        W(KK) = A(K)
   20 CONTINUE
      IF (M.GT.0) THEN
        call cpu_time(T)
        CALL LA15BD(BA,INDB,IB,M,IPB,IWB,Y,RINFO(3),W,.FALSE.,I15,
     +              I15(K15+1))
        call cpu_time(DUM)
        RINFO(6) = RINFO(6) + DUM - T

        IF (CNTL(5).GE.1) THEN
          CALL LA04ED(A,LA,IRN,IP,M,N,B,W,DW,IX,JIN,RINFO,BA,INDB,IPB,
     +                IWB,YY,V,IB,.FALSE.,I15,C15)
        END IF
      ELSE
        DW = ZERO
      END IF

C     Compute Z(JIN) afresh
      ZQ = ZERO
      IF (.NOT.FEAS) THEN
        DO 30 I = 1,M
          IF (IX(I).GT.NPM) THEN
            IF (IX(I).GT.NPM+NPM) THEN
              ZQ = ZQ - W(I)
            ELSE
              ZQ = ZQ + W(I)
            END IF
          END IF
   30   CONTINUE
      ELSE
        ZQ = C(JIN)
        DO 40 I = 1,M
          J = IX(I)
          IF (J.LE.N) THEN
            ZQ = ZQ - C(J)*W(I)
          END IF
   40   CONTINUE
      END IF

      ZI = Z(JIN)
      IF (ABS(ZQ-ZI).GE.MAX(ABS(ZQ),ABS(ZI))*0.1) THEN
        IF (MP.GE.0) WRITE (MP,FMT='(/,/,A,1P,2E12.4)')
     +      ' Approximate and newly computed reduced costs are',Z(JIN),
     +      ZQ
        IF (JOB.EQ.3) THEN
          JOB = 5
        ELSE
          JOB = -1
          IF (FEAS) JOB = 0
        END IF
        GO TO 280
      END IF

Compute G(IQ) afresh.
      GQ = ONE
      DO 50 I = 1,M
        GQ = GQ + W(I)**2
   50 CONTINUE

C     First scan using row selection scheme of P.Harris. Check for
C     unbounded solution.
      RPMIN = SIGMA
      RMIN = SIGMA
      PVMIN = CNTL(11)*SQRT(GQ)
      IF (CNTL(5).GE.1) PVMIN = MAX(DW,PVMIN)
      IRP = 0
      DO 60 I = 1,M
        WW = W(I)*SIGN
        J = IX(I)
        IF (J.GT.NPM) THEN
          IF (J.GT.NPM+NPM) THEN
            XJ = X(J-NPM-NPM)
          ELSE
            XJ = X(J-NPM)
          END IF
        ELSE
          XJ = X(J)
        END IF
        RP = ZERO
        IF (WW.GT.ZERO) THEN
          IF (J.LE.KB) THEN
            IF (BND(1,J).NE.-SIGMA) THEN
C                 Explicit lower bound may become active
              R = (XJ-BND(1,J))/WW
              RP = R + CNTL(3)/WW
            END IF
          ELSE IF (J.GT.NPM+NPM) THEN
C              Artificial variable bounded below
            R = XJ/WW
            RP = R + CNTL(3)/WW
          ELSE IF (J.GE.LB) THEN
            R = XJ/WW
            RP = R + CNTL(3)/WW
          END IF
        ELSE IF (WW.LT.ZERO) THEN
          IF (J.LE.KB) THEN
            IF (BND(2,J).NE.SIGMA) THEN
C                 Explicit upper bound may become active
              R = (XJ-BND(2,J))/WW
              RP = R - CNTL(3)/WW
            END IF
          ELSE IF (J.GT.N .AND. J.LE.NPM+NPM) THEN
C              Artificial variable bounded above
            R = XJ/WW
            RP = R - CNTL(3)/WW
          END IF
        END IF
        IF (RP.GT.ZERO) THEN
          IF (R.LE.RPMIN) THEN
            IF (RP.LT.RMIN) IRP = 0
            RPMIN = MIN(RP,RPMIN)
            RMIN = MIN(R,RMIN)
            IF (ABS(WW).GT.PVMIN) THEN
              IRP = IRP + 1
              IR(IRP) = I
            END IF
          END IF
        END IF
   60 CONTINUE

C Look for a swap of explicit bounds
      IF (JIN.LE.KB) THEN
        R = ZERO
        IF (SIGN.GT.0) THEN
          IF (BND(2,JIN).NE.SIGMA) THEN
            IF (JX(JIN).EQ.1) THEN
              R = BND(2,JIN) - X(JIN) - BND(1,JIN)
            ELSE IF (JX(JIN).EQ.0) THEN
              R = BND(2,JIN) - X(JIN)
            END IF
          END IF
        ELSE
          IF (BND(1,JIN).NE.-SIGMA) THEN
            IF (JX(JIN).EQ.2) THEN
              R = X(JIN) + BND(2,JIN) - BND(1,JIN)
            ELSE IF (JX(JIN).EQ.0) THEN
              R = X(JIN) - BND(1,JIN)
            END IF
          END IF
        END IF
        IF (R.GT.ZERO) THEN
          IF (R.LE.RPMIN) THEN
            RMIN = R
            JOUT = JIN
            PIVOT = ZERO
C   Favour this choice of JOUT
            GO TO 90
          END IF
        END IF
      END IF

C     Second row scan to resolve any ambiguity.
      PIVOT = ZERO
      DO 70 K = 1,IRP
        I = IR(K)
        WW = W(I)*SIGN
        J = IX(I)
        IF (J.GT.NPM) THEN
          IF (J.GT.NPM+NPM) THEN
            XJ = X(J-NPM-NPM)
          ELSE
            XJ = X(J-NPM)
          END IF
        ELSE
          XJ = X(J)
        END IF
        IF (J.GT.N) THEN
          R = XJ/WW
C            Favour removal of artificials in phase 1
          IF (.NOT.FEAS) WW = WW*10000.0
        ELSE IF (WW.LE.ZERO) THEN
          R = (XJ-BND(2,J))/WW
        ELSE IF (J.GE.LB) THEN
          R = XJ/WW
        ELSE
          R = (XJ-BND(1,J))/WW
        END IF
        IF (R.GT.RPMIN) GO TO 70
        IF (ABS(WW).LE.ABS(PIVOT)) GO TO 70
        PIVOT = WW
        RMIN = R
        JP = I
        JOUT = J
   70 CONTINUE

      IF (PIVOT.EQ.ZERO) THEN
C     The solution is unbounded
        IF (MP.GE.0) WRITE (MP,FMT='(/,/,A,1P,2E12.4)')
     +      ' LA04: approximate and newly computed reduced costs are',
     +      Z(JIN),ZQ
        IF (JOB.NE.3) THEN
          JOB = -2
C Replace small entries in W by zero
          DO 80 I = 1,M
            IF (ABS(W(I)).LE.PVMIN) W(I) = ZERO
   80     CONTINUE
        ELSE
          IF (MP.GE.0) WRITE (MP,FMT='(/,/,A)')
     +        ' LA04: no row is acceptable as the pivot row'
          JOB = 5
        END IF
        GO TO 280
      END IF

      PIVOT = W(JP)

   90 SIGNX = 0
      IF (JOUT.GT.NPM) THEN
        SIGNX = -1
        JOUT = JOUT - NPM
        IF (JOUT.GT.NPM) THEN
          SIGNX = 1
          JOUT = JOUT - NPM
        END IF
        NINF = NINF - 1
      END IF

      RMIN = MAX(ZERO,RMIN)*SIGN
C     Revise solution.
      DO 100 I = 1,M
        J = IX(I)
        IF (J.GT.NPM) THEN
          J = J - NPM
          IF (J.GT.NPM) J = J - NPM
        END IF
        X(J) = X(J) - W(I)*RMIN
  100 CONTINUE

      IF (JIN.NE.JOUT) THEN
C     Modify the LU factorization of  B corresponding to the
C     replacement of col. IX(JP) by col. IQ
        IF (M.GT.0) THEN
          call cpu_time(T)
          IF (CNTL(4).EQ.ONE) CALL LA15BD(BA,INDB,IB,M,IPB,IWB,V,
     +                             RINFO(3),W,.TRUE.,I15,I15(K15+1))
          call cpu_time(DUM)
          RINFO(6) = DUM - T + RINFO(6)
          T = DUM
          G5 = RINFO(3)
          C15(1) = CNTL(10)*RINFO(23)
          CALL LA15CD(BA,INDB,IB,M,IPB,IWB,Y,G5,RINFO(2),JP,I15,
     +                C15,I15(K15+1))
          RINFO(3) = G5
          call cpu_time(DUM)
          RINFO(7) = DUM - T + RINFO(7)
          IF (G5.GE.RINFO(24)) THEN
            IF (MP.GE.0) WRITE (MP,FMT='(/,/,A)')
     +          ' LA04: refactorization because of LA15 growth'
            JOB = 8
            GO TO 280
          ELSE IF (G5.LT.-6.5) THEN
            IF (IK.LE.5) THEN
              JOB = -4
              RINFO(35) = NINT(IB*1.5)
              GO TO 280
            END IF
            IF (MP.GE.0) WRITE (MP,FMT='(/,/,A,F7.1)')
     +          ' LA04: refactorization because of LA15 fill-in'
            JOB = 8
            GO TO 280
          ELSE IF (G5.LT.ZERO) THEN
C Increase CNTL(11) to stop this happening again
            CNTL(11) = ABS(PIVOT)*1.1/SQRT(GQ)
            IF (MP.GE.0) WRITE (MP,FMT='(/,/,A/,A,1PE10.2)')
     +          ' LA04: refactorization to avoid singular basis.',
     +          ' CNTL(11) increased to ', CNTL(11)
            JOB = 8
            GO TO 280
          END IF
        END IF
C     Calculate X(JIN)
        X(JIN) = X(JIN) + RMIN
        IF (JIN.LE.KB) THEN
          IF (JX(JIN).EQ.1) THEN
            X(JIN) = X(JIN) + BND(1,JIN)
          ELSE IF (JX(JIN).EQ.2) THEN
            X(JIN) = X(JIN) + BND(2,JIN)
          END IF
        END IF
      END IF

C Revise the value of the objective
      RINFO(1) = RINFO(1) + ZQ*RMIN

C     Increment iteration count and write summary
      ITER = ITER + 1
      call cpu_time(DUM)
      T = DUM - RINFO(9) + RINFO(4)
      IF (MP.GE.0) THEN
        IF (MOD(IK,10).EQ.0) WRITE (MP,FMT='(/,A/,A,A,/,A,A)') ' LA04:',
     +      '   Iter LENL   JIN     NINF       G        Z      ',
     +      '   GAMA       Time',
     +      '        LENU  JOUT  Obj. Fun     NCP    Aprox Z   ',
     +      'Aprox GAMA    Pivot'
        WRITE (MP,FMT='(3I6,I7,1P,E12.1,2E11.3,0P,F11.3)') ITER,
     +    I15(K15+1),JIN,NINF,RINFO(3),ZQ,GQ,T
        WRITE (MP,FMT='(I12,I6,1P,E14.6,I4,E12.3,E11.3,E11.1)')
     +    I15(K15+2),JOUT,RINFO(1),I15(K15+3),Z(JIN),G(JIN),PIVOT
      END IF

      IF (JIN.EQ.JOUT) THEN
        IF (SIGN.LT.ZERO) THEN
          JX(JOUT) = 1
        ELSE
          JX(JOUT) = 2
        END IF
        X(JOUT) = ZERO
        GO TO 170
      END IF

      IF (JOUT.LE.KB) THEN
        IF (SIGNX.EQ.ZERO) THEN
          IF (PIVOT*SIGN.GE.ZERO) THEN
            JX(JOUT) = 1
            X(JOUT) = X(JOUT) - BND(1,JOUT)
          ELSE
            JX(JOUT) = 2
            X(JOUT) = X(JOUT) - BND(2,JOUT)
          END IF
        END IF
      END IF
      IF (JIN.LE.KB) JX(JIN) = -1

C     Reclassify JIN as basic and IX(JP) as nonbasic
      IX(JP) = JIN
      G(JIN) = -JP

      DO 110 I = 1,M
        Y(I) = ZERO
  110 CONTINUE
      Y(JP) = ONE
      IF (M.GT.0) THEN
        call cpu_time(T)
        CALL LA15BD(BA,INDB,IB,M,IPB,IWB,V,RINFO(3),Y,.TRUE.,
     +              I15,I15(K15+1))
        call cpu_time(DUM)
        RINFO(6) = DUM - T + RINFO(6)
      END IF

      IF (CNTL(4).EQ.1) THEN
C     Revise column weights and reduced costs
        DO 140 I = 1,N
          IF (G(I).LE.ZERO) GO TO 140
          ALPHA = ZERO
          K1 = IP(I)
          K2 = IP(I+1) - 1
          DO 120 K = K1,K2
            KK = IRN(K)
            ALPHA = ALPHA + A(K)*Y(KK)
  120     CONTINUE
          IF (ALPHA.EQ.ZERO) GO TO 140
          BETA = ZERO
          DO 130 K = K1,K2
            KK = IRN(K)
            BETA = BETA + A(K)*W(KK)
  130     CONTINUE
          G(I) = MAX(G(I)+ALPHA* (ALPHA*GQ-2.*BETA),1.+ALPHA*ALPHA)
          Z(I) = Z(I) - ALPHA*ZQ
  140   CONTINUE
      ELSE
C     Revise reduced costs
        DO 160 I = 1,N
          IF (G(I).LE.ZERO) GO TO 160
          ALPHA = ZERO
          K1 = IP(I)
          K2 = IP(I+1) - 1
          DO 150 K = K1,K2
            KK = IRN(K)
            ALPHA = ALPHA + A(K)*Y(KK)
  150     CONTINUE
          Z(I) = Z(I) - ALPHA*ZQ
  160   CONTINUE
      END IF
      IF (JOUT.LE.N) THEN
        ZQ = -ZQ/PIVOT - SIGNX
        Z(JOUT) = ZQ
        G(JOUT) = MAX(GQ/ (PIVOT*PIVOT),ONE)
      END IF

C Check if it is time to refactorize
  170 IK = IK + 1
      call cpu_time(DUM)
      T = DUM - RINFO(9) + RINFO(4)
      J = MOD(IK,10) + 1
      RINFO(10+J) = (RINFO(8)+T-RINFO(10))/IK
      IF (IK.EQ.1) RINFO(21) = RINFO(10+J)
      RINFO(21) = MIN(RINFO(21),RINFO(10+J))
      IF (IK.GE.10) THEN
        ERR = ZERO
        DO 180 I = 2,10
          IF (J.NE.I-1) ERR = MAX(ERR,RINFO(10+I-1)-RINFO(10+I))
  180   CONTINUE
        IF (J.NE.10) ERR = MAX(ERR,RINFO(10+10)-RINFO(10+1))
        IF (RINFO(10+J).GT. (RINFO(21)+ERR)*1.05) THEN
          IF (MP.GE.0) WRITE (MP,FMT='(/,/,A)')
     +        ' LA04: refactorization to improve average iteration time'
          JOB = 4
          GO TO 280
        END IF
      END IF


C     Refine X, if this is required
      IF (CNTL(5).GE.1) THEN
        DO 190 I = 1,M
          Y(I) = X(N+I)
  190   CONTINUE
        DO 200 I = 1,M
          W(I) = B(I)
          J = IX(I)
          IF (J.GT.NPM) THEN
            J = J - NPM
            IF (J.GT.NPM) J = J - NPM
          END IF
          V(I) = X(J)
          IF (J.GT.N) Y(J-N) = ZERO
  200   CONTINUE
C Allow for variables at nonzero bounds
        DO 220 J = 1,N
          XJ = ZERO
          IF (G(J).GT.ZERO) XJ = X(J)
          IF (J.LE.KB) THEN
            IF (JX(J).EQ.1) THEN
              XJ = XJ + BND(1,J)
            ELSE IF (JX(J).EQ.2) THEN
              XJ = XJ + BND(2,J)
            END IF
          END IF
          IF (XJ.NE.ZERO) THEN
            DO 210 K = IP(J),IP(J+1) - 1
              I = IRN(K)
              W(I) = W(I) - A(K)*XJ
  210       CONTINUE
          END IF
  220   CONTINUE
        DO 230 J = 1,M
          W(J) = W(J) - Y(J)
  230   CONTINUE

        CALL LA04ED(A,LA,IRN,IP,M,N,W,V,DV,IX,0,RINFO,BA,INDB,IPB,IWB,
     +              YY,Y,IB,.FALSE.,I15,C15)

        DO 240 I = 1,M
          J = IX(I)
          IF (J.GT.NPM) THEN
            J = J - NPM
            IF (J.GT.NPM) J = J - NPM
          END IF
          X(J) = V(I)
  240   CONTINUE

      END IF

C Determine next pivot column
      JIN = 0
      JOB = 3
C Test for termination of phase 1
      IF (.NOT.FEAS) THEN
        IF (RINFO(1).LT.CNTL(3)) GO TO 280
      END IF
      ZGMAX = ZERO
      DO 250 I = 1,KB
        IF (G(I).LE.ZERO) GO TO 250
        ZI = Z(I)
        IF (JX(I).EQ.1) THEN
          IF (ZI.GT.ZERO) GO TO 250
        ELSE IF (JX(I).EQ.2) THEN
          IF (ZI.LT.ZERO) GO TO 250
        END IF
        ZI = ZI*ZI/G(I)
        IF (ZGMAX.LT.ZI) THEN
          ZGMAX = ZI
          JIN = I
        END IF
  250 CONTINUE
      DO 260 I = KB + 1,LB - 1
        IF (G(I).LE.ZERO) GO TO 260
        ZI = Z(I)
        ZI = ZI*ZI/G(I)
        IF (ZGMAX.LT.ZI) THEN
          ZGMAX = ZI
          JIN = I
        END IF
  260 CONTINUE
      DO 270 I = LB,N
        IF (G(I).LE.ZERO) GO TO 270
        ZI = Z(I)
        IF (ZI.GT.ZERO) GO TO 270
        ZI = ZI*ZI/G(I)
        IF (ZGMAX.LT.ZI) THEN
          ZGMAX = ZI
          JIN = I
        END IF
  270 CONTINUE

  280 RINFO(25) = JOUT
      RINFO(26) = JIN
      RINFO(27) = NINF
      RINFO(28) = IK
      RINFO(29) = ITER
      RINFO(32) = I15(K15+1) + I15(K15+2)
      RINFO(33) = I15(K15+3)

      END

      SUBROUTINE LA04ED(A,LA,IRN,IP,M,N,B,X,DX,IX,JIN,RINFO,BA,INDB,IPB,
     +                  IWB,W,Y,IB,TRANS,I15,C15)
C Perform iterative refinememnt.

C     .. Scalar Arguments ..
      INTEGER LA,IB,M,N
      LOGICAL TRANS
      DOUBLE PRECISION DX

C     .. Array Arguments ..
      DOUBLE PRECISION
     +                A(LA),BA(IB),B(M),RINFO(40),W(M),X(M),Y(M),C15(3)
      INTEGER INDB(IB,2),IP(N+1),IPB(M,2),IRN(LA),IWB(M,4),IX(M),JIN,
     +        I15(10)

C A,LA,...,N (IN) specify the matrix.
C B (IN) RHS if JIN = 0.
C X (INOUT) holds the current solution.
C DX (OUT) infinity norm of the final correction.
C IX (IN) holds the basis.
C JIN (IN) index of column of A to be the RHS or zero if B is RHS
C RINFO(3) (IN) contains the parameter G of LA15A/C.
C RINFO(6) (INOUT) holds the c.p.u. time spent in LA15B
C      since the last initial entry.
C BA,...,IWB are passed to LA15 and used to contain the factorized
C      basis.
C W is a workspace array.
C Y is a workspace array.
C IB (IN) holds the size of array B.
C TRANS (IN) true if the system involves the basis transpose
C I15(J),J=1,3 used for LA15's ICNTL argument
C I15(J),J=4,10 used for LA15's KEEP argument
C C15(J),J=1,3 used for LA15's CNTL argument

C     .. External Subroutines ..
      EXTERNAL LA15BD

C     .. Local Scalars ..
      INTEGER I,J,K,NPM,K15
      REAL SA,SX,T
      DOUBLE PRECISION
     +                DUM,OLDSUM,RA,RX,SUM,SUM2,TERM,ZERO
      PARAMETER (ZERO=0)
C DUM  Used in cpu_time call
C I    Row index
C J    Column index
C K    Index in array A
C K15  set to 3 where I15(K15+J),J=1,7 used for LA15's KEEP argument
C NPM  N+M
C OLDSUM Old value of SUM
C RA   Remaining part of an entry of A
C RX   Remaining part of an entry of X
C SA   Single precision part of an entry of A
C SUM  Temporary for main part of a sum
C SUM2 Temporary for second part of a sum
C SX   Single precision part of an entry of X
C T    Temporary for holding C.P.U. time
C TERM Term in acumulation

      K15 = 3

      NPM = N + M

C Compute residual
      IF (JIN.EQ.0) THEN
        DO 10 I = 1,M
          W(I) = B(I)
   10   CONTINUE
      ELSE
        DO 20 I = 1,M
          W(I) = ZERO
   20   CONTINUE
        DO 30 K = IP(JIN),IP(JIN+1) - 1
          W(IRN(K)) = A(K)
   30   CONTINUE
      END IF
      IF (TRANS) THEN
        DO 50 I = 1,M
          J = IX(I)
          IF (J.GT.N+M) THEN
            J = J - N - M
            IF (J.GT.N+M) J = J - N - M
          END IF
          IF (J.LE.N) THEN
C  Ordinary column
            SA = W(I)
            SUM = SA
            SUM2 = W(I) - SUM
            DO 40 K = IP(J),IP(J+1) - 1
C  Equivalent to SUM = SUM - A(K)*X(INRN(K))
              RX = -X(IRN(K))
              SX = RX
              RX = RX - SX
              RA = A(K)
              SA = RA
              RA = RA - SA
              TERM = DBLE(SX)*DBLE(SA)
              OLDSUM = SUM
              SUM = SUM + TERM
              TERM = (OLDSUM-SUM) + TERM
              SUM2 = SUM2 + RX*SA + SX*RA + RX*RA  + TERM
   40       CONTINUE
            W(I) = SUM + SUM2
          ELSE
C  Artificial column
            W(I) = W(I) - X(J-N)
          END IF
   50   CONTINUE
      ELSE
C Y is used for the second part of the residual
        DO 60 I = 1,M
          Y(I) = ZERO
   60   CONTINUE
        DO 80 I = 1,M
          IF (X(I).EQ.ZERO) GO TO 80
          RX = X(I)
          SX = RX
          RX = RX - SX
          J = IX(I)
          IF (J.GT.NPM) THEN
            J = J - NPM
            IF (J.GT.NPM) J = J - NPM
          END IF
          IF (J.LE.N) THEN
C  Ordinary column
            DO 70 K = IP(J),IP(J+1) - 1
C Equivalent to W(IRN(K)) = W(IRN(K)) - X(I)*A(K)
              RA = -A(K)
              SA = RA
              RA = RA - SA
              TERM = DBLE(SX)*DBLE(SA)
              OLDSUM = W(IRN(K))
              W(IRN(K)) = OLDSUM + TERM
              TERM = (OLDSUM-W(IRN(K))) + TERM
              Y(IRN(K)) = Y(IRN(K)) + RX*SA + SX*RA + RX*RA + TERM
   70       CONTINUE
          ELSE
C  Artificial column
C  Equivalent to W(J-N) = W(J-N) - X(I)
            TERM = -X(I)
            OLDSUM = W(J-N)
            W(J-N) = OLDSUM + TERM
            TERM = (OLDSUM-W(J-N)) + TERM
            Y(J-N) = Y(J-N) + TERM
          END IF
   80   CONTINUE
        DO 90 I = 1,M
          W(I) = W(I) + Y(I)
   90   CONTINUE

      END IF
      call cpu_time(T)
      CALL LA15BD(BA,INDB,IB,M,IPB,IWB,Y,RINFO(3),W,TRANS,I15,
     +            I15(K15+1))
      call cpu_time(DUM)
      RINFO(6) = RINFO(6) + DUM - T
      DX = ZERO
      DO 100 I = 1,M
        DX = MAX(DX,ABS(W(I)))
        X(I) = X(I) + W(I)
  100 CONTINUE

      END





      SUBROUTINE LA04SD(A,LA,IRN,IP,M,N,B,C,BND,KB,SIGMA,RS,CS,WS,JCN)

C Scale an LP problem prior to using LA04A

C     .. Scalar Arguments ..
      INTEGER LA,M,N,KB
      DOUBLE PRECISION
     +                SIGMA
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION
     +                 A(LA),B(M),BND(2,KB),C(N),CS(N+1),RS(M+1),
     +                 WS(2*M+3*N+5)
      INTEGER IP(N+1),IRN(LA),JCN(LA)

C A   (INOUT) holds the matrix entries by columns.
C LA  (IN) length of A, IRN, JCN.
C IRN (INOUT) holds the row indices of the matrix entries in A.
C IP  (IN) IP(J) holds the position in A and IRN of the first
C     entry of column J and IP(N+1) holds the first unused position.
C M   (IN) Number of rows.
C N   (IN) Number of columns.
C B   (INOUT) RHS vector.
C C   (INOUT) Objective vector.
C BND (INOUT) BND(1,J) and BND(2,J) hold the bounds for variable J,
C     J = 1,2,...,KB.
C KB  (IN) Second extent of BND.
C SIGMA (IN) Value that represents an infinite bound.
C RS (OUT) Row scaling factors
C CS (OUT) Column scaling factors
C WS   workarray
C JCN  workarray used to hold the column indices of the entries.

C     .. Local Scalars ..
      INTEGER I,J,K,NZ
      DOUBLE PRECISION
     +                ONE,ZERO
      PARAMETER (ONE=1,ZERO=0)
C I     Row index
C J     Column index
C K     Index in A and IRN.
C NZ    No. of entries in A.

C     .. Intrinsic Functions ..
      INTRINSIC ABS,EXP

C     .. External Subroutines ..
      EXTERNAL MC29AD

      NZ = IP(N+1) - 1

C Set up column indices
      DO 20 J = 1,N
        DO 10 K = IP(J),IP(J+1) - 1
          JCN(K) = J
   10   CONTINUE
   20 CONTINUE

C Insert b into A.
      DO 30 I = 1,M
        IF (B(I).NE.ZERO) THEN
          NZ = NZ + 1
          A(NZ) = B(I)
          IRN(NZ) = I
          JCN(NZ) = N + 1
        END IF
   30 CONTINUE

C Insert c into A.
      DO 40 J = 1,N
        IF (C(J).NE.ZERO) THEN
          NZ = NZ + 1
          A(NZ) = C(J)
          IRN(NZ) = M + 1
          JCN(NZ) = J
        END IF
   40 CONTINUE

C Find scaling factors
      IF (NZ.GT.0) THEN
        CALL MC29AD(M+1,N+1,NZ,A,IRN,JCN,RS,CS,WS,0,I)
        DO 50 I = 1,M + 1
          RS(I) = EXP(RS(I))
   50   CONTINUE
        DO 60 J = 1,N + 1
          CS(J) = EXP(CS(J))
   60   CONTINUE
        RS(M+1) = ONE/CS(N+1)
      ELSE
        DO 70 I = 1,M + 1
          RS(I) = ONE
   70   CONTINUE
        DO 80 J = 1,N + 1
          CS(J) = ONE
   80   CONTINUE
      END IF

C Scale the problem
      DO 90 I = 1,M
        B(I) = B(I)*RS(I)*CS(N+1)
   90 CONTINUE
      DO 100 J = 1,N
        C(J) = C(J)*CS(J)*RS(M+1)
  100 CONTINUE
      DO 110 I = 1,KB
        IF (ABS(BND(1,I)).NE.SIGMA) THEN
          BND(1,I) = BND(1,I)*CS(N+1)/CS(I)
        END IF
        IF (ABS(BND(2,I)).NE.SIGMA) THEN
          BND(2,I) = BND(2,I)*CS(N+1)/CS(I)
        END IF
  110 CONTINUE
      DO 120 K = 1,NZ
        A(K) = A(K)*RS(IRN(K))*CS(JCN(K))
  120 CONTINUE

      END


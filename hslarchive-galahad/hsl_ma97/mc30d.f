* COPYRIGHT (c) 1993 Council for the Central Laboratory
*                    of the Research Councils
* Original date 15 March 1993

* 12th July 2004 Version 1.0.0. Version numbering added.

      SUBROUTINE MC30AD(N,NE,A,IRN,ICN,S,W,LP,IFAIL)
      INTEGER N,NE
      DOUBLE PRECISION A(NE)
      INTEGER IRN(NE),ICN(NE)
      DOUBLE PRECISION S(N),W(N,4)
      INTEGER LP,IFAIL
C N is an integer variable that must be set to the matrix order.
C      It is not altered by the subroutine.
C NE is an integer variable that must be set to the number of entries.
C      It is not altered by the subroutine.
C A is an array that holds the values of the entries.
C IRN  is an integer array that must be set to the row indices of the
C      entries. It is not altered by the subroutine.
C ICN  is an integer array that must be set to the column indices of the
C      entries. It is not altered by the subroutine.
C S is an array that need not be be on entry. On return, it holds the
C      logarithms of the scaling factors.
C W is a workarray.
C      W(:,1)  holds row non-zero counts (diagonal matrix M).
C      W(:,2)  holds the residual vector r.
C      W(:,3)  holds the cg vector p.
C      W(:,4)  holds the cg vector (M+E)p.
C LP must be set to the unit number for messages.
C      It is not altered by the subroutine.
C IFAIL need not be set by the user. On return it has one of the
C     following values:
C     0 successful entry.
C     -1 N < 1.
C     -2 NE < 1.

      INTRINSIC LOG,ABS,MAX,MIN

C Constants
      INTEGER M,MAXIT,R,P,MP
      PARAMETER (M=1,MAXIT=10,MP=4,P=3,R=2)
      DOUBLE PRECISION ONE,RMIN,ZERO
      PARAMETER (ONE=1D0,RMIN=0.1,ZERO=0D0)
C M     W(:,M)  holds row non-zero counts (diagonal matrix M).
C MAXIT is the maximal permitted number of iterations.
C MP    W(:,MP)  holds the cg vector (M+E)p.
C P     W(:,P)  holds the cg vector p.
C R     W(:,R)  holds the residual vector.
C RMIN is used in a convergence test on (residual norm)**2

C Local variables
      DOUBLE PRECISION AK,BK
      INTEGER I,ITER,J,K
      DOUBLE PRECISION PP,RM,RR,RRL,U
C AK Scalar of cg iteration.
C BK Scalar of cg iteration.
C I Row index.
C ITER Iteration index.
C J Column index.
C K Entry number.
C PP Scalar p'(M+E)p of cg iteration.
C RM Threshold for RR.
C RR Scalar r'(inv M)r of cg iteration.
C RRL Previous value of RR.
C U abs(A(K)).

C Check N and NE.
      IFAIL = 0
      IF (N.LT.1) THEN
         IFAIL = -1
         GO TO 130
      ELSE IF (NE.LE.0) THEN
         IFAIL = -2
         GO TO 130
      END IF
C
C     Initialise for accumulation of sums and products.
      DO 10 I = 1,N
         S(I) = ZERO
         W(I,M) = ZERO
         W(I,R) = ZERO
   10 CONTINUE

C     Count non-zeros in the rows, and compute rhs vectors.
      DO 40 K = 1,NE
         U = ABS(A(K))
         IF (U.EQ.ZERO) GO TO 40
         I = IRN(K)
         J = ICN(K)
         IF ( MIN(I,J).LT.1 .OR. MAX(I,J).GT.N ) GO TO 40
         U = LOG(U)
         W(I,M) = W(I,M) + ONE
         W(I,R) = W(I,R) - U
         W(J,M) = W(J,M) + ONE
         IF(I.EQ.J) GO TO 40
         W(J,R) = W(J,R) - U
   40 CONTINUE

C     Find the initial vectors
      RR = ZERO
      DO 50 I = 1,N
         IF (W(I,M).EQ.ZERO) W(I,M) = ONE
         W(I,P) = W(I,R)/W(I,M)
         W(I,MP) = W(I,R)
         RR = RR + W(I,R)**2/W(I,M)
   50 CONTINUE

      RM = RMIN*NE
      IF (RR.LE.RM) RETURN
C
C     Iteration loop
      DO 120 ITER = 1,MAXIT
C    Sweep through matrix to add Ep to Mp
         DO 80 K = 1,NE
            IF (A(K).EQ.ZERO) GO TO 80
            I = ICN(K)
            J = IRN(K)
            IF(I.EQ.J) GO TO 80
            IF ( MIN(I,J).LT.1 .OR. MAX(I,J).GT.N ) GO TO 80
            W(I,MP) = W(I,MP) + W(J,P)
            W(J,MP) = W(J,MP) + W(I,P)
   80    CONTINUE
         PP = ZERO
         DO 90 I = 1,N
            PP = PP + W(I,P)*W(I,MP)
   90    CONTINUE
         AK = RR/PP
C     Update solution and residual
         RRL = RR
         RR = ZERO
         DO 100 I = 1,N
            S(I) = S(I) + AK*W(I,P)
            W(I,R) = W(I,R) - AK*W(I,MP)
            RR = RR + W(I,R)**2/W(I,M)
  100    CONTINUE
         IF (RR.LE.RM) RETURN
C Update vector P.
         BK = RR/RRL
         DO 110 I = 1,N
            W(I,P) = W(I,R)/W(I,M) + BK*W(I,P)
            W(I,MP) = W(I,P)*W(I,M)
  110    CONTINUE
  120 CONTINUE

C Error returns
  130 IF (LP.GT.0) WRITE (LP,'(/A/A,I3)')
     +    ' **** Error return from MC30AD ****',' IFAIL =',IFAIL
      END

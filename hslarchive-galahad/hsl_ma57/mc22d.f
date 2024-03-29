* COPYRIGHT (c) 1976 AEA Technology
* Original date 21 Jan 1993
C 8 August 2000: CONTINUEs given to DOs.
C 20/2/02 Cosmetic changes applied to reduce single/double differences

C
C 12th July 2004 Version 1.0.0. Version numbering added.

      SUBROUTINE MC22AD(N,ICN,A,NZ,LENROW,IP,IQ,IW,IW1)
C     .. Scalar Arguments ..
      INTEGER N,NZ
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(NZ)
      INTEGER ICN(NZ),IP(N),IQ(N),IW(N,2),IW1(NZ),LENROW(N)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION AVAL
      INTEGER I,ICHAIN,IOLD,IPOS,J,J2,JJ,JNUM,JVAL,LENGTH,NEWPOS
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC IABS
C     ..
C     .. Executable Statements ..
      IF (NZ.LE.0) GO TO 1000
      IF (N.LE.0) GO TO 1000
C SET START OF ROW I IN IW(I,1) AND LENROW(I) IN IW(I,2)
      IW(1,1) = 1
      IW(1,2) = LENROW(1)
      DO 10 I = 2,N
        IW(I,1) = IW(I-1,1) + LENROW(I-1)
        IW(I,2) = LENROW(I)
   10 CONTINUE
C PERMUTE LENROW ACCORDING TO IP.  SET OFF-SETS FOR NEW POSITION
C     OF ROW IOLD IN IW(IOLD,1) AND PUT OLD ROW INDICES IN IW1 IN
C     POSITIONS CORRESPONDING TO THE NEW POSITION OF THIS ROW IN A/ICN.
      JJ = 1
      DO 20 I = 1,N
        IOLD = IP(I)
        IOLD = IABS(IOLD)
        LENGTH = IW(IOLD,2)
        LENROW(I) = LENGTH
        IF (LENGTH.EQ.0) GO TO 20
        IW(IOLD,1) = IW(IOLD,1) - JJ
        J2 = JJ + LENGTH - 1
        DO 15 J = JJ,J2
          IW1(J) = IOLD
   15   CONTINUE
        JJ = J2 + 1
   20 CONTINUE
C SET INVERSE PERMUTATION TO IQ IN IW(.,2).
      DO 30 I = 1,N
        IOLD = IQ(I)
        IOLD = IABS(IOLD)
        IW(IOLD,2) = I
   30 CONTINUE
C PERMUTE A AND ICN IN PLACE, CHANGING TO NEW COLUMN NUMBERS.
C
C ***   MAIN LOOP   ***
C EACH PASS THROUGH THIS LOOP PLACES A CLOSED CHAIN OF COLUMN INDICES
C     IN THEIR NEW (AND FINAL) POSITIONS ... THIS IS RECORDED BY
C     SETTING THE IW1 ENTRY TO ZERO SO THAT ANY WHICH ARE SUBSEQUENTLY
C     ENCOUNTERED DURING THIS MAJOR SCAN CAN BE BYPASSED.
      DO 200 I = 1,NZ
        IOLD = IW1(I)
        IF (IOLD.EQ.0) GO TO 200
        IPOS = I
        JVAL = ICN(I)
C IF ROW IOLD IS IN SAME POSITIONS AFTER PERMUTATION GO TO 150.
        IF (IW(IOLD,1).EQ.0) GO TO 150
        AVAL = A(I)
C **  CHAIN LOOP  **
C EACH PASS THROUGH THIS LOOP PLACES ONE (PERMUTED) COLUMN INDEX
C     IN ITS FINAL POSITION  .. VIZ. IPOS.
        DO 100 ICHAIN = 1,NZ
C NEWPOS IS THE ORIGINAL POSITION IN A/ICN OF THE ELEMENT TO BE PLACED
C IN POSITION IPOS.  IT IS ALSO THE POSITION OF THE NEXT ELEMENT IN
C     THE CHAIN.
          NEWPOS = IPOS + IW(IOLD,1)
C IS CHAIN COMPLETE ?
          IF (NEWPOS.EQ.I) GO TO 130
          A(IPOS) = A(NEWPOS)
          JNUM = ICN(NEWPOS)
          ICN(IPOS) = IW(JNUM,2)
          IPOS = NEWPOS
          IOLD = IW1(IPOS)
          IW1(IPOS) = 0
C **  END OF CHAIN LOOP  **
  100   CONTINUE
  130   A(IPOS) = AVAL
  150   ICN(IPOS) = IW(JVAL,2)
C ***   END OF MAIN LOOP   ***
  200 CONTINUE
C
 1000 RETURN

      END

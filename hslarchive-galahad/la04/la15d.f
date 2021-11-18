! COPYRIGHT (c) 1975 AEA Technology and
! Council for the Central Laboratory of the Research Councils
!
! Version 1.3.0
! See ChangeLog for version history.
!
      SUBROUTINE LA15ID(ICNTL,CNTL,KEEP)
C     .. Array Arguments ..
      DOUBLE PRECISION CNTL(3)
      INTEGER ICNTL(3)
      INTEGER KEEP(7)
C     ..
C     .. Local Scalars ..
      INTEGER I
C     ..
      ICNTL(1) = 6
      ICNTL(2) = 0
      ICNTL(3) = 0

      CNTL(1) = 0.0
      CNTL(2) = 0.0
      CNTL(3) = 0.0
C
C  Initialize kept data to avoid undefined assignment
      DO 10 I = 1,7
      KEEP(I) = 0
   10 CONTINUE

      RETURN
      END
      SUBROUTINE LA15AD(A,IND,NZERO,IA,N,IP,IW,W,G,U,ICNTL,CNTL,KEEP)
C     .. Scalar Arguments ..
      DOUBLE PRECISION G,U
      INTEGER IA,N,NZERO
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(IA),W(N)
      INTEGER IND(IA,2),IP(N,2),IW(N,8)
      DOUBLE PRECISION CNTL(3)
      INTEGER ICNTL(3)
      INTEGER KEEP(7)
C IP(I,1),IP(I,2) POINT TO THE START OF ROW/COL I.
C IW(I,1),IW(I,2) HOLD THE NUMBER OF NON-ZEROS IN ROW/COL I.
C DURING THE MAIN BODY OF THIS SUBROUTINE THE VECTORS IW(.,3),IW(.,5),
C     IW(.,7) ARE USED TO HOLD DOUBLY LINKED LISTS OF ROWS THAT HAVE
C     NOT BEEN PIVOTAL AND HAVE EQUAL NUMBERS OF NON-ZEROS.
C IW(.,4),IW(.,6),IW(.,8) HOLD SIMILAR LISTS FOR THE COLUMNS.
C IW(I,3),IW(I,4) HOLD FIRST ROW/COLUMN TO HAVE I NON-ZEROS
C     OR ZERO IF THERE ARE NONE.
C IW(I,5), IW(I,6) HOLD ROW/COL NUMBER OF ROW/COL PRIOR TO ROW/COL I
C     IN ITS LIST, OR ZERO IF NONE.
C IW(I,7), IW(I,8) HOLD ROW/COL NUMBER OF ROW/COL AFTER ROW/COL I
C     IN ITS LIST, OR ZERO IF NONE.
C FOR ROWS/COLS THAT HAVE BEEN PIVOTAL IW(I,5),IW(I,6) HOLD NEGATION OF
C     POSITION OF ROW/COL I IN THE PIVOTAL ORDERING.
C ICNT59 CONTROLS FOR CALL TO MC59
C INFO59 INFORMATION RETURNED BY MC59 (MOSTLY IGNORED)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION AM,AMAX,AU,EPS,JCOST,KCOST,NZ1
      INTEGER I,IDUMMY,II,IL,IN,IPP,IPV,IR,J,JP,K,K1,K2,KC,KJ,KK,KL,KLC,
     +        KN,KNP,KP,KPC,KPL,KQ,KR,KRL,KS,L,MCP,NC,NZ,NZC
      INTEGER LP,LENL,LENU,NCP,LROW,LCOL,ICNT59(10),INFO59(10)
      REAL SMALL
C     ..
C     .. External Functions ..
      DOUBLE PRECISION FD15AD
      EXTERNAL FD15AD
C     ..
C     .. External Subroutines ..
      EXTERNAL LA15ED,MC59AD
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MAX
C     ..
C     .. Executable Statements ..
C
C  Restore kept data
C
      LP    = ICNTL(1)
      SMALL = CNTL(1)
      LENL  = KEEP(1)
      LENU  = KEEP(2)
      NCP   = KEEP(3)
      LROW  = KEEP(4)
      LCOL  = KEEP(5)
C
      EPS = FD15AD('E')
C EPS IS THE RELATIVE ACCURACY OF FLOATING-POINT COMPUTATION
      IF (U.GT.1.0D0) U = 1.0D0
      IF (U.LT.EPS) U = EPS
      IF (N.LT.1) GO TO 640
      IF (NZERO.LT.0) GO TO 700
      G = 0.
      DO 20 I = 1,N
        W(I) = 0.
        DO 10 J = 1,5
          IW(I,J) = 0
   10   CONTINUE
   20 CONTINUE
C
C FLUSH OUT SMALL ENTRIES, COUNT ELEMENTS IN ROWS AND COLUMNS
      L = 1
      LENU = NZERO
      DO 50 IDUMMY = 1,NZERO
        IF (L.GT.LENU) GO TO 60
        DO 30 K = L,LENU
          IF (ABS(A(K)).LE.SMALL) GO TO 40
          I = IND(K,1)
          J = IND(K,2)
          G = MAX(ABS(A(K)),G)
          IF (I.LT.1 .OR. I.GT.N) GO TO 650
          IF (J.LT.1 .OR. J.GT.N) GO TO 650
          IW(I,1) = IW(I,1) + 1
          IW(J,2) = IW(J,2) + 1
   30   CONTINUE
        GO TO 60

   40   L = K
        A(L) = A(LENU)
        IND(L,1) = IND(LENU,1)
        IND(L,2) = IND(LENU,2)
        LENU = LENU - 1
   50 CONTINUE
C
   60 LENL = 0
      LROW = LENU
      LCOL = LROW
C MCP IS THE MAXIMUM NUMBER OF COMPRESSES PERMITTED BEFORE AN
C     ERROR RETURN RESULTS.
      MCP = MAX(N/10,20)
      NCP = 0

      IF(LENU.GT.0)THEN
C Reorder by rows. Note that MC59 alters IP(1,2).
        ICNT59(1) = 1
        ICNT59(2) = 0
        ICNT59(3) = 0
        ICNT59(4) = ICNTL(1)
        IF ( ICNT59(4).EQ.0 ) ICNT59(4) = -1
        ICNT59(5) = ICNT59(4)
        ICNT59(6) = 0
        CALL MC59AD(ICNT59,N,N,LENU,IND(1,2),IA,IND(1,1),IA,A,N+1,IP,
     +            N+1,IW(1,7),INFO59)
      ELSE
        DO 70 IR = 1,N
          IP(IR,1) = 1
   70   CONTINUE
      END IF
C INITIALIZE IP(I,2) TO POINT JUST BEYOND WHERE THE LAST COMPONENT
C     OF COLUMN I OF A WILL BE STORED.
      K = 1
      DO 80 IR = 1,N
        K = K + IW(IR,2)
        IP(IR,2) = K
   80 CONTINUE
C CHECK FOR DOUBLE ENTRIES WHILE USING THE NEWLY CONSTRUCTED
C     ROW FILE TO CONSTRUCT THE COLUMN FILE. NOTE THAT BY PUTTING
C    THE ENTRIES IN BACKWARDS AND DECREASING IP(J,2) EACH TIME IT
C     IS USED WE AUTOMATICALLY LEAVE IT POINTING TO THE FIRST ELEMENT.
      KL = LENU
      DO 100 II = 1,N
        IR = N + 1 - II
        KP = IP(IR,1)
        DO 90 K = KP,KL
          J = IND(K,2)
          IF (IW(J,5).EQ.IR) GO TO 630
          IW(J,5) = IR
          KR = IP(J,2) - 1
          IP(J,2) = KR
          IND(KR,1) = IR
   90   CONTINUE
        KL = KP - 1
  100 CONTINUE
C
C SET UP LINKED LISTS OF ROWS AND COLS WITH EQUAL NUMBERS OF NON-ZEROS.
      DO 120 L = 1,2
        DO 110 I = 1,N
          NZ = IW(I,L)
          IF(NZ.GT.0)THEN
            IN = IW(NZ,L+2)
            IW(NZ,L+2) = I
            IW(I,L+6) = IN
            IF (IN.NE.0) IW(IN,L+4) = I
          END IF
          IW(I,L+4) = 0
  110   CONTINUE
  120 CONTINUE
C
C
C START OF MAIN ELIMINATION LOOP.
      DO 560 IPV = 1,N
C FIND PIVOT. JCOST IS MARKOWITZ COST OF CHEAPEST PIVOT FOUND SO FAR,
C     WHICH IS IN ROW IPP AND COLUMN JP.
C NZ1 holds NZ-1
        JCOST = N
        JCOST = JCOST*JCOST
C LOOP ON LENGTH OF COLUMN TO BE SEARCHED
        DO 210 NZ = 1,N
          NZ1 = NZ - 1
          IF (JCOST.LE.NZ1**2) GO TO 220
          J = IW(NZ,4)
C SEARCH COLUMNS WITH NZ NON-ZEROS.
          DO 160 IDUMMY = 1,N
            IF (J.LE.0) GO TO 170
            KP = IP(J,2)
            KL = KP + IW(J,2) - 1
            DO 150 K = KP,KL
              I = IND(K,1)
              KCOST = NZ1* (IW(I,1)-1)
              IF (KCOST.GE.JCOST) GO TO 150
              IF (NZ.EQ.1) GO TO 140
C FIND LARGEST ELEMENT IN ROW OF POTENTIAL PIVOT.
              AMAX = 0.D0
              K1 = IP(I,1)
              K2 = IW(I,1) + K1 - 1
              DO 130 KK = K1,K2
                AMAX = MAX(AMAX,ABS(A(KK)))
                IF (IND(KK,2).EQ.J) KJ = KK
  130         CONTINUE
C PERFORM STABILITY TEST.
              IF (ABS(A(KJ)).LT.AMAX*U) GO TO 150
  140         JCOST = KCOST
              IPP = I
              JP = J
              IF (JCOST.LE.NZ1**2) GO TO 220
  150       CONTINUE
            J = IW(J,8)
  160     CONTINUE
C SEARCH ROWS WITH NZ NON-ZEROS.
  170     I = IW(NZ,3)
          DO 200 IDUMMY = 1,N
            IF (I.LE.0) GO TO 210
            AMAX = 0.D0
            KP = IP(I,1)
            KL = KP + IW(I,1) - 1
C FIND LARGEST ELEMENT IN THE ROW
            DO 180 K = KP,KL
              AMAX = MAX(ABS(A(K)),AMAX)
  180       CONTINUE
            AU = AMAX*U
            DO 190 K = KP,KL
C PERFORM STABILITY TEST.
              IF (ABS(A(K)).LT.AU) GO TO 190
              J = IND(K,2)
              KCOST = NZ1* (IW(J,2)-1)
              IF (KCOST.GE.JCOST) GO TO 190
              JCOST = KCOST
              IPP = I
              JP = J
              IF (JCOST.LE.NZ1**2) GO TO 220
  190       CONTINUE
            I = IW(I,7)
  200     CONTINUE
  210   CONTINUE
        GO TO 680
C
C PIVOT FOUND.
C REMOVE ROWS AND COLUMNS INVOLVED IN ELIMINATION FROM ORDERING VECTORS.
  220   KP = IP(JP,2)
        KL = IW(JP,2) + KP - 1
        DO 260 L = 1,2
          DO 250 K = KP,KL
            I = IND(K,L)
            IL = IW(I,L+4)
            IN = IW(I,L+6)
            IF (IL.EQ.0) GO TO 230
            IW(IL,L+6) = IN
            GO TO 240

  230       NZ = IW(I,L)
            IW(NZ,L+2) = IN
  240       IF (IN.GT.0) IW(IN,L+4) = IL
  250     CONTINUE
          KP = IP(IPP,1)
          KL = KP + IW(IPP,1) - 1
  260   CONTINUE
C STORE PIVOT
        IW(IPP,5) = -IPV
        IW(JP,6) = -IPV
C ELIMINATE PIVOTAL ROW FROM COLUMN FILE AND FIND PIVOT IN ROW FILE.
        DO 290 K = KP,KL
          J = IND(K,2)
          KPC = IP(J,2)
          IW(J,2) = IW(J,2) - 1
          KLC = KPC + IW(J,2)
          DO 270 KC = KPC,KLC
            IF (IPP.EQ.IND(KC,1)) GO TO 280
  270     CONTINUE
  280     IND(KC,1) = IND(KLC,1)
          IND(KLC,1) = 0
          IF (J.EQ.JP) KR = K
  290   CONTINUE
C BRING PIVOT TO FRONT OF PIVOTAL ROW.
        AU = A(KR)
        A(KR) = A(KP)
        A(KP) = AU
        IND(KR,2) = IND(KP,2)
        IND(KP,2) = JP
C
C PERFORM ELIMINATION ITSELF, LOOPING ON NON-ZEROS IN PIVOT COLUMN.
        NZC = IW(JP,2)
        IF (NZC.EQ.0) GO TO 520
        DO 510 NC = 1,NZC
          KC = IP(JP,2) + NC - 1
          IR = IND(KC,1)
C SEARCH NON-PIVOT ROW FOR ELEMENT TO BE ELIMINATED.
          KR = IP(IR,1)
          KRL = KR + IW(IR,1) - 1
          DO 300 KNP = KR,KRL
            IF (JP.EQ.IND(KNP,2)) GO TO 310
  300     CONTINUE
C BRING ELEMENT TO BE ELIMINATED TO FRONT OF ITS ROW.
  310     AM = A(KNP)
          A(KNP) = A(KR)
          A(KR) = AM
          IND(KNP,2) = IND(KR,2)
          IND(KR,2) = JP
          AM = -A(KR)/A(KP)
C COMPRESS ROW FILE UNLESS IT IS CERTAIN THAT THERE IS ROOM FOR NEW ROW.
          IF (LROW+IW(IR,1)+IW(IPP,1)+LENL.LE.IA) GO TO 320
          IF (NCP.GE.MCP .OR. LENU+IW(IR,1)+IW(IPP,1)+LENL.GT.
     +        IA) GO TO 670
          CALL LA15ED(A,IND(1,2),IP,N,IW,IA,.TRUE.,NCP,LROW,LCOL)
          KP = IP(IPP,1)
          KR = IP(IR,1)
  320     KRL = KR + IW(IR,1) - 1
          KQ = KP + 1
          KPL = KP + IW(IPP,1) - 1
C PLACE PIVOT ROW (EXCLUDING PIVOT ITSELF) IN W.
          IF (KQ.GT.KPL) GO TO 340
          DO 330 K = KQ,KPL
            J = IND(K,2)
            W(J) = A(K)
  330     CONTINUE
  340     IP(IR,1) = LROW + 1
C
C TRANSFER MODIFIED ELEMENTS.
          IND(KR,2) = 0
          KR = KR + 1
          IF (KR.GT.KRL) GO TO 400
          DO 390 KS = KR,KRL
            J = IND(KS,2)
            AU = A(KS) + AM*W(J)
            IND(KS,2) = 0
C IF ELEMENT IS VERY SMALL REMOVE IT FROM U.
            IF (ABS(AU).LE.SMALL) GO TO 350
            G = MAX(G,ABS(AU))
            LROW = LROW + 1
            A(LROW) = AU
            IND(LROW,2) = J
            GO TO 380

  350       LENU = LENU - 1
C REMOVE ELEMENT FROM COL FILE.
            K = IP(J,2)
            KL = K + IW(J,2) - 1
            IW(J,2) = KL - K
            DO 360 KK = K,KL
              IF (IND(KK,1).EQ.IR) GO TO 370
  360       CONTINUE
  370       IND(KK,1) = IND(KL,1)
            IND(KL,1) = 0
  380       W(J) = 0.D0
  390     CONTINUE
C
C SCAN PIVOT ROW FOR FILLS.
  400     IF (KQ.GT.KPL) GO TO 490
          DO 480 KS = KQ,KPL
            J = IND(KS,2)
            AU = AM*W(J)
            IF (ABS(AU).LE.SMALL) GO TO 470
            LROW = LROW + 1
            A(LROW) = AU
            IND(LROW,2) = J
            LENU = LENU + 1
C
C CREATE FILL IN COLUMN FILE.
            NZ = IW(J,2)
            K = IP(J,2)
            KL = K + NZ - 1
C IF POSSIBLE PLACE NEW ELEMENT AT END OF PRESENT ENTRY.
            IF (KL.NE.LCOL) GO TO 410
            IF (LCOL+LENL.GE.IA) GO TO 430
            LCOL = LCOL + 1
            GO TO 420

  410       IF (IND(KL+1,1).NE.0) GO TO 430
  420       IND(KL+1,1) = IR
            GO TO 460
C NEW ENTRY HAS TO BE CREATED.
  430       IF (LCOL+LENL+NZ+1.LT.IA) GO TO 440
C COMPRESS COLUMN FILE IF THERE IS NOT ROOM FOR NEW ENTRY.
            IF (NCP.GE.MCP .OR. LENU+LENL+NZ+1.GE.IA) GO TO 670
            CALL LA15ED(A,IND,IP(1,2),N,IW(1,2),IA,.FALSE.,NCP,
     +                  LROW,LCOL)
            K = IP(J,2)
            KL = K + NZ - 1
C TRANSFER OLD ENTRY INTO NEW.
  440       IP(J,2) = LCOL + 1
            DO 450 KK = K,KL
              LCOL = LCOL + 1
              IND(LCOL,1) = IND(KK,1)
              IND(KK,1) = 0
  450       CONTINUE
C ADD NEW ELEMENT.
            LCOL = LCOL + 1
            IND(LCOL,1) = IR
  460       G = MAX(G,ABS(AU))
            IW(J,2) = NZ + 1
  470       W(J) = 0.D0
  480     CONTINUE
  490     IW(IR,1) = LROW + 1 - IP(IR,1)
C
C STORE MULTIPLIER
          IF (LENL+LCOL+1.LE.IA) GO TO 500
C COMPRESS COL FILE IF NECESSARY.
          IF (NCP.GE.MCP) GO TO 670
          CALL LA15ED(A,IND,IP(1,2),N,IW(1,2),IA,.FALSE.,NCP,
     +                LROW,LCOL)
  500     K = IA - LENL
          LENL = LENL + 1
          A(K) = AM
          IND(K,1) = IPP
          IND(K,2) = IR
          LENU = LENU - 1
  510   CONTINUE
C
C INSERT ROWS AND COLUMNS INVOLVED IN ELIMINATION IN LINKED LISTS
C     OF EQUAL NUMBERS OF NON-ZEROS.
  520   K1 = IP(JP,2)
        K2 = IW(JP,2) + K1 - 1
        IW(JP,2) = 0
        DO 550 L = 1,2
          IF (K2.LT.K1) GO TO 540
          DO 530 K = K1,K2
            IR = IND(K,L)
            IF (L.EQ.1) IND(K,L) = 0
            NZ = IW(IR,L)
            IF (NZ.GT.0) THEN
              IN = IW(NZ,L+2)
              IW(IR,L+6) = IN
              IW(NZ,L+2) = IR
              IF (IN.NE.0) IW(IN,L+4) = IR
            END IF
            IW(IR,L+4) = 0
  530     CONTINUE
  540     K1 = IP(IPP,1) + 1
          K2 = IW(IPP,1) + K1 - 2
  550   CONTINUE
  560 CONTINUE
C
C RESET COLUMN FILE TO REFER TO U AND STORE ROW/COL NUMBERS IN
C     PIVOTAL ORDER IN IW(.,3),IW(.,4)
      DO 570 I = 1,N
        J = -IW(I,5)
        IW(J,3) = I
        J = -IW(I,6)
        IW(J,4) = I
        IW(I,2) = 0
  570 CONTINUE
      DO 590 I = 1,N
        KP = IP(I,1)
        KL = IW(I,1) + KP - 1
        DO 580 K = KP,KL
          J = IND(K,2)
          IW(J,2) = IW(J,2) + 1
  580   CONTINUE
  590 CONTINUE
      K = 1
      DO 600 I = 1,N
        K = K + IW(I,2)
        IP(I,2) = K
  600 CONTINUE
      LCOL = K - 1
      DO 620 II = 1,N
        I = IW(II,3)
        KP = IP(I,1)
        KL = IW(I,1) + KP - 1
        DO 610 K = KP,KL
          J = IND(K,2)
          KN = IP(J,2) - 1
          IP(J,2) = KN
          IND(KN,1) = I
  610   CONTINUE
  620 CONTINUE
      GO TO 710
C
C     THE FOLLOWING INSTRUCTIONS IMPLEMENT THE FAILURE EXITS.
  630 IF (LP.GT.0) WRITE (LP,FMT=9000) IR,J

 9000 FORMAT (/,/,' ERROR RETURN FROM LA15AD BECAUSE',/,' THERE IS MOR',
     +       'E THAN ONE ENTRY IN ROW',I5,' AND COLUMN',I5)

      G = -4.
      GO TO 710

  640 IF (LP.GT.0) WRITE (LP,FMT=9010)

 9010 FORMAT (/,/,' ERROR RETURN FROM LA15AD BECAUSE N IS NOT POSITIVE')

      G = -1.0D0
      GO TO 710

  650 IF (LP.GT.0) WRITE (LP,FMT=9020) K,I,J

 9020 FORMAT (/,/,' ERROR RETURN FROM LA15AD BECAUSE',/,' ELEMENT',I7,
     +       ' IS IN ROW',I5,' AND COLUMN',I5)

      G = -3.D0
      GO TO 710

  670 IF (LP.GT.0) WRITE (LP,FMT=9040)

 9040 FORMAT (/,/,' ERROR RETURN FROM LA15AD BECAUSE IA IS TOO SMALL')

      G = -7.D0
      GO TO 710

  680 IPV = IPV - 1
      IF (LP.GT.0) WRITE (LP,FMT=9050) IPV
 9050 FORMAT (/,/,' ERROR RETURN FROM LA15AD BECAUSE THE MATRIX IS ',
     +   'SINGULAR WITH RANK ',I7)
      G = -5.D0
      GO TO 710

  700 IF (LP.GT.0) WRITE (LP,FMT=9060)

 9060 FORMAT (/,/,' ERROR RETURN FROM LA15AD BECAUSE NZ IS NEGATIVE')

      G = -8.
C
C  Save kept data and return to caller
C
  710 CONTINUE
      KEEP(1)  = LENL
      KEEP(2)  = LENU
      KEEP(3)  = NCP
      KEEP(4)  = LROW
      KEEP(5)  = LCOL
      RETURN
      END
      SUBROUTINE LA15BD(A,IND,IA,N,IP,IW,W,G,B,TRANS,ICNTL,KEEP)
C IP(I,1),IP(I,2) POINT TO START OF ROW/COLUMN I OF U.
C IW(I,1),IW(I,2) ARE LENGTHS OF ROW/COL I OF U.
C IW(.,3),IW(.,4) HOLD ROW/COL NUMBERS IN PIVOTAL ORDER.
C     .. Scalar Arguments ..
      DOUBLE PRECISION G
      INTEGER IA,N
      LOGICAL TRANS
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(IA),B(N),W(N)
      INTEGER IND(IA,2),IP(N,2),IW(N,4)
      INTEGER ICNTL(3)
      INTEGER KEEP(7)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION AM
      INTEGER I,II,J,K,K2,KK,KL,KLL,KP,KPC,L1,N1,NZ,LP,LENL
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
C     .. Executable Statements ..
C
C  Restore kept data (not altered by LA15BD)
C
      LP    = ICNTL(1)
      LENL  = KEEP(1)
C
      IF (G.LT.0.D0) GO TO 130
      KLL = IA - LENL + 1
      IF (TRANS) GO TO 80
C
C     MULTIPLY VECTOR BY INVERSE OF L
      IF (LENL.LE.0) GO TO 20
      L1 = IA + 1
      DO 10 KK = 1,LENL
        K = L1 - KK
        I = IND(K,1)
        IF (B(I).EQ.0.D0) GO TO 10
        J = IND(K,2)
        B(J) = B(J) + A(K)*B(I)
   10 CONTINUE
   20 DO 30 I = 1,N
        W(I) = B(I)
        B(I) = 0.D0
   30 CONTINUE
C
C     MULTIPLY VECTOR BY INVERSE OF U
      N1 = N + 1
      DO 70 II = 1,N
        I = N1 - II
        I = IW(I,3)
        AM = W(I)
        KP = IP(I,1)
        IF (KP.GT.0) GO TO 50
        KP = -KP
        IP(I,1) = KP
        NZ = IW(I,1)
        KL = KP - 1 + NZ
        K2 = KP + 1
        DO 40 K = K2,KL
          J = IND(K,2)
          AM = AM - A(K)*B(J)
   40   CONTINUE
   50   IF (AM.EQ.0.D0) GO TO 70
        J = IND(KP,2)
        B(J) = AM/A(KP)
        KPC = IP(J,2)
        KL = IW(J,2) + KPC - 1
        IF (KL.EQ.KPC) GO TO 70
        K2 = KPC + 1
        DO 60 K = K2,KL
          I = IND(K,1)
          IP(I,1) = -ABS(IP(I,1))
   60   CONTINUE
   70 CONTINUE
      GO TO 140
C
C     MULTIPLY VECTOR BY INVERSE OF TRANSPOSE OF U
   80 DO 90 I = 1,N
        W(I) = B(I)
        B(I) = 0.D0
   90 CONTINUE
      DO 110 II = 1,N
        I = IW(II,4)
        AM = W(I)
        IF (AM.EQ.0.D0) GO TO 110
        J = IW(II,3)
        KP = IP(J,1)
        AM = AM/A(KP)
        B(J) = AM
        KL = IW(J,1) + KP - 1
        IF (KP.EQ.KL) GO TO 110
        K2 = KP + 1
        DO 100 K = K2,KL
          I = IND(K,2)
          W(I) = W(I) - AM*A(K)
  100   CONTINUE
  110 CONTINUE
C
C     MULTIPLY VECTOR BY INVERSE OF TRANSPOSE OF L
      IF (KLL.GT.IA) RETURN
      DO 120 K = KLL,IA
        J = IND(K,2)
        IF (B(J).EQ.0.D0) GO TO 120
        I = IND(K,1)
        B(I) = B(I) + A(K)*B(J)
  120 CONTINUE
      GO TO 140
C
  130 IF (LP.GT.0) WRITE (LP,FMT=9000)

 9000 FORMAT (/,/,' ERROR RETURN FROM LA15BD BECAUSE EARLIER ENTRY',
     +       ' GAVE ERROR RETURN')
  140 RETURN
      END
      SUBROUTINE LA15CD(A,IND,IA,N,IP,IW,W,G,U,MM,ICNTL,CNTL,KEEP)
C     .. Scalar Arguments ..
      DOUBLE PRECISION G,U
      INTEGER IA,MM,N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(IA),W(N)
      INTEGER IND(IA,2),IP(N,2),IW(N,4)
      DOUBLE PRECISION CNTL(3)
      INTEGER ICNTL(3)
      INTEGER KEEP(7)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION AM,AU
      INTEGER I,II,IJ,IM,IN,INS,IPP,IR,IS,J,JM,JNS,JP,K,KJ,KK,KL,KM,KNP,
     +        KP,KPL,KQ,KR,KRL,KS,L,LAST,LAST1,LAST2,M,M1,MCP,NZ
      INTEGER LP,LENL,LENU,NCP,LROW,LCOL
      REAL SMALL
C     ..
C     .. External Subroutines ..
      EXTERNAL LA15ED
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MAX
C     ..
C     .. Executable Statements ..
C
C  Restore kept data
C
      LP    = ICNTL(1)
      SMALL = CNTL(1)
      LENL  = KEEP(1)
      LENU  = KEEP(2)
      NCP   = KEEP(3)
      LROW  = KEEP(4)
      LCOL  = KEEP(5)
C
      IF (G.LT.0.D0) GO TO 620
      IF (MM.LT.1 .OR. MM.GT.N) GO TO 630
      JM = MM
C MCP LIMITS THE VALUE OF NCP PERMITTED BEFORE AN ERROR RETURN RESULTS.
      MCP = NCP + 20
C REMOVE OLD COLUMN
      LENU = LENU - IW(JM,2)
      KP = IP(JM,2)
      IM = IND(KP,1)
      KL = KP + IW(JM,2) - 1
      IW(JM,2) = 0
      DO 30 K = KP,KL
        I = IND(K,1)
        IND(K,1) = 0
        KR = IP(I,1)
        NZ = IW(I,1) - 1
        IW(I,1) = NZ
        KRL = KR + NZ
        DO 10 KM = KR,KRL
          IF (IND(KM,2).EQ.JM) GO TO 20
   10   CONTINUE
   20   A(KM) = A(KRL)
        IND(KM,2) = IND(KRL,2)
        IND(KRL,2) = 0
   30 CONTINUE
C
C INSERT NEW COLUMN
C JDH 09/04/13: Added initialization of LAST. I /think/ 0 is the correct
C value, but can't be certain. It seems to work.
      LAST = 0
      DO 110 II = 1,N
        I = IW(II,3)
        IF (I.EQ.IM) M = II
        IF (ABS(W(I)).LE.SMALL) GO TO 100
        LENU = LENU + 1
        LAST = II
        IF (LCOL+LENL.LT.IA) GO TO 40
C COMPRESS COLUMN FILE IF NECESSARY.
        IF (NCP.GE.MCP .OR. LENL+LENU.GE.IA) GO TO 610
        CALL LA15ED(A,IND,IP(1,2),N,IW(1,2),IA,.FALSE.,NCP,LROW,
     +              LCOL)
   40   LCOL = LCOL + 1
        NZ = IW(JM,2)
        IF (NZ.EQ.0) IP(JM,2) = LCOL
        IW(JM,2) = NZ + 1
        IND(LCOL,1) = I
        NZ = IW(I,1)
        KPL = IP(I,1) + NZ
        IF (KPL.GT.LROW) GO TO 50
        IF (IND(KPL,2).EQ.0) GO TO 90
C NEW ENTRY HAS TO BE CREATED.
   50   IF (LENL+LROW+NZ.LT.IA) GO TO 60
        IF (NCP.GE.MCP .OR. LENL+LENU+NZ.GE.IA) GO TO 610
C COMPRESS ROW FILE IF NECESSARY.
        CALL LA15ED(A,IND(1,2),IP,N,IW,IA,.TRUE.,NCP,LROW,LCOL)
   60   KP = IP(I,1)
        IP(I,1) = LROW + 1
        IF (NZ.EQ.0) GO TO 80
        KPL = KP + NZ - 1
        DO 70 K = KP,KPL
          LROW = LROW + 1
          A(LROW) = A(K)
          IND(LROW,2) = IND(K,2)
          IND(K,2) = 0
   70   CONTINUE
   80   LROW = LROW + 1
        KPL = LROW
C PLACE NEW ELEMENT AT END OF ROW.
   90   IW(I,1) = NZ + 1
        A(KPL) = W(I)
        IND(KPL,2) = JM
  100   W(I) = 0.D0
  110 CONTINUE
      IF (IW(IM,1).EQ.0 .OR. IW(JM,2).EQ.0 .OR. M.GT.LAST) GO TO 600
C
C FIND COLUMN SINGLETONS, OTHER THAN THE SPIKE. NON-SINGLETONS ARE
C     MARKED WITH W(J)=1. ONLY IW(.,3) IS REVISED AND IW(.,4) IS USED
C     FOR WORKSPACE.
      INS = M
      M1 = M
      W(JM) = 1.0D0
      DO 140 II = M,LAST
        I = IW(II,3)
        J = IW(II,4)
        IF (W(J).EQ.0.D0) GO TO 130
        KP = IP(I,1)
        KL = KP + IW(I,1) - 1
        DO 120 K = KP,KL
          J = IND(K,2)
          W(J) = 1.0D0
  120   CONTINUE
        IW(INS,4) = I
        INS = INS + 1
        GO TO 140
C PLACE SINGLETONS IN NEW POSITION.
  130   IW(M1,3) = I
        M1 = M1 + 1
  140 CONTINUE
C PLACE NON-SINGLETONS IN NEW POSITION.
      IJ = M + 1
      IF (M1.GE.LAST) GO TO 160
      LAST1 = LAST - 1
      DO 150 II = M1,LAST1
        IW(II,3) = IW(IJ,4)
        IJ = IJ + 1
  150 CONTINUE
C PLACE SPIKE AT END.
  160 IW(LAST,3) = IM
C
C FIND ROW SINGLETONS, APART FROM SPIKE ROW. NON-SINGLETONS ARE MARKED
C     WITH W(I)=2.D0 AGAIN ONLY IW(.,3) IS REVISED AND IW(.,4) IS USED
C     FOR WORKSPACE.
      LAST1 = LAST
      JNS = LAST
      W(IM) = 2.D0
      J = JM
      DO 190 IJ = M1,LAST
        II = LAST + M1 - IJ
        I = IW(II,3)
        IF (W(I).NE.2.D0) GO TO 180
        K = IP(I,1)
        IF (II.NE.LAST) J = IND(K,2)
        KP = IP(J,2)
        KL = KP + IW(J,2) - 1
        IW(JNS,4) = I
        JNS = JNS - 1
        DO 170 K = KP,KL
          I = IND(K,1)
          W(I) = 2.D0
  170   CONTINUE
        GO TO 190

  180   IW(LAST1,3) = I
        LAST1 = LAST1 - 1
  190 CONTINUE
      DO 200 II = M1,LAST1
        JNS = JNS + 1
        I = IW(JNS,4)
        W(I) = 3.D0
        IW(II,3) = I
  200 CONTINUE
C
C DEAL WITH SINGLETON SPIKE COLUMN. NOTE THAT BUMP ROWS ARE MARKED BY
C    W(I)=3.D0
      DO 240 II = M1,LAST1
        KP = IP(JM,2)
        KL = KP + IW(JM,2) - 1
        IS = 0
        DO 210 K = KP,KL
          L = IND(K,1)
          IF (W(L).NE.3.D0) GO TO 210
          IF (IS.NE.0) GO TO 250
          I = L
          KNP = K
          IS = 1
  210   CONTINUE
        IF (IS.EQ.0) GO TO 600
C MAKE A(I,JM) A PIVOT.
        IND(KNP,1) = IND(KP,1)
        IND(KP,1) = I
        KP = IP(I,1)
        DO 220 K = KP,IA
          IF (IND(K,2).EQ.JM) GO TO 230
  220   CONTINUE
  230   AM = A(KP)
        A(KP) = A(K)
        A(K) = AM
        IND(K,2) = IND(KP,2)
        IND(KP,2) = JM
        JM = IND(K,2)
        IW(II,4) = I
        W(I) = 2.D0
  240 CONTINUE
      II = LAST1
      GO TO 270

  250 IN = M1
      DO 260 IJ = II,LAST1
        IW(IJ,4) = IW(IN,3)
        IN = IN + 1
  260 CONTINUE
  270 LAST2 = LAST1 - 1
      IF (M1.EQ.LAST1) GO TO 580
      DO 280 I = M1,LAST2
        IW(I,3) = IW(I,4)
  280 CONTINUE
      M1 = II
      IF (M1.EQ.LAST1) GO TO 580
C
C CLEAR W
      DO 290 I = 1,N
        W(I) = 0.D0
  290 CONTINUE
C
C PERFORM ELIMINATION
      IR = IW(LAST1,3)
      DO 570 II = M1,LAST1
        IPP = IW(II,3)
        KP = IP(IPP,1)
        KR = IP(IR,1)
        JP = IND(KP,2)
        IF (II.EQ.LAST1) JP = JM
C SEARCH NON-PIVOT ROW FOR ELEMENT TO BE ELIMINATED.
C  AND BRING IT TO FRONT OF ITS ROW
        KRL = KR + IW(IR,1) - 1
        DO 300 KNP = KR,KRL
          IF (JP.EQ.IND(KNP,2)) GO TO 310
  300   CONTINUE
        IF(II-LAST1.NE.0) GO TO 570
        IF(II-LAST1.EQ.0) GO TO 600
C BRING ELEMENT TO BE ELIMINATED TO FRONT OF ITS ROW.
  310   AM = A(KNP)
        A(KNP) = A(KR)
        A(KR) = AM
        IND(KNP,2) = IND(KR,2)
        IND(KR,2) = JP
        IF (II.EQ.LAST1) GO TO 320
        IF (ABS(A(KP)).LT.U*ABS(AM)) GO TO 320
        IF (ABS(AM).LT.U*ABS(A(KP))) GO TO 350
        IF (IW(IPP,1).LE.IW(IR,1)) GO TO 350
C PERFORM INTERCHANGE
  320   IW(LAST1,3) = IPP
        IW(II,3) = IR
        IR = IPP
        IPP = IW(II,3)
        K = KR
        KR = KP
        KP = K
        KJ = IP(JP,2)
        DO 330 K = KJ,IA
          IF (IND(K,1).EQ.IPP) GO TO 340
  330   CONTINUE
  340   IND(K,1) = IND(KJ,1)
        IND(KJ,1) = IPP
  350   IF (A(KP).EQ.0.D0) GO TO 600
        IF (II.EQ.LAST1) GO TO 570
        AM = -A(KR)/A(KP)
C COMPRESS ROW FILE UNLESS IT IS CERTAIN THAT THERE IS ROOM FOR NEW ROW.
        IF (LROW+IW(IR,1)+IW(IPP,1)+LENL.LE.IA) GO TO 360
        IF (NCP.GE.MCP .OR. LENU+IW(IR,1)+IW(IPP,1)+LENL.GT.
     +      IA) GO TO 610
        CALL LA15ED(A,IND(1,2),IP,N,IW,IA,.TRUE.,NCP,LROW,LCOL)
        KP = IP(IPP,1)
        KR = IP(IR,1)
  360   KRL = KR + IW(IR,1) - 1
        KQ = KP + 1
        KPL = KP + IW(IPP,1) - 1
C PLACE PIVOT ROW (EXCLUDING PIVOT ITSELF) IN W.
        IF (KQ.GT.KPL) GO TO 380
        DO 370 K = KQ,KPL
          J = IND(K,2)
          W(J) = A(K)
  370   CONTINUE
  380   IP(IR,1) = LROW + 1
C
C TRANSFER MODIFIED ELEMENTS.
        IND(KR,2) = 0
        KR = KR + 1
        IF (KR.GT.KRL) GO TO 440
        DO 430 KS = KR,KRL
          J = IND(KS,2)
          AU = A(KS) + AM*W(J)
          IND(KS,2) = 0
C IF ELEMENT IS VERY SMALL REMOVE IT FROM U.
          IF (ABS(AU).LE.SMALL) GO TO 390
          G = MAX(G,ABS(AU))
          LROW = LROW + 1
          A(LROW) = AU
          IND(LROW,2) = J
          GO TO 420

  390     LENU = LENU - 1
C REMOVE ELEMENT FROM COL FILE.
          K = IP(J,2)
          KL = K + IW(J,2) - 1
          IW(J,2) = KL - K
          DO 400 KK = K,KL
            IF (IND(KK,1).EQ.IR) GO TO 410
  400     CONTINUE
  410     IND(KK,1) = IND(KL,1)
          IND(KL,1) = 0
  420     W(J) = 0.D0
  430   CONTINUE
C
C SCAN PIVOT ROW FOR FILLS.
  440   IF (KQ.GT.KPL) GO TO 530
        DO 520 KS = KQ,KPL
          J = IND(KS,2)
          AU = AM*W(J)
          IF (ABS(AU).LE.SMALL) GO TO 510
          LROW = LROW + 1
          A(LROW) = AU
          IND(LROW,2) = J
          LENU = LENU + 1
C
C CREATE FILL IN COLUMN FILE.
          NZ = IW(J,2)
          K = IP(J,2)
          KL = K + NZ - 1
C IF POSSIBLE PLACE NEW ELEMENT AT END OF PRESENT ENTRY.
          IF (KL.NE.LCOL) GO TO 450
          IF (LCOL+LENL.GE.IA) GO TO 470
          LCOL = LCOL + 1
          GO TO 460

  450     IF (IND(KL+1,1).NE.0) GO TO 470
  460     IND(KL+1,1) = IR
          GO TO 500
C NEW ENTRY HAS TO BE CREATED.
  470     IF (LCOL+LENL+NZ+1.LT.IA) GO TO 480
C COMPRESS COLUMN FILE IF THERE IS NOT ROOM FOR NEW ENTRY.
          IF (NCP.GE.MCP .OR. LENU+LENL+NZ+1.GE.IA) GO TO 610
          CALL LA15ED(A,IND,IP(1,2),N,IW(1,2),IA,.FALSE.,NCP,
     +                LROW,LCOL)
          K = IP(J,2)
          KL = K + NZ - 1
C TRANSFER OLD ENTRY INTO NEW.
  480     IP(J,2) = LCOL + 1
          DO 490 KK = K,KL
            LCOL = LCOL + 1
            IND(LCOL,1) = IND(KK,1)
            IND(KK,1) = 0
  490     CONTINUE
C ADD NEW ELEMENT.
          LCOL = LCOL + 1
          IND(LCOL,1) = IR
  500     G = MAX(G,ABS(AU))
          IW(J,2) = NZ + 1
  510     W(J) = 0.D0
  520   CONTINUE
  530   IW(IR,1) = LROW + 1 - IP(IR,1)
C
C STORE MULTIPLIER
        IF (LENL+LCOL+1.LE.IA) GO TO 540
C COMPRESS COL FILE IF NECESSARY.
        IF (NCP.GE.MCP) GO TO 610
        CALL LA15ED(A,IND,IP(1,2),N,IW(1,2),IA,.FALSE.,NCP,LROW,
     +              LCOL)
  540   K = IA - LENL
        LENL = LENL + 1
        A(K) = AM
        IND(K,1) = IPP
        IND(K,2) = IR
C CREATE BLANK IN PIVOTAL COLUMN.
        KP = IP(JP,2)
        NZ = IW(JP,2) - 1
        KL = KP + NZ
        DO 550 K = KP,KL
          IF (IND(K,1).EQ.IR) GO TO 560
  550   CONTINUE
  560   IND(K,1) = IND(KL,1)
        IW(JP,2) = NZ
        IND(KL,1) = 0
        LENU = LENU - 1
  570 CONTINUE
C
C CONSTRUCT COLUMN PERMUTATION AND STORE IT IN IW(.,4)
  580 DO 590 II = M,LAST
        I = IW(II,3)
        K = IP(I,1)
        J = IND(K,2)
        IW(II,4) = J
  590 CONTINUE
      GO TO 640
C
C     THE FOLLOWING INSTRUCTIONS IMPLEMENT THE FAILURE EXITS.
  600 IF (LP.NE.0) WRITE (LP,FMT=9000) MM

 9000 FORMAT (/,/,' ERROR RETURN FROM LA15CD BECAUSE',/,' SINGULAR MAT',
     +       'RIX CREATED BY REPLACEMENT OF COL',I5)

      G = -6.D0
      GO TO 640

  610 IF (LP.NE.0) WRITE (LP,FMT=9010)

 9010 FORMAT (/,/,' ERROR RETURN FROM LA15CD BECAUSE IA IS TOO SMALL')

      G = -7.D0
      GO TO 640

  620 IF (LP.NE.0) WRITE (LP,FMT=9020)

 9020 FORMAT (/,/,' ERROR RETURN FROM LA15CD BECAUSE EARLIER ENTRY',
     +       ' GAVE ERROR RETURN')

      GO TO 640

  630 IF (LP.NE.0) WRITE (LP,FMT=9030) MM

 9030 FORMAT (/,/,' ERROR RETURN FROM LA15CD BECAUSE M HAS THE VALUE',
     +       I8)

      G = -9.
C
C  Save kept data and return to caller
C
  640 CONTINUE
      KEEP(1)  = LENL
      KEEP(2)  = LENU
      KEEP(3)  = NCP
      KEEP(4)  = LROW
      KEEP(5)  = LCOL
      RETURN
      END
      SUBROUTINE LA15ED(A,IRN,IP,N,IW,IA,REALS,NCP,LROW,LCOL)
C     .. Scalar Arguments ..
      INTEGER IA,N,NCP,LROW,LCOL
      LOGICAL REALS
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(IA)
      INTEGER IP(N),IRN(IA),IW(N)
C     ..
C     .. Local Scalars ..
      INTEGER IPI,J,K,KL,KN,NZ
C     ..
C     .. Executable Statements ..
      NCP = NCP + 1
C     COMPRESS FILE OF POSITIVE INTEGERS. ENTRY J STARTS AT IRN(IP(J))
C  AND CONTAINS IW(J) INTEGERS,J=1,N. OTHER COMPONENTS OF IRN ARE ZERO.
C  LENGTH OF COMPRESSED FILE PLACED IN LROW IF REALS IS .TRUE. OR LCOL
C  OTHERWISE.
C  IF REALS IS .TRUE. ARRAY A CONTAINS A REAL FILE ASSOCIATED WITH IRN
C  AND THIS IS COMPRESSED TOO.
C  A,IRN,IP,IW,IA ARE INPUT/OUTPUT VARIABLES.
C  N,REALS ARE INPUT/UNCHANGED VARIABLES.
C
      DO 10 J = 1,N
C STORE THE LAST ELEMENT OF ENTRY J IN IW(J) THEN OVERWRITE IT BY -J.
        NZ = IW(J)
        IF (NZ.LE.0) GO TO 10
        K = IP(J) + NZ - 1
        IW(J) = IRN(K)
        IRN(K) = -J
   10 CONTINUE
C KN IS THE POSITION OF NEXT ENTRY IN COMPRESSED FILE.
      KN = 0
      IPI = 0
      KL = LCOL
      IF (REALS) KL = LROW
C LOOP THROUGH THE OLD FILE SKIPPING ZERO (DUMMY) ELEMENTS AND
C     MOVING GENUINE ELEMENTS FORWARD. THE ENTRY NUMBER BECOMES
C     KNOWN ONLY WHEN ITS END IS DETECTED BY THE PRESENCE OF A NEGATIVE
C     INTEGER.
      DO 30 K = 1,KL
        IF (IRN(K).EQ.0) GO TO 30
        KN = KN + 1
        IF (REALS) A(KN) = A(K)
        IF (IRN(K).GE.0) GO TO 20
C END OF ENTRY. RESTORE IRN(K), SET POINTER TO START OF ENTRY AND
C     STORE CURRENT KN IN IPI READY FOR USE WHEN NEXT LAST ENTRY
C     IS DETECTED.
        J = -IRN(K)
        IRN(K) = IW(J)
        IP(J) = IPI + 1
        IW(J) = KN - IPI
        IPI = KN
   20   IRN(KN) = IRN(K)
   30 CONTINUE
      IF (REALS) LROW = KN
      IF (.NOT.REALS) LCOL = KN
      RETURN
      END

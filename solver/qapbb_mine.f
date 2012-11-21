      PROGRAM QAPBB
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IMPLICIT INTEGER (A - Z)
      parameter (ndim = 33)
c auslegung fuer ndim = 33
      DIMENSION A(ndim,ndim),B(ndim,ndim),C(ndim,ndim),ZUL(ndim,ndim)
      dimension UMSPEI(ndim*ndim),LOESG(ndim),U(ndim), V(ndim)
     *          , DD(ndim), PARTPE(ndim), Y(ndim),LAB(ndim),Z1(ndim)
     *          ,H1(ndim),MENG(ndim-2),PHIOFM(ndim-2),MENGE(ndim-2),
     *          VEKSUM(ndim-2), VEKQU(ndim-2), ALTER(ndim-2)
      dimension ASPEI(11968),BSPEI(11968), CSPEI(12528)
      LOGICAL BOOL(ndim), BOOL1(ndim), HL1(ndim)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      UNENDL = 1000000000
      call ilese(n, a, b, c, ndim)
c      print *,'0=no feasible sol., else input known value'
c      read *,olwert
        olwert = 0
        if (olwert .eq. 0) olwert = unendl
c      WRITE(6,1000) N
c      DO 12 I=1,N
c 12   WRITE(6,1010) (A(I,J),J=1,N)
c      WRITE(6,1015)
c      DO 13 I=1,N
c 13   WRITE(6,1010) (B(I,J),J=1,N)
c      write(6,1016)
c      do 14 i=1,n
c14    write(6,1010) (c(i,j),j=1,n)
c        write( 6,1011) olwert
      CALL QAP (N,A,B,UNENDL,LOESG,OLWERT,C,UMSPEI,ZUL,U,V,DD,
     *          PARTPE,Y,LAB,Z1,H1,MENG,PHIOFM,MENGE,VEKSUM,VEKQU,
     *          ALTER,ASPEI,BSPEI,CSPEI,BOOL,BOOL1,HL1,NDIM)
      WRITE(6,2000)
      DO 200 I=1,N
 200  WRITE(6,2005) I, LOESG(I)
      WRITE(6,2010) OLWERT
      print *,OLWERT
      print *,(loesg(i),i=1,n)      
 102  FORMAT(I2)
 103  FORMAT(20I3)
 1000 FORMAT(20X,'N =',I3,///, 15X,17HDISTANCE MATRIX A,/)
 1010 FORMAT(6X,20I4)
 1011 format(6x,'a feasible solution with value=',i10,4x,'was supplied')
 1015 FORMAT(//15X,19HCONNECTION MATRIX B,/)
 1016 format(//15x,'LINEAR TERM:',/)
 2000 FORMAT(///9X,15HEXACT SOLUTION:)
 2005 FORMAT(9X,I5,4H -->,I2)
 2010 FORMAT(/9X,25HOBJECTIVE FUNCTION VALUE:,I6)

      END

      SUBROUTINE QAP (N,A,B,UNENDL,LOESG,OLWERT,C,UMSPEI,ZUL,U,V,DD,
     *                PARTPE,Y,LAB,Z1,H1,MENG,PHIOFM,MENGE,VEKSUM,
     *         VEKQU,ALTER,ASPEI,BSPEI,CSPEI,BOOL,BOOL1,HL1,NDIM)
C *** ****************************************************************
C     *                                                              *
C     *          PROCEDURE FOR SOLVING QUADRATIC ASSIGNMENT          *
C     *                          PROBLEMS                            *
C     *                                                              *
C *** ****************************************************************
C     *                                                              *
C     *  1. CALL:                                                    *
C     *     CALL QAP (N,A,B,UNENDL,LOESG,OLWERT,C,UMSPEI,ZUL,U,V,DD, *
C     *               PARTPE,Y,LAB,Z1,H1,MENG,PHIOFM,MENGE,VEKSUM,   *
C     *               VEKQU,ALTER,ASPEI,BSPEI,CSPEI,BOOL,BOOL1,HL1)  *
C     *                                                              *
C     *  2. COMPUTER CODE:                                           *
C     *     FORTRAN IV                                               *
C     *                                                              *
C     *  3. METHOD:                                                  *
C     *     PERTURBATION METHOD FOR THE OPTIMAL SOLUTION OF          *
C     *     QUADRATIC ASSIGNMENT PROBLEMS                            *
C     *                                                              *
C     *  4. PARAMETERS:                                              *
C     *     INPUT:                                                   *
C     *        N        DIMENSION OF THE PROBLEM                     *
C     *        NDIM     DIMENSION OF CALLING PROGRAM
C     *        A(I,J)   DISTANCE MATRIX (INTEGER)  I,J=1,...,N       *
C     *        B(I,J)   CONNECTION MATRIX (INTEGER)  I,J=1,...,N     *
C     *        UNENDL   LARGE MACHINE NUMBER (INTEGER)               *
C     *     OUTPUT:                                                  *
C     *        LOESG(I) OPTIMAL ASSIGNMENT (INTEGER)                 *
C     *                  I --> LOESG(I)  I=1,...,N                   *
C     *        OLWERT   OPTIMAL OBJECTIVE FUNCTION VALUE (INTEGER)   *
C     *     INTEGER ARRAYS OF DIMENSION (N,N)                        *
C     *        C(I,J), ZUL(I,J)                                      *
C     *     INTEGER ARRAYS OF DIMENSION N                            *
C     *        U(I), V(I), DD(I), PARTPE(I), Y(I), LAB(I), Z1(I),    *
C     *        H1(I)                                                 *
C     *     INTEGER ARRAY OF DIMENSION N*N                           *
C     *        UMSPEI(I)                                             *
C     *     INTEGER ARRAYS OF DIMENSION N-2                          *
C     *        MENG(I), PHIOFM(I), MENGE(I), VEKSUM(I), VEKQU(I),    *
C     *        ALTER(I)                                              *
C     *     INTEGER ARRAYS OF DIMENSION N*(N+1)*(2*N-2)/6            *
C     *        ASPEI(I), BSPEI(I)                                    *
C     *     INTEGER ARRAY OF DIMENSION (N*(N+1)*(2*N+1)/6)-1         *
C     *        CSPEI(I)                                              *
C     *     LOGICAL ARRAYS OF DIMENSION N                            *
C     *        BOOL(I), BOOL1(I), HL1(I)                             *
C     *                                                              *
C     *  5. EXTERNAL SUBROUTINES:                                    *
C     *     SUBROUTINE ALTKOS                                        *
C     *     SUBROUTINE WEGSPE                                        *
C     *     SUBROUTINE PROGNO                                        *
C     *     SUBROUTINE LSAP                                          *
C     *     SUBROUTINE SSORT                                         *
C     *                                                              *
C<<
C     *  6. COMMENT:                                                 *
C     *     THE ELEMENTS OF THE MATRICES A AND B HAVE TO BE          *
C     *     NONNEGATIVE INTEGERS. THEY ARE CHANGED DURING THE        *
C     *     PROCEDURE.                                               *
C     *                                                              *
C     *  7. AUTHORS:                                                 *
C     *     T. BOENNIGER                                             *
C     *     R.E. BURKARD                                             *
C     *     K.-H. STRATMANN                                          *
C     *                                                              *
C *** ****************************************************************
C
C<<
      IMPLICIT INTEGER (A - Z)
      INTEGER A(NDIM,1), B(NDIM,1), C(NDIM,1), UMSPEI(1), U(1), V(1),
     *        DD(1), ZUL(NDIM,1), MENG(1), PHIOFM(1), PARTPE(1),
     *        LOESG(1), Y(1), LAB(1), Z1(1), MENGE(1), ASPEI(1),
     *        BSPEI(1), CSPEI(1), VEKSUM(1), VEKQU(1), ALTER(1),
     *        H1(1)
      LOGICAL BOOL(1), BOOL1(1), HL1(1), WEITER
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C *** INITIALIZING OF PARAMETERS
C
      K = 0
c ***      OLWERT = UNENDL (here a user defined value can come in)
      NM2 = N-2
      I2 = N + 1
      J = 0
      J1 = 0
      DO 5035 I = 1,NM2
      NMK = I2 - I
      J2 = NMK*NMK
      J = J + J2 - NMK
      VEKSUM(I) = J
      J1 = J1 + J2
 5035 VEKQU(I) = J1
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C *** 1. REDUCTION: C(I,J) = A(I,I) * B(J,J)
C
      DO 5010 I = 1,N
      BOOL(I) = .FALSE.
      BOOL1(I) = .FALSE.
      ZSTERN = A(I,I)
      DO 5010 J = 1,N
      ZUL(I,J) = -1
c consider linear term c
 5010 C(I,J) = ZSTERN * B(J,J) + c( i, j)
      DO 5030 J = 1,N
      A(J,J) = UNENDL
 5030 B(J,J) = UNENDL
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C *** REDUCTION OF MATRICES A AND B
C
      DO 5103 I=1,N
      RA = A(I,1)
      RB = B(I,1)
      DO 5101 J=2,N
      RAA = A(I,J)
      RBB = B(I,J)
      IF(RAA .LT. RA) RA = RAA
      IF(RBB .LT. RB) RB = RBB
 5101 CONTINUE
      DO 5102 J=1,N
      A(I,J) = A(I,J)-RA
 5102 B(I,J) = B(I,J)-RB
      U(I) = RA
 5103 V(I) = RB
C<<
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C *** DETERMINATION OF THE MATRIX C
C     ROWWISE REDUCTION OF MATRICES A AND B
C
      DO 5106 I=1,N
      CA1 = U(I)
      DO 5106 J=1,N
      BMJ = V(J)
      BM = (N-1)*BMJ
      DO 5104 KK=1,N
      IF(KK .EQ. J) GOTO 5104
      BM = BM+B(J,KK)
 5104 CONTINUE
      CA = CA1*BM
      AM = 0
      DO 5105 KK=1,N
      IF(KK .EQ. I) GOTO 5105
      AM = AM+A(I,KK)
 5105 CONTINUE
 5106 C(I,J) = CA+BMJ*AM+C(I,J)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C *** COLUMNWISE REDUCTION OF MATRICES A AND B
C
      DO 5109 I=1,N
      RA = A(1,I)
      RB = B(1,I)
      DO 5107 J=2,N
      RAA = A(J,I)
      RBB = B(J,I)
      IF(RAA .LT. RA) RA = RAA
      IF(RBB .LT. RB) RB = RBB
 5107 CONTINUE
      DO 5108 J=1,N
      A(J,I) = A(J,I)-RA
 5108 B(J,I) = B(J,I)-RB
      U(I) = RA
 5109 V(I) = RB
      DO 5112 I=1,N
      A(I,I) = 0
      B(I,I) = 0
      CA1 = U(I)
      DO 5112 J=1,N
      BMJ = V(J)
      BM = (N-1)*BMJ
      DO 5110 KK=1,N
      IF(KK .EQ. J) GOTO 5110
      BM = BM+B(KK,J)
 5110 CONTINUE
      CA = CA1*BM
      AM = 0
      DO 5111 KK=1,N
      IF(KK .EQ. I) GOTO 5111
      AM = AM+A(KK,I)
 5111 CONTINUE
      CCC = C(I,J)+CA+BMJ*AM
      C(I,J) = CCC
 5112 CONTINUE
      ZPART = 0
      NMK = N
      WEITER = .TRUE.
C<<
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C *** IMPROVEMENT OF THE BOUNDS C(I,J) BY ADDING MINIMAL SCALAR
C     PRODUCTS (GILMORE BOUNDS). THE SUBROUTINES WEGSPE AND PROGNO
C     COMPUTE THESE MINIMAL SCALAR PRODUCTS.
C
      CALL WEGSPE (N,K,NMK,IZAEHL,JZAEHL,A,B,VEKSUM,U,BOOL,BOOL1,
     *             ASPEI,BSPEI,H1,NDIM)
      CALL PROGNO (N,K,NMK,UMSPEI,VEKSUM,VEKQU,WEITER,ASPEI,BSPEI,
     *             CSPEI)
      WEITER = .FALSE.
 5210 IZ = 0
      I2 = 0
      J2 = 0
      IKAP = 0
      DO 5120 I = 1,N
      IF (BOOL(I)) GOTO 5120
      JJ = IZ*NMK
      IZ = IZ + 1
      JZ = 0
      DO 5125 J = 1,N
      IF (BOOL1(J)) GOTO 5125
      JZ = JZ + 1
      JJ = JJ+1
      IF (ZUL(I,J) .LT. 0) GOTO 5145
      UMSPEI(JJ) = UNENDL
      IKAP = IKAP + 1
      IF (IKAP - 2) 5375, 5380, 5125
 5375 IHALT = IZ
      JHALT = JZ
      GOTO 5125
 5380 IF (IZ .EQ. IHALT) I2 = IZ
      IF (JZ .EQ. JHALT) J2 = JZ
      GOTO 5125
 5145 UMSPEI(JJ) = C(I,J) + UMSPEI(JJ)
 5125 CONTINUE
 5120 CONTINUE
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C *** COMPUTATION OF A BOUND BY SOLVING LSAP WITH COSTMATRIX UMSPEI
C
 5180 CALL LSAP (NMK,UNENDL,UMSPEI,ZSTERN,Y,Z1,LAB,DD,U,V,H1,HL1)
      JJ = 0
      DO 5182 I=1,NMK
      DO 5181 J=1,NMK
      JJ = JJ+1
      UMSPEI(JJ) = UMSPEI(JJ)-U(I)-V(J)
 5181 CONTINUE
 5182 CONTINUE
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C *** ZPART IS THE FIXED FRACTION OF THE OBJECTIVE FUNCTION VALUE
C     IMPLIED BY THE PRESENT PARTIAL PERMUTATION.
C     THE PRESENT BOUND IS ZPART+ZSTERN.
C
      IF (ZPART+ZSTERN.GE.OLWERT) GOTO 5250
      if (lbd .lt. zpart+zstern) then
                lbd = zpart + zstern
c                print 998, lbd
                        endif
998     format(3x,'lb=',i12)
      IF (WEITER) GOTO 5220
 5135 IF (IKAP - 1) 5400, 5410, 5420
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C *** COMPUTATION OF THE ALTERNATIVE COSTS
C
 5400 CALL ALTKOS (NMK,UMSPEI,Z1,UNENDL,IZAEHL,JZAEHL,ALTERK)
      GOTO 5490
C<<
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C *** COMPUTATION OF THE NEXT SINGLE ASSIGNMENT IN A FIXED ROW
C     OR COLUMN
C
 5410 MIN = UNENDL
      J1 = Z1(IHALT)
      JJ = (IHALT-1)*NMK
      DO 5430 J = 1,NMK
      JJ = JJ+1
      T = UMSPEI(JJ)
      IF (MIN .LE. T .OR. J .EQ. J1)  GOTO 5430
      MIN = T
 5430 CONTINUE
      MIN1 = MIN
      MIN = UNENDL
      JJ = J1
      DO 5510 I = 1,NMK
      T = UMSPEI(JJ)
      JJ = JJ+NMK
      IF (MIN .LE. T .OR. I .EQ. IHALT) GOTO 5510
      MIN = T
 5510 CONTINUE
      MIN1 = MIN1 + MIN
      MIN = UNENDL
      JJ = JHALT
      DO 5440 I = 1,NMK
      T = UMSPEI(JJ)
      JJ = JJ+NMK
      IF (T .GE. MIN .OR. JHALT .EQ.Z1(I)) GOTO 5440
      MIN = T
 5440 CONTINUE
      MIN2 = MIN
      DO 5530 I = 1,NMK
      IF (Z1(I) .EQ. JHALT) GOTO 5540
 5530 CONTINUE
 5540 I1 = I
      JJ = (I1-1)*NMK
      MIN = UNENDL
      DO 5520 J = 1,NMK
      JJ = JJ+1
      T = UMSPEI(JJ)
      IF (T .GE. MIN .OR. J .EQ. JHALT) GOTO 5520
      MIN = T
 5520 CONTINUE
      IF ((MIN+MIN2) .LT. MIN1) GOTO 5450
      IZAEHL = I1
      JZAEHL = JHALT
      ALTERK = MIN + MIN2
      GOTO 5490
 5450 IZAEHL = IHALT
      JZAEHL = J1
      ALTERK = MIN1
      GOTO 5490
C<<
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C *** COMPUTATION OF THE NEXT SINGLE ASSIGNMENT IN THE PREVIOUSLY
C     FIXED ROW (RESP. PREVIOUSLY FIXED COLUMN)
C
 5420 IF (I2 .EQ. 0) GOTO 5480
      IZAEHL = I2
      JZAEHL = Z1(IZAEHL)
      GOTO 5370
 5480 JZAEHL = J2
      DO 5470 I = 1,NMK
      IF (Z1(I) .EQ. JZAEHL) GOTO 5485
 5470 CONTINUE
 5485 IZAEHL = I
 5370 MIN = UNENDL
      JJ = (IZAEHL-1)*NMK
      DO 5350 I = 1,NMK
      JJ = JJ+1
      T = UMSPEI(JJ)
      IF (T .GE. MIN .OR. I .EQ. JZAEHL) GOTO 5350
      MIN = T
 5350 CONTINUE
      ALTERK = MIN
      MIN = UNENDL
      JJ = JZAEHL
      DO 5360 J = 1,NMK
      T = UMSPEI(JJ)
      JJ = JJ+NMK
      IF (T .GE. MIN .OR. J .EQ. IZAEHL) GOTO 5360
      MIN = T
 5360 CONTINUE
      ALTERK = ALTERK + MIN
 5490 IZ = 0
      ALTER(K+1) = ALTERK + ZSTERN
      DO 5155 I = 1,N
      IF (BOOL(I)) GOTO 5155
      IZ = IZ + 1
      IF (IZAEHL .EQ. IZ) GOTO 5165
 5155 CONTINUE
 5165 IZAEHL = I
      IZ = 0
      DO 5175 J = 1,N
      IF (BOOL1(J)) GOTO 5175
      IZ = IZ + 1
      IF (JZAEHL .EQ. IZ) GOTO 5185
 5175 CONTINUE
 5185 JZAEHL = J
      ZUL(IZAEHL,JZAEHL) =K
      WEITER = .TRUE.
      BOOL(IZAEHL) = .TRUE.
      BOOL1(JZAEHL) = .TRUE.
      NMK = N - K - 1
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C *** COMPUTATION OF THE COST MATRIX CORRESPONDING TO THE NEW PARTIAL
C     PERMUTATION
C
      CALL WEGSPE (N,K,NMK,IZAEHL,JZAEHL,A,B,VEKSUM,U,BOOL,BOOL1,
     *             ASPEI,BSPEI,H1,NDIM)
      CALL PROGNO (N,K,NMK,UMSPEI,VEKSUM,VEKQU,WEITER,ASPEI,BSPEI,
     *             CSPEI)
C<<
      IZ = 0
      DO 5130 I = 1,N
      IF (BOOL(I)) GOTO 5130
      T1 = A(I,IZAEHL)
      T2 = A(IZAEHL,I)
      JZ = IZ
      DO 5140 J = 1,N
      IF (BOOL1(J)) GOTO 5140
      JZ = JZ+1
      UMSPEI(JZ) = C(I,J)+T1*B(J,JZAEHL)+T2*B(JZAEHL,J)+UMSPEI(JZ)
 5140 CONTINUE
      IZ = IZ+NMK
 5130 CONTINUE
      ZPART = ZPART + C(IZAEHL,JZAEHL)
      GOTO 5180
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C *** THE BOUND FOR THE NEW PARTIAL PERMUTATION IS NOT LESS THAN A
C     PREVIOUSLY FOUND OBJECTIVE FUNCTION VALUE.
C     - BACKTRACKING -
C
 5250 IF (.NOT. WEITER) GOTO 5255
      WEITER = .FALSE.
      K=K+1
      GOTO 5230
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C *** EXIT: THE SOLUTION TREE IS COMPLETELY FATHOMED.
C
 5255 IF (K. EQ. 0) GOTO 5600
      IZAEHL = MENGE(K)
      JZAEHL = PHIOFM(K)
 5220 DO 5150 I = 1,N
      IF (BOOL(I)) GOTO 5150
      T1 = A(IZAEHL,I)
      T2 = A(I,IZAEHL)
      DO 5160 J = 1,N
      IF (BOOL1(J)) GOTO 5160
      T = T1 * B(JZAEHL,J) + T2 * B(J,JZAEHL)
      IF (.NOT. WEITER) T = -T
      C(I,J) = C(I,J) + T
 5160 CONTINUE
 5150 CONTINUE
      IF (.NOT. WEITER) GOTO 5230
      PARTPE(IZAEHL) = JZAEHL
      K = K + 1
      MENGE(K) = IZAEHL
      PHIOFM(K) = JZAEHL
      IF (K .EQ. (N-2)) GOTO 5270
      IKAP = 0
      GOTO 5135
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C *** CANCELLATION OF THE LAST SINGLE ASSIGNMENT
C
 5230 DO 5260 I = 1,N
      IF (BOOL(I)) GOTO 5260
      DO 5240 J = 1,N
      IF (.NOT. BOOL1(J) .AND. ZUL(I,J) .EQ. K) ZUL(I,J) = -1
 5240 CONTINUE
 5260 CONTINUE
C<<
      ZPART = ZPART - C(IZAEHL,JZAEHL)
      BOOL(IZAEHL) = .FALSE.
      BOOL1(JZAEHL) = .FALSE.
      K = K - 1
      NMK = N - K
      IF (ALTER(K+1)+ZPART.GE.OLWERT) GOTO 5255
      CALL PROGNO (N,K,NMK,UMSPEI,VEKSUM,VEKQU,WEITER,ASPEI,BSPEI,
     *             CSPEI)
      GOTO 5210
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C *** COMPUTATION OF THE OBJECTIVE FUNCTION VALUES FOR THE REMAINING
C     TWO COMPLETE PERMUTATIONS
C
 5270 DO 5290 I = 1,N
      IF (BOOL(I)) GOTO 5290
      IZ = I
      GOTO 5285
 5290 CONTINUE
 5285 DO 5280 I = 1,N
      IF (BOOL1(I)) GOTO 5280
      J = I
      GOTO 5295
 5280 CONTINUE
 5295 BOOL(IZ) = .TRUE.
      BOOL1(J) = .TRUE.
      DO 5310 I = 1,N
      IF (BOOL(I)) GOTO 5310
      I1 = I
      GOTO 5305
 5310 CONTINUE
 5305 DO 5300 I = 1,N
      IF(BOOL1(I)) GOTO 5300
      J1 = I
      GOTO 5315
 5300 CONTINUE
 5315 WEITER = .FALSE.
      IZE = 0
 5330 ZSTERN = C(IZ,J) + C(I1,J1) + A(IZ,I1)*B(J,J1) + A(I1,IZ)*B(J1,J)
      BOOL(IZ) = .FALSE.
      BOOL1(J) = .FALSE.
      IF ((ZSTERN+ZPART) .GE. OLWERT) GOTO 5340
      OLWERT=ZSTERN+ZPART
c      print 999,olwert
      lbd = 0
999     format(20x,'new solution=',i12)
      DO 5320 I=1,N
      IF (BOOL(I)) LOESG(I) = PARTPE(I)
 5320 CONTINUE
      DO 5335 I = 1,NM2
 5335 MENG(I) = MENGE(I)
      LOESG(IZ) = J
      LOESG(I1) = J1
 5340 IF(IZE.NE.0) GOTO 5255
      IZE = IZ
      IZ = I1
      I1 = IZE
      GOTO 5330
 5600 CONTINUE
      RETURN
      END
C<<
      SUBROUTINE ALTKOS (N,UMSPEI,Z1,UNENDL,IZAEHL,JZAEHL,ALTERK)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      INTEGER UMSPEI(1), Z1(1)
      INTEGER A, ALTERK, UNENDL
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C *** DETERMINATION OF ALTERNATIVE COSTS AND THE RESULTING SINGLE
C     ASSIGNMENT.
C
      ALTERK = -1
      JJ = 0
      DO 4110 I = 1,N
      J = Z1(I)
      J1 = J-N
      MIN = UNENDL
      DO 4120 I1 = 1,N
      J1 = J1+N
      IF (I1 .EQ. I) GOTO 4120
      A = UMSPEI(J1)
      IF (A .LT. MIN) MIN = A
 4120 CONTINUE
      MIN1 = MIN
      MIN = UNENDL
      DO 4130 I1 = 1,N
      JJ = JJ+1
      IF(I1 .EQ. J) GOTO 4130
      A = UMSPEI(JJ)
      IF (A .LT. MIN) MIN = A
 4130 CONTINUE
      MIN = MIN + MIN1
      IF (MIN .LE. ALTERK) GOTO 4110
      IZAEHL = I
      JZAEHL = J
      ALTERK = MIN
 4110 CONTINUE
      RETURN
      END
C<<
      SUBROUTINE WEGSPE (N,K,NMK,IZAEHL,JZAEHL,A,B,VEKSUM,VEKT,BOOL,
     *                   BOOL1,ASPEI,BSPEI,H1,NDIM)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      INTEGER A(NDIM,1), B(NDIM,1), VEKT(1), VEKSUM(1), ASPEI(1)
      INTEGER BSPEI(1), H1(1), T
      LOGICAL BOOL(1), BOOL1(1), LOGI
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C *** ORDERING OF THE ELEMENTS IN THE ROWS OF THE MATRIX A
C     DECREASINGLY AND IN MATRIX B INCREASINGLY. THE RESULTING
C     MATRICES ARE ROWWISE STORED ON THE VECTORS ASPEI AND BSPEI.
C
      NMKM2 = NMK - 1
      IF (NMK .NE. N) GOTO 2900
      J1 = 1
      DO 2810 I = 1,N
      IZ = 0
      DO 2820 J = 1,N
      IF (J .EQ. I) GOTO 2820
      IZ = IZ + 1
      VEKT(IZ) = A(I,J)
 2820 CONTINUE
      CALL SSORT (VEKT,H1,NMKM2)
      DO 2810 J = 1,NMKM2
      NMKJ = NMK-J
      ASPEI(J1) = VEKT(NMKJ)
 2810 J1 = J1 + 1
      J1 = 1
      DO 2840 I = 1,N
      IZ = 0
      DO 2850 J = 1,N
      IF (J .EQ. I) GOTO 2850
      IZ = IZ + 1
      VEKT(IZ) = B(I,J)
 2850 CONTINUE
      CALL SSORT (VEKT,H1,NMKM2)
      DO 2840 J = 1,NMKM2
      BSPEI(J1) = VEKT(J)
 2840 J1 = J1 + 1
      RETURN
 2900 NSUM1 = VEKSUM(K+1)
      J1 = NSUM1
      J2 = NSUM1 - (NMK+1) * NMK
      DO 2910 I = 1,N
      IF (BOOL(I)) GOTO 2930
      T = A(I,IZAEHL)
      LOGI = .TRUE.
      DO 2920 J = 1,NMK
      J2 = J2 + 1
      IF (ASPEI(J2) .EQ. T .AND. LOGI) GOTO 2980
      J1 = J1 + 1
      ASPEI(J1) = ASPEI(J2)
      GOTO 2920
 2980 LOGI = .FALSE.
 2920 CONTINUE
      GOTO 2910
C<<
 2930 IF (I .EQ. IZAEHL) J2 = J2 + NMK
 2910 CONTINUE
      J1 = NSUM1
      J2 = NSUM1 - (NMK+1) * NMK
      DO 2940 I = 1,N
      IF (BOOL1(I)) GOTO 2960
      T = B(I,JZAEHL)
      LOGI = .TRUE.
      DO 2950 J = 1,NMK
      J2 = J2 + 1
      IF (BSPEI(J2) .EQ. T .AND. LOGI) GOTO 2970
      J1 = J1 + 1
      BSPEI(J1) = BSPEI(J2)
      GOTO 2950
 2970 LOGI = .FALSE.
 2950 CONTINUE
      GOTO 2940
 2960 IF (I .EQ. JZAEHL) J2 = J2 + NMK
 2940 CONTINUE
      RETURN
      END
C<<
      SUBROUTINE PROGNO (N,K,NMK,UMSPEI,VEKSUM,VEKQU,WEITER,ASPEI,
     *                   BSPEI,CSPEI)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      INTEGER UMSPEI(1), VEKSUM(1), VEKQU(1), ASPEI(1), BSPEI(1)
      INTEGER CSPEI(1), T
      LOGICAL WEITER
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C *** COMPUTATION OF THE FINAL COST MATRIX UMSPEI FOR THE LINEAR
C     SUBPROBLEM LSAP ON THE K-TH STAGE OF THE DECISION TREE
C
      JU = 0
      IF (NMK .NE. N) GOTO 3600
      NSUM2 = 0
      J1 = 0
      IF (WEITER) GOTO 3680
      GOTO 3690
 3600 IF (WEITER) GOTO 3630
      J1 = VEKQU(K)
 3690 DO 3640 I = 1,NMK
      DO 3640 J = 1,NMK
      J1 = J1 + 1
      JU = JU + 1
 3640 UMSPEI(JU) = CSPEI(J1)
      RETURN
 3630 NSUM2 = VEKSUM(K+1)
      J1 = VEKQU(K+1)
 3680 J2 = NSUM2
      NMKM1 = NMK - 1
      DO 3650 I = 1,NMK
      J3 = NSUM2
      DO 3660 J = 1,NMK
      J1 = J1 + 1
      JU = JU + 1
      T = 0
      DO 3670 IZ = 1,NMKM1
      J2IZ = J2+IZ
      J3IZ = J3+IZ
 3670 T = T+ASPEI(J2IZ)*BSPEI(J3IZ)
      J3 = J3 + NMKM1
      CSPEI(J1) = T
 3660 UMSPEI(JU) = T
 3650 J2 = J2 + NMKM1
      RETURN
      END

      SUBROUTINE LSAP(N,SUP,C,Z,ZEILE,SPALTE,DMINUS,DPLUS,YS,YT,VOR,
     1                LABEL)
C *** ****************************************************************
C     *                                                              *
C     *                                                              *
C     *            LINEAR SUM ASSIGNMENT PROBLEM                     *
C     *                                                              *
C *** ****************************************************************
C     *  1.  CALL:                                                   *
C     *      CALL LSAP(N,SUP,C,Z,ZEILE,SPALTE,DMINUS,DPLUS,YS,YT,VOR,*
C     *                LABEL)                                        *
C     *                                                              *
C     *  2.  COMPUTER CODE:                                          *
C     *      FORTRAN IV                                              *
C     *                                                              *
C     *  3.  METHOD:                                                 *
C     *      SHORTEST AUGMENTING PATH METHOD                         *
C     *                                                              *
C     *  4.  PARAMETERS:                                             *
C     *      INPUT:                                                  *
C     *         N     DIMENSION OF THE COST MATRIX  C                *
C     *         C(I)  COST MATRIX (I=1,...N*N) , ROWWISE STORED      *
C     *         SUP   LARGE MACHINE NUMBER                           *
C     *      OUTPUT:                                                 *
C     *         SPALTE(I)   OPTIMAL ASSIGNMENT (I=1,...,N)           *
C     *         Z           OPTIMAL VALUE                            *
C     *         YS(I)     OPTIMAL DUAL (ROW) VARIABLES    (I=1,...,N)*
C     *         YT(I)     OPTIMAL DUAL (COLUMN) VARIABLES (I=1,...,N)*
C     *      INTEGER ARRAYS OF LENGHTS N:                            *
C     *         ZEILE(I),DMINUS(I),DPLUS(I),YS(I),YT(I),VOR(1)       *
C     *      LOGICAL ARRAY OF LENGTH N:                              *
C     *         LABEL(I)                                             *
C     *                                                              *
C     *  5.  EXTERNAL SUBROUTINES:                                   *
C     *      NONE                                                    *
C     *                                                              *
C     *  6.  AUTHOR:                                                 *
C     *      U.DERIGS                                                *
C     *                                                              *
C *** ****************************************************************
C
C<<
      IMPLICIT INTEGER (A-Z)
      DIMENSION YS(1),YT(1),VOR(1)
      DIMENSION C(1),ZEILE(1),SPALTE(1),DMINUS(1),DPLUS(1)
      LOGICAL LABEL(1)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C *** STARTPROCEDURE
C     CONSTRUCTION OF AN INITIAL (PARTIAL) ASSIGNMENT
C
      DO 50 I=1,N
      ZEILE(I)=0
      SPALTE(I)=0
      VOR(I)=0
      YS(I)=0
   50 YT(I)=0
      IK=0
      DO 2 I=1,N
      DO 3 J=1,N
      IK=IK+1
      CC=C(IK)
      IF(J.EQ.1) GOTO 4
      IF((CC-UI).GE.0) GOTO 3
    4 UI=CC
      JO=J
    3 CONTINUE
      YS(I)=UI
      IF(ZEILE(JO).NE.0) GOTO 2
      ZEILE(JO)=I
      SPALTE(I)=JO
    2 CONTINUE
      DO 5 J=1,N
      YT(J)=0
      IF(ZEILE(J).EQ.0) YT(J)=SUP
    5 CONTINUE
      IK=0
      DO 6 I=1,N
      UI=YS(I)
      DO 7 J=1,N
      IK=IK+1
      VJ=YT(J)
      IF(VJ.LE.0) GOTO 7
      CC=C(IK)-UI
      IF(CC.GE.VJ) GOTO 7
      YT(J)=CC
      VOR(J)=I
    7 CONTINUE
    6 CONTINUE
      DO 8 J=1,N
      I=VOR(J)
      IF(I.EQ.0) GOTO 8
      IF(SPALTE(I).NE.0) GOTO 8
      SPALTE(I)=J
      ZEILE(J)=I
    8 CONTINUE
C<<
      DO 9 I=1,N
      IF(SPALTE(I).NE.0) GOTO 9
      UI=YS(I)
      IK=(I-1)*N
      DO 10 J=1,N
      IK=IK+1
      IF(ZEILE(J).NE.0) GOTO 10
      CC=C(IK)
      IF((CC-UI-YT(J)).GT.0) GOTO 10
      SPALTE(I)=J
      ZEILE(J)=I
      GOTO 9
   10 CONTINUE
    9 CONTINUE
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C *** CONSTRUCTION OF THE OPTIMAL ASSIGNMENT
C
      DO 1000 U=1,N
      IF(SPALTE(U).GT.0) GOTO 1000
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C *** SHORTEST PATH COMPUTATION
C
      US=(U-1)*N
      DO 100 I=1,N
      VOR(I)=U
      LABEL(I)=.FALSE.
      DPLUS(I)=SUP
      USI=US+I
  100 DMINUS(I)=C(USI)-YS(U)-YT(I)
      DPLUS(U)=0
  105 D=SUP
      DO 110 I=1,N
      IF(LABEL(I)) GOTO 110
      IF(DMINUS(I).GE.D) GOTO 110
      D=DMINUS(I)
      INDEX=I
  110 CONTINUE
      IF(ZEILE(INDEX).LE.0) GOTO 400
      LABEL(INDEX)=.TRUE.
      W=ZEILE(INDEX)
      WS=(W-1)*N
      DPLUS(W)=D
      DO 130 I=1,N
      IF(LABEL(I)) GOTO 130
      WSI=WS+I
      VGL=D+C(WSI)-YS(W)-YT(I)
      IF(DMINUS(I).LE.VGL) GOTO 130
      DMINUS(I)=VGL
      VOR(I)=W
  130 CONTINUE
      GOTO 105
C<<
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C *** AUGMENTATION
C
  400 W=VOR(INDEX)
      ZEILE(INDEX)=W
      IND=SPALTE(W)
      SPALTE(W)=INDEX
      IF(W.EQ.U) GOTO 500
      INDEX=IND
      GOTO 400
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C *** TRANSFORMATION
C
  500 DO 510 I=1,N
      IF(DPLUS(I).EQ.SUP) GOTO 505
      YS(I)=YS(I)+D-DPLUS(I)
  505 IF(DMINUS(I).GE.D) GOTO 510
      YT(I)=YT(I)+DMINUS(I)-D
  510 CONTINUE
 1000 CONTINUE
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C *** COMPUTATION OF THE OPTIMAL VALUE
C
      Z=0
      DO 2000 I=1,N
      IS=(I-1)*N
      J=SPALTE(I)
      ISJ=IS+J
      Z=Z+C(ISJ)
 2000 CONTINUE
      RETURN
      END

      SUBROUTINE SSORT (A,B,L)
C *** ****************************************************************
C     *                                                              *
C     *                       SORT - PROGRAM                         *
C     *                                                              *
C *** ****************************************************************
C     *                                                              *
C     *  1. CALL:                                                    *
C     *     CALL SSORT (A,B,N)                                       *
C     *                                                              *
C     *  2. COMPUTER CODE:                                           *
C     *     FORTRAN IV                                               *
C     *                                                              *
C     *  3. METHOD:                                                  *
C     *     SORTING OF THE VECTOR A(I) IN INCREASING ORDER BY THE    *
C     *     SHELL-ALGORITHM.                                         *
C     *                                                              *
C     *  4. PARAMETERS:                                              *
C     *     INPUT:                                                   *
C     *        N        DIMENSION OF THE VECTOR                      *
C     *        A(I)     VECTOR TO SORT (INTEGER)  I=1,...,N          *
C     *        B(I)     = I    (INTEGER)  I=1,...,N                  *
C     *     OUTPUT:                                                  *
C     *        A(I)     THE SORTED VECTOR                            *
C     *        B(I)     PERMUTATION VECTOR OF THE SORTED VECTOR      *
C     *                                                              *
C     *  5. EXTERNAL SUBROUTINES:                                    *
C     *     NONE                                                     *
C     *                                                              *
C *** ****************************************************************
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IMPLICIT INTEGER (A-Z)
      DIMENSION A(1), B(1)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      F=1
      IF(L.LE.F) RETURN
      N2 =(L-F+1)/2
      S =1023
      DO 100 T=1,10
      IF(S.GT.N2) GOTO 90
      LS=L-S
      DO 20 I=F,LS
      IS =I+S
      AH =A(IS)
      BH = B(IS)
      J=I
      JS =IS
    5 IF(AH.GE.A(J)) GOTO 10
      A(JS) =A(J)
      B(JS) = B(J)
      JS =J
      J =J-S
      IF(J.GE.F) GOTO 5
   10 A(JS) =AH
      B(JS) = BH
   20 CONTINUE
   90 S=S/2
  100 CONTINUE
      RETURN
      END

        subroutine ilese(n, a, b, c, ndim)
        integer a( ndim, ndim), b(ndim, ndim), c(ndim, ndim)
        read (5,*) n
        if (n .gt. ndim) then
                print *,'n zu gross: n, nmax=', n, ndim
                stop
                        endif
        do 10 i=1, n
10              read (5,*) (a(i, j), j=1, n)
        do 20 i=1, n
20              read (5,*) (b(i, j), j=1, n)
        do 30 i=1, n
30              read (5,*, end = 40) (c(i, j), j=1, n)
        return
40      continue
        print *,'no linear term.'
        do 50 i=1,n
        do 50 j=1,n
50              c(i,j) = 0
        return
        end

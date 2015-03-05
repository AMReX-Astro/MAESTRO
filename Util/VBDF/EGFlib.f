C-----------------------------------------------------------------------
      SUBROUTINE EGFPAR ( NP, T, X, Y, CP, WEG, IWEG )
C-----------------------------------------------------------------------
C
C     This subroutine initializes the thermomolecular
C     parameters that are needed in order to evaluate
C     the transport linear systems.
C     The parameters that have to be evaluated depend on the transport
C     coefficients that will be subsequently computed, as indicated
C     by JFLAG. This flag has been initialized when calling EGINI.
C
C     Input
C     -----
C        NP           number of nodes
C        T(NP)        temperature
C        X(NP,NS)     species mole fractions
C        Y(NP,NS)     species mass fractions
C        CP(NP,NS)    species heat capacities at constant pressure
C                     per unit mass
C        WEG          double precision work array for EGLIB
C        IWEG         integer work array for EGLIB
C
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      DIMENSION WEG(*), IWEG(*)
C-----------------------------------------------------------------------
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LFPAR ( NP, JFLAG, WEG(IEGPA), 
     &       T, WEG(IDLT1), WEG(IDLT2), WEG(IDLT3), 
     &       WEG(IDLT4), WEG(IDLT5), WEG(IDLT6),
     &       WEG(IEGRU), NS, WEG(IAAA), WEG(IBBB),
     &       WEG(IEGWT), WEG(IBIN), WEG(IETA), 
     &       WEG(IETALG), WEG(ICTAIJ), WEG(IFITA), 
     &       WEG(ICINT), CP, WEG(ICXI), 
     &       WEG(IEPSIJ), WEG(IEGEPS), WEG(IEGCFD), WEG(IEGCFE), 
     &       WEG(IEGZRT), WEG(IEGDIP), IWEG(IEGLIN) )
C-----------------------------------------------------------------------
      CALL EGFXY ( NP, NS, WEG(IEGWT), WEG(IXTR), WEG(IYTR), WEG(IAUX), 
     &             X, Y, WEG(IWWTR), WEG(ISUMTR) )
C-----------------------------------------------------------------------
      RETURN
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE LFPAR ( NP, IFLAG, PATMOS, 
     &                   T, DLT1, DLT2, DLT3, DLT4, DLT5, DLT6,
     &                   RU, NS, AAA, CROT, WT, BIN, ETA,
     &                   ETALG, CTAIJ, FITA, 
     &                   CINT, CPMS, CXI, EPSIJ, 
     &                   EPS, COFD, COFE, ZROT, DIP, LIN )
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      DIMENSION BIN(NP,*), ETA(NP,NS), ETALG(NP,NS),
     &          CTAIJ(*), WT(*), 
     &          EPSIJ(*), LIN(NS), ZROT(NS), EPS(NS), DIP(NS)
      DIMENSION T(NP), DLT1(NP), DLT2(NP), DLT3(NP), DLT4(NP),
     &          DLT5(NP), DLT6(NP), AAA(NP), CROT(NP)
      DIMENSION CPMS(NP,NS), CINT(NP,NS), CXI(NP,NS)
      DIMENSION FITA(7,NS,NS)
C-----------------------------------------------------------------------
      DATA PI/3.1415926535D0/
      DATA PI32O2/2.7842D+00/, P2O4P2/4.4674D+00/, PI32/5.5683D+00/
C-----------------------------------------------------------------------
C     logarithm of the temperature and its powers
C-----------------------------------------------------------------------
      do nn = 1, np
         dlt1(nn) = dlog ( t(nn) )
         dlt2(nn) = dlt1(nn) * dlt1(nn)
         dlt3(nn) = dlt2(nn) * dlt1(nn)
         dlt4(nn) = dlt3(nn) * dlt1(nn)
         dlt5(nn) = dlt4(nn) * dlt1(nn)
         dlt6(nn) = dlt5(nn) * dlt1(nn)
      enddo
C-----------------------------------------------------------------------
      CALL EGFCOFE ( NP, COFE, NS, DLT1, DLT2, DLT3, ETA, ETALG )
      IF ( IFLAG .LE. 1 ) RETURN
C-----------------------------------------------------------------------
      CALL EGFCOFD ( NP, COFD, NS, DLT1, DLT2, DLT3, BIN )
      IF (IFLAG .LE. 2) RETURN
C-----------------------------------------------------------------------
C     Test for IFLAG = 3
C-----------------------------------------------------------------------
      IF (IFLAG .EQ. 3) RETURN
C-----------------------------------------------------------------------
C         COMPUTE PARKER CORRECTION FOR ZROT
C         AND ALSO THE ROTATIONAL AND INTERNAL PARTS OF SPECIFIC HEAT
C-----------------------------------------------------------------------
      DO 400 K = 1, NS
         IF (LIN(K) .EQ. 0) THEN
            do nn = 1, np
               CROT(nn)   = 0.0D0
               CINT(nn,K) = 0.0D0
            enddo
         ELSEIF (LIN(K) .EQ. 1) THEN
            WRU = WT(K) / RU
            do nn = 1, np
               CROT(nn)   = 1.0D0
               CINT(nn,K) = CPMS(nn,K) * WRU - 2.5D0
            enddo
         ELSEIF (LIN(K) .EQ. 2) THEN
            WRU = WT(K) / RU
            do nn = 1, np
               CROT(nn)   = 1.5D0
               CINT(nn,K) = CPMS(nn,K) * WRU - 2.5D0
            enddo
         ENDIF
C........
         DR   = EPS(K) / 298.0D0
         SQDR = DSQRT(DR)
         DR32 = SQDR*DR
         AAAA = (1.0D0 + PI32O2*SQDR + P2O4P2*DR + PI32*DR32) 
     &          * MAX(1.0D0, ZROT(K))
C........
         do nn = 1, np
            DD   = EPS(K) / T(nn)
            SQDD = DSQRT(DD)
            DD32 = SQDD*DD
            BBBB = (1.0D0 + PI32O2*SQDD + P2O4P2*DD + PI32*DD32) 
C........
            XI = ( AAAA / BBBB ) 
            CXI(nn,K) = CROT(nn) / ( XI * PI )
         enddo
  400 CONTINUE
C-----------------------------------------------------------------------
C     Test for IFLAG = 4 or 5
C-----------------------------------------------------------------------
      IF (IFLAG .LE. 5) RETURN
C-----------------------------------------------------------------------
C        EVALUATE THE BINARY INTERNAL DIFFUSION COEFFICIENTS
C                     D_int,ii
C-----------------------------------------------------------------------
      DO I = 1, NS
         II = (I-1) * NS + I
         IF ( IFLAG .EQ. 7 ) THEN
            FITA0 = FITA(1,I,I)
            FITA1 = FITA(2,I,I)
            FITA2 = FITA(3,I,I)
            FITA3 = FITA(4,I,I)
            FITA4 = FITA(5,I,I)
            FITA5 = FITA(6,I,I)
            FITA6 = FITA(7,I,I)
            do nn = 1, np
               AAA(nn) = FITA0          + FITA1*DLT1(nn) 
     &                 + FITA2*DLT2(nn) + FITA3*DLT3(nn)
     &                 + FITA4*DLT4(nn) + FITA5*DLT5(nn)
     &                 + FITA6*DLT6(nn) 
            enddo
         ELSE
            do nn = 1, np
               AAA(nn) = CTAIJ(II)
            enddo
         ENDIF
         IBIN = NS*(I-1) - (I*(I-1))/2 + I
         do nn = 1, np
            BIN(nn,IBIN) = 5.0D0 * PATMOS * WT(I) / 
     &                   ( 6.0D0 * RU * T(nn) * AAA(nn) * ETA(nn,I) )
         enddo
         IF ( DIP(I) .NE. 0.0D0 ) THEN
            do nn = 1, np
               BIN(nn,IBIN) = BIN(nn,IBIN) 
     &                      * ( 1.0D0 + 2.985D3 / T(nn) ** 1.5D0 )
            enddo
         ENDIF
      ENDDO
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGFCOFD ( NP, COF, NS, DLT, DLT2, DLT3, BIN )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION COF(4,NS,*), BIN(NP,*)
      DIMENSION DLT(NP), DLT2(NP), DLT3(NP)
C-----------------------------------------------------------------------
C     This subroutine returns the array BIN formed by the
C     reciprocals of the binary diffusion coefficients at
C     atmospheric pressure.
C-----------------------------------------------------------------------
      DO 120 K=1, NS
         KP1 = K + 1
         DO 110 L=KP1, NS
            COF1 = -COF(1,K,L)
            COF2 = -COF(2,K,L)
            COF3 = -COF(3,K,L)
            COF4 = -COF(4,K,L)
            KL   = NS*(K-1) - (K*(K-1))/2 + L
            do nn = 1, np
               DLK = DEXP( COF1 + COF2*DLT(nn) + COF3*DLT2(nn)
     1                          + COF4*DLT3(nn) )
               BIN(nn,KL) = DLK
            enddo
110      CONTINUE
120   CONTINUE
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGFCOFE ( NP, COF, NS, DLT, DLT2, DLT3, ETA, ETALG )
C-----------------------------------------------------------------------
C     This subroutine returns the arrays ETA and ETALG formed by the
C     pure species shear viscosities and their logs'.
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION COF(4,NS), ETA(NP,NS), ETALG(NP,NS)
      DIMENSION DLT(NP), DLT2(NP), DLT3(NP)
C-----------------------------------------------------------------------
      DO 120 K=1, NS
            COF1 = COF(1,K)
            COF2 = COF(2,K)
            COF3 = COF(3,K)
            COF4 = COF(4,K)
            do nn = 1, np
               ETALG(nn,K) = COF1 + COF2*DLT(nn) + COF3*DLT2(nn) 
     &                            + COF4*DLT3(nn) 
               ETA(nn,K)   = DEXP( ETALG(nn,K) )
            enddo
120   CONTINUE
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGFXY ( NP, NS, WT, XTR, YTR, AUX, X, Y, WWTR, SUMTR )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      DIMENSION WT(NS), XTR(NP,NS), YTR(NP,NS), AUX(NP,NS), X(NP,NS), 
     &          Y(NP,NS), WWTR(NP), SUMTR(NP)
C-----------------------------------------------------------------------
C     Add a small constant to the mole and mass fractions
C-----------------------------------------------------------------------
      CALL EGZERO ( np, SUMTR )
      DO I = 1, NS
         do nn = 1, np
            SUMTR(nn) = SUMTR(nn) + X(nn,I)
         enddo
      ENDDO
      SSS = 1.0D-16
      AAA = 1.0D0 / DFLOAT(NS)
      CALL EGZERO ( np, WWTR )
      DO I = 1, NS
         do nn = 1, np
            XTR(nn,I) = X(nn,I) + SSS * ( SUMTR(nn)*AAA - X(nn,I) )
            WWTR(nn)  = WWTR(nn) + XTR(nn,I) * WT(I)
         enddo
      ENDDO
      DO I = 1, NS
         do nn = 1, np
            YTR(nn,I) = XTR(nn,I) * WT(I) / WWTR(nn)
         enddo
      ENDDO
C-----------------------------------------------------------------------
C     AUX(nn,i) = \sum_{j .ne. i} YTR(nn,j)
C-----------------------------------------------------------------------
      CALL EGZERO ( np*NS, AUX )
      DO I = 1, NS
         DO J = I+1, NS
            do nn = 1, np
               AUX(nn,I) = AUX(nn,I) + YTR(nn,J)
               AUX(nn,J) = AUX(nn,J) + YTR(nn,I)
            enddo
         ENDDO
      ENDDO
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGFBIN (NP, T, WEG, BIN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 2
C        ITLS  = 0
C
C     Input
C     -----
C        NP             number of nodes
C        T(NP)          temperature
C        WEG            double precision work array for EGLIB
C
C     Output
C     ------
C        BIN(NP,NS,NS)  binary diffusion coefficients at atmospheric 
C                       pressure for species pairs (i,j) with i<>j.
C
C-----------------------------------------------------------------------
      INCLUDE 'eg.cmn'
      DIMENSION WEG(*), BIN(NP,NS,NS)
C-----------------------------------------------------------------------
      IND = 0
      DO I = 1, NS
         IND = IND + NP
         DO J = I+1, NS
            do nn = 1, np
               BIN(nn,I,J) = 1.0D0 / WEG ( IBIN + IND + nn - 1 )
               BIN(nn,J,I) = BIN(nn,I,J)
            enddo
            IND = IND + NP
         ENDDO
      ENDDO
C-----------------------------------------------------------------------      
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGFD1 (NP, PRES, IPRES, T, Y, WW, WEG, D)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 2
C        ITLS  = 1
C
C     Input
C     -----
C        NP           number of nodes
C        PRES(*)      pressure
C        IPRES        flag for pressure (0:scalar, .ne.0:node dependent)
C        T(NP)        temperature
C        Y(NP,NS)     species mass fractions
C        WW(NP)       mean molecular weight
C        WEG          double precision work array for EGLIB
C
C     Output
C     ------
C        D(NP,NS,NS)  flux diffusion matrix
C
C     
C     Two standard iterations are performed on the matrix L_[00]
C
C     Note
C     ----
C        D satisfies the mass conservation constraints
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGFD1 ( NP, NS, PRES, IPRES, WEG(IXTR), WEG(IYTR), 
     &              WEG(ISUMTR), WEG(IEGPA), D, WEG(IBIN), WEG(IG), 
     &              WEG(IDMI), WEG(IAAA), WEG(ITEMP), WEG(IAUX) )
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGFD1 ( NP, NS, PRES, IPRES, XTR, YTR, SUMTR, 
     &                    PATMOS, D, BIN, G, DMI, AAA, TEMP, AUX )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      CALL EGFEML00 ( NP, NS, XTR, BIN, G )
C-----------------------------------------------------------------------
C        Call the iterative procedures
C-----------------------------------------------------------------------
      ITERMX = 2
      CALL EGFSI1 ( NP, NS, G, DMI, TEMP, YTR, SUMTR, D, ITERMX, 
     &              AUX, AAA )
C-----------------------------------------------------------------------
C        Correct for the pressure dependence of D
C-----------------------------------------------------------------------
      IF ( IPRES .EQ. 0 ) THEN
         AA = PATMOS / PRES
         CALL DSCAL ( np*NS*NS, AA, D, 1 )
      ELSE
         CALL EGFSCAD ( NP, NS, PATMOS, PRES, D, AAA )
      ENDIF
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGFD2 (NP, PRES, IPRES, T, Y, WW, WEG, D)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 2
C        ITLS  = 1
C
C     Input
C     -----
C        NP           number of nodes
C        PRES(*)      pressure
C        IPRES        flag for pressure (0:scalar, .ne.0:node dependent)
C        T(NP)        temperature
C        Y(NP,NS)     species mass fractions
C        WW(NP)       mean molecular weight
C        WEG          double precision work array for EGLIB
C
C     Output
C     ------
C        D(NP,NS,NS)  flux diffusion matrix
C
C     
C     We form the Choleski decomposition of matrix L_[00] 
C
C     Note
C     ----
C        D satisfies the mass conservation constraints
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGFD2 ( NP, NS, PRES, IPRES, WEG(IXTR), WEG(IYTR), 
     &              WEG(ISUMTR), WEG(IEGPA), D, WEG(IBIN), 
     &              WEG(IG), WEG(IAN), WEG(IAAA), WEG(ITEMP) )
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGFD2 ( NP, NS, PRES, IPRES, XTR, YTR, SUMTR, 
     &                    PATMOS, D, BIN, G, AN, AAA, TEMP )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C     Form the linear system
C-----------------------------------------------------------------------
      CALL EGFDDEC ( NP, NS, XTR, YTR, G, BIN, TEMP )
C-----------------------------------------------------------------------
C     Form the matrix D
C-----------------------------------------------------------------------
      CALL EGFEVD ( NP, NS, NS, AN, YTR, G, SUMTR, D, TEMP, AAA )
C-----------------------------------------------------------------------
C     Project the matrix D
C-----------------------------------------------------------------------
      CALL EGFPDP ( NP, NS, YTR, TEMP, SUMTR, D, AAA )
C-----------------------------------------------------------------------
C        Correct for the pressure dependence of D
C-----------------------------------------------------------------------
      IF ( IPRES .EQ. 0 ) THEN
         AA = PATMOS / PRES
         CALL DSCAL ( np*NS*NS, AA, D, 1 )
      ELSE
         CALL EGFSCAD ( NP, NS, PATMOS, PRES, D, AAA )
      ENDIF
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGFDR1 (NP, T, Y, WW, WEG, D)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 2
C        ITLS  = 1
C
C     Input
C     -----
C        NP           number of nodes
C        T(NP)        temperature
C        Y(NP,NS)     species mass fractions
C        WW(NP)       mean molecular weight
C        WEG          double precision work array for EGLIB
C
C     Output
C     ------
C        D(NP,NS,NS)  rho * flux diffusion matrix
C
C        where rho is the density.
C
C     
C     Two standard iterations are performed on the matrix L_[00]
C
C     Note
C     ----
C        D satisfies the mass conservation constraints
C        D is also PRESSURE INDEPENDENT.
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGFDR1 ( NP, NS, T, WEG(IXTR), WEG(IYTR), WEG(IWWTR), 
     &               WEG(ISUMTR), WEG(IEGRU), WEG(IEGPA), D, WEG(IBIN), 
     &               WEG(IG), WEG(IDMI), WEG(ITEMP), WEG(IAUX),
     &               WEG(IAAA) )
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGFDR1 ( NP, NS, TEMPER, XTR, YTR, WWTR, SUMTR,  
     &                     RU, PATMOS, D, BIN, G, DMI, TEMP, AUX, AAA )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      CALL EGFEML00 ( NP, NS, XTR, BIN, G )
C-----------------------------------------------------------------------
C        Call the iterative procedures
C-----------------------------------------------------------------------
      ITERMX = 2
      CALL EGFSI1 ( NP, NS, G, DMI, TEMP, YTR, SUMTR, D, ITERMX, 
     &              AUX, AAA )
C-----------------------------------------------------------------------
      CALL EGFSCADR ( NP, NS, PATMOS, WWTR, RU, TEMPER, D, AAA )
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGFDR2 (NP, T, Y, WW, WEG, D)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 2
C        ITLS  = 1
C
C     Input
C     -----
C        NP           number of nodes
C        T(NP)        temperature
C        Y(NP,NS)     species mass fractions
C        WW(NP)       mean molecular weight
C        WEG          double precision work array for EGLIB
C
C     Output
C     ------
C        D(NP,NS,NS)  rho * flux diffusion matrix
C
C        where rho is the density.
C
C     
C     We form the Choleski decomposition of matrix L_[00] 
C
C     Note
C     ----
C        D satisfies the mass conservation constraints
C        D is also PRESSURE INDEPENDENT.
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGFDR2 ( NP, NS, T, WEG(IXTR), WEG(IYTR), WEG(ISUMTR), 
     &               WEG(IWWTR), WEG(IEGRU), WEG(IEGPA), D, WEG(IBIN), 
     &               WEG(IG), WEG(IAN), WEG(ITEMP), WEG(IAAA) )
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGFDR2 ( NP, NS, TEMPER, XTR, YTR, SUMTR, WWTR, 
     &                     RU, PATMOS, D, BIN, G, AN, TEMP, AAA )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C     Form the linear system
C-----------------------------------------------------------------------
      CALL EGFDDEC ( NP, NS, XTR, YTR, G, BIN, TEMP )
C-----------------------------------------------------------------------
C     Form the matrix D
C-----------------------------------------------------------------------
      CALL EGFEVD ( NP, NS, NS, AN, YTR, G, SUMTR, D, TEMP, AAA )
C-----------------------------------------------------------------------
C     Project the matrix D
C-----------------------------------------------------------------------
      CALL EGFPDP ( NP, NS, YTR, TEMP, SUMTR, D, AAA )
C-----------------------------------------------------------------------
      CALL EGFSCADR ( NP, NS, PATMOS, WWTR, RU, TEMPER, D, AAA )
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGFE1 ( NP, ALPHA, T, X, WEG, ETA )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 1
C        ITLS  = 0
C
C     Input
C     -----
C        NP           number of nodes
C        ALPHA        parameter for the averaging formula
C        T(NP)        temperature
C        X(NP,NS)     species mole fractions
C        WEG          double precision work array for EGLIB
C
C     Output
C     ------
C        ETA(NP)      shear viscosity 
C
C     
C     The empirical, mixture average formula of order ALPHA is used
C     to evaluate ETA.
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGFE1 ( NP, NS, X, ETA, WEG(IETALG), ALPHA, WEG(IAAA) ) 
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGFE1 ( NP, NS, X, ETAMA, ETALG, ALPHA, SUM )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      DIMENSION X(NP,NS), ETALG(NP,NS), ETAMA(NP), SUM(NP)
C-----------------------------------------------------------------------
      call egzero ( np, sum )
      IF ( ALPHA .EQ. 0.0D0 ) THEN
         DO I = 1, NS
            do nn = 1, np
               SUM(nn) = SUM(nn) + X(nn,I) * ETALG(nn,I)
            enddo
         ENDDO
         do nn = 1, np
            ETAMA(nn) = DEXP(SUM(nn))
         enddo
      ELSE
         DO I = 1, NS
            do nn = 1, np
               SUM(nn) = SUM(nn) + X(nn,I) * DEXP ( ALPHA*ETALG(nn,I) ) 
            enddo
         ENDDO
         RALPHA = (1.0D0/ALPHA)
         do nn = 1, np
            ETAMA(nn) = SUM(nn) ** RALPHA
         enddo
      ENDIF
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGFE2 ( NP, T, Y, WEG, ETA )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 2
C        ITLS  = 1
C
C     Input
C     -----
C        NP           number of nodes
C        T(NP)        temperature
C        Y(NP,NS)     species mass fractions
C        WEG          double precision work array for EGLIB
C
C     Output
C     ------
C        ETA(NP)      shear viscosity 
C
C     
C     One CG iteration is performed on the matrix H
C     The collision integral A_ij is assumed to be temperature
C     independent.
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGFE2 ( NP, NS, T, WEG(IXTR), WEG(IEGWT), ETA, 
     &            WEG(IEGRU), WEG(IEGPA),
     &            WEG(IBETA), WEG(IAN), WEG(IZN), WEG(IRN), 
     &            WEG(IG), WEG(IDMI), WEG(ITEMP), WEG(IAAA),
     &            WEG(IBBB), WEG(ICTAIJ), WEG(IETA), WEG(IBIN) )
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGFE2 ( NP, NS, T, X, WT, ETACG, RU, PATMOS, 
     &                  BETA, AN, ZN, RN, G, DMI, TEMP, AAA, BBB, 
     &                  AIJ, ETA, BIN )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      CALL EGFEMHA ( NP, NS, T, X, WT, RU, PATMOS, BETA, G, AIJ, 
     &               ETA, BIN )
C-----------------------------------------------------------------------
      ITERMX = 1
      CALL DCOPY  ( np*NS, BETA, 1, RN, 1 )
      CALL EGFCG1 ( NP, NS, G, DMI, AN, ZN, RN, TEMP, 
     &              AAA, BBB, ITERMX )
C-----------------------------------------------------------------------
      CALL EGFDOT ( NP, NS, ETACG, AN, BETA )
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGFE3 ( NP, T, Y, WEG, ETA )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 3
C        ITLS  = 1
C
C     Input
C     -----
C        NP           number of nodes
C        T(NP)        temperature
C        Y(NP,NS)     species mass fractions
C        WEG          double precision work array for EGLIB
C
C     Output
C     ------
C        ETA(NP)      shear viscosity 
C
C     
C     One CG iteration is performed on the matrix H
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGFE3 ( NP, NS, T, WEG(IDLT1), WEG(IDLT2), WEG(IDLT3),
     &           WEG(IDLT4), WEG(IDLT5), WEG(IDLT6),
     &           WEG(IXTR), WEG(IEGWT), ETA, 
     &           WEG(IEGRU), WEG(IEGPA),
     &           WEG(IBETA), WEG(IAN), WEG(IZN), WEG(IRN), 
     &           WEG(IG), WEG(IDMI), WEG(ITEMP), WEG(IAAA),
     &           WEG(IBBB), WEG(IAIJ), WEG(IFITA),
     &           WEG(IETA), WEG(IBIN) )
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGFE3 ( NP, NS, T, DLT1, DLT2, DLT3, DLT4, DLT5, DLT6,
     &                 X, WT, ETACG, RU, PATMOS, 
     &                 BETA, AN, ZN, RN, G, DMI, TEMP, AAA, BBB, 
     &                 AIJ, FITA, ETA, BIN )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      CALL EGFEMH ( NP, NS, T, DLT1, DLT2, DLT3, DLT4, DLT5, DLT6,
     &              X, WT, RU, PATMOS, BETA, G, AIJ, FITA, 
     &              ETA, BIN )
C-----------------------------------------------------------------------
      ITERMX = 1
      CALL DCOPY  ( np*NS, BETA, 1, RN, 1 )
      CALL EGFCG1 ( NP, NS, G, DMI, AN, ZN, RN, TEMP, 
     &              AAA, BBB, ITERMX )
C-----------------------------------------------------------------------
      CALL EGFDOT ( NP, NS, ETACG, AN, BETA )
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGFE4 ( NP, T, Y, WEG, ETA )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 3
C        ITLS  = 1
C
C     Input
C     -----
C        NP           number of nodes
C        T(NP)        temperature
C        Y(NP,NS)     species mass fractions
C        WEG          double precision work array for EGLIB
C
C     Output
C     ------
C        ETA(NP)      shear viscosity 
C
C     
C     We form the Choleski decomposition of matrix H
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGFE4 ( NP, NS, T, WEG(IDLT1), WEG(IDLT2), WEG(IDLT3),
     &             WEG(IDLT4), WEG(IDLT5), WEG(IDLT6),
     &             WEG(IXTR), WEG(IEGWT), ETA, 
     &             WEG(IEGRU), WEG(IEGPA),
     &             WEG(IBETA), WEG(IAN), WEG(IG), WEG(ITEMP),
     &             WEG(IAIJ), WEG(IFITA), WEG(IETA), WEG(IBIN) )
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGFE4 ( NP, NS, T, DLT1, DLT2, DLT3, DLT4, DLT5, DLT6,
     &                 X, WT, ETACG, RU, PATMOS, 
     &                 BETA, AN, G, TEMP, 
     &                 AIJ, FITA, ETA, BIN )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      CALL EGFEMH ( NP, NS, T, DLT1, DLT2, DLT3, DLT4, DLT5, DLT6,
     &              X, WT, RU, PATMOS, BETA, G, AIJ, FITA,
     &              ETA, BIN )
C-----------------------------------------------------------------------
      CALL EGFDEC ( NP, NS, G, TEMP, IER )
      CALL DCOPY  ( np*NS, BETA, 1, AN, 1 )
      CALL EGFSOL ( NP, NS, G, AN, TEMP )
C-----------------------------------------------------------------------
      CALL EGFDOT ( NP, NS, ETACG, AN, BETA )
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGFK1 ( NP, ALPHA, T, X, WEG, VOLVIS )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 4
C        ITLS  = 0
C
C     Input
C     -----
C        NP           number of nodes
C        ALPHA        parameter for the averaging formula
C        T(NP)        temperature
C        X(NP,NS)     species mole fractions
C        WEG          double precision work array for EGLIB
C
C     Output
C     ------
C        VOLVIS(NP)   volume viscosity 
C
C     
C     The empirical, mixture average formula of order ALPHA is used
C     to evaluate VOLVIS.
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGFK1 ( NP, NS, X, VOLVIS, WEG(IETALG), WEG(ICXI),
     &              WEG(ICINT), ALPHA, WEG(IAAA), WEG(IBBB) ) 
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGFK1 ( NP, NS, X, VOLVIS, ETALG, CXI, CINT, ALPHA,
     &                  SXP, SUM )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      DIMENSION X(NP,NS), ETALG(NP,NS), CXI(NP,NS), CINT(NP,NS)
      DIMENSION VOLVIS(NP), SXP(NP), SUM(NP)
C-----------------------------------------------------------------------
      call egzero ( np, sxp )
      DO I = 1, NS
         IF ( CXI(1,I) .NE. 0.0D0 ) THEN
            do nn = 1, np
               SXP(nn) = SXP(nn) + X(nn,I)
            enddo
         ENDIF
      ENDDO
      do nn = 1, np
         SXP(nn) = 1.0D0 / SXP(nn)
      enddo
      call egzero ( np, sum )
      IF ( ALPHA .EQ. 0.0D0 ) THEN
         DO I = 1, NS
            IF ( CXI(1,I) .NE. 0.0D0 ) THEN
               do nn = 1, np
                  CCC = CINT(nn,I) / ( 1.5D0 + CINT(nn,I) )
                  VVV = ETALG(nn,I) + DLOG ( 0.25D0*CCC*CCC/CXI(nn,I) )
                  SUM(nn) = SUM(nn) + SXP(nn) * X(nn,I) * VVV
               enddo
            ENDIF
         ENDDO
         do nn = 1, np
            VOLVIS(nn) = DEXP(SUM(nn))
         enddo
      ELSE
         DO I = 1, NS
            IF ( CXI(1,I) .NE. 0.0D0 ) THEN
               do nn = 1, np
                  CCC = CINT(nn,I) / ( 1.5D0 + CINT(nn,I) )
                  VVV = ETALG(nn,I) + DLOG ( 0.25D0*CCC*CCC/CXI(nn,I) )
                  SUM(nn) = SUM(nn) 
     &                    + SXP(nn) * X(nn,I) * DEXP ( ALPHA*VVV )
               enddo
            ENDIF
         ENDDO
         RALPHA = (1.0D0/ALPHA)
         do nn = 1, np
            VOLVIS(nn) = SUM(nn) ** RALPHA
         enddo
      ENDIF
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGFK2 (NP, T, Y, WEG, VV01)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 4
C        ITLS  = 1
C
C     Input
C     -----
C        NP           number of nodes
C        T(NP)        temperature
C        Y(NP,NS)     species mass fractions
C        WEG          double precision work array for EGLIB
C
C     Output
C     ------
C        VV01(NP)     volume viscosity 
C
C     
C     This subroutine evaluates the volume viscosity \kappa_[01]
C
C     The collision integral A_ij is assumed to be temperature
C     independent.
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGFK2 ( NP, NS, T, WEG(IXTR), WEG(IEGWT), WEG(IEGRU), 
     &              WEG(IEGPA), VV01, 
     &              WEG(IBIN), WEG(ICTAIJ), WEG(ICINT),
     &              WEG(IETA), WEG(ICXI), WEG(IG), WEG(IBETA),
     &              WEG(IAAA), WEG(IBBB) )
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGFK2 ( NP, NS, TEMPER, XTR, WT, RU, PATMOS, VV01, 
     &                    BIN, AIJ, CINT, ETA, CXI, G, BETA, AAA, BBB )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      DIMENSION G(NP,*), BETA(NP,*), VV01(NP)
C-----------------------------------------------------------------------
      CALL EGFEMKA01 ( NP, NS, TEMPER, XTR, WT, RU, PATMOS,
     &                 BIN, AIJ, CINT, ETA, CXI, G, BETA, AAA, BBB )
C-----------------------------------------------------------------------
C        Evaluate the transport coefficient
C-----------------------------------------------------------------------
      call egzero ( np, vv01 )
      DO I = 1, NS
         III  = NS*(I-1) - (I*(I-1))/2 + I
         do nn = 1, np
            VV01(nn) = VV01(nn) + BETA(nn,I) * BETA(nn,I) / G(nn,III)
         enddo
      ENDDO
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGFK3 (NP, T, Y, WEG, VV01)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 5
C        ITLS  = 1
C
C     Input
C     -----
C        NP           number of nodes
C        T(NP)        temperature
C        Y(NP,NS)     species mass fractions
C        WEG          double precision work array for EGLIB
C
C     Output
C     ------
C        VV01(NP)     volume viscosity 
C
C     
C     This subroutine evaluates the volume viscosity \kappa_[01]
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGFK3 ( NP, NS, T, WEG(IDLT1), WEG(IDLT2), WEG(IDLT3),
     &            WEG(IDLT4), WEG(IDLT5), WEG(IDLT6),
     &            WEG(IXTR), WEG(IEGWT), WEG(IEGRU), 
     &            WEG(IEGPA), VV01, 
     &            WEG(IBIN), WEG(IAIJ), WEG(IFITA), WEG(ICINT), 
     &            WEG(IETA), WEG(ICXI), WEG(IG), WEG(IBETA),
     &            WEG(IAAA), WEG(IBBB) )
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGFK3 ( NP, NS, TEMPER, DLT1, DLT2, DLT3, DLT4, 
     &                    DLT5, DLT6, XTR, WT, RU, PATMOS, VV01, 
     &                    BIN, AIJ, FITA, CINT, ETA, CXI, G, BETA, 
     &                    AAA, BBB )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      DIMENSION G(NP,*), BETA(NP,*), VV01(NP)
C-----------------------------------------------------------------------
      CALL EGFEMK01 ( NP, NS, TEMPER, DLT1, DLT2, DLT3, DLT4, 
     &                DLT5, DLT6, XTR, WT, RU, PATMOS,
     &                BIN, AIJ, FITA, CINT, ETA, CXI, G, BETA, 
     &                AAA, BBB )
C-----------------------------------------------------------------------
C        Evaluate the transport coefficient
C-----------------------------------------------------------------------
      call egzero ( np, vv01 )
      DO I = 1, NS
         III  = NS*(I-1) - (I*(I-1))/2 + I
         do nn = 1, np
            VV01(nn) = VV01(nn) + BETA(nn,I) * BETA(nn,I) / G(nn,III)
         enddo
      ENDDO
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGFK4 (NP, T, Y, WEG, VV)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 4
C        ITLS  = 2
C
C     Input
C     -----
C        NP           number of nodes
C        T(NP)        temperature
C        Y(NP,NS)     species mass fractions
C        WEG          double precision work array for EGLIB
C
C     Output
C     ------
C        VV(NP)       volume viscosity 
C
C     
C     This subroutine evaluates the volume viscosity 
C                  VV = \kappa_[s]^1 + \kappa_[01]
C     where \kappa_[s]^1 is obtained after one standard iteration
C     on the Schur complement K_[s].
C
C     The collision integral A_ij is assumed to be temperature
C     independent.
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGFK4 ( NP, NS, T, WEG(IXTR), WEG(IEGWT), WEG(IEGRU), 
     &              WEG(IEGPA), VV, 
     &              WEG(IBIN), WEG(ICTAIJ), WEG(ICINT),
     &              WEG(IETA),  WEG(ICXI), WEG(IG), WEG(IBETA),
     &              WEG(IAN), WEG(ITEMP), WEG(IAAA), WEG(IBBB) )
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGFK4 ( NP, NS, TEMPER, XTR, WT, RU, PATMOS, VV, 
     &           BIN, AIJ, CINT, ETA, CXI, G, BETA, AN, TEMP, AAA, BBB )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      CALL EGFEMKA ( NP, NS, TEMPER, XTR, WT, RU, PATMOS,
     &               BIN, AIJ, CINT, ETA, CXI, G, BETA, AAA, BBB )
C-----------------------------------------------------------------------
C        Evaluate the transport coefficient
C-----------------------------------------------------------------------
      CALL EGFEVK ( NP, NS, VV, G, AN, BETA, TEMP, AAA, BBB )
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGFK5 (NP, T, Y, WEG, VV)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 5
C        ITLS  = 2
C
C     Input
C     -----
C        NP           number of nodes
C        T(NP)        temperature
C        Y(NP,NS)     species mass fractions
C        WEG          double precision work array for EGLIB
C
C     Output
C     ------
C        VV(NP)       volume viscosity 
C
C     
C     This subroutine evaluates the volume viscosity 
C                  VV = \kappa_[s]^1 + \kappa_[01]
C     where \kappa_[s]^1 is obtained after one standard iteration
C     on the Schur complement K_[s].
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGFK5 ( NP, NS, T, WEG(IDLT1), WEG(IDLT2), WEG(IDLT3),
     &           WEG(IDLT4), WEG(IDLT5), WEG(IDLT6),
     &           WEG(IXTR), WEG(IEGWT), WEG(IEGRU), 
     &           WEG(IEGPA), VV, 
     &           WEG(IBIN), WEG(IAIJ), WEG(IFITA), WEG(ICINT), 
     &           WEG(IETA), WEG(ICXI), WEG(IG), WEG(IBETA),
     &           WEG(IAN), WEG(ITEMP), WEG(IAAA), WEG(IBBB) )
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGFK5 ( NP, NS, TEMPER, DLT1, DLT2, DLT3, DLT4, 
     &           DLT5, DLT6, XTR, WT, RU, PATMOS, VV, 
     &           BIN, AIJ, FITA, CINT, ETA, CXI, G, BETA, AN, TEMP, 
     &           AAA, BBB )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      CALL EGFEMK ( NP, NS, TEMPER, DLT1, DLT2, DLT3, DLT4, 
     &           DLT5, DLT6, XTR, WT, RU, PATMOS,
     &           BIN, AIJ, FITA, CINT, ETA, CXI, G, BETA, AAA, BBB )
C-----------------------------------------------------------------------
C        Evaluate the transport coefficient
C-----------------------------------------------------------------------
      CALL EGFEVK ( NP, NS, VV, G, AN, BETA, TEMP, AAA, BBB )
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGFK6 ( NP, T, Y, WEG, VV )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 5
C        ITLS  = 2
C
C     Input
C     -----
C        NP           number of nodes
C        T(NP)        temperature
C        Y(NP,NS)     species mass fractions
C        WEG          double precision work array for EGLIB
C
C     Output
C     ------
C        VV(NP)       volume viscosity 
C
C     
C     We form the Choleski decomposition of matrix K
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGFK6 ( NP, NS, T, WEG(IDLT1), WEG(IDLT2), WEG(IDLT3),
     &           WEG(IDLT4), WEG(IDLT5), WEG(IDLT6),
     &           WEG(IXTR), WEG(IEGWT), WEG(IEGRU), 
     &           WEG(IEGPA), VV, 
     &           WEG(IBIN), WEG(IAIJ), WEG(IFITA), WEG(ICINT), 
     &           WEG(IETA), WEG(ICXI), WEG(IG), WEG(IBETA),
     &           WEG(IAN), WEG(ITEMP), WEG(IAAA), WEG(IBBB) )
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGFK6 ( NP, NS, TEMPER, DLT1, DLT2, DLT3, DLT4, 
     &           DLT5, DLT6, XTR, WT, RU, PATMOS, VV, 
     &           BIN, AIJ, FITA, CINT, ETA, CXI, G, BETA, AN, TEMP, 
     &           AAA, BBB )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      DIMENSION G(NP,*)
      DIMENSION CINT(NP,NS), XTR(NP,NS) 
C-----------------------------------------------------------------------
      CALL EGFEMK ( NP, NS, TEMPER, DLT1, DLT2, DLT3, DLT4, 
     &           DLT5, DLT6, XTR, WT, RU, PATMOS,
     &           BIN, AIJ, FITA, CINT, ETA, CXI, G, BETA, AAA, BBB )
C-----------------------------------------------------------------------
C     Make the matrix G positive definite
C-----------------------------------------------------------------------
      NG = 2 * NS
      DO I = 1, NS
         II1 = NG*(I-1) - (I*(I-1))/2 + I
         II2 = NG*(I-1) - (I*(I-1))/2 + I + NS
         II3 = NG*(I+NS-1) - ((I+NS)*(I+NS-1))/2 + I + NS
C........
         do nn = 1, np
            G(nn,II1) = G(nn,II1) + 1.5D0      * 1.5D0      
     &                                         * XTR(nn,I) * XTR(nn,I)
            G(nn,II2) = G(nn,II2) + 1.5D0      * CINT(nn,I) 
     &                                         * XTR(nn,I) * XTR(nn,I)
            G(nn,II3) = G(nn,II3) + CINT(nn,I) * CINT(nn,I) 
     &                                         * XTR(nn,I) * XTR(nn,I)
         enddo
         DO J = I+1, NS
            IJ1 = NG*(I-1) - (I*(I-1))/2 + J
            IJ2 = NG*(I-1) - (I*(I-1))/2 + J + NS
            JI2 = NG*(J-1) - (J*(J-1))/2 + I + NS
            IJ3 = NG*(I+NS-1) - ((I+NS)*(I+NS-1))/2 + J + NS
C...........
            do nn = 1, np
               FAC = XTR(nn,I) * XTR(nn,J)
               G(nn,IJ1) = G(nn,IJ1) + 1.5D0 * 1.5D0 * FAC      
               G(nn,IJ2) = G(nn,IJ2) + 1.5D0 * CINT(nn,J) * FAC 
               G(nn,JI2) = G(nn,JI2) + 1.5D0 * CINT(nn,I) * FAC 
               G(nn,IJ3) = G(nn,IJ3) + CINT(nn,I) * CINT(nn,J) * FAC 
            enddo
         ENDDO
      ENDDO
C-----------------------------------------------------------------------
C        Evaluate the transport coefficient
C-----------------------------------------------------------------------
      CALL EGFDEC ( NP, NG, G, TEMP, IER )
      CALL DCOPY  ( np*NG, BETA, 1, AN, 1 )
      CALL EGFSOL ( NP, NG, G, AN, TEMP )
C-----------------------------------------------------------------------
      CALL EGFDOT ( NP, NG, VV, AN, BETA )
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGFL1 ( NP, ALPHA, T, X, WEG, CON )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 1
C        ITLS  = 0
C
C     Input
C     -----
C        NP           number of nodes
C        ALPHA        parameter for the averaging formula
C        T(NP)        temperature
C        X(NP,NS)     species mole fractions
C        WEG          double precision work array for EGLIB
C
C     Output
C     ------
C        CON(NP)      thermal conductivity 
C
C     
C     The empirical, mixture average formula of order ALPHA is used
C     to evaluate CON.
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGFL1 ( NP, NS, 
     &              X, T, WEG(IDLT1), WEG(IDLT2), WEG(IDLT3),
     &              CON, WEG(IEGCFL), ALPHA, WEG(IAAA) ) 
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGFL1 ( NP, NS, X, T, DLT, DLT2, DLT3, 
     &                  CONMA, COF, ALPHA, SUM )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      DIMENSION X(NP,NS), COF(4,NS), CONMA(NP)
      DIMENSION T(NP), DLT(NP), DLT2(NP), DLT3(NP), SUM(NP)
C-----------------------------------------------------------------------
      call egzero ( np, sum )
      IF ( ALPHA .EQ. 0.0D0 ) THEN
         DO I = 1, NS
            COF1 = COF(1,I)
            COF2 = COF(2,I)
            COF3 = COF(3,I)
            COF4 = COF(4,I)
            do nn = 1, np
               CONLG = COF1 + COF2*DLT(nn) + COF3*DLT2(nn) 
     &                      + COF4*DLT3(nn) 
               SUM(nn) = SUM(nn) + X(nn,I) * CONLG
            enddo
         ENDDO
         do nn = 1, np
            CONMA(nn) = DEXP(SUM(nn))
         enddo
      ELSE
         DO I = 1, NS
            COF1 = COF(1,I)
            COF2 = COF(2,I)
            COF3 = COF(3,I)
            COF4 = COF(4,I)
            do nn = 1, np
               CONLG = COF1 + COF2*DLT(nn) + COF3*DLT2(nn) 
     &                      + COF4*DLT3(nn) 
               SUM(nn) = SUM(nn) + X(nn,I) * DEXP ( ALPHA*CONLG )
            enddo
         ENDDO
         RALPHA = (1.0D0/ALPHA)
         do nn = 1, np
            CONMA(nn) = SUM(nn) ** RALPHA
         enddo
      ENDIF
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGFL2 (NP, T, Y, WEG, IWEG, TC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 6
C        ITLS  = 1
C
C     Input
C     -----
C        NP           number of nodes
C        T(NP)        temperature
C        Y(NP,NS)     species mass fractions
C        WEG          double precision work array for EGLIB
C        IWEG         integer work array for EGLIB
C
C     Output
C     ------
C        TC(NP)       thermal conductivity 
C
C     
C     One CG iteration is performed on the matrix \Lambda_[e].
C
C     The collision integrals A_ij, B_ij and C_ij are assumed to be
C     temperature independent.
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*), IWEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGFL2 ( NP, NS, T, WEG(IXTR), WEG(IEGWT), WEG(IEGRU), 
     &              WEG(IEGPA), TC, WEG(IBIN), 
     &              WEG(ICTAIJ), WEG(ICTBIJ), 
     &              WEG(ICINT), WEG(IETA), WEG(ICXI), 
     &              IWEG(IEGLIN),
     &              WEG(IG), WEG(IDMI), WEG(IAN), WEG(IZN), 
     &              WEG(IRN), WEG(ITEMP), WEG(IBETA), 
     &              WEG(IAAA), WEG(IBBB) )
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGFL2 ( NP, NS, TEMPER, XTR, WT, RU, PATMOS,
     &           TC, 
     &           BIN, AIJ, BIJ, CINT, ETA, CXI, LIN,
     &           G, DMI, AN, ZN, RN, TEMP, BETA, AAA, BBB )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      CALL EGFEMAAE ( NP, NS, TEMPER, XTR, WT, RU, PATMOS,
     &                BIN, AIJ, BIJ, CINT, ETA, CXI,
     &                LIN, G, BETA )
C-----------------------------------------------------------------------
C        Call the iterative procedures
C-----------------------------------------------------------------------
      NG = NS
C.....
      ITERMX = 1
      CALL DCOPY  ( np*NG, BETA, 1, RN, 1 )
      CALL EGFCG1 ( NP, NG, G, DMI, AN, ZN, RN, TEMP, 
     &              AAA, BBB, ITERMX )
C-----------------------------------------------------------------------
C        Evaluate the transport coefficients
C-----------------------------------------------------------------------
      CALL EGFEVL ( NP, NG, TC, AN, BETA, PATMOS, TEMPER )
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGFL3 (NP, T, Y, WEG, IWEG, TC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 7
C        ITLS  = 1
C
C     Input
C     -----
C        NP           number of nodes
C        T(NP)        temperature
C        Y(NP,NS)     species mass fractions
C        WEG          double precision work array for EGLIB
C        IWEG         integer work array for EGLIB
C
C     Output
C     ------
C        TC(NP)       thermal conductivity 
C
C     
C     One CG iteration is performed on the matrix \Lambda_[e].
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*), IWEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGFL3 ( NP, NS, T, WEG(IDLT1), WEG(IDLT2), WEG(IDLT3),
     &             WEG(IDLT4), WEG(IDLT5), WEG(IDLT6),
     &             WEG(IXTR), WEG(IEGWT), WEG(IEGRU), 
     &             WEG(IEGPA), TC, 
     &             WEG(IBIN), WEG(IAIJ), WEG(IBIJ), 
     &             WEG(IFITA), WEG(IFITB), 
     &             WEG(ICINT), WEG(IETA), WEG(ICXI), 
     &             IWEG(IEGLIN),
     &             WEG(IG), WEG(IDMI), WEG(IAN), WEG(IZN), 
     &             WEG(IRN), WEG(ITEMP), WEG(IBETA), 
     &             WEG(IAAA), WEG(IBBB) )
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGFL3 ( NP, NS, TEMPER, DLT1, DLT2, DLT3, DLT4, 
     &           DLT5, DLT6, XTR, WT, RU, PATMOS, TC,
     &           BIN, AIJ, BIJ, FITA, FITB, 
     &           CINT, ETA, CXI, LIN,
     &           G, DMI, AN, ZN, RN, TEMP, BETA, AAA, BBB )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      CALL EGFEMAE ( NP, NS, TEMPER, DLT1, DLT2, DLT3, DLT4, 
     &           DLT5, DLT6, XTR, WT, RU, PATMOS,
     &           BIN, AIJ, BIJ, FITA, FITB, CINT, ETA, CXI,
     &           LIN, G, BETA )
C-----------------------------------------------------------------------
C        Call the iterative procedures
C-----------------------------------------------------------------------
      NG = NS
C.....
      ITERMX = 1
      CALL DCOPY  ( np*NG, BETA, 1, RN, 1 )
      CALL EGFCG1 ( NP, NG, G, DMI, AN, ZN, RN, TEMP, 
     &              AAA, BBB, ITERMX )
C-----------------------------------------------------------------------
C        Evaluate the transport coefficients
C-----------------------------------------------------------------------
      CALL EGFEVL ( NP, NG, TC, AN, BETA, PATMOS, TEMPER )
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGFL4 (NP, T, Y, WEG, IWEG, TC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 6
C        ITLS  = 2
C
C     Input
C     -----
C        NP           number of nodes
C        T(NP)        temperature
C        Y(NP,NS)     species mass fractions
C        WEG          double precision work array for EGLIB
C        IWEG         integer work array for EGLIB
C
C     Output
C     ------
C        TC(NP)       thermal conductivity 
C
C     
C     One CG iteration is performed on the matrix \Lambda.
C
C     The collision integrals A_ij, B_ij and C_ij are assumed to be
C     temperature independent.
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*), IWEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGFL4 ( NP, NS, T, WEG(IXTR), WEG(IEGWT), WEG(IEGRU), 
     &             WEG(IEGPA), TC, WEG(IBIN),
     &             WEG(ICTAIJ), WEG(ICTBIJ), 
     &             WEG(ICINT), WEG(IETA), WEG(ICXI),
     &             IWEG(IEGLIN), WEG(IG), WEG(IDMI), WEG(IAN),
     &             WEG(IZN),  WEG(IRN), WEG(ITEMP), WEG(IBETA),
     &             WEG(IAAA), WEG(IBBB) )
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGFL4 ( NP, NS, TEMPER, XTR, WT, RU, PATMOS, TC, 
     &           BIN, AIJ, BIJ, CINT, ETA, CXI, LIN,
     &           G, DMI, AN, ZN, RN, TEMP, BETA, AAA, BBB )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      CALL EGFEMAA ( NP, NS, TEMPER, XTR, WT, RU, PATMOS,
     &               BIN, AIJ, BIJ, CINT, ETA, CXI,
     &               LIN, G, BETA )
C-----------------------------------------------------------------------
C        Call the iterative procedures
C-----------------------------------------------------------------------
      NG = 2 * NS
C.....
      ITERMX = 1
      CALL DCOPY  ( np*NG, BETA, 1, RN, 1 )
      CALL EGFCG2 ( NP, NS, NG, G, DMI, AN, ZN, RN, TEMP, 
     &              AAA, BBB, ITERMX )
C-----------------------------------------------------------------------
C        Evaluate the transport coefficients
C-----------------------------------------------------------------------
      CALL EGFEVL ( NP, NG, TC, AN, BETA, PATMOS, TEMPER )
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGFL5 (NP, T, Y, WEG, IWEG, TC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 7
C        ITLS  = 2
C
C     Input
C     -----
C        NP           number of nodes
C        T(NP)        temperature
C        Y(NP,NS)     species mass fractions
C        WEG          double precision work array for EGLIB
C        IWEG         integer work array for EGLIB
C
C     Output
C     ------
C        TC(NP)       thermal conductivity 
C
C     
C     One CG iteration is performed on the matrix \Lambda.
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*), IWEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGFL5 ( NP, NS, T, WEG(IDLT1), WEG(IDLT2), WEG(IDLT3),
     &            WEG(IDLT4), WEG(IDLT5), WEG(IDLT6),
     &            WEG(IXTR), WEG(IEGWT), WEG(IEGRU), 
     &            WEG(IEGPA), TC, 
     &            WEG(IBIN), WEG(IAIJ), WEG(IBIJ), 
     &            WEG(IFITA), WEG(IFITB), 
     &            WEG(ICINT), WEG(IETA), WEG(ICXI), IWEG(IEGLIN),
     &            WEG(IG), WEG(IDMI), WEG(IAN), WEG(IZN), 
     &            WEG(IRN), WEG(ITEMP), WEG(IBETA),
     &            WEG(IAAA), WEG(IBBB) )
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGFL5 ( NP, NS, TEMPER, DLT1, DLT2, DLT3, DLT4, 
     &           DLT5, DLT6, XTR, WT, RU, PATMOS, TC, 
     &           BIN, AIJ, BIJ, FITA, FITB, 
     &           CINT, ETA, CXI, LIN,
     &           G, DMI, AN, ZN, RN, TEMP, BETA, AAA, BBB )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      dimension an(np,*), beta(np,*), temper(np), tc(np)
C-----------------------------------------------------------------------
      CALL EGFEMA ( NP, NS, TEMPER, DLT1, DLT2, DLT3, DLT4, 
     &           DLT5, DLT6, XTR, WT, RU, PATMOS,
     &           BIN, AIJ, BIJ, FITA, FITB, CINT, ETA, CXI,
     &           LIN, G, BETA )
C-----------------------------------------------------------------------
C        Call the iterative procedures
C-----------------------------------------------------------------------
      NG = 2 * NS
C.....
      ITERMX = 1
      CALL DCOPY  ( np*NG, BETA, 1, RN, 1 )
      CALL EGFCG2 ( NP, NS, NG, G, DMI, AN, ZN, RN, TEMP, 
     &              AAA, BBB, ITERMX )
C-----------------------------------------------------------------------
C        Evaluate the transport coefficients
C-----------------------------------------------------------------------
      CALL EGFEVL ( NP, NG, TC, AN, BETA, PATMOS, TEMPER )
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGFLCT1 (NP, T, Y, WEG, IWEG, TC, CHIT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 7
C        ITLS  = 1
C
C     Input
C     -----
C        NP           number of nodes
C        T(NP)        temperature
C        Y(NP,NS)     species mass fractions
C        WEG          double precision work array for EGLIB
C        IWEG         integer work array for EGLIB
C
C     Output
C     ------
C        TC(NP)        thermal conductivity 
C        CHIT(NP,NS)   rescaled thermal diffusion ratios
C
C     
C     Three CG iterations are performed on the matrix \Lambda_[e].
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*), IWEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGFLCT1 ( NP, NS, T, WEG(IDLT1), WEG(IDLT2), 
     &         WEG(IDLT3), WEG(IDLT4), WEG(IDLT5), WEG(IDLT6),
     &         WEG(IXTR), WEG(IEGWT), WEG(IEGRU), 
     &         WEG(IEGPA), TC, CHIT,
     &         WEG(IBIN), WEG(IAIJ), WEG(IBIJ), WEG(ICIJ), 
     &         WEG(IFITA), WEG(IFITB), WEG(IFITC),
     &         WEG(ICINT), WEG(IETA), WEG(ICXI), 
     &         IWEG(IEGLIN),
     &         WEG(IG), WEG(IDMI), WEG(IAN), WEG(IZN), 
     &         WEG(IRN), WEG(ITEMP), WEG(IBETA),
     &         WEG(IAAA), WEG(IBBB) )
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGFLCT1 ( NP, NS, TEMPER, DLT1, DLT2, DLT3, DLT4, 
     &           DLT5, DLT6, XTR, WT, RU, PATMOS,
     &           TC, CHIT,
     &           BIN, AIJ, BIJ, CIJ, FITA, FITB, FITC,
     &           CINT, ETA, CXI, LIN,
     &           G, DMI, AN, ZN, RN, TEMP, BETA, AAA, BBB )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      CALL EGFEMAE ( NP, NS, TEMPER, DLT1, DLT2, DLT3, DLT4, 
     &           DLT5, DLT6, XTR, WT, RU, PATMOS,
     &           BIN, AIJ, BIJ, FITA, FITB, CINT, ETA, CXI,
     &           LIN, G, BETA )
C-----------------------------------------------------------------------
C        Call the iterative procedures
C-----------------------------------------------------------------------
      NG = NS
C.....
      ITERMX = 3
      CALL DCOPY  ( np*NG, BETA, 1, RN, 1 )
      CALL EGFCG1 ( NP, NG, G, DMI, AN, ZN, RN, TEMP, 
     &              AAA, BBB, ITERMX )
C-----------------------------------------------------------------------
C        Evaluate the transport coefficients
C-----------------------------------------------------------------------
      CALL EGFEVL  ( NP, NG, TC, AN, BETA, PATMOS, TEMPER )
      CALL EGFEVCT ( NP, NS, CHIT, AN, XTR, BIN, CIJ, FITC,
     &               DLT1, DLT2, DLT3, DLT4, DLT5, DLT6, WT )
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGFLCT2 ( NP, T, Y, WEG, IWEG, TC, CHIT )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 7
C        ITLS  = 1
C
C     Input
C     -----
C        NP           number of nodes
C        T(NP)        temperature
C        Y(NP,NS)     species mass fractions
C        WEG          double precision work array for EGLIB
C        IWEG         integer work array for EGLIB
C
C     Output
C     ------
C        TC(NP)        thermal conductivity 
C        CHIT(NP,NS)   rescaled thermal diffusion ratios
C
C     
C     We form the Choleski decomposition of matrix \Lambda_[e] 
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*), IWEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGFLCT2 ( NP, NS, T, WEG(IDLT1), WEG(IDLT2),
     &         WEG(IDLT3), WEG(IDLT4), WEG(IDLT5), WEG(IDLT6),
     &         WEG(IXTR), WEG(IEGWT), WEG(IEGRU), 
     &         WEG(IEGPA), TC, CHIT,
     &         WEG(IBIN), WEG(IAIJ), WEG(IBIJ), WEG(ICIJ), 
     &         WEG(IFITA), WEG(IFITB), WEG(IFITC),
     &         WEG(ICINT), WEG(IETA), WEG(ICXI), 
     &         IWEG(IEGLIN),
     &         WEG(IG), WEG(IAN), WEG(ITEMP), WEG(IBETA) )
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGFLCT2 ( NP, NS, TEMPER, DLT1, DLT2, DLT3, DLT4, 
     &           DLT5, DLT6, XTR, WT, RU, PATMOS,
     &           TC, CHIT,
     &           BIN, AIJ, BIJ, CIJ, FITA, FITB, FITC,
     &           CINT, ETA, CXI, LIN,
     &           G, AN, TEMP, BETA )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      CALL EGFEMAE ( NP, NS, TEMPER, DLT1, DLT2, DLT3, DLT4, 
     &           DLT5, DLT6, XTR, WT, RU, PATMOS,
     &           BIN, AIJ, BIJ, FITA, FITB, CINT, ETA, CXI,
     &           LIN, G, BETA )
C-----------------------------------------------------------------------
C        Call the iterative procedures
C-----------------------------------------------------------------------
      NG = NS
      CALL EGFDEC ( NP, NG, G, TEMP, IER )
      CALL DCOPY  ( np*NG, BETA, 1, AN, 1 )
      CALL EGFSOL ( NP, NG, G, AN, TEMP )
C-----------------------------------------------------------------------
C        Evaluate the transport coefficients
C-----------------------------------------------------------------------
      CALL EGFEVL  ( NP, NG, TC, AN, BETA, PATMOS, TEMPER )
      CALL EGFEVCT ( NP, NS, CHIT, AN, XTR, BIN, CIJ, FITC,
     &               DLT1, DLT2, DLT3, DLT4, DLT5, DLT6, WT )
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGFLCT3 (NP, T, Y, WEG, IWEG, TC, CHIT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 7
C        ITLS  = 2
C
C     Input
C     -----
C        NP           number of nodes
C        T(NP)        temperature
C        Y(NP,NS)     species mass fractions
C        WEG          double precision work array for EGLIB
C        IWEG         integer work array for EGLIB
C
C     Output
C     ------
C        TC(NP)        thermal conductivity 
C        CHIT(NP,NS)   rescaled thermal diffusion ratios
C
C     
C     Three CG iterations are performed on the matrix \Lambda.
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*), IWEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGFLCT3 ( NP, NS, T, WEG(IDLT1), WEG(IDLT2),
     &         WEG(IDLT3), WEG(IDLT4), WEG(IDLT5), WEG(IDLT6),
     &         WEG(IXTR), WEG(IEGWT), WEG(IEGRU), 
     &         WEG(IEGPA), TC, CHIT,
     &         WEG(IBIN), WEG(IAIJ), WEG(IBIJ), WEG(ICIJ), 
     &         WEG(IFITA), WEG(IFITB), WEG(IFITC),
     &         WEG(ICINT), WEG(IETA), WEG(ICXI), 
     &         IWEG(IEGLIN),WEG(IG), WEG(IDMI), WEG(IAN), 
     &         WEG(IZN), WEG(IRN), WEG(ITEMP), WEG(IBETA),
     &         WEG(IAAA), WEG(IBBB) )
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGFLCT3 ( NP, NS, TEMPER, DLT1, DLT2, DLT3, DLT4, 
     &           DLT5, DLT6, XTR, WT, RU, PATMOS,
     &           TC, CHIT,
     &           BIN, AIJ, BIJ, CIJ, FITA, FITB, FITC, 
     &           CINT, ETA, CXI, LIN,
     &           G, DMI, AN, ZN, RN, TEMP, BETA, AAA, BBB )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      CALL EGFEMA ( NP, NS, TEMPER, DLT1, DLT2, DLT3, DLT4, 
     &           DLT5, DLT6, XTR, WT, RU, PATMOS,
     &           BIN, AIJ, BIJ, FITA, FITB, CINT, ETA, CXI,
     &           LIN, G, BETA )
C-----------------------------------------------------------------------
C        Call the iterative procedures
C-----------------------------------------------------------------------
      NG = 2 * NS
C.....
      ITERMX = 3
      CALL DCOPY  ( np*NG, BETA, 1, RN, 1 )
      CALL EGFCG2 ( NP, NS, NG, G, DMI, AN, ZN, RN, TEMP, 
     &              AAA, BBB, ITERMX )
C-----------------------------------------------------------------------
C        Evaluate the transport coefficients
C-----------------------------------------------------------------------
      CALL EGFEVL  ( NP, NG, TC, AN, BETA, PATMOS, TEMPER )
      CALL EGFEVCT ( NP, NS, CHIT, AN, XTR, BIN, CIJ, FITC,
     &               DLT1, DLT2, DLT3, DLT4, DLT5, DLT6, WT )
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGFLCT4 (NP, T, Y, WEG, IWEG, TC, CHIT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 7
C        ITLS  = 2
C
C     Input
C     -----
C        NP           number of nodes
C        T(NP)        temperature
C        Y(NP,NS)     species mass fractions
C        WEG          double precision work array for EGLIB
C        IWEG         integer work array for EGLIB
C
C     Output
C     ------
C        TC(NP)        thermal conductivity 
C        CHIT(NP,NS)   rescaled thermal diffusion ratios
C
C     
C     We form the Choleski decomposition of matrix \Lambda 
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*), IWEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGFLCT4 ( NP, NS, T, WEG(IDLT1), WEG(IDLT2),
     &         WEG(IDLT3), WEG(IDLT4), WEG(IDLT5), WEG(IDLT6),
     &         WEG(IXTR), WEG(IEGWT), WEG(IEGRU), 
     &         WEG(IEGPA), TC, CHIT,
     &         WEG(IBIN), WEG(IAIJ), WEG(IBIJ), WEG(ICIJ), 
     &         WEG(IFITA), WEG(IFITB), WEG(IFITC),
     &         WEG(ICINT), WEG(IETA), WEG(ICXI), IWEG(IEGLIN),
     &         WEG(IG), WEG(IAN), WEG(ITEMP), WEG(IBETA) )
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGFLCT4 ( NP, NS, TEMPER, DLT1, DLT2, DLT3, DLT4, 
     &           DLT5, DLT6, XTR, WT, RU, PATMOS,
     &           TC, CHIT,
     &           BIN, AIJ, BIJ, CIJ, FITA, FITB, FITC,
     &           CINT, ETA, CXI, LIN,
     &           G, AN, TEMP, BETA )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      CALL EGFEMA ( NP, NS, TEMPER, DLT1, DLT2, DLT3, DLT4, 
     &           DLT5, DLT6, XTR, WT, RU, PATMOS,
     &           BIN, AIJ, BIJ, FITA, FITB, CINT, ETA, CXI,
     &           LIN, G, BETA )
C-----------------------------------------------------------------------
      NG = 2 * NS
      CALL EGFDEC ( NP, NG, G, TEMP, IER )
      CALL DCOPY  ( np*NG, BETA, 1, AN, 1 )
      CALL EGFSOL ( NP, NG, G, AN, TEMP )
C-----------------------------------------------------------------------
C        Evaluate the transport coefficients
C-----------------------------------------------------------------------
      CALL EGFEVL  ( NP, NG, TC, AN, BETA, PATMOS, TEMPER )
      CALL EGFEVCT ( NP, NS, CHIT, AN, XTR, BIN, CIJ, FITC,
     &               DLT1, DLT2, DLT3, DLT4, DLT5, DLT6, WT )
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGFLTD1 (NP, PRES, IPRES, T, Y, WW, WEG, IWEG, 
     &                    PTC, THETA, D)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 7
C        ITLS  = 2
C
C     Input
C     -----
C        NP           number of nodes
C        PRES(*)      pressure
C        IPRES        flag for pressure (0:scalar, .ne.0:node dependent)
C        T(NP)        temperature
C        Y(NP,NS)     species mass fractions
C        WW(NP)       mean molecular weight
C        WEG          double precision work array for EGLIB
C        IWEG         integer work array for EGLIB
C
C     Output
C     ------
C        PTC(NP)      partial thermal conductivity
C        THETA(NP,NS) thermal diffusion vector
C        D(NP,NS,NS)  flux diffusion matrix
C
C
C     Two CG iterations are performed on the matrix L_[e] in order
C     to evaluate PTC and THETA.
C     Two projected standard iterations are performed on the matrix
C     L_[e] in order to evaluate D.
C
C     Note
C     ----
C        THETA and D satisfy the mass conservation constraints.
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*), IWEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGFLTD1 ( NP, NS, PRES, IPRES, T, WEG(IDLT1), WEG(IDLT2),
     &        WEG(IDLT3), WEG(IDLT4), WEG(IDLT5), WEG(IDLT6),
     &        WEG(IXTR), WEG(IYTR), WEG(IEGWT), WEG(ISUMTR),  
     &        WEG(IEGRU), WEG(IEGPA), PTC, THETA, D,
     &        WEG(IBIN), WEG(IAIJ), WEG(IBIJ), WEG(ICIJ), 
     &        WEG(IFITA), WEG(IFITB), WEG(IFITC),
     &        WEG(ICINT), WEG(IETA), WEG(ICXI), IWEG(IEGLIN),
     &        WEG(IG), WEG(IDMI), WEG(IAN), WEG(IZN), 
     &        WEG(IRN), WEG(ITEMP), WEG(IBETA), WEG(IAUX),
     &        WEG(IAAA), WEG(IBBB) )
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGFLTD1 ( NP, NS, PRES, IPRES, 
     &           TEMPER, DLT1, DLT2, DLT3, DLT4, DLT5, DLT6,  
     &           XTR, YTR, WT, SUMTR, RU, PATMOS, PTC, THETA, D,
     &           BIN, AIJ, BIJ, CIJ, FITA, FITB, FITC, CINT, ETA, CXI,
     &           LIN, G, DMI, AN, ZN, RN, TEMP, BETA, AUX, AAA, BBB )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      CALL EGFEMLE ( NP, NS, TEMPER, DLT1, DLT2, DLT3,  
     &           DLT4, DLT5, DLT6, XTR, WT, RU, PATMOS,
     &           BIN, AIJ, BIJ, CIJ, FITA, FITB, FITC, CINT, ETA, CXI,
     &           LIN, G, BETA )
C-----------------------------------------------------------------------
C        Call the iterative procedures
C-----------------------------------------------------------------------
      NG = 2 * NS
C.....
      ITERMX = 2
      CALL DCOPY  ( np*NG, BETA, 1, RN, 1 )
      CALL EGFCG2 ( NP, NS, NG, G, DMI, AN, ZN, RN, TEMP, 
     &              AAA, BBB, ITERMX )
C.....
      ITERMX = 2
      CALL EGFSI2 ( NP, NS, NG, G, DMI, TEMP, YTR, SUMTR, D, ITERMX, 
     &              AUX, AAA )
C-----------------------------------------------------------------------
C        Evaluate the transport coefficients
C-----------------------------------------------------------------------
      CALL EGFEVL ( NP, NG, PTC, AN, BETA, PATMOS, TEMPER )
      CALL EGFEVT ( NP, NS, THETA, AN, YTR, SUMTR, AAA )
C-----------------------------------------------------------------------
C        Correct for the pressure dependence of theta and D
C-----------------------------------------------------------------------
      IF ( IPRES .EQ. 0 ) THEN
         AA = PATMOS / PRES
         CALL DSCAL ( np*NS, AA, THETA, 1 )
         CALL DSCAL ( np*NS*NS, AA, D, 1 )
      ELSE
         CALL EGFSCATD ( NP, NS, PATMOS, PRES, THETA, D, AAA )
      ENDIF
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGFLTD2 (NP, PRES, IPRES, T, Y, WW, WEG, IWEG, 
     &                    PTC, THETA, D)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 7
C        ITLS  = 2
C
C     Input
C     -----
C        NP           number of nodes
C        PRES(*)      pressure
C        IPRES        flag for pressure (0:scalar, .ne.0:node dependent)
C        T(NP)        temperature
C        Y(NP,NS)     species mass fractions
C        WW(NP)       mean molecular weight
C        WEG          double precision work array for EGLIB
C        IWEG         integer work array for EGLIB
C
C     Output
C     ------
C        PTC(NP)      partial thermal conductivity
C        THETA(NP,NS) thermal diffusion vector
C        D(NP,NS,NS)  flux diffusion matrix
C
C     
C     We form the Choleski decomposition of matrix L_[e] 
C
C     Note
C     ----
C        THETA and D satisfy the mass conservation constraints
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*), IWEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGFLTD2 ( NP, NS, PRES, IPRES, T, WEG(IDLT1), WEG(IDLT2),
     &        WEG(IDLT3), WEG(IDLT4), WEG(IDLT5), WEG(IDLT6),
     &        WEG(IXTR), WEG(IYTR), WEG(IEGWT), WEG(ISUMTR), 
     &        WEG(IEGRU), WEG(IEGPA), PTC, THETA, D,
     &        WEG(IBIN), WEG(IAIJ), WEG(IBIJ), WEG(ICIJ), 
     &        WEG(IFITA), WEG(IFITB), WEG(IFITC),
     &        WEG(ICINT), WEG(IETA), WEG(ICXI), IWEG(IEGLIN),
     &        WEG(IG), WEG(IAN), WEG(IAAA), WEG(ITEMP), WEG(IBETA) )
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGFLTD2 ( NP, NS, PRES, IPRES, 
     &           TEMPER, DLT1, DLT2, DLT3, DLT4, DLT5, DLT6,  
     &           XTR, YTR, WT, SUMTR, RU, PATMOS, PTC, THETA, D,
     &           BIN, AIJ, BIJ, CIJ, FITA, FITB, FITC, CINT, ETA, CXI,
     &           LIN, G, AN, AAA, TEMP, BETA )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      DIMENSION G(NP,*), YTR(NP,NS)
C-----------------------------------------------------------------------
      CALL EGFEMLE ( NP, NS, TEMPER, DLT1, DLT2, DLT3,  
     &           DLT4, DLT5, DLT6, XTR, WT, RU, PATMOS,
     &           BIN, AIJ, BIJ, CIJ, FITA, FITB, FITC, CINT, ETA, CXI,
     &           LIN, G, BETA )
C-----------------------------------------------------------------------
C     Make the matrix G positive definite
C-----------------------------------------------------------------------
      NG = 2 * NS
      DO I = 1, NS
         II1 = NG*(I-1) - (I*(I-1))/2 + I
         do nn = 1, np
            G(nn,II1) = G(nn,II1) + YTR(nn,I) * YTR(nn,I)
         enddo
         DO J = I+1, NS
            IJ1 = NG*(I-1) - (I*(I-1))/2 + J
            do nn = 1, np
               G(nn,IJ1) = G(nn,IJ1) + YTR(nn,I) * YTR(nn,J)
            enddo
         ENDDO
      ENDDO
C-----------------------------------------------------------------------
      CALL EGFDEC ( NP, NG, G, TEMP, IER )
      CALL DCOPY  ( np*NG, BETA, 1, AN, 1 )
      CALL EGFSOL ( NP, NG, G, AN, TEMP )
C-----------------------------------------------------------------------
      CALL EGFEVL ( NP, NG, PTC, AN, BETA, PATMOS, TEMPER )
      CALL EGFEVT ( NP, NS, THETA, AN, YTR, SUMTR, AAA )
C-----------------------------------------------------------------------
      CALL EGFEVD ( NP, NS, NG, AN, YTR, G, SUMTR, D, TEMP, AAA )
C-----------------------------------------------------------------------
C     Project the matrix D
C-----------------------------------------------------------------------
      CALL EGFPDP ( NP, NS, YTR, TEMP, SUMTR, D, AAA )
C-----------------------------------------------------------------------
C        Correct for the pressure dependence of theta and D
C-----------------------------------------------------------------------
      IF ( IPRES .EQ. 0 ) THEN
         AA = PATMOS / PRES
         CALL DSCAL ( np*NS, AA, THETA, 1 )
         CALL DSCAL ( np*NS*NS, AA, D, 1 )
      ELSE
         CALL EGFSCATD ( NP, NS, PATMOS, PRES, THETA, D, AAA )
      ENDIF
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGFLTD3 (NP, PRES, IPRES, T, Y, WW, WEG, IWEG, 
     &                    PTC, THETA, D)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 7
C        ITLS  = 3
C
C     Input
C     -----
C        NP           number of nodes
C        PRES(*)      pressure
C        IPRES        flag for pressure (0:scalar, .ne.0:node dependent)
C        T(NP)        temperature
C        Y(NP,NS)     species mass fractions
C        WW(NP)       mean molecular weight
C        WEG          double precision work array for EGLIB
C        IWEG         integer work array for EGLIB
C
C     Output
C     ------
C        PTC(NP)      partial thermal conductivity
C        THETA(NP,NS) thermal diffusion vector
C        D(NP,NS,NS)  flux diffusion matrix
C     
C
C     Three CG iterations are performed on the matrix L in order
C     to evaluate PTC and THETA.
C     Two projected standard iterations are performed on the matrix
C     L_[00] in order to evaluate D.
C
C     Note
C     ----
C        THETA and D satisfy the mass conservation constraints.
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*), IWEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGFLTD3 ( NP, NS, PRES, IPRES, T, WEG(IDLT1), WEG(IDLT2),
     &         WEG(IDLT3), WEG(IDLT4), WEG(IDLT5), WEG(IDLT6),
     &         WEG(IXTR), WEG(IYTR), WEG(IEGWT), WEG(ISUMTR), 
     &         WEG(IEGRU), WEG(IEGPA), PTC, THETA, D,
     &         WEG(IBIN), WEG(IAIJ), WEG(IBIJ), WEG(ICIJ), 
     &         WEG(IFITA), WEG(IFITB), WEG(IFITC),
     &         WEG(ICINT), WEG(IETA), WEG(ICXI), IWEG(IEGLIN),
     &         WEG(IG), WEG(IDMI), WEG(IAN), WEG(IZN), 
     &         WEG(IRN), WEG(ITEMP), WEG(IBETA), WEG(IAUX),
     &         WEG(IAAA), WEG(IBBB) )
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGFLTD3 ( NP, NS, PRES, IPRES, 
     &           TEMPER, DLT1, DLT2, DLT3, DLT4, DLT5, DLT6,  
     &           XTR, YTR, WT, SUMTR, RU, PATMOS, PTC, THETA, D,
     &           BIN, AIJ, BIJ, CIJ, FITA, FITB, FITC, CINT, ETA, CXI,
     &           LIN, G, DMI, AN, ZN, RN, TEMP, BETA, AUX, AAA, BBB )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      CALL EGFEML ( NP, NS, TEMPER, DLT1, DLT2, DLT3,  
     &           DLT4, DLT5, DLT6, XTR, WT, RU, PATMOS,
     &           BIN, AIJ, BIJ, CIJ, FITA, FITB, FITC, CINT, ETA, CXI,
     &           LIN, G, BETA )
C-----------------------------------------------------------------------
C        Call the iterative procedures
C-----------------------------------------------------------------------
      NG = 3 * NS
C.....
      ITERMX = 3
      CALL DCOPY  ( np*NG, BETA, 1, RN, 1 )
      CALL EGFCG3 ( NP, NS, NG, G, DMI, AN, ZN, RN, TEMP, 
     &              AAA, BBB, ITERMX )
C-----------------------------------------------------------------------
C        Evaluate the transport coefficients
C-----------------------------------------------------------------------
      CALL EGFEVL ( NP, NG, PTC, AN, BETA, PATMOS, TEMPER )
      CALL EGFEVT ( NP, NS, THETA, AN, YTR, SUMTR, AAA )
C-----------------------------------------------------------------------
      CALL EGFEML00 ( NP, NS, XTR, BIN, G )
C.....
      ITERMX = 2
      CALL EGFSI1 ( NP, NS, G, DMI, TEMP, YTR, SUMTR, D, ITERMX, 
     &              AUX, AAA )
C-----------------------------------------------------------------------
C        Correct for the pressure dependence of theta and D
C-----------------------------------------------------------------------
      IF ( IPRES .EQ. 0 ) THEN
         AA = PATMOS / PRES
         CALL DSCAL ( np*NS, AA, THETA, 1 )
         CALL DSCAL ( np*NS*NS, AA, D, 1 )
      ELSE
         CALL EGFSCATD ( NP, NS, PATMOS, PRES, THETA, D, AAA )
      ENDIF
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGFLTD4 ( NP, PRES, IPRES, T, Y, WW, WEG, IWEG, 
     &                     PTC, THETA, D )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 7
C        ITLS  = 3
C
C     Input
C     -----
C        NP           number of nodes
C        PRES(*)      pressure
C        IPRES        flag for pressure (0:scalar, .ne.0:node dependent)
C        T(NP)        temperature
C        Y(NP,NS)     species mass fractions
C        WW(NP)       mean molecular weight
C        WEG          double precision work array for EGLIB
C        IWEG         integer work array for EGLIB
C
C     Output
C     ------
C        PTC(NP)      partial thermal conductivity
C        THETA(NP,NS) thermal diffusion vector
C        D(NP,NS,NS)  flux diffusion matrix
C
C     
C     We form the Choleski decomposition of matrices L and L_[00] 
C
C     Note
C     ----
C        THETA and D satisfy the mass conservation constraints
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*), IWEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGFLTD4 ( NP, NS, PRES, IPRES, T, WEG(IDLT1), WEG(IDLT2),
     &         WEG(IDLT3), WEG(IDLT4), WEG(IDLT5), WEG(IDLT6),
     &         WEG(IXTR), WEG(IYTR), WEG(IEGWT), WEG(ISUMTR), 
     &         WEG(IEGRU), WEG(IEGPA), PTC, THETA, D,
     &         WEG(IBIN), WEG(IAIJ), WEG(IBIJ), WEG(ICIJ), 
     &         WEG(IFITA), WEG(IFITB), WEG(IFITC),
     &         WEG(ICINT), WEG(IETA), WEG(ICXI), IWEG(IEGLIN),
     &         WEG(IG), WEG(IAN), WEG(IAAA), WEG(ITEMP), WEG(IBETA) )
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGFLTD4 ( NP, NS, PRES, IPRES, 
     &           TEMPER, DLT1, DLT2, DLT3, DLT4, DLT5, DLT6,  
     &           XTR, YTR, WT, SUMTR, RU, PATMOS, PTC, THETA, D,
     &           BIN, AIJ, BIJ, CIJ, FITA, FITB, FITC, CINT, ETA, CXI,
     &           LIN, G, AN, AAA, TEMP, BETA )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      DIMENSION G(NP,*), YTR(NP,NS)
C-----------------------------------------------------------------------
      CALL EGFEML ( NP, NS, TEMPER, DLT1, DLT2, DLT3,  
     &              DLT4, DLT5, DLT6, XTR, WT, RU, PATMOS,
     &              BIN, AIJ, BIJ, CIJ, FITA, FITB, FITC, 
     &              CINT, ETA, CXI, LIN, G, BETA )
C-----------------------------------------------------------------------
C     Make the matrix G positive definite
C-----------------------------------------------------------------------
      NG = 3 * NS
      DO I = 1, NS
         II1 = NG*(I-1) - (I*(I-1))/2 + I
         do nn = 1, np
            G(nn,II1) = G(nn,II1) + YTR(nn,I) * YTR(nn,I)
         enddo
         DO J = I+1, NS
            IJ1 = NG*(I-1) - (I*(I-1))/2 + J
            do nn = 1, np
               G(nn,IJ1) = G(nn,IJ1) + YTR(nn,I) * YTR(nn,J)
            enddo
         ENDDO
      ENDDO
C-----------------------------------------------------------------------
      CALL EGFDEC ( NP, NG, G, TEMP, IER )
      CALL DCOPY  ( np*NG, BETA, 1, AN, 1 )
      CALL EGFSOL ( NP, NG, G, AN, TEMP )
C-----------------------------------------------------------------------
      CALL EGFEVL ( NP, NG, PTC, AN, BETA, PATMOS, TEMPER )
      CALL EGFEVT ( NP, NS, THETA, AN, YTR, SUMTR, AAA )
C-----------------------------------------------------------------------
      CALL EGFEML00 ( NP, NS, XTR, BIN, G )
C-----------------------------------------------------------------------
C     Make the matrix G positive definite
C-----------------------------------------------------------------------
      NG = NS
      DO I = 1, NS
         II1 = NG*(I-1) - (I*(I-1))/2 + I
         do nn = 1, np
            G(nn,II1) = G(nn,II1) + YTR(nn,I) * YTR(nn,I)
         enddo
         DO J = I+1, NS
            IJ1 = NG*(I-1) - (I*(I-1))/2 + J
            do nn = 1, np
               G(nn,IJ1) = G(nn,IJ1) + YTR(nn,I) * YTR(nn,J)
            enddo
         ENDDO
      ENDDO
C-----------------------------------------------------------------------
      CALL EGFDEC ( NP, NG, G, TEMP, IER )
C-----------------------------------------------------------------------
      CALL EGFEVD ( NP, NS, NG, AN, YTR, G, SUMTR, D, TEMP, AAA )
C-----------------------------------------------------------------------
C     Project the matrix D
C-----------------------------------------------------------------------
      CALL EGFPDP ( NP, NS, YTR, TEMP, SUMTR, D, AAA )
C-----------------------------------------------------------------------
C        Correct for the pressure dependence of theta and D
C-----------------------------------------------------------------------
      IF ( IPRES .EQ. 0 ) THEN
         AA = PATMOS / PRES
         CALL DSCAL ( np*NS, AA, THETA, 1 )
         CALL DSCAL ( np*NS*NS, AA, D, 1 )
      ELSE
         CALL EGFSCATD ( NP, NS, PATMOS, PRES, THETA, D, AAA )
      ENDIF
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGFLTD5 (NP, PRES, IPRES, T, Y, WW, WEG, IWEG, 
     &                    PTC, THETA, D)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 7
C        ITLS  = 3
C
C     Input
C     -----
C        NP           number of nodes
C        PRES(*)      pressure
C        IPRES        flag for pressure (0:scalar, .ne.0:node dependent)
C        T(NP)        temperature
C        Y(NP,NS)     species mass fractions
C        WW(NP)       mean molecular weight
C        WEG          double precision work array for EGLIB
C        IWEG         integer work array for EGLIB
C
C     Output
C     ------
C        PTC(NP)      partial thermal conductivity
C        THETA(NP,NS) thermal diffusion vector
C        D(NP,NS,NS)  flux diffusion matrix
C     
C
C     Three CG iterations are performed on the matrix L in order
C     to evaluate PTC and THETA.
C     Two projected standard iterations are performed on the matrix
C     L in order to evaluate D.
C
C     Note
C     ----
C        THETA and D satisfy the mass conservation constraints.
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*), IWEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGFLTD5 ( NP, NS, PRES, IPRES, T, WEG(IDLT1), WEG(IDLT2),
     &         WEG(IDLT3), WEG(IDLT4), WEG(IDLT5), WEG(IDLT6),
     &         WEG(IXTR), WEG(IYTR), WEG(IEGWT), WEG(ISUMTR), 
     &         WEG(IEGRU), WEG(IEGPA), PTC, THETA, D,
     &         WEG(IBIN), WEG(IAIJ), WEG(IBIJ), WEG(ICIJ), 
     &         WEG(IFITA), WEG(IFITB), WEG(IFITC),
     &         WEG(ICINT), WEG(IETA), WEG(ICXI), IWEG(IEGLIN),
     &         WEG(IG), WEG(IDMI), WEG(IAN), WEG(IZN), WEG(IRN), 
     &         WEG(ITEMP), WEG(IBETA), WEG(IAUX), WEG(IAAA), WEG(IBBB) )
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGFLTD5 ( NP, NS, PRES, IPRES, 
     &           TEMPER, DLT1, DLT2, DLT3, DLT4, DLT5, DLT6,  
     &           XTR, YTR, WT, SUMTR, RU, PATMOS, PTC, THETA, D,
     &           BIN, AIJ, BIJ, CIJ, FITA, FITB, FITC, CINT, ETA, CXI,
     &           LIN, G, DMI, AN, ZN, RN, TEMP, BETA, AUX, AAA, BBB )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      CALL EGFEML ( NP, NS, TEMPER, DLT1, DLT2, DLT3,  
     &           DLT4, DLT5, DLT6, XTR, WT, RU, PATMOS,
     &           BIN, AIJ, BIJ, CIJ, FITA, FITB, FITC, CINT, ETA, CXI,
     &           LIN, G, BETA )
C-----------------------------------------------------------------------
C        Call the iterative procedures
C-----------------------------------------------------------------------
      NG = 3 * NS
C.....
      ITERMX = 3
      CALL DCOPY  ( np*NG, BETA, 1, RN, 1 )
      CALL EGFCG3 ( NP, NS, NG, G, DMI, AN, ZN, RN, TEMP, 
     &              AAA, BBB, ITERMX )
C.....
      ITERMX = 2
      CALL EGFSI3 ( NP, NS, NG, G, DMI, TEMP, YTR, SUMTR, D, ITERMX, 
     &              AUX, AAA )
C-----------------------------------------------------------------------
C        Evaluate the transport coefficients
C-----------------------------------------------------------------------
      CALL EGFEVL ( NP, NG, PTC, AN, BETA, PATMOS, TEMPER )
      CALL EGFEVT ( NP, NS, THETA, AN, YTR, SUMTR, AAA )
C-----------------------------------------------------------------------
C        Correct for the pressure dependence of theta and D
C-----------------------------------------------------------------------
      IF ( IPRES .EQ. 0 ) THEN
         AA = PATMOS / PRES
         CALL DSCAL ( np*NS, AA, THETA, 1 )
         CALL DSCAL ( np*NS*NS, AA, D, 1 )
      ELSE
         CALL EGFSCATD ( NP, NS, PATMOS, PRES, THETA, D, AAA )
      ENDIF
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGFLTD6 ( NP, PRES, IPRES, T, Y, WW, WEG, IWEG, 
     &                     PTC, THETA, D )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 7
C        ITLS  = 3
C
C     Input
C     -----
C        NP           number of nodes
C        PRES(*)      pressure
C        IPRES        flag for pressure (0:scalar, .ne.0:node dependent)
C        T(NP)        temperature
C        Y(NP,NS)     species mass fractions
C        WW(NP)       mean molecular weight
C        WEG          double precision work array for EGLIB
C        IWEG         integer work array for EGLIB
C
C     Output
C     ------
C        PTC(NP)      partial thermal conductivity
C        THETA(NP,NS) thermal diffusion vector
C        D(NP,NS,NS)  flux diffusion matrix
C
C     
C     We form the Choleski decomposition of matrix L
C
C     Note
C     ----
C        THETA and D satisfy the mass conservation constraints
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*), IWEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGFLTD6 ( NP, NS, PRES, IPRES, T, WEG(IDLT1), WEG(IDLT2),
     &         WEG(IDLT3), WEG(IDLT4), WEG(IDLT5), WEG(IDLT6),
     &         WEG(IXTR), WEG(IYTR), WEG(IEGWT), WEG(ISUMTR), 
     &         WEG(IEGRU), WEG(IEGPA), PTC, THETA, D,
     &         WEG(IBIN), WEG(IAIJ), WEG(IBIJ), WEG(ICIJ), 
     &         WEG(IFITA), WEG(IFITB), WEG(IFITC),
     &         WEG(ICINT), WEG(IETA), WEG(ICXI), IWEG(IEGLIN),
     &         WEG(IG), WEG(IAN), WEG(IAAA), WEG(ITEMP), WEG(IBETA) )
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGFLTD6 ( NP, NS, PRES, IPRES, 
     &           TEMPER, DLT1, DLT2, DLT3, DLT4, DLT5, DLT6,  
     &           XTR, YTR, WT, SUMTR, RU, PATMOS, PTC, THETA, D,
     &           BIN, AIJ, BIJ, CIJ, FITA, FITB, FITC, CINT, ETA, CXI,
     &           LIN, G, AN, AAA, TEMP, BETA )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      DIMENSION G(NP,*), YTR(NP,NS)
C-----------------------------------------------------------------------
      CALL EGFEML ( NP, NS, TEMPER, DLT1, DLT2, DLT3,  
     &              DLT4, DLT5, DLT6, XTR, WT, RU, PATMOS,
     &              BIN, AIJ, BIJ, CIJ, FITA, FITB, FITC, 
     &              CINT, ETA, CXI, LIN, G, BETA )
C-----------------------------------------------------------------------
C     Make the matrix G positive definite
C-----------------------------------------------------------------------
      NG = 3 * NS
      DO I = 1, NS
         II1 = NG*(I-1) - (I*(I-1))/2 + I
         do nn = 1, np
            G(nn,II1) = G(nn,II1) + YTR(nn,I) * YTR(nn,I)
         enddo
         DO J = I+1, NS
            IJ1 = NG*(I-1) - (I*(I-1))/2 + J
            do nn = 1, np
               G(nn,IJ1) = G(nn,IJ1) + YTR(nn,I) * YTR(nn,J)
            enddo
         ENDDO
      ENDDO
C-----------------------------------------------------------------------
      CALL EGFDEC ( NP, NG, G, TEMP, IER )
      CALL DCOPY  ( np*NG, BETA, 1, AN, 1 )
      CALL EGFSOL ( NP, NG, G, AN, TEMP )
C-----------------------------------------------------------------------
      CALL EGFEVL ( NP, NG, PTC, AN, BETA, PATMOS, TEMPER )
      CALL EGFEVT ( NP, NS, THETA, AN, YTR, SUMTR, AAA )
C-----------------------------------------------------------------------
      CALL EGFEVD ( NP, NS, NG, AN, YTR, G, SUMTR, D, TEMP, AAA )
C-----------------------------------------------------------------------
C     Project the matrix D
C-----------------------------------------------------------------------
      CALL EGFPDP ( NP, NS, YTR, TEMP, SUMTR, D, AAA )
C-----------------------------------------------------------------------
C        Correct for the pressure dependence of theta and D
C-----------------------------------------------------------------------
      IF ( IPRES .EQ. 0 ) THEN
         AA = PATMOS / PRES
         CALL DSCAL ( np*NS, AA, THETA, 1 )
         CALL DSCAL ( np*NS*NS, AA, D, 1 )
      ELSE
         CALL EGFSCATD ( NP, NS, PATMOS, PRES, THETA, D, AAA )
      ENDIF
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGFLTDR1 (NP, T, Y, WW, WEG, IWEG, PTC, THETA, D)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 7
C        ITLS  = 2
C
C     Input
C     -----
C        NP           number of nodes
C        T(NP)        temperature
C        Y(NP,NS)     species mass fractions
C        WW(NP)       mean molecular weight
C        WEG          double precision work array for EGLIB
C        IWEG         integer work array for EGLIB
C
C     Output
C     ------
C        PTC(NP)      partial thermal conductivity
C        THETA(NP,NS) rho * thermal diffusion vector
C        D(NP,NS,NS)  rho * flux diffusion matrix
C
C        where rho is the density.
C
C     Two CG iterations are performed on the matrix L_[e] in order
C     to evaluate PTC and THETA.
C     Two projected standard iterations are performed on the matrix
C     L_[e] in order to evaluate D.
C
C     Note
C     ----
C        THETA and D satisfy the mass conservation constraints.
C        These transport coefficients are also PRESSURE INDEPENDENT.
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*), IWEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGFLTDR1 ( NP, NS, T, WEG(IDLT1), WEG(IDLT2),
     &        WEG(IDLT3), WEG(IDLT4), WEG(IDLT5), WEG(IDLT6),
     &        WEG(IXTR), WEG(IYTR), WEG(IEGWT), WEG(IWWTR), WEG(ISUMTR), 
     &        WEG(IEGRU), WEG(IEGPA), PTC, THETA, D,
     &        WEG(IBIN), WEG(IAIJ), WEG(IBIJ), WEG(ICIJ), 
     &        WEG(IFITA), WEG(IFITB), WEG(IFITC),
     &        WEG(ICINT), WEG(IETA), WEG(ICXI), IWEG(IEGLIN),
     &        WEG(IG), WEG(IDMI), WEG(IAN), WEG(IZN), 
     &        WEG(IRN), WEG(ITEMP), WEG(IBETA), WEG(IAUX),
     &        WEG(IAAA), WEG(IBBB) )
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGFLTDR1 ( NP, NS, TEMPER, DLT1, DLT2, DLT3,  
     &           DLT4, DLT5, DLT6, XTR, YTR, WT, WWTR, SUMTR,
     &           RU, PATMOS, PTC, THETA, D,
     &           BIN, AIJ, BIJ, CIJ, FITA, FITB, FITC, CINT, ETA, CXI,
     &           LIN, G, DMI, AN, ZN, RN, TEMP, BETA, AUX, AAA, BBB )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      CALL EGFEMLE ( NP, NS, TEMPER, DLT1, DLT2, DLT3,  
     &           DLT4, DLT5, DLT6, XTR, WT, RU, PATMOS,
     &           BIN, AIJ, BIJ, CIJ, FITA, FITB, FITC, CINT, ETA, CXI,
     &           LIN, G, BETA )
C-----------------------------------------------------------------------
C        Call the iterative procedures
C-----------------------------------------------------------------------
      NG = 2 * NS
C.....
      ITERMX = 2
      CALL DCOPY  ( np*NG, BETA, 1, RN, 1 )
      CALL EGFCG2 ( NP, NS, NG, G, DMI, AN, ZN, RN, TEMP, 
     &              AAA, BBB, ITERMX )
C.....
      ITERMX = 2
      CALL EGFSI2 ( NP, NS, NG, G, DMI, TEMP, YTR, SUMTR, D, ITERMX, 
     &              AUX, AAA )
C-----------------------------------------------------------------------
C        Evaluate the transport coefficients
C-----------------------------------------------------------------------
      CALL EGFEVL ( NP, NG, PTC, AN, BETA, PATMOS, TEMPER )
      CALL EGFEVT ( NP, NS, THETA, AN, YTR, SUMTR, AAA )
C-----------------------------------------------------------------------
      CALL EGFSCATDR ( NP, NS, PATMOS, WWTR, RU, TEMPER, THETA, D, AAA )
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGFLTDR2 (NP, T, Y, WW, WEG, IWEG, PTC, THETA, D)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 7
C        ITLS  = 2
C
C     Input
C     -----
C        NP           number of nodes
C        T(NP)        temperature
C        Y(NP,NS)     species mass fractions
C        WW(NP)       mean molecular weight
C        WEG          double precision work array for EGLIB
C        IWEG         integer work array for EGLIB
C
C     Output
C     ------
C        PTC(NP)      partial thermal conductivity
C        THETA(NP,NS) rho * thermal diffusion vector
C        D(NP,NS,NS)  rho * flux diffusion matrix
C
C        where rho is the density.
C
C     
C     We form the Choleski decomposition of matrix L_[e] 
C
C     Note
C     ----
C        THETA and D satisfy the mass conservation constraints
C        These transport coefficients are also PRESSURE INDEPENDENT.
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*), IWEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGFLTDR2 ( NP, NS, T, WEG(IDLT1), WEG(IDLT2),
     &        WEG(IDLT3), WEG(IDLT4), WEG(IDLT5), WEG(IDLT6),
     &        WEG(IXTR), WEG(IYTR), WEG(IEGWT), WEG(IWWTR), WEG(ISUMTR), 
     &        WEG(IEGRU), WEG(IEGPA), PTC, THETA, D,
     &        WEG(IBIN), WEG(IAIJ), WEG(IBIJ), WEG(ICIJ), 
     &        WEG(IFITA), WEG(IFITB), WEG(IFITC),
     &        WEG(ICINT), WEG(IETA), WEG(ICXI), IWEG(IEGLIN),
     &        WEG(IG), WEG(IAN), WEG(IAAA), WEG(ITEMP), WEG(IBETA) )
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGFLTDR2 ( NP, NS, TEMPER, DLT1, DLT2, DLT3,  
     &           DLT4, DLT5, DLT6, XTR, YTR, WT, WWTR, SUMTR,
     &           RU, PATMOS, PTC, THETA, D,
     &           BIN, AIJ, BIJ, CIJ, FITA, FITB, FITC, CINT, ETA, CXI,
     &           LIN, G, AN, AAA, TEMP, BETA )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      DIMENSION G(NP,*), YTR(NP,NS)
C-----------------------------------------------------------------------
      CALL EGFEMLE ( NP, NS, TEMPER, DLT1, DLT2, DLT3,  
     &           DLT4, DLT5, DLT6, XTR, WT, RU, PATMOS,
     &           BIN, AIJ, BIJ, CIJ, FITA, FITB, FITC, CINT, ETA, CXI,
     &           LIN, G, BETA )
C-----------------------------------------------------------------------
C     Make the matrix G positive definite
C-----------------------------------------------------------------------
      NG = 2 * NS
      DO I = 1, NS
         II1 = NG*(I-1) - (I*(I-1))/2 + I
         do nn = 1, np
            G(nn,II1) = G(nn,II1) + YTR(nn,I) * YTR(nn,I)
         enddo
         DO J = I+1, NS
            IJ1 = NG*(I-1) - (I*(I-1))/2 + J
            do nn = 1, np
               G(nn,IJ1) = G(nn,IJ1) + YTR(nn,I) * YTR(nn,J)
            enddo
         ENDDO
      ENDDO
C-----------------------------------------------------------------------
      CALL EGFDEC ( NP, NG, G, TEMP, IER )
      CALL DCOPY  ( np*NG, BETA, 1, AN, 1 )
      CALL EGFSOL ( NP, NG, G, AN, TEMP )
C-----------------------------------------------------------------------
      CALL EGFEVL ( NP, NG, PTC, AN, BETA, PATMOS, TEMPER )
      CALL EGFEVT ( NP, NS, THETA, AN, YTR, SUMTR, AAA )
C-----------------------------------------------------------------------
      CALL EGFEVD ( NP, NS, NG, AN, YTR, G, SUMTR, D, TEMP, AAA )
C-----------------------------------------------------------------------
C     Project the matrix D
C-----------------------------------------------------------------------
      CALL EGFPDP ( NP, NS, YTR, TEMP, SUMTR, D, AAA )
C-----------------------------------------------------------------------
      CALL EGFSCATDR ( NP, NS, PATMOS, WWTR, RU, TEMPER, THETA, D, AAA )
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGFLTDR3 (NP, T, Y, WW, WEG, IWEG, PTC, THETA, D)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 7
C        ITLS  = 3
C
C     Input
C     -----
C        NP           number of nodes
C        T(NP)        temperature
C        Y(NP,NS)     species mass fractions
C        WW(NP)       mean molecular weight
C        WEG          double precision work array for EGLIB
C        IWEG         integer work array for EGLIB
C
C     Output
C     ------
C        PTC(NP)      partial thermal conductivity
C        THETA(NP,NS) rho * thermal diffusion vector
C        D(NP,NS,NS)  rho * flux diffusion matrix
C
C        where rho is the density.
C     
C
C     Three CG iterations are performed on the matrix L in order
C     to evaluate PTC and THETA.
C     Two projected standard iterations are performed on the matrix
C     L_[00] in order to evaluate D.
C
C     Note
C     ----
C        THETA and D satisfy the mass conservation constraints.
C        These transport coefficients are also PRESSURE INDEPENDENT.
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*), IWEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGFLTDR3 ( NP, NS, T, WEG(IDLT1), WEG(IDLT2),
     &        WEG(IDLT3), WEG(IDLT4), WEG(IDLT5), WEG(IDLT6),
     &        WEG(IXTR), WEG(IYTR), WEG(IEGWT), WEG(IWWTR), WEG(ISUMTR), 
     &        WEG(IEGRU), WEG(IEGPA), PTC, THETA, D,
     &        WEG(IBIN), WEG(IAIJ), WEG(IBIJ), WEG(ICIJ), 
     &        WEG(IFITA), WEG(IFITB), WEG(IFITC),
     &        WEG(ICINT), WEG(IETA), WEG(ICXI), IWEG(IEGLIN),
     &        WEG(IG), WEG(IDMI), WEG(IAN), WEG(IZN), 
     &        WEG(IRN), WEG(ITEMP), WEG(IBETA), WEG(IAUX),
     &        WEG(IAAA), WEG(IBBB) )
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGFLTDR3 ( NP, NS, TEMPER, DLT1, DLT2, DLT3,  
     &           DLT4, DLT5, DLT6, XTR, YTR, WT, WWTR, SUMTR,
     &           RU, PATMOS, PTC, THETA, D,
     &           BIN, AIJ, BIJ, CIJ, FITA, FITB, FITC, CINT, ETA, CXI,
     &           LIN, G, DMI, AN, ZN, RN, TEMP, BETA, AUX, AAA, BBB )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      CALL EGFEML ( NP, NS, TEMPER, DLT1, DLT2, DLT3,  
     &           DLT4, DLT5, DLT6, XTR, WT, RU, PATMOS,
     &           BIN, AIJ, BIJ, CIJ, FITA, FITB, FITC, CINT, ETA, CXI,
     &           LIN, G, BETA )
C-----------------------------------------------------------------------
C        Call the iterative procedures
C-----------------------------------------------------------------------
      NG = 3 * NS
C.....
      ITERMX = 3
      CALL DCOPY  ( np*NG, BETA, 1, RN, 1 )
      CALL EGFCG3 ( NP, NS, NG, G, DMI, AN, ZN, RN, TEMP, 
     &              AAA, BBB, ITERMX )
C-----------------------------------------------------------------------
C        Evaluate the transport coefficients
C-----------------------------------------------------------------------
      CALL EGFEVL ( NP, NG, PTC, AN, BETA, PATMOS, TEMPER )
      CALL EGFEVT ( NP, NS, THETA, AN, YTR, SUMTR, AAA )
C-----------------------------------------------------------------------
      CALL EGFEML00 ( NP, NS, XTR, BIN, G )
C.....
      ITERMX = 2
      CALL EGFSI1 ( NP, NS, G, DMI, TEMP, YTR, SUMTR, D, ITERMX, 
     &              AUX, AAA )
C-----------------------------------------------------------------------
      CALL EGFSCATDR ( NP, NS, PATMOS, WWTR, RU, TEMPER, THETA, D, AAA )
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGFLTDR4 ( NP, T, Y, WW, WEG, IWEG, PTC, THETA, D )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 7
C        ITLS  = 3
C
C     Input
C     -----
C        NP           number of nodes
C        T(NP)        temperature
C        Y(NP,NS)     species mass fractions
C        WW(NP)       mean molecular weight
C        WEG          double precision work array for EGLIB
C        IWEG         integer work array for EGLIB
C
C     Output
C     ------
C        PTC(NP)      partial thermal conductivity
C        THETA(NP,NS) rho * thermal diffusion vector
C        D(NP,NS,NS)  rho * flux diffusion matrix
C
C        where rho is the density.
C
C     
C     We form the Choleski decomposition of matrices L and L_[00] 
C
C     Note
C     ----
C        THETA and D satisfy the mass conservation constraints
C        These transport coefficients are also PRESSURE INDEPENDENT.
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*), IWEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGFLTDR4 ( NP, NS, T, WEG(IDLT1), WEG(IDLT2),
     &        WEG(IDLT3), WEG(IDLT4), WEG(IDLT5), WEG(IDLT6),
     &        WEG(IXTR), WEG(IYTR), WEG(IEGWT), WEG(IWWTR), WEG(ISUMTR),  
     &        WEG(IEGRU), WEG(IEGPA), PTC, THETA, D,
     &        WEG(IBIN), WEG(IAIJ), WEG(IBIJ), WEG(ICIJ), 
     &        WEG(IFITA), WEG(IFITB), WEG(IFITC),
     &        WEG(ICINT), WEG(IETA), WEG(ICXI), IWEG(IEGLIN),
     &        WEG(IG), WEG(IAN), WEG(IAAA), WEG(ITEMP), WEG(IBETA) )
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGFLTDR4 ( NP, NS, TEMPER, DLT1, DLT2, DLT3,  
     &           DLT4, DLT5, DLT6, XTR, YTR, WT, WWTR, SUMTR,
     &           RU, PATMOS, PTC, THETA, D,
     &           BIN, AIJ, BIJ, CIJ, FITA, FITB, FITC, CINT, ETA, CXI,
     &           LIN, G, AN, AAA, TEMP, BETA )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      DIMENSION G(NP,*), YTR(NP,NS)
C-----------------------------------------------------------------------
      CALL EGFEML ( NP, NS, TEMPER, DLT1, DLT2, DLT3,  
     &              DLT4, DLT5, DLT6, XTR, WT, RU, PATMOS,
     &              BIN, AIJ, BIJ, CIJ, FITA, FITB, FITC, 
     &              CINT, ETA, CXI, LIN, G, BETA )
C-----------------------------------------------------------------------
C     Make the matrix G positive definite
C-----------------------------------------------------------------------
      NG = 3 * NS
      DO I = 1, NS
         II1 = NG*(I-1) - (I*(I-1))/2 + I
         do nn = 1, np
            G(nn,II1) = G(nn,II1) + YTR(nn,I) * YTR(nn,I)
         enddo
         DO J = I+1, NS
            IJ1 = NG*(I-1) - (I*(I-1))/2 + J
            do nn = 1, np
               G(nn,IJ1) = G(nn,IJ1) + YTR(nn,I) * YTR(nn,J)
            enddo
         ENDDO
      ENDDO
C-----------------------------------------------------------------------
      CALL EGFDEC ( NP, NG, G, TEMP, IER )
      CALL DCOPY  ( np*NG, BETA, 1, AN, 1 )
      CALL EGFSOL ( NP, NG, G, AN, TEMP )
C-----------------------------------------------------------------------
      CALL EGFEVL ( NP, NG, PTC, AN, BETA, PATMOS, TEMPER )
      CALL EGFEVT ( NP, NS, THETA, AN, YTR, SUMTR, AAA )
C-----------------------------------------------------------------------
      CALL EGFEML00 ( NP, NS, XTR, BIN, G )
C-----------------------------------------------------------------------
C     Make the matrix G positive definite
C-----------------------------------------------------------------------
      NG = NS
      DO I = 1, NS
         II1 = NG*(I-1) - (I*(I-1))/2 + I
         do nn = 1, np
            G(nn,II1) = G(nn,II1) + YTR(nn,I) * YTR(nn,I)
         enddo
         DO J = I+1, NS
            IJ1 = NG*(I-1) - (I*(I-1))/2 + J
            do nn = 1, np
               G(nn,IJ1) = G(nn,IJ1) + YTR(nn,I) * YTR(nn,J)
            enddo
         ENDDO
      ENDDO
C-----------------------------------------------------------------------
      CALL EGFDEC ( NP, NG, G, TEMP, IER )
C-----------------------------------------------------------------------
      CALL EGFEVD ( NP, NS, NG, AN, YTR, G, SUMTR, D, TEMP, AAA )
C-----------------------------------------------------------------------
C     Project the matrix D
C-----------------------------------------------------------------------
      CALL EGFPDP ( NP, NS, YTR, TEMP, SUMTR, D, AAA )
C-----------------------------------------------------------------------
      CALL EGFSCATDR ( NP, NS, PATMOS, WWTR, RU, TEMPER, THETA, D, AAA )
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGFLTDR5 (NP, T, Y, WW, WEG, IWEG, PTC, THETA, D)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 7
C        ITLS  = 3
C
C     Input
C     -----
C        NP           number of nodes
C        T(NP)        temperature
C        Y(NP,NS)     species mass fractions
C        WW(NP)       mean molecular weight
C        WEG          double precision work array for EGLIB
C        IWEG         integer work array for EGLIB
C
C     Output
C     ------
C        PTC(NP)      partial thermal conductivity
C        THETA(NP,NS) rho * thermal diffusion vector
C        D(NP,NS,NS)  rho * flux diffusion matrix
C
C        where rho is the density.
C     
C
C     Three CG iterations are performed on the matrix L in order
C     to evaluate PTC and THETA.
C     Two projected standard iterations are performed on the matrix
C     L in order to evaluate D.
C
C     Note
C     ----
C        THETA and D satisfy the mass conservation constraints.
C        These transport coefficients are also PRESSURE INDEPENDENT.
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*), IWEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGFLTDR5 ( NP, NS, T, WEG(IDLT1), WEG(IDLT2),
     &        WEG(IDLT3), WEG(IDLT4), WEG(IDLT5), WEG(IDLT6),
     &        WEG(IXTR), WEG(IYTR), WEG(IEGWT), WEG(IWWTR), WEG(ISUMTR), 
     &        WEG(IEGRU), WEG(IEGPA), PTC, THETA, D,
     &        WEG(IBIN), WEG(IAIJ), WEG(IBIJ), WEG(ICIJ), 
     &        WEG(IFITA), WEG(IFITB), WEG(IFITC),
     &        WEG(ICINT), WEG(IETA), WEG(ICXI), IWEG(IEGLIN),
     &        WEG(IG), WEG(IDMI), WEG(IAN), WEG(IZN), 
     &        WEG(IRN), WEG(ITEMP), WEG(IBETA), WEG(IAUX),
     &        WEG(IAAA), WEG(IBBB) )
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGFLTDR5 ( NP, NS, TEMPER, DLT1, DLT2, DLT3,  
     &           DLT4, DLT5, DLT6, XTR, YTR, WT, WWTR, SUMTR,
     &           RU, PATMOS, PTC, THETA, D,
     &           BIN, AIJ, BIJ, CIJ, FITA, FITB, FITC, CINT, ETA, CXI,
     &           LIN, G, DMI, AN, ZN, RN, TEMP, BETA, AUX, AAA, BBB )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      CALL EGFEML ( NP, NS, TEMPER, DLT1, DLT2, DLT3,  
     &           DLT4, DLT5, DLT6, XTR, WT, RU, PATMOS,
     &           BIN, AIJ, BIJ, CIJ, FITA, FITB, FITC, CINT, ETA, CXI,
     &           LIN, G, BETA )
C-----------------------------------------------------------------------
C        Call the iterative procedures
C-----------------------------------------------------------------------
      NG = 3 * NS
C.....
      ITERMX = 3
      CALL DCOPY  ( np*NG, BETA, 1, RN, 1 )
      CALL EGFCG3 ( NP, NS, NG, G, DMI, AN, ZN, RN, TEMP, 
     &              AAA, BBB, ITERMX )
C.....
      ITERMX = 2
      CALL EGFSI3 ( NP, NS, NG, G, DMI, TEMP, YTR, SUMTR, D, ITERMX, 
     &              AUX, AAA )
C-----------------------------------------------------------------------
C        Evaluate the transport coefficients
C-----------------------------------------------------------------------
      CALL EGFEVL ( NP, NG, PTC, AN, BETA, PATMOS, TEMPER )
      CALL EGFEVT ( NP, NS, THETA, AN, YTR, SUMTR, AAA )
C-----------------------------------------------------------------------
      CALL EGFSCATDR ( NP, NS, PATMOS, WWTR, RU, TEMPER, THETA, D, AAA )
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGFLTDR6 ( NP, T, Y, WW, WEG, IWEG, PTC, THETA, D )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 7
C        ITLS  = 3
C
C     Input
C     -----
C        NP           number of nodes
C        T(NP)        temperature
C        Y(NP,NS)     species mass fractions
C        WW(NP)       mean molecular weight
C        WEG          double precision work array for EGLIB
C        IWEG         integer work array for EGLIB
C
C     Output
C     ------
C        PTC(NP)      partial thermal conductivity
C        THETA(NP,NS) rho * thermal diffusion vector
C        D(NP,NS,NS)  rho * flux diffusion matrix
C
C        where rho is the density.
C
C     
C     We form the Choleski decomposition of matrix L
C
C     Note
C     ----
C        THETA and D satisfy the mass conservation constraints
C        These transport coefficients are also PRESSURE INDEPENDENT.
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*), IWEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGFLTDR6 ( NP, NS, T, WEG(IDLT1), WEG(IDLT2),
     &        WEG(IDLT3), WEG(IDLT4), WEG(IDLT5), WEG(IDLT6),
     &        WEG(IXTR), WEG(IYTR), WEG(IEGWT), WEG(IWWTR), WEG(ISUMTR),  
     &        WEG(IEGRU), WEG(IEGPA), PTC, THETA, D,
     &        WEG(IBIN), WEG(IAIJ), WEG(IBIJ), WEG(ICIJ), 
     &        WEG(IFITA), WEG(IFITB), WEG(IFITC),
     &        WEG(ICINT), WEG(IETA), WEG(ICXI), IWEG(IEGLIN),
     &        WEG(IG), WEG(IAN), WEG(IAAA), WEG(ITEMP), WEG(IBETA) )
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGFLTDR6 ( NP, NS, TEMPER, DLT1, DLT2, DLT3,  
     &           DLT4, DLT5, DLT6, XTR, YTR, WT, WWTR, SUMTR, 
     &           RU, PATMOS, PTC, THETA, D,
     &           BIN, AIJ, BIJ, CIJ, FITA, FITB, FITC, CINT, ETA, CXI,
     &           LIN, G, AN, AAA, TEMP, BETA )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      DIMENSION G(NP,*), YTR(NP,NS)
C-----------------------------------------------------------------------
      CALL EGFEML ( NP, NS, TEMPER, DLT1, DLT2, DLT3,  
     &              DLT4, DLT5, DLT6, XTR, WT, RU, PATMOS,
     &              BIN, AIJ, BIJ, CIJ, FITA, FITB, FITC, 
     &              CINT, ETA, CXI, LIN, G, BETA )
C-----------------------------------------------------------------------
C     Make the matrix G positive definite
C-----------------------------------------------------------------------
      NG = 3 * NS
      DO I = 1, NS
         II1 = NG*(I-1) - (I*(I-1))/2 + I
         do nn = 1, np
            G(nn,II1) = G(nn,II1) + YTR(nn,I) * YTR(nn,I)
         enddo
         DO J = I+1, NS
            IJ1 = NG*(I-1) - (I*(I-1))/2 + J
            do nn = 1, np
               G(nn,IJ1) = G(nn,IJ1) + YTR(nn,I) * YTR(nn,J)
            enddo
         ENDDO
      ENDDO
C-----------------------------------------------------------------------
      CALL EGFDEC ( NP, NG, G, TEMP, IER )
      CALL DCOPY  ( np*NG, BETA, 1, AN, 1 )
      CALL EGFSOL ( NP, NG, G, AN, TEMP )
C-----------------------------------------------------------------------
      CALL EGFEVL ( NP, NG, PTC, AN, BETA, PATMOS, TEMPER )
      CALL EGFEVT ( NP, NS, THETA, AN, YTR, SUMTR, AAA )
C-----------------------------------------------------------------------
      CALL EGFEVD ( NP, NS, NG, AN, YTR, G, SUMTR, D, TEMP, AAA )
C-----------------------------------------------------------------------
C     Project the matrix D
C-----------------------------------------------------------------------
      CALL EGFPDP ( NP, NS, YTR, TEMP, SUMTR, D, AAA )
C-----------------------------------------------------------------------
      CALL EGFSCATDR ( NP, NS, PATMOS, WWTR, RU, TEMPER, THETA, D, AAA )
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGFTD1 (NP, PRES, IPRES, T, Y, WW, WEG, IWEG, 
     &                   THETA, D)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 7
C        ITLS  = 3
C
C     Input
C     -----
C        NP           number of nodes
C        PRES(*)      pressure
C        IPRES        flag for pressure (0:scalar, .ne.0:node dependent)
C        T(NP)        temperature
C        Y(NP,NS)     species mass fractions
C        WW(NP)       mean molecular weight
C        WEG          double precision work array for EGLIB
C        IWEG         integer work array for EGLIB
C
C     Output
C     ------
C        THETA(NP,NS) thermal diffusion vector
C        D(NP,NS,NS)  flux diffusion matrix
C
C        where rho is the density.
C     
C
C     Three CG iterations are performed on the matrix L in order
C     to evaluate THETA.
C     Two projected standard iterations are performed on the matrix
C     L_[00] in order to evaluate D.
C
C     Note
C     ----
C        THETA and D satisfy the mass conservation constraints.
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*), IWEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGFTD1 ( NP, NS, PRES, IPRES, T, WEG(IDLT1), WEG(IDLT2),
     &        WEG(IDLT3), WEG(IDLT4), WEG(IDLT5), WEG(IDLT6),
     &        WEG(IXTR), WEG(IYTR), WEG(IEGWT), WEG(IWWTR), WEG(ISUMTR), 
     &        WEG(IEGRU), WEG(IEGPA), THETA, D,
     &        WEG(IBIN), WEG(IAIJ), WEG(IBIJ), WEG(ICIJ), 
     &        WEG(IFITA), WEG(IFITB), WEG(IFITC),
     &        WEG(ICINT), WEG(IETA), WEG(ICXI), IWEG(IEGLIN), WEG(IG), 
     &        WEG(IDMI), WEG(IAN), WEG(IZN), WEG(IRN), 
     &        WEG(ITEMP), WEG(IBETA), WEG(IAUX),
     &        WEG(IAAA), WEG(IBBB) )
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGFTD1 ( NP, NS, PRES, IPRES, 
     &           TEMPER, DLT1, DLT2, DLT3, DLT4, DLT5, DLT6,  
     &           XTR, YTR, WT, WWTR, SUMTR, RU, PATMOS, THETA, D,
     &           BIN, AIJ, BIJ, CIJ, FITA, FITB, FITC,
     &           CINT, ETA, CXI, LIN,
     &           G, DMI, AN, ZN, RN, TEMP, BETA, AUX, AAA, BBB )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      CALL EGFEML ( NP, NS, TEMPER, DLT1, DLT2, DLT3,  
     &           DLT4, DLT5, DLT6, XTR, WT, RU, PATMOS,
     &           BIN, AIJ, BIJ, CIJ, FITA, FITB, FITC, CINT, ETA, CXI,
     &           LIN, G, BETA )
C-----------------------------------------------------------------------
C       Evaluate theta 
C-----------------------------------------------------------------------
      ITERMX = 3
      NG = 3 * NS
      CALL DCOPY  ( np*NG, BETA, 1, RN, 1 )
      CALL EGFCG3 ( NP, NS, NG, G, DMI, AN, ZN, RN, TEMP, 
     &              AAA, BBB, ITERMX )
C-----------------------------------------------------------------------
      CALL EGFEVT ( NP, NS, THETA, AN, YTR, SUMTR, AAA )
      CALL EGFEVDI1 ( NP, NS, D, WT, WWTR, TEMP, XTR, YTR, BIN, AUX )
C-----------------------------------------------------------------------
C     Project the matrix D
C-----------------------------------------------------------------------
      CALL EGFPDP ( NP, NS, YTR, TEMP, SUMTR, D, AAA )
C-----------------------------------------------------------------------
C        Correct for the pressure dependence of theta and D
C-----------------------------------------------------------------------
      IF ( IPRES .EQ. 0 ) THEN
         AA = PATMOS / PRES
         CALL DSCAL ( np*NS, AA, THETA, 1 )
         CALL DSCAL ( np*NS*NS, AA, D, 1 )
      ELSE
         CALL EGFSCATD ( NP, NS, PATMOS, PRES, THETA, D, AAA )
      ENDIF
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGFTDR1 (NP, T, Y, WW, WEG, IWEG, THETA, D)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 7
C        ITLS  = 3
C
C     Input
C     -----
C        NP           number of nodes
C        T(NP)        temperature
C        Y(NP,NS)     species mass fractions
C        WW(NP)       mean molecular weight
C        WEG          double precision work array for EGLIB
C        IWEG         integer work array for EGLIB
C
C     Output
C     ------
C        THETA(NP,NS) rho * thermal diffusion vector
C        D(NP,NS,NS)  rho * flux diffusion matrix
C
C        where rho is the density.
C     
C
C     Three CG iterations are performed on the matrix L in order
C     to evaluate THETA.
C     Two projected standard iterations are performed on the matrix
C     L_[00] in order to evaluate D.
C
C     Note
C     ----
C        THETA and D satisfy the mass conservation constraints.
C        These transport coefficients are also PRESSURE INDEPENDENT.
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*), IWEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGFTDR1 ( NP, NS, T, WEG(IDLT1), WEG(IDLT2),
     &       WEG(IDLT3), WEG(IDLT4), WEG(IDLT5), WEG(IDLT6),
     &       WEG(IXTR), WEG(IYTR), WEG(IEGWT), WEG(IWWTR), WEG(ISUMTR), 
     &       WEG(IEGRU), WEG(IEGPA), THETA, D,
     &       WEG(IBIN), WEG(IAIJ), WEG(IBIJ), WEG(ICIJ), 
     &       WEG(IFITA), WEG(IFITB), WEG(IFITC),
     &       WEG(ICINT), WEG(IETA), WEG(ICXI), IWEG(IEGLIN),
     &       WEG(IG), WEG(IDMI), WEG(IAN), WEG(IZN), 
     &       WEG(IRN), WEG(ITEMP), WEG(IBETA), WEG(IAUX),
     &       WEG(IAAA), WEG(IBBB) )
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGFTDR1 ( NP, NS, TEMPER, DLT1, DLT2, DLT3,  
     &           DLT4, DLT5, DLT6, XTR, YTR, WT, WWTR, SUMTR, 
     &           RU, PATMOS, THETA, D,
     &           BIN, AIJ, BIJ, CIJ, FITA, FITB, FITC,
     &           CINT, ETA, CXI, LIN,
     &           G, DMI, AN, ZN, RN, TEMP, BETA, AUX, AAA, BBB )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      CALL EGFEML ( NP, NS, TEMPER, DLT1, DLT2, DLT3,  
     &           DLT4, DLT5, DLT6, XTR, WT, RU, PATMOS,
     &           BIN, AIJ, BIJ, CIJ, FITA, FITB, FITC, CINT, ETA, CXI,
     &           LIN, G, BETA )
C-----------------------------------------------------------------------
C       Evaluate theta 
C-----------------------------------------------------------------------
      ITERMX = 3
      NG = 3 * NS
      CALL DCOPY  ( np*NG, BETA, 1, RN, 1 )
      CALL EGFCG3 ( NP, NS, NG, G, DMI, AN, ZN, RN, TEMP, 
     &              AAA, BBB, ITERMX )
C-----------------------------------------------------------------------
      CALL EGFEVT ( NP, NS, THETA, AN, YTR, SUMTR, AAA )
      CALL EGFEVDI1 ( NP, NS, D, WT, WWTR, TEMP, XTR, YTR, BIN, AUX )
C-----------------------------------------------------------------------
C     Project the matrix D
C-----------------------------------------------------------------------
      CALL EGFPDP ( NP, NS, YTR, TEMP, SUMTR, D, AAA )
C-----------------------------------------------------------------------
C        Correct for the pressure dependence of theta and D
C-----------------------------------------------------------------------
      CALL EGFSCATDR ( NP, NS, PATMOS, WWTR, RU, TEMPER, THETA, D, AAA )
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGFV1 (NP, PRES, IPRES, T, Y, WW, WEG, D)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 2
C        ITLS  = 0
C
C     Input
C     -----
C        NP           number of nodes
C        PRES(*)      pressure
C        IPRES        flag for pressure (0:scalar, .ne.0:node dependent)
C        T(NP)        temperature
C        Y(NP,NS)     species mass fractions
C        WW(NP)       mean molecular weight
C        WEG          double precision work array for EGLIB
C
C     Output
C     ------
C        D(NP,NS)     flux diffusion coefficients
C
C     
C     The array D corresponds to the diagonal of the matrix
C     ~D_[00]^1 before projection.
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGFV1 ( NP, NS, PRES, IPRES, WEG(IXTR), WEG(IYTR), 
     &              WEG(IEGWT), WEG(IWWTR), WEG(IEGPA), D, WEG(IBIN), 
     &              WEG(IAUX), WEG(IAAA) )
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGFV1 ( NP, NS, PRES, IPRES, XTR, YTR, WT, WWTR, 
     &                    PATMOS, D, BIN, AUX, AAA )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XTR(NP,NS), YTR(NP,NS), WT(NS), D(NP,NS), 
     &          BIN(NP,*), AUX(NP,NS), WWTR(NP)
C-----------------------------------------------------------------------
      CALL EGZERO ( np*NS, D )
      DO I = 1, NS
         IP1 = I + 1
         DO J = IP1, NS
            IJ = NS*(I-1) - (I*(I-1))/2 + J
            do nn = 1, np
               AAI = XTR(nn,I) * BIN(nn,IJ)
               AAJ = XTR(nn,J) * BIN(nn,IJ)
               D(nn,I) = D(nn,I) + AAJ
               D(nn,J) = D(nn,J) + AAI
            enddo
         ENDDO
      ENDDO
C-----------------------------------------------------------------------
      DO I = 1, NS
         do nn = 1, np
            D(nn,I) = WT(I) * AUX(nn,I) / ( D(nn,I) * WWTR(nn) )
         enddo
      ENDDO
C-----------------------------------------------------------------------
C        Correct for the pressure dependence of D
C-----------------------------------------------------------------------
      IF ( IPRES .EQ. 0 ) THEN
         AA = PATMOS / PRES
         CALL DSCAL ( np*NS, AA, D, 1 )
      ELSE
         CALL EGFSCAV ( NP, NS, PATMOS, PRES, D, AAA )
      ENDIF
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGFVR1 (NP, T, Y, WEG, D)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 2
C        ITLS  = 0
C
C     Input
C     -----
C        NP           number of nodes
C        T(NP)        temperature
C        Y(NP,NS)     species mass fractions
C        WEG          double precision work array for EGLIB
C
C     Output
C     ------
C        D(NP,NS)     rho * flux diffusion coefficients
C
C        where rho is the density.
C
C     
C     The array D corresponds to the diagonal of the matrix
C     rho * ~D_[00]^1 before projection.
C     
C
C     Note
C     ----
C        The coefficients D are PRESSURE INDEPENDENT.
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGFVR1 ( NP, NS, T, WEG(IXTR), WEG(IEGWT), 
     &               WEG(IEGRU), WEG(IEGPA), D, WEG(IBIN), 
     &               WEG(IAUX) )
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGFVR1 ( NP, NS, TEMP, XTR, WT, RU, PATMOS, 
     &                     D, BIN, AUX )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XTR(NP,NS), WT(NS), D(NP,NS), BIN(NP,*), AUX(NP,NS)
      DIMENSION TEMP(NP)
C-----------------------------------------------------------------------
      CALL EGZERO ( np*NS, D )
      DO I = 1, NS
         IP1 = I + 1
         DO J = IP1, NS
            IJ = NS*(I-1) - (I*(I-1))/2 + J
            do nn = 1, np
               AAI = XTR(nn,I) * BIN(nn,IJ)
               AAJ = XTR(nn,J) * BIN(nn,IJ)
               D(nn,I) = D(nn,I) + AAJ
               D(nn,J) = D(nn,J) + AAI
            enddo
         ENDDO
      ENDDO
C-----------------------------------------------------------------------
      DO I = 1, NS
         FAC = WT(I) * PATMOS / RU
         do nn = 1, np
            D(nn,I) = FAC * AUX(nn,I) / ( D(nn,I) * TEMP(nn) )
         enddo
      ENDDO
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGFYV (NP, PRES, IPRES, Y, WW, WEG, F, IDDEC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 2
C        ITLS  = 1
C
C     Input
C     -----
C        NP           number of nodes
C        PRES         pressure
C        IPRES        flag for pressure (0:scalar, .ne.0:node dependent)
C        Y(NP,NS)     species mass fractions
C        WW           mean molecular weight
C        WEG          double precision work array for EGLIB
C        F(NP,NS)     species diffusion driving forces
C        IDDEC        flag indicating whether the diffusion matrix
C                     needs first to be assembled and decomposed
C                     0: no, .ne.0: yes
C
C     Output
C     ------
C        F(NP,NS)     species flux diffusion velocities
C                     F_i <-- - sum_j \widetilde D_{ij} F_j  
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGFYV ( NP, NS, PRES, IPRES, WEG(IXTR), WEG(IYTR), 
     &              WEG(IEGPA), WEG(IG), WEG(ITEMP), 
     &              WEG(IBIN), WEG(IAAA), F, IDDEC )
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGFYV ( NP, NS, PRES, IPRES, XTR, YTR, 
     &                    PATMOS, G, TEMP, BIN, AAA, F, IDDEC )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION F(NP,NS), YTR(NP,NS), PRES(*), AAA(NP)
C-----------------------------------------------------------------------
      IF ( IDDEC .NE. 0 ) THEN
         CALL EGFDDEC  ( NP, NS, XTR, YTR, G, BIN, TEMP )
      ENDIF
      CALL EGFSOL  ( NP, NS, G, F, AAA )
      IF ( IPRES .EQ. 0 ) THEN
         PPP = - PATMOS / PRES(1)
         DO I = 1, NS
            do nn = 1, np
               F(nn,I) = F(nn,I) * YTR(nn,I) * PPP
            enddo
         ENDDO
      ELSE
         do nn = 1, np
            AAA(nn) = - PATMOS / PRES(nn) 
         enddo
         DO I = 1, NS
            do nn = 1, np
               F(nn,I) = F(nn,I) * YTR(nn,I) * AAA(nn)
            enddo
         ENDDO
      ENDIF
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGFRYV (NP, T, Y, WEG, F, IDDEC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Minimum flags required in EGINI
C     -------------------------------
C        IFLAG = 2
C        ITLS  = 1
C
C     Input
C     -----
C        NP           number of nodes
C        T            temperature
C        Y(NP,NS)     species mass fractions
C        WEG          double precision work array for EGLIB
C        F(NP,NS)     species diffusion driving forces
C        IDDEC        flag indicating whether the diffusion matrix
C                     needs first to be assembled and decomposed
C                     0: no, .ne.0: yes
C
C     Output
C     ------
C        F(NP,NS)     rescaled species flux diffusion velocities
C                     F_i <-- - rho * sum_j \widetilde D_{ij} F_j  
C
C-----------------------------------------------------------------------
      DIMENSION WEG(*)
      INCLUDE 'eg.cmn'
C-----------------------------------------------------------------------
      CALL LEGFRYV ( NP, NS, WEG(IXTR), WEG(IYTR), WEG(IEGWT),
     &               WEG(IEGPA), WEG(IEGRU), T, WEG(IG), WEG(ITEMP), 
     &               WEG(IBIN), WEG(IAAA), F, IDDEC )
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE LEGFRYV ( NP, NS, XTR, YTR, WT, PATMOS, RU, TEMPER, 
     &                     G, TEMP, BIN, AAA, F, IDDEC )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION TEMPER(NP), F(NP,NS), WT(NS), AAA(NP), XTR(NP,NS)
C-----------------------------------------------------------------------
      IF ( IDDEC .NE. 0 ) THEN
         CALL EGFDDEC  ( NP, NS, XTR, YTR, G, BIN, TEMP )
      ENDIF
      CALL EGFSOL   ( NP, NS, G, F, AAA )
C-----------------------------------------------------------------------
      do nn = 1, np
         AAA(nn) = - PATMOS / RU / TEMPER(nn)
      enddo
      DO I = 1, NS
         do nn = 1, np
            F(nn,I) = F(nn,I) * WT(I) * XTR(nn,I) * AAA(nn)
         enddo
      ENDDO
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGFDOT ( NP, NG, DOT, AN, BETA )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      DIMENSION DOT(NP), AN(NP,NG), BETA(NP,NG)
C-----------------------------------------------------------------------
      call egzero ( np, dot )
      DO I = 1, ng
         do nn = 1, np
            dot(nn) = dot(nn) + AN(nn,I) * BETA(nn,I)
         enddo
      ENDDO
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGFEVCT ( NP, NS, CHIT, AN, XTR, BIN, CIJ, FITC,
     &                     DLT1, DLT2, DLT3, DLT4, DLT5, DLT6, WT )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      DIMENSION CHIT(NP,NS), AN(NP,NS), XTR(NP,NS), BIN(NP,*),
     &          CIJ(NP), FITC(7,NS,NS), WT(NS) 
      DIMENSION DLT1(NP), DLT2(NP), DLT3(NP), DLT4(NP), DLT5(NP),
     &          DLT6(NP)
C-----------------------------------------------------------------------
      CALL EGZERO ( np*NS, CHIT )
      DO I = 1, NS
         DO J = I+1, NS
            IJ = NS*(I-1) - (I*(I-1))/2 + J
C........
            FITC0 = FITC(1,I,J)
            FITC1 = FITC(2,I,J)
            FITC2 = FITC(3,I,J)
            FITC3 = FITC(4,I,J)
            FITC4 = FITC(5,I,J)
            FITC5 = FITC(6,I,J)
            FITC6 = FITC(7,I,J)
            do nn = 1, np
               CIJ(nn) = FITC0          + FITC1*DLT1(nn) 
     &                 + FITC2*DLT2(nn) + FITC3*DLT3(nn)
     &                 + FITC4*DLT4(nn) + FITC5*DLT5(nn)
     &                 + FITC6*DLT6(nn) 
            enddo
C........
            do nn = 1, np
               TERM = 0.5D0 * BIN(nn,IJ) * (6.0D0*CIJ(nn) - 5.0D0 )
     &             * ( WT(I)*AN(nn,J) - WT(J)*AN(nn,I) ) 
     &             / ( WT(I)+WT(J) )
               CHIT(nn,I) = CHIT(nn,I) + XTR(nn,J) * TERM
               CHIT(nn,J) = CHIT(nn,J) - XTR(nn,I) * TERM
            enddo
         ENDDO
      ENDDO
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGFEVD ( NP, NS, NG, AN, YTR, G, SUMTR, D, TEMP, AAA )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      DIMENSION AN(NP,NG), YTR(NP,NS), G(NP,*), TEMP(NP,NG), AAA(NP),
     &          SUMTR(NP), D(NP,NS,NS)
C-----------------------------------------------------------------------
      do nn = 1, np
         aaa(nn) = 1.d0 / sumtr(nn)
      enddo
      DO I = 1, NS
         CALL EGZERO ( np*NG, AN )
         CALL DCOPY ( np*NS, YTR, 1, AN, 1 )
         DO J = 1, NS
            do nn = 1, np
               AN(nn,J) = - aaa(nn) * AN(nn,J)
            enddo
         ENDDO
         do nn = 1, np
            AN(nn,I) = AN(nn,I) + 1.0D0
         enddo
         CALL EGFSOL ( NP, NG, G, AN, TEMP )
         DO J = 1, NS
            do nn = 1, np
               D(nn,I,J) = YTR(nn,I) * AN(nn,J)
            enddo
         ENDDO
      ENDDO
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGFEVDI1 ( NP, NS, D, WT, WWTR, TEMP, 
     &                      XTR, YTR, BIN, AUX )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      DIMENSION WT(NS), WWTR(NP), TEMP(NP,NS), XTR(NP,NS), YTR(NP,NS),
     &          BIN(NP,*), AUX(NP,NS)
      DIMENSION D(NP,NS,NS)
C-----------------------------------------------------------------------
      CALL EGZERO ( np*NS, TEMP )
      DO I = 1, NS
         DO J = I+1, NS
            IJ = NS*(I-1) - (I*(I-1))/2 + J
            do nn = 1, np
               TEMP(nn,I) = TEMP(nn,I) + XTR(nn,J) * BIN(nn,IJ)
               TEMP(nn,J) = TEMP(nn,J) + XTR(nn,I) * BIN(nn,IJ)
            enddo
         ENDDO
      ENDDO
      DO I = 1, NS
         do nn = 1, np
            TEMP(nn,I) = AUX(nn,I) / TEMP(nn,I)
         enddo
      ENDDO
      DO I = 1, NS
         DO J = I+1, NS
            IJ = NS*(I-1) - (I*(I-1))/2 + J
            do nn = 1, np
               D(nn,I,J) = WT(I) / WWTR(nn) * TEMP(nn,I) * TEMP(nn,J) 
     &                                      * XTR(nn,I) * BIN(nn,IJ)
               D(nn,J,I) = WT(J) / WWTR(nn) * TEMP(nn,I) * TEMP(nn,J) 
     &                                      * XTR(nn,J) * BIN(nn,IJ)
            enddo
         ENDDO
         do nn = 1, np
            D(nn,I,I) = WT(I) / WWTR(nn) * TEMP(nn,I) 
     &                                   * ( 1.0D0 + YTR(nn,I) )
         enddo
      ENDDO
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGFEVK ( NP, NS, VV, G, AN, BETA, TEMP, VV01, VVS )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      DIMENSION G(NP,*), BETA(NP,*), AN(NP,*), TEMP(NP,*)
      DIMENSION VV(NP), VV01(NP), VVS(NP)
C-----------------------------------------------------------------------
      NS2  = 2 * NS
      DO I = NS+1, NS2
         III  = NS2*(I-1) - (I*(I-1))/2 + I
         do nn = 1, np
            AN(nn,I-NS) = BETA(nn,I) / G(nn,III)
         enddo
      ENDDO
      DO I = 1, NS
         call dcopy ( np, BETA(1,i), 1, TEMP(1,i), 1 )
         DO J = 1, NS
            IJ = NS2*(I-1) - (I*(I-1))/2 + J + NS
            do nn = 1, np
               TEMP(nn,I) = TEMP(nn,I) - G(nn,IJ) * AN(nn,J)
            enddo
         ENDDO
      ENDDO
      call egzero ( np, vv01 )
      call egzero ( np, vvs )
      DO I = 1, NS
         II = NS2*(I-1) - (I*(I-1))/2 + I 
         do nn = 1, np
            VV01(nn) = VV01(nn) + AN(nn,I) * BETA(nn,I+NS)
            VVS(nn)  = VVS(nn)  + TEMP(nn,I) * TEMP(nn,I) / G(nn,II)
         enddo
      ENDDO
      do nn = 1, np
         VV(nn) = VV01(nn) + VVS(nn)
      enddo
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGFEVL ( NP, NG, TC, AN, BETA, PATMOS, TEMPER )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      DIMENSION TC(NP), AN(NP,NG), BETA(NP,NG), TEMPER(NP)
C-----------------------------------------------------------------------
      call egzero ( np, tc )
      DO I = 1, NG
         do nn = 1, np
            TC(nn) = TC(nn) + AN(nn,I) * BETA(nn,I)
         enddo
      ENDDO
      do nn = 1, np
         tc(nn) = tc(nn) * patmos / temper(nn)
      enddo
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGFEVT ( NP, NS, THETA, AN, YTR, SUMTR, AAA )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      DIMENSION THETA(NP,NS), AN(NP,NS), YTR(NP,NS), AAA(NP), SUMTR(NP)
C-----------------------------------------------------------------------
      call egzero ( np, aaa )
      DO I = 1, NS
         do nn = 1, np
            aaa(nn) = aaa(nn) + YTR(nn,I) * AN(nn,I)
         enddo
      ENDDO
      DO I = 1, NS
         do nn = 1, np
            THETA(nn,I) =  aaa(nn) / sumtr(nn) - AN(nn,I)
         enddo
      ENDDO
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGFSCAD ( NP, NS, PATMOS, PRES, D, AAA )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      DIMENSION D(NP,NS,NS), PRES(NP), AAA(NP)
C-----------------------------------------------------------------------
      do nn = 1, np
         aaa(nn) = patmos / pres(nn)
      enddo
      do i = 1, ns
         do j = 1, ns
            do nn = 1, np
               d(nn,i,j) = aaa(nn) * d(nn,i,j)
            enddo
         enddo
      enddo
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGFSCAV ( NP, NS, PATMOS, PRES, D, AAA )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      DIMENSION D(NP,NS), PRES(NP), AAA(NP)
C-----------------------------------------------------------------------
      do nn = 1, np
         aaa(nn) = patmos / pres(nn)
      enddo
      do i = 1, ns
         do nn = 1, np
            d(nn,i) = aaa(nn) * d(nn,i)
         enddo
      enddo
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGFSCATD ( NP, NS, PATMOS, PRES, THETA, D, AAA )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      DIMENSION THETA(NP,NS), D(NP,NS,NS), PRES(NP), AAA(NP)
C-----------------------------------------------------------------------
      do nn = 1, np
         aaa(nn) = patmos / pres(nn)
      enddo
      do i = 1, ns
         do nn = 1, np
            theta(nn,i) = aaa(nn) * theta(nn,i)
         enddo
      enddo
      do i = 1, ns
         do j = 1, ns
            do nn = 1, np
               d(nn,i,j) = aaa(nn) * d(nn,i,j)
            enddo
         enddo
      enddo
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGFSCADR ( NP, NS, PATMOS, WWTR, RU, TEMPER, D, AAA )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      DIMENSION WWTR(NP), TEMPER(NP), D(NP,NS,NS), AAA(NP)
C-----------------------------------------------------------------------
      do nn = 1, np
         aaa(nn) = PATMOS * WWTR(nn) / ( RU * TEMPER(nn) )
      enddo
      do i = 1, ns
         do j = 1, ns
            do nn = 1, np
               d(nn,i,j) = aaa(nn) * d(nn,i,j)
            enddo
         enddo
      enddo
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGFSCATDR ( NP, NS, PATMOS, WWTR, RU, TEMPER, 
     &                       THETA, D, AAA )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      DIMENSION WWTR(NP), TEMPER(NP), THETA(NP,NS), D(NP,NS,NS), AAA(NP)
C-----------------------------------------------------------------------
      do nn = 1, np
         aaa(nn) = PATMOS * WWTR(nn) / ( RU * TEMPER(nn) )
      enddo
      do i = 1, ns
         do nn = 1, np
            theta(nn,i) = aaa(nn) * theta(nn,i)
         enddo
      enddo
      do i = 1, ns
         do j = 1, ns
            do nn = 1, np
               d(nn,i,j) = aaa(nn) * d(nn,i,j)
            enddo
         enddo
      enddo
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGFAXS ( NP, N, A, X, B )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      DIMENSION A(NP,*), X(NP,N), B(NP,N)
C-----------------------------------------------------------------------
C                             B == A.X
C
C     The matrix A is stored in symmetric form, i.e., instead of
C     A(nn,i,j) we store A(nn,n*(i-1) - i*(i-1)/2 + j)
C-----------------------------------------------------------------------
      DO 10 I = 1, N
         III = N*(I-1) - (I*(I-1))/2 + I
         do nn = 1, np
            B(nn,I) = A(nn,III) * X(nn,I)
         enddo
10    CONTINUE
      DO 30 I = 1, N
         DO 20 J = I+1, N
            III = N*(I-1) - (I*(I-1))/2 + J
            do nn = 1, np
               B(nn,I) = B(nn,I) + A(nn,III)*X(nn,J)
               B(nn,J) = B(nn,J) + A(nn,III)*X(nn,I)
            enddo
20       CONTINUE
30    CONTINUE
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGFCG1 ( NP, NG, G, DMI, AN, ZN, RN, TEMP, 
     &                    AAA, BBB, ITERMX )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      DIMENSION G(NP,*), DMI(NP,NG),  AN(NP,NG), ZN(NP,NG),  RN(NP,NG),
     &          TEMP(NP,NG)
      DIMENSION AAA(NP), BBB(NP)
C-----------------------------------------------------------------------
C           INITIALISATION DES ITERATIONS
C-----------------------------------------------------------------------
      NITER = 0
      DO I = 1, NG
         III = NG*(I-1) - (I*(I-1))/2 + I
         do nn = 1, np
            AN(nn,I) = 0.0D0
            ZN(nn,I) = 0.0D0
            DMI(nn,I) = 1.0D0 / G(nn,III)
         enddo
      ENDDO
      call egzero ( np, aaa )
      call egzero ( np, bbb )
      DO I = 1, NG
         do nn = 1, np
            AAA(nn) = AAA(nn) + DMI(nn,I) * RN(nn,I)*RN(nn,I)
         enddo
      ENDDO
C-----------------------------------------------------------------------
 100  CONTINUE
      NITER = NITER + 1
      DO I = 1, NG
         do nn = 1, np
            ZN(nn,I) = DMI(nn,I)*RN(nn,I) + BBB(nn)*ZN(nn,I)
         enddo
      ENDDO
      CALL EGFAXS(np, NG, G, ZN, TEMP)
      call egzero ( np, bbb )
      do i = 1, ng
         do nn = 1, np
            bbb(nn) = bbb(nn) + ZN(nn,I) * TEMP(nn,I)
         enddo
      enddo
      DO I = 1, NG
         do nn = 1, np
            AN(nn,I) = AN(nn,I) + AAA(nn)/BBB(nn)*ZN(nn,I)
            RN(nn,I) = RN(nn,I) - AAA(nn)/BBB(nn)*TEMP(nn,I)
         enddo
      ENDDO
      call egzero ( np, bbb )
      DO I = 1, NG
         do nn = 1, np
            BBB(nn) = BBB(nn) + DMI(nn,I) * RN(nn,I) * RN(nn,I)
         enddo
      ENDDO
      do nn = 1, np
         save    = aaa(nn)
         aaa(nn) = bbb(nn)
         bbb(nn) = bbb(nn) / save
      enddo
      IF ( NITER .LT. ITERMX ) GO TO 100
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGFCG2 ( NP, NS, NG, G, DMI, AN, ZN, RN, TEMP, 
     &                    AAA, BBB, ITERMX )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      DIMENSION G(NP,*), DMI(NP,3,NS),  AN(NP,NG), ZN(NP,NG),  
     &          RN(NP,NG), TEMP(NP,NG)
      DIMENSION AAA(NP), BBB(NP)
C-----------------------------------------------------------------------
C           INITIALISATION DES ITERATIONS
C
C     The preconditioner matrix DMI is stored in compact form
C     DMI(nn,n,i) with n referring to
C                   1  2
C            DMI =  2  3
C-----------------------------------------------------------------------
      NITER = 0
      DO I = 1, NG
         do nn = 1, np
            AN(nn,I) = 0.0D0
            ZN(nn,I) = 0.0D0
         enddo
      ENDDO
      DO I = 1, NS
         II1 = NG*(I-1) - (I*(I-1))/2 + I
         II2 = NG*(I-1) - (I*(I-1))/2 + I + NS
         II3 = NG*(I+NS-1) - ((I+NS)*(I+NS-1))/2 + I + NS
         do nn = 1, np
            DETM  = 1.0D0 / 
     &              ( (G(nn,II1)*G(nn,II3) - G(nn,II2)*G(nn,II2)) )
C........
            DMI(nn,1,I) =   G(nn,II3) * DETM
            DMI(nn,2,I) = - G(nn,II2) * DETM
            DMI(nn,3,I) =   G(nn,II1) * DETM
C........
         enddo
      ENDDO
      call egzero ( np, aaa )
      call egzero ( np, bbb )
      DO I = 1, NS
         do nn = 1, np
            AAA(nn) = AAA(nn) 
     &              + RN(nn,I) *    ( DMI(nn,1,I) * RN(nn,I) 
     &                              + DMI(nn,2,I) * RN(nn,I+NS) )
     &              + RN(nn,I+NS) * ( DMI(nn,2,I) * RN(nn,I) 
     &                              + DMI(nn,3,I) * RN(nn,I+NS) )
         enddo
      ENDDO
C-----------------------------------------------------------------------
 100  CONTINUE
      NITER = NITER + 1
      DO I = 1, NS
         do nn = 1, np
            ZN(nn,I)    = DMI(nn,1,I)*RN(nn,I) 
     &                  + DMI(nn,2,I)*RN(nn,I+NS) + BBB(nn)*ZN(nn,I)
            ZN(nn,I+NS) = DMI(nn,2,I)*RN(nn,I) 
     &                  + DMI(nn,3,I)*RN(nn,I+NS) + BBB(nn)*ZN(nn,I+NS)
         enddo
      ENDDO
      CALL EGFAXS(np, NG, G, ZN, TEMP)
      call egzero ( np, bbb )
      do i = 1, ng
         do nn = 1, np
            bbb(nn) = bbb(nn) + ZN(nn,I) * TEMP(nn,I)
         enddo
      enddo
      DO I = 1, NG
         do nn = 1, np
            AN(nn,I) = AN(nn,I) + AAA(nn)/BBB(nn)*ZN(nn,I)
            RN(nn,I) = RN(nn,I) - AAA(nn)/BBB(nn)*TEMP(nn,I)
         enddo
      ENDDO
      call egzero ( np, bbb )
      DO I = 1, NS
         do nn = 1, np
            BBB(nn) = BBB(nn) 
     &              + RN(nn,I) *    ( DMI(nn,1,I) * RN(nn,I) 
     &                              + DMI(nn,2,I) * RN(nn,I+NS) )
     &              + RN(nn,I+NS) * ( DMI(nn,2,I) * RN(nn,I) 
     &                              + DMI(nn,3,I) * RN(nn,I+NS) )
         enddo
      ENDDO
      do nn = 1, np
         save    = aaa(nn)
         aaa(nn) = bbb(nn)
         bbb(nn) = bbb(nn) / save
      enddo
      IF ( NITER .LT. ITERMX ) GO TO 100
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGFCG3 ( NP, NS, NG, G, DMI, AN, ZN, RN, TEMP, 
     &                    AAA, BBB, ITERMX )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      DIMENSION G(NP,*), DMI(NP,6,NS),  AN(NP,NG), ZN(NP,NG),  
     &          RN(NP,NG), TEMP(NP,NG)
      DIMENSION AAA(NP), BBB(NP)
C-----------------------------------------------------------------------
C           INITIALISATION DES ITERATIONS
C
C     The preconditioner matrix DMI is stored in compact form
C     DMI(nn,n,i) with n referring to
C                   1  2  3
C            DMI =  2  4  5
C                   3  5  6
C-----------------------------------------------------------------------
      NS2 = 2 * NS
      NS3 = 3 * NS
      NITER = 0
      DO I = 1, NG
         do nn = 1, np
            AN(nn,I) = 0.0D0
            ZN(nn,I) = 0.0D0
         enddo
      ENDDO
      DO I = 1, NS
         II1 = NS3*(I-1) - (I*(I-1))/2 + I
         II2 = NS3*(I-1) - (I*(I-1))/2 + I + NS
         II3 = NS3*(I-1) - (I*(I-1))/2 + I + NS2
         II4 = NS3*(I+NS-1) - ((I+NS)*(I+NS-1))/2 + I + NS
         II5 = NS3*(I+NS-1) - ((I+NS)*(I+NS-1))/2 + I + NS2
         II6 = NS3*(I+NS2-1) - ((I+NS2)*(I+NS2-1))/2 + I + NS2
         do nn = 1, np
            DETM = 1.0D0/( G(nn,II1) * G(nn,II4) * G(nn,II6)
     &                  -  G(nn,II1) * G(nn,II5) * G(nn,II5)
     &                  -  G(nn,II2) * G(nn,II2) * G(nn,II6) )
C...........
            DMI(nn,1,I) = ( G(nn,II4) * G(nn,II6) 
     &                    - G(nn,II5) * G(nn,II5) ) * DETM
            DMI(nn,2,I) = - G(nn,II2) * G(nn,II6) * DETM
            DMI(nn,3,I) =   G(nn,II2) * G(nn,II5) * DETM
            DMI(nn,4,I) =   G(nn,II1) * G(nn,II6) * DETM
            DMI(nn,5,I) = - G(nn,II5) * G(nn,II1) * DETM
            DMI(nn,6,I) = ( G(nn,II1) * G(nn,II4) 
     &                    - G(nn,II2) * G(nn,II2) ) * DETM
         enddo
      ENDDO
C-----------------------------------------------------------------------
      call egzero ( np, aaa )
      call egzero ( np, bbb )
      DO I = 1, NS
         do nn = 1, np
            FAC = RN(nn,I) *     ( DMI(nn,1,I) * RN(nn,I) 
     &                           + DMI(nn,2,I) * RN(nn,I+NS) 
     &                           + DMI(nn,3,I) * RN(nn,I+NS2) )
     &          + RN(nn,I+NS) *  ( DMI(nn,2,I) * RN(nn,I) 
     &                           + DMI(nn,4,I) * RN(nn,I+NS) 
     &                           + DMI(nn,5,I) * RN(nn,I+NS2) )
     &          + RN(nn,I+NS2) * ( DMI(nn,3,I) * RN(nn,I) 
     &                           + DMI(nn,5,I) * RN(nn,I+NS) 
     &                           + DMI(nn,6,I) * RN(nn,I+NS2) ) 
            AAA(nn) = AAA(nn) + FAC
         enddo
      ENDDO
C-----------------------------------------------------------------------
 100  CONTINUE
      NITER = NITER + 1
      DO I = 1, NS
         do nn = 1, np
            ZN(nn,I)     = DMI(nn,1,I) * RN(nn,I) 
     &                   + DMI(nn,2,I) * RN(nn,I+NS) 
     &                   + DMI(nn,3,I) * RN(nn,I+NS2) 
     &                   + BBB(nn) * ZN(nn,I)
            ZN(nn,I+NS)  = DMI(nn,2,I) * RN(nn,I) 
     &                   + DMI(nn,4,I) * RN(nn,I+NS) 
     &                   + DMI(nn,5,I) * RN(nn,I+NS2) 
     &                   + BBB(nn) * ZN(nn,I+NS)
            ZN(nn,I+NS2) = DMI(nn,3,I) * RN(nn,I) 
     &                   + DMI(nn,5,I) * RN(nn,I+NS) 
     &                   + DMI(nn,6,I) * RN(nn,I+NS2) 
     &                   + BBB(nn) * ZN(nn,I+NS2)
         enddo
      ENDDO
      CALL EGFAXS(np, NG, G, ZN, TEMP)
      call egzero ( np, bbb )
      do i = 1, ng
         do nn = 1, np
            bbb(nn) = bbb(nn) + ZN(nn,I) * TEMP(nn,I)
         enddo
      enddo
      DO I = 1, NG
         do nn = 1, np
            AN(nn,I) = AN(nn,I) + AAA(nn)/BBB(nn)*ZN(nn,I)
            RN(nn,I) = RN(nn,I) - AAA(nn)/BBB(nn)*TEMP(nn,I)
         enddo
      ENDDO
      call egzero ( np, bbb )
      DO I = 1, NS
         do nn = 1, np
            FAC = RN(nn,I) *     ( DMI(nn,1,I) * RN(nn,I) 
     &                           + DMI(nn,2,I) * RN(nn,I+NS) 
     &                           + DMI(nn,3,I) * RN(nn,I+NS2) )
     &          + RN(nn,I+NS) *  ( DMI(nn,2,I) * RN(nn,I) 
     &                           + DMI(nn,4,I) * RN(nn,I+NS) 
     &                           + DMI(nn,5,I) * RN(nn,I+NS2) )
     &          + RN(nn,I+NS2) * ( DMI(nn,3,I) * RN(nn,I) 
     &                           + DMI(nn,5,I) * RN(nn,I+NS) 
     &                           + DMI(nn,6,I) * RN(nn,I+NS2) ) 
            BBB(nn) = BBB(nn) + FAC
         enddo
      ENDDO
      do nn = 1, np
         save    = aaa(nn)
         aaa(nn) = bbb(nn)
         bbb(nn) = bbb(nn) / save
      enddo
      IF ( NITER .LT. ITERMX ) GO TO 100
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGFSI1 ( NP, NS, G, DMI, TEMP, YTR, SUMTR,
     &                    D, ITERMX, AUX, AAA )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION G(NP,*), DMI(NP,NS), TEMP(NP,NS), YTR(NP,NS),  
     &          D(NP,NS,NS), AAA(NP), SUMTR(NP), AUX(NP,NS)
C-----------------------------------------------------------------------
C     At most two iterations are performed.
C-----------------------------------------------------------------------
      CALL EGZERO ( NP * NS * NS, D )
      DO I = 1, NS
         III = NS*(I-1) - (I*(I-1))/2 + I
         do nn = 1, np
            DMI(nn,I) = AUX(nn,I) / G(nn,III)
            D(nn,I,I) = DMI(nn,I)
         enddo
      ENDDO
C-----------------------------------------------------------------------
      IF (ITERMX .GT. 1) THEN
         DO I = 1, NS
            do nn = 1, np
               D(nn,I,I) = DMI(nn,I) * ( 1.0D0 + YTR(nn,I) )
            enddo
            DO J = I+1, NS
               IND = NS*(I-1) - (I*(I-1))/2 + J
               do nn = 1, np
                  D(nn,I,J) = - DMI(nn,I) * DMI(nn,J) * G(nn,IND)
                  D(nn,J,I) = D(nn,I,J)
               enddo
            ENDDO
         ENDDO
      ENDIF
C-----------------------------------------------------------------------
C     Matrix tilde D
C-----------------------------------------------------------------------
      DO I = 1, NS
         DO J = 1, NS
            do nn = 1, np
               D(nn,I,J) = YTR(nn,i) * D(nn,I,J)
            enddo
         ENDDO
      ENDDO
C-----------------------------------------------------------------------
C     Project tilde D
C-----------------------------------------------------------------------
      CALL EGFPDP ( np, NS, YTR, TEMP, SUMTR, D, AAA )
C-----------------------------------------------------------------------
      RETURN                                                                    
      END                                                                       
C=======================================================================
C=======================================================================
      SUBROUTINE EGFSI2 ( NP, NS, NG, G, DMI, TEMP, YTR, SUMTR,
     &                    D, ITERMX, AUX, AAA )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION G(NP,*), DMI(NP,3,NS), TEMP(NP,NG), YTR(NP,NS), 
     &          SUMTR(NP), D(NP,NS,NS), AAA(NP), AUX(NP,NS)
C-----------------------------------------------------------------------
C     At most two iterations are performed
C-----------------------------------------------------------------------
      NS2 = 2 * NS
      CALL EGZERO ( NP * NS * NS, D )
      DO I = 1, NS
         II1 = NS2*(I-1) - (I*(I-1))/2 + I
         II2 = NS2*(I-1) - (I*(I-1))/2 + I + NS
         II4 = NS2*(I+NS-1) - ((I+NS)*(I+NS-1))/2 + I + NS
         do nn = 1, np
            FAC  = 1.0D0 / AUX(nn,I)
            DETM = 1.0D0 / 
     &             ( FAC*G(nn,II1)*G(nn,II4) - G(nn,II2)*G(nn,II2) )
            DMI(nn,1,I) =   G(nn,II4) * DETM
            DMI(nn,2,I) = - G(nn,II2) * DETM
            D(nn,I,I) = DMI(nn,1,I)
         enddo
      ENDDO
C-----------------------------------------------------------------------
      IF (ITERMX .GT. 1) THEN
         DO I = 1, NS
            II1 = NS2*(I-1) - (I*(I-1))/2 + I
            II2 = NS2*(I-1) - (I*(I-1))/2 + I + NS
            II3 = NS2*(I+NS-1) - ((I+NS)*(I+NS-1))/2 + I + NS
            do nn = 1, np
               D(nn,I,I) = 2.0D0 * DMI(nn,1,I) 
     &                        - DMI(nn,1,I)*DMI(nn,1,I)*G(nn,II1)
     &                        - DMI(nn,1,I)*DMI(nn,2,I)*G(nn,II2)*2.0D0
     &                        - DMI(nn,2,I)*DMI(nn,2,I)*G(nn,II3)
            enddo
            DO J = I+1, NS
               IJ1 = NS2*(I-1) - (I*(I-1))/2 + J
               IJ2 = NS2*(I-1) - (I*(I-1))/2 + J + NS
               JI2 = NS2*(J-1) - (J*(J-1))/2 + I + NS
               IJ3 = NS2*(I+NS-1) - ((I+NS)*(I+NS-1))/2 + J + NS
C..............
               do nn = 1, np
                  D(nn,I,J) = 
     &                  - ( DMI(nn,1,I) * DMI(nn,1,J) * G(nn,IJ1)
     &                    + DMI(nn,1,I) * DMI(nn,2,J) * G(nn,IJ2)
     &                    + DMI(nn,2,I) * DMI(nn,1,J) * G(nn,JI2)
     &                    + DMI(nn,2,I) * DMI(nn,2,J) * G(nn,IJ3) )
                  D(nn,J,I) = D(nn,I,J)
               enddo
            ENDDO
         ENDDO
      ENDIF
C-----------------------------------------------------------------------
C     Matrix tilde D
C-----------------------------------------------------------------------
      DO I = 1, NS
         DO J = 1, NS
            do nn = 1, np
               D(nn,I,J) = YTR(nn,i) * D(nn,I,J)
            enddo
         ENDDO
      ENDDO
C-----------------------------------------------------------------------
C     Project tilde D
C-----------------------------------------------------------------------
      CALL EGFPDP ( np, NS, YTR, TEMP, SUMTR, D, AAA )
C-----------------------------------------------------------------------
      RETURN                                                                    
      END                                                                       
C=======================================================================
C=======================================================================
      SUBROUTINE EGFSI3 ( NP, NS, NG, G, DMI, TEMP, YTR, SUMTR,
     &                    D, ITERMX, AUX, AAA )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION G(NP,*), DMI(NP,6,NS), TEMP(NP,NG), YTR(NP,NS), 
     &          SUMTR(NP), D(NP,NS,NS), AAA(NP), AUX(NP,NS)
C-----------------------------------------------------------------------
C     At most two iterations are performed
C-----------------------------------------------------------------------
      NS2 = 2 * NS
      NS3 = 3 * NS
      CALL EGZERO ( NP * NS * NS, D )
      DO I = 1, NS
         II1 = NS3*(I-1) - (I*(I-1))/2 + I
         II2 = NS3*(I-1) - (I*(I-1))/2 + I + NS
         II3 = NS3*(I-1) - (I*(I-1))/2 + I + NS2
         II4 = NS3*(I+NS-1) - ((I+NS)*(I+NS-1))/2 + I + NS
         II5 = NS3*(I+NS-1) - ((I+NS)*(I+NS-1))/2 + I + NS2
         II6 = NS3*(I+NS2-1) - ((I+NS2)*(I+NS2-1))/2 + I + NS2
         do nn = 1, np
            FAC  = 1.0D0 / AUX(nn,I)
            DETM = 1.0D0/(FAC*G(nn,II1)*G(nn,II4)*G(nn,II6)
     &               - FAC*G(nn,II1)*G(nn,II5)*G(nn,II5)
     &               - G(nn,II2)*G(nn,II2)*G(nn,II6))
            DMI(nn,1,I) = ( G(nn,II4) * G(nn,II6) 
     &                    - G(nn,II5) * G(nn,II5) ) * DETM
            DMI(nn,2,I) = - G(nn,II2) * G(nn,II6) * DETM
            DMI(nn,3,I) =   G(nn,II2) * G(nn,II5) * DETM
            D(nn,I,I) = DMI(nn,1,I)
         enddo
      ENDDO
C-----------------------------------------------------------------------
      IF (ITERMX .GT. 1) THEN
         DO I = 1, NS
            III = NS3*(I-1) - (I*(I-1))/2 + I
            do nn = 1, np
               D(nn,I,I) = DMI(nn,1,I) 
     &                   * ( 1.0D0 + G(nn,III) * DMI(nn,1,I)
     &                             * YTR(nn,I) / AUX(nn,I) )
            enddo
            DO J = I+1, NS
               IJ1 = NS3*(I-1) - (I*(I-1))/2 + J
               IJ2 = NS3*(I-1) - (I*(I-1))/2 + J + NS
               IJ3 = NS3*(I-1) - (I*(I-1))/2 + J + NS2
               IJ4 = NS3*(I+NS-1) - ((I+NS)*(I+NS-1))/2 + J + NS
               IJ5 = NS3*(I+NS-1) - ((I+NS)*(I+NS-1))/2 + J + NS2
               IJ6 = NS3*(I+NS2-1) - ((I+NS2)*(I+NS2-1))/2 + J + NS2
               JI2 = NS3*(J-1) - (J*(J-1))/2 + I + NS
               JI3 = NS3*(J-1) - (J*(J-1))/2 + I + NS2
               JI5 = NS3*(J+NS-1) - ((J+NS)*(J+NS-1))/2 + I + NS2
C..............
               do nn = 1, np
                  D(nn,I,J) = 
     &                  - ( DMI(nn,1,I) * DMI(nn,1,J) * G(nn,IJ1)
     &                    + DMI(nn,1,I) * DMI(nn,2,J) * G(nn,IJ2)
     &                    + DMI(nn,1,I) * DMI(nn,3,J) * G(nn,IJ3)
     &                    + DMI(nn,2,I) * DMI(nn,1,J) * G(nn,JI2)
     &                    + DMI(nn,2,I) * DMI(nn,2,J) * G(nn,IJ4)
     &                    + DMI(nn,2,I) * DMI(nn,3,J) * G(nn,IJ5)
     &                    + DMI(nn,3,I) * DMI(nn,1,J) * G(nn,JI3)
     &                    + DMI(nn,3,I) * DMI(nn,2,J) * G(nn,JI5)
     &                    + DMI(nn,3,I) * DMI(nn,3,J) * G(nn,IJ6) )
                  D(nn,J,I) = D(nn,I,J)
               enddo
            ENDDO
         ENDDO
      ENDIF
C-----------------------------------------------------------------------
C     Matrix tilde D
C-----------------------------------------------------------------------
      DO I = 1, NS
         DO J = 1, NS
            do nn = 1, np
               D(nn,I,J) = YTR(nn,i) * D(nn,I,J)
            enddo
         ENDDO
      ENDDO
C-----------------------------------------------------------------------
C     Project tilde D
C-----------------------------------------------------------------------
      CALL EGFPDP ( np, NS, YTR, TEMP, SUMTR, D, AAA )
C-----------------------------------------------------------------------
      RETURN                                                                    
      END                                                                       
C=======================================================================
C=======================================================================
      SUBROUTINE EGFPDP ( NP, NS, YTR, TEMP, SUMTR, D, AAA )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION YTR(np,NS), TEMP(np,NS), D(np,NS,NS), SUMTR(np), AAA(np)
C-----------------------------------------------------------------------
C     Input : matrix D
C     Output (overwritten on D) : PDP with P_ij = Id - YTR*U / sumtr
C-----------------------------------------------------------------------
      CALL EGZERO ( np*NS, TEMP )
      DO J = 1, NS
         DO I = 1, NS
            do nn = 1, np
               TEMP(nn,I) = TEMP(nn,I) + D(nn,I,J) * YTR(nn,j)
            enddo
         ENDDO
      ENDDO
C.....
      DO I = 1, NS
         do nn = 1, np
            AAA(nn) = TEMP(nn,I) / SUMTR(nn)
         enddo
         DO J = 1, NS
            do nn = 1, np
               D(nn,I,J) = D(nn,I,J) - AAA(nn)
            enddo
         ENDDO
      ENDDO
C.....
      CALL EGZERO ( np*NS, TEMP )
      DO J = 1, NS
         DO I = 1, NS
            do nn = 1, np
               TEMP(nn,I) = TEMP(nn,I) + D(nn,J,I) 
            enddo
         ENDDO
      ENDDO
      DO I = 1, NS
         do nn = 1, np
            AAA(nn) = YTR(nn,I) / SUMTR(nn)
         enddo
         DO J = 1, NS
            do nn = 1, np
               D(nn,I,J) = D(nn,I,J) - AAA(nn) * TEMP(nn,J)
            enddo
         ENDDO
      ENDDO
C-----------------------------------------------------------------------
      RETURN                                                                    
      END                                                                       
C=======================================================================
C=======================================================================
      SUBROUTINE EGFDDEC ( NP, NS, XTR, YTR, G, BIN, TEMP )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C
C     Output
C     ------
C        G(NP,*)   Choleski decomposition of the positive
C                  definite version of matrix L_[00]
C                  dimension G(NP,*) at least NP * NS*(NS+1)/2 
C
C-----------------------------------------------------------------------
      DIMENSION G(NP,*), YTR(NP,NS)
C-----------------------------------------------------------------------
      CALL EGFEML00 ( NP, NS, XTR, BIN, G )
C-----------------------------------------------------------------------
C     Make the matrix G positive definite
C-----------------------------------------------------------------------
      NG = NS
      DO I = 1, NS
         II1 = NG*(I-1) - (I*(I-1))/2 + I
         do nn = 1, np
            G(nn,II1) = G(nn,II1) + YTR(nn,I) * YTR(nn,I)
         enddo
         DO J = I+1, NS
            IJ1 = NG*(I-1) - (I*(I-1))/2 + J
            do nn = 1, np
               G(nn,IJ1) = G(nn,IJ1) + YTR(nn,I) * YTR(nn,J)
            enddo
         ENDDO
      ENDDO
C-----------------------------------------------------------------------
      CALL EGFDEC ( NP, NG, G, TEMP, IER )
      IF ( IER .NE. 0 ) THEN
         WRITE(*,'(1X,''stopping in EGFDDEC with IER = '',i3)') IER
         STOP
      ENDIF
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGFDEC ( NP, N, A, W, IER )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C     A is stored in symmetric form
C     A(nn,i,j) --> A(nn,ij) with ij = N*(i-1) - i(i-1)/2 + j for j >= i
C-----------------------------------------------------------------------
      DIMENSION A(NP,*), W(NP,N)
      IER = 0
      DO  K = 1,N
         KM1 = K - 1
         KP1 = K + 1
         DO  J = 1, KM1
            jj = n*(j-1) - (j*(j-1))/2 + j
            kj = n*(j-1) - (j*(j-1))/2 + k
            do nn = 1, np
               W(nn,J)= A(nn,JJ)*A(nn,KJ)
            enddo
         ENDDO
         kk = n*(k-1) - (k*(k-1))/2 + k
         DO J = 1, KM1
            kj = n*(j-1) - (j*(j-1))/2 + k
            do nn = 1, np
               A(nn,KK) = A(nn,KK) - W(nn,J)*A(nn,KJ)
            enddo
         ENDDO
         do nn = 1, np
            IF (A(nn,KK) .EQ. 0.0D0) THEN
               WRITE(6, '(''SINGULAR MATRIX IN EGFDEC'')' )
               WRITE(6, '(''NODE = '',I7)' ) nn
               IER = K
               stop
            ENDIF
         enddo
         DO J = 1, KM1
            DO I = KP1, N
               ik = n*(k-1) - (k*(k-1))/2 + i
               ij = n*(j-1) - (j*(j-1))/2 + i
               do nn = 1, np
                  A(nn,IK) = A(nn,IK) - A(nn,IJ)*W(nn,J)
               enddo
            ENDDO
         ENDDO
         DO I = KP1, N
            ik = n*(k-1) - (k*(k-1))/2 + i
            do nn = 1, np
               A(nn,IK) = A(nn,IK)/A(nn,KK)
            enddo
         ENDDO
      ENDDO
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGFSOL ( NP, N, A, B, FAC )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(NP,*), B(NP,N), FAC(NP)
      NM1 = N - 1
      DO J = 1, NM1
         JP1 = J + 1
         do nn = 1, np
            FAC(nn) = -B(nn,J)
         enddo
         DO  K = JP1, N
            kj = n*(j-1) - (j*(j-1))/2 + k
            do nn = 1, np
               B(nn,K) = B(nn,K) + A(nn,KJ)*FAC(nn)
            enddo
         ENDDO
      ENDDO
      DO J = 1, N
         jj = n*(j-1) - (j*(j-1))/2 + j
         do nn = 1, np
            B(nn,J) = B(nn,J)/A(nn,JJ)
         enddo
      ENDDO
      DO JB = 1, NM1
         J = N + 1 - JB
         JM1 = J - 1
         do nn = 1, np
            FAC(nn) = -B(nn,J)
         enddo
         DO K = 1, JM1
            jk = n*(k-1) - (k*(k-1))/2 + j
            do nn = 1, np
               B(nn,K) = B(nn,K) + A(nn,JK)*FAC(nn)
            enddo
         ENDDO
      ENDDO
      RETURN
      END
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
      SUBROUTINE EGFEMH ( NP, NS, 
     &                    T, DLT1, DLT2, DLT3, DLT4, DLT5, DLT6, 
     &                    X, WT, RU, PATMOS, 
     &                    BETA, G, AIJ, FITA, ETA, BIN )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      DIMENSION X(NP,NS), BIN(NP,*), WT(NS), AIJ(NP), ETA(NP,NS),
     &          FITA(7,NS,NS) 
      DIMENSION T(NP), DLT1(NP), DLT2(NP), DLT3(NP), DLT4(NP), DLT5(NP),
     &          DLT6(NP)
      DIMENSION G(NP,*), BETA(NP,NS)
C-----------------------------------------------------------------------
C         EVALUATE THE RHS BETA
C-----------------------------------------------------------------------
      DO K = 1, NS
         do nn = 1, np
            BETA(nn,K) = X(nn,K)
         enddo
      ENDDO
C-----------------------------------------------------------------------
C         EVALUATE THE MATRIX H
C-----------------------------------------------------------------------
      DO I = 1, NS
         III = NS*(I-1) - (I*(I-1))/2 + I
         do nn = 1, np
            G(nn,III) = X(nn,I) * X(nn,I) / ETA(nn,I)
         enddo
      ENDDO
      FFF = 6.0D0 * RU / ( 5.0D0 * PATMOS )
      ZZZ = 5.0D0 / 3.0D0
      DO I = 1, NS
         IP1 = I + 1
         III = NS*(I-1) - (I*(I-1))/2 + I
         DO J = IP1, NS
            IND = NS*(I-1) - (I*(I-1))/2 + J
            JJJ = NS*(J-1) - (J*(J-1))/2 + J
            FAC = FFF / ( WT(I)+WT(J) )
            WIJ = WT(I) / WT(J)
            WJI = WT(J) / WT(I)
C...........
            FITA0 = FITA(1,I,J)
            FITA1 = FITA(2,I,J)
            FITA2 = FITA(3,I,J)
            FITA3 = FITA(4,I,J)
            FITA4 = FITA(5,I,J)
            FITA5 = FITA(6,I,J)
            FITA6 = FITA(7,I,J)
            do nn = 1, np
               AIJ(nn) = FITA0          + FITA1*DLT1(nn) 
     &                 + FITA2*DLT2(nn) + FITA3*DLT3(nn)
     &                 + FITA4*DLT4(nn) + FITA5*DLT5(nn)
     &                 + FITA6*DLT6(nn) 
            enddo
C...........
            do nn = 1, np
               AAA = X(nn,I) * X(nn,J) * FAC * BIN(nn,IND) * T(nn)
               G(nn,IND) = AAA * ( AIJ(nn) - ZZZ )
               G(nn,III) = G(nn,III) + AAA * 
     &                      ( AIJ(nn)*WJI + ZZZ )
               G(nn,JJJ) = G(nn,JJJ) + AAA * 
     &                      ( AIJ(nn)*WIJ + ZZZ )
            enddo
         ENDDO
      ENDDO
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGFEMHA ( NP, NS, T, X, WT, RU, PATMOS, 
     &                     BETA, G, AIJ, ETA, BIN )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      DIMENSION X(NP,NS), BIN(NP,*), WT(NS), AIJ(NS,NS), 
     &          ETA(NP,NS), T(NP)
      DIMENSION G(NP,*), BETA(NP,NS)
C-----------------------------------------------------------------------
C         EVALUATE THE RHS BETA
C-----------------------------------------------------------------------
      DO K = 1, NS
         do nn = 1, np
            BETA(nn,K) = X(nn,K)
         enddo
      ENDDO
C-----------------------------------------------------------------------
C         EVALUATE THE MATRIX H
C         AIJ is node independent
C-----------------------------------------------------------------------
      DO I = 1, NS
         III = NS*(I-1) - (I*(I-1))/2 + I
         do nn = 1, np
            G(nn,III) = X(nn,I) * X(nn,I) / ETA(nn,I)
         enddo
      ENDDO
      FFF = 6.0D0 * RU / ( 5.0D0 * PATMOS )
      ZZZ = 5.0D0 / 3.0D0
      DO I = 1, NS
         IP1 = I + 1
         III = NS*(I-1) - (I*(I-1))/2 + I
         DO J = IP1, NS
            IND = NS*(I-1) - (I*(I-1))/2 + J
            JJJ = NS*(J-1) - (J*(J-1))/2 + J
            FAC = FFF / ( WT(I)+WT(J) )
            WIJ = WT(I) / WT(J)
            WJI = WT(J) / WT(I)
            A1  = ( AIJ(I,J) - ZZZ )
            A2  = ( AIJ(I,J)*WJI + ZZZ )
            A3  = ( AIJ(I,J)*WIJ + ZZZ )
            do nn = 1, np
               AAA = X(nn,I) * X(nn,J) * FAC * BIN(nn,IND) * T(nn)
               G(nn,IND) = AAA * A1
               G(nn,III) = G(nn,III) + AAA * A2
               G(nn,JJJ) = G(nn,JJJ) + AAA * A3
            enddo
         ENDDO
      ENDDO
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGFEMK ( NP, NS, 
     &                    TEMPER, DLT1, DLT2, DLT3, DLT4, DLT5, DLT6, 
     &                    X, WT, RU, PATMOS,
     &                    BIN, AIJ, FITA, CINT, ETA, CXI,
     &                    G, BETA, CCC, RCV )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      DIMENSION G(NP,*), BETA(NP,*), CCC(NP), RCV(NP)
      DIMENSION BIN(NP,*), AIJ(NP), CINT(NP,NS), ETA(NP,NS), 
     &          X(NP,NS), WT(NS), CXI(NP,NS), FITA(7,NS,NS)
      DIMENSION TEMPER(NP), DLT1(NP), DLT2(NP), DLT3(NP), DLT4(NP), 
     &          DLT5(NP), DLT6(NP)
C-----------------------------------------------------------------------
C         Form the RHS 
C-----------------------------------------------------------------------
      call egzero ( np, ccc )
      DO I = 1, NS
         do nn = 1, np
            CCC(nn) = CCC(nn) + X(nn,I) * CINT(nn,I)
         enddo
      ENDDO
      do nn = 1, np
         rcv(nn) = 1.0D0 / ( 1.5D0 + CCC(nn) )
      enddo
      DO I = 1, NS
         do nn = 1, np
            BETA(nn,I)      = X(nn,I) * CCC(nn) * RCV(nn)
            BETA(nn,I+NS)   = - X(nn,I) * CINT(nn,I) * RCV(nn)
         enddo
      ENDDO
C-----------------------------------------------------------------------
C         Form the transport linear system matrix K
C-----------------------------------------------------------------------
      NS2 = 2 * NS
C ----- Initialize the diagonals
      DO I = 1, NS
         II1 = NS2*(I-1)     - (I*(I-1))/2           + I
         II2 = NS2*(I-1)     - (I*(I-1))/2           + I + NS
         II3 = NS2*(I+NS-1)  - ((I+NS)*(I+NS-1))/2   + I + NS
         do nn = 1, np
            FAC = 4.0D0 * CXI(nn,I) / ETA(nn,I) * X(nn,I) * X(nn,I)
            G(nn,II1) = FAC 
            G(nn,II2) = - FAC 
            G(nn,II3) = FAC 
         enddo
      ENDDO
C-----------------------------------------------------------------------
      FFF = 2.4D0 * RU / PATMOS
      DO I = 1, NS
         II1 = NS2*(I-1) - (I*(I-1))/2 + I
         II2 = NS2*(I-1) - (I*(I-1))/2 + I + NS
         II3 = NS2*(I+NS-1) - ((I+NS)*(I+NS-1))/2 + I + NS
C........
         IP1 = I + 1
         DO J = IP1, NS
            IND = NS*(I-1) - (I*(I-1))/2 + J
C
            JJ1 = NS2*(J-1) - (J*(J-1))/2 + J
            JJ2 = NS2*(J-1) - (J*(J-1))/2 + J + NS
            JJ3 = NS2*(J+NS-1) - ((J+NS)*(J+NS-1))/2 + J + NS
C...........
            IJ1 = NS2*(I-1) - (I*(I-1))/2 + J
            IJ2 = NS2*(I-1) - (I*(I-1))/2 + J + NS
            IJ3 = NS2*(I+NS-1) - ((I+NS)*(I+NS-1))/2 + J + NS
C...........
            JI2 = NS2*(J-1) - (J*(J-1))/2 + I + NS
C...........
            Z2  = FFF / ( WT(I)+WT(J) )
            Z1  = 1.25D0 * Z2
            ZI  = FFF / WT(I) 
            ZJ  = FFF / WT(J) 
            ZIJ = ZI + ZJ
            WJI = WT(J) / WT(I)
            WIJ = WT(I) / WT(J)
C...........
            FITA0 = FITA(1,I,J)
            FITA1 = FITA(2,I,J)
            FITA2 = FITA(3,I,J)
            FITA3 = FITA(4,I,J)
            FITA4 = FITA(5,I,J)
            FITA5 = FITA(6,I,J)
            FITA6 = FITA(7,I,J)
            do nn = 1, np
               AIJ(nn) = FITA0          + FITA1*DLT1(nn) 
     &                 + FITA2*DLT2(nn) + FITA3*DLT3(nn)
     &                 + FITA4*DLT4(nn) + FITA5*DLT5(nn)
     &                 + FITA6*DLT6(nn) 
            enddo
C...........
            do nn = 1, np
               AAA = X(nn,I) * X(nn,J) * BIN(nn,IND) * TEMPER(nn) * Z1
               BB1 = X(nn,I) * X(nn,J) * BIN(nn,IND) * AIJ(nn)
     &             * TEMPER(nn) 
               BBB = BB1 * Z2
               BBJ = BBB * WJI
               BBI = BBB * WIJ
               BB2 = BB1 * ZIJ
               CIJ = ( CXI(nn,I) + CXI(nn,J) )
C ----- 10,10
               G(nn,IJ1) = - AAA + BBB * CIJ
               G(nn,II1) = G(nn,II1) + AAA + BBJ * CIJ
               G(nn,JJ1) = G(nn,JJ1) + AAA + BBI * CIJ
C ----- 10,01 and 01,10
               G(nn,IJ2) = - BB1 * CXI(nn,J) * ZI
               G(nn,JI2) = - BB1 * CXI(nn,I) * ZJ
               G(nn,II2) = G(nn,II2) - BB1 * CXI(nn,I) * ZI
               G(nn,JJ2) = G(nn,JJ2) - BB1 * CXI(nn,J) * ZJ
C ----- 01,01
               G(nn,IJ3) = 0.0D0
               G(nn,II3) = G(nn,II3) + BB2 * CXI(nn,I)
               G(nn,JJ3) = G(nn,JJ3) + BB2 * CXI(nn,J)
            enddo
         ENDDO
      ENDDO
      DO I = NS+1, NS2
         II3 = NS2*(I-1) - (I*(I-1))/2 + I
         IF ( CXI(1,I-NS) .EQ. 0.0D0) THEN
            do nn = 1, np
               G(nn,II3) = 1.0D0
            enddo
         ENDIF
      ENDDO
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGFEMK01 ( NP, NS, 
     &                      TEMPER, DLT1, DLT2, DLT3, DLT4, DLT5, DLT6, 
     &                      X, WT, RU, PATMOS,
     &                      BIN, AIJ, FITA, CINT, ETA, CXI,
     &                      G, BETA, CCC, RCV )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      DIMENSION G(NP,*), BETA(NP,*), CCC(NP), RCV(NP)
      DIMENSION BIN(NP,*), AIJ(NP), CINT(NP,NS), ETA(NP,NS), 
     &          X(NP,NS), WT(NS), CXI(NP,NS), FITA(7,NS,NS)
      DIMENSION TEMPER(NP), DLT1(NP), DLT2(NP), DLT3(NP), DLT4(NP), 
     &          DLT5(NP), DLT6(NP)
C-----------------------------------------------------------------------
C         Form the RHS 
C-----------------------------------------------------------------------
      call egzero ( np, ccc )
      DO I = 1, NS
         do nn = 1, np
            CCC(nn) = CCC(nn) + X(nn,I) * CINT(nn,I)
         enddo
      ENDDO
      do nn = 1, np
         rcv(nn) = 1.0D0 / ( 1.5D0 + CCC(nn) )
      enddo
      DO I = 1, NS
         do nn = 1, np
            BETA(nn,I)   = - X(nn,I) * CINT(nn,I) * RCV(nn)
         enddo
      ENDDO
C-----------------------------------------------------------------------
C         Form the transport linear system matrix K_[01]
C-----------------------------------------------------------------------
C ----- Initialize the diagonals
      DO I = 1, NS
         II1 = NS*(I-1) - (I*(I-1))/2 + I
         do nn = 1, np
            FAC = 4.0D0 * CXI(nn,I) / ETA(nn,I) * X(nn,I) * X(nn,I)
            G(nn,II1) = FAC 
         enddo
      ENDDO
C-----------------------------------------------------------------------
      FFF = 2.4D0 * RU / PATMOS
      DO I = 1, NS
         II1 = NS*(I-1) - (I*(I-1))/2 + I
C........
         IP1 = I + 1
         DO J = IP1, NS
            JJ1 = NS*(J-1) - (J*(J-1))/2 + J
            IJ1 = NS*(I-1) - (I*(I-1))/2 + J
C...........
            ZZZ = FFF * ( WT(I) + WT(J) ) / ( WT(I) * WT(J) )
C...........
            FITA0 = FITA(1,I,J)
            FITA1 = FITA(2,I,J)
            FITA2 = FITA(3,I,J)
            FITA3 = FITA(4,I,J)
            FITA4 = FITA(5,I,J)
            FITA5 = FITA(6,I,J)
            FITA6 = FITA(7,I,J)
            do nn = 1, np
               AIJ(nn) = FITA0          + FITA1*DLT1(nn) 
     &                 + FITA2*DLT2(nn) + FITA3*DLT3(nn)
     &                 + FITA4*DLT4(nn) + FITA5*DLT5(nn)
     &                 + FITA6*DLT6(nn) 
            enddo
C...........
            do nn = 1, np
               BB2  = X(nn,I) * X(nn,J) * BIN(nn,IJ1) * TEMPER(nn)
     &              * AIJ(nn) * ZZZ
C ----- 01,01
               G(nn,IJ1) = 0.0D0
               G(nn,II1) = G(nn,II1) + BB2 * CXI(nn,I)
               G(nn,JJ1) = G(nn,JJ1) + BB2 * CXI(nn,J)
            enddo
         ENDDO
      ENDDO
      DO I = 1, NS
         II1 = NS*(I-1) - (I*(I-1))/2 + I
         IF ( CXI(1,I) .EQ. 0.0D0) THEN
            do nn = 1, np
               G(nn,II1) = 1.0D0
            enddo
         ENDIF
      ENDDO
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGFEMKA ( NP, NS, TEMPER, X, WT, RU, PATMOS,
     &                     BIN, AIJ, CINT, ETA, CXI,
     &                     G, BETA, CCC, RCV )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      DIMENSION G(NP,*), BETA(NP,*), CCC(NP), RCV(NP)
      DIMENSION BIN(NP,*), AIJ(NS,NS), CINT(NP,NS), ETA(NP,NS), 
     &          X(NP,NS), WT(NS), CXI(NP,NS), TEMPER(NP)
C-----------------------------------------------------------------------
C         Form the RHS 
C-----------------------------------------------------------------------
      call egzero ( np, ccc )
      DO I = 1, NS
         do nn = 1, np
            CCC(nn) = CCC(nn) + X(nn,I) * CINT(nn,I)
         enddo
      ENDDO
      do nn = 1, np
         rcv(nn) = 1.0D0 / ( 1.5D0 + CCC(nn) )
      enddo
      DO I = 1, NS
         do nn = 1, np
            BETA(nn,I)      = X(nn,I) * CCC(nn) * RCV(nn)
            BETA(nn,I+NS)   = - X(nn,I) * CINT(nn,I) * RCV(nn)
         enddo
      ENDDO
C-----------------------------------------------------------------------
C         Form the transport linear system matrix K
C         A_ij is node independent
C-----------------------------------------------------------------------
      NS2 = 2 * NS
C ----- Initialize the diagonals
      DO I = 1, NS
         II1 = NS2*(I-1)     - (I*(I-1))/2           + I
         II2 = NS2*(I-1)     - (I*(I-1))/2           + I + NS
         II3 = NS2*(I+NS-1)  - ((I+NS)*(I+NS-1))/2   + I + NS
         do nn = 1, np
            FAC = 4.0D0 * CXI(nn,I) / ETA(nn,I) * X(nn,I) * X(nn,I)
            G(nn,II1) = FAC 
            G(nn,II2) = - FAC 
            G(nn,II3) = FAC 
         enddo
      ENDDO
C-----------------------------------------------------------------------
      FFF = 2.4D0 * RU / PATMOS
      DO I = 1, NS
         II1 = NS2*(I-1) - (I*(I-1))/2 + I
         II2 = NS2*(I-1) - (I*(I-1))/2 + I + NS
         II3 = NS2*(I+NS-1) - ((I+NS)*(I+NS-1))/2 + I + NS
C........
         IP1 = I + 1
         DO J = IP1, NS
            IND = NS*(I-1) - (I*(I-1))/2 + J
C
            JJ1 = NS2*(J-1) - (J*(J-1))/2 + J
            JJ2 = NS2*(J-1) - (J*(J-1))/2 + J + NS
            JJ3 = NS2*(J+NS-1) - ((J+NS)*(J+NS-1))/2 + J + NS
C...........
            IJ1 = NS2*(I-1) - (I*(I-1))/2 + J
            IJ2 = NS2*(I-1) - (I*(I-1))/2 + J + NS
            IJ3 = NS2*(I+NS-1) - ((I+NS)*(I+NS-1))/2 + J + NS
C...........
            JI2 = NS2*(J-1) - (J*(J-1))/2 + I + NS
C...........
            Z1  = 1.25D0 * FFF / ( WT(I)+WT(J) )
            Z2  = FFF * AIJ(I,J) / ( WT(I)+WT(J) )
            ZI  = FFF * AIJ(I,J) / WT(I)
            ZJ  = FFF * AIJ(I,J) / WT(J)
            ZIJ = ZI + ZJ
            WJI = WT(J) / WT(I)
            WIJ = WT(I) / WT(J)
            do nn = 1, np
               BB1 = X(nn,I) * X(nn,J) * BIN(nn,IND) * TEMPER(nn) 
               AAA = BB1 * Z1
               BBB = BB1 * Z2
               BBJ = BBB * WJI
               BBI = BBB * WIJ
               BB2 = BB1 * ZIJ
               CIJ = ( CXI(nn,I) + CXI(nn,J) )
C ----- 10,10
               G(nn,IJ1) = - AAA + BBB * CIJ
               G(nn,II1) = G(nn,II1) + AAA + BBJ * CIJ
               G(nn,JJ1) = G(nn,JJ1) + AAA + BBI * CIJ
C ----- 10,01 and 01,10
               G(nn,IJ2) = - BB1 * CXI(nn,J) * ZI
               G(nn,JI2) = - BB1 * CXI(nn,I) * ZJ
               G(nn,II2) = G(nn,II2) - BB1 * CXI(nn,I) * ZI
               G(nn,JJ2) = G(nn,JJ2) - BB1 * CXI(nn,J) * ZJ
C ----- 01,01
               G(nn,IJ3) = 0.0D0
               G(nn,II3) = G(nn,II3) + BB2 * CXI(nn,I)
               G(nn,JJ3) = G(nn,JJ3) + BB2 * CXI(nn,J)
            enddo
         ENDDO
      ENDDO
      DO I = NS+1, NS2
         II3 = NS2*(I-1) - (I*(I-1))/2 + I
         IF ( CXI(1,I-NS) .EQ. 0.0D0) THEN
            do nn = 1, np
               G(nn,II3) = 1.0D0
            enddo
         ENDIF
      ENDDO
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGFEMKA01 ( NP, NS, TEMPER, X, WT, RU, PATMOS,
     &                       BIN, AIJ, CINT, ETA, CXI,
     &                       G, BETA, CCC, RCV )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      DIMENSION G(NP,*), BETA(NP,*), CCC(NP), RCV(NP)
      DIMENSION BIN(NP,*), AIJ(NS,NS), CINT(NP,NS), ETA(NP,NS), 
     &          X(NP,NS), WT(NS), CXI(NP,NS), TEMPER(NP)
C-----------------------------------------------------------------------
C         Form the RHS 
C-----------------------------------------------------------------------
      call egzero ( np, ccc )
      DO I = 1, NS
         do nn = 1, np
            CCC(nn) = CCC(nn) + X(nn,I) * CINT(nn,I)
         enddo
      ENDDO
      do nn = 1, np
         rcv(nn) = 1.0D0 / ( 1.5D0 + CCC(nn) )
      enddo
      DO I = 1, NS
         do nn = 1, np
            BETA(nn,I)   = - X(nn,I) * CINT(nn,I) * RCV(nn)
         enddo
      ENDDO
C-----------------------------------------------------------------------
C         Form the transport linear system matrix K_[01]
C-----------------------------------------------------------------------
C ----- Initialize the diagonals
      DO I = 1, NS
         II1 = NS*(I-1) - (I*(I-1))/2 + I
         do nn = 1, np
            FAC = 4.0D0 * CXI(nn,I) / ETA(nn,I) * X(nn,I) * X(nn,I)
            G(nn,II1) = FAC 
         enddo
      ENDDO
C-----------------------------------------------------------------------
      FFF = 2.4D0 * RU / PATMOS
      DO I = 1, NS
         II1 = NS*(I-1) - (I*(I-1))/2 + I
C........
         IP1 = I + 1
         DO J = IP1, NS
            JJ1 = NS*(J-1) - (J*(J-1))/2 + J
            IJ1 = NS*(I-1) - (I*(I-1))/2 + J
C...........
            ZZZ = FFF * AIJ(I,J) * ( WT(I) + WT(J) ) / ( WT(I) * WT(J) )
            do nn = 1, np
               BB2  = X(nn,I) * X(nn,J) * BIN(nn,IJ1) * TEMPER(nn) * ZZZ
C ----- 01,01
               G(nn,IJ1) = 0.0D0
               G(nn,II1) = G(nn,II1) + BB2 * CXI(nn,I)
               G(nn,JJ1) = G(nn,JJ1) + BB2 * CXI(nn,J)
            enddo
         ENDDO
      ENDDO
      DO I = 1, NS
         II1 = NS*(I-1) - (I*(I-1))/2 + I
         IF ( CXI(1,I) .EQ. 0.0D0) THEN
            do nn = 1, np
               G(nn,II1) = 1.0D0
            enddo
         ENDIF
      ENDDO
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGFEML ( NP, NS, 
     &                    TEMPER, DLT1, DLT2, DLT3, DLT4, DLT5, DLT6, 
     &                    X, WT, RU, PATMOS,
     &                    BIN, AIJ, BIJ, CIJ, FITA, FITB, FITC,
     &                    CINT, ETA, CXI, LIN, G, BETA )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      DIMENSION G(NP,*), BETA(NP,*)
      DIMENSION BIN(NP,*), AIJ(NP), BIJ(NP), 
     &          CIJ(NP), CINT(NP,NS), ETA(NP,NS), X(NP,NS), 
     &          WT(NS), CXI(NP,NS), LIN(NS)
      DIMENSION TEMPER(NP), DLT1(NP), DLT2(NP), DLT3(NP), DLT4(NP), 
     &          DLT5(NP), DLT6(NP)
      DIMENSION FITA(7,NS,NS), FITB(7,NS,NS), FITC(7,NS,NS)
C-----------------------------------------------------------------------
C         Form the RHS 
C-----------------------------------------------------------------------
      DO I = 1, NS
         do nn = 1, np
            BETA(nn,I)      = 0.0D0
            BETA(nn,I+NS)   = 2.5D0 * X(nn,I)
            BETA(nn,I+2*NS) = CINT(nn,I) * X(nn,I)
         enddo
      ENDDO
C-----------------------------------------------------------------------
C         Form the transport linear system matrix L
C
C     Note: the system matrix is evaluated at atmospheric pressure
C-----------------------------------------------------------------------
      NS2 = 2 * NS
      NS3 = 3 * NS
C ----- Initialize the diagonals
      DO I = 1, NS
         IND = NS*(I-1) - (I*(I-1))/2 + I
C
         II1 = NS3*(I-1)     - (I*(I-1))/2           + I
         II2 = NS3*(I-1)     - (I*(I-1))/2           + I + NS
         II3 = NS3*(I-1)     - (I*(I-1))/2           + I + NS2
         II4 = NS3*(I+NS-1)  - ((I+NS)*(I+NS-1))/2   + I + NS
         II5 = NS3*(I+NS-1)  - ((I+NS)*(I+NS-1))/2   + I + NS2
         II6 = NS3*(I+NS2-1) - ((I+NS2)*(I+NS2-1))/2 + I + NS2
         ZZZ = 5.0D0 * PATMOS * WT(I) / ( 3.0D0 * RU )
         do nn = 1, np
            G(nn,II1) = 0.0D0
            G(nn,II2) = 0.0D0
            G(nn,II3) = 0.0D0
C-----------------------------------------------------------------------
C        FAC = 2 * A_{ii} / {\cal D}_{ii}
C-----------------------------------------------------------------------
            FAC = ZZZ / ( TEMPER(nn) * ETA(nn,I) )
            G(nn,II4) = FAC * X(nn,I) * X(nn,I) 
     &                   * ( 1.0D0 + 1.0D1 * CXI(nn,I) / 3.0D0 )
            G(nn,II5) = - FAC * X(nn,I) * X(nn,I) 
     &                   * 2.0D0 * CXI(nn,I) 
            G(nn,II6) = FAC * X(nn,I) * X(nn,I) 
     &                   * 1.2D0 * CXI(nn,I)
     &                   + X(nn,I) * X(nn,I) * CINT(nn,I) * BIN(nn,IND)
         enddo
      ENDDO
C-----------------------------------------------------------------------
      FFF = 2.0D1 / 3.0D0
      DO I = 1, NS
         II1 = NS3*(I-1) - (I*(I-1))/2 + I
         II2 = NS3*(I-1) - (I*(I-1))/2 + I + NS
         II3 = NS3*(I-1) - (I*(I-1))/2 + I + NS2
         II4 = NS3*(I+NS-1) - ((I+NS)*(I+NS-1))/2 + I + NS
         II5 = NS3*(I+NS-1) - ((I+NS)*(I+NS-1))/2 + I + NS2
         II6 = NS3*(I+NS2-1) - ((I+NS2)*(I+NS2-1))/2 + I + NS2
C........
         IP1 = I + 1
         DO J = IP1, NS
         IND = NS*(I-1) - (I*(I-1))/2 + J
C
            JJ1 = NS3*(J-1) - (J*(J-1))/2 + J
            JJ2 = NS3*(J-1) - (J*(J-1))/2 + J + NS
            JJ3 = NS3*(J-1) - (J*(J-1))/2 + J + NS2
            JJ4 = NS3*(J+NS-1) - ((J+NS)*(J+NS-1))/2 + J + NS
            JJ5 = NS3*(J+NS-1) - ((J+NS)*(J+NS-1))/2 + J + NS2
            JJ6 = NS3*(J+NS2-1) - ((J+NS2)*(J+NS2-1))/2 + J + NS2
C...........
            IJ1 = NS3*(I-1) - (I*(I-1))/2 + J
            IJ2 = NS3*(I-1) - (I*(I-1))/2 + J + NS
            IJ3 = NS3*(I-1) - (I*(I-1))/2 + J + NS2
            IJ4 = NS3*(I+NS-1) - ((I+NS)*(I+NS-1))/2 + J + NS
            IJ5 = NS3*(I+NS-1) - ((I+NS)*(I+NS-1))/2 + J + NS2
            IJ6 = NS3*(I+NS2-1) - ((I+NS2)*(I+NS2-1))/2 + J + NS2
C...........
            JI2 = NS3*(J-1) - (J*(J-1))/2 + I + NS
            JI3 = NS3*(J-1) - (J*(J-1))/2 + I + NS2
            JI5 = NS3*(J+NS-1) - ((J+NS)*(J+NS-1))/2 + I + NS2
C...........
            WWI = WT(I) / ( WT(I) + WT(J) )
            WWJ = WT(J) / ( WT(I) + WT(J) )
            WIJ = WT(I) / WT(J)
            WJI = WT(J) / WT(I)
C..............
            FITA0 = FITA(1,I,J)
            FITA1 = FITA(2,I,J)
            FITA2 = FITA(3,I,J)
            FITA3 = FITA(4,I,J)
            FITA4 = FITA(5,I,J)
            FITA5 = FITA(6,I,J)
            FITA6 = FITA(7,I,J)
C..............
            FITB0 = FITB(1,I,J)
            FITB1 = FITB(2,I,J)
            FITB2 = FITB(3,I,J)
            FITB3 = FITB(4,I,J)
            FITB4 = FITB(5,I,J)
            FITB5 = FITB(6,I,J)
            FITB6 = FITB(7,I,J)
C..............
            FITC0 = FITC(1,I,J)
            FITC1 = FITC(2,I,J)
            FITC2 = FITC(3,I,J)
            FITC3 = FITC(4,I,J)
            FITC4 = FITC(5,I,J)
            FITC5 = FITC(6,I,J)
            FITC6 = FITC(7,I,J)
            do nn = 1, np
               AIJ(nn) = FITA0          + FITA1*DLT1(nn) 
     &                 + FITA2*DLT2(nn) + FITA3*DLT3(nn)
     &                 + FITA4*DLT4(nn) + FITA5*DLT5(nn)
     &                 + FITA6*DLT6(nn) 
               BIJ(nn) = FITB0          + FITB1*DLT1(nn) 
     &                 + FITB2*DLT2(nn) + FITB3*DLT3(nn)
     &                 + FITB4*DLT4(nn) + FITB5*DLT5(nn)
     &                 + FITB6*DLT6(nn) 
               CIJ(nn) = FITC0          + FITC1*DLT1(nn) 
     &                 + FITC2*DLT2(nn) + FITC3*DLT3(nn)
     &                 + FITC4*DLT4(nn) + FITC5*DLT5(nn)
     &                 + FITC6*DLT6(nn) 
            enddo
C..............
            do nn = 1, np
C ----- 00,00
               AAA = X(nn,I) * X(nn,J) * BIN(nn,IND)
               G(nn,IJ1) = - AAA
               G(nn,II1) = G(nn,II1) + AAA
               G(nn,JJ1) = G(nn,JJ1) + AAA
C ----- 00,10 and 10,00
               AAI = 0.5D0 * WWI * (6.0D0*CIJ(nn) - 5.0D0)
               AAJ = 0.5D0 * WWJ * (6.0D0*CIJ(nn) - 5.0D0)
               G(nn,IJ2) = AAI * AAA
               G(nn,JI2) = AAJ * AAA
               G(nn,II2) = G(nn,II2) - AAJ * AAA
               G(nn,JJ2) = G(nn,JJ2) - AAI * AAA
C ----- 00,01 and 01,00
               G(nn,IJ3) = 0.0D0
               G(nn,JI3) = 0.0D0
C ----- 10,10
               BB1 = FFF * AIJ(nn) 
     &                              * (CXI(nn,I)+CXI(nn,J))
               BBB = WWI * WWJ * (
     &               13.75D0 - 3.0D0*BIJ(nn) 
     &                - 4.0D0*AIJ(nn) - BB1 )
               BBJ = ( 7.5D0 * WWI * WWI + 
     &               (6.25D0 - 3.0D0*BIJ(nn)) * WWJ * WWJ 
     &                + 4.0D0*AIJ(nn) * WWI * WWJ
     &                + WWI * WWJ * BB1 )
               BBI = ( 7.5D0 * WWJ * WWJ + 
     &               (6.25D0 - 3.0D0*BIJ(nn)) * WWI * WWI 
     &                + 4.0D0*AIJ(nn) * WWI * WWJ
     &                + WWI * WWJ * BB1 )
               G(nn,IJ4) = - BBB * AAA
               G(nn,II4) = G(nn,II4) + BBJ * AAA
               G(nn,JJ4) = G(nn,JJ4) + BBI * AAA
C ----- 10,01 and 01,10
               CCI = WWI * 4.0D0 * AIJ(nn) * CXI(nn,I)
               CCJ = WWJ * 4.0D0 * AIJ(nn) * CXI(nn,J)
               G(nn,IJ5) = - CCJ * AAA
               G(nn,JI5) = - CCI * AAA
               G(nn,II5) = G(nn,II5) - CCI * AAA
               G(nn,JJ5) = G(nn,JJ5) - CCJ * AAA
C ----- 01,01
               DDJ = AAA * CINT(nn,I) 
     &            +  AAA * 2.4D0 * AIJ(nn) * WIJ * CXI(nn,I)
               DDI = AAA * CINT(nn,J) 
     &            +  AAA * 2.4D0 * AIJ(nn) * WJI * CXI(nn,J)
               G(nn,IJ6) = 0.0D0
               G(nn,II6) = G(nn,II6) + DDJ
               G(nn,JJ6) = G(nn,JJ6) + DDI
            enddo
         ENDDO
      ENDDO
      DO I = NS2+1, NS3
         II6 = NS3*(I-1) - (I*(I-1))/2 + I
         IF ( LIN(I-NS2) .EQ. 0) THEN
            do nn = 1, np
               G(nn,II6) = 1.0D0
            enddo
         ENDIF
      ENDDO
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGFEMA ( NP, NS, 
     &                    TEMPER, DLT1, DLT2, DLT3, DLT4, DLT5, DLT6, 
     &                    X, WT, RU, PATMOS,
     &                    BIN, AIJ, BIJ, FITA, FITB, 
     &                    CINT, ETA, CXI, LIN, G, BETA )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      DIMENSION G(NP,*), BETA(NP,*)
      DIMENSION BIN(NP,*), AIJ(NP), BIJ(NP), 
     &          CINT(NP,NS), ETA(NP,NS), X(NP,NS), 
     &          WT(NS), CXI(NP,NS), LIN(NS)
      DIMENSION TEMPER(NP), DLT1(NP), DLT2(NP), DLT3(NP), DLT4(NP), 
     &          DLT5(NP), DLT6(NP)
      DIMENSION FITA(7,NS,NS), FITB(7,NS,NS)
C-----------------------------------------------------------------------
C         Form the RHS 
C-----------------------------------------------------------------------
      DO I = 1, NS
         do nn = 1, np
            BETA(nn,I)    = 2.5D0 * X(nn,I)
            BETA(nn,I+NS) = CINT(nn,I) * X(nn,I)
         enddo
      ENDDO
C-----------------------------------------------------------------------
C         Form the transport linear system matrix \Lambda
C
C     Note: the system matrix is evaluated at atmospheric pressure
C-----------------------------------------------------------------------
      NS2 = 2 * NS
C ----- Initialize the diagonals
      DO I = 1, NS
         IND = NS*(I-1) - (I*(I-1))/2 + I
C
         II4 = NS2*(I-1)    - (I*(I-1))/2         + I 
         II5 = NS2*(I-1)    - (I*(I-1))/2         + I + NS
         II6 = NS2*(I+NS-1) - ((I+NS)*(I+NS-1))/2 + I + NS
         ZZZ = 5.0D0 * PATMOS * WT(I) / ( 3.0D0 * RU )
C-----------------------------------------------------------------------
C        FAC = 2 * A_{ii} / {\cal D}_{ii}
C-----------------------------------------------------------------------
         do nn = 1, np
            FAC = ZZZ / ( TEMPER(nn) * ETA(nn,I) )
            G(nn,II4) = FAC * X(nn,I) * X(nn,I) 
     &                   * ( 1.0D0 + 1.0D1 * CXI(nn,I) / 3.0D0 )
            G(nn,II5) = - FAC * X(nn,I) * X(nn,I) 
     &                   * 2.0D0 * CXI(nn,I) 
            G(nn,II6) = FAC * X(nn,I) * X(nn,I) 
     &                   * 1.2D0 * CXI(nn,I)
     &                   + X(nn,I) * X(nn,I) * CINT(nn,I) * BIN(nn,IND)
         enddo
      ENDDO
C-----------------------------------------------------------------------
      FFF = 2.0D1 / 3.0D0
      DO I = 1, NS
         II4 = NS2*(I-1)    - (I*(I-1))/2         + I 
         II5 = NS2*(I-1)    - (I*(I-1))/2         + I + NS
         II6 = NS2*(I+NS-1) - ((I+NS)*(I+NS-1))/2 + I + NS
C........
         IP1 = I + 1
         DO J = IP1, NS
            IND = NS*(I-1) - (I*(I-1))/2 + J
C
            JJ4 = NS2*(J-1)    - (J*(J-1))/2         + J 
            JJ5 = NS2*(J-1)    - (J*(J-1))/2         + J + NS
            JJ6 = NS2*(J+NS-1) - ((J+NS)*(J+NS-1))/2 + J + NS
C...........
            IJ4 = NS2*(I-1)    - (I*(I-1))/2         + J 
            IJ5 = NS2*(I-1)    - (I*(I-1))/2         + J + NS
            IJ6 = NS2*(I+NS-1) - ((I+NS)*(I+NS-1))/2 + J + NS
C...........
            JI5 = NS2*(J-1) - (J*(J-1))/2 + I + NS
C...........
            WWI = WT(I) / ( WT(I) + WT(J) )
            WWJ = WT(J) / ( WT(I) + WT(J) )
            WIJ = WT(I) / WT(J)
            WJI = WT(J) / WT(I)
C..............
            FITA0 = FITA(1,I,J)
            FITA1 = FITA(2,I,J)
            FITA2 = FITA(3,I,J)
            FITA3 = FITA(4,I,J)
            FITA4 = FITA(5,I,J)
            FITA5 = FITA(6,I,J)
            FITA6 = FITA(7,I,J)
C..............
            FITB0 = FITB(1,I,J)
            FITB1 = FITB(2,I,J)
            FITB2 = FITB(3,I,J)
            FITB3 = FITB(4,I,J)
            FITB4 = FITB(5,I,J)
            FITB5 = FITB(6,I,J)
            FITB6 = FITB(7,I,J)
            do nn = 1, np
               AIJ(nn) = FITA0          + FITA1*DLT1(nn) 
     &                 + FITA2*DLT2(nn) + FITA3*DLT3(nn)
     &                 + FITA4*DLT4(nn) + FITA5*DLT5(nn)
     &                 + FITA6*DLT6(nn) 
               BIJ(nn) = FITB0          + FITB1*DLT1(nn) 
     &                 + FITB2*DLT2(nn) + FITB3*DLT3(nn)
     &                 + FITB4*DLT4(nn) + FITB5*DLT5(nn)
     &                 + FITB6*DLT6(nn) 
            enddo
C..............
            do nn = 1, np
               AAA = X(nn,I) * X(nn,J) * BIN(nn,IND)
C ----- 10,10
               BB1 = FFF * AIJ(nn) * (CXI(nn,I)+CXI(nn,J))
               BBB = WWI * WWJ * (
     &               13.75D0 - 3.0D0*BIJ(nn) 
     &                - 4.0D0*AIJ(nn) - BB1 )
               BBJ = ( 7.5D0 * WWI * WWI + 
     &               (6.25D0 - 3.0D0*BIJ(nn)) * WWJ * WWJ 
     &                + 4.0D0*AIJ(nn) * WWI * WWJ
     &                + WWI * WWJ * BB1 )
               BBI = ( 7.5D0 * WWJ * WWJ + 
     &               (6.25D0 - 3.0D0*BIJ(nn)) * WWI * WWI 
     &                + 4.0D0*AIJ(nn) * WWI * WWJ
     &                + WWI * WWJ * BB1 )
               G(nn,IJ4) = - BBB * AAA
               G(nn,II4) = G(nn,II4) + BBJ * AAA
               G(nn,JJ4) = G(nn,JJ4) + BBI * AAA
C ----- 10,01 and 01,10
               CCI = WWI * 4.0D0 * AIJ(nn) * CXI(nn,I)
               CCJ = WWJ * 4.0D0 * AIJ(nn) * CXI(nn,J)
               G(nn,IJ5) = - CCJ * AAA
               G(nn,JI5) = - CCI * AAA
               G(nn,II5) = G(nn,II5) - CCI * AAA
               G(nn,JJ5) = G(nn,JJ5) - CCJ * AAA
C ----- 01,01
               DDJ = AAA * CINT(nn,I) 
     &            +  AAA * 2.4D0 * AIJ(nn) * WIJ * CXI(nn,I)
               DDI = AAA * CINT(nn,J) 
     &            +  AAA * 2.4D0 * AIJ(nn) * WJI * CXI(nn,J)
               G(nn,IJ6) = 0.0D0
               G(nn,II6) = G(nn,II6) + DDJ
               G(nn,JJ6) = G(nn,JJ6) + DDI
            enddo
         ENDDO
      ENDDO
      DO I = NS+1, NS2
         II6 = NS2*(I-1) - (I*(I-1))/2 + I
         IF ( LIN(I-NS) .EQ. 0) THEN
            do nn = 1, np
               G(nn,II6) = 1.0D0
            enddo
         ENDIF
      ENDDO
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGFEMAA ( NP, NS, TEMPER, X, WT, RU, PATMOS,
     &                     BIN, AIJ, BIJ, CINT, ETA, CXI,
     &                     LIN, G, BETA )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      DIMENSION G(NP,*), BETA(NP,*)
      DIMENSION BIN(NP,*), AIJ(NS,NS), BIJ(NS,NS), 
     &          CINT(NP,NS), ETA(NP,NS), X(NP,NS), 
     &          WT(NS), CXI(NP,NS), LIN(NS), TEMPER(NP)
C-----------------------------------------------------------------------
C         Form the RHS 
C-----------------------------------------------------------------------
      DO I = 1, NS
         do nn = 1, np
            BETA(nn,I)    = 2.5D0 * X(nn,I)
            BETA(nn,I+NS) = CINT(nn,I) * X(nn,I)
         enddo
      ENDDO
C-----------------------------------------------------------------------
C         Form the transport linear system matrix \Lambda
C
C     Note: the system matrix is evaluated at atmospheric pressure
C     A_ij, B_ij, and C_ij are node independent
C-----------------------------------------------------------------------
      NS2 = 2 * NS
C ----- Initialize the diagonals
      DO I = 1, NS
         IND = NS*(I-1) - (I*(I-1))/2 + I
C
         II4 = NS2*(I-1)    - (I*(I-1))/2         + I 
         II5 = NS2*(I-1)    - (I*(I-1))/2         + I + NS
         II6 = NS2*(I+NS-1) - ((I+NS)*(I+NS-1))/2 + I + NS
         ZZZ = 5.0D0 * PATMOS * WT(I) / ( 3.0D0 * RU )
C-----------------------------------------------------------------------
C        FAC = 2 * A_{ii} / {\cal D}_{ii}
C-----------------------------------------------------------------------
         do nn = 1, np
            FAC = ZZZ / ( TEMPER(nn) * ETA(nn,I) )
            G(nn,II4) = FAC * X(nn,I) * X(nn,I) 
     &                   * ( 1.0D0 + 1.0D1 * CXI(nn,I) / 3.0D0 )
            G(nn,II5) = - FAC * X(nn,I) * X(nn,I) 
     &                   * 2.0D0 * CXI(nn,I) 
            G(nn,II6) = FAC * X(nn,I) * X(nn,I) 
     &                   * 1.2D0 * CXI(nn,I)
     &                   + X(nn,I) * X(nn,I) * CINT(nn,I) * BIN(nn,IND)
         enddo
      ENDDO
C-----------------------------------------------------------------------
      FFF = 2.0D1 / 3.0D0
      DO I = 1, NS
         II4 = NS2*(I-1)    - (I*(I-1))/2         + I 
         II5 = NS2*(I-1)    - (I*(I-1))/2         + I + NS
         II6 = NS2*(I+NS-1) - ((I+NS)*(I+NS-1))/2 + I + NS
C........
         IP1 = I + 1
         DO J = IP1, NS
            IND = NS*(I-1) - (I*(I-1))/2 + J
C
            JJ4 = NS2*(J-1)    - (J*(J-1))/2         + J 
            JJ5 = NS2*(J-1)    - (J*(J-1))/2         + J + NS
            JJ6 = NS2*(J+NS-1) - ((J+NS)*(J+NS-1))/2 + J + NS
C...........
            IJ4 = NS2*(I-1)    - (I*(I-1))/2         + J 
            IJ5 = NS2*(I-1)    - (I*(I-1))/2         + J + NS
            IJ6 = NS2*(I+NS-1) - ((I+NS)*(I+NS-1))/2 + J + NS
C...........
            JI5 = NS2*(J-1) - (J*(J-1))/2 + I + NS
C...........
            WWI = WT(I) / ( WT(I) + WT(J) )
            WWJ = WT(J) / ( WT(I) + WT(J) )
            WIJ = WT(I) / WT(J)
            WJI = WT(J) / WT(I)
            ZZ1 = FFF * AIJ(I,J) 
            ZZ2 = 13.75D0 - 3.0D0*BIJ(I,J) - 4.0D0*AIJ(I,J)
            ZZ3 = 7.5D0 * WWI * WWI 
     &            + ( 6.25D0 - 3.0D0*BIJ(I,J) ) * WWJ * WWJ 
     &            + 4.0D0*AIJ(I,J) * WWI * WWJ
            ZZ4 = 7.5D0 * WWJ * WWJ 
     &            + ( 6.25D0 - 3.0D0*BIJ(I,J) ) * WWI * WWI 
     &            + 4.0D0*AIJ(I,J) * WWI * WWJ
            ZZ5 = WWI * 4.0D0 * AIJ(I,J) 
            ZZ6 = WWJ * 4.0D0 * AIJ(I,J) 
            ZZ7 = 2.4D0 * AIJ(I,J) * WIJ 
            ZZ8 = 2.4D0 * AIJ(I,J) * WJI 
            do nn = 1, np
               AAA = X(nn,I) * X(nn,J) * BIN(nn,IND)
C ----- 10,10
               BB1 = ZZ1 * (CXI(nn,I)+CXI(nn,J))
               BBB = WWI * WWJ * ( ZZ2 - BB1 )
               BBJ = ZZ3 + WWI * WWJ * BB1 
               BBI = ZZ4 + WWI * WWJ * BB1 
               G(nn,IJ4) = - BBB * AAA
               G(nn,II4) = G(nn,II4) + BBJ * AAA
               G(nn,JJ4) = G(nn,JJ4) + BBI * AAA
C ----- 10,01 and 01,10
               CCI = ZZ5 * CXI(nn,I)
               CCJ = ZZ6 * CXI(nn,J)
               G(nn,IJ5) = - CCJ * AAA
               G(nn,JI5) = - CCI * AAA
               G(nn,II5) = G(nn,II5) - CCI * AAA
               G(nn,JJ5) = G(nn,JJ5) - CCJ * AAA
C ----- 01,01
               DDJ = AAA * CINT(nn,I) 
     &            +  AAA * ZZ7 * CXI(nn,I)
               DDI = AAA * CINT(nn,J)
     &            +  AAA * ZZ8 * CXI(nn,J)
               G(nn,IJ6) = 0.0D0
               G(nn,II6) = G(nn,II6) + DDJ
               G(nn,JJ6) = G(nn,JJ6) + DDI
            enddo
         ENDDO
      ENDDO
      DO I = NS+1, NS2
         II6 = NS2*(I-1) - (I*(I-1))/2 + I
         IF ( LIN(I-NS) .EQ. 0) THEN
            do nn = 1, np
               G(nn,II6) = 1.0D0
            enddo
         ENDIF
      ENDDO
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGFEML00 ( NP, NS, X, BIN, G )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      DIMENSION G(NP,*), BIN(NP,*), X(NP,NS)
C-----------------------------------------------------------------------
C         Form the transport linear system matrix L_[00]
C
C     Note: the system matrix is evaluated at atmospheric pressure.
C-----------------------------------------------------------------------
C ----- Initialize the diagonals
      DO I = 1, NS
         III = NS*(I-1) - (I*(I-1))/2 + I
         do nn = 1, np
            G(nn,III) = 0.0D0
         enddo
      ENDDO
C-----------------------------------------------------------------------
      DO I = 1, NS
         IP1 = I + 1
         III = NS*(I-1) - (I*(I-1))/2 + I
         DO J = IP1, NS
            IND = NS*(I-1) - (I*(I-1))/2 + J
            JJJ = NS*(J-1) - (J*(J-1))/2 + J
            do nn = 1, np
               AAA = X(nn,I) * X(nn,J) * BIN(nn,IND)
               G(nn,IND) = - AAA
               G(nn,III) = G(nn,III) + AAA
               G(nn,JJJ) = G(nn,JJJ) + AAA
            enddo
         ENDDO
      ENDDO
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGFEMLE ( NP, NS, 
     &                     TEMPER, DLT1, DLT2, DLT3, DLT4, DLT5, DLT6, 
     &                     X, WT, RU, PATMOS,
     &                     BIN, AIJ, BIJ, CIJ, FITA, FITB, FITC,
     &                     CINT, ETA, CXI, LIN, G, BETA )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      DIMENSION G(NP,*), BETA(NP,*)
      DIMENSION BIN(NP,*), AIJ(NP), BIJ(NP), 
     &          CIJ(NP), CINT(NP,NS), ETA(NP,NS), X(NP,NS), 
     &          WT(NS), CXI(NP,NS), LIN(NS)
      DIMENSION TEMPER(NP), DLT1(NP), DLT2(NP), DLT3(NP), DLT4(NP), 
     &          DLT5(NP), DLT6(NP)
      DIMENSION FITA(7,NS,NS), FITB(7,NS,NS), FITC(7,NS,NS)
C-----------------------------------------------------------------------
C         Form the RHS 
C-----------------------------------------------------------------------
      DO I = 1, NS
         do nn = 1, np
            BETA(nn,I)      = 0.0D0
            BETA(nn,I+NS)   = ( 2.5D0 + CINT(nn,I) ) * X(nn,I)
         enddo
      ENDDO
C-----------------------------------------------------------------------
C         Form the transport linear system matrix L_[e]
C
C     Note: the system matrix is evaluated at atmospheric pressure
C-----------------------------------------------------------------------
      NS2 = 2 * NS
C ----- Initialize the diagonals
      DO I = 1, NS
         IND = NS*(I-1) - (I*(I-1))/2 + I
C
         II1 = NS2*(I-1)     - (I*(I-1))/2           + I
         II2 = NS2*(I-1)     - (I*(I-1))/2           + I + NS
         II4 = NS2*(I+NS-1)  - ((I+NS)*(I+NS-1))/2   + I + NS
         ZZZ = 5.0D0 * PATMOS * WT(I) / ( 3.0D0 * RU )
         do nn = 1, np
            G(nn,II1) = 0.0D0
            G(nn,II2) = 0.0D0
C-----------------------------------------------------------------------
C        FAC = 2 * A_{ii} / {\cal D}_{ii}
C-----------------------------------------------------------------------
            FAC = ZZZ / ( TEMPER(nn) * ETA(nn,I) )
            G(nn,II4) = FAC * X(nn,I) * X(nn,I) 
     &                   * ( 1.0D0 + 8.0D0 * CXI(nn,I) / 1.5D1 )
     &                + X(nn,I) * X(nn,I) * CINT(nn,I) * BIN(nn,IND)
         enddo
      ENDDO
C-----------------------------------------------------------------------
      DDD = 4.0D0 / 3.0D0
      EEE = 4.0D0 / 1.5D1
      FFF = 2.0D1 / 3.0D0
      DO I = 1, NS
         II1 = NS2*(I-1) - (I*(I-1))/2 + I
         II2 = NS2*(I-1) - (I*(I-1))/2 + I + NS
         II4 = NS2*(I+NS-1) - ((I+NS)*(I+NS-1))/2 + I + NS
C........
         IP1 = I + 1
         DO J = IP1, NS
            IND = NS*(I-1) - (I*(I-1))/2 + J
C
            JJ1 = NS2*(J-1) - (J*(J-1))/2 + J
            JJ2 = NS2*(J-1) - (J*(J-1))/2 + J + NS
            JJ4 = NS2*(J+NS-1) - ((J+NS)*(J+NS-1))/2 + J + NS
C...........
            IJ1 = NS2*(I-1) - (I*(I-1))/2 + J
            IJ2 = NS2*(I-1) - (I*(I-1))/2 + J + NS
            IJ4 = NS2*(I+NS-1) - ((I+NS)*(I+NS-1))/2 + J + NS
C...........
            JI2 = NS2*(J-1) - (J*(J-1))/2 + I + NS
C...........
            WWI  = WT(I) / ( WT(I) + WT(J) )
            WWJ  = WT(J) / ( WT(I) + WT(J) )
            WIOJ = 3.0D0 * WT(I) / WT(J) - 2.0D0
            WJOI = 3.0D0 * WT(J) / WT(I) - 2.0D0
C..............
            FITA0 = FITA(1,I,J)
            FITA1 = FITA(2,I,J)
            FITA2 = FITA(3,I,J)
            FITA3 = FITA(4,I,J)
            FITA4 = FITA(5,I,J)
            FITA5 = FITA(6,I,J)
            FITA6 = FITA(7,I,J)
C..............
            FITB0 = FITB(1,I,J)
            FITB1 = FITB(2,I,J)
            FITB2 = FITB(3,I,J)
            FITB3 = FITB(4,I,J)
            FITB4 = FITB(5,I,J)
            FITB5 = FITB(6,I,J)
            FITB6 = FITB(7,I,J)
C..............
            FITC0 = FITC(1,I,J)
            FITC1 = FITC(2,I,J)
            FITC2 = FITC(3,I,J)
            FITC3 = FITC(4,I,J)
            FITC4 = FITC(5,I,J)
            FITC5 = FITC(6,I,J)
            FITC6 = FITC(7,I,J)
            do nn = 1, np
               AIJ(nn) = FITA0          + FITA1*DLT1(nn) 
     &                 + FITA2*DLT2(nn) + FITA3*DLT3(nn)
     &                 + FITA4*DLT4(nn) + FITA5*DLT5(nn)
     &                 + FITA6*DLT6(nn) 
               BIJ(nn) = FITB0          + FITB1*DLT1(nn) 
     &                 + FITB2*DLT2(nn) + FITB3*DLT3(nn)
     &                 + FITB4*DLT4(nn) + FITB5*DLT5(nn)
     &                 + FITB6*DLT6(nn) 
               CIJ(nn) = FITC0          + FITC1*DLT1(nn) 
     &                 + FITC2*DLT2(nn) + FITC3*DLT3(nn)
     &                 + FITC4*DLT4(nn) + FITC5*DLT5(nn)
     &                 + FITC6*DLT6(nn) 
            enddo
C..............
            do nn = 1, np
C ----- 00,00
               AAA = X(nn,I) * X(nn,J) * BIN(nn,IND)
               G(nn,IJ1) = - AAA
               G(nn,II1) = G(nn,II1) + AAA
               G(nn,JJ1) = G(nn,JJ1) + AAA
C ----- 00,e and e,00
               AAI = 0.5D0 * WWI * (6.0D0*CIJ(nn) - 5.0D0)
               AAJ = 0.5D0 * WWJ * (6.0D0*CIJ(nn) - 5.0D0)
               G(nn,IJ2) = AAI * AAA
               G(nn,JI2) = AAJ * AAA
               G(nn,II2) = G(nn,II2) - AAJ * AAA
               G(nn,JJ2) = G(nn,JJ2) - AAI * AAA
C ----- e,e
               TTI = DDD * AIJ(nn) * WIOJ * CXI(nn,I)
               TTJ = DDD * AIJ(nn) * WJOI * CXI(nn,J)
               BBB = WWI * WWJ * (
     &               13.75D0 - 3.0D0*BIJ(nn) 
     &                - 4.0D0*AIJ(nn) + TTI + TTJ )
               BBJ = 7.5D0 * WWI * WWI  
     &            + (6.25D0 - 3.0D0*BIJ(nn)) * WWJ * WWJ 
     &            + 4.0D0*AIJ(nn) * WWI * WWJ
     &            + EEE * AIJ(nn) * CXI(nn,I)
     &                          * WWI * WWJ * WIOJ * WIOJ
     &            + FFF * AIJ(nn) * CXI(nn,J) * WWI * WWJ
               BBI = 7.5D0 * WWJ * WWJ  
     &            + (6.25D0 - 3.0D0*BIJ(nn)) * WWI * WWI 
     &            + 4.0D0*AIJ(nn) * WWI * WWJ
     &            + EEE * AIJ(nn) * CXI(nn,J)
     &                          * WWI * WWJ * WJOI * WJOI
     &            + FFF * AIJ(nn) * CXI(nn,I) * WWI * WWJ
               DDJ = AAA * CINT(nn,I) 
               DDI = AAA * CINT(nn,J) 
               G(nn,IJ4) = - BBB * AAA
               G(nn,II4) = G(nn,II4) + BBJ * AAA + DDJ
               G(nn,JJ4) = G(nn,JJ4) + BBI * AAA + DDI
            enddo
         ENDDO
      ENDDO
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGFEMAE ( NP, NS, 
     &                     TEMPER, DLT1, DLT2, DLT3, DLT4, DLT5, DLT6, 
     &                     X, WT, RU, PATMOS,
     &                     BIN, AIJ, BIJ, FITA, FITB, 
     &                     CINT, ETA, CXI, LIN, G, BETA )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      DIMENSION G(NP,*), BETA(NP,*)
      DIMENSION BIN(NP,*), AIJ(NP), BIJ(NP), 
     &          CINT(NP,NS), ETA(NP,NS), X(NP,NS), 
     &          WT(NS), CXI(NP,NS), LIN(NS)
      DIMENSION TEMPER(NP), DLT1(NP), DLT2(NP), DLT3(NP), DLT4(NP), 
     &          DLT5(NP), DLT6(NP)
      DIMENSION FITA(7,NS,NS), FITB(7,NS,NS)
C-----------------------------------------------------------------------
C         Form the RHS 
C-----------------------------------------------------------------------
      DO I = 1, NS
         do nn = 1, np
            BETA(nn,I) = ( 2.5D0 + CINT(nn,I) ) * X(nn,I)
         enddo
      ENDDO
C-----------------------------------------------------------------------
C         Form the transport linear system matrix \Lambda_[e]
C
C     Note: the system matrix is evaluated at atmospheric pressure
C-----------------------------------------------------------------------
C ----- Initialize the diagonals
      DO I = 1, NS
         II4 = NS*(I-1)  - (I*(I-1))/2   + I 
         ZZZ = 5.0D0 * PATMOS * WT(I) / ( 3.0D0 * RU )
C-----------------------------------------------------------------------
C        FAC = 2 * A_{ii} / {\cal D}_{ii}
C-----------------------------------------------------------------------
         do nn = 1, np
            FAC = ZZZ / ( TEMPER(nn) * ETA(nn,I) )
            G(nn,II4) = FAC * X(nn,I) * X(nn,I) 
     &                   * ( 1.0D0 + 8.0D0 * CXI(nn,I) / 1.5D1 )
     &                + X(nn,I) * X(nn,I) * CINT(nn,I) * BIN(nn,II4)
         enddo
      ENDDO
C-----------------------------------------------------------------------
      DDD = 4.0D0 / 3.0D0
      EEE = 4.0D0 / 1.5D1
      FFF = 2.0D1 / 3.0D0
      DO I = 1, NS
         II4 = NS*(I-1) - (I*(I-1))/2 + I
C........
         IP1 = I + 1
         DO J = IP1, NS
            JJ4 = NS*(J-1) - (J*(J-1))/2 + J
            IJ4 = NS*(I-1) - (I*(I-1))/2 + J
C...........
            WWI  = WT(I) / ( WT(I) + WT(J) )
            WWJ  = WT(J) / ( WT(I) + WT(J) )
            WIOJ = 3.0D0 * WT(I) / WT(J) - 2.0D0
            WJOI = 3.0D0 * WT(J) / WT(I) - 2.0D0
C..............
            FITA0 = FITA(1,I,J)
            FITA1 = FITA(2,I,J)
            FITA2 = FITA(3,I,J)
            FITA3 = FITA(4,I,J)
            FITA4 = FITA(5,I,J)
            FITA5 = FITA(6,I,J)
            FITA6 = FITA(7,I,J)
C..............
            FITB0 = FITB(1,I,J)
            FITB1 = FITB(2,I,J)
            FITB2 = FITB(3,I,J)
            FITB3 = FITB(4,I,J)
            FITB4 = FITB(5,I,J)
            FITB5 = FITB(6,I,J)
            FITB6 = FITB(7,I,J)
            do nn = 1, np
               AIJ(nn) = FITA0          + FITA1*DLT1(nn) 
     &                 + FITA2*DLT2(nn) + FITA3*DLT3(nn)
     &                 + FITA4*DLT4(nn) + FITA5*DLT5(nn)
     &                 + FITA6*DLT6(nn) 
               BIJ(nn) = FITB0          + FITB1*DLT1(nn) 
     &                 + FITB2*DLT2(nn) + FITB3*DLT3(nn)
     &                 + FITB4*DLT4(nn) + FITB5*DLT5(nn)
     &                 + FITB6*DLT6(nn) 
            enddo
C..............
            do nn = 1, np
               TTI = DDD * AIJ(nn) * WIOJ * CXI(nn,I)
               TTJ = DDD * AIJ(nn) * WJOI * CXI(nn,J)
               AAA = X(nn,I) * X(nn,J) * BIN(nn,IJ4)
               BBB = WWI * WWJ * (
     &               13.75D0 - 3.0D0*BIJ(nn) 
     &                - 4.0D0*AIJ(nn) + TTI + TTJ )
               BBJ = 7.5D0 * WWI * WWI  
     &            + (6.25D0 - 3.0D0*BIJ(nn)) * WWJ * WWJ 
     &            + 4.0D0*AIJ(nn) * WWI * WWJ
     &            + EEE * AIJ(nn) * CXI(nn,I)
     &                          * WWI * WWJ * WIOJ * WIOJ
     &            + FFF * AIJ(nn) * CXI(nn,J) * WWI * WWJ
               BBI = 7.5D0 * WWJ * WWJ  
     &            + (6.25D0 - 3.0D0*BIJ(nn)) * WWI * WWI 
     &            + 4.0D0*AIJ(nn) * WWI * WWJ
     &            + EEE * AIJ(nn) * CXI(nn,J)
     &                          * WWI * WWJ * WJOI * WJOI
     &            + FFF * AIJ(nn) * CXI(nn,I) * WWI * WWJ
               DDJ = AAA * CINT(nn,I) 
               DDI = AAA * CINT(nn,J) 
C...........
               G(nn,IJ4) = - BBB * AAA
               G(nn,II4) = G(nn,II4) + BBJ * AAA + DDJ
               G(nn,JJ4) = G(nn,JJ4) + BBI * AAA + DDI
            enddo
         ENDDO
      ENDDO
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE EGFEMAAE ( NP, NS, TEMPER, X, WT, RU, PATMOS,
     &                      BIN, AIJ, BIJ, CINT, ETA, CXI,
     &                      LIN, G, BETA )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
      DIMENSION G(NP,*), BETA(NP,*)
      DIMENSION BIN(NP,*), AIJ(NS,NS), BIJ(NS,NS), 
     &          CINT(NP,NS), ETA(NP,NS), X(NP,NS), 
     &          WT(NS), CXI(NP,NS), LIN(NS), TEMPER(NP)
C-----------------------------------------------------------------------
C         Form the RHS 
C-----------------------------------------------------------------------
      DO I = 1, NS
         do nn = 1, np
            BETA(nn,I) = ( 2.5D0 + CINT(nn,I) ) * X(nn,I)
         enddo
      ENDDO
C-----------------------------------------------------------------------
C         Form the transport linear system matrix \Lambda_[e]
C
C     Note: the system matrix is evaluated at atmospheric pressure
C     A_ij, B_ij, and C_ij are node independent
C-----------------------------------------------------------------------
C ----- Initialize the diagonals
      DO I = 1, NS
         II4 = NS*(I-1)  - (I*(I-1))/2   + I 
         ZZZ = 5.0D0 * PATMOS * WT(I) / ( 3.0D0 * RU )
C-----------------------------------------------------------------------
C        FAC = 2 * A_{ii} / {\cal D}_{ii}
C-----------------------------------------------------------------------
         do nn = 1, np
            FAC = ZZZ / ( TEMPER(nn) * ETA(nn,I) )
            G(nn,II4) = FAC * X(nn,I) * X(nn,I) 
     &                   * ( 1.0D0 + 8.0D0 * CXI(nn,I) / 1.5D1 )
     &                + X(nn,I) * X(nn,I) * CINT(nn,I) * BIN(nn,II4)
         enddo
      ENDDO
C-----------------------------------------------------------------------
      DDD = 4.0D0 / 3.0D0
      EEE = 4.0D0 / 1.5D1
      FFF = 2.0D1 / 3.0D0
      DO I = 1, NS
         II4 = NS*(I-1) - (I*(I-1))/2 + I
C........
         IP1 = I + 1
         DO J = IP1, NS
            JJ4 = NS*(J-1) - (J*(J-1))/2 + J
            IJ4 = NS*(I-1) - (I*(I-1))/2 + J
C...........
            WWI  = WT(I) / ( WT(I) + WT(J) )
            WWJ  = WT(J) / ( WT(I) + WT(J) )
            WIOJ = 3.0D0 * WT(I) / WT(J) - 2.0D0
            WJOI = 3.0D0 * WT(J) / WT(I) - 2.0D0
            ZZ1  = DDD * AIJ(I,J) * WIOJ 
            ZZ2  = DDD * AIJ(I,J) * WJOI 
            ZZ3  = 13.75D0 - 3.0D0*BIJ(I,J) - 4.0D0*AIJ(I,J) 
            ZZ4  = 7.5D0 * WWI * WWI  
     &            + (6.25D0 - 3.0D0*BIJ(I,J)) * WWJ * WWJ 
     &            + 4.0D0*AIJ(I,J) * WWI * WWJ
            ZZ5  = 7.5D0 * WWJ * WWJ  
     &            + (6.25D0 - 3.0D0*BIJ(I,J)) * WWI * WWI 
     &            + 4.0D0*AIJ(I,J) * WWI * WWJ
            ZZ6  = EEE * AIJ(I,J) * WWI * WWJ * WIOJ * WIOJ
            ZZ7  = EEE * AIJ(I,J) * WWI * WWJ * WJOI * WJOI
            ZZ8  = FFF * AIJ(I,J) * WWI * WWJ 
            do nn = 1, np
               TTI = ZZ1 * CXI(nn,I)
               TTJ = ZZ2 * CXI(nn,J)
               AAA = X(nn,I) * X(nn,J) * BIN(nn,IJ4)
               BBB = WWI * WWJ * ( ZZ3 + TTI + TTJ )
               BBJ = ZZ4 + ZZ6 * CXI(nn,I) + ZZ8 * CXI(nn,J)
               BBI = ZZ5 + ZZ7 * CXI(nn,J) + ZZ8 * CXI(nn,I)
               DDJ = AAA * CINT(nn,I) 
               DDI = AAA * CINT(nn,J) 
C...........
               G(nn,IJ4) = - BBB * AAA
               G(nn,II4) = G(nn,II4) + BBJ * AAA + DDJ
               G(nn,JJ4) = G(nn,JJ4) + BBI * AAA + DDI
            enddo
         ENDDO
      ENDDO
C-----------------------------------------------------------------------
      RETURN
      END

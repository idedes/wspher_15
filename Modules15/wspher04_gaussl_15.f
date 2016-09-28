C FILE NAME = wspher04_gaussl_15.f ! Keep this symbol:    $ident@string$
C
C=======================================================================
C=======================================================================
C============      G A U S S    R O U T I N E   S E T   ================
C=======================================================================
C=======================================================================
C
      SUBROUTINE GENBAS(NDIM_N,NDIM_L,NDGAUS,NGAUSS,ANODES,AWEIGH,
     *                  NSHELL,BAUXIL,CAUXIL,EPSLAG,QRNODE,QRWEIG,
     *                         QNORNL,QDFLAG,QFLAGN,QPOLYN,QDPOLY,
     *                         IALPHA,QPOLYN_DENSIT,QDPOLY_DENSIT)
C
      IMPLICIT  REAL*16 (A-H,O-Z)
C
      DIMENSION
     *          QPOLYN(1:NDGAUS,0:NDIM_N,0:NDIM_L), 
     *          QDPOLY(1:NDGAUS,0:NDIM_N,0:NDIM_L)  
      DIMENSION
     *          QFLAGN(0:NDIM_N,0:NDIM_L),
     *          QDFLAG(0:NDIM_N,0:NDIM_L)
      DIMENSION
     *          QRNODE(1:NDGAUS,0:NDIM_L),
     *          QRWEIG(1:NDGAUS,0:NDIM_L)
      DIMENSION 
     *          QNORNL(0:NDIM_N,0:NDIM_L)
      DIMENSION
     *          ANODES(1:NDGAUS),
     *          AWEIGH(1:NDGAUS)
      DIMENSION
     *          BAUXIL(1:NDGAUS),
     *          CAUXIL(1:NDGAUS)
      DIMENSION
     *          QFLEPS(0:NDIM_N,0:NDIM_L),
     *          QDFEPS(0:NDIM_N,0:NDIM_L) 
      DIMENSION
     *          QPOLYN_DENSIT(1:NDGAUS,0:NDIM_N,0:NDIM_L,0:NDIM_L), 
     *          QDPOLY_DENSIT(1:NDGAUS,0:NDIM_N,0:NDIM_L,0:NDIM_L)
C
      COMMON
     *         /PRFILE/ ISCREN,LOGFIL,LOGCHP,LOGCHN,LOGHIP,LOGHIN,
     *                  LOGPRO,LOGNEU,LOGWRI,LSCREN,IRESUL,ICONVE,
     *                  NOUTPT,LOGAUX,LOGENE,LOGRAD,LOGGRD,LOGMSH,
     *                                                     LOGBIS                                                 
C_______________________________________________________________________
C
      DATA
     *        EPSORT / 1.0Q-30 /     
C
C=======================================================================
C             Harmonic-Oscillator Basis-Generating Routine
C=======================================================================
C
C     Subroutine GENBAS:  Provides the Gauss-Laguerre integration 
C     points and weights with the indices alpha=LQNUMB-1/2 rather
C     than usual <alpha=LQNUMB+1/2> for the kinetic energy matrix
C     elements integration.  This part of the result  is provided
C     through call to subroutine LAGAUS for LQNUMB in [0,NDIM_L].
C
C     This subroutine  also provides the numerical values  of the 
C     generalised  Laguerre polynomials  and corresponding  first 
C     order derivatives. This part comes through  <<call LAGUER>>
C
C     Finally it calculates the norms of the Laguerre polynomials
C     and verifies the functioning  of the integration by testing
C     the Laguerre orthogonality relations.
C
C     Meaning of some parameters:
C
C     NSHELL - The maximum value of the LQNUMB index  of Laguerre 
C              functions to be  calculated.  At the end, we store 
C              the polynomials with indices LQNUMB in [0,NSHELL]
C
C              The NQNUMB indices run between 0 and NSHELL/2 
C
C=======================================================================
C
      IF (NGAUSS.GT.NDGAUS) THEN
          WRITE(LSCREN,'(/,''NGAUSS='',I2,'' exceeds NDGAUS='',I2,/)')
     *                       NGAUSS,                 NDGAUS
          STOP 'Stop! NGAUSS > NDGAUS in GENBAS'
      END IF
C
C=======================================================================
C     Proceeding to define the Gauss quadrature nodes and weights
C     for the series of the L-values; L are integer, ALPHAL=L+0.5
C=======================================================================
C
      MAXIML=NSHELL
C
      DO LQNUMB=0,MAXIML      
C
C        Attention: Below we define the weight factor 
C                   as proportional to   z^{\ell-1/2}
C                   rather then the standard exponent z^{\ell+1/2}
C
         IF (IALPHA.EQ.0) ALPHAL=QFLOAT(LQNUMB)-0.5Q0
         IF (IALPHA.EQ.1) ALPHAL=QFLOAT(LQNUMB)+0.5Q0
C
         DO NPOINT=1,NGAUSS
            BAUXIL(NPOINT)= ALPHAL+QFLOAT(2*NPOINT-1)
            CAUXIL(NPOINT)=(ALPHAL+QFLOAT(  NPOINT-1))*QFLOAT(NPOINT-1)
         END DO
C
         IF (LOGWRI.GT.4) THEN
             WRITE(LOGFIL,'(9X,''Entering LAGAUS from GENBAS, '',
     *                         ''Orbital momentum L='',i2, 1X,
     *                         ''Weight function z^{\ell-1/2}'')') 
     *                                            LQNUMB
         END IF
C
         CALL LAGAUS(NGAUSS,ANODES,AWEIGH,ALPHAL,BAUXIL,CAUXIL,
     *                                           EPSLAG,NDGAUS)
         DO IPOINT=1,NGAUSS
C
            QRNODE(IPOINT,LQNUMB)=ANODES(IPOINT)
            QRWEIG(IPOINT,LQNUMB)=AWEIGH(IPOINT)
C
         END DO
C
C=======================================================================
C @@@ IRENE
C     ATTENTION: Le test suivant demande a reduire NGAUSS, dans le
C                fichier de donnees, suffisamment pour ne pas depasser
C                NDIM_N (qui est ici de 20). Dans le cas contraire,
C                le test est evite via le go to 1.
C
C=======================================================================
C @@@ IRENE VERIFY THIS NONSENSE
C
         IF (LOGWRI.GT.1) THEN
             WRITE(LOGFIL,'(12x,''LQNUMB='',I3,''  NGAUSS='',I3,2X,
     *                         ''NDIM_N='',I3,'' check the nonsense'')')
     *                           LQNUMB,          NGAUSS,NDIM_N
         END IF
C
         IF (NGAUSS.GT.NDIM_N) GO TO 1
C
C=======================================================================
C
         KIND_L=1
C
         NQNUMB=NGAUSS
         LORDER=LQNUMB
C
C=======================================================================
C
         IF (IALPHA.EQ.1) THEN
C
             DO IPOINT=1,NQNUMB
C
                XARGUM=QRNODE(IPOINT,LQNUMB)
C
                CALL LAGUER(XARGUM,QFLAGN,QDFLAG,NQNUMB,LORDER,KIND_L,
     *                                                  NDIM_L,NDIM_N)
             END DO
C
         END IF
C
C=======================================================================
C
         IF (IALPHA.NE.1) THEN
C
             DO IPOINT=1,NQNUMB
C
                XARGUM=QRNODE(IPOINT,LQNUMB)
C
                CALL LAGUER(XARGUM,QFLAGN,QDFLAG,NQNUMB,LORDER,KIND_L,
     *                                                  NDIM_L,NDIM_N)    
             END DO
C
         END IF
C
C=======================================================================
C
   1  CONTINUE
C
C=======================================================================
C
      END DO
C
C=======================================================================
C     Calculating the normalisation constants for the HO functions
C=======================================================================
C
      NORDER=NSHELL/2
      LORDER=NSHELL
      KIND_L=1
C
      CALL NORMNL(QNORNL,NORDER,LORDER,NDIM_N,NDIM_L,KIND_L)
C
C=======================================================================
C     Below we calculate the values of the generalised Laguerre
C     polynomials and their first derivatives using recurrences;
C     In other words: we generate the spatial representation of
C     the harmonic oscillator basis in terms of n & l q-numbers
C=======================================================================
C       
       DO N_MAIN=NSHELL,0,-1
C
          DO LQNUMB=N_MAIN,0,-2
             NQNUMB=(N_MAIN-LQNUMB)/2
C
             DO IPOINT=1,NGAUSS
C
                RARGUM=QRNODE(IPOINT,LQNUMB)
C
                CALL LAGUER(RARGUM,QFLAGN,QDFLAG,NQNUMB,LQNUMB,KIND_L,
     *                                                  NDIM_L,NDIM_N)
C
                QPOLYN(IPOINT,NQNUMB,LQNUMB)=QFLAGN(NQNUMB,LQNUMB)
     *                                      *QNORNL(NQNUMB,LQNUMB)
C
                QDPOLY(IPOINT,NQNUMB,LQNUMB)=QDFLAG(NQNUMB,LQNUMB)
     *                                      *QNORNL(NQNUMB,LQNUMB)
            END DO
C
         END DO
      END DO
C      
      DO LQNUMB=NSHELL,0,-1
         DO IPOINT=1,NGAUSS
C            
            RARGUM=QRNODE(IPOINT,LQNUMB)
C
            CALL LAGUER(RARGUM,QFLAGN,QDFLAG,NORDER,LORDER,KIND_L,
     *                                              NDIM_L,NDIM_N)
            DO NAUXIL=0,NORDER
               DO LAUXIL=0,LORDER
C                       
                  QPOLYN_DENSIT(IPOINT,NAUXIL,LAUXIL,LQNUMB)
     *            =
     *            QFLAGN(NAUXIL,LAUXIL)*QNORNL(NAUXIL,LAUXIL)
C     
                  QDPOLY_DENSIT(IPOINT,NAUXIL,LAUXIL,LQNUMB)
     *            =
     *            QDFLAG(NAUXIL,LAUXIL)*QNORNL(NAUXIL,LAUXIL)
C     
               END DO
            END DO
             
         END DO
      END DO
C
C=======================================================================
C     Below we are verifying correctness of the numerical result
C     by testing the precision  of the orthonormality relations 
C=======================================================================
C
      IF (LOGWRI.GT.5) THEN
          WRITE(LOGFIL,'(/,12X,''Testing precision of orthogonality'',
     *                                                             /)')
      END IF
C
      DO LQNUMB=NSHELL,0,-1
C
C        Here we are going to verify the orthonormality relations
C        at fixed  L <<==>> LQNUMB values for L=0 to NSHELL
C_______________________________________________________________________
C
         IF (LOGWRI.GT.5) THEN
             WRITE(LOGFIL,'(/,9X,''LGNUMB='',I3)') LQNUMB
             WRITE(LOGFIL,'(36X,''.123456789 123456789'',
     *                          '' 123456789 123'')')
         END IF
C_______________________________________________________________________
C
         ORTMAX=-999.0
C
         DO N1NUMB=0,(NSHELL-LQNUMB)/2
            DO N2NUMB=0,(NSHELL-LQNUMB)/2       
C  
               ORTOGO=0.0Q0
C       
               DO I_PNTS=1,NGAUSS
C
                  IF (IALPHA.EQ.0) THEN
C
C                     We are using the weight factor e^{-z} z^{l-1/2}
C                     and THEREFORE the polynomials will be multiplied 
C                     by an extra factor "z"...
C
                      ZNODES=QRNODE(I_PNTS,LQNUMB)
C
                      ORTOGO=ORTOGO+QPOLYN(I_PNTS,N1NUMB,LQNUMB)
     *                             *QPOLYN(I_PNTS,N2NUMB,LQNUMB)
     *                             *QRWEIG(I_PNTS,LQNUMB)
     *                             *ZNODES                     
                  END IF
C
                  IF (IALPHA.EQ.1) THEN
C
C                     ... whereas: 
C                         the usual weight function is e^{-z} z^{l+1/2}
C
                      ORTOGO=ORTOGO+QPOLYN(I_PNTS,N1NUMB,LQNUMB)
     *                             *QPOLYN(I_PNTS,N2NUMB,LQNUMB)
     *                             *QRWEIG(I_PNTS,LQNUMB)
                  END IF
C
               END DO
C
               IF (N1NUMB.NE.N2NUMB) THEN
                   IF (ORTMAX.LT.ORTOGO) ORTMAX=ORTOGO
               END IF
C
               IF (LOGWRI.GT.5) THEN
C
               WRITE(LOGFIL,'(9X,''N1='',I2,'' N2='',I2,3X,''ORTOGO='',
     *                                                       F40.33)')
     *                             N1NUMB,     N2NUMB,       ORTOGO
               END IF
C______________________________________________________________________
C
                                     TESVAL=ORTOGO
               IF (N1NUMB.EQ.N2NUMB) TESVAL=ORTOGO-1.0Q0
C
               IF (ABS(TESVAL).GT.EPSORT) THEN
C
                   WRITE(LSCREN,'(/,''You requested the quadruple '',
     *                         ''precision; found non-orthogonality '',
     *                         ''on the level of EPSORT='',e12.5)') 
     *                                                      EPSORT
                   WRITE(LSCREN,'(''with LQNUMB='',i2,'' NGAUSS='',i2,
     *                          /,''     N1NUMB='',i2,'' N2NUMB='',i2)')
     *                                   LQNUMB,         NGAUSS,
     *                                   N1NUMB,         N2NUMB
                   WRITE(LSCREN,'(''TESVAL='',E32.24)') TESVAL
C
                   STOP 'Stop! - from integration precision test'
C
               END IF
C______________________________________________________________________
C      
            END DO ! N2NUMB
         END DO ! N1NUMB
C
         IF (LOGWRI.GT.5 .AND. LQNUMB.LT.(NSHELL-1)) THEN
             WRITE(LOGFIL,'(13X,''LQNUMB='',I2,'' ORTMAX='',E12.5)')
     *                            LQNUMB,         ORTMAX
         END IF
C      
      END DO ! LGNUMB
C
      IF (LOGWRI.GT.5) THEN
          WRITE(LOGFIL,'(T24,''EPSLAG='',E12.5)') EPSLAG
      END IF
C
C=======================================================================
C
      CALL SET_FACTORIAL      
C
C=======================================================================
C
      IF (LOGWRI.GT.0) THEN
          WRITE(LOGFIL,'(9X,''Exiting  GENBAS'')')
      END IF
C
C=======================================================================
C      
      RETURN
      END      
C
C=======================================================================
C=======================================================================
C
      SUBROUTINE NORMNL(XNORNL,NORDER,LORDER,NDIM_N,NDIM_L,KIND_L)
C
C                Here QUADRUPLE precision; its impact begins for L>30
      IMPLICIT
     *          REAL*16 (A-H,O-Z)
C
      DIMENSION 
     *          XNORNL(0:NDIM_N,0:NDIM_L)
      DIMENSION 
     *          FACTOR(0:40)
C       
C=======================================================================
C     This subroutine  calculates the normalisation constants
C     for the radial wave-functions of the spherical harmonic
C     oscillator solutions. The stretching factor is ignored.
C-----------------------------------------------------------------------
CID:
C     The formula is (for KIND_L=1):
C
C     N_{nl} = SQRT[ ( 2^{n+l} * n! ) / (( 2n + 2l + 1 )!! * sqrt(pi)) ]
CID.
C=======================================================================
C
      FACTOR(0)=1.0Q0
      FACTOR(1)=1.0Q0
C
      DO L=2,40
         FACTOR(L)=QSQRT(QFLOAT(L))*FACTOR(L-1)
C
C        WRITE(65,'(''L='',I2,3X,''FACTOR(L)='',F40.2)') L,FACTOR(L)
C
      END DO
C

      IF (KIND_L.EQ.0) THEN
C
          IF (NORDER+LORDER.GT.40) STOP 'Maupa oraz Bauvann in NORML'
C                
C                          Normalisation for full-integer L-values
          DO N=0,NORDER
             DO L=1,LORDER
                XNORNL(N,L)=1.0Q0/(FACTOR(L+N)/FACTOR(N))
             END DO
          END DO
C
      ELSE
C                          Normalisation for half-integer L-values
C
          PINUMB=4.Q0*QATAN(1.0Q0)
          SQRTPI=QSQRT(PINUMB)
          FRACTN=1.Q0/SQRTPI
C
          FRACT1=FRACTN
C
          DO J=0,LORDER
             FRACT1=FRACT1*2.0Q0/(2.0Q0*QFLOAT(J)+1.0Q0)
             XNORNL(0,J)=QSQRT(FRACT1)
          END DO
C
          DO I=1,NORDER
             FRACTN=FRACTN*2.0Q0*QFLOAT(I)/(2.0Q0*QFLOAT(I)-1.0Q0)
             FRACT2=1.0Q0
C
             DO K=0,LORDER
                FRACT2=FRACT2*2.0Q0/(2.0Q0*QFLOAT(I)
     *		          +2.0Q0*QFLOAT(K)+1.0Q0)
                XNORNL(I,K)=QSQRT(FRACTN*FRACT2)           
             END DO
          END DO
C
      END IF
C
C=======================================================================
C
      RETURN
      END  
C
C=======================================================================
C
      SUBROUTINE SET_FACTORIAL

      INTEGER MAX_FACTORIAL
      PARAMETER(MAX_FACTORIAL=150)

      REAL*16 FACLOG(0:MAX_FACTORIAL)
      COMMON  /FACLOG_TABLE/ FACLOG

      REAL*16 FN
      INTEGER N

      CALL CPUTIM('FACTOR',1)

      FACLOG(0)=0.Q0
      FN=0.Q0
      DO N=1,MAX_FACTORIAL
         FN=FN+1.Q0
         FACLOG(N)=FACLOG(N-1)+QLOG(FN)
      END DO

      CALL CPUTIM('FACTOR',0)

      RETURN
      END
C
C=======================================================================
C=======================================================================
C
      SUBROUTINE LAGUER(XARGUM,FLAGNL,DFLAGN,NORDER,LORDER,KIND_L,
     *                                              NDIM_L,NDIM_N)
      IMPLICIT REAL*16 (A-H,O-Z)
C
      DIMENSION
     *          FLAGNL(0:NDIM_N,0:NDIM_L),
     *          DFLAGN(0:NDIM_N,0:NDIM_L)
C
C=======================================================================
C
C     This subroutine  calculates the values  of the generalised
C     Laguerre functions at the argument  <xargum>; we calculate
C     simultaneously all the Laguerre functions with the indices
C     n=0,1, ... NORDER, for the largest L-index, and next apply
C     the recurrence 
C                      L^{a-1}_n = L^a_n - L^a_{n-1}
C     Attention!
C               Here we calculate only those Laguerre functions
C               that are needed as solutions of the Schrodinger
C               equation for the spherically-symmetric harm-osc
C               potential i.e. alpha=L+1/2  
C
C     KIND_L - either 0 (integer L-indices) or 1 (half-integer)
C
C=======================================================================
C
      CALL CPUTIM('LAGUER',1)
C
C=======================================================================
C
      IF (NORDER.GT.NDIM_N) THEN
          WRITE(0,'(''NORDER='',I3,3X,''NDIM_N='',I3)') NORDER,NDIM_N
          STOP 'NORDER > NDIM_N in LAGUER'
      END IF
C
      IF (NORDER.LT.000000) STOP 'NORDER negative in LAGUER'
C
      IF (LORDER.GE.NDIM_L) STOP 'LORDER geNDIM_L in LAGUER'
      IF (LORDER.LT.000000) STOP 'LORDER negative in LAGUER'
C
      IF (.NOT.(KIND_L.EQ.0.OR.KIND_L.EQ.1)) THEN
          STOP 'Param. KIND_L must be zero or one in LAGUER'
      END IF
C
C=======================================================================
C
C     Attention: Below, removing +0.5 allows to calculate the Lag.
C                functions for the integer order-L polynomials ...
C
      L_HELP=LORDER+1
C
      IF (KIND_L.EQ.0) THEN
          XALPHA=QFLOAT(L_HELP)
      ELSE
          XALPHA=QFLOAT(L_HELP)+0.5Q0
      END IF
C
      X_LAG0=1.0Q0
      X_LAG1=1.0Q0+XALPHA-XARGUM
C
      FLAGNL(0,L_HELP)=X_LAG0
      FLAGNL(1,L_HELP)=X_LAG1
C
C=======================================================================
C
C            In what follows we are using the recurrence:
C
C     L^a_{n+1}=\frac{1}{n+1}[(2n+a+1-x) L^a_n - (n+a) L^a_{n-1}]
C
C     rewritten as
C
C     L^a_{n} = \frac{1}{n}[(2n+a-1-x) L^a_{n-1} - (n-1+a) L^a_{n-2}] 
C
C=======================================================================
C
      DO N=2,NORDER
C
         X_LAG2=((QFLOAT(2*N-1)+XALPHA-XARGUM) *X_LAG1
     *         - (QFLOAT(N-1)  +XALPHA       ) *X_LAG0)/QFLOAT(N)
C
         X_LAG0=X_LAG1
         X_LAG1=X_LAG2
C
         FLAGNL(N,L_HELP)=X_LAG2
C
      END DO        
C_______________________________________________________________________
C
C     Defining generalissd Laguerre polynomials n=0 for all L-indices
C
      DO L=0,L_HELP
         FLAGNL(0,L)=1.0Q0
         DFLAGN(0,L)=0.0Q0
      END DO
C_______________________________________________________________________
C
C     Applying the recurrence: L^{a-1}_n = L^a_n - L^a_{n-1}; n>0
C
      DO L=L_HELP,1,-1
         DO N=NORDER,1,-1
            FLAGNL(N,L-1)=FLAGNL(N,L)-FLAGNL(N-1,L)
         END DO
      END DO
C
C_______________________________________________________________________
C
C     Calculating derivatives of the generld. Laguerre functions
C
      DO L=LORDER,0,-1
         DO N=NORDER,1,-1
            DFLAGN(N,L)=-FLAGNL(N-1,L+1) ! dL^a_n=-L^{a+1}_{n-1}
         END DO
      END DO
C
C=======================================================================
C
      CALL CPUTIM('LAGUER',0)
C
C=======================================================================
C      
      RETURN 
      END 
C
C=======================================================================
C=======================================================================
C
      SUBROUTINE LAGAUS(NGAUSS,XNODES,WEIGHT,ALPHAG,BAUXIL,CAUXIL,
     *                                              EPSLAG,NDLAGR)
      IMPLICIT  REAL*16 (A-H,O-Z)
C
      DIMENSION
     *          XNODES(NDLAGR),WEIGHT(NDLAGR)
      DIMENSION
     *          BAUXIL(NDLAGR),CAUXIL(NDLAGR)
C
      COMMON
     *         /PRFILE/ ISCREN,LOGFIL,LOGCHP,LOGCHN,LOGHIP,LOGHIN,
     *                  LOGPRO,LOGNEU,LOGWRI,LSCREN,IRESUL,ICONVE,
     *                  NOUTPT,LOGAUX,LOGENE,LOGRAD,LOGGRD,LOGMSH,
     *                                                     LOGBIS                                                  
C
C=======================================================================
C
C     LAGAUS - calculates zeros and weights of the Gauss-Laguerre
C              integration formula with the L-dependent weights
C
C=======================================================================
C
C     Karolina ???
C @@@
C ### attention - here we fix LINDEX from alpha=l+1/2
C
      FLOATN=QFLOAT(NGAUSS)
C
      LINDEX=ALPHAG-0.5Q0
C
      CSXSUM=0.0Q0
      CSWSUM=0.0Q0
C
      PINUMB=4.Q0*QATAN(1.0Q0)
      SQRTPI=QSQRT(PINUMB)
C
      IF (LINDEX.EQ.-1) THEN
          FRECUR=SQRTPI
      ELSE
          FRECUR=SQRTPI/2.000Q0
      END IF
C
      DO I=1,LINDEX
         FRECUR=FRECUR*(2.0Q0*QFLOAT(I)+1.0Q0)/2.0Q0
      END DO
C
      DO J=2,NGAUSS
         FRECUR=FRECUR*CAUXIL(J)
      END DO
C
C=======================================================================
C
      DO I=1,NGAUSS
C
         IF (I-1) 6,2,3
C
C            Smallest zero
C
   2     CONTINUE    
C
         X_TERM=(1.0Q0+ALPHAG)*(3.000Q0+.920Q0*ALPHAG)
     *         /(1.0Q0+2.40000Q0*FLOATN+1.80Q0*ALPHAG)
C
         GO TO 6
C
   3     CONTINUE
C   
         IF (I-2) 6,4,5
C
C            Second zero
C
   4     CONTINUE
C   
         X_TERM=X_TERM+(15.0Q0+6.2500Q0*ALPHAG)
     *         /(1.000Q0+.9000Q0*ALPHAG+2.5000Q0*FLOATN)
C
         GO TO 6
C
C            All other zeros
C
   5     CONTINUE
C   
         FLOATI=QFLOAT(I-2)
C   
         R1TERM=(1.00Q0+2.5500Q0*FLOATI)/(1.9Q0*FLOATI)
         R2TERM= 1.26Q0*FLOATI*ALPHAG /(1.0Q0+3.5000Q0*FLOATI)
C  
         RATIOS=(R1TERM+R2TERM)/(1.000Q0+0.3000Q0*ALPHAG)
C
         X_TERM= X_TERM+RATIOS*(X_TERM-XNODES(I-2))
C
   6     CONTINUE
   
         XTITER=X_TERM
         ALITER=ALPHAG
         EPSITE=EPSLAG
C
         IF (LOGWRI.GT.5) THEN
C
             WRITE(LOGFIL,'(12X,''Entering LGROOT from LAGAUS, '',
     *                          ''I_NODE='',I3)') I    
         END IF
C
         CALL LGROOT(XTITER,NGAUSS,ALITER,DPNAUX,XPOLYN,BAUXIL,CAUXIL,
     *                                                  EPSITE,NDLAGR)
         X_TERM=XTITER
C
         XNODES(I)=X_TERM
         WEIGHT(I)=FRECUR/DPNAUX/XPOLYN
C
         CSXSUM=CSXSUM+X_TERM
         CSWSUM=CSWSUM+WEIGHT(I)
C
      END DO
C
C=======================================================================
C
      IF (LOGWRI.GT.4) THEN
          WRITE(LOGFIL,'(9X,''Exiting  LAGAUS'')')
      END IF
C
C=======================================================================
C
      RETURN
      END
C
C=======================================================================
C=======================================================================
C
      SUBROUTINE LGROOT(XARGUM,NGAUSS,ALPHAG,DPNAUX,XPOLYN,
     *                         BAUXIL,CAUXIL,EPSLAG,NDLAGR)
C
      IMPLICIT  REAL*16 (A-H,O-Z)
C
      DIMENSION
     *          BAUXIL(NDLAGR),CAUXIL(NDLAGR)
C
      COMMON
     *         /PRFILE/ ISCREN,LOGFIL,LOGCHP,LOGCHN,LOGHIP,LOGHIN,
     *                  LOGPRO,LOGNEU,LOGWRI,LSCREN,IRESUL,ICONVE,
     *                  NOUTPT,LOGAUX,LOGENE,LOGRAD,LOGGRD,LOGMSH,
     *                                                     LOGBIS                                                  
C
C=======================================================================
C
C     This SUBROUTINE improves approximation root XARGUM
C     In addition we also obtain:
C
C                    DPNAUX=derivative of P(N) at XARGUM
C                and XPOLYN=value of P(N-1) at XARGUM
C
C=======================================================================
C
      ITERAT=0
C
   1  CONTINUE
C
C=======================================================================
C
      ITERAT=ITERAT+1
C
      CALL LGRECR(PNITER,DPITER,XPOLYN,XARGUM,NGAUSS,ALPHAG,
     *                                 BAUXIL,CAUXIL,NDLAGR)
C
      D=PNITER/DPITER
      XARGUM=XARGUM-D
C
      IF (ABS(D/XARGUM)-EPSLAG) 3,3,2
C
   2  IF (ITERAT-10) 1,3,3
C
   3  CONTINUE
C
C=======================================================================
C
      IF (LOGWRI.GT.5) THEN
C
          IF (ITERAT.GE.10) THEN
              WRITE(LOGFIL,'(16X,''From LGROOT ITERAT='',i3,1X,
     *                           ''maximum 10 are allowed, '',
     *                           ''precision EPSLAG='',D12.5,1X,
     *                           ''not respected'')') ITERAT,EPSLAG
          ELSE
              WRITE(LOGFIL,'(16X,''From LGROOT ITERAT='',i3,1X,
     *                           ''maximum 10 are allowed, '',
     *                           ''precision EPSLAG='',D12.5,1X,
     *                           ''satisfied'')') ITERAT,EPSLAG
          END IF
C
      END IF
C
C=======================================================================
C   
      DPNAUX=DPITER
C
C=======================================================================
C
      RETURN
      END
C
C=======================================================================
C=======================================================================
C
      SUBROUTINE LGRECR(PNITER,DPITER,XPOLYN,XARGUM,NGAUSS,ALPHAG,
     *                                       BAUXIL,CAUXIL,NDLAGR)
      IMPLICIT  REAL*16 (A-H,O-Z)
C
      DIMENSION
     *          BAUXIL(NDLAGR),CAUXIL(NDLAGR)
C
C=======================================================================
C     Auxiliary routine - recurrences (called from LGROOT)
C=======================================================================
C
      PLAGUE=1.0000Q0
C
      P_ITER=XARGUM-ALPHAG-1.0Q0
C
      DP1AUX=0.0Q0
      DP_AUX=1.0Q0
C
      DO J=2,NGAUSS
C
         Q_AUXI=(XARGUM-BAUXIL(J))*P_ITER-CAUXIL(J)*PLAGUE
         DQAUXI=(XARGUM-BAUXIL(J))*DP_AUX-CAUXIL(J)*DP1AUX+P_ITER
C
         PLAGUE=P_ITER
         P_ITER=Q_AUXI
         DP1AUX=DP_AUX
         DP_AUX=DQAUXI
C
      END DO
C
      PNITER=P_ITER
      DPITER=DP_AUX
      XPOLYN=PLAGUE
C
C=======================================================================
C
      RETURN
      END
C       
C=======================================================================
C=======================================================================
C       
      SUBROUTINE ORDHEA(VECTOR,INDEXV,MAXVEC,NDVECT)
      DIMENSION
     *          VECTOR(1:NDVECT),INDEXV(1:NDVECT)
C
C=======================================================================
C
C        This routine rearranges vector 'VECTOR' of dimension 'MAXVEC'
C        so that it is in ascending order on return.   Vector 'INDEXV'
C        contains on return,  in its k-th element, the original index
C        of resulting VECTOR(K).
C
C        On input 'INDEXV' must contain original indices.
C
C=======================================================================
C                                   Method used: heapsort of Williams
C=======================================================================
C                  Based on  "numerical  recipes"
C=======================================================================
C
      IF (MAXVEC.LT.1) THEN
          STOP 'ORDHEA called with negative number of vector elements'
      END IF
C
      IF (MAXVEC.EQ.1) RETURN
C
      IF (MAXVEC.GT.NDVECT) STOP 'MAXVEC > NDVECT in ORDHEA'
C
C=======================================================================
C
      LACTIV=MAXVEC/2+1
      INDACT=MAXVEC
C
C=======================================================================
C
    1 CONTINUE
C
C=======================================================================
C
      IF (LACTIV.GT.1) THEN
C
          LACTIV=LACTIV-1
C
          VELEMT=VECTOR(LACTIV)
          IELEMT=INDEXV(LACTIV)
C
      ELSE
C
          VELEMT=VECTOR(INDACT)
          IELEMT=INDEXV(INDACT)
C
          VECTOR(INDACT)=VECTOR(1)
          INDEXV(INDACT)=INDEXV(1)
C
          INDACT=INDACT-1
C
          IF (INDACT.EQ.1) THEN
              VECTOR(1)=VELEMT
              INDEXV(1)=IELEMT
C
              RETURN
C
          END IF
C
      END IF
C
C=======================================================================
C
      I=LACTIV
      J=LACTIV+LACTIV
C
C=======================================================================
C
    2 CONTINUE
C
C=======================================================================
C
      IF (J.LE.INDACT) THEN
C
          IF (J.LT.INDACT.AND.VECTOR(J).LT.VECTOR(J+1)) J=J+1
C
          IF (VELEMT.LT.VECTOR(J)) THEN
C
              VECTOR(I)=VECTOR(J)
              INDEXV(I)=INDEXV(J)
C
              I=J
              J=J+J
C
          ELSE
              J=INDACT+1
          END IF
C
          GO TO 2
C
      END IF
C
C=======================================================================
C
      VECTOR(I)=VELEMT
      INDEXV(I)=IELEMT
C
C=======================================================================
C
      GO TO 1
C
C=======================================================================
C
      END
C      
C=======================================================================
C=======================================================================
C
      SUBROUTINE ORDLAB(NDIM_N,NDIM_L,LABORD_PROTON,LABORD_NEUTRS,
     *                                LABSYM_PROTON,LABSYM_NEUTRS,
     *                                              NDNUCL,INUCLI)
C      
      INCLUDE  'MATDIM/NDIMaN.f'
      
      PARAMETER
     *         (LEVDIM=30,NDIMaL=NDIMaN)
C
      CHARACTER 
     *          LABSYM_PROTON*6,LABSYM_NEUTRS*6
      CHARACTER 
     *          KEYWOR*7,STRING*6,FILEXP*256
C
      DIMENSION
     *          LABORD_PROTON(1:NDNUCL,0:NDIM_N,0:NDIM_L,0:2*NDIM_L+1),
     *          LABORD_NEUTRS(1:NDNUCL,0:NDIM_N,0:NDIM_L,0:2*NDIM_L+1)
      DIMENSION
     *          LABSYM_PROTON(1:NDNUCL,0:NDIM_N,0:NDIM_L,0:2*NDIM_L+1),
     *          LABSYM_NEUTRS(1:NDNUCL,0:NDIM_N,0:NDIM_L,0:2*NDIM_L+1)
      DIMENSION
     *          LABORD(0:NDIMaN,0:NDIMaL,0:2*NDIMaL+1)
C
      COMMON
     *       /PRFILE/ ISCREN,LOGFIL,LOGCHP,LOGCHN,LOGHIP,LOGHIN,
     *                LOGPRO,LOGNEU,LOGWRI,LSCREN,IRESUL,ICONVE,
     *                NOUTPT,LOGAUX,LOGENE,LOGRAD,LOGGRD,LOGMSH,
     *                                                   LOGBIS                                                  
      COMMON
     *       /NUCLEU/ IZ_FIX,IN_FIX,ISOSPI    
C      
      DATA
     *        N_UNIT / 10 /
C
C=======================================================================     
C     The purpose of this routine is ??? byc moze do wywalenia
C=======================================================================     
C
      LABORD(0,0, 1)=1 ! 1s1/2
C                      !        N=0 full
      LABORD(0,1, 3)=1 ! 1p3/2
      LABORD(0,1, 1)=1 ! 1p1/2
C                      !        N=1 full
      LABORD(0,2, 5)=1 ! 1d5/2
C                      !        Intruder from N=2
      LABORD(1,0, 1)=1 ! 2s2/2
      LABORD(0,2, 3)=1 ! 1d3/2
C                      !        N=2 full
      LABORD(0,3, 7)=1 ! 1f7/2
C                      !        Intruder from N=3
      LABORD(1,1, 3)=1 ! 2p3/2 
      LABORD(0,3, 5)=1 ! 1f5/2  
      LABORD(1,1, 1)=1 ! 2p1/2
C                      !        N=3 full
      LABORD(0,4, 9)=1 ! 1g9/2
C                      !        Intruder from N=4  
      LABORD(0,4, 7)=1 ! 1g7/2 
      LABORD(1,2, 5)=1 ! 2d5/2
      LABORD(1,2, 3)=1 ! 2d3/2 
      LABORD(2,0, 1)=1 ! 3s1/2 
C                      !        N=4 full
      LABORD(0,5,11)=1 ! 1h11/2
C                      !        Intruder from N=5  
      LABORD(0,5, 9)=1 ! 1h9/2 
      LABORD(1,3, 7)=1 ! 2f7/2 
      LABORD(1,3, 5)=1 ! 2f5/2 
      LABORD(2,1, 3)=1 ! 3p3/2 
      LABORD(2,1, 1)=1 ! 3p1/2
C                      !        N=5 full 
      LABORD(0,6,13)=1 ! 1i13/2
C                      !        Intruder from N=6 
      LABORD(1,4, 9)=1 ! 2g9/2 
      LABORD(2,2, 5)=1 ! 3d5/2 
      LABORD(0,6,11)=1 ! 1i11/2 
      LABORD(1,4, 7)=1 ! 2g7/2 
      LABORD(3,0, 1)=1 ! 4s1/2 
      LABORD(2,2, 3)=1 ! 3d3/2 
C                      !        N=6 full 
      LABORD(0,7,15)=1 ! 1j15/2 
C                      !        Intruder from N=7    
C
C=======================================================================
C      
      FILEXP='ws_exper/levels.new.d'
C      
      OPEN (N_UNIT,FILE=FILEXP,STATUS='UNKNOWN',FORM='FORMATTED')
C
      IF (LOGWRI.GT.5) THEN
C
          WRITE(LOGFIL,'(/,9X,''Opening UNIT='',I2,1x,
     *                        ''File=ws_exper/levels.new.d'',/)') N_UNIT
C
          WRITE(LOGFIL,'(  9X,''IZ_FIX='',I3,'' IN_FIX='',I3,1X,
     *                        ''ISOSPI='',I1,/)')                   
     *                          IZ_FIX,         IN_FIX,ISOSPI
      END IF
C      
      DO N=0,NDIM_N
         DO L=0,NDIM_L
            DO J=0,2*NDIM_L+1
               LABORD_PROTON(INUCLI,N,L,J)=0
               LABORD_NEUTRS(INUCLI,N,L,J)=0	       
            END DO
         END DO
      END DO
C
C=======================================================================
C
   2  CONTINUE
C
C=======================================================================
C   
      READ (N_UNIT,'(A7)',END=3)  KEYWOR
C
C=======================================================================
C      
      IF (KEYWOR.EQ.'n>>l>>j') THEN
C      
          NUMOCC_PROTON=0
          NUMOCC_NEUTRS=0
C	  
          DO I=1,LEVDIM
C
C             WRITE(0,*)I,IN_FIX,NUMOCC_NEUTRS
C             WRITE(0,*)I,IZ_FIX,NUMOCC_PROTON
C
             IF (NUMOCC_NEUTRS.LT.IN_FIX) THEN
C             
                 READ (N_UNIT,*) NNUMBR,LNUMBR,JNUMBR,STRING
C
                 IF (LOGWRI.GT.5)                 
     *               WRITE(LOGFIL,'(7X,3I3,2X,A)') NNUMBR,LNUMBR,JNUMBR,
     *                                                           STRING            
                 LABORD_NEUTRS(INUCLI,NNUMBR-1,LNUMBR,JNUMBR)=1
                 LABSYM_NEUTRS(INUCLI,NNUMBR-1,LNUMBR,JNUMBR)=STRING
C
                 IF (LABORD_NEUTRS(INUCLI,NNUMBR-1,LNUMBR,JNUMBR).NE.
     *                      LABORD(NNUMBR-1,LNUMBR,JNUMBR)) THEN
C
                     IF (ISCREN.NE.0) THEN
C
                         WRITE(LSCREN,'(9X,''LABORD_N('',I1,'','',I1,
     *                            '','',I2,'')='',I1,'' LABORD='',I1)')
     *                                  NNUMBR-1,LNUMBR,JNUMBR,
     *                    LABORD_NEUTRS(INUCLI,NNUMBR-1,LNUMBR,JNUMBR),
     *                                  LABORD(NNUMBR-1,LNUMBR,JNUMBR)
                     END IF
C
                     STOP 'Occupation structure inconsistency in LABORD'            
C
                 END IF             
C             
                 NUMOCC_NEUTRS=NUMOCC_NEUTRS+JNUMBR+1
C          
                 IF (NUMOCC_PROTON.LT.IZ_FIX) THEN
C             
                     LABORD_PROTON(INUCLI,NNUMBR-1,LNUMBR,JNUMBR)=1
                     LABSYM_PROTON(INUCLI,NNUMBR-1,LNUMBR,JNUMBR)=STRING
C		     
                     NUMOCC_PROTON=NUMOCC_PROTON+JNUMBR+1
C
                 END IF
             ELSE
                 GO TO 1
             END IF
          END DO
      END IF
C     
C=======================================================================
C     
   1  CONTINUE   
C     
C=======================================================================
C       
      IF (KEYWOR.EQ.'ENDFILE') THEN
C
          IF (LOGWRI.GT.3) WRITE(LOGFIL,'()')
C 
          LBPCNT=0
          LBNCNT=0
C
          DO N=0,NDIM_N
             DO L=0,NDIM_L
                DO J=1,2*NDIM_L+1,2
C
                   IF (LABORD_PROTON(INUCLI,N,L,J).EQ.1) THEN
C
                       LBPCNT=LBPCNT+1                   
C
                       IF (LOGWRI.GT.5)
     *                       
     *                     WRITE(LOGFIL,'(9X,''p'',3I3,I5,2x,a)') 
     *                           N,L,J,LABORD_PROTON(INUCLI,N,L,J),
     *                                 LABSYM_PROTON(INUCLI,N,L,J)
                   END IF
C
                   IF (LABORD_NEUTRS(INUCLI,N,L,J).EQ.1) THEN
C
                       LBNCNT=LBNCNT+1                   
C                   
                       IF (LOGWRI.GT.5)
     *                       
     *                     WRITE(LOGFIL,'(9X,''n'',3I3,I5,2x,a)') 
     *                           N,L,J,LABORD_NEUTRS(INUCLI,N,L,J),
     *                                 LABSYM_NEUTRS(INUCLI,N,L,J)
                   END IF
C                   
                END DO ! end N
             END DO ! end L
          END DO ! end J
C       
          CLOSE(N_UNIT)
C       
          IF (LOGWRI.GT.5) THEN
C
              WRITE(LOGFIL,'()')
              WRITE(LOGFIL,'(9X,''In total LBPCNT='',I2,1x,
     *                          ''proton  labels and'',/,9X, 
     *                          ''         LBNCNT='',I2,1x,
     *                          ''neutron labels'')')LBPCNT,LBNCNT
              WRITE(LOGFIL,'(/,9X,''Closing UNIT='',i2,/)') N_UNIT
C
              WRITE(LOGFIL,'(/,''Exiting  ORDLAB'')')
C
          END IF
C
          RETURN
C
      END IF
C               
C=======================================================================
C      
      GO TO 2
C               
C=======================================================================
C
   3  CONTINUE
C
       WRITE(100,'()')
       
       DO N=0,NDIM_N
          DO L=0,NDIM_L
             DO J=0,2*NDIM_L+1
 	              IF (LABORD_PROTON(INUCLI,N,L,J).EQ.1) THEN
C
                    WRITE(LSCREN,'(''p'',4I3)')
     *                               N,L,J,LABORD_PROTON(INUCLI,N,L,J)
                END IF
             END DO
          END DO
       END DO
       
       WRITE(100,'()')
C_______________________________________________________________________
C       
       DO N=0,NDIM_N
          DO L=0,NDIM_L
             DO J=0,2*NDIM_L+1
 	              IF (LABORD_NEUTRS(INUCLI,N,L,J).EQ.1) THEN
C
                    WRITE(LSCREN,'(''n'',4I3)')
     *                               N,L,J,LABORD_NEUTRS(INUCLI,N,L,J)
                END IF
             END DO
          END DO
       END DO
C      
C=======================================================================
C
      IF (LOGWRI.GT.5) THEN
          WRITE(LOGFIL,'(/,''Exiting  ORDLAB'')')
      END IF
C      
C=======================================================================
C  
      RETURN
      END
C
C=======================================================================
C=======================================================================
C
      SUBROUTINE CHEKHO(POLYNL,DPOLYN,POLYNX,RNODES,RWEIGH,HOMEGA,
     *                  NDBASE,NDIM_N,NDIM_L,NSHELL,NDGAUS,NGAUSS,
     *                                       IALPHA,NDNUCL,INUCLI)
C
      DIMENSION
     *          ENKIN4(1:NDBASE,1:NDBASE,0:NDIM_L),
     *          ENPHO4(1:NDBASE,1:NDBASE,0:NDIM_L),
     *          ENPOHO(1:NDBASE,1:NDBASE,0:NDIM_L)
      DIMENSION
     *          POLYNL(1:NDGAUS,0:NDIM_N,0:NDIM_L), 
     *          POLYNX(1:NDGAUS,0:NDIM_N,0:NDIM_L), 
     *          DPOLYN(1:NDGAUS,0:NDIM_N,0:NDIM_L) 
      DIMENSION
     *          RNODES(1:NDGAUS,0:NDIM_L),
     *          RWEIGH(1:NDGAUS,0:NDIM_L)
      DIMENSION
     *          HOMEGA(1:NDNUCL)
C
      COMMON
     *         /PRFILE/ ISCREN,LOGFIL,LOGCHP,LOGCHN,LOGHIP,LOGHIN,
     *                  LOGPRO,LOGNEU,LOGWRI,LSCREN,IRESUL,ICONVE,
     *                  NOUTPT,LOGAUX,LOGENE,LOGRAD,LOGGRD,LOGMSH,
     *                                                     LOGBIS  
C
C=======================================================================
C
C     This an intrinsic-test routine checking the fact that the total 
C     HO hamiltonian  (the sum of the kinetic and potential energies) 
C     is diagonal in the  {|n l s j m_j>}  basis,  the diagonal terms 
C     being the eigen-energies E_nl=(2n+l+3/2)*homega.
C
C     Furthermore according to the virial theorem, the diagonal terms 
C     of the kinetic and potential energies must be equal, both being 
C     equal to 0.5*E_nl.  In order to perform the accuracy tests,  it 
C     will be of advantage to express the matrix elements in units of 
C     HOMEGA.
C
C     This subroutine is identical to KINMAT but here the calculation 
C     of the HO potential has been included,  and the discussed check 
C     performed  together  with the calculation  of the  HO potential 
C     matrix elements, and therefore, the kinetic energy test  (valid 
C     then for all potentials) is performed.
C                                                           H.M. 2008 
C=======================================================================     
C
C=======================================================================     
C     This subroutine  calculates  the matrix elements of the kinetic
C     energy term  [the H.O. scaling factors  may be dependent on the 
C     nucleonic mass]
C=======================================================================     
C
C @@@ ???
CID      EPSTST=1.Q-29 ! These two definitions of EPSTST serve
CID      EPSTST=1.Q-28 ! when working with quadruple precision
C
      EPSTST=1.E-12 ! Working in double precision
C
      HOMAUX=HOMEGA(INUCLI)
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C=======================================================================     
C
      DO LQNUMB=NSHELL,0,-1
         DO NQNUMB=0,(NSHELL-LQNUMB)/2
            DO I_PNTS=1,NGAUSS
C
               ZNODES=RNODES(I_PNTS,LQNUMB)
C
               TERM_1=0.50*(FLOAT(LQNUMB+1)-ZNODES)
C
               POLYNX(I_PNTS,NQNUMB,LQNUMB)=TERM_1
     *	                                    *POLYNL(I_PNTS,NQNUMB,LQNUMB)
     *                                     +ZNODES
     *                                     *DPOLYN(I_PNTS,NQNUMB,LQNUMB)
            END DO
         END DO
      END DO
C_______________________________________________________________________
C
      DO LQNUMB=NSHELL,0,-1
C
C        Here we are going to construct the T-sub-block at given L
C
         ETOTA4_MAXABS=-99999.0
         ENEDI4_MAXABS=-99999.0
C
         DO N1NUMB=0,(NSHELL-LQNUMB)/2
C
            IF (N1NUMB+1.GT.NDBASE) THEN
                STOP 'N1NUMB+1.GT.NDBASE in CHEKHO'
            END IF
C
            DO N2NUMB=0,(NSHELL-LQNUMB)/2       
C       
               T_ELEM=0.00 ! Kinetic Energy
               HOELEM=0.00 ! Harmonic Oscillator Potential
C       
               DO I_PNTS=1,NGAUSS
C
                  IF (IALPHA.EQ.0) THEN
C
                      ZNODES=RNODES(I_PNTS,LQNUMB)
C
                      TERM_1=POLYNX(I_PNTS,N1NUMB,LQNUMB)
     *                      *POLYNX(I_PNTS,N2NUMB,LQNUMB)
C
                      TERM_2=0.250*FLOAT(LQNUMB*(LQNUMB+1))
     *		                    *POLYNL(I_PNTS,N1NUMB,LQNUMB)
     *                      *POLYNL(I_PNTS,N2NUMB,LQNUMB)
C
                      TERMHO=POLYNL(I_PNTS,N1NUMB,LQNUMB)
     *                      *ZNODES*ZNODES
     *                      *POLYNL(I_PNTS,N2NUMB,LQNUMB)
C
                      T_ELEM=T_ELEM+(TERM_1+TERM_2)
     *                      *RWEIGH(I_PNTS,LQNUMB)
C                 
                      HOELEM=HOELEM+TERMHO*RWEIGH(I_PNTS,LQNUMB)
C
                  END IF
C_______________________________________________________________________
C
                  IF (IALPHA.EQ.1) THEN
C
                      ZNODES=RNODES(I_PNTS,LQNUMB)
C 
                      TERM_1=POLYNX(I_PNTS,N1NUMB,LQNUMB)
     *                      *POLYNX(I_PNTS,N2NUMB,LQNUMB)/ZNODES
C
                      TERM_2=0.250*FLOAT(LQNUMB*(LQNUMB+1))
     *		                    *POLYNL(I_PNTS,N1NUMB,LQNUMB)
     *                      *POLYNL(I_PNTS,N2NUMB,LQNUMB)/ZNODES
C
                      TERMHO=POLYNL(I_PNTS,N1NUMB,LQNUMB)
     *                      *POLYNL(I_PNTS,N2NUMB,LQNUMB)*ZNODES
C 
                      T_ELEM=T_ELEM
     *                      +
     *                      (TERM_1+TERM_2)*RWEIGH(I_PNTS,LQNUMB)
C                 
                      HOELEM=HOELEM+TERMHO*RWEIGH(I_PNTS,LQNUMB)
C
                  END IF ! IALPHA=1
C
               END DO ! I_PNTS =>> NGAUSS
C
C              The factor of 2.0 in the following line takes into
C              account the fact that  XNORNL  are defined without
C              the factor sqrt(2). XNORNL are defined without the 
C              stretching factor "a" either, but this one cancels 
C              out in the calculations.
C
               ENKIN4(N1NUMB+1,N2NUMB+1,LQNUMB)=2.00*QEXT(HOMAUX)
     *                                              *T_ELEM
C
               ENPHO4(N1NUMB+1,N2NUMB+1,LQNUMB)=2.00*0.50*0.50
     *                                         *HOMAUX
     *                                         *HOELEM
C
C=======================================================================
C
               ETOTA4=ENKIN4(N1NUMB+1,N2NUMB+1,LQNUMB)
     *	              +ENPHO4(N1NUMB+1,N2NUMB+1,LQNUMB)
C
C=======================================================================
C
               IF (N1NUMB.NE.N2NUMB) THEN
                   ABSET4=ABS(ETOTA4)
                   IF (ETOTA4_MAXABS.LT.ABSET4) ETOTA4_MAXABS=ABSET4
                   IF (ABSET4.GT.EPSTST) THEN
                       STOP 'HO oscillator - tot in CHEKHO'
                   END IF
               END IF
C
C=======================================================================
C
               IF (N1NUMB.EQ.N2NUMB) THEN
C
                   EVERI4=(2.0000*FLOAT(N1NUMB)+FLOAT(LQNUMB)+1.50)
     *                   * HOMAUX
C
                   ENEDI4=ABS(ENKIN4(N1NUMB+1,N2NUMB+1,LQNUMB)
     * 		                -0.5*EVERI4)
C
                   IF (ENEDI4_ABSMAX.LT.ENEDI4) ENEDI4_ABSMAX=ENEDI4 
C
                   IF (ENEDI4.GT.EPSTST) THEN
C
                       WRITE(LSCREN ,'(''ENEDI4='',F38.33,3X,
     *		                               ''EPSTST='',F38.33)')
     *                                   ENEDI4,EPSTST
C
                       STOP 'HO oscillator - kinetic energy in CHEKHO'
C
                   END IF
C
                   ENEDI4=ABS(ENPHO4(N1NUMB+1,N2NUMB+1,LQNUMB)
     * 		                -0.5*EVERI4)
C
                   IF (ENEDI4.GT.EPSTST) THEN
C
                       WRITE(0,'(''ENEDI4='',F38.33,3X,
     *		                         ''EPSTST='',F38.33)')
     *                             ENEDI4,EPSTST
C
                       STOP 'HO oscillator - potential energy in CHEKHO'
C
                    END IF
C
                END IF
C
C=======================================================================
C
            END DO
C
         END DO
C
         IF (LOGWRI.GT.0 .AND. LQNUMB.LT.(NSHELL-1)) THEN
             WRITE(LOGFIL,'(9X,''LQNUMB='',I2,'' ETOTA4='',F20.15,1X,
     *                                        '' ENEDI4='',F20.15)')
     *                           LQNUMB,         ETOTA4_MAXABS,
     *                                           ENEDI4_ABSMAX
         END IF
C      
      END DO
C
      IF (LOGWRI.GT.0) THEN
          WRITE(LOGFIL,'(19X,''EPSTST='',F20.15,''  EPSTST='',F20.15)')
     *                         EPSTST,              EPSTST
      END IF
C
C=======================================================================
C
C     Evaluating directly the matrix elements for the HO potential:
C
C     Diagonal matrix elements  of the kinetic energy operator must 
C     always be equal, in the present basis, to the diagonal matrix 
C     elements of the HO potential, both being equal to half of the
C     total energy (virial theorem).
C
C     It follows  that off-diagonal matrix elements  of the kinetic 
C     energy operator are exactly opposite in sign  to those of the 
C     HO potential.
C
C=======================================================================
C
      DO LQNUMB=NSHELL,0,-1
C
         DO N1NUMB=0,(NSHELL-LQNUMB)/2
C
            DO N2NUMB=0,(NSHELL-LQNUMB)/2       
C
               ENPOHO(N1NUMB+1,N2NUMB+1,LQNUMB)=0.00
C
               IF (N1NUMB.EQ.N2NUMB) THEN
C
                   ENPOHO(N1NUMB+1,N1NUMB+1,LQNUMB)
     *		          =0.500*(2.0000*FLOAT(N1NUMB)+FLOAT(LQNUMB)+1.50)
     *            *HOMAUX

               END IF
C
               IF (N1NUMB.EQ.(N2NUMB-1)) THEN
C
                   ENPOHO(N1NUMB+1,N2NUMB+1,LQNUMB)
     *		          =-0.500*HOMAUX
     *            *SQRT(FLOAT(N2NUMB)*(FLOAT(N2NUMB+LQNUMB)+0.50))
C
	              END IF
C
               IF (N1NUMB.EQ.(N2NUMB+1)) THEN
C
                   ENPOHO(N1NUMB+1,N2NUMB+1,LQNUMB)
     *		          =-0.500*HOMAUX
     *                   *SQRT(FLOAT(N2NUMB+1)
     *                       *(FLOAT(N2NUMB+LQNUMB)+1.500))
C
	              END IF
C
C=======================================================================     
C
               IF (N1NUMB.EQ.N2NUMB) THEN
C
                   DIFFMA=ENPOHO(N1NUMB+1,N2NUMB+1,LQNUMB)
     *                   -ENKIN4(N1NUMB+1,N2NUMB+1,LQNUMB)
               ELSE
C
                   DIFFMA=ENPOHO(N1NUMB+1,N2NUMB+1,LQNUMB)
     *                 -(-ENKIN4(N1NUMB+1,N2NUMB+1,LQNUMB))
	           END IF
C
            IF (ABS(DIFFMA).GT.EPSTST) THEN
                STOP 'ABS(DIFFMA).GT.EPSTST in CHEKHO'
            END IF
C
C=======================================================================     
C
            END DO
C
         END DO
C      
      END DO
C
C=======================================================================
C
      IF (LOGWRI.GT.3) THEN
          WRITE(LOGFIL,'(''Exiting  CHEKHO'')')
      END IF
C
C=======================================================================
C
      RETURN
      END
C      
C=======================================================================
C=======================================================================
C 

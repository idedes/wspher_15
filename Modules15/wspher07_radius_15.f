C FILE NAME = wspher07_radius_15.f ! Keep this symbol:    $ident@string$
C
C=======================================================================     
C=======================================================================     
C                R.M.S. RADIUS CALCULATION RELATED ROUTINES
C=======================================================================     
C=======================================================================     
C      
      SUBROUTINE NUCRAD(NDSPEC,LDSPEC,NDGAUS,NGAUSS,NDBASE,POLYNL,
     *                  NDIM_L,NDIM_N,RAD_UP,RAD_DN,RNODES,RWEIGH,
     *                  AOSCIL,RFUNUP,RFUNDN,NSHELL,ENERUP,ENERDN,
     *                  NWSSPH,LWSSPH,JWSSPH,RMSTHE,RPAIRI,RNOPAI,
     *                  LABORD,CMATUP,CMATDN,ENETHE,VCOEUP,VCOEDN,
     *                                       ISOSPI,N_PART,INUCLI)
C
      INCLUDE   'MATDIM/NDNUCL.f'
C
      DIMENSION
     *          LABORD(0:NDIM_N,0:NDIM_L,0:2*NDIM_L+1)
      DIMENSION
     *          RFUNUP(1:NDGAUS,0:NDIM_N,0:NDIM_L),
     *          RFUNDN(1:NDGAUS,0:NDIM_N,0:NDIM_L)
      DIMENSION
     *          CMATUP(0:NDIM_L,1:NDBASE,1:NDBASE),
     *          CMATDN(0:NDIM_L,1:NDBASE,1:NDBASE)
      DIMENSION     
     *          POLYNL(1:NDGAUS,0:NDIM_N,0:NDIM_L)
      DIMENSION
     *          RNODES(1:NDGAUS,0:NDIM_L),
     *          RWEIGH(1:NDGAUS,0:NDIM_L)
      DIMENSION
     *          ENERUP(0:NDIM_L,1:NDBASE),
     *          ENERDN(0:NDIM_L,1:NDBASE)
      DIMENSION
     *          VCOEUP(0:NDIM_L,1:NDBASE),
     *          VCOEDN(0:NDIM_L,1:NDBASE)
      DIMENSION
     *          RAD_UP(0:NDIM_L,0:NDIM_N),
     *          RAD_DN(0:NDIM_L,0:NDIM_N)
      DIMENSION
     *          ENETHE(1:NDSPEC)
      DIMENSION
     *          NWSSPH(1:NDSPEC),
     *          LWSSPH(1:NDSPEC),
     *          JWSSPH(1:NDSPEC)
      DIMENSION
     *          AOSCIL(1:NDNUCL) 
C      
      COMMON
     *         /IFPAIR/ IF_PAI,IFGPAI,IFDPAI
      COMMON
     *         /BCSBCS_ENEMAX/ ENEMAX
      COMMON
     *         /PRFILE/ ISCREN,LOGFIL,LOGCHP,LOGCHN,LOGHIP,LOGHIN,
     *                  LOGPRO,LOGNEU,LOGWRI,LSCREN,IRESUL,ICONVE,
     *                  NOUTPT,LOGAUX,LOGENE,LOGRAD,LOGGRD,LOGMSH,
     *                                                     LOGBIS  
C
      DATA
     *          EPS_V2 / 1.0e-6 /     
C
C=======================================================================     
C     This subroutine calculates the nuclear r.m.s. radius for
C     protons or neutrons,  including (optionally) the pairing  
C=======================================================================     
C
      FACTOR=2.0*AOSCIL(INUCLI)**2
C
C=======================================================================
C     Calculating the single-particle nucleonic wave functions
C     in their spatial representations to be used to integrate         
C=======================================================================
C
      CALL FUNRAD(CMATUP,CMATDN,NGAUSS,NSHELL,NDIM_L,NDBASE,
     *                   NDGAUS,NDIM_N,RFUNUP,RFUNDN,POLYNL)
C           
C=======================================================================
C     Calculating the root-mean-square radius, no pairing case         
C=======================================================================
C           
      RUPAUX=0.
      RDNAUX=0.
      ICUMUL=0
C           
      DO I=1,LDSPEC  !  Do-loop over the spherically-degenerate states
C
         LQNUMB=LWSSPH(I)
         NNUMBR=NWSSPH(I)
         JNUMBR=JWSSPH(I)
C
         JNUMUP=LQNUMB+LQNUMB+1
         JNUMDN=LQNUMB+LQNUMB-1
C
C        LABORD contains the information about the nuclear s.p. levels
C        that should be taken into account in no-pairing case. Even if
C        the levels cross at variance with experiment the right states
C        will be taken  for the summations of the radius contributions 
C
         IF (LABORD(NNUMBR,LQNUMB,JNUMBR).EQ.1) THEN
C
             ICUMUL=ICUMUL+(JWSSPH(I)+1) ! Formally number of particles
C                                          on occupied WS states with (2j+1) dgs
C
C            Below we Gauss-integrate the expectation value of r^2
C                     one WS single-particle state after the other
	            UPPAUX=0.
	            DWNAUX=0.
	 
	            DO I_PNTS=1,NGAUSS
C
                IF (JWSSPH(I).EQ.JNUMUP) THEN 
C
                    UPPAUX=UPPAUX+(RFUNUP(I_PNTS,NNUMBR,LQNUMB)
     *                            *RNODES(I_PNTS,LQNUMB))**2
     *                            *RWEIGH(I_PNTS,LQNUMB)
                ELSE
                    UPPAUX=0.		  
                END IF
C
                IF (JWSSPH(I).EQ.JNUMDN) THEN
C
                    DWNAUX=DWNAUX+(RFUNDN(I_PNTS,NNUMBR,LQNUMB)
     *                            *RNODES(I_PNTS,LQNUMB))**2
     *                            *RWEIGH(I_PNTS,LQNUMB)
                ELSE
                    DWNAUX=0.  
                END IF
C 
             END DO
C
             RUPAUX=RUPAUX+UPPAUX*(LQNUMB+1)
             RDNAUX=RDNAUX+DWNAUX*LQNUMB
C
         END IF
C
      END DO
C
      IF (ICUMUL.NE.N_PART) THEN
C
          WRITE(0,'(/,''No-pairing, cumulative sum must give N_PART='',
     *                                             i3,'' ICUMUL='',i3)')
     *                                             N_PART,ICUMUL
          IF (ISOSPI.EQ.1) THEN
              STOP 'Ballaggane with the radius for protons  in NUCRAD'
          ELSE
              STOP 'Ballaggane with the radius for neutrons in NUCRAD'
          END IF 
C
      END IF
C  
      RADIUS=RUPAUX+RDNAUX  
      RMSRAA=SQRT((RADIUS*FACTOR/N_PART))
      RNOPAI=RMSRAA
C  
      RMSTHE=RMSRAA
C           
C=======================================================================
C
      IF (IF_PAI.EQ.0) THEN
          IF (LOGWRI.GT.5) THEN
              WRITE(LOGFIL,'(15X,''Exiting  NUCRAD, no pairing option''
     *                                                              )')
          END IF 
          RETURN
      END IF           
C           
C=======================================================================
C=======================================================================
C     Calculating the  rms  radius within the BCS pairing      
C=======================================================================
C=======================================================================
C
C     In the case of monopole pairing we take into account 
C     the following two possibilities:
C
C     a. The usual BCS system of 2 equations with given G
C
      IF (IFGPAI.EQ.1) THEN
          CALL PAIRIN(ENERUP,ENERDN,NSHELL,NGAUSS,ENEMAX,
     *                       VCOEUP,VCOEDN,NDIM_L,NDBASE)
      END IF
C
C     b. Calculating at empirical Delta -> finding Lambda
C
      IF (IFDPAI.EQ.1) THEN
          CALL PAIDEL(ENERUP,ENERDN,NSHELL,ENEMAX,
     *                VCOEUP,VCOEDN,NDIM_L,NDBASE,INUCLI)
      END IF
C
C     c. Solving particle number and Delta^(3) equations
C
      IFDEL3=1
      IFDEL3=0
C
      IF (IFDEL3.EQ.1) THEN
C
          CALL PAIDL3(ENERUP,ENERDN,NSHELL,NGAUSS,ENEMAX,
     *                       VCOEUP,VCOEDN,NDIM_L,NDBASE)
      END IF
C           
C=======================================================================
C     The above calls provide the occupation coefficients 
C     calculated  using two alternative approaches;  what
C     we do below,  is independent of the algorithm above        
C=======================================================================
C
      ICUMUL=0
      ACUMUL=0.
      RUPAUX=0.
      RDNAUX=0.
C
      DO I=1,LDSPEC ! The do-loop over the spher.-degenerate states
C
         IF (ENETHE(I).GT.ENEMAX) GO TO 1
C
         LQNUMB=LWSSPH(I)
       	 NNUMBR=NWSSPH(I)
	        JNUMBR=JWSSPH(I)
C
         JNUMUP=LQNUMB+LQNUMB+1
	        JNUMDN=LQNUMB+LQNUMB-1
C
         ICUMUL=ICUMUL+(JWSSPH(I)+1) ! Formally the number of particles
C                                      on all WS states with (2j+1) dgs
         IF (JWSSPH(I).EQ.JNUMUP) THEN 
             ACUMUL=ACUMUL+(JWSSPH(I)+1)*(VCOEUP(LQNUMB,NNUMBR+1))
         END IF
C
         IF (JWSSPH(I).EQ.JNUMDN) THEN 
             ACUMUL=ACUMUL+(JWSSPH(I)+1)*(VCOEDN(LQNUMB,NNUMBR+1))
         END IF
C
C        Below we Gauss-integrate the expectation value of r^2
C                 one WS single-particle state after the other
	        UPPAUX=0.
	        DWNAUX=0.
	 
	        DO I_PNTS=1,NGAUSS
C
            IF (JWSSPH(I).EQ.JNUMUP) THEN 
C
                UPPAUX=UPPAUX+(RFUNUP(I_PNTS,NNUMBR,LQNUMB)
     *                        *RNODES(I_PNTS,LQNUMB))**2
     *                        *RWEIGH(I_PNTS,LQNUMB)
            ELSE
                UPPAUX=0.		  
            END IF
C
            IF (JWSSPH(I).EQ.JNUMDN) THEN
C
                DWNAUX=DWNAUX+(RFUNDN(I_PNTS,NNUMBR,LQNUMB)
     *                        *RNODES(I_PNTS,LQNUMB))**2
     *                        *RWEIGH(I_PNTS,LQNUMB)
            ELSE
                DWNAUX=0.  
            END IF
C 
         END DO
C   
         UPPAUX=UPPAUX*VCOEUP(LQNUMB,NNUMBR+1)
         DWNAUX=DWNAUX*VCOEDN(LQNUMB,NNUMBR+1)
C
	        RUPAUX=RUPAUX+UPPAUX*(LQNUMB+1)
         RDNAUX=RDNAUX+DWNAUX*LQNUMB
C
   1     CONTINUE
C
      END DO
C  
      RADIUS=RUPAUX+RDNAUX  
      RMSRAA=SQRT((RADIUS*FACTOR/N_PART))
      RPAIRI=RMSRAA
C  
      RMSTHE=RMSRAA
C           
C=======================================================================
C
      IF (ABS(ACUMUL-REAL(N_PART)).GT.EPS_V2) THEN
C
          IF (LOGWRI.GT.0) THEN
              WRITE(LOGFIL,'(/,''ACUMUL='',f16.12,'' should be equal '',
     *                         ''to N_PART='',I3,'' with precision '',
     *                         ''EPS_V2='',e10.3)')
     *                           ACUMUL,N_PART,EPS_V2
          END IF
C
      ELSE
C
          IF (LOGWRI.GT.0) THEN
              WRITE(LOGFIL,'(/,''ACUMUL='',f16.12,'' compared '',
     *                         ''to N_PART='',I3,'' with precision '',
     *                         ''EPS_V2='',e10.3)')
     *                           ACUMUL,N_PART,EPS_V2
          END IF
C
      END IF
C           
C=======================================================================
C
      IF (LOGWRI.GT.5) THEN
          WRITE(LOGFIL,'(12X,''Exiting  NUCRAD, after pairing option''
     *                                                             )')
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
      SUBROUTINE FUNRAD(CMATUP,CMATDN,NGAUSS,NSHELL,NDIM_L,NDBASE,
     *                         NDGAUS,NDIM_N,RFUNUP,RFUNDN,POLYNL)
C      
      DIMENSION
     *          CMATUP(0:NDIM_L,1:NDBASE,1:NDBASE),
     *          CMATDN(0:NDIM_L,1:NDBASE,1:NDBASE)
      DIMENSION     
     *          POLYNL(1:NDGAUS,0:NDIM_N,0:NDIM_L)
      DIMENSION
     *          RFUNUP(1:NDGAUS,0:NDIM_N,0:NDIM_L),
     *          RFUNDN(1:NDGAUS,0:NDIM_N,0:NDIM_L)
C
C=======================================================================     
C     This subroutine calculates the spatial representation of the
C     radial wave-functions, solutions to the WS problem, at Gauss
C     nodes - the same that have been used for the matrix-elements
C=======================================================================     
C
      DO I_PNTS=1,NGAUSS
      
         DO LQNUMB=NSHELL,0,-1
C
            NBBASE=(NSHELL-LQNUMB)/2
C
            DO NNUMBR=0,NBBASE
C 
               AUXUPP=0. ! Below: NNUMBR+1 denumerates eigenvectors
               AUXDWN=0. ! while  N1NUMB+1 its expansion components
C
               DO N1NUMB=0,NBBASE
	       
                  AUXUPP=AUXUPP+POLYNL(I_PNTS,N1NUMB,LQNUMB)
     *                         *CMATUP(LQNUMB,N1NUMB+1,NNUMBR+1)
C     
                  AUXDWN=AUXDWN+POLYNL(I_PNTS,N1NUMB,LQNUMB)
     *                         *CMATDN(LQNUMB,N1NUMB+1,NNUMBR+1)

     	       END DO
C
	           RFUNUP(I_PNTS,NNUMBR,LQNUMB)=AUXUPP ! Spatial represent.
	           RFUNDN(I_PNTS,NNUMBR,LQNUMB)=AUXDWN ! for up's & down's
C	       
            END DO
         END DO

      END DO
C
C=======================================================================
C
      RETURN
      END
C      
C=======================================================================
C=======================================================================
C 

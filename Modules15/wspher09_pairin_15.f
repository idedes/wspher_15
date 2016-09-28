C FILE NAME = wspher09_pairin_15.f ! Keep this symbol:    $ident@string$
C      
C=======================================================================
C=======================================================================
C                       PAIRING SUBROUTINE PACKAGE
C=======================================================================
C=======================================================================
C
      SUBROUTINE EQSBCS(BCSLMB,DELTA2,DNOVDL,DNOVDD,DGOVDL,DGOVDD,
     *                                              AN_BCS,AG_BCS)
C
      INCLUDE  'MATDIM/NDSPEC.f'
C
      DIMENSION
     *          NBLOCK(1:NDSPEC)
C
      COMMON
     *         /BCSBCS_RES_WS/ ESINGP(1:NDSPEC)
      COMMON
     *         /BCSBCS_BLCKLV/ NBLOCK
      COMMON
     *         /BCSBCS_LIMPAI/ NOLEVT,INDEXX_MAXSPH
C
C=======================================================================
C
C     This subroutine calculates the rhs of BCS-equations and
C     their derivatives.  It is called  in subroutine  SYSPAI
C
C     Left hand sides are:  The particle number  and 2/G_PAIR
C
C=======================================================================
C
C     BCSLMB - Denotes the pairing BCS - lambda
C
C=======================================================================
C
      DELTA2=MAX(DELTA2,1.0E-9)
C
      AN_BCS=0
      AG_BCS=0
C
      SUMM_0=0
      SUMM_1=0
C
      DO LEVELS=1,NOLEVT
C
         IF (NBLOCK(LEVELS).EQ.1) GO TO 1
C
         ENEDIF=ESINGP(LEVELS)-BCSLMB
         ENER_2=ENEDIF*ENEDIF +DELTA2
C
         ONEOSQ=1/SQRT(ENER_2)
         ONEOE2=ONEOSQ/ENER_2
C
         AN_BCS=AN_BCS-ONEOSQ*ENEDIF+1.0
         AG_BCS=AG_BCS+ONEOSQ
C
         SUMM_0=SUMM_0+ONEOE2
         SUMM_1=SUMM_1+ONEOE2*ENEDIF
C
   1     CONTINUE
C
      END DO
C
      DNOVDL=+SUMM_0*DELTA2
      DNOVDD=+SUMM_1*0.5000
C
      DGOVDL=+SUMM_1
      DGOVDD=-SUMM_0*0.5000
C
C=======================================================================
C
      RETURN
      END
C
C=======================================================================
C=======================================================================
C
      SUBROUTINE EQSBC3(BCSLMB,DELTA2,DNOVDL,DNOVDD,DGOVDL,DGOVDD,
     *                                              AN_BCS,AG_BCS)
C
      INCLUDE  'MATDIM/NDSPEC.f'
C
      DIMENSION
     *          NBLOCK(1:NDSPEC)
C
      COMMON
     *         /BCSBCS_RES_WS/ ESINGP(1:NDSPEC)
      COMMON
     *         /BCSBCS_DELTA3/ DELTA3_INDICT
      COMMON
     *         /BCSBCS_BLCKLV/ NBLOCK
      COMMON
     *         /BCSBCS_LIMPAI/ NOLEVT,INDEXX_MAXSPH
C
C=======================================================================
C
C     This subroutine calculates the rhs  of BCS-LIKE equations and
C     and their derivatives.  It is called  from subroutine  SYSPAI
C
C     Left hand sides are:  The particle number  and  DELTA3_INDICT
C
C=======================================================================
C
C     BCSLMB - Denotes the pairing BCS - lambda
C
C=======================================================================
C
      DELTA2=MAX(DELTA2,1.0E-9)
C
      AN_BCS=0
      AG_BCS=0
C
      SUMM_0=0
      SUMM_1=0
      LEVCNT=0
C
      DO LEVELS=1,NOLEVT
C
         IF (NBLOCK(LEVELS).EQ.1) THEN  ! Only one level must be 
             LEVBLK=LEVELS              !          blocked  here
             LEVCNT=LEVCNT+1
             IF (LEVCNT.GT.1) THEN
                 STOP 'Level count > 1 not allowed here <<= EQSBC3'
             END IF
             GO TO 1
         END IF
C         
         ENEDIF=ESINGP(LEVELS)-BCSLMB
         ENER_2=ENEDIF*ENEDIF +DELTA2   ! Quasiparticle squared
C
         ONEOSQ=1/SQRT(ENER_2)          ! Term  in 2/G equation
         ONEOE2=ONEOSQ/ENER_2           ! 1.0 / quasiparticle^3
C
         AN_BCS=AN_BCS-ONEOSQ*ENEDIF+1.0! sum_k [1-(e-l)/quasip]
C                             Proceeding to calculate derivatives
         SUMM_0=SUMM_0+ONEOE2
         SUMM_1=SUMM_1+ONEOE2*ENEDIF
C
   1     CONTINUE
C
      END DO
C
      ENEDIF=ESINGP(LEVBLK)-BCSLMB
C
      AG_BCS=SQRT(ENEDIF**2+DELTA2)
C
      DNOVDL=SUMM_0*DELTA2
      DNOVDD=SUMM_1*0.5000
C
      DGOVDL=ENEDIF/SQRT(ENEDIF**2+DELTA2) 
C
      DGOVDD=0.5000/SQRT(ENEDIF**2+DELTA2)
      WRITE(66,'(''LEVBLK='',I3,'' NOLEVT='',I3)') LEVBLK,NOLEVT
      WRITE(66,'(''BCSLMB,DELTA2,ENEDIF,='',3F20.12,
     *           '' AG,AN='',2F20.12)') BCSLMB,SQRT(DELTA2),ENEDIF,
     *                                  AG_BCS,AN_BCS
C
C=======================================================================
C
      RETURN
      END
C
C=======================================================================
C=======================================================================
C
      SUBROUTINE SYSPAI(X_INIT,Y_INIT,FFUNCT,GFUNCT,EPS_XY,EPS_FG,
     *                  IERROR,MXITER,ITERAC,DIFF_F,DIFF_G,FNNAME)
C
      PARAMETER
     *         (SHIUPP=1.001,SHIDWN=0.999,DETMIN=1.0E-30)
C
      EXTERNAL
     *          FNNAME
C
      DIMENSION
     *          DXPUSH(0:7),DYPUSH(0:7)
C
C=======================================================================
C
      DATA
     *     DXPUSH /  SHIUPP, SHIUPP, 0.00E0, SHIDWN,
     *               SHIDWN, SHIDWN, 0.00E0, SHIUPP /
      DATA
     *     DYPUSH /  0.00E0, SHIDWN, SHIDWN, SHIDWN,
     *               0.00E0, SHIUPP, SHIUPP, SHIUPP /
C
C=======================================================================
C
C     This subroutine solves system of two non-linear equations of
C     the form:
C                     FN_ACT(XARGUM,YARGUM)=FFUNCT
C                     GN_ACT(XARGUM,YARGUM)=GFUNCT
C
C     using Newton method.
C
C     FFUNCT and GFUNCT are input constants, whereas the functions 
C     FN_ACT and GN_ACT are provided by the external subroutine:
C
C              FNNAME(XARGUM,YARGUM,DFOVDX,DFOVDY,DGOVDX,DGOVDY,
C                                                 FN_ACT,GN_ACT)
C
C     XARGUM - The actual value of the x-argument
C     YARGUM - The actual value of the x-argument
C
C     FN_ACT - The actual value of the f-function
C     GN_ACT - The actual value of the g-function
C
C     X_INIT - on entry => an initial guess for x
C              on exit  => the f i n a l solution
C
C     Y_INIT - on entry => an initial guess for y
C              on exit  => the f i n a l solution
C
C     ITERAC - on exit: if 0 no solution found within the iteration
C                                                             limit
C                       if non zer0 => actual no of iterations used
C
C     MXITER - Maximum number of iterations allowed
C
C     If convergence has not been achieved, ITERAC  is set to 0. In
C     any case  X_INIT and  Y_INIT  contain the best solution found
C     so far, while DIFF_F  and DIFF_G contain  deviation of FN_ACT
C     and GN_ACT from FFUNCT and GFUNCT at this point.  The routine
C     terminates successfully when:
C
C           ABS(FN_ACT-FFUNCT) + ABS(GN_ACT-GFUNCT) .LE. EPS_FG,
C     or
C                  ABS(DELTAX) + ABS(DELTAY) .LE. EPS_XY
C
C=======================================================================
C
      DEVIAT=1.0E+37
      IERROR=0
C
C=======================================================================
C
      DO I=1,MXITER
C
         ITERAC=I
C
         CALL FNNAME(X_INIT,Y_INIT,DFOVDX,DFOVDY,DGOVDX,DGOVDY,
     *                                           FN_ACT,GN_ACT)
C
         F_DIFF=FN_ACT-FFUNCT
         G_DIFF=GN_ACT-GFUNCT
C
         ABSFDI=ABS(F_DIFF)
         ABSGDI=ABS(G_DIFF)
C
         SUMABS=ABSFDI+ABSGDI
C
C=======================================================================
C        Minimising the deviation in terms of functions
C=======================================================================
C
         IF (SUMABS.LT.DEVIAT) THEN
C
             DIFF_F=ABSFDI
             DIFF_G=ABSGDI
             DEVIAT=SUMABS
             X_BEST=X_INIT
             Y_BEST=Y_INIT
C
         END IF
C
C=======================================================================
C        Checking the precision in terms of functions
C=======================================================================
C
         IF (SUMABS.LE.EPS_FG) GO TO 2
C
C=======================================================================
C
         DETERM=DFOVDX*DGOVDY-DFOVDY*DGOVDX
C
C=======================================================================
C        Checking  for possibly vanishing determinant
C=======================================================================
C
         IF (ABS(DETERM).LE.DETMIN) THEN
C
             IAUXIL=MOD(ITERAC-1,8)
C
             X_INIT=DXPUSH(IAUXIL)*X_INIT
             Y_INIT=DYPUSH(IAUXIL)*Y_INIT
C
             GO TO 1
C
         END IF
C
C=======================================================================
C
         DELTAX=(F_DIFF*DGOVDY-G_DIFF*DFOVDY)/DETERM
         DELTAY=(G_DIFF*DFOVDX-F_DIFF*DGOVDX)/DETERM
C
C=======================================================================
C        Checking the precision in arguments (increments small enough?)
C=======================================================================
C
         IF (ABS(DELTAX)+ABS(DELTAY).LE.EPS_XY) GO TO 2
C
C=======================================================================
C        Correct the running arguments XARGUM and YARGUM and proceed
C=======================================================================
C
         X_INIT=X_INIT-DELTAX
         Y_INIT=Y_INIT-DELTAY
C
   1     CONTINUE
C
      END DO
C
C=======================================================================
C
      IERROR=1
      RETURN
C
C=======================================================================
C
   2  CONTINUE
C
      X_INIT=X_BEST
      Y_INIT=Y_BEST
C
C=======================================================================
C
      RETURN
      END
C      
C=======================================================================
C=======================================================================
C
      SUBROUTINE PAIRIN(ENERUP,ENERDN,NSHELL,NGAUSS,ENEMAX,
     *                         VCOEUP,VCOEDN,NDIM_L,NDBASE)
C
      INCLUDE   'MATDIM/NDSPEC.f'
      INCLUDE   'MATDIM/NDIM_N.f'
      INCLUDE   'MATDIM/NDNUCL.f'
      INCLUDE   'MATDIM/NDGAUS.f'
C
      EXTERNAL
     *          EQSBCS
C
      DIMENSION
     *          ENERUP(0:NDIM_L,1:NDBASE),
     *          ENERDN(0:NDIM_L,1:NDBASE)
      DIMENSION
     *          VCOEUP(0:NDIM_L,1:NDBASE),
     *          VCOEDN(0:NDIM_L,1:NDBASE)
C
      COMMON
     *       /DIMSPH/ LDSING_PROTON,LEVTHE_PROTON,
     *                LEVEXP_PROTON(1:NDNUCL),
     *      	      LDSING_NEUTRS,LEVTHE_NEUTRS,
     *                LEVEXP_NEUTRS(1:NDNUCL)
      COMMON
     *       /SPECTR/ ENETHE_PROTON(1:NDSPEC),
     *      	      ENETHE_NEUTRS(1:NDSPEC)
      COMMON
     *       /DEGENE/ IDGLEV_PROTON(1:NDSPEC),
     *                ICUMUL_PROTON(0:NDSPEC),
     *                IDGLEV_NEUTRS(1:NDSPEC),
     *                ICUMUL_NEUTRS(0:NDSPEC)
C
      COMMON
     *       /BCSBCS_ENETHE/ ENETHE(1:NDSPEC)
     *       /BCSBCS_IDGLEV/ IDGLEV(1:NDSPEC)
      COMMON
     *       /BCSBCS_BLCKLV/ NBLOCK(1:NDSPEC)
      COMMON
     *       /BCSBCS_RES_WS/ ESINGP(1:NDSPEC)
      COMMON
     *       /BCSBCS_LIMPAI/ NOLEVT,INDEXX_MAXSPH
C
      COMMON
     *       /NUCLEU/ IZ_FIX,IN_FIX,ISOSPI    
      COMMON
     *       /BCSVAL/ XLAMBD,DELTA2
C
      DATA
     *        EPS_XY / 5.0E-9 /,
     *        EPS_FG / 5.0E-9 /,
     *        EPSILO / 5.0E-8 /
      DATA
     *        MXITER /  100 /
C
C=======================================================================     
C     This subroutine calculates the BCS pairing Delta and Lambda
C     parameters as well as the pairing v2 occupation coefficients
C     assuming the G-pairing strength known
C=======================================================================
C
      IF (ISOSPI.EQ.1) THEN
C
          N_PART=IZ_FIX
          LDSING=LDSING_PROTON      
C                     LDSING = The number of degenerate 
          DO INDEXX=1,LDSING ! spherical  proton-levels 
	     ENETHE(INDEXX)=ENETHE_PROTON(INDEXX)
             IDGLEV(INDEXX)=IDGLEV_PROTON(INDEXX)
          END DO
C
      END IF
C
      IF (ISOSPI.EQ.0) THEN
C
          N_PART=IN_FIX
          LDSING=LDSING_NEUTRS
C                     LDSING = The number of degenerate 
          DO INDEXX=1,LDSING ! spherical neutron-levels
	     ENETHE(INDEXX)=ENETHE_NEUTRS(INDEXX)
             IDGLEV(INDEXX)=IDGLEV_NEUTRS(INDEXX)
          END DO
C
      END IF
C
C=======================================================================
C
      NOLEVT=0
      INDEKS=1
      EPSENE=1.0E-11
C      
C     We artificially remove the spherical degeneracy in order
C     to use the standard BCS equations-solving routine EQSBCS     
C 
      DO INDEXX=1,LDSING
C
         IF (INDEKS.LE.N_PART) THEN
C
             IDEGEN=IDGLEV(INDEXX)/2
             IDEGE2=IDEGEN/2
C
	     DO I=0,IDEGEN-1
                ESINGP(INDEKS)=ENETHE(INDEXX)+(REAL(I)-IDEGE2)*EPSENE
                INDEKS=INDEKS+1
	     END DO
C
	 END IF 
C
      END DO 
C
      NOLEVT=INDEKS-1
C
      ENEMAX=ESINGP(NOLEVT)+EPSENE ! Energy limit for summations
C
      GACTIV=G_FITT(IN_FIX,IZ_FIX,ISOSPI)
C
      A_MASS=REAL(IZ_FIX+IN_FIX)
C
C=======================================================================
C     Blocking of levels is not taken into account
C=======================================================================
C
      DO LEVELS=1,NDSPEC
         NBLOCK(LEVELS)=0
      END DO
C
C=======================================================================
C     BCS (to start with)
C=======================================================================
C
      DELTA2=(12.0/(SQRT(A_MASS)))**2
C
      XLAMBD=0.5*(ESINGP(N_PART/2)+ESINGP(N_PART/2+1))
C
      TWOVRG=2.0/GACTIV
C
      PARTIC=N_PART
C
C=======================================================================
C     Solving the usual BCS-equations
C=======================================================================
C
      CALL SYSPAI(XLAMBD,DELTA2,PARTIC,TWOVRG,EPS_XY,EPS_FG,
     *            IERROR,MXITER,ITERAC,DIFF_F,DIFF_G,EQSBCS)
C
      IF (IERROR.EQ.1) THEN
          XLAMBD=0.5*(ESINGP(N_PART/2)+ESINGP(N_PART/2+1))
          DELTA2=0.01
      END IF
C
C=======================================================================
C      
      DO LQNUMB=NSHELL,0,-1
C
         DO NNUMBR=0,(NSHELL-LQNUMB)/2
C	
	    EMLAUP=ENERUP(LQNUMB,NNUMBR+1)-XLAMBD
	    EMLADN=ENERDN(LQNUMB,NNUMBR+1)-XLAMBD
C	 
            VCOEUP(LQNUMB,NNUMBR+1)
     *     =
     *	    0.5*(1.-EMLAUP/(SQRT(EMLAUP**2+DELTA2)))
C	 
            IF (LQNUMB.NE.0) THEN
C	 
                VCOEDN(LQNUMB,NNUMBR+1)
     *         =
     *          0.5*(1.-EMLADN/(SQRT(EMLADN**2+DELTA2)))
            ELSE
	        VCOEDN(LQNUMB,NNUMBR+1)=0.0
            END IF
C
         END DO
      END DO     
C
C=======================================================================
C     Testing the v2 coefficients...
C=======================================================================
C
      IF (IERROR.EQ.0) THEN
C
          TESTV2=0.
C
          DO INDEXX=1,NOLEVT
      	     EMLAMB=ESINGP(INDEXX)-XLAMBD	       
             V2COEF=0.5*(1.-EMLAMB/(SQRT(EMLAMB**2+DELTA2)))
	     TESTV2=TESTV2+V2COEF
          END DO
C
          TESTV2=2.*TESTV2
          COMPV2=TESTV2-N_PART
C
          IF (ABS(COMPV2).GT.EPSILO) THEN ! STOP 'V2 coefficients!!!'
              WRITE(0,'(''Test v2='',f30.16)')TESTV2
          END IF
C
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
      SUBROUTINE PAIDL3(ENERUP,ENERDN,NSHELL,NGAUSS,ENEMAX,
     *                         VCOEUP,VCOEDN,NDIM_L,NDBASE)
C
      INCLUDE   'MATDIM/NDSPEC.f'
      INCLUDE   'MATDIM/NDIM_N.f'
      INCLUDE   'MATDIM/NDNUCL.f'
      INCLUDE   'MATDIM/NDGAUS.f'
C
      EXTERNAL
     *          EQSBC3
C
      DIMENSION
     *          ENERUP(0:NDIM_L,1:NDBASE),
     *          ENERDN(0:NDIM_L,1:NDBASE)
      DIMENSION
     *          VCOEUP(0:NDIM_L,1:NDBASE),
     *          VCOEDN(0:NDIM_L,1:NDBASE)
C
      COMMON
     *       /DIMSPH/ LDSING_PROTON,LEVTHE_PROTON,
     *                LEVEXP_PROTON(1:NDNUCL),
     *      	      LDSING_NEUTRS,LEVTHE_NEUTRS,
     *                LEVEXP_NEUTRS(1:NDNUCL)
      COMMON
     *       /SPECTR/ ENETHE_PROTON(1:NDSPEC),
     *      	      ENETHE_NEUTRS(1:NDSPEC)
      COMMON
     *       /DEGENE/ IDGLEV_PROTON(1:NDSPEC),
     *                ICUMUL_PROTON(0:NDSPEC),
     *                IDGLEV_NEUTRS(1:NDSPEC),
     *                ICUMUL_NEUTRS(0:NDSPEC)
C
      COMMON
     *       /BCSBCS_ENETHE/ ENETHE(1:NDSPEC)
     *       /BCSBCS_IDGLEV/ IDGLEV(1:NDSPEC)
      COMMON
     *       /BCSBCS_BLCKLV/ NBLOCK(1:NDSPEC)
      COMMON
     *       /BCSBCS_RES_WS/ ESINGP(1:NDSPEC)
C
      COMMON
     *       /NUCLEU/ IZ_FIX,IN_FIX,ISOSPI    
      COMMON
     *       /BCSBCS_DELTA3/ DELTA3_INDICT
      COMMON
     *       /BCSBCS_LIMPAI/ NOLEVT,INDEXX_MAXSPH
      COMMON
     *       /BCSVAL/ XLAMBD,DELTA2
C
      DATA
     *        EPS_XY / 5.0E-9 /,
     *        EPS_FG / 5.0E-9 /,
     *        EPSILO / 5.0E-8 /
      DATA
     *        MXITER /  100 /
C
C=======================================================================     
C     This subroutine calculates the BCS pairing Delta and Lambda
C     parameters as well as the pairing v2 occupation coefficients
C     assuming the pairing Delta^(3) indicator is known (odd n_o)
C=======================================================================
C
      IZ_FIX=19
      DELTA3_INDICT=3.1
C
      IF (ISOSPI.EQ.1) THEN
C
          N_PART=IZ_FIX
          LDSING=LDSING_PROTON      
C                     LDSING = The number of degenerate 
          DO INDEXX=1,LDSING ! spherical  proton-levels 
	     ENETHE(INDEXX)=ENETHE_PROTON(INDEXX)
             IDGLEV(INDEXX)=IDGLEV_PROTON(INDEXX)
          END DO
C
      END IF
C
      IF (ISOSPI.EQ.0) THEN
C
          N_PART=IN_FIX
          LDSING=LDSING_NEUTRS
C                     LDSING = The number of degenerate 
          DO INDEXX=1,LDSING ! spherical neutron-levels
	     ENETHE(INDEXX)=ENETHE_NEUTRS(INDEXX)
             IDGLEV(INDEXX)=IDGLEV_NEUTRS(INDEXX)
          END DO
C
      END IF
C
C=======================================================================
C
      NOLEVT=0
      INDEKS=1
      EPSENE=1.0E-11
C      
C     We artificially remove the spherical degeneracy in order
C     to use the standard BCS equations-solving routine EQSBCS     
C 
      DO INDEXX=1,LDSING
C
         IF (INDEKS.LE.N_PART) THEN
C
             IDEGEN=IDGLEV(INDEXX)/2
             IDEGE2=IDEGEN/2
C
	     DO I=0,IDEGEN-1
                ESINGP(INDEKS)=ENETHE(INDEXX)+(REAL(I)-IDEGE2)*EPSENE
                INDEKS=INDEKS+1
                WRITE(67,'(I3,'' ESINGP='',F14.8)')INDEKS,ESINGP(INDEKS)
	     END DO
C
	 END IF 
C
      END DO 
C
      NOLEVT=INDEKS-1
      NOLEVT=(N_PART/2)*2
      WRITE(0,'(''N_PART,NOLEVT='',2I3)')N_PART,NOLEVT
      WRITE(66,'(''N_PART,NOLEVT='',2I3)')N_PART,NOLEVT
C
      ENEMAX=ESINGP(NOLEVT)+EPSENE ! Energy limit for summations
C
C=======================================================================
C     Blocking of levels is not taken into account
C=======================================================================
C
      DO LEVELS=1,NDSPEC
         NBLOCK(LEVELS)=0
      END DO
C
      LEVELS=N_PART/2
      NBLOCK(LEVELS)=1
C
C=======================================================================
C     Intitial values for the solutions (to be)
C=======================================================================
C
      DELTA2=DELTA3_INDICT**2
C
      XLAMBD=0.5*(ESINGP(N_PART/2)+ESINGP(N_PART/2+1))
C
      PARTIC=(N_PART/2)*2
C
C=======================================================================
C     Solving the system of Delta^(3) & Particle numbe equations
C=======================================================================
C
      CALL SYSPAI(XLAMBD,DELTA2,PARTIC,DELTA3_INDICT,EPS_XY,EPS_FG,
     *            IERROR,MXITER,ITERAC,DIFF_F,DIFF_G,EQSBC3)
C
      write(0,'(''XLAMBD='',f12.6,'' DELTA2='',f16.12)')
     *   XLAMBD,sqrt(DELTA2)
      IF (IERROR.EQ.1) THEN
          XLAMBD=0.5*(ESINGP(N_PART/2)+ESINGP(N_PART/2+1))
          DELTA2=0.01
      END IF
C
C=======================================================================
C      
      DO LQNUMB=NSHELL,0,-1
C
         DO NNUMBR=0,(NSHELL-LQNUMB)/2
C	
	    EMLAUP=ENERUP(LQNUMB,NNUMBR+1)-XLAMBD
	    EMLADN=ENERDN(LQNUMB,NNUMBR+1)-XLAMBD
C	 
            VCOEUP(LQNUMB,NNUMBR+1)
     *     =
     *	    0.5*(1.-EMLAUP/(SQRT(EMLAUP**2+DELTA2)))
C	 
            IF (LQNUMB.NE.0) THEN
C	 
                VCOEDN(LQNUMB,NNUMBR+1)
     *         =
     *          0.5*(1.-EMLADN/(SQRT(EMLADN**2+DELTA2)))
            ELSE
	        VCOEDN(LQNUMB,NNUMBR+1)=0.0
            END IF
C
         END DO
      END DO     
C
C=======================================================================
C     Testing the v2 coefficients...
C=======================================================================
C
      IF (IERROR.EQ.0) THEN
C
          TESTV2=0.
C
          DO INDEXX=1,NOLEVT
      	     EMLAMB=ESINGP(INDEXX)-XLAMBD	       
             V2COEF=0.5*(1.-EMLAMB/(SQRT(EMLAMB**2+DELTA2)))
	     TESTV2=TESTV2+V2COEF
          END DO
C
          TESTV2=2.*TESTV2
          COMPV2=TESTV2-N_PART
C
          IF (ABS(COMPV2).GT.EPSILO) THEN ! STOP 'V2 coefficients!!!'
              WRITE(0,'(''Test v2='',f30.16)')TESTV2
          END IF
C
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
      SUBROUTINE PAIDEL(ENERUP,ENERDN,NSHELL,ENEMAX,
     *                  VCOEUP,VCOEDN,NDIM_L,NDBASE,INUCLI)
C
      INCLUDE   'MATDIM/NDSPEC.f'
      INCLUDE   'MATDIM/NDIM_N.f'
      INCLUDE   'MATDIM/NDNUCL.f'
C
      EXTERNAL
     *          ONEBCS
C
      DIMENSION
     *          ENERUP(0:NDIM_L,1:NDBASE),
     *          ENERDN(0:NDIM_L,1:NDBASE)
      DIMENSION
     *          VCOEUP(0:NDIM_L,1:NDBASE),
     *          VCOEDN(0:NDIM_L,1:NDBASE)
C
      COMMON
     *       /DIMSPH/ LDSING_PROTON,LEVTHE_PROTON,
     *                LEVEXP_PROTON(1:NDNUCL),
     *      	      LDSING_NEUTRS,LEVTHE_NEUTRS,
     *                LEVEXP_NEUTRS(1:NDNUCL)
      COMMON
     *       /SPECTR/ ENETHE_PROTON(1:NDSPEC),
     *      	      ENETHE_NEUTRS(1:NDSPEC)
      COMMON
     *       /DEGENE/ IDGLEV_PROTON(1:NDSPEC),
     *                ICUMUL_PROTON(0:NDSPEC),
     *                IDGLEV_NEUTRS(1:NDSPEC),
     *                ICUMUL_NEUTRS(0:NDSPEC)
C
      COMMON
     *       /BCSBCS_ENETHE/ ENETHE(1:NDSPEC)
     *       /BCSBCS_IDGLEV/ IDGLEV(1:NDSPEC)
      COMMON
     *       /BCSBCS_BLCKLV/ NBLOCK(1:NDSPEC)
      COMMON
     *       /BCSBCS_RES_WS/ ESINGP(1:NDSPEC)
      COMMON
     *       /BCSBCS_LIMPAI/ NOLEVT,INDEXX_MAXSPH
C
      COMMON
     *       /NUCLEU/ IZ_FIX,IN_FIX,ISOSPI    
      COMMON
     *       /BCSVAL/ XLAMBD,DELTA2
      COMMON
     *       /BCSDEL/ DELT_2(1:NDNUCL,1:2)
      COMMON
     *       /BCSLAM/ XLAMBD_PROTON(1:NDNUCL),
     *                XLAMBD_NEUTRS(1:NDNUCL)
      COMMON
     *       /STEPBI/ BISTEP
      COMMON
     *       /PRFILE/ ISCREN,LOGFIL,LOGCHP,LOGCHN,LOGHIP,LOGHIN,
     *                LOGPRO,LOGNEU,LOGWRI,LSCREN,IRESUL,ICONVE,
     *                NOUTPT,LOGAUX,LOGENE,LOGRAD,LOGGRD,LOGMSH,
     *                                                   LOGBIS
C
      DATA
     *        EPSBIS / 1.0E-14 /,
     *        EPSILO / 5.0E-08 /
      DATA
     *        MXITER /  150 /
C
C=======================================================================     
C     This subroutine calculates the BCS pairing Delta and Lambda
C     parameters as well as the paring v2 occupation coefficients
C=======================================================================
C
      IF (LOGWRI.GT.0) THEN
          WRITE(LOGFIL,'(/,''Entering PAIDEL with ISOSPI= '',I1,2X,
     *                     ''IZ_FIX= '',I3,2X,''and IN_FIX= '',I3)')
     *                                  ISOSPI,IZ_FIX,IN_FIX
      END IF
C
C=======================================================================     
C
      IF (ISOSPI.EQ.1) THEN
C
          N_PART=IZ_FIX
          LDSING=LDSING_PROTON      
C                     LDSING = The number of degenerate 
          DO INDEXX=1,LDSING ! spherical  proton-levels 
	     ENETHE(INDEXX)=ENETHE_PROTON(INDEXX)
             IDGLEV(INDEXX)=IDGLEV_PROTON(INDEXX)
          END DO
C
      END IF
C
      IF (ISOSPI.EQ.0) THEN
C
          N_PART=IN_FIX
          LDSING=LDSING_NEUTRS
C                     LDSING = The number of degenerate 
          DO INDEXX=1,LDSING ! spherical neutron-levels
	     ENETHE(INDEXX)=ENETHE_NEUTRS(INDEXX)
             IDGLEV(INDEXX)=IDGLEV_NEUTRS(INDEXX)
          END DO
C
      END IF
C
C=======================================================================
C
      NOLEVT=0
      INDEKS=0
      EPSENE=1.0E-11
C      
C     We artificially remove the spherical degeneracy in order
C     to use particle-number equation  solving  routine ONEBCS  
C
C     IMPORTANT: we are taking the convention to split the le-
C                vels in (2j+1) "new" levels. With this choice
C                the particle equation to solve is N=\sum v^2   
C 
      DO INDEXX=1,LDSING
C
         IF (INDEKS.LE.2*N_PART) THEN
C
             IDEGEN=IDGLEV(INDEXX)   !/2
             IDEGE2=IDEGEN           !/2
C
	     DO I=0,IDEGEN-1
C
                INDEKS=INDEKS+1
                ESINGP(INDEKS)=ENETHE(INDEXX)+(REAL(I)-IDEGE2)*EPSENE
C
                IF (LOGWRI.GE.5) THEN
                    WRITE(LOGFIL,'('' INDEXX= '',I4,'' I= '',I4,
     *                         '' INDEKS= '',I4,'' ENETHE(INDEXX)= '',
     *                            F16.12,'' ESINGP(INDEKS)= '',F16.12)')
     *                            INDEXX,I,INDEKS,ENETHE(INDEXX),
     *                                            ESINGP(INDEKS)
                END IF
C                
	     END DO
C             
             INDEXX_MAXSPH=INDEXX
C
         END IF 
C
      END DO 
C
      NOLEVT=INDEKS
C
      ENEMAX=ESINGP(NOLEVT)+EPSENE ! Energy limit for summations
C      
      IF (LOGWRI.GT.0) THEN
          WRITE(LOGFIL,'(/,24X,''Before entering BISECT: N_PART= '',I5,
     *                     '' NOLEVT= '',I5,'' ENEMAX= '',F20.9,/)')
     *                        N_PART,NOLEVT,ENEMAX
      END IF
C
C=======================================================================
C     Blocking of levels is not taken into account
C=======================================================================
C
      DO LEVELS=1,NDSPEC
         NBLOCK(LEVELS)=0
      END DO
C
C=======================================================================
C     Preparing empirical Delta's for one BCS equation run
C=======================================================================
C
CID      DELTAP=DELEXP(IZ_FIX,IN_FIX,ISOSPI)
      IF (ISOSPI.EQ.1) DELTA2=DELT_2(INUCLI,1)
      IF (ISOSPI.EQ.0) DELTA2=DELT_2(INUCLI,2)
          
      DELTAP=SQRT(DELTA2)
C
CID      XLAMBD=0.5*(ESINGP(N_PART/2)+ESINGP(N_PART/2+1))
      XLAMBD=0.5*(ESINGP(N_PART)+ESINGP(N_PART+1))
C
      ALIMIT=XLAMBD-BISTEP
      BLIMIT=XLAMBD+BISTEP
C
C
C=======================================================================
C
      IF (LOGWRI.GT.0) THEN
          WRITE(LOGFIL,'(24X,''DELTA2= '',F20.13,'' XLAMBD= '',F20.13)')
     *                                          DELTA2,XLAMBD
          WRITE(LOGFIL,'(24X,''ALIMIT= '',F20.13,'' BLIMIT= '',F20.13)')
     *                                          ALIMIT,BLIMIT
      END IF
C
C=======================================================================
C
C     Defining MXITER from EPSBIS using the mathematical theorem
C
      ABS_AB=ABS(ALIMIT-BLIMIT)
    
      MXITER=NINT(LOG(ABS_AB/EPSBIS)/LOG(2.0))
C
      IF (LOGWRI.GT.0) THEN
          WRITE(LOGFIL,'(/,24X,''MXITER= '',I3)')MXITER
      END IF
C
C=======================================================================
C     Solving particle number BCS-equation at fixed Delta
C=======================================================================
C
      CALL BISECT(ALIMIT,BLIMIT,XLAMB1,ONEBCS,EPSBIS,MXITER,
     *                                 BISTEP,IT_ACT,IERBIS)
C
      IF (IERBIS.NE.0) THEN
C
          XLAMBD=0.5*(ESINGP(N_PART)+ESINGP(N_PART+1))
C
          WRITE(LOGFIL,'(''BISECT: IERBIS='',I3,'' MXITER='',I3,
     *              '' IT_ACT= '',i3)')IERBIS,MXITER,IT_ACT
          WRITE(LOGFIL,'(''ALIMIT= '',F20.13,'' BLIMIT= '',F20.13,
     *                   '' XLAMBD= '',F20.13)') ALIMIT,BLIMIT,XLAMBD
     
          ONEBCS_FUNCTN=ONEBCS(XLAMB1)
          
          WRITE(LOGFIL,'(''AFTER BISECT: ONEBCS_FUNCTN= '',F20.13)')
     *                                   ONEBCS_FUNCTN
          
          WRITE(0,'(''BISECT: IERBIS='',I3,'' MXITER='',I3,
     *              '' IT_ACT= '',i3)')IERBIS,MXITER,IT_ACT
          STOP 'STOP from PAIDEL: IERBIS.NE.0 !!'
      ELSE
      
          XLAMBD=XLAMB1
          IF (LOGWRI.GT.0) THEN
              WRITE(LOGFIL,'(/,24X,''ONEBCS Solution: DELTAP= '',f20.13,
     *               '' XLAMBD= '',F20.13,/)')DELTAP,XLAMBD
          END IF
     
          IF (ISOSPI.EQ.1) XLAMBD_PROTON(INUCLI)=XLAMBD
          IF (ISOSPI.EQ.0) XLAMBD_NEUTRS(INUCLI)=XLAMBD
          
      END IF
C
C=======================================================================
C      
      DO LQNUMB=NSHELL,0,-1
C
         DO NNUMBR=0,(NSHELL-LQNUMB)/2
C	
	    EMLAUP=ENERUP(LQNUMB,NNUMBR+1)-XLAMBD
	    EMLADN=ENERDN(LQNUMB,NNUMBR+1)-XLAMBD
C	 
            VCOEUP(LQNUMB,NNUMBR+1)
     *     =
     *	    0.5*(1.-EMLAUP/(SQRT(EMLAUP**2+DELTA2)))
C	 
            IF (LQNUMB.NE.0) THEN
C	 
                VCOEDN(LQNUMB,NNUMBR+1)
     *         =
     *          0.5*(1.-EMLADN/(SQRT(EMLADN**2+DELTA2)))
            ELSE
	        VCOEDN(LQNUMB,NNUMBR+1)=0.0
            END IF
C
         END DO
      END DO     
C
C=======================================================================
C     Testing the v2 coefficients...
C=======================================================================
C
      TESTV2=0.
C
      DO INDEXX=1,NOLEVT
         EMLAMB=ESINGP(INDEXX)-XLAMBD	       
         V2COEF=0.5*(1.-EMLAMB/(SQRT(EMLAMB**2+DELTA2)))
	 TESTV2=TESTV2+V2COEF
      END DO
C
      TESTV2=TESTV2    !*2
      COMPV2=TESTV2-N_PART
C
      IF (ABS(COMPV2).GT.EPSILO) THEN ! 
          WRITE(LOGFIL,'(''Test ONE-EQ v2='',f30.16)')TESTV2
          STOP 'STOP from PAIDEL: V2 coefficients!!!'
      END IF
C
C=======================================================================
C
      IF (LOGWRI.GT.0) THEN
          WRITE(LOGFIL,'(/,''Exiting PAIDEL'')')
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
      FUNCTION ONEBCS(XLAMBD)
C
      INCLUDE   'MATDIM/NDSPEC.f'
C
      COMMON
     *       /BCSBCS_RES_WS/ ESINGP(1:NDSPEC)
      COMMON
     *       /BCSBCS_LIMPAI/ NOLEVT,INDEXX_MAXSPH
C
      COMMON
     *       /NUCLEU/ IZ_FIX,IN_FIX,ISOSPI    
      COMMON
     *       /BCSVAL/ XDUMMY,DELTA2
C
C=======================================================================
C     Preparation for solving the particle-number equation 
C                         N-sum_j v_j^2 = 0
C     for Lambda 
C=======================================================================
C
      PARTEQ=0.
                       N_PART=IZ_FIX
      IF (ISOSPI.EQ.0) N_PART=IN_FIX
C
      DO INDEXX=1,NOLEVT
      	 EMLAMB=ESINGP(INDEXX)-XLAMBD	       
         V2COEF=0.5*(1.-EMLAMB/(SQRT(EMLAMB**2+DELTA2)))
	 PARTEQ=PARTEQ+V2COEF
      END DO
C
      PARTEQ=PARTEQ         !*2
      ONEBCS=PARTEQ-N_PART
C
C=======================================================================     
C
      RETURN
      END      
C
C=======================================================================
C=======================================================================
C
      FUNCTION DELEXP(IZ_FIX,IN_FIX,ISOSPI)
      PARAMETER
     *          (NDDATA=8)
      DIMENSION
     *           I_PROT(1:NDDATA),
     *           I_NEUT(1:NDDATA)
      DIMENSION
     *           DELPRO(1:NDDATA),
     *           DELNEU(1:NDDATA)
C
C=======================================================================
C
C     This function provides the results of the Delta^3 filter 
C     calculations, based on the mass tables of Wapstra et al.
C     For the meaning of this type of filters that approximate
C     the experimental pairing gaps cf.:
C
C     arXiv:nucl-th/0003019 
C
C     Title: Odd-even staggering of binding energies as 
C            a consequence of pairing and mean-field effects
C
C     Authors: J. Dobaczewski, P. Magierski, W. Nazarewicz, 
C              W. Satula, Z. Szymanski
C
C     Journal-ref: Phys.Rev. C63 (2001) 024308
C                                                          AND
C     arXiv:nucl-th/9804060
C
C     Title: Odd-Even Staggering of Nuclear Masses: Pairing or 
C            Shape Effect?
C
C     Authors: W. Satula, J. Dobaczewski, W. Nazarewicz
C     Journal-ref: Phys.Rev.Lett.81:3599-3602,1998
C
C=======================================================================
C
      I_PROT(1)=8      ! 16O
      I_NEUT(1)=8      ! 16O
C
      DELPRO(1)=1.6616 ! Delta in MeV
      DELNEU(1)=1.9505
C
      I_PROT(2)=20     ! 40Ca
      I_NEUT(2)=20     ! 40Ca
C
      DELPRO(2)=1.3415 ! Delta in MeV
      DELNEU(2)=1.5590
C
      I_PROT(3)=20     ! 48Ca
      I_NEUT(3)=28     ! 48Ca
C
      DELPRO(3)=1.2680 ! Delta in MeV
      DELNEU(3)=0.6031
C
      I_PROT(4)=28     ! 56Ni
      I_NEUT(4)=28     ! 56Ni
C
      DELPRO(4)=0.7911 ! Delta in MeV
      DELNEU(4)=0.9859
C
      I_PROT(5)=40     ! 90Zr
      I_NEUT(5)=50     ! 90Zr
C
      DELPRO(5)=1.1537 ! Delta in MeV
      DELNEU(5)=0.7201 
C
      I_PROT(6)=50     ! 132Sn
      I_NEUT(6)=82     ! 132Sn
C
      DELPRO(6)=0.6143 ! Delta in MeV
      DELNEU(6)=0.6614 
C
      I_PROT(7)=64     ! 146Gd
      I_NEUT(7)=82     ! 146Gd      
C
      DELPRO(7)=1.2243 ! Delta in MeV
      DELNEU(7)=0.8215 
C
      I_PROT(8)=82     ! 208Pb
      I_NEUT(8)=126    ! 208Pb      
C
      DELPRO(8)=0.5924 ! Delta in MeV
      DELNEU(8)=0.6244 
C
C=======================================================================
C
      DO I=1,NDDATA
         IF (I_PROT(I).EQ.IZ_FIX.AND.I_NEUT(I).EQ.IN_FIX) THEN
             IF (ISOSPI.EQ.1) DELEXP=DELPRO(I)
             IF (ISOSPI.EQ.0) DELEXP=DELNEU(I)
         END IF
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
      SUBROUTINE BISECT(ALIMIT,BLIMIT,XARGUM,FUNBIS,EPSBIS,ITEMAX,
     *                                       BISTEP,IT_ACT,IERBIS)
      EXTERNAL
     *         FUNBIS
      COMMON
     *        /PRFILE/ ISCREN,LOGFIL,LOGCHP,LOGCHN,LOGHIP,LOGHIN,
     *                 LOGPRO,LOGNEU,LOGWRI,LSCREN,IRESUL,ICONVE,
     *                 NOUTPT,LOGAUX,LOGENE,LOGRAD,LOGGRD,LOGMSH,
     *                                                    LOGBIS
C
C=======================================================================
C
C     BISECTION METHOD FOR THE SOLUTION OF NONLINEAR EQUATIONS
C
C=======================================================================
C
C     ALIMIT = Lower limit for the search interval (input)
C     BLIMIT = Upper limit for the search interval (input)
C
C     XARGUM = On exit ===>>>> the solution of the problem
C
C     EPSBIS = Solution tolerance: if the function studied
C              satisfies F(x)< EPSBIS iterations will stop
C
C     ITEMAX = Maximum number of iterations allowed
C
C     FUNBIS = External function name; this FUNBIS defines
C              the function to be studied
C
C=======================================================================
C
C     IERBIS = Error index; If solution found  IERBIS=0
C
C                           If the function does not 
C                           change sign in the original
C                           interval          IERBIS=10
C
C                           If the number of iterations
C                           exceeded         IERBIS=100
C
C=======================================================================
C
      XA_LIM=ALIMIT
      XB_LIM=BLIMIT
C
      IT_ACT=0
C
C=======================================================================
C
      FA_BIS=FUNBIS(XA_LIM)
C
      FB_BIS=FUNBIS(XB_LIM)
C
C=======================================================================
C
      IF (LOGWRI.GE.5) THEN
C
          WRITE(LOGBIS,'(''Routine Parameters: EPSBIS= '',E10.2,
     *                   '' ITEMAX= '',I4,'' BISTEP= '',F10.4,/)')
     *                                       EPSBIS,ITEMAX,BISTEP
          
          WRITE(LOGBIS,'(''Initial guess: XA_LIM= '',F20.13,
     *                                 '' FA_BIS= '',F20.13,/,14X,
     *                                 '' XB_LIM= '',F20.13,
     *                                 '' FB_BIS= '',F20.13,/)') XA_LIM,
     *                                             FA_BIS,XB_LIM,FB_BIS
      END IF
C
C=======================================================================
C     Verifying whether the solution is not lying on the interval
C                                                          limits
C=======================================================================
C
      IF (LOGWRI.GE.5) THEN
          WRITE(LOGBIS,'(''Verifying whether the solution is not '',
     *                   ''lying on the interval limits...'')')
      END IF
      
      ABS_FA=ABS(FA_BIS)
      XARGUM=ALIMIT
C
      IF (ABS_FA.LT.EPSBIS) THEN
          IERBIS=0
          IF (LOGWRI.GE.5) THEN
              WRITE(LOGBIS,'(''ABS_FA= '',F20.13,'' is less than '',
     *                       ''EPSBIS= '',E10.2,'' ==> IERBIS=0 ==>'',
     *                       '' the solution is XARGUM= '',F20.13,/)')
     *                         ABS_FA,EPSBIS,XARGUM
              WRITE(LOGBIS,'(''RETURN'',/,/,''Exiting BISECT'')')
          END IF
          RETURN
      END IF
C
      ABS_FB=ABS(FB_BIS)
      XARGUM=BLIMIT
C
      IF (ABS_FB.LT.EPSBIS) THEN
          IERBIS=0
          IF (LOGWRI.GE.5) THEN
              WRITE(LOGBIS,'(''ABS_FB= '',F20.13,'' is less than '',
     *                       ''EPSBIS= '',E10.2,'' ==> IERROR=0 ==>'',
     *                       '' the solution is XARGUM= '',F20.13,/)')
     *                         ABS_FB,EPSBIS,XARGUM
              WRITE(LOGBIS,'(''RETURN'',/,/,''Exiting BISECT'')')
          END IF
          RETURN
      END IF
      
      IF (LOGWRI.GE.5) THEN
          WRITE(LOGFIL,'(''No solution in the interval limits. '',
     *                   ''We continue...'',/)')
      END IF
C
C=======================================================================
C         The two IF's below verify whether we are still within
C         the correcct interval (function changes sign) and we
C         constrain the interval 
C=======================================================================
C
      IF (LOGWRI.GE.5) THEN
          WRITE(LOGBIS,'(''Verifying whether we are in the correct '',
     *                   ''interval...'',/)')
      END IF
C      
      IF (FA_BIS*FB_BIS.GT.0.0) THEN
C          
          IF (LOGWRI.GE.5) THEN 
              WRITE(LOGBIS,'(''1st Time: FA_BIS*FB_BIS.GT.0.0'',/)')
          END IF
C
          IF (FB_BIS.GT.0.0) THEN
C
              XB_LIM=XB_LIM-BISTEP
C
              FB_BIS=FUNBIS(XB_LIM)
C
              IF (LOGWRI.GE.5) THEN
                  WRITE(LOGFIL,'(''FB_BIS.GT.0.0: new '',
     *                           ''XB_NEW=XB_OLD-BISTEP= '',F20.13,
     *                           '' FB_BIS= '',F20.13,/)')XB_LIM,FB_BIS
              END IF
C
          ELSE
C
              XB_LIM=XB_LIM+BISTEP
C
              FB_BIS=FUNBIS(XB_LIM)
C              
              IF (LOGWRI.GE.5) THEN
                  WRITE(LOGFIL,'(''FB_BIS.LT.0.0: new '',
     *                           ''XB_NEW=XB_OLD+BISTEP= '',F20.13,
     *                           '' FB_BIS= '',F20.13,/)')XB_LIM,FB_BIS
              END IF
C
          END IF
C
      END IF
C
C=======================================================================
C
      IF (FA_BIS*FB_BIS.GT.0.0) THEN
C
          IF (LOGWRI.GE.5) THEN 
              WRITE(LOGBIS,'(''2nd Time: FA_BIS*FB_BIS.GT.0.0'',/)')
          END IF
C
          IF (FB_BIS.GT.0.0) THEN
C
              XB_LIM=XB_LIM-2.0*BISTEP
C
              FB_BIS=FUNBIS(XB_LIM)
C              
              IF (LOGWRI.GE.5) THEN
                  WRITE(LOGFIL,'(''FB_BIS.GT.0.0: new '',
     *                           ''XB_NEW=XB_OLD-2.0*BISTEP= '',F20.13,
     *                           '' FB_BIS= '',F20.13,/)')XB_LIM,FB_BIS
              END IF
C
          ELSE
C
              XB_LIM=XB_LIM+2.0*BISTEP
C
              FB_BIS=FUNBIS(XB_LIM)
C              
              IF (LOGWRI.GE.5) THEN
                  WRITE(LOGFIL,'(''FB_BIS.GT.0.0: new '',
     *                           ''XB_NEW=XB_OLD+2.0*BISTEP= '',F20.13,
     *                           '' FB_BIS= '',F20.13,/)')XB_LIM,FB_BIS
              END IF
C
          END IF
C
      END IF
C
C=======================================================================
C
      IF (FA_BIS*FB_BIS.GT.0.0) THEN
C
          IERBIS=10
C
          WRITE(NOUTPT,'(/,80(1H=),/,30X,
     *              ''BISECT OUT OF RANGE: '',/,
     *              ''XA_LIM='',E11.4,1X, ''XB_LIM='',E11.4,1X,
     *              '' F(XA_LIM)='',E10.3,'' F(XB_LIM)='',E10.3,/,
     *                                                80(1H=),/)')
     *
     *                                  XA_LIM,XB_LIM,FA_BIS,FB_BIS
          
          IF (LOGWRI.GE.5) THEN
C          
              WRITE(LOGBIS,'(''IERBIS= '',I2)')IERBIS
              WRITE(LOGBIS,'(/,80(1H=),/,30X,
     *              ''BISECT OUT OF RANGE: '',/,
     *              ''XA_LIM='',E11.4,1X, ''XB_LIM='',E11.4,1X,
     *              '' F(XA_LIM)='',E10.3,'' F(XB_LIM)='',E10.3,/,
     *                                                80(1H=),/)')
     *
     *                                  XA_LIM,XB_LIM,FA_BIS,FB_BIS
C              
              WRITE(LOGBIS,'(''RETURN'',/,/,''Exiting BISECT'')')
C          
          END IF
C          
          RETURN
C
C=======================================================================
C
      END IF
C
C=======================================================================
C
      IF (LOGWRI.GE.5) THEN
          WRITE(LOGBIS,'(/,''...everything seems ok, we begin with '',
     *                   ''the real algorithm'',/,/)')
      END IF
C
C=======================================================================
C                   Bisection algorithm begins here ...
C=======================================================================
C
   1  CONTINUE
C
      IT_ACT=IT_ACT+1
C
C      IF (ABS(XB_LIM-XA_LIM).LT.1.0E-14) THEN
C
C          WRITE(NOUTPT,'(1H*,78(1H=),1H*,/,1H*,
C     *              19X,''FURTHER ITERATIONS BELOW 8 BYTE ACCURACY'',
C     *              19X,1H*,/,1H*,78(1H=),1H*)')
C          FX_BIS=0.
C          GO TO  8
C
C      END IF
C
C=======================================================================
C
      XARGUM=0.5*(XA_LIM+XB_LIM)
C
      FX_BIS=FUNBIS(XARGUM)
C
C   8  CONTINUE
C
      ABS_FX=ABS(FX_BIS)
C      
      IF (LOGWRI.GE.5) THEN
          WRITE(LOGBIS,'(''IT_ACT= '',I4,
     *                   '' XA_LIM= '',F20.15,'' XB_LIM= '',F20.15,/)')
     *                                         IT_ACT,XA_LIM,XB_LIM
C          
          WRITE(LOGBIS,'(13X,''ABS_FX= '',F40.35,
     *                           '' XARGUM= '',F40.35,/,
     *                        13X,''EPSBIS= '',F40.35,/)')ABS_FX,
     *                                                    XARGUM,EPSBIS 
      END IF 
C
      IF (ABS_FX.LT.EPSBIS) THEN
          IERBIS=0
          IF (LOGWRI.GE.5) THEN
              WRITE(LOGBIS,'(''ABS_FX= '',F20.15,'' is less than '',/,
     *                       ''EPSBIS= '',F20.15,''  ==> IERBIS=0'')')
     *                         ABS_FX,EPSBIS
              WRITE(LOGBIS,'(/,''===>>> XARGUM= '',F20.15,/)')XARGUM
              WRITE(LOGBIS,'(''RETURN'',/,/,''Exiting BISECT'')')
          END IF
          RETURN
      END IF
C
      IF (IT_ACT.GT.ITEMAX) THEN
          IERBIS=0  !100
          IF (LOGWRI.GE.5) THEN
              WRITE(LOGBIS,'(''IT_ACT= '',I4,'' is bigger than '',
     *                       ''ITEMAX= '',i4,'' ==> IERBIS=100'')')
     *                         IT_ACT,ITEMAX
              WRITE(LOGBIS,'(''RETURN'',/,/,''Exiting BISECT'')')
          END IF
          RETURN
      END IF
C
C=======================================================================
C         Function changes sign in the left subinterval
C=======================================================================
C
      IF (FA_BIS*FX_BIS.LT.0.0) THEN
C
          XB_LIM=XARGUM
          
          IF (LOGWRI.GE.5) THEN
              WRITE(LOGBIS,'(''Function changes sign in the '',
     *                       ''left subinterval: FA_BIS*FX_BIS.LT.0.0'',
     *                       '': new XB_LIM= '',F20.13)') XB_LIM
              WRITE(LOGBIS,'(''...we go to the next iteration...'',/)')
          END IF
C
          GO TO 1
C
      END IF
C
C=======================================================================
C         Function changes sign in the right subinterval
C=======================================================================
C
      IF (FA_BIS*FX_BIS.GT.0.0) THEN
C
          XA_LIM=XARGUM
          
          IF (LOGWRI.GE.5) THEN
              WRITE(LOGBIS,'(''Function changes sign in the '',
     *                       ''left subinterval: FA_BIS*FX_BIS.GT.0.0'',
     *                       '': new XA_LIM= '',F20.13)') XA_LIM
              WRITE(LOGBIS,'(''...we go to the next iteration...'',/)')
          END IF
C
          GO TO 1
C
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
      FUNCTION G_FITT(INCALC,IZCALC,ISOSPI)
C
C=======================================================================
C     P A I R I N G    P A R A M E T R I Z A T I O N    (empirical)
C=======================================================================
C
      A_MASS=INCALC+IZCALC
C
      IF (ISOSPI.EQ.0) THEN
C
          IF (IZCALC.GE.88) THEN
C
              G_FITT=(19.3-0.084*(INCALC-IZCALC))/A_MASS
C
          ELSE
C
              G_FITT=(18.95-0.078*(INCALC-IZCALC))/A_MASS
C
          END IF
C
      ELSE
C
          IF (IZCALC.GE.88) THEN
C
              G_FITT=(13.3+0.217*(INCALC-IZCALC))/A_MASS
C
          ELSE
C
              G_FITT=(17.90+0.176*(INCALC-IZCALC))/A_MASS
C
          END IF
C
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

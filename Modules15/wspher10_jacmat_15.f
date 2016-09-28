C FILE NAME = wspher10_jacmat_15.f ! Keep this symbol:    $ident@string$
C
C=======================================================================
C=======================================================================
C             JACOBIAN AND PARTIAL DERIVATIVES FOR MINIMISATION
C=======================================================================
C=======================================================================
C
      SUBROUTINE JMATRX(LDFUNC,ARGPAR)
C      
      INCLUDE   'MATDIM/NDSPEC.f'
      INCLUDE   'MATDIM/NDPARS.f'
      INCLUDE   'MATDIM/NDNUCL.f'
      INCLUDE   'MATDIM/NDLEXP.f'
      INCLUDE   'MATDIM/NDIM_M.f'
      INCLUDE   'MATDIM/NDFUNC.f'
      INCLUDE   'MATDIM/ND_RHO.f'
      INCLUDE   'MATDIM/N_NOYX.f'
C
      CHARACTER
     *          LABTHE_PROTON*6,LABTHE_NEUTRS*6
      CHARACTER
     *          LABEXP_PROTON*6,LABEXP_NEUTRS*6
      CHARACTER
     *          LABPRO_PRINTG*8,LABNEU_PRINTG*8,
     *          LABTOT_PRINTG*8
      CHARACTER
     *          ACTION*1,CALLED*6
      CHARACTER
     *          JOBU*1,JOBVT*1
      CHARACTER
     *          NUCSYM*6,INPSYM*6,TYPCHI*6,SYMBNU*5
C
      DIMENSION
     *          ARGPAR(1:NDPARS),DERIVF(1:NDPARS)
      DIMENSION
     *          IDPROT(1:NDIM_M),IDNEUT(1:NDIM_M)
      DIMENSION
     *          DEGENE(1:NDIM_M)
      DIMENSION  
     *                 XJMATR(1:NDIM_M,1:NDPARS),
     *          XJMATR_PROTON(1:NDIM_M,1:NDPARS),
     *          XJMATR_NEUTRS(1:NDIM_M,1:NDPARS)
      DIMENSION  
     *                 CONVER(1:NDPARS,1:NDPARS),
     *                 CONVOL(1:NDPARS,1:NDPARS),
     *                 YJTJMX(1:NDPARS,1:NDPARS),
     *          CONVER_PROTON(1:NDPARS,1:NDPARS),
     *          CONVOL_PROTON(1:NDPARS,1:NDPARS),
     *          YJTJMX_PROTON(1:NDPARS,1:NDPARS),
     *          CONVER_NEUTRS(1:NDPARS,1:NDPARS),
     *          CONVOL_NEUTRS(1:NDPARS,1:NDPARS),
     *          YJTJMX_NEUTRS(1:NDPARS,1:NDPARS),
     *                   TEST(1:NDPARS,1:NDPARS),
     *            TEST_PROTON(1:NDPARS,1:NDPARS),
     *            TEST_NEUTRS(1:NDPARS,1:NDPARS)
      DIMENSION
     *          SINGLE_PROTON(1:NDPARS),
     *          UMATRI_PROTON(1:NDIM_M,1:NDIM_M),
     *          VTMATR_PROTON(1:NDPARS,1:NDPARS),
     *          SINGLE_NEUTRS(1:NDPARS),
     *          UMATRI_NEUTRS(1:NDIM_M,1:NDIM_M),
     *          VTMATR_NEUTRS(1:NDPARS,1:NDPARS),
     *                 SINGLE(1:NDPARS),
     *                 UMATRI(1:NDIM_M,1:NDIM_M),
     *                 VTMATR(1:NDPARS,1:NDPARS)
      DIMENSION  
     *          DERIVV(1:NDPARS,1:NDIM_M)
      DIMENSION
     *          FUNCT0(1:NDIM_M),
     *          FUNCT1(1:NDIM_M),
     *          LEVELS(1:NDIM_M),
     *          LEXPLS(1:NDIM_M),
     *          WEISQR_PROTON(1:NDIM_M),
     *          WEISQR_NEUTRS(1:NDIM_M),
     *          WEISQR(1:NDIM_M),
     *          FUNEXP(1:NDIM_M)
      DIMENSION
     *          PARNEW_PROTON(1:NDPARS),
     *          PARNEW_NEUTRS(1:NDPARS),
     *                 PARNEW(1:NDPARS)
      DIMENSION
     *          LABPRO_PRINTG(1:NDIM_M),LABNEU_PRINTG(1:NDIM_M),
     *          LABTOT_PRINTG(1:NDIM_M)
      DIMENSION
     *          WORK(1:NDSPEC),IPIV(1:NDPARS)
      DIMENSION
     *          CHISQU_GRADNT(1:NDPARS),
     *          FUNAUX(1:NDIM_M)
C     
      COMMON
     *       /DIMSPH/ LDSING_PROTON,LEVTHE_PROTON,
     *                LEVEXP_PROTON(1:NDNUCL),
     *                LDSING_NEUTRS,LEVTHE_NEUTRS,
     *                LEVEXP_NEUTRS(1:NDNUCL)
      COMMON
     *       /CHIVAL_PROTON/ CHICOR_PROTON,RADDIF_PROTON,CHIINV_PROTON,
     *                       FERDIF_PROTON,GAPDIF_PROTON,CHIWEI_PROTON,
     *                       EABSWD_PROTON,EABSAV_PROTON,ERRMAX_PROTON,
     *                       DIFFUP_PROTON,DIFFDW_PROTON,CHIRHO_PROTON,
     *                                                   CHIDEG_PROTON
      COMMON
     *       /CHIVAL_NEUTRS/ CHICOR_NEUTRS,RADDIF_NEUTRS,CHIINV_NEUTRS,
     *                       FERDIF_NEUTRS,GAPDIF_NEUTRS,CHIWEI_NEUTRS,
     *                       EABSWD_NEUTRS,EABSAV_NEUTRS,ERRMAX_NEUTRS,
     *                       DIFFUP_NEUTRS,DIFFDW_NEUTRS,CHIDEG_NEUTRS
      COMMON
     *       /NUCWEI/ WEINUC_PROTON(1:NDNUCL),
     *                WEINUC_NEUTRS(1:NDNUCL),
     *                WEIRAD_PROTON(1:NDNUCL),
     *                WEIRAD_NEUTRS(1:NDNUCL)
      COMMON
     *       /WEIGHT/ WEIGHT_CORREL,WEIGHT_RADIUS,WEIGHT_INVERT,
     *                WEIGHT_EFERMI,WEIGHT_ENEGAP,WEIWEI,
     *                WEIGHT_ERRABS,WEIGHT_EABSAV,WEIGHT_ERRMAX,
     *                WEIGHT_DENSUP,WEIGHT_DENSDW,WEIGHT_RHODEN
      COMMON
     *       /BETDAT/ APARAM,BPARAM,PROBAB,FREDEG
      COMMON
     *       /MODIFR/ RMSEXP_PROTON(1:NDNUCL),RMSEXP_NEUTRS(1:NDNUCL),
     *                RMSTHE_PROTON,RMSTHE_NEUTRS,
     *                RMODIF_PROTON,RMODIF_NEUTRS,
     *                              PUSHIN,DELTAR
       COMMON
     *       /THEORF/ FERTHE_PROTON,GAPTHE_PROTON,
     *                DENUPP_PROTON,DENLOW_PROTON,
     *
     *                FERTHE_NEUTRS,GAPTHE_NEUTRS,
     *                DENUPP_NEUTRS,DENLOW_NEUTRS 
      COMMON
     *       /FERGPS/ FERMEX_PROTON(1:NDNUCL),GAPEXP_PROTON(1:NDNUCL),
     *                DENSUP_PROTON(1:NDNUCL),DENSDW_PROTON(1:NDNUCL),
     *
     *                FERMEX_NEUTRS(1:NDNUCL),GAPEXP_NEUTRS(1:NDNUCL),
     *                DENSUP_NEUTRS(1:NDNUCL),DENSDW_NEUTRS(1:NDNUCL) 
      COMMON
     *       /EX_LAB/ LABEXP_PROTON(1:NDNUCL,1:NDLEXP),
     *                LABEXP_NEUTRS(1:NDNUCL,1:NDLEXP)
      COMMON
     *       /EX_ENE/ EXPEXP_PROTON(1:NDNUCL,1:NDLEXP),
     *                IDEGEX_PROTON(1:NDNUCL,1:NDLEXP),
     *
     *                EXPEXP_NEUTRS(1:NDNUCL,1:NDLEXP),
     *                IDEGEX_NEUTRS(1:NDNUCL,1:NDLEXP)
      COMMON
     *       /SPECTR/ ENETHE_PROTON(1:NDSPEC),
     *                ENETHE_NEUTRS(1:NDSPEC)
      COMMON
     *       /THELAB/ LABTHE_PROTON(1:NDSPEC),
     *                LABTHE_NEUTRS(1:NDSPEC)
      COMMON
     *       /SELNUC/ NUCACT,
     *                ITAKNU(1:NDNUCL),
     *                INPUTZ(1:NDNUCL),
     *                INPUTN(1:NDNUCL),
     *                INPSYM(1:NDNUCL)
      COMMON
     *       /NUCINF/ LDNUCL,
     *                NUCSYM(1:NDNUCL),
     *                NUMB_Z(1:NDNUCL),
     *                NUMB_N(1:NDNUCL)
      COMMON
     *       /MESHIN/ I_MESH(1:NDPARS),
     *                XMIN_I(1:NDPARS),
     *                XMAX_I(1:NDPARS),
     *                MAXPAR(1:NDPARS)
      COMMON
     *       /TESFIT/ IFTEST,IFFITS,IFMESH,ISPE_P,IFMOCA
      COMMON
     *       /EXPDEN/ LABELN(1:N_NOYX,1:2),
     *                R_MESH(1:ND_RHO,1:N_NOYX),
     *                RHOEXP(1:ND_RHO,1:N_NOYX)
     *
     *       /EXPSYM/ SYMBNU(1:N_NOYX)
      COMMON
     *       /THEDEN/ N_NUCL,I_NUCL(1:N_NOYX),
     *                RHOTHE(1:ND_RHO,1:N_NOYX)
      COMMON
     *       /SEARCH/ IFTAKE(1:NDPARS)
      COMMON
     *       /ACTACT/ ACTION(1:NDPARS)
      COMMON
     *       /CNTROL/ PARAMT(1:NDPARS),
     *                DPARAM(1:NDPARS)     
      COMMON
     *       /POTPOT/ PARPOT(1:NDPARS)
      COMMON
     *       /ARGINI/ XINITS(1:NDPARS)
      COMMON
     *       /MINITE/ IDEFCN,ITECHI 
      COMMON
     *       /CHI_PN/ CHISQU_PROTON,CHISQU_NEUTRS
      COMMON
     *       /CHICHO/ VALCHI
      COMMON
     *       /CHIGRD/ CHISQU_GRDNRM,
     *                CHISQU_GRADIE(1:NDPARS)
      COMMON
     *       /DENSTR/ IFDENS
CID      COMMON
CID     *       /SVDSVC/ SVDCUT
      COMMON
     *       /P_N_LM/ IFPROT,IFNEUT,IFBOTH
      COMMON 
     *       /CHITAK/ IF_SPE,IF_RAD,IF_GAP,IF_FER,IF_DEN,IF_RHO,IF_INV,
     *                LDCHI2,TYPCHI(1:NDNUCL),ITKCHI(1:NDNUCL)
      COMMON
     *       /PRFILE/ ISCREN,LOGFIL,LOGCHP,LOGCHN,LOGHIP,LOGHIN,
     *                LOGPRO,LOGNEU,LOGWRI,LSCREN,IRESUL,ICONVE,
     *                NOUTPT,LOGAUX,LOGENE,LOGRAD,LOGGRD,LOGMSH,
     *                                                   LOGBIS 
C
      DATA    EPSTES /1.0E-7/
C
C=======================================================================
C
C     This subroutine calculates first the JACOBIAN matrix of the form
C
C                           J(k,j)=dF(k)/dp(j) 
C
C     where F denotes all the calculated theoretical quantities  which 
C     are compared with experimental values, mainly sp levels but also 
C     r.m.s. radii and optionally, Fermi level, energy gaps and others 
C     All the parameters which are used for the minimisation are deno-
C     det "p" for short.
C @@@
C     The J matrix gives the covariance matrix (on output --> CONVER):
C
C                            [C]=[J^T J]^(-1)
C
C     The choice of the physical quantities for the Jacobian depends 
C     on the definition of the chi^2 function to be minimised.
C
C=======================================================================
C     
CID      CALL CPUTIM('JMATRX',1)
C
C=======================================================================
C      
      If (ISCREN.GT.0) WRITE(LSCREN,'(''entering JMATRX...'')')
C
      IF (LOGWRI.GT.0) THEN
          WRITE(LOGFIL,'(12X,''entering JMATRX...'')')
          WRITE(LOGAUX,'(''entering JMATRX...'')')
      END IF
C
C=======================================================================
C
      LDSEAR=0
      IACTIV=0
C      
      DO IPARAM=1,NDPARS
         ACTION(IPARAM)=' '
         IF (IFTAKE(IPARAM).EQ.1) THEN
             IACTIV=IACTIV+1
             PARPOT(IPARAM)=ARGPAR(IACTIV)
         END IF
C
         IF (LOGWRI.GT.4) THEN
             WRITE(LOGFIL,'(20X,I3,F10.4)')IPARAM,PARPOT(IPARAM)
         END IF
C         
         IF (IFMESH.EQ.1) THEN
             IF (I_MESH(IPARAM).EQ.1) THEN
                 PARPOT(IPARAM)=XINITS(IPARAM)
             END IF
         END IF
C
      END DO
C
      LDSEAR=IACTIV
C_______________________________________________________________________
C
C     Initialising auxiliary variables
C      
      DO I_NDEX=1,NDIM_M
         FUNCT0(I_NDEX)=0.
         FUNCT1(I_NDEX)=0.
         LEVELS(I_NDEX)=0
         LEXPLS(I_NDEX)=0
	 WEISQR(I_NDEX)=0.
      END DO
C      
      DO IPARAM=1,NDPARS
         DO INDEXX=1,NDFUNC
            XJMATR(INDEXX,IPARAM)=0.
            XJMATR_PROTON(INDEXX,IPARAM)=0.
            XJMATR_NEUTRS(INDEXX,IPARAM)=0.
         END DO
      END DO
C      
      CHIENE_WEIPRO=0.0
      CHIENE_WEINEU=0.0
      CHIRAD_WEIPRO=0.0
      CHIRAD_WEINEU=0.0
      CHIRHO_WEIPRO=0.0
C      
      ILEVEL=0
      LEVPRO=0
      LEVNEU=0
      
      ILEVEL_PROTON=0
      ILEVEL_NEUTRS=0
       
      INITI1=0
      INITI2=0
      
      I_JMAT=0
      I_JMAT_PROTON=0
      I_JMAT_NEUTRS=0
C
C     The meaning of the variables below:
C
C     CHISQU_PROTON and CHISQU_NEUTRS =>> total chi^2 summed on all
C                                         the components and nuclei
C
      CHISQU_PROTON=0.0 ! As above, for protons
      CHISQU_NEUTRS=0.0 ! As above, for neutrosn
      CHISQU_TOTALS=0.0 ! As above, protons and neutrons togehter
C
C=======================================================================
C
C     First, calculating the total weights (protons and neutrons)
C     that will be used as the normalizing constants of the \chi^2
C
      WEIGHT_PROTON=0.0
      WEIGHT_NEUTRS=0.0
C
      DEGPRO_NOWEIG=0.0
      DEGNEU_NOWEIG=0.0
C
      SUMWEI_PROTON=0.0
      SUMWEI_NEUTRS=0.0
C
      DO INUCLI=1,LDNUCL
         IF (ITAKNU(INUCLI).EQ.1) THEN
C_______________________________________________________________________
C
             IF (IF_SPE.EQ.1) THEN
C
                 DEGSUM_PROTON=0.0
C
                 DO IEXPER=1,LEVEXP_PROTON(INUCLI)
                    DEGAUX=REAL(IDEGEX_PROTON(INUCLI,IEXPER))
                    DEGSUM_PROTON=DEGSUM_PROTON+DEGAUX
                 END DO
C
                 DEGPRO_NOWEIG=DEGPRO_NOWEIG+DEGSUM_PROTON
C
                 WEIGHT_PROTON=WEIGHT_PROTON+DEGSUM_PROTON
     *                                      *WEINUC_PROTON(INUCLI)
C
                 DEGSUM_NEUTRS=0.0
C
                 DO IEXPER=1,LEVEXP_NEUTRS(INUCLI)
                    DEGAUX=REAL(IDEGEX_NEUTRS(INUCLI,IEXPER))
                    DEGSUM_NEUTRS=DEGSUM_NEUTRS+DEGAUX
                 END DO
C
                 DEGNEU_NOWEIG=DEGNEU_NOWEIG+DEGSUM_NEUTRS
C
                 WEIGHT_NEUTRS=WEIGHT_NEUTRS+DEGSUM_NEUTRS
     *                                      *WEINUC_NEUTRS(INUCLI)
C
             END IF
C_______________________________________________________________________
C
             IF (IF_RAD.EQ.1) THEN
C
                 WEIGHT_PROTON=WEIGHT_PROTON+WEIRAD_PROTON(INUCLI)
                 WEIGHT_NEUTRS=WEIGHT_NEUTRS+WEIRAD_NEUTRS(INUCLI)
C
             END IF
C_______________________________________________________________________
C
             SUMWEI_PROTON=SUMWEI_PROTON+WEINUC_PROTON(INUCLI)
             SUMWEI_NEUTRS=SUMWEI_NEUTRS+WEINUC_NEUTRS(INUCLI)
C
         END IF
      END DO
C
C=======================================================================
C    
      DO INUCLI=1,LDNUCL
         
         IF (ITAKNU(INUCLI).EQ.1) THEN
             
             IZ_FIX=NUMB_Z(INUCLI)
             IN_FIX=NUMB_N(INUCLI)
C
C=======================================================================
C
             IF (LOGWRI.GT.4) THEN
                 WRITE(LOGFIL,'(12X,''Entering WS_RUN from JMATRX '',
     *                              ''[1] IZ_FIX='',I3,
     *                                 '' IN_FIX='',I3)')
     *                                    IZ_FIX,IN_FIX
             END IF
                                                     I_MODE=2
                                                            I_FLAG=1
             CALL WS_RUN(INUCLI,CHISQU_CHOICE,IABORT,I_MODE,I_FLAG,
     *                                 CHISQU_AUXIL1,CHISQU_AUXIL2)      
C
             IF (IDEFCN.GE.ITECHI) THEN
                 WRITE(LSCREN,'(/,''Omitting JMATRX, no-convergence '',
     *                       ''case -> RETURN from JMATRX'',/)')
                 GO TO 100
CID                 RETURN
             END IF
C
C=======================================================================                 
C
C            The following quantities contain ALL \chi^2 contributions 
C            They already contain their corresponding weight factors
C
C            These quantites are the ones printed in the LOGAUX files
C        
             CHISQU_PROTON=CHISQU_PROTON+(CHISQU_AUXIL1/SUMWEI_PROTON)
C
             CHISQU_NEUTRS=CHISQU_NEUTRS+(CHISQU_AUXIL2/SUMWEI_NEUTRS)
C
             CHISQU_TOTALS=CHISQU_TOTALS+(CHISQU_AUXIL1/SUMWEI_PROTON)
     *                                  +(CHISQU_AUXIL2/SUMWEI_NEUTRS)
C_______________________________________________________________________
C
C            These quantites are the ones printed in the LOGAUX files
C
             CHIENE_WEIPRO=CHIENE_WEIPRO+CHIDEG_PROTON
             CHIENE_WEINEU=CHIENE_WEINEU+CHIDEG_NEUTRS
C
             CHIRAD_WEIPRO=CHIRAD_WEIPRO+RADDIF_PROTON
     *                                  *WEIRAD_PROTON(INUCLI)
             CHIRAD_WEINEU=CHIRAD_WEINEU+RADDIF_NEUTRS
     *                                  *WEIRAD_NEUTRS(INUCLI)
C
             CHIRHO_WEIPRO=CHIRHO_WEIPRO+CHIRHO_PROTON*WEIGHT_RHODEN
C
C=======================================================================                 
C            Beginning the IF sequence to choose the selected
C            components of the total chi^2 = chi_1^2 + ch_2^2 + ...  
C=======================================================================                 
C         
             IF (IF_SPE.EQ.1) THEN
C            
                 DO IEXPER=1,LEVEXP_PROTON(INUCLI)
                    DO ITHEOR=1,LEVTHE_PROTON
C                   
                       IF (LABEXP_PROTON(INUCLI,IEXPER).EQ.
     *                     LABTHE_PROTON(ITHEOR)) THEN
C                      
                           ILEVEL=ILEVEL+1
                           LEVPRO=ILEVEL
C                      
                           ILEVEL_PROTON=ILEVEL_PROTON+1
                           IDPROT(ILEVEL_PROTON)=ILEVEL
C                      
                           IF (IEXPER.EQ.1) INITI1=ILEVEL
C                      
                           LEVELS(ILEVEL)=ITHEOR
                           LEXPLS(ILEVEL)=IEXPER
C                      
                           DEGENE(ILEVEL)
     *                    =
     *                     REAL(IDEGEX_PROTON(INUCLI,IEXPER))
C                      
                           FUNCT0(ILEVEL)=ENETHE_PROTON(ITHEOR)
                           FUNEXP(ILEVEL)=EXPEXP_PROTON(INUCLI,IEXPER)
C                      
                           WEISQR(ILEVEL)
     *                    =
     *                     SQRT(DEGENE(ILEVEL)*WEINUC_PROTON(INUCLI)
     *                                        /SUMWEI_PROTON        ) 
                       END IF
C                   
                    END DO
                 END DO
C
                 DO IEXPER=1,LEVEXP_NEUTRS(INUCLI)
                    DO ITHEOR=1,LEVTHE_NEUTRS
C                   
                       IF (LABEXP_NEUTRS(INUCLI,IEXPER).EQ.
     *                     LABTHE_NEUTRS(ITHEOR)) THEN
C                      
                           ILEVEL=ILEVEL+1
                           LEVNEU=ILEVEL
C                      
                           ILEVEL_NEUTRS=ILEVEL_NEUTRS+1
                           IDNEUT(ILEVEL_NEUTRS)=ILEVEL
C                      
                           IF (IEXPER.EQ.1) INITI2=ILEVEL
C                      
                           LEVELS(ILEVEL)=ITHEOR
                           LEXPLS(ILEVEL)=IEXPER
C                      
                           DEGENE(ILEVEL)
     *                    =
     *                     REAL(IDEGEX_NEUTRS(INUCLI,IEXPER))
C                      
                           FUNCT0(ILEVEL)=ENETHE_NEUTRS(ITHEOR)
                           FUNEXP(ILEVEL)=EXPEXP_NEUTRS(INUCLI,IEXPER)
C                      
                           WEISQR(ILEVEL)
     *                    =
     *                     SQRT(DEGENE(ILEVEL)*WEINUC_NEUTRS(INUCLI)
     *                                        /SUMWEI_NEUTRS        )
C                      
                       END IF
C                   
                    END DO
                 END DO
C            
             END IF ! IF_SPE=1
C         
C=======================================================================
C 
             IF (IF_RAD.EQ.1) THEN
C
                 ILEVEL=ILEVEL+1
C
                 ILEVEL_PROTON=ILEVEL_PROTON+1
                 IDPROT(ILEVEL_PROTON)=ILEVEL
C
                 FUNCT0(ILEVEL)=RMSTHE_PROTON
                 FUNEXP(ILEVEL)=RMSEXP_PROTON(INUCLI)
                 WEISQR(ILEVEL)=SQRT(WEIRAD_PROTON(INUCLI)
     *                         *     WEINUC_PROTON(INUCLI)
     *                         /     SUMWEI_PROTON) 
C
                 ILEVEL=ILEVEL+1
C
                 ILEVEL_NEUTRS=ILEVEL_NEUTRS+1
                 IDNEUT(ILEVEL_NEUTRS)=ILEVEL
C
                 FUNCT0(ILEVEL)=RMSTHE_NEUTRS
                 FUNEXP(ILEVEL)=RMSEXP_NEUTRS(INUCLI)
                 WEISQR(ILEVEL)=SQRT(WEIRAD_NEUTRS(INUCLI)
     *                         *     WEINUC_NEUTRS(INUCLI)
     *                         /     SUMWEI_NEUTRS) 
C
             END IF ! IF_RAD=1
C         
C=======================================================================
C
             IF (IF_GAP.EQ.1) THEN
C
                 ILEVEL=ILEVEL+1
C
                 ILEVEL_PROTON=ILEVEL_PROTON+1
                 IDPROT(ILEVEL_PROTON)=ILEVEL
C
                 FUNCT0(ILEVEL)=GAPTHE_PROTON
                 FUNEXP(ILEVEL)=GAPEXP_PROTON(INUCLI)
                 WEISQR(ILEVEL)=SQRT(WEIGHT_ENEGAP/WEIGHT_PROTON) 
C
                 ILEVEL=ILEVEL+1
C
                 ILEVEL_NEUTRS=ILEVEL_NEUTRS+1
                 IDNEUT(ILEVEL_NEUTRS)=ILEVEL
C
                 FUNCT0(ILEVEL)=GAPTHE_NEUTRS
                 FUNEXP(ILEVEL)=GAPEXP_NEUTRS(INUCLI)
                 WEISQR(ILEVEL)=SQRT(WEIGHT_ENEGAP/WEIGHT_NEUTRS)
C
             END IF ! IF_GAP=1
C         
C=======================================================================
C
             IF (IF_FER.EQ.1) THEN
C
                 ILEVEL=ILEVEL+1
C
                 ILEVEL_PROTON=ILEVEL_PROTON+1
                 IDPROT(ILEVEL_PROTON)=ILEVEL
C
                 FUNCT0(ILEVEL)=FERTHE_PROTON
                 FUNEXP(ILEVEL)=FERMEX_PROTON(INUCLI)
                 WEISQR(ILEVEL)=SQRT(WEIGHT_EFERMI/WEIGHT_PROTON)
C
                 ILEVEL=ILEVEL+1
C
                 ILEVEL_NEUTRS=ILEVEL_NEUTRS+1
                 IDNEUT(ILEVEL_NEUTRS)=ILEVEL
C
                 FUNCT0(ILEVEL)=FERTHE_NEUTRS
                 FUNEXP(ILEVEL)=FERMEX_NEUTRS(INUCLI)
                 WEISQR(ILEVEL)=SQRT(WEIGHT_EFERMI/WEIGHT_NEUTRS)
C
             END IF ! IF_FER=1
C         
C=======================================================================
C
             IF (IF_DEN.EQ.1) THEN
C
                 ILEVEL=ILEVEL+1
C
                 ILEVEL_PROTON=ILEVEL_PROTON+1
                 IDPROT(ILEVEL_PROTON)=ILEVEL
C
                 FUNCT0(ILEVEL)=DENLOW_PROTON
                 FUNEXP(ILEVEL)=DENSDW_PROTON(INUCLI)
                 WEISQR(ILEVEL)=SQRT(WEIGHT_DENSDW/WEIGHT_PROTON)
C
                 ILEVEL=ILEVEL+1
C
                 ILEVEL_PROTON=ILEVEL_PROTON+1
                 IDPROT(ILEVEL_PROTON)=ILEVEL
C
                 FUNCT0(ILEVEL)=DENUPP_PROTON
                 FUNEXP(ILEVEL)=DENSUP_PROTON(INUCLI)
                 WEISQR(ILEVEL)=SQRT(WEIGHT_DENSUP/WEIGHT_PROTON)
C_______________________________________________________________________
C            
                 ILEVEL=ILEVEL+1
C
                 ILEVEL_NEUTRS=ILEVEL_NEUTRS+1
                 IDNEUT(ILEVEL_NEUTRS)=ILEVEL
C
                 FUNCT0(ILEVEL)=DENLOW_NEUTRS
                 FUNEXP(ILEVEL)=DENSDW_NEUTRS(INUCLI)
                 WEISQR(ILEVEL)=SQRT(WEIGHT_DENSDW/WEIGHT_NEUTRS)
C
                 ILEVEL=ILEVEL+1
C
                 ILEVEL_NEUTRS=ILEVEL_NEUTRS+1
                 IDNEUT(ILEVEL_NEUTRS)=ILEVEL
C
                 FUNCT0(ILEVEL)=DENUPP_NEUTRS
                 FUNEXP(ILEVEL)=DENSUP_NEUTRS(INUCLI)
                 WEISQR(ILEVEL)=SQRT(WEIGHT_DENSUP/WEIGHT_NEUTRS)
C
             END IF ! IF_DEN=1  
C         
C=======================================================================
C  @@@ irene - what is this?
             IF (IF_RHO.EQ.1) THEN
*C            
*                CALL RHOTHE_DENSIT(INUCLI,IZ_FIX,IN_FIX)
*C             
*                DO K_NOYX=1,N_NUCL
*C                 
*                   I_NOYX=I_NUCL(K_NOYX)
*C                 
*                   DO ID_RHO=1,ND_RHO
*C                   
*                      ILEVEL=ILEVEL+1
*                      ILEVEL_PROTON=ILEVEL_PROTON+1
*                      IDPROT(ILEVEL_PROTON)=ILEVEL
*C                   
*                      FUNCT0(ILEVEL)=RHOTHE(ID_RHO,I_NOYX)
*                      FUNEXP(ILEVEL)=RHOEXP(ID_RHO,I_NOYX)
*                      WEISQR(ILEVEL)=SQRT(WEIGHT_RHODEN)
*C                   
*                   END DO
*C            
*                END DO 
C             
             END IF ! IF_RHO=1
C         
C=======================================================================
C
             IPARTL=ILEVEL ! partial dimension of FUNCT0 AND FUNEXP
C         
C=======================================================================
C
C            Calculating partial derivatives of the energies and other
C            observables taken into account in the minimisation
C
C=======================================================================
C
             IACTIV=0
             IACPRO=0
             IACNEU=0
C
             DO IPARAM=1,NDPARS
                DERIVF(IPARAM)=0.0
             END DO
C
             IAUXI1=I_JMAT_PROTON
             IAUXI2=IAUXI1
C
             IAUXI3=I_JMAT_NEUTRS
             IAUXI4=IAUXI3
C         
C=======================================================================
C=======================================================================
C            Main DO LOOP over the parameters of the fit    
C=======================================================================
C=======================================================================
C
             DO IPARAM=1,NDPARS
C            
                IF (IFTAKE(IPARAM).EQ.1) THEN
C                
                    IACTIV=IACTIV+1
C                
                    IF (IFPROT.EQ.1) THEN
                        IACPRO=IACPRO+1
                    END IF
C                
                    IF (IFNEUT.EQ.1) THEN
                        IACNEU=IACNEU+1
                    END IF
C                
                    PARPOT(IPARAM)=PARPOT(IPARAM)+DPARAM(IPARAM)
                    ACTION(IPARAM)='*'
C
                    IF (LOGWRI.GT.4) THEN
                        WRITE(LOGFIL,'(12X,''Entering WS_RUN from '',
     *                                     ''JMATRX [2] IPARAM='',I2,
     *                                  1X,''IZ_FIX='',I3,1X,
     *                                     ''IN_FIX='',I3)') IPARAM,
     *                                       IZ_FIX,IN_FIX
                    END IF
                                                            I_MODE=1
                                                            I_FLAG=5
                    CALL WS_RUN(INUCLI,CHISQU_CHOICE,IABORT,I_MODE,
     *                                        I_FLAG,AUXIL1,AUXIL2)
C         
C=======================================================================
C                
                    IF (IF_SPE.EQ.1) THEN
C      
                        DO INDEXX=INITI1,LEVPRO
                        
                           IAUXI1=IAUXI1+1
C                       
                           IWHICH=LEVELS(INDEXX)
C                       
                           FUNCT1(INDEXX)=ENETHE_PROTON(IWHICH)
C                       
                           DERIVF(IACTIV)
     *                    =
     *                     (FUNCT1(INDEXX)-FUNCT0(INDEXX))
     *                                    /DPARAM(IPARAM)
C     
                           IF (IFDENS.EQ.0 .AND. IFPROT.EQ.1) THEN
C
                               XJMATR_PROTON(IAUXI1,IACPRO)
     *                        =
     *                         DERIVF(IACTIV)*WEISQR(INDEXX)
                               LABPRO_PRINTG(IAUXI1)
     *                        =
     *                         'p '//LABTHE_NEUTRS(IWHICH)
                           END IF
C    
                           IF (IFDENS.EQ.1 .OR. 
     *                        (IFDENS.EQ.0 .AND. IFBOTH.EQ.1)) THEN
C
                               XJMATR(INDEXX,IACTIV)=DERIVF(IACTIV)
     *                                              *WEISQR(INDEXX)
                               LABTOT_PRINTG(INDEXX)
     *                        =
     *                        'p '//LABTHE_NEUTRS(IWHICH)
                           END IF
C                       
                        END DO ! of INDEXX
C
                        I_JMAT_PROTON=IAUXI1
                        IAUXI1=IAUXI2
C                    
                        DO INDEXX=INITI2,LEVNEU
C
                           IAUXI3=IAUXI3+1
C                       
                           IWHICH=LEVELS(INDEXX)
C                       
                           FUNCT1(INDEXX)=ENETHE_NEUTRS(IWHICH)
C                       
                           DERIVF(IACTIV)
     *                    =
     *                     (FUNCT1(INDEXX)-FUNCT0(INDEXX))
     *                                    /DPARAM(IPARAM)
C     
                           IF (IFDENS.EQ.0 .AND. IFNEUT.EQ.1) THEN
C
                               XJMATR_NEUTRS(IAUXI3,IACNEU)
     *                        =
     *                         DERIVF(IACTIV)*WEISQR(INDEXX)
                               LABNEU_PRINTG(IAUXI3)
     *                        =
     *                         'n '//LABTHE_NEUTRS(IWHICH)
                           END IF
C                       
                           IF (IFDENS.EQ.1 .OR. 
     *                        (IFDENS.EQ.0 .AND. IFBOTH.EQ.1)) THEN
C
                               XJMATR(INDEXX,IACTIV)=DERIVF(IACTIV)
     *                                             *WEISQR(INDEXX)
                               LABTOT_PRINTG(INDEXX)
     *                        =
     *                         'n '//LABTHE_NEUTRS(IWHICH)
                           END IF
C                       
                        END DO
C
                        I_JMAT_NEUTRS=IAUXI3
                        IAUXI3=IAUXI4
C
                        I_JMAT=LEVNEU
C                    
                    END IF ! IF_SPE=1
C         
C=======================================================================
C
                    IF (IF_RAD.EQ.1) THEN
C                    
                        I_JMAT=I_JMAT+1
                        I_JMAT_PROTON=I_JMAT_PROTON+1
C                    
                        FUNCT1(I_JMAT)=RMSTHE_PROTON
C                    
                        IF (IFDENS.EQ.0 .AND. IFPROT.EQ.1) THEN
C
                            XJMATR_PROTON(I_JMAT_PROTON,IACPRO)
     *                                         =(FUNCT1(I_JMAT)
     *                                         - FUNCT0(I_JMAT))
     *                                         / DPARAM(IPARAM)
     *                                         * WEISQR(I_JMAT)
                            LABPRO_PRINTG(I_JMAT_PROTON)='p radius'
                        END IF
C                       
                        IF (IFDENS.EQ.1 .OR. 
     *                     (IFDENS.EQ.0 .AND. IFBOTH.EQ.1)) THEN
C
                            XJMATR(I_JMAT,IACTIV)=(FUNCT1(I_JMAT)
     *                                           - FUNCT0(I_JMAT))
     *                                           / DPARAM(IPARAM)
     *                                           * WEISQR(I_JMAT)
                            LABTOT_PRINTG(I_JMAT)='p radius'
                        END IF
C_______________________________________________________________________
C                    
                        I_JMAT=I_JMAT+1
                        I_JMAT_NEUTRS=I_JMAT_NEUTRS+1
C
                        FUNCT1(I_JMAT)=RMSTHE_NEUTRS
C
                        IF (IFDENS.EQ.0 .AND. IFNEUT.EQ.1) THEN
C
                            XJMATR_NEUTRS(I_JMAT_NEUTRS,IACNEU)
     *                                         =(FUNCT1(I_JMAT)
     *                                         - FUNCT0(I_JMAT))
     *                                         / DPARAM(IPARAM)
     *                                         * WEISQR(I_JMAT)
                            LABNEU_PRINTG(I_JMAT_NEUTRS)='n radius'
                        END IF
C                       
                        IF (IFDENS.EQ.1 .OR. 
     *                     (IFDENS.EQ.0 .AND. IFBOTH.EQ.1)) THEN
C
                            XJMATR(I_JMAT,IACTIV)=(FUNCT1(I_JMAT)
     *                                           - FUNCT0(I_JMAT))
     *                                           / DPARAM(IPARAM)
     *                                           * WEISQR(I_JMAT)
                            LABTOT_PRINTG(I_JMAT)='n radius'
                        END IF
C
                    END IF ! IF_RAD=1
C         
C=======================================================================
C
                    IF (IF_GAP.EQ.1) THEN
C
                        I_JMAT=I_JMAT+1
                        I_JMAT_PROTON=I_JMAT_PROTON+1
C
                        FUNCT1(I_JMAT)=GAPTHE_PROTON
C
                        IF (IFDENS.EQ.0 .AND. IFPROT.EQ.1) THEN
C
                            XJMATR_PROTON(I_JMAT_PROTON,IACPRO)
     *                                         =(FUNCT1(I_JMAT)
     *                                         - FUNCT0(I_JMAT))
     *                                         / DPARAM(IPARAM)
     *                                         * WEISQR(I_JMAT)
                        END IF
C                       
                        IF (IFDENS.EQ.1 .OR.
     *                     (IFDENS.EQ.0 .AND. IFBOTH.EQ.1)) THEN
C
                            XJMATR(I_JMAT,IACTIV)=(FUNCT1(I_JMAT)
     *                                           - FUNCT0(I_JMAT))
     *                                           / DPARAM(IPARAM)
     *                                           * WEISQR(I_JMAT)
                        END IF
C_______________________________________________________________________
C                    
                        I_JMAT=I_JMAT+1
                        I_JMAT_NEUTRS=I_JMAT_NEUTRS+1
C
                        FUNCT1(I_JMAT)=GAPTHE_NEUTRS
C
                        IF (IFDENS.EQ.0 .AND. IFNEUT.EQ.1) THEN
C
                            XJMATR_NEUTRS(I_JMAT_NEUTRS,IACNEU)
     *                                         =(FUNCT1(I_JMAT)
     *                                         - FUNCT0(I_JMAT))
     *                                         / DPARAM(IPARAM)
     *                                         * WEISQR(I_JMAT)
                        END IF
C                       
                        IF (IFDENS.EQ.1 .OR.
     *                     (IFDENS.EQ.0 .AND. IFBOTH.EQ.1)) THEN
C
                            XJMATR(I_JMAT,IACTIV)=(FUNCT1(I_JMAT)
     *                                           - FUNCT0(I_JMAT))
     *                                           / DPARAM(IPARAM)
     *                                           * WEISQR(I_JMAT)
                        END IF
C
                    END IF ! IF_GAP=1
C         
C=======================================================================
C
                    IF (IF_FER.EQ.1) THEN
C
                        I_JMAT=I_JMAT+1
                        I_JMAT_PROTON=I_JMAT_PROTON
C
                        FUNCT1(I_JMAT)=FERTHE_PROTON
C
                        IF (IFDENS.EQ.0 .AND. IFPROT.EQ.1) THEN
C
                            XJMATR_PROTON(I_JMAT_PROTON,IACPRO)
     *                                         =(FUNCT1(I_JMAT)
     *                                         - FUNCT0(I_JMAT))
     *                                         / DPARAM(IPARAM)
     *                                         * WEISQR(I_JMAT)
                        END IF
C                       
                        IF (IFDENS.EQ.1 .OR.
     *                     (IFDENS.EQ.0 .AND. IFBOTH.EQ.1)) THEN
C
                            XJMATR(I_JMAT,IACTIV)=(FUNCT1(I_JMAT)
     *                                           - FUNCT0(I_JMAT))
     *                                           / DPARAM(IPARAM)
     *                                           * WEISQR(I_JMAT)
                        END IF
C_______________________________________________________________________
C                    
                        I_JMAT=I_JMAT+1
                        I_JMAT_NEUTRS=I_JMAT_NEUTRS+1
C
                        FUNCT1(I_JMAT)=FERTHE_NEUTRS
C
                        IF (IFDENS.EQ.0 .AND. IFNEUT.EQ.1) THEN
C
                            XJMATR_NEUTRS(I_JMAT_NEUTRS,IACNEU)
     *                                         =(FUNCT1(I_JMAT)
     *                                         - FUNCT0(I_JMAT))
     *                                         / DPARAM(IPARAM)
     *                                         * WEISQR(I_JMAT)
                        END IF
C                       
                        IF (IFDENS.EQ.1 .OR.
     *                     (IFDENS.EQ.0 .AND. IFBOTH.EQ.1)) THEN
C
                            XJMATR(I_JMAT,IACTIV)=(FUNCT1(I_JMAT)
     *                                           - FUNCT0(I_JMAT))
     *                                           / DPARAM(IPARAM)
     *                                           * WEISQR(I_JMAT)
                        END IF
C
                    END IF ! IF_FER=1
C
C=======================================================================
C
                    IF (IF_DEN.EQ.1) THEN
C
                        I_JMAT=I_JMAT+1
                        I_JMAT_PROTON=I_JMAT_PROTON
C
                        FUNCT1(I_JMAT)=DENLOW_PROTON
C
                        IF (IFDENS.EQ.0 .AND. IFPROT.EQ.1) THEN
C
                            XJMATR_PROTON(I_JMAT_PROTON,IACPRO)
     *                                         =(FUNCT1(I_JMAT)
     *                                         - FUNCT0(I_JMAT))
     *                                         / DPARAM(IPARAM)
     *                                         * WEISQR(I_JMAT)
                        END IF
C                       
                        IF (IFDENS.EQ.1 .OR.
     *                     (IFDENS.EQ.0 .AND. IFBOTH.EQ.1)) THEN
C
                            XJMATR(I_JMAT,IACTIV)=(FUNCT1(I_JMAT)
     *                                           - FUNCT0(I_JMAT))
     *                                           / DPARAM(IPARAM)
     *                                           * WEISQR(I_JMAT)
                        END IF
C_______________________________________________________________________
C                    
                        I_JMAT=I_JMAT+1
                        I_JMAT_PROTON=I_JMAT_PROTON
C
                        FUNCT1(I_JMAT)=DENUPP_PROTON
C
                        IF (IFDENS.EQ.0 .AND. IFPROT.EQ.1) THEN
C
                            XJMATR_PROTON(I_JMAT_PROTON,IACPRO)
     *                                         =(FUNCT1(I_JMAT)
     *                                         - FUNCT0(I_JMAT))
     *                                         / DPARAM(IPARAM)
     *                                         * WEISQR(I_JMAT)
                        END IF
C                       
                        IF (IFDENS.EQ.1 .OR.
     *                     (IFDENS.EQ.0 .AND. IFBOTH.EQ.1)) THEN
C
                            XJMATR(I_JMAT,IACTIV)=(FUNCT1(I_JMAT)
     *                                           - FUNCT0(I_JMAT))
     *                                           / DPARAM(IPARAM)
     *                                           * WEISQR(I_JMAT)
                        END IF
C_______________________________________________________________________
C
                        I_JMAT=I_JMAT+1
                        I_JMAT_NEUTRS=I_JMAT_NEUTRS+1
C  
                        FUNCT1(I_JMAT)=DENLOW_NEUTRS
C
                        IF (IFDENS.EQ.0 .AND. IFNEUT.EQ.1) THEN
C
                            XJMATR_NEUTRS(I_JMAT_NEUTRS,IACNEU)
     *                                         =(FUNCT1(I_JMAT)
     *                                         - FUNCT0(I_JMAT))
     *                                         / DPARAM(IPARAM)
     *                                         * WEISQR(I_JMAT)
                        END IF
C                       
                        IF (IFDENS.EQ.1 .OR.
     *                     (IFDENS.EQ.0 .AND. IFBOTH.EQ.1)) THEN
C
                            XJMATR(I_JMAT,IACTIV)=(FUNCT1(I_JMAT)
     *                                           - FUNCT0(I_JMAT))
     *                                           / DPARAM(IPARAM)
     *                                           * WEISQR(I_JMAT)
                        END IF
C_______________________________________________________________________
C
                        I_JMAT=I_JMAT+1
                        I_JMAT_NEUTRS=I_JMAT_NEUTRS+1
C
                        FUNCT1(I_JMAT)=DENUPP_NEUTRS
C
                        IF (IFDENS.EQ.0 .AND. IFNEUT.EQ.1) THEN
C
                            XJMATR_NEUTRS(I_JMAT_NEUTRS,IACNEU)
     *                                         =(FUNCT1(I_JMAT)
     *                                         - FUNCT0(I_JMAT))
     *                                         / DPARAM(IPARAM)
     *                                         * WEISQR(I_JMAT)
                        END IF
C                       
                        IF (IFDENS.EQ.1 .OR.
     *                     (IFDENS.EQ.0 .AND. IFBOTH.EQ.1)) THEN
C
                            XJMATR(I_JMAT,IACTIV)=(FUNCT1(I_JMAT)
     *                                           - FUNCT0(I_JMAT))
     *                                           / DPARAM(IPARAM)
     *                                           * WEISQR(I_JMAT)
                        END IF
C                
                    END IF ! IF_DEN=1
C         
C=======================================================================
C
C  @@@ IRENE: WHAT DOES THIS MEAN?
                    IF (IF_RHO.EQ.1) THEN
*C            
*                       CALL RHOTHE_DENSIT(INUCLI,IZ_FIX,IN_FIX)
*C             
*                       DO K_NOYX=1,N_NUCL
*C                 
*                          I_NOYX=I_NUCL(K_NOYX)
*C                 
*                          DO ID_RHO=1,ND_RHO
*C                   
*                             I_JMAT=I_JMAT+1
*                             I_JMAT_PROTON=I_JMAT_PROTON+1
                    
*                             FUNCT1(I_JMAT)=RHOTHE(ID_RHO,I_NOYX)
                          
*                             IF (IFDENS.EQ.0 .AND. IPARAM.LE.20) THEN
C
*                                 XJMATR_PROTON(I_JMAT_PROTON,IACPRO)
*     *                                              =(FUNCT1(I_JMAT)
*     *                                              - FUNCT0(I_JMAT))
*     *                                              / DPARAM(IPARAM)
*     *                                              * WEISQR(I_JMAT)
*                             END IF
*C                       
*                             IF (IFDENS.EQ.1) THEN
C
*                                 XJMATR(I_JMAT,IACTIV)=(FUNCT1(I_JMAT)
*     *                                                - FUNCT0(I_JMAT))
*     *                                                / DPARAM(IPARAM)
*     *                                                * WEISQR(I_JMAT)
*                             END IF
*C                   
*                          END DO
*                       END DO
*C                      
                    END IF ! of IF_RHO
C                
C                   IF (IF_INV.EQ.1) THEN
CID                 END IF
C         
C=======================================================================
C                
                    ACTION(IPARAM)=' ' 
                    PARPOT(IPARAM)=PARPOT(IPARAM)-DPARAM(IPARAM)
            
                END IF ! (IFTAKE(IPARAM).EQ.1)
C            
             END DO ! IPARAM
C         
C=======================================================================
C         
             IF (IPARTL.NE.I_JMAT) THEN
                 WRITE(0,'(/,''Alarm in JMATRX: IPARTL= '',I3,
     *                 '' is different from I_JMAT= '',I3,
     *                 '' -> They should be equal!'',/)') IPARTL,I_JMAT
                 STOP 'STOP in JMATRX: IPARTL.NE.I_JMAT'
             END IF
C         
C=======================================================================
C     
         END IF ! Over (ITAKNU(INUCLI).EQ.1)
      
      END DO ! Over INUCLI
C
C=======================================================================
C=======================================================================
C=======================================================================
C      
      LEVNUM=ILEVEL
C      
      ITOTAL=ILEVEL_PROTON+ILEVEL_NEUTRS
C
      VALCHI=CHISQU_TOTALS
C
C=======================================================================
C                       Some 'alarm' prints 
C=======================================================================
C     
      IF (IFDENS.EQ.1 .OR. (IFDENS.EQ.0 .AND. IFBOTH.EQ.1)) THEN
          
          IF (LEVNUM.NE.LDFUNC) THEN
              WRITE(LSCREN,'(/,''Alarm in JMATRX: LEVNUM= '',I3,
     *                 '' is different from LDFUNC= '',I3,
     *                 '' -> They should be equal!'',/)') LEVNUM,LDFUNC
              STOP 'STOP in JMATRX: LEVNUM.NE.LDFUNC'
          END IF
C_______________________________________________________________________
C      
          IF (I_JMAT.NE.LDFUNC) THEN
              WRITE(LSCREN,'(/,''Alarm in JMATRX: I_JMAT= '',I3,
     *                '' is different from LDFUNC= '',I3,
     *                '' -> They should be equal!'',/)') I_JMAT,LDFUNC
              STOP 'STOP in JMATRX: I_JMAT.NE.LDFUNC'
          END IF
C_______________________________________________________________________
C      
          IF (ITOTAL.NE.LDFUNC) THEN
              WRITE(LSCREN,'(/,''Alarm in JMATRX: '',
     *                ''ITOTAL=ILEVEL_PROTON+ILEVEL_NEUTRS= '',
     *                  I3,'' = '',I3,'' + '',I3,
     *                '' is different from LDFUNC= '',I3,
     *                '' -> They should be equal!'',/)') ITOTAL,
     *                                                   ILEVEL_PROTON,
     *                                                   ILEVEL_NEUTRS,
     *                                                          LDFUNC
             STOP 'STOP in JMATRX: ITOTAL.NE.LDFUNC'
          END IF
      
      END IF
C
C=======================================================================
C      
      IF (IFDENS.EQ.0 .AND. IFPROT.EQ.1) THEN
C          
          IF (ILEVEL_PROTON.NE.LDFUNC) THEN
              WRITE(LSCREN,'(/,''Alarm in JMATRX: ILEVEL_PROTON= '',I3,
     *                 '' is different from LDFUNC= '',I3,
     *                 '' -> They should be equal!'',/)')ILEVEL_PROTON,
     *                                                          LDFUNC
              STOP 'STOP in JMATRX: ILEVEL_PROTON.NE.LDFUNC'
          END IF
C_______________________________________________________________________
C      
          IF (I_JMAT_PROTON.NE.LDFUNC) THEN
              WRITE(LSCREN,'(/,''Alarm in JMATRX: I_JMAT_PROTON= '',I3,
     *                '' is different from LDFUNC= '',I3,
     *                '' -> They should be equal!'',/)') I_JMAT_PROTON,
     *                                                          LDFUNC
              STOP 'STOP in JMATRX: I_JMAT.NE.LDFUNC'
          END IF
C      
      END IF      
C
C=======================================================================
C      
      IF (IFDENS.EQ.0 .AND. IFNEUT.EQ.1) THEN
C
          IF (ILEVEL_NEUTRS.NE.LDFUNC) THEN
              WRITE(LSCREN,'(/,''Alarm in JMATRX: ILEVEL_NEUTRS= '',I3,
     *                 '' is different from LDFUNC= '',I3,
     *                 '' -> They should be equal!'',/)')ILEVEL_NEUTRS,
     *                                                          LDFUNC
              STOP 'STOP in JMATRX: ILEVEL_NEUTRS.NE.LDFUNC'
          END IF
C_______________________________________________________________________
C      
          IF (I_JMAT_NEUTRS.NE.LDFUNC) THEN
              WRITE(LSCREN,'(/,''Alarm in JMATRX: I_JMAT_NEUTRS= '',I3,
     *                '' is different from LDFUNC= '',I3,
     *                '' -> They should be equal!'',/)') I_JMAT_NEUTRS,
     *                                                          LDFUNC
              STOP 'STOP in JMATRX: I_JMAT.NE.LDFUNC'
          END IF
C      
      END IF
C
C=======================================================================
C      
                                I_MODE=2
      CALL INPRIN(IDEFCN,IACTIV,I_MODE,CHISQU_PROTON,CHISQU_NEUTRS,
     *                                 CHIENE_WEIPRO,CHIENE_WEINEU,
     *                                 CHIRAD_WEIPRO,CHIRAD_WEINEU,
     *                                               CHIRHO_WEIPRO)  
C
C=======================================================================
C          Calculating the chi^2 gradient
C=======================================================================
C     
      IF (IFDENS.EQ.0 .AND. IFPROT.EQ.1) THEN
C 
           CHISQU_GRDNRM=0.0 ! \chi^2 gradient norm
C
           DO INDEXP=1,IACPRO
               CHISQU_GRADNT(INDEXP)=0.0
           END DO
C_______________________________________________________________________
C        
           DO INDEXP=1,IACPRO ! parameters
C        
              DO INDEXF=1,ILEVEL_PROTON ! functions
C
                 ILEVEL=IDPROT(INDEXF)
C            
                 FUNAUX(INDEXF)=(FUNCT0(ILEVEL)-FUNEXP(ILEVEL))
     *                         * WEISQR(ILEVEL)
C
                 CHISQU_GRADNT(INDEXP)=CHISQU_GRADNT(INDEXP)
     *                                +FUNAUX(INDEXF)
     *                                *XJMATR_PROTON(INDEXF,INDEXP)
C
             END DO
C         
             CHISQU_GRADNT(INDEXP)=2*CHISQU_GRADNT(INDEXP)
             CHISQU_GRDNRM=CHISQU_GRDNRM+(CHISQU_GRADNT(INDEXP))**2
C         
          END DO 
C
          CHISQU_GRDNRM=SQRT(CHISQU_GRDNRM)
C
      END IF
C
C=======================================================================
C
      IF (IFDENS.EQ.0 .AND. IFNEUT.EQ.1) THEN
C 
           CHISQU_GRDNRM=0.0 ! \chi^2 gradient norm
C
           DO INDEXP=1,IACNEU
               CHISQU_GRADNT(INDEXP)=0.0
           END DO
C________________________________________________________________________
C        
           DO INDEXP=1,IACNEU ! parameters
C        
              DO INDEXF=1,ILEVEL_NEUTRS ! functions
C
                 ILEVEL=IDNEUT(INDEXF)
C            
                 FUNAUX(INDEXF)=(FUNCT0(ILEVEL)-FUNEXP(ILEVEL))
     *                         * WEISQR(ILEVEL)
C
                 CHISQU_GRADNT(INDEXP)=CHISQU_GRADNT(INDEXP)
     *                                +FUNAUX(INDEXF)
     *                                *XJMATR_NEUTRS(INDEXF,INDEXP)
C
             END DO
C         
             CHISQU_GRADNT(INDEXP)=2*CHISQU_GRADNT(INDEXP)
             CHISQU_GRDNRM=CHISQU_GRDNRM+(CHISQU_GRADNT(INDEXP))**2
C         
          END DO 
C
          CHISQU_GRDNRM=SQRT(CHISQU_GRDNRM)
C
      END IF
C
C=======================================================================
C
      IF (IFDENS.EQ.1 .OR. (IFDENS.EQ.0 .AND. IFBOTH.EQ.1)) THEN
C 
           CHISQU_GRDNRM=0.0 ! \chi^2 gradient norm
C
           DO INDEXP=1,IACTIV
               CHISQU_GRADNT(INDEXP)=0.0
           END DO
C________________________________________________________________________
C        
           DO INDEXP=1,IACTIV ! parameters
C        
              DO INDEXF=1,LDFUNC ! functions
C            
                 FUNAUX(INDEXF)=(FUNCT0(INDEXF)-FUNEXP(INDEXF))
     *                         * WEISQR(INDEXF)
C
                 CHISQU_GRADNT(INDEXP)=CHISQU_GRADNT(INDEXP)
     *                                +FUNAUX(INDEXF)
     *                                *XJMATR(INDEXF,INDEXP)
C
             END DO
C         
             CHISQU_GRADNT(INDEXP)=2*CHISQU_GRADNT(INDEXP)
             CHISQU_GRDNRM=CHISQU_GRDNRM+(CHISQU_GRADNT(INDEXP))**2
C         
          END DO 
C
          CHISQU_GRDNRM=SQRT(CHISQU_GRDNRM)
C
      END IF 
C
C=======================================================================
C     Checking whether we proceed with the estimation of the 
C     confidence intervals for parameters.
C     If number of parameters is greater than number of experimental
C     levels (plus the rms radii) further calculations are pointless
C=======================================================================
C
      IF (IFDENS.EQ.0) THEN
C
          IF (IACPRO.GE.ILEVEL_PROTON) THEN
C
              IF (LOGWRI.GT.0) 
     *           WRITE(IRESUL,'(''The confidence intervals CANNOT be '',
     *                       ''calculated'',/,''NUMBER OF PARAMETERS '',
     *                       ''> NUMBER OF LEVELS --> PROTONS'')')
              WRITE(LSCREN,'(/,''RETURN in JMATRX: '',
     *                  ''The confidence intervals CANNOT be '',
     *                  ''calculated'',/,''NUMBER OF PARAMETERS '',
     *                  ''> NUMBER OF LEVELS --> PROTONS'',/)')
C
               GO TO 100
C              CALL CPUTIM('JMATRX',0)
C              RETURN
          END IF
C
          IF (IACNEU.GE.ILEVEL_NEUTRS) THEN
C
              IF (LOGWRI.GT.0) 
     *           WRITE(IRESUL,'(''The confidence intervals CANNOT be '',
     *                       ''calculated'',/,''NUMBER OF PARAMETERS '',
     *                       ''> NUMBER OF LEVELS --> NEUTRONS'')')
              WRITE(LSCREN,'(/,''RETURN in JMATRX: '',
     *                  ''The confidence intervals CANNOT be '',
     *                  ''calculated'',/,''NUMBER OF PARAMETERS '',
     *                  ''> NUMBER OF LEVELS --> NEUTRS'',/)')
C
               GO TO 100
C              CALL CPUTIM('JMATRX',0)
C              RETURN
          END IF
      ELSE
          IF (IACTIV.GE.LEVNUM) THEN
C
              IF (LOGWRI.GT.0) 
     *           WRITE(IRESUL,'(''The confidence intervals CANNOT be '',
     *                       ''calculated'',/,''NUMBER OF PARAMETERS '',
     *                       ''> NUMBER OF LEVELS'')')
              WRITE(LSCREN,'(/,''RETURN in JMATRX: '',
     *                  ''The confidence intervals CANNOT be '',
     *                  ''calculated'',/,''NUMBER OF PARAMETERS '',
     *                  ''> NUMBER OF LEVELS'',/)')
C
               GO TO 100
C              CALL CPUTIM('JMATRX',0)
C              RETURN
          END IF
      END IF
C
C=======================================================================
C     Calculating the S^2 and the critical t...
C=======================================================================
C
      PROBAB=0.975
C
      IF (IFDENS.EQ.1) THEN
C
          FREDEG=FLOAT(LEVNUM-IACTIV+1) !!!check definition
          SUMSQA=0.
C
          DO INDEXX=1,LEVNUM
             AUXILL=WEISQR(INDEXX)*(FUNCT0(INDEXX)-FUNEXP(INDEXX))
             SUMSQA=SUMSQA+AUXILL*AUXILL
          END DO
C
          SUMSQA=SUMSQA/FREDEG
C
          CALL TCRITC(CRITIC)
C
          DELTAP=SQRT(SUMSQA)*CRITIC!/SQRT(FREDEG)
C
      END IF
C
C=======================================================================
C
      IF (IFDENS.EQ.0) THEN
C
          FREDEG=FLOAT(ILEVEL_PROTON-IACPRO+1)
          SUMSQA_PROTON=0.
C
          DO INDEXX=1,ILEVEL_PROTON
             INDEX2=IDPROT(INDEXX)
             AUXILL=WEISQR(INDEX2)
     *             *(FUNCT0(INDEX2)-FUNEXP(INDEX2))
             SUMSQA_PROTON=SUMSQA_PROTON+AUXILL*AUXILL
          END DO
C
          SUMSQA_PROTON=SUMSQA_PROTON/FREDEG
C
          CALL TCRITC(CRITIC)
C
          DELTAP_PROTON=SQRT(SUMSQA_PROTON)*CRITIC!/SQRT(FREDEG)
C_______________________________________________________________________
C
          FREDEG=FLOAT(ILEVEL_NEUTRS-IACNEU+1)
          SUMSQA_NEUTRS=0.
C
          DO INDEXX=1,ILEVEL_NEUTRS
             INDEX3=IDNEUT(INDEXX)
             AUXILL=WEISQR(INDEX3)
     *             *(FUNCT0(INDEX3)-FUNEXP(INDEX3))
             SUMSQA_NEUTRS=SUMSQA_NEUTRS+AUXILL*AUXILL
          END DO
C
          SUMSQA_NEUTRS=SUMSQA_NEUTRS/FREDEG
C
          CALL TCRITC(CRITIC)
C
          DELTAP_NEUTRS=SQRT(SUMSQA_NEUTRS)*CRITIC!/SQRT(FREDEG)
C
      END IF
C
C=======================================================================
C
C     Printing the J(k,l) matrix...
C
C=======================================================================
C
      IF (LOGWRI.GT.4) THEN
          WRITE(IRESUL,'(/,80(''*''),/,80(''*''),/,''*'',78X,''*'',/,
     *                 ''*   RESULTS OF THE ERROR ANALYSIS'',
     *                 T80,''*'',/,''*'',78X,''*'',/,80(''*''))')
C
          WRITE(IRESUL,'(80(''*''),/,''*'',78X,''*'',/,
     *                 ''*   Jacobian for '',I2,'' active parameters '',
     *                 ''and '',I2,'' proton and/or neutron levels:'',
     *                 T80,''*'',/,''*'',78X,''*'')')IACTIV,LEVNUM-2
      END IF
C
      IF (IFDENS.EQ.0) THEN
C
          IF (IFPROT.EQ.1) THEN
C          
              WRITE(LOGAUX,'(/,''<<JACOBMATRX>>    '',
     *                       <IACPRO>(''Par'',i2.2,4X))')(I,I=1,IACPRO)
C
              DO I_LEVL=1,ILEVEL_PROTON
C                WRITE(111,'(<IACPRO>(F15.8))')
C     *           (XJMATR_PROTON(I_LEVL,I_PARM),I_PARM=1,IACPRO)
                 IF (LOGWRI.GT.4)
     *               WRITE(IRESUL,'(''*'',<IACPRO>(F9.3),T80,''*'')')
     *                 (XJMATR_PROTON(I_LEVL,I_PARM),I_PARM=1,IACPRO)
                 WRITE(LOGAUX,'(4x,''Lev'',i3.3,4x,<IACPRO>(F9.3))')
     *             I_LEVL,
     *            (XJMATR_PROTON(I_LEVL,I_PARM),I_PARM=1,IACPRO)
              END DO
C
              IF (LOGWRI.GT.4)
     *               WRITE(IRESUL,'(''*'',78X,''*'')')
          
          END IF
C
          IF (IFNEUT.EQ.1) THEN
              
              WRITE(LOGAUX,'(/,''<<JACOBMATRX>>    '',
     *                       <IACNEU>(''Par'',i2.2,4X))')(I,I=1,IACNEU)
C
              DO I_LEVL=1,ILEVEL_NEUTRS
C                WRITE(111,'(<IACNEU>(F15.8))')
C     *           (XJMATR_NEUTRS(I_LEVL,I_PARM),I_PARM=1,IACNEU)
                 IF (LOGWRI.GT.4)
     *               WRITE(IRESUL,'(''*'',<IACNEU>(F9.3),T80,''*'')')
     *            (XJMATR_NEUTRS(I_LEVL,I_PARM),I_PARM=1,IACNEU)
                 WRITE(LOGAUX,'(4x,''Lev'',i3.3,4x,<IACNEU>(F9.3))')
     *             I_LEVL,
     *            (XJMATR_NEUTRS(I_LEVL,I_PARM),I_PARM=1,IACNEU)
              END DO
          
          END IF
C
      ELSE ! ( IFDENS.EQ.1 )
      
          WRITE(LOGAUX,'(/,''<<JACOBMATRX>>    '',
     *                       <IACTIV>(''Par'',i2.2,4X))')(I,I=1,IACTIV)
C
          DO I_LEVL=1,LEVNUM
C            WRITE(111,'(<IACTIV>(F15.8))')
C     *           (XJMATR(I_LEVL,I_PARM),I_PARM=1,IACTIV)
             IF (LOGWRI.GT.4)
     *               WRITE(IRESUL,'(''*'',<IACTIV>(F9.3),T80,''*'')')
     *           (XJMATR(I_LEVL,I_PARM),I_PARM=1,IACTIV)
             WRITE(LOGAUX,'(''Lev'',i2.2,1X,a8,<IACTIV>(F9.3))')
     *            I_LEVL,LABTOT_PRINTG(I_LEVL),
     *           (XJMATR(I_LEVL,I_PARM),I_PARM=1,IACTIV)
          END DO
C
      END IF

      IF (LOGWRI.GT.4) WRITE(IRESUL,'(''*'',78X,''*'',/,80(''*''))')
C
C=======================================================================
C     Calculating the inverse covariance matrix (CONVER)
C     [C]^(-1)=[J^T J] --> CONVER(k,l)=\sum_i J(i,k)J(i,l)
C=======================================================================
C
      IF (IFDENS.EQ.0) THEN
C
          IF (IFPROT.EQ.1) THEN
C
          DO I_PARL=1,IACPRO
             DO I_PARK=1,IACPRO
C         
                CONVER_PROTON(I_PARL,I_PARK)=0.
C
                DO I_LEVL=1,ILEVEL_PROTON
C
                   CONVER_PROTON(I_PARL,I_PARK)
     *            =                   
     *             CONVER_PROTON(I_PARL,I_PARK) 
     *            +
     *             XJMATR_PROTON(I_LEVL,I_PARK)
     *            *
     *             XJMATR_PROTON(I_LEVL,I_PARL)
C
               END DO
C
               CONVOL_PROTON(I_PARL,I_PARK)
     *        = 
     *         CONVER_PROTON(I_PARL,I_PARK)
C
            END DO
         END DO
C
         END IF
C  
         IF (IFNEUT.EQ.1) THEN
C
         DO I_PARL=1,IACNEU
            DO I_PARK=1,IACNEU
C         
               CONVER_NEUTRS(I_PARL,I_PARK)=0.
C
               DO I_LEVL=1,ILEVEL_NEUTRS
C
                  CONVER_NEUTRS(I_PARL,I_PARK)=
     *            CONVER_NEUTRS(I_PARL,I_PARK)+
     *            XJMATR_NEUTRS(I_LEVL,I_PARK)*
     *            XJMATR_NEUTRS(I_LEVL,I_PARL)
C
               END DO
C
               CONVOL_NEUTRS(I_PARL,I_PARK)
     *        =
     *         CONVER_NEUTRS(I_PARL,I_PARK)
C
            END DO
         END DO
C
         END IF
C     
C         IF (IFPROT.EQ.0) THEN 
C     
C             DO I_PARL=1,IACPRO
C                WRITE(112,'(25(F15.8))')
C     *               (CONVER_PROTON(I_PARL,I_PARM),I_PARM=1,IACPRO)
C             END DO
C             WRITE(112,'( )')
C         END IF
C     
C         IF (IFNEUT.EQ.0) THEN
C             DO I_PARL=1,IACNEU
C                WRITE(112,'(25(F15.8))')
C     *               (CONVER_NEUTRS(I_PARL,I_PARM),I_PARM=1,IACNEU)
C             END DO
C         END IF
C
C=======================================================================
C        Inverting the CONVER matrix; at the output CONVER stores
C        the covariance matrix
C=======================================================================
C 
C        DGETRF calculates LU decomposition  of the CONVER matrix 
C        for subroutine DGETRI...
C
         IF (IFPROT.EQ.0) GO TO 5
C
         CALL DGETRF(IACPRO,IACPRO,CONVER_PROTON,NDPARS,IPIV,INFOP)
C
C        If INFOP(N) is not equal to 0, the matrix will not be inverted!
C
         IF (INFOP.NE.0) THEN
             IF (LOGWRI.GT.4)
     *           WRITE(IRESUL,'(''* ...in DGETRF... INFOP='',I3,T80,
     *                          ''*'')')INFOP
         END IF
C
         CALL DGETRI(IACPRO,CONVER_PROTON,NDPARS,IPIV,WORK,NDSPEC,INFOP) 
C
         IF (INFOP.NE.0) THEN
             IF (LOGWRI.GT.4)
     *           WRITE(IRESUL,'(''* ...in DGETRF... INFOP='',I3,T80,
     *                          ''*'')')INFOP
         END IF
C
  5      IF (IFNEUT.EQ.0) GO TO 6
C
         CALL DGETRF(IACNEU,IACNEU,CONVER_NEUTRS,NDPARS,IPIV,INFON)
C
         IF (INFON.NE.0) THEN
             IF (LOGWRI.GT.4)
     *           WRITE(IRESUL,'(''* ...in DGETRF... INFON='',I3,T80,
     *                          ''*'')')INFON
         END IF
C
         CALL DGETRI(IACNEU,CONVER_NEUTRS,NDPARS,IPIV,WORK,NDSPEC,INFON) 
C
         IF (INFON.NE.0) THEN
             IF (LOGWRI.GT.4)
     *           WRITE(IRESUL,'(''* ...in DGETRF... INFON='',I3,T80,
     *                          ''*'')')INFON
         END IF
C
  6      CONTINUE
C
C-----------------------------------------------------------------------
C        Printing the covariance matrix...
C-----------------------------------------------------------------------
C
C         IF (IFPROT.EQ.0) THEN 
C             DO I_PARL=1,IACPRO
C                WRITE(112,'(<IACPRO>(F30.8))')
C     *               (CONVER_PROTON(I_PARL,I_PARM),I_PARM=1,IACPRO)
C             END DO
C             WRITE(112,'( )')
C         END IF
C
C         IF (IFNEUT.EQ.0) THEN
C             DO I_PARL=1,IACNEU
C                WRITE(112,'(<IACNEU>(F30.8))')
C     *               (CONVER_NEUTRS(I_PARL,I_PARM),I_PARM=1,IACNEU)
C             END DO
C             WRITE(112,'( )')
C         END IF
C
C-----------------------------------------------------------------------
C        Testing the result...
C-----------------------------------------------------------------------
C
         IERROP=0
	 IERRON=0
C
C        Proton part
C
         IF (IFPROT.EQ.0) GO TO 7
C
         DO L=1,IACPRO
            DO K=1,IACPRO
               TEST_PROTON(L,K)=0.
               DO I=1,IACPRO
                  TEST_PROTON(L,K)=TEST_PROTON(L,K)+
     *                             CONVER_PROTON(L,I)*CONVOL_PROTON(I,K)
               END DO
            END DO
         END DO
C
         DO L=1,IACPRO
            DO K=1,IACPRO
C
               IF (L.NE.K.AND.TEST_PROTON(L,K).GT.EPSTES) THEN
	           IERROP=1
	       END IF
C              
            END DO
C
	    DIFFVAL=TEST_PROTON(L,L)-1.0
C
            IF (DIFFVAL.GT.EPSTES) THEN
	        IERROP=1
	    END IF
C
         END DO
C
C        Neutron part
C
  7      IF (IFNEUT.EQ.0) GO TO 8
C
         DO L=1,IACNEU
            DO K=1,IACNEU
               TEST_NEUTRS(L,K)=0.
               DO I=1,IACTIV
                  TEST_NEUTRS(L,K)=TEST_NEUTRS(L,K)+
     *                             CONVER_NEUTRS(L,I)*CONVOL_NEUTRS(I,K)
               END DO
            END DO
         END DO
C
         DO L=1,IACPRO
            DO K=1,IACPRO
C
               IF (L.NE.K.AND.TEST_NEUTRS(L,K).GT.EPSTES) THEN
	           IERRON=1
	       END IF
C              
            END DO
C
	    DIFFVAL=TEST_NEUTRS(L,L)-1.0
C
            IF (DIFFVAL.GT.EPSTES) THEN
	        IERRON=1
            END IF
C
         END DO
C	 
  8      CONTINUE
C
C        Printing the information about unsuccesful/succesful inversion
C
         IF (LOGWRI.GT.4)
     *       WRITE(IRESUL,'(''*'',78X,''*'',/,
     *                    ''*   Covariance matrix:'',
     *                    T80,''*'',/,''*'',78X,''*'')')
C
         IF (IERROP.EQ.1) THEN
             IF (LOGWRI.GT.4)
     *         WRITE(IRESUL,'(''*   Covariance matrix is incorrectly'',
     *                      '' calculated!!! PROTONS'',T80,''*'')')
         ELSE
             IF (IFPROT.EQ.0) GO TO 9
             DO I_PARL=1,IACPRO
                IF (LOGWRI.GT.4)
     *              WRITE(IRESUL,'(''*'',<IACPRO>(F10.2),T80,''*'')')
     *               (CONVER_PROTON(I_PARL,I_PARM),I_PARM=1,IACPRO)
             END DO
	 END IF
C
  9      CONTINUE
C
         IF (LOGWRI.GT.4) WRITE(IRESUL,'(''*'',78X,''*'')')
C
         IF (IERRON.EQ.1) THEN
             IF (LOGWRI.GT.4)
     *         WRITE(IRESUL,'(''*   Covariance matrix is incorrectly'',
     *                      '' calculated!!! NEUTRONS'',T80,''*'')')
         ELSE
            IF (IFNEUT.EQ.0) GO TO 10
            DO I_PARL=1,IACNEU
               IF (LOGWRI.GT.4)
     *             WRITE(IRESUL,'(''*'',<IACNEU>(F10.2),T80,''*'')')
     *              (CONVER_NEUTRS(I_PARL,I_PARM),I_PARM=1,IACNEU)
            END DO
	 END IF
C
  10     CONTINUE
C
         IF (LOGWRI.GT.4) WRITE(IRESUL,'(''*'',78X,''*'')')
C
         IF (LOGWRI.GT.4) WRITE(IRESUL,'(80(''*''),/,''*'',78X,''*'',/,
     *                ''*   Confidence Intervals for parameters:'',
     *                T80,''*'',/,''*'',78X,''*'')')
C
         IF (IERROP.EQ.0) THEN
C
             IF (IFPROT.EQ.0) GO TO 11
C
             DO I_PARL=1,IACPRO
                IF (LOGWRI.GT.4)
     *              WRITE(IRESUL,'(''*'',F9.2,T80,''*'')')
     *                SQRT(CONVER_PROTON(I_PARL,I_PARL)*DELTAP_PROTON)
             END DO
	 ELSE
	     IF (LOGWRI.GT.4)
     *       WRITE(IRESUL,'(''*  NOT CALCULATED (PROTONS)'',T80,''*'')')
	 END IF
C
  11     CONTINUE
C
         IF (LOGWRI.GT.4) WRITE(IRESUL,'(''*'',78X,''*'')')
C
         IF (IERRON.EQ.0) THEN
            IF (IFNEUT.EQ.0) GO TO 12
            DO I_PARL=1,IACNEU
               IF (LOGWRI.GT.4)
     *             WRITE(IRESUL,'(''*'',F9.2,T80,''*'')')
     *              SQRT(CONVER_NEUTRS(I_PARL,I_PARL)*DELTAP_NEUTRS)
            END DO
	 ELSE
	    IF (LOGWRI.GT.4)
     *          WRITE(IRESUL,'(''*   NOT CALCULATED (NEUTRONS)'',T80,
     *                     ''*'')')
	    
	 END IF
         IF (LOGWRI.GT.4) WRITE(IRESUL,'(''*'',78X,''*'')')
C
  12     CONTINUE
C         
C=======================================================================
C
      ELSE
C
C=======================================================================
C
         DO I_PARL=1,IACTIV
            DO I_PARK=1,IACTIV
C         
               CONVER(I_PARL,I_PARK)=0.
C
               DO I_LEVL=1,LEVNUM
C
                  CONVER(I_PARL,I_PARK) =CONVER(I_PARL,I_PARK)+
     *                                   XJMATR(I_LEVL,I_PARK)*
     *                                   XJMATR(I_LEVL,I_PARL)
C
               END DO
               CONVOL(I_PARL,I_PARK)= CONVER(I_PARL,I_PARK)
            END DO
C
         END DO
C
         WRITE(112,'(''BEFORE GAUSSJ'')')
C
         DO I_PARL=1,IACTIV
            WRITE(112,'(<IACTIV>(F15.8))')
     *           (CONVER(I_PARL,I_PARM),I_PARM=1,IACTIV)
         END DO
C
C-----------------------------------------------------------------------
C        Inverting the CONVER matrix, at the output CONVER stores
C        the covariance matrix
C-----------------------------------------------------------------------
C
         CALL DGETRF(IACTIV,IACTIV,CONVER,NDPARS,IPIV,INFO)
         WRITE(0,*)INFO
C
         IF (INFO.NE.0) THEN
             IF (LOGWRI.GT.4)
     *       WRITE(IRESUL,'(''* ...in DGETRF... INFO='',I3,T80,''*'')')
     *             INFO
         END IF
C
         CALL DGETRI(IACTIV,CONVER,NDPARS,IPIV,WORK,NDSPEC,INFO)
C
         IF (INFO.NE.0) THEN
             IF (LOGWRI.GT.4)
     *       WRITE(IRESUL,'(''* ...in DGETRI... INFO='',I3,T80,''*'')')
     *             INFO
         END IF
         
C
C-----------------------------------------------------------------------
C        Printing the covariance matrix...
C-----------------------------------------------------------------------
C
         WRITE(112,'(''AFTER GAUSSJ'')')
         DO I_PARL=1,IACTIV
            WRITE(112,'(<IACTIV>(F30.8))')
     *           (CONVER(I_PARL,I_PARM),I_PARM=1,IACTIV)
         END DO
C      
C-----------------------------------------------------------------------
C        Testing the result...
C-----------------------------------------------------------------------
C
         I_ERRO=0
         DO L=1,IACTIV
            DO K=1,IACTIV
               TEST(L,K)=0.
               DO I=1,IACTIV
                  TEST(L,K)=TEST(L,K)+CONVER(L,I)*CONVOL(I,K)
               END DO
            END DO
         END DO
C
         DO I_PARL=1,IACTIV
            WRITE(113,'(<IACTIV>(F15.8))')
     *           (TEST(I_PARL,I_PARM),I_PARM=1,IACTIV)
         END DO
      
         DO L=1,IACTIV
            DO K=1,IACTIV
C
               IF (L.NE.K.AND.TEST(L,K).GT.EPSTES) THEN
	           I_ERRO=1
	       END IF
C              
            END DO
	    DIFFVAL=TEST(L,L)-1.0
            IF (DIFFVAL.GT.EPSTES) THEN
                I_ERRO=1
	    END IF
         END DO
C
C        Printing the results...
C
         IF (I_ERRO.EQ.1) THEN
             IF (LOGWRI.GT.4) THEN
                 WRITE(IRESUL,'(''*'',78X,''*'',/,
     *                     ''*   Covariance matrix is incorrectly'',
     *                     '' calculated!!!'',T80,''*'')')
                 WRITE(IRESUL,'(''*'',78X,''*'')')
	         WRITE(IRESUL,'(80(''*''),/,''*'',78X,''*'',/,
     *                        ''*   Confidence intervals '',
     *                        ''NOT CALCULATED'',T80,''*'')')
                 WRITE(IRESUL,'(''*'',78X,''*'')')
             END IF
         ELSE
             IF (LOGWRI.GT.4) THEN
	         WRITE(IRESUL,'(''*'',78X,''*'',/,
     *                    ''*   Covariance matrix:'',
     *                    T80,''*'',/,''*'',78X,''*'')')
                 DO I_PARL=1,IACTIV
                    WRITE(IRESUL,'(''*'',<IACTIV>(F9.2),T80,''*'')')
     *              (CONVER(I_PARL,I_PARM),I_PARM=1,IACTIV)
                 END DO
                 WRITE(IRESUL,'(''*'',78X,''*'')')
                 WRITE(IRESUL,'(80(''*''),/,''*'',78X,''*'',/,
     *                    ''*   Confidence Intervals for parameters:'',
     *                    T80,''*'',/,''*'',78X,''*'')')
	         DO I_PARL=1,IACTIV
                    WRITE(IRESUL,'(''*'',F9.2,T80,''*'')')
     *              SQRT(CONVER(I_PARL,I_PARL)*DELTAP)
                 END DO
                 WRITE(IRESUL,'(''*'',78X,''*'')')
             END IF
         END IF
C
      END IF
C
C
C
C
C
C
C
C
C
C
C
C=======================================================================
C
C     Now we test the J matrix using the Singular Value Decomposition 
C     routine from LAPACK.  The results are stored in SINGLE, UMATRI 
C     and VTMATR.  We chose  to call  the routines  with 'A' option, 
C     which gives full V and U matrix. 
C
C     REMARK: If one chooses to change this option, size declaration 
C     of these matrices has to be changed.
C
C=======================================================================
C
C     Choice of the method of storing the U and V - 'FULL'
C
*      INFO_P=0
*      INFO_N=0
*      JOBU ='A'
*      JOBVT='A'
*C
*      IF (IFDENS.EQ.0) THEN
*C
*         IF (IFPROT.EQ.0) GO TO 13
*C
*         NUMROW=LEVPRO+1
*C
*         CALL DGESVD(  JOBU, JOBVT,NUMROW,IACPRO,XJMATR_PROTON,NDIM_M,
*     *        NDPARS,SINGLE_PROTON,UMATRI_PROTON,NDIM_M,VTMATR_PROTON,
*     *        NDPARS, WORK ,NDSPEC, INFO_P)
*C
*  13     IF (IFNEUT.EQ.0) GO TO 14
*C
*         NUMROW=LEVNEU+1
*C
*         CALL DGESVD(  JOBU, JOBVT,NUMROW,IACNEU,XJMATR_NEUTRS,NDIM_M,
*     *        NDPARS,SINGLE_NEUTRS,UMATRI_NEUTRS,NDIM_M,VTMATR_NEUTRS,
*     *        NDPARS, WORK ,NDSPEC, INFO_N)
*C
*C If something goes wrong INFO values are not 0, and we print them
*C
*  14     CONTINUE
*C
*         IF (INFO_P.NE.0) THEN
*             WRITE(IRESUL,'(''* ...in DGESVD... INFO_P='',T80,''*'')')
*     *             INFO_P
*         END IF
*         IF (INFO_N.NE.0) THEN
*             WRITE(IRESUL,'(''* ...in DGESVD... INFO_N='',T80,''*'')')
*     *             INFO_N
*         END IF
*C
*C Printing the results to IRESUL file...
*C
*         WRITE(IRESUL,'(80(''*''))')
*C
*         IF (IFPROT.EQ.0) GO TO 15
*C
*         WRITE(IRESUL,'(''*'',T80,''*'')')
*         WRITE(IRESUL,'(''*   Proton  singular values:'',T80,''*'')')
*         WRITE(IRESUL,'(''*'',T80,''*'')')
*         DO INUM=1,IACPRO
*            WRITE(IRESUL,'(''*'',F11.6,T80,''*'')')SINGLE_PROTON(INUM)
*	 END DO
*C
* 15      IF (IFNEUT.EQ.0) GO TO 16
*C
*         WRITE(IRESUL,'(''*'',T80,''*'')')
*         WRITE(IRESUL,'(''*   Neutron singular values:'',T80,''*'')')
*         WRITE(IRESUL,'(''*'',T80,''*'')')
*	 DO INUM=1,IACNEU
*            WRITE(IRESUL,'(''*'',F11.6,T80,''*'')')SINGLE_NEUTRS(INUM)
*         END DO
*C
*  16     IF (IFPROT.EQ.0) GO TO 17
*C
*         WRITE(IRESUL,'(''*'',T80,''*'')')
*         WRITE(IRESUL,'(''*   Proton  VT matrix:'',T80,''*'')')
*         WRITE(IRESUL,'(''*'',T80,''*'')')
*         DO INUL=1,IACPRO
*            WRITE(IRESUL,'(''*'',T5,<IACPRO>F9.4,T80,''*'')')
*     *           (VTMATR_PROTON(INUL,INUM),INUM=1,IACPRO)
*         END DO
*C
*  17     IF (IFNEUT.EQ.0) GO TO 18
*C
*         WRITE(IRESUL,'(''*'',T80,''*'')')
*         WRITE(IRESUL,'(''*   Neutron VT matrix:'',T80,''*'')')
*         WRITE(IRESUL,'(''*'',T80,''*'')')
*         DO INUL=1,IACNEU
*            WRITE(IRESUL,'(''*'',T5,<IACNEU>F9.4,T80,''*'')')
*     *           (VTMATR_NEUTRS(INUL,INUM),INUM=1,IACNEU)
*         END DO
*C
*  18     IF (IFPROT.EQ.0) GO TO 19
*C
*         WRITE(IRESUL,'(''*'',T80,''*'')')
*         WRITE(IRESUL,'(''*   Proton  U matrix:'',T80,''*'')')
*         WRITE(IRESUL,'(''*'',T80,''*'')')
*C      
*         DO INUL=1,LEVPRO+1
*            WRITE(IRESUL,'(''*'',T2,<LEVPRO+3>F9.4,T80,''*'')')
*     *           (UMATRI_PROTON(INUL,INUM),INUM=1,LEVPRO+1)
*         END DO
*C
*  19     IF (IFNEUT.EQ.0) GO TO 20
*C
*         WRITE(IRESUL,'(''*'',T80,''*'')')
*         WRITE(IRESUL,'(''*   Neutron U matrix:'',T80,''*'')')
*         WRITE(IRESUL,'(''*'',T80,''*'')')
*         DO INUL=1,LEVNEU+1
*            WRITE(IRESUL,'(''*'',T2,<LEVNEU+3>F9.4,T80,''*'')')
*     *           (UMATRI_NEUTRS(INUL,INUM),INUM=1,LEVNEU+1)
*         END DO
*         WRITE(IRESUL,'(''*'',T80,''*'',/,80(''*''))')
*C
*C Calculating new set of uncorrelated parameters for protons..
*C both stored in PARNEW()
*C
*  20     IF (IFPROT.EQ.0) GO TO 21
*C
*         DO JPARNB=1,IACPRO
*	    SUM_VS=0.
*	    DO LPARNB=1,IACPRO
*	       SUMUWF=0.
*	       DO KLEVNB=1,LEVPRO+1
*	          SUMUWF=SUMUWF+UMATRI_PROTON(KLEVNB,LPARNB)
*     *			       *WEISQR_PROTON(KLEVNB)
*     *                         *(FUNEXP(KLEVNB)-FUNCT0(KLEVNB))
*	       END DO
*               IF (SINGLE_PROTON(LPARNB).GT.SVDCUT) THEN
*	           SUM_VS=SUM_VS+VTMATR_PROTON(LPARNB,JPARNB)
*     *                          *(1./SINGLE_PROTON(LPARNB))
*     *                          *SUMUWF
*C	       PRINT*,'SUMUWF',SUMUWF
*C	       PRINT*,'     v',VTMATR_PROTON(LPARNB,JPARNB)
*C	       PRINT*,'  S^-1',1./SINGLE_PROTON(LPARNB)
*               ELSE
*	           SUM_VS=SUM_VS
*	       END IF
*	    END DO
*C               PRINT*,'SUM_VS=',SUM_VS
*	    PARNEW(JPARNB)=ARGPAR(JPARNB)+SUM_VS
*C	    PRINT*,JPARNB,ARGPAR(JPARNB),PARNEW(JPARNB)
*	 END DO
*C      
*C Calculating new set of uncorrelated parameters for neutrons..
*C
*  21     IF (IFNEUT.EQ.0) GO TO 22
*C
*         DO JPARNB=1,IACNEU
*	    SUM_VS=0.
*	    DO LPARNB=1,IACNEU
*	       SUMUWF=0.
*	       DO KLEVNB=LEVPRO+2,LEVNUM
*	          SUMUWF=SUMUWF+UMATRI_NEUTRS(KLEVNB-LEVPRO-1,LPARNB)
*     *			       *WEISQR_NEUTRS(KLEVNB)
*     *                         *(FUNEXP(KLEVNB)-FUNCT0(KLEVNB))
*	       END DO
*               IF (SINGLE_NEUTRS(LPARNB).GT.SVDCUT) THEN
*	           SUM_VS=SUM_VS+VTMATR_NEUTRS(LPARNB,JPARNB)
*     *                          *(1./SINGLE_NEUTRS(LPARNB))
*     *                          *SUMUWF
*C	       PRINT*,'SUMUWF',SUMUWF
*C	       PRINT*,'     v',VTMATR_NEUTRS(LPARNB,JPARNB)
*C	       PRINT*,'  S^-1',1./SINGLE_neutrs(LPARNB)
*               ELSE
*	           SUM_VS=SUM_VS
*	       END IF
*	    END DO
*C            PRINT*,'SUM_VS=',SUM_VS
*	    PARNEW(JPARNB+IACPRO)=ARGPAR(JPARNB+IACPRO)+SUM_VS
*C	    PRINT*,JPARNB,ARGPAR(JPARNB+IACPRO),PARNEW(JPARNB+IACPRO)
*	 END DO
*C
*  22     CONTINUE
*C-----------------------------------------------------------------------
*      ELSE
*C-----------------------------------------------------------------------
*         NUMROW=LEVNUM
*C
*         CALL DGESVD(  JOBU, JOBVT,NUMROW,IACTIV,XJMATR,NDIM_M,
*     *        NDPARS,SINGLE,UMATRI,NDIM_M,VTMATR,NDPARS, WORK ,NDSPEC,
*     *                                                          INFO )
*C
*         IF (INFO.NE.0) THEN
*             WRITE(IRESUL,'(''* ...in DGESVD... INFO='',I3,T80,''*'')')
*     *             INFO
*         END IF
*C
*C Printing the results to IRESUL file...
*C
*         WRITE(IRESUL,'(80(''*''))')
*         WRITE(IRESUL,'(''*'',T80,''*'')')
*         WRITE(IRESUL,'(''*   Singular values:'',T80,''*'')')
*         WRITE(IRESUL,'(''*'',T80,''*'')')
*         DO INUM=1,IACTIV
*            WRITE(IRESUL,'(''*'',F11.6,T80,''*'')')SINGLE(INUM)
*	 END DO
*         WRITE(IRESUL,'(''*'',T80,''*'')')
*         WRITE(IRESUL,'(''*   VT matrix:'',T80,''*'')')
*         WRITE(IRESUL,'(''*'',T80,''*'')')
*         DO INUL=1,IACTIV
*            WRITE(IRESUL,'(''*'',T5,<IACTIV>F9.4,T80,''*'')')
*     *           (VTMATR(INUL,INUM),INUM=1,IACTIV)
*         END DO
*         WRITE(IRESUL,'(''*'',T80,''*'')')
*         WRITE(IRESUL,'(''*   U matrix:'',T80,''*'')')
*         WRITE(IRESUL,'(''*'',T80,''*'')')
*C      
*         DO INUL=1,LEVNUM
*            WRITE(IRESUL,'(''*'',T2,<LEVNUM>F9.4,T80,''*'')')
*     *           (UMATRI_PROTON(INUL,INUM),INUM=1,LEVNUM)
*         END DO
*         WRITE(IRESUL,'(''*'',T80,''*'',/,80(''*''))')
*C
*         DO JPARNB=1,IACTIV
*	    SUM_VS=0.
*	    DO LPARNB=1,IACTIV
*	       SUMUWF=0.
*	       DO KLEVNB=1,LEVNUM
*	          SUMUWF=SUMUWF+UMATRI(KLEVNB,LPARNB)
*     *			       *WEISQR(KLEVNB)
*     *                         *(FUNEXP(KLEVNB)-FUNCT0(KLEVNB))
*	       END DO
*               IF (SINGLE(LPARNB).GT.SVDCUT) THEN
*	           SUM_VS=SUM_VS+VTMATR(LPARNB,JPARNB)
*     *                          *(1./SINGLE(LPARNB))
*     *                          *SUMUWF
*               ELSE
*	           SUM_VS=SUM_VS
*	       END IF
*	    END DO
*	    PARNEW(JPARNB)=ARGPAR(JPARNB)+SUM_VS
*	 END DO
*      END IF
*C
*C Printing old & new parameters to the IRESUL file
*C
*      WRITE(IRESUL,'(''*'',T80,''*'',/,''*   New set of parameters'',
*     *               '' from the SVD'',T80,''*'',/,''*'',T80,''*'')')
*      WRITE(IRESUL,'(''*'',T12,''OLD'',17X,''NEW'',T80,''*'')')
*      WRITE(IRESUL,'(''*'',T80,''*'')')
*      DO I=1,IACTIV
*         WRITE(IRESUL,'(''*'',2F20.8,T80,''*'')') ARGPAR(I),PARNEW(I)
*      END DO
*      WRITE(IRESUL,'(''*'',T80,''*'',/,80(''*''))')

*C-----------------------------------------------------------------------
*C  Again calculating the covariance matrix, this time using the SVD 
*C  results with singular values cut at the user defined level SVDCUT
*C
*C                        (J^TJ)^-1=VS^-2V^T
*C-----------------------------------------------------------------------
*      WRITE(IRESUL,'(''*'',78X,''*'',/,
*     *               ''*   Covariance matrix:'',
*     *           T80,''*'',/,''*'',78X,''*'')')

*      IF (IFDENS.EQ.0) THEN

*         DO LPARNB=1,IACPRO
*            DO KPARNB=1,IACPRO
*	       SUMSUM=0.
*	       DO IPARNB=1,IACPRO
*	          IF (SINGLE_PROTON(IPARNB).GT.SVDCUT) THEN
*	             SUMSUM=SUMSUM+VTMATR_PROTON(IPARNB,LPARNB)
*     *                            *(SINGLE_PROTON(IPARNB)**(-2))
*     *                            *VTMATR_PROTON(IPARNB,KPARNB)
*                  ELSE
*	             SUMSUM=SUMSUM
*                  END IF
*	       END DO
*	       YJTJMX_PROTON(LPARNB,KPARNB)=SUMSUM
*            END DO
*            WRITE(IRESUL,'(''*'',<IACPRO>(F10.4),T80,''*'')')
*     *                    (YJTJMX_PROTON(LPARNB,I),I=1,IACPRO)	 
*         END DO
*         WRITE(IRESUL,'(''*'',T80,''*'')')
*         DO LPARNB=1,IACNEU
*            DO KPARNB=1,IACNEU
*	       SUMSUM=0.
*	       DO IPARNB=1,IACNEU
*	          IF (SINGLE_NEUTRS(IPARNB).GT.SVDCUT) THEN
*	             SUMSUM=SUMSUM+VTMATR_NEUTRS(IPARNB,LPARNB)
*     *                            *(SINGLE_NEUTRS(IPARNB)**(-2))
*     *                            *VTMATR_NEUTRS(IPARNB,KPARNB)
*                  ELSE
*	             SUMSUM=SUMSUM
*	          END IF
*	       END DO
*               YJTJMX_NEUTRS(LPARNB,KPARNB)=SUMSUM
*            END DO
*            WRITE(IRESUL,'(''*'',<IACNEU>(F10.4),T80,''*'')')
*     *                    (YJTJMX_NEUTRS(LPARNB,I),I=1,IACNEU)	 
*         END DO
*C
*         WRITE(IRESUL,'(''*'',T80,''*'')')
*         WRITE(IRESUL,'(80(''*''),/,''*'',78X,''*'',/,
*     *                ''*   Confidence Intervals for parameters:'',
*     *                T80,''*'',/,''*'',78X,''*'')')
*         DO I_PARL=1,IACPRO
*            WRITE(IRESUL,'(''*'',F9.3,T80,''*'')')
*     *           SQRT(YJTJMX_PROTON(I_PARL,I_PARL)*DELTAP_PROTON)
*         END DO
*         WRITE(IRESUL,'(''*'',T80,''*'')')
*         DO I_PARL=1,IACNEU
*            WRITE(IRESUL,'(''*'',F9.3,T80,''*'')')
*     *           SQRT(YJTJMX_NEUTRS(I_PARL,I_PARL)*DELTAP_NEUTRS)
*         END DO
*         WRITE(IRESUL,'(''*'',T80,''*'',/,80(''*''))')

*      ELSE
*         DO LPARNB=1,IACTIV
*            DO KPARNB=1,IACTIV
*	       SUMSUM=0.
*	       DO IPARNB=1,IACTIV
*	          IF (SINGLE(IPARNB).GT.SVDCUT) THEN
*	             SUMSUM=SUMSUM+VTMATR(IPARNB,LPARNB)
*     *                            *(SINGLE(IPARNB)**(-2))
*     *                            *VTMATR(IPARNB,KPARNB)
*                  ELSE
*	             SUMSUM=SUMSUM
*                  END IF
*	       END DO
*	       YJTJMX(LPARNB,KPARNB)=SUMSUM
*            END DO
*            WRITE(IRESUL,'(''*'',<IACTIV>(F10.4),T80,''*'')')
*     *                    (YJTJMX(LPARNB,I),I=1,IACTIV)	 
*         END DO
*         WRITE(IRESUL,'(''*'',T80,''*'')')

*         WRITE(IRESUL,'(80(''*''),/,''*'',78X,''*'',/,
*     *                ''*   Confidence Intervals for parameters:'',
*     *                T80,''*'',/,''*'',78X,''*'')')
*         DO I_PARL=1,IACTIV
*            WRITE(IRESUL,'(''*'',F9.3,T80,''*'')')
*     *           SQRT(YJTJMX(I_PARL,I_PARL)*DELTAP)
*         END DO
*         WRITE(IRESUL,'(''*'',T80,''*'',/,80(''*''))')

*      END IF
C
C=======================================================================
C
  100 CONTINUE
C
      IF (LOGWRI.GT.0) THEN
          WRITE(LOGFIL,'(/,9X,''Exiting  JMATRX'')')
          WRITE(LOGAUX,'(/,''... leaving JMATRX'',/)')
      END IF
C
      IF (ISCREN.GT.0) WRITE(LSCREN,'(/,''... leaving JMATRX'',/)')
C
C=======================================================================
C     
CID      CALL CPUTIM('JMATRX',0)
C
C=======================================================================
C 
      RETURN
      END
C
C=======================================================================
C=======================================================================
C
      SUBROUTINE TCRITC(CRITIC)
C
      PARAMETER (ITERNB=50)
      COMMON
     *     /BETDAT/ APARAM,BPARAM,PROBAB,FREDEG
      COMMON
     *     /XARGUM/ X_ARGU
      EXTERNAL
     *      BISEC,BETFUN,TFUNEV
      DATA
     *    AINTER /0.0/,
     *    BINTER /1.0/,
     *    ALIMTT /0.0/,
     *    BLIMTT /700.0/,
     *    EPSLOO /1.0E-13/
C
C=======================================================================
C
C  This routine calculates the critical t values
C  (student's distribution) for given ONE-SIDED (!)probability:
C  Pr(T<A)=1-p with A=t(p,n) where n is the number of 
C  degrees of freedom. We also have Pr(-A<T<A)=1-2p.
C  The confidence interval is given by: D=A*Sn/sqrt(n)
C
C  The probability function is given by incomplete beta functions:
C  P=Ix(a,b) with a=b=n/2 and x=(t+sqrt(t^2+n))/(2*sqrt(t^2+n))
C  
C
C  Meaning of the input parameters (common BETDAT):
C
C  PROBAB        - given propability (1-p)
C
C  FREDEG        - number of degrees of freedom (n)
C
C  APARAM,BPARAM - parameters for the evaluation of the inc. beta func.
C                  (a,b)
C  TCRITV        - calculated value of the critical t
C   
C      TT=1.37638191977842
c      TT=1.37638191977869
C      FREDEG=1.
C      Y=(TT+SQRT(TT**2+FREDEG))/(2.*SQRT(TT**2+FREDEG))
c      PRINT*,Y
      APARAM=0.5*FREDEG
      BPARAM=0.5*FREDEG
c      PROBAB=0.80
c     WRITE(666,'(3F12.5)')FREDEG,PROBAB
      
      CALL BISEC(AINTER,BINTER,XARGUM,BETFUN,EPSLOO,ITERNB,IERROR)
c     WRITE(666,'(F11.5,I3)')XARGUM,IERROR
      X_ARGU=XARGUM
      CALL BISEC(ALIMTT,BLIMTT,TARGUM,TFUNEV,EPSLOO,ITERNB,IERROR)       
c     WRITE(666,'(F11.5,I3)')TARGUM,IERROR
C
      CRITIC=TARGUM
C
C=======================================================================
C
      RETURN
      END 
C
C=======================================================================
C=======================================================================
C
      FUNCTION BETFUN(X_ARGU)
C      
      COMMON
     *      /BETDAT/ APARAM,BPARAM,PROBAB,FREDEG
C
C=======================================================================
C     This routine prepares the function for the bisection routine
C                  f(x)=BETAI(x)-P
C=======================================================================
C
      BETFUN=BETAI(APARAM,BPARAM,X_ARGU)-PROBAB
C
C=======================================================================
C
      RETURN
      END 
C
C=======================================================================
C=======================================================================
C
      FUNCTION TFUNEV(T_ARGU)
C      
      COMMON
     *       /BETDAT/ APARAM,BPARAM,PROBAB,FREDEG
      COMMON
     *       /XARGUM/ X_ARGU
C
C=======================================================================
C     This routine prepares the function for the bisection routine
C                  f(x)=sqrt(t^2+n)*(2*x-1)-t
C=======================================================================
C      
      TFUNEV=SQRT(T_ARGU**2+FREDEG)*(2.*X_ARGU-1.)-T_ARGU
C
C=======================================================================
C 
      RETURN     
      END
C
C=======================================================================
C
      SUBROUTINE BISEC(A,B,X,F,EPS,ITER,IER)
      EXTERNAL
     *         F
C
C=======================================================================
C         SOLVES THE EQUATION F(X)=0 USING A BISECTION METHOD
C=======================================================================
C
      XA=A
      XB=B
      IT=0
C
      FA=F(XA)
      FB=F(XB)
C
      ABFA=ABS(FA)
      X=A
C
      IF (ABFA.LT.EPS) GO TO 10
C
      ABFB=ABS(FB)
      X=B
C
      IF (ABFB.LT.EPS) GO TO 10
      IF (FA*FB.GT.0.0) GO TO 11
C
   1  CONTINUE
C
      IT=IT+1
      X=0.5*(XA+XB)
      FX=F(X)
      ABFF=ABS(FX)
C
      IF (ABFF.LT.EPS) GO TO 10
      IF (IT.GT.ITER) GO TO 20
      IF (FA*FX) 2,10,4
C
   2  CONTINUE
C
C=======================================================================
C         FUNCTION CHANGES SIGN IN THE LEFT SUBINTERVAL
C=======================================================================
C
      XB=X
C
      GO TO 1
C
   4  CONTINUE
C
C=======================================================================
C         FUNCTION CHANGES SIGN IN THE RIGHT SUBINTERVAL
C=======================================================================
C
      XA=X
C
      GO TO 1
C
   10 CONTINUE
C
      IER=0
C
      RETURN
C
   11 CONTINUE
C
      PRINT 30
   30 FORMAT(///,1X,18HBISEC OUT OF RANGE ,///)
C
      STOP 'BISEC'
C
   20 CONTINUE
C
      IER=10
C
      RETURN
      END
C
C=======================================================================
C=======================================================================
C
      FUNCTION betai(a,b,x)
      REAL betai,a,b,x
CU    USES betacf,gammln
      REAL bt,betacf,gammln
      if(x.lt.0..or.x.gt.1.)pause 'bad argument x in betai'
      if(x.eq.0..or.x.eq.1.)then
        bt=0.
      else
        bt=exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.-x))
      endif
      if(x.lt.(a+1.)/(a+b+2.))then
        betai=bt*betacf(a,b,x)/a
        return
      else
        betai=1.-bt*betacf(b,a,1.-x)/b
        return
      endif
      END
C
C=======================================================================
C=======================================================================
C
      FUNCTION gammln(xx)
      REAL gammln,xx
      INTEGER j
CID      DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
      DIMENSION cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     *24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     *-.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do 11 j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
11    continue
      gammln=tmp+log(stp*ser/x)
      return
      END
C
C=======================================================================
C=======================================================================
C
      FUNCTION betacf(a,b,x)
      INTEGER MAXIT
      REAL betacf,a,b,x,EPS,FPMIN
      PARAMETER (MAXIT=100,EPS=3.e-7,FPMIN=1.e-30)
      INTEGER m,m2
      REAL aa,c,d,del,h,qab,qam,qap
      qab=a+b
      qap=a+1.
      qam=a-1.
      c=1.
      d=1.-qab*x/qap
      if(abs(d).lt.FPMIN)d=FPMIN
      d=1./d
      h=d
      do 11 m=1,MAXIT
        m2=2*m
        aa=m*(b-m)*x/((qam+m2)*(a+m2))
        d=1.+aa*d
        if(abs(d).lt.FPMIN)d=FPMIN
        c=1.+aa/c
        if(abs(c).lt.FPMIN)c=FPMIN
        d=1./d
        h=h*d*c
        aa=-(a+m)*(qab+m)*x/((a+m2)*(qap+m2))
        d=1.+aa*d
        if(abs(d).lt.FPMIN)d=FPMIN
        c=1.+aa/c
        if(abs(c).lt.FPMIN)c=FPMIN
        d=1./d
        del=d*c
        h=h*del
        if(abs(del-1.).lt.EPS)goto 1
11    continue
      pause 'a or b too big, or MAXIT too small in betacf'
1     betacf=h
      return
      END
C
C=======================================================================
C=======================================================================
C

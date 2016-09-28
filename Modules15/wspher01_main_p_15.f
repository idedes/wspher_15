C FILE NAME = wspher01_main_p_15.f ! Keep this symbol:    $ident@string$
C=======================================================================
C=======================================================================
C=======================================================================
C     
      PROGRAM WSPHER15
C      
      EXTERNAL 
     *          FUNMIN
C
      INCLUDE  'MATDIM/MAXFEV.f'
      INCLUDE  'MATDIM/NDFUNC.f'
      INCLUDE  'MATDIM/NDIM_M.f'
      INCLUDE  'MATDIM/NDMESH.f'
      INCLUDE  'MATDIM/NDSPEC.f'
      INCLUDE  'MATDIM/NDGAUS.f'
      INCLUDE  'MATDIM/NDPARS.f'
      INCLUDE  'MATDIM/NDIM_N.f'
      INCLUDE  'MATDIM/NDNUCL.f'
      INCLUDE  'MATDIM/NDLEXP.f'
      INCLUDE  'MATDIM/NDRADZ.f'
      INCLUDE  'MATDIM/NDRADN.f'
      INCLUDE  'MATDIM/NDMAIN.f'
      INCLUDE  'MATDIM/NDLAST.f'
      INCLUDE  'MATDIM/NDIM_P.f'
      INCLUDE  'MATDIM/ND_RHO.f'
      INCLUDE  'MATDIM/N_NOYX.f'
      INCLUDE  'MATDIM/NDTITL.f'
C
      PARAMETER
     *         (NDIM_L=NDIM_N,NDBASE=NDIM_N+1,NDMES2=NDMESH*NDMESH)
      PARAMETER 
     *         (IFMODE=1,NPRINT=5,NDWEIG=4000)
C
      REAL*16
     *          ANODES,AWEIGH,BAUXIL,CAUXIL,EPSLAG,QRNODE,
     *          QRWEIG,QNORNL,QDFLAG,QFLAGN,QPOLYN,QDPOLY,
     *          QPOLYN_DENSIT,QDPOLY_DENSIT
C
      CHARACTER
     *          PRONAM*10,WHATEX*06,LABSYM_PROTON*6,LABSYM_NEUTRS*6
      CHARACTER
     *          FILNAM*256,TITPAR*13,TITMSH*8,TITLES*12,
     *                                        TITLES_EXPORT*12
      CHARACTER
     *          INPSYM*6,TYPCHI*6,TAKCHI*3,NUCSYM*6
      CHARACTER
     *          DIRNAM*256,GENAME*256,STRING*126
      CHARACTER
     *          LABTHE_PROTON*6,LABTHE_NEUTRS*6,
     *          LABEXP_PROTON*6,LABEXP_NEUTRS*6
      CHARACTER
     *          TITAUX*7,AUXLAB*6,LABTEX*11
      CHARACTER
     *          FILNAM_ENDING*3,TAKPAR*3
      CHARACTER
     *          CALLED*6,VERSIO*3
      CHARACTER
     *          TITLES_LATEXS*050,NUCNAM_LATEXS*010,PARPOT_UNITSS*040
C
C=======================================================================
C
C               Basis related, Laguerre functions, Gauss integration
C
      DIMENSION
     *          QPOLYN_DENSIT(1:NDGAUS,0:NDIM_N,0:NDIM_L,0:NDIM_L), 
     *          QDPOLY_DENSIT(1:NDGAUS,0:NDIM_N,0:NDIM_L,0:NDIM_L)
      DIMENSION
     *          QPOLYN(1:NDGAUS,0:NDIM_N,0:NDIM_L), 
     *          QDPOLY(1:NDGAUS,0:NDIM_N,0:NDIM_L)
      DIMENSION
     *          QFLAGN(0:NDIM_N,0:NDIM_L),
     *          QDFLAG(0:NDIM_N,0:NDIM_L)
      DIMENSION
     *          QRNODE(1:NDGAUS,0:NDIM_L),
     *          QRWEIG(1:NDGAUS,0:NDIM_L)
C
      DIMENSION
     *          DENSIT_PROTON(1:NDGAUS,0:NDIM_L),
     *          DENSIT_NEUTRS(1:NDGAUS,0:NDIM_L)
C      
      DIMENSION 
     *          DIAGSC(1:NDPARS),QTRANF(1:NDPARS),
     *          WORKA1(1:NDPARS),WORKA2(1:NDPARS),
     *          WORKA3(1:NDPARS),WORKA4(1:NDFUNC)
      DIMENSION
     *          FJACOB(1:NDFUNC,1:NDPARS)
      DIMENSION 
     *          I_PERM(1:NDPARS)
      DIMENSION 
     *          FUNJAC(1:NDIM_M,1:NDPARS)
      DIMENSION 
     *          IPIVOT(1:NDPARS)
      DIMENSION
     *          INDEXI(1:NDIM_M),MCOUNT(1:NDIM_M),STPMSH(1:NDIM_M),
     *                                            MAXMSH(1:NDIM_P)
      DIMENSION
     *          PARMSH(1:NDIM_M,1:NDMESH),
     *          AXSMAP(1:NDMESH,1:NDMESH)
      DIMENSION
     *          TITPAR(1:NDPARS),TITMSH(1:NDIM_P),TAKCHI(1:NDNUCL),
     *                           TITAUX(1:NDIM_M),TAKPAR(1:NDIM_M)
      DIMENSION
     *          CHITOT_OFMESH(1:NDMESH,1:NDMESH),
     *          CHIPRO_OFMESH(1:NDMESH,1:NDMESH),
     *          CHINEU_OFMESH(1:NDMESH,1:NDMESH)
      DIMENSION
     *          PARPOT_OFMESH(1:NDMESH,1:NDMESH,1:NDPARS),
     *          PARPOT_PRINTI(1:NDPARS,1:NDMES2)
C
      DIMENSION
     *          RFUNUP_PROTON(1:NDGAUS,0:NDIM_N,0:NDIM_L),
     *          RFUNDN_PROTON(1:NDGAUS,0:NDIM_N,0:NDIM_L),
     *          RFUNUP_NEUTRS(1:NDGAUS,0:NDIM_N,0:NDIM_L),
     *          RFUNDN_NEUTRS(1:NDGAUS,0:NDIM_N,0:NDIM_L)
      DIMENSION
     *          LABORD_AUXILP(0:NDIM_N,0:NDIM_L,0:2*NDIM_L+1),
     *          LABORD_AUXILN(0:NDIM_N,0:NDIM_L,0:2*NDIM_L+1)
      DIMENSION 
     *          POLYNX(1:NDGAUS,0:NDIM_N,0:NDIM_L)
      DIMENSION 
     *          QNORNL(0:NDIM_N,0:NDIM_L)
      DIMENSION
     *          RAD_UP(0:NDIM_L,0:NDIM_N),
     *          RAD_DN(0:NDIM_L,0:NDIM_N)
      DIMENSION
     *          HSORUP(1:NDBASE,1:NDBASE),HSORDN(1:NDBASE,1:NDBASE)
      DIMENSION
     *          SPENUP(1:NDBASE),SPENDN(1:NDBASE),AUXVEC(1:NDBASE)     
      DIMENSION
     *          ANODES(1:NDGAUS),AWEIGH(1:NDGAUS),BAUXIL(1:NDGAUS),
     *                                            CAUXIL(1:NDGAUS)
      DIMENSION
     *          ARGMNT(1:NDPARS)
      DIMENSION
     *          RMSPRO_WEIDEP(1:NDWEIG),
     *          RMSNEU_WEIDEP(1:NDWEIG),
     *          RMSPRO_WEINUC(1:NDWEIG,1:NDNUCL),
     *          RMSNEU_WEINUC(1:NDWEIG,1:NDNUCL),
     *          RMSVAL_RATPRO(1:NDWEIG,1:NDNUCL),
     *          RMSVAL_RATNEU(1:NDWEIG,1:NDNUCL),
     *          RMSGLO_RATPRO(1:NDWEIG),
     *          RMSGLO_RATNEU(1:NDWEIG),
     *          IND_V1(1:NDWEIG),
     *          IND_V2(1:NDWEIG),
     *          IND_V3(1:NDWEIG),
     *          IND_V4(1:NDWEIG),
     *          IND_V5(1:NDWEIG)
      DIMENSION
     *          WEIG1P(1:NDWEIG),
     *          WEIG2P(1:NDWEIG),
     *          WEIG3P(1:NDWEIG),
     *          WEIG6P(1:NDWEIG),
     *          WEIG8P(1:NDWEIG),
     *          WEIG1N(1:NDWEIG),
     *          WEIG2N(1:NDWEIG),
     *          WEIG3N(1:NDWEIG),
     *          WEIG6N(1:NDWEIG),
     *          WEIG8N(1:NDWEIG)
      DIMENSION
     *          IASTER(1:NDWEIG),JASTER(1:NDWEIG)
C
      COMMON
     *       /TITEXP/ TITLES_EXPORT(1:NDTITL)
      COMMON
     *       /POTPOT/ PARPOT(1:NDPARS)
      COMMON
     *       /ARGINI/ XINITS(1:NDIM_P)
C @@@ IRENE WRONG ALIGNMENT =>> PLEASE SPLIT INTO TWO COMMONS
      COMMON
     *       /MESHIN/ I_MESH(1:NDIM_P),
     *                XMIN_I(1:NDIM_P),
     *                XMAX_I(1:NDIM_P),
     *                MAXPAR(1:NDIM_P)
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
C     
      COMMON
     *       /LEVORD/ LABORD_PROTON(1:NDNUCL,
     *                              0:NDIM_N,0:NDIM_L,0:2*NDIM_L+1),
     *                LABORD_NEUTRS(1:NDNUCL,
     *                              0:NDIM_N,0:NDIM_L,0:2*NDIM_L+1)
      COMMON
     *       /LEVSYM/ LABSYM_PROTON(1:NDNUCL,
     *                              0:NDIM_N,0:NDIM_L,0:2*NDIM_L+1),
     *                LABSYM_NEUTRS(1:NDNUCL,
     *                              0:NDIM_N,0:NDIM_L,0:2*NDIM_L+1)
      COMMON
     *       /DIMSPH/ LDSING_PROTON,LEVTHE_PROTON,
     *                              LEVEXP_PROTON(1:NDNUCL),
     *                LDSING_NEUTRS,LEVTHE_NEUTRS,
     *                              LEVEXP_NEUTRS(1:NDNUCL)
      COMMON
     *       /EKINET/ ENEKIN_PROTON(1:NDNUCL,1:NDBASE,
     *                                       1:NDBASE,0:NDIM_L),
     *                ENEKIN_NEUTRS(1:NDNUCL,1:NDBASE,
     *                                       1:NDBASE,0:NDIM_L)
      COMMON
     *       /CFACTO/ CMATUP_PROTON(0:NDIM_L,1:NDBASE,1:NDBASE),
     *                CMATDN_PROTON(0:NDIM_L,1:NDBASE,1:NDBASE),
     *
     *                CMATUP_NEUTRS(0:NDIM_L,1:NDBASE,1:NDBASE),
     *                CMATDN_NEUTRS(0:NDIM_L,1:NDBASE,1:NDBASE)
      COMMON
     *       /POLPOL/ XPOLYN(1:NDGAUS,0:NDIM_N,0:NDIM_L), 
     *                XDPOLY(1:NDGAUS,0:NDIM_N,0:NDIM_L)
      COMMON
     *       /POLDEN/ XPOLYN_DENSIT(1:NDGAUS,0:NDIM_N,0:NDIM_L,
     *                                                0:NDIM_L), 
     *                XDPOLY_DENSIT(1:NDGAUS,0:NDIM_N,0:NDIM_L,
     *                                                0:NDIM_L)
      COMMON
     *       /ENUPDN/ ENERUP_PROTON(0:NDIM_L,1:NDBASE),
     *                ENERDN_PROTON(0:NDIM_L,1:NDBASE),
     *
     *                ENERUP_NEUTRS(0:NDIM_L,1:NDBASE),
     *                ENERDN_NEUTRS(0:NDIM_L,1:NDBASE)
      COMMON
     *       /NMLCNT/ XNORNL(0:NDIM_N,0:NDIM_L)
      COMMON 
     *       /LAGLAG/ XFLAGN(0:NDIM_N,0:NDIM_L),
     *                XDFLAG(0:NDIM_N,0:NDIM_L)
      COMMON
     *       /NODWEI/ XRNODE(1:NDGAUS,0:NDIM_L),
     *                XRWEIG(1:NDGAUS,0:NDIM_L)
      COMMON
     *       /VLIMIT/ VMISTR(1:NDPARS),
     *                VMIMIN(1:NDPARS),
     *                VMIMAX(1:NDPARS),
     *                VMISTP(1:NDPARS)
      COMMON
     *       /SEARCH/ IFTAKE(1:NDPARS)
C 
      COMMON
     *       /V2COEF/ VCOEUP_PROTON(0:NDIM_L,1:NDBASE),
     *                VCOEDN_PROTON(0:NDIM_L,1:NDBASE),
     *
     *                VCOEUP_NEUTRS(0:NDIM_L,1:NDBASE),
     *                VCOEDN_NEUTRS(0:NDIM_L,1:NDBASE)
      COMMON
     *       /MODIFR/ RMSEXP_PROTON(1:NDNUCL),
     *                RMSEXP_NEUTRS(1:NDNUCL),
     *                RMSTHE_PROTON,RMSTHE_NEUTRS,
     *                RMODIF_PROTON,RMODIF_NEUTRS,
     *                              PUSHIN,DELTAR
      COMMON
     *       /ER_RMS/ RMSERR_PROTON(1:NDNUCL),
     *                RMSERR_NEUTRS(1:NDNUCL)
      COMMON
     *       /FERGPS/ FERMEX_PROTON(1:NDNUCL),GAPEXP_PROTON(1:NDNUCL),
     *                DENSUP_PROTON(1:NDNUCL),DENSDW_PROTON(1:NDNUCL),
     *
     *                FERMEX_NEUTRS(1:NDNUCL),GAPEXP_NEUTRS(1:NDNUCL),
     *                DENSUP_NEUTRS(1:NDNUCL),DENSDW_NEUTRS(1:NDNUCL)
      COMMON
     *       /DATA_V/ AOSCIL_PROTON(1:NDNUCL),
     *                AOSCIL_NEUTRS(1:NDNUCL),
     *                HOMEGA(1:NDNUCL)    
      COMMON
     *       /QUANTN/ NWSSPH_PROTON(1:NDSPEC),
     *                LWSSPH_PROTON(1:NDSPEC),
     *                JWSSPH_PROTON(1:NDSPEC),
     *
     *                NWSSPH_NEUTRS(1:NDSPEC),
     *                LWSSPH_NEUTRS(1:NDSPEC),
     *                JWSSPH_NEUTRS(1:NDSPEC)
      COMMON
     *       /THELAB/ LABTHE_PROTON(1:NDSPEC),
     *                LABTHE_NEUTRS(1:NDSPEC)
      COMMON
     *       /SPECTR/ ENETHE_PROTON(1:NDSPEC),
     *                ENETHE_NEUTRS(1:NDSPEC)
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
     *       /NUCLEU/ IZ_FIX,IN_FIX,ISOSPI    
      COMMON
     *       /CHIEXT/ CHIMI1,CHIMI2_PROTON,CHIMI3_PROTON,
     *                       CHIMI2_NEUTRS,CHIMI3_NEUTRS
      COMMON
     *       /GAUSHE/ NGAUSS,NSHELL_PROTON,NSHELL_NEUTRS
      COMMON
     *       /CHI_PN/ CHISQU_PROTON,CHISQU_NEUTRS
      COMMON
     *       /CHICHO/ VALCHI
     *       /DENSTR/ IFDENS
     *       /TENSOR/ IFTENS
      COMMON
     *       /ENEPRI/ ENEMAX
     *       /INPSEL/ WHATEX
      COMMON
     *       /TENOPT/ ICENTT,ISORBT,ITENSR
     *       /TESFIT/ IFTEST,IFFITS,IFMESH,ISPE_P,IFMOCA
      COMMON
     *       /MONTEC/ IFPSEU,IFPARA,LDMONT,LDBINS
      COMMON
     *       /FEVALS/ IFUNCT_EVALUS
     *       /P_N_LM/ IFPROT,IFNEUT,IFBOTH
     *       /MINITE/ IDEFCN,ITECHI
      COMMON
     *       /CNTCNT/ ICOUNT_DENGRD,ICOUNT_FUNMIN,ICOUNT_HAMMAT,
     *                ICOUNT_EXPTHE,ICOUNT_CHOICE,ICOUNT_DIAMAT
      COMMON 
     *       /CHITAK/ IF_SPE,IF_RAD,IF_GAP,IF_FER,IF_DEN,IF_RHO,IF_INV,
     *                LDCHI2,TYPCHI(1:NDNUCL),ITKCHI(1:NDNUCL)
      COMMON  
     *       /VERSIN/ VERSIO
      COMMON
     *       /PARAMT_TITLES/ TITLES(1:NDTITL),
     *                       TITLES_LATEXS(1:NDTITL),
     *                       NUCNAM_LATEXS(1:NDNUCL)
      COMMON
     *       /PARUNI/ PARPOT_UNITSS(1:NDTITL)
      COMMON
     *       /NUCWEI/ WEINUC_PROTON(1:NDNUCL),
     *                WEINUC_NEUTRS(1:NDNUCL),
     *                WEIRAD_PROTON(1:NDNUCL),
     *                WEIRAD_NEUTRS(1:NDNUCL)
      COMMON
     *       /RMSGLO_NUCLEU/ RMSGLO_TAKPRO,
     *                       RMSGLO_TAKNEU,
     *                       RMSVAL_TAKPRO(1:NDNUCL),
     *                       RMSVAL_TAKNEU(1:NDNUCL)
      COMMON
     *       /RMSREF/ RMSVAL_REFPRO(1:NDNUCL),
     *                RMSVAL_REFNEU(1:NDNUCL),
     *                RMSGLO_REFPRO,
     *                RMSGLO_REFNEU
      COMMON
     *       /RMSIND/ RMSIND_PROTON(1:NDNUCL),
     *                RMSIND_NEUTRS(1:NDNUCL)
C
      COMMON
     *       /PRFILE/ ISCREN,LOGFIL,LOGCHP,LOGCHN,LOGHIP,LOGHIN,
     *                LOGPRO,LOGNEU,LOGWRI,LSCREN,IRESUL,ICONVE,
     *                NOUTPT,LOGAUX,LOGENE,LOGRAD,LOGGRD,LOGMSH,
     *                                                   LOGBIS 
C
      DATA
     *       EPSLAG / 1.0Q-32 /
C
C=======================================================================
C=======================================================================
C=======================================================================
C
      CALL CPUTIM('WSPH15',1)
C
C=======================================================================
C=======================================================================
C=======================================================================
C
                       PRONAM='wspher_15W'
C
C=======================================================================
C
      VERSIO='15W'
C
C=======================================================================
C     
C     Creating the necessary titles
C
      CALL CREATG_TITLES(NDTITL,NDNUCL,TITLES,TITLES_LATEXS,
     *                          PARPOT_UNITSS,NUCNAM_LATEXS)
C
      DO I=1,NDTITL
         TITLES_EXPORT(I)=TITLES(I)
      END DO
C
C=======================================================================
C
C     Creating the necessary folders to store the results.
C @@@ Irene?
C     The indication -p means that if the directory is already 
C     created the machine will not complain, i.e. the folders
C     will not be over-written (but yes the files if inside the
C     code are created again with the same name).
C
C      THE OLD FILES WILL BE OVERWRITTEN? YES, THEY WILL
C
      CALL SYSTEM('mkdir -pv Results')
C
      CALL SYSTEM('mkdir -pv LogFiles/IFDENS-0')
      CALL SYSTEM('mkdir -pv LogFiles/IFDENS-1_IFTENS-0')
      CALL SYSTEM('mkdir -pv LogFiles/IFDENS-1_IFTENS-1')
C
      CALL SYSTEM('mkdir -pv Fitting-Results')
      CALL SYSTEM('mkdir -pv Fitting-Results/Global-Fit')
      CALL SYSTEM('mkdir -pv Fitting-Results/Fortran-Readable')
      CALL SYSTEM('mkdir -pv Fitting-Results/Latex-Readable')
C
      CALL SYSTEM('mkdir -pv EnergyLevels/IFDENS-0/Fitted')
      CALL SYSTEM('mkdir -pv EnergyLevels/IFDENS-0/Predic')
      CALL SYSTEM('mkdir -pv EnergyLevels/IFDENS-0/OneRun')
      CALL SYSTEM('mkdir -pv EnergyLevels/IFDENS-0/WSUniv')
C
      CALL SYSTEM('mkdir -pv EnergyLevels/IFDENS-1_IFTENS-0/Fitted')
      CALL SYSTEM('mkdir -pv EnergyLevels/IFDENS-1_IFTENS-0/Predic')
      CALL SYSTEM('mkdir -pv EnergyLevels/IFDENS-1_IFTENS-0/OneRun')
C
      CALL SYSTEM('mkdir -pv EnergyLevels/IFDENS-1_IFTENS-1/Fitted')
      CALL SYSTEM('mkdir -pv EnergyLevels/IFDENS-1_IFTENS-1/Predic')
      CALL SYSTEM('mkdir -pv EnergyLevels/IFDENS-1_IFTENS-1/OneRun')
C
      CALL SYSTEM('mkdir -pv DensitCurves/IFDENS-0/Fitted')
      CALL SYSTEM('mkdir -pv DensitCurves/IFDENS-0/Predic')
      CALL SYSTEM('mkdir -pv DensitCurves/IFDENS-0/OneRun')
C
      CALL SYSTEM('mkdir -pv DensitCurves/IFDENS-1_IFTENS-0/Fitted')
      CALL SYSTEM('mkdir -pv DensitCurves/IFDENS-1_IFTENS-0/Predic')
      CALL SYSTEM('mkdir -pv DensitCurves/IFDENS-1_IFTENS-0/OneRun')
C
      CALL SYSTEM('mkdir -pv DensitCurves/IFDENS-1_IFTENS-1/Fitted')
      CALL SYSTEM('mkdir -pv DensitCurves/IFDENS-1_IFTENS-1/Predic')
      CALL SYSTEM('mkdir -pv DensitCurves/IFDENS-1_IFTENS-1/OneRun')
C
      CALL SYSTEM('mkdir -pv LogFMesh/IFDENS-0')
      CALL SYSTEM('mkdir -pv LogFMesh/IFDENS-1_IFTENS-0')
      CALL SYSTEM('mkdir -pv LogFMesh/IFDENS-1_IFTENS-1')
C
      CALL SYSTEM('mkdir -pv MeshTables/IFDENS-0')
      CALL SYSTEM('mkdir -pv MeshTables/IFDENS-1_IFTENS-0')
      CALL SYSTEM('mkdir -pv MeshTables/IFDENS-1_IFTENS-1')
C
      CALL SYSTEM('mkdir -pv Chi2Maps/IFDENS-0')
      CALL SYSTEM('mkdir -pv Chi2Maps/IFDENS-1_IFTENS-0')
      CALL SYSTEM('mkdir -pv Chi2Maps/IFDENS-1_IFTENS-1')
C
      CALL SYSTEM('mkdir -pv SPECurves/IFDENS-0')
      CALL SYSTEM('mkdir -pv SPECurves/IFDENS-1_IFTENS-0')
      CALL SYSTEM('mkdir -pv SPECurves/IFDENS-1_IFTENS-1')
C
      CALL SYSTEM('mkdir -pv LogFMont/IFDENS-0')
      CALL SYSTEM('mkdir -pv LogFMont/IFDENS-1_IFTENS-0')
      CALL SYSTEM('mkdir -pv LogFMont/IFDENS-1_IFTENS-1')
C
      CALL SYSTEM('mkdir -pv LogMonte/IFDENS-0')
      CALL SYSTEM('mkdir -pv LogMonte/IFDENS-1_IFTENS-0')
      CALL SYSTEM('mkdir -pv LogMonte/IFDENS-1_IFTENS-1')
C
      CALL SYSTEM('mkdir -pv SPEMontC/IFDENS-0')
      CALL SYSTEM('mkdir -pv SPEMontC/IFDENS-1_IFTENS-0')
      CALL SYSTEM('mkdir -pv SPEMontC/IFDENS-1_IFTENS-1')
C
C=======================================================================
C
C     Defining the output unit numbers for the log-files
C
      LOGCHP=22 ! Log file to store chi^2 for protons
      LOGCHN=23 ! Log file to store chi^2 for neutrons
C
      LOGPRO=82
      LOGNEU=83
C
      IRESUL=10
      ICONVE=20

      NOUTPT=6
      N_OUTP=6
      LSCREN=0
C
      IORTHO=1
C
C=======================================================================
C
      LOGBIS=41
C
      LOGAUX=50
      LOGENE=550
      LOGRAD=551
C
C=======================================================================
C
C                 Opening the central logfile called "total" =>> _t
      LOGFIL=21
      LOGWRI=01 ! This will be overwritten in SUBROUTINE NAMELI
C
      IALPHA=00 ! Integration weight proportional to z^{\ell-1/2}
C               ! This is the physicist's choice for controlled reasons
C
      IF (LOGWRI.GT.0) THEN
C
          STRING(1:16)=PRONAM//'_t.log'
C
          OPEN(LOGFIL,FILE=STRING(1:16),STATUS='UNKNOWN',
     *                                  FORM='FORMATTED')
C
          WRITE(LOGFIL,'(/,''Entering the MAIN program '',A10)')
     *                                                   PRONAM
C
          WRITE(LOGFIL,'(/,''ATTENTION: VERSION IALPHA='',I1,'' IS ''
     *                     ''USED =>>'')')      IALPHA
          WRITE(LOGFIL,'(11X,''Weight = z^{\ell-1/2}'')')
C
      END IF
C
C=======================================================================
C=======================================================================
C     Opening the log-files for the output-control information
C=======================================================================
C
      IF (LOGWRI.GT.0) THEN
C
          WRITE(LOGFIL,'(/,''Finished the do-loop over nuclei '',
     *                        ''<<= experimental input'')')
C
          WRITE(LOGFIL,'(/,''In MAIN, LOG-TYPE files:'',/)')
C_______________________________________________________________________
C
          STRING(1:16)=PRONAM//'_p.log'
C
          OPEN(UNIT=LOGPRO,FILE=STRING(1:16),STATUS='UNKNOWN',
     *                                       FORM='FORMATTED')
C
          WRITE(LOGFIL,'(9X,''Opening the file: '',A,
     *                   1X,''UNIT='',I2)') STRING(1:16),LOGPRO
C_______________________________________________________________________
C
          STRING(1:16)=PRONAM//'_n.log'
C
          OPEN(UNIT=LOGNEU,FILE=STRING(1:16),STATUS='UNKNOWN',
     *                                       FORM='FORMATTED')
C
          WRITE(LOGFIL,'(9X,''Opening the file: '',A,
     *                   1X,''UNIT='',I2)') STRING(1:16),LOGNEU
C_______________________________________________________________________
C
          STRING(1:21)=PRONAM//'_iresul.log'
          OPEN(UNIT=IRESUL,FILE=STRING(1:21),STATUS='UNKNOWN',
     *                                       FORM='FORMATTED')
          WRITE(LOGFIL,'(9X,''Opening the file: '',A,
     *                   1X,''UNIT='',I2)') STRING(1:21),IRESUL
C_______________________________________________________________________
C
          IF (LOGWRI.GE.5) THEN
C
              STRING(1:21)=PRONAM//'_bisect.log'
              OPEN(UNIT=LOGBIS,FILE=STRING(1:21),STATUS='UNKNOWN',
     *                                           FORM='FORMATTED')
              WRITE(LOGFIL,'(/,''Opening BISECT  logfile from '',A,/)')
     *                                           PRONAM
          END IF
C
C=======================================================================
C
          OPEN(LOGCHN,FILE='Results/Log_file_CHIN.d',STATUS='UNKNOWN',
     *                                               FORM='FORMATTED')
          WRITE(LOGFIL,'(9X,''Opening the file:'',
     *                   1X,''Results/Log_file_CHIN.d'',
     *                   1X,''UNIT='',I2)') LOGHIN
C_______________________________________________________________________
C
          OPEN(LOGCHP,FILE='Results/Log_file_CHIP.d',STATUS='UNKNOWN',
     *                                               FORM='FORMATTED')
          WRITE(LOGFIL,'(  9X,''Opening the file:'',
     *                     1X,''Results/Log_file_CHIP.d'',
     *                     1X,''UNIT='',I2)') LOGHIP
C_______________________________________________________________________
C
          WRITE(LOGFIL,'()')
C
          WRITE(LOGFIL,'(/,''Entering NAMELI from MAIN'')')
C
          WRITE(N_OUTP,'(///,''ATTENTION: VERSION IALPHA='',I1,'' IS ''
     *                       ''USED !!!'',/)')    IALPHA
C
      END IF
C
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C  @@@ IRENE, check this: IS IT A RIGHT PLACE for "this" to be HERE?
      CHIMI1       =9999. ! used only in printing routine
      CHIMI2_PROTON=9999. ! used only in printing routine      
      CHIMI3_PROTON=9999. ! used only in printing routine      
      CHIMI2_NEUTRS=9999. ! used only in printing routine      
      CHIMI3_NEUTRS=9999. ! used only in printing routine      
C
C=======================================================================
C     Routine reading the input data for this particular run
C=======================================================================
C
C     NGAUSS     !  The number of nodes in the Gauss formula
C     NSHELL     !  Maximum value  of the possible L-indices
C
      CALL NAMELI(NUCACT,ITAKNU,INPUTZ,INPUTN,INPSYM,NSHELL_PROTON,
     *                                               NSHELL_NEUTRS,
     *                   LDRAND,NGAUSS,ENEMAX,I_SEED,ISIMPL,NEWSED,
     *                          AOSCIL_PROTON,AOSCIL_NEUTRS,HOMEGA)
C
C=======================================================================
C     Reading and preparing the experimental charge-densities 
C     for a l l the nuclei which are experimentally available
C=======================================================================
C     
      IF (LOGWRI.GT.4) THEN
C
          WRITE(LOGFIL,'()')
C
          DO I=1,3
             WRITE(LOGFIL,'(22(''=''))')
          END DO
C
          WRITE(LOGFIL,'(/,''Entering PREPAR_DENEXP from MAIN'')')
C
      END IF
C
      CALL PREPAR_DENEXP
C
C=======================================================================
C     Reading the experimental information  for EACH nucleus 
C     spherical energy levels, radii, Fermi and gap energies 
C=======================================================================
C
      DO INUCLI=1,LDNUCL
         DO ILEVEL=1,NDLEXP
C
            LABEXP_PROTON(INUCLI,ILEVEL)='XXXXXX'
            LABEXP_NEUTRS(INUCLI,ILEVEL)='XXXXXX'
C
            EXPEXP_PROTON(INUCLI,ILEVEL)=999.999
            EXPEXP_NEUTRS(INUCLI,ILEVEL)=999.999
C
            IDEGEX_PROTON(INUCLI,ILEVEL)=999
            IDEGEX_NEUTRS(INUCLI,ILEVEL)=999
C
         END DO
      END DO
C
      IF (LOGWRI.GT.4) THEN
          WRITE(LOGFIL,'(/,''Entering the do-loop over nuclei from '',
     *                     ''MAIN, <<= experimental input'')')
      END IF
C
      LDNUCL_AUXILI=8
C
      DO INUCLI=1,LDNUCL_AUXILI
C      
         IZ_FIX=NUMB_Z(INUCLI)
         IN_FIX=NUMB_N(INUCLI)
C 
         IDEFCN=1
C
C=======================================================================
C        Optionally reading the experimental radii: protons...
C=======================================================================
C
         IF (LOGWRI.GT.4) THEN
C
             WRITE(LOGFIL,'()')
C
             DO I=1,3
                WRITE(LOGFIL,'(9X,22(''=''))')
             END DO
C
             WRITE(LOGFIL,'(/,9X,''Entering RDPROT_RMSRAD from MAIN,'',
     *                        1X,''INUCLI='',I1,'' LDNUCL='',I1,
     *                        1X,''IZ_FIX='',I2,'' IN_FIX='',I3)')
     *                             INUCLI,         LDNUCL,
     *                             IZ_FIX,         IN_FIX
         END IF
                                                               IFEXPE=1
         CALL RDPROT_RMSRAD(NDRADZ,RADPRO,CHGRMS,RERROR,IZ_FIX,IN_FIX,
     *                                                  IFEXPE,IPRIEX)
C
C=======================================================================
C        ... and then neutrons
C=======================================================================
C
         IF (LOGWRI.GT.4 .AND. IDEFCN.EQ.1) THEN
             WRITE(LOGFIL,'(/,9X,''Entering PROVID_RMS_NR from MAIN,'',
     *                        1X,''INUCLI='',I1,'' LDNUCL='',I1,
     *                        1X,''IZ_FIX='',I2,'' IN_FIX='',I3)')
     *                             INUCLI,         LDNUCL,
     *                             IZ_FIX,         IN_FIX
         END IF
C
         CALL PROVID_RMS_NR(NDRADN,RADNEU,ERRNEU,IZ_FIX,IN_FIX,IDEFCN,
     *                                                  IFEXPE,IPRIEX)
         RMSEXP_PROTON(INUCLI)=RADPRO
         RMSERR_PROTON(INUCLI)=RERROR
C      
         RMSEXP_NEUTRS(INUCLI)=RADNEU 
         RMSERR_NEUTRS(INUCLI)=ERRNEU   
C
C=======================================================================
C        Optionally reading experimental  single particle
C        levels; option for doubly magic spherical nuclei
C=======================================================================
C
         IF (LOGWRI.GT.4) THEN
             WRITE(LOGFIL,'(/,9X,''Entering READEX_LEVELS from MAIN, '',
     *                           ''INUCLI='',I1,'' LDNUCL='',I1,
     *                        1X,''IZ_FIX='',I2,'' IN_FIX='',I3)')
     *                             INUCLI,         LDNUCL,
     *                             IZ_FIX,         IN_FIX
         END IF
C     
         CALL READEX_LEVELS(IFEXPE,IDEFCN)
C
C=======================================================================     
C        Reading the experimental single particle energies p/n
C=======================================================================     
C
         ISOSPI=1
C
         IF (LOGWRI.GT.4) THEN
             WRITE(LOGFIL,'(/,9X,''Entering PREPAR_EXPLEV from MAIN,'',
     *                        1X,''ISOSPI='',I1,
     *                        1X,''IZ_FIX='',I2,'' IN_FIX='',I3)') 
     *                             ISOSPI,
     *                             IZ_FIX,         IN_FIX
         END IF
C
         CALL PREPAR_EXPLEV(IZ_FIX,IN_FIX,FERMEX,GAPEXP,DENSUP,DENSDW,
     *                             IDEFCN,WHATEX,ISOSPI,IFEXPE,INUCLI)
          
         FERMEX_PROTON(INUCLI)=FERMEX
         GAPEXP_PROTON(INUCLI)=GAPEXP
         DENSUP_PROTON(INUCLI)=DENSUP
         DENSDW_PROTON(INUCLI)=DENSDW
C
         ISOSPI=0      
C
         IF (LOGWRI.GT.4) THEN
             WRITE(LOGFIL,'(/,9X,''Entering PREPAR_EXPLEV from MAIN,'',
     *                        1X,''ISOSPI='',I1, 
     *                        1X,''IZ_FIX='',I2,'' IN_FIX='',I3)') 
     *                             ISOSPI,
     *                             IZ_FIX,         IN_FIX
         END IF
C
         CALL PREPAR_EXPLEV(IZ_FIX,IN_FIX,FERMEX,GAPEXP,DENSUP,DENSDW,
     *                             IDEFCN,WHATEX,ISOSPI,IFEXPE,INUCLI)
     
         FERMEX_NEUTRS(INUCLI)=FERMEX
         GAPEXP_NEUTRS(INUCLI)=GAPEXP
         DENSUP_NEUTRS(INUCLI)=DENSUP
         DENSDW_NEUTRS(INUCLI)=DENSDW 
C
C=======================================================================
C          
      END DO  ! Over INUCLI
C
C=======================================================================
C=======================================================================
C     Defining  the Gauss integration parameters  as well as
C     generalised Laguerre basis polynomials and derivatives
C=======================================================================
C=======================================================================
C   
      IF (LOGWRI.GT.0) THEN
          WRITE(LOGFIL,'(/,9X,''Entering GENBAS from MAIN'',/)')
      END IF
C
      CALL GENBAS(NDIM_N,NDIM_L,NDGAUS,NGAUSS,ANODES,AWEIGH,
     *            NSHELL_NEUTRS,BAUXIL,CAUXIL,EPSLAG,QRNODE,
     *            QRWEIG,QNORNL,QDFLAG,QFLAGN,QPOLYN,QDPOLY,
     *                   IALPHA,QPOLYN_DENSIT,QDPOLY_DENSIT)
C
C=======================================================================
C     Changing Nodes, Weights and Polynomials from REAL*16 to REAL*8
C=======================================================================
C
      MAXIML=NSHELL_NEUTRS
      NORDER=NSHELL_NEUTRS/2
      LORDER=NSHELL_NEUTRS
C
C=======================================================================
C  
C     Testing the precision of nodes and weights / printout neutralised
C     
C     WRITE(90,'(80(''#''),/,''#'',78X,''#'')')
C     WRITE(90,'(''#'',9X,A,9X,''#'',/,''#'',78X,''#'')')
C    *    'Comparing Real*16 with Real*8 obtained by applying R=REAL(Q)'
C     WRITE(90,'(80(''#''),/)')
C
C     WRITE(90,'(4X,''L'',4X,''IG'',15X,''NODE-Q/R'',28X,''WEIGHT-Q/R''
C    *                                                              )')
C       
      DO LQNUMB=0,MAXIML
C
C        WRITE(90,'(18X,''123456789 123456789 123456'',
C    *              10X,''123456789 123456789 123456'')') 
C      
         DO IPOINT=1,NGAUSS
C            
            XRNODE(IPOINT,LQNUMB)=REAL(QRNODE(IPOINT,LQNUMB))
            XRWEIG(IPOINT,LQNUMB)=REAL(QRWEIG(IPOINT,LQNUMB))
C           
C           WRITE(90,'(2I6.3,2E36.26,/,12x,2e36.26)') 
C    *                                  LQNUMB,IPOINT,
C    *                                  QRNODE(IPOINT,LQNUMB),
C    *                                  QRWEIG(IPOINT,LQNUMB),
C    *                                  XRNODE(IPOINT,LQNUMB),
C    *                                  XRWEIG(IPOINT,LQNUMB)
C           
         END DO
      END DO
C
C=======================================================================
C     Conclusion: There remain up to 16 coinciding decimals from this
C                 test. To be compared with 32 digits precision, when
C                 testing the orthogonality in quadruple precision =>
C     Conclusion: Double precision retains 16 decimal places??????
C=======================================================================
C=======================================================================
C
C     Similar as above but: for precision of normalisation constants,
C     Laguerre functions and their derivatives / printout neutralised
C
C       
C     WRITE(91,'(80(''#''),/,''#'',78X,''#'')')
C     WRITE(91,'(''#'',20X,A,29X,''#'',/,''#'',78X,''#'')')
C    *               'Comparing Real*16 with Real*8'
C     WRITE(91,'(80(''#''),/)')
C      
C     WRITE(91,'(4X,''N'',5X,''L'',20X,''NORM-Q/R'',
C    *           32X,''Lag^N_L'',32X,''D Lag^N_L'')')
C       
      DO NQNUMB=0,NORDER
C
C        WRITE(91,'(20X,''123456789 123456789 12345678'',
C    *              12X,''123456789 123456789 12345678'', 
C    *              12X,''123456789 123456789 12345678'')') 
C
         DO LQNUMB=0,LORDER
C            
            XNORNL(NQNUMB,LQNUMB)=REAL(QNORNL(NQNUMB,LQNUMB))
C            
            XFLAGN(NQNUMB,LQNUMB)=REAL(QFLAGN(NQNUMB,LQNUMB))
            XDFLAG(NQNUMB,LQNUMB)=REAL(QDFLAG(NQNUMB,LQNUMB))
C            
C           WRITE(91,'(2I6.3,3E40.28,/,12x,3E40.28)')
C    *                                 NQNUMB,LQNUMB,
C    *                                 QNORNL(NQNUMB,LQNUMB),
C    *                                 QFLAGN(NQNUMB,LQNUMB),
C    *                                 QDFLAG(NQNUMB,LQNUMB),
C    *                                 XNORNL(NQNUMB,LQNUMB),
C    *                                 XFLAGN(NQNUMB,LQNUMB),
C    *                                 XDFLAG(NQNUMB,LQNUMB)
C            
         END DO
      END DO
C
C=======================================================================
C
C                 N I H I L   N O V I   S U B   S O L E
C
C     Conclusion: There remain up to 16 coinciding decimals from this
C                 test. To be compared with 32 digits precision, when
C                 testing the orthogonality in quadruple precision =>
C=======================================================================
C=======================================================================
C
C     Similar to the above, but for the normalised Laguerre functions
C                                           with prints commented out
C
C     This is because the Great Artist expected to discover differ...
C       
C     WRITE(92,'(80(''#''),/,''#'',78X,''#'')')
C     WRITE(92,'(''#'',20X,A,29X,''#'',/,''#'',78X,''#'')')
C    *               'Comparing Real*16 with Real*8'
C     WRITE(92,'(80(''#''),/)')
C       
C     WRITE(92,'(4X,''L'',5X,''N'',4X,''IG'',18X,''Norm-L:Q/R'',
C    *                                       30X,''Norm-DL:Q/R'')')
C       
      DO N_MAIN=NSHELL_NEUTRS,0,-1
C         
         DO LQNUMB=N_MAIN,0,-2
            NQNUMB=(N_MAIN-LQNUMB)/2
C
C           WRITE(92,'(26X,''123456789 123456789 12345678'',
C    *                 12X,''123456789 123456789 12345678'')') 
C            
            DO IPOINT=1,NGAUSS
C                
               XPOLYN(IPOINT,NQNUMB,LQNUMB)=
     *         REAL(QPOLYN(IPOINT,NQNUMB,LQNUMB))
C     
               XDPOLY(IPOINT,NQNUMB,LQNUMB)=
     *         REAL(QDPOLY(IPOINT,NQNUMB,LQNUMB))
C     
C              WRITE(92,'(3I6.3,2E40.28,/,18x,2E40.28)')
C    *                                    LQNUMB,NQNUMB,IPOINT,
C    *                                    QPOLYN(IPOINT,NQNUMB,LQNUMB),
C    *                                    QDPOLY(IPOINT,NQNUMB,LQNUMB),
C    *                                    XPOLYN(IPOINT,NQNUMB,LQNUMB),
C    *                                    XDPOLY(IPOINT,NQNUMB,LQNUMB)
C             
               DO NAUXIL=0,NORDER
                  DO LAUXIL=0,LORDER
C                   
                     XPOLYN_DENSIT(IPOINT,NAUXIL,LAUXIL,LQNUMB)=
     *               REAL(QPOLYN_DENSIT(IPOINT,NAUXIL,LAUXIL,LQNUMB))
C                   
                     XDPOLY_DENSIT(IPOINT,NAUXIL,LAUXIL,LQNUMB)=
     *               REAL(QDPOLY_DENSIT(IPOINT,NAUXIL,LAUXIL,LQNUMB))
C                   
                  END DO
               END DO
C            
            END DO ! IPOINT of Gauss integration
         END DO ! LQNUMB   
      END DO ! N_MAIN
C
C=======================================================================
C     Here we wish to calm down the passionate reader who at this point
C     is dying of curiosity about results of these profound efforts ...
C     And the conclusion is … surprise / surprise =>> the same as above
C=======================================================================
C=======================================================================
C
      DO INUCLI=1,LDNUCL
C          
         IZ_FIX=NUMB_Z(INUCLI)
         IN_FIX=NUMB_N(INUCLI)      
C
C=======================================================================
C=======================================================================
C        Preparing the kinetic energy term; this part is needed
C        once per run of the code;  we avoid repeating the same
C        operations millions of times ...
C=======================================================================
C=======================================================================
C
         IF (LOGWRI.GT.5) THEN
             WRITE(LOGFIL,'(/,''Entering KINMAT from MAIN,  protons, '',
     *                        ''INUCLI='',I2,'' LDNUCL='',I2)')
     *                          INUCLI,         LDNUCL
         END IF
C
         CALL KINMAT(XPOLYN,XDPOLY,POLYNX,XRNODE,XRWEIG,ENEKIN_PROTON,
     *               HOMEGA,NDBASE,NDIM_N,NDIM_L,NSHELL_PROTON,NDGAUS,
     *                                           NGAUSS,NDNUCL,INUCLI)
C
         IF (LOGWRI.GT.5) THEN
             WRITE(LOGFIL,'(/,''Entering KINMAT from MAIN, neutrons, '',
     *                        ''INUCLI='',I2,'' LDNUCL='',I2)')
     *                          INUCLI,         LDNUCL
         END IF
C
         CALL KINMAT(XPOLYN,XDPOLY,POLYNX,XRNODE,XRWEIG,ENEKIN_NEUTRS,
     *               HOMEGA,NDBASE,NDIM_N,NDIM_L,NSHELL_NEUTRS,NDGAUS,
     *                                           NGAUSS,NDNUCL,INUCLI)
C
C=======================================================================     
C=======================================================================     
C
C        This an intrinsic-test routine, checking the fact that the 
C        total HO hamiltonian (the sum of the kinetic and potential
C        energies) is diagonal in the  {|n l s j m_j>}  basis,  the
C        diagonal terms being the eigen-energies 
C
C                         E_nl=(2n+l+3/2)*homega
C
C        Since the call below is testing the HO properties using HO 
C        wave functions,  rather than the actual physical potential,
C        it is of no interest [we verified the machine precision of
C        the numerical algorithm calculating the Laguerre-type wave
C        functions already]. 
C
C        IF (LOGWRI.GT.0) THEN
C            WRITE(LOGFIL,'(/,''Entering CHEKHO from MAIN,  '',
C    *                        ''NSHELL_PROTON only, IALPHA='',I2,1X,
C    *                                            ''INUCLI='',I2)')
C    *                                              IALPHA,INUCLI
C        END IF
C
C        CALL CHEKHO(XPOLYN,XDPOLY,POLYNX,XRNODE,XRWEIG,HOMEGA,
C    *               NDBASE,NDIM_N,NDIM_L,NSHELL_PROTON,NDGAUS,
C    *                             NGAUSS,IALPHA,NDNUCL,INUCLI)
C
C        The conclusions of these equally profound efforts is that...
C        harmonic oscillator wave functions,  in quadrupole precision
C        when multiplied by r^2 and integrated preserve the precision
C        So do the first order derivatives...
C
C=======================================================================     
C        Reading the order of the theoretical single-particle levels 
C=======================================================================     
C
         IF (LOGWRI.GT.5) THEN
             WRITE(LOGFIL,'(/,''Entering ORDLAB from MAIN, '',
     *                        ''NSHELL_PROTON only, IALPHA='',I1)')
     *                                              IALPHA
         END IF
C
         IF (IZ_FIX.LE.82 .AND. IN_FIX.LE.126) THEN
             CALL ORDLAB(NDIM_N,NDIM_L,LABORD_PROTON,LABORD_NEUTRS,
     *                                 LABSYM_PROTON,LABSYM_NEUTRS,
     *                                               NDNUCL,INUCLI)
         END IF
C
      END DO ! INUCLI
C
C=======================================================================     
C=======================================================================     
C=======================================================================     
C     The following is a single run for the standard 
C     WS Hamiltonian with universal parameters (WSSTAN)  
C=======================================================================
C=======================================================================     
C=======================================================================     
C
      IF (ISIMPL.EQ.1) THEN
C
          WRITE(LSCREN,'(''Simple run with the spherical WS '',
     *                          ''Hamiltonian'')')
C
          IF (LOGWRI.GT.0) THEN
              WRITE(LOGFIL,'(/,''Entering WSUNIV, SINGLE RUN'')')
          END IF 
C
          CALL ONERUN_WSUNIV
          
          IF (LOGWRI.GT.0) THEN
              WRITE(LOGFIL,'(/,''Exiting WSUNIV - SINGLE RUN'')')
          END IF
C
C=======================================================================
C       
          STOP 'Single-run with the standard WS Hamiltonian successful'
C      
      END IF
C
C=======================================================================     
C=======================================================================     
C=======================================================================     
C     The following is a single run for a user-choice of 
C     parameters and options.
C=======================================================================
C=======================================================================     
C=======================================================================     
C
      IF (IFTEST.EQ.1) THEN 
C
          IF (LOGWRI.GT.0) THEN
              WRITE(LOGFIL,'(/,''Entering ONERUN, IFTEST='',I1,1X,
     *                            ''from MAIN'')')   IFTEST
          END IF
C
C=======================================================================     
C
          CALL ONERUN
C
C=======================================================================     
C      
       STOP 'Stop! Single-run version successful'
C
      END IF
C
C=======================================================================     
C=======================================================================     
C=======================================================================     
C     B e l o w   t h e   p a r a m e t e r - f i t t i n g 
C=======================================================================     
C=======================================================================     
C=======================================================================     
C
      IF (IFFITS.EQ.1) THEN
C
          IDEFCN=1
C
C=======================================================================
C
          IFUNCT_EVALUS=0
C
          ICOUNT=0
          DO I=1,NDPARS
             IF (IFTAKE(I).EQ.1) THEN
                  ICOUNT=ICOUNT+1
                  ARGMNT(ICOUNT)=VMISTR(I)
                  TITPAR(ICOUNT)=TITLES(I)(2:8)//'PROTON'
             END IF
          END DO
C
C=======================================================================
C
          IF (LOGWRI.NE.0) THEN
              WRITE(LOGFIL,'()')
              DO I=1,5
                 WRITE(LOGFIL,'(''Entering LMMINI '',
     *                          ''from MAIN, IFFITS=1'')')
              END DO
          END IF
C          
          CALL LMMINI(ARGMNT,I_SEED,LDRAND,NEWSED,NDLAST)
C
          IF (LOGWRI.NE.0) THEN
              WRITE(LOGFIL,'()')
              DO I=1,5
                 WRITE(LOGFIL,'(''END OF RUN - '',
     *                          ''from MAIN, IFFITS=1'')')
              END DO
          END IF
C           
C=======================================================================
C          
          CALL CPUTIM('WSPH15',0)
          CALL TIMPRI
C           
C=======================================================================
C      
          STOP 'Stop! Parameter-fitting version successful'
C           
C=======================================================================
C
      END IF    
C
C=======================================================================     
C=======================================================================     
C=======================================================================
C
C                         M E S H - O P T I O N    
C 
C     Below we run an option of minimising over N variables (N chosen 
C     by the user in the input file), whereas all other variables are
C     fixed by the user. The results are stored on a disc for plotting 
C     analysis purposes ( M A P   P L O T S ) .
C     
C     For making easier the code, the user can only choose TWO varia-
C     bles in the input file with I_MESH(I) = 1
C
C=======================================================================     
C=======================================================================     
C=======================================================================     
C   
      IF (IFMESH.EQ.1) THEN
          
          IF (LOGWRI.GT.0) THEN
              WRITE(LOGFIL,'(/,''Running the option of a'',
     *                         '' Minimisation on a Rectangular'',
     *                         '' Mesh of points'')')
C
              WRITE(LOGFIL,'(/,9X,''Entering MESHIN_MAPING '',
     *                            ''from MAIN'')')
          END IF
C
          IF (LOGWRI.GT.0) THEN
              WRITE(LSCREN,'(''Setting LOGWRI=0 to avoid a '',
     *                       ''HUGE output file'')')
              LOGWRI=0
          END IF
C
          CALL MESHIN_MAPING(LDRAND)
C
          CALL CPUTIM('WSPH15',0)
          CALL TIMPRI
C           
C=======================================================================
C      
          STOP 'From MAIN: Mesh-option version successful'
C           
C=======================================================================
C      
      END IF ! IFMESH = 1
C
C=======================================================================     
C=======================================================================     
C=======================================================================
C
C      SINGLE PARTICLE ENERGIES AS A FUNCTION OF ONE PARAMETER 
C
C=======================================================================     
C=======================================================================     
C=======================================================================
C
      IF (ISPE_P.EQ.1) THEN
C
          IF (LOGWRI.GT.0) THEN
              WRITE(LOGFIL,'(/,''Running the option of '',
     *                         '' SPE as a function of '',
     *                         '' one parameter'')')
C
              WRITE(LOGFIL,'(/,9X,''Entering SPENER_OFPARA '',
     *                            ''from MAIN'')')
          END IF
C
C          IF (LOGWRI.GT.0) THEN
C              WRITE(LSCREN,'(''Setting LOGWRI=0 to avoid a '',
C     *                       ''HUGE output file'')')
C              LOGWRI=0
C          END IF
C
          CALL SPENER_OFPARA(LDRAND)
C
          CALL CPUTIM('WSPH15',0)
          CALL TIMPRI
C           
C=======================================================================
C      
          STOP 'From MAIN: SPE of parameter version successful'
C           
C=======================================================================
C
      END IF ! ISPE_P=1
C
C=======================================================================     
C=======================================================================     
C=======================================================================
C
C      MONTE-CARLO RUN: S.P.E. LEVELS IN HISTOGRAM FORMS
C
C=======================================================================     
C=======================================================================     
C=======================================================================
C
      IF (IFMOCA.EQ.1) THEN
C
          IF (LOGWRI.GT.0) THEN
              WRITE(LOGFIL,'(/,''Running the option of Monte-Carlo'')')
          END IF
C
          IF (IFPARA.EQ.1) THEN
              IF (LOGWRI.GT.0) WRITE(LOGFIL,'(/,9X,''Entering '',
     *                             ''SPENER_MONTEC from MAIN'')')
              CALL SPENER_MONTEC
          END IF
C
          IF (IFPSEU.EQ.1) THEN
              IF (LOGWRI.GT.0) WRITE(LOGFIL,'(/,9X,''Entering '',
     *                             ''FITTIN_MONTEC from MAIN'')')
              CALL FITTIN_MONTEC(I_SEED,NEWSED)
          END IF
C
          STOP 'From MAIN: Monte-Carlo SPE Histogram version successful'
C
      END IF
C           
C=======================================================================
C
      STOP 'ALLRIGHT'      
      END 
C
C=======================================================================
C=======================================================================
C=======================================================================

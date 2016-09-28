C FILE NAME = wspher02_nameli_15.f ! Keep this symbol:    $ident@string$
C=======================================================================
C=======================================================================
C=======================================================================
C      
C=======================================================================
C=======================================================================
C                     INPUT SUBROUTINE NAMELIST STYLE
C=======================================================================
C=======================================================================
C      
      SUBROUTINE NAMELI(NUCACT,ITAKNU,INPUTZ,INPUTN,INPSYM,
     *                                NSHELL_PROTON,NSHELL_NEUTRS,
     *                  LDRAND,NGAUSS,ENEMAX,I_SEED,ISIMPL,NEWSED,
     *                         AOSCIL_PROTON,AOSCIL_NEUTRS,HOMEGA)
C      
      INCLUDE   'MATDIM/NDPARS.f'
      INCLUDE   'MATDIM/NDRAND.f'
      INCLUDE   'MATDIM/NDRAUS.f'
      INCLUDE   'MATDIM/NDNUCL.f'
      INCLUDE   'MATDIM/NDLEXP.f'
      INCLUDE   'MATDIM/NDLAST.f'
      INCLUDE   'MATDIM/NDCOLO.f'
      INCLUDE   'MATDIM/NDMONT.f'
      INCLUDE   'MATDIM/NDBINS.f'
C      
      PARAMETER
     *         (N_READ=5)
C
      CHARACTER 
     *          KEYWOR*10,SPACIN*10,CHIDEF*6,WHATEX*6,STRING*256
      CHARACTER 
     *          EXTENT_PROTON*40,EXTENT_NEUTRS*40
      CHARACTER 
     *          LABPRO_REMOVE*6,LABNEU_REMOVE*6,
     *          REMOVE_PROTON*3,REMOVE_NEUTRS*3
      CHARACTER
     *          INPSYM*6,TYPCHI*6,NUCSYM*6
      DIMENSION
     *          ITAKNU(1:NDNUCL),
     *          INPUTZ(1:NDNUCL),
     *          INPUTN(1:NDNUCL),
     *          INPSYM(1:NDNUCL)
      DIMENSION
     *          AOSCIL_PROTON(1:NDNUCL),
     *          AOSCIL_NEUTRS(1:NDNUCL),
     *          HOMEGA(1:NDNUCL)
C
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
     *       /EXTREM/ V0CMAX,V0CMIN,XL_MAX,XL_MIN,XEFMAX,XEFMIN,
     *                R0CMAX,R0CMIN,R0SMAX,R0SMIN,REFMAX,REFMIN,
     *                A0CMAX,A0CMIN,A0SMAX,A0SMIN,AEFMAX,AEFMIN,
     *                                            CouMAX,CouMIN,
     *                              XPPMAX,XPPMIN,XPNMAX,XPNMIN,
     *                              XNPMAX,XNPMIN,XNNMAX,XNNMIN,
     *                              YPPMAX,YPPMIN,YPNMAX,YPNMIN,
     *                              YNPMAX,YNPMIN,YNNMAX,YNNMIN,
     *                              CPPMAX,CPPMIN,CPNMAX,CPNMIN,
     *                              CNPMAX,CNPMIN,CNNMAX,CNNMIN,
     *
     *                                            UPFACT,DWFACT
      COMMON
     *       /KAPPAS/ XKCMIN,XKCMAX,XACMIN,XACMAX,XRCMIN,XRCMAX,
     *                XKSMIN,XKSMAX,XASMIN,XASMAX,XRSMIN,XRSMAX,
     *                XKEMIN,XKEMAX,XAEMIN,XAEMAX,XREMIN,XREMAX,
     *                                            XCoMIN,XCoMAX
      COMMON
     *       /EXTKAP/ V0CMIN_KAPPAR,V0CMAX_KAPPAR,XKCMIN_KAPPAR,
     *                XKCMAX_KAPPAR,A0CMIN_KAPPAR,A0CMAX_KAPPAR,
     *                XACMIN_KAPPAR,XACMAX_KAPPAR,R0CMIN_KAPPAR,
     *                R0CMAX_KAPPAR,XRCMIN_KAPPAR,XRCMAX_KAPPAR,
     *                XL_MIN_KAPPAR,XL_MAX_KAPPAR,XKSMIN_KAPPAR,
     *                XKSMAX_KAPPAR,A0SMIN_KAPPAR,A0SMAX_KAPPAR,
     *                XASMIN_KAPPAR,XASMAX_KAPPAR,R0SMIN_KAPPAR,
     *                R0SMAX_KAPPAR,XRSMIN_KAPPAR,XRSMAX_KAPPAR
      COMMON
     *       /HBAR_V/ HBAR_C    
C 
      COMMON
     *       /NUCINF/ LDNUCL,
     *                NUCSYM(1:NDNUCL),
     *                NUMB_Z(1:NDNUCL),
     *                NUMB_N(1:NDNUCL)
      COMMON
     *       /FILEXT/ EXTENT_PROTON,EXTENT_NEUTRS
      COMMON
     *       /MODIFR/ RMSEXP_PROTON(1:NDNUCL),
     *                RMSEXP_NEUTRS(1:NDNUCL),
     *                RMSTHE_PROTON,RMSTHE_NEUTRS,
     *                RMODIF_PROTON,RMODIF_NEUTRS,
     *                              PUSHIN,DELTAR
      COMMON
     *       /REMLAB/ REMOVE_PROTON,REMOVE_NEUTRS
      COMMON
     *       /OUTLAB/ LABPRO_REMOVE(1:NDRAUS),
     *                LABNEU_REMOVE(1:NDRAUS)
      COMMON
     *       /OUTIND/ INDPRO_LETOUT(1:NDRAUS),
     *                INDNEU_LETOUT(1:NDRAUS)
      COMMON
     *       /RANRAN/ IRMFST,I_RAND,IRANDO
      COMMON
     *       /IFPAIR/ IF_PAI,IFGPAI,IFDPAI
      COMMON
     *       /SEARCH/ IFTAKE(1:NDPARS)
      COMMON
     *       /VLIMIT/ VMISTR(1:NDPARS),
     *                VMIMIN(1:NDPARS),
     *                VMIMAX(1:NDPARS),
     *                VMISTP(1:NDPARS)
      COMMON
     *       /CNTROL/ PARAMT(1:NDPARS),
     *                DPARAM(1:NDPARS)     
      COMMON
     *       /POTPOT/ PARPOT(1:NDPARS)
      COMMON
     *       /MESHIN/ I_MESH(1:NDPARS),
     *                XMIN_I(1:NDPARS),
     *                XMAX_I(1:NDPARS),
     *                MAXPAR(1:NDPARS)
      COMMON
     *       /DELTAF/ DFACTO(1:NDNUCL,1:2),
     *       /BCSDEL/ DELT_2(1:NDNUCL,1:2)
      COMMON 
     *       /CHITAK/ IF_SPE,IF_RAD,IF_GAP,IF_FER,IF_DEN,IF_RHO,IF_INV,
     *                LDCHI2,TYPCHI(1:NDNUCL),ITKCHI(1:NDNUCL)
      COMMON
     *       /TESFIT/ IFTEST,IFFITS,IFMESH,ISPE_P,IFMOCA
      COMMON
     *       /MINITE/ IDEFCN,ITECHI
      COMMON
     *       /WHERBE/ IRUNMI,NEWOLD 
      COMMON
     *       /PRINAL/ IMTALK
     *       /INPSEL/ WHATEX
     *       /WHICHI/ CHIDEF
     *       /DENSTR/ IFDENS
     *       /TENSOR/ IFTENS
      COMMON
     *       /TENOPT/ ICENTT,ISORBT,ITENSR
      COMMON
     *       /DEPPAI/ IFDEEP,IFPRON
      COMMON
     *       /SVDSVC/ SVDCUT
      COMMON
     *       /IFLAMB/ IFPAR1,IFPAR2,IFPAR3,IFPAR4,
     *                IFPAR5,IFPAR6,IFPAR7,IFPAR8,
     *                IFPA09,IFPA10,IFPA11,IFPA12
      COMMON
     *       /XLENGT/ XMIN_T,XMAX_T,XMIN_E,XMAX_E
      COMMON
     *       /ATTRIB/ ISTYLE,I_TYPE,ICOLOR(1:NDCOLO),ITHICK(1:NDCOLO)
      COMMON
     *       /TOLERS/ TOLERF,TOLERX,TOLERG,FACTOR
      COMMON
     *       /STOPAR/ EPSLAS,LDLAST
      COMMON
     *       /LAMUNI/ UNITLA
      COMMON
     *       /IFKAPP/ IFK_VC,IFK_VS,IFK_RC,IFK_RS,IFK_AC,IFK_AS
      COMMON
     *       /PARCOR/ IFCORR,IFRCVC,IFRCAC,IFVCAC,IFRSVS,IFRSAS,IFVSAS
      COMMON
     *       /STEPBI/ BISTEP
      COMMON
     *       /FACCOU/ COUFAC
      COMMON
     *       /CTSHIF/ V0SHIF_PROTON(1:NDNUCL),
     *                V0SHIF_NEUTRS(1:NDNUCL)
      COMMON
     *       /RMSREF/ RMSVAL_REFPRO(1:NDNUCL),
     *                RMSVAL_REFNEU(1:NDNUCL),
     *                RMSGLO_REFPRO,
     *                RMSGLO_REFNEU
      COMMON
     *       /RMSIND/ RMSIND_PROTON(1:NDNUCL),
     *                RMSIND_NEUTRS(1:NDNUCL)
      COMMON
     *       /MONTEC/ IFPSEU,IFPARA,LDMONT,LDBINS
      COMMON
     *       /GAUSSI/ PARPOT_XMEANS(1:NDPARS),
     *                PARPOT_SIGMAS(1:NDPARS)
      COMMON
     *       /PRFILE/ ISCREN,LOGFIL,LOGCHP,LOGCHN,LOGHIP,LOGHIN,
     *                LOGPRO,LOGNEU,LOGWRI,LSCREN,IRESUL,ICONVE,
     *                NOUTPT,LOGAUX,LOGENE,LOGRAD,LOGGRD,LOGMSH,
     *                                                   LOGBIS                                                 
C_______________________________________________________________________
C
      DATA
     *     HBAR_C /197.327053/,! Units: MeV*fm;   ERROR: (+/-)0.000059
     *     XMASSP /938.27231 /,! Units: MeV/c**2; ERROR: (+/-)0.00028
     *     XMASSN /939.56563 / ! Units: MeV/c**2; ERROR: (+/-)0.00028
C @@@ DISCUSS THAT FINALLY...
      DATA
     *     NGSMAX /80/  ! Maximum order in the Gauss integrations
     *     NSHMAX /20/  ! Effective maximum of shells guaranteeing @@@
C                         that  the quadruple precision integrals
C                         will be possible to calculate
C     
C=======================================================================
C     In this subroutine we read input data in the NAMELIST style
C=======================================================================
C
C=======================================================================
C
   1  CONTINUE     
C     
C=======================================================================
C     
      READ (N_READ,'(A10)',END=2)  KEYWOR
C
C=======================================================================
C=======================================================================
C     Here the proton and neutron numbers of the nucleus in question...
C=======================================================================
C=======================================================================
C
      IF (KEYWOR.EQ.'IFHOWTORUN') THEN
C
          READ (N_READ,*) ISIMPL,IFTEST,IFFITS,IFMESH,ISPE_P,IFMOCA,
     *                           IF_DIR,IFDENS,IFTENS,ISCREN,LOGWRI
C
C         In the case of a single-nucleus run IZ_FIX and IN_FIX will
C         define the nucleus in question, see however the option for
C         fitting the parameters for several nuclei at the same time
C
C         ISIMPL=1 =>> A single WS run with the standard Hamiltonian
C         IFTEST=1 =>> Testing: The WS_RUN routine with the standard
C                                                         parameters
C         IFFITS=1 =>> Fitting the WS-Hamiltonian parameters various
C                                           variants described below
C         IF_DIR=1 =>> Option with the Dirac form of the Hamiltonian
C         IFDENS=1 =>> Option with the density-dependent  spin-orbit
C         IFTENS=1 =>> Mean-field  tensor-interaction term  included
C
C         ISCREN=1 =>> Partial output on the screen for online cntrl
C         LOGWRI>0 =>> Controlling the format of the log-file output,
C                      If 0, no logfile, If > 0, various variants of
C                      the details of the information.
C_______________________________________________________________________
C
          WRITE(LSCREN,'(/,9X,''From NAMELI, LOGWRI='',i1)') LOGWRI
C
          WRITE(LOGFIL,'(/,9X,''From NAMELI, LOGWRI='',i1)') LOGWRI
C
          WRITE(LOGFIL,'(/,A)') KEYWOR
C_______________________________________________________________________
C
          IF (ISIMPL.EQ.1) THEN
C
              IF (ISCREN.GT.0) THEN
                  WRITE(LSCREN,'(''Single-run with Universal'',1X, 
     *                           ''Parameters Variant Selected'')')
              END IF
C
              WRITE(LOGFIL,'( 9X,''Single-run with Universal'',1X, 
     *                           ''Parameters Variant Selected'')')
          END IF
C_______________________________________________________________________
C
          IF (IFTEST.EQ.1) THEN
C
              IF (ISCREN.GT.0) THEN
                  WRITE(LSCREN,'(''Single-run with User'',1X, 
     *                           ''Parameters Variant Selected'')')
              END IF
C
              WRITE(LOGFIL,'( 9X,''Single-run with User'',1X, 
     *                           ''Parameters Variant Selected'')')
          END IF
C_______________________________________________________________________
C
          IF (IFFITS.EQ.1) THEN
C
              IF (ISCREN.GT.0) THEN
                  WRITE(LSCREN,'(''Parameter Fitting Variant '',
     *                           ''Selected'')')
              END IF
C
              WRITE(LOGFIL,'( 9X,''Parameter Fitting Variant '',
     *                           ''Selected'')')
          END IF
C_______________________________________________________________________
C
          IF (IFMESH.EQ.1) THEN
C
              IF (ISCREN.GT.0) THEN
                  WRITE(LSCREN,'(''Parameter Mesh Variant '',
     *                           ''Selected'')')
              END IF
C
              WRITE(LOGFIL,'( 9X,''Parameter Mesh Variant '',
     *                           ''Selected'')')
          END IF
C_______________________________________________________________________      
C
          IF (IF_DIR.EQ.1) THEN
C
              IF (ISCREN.NE.0) THEN
                  WRITE(LSCREN,'(''Dirac Equation Variant '',
     *                           ''Selected'')')
              END IF
C
              WRITE(LOGFIL,'( 9X,''Dirac Equation Variant '',
     *                           ''Selected'')')
C
              STOP 'DIRAC option not implemented yet'
C
          END IF
C_______________________________________________________________________
C
          IF (IFDENS.EQ.0) THEN
C
              IF (ISCREN.GT.0) THEN
                  WRITE(LSCREN,'(''Calculation WITHOUT density-'',
     *                           ''dependent spin-orbit'')')
              END IF
C
              WRITE(LOGFIL,'( 9X,''Calculation WITHOUT density-'',
     *                           ''dependent spin-orbit'')')
          ELSE
C
              IF (ISCREN.GT.0) THEN
                  WRITE(LSCREN,'(''Calculating WITH density-'',
     *                           ''dependent spin-orbit'')')
              END IF
C
              WRITE(LOGFIL,'( 9X,''Calculating WITH density-'',
     *                           ''dependent spin-orbit'')')
          END IF
C_______________________________________________________________________
C
          IF (IFTENS.EQ.0) THEN
C
              IF (ISCREN.GT.0) THEN
                  WRITE(LSCREN,'(''Calculation WITHOUT tensor '',
     *                           ''contribution'')')
              END IF
C
              WRITE(LOGFIL,'( 9X,''Calculation WITHOUT tensor '',
     *                           ''contribution'')')
          ELSE
C
              IF (ISCREN.GT.0) THEN
                  WRITE(LSCREN,'(''Calculating WITH tensor '',
     *                           ''contribution'')')
              END IF
C
              WRITE(LOGFIL,'( 9X,''Calculating WITH tensor '',
     *                           ''contribution'')')
          END IF
C
      END IF
C
C=======================================================================
C=======================================================================
C=======================================================================
C     Reading what contributions to chi^2 which we take into account
C=======================================================================
C=======================================================================
C=======================================================================
C  
      IF (KEYWOR.EQ.'IFTAKECHI2') THEN
C
          READ(N_READ,*) IF_SPE,IF_RAD,IF_GAP,IF_FER,IF_DEN,IF_RHO,
     *                                                      IF_INV
C
          WRITE(LOGFIL,'(A)') KEYWOR
          WRITE(LOGFIL,'(9X,''Fitting requests:'')')
C
          IF (IF_SPE.EQ.1) THEN
              WRITE(LOGFIL,'(27X,''Single Particle levels - YES'')')
              ITKCHI(1)=1
          ELSE
              WRITE(LOGFIL,'(27X,''Single Particle levels - NO'')')
              ITKCHI(1)=0
          END IF
C 
          IF (IF_RAD.EQ.1) THEN
              WRITE(LOGFIL,'(27X,''R.M.S. radii ......... - YES'')')
              ITKCHI(2)=1
          ELSE
              WRITE(LOGFIL,'(27X,''R.M.S. radii ......... - NO'')')
              ITKCHI(2)=0
          END IF
C
          IF (IF_GAP.EQ.1) THEN
              WRITE(LOGFIL,'(27X,''Main gap sizes ....... - YES'')')
              ITKCHI(3)=1
          ELSE
              WRITE(LOGFIL,'(27X,''Main gap sizes ....... - NO'')')
              ITKCHI(3)=0
          END IF
C
          IF (IF_FER.EQ.1) THEN
              WRITE(LOGFIL,'(27X,''Average Fermi levels   - YES'')')
              ITKCHI(4)=1
          ELSE
              WRITE(LOGFIL,'(27X,''Average Fermi levels   - NO'')')
              ITKCHI(4)=0
          END IF
C
          IF (IF_DEN.EQ.1) THEN
              WRITE(LOGFIL,'(27X,''Densities =/- shells   - YES'')')
              ITKCHI(5)=1
          ELSE
              WRITE(LOGFIL,'(27X,''Densities =/- shells   - NO'')')
              ITKCHI(5)=0
          END IF
C
          IF (IF_RHO.EQ.1) THEN
              WRITE(LOGFIL,'(27X,''Densities rho(r) ..... - YES'')')
              ITKCHI(6)=1
          ELSE
              WRITE(LOGFIL,'(27X,''Densities rho(r) ..... - NO'')')
              ITKCHI(6)=0
          END IF
C
          IF (IF_INV.EQ.1) THEN
              WRITE(LOGFIL,'(27X,''Level inversions ..... - YES'')')
              ITKCHI(7)=1
          ELSE
              WRITE(LOGFIL,'(27X,''Level inversions ..... - NO'')')
              ITKCHI(7)=0
          END IF
C
          TYPCHI(1)='IF_SPE'
          TYPCHI(2)='IF_RAD'
          TYPCHI(3)='IF_GAP'
          TYPCHI(4)='IF_FER'
          TYPCHI(5)='IF_DEN'
          TYPCHI(6)='IF_RHO'
          TYPCHI(7)='IF_INV'
C
          LDCHI2=7
C
      END IF
C
C=======================================================================
C=======================================================================
C     Nuclei chosen for the minimisation
C=======================================================================
C=======================================================================
C
      IF (KEYWOR.EQ.'EXPER_FILE') THEN
C
          READ (N_READ,*) IFDEEP,IFPRON
C
          WRITE(LOGFIL,'(A)') KEYWOR
          WRITE(LOGFIL,'(9X,''Experimental file:'')')
C
          IF (IFDEEP.EQ.1) THEN
              WRITE(LOGFIL,'(27X,''Deeply-bound states... - YES'')')
          ELSE
              WRITE(LOGFIL,'(27X,''Deeply-bound states... - NO'')')
          END IF
C
          IF (IFPRON.EQ.1) THEN
              WRITE(LOGFIL,'(27X,''P-N pairing N=Z nuclei - YES'')')
          ELSE
              WRITE(LOGFIL,'(27X,''P-N pairing N=Z nuclei - NO'')')
          END IF
C
      END IF
C
C=======================================================================
C=======================================================================
C     How many nuclei to chose from for the minimisation
C=======================================================================
C=======================================================================
C
      IF (KEYWOR.EQ.'HOWMANYNUC') THEN
C
          READ(N_READ,*) LDNUCL
C
          IF (LDNUCL.GT.NDNUCL) THEN
              WRITE(0,'(/,''Alarm in NAMELI: LDNUCL= '',I3,
     *                    ''is .GT. NDNUCL= '',I3,/)') LDNUCL,NDNUCL
              STOP 'STOP in NAMELI: LDNUCL.GT.NDNUCL'
          END IF
C 
          WRITE(LOGFIL,'(A)') KEYWOR
          WRITE(LOGFIL,'(9X,''Scanning data base for LDNUCL='',I2,
     *                      '' nuclei'')') LDNUCL
C
      END IF
C
C=======================================================================
C=======================================================================
C     Which nuclei really chosen?
C=======================================================================
C=======================================================================
C
      IF (KEYWOR.EQ.'ITAKENUCLS') THEN
C        
          READ(N_READ,*) (ITAKNU(I),I=1,LDNUCL)
C        
          NUCACT=0
C        
          DO I=1,LDNUCL
C          
             IF (ITAKNU(I).EQ.1) THEN
                 NUCACT=NUCACT+1
             END IF
C        
          END DO
C_______________________________________________________________________
C
          WRITE(LOGFIL,'(A)') KEYWOR
C_______________________________________________________________________
C        
          IF (NUCACT.EQ.0) THEN
              WRITE(LSCREN,'(/,''There is NO active nucleus '',
     *                         ''for the calculation!'',/)')
C
              STOP 'STOP in NAMELI: No active nucleus in the input file'
C
          END IF
C
          IF (NUCACT.GT.LDNUCL) THEN
              WRITE(0,'(/,''Alarm in NAMELI: NUCACT= '',I3,''is .GT. ''
     *                    ''LDNUCL= '',i3,/)') NUCACT,LDNUCL
              STOP 'STOP: NUCACT.GT.LDNUCL'
          END IF
C_______________________________________________________________________
C        
          LDACTV=0
C_______________________________________________________________________
C        
          IF (ITAKNU(1).EQ.1) THEN
C            
              LDACTV=LDACTV+1
C            
              INPUTZ(LDACTV)=8
              INPUTN(LDACTV)=8
              INPSYM(LDACTV)='   16O'
C
              WRITE(LOGFIL,'(9X,''Actually taken: '',A6)')
     *                                      INPSYM(LDACTV)
          END IF
C_______________________________________________________________________
C        
          IF (ITAKNU(2).EQ.1) THEN
C            
              LDACTV=LDACTV+1
C            
              INPUTZ(LDACTV)=20
              INPUTN(LDACTV)=20
              INPSYM(LDACTV)='  40Ca'
C
              WRITE(LOGFIL,'(9X,''Actually taken: '',A6)')
     *                                      INPSYM(LDACTV)
          END IF
C_______________________________________________________________________
C        
          IF (ITAKNU(3).EQ.1) THEN
C            
              LDACTV=LDACTV+1
C            
              INPUTZ(LDACTV)=20
              INPUTN(LDACTV)=28
              INPSYM(LDACTV)='  48Ca'
C
              WRITE(LOGFIL,'(9X,''Actually taken: '',A6)')
     *                                      INPSYM(LDACTV)
          END IF
C_______________________________________________________________________
C        
          IF (ITAKNU(4).EQ.1) THEN
C            
              LDACTV=LDACTV+1
C            
              INPUTZ(LDACTV)=28
              INPUTN(LDACTV)=28
              INPSYM(LDACTV)='  56Ni'
C
              WRITE(LOGFIL,'(9X,''Actually taken: '',A6)')
     *                                      INPSYM(LDACTV)
          END IF
C_______________________________________________________________________
C        
          IF (ITAKNU(5).EQ.1) THEN
C            
              LDACTV=LDACTV+1
C            
              INPUTZ(LDACTV)=40
              INPUTN(LDACTV)=50
              INPSYM(LDACTV)='  90Zr'
C
              WRITE(LOGFIL,'(9X,''Actually taken: '',A6)')
     *                                      INPSYM(LDACTV)
          END IF
C_______________________________________________________________________
C        
          IF (ITAKNU(6).EQ.1) THEN
C            
              LDACTV=LDACTV+1
C            
              INPUTZ(LDACTV)=50
              INPUTN(LDACTV)=82
              INPSYM(LDACTV)=' 132Sn'
C
              WRITE(LOGFIL,'(9X,''Actually taken: '',A6)')
     *                                      INPSYM(LDACTV)
          END IF
C_______________________________________________________________________
C        
          IF (ITAKNU(7).EQ.1) THEN
C            
              LDACTV=LDACTV+1
C            
              INPUTZ(LDACTV)=64
              INPUTN(LDACTV)=82
              INPSYM(LDACTV)=' 146Gd'
C
              WRITE(LOGFIL,'(9X,''Actually taken: '',A6)')
     *                                      INPSYM(LDACTV)
          END IF
C_______________________________________________________________________
C        
          IF (ITAKNU(8).EQ.1) THEN
C            
              LDACTV=LDACTV+1
C            
              INPUTZ(LDACTV)=82
              INPUTN(LDACTV)=126
              INPSYM(LDACTV)=' 208Pb'
C
              WRITE(LOGFIL,'(9X,''Actually taken: '',A6)')
     *                                      INPSYM(LDACTV)
          END IF
C_______________________________________________________________________
C
          IF (ITAKNU(9).EQ.1) THEN
C
              LDACTV=LDACTV+1
C
              INPUTZ(LDACTV)=114
              INPUTN(LDACTV)=164
              INPSYM(LDACTV)=' 278Fl'
C
              WRITE(LOGFIL,'(9X,''Actually taken: '',A6)')
     *                                      INPSYM(LDACTV)
C
          END IF
C_______________________________________________________________________
C
          IF (ITAKNU(10).EQ.1) THEN
C
              LDACTV=LDACTV+1
C
              INPUTZ(LDACTV)=114
              INPUTN(LDACTV)=186
              INPSYM(LDACTV)=' 300Fl'
C
              WRITE(LOGFIL,'(9X,''Actually taken: '',A6)')
     *                                      INPSYM(LDACTV)
C
          END IF
C_______________________________________________________________________
C
          IF (ITAKNU(11).EQ.1) THEN
C
              LDACTV=LDACTV+1
C
              INPUTZ(LDACTV)=114
              INPUTN(LDACTV)=220
              INPSYM(LDACTV)=' 334Fl'
C
              WRITE(LOGFIL,'(9X,''Actually taken: '',A6)')
     *                                      INPSYM(LDACTV)
C
          END IF
C_______________________________________________________________________
C
          IF (ITAKNU(12).EQ.1) THEN
C
              LDACTV=LDACTV+1
C
              INPUTZ(LDACTV)=126
              INPUTN(LDACTV)=164
              INPSYM(LDACTV)='290Ubh'
C
              WRITE(LOGFIL,'(9X,''Actually taken: '',A6)')
     *                                      INPSYM(LDACTV)
C
          END IF
C_______________________________________________________________________
C
          IF (ITAKNU(13).EQ.1) THEN
C
              LDACTV=LDACTV+1
C
              INPUTZ(LDACTV)=126
              INPUTN(LDACTV)=186
              INPSYM(LDACTV)='312Ubh'
C
              WRITE(LOGFIL,'(9X,''Actually taken: '',A6)')
     *                                      INPSYM(LDACTV)
C
          END IF
C_______________________________________________________________________
C
          IF (ITAKNU(14).EQ.1) THEN
C
              LDACTV=LDACTV+1
C
              INPUTZ(LDACTV)=126
              INPUTN(LDACTV)=220
              INPSYM(LDACTV)='346Ubh'
C
              WRITE(LOGFIL,'(9X,''Actually taken: '',A6)')
     *                                      INPSYM(LDACTV)
C
          END IF
C_______________________________________________________________________
C
          IF (LDACTV.NE.NUCACT) THEN
C
              WRITE(LSCREN,'(''WARNING!: LDACTV='',I3,
     *                       '' is NOT equal to NUCACT='',/)')
     *                                   LDACTV,NUCACT
C
              STOP 'Problem in KEYWOR=ITAKENUCLS in NAMELI'
C
          END IF
C_______________________________________________________________________
C      
          IF (NUCACT.EQ.1) THEN
C
              WRITE(LSCREN,'(''Calculations with '',I3,
     *                       '' nucleus: '',(A6))')
     *                          NUCACT,INPSYM(1)
          ELSE
C
              WRITE(LSCREN,'(/,9X,''Calculations with '',I3,
     *                         '' nuclei: '',<NUCACT>(A6))')
     *                         NUCACT,(INPSYM(I),I=1,NUCACT)
          END IF
C
          IF (NUCACT.EQ.1) THEN
      
              WRITE(LOGFIL,'(9X,''Calculations with'',I3,1X,
     *                          ''nucleus: '',(A6),'' --> Z='',
     *                                        I3,''; N='',I3)')
     *                     NUCACT,INPSYM(1),INPUTZ(1),INPUTN(1)
          ELSE
C
              WRITE(LOGFIL,'(9X,''Calculations with'',I3,1X,
     *                          ''nuclei:  '',(A6),'' --> Z='',
     *                                        I3,''; N='',I3)')
     *                     NUCACT,INPSYM(1),INPUTZ(1),INPUTN(1)
              DO I=2,NUCACT
C
                 WRITE(LOGFIL,'(39X,A6,'' --> Z='',I3,''; N='',I3)')
     *                                 INPSYM(I),INPUTZ(I),INPUTN(I)
              END DO
C
          END IF
C
      END IF
C
C=======================================================================
C
      IF (KEYWOR.EQ.'NUCLEIDATA') THEN
C          
          READ(N_READ,'(12X,<LDNUCL>(A6,1X))') (NUCSYM(I),I=1,LDNUCL)
C                                        Symbols of nuclei; 16O, etc.
          READ(N_READ,*) (NUMB_Z(I),I=1,LDNUCL)
          READ(N_READ,*) (NUMB_N(I),I=1,LDNUCL)
C_______________________________________________________________________
C
          WRITE(LOGFIL,'(A)') KEYWOR
C
          WRITE(LOGFIL,'(9X,''READING:'',<LDNUCL>(1X,A6))')
     *                  (NUCSYM(I),I=1,LDNUCL)
          WRITE(LOGFIL,'(17X,<LDNUCL>(1X,I6))')(NUMB_Z(I),I=1,LDNUCL)
          WRITE(LOGFIL,'(17X,<LDNUCL>(1X,I6))')(NUMB_N(I),I=1,LDNUCL)
C_______________________________________________________________________
C
          IF (LOGWRI.GE.5) THEN
C          
              WRITE(LOGFIL,'(A)') KEYWOR
              WRITE(LOGFIL,'(/,9X,''Our data base contains '',
     *                            ''the following'',/,9X,
     *                            ''doubly-magic spherical nuclei:'',
     *                                                            /)')
              DO I=1,LDNUCL
                 WRITE(LOGFIL,'(12X,A6,3X,I4,I4)') NUCSYM(I),NUMB_Z(I),
     *                                                       NUMB_N(I)
              END DO
C
          END IF
C          
      END IF
C
C=======================================================================
C=======================================================================
C     Pairing options and yes/no for the including of the density
C=======================================================================
C=======================================================================
C
      IF (KEYWOR.EQ.'IF-PAIRING') THEN
C
          READ (N_READ,*) IFGPAI,IFDPAI
C

          IF (IFGPAI.NE.0 .AND. IFDPAI.NE.0) THEN
C
              WRITE(LSCREN,'(''Conflicting option IFGPAI='',i1,1x,
     *                  ''and IFDPAI='',i1,'' forcing Delta-option'')')
     *                    IFGPAI,IFDPAI
              IFGPAI=0 
          END IF
C
          IF (IFGPAI.NE.0 .OR. IFDPAI.NE.0) THEN
              IF_PAI=1
          END IF          
C_______________________________________________________________________
C
          WRITE(LOGFIL,'(A)') KEYWOR
C_______________________________________________________________________
C
          IF (IFGPAI.NE.0) THEN
C
              WRITE(LSCREN,
     *            '(''Pairing from empirical G-parametrisation'')')
C
              WRITE(LOGFIL,
     *            '(9X,''Pairing from empirical G-parametrisation'')')
          END IF          
C_______________________________________________________________________
C
          IF (IFDPAI.NE.0) THEN
C
              WRITE(LSCREN,
     *       '(''Pairing from empirical-delta [mass-differences]'')')
C
              WRITE(LOGFIL,
     *       '(9X,''Pairing from empirical-delta [mass-differences]''
     *                                                            )')
          END IF


C  @@@ IRENE: WHAT DOES THIS "IF" HAVE TO LOOK FOR UNDER PAIRING? MISLEADING
C  @@@ IRENE: WHAT DOES THIS "IF" HAVE TO LOOK FOR UNDER PAIRING? MISLEADING
C  @@@ IRENE: WHAT DOES THIS "IF" HAVE TO LOOK FOR UNDER PAIRING? MISLEADING
C  @@@ IRENE: WHAT DOES THIS "IF" HAVE TO LOOK FOR UNDER PAIRING? MISLEADING
C  @@@ IRENE: WHAT DOES THIS "IF" HAVE TO LOOK FOR UNDER PAIRING? MISLEADING
C  @@@ IRENE: WHAT DOES THIS "IF" HAVE TO LOOK FOR UNDER PAIRING? MISLEADING
C  @@@ IRENE: WHAT DOES THIS "IF" HAVE TO LOOK FOR UNDER PAIRING? MISLEADING
C  @@@ IRENE: WHAT DOES THIS "IF" HAVE TO LOOK FOR UNDER PAIRING? MISLEADING
C  @@@ IRENE: WHAT DOES THIS "IF" HAVE TO LOOK FOR UNDER PAIRING? MISLEADING
C  @@@ IRENE: WHAT DOES THIS "IF" HAVE TO LOOK FOR UNDER PAIRING? MISLEADING


C_______________________________________________________________________
C
          IF (IFDENS.EQ.0) THEN
C
              WRITE(LSCREN,'(''Calculating WITHOUT density-dependent '',
     *                       ''Spin-Orbit'')')
C
              WRITE(LOGFIL,'(9X,''Calculating WITHOUT '',
     *                          ''density-dependent Spin-Orbit'')')
          ELSE
C
              WRITE(LSCREN,'(''Calculating WITH density-dependent '',
     *                       ''Spin-Orbit'')')
C
              WRITE(LOGFIL,'(9X,''Calculating WITH '',
     *                          ''density dependent Spin-Orbit'')')
C
          END IF
C_______________________________________________________________________
C
      END IF
C
C=======================================================================
C=======================================================================
C     Pairing BISECTION options
C=======================================================================
C=======================================================================
C 
      IF (KEYWOR.EQ.'BISECTSTEP') THEN
C
          READ(N_READ,*) BISTEP
C
          IF (IFDPAI.EQ.1) THEN
              WRITE(LOGFIL,'(A)') KEYWOR
              WRITE(LOGFIL,'(9X,''1 BCS equation, BISTEP= '',F6.3)') 
     *                                            BISTEP
          ELSE
              WRITE(LOGFIL,'(A)') KEYWOR
              WRITE(LOGFIL,'(9X,''IFDPAI=0: '',
     *                          ''NO Delta-parameterised pairing '',
     *                          ''was activated'')')
          END IF
C
      END IF
C
C=======================================================================
C          
      IF (KEYWOR.EQ.'DELTAFACTO') THEN
C
          READ(N_READ,*) (DFACTO(I,1),I=1,LDNUCL) ! Protons
          READ(N_READ,*) (DFACTO(I,2),I=1,LDNUCL) ! Neutrons              
C
          IF (IFDPAI.EQ.1) THEN
C
              WRITE(LOGFIL,'(A)') KEYWOR
              WRITE(LOGFIL,'(9X,''INUCLI | PROTON:  DFACTO  '',
     *                          ''DELEXP  PRODUC'',
     *                          '' | NEUTRS:  DFACTO  DELEXP'',
     *                                            ''  PRODUC'')')
              DO INUCLI=1,LDNUCL
C                
                 IZ_FIX=NUMB_Z(INUCLI)
                 IN_FIX=NUMB_N(INUCLI)
C                 
                 DELTAP=DELEXP(IZ_FIX,IN_FIX,1)
                 DELTAN=DELEXP(IZ_FIX,IN_FIX,0)
C
                 DELT_2(INUCLI,1)=(DFACTO(INUCLI,1)*DELTAP)**2
                 DELT_2(INUCLI,2)=(DFACTO(INUCLI,2)*DELTAN)**2
C
                 WRITE(LOGFIL,'(9X,I4,14X,F6.2,2X,F6.2,2X,F6.2,11X,
     *                                       F6.2,2X,F6.2,2X,F6.2)') 
     *                              INUCLI,DFACTO(INUCLI,1),DELTAP,
     *                                     DFACTO(INUCLI,1)*DELTAP,
     *                                     DFACTO(INUCLI,2),DELTAN,
     *                                     DFACTO(INUCLI,2)*DELTAN
              END DO
C
          ELSE
C
              WRITE(LOGFIL,'(A)') KEYWOR
              WRITE(LOGFIL,'(9X,''IFDPAI=0: NO Delta-parameterised '',
     *                          ''pairing was activated'')')
          END IF
C
      END IF ! DELTAFACTO
C
C=======================================================================
C=======================================================================
C     Introducing (or not) the tensor interaction to the central, 
C     spin-orbit and tensor interaction as such to the hamiltonian
C=======================================================================
C=======================================================================
C
      IF (KEYWOR.EQ.'TENSOPTION') THEN
C
          READ (N_READ,*) ICENTT,ISORBT,ITENSR
CID:
C         ICENTT - add tensor to the central potential
C         ISORBT - add tensor to the spin-orbit potential
C         ITENSR - add tensor as such (to be clarified...)
CID. @@@ IRENE, LET US CLARIFY @@@
          WRITE(LOGFIL,'(A)') KEYWOR
C_______________________________________________________________________
C
          IF (ICENTT.EQ.0) THEN
C
              WRITE(LSCREN,'(''Calculation WITHOUT tensor '',
     *                       ''contribution in the central '',
     *                                        ''potential'')')
C
              WRITE(LOGFIL,'(9X,''Calculating WITHOUT tensor '',
     *                          ''contribution in the central '',
     *                                           ''potential'')')
          ELSE
C
              WRITE(LSCREN,'(''Calculating WITH tensor contribution '',
     *                       ''in the central potential'')')
C
              WRITE(LOGFIL,'(9X,''Calculating WITH tensor '',
     *                     ''contribution in the central potential'')')
C
          END IF 
C_______________________________________________________________________
C          
          IF (ISORBT.EQ.0) THEN
C
              WRITE(LSCREN,'(''Calculating WITHOUT tensor '',
     *                       ''contribution '',
     *                       ''in the spin-orbit potential'')')
C
              WRITE(LOGFIL,'(9X,''Calculating WITHOUT tensor '',
     *                 ''contribution in the spin-orbit potential'')')
C
          ELSE
C
              WRITE(LSCREN,'(''Calculating WITH tensor contribution '',
     *                       ''in the spin-orbit potential'')')
C
              WRITE(LOGFIL,'(9X,''Calculating WITH tensor '',
     *                  ''contribution in the spin-orbit potential'')')
C
          END IF
C_______________________________________________________________________
C 
          IF (ITENSR.EQ.0) THEN
C
              WRITE(LSCREN,'(''Calculating WITHOUT tensor term'')')
C
              WRITE(LOGFIL,'(9X,''Calculating WITHOUT tensor '',
     *                                                      ''term'')')
C
          ELSE
C
              WRITE(LSCREN,'(''Calculating WITH tensor term'')')
C
              WRITE(LOGFIL,'(9X,''Calculating WITH tensor term'')')
C
          END IF
C
      END IF
C
C=======================================================================
C=======================================================================
C     Parameters of the basis and the Gauss integration control
C=======================================================================
C=======================================================================
C
      IF (KEYWOR.EQ.'PARAMSBASE') THEN
C
          READ (N_READ,*) NSHELL_PROTON,NSHELL_NEUTRS,NGAUSS,HOMEG0
C
          WRITE(LOGFIL,'(A)') KEYWOR
          WRITE(LOGFIL,'(9X,''NSHELL_PROTON,NSHELL_NEUTRS,NGAUSS,'',
     *                                                   ''HOMEG0'')')
          WRITE(LOGFIL,'(9X,I7,7X,I7,8X,I3,3X,F6.3)')
     *                       NSHELL_PROTON,NSHELL_NEUTRS,NGAUSS,HOMEG0
C
          WRITE(LOGFIL,'(9X,''The following quantities are '',
     *                             ''calculated here ONCE FOREVER '',
     *                              ''including the reduced masses'')')
C
          WRITE(LOGFIL,'(9X,''INUCLI  IZ_FIX  IN_FIX  A_MASS  '',
     *                      ''X_MASS_PROTON  X_MASS_NEUTRS  '',
     *                      ''HOMEGA_VALUES  AOSCIL_PROTON  '',
     *                                     ''AOSCIL_NEUTRS'')')
C_______________________________________________________________________
C          
          DO I=1,LDNUCL
C           
             IZ_FIX=NUMB_Z(I)
             IN_FIX=NUMB_N(I)
C           
             A_MASS=IZ_FIX+IN_FIX
C          
             X_MASS_PROTON=XMASSP*((IZ_FIX-1)*XMASSP+IN_FIX*XMASSN)
     *                               /(IZ_FIX*XMASSP+IN_FIX*XMASSN)
C
             X_MASS_NEUTRS=XMASSN*(IZ_FIX*XMASSP+(IN_FIX-1)*XMASSN)
     *                               /(IZ_FIX*XMASSP+IN_FIX*XMASSN)
C
             HOMEGA(I)=HOMEG0/A_MASS**(1./3.)
             HOMEGA(I)=1.2000*HOMEGA(I)         
C
             AOSCIL_PROTON(I)=HBAR_C**2/HOMEGA(I)/X_MASS_PROTON
             AOSCIL_PROTON(I)=SQRT(AOSCIL_PROTON(I))
C
             AOSCIL_NEUTRS(I)=HBAR_C**2/HOMEGA(I)/X_MASS_NEUTRS
             AOSCIL_NEUTRS(I)=SQRT(AOSCIL_NEUTRS(I))
C
             XPCENT=100.0*X_MASS_PROTON/XMASSP
             XNCENT=100.0*X_MASS_NEUTRS/XMASSN  
C         
             WRITE(LOGFIL,'(9X,3(I6,2X),F6.0,2X,F7.3,F5.1,''%'',2X,
     *                                          F7.3,F5.1,''%'',
     *                                               3(F13.4,2X))')
     *                                     I,IZ_FIX,IN_FIX,A_MASS,
     *                                       X_MASS_PROTON,XPCENT,
     *                                       X_MASS_NEUTRS,XNCENT,
     *                                 HOMEGA(I),AOSCIL_PROTON(I),
     *                                           AOSCIL_NEUTRS(I)
C  
          END DO
C_______________________________________________________________________
C
          IF (NGAUSS.GT.NGSMAX) THEN
C
              WRITE(LSCREN,'(/,''NGAUSS='',I3,'' in Gauss-Laguerre '',
     *                         ''exceeds NGMAXX='',I2,/)') NGAUSS,NGSMAX
              STOP 'Stop! NGAUSS > NGSMAX in NAMELI'
C
          END IF
C_______________________________________________________________________
C
          IF (NSHELL_NEUTRS.GT.2*NGAUSS) THEN
C
              WRITE(LSCREN,'(/,''NSHELL_NEUTRS='',I2,'' exceeds '',
     *                         ''2*NGAUSS='',I3,/)')
     *                           NSHELL_NEUTRS,2*NGAUSS
              STOP 'Stop! NSHELL_NEUTRS > 2*NGAUSS in NAMELI'
C
          END IF
C
          IF (NSHELL_PROTON.GT.2*NGAUSS) THEN
C
              WRITE(LSCREN,'(/,''NSHELL_PROTON='',I2,'' exceeds '',
     *                         ''2*NGAUSS='',I3,/)')
     *                           NSHELL_PROTON,2*NGAUSS
              STOP 'Stop! NSHELL_PROTON > 2*NGAUSS in NAMELI'
C
          END IF
C_______________________________________________________________________
C
          IF (NSHELL_NEUTRS.GT.NSHMAX) THEN
C
              WRITE(LOGFIL,'(''NSHELL_NEUTRS='',I2,'' exceeds '',
     *                       ''NSHMAX='',I2)') NSHELL_NEUTRS,NSHMAX
          END IF
C
          IF (NSHELL_PROTON.GT.NSHMAX) THEN
C
              WRITE(LOGFIL,'(''NSHELL_PROTON='',I2,'' exceeds '',
     *                       ''NSHMAX='',I2)') NSHELL_PROTON,NSHMAX
          END IF
C_______________________________________________________________________
C
      END IF
C
C=======================================================================
C=======================================================================
C     One of the options for the experimental input
C=======================================================================
C=======================================================================
C
      IF (KEYWOR.EQ.'EXPOPTIONS') THEN
C
          READ (N_READ,*) WHATEX
C
          WRITE(LOGFIL,'(A)') KEYWOR
          WRITE(LOGFIL,'(  9X,''Options for the experimental '',
     *            ''input:'',/,9X,''Among: EXTRUE, EXOROS, SHELLM, '',
     *                         ''EXPNEW you selected: '',a6)') WHATEX
C
          WRITE(LSCREN,'(/,''Options for the experimental input:'',/,
     *                     ''Among: EXTRUE, EXOROS, SHELLM, EXPNEW '',
     *                     ''you selected: '',a6)') WHATEX
C
      END IF
C
C=======================================================================
C=======================================================================
C     Options related to the print control
C=======================================================================
C=======================================================================
C
      IF (KEYWOR.EQ.'WHATPRINTS') THEN
C
          READ (N_READ,*) ENEMAX
C
          WRITE(LOGFIL,'(A)') KEYWOR
          WRITE(LOGFIL,'(9X,''ENEMAX'',/,9X,F6.2)') ENEMAX
C
      END IF
C
C=======================================================================
C=======================================================================
C     Options related to the SVD routine (singular values cutoff)
C=======================================================================
C=======================================================================
C
C @@@ VERIFY THE FUNCTIONING
      IF (KEYWOR.EQ.'SVDWHATCUT') THEN
C
          READ (N_READ,*) SVDCUT
C
          WRITE(LOGFIL,'(A)') KEYWOR
          WRITE(LOGFIL,'(9X,''Singular Value Decomposition cut '',
     *                      ''SVDCUT='',e12.5)')
     *                        SVDCUT
C
          WRITE(LSCREN,'(/,''Singular Value Decomposition cut '',
     *                     ''SVDCUT='',e12.5)')
     *                       SVDCUT
      END IF
C
C=======================================================================
C=======================================================================
C     Which conditions we want for the Spin-Orbit with density
C=======================================================================
C=======================================================================
C
      IF (KEYWOR.EQ.'IFLAMBVSOD') THEN
C
C         With these parameters we decide how we want to set the 
C         SO-density \lambda's. 
C
C         IFPAR1 - \lambda_pp
C         IFPAR2 - \lambda_pn
C         IFPAR3 - \lambda_np
C         IFPAR4 - \lambda_nn
C
C         Example: if IFPAR1=+1 and IFPAR2=IFPAR3=IFPAR4=-1, then:
C
C         \lambda_pp=\lambda_pn=\lambda_np=\lambda_nn
C
C         If IFPAR1=+1 and IFPAR4=-1 and IFPAR2=+2 and IFPAR3=-2
C
C         \lambda_pp=\lambda_nn
C         \lambda_np=\lambda_pn
C
C         If IFPAR1=IFPAR2=IFPAR3=IFPAR4=+1 the 4 \lambda parameters
C                                                    are independent
C
          READ (N_READ,*) IFPAR1,IFPAR2,IFPAR3,IFPAR4
C
          WRITE(LOGFIL,'(A)') KEYWOR
C
          WRITE(LOGFIL,'(  9X,''Spin-Orbit density dependent '',
     *                        ''conditions: '',
     *                        ''IFPAR1= '',I2,'' IFPAR2= '',I2,1X,
     *                        ''IFPAR3= '',I2,'' IFPAR4= '',I2)')
     *                          IFPAR1,          IFPAR2,
     *                          IFPAR3,          IFPAR4
C
      END IF
C
C=======================================================================
C=======================================================================
C     Which conditions we want for the TENSOR Spin-Orbit with density
C=======================================================================
C=======================================================================
C
      IF (KEYWOR.EQ.'IFLAMBVSOT') THEN
C  @@@ IRENE, DO WE HAVE HERE A LIBERTY OF HAVING INDEPENDENT LAMBDA ???
C         With these parameters we decide how we want to set the 
C         SO-density TENSOR \lambda's. 
C
C         IFPAR5 - \lambda_pp^{t,so}
C         IFPAR6 - \lambda_pn^{t,so}
C         IFPAR7 - \lambda_np^{t,so}
C         IFPAR8 - \lambda_nn^{t,so}
C
C         The functioning is the same as IFLAMBVSOD
C
          READ (N_READ,*) IFPAR5,IFPAR6,IFPAR7,IFPAR8
C
          WRITE(LOGFIL,'(A)') KEYWOR
C
          WRITE(LOGFIL,'(  9X,''Spin-Orbit tensor  dependent '',
     *                        ''conditions: '',
     *                        ''IFPAR5= '',I2,'' IFPAR6= '',I2,1X,
     *                        ''IFPAR7= '',I2,'' IFPAR8= '',I2,''??'')')
     *                          IFPAR5,          IFPAR6,
     *                          IFPAR7,          IFPAR8
C
      END IF
C
C=======================================================================
C=======================================================================
C     Which conditions we want for the TENSOR contribution to central V 
C=======================================================================
C=======================================================================
C  @@@ IRENE, DO WE HAVE HERE A LIBERTY OF HAVING INDEPENDENT LAMBDA ???
      IF (KEYWOR.EQ.'IFLAMBVCNT') THEN
C
C         With these parameters we decide how we want to set the 
C         Central TENSOR \lambda's. 
C
C         IFPA09 - \lambda_pp^{t,cnt}
C         IFPA10 - \lambda_pn^{t,cnt}
C         IFPA11 - \lambda_np^{t,cnt}
C         IFPA12 - \lambda_nn^{t,cnt}
C
C         The functioning is the same as IFLAMBVSOD
C
          READ (N_READ,*) IFPA09,IFPA10,IFPA11,IFPA12
C
          WRITE(LOGFIL,'(A)') KEYWOR
C
          WRITE(LOGFIL,'(  9X,''Central-V  Density dependent '',
     *                        ''conditions: '',
     *                        ''IFPA09= '',I2,'' IFPA10= '',I2,1X,
     *                        ''IFPA11= '',I2,'' IFPA12= '',I2)')
     *                          IFPA09,          IFPA10,
     *                          IFPA11,          IFPA12
C
      END IF   
C
C=======================================================================
C=======================================================================
C     Reading the unit over which we divide the \lambda's for
C     printing/plotting reasons.
C=======================================================================
C=======================================================================
C 
      IF (KEYWOR.EQ.'UNITLAMBDA') THEN
C
          READ(*,*) UNITLA
C
          WRITE(LOGFIL,'(A)') KEYWOR
          WRITE(LOGFIL,'(9X,''Unit over which we divide the '',
     *                      ''\lambdas writing purposes only: UNITLA'',
     *                       /,71X,F6.2)') UNITLA
          DO I=1,5
             WRITE(LOGFIL,'(9X,''IS IT STILL ACTUAL?????'')')
          END DO
      END IF    
C
C=======================================================================
C=======================================================================
C     Reading if we want to use the 'kappa-parametrization' for 
C     depth, radii and diffusness
C=======================================================================
C=======================================================================
C 
      IF (KEYWOR.EQ.'IFKAPPAPAR') THEN
          READ(*,*) IFK_VC,IFK_VS,IFK_RC,IFK_RS,IFK_AC,IFK_AS
          IF (LOGWRI.GT.0) THEN
              WRITE(LOGFIL,'(A)')KEYWOR
              WRITE(LOGFIL,'(9X,''IFK_VC IFK_VS IFK_RC IFK_RS '',
     *                                        ''IFK_AC IFK_AS'')')
              WRITE(LOGFIL,'(9X,6(I6,1X))')IFK_VC,IFK_VS,IFK_RC,
     *                                     IFK_RS,IFK_AC,IFK_AS
          END IF
      END IF
C 
C=======================================================================
C=======================================================================
C     Reading if we want to introduce parametric correlations
C=======================================================================
C=======================================================================
C
      IF (KEYWOR.EQ.'IFPARAMCOR') THEN
          READ(N_READ,*) IFCORR
          IF (LOGWRI.GT.0) THEN
              WRITE(LOGFIL,'(A)')KEYWOR
              WRITE(LOGFIL,'(9X,''IFCORR'')')
              WRITE(LOGFIL,'(9X,I6)')IFCORR
          END IF
      END IF
C 
C=======================================================================
C=======================================================================
C     Reading the indices that activate the parametric correlations
C     that we want to eliminate
C=======================================================================
C=======================================================================
C
      IF (KEYWOR.EQ.'PARAMCORRE') THEN
          READ(N_READ,*) IFRCVC,IFRCAC,IFVCAC,IFRSVS,IFRSAS,IFVSAS
          IF (LOGWRI.GT.0) THEN
              WRITE(LOGFIL,'(A)')KEYWOR
              WRITE(LOGFIL,'(9X,''IFRCVC IFRCAC IFVCAC IFRSVS '',
     *                                        ''IFRSAS IFVSAS'')')
              WRITE(LOGFIL,'(9X,6(I6,1X))')IFRCVC,IFRCAC,IFVCAC,
     *                                     IFRSVS,IFRSAS,IFVSAS
          END IF
      END IF
C 
C=======================================================================
C=======================================================================
C     Central potential parameters - either values to be used in
C     the single run option or starting values, if the parameter
C     has been selected by choosing <IFFITS=1>. If random number
C     option used  ->> the initial parameter values are modified
C=======================================================================
C=======================================================================
C
      IF (KEYWOR.EQ.'CENTRALPOT') THEN
C
C         User's Central Potential Parameter Values
C
          READ (N_READ,*) V0CENT_PROTON,R0CENT_PROTON,A0CENT_PROTON,
     *                    XK_V0C_PROTON,XK_R0C_PROTON,XK_A0C_PROTON,
     *                                         R0COUL,XK_COU,COUFAC,
     *                    V0CENT_NEUTRS,R0CENT_NEUTRS,A0CENT_NEUTRS,
     *                    XK_V0C_NEUTRS,XK_R0C_NEUTRS,XK_A0C_NEUTRS
C 
          WRITE(LOGFIL,'(A)') KEYWOR
C
          WRITE(LOGFIL,'(9X,''User Defined Central Potential '',
     *                      ''Parameter Values:'',/)')
C
          WRITE(LOGFIL,'(9X,''V0CENT  R0CENT  A0CENT  XK_V0C  '',
     *                      ''XK_R0C  XK_A0C  R0COUL  XK_COU'')')
C
          WRITE(LOGFIL,'(8X,F8.4,1X,7(F6.4,2X),2X,''protons'')')
     *                                           V0CENT_PROTON,
     *                             R0CENT_PROTON,A0CENT_PROTON,
     *                             XK_V0C_PROTON,XK_R0C_PROTON,
     *                             XK_A0C_PROTON,R0COUL,XK_COU
C
          WRITE(LOGFIL,'(8X,F8.4,1X,5(F6.4,2X),18X,''neutrons'')')
     *                             V0CENT_NEUTRS,R0CENT_NEUTRS,
     *                             A0CENT_NEUTRS,XK_V0C_NEUTRS,
     *                             XK_R0C_NEUTRS,XK_A0C_NEUTRS          
      END IF          
C 
C=======================================================================
C=======================================================================
C     Spin-Orbit potential parameters
C=======================================================================
C=======================================================================
C
      IF (KEYWOR.EQ.'SPIN-ORBIT') THEN
C
C         User's Traditional WS Spin-Orbit Parameter Values
C
          READ (N_READ,*) V0SORB_PROTON,R0SORB_PROTON,A0SORB_PROTON,
     *                    XK_LAM_PROTON,XK_RSO_PROTON,XK_ASO_PROTON,
     *
     *                    V0SORB_NEUTRS,R0SORB_NEUTRS,A0SORB_NEUTRS,
     *                    XK_LAM_NEUTRS,XK_RSO_NEUTRS,XK_ASO_NEUTRS
C
          IF (IFDENS.EQ.0) THEN
C
              WRITE(LOGFIL,'(A)') KEYWOR
C
              WRITE(LOGFIL,'(9X,''User Defined Traditional WS '',
     *                      ''Spin-Orbit Parameter Values:'',/)')
C
              WRITE(LOGFIL,'(9X,''V0SORB  R0SORB  A0SORB  XK_LAM  '',
     *                                       ''XK_RSO  XK_ASO'')')
C
              WRITE(LOGFIL,'(8X,F8.4,1X,5(F6.4,2X),18X,''protons'')')
     *                    V0SORB_PROTON,R0SORB_PROTON,A0SORB_PROTON,
     *                    XK_LAM_PROTON,XK_RSO_PROTON,XK_ASO_PROTON
C
              WRITE(LOGFIL,'(8X,F8.4,1X,5(F6.4,2X),18X,''neutrons'')')
     *                    V0SORB_NEUTRS,R0SORB_NEUTRS,A0SORB_NEUTRS,
     *                    XK_LAM_NEUTRS,XK_RSO_NEUTRS,XK_ASO_NEUTRS
          END IF
C
      END IF
C 
C=======================================================================
C=======================================================================
C     Kappa parametrization values for Central Potential Parameters
C=======================================================================
C=======================================================================
C
      IF (KEYWOR.EQ.'KAPPA-CENT') THEN
C
          READ(N_READ,*) V0CENT_KAPPAR,XK_V0C_KAPPAR,
     *                   R0CENT_KAPPAR,XK_R0C_KAPPAR,
     *                   A0CENT_KAPPAR,XK_A0C_KAPPAR
C
          WRITE(LOGFIL,'(A)') KEYWOR
C
          IF (IFK_VC.EQ.1) WRITE(LOGFIL,'(9X,''V0CENT  XK_V0C  '',$)')
          IF (IFK_RC.EQ.1) WRITE(LOGFIL,'(''R0CENT  XK_R0C  '',$)')
          IF (IFK_AC.EQ.1) WRITE(LOGFIL,'(''A0CENT  XK_A0C  '')')
C
          IF (IFK_VC.EQ.1) WRITE(LOGFIL,'(9X,2(F6.2,2X),$)')
     *                                                   V0CENT_KAPPAR,
     *                                                   XK_V0C_KAPPAR
C
          IF (IFK_RC.EQ.1) WRITE(LOGFIL,'(2(F6.2,2X),$)')R0CENT_KAPPAR,
     *                                                   XK_R0C_KAPPAR
C
          IF (IFK_AC.EQ.1) WRITE(LOGFIL,'(2(F6.2,2X))')A0CENT_KAPPAR,
     *                                                 XK_A0C_KAPPAR
C
      END IF
C 
C=======================================================================
C=======================================================================
C     Kappa parametrization values for Pure WS Potential Parameters
C=======================================================================
C=======================================================================
C
      IF (KEYWOR.EQ.'KAPPA-SORB') THEN
C
          READ(N_READ,*) V0SORB_KAPPAR,XK_V0S_KAPPAR,
     *                   R0SORB_KAPPAR,XK_R0S_KAPPAR,
     *                   A0SORB_KAPPAR,XK_A0S_KAPPAR
C
          WRITE(LOGFIL,'(A)') KEYWOR
C
          IF (IFK_VS.EQ.1) WRITE(LOGFIL,'(9X,''V0SORB  XK_V0S  '',$)')
          IF (IFK_RS.EQ.1) WRITE(LOGFIL,'(''R0SORB  XK_R0S  '',$)')
          IF (IFK_AS.EQ.1) WRITE(LOGFIL,'(''A0SORB  XK_A0S  '')')
C
          IF (IFK_VS.EQ.1) WRITE(LOGFIL,'(9X,2(F6.2,2X),$)')
     *                                                   V0SORB_KAPPAR,
     *                                                   XK_V0S_KAPPAR
C
          IF (IFK_RS.EQ.1) WRITE(LOGFIL,'(2(F6.2,2X),$)')R0SORB_KAPPAR,
     *                                                   XK_R0S_KAPPAR
C
          IF (IFK_AS.EQ.1) WRITE(LOGFIL,'(2(F6.2,2X))')A0SORB_KAPPAR,
     *                                                 XK_A0S_KAPPAR
C
      END IF
C 
C=======================================================================
C=======================================================================
C     Optional: effective-mass term for the Dirac equation
C=======================================================================
C=======================================================================
C
      IF (KEYWOR.EQ.'EFFECTMASS') THEN
C
C         Effective Mass Parameters
C
          READ (N_READ,*) V0EFFM_PROTON,R0EFFM_PROTON,A0EFFM_PROTON,
     *                    XK_LEF_PROTON,XK_REF_PROTON,XK_AEF_PROTON,
     *                    V0EFFM_NEUTRS,R0EFFM_NEUTRS,A0EFFM_NEUTRS,
     *                    XK_LEF_NEUTRS,XK_REF_NEUTRS,XK_AEF_NEUTRS
C
          WRITE(LOGFIL,'(A)') KEYWOR
C
          WRITE(LOGFIL,'(9X,''User Defined Effective Mass '',
     *                      ''Parameter Values:'',/)')
C
          WRITE(LOGFIL,'(9X,''V0EFFM  R0EFFM  A0EFFM  XK_LEF  '',
     *                                       ''XK_REF  XK_AEF'')')
C
          WRITE(LOGFIL,'(8X,F8.4,1X,5(F6.4,2X),18X,''protons'')')
     *                    V0EFFM_PROTON,R0EFFM_PROTON,A0EFFM_PROTON,
     *                    XK_LEF_PROTON,XK_REF_PROTON,XK_AEF_PROTON
C
          WRITE(LOGFIL,'(8X,F8.4,1X,5(F6.4,2X),18X,''neutrons'')')
     *                    V0EFFM_NEUTRS,R0EFFM_NEUTRS,A0EFFM_NEUTRS,
     *                    XK_LEF_NEUTRS,XK_REF_NEUTRS,XK_AEF_NEUTRS
C
      END IF        
C 
C=======================================================================
C=======================================================================
C     Optional: spin-orbit with density
C=======================================================================
C=======================================================================
C 
      IF (KEYWOR.EQ.'DENSTSORBI') THEN
C
C         Density Spin-Orbit Potential Parameters
C
          READ (N_READ,*) ALAMPP,ALAMPN,ALAMNP,ALAMNN   
C
          IF (IFDENS.EQ.1) THEN
              WRITE(LOGFIL,'(A)')KEYWOR
              WRITE(LOGFIL,'(9X,''User Defined Density Spin-Orbit '',
     *                      ''Potential Parameter Values:'',/)')
              WRITE(LOGFIL,'(9X,''ALAMPP  ALAMPN  ALAMNP  ALAMNN'')')
              WRITE(LOGFIL,'(9X,4(F6.2,2X))')
     *                    ALAMPP,ALAMPN,ALAMNP,ALAMNN
          END IF
C
      END IF    
C 
C=======================================================================
C=======================================================================
C     Optional: spin-orbit with TENSOR density
C=======================================================================
C=======================================================================
C 
      IF (KEYWOR.EQ.'TENSOSORBI') THEN
C
C         Density Spin-Orbit TENSOR Potential Parameters
C
          READ (N_READ,*) TLAMPP,TLAMPN,TLAMNP,TLAMNN 
C
          IF (ISORBT.EQ.1) THEN
              WRITE(LOGFIL,'(A)')KEYWOR
              WRITE(LOGFIL,'(9X,''User Defined Density Spin-Orbit '',
     *                      ''TENSOR Potential Parameter Values:'',/)')
              WRITE(LOGFIL,'(9X,''TLAMPP  TLAMPN  TLAMNP  TLAMNN'')')
              WRITE(LOGFIL,'(9X,4(F6.2,2X))')
     *                    TLAMPP,TLAMPN,TLAMNP,TLAMNN
          END IF
C
      END IF 
C 
C=======================================================================
C=======================================================================
C     Optional: spin-orbit with TENSOR density
C=======================================================================
C=======================================================================
C 
      IF (KEYWOR.EQ.'TENSOCENTR') THEN
C
C         Central TENSOR Potential Parameters
C
          READ (N_READ,*) CLAMPP,CLAMPN,CLAMNP,CLAMNN
C
          IF (ICENTT.EQ.1) THEN
              WRITE(LOGFIL,'(A)')KEYWOR
              WRITE(LOGFIL,'(9X,''User Defined Central TENSOR '',
     *                      ''Potential Parameter Values:'',/)')
              WRITE(LOGFIL,'(9X,''CLAMPP  CLAMPN  CLAMNP  CLAMNN'')')
              WRITE(LOGFIL,'(9X,4(F6.2,2X))')
     *                    CLAMPP,CLAMPN,CLAMNP,CLAMNN
          END IF
      END IF 
C 
C=======================================================================
C=======================================================================
C     Reading the shift applied to V_0_central (only applied for ONERUN)
C=======================================================================
C=======================================================================
C 
      IF (KEYWOR.EQ.'SHIFTCENTR') THEN
C
          READ(N_READ,*) (V0SHIF_PROTON(I),I=1,LDNUCL)
          READ(N_READ,*) (V0SHIF_NEUTRS(I),I=1,LDNUCL)
C
          IF (IFTEST.EQ.1 .AND. LOGWRI.GT.0) THEN
              WRITE(LOGFIL,'(A)')KEYWOR
              WRITE(LOGFIL,'(9X,''SHIFOx  SHIFC0  SHIFC8  SHIFNi  '',
     *                          ''SHIFZr  SHIFSn  SHIFGd  SHIFPb'')')
              WRITE(LOGFIL,'(9X,<LDNUCL>(F6.2,2X),/,
     *                       9X,<LDNUCL>(F6.2,2X))')
     *                      (V0SHIF_PROTON(I),I=1,LDNUCL),
     *                      (V0SHIF_NEUTRS(I),I=1,LDNUCL)
          END IF
C
      END IF
C      
C=======================================================================
C      < Reading the control parameters for the minimisation >
C=======================================================================
C
C      Defining the Dirac parameters valid for the whole run:
C
C      IRMFST (I - Relativistic Mean Field - Strict)  controls 
C              the equality  between spin-orbit parameters and 
C              effective mass parameters. According to the RMF, 
C              both sets should be equal.  This relation is to 
C              forced when IRMFST=1. In such a case, we cannot
C              vary the effective mass parameters  (see NAMELI, 
C              where  independent  minimisation variables  are 
C                                                      defined
C
C      I_RAND - It should be equal to 1 if minimisation should
C               start from randomly chosen initial values
C
C      LDRAND - The number  of restarts  of the  random-number
C               generating routine and thus also the number of
C               independent minimisation - check LDRAND
C
C      I_SEED - First we initialise the "seed" by calling 
C
C                        RANDIN(I_SEED), 
C             
C               where I_SEED can be set to any value.
C 
C      Next generate your random numbers (uniformly distributed 
C      within 0 and 1) by calling ZUFALL(N,A);  N is the actual 
C      number of the random numbers needed  and A is the vector 
C      that contains them
C
C=======================================================================
C 
      IF (KEYWOR.EQ.'WHICHMINIM') THEN
C
          READ (*,*) IRMFST,I_RAND,LDRAND,I_SEED,NEWSED,IMTALK,IRUNMI
C
          WRITE(LOGFIL,'(A)') KEYWOR
C
          WRITE(LOGFIL,'(9X,''IRMFST  I_RAND  LDRAND  '',
     *                      ''I_SEED  NEWSED  IMTALK  IRUNMI'')')
C
          WRITE(LOGFIL,'(9X,7(I6,2X))') IRMFST,I_RAND,LDRAND,I_SEED,
     *                                         NEWSED,IMTALK,IRUNMI
          IF (LDRAND.GT.NDRAND) THEN
              WRITE(NOUTPT,'(''LDRAND='',I6,'' exceeds NDRAND='',I6)')
     *                         LDRAND,                 NDRAND
              STOP 'LDRAND > NDRAND in NAMELI'
          END IF
      END IF
C
C=======================================================================
C=======================================================================
C     Reading ITECHI - the maximum no. of chi^2 evaluations 
C=======================================================================
C=======================================================================
C 
      IF (KEYWOR.EQ.'MINIMEVALS') THEN
C
          READ (N_READ,*) ITECHI ! Maximum number of chi^2 evaluations
C
          IF (IFFITS.EQ.1 .OR. IFMESH.EQ.1) THEN
C
              WRITE(LOGFIL,'(A)') KEYWOR
              WRITE(LOGFIL,'(9X,''Maximum number of chi^2 evaluations'',
     *                          '' ITECHI= '',I5)') ITECHI
          END IF
C
      END IF
C
C=======================================================================
C     ... and also the tolerance factors
C
      IF (KEYWOR.EQ.'TOLERANCES') THEN
C
          READ (*,*) DWFACT,UPFACT,PUSHIN
C
          IF (IFFITS.EQ.1 .OR. IFMESH.EQ.1) THEN
C
              WRITE(LOGFIL,'(A)') KEYWOR
C
              WRITE(LOGFIL,'(9X,''chi^2 Tolerances: '',
     *                  ''DWFACT  UPFACT  PUSHIN'',/,28X,3(f6.4,2x))') 
     *                    DWFACT, UPFACT, PUSHIN
          END IF
      END IF
C
C=======================================================================
C     ... and the conditions for stoping the minimisation (for LEVMAR)
C
      IF (KEYWOR.EQ.'STOPTOLERS') THEN
C  @@@ IRENE: WHAT IS FACTOR ???
          READ(*,*) TOLERF,TOLERX,TOLERG,FACTOR
C
          IF (IFFITS.EQ.1 .OR. IFMESH.EQ.1) THEN
C
              WRITE(LOGFIL,'(A)') KEYWOR
C
              WRITE(LOGFIL,'(9X,''TOLERF  TOLERX  TOLERG  FACTOR'')')
C
              WRITE(LOGFIL,'(9X,3(ES6.0,2X),F6.2)') TOLERF,TOLERX,
     *                                              TOLERG,FACTOR
          END IF
C
      END IF
C
C=======================================================================
C     ... and the conditions for stoping the minimisation (for LEVMAR)
C 
      IF (KEYWOR.EQ.'STOPINCOND') THEN
C
C         If the gradient condition does not impo
C         LDLAST - We continue iterating at most LDSAT times verifying
C                  whether 
C
          READ(*,*) LDLAST,EPSLAS
C
          IF (IFFITS.EQ.1 .OR. IFMESH.EQ.1) THEN
C
              WRITE(LOGFIL,'(A)') KEYWOR
C
              WRITE(LOGFIL,'(9X,''Stopping conditions for LEVMAR: '',
     *                      ''LDLAST  EPSLAS'',/,41X,I6,2X,ES6.0)')
     *                        LDLAST, EPSLAS
          END IF
C
          IF (LDLAST.GT.NDLAST) THEN
C
              WRITE(LOGFIL,'(/,''Alarm in NAMELI: LDLAST= '',I3,
     *                  '' is .GT. NDLAST= '',I3,/)') LDLAST,NDLAST

              WRITE(0,'(/,''Alarm in NAMELI: LDLAST= '',I3,
     *                  '' is .GT. NDLAST= '',I3,/)') LDLAST,NDLAST
C
              STOP 'STOP in NAMELI: LDLAST.GT.NDLAST'
C
          END IF
C
      END IF
C
C=======================================================================
C=======================================================================
C     Reading the weight-factors to define the minimised function
C=======================================================================
C=======================================================================
C 
      IF (KEYWOR.EQ.'MINWEIGHTS') THEN
C  @@@ IRENE; EVEN IF ONE CAN GUESS THE MEANING OF MOST OF THE WEIGHTS,
C      PLEASE INTRODUCE THE DEFINITIONS (IN A FEW WORDS) OF ALL OF THEM 
          READ (*,*) WEICOR,WEIRAD,WEIINV,WEIFER,WEIGAP,WEIWEI,
     *               WEIDIF,WEIABS,WEIMAX,WEIDUP,WEIDDW,WEIRHO
C
C ###     Renaming the WEIGHTS
C       
          WEIGHT_CORREL=WEICOR
          WEIGHT_RADIUS=WEIRAD
          WEIGHT_INVERT=WEIINV
          WEIGHT_EFERMI=WEIFER
          WEIGHT_ENEGAP=WEIGAP
C
C ###     WEIWEI rests with the same name
C
          WEIGHT_ERRABS=WEIDIF
          WEIGHT_EABSAV=WEIABS
          WEIGHT_ERRMAX=WEIMAX
          WEIGHT_DENSUP=WEIDUP
          WEIGHT_DENSDW=WEIDDW   
          WEIGHT_RHODEN=WEIRHO   
C
          WRITE(LOGFIL,'(A)')KEYWOR
C
          WRITE(LOGFIL,'(9X,''WEICOR  WEIRAD  WEIINV  WEIFER  '',
     *                      ''WEIGAP  WEIWEI  WEIDIF  WEIABS  '',
     *                      ''WEIMAX  WEIDUP  WEIDDW  WEIRHO'')')
C
          WRITE(LOGFIL,'(9X,12(F6.2,2X))') WEICOR,WEIRAD,WEIINV,
     *                                     WEIFER,WEIGAP,WEIWEI,
     *                WEIDIF,WEIABS,WEIMAX,WEIDUP,WEIDDW,WEIRHO
     
      END IF
C
C=======================================================================
C=======================================================================
C     Reading the choice parameter that will define 'chi^2'
C=======================================================================
C=======================================================================
C 
      IF (KEYWOR.EQ.'CHOICEOFHI') THEN
C  @@@ IRENE: PLEAS PUT THE COMMENTS ABOUT ALL POSSIBLE CHOICES
C      AND TEHIR MEANING
          READ (*,*) CHIDEF
C
          WRITE(LOGFIL,'(A)') KEYWOR
          WRITE(LOGFIL,'(9X,''CHIDEF= '',A)') CHIDEF
C
      END IF
C
C=======================================================================
C=======================================================================
C     Reading the weight factors for energies and radii for each nucleus
C=======================================================================
C=======================================================================
C 
      IF (KEYWOR.EQ.'CHIWEI_PRO') THEN
                                           LDNUCL_AUXILI=8
          READ (*,*) (WEINUC_PROTON(I),I=1,LDNUCL_AUXILI)
C
          WRITE(LOGFIL,'(A)') KEYWOR
C
          WRITE(LOGFIL,'(9X,''WEI_Ox  WEI_Ca  WEI_Ca  WEI_Ni  '',
     *                      ''WEI_Zr  WEI_Sn  WEI_Gd  WEI_Pb'')')
C
          WRITE(LOGFIL,'(9X,<LDNUCL_AUXILI>(F6.2,2X))')
     *                              (WEINUC_PROTON(I),I=1,LDNUCL_AUXILI)
C          
      END IF
C
C=======================================================================
C    
      IF (KEYWOR.EQ.'CHIWEI_NEU') THEN
                                           LDNUCL_AUXILI=8
          READ (*,*) (WEINUC_NEUTRS(I),I=1,LDNUCL_AUXILI)
C
          WRITE(LOGFIL,'(A)') KEYWOR
C
          WRITE(LOGFIL,'(9X,''WEI_Ox  WEI_Ca  WEI_Ca  WEI_Ni  '',
     *                      ''WEI_Zr  WEI_Sn  WEI_Gd  WEI_Pb'')')
C
          WRITE(LOGFIL,'(9X,<LDNUCL_AUXILI>(F6.2,2X))')
     *                              (WEINUC_NEUTRS(I),I=1,LDNUCL_AUXILI)
C          
      END IF
C
C=======================================================================
C      
      IF (KEYWOR.EQ.'RADWEI_PRO') THEN
                                           LDNUCL_AUXILI=8
          READ (*,*) (WEIRAD_PROTON(I),I=1,LDNUCL_AUXILI)
C
          WRITE(LOGFIL,'(A)') KEYWOR
C
          WRITE(LOGFIL,'(9X,''WEI_Ox  WEI_Ca  WEI_Ca  WEI_Ni  '',
     *                      ''WEI_Zr  WEI_Sn  WEI_Gd  WEI_Pb'')')
C
          WRITE(LOGFIL,'(9X,<LDNUCL_AUXILI>(F6.2,2X))')
     *                              (WEIRAD_PROTON(I),I=1,LDNUCL_AUXILI)
C          
      END IF
C
C=======================================================================
C   
      IF (KEYWOR.EQ.'RADWEI_NEU') THEN
                                           LDNUCL_AUXILI=8
          READ (*,*) (WEIRAD_NEUTRS(I),I=1,LDNUCL_AUXILI)
C
          WRITE(LOGFIL,'(A)') KEYWOR
C
          WRITE(LOGFIL,'(9X,''WEI_Ox  WEI_Ca  WEI_Ca  WEI_Ni  '',
     *                      ''WEI_Zr  WEI_Sn  WEI_Gd  WEI_Pb'')')

          WRITE(LOGFIL,'(9X,<LDNUCL_AUXILI>(F6.2,2X))')
     *                              (WEIRAD_NEUTRS(I),I=1,LDNUCL_AUXILI)
C          
      END IF
C 
C=======================================================================
C=======================================================================
C
      IF (KEYWOR.EQ.'RMSVAL_PRO') THEN
                                          LDNUCL_AUXILI=8
          READ(*,*) (RMSVAL_REFPRO(I),I=1,LDNUCL_AUXILI)
C 
          WRITE(LOGFIL,'(A)') KEYWOR
C
          WRITE(LOGFIL,'(9X,''Pro_Ox  Pro_Ca  Pro_Ca  Pro_Ni  '',
     *                      ''Pro_Zr  Pro_Sn  Pro_Gd  Pro_Pb'')')

          WRITE(LOGFIL,'(9X,<LDNUCL_AUXILI>(F6.2,2X))')
     *                              (RMSVAL_REFPRO(I),I=1,LDNUCL_AUXILI)
      END IF
C
C=======================================================================
C
      IF (KEYWOR.EQ.'RMSVAL_NEU') THEN
                                          LDNUCL_AUXILI=8
          READ(*,*) (RMSVAL_REFNEU(I),I=1,LDNUCL_AUXILI)
C 
          WRITE(LOGFIL,'(A)') KEYWOR
C
          WRITE(LOGFIL,'(9X,''Neu_Ox  Neu_Ca  Neu_Ca  Neu_Ni  '',
     *                      ''Neu_Zr  Neu_Sn  Neu_Gd  Neu_Pb'')')

          WRITE(LOGFIL,'(9X,<LDNUCL_AUXILI>(F6.2,2X))')
     *                              (RMSVAL_REFNEU(I),I=1,LDNUCL_AUXILI)
      END IF
C
C=======================================================================
C
      IF (KEYWOR.EQ.'RMSGLO_REF') THEN
C
          READ(*,*) RMSGLO_REFPRO,RMSGLO_REFNEU
C 
          WRITE(LOGFIL,'(A)') KEYWOR
C
          WRITE(LOGFIL,'(9X,''RMSPRO  RMSNEU'')')

          WRITE(LOGFIL,'(9X,<LDNUCL>(F6.2,2X))')RMSGLO_REFPRO, 
     *                                          RMSGLO_REFNEU
      END IF
C 
C=======================================================================
C=======================================================================
C
      IF (KEYWOR.EQ.'RMSIND_PRO') THEN
                                          LDNUCL_AUXILI=8
          READ(*,*) (RMSIND_PROTON(I),I=1,LDNUCL_AUXILI)
C 
          WRITE(LOGFIL,'(A)') KEYWOR
C
          WRITE(LOGFIL,'(9X,''Pro_Ox  Pro_Ca  Pro_Ca  Pro_Ni  '',
     *                      ''Pro_Zr  Pro_Sn  Pro_Gd  Pro_Pb'')')

          WRITE(LOGFIL,'(9X,<LDNUCL_AUXILI>(F6.2,2X))')
     *                              (RMSIND_PROTON(I),I=1,LDNUCL_AUXILI)
      END IF
C
C=======================================================================
C
      IF (KEYWOR.EQ.'RMSIND_NEU') THEN
                                          LDNUCL_AUXILI=8
          READ(*,*) (RMSIND_NEUTRS(I),I=1,LDNUCL_AUXILI)
C 
          WRITE(LOGFIL,'(A)') KEYWOR
C
          WRITE(LOGFIL,'(9X,''Neu_Ox  Neu_Ca  Neu_Ca  Neu_Ni  '',
     *                      ''Neu_Zr  Neu_Sn  Neu_Gd  Neu_Pb'')')

          WRITE(LOGFIL,'(9X,<LDNUCL_AUXILI>(F6.2,2X))')
     *                              (RMSIND_PROTON(I),I=1,LDNUCL_AUXILI)
      END IF
C 
C=======================================================================
C=======================================================================
C     Central  potential  parameters - ranges  for minimisation
C=======================================================================
C=======================================================================
C
      IF (KEYWOR.EQ.'CENTLIMITS') THEN
C
          READ (N_READ,*) V0CMIN,V0CMAX,A0CMIN,A0CMAX,R0CMIN,R0CMAX,
     *                                                CouMIN,CouMAX
          WRITE(LOGFIL,'(A)') KEYWOR
C
          WRITE(LOGFIL,'(9X,''V0CMIN  V0CMAX  '',
     *                      ''A0CMIN  A0CMAX  '',
     *                      ''R0CMIN  R0CMAX'',/,
     *                   9X,  F6.2,2X,F6.2,2X,
     *                        F6.4,2X,F6.4,2X,
     *                        F6.4,2X,F6.4,2X)')
     *
     *                        V0CMIN,V0CMAX,A0CMIN,A0CMAX,R0CMIN,R0CMAX
C
      END IF          
C 
C=======================================================================
C=======================================================================
C     Spin-orbit potential parameters - ranges for minimisation
C=======================================================================
C=======================================================================
C
      IF (KEYWOR.EQ.'SORBLIMITS') THEN
C
          READ (N_READ,*) XL_MIN,XL_MAX,A0SMIN,A0SMAX,R0SMIN,R0SMAX
C
          IF (IFDENS.EQ.0) THEN
C
             WRITE(LOGFIL,'(A)') KEYWOR
C
             WRITE(LOGFIL,'(9X,''XL_MIN  XL_MAX  '',
     *                         ''A0SMIN  A0SMAX  '',
     *                         ''R0SMIN  R0SMAX'',/,9X,
     *                           F6.2,2X,F6.2,2X,
     *                           F6.4,2X,F6.4,2X,
     *                           F6.4,2X,F6.4,2X)')
     *
     *                    XL_MIN,XL_MAX,A0SMIN,A0SMAX,R0SMIN,R0SMAX
          END IF
C
      END IF         
C 
C=======================================================================
C=======================================================================
C     Effective-mass   term  parameters - ranges for minimisation
C=======================================================================
C=======================================================================
C
      IF (KEYWOR.EQ.'EFFMLIMITS') THEN
C
          READ (N_READ,*) XEFMIN,XEFMAX,AEFMIN,AEFMAX,REFMIN,REFMAX
C
          WRITE(LOGFIL,'(A)')KEYWOR
C               
          WRITE(LOGFIL,'(9X,''XEFMIN  XEFMAX  '',
     *                      ''AEFMIN  AEFMAX  '',
     *                      ''REFMIN  REFMAX'',/,9X,
     *                        F6.2,2X,F6.2,2X,
     *                        F6.4,2X,F6.4,2X,
     *                        F6.4,2X,F6.4,2X)')
     *
     *                        XEFMIN,XEFMAX,AEFMIN,AEFMAX,REFMIN,REFMAX
C
      END IF          
C 
C=======================================================================
C=======================================================================
C     Density dependent S-O pot. parameters - ranges for minimisation
C=======================================================================
C=======================================================================
C 
      IF (KEYWOR.EQ.'DENSLIMITS') THEN
C
          READ (N_READ,*) XPPMIN,XPPMAX,XPNMIN,XPNMAX,
     *                    XNPMIN,XNPMAX,XNNMIN,XNNMAX
C
          IF (IFDENS.EQ.1) THEN
C
             WRITE(LOGFIL,'(A)') KEYWOR
C               
             WRITE(LOGFIL,'(9X,''XPPMIN  XPPMAX  '',
     *                         ''XPNMIN  XPNMAX'',/,9X,
     *                           F6.1,2X,F6.1,2X,
     *                           F6.1,2X,F6.1,2X,//,9X,
     *                         ''XNPMIN  XNPMAX  '',
     *                         ''XNNMIN  XNNMAX'',/,9X,
     *                           F6.1,2X,F6.1,2X,
     *                           F6.1,2X,F6.1,2X)')
     *
     *                           XPPMIN,XPPMAX,XPNMIN,XPNMAX,
     *                           XNPMIN,XNPMAX,XNNMIN,XNNMAX
          END IF
C
      END IF
C 
C=======================================================================
C=======================================================================
C     Tensor parameters - ranges for minimisation
C=======================================================================
C=======================================================================
C 
      IF (KEYWOR.EQ.'TENSLIMITS') THEN
C
          READ (N_READ,*) YPPMIN,YPPMAX,YPNMIN,YPNMAX,
     *                    YNPMIN,YNPMAX,YNNMIN,YNNMAX
C
          IF (ISORBT.EQ.1 .AND. IFTENS.EQ.1) THEN
C
              WRITE(LOGFIL,'(A)')KEYWOR
C               
              WRITE(LOGFIL,'(9X,''YPPMIN  YPPMAX  '',
     *                          ''YPNMIN  YPNMAX'',/,9X,
     *                            F6.1,2X,F6.1,2X,
     *                            F6.1,2X,F6.1,2X,//,9X,
     *                          ''YNPMIN  YNPMAX  '',
     *                          ''YNNMIN  YNNMAX'',/,9X,
     *                            F6.1,2X,F6.1,2X,
     *                            F6.1,2X,F6.1,2X)')
     *
     *                            YPPMIN,YPPMAX,YPNMIN,YPNMAX,
     *                            YNPMIN,YNPMAX,YNNMIN,YNNMAX  
          END IF    
      END IF          
C 
C=======================================================================
C=======================================================================
C     Tensor parameters - ranges for minimisation
C=======================================================================
C=======================================================================
C
      IF (KEYWOR.EQ.'TNCNLIMITS') THEN
C
          READ (N_READ,*) CPPMIN,CPPMAX,CPNMIN,CPNMAX,
     *                    CNPMIN,CNPMAX,CNNMIN,CNNMAX
C
          IF (ICENTT.EQ.1 .AND. IFTENS.EQ.1) THEN
C
              WRITE(LOGFIL,'(A)') KEYWOR
C
              WRITE(LOGFIL,'(9X,''CPPMIN  CPPMAX  '',
     *                          ''CPNMIN  CPNMAX'',/,9X,
     *                            F6.1,2X,F6.1,2X,
     *                            F6.1,2X,F6.1,2X,//,9X,
     *                          ''CNPMIN  CNPMAX  '',
     *                          ''CNNMIN  CNNMAX'',/,9X,
     *                            F6.1,2X,F6.1,2X,
     *                            F6.1,2X,F6.1,2X)')
     *
     *                            CPPMIN,CPPMAX,CPNMIN,CPNMAX,
     *                            CNPMIN,CNPMAX,CNNMIN,CNNMAX
          END IF    
C
      END IF  
C 
C=======================================================================
C=======================================================================
C     Choice of the argument increments for the derivatives
C=======================================================================
C=======================================================================
C
      IF (KEYWOR.EQ.'DERPARSTPS') THEN
C
          READ (N_READ,*) DL_V0C,DL_A0C,DL_R0C,DL_XLA,DL_A0S,DL_R0S,
     *                                         DL_XRM,DL_ARM,DL_RRM,
     *                                                       DL_COU,
     *                                  DL_XPP,DL_XPN,DL_XNP,DL_XNN,
     *                                  DL_YPP,DL_YPN,DL_YNP,DL_YNN,
     *                                  DL_CPP,DL_CPN,DL_CNP,DL_CNN
C
          WRITE(LOGFIL,'(A)') KEYWOR
C
          WRITE(LOGFIL,'(9X,''DL_V0C  DL_A0C  DL_R0C'')')
          WRITE(LOGFIL,'(9X,3(F6.4,2X),/)')DL_V0C,DL_A0C,DL_R0C
C
          IF (IFDENS.EQ.0) THEN
              WRITE(LOGFIL,'(9X,''DL_XLA  DL_A0S  DL_R0S'')')
              WRITE(LOGFIL,'(9X,3(F6.4,2X),/)')DL_XLA,DL_A0S,DL_R0S
          END IF
C
          WRITE(LOGFIL,'(9X,''DL_XRM  DL_ARM  DL_RRM  DL_COU'')')
          WRITE(LOGFIL,'(9X,4(F6.4,2X),/)')DL_XRM,DL_ARM,DL_RRM,DL_COU
C
          IF (IFDENS.EQ.1) THEN
C
              WRITE(LOGFIL,'(9X,''DL_XPP  DL_XPN  DL_XNP  DL_XNN'')')
              WRITE(LOGFIL,'(9X,4(F6.4,2X))')
     *                                    DL_XPP,DL_XPN,DL_XNP,DL_XNN
C
              IF (ISORBT.EQ.1 .AND. IFTENS.EQ.1) THEN
C
                  WRITE(LOGFIL,'(9X,''DL_YPP  DL_YPN  DL_YNP  '',
     *                                              ''DL_YNN'')')
                  WRITE(LOGFIL,'(9X,4(F6.4,2X))')
     *                                   DL_YPP,DL_YPN,DL_YNP,DL_YNN
              END IF
C
              IF (ICENTT.EQ.1 .AND. IFTENS.EQ.1) THEN
C
                  WRITE(LOGFIL,'(9X,''DL_CPP  DL_CPN  DL_CNP  '',
     *                                              ''DL_CNN'')')
                  WRITE(LOGFIL,'(9X,4(F6.4,2X))')
     *                                  DL_CPP,DL_CPN,DL_CNP,DL_CNN
              END IF
C
          END IF   
C_______________________________________________________________________
C
          WRITE(NOUTPT,'(/,80(''#''),/,''#'',T80,''#'',/,
     *                ''#  V0CMIN  V0CMAX  DL_V0C    '',
     *                ''A0CMIN  A0CMAX  DL_A0C    '',
     *                ''R0CMIN  R0CMAX  DL_R0C  #'',/,
     *                ''#  '',F6.2,2X,F6.2,2X,F6.4,1X,
     *                ''   '',F6.4,2X,F6.4,2X,F6.4,1X,
     *                ''   '',F6.4,2X,F6.4,2X,F6.4,''  #'')')
     *
     *                V0CMIN,V0CMAX,DL_V0C,A0CMIN,A0CMAX,DL_A0C,
     *                                     R0CMIN,R0CMAX,DL_R0C
C_______________________________________________________________________
C
          WRITE(NOUTPT,'(  ''#'',T80,''#'',/,
     *                ''#  XL_MIN  XL_MAX  DL_XLA    '',
     *                ''A0SMIN  A0SMAX  DL_A0S    '',
     *                ''R0SMIN  R0SMAX  DL_R0S  #'',/,
     *                ''#  '',F6.2,2X,F6.2,2X,F6.4,1X,
     *                ''   '',F6.4,2X,F6.4,2X,F6.4,1X,
     *                ''   '',F6.4,2X,F6.4,2X,F6.4,''  #'')')
     *
     *                XL_MIN,XL_MAX,DL_XLA,A0SMIN,A0SMAX,DL_A0S,
     *                R0SMIN,R0SMAX,DL_R0S
C_______________________________________________________________________
C
          WRITE(NOUTPT,'(  ''#'',T80,''#'',/,
     *                ''#  XEFMIN  XEFMAX  DL_XRM    '',
     *                ''AEFMIN  AEFMAX  DL_ARM    '',
     *                ''REFMIN  REFMAX  DL_RRM  #'',/,
     *                ''#  '',F6.2,2X,F6.2,2X,F6.4,1X,
     *                ''   '',F6.4,2X,F6.4,2X,F6.4,1X,
     *                ''   '',F6.4,2X,F6.4,2X,F6.4,''  #'',/,
     *                ''#'',T80,''#'')')
     *
     *                XEFMIN,XEFMAX,DL_XRM,AEFMIN,AEFMAX,DL_ARM,
     *                REFMIN,REFMAX,DL_RRM
C_______________________________________________________________________
C
          IF (IFDENS.EQ.1) THEN
C
              WRITE(NOUTPT,'(  ''#'',T80,''#'',/,
     *                ''#  XPPMIN  XPPMAX  DL_XPP    '',
     *                ''XPNMIN  XPNMAX  DL_XPN'',T80,''#'',/,
     *                ''#  '',F6.1,2X,F6.2,2X,F6.4,1X,
     *                ''   '',F6.1,2X,F6.2,2X,F6.4,1X,
     *                T80,''#'',/,''#'',T80,''#'',/,
     *                ''#  XNPMIN  XNPMAX  DL_XNP    '',
     *                ''XNNMIN  XNNMAX  DL_XNN'',T80,''#'',
     *                /,
     *                ''#  '',F6.1,2X,F6.2,2X,F6.4,1X,
     *                ''   '',F6.1,2X,F6.2,2X,F6.4,1X,
     *                T80,''#'',/,
     *                ''#'',T80,''#'',/,80(''#''))')
     *
     *                XPPMIN,XPPMAX,DL_XPP,XPNMIN,XPNMAX,DL_XPN,
     *                XNPMIN,XNPMAX,DL_XNP,XNNMIN,XNNMAX,DL_XNN  
C_______________________________________________________________________
C
              WRITE(NOUTPT,'(  ''#'',T80,''#'',/,
     *                ''#  YPPMIN  YPPMAX  DL_YPP    '',
     *                ''YPNMIN  YPNMAX  DL_YPN'',T80,''#'',/,
     *                ''#  '',F6.1,2X,F6.2,2X,F6.4,1X,
     *                ''   '',F6.1,2X,F6.2,2X,F6.4,1X,
     *                T80,''#'',/,''#'',T80,''#'',/,
     *                ''#  YNPMIN  YNPMAX  DL_YNP    '',
     *                ''YNNMIN  YNNMAX  DL_YNN'',T80,''#'',
     *                /,
     *                ''#  '',F6.1,2X,F6.2,2X,F6.4,1X,
     *                ''   '',F6.1,2X,F6.2,2X,F6.4,1X,
     *                T80,''#'',/,
     *                ''#'',T80,''#'',/,80(''#''))')
     *
     *                YPPMIN,YPPMAX,DL_YPP,YPNMIN,YPNMAX,DL_YPN,
     *                YNPMIN,YNPMAX,DL_YNP,YNNMIN,YNNMAX,DL_YNN 
C_______________________________________________________________________
C
              WRITE(NOUTPT,'(  ''#'',T80,''#'',/,
     *                ''#  CPPMIN  CPPMAX  DL_CPP    '',
     *                ''CPNMIN  CPNMAX  DL_CPN'',T80,''#'',/,
     *                ''#  '',F6.1,2X,F6.1,2X,F6.4,1X,
     *                ''   '',F6.1,2X,F6.1,2X,F6.4,1X,
     *                T80,''#'',/,''#'',T80,''#'',/,
     *                ''#  CNPMIN  CNPMAX  DL_CNP    '',
     *                ''CNNMIN  CNNMAX  DL_CNN'',T80,''#'',
     *                /,
     *                ''#  '',F6.1,2X,F6.1,2X,F6.4,1X,
     *                ''   '',F6.1,2X,F6.1,2X,F6.4,1X,
     *                T80,''#'',/,
     *                ''#'',T80,''#'',/,80(''#''))')
     *
     *                CPPMIN,CPPMAX,DL_CPP,CPNMIN,CPNMAX,DL_CPN,
     *                CNPMIN,CNPMAX,DL_CNP,CNNMIN,CNNMAX,DL_CNN  
C
          END IF
C  
      END IF        
C 
C=======================================================================
C=======================================================================
C     Central potential KAPPA parameters - ranges for minimisation
C=======================================================================
C=======================================================================
C 
      IF (KEYWOR.EQ.'CENTKAPPAS') THEN
C
          READ (N_READ,*) XKCMIN,XKCMAX,XACMIN,XACMAX,XRCMIN,XRCMAX,
     *                                                XCoMIN,XCoMAX
C
          WRITE(LOGFIL,'(A)') KEYWOR
C
          WRITE(LOGFIL,'(9X,''XKCMIN  XKCMAX  '',
     *                      ''XACMIN  XACMAX  '',
     *                      ''XRCMIN  XRCMAX  '',
     *                      ''XCoMIN  XCoMAX'',/,
     *                        9X,F6.2,2X,F6.2,2X,
     *                           F6.3,2X,F6.4,2X,
     *                           F6.3,2X,F6.4,2X,
     *                           F6.3,2X,F6.4)')
     *
     *                        XKCMIN,XKCMAX,XACMIN,XACMAX,XRCMIN,XRCMAX,
     *                                                    XCoMIN,XCoMAX
C
      END IF          
C 
C=======================================================================
C=======================================================================
C     Spin-orbit potential KAPPA parameters - ranges for minimisation
C=======================================================================
C=======================================================================
C 
      IF (KEYWOR.EQ.'SORBKAPPAS') THEN
C
          READ (N_READ,*) XKSMIN,XKSMAX,XASMIN,XASMAX,XRSMIN,XRSMAX
C
          IF (IFDENS.EQ.0) THEN

              WRITE(LOGFIL,'(A)') KEYWOR

              WRITE(LOGFIL,'(9X,''XKSMIN  XKSMAX  '',
     *                          ''XASMIN  XASMAX  '',
     *                          ''XRSMIN  XRSMAX  '',/,
     *                            9X,F6.2,2X,F6.2,2X,
     *                               F6.3,2X,F6.4,2X,
     *                               F6.3,2X,F6.4,2X)')
     *
     *                      XKSMIN,XKSMAX,XASMIN,XASMAX,XRSMIN,XRSMAX
          END IF    
C
      END IF          
C 
C=======================================================================
C=======================================================================
C     Reduced-mass term KAPPA parameters - ranges for minimisation
C=======================================================================
C=======================================================================
C 
      IF (KEYWOR.EQ.'EFFMKAPPAS') THEN
C
          READ (N_READ,*) XKEMIN,XKEMAX,XAEMIN,XAEMAX,XREMIN,XREMAX
C
          WRITE(LOGFIL,'(A)') KEYWOR
C
          WRITE(LOGFIL,'(9X,''XKEMIN  XKEMAX  '',
     *                      ''XAEMIN  XAEMAX  '',
     *                      ''XREMIN  XREMAX  '',/,
     *                        9X,F6.2,2X,F6.2,2X,
     *                           F6.3,2X,F6.4,2X,
     *                           F6.3,2X,F6.4,2X)')
     *
     *                        XKEMIN,XKEMAX,XAEMIN,XAEMAX,XREMIN,XREMAX
C
      END IF          
C 
C=======================================================================
C=======================================================================
C     Choice of the argument increments (KAPPA-old) for the derivatives
C=======================================================================
C=======================================================================
C
      IF (KEYWOR.EQ.'DERKAPSTPS') THEN
C
          READ (N_READ,*) DLKV0C,DLKA0C,DLKR0C,DLKXLA,DLKA0S,DLKR0S,
     *                                         DLKXRM,DLKARM,DLKRRM,
     *                                                       DLKCOU
          WRITE(LOGFIL,'(A)') KEYWOR
C
          WRITE(LOGFIL,'(9X,''DLKV0C  DLKA0C  DLKR0C  DLKXLA  '',
     *                      ''DLKA0S  DLKR0S  DLKXRM  DLKARM  '',
     *                      ''DLKRRM  DLKCOU'')')
C
          WRITE(LOGFIL,'(9X,10(F6.4,2X))') DLKV0C,DLKA0C,DLKR0C,DLKXLA,
     *                                     DLKA0S,DLKR0S,DLKXRM,DLKARM,
     *                                                   DLKRRM,DLKCOU
C_______________________________________________________________________
C
          WRITE(NOUTPT,'(/,80(''#''),/,''#'',T80,''#'',/,
     *                   ''#  XKCMIN  XKCMAX  DLKV0C    '',
     *                   ''XACMIN  XACMAX  DLKA0C    '',
     *                   ''XRCMIN  XRCMAX  DLKR0C  #'',/,
     *                   ''#  '',F6.2,2X,F6.2,2X,F6.4,1X,
     *                   ''   '',F6.3,2X,F6.4,2X,F6.4,1X,
     *                   ''   '',F6.3,2X,F6.4,2X,F6.4,''  #'')')
     *
     *                   XKCMIN,XKCMAX,DLKV0C,XACMIN,XACMAX,DLKA0C,
     *                                        XRCMIN,XRCMAX,DLKR0C
C_______________________________________________________________________
C
          WRITE(NOUTPT,'(''#'',T80,''#'',/,
     *                ''#  XKSMIN  XKSMAX  DLKXLA    '',
     *                ''XASMIN  XASMAX  DLKA0S    '',
     *                ''XRSMIN  XRSMAX  DLKR0S  #'',/,
     *                ''#  '',F6.2,2X,F6.2,2X,F6.4,1X,
     *                ''   '',F6.3,2X,F6.4,2X,F6.4,1X,
     *                ''   '',F6.3,2X,F6.4,2X,F6.4,''  #'')')
     *
     *                XKSMIN,XKSMAX,DLKXLA,XASMIN,XASMAX,DLKA0S,
     *                                     XRSMIN,XRSMAX,DLKR0S
C_______________________________________________________________________
C 
          WRITE(NOUTPT,'(  ''#'',T80,''#'',/,
     *                ''#  XKEMIN  XKEMAX  DLKXRM    '',
     *                ''XAEMIN  XAEMAX  DLKARM    '',
     *                ''XREMIN  XREMAX  DLKRRM  #'',/,
     *                ''#  '',F6.2,2X,F6.2,2X,F6.4,1X,
     *                ''   '',F6.3,2X,F6.4,2X,F6.4,1X,
     *                ''   '',F6.3,2X,F6.4,2X,F6.4,''  #'',/,
     *                ''#'',T80,''#'',/,
     *                ''#'',T56,''XCoMIN  XCoMAX  DLKCOU  #'',/,
     *                ''#'',T56,F6.3,2X,F6.4,2X,F6.4,''  #'',/,
     *                ''#'',T80,''#'',/,80(''#''),/)')
     *
     *                XKEMIN,XKEMAX,DLKXRM,XAEMIN,XAEMAX,DLKARM,
     *                                     XREMIN,XREMAX,DLKRRM,
     *                                     XCoMIN,XCoMAX,DLKCOU
C
      END IF           
C 
C=======================================================================
C=======================================================================
C     Choice of the argument limits and increments (KAPPA-new) 
C     random starts and for the derivatives, respectively
C=======================================================================
C=======================================================================
C
      IF (KEYWOR.EQ.'KAPA-C-LIM') THEN
          READ(*,*) V0CMIN_KAPPAR,V0CMAX_KAPPAR,XKCMIN_KAPPAR,
     *              XKCMAX_KAPPAR,A0CMIN_KAPPAR,A0CMAX_KAPPAR,
     *              XACMIN_KAPPAR,XACMAX_KAPPAR,R0CMIN_KAPPAR,
     *              R0CMAX_KAPPAR,XRCMIN_KAPPAR,XRCMAX_KAPPAR
          WRITE(LOGFIL,'(A)')KEYWOR
          WRITE(LOGFIL,'(9X,''V0CMIN  V0CMAX  XKCMIN  XKCMAX  '',
     *                      ''A0CMIN  A0CMAX  XACMIN  XACMAX  '',
     *                      ''R0CMIN  R0CMAX  XRCMIN  XRCMAX'')')
          WRITE(LOGFIL,'(9X,12(F6.2,2X))')
     *              V0CMIN_KAPPAR,V0CMAX_KAPPAR,XKCMIN_KAPPAR,
     *              XKCMAX_KAPPAR,A0CMIN_KAPPAR,A0CMAX_KAPPAR,
     *              XACMIN_KAPPAR,XACMAX_KAPPAR,R0CMIN_KAPPAR,
     *              R0CMAX_KAPPAR,XRCMIN_KAPPAR,XRCMAX_KAPPAR
      END IF
C
C=======================================================================
C
      IF (KEYWOR.EQ.'KAPASO-LIM') THEN
          READ(*,*) XL_MIN_KAPPAR,XL_MAX_KAPPAR,XKSMIN_KAPPAR,
     *              XKSMAX_KAPPAR,A0SMIN_KAPPAR,A0SMAX_KAPPAR,
     *              XASMIN_KAPPAR,XASMAX_KAPPAR,R0SMIN_KAPPAR,
     *              R0SMAX_KAPPAR,XRSMIN_KAPPAR,XRSMAX_KAPPAR

          WRITE(LOGFIL,'(A)')KEYWOR
          WRITE(LOGFIL,'(9X,''XL_MIN  XL_MAX  XKSMIN  XKSMAX  '',
     *                      ''A0SMIN  A0SMAX  XASMIN  XASMAX  '',
     *                      ''R0SMIN  R0SMAX  XRSMIN  XRSMAX'')')
          WRITE(LOGFIL,'(9X,12(F6.2,2X))')
     *              XL_MIN_KAPPAR,XL_MAX_KAPPAR,XKSMIN_KAPPAR,
     *              XKSMAX_KAPPAR,A0SMIN_KAPPAR,A0SMAX_KAPPAR,
     *              XASMIN_KAPPAR,XASMAX_KAPPAR,R0SMIN_KAPPAR,
     *              R0SMAX_KAPPAR,XRSMIN_KAPPAR,XRSMAX_KAPPAR
      END IF
C
C=======================================================================
C
      IF (KEYWOR.EQ.'STEPSKAPPA') THEN
          READ(*,*) DL_V0C_KAPPAR,DLKV0C_KAPPAR,DL_A0C_KAPPAR,
     *              DLKA0C_KAPPAR,DL_R0C_KAPPAR,DLKR0C_KAPPAR,
     *              DL_VSO_KAPPAR,DLKXLA_KAPPAR,DL_A0S_KAPPAR,
     *              DLKA0S_KAPPAR,DL_R0S_KAPPAR,DLKR0S_KAPPAR

          WRITE(LOGFIL,'(A)')KEYWOR
          WRITE(LOGFIL,'(9X,''DL_V0C  DLKV0C  DL_A0C  DLKA0C  '',
     *                      ''DL_R0C  DLKR0C  DL_VSO  DLKXLA  '',
     *                      ''DL_A0S  DLKA0S  DL_R0S  DLKR0S'')')
          WRITE(LOGFIL,'(9X,14(F6.4,2X))')
     *              DL_V0C_KAPPAR,DLKV0C_KAPPAR,DL_A0C_KAPPAR,
     *              DLKA0C_KAPPAR,DL_R0C_KAPPAR,DLKR0C_KAPPAR,
     *              DL_VSO_KAPPAR,DLKXLA_KAPPAR,DL_A0S_KAPPAR,
     *              DLKA0S_KAPPAR,DL_R0S_KAPPAR,DLKR0S_KAPPAR
      END IF
C
C=======================================================================
C     Defining the starting parameter values for the chi^2
C     minimisation as well as parameter selection
C=======================================================================
C=======================================================================
C     Protons
C=======================================================================
C=======================================================================
C 
      IF (KEYWOR.EQ.'IFTAKEPARP') THEN
C
          READ (*,*) IF_V0C_PROTON,IF_A0C_PROTON,IF_R0C_PROTON,
     *               IF_XLA_PROTON,IF_A0S_PROTON,IF_R0S_PROTON,
     *               IF_XRM_PROTON,IF_ARM_PROTON,IF_RRM_PROTON,
     *                                                  IF_COU
C
          WRITE(LOGFIL,'(A)') KEYWOR
C
          WRITE(LOGFIL,'(9X,''IF_V0C  IF_A0C  IF_R0C  IF_XLA  '',
     *                      ''IF_A0S  IF_R0S  IF_XRM  IF_ARM  '',
     *                      ''IF_RRM  IF_COU'')')
C
          WRITE(LOGFIL,'(9X,10(I6,2X))')
     *                   IF_V0C_PROTON,IF_A0C_PROTON,IF_R0C_PROTON,
     *                   IF_XLA_PROTON,IF_A0S_PROTON,IF_R0S_PROTON,
     *                   IF_XRM_PROTON,IF_ARM_PROTON,IF_RRM_PROTON,
     *                                                      IF_COU
C_______________________________________________________________________
C 
C         The values of the vector index below are associated
C         with the potential parameters. We have:
C_______________________________________________________________________
C
C         1 -> central potential depth (p)
C         2 -> central potential radius (p)
C         3 -> central potential diffuseness (p)
C
C         4 -> spin-orbit strength (p)
C         5 -> spin-orbit potential radius (p)
C         6 -> spin-orbit potential diffuseness (p)
C
C         7 -> effective mass strength (p)
C         8 -> effective mass radius (p) 
C         9 -> effective mass diffuseness (p) 
C
C         10-> Coulomb radius parameter 
C_______________________________________________________________________
C 
C         Below, the we use the following coding: 
C 
C         "V" - variable, "STR" - starting value,
C         "MIN" -> minimal value of the interval,
C         "MAX" -> maximal value of the interval,
C         "STP" -> step to calculate derivatives
C_______________________________________________________________________
C 
C         Default starting values of parameters 
C 
          VMISTR(1)=V0CENT_PROTON
          VMISTR(2)=R0CENT_PROTON
          VMISTR(3)=A0CENT_PROTON
C
          VMISTR(4)=V0SORB_PROTON
          VMISTR(5)=R0SORB_PROTON
          VMISTR(6)=A0SORB_PROTON
C
          VMISTR(7)=V0EFFM_PROTON
          VMISTR(8)=R0EFFM_PROTON
          VMISTR(9)=A0EFFM_PROTON
C
          VMISTR(10)=R0COUL
C
          V0COLD=V0CENT_PROTON
          R0COLD=R0CENT_PROTON
          A0COLD=A0CENT_PROTON
C
          V0SOLD=V0SORB_PROTON
          R0SOLD=R0SORB_PROTON
          A0SOLD=A0SORB_PROTON
C      
          V0EOLD=V0EFFM_PROTON
          R0EOLD=R0EFFM_PROTON
          A0EOLD=A0EFFM_PROTON
C      
          C0COLD=R0COUL
C_______________________________________________________________________
C
C         Minimum values (lower bounds of the intervals)
C
          VMIMIN(1)=V0CMIN
          VMIMIN(2)=R0CMIN
          VMIMIN(3)=A0CMIN
C          
          VMIMIN(4)=XL_MIN
          VMIMIN(5)=R0SMIN
          VMIMIN(6)=A0SMIN
C
          VMIMIN(7)=XEFMIN
          VMIMIN(8)=REFMIN
          VMIMIN(9)=AEFMIN
C_______________________________________________________________________
C
          VMIMIN(10)=CouMIN
C_______________________________________________________________________
C
C         Maximum values (upper bounds of the intervals)
C
          VMIMAX(1)=V0CMAX
          VMIMAX(2)=R0CMAX
          VMIMAX(3)=A0CMAX
C
          VMIMAX(4)=XL_MAX
          VMIMAX(5)=R0SMAX
          VMIMAX(6)=A0SMAX
C
          VMIMAX(7)=XEFMAX
          VMIMAX(8)=REFMAX
          VMIMAX(9)=AEFMAX
C_______________________________________________________________________
C
          VMIMAX(10)=CouMAX
C_______________________________________________________________________
C
C         Step values (for calculating the derivatives)
C
          VMISTP(1)=DL_V0C
          VMISTP(2)=DL_R0C
          VMISTP(3)=DL_A0C
C
          VMISTP(4)=DL_XLA
          VMISTP(5)=DL_R0S
          VMISTP(6)=DL_A0S
C
          VMISTP(7)=DL_XRM
          VMISTP(8)=DL_RRM
          VMISTP(9)=DL_ARM
C_______________________________________________________________________
C
          VMISTP(10)=DL_COU
C_______________________________________________________________________
C 
C         Yes/No decision: take or not to take this is the question 
C 
          IFTAKE(1)=IF_V0C_PROTON
          IFTAKE(2)=IF_R0C_PROTON
          IFTAKE(3)=IF_A0C_PROTON
C          
          IFTAKE(4)=IF_XLA_PROTON
          IFTAKE(5)=IF_R0S_PROTON
          IFTAKE(6)=IF_A0S_PROTON
C          
          IFTAKE(7)=IF_XRM_PROTON
          IFTAKE(8)=IF_RRM_PROTON
          IFTAKE(9)=IF_ARM_PROTON
C_______________________________________________________________________
C          
          IFTAKE(10)=IF_COU
C_______________________________________________________________________
C
      END IF
C
C=======================================================================
C=======================================================================
C 
      IF (KEYWOR.EQ.'IFTAKEKAPP') THEN
C
          READ (*,*) IFKV0C_PROTON,IFKA0C_PROTON,IFKR0C_PROTON,
     *	             IFKXLA_PROTON,IFKA0S_PROTON,IFKR0S_PROTON,
     *               IFKXRM_PROTON,IFKARM_PROTON,IFKRRM_PROTON,
     *                                                  IFKCOU
C
          WRITE(LOGFIL,'(A)') KEYWOR

          WRITE(LOGFIL,'(9X,''IFKV0C  IFKA0C  IFKR0C  IFKXLA  '',
     *                      ''IFKA0S  IFKR0S  IFKXRM  IFKARM  '',
     *                      ''IFKRRM  IFKCOU'')')
C
          WRITE(LOGFIL,'(9X,10(I6,2X))')
     *                   IFKV0C_PROTON,IFKA0C_PROTON,IFKR0C_PROTON,
     *	                 IFKXLA_PROTON,IFKA0S_PROTON,IFKR0S_PROTON,
     *                   IFKXRM_PROTON,IFKARM_PROTON,IFKRRM_PROTON,
     *                                                      IFKCOU
C_______________________________________________________________________
C 
C         The values of the vector index below are associated
C         with the potential parameters. We have:
C_______________________________________________________________________
C
C         11-> central potential depth (p)
C         12-> central potential radius (p)
C         13-> central potential diffuseness (p)
C
C         14-> spin-orbit strength (p)
C         15-> spin-orbit potential radius (p)
C         16-> spin-orbit potential diffuseness (p)
C
C         17-> effective mass strength (p)
C         18-> effective mass radius (p) 
C         19-> effective mass diffuseness (p) 
C
C         20-> Coulomb radius parameter 
C_______________________________________________________________________
C 
C         Below, the we use the following coding: 
C 
C         "V" - variable, "STR" - starting value,
C         "MIN" -> minimal value of the interval,
C         "MAX" -> maximal value of the interval,
C         "STP" -> step to calculate derivatives
C_______________________________________________________________________
C 
C         Default starting values of parameters 
C 
          VMISTR(11)=XK_V0C_PROTON
          VMISTR(12)=XK_R0C_PROTON
          VMISTR(13)=XK_A0C_PROTON
C
          VMISTR(14)=XK_LAM_PROTON
          VMISTR(15)=XK_RSO_PROTON
          VMISTR(16)=XK_ASO_PROTON
C
          VMISTR(17)=XK_LEF_PROTON
          VMISTR(18)=XK_REF_PROTON
          VMISTR(19)=XK_AEF_PROTON
C
          VMISTR(20)=XK_COU
C      
          OK_V0C=XK_V0C_PROTON
          OK_R0C=XK_R0C_PROTON
          OK_A0C=XK_A0C_PROTON
C      
          OK_LAM=XK_LAM_PROTON
          OK_RSO=XK_RSO_PROTON
          OK_ASO=XK_ASO_PROTON
C      
          OK_LEF=XK_LEF_PROTON
          OK_REF=XK_REF_PROTON
          OK_AEF=XK_AEF_PROTON
C      
          OK_COU=XK_COU
C_______________________________________________________________________
C
C         Minimum values (lower bounds of the intervals)
C
          VMIMIN(11)=XKCMIN
          VMIMIN(12)=XRCMIN
          VMIMIN(13)=XACMIN
C          
          VMIMIN(14)=XKSMIN
          VMIMIN(15)=XRSMIN
          VMIMIN(16)=XASMIN
C
          VMIMIN(17)=XKEMIN
          VMIMIN(18)=XREMIN
          VMIMIN(19)=XAEMIN
C_______________________________________________________________________
C
          VMIMIN(20)=XCoMIN
C_______________________________________________________________________
C
C         Maximum values (upper bounds of the intervals)
C
          VMIMAX(11)=XKCMAX
          VMIMAX(12)=XRCMAX
          VMIMAX(13)=XACMAX
C
          VMIMAX(14)=XKSMAX
          VMIMAX(15)=XRSMAX
          VMIMAX(16)=XASMAX
C
          VMIMAX(17)=XKEMAX
          VMIMAX(18)=XREMAX
          VMIMAX(19)=XAEMAX
C_______________________________________________________________________
C
          VMIMAX(20)=XCoMAX
C_______________________________________________________________________
C
C         Step values (for calculating the derivatives - KAPPA)
C
          VMISTP(11)=DLKV0C
          VMISTP(12)=DLKR0C
          VMISTP(13)=DLKA0C
C
          VMISTP(14)=DLKXLA
          VMISTP(15)=DLKR0S
          VMISTP(16)=DLKA0S
C
          VMISTP(17)=DLKXRM
          VMISTP(18)=DLKRRM
          VMISTP(19)=DLKARM
C_______________________________________________________________________
C
          VMISTP(20)=DLKCOU
C_______________________________________________________________________
C 
C         Yes/No decision: take or not to take this KAPPA 
C 
          IFTAKE(11)=IFKV0C_PROTON
          IFTAKE(12)=IFKR0C_PROTON
          IFTAKE(13)=IFKA0C_PROTON
C          
          IFTAKE(14)=IFKXLA_PROTON
          IFTAKE(15)=IFKR0S_PROTON
          IFTAKE(16)=IFKA0S_PROTON
C          
          IFTAKE(17)=IFKXRM_PROTON
          IFTAKE(18)=IFKRRM_PROTON
          IFTAKE(19)=IFKARM_PROTON
C_______________________________________________________________________
C          
          IFTAKE(20)=IFKCOU
C_______________________________________________________________________
C
          DO I=1,NDPARS
             DPARAM(I)=VMISTP(I)
          END DO
C_______________________________________________________________________
C
      END IF
C
C=======================================================================
C=======================================================================
C     Neutrons
C=======================================================================
C=======================================================================
C
      IF (KEYWOR.EQ.'IFTAKEPARN') THEN
C
          READ (*,*) IF_V0C_NEUTRS,IF_A0C_NEUTRS,IF_R0C_NEUTRS,
     *               IF_XLA_NEUTRS,IF_A0S_NEUTRS,IF_R0S_NEUTRS,
     *               IF_XRM_NEUTRS,IF_ARM_NEUTRS,IF_RRM_NEUTRS
C
          WRITE(LOGFIL,'(A)') KEYWOR
C
          WRITE(LOGFIL,'(9X,''IF_V0C  IF_A0C  IF_R0C  IF_XLA  '',
     *                      ''IF_A0S  IF_R0S  IF_XRM  IF_ARM  '',
     *                      ''IF_RRM '')')
C
          WRITE(LOGFIL,'(9X,9(I6,2X))')
     *                   IF_V0C_NEUTRS,IF_A0C_NEUTRS,IF_R0C_NEUTRS,
     *                   IF_XLA_NEUTRS,IF_A0S_NEUTRS,IF_R0S_NEUTRS,
     *                   IF_XRM_NEUTRS,IF_ARM_NEUTRS,IF_RRM_NEUTRS
C_______________________________________________________________________
C 
C         The values of the vector index below are associated
C         with the potential parameters. We have:
C_______________________________________________________________________
C
C         21 -> central potential depth (n)
C         22 -> central potential radius (n)
C         23 -> central potential diffuseness (n)
C
C         24 -> spin-orbit strength (n)
C         25 -> spin-orbit potential radius (n)
C         26 -> spin-orbit potential diffuseness (n)
C
C         27 -> effective mass strength (n)
C         28 -> effective mass radius (n) 
C         29 -> effective mass diffuseness (n) 
C_______________________________________________________________________
C 
C         Below, the we use  the followng coding: 
C 
C         "V" - variable, "STR" - starting value,
C         "MIN" -> minimal value of the interval,
C         "MAX" -> maximal value of the interval,
C         "STP" -> step to calculate derivatives
C_______________________________________________________________________
C 
C         Default starting values of parameters 
C 
          VMISTR(21)=V0CENT_NEUTRS
          VMISTR(22)=R0CENT_NEUTRS
          VMISTR(23)=A0CENT_NEUTRS
C
          VMISTR(24)=V0SORB_NEUTRS
          VMISTR(25)=R0SORB_NEUTRS
          VMISTR(26)=A0SORB_NEUTRS
C
          VMISTR(27)=V0EFFM_NEUTRS
          VMISTR(28)=R0EFFM_NEUTRS
          VMISTR(29)=A0EFFM_NEUTRS
C
          V0COLD=V0CENT_PROTON
          R0COLD=R0CENT_PROTON
          A0COLD=A0CENT_PROTON
C
          V0SOLD=V0SORB_PROTON
          R0SOLD=R0SORB_PROTON
          A0SOLD=A0SORB_PROTON
C      
          V0EOLD=V0EFFM_PROTON
          R0EOLD=R0EFFM_PROTON
          A0EOLD=A0EFFM_PROTON
C_______________________________________________________________________
C
C         Minimum values (lower bounds of the intervals)
C
          VMIMIN(21)=V0CMIN
          VMIMIN(22)=R0CMIN
          VMIMIN(23)=A0CMIN
C          
          VMIMIN(24)=XL_MIN
          VMIMIN(25)=R0SMIN
          VMIMIN(26)=A0SMIN
C
          VMIMIN(27)=XEFMIN
          VMIMIN(28)=REFMIN
          VMIMIN(29)=AEFMIN
C_______________________________________________________________________
C
C         Maximum values (upper bounds of the intervals)
C
          VMIMAX(21)=V0CMAX
          VMIMAX(22)=R0CMAX
          VMIMAX(23)=A0CMAX
C
          VMIMAX(24)=XL_MAX
          VMIMAX(25)=R0SMAX
          VMIMAX(26)=A0SMAX
C
          VMIMAX(27)=XEFMAX
          VMIMAX(28)=REFMAX
          VMIMAX(29)=AEFMAX
C_______________________________________________________________________
C
C         Step values (for calculating the derivatives)
C
          VMISTP(21)=DL_V0C
          VMISTP(22)=DL_R0C
          VMISTP(23)=DL_A0C
C
          VMISTP(24)=DL_XLA
          VMISTP(25)=DL_R0S
          VMISTP(26)=DL_A0S
C
          VMISTP(27)=DL_XRM
          VMISTP(28)=DL_RRM
          VMISTP(29)=DL_ARM
C_______________________________________________________________________
C 
C         Yes/No decision: take or not to take this is the question 
C 
          IFTAKE(21)=IF_V0C_NEUTRS
          IFTAKE(22)=IF_R0C_NEUTRS
          IFTAKE(23)=IF_A0C_NEUTRS
C          
          IFTAKE(24)=IF_XLA_NEUTRS
          IFTAKE(25)=IF_R0S_NEUTRS
          IFTAKE(26)=IF_A0S_NEUTRS
C          
          IFTAKE(27)=IF_XRM_NEUTRS
          IFTAKE(28)=IF_RRM_NEUTRS
          IFTAKE(29)=IF_ARM_NEUTRS
C_______________________________________________________________________
C
      END IF
C
C=======================================================================
C=======================================================================
C 
      IF (KEYWOR.EQ.'IFTAKEKAPN') THEN
C
          READ (*,*) IFKV0C_NEUTRS,IFKA0C_NEUTRS,IFKR0C_NEUTRS,
     *	             IFKXLA_NEUTRS,IFKA0S_NEUTRS,IFKR0S_NEUTRS,
     *               IFKXRM_NEUTRS,IFKARM_NEUTRS,IFKRRM_NEUTRS
C
          WRITE(LOGFIL,'(A)') KEYWOR
C
          WRITE(LOGFIL,'(9X,''IFKV0C  IFKA0C  IFKR0C  IFKXLA  '',
     *                      ''IFKA0S  IFKR0S  IFKXRM  IFKARM  '',
     *                      ''IFKRRM  '')')
C
          WRITE(LOGFIL,'(9X,9(I6,2X))')
     *                   IFKV0C_NEUTRS,IFKA0C_NEUTRS,IFKR0C_NEUTRS,
     *	                 IFKXLA_NEUTRS,IFKA0S_NEUTRS,IFKR0S_NEUTRS,
     *                   IFKXRM_NEUTRS,IFKARM_NEUTRS,IFKRRM_NEUTRS
C_______________________________________________________________________
C 
C         The values of the vector index below are associated
C         with the potential parameters. We have:
C_______________________________________________________________________
C
C         30-> central potential depth (p)
C         31-> central potential radius (p)
C         32-> central potential diffuseness (p)
C
C         33-> spin-orbit strength (p)
C         34-> spin-orbit potential radius (p)
C         35-> spin-orbit potential diffuseness (p)
C
C         36-> effective mass strength (p)
C         37-> effective mass radius (p) 
C         38-> effective mass diffuseness (p) 
C_______________________________________________________________________
C 
C         Below, the we use the following coding: 
C 
C         "V" - variable, "STR" - starting value,
C         "MIN" -> minimal value of the interval,
C         "MAX" -> maximal value of the interval,
C         "STP" -> step to calculate derivatives
C_______________________________________________________________________
C 
C         Default starting values of parameters 
C 
          VMISTR(30)=XK_V0C_PROTON
          VMISTR(31)=XK_R0C_PROTON
          VMISTR(32)=XK_A0C_PROTON
C
          VMISTR(33)=XK_LAM_PROTON
          VMISTR(34)=XK_RSO_PROTON
          VMISTR(35)=XK_ASO_PROTON
C
          VMISTR(36)=XK_LEF_PROTON
          VMISTR(37)=XK_REF_PROTON
          VMISTR(37)=XK_AEF_PROTON
C      
          OK_V0C=XK_V0C_PROTON
          OK_R0C=XK_R0C_PROTON
          OK_A0C=XK_A0C_PROTON
C      
          OK_LAM=XK_LAM_PROTON
          OK_RSO=XK_RSO_PROTON
          OK_ASO=XK_ASO_PROTON
C      
          OK_LEF=XK_LEF_PROTON
          OK_REF=XK_REF_PROTON
          OK_AEF=XK_AEF_PROTON
C      
          OK_COU=XK_COU
C_______________________________________________________________________
C
C         Minimum values (lower bounds of the intervals)
C
          VMIMIN(30)=XKCMIN
          VMIMIN(31)=XRCMIN
          VMIMIN(32)=XACMIN
C          
          VMIMIN(33)=XKSMIN
          VMIMIN(34)=XRSMIN
          VMIMIN(35)=XASMIN
C
          VMIMIN(36)=XKEMIN
          VMIMIN(37)=XREMIN
          VMIMIN(38)=XAEMIN
C_______________________________________________________________________
C
C         Maximum values (upper bounds of the intervals)
C
          VMIMAX(30)=XKCMAX
          VMIMAX(31)=XRCMAX
          VMIMAX(32)=XACMAX
C
          VMIMAX(33)=XKSMAX
          VMIMAX(34)=XRSMAX
          VMIMAX(35)=XASMAX
C
          VMIMAX(36)=XKEMAX
          VMIMAX(37)=XREMAX
          VMIMAX(38)=XAEMAX
C_______________________________________________________________________
C
C         Step values (for calculating the derivatives - KAPPA)
C
          VMISTP(30)=DLKV0C
          VMISTP(31)=DLKR0C
          VMISTP(32)=DLKA0C
C
          VMISTP(33)=DLKXLA
          VMISTP(34)=DLKR0S
          VMISTP(35)=DLKA0S
C
          VMISTP(36)=DLKXRM
          VMISTP(37)=DLKRRM
          VMISTP(38)=DLKARM
C_______________________________________________________________________
C 
C         Yes/No decision: take or not to take this KAPPA 
C 
          IFTAKE(30)=IFKV0C_NEUTRS
          IFTAKE(31)=IFKR0C_NEUTRS
          IFTAKE(32)=IFKA0C_NEUTRS
C          
          IFTAKE(33)=IFKXLA_NEUTRS
          IFTAKE(34)=IFKR0S_NEUTRS
          IFTAKE(35)=IFKA0S_NEUTRS
C          
          IFTAKE(36)=IFKXRM_NEUTRS
          IFTAKE(37)=IFKRRM_NEUTRS
          IFTAKE(38)=IFKARM_NEUTRS
C_______________________________________________________________________
C
          DO I=1,NDPARS
             DPARAM(I)=VMISTP(I)
          END DO
C_______________________________________________________________________
C
      END IF
C
C=======================================================================
C=======================================================================
C     Density-dependent spin-orbit Hamiltonian
C=======================================================================
C=======================================================================
C 
      IF (KEYWOR.EQ.'IFTAKEDENS') THEN
C
          READ (*,*) IF_XPP,IF_XPN,IF_XNP,IF_XNN
C
          IF (IFDENS.EQ.1) THEN
C
              WRITE(LOGFIL,'(A)')KEYWOR
C
              WRITE(LOGFIL,'(9X,''IF_XPP  IF_XPN  IF_XNP  IF_XNN'')')
              WRITE(LOGFIL,'(9X,4(I6,2X))')IF_XPP,IF_XPN,IF_XNP,IF_XNN
C
          END IF
C_______________________________________________________________________
C 
C         The values of the vector index below are associated
C         with the potential parameters. We have:
C_______________________________________________________________________
C
C         39 -> lambda  proton-proton
C         40 -> lambda  proton-neutron
C         41 -> lambda neutron-proton
C         42 -> lambda neutron-neutron
C_______________________________________________________________________
C 
C         Below, the we use the following coding: 
C 
C         "V" - variable, "STR" - starting value,
C         "MIN" -> minimal value of the interval,
C         "MAX" -> maximal value of the interval,
C         "STP" -> step to calculate derivatives
C_______________________________________________________________________
C 
C         Default starting values of parameters 
C 
          VMISTR(39)=ALAMPP
          VMISTR(40)=ALAMPN
          VMISTR(41)=ALAMNP
          VMISTR(42)=ALAMNN
C_______________________________________________________________________
C
C         Minimum values (lower bounds of the intervals)
C
          VMIMIN(39)=XPPMIN
          VMIMIN(40)=XPNMIN
          VMIMIN(41)=XNPMIN
          VMIMIN(42)=XNNMIN
C_______________________________________________________________________
C
C         Maximum values (upper bounds of the intervals)
C
          VMIMAX(39)=XPPMAX
          VMIMAX(40)=XPNMAX
          VMIMAX(41)=XNPMAX
          VMIMAX(42)=XNNMAX
C_______________________________________________________________________
C
C         Step values (for calculating the derivatives)
C
          VMISTP(39)=DL_XPP
          VMISTP(40)=DL_XPN
          VMISTP(41)=DL_XNP
          VMISTP(42)=DL_XNN
C_______________________________________________________________________
C 
C         Yes/No decision: take or not to take this is the question 
C 
          IFTAKE(39)=IF_XPP
          IFTAKE(40)=IF_XPN
          IFTAKE(41)=IF_XNP
          IFTAKE(42)=IF_XNN
C_______________________________________________________________________
C
      END IF
C
C=======================================================================
C=======================================================================
C     Tensor part for SO-Potential
C=======================================================================
C=======================================================================
C 
      IF (KEYWOR.EQ.'IFTAKETENS') THEN
C
          READ (*,*) IF_YPP,IF_YPN,IF_YNP,IF_YNN
C
          IF (IFTENS.EQ.1 .AND. ISORBT.EQ.1) THEN
C
              WRITE(LOGFIL,'(A)') KEYWOR
C
              WRITE(LOGFIL,'(9X,''IF_YPP  IF_YPN  IF_YNP  IF_YNN'')')
              WRITE(LOGFIL,'(9X,4(I6,2X))')IF_YPP,IF_YPN,IF_YNP,IF_YNN
C
          END IF
C_______________________________________________________________________
C 
C         The values of the vector index below are associated
C         with the potential parameters. We have:
C_______________________________________________________________________
C
C         43 -> tensor lambda  proton-proton
C         44 -> tensor lambda  proton-neutron
C         45 -> tensor lambda neutron-proton
C         46 -> tensor lambda neutron-neutron
C_______________________________________________________________________
C 
C         Below, the we use the following coding: 
C 
C         "V" - variable, "STR" - starting value,
C         "MIN" -> minimal value of the interval,
C         "MAX" -> maximal value of the interval,
C         "STP" -> step to calculate derivatives
C_______________________________________________________________________
C 
C         Default starting values of parameters 
C 
          VMISTR(43)=TLAMPP
          VMISTR(44)=TLAMPN
          VMISTR(45)=TLAMNP
          VMISTR(46)=TLAMNN
C_______________________________________________________________________
C
C         Minimum values (lower bounds of the intervals)
C
          VMIMIN(43)=YPPMIN
          VMIMIN(44)=YPNMIN
          VMIMIN(45)=YNPMIN
          VMIMIN(46)=YNNMIN
C_______________________________________________________________________
C
C         Maximum values (upper bounds of the intervals)
C
          VMIMAX(43)=YPPMAX
          VMIMAX(44)=YPNMAX
          VMIMAX(45)=YNPMAX
          VMIMAX(46)=YNNMAX
C_______________________________________________________________________
C
C         Step values (for calculating the derivatives)
C
          VMISTP(43)=DL_YPP
          VMISTP(44)=DL_YPN
          VMISTP(45)=DL_YNP
          VMISTP(46)=DL_YNN
C_______________________________________________________________________
C 
C         Yes/No decision: take or not to take this is the question 
C 
          IFTAKE(43)=IF_YPP
          IFTAKE(44)=IF_YPN
          IFTAKE(45)=IF_YNP
          IFTAKE(46)=IF_YNN
C_______________________________________________________________________
C
          DO I=1,NDPARS
             DPARAM(I)=VMISTP(I)
          END DO
C_______________________________________________________________________
C
      END IF
C
C=======================================================================
C=======================================================================
C     Tensor part for Central-Potential
C=======================================================================
C=======================================================================
C 
      IF (KEYWOR.EQ.'IFTAKECNTN') THEN
C
          READ (*,*) IF_CPP,IF_CPN,IF_CNP,IF_CNN
C
          IF (IFTENS.EQ.1 .AND. ICENTT.EQ.1) THEN
C
              WRITE(LOGFIL,'(A)') KEYWOR
C
              WRITE(LOGFIL,'(9X,''IF_CPP  IF_CPN  IF_CNP  IF_CNN'')')
              WRITE(LOGFIL,'(9X,4(I6,2X))')IF_CPP,IF_CPN,IF_CNP,IF_CNN
C
          END IF
C_______________________________________________________________________
C 
C         The values of the vector index below are associated
C         with the potential parameters. We have:
C_______________________________________________________________________
C
C         47 -> central tensor lambda  proton-proton
C         48 -> central tensor lambda  proton-neutron
C         49 -> central tensor lambda neutron-proton
C         50 -> central tensor lambda neutron-neutron
C_______________________________________________________________________
C 
C         Below, the we use the following coding: 
C 
C         "V" - variable, "STR" - starting value,
C         "MIN" -> minimal value of the interval,
C         "MAX" -> maximal value of the interval,
C         "STP" -> step to calculate derivatives
C_______________________________________________________________________
C 
C         Default starting values of parameters 
C 
          VMISTR(47)=CLAMPP
          VMISTR(48)=CLAMPN
          VMISTR(49)=CLAMNP
          VMISTR(50)=CLAMNN
C_______________________________________________________________________
C
C         Minimum values (lower bounds of the intervals)
C
          VMIMIN(47)=CPPMIN
          VMIMIN(48)=CPNMIN
          VMIMIN(49)=CNPMIN
          VMIMIN(50)=CNNMIN
C_______________________________________________________________________
C
C         Maximum values (upper bounds of the intervals)
C
          VMIMAX(47)=CPPMAX
          VMIMAX(48)=CPNMAX
          VMIMAX(49)=CNPMAX
          VMIMAX(50)=CNNMAX
C_______________________________________________________________________
C
C         Step values (for calculating the derivatives)
C
          VMISTP(47)=DL_CPP
          VMISTP(48)=DL_CPN
          VMISTP(49)=DL_CNP
          VMISTP(50)=DL_CNN
C_______________________________________________________________________
C 
C         Yes/No decision: take or not to take this is the question 
C 
          IFTAKE(47)=IF_CPP
          IFTAKE(48)=IF_CPN
          IFTAKE(49)=IF_CNP
          IFTAKE(50)=IF_CNN
C_______________________________________________________________________
C
          DO I=1,NDPARS
             DPARAM(I)=VMISTP(I)
          END DO
C_______________________________________________________________________
C
      END IF
C
C=======================================================================
C=======================================================================
C     Kappa parametrization: Central Potential
C=======================================================================
C=======================================================================
C
      IF (KEYWOR.EQ.'IFTAKEKAPC') THEN
C
          READ(*,*) IFTV0C,IFTKVC,IFTR0C,IFTKRC,IFTA0C,IFTKAC
C
          WRITE(LOGFIL,'(A)') KEYWOR
          WRITE(LOGFIL,'(9X,''IFTV0C  IFTKVC  IFTR0C  IFTKRC  '',
     *                      ''IFTA0C  IFTKAC'')')
          WRITE(LOGFIL,'(9X,6(I6,2X))')IFTV0C,IFTKVC,IFTR0C,
     *                                 IFTKRC,IFTA0C,IFTKAC
C_______________________________________________________________________
C 
C         The values of the vector index below are associated
C         with the potential parameters. We have:
C_______________________________________________________________________
C
C         51 -> V_0
C         52 -> kappa for V_0
C         53 -> r_0
C         54 -> kappa for r_0
C         55 -> a_0
C         56 -> kappa for a_0
C_______________________________________________________________________
C 
C         Below, the we use the following coding: 
C 
C         "V" - variable, "STR" - starting value,
C         "MIN" -> minimal value of the interval,
C         "MAX" -> maximal value of the interval,
C         "STP" -> step to calculate derivatives
C_______________________________________________________________________
C 
C         Default starting values of parameters 
C 
          VMISTR(51)=V0CENT_KAPPAR
          VMISTR(52)=XK_V0C_KAPPAR
          VMISTR(53)=R0CENT_KAPPAR
          VMISTR(54)=XK_R0C_KAPPAR
          VMISTR(55)=A0CENT_KAPPAR
          VMISTR(56)=XK_A0C_KAPPAR
C_______________________________________________________________________
C
C         Minimum values (lower bounds of the intervals)
C
          VMIMIN(51)=V0CMIN_KAPPAR
          VMIMIN(52)=XKCMIN_KAPPAR
          VMIMIN(53)=R0CMIN_KAPPAR
          VMIMIN(54)=XRCMIN_KAPPAR
          VMIMIN(55)=A0CMIN_KAPPAR
          VMIMIN(56)=XACMIN_KAPPAR
C_______________________________________________________________________
C
C         Maximum values (upper bounds of the intervals)
C
          VMIMAX(51)=V0CMAX_KAPPAR
          VMIMAX(52)=XKCMAX_KAPPAR
          VMIMAX(53)=R0CMAX_KAPPAR
          VMIMAX(54)=XRCMAX_KAPPAR
          VMIMAX(55)=A0CMAX_KAPPAR
          VMIMAX(56)=XACMAX_KAPPAR
C_______________________________________________________________________
C
C         Step values (for calculating the derivatives)
C
          VMISTP(51)=DL_V0C_KAPPAR
          VMISTP(52)=DLKV0C_KAPPAR
          VMISTP(53)=DL_R0C_KAPPAR
          VMISTP(54)=DLKR0C_KAPPAR
          VMISTP(55)=DL_A0C_KAPPAR
          VMISTP(56)=DLKA0C_KAPPAR
C_______________________________________________________________________
C 
C         Yes/No decision: take or not to take this is the question 
C 
          IFTAKE(51)=IFTV0C
          IFTAKE(52)=IFTKVC
          IFTAKE(53)=IFTR0C
          IFTAKE(54)=IFTKRC
          IFTAKE(55)=IFTA0C
          IFTAKE(56)=IFTKAC
C_______________________________________________________________________
C
          DO I=1,NDPARS
             DPARAM(I)=VMISTP(I)
          END DO
C_______________________________________________________________________
C
      END IF
C
C=======================================================================
C=======================================================================
C     Kappa parametrization: Pure SO WS Potential
C=======================================================================
C=======================================================================
C
      IF (KEYWOR.EQ.'IFTAKEKAPS') THEN
C
          READ(*,*) IFTV0S,IFTKVS,IFTR0S,IFTKRS,IFTA0S,IFTKAS
C
          WRITE(LOGFIL,'(A)') KEYWOR
          WRITE(LOGFIL,'(9X,''IFTV0S  IFTKVS  IFTR0S  IFTKRS  '',
     *                      ''IFTA0S  IFTKAS'')')
          WRITE(LOGFIL,'(9X,6(I6,2X))')IFTV0S,IFTKVS,IFTR0S,
     *                                 IFTKRS,IFTA0S,IFTKAS
C_______________________________________________________________________
C 
C         The values of the vector index below are associated
C         with the potential parameters. We have:
C_______________________________________________________________________
C
C         57 -> lambda_0
C         58 -> kappa for lambda_0
C         59 -> r_0
C         60 -> kappa for r_0
C         61 -> a_0
C         62 -> kappa for a_0
C_______________________________________________________________________
C 
C         Below, the we use the following coding: 
C 
C         "V" - variable, "STR" - starting value,
C         "MIN" -> minimal value of the interval,
C         "MAX" -> maximal value of the interval,
C         "STP" -> step to calculate derivatives
C_______________________________________________________________________
C 
C         Default starting values of parameters 
C 
          VMISTR(57)=V0SORB_KAPPAR
          VMISTR(58)=XK_V0S_KAPPAR
          VMISTR(59)=R0SORB_KAPPAR
          VMISTR(60)=XK_R0S_KAPPAR
          VMISTR(61)=A0SORB_KAPPAR
          VMISTR(62)=XK_A0S_KAPPAR
C_______________________________________________________________________
C
C         Minimum values (lower bounds of the intervals)
C
          VMIMIN(57)=XL_MIN_KAPPAR
          VMIMIN(58)=XKSMIN_KAPPAR
          VMIMIN(59)=R0SMIN_KAPPAR
          VMIMIN(60)=XRSMIN_KAPPAR
          VMIMIN(61)=A0SMIN_KAPPAR
          VMIMIN(62)=XASMIN_KAPPAR
C_______________________________________________________________________
C
C         Maximum values (upper bounds of the intervals)
C
          VMIMAX(57)=XL_MAX_KAPPAR
          VMIMAX(58)=XKSMAX_KAPPAR
          VMIMAX(59)=R0SMAX_KAPPAR
          VMIMAX(60)=XRSMAX_KAPPAR
          VMIMAX(61)=A0SMAX_KAPPAR
          VMIMAX(62)=XASMAX_KAPPAR
C_______________________________________________________________________
C
C         Step values (for calculating the derivatives)
C
          VMISTP(57)=DL_VSO_KAPPAR
          VMISTP(58)=DLKXLA_KAPPAR
          VMISTP(59)=DL_R0S_KAPPAR
          VMISTP(60)=DLKR0S_KAPPAR
          VMISTP(61)=DL_A0S_KAPPAR
          VMISTP(62)=DLKA0S_KAPPAR
C_______________________________________________________________________
C 
C         Yes/No decision: take or not to take this is the question 
C 
          IFTAKE(57)=IFTV0S
          IFTAKE(58)=IFTKVS
          IFTAKE(59)=IFTR0S
          IFTAKE(60)=IFTKRS
          IFTAKE(61)=IFTA0S
          IFTAKE(62)=IFTKAS
C_______________________________________________________________________
C
          DO I=1,NDPARS
             DPARAM(I)=VMISTP(I)
          END DO
C_______________________________________________________________________
C
      END IF
C 
C=======================================================================
C=======================================================================
C     Central potential parameters  - ranges for mesh option
C=======================================================================
C=======================================================================
C 
      IF (KEYWOR.EQ.'CENTMESHLM') THEN
C
          READ (N_READ,*) V0CMIN_IFMESH,V0CMAX_IFMESH,NOPV0C,
     *                    A0CMIN_IFMESH,A0CMAX_IFMESH,NOPA0C,
     *                    R0CMIN_IFMESH,R0CMAX_IFMESH,NOPR0C,
     *                    CouMIN_IFMESH,CouMAX_IFMESH,NOPCou
C
          IF (IFMESH.EQ.1) THEN
C
              WRITE(LOGFIL,'(A)') KEYWOR
C
              WRITE(LOGFIL,'(9X,''V0CMIN  V0CMAX  NOPV0C  '',
     *                          ''A0CMIN  A0CMAX  NOPA0C  '',
     *                          ''R0CMIN  R0CMAX  NOPR0C  '',
     *                          ''CouMIN  CouMAX  NPOCou'',/,9X,
     *                            F6.1,2X,F6.2,2X,I6,2X,
     *                            F6.2,2X,F6.2,2X,I6,2X,
     *                            F6.2,2X,F6.2,2X,I6,2X,
     *                            F6.2,2X,F6.2,2X,I6,2X)')
     *
     *                            V0CMIN_IFMESH,V0CMAX_IFMESH,NOPV0C,
     *                            A0CMIN_IFMESH,A0CMAX_IFMESH,NOPA0C,
     *                            R0CMIN_IFMESH,R0CMAX_IFMESH,NOPR0C,
     *                            CouMIN_IFMESH,CouMAX_IFMESH,NOPCou
          END IF
C
      END IF          
C 
C=======================================================================
C=======================================================================
C     Spin-orbit potential parameter  ranges for mesh option
C=======================================================================
C=======================================================================
C 
      IF (KEYWOR.EQ.'SORBMESHLM') THEN
C
          READ (N_READ,*) XL_MIN_IFMESH,XL_MAX_IFMESH,NOP_XL,
     *                    A0SMIN_IFMESH,A0SMAX_IFMESH,NOPA0S,
     *                    R0SMIN_IFMESH,R0SMAX_IFMESH,NOPR0S
C
          IF (IFDENS.EQ.0 .AND. IFMESH.EQ.1) THEN
C
              WRITE(LOGFIL,'(A)') KEYWOR
C
              WRITE(LOGFIL,'(9X,''XL_MIN  XL_MAX  NOP_XL  '',
     *                          ''A0SMIN  A0SMAX  NOPA0S  '',
     *                          ''R0SMIN  R0SMAX  NOPR0S'',/,9X,
     *                            F6.2,2X,F6.2,2X,I6,2X,
     *                            F6.2,2X,F6.2,2X,I6,2X,
     *                            F6.2,2X,F6.2,2X,I6,2X)')
     *
     *                            XL_MIN_IFMESH,XL_MAX_IFMESH,NOP_XL,
     *                            A0SMIN_IFMESH,A0SMAX_IFMESH,NOPA0S,
     *                            R0SMIN_IFMESH,R0SMAX_IFMESH,NOPR0S
          END IF
C
      END IF          
C 
C=======================================================================
C=======================================================================
C     Effective-mass term parameters, ranges for mesh option
C=======================================================================
C=======================================================================
C 
      IF (KEYWOR.EQ.'EFFMMESHLM') THEN
C
          READ (N_READ,*) XEFMIN_IFMESH,XEFMAX_IFMESH,NOPXEF,
     *                    AEFMIN_IFMESH,AEFMAX_IFMESH,NOPAEF,
     *                    REFMIN_IFMESH,REFMAX_IFMESH,NOPREF
C
          IF (IFMESH.EQ.1) THEN
C
              WRITE(LOGFIL,'(A)') KEYWOR
C
              WRITE(LOGFIL,'(9X,''XEFMIN  XEFMAX  NOPXEF  '',
     *                          ''AEFMIN  AEFMAX  NOPAEF  '',
     *                          ''REFMIN  REFMAX  NOPREF'',/,9X,
     *                            F6.2,2X,F6.2,2X,I6,2X,
     *                            F6.2,2X,F6.2,2X,I6,2X,
     *                            F6.2,2X,F6.2,2X,I6,2X)')
     *
     *                            XEFMIN_IFMESH,XEFMAX_IFMESH,NOPXEF,
     *                            AEFMIN_IFMESH,AEFMAX_IFMESH,NOPAEF,
     *                            REFMIN_IFMESH,REFMAX_IFMESH,NOPREF
          END IF
C
      END IF          
C 
C=======================================================================
C=======================================================================
C     Density dependent S-O pot. parameters - ranges for mesh option
C=======================================================================
C=======================================================================
C 
      IF (KEYWOR.EQ.'DENSMESHLM') THEN
C
          READ (N_READ,*) XPPMIN_IFMESH,XPPMAX_IFMESH,NOPXPP,
     *                    XPNMIN_IFMESH,XPNMAX_IFMESH,NOPXPN,
     *                    XNPMIN_IFMESH,XNPMAX_IFMESH,NOPXNP,
     *                    XNNMIN_IFMESH,XNNMAX_IFMESH,NOPXNN
C
          IF (IFDENS.EQ.1 .AND. IFMESH.EQ.1) THEN
C
              WRITE(LOGFIL,'(A)') KEYWOR
C
              WRITE(LOGFIL,'(9X,''XPPMIN  XPPMAX  NOPXPP  '',
     *                          ''XPNMIN  XPNMAX  NOPXPN  '',
     *                          ''XNPMIN  XNPMAX  NOPXNP  '',
     *                          ''XNNMIN  XNNMAX  NOPXNN'',/,9X,
     *                            F6.2,2X,F6.2,2X,I6,2X,
     *                            F6.2,2X,F6.2,2X,I6,2X,
     *                            F6.2,2X,F6.2,2X,I6,2X,
     *                            F6.2,2X,F6.2,2X,I6,2X)')
     *
     *                            XPPMIN_IFMESH,XPPMAX_IFMESH,NOPXPP,
     *                            XPNMIN_IFMESH,XPNMAX_IFMESH,NOPXPN,
     *                            XNPMIN_IFMESH,XNPMAX_IFMESH,NOPXNP,
     *                            XNNMIN_IFMESH,XNNMAX_IFMESH,NOPXNN
          END IF
C
      END IF          
C 
C=======================================================================
C=======================================================================
C     SO-Tensor parameters - ranges for mesh option
C=======================================================================
C=======================================================================
C 
      IF (KEYWOR.EQ.'TENSMESHLM') THEN
C
          READ (N_READ,*) YPPMIN_IFMESH,YPPMAX_IFMESH,NOPYPP,
     *                    YPNMIN_IFMESH,YPNMAX_IFMESH,NOPYPN,
     *                    YNPMIN_IFMESH,YNPMAX_IFMESH,NOPYNP,
     *                    YNNMIN_IFMESH,YNNMAX_IFMESH,NOPYNN
C
          IF (IFTENS.EQ.1 .AND. ISORBT.EQ.1 .AND. IFMESH.EQ.1) THEN
C
              WRITE(LOGFIL,'(A)') KEYWOR
C
              WRITE(LOGFIL,'(9X,''YPPMIN  YPPMAX  NOPYPP  '',
     *                          ''YPNMIN  YPNMAX  NOPYPN  '',
     *                          ''YNPMIN  YNPMAX  NOPYNP  '',
     *                          ''YNNMIN  YNNMAX  NOPYNN'',/,9X,
     *                            F6.2,2X,F6.2,2X,I6,2X,
     *                            F6.2,2X,F6.2,2X,I6,2X,
     *                            F6.2,2X,F6.2,2X,I6,2X,
     *                            F6.2,2X,F6.2,2X,I6,2X)')
     *
     *                            YPPMIN_IFMESH,YPPMAX_IFMESH,NOPYPP,
     *                            YPNMIN_IFMESH,YPNMAX_IFMESH,NOPYPN,
     *                            YNPMIN_IFMESH,YNPMAX_IFMESH,NOPYNP,
     *                            YNNMIN_IFMESH,YNNMAX_IFMESH,NOPYNN
          END IF
C
      END IF 
C 
C=======================================================================
C=======================================================================
C     Central-Tensor parameters - ranges for mesh option
C=======================================================================
C=======================================================================
C 
      IF (KEYWOR.EQ.'CNTNMESHLM') THEN
C
          READ (N_READ,*) CPPMIN_IFMESH,CPPMAX_IFMESH,NOPCPP,
     *                    CPNMIN_IFMESH,CPNMAX_IFMESH,NOPCPN,
     *                    CNPMIN_IFMESH,CNPMAX_IFMESH,NOPCNP,
     *                    CNNMIN_IFMESH,CNNMAX_IFMESH,NOPCNN
C
          IF (IFTENS.EQ.1 .AND. ICENTT.EQ.1 .AND. IFMESH.EQ.1) THEN
C
              WRITE(LOGFIL,'(A)') KEYWOR
C
              WRITE(LOGFIL,'(9X,''CPPMIN  CPPMAX  NOPCPP  '',
     *                          ''CPNMIN  CPNMAX  NOPCPN  '',
     *                          ''CNPMIN  CNPMAX  NOPCNP  '',
     *                          ''CNNMIN  CNNMAX  NOPCNN'',/,9X,
     *                            F6.2,2X,F6.2,2X,I6,2X,
     *                            F6.2,2X,F6.2,2X,I6,2X,
     *                            F6.2,2X,F6.2,2X,I6,2X,
     *                            F6.2,2X,F6.2,2X,I6,2X)')
     *
     *                            CPPMIN_IFMESH,CPPMAX_IFMESH,NOPCPP,
     *                            CPNMIN_IFMESH,CPNMAX_IFMESH,NOPCPN,
     *                            CNPMIN_IFMESH,CNPMAX_IFMESH,NOPCNP,
     *                            CNNMIN_IFMESH,CNNMAX_IFMESH,NOPCNN
          END IF
C
      END IF 
C 
C=======================================================================
C=======================================================================
C     Central potential KAPPA parameters - ranges for mesh option
C=======================================================================
C=======================================================================
C 
      IF (KEYWOR.EQ.'CENTKAPMSH') THEN
C
          READ (N_READ,*) XKCMIN_IFMESH,XKCMAX_IFMESH,NOPXKC,            
     *                    XACMIN_IFMESH,XACMAX_IFMESH,NOPXAC,
     *                    XRCMIN_IFMESH,XRCMAX_IFMESH,NOPXRC,
     *                    XCoMIN_IFMESH,XCoMAX_IFMESH,NOPXCo
C
          IF (IFMESH.EQ.1) THEN
C
              WRITE(LOGFIL,'(A)') KEYWOR
C
              WRITE(LOGFIL,'(9X,''XKCMIN  XKCMAX  NOPXKC  '',
     *                          ''XACMIN  XACMAX  NOPXAC  '',
     *                          ''XRCMIN  XRCMAX  NOPXRC  '',
     *                          ''XCoMIN  XCoMAX  NOPXCo'',/,9X,
     *                            F6.2,2X,F6.2,2X,I6,2X,
     *                            F6.2,2X,F6.2,2X,I6,2X,
     *                            F6.2,2X,F6.2,2X,I6,2X,
     *                            F6.2,2X,F6.2,2X,I6,2X)')
     *
     *                            XKCMIN_IFMESH,XKCMAX_IFMESH,NOPXKC,            
     *                            XACMIN_IFMESH,XACMAX_IFMESH,NOPXAC,
     *                            XRCMIN_IFMESH,XRCMAX_IFMESH,NOPXRC,
     *                            XCoMIN_IFMESH,XCoMAX_IFMESH,NOPXCo
          END IF
C
      END IF          
C 
C=======================================================================
C=======================================================================
C     Spin-orbit potential KAPPA parameters - ranges for mesh option
C=======================================================================
C=======================================================================
C 
      IF (KEYWOR.EQ.'SORBKAPMSH') THEN
C
          READ (N_READ,*) XKSMIN_IFMESH,XKSMAX_IFMESH,NOPXKS,
     *                    XASMIN_IFMESH,XASMAX_IFMESH,NOPXAS,
     *                    XRSMIN_IFMESH,XRSMAX_IFMESH,NOPXRS
C
          IF (IFMESH.EQ.1) THEN
C
              WRITE(LOGFIL,'(A)') KEYWOR
C
              WRITE(LOGFIL,'(9X,''XKSMIN  XKSMAX  NOPXKS  '',
     *                          ''XASMIN  XASMAX  NOPXAS  '',
     *                          ''XRSMIN  XRSMAX  NOPXRS'',/,9X,
     *                            F6.2,2X,F6.2,2X,I6,2X,
     *                            F6.2,2X,F6.2,2X,I6,2X,
     *                            F6.2,2X,F6.2,2X,I6,2X)')
     *
     *                            XKSMIN_IFMESH,XKSMAX_IFMESH,NOPXKS,
     *                            XASMIN_IFMESH,XASMAX_IFMESH,NOPXAS,
     *                            XRSMIN_IFMESH,XRSMAX_IFMESH,NOPXRS
          END IF
C
      END IF          
C 
C=======================================================================
C=======================================================================
C     Reduced-mass term KAPPA parameters - ranges for mesh option
C=======================================================================
C=======================================================================
C 
      IF (KEYWOR.EQ.'EFFMKAPMSH') THEN
C
          READ (N_READ,*) XKEMIN_IFMESH,XKEMAX_IFMESH,NOPXKE,
     *                    XAEMIN_IFMESH,XAEMAX_IFMESH,NOPXAE,
     *                    XREMIN_IFMESH,XREMAX_IFMESH,NOPXRE
C
          IF (IFMESH.EQ.1) THEN
C
              WRITE(LOGFIL,'(A)') KEYWOR
C
              WRITE(LOGFIL,'(9X,''XKEMIN  XKEMAX  NOPXKE  '',
     *                          ''XAEMIN  XAEMAX  NOPXAE  '',
     *                          ''XREMIN  XREMAX  NOPXRE'',/,9X,
     *                            F6.2,2X,F6.2,2X,I6,2X,
     *                            F6.2,2X,F6.2,2X,I6,2X,
     *                            F6.2,2X,F6.2,2X,I6,2X)')
     *
     *                            XKEMIN_IFMESH,XKEMAX_IFMESH,NOPXKE,
     *                            XAEMIN_IFMESH,XAEMAX_IFMESH,NOPXAE,
     *                            XREMIN_IFMESH,XREMAX_IFMESH,NOPXRE
          END IF
C
      END IF     
          
C 
C=======================================================================
C=======================================================================
C     V_o, r_o and a_o for KAPPA parameters - ranges for mesh option
C=======================================================================
C=======================================================================
C 
      IF (KEYWOR.EQ.'CNTVRAMESH') THEN
C
          READ (N_READ,*) VC_MIN_IFMESH,VC_MAX_IFMESH,NOP_VC,
     *                    AC_MIN_IFMESH,AC_MAX_IFMESH,NOP_AC,
     *                    RC_MIN_IFMESH,RC_MAX_IFMESH,NOP_RC
C
          IF (IFMESH.EQ.1) THEN
C
              WRITE(LOGFIL,'(A)') KEYWOR
C
              WRITE(LOGFIL,'(9X,''V0CMIN  V0CMAX  NOPV0C  '',
     *                          ''A0CMIN  A0CMAX  NOPA0C  '',
     *                          ''R0CMIN  R0CMAX  NOPR0C  '',/,9X,
     *                            F6.1,2X,F6.2,2X,I6,2X,
     *                            F6.2,2X,F6.2,2X,I6,2X,
     *                            F6.2,2X,F6.2,2X,I6,2X)')
     *
     *                    VC_MIN_IFMESH,VC_MAX_IFMESH,NOP_VC,
     *                    AC_MIN_IFMESH,AC_MAX_IFMESH,NOP_AC,
     *                    RC_MIN_IFMESH,RC_MAX_IFMESH,NOP_RC
          END IF
C
      END IF    
C 
C=======================================================================
C=======================================================================
C  kappa V_o, r_o and a_o for KAPPA parameters - ranges for mesh option
C=======================================================================
C=======================================================================
C 
      IF (KEYWOR.EQ.'CNTKAPMESH') THEN
C
          READ (N_READ,*) VCKMIN_IFMESH,VCKMAX_IFMESH,NOPVCK,
     *                    ACKMIN_IFMESH,ACKMAX_IFMESH,NOPACK,
     *                    RCKMIN_IFMESH,RCKMAX_IFMESH,NOPRCK
C
          IF (IFMESH.EQ.1) THEN
C
              WRITE(LOGFIL,'(A)') KEYWOR
C
              WRITE(LOGFIL,'(9X,''V0CMIN  V0CMAX  NOPV0C  '',
     *                          ''A0CMIN  A0CMAX  NOPA0C  '',
     *                          ''R0CMIN  R0CMAX  NOPR0C  '',/,9X,
     *                            F6.1,2X,F6.2,2X,I6,2X,
     *                            F6.2,2X,F6.2,2X,I6,2X,
     *                            F6.2,2X,F6.2,2X,I6,2X)')
     *
     *                    VCKMIN_IFMESH,VCKMAX_IFMESH,NOPVCK,
     *                    ACKMIN_IFMESH,ACKMAX_IFMESH,NOPACK,
     *                    RCKMIN_IFMESH,RCKMAX_IFMESH,NOPRCK
          END IF
C
      END IF            
C 
C=======================================================================
C=======================================================================
C     V_o, r_o and a_o for KAPPA parameters - ranges for mesh option
C=======================================================================
C=======================================================================
C 
      IF (KEYWOR.EQ.'SORLRAMESH') THEN
C
          READ (N_READ,*) VS_MIN_IFMESH,VS_MAX_IFMESH,NOP_VS,
     *                    AS_MIN_IFMESH,AS_MAX_IFMESH,NOP_AS,
     *                    RS_MIN_IFMESH,RS_MAX_IFMESH,NOP_RS
C
          IF (IFMESH.EQ.1) THEN
C
              WRITE(LOGFIL,'(A)') KEYWOR
C
              WRITE(LOGFIL,'(9X,''V0CMIN  V0CMAX  NOPV0C  '',
     *                          ''A0CMIN  A0CMAX  NOPA0C  '',
     *                          ''R0CMIN  R0CMAX  NOPR0C  '',/,9X,
     *                            F6.1,2X,F6.2,2X,I6,2X,
     *                            F6.2,2X,F6.2,2X,I6,2X,
     *                            F6.2,2X,F6.2,2X,I6,2X)')
     *
     *                    VS_MIN_IFMESH,VS_MAX_IFMESH,NOP_VS,
     *                    AS_MIN_IFMESH,AS_MAX_IFMESH,NOP_AS,
     *                    RS_MIN_IFMESH,RS_MAX_IFMESH,NOP_RS
          END IF
C
      END IF    
C 
C=======================================================================
C=======================================================================
C  kappa V_o, r_o and a_o for KAPPA parameters - ranges for mesh option
C=======================================================================
C=======================================================================
C 
      IF (KEYWOR.EQ.'SORKAPMESH') THEN
C
          READ (N_READ,*) VSKMIN_IFMESH,VSKMAX_IFMESH,NOPVSK,
     *                    ASKMIN_IFMESH,ASKMAX_IFMESH,NOPASK,
     *                    RSKMIN_IFMESH,RSKMAX_IFMESH,NOPRSK
C
          IF (IFMESH.EQ.1) THEN
C
              WRITE(LOGFIL,'(A)') KEYWOR
C
              WRITE(LOGFIL,'(9X,''V0CMIN  V0CMAX  NOPV0C  '',
     *                          ''A0CMIN  A0CMAX  NOPA0C  '',
     *                          ''R0CMIN  R0CMAX  NOPR0C  '',/,9X,
     *                            F6.1,2X,F6.2,2X,I6,2X,
     *                            F6.2,2X,F6.2,2X,I6,2X,
     *                            F6.2,2X,F6.2,2X,I6,2X)')
     *
     *                    VSKMIN_IFMESH,VSKMAX_IFMESH,NOPVSK,
     *                    ASKMIN_IFMESH,ASKMAX_IFMESH,NOPASK,
     *                    RSKMIN_IFMESH,RSKMAX_IFMESH,NOPRSK
          END IF
C
      END IF            
C=======================================================================
C=======================================================================
C     Protons
C=======================================================================
C=======================================================================
C 
      IF (KEYWOR.EQ.'IFMESHPARP') THEN
C
          READ (*,*) IM_V0C_PROTON,IM_A0C_PROTON,IM_R0C_PROTON,
     *               IM_XLA_PROTON,IM_A0S_PROTON,IM_R0S_PROTON,
     *               IM_XRM_PROTON,IM_ARM_PROTON,IM_RRM_PROTON,
     *                                                  IM_COU
C
          IF (IFMESH.EQ.1) THEN
C
              WRITE(LOGFIL,'(A)') KEYWOR
C
              WRITE(LOGFIL,'(9X,''IM_V0C  IM_A0C  IM_R0C  IM_XLA  '',
     *                          ''IM_A0S  IM_R0S  IM_XRM  IM_ARM  '',
     *                          ''IM_RRM  IM_COU'')')
C
              WRITE(LOGFIL,'(9X,10(I6,2X))')
     *                       IM_V0C_PROTON,IM_A0C_PROTON,IM_R0C_PROTON,
     *                       IM_XLA_PROTON,IM_A0S_PROTON,IM_R0S_PROTON,
     *                       IM_XRM_PROTON,IM_ARM_PROTON,IM_RRM_PROTON,
     *                                                          IM_COU
          END IF
C_______________________________________________________________________
C 
C         The values of the vector index below are associated
C         with the potential parameters. We have:
C_______________________________________________________________________
C
C         1 -> central potential depth (p)
C         2 -> central potential radius (p)
C         3 -> central potential diffuseness (p)
C
C         4 -> spin-orbit strength (p)
C         5 -> spin-orbit potential radius (p)
C         6 -> spin-orbit potential diffuseness (p)
C
C         7 -> effective mass strength (p)
C         8 -> effective mass radius (p) 
C         9 -> effective mass diffuseness (p) 
C
C         10-> Coulomb radius parameter 
C_______________________________________________________________________
C 
C         Below, the we use the following coding: 
C 
C         "V" - variable, "STR" - starting value,
C         "MIN" -> minimal value of the interval,
C         "MAX" -> maximal value of the interval,
C         "STP" -> step to calculate derivatives
C_______________________________________________________________________
C
C         Minimum values (lower bounds of the intervals)
C
          XMIN_I(1)=V0CMIN_IFMESH
          XMIN_I(2)=R0CMIN_IFMESH
          XMIN_I(3)=A0CMIN_IFMESH
C          
          XMIN_I(4)=XL_MIN_IFMESH
          XMIN_I(5)=R0SMIN_IFMESH
          XMIN_I(6)=A0SMIN_IFMESH
C
          XMIN_I(7)=XEFMIN_IFMESH
          XMIN_I(8)=REFMIN_IFMESH
          XMIN_I(9)=AEFMIN_IFMESH
C_______________________________________________________________________
C
          XMIN_I(10)=CouMIN_IFMESH
C_______________________________________________________________________
C
C         Maximum values (upper bounds of the intervals)
C
          XMAX_I(1)=V0CMAX_IFMESH
          XMAX_I(2)=R0CMAX_IFMESH
          XMAX_I(3)=A0CMAX_IFMESH
C
          XMAX_I(4)=XL_MAX_IFMESH
          XMAX_I(5)=R0SMAX_IFMESH
          XMAX_I(6)=A0SMAX_IFMESH
C
          XMAX_I(7)=XEFMAX_IFMESH
          XMAX_I(8)=REFMAX_IFMESH
          XMAX_I(9)=AEFMAX_IFMESH
C_______________________________________________________________________
C
          XMAX_I(10)=CouMAX_IFMESH
C_______________________________________________________________________
C
C         Maximum number of points (user's choice)
C
          MAXPAR(1)=NOPV0C
          MAXPAR(2)=NOPR0C
          MAXPAR(3)=NOPA0C
C          
          MAXPAR(4)=NOP_XL
          MAXPAR(5)=NOPR0S
          MAXPAR(6)=NOPA0S
C
          MAXPAR(7)=NOPXEF
          MAXPAR(8)=NOPREF
          MAXPAR(9)=NOPAEF
C_______________________________________________________________________
C
          MAXPAR(10)=NOPCou
C_______________________________________________________________________
C 
C         Yes/No decision: take or not to take this is the question 
C 
          I_MESH(1)=IM_V0C_PROTON
          I_MESH(2)=IM_R0C_PROTON
          I_MESH(3)=IM_A0C_PROTON
C          
          I_MESH(4)=IM_XLA_PROTON
          I_MESH(5)=IM_R0S_PROTON
          I_MESH(6)=IM_A0S_PROTON
C          
          I_MESH(7)=IM_XRM_PROTON
          I_MESH(8)=IM_RRM_PROTON
          I_MESH(9)=IM_ARM_PROTON
C_______________________________________________________________________
C          
          I_MESH(10)=IM_COU
C_______________________________________________________________________
C
      END IF
C
C=======================================================================
C=======================================================================
C
      IF (KEYWOR.EQ.'IFMESHKAPP') THEN
C
          READ (*,*) IMKV0C_PROTON,IMKA0C_PROTON,IMKR0C_PROTON,
     *	             IMKXLA_PROTON,IMKA0S_PROTON,IMKR0S_PROTON,
     *               IMKXRM_PROTON,IMKARM_PROTON,IMKRRM_PROTON,
     *                                                  IMKCOU
C
          IF (IFMESH.EQ.1) THEN
C
              WRITE(LOGFIL,'(A)') KEYWOR
C
              WRITE(LOGFIL,'(9X,''IMKV0C  IMKA0C  IMKR0C  IMKXLA  '',
     *                          ''IMKA0S  IMKR0S  IMKXRM  IMKARM  '',
     *                          ''IMKRRM  IMKCOU'')')
C
              WRITE(LOGFIL,'(9X,10(I6,2X))')
     *                       IMKV0C_PROTON,IMKA0C_PROTON,IMKR0C_PROTON,
     *	                     IMKXLA_PROTON,IMKA0S_PROTON,IMKR0S_PROTON,
     *                       IMKXRM_PROTON,IMKARM_PROTON,IMKRRM_PROTON,
     *                                                          IMKCOU
          END IF
C_______________________________________________________________________
C 
C         The values of the vector index below are associated
C         with the potential parameters. We have:
C_______________________________________________________________________
C
C         11-> central potential depth       (p)
C         12-> central potential radius      (p)
C         13-> central potential diffuseness (p)
C
C         14-> spin-orbit strength              (p)
C         15-> spin-orbit potential radius      (p)
C         16-> spin-orbit potential diffuseness (p)
C
C         17-> effective mass strength    (p)
C         18-> effective mass radius      (p) 
C         19-> effective mass diffuseness (p) 
C
C         20-> Coulomb radius parameter 
C_______________________________________________________________________
C 
C         Below, the we use the following coding: 
C 
C         "V" - variable, "STR" - starting value,
C         "MIN" -> minimal value of the interval,
C         "MAX" -> maximal value of the interval,
C         "STP" -> step to calculate derivatives
C_______________________________________________________________________
C
C         Minimum values (lower bounds of the intervals)
C
          XMIN_I(11)=XKCMIN_IFMESH
          XMIN_I(12)=XRCMIN_IFMESH
          XMIN_I(13)=XACMIN_IFMESH
C          
          XMIN_I(14)=XKSMIN_IFMESH
          XMIN_I(15)=XRSMIN_IFMESH
          XMIN_I(16)=XASMIN_IFMESH
C
          XMIN_I(17)=XKEMIN_IFMESH
          XMIN_I(18)=XREMIN_IFMESH
          XMIN_I(19)=XAEMIN_IFMESH
C_______________________________________________________________________
C
          XMIN_I(20)=XCoMIN_IFMESH
C_______________________________________________________________________
C
C         Maximum values (upper bounds of the intervals)
C
          XMAX_I(11)=XKCMAX_IFMESH
          XMAX_I(12)=XRCMAX_IFMESH
          XMAX_I(13)=XACMAX_IFMESH
C
          XMAX_I(14)=XKSMAX_IFMESH
          XMAX_I(15)=XRSMAX_IFMESH
          XMAX_I(16)=XASMAX_IFMESH
C
          XMAX_I(17)=XKEMAX_IFMESH
          XMAX_I(18)=XREMAX_IFMESH
          XMAX_I(19)=XAEMAX_IFMESH
C_______________________________________________________________________
C
          XMAX_I(20)=XCoMAX_IFMESH
C_______________________________________________________________________
C
C         Maximum number of points (user's choice)
C
          MAXPAR(11)=NOPXKC
          MAXPAR(12)=NOPXRC
          MAXPAR(13)=NOPXAC
C          
          MAXPAR(14)=NOPXKS
          MAXPAR(15)=NOPXRS
          MAXPAR(16)=NOPXAS
C
          MAXPAR(17)=NOPXKE
          MAXPAR(18)=NOPXRE
          MAXPAR(19)=NOPXAE
C_______________________________________________________________________
C
          MAXPAR(20)=NOPXCo
C_______________________________________________________________________
C
C         Yes/No decision: take or not to take this KAPPA 
C 
          I_MESH(11)=IMKV0C_PROTON
          I_MESH(12)=IMKR0C_PROTON
          I_MESH(13)=IMKA0C_PROTON
C          
          I_MESH(14)=IMKXLA_PROTON
          I_MESH(15)=IMKR0S_PROTON
          I_MESH(16)=IMKA0S_PROTON
C          
          I_MESH(17)=IMKXRM_PROTON
          I_MESH(18)=IMKRRM_PROTON
          I_MESH(19)=IMKARM_PROTON
C_______________________________________________________________________
C          
          I_MESH(20)=IMKCOU
C_______________________________________________________________________
C
      END IF
C
C=======================================================================
C=======================================================================
C     Neutrons
C=======================================================================
C=======================================================================
C 
      IF (KEYWOR.EQ.'IFMESHPARN') THEN
C
          READ (*,*) IM_V0C_NEUTRS,IM_A0C_NEUTRS,IM_R0C_NEUTRS,
     *               IM_XLA_NEUTRS,IM_A0S_NEUTRS,IM_R0S_NEUTRS,
     *               IM_XRM_NEUTRS,IM_ARM_NEUTRS,IM_RRM_NEUTRS
C
          IF (IFMESH.EQ.1) THEN
C
              WRITE(LOGFIL,'(A)') KEYWOR
C
              WRITE(LOGFIL,'(9X,''IM_V0C  IM_A0C  IM_R0C  IM_XLA  '',
     *                          ''IM_A0S  IM_R0S  IM_XRM  IM_ARM  '',
     *                          ''IM_RRM  IM_COU'')')
C
              WRITE(LOGFIL,'(9X,9(I6,2X))')
     *                       IM_V0C_NEUTRS,IM_A0C_NEUTRS,IM_R0C_NEUTRS,
     *                       IM_XLA_NEUTRS,IM_A0S_NEUTRS,IM_R0S_NEUTRS,
     *                       IM_XRM_NEUTRS,IM_ARM_NEUTRS,IM_RRM_NEUTRS
          END IF
C_______________________________________________________________________
C 
C         The values of the vector index below are associated
C         with the potential parameters. We have:
C_______________________________________________________________________
C
C         21 -> central potential depth       (n)
C         22 -> central potential radius      (n)
C         23 -> central potential diffuseness (n)
C
C         24 -> spin-orbit strength              (n)
C         25 -> spin-orbit potential radius      (n)
C         26 -> spin-orbit potential diffuseness (n)
C
C         27 -> effective mass strength    (n)
C         28 -> effective mass radius      (n) 
C         29 -> effective mass diffuseness (n) 
C_______________________________________________________________________
C 
C         Below, the we use the following coding: 
C 
C         "V" - variable, "STR" - starting value,
C         "MIN" -> minimal value of the interval,
C         "MAX" -> maximal value of the interval,
C         "STP" -> step to calculate derivatives
C_______________________________________________________________________
C
C         Minimum values (lower bounds of the intervals)
C
          XMIN_I(21)=V0CMIN_IFMESH
          XMIN_I(22)=R0CMIN_IFMESH
          XMIN_I(23)=A0CMIN_IFMESH
C          
          XMIN_I(24)=XL_MIN_IFMESH
          XMIN_I(25)=R0SMIN_IFMESH
          XMIN_I(26)=A0SMIN_IFMESH
C
          XMIN_I(27)=XEFMIN_IFMESH
          XMIN_I(28)=REFMIN_IFMESH
          XMIN_I(29)=AEFMIN_IFMESH
C_______________________________________________________________________
C
C         Maximum values (upper bounds of the intervals)
C
          XMAX_I(21)=V0CMAX_IFMESH
          XMAX_I(22)=R0CMAX_IFMESH
          XMAX_I(23)=A0CMAX_IFMESH
C
          XMAX_I(24)=XL_MAX_IFMESH
          XMAX_I(25)=R0SMAX_IFMESH
          XMAX_I(26)=A0SMAX_IFMESH
C
          XMAX_I(27)=XEFMAX_IFMESH
          XMAX_I(28)=REFMAX_IFMESH
          XMAX_I(29)=AEFMAX_IFMESH
C_______________________________________________________________________
C
C         Maximum number of points (user's choice)
C
          MAXPAR(21)=NOPV0C
          MAXPAR(22)=NOPR0C
          MAXPAR(23)=NOPA0C
C          
          MAXPAR(24)=NOP_XL
          MAXPAR(25)=NOPR0S
          MAXPAR(26)=NOPA0S
C
          MAXPAR(27)=NOPXEF
          MAXPAR(28)=NOPREF
          MAXPAR(29)=NOPAEF
C_______________________________________________________________________
C 
C         Yes/No decision: take or not to take this is the question 
C 
          I_MESH(21)=IM_V0C_NEUTRS
          I_MESH(22)=IM_R0C_NEUTRS
          I_MESH(23)=IM_A0C_NEUTRS
C          
          I_MESH(24)=IM_XLA_NEUTRS
          I_MESH(25)=IM_R0S_NEUTRS
          I_MESH(26)=IM_A0S_NEUTRS
C          
          I_MESH(27)=IM_XRM_NEUTRS
          I_MESH(28)=IM_RRM_NEUTRS
          I_MESH(29)=IM_ARM_NEUTRS
C_______________________________________________________________________
C
      END IF
C
C=======================================================================
C=======================================================================
C 
      IF (KEYWOR.EQ.'IFMESHKAPN') THEN
C
          READ (*,*) IMKV0C_NEUTRS,IMKA0C_NEUTRS,IMKR0C_NEUTRS,
     *	             IMKXLA_NEUTRS,IMKA0S_NEUTRS,IMKR0S_NEUTRS,
     *               IMKXRM_NEUTRS,IMKARM_NEUTRS,IMKRRM_NEUTRS 
C
          IF (IFMESH.EQ.1) THEN

              WRITE(LOGFIL,'(A)') KEYWOR

              WRITE(LOGFIL,'(9X,''IMKV0C  IMKA0C  IMKR0C  IMKXLA  '',
     *                          ''IMKA0S  IMKR0S  IMKXRM  IMKARM  '',
     *                          ''IMKRRM  '')')
C
              WRITE(LOGFIL,'(9X,9(I6,2X))')
     *                       IMKV0C_NEUTRS,IMKA0C_NEUTRS,IMKR0C_NEUTRS,
     *	                     IMKXLA_NEUTRS,IMKA0S_NEUTRS,IMKR0S_NEUTRS,
     *                       IMKXRM_NEUTRS,IMKARM_NEUTRS,IMKRRM_NEUTRS 
          END IF
C_______________________________________________________________________
C 
C         The values of the vector index below are associated
C         with the potential parameters. We have:
C_______________________________________________________________________
C
C         30-> central potential depth       (p)
C         31-> central potential radius      (p)
C         32-> central potential diffuseness (p)
C
C         33-> spin-orbit strength              (p)
C         34-> spin-orbit potential radius      (p)
C         35-> spin-orbit potential diffuseness (p)
C
C         36-> effective mass strength    (p)
C         37-> effective mass radius      (p) 
C         38-> effective mass diffuseness (p) 
C_______________________________________________________________________
C 
C         Below, the we use the following coding: 
C 
C         "V" - variable, "STR" - starting value,
C         "MIN" -> minimal value of the interval,
C         "MAX" -> maximal value of the interval,
C         "STP" -> step to calculate derivatives
C_______________________________________________________________________
C
C         Minimum values (lower bounds of the intervals)
C
          XMIN_I(30)=XKCMIN_IFMESH
          XMIN_I(31)=XRCMIN_IFMESH
          XMIN_I(32)=XACMIN_IFMESH
C          
          XMIN_I(33)=XKSMIN_IFMESH
          XMIN_I(34)=XRSMIN_IFMESH
          XMIN_I(35)=XASMIN_IFMESH
C
          XMIN_I(36)=XKEMIN_IFMESH
          XMIN_I(37)=XREMIN_IFMESH
          XMIN_I(38)=XAEMIN_IFMESH
C_______________________________________________________________________
C
C         Maximum values (upper bounds of the intervals)
C
          XMAX_I(30)=XKCMAX_IFMESH
          XMAX_I(31)=XRCMAX_IFMESH
          XMAX_I(32)=XACMAX_IFMESH
C
          XMAX_I(33)=XKSMAX_IFMESH
          XMAX_I(34)=XRSMAX_IFMESH
          XMAX_I(35)=XASMAX_IFMESH
C
          XMAX_I(36)=XKEMAX_IFMESH
          XMAX_I(37)=XREMAX_IFMESH
          XMAX_I(38)=XAEMAX_IFMESH
C_______________________________________________________________________
C
C         Maximum number of points (user's choice)
C
          MAXPAR(30)=NOPXKC
          MAXPAR(31)=NOPXRC
          MAXPAR(32)=NOPXAC
C          
          MAXPAR(33)=NOPXKS
          MAXPAR(34)=NOPXRS
          MAXPAR(35)=NOPXAS
C
          MAXPAR(36)=NOPXKE
          MAXPAR(37)=NOPXRE
          MAXPAR(38)=NOPXAE
C_______________________________________________________________________
C 
C         Yes/No decision: take or not to take this KAPPA 
C 
          I_MESH(30)=IMKV0C_NEUTRS
          I_MESH(31)=IMKR0C_NEUTRS
          I_MESH(32)=IMKA0C_NEUTRS
C          
          I_MESH(33)=IMKXLA_NEUTRS
          I_MESH(34)=IMKR0S_NEUTRS
          I_MESH(35)=IMKA0S_NEUTRS
C          
          I_MESH(36)=IMKXRM_NEUTRS
          I_MESH(37)=IMKRRM_NEUTRS
          I_MESH(38)=IMKARM_NEUTRS
C          
      END IF 
C
C=======================================================================
C=======================================================================
C     Density-dependent spin-orbit Hamiltonian
C=======================================================================
C=======================================================================
C 
      IF (KEYWOR.EQ.'IFMESHDENS') THEN
C
          READ (*,*) IM_XPP,IM_XPN,IM_XNP,IM_XNN
C
          IF (IFDENS.EQ.1 .AND. IFMESH.EQ.1) THEN
C
              WRITE(LOGFIL,'(A)') KEYWOR
C
              WRITE(LOGFIL,'(9X,''IM_XPP  IM_XPN  IM_XNP  IM_XNN'')')
              WRITE(LOGFIL,'(9X,4(I6,2X))')IM_XPP,IM_XPN,IM_XNP,IM_XNN
C
          END IF
C_______________________________________________________________________
C 
C         The values of the vector index below are associated
C         with the potential parameters. We have:
C_______________________________________________________________________
C
C         39 -> lambda  proton-proton
C         40 -> lambda  proton-neutron
C         41 -> lambda neutron-proton
C         42 -> lambda neutron-neutron
C ###     43 -> alpha spin-orbit
C ###     44 -> beta spin-orbit
C_______________________________________________________________________
C 
C         Below, the we use the following coding: 
C 
C         "V" - variable, "STR" - starting value,
C         "MIN" -> minimal value of the interval,
C         "MAX" -> maximal value of the interval,
C         "STP" -> step to calculate derivatives
C_______________________________________________________________________
C
C         Minimum values (lower bounds of the intervals)
C
          XMIN_I(39)=XPPMIN_IFMESH
          XMIN_I(40)=XPNMIN_IFMESH
          XMIN_I(41)=XNPMIN_IFMESH
          XMIN_I(42)=XNNMIN_IFMESH
C_______________________________________________________________________
C
C         Maximum values (upper bounds of the intervals)
C
          XMAX_I(39)=XPPMAX_IFMESH
          XMAX_I(40)=XPNMAX_IFMESH
          XMAX_I(41)=XNPMAX_IFMESH
          XMAX_I(42)=XNNMAX_IFMESH
C_______________________________________________________________________
C
C         Maximum number of points (user's choice)
C
          MAXPAR(39)=NOPXPP
          MAXPAR(40)=NOPXPN
          MAXPAR(41)=NOPXNP
          MAXPAR(42)=NOPXNN          
C_______________________________________________________________________
C 
C         Yes/No decision: take or not to take this is the question 
C 
          I_MESH(39)=IM_XPP
          I_MESH(40)=IM_XPN
          I_MESH(41)=IM_XNP
          I_MESH(42)=IM_XNN
C_______________________________________________________________________
C
      END IF
C
C=======================================================================
C=======================================================================
C     SO-Tensor part
C=======================================================================
C=======================================================================
C 
      IF (KEYWOR.EQ.'IFMESHTENS') THEN
C
          READ (*,*) IM_YPP,IM_YPN,IM_YNP,IM_YNN
C
          IF (IFTENS.EQ.1 .AND. ISORBT.EQ.1 .AND. IFMESH.EQ.1) THEN
C
              WRITE(LOGFIL,'(A)') KEYWOR

              WRITE(LOGFIL,'(9X,''IM_YPP  IM_YPN  IM_YNP  IM_YNN'')')
              WRITE(LOGFIL,'(9X,4(I6,2X))') IM_YPP,IM_YPN,IM_YNP,IM_YNN
C
          END IF
C_______________________________________________________________________
C 
C         The values of the vector index below are associated
C         with the potential parameters. We have:
C_______________________________________________________________________
C
C         45 -> tensor lambda  proton-proton
C         46 -> tensor lambda  proton-neutron
C         47 -> tensor lambda neutron-proton
C         48 -> tensor lambda neutron-neutron
C_______________________________________________________________________
C 
C         Below, the we use the following coding: 
C 
C         "V" - variable, "STR" - starting value,
C         "MIN" -> minimal value of the interval,
C         "MAX" -> maximal value of the interval,
C         "STP" -> step to calculate derivatives
C_______________________________________________________________________
C
C         Minimum values (lower bounds of the intervals)
C
          XMIN_I(43)=YPPMIN_IFMESH
          XMIN_I(44)=YPNMIN_IFMESH
          XMIN_I(45)=YNPMIN_IFMESH
          XMIN_I(46)=YNNMIN_IFMESH
C_______________________________________________________________________
C
C         Maximum values (upper bounds of the intervals)
C
          XMAX_I(43)=YPPMAX_IFMESH
          XMAX_I(44)=YPNMAX_IFMESH
          XMAX_I(45)=YNPMAX_IFMESH
          XMAX_I(46)=YNNMAX_IFMESH
C_______________________________________________________________________
C
C         Maximum number of points (user's choice)
C
          MAXPAR(43)=NOPYPP
          MAXPAR(44)=NOPYPN
          MAXPAR(45)=NOPYNP
          MAXPAR(46)=NOPYNN            
C_______________________________________________________________________
C 
C         Yes/No decision: take or not to take this is the question 
C 
          I_MESH(43)=IM_YPP
          I_MESH(44)=IM_YPN
          I_MESH(45)=IM_YNP
          I_MESH(46)=IM_YNN
C_______________________________________________________________________
C
      END IF
C
C=======================================================================
C=======================================================================
C     Central-potential tensor part
C=======================================================================
C=======================================================================
C 
      IF (KEYWOR.EQ.'IFMESHCNTN') THEN
C
          READ (*,*) IM_CPP,IM_CPN,IM_CNP,IM_CNN
C
          IF (IFTENS.EQ.1 .AND. ICENTT.EQ.1 .AND. IFMESH.EQ.1) THEN
C
              WRITE(LOGFIL,'(A)') KEYWOR
C
              WRITE(LOGFIL,'(9X,''IM_CPP  IM_CPN  IM_CNP  IM_CNN'')')
              WRITE(LOGFIL,'(9X,4(I6,2X))')IM_CPP,IM_CPN,IM_CNP,IM_CNN
C
          END IF
C_______________________________________________________________________
C 
C         The values of the vector index below are associated
C         with the potential parameters. We have:
C_______________________________________________________________________
C
C         45 -> tensor lambda  proton-proton
C         46 -> tensor lambda  proton-neutron
C         47 -> tensor lambda neutron-proton
C         48 -> tensor lambda neutron-neutron
C_______________________________________________________________________
C 
C         Below, the we use the following coding: 
C 
C         "V" - variable, "STR" - starting value,
C         "MIN" -> minimal value of the interval,
C         "MAX" -> maximal value of the interval,
C         "STP" -> step to calculate derivatives
C_______________________________________________________________________
C
C         Minimum values (lower bounds of the intervals)
C
          XMIN_I(47)=CPPMIN_IFMESH
          XMIN_I(48)=CPNMIN_IFMESH
          XMIN_I(49)=CNPMIN_IFMESH
          XMIN_I(50)=CNNMIN_IFMESH
C_______________________________________________________________________
C
C         Maximum values (upper bounds of the intervals)
C
          XMAX_I(47)=CPPMAX_IFMESH
          XMAX_I(48)=CPNMAX_IFMESH
          XMAX_I(49)=CNPMAX_IFMESH
          XMAX_I(50)=CNNMAX_IFMESH
C_______________________________________________________________________
C
C         Maximum number of points (user's choice)
C
          MAXPAR(47)=NOPCPP
          MAXPAR(48)=NOPCPN
          MAXPAR(49)=NOPCNP
          MAXPAR(50)=NOPCNN            
C_______________________________________________________________________
C 
C         Yes/No decision: take or not to take this is the question 
C 
          I_MESH(47)=IM_CPP
          I_MESH(48)=IM_CPN
          I_MESH(49)=IM_CNP
          I_MESH(50)=IM_CNN
C_______________________________________________________________________
C
      END IF
C
C=======================================================================
C=======================================================================
C     Central-potential KAPPA parametrization
C=======================================================================
C=======================================================================
C 
      IF (KEYWOR.EQ.'IFMESHCNTK') THEN
C
          READ (*,*) IM_V0C,IM_KVC,IM_R0C,IM_KRC,IM_A0C,IM_KAC
C
          IF (IFMESH.EQ.1) THEN
C
              WRITE(LOGFIL,'(A)') KEYWOR
C
              WRITE(LOGFIL,'(9X,''IM_V0C  IM_KVC  IM_R0C  IM_KRC  '',
     *                          ''IM_A0C  IM_KAC'')')
              WRITE(LOGFIL,'(9X,4(I6,2X))')IM_V0C,IM_KVC,IM_R0C,
     *                                     IM_KRC,IM_A0C,IM_KAC
C
          END IF
C_______________________________________________________________________
C 
C         The values of the vector index below are associated
C         with the potential parameters. We have:
C_______________________________________________________________________
C
C         51 -> Vo central
C         52 -> kappa Vo central
C         53 -> ro central
C         54 -> kappa ro central
C         55 -> ao central
C         56 -> kappa ao central
C_______________________________________________________________________
C 
C         Below, the we use the following coding: 
C 
C         "V" - variable, "STR" - starting value,
C         "MIN" -> minimal value of the interval,
C         "MAX" -> maximal value of the interval,
C         "STP" -> step to calculate derivatives
C_______________________________________________________________________
C
C         Minimum values (lower bounds of the intervals)
C
          XMIN_I(51)=VC_MIN_IFMESH
          XMIN_I(52)=VCKMIN_IFMESH
          XMIN_I(53)=RC_MIN_IFMESH
          XMIN_I(54)=RCKMIN_IFMESH
          XMIN_I(55)=AC_MIN_IFMESH
          XMIN_I(56)=ACKMIN_IFMESH
C_______________________________________________________________________
C
C         Maximum values (upper bounds of the intervals)
C
          XMAX_I(51)=VC_MAX_IFMESH
          XMAX_I(52)=VCKMAX_IFMESH
          XMAX_I(53)=RC_MAX_IFMESH
          XMAX_I(54)=RCKMAX_IFMESH
          XMAX_I(55)=AC_MAX_IFMESH
          XMAX_I(56)=ACKMAX_IFMESH
C_______________________________________________________________________
C
C         Maximum number of points (user's choice)
C
          MAXPAR(51)=NOP_VC
          MAXPAR(52)=NOPVCK
          MAXPAR(53)=NOP_RC
          MAXPAR(54)=NOPRCK
          MAXPAR(55)=NOP_AC
          MAXPAR(56)=NOPACK      
C_______________________________________________________________________
C 
C         Yes/No decision: take or not to take this is the question 
C 
          I_MESH(51)=IM_V0C
          I_MESH(52)=IM_KVC
          I_MESH(53)=IM_R0C
          I_MESH(54)=IM_KRC
          I_MESH(55)=IM_A0C
          I_MESH(56)=IM_KAC
C_______________________________________________________________________
C
      END IF
C
C=======================================================================
C=======================================================================
C     Central-potential KAPPA parametrization
C=======================================================================
C=======================================================================
C 
      IF (KEYWOR.EQ.'IFMESHSORK') THEN
C
          READ (*,*) IM_VSO,IM_KVS,IM_RSO,IM_KRS,IM_ASO,IM_KAS
C
          IF (IFMESH.EQ.1) THEN
C
              WRITE(LOGFIL,'(A)') KEYWOR
C
              WRITE(LOGFIL,'(9X,''IM_VSO  IM_KVS  IM_RSO  IM_KRS  '',
     *                          ''IM_ASO  IM_KAS  '')')
              WRITE(LOGFIL,'(9X,4(I6,2X))')IM_VSO,IM_KVS,IM_RSO,
     *                                     IM_KRS,IM_ASO,IM_KAS
C
          END IF
C_______________________________________________________________________
C 
C         The values of the vector index below are associated
C         with the potential parameters. We have:
C_______________________________________________________________________
C
C         57 -> Vo spin-orbit
C         58 -> kappa Vo spin-orbit
C         59 -> ro spin-orbit
C         60 -> kappa ro spin-orbit
C         61 -> ao spin-orbit
C         62 -> kappa ao spin-orbit
C_______________________________________________________________________
C 
C         Below, the we use the following coding: 
C 
C         "V" - variable, "STR" - starting value,
C         "MIN" -> minimal value of the interval,
C         "MAX" -> maximal value of the interval,
C         "STP" -> step to calculate derivatives
C_______________________________________________________________________
C
C         Minimum values (lower bounds of the intervals)
C
          XMIN_I(57)=VS_MIN_IFMESH
          XMIN_I(58)=VSKMIN_IFMESH
          XMIN_I(59)=RS_MIN_IFMESH
          XMIN_I(60)=RSKMIN_IFMESH
          XMIN_I(61)=AS_MIN_IFMESH
          XMIN_I(62)=ASKMIN_IFMESH
C_______________________________________________________________________
C
C         Maximum values (upper bounds of the intervals)
C
          XMAX_I(57)=VS_MAX_IFMESH
          XMAX_I(58)=VSKMAX_IFMESH
          XMAX_I(59)=RS_MAX_IFMESH
          XMAX_I(60)=RSKMAX_IFMESH
          XMAX_I(61)=AS_MAX_IFMESH
          XMAX_I(62)=ASKMAX_IFMESH
C_______________________________________________________________________
C
C         Maximum number of points (user's choice)
C
          MAXPAR(57)=NOP_VS
          MAXPAR(58)=NOPVSK
          MAXPAR(59)=NOP_RS
          MAXPAR(60)=NOPRSK
          MAXPAR(61)=NOP_AS
          MAXPAR(62)=NOPASK      
C_______________________________________________________________________
C 
C         Yes/No decision: take or not to take this is the question 
C 
          I_MESH(57)=IM_VSO
          I_MESH(58)=IM_KVS
          I_MESH(59)=IM_RSO
          I_MESH(60)=IM_KRS
          I_MESH(61)=IM_ASO
          I_MESH(62)=IM_KAS
C_______________________________________________________________________
C
      END IF
C
C=======================================================================
C=======================================================================
C     Entering to the Monte-Carlo Part. 
C
C     Reading the number of M-C restarts and the number of bins
C     that we want for the future histograms.
C=======================================================================
C=======================================================================
C
      IF (KEYWOR.EQ.'MONTECARLO') THEN
C
          READ(*,*) IFPSEU,IFPARA,LDMONT,LDBINS
C
          IF (IFMOCA.EQ.1) THEN
              WRITE(LOGFIL,'(A)') KEYWOR
              WRITE(LOGFIL,'(9X,''LDMONT  LDBINS'',/,9X,I6,2X,I6)')
     *                            LDMONT, LDBINS
          END IF
C
          IF (LDMONT.GT.NDMONT) THEN
              WRITE(LOGFIL,'(''Stop in NAMELI: '',
     *                       ''LDMONT= '',I6,'' is GT NDMONT= '',I6)')
     *                         LDMONT,                NDMONT
              STOP 'STOP in NAMELI: LDMONT.GT.NDMONT'
          END IF
C
          IF (LDMONT.GT.NDRAND) THEN
              WRITE(NOUTPT,'(''LDMONT='',I6,'' exceeds NDRAND='',I6)')
     *                         LDMONT,                 NDRAND
              STOP 'STOP in NAMELI: LDMONT.GT.NDRAND'
          END IF
C
          IF (LDBINS.GT.NDBINS) THEN
              WRITE(LOGFIL,'(''Stop in NAMELI: '',
     *                       ''LDBINS= '',I6,'' is GT NDBINS= '',I6)')
     *                         LDBINS,                NDBINS
              STOP 'STOP in NAMELI: LDBINS.GT.NDBINS'
          END IF
C
      END IF
C
C=======================================================================
C=======================================================================
C     Reading the mean values of the parameters
C=======================================================================
C=======================================================================
C
      IF (KEYWOR.EQ.'CENTPROT_M') THEN
          READ(*,*) VPCENT_XMEANS,RPCENT_XMEANS,APCENT_XMEANS,
     *              XK_VPC_XMEANS,XK_RPC_XMEANS,XK_APC_XMEANS,
     *              R0COUL_XMEANS,XK_COU_XMEANS
          IF (IFMOCA.EQ.1) THEN
              WRITE(LOGFIL,'(A)') KEYWOR
              WRITE(LOGFIL,'(9X,''V0CENT  R0CENT  A0CENT  '',
     *                          ''XK_V0C  XK_R0C  XK_A0C  '',
     *                          ''R0COUL  XK_COU'')')
              WRITE(LOGFIL,'(9X,8(F6.2,2X))')
     *                            VPCENT_XMEANS,RPCENT_XMEANS,
     *                                          APCENT_XMEANS,
     *                            XK_VPC_XMEANS,XK_RPC_XMEANS,
     *                                          XK_APC_XMEANS,
     *                            R0COUL_XMEANS,XK_COU_XMEANS
          END IF
      END IF
C
C=======================================================================
C
      IF (KEYWOR.EQ.'CENTNEUT_M') THEN
          READ(*,*) VNCENT_XMEANS,RNCENT_XMEANS,ANCENT_XMEANS,
     *              XK_VNC_XMEANS,XK_RNC_XMEANS,XK_ANC_XMEANS
          IF (IFMOCA.EQ.1) THEN
              WRITE(LOGFIL,'(A)') KEYWOR
              WRITE(LOGFIL,'(9X,''V0CENT  R0CENT  A0CENT  '',
     *                          ''XK_V0C  XK_R0C  XK_A0C  '',/,
     *                           9X,6(F6.2,2X))')
     *                            VNCENT_XMEANS,RNCENT_XMEANS,
     *                                          ANCENT_XMEANS,
     *                            XK_VNC_XMEANS,XK_RNC_XMEANS,
     *                                          XK_ANC_XMEANS
          END IF
      END IF
C
C=======================================================================
C
      IF (KEYWOR.EQ.'SORBPROT_M') THEN
          READ(*,*) VPSORB_XMEANS,RPSORB_XMEANS,APSORB_XMEANS,
     *              XK_VPS_XMEANS,XK_RPS_XMEANS,XK_APS_XMEANS
          IF (IFMOCA.EQ.1) THEN
              WRITE(LOGFIL,'(A)') KEYWOR
              WRITE(LOGFIL,'(9X,''V0SORB  R0SORB  A0SORB  '',
     *                          ''XK_V0S  XK_R0S  XK_A0S  '',/,
     *                          9X,6(F6.2,2X))')
     *                     VPSORB_XMEANS,RPSORB_XMEANS,APSORB_XMEANS,
     *                     XK_VPS_XMEANS,XK_RPS_XMEANS,XK_APS_XMEANS
          END IF
      END IF
C
C=======================================================================
C
      IF (KEYWOR.EQ.'SORBNEUT_M') THEN
          READ(*,*) VNSORB_XMEANS,RNSORB_XMEANS,ANSORB_XMEANS,
     *              XK_VNS_XMEANS,XK_RNS_XMEANS,XK_ANS_XMEANS
          IF (IFMOCA.EQ.1) THEN
              WRITE(LOGFIL,'(A)') KEYWOR
              WRITE(LOGFIL,'(9X,''V0SORB  R0SORB  A0SORB  '',
     *                          ''XK_V0S  XK_R0S  XK_A0S  '',/,
     *                          9X,6(F6.2,2X))')
     *                     VNSORB_XMEANS,RNSORB_XMEANS,ANSORB_XMEANS,
     *                     XK_VNS_XMEANS,XK_RNS_XMEANS,XK_ANS_XMEANS
          END IF
      END IF
C
C=======================================================================
C
      IF (KEYWOR.EQ.'EFFEPROT_M') THEN
          READ(*,*) VPEFFM_XMEANS,RPEFFM_XMEANS,APEFFM_XMEANS,
     *              XK_VPE_XMEANS,XK_RPE_XMEANS,XK_APE_XMEANS
          IF (IFMOCA.EQ.1) THEN
              WRITE(LOGFIL,'(A)') KEYWOR
              WRITE(LOGFIL,'(9X,''V0EFFM  R0EFFM  A0EFFM  '',
     *                          ''XK_LEF  XK_REF  XK_AEF  '',/,
     *                          9X,6(F6.2,2X))')
     *                     VPEFFM_XMEANS,RPEFFM_XMEANS,APEFFM_XMEANS,
     *                     XK_VPE_XMEANS,XK_RPE_XMEANS,XK_APE_XMEANS
          END IF
      END IF
C
C=======================================================================
C
      IF (KEYWOR.EQ.'EFFENEUT_M') THEN
          READ(*,*) VNEFFM_XMEANS,RNEFFM_XMEANS,ANEFFM_XMEANS,
     *              XK_VNE_XMEANS,XK_RNE_XMEANS,XK_ANE_XMEANS
          IF (IFMOCA.EQ.1) THEN
              WRITE(LOGFIL,'(A)') KEYWOR
              WRITE(LOGFIL,'(9X,''V0EFFM  R0EFFM  A0EFFM  '',
     *                          ''XK_LEF  XK_REF  XK_AEF  '',/,
     *                          9X,6(F6.2,2X))')
     *                     VNEFFM_XMEANS,RNEFFM_XMEANS,ANEFFM_XMEANS,
     *                     XK_VNE_XMEANS,XK_RNE_XMEANS,XK_ANE_XMEANS
          END IF
      END IF
C
C=======================================================================
C
      IF (KEYWOR.EQ.'SORBDENS_M') THEN
          READ(*,*) ALAMPP_XMEANS,ALAMPN_XMEANS,ALAMNP_XMEANS,
     *                                          ALAMNN_XMEANS
          IF (IFMOCA.EQ.1) THEN
              WRITE(LOGFIL,'(A)') KEYWOR
              WRITE(LOGFIL,'(9X,''ALAMPP  ALAMPN  ALAMNP  ALAMNN'',/,
     *                       9X,4(F6.2,2X))')
     *                     ALAMPP_XMEANS,ALAMPN_XMEANS,
     *                     ALAMNP_XMEANS,ALAMNN_XMEANS
          END IF
      END IF
C
C=======================================================================
C
      IF (KEYWOR.EQ.'SORBTENS_M') THEN
          READ(*,*) TLAMPP_XMEANS,TLAMPN_XMEANS,TLAMNP_XMEANS,
     *                                          TLAMNN_XMEANS
          IF (IFMOCA.EQ.1) THEN
              WRITE(LOGFIL,'(A)') KEYWOR
              WRITE(LOGFIL,'(9X,''TLAMPP  TLAMPN  TLAMNP  TLAMNN'',/,
     *                       9X,4(F6.2,2X))')
     *                     TLAMPP_XMEANS,TLAMPN_XMEANS,
     *                     TLAMNP_XMEANS,TLAMNN_XMEANS
          END IF
      END IF
C
C=======================================================================
C
      IF (KEYWOR.EQ.'CENTTENS_M') THEN
          READ(*,*) CLAMPP_XMEANS,CLAMPN_XMEANS,CLAMNP_XMEANS,
     *                                          CLAMNN_XMEANS
          IF (IFMOCA.EQ.1) THEN
              WRITE(LOGFIL,'(A)') KEYWOR
              WRITE(LOGFIL,'(9X,''TLAMPP  TLAMPN  TLAMNP  TLAMNN'',/,
     *                       9X,4(F6.2,2X))')
     *                     CLAMPP_XMEANS,CLAMPN_XMEANS,
     *                     CLAMNP_XMEANS,CLAMNN_XMEANS
          END IF
      END IF
C
C=======================================================================
C
      IF (KEYWOR.EQ.'KAPACENT_M') THEN
          READ(*,*) V0CENT_XMEANS,XK_V0C_XMEANS,R0CENT_XMEANS,
     *              XK_R0C_XMEANS,A0CENT_XMEANS,XK_A0C_XMEANS
          IF (IFMOCA.EQ.1) THEN
              WRITE(LOGFIL,'(A)') KEYWOR
              WRITE(LOGFIL,'(9X,''V0CENT  XK_V0C  R0CENT  XK_R0C  '',
     *                          ''A0CENT  XK_A0C'',/,
     *                       9X,6(F6.2,2X))')
     *              V0CENT_XMEANS,XK_V0C_XMEANS,R0CENT_XMEANS,
     *              XK_R0C_XMEANS,A0CENT_XMEANS,XK_A0C_XMEANS
          END IF
      END IF
C
C=======================================================================
C
      IF (KEYWOR.EQ.'KAPASORB_M') THEN
          READ(*,*) V0SORB_XMEANS,XK_V0S_XMEANS,R0SORB_XMEANS,
     *              XK_R0S_XMEANS,A0SORB_XMEANS,XK_A0S_XMEANS
          IF (IFMOCA.EQ.1) THEN
              WRITE(LOGFIL,'(A)') KEYWOR
              WRITE(LOGFIL,'(9X,''V0CENT  XK_V0C  R0CENT  XK_R0C  '',
     *                          ''A0CENT  XK_A0C'',/,
     *                       9X,6(F6.2,2X))')
     *              V0SORB_XMEANS,XK_V0S_XMEANS,R0SORB_XMEANS,
     *              XK_R0s_XMEANS,A0SORB_XMEANS,XK_A0S_XMEANS
          END IF
      END IF
C
C=======================================================================
C=======================================================================
C     Reading the sigma values of the parameters
C=======================================================================
C=======================================================================
C
      IF (KEYWOR.EQ.'CENTPROT_S') THEN
          READ(*,*) VPCENT_SIGMAS,RPCENT_SIGMAS,APCENT_SIGMAS,
     *              XK_VPC_SIGMAS,XK_RPC_SIGMAS,XK_APC_SIGMAS,
     *              R0COUL_SIGMAS,XK_COU_SIGMAS
          IF (IFMOCA.EQ.1) THEN
              WRITE(LOGFIL,'(A)') KEYWOR
              WRITE(LOGFIL,'(9X,''V0CENT  R0CENT  A0CENT  '',
     *                          ''XK_V0C  XK_R0C  XK_A0C  '',
     *                          ''R0COUL  XK_COU'',/,9X,8(F6.2,2X))')
     *                            VPCENT_SIGMAS,RPCENT_SIGMAS,
     *                                          APCENT_SIGMAS,
     *                            XK_VPC_SIGMAS,XK_RPC_SIGMAS,
     *                                          XK_APC_SIGMAS,
     *                            R0COUL_SIGMAS,XK_COU_SIGMAS
          END IF
      END IF
C
C=======================================================================
C
      IF (KEYWOR.EQ.'CENTNEUT_S') THEN
          READ(*,*) VNCENT_SIGMAS,RNCENT_SIGMAS,ANCENT_SIGMAS,
     *              XK_VNC_SIGMAS,XK_RNC_SIGMAS,XK_ANC_SIGMAS
          IF (IFMOCA.EQ.1) THEN
              WRITE(LOGFIL,'(A)') KEYWOR
              WRITE(LOGFIL,'(9X,''V0CENT  R0CENT  A0CENT  '',
     *                          ''XK_V0C  XK_R0C  XK_A0C  '',/,
     *                       9X,6(F6.2,2X))')
     *                            VNCENT_SIGMAS,RNCENT_SIGMAS,
     *                                          ANCENT_SIGMAS,
     *                            XK_VNC_SIGMAS,XK_RNC_SIGMAS,
     *                                          XK_ANC_SIGMAS
          END IF
      END IF
C
C=======================================================================
C
      IF (KEYWOR.EQ.'SORBPROT_S') THEN
          READ(*,*) VPSORB_SIGMAS,RPSORB_SIGMAS,APSORB_SIGMAS,
     *              XK_VPS_SIGMAS,XK_RPS_SIGMAS,XK_APS_SIGMAS
          IF (IFMOCA.EQ.1) THEN
              WRITE(LOGFIL,'(A)') KEYWOR
              WRITE(LOGFIL,'(9X,''V0SORB  R0SORB  A0SORB  '',
     *                          ''XK_V0S  XK_R0S  XK_A0S  '',/,
     *                       9X,6(F6.2,2X))')
     *                     VPSORB_SIGMAS,RPSORB_SIGMAS,APSORB_SIGMAS,
     *                     XK_VPS_SIGMAS,XK_RPS_SIGMAS,XK_APS_SIGMAS
          END IF
      END IF
C
C=======================================================================
C
      IF (KEYWOR.EQ.'SORBNEUT_S') THEN
          READ(*,*) VNSORB_SIGMAS,RNSORB_SIGMAS,ANSORB_SIGMAS,
     *              XK_VNS_SIGMAS,XK_RNS_SIGMAS,XK_ANS_SIGMAS
          IF (IFMOCA.EQ.1) THEN
              WRITE(LOGFIL,'(A)') KEYWOR
              WRITE(LOGFIL,'(9X,''V0SORB  R0SORB  A0SORB  '',
     *                          ''XK_V0S  XK_R0S  XK_A0S  '',/,
     *                       9X,6(F6.2,2X))')
     *                     VNSORB_SIGMAS,RNSORB_SIGMAS,ANSORB_SIGMAS,
     *                     XK_VNS_SIGMAS,XK_RNS_SIGMAS,XK_ANS_SIGMAS
          END IF
      END IF
C
C=======================================================================
C
      IF (KEYWOR.EQ.'EFFEPROT_S') THEN
          READ(*,*) VPEFFM_SIGMAS,RPEFFM_SIGMAS,APEFFM_SIGMAS,
     *              XK_VPE_SIGMAS,XK_RPE_SIGMAS,XK_APE_SIGMAS
          IF (IFMOCA.EQ.1) THEN
              WRITE(LOGFIL,'(A)') KEYWOR
              WRITE(LOGFIL,'(9X,''V0EFFM  R0EFFM  A0EFFM  '',
     *                          ''XK_LEF  XK_REF  XK_AEF  '',/,
     *                       9X,6(F6.2,2X))')
     *                     VPEFFM_SIGMAS,RPEFFM_SIGMAS,APEFFM_SIGMAS,
     *                     XK_VPE_SIGMAS,XK_RPE_SIGMAS,XK_APE_SIGMAS
          END IF
      END IF
C
C=======================================================================
C
      IF (KEYWOR.EQ.'EFFENEUT_S') THEN
          READ(*,*) VNEFFM_SIGMAS,RNEFFM_SIGMAS,ANEFFM_SIGMAS,
     *              XK_VNE_SIGMAS,XK_RNE_SIGMAS,XK_ANE_SIGMAS
          IF (IFMOCA.EQ.1) THEN
              WRITE(LOGFIL,'(A)') KEYWOR
              WRITE(LOGFIL,'(9X,''V0EFFM  R0EFFM  A0EFFM  '',
     *                          ''XK_LEF  XK_REF  XK_AEF  '',/,
     *                       9X,6(F6.2,2X))')
     *                     VNEFFM_SIGMAS,RNEFFM_SIGMAS,ANEFFM_SIGMAS,
     *                     XK_VNE_SIGMAS,XK_RNE_SIGMAS,XK_ANE_SIGMAS
          END IF
      END IF
C
C=======================================================================
C
      IF (KEYWOR.EQ.'SORBDENS_S') THEN
          READ(*,*) ALAMPP_SIGMAS,ALAMPN_SIGMAS,ALAMNP_SIGMAS,
     *                                          ALAMNN_SIGMAS
          IF (IFMOCA.EQ.1) THEN
              WRITE(LOGFIL,'(A)') KEYWOR
              WRITE(LOGFIL,'(9X,''ALAMPP  ALAMPN  ALAMNP  ALAMNN'',/,
     *                       9X,4(F6.2,2X))')
     *                     ALAMPP_SIGMAS,ALAMPN_SIGMAS,
     *                     ALAMNP_SIGMAS,ALAMNN_SIGMAS
          END IF
      END IF
C
C=======================================================================
C
      IF (KEYWOR.EQ.'SORBTENS_S') THEN
          READ(*,*) TLAMPP_SIGMAS,TLAMPN_SIGMAS,TLAMNP_SIGMAS,
     *                                          TLAMNN_SIGMAS
          IF (IFMOCA.EQ.1) THEN
              WRITE(LOGFIL,'(A)') KEYWOR
              WRITE(LOGFIL,'(9X,''TLAMPP  TLAMPN  TLAMNP  TLAMNN'',/,
     *                       9X,4(F6.2,2X))')
     *                     TLAMPP_SIGMAS,TLAMPN_SIGMAS,
     *                     TLAMNP_SIGMAS,TLAMNN_SIGMAS
          END IF
      END IF
C
C=======================================================================
C
      IF (KEYWOR.EQ.'CENTTENS_S') THEN
          READ(*,*) CLAMPP_SIGMAS,CLAMPN_SIGMAS,CLAMNP_SIGMAS,
     *                                          CLAMNN_SIGMAS
          IF (IFMOCA.EQ.1) THEN
              WRITE(LOGFIL,'(A)') KEYWOR
              WRITE(LOGFIL,'(9X,''TLAMPP  TLAMPN  TLAMNP  TLAMNN'',/,
     *                       9X,4(F6.2,2X))')
     *                     CLAMPP_SIGMAS,CLAMPN_SIGMAS,
     *                     CLAMNP_SIGMAS,CLAMNN_SIGMAS
          END IF
      END IF
C
C=======================================================================
C
      IF (KEYWOR.EQ.'KAPACENT_S') THEN
          READ(*,*) V0CENT_SIGMAS,XK_V0C_SIGMAS,R0CENT_SIGMAS,
     *              XK_R0C_SIGMAS,A0CENT_SIGMAS,XK_A0C_SIGMAS
          IF (IFMOCA.EQ.1) THEN
              WRITE(LOGFIL,'(A)') KEYWOR
              WRITE(LOGFIL,'(9X,''V0CENT  XK_V0C  R0CENT  XK_R0C  '',
     *                          ''A0CENT  XK_A0C'',/,
     *                       9X,6(F6.2,2X))')
     *              V0CENT_SIGMAS,XK_V0C_SIGMAS,R0CENT_SIGMAS,
     *              XK_R0C_SIGMAS,A0CENT_SIGMAS,XK_A0C_SIGMAS
          END IF
      END IF
C
C=======================================================================
C
      IF (KEYWOR.EQ.'KAPASORB_S') THEN
          READ(*,*) V0SORB_SIGMAS,XK_V0S_SIGMAS,R0SORB_SIGMAS,
     *              XK_R0S_SIGMAS,A0SORB_SIGMAS,XK_A0S_SIGMAS
          IF (IFMOCA.EQ.1) THEN
              WRITE(LOGFIL,'(A)') KEYWOR
              WRITE(LOGFIL,'(9X,''V0CENT  XK_V0C  R0CENT  XK_R0C  '',
     *                          ''A0CENT  XK_A0C'',/,
     *                       9X,6(F6.2,2X))')
     *              V0SORB_SIGMAS,XK_V0S_SIGMAS,R0SORB_SIGMAS,
     *              XK_R0S_SIGMAS,A0SORB_SIGMAS,XK_A0S_SIGMAS
          END IF
      END IF
C
C=======================================================================
C=======================================================================
C     Reading the dangerous information about optional NOT
C     taking into account certain experimental levels
C=======================================================================
C=======================================================================
C 
      IF (KEYWOR.EQ.'WHATNOTAKE') THEN
C
          REMOVE_PROTON='NOT'
          REMOVE_NEUTRS='NOT'
C
          DO I=1,NDRAUS
             INDNEU_LETOUT(I)=1
             INDPRO_LETOUT(I)=1
          END DO
C
          DO I=1,NDRAUS
C
             READ (*,'(3X,I10,2X,A6,A10,I10,2X,A6)') 
     *                 INDNEU_LETOUT(I),LABNEU_REMOVE(I),SPACIN,
     *                 INDPRO_LETOUT(I),LABPRO_REMOVE(I)             
          END DO
C
          WRITE(LOGFIL,'(A)') KEYWOR
C
          DO I=1,NDRAUS
C
             WRITE(LOGFIL,'(9X,I1,2X,A6,10X,I2,2X,A6)') 
     *                      INDNEU_LETOUT(I),LABNEU_REMOVE(I),
     *                      INDPRO_LETOUT(I),LABPRO_REMOVE(I) 
          END DO
C_______________________________________________________________________
C
          NOTNEU=0
          NOTPRO=0
C
          DO I=1,NDRAUS
             IF (INDNEU_LETOUT(I).EQ.1) THEN
                 NOTNEU=NOTNEU+1
             END IF
             IF (INDPRO_LETOUT(I).EQ.1) THEN
                 NOTPRO=NOTPRO+1
             END IF
          END DO
C
C=======================================================================
C=======================================================================
C
          IF (NOTNEU.NE.0) THEN
C_______________________________________________________________________
C
              REMOVE_NEUTRS='YES'
C
              IF (ISCREN.EQ.1) THEN
                  WRITE(LSCREN,'()')
                  DO I=1,15
                     WRITE(LSCREN,'(9X,''YOU REMOVED'',I2,
     *                                 '' VALID EXPERIMENTAL '',
     *                                ''NEUTRON LEVELS FROM THE FIT'')')
     *                                                      NOTNEU
                  END DO
                  WRITE(LSCREN,'()')
              END IF
C_______________________________________________________________________
C
              WRITE(LOGFIL,'()')
C
              DO I=1,15
                 WRITE(LOGFIL,'(9X,''YOU REMOVED'',I2,
     *                             '' VALID EXPERIMENTAL '',
     *                             ''NEUTRON LEVELS FROM THE FIT'')')
     *                                                    NOTNEU
              END DO
C
              WRITE(LOGFIL,'()')
C
              WRITE(LOGFIL,'()')
              DO I=1,15
                 WRITE(LOGFIL,'(9X,''YOU REMOVED'',I2,1X,
     *                             ''VALID EXPERIMENTAL '',
     *                             ''NEUTRON LEVELS FROM THE FIT'')')
     *                                                    NOTNEU
              END DO
C
              WRITE(LOGFIL,'()')
C
              DO I=1,NDRAUS
                 IF (INDNEU_LETOUT(I).EQ.1) THEN
                     WRITE(LOGFIL,'(''Removed neutron '',a6,'' level'')
     *                                                               ')
     *                                LABNEU_REMOVE(I)
                 END IF
              END DO
C_______________________________________________________________________
C
              DO I=1,NDRAUS
C
                 IF (ISCREN.EQ.1 .AND. INDNEU_LETOUT(I).EQ.1) THEN
C
                     WRITE(LSCREN,'(''Removed neutron '',a6,'' level'')
     *                                                               ')
     *                                LABNEU_REMOVE(I)
                 END IF
              END DO
C
          END IF
C
C=======================================================================
C=======================================================================
C
          IF (NOTPRO.NE.0) THEN
C
              REMOVE_PROTON='YES'
C
              IF (ISCREN.EQ.1) THEN
C
                  WRITE(LSCREN,'()')
C
                  DO I=1,15
                     WRITE(LSCREN,'(9X,''YOU REMOVED'',I2,1X,
     *                                 ''VALID EXPERIMENTAL '',
     *                                ''PROTON  LEVELS FROM THE FIT'')')
     *                                                      NOTPRO
                  END DO
                  WRITE(LSCREN,'()')
              END IF
C
              WRITE(LOGFIL,'()')
              DO I=1,15
                 WRITE(LOGFIL,'(9X,''YOU REMOVED'',I2,1X,
     *                             ''VALID EXPERIMENTAL '',
     *                             ''PROTON  LEVELS FROM THE FIT'')')
     *                                                    NOTPRO
              END DO
C
              WRITE(LOGFIL,'()')
              WRITE(LOGFIL,'()')
C
              DO I=1,15
                 WRITE(LOGFIL,'(9X,''YOU REMOVED'',I2,1X,
     *                             ''VALID EXPERIMENTAL '',
     *                             ''PROTON  LEVELS FROM THE FIT'')')
     *                                                    NOTPRO
              END DO
C
              WRITE(LOGFIL,'()')
C
              DO I=1,NDRAUS
                 IF (INDPRO_LETOUT(I).EQ.1) THEN
                     WRITE(LSCREN,'(''Removed proton '',a6,'' level'')')
     *                     LABPRO_REMOVE(I)                            
                 END IF
              END DO
C
              DO I=1,NDRAUS
                 IF (ISCREN.EQ.1 .AND. INDPRO_LETOUT(I).EQ.1) THEN
                     WRITE(LSCREN,'(''Removed proton '',a6,'' level'')')
     *                     LABPRO_REMOVE(I)                            
                 END IF
              END DO
C
          END IF
C
      END IF 
C
C=======================================================================
C=======================================================================
C=======================================================================
C     Reading the dangerous information about optional NOT taking
C                        into account certain experimental levels
C=======================================================================
C=======================================================================
C=======================================================================
C 
      IF (KEYWOR.EQ.'IDENTIFIER') THEN
C
          READ (*,*) EXTENT_NEUTRS
          READ (*,*) EXTENT_PROTON
C
          WRITE(LOGFIL,'(A)') KEYWOR
C
          LENGTH=0
          DO I=1,40
             IF (EXTENT_NEUTRS(I:I).NE.' ') LENGTH=LENGTH+1
          END DO
C
          WRITE(LOGFIL,'(9X,''Summary of the results '',
     *                       ''goes to the files with extensions:'',/,
     *                    9X,''EXTENT_NEUTRS='',A,''.dat'')') 
     *                         EXTENT_NEUTRS(1:LENGTH)
          LENGTH=0
          DO I=1,40
             IF (EXTENT_PROTON(I:I).NE.' ') LENGTH=LENGTH+1
          END DO
C
          WRITE(LOGFIL,'(9X,''EXTENT_PROTON='',A,''.dat'')') 
     *                        EXTENT_PROTON(1:LENGTH)
C
          LENGTH=0
          DO I=1,40
             IF (EXTENT_NEUTRS(I:I).NE.' ') LENGTH=LENGTH+1
          END DO
C
          WRITE(LSCREN,'()')
          WRITE(LSCREN,'(''Summary of the results '',
     *                   ''goes to the files with extensions:'',/,
     *                   ''EXTENT_NEUTRS='',A,''.dat'')') 
     *                     EXTENT_NEUTRS(1:LENGTH)
          LENGTH=0
          DO I=1,40
             IF (EXTENT_PROTON(I:I).NE.' ') LENGTH=LENGTH+1
          END DO
          WRITE(LSCREN,'(''EXTENT_PROTON='',A,''.dat'')') 
     *                     EXTENT_PROTON(1:LENGTH)
C
      END IF
C
C=======================================================================
C=======================================================================
C=======================================================================
C     Reading the length of the energy levels (plotting reasons)
C=======================================================================
C=======================================================================
C=======================================================================
C 
      IF (KEYWOR.EQ.'LEV_LENGTH') THEN
C
          READ(*,*) XMIN_T,XMAX_T,XMIN_E,XMAX_E
C
          WRITE(LOGFIL,'(A)') KEYWOR
C
          WRITE(LOGFIL,'(9X,''XMIN_T  XMAX_T  XMIN_E  XMAX_E'')')
          WRITE(LOGFIL,'(9X,4F6.2)')XMIN_T,XMAX_T,XMIN_E,XMAX_E
C
      END IF
C
C=======================================================================
C      
      IF (KEYWOR.EQ.'ATTRIBUTES') THEN
C
          READ(*,*) ISTYLE,I_TYPE
          READ(*,*) STRING
C
          READ(*,*) (ICOLOR(I),I=1,NDCOLO)
          READ(*,*) STRING
C
          READ(*,*) (ITHICK(I),I=1,NDCOLO)
C
          WRITE(LOGFIL,'(A)') KEYWOR
C
          WRITE(LOGFIL,'(9X,''STYLE1  TYPEPN'')')
C
          WRITE(LOGFIL,'(9X,I6,2X,I6)') ISTYLE,I_TYPE
C
          WRITE(LOGFIL,'(9X,<NDCOLO>(''COLO'',I2.2,2X))') (I,I=1,NDCOLO)
          WRITE(LOGFIL,'(9X,<NDCOLO>(I6,2X))')    (ICOLOR(I),I=1,NDCOLO)
          WRITE(LOGFIL,'(9X,<NDCOLO>(''ITHI'',I2.2,2X))') (I,I=1,NDCOLO)
          WRITE(LOGFIL,'(9X,<NDCOLO>(I6,2X))')    (ITHICK(I),I=1,NDCOLO)
C
      END IF
C
C=======================================================================
C=======================================================================
C                         Ready to run ?
C=======================================================================
C=======================================================================
C @@@ WE WILL NEED INTELLIGENT PRINT OF THIS INFORMATION - LET'S DISCUSS
      IF (KEYWOR.EQ.'EXECUTE!!!') THEN
C
          PARPOT( 1)=V0CENT_PROTON
          PARPOT( 2)=R0CENT_PROTON
          PARPOT( 3)=A0CENT_PROTON
          PARPOT( 4)=V0SORB_PROTON
          PARPOT( 5)=R0SORB_PROTON
          PARPOT( 6)=A0SORB_PROTON
          PARPOT( 7)=V0EFFM_PROTON
          PARPOT( 8)=R0EFFM_PROTON
          PARPOT( 9)=A0EFFM_PROTON
          PARPOT(10)=R0COUL
     
          PARPOT(11)=XK_V0C_PROTON
          PARPOT(12)=XK_R0C_PROTON
          PARPOT(13)=XK_A0C_PROTON
          PARPOT(14)=XK_LAM_PROTON
          PARPOT(15)=XK_RSO_PROTON
          PARPOT(16)=XK_ASO_PROTON
          PARPOT(17)=XK_LEF_PROTON
          PARPOT(18)=XK_REF_PROTON
          PARPOT(19)=XK_AEF_PROTON
          PARPOT(20)=XK_COU
     
          PARPOT(21)=V0CENT_NEUTRS
          PARPOT(22)=R0CENT_NEUTRS
          PARPOT(23)=A0CENT_NEUTRS
          PARPOT(24)=V0SORB_NEUTRS
          PARPOT(25)=R0SORB_NEUTRS
          PARPOT(26)=A0SORB_NEUTRS
          PARPOT(27)=V0EFFM_NEUTRS
          PARPOT(28)=R0EFFM_NEUTRS
          PARPOT(29)=A0EFFM_NEUTRS
     
          PARPOT(30)=XK_V0C_NEUTRS
          PARPOT(31)=XK_R0C_NEUTRS
          PARPOT(32)=XK_A0C_NEUTRS
          PARPOT(33)=XK_LAM_NEUTRS
          PARPOT(34)=XK_RSO_NEUTRS
          PARPOT(35)=XK_ASO_NEUTRS
          PARPOT(36)=XK_LEF_NEUTRS
          PARPOT(37)=XK_REF_NEUTRS
          PARPOT(38)=XK_AEF_NEUTRS
     
          PARPOT(39)=ALAMPP
          PARPOT(40)=ALAMPN
          PARPOT(41)=ALAMNP
          PARPOT(42)=ALAMNN

          PARPOT(43)=TLAMPP
          PARPOT(44)=TLAMPN
          PARPOT(45)=TLAMNP
          PARPOT(46)=TLAMNN
          
          PARPOT(47)=CLAMPP
          PARPOT(48)=CLAMPN
          PARPOT(49)=CLAMNP
          PARPOT(50)=CLAMNN
C
          PARPOT(51)=V0CENT_KAPPAR
          PARPOT(52)=XK_V0C_KAPPAR
          PARPOT(53)=R0CENT_KAPPAR
          PARPOT(54)=XK_R0C_KAPPAR
          PARPOT(55)=A0CENT_KAPPAR
          PARPOT(56)=XK_A0C_KAPPAR
C
          PARPOT(57)=V0SORB_KAPPAR
          PARPOT(58)=XK_V0S_KAPPAR
          PARPOT(59)=R0SORB_KAPPAR
          PARPOT(60)=XK_R0S_KAPPAR
          PARPOT(61)=A0SORB_KAPPAR
          PARPOT(62)=XK_A0S_KAPPAR
C
C=======================================================================
C
          PARPOT_XMEANS( 1)=VPCENT_XMEANS
          PARPOT_XMEANS( 2)=RPCENT_XMEANS
          PARPOT_XMEANS( 3)=APCENT_XMEANS
          PARPOT_XMEANS( 4)=VPSORB_XMEANS
          PARPOT_XMEANS( 5)=RPSORB_XMEANS
          PARPOT_XMEANS( 6)=APSORB_XMEANS
          PARPOT_XMEANS( 7)=VPEFFM_XMEANS
          PARPOT_XMEANS( 8)=RPEFFM_XMEANS
          PARPOT_XMEANS( 9)=APEFFM_XMEANS
          PARPOT_XMEANS(10)=R0COUL_XMEANS
     
          PARPOT_XMEANS(11)=XK_VPC_XMEANS
          PARPOT_XMEANS(12)=XK_RPC_XMEANS
          PARPOT_XMEANS(13)=XK_APC_XMEANS
          PARPOT_XMEANS(14)=XK_VPS_XMEANS
          PARPOT_XMEANS(15)=XK_RPS_XMEANS
          PARPOT_XMEANS(16)=XK_APS_XMEANS
          PARPOT_XMEANS(17)=XK_VPE_XMEANS
          PARPOT_XMEANS(18)=XK_RPE_XMEANS
          PARPOT_XMEANS(19)=XK_APE_XMEANS
          PARPOT_XMEANS(20)=XK_COU_XMEANS
     
          PARPOT_XMEANS(21)=VNCENT_XMEANS
          PARPOT_XMEANS(22)=RNCENT_XMEANS
          PARPOT_XMEANS(23)=ANCENT_XMEANS
          PARPOT_XMEANS(24)=VNSORB_XMEANS
          PARPOT_XMEANS(25)=RNSORB_XMEANS
          PARPOT_XMEANS(26)=ANSORB_XMEANS
          PARPOT_XMEANS(27)=VNEFFM_XMEANS
          PARPOT_XMEANS(28)=RNEFFM_XMEANS
          PARPOT_XMEANS(29)=ANEFFM_XMEANS
     
          PARPOT_XMEANS(30)=XK_VNC_XMEANS
          PARPOT_XMEANS(31)=XK_RNC_XMEANS
          PARPOT_XMEANS(32)=XK_ANC_XMEANS
          PARPOT_XMEANS(33)=XK_VNS_XMEANS
          PARPOT_XMEANS(34)=XK_RNS_XMEANS
          PARPOT_XMEANS(35)=XK_ANS_XMEANS
          PARPOT_XMEANS(36)=XK_VNE_XMEANS
          PARPOT_XMEANS(37)=XK_RNE_XMEANS
          PARPOT_XMEANS(38)=XK_ANE_XMEANS
     
          PARPOT_XMEANS(39)=ALAMPP_XMEANS
          PARPOT_XMEANS(40)=ALAMPN_XMEANS
          PARPOT_XMEANS(41)=ALAMNP_XMEANS
          PARPOT_XMEANS(42)=ALAMNN_XMEANS

          PARPOT_XMEANS(43)=TLAMPP_XMEANS
          PARPOT_XMEANS(44)=TLAMPN_XMEANS
          PARPOT_XMEANS(45)=TLAMNP_XMEANS
          PARPOT_XMEANS(46)=TLAMNN_XMEANS
          
          PARPOT_XMEANS(47)=CLAMPP_XMEANS
          PARPOT_XMEANS(48)=CLAMPN_XMEANS
          PARPOT_XMEANS(49)=CLAMNP_XMEANS
          PARPOT_XMEANS(50)=CLAMNN_XMEANS
C
          PARPOT_XMEANS(51)=V0CENT_XMEANS
          PARPOT_XMEANS(52)=XK_V0C_XMEANS
          PARPOT_XMEANS(53)=R0CENT_XMEANS
          PARPOT_XMEANS(54)=XK_R0C_XMEANS
          PARPOT_XMEANS(55)=A0CENT_XMEANS
          PARPOT_XMEANS(56)=XK_A0C_XMEANS
C
          PARPOT_XMEANS(57)=V0SORB_XMEANS
          PARPOT_XMEANS(58)=XK_V0S_XMEANS
          PARPOT_XMEANS(59)=R0SORB_XMEANS
          PARPOT_XMEANS(60)=XK_R0S_XMEANS
          PARPOT_XMEANS(61)=A0SORB_XMEANS
          PARPOT_XMEANS(62)=XK_A0S_XMEANS
C
C=======================================================================
C
          PARPOT_SIGMAS( 1)=VPCENT_SIGMAS
          PARPOT_SIGMAS( 2)=RPCENT_SIGMAS
          PARPOT_SIGMAS( 3)=APCENT_SIGMAS
          PARPOT_SIGMAS( 4)=VPSORB_SIGMAS
          PARPOT_SIGMAS( 5)=RPSORB_SIGMAS
          PARPOT_SIGMAS( 6)=APSORB_SIGMAS
          PARPOT_SIGMAS( 7)=VPEFFM_SIGMAS
          PARPOT_SIGMAS( 8)=RPEFFM_SIGMAS
          PARPOT_SIGMAS( 9)=APEFFM_SIGMAS
          PARPOT_SIGMAS(10)=R0COUL_SIGMAS
     
          PARPOT_SIGMAS(11)=XK_VPC_SIGMAS
          PARPOT_SIGMAS(12)=XK_RPC_SIGMAS
          PARPOT_SIGMAS(13)=XK_APC_SIGMAS
          PARPOT_SIGMAS(14)=XK_VPS_SIGMAS
          PARPOT_SIGMAS(15)=XK_RPS_SIGMAS
          PARPOT_SIGMAS(16)=XK_APS_SIGMAS
          PARPOT_SIGMAS(17)=XK_VPE_SIGMAS
          PARPOT_SIGMAS(18)=XK_RPE_SIGMAS
          PARPOT_SIGMAS(19)=XK_APE_SIGMAS
          PARPOT_SIGMAS(20)=XK_COU_SIGMAS
     
          PARPOT_SIGMAS(21)=VNCENT_SIGMAS
          PARPOT_SIGMAS(22)=RNCENT_SIGMAS
          PARPOT_SIGMAS(23)=ANCENT_SIGMAS
          PARPOT_SIGMAS(24)=VNSORB_SIGMAS
          PARPOT_SIGMAS(25)=RNSORB_SIGMAS
          PARPOT_SIGMAS(26)=ANSORB_SIGMAS
          PARPOT_SIGMAS(27)=VNEFFM_SIGMAS
          PARPOT_SIGMAS(28)=RNEFFM_SIGMAS
          PARPOT_SIGMAS(29)=ANEFFM_SIGMAS
     
          PARPOT_SIGMAS(30)=XK_VNC_SIGMAS
          PARPOT_SIGMAS(31)=XK_RNC_SIGMAS
          PARPOT_SIGMAS(32)=XK_ANC_SIGMAS
          PARPOT_SIGMAS(33)=XK_VNS_SIGMAS
          PARPOT_SIGMAS(34)=XK_RNS_SIGMAS
          PARPOT_SIGMAS(35)=XK_ANS_SIGMAS
          PARPOT_SIGMAS(36)=XK_VNE_SIGMAS
          PARPOT_SIGMAS(37)=XK_RNE_SIGMAS
          PARPOT_SIGMAS(38)=XK_ANE_SIGMAS
     
          PARPOT_SIGMAS(39)=ALAMPP_SIGMAS
          PARPOT_SIGMAS(40)=ALAMPN_SIGMAS
          PARPOT_SIGMAS(41)=ALAMNP_SIGMAS
          PARPOT_SIGMAS(42)=ALAMNN_SIGMAS

          PARPOT_SIGMAS(43)=TLAMPP_SIGMAS
          PARPOT_SIGMAS(44)=TLAMPN_SIGMAS
          PARPOT_SIGMAS(45)=TLAMNP_SIGMAS
          PARPOT_SIGMAS(46)=TLAMNN_SIGMAS
          
          PARPOT_SIGMAS(47)=CLAMPP_SIGMAS
          PARPOT_SIGMAS(48)=CLAMPN_SIGMAS
          PARPOT_SIGMAS(49)=CLAMNP_SIGMAS
          PARPOT_SIGMAS(50)=CLAMNN_SIGMAS
C
          PARPOT_SIGMAS(51)=V0CENT_SIGMAS
          PARPOT_SIGMAS(52)=XK_V0C_SIGMAS
          PARPOT_SIGMAS(53)=R0CENT_SIGMAS
          PARPOT_SIGMAS(54)=XK_R0C_SIGMAS
          PARPOT_SIGMAS(55)=A0CENT_SIGMAS
          PARPOT_SIGMAS(56)=XK_A0C_SIGMAS
C
          PARPOT_SIGMAS(57)=V0SORB_SIGMAS
          PARPOT_SIGMAS(58)=XK_V0S_SIGMAS
          PARPOT_SIGMAS(59)=R0SORB_SIGMAS
          PARPOT_SIGMAS(60)=XK_R0S_SIGMAS
          PARPOT_SIGMAS(61)=A0SORB_SIGMAS
          PARPOT_SIGMAS(62)=XK_A0S_SIGMAS
C
C=======================================================================
C
          IF (LOGWRI.GT.0) THEN
              WRITE(LOGFIL,'(/,''Exiting  NAMELI'')')
          END IF
C
          RETURN
C
      END IF
C
C=======================================================================
C                           Scanning continues
C=======================================================================
C
      GO TO  1
C
C=======================================================================
C
   2  CONTINUE
C
C=======================================================================
C     Illegal end of file in the input data stream
C=======================================================================
C
      WRITE(NOUTPT,'(//,''Something wrong with the structure of the '',
     *             ''input file - perhaps forgotten <EXECUTESET> ?'')')
C
C=======================================================================
C                    Forgotten 'EXECUTESET' statement ? ?
C=======================================================================
C
      WRITE(NOUTPT,'(''No "EXECUTE!!!" in the input data stream '',
     *               ''(in NAMELI)'')')
C   
      STOP
     *
     *   'Wrong input structure, eof not allowed here, STOP from NAMELI'
C
      END
C
C=======================================================================
C=======================================================================
C=======================================================================

C FILE NAME = wspher05_hamilt_15.f ! Keep this symbol:    $ident@string$
C      
C=======================================================================
C=======================================================================
C          CONSTRUCTION OF THE HAMILTONIAN MATRIX AND RELATED
C=======================================================================
C=======================================================================
C
      SUBROUTINE WS_UNI(INUCLI,CHISQU_CHOICE,IABORT,I_MODE,I_FLAG,
     *                                CHISQU_PROTON,CHISQU_NEUTRS)
C
      INCLUDE   'MATDIM/NDSPEC.f'
      INCLUDE   'MATDIM/NDGAUS.f'
      INCLUDE   'MATDIM/NDIM_N.f'
      INCLUDE   'MATDIM/NDNUCL.f'
      INCLUDE   'MATDIM/NDLEXP.f'
      INCLUDE   'MATDIM/NDMAIN.f'
C      
      PARAMETER
     *         (NDIM_L=NDIM_N,NDBASE=NDIM_N+1)
C
      CHARACTER
     *          WHATEX*6,TYPCHI*6,INPSYM*6,NUCSYM*6
      CHARACTER
     *          LABTHE_PROTON*6,LABTHE_NEUTRS*6
      CHARACTER
     *          LABEXP_PROTON*6,LABEXP_NEUTRS*6
C
      DIMENSION
     *          LABORD_AUXILP(0:NDIM_N,0:NDIM_L,0:2*NDIM_L+1),
     *          LABORD_AUXILN(0:NDIM_N,0:NDIM_L,0:2*NDIM_L+1)
      DIMENSION
     *          RFUNUP_PROTON(1:NDGAUS,0:NDIM_N,0:NDIM_L),
     *          RFUNDN_PROTON(1:NDGAUS,0:NDIM_N,0:NDIM_L),
     *          RFUNUP_NEUTRS(1:NDGAUS,0:NDIM_N,0:NDIM_L),
     *          RFUNDN_NEUTRS(1:NDGAUS,0:NDIM_N,0:NDIM_L)
      DIMENSION
     *          RAD_UP(0:NDIM_L,0:NDIM_N),
     *         	RAD_DN(0:NDIM_L,0:NDIM_N)
      DIMENSION
     *          HSORUP(1:NDBASE,1:NDBASE),
     *          HSORDN(1:NDBASE,1:NDBASE)
      DIMENSION
     *          SPENUP(1:NDBASE),
     *          SPENDN(1:NDBASE),     
     *          AUXVEC(1:NDBASE)  
C
      COMMON
     *       /LEVORD/ LABORD_PROTON(1:NDNUCL,
     *                              0:NDIM_N,0:NDIM_L,0:2*NDIM_L+1),
     *                LABORD_NEUTRS(1:NDNUCL,
     *                              0:NDIM_N,0:NDIM_L,0:2*NDIM_L+1)
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
     *       /DIMSPH/ LDSING_PROTON,LEVTHE_PROTON,
     *                              LEVEXP_PROTON(1:NDNUCL),
     *                LDSING_NEUTRS,LEVTHE_NEUTRS,
     *                              LEVEXP_NEUTRS(1:NDNUCL)
      COMMON
     *       /NODWEI/ XRNODE(1:NDGAUS,0:NDIM_L),
     *                XRWEIG(1:NDGAUS,0:NDIM_L)
      COMMON
     *       /POLPOL/ XPOLYN(1:NDGAUS,0:NDIM_N,0:NDIM_L), 
     *                XDPOLY(1:NDGAUS,0:NDIM_N,0:NDIM_L) 
      COMMON
     *       /V2COEF/ VCOEUP_PROTON(0:NDIM_L,1:NDBASE),
     *                VCOEDN_PROTON(0:NDIM_L,1:NDBASE),
     *
     *                VCOEUP_NEUTRS(0:NDIM_L,1:NDBASE),
     *                VCOEDN_NEUTRS(0:NDIM_L,1:NDBASE)
      COMMON
     *       /ENUPDN/ ENERUP_PROTON(0:NDIM_L,1:NDBASE),
     *                ENERDN_PROTON(0:NDIM_L,1:NDBASE),

     *
     *                ENERUP_NEUTRS(0:NDIM_L,1:NDBASE),
     *                ENERDN_NEUTRS(0:NDIM_L,1:NDBASE)
      COMMON
     *       /MODIFR/ RMSEXP_PROTON(1:NDNUCL),
     *                RMSEXP_NEUTRS(1:NDNUCL),
     *                RMSTHE_PROTON,RMSTHE_NEUTRS,
     *                RMODIF_PROTON,RMODIF_NEUTRS,
     *                              PUSHIN,DELTAR
      COMMON
     *       /FERGPS/ FERMEX_PROTON(1:NDNUCL),GAPEXP_PROTON(1:NDNUCL),
     *                DENSUP_PROTON(1:NDNUCL),DENSDW_PROTON(1:NDNUCL),
     *
     *                FERMEX_NEUTRS(1:NDNUCL),GAPEXP_NEUTRS(1:NDNUCL),
     *                DENSUP_NEUTRS(1:NDNUCL),DENSDW_NEUTRS(1:NDNUCL) 
      COMMON
     *       /QUANTN/ NWSSPH_PROTON(1:NDSPEC),
     *                LWSSPH_PROTON(1:NDSPEC),
     *                JWSSPH_PROTON(1:NDSPEC),
     *
     *                NWSSPH_NEUTRS(1:NDSPEC),
     *                LWSSPH_NEUTRS(1:NDSPEC),
     *                JWSSPH_NEUTRS(1:NDSPEC)
      COMMON
     *       /SPECTR/ ENETHE_PROTON(1:NDSPEC),
     *                ENETHE_NEUTRS(1:NDSPEC)
      COMMON
     *       /THELAB/ LABTHE_PROTON(1:NDSPEC),
     *                LABTHE_NEUTRS(1:NDSPEC)
      COMMON
     *       /EX_ENE/ EXPEXP_PROTON(1:NDNUCL,1:NDLEXP),
     *                IDEGEX_PROTON(1:NDNUCL,1:NDLEXP),
     *
     *                EXPEXP_NEUTRS(1:NDNUCL,1:NDLEXP),
     *                IDEGEX_NEUTRS(1:NDNUCL,1:NDLEXP)
      COMMON
     *       /EX_LAB/ LABEXP_PROTON(1:NDNUCL,1:NDLEXP),
     *                LABEXP_NEUTRS(1:NDNUCL,1:NDLEXP)
C     
      COMMON
     *       /GAUSHE/ NGAUSS,NSHELL_PROTON,NSHELL_NEUTRS
      COMMON
     *       /DATA_V/ AOSCIL_PROTON(1:NDNUCL),
     *                AOSCIL_NEUTRS(1:NDNUCL),
     *                HOMEGA(1:NDNUCL)   
      COMMON
     *       /NUCLEU/ IZ_FIX,IN_FIX,ISOSPI   
      COMMON
     *       /NUCINF/ LDNUCL,
     *                NUCSYM(1:NDNUCL),
     *                NUMB_Z(1:NDNUCL),
     *                NUMB_N(1:NDNUCL)
      COMMON
     *       /MINITE/ IDEFCN,ITECHI
      COMMON
     *       /TESFIT/ IFTEST,IFFITS,IFMESH,ISPE_P,IFMOCA
      COMMON
     *       /DENSIF/ ITRNUM
     *       /WHENST/ DIFFVA
     *       /DENSTR/ IFDENS
     *       /GNUPLO/ IGNUPL
      COMMON
     *       /ACTIVE/ IACTIV,IwMODE
     *       /FEVALS/ IFUNCT_EVALUS
     *       /MASSIV/ IMASIV_PRODUC
      COMMON
     *       /PRFILE/ ISCREN,LOGFIL,LOGCHP,LOGCHN,LOGHIP,LOGHIN,
     *                LOGPRO,LOGNEU,LOGWRI,LSCREN,IRESUL,ICONVE,
     *                NOUTPT,LOGAUX,LOGENE,LOGRAD,LOGGRD,LOGMSH,
     *                                                   LOGBIS
C  @@@ IRENE - WHAT IS THE ROLE OF THIS EPSSIL ???
      DATA 
     *        EPSSIL / 1.0e-3 /
      DATA
     *        SPENER_MXPRIN / 10.000 / ! No print above that energy
C
C=======================================================================
C
C     This subroutine functions as the standard WS-ODD code. Its form
C     has been adapted to enable the connection with the minimisation
C     subroutines.
C
C     Some extra routines, allowing a comparison with experiment, and 
C     facilitating the definition of the ch^2 type routines are:
C
C     RDPROT_RMSRAD - reads the experimental  proton  mean square radii;
C     PROVID_RMS_NR - gives the experimental  neutron mean square radii;
C     READEX_LEVELS - reads the experimental spherical levels for doubly
C                     magic nuclei;
C
C     PREPAR_EXPLEV - transforms single particle excitation energies of
C                     one-hole and one-particle nuclei,  which are read 
C                     in SUBROUTINE READEX_LEVELS, into format suitable
C                     for comparisons with the calculations;
C
C     SERVIS - internal control subroutine
C     EXPTHE - internal control subroutine
C
C=======================================================================
C
      CALL CPUTIM('WS_UNI',1)
C
C=======================================================================
C
      IwMODE=I_MODE
C
      IF (I_MODE.EQ.0.OR.I_MODE.EQ.2) THEN
          IF (LOGWRI.GT.4) WRITE(IRESUL,'(/,''IwMODE='',I1)') IwMODE
      ELSE
          IF (LOGWRI.GT.4) WRITE(IRESUL,'(  ''IwMODE='',I1)') IwMODE
      END IF
C      
C     IFEXPE=IFFITS
C
C @@@ IRENE - what is is nonsense? NOUTPT is Standard output = 6 always
C      NOUTPT=IRESUL
      IABORT=0
C
C      IDEFCN=IDEFCN+1
C
C=======================================================================
C
      IF (IDEFCN.GT.ITECHI) THEN
C
CID          WRITE(ICONVE,'(/,80(''!''),/,''!WS_UNI'',72x,''!''/,
C     *                 ''!   The allowed maximum number of function'',
C     *                 '' evaluations ITECHI='',i4,''  exceeded   !'',/,
C     *                 ''!'',78x,''!''/,80(''!''))')
C     *                                                           ITECHI
C          
          WRITE(LOGAUX,'(/,80(''!''),/,''!WS_UNI'',72x,''!''/,
     *                 ''!   The allowed maximum number of function'',
     *                 '' evaluations ITECHI='',i4,''  exceeded   !'',/,
     *                 ''!'',78x,''!''/,80(''!''))')
     *                                                           ITECHI
          IABORT=1
C          
CID          WRITE(ICONVE,'(''RETURN from WS_UNI - before WSSTAN'')')
C          
          WRITE(LOGAUX,'(''RETURN from WS_UNI - before WSSTAN'')')
C          
          RETURN
C
      END IF
C
C=======================================================================
C     Solving the Schroedinger Equation (params in POTPOT)
C=======================================================================
C
      IZ_FIX=NUMB_Z(INUCLI)
      IN_FIX=NUMB_N(INUCLI)
C
      IF (IFDENS.EQ.0) THEN
C
          IF (LOGWRI.GT.4) THEN
              WRITE(LOGFIL,'(15X,''Entering WSSTAN from WS_UNI [1], '',
     *                           ''IFDENS=0 IZ_FIX='',I2,
     *                                   '' IN_FIX='',I3)') 
     *                                      IZ_FIX,IN_FIX
          END IF 
C
          CALL WSSTAN(INUCLI,
     *                XPOLYN,XRNODE,XRWEIG,NGAUSS,NSHELL_PROTON,
     *                NDGAUS,NDIM_N,NDIM_L,NDBASE,NSHELL_NEUTRS,
     *                       HSORUP,HSORDN,SPENUP,SPENDN,AUXVEC,IABORT,
     *                                     CMATUP_PROTON,CMATDN_PROTON,
     *                                     CMATUP_NEUTRS,CMATDN_NEUTRS,
     *                                     ENERUP_PROTON,ENERDN_PROTON,
     *                                     ENERUP_NEUTRS,ENERDN_NEUTRS,
     *                                     ENEKIN_PROTON,ENEKIN_NEUTRS)
C
      END IF
C
C=======================================================================
C      
      DO N=0,NDIM_N
         DO L=0,NDIM_L
            DO J=0,2*NDIM_L+1
               LABORD_AUXILP(N,L,J)=LABORD_PROTON(INUCLI,N,L,J)
               LABORD_AUXILN(N,L,J)=LABORD_NEUTRS(INUCLI,N,L,J)
            END DO
         END DO
      END DO
C
      RMSTHE_PROTON=0.
      RMSTHE_NEUTRS=0.
C
      ISOSPI=1      
C
      IF (LOGWRI.GT.5) THEN
          WRITE(LOGFIL,'(/,15X,''Entering NUCRAD from WS_UNI'',1X,
     *                                             ''ISOSPI='',I1)')
     *                                               ISOSPI
      END IF  
C
                         LDSPEC=LDSING_PROTON
      CALL NUCRAD(NDSPEC,LDSPEC,NDGAUS,NGAUSS,NDBASE,XPOLYN,
     *            NDIM_L,NDIM_N,RAD_UP,RAD_DN,XRNODE,XRWEIG,
     *            AOSCIL_PROTON,RFUNUP_PROTON,RFUNDN_PROTON,
     *            NSHELL_PROTON,ENERUP_PROTON,ENERDN_PROTON,
     *            NWSSPH_PROTON,LWSSPH_PROTON,JWSSPH_PROTON,
     *            RMSTHE_PROTON,RMSPAI_PROTON,RMSNOP_PROTON,
     *            LABORD_AUXILP,CMATUP_PROTON,CMATDN_PROTON,  
     *            ENETHE_PROTON,VCOEUP_PROTON,VCOEDN_PROTON,
     *                                 ISOSPI,IZ_FIX,INUCLI)
C
      ISOSPI=0      
C
      IF (LOGWRI.GT.5) THEN
          WRITE(LOGFIL,'(/,15X,''Entering NUCRAD from WS_UNI'',1X,
     *                                             ''ISOSPI='',I1)')
     *                                               ISOSPI
      END IF  
                         LDSPEC=LDSING_NEUTRS
      CALL NUCRAD(NDSPEC,LDSPEC,NDGAUS,NGAUSS,NDBASE,XPOLYN,
     *            NDIM_L,NDIM_N,RAD_UP,RAD_DN,XRNODE,XRWEIG,
     *            AOSCIL_NEUTRS,RFUNUP_NEUTRS,RFUNDN_NEUTRS,
     *            NSHELL_NEUTRS,ENERUP_NEUTRS,ENERDN_NEUTRS,
     *            NWSSPH_NEUTRS,LWSSPH_NEUTRS,JWSSPH_NEUTRS,
     *            RMSTHE_NEUTRS,RMSPAI_NEUTRS,RMSNOP_NEUTRS,
     *            LABORD_AUXILN,CMATUP_NEUTRS,CMATDN_NEUTRS,  
     *            ENETHE_NEUTRS,VCOEUP_NEUTRS,VCOEDN_NEUTRS,
     *                                 ISOSPI,IN_FIX,INUCLI)
C
C=======================================================================
C     Testing the quality of the experimental spectra for
C     doubly-magic, spherical nuclei           (IFEXPE=1)
C=======================================================================
C
      ISOSPI=1
C
      IF (LOGWRI.GT.5) THEN
          WRITE(LOGFIL,'(/,15X,''Entering EXPTHE from WS_UNI'',1X,
     *                                             ''ISOSPI='',I1)')
     *                                               ISOSPI
      END IF  
C      
      CALL EXPTHE(IZ_FIX,IN_FIX,DENSUP_PROTON,DENSDW_PROTON,
     *            FERMEX_PROTON,GAPEXP_PROTON,IDEFCN,ISOSPI,
     *                   IPARAM,CHISQU_CHOICE,I_FLAG,INUCLI)
C
      CHISQU_PROTON=CHISQU_CHOICE
C
      ISOSPI=0
C
      IF (LOGWRI.GT.5) THEN
          WRITE(LOGFIL,'(/,15X,''Entering EXPTHE from WS_UNI'',1X,
     *                                             ''ISOSPI='',I1)')
     *                                               ISOSPI
      END IF  
C
      CALL EXPTHE(IZ_FIX,IN_FIX,DENSUP_NEUTRS,DENSDW_NEUTRS,
     *            FERMEX_NEUTRS,GAPEXP_NEUTRS,IDEFCN,ISOSPI,
     *                   IPARAM,CHISQU_CHOICE,I_FLAG,INUCLI)
C
      CHISQU_NEUTRS=CHISQU_CHOICE
C 
      IFUNCT_EVALUS=IFUNCT_EVALUS+1
C
C=======================================================================     
C
      WRITE(NOUTPT,'()')
C
      WRITE(NOUTPT,'(80(''#''),/,''#'',78X,''#'')')
      WRITE(NOUTPT,'(''#'',20X,''IZ_FIX='',I3,20X,
     *                                    ''IN_FIX='',I3,
     *                                   18X,''#'',/,''#'',78X,''#'')')
     *                                 IZ_FIX,            IN_FIX
       
      WRITE(NOUTPT,'(80(''#''),/,''#'',78X,''#'')')
      WRITE(NOUTPT,'(''#'',2X,''NEUTRONS: Theoretical vs. '',
     *                            ''Experimental Energies'',29x,''#'')')
      WRITE(NOUTPT,'(''#'',78X,''#'',/,80(''#''))')
      WRITE(NOUTPT,'(''#'',78X,''#'')')
       
      WRITE(NOUTPT,'(''#'',3X,''No)'',7x,''Labels'',7x,
     *                                              ''EneThe'',
     *                      14X,''No)'',7x,''Labels'',7x,''EneExp'',3X,
     *                                                         ''#'')')
C
      DO ITHEOR=1,LDSING_NEUTRS
C
         IF (ENETHE_NEUTRS(ITHEOR+1).GT.SPENER_MXPRIN) THEN
             GO TO 7
         END IF
C
         DO IEXPER=1,LEVEXP_NEUTRS(INUCLI)
C
            IF (LABEXP_NEUTRS(INUCLI,IEXPER) .EQ.
     *          LABTHE_NEUTRS(ITHEOR))       THEN
C
                WRITE(NOUTPT,'(''#'',i5,'')'',7x,a6,f20.13,
     *                              4X,i5,'')'',7x,a6,f15.4,1X,''#'')')
     *                                  ITHEOR,
     *                                  LABTHE_NEUTRS(ITHEOR),
     *                                  ENETHE_NEUTRS(ITHEOR),
     *                                  IEXPER,
     *                                  LABEXP_NEUTRS(INUCLI,IEXPER),
     *                                  EXPEXP_NEUTRS(INUCLI,IEXPER)
                 GO TO 5
            END IF
         END DO
           
         WRITE(NOUTPT,'(''#'',i5,'')'',7x,a6,f20.13,39X,    
     *                            ''#'')')ITHEOR,LABTHE_NEUTRS(ITHEOR),
     *                                           ENETHE_NEUTRS(ITHEOR)
  5      CONTINUE
C
      END DO
C
  7   CONTINUE
C
      WRITE(NOUTPT,'(''#'',78X,''#'',/,80(''#''),/,''#'',78X,''#'')')
C
      WRITE(NOUTPT,'(''#'',2X,''PROTONS: Theoretical vs. '',
     *                            ''Experimental Energies'',30x,''#'')')
      WRITE(NOUTPT,'(''#'',78X,''#'',/,80(''#''))')
      WRITE(NOUTPT,'(''#'',78X,''#'')')
      
      WRITE(NOUTPT,'(''#'',3X,''No)'',7x,''Labels'',7x,
     *                                              ''EneThe'',
     *                      14X,''No)'',7x,''Labels'',7x,''EneExp'',3X,
     *                                                         ''#'')')
      DO ITHEOR=1,LDSING_PROTON
C
         IF (ENETHE_NEUTRS(ITHEOR+1).GT.SPENER_MXPRIN) THEN
             GO TO 8
         END IF
C
         DO IEXPER=1,LEVEXP_PROTON(INUCLI)
C
            IF (LABEXP_PROTON(INUCLI,IEXPER) .EQ.
     *          LABTHE_PROTON(ITHEOR))       THEN
C
                WRITE(NOUTPT,'(''#'',i5,'')'',7x,a6,f20.13,
     *                           4X,i5,'')'',7x,a6,f15.4,1X,''#'')')
     *                                  ITHEOR,
     *                                  LABTHE_PROTON(ITHEOR),
     *                                  ENETHE_PROTON(ITHEOR),
     *                                  IEXPER,
     *                                  LABEXP_PROTON(INUCLI,IEXPER),
     *                                  EXPEXP_PROTON(INUCLI,IEXPER)
                GO TO 6
            END IF
         END DO
           
         WRITE(NOUTPT,'(''#'',i5,'')'',7x,a6,f20.13,39X, 
     *                            ''#'')')ITHEOR,LABTHE_PROTON(ITHEOR),
     *                                           ENETHE_PROTON(ITHEOR)
  6      CONTINUE
C
      END DO
C
  8   CONTINUE
C
      WRITE(NOUTPT,'(''#'',78X,''#'',/,80(''#''))')
C
C=======================================================================
C
      IF (LOGWRI.GT.4) THEN
          WRITE(LOGFIL,'(12X,''Exiting  WS_UNI'')')
      END IF
C      
C=======================================================================
C
      CALL CPUTIM('WS_UNI',0)   
C      
C=======================================================================
C      
      RETURN
      END
C      
C=======================================================================
C=======================================================================
C
      SUBROUTINE WS_RUN(INUCLI,CHISQU_CHOICE,IABORT,I_MODE,I_FLAG,
     *                                CHISQU_PROTON,CHISQU_NEUTRS)
C
      INCLUDE   'MATDIM/NDSPEC.f'
      INCLUDE   'MATDIM/NDGAUS.f'
      INCLUDE   'MATDIM/NDIM_N.f'
      INCLUDE   'MATDIM/NDNUCL.f'
      INCLUDE   'MATDIM/NDLEXP.f'
      INCLUDE   'MATDIM/NDMAIN.f'
C      
      PARAMETER
     *         (NDIM_L=NDIM_N,NDBASE=NDIM_N+1)
C
      CHARACTER
     *          WHATEX*6,TYPCHI*6,INPSYM*6,NUCSYM*6
      CHARACTER
     *          LABEXP_PROTON*6,LABEXP_NEUTRS*6
C
      DIMENSION
     *          LABORD_AUXILP(0:NDIM_N,0:NDIM_L,0:2*NDIM_L+1),
     *          LABORD_AUXILN(0:NDIM_N,0:NDIM_L,0:2*NDIM_L+1)
      DIMENSION
     *          RFUNUP_PROTON(1:NDGAUS,0:NDIM_N,0:NDIM_L),
     *          RFUNDN_PROTON(1:NDGAUS,0:NDIM_N,0:NDIM_L),
     *          RFUNUP_NEUTRS(1:NDGAUS,0:NDIM_N,0:NDIM_L),
     *          RFUNDN_NEUTRS(1:NDGAUS,0:NDIM_N,0:NDIM_L)
      DIMENSION
     *          RAD_UP(0:NDIM_L,0:NDIM_N),
     *         	RAD_DN(0:NDIM_L,0:NDIM_N)
      DIMENSION
     *          HSORUP(1:NDBASE,1:NDBASE),
     *          HSORDN(1:NDBASE,1:NDBASE)
      DIMENSION
     *          SPENUP(1:NDBASE),
     *          SPENDN(1:NDBASE),     
     *          AUXVEC(1:NDBASE)  
C
      COMMON
     *       /LEVORD/ LABORD_PROTON(1:NDNUCL,
     *                              0:NDIM_N,0:NDIM_L,0:2*NDIM_L+1),
     *                LABORD_NEUTRS(1:NDNUCL,
     *                              0:NDIM_N,0:NDIM_L,0:2*NDIM_L+1)
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
     *       /DIMSPH/ LDSING_PROTON,LEVTHE_PROTON,
     *                              LEVEXP_PROTON(1:NDNUCL),
     *                LDSING_NEUTRS,LEVTHE_NEUTRS,
     *                              LEVEXP_NEUTRS(1:NDNUCL)
      COMMON
     *       /NODWEI/ XRNODE(1:NDGAUS,0:NDIM_L),
     *                XRWEIG(1:NDGAUS,0:NDIM_L)
      COMMON
     *       /POLPOL/ XPOLYN(1:NDGAUS,0:NDIM_N,0:NDIM_L), 
     *                XDPOLY(1:NDGAUS,0:NDIM_N,0:NDIM_L) 
      COMMON
     *       /V2COEF/ VCOEUP_PROTON(0:NDIM_L,1:NDBASE),
     *                VCOEDN_PROTON(0:NDIM_L,1:NDBASE),
     *
     *                VCOEUP_NEUTRS(0:NDIM_L,1:NDBASE),
     *                VCOEDN_NEUTRS(0:NDIM_L,1:NDBASE)
      COMMON
     *       /ENUPDN/ ENERUP_PROTON(0:NDIM_L,1:NDBASE),
     *                ENERDN_PROTON(0:NDIM_L,1:NDBASE),
     *
     *                ENERUP_NEUTRS(0:NDIM_L,1:NDBASE),
     *                ENERDN_NEUTRS(0:NDIM_L,1:NDBASE)
      COMMON
     *       /MODIFR/ RMSEXP_PROTON(1:NDNUCL),
     *                RMSEXP_NEUTRS(1:NDNUCL),
     *                RMSTHE_PROTON,RMSTHE_NEUTRS,
     *                RMODIF_PROTON,RMODIF_NEUTRS,
     *                              PUSHIN,DELTAR
      COMMON
     *       /FERGPS/ FERMEX_PROTON(1:NDNUCL),GAPEXP_PROTON(1:NDNUCL),
     *                DENSUP_PROTON(1:NDNUCL),DENSDW_PROTON(1:NDNUCL),
     *
     *                FERMEX_NEUTRS(1:NDNUCL),GAPEXP_NEUTRS(1:NDNUCL),
     *                DENSUP_NEUTRS(1:NDNUCL),DENSDW_NEUTRS(1:NDNUCL) 
      COMMON
     *       /QUANTN/ NWSSPH_PROTON(1:NDSPEC),
     *                LWSSPH_PROTON(1:NDSPEC),
     *                JWSSPH_PROTON(1:NDSPEC),
     *
     *                NWSSPH_NEUTRS(1:NDSPEC),
     *                LWSSPH_NEUTRS(1:NDSPEC),
     *                JWSSPH_NEUTRS(1:NDSPEC)
      COMMON
     *       /SPECTR/ ENETHE_PROTON(1:NDSPEC),
     *                ENETHE_NEUTRS(1:NDSPEC)
C     
      COMMON
     *       /GAUSHE/ NGAUSS,NSHELL_PROTON,NSHELL_NEUTRS
      COMMON
     *       /DATA_V/ AOSCIL_PROTON(1:NDNUCL),
     *                AOSCIL_NEUTRS(1:NDNUCL),
     *                HOMEGA(1:NDNUCL)   
      COMMON
     *       /NUCLEU/ IZ_FIX,IN_FIX,ISOSPI   
      COMMON
     *       /NUCINF/ LDNUCL,
     *                NUCSYM(1:NDNUCL),
     *                NUMB_Z(1:NDNUCL),
     *                NUMB_N(1:NDNUCL)
      COMMON
     *       /MINITE/ IDEFCN,ITECHI
      COMMON
     *       /TESFIT/ IFTEST,IFFITS,IFMESH,ISPE_P,IFMOCA
      COMMON
     *       /DENSIF/ ITRNUM
     *       /WHENST/ DIFFVA
     *       /DENSTR/ IFDENS
     *       /GNUPLO/ IGNUPL
      COMMON
     *       /ACTIVE/ IACTIV,IwMODE
     *       /FEVALS/ IFUNCT_EVALUS
     *       /MASSIV/ IMASIV_PRODUC
      COMMON
     *       /PRFILE/ ISCREN,LOGFIL,LOGCHP,LOGCHN,LOGHIP,LOGHIN,
     *                LOGPRO,LOGNEU,LOGWRI,LSCREN,IRESUL,ICONVE,
     *                NOUTPT,LOGAUX,LOGENE,LOGRAD,LOGGRD,LOGMSH,
     *                                                   LOGBIS                                                  
C  @@@ IRENE - WHAT IS THE ROLE OF THIS EPSSIL ???
      DATA 
     *        EPSSIL / 1.0e-3 /
C
C=======================================================================
C
C     This subroutine functions as the standard WS-ODD code. Its form
C     has been adapted to enable the connection with the minimisation
C     subroutines.
C
C     Some extra routines, allowing a comparison with experiment, and 
C     facilitating the definition of the ch^2 type routines are:
C
C     RDPROT_RMSRAD - reads the experimental  proton  mean square radii;
C     PROVID_RMS_NR - gives the experimental  neutron mean square radii;
C     READEX_LEVELS - reads the experimental spherical levels for doubly
C                     magic nuclei;
C
C     PREPAR_EXPLEV - transforms single particle excitation energies of
C                     one-hole and one-particle nuclei,  which are read 
C                     in SUBROUTINE READEX_LEVELS, into format suitable
C                     for comparisons with the calculations;
C
C     SERVIS - internal control subroutine
C     EXPTHE - internal control subroutine
C
C=======================================================================
C
      CALL CPUTIM('WS_RUN',1)
C
C=======================================================================
C
      IwMODE=I_MODE
C
      IF (I_MODE.EQ.0.OR.I_MODE.EQ.2) THEN
          IF (LOGWRI.GT.4) WRITE(IRESUL,'(/,''IwMODE='',I1)') IwMODE
      ELSE
          IF (LOGWRI.GT.4) WRITE(IRESUL,'(  ''IwMODE='',I1)') IwMODE
      END IF
C      
C     IFEXPE=IFFITS
C
C @@@ IRENE - what is is nonsense? NOUTPT is Standard output = 6 always
C      NOUTPT=IRESUL
      IABORT=0
C
C      IDEFCN=IDEFCN+1
C
C=======================================================================
C
      IF (IDEFCN.GT.ITECHI) THEN
C
CID          WRITE(ICONVE,'(/,80(''!''),/,''!WS_RUN'',72x,''!''/,
C     *                 ''!   The allowed maximum number of function'',
C     *                 '' evaluations ITECHI='',i4,''  exceeded   !'',/,
C     *                 ''!'',78x,''!''/,80(''!''))')
C     *                                                           ITECHI
C          
          WRITE(LOGAUX,'(/,80(''!''),/,''!WS_RUN'',72x,''!''/,
     *                 ''!   The allowed maximum number of function'',
     *                 '' evaluations ITECHI='',i4,''  exceeded   !'',/,
     *                 ''!'',78x,''!''/,80(''!''))')
     *                                                           ITECHI
          IABORT=1
C          
CID          WRITE(ICONVE,'(''RETURN from WS_RUN - before HAMMAT'')')
C          
          WRITE(LOGAUX,'(''RETURN from WS_RUN - before HAMMAT'')')
C
          WRITE(000000,'(''RETURN from WS_RUN - before HAMMAT'')')
C
          CALL CPUTIM('WS_RUN',0)
C          
          RETURN
C
      END IF
C
C=======================================================================
C     Solving the Schroedinger Equation (params in POTPOT)
C=======================================================================
C
      DIFFVA=999.
      ITRNUM=0
C
      IZ_FIX=NUMB_Z(INUCLI)
      IN_FIX=NUMB_N(INUCLI)
C
      IF (IFDENS.EQ.0) THEN
C
          IF (LOGWRI.GT.4) THEN
              WRITE(LOGFIL,'(15X,''Entering HAMMAT from WS_RUN [1], '',
     *                           ''IFDENS=0 IZ_FIX='',I2,
     *                                   '' IN_FIX='',I3)') 
     *                                      IZ_FIX,IN_FIX
          END IF 
C
          CALL HAMMAT(INUCLI,XPOLYN,XDPOLY,XRNODE,XRWEIG,NGAUSS,
     *                       NSHELL_PROTON,NDGAUS,NSHELL_NEUTRS,
     *                HSORUP,HSORDN,SPENUP,SPENDN,AUXVEC,IABORT,
     *                              CMATUP_PROTON,CMATDN_PROTON,
     *                              CMATUP_NEUTRS,CMATDN_NEUTRS,
     *                              ENEKIN_PROTON,ENEKIN_NEUTRS,
     *                              LABORD_PROTON,LABORD_NEUTRS)
C
      END IF
C
C=======================================================================
C      
      IF (IFDENS.EQ.1) THEN
C      
   1      CONTINUE 
C
          IF (SQRT(DIFFVA).GE.EPSSIL) THEN 
C
              IF (LOGWRI.GT.4) THEN
              WRITE(LOGFIL,'(15X,''Entering HAMMAT from WS_RUN [2], '',
     *                           ''DIFFVA= '',F22.13,/,45X,
     *                           ''EPSSIL= '',F22.13)') DIFFVA,EPSSIL
              END IF 
C
              CALL HAMMAT(INUCLI,XPOLYN,XDPOLY,XRNODE,XRWEIG,NGAUSS,
     *                           NSHELL_PROTON,NDGAUS,NSHELL_NEUTRS,
     *                    HSORUP,HSORDN,SPENUP,SPENDN,AUXVEC,IABORT,
     *                                  CMATUP_PROTON,CMATDN_PROTON,
     *                                  CMATUP_NEUTRS,CMATDN_NEUTRS,
     *                                  ENEKIN_PROTON,ENEKIN_NEUTRS,
     *                                  LABORD_PROTON,LABORD_NEUTRS)
              IF (IFTEST.EQ.1) THEN
C              
                  WRITE(21,'()')
                  WRITE(21,'(''ITRUM='',I3.2,''  DIFFVA='',ES15.5,
     *                       ''  SQRT(DIFFVA)='',ES15.5)')
     *                         ITRNUM,DIFFVA,SQRT(DIFFVA)
              END IF
               
              ITRNUM=ITRNUM+1
C
	      IF (IABORT.EQ.1) THEN
C
                  WRITE(LOGAUX,'(''RETURN in WS_RUN after HAMMAT: '',
     *                           ''IABORT= '',I1,'' and ITRNUM= '',I5)')
     *                             IABORT,              ITRNUM
                  WRITE(LOGFIL,'(''RETURN in WS_RUN after HAMMAT: '',
     *                           ''IABORT= '',I1,'' and ITRNUM= '',I5)')
     *                             IABORT,              ITRNUM
                  WRITE(000000,'(''RETURN in WS_RUN after HAMMAT: '',
     *                           ''IABORT= '',I1,'' and ITRNUM= '',I5)')
     *                             IABORT,              ITRNUM 
C
                  CALL CPUTIM('WS_RUN',0)
C
                  RETURN
C
              END IF
          ELSE
              GO TO 2
          END IF
C
          GO TO 1
C
   2      CONTINUE 
C
          IF (LOGWRI.GT.0) THEN
              WRITE(LOGFIL,'(15X,''WS_RUN: ITRNUM='',I3)') ITRNUM
          END IF
C
      END IF
C
C=======================================================================
C      
      DO N=0,NDIM_N
         DO L=0,NDIM_L
            DO J=0,2*NDIM_L+1
               LABORD_AUXILP(N,L,J)=LABORD_PROTON(INUCLI,N,L,J)
               LABORD_AUXILN(N,L,J)=LABORD_NEUTRS(INUCLI,N,L,J)
            END DO
         END DO
      END DO
C
      RMSTHE_PROTON=0.
      RMSTHE_NEUTRS=0.
C
      ISOSPI=1      
C
      IF (LOGWRI.GT.5) THEN
          WRITE(LOGFIL,'(/,15X,''Entering NUCRAD from WS_RUN'',1X,
     *                                             ''ISOSPI='',I1)')
     *                                               ISOSPI
      END IF  
C
C=======================================================================
C
      IF (IZ_FIX.GT.82) GO TO 100
C
C=======================================================================
C
                         LDSPEC=LDSING_PROTON
      CALL NUCRAD(NDSPEC,LDSPEC,NDGAUS,NGAUSS,NDBASE,XPOLYN,
     *            NDIM_L,NDIM_N,RAD_UP,RAD_DN,XRNODE,XRWEIG,
     *            AOSCIL_PROTON,RFUNUP_PROTON,RFUNDN_PROTON,
     *            NSHELL_PROTON,ENERUP_PROTON,ENERDN_PROTON,                                           
     *            NWSSPH_PROTON,LWSSPH_PROTON,JWSSPH_PROTON,
     *            RMSTHE_PROTON,RMSPAI_PROTON,RMSNOP_PROTON,
     *            LABORD_AUXILP,CMATUP_PROTON,CMATDN_PROTON,  
     *            ENETHE_PROTON,VCOEUP_PROTON,VCOEDN_PROTON,
     *                                 ISOSPI,IZ_FIX,INUCLI)
C
      ISOSPI=0      
C
      IF (LOGWRI.GT.5) THEN
          WRITE(LOGFIL,'(/,15X,''Entering NUCRAD from WS_RUN'',1X,
     *                                             ''ISOSPI='',I1)')
     *                                               ISOSPI
      END IF  
                         LDSPEC=LDSING_NEUTRS
      CALL NUCRAD(NDSPEC,LDSPEC,NDGAUS,NGAUSS,NDBASE,XPOLYN,
     *            NDIM_L,NDIM_N,RAD_UP,RAD_DN,XRNODE,XRWEIG,
     *            AOSCIL_NEUTRS,RFUNUP_NEUTRS,RFUNDN_NEUTRS,
     *            NSHELL_NEUTRS,ENERUP_NEUTRS,ENERDN_NEUTRS,                                           
     *            NWSSPH_NEUTRS,LWSSPH_NEUTRS,JWSSPH_NEUTRS,
     *            RMSTHE_NEUTRS,RMSPAI_NEUTRS,RMSNOP_NEUTRS,
     *            LABORD_AUXILN,CMATUP_NEUTRS,CMATDN_NEUTRS,  
     *            ENETHE_NEUTRS,VCOEUP_NEUTRS,VCOEDN_NEUTRS,
     *                                 ISOSPI,IN_FIX,INUCLI)
C
C=======================================================================
C
      IF (IZ_FIX.GT.82) GO TO 100
C
C=======================================================================
C     Testing the quality of the experimental spectra for
C     doubly-magic, spherical nuclei           (IFEXPE=1)
C=======================================================================
C
      ISOSPI=1
C
      IF (LOGWRI.GT.5) THEN
          WRITE(LOGFIL,'(/,15X,''Entering EXPTHE from WS_RUN'',1X,
     *                                             ''ISOSPI='',I1)')
     *                                               ISOSPI
      END IF  
C      
      CALL EXPTHE(IZ_FIX,IN_FIX,DENSUP_PROTON,DENSDW_PROTON,
     *            FERMEX_PROTON,GAPEXP_PROTON,IDEFCN,ISOSPI,
     *                   IPARAM,CHISQU_CHOICE,I_FLAG,INUCLI)
C
      CHISQU_PROTON=CHISQU_CHOICE
C
      ISOSPI=0
C
      IF (LOGWRI.GT.5) THEN
          WRITE(LOGFIL,'(/,15X,''Entering EXPTHE from WS_RUN'',1X,
     *                                             ''ISOSPI='',I1)')
     *                                               ISOSPI
      END IF  
C
      CALL EXPTHE(IZ_FIX,IN_FIX,DENSUP_NEUTRS,DENSDW_NEUTRS,
     *            FERMEX_NEUTRS,GAPEXP_NEUTRS,IDEFCN,ISOSPI,
     *                   IPARAM,CHISQU_CHOICE,I_FLAG,INUCLI)
C
      CHISQU_NEUTRS=CHISQU_CHOICE
C 
      IFUNCT_EVALUS=IFUNCT_EVALUS+1
C
C=======================================================================
C
 100  CONTINUE
C
C=======================================================================
C     Printing actual parameters and corresponding Chi^2 (AS ORIG.)
C=======================================================================
C
c ###      IF (I_FLAG.EQ.1) THEN 
c ###          CALL INPRIN(IDEFCN,IACTIV,I_MODE,CHISQU_PROTON,
c ###     *                                     CHISQU_NEUTRS)
c ###      END IF
C ###
C ###    !!!!  NOW THE 'CALL' IS IN SUBROUTINE FUNMIN   !!!!
C=======================================================================
C
      IF (LOGWRI.GT.4) THEN
          WRITE(LOGFIL,'(12X,''Exiting  WS_RUN'')')
      END IF
C      
C=======================================================================
C
      CALL CPUTIM('WS_RUN',0)   
C      
C=======================================================================
C      
      RETURN
      END
C      
C=======================================================================
C=======================================================================
C
      SUBROUTINE KINMAT(POLYNL,DPOLYN,POLYNX,RNODES,RWEIGH,ENEKIN,
     *                  HOMEGA,NDBASE,NDIM_N,NDIM_L,NSHELL,NDGAUS,
     *                                       NGAUSS,NDNUCL,INUCLI)
C
      DIMENSION
     *          HOMEGA(1:NDNUCL)
      DIMENSION
     *          ENEKIN(1:NDNUCL,1:NDBASE,1:NDBASE,0:NDIM_L)
      DIMENSION
     *          POLYNL(1:NDGAUS,0:NDIM_N,0:NDIM_L), 
     *          POLYNX(1:NDGAUS,0:NDIM_N,0:NDIM_L), 
     *          DPOLYN(1:NDGAUS,0:NDIM_N,0:NDIM_L) 
      DIMENSION
     *          RNODES(1:NDGAUS,0:NDIM_L),
     *          RWEIGH(1:NDGAUS,0:NDIM_L)
C
      COMMON
     *       /PRFILE/ ISCREN,LOGFIL,LOGCHP,LOGCHN,LOGHIP,LOGHIN,
     *                LOGPRO,LOGNEU,LOGWRI,LSCREN,IRESUL,ICONVE,
     *                NOUTPT,LOGAUX,LOGENE,LOGRAD,LOGGRD,LOGMSH,
     *                                                   LOGBIS  
C
C @@@ How does this routine distinguish
C     between the protons and the neutrons? MASS DIFFEENCES?
C=======================================================================     
C     This subroutine calculates the matrix elements of the kinetic
C     energy term [the H.O. scaling factors may be dependent on the 
C     nucleonic mass]. Formally the difference between the neutrons
C     and protons comes only through the numbers of shells, because
C     the matrix elements take the same form.
C=======================================================================     
C
      DO LQNUMB=NSHELL,0,-1
         DO NQNUMB=0,(NSHELL-LQNUMB)/2
            DO I_PNTS=1,NGAUSS
C
            TERM_1=0.50*(FLOAT(LQNUMB+1)-RNODES(I_PNTS,LQNUMB))
C
            POLYNX(I_PNTS,NQNUMB,LQNUMB)=POLYNL(I_PNTS,NQNUMB,LQNUMB)
     *                                  *TERM_1
     *                                  +DPOLYN(I_PNTS,NQNUMB,LQNUMB)
     *                                  *RNODES(I_PNTS,LQNUMB)
            END DO
         END DO
      END DO
C_______________________________________________________________________
C
      DO LQNUMB=NSHELL,0,-1
C
C        Here we are going to construct the T-sub-block at given L
C
         DO N1NUMB=0,(NSHELL-LQNUMB)/2
            DO N2NUMB=0,(NSHELL-LQNUMB)/2       
C       
               T_ELEM=0.00
C       
               DO I_PNTS=1,NGAUSS
C
                  TERM_1=POLYNX(I_PNTS,N1NUMB,LQNUMB)
     *                  *POLYNX(I_PNTS,N2NUMB,LQNUMB)
C
                  TERM_2=POLYNL(I_PNTS,N1NUMB,LQNUMB)*0.250
     *                  *POLYNL(I_PNTS,N2NUMB,LQNUMB)
     *                  *FLOAT(LQNUMB*(LQNUMB+1))
C       
                  T_ELEM=T_ELEM+(TERM_1+TERM_2)*RWEIGH(I_PNTS,LQNUMB)    
C                  
               END DO
C
C              The factor of 2.0 in the following line takes into
C              account the fact that  XNORNL  are defined without
C              the factor sqrt(2). XNORNL are defined without the 
C              stretching factor "a" either, but this one cancels 
C              out in the calculations.
C
               ENEKIN(INUCLI,N1NUMB+1,N2NUMB+1,LQNUMB)=2.00
     *	                                               *HOMEGA(INUCLI)
     *                                                *T_ELEM
C
C=======================================================================
C      
            END DO
         END DO
C      
      END DO
C
C=======================================================================
C
      IF (LOGWRI.GT.5) THEN
          WRITE(LOGFIL,'(''Exiting  KINMAT'')')
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
      SUBROUTINE HAMMAT(INUCLI,POLYNL,DPOLYN,RNODES,RWEIGH,NGAUSS,
     *                         NSHELL_PROTON,NDGAUS,NSHELL_NEUTRS,
     *                  HSORUP,HSORDN,SPENUP,SPENDN,AUXVEC,IABORT,
     *                                CMATUP_PROTON,CMATDN_PROTON,
     *                                CMATUP_NEUTRS,CMATDN_NEUTRS,
     *                                ENEKIN_PROTON,ENEKIN_NEUTRS,
     *                                LABORD_PROTON,LABORD_NEUTRS)
C     
      INCLUDE  'MATDIM/NDSPEC.f'
      INCLUDE  'MATDIM/NITERA.f'
      INCLUDE  'MATDIM/NDNUCL.f'
      INCLUDE  'MATDIM/NDPARS.f'
      INCLUDE  'MATDIM/NDIMaN.f'
      INCLUDE  'MATDIM/LDGAUS.f'
      INCLUDE  'MATDIM/NSIMPS.f'
      INCLUDE  'MATDIM/NDIM_N.f'
C      
      PARAMETER
     *         (NDIMaL=NDIMaN,NDIM_L=NDIM_N,NDBASE=NDIM_N+1)
C
      CHARACTER
     *          STRNG1*1,STRNG2*23
      CHARACTER
     *          LABTHE_PROTON*6,LABTHE_NEUTRS*6
      DIMENSION
     *          ENEKIN_PROTON(1:NDNUCL,1:NDBASE,1:NDBASE,0:NDIM_L),
     *          ENEKIN_NEUTRS(1:NDNUCL,1:NDBASE,1:NDBASE,0:NDIM_L)
      DIMENSION
     *          CAUXUP_PROTON(0:NDIM_N,1:NDBASE),
     *          CAUXDN_PROTON(0:NDIM_N,1:NDBASE),
     *
     *          CAUXUP_NEUTRS(0:NDIM_N,1:NDBASE),
     *          CAUXDN_NEUTRS(0:NDIM_N,1:NDBASE)
      DIMENSION
     *          CMATUP_PROTON(0:NDIM_L,1:NDBASE,1:NDBASE),
     *          CMATDN_PROTON(0:NDIM_L,1:NDBASE,1:NDBASE),
     *
     *          CMATUP_NEUTRS(0:NDIM_L,1:NDBASE,1:NDBASE),
     *          CMATDN_NEUTRS(0:NDIM_L,1:NDBASE,1:NDBASE)
      DIMENSION
     *          POLYNL(1:NDGAUS,0:NDIM_N,0:NDIM_L),
     *          DPOLYN(1:NDGAUS,0:NDIM_N,0:NDIM_L)
      DIMENSION
     *          LAUXIL(1:NDSPEC),JAUXIL(1:NDSPEC),
     *                           NAUXIL(1:NDSPEC)
      DIMENSION
     *          HSORUP(1:NDBASE,1:NDBASE),
     *          HSORDN(1:NDBASE,1:NDBASE)
      DIMENSION
     *          VPOTUP(1:NDBASE,1:NDBASE),
     *          VPOTDN(1:NDBASE,1:NDBASE)
      DIMENSION
     *          VCOEUP(0:NDIM_L,1:NDBASE),
     *          VCOEDN(0:NDIM_L,1:NDBASE)
      DIMENSION
     *          RNODES(1:NDGAUS,0:NDIM_L),
     *          RWEIGH(1:NDGAUS,0:NDIM_L)
      DIMENSION
     *          INDEXS(1:NDSPEC)
      DIMENSION
     *          SPENUP(1:NDBASE),
     *          SPENDN(1:NDBASE),     
     *          AUXVEC(1:NDBASE)     
      DIMENSION
     *          ENETHE(1:NDSPEC)
      DIMENSION
     *          LABORD_PROTON(1:NDNUCL,
     *                        0:NDIM_N,0:NDIM_L,0:2*NDIM_L+1),
     *          LABORD_NEUTRS(1:NDNUCL,
     *                        0:NDIM_N,0:NDIM_L,0:2*NDIM_L+1)
C
      COMMON
     *       /DATA_V/ AOSCIL_PROTON(1:NDNUCL),
     *                AOSCIL_NEUTRS(1:NDNUCL),
     *                HOMEGA(1:NDNUCL)
      COMMON
     *       /SEARCH/ IFTAKE(1:NDPARS)
      COMMON
     *       /DIMSPH/ LDSING_PROTON,LEVTHE_PROTON,
     *                              LEVEXP_PROTON(1:NDNUCL),
     *                LDSING_NEUTRS,LEVTHE_NEUTRS,
     *                              LEVEXP_NEUTRS(1:NDNUCL)
      COMMON
     *       /ENUPDN/ ENERUP_PROTON(0:NDIM_L,1:NDBASE),
     *                ENERDN_PROTON(0:NDIM_L,1:NDBASE),
     *
     *                ENERUP_NEUTRS(0:NDIM_L,1:NDBASE),
     *                ENERDN_NEUTRS(0:NDIM_L,1:NDBASE)
      COMMON
     *       /EUDAUX/ ENERUP_AUXPRO(0:NDIM_L,1:NDBASE),
     *                ENERDN_AUXPRO(0:NDIM_L,1:NDBASE),
     *
     *                ENERUP_AUXNEU(0:NDIM_L,1:NDBASE),
     *                ENERDN_AUXNEU(0:NDIM_L,1:NDBASE)
      COMMON
     *       /CUDAUX/ CMATUP_AUXPRO(0:NDIM_L,1:NDBASE,1:NDBASE),
     *                CMATDN_AUXPRO(0:NDIM_L,1:NDBASE,1:NDBASE),
     *
     *                CMATUP_AUXNEU(0:NDIM_L,1:NDBASE,1:NDBASE),
     *                CMATDN_AUXNEU(0:NDIM_L,1:NDBASE,1:NDBASE)
      COMMON
     *       /ENRGYY/ ENERGY_PROTON(1:NDSPEC,0:NITERA), 
     *                ENERGY_NEUTRS(1:NDSPEC,0:NITERA)
      COMMON
     *       /SPECTR/ ENETHE_PROTON(1:NDSPEC),
     *                ENETHE_NEUTRS(1:NDSPEC)
      COMMON
     *       /THELAB/ LABTHE_PROTON(1:NDSPEC),
     *                LABTHE_NEUTRS(1:NDSPEC)  
      COMMON
     *       /QUANTN/ NWSSPH_PROTON(1:NDSPEC),
     *                LWSSPH_PROTON(1:NDSPEC),
     *                JWSSPH_PROTON(1:NDSPEC),
     *
     *                NWSSPH_NEUTRS(1:NDSPEC),
     *                LWSSPH_NEUTRS(1:NDSPEC),
     *                JWSSPH_NEUTRS(1:NDSPEC)
      COMMON
     *       /V2COEF/ VCOEUP_PROTON(0:NDIM_L,1:NDBASE),
     *                VCOEDN_PROTON(0:NDIM_L,1:NDBASE),
     *
     *                VCOEUP_NEUTRS(0:NDIM_L,1:NDBASE),
     *                VCOEDN_NEUTRS(0:NDIM_L,1:NDBASE)
      COMMON
     *       /V2_COE/ V2_PUP_VECTOR(0:NDIM_L,1:NDBASE),
     *                V2_PDW_VECTOR(0:NDIM_L,1:NDBASE),
     *                DUP_LN(0:NDIM_L,1:NDBASE),
     *                DDW_LN(0:NDIM_L,1:NDBASE),
     *                ENERGY,XLAMBD,DELTA2
      COMMON
     *       /POTPOT/ PARPOT(1:NDPARS)
      COMMON
     *       /BCSENE/ ENEMAX_BCSPAI(1:NDNUCL,1:2)
      COMMON
     *       /NUCLEU/ IZ_FIX,IN_FIX,ISOSPI    
      COMMON
     *       /DENSIF/ ITRNUM
      COMMON
     *       /DENSTR/ IFDENS
      COMMON
     *       /IFPAIR/ IF_PAI,IFGPAI,IFDPAI
      COMMON
     *       /TESFIT/ IFTEST,IFFITS,IFMESH,ISPE_P,IFMOCA
      COMMON
     *       /MASSIV/ IMASIV_PRODUC
      COMMON
     *       /CNTCNT/ ICOUNT_DENGRD,ICOUNT_FUNMIN,ICOUNT_HAMMAT,
     *                ICOUNT_EXPTHE,ICOUNT_CHOICE,ICOUNT_DIAMAT
      COMMON
     *       /PRFILE/ ISCREN,LOGFIL,LOGCHP,LOGCHN,LOGHIP,LOGHIN,
     *                LOGPRO,LOGNEU,LOGWRI,LSCREN,IRESUL,ICONVE,
     *                NOUTPT,LOGAUX,LOGENE,LOGRAD,LOGGRD,LOGMSH,
     *                                                   LOGBIS  
C
      DATA
     *        EPSTSH        / 1.0E-4 / ! Hermiticity test
      DATA
     *        SPENER_MXPRIN / 10.000 / ! No print above that energy
C
C=======================================================================     
C     This subroutine calculates the matrix elements of the Hamiltonian
C     and diagonalises (through DIAMAT) the corresponding blocks of the
C     Hamiltonian with parallel, and anti-parallel "l-s" configurations
C=======================================================================     
C
C     We refer to the programmed mathematical expressions using Chapter
C     and equation numbers in the following documents:
C
C     Document 1: "Interactions and mean field Hamiltonians" =>> D1
C
C     Document 2: "Nuclear Pairing, Nuclear Superfludity..." =>> D2
C
      CALL CPUTIM('HAMMAT',1)
C
C=======================================================================     
C      
      ICOUNT_HAMMAT=ICOUNT_HAMMAT+1
C
C=======================================================================     
C
      IF (IMASIV_PRODUC.EQ.1 .AND. LOGWRI.GT.4) THEN
C
          WRITE(LOGFIL,'(/,15X,''Entering HAMMAT '',
     *                         ''for massive production'')')
C
          WRITE(LOGFIL,'(15X,''IZ_FIX= '',I3,'' IN_FIX= '',I3,1X,
     *                       ''AOSCIL_PROTON= '',F8.4,1X,
     *                       ''AOSCIL_NEUTRS= '',F8.4,/)')
     *                         IZ_FIX,IN_FIX,AOSCIL_PROTON(INUCLI),
     *                                       AOSCIL_NEUTRS(INUCLI)
C
          WRITE(LOGFIL,'(12X,''Parameters: '',/)')
C
          DO IPARAM=1,NDPARS
              IF (IFTAKE(IPARAM).EQ.1) THEN
                  WRITE(LOGFIL,'(12X,I3,f9.4)') IPARAM,PARPOT(IPARAM)
              END IF
          END DO
C         
      END IF
C
C=======================================================================     
C
      IABORT=0
C
C=======================================================================     
C
C     Begin by introducing all the constant constants ...
C
      CALL INTROD(INUCLI)
C
C=======================================================================
C=======================================================================
C     Beginning with the algorithm for the  p r o t o n s
C=======================================================================
C=======================================================================
C
      IF (LOGWRI.GT.0) WRITE(NOUTPT,'(/)')
C
      NCOUNT=0
      ISOSPI=1
C
      IF (IFFITS.EQ.1 .OR. IFMESH.EQ.1) THEN
C
          IF (LOGWRI.GT.5) THEN
              WRITE(LOGFIL,'(15X,''Entering SELECT_VSODEN from '',
     *                                                 ''HAMMAT'')')
          END IF
C
          CALL SELECT_VSODEN
C
      END IF
C
C=======================================================================
C
C     Product L * S :
C
C     FACTUP - spin parallel configuration s=1/2 
C     
C     FACTDN - spin antiparallel anti-configuration s=-1/2 
C
C=======================================================================
C  
      DO LQNUMB=NSHELL_PROTON,0,-1
C          
         FACTUP=+0.5*LQNUMB
         FACTDN=-0.5*(LQNUMB+1)     
C
C        Here we are going to construct the H-sub-block at given L
C 
         DO N1NUMB=0,(NSHELL_PROTON-LQNUMB)/2
            SPENUP(N1NUMB+1)=0.0
            SPENDN(N1NUMB+1)=0.0
            DO N2NUMB=0,(NSHELL_PROTON-LQNUMB)/2       
               HSORUP(N1NUMB+1,N2NUMB+1)=0.0
               HSORDN(N1NUMB+1,N2NUMB+1)=0.0
               VPOTUP(N1NUMB+1,N2NUMB+1)=0.0
               VPOTDN(N1NUMB+1,N2NUMB+1)=0.0
            END DO
         END DO
C
C        We are within the l-do-loop
C_______________________________________________________________________
C       
         VCELEM=0.0
         VSELEM=0.0
C_______________________________________________________________________
C                                                            Integration
         DO I_PNTS=1,NGAUSS
C
C           We are using the weight factor e^{-z} z^{l-1/2} and
C           consequently the potential is multiplied by "z" for
C           details, cf. D1, eqs.(4.8.8) - (4.8.10)
C
C           We selected to program Eq.(4.8.10) because this one
C           corresponds directly to the form of the kinetic energy
C
            ZNODES=RNODES(I_PNTS,LQNUMB) 
            ZWEIGH=RWEIGH(I_PNTS,LQNUMB)
         
            VCELEM=VCENTR(INUCLI,ZNODES,I_PNTS,LQNUMB)
         
            IF (ITRNUM.EQ.0) THEN ! Using the analytical SO-Potential
C                
                VSELEM=V_SORB(INUCLI,ZNODES)
C
            ELSE   ! Using the density-dependent SO-Potential
C             
                VSELEM=V_SORB_DENSIT(INUCLI,ZNODES,I_PNTS,LQNUMB)
C
            END IF
C_______________________________________________________________________
C         
            DO N1NUMB=0,(NSHELL_PROTON-LQNUMB)/2
               DO N2NUMB=0,(NSHELL_PROTON-LQNUMB)/2       
C
                  FACTOR=ZNODES*ZWEIGH
     *                  *POLYNL(I_PNTS,N1NUMB,LQNUMB)
     *                  *POLYNL(I_PNTS,N2NUMB,LQNUMB)
C
                  VPOTUP(N1NUMB+1,N2NUMB+1)=VPOTUP(N1NUMB+1,N2NUMB+1)
     *                                     +VCELEM       *FACTOR
     *                                     +VSELEM*FACTUP*FACTOR
C     
                  VPOTDN(N1NUMB+1,N2NUMB+1)=VPOTDN(N1NUMB+1,N2NUMB+1)
     *                                     +VCELEM       *FACTOR
     *                                     +VSELEM*FACTDN*FACTOR
C
               END DO ! End of the n2 do-loop 
C            
            END DO ! End of the n1 do-loop  
C_______________________________________________________________________
C             
         END DO ! Ending Gauss-Laguerre integration for the potentials
C_______________________________________________________________________
C   
C        Adding the kinetic energy matrix elements already calculated 
C        in KINMAT to the matrix elements of the potential       
C         
         DO N1NUMB=0,(NSHELL_PROTON-LQNUMB)/2
            DO N2NUMB=0,(NSHELL_PROTON-LQNUMB)/2  
C         
               EKINET=ENEKIN_PROTON(INUCLI,N1NUMB+1,N2NUMB+1,LQNUMB)
C
               HSORUP(N1NUMB+1,N2NUMB+1)
     *        =
     *         EKINET+VPOTUP(N1NUMB+1,N2NUMB+1)
C               
               HSORDN(N1NUMB+1,N2NUMB+1)
     *        =
     *         EKINET+VPOTDN(N1NUMB+1,N2NUMB+1)
C                      
             END DO ! End of the n2 do-loop 
         END DO ! End of the n1 do-loop  
C
C=======================================================================
C
         LDBASE=(NSHELL_PROTON-LQNUMB)/2+1       
C
C=======================================================================
C        Symmetry test for the Hamiltonian matrix under the l do-loop
C=======================================================================
C
         DO N1NUMB=0,(NSHELL_PROTON-LQNUMB)/2
            DO N2NUMB=N1NUMB,(NSHELL_PROTON-LQNUMB)/2       
C
               IF (ABS(HSORUP(N1NUMB+1,N2NUMB+1)
     *                -HSORUP(N2NUMB+1,N1NUMB+1)).GT.EPSTSH) THEN
C
                   WRITE(NOUTPT,'(''ITRNUM='',I3,3X,''N1='',I3,3X,
     *                                              ''N2='',I3,3X,
     *                            ''HSORUP(N1,N2)='',F30.25,3X,
     *                            ''HSORUP(N2,N1)='',F30.25)')
     *
     *                              ITRNUM,N1NUMB,N2NUMB,
     *                              HSORUP(N1NUMB+1,N2NUMB+1),
     *                              HSORUP(N2NUMB+1,N1NUMB+1)
C
                   STOP 'Non-hermiticity in matrix HSORUP -PROTONS- !'
               END IF
C
               IF (ABS(HSORDN(N1NUMB+1,N2NUMB+1)
     *                -HSORDN(N2NUMB+1,N1NUMB+1)).GT.EPSTSH) THEN
C
                   WRITE(NOUTPT,'(''ITRNUM='',I3,3X,''N1='',I3,3X,
     *                                              ''N2='',I3,3X,
     *                            ''HSORDN(N1,N2)='',F30.25,3X,
     *                            ''HSORDN(N2,N1)='',F30.25)')
     *
     *                              ITRNUM,N1NUMB,N2NUMB,
     *                              HSORDN(N1NUMB+1,N2NUMB+1),
     *                              HSORDN(N2NUMB+1,N1NUMB+1)
C
                   STOP 'Non hermiticity in matrix HSORDN -PROTONS- !'
               END IF
C
            END DO
         END DO
C
C=======================================================================
C
C        Diagonalising  the real-symmetric Hamiltonian-matrix
C        including the calculation of eigenvectors by setting:
C                                                EIGVEC=.TRUE. 
C                                                                
         CALL DIAMAT(HSORUP,SPENUP,AUXVEC,NDBASE,LDBASE,.TRUE.)
C
C        Normalised eigenvectors are stored in columns HSORUP
C
C=======================================================================
C
C        For every l-value we have there is a spin-parallel-L
C        configuration, say spin-up,  and another one that is
C        anti-parallel, say spin-down, the latter exclusively
C        for L>0.  We proceed by storing the spin-up results:
C
         DO NWSNUM=1,(NSHELL_PROTON-LQNUMB)/2+1
            DO N1NUMB=0,(NSHELL_PROTON-LQNUMB)/2
C
               CAUXUP_PROTON(N1NUMB,NWSNUM)=HSORUP(N1NUMB+1,NWSNUM)
C               
               CMATUP_PROTON(LQNUMB,N1NUMB+1,NWSNUM)
     *        =              HSORUP(N1NUMB+1,NWSNUM)
C 
               ENERUP_PROTON(LQNUMB,NWSNUM)=SPENUP(NWSNUM)
C
            END DO
            
         END DO  
C
C=======================================================================
C
         DO I=1,LDBASE
                   NCOUNT =NCOUNT+1
            ENETHE(NCOUNT)=SPENUP(I)
            LAUXIL(NCOUNT)=LQNUMB
            JAUXIL(NCOUNT)=LQNUMB+LQNUMB+1
            NAUXIL(NCOUNT)=I-1
            INDEXS(NCOUNT)=NCOUNT
         END DO    
C
C=======================================================================
C
C        If L>0, we proceed to treat the spin-antiparallel-L case
C
         IF (LQNUMB.NE.0) THEN  ! Diagonalisation for L > 0, DOWN
C                                                    EIGVEC=.TRUE.
             CALL DIAMAT(HSORDN,SPENDN,AUXVEC,NDBASE,LDBASE,.TRUE.)
C
             DO NWSNUM=1,(NSHELL_PROTON-LQNUMB)/2+1
                DO N1NUMB=0,(NSHELL_PROTON-LQNUMB)/2 
C
                   CAUXDN_PROTON(N1NUMB,NWSNUM)=HSORDN(N1NUMB+1,NWSNUM)
                   
                   CMATDN_PROTON(LQNUMB,N1NUMB+1,NWSNUM)
     *        =                  HSORDN(N1NUMB+1,NWSNUM)
C
                   ENERDN_PROTON(LQNUMB,NWSNUM)=SPENDN(NWSNUM)
C
                END DO
             END DO
C
             DO I=1,LDBASE
                       NCOUNT =NCOUNT+1
                ENETHE(NCOUNT)=SPENDN(I)
                LAUXIL(NCOUNT)=LQNUMB
                JAUXIL(NCOUNT)=LQNUMB+LQNUMB-1
                NAUXIL(NCOUNT)=I-1
                INDEXS(NCOUNT)=NCOUNT
             END DO 
C
         END IF 
C______________________________________________________________________
C      
      END DO   !LQNUMB   
C
C=======================================================================
C
      DO I=1,NCOUNT
          INDEXS(I)=I
      END DO
      
      CALL ORDHEA(ENETHE,INDEXS,NCOUNT,NDSPEC)
C           
C=======================================================================
C
      LDSING_PROTON=NCOUNT
C
      NOCCUP_PROTON=0
C
      DO I=1,NCOUNT
C
         INDOLD=INDEXS(I)
C
         NWSSPH_PROTON(I)=NAUXIL(INDOLD)
         LWSSPH_PROTON(I)=LAUXIL(INDOLD)
         JWSSPH_PROTON(I)=JAUXIL(INDOLD)
C
         ENERGY_PROTON(I,ITRNUM)=ENETHE(I)
         ENETHE_PROTON(I       )=ENETHE(I)
C
         IF (IZ_FIX.GT.82 .AND. NOCCUP_PROTON.LT.IZ_FIX) THEN
C
            NOCCUP_PROTON=NOCCUP_PROTON+JWSSPH_PROTON(I)+1
            LABORD_PROTON(INUCLI,NWSSPH_PROTON(I),LWSSPH_PROTON(I),
     *                                            JWSSPH_PROTON(I))=1
C
         END IF
C      
      END DO
C
C=======================================================================
C=======================================================================
C     Repeating the same algorithm  for   n e u t r o n s
C=======================================================================
C=======================================================================
C
      NCOUNT=0
      ISOSPI=0
C      
      XNUMBN=0.0
C
C=======================================================================
C
C     Product L * S :
C
C     FACTUP - spin parallel configuration s=1/2 
C     
C     FACTDN - spin antiparallel anti-configuration s=-1/2 
C
C=======================================================================
C  
      DO LQNUMB=NSHELL_NEUTRS,0,-1
C          
         FACTUP=+0.5*LQNUMB
         FACTDN=-0.5*(LQNUMB+1)         
C
C       Here we are going to construct the H-sub-block at given L
C
         DO N1NUMB=0,(NSHELL_NEUTRS-LQNUMB)/2
C
            SPENUP(N1NUMB+1)=0.0
            SPENDN(N1NUMB+1)=0.0
C
            DO N2NUMB=0,(NSHELL_NEUTRS-LQNUMB)/2       
               HSORUP(N1NUMB+1,N2NUMB+1)=0.0
               HSORDN(N1NUMB+1,N2NUMB+1)=0.0
               VPOTUP(N1NUMB+1,N2NUMB+1)=0.0
               VPOTDN(N1NUMB+1,N2NUMB+1)=0.0
            END DO
C
         END DO
C
C        We are within the l-do-loop
C_______________________________________________________________________
C       
         VCELEM=0.0
         VSELEM=0.0
C_______________________________________________________________________
C                                                         Integration
         DO I_PNTS=1,NGAUSS
C
C        We are using the weight factor e^{-z} z^{l-1/2}
C        and consequently the potential is multiplied by "z"
C
            ZNODES=RNODES(I_PNTS,LQNUMB) 
            ZWEIGH=RWEIGH(I_PNTS,LQNUMB)
         
            VCELEM=VCENTR(INUCLI,ZNODES,I_PNTS,LQNUMB)
         
            IF (ITRNUM.EQ.0) THEN ! Using the analytical SO-Potential
C                
                VSELEM=V_SORB(INUCLI,ZNODES)
C
            ELSE   ! Using the density-dependent SO-Potential
C                
                VSELEM=V_SORB_DENSIT(INUCLI,ZNODES,I_PNTS,LQNUMB)
C
            END IF
C_______________________________________________________________________
C         
            DO N1NUMB=0,(NSHELL_NEUTRS-LQNUMB)/2
               DO N2NUMB=0,(NSHELL_NEUTRS-LQNUMB)/2       
C
                  FACTOR=ZNODES*ZWEIGH
     *                  *POLYNL(I_PNTS,N1NUMB,LQNUMB)
     *                  *POLYNL(I_PNTS,N2NUMB,LQNUMB)
C
                  VPOTUP(N1NUMB+1,N2NUMB+1)=VPOTUP(N1NUMB+1,N2NUMB+1)
     *                                     +VCELEM*FACTOR
     *                                     +VSELEM*FACTUP*FACTOR
C     
                  VPOTDN(N1NUMB+1,N2NUMB+1)=VPOTDN(N1NUMB+1,N2NUMB+1)
     *                                     +VCELEM*FACTOR
     *                                     +VSELEM*FACTDN*FACTOR
C
               END DO ! End of the n2 do-loop 
C            
            END DO ! End of the n1 do-loop  
C_______________________________________________________________________
C             
         END DO ! Ending Gauss-Laguerre integration for the potentials
C_______________________________________________________________________
C   
C        Adding the kinetic energy matrix elements already calculated 
C        in KINMAT to the matrix elements of the potential       
C
         DO N1NUMB=0,(NSHELL_NEUTRS-LQNUMB)/2
            DO N2NUMB=0,(NSHELL_NEUTRS-LQNUMB)/2  
C         
               EKINET=ENEKIN_NEUTRS(INUCLI,N1NUMB+1,N2NUMB+1,LQNUMB)
C
               HSORUP(N1NUMB+1,N2NUMB+1)
     *        =
     *         EKINET+VPOTUP(N1NUMB+1,N2NUMB+1)
C               
               HSORDN(N1NUMB+1,N2NUMB+1)
     *        =
     *         EKINET+VPOTDN(N1NUMB+1,N2NUMB+1)
C         
             END DO ! End of the n2 do-loop 
         END DO ! End of the n1 do-loop  
C
C=======================================================================
C      
         LDBASE=(NSHELL_NEUTRS-LQNUMB)/2+1
C
C=======================================================================
C        Symmetry test for the Hamiltonian matrix
C=======================================================================
C
         DO N1NUMB=0,(NSHELL_NEUTRS-LQNUMB)/2
            DO N2NUMB=N1NUMB,(NSHELL_NEUTRS-LQNUMB)/2       
C
               IF (ABS(HSORUP(N1NUMB+1,N2NUMB+1)
     *                -HSORUP(N2NUMB+1,N1NUMB+1)).GT.EPSTSH) THEN
C
                   WRITE(NOUTPT,'(''N1='',I3,3X,''N2='',I3,3X,
     *                            ''HSORUP(N1,N2)='',F30.25,3X,
     *                            ''HSORUP(N2,N1)='',F30.25)')
     *
     *                   N1NUMB,N2NUMB,
     *                   HSORUP(N1NUMB+1,N2NUMB+1),
     *                   HSORUP(N2NUMB+1,N1NUMB+1)
C
                   STOP 'Non hermiticity in matrix HSORUP -NEUTRONS- !'
               END IF
C
               IF (ABS(HSORDN(N1NUMB+1,N2NUMB+1)
     *                -HSORDN(N2NUMB+1,N1NUMB+1)).GT.EPSTSH) THEN
C
                   WRITE(NOUTPT,'(''N1='',I3,3X,''N2='',I3,3X,
     *                            ''HSORDN(N1,N2)='',F30.25,3X,
     *                            ''HSORDN(N2,N1)='',F30.25)')
     *
     *                   N1NUMB,N2NUMB,
     *                   HSORDN(N1NUMB+1,N2NUMB+1),
     *                   HSORDN(N2NUMB+1,N1NUMB+1)
C
                   STOP 'Non-hermiticity in matrix HSORDN -NEUTRONS- !'
               END IF
C
            END DO
         END DO
C
C=======================================================================
C        Diagonalising the real-symmetric hamiltonian matrix
C=======================================================================
C                                                EIGVEC=.TRUE.
         CALL DIAMAT(HSORUP,SPENUP,AUXVEC,NDBASE,LDBASE,.TRUE.)
C_____________________________________________________________________
C
         DO NWSNUM=1,(NSHELL_NEUTRS-LQNUMB)/2+1            
            DO N1NUMB=0,(NSHELL_NEUTRS-LQNUMB)/2
C
               CAUXUP_NEUTRS(N1NUMB,NWSNUM)=HSORUP(N1NUMB+1,NWSNUM)
               
               CMATUP_NEUTRS(LQNUMB,N1NUMB+1,NWSNUM)
     *        =              HSORUP(N1NUMB+1,NWSNUM)
C
               ENERUP_NEUTRS(LQNUMB,NWSNUM)=SPENUP(NWSNUM)
C
            END DO
         END DO
C
         DO I=1,LDBASE
                   NCOUNT=NCOUNT+1
            ENETHE(NCOUNT)=SPENUP(I)
            LAUXIL(NCOUNT)=LQNUMB
            JAUXIL(NCOUNT)=LQNUMB+LQNUMB+1
            NAUXIL(NCOUNT)=I-1
            INDEXS(NCOUNT)=NCOUNT
         END DO    
C_____________________________________________________________________
C
         IF (LQNUMB.NE.0) THEN        ! Diagonalisation for L > 0
C                                                    EIGVEC=.TRUE.
             CALL DIAMAT(HSORDN,SPENDN,AUXVEC,NDBASE,LDBASE,.TRUE.)
C
             DO N1NUMB=0,(NSHELL_NEUTRS-LQNUMB)/2
                DO NWSNUM=1,(NSHELL_NEUTRS-LQNUMB)/2+1
C
                   CAUXDN_NEUTRS(N1NUMB,NWSNUM)=HSORDN(N1NUMB+1,NWSNUM)
                   
                   CMATDN_NEUTRS(LQNUMB,N1NUMB+1,NWSNUM)
     *        =                  HSORDN(N1NUMB+1,NWSNUM)
C
                   ENERDN_NEUTRS(LQNUMB,NWSNUM)=SPENDN(NWSNUM)
C
                END DO
             END DO
C
             DO I=1,LDBASE
                       NCOUNT =NCOUNT+1
                ENETHE(NCOUNT)=SPENDN(I)
                LAUXIL(NCOUNT)=LQNUMB
                JAUXIL(NCOUNT)=LQNUMB+LQNUMB-1
                NAUXIL(NCOUNT)=I-1
                INDEXS(NCOUNT)=NCOUNT
             END DO 
C
         END IF 
C______________________________________________________________________
C      
      END DO  ! LQNUMB   
C
C=======================================================================
C
      DO I=1,NCOUNT
          INDEXS(I)=I
      END DO
      
      CALL ORDHEA(ENETHE,INDEXS,NCOUNT,NDSPEC)
C           
C=======================================================================
C
      LDSING_NEUTRS=NCOUNT
C
      NOCCUP_NEUTRS=0
C
      DO I=1,NCOUNT
C
         INDOLD=INDEXS(I)
C
         NWSSPH_NEUTRS(I)=NAUXIL(INDOLD)
         LWSSPH_NEUTRS(I)=LAUXIL(INDOLD)
         JWSSPH_NEUTRS(I)=JAUXIL(INDOLD)
C
         ENERGY_NEUTRS(I,ITRNUM)=ENETHE(I)
         ENETHE_NEUTRS(I       )=ENETHE(I)
C
         IF (IN_FIX.GT.126 .AND. NOCCUP_NEUTRS.LT.IN_FIX) THEN
C
             NOCCUP_NEUTRS=NOCCUP_NEUTRS+JWSSPH_NEUTRS(I)+1
             LABORD_NEUTRS(INUCLI,NWSSPH_NEUTRS(I),LWSSPH_NEUTRS(I),
     *                                             JWSSPH_NEUTRS(I))=1
C
         END IF
C
      END DO
C
C=======================================================================
C 
C     Verifying the self-consistency condition, while iterating in 
C     the case of the solution-dependent Hamiltonian spin-orbit or
C                                                           tensor
C
      CALL CONSIS_VERIFS(IABORT)
C       
C=======================================================================
C     Below we prepare the spherical labels, e.g. for the plotting 
C     system and for the chi^2 minimisation, and in order to print 
C     the corresponding tables
C
      CALL SPHLAB(IZ_FIX,IN_FIX,EFERMI_PROTON,EFERMI_NEUTRS)
C       
C=======================================================================
C
C     Printing the energies and parameters in the standard output
C
      IF (LOGWRI.GT.0) THEN
          WRITE(NOUTPT,'(80(''#''),/,''#'',T80,''#'',/,
     *             ''#  ITRNUM= '',I2,5X,''IZ_FIX= '',I3,
     *             5X,''IN_FIX= '',I3,T80,''#'',/,''#'',T80,''#'',/,
     *             ''#  Theoretical PROTON Spectrum'',T80,''#'',/,
     *             ''#'',78X,''#'',/,80(''#''),/,''#'',T80,''#'')')
     *                                        ITRNUM,IZ_FIX,IN_FIX
          I_P=0
C
          DO INDEXX=1,LDSING_PROTON
C
             I_P=I_P+JWSSPH_PROTON(INDEXX)+1
C
             WRITE(NOUTPT,'(''# '',I4,'')'',3X,F20.13,6X,A,6X,I4,27X,
     *                                                         ''#'')')
     *         INDEXX,ENETHE_PROTON(INDEXX),LABTHE_PROTON(INDEXX),I_P
C
               IF (ENETHE_PROTON(INDEXX+1).GT.SPENER_MXPRIN) THEN
                  GO TO 1
               END IF
          END DO
C
          WRITE(NOUTPT,'(''#'',78X,''#'',/,80(''#''),/)')
C
   1      CONTINUE
C_______________________________________________________________________
C      
          WRITE(NOUTPT,'(80(''#''),/,''#'',T80,''#'',/,
     *             ''#  ITRNUM= '',I2,5X,''IZ_FIX= '',I3,
     *             5X,''IN_FIX= '',I3,T80,''#'',/,''#'',T80,''#'',/,
     *             ''#  Theoretical NEUTRON Spectrum'',T80,''#'',/,
     *             ''#'',78X,''#'',/,80(''#''),/,''#'',T80,''#'')')
     *                                        ITRNUM,IZ_FIX,IN_FIX
C
          I_N=0
          DO INDEXX=1,LDSING_NEUTRS
C
             I_N=I_N+JWSSPH_NEUTRS(INDEXX)+1
C
             WRITE(NOUTPT,'(''# '',I4,'')'',3X,F20.13,6X,A,6X,I4,27X,
     *                                                         ''#'')')
     *         INDEXX,ENETHE_NEUTRS(INDEXX),LABTHE_NEUTRS(INDEXX),I_N
C
             IF (ENETHE_NEUTRS(INDEXX+1).GT.SPENER_MXPRIN) THEN
                 GO TO 2
             END IF
          END DO
C
          WRITE(NOUTPT,'(''#'',78X,''#'',/,80(''#''),/)')
C
   2      CONTINUE
C
      END IF
C           
C=======================================================================
C
      DO LQNUMB=NSHELL_PROTON,0,-1
         DO NWSNUM=1,(NSHELL_PROTON-LQNUMB)/2+1
             
            ENERUP_AUXPRO(LQNUMB,NWSNUM)=ENERUP_PROTON(LQNUMB,NWSNUM)
            ENERDN_AUXPRO(LQNUMB,NWSNUM)=ENERDN_PROTON(LQNUMB,NWSNUM)
            
            DO N1NUMB=0,(NSHELL_PROTON-LQNUMB)/2
                
               CMATUP_AUXPRO(LQNUMB,N1NUMB+1,NWSNUM)
     *        =CMATUP_PROTON(LQNUMB,N1NUMB+1,NWSNUM)
               
               CMATDN_AUXPRO(LQNUMB,N1NUMB+1,NWSNUM)
     *        =CMATDN_PROTON(LQNUMB,N1NUMB+1,NWSNUM)
     
            END DO
            
         END DO
      END DO 
C           
C=======================================================================
C
      DO LQNUMB=NSHELL_NEUTRS,0,-1
         DO NWSNUM=1,(NSHELL_NEUTRS-LQNUMB)/2+1
             
            ENERUP_AUXNEU(LQNUMB,NWSNUM)=ENERUP_NEUTRS(LQNUMB,NWSNUM)
            ENERDN_AUXNEU(LQNUMB,NWSNUM)=ENERDN_NEUTRS(LQNUMB,NWSNUM)
            
            DO N1NUMB=0,(NSHELL_NEUTRS-LQNUMB)/2
                
               CMATUP_AUXNEU(LQNUMB,N1NUMB+1,NWSNUM)
     *        =CMATUP_NEUTRS(LQNUMB,N1NUMB+1,NWSNUM)
               
               CMATDN_AUXNEU(LQNUMB,N1NUMB+1,NWSNUM)
     *        =CMATDN_NEUTRS(LQNUMB,N1NUMB+1,NWSNUM)
     
            END DO
            
         END DO
      END DO   
C
C=======================================================================
C     Calling the pairing subroutine (optional)
C=======================================================================
C
      IF (IF_PAI.EQ.1) THEN 
C          
          ISOSPI=1
C          
          CALL PAIDEL(ENERUP_AUXPRO,ENERDN_AUXPRO,NSHELL_PROTON,
     *                ENEMAX,VCOEUP_PROTON,VCOEDN_PROTON,NDIM_L,
     *                                            NDBASE,INUCLI)
C          
          ENEMAX_BCSPAI(INUCLI,1)=ENEMAX
C_______________________________________________________________________
C         
          ISOSPI=0
C          
          CALL PAIDEL(ENERUP_AUXNEU,ENERDN_AUXNEU,NSHELL_NEUTRS,
     *                ENEMAX,VCOEUP_NEUTRS,VCOEDN_NEUTRS,NDIM_L,
     *                                            NDBASE,INUCLI)
C          
          ENEMAX_BCSPAI(INUCLI,2)=ENEMAX
C      
      END IF
C
C=======================================================================
C 
C     Checking the proton-density integral
C  
      XNUMBZ=0.0
C
      DO LQNUMB=NSHELL_PROTON,0,-1
         DO I_PNTS=1,NGAUSS
C
            ZNODES=RNODES(I_PNTS,LQNUMB)
C
            XNUMBZ=XNUMBZ+DENSIT_LAGUER(INUCLI,1,ZNODES,I_PNTS,LQNUMB)
     *                   *ZNODES*RWEIGH(I_PNTS,LQNUMB)*2
         END DO
      END DO
C______________________________________________________________________
C  
C     Pairing tests and checks
C        
      IF (IF_PAI.EQ.1) THEN
C      
          V2_SUM=0.0
C
          DO LQNUMB=NSHELL_PROTON,0,-1
C
             JQNUMB=2*LQNUMB+1
C
             DO NWSNUM=1,(NSHELL_PROTON-LQNUMB)/2+1
C
                XNUMBZ_STATUP=0.0
C
                DO I_PNTS=1,NGAUSS ! Integration over "z"
C
                   ZNODES=RNODES(I_PNTS,LQNUMB)
C
                   XNUMBZ_STATUP=XNUMBZ_STATUP
     *                          +DENSIT_STATES(INUCLI,1,ZNODES,I_PNTS,
     *                                           LQNUMB,NWSNUM,JQNUMB)
     *                          *ZNODES*RWEIGH(I_PNTS,LQNUMB)*2
C
                END DO
C
                V2_SUM=V2_SUM+V2_PUP_VECTOR(LQNUMB,NWSNUM)*(JQNUMB+1)
C
             END DO
C
          END DO
C
          DO LQNUMB=NSHELL_PROTON,0,-1
C
             JQNUMB=2*LQNUMB-1
C
             DO NWSNUM=1,(NSHELL_PROTON-LQNUMB)/2+1
C
                XNUMBZ_STATDN=0.0
C
                DO I_PNTS=1,NGAUSS
C
                   ZNODES=RNODES(I_PNTS,LQNUMB)
C
                   XNUMBZ_STATDN=XNUMBZ_STATDN
     *                          +DENSIT_STATES(INUCLI,1,ZNODES,I_PNTS,
     *                                           LQNUMB,NWSNUM,JQNUMB)
     *                          *ZNODES*RWEIGH(I_PNTS,LQNUMB)*2
C
                END DO
C
                V2_SUM=V2_SUM+V2_PDW_VECTOR(LQNUMB,NWSNUM)*(JQNUMB+1)
C
             END DO
C
          END DO
C
      END IF
C______________________________________________________________________
C
      IF (LOGWRI.GT.4) THEN
          WRITE(LOGFIL,'(18X,''ITRNUM= '',I2,'' XNUMBZ='',F20.13,
     *                                       '' IZ_FIX='',I2)')
     *                         ITRNUM,XNUMBZ,IZ_FIX
      END IF
C
      IF (ISCREN.GT.5) THEN
          WRITE(0,'(''ITRNUM= '',I2)')ITRNUM
          WRITE(0,'(''In HAMMAT: XNUMBZ_PROTON='',F20.13)')XNUMBZ
          WRITE(0,'(11X,''IZ_FIX='',I6)')IZ_FIX
      END IF
C
      IF (ABS(XNUMBZ-IZ_FIX).GT.1.E-4) THEN
          DO I=1,15
             WRITE(LOGFIL,'(45X,''CHECK WHAT IS WRONG WITH THE '',
     *                          ''PROTON DENSITY'')')
          END DO
CID       STOP 'Stop in HAMMAT: Proton density not correct!'
      END IF
C
C=======================================================================
C 
C      Checking the neutron-density integral
C       
       XNUMBN=0.0
C       
       DO LQNUMB=NSHELL_PROTON,0,-1
          DO I_PNTS=1,NGAUSS
C           
             ZNODES=RNODES(I_PNTS,LQNUMB)
C                
             XNUMBN=XNUMBN+DENSIT_LAGUER(INUCLI,0,ZNODES,I_PNTS,LQNUMB)
     *                     *ZNODES*RWEIGH(I_PNTS,LQNUMB)*2
          END DO
      END DO
C_______________________________________________________________________
C
      IF (LOGWRI.GT.4) THEN
          WRITE(LOGFIL,'(18X,''ITRNUM= '',I2,'' XNUMBN='',F20.13,
     *                                       '' IN_FIX='',I2)')
     *                         ITRNUM,XNUMBN,IN_FIX
      END IF
C
      IF (ISCREN.GT.5) THEN
          WRITE(0,'(''ITRNUM= '',I2)')ITRNUM
          WRITE(0,'(''In HAMMAT: XNUMBN='',F20.13)')XNUMBN
          WRITE(0,'(11X,''IN_FIX='',I6)')IN_FIX
      END IF
C
      IF (ABS(XNUMBN-IN_FIX).GT.1.E-4) THEN
          DO I=1,15
             WRITE(LOGFIL,'(45X,''CHECK WHAT IS WRONG WITH THE '',
     *                          ''NEUTRON DENSITY'')')
          END DO
CID       STOP 'Stop in HAMMAT: Proton density not correct!'
      END IF
C       
C=======================================================================
C
C      CALL PRTNG(INUCLI,1)  ! ID: Creating the Y/N file     
C       
C=======================================================================
C
      IF (LOGWRI.GT.4) THEN
          WRITE(LOGFIL,'(15X,''Exiting  HAMMAT'')')
      END IF   
C    
C=======================================================================
C
      CALL CPUTIM('HAMMAT',0)   
C
C=======================================================================
C      
      RETURN
      END   
C
C=======================================================================     
C=======================================================================     
C      
      SUBROUTINE WSSTAN(INUCLI,
     *                  POLYNL,RNODES,RWEIGH,NGAUSS,NSHELL_PROTON,
     *                  NDGAUS,NDIM_N,NDIM_L,NDBASE,NSHELL_NEUTRS,
     *                  HSORUP,HSORDN,SPENUP,SPENDN,AUXVEC,IABORT,
     *                                CMATUP_PROTON,CMATDN_PROTON,
     *                                CMATUP_NEUTRS,CMATDN_NEUTRS,
     *                                ENERUP_PROTON,ENERDN_PROTON,
     *                                ENERUP_NEUTRS,ENERDN_NEUTRS,
     *                                ENEKIN_PROTON,ENEKIN_NEUTRS)
C
      INCLUDE  'MATDIM/NDSPEC.f'
      INCLUDE  'MATDIM/NITERA.f'
      INCLUDE  'MATDIM/NDNUCL.f'
      INCLUDE  'MATDIM/NDPARS.f'
CID  
C      REAL*16   
C     *          POLYNL,RNODES,RWEIGH
CID
      DIMENSION
     *          POLYNL(1:NDGAUS,0:NDIM_N,0:NDIM_L)
      DIMENSION
     *          INDEXS(1:NDSPEC)
      DIMENSION
     *          LAUXIL(1:NDSPEC),JAUXIL(1:NDSPEC),
     *                           NAUXIL(1:NDSPEC)
      DIMENSION
     *          ENEKIN_PROTON(1:NDNUCL,1:NDBASE,1:NDBASE,0:NDIM_L),
     *          ENEKIN_NEUTRS(1:NDNUCL,1:NDBASE,1:NDBASE,0:NDIM_L)
      DIMENSION
     *          CMATUP_PROTON(0:NDIM_L,1:NDBASE,1:NDBASE),
     *          CMATDN_PROTON(0:NDIM_L,1:NDBASE,1:NDBASE),
     *
     *          CMATUP_NEUTRS(0:NDIM_L,1:NDBASE,1:NDBASE),
     *          CMATDN_NEUTRS(0:NDIM_L,1:NDBASE,1:NDBASE)
      DIMENSION
     *          ENERUP_PROTON(0:NDIM_L,1:NDBASE),
     *          ENERDN_PROTON(0:NDIM_L,1:NDBASE),
     *
     *          ENERUP_NEUTRS(0:NDIM_L,1:NDBASE),
     *          ENERDN_NEUTRS(0:NDIM_L,1:NDBASE)
      DIMENSION
     *          HSORUP(1:NDBASE,1:NDBASE),
     *          HSORDN(1:NDBASE,1:NDBASE)
      DIMENSION
     *          RNODES(1:NDGAUS,0:NDIM_L),
     *          RWEIGH(1:NDGAUS,0:NDIM_L)
      DIMENSION
     *          SPENUP(1:NDBASE),
     *          SPENDN(1:NDBASE),     
     *          AUXVEC(1:NDBASE)     
      DIMENSION
     *          ENETHE(1:NDSPEC)
C     
      COMMON
     *       /POTPOT/ PARPOT(1:NDPARS) 
      COMMON
     *       /AUXPOT/ VKACEN_PROTON,RKACEN_PROTON,AKACEN_PROTON,
     *                VKASOR_PROTON,RKASOR_PROTON,AKASOR_PROTON,
     *                VKAEFM_PROTON,RKAEFM_PROTON,AKAEFM_PROTON,
     *                                            V0UNIT_PROTON,
     *                VKACEN_NEUTRS,RKACEN_NEUTRS,AKACEN_NEUTRS,
     *                VKASOR_NEUTRS,RKASOR_NEUTRS,AKASOR_NEUTRS,
     *                VKAEFM_NEUTRS,RKAEFM_NEUTRS,AKAEFM_NEUTRS,
     *                                            V0UNIT_NEUTRS,
     *                                                   RKACOU,
     *                              ALAMPP,ALAMPN,ALAMNP,ALAMNN,
     *                              TLAMPP,TLAMPN,TLAMNP,TLAMNN,
     *                              CLAMPP,CLAMPN,CLAMNP,CLAMNN
      COMMON
     *       /KAPPAR/ V0CENT_KAPPAR,XK_V0C_KAPPAR,
     *                R0CENT_KAPPAR,XK_R0C_KAPPAR,
     *                A0CENT_KAPPAR,XK_A0C_KAPPAR,
     *                V0SORB_KAPPAR,XK_V0S_KAPPAR,
     *                R0SORB_KAPPAR,XK_R0S_KAPPAR,
     *                A0SORB_KAPPAR,XK_A0S_KAPPAR
      COMMON
     *       /NUCLEU/ IZ_FIX,IN_FIX,ISOSPI    
      COMMON
     *       /HBAR_V/ HBAR_C
      COMMON
     *       /DATA_V/ AOSCIL_PROTON(1:NDNUCL),
     *                AOSCIL_NEUTRS(1:NDNUCL),
     *                HOMEGA(1:NDNUCL)      
      COMMON
     *       /DIMSPH/ LDSING_PROTON,LEVTHE_PROTON,
     *                LEVEXP_PROTON(1:NDNUCL),
     *                LDSING_NEUTRS,LEVTHE_NEUTRS,
     *                LEVEXP_NEUTRS(1:NDNUCL)
      COMMON
     *       /SPECTR/ ENETHE_PROTON(1:NDSPEC),
     *                ENETHE_NEUTRS(1:NDSPEC)
      COMMON
     *       /QUANTN/ NWSSPH_PROTON(1:NDSPEC),
     *                LWSSPH_PROTON(1:NDSPEC),
     *                JWSSPH_PROTON(1:NDSPEC),
     *
     *                NWSSPH_NEUTRS(1:NDSPEC),
     *                LWSSPH_NEUTRS(1:NDSPEC),
     *                JWSSPH_NEUTRS(1:NDSPEC)
      COMMON
     *       /ENRGYY/ ENERGY_PROTON(1:NDSPEC,0:NITERA), 
     *                ENERGY_NEUTRS(1:NDSPEC,0:NITERA)
      COMMON
     *       /PRFILE/ ISCREN,LOGFIL,LOGCHP,LOGCHN,LOGHIP,LOGHIN,
     *                LOGPRO,LOGNEU,LOGWRI,LSCREN,IRESUL,ICONVE,
     *                NOUTPT,LOGAUX,LOGENE,LOGRAD,LOGGRD,LOGMSH,
     *                                                   LOGBIS
      DATA
     *        EPSTSH / 1.0000E-6 /
      DATA
     *        ALFCON / 137.03602 /  
C
C=======================================================================     
C=======================================================================     
C
C          -> BASED ON FOLLOWING SUBROUTINES: HAMMAT, INTROD
C
C          -> USED IF ONE ONLY WANTS THE STANDARD WS HAMILTONIAN !
C
C=======================================================================     
C=======================================================================     
C
      IABORT=0
C
C=======================================================================     
C                    Universal Woods-Saxon parameters
C=======================================================================     
C
      X_MASS=938.9059
      HBAR_C=197.3289
C
      V0CENT=-49.60            
C
      AKACEN=0.7
      AKASOR=0.7
C
      XKAPPA=0.860000
C
      R0CENT_PROTON=1.275
      R0CENT_NEUTRS=1.347
C
      R0SORB_PROTON=1.320
      R0SORB_NEUTRS=1.310
C
      V0SORB_PROTON=36.0 ! Unitless
      V0SORB_NEUTRS=35.0 ! Unitless
C
      R0COUL=1.275
C
C=======================================================================     
C
      AOSCIL=HBAR_C**2/HOMEGA(INUCLI)/X_MASS
      AOSCIL=SQRT(AOSCIL)
C
      A_MASS=IN_FIX+IZ_FIX
C
      RKACOU=R0COUL*A_MASS**(1./3.)
C
      VKASOR=1.0
C
C=======================================================================
C
      V0CENT_KAPPAR=V0CENT
      XK_V0C_KAPPAR=XKAPPA
      R0CENT_KAPPAR=0.0000
      XK_R0C_KAPPAR=0.0000
      A0CENT_KAPPAR=0.0000
      XK_A0C_KAPPAR=0.0000
C
      V0SORB_KAPPAR=0.0000
      XK_V0S_KAPPAR=0.0000
      R0SORB_KAPPAR=0.0000
      XK_R0S_KAPPAR=0.0000
      A0SORB_KAPPAR=0.0000
      XK_A0S_KAPPAR=0.0000
C
C=======================================================================
C=======================================================================
C              Beginning with the algorithm for the protons
C=======================================================================
C=======================================================================
C
      WRITE(NOUTPT,'(/)')
C
      NCOUNT=0
      ISOSPI=1
C
C=======================================================================
C
      VKACEN=V0CENT*(1+XKAPPA*(IN_FIX-IZ_FIX)/A_MASS)
C
      RKACEN=R0CENT_PROTON*A_MASS**(1./3.)
C
      RKASOR=R0SORB_PROTON*A_MASS**(1./3.)
C
      V0UNIT=-0.25*V0SORB_PROTON*VKACEN*(HBAR_C/X_MASS)**2
     *                             /((A_MASS-1)/A_MASS)**2
C
C=======================================================================
C
      VKACEN_PROTON=VKACEN
      RKACEN_PROTON=RKACEN
      AKACEN_PROTON=AKACEN
C
      VKASOR_PROTON=V0UNIT
      RKASOR_PROTON=RKASOR
      AKASOR_PROTON=AKASOR
C
C=======================================================================
C     Below, we use the convention involving Pauli matrices
C     rather then  s = 1/2*sigma;  this implies eliminating
C     the factor of 1/2
C     Spin-Orbit Interaction: L * \sigma
C=======================================================================
C
      DO LQNUMB=NSHELL_PROTON,0,-1
C
C-----------------------------------------------------------------------
C        Spin-orbit parallel and anti-parallel configurations
C-----------------------------------------------------------------------
C
         FACTUP=+LQNUMB
C         
         FACTDN=-(LQNUMB+1)
C
C-----------------------------------------------------------------------
C        Here we are going to construct the H-sub-block at given L
C-----------------------------------------------------------------------
C
         DO N1NUMB=0,(NSHELL_PROTON-LQNUMB)/2
            SPENUP(N1NUMB+1)=0.0
            SPENDN(N1NUMB+1)=0.0
            DO N2NUMB=0,(NSHELL_PROTON-LQNUMB)/2       
               HSORUP(N1NUMB+1,N2NUMB+1)=0.0
               HSORDN(N1NUMB+1,N2NUMB+1)=0.0
            END DO
         END DO
C
         DO N1NUMB=0,(NSHELL_PROTON-LQNUMB)/2
            DO N2NUMB=0,(NSHELL_PROTON-LQNUMB)/2       
C       
               VCELEM=0.0
               VSELEM=0.0
C       
               DO I_PNTS=1,NGAUSS
C
C-----------------------------------------------------------------------
C                 We are using the weight factor e^{-z} z^{l-1/2}
C                 and consequently the potential is multiplied by "z"
C-----------------------------------------------------------------------
C
                  ZNODES=RNODES(I_PNTS,LQNUMB)
C
                  FACTOR=ZNODES*RWEIGH(I_PNTS,LQNUMB)
     *                  *POLYNL(I_PNTS,N1NUMB,LQNUMB)
     *                  *POLYNL(I_PNTS,N2NUMB,LQNUMB)
C
C=======================================================================
C                 Following structure taken from function VCENTR(ZNODES)
C=======================================================================
C 
                  RADIUS=AOSCIL*SQRT(ZNODES) ! Radius in Fermi
C
                  VCENWS=VKACEN/(1.0+EXP((RADIUS-RKACEN)/AKACEN))
C            
                  IF (RADIUS.LE.RKACOU) THEN
                      V_COUL=HBAR_C*(IZ_FIX-1)/ALFCON/RKACOU
     *                      *(1.500-0.5000000*(RADIUS/RKACOU)**2)   
                  ELSE
                      V_COUL=HBAR_C*(IZ_FIX-1)/(ALFCON*RADIUS)
                  END IF
C
                  VCELEM=VCELEM+(VCENWS+V_COUL)*FACTOR
C
C=======================================================================
C                 Following structure taken from function V_SORB(ZNODES)
C=======================================================================
C      
                  VSOELM=-VKASOR
     *		        *V0UNIT/AKASOR*EXP((RADIUS-RKASOR)/AKASOR)
     *                  /(1.000+EXP((RADIUS-RKASOR)/AKASOR))**2
     *                  /RADIUS
C
                  VSELEM=VSELEM+VSOELM*FACTOR
C
C=======================================================================
C
               END DO
C
C=======================================================================
C
               EKINET=ENEKIN_PROTON(INUCLI,N1NUMB+1,N2NUMB+1,LQNUMB)
C
C=======================================================================
C
               HSORUP(N1NUMB+1,N2NUMB+1)=EKINET+VCELEM+VSELEM*FACTUP
C
               HSORDN(N1NUMB+1,N2NUMB+1)=EKINET+VCELEM+VSELEM*FACTDN
C
C=======================================================================
C
            END DO
         END DO
C
C=======================================================================
C      
         LDBASE=(NSHELL_PROTON-LQNUMB)/2+1
C
C=======================================================================
C        Testing hermiticity of proton up and down blocs
C=======================================================================
C
         DO N1NUMB=0,(NSHELL_PROTON-LQNUMB)/2
            DO N2NUMB=N1NUMB,(NSHELL_PROTON-LQNUMB)/2       
C
               IF (ABS(HSORUP(N1NUMB+1,N2NUMB+1)
     *                -HSORUP(N2NUMB+1,N1NUMB+1)).GT.EPSTSH) THEN
C
                  WRITE(NOUTPT,'(''N1='',I3,3X,''N2='',I3,3X,
     *                      ''HSORUP(N1,N2)='',F30.25,3X,
     *                      ''HSORUP(N2,N1)='',F30.25)')
     *
     *                  N1NUMB,N2NUMB,
     *                  HSORUP(N1NUMB+1,N2NUMB+1),
     *                  HSORUP(N2NUMB+1,N1NUMB+1)
C
                  STOP 'Non hermiticity in HSORUP -PROTONS- WSSTAN'
C
               END IF
C
               IF (ABS(HSORDN(N1NUMB+1,N2NUMB+1)
     *                -HSORDN(N2NUMB+1,N1NUMB+1)).GT.EPSTSH) THEN
C
                  WRITE(NOUTPT,'(''N1='',I3,3X,''N2='',I3,3X,
     *                      ''HSORDN(N1,N2)='',F30.25,3X,
     *                      ''HSORDN(N2,N1)='',F30.25)')
     *
     *                  N1NUMB,N2NUMB,
     *                  HSORDN(N1NUMB+1,N2NUMB+1),
     *                  HSORDN(N2NUMB+1,N1NUMB+1)
C
                  STOP 'Non hermiticity in HSORDN -PROTONS- WSSTAN'
C
               END IF
C
            END DO
         END DO
C
C=======================================================================
C        Diagonalizing the proton spin-up matrix (l parralel s)
C=======================================================================
C
         CALL DIAMAT(HSORUP,SPENUP,AUXVEC,NDBASE,LDBASE,.TRUE.)
C
C=======================================================================
C
         DO N1NUMB=0,(NSHELL_PROTON-LQNUMB)/2
            DO N2NUMB=0,(NSHELL_PROTON-LQNUMB)/2 
C
               CMATUP_PROTON(LQNUMB,N1NUMB+1,N2NUMB+1) 
     *         =             HSORUP(N1NUMB+1,N2NUMB+1)
C
               ENERUP_PROTON(LQNUMB,N1NUMB+1)=SPENUP(N1NUMB+1)
C
            END DO
         END DO
C
         DO I=1,LDBASE
            NCOUNT=NCOUNT+1
            ENETHE(NCOUNT)=SPENUP(I)
            LAUXIL(NCOUNT)=LQNUMB
            JAUXIL(NCOUNT)=LQNUMB+LQNUMB+1
            NAUXIL(NCOUNT)=I-1
            INDEXS(NCOUNT)=NCOUNT
         END DO    
C
C=======================================================================
C        Diagonalizing the proton spin-down matrix (l anti-parralel s)
C=======================================================================
C
         IF (LQNUMB.NE.0) THEN
C
C=======================================================================
C
             CALL DIAMAT(HSORDN,SPENDN,AUXVEC,NDBASE,LDBASE,.TRUE.)
C
C=======================================================================
C
             DO N1NUMB=0,(NSHELL_PROTON-LQNUMB)/2
                DO N2NUMB=0,(NSHELL_PROTON-LQNUMB)/2 
C
                   CMATDN_PROTON(LQNUMB,N1NUMB+1,N2NUMB+1) 
     *             =             HSORDN(N1NUMB+1,N2NUMB+1)
C
                   ENERDN_PROTON(LQNUMB,N1NUMB+1)=SPENDN(N1NUMB+1)
C
                END DO
             END DO
C
             DO I=1,LDBASE
                NCOUNT=NCOUNT+1
                ENETHE(NCOUNT)=SPENDN(I)
                LAUXIL(NCOUNT)=LQNUMB
                JAUXIL(NCOUNT)=LQNUMB+LQNUMB-1
                NAUXIL(NCOUNT)=I-1
                INDEXS(NCOUNT)=NCOUNT
             END DO 
C
         END IF             
C
      END DO            
C
C=======================================================================
C     Ordering the mixed spin-up/spin-down spectra
C=======================================================================
C
      CALL ORDHEA(ENETHE,INDEXS,NCOUNT,NDSPEC)
C           
C=======================================================================
C
      LDSING_PROTON=NCOUNT
C
      DO I=1,NCOUNT
C
         INDOLD=INDEXS(I)
C
         NWSSPH_PROTON(I)=NAUXIL(INDOLD)
         LWSSPH_PROTON(I)=LAUXIL(INDOLD)
         JWSSPH_PROTON(I)=JAUXIL(INDOLD)
C
         ENERGY_PROTON(I,0)=ENETHE(I)
         ENETHE_PROTON(I  )=ENETHE(I)
C
      END DO
C
C=======================================================================
C=======================================================================
C            Repeating the same algorithm for the neutrons
C=======================================================================
C=======================================================================
C
      WRITE(NOUTPT,'(/)')
C
      NCOUNT=0
      ISOSPI=0
C
C=======================================================================
C
      VKACEN=V0CENT*(1-XKAPPA*(IN_FIX-IZ_FIX)/A_MASS)
C
      RKACEN=R0CENT_NEUTRS*A_MASS**(1./3.)
C
      RKASOR=R0SORB_NEUTRS*A_MASS**(1./3.)
C
      V0UNIT=-0.25*V0SORB_NEUTRS*VKACEN*(HBAR_C/X_MASS)**2
     *                             /((A_MASS-1)/A_MASS)**2
C
C=======================================================================
C
      VKACEN_NEUTRS=VKACEN
      RKACEN_NEUTRS=RKACEN
      AKACEN_NEUTRS=AKACEN
C
      VKASOR_NEUTRS=V0UNIT
      RKASOR_NEUTRS=RKASOR
      AKASOR_NEUTRS=AKASOR
C
C-----------------------------------------------------------------------
C     Below, we use the convention involving Pauli matrices
C     rather then  s = 1/2*sigma;  this implies eliminating
C     the factor of 1/2
C     Spin-Orbit Interaction: L * \sigma
C-----------------------------------------------------------------------
C
      DO LQNUMB=NSHELL_NEUTRS,0,-1
C
C-----------------------------------------------------------------------
C        Spin-orbit parallel configuration
C-----------------------------------------------------------------------
C
         FACTUP=     LQNUMB
C
C-----------------------------------------------------------------------
C        Spin-orbit anti-parallel configuration
C-----------------------------------------------------------------------
C
         FACTDN=-    (LQNUMB+1)
C
C-----------------------------------------------------------------------
C        Here we are going to construct the H-sub-block at given L
C-----------------------------------------------------------------------
C
         DO N1NUMB=0,(NSHELL_NEUTRS-LQNUMB)/2
            SPENUP(N1NUMB+1)=0.0
            SPENDN(N1NUMB+1)=0.0
            DO N2NUMB=0,(NSHELL_NEUTRS-LQNUMB)/2       
               HSORUP(N1NUMB+1,N2NUMB+1)=0.0
               HSORDN(N1NUMB+1,N2NUMB+1)=0.0
            END DO
         END DO
C
         DO N1NUMB=0,(NSHELL_NEUTRS-LQNUMB)/2
            DO N2NUMB=0,(NSHELL_NEUTRS-LQNUMB)/2       
C       
               VCELEM=0.0
               VSELEM=0.0
C       
               DO I_PNTS=1,NGAUSS
C
C-----------------------------------------------------------------------
C                 We are using the weight factor e^{-z} z^{l-1/2}
C                 and consequently the potential is multiplied by "z"
C-----------------------------------------------------------------------
C
                  ZNODES=RNODES(I_PNTS,LQNUMB)
C
                  FACTOR=ZNODES*RWEIGH(I_PNTS,LQNUMB)
     *                  *POLYNL(I_PNTS,N1NUMB,LQNUMB)
     *                  *POLYNL(I_PNTS,N2NUMB,LQNUMB)
C
C=======================================================================
C                 Following structure taken from function VCENTR(ZNODES)
C=======================================================================
C 
                  RADIUS=AOSCIL*SQRT(ZNODES) ! Radius in Fermi
C
                  VCENWS=VKACEN/(1.0+EXP((RADIUS-RKACEN)/AKACEN))
C
                  VCELEM=VCELEM+VCENWS*FACTOR
C
C=======================================================================
C                 Following structure taken from function V_SORB(ZNODES)
C=======================================================================
C      
                  VSOELM=-VKASOR
     *		        *V0UNIT/AKASOR*EXP((RADIUS-RKASOR)/AKASOR)
     *                  /(1.000+EXP((RADIUS-RKASOR)/AKASOR))**2
     *                  /RADIUS
C
                  VSELEM=VSELEM+VSOELM*FACTOR
C
C=======================================================================
C                 
               END DO
C
C=======================================================================
C
               EKINET=ENEKIN_NEUTRS(INUCLI,N1NUMB+1,N2NUMB+1,LQNUMB)
C
C=======================================================================
C
               HSORUP(N1NUMB+1,N2NUMB+1)=EKINET+VCELEM+VSELEM*FACTUP
               HSORDN(N1NUMB+1,N2NUMB+1)=EKINET+VCELEM+VSELEM*FACTDN
C
C=======================================================================
C
            END DO
         END DO
C      
         LDBASE=(NSHELL_NEUTRS-LQNUMB)/2+1
C
C=======================================================================
C        Testing hermiticity of neutron up and down blocs
C=======================================================================
C
         DO N1NUMB=0,(NSHELL_PROTON-LQNUMB)/2
            DO N2NUMB=N1NUMB,(NSHELL_PROTON-LQNUMB)/2       
C
               IF (ABS(HSORUP(N1NUMB+1,N2NUMB+1)
     *                -HSORUP(N2NUMB+1,N1NUMB+1)).GT.EPSTSH) THEN
C
                  WRITE(NOUTPT,'(''N1='',I3,3X,''N2='',I3,3X,
     *                      ''HSORUP(N1,N2)='',F30.25,3X,
     *                      ''HSORUP(N2,N1)='',F30.25)')
     *
     *                  N1NUMB,N2NUMB,
     *                  HSORUP(N1NUMB+1,N2NUMB+1),
     *                  HSORUP(N2NUMB+1,N1NUMB+1)
C
                  STOP 'NON HERMITICITY IN MATRIX HSORUP -NEUTRONS- !'
C
               END IF
C
               IF (ABS(HSORDN(N1NUMB+1,N2NUMB+1)
     *                -HSORDN(N2NUMB+1,N1NUMB+1)).GT.EPSTSH) THEN
C
                  WRITE(NOUTPT,'(''N1='',I3,3X,''N2='',I3,3X,
     *                      ''HSORDN(N1,N2)='',F30.25,3X,
     *                      ''HSORDN(N2,N1)='',F30.25)')
     *
     *                  N1NUMB,N2NUMB,
     *                  HSORDN(N1NUMB+1,N2NUMB+1),
     *                  HSORDN(N2NUMB+1,N1NUMB+1)
C
                  STOP 'NON HERMITICITY IN MATRIX HSORDN -NEUTRONS- !'
C
               END IF
C
            END DO
         END DO
C
C=======================================================================
C        Diagonalising the neutron spin-up matrix (l parallel s)
C=======================================================================
C
         CALL DIAMAT(HSORUP,SPENUP,AUXVEC,NDBASE,LDBASE,.TRUE.)
C
C=======================================================================
C
         DO N1NUMB=0,(NSHELL_NEUTRS-LQNUMB)/2
            DO N2NUMB=0,(NSHELL_NEUTRS-LQNUMB)/2 
C
               CMATUP_NEUTRS(LQNUMB,N1NUMB+1,N2NUMB+1) 
     *         =             HSORUP(N1NUMB+1,N2NUMB+1)
C
               ENERUP_NEUTRS(LQNUMB,N1NUMB+1)=SPENUP(N1NUMB+1)
C
            END DO
         END DO
C
         DO I=1,LDBASE
            NCOUNT=NCOUNT+1
            ENETHE(NCOUNT)=SPENUP(I)
            LAUXIL(NCOUNT)=LQNUMB
            JAUXIL(NCOUNT)=LQNUMB+LQNUMB+1
            NAUXIL(NCOUNT)=I-1
            INDEXS(NCOUNT)=NCOUNT
         END DO    
C
C=======================================================================
C        Diagonalising the neutron spin-down matrix (l anti-parallel s)
C=======================================================================
C
         IF (LQNUMB.NE.0) THEN
C
C=======================================================================
C
             CALL DIAMAT(HSORDN,SPENDN,AUXVEC,NDBASE,LDBASE,.TRUE.)
C
C=======================================================================
C
             DO N1NUMB=0,(NSHELL_NEUTRS-LQNUMB)/2
                DO N2NUMB=0,(NSHELL_NEUTRS-LQNUMB)/2 
C
                   CMATDN_NEUTRS(LQNUMB,N1NUMB+1,N2NUMB+1) 
     *             =             HSORDN(N1NUMB+1,N2NUMB+1)
C
                   ENERDN_NEUTRS(LQNUMB,N1NUMB+1)=SPENDN(N1NUMB+1)
C
                END DO
             END DO
C
             DO I=1,LDBASE
                NCOUNT=NCOUNT+1
                ENETHE(NCOUNT)=SPENDN(I)
                LAUXIL(NCOUNT)=LQNUMB
                JAUXIL(NCOUNT)=LQNUMB+LQNUMB-1
                NAUXIL(NCOUNT)=I-1
                INDEXS(NCOUNT)=NCOUNT
             END DO 
C
         END IF             
C
      END DO            
C
C=======================================================================
C
      CALL ORDHEA(ENETHE,INDEXS,NCOUNT,NDSPEC)
C           
C=======================================================================
C
      LDSING_NEUTRS=NCOUNT
C
      DO I=1,NCOUNT
C
         INDOLD=INDEXS(I)
C
         NWSSPH_NEUTRS(I)=NAUXIL(INDOLD)
         LWSSPH_NEUTRS(I)=LAUXIL(INDOLD)
         JWSSPH_NEUTRS(I)=JAUXIL(INDOLD)
C
         ENERGY_NEUTRS(I,0)=ENETHE(I)
         ENETHE_NEUTRS(I  )=ENETHE(I)
C
      END DO
C       
C=======================================================================
C     Below we prepare the spherical labels e.g. for the plotting 
C     system and for the chi^2 minimisation, and print the tables 
C=======================================================================
C
      CALL SPHLAB(IZ_FIX,IN_FIX,EFERMI_PROTON,EFERMI_NEUTRS)
C       
C=======================================================================
C
      RETURN
      END
C
C=======================================================================
C=======================================================================
C
      SUBROUTINE ONERUN_WSUNIV
C
      INCLUDE   'MATDIM/NDSPEC.f'
      INCLUDE   'MATDIM/NDPARS.f'
      INCLUDE   'MATDIM/NDNUCL.f'
      INCLUDE   'MATDIM/NDLEXP.f'
      INCLUDE   'MATDIM/NDMESH.f'
      INCLUDE   'MATDIM/ND_RHO.f'
      INCLUDE   'MATDIM/N_NOYX.f'
      INCLUDE   'MATDIM/NDIM_N.f'
      INCLUDE   'MATDIM/NDTITL.f'
C      
      EXTERNAL
     *         DENSIT,DENTHE_INTEGR,DENEXP,DENEXP_INTEGR,
     *         DENTHE_DENEXP
c
      PARAMETER
     *         (NDIM_L=NDIM_N,NDBASE=NDIM_N+1)
C
      CHARACTER
     *          NUCSYM*6,INPSYM*6,TYPCHI*6
      CHARACTER
     *          LABEXP_PROTON*6,LABEXP_NEUTRS*6,
     *          LABTHE_PROTON*6,LABTHE_NEUTRS*6,
     *          LABTHE_PRONUC*6,LABTHE_NEUNUC*6,
     *          LABTHE_PRINTG*6
      CHARACTER
     *          L_CURV_PROTON*6,L_CURV_NEUTRS*6
      CHARACTER
     *          CURNAM*020,L_CURV*100,CURLAB*100,MAIN_T*256,
     *          XTITLE*256,YTITLE*256,TITLES*012,
     *          TITLES_LATEXS*050,NUCNAM_LATEXS*010
C
      DIMENSION
     *          PARPOT_PRINTG(1:NDPARS)
      DIMENSION
     *          RMSTHE_PRONUC(1:NDNUCL),
     *          LEVTHE_PRONUC(1:NDNUCL),
     *          ENETHE_PRONUC(1:NDSPEC,1:NDNUCL),
     *          LABTHE_PRONUC(1:NDSPEC,1:NDNUCL)
      DIMENSION
     *          RMSTHE_NEUNUC(1:NDNUCL),
     *          LEVTHE_NEUNUC(1:NDNUCL),
     *          ENETHE_NEUNUC(1:NDSPEC,1:NDNUCL),
     *          LABTHE_NEUNUC(1:NDSPEC,1:NDNUCL)
      DIMENSION
     *          ENETHE_PRINTG(1:NDSPEC),
     *          LABTHE_PRINTG(1:NDSPEC)
      DIMENSION
     *          EABSAV_PRONUC(1:NDNUCL),
     *          EABSAV_NEUNUC(1:NDNUCL),
     *          RMSVAL_PRONUC(1:NDNUCL),
     *          RMSVAL_NEUNUC(1:NDNUCL),
     *          ERRMAX_PRONUC(1:NDNUCL),
     *          ERRMAX_NEUNUC(1:NDNUCL)
      DIMENSION
     *          X_CURV_PROTON(1:NDMESH,1:NDSPEC),
     *          Y_CURV_PROTON(1:NDMESH,1:NDSPEC),
     *          L_CURV_PROTON(1:NDMESH,1:NDSPEC)
      DIMENSION
     *          X_CURV_NEUTRS(1:NDMESH,1:NDSPEC),
     *          Y_CURV_NEUTRS(1:NDMESH,1:NDSPEC),
     *          L_CURV_NEUTRS(1:NDMESH,1:NDSPEC)
      DIMENSION
     *          ND_MAX(1:NDSPEC),
     *          X_CURV(1:NDMESH,1:NDSPEC),
     *          Y_CURV(1:NDMESH,1:NDSPEC),
     *          L_CURV(1:NDMESH,1:NDSPEC)
      DIMENSION 
     *          ENERUP_PRONUC(0:NDIM_L,1:NDBASE,1:NDNUCL),
     *          ENERDN_PRONUC(0:NDIM_L,1:NDBASE,1:NDNUCL),
     *
     *          ENERUP_NEUNUC(0:NDIM_L,1:NDBASE,1:NDNUCL),
     *          ENERDN_NEUNUC(0:NDIM_L,1:NDBASE,1:NDNUCL)
      DIMENSION  
     *          CMATUP_PRONUC(0:NDIM_L,1:NDBASE,1:NDBASE,1:NDNUCL),
     *          CMATDN_PRONUC(0:NDIM_L,1:NDBASE,1:NDBASE,1:NDNUCL),
     *
     *          CMATUP_NEUNUC(0:NDIM_L,1:NDBASE,1:NDBASE,1:NDNUCL),
     *          CMATDN_NEUNUC(0:NDIM_L,1:NDBASE,1:NDBASE,1:NDNUCL)
C          
      COMMON
     *       /GAUSHE/ NGAUSS,NSHELL_PROTON,NSHELL_NEUTRS
      COMMON
     *       /POTPOT/ PARPOT(1:NDPARS)
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
     *       /DATA_V/ AOSCIL_PROTON(1:NDNUCL),
     *                AOSCIL_NEUTRS(1:NDNUCL),
     *                HOMEGA(1:NDNUCL) 
      COMMON
     *       /EXPDEN/ LABELN(1:N_NOYX,1:2),
     *                R_MESH(1:ND_RHO,1:N_NOYX),
     *                RHOEXP(1:ND_RHO,1:N_NOYX)
      COMMON
     *       /DIMSPH/ LDSING_PROTON,LEVTHE_PROTON,
     *                              LEVEXP_PROTON(1:NDNUCL),
     *                LDSING_NEUTRS,LEVTHE_NEUTRS,
     *                              LEVEXP_NEUTRS(1:NDNUCL)
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
     *       /EUDAUX/ ENERUP_AUXPRO(0:NDIM_L,1:NDBASE),
     *                ENERDN_AUXPRO(0:NDIM_L,1:NDBASE),
     *
     *                ENERUP_AUXNEU(0:NDIM_L,1:NDBASE),
     *                ENERDN_AUXNEU(0:NDIM_L,1:NDBASE)
      COMMON
     *       /CUDAUX/ CMATUP_AUXPRO(0:NDIM_L,1:NDBASE,1:NDBASE),
     *                CMATDN_AUXPRO(0:NDIM_L,1:NDBASE,1:NDBASE),
     *
     *                CMATUP_AUXNEU(0:NDIM_L,1:NDBASE,1:NDBASE),
     *                CMATDN_AUXNEU(0:NDIM_L,1:NDBASE,1:NDBASE)
      COMMON
     *       /MODIFR/ RMSEXP_PROTON(1:NDNUCL),
     *                RMSEXP_NEUTRS(1:NDNUCL),
     *                RMSTHE_PROTON,RMSTHE_NEUTRS,
     *                RMODIF_PROTON,RMODIF_NEUTRS,
     *                              PUSHIN,DELTAR    
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
     *       /NUCLEU/ IZ_FIX,IN_FIX,ISOSPI
      COMMON
     *       /ARGUMN/ IARG_R,IARG_Z
      COMMON
     *       /CHI_PN/ CHISQU_PROTON,CHISQU_NEUTRS
      COMMON
     *       /CHICHO/ VALCHI
     *       /DENSTR/ IFDENS
     *       /TENSOR/ IFTENS
      COMMON
     *       /TENOPT/ ICENTT,ISORBT,ITENSR
     *       /TESFIT/ IFTEST,IFFITS,IFMESH,ISPE_P,IFMOCA
     *       /P_N_LM/ IFPROT,IFNEUT,IFBOTH
     *       /MINITE/ IDEFCN,ITECHI
      COMMON
     *       /FITTED/ IFITED
     *       /MASSIV/ IMASIV_PRODUC
      COMMON
     *       /IFPAIR/ IF_PAI,IFGPAI,IFDPAI
      COMMON
     *       /SIMUSE/ AOSCIL,RADCUT,INUCLI
      COMMON
     *       /CUTOFF/ RADIUS_CUTOFF(1:NDNUCL)
      COMMON 
     *       /CHITAK/ IF_SPE,IF_RAD,IF_GAP,IF_FER,IF_DEN,IF_RHO,IF_INV,
     *                LDCHI2,TYPCHI(1:NDNUCL),ITKCHI(1:NDNUCL)
      COMMON
     *       /PARAMT_TITLES/ TITLES(1:NDTITL),
     *                       TITLES_LATEXS(1:NDTITL),
     *                       NUCNAM_LATEXS(1:NDNUCL)
      COMMON
     *       /AUXPOT/ VKACEN_PROTON,RKACEN_PROTON,AKACEN_PROTON,
     *                VKASOR_PROTON,RKASOR_PROTON,AKASOR_PROTON,
     *                VKAEFM_PROTON,RKAEFM_PROTON,AKAEFM_PROTON,
     *                                            V0UNIT_PROTON,
     *                VKACEN_NEUTRS,RKACEN_NEUTRS,AKACEN_NEUTRS,
     *                VKASOR_NEUTRS,RKASOR_NEUTRS,AKASOR_NEUTRS,
     *                VKAEFM_NEUTRS,RKAEFM_NEUTRS,AKAEFM_NEUTRS,
     *                                            V0UNIT_NEUTRS,
     *                                                   RKACOU,
     *                              ALAMPP,ALAMPN,ALAMNP,ALAMNN,
     *                              TLAMPP,TLAMPN,TLAMNP,TLAMNN,
     *                              CLAMPP,CLAMPN,CLAMNP,CLAMNN
      COMMON
     *       /KAPPAR/ V0CENT_KAPPAR,XK_V0C_KAPPAR,
     *                R0CENT_KAPPAR,XK_R0C_KAPPAR,
     *                A0CENT_KAPPAR,XK_A0C_KAPPAR,
     *                V0SORB_KAPPAR,XK_V0S_KAPPAR,
     *                R0SORB_KAPPAR,XK_R0S_KAPPAR,
     *                A0SORB_KAPPAR,XK_A0S_KAPPAR
C
      COMMON
     *       /PRFILE/ ISCREN,LOGFIL,LOGCHP,LOGCHN,LOGHIP,LOGHIN,
     *                LOGPRO,LOGNEU,LOGWRI,LSCREN,IRESUL,ICONVE,
     *                NOUTPT,LOGAUX,LOGENE,LOGRAD,LOGGRD,LOGMSH,
     *                                                   LOGBIS 
C       
C=======================================================================
C
C     This subroutine calculates the energy levels, densities... for
C     one single run with the parameter values desired for the user.
C
C     It is only activated if IFTEST=1 for the input file
C
C     The desired parameter values come from the input file, too
C       
C=======================================================================
C
      IF (LOGWRI.GT.0) THEN
          WRITE(LOGFIL,'(/,9X,''Entering ONERUN_WSUNIV'',/)')
      END IF
C
C=======================================================================
C
      IFPROT=0
      IFNEUT=0
      IFBOTH=0
C
C=======================================================================
C
C     First, calculating the total weights (protons and neutrons)
C     that will be used as the normalizing constants of the \chi^2
C
      WEIGHT_PROTON=0.0
      WEIGHT_NEUTRS=0.0
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
                 WEIGHT_PROTON=WEINUC_PROTON(INUCLI)
C
                 DEGSUM_NEUTRS=0.0
C
                 DO IEXPER=1,LEVEXP_NEUTRS(INUCLI)
                    DEGAUX=REAL(IDEGEX_NEUTRS(INUCLI,IEXPER))
                    DEGSUM_NEUTRS=DEGSUM_NEUTRS+DEGAUX
                 END DO
C
                 WEIGHT_NEUTRS=WEINUC_NEUTRS(INUCLI)
C
             END IF
C_______________________________________________________________________
C
         END IF
      END DO
C       
C=======================================================================
C
      IMASIV_PRODUC=1
      IRANDO=0
      LDRAND=0
      IEVALU_PRINTG=1
      CHIGRD_PRINTG=999.
      CHISQU_MINIML=999.
C       
C=======================================================================
C
      DO IPARAM=1,NDPARS
         PARPOT_PRINTG(IPARAM)=PARPOT(IPARAM)
      END DO
C
      IF (NUCACT.LT.LDNUCL) RMSGLO_PRINTG=9.9999
C          
      RMSGLO_PROTON=0.
      RMSGLO_NEUTRS=0.
C
      RMSGLO_PROAUX=0.
      RMSGLO_NEUAUX=0.
C
      IDEG_P=0
      IDEG_N=0
C       
C=======================================================================
C
      DO INUCLI=1,LDNUCL
C
         IZ_FIX=NUMB_Z(INUCLI)
         IN_FIX=NUMB_N(INUCLI)
C
C=======================================================================     
C
C        Diagonalising the Hamiltonian for the single-run mode
C
         IF (LOGWRI.GT.4) THEN
             WRITE(LOGFIL,'(12x,''Entering WS_RUN from ONERUN with: '',
     *                          ''IZ_FIX= '',I3,'' and IN_FIX= '',I3)')
     *                            IZ_FIX,              IN_FIX
         END IF
C
         IDEFCN=0
                                                 I_MODE=0
                                                        I_FLAG=1
         CALL WS_UNI(INUCLI,CHISQU_CHOICE,IABORT,I_MODE,I_FLAG,
     *                             CHISQU_PROTON,CHISQU_NEUTRS)
C
         VALCHI=CHISQU_PROTON/WEIGHT_PROTON
     *         +CHISQU_NEUTRS/WEIGHT_NEUTRS 
C
         IF (ITAKNU(INUCLI).EQ.1) CHISQU_MINIML=VALCHI
C_______________________________________________________________________
C     
         RMSTHE_PRONUC(INUCLI)=RMSTHE_PROTON
         LEVTHE_PRONUC(INUCLI)=LEVTHE_PROTON
         
         DO ITHEOR=1,LEVTHE_PROTON
            ENETHE_PRONUC(ITHEOR,INUCLI)=ENETHE_PROTON(ITHEOR)
            LABTHE_PRONUC(ITHEOR,INUCLI)=LABTHE_PROTON(ITHEOR)
         END DO
C_______________________________________________________________________
C         
         RMSTHE_NEUNUC(INUCLI)=RMSTHE_NEUTRS
         LEVTHE_NEUNUC(INUCLI)=LEVTHE_NEUTRS
         
         DO ITHEOR=1,LEVTHE_PROTON
            ENETHE_NEUNUC(ITHEOR,INUCLI)=ENETHE_NEUTRS(ITHEOR)
            LABTHE_NEUNUC(ITHEOR,INUCLI)=LABTHE_NEUTRS(ITHEOR)
         END DO
C_______________________________________________________________________
C
C        Calculating the TOTAL degenerancy of all the nuclei
C
         DO IEXPER=1,LEVEXP_PROTON(INUCLI)
            IDEG_P=IDEG_P+IDEGEX_PROTON(INUCLI,IEXPER)
         END DO
C
         DO IEXPER=1,LEVEXP_NEUTRS(INUCLI)
            IDEG_N=IDEG_N+IDEGEX_NEUTRS(INUCLI,IEXPER)
         END DO
C_______________________________________________________________________
C                  
         EABSAV_PRONUC(INUCLI)=EABSAV_PROTON
         EABSAV_NEUNUC(INUCLI)=EABSAV_NEUTRS
         
         ERRMAX_PRONUC(INUCLI)=ERRMAX_PROTON
         ERRMAX_NEUNUC(INUCLI)=ERRMAX_NEUTRS
         
         RMSVAL_PRONUC(INUCLI)=SQRT(CHIWEI_PROTON)
         RMSVAL_NEUNUC(INUCLI)=SQRT(CHIWEI_NEUTRS)
         
         RMSGLO_PROTON=RMSGLO_PROTON+CHIWEI_PROTON
         RMSGLO_NEUTRS=RMSGLO_NEUTRS+CHIWEI_NEUTRS

         RMSGLO_PROAUX=RMSGLO_PROAUX+CHIDEG_PROTON
         RMSGLO_NEUAUX=RMSGLO_NEUAUX+CHIDEG_NEUTRS
C_______________________________________________________________________
C
         DO LQNUMB=NSHELL_PROTON,0,-1
            DO NWSNUM=1,(NSHELL_PROTON-LQNUMB)/2+1
C
               ENERUP_PRONUC(LQNUMB,NWSNUM,INUCLI)
     *        =
     *         ENERUP_AUXPRO(LQNUMB,NWSNUM)
C
               ENERDN_PRONUC(LQNUMB,NWSNUM,INUCLI)
     *        =
     *         ENERDN_AUXPRO(LQNUMB,NWSNUM)
C
               DO N1NUMB=0,(NSHELL_PROTON-LQNUMB)/2
C
                  CMATUP_PRONUC(LQNUMB,N1NUMB+1,NWSNUM,INUCLI)
     *           =CMATUP_AUXPRO(LQNUMB,N1NUMB+1,NWSNUM)
C
                  CMATDN_PRONUC(LQNUMB,N1NUMB+1,NWSNUM,INUCLI)
     *           =CMATDN_AUXPRO(LQNUMB,N1NUMB+1,NWSNUM)
C
               END DO
C
            END DO
         END DO  
C_______________________________________________________________________
C
         DO LQNUMB=NSHELL_NEUTRS,0,-1
            DO NWSNUM=1,(NSHELL_NEUTRS-LQNUMB)/2+1
C
               ENERUP_NEUNUC(LQNUMB,NWSNUM,INUCLI)
     *        =
     *         ENERUP_AUXNEU(LQNUMB,NWSNUM)
C
               ENERDN_NEUNUC(LQNUMB,NWSNUM,INUCLI)
     *        =
     *         ENERDN_AUXNEU(LQNUMB,NWSNUM)
C
               DO N1NUMB=0,(NSHELL_NEUTRS-LQNUMB)/2
C
                  CMATUP_NEUNUC(LQNUMB,N1NUMB+1,NWSNUM,INUCLI)
     *           =CMATUP_AUXNEU(LQNUMB,N1NUMB+1,NWSNUM)
C
                  CMATDN_NEUNUC(LQNUMB,N1NUMB+1,NWSNUM,INUCLI)
     *           =CMATDN_AUXNEU(LQNUMB,N1NUMB+1,NWSNUM)
C
               END DO
C
            END DO
         END DO   
C_______________________________________________________________________
C          
      END DO !INUCLI
C_______________________________________________________________________
C      
      RMSGLO_PROTON=SQRT(RMSGLO_PROTON)
      RMSGLO_NEUTRS=SQRT(RMSGLO_NEUTRS)
C
      RMSGLO_PROAUX=SQRT(RMSGLO_PROAUX/IDEG_P)
      RMSGLO_NEUAUX=SQRT(RMSGLO_NEUAUX/IDEG_N)
C
C=======================================================================
C
      DO INUCLI=1,LDNUCL
C
         IZ_FIX=NUMB_Z(INUCLI)
         IN_FIX=NUMB_N(INUCLI)
C
C=======================================================================     
C                    Universal Woods-Saxon parameters
C=======================================================================     
C
         X_MASS=938.9059
         HBAR_C=197.3289
C
         V0CENT=-49.60            
C
         AKACEN=0.7
         AKASOR=0.7
C
         XKAPPA=0.860000
C
         R0CENT_PROTON=1.275
         R0CENT_NEUTRS=1.347
C
         R0SORB_PROTON=1.320
         R0SORB_NEUTRS=1.310
C
         V0SORB_PROTON=36.0 ! Unitless
         V0SORB_NEUTRS=35.0 ! Unitless
C
         R0COUL=1.275
C
C=======================================================================     
C
         AOSCIL=HBAR_C**2/HOMEGA(INUCLI)/X_MASS
         AOSCIL=SQRT(AOSCIL)
C
         A_MASS=IN_FIX+IZ_FIX
C
         RKACOU=R0COUL*A_MASS**(1./3.)
C
         VKASOR=1.0
C
C=======================================================================
C
         V0CENT_KAPPAR=V0CENT
         XK_V0C_KAPPAR=XKAPPA
         R0CENT_KAPPAR=0.0000
         XK_R0C_KAPPAR=0.0000
         A0CENT_KAPPAR=0.0000
         XK_A0C_KAPPAR=0.0000
C
         V0SORB_KAPPAR=0.0000
         XK_V0S_KAPPAR=0.0000
         R0SORB_KAPPAR=0.0000
         XK_R0S_KAPPAR=0.0000
         A0SORB_KAPPAR=0.0000
         XK_A0S_KAPPAR=0.0000
C
C=======================================================================
C
         VKACEN_PROTON=V0CENT*(1+XKAPPA*(IN_FIX-IZ_FIX)/A_MASS)
C
         RKACEN_PROTON=R0CENT_PROTON*A_MASS**(1./3.)
C
         RKASOR_PROTON=R0SORB_PROTON*A_MASS**(1./3.)
C
         V0UNIT=-0.25*V0SORB_PROTON*VKACEN_PROTON*(HBAR_C/X_MASS)**2
     *                             /((A_MASS-1)/A_MASS)**2
C
C=======================================================================
C
         AKACEN_PROTON=AKACEN
         VKASOR_PROTON=V0UNIT
         AKASOR_PROTON=AKASOR
C
C=======================================================================
C
         VKACEN_NEUTRS=V0CENT*(1-XKAPPA*(IN_FIX-IZ_FIX)/A_MASS)
C
         RKACEN_NEUTRS=R0CENT_NEUTRS*A_MASS**(1./3.)
C
         RKASOR_NEUTRS=R0SORB_NEUTRS*A_MASS**(1./3.)
C
         V0UNIT=-0.25*V0SORB_NEUTRS*VKACEN_NEUTRS*(HBAR_C/X_MASS)**2
     *                             /((A_MASS-1)/A_MASS)**2
C
C=======================================================================
C
         AKACEN_NEUTRS=AKACEN
         VKASOR_NEUTRS=V0UNIT
         AKASOR_NEUTRS=AKASOR
C
C=======================================================================
C
         DO LQNUMB=NSHELL_PROTON,0,-1
            DO NWSNUM=1,(NSHELL_PROTON-LQNUMB)/2+1
C
               ENERUP_AUXPRO(LQNUMB,NWSNUM)
     *        =
     *         ENERUP_PRONUC(LQNUMB,NWSNUM,INUCLI)
C
               ENERDN_AUXPRO(LQNUMB,NWSNUM)
     *        =
     *         ENERDN_PRONUC(LQNUMB,NWSNUM,INUCLI)
C
               DO N1NUMB=0,(NSHELL_PROTON-LQNUMB)/2
C
                  CMATUP_AUXPRO(LQNUMB,N1NUMB+1,NWSNUM)
     *           =CMATUP_PRONUC(LQNUMB,N1NUMB+1,NWSNUM,INUCLI)
C
                  CMATDN_AUXPRO(LQNUMB,N1NUMB+1,NWSNUM)
     *           =CMATDN_PRONUC(LQNUMB,N1NUMB+1,NWSNUM,INUCLI)
C
               END DO
C
            END DO
         END DO  
C_______________________________________________________________________
C
         DO LQNUMB=NSHELL_NEUTRS,0,-1
            DO NWSNUM=1,(NSHELL_NEUTRS-LQNUMB)/2+1
C
               ENERUP_AUXNEU(LQNUMB,NWSNUM)
     *        =
     *         ENERUP_NEUNUC(LQNUMB,NWSNUM,INUCLI)
C
               ENERDN_AUXNEU(LQNUMB,NWSNUM)
     *        =
     *         ENERDN_NEUNUC(LQNUMB,NWSNUM,INUCLI)
C
               DO N1NUMB=0,(NSHELL_NEUTRS-LQNUMB)/2
C
                  CMATUP_AUXNEU(LQNUMB,N1NUMB+1,NWSNUM)
     *           =CMATUP_NEUNUC(LQNUMB,N1NUMB+1,NWSNUM,INUCLI)
C
                  CMATDN_AUXNEU(LQNUMB,N1NUMB+1,NWSNUM)
     *           =CMATDN_NEUNUC(LQNUMB,N1NUMB+1,NWSNUM,INUCLI)
C
               END DO
C
            END DO
         END DO   
C_______________________________________________________________________
C 
         IFITED=3 ! WS UNIVERSAL single run
C_______________________________________________________________________
C             
         IF (IFDENS.EQ.0) THEN
C
             ISOSPI=1
             IFPROT=1
             IFNEUT=0
C
             ERRMAX_PRINTG=ERRMAX_PRONUC(INUCLI)
C            EABSAV_PRINTG=EABSAV_PRONUC(INUCLI)
             RMSVAL_PRINTG=RMSVAL_PRONUC(INUCLI)
             
             RMSEXP_PRINTG=RMSEXP_PROTON(INUCLI)
             RMSTHE_PRINTG=RMSTHE_PRONUC(INUCLI)
             LEVTHE_PRINTG=LEVTHE_PRONUC(INUCLI)
C
             DO ITHEOR=1,LEVTHE_PRINTG
                ENETHE_PRINTG(ITHEOR)=ENETHE_PRONUC(ITHEOR,INUCLI)
                LABTHE_PRINTG(ITHEOR)=LABTHE_PRONUC(ITHEOR,INUCLI)
             END DO
C
             CALL WRITIN_ENELEV(ISOSPI,IRANDO,LDRAND,IEVALU_PRINTG,
     *                                 PARPOT_PRINTG,CHISQU_MINIML,
     *                                 CHIGRD_PRINTG,ERRMAX_PRINTG,
     *                                 RMSVAL_PRINTG,RMSGLO_PROAUX,
     *                                 RMSEXP_PRINTG,RMSTHE_PRINTG,
     *                                 LEVTHE_PRINTG,ENETHE_PRINTG,
     *                                 LABTHE_PRINTG,INUCLI)
C_______________________________________________________________________
C         
             ISOSPI=0
             IFNEUT=1
             IFPROT=0
C
             ERRMAX_PRINTG=ERRMAX_NEUNUC(INUCLI)
C             EABSAV_PRINTG=EABSAV_NEUNUC(INUCLI)
             RMSVAL_PRINTG=RMSVAL_NEUNUC(INUCLI)
C
             RMSEXP_PRINTG=RMSEXP_NEUTRS(INUCLI)
             RMSTHE_PRINTG=RMSTHE_NEUNUC(INUCLI)
             LEVTHE_PRINTG=LEVTHE_NEUNUC(INUCLI)
C
             DO ITHEOR=1,LEVTHE_PRINTG
                ENETHE_PRINTG(ITHEOR)=ENETHE_NEUNUC(ITHEOR,INUCLI)
                LABTHE_PRINTG(ITHEOR)=LABTHE_NEUNUC(ITHEOR,INUCLI)
             END DO
C
             CALL WRITIN_ENELEV(ISOSPI,IRANDO,LDRAND,IEVALU_PRINTG,
     *                                 PARPOT_PRINTG,CHISQU_MINIML,
     *                                 CHIGRD_PRINTG,ERRMAX_PRINTG,
     *                                 RMSVAL_PRINTG,RMSGLO_NEUAUX,
     *                                 RMSEXP_PRINTG,RMSTHE_PRINTG,
     *                                 LEVTHE_PRINTG,ENETHE_PRINTG,
     *                                 LABTHE_PRINTG,INUCLI)
         END IF
C
      END DO  !INUCLI
C       
C=======================================================================
C
      IF (LOGWRI.GT.0) THEN
          WRITE(LOGFIL,'(/,''Exiting ONERUN_WSUNIV'')')
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
      SUBROUTINE ONERUN
C
      INCLUDE   'MATDIM/NDSPEC.f'
      INCLUDE   'MATDIM/NDPARS.f'
      INCLUDE   'MATDIM/NDNUCL.f'
      INCLUDE   'MATDIM/NDLEXP.f'
      INCLUDE   'MATDIM/NDMESH.f'
      INCLUDE   'MATDIM/ND_RHO.f'
      INCLUDE   'MATDIM/N_NOYX.f'
      INCLUDE   'MATDIM/NDIM_N.f'
      INCLUDE   'MATDIM/NDTITL.f'
C      
      EXTERNAL
     *         DENSIT,DENTHE_INTEGR,DENEXP,DENEXP_INTEGR,
     *         DENTHE_DENEXP
c
      PARAMETER
     *         (NDIM_L=NDIM_N,NDBASE=NDIM_N+1)
C
      CHARACTER
     *          NUCSYM*6,INPSYM*6,TYPCHI*6
      CHARACTER
     *          LABEXP_PROTON*6,LABEXP_NEUTRS*6,
     *          LABTHE_PROTON*6,LABTHE_NEUTRS*6,
     *          LABTHE_PRONUC*6,LABTHE_NEUNUC*6,
     *          LABTHE_PRINTG*6
      CHARACTER
     *          L_CURV_PROTON*6,L_CURV_NEUTRS*6
      CHARACTER
     *          CURNAM*020,L_CURV*100,CURLAB*100,MAIN_T*256,
     *          XTITLE*256,YTITLE*256,TITLES*012,
     *          TITLES_LATEXS*050,NUCNAM_TITLES*010
C
      DIMENSION
     *          PARPOT_PRINTG(1:NDPARS)
      DIMENSION
     *          RMSTHE_PRONUC(1:NDNUCL),
     *          LEVTHE_PRONUC(1:NDNUCL),
     *          ENETHE_PRONUC(1:NDSPEC,1:NDNUCL),
     *          LABTHE_PRONUC(1:NDSPEC,1:NDNUCL)
      DIMENSION
     *          RMSTHE_NEUNUC(1:NDNUCL),
     *          LEVTHE_NEUNUC(1:NDNUCL),
     *          ENETHE_NEUNUC(1:NDSPEC,1:NDNUCL),
     *          LABTHE_NEUNUC(1:NDSPEC,1:NDNUCL)
      DIMENSION
     *          ENETHE_PRINTG(1:NDSPEC),
     *          LABTHE_PRINTG(1:NDSPEC)
      DIMENSION
     *          EABSAV_PRONUC(1:NDNUCL),
     *          EABSAV_NEUNUC(1:NDNUCL),
     *          RMSVAL_PRONUC(1:NDNUCL),
     *          RMSVAL_NEUNUC(1:NDNUCL),
     *          ERRMAX_PRONUC(1:NDNUCL),
     *          ERRMAX_NEUNUC(1:NDNUCL)
      DIMENSION
     *          X_CURV_PROTON(1:NDMESH,1:NDSPEC),
     *          Y_CURV_PROTON(1:NDMESH,1:NDSPEC),
     *          L_CURV_PROTON(1:NDMESH,1:NDSPEC)
      DIMENSION
     *          X_CURV_NEUTRS(1:NDMESH,1:NDSPEC),
     *          Y_CURV_NEUTRS(1:NDMESH,1:NDSPEC),
     *          L_CURV_NEUTRS(1:NDMESH,1:NDSPEC)
      DIMENSION
     *          ND_MAX(1:NDSPEC),
     *          X_CURV(1:NDMESH,1:NDSPEC),
     *          Y_CURV(1:NDMESH,1:NDSPEC),
     *          L_CURV(1:NDMESH,1:NDSPEC)
      DIMENSION 
     *          ENERUP_PRONUC(0:NDIM_L,1:NDBASE,1:NDNUCL),
     *          ENERDN_PRONUC(0:NDIM_L,1:NDBASE,1:NDNUCL),
     *
     *          ENERUP_NEUNUC(0:NDIM_L,1:NDBASE,1:NDNUCL),
     *          ENERDN_NEUNUC(0:NDIM_L,1:NDBASE,1:NDNUCL)
      DIMENSION  
     *          CMATUP_PRONUC(0:NDIM_L,1:NDBASE,1:NDBASE,1:NDNUCL),
     *          CMATDN_PRONUC(0:NDIM_L,1:NDBASE,1:NDBASE,1:NDNUCL),
     *
     *          CMATUP_NEUNUC(0:NDIM_L,1:NDBASE,1:NDBASE,1:NDNUCL),
     *          CMATDN_NEUNUC(0:NDIM_L,1:NDBASE,1:NDBASE,1:NDNUCL)
C          
      COMMON
     *       /GAUSHE/ NGAUSS,NSHELL_PROTON,NSHELL_NEUTRS
      COMMON
     *       /POTPOT/ PARPOT(1:NDPARS)
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
     *       /DATA_V/ AOSCIL_PROTON(1:NDNUCL),
     *                AOSCIL_NEUTRS(1:NDNUCL),
     *                HOMEGA(1:NDNUCL) 
      COMMON
     *       /EXPDEN/ LABELN(1:N_NOYX,1:2),
     *                R_MESH(1:ND_RHO,1:N_NOYX),
     *                RHOEXP(1:ND_RHO,1:N_NOYX)
      COMMON
     *       /DIMSPH/ LDSING_PROTON,LEVTHE_PROTON,
     *                              LEVEXP_PROTON(1:NDNUCL),
     *                LDSING_NEUTRS,LEVTHE_NEUTRS,
     *                              LEVEXP_NEUTRS(1:NDNUCL)
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
     *       /EUDAUX/ ENERUP_AUXPRO(0:NDIM_L,1:NDBASE),
     *                ENERDN_AUXPRO(0:NDIM_L,1:NDBASE),
     *
     *                ENERUP_AUXNEU(0:NDIM_L,1:NDBASE),
     *                ENERDN_AUXNEU(0:NDIM_L,1:NDBASE)
      COMMON
     *       /CUDAUX/ CMATUP_AUXPRO(0:NDIM_L,1:NDBASE,1:NDBASE),
     *                CMATDN_AUXPRO(0:NDIM_L,1:NDBASE,1:NDBASE),
     *
     *                CMATUP_AUXNEU(0:NDIM_L,1:NDBASE,1:NDBASE),
     *                CMATDN_AUXNEU(0:NDIM_L,1:NDBASE,1:NDBASE)
      COMMON
     *       /MODIFR/ RMSEXP_PROTON(1:NDNUCL),
     *                RMSEXP_NEUTRS(1:NDNUCL),
     *                RMSTHE_PROTON,RMSTHE_NEUTRS,
     *                RMODIF_PROTON,RMODIF_NEUTRS,
     *                              PUSHIN,DELTAR    
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
     *       /NUCLEU/ IZ_FIX,IN_FIX,ISOSPI
      COMMON
     *       /ARGUMN/ IARG_R,IARG_Z
      COMMON
     *       /CHI_PN/ CHISQU_PROTON,CHISQU_NEUTRS
      COMMON
     *       /CHICHO/ VALCHI
     *       /DENSTR/ IFDENS
     *       /TENSOR/ IFTENS
      COMMON
     *       /TENOPT/ ICENTT,ISORBT,ITENSR
     *       /TESFIT/ IFTEST,IFFITS,IFMESH,ISPE_P,IFMOCA
     *       /P_N_LM/ IFPROT,IFNEUT,IFBOTH
     *       /MINITE/ IDEFCN,ITECHI
      COMMON
     *       /FITTED/ IFITED
     *       /MASSIV/ IMASIV_PRODUC
      COMMON
     *       /IFPAIR/ IF_PAI,IFGPAI,IFDPAI
      COMMON
     *       /SIMUSE/ AOSCIL,RADCUT,INUCLI
      COMMON
     *       /CUTOFF/ RADIUS_CUTOFF(1:NDNUCL)
      COMMON 
     *       /CHITAK/ IF_SPE,IF_RAD,IF_GAP,IF_FER,IF_DEN,IF_RHO,IF_INV,
     *                LDCHI2,TYPCHI(1:NDNUCL),ITKCHI(1:NDNUCL)
      COMMON
     *       /PARAMT_TITLES/ TITLES(1:NDTITL),
     *                       TITLES_LATEXS(1:NDTITL),
     *                       NUCNAM_LATEXS(1:NDNUCL)
C
      COMMON
     *       /PRFILE/ ISCREN,LOGFIL,LOGCHP,LOGCHN,LOGHIP,LOGHIN,
     *                LOGPRO,LOGNEU,LOGWRI,LSCREN,IRESUL,ICONVE,
     *                NOUTPT,LOGAUX,LOGENE,LOGRAD,LOGGRD,LOGMSH,
     *                                                   LOGBIS 
C       
C=======================================================================
C
C     This subroutine calculates the energy levels, densities... for
C     one single run with the parameter values desired for the user.
C
C     It is only activated if IFTEST=1 for the input file
C
C     The desired parameter values come from the input file, too
C       
C=======================================================================
C
      IF (LOGWRI.GT.0) THEN
          WRITE(LOGFIL,'(/,9X,''Entering ONERUN'',/)')
      END IF
C
C=======================================================================
C
      IFPROT=0
      IFNEUT=0
      IFBOTH=0
C
C=======================================================================
C
C     First, calculating the total weights (protons and neutrons)
C     that will be used as the normalizing constants of the \chi^2
C
      WEIGHT_PROTON=0.0
      WEIGHT_NEUTRS=0.0
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
         END IF
      END DO
C       
C=======================================================================
C
      IMASIV_PRODUC=1
      IRANDO=0
      LDRAND=0
      IEVALU_PRINTG=1
      CHIGRD_PRINTG=999.
      CHISQU_MINIML=999.
C       
C=======================================================================
C
      DO IPARAM=1,NDPARS
         PARPOT_PRINTG(IPARAM)=PARPOT(IPARAM)
      END DO
C
      IF (NUCACT.LT.LDNUCL) RMSGLO_PRINTG=9.9999
C          
      RMSGLO_PROTON=0.
      RMSGLO_NEUTRS=0.
C
      RMSGLO_PROAUX=0.
      RMSGLO_NEUAUX=0.
C
      IDEG_P=0
      IDEG_N=0
C       
C=======================================================================
C
      DO INUCLI=1,LDNUCL
C
         IZ_FIX=NUMB_Z(INUCLI)
         IN_FIX=NUMB_N(INUCLI)
C
         WRITE(0,*)INUCLI,LDNUCL,IZ_FIX,IN_FIX
C
C=======================================================================     
C
C        Diagonalising the Hamiltonian for the single-run mode
C
         IF (LOGWRI.GT.4) THEN
             WRITE(LOGFIL,'(12x,''Entering WS_RUN from ONERUN with: '',
     *                          ''IZ_FIX= '',I3,'' and IN_FIX= '',I3)')
     *                            IZ_FIX,              IN_FIX
         END IF
C
         IDEFCN=0
                                                 I_MODE=0
                                                        I_FLAG=1
         CALL WS_RUN(INUCLI,CHISQU_CHOICE,IABORT,I_MODE,I_FLAG,
     *                             CHISQU_PROTON,CHISQU_NEUTRS)
C
         VALCHI=CHISQU_PROTON/WEIGHT_PROTON
     *         +CHISQU_NEUTRS/WEIGHT_NEUTRS 
C
         IF (ITAKNU(INUCLI).EQ.1) CHISQU_MINIML=VALCHI
C_______________________________________________________________________
C     
         RMSTHE_PRONUC(INUCLI)=RMSTHE_PROTON
         LEVTHE_PRONUC(INUCLI)=LEVTHE_PROTON
         
         DO ITHEOR=1,LEVTHE_PROTON
            ENETHE_PRONUC(ITHEOR,INUCLI)=ENETHE_PROTON(ITHEOR)
            LABTHE_PRONUC(ITHEOR,INUCLI)=LABTHE_PROTON(ITHEOR)
         END DO
C_______________________________________________________________________
C         
         RMSTHE_NEUNUC(INUCLI)=RMSTHE_NEUTRS
         LEVTHE_NEUNUC(INUCLI)=LEVTHE_NEUTRS
         
         DO ITHEOR=1,LEVTHE_PROTON
            ENETHE_NEUNUC(ITHEOR,INUCLI)=ENETHE_NEUTRS(ITHEOR)
            LABTHE_NEUNUC(ITHEOR,INUCLI)=LABTHE_NEUTRS(ITHEOR)
         END DO
C_______________________________________________________________________
C
C        Calculating the TOTAL degenerancy of all the nuclei
C
         DO IEXPER=1,LEVEXP_PROTON(INUCLI)
            IDEG_P=IDEG_P+IDEGEX_PROTON(INUCLI,IEXPER)
         END DO
C
         DO IEXPER=1,LEVEXP_NEUTRS(INUCLI)
            IDEG_N=IDEG_N+IDEGEX_NEUTRS(INUCLI,IEXPER)
         END DO
C_______________________________________________________________________
C                  
         EABSAV_PRONUC(INUCLI)=EABSAV_PROTON
         EABSAV_NEUNUC(INUCLI)=EABSAV_NEUTRS
         
         ERRMAX_PRONUC(INUCLI)=ERRMAX_PROTON
         ERRMAX_NEUNUC(INUCLI)=ERRMAX_NEUTRS
         
         RMSVAL_PRONUC(INUCLI)=SQRT(CHIWEI_PROTON)
         RMSVAL_NEUNUC(INUCLI)=SQRT(CHIWEI_NEUTRS)
         
         RMSGLO_PROTON=RMSGLO_PROTON+CHIWEI_PROTON
         RMSGLO_NEUTRS=RMSGLO_NEUTRS+CHIWEI_NEUTRS

         RMSGLO_PROAUX=RMSGLO_PROAUX+CHIDEG_PROTON
         RMSGLO_NEUAUX=RMSGLO_NEUAUX+CHIDEG_NEUTRS
C_______________________________________________________________________
C
         DO LQNUMB=NSHELL_PROTON,0,-1
            DO NWSNUM=1,(NSHELL_PROTON-LQNUMB)/2+1
C
               ENERUP_PRONUC(LQNUMB,NWSNUM,INUCLI)
     *        =
     *         ENERUP_AUXPRO(LQNUMB,NWSNUM)
C
               ENERDN_PRONUC(LQNUMB,NWSNUM,INUCLI)
     *        =
     *         ENERDN_AUXPRO(LQNUMB,NWSNUM)
C
               DO N1NUMB=0,(NSHELL_PROTON-LQNUMB)/2
C
                  CMATUP_PRONUC(LQNUMB,N1NUMB+1,NWSNUM,INUCLI)
     *           =CMATUP_AUXPRO(LQNUMB,N1NUMB+1,NWSNUM)
C
                  CMATDN_PRONUC(LQNUMB,N1NUMB+1,NWSNUM,INUCLI)
     *           =CMATDN_AUXPRO(LQNUMB,N1NUMB+1,NWSNUM)
C
               END DO
C
            END DO
         END DO  
C_______________________________________________________________________
C
         DO LQNUMB=NSHELL_NEUTRS,0,-1
            DO NWSNUM=1,(NSHELL_NEUTRS-LQNUMB)/2+1
C
               ENERUP_NEUNUC(LQNUMB,NWSNUM,INUCLI)
     *        =
     *         ENERUP_AUXNEU(LQNUMB,NWSNUM)
C
               ENERDN_NEUNUC(LQNUMB,NWSNUM,INUCLI)
     *        =
     *         ENERDN_AUXNEU(LQNUMB,NWSNUM)
C
               DO N1NUMB=0,(NSHELL_NEUTRS-LQNUMB)/2
C
                  CMATUP_NEUNUC(LQNUMB,N1NUMB+1,NWSNUM,INUCLI)
     *           =CMATUP_AUXNEU(LQNUMB,N1NUMB+1,NWSNUM)
C
                  CMATDN_NEUNUC(LQNUMB,N1NUMB+1,NWSNUM,INUCLI)
     *           =CMATDN_AUXNEU(LQNUMB,N1NUMB+1,NWSNUM)
C
               END DO
C
            END DO
         END DO   
C_______________________________________________________________________
C          
      END DO !INUCLI
C_______________________________________________________________________
C      
      RMSGLO_PROTON=SQRT(RMSGLO_PROTON)
      RMSGLO_NEUTRS=SQRT(RMSGLO_NEUTRS)
C
      RMSGLO_PROAUX=SQRT(RMSGLO_PROAUX/IDEG_P)
      RMSGLO_NEUAUX=SQRT(RMSGLO_NEUAUX/IDEG_N)
C
C=======================================================================
C
      DO INUCLI=1,LDNUCL
C
         IZ_FIX=NUMB_Z(INUCLI)
         IN_FIX=NUMB_N(INUCLI)
C
         CALL INTROD(INUCLI)
C_______________________________________________________________________
C
         DO LQNUMB=NSHELL_PROTON,0,-1
            DO NWSNUM=1,(NSHELL_PROTON-LQNUMB)/2+1
C
               ENERUP_AUXPRO(LQNUMB,NWSNUM)
     *        =
     *         ENERUP_PRONUC(LQNUMB,NWSNUM,INUCLI)
C
               ENERDN_AUXPRO(LQNUMB,NWSNUM)
     *        =
     *         ENERDN_PRONUC(LQNUMB,NWSNUM,INUCLI)
C
               DO N1NUMB=0,(NSHELL_PROTON-LQNUMB)/2
C
                  CMATUP_AUXPRO(LQNUMB,N1NUMB+1,NWSNUM)
     *           =CMATUP_PRONUC(LQNUMB,N1NUMB+1,NWSNUM,INUCLI)
C
                  CMATDN_AUXPRO(LQNUMB,N1NUMB+1,NWSNUM)
     *           =CMATDN_PRONUC(LQNUMB,N1NUMB+1,NWSNUM,INUCLI)
C
               END DO
C
            END DO
         END DO  
C_______________________________________________________________________
C
         DO LQNUMB=NSHELL_NEUTRS,0,-1
            DO NWSNUM=1,(NSHELL_NEUTRS-LQNUMB)/2+1
C
               ENERUP_AUXNEU(LQNUMB,NWSNUM)
     *        =
     *         ENERUP_NEUNUC(LQNUMB,NWSNUM,INUCLI)
C
               ENERDN_AUXNEU(LQNUMB,NWSNUM)
     *        =
     *         ENERDN_NEUNUC(LQNUMB,NWSNUM,INUCLI)
C
               DO N1NUMB=0,(NSHELL_NEUTRS-LQNUMB)/2
C
                  CMATUP_AUXNEU(LQNUMB,N1NUMB+1,NWSNUM)
     *           =CMATUP_NEUNUC(LQNUMB,N1NUMB+1,NWSNUM,INUCLI)
C
                  CMATDN_AUXNEU(LQNUMB,N1NUMB+1,NWSNUM)
     *           =CMATDN_NEUNUC(LQNUMB,N1NUMB+1,NWSNUM,INUCLI)
C
               END DO
C
            END DO
         END DO   
C_______________________________________________________________________
C 
         IFITED=2 ! single run
C_______________________________________________________________________
C             
         IF (IFDENS.EQ.0) THEN
C
             ISOSPI=1
             IFPROT=1
             IFNEUT=0
C
             ERRMAX_PRINTG=ERRMAX_PRONUC(INUCLI)
C            EABSAV_PRINTG=EABSAV_PRONUC(INUCLI)
             RMSVAL_PRINTG=RMSVAL_PRONUC(INUCLI)
             
             RMSEXP_PRINTG=RMSEXP_PROTON(INUCLI)
             RMSTHE_PRINTG=RMSTHE_PRONUC(INUCLI)
             LEVTHE_PRINTG=LEVTHE_PRONUC(INUCLI)
C
             DO ITHEOR=1,LEVTHE_PRINTG
                ENETHE_PRINTG(ITHEOR)=ENETHE_PRONUC(ITHEOR,INUCLI)
                LABTHE_PRINTG(ITHEOR)=LABTHE_PRONUC(ITHEOR,INUCLI)
             END DO
C
             CALL WRITIN_ENELEV(ISOSPI,IRANDO,LDRAND,IEVALU_PRINTG,
     *                                 PARPOT_PRINTG,CHISQU_MINIML,
     *                                 CHIGRD_PRINTG,ERRMAX_PRINTG,
     *                                 RMSVAL_PRINTG,RMSGLO_PROAUX,
     *                                 RMSEXP_PRINTG,RMSTHE_PRINTG,
     *                                 LEVTHE_PRINTG,ENETHE_PRINTG,
     *                                 LABTHE_PRINTG,INUCLI)
C_______________________________________________________________________
C         
             ISOSPI=0
             IFNEUT=1
             IFPROT=0
C
             ERRMAX_PRINTG=ERRMAX_NEUNUC(INUCLI)
C             EABSAV_PRINTG=EABSAV_NEUNUC(INUCLI)
             RMSVAL_PRINTG=RMSVAL_NEUNUC(INUCLI)
C
             RMSEXP_PRINTG=RMSEXP_NEUTRS(INUCLI)
             RMSTHE_PRINTG=RMSTHE_NEUNUC(INUCLI)
             LEVTHE_PRINTG=LEVTHE_NEUNUC(INUCLI)
C
             DO ITHEOR=1,LEVTHE_PRINTG
                ENETHE_PRINTG(ITHEOR)=ENETHE_NEUNUC(ITHEOR,INUCLI)
                LABTHE_PRINTG(ITHEOR)=LABTHE_NEUNUC(ITHEOR,INUCLI)
             END DO
C
             CALL WRITIN_ENELEV(ISOSPI,IRANDO,LDRAND,IEVALU_PRINTG,
     *                                 PARPOT_PRINTG,CHISQU_MINIML,
     *                                 CHIGRD_PRINTG,ERRMAX_PRINTG,
     *                                 RMSVAL_PRINTG,RMSGLO_NEUAUX,
     *                                 RMSEXP_PRINTG,RMSTHE_PRINTG,
     *                                 LEVTHE_PRINTG,ENETHE_PRINTG,
     *                                 LABTHE_PRINTG,INUCLI)
         END IF
C_______________________________________________________________________
C         
         IF (IFDENS.EQ.1) THEN
C
             ISOSPI=1
C
             ERRMAX_PRINTG=ERRMAX_PRONUC(INUCLI)
C             EABSAV_PRINTG=EABSAV_PRONUC(INUCLI)
             RMSVAL_PRINTG=RMSVAL_PRONUC(INUCLI)
C
             RMSEXP_PRINTG=RMSEXP_PROTON(INUCLI)
             RMSTHE_PRINTG=RMSTHE_PRONUC(INUCLI)
             LEVTHE_PRINTG=LEVTHE_PRONUC(INUCLI)
C
             DO ITHEOR=1,LEVTHE_PRINTG
                ENETHE_PRINTG(ITHEOR)=ENETHE_PRONUC(ITHEOR,INUCLI)
                LABTHE_PRINTG(ITHEOR)=LABTHE_PRONUC(ITHEOR,INUCLI)
             END DO
C
             CALL WRITIN_ENELEV(ISOSPI,IRANDO,LDRAND,IEVALU_PRINTG,
     *                                 PARPOT_PRINTG,CHISQU_MINIML,
     *                                 CHIGRD_PRINTG,ERRMAX_PRINTG,
     *                                 RMSVAL_PRINTG,RMSGLO_PROAUX,
     *                                 RMSEXP_PRINTG,RMSTHE_PRINTG,
     *                                 LEVTHE_PRINTG,ENETHE_PRINTG,
     *                                 LABTHE_PRINTG,INUCLI)
C_______________________________________________________________________
C     
             ISOSPI=0
C
             ERRMAX_PRINTG=ERRMAX_NEUNUC(INUCLI)
C             EABSAV_PRINTG=EABSAV_NEUNUC(INUCLI)
             RMSVAL_PRINTG=RMSVAL_NEUNUC(INUCLI)
C
             RMSEXP_PRINTG=RMSEXP_NEUTRS(INUCLI)
             RMSTHE_PRINTG=RMSTHE_NEUNUC(INUCLI)
             LEVTHE_PRINTG=LEVTHE_NEUNUC(INUCLI)
C
             DO ITHEOR=1,LEVTHE_PRINTG
                ENETHE_PRINTG(ITHEOR)=ENETHE_NEUNUC(ITHEOR,INUCLI)
                LABTHE_PRINTG(ITHEOR)=LABTHE_NEUNUC(ITHEOR,INUCLI)
             END DO
C
             CALL WRITIN_ENELEV(ISOSPI,IRANDO,LDRAND,IEVALU_PRINTG,
     *                                 PARPOT_PRINTG,CHISQU_MINIML,
     *                                 CHIGRD_PRINTG,ERRMAX_PRINTG,
     *                                 RMSVAL_PRINTG,RMSGLO_NEUAUX,
     *                                 RMSEXP_PRINTG,RMSTHE_PRINTG,
     *                                 LEVTHE_PRINTG,ENETHE_PRINTG,
     *                                 LABTHE_PRINTG,INUCLI)
         END IF
C
C=======================================================================
C
C        Now, preparing outputs for density-dependent functions
C
C=======================================================================
C
         ISPLIT=0
         CURNAM='Dens-'
C_______________________________________________________________________
C
C        Total Proton Densities: Theory vs. Experiment
C
         ISOSPI=1
C
         N_CURV=2
C                        
         DO I_CURV=1,N_CURV
            ND_MAX(I_CURV)=60
         END DO
C
         IARG_R=1
         IARG_Z=0
C                        
         AOSCIL=AOSCIL_PROTON(INUCLI)
         RADCUT=RADIUS_CUTOFF(INUCLI)
C
         IF (IF_PAI.EQ.0) CURLAB='$\\rho^{th}_{\\pi}(r)$'
         IF (IF_PAI.EQ.1) CURLAB='$\\rho^{pair}_{\\pi}(r)$'
C                        
         DO IPOINT=1,ND_MAX(N_CURV)
C           
            RARGUM=0.0001+(IPOINT-1)*0.2 
C           
            X_CURV(IPOINT,1)=RARGUM
            Y_CURV(IPOINT,1)=DENSIT(RARGUM)
            L_CURV(IPOINT,1)=CURLAB
C
            X_CURV(IPOINT,2)=RARGUM
            Y_CURV(IPOINT,2)=DENEXP(RARGUM)
            L_CURV(IPOINT,2)='$\\rho^{exp}_{\\pi}(r)$'
C                        
         END DO
C                     
         MAIN_T='Proton Density Profiles'
         XTITLE='r (fm)'
         YTITLE='Density (N/fm$^3$)'
         ISOSPI=1 
C                     
         CALL DENSIT_OUTPUT(CURNAM,N_CURV,ND_MAX,X_CURV,Y_CURV,L_CURV,
     *                      ISOSPI,IZ_FIX,IN_FIX,MAIN_T,XTITLE,YTITLE,
     *                                           ISPLIT,INUCLI,LDRAND) 
C_______________________________________________________________________
C
C        Total Neutron Density: Theory
C
         ISOSPI=0
C
         N_CURV=1
         ND_MAX(I_CURV)=60
C
         IARG_R=1
         IARG_Z=0
C
         AOSCIL=AOSCIL_NEUTRS(INUCLI)
         IF (IF_PAI.EQ.0) CURLAB='$\\rho^{th}_{\\nu}(r)$'
         IF (IF_PAI.EQ.1) CURLAB='$\\rho^{pair}_{\\nu}(r)$'
C                        
         DO IPOINT=1,ND_MAX(I_CURV)
C           
            RARGUM=0.0001+(IPOINT-1)*0.2 
C           
            X_CURV(IPOINT,1)=RARGUM
            Y_CURV(IPOINT,1)=DENSIT(RARGUM)
            L_CURV(IPOINT,1)=CURLAB
C                        
         END DO
C                     
         MAIN_T='Neutron Density Profile'
         XTITLE='r (fm)'
         YTITLE='Density (N/fm$^3$)'
         ISOSPI=0 
C                     
         CALL DENSIT_OUTPUT(CURNAM,N_CURV,ND_MAX,X_CURV,Y_CURV,L_CURV,
     *                      ISOSPI,IZ_FIX,IN_FIX,MAIN_T,XTITLE,YTITLE,
     *                                           ISPLIT,INUCLI,LDRAND)
C
C=======================================================================
C
C        Level densities summed one after the other
C        First protons, then neutrons
C_______________________________________________________________________
C
         ISPLIT=1                
         CURNAM='RhoS-'
         ND_MAX_AUXIL=60
C_______________________________________________________________________
C                    
         ISOSPI=1
C                     
         DO IPOINT=1,ND_MAX_AUXIL
C                        
            XPOINT=0.0001+(IPOINT-1)*0.2
C                        
            ZARGUM=(XPOINT/AOSCIL_PROTON(INUCLI))**2
C
            CALL RHONLJ_SUMMED(INUCLI,ISOSPI,ZARGUM,
     *                         IPOINT,N_CURV,Y_CURV,L_CURV)
C     
            DO I_CURV=1,N_CURV
C                           
               X_CURV(IPOINT,I_CURV)=XPOINT
               ND_MAX(I_CURV)=ND_MAX_AUXIL
C                           
            END DO
C     
         END DO
C                     
         WRITE(MAIN_T,'(''Calculated Proton Density Profiles'')')
         XTITLE='r (fm)'
         YTITLE='Density (N/fm$^3$)'
C                     
         CALL DENSIT_OUTPUT(CURNAM,N_CURV,ND_MAX,X_CURV,Y_CURV,L_CURV,
     *                      ISOSPI,IZ_FIX,IN_FIX,MAIN_T,XTITLE,YTITLE,
     *                                           ISPLIT,INUCLI,LDRAND)
C_______________________________________________________________________
C                     
         ISOSPI=0
C                     
         DO IPOINT=1,ND_MAX_AUXIL
C                        
            XPOINT=0.0001+(IPOINT-1)*0.2
C                        
            ZARGUM=(XPOINT/AOSCIL_NEUTRS(INUCLI))**2
C
            CALL RHONLJ_SUMMED(INUCLI,ISOSPI,ZARGUM,
     *                         IPOINT,N_CURV,Y_CURV,L_CURV)
C     
            DO I_CURV=1,N_CURV
C                           
               X_CURV(IPOINT,I_CURV)=XPOINT
               ND_MAX(I_CURV)=ND_MAX_AUXIL
C                           
            END DO
C     
         END DO
                     
         WRITE(MAIN_T,'(''Calculated Neutron Density Profiles'')')
         XTITLE='r (fm)'
         YTITLE='Density (N/fm$^3$)'
C                     
         CALL DENSIT_OUTPUT(CURNAM,N_CURV,ND_MAX,X_CURV,Y_CURV,L_CURV,
     *                      ISOSPI,IZ_FIX,IN_FIX,MAIN_T,XTITLE,YTITLE,
     *                                           ISPLIT,INUCLI,LDRAND)
C
C=======================================================================
C        
      END DO  !INUCLI
C       
C=======================================================================
C
      IF (LOGWRI.GT.0) THEN
          WRITE(LOGFIL,'(/,''Exiting ONERUN'')')
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
      SUBROUTINE SPENER_OFPARA(LDRAND)
C
      INCLUDE 'MATDIM/NDMESH.f'
      INCLUDE 'MATDIM/NDSPEC.f'
      INCLUDE 'MATDIM/NDPARS.f'
      INCLUDE 'MATDIM/NDNUCL.f'
      INCLUDE 'MATDIM/NDLEXP.f'
      INCLUDE 'MATDIM/NDIM_M.f'
      INCLUDE 'MATDIM/NDTITL.f'
      INCLUDE 'MATDIM/NDLAST.f'
C
      PARAMETER
     *          (NDMES2=NDMESH*NDMESH)
C
      CHARACTER
     *          NUCSYM*6,INPSYM*6,TYPCHI*6,VERSIO*3
      CHARACTER
     *          FILNAM*256,L_CURV*100,NUCNAM*05,TITLES*12,
     *          ISONAM*002,TEXLAM*020,TITPAR*13,FITNUC*02,
     *          TITLES_LATEXS*050,NUCNAM_LATEXS*010
      CHARACTER
     *          LABTHE_PROTON*6,LABTHE_NEUTRS*6,
     *          LABEXP_PROTON*6,LABEXP_NEUTRS*6
      CHARACTER
     *          LABREF*6,LABTHE*6,LABEXP*6,LABTEX*11,L_CURV_AUXILI*100,
     *          L_CURV_AUXPRO*100,L_CURV_AUXNEU*100,
     *          LABREF_PROTON*6,LABREF_NEUTRS*6
      CHARACTER
     *          MAIN_T*256,XTITLE*256,YTITLE*256,PARNAM*8,TITMSH*8
C
      DIMENSION
     *          STPMSH(1:NDIM_M),
     *          XINITS(1:NDPARS),
     *          ARGMNT(1:NDPARS)
      DIMENSION
     *          TITMSH(1:NDPARS),
     *          TITPAR(1:NDPARS),
     *          PARPOT_PRINTI(1:NDPARS,1:NDMES2)
      DIMENSION
     *          NDPOIN(1:NDSPEC),
     *          X_CURV(1:NDMESH,1:NDSPEC),
     *          Y_CURV(1:NDMESH,1:NDSPEC),
     *          L_CURV(1:NDMESH,1:NDSPEC),
     *          ICOLOR(1:NDSPEC),
     *          ITHICK(1:NDSPEC)
      DIMENSION
     *          X_CURV_AUXILI(1:NDMESH,1:NDSPEC,1:NDNUCL),
     *          Y_CURV_AUXILI(1:NDMESH,1:NDSPEC,1:NDNUCL),
     *          L_CURV_AUXILI(1:NDMESH,1:NDSPEC,1:NDNUCL),
     *          ICOLOR_AUXILI(1:NDSPEC,1:NDNUCL),
     *          ITHICK_AUXILI(1:NDSPEC,1:NDNUCL)
      DIMENSION
     *          NDPOIN_PROTON(1:NDSPEC),
     *          X_CURV_AUXPRO(1:NDMESH,1:NDSPEC,1:NDNUCL),
     *          Y_CURV_AUXPRO(1:NDMESH,1:NDSPEC,1:NDNUCL),
     *          L_CURV_AUXPRO(1:NDMESH,1:NDSPEC,1:NDNUCL),
     *          ICOLOR_AUXPRO(1:NDSPEC,1:NDNUCL),
     *          ITHICK_AUXPRO(1:NDSPEC,1:NDNUCL)
      DIMENSION
     *          NDPOIN_NEUTRS(1:NDSPEC),
     *          X_CURV_AUXNEU(1:NDMESH,1:NDSPEC,1:NDNUCL),
     *          Y_CURV_AUXNEU(1:NDMESH,1:NDSPEC,1:NDNUCL),
     *          L_CURV_AUXNEU(1:NDMESH,1:NDSPEC,1:NDNUCL),
     *          ICOLOR_AUXNEU(1:NDSPEC,1:NDNUCL),
     *          ITHICK_AUXNEU(1:NDSPEC,1:NDNUCL)
      DIMENSION
     *          RMSVAL(1:NDSPEC,1:NDNUCL),
     *          RMSVAL_PROTON(1:NDSPEC,1:NDNUCL),
     *          RMSVAL_NEUTRS(1:NDSPEC,1:NDNUCL)
      DIMENSION
     *          LABREF(1:NDSPEC),
     *          LABREF_PROTON(1:NDSPEC),
     *          LABREF_NEUTRS(1:NDSPEC)
      COMMON
     *       /DENSTR/ IFDENS
      COMMON
     *       /TENSOR/ IFTENS
      COMMON
     *       /TENOPT/ ICENTT,ISORBT,ITENSR 
      COMMON
     *       /DEPPAI/ IFDEEP,IFPRON
      COMMON
     *       /IFPAIR/ IF_PAI,IFGPAI,IFDPAI
      COMMON 
     *       /CHITAK/ IF_SPE,IF_RAD,IF_GAP,IF_FER,IF_DEN,IF_RHO,IF_INV,
     *                LDCHI2,TYPCHI(1:NDNUCL),ITKCHI(1:NDNUCL)
      COMMON
     *       /IFLAMB/ IFPAR1,IFPAR2,IFPAR3,IFPAR4,
     *                IFPAR5,IFPAR6,IFPAR7,IFPAR8,
     *                IFPA09,IFPA10,IFPA11,IFPA12
      COMMON
     *       /PARCOR/ IFCORR,IFRCVC,IFRCAC,IFVCAC,IFRSVS,IFRSAS,IFVSAS
      COMMON
     *       /POTPOT/ PARPOT(1:NDPARS)
      COMMON
     *       /MESHIN/ I_MESH(1:NDPARS),
     *                XMIN_I(1:NDPARS),
     *                XMAX_I(1:NDPARS),
     *                MAXPAR(1:NDPARS)
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
     *       /VLIMIT/ VMISTR(1:NDPARS),
     *                VMIMIN(1:NDPARS),
     *                VMIMAX(1:NDPARS),
     *                VMISTP(1:NDPARS)
      COMMON
     *       /SEARCH/ IFTAKE(1:NDPARS)
      COMMON
     *       /SPECTR/ ENETHE_PROTON(1:NDSPEC),
     *                ENETHE_NEUTRS(1:NDSPEC)
      COMMON
     *       /THELAB/ LABTHE_PROTON(1:NDSPEC),
     *                LABTHE_NEUTRS(1:NDSPEC)  
      COMMON
     *       /QUANTN/ NWSSPH_PROTON(1:NDSPEC),
     *                LWSSPH_PROTON(1:NDSPEC),
     *                JWSSPH_PROTON(1:NDSPEC),
     *
     *                NWSSPH_NEUTRS(1:NDSPEC),
     *                LWSSPH_NEUTRS(1:NDSPEC),
     *                JWSSPH_NEUTRS(1:NDSPEC)
      COMMON
     *       /EX_ENE/ EXPEXP_PROTON(1:NDNUCL,1:NDLEXP),
     *                IDEGEX_PROTON(1:NDNUCL,1:NDLEXP),
     *
     *                EXPEXP_NEUTRS(1:NDNUCL,1:NDLEXP),
     *                IDEGEX_NEUTRS(1:NDNUCL,1:NDLEXP)
      COMMON
     *       /EX_LAB/ LABEXP_PROTON(1:NDNUCL,1:NDLEXP),
     *                LABEXP_NEUTRS(1:NDNUCL,1:NDLEXP)
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
     *       /AUXPOT/ VKACEN_PROTON,RKACEN_PROTON,AKACEN_PROTON,
     *                VKASOR_PROTON,RKASOR_PROTON,AKASOR_PROTON,
     *                VKAEFM_PROTON,RKAEFM_PROTON,AKAEFM_PROTON,
     *                                            V0UNIT_PROTON,
     *                VKACEN_NEUTRS,RKACEN_NEUTRS,AKACEN_NEUTRS,
     *                VKASOR_NEUTRS,RKASOR_NEUTRS,AKASOR_NEUTRS,
     *                VKAEFM_NEUTRS,RKAEFM_NEUTRS,AKAEFM_NEUTRS,
     *                                            V0UNIT_NEUTRS,
     *                                                   RKACOU,
     *                              ALAMPP,ALAMPN,ALAMNP,ALAMNN,
     *                              TLAMPP,TLAMPN,TLAMNP,TLAMNN,
     *                              CLAMPP,CLAMPN,CLAMNP,CLAMNN
      COMMON
     *       /PARAMT_TITLES/ TITLES(1:NDTITL),
     *                       TITLES_LATEXS(1:NDTITL),
     *                       NUCNAM_LATEXS(1:NDNUCL)
      COMMON
     *       /VERSIN/ VERSIO
      COMMON
     *       /PRFILE/ ISCREN,LOGFIL,LOGCHP,LOGCHN,LOGHIP,LOGHIN,
     *                LOGPRO,LOGNEU,LOGWRI,LSCREN,IRESUL,ICONVE,
     *                NOUTPT,LOGAUX,LOGENE,LOGRAD,LOGGRD,LOGMSH,
     *                                                   LOGBIS 
C       
C=======================================================================
C
C     This subroutine prepares all the input parapeters for WS_RUN
C     in order to be able to plot the single particle energies  as 
C     a function of O N E parameter, that the user has choosen via 
C     NAMELI. 
C
C     ATTENTION: this subroutine uses some vectors that are also used
C     in MESHIN_MAPING. This has been deliberately programmed.
C
C=======================================================================
C
C     IFTAKE .EQ. -1: parameter eliminated via parametric correlation
C     IFTAKE .EQ. +1: parameter entering to the minimization
C     IFTAKE .EQ. +4: parameter considered constant but printed anyway
C
C=======================================================================     
C
      WRITE(LOGFIL,'(/,''Entering to SPENER_OFPARA'',/)')
C
C=======================================================================     
C     Over which parameter the levels are going to be plotted
C=======================================================================
C
      ICOUNT=0
C
      DO IPARAM=1,NDPARS
         IF (I_MESH(IPARAM).EQ.1) THEN
C
             ICOUNT=ICOUNT+1
C
             STPMSH(IPARAM)=(XMAX_I(IPARAM)-XMIN_I(IPARAM))
     *                     /(MAXPAR(IPARAM)-1)
C
             MAXMSH=MAXPAR(IPARAM)
C
             INDEXI=IPARAM
C
             PARNAM=TITLES(IPARAM)(3:10)
             TITMSH(ICOUNT)=PARNAM
C
         END IF
      END DO
C
      IF (ICOUNT.NE.1) THEN
C
          WRITE(LSCREN,'(''ALARM from SPENER_OFPARA: '',
     *                   ''ICOUNT='',I2,''when it has to be '',
     *                   ''equal to 1!! - Change the input file '')')
     *                    ICOUNT
          STOP 'Stop - Change the input file in ISPE_P option'
C
      END IF
C
      ICMESH=ICOUNT
C
C=======================================================================
C         Over which parameters we are going to minimise --> ARGMNT(I)
C=======================================================================
C        
      ICOUNT_OFPARS=0
C
      IFPROT=0
      IFNEUT=0
      IFBOTH=0
C
      DO I=1,NDPARS
         IF (IFTAKE(I).EQ.1) THEN
C
             ICOUNT_OFPARS=ICOUNT_OFPARS+1
             ARGMNT(ICOUNT_OFPARS)=VMISTR(I)
             TITPAR(ICOUNT_OFPARS)=TITLES(I)(3:10)
C
             IF (I.LE.20 .AND. IFDENS.EQ.0) THEN
                 IFPROT=1
             END IF
             IF (I.GT.20 .AND. IFDENS.EQ.0) THEN
                 IFNEUT=1
             END IF
             IF (I.GE.51 .AND. IFDENS.EQ.0) THEN
                 IFPROT=0
                 IFNEUT=0
                 IFBOTH=1
             END IF
C
         END IF
      END DO
C
      IF (LOGWRI.GT.0) THEN
          IF (IFPROT.EQ.1) 
     *        WRITE(LOGFIL,'(9X,''SPENER_OFPARA for protons'',/)')
          IF (IFNEUT.EQ.1) 
     *        WRITE(LOGFIL,'(9X,''SPENER_OFPARA for neutrons'',/)')
      END IF
C
C=======================================================================
C
      IPRINT_PROTON=0
      IPRINT_NEUTRS=0
C
      IF ((IFDENS.EQ.0 .AND. IFNEUT.EQ.1) .OR.
     *    (IFDENS.EQ.0 .AND. IFBOTH.EQ.1) .OR.
     *     IFDENS.EQ.1) THEN
           IPRINT_NEUTRS=1
      END IF
C
      IF ((IFDENS.EQ.0 .AND. IFPROT.EQ.1) .OR.
     *    (IFDENS.EQ.0 .AND. IFBOTH.EQ.1) .OR.
     *     IFDENS.EQ.1) THEN
           IPRINT_PROTON=1
      END IF
C
C=======================================================================
C
      IF (ICOUNT_OFPARS.NE.0) GO TO 1 ! ==>> Minimization required
C
C=======================================================================
C
      IF (ICOUNT_OFPARS.EQ.0) THEN
C
          IFPROT=0
          IFNEUT=0
          IFBOTH=0
C
          DO I=1,NDPARS
             IF (I_MESH(I).EQ.1) THEN
C
                 IF (I.LE.20 .AND. IFDENS.EQ.0) THEN
                     IFPROT=1
                 END IF
                 IF (I.GT.20 .AND. IFDENS.EQ.0) THEN
                     IFNEUT=1
                 END IF
                 IF (I.GE.51 .AND. IFDENS.EQ.0) THEN
                     IFPROT=0
                     IFNEUT=0
                     IFBOTH=1
                 END IF
C
             END IF
          END DO
C
      END IF
C
      IF (LOGWRI.GT.0) THEN
          IF (IFPROT.EQ.1) 
     *        WRITE(LOGFIL,'(9X,''SPENER_OFPARA for protons'',/)')
          IF (IFNEUT.EQ.1) 
     *        WRITE(LOGFIL,'(9X,''SPENER_OFPARA for neutrons'',/)')
      END IF
C
C=======================================================================
C
      IPRINT_PROTON=0
      IPRINT_NEUTRS=0
C
      IF ((IFDENS.EQ.0 .AND. IFNEUT.EQ.1) .OR.
     *    (IFDENS.EQ.0 .AND. IFBOTH.EQ.1) .OR.
     *     IFDENS.EQ.1) THEN
           IPRINT_NEUTRS=1
      END IF
C
      IF ((IFDENS.EQ.0 .AND. IFPROT.EQ.1) .OR.
     *    (IFDENS.EQ.0 .AND. IFBOTH.EQ.1) .OR.
     *     IFDENS.EQ.1) THEN
           IPRINT_PROTON=1
      END IF
C
C=======================================================================
C
      IF (IPRINT_PROTON.EQ.1) THEN
          CONTINUE
      ELSE
          GO TO 3
      END IF
C
C=======================================================================
C
C     NO MINIMIZATION - PROTONS
C
      IF (LOGWRI.GT.0) 
     *    WRITE(LOGFIL,'(9x,''NO Minimization option'',/)')
C
      JPARAM=INDEXI
C
      DO INUCLI=1,LDNUCL
         IF (ITAKNU(INUCLI).EQ.1) THEN
C
             DO INDEX1=1,MAXMSH
C
                IF (LOGWRI.GT.0) THEN
                    WRITE(NOUTPT,'(/,80(''#''),/,''#'',T80,''#'',/,
     *                             ''#  INDEX1= '',I6,T80,''#'',/,
     *                             ''#'',T80,''#'',/,80(''#''))')INDEX1
                END IF
C
                XINITS(JPARAM)=XMIN_I(JPARAM)+(INDEX1-1)*STPMSH(JPARAM)
                PARPOT(JPARAM)=XINITS(JPARAM)
C
                IDEFCN=0
                I_MODE=0
                I_FLAG=1
C
                CALL WS_RUN(INUCLI,CHISQU_CHOICE,IABORT,I_MODE,I_FLAG,
     *                                    CHISQU_PROTON,CHISQU_NEUTRS)
C
C               Updating the parameter vector for printing purposes
C
                DO KPARAM=1,NDPARS
                   PARPOT_PRINTI(KPARAM,INDEX1)=PARPOT(KPARAM)
                END DO
C
C               storing the rms value for protons
C
                RMSVAL(INDEX1,INUCLI)=SQRT(CHIWEI_PROTON)
C
C               Storing the values of the levels for each nucleus 
C               for the plotting output (depending on the reference,
C                                                         see below)
C
                NDCURV=LEVTHE_PROTON
C
C               Storing the reference order of the levels
C
                IF (INDEX1.EQ.1) THEN
                    DO IREFER=1,NDCURV
                       LABREF(IREFER)=LABTHE_PROTON(IREFER)
                    END DO
                END IF
C
                DO ICURVE=1,NDCURV
C
                   X_CURV(INDEX1,ICURVE)=XINITS(JPARAM)
C
                   DO ITHEOR=1,NDCURV
                      IF (LABTHE_PROTON(ITHEOR).EQ.LABREF(ICURVE))THEN
C
                          Y_CURV(INDEX1,ICURVE)
     *                   =ENETHE_PROTON(ITHEOR)
C
                          LABTHE=LABTHE_PROTON(ITHEOR)
                          CALL LATEXS_LABELS(LABTHE,LABTEX)
                          L_CURV(INDEX1,ICURVE)=LABTEX
C
                          LQNUMB=LWSSPH_PROTON(ITHEOR)
                          LPARIT=(-1)**LQNUMB
C
                          IF (LPARIT.EQ.+1) ICOLOR(ICURVE)=44
C
                          IF (LPARIT.EQ.-1) ICOLOR(ICURVE)=53
C
                          ITHICK(ICURVE)=4
C
                      END IF
                   END DO
C                   
                END DO !ICURVE
C
                DO I=1,NDCURV
                   NDPOIN(I)=MAXMSH
                END DO
C
             END DO !INDEX1
C
C            Looking the correspondance with the experiment
C
             IEXTRA=0
             DO IEXPER=1,LEVEXP_PROTON(INUCLI)
                DO ITHEOR=1,NDCURV
C
                   LABEXP=LABEXP_PROTON(INUCLI,IEXPER)
                   CALL LATEXS_LABELS(LABEXP,LABTEX)
C
                   IF (LABTEX.EQ.L_CURV(1,ITHEOR)) THEN
C
                       DO IPOINT=1,NDPOIN(ITHEOR)-1
C
                          IF ((EXPEXP_PROTON(INUCLI,IEXPER)
     *                   .GE.  Y_CURV(IPOINT,ITHEOR))
     *                   .AND.(EXPEXP_PROTON(INUCLI,IEXPER)
     *                   .LT.  Y_CURV(IPOINT+1,ITHEOR)))
     *                    THEN
C
                             IEXTRA=IEXTRA+1
C
                             IF (IPOINT.EQ.1) THEN
                                 X_CURV(1,NDCURV+IEXTRA)
     *                          =X_CURV(IPOINT,ITHEOR)
                             ELSE
                                 X_CURV(1,NDCURV+IEXTRA)
     *                          =X_CURV(IPOINT-1,ITHEOR)
                             END IF
C
                             X_CURV(2,NDCURV+IEXTRA)
     *                      =X_CURV(IPOINT+1,ITHEOR)

C
                             Y_CURV(1,NDCURV+IEXTRA)
     *                      =EXPEXP_PROTON(INUCLI,IEXPER)
C
                             Y_CURV(2,NDCURV+IEXTRA)
     *                      =EXPEXP_PROTON(INUCLI,IEXPER)
C
                             L_CURV(1,NDCURV+IEXTRA)
     *                      =LABTEX
C
                             L_CURV(2,NDCURV+IEXTRA)
     *                      =LABTEX
C
                             NDPOIN(NDCURV+IEXTRA)=2
                             ITHICK(NDCURV+IEXTRA)=6
C
                             IF (ICOLOR(ITHEOR).EQ.44)
     *                           ICOLOR(NDCURV+IEXTRA)=71
C
                             IF (ICOLOR(ITHEOR).EQ.53)
     *                          ICOLOR(NDCURV+IEXTRA)=72
C
                          END IF
C
                          IF ((EXPEXP_PROTON(INUCLI,IEXPER)
     *                   .LE.  Y_CURV(IPOINT,ITHEOR))
     *                   .AND.(EXPEXP_PROTON(INUCLI,IEXPER)
     *                   .GT.  Y_CURV(IPOINT+1,ITHEOR)))
     *                    THEN
C
                             IEXTRA=IEXTRA+1
C
                             IF (IPOINT.EQ.1) THEN
                                 X_CURV(1,NDCURV+IEXTRA)
     *                          =X_CURV(IPOINT,ITHEOR)
                             ELSE
                                 X_CURV(1,NDCURV+IEXTRA)
     *                          =X_CURV(IPOINT-1,ITHEOR)
                             END IF
C
                             X_CURV(2,NDCURV+IEXTRA)
     *                      =X_CURV(IPOINT+1,ITHEOR)
C
                             Y_CURV(1,NDCURV+IEXTRA)
     *                      =EXPEXP_PROTON(INUCLI,IEXPER)
C
                             Y_CURV(2,NDCURV+IEXTRA)
     *                      =EXPEXP_PROTON(INUCLI,IEXPER)
C
                             L_CURV(1,NDCURV+IEXTRA)
     *                      =LABTEX
C
                             L_CURV(2,NDCURV+IEXTRA)
     *                      =LABTEX
C
                             NDPOIN(NDCURV+IEXTRA)=2
                             ITHICK(NDCURV+IEXTRA)=6
C
                             IF (ICOLOR(ITHEOR).EQ.44)
     *                           ICOLOR(NDCURV+IEXTRA)=71
C
                             IF (ICOLOR(ITHEOR).EQ.53)
     *                           ICOLOR(NDCURV+IEXTRA)=72
C
                          END IF
C
                       END DO !IPOINT=1,NDPOIN(ITHEOR)-1
C
                   END IF
C
                END DO
             END DO
C
             IF (IEXTRA.NE.LEVEXP_PROTON(INUCLI)) THEN
                 WRITE(0,'(''IEXTRA= '',I4,''LEVEXP= '',I4)')IEXTRA,
     *                                         LEVEXP_PROTON(INUCLI)
CID                 STOP 'IEXTRA.NE.LEVEXP_PROTON(INUCLI)'
             END IF
C
             NDCURV=NDCURV+IEXTRA   !LEVEXP_PROTON(INUCLI)
C
             ISOSPI=1
             IZ_FIX=NUMB_Z(INUCLI)
             IN_FIX=NUMB_N(INUCLI)
C
             MAIN_T='Energy Level Parametric Dependence'
             XTITLE='\\boldmath'//TITLES_LATEXS(JPARAM)
             YTITLE='Single Particle Energies (MeV)'
C
             IF (LOGWRI.GT.0) 
     *           WRITE(LOGFIL,'(9X,''Entering SPENER_PLOTIN '',
     *                                         ''for protons'')')
C
             CALL SPENER_PLOTIN(NDCURV,NDPOIN,X_CURV,Y_CURV,L_CURV,
     *                          MAIN_T,XTITLE,YTITLE,ICOLOR,ITHICK,
     *                          INUCLI,ISOSPI,IZ_FIX,IN_FIX,PARNAM,
     *                                 JPARAM,PARPOT_PRINTI,RMSVAL)
C
         END IF !ITAKNU(INUCLI).EQ.1
      END DO !INUCLI
C
C=======================================================================
C
   3  CONTINUE
C
      IF (IPRINT_NEUTRS.EQ.1) THEN
          CONTINUE   ! ==>> WE DO THE NEUTRONS
      ELSE
          GO TO 2    ! ==>> WE JUMP TO THE END
      END IF
C
C=======================================================================
C
C     NO MINIMIZATION - NEUTRONS
C
      IF (LOGWRI.GT.0) 
     *    WRITE(LOGFIL,'(9x,''NO Minimization option'',/)')
C
      JPARAM=INDEXI
C
      DO INUCLI=1,LDNUCL
         IF (ITAKNU(INUCLI).EQ.1) THEN
C
             DO INDEX1=1,MAXMSH
C
                IF (LOGWRI.GT.0) THEN
                    WRITE(NOUTPT,'(/,80(''#''),/,''#'',T80,''#'',/,
     *                             ''#  INDEX1= '',I6,T80,''#'',/,
     *                             ''#'',T80,''#'',/,80(''#''))')INDEX1
                END IF
C
                XINITS(JPARAM)=XMIN_I(JPARAM)+(INDEX1-1)*STPMSH(JPARAM)
                PARPOT(JPARAM)=XINITS(JPARAM)
C
                IDEFCN=0
                I_MODE=0
                I_FLAG=1
C
                CALL WS_RUN(INUCLI,CHISQU_CHOICE,IABORT,I_MODE,I_FLAG,
     *                                    CHISQU_PROTON,CHISQU_NEUTRS)
C
C               Updating the parameter vector for printing purposes
C
                DO KPARAM=1,NDPARS
                   PARPOT_PRINTI(KPARAM,INDEX1)=PARPOT(KPARAM)
                END DO
C
                RMSVAL(INDEX1,INUCLI)=SQRT(CHIWEI_NEUTRS)
C
C               Storing the values of the levels for each nucleus 
C               for the plotting output (depending on the reference,
C                                                         see below)
C
                NDCURV=LEVTHE_NEUTRS
C
C               Storing the reference order of the levels
C
                IF (INDEX1.EQ.1) THEN
                    DO IREFER=1,NDCURV
                       LABREF(IREFER)=LABTHE_NEUTRS(IREFER)
                    END DO
                END IF
C
                DO ICURVE=1,NDCURV
C
                   X_CURV(INDEX1,ICURVE)=XINITS(JPARAM)
C
                   DO ITHEOR=1,NDCURV
                      IF (LABTHE_NEUTRS(ITHEOR).EQ.LABREF(ICURVE))THEN
C
                          Y_CURV(INDEX1,ICURVE)
     *                   =ENETHE_NEUTRS(ITHEOR)
C
                          LABTHE=LABTHE_NEUTRS(ITHEOR)
                          CALL LATEXS_LABELS(LABTHE,LABTEX)
                          L_CURV(INDEX1,ICURVE)=LABTEX
C
                          LQNUMB=LWSSPH_NEUTRS(ITHEOR)
                          LPARIT=(-1)**LQNUMB
C
                          IF (LPARIT.EQ.+1) ICOLOR(ICURVE)=44
C
                          IF (LPARIT.EQ.-1) ICOLOR(ICURVE)=53
C
                          ITHICK(ICURVE)=4
C
                      END IF
                   END DO
C                   
                END DO !ICURVE
C
                DO I=1,NDCURV
                   NDPOIN(I)=MAXMSH
                END DO
C
             END DO !INDEX1
C
C            Looking the correspondance with the experiment
C
             IEXTRA=0
             DO IEXPER=1,LEVEXP_NEUTRS(INUCLI)
                DO ITHEOR=1,NDCURV
C
                   LABEXP=LABEXP_NEUTRS(INUCLI,IEXPER)
                   CALL LATEXS_LABELS(LABEXP,LABTEX)
C
                   IF (LABTEX.EQ.L_CURV(1,ITHEOR)) THEN
C
                       DO IPOINT=1,NDPOIN(ITHEOR)-1
C
                          IF ((EXPEXP_NEUTRS(INUCLI,IEXPER)
     *                   .GE.  Y_CURV(IPOINT,ITHEOR))
     *                   .AND.(EXPEXP_NEUTRS(INUCLI,IEXPER)
     *                   .LT.  Y_CURV(IPOINT+1,ITHEOR)))
     *                    THEN
C
                             IEXTRA=IEXTRA+1
C
                             IF (IPOINT.EQ.1) THEN
                                 X_CURV(1,NDCURV+IEXTRA)
     *                          =X_CURV(IPOINT,ITHEOR)
                             ELSE
                                 X_CURV(1,NDCURV+IEXTRA)
     *                          =X_CURV(IPOINT-1,ITHEOR)
                             END IF
C
                             X_CURV(2,NDCURV+IEXTRA)
     *                      =X_CURV(IPOINT+1,ITHEOR)
C
                             Y_CURV(1,NDCURV+IEXTRA)
     *                      =EXPEXP_NEUTRS(INUCLI,IEXPER)
C
                             Y_CURV(2,NDCURV+IEXTRA)
     *                      =EXPEXP_NEUTRS(INUCLI,IEXPER)
C
                             L_CURV(1,NDCURV+IEXTRA)
     *                      =LABTEX
C
                             L_CURV(2,NDCURV+IEXTRA)
     *                      =LABTEX
C
                             NDPOIN(NDCURV+IEXTRA)=2
                             ITHICK(NDCURV+IEXTRA)=6
C
                             IF (ICOLOR(ITHEOR).EQ.44)
     *                           ICOLOR(NDCURV+IEXTRA)=71
C
                             IF (ICOLOR(ITHEOR).EQ.53)
     *                          ICOLOR(NDCURV+IEXTRA)=72
C
                          END IF
C
                          IF ((EXPEXP_NEUTRS(INUCLI,IEXPER)
     *                   .LE.  Y_CURV(IPOINT,ITHEOR))
     *                   .AND.(EXPEXP_NEUTRS(INUCLI,IEXPER)
     *                   .GT.  Y_CURV(IPOINT+1,ITHEOR)))
     *                    THEN
C
                             IEXTRA=IEXTRA+1
C
                             IF (IPOINT.EQ.1) THEN
                                 X_CURV(1,NDCURV+IEXTRA)
     *                          =X_CURV(IPOINT,ITHEOR)
                             ELSE
                                 X_CURV(1,NDCURV+IEXTRA)
     *                          =X_CURV(IPOINT-1,ITHEOR)
                             END IF
C
                             X_CURV(2,NDCURV+IEXTRA)
     *                      =X_CURV(IPOINT+1,ITHEOR)
C
                             Y_CURV(1,NDCURV+IEXTRA)
     *                      =EXPEXP_NEUTRS(INUCLI,IEXPER)
C
                             Y_CURV(2,NDCURV+IEXTRA)
     *                      =EXPEXP_NEUTRS(INUCLI,IEXPER)
C
                             L_CURV(1,NDCURV+IEXTRA)
     *                      =LABTEX
C
                             L_CURV(2,NDCURV+IEXTRA)
     *                      =LABTEX
C
                             NDPOIN(NDCURV+IEXTRA)=2
                             ITHICK(NDCURV+IEXTRA)=6
C
                             IF (ICOLOR(ITHEOR).EQ.44)
     *                           ICOLOR(NDCURV+IEXTRA)=71
C
                             IF (ICOLOR(ITHEOR).EQ.53)
     *                           ICOLOR(NDCURV+IEXTRA)=72
C
                          END IF
C
                       END DO !IPOINT=1,NDPOIN(ITHEOR)-1
C
                   END IF
                END DO
             END DO
C
             IF (IEXTRA.NE.LEVEXP_NEUTRS(INUCLI)) THEN
                 WRITE(0,'(''IEXTRA= '',I4,''LEVEXP= '',I4)')IEXTRA,
     *                                         LEVEXP_NEUTRS(INUCLI)
CID                 STOP 'IEXTRA.NE.LEVEXP_NEUTRS(INUCLI)'
             END IF
C
             NDCURV=NDCURV+IEXTRA   !LEVEXP_NEUTRS(INUCLI)
C
             ISOSPI=0
             IZ_FIX=NUMB_Z(INUCLI)
             IN_FIX=NUMB_N(INUCLI)
C
             MAIN_T='Energy Level Parametric Dependence'
             XTITLE='\\boldmath'//TITLES_LATEXS(JPARAM)
             YTITLE='Single Particle Energies (MeV)'
C
             IF (LOGWRI.GT.0) 
     *           WRITE(LOGFIL,'(9X,''Entering SPENER_PLOTIN '',
     *                                         ''for neutrons'')')
C
             CALL SPENER_PLOTIN(NDCURV,NDPOIN,X_CURV,Y_CURV,L_CURV,
     *                          MAIN_T,XTITLE,YTITLE,ICOLOR,ITHICK,
     *                          INUCLI,ISOSPI,IZ_FIX,IN_FIX,PARNAM,
     *                                 JPARAM,PARPOT_PRINTI,RMSVAL)
C
         END IF !ITAKNU(INUCLI).EQ.1
      END DO !INUCLI
C
C=======================================================================
C
      IF (ICOUNT_OFPARS.EQ.0) GO TO 2 ! ==>> NO Minimization required
C                                               => we jump to the end
C
C=======================================================================
C
   1  CONTINUE
C
C=======================================================================
C
C     MINIMIZATION - PROTONS / NEUTRONS
C
      IF (LOGWRI.GT.0) 
     *    WRITE(LOGFIL,'(9x,''Minimization option'',/)')
C
      JPARAM=INDEXI
C
      NCOUNT_OFMESH=0
      NTOTAL_OFMESH=MAXMSH
C
      DO INDEX1=1,MAXMSH
C
         IF (LOGWRI.GT.0) THEN
             WRITE(NOUTPT,'(/,80(''#''),/,''#'',T80,''#'',/,
     *                       ''#  INDEX1= '',I6,T80,''#'',/,
     *                       ''#'',T80,''#'',/,80(''#''))')INDEX1
         END IF
C
         XINITS(JPARAM)=XMIN_I(JPARAM)+(INDEX1-1)*STPMSH(JPARAM)
         PARPOT(JPARAM)=XINITS(JPARAM)
C
         NCOUNT_OFMESH=NCOUNT_OFMESH+1
C
         CALL LMMESH(NEWSED,ARGMNT,I_SEED,LDRAND,LEVNUM,
     *               NTOTAL_OFMESH,NCOUNT_OFMESH,TITMSH,
     *               ICMESH,ICOUNT_OFPARS,CHITOT_MINAUX,
     *                      CHISQU_PROAUX,CHISQU_NEUAUX,
     *                             PARPOT_PRINTI,NDLAST)
C
C        Updating the parameter vector        
C
         DO KPARAM=1,NDPARS
            PARPOT(KPARAM)=PARPOT_PRINTI(KPARAM,NCOUNT_OFMESH)
         END DO
C_______________________________________________________________________
C
C        Once we have the vector of parameters, we enter to WS_RUN
C        to calculate the S. P. Energy Levels that will be plotted
C
         WRITE(NOUTPT,'(/,''Entering to the final calculation'',/)')
C
         DO INUCLI=1,LDNUCL
CID            IF (ITAKNU(INUCLI).EQ.1) THEN
C
                IZ_FIX=NUMB_Z(INUCLI)
                IN_FIX=NUMB_N(INUCLI)
C
                IF (LOGWRI.GT.4) THEN
                    WRITE(LOGFIL,'(12X,''Entering WS_RUN from '',
     *                                 ''SPENER_OFPARA '',
     *                                 ''INUCLI='',I2)') INUCLI
                END IF
C
                I_MODE=0
                I_FLAG=1
C
                CALL WS_RUN(INUCLI,CHISQU_CHOICE,IABORT,I_MODE,
     *                      I_FLAG,CHISQU_PROTON,CHISQU_NEUTRS)
C
C               Storing the values of the levels for each nucleus 
C               for the plotting output (depending on the reference,
C                                                         see below)
C
                NDCURV_PROTON=LEVTHE_PROTON
                NDCURV_NEUTRS=LEVTHE_NEUTRS
C
                RMSVAL_PROTON(INDEX1,INUCLI)=SQRT(CHIWEI_PROTON)
                RMSVAL_NEUTRS(INDEX1,INUCLI)=SQRT(CHIWEI_NEUTRS)
C
C               Storing the reference order of the levels for protons
C
                IF (INDEX1.EQ.1) THEN
                    DO IREFER=1,NDCURV_PROTON
                       LABREF_PROTON(IREFER)=LABTHE_PROTON(IREFER)
                    END DO
                END IF
C
                DO ICURVE=1,NDCURV_PROTON
C
                   X_CURV_AUXPRO(INDEX1,ICURVE,INUCLI)=XINITS(JPARAM)
C
                   DO ITHEOR=1,NDCURV_PROTON
                      IF (LABTHE_PROTON(ITHEOR)
     *               .EQ. LABREF_PROTON(ICURVE))THEN
C
                          Y_CURV_AUXPRO(INDEX1,ICURVE,INUCLI)
     *                   =ENETHE_PROTON(ITHEOR)
C
                          LABTHE=LABTHE_PROTON(ITHEOR)
                          CALL LATEXS_LABELS(LABTHE,LABTEX)
                          L_CURV_AUXPRO(INDEX1,ICURVE,INUCLI)=LABTEX
C
                          LQNUMB=LWSSPH_PROTON(ITHEOR)
                          LPARIT=(-1)**LQNUMB
C
                          IF (LPARIT.EQ.+1) 
     *                        ICOLOR_AUXPRO(ICURVE,INUCLI)=44
C
                          IF (LPARIT.EQ.-1) 
     *                        ICOLOR_AUXPRO(ICURVE,INUCLI)=53
C
                          ITHICK_AUXPRO(ICURVE,INUCLI)=4
C
                      END IF
                   END DO
C                   
                END DO !ICURVE
C
                DO I=1,NDCURV_PROTON
                   NDPOIN_PROTON(I)=MAXMSH
                END DO
C
C               the same for neutrons
C
                IF (INDEX1.EQ.1) THEN
                    DO IREFER=1,NDCURV_NEUTRS
                       LABREF_NEUTRS(IREFER)=LABTHE_NEUTRS(IREFER)
                    END DO
                END IF
C
                DO ICURVE=1,NDCURV_NEUTRS
C
                   X_CURV_AUXNEU(INDEX1,ICURVE,INUCLI)=XINITS(JPARAM)
C
                   DO ITHEOR=1,NDCURV_NEUTRS
                      IF (LABTHE_NEUTRS(ITHEOR)
     *               .EQ. LABREF_NEUTRS(ICURVE))THEN
C
                          Y_CURV_AUXNEU(INDEX1,ICURVE,INUCLI)
     *                   =ENETHE_NEUTRS(ITHEOR)
C
                          LABTHE=LABTHE_NEUTRS(ITHEOR)
                          CALL LATEXS_LABELS(LABTHE,LABTEX)
                          L_CURV_AUXNEU(INDEX1,ICURVE,INUCLI)=LABTEX
C
                          LQNUMB=LWSSPH_NEUTRS(ITHEOR)
                          LPARIT=(-1)**LQNUMB
C
                          IF (LPARIT.EQ.+1)
     *                        ICOLOR_AUXNEU(ICURVE,INUCLI)=44
C
                          IF (LPARIT.EQ.-1)
     *                        ICOLOR_AUXNEU(ICURVE,INUCLI)=53
C
                          ITHICK_AUXNEU(ICURVE,INUCLI)=4
C
                      END IF
                   END DO
C
                END DO !ICURVE
C
                DO I=1,NDCURV_NEUTRS
                   NDPOIN_NEUTRS(I)=MAXMSH
                END DO
C
CID             END IF !ITAKNU(INUCLI).EQ.1
         END DO !INUCLI
C
      END DO !INDEX1
C
C=======================================================================
C
C     PRINTING PROTONS
C
      IF (IPRINT_PROTON.EQ.1) THEN
          CONTINUE
      ELSE
          GO TO 14
      END IF
C
      DO INUCLI=1,LDNUCL
CID         IF (ITAKNU(INUCLI).EQ.1) THEN
C
C            Looking the correspondance with the experiment
C
             NDCURV=NDCURV_PROTON
C
             IEXTRA=0
             DO IEXPER=1,LEVEXP_PROTON(INUCLI)
                DO ITHEOR=1,NDCURV
C
                   LABEXP=LABEXP_PROTON(INUCLI,IEXPER)
                   CALL LATEXS_LABELS(LABEXP,LABTEX)
C
                   IF (LABTEX.EQ.L_CURV_AUXPRO(1,ITHEOR,INUCLI)) THEN
C
                       DO IPOINT=1,NDPOIN_PROTON(ITHEOR)-1
C
                          IF ((EXPEXP_PROTON(INUCLI,IEXPER)
     *                   .GE.  Y_CURV_AUXPRO(IPOINT,ITHEOR,INUCLI))
     *                   .AND.(EXPEXP_PROTON(INUCLI,IEXPER)
     *                   .LT.  Y_CURV_AUXPRO(IPOINT+1,ITHEOR,INUCLI)))
     *                    THEN
C
                             IEXTRA=IEXTRA+1
C
                             IF (IPOINT.EQ.1) THEN
                                 X_CURV_AUXPRO(1,NDCURV+IEXTRA,INUCLI)
     *                          =X_CURV_AUXPRO(IPOINT,ITHEOR,INUCLI)
                             ELSE
                                 X_CURV_AUXPRO(1,NDCURV+IEXTRA,INUCLI)
     *                          =X_CURV_AUXPRO(IPOINT-1,ITHEOR,INUCLI)
                             END IF
C
                             X_CURV_AUXPRO(2,NDCURV+IEXTRA,INUCLI)
     *                      =X_CURV_AUXPRO(IPOINT+1,ITHEOR,INUCLI)
C
                             Y_CURV_AUXPRO(1,NDCURV+IEXTRA,INUCLI)
     *                      =EXPEXP_PROTON(INUCLI,IEXPER)
C
                             Y_CURV_AUXPRO(2,NDCURV+IEXTRA,INUCLI)
     *                      =EXPEXP_PROTON(INUCLI,IEXPER)
C
                             L_CURV_AUXPRO(1,NDCURV+IEXTRA,INUCLI)
     *                      =LABTEX
C
                             L_CURV_AUXPRO(2,NDCURV+IEXTRA,INUCLI)
     *                      =LABTEX
C
                             NDPOIN_PROTON(NDCURV+IEXTRA)=2
                             ITHICK_AUXPRO(NDCURV+IEXTRA,INUCLI)=6
C
                             IF (ICOLOR_AUXPRO(ITHEOR,INUCLI).EQ.44)
     *                           ICOLOR_AUXPRO(NDCURV+IEXTRA,INUCLI)=71
C
                             IF (ICOLOR_AUXPRO(ITHEOR,INUCLI).EQ.53)
     *                          ICOLOR_AUXPRO(NDCURV+IEXTRA,INUCLI)=72
C
                          END IF
C
                          IF ((EXPEXP_PROTON(INUCLI,IEXPER)
     *                   .LE.  Y_CURV_AUXPRO(IPOINT,ITHEOR,INUCLI))
     *                   .AND.(EXPEXP_PROTON(INUCLI,IEXPER)
     *                   .GT.  Y_CURV_AUXPRO(IPOINT+1,ITHEOR,INUCLI)))
     *                    THEN
C
                             IEXTRA=IEXTRA+1
C
                             IF (IPOINT.EQ.1) THEN
                                 X_CURV_AUXPRO(1,NDCURV+IEXTRA,INUCLI)
     *                          =X_CURV_AUXPRO(IPOINT,ITHEOR,INUCLI)
                             ELSE
                                 X_CURV_AUXPRO(1,NDCURV+IEXTRA,INUCLI)
     *                          =X_CURV_AUXPRO(IPOINT-1,ITHEOR,INUCLI)
                             END IF
C
                             X_CURV_AUXPRO(2,NDCURV+IEXTRA,INUCLI)
     *                      =X_CURV_AUXPRO(IPOINT+1,ITHEOR,INUCLI)
C
                             Y_CURV_AUXPRO(1,NDCURV+IEXTRA,INUCLI)
     *                      =EXPEXP_PROTON(INUCLI,IEXPER)
C
                             Y_CURV_AUXPRO(2,NDCURV+IEXTRA,INUCLI)
     *                      =EXPEXP_PROTON(INUCLI,IEXPER)
C
                             L_CURV_AUXPRO(1,NDCURV+IEXTRA,INUCLI)
     *                      =LABTEX
C
                             L_CURV_AUXPRO(2,NDCURV+IEXTRA,INUCLI)
     *                      =LABTEX
C
                             NDPOIN_PROTON(NDCURV+IEXTRA)=2
                             ITHICK_AUXPRO(NDCURV+IEXTRA,INUCLI)=6
C
                             IF (ICOLOR_AUXPRO(ITHEOR,INUCLI).EQ.44)
     *                           ICOLOR_AUXPRO(NDCURV+IEXTRA,INUCLI)=71
C
                             IF (ICOLOR_AUXPRO(ITHEOR,INUCLI).EQ.53)
     *                          ICOLOR_AUXPRO(NDCURV+IEXTRA,INUCLI)=72
C
                          END IF
C
                       END DO !IPOINT=1,NDPOIN(ITHEOR)-1
C
                   END IF
C
                END DO
             END DO
C
             IF (IEXTRA.NE.LEVEXP_PROTON(INUCLI)) THEN
                 WRITE(0,'(''IEXTRA= '',I4,''LEVEXP_PROTON= '',I4)')
     *                                IEXTRA,LEVEXP_PROTON(INUCLI)
CID                 STOP 'IEXTRA.NE.LEVEXP_PROTON(INUCLI)'
             END IF
C
             NDCURV=NDCURV+IEXTRA   !LEVEXP_PROTON(INUCLI)
C_______________________________________________________________________
C
C            And now, printing the results for each nucleus
C
             DO ICURVE=1,NDCURV
                DO IPOINT=1,NDPOIN_PROTON(ICURVE)
C
                   X_CURV(IPOINT,ICURVE)
     *            =X_CURV_AUXPRO(IPOINT,ICURVE,INUCLI)
C
                   Y_CURV(IPOINT,ICURVE)
     *            =Y_CURV_AUXPRO(IPOINT,ICURVE,INUCLI)
C
                   L_CURV(IPOINT,ICURVE)
     *            =L_CURV_AUXPRO(IPOINT,ICURVE,INUCLI)
C
                END DO
C
                ICOLOR(ICURVE)=ICOLOR_AUXPRO(ICURVE,INUCLI)
                ITHICK(ICURVE)=ITHICK_AUXPRO(ICURVE,INUCLI)
C
             END DO
C
             ISOSPI=1
             IZ_FIX=NUMB_Z(INUCLI)
             IN_FIX=NUMB_N(INUCLI)
C
             MAIN_T='Energy Level Parametric Dependence'
             XTITLE='\\boldmath'//TITLES_LATEXS(JPARAM)
             YTITLE='Single Particle Energies (MeV)'
C
             CALL SPENER_PLOTIN(NDCURV,NDPOIN_PROTON,X_CURV,Y_CURV,
     *                          L_CURV,MAIN_T,XTITLE,YTITLE,ICOLOR,
     *                          ITHICK,INUCLI,ISOSPI,IZ_FIX,IN_FIX,
     *                                 PARNAM,JPARAM,PARPOT_PRINTI,
     *                                               RMSVAL_PROTON)
C
CID         END IF
      END DO
C
C=======================================================================
C
C     PRINTING NEUTRONS
C
  14  CONTINUE
C
      IF (IPRINT_NEUTRS.EQ.1) THEN
          CONTINUE
      ELSE
          GO TO 2
      END IF
C
      DO INUCLI=1,LDNUCL
CID         IF (ITAKNU(INUCLI).EQ.1) THEN
C
C            Looking the correspondance with the experiment
C
             NDCURV=NDCURV_NEUTRS
C
             IEXTRA=0
             DO IEXPER=1,LEVEXP_NEUTRS(INUCLI)
                DO ITHEOR=1,NDCURV
C
                   LABEXP=LABEXP_NEUTRS(INUCLI,IEXPER)
                   CALL LATEXS_LABELS(LABEXP,LABTEX)
C
                   IF (LABTEX.EQ.L_CURV_AUXNEU(1,ITHEOR,INUCLI)) THEN
C
                       DO IPOINT=1,NDPOIN_NEUTRS(ITHEOR)-1
C
                          IF ((EXPEXP_NEUTRS(INUCLI,IEXPER)
     *                   .GE.  Y_CURV_AUXNEU(IPOINT,ITHEOR,INUCLI))
     *                   .AND.(EXPEXP_NEUTRS(INUCLI,IEXPER)
     *                   .LT.  Y_CURV_AUXNEU(IPOINT+1,ITHEOR,INUCLI)))
     *                    THEN
C
                             IEXTRA=IEXTRA+1
C
                             IF (IPOINT.EQ.1) THEN
                                 X_CURV_AUXNEU(1,NDCURV+IEXTRA,INUCLI)
     *                          =X_CURV_AUXNEU(IPOINT,ITHEOR,INUCLI)
                             ELSE
                                 X_CURV_AUXNEU(1,NDCURV+IEXTRA,INUCLI)
     *                          =X_CURV_AUXNEU(IPOINT-1,ITHEOR,INUCLI)
                             END IF
C
                             X_CURV_AUXNEU(2,NDCURV+IEXTRA,INUCLI)
     *                      =X_CURV_AUXNEU(IPOINT+1,ITHEOR,INUCLI)
C
                             Y_CURV_AUXNEU(1,NDCURV+IEXTRA,INUCLI)
     *                      =EXPEXP_NEUTRS(INUCLI,IEXPER)
C
                             Y_CURV_AUXNEU(2,NDCURV+IEXTRA,INUCLI)
     *                      =EXPEXP_NEUTRS(INUCLI,IEXPER)
C
                             L_CURV_AUXNEU(1,NDCURV+IEXTRA,INUCLI)
     *                      =LABTEX
C
                             L_CURV_AUXNEU(2,NDCURV+IEXTRA,INUCLI)
     *                      =LABTEX
C
                             NDPOIN(NDCURV+IEXTRA)=2
                             ITHICK_AUXNEU(NDCURV+IEXTRA,INUCLI)=6
C
                             IF (ICOLOR_AUXNEU(ITHEOR,INUCLI).EQ.44)
     *                           ICOLOR_AUXNEU(NDCURV+IEXTRA,INUCLI)=71
C
                             IF (ICOLOR_AUXNEU(ITHEOR,INUCLI).EQ.53)
     *                          ICOLOR_AUXNEU(NDCURV+IEXTRA,INUCLI)=72
C
                          END IF
C
                          IF ((EXPEXP_NEUTRS(INUCLI,IEXPER)
     *                   .LE.  Y_CURV_AUXNEU(IPOINT,ITHEOR,INUCLI))
     *                   .AND.(EXPEXP_NEUTRS(INUCLI,IEXPER)
     *                   .GT.  Y_CURV_AUXNEU(IPOINT+1,ITHEOR,INUCLI)))
     *                    THEN
C
                             IEXTRA=IEXTRA+1
C
                             IF (IPOINT.EQ.1) THEN
                                 X_CURV_AUXNEU(1,NDCURV+IEXTRA,INUCLI)
     *                          =X_CURV_AUXNEU(IPOINT,ITHEOR,INUCLI)
                             ELSE
                                 X_CURV_AUXNEU(1,NDCURV+IEXTRA,INUCLI)
     *                          =X_CURV_AUXNEU(IPOINT-1,ITHEOR,INUCLI)
                             END IF
C
                             X_CURV_AUXNEU(2,NDCURV+IEXTRA,INUCLI)
     *                      =X_CURV_AUXNEU(IPOINT+1,ITHEOR,INUCLI)
C
                             Y_CURV_AUXNEU(1,NDCURV+IEXTRA,INUCLI)
     *                      =EXPEXP_NEUTRS(INUCLI,IEXPER)
C
                             Y_CURV_AUXNEU(2,NDCURV+IEXTRA,INUCLI)
     *                      =EXPEXP_NEUTRS(INUCLI,IEXPER)
C
                             L_CURV_AUXNEU(1,NDCURV+IEXTRA,INUCLI)
     *                      =LABTEX
C
                             L_CURV_AUXNEU(2,NDCURV+IEXTRA,INUCLI)
     *                      =LABTEX
C
                             NDPOIN_NEUTRS(NDCURV+IEXTRA)=2
                             ITHICK_AUXNEU(NDCURV+IEXTRA,INUCLI)=6
C
                             IF (ICOLOR_AUXNEU(ITHEOR,INUCLI).EQ.44)
     *                           ICOLOR_AUXNEU(NDCURV+IEXTRA,INUCLI)=71
C
                             IF (ICOLOR_AUXNEU(ITHEOR,INUCLI).EQ.53)
     *                          ICOLOR_AUXNEU(NDCURV+IEXTRA,INUCLI)=72
C
                          END IF
C
                       END DO !IPOINT=1,NDPOIN(ITHEOR)-1
C
                   END IF
C
                END DO
             END DO
C
             IF (IEXTRA.NE.LEVEXP_NEUTRS(INUCLI)) THEN
                 WRITE(0,'(''IEXTRA= '',I4,''LEVEXP_NEUTRS= '',I4)')
     *                                IEXTRA,LEVEXP_NEUTRS(INUCLI)
CID                 STOP 'IEXTRA.NE.LEVEXP_NEUTRS(INUCLI)'
             END IF
C
             NDCURV=NDCURV+IEXTRA   !LEVEXP_NEUTRS(INUCLI)
C_______________________________________________________________________
C
C            And now, printing the results for each nucleus
C
             DO ICURVE=1,NDCURV
                DO IPOINT=1,NDPOIN_NEUTRS(ICURVE)
C
                   X_CURV(IPOINT,ICURVE)
     *            =X_CURV_AUXNEU(IPOINT,ICURVE,INUCLI)
C
                   Y_CURV(IPOINT,ICURVE)
     *            =Y_CURV_AUXNEU(IPOINT,ICURVE,INUCLI)
C
                   L_CURV(IPOINT,ICURVE)
     *            =L_CURV_AUXNEU(IPOINT,ICURVE,INUCLI)
C
                END DO
C
                ICOLOR(ICURVE)=ICOLOR_AUXNEU(ICURVE,INUCLI)
                ITHICK(ICURVE)=ITHICK_AUXNEU(ICURVE,INUCLI)
C
             END DO
C
             ISOSPI=0
             IZ_FIX=NUMB_Z(INUCLI)
             IN_FIX=NUMB_N(INUCLI)
C
             MAIN_T='Energy Level Parametric Dependence'
             XTITLE='\\boldmath'//TITLES_LATEXS(JPARAM)
             YTITLE='Single Particle Energies (MeV)'
C
             CALL SPENER_PLOTIN(NDCURV,NDPOIN_NEUTRS,X_CURV,Y_CURV,
     *                          L_CURV,MAIN_T,XTITLE,YTITLE,ICOLOR,
     *                          ITHICK,INUCLI,ISOSPI,IZ_FIX,IN_FIX,
     *                                 PARNAM,JPARAM,PARPOT_PRINTI,
     *                                               RMSVAL_NEUTRS)
C
CID         END IF
      END DO
C
C=======================================================================
C
   2  CONTINUE
C
C=======================================================================
C
      WRITE(LOGFIL,'(/,''Exiting SPENER_OFPARA'',/)')
C
C=======================================================================
C
      RETURN
      END
C
C=======================================================================
C=======================================================================
C
      SUBROUTINE SPENER_PLOTIN(NDCURV,NDPOIN,X_CURV,Y_CURV,L_CURV,
     *                         MAIN_T,XTITLE,YTITLE,ICOLOR,ITHICK,
     *                         INUCLI,ISOSPI,IZ_FIX,IN_FIX,PARNAM,
     *                         JPARAM,PARPOT_PRINTI,RMSVAL)
C
      INCLUDE 'MATDIM/NDMESH.f'
      INCLUDE 'MATDIM/NDSPEC.f'
      INCLUDE 'MATDIM/NDPARS.f'
      INCLUDE 'MATDIM/NDNUCL.f'
      INCLUDE 'MATDIM/NDCOLO.f'
      INCLUDE 'MATDIM/NDTITL.f'
C
      PARAMETER
     *         (NDMES2=NDMESH*NDMESH)
      CHARACTER
     *          NUCSYM*6,INPSYM*6,TYPCHI*6,VERSIO*3
      CHARACTER
     *          FILNAM*256,L_CURV*100,NUCNAM*06,TITLES*12,
     *          ISONAM*002,TEXLAM*020,TITPAR*13,FITNUC*02,
     *          TITLES_LATEXS*050,NUCNAM_LATEXS*010
      CHARACTER
     *          MAIN_T*256,XTITLE*256,YTITLE*256
      CHARACTER
     *          FILNAM_IFITED*14,FILNAM_NUCFIT*30,
     *          FITNUC_AUXILI*02,PARNAM*8
      CHARACTER
     *          AUXTX1*30,AUXTX2*30,AUXTX3*30,AUXTX4*30,
     *          AUXTX5*30,AUXTX6*30
C
      DIMENSION
     *          NDPOIN(1:NDSPEC),
     *          X_CURV(1:NDMESH,1:NDSPEC),
     *          Y_CURV(1:NDMESH,1:NDSPEC),
     *          L_CURV(1:NDMESH,1:NDSPEC),
     *          ICOLOR(1:NDSPEC),
     *          ITHICK(1:NDSPEC)
      DIMENSION
     *          FITNUC(1:NDNUCL),
     *          FITNUC_AUXILI(1:NDNUCL)
      DIMENSION
     *          PARPOT_PRINTI(1:NDPARS,1:NDMES2),
     *          RMSVAL(1:NDSPEC,1:NDNUCL)
      COMMON
     *       /DENSTR/ IFDENS
      COMMON
     *       /TENSOR/ IFTENS
      COMMON
     *       /TENOPT/ ICENTT,ISORBT,ITENSR 
      COMMON
     *       /DEPPAI/ IFDEEP,IFPRON
      COMMON
     *       /IFPAIR/ IF_PAI,IFGPAI,IFDPAI
      COMMON 
     *       /CHITAK/ IF_SPE,IF_RAD,IF_GAP,IF_FER,IF_DEN,IF_RHO,IF_INV,
     *                LDCHI2,TYPCHI(1:NDNUCL),ITKCHI(1:NDNUCL)
      COMMON
     *       /IFLAMB/ IFPAR1,IFPAR2,IFPAR3,IFPAR4,
     *                IFPAR5,IFPAR6,IFPAR7,IFPAR8,
     *                IFPA09,IFPA10,IFPA11,IFPA12
      COMMON
     *       /PARCOR/ IFCORR,IFRCVC,IFRCAC,IFVCAC,IFRSVS,IFRSAS,IFVSAS
      COMMON
     *       /IFKAPP/ IFK_VC,IFK_VS,IFK_RC,IFK_RS,IFK_AC,IFK_AS
      COMMON
     *       /POTPOT/ PARPOT(1:NDPARS)
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
     *       /SEARCH/ IFTAKE(1:NDPARS)
      COMMON
     *       /VERSIN/ VERSIO
      COMMON
     *       /PARAMT_TITLES/ TITLES(1:NDTITL),
     *                       TITLES_LATEXS(1:NDTITL),
     *                       NUCNAM_LATEXS(1:NDNUCL)
      COMMON
     *       /PRFILE/ ISCREN,LOGFIL,LOGCHP,LOGCHN,LOGHIP,LOGHIN,
     *                LOGPRO,LOGNEU,LOGWRI,LSCREN,IRESUL,ICONVE,
     *                NOUTPT,LOGAUX,LOGENE,LOGRAD,LOGGRD,LOGMSH,
     *                                                   LOGBIS
      DATA
     *       FITNUC(01) / 'Ox' /,
     *       FITNUC(02) / 'C0' /,
     *       FITNUC(03) / 'C8' /,
     *       FITNUC(04) / 'Ni' /,
     *       FITNUC(05) / 'Zr' /,
     *       FITNUC(06) / 'Sn' /,
     *       FITNUC(07) / 'Gd' /,
     *       FITNUC(08) / 'Pb' /
C
C=======================================================================
C
      IF (LOGWRI.GT.0) WRITE(LOGFIL,'(/,''Entering SPENER_PLOTIN'')')
      IF (LOGWRI.GT.0) WRITE(LOGFIL,'(''INUCLI= '',I4)')INUCLI
      IF (LOGWRI.GT.0) 
     *    WRITE(LOGFIL,'(/,'' IZ= '',I3,'' IN= '',I3,/)') 
     *                  NUMB_Z(INUCLI),NUMB_N(INUCLI)
C
C=======================================================================
C
      WRITE(0,'(/,''Entering SPENER_PLOTIN'')')
      WRITE(0,'(''INUCLI= '',I4,'' ISOSPI= '',I2)')INUCLI,ISOSPI
      WRITE(0,'(/,'' IZ= '',I3,'' IN= '',I3,/)')
     *                  NUMB_Z(INUCLI),NUMB_N(INUCLI)
C
C=======================================================================
C
      YMINIM=+1.0E+10
      YMAXIM=-1.0E+10
      
      DO I_CURV=1,NDCURV
         DO IPOINT=1,NDPOIN(I_CURV)
            IF (Y_CURV(IPOINT,I_CURV).LT.YMINIM) THEN
                YMINIM=Y_CURV(IPOINT,I_CURV)
            END IF
            IF (Y_CURV(IPOINT,I_CURV).GT.YMAXIM) THEN
                YMAXIM=Y_CURV(IPOINT,I_CURV)
            END IF
         END DO
      END DO
C
C=======================================================================
C
      XMINIM=+1.0E+10
      XMAXIM=-1.0E+10
C
      DO I_CURV=1,NDCURV
         DO IPOINT=1,NDPOIN(I_CURV)
            IF (X_CURV(IPOINT,I_CURV).LT.XMINIM) THEN
                XMINIM=X_CURV(IPOINT,I_CURV)
            END IF
            IF (X_CURV(IPOINT,I_CURV).GT.XMAXIM) THEN
                XMAXIM=X_CURV(IPOINT,I_CURV)
            END IF
         END DO
      END DO
C
C=======================================================================
C
      NRESUL=60
C      
      IF (ISOSPI.EQ.1) ISONAM='-P'
      IF (ISOSPI.EQ.0) ISONAM='-N'
C
C=======================================================================
C 
      KACTIV=0
      DO JNUCLI=1,LDNUCL
         IF (ITAKNU(JNUCLI).EQ.1) THEN
             KACTIV=KACTIV+1
             FITNUC_AUXILI(KACTIV)=FITNUC(JNUCLI)
         END IF
      END DO
C
      I1=3*KACTIV
C
      WRITE(FILNAM_NUCFIT,'(''_'',<NUCACT>(A2,''-''))')
     *                      (FITNUC_AUXILI(I),I=1,NUCACT)
C
C=======================================================================
C
      WRITE(NUCNAM,'(A6)') NUCSYM(INUCLI)
C
      J1=1
      J2=6
      IF (NUCNAM(1:1).EQ.' ') J1=2
      IF (NUCNAM(2:2).EQ.' ') J1=3
      IF (NUCNAM(3:3).EQ.' ') J1=4
C
      WRITE(0,*)INUCLI,NUCNAM(J1:J2)
C
C=======================================================================
C
      IF (IFDENS.EQ.0) THEN
          WRITE(FILNAM,'(''SPECurves/IFDENS-0/SPECur_'',A,A2,A,''_'',A8,
     *                  ''_IFDENS-'',I1,''_IFTENS-'',I1,''_IF-PAI-'',I1,
     *                  ''_IF-RAD-'',I1,''_IF-ORD-'',I1,
     *                  ''_IF-RHO-'',I1,''_IFDEEP-'',I1,
     *                  ''_IFPRON-'',I1,
     *                  ''_IFK-VC-'',I1,''_IFK-VS-'',I1,
     *                  ''_IFK-RC-'',I1,''_IFK-RS-'',I1,
     *                  ''_IFK-AC-'',I1,''_IFK-AS-'',I1,
     *                                   ''_'',A3,''.dat'')')
     *
     *               NUCNAM(J1:J2),ISONAM,FILNAM_NUCFIT(01:I1),PARNAM,
     *               IFDENS,IFTENS,IF_PAI,IF_RAD,IF_INV,IF_RHO,
     *               IFDEEP,IFPRON,IFK_VC,IFK_VS,IFK_RC,IFK_RS,IFK_AC,
     *                                                  IFK_AS,VERSIO
      END IF
C
      WRITE(0,*)FILNAM
C
C=======================================================================
C
      IF (IFDENS.EQ.1 .AND. IFTENS.EQ.0) THEN
C
          WRITE(TEXLAM,'(''_LAMBDA'',SP,I2,I2,I2,I2)')
     *                            IFPAR1,IFPAR2,IFPAR3,IFPAR4
C
          WRITE(FILNAM,'(''SPECurves/IFDENS-1_IFTENS-0/SPECur_'',
     *                                          A,A2,A,''_'',A8,
     *                  ''_IFDENS-'',I1,''_IFTENS-'',I1,''_IF-PAI-'',I1,
     *                  ''_IF-RAD-'',I1,''_IF-ORD-'',I1,''_IF-RHO-'',I1,
     *                  ''_IFDEEP-'',I1,''_IFPRON-'',I1,
     *                  ''_IFK-VC-'',I1,''_IFK-VS-'',I1,
     *                  ''_IFK-RC-'',I1,''_IFK-RS-'',I1,
     *                  ''_IFK-AC-'',I1,''_IFK-AS-'',I1,
     *                                          A15,''_'',A3,''.dat'')')
     *
     *               NUCNAM(J1:J2),ISONAM,FILNAM_NUCFIT(01:I1),PARNAM,
     *               IFDENS,IFTENS,IF_PAI,IF_RAD,IF_INV,IF_RHO,
     *               IFDEEP,IFPRON,IFK_VC,IFK_VS,IFK_RC,IFK_RS,IFK_AC,
     *                                           IFK_AS,TEXLAM,VERSIO
      END IF
C
C=======================================================================
C
      IF (IFDENS.EQ.1 .AND. IFTENS.EQ.1) THEN
C
          WRITE(TEXLAM,'(''_LAMBDA'',SP,I2,I2,I2,I2)')
     *                            IFPAR1,IFPAR2,IFPAR3,IFPAR4
C
          WRITE(FILNAM,'(''SPECurves/IFDENS-1_IFTENS-1/SPECur_'',
     *                                          A,A2,A,''_'',A8,
     *                  ''_IFDENS-'',I1,''_IFTENS-'',I1,''_IF-PAI-'',I1,
     *                  ''_ICENTT-'',I1,''_ISORBT-'',I1,''_ITENSR-'',I1,
     *                  ''_IF-RAD-'',I1,''_IF-ORD-'',I1,''_IF-RHO-'',I1,
     *                  ''_IFDEEP-'',I1,''_IFPRON-'',I1,
     *                  ''_IFK-VC-'',I1,''_IFK-VS-'',I1,
     *                  ''_IFK-RC-'',I1,''_IFK-RS-'',I1,
     *                  ''_IFK-AC-'',I1,''_IFK-AS-'',I1,
     *                      A15,''_'',A3,''.dat'')')
     *
     *               NUCNAM(J1:J2),ISONAM,FILNAM_NUCFIT(01:I1),PARNAM,
     *               IFDENS,IFTENS,IF_PAI,ICENTT,ISORBT,ITENSR,
     *               IF_RAD,IF_INV,IF_RHO,IFDEEP,IFPRON,
     *               IFK_VC,IFK_VS,IFK_RC,IFK_RS,IFK_AC,IFK_AS,
     *                                           TEXLAM,VERSIO
      END IF
C
C=======================================================================
C
      OPEN(NRESUL,FILE=FILNAM,STATUS='UNKNOWN')
C
C=======================================================================
C
      WRITE(NRESUL,'(/,''<<CODEVERSIO>> VERSIO'')')
C
      WRITE(NRESUL,'(16X,A3)') VERSIO
C
C=======================================================================
C
      WRITE(NRESUL,'(/,''<<MAIN_TITLE>>'')')
C      
      WRITE(NRESUL,'(15X,A,/)') MAIN_T
C
C=======================================================================
C     
      WRITE(NRESUL,'(''<<NUCLIDINFO>>'',1X,''NUMB_Z'',2X,''NUMB_N'')')
C
      WRITE(NRESUL,'(15X,I3,6X,I3,/)') IZ_FIX,IN_FIX
C
      A_MASS=IZ_FIX+IN_FIX
C
C=======================================================================
C
      WRITE(NRESUL,'(''<<NoNUCACTIV>>'',1X,''NUCACT'')')
C      
      WRITE(NRESUL,'(15X,I3,/)') NUCACT
C
C=======================================================================
C
      WRITE(NRESUL,'(''<<EXPER_FILE>>'',1X,''IFDEEP'',2X
     *                                     ''IFPRON'')')
      
      WRITE(NRESUL,'(15X,2(I3,6X),/)') IFDEEP,IFPRON
C
C=======================================================================
C
      WRITE(NRESUL,'(''<<POTEN_INFO>>'',1X,''ISOSPI'',2X,''IFDENS'',
     *                                  2X,''IFTENS'',2X,''IF_PAI'')')
C      
      WRITE(NRESUL,'(15X,4(I3,5X),/)') ISOSPI,IFDENS,IFTENS,IF_PAI
C
C=======================================================================
C
      WRITE(NRESUL,'(''<<TENSR_INFO>>'',1X,''ICENTT'',2X
     *                                     ''ISORBT'',2X,''ITENSR'')')
C      
      WRITE(NRESUL,'(15X,3(I3,5X),/)') ICENTT,ISORBT,ITENSR
C
C=======================================================================
C
      WRITE(NRESUL,'(''<<TAKEN_CHI2>> IF_SPE  IF_RAD  IF_GAP  '',
     *               ''IF_FER  IF_DEN  IF_RHO  IF_INV'')')
      WRITE(NRESUL,'(15X,7(I3,5X),/)')IF_SPE,IF_RAD,IF_GAP,IF_FER,
     *                                       IF_DEN,IF_RHO,IF_INV
C
C=======================================================================
C
C     IFTAKE .EQ. -1: parameter eliminated via parametric correlation
C     IFTAKE .EQ. +1: parameter entering to the minimization
C     IFTAKE .EQ. +4: parameter considered constant but printed anyway
C
      INFTOP=0
      WRITE(NRESUL,'(''<<PARAM_INFO>> '',A8,$)')TITLES(JPARAM)(3:10)
      DO KPARAM=1,NDPARS
         IF (IFTAKE(KPARAM).EQ.-1 .OR. IFTAKE(KPARAM).EQ.1
     *  .OR. IFTAKE(KPARAM).EQ.4) THEN
             INFTOP=INFTOP+1
             WRITE(NRESUL,'(2X,A8,$)')TITLES(KPARAM)(3:10)
         END IF
      END DO
      WRITE(NRESUL,'()')
C
      WRITE(NRESUL,'(15X,I8,2X,I8)')NDPOIN(1),INFTOP
      DO I=1,NDPOIN(1)
         WRITE(NRESUL,'(15X,F8.4,$)')PARPOT_PRINTI(JPARAM,I)
         DO KPARAM=1,NDPARS
            IF (IFTAKE(KPARAM).EQ.-1 .OR. IFTAKE(KPARAM).EQ.1
     *     .OR. IFTAKE(KPARAM).EQ.4) THEN
                WRITE(NRESUL,'(2X,F8.4,$)')PARPOT_PRINTI(KPARAM,I)
            END IF
         END DO
         WRITE(NRESUL,'()')
      END DO
C
      WRITE(NRESUL,'()')
C
C=======================================================================
C
      WRITE(NRESUL,'(''<<PARAM_TITL>>'',/,14X,'' '',$)')
      DO KPARAM=1,NDPARS
         IF (IFTAKE(KPARAM).EQ.-1 .OR. IFTAKE(KPARAM).EQ.1
     *  .OR. IFTAKE(KPARAM).EQ.4) THEN
             WRITE(NRESUL,'(A50,''  '',$)')TITLES_LATEXS(KPARAM)
         END IF
      END DO
      WRITE(NRESUL,'()')
C
      WRITE(NRESUL,'()')
C
C=======================================================================
C
      WRITE(NRESUL,'(''<<RMS_VALUES>>'')')
      WRITE(NRESUL,'(15X,I6)')NDPOIN(1)
      DO I=1,NDPOIN(1)
         WRITE(NRESUL,'(15X,F8.4,2X,F8.4)')
     *         PARPOT_PRINTI(JPARAM,I),RMSVAL(I,INUCLI)
      END DO
C
C=======================================================================
C      
      WRITE(NRESUL,'(''<<XAXIS_TEXT>>'',/,15x,A,/)') XTITLE
C      
      WRITE(NRESUL,'(''<<YAXIS_TEXT>>'',/,15x,A,/)') YTITLE
C
C=======================================================================
C
      WRITE(NRESUL,'(''<<SIDE_TEXTS>>'')')
C
      AUXTX1=' '
      AUXTX2=' '
      AUXTX3=' '
      AUXTX4=' '
      AUXTX5=' '
      AUXTX6=' '
C
      IF (IFDENS.EQ.1) AUXTX1='${\\rm DENS}\\,\\vert\\,$'
      IF (IFTENS.EQ.1) AUXTX2='${\\rm TENS}\\,\\vert\\,$'
      IF (IF_PAI.EQ.1) AUXTX3='${\\rm PAIR}\\,\\vert\\,$'
      IF (IF_RAD.EQ.1) AUXTX4='${\\rm RADI}\\,\\vert\\,$'
      IF (IFDEEP.EQ.1) AUXTX5='${\\rm DEEP}\\,\\vert\\,$'
      IF (IFPRON.EQ.1) AUXTX6='${\\rm PRON}\\,\\vert\\,$'
C
      WRITE(NRESUL,'(15X,''~~~Id:'',6A)')AUXTX1,AUXTX2,AUXTX3,AUXTX4,
     *                                   AUXTX5,AUXTX6
C
      WRITE(NRESUL,'(15X,''Fitted: '',$)')
      DO KNUCLI=1,LDNUCL
         IF (ITAKNU(KNUCLI).EQ.1) THEN
             WRITE(NRESUL,'(A10,'','',$)') NUCNAM_LATEXS(KNUCLI)
         END IF
      END DO
      WRITE(NRESUL,'()')
C
      WRITE(NRESUL,'(15X,''Fixed: '',$)')
      DO KPARAM=1,NDPARS
         IF (IFTAKE(KPARAM).EQ.4) THEN
             WRITE(NRESUL,'(A,'' $='',F10.2,''$, '',$)') 
     *             TITLES_LATEXS(KPARAM),PARPOT_PRINTI(KPARAM,1)
         END IF
      END DO
      WRITE(NRESUL,'()')
C
C=======================================================================
C     
      WRITE(NRESUL,'(/,15X,30(''=''),/)')
          
      WRITE(NRESUL,'(''<<X_AX_PARAM>> XMINIM_FIG  XMAXIM_FIG  '',
     *                   ''X_STEP_FIG'')')
          
      WRITE(NRESUL,'(15X,3(F10.2,2X),/)')XMINIM,XMAXIM,(XMAXIM-XMINIM)
     
      WRITE(NRESUL,'(''<<Y_AX_PARAM>> YMINIM_FIG  YMAXIM_FIG  '',
     *                   ''Y_STEP_FIG'')')
          
      WRITE(NRESUL,'(15X,3(F10.2,2X),/)')
     *                         YMINIM,YMAXIM,(YMAXIM-YMINIM)
     
      WRITE(NRESUL,'(15X,30(''=''),/)')
C
C=======================================================================
C
      I_COLM=0
      I_LEFT=1
      I_RIGH=1
C
      WRITE(NRESUL,'(''<<WHATLEGEND>>'',1X,''COLM_LABEL'',2X,
     *               ''LEFT_LABEL'',2X,''RIGH_LABEL'')')
      WRITE(NRESUL,'(15X,I5,7X,I5,7X,I5)')I_COLM,I_LEFT,I_RIGH
C
C=======================================================================
C      
      IFBULL=0
      IFHIST=0
      IFOVER=0
      ISPLIT=0
C      
      WRITE(NRESUL,'(/,''<<NO_OF_CURV>>'',1X,''NUMB_CURVE'',2X,
     *                 ''IF_SPLITIN'',2X,''IF_BULLETS'',2X,
     *                 ''IF_HISTOGS'',2X,''IF_OVERFITING'')')
      WRITE(NRESUL,'(15X,I5,7X,I5,7X,I5,7X,I5,7X,I5)')NDCURV,ISPLIT,
     *                                         IFBULL,IFHIST,IFOVER
C
C=======================================================================
C
      ISTYLE=0
      I_TYPE=3
C
      DO I_CURV=1,NDCURV
C
         WRITE(NRESUL,'(''<<>>'')')
         WRITE(NRESUL,'(''<<CURVE_'',I4.4,''>>'',1X,
     *                  ''THICK_LINE'',2X,''STYLE_LINE'',2X,
     *                  ''COLOR_LINE'',2X,''TYPE_POINT'')')I_CURV
         WRITE(NRESUL,'(15X,I5,7X,I5,7X,I5,7X,I5)')ITHICK(I_CURV),
     *                                             ISTYLE,
     *                                             ICOLOR(I_CURV),
     *                                             I_TYPE
         WRITE(NRESUL,'(15X,''NUMB_POINT'',/,15X,I5)')NDPOIN(I_CURV)
C
         DO IPOINT=1,NDPOIN(I_CURV)
C
            WRITE(NRESUL,'(15X,F10.4,4X,E12.4,4X,A100)') 
     *                           X_CURV(IPOINT,I_CURV),
     *                           Y_CURV(IPOINT,I_CURV),
     *                           L_CURV(IPOINT,I_CURV)
C
         END DO
C
      END DO
C
C=======================================================================
C      
      WRITE(NRESUL,'(/,''<<GO_GETTHEM>>'',/)')
C
C=======================================================================
C
      RETURN
      END
C
C=======================================================================
C=======================================================================
C
      SUBROUTINE FITTIN_MONTEC(I_SEED,NEWSED)
C
      INCLUDE 'MATDIM/NDSPEC.f'
      INCLUDE 'MATDIM/NDPARS.f'
      INCLUDE 'MATDIM/NDNUCL.f'
      INCLUDE 'MATDIM/NDLEXP.f'
      INCLUDE 'MATDIM/NDIM_M.f'
      INCLUDE 'MATDIM/MAXFEV.f'
      INCLUDE 'MATDIM/NDLAST.f'
      INCLUDE 'MATDIM/NDBINS.f'
      INCLUDE 'MATDIM/NDMONT.f'
      INCLUDE 'MATDIM/NDCOLO.f'
      INCLUDE 'MATDIM/NDTITL.f'
C
      PARAMETER 
     *         (NDCURV=NDSPEC,LENGTH=256,INMODE=1,NPRINT=-1,ND_RMS=2)
C
      CHARACTER
     *          LABPRO_PSEUDO*6,LABNEU_PSEUDO*6,
     *          LABEXP_PROTON*6,LABEXP_NEUTRS*6,
     *          LABEXP_PROREA*6,LABEXP_NEUREA*6,
     *          LABTHE_PROTON*6,LABTHE_NEUTRS*6
      CHARACTER
     *          LABVEC_PROTON*6,LABVEC_NEUTRS*6,
     *          LABREF_PROTON*6,LABREF_NEUTRS*6
      CHARACTER
     *          INPSYM*6,NUCSYM*6,TYPCHI*6
      CHARACTER
     *          FILNAM*256,FILNAM_ENDING*13,
     *          LABTHE_PRINTG*06,TEXLAM*20,
     *          FILNA2*256,FITNUC*2,FILNAM_NUCFIT*30,
     *          FITNUC_AUXILI*30,FILNA3*256,VERSIO*3
      CHARACTER
     *          TITLES*12,TITLES_LATEXS*50,NUCNAM_LATEXS*10,LABELS*100,
     *          PARPOT_UNITSS*40,XUNITS*40
      CHARACTER
     *          FILPAR*300,PARNAM*8,FILRMS*300,FILSPE*300
      CHARACTER
     *          HISTIT*256,HX_TIT*256,HY_TIT*256,LABTIT*100,
     *          LABTIT_PARPOT*100,LABTIT_PROTON*100,LABTIT_NEUTRS*100
      CHARACTER
     *          SPELAB*6,CDUMMY*6,LABTHE*6,LABTEX*11
C
      DIMENSION
     *          LABREF_PROTON(1:NDSPEC),
     *          LABREF_NEUTRS(1:NDSPEC)
      DIMENSION
     *          ARGPAR(1:NDPARS),
     *          ARGPAR_TITLES(1:NDPARS),
     *          ARGPAR_LATEXS(1:NDPARS)
      DIMENSION
     *          FITNUC(1:NDNUCL),
     *          FITNUC_AUXILI(1:NDNUCL)
      DIMENSION
     *          LABELS(1:NDBINS,1:NDCURV),
     *          XHISTO(1:NDBINS,1:NDCURV),
     *          XHISTO_MINIMS(1:NDCURV),
     *          XHISTO_MAXIMS(1:NDCURV)
      DIMENSION
     *          PARPOT_MAXIMA(1:NDPARS),
     *          PARPOT_MINIMA(1:NDPARS)
      DIMENSION
     *          SPELAB(1:NDSPEC)
      DIMENSION
     *          BINSIZ(1:NDCURV),
     *          YHISTO_UNNORM(1:NDBINS,1:NDCURV),
     *          YHISTO_INTEGR(1:NDBINS,1:NDCURV),
     *          YHISTO_MAXIMS(1:NDBINS,1:NDCURV)
      DIMENSION
     *          XHISTO_PRINTI(1:NDBINS,1:NDCURV),
     *          YHISTO_PRINTI(1:NDBINS,1:NDCURV)
      DIMENSION
     *          YEXPEC(1:NDCURV),
     *          ICOLOR_HISTOG(1:NDCURV)
      DIMENSION
     *          LABTIT(1:ND_RMS),
     *          LABTIT_PARPOT(1:NDPARS),
     *          LABTIT_PROTON(1:NDSPEC),
     *          LABTIT_NEUTRS(1:NDSPEC)
      DIMENSION
     *          VECTOR_RMSVAL(1:ND_RMS),
     *          VECTOR_PARPOT(1:NDPARS),
     *          VECTOR_ENEPRO(1:NDSPEC),
     *          VECTOR_ENENEU(1:NDSPEC)
C
      COMMON
     *       /GAUSHE/ NGAUSS,NSHELL_PROTON,NSHELL_NEUTRS
      COMMON
     *       /CHITAK/ IF_SPE,IF_RAD,IF_GAP,IF_FER,IF_DEN,IF_RHO,IF_INV,
     *                LDCHI2,TYPCHI(1:NDNUCL),ITKCHI(1:NDNUCL)
      COMMON
     *       /IFLAMB/ IFPAR1,IFPAR2,IFPAR3,IFPAR4,
     *                IFPAR5,IFPAR6,IFPAR7,IFPAR8,
     *                IFPA09,IFPA10,IFPA11,IFPA12
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
     *       /RANRAN/ IRMFST,I_RAND,IRANDO
      COMMON
     *       /NUCLEU/ IZ_FIX,IN_FIX,ISOSPI 
      COMMON
     *       /SEARCH/ IFTAKE(1:NDPARS)    
      COMMON
     *       /POTPOT/ PARPOT(1:NDPARS)
      COMMON
     *       /VLIMIT/ VMISTR(1:NDPARS),
     *                VMIMIN(1:NDPARS),
     *                VMIMAX(1:NDPARS),
     *                VMISTP(1:NDPARS)
      COMMON
     *       /TENOPT/ ICENTT,ISORBT,ITENSR
     *       /IFPAIR/ IF_PAI,IFGPAI,IFDPAI
      COMMON
     *       /CHIGRD/ CHISQU_GRDNRM,
     *                CHISQU_GRADIE(1:NDPARS)
      COMMON
     *       /MINITE/ IDEFCN,ITECHI
     *       /P_N_LM/ IFPROT,IFNEUT,IFBOTH
     *       /FEVALS/ IFUNCT_EVALUS
      COMMON
     *       /STOPAR/ EPSLAS,LDLAST
     *       /MASSIV/ IMASIV_PRODUC
     *       /DEPPAI/ IFDEEP,IFPRON
      COMMON
     *       /DENSTR/ IFDENS
     *       /TENSOR/ IFTENS
      COMMON
     *       /IFKAPP/ IFK_VC,IFK_VS,IFK_RC,IFK_RS,IFK_AC,IFK_AS
      COMMON
     *       /MONTEC/ IFPSEU,IFPARA,LDMONT,LDBINS
      COMMON
     *       /PSEUDO_ZN_NUM/ IZ_PSE(1:NDNUCL),
     *                       IN_PSE(1:NDNUCL)
      COMMON
     *       /PSEUDO_LEVELS/ LEVNEU_PSEUDO(1:NDNUCL),
     *                       LEVPRO_PSEUDO(1:NDNUCL)
      COMMON
     *       /PSEUDO_PROTON/ ENEPRO_PSEUDO(1:NDLEXP,1:NDNUCL),
     *                       LABPRO_PSEUDO(1:NDLEXP,1:NDNUCL),
     *                       SIGPRO_PSEUDO(1:NDLEXP,1:NDNUCL),
     *                       IDEPRO_PSEUDO(1:NDLEXP,1:NDNUCL),
     *                       LEVTAK_PROTON(1:NDLEXP,1:NDNUCL)
      COMMON
     *       /PSEUDO_NEUTRS/ ENENEU_PSEUDO(1:NDLEXP,1:NDNUCL),
     *                       LABNEU_PSEUDO(1:NDLEXP,1:NDNUCL),
     *                       SIGNEU_PSEUDO(1:NDLEXP,1:NDNUCL),
     *                       IDENEU_PSEUDO(1:NDLEXP,1:NDNUCL),
     *                       LEVTAK_NEUTRS(1:NDLEXP,1:NDNUCL)
      COMMON
     *       /SIGEXP_PRONEU/ SIGEXP_PROTON(1:NDLEXP,1:NDNUCL),
     *                       SIGEXP_NEUTRS(1:NDLEXP,1:NDNUCL)
      COMMON
     *       /EX_ENE/ EXPEXP_PROTON(1:NDNUCL,1:NDLEXP),
     *                IDEGEX_PROTON(1:NDNUCL,1:NDLEXP),
     *
     *                EXPEXP_NEUTRS(1:NDNUCL,1:NDLEXP),
     *                IDEGEX_NEUTRS(1:NDNUCL,1:NDLEXP)
      COMMON
     *       /EX_LAB/ LABEXP_PROTON(1:NDNUCL,1:NDLEXP),
     *                LABEXP_NEUTRS(1:NDNUCL,1:NDLEXP)
      COMMON
     *       /DIMSPH/ LDSING_PROTON,LEVTHE_PROTON,
     *                LEVEXP_PROTON(1:NDNUCL),
     *                LDSING_NEUTRS,LEVTHE_NEUTRS,
     *                LEVEXP_NEUTRS(1:NDNUCL)
      COMMON
     *       /REALEX_PROTON/ EXPEXP_PROREA(1:NDNUCL,1:NDLEXP),
     *                       LABEXP_PROREA(1:NDNUCL,1:NDLEXP),
     *                       IDEGEX_PROREA(1:NDNUCL,1:NDLEXP),
     *                       LEVEXP_PROREA(1:NDNUCL)
      COMMON
     *       /REALEX_NEUTRS/ EXPEXP_NEUREA(1:NDNUCL,1:NDLEXP),
     *                       LABEXP_NEUREA(1:NDNUCL,1:NDLEXP),
     *                       IDEGEX_NEUREA(1:NDNUCL,1:NDLEXP),
     *                       LEVEXP_NEUREA(1:NDNUCL)
      COMMON
     *       /SPECTR/ ENETHE_PROTON(1:NDSPEC),
     *                ENETHE_NEUTRS(1:NDSPEC)
      COMMON
     *       /THELAB/ LABTHE_PROTON(1:NDSPEC),
     *                LABTHE_NEUTRS(1:NDSPEC)  
      COMMON
     *       /QUANTN/ NWSSPH_PROTON(1:NDSPEC),
     *                LWSSPH_PROTON(1:NDSPEC),
     *                JWSSPH_PROTON(1:NDSPEC),
     *
     *                NWSSPH_NEUTRS(1:NDSPEC),
     *                LWSSPH_NEUTRS(1:NDSPEC),
     *                JWSSPH_NEUTRS(1:NDSPEC)
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
     *       /VERSIN/ VERSIO
      COMMON
     *       /PARAMT_TITLES/ TITLES(1:NDTITL),
     *                       TITLES_LATEXS(1:NDTITL),
     *                       NUCNAM_LATEXS(1:NDNUCL)
      COMMON
     *       /PARUNI/ PARPOT_UNITSS(1:NDTITL)
C
      COMMON
     *       /PRFILE/ ISCREN,LOGFIL,LOGCHP,LOGCHN,LOGHIP,LOGHIN,
     *                LOGPRO,LOGNEU,LOGWRI,LSCREN,IRESUL,ICONVE,
     *                NOUTPT,LOGAUX,LOGENE,LOGRAD,LOGGRD,LOGMSH,
     *                                                   LOGBIS 
      COMMON
     *       /PRINAL/ IMTALK
      COMMON
     *       /OCCREC/ N_ACTU_PROTON,N_CORR_PROTON,
     *                N_ACTU_NEUTRS,N_CORR_NEUTRS
C
      DATA
     *       FITNUC(01) / 'Ox' /,
     *       FITNUC(02) / 'C0' /,
     *       FITNUC(03) / 'C8' /,
     *       FITNUC(04) / 'Ni' /,
     *       FITNUC(05) / 'Zr' /,
     *       FITNUC(06) / 'Sn' /,
     *       FITNUC(07) / 'Gd' /,
     *       FITNUC(08) / 'Pb' /
C
C=======================================================================
C
C     This subroutine does the minimization of the \chi^2 using the 
C     pseudo-experimental levels. Using the Monte-Carlo techniques, 
C     we give to these levels and error, normally distributed with 
C     the \sigma of our choice.
C
C     The mean and the sigma values of the gaussian over which the 
C     random numbers are generated are read from the input file:
C
C                          ws15_pseudo-exper.d
C
C     They depend on the nucleus, particle and energy level.
C____________________
C
C     BIG WARNING:
C
C     Since we did not want to complicate the minimization subroutines
C     of the rest of our code, we redifine the variables EXPEXP with
C     the pseudo-experimental levels. Anyway, the real experimental
C     values are kept in someother vectors under the common blocks:
C
C                   /REALEX_PROTON/ and /REALEX_NEUTRS/
C
C=======================================================================
C
C     First, reading the pseudo-experimental levels
C
      IF (LOGWRI.GT.0) WRITE(LOGFIL,'(''Entering READIN_PSEUDO from '',
     *                                ''FITTIN_MONTEC'')')
      CALL READIN_PSEUDO
C
C=======================================================================
C
C     Here the 'dangerous' part: transforming the real experimental
C     results to the pseudo-experimental
C
      IF (LOGWRI.GT.0) WRITE(LOGFIL,'(/,''Entering EX_TO_PSEUDO from '',
     *                                  ''FITTIN_MONTEC'')')
      CALL EX_TO_PSEUDO
C
      LEVACT_PROTON=0
      LEVACT_NEUTRS=0
      DO INUCLI=1,LDNUCL
         IF (ITAKNU(INUCLI).EQ.1) THEN
             LEVACT_PROTON=LEVACT_PROTON+LEVEXP_PROTON(INUCLI)
             LEVACT_NEUTRS=LEVACT_NEUTRS+LEVEXP_NEUTRS(INUCLI)
             SIGACT_PROTON=SIGEXP_PROTON(1,INUCLI)
             SIGACT_NEUTRS=SIGEXP_NEUTRS(1,INUCLI)
         END IF
      END DO
C
C=======================================================================
C
C     Identifying the parameters over we are going to minimize, 
C     creating the logfiles and counting how many contributions
C     we want to have to the \chi^2
C
C=======================================================================
C
      IFUNCT_EVALUS=0
C
      IFIRST=1
      IFPROT=0
      IFNEUT=0
      IFBOTH=0
      IACTIV=0
C
      DO IPARAM=1,NDPARS ! Fixed as NDPARS=48
C
         IF (IFTAKE(IPARAM).EQ.1) THEN
C
             IACTIV=IACTIV+1
             ARGPAR(IACTIV)=VMISTR(IPARAM)
C  
             IF (IPARAM.LE.20 .AND. IFDENS.EQ.0) THEN
                 IFPROT=1
                 LEVACT=LEVACT_PROTON
                 SIGACT=SIGACT_PROTON
             END IF
C
             IF (IPARAM.GT.20 .AND. IFDENS.EQ.0) THEN
                 IFNEUT=1
                 LEVACT=LEVACT_NEUTRS
                 SIGACT=SIGACT_NEUTRS
             END IF
C
             IF (IPARAM.GE.51 .AND. IFDENS.EQ.0) THEN
                 IFPROT=0
                 IFNEUT=0
                 IFBOTH=1
                 LEVACT=LEVACT_NEUTRS+LEVACT_PROTON
                 SIGACT=SIGACT_PROTON !SIGACT_NEUTRS
             END IF
C
         END IF
C
      END DO
C_______________________________________________________________________
C
      IF ((IFDENS.EQ.0).AND.(IFNEUT.EQ.1).AND.(IFPROT.EQ.1)) THEN
C  
           IFPROT=0
           IFNEUT=0
           IFBOTH=1
C
          WRITE(0,'(/,''ALERT in LMMINI!! Density is NOT activated'',1X,
     *                ''and IFNEUT and IFPROT are equal to 1 !!!!'',/)')
C
C          STOP 'Chose proton OR neutron parameters, stop FROM LMMINI'
C
      END IF
C
C=======================================================================
C
      IF (LOGWRI.GT.0) THEN
C
          WRITE(LOGFIL,'(/,9X,''Number of active parameters '',
     *                        ''IACTIV='',I2)') IACTIV
C
          IF (IFPROT.EQ.1) THEN
              WRITE(LOGFIL,'(9X,''Minimisation for protons'')')
          END IF
C
          IF (IFNEUT.EQ.1) THEN
              WRITE(LOGFIL,'(9X,''Minimisation for neutrons'')')
          END IF
C
          IF (IFBOTH.EQ.1) THEN
              WRITE(LOGFIL,'(9X,''Minimisation for protons '',
     *                          ''and neutrons'')')
          END IF
C
      END IF
C
C=======================================================================
C
      IF (IFDENS.EQ.0 .AND. IFPROT.EQ.1) FILNAM_ENDING='P'
      IF (IFDENS.EQ.0 .AND. IFNEUT.EQ.1) FILNAM_ENDING='N'
      IF (IFDENS.EQ.0 .AND. IFBOTH.EQ.1) FILNAM_ENDING='B'
C
      IF (IFDENS.EQ.1) THEN
C
          LEVACT=LEVACT_NEUTRS+LEVACT_PROTON
          SIGACT=SIGACT_PROTON !SIGACT_NEUTRS
C              
          FILNAM_ENDING='B'
C              
          WRITE(TEXLAM,'(''_LAMBDA'',SP,I2,I2,I2,I2)')
     *                               IFPAR1,IFPAR2,IFPAR3,IFPAR4
C          
      END IF
C
C=======================================================================
C 
      KACTIV=0
      DO JNUCLI=1,LDNUCL
         IF (ITAKNU(JNUCLI).EQ.1) THEN
             KACTIV=KACTIV+1
             FITNUC_AUXILI(KACTIV)=FITNUC(JNUCLI)
         END IF
      END DO
      
      I1=3*KACTIV-1
      
      WRITE(FILNAM_NUCFIT,'(<NUCACT>(A2,''-''))')
     *                     (FITNUC_AUXILI(I),I=1,NUCACT)
C
C=======================================================================
C     
      IF (IFDENS.EQ.0) THEN
C      
          WRITE(FILNAM,'(''LogMonte/IFDENS-0/'',A,''_'',A1,
     *                   ''_IFDENS-'',I1,''_IFTENS-'',I1,
     *                   ''_IF-RAD-'',I1,''_IF-ORD-'',I1,
     *                   ''_IF-RHO-'',I1,
     *                   ''_IFDEEP-'',I1,''_IFPRON-'',I1,
     *                   ''_IFK-VC-'',I1,''_IFK-VS-'',I1,
     *                   ''_IFK-RC-'',I1,''_IFK-RS-'',I1,
     *                   ''_IFK-AC-'',I1,''_IFK-AS-'',I1,
     *                   ''_LEVACT-'',I2.2,''_SIGMAS-'',F5.3,
     *                   ''_LDMONT-'',I5.5,''_'',A3)')
     *
     *                      FILNAM_NUCFIT(1:I1),FILNAM_ENDING,
     *                      IFDENS,IFTENS,IF_RAD,IF_INV,IF_RHO,
     *                      IFDEEP,IFPRON,IFK_VC,IFK_VS,IFK_RC,
     *                      IFK_RS,IFK_AC,IFK_AS,LEVACT,SIGACT,
     *                      LDMONT,VERSIO
C          
          IF (NUCACT.GT.1 .AND. NUCACT.LT.LDNUCL) THEN
C              
              WRITE(FILNA2,'(''LogMonte/IFDENS-0/'',A,''_'',
     *                                      A1,''_Energies'',
     *                       ''_IFDENS-'',I1,''_IFTENS-'',I1,
     *                       ''_IF-RAD-'',I1,''_IF-ORD-'',I1,
     *                       ''_IF-RHO-'',I1,
     *                       ''_IFDEEP-'',I1,''_IFPRON-'',I1,
     *                       ''_IFK-VC-'',I1,''_IFK-VS-'',I1,
     *                       ''_IFK-RC-'',I1,''_IFK-RS-'',I1,
     *                       ''_IFK-AC-'',I1,''_IFK-AS-'',I1,
     *                       ''_LEVACT-'',I2.2,''_SIGMAS-'',F5.3,
     *                       ''_LDMONT-'',I5.5,''_'',A3)')
     * 
     *                          FILNAM_NUCFIT(1:I1),FILNAM_ENDING,
     *                          IFDENS,IFTENS,IF_RAD,IF_INV,IF_RHO,
     *                          IFDEEP,IFPRON,IFK_VC,IFK_VS,IFK_RC,
     *                          IFK_RS,IFK_AC,IFK_AS,LEVACT,SIGACT,
     *                          LDMONT,VERSIO
             
              WRITE(FILNA3,'(''LogMonte/IFDENS-0/'',A,''_'',
     *                                         A1,''_Radii'',
     *                       ''_IFDENS-'',I1,''_IFTENS-'',I1,
     *                       ''_IF-RAD-'',I1,''_IF-ORD-'',I1,
     *                       ''_IF-RHO-'',I1,
     *                       ''_IFDEEP-'',I1,''_IFPRON-'',I1,
     *                       ''_IFK-VC-'',I1,''_IFK-VS-'',I1,
     *                       ''_IFK-RC-'',I1,''_IFK-RS-'',I1,
     *                       ''_IFK-AC-'',I1,''_IFK-AS-'',I1,
     *                       ''_LEVACT-'',I2.2,''_SIGMAS-'',F5.3,
     *                       ''_LDMONT-'',I5.5,''_'',A3)')
     * 
     *                          FILNAM_NUCFIT(1:I1),FILNAM_ENDING,
     *                          IFDENS,IFTENS,IF_RAD,IF_INV,IF_RHO,
     *                          IFDEEP,IFPRON,IFK_VC,IFK_VS,IFK_RC,
     *                          IFK_RS,IFK_AC,IFK_AS,LEVACT,SIGACT,
     *                          LDMONT,VERSIO
C
              DO K=1,256
                 IF (FILNA2(K:K).EQ.'.') THEN
                     FILNA2(K:K)='_'
                 END IF
                 IF (FILNA2(K:K).EQ.' ') THEN
                     KLENG1=K-1
                     GO TO 1
                 END IF
              END DO
C
   1          CONTINUE
C
              FILNA2=FILNA2(1:KLENG1)//'.log'
C
              DO K=1,256
                 IF (FILNA3(K:K).EQ.'.') THEN
                     FILNA3(K:K)='_'
                 END IF
                 IF (FILNA3(K:K).EQ.' ') THEN
                     KLENG2=K-1
                     GO TO 2
                 END IF
              END DO
C
   2          CONTINUE
C
              FILNA3=FILNA3(1:KLENG2)//'.log'
C
              OPEN(LOGENE,FILE=FILNA2,STATUS='UNKNOWN',FORM='FORMATTED')
              OPEN(LOGRAD,FILE=FILNA3,STATUS='UNKNOWN',FORM='FORMATTED')
C
          END IF
C      
      END IF
C
C=======================================================================
C     
      IF (IFDENS.EQ.0) THEN
C      
          WRITE(FILNAM,'(''LogMonte/IFDENS-0/'',A,''_'',A1,
     *                   ''_IFDENS-'',I1,''_IFTENS-'',I1,
     *                   ''_IF-PAI-'',I1,
     *                   ''_IF-RAD-'',I1,''_IF-ORD-'',I1,
     *                   ''_IF-RHO-'',I1,
     *                   ''_IFDEEP-'',I1,''_IFPRON-'',I1,
     *                   ''_IFK-VC-'',I1,''_IFK-VS-'',I1,
     *                   ''_IFK-RC-'',I1,''_IFK-RS-'',I1,
     *                   ''_IFK-AC-'',I1,''_IFK-AS-'',I1,
     *                   ''_LEVACT-'',I2.2,''_SIGMAS-'',F5.3,
     *                   ''_LDMONT-'',I5.5,''_'',A3)')
     *
     *                          FILNAM_NUCFIT(1:I1),FILNAM_ENDING,
     *                          IFDENS,IFTENS,IF_PAI,IF_RAD,IF_INV,
     *                          IF_RHO,IFDEEP,IFPRON,IFK_VC,IFK_VS,
     *                          IFK_RC,IFK_RS,IFK_AC,IFK_AS,
     *                          LEVACT,SIGACT,LDMONT,VERSIO
C          
          IF (NUCACT.GT.1 .AND. NUCACT.LT.LDNUCL) THEN
C              
              WRITE(FILNA2,'(''LogMonte/IFDENS-0/'',A,''_'',
     *                                      A1,''_Energies'',
     *                       ''_IFDENS-'',I1,''_IFTENS-'',I1,
     *                       ''_IF-PAI-'',I1,
     *                       ''_IF-RAD-'',I1,''_IF-ORD-'',I1,
     *                       ''_IF-RHO-'',I1,
     *                       ''_IFDEEP-'',I1,''_IFPRON-'',I1,
     *                       ''_IFK-VC-'',I1,''_IFK-VS-'',I1,
     *                       ''_IFK-RC-'',I1,''_IFK-RS-'',I1,
     *                       ''_IFK-AC-'',I1,''_IFK-AS-'',I1,
     *                       ''_LEVACT-'',I2.2,''_SIGMAS-'',F5.3,
     *                       ''_LDMONT-'',I5.5,''_'',A3)')
     * 
     *                          FILNAM_NUCFIT(1:I1),FILNAM_ENDING,
     *                          IFDENS,IFTENS,IF_PAI,IF_RAD,IF_INV,
     *                          IF_RHO,IFDEEP,IFPRON,IFK_VC,IFK_VS,
     *                          IFK_RC,IFK_RS,IFK_AC,IFK_AS,
     *                          LEVACT,SIGACT,LDMONT,VERSIO
             
              WRITE(FILNA3,'(''LogMonte/IFDENS-0/'',A,''_'',
     *                                         A1,''_Radii'',
     *                       ''_IFDENS-'',I1,''_IFTENS-'',I1,
     *                       ''_IF-PAI-'',I1,
     *                       ''_IF-RAD-'',I1,''_IF-ORD-'',I1,
     *                       ''_IF-RHO-'',I1,
     *                       ''_IFDEEP-'',I1,''_IFPRON-'',I1,
     *                       ''_IFK-VC-'',I1,''_IFK-VS-'',I1,
     *                       ''_IFK-RC-'',I1,''_IFK-RS-'',I1,
     *                       ''_IFK-AC-'',I1,''_IFK-AS-'',I1,
     *                       ''_LEVACT-'',I2.2,''_SIGMAS-'',F5.3,
     *                       ''_LDMONT-'',I5.5,''_'',A3)')
     * 
     *                          FILNAM_NUCFIT(1:I1),FILNAM_ENDING,
     *                          IFDENS,IFTENS,IF_PAI,IF_RAD,IF_INV,
     *                          IF_RHO,IFDEEP,IFPRON,IFK_VC,IFK_VS,
     *                          IFK_RC,IFK_RS,IFK_AC,IFK_AS,
     *                          LEVACT,SIGACT,LDMONT,VERSIO
C
              DO K=1,256
                 IF (FILNA2(K:K).EQ.'.') THEN
                     FILNA2(K:K)='_'
                 END IF
                 IF (FILNA2(K:K).EQ.' ') THEN
                     KLENG1=K-1
                     GO TO 3
                 END IF
              END DO
C
   3          CONTINUE
C
              FILNA2=FILNA2(1:KLENG1)//'.log'
C
              DO K=1,256
                 IF (FILNA3(K:K).EQ.'.') THEN
                     FILNA3(K:K)='_'
                 END IF
                 IF (FILNA3(K:K).EQ.' ') THEN
                     KLENG2=K-1
                     GO TO 4
                 END IF
              END DO
C
   4          CONTINUE
C
              FILNA3=FILNA3(1:KLENG2)//'.log'
C
              OPEN(LOGENE,FILE=FILNA2,STATUS='UNKNOWN',FORM='FORMATTED')
              OPEN(LOGRAD,FILE=FILNA3,STATUS='UNKNOWN',FORM='FORMATTED')
C
          END IF
C
      END IF
C
C=======================================================================
C
      IF (IFDENS.EQ.1) THEN
C
          WRITE(FILNAM,'(''LogMonte/IFDENS-1_IFTENS-'',I1,''/'',
     *                                        A,''_'',A1,
     *                   ''_IFDENS-'',I1,''_IFTENS-'',I1,
     *                   ''_IF-PAI-'',I1,
     *                   ''_IF-RAD-'',I1,''_IF-ORD-'',I1,
     *                   ''_IFDEEP-'',I1,''_IFPRON-'',I1,
     *                   ''_IFK-VC-'',I1,''_IFK-VS-'',I1,
     *                   ''_IFK-RC-'',I1,''_IFK-RS-'',I1,
     *                   ''_IFK-AC-'',I1,''_IFK-AS-'',I1,
     *                   ''_LEVACp-'',I2.2,''_LEVACn-'',I2.2,
     *                                     ''_SIGMAS-'',F5.3,
     *                   ''_LDMONT-'',I5.5,A15,''_'',A3)')
     *
     *                      IFTENS,FILNAM_NUCFIT(1:I1),FILNAM_ENDING,
     *                      IFDENS,IFTENS,IF_PAI,IF_RAD,IF_INV,
     *                      IFDEEP,IFPRON,IFK_VC,IFK_VS,
     *                      IFK_RC,IFK_RS,IFK_AC,IFK_AS,
     *                      LEVACT_PROTON,LEVACT_NEUTRS,SIGACT,
     *                                    LDMONT,TEXLAM,VERSIO
C______________________________________________________________________
C
          IF (IFTENS.EQ.1) THEN ! IFDENS=1
C
              WRITE(FILNAM,'(''LogMonte/IFDENS-1_IFTENS-'',I1,''/'',
     *                                        A,''_'',A1,
     *                   ''_IFDENS-'',I1,''_IFTENS-'',I1,
     *                   ''_IF-PAI-'',I1,
     *                   ''_ICENTT-'',I1,''_ISORBT-'',I1,
     *                                   ''_ITENSR-'',I1,
     *                   ''_IF-RAD-'',I1,''_IF-ORD-'',I1,
     *                   ''_IFDEEP-'',I1,''_IFPRON-'',I1,
     *                   ''_IFK-VC-'',I1,''_IFK-VS-'',I1,
     *                   ''_IFK-RC-'',I1,''_IFK-RS-'',I1,
     *                   ''_IFK-AC-'',I1,''_IFK-AS-'',I1,
     *                   ''_LEVACp-'',I2.2,''_LEVACn-'',I2.2,
     *                                     ''_SIGMAS-'',F5.3,
     *                   ''_LDMONT-'',I5.5,A15,''_'',A3)')
     *
     *                      IFTENS,FILNAM_NUCFIT(1:I1),FILNAM_ENDING,
     *                      IFDENS,IFTENS,IF_PAI,ICENTT,ISORBT,ITENSR,
     *                      IF_RAD,IF_INV,IFDEEP,IFPRON,IFK_VC,IFK_VS,
     *                      IFK_RC,IFK_RS,IFK_AC,IFK_AS,
     *                      LEVACT_PROTON,LEVACT_NEUTRS,SIGACT,
     *                                    LDMONT,TEXLAM,VERSIO
          END IF
C______________________________________________________________________
C
          IF (NUCACT.GT.1 .AND. NUCACT.LT.LDNUCL) THEN ! IFDENS=1
C
              WRITE(FILNA2,'(''LogMonte/IFDENS-1_IFTENS-'',I1,''/'',
     *                              A,''_'',A1,''_Energies'',
     *                       ''_IFDENS-'',I1,''_IFTENS-'',I1,
     *                       ''_IF-PAI-'',I1,
     *                       ''_IF-RAD-'',I1,''_IF-ORD-'',I1,
     *                       ''_IFDEEP-'',I1,''_IFPRON-'',I1,
     *                       ''_IFK-VC-'',I1,''_IFK-VS-'',I1,
     *                       ''_IFK-RC-'',I1,''_IFK-RS-'',I1,
     *                       ''_IFK-AC-'',I1,''_IFK-AS-'',I1,
     *                       ''_LEVACp-'',I2.2,''_LEVACn-'',I2.2,
     *                                     ''_SIGMAS-'',F5.3,
     *                       ''_LDMONT-'',I5.5,A15,''_'',A3)')
     *
     *                      IFTENS,FILNAM_NUCFIT(1:I1),FILNAM_ENDING,
     *                      IFDENS,IFTENS,IF_PAI,
     *                      IF_RAD,IF_INV,IFDEEP,IFPRON,IFK_VC,IFK_VS,
     *                      IFK_RC,IFK_RS,IFK_AC,IFK_AS,
     *                      LEVACT_PROTON,LEVACT_NEUTRS,SIGACT,
     *                                    LDMONT,TEXLAM,VERSIO
C______________________________________________________________________
C
              WRITE(FILNA3,'(''LogMonte/IFDENS-1_IFTENS-'',I1,''/'',
     *                                 A,''_'',A1,''_Radii'',
     *                       ''_IFDENS-'',I1,''_IFTENS-'',I1,
     *                       ''_IF-PAI-'',I1,
     *                       ''_IF-RAD-'',I1,''_IF-ORD-'',I1,
     *                       ''_IFDEEP-'',I1,''_IFPRON-'',I1,
     *                       ''_IFK-VC-'',I1,''_IFK-VS-'',I1,
     *                       ''_IFK-RC-'',I1,''_IFK-RS-'',I1,
     *                       ''_IFK-AC-'',I1,''_IFK-AS-'',I1,
     *                       ''_LEVACp-'',I2.2,''_LEVACn-'',I2.2,
     *                                     ''_SIGMAS-'',F5.3,
     *                       ''_LDMONT-'',I5.5,A15,''_'',A3)')
     * 
     *                      IFTENS,FILNAM_NUCFIT(1:I1),FILNAM_ENDING,
     *                      IFDENS,IFTENS,IF_PAI,
     *                      IF_RAD,IF_INV,IFDEEP,IFPRON,IFK_VC,IFK_VS,
     *                      IFK_RC,IFK_RS,IFK_AC,IFK_AS,
     *                      LEVACT_PROTON,LEVACT_NEUTRS,SIGACT,
     *                                    LDMONT,TEXLAM,VERSIO
C______________________________________________________________________
C
              IF (IFTENS.EQ.1) THEN ! IFDENS=1
C
                  WRITE(FILNA2,'(''LogMonte/IFDENS-1_IFTENS-'',I1,''/'',
     *                              A,''_'',A1,''_Energies'',
     *                           ''_IFDENS-'',I1,''_IFTENS-'',I1,
     *                           ''_IF-PAI-'',I1,
     *                           ''_ICENTT-'',I1,''_ISORBT-'',I1,
     *                                           ''_ITENSR-'',I1,
     *                           ''_IF-RAD-'',I1,''_IF-ORD-'',I1,
     *                           ''_IFDEEP-'',I1,''_IFPRON-'',I1,
     *                           ''_IFK-VC-'',I1,''_IFK-VS-'',I1,
     *                           ''_IFK-RC-'',I1,''_IFK-RS-'',I1,
     *                           ''_IFK-AC-'',I1,''_IFK-AS-'',I1,
     *                           ''_LEVACp-'',I2.2,''_LEVACn-'',I2.2,
     *                                     ''_SIGMAS-'',F5.3,
     *                           ''_LDMONT-'',I5.5,A15,''_'',A3)')
     *
     *                      IFTENS,FILNAM_NUCFIT(1:I1),FILNAM_ENDING,
     *                      IFDENS,IFTENS,IF_PAI,ICENTT,ISORBT,ITENSR,
     *                      IF_RAD,IF_INV,IFDEEP,IFPRON,IFK_VC,IFK_VS,
     *                      IFK_RC,IFK_RS,IFK_AC,IFK_AS,
     *                      LEVACT_PROTON,LEVACT_NEUTRS,SIGACT,
     *                                    LDMONT,TEXLAM,VERSIO
C
C______________________________________________________________________
C
                  WRITE(FILNA3,'(''LogMonte/IFDENS-1_IFTENS-'',I1,''/'',
     *                                 A,''_'',A1,''_Radii'',
     *                           ''_IFDENS-'',I1,''_IFTENS-'',I1,
     *                           ''_IF-PAI-'',I1,
     *                           ''_ICENTT-'',I1,''_ISORBT-'',I1,
     *                                           ''_ITENSR-'',I1,
     *                           ''_IF-RAD-'',I1,''_IF-ORD-'',I1,
     *                           ''_IFDEEP-'',I1,''_IFPRON-'',I1,
     *                           ''_IFK-VC-'',I1,''_IFK-VS-'',I1,
     *                           ''_IFK-RC-'',I1,''_IFK-RS-'',I1,
     *                           ''_IFK-AC-'',I1,''_IFK-AS-'',I1,
     *                           ''_LEVACp-'',I2.2,''_LEVACn-'',I2.2,
     *                                     ''_SIGMAS-'',F5.3,
     *                           ''_LDMONT-'',I5.5,A15,''_'',A3)')
     *
     *                      IFTENS,FILNAM_NUCFIT(1:I1),FILNAM_ENDING,
     *                      IFDENS,IFTENS,IF_PAI,ICENTT,ISORBT,ITENSR,
     *                      IF_RAD,IF_INV,IFDEEP,IFPRON,IFK_VC,IFK_VS,
     *                      IFK_RC,IFK_RS,IFK_AC,IFK_AS,
     *                      LEVACT_PROTON,LEVACT_NEUTRS,SIGACT,
     *                                    LDMONT,TEXLAM,VERSIO
C
              END IF
C
              DO K=1,256
                 IF (FILNA2(K:K).EQ.'.') THEN
                     FILNA2(K:K)='_'
                 END IF
                 IF (FILNA2(K:K).EQ.' ') THEN
                     KLENG1=K-1
                     GO TO 5
                 END IF
              END DO
C
   5          CONTINUE
C
              FILNA2=FILNA2(1:KLENG1)//'.log'
C
              DO K=1,256
                 IF (FILNA3(K:K).EQ.'.') THEN
                     FILNA3(K:K)='_'
                 END IF
                 IF (FILNA3(K:K).EQ.' ') THEN
                     KLENG2=K-1
                     GO TO 6
                 END IF
              END DO
C
   6          CONTINUE
C
              FILNA3=FILNA3(1:KLENG2)//'.log'
C
              OPEN(LOGENE,FILE=FILNA2,STATUS='UNKNOWN',FORM='FORMATTED')
              OPEN(LOGRAD,FILE=FILNA3,STATUS='UNKNOWN',FORM='FORMATTED')
C
          END IF
C
      END IF
C
C=======================================================================
C
      DO K=1,256
         IF (FILNAM(K:K).EQ.'.') THEN
             FILNAM(K:K)='_'
         END IF
         IF (FILNAM(K:K).EQ.' ') THEN
             KLENG3=K-1
             GO TO 7
         END IF
      END DO
   7  CONTINUE
C
      FILNAM=FILNAM(1:KLENG3)//'.log'
C
C=======================================================================
C
      OPEN(LOGAUX,FILE=FILNAM,STATUS='UNKNOWN',FORM='FORMATTED')
C
C=======================================================================
C
C     Defining the dimension of the vector-function to be minimised
C
      LEVNUM=0
C
      IF (IFDENS.EQ.1 .OR. (IFDENS.EQ.0 .AND. IFBOTH.EQ.1)) THEN
C
          DO INUCLI=1,LDNUCL
C
             IF (ITAKNU(INUCLI).EQ.1) THEN
C
                 IF (IF_SPE.EQ.1) LEVNUM=LEVNUM+LEVEXP_PROTON(INUCLI)
     *                                         +LEVEXP_NEUTRS(INUCLI)
     
                 IF (IF_RAD.EQ.1) LEVNUM=LEVNUM+2 ! Two radii for each 
C                                                       nucleus, (P+N)
C                                                              
                 IF (IF_GAP.EQ.1) LEVNUM=LEVNUM+2 ! Two Gap Energies 
C                                                   for each nucleus
C                                                                (P+N)
C                                                     
                 IF (IF_FER.EQ.1) LEVNUM=LEVNUM+2 ! Two Fermi Energies 
C                                                   for  each  nucleus
C                                                                (P+N)
C
                 IF (IF_DEN.EQ.1) LEVNUM=LEVNUM+4 ! Four densities for 
C                                                   each nucleus: 
C                                                   below and above the 
C                                                   shell closure (P+N)
                 IF (IF_RHO.EQ.1) THEN
C             
c                     NRHOEX=0
C             
c                     CALL COUNTI_RHOEXP(INUCLI,NRHOEX)
C             
c                     LEVNUM=LEVNUM+NRHOEX
C             
                 END IF
C            
             END IF
C        
          END DO ! of LDNUCL
C       
       END IF ! IFDENS=1
C
C=======================================================================
C       
       IF ((IFDENS.EQ.0).AND.(IFPROT.EQ.1)) THEN 
C           
           DO INUCLI=1,LDNUCL
C        
              IF (ITAKNU(INUCLI).EQ.1) THEN
C         
                  IF (IF_SPE.EQ.1) LEVNUM=LEVNUM+LEVEXP_PROTON(INUCLI)
C     
                  IF (IF_RAD.EQ.1) LEVNUM=LEVNUM+1 ! Proton radius
C                                                              
                  IF (IF_GAP.EQ.1) LEVNUM=LEVNUM+1 ! Proton Gap energy
C                                                     
                  IF (IF_FER.EQ.1) LEVNUM=LEVNUM+1 ! Proton Fermi energy
C        
                  IF (IF_DEN.EQ.1) LEVNUM=LEVNUM+2 ! Proton up and down 
C                                                    density energies
C
                  IF (IF_RHO.EQ.1) THEN
*C             
*                      NRHOEX=0
*C             
*                      CALL COUNTI_RHOEXP(INUCLI,NRHOEX)
*C             
*                      LEVNUM=LEVNUM+NRHOEX
C             
                  END IF
C            
              END IF
C        
           END DO
C       
       END IF
C
C=======================================================================
C       
       IF ((IFDENS.EQ.0).AND.(IFNEUT.EQ.1)) THEN 
C
           DO INUCLI=1,LDNUCL
C
              IF (ITAKNU(INUCLI).EQ.1) THEN
C
                  IF (IF_SPE.EQ.1) LEVNUM=LEVNUM+LEVEXP_NEUTRS(INUCLI)
     
                  IF (IF_RAD.EQ.1) LEVNUM=LEVNUM+1 ! Neutron radius
C                                                              
                  IF (IF_GAP.EQ.1) LEVNUM=LEVNUM+1 ! Neutron Gap energy
C                                                     
                  IF (IF_FER.EQ.1) LEVNUM=LEVNUM+1 ! Neutron Fermi 
C                                                                energy
C
                  IF (IF_DEN.EQ.1) LEVNUM=LEVNUM+2 ! Neutron up and 
C                                                    down density 
C                                                              energies
C                 IF (IF_RHO.EQ.1)
C                 IF (IF_INV.EQ.1)
            
              END IF
         
           END DO
       
       END IF       
C
C=======================================================================
C
       IF (LEVNUM.GT.NDIM_M) THEN
C
           WRITE(LSCREN,'(/,''LEVNUM= '',I4,''and NDIM_M= '',I4,
     *                  ''--> LEVNUM should be LT NDIM_M!!'',/)')
C
           STOP 'STOP in LMMINI: LEVNUM.GT.NDIM_M!!'
C
       END IF
C
C=======================================================================
C 
      IF (LOGWRI.GT.0) THEN
C 
          WRITE(LOGFIL,'(9X,''LEVNUM='',I2,'' including: '')') LEVNUM
C
          IF (IF_SPE.EQ.1) WRITE(LOGFIL,'(30X,''Single Particle'',
     *                                        '' Energies'')')
C
          IF (IF_RAD.EQ.1) WRITE(LOGFIL,'(30X,''RMS Radii'')')
C
          IF (IF_GAP.EQ.1) WRITE(LOGFIL,'(30X,''Gap Energy'')')   
C       
          IF (IF_FER.EQ.1) WRITE(LOGFIL,'(30X,''Fermi Energy'')')  
       
          IF (IF_DEN.EQ.1) WRITE(LOGFIL,'(''Densities below and '',
     *                                    ''above the shell closure'')') 
C         IF (IF_RHO.EQ.1)
C         IF (IF_INV.EQ.1)      
C 
      END IF
C
C=======================================================================
C
C     To avoid huge log/output files
C
      LOGWRI=0
      IMTALK=0
C
      WRITE(LOGFIL,'(/,''Setting LOGWRI=0 and IMTALK=0 to avoid '',
     *                 ''huge log/output files'')')
C
C=======================================================================
C
C     Opening the temporary files where we will store the results to
C     create the histograms. This is done to avoid the constructions
C     of huge matrices.
C
C
C     Opening files for the parameters
C
      NWRITE_PARAMS=1000
C
      DO IPARAM=1,NDPARS
         IF (IFTAKE(IPARAM).EQ.1) THEN
C
             NWRITE_PARAMS=NWRITE_PARAMS+1
             NRESUL=NWRITE_PARAMS
C
             PARNAM=TITLES(IPARAM)(3:10)
             CALL HISTOG_PARNAM(PARNAM,LEVACT,LEVACT_PROTON,
     *                          LEVACT_NEUTRS,SIGACT,FILPAR)
C
             WRITE(0,'(''Opening file: '',A)') FILPAR
C
             OPEN(UNIT=NWRITE_PARAMS,FILE=FILPAR,STATUS='UNKNOWN')
C
             LDCURV=1
C
             IF (IPARAM.LE.20) ISOSPI=1
             IF (IPARAM.GT.20) ISOSPI=0
             IF (IPARAM.GT.38) ISOSPI=2 !both particles
C
             YEXPEC=VMISTR(IPARAM)
C
             LENGT1=0
             DO K=50,1,-1
                IF (TITLES_LATEXS(IPARAM)(K:K).NE.' ') THEN
                    LENGT1=K
                    GO TO 8
                END IF
             END DO
   8         CONTINUE
C
             LENGT2=0
             DO K=40,1,-1
                IF (PARPOT_UNITSS(IPARAM)(K:K).NE.' ') THEN
                    LENGT2=K
                    GO TO 9
                END IF
             END DO
   9         CONTINUE
C
             HISTIT='Parameter Distribution'
             HX_TIT='\\boldmath'//TITLES_LATEXS(IPARAM)(1:LENGT1)//
     *                      ' ('//PARPOT_UNITSS(IPARAM)(1:LENGT2)//')'
             HY_TIT='Probability Distribution'
             XUNITS=PARPOT_UNITSS(IPARAM)
C
             CALL PLOTINFOPRINT(NRESUL,HISTIT,HX_TIT,HY_TIT,
     *                          LEVACT_PROTON,LEVACT_NEUTRS,
     *                          SIGACT_PROTON,SIGACT_NEUTRS,
     *                                 LDCURV,YEXPEC,XUNITS)
C
         END IF
      END DO
C
C     Opening files for rms
C
      IPRINT_PROTON=0
      IPRINT_NEUTRS=0
C
      IF ((IFDENS.EQ.0 .AND. IFPROT.EQ.1) .OR.
     *    (IFDENS.EQ.0 .AND. IFBOTH.EQ.1) .OR.
     *     IFDENS.EQ.1) THEN
C
          NWRITE_RMSPRO=2001
          NRESUL=NWRITE_RMSPRO
          ISOSPI=1
          CALL HISTOG_RMSNAM(ISOSPI,LEVACT,LEVACT_PROTON,
     *                       LEVACT_NEUTRS,SIGACT,FILRMS)
C
          WRITE(0,'(''Opening file: '',A)') FILRMS
C
          OPEN(UNIT=NWRITE_RMSPRO,FILE=FILRMS,STATUS='UNKNOWN')
C
          LDCURV=1
C
          YEXPEC=0.0
C
          HISTIT='\\boldmath$\\chi$ Distribution'
          HX_TIT='\\boldmath${\\rm rms}_p$ (MeV)'
          HY_TIT='Probability Distribution'
          XUNITS='MeV'
C
          CALL PLOTINFOPRINT(NRESUL,HISTIT,HX_TIT,HY_TIT,
     *                       LEVACT_PROTON,LEVACT_NEUTRS,
     *                       SIGACT_PROTON,SIGACT_NEUTRS,
     *                              LDCURV,YEXPEC,XUNITS)
C
          IPRINT_PROTON=1
C
      END IF
C
      IF ((IFDENS.EQ.0 .AND. IFNEUT.EQ.1) .OR.
     *    (IFDENS.EQ.0 .AND. IFBOTH.EQ.1) .OR.
     *     IFDENS.EQ.1) THEN
C
          NWRITE_RMSNEU=2000
          NRESUL=NWRITE_RMSNEU
          ISOSPI=0
          CALL HISTOG_RMSNAM(ISOSPI,LEVACT,LEVACT_PROTON,
     *                       LEVACT_NEUTRS,SIGACT,FILRMS)
C
          WRITE(0,'(''Opening file: '',A)') FILRMS
C
          OPEN(UNIT=NWRITE_RMSNEU,FILE=FILRMS,STATUS='UNKNOWN')
C
          LDCURV=1
C
          YEXPEC=0.0
C
          HISTIT='\\boldmath$\\chi$ Distribution'
          HX_TIT='\\boldmath${\\rm rms}_p$ (MeV)'
          HY_TIT='Probability Distribution'
          XUNITS='MeV'
C
          CALL PLOTINFOPRINT(NRESUL,HISTIT,HX_TIT,HY_TIT,
     *                       LEVACT_PROTON,LEVACT_NEUTRS,
     *                       SIGACT_PROTON,SIGACT_NEUTRS,
     *                              LDCURV,YEXPEC,XUNITS)
C
          IPRINT_NEUTRS=1
C
      END IF
C
C     Opening files for the spe
C
      IF ((IFDENS.EQ.0 .AND. IFPROT.EQ.1) .OR.
     *    (IFDENS.EQ.0 .AND. IFBOTH.EQ.1) .OR.
     *     IFDENS.EQ.1) THEN
C
          NWRITE_ENEPRO=3001
          NRESUL=NWRITE_ENEPRO
C
          ISOSPI=1
          CALL HISTOG_SPENAM(ISOSPI,LEVACT,LEVACT_PROTON,
     *                       LEVACT_NEUTRS,SIGACT,FILSPE)
C
          WRITE(0,'(''Opening file: '',A)') FILSPE
C
          OPEN(UNIT=NRESUL,FILE=FILSPE,STATUS='UNKNOWN')
C
          LDCURV=(NSHELL_PROTON+1)*(NSHELL_PROTON/2+1)
C
          YEXPEC=99.9999
C
          WRITE(HISTIT,'(''Uncertainty Distributions for '',
     *                   ''\\boldmath$N_p= '',i3,
     *                   ''$'')')LEVACT_PROTON
          HX_TIT='\\boldmath${\\rm Energy (MeV)}$'
          HY_TIT='Proton Histograms'
          XUNITS='MeV'
C
          CALL PLOTINFOPRINT(NRESUL,HISTIT,HX_TIT,HY_TIT,
     *                       LEVACT_PROTON,LEVACT_NEUTRS,
     *                       SIGACT_PROTON,SIGACT_NEUTRS,
     *                              LDCURV,YEXPEC,XUNITS)
C
      END IF
C
      IF ((IFDENS.EQ.0 .AND. IFNEUT.EQ.1) .OR.
     *    (IFDENS.EQ.0 .AND. IFBOTH.EQ.1) .OR.
     *     IFDENS.EQ.1) THEN
C
          NWRITE_ENENEU=3000
          NRESUL=NWRITE_ENENEU
C
          ISOSPI=0
          CALL HISTOG_SPENAM(ISOSPI,LEVACT,LEVACT_PROTON,
     *                       LEVACT_NEUTRS,SIGACT,FILSPE)
C
          OPEN(UNIT=NRESUL,FILE=FILSPE,STATUS='UNKNOWN')
C
          WRITE(0,'(''Opening file: '',A)') FILSPE
C
          LDCURV=(NSHELL_NEUTRS+1)*(NSHELL_NEUTRS/2+1)
C
          YEXPEC=99.9999
C
          WRITE(HISTIT,'(''Uncertainty Distributions for '',
     *                   ''\\boldmath$N_n= '',i3,
     *                   ''$'')')LEVACT_NEUTRS
          HX_TIT='\\boldmath${\\rm Energy (MeV)}$'
          HY_TIT='Neutron Histograms'
          XUNITS='MeV'
C
          CALL PLOTINFOPRINT(NRESUL,HISTIT,HX_TIT,HY_TIT,
     *                       LEVACT_PROTON,LEVACT_NEUTRS,
     *                       SIGACT_PROTON,SIGACT_NEUTRS,
     *                              LDCURV,YEXPEC,XUNITS)
C
      END IF
C
C=======================================================================
C
C     Initializing the maxima and the minima for the parameters,
C     rms and s.p.e.
C
      DO IPARAM=1,NDPARS
         PARPOT_MAXIMA(IPARAM)=-1.0E+10
         PARPOT_MINIMA(IPARAM)=+1.0E+10
      END DO
C
      RMSVAL_PROMAX=-1.0E+10
      RMSVAL_PROMIN=+1.0E+10
C
      RMSVAL_NEUMAX=-1.0E+10
      RMSVAL_NEUMIN=+1.0E+10
C
      ENETHE_PROMAX=-1.0E+10
      ENETHE_PROMIN=+1.0E+10
C
      ENETHE_NEUMAX=-1.0E+10
      ENETHE_NEUMIN=+1.0E+10
C
C=======================================================================
C
C     Beginning the real calculations: generating the random numbers
C     for the noise for the psuedo-experimental levels and minimizing
C     the \chi^2
C
      INSEED=0!7893
      CALL ZBQLINI(INSEED)
C
      IERROR_PROTON=0  !To count how many times we did not get the
      IERROR_NEUTRS=0  !corresponding particle number
C
      DO IMONTE=1,LDMONT
C_______________________________________________________________________
C
C        Randomly generating the noise for the pseudo-levels
C
         DO INUCLI=1,LDNUCL
            IF (ITAKNU(INUCLI).EQ.1) THEN
C
                IEXPER=0
                DO I=1,LEVPRO_PSEUDO(INUCLI)          !Protons
C
                   IF (LEVTAK_PROTON(I,INUCLI).EQ.1) THEN
C
                       IEXPER=IEXPER+1
C
                       GAUSIG=SIGPRO_PSEUDO(I,INUCLI)
                       GAUSMU=ENEPRO_PSEUDO(I,INUCLI)
C
                       RANDOM=ZBQLNOR(GAUSMU,GAUSIG)
                       EXPEXP_PROTON(INUCLI,IEXPER)=RANDOM
C
                    END IF
C
                END DO
C
                IEXPER=0
                DO I=1,LEVNEU_PSEUDO(INUCLI)         !Neutrons
C
                   IF (LEVTAK_NEUTRS(I,INUCLI).EQ.1) THEN
C
                       IEXPER=IEXPER+1
C
                       GAUSIG=SIGNEU_PSEUDO(I,INUCLI)
                       GAUSMU=ENENEU_PSEUDO(I,INUCLI)
C
                       RANDOM=ZBQLNOR(GAUSMU,GAUSIG)
                       EXPEXP_NEUTRS(INUCLI,IEXPER)=RANDOM
C
                   END IF
C
                END DO
C
            END IF
         END DO
C_______________________________________________________________________
C
C        Calling the minimization subroutine
C
 100     CONTINUE
C
         CALL LMMINI_MONTEC(IMONTE,LDMONT,IACTIV,ARGPAR,
     *                      LEVNUM,NDLAST,NEWSED,I_SEED)
C_______________________________________________________________________
C
C        Checking if we arrived to the correct particle number.
C        If not, we repeat the iteration -->> GO TO 100
C
         IF (IPRINT_PROTON.EQ.1) THEN
             IF (N_ACTU_PROTON.NE.N_CORR_PROTON) THEN
C
                 WRITE(LOGAUX,'()')
                 WRITE(000000,'()')
                 DO I=1,3
                    WRITE(LOGAUX,'(''PROTON PARTICLE NUMBER NOT '',
     *                             ''SATISFIED - WE REPEAT THIS '',
     *                             ''ITERATION'')')
                    WRITE(000000,'(''PROTON PARTICLE NUMBER NOT '',
     *                             ''SATISFIED - WE REPEAT THIS '',
     *                             ''ITERATION'')')
                 END DO
C
                 IERROR_PROTON=IERROR_PROTON+1
C
                 WRITE(LOGAUX,'()')
                 WRITE(000000,'()')
                 DO I=1,3
                    WRITE(LOGAUX,'(''It already happened for '',
     *                             ''IERROR_PROTON= '',i6,'' times'')')
     *                               IERROR_PROTON
                    WRITE(000000,'(''It already happened for '',
     *                             ''IERROR_PROTON= '',i6,'' times'')')
     *                               IERROR_PROTON
                 END DO
                 WRITE(LOGAUX,'()')
                 WRITE(000000,'()')
C
                 I_SEED=NEWSED
C
                 GO TO 100
C
             END IF
         END IF
C
         IF (IPRINT_NEUTRS.EQ.1) THEN
             IF (N_ACTU_PROTON.NE.N_CORR_PROTON) THEN
C
                 WRITE(LOGAUX,'()')
                 WRITE(000000,'()')
                 DO I=1,3
                    WRITE(LOGAUX,'(/,''NEUTRON PARTICLE NUMBER NOT '',
     *                               ''SATISFIED - WE REPEAT THIS '',
     *                               ''ITERATION'',/)')
                    WRITE(000000,'(/,''NEUTRON PARTICLE NUMBER NOT '',
     *                               ''SATISFIED - WE REPEAT THIS '',
     *                               ''ITERATION'',/)')
                 END DO
C
                 IERROR_NEUTRS=IERROR_NEUTRS+1
C
                 WRITE(LOGAUX,'()')
                 WRITE(000000,'()')
                 DO I=1,3
                    WRITE(LOGAUX,'(''It already happened for '',
     *                             ''IERROR_NEUTRS= '',i6,'' times'')')
     *                               IERROR_NEUTRS
                    WRITE(LOGAUX,'(''It already happened for '',
     *                             ''IERROR_NEUTRS= '',i6,'' times'')')
     *                               IERROR_NEUTRS
                 END DO
                 WRITE(LOGAUX,'()')
                 WRITE(000000,'()')
C
                 I_SEED=NEWSED
C
                 GO TO 100
C
             END IF
         END IF
C_______________________________________________________________________
C
C        Calculating the s.p.e. and rms values after each minimization
C
         DO INUCLI=1,LDNUCL
            IF (ITAKNU(INUCLI).EQ.1) THEN
C
                I_FLAG=1
                I_MODE=0
C
                CALL WS_RUN(INUCLI,CHISQU_CHOICE,IABORT,I_MODE,
     *                      I_FLAG,CHISQU_PROTON,CHISQU_NEUTRS)
C
                LTPROT=LEVTHE_PROTON
                LTNEUT=LEVTHE_NEUTRS
C_______________________________________________________________________
C
C               Printing the rms values at each IMONTE
C
                IF (IPRINT_PROTON.EQ.1) THEN
                    LDCURV=1
                    LABTIT(LDCURV)='${\\rm rms}_p$'
                    VECTOR_RMSVAL(LDCURV)=SQRT(CHIWEI_PROTON)
C
                    CALL PLOTCURVPRINT(NWRITE_RMSPRO,ND_RMS,LDCURV,
     *                                 IMONTE,LABTIT,VECTOR_RMSVAL)
                END IF
C
                IF (IPRINT_NEUTRS.EQ.1) THEN
                    LDCURV=1
                    LABTIT(LDCURV)='${\\rm rms}_n$'
                    VECTOR_RMSVAL(LDCURV)=SQRT(CHIWEI_NEUTRS)
C
                    CALL PLOTCURVPRINT(NWRITE_RMSNEU,ND_RMS,LDCURV,
     *                                 IMONTE,LABTIT,VECTOR_RMSVAL)
                END IF
C
C               MAX - MIN for rms proton
C
                IF (SQRT(CHIWEI_PROTON).LT.RMSVAL_PROMIN) THEN
                    RMSVAL_PROMIN=SQRT(CHIWEI_PROTON)
                END IF
C
                IF (SQRT(CHIWEI_PROTON).GT.RMSVAL_PROMAX) THEN
                    RMSVAL_PROMAX=SQRT(CHIWEI_PROTON)
                END IF
C
C               MAX - MIN for rms neutrons
C
                IF (SQRT(CHIWEI_NEUTRS).LT.RMSVAL_NEUMIN) THEN
                    RMSVAL_NEUMIN=SQRT(CHIWEI_NEUTRS)
                END IF
C
                IF (SQRT(CHIWEI_NEUTRS).GT.RMSVAL_NEUMAX) THEN
                    RMSVAL_NEUMAX=SQRT(CHIWEI_NEUTRS)
                END IF
C_______________________________________________________________________
C
C               Storing the reference label
C
                IF (IMONTE.EQ.1) THEN
C
                    DO IREFER=1,LTPROT
                       LABREF_PROTON(IREFER)=LABTHE_PROTON(IREFER)
                    END DO
C
                    DO IREFER=1,LTNEUT
                       LABREF_NEUTRS(IREFER)=LABTHE_NEUTRS(IREFER)
                    END DO
C
                END IF
C
C               Printing the results for each IMONTE (=restart)
C
                IF (IPRINT_PROTON.EQ.1) THEN
C
                    LDCURV=LTPROT
C
                    DO IREFER=1,LTPROT
C
                       LABTIT_PROTON(IREFER)=LABREF_PROTON(IREFER)
C
                       DO ITHEOR=1,LTPROT
                          IF (LABTHE_PROTON(ITHEOR)
     *                   .EQ. LABREF_PROTON(IREFER))THEN
C
                              VECTOR_ENEPRO(IREFER)
     *                       =ENETHE_PROTON(ITHEOR)
C
                              IF (ENETHE_PROTON(ITHEOR)
     *                       .LT. ENETHE_PROMIN) THEN
                                  ENETHE_PROMIN
     *                           =ENETHE_PROTON(ITHEOR)
                              END IF
C
                              IF (ENETHE_PROTON(ITHEOR)
     *                       .GT. ENETHE_PROMAX) THEN
                                  ENETHE_PROMAX
     *                           =ENETHE_PROTON(ITHEOR)
                              END IF
C
                          END IF
                       END DO
                    END DO
C
                    CALL PLOTCURVPRINT(NWRITE_ENEPRO,NDSPEC,LDCURV,
     *                                        IMONTE,LABTIT_PROTON,
     *                                               VECTOR_ENEPRO)
C
                END IF
C
                IF (IPRINT_NEUTRS.EQ.1) THEN
C
                    LDCURV=LTNEUT
C
                    DO IREFER=1,LTNEUT
C
                       LABTIT_NEUTRS(IREFER)=LABREF_NEUTRS(IREFER)
C
                       DO ITHEOR=1,LTNEUT
                          IF (LABTHE_NEUTRS(ITHEOR)
     *                   .EQ. LABREF_NEUTRS(IREFER)) THEN
C
                              VECTOR_ENENEU(IREFER)
     *                       =ENETHE_NEUTRS(ITHEOR)
C
                              IF (ENETHE_NEUTRS(ITHEOR)
     *                       .LT. ENETHE_NEUMIN) THEN
                                  ENETHE_NEUMIN
     *                           =ENETHE_NEUTRS(ITHEOR)
                              END IF
C
                              IF (ENETHE_NEUTRS(ITHEOR)
     *                       .GT. ENETHE_NEUMAX) THEN
                                  ENETHE_NEUMAX
     *                           =ENETHE_NEUTRS(ITHEOR)
                              END IF
C
                          END IF
                       END DO
                    END DO
C
                    CALL PLOTCURVPRINT(NWRITE_ENENEU,NDSPEC,LDCURV,
     *                                        IMONTE,LABTIT_NEUTRS,
     *                                               VECTOR_ENENEU)
C
                END IF
C
            END IF !ITAKNU
         END DO !INUCLI
C_______________________________________________________________________
C
C        Looking for the parameter extrems
C
         DO IPARAM=1,NDPARS
            IF (PARPOT_MAXIMA(IPARAM).LT.PARPOT(IPARAM)) THEN
                PARPOT_MAXIMA(IPARAM)=PARPOT(IPARAM)
            END IF
            IF (PARPOT_MINIMA(IPARAM).GT.PARPOT(IPARAM)) THEN
                PARPOT_MINIMA(IPARAM)=PARPOT(IPARAM)
            END IF
         END DO
C
C        Printing the parameters at each IMONTE
C
         NWRITE_PARAMS=1000
C
         DO IPARAM=1,NDPARS
            IF (IFTAKE(IPARAM).EQ.1) THEN
C
                NWRITE_PARAMS=NWRITE_PARAMS+1
                LDCURV=1
C
                LABTIT_PARPOT(LDCURV)=TITLES_LATEXS(IPARAM)
                VECTOR_PARPOT(LDCURV)=PARPOT(IPARAM)
C
                CALL PLOTCURVPRINT(NWRITE_PARAMS,NDPARS,LDCURV,
     *                                    IMONTE,LABTIT_PARPOT,
     *                                           VECTOR_PARPOT)
C
            END IF
         END DO
C
         I_SEED=NEWSED
C_______________________________________________________________________
C
      END DO !IMONTE
C
C=======================================================================
C=======================================================================
C
C     Printing RMS and SPE max and min values
C
      YMINIM=1.0
      YMAXIM=REAL(LDMONT)
C
      IF (IPRINT_PROTON.EQ.1) THEN
C
          CALL PLOTEXTRPRINT(NWRITE_RMSPRO,RMSVAL_PROMIN,
     *                       RMSVAL_PROMAX,YMINIM,YMAXIM)
C
          CALL PLOTEXTRPRINT(NWRITE_ENEPRO,ENETHE_PROMIN,
     *                       ENETHE_PROMAX,YMINIM,YMAXIM)
C
      END IF
C
      IF (IPRINT_NEUTRS.EQ.1) THEN
C
          CALL PLOTEXTRPRINT(NWRITE_RMSNEU,RMSVAL_NEUMIN,
     *                       RMSVAL_NEUMAX,YMINIM,YMAXIM)
C
          CALL PLOTEXTRPRINT(NWRITE_ENENEU,ENETHE_NEUMIN,
     *                       ENETHE_NEUMAX,YMINIM,YMAXIM)
C
      END IF
C
C     And now the parameters
C
      NWRITE_PARAMS=1000
C
      DO IPARAM=1,NDPARS
         IF (IFTAKE(IPARAM).EQ.1) THEN
C
             NWRITE_PARAMS=NWRITE_PARAMS+1
C
             XMAXIM=PARPOT_MAXIMA(IPARAM)
             XMINIM=PARPOT_MINIMA(IPARAM)
C
             CALL PLOTEXTRPRINT(NWRITE_PARAMS,XMINIM,XMAXIM,
     *                                        YMINIM,YMAXIM)
C
         END IF
      END DO
C
C=======================================================================
C
      CLOSE(NWRITE_PARAMS)
      CLOSE(NWRITE_ENENEU)
      CLOSE(NWRITE_ENEPRO)
      CLOSE(NWRITE_RMSNEU)
      CLOSE(NWRITE_RMSPRO)
C
C=======================================================================
C=======================================================================
C=======================================================================
C=======================================================================
C
      REWIND LOGPRO
   10 CONTINUE
      READ (LOGPRO,'(A120)',END=11) STRING
      IF (LOGWRI.GT.4) WRITE(NpUNIT,'(A80)'        ) STRING
C
      GO TO 10
C
   11 CONTINUE  ! Ready with the copy for the protons      
C      
C=======================================================================
C 
      REWIND LOGNEU    
   20 CONTINUE
      READ (LOGNEU,'(A120)',END=21) STRING
      IF (LOGWRI.GT.4) WRITE(NnUNIT,'(A80)'        ) STRING
C      
      GO TO 20
C
   21 CONTINUE  ! Ready with the copy for the neutrons     
C
C=======================================================================      
C=======================================================================
C
      WRITE(LOGFIL,'(/,9X,''Exting FITTIN_MONTEC'')')
C
C=======================================================================
C
      RETURN
      END
C
C=======================================================================
C=======================================================================
C
      SUBROUTINE EX_TO_PSEUDO
C
      INCLUDE 'MATDIM/NDNUCL.f'
      INCLUDE 'MATDIM/NDLEXP.f'
C
      CHARACTER
     *          NUCSYM*6,INPSYM*6
      CHARACTER
     *          LABPRO_PSEUDO*6,LABNEU_PSEUDO*6,
     *          LABEXP_PROTON*6,LABEXP_NEUTRS*6,
     *          LABEXP_PROREA*6,LABEXP_NEUREA*6
C
      COMMON
     *       /NUCINF/ LDNUCL,
     *                NUCSYM(1:NDNUCL),
     *                NUMB_Z(1:NDNUCL),
     *                NUMB_N(1:NDNUCL)
      COMMON
     *       /SELNUC/ NUCACT,
     *                ITAKNU(1:NDNUCL),
     *                INPUTZ(1:NDNUCL),
     *                INPUTN(1:NDNUCL),
     *                INPSYM(1:NDNUCL)
      COMMON
     *       /MONTEC/ IFPSEU,IFPARA,LDMONT,LDBINS
      COMMON
     *       /PSEUDO_ZN_NUM/ IZ_PSE(1:NDNUCL),
     *                       IN_PSE(1:NDNUCL)
      COMMON
     *       /PSEUDO_LEVELS/ LEVNEU_PSEUDO(1:NDNUCL),
     *                       LEVPRO_PSEUDO(1:NDNUCL)
      COMMON
     *       /PSEUDO_PROTON/ ENEPRO_PSEUDO(1:NDLEXP,1:NDNUCL),
     *                       LABPRO_PSEUDO(1:NDLEXP,1:NDNUCL),
     *                       SIGPRO_PSEUDO(1:NDLEXP,1:NDNUCL),
     *                       IDEPRO_PSEUDO(1:NDLEXP,1:NDNUCL),
     *                       LEVTAK_PROTON(1:NDLEXP,1:NDNUCL)
      COMMON
     *       /PSEUDO_NEUTRS/ ENENEU_PSEUDO(1:NDLEXP,1:NDNUCL),
     *                       LABNEU_PSEUDO(1:NDLEXP,1:NDNUCL),
     *                       SIGNEU_PSEUDO(1:NDLEXP,1:NDNUCL),
     *                       IDENEU_PSEUDO(1:NDLEXP,1:NDNUCL),
     *                       LEVTAK_NEUTRS(1:NDLEXP,1:NDNUCL)
      COMMON
     *       /SIGEXP_PRONEU/ SIGEXP_PROTON(1:NDLEXP,1:NDNUCL),
     *                       SIGEXP_NEUTRS(1:NDLEXP,1:NDNUCL)
      COMMON
     *       /EX_ENE/ EXPEXP_PROTON(1:NDNUCL,1:NDLEXP),
     *                IDEGEX_PROTON(1:NDNUCL,1:NDLEXP),
     *
     *                EXPEXP_NEUTRS(1:NDNUCL,1:NDLEXP),
     *                IDEGEX_NEUTRS(1:NDNUCL,1:NDLEXP)
      COMMON
     *       /EX_LAB/ LABEXP_PROTON(1:NDNUCL,1:NDLEXP),
     *                LABEXP_NEUTRS(1:NDNUCL,1:NDLEXP)
      COMMON
     *       /DIMSPH/ LDSING_PROTON,LEVTHE_PROTON,
     *                LEVEXP_PROTON(1:NDNUCL),
     *                LDSING_NEUTRS,LEVTHE_NEUTRS,
     *                LEVEXP_NEUTRS(1:NDNUCL)
      COMMON
     *       /REALEX_PROTON/ EXPEXP_PROREA(1:NDNUCL,1:NDLEXP),
     *                       LABEXP_PROREA(1:NDNUCL,1:NDLEXP),
     *                       IDEGEX_PROREA(1:NDNUCL,1:NDLEXP),
     *                       LEVEXP_PROREA(1:NDNUCL)
      COMMON
     *       /REALEX_NEUTRS/ EXPEXP_NEUREA(1:NDNUCL,1:NDLEXP),
     *                       LABEXP_NEUREA(1:NDNUCL,1:NDLEXP),
     *                       IDEGEX_NEUREA(1:NDNUCL,1:NDLEXP),
     *                       LEVEXP_NEUREA(1:NDNUCL)
C
      COMMON
     *       /PRFILE/ ISCREN,LOGFIL,LOGCHP,LOGCHN,LOGHIP,LOGHIN,
     *                LOGPRO,LOGNEU,LOGWRI,LSCREN,IRESUL,ICONVE,
     *                NOUTPT,LOGAUX,LOGENE,LOGRAD,LOGGRD,LOGMSH,
     *                                                   LOGBIS 
C
C=======================================================================
C
C     This subroutine transforms the real experimental single particle
C     energy values to the pseudo-experimental ones. It is called from
C     FITTIN_MONTEC.
C
C=======================================================================
C
      WRITE(LOGFIL,'(/,''Entering EX_TO_PSEUDO'')')
C
C=======================================================================
C
C     First, we store the real-experimental values 
C
      DO I=1,LDNUCL
C
         LEVEXP_PROREA(I)=LEVEXP_PROTON(I)
         DO J=1,LEVEXP_PROTON(I)
            EXPEXP_PROREA(I,J)=EXPEXP_PROTON(I,J)
            LABEXP_PROREA(I,J)=LABEXP_PROTON(I,J)
            IDEGEX_PROREA(I,J)=IDEGEX_PROTON(I,J)
         END DO
C
         LEVEXP_NEUREA(I)=LEVEXP_NEUTRS(I)
         DO J=1,LEVEXP_NEUTRS(I)
            EXPEXP_NEUREA(I,J)=EXPEXP_NEUTRS(I,J)
            LABEXP_NEUREA(I,J)=LABEXP_NEUTRS(I,J)
            IDEGEX_NEUREA(I,J)=IDEGEX_NEUTRS(I,J)
         END DO
C
      END DO
C
C=======================================================================
C
      DO I=1,LDNUCL
C
         IF (NUMB_Z(I).NE.IZ_PSE(I) .OR. NUMB_N(I).NE.IN_PSE(I)) THEN
             WRITE(LOGFIL,'(''WARNING in EX_TO_PSEUDO - not the same'',
     *                      ''nucleus: NUMB_Z= '',i3,'' and IZ_PSE= '',
     *                       I3,'' and NUMB_N= '',I3,'' and IN_PSE= '',
     *                       i3)')NUMB_Z(I),IZ_PSE(I),NUMB_N(I),
     *                                                IN_PSE(I)
             STOP 'STOP in EX_TO_PSEUDO: Wrong nucleus order'
         END IF
C
         IEXPER=0
         DO J=1,LEVPRO_PSEUDO(I)
            IF (LEVTAK_PROTON(J,I).EQ.1) THEN
                IEXPER=IEXPER+1
                EXPEXP_PROTON(I,IEXPER)=ENEPRO_PSEUDO(J,I)
                LABEXP_PROTON(I,IEXPER)=LABPRO_PSEUDO(J,I)
                IDEGEX_PROTON(I,IEXPER)=IDEPRO_PSEUDO(J,I)
                SIGEXP_PROTON(IEXPER,I)=SIGPRO_PSEUDO(J,I)
            END IF
         END DO
         LEVEXP_PROTON(I)=IEXPER
C
         IEXPER=0
         DO J=1,LEVNEU_PSEUDO(I)
            IF (LEVTAK_NEUTRS(J,I).EQ.1) THEN 
                IEXPER=IEXPER+1
                EXPEXP_NEUTRS(I,IEXPER)=ENENEU_PSEUDO(J,I)
                LABEXP_NEUTRS(I,IEXPER)=LABNEU_PSEUDO(J,I)
                IDEGEX_NEUTRS(I,IEXPER)=IDENEU_PSEUDO(J,I)
                SIGEXP_NEUTRS(IEXPER,I)=SIGNEU_PSEUDO(J,I)
            END IF
         END DO
         LEVEXP_NEUTRS(I)=IEXPER
C
         IF (LOGWRI.GT.0) THEN
             WRITE(LOGFIL,'(/,''For INUCLI= '',I1,'' we have '',
     *                        ''LEVEXP_PROTON= '',I3,'' and '',
     *                        ''LEVEXP_NEUTRS= '',I3,/)')
     *                        I,LEVEXP_PROTON(I),LEVEXP_NEUTRS(I)
         END IF
C
      END DO
C
C=======================================================================
C
      WRITE(LOGFIL,'(/,''We are working with the following '',
     *                                     ''pseudo-levels:'',/)')
C
      WRITE(LOGFIL,'(9X,''Protons:'',/)')
      DO INUCLI=1,LDNUCL
         IF (ITAKNU(INUCLI).EQ.1) THEN
C
             DO IEXPER=1,LEVEXP_PROTON(INUCLI)
                WRITE(LOGFIL,'(15X,2I4,F15.3,3X,A6,3X,I4,F15.4)')
     *                            INUCLI,IEXPER,
     *                            EXPEXP_PROTON(INUCLI,IEXPER),
     *                            LABEXP_PROTON(INUCLI,IEXPER),
     *                            IDEGEX_PROTON(INUCLI,IEXPER),
     *                            SIGEXP_PROTON(IEXPER,INUCLI)
             END DO
C
             WRITE(LOGFIL,'()')
C
         END IF
      END DO
C
      WRITE(LOGFIL,'(9X,''Neutrons:'',/)')
      DO INUCLI=1,LDNUCL
         IF (ITAKNU(INUCLI).EQ.1) THEN
C
             DO IEXPER=1,LEVEXP_NEUTRS(INUCLI)
                WRITE(LOGFIL,'(15X,2I4,F15.3,3X,A6,3X,I4,F15.4)')
     *                            INUCLI,IEXPER,
     *                            EXPEXP_NEUTRS(INUCLI,IEXPER),
     *                            LABEXP_NEUTRS(INUCLI,IEXPER),
     *                            IDEGEX_NEUTRS(INUCLI,IEXPER),
     *                            SIGEXP_NEUTRS(IEXPER,INUCLI)
             END DO
C
             WRITE(LOGFIL,'()')
C
         END IF
      END DO
C
C=======================================================================
C
      IF (LOGWRI.GT.0) WRITE(LOGFIL,'(''Exiting EX_TO_PSEUDO'')')
C
C=======================================================================
C
      RETURN
      END
C
C=======================================================================
C=======================================================================
C
      SUBROUTINE HISTOG_PARNAM(PARNAM,LEVACT,LEVACT_PROTON,
     *                         LEVACT_NEUTRS,SIGMAS,FILPAR)
C
      INCLUDE 'MATDIM/NDNUCL.f'
C
      CHARACTER
     *          PARNAM*8,FILPAR*300,VERSIO*3,TYPCHI*6,
     *          NUCSYM*6,INPSYM*6,FITNUC*2,
     *          FITNUC_AUXILI*2,FILNAM_NUCFIT*30,TEXLAM*20
C
      DIMENSION
     *          FITNUC(1:NDNUCL),
     *          FITNUC_AUXILI(1:NDNUCL)
C
      COMMON
     *       /NUCINF/ LDNUCL,
     *                NUCSYM(1:NDNUCL),
     *                NUMB_Z(1:NDNUCL),
     *                NUMB_N(1:NDNUCL)
       COMMON
     *       /SELNUC/ NUCACT,
     *                ITAKNU(1:NDNUCL),
     *                INPUTZ(1:NDNUCL),
     *                INPUTN(1:NDNUCL),
     *                INPSYM(1:NDNUCL)
      COMMON
     *       /DENSTR/ IFDENS
      COMMON
     *       /TENSOR/ IFTENS
      COMMON
     *       /TENOPT/ ICENTT,ISORBT,ITENSR 
      COMMON
     *       /DEPPAI/ IFDEEP,IFPRON
      COMMON
     *       /IFPAIR/ IF_PAI,IFGPAI,IFDPAI
      COMMON
     *       /IFLAMB/ IFPAR1,IFPAR2,IFPAR3,IFPAR4,
     *                IFPAR5,IFPAR6,IFPAR7,IFPAR8,
     *                IFPA09,IFPA10,IFPA11,IFPA12
      COMMON
     *       /MONTEC/ IFPSEU,IFPARA,LDMONT,LDBINS
      COMMON
     *       /IFKAPP/ IFK_VC,IFK_VS,IFK_RC,IFK_RS,IFK_AC,IFK_AS
      COMMON 
     *       /CHITAK/ IF_SPE,IF_RAD,IF_GAP,IF_FER,IF_DEN,IF_RHO,IF_INV,
     *                LDCHI2,TYPCHI(1:NDNUCL),ITKCHI(1:NDNUCL)
      COMMON
     *       /VERSIN/ VERSIO
C
      DATA
     *       FITNUC(01) / 'Ox' /,
     *       FITNUC(02) / 'C0' /,
     *       FITNUC(03) / 'C8' /,
     *       FITNUC(04) / 'Ni' /,
     *       FITNUC(05) / 'Zr' /,
     *       FITNUC(06) / 'Sn' /,
     *       FITNUC(07) / 'Gd' /,
     *       FITNUC(08) / 'Pb' /
C
C=======================================================================
C
      IF (IFDENS.EQ.1) THEN
         WRITE(TEXLAM,'(''_LAMBDA'',SP,I2,I2,I2,I2)')
     *                               IFPAR1,IFPAR2,IFPAR3,IFPAR4    
      END IF
C
C=======================================================================
C 
      KACTIV=0
      DO JNUCLI=1,LDNUCL
         IF (ITAKNU(JNUCLI).EQ.1) THEN
             KACTIV=KACTIV+1
             FITNUC_AUXILI(KACTIV)=FITNUC(JNUCLI)
         END IF
      END DO
      
      I1=3*KACTIV-1
      
      WRITE(FILNAM_NUCFIT,'(<NUCACT>(A2,''-''))')
     *                     (FITNUC_AUXILI(I),I=1,NUCACT)
C
C=======================================================================
C
      CALL SYSTEM('mkdir -pv ParHists/IFDENS-0')
      CALL SYSTEM('mkdir -pv ParHists/IFDENS-1_IFTENS-0')
      CALL SYSTEM('mkdir -pv ParHists/IFDENS-1_IFTENS-1')
C
C=======================================================================
C     
      IF (IFDENS.EQ.0) THEN
C      
          WRITE(FILPAR,'(''ParHists/IFDENS-0/'',A8,''_'',A,
     *                   ''_IFDENS-'',I1,''_IFTENS-'',I1,
     *                   ''_IF-PAI-'',I1,
     *                   ''_IF-RAD-'',I1,''_IF-ORD-'',I1,
     *                   ''_IF-RHO-'',I1,
     *                   ''_IFDEEP-'',I1,''_IFPRON-'',I1,
     *                   ''_IFK-VC-'',I1,''_IFK-VS-'',I1,
     *                   ''_IFK-RC-'',I1,''_IFK-RS-'',I1,
     *                   ''_IFK-AC-'',I1,''_IFK-AS-'',I1,
     *                   ''_LEVACT-'',I2.2,''_SIGMAS-'',F5.3,
     *                   ''_LDMONT-'',I5.5,''_'',A3)')
     *
     *                          PARNAM,FILNAM_NUCFIT(1:I1),
     *                          IFDENS,IFTENS,IF_PAI,IF_RAD,IF_INV,
     *                          IF_RHO,IFDEEP,IFPRON,IFK_VC,IFK_VS,
     *                          IFK_RC,IFK_RS,IFK_AC,IFK_AS,
     *                          LEVACT,SIGMAS,LDMONT,VERSIO
C      
      END IF
C
C=======================================================================
C      
      IF (IFDENS.EQ.1) THEN
C      
          WRITE(FILPAR,'(''ParHists/IFDENS-1_IFTENS-'',I1,''/'',
     *                                        A8,''_'',A,
     *                   ''_IFDENS-'',I1,''_IFTENS-'',I1,
     *                   ''_IF-PAI-'',I1,
     *                   ''_IF-RAD-'',I1,''_IF-ORD-'',I1,
     *                   ''_IFDEEP-'',I1,''_IFPRON-'',I1,
     *                   ''_IFK-VC-'',I1,''_IFK-VS-'',I1,
     *                   ''_IFK-RC-'',I1,''_IFK-RS-'',I1,
     *                   ''_IFK-AC-'',I1,''_IFK-AS-'',I1,
     *                   ''_LEVACp-'',I2.2,''_LEVACn-'',I2.2,
     *                                     ''_SIGMAS-'',F5.3,
     *                   ''_LDMONT-'',I5.5,A15,''_'',A3)')
     *
     *                      IFTENS,PARNAM,FILNAM_NUCFIT(1:I1),
     *                      IFDENS,IFTENS,IF_PAI,IF_RAD,IF_INV,
     *                      IFDEEP,IFPRON,IFK_VC,IFK_VS,
     *                      IFK_RC,IFK_RS,IFK_AC,IFK_AS,
     *                      LEVACT_PROTON,LEVACT_NEUTRS,SIGMAS,
     *                                    LDMONT,TEXLAM,VERSIO
C______________________________________________________________________
C
          IF (IFTENS.EQ.1) THEN ! IFDENS=1
C
              WRITE(FILPAR,'(''ParHists/IFDENS-1_IFTENS-'',I1,''/'',
     *                                        A8,''_'',A,
     *                   ''_IFDENS-'',I1,''_IFTENS-'',I1,
     *                   ''_IF-PAI-'',I1,
     *                   ''_ICENTT-'',I1,''_ISORBT-'',I1,
     *                                   ''_ITENSR-'',I1,
     *                   ''_IF-RAD-'',I1,''_IF-ORD-'',I1,
     *                   ''_IFDEEP-'',I1,''_IFPRON-'',I1,
     *                   ''_IFK-VC-'',I1,''_IFK-VS-'',I1,
     *                   ''_IFK-RC-'',I1,''_IFK-RS-'',I1,
     *                   ''_IFK-AC-'',I1,''_IFK-AS-'',I1,
     *                   ''_LEVACp-'',I2.2,''_LEVACn-'',I2.2,
     *                                     ''_SIGMAS-'',F5.3,
     *                   ''_LDMONT-'',I5.5,A15,''_'',A3)')
     *
     *                      IFTENS,PARNAM,FILNAM_NUCFIT(1:I1),
     *                      IFDENS,IFTENS,IF_PAI,ICENTT,ISORBT,ITENSR,
     *                      IF_RAD,IF_INV,IFDEEP,IFPRON,IFK_VC,IFK_VS,
     *                      IFK_RC,IFK_RS,IFK_AC,IFK_AS,
     *                      LEVACT_PROTON,LEVACT_NEUTRS,SIGMAS,
     *                                    LDMONT,TEXLAM,VERSIO
          END IF
C
      END IF
C
      DO K=1,300
         IF (FILPAR(K:K).EQ.'.') THEN
             FILPAR(K:K)='_'
         END IF
         IF (FILPAR(K:K).EQ.' ') THEN
             KLENGT=K-1
             GO TO 1
         END IF
      END DO
C
   1  CONTINUE
C
      FILPAR=FILPAR(1:KLENGT)//'.dat'
C
C=======================================================================
C
      RETURN
      END
C
C=======================================================================
C=======================================================================
C
      SUBROUTINE HISTOG_RMSNAM(ISOSPI,LEVACT,LEVACT_PROTON,
     *                         LEVACT_NEUTRS,SIGMAS,FILRMS)
C
      INCLUDE 'MATDIM/NDNUCL.f'
C
      CHARACTER
     *          FILISO*1,FILRMS*300,VERSIO*3,TYPCHI*6,
     *          FITNUC*2,FITNUC_AUXILI*2,FILNAM_NUCFIT*30,TEXLAM*20
      CHARACTER
     *          NUCSYM*6,INPSYM*6
C
      DIMENSION
     *          FITNUC(1:NDNUCL),
     *          FITNUC_AUXILI(1:NDNUCL)
C
      COMMON
     *       /NUCINF/ LDNUCL,
     *                NUCSYM(1:NDNUCL),
     *                NUMB_Z(1:NDNUCL),
     *                NUMB_N(1:NDNUCL)
      COMMON
     *       /SELNUC/ NUCACT,
     *                ITAKNU(1:NDNUCL),
     *                INPUTZ(1:NDNUCL),
     *                INPUTN(1:NDNUCL),
     *                INPSYM(1:NDNUCL)
      COMMON
     *       /DENSTR/ IFDENS
      COMMON
     *       /TENSOR/ IFTENS
      COMMON
     *       /TENOPT/ ICENTT,ISORBT,ITENSR 
      COMMON
     *       /DEPPAI/ IFDEEP,IFPRON
      COMMON
     *       /IFPAIR/ IF_PAI,IFGPAI,IFDPAI
      COMMON
     *       /IFLAMB/ IFPAR1,IFPAR2,IFPAR3,IFPAR4,
     *                IFPAR5,IFPAR6,IFPAR7,IFPAR8,
     *                IFPA09,IFPA10,IFPA11,IFPA12
      COMMON
     *       /MONTEC/ IFPSEU,IFPARA,LDMONT,LDBINS
      COMMON
     *       /IFKAPP/ IFK_VC,IFK_VS,IFK_RC,IFK_RS,IFK_AC,IFK_AS
      COMMON 
     *       /CHITAK/ IF_SPE,IF_RAD,IF_GAP,IF_FER,IF_DEN,IF_RHO,IF_INV,
     *                LDCHI2,TYPCHI(1:NDNUCL),ITKCHI(1:NDNUCL)
      COMMON
     *       /VERSIN/ VERSIO
C
      DATA
     *       FITNUC(01) / 'Ox' /,
     *       FITNUC(02) / 'C0' /,
     *       FITNUC(03) / 'C8' /,
     *       FITNUC(04) / 'Ni' /,
     *       FITNUC(05) / 'Zr' /,
     *       FITNUC(06) / 'Sn' /,
     *       FITNUC(07) / 'Gd' /,
     *       FITNUC(08) / 'Pb' /
C
C=======================================================================
C
      IF (ISOSPI.EQ.1) FILISO='P'
      IF (ISOSPI.EQ.0) FILISO='N'
C
C=======================================================================
C
      IF (IFDENS.EQ.1) THEN
         WRITE(TEXLAM,'(''_LAMBDA'',SP,I2,I2,I2,I2)')
     *                               IFPAR1,IFPAR2,IFPAR3,IFPAR4    
      END IF
C
C=======================================================================
C 
      KACTIV=0
      DO JNUCLI=1,LDNUCL
         IF (ITAKNU(JNUCLI).EQ.1) THEN
             KACTIV=KACTIV+1
             FITNUC_AUXILI(KACTIV)=FITNUC(JNUCLI)
         END IF
      END DO
      
      I1=3*KACTIV-1
      
      WRITE(FILNAM_NUCFIT,'(<NUCACT>(A2,''-''))')
     *                     (FITNUC_AUXILI(I),I=1,NUCACT)
C
C=======================================================================
C
      CALL SYSTEM('mkdir -pv RMSHists/IFDENS-0')
      CALL SYSTEM('mkdir -pv RMSHists/IFDENS-1_IFTENS-0')
      CALL SYSTEM('mkdir -pv RMSHists/IFDENS-1_IFTENS-1')
C
C=======================================================================
C
      IF (IFDENS.EQ.0) THEN
C      
          WRITE(FILRMS,'(''RMSHists/IFDENS-0/RMS-'',A1,''_'',A,
     *                   ''_IFDENS-'',I1,''_IFTENS-'',I1,
     *                   ''_IF-PAI-'',I1,
     *                   ''_IF-RAD-'',I1,''_IF-ORD-'',I1,
     *                   ''_IF-RHO-'',I1,
     *                   ''_IFDEEP-'',I1,''_IFPRON-'',I1,
     *                   ''_IFK-VC-'',I1,''_IFK-VS-'',I1,
     *                   ''_IFK-RC-'',I1,''_IFK-RS-'',I1,
     *                   ''_IFK-AC-'',I1,''_IFK-AS-'',I1,
     *                   ''_LEVACT-'',I2.2,''_SIGMAS-'',F5.3,
     *                   ''_LDMONT-'',I5.5,''_'',A3)')
     *
     *                          FILISO,FILNAM_NUCFIT(1:I1),
     *                          IFDENS,IFTENS,IF_PAI,IF_RAD,IF_INV,
     *                          IF_RHO,IFDEEP,IFPRON,IFK_VC,IFK_VS,
     *                          IFK_RC,IFK_RS,IFK_AC,IFK_AS,
     *                          LEVACT,SIGMAS,LDMONT,VERSIO
C
      END IF
C
C=======================================================================
C      
      IF (IFDENS.EQ.1) THEN
C      
          WRITE(FILRMS,'(''RMSHists/IFDENS-1_IFTENS-'',I1,''/RMS-'',
     *                                        A1,''_'',A,
     *                   ''_IFDENS-'',I1,''_IFTENS-'',I1,
     *                   ''_IF-PAI-'',I1,
     *                   ''_IF-RAD-'',I1,''_IF-ORD-'',I1,
     *                   ''_IFDEEP-'',I1,''_IFPRON-'',I1,
     *                   ''_IFK-VC-'',I1,''_IFK-VS-'',I1,
     *                   ''_IFK-RC-'',I1,''_IFK-RS-'',I1,
     *                   ''_IFK-AC-'',I1,''_IFK-AS-'',I1,
     *                   ''_LEVACp-'',I2.2,''_LEVACn-'',I2.2,
     *                                     ''_SIGMAS-'',F5.3,
     *                   ''_LDMONT-'',I5.5,A15,''_'',A3)')
     *
     *                      IFTENS,FILISO,FILNAM_NUCFIT(1:I1),
     *                      IFDENS,IFTENS,IF_PAI,IF_RAD,IF_INV,
     *                      IFDEEP,IFPRON,IFK_VC,IFK_VS,
     *                      IFK_RC,IFK_RS,IFK_AC,IFK_AS,
     *                      LEVACT_PROTON,LEVACT_NEUTRS,SIGMAS,
     *                                    LDMONT,TEXLAM,VERSIO
C______________________________________________________________________
C
          IF (IFTENS.EQ.1) THEN ! IFDENS=1
C
              WRITE(FILRMS,'(''RMSHists/IFDENS-1_IFTENS-'',I1,''/RMS-'',
     *                                        A1,''_'',A,
     *                   ''_IFDENS-'',I1,''_IFTENS-'',I1,
     *                   ''_IF-PAI-'',I1,
     *                   ''_ICENTT-'',I1,''_ISORBT-'',I1,
     *                                   ''_ITENSR-'',I1,
     *                   ''_IF-RAD-'',I1,''_IF-ORD-'',I1,
     *                   ''_IFDEEP-'',I1,''_IFPRON-'',I1,
     *                   ''_IFK-VC-'',I1,''_IFK-VS-'',I1,
     *                   ''_IFK-RC-'',I1,''_IFK-RS-'',I1,
     *                   ''_IFK-AC-'',I1,''_IFK-AS-'',I1,
     *                   ''_LEVACp-'',I2.2,''_LEVACn-'',I2.2,
     *                                     ''_SIGMAS-'',F5.3,
     *                   ''_LDMONT-'',I5.5,A15,''_'',A3)')
     *
     *                      IFTENS,FILISO,FILNAM_NUCFIT(1:I1),
     *                      IFDENS,IFTENS,IF_PAI,ICENTT,ISORBT,ITENSR,
     *                      IF_RAD,IF_INV,IFDEEP,IFPRON,IFK_VC,IFK_VS,
     *                      IFK_RC,IFK_RS,IFK_AC,IFK_AS,
     *                      LEVACT_PROTON,LEVACT_NEUTRS,SIGMAS,
     *                                    LDMONT,TEXLAM,VERSIO
          END IF
C
      END IF
C
      DO K=1,300
         IF (FILRMS(K:K).EQ.'.') THEN
             FILRMS(K:K)='_'
         END IF
         IF (FILRMS(K:K).EQ.' ') THEN
             KLENGT=K-1
             GO TO 1
         END IF
      END DO
C
   1  CONTINUE
C
      FILRMS=FILRMS(1:KLENGT)//'.dat'
C
C=======================================================================
C
      RETURN
      END
C
C=======================================================================
C=======================================================================
C
      SUBROUTINE HISTOG_SPENAM(ISOSPI,LEVACT,LEVACT_PROTON,
     *                         LEVACT_NEUTRS,SIGMAS,FILSPE)
C
      INCLUDE 'MATDIM/NDNUCL.f'
C
      CHARACTER
     *          FILISO*1,FILSPE*300,VERSIO*3,TYPCHI*6,
     *          FITNUC*2,FITNUC_AUXILI*2,FILNAM_NUCFIT*30,TEXLAM*20
      CHARACTER
     *          NUCSYM*6,INPSYM*6
C
      DIMENSION
     *          FITNUC(1:NDNUCL),
     *          FITNUC_AUXILI(1:NDNUCL)
C
      COMMON
     *       /NUCINF/ LDNUCL,
     *                NUCSYM(1:NDNUCL),
     *                NUMB_Z(1:NDNUCL),
     *                NUMB_N(1:NDNUCL)
      COMMON
     *       /SELNUC/ NUCACT,
     *                ITAKNU(1:NDNUCL),
     *                INPUTZ(1:NDNUCL),
     *                INPUTN(1:NDNUCL),
     *                INPSYM(1:NDNUCL)
      COMMON
     *       /DENSTR/ IFDENS
      COMMON
     *       /TENSOR/ IFTENS
      COMMON
     *       /TENOPT/ ICENTT,ISORBT,ITENSR
      COMMON
     *       /DEPPAI/ IFDEEP,IFPRON
      COMMON
     *       /IFPAIR/ IF_PAI,IFGPAI,IFDPAI
      COMMON
     *       /IFLAMB/ IFPAR1,IFPAR2,IFPAR3,IFPAR4,
     *                IFPAR5,IFPAR6,IFPAR7,IFPAR8,
     *                IFPA09,IFPA10,IFPA11,IFPA12
      COMMON
     *       /MONTEC/ IFPSEU,IFPARA,LDMONT,LDBINS
      COMMON
     *       /IFKAPP/ IFK_VC,IFK_VS,IFK_RC,IFK_RS,IFK_AC,IFK_AS
      COMMON
     *       /CHITAK/ IF_SPE,IF_RAD,IF_GAP,IF_FER,IF_DEN,IF_RHO,IF_INV,
     *                LDCHI2,TYPCHI(1:NDNUCL),ITKCHI(1:NDNUCL)
      COMMON
     *       /VERSIN/ VERSIO
C
      DATA
     *       FITNUC(01) / 'Ox' /,
     *       FITNUC(02) / 'C0' /,
     *       FITNUC(03) / 'C8' /,
     *       FITNUC(04) / 'Ni' /,
     *       FITNUC(05) / 'Zr' /,
     *       FITNUC(06) / 'Sn' /,
     *       FITNUC(07) / 'Gd' /,
     *       FITNUC(08) / 'Pb' /
C
C=======================================================================
C
      IF (ISOSPI.EQ.1) FILISO='P'
      IF (ISOSPI.EQ.0) FILISO='N'
C
C=======================================================================
C
      IF (IFDENS.EQ.1) THEN
          WRITE(TEXLAM,'(''_LAMBDA'',SP,I2,I2,I2,I2)')
     *                               IFPAR1,IFPAR2,IFPAR3,IFPAR4
      END IF
C
C=======================================================================
C
      KACTIV=0
      DO JNUCLI=1,LDNUCL
         IF (ITAKNU(JNUCLI).EQ.1) THEN
             KACTIV=KACTIV+1
             FITNUC_AUXILI(KACTIV)=FITNUC(JNUCLI)
          END IF
      END DO

      I1=3*KACTIV-1

      WRITE(FILNAM_NUCFIT,'(<NUCACT>(A2,''-''))')
     *                     (FITNUC_AUXILI(I),I=1,NUCACT)
C
C=======================================================================
C
      CALL SYSTEM('mkdir -pv SPEHists/IFDENS-0')
      CALL SYSTEM('mkdir -pv SPEHists/IFDENS-1_IFTENS-0')
      CALL SYSTEM('mkdir -pv SPEHists/IFDENS-1_IFTENS-1')
C
C=======================================================================
C
      IF (IFDENS.EQ.0) THEN
C
          WRITE(FILSPE,'(''SPEHists/IFDENS-0/SPE-'',A1,''_'',A,
     *                   ''_IFDENS-'',I1,''_IFTENS-'',I1,
     *                   ''_IF-PAI-'',I1,
     *                   ''_IF-RAD-'',I1,''_IF-ORD-'',I1,
     *                   ''_IF-RHO-'',I1,
     *                   ''_IFDEEP-'',I1,''_IFPRON-'',I1,
     *                   ''_IFK-VC-'',I1,''_IFK-VS-'',I1,
     *                   ''_IFK-RC-'',I1,''_IFK-RS-'',I1,
     *                   ''_IFK-AC-'',I1,''_IFK-AS-'',I1,
     *                   ''_LEVACT-'',I2.2,''_SIGMAS-'',F5.3,
     *                   ''_LDMONT-'',I5.5,''_'',A3)')
     *
     *                          FILISO,FILNAM_NUCFIT(1:I1),
     *                          IFDENS,IFTENS,IF_PAI,IF_RAD,IF_INV,
     *                          IF_RHO,IFDEEP,IFPRON,IFK_VC,IFK_VS,
     *                          IFK_RC,IFK_RS,IFK_AC,IFK_AS,
     *                          LEVACT,SIGMAS,LDMONT,VERSIO
C
      END IF
C
C=======================================================================
C
      IF (IFDENS.EQ.1) THEN
C
          WRITE(FILSPE,'(''SPEHists/IFDENS-1_IFTENS-'',I1,''/SPE-'',
     *                                        A1,''_'',A,
     *                   ''_IFDENS-'',I1,''_IFTENS-'',I1,
     *                   ''_IF-PAI-'',I1,
     *                   ''_IF-RAD-'',I1,''_IF-ORD-'',I1,
     *                   ''_IFDEEP-'',I1,''_IFPRON-'',I1,
     *                   ''_IFK-VC-'',I1,''_IFK-VS-'',I1,
     *                   ''_IFK-RC-'',I1,''_IFK-RS-'',I1,
     *                   ''_IFK-AC-'',I1,''_IFK-AS-'',I1,
     *                   ''_LEVACp-'',I2.2,''_LEVACn-'',I2.2,
     *                                     ''_SIGMAS-'',F5.3,
     *                   ''_LDMONT-'',I5.5,A15,''_'',A3)')
     *
     *                      IFTENS,FILISO,FILNAM_NUCFIT(1:I1),
     *                      IFDENS,IFTENS,IF_PAI,IF_RAD,IF_INV,
     *                      IFDEEP,IFPRON,IFK_VC,IFK_VS,
     *                      IFK_RC,IFK_RS,IFK_AC,IFK_AS,
     *                      LEVACT_PROTON,LEVACT_NEUTRS,SIGMAS,
     *                                    LDMONT,TEXLAM,VERSIO
C______________________________________________________________________
C
          IF (IFTENS.EQ.1) THEN ! IFDENS=1
C
              WRITE(FILSPE,'(''SPEHists/IFDENS-1_IFTENS-'',I1,''/SPE-'',
     *                                        A1,''_'',A,
     *                   ''_IFDENS-'',I1,''_IFTENS-'',I1,
     *                   ''_IF-PAI-'',I1,
     *                   ''_ICENTT-'',I1,''_ISORBT-'',I1,
     *                                   ''_ITENSR-'',I1,
     *                   ''_IF-RAD-'',I1,''_IF-ORD-'',I1,
     *                   ''_IFDEEP-'',I1,''_IFPRON-'',I1,
     *                   ''_IFK-VC-'',I1,''_IFK-VS-'',I1,
     *                   ''_IFK-RC-'',I1,''_IFK-RS-'',I1,
     *                   ''_IFK-AC-'',I1,''_IFK-AS-'',I1,
     *                   ''_LEVACp-'',I2.2,''_LEVACn-'',I2.2,
     *                                     ''_SIGMAS-'',F5.3,
     *                   ''_LDMONT-'',I5.5,A15,''_'',A3)')
     *
     *                      IFTENS,FILISO,FILNAM_NUCFIT(1:I1),
     *                      IFDENS,IFTENS,IF_PAI,ICENTT,ISORBT,ITENSR,
     *                      IF_RAD,IF_INV,IFDEEP,IFPRON,IFK_VC,IFK_VS,
     *                      IFK_RC,IFK_RS,IFK_AC,IFK_AS,
     *                      LEVACT_PROTON,LEVACT_NEUTRS,SIGMAS,
     *                                    LDMONT,TEXLAM,VERSIO
          END IF
C
      END IF
C
      DO K=1,300
         IF (FILSPE(K:K).EQ.'.') THEN
             FILSPE(K:K)='_'
         END IF
         IF (FILSPE(K:K).EQ.' ') THEN
             KLENGT=K-1
             GO TO 1
         END IF
      END DO
C
   1  CONTINUE
C
      FILSPE=FILSPE(1:KLENGT)//'.dat'
C
C=======================================================================
C
      RETURN
      END
C
C=======================================================================
C=======================================================================
C
      SUBROUTINE SPENER_MONTEC
C
      INCLUDE 'MATDIM/NDMESH.f'
      INCLUDE 'MATDIM/NDSPEC.f'
      INCLUDE 'MATDIM/NDPARS.f'
      INCLUDE 'MATDIM/NDNUCL.f'
      INCLUDE 'MATDIM/NDLEXP.f'
      INCLUDE 'MATDIM/NDIM_M.f'
      INCLUDE 'MATDIM/NDTITL.f'
      INCLUDE 'MATDIM/NDLAST.f'
      INCLUDE 'MATDIM/NDBINS.f'
      INCLUDE 'MATDIM/NDMONT.f'
      INCLUDE 'MATDIM/NDCOLO.f'
C
      PARAMETER
     *          (NDMES2=NDMESH*NDMESH)
C
      CHARACTER
     *          NUCSYM*6,INPSYM*6,TYPCHI*6,VERSIO*3,FILNAM*256
      CHARACTER
     *          TITLES*12,TITLES_LATEXS*050,NUCNAM_LATEXS*010
      CHARACTER
     *          LABTHE_PROTON*6,LABTHE_NEUTRS*6,
     *          LABEXP_PROTON*6,LABEXP_NEUTRS*6
      CHARACTER
     *          LABREF*6,LABTHE*6,LABEXP*6,LABTEX*11,LABELS*100
      CHARACTER
     *          LABREF_PROTON*6,LABREF_NEUTRS*6
      CHARACTER
     *          HISTIT*256,HX_TIT*256,HY_TIT*256
C
      DIMENSION
     *          ARGMNT(1:NDPARS)
      DIMENSION
     *          SPEVEC_PROTON(1:NDMONT,1:NDSPEC),
     *          SPEVEC_NEUTRS(1:NDMONT,1:NDSPEC)
      DIMENSION
     *          LABREF_PROTON(1:NDSPEC),
     *          LABREF_NEUTRS(1:NDSPEC)
      DIMENSION
     *          LABELS(1:NDBINS,1:NDSPEC),
     *          XHISTO(1:NDBINS,1:NDSPEC),
     *          XHISTO_MINIMS(1:NDSPEC),
     *          XHISTO_MAXIMS(1:NDSPEC)
      DIMENSION
     *          BINSIZ(1:NDSPEC)
      DIMENSION
     *          YHISTO_UNNORM(1:NDBINS,1:NDSPEC),
     *          YHISTO_INTEGR(1:NDBINS,1:NDSPEC),
     *          YHISTO_MAXIMS(1:NDBINS,1:NDSPEC)
      DIMENSION
     *          YEXPEC(1:NDSPEC),
     *          ICOLOR_HISTOG(1:NDSPEC)
C      
      COMMON
     *       /NUCLEU/ IZ_FIX,IN_FIX,ISOSPI 
      COMMON
     *       /DENSTR/ IFDENS
      COMMON
     *       /TENSOR/ IFTENS
      COMMON
     *       /TENOPT/ ICENTT,ISORBT,ITENSR 
      COMMON
     *       /DEPPAI/ IFDEEP,IFPRON
      COMMON
     *       /IFPAIR/ IF_PAI,IFGPAI,IFDPAI
      COMMON 
     *       /CHITAK/ IF_SPE,IF_RAD,IF_GAP,IF_FER,IF_DEN,IF_RHO,IF_INV,
     *                LDCHI2,TYPCHI(1:NDNUCL),ITKCHI(1:NDNUCL)
      COMMON
     *       /IFLAMB/ IFPAR1,IFPAR2,IFPAR3,IFPAR4,
     *                IFPAR5,IFPAR6,IFPAR7,IFPAR8,
     *                IFPA09,IFPA10,IFPA11,IFPA12
      COMMON
     *       /PARCOR/ IFCORR,IFRCVC,IFRCAC,IFVCAC,IFRSVS,IFRSAS,IFVSAS

      COMMON
     *       /POTPOT/ PARPOT(1:NDPARS)
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
     *       /VLIMIT/ VMISTR(1:NDPARS),
     *                VMIMIN(1:NDPARS),
     *                VMIMAX(1:NDPARS),
     *                VMISTP(1:NDPARS)
      COMMON
     *       /SEARCH/ IFTAKE(1:NDPARS)
      COMMON
     *       /SPECTR/ ENETHE_PROTON(1:NDSPEC),
     *                ENETHE_NEUTRS(1:NDSPEC)
      COMMON
     *       /THELAB/ LABTHE_PROTON(1:NDSPEC),
     *                LABTHE_NEUTRS(1:NDSPEC)  
      COMMON
     *       /QUANTN/ NWSSPH_PROTON(1:NDSPEC),
     *                LWSSPH_PROTON(1:NDSPEC),
     *                JWSSPH_PROTON(1:NDSPEC),
     *
     *                NWSSPH_NEUTRS(1:NDSPEC),
     *                LWSSPH_NEUTRS(1:NDSPEC),
     *                JWSSPH_NEUTRS(1:NDSPEC)
      COMMON
     *       /EX_ENE/ EXPEXP_PROTON(1:NDNUCL,1:NDLEXP),
     *                IDEGEX_PROTON(1:NDNUCL,1:NDLEXP),
     *
     *                EXPEXP_NEUTRS(1:NDNUCL,1:NDLEXP),
     *                IDEGEX_NEUTRS(1:NDNUCL,1:NDLEXP)
      COMMON
     *       /EX_LAB/ LABEXP_PROTON(1:NDNUCL,1:NDLEXP),
     *                LABEXP_NEUTRS(1:NDNUCL,1:NDLEXP)
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
     *       /AUXPOT/ VKACEN_PROTON,RKACEN_PROTON,AKACEN_PROTON,
     *                VKASOR_PROTON,RKASOR_PROTON,AKASOR_PROTON,
     *                VKAEFM_PROTON,RKAEFM_PROTON,AKAEFM_PROTON,
     *                                            V0UNIT_PROTON,
     *                VKACEN_NEUTRS,RKACEN_NEUTRS,AKACEN_NEUTRS,
     *                VKASOR_NEUTRS,RKASOR_NEUTRS,AKASOR_NEUTRS,
     *                VKAEFM_NEUTRS,RKAEFM_NEUTRS,AKAEFM_NEUTRS,
     *                                            V0UNIT_NEUTRS,
     *                                                   RKACOU,
     *                              ALAMPP,ALAMPN,ALAMNP,ALAMNN,
     *                              TLAMPP,TLAMPN,TLAMNP,TLAMNN,
     *                              CLAMPP,CLAMPN,CLAMNP,CLAMNN
      COMMON
     *       /MONTEC/ IFPSEU,IFPARA,LDMONT,LDBINS
      COMMON
     *       /GAUSSI/ PARPOT_XMEANS(1:NDPARS),
     *                PARPOT_SIGMAS(1:NDPARS)
      COMMON
     *       /PARAMT_TITLES/ TITLES(1:NDTITL),
     *                       TITLES_LATEXS(1:NDTITL),
     *                       NUCNAM_LATEXS(1:NDNUCL)
      COMMON
     *       /ATTRIB/ ISTYLE,I_TYPE,ICOLOR(1:NDCOLO),ITHICK(1:NDCOLO)
      COMMON
     *       /VERSIN/ VERSIO
      COMMON
     *       /PRFILE/ ISCREN,LOGFIL,LOGCHP,LOGCHN,LOGHIP,LOGHIN,
     *                LOGPRO,LOGNEU,LOGWRI,LSCREN,IRESUL,ICONVE,
     *                NOUTPT,LOGAUX,LOGENE,LOGRAD,LOGGRD,LOGMSH,
     *                                                   LOGBIS 
C
C=======================================================================
C
C     This subroutine creates the histograms of the SPE levels .
C     It uses Monte-Carlo random generator to create LDMONT random
C     restarts. 
C
C     The random numbers are choosen using PARPOT_XMEANS
C     and PARPOT_SIGMAS distributed in a Gaussian Distribution.
C
C=======================================================================
C
                   INSEED=0!7893
      CALL ZBQLINI(INSEED)
C
      DO INUCLI=1,LDNUCL
         IF (ITAKNU(INUCLI).EQ.1) THEN
C
             IZ_FIX=NUMB_Z(INUCLI)
             IN_FIX=NUMB_N(INUCLI)
C
             CALL LOGMOC_CREATN(INUCLI,FILNAM)
             LOGMOC=44
             OPEN(UNIT=LOGMOC,FILE=FILNAM,STATUS='UNKNOWN')
C
             WRITE(LOGMOC,'(''  IMONTE'',$)')
             DO IPARAM=1,NDPARS
                IF (IFTAKE(IPARAM).EQ.1) THEN
                    WRITE(LOGMOC,'(3X,A8,$)')TITLES(IPARAM)(3:10)
                END IF
             END DO
             WRITE(LOGMOC,'(3X,''RMSVAL_p'',3x,''RMSVAL_n'')')
C
             DO IMONTE=1,LDMONT
C
                IF (LOGWRI.GT.0) THEN
                    WRITE(NOUTPT,'(/,80(''#''),/,''#'',T80,''#'',/,
     *                             ''#  IMONTE= '',I6,T80,''#'',/,
     *                             ''#'',T80,''#'',/,80(''#''))')IMONTE
                END IF
C
                IACTIV=0
                DO IPARAM=1,NDPARS
C
                   IF (IFTAKE(IPARAM).EQ.1) THEN
C
                       IACTIV=IACTIV+1
C
                       GAUSIG=PARPOT_SIGMAS(IPARAM)
                       GAUSMU=PARPOT_XMEANS(IPARAM)
C
                       RANDOM=ZBQLNOR(GAUSMU,GAUSIG)
                       ARGMNT(IACTIV)=RANDOM
                       PARPOT(IPARAM)=RANDOM
C
                   END IF
C
                END DO !IPARAM
C
                WRITE(LOGMOC,'(I8,$)')IMONTE
                IACTIV=0
                DO IPARAM=1,NDPARS
                   IF (IFTAKE(IPARAM).EQ.1) THEN
                       IACTIV=IACTIV+1
                       WRITE(LOGMOC,'(3X,F8.4,$)')PARPOT(IPARAM)
                   END IF
                END DO
                WRITE(LOGMOC,'()')
C
C               Being sure that the random numbers have physical sense
C
                CALL VERIFY(ARGMNT,IABORT)
C
                IACTIV=0
                DO IPARAM=1,NDPARS
                   IF (IFTAKE(IPARAM).EQ.1) THEN
                       IACTIV=IACTIV+1
                       PARPOT(IPARAM)=ARGMNT(IACTIV)
                   END IF
                END DO
C
                WRITE(LOGMOC,'(I8,$)')IMONTE
                DO IPARAM=1,NDPARS
                   IF (IFTAKE(IPARAM).EQ.1) THEN
                       WRITE(LOGMOC,'(3X,F8.4,$)')PARPOT(IPARAM)
                   END IF
                END DO
C
C               Levels Calculation
C
                IDEFCN=0
                I_MODE=0
                I_FLAG=1
C
                CALL WS_RUN(INUCLI,CHISQU_CHOICE,IABORT,I_MODE,
     *                      I_FLAG,CHISQU_PROTON,CHISQU_NEUTRS)
C
                WRITE(LOGMOC,'(3X,F8.4,3X,F8.4)')SQRT(CHIWEI_PROTON),
     *                                           SQRT(CHIWEI_NEUTRS)
C
                LTPROT=LEVTHE_PROTON
                LTNEUT=LEVTHE_NEUTRS
C
C               Storing the reference label
C
                IF (IMONTE.EQ.1) THEN
                    DO IREFER=1,LTPROT
                       LABREF_PROTON(IREFER)=LABTHE_PROTON(IREFER)
                    END DO
                    DO IREFER=1,LTNEUT
                       LABREF_NEUTRS(IREFER)=LABTHE_NEUTRS(IREFER)
                    END DO
                END IF
C
C               Storing the results for each IMONTE (=restart)
C
                DO ITHEOR=1,LTPROT
                   DO IREFER=1,LTPROT
                      IF (LABTHE_PROTON(ITHEOR)
     *               .EQ. LABREF_PROTON(IREFER))THEN
                          SPEVEC_PROTON(IMONTE,ITHEOR)
     *                   =ENETHE_PROTON(ITHEOR)
                      END IF  
                   END DO       
                END DO
C
                DO ITHEOR=1,LTNEUT
                   DO IREFER=1,LTNEUT
                      IF (LABTHE_NEUTRS(ITHEOR)
     *               .EQ. LABREF_NEUTRS(IREFER)) THEN
                          SPEVEC_NEUTRS(IMONTE,ITHEOR)
     *                   =ENETHE_NEUTRS(ITHEOR)
                      END IF
                   END DO
                END DO
C 
             END DO !IMONTE
C_______________________________________________________________________
C
C            Once we have the LDMONT different results, 
C            we construct the histograms
C_______________________________________________________________________
C
C            Protons
C
             ISOSPI=1
C
             CALL HISTOG_FABRIC(LTPROT,SPEVEC_PROTON,LDHIST,BINSIZ,
     *                          XHISTO,XHISTO_MINIMS,XHISTO_MAXIMS,
     *                                 YHISTO_UNNORM,YHISTO_INTEGR,
     *                                               YHISTO_MAXIMS)
C
             DO ITHEOR=1,LTPROT
                DO IHISTO=1,LDHIST
                   LABTHE=LABREF_PROTON(ITHEOR)
                   CALL LATEXS_LABELS(LABTHE,LABTEX)
                   LABELS(IHISTO,ITHEOR)=LABTEX
                END DO
             END DO
C
C            Identifying the EXPECTED values (for plotting purposes)
C
             II=0
             DO ITHEOR=1,LTPROT
                YEXPEC(ITHEOR)=999.0
                ICOLOR_HISTOG(ITHEOR)=70
                DO IEXPER=1,LEVEXP_PROTON(INUCLI)
                   IF (LABREF_PROTON(ITHEOR)
     *            .EQ. LABEXP_PROTON(INUCLI,IEXPER)) THEN
                       YEXPEC(ITHEOR)=EXPEXP_PROTON(INUCLI,IEXPER)
                       II=II+1
                       ICOLOR_HISTOG(ITHEOR)=ICOLOR(II)
                   END IF
                END DO
             END DO
C
             HISTIT='Proton Energy Distributions'
             HX_TIT='Energy (MeV)'
             HY_TIT='Probability Distribution'
C
             CALL HISTOG_PRINTI(LTPROT,LDHIST,XHISTO,YHISTO_MAXIMS,
     *                          XHISTO_MINIMS,XHISTO_MAXIMS,YEXPEC,
     *                          HISTIT,HX_TIT,HY_TIT,LABELS,INUCLI,
     *                                               ICOLOR_HISTOG)
C_______________________________________________________________________
C
C            Neutrons
C
             ISOSPI=0
C
             CALL HISTOG_FABRIC(LTNEUT,SPEVEC_NEUTRS,LDHIST,BINSIZ,
     *                          XHISTO,XHISTO_MINIMS,XHISTO_MAXIMS,
     *                                 YHISTO_UNNORM,YHISTO_INTEGR,
     *                                               YHISTO_MAXIMS)
C
             DO ITHEOR=1,LTNEUT
                DO IHISTO=1,LDHIST
                   LABTHE=LABREF_NEUTRS(ITHEOR)
                   CALL LATEXS_LABELS(LABTHE,LABTEX)
                   LABELS(IHISTO,ITHEOR)=LABTEX
                END DO
             END DO
C
C            Identifying the EXPECTED values (for plotting purposes)
C
             II=0
             DO ITHEOR=1,LTNEUT
                YEXPEC(ITHEOR)=999.0
                ICOLOR_HISTOG(ITHEOR)=70
                DO IEXPER=1,LEVEXP_NEUTRS(INUCLI)
                   IF (LABREF_NEUTRS(ITHEOR)
     *            .EQ. LABEXP_NEUTRS(INUCLI,IEXPER)) THEN
                       YEXPEC(ITHEOR)=EXPEXP_NEUTRS(INUCLI,IEXPER)
                       II=II+1
                       ICOLOR_HISTOG(ITHEOR)=ICOLOR(II)
                   END IF
                END DO
             END DO
C
             HISTIT='Neutron Energy Distributions'
             HX_TIT='Energy (MeV)'
             HY_TIT='Probability Distribution'
C
             CALL HISTOG_PRINTI(LTNEUT,LDHIST,XHISTO,YHISTO_MAXIMS,
     *                          XHISTO_MINIMS,XHISTO_MAXIMS,YEXPEC,
     *                          HISTIT,HX_TIT,HY_TIT,LABELS,INUCLI,
     *                                               ICOLOR_HISTOG)
C
         END IF !ITAKNU
      END DO !INUCLI
C
C=======================================================================
C
      RETURN
      END
C
C=======================================================================
C=======================================================================
C
      SUBROUTINE LOGMOC_CREATN(INUCLI,FILNAM)
C
      INCLUDE 'MATDIM/NDNUCL.f'
C
      CHARACTER
     *          FILNAM*256,NUCSYM*6,TYPCHI*6,TEXLAM*20,NUCNAM*6,
     *                                                 VERSIO*3
C      
      COMMON
     *       /NUCINF/ LDNUCL,
     *                NUCSYM(1:NDNUCL),
     *                NUMB_Z(1:NDNUCL),
     *                NUMB_N(1:NDNUCL)
      COMMON
     *       /NUCLEU/ IZ_FIX,IN_FIX,ISOSPI 
      COMMON
     *       /DENSTR/ IFDENS
      COMMON
     *       /TENSOR/ IFTENS
      COMMON
     *       /TENOPT/ ICENTT,ISORBT,ITENSR 
      COMMON
     *       /DEPPAI/ IFDEEP,IFPRON
      COMMON
     *       /IFPAIR/ IF_PAI,IFGPAI,IFDPAI
      COMMON 
     *       /CHITAK/ IF_SPE,IF_RAD,IF_GAP,IF_FER,IF_DEN,IF_RHO,IF_INV,
     *                LDCHI2,TYPCHI(1:NDNUCL),ITKCHI(1:NDNUCL)
      COMMON
     *       /IFLAMB/ IFPAR1,IFPAR2,IFPAR3,IFPAR4,
     *                IFPAR5,IFPAR6,IFPAR7,IFPAR8,
     *                IFPA09,IFPA10,IFPA11,IFPA12
      COMMON
     *       /PARCOR/ IFCORR,IFRCVC,IFRCAC,IFVCAC,IFRSVS,IFRSAS,IFVSAS
      COMMON
     *       /VERSIN/ VERSIO
C
C=======================================================================
C
      WRITE(NUCNAM,'(A6)') NUCSYM(INUCLI)
C
      J1=1
      J2=6
      IF (NUCNAM(1:1).EQ.' ') J1=2
      IF (NUCNAM(2:2).EQ.' ') J1=3
      IF (NUCNAM(3:3).EQ.' ') J1=4
C
C=======================================================================
C
      IF (IFDENS.EQ.0) THEN
          WRITE(FILNAM,'(''LogFMont/IFDENS-0/SPE-MC_'',A,
     *                  ''_IFDENS-'',I1,''_IFTENS-'',I1,''_IF-PAI-'',I1,
     *                  ''_IF-RAD-'',I1,''_IF-ORD-'',I1,''_IFDEEP-'',I1,
     *                  ''_IFPRON-'',I1,
     *                  ''_IFCORR-'',I1,''_IFRCVC-'',I1,''_IFRCAC-'',I1,
     *                  ''_IFVCAC-'',I1,''_IFRSVS-'',I1,''_IFRSAS-'',I1,
     *                  ''_IFVSAS-'',I1,
     *                  ''_'',A3,''.log'')')
     *
     *               NUCNAM(J1:J2),IFDENS,IFTENS,IF_PAI,
     *               IF_RAD,IF_INV,IFDEEP,IFPRON,IFCORR,IFRCVC,
     *               IFRCAC,IFVCAC,IFRSVS,IFRSAS,IFVSAS,VERSIO
      END IF
C
C=======================================================================
C
      IF (IFDENS.EQ.1 .AND. IFTENS.EQ.0) THEN
C
          WRITE(TEXLAM,'(''_LAMBDA'',SP,I2,I2,I2,I2)')
     *                            IFPAR1,IFPAR2,IFPAR3,IFPAR4
C
          WRITE(FILNAM,'(''LogFMont/IFDENS-1_IFTENS-0/SPE-MC_'',A,
     *                  ''_IFDENS-'',I1,''_IFTENS-'',I1,''_IF-PAI-'',I1,
     *                  ''_IF-RAD-'',I1,''_IF-ORD-'',I1,''_IFDEEP-'',I1,
     *                  ''_IFPRON-'',I1,
     *                  ''_IFCORR-'',I1,''_IFRCVC-'',I1,''_IFRCAC-'',I1,
     *                  ''_IFVCAC-'',I1,''_IFRSVS-'',I1,''_IFRSAS-'',I1,
     *                  ''_IFVSAS-'',I1,
     *                  A15,''_'',A3,''.log'')')
     *
     *               NUCNAM(J1:J2),
     *               IFDENS,IFTENS,IF_PAI,IF_RAD,IF_INV,IFDEEP,IFPRON,
     *               IFCORR,IFRCVC,IFRCAC,IFVCAC,IFRSVS,IFRSAS,IFVSAS,
     *                                                  TEXLAM,VERSIO
      END IF
C
C=======================================================================
C
      IF (IFDENS.EQ.1 .AND. IFTENS.EQ.1) THEN
C
          WRITE(TEXLAM,'(''_LAMBDA'',SP,I2,I2,I2,I2)')
     *                            IFPAR1,IFPAR2,IFPAR3,IFPAR4
C
          WRITE(FILNAM,'(''LogFMont/IFDENS-1_IFTENS-0/SPE-MC_'',A,
     *                  ''_IFDENS-'',I1,''_IFTENS-'',I1,''_IF-PAI-'',I1,
     *                  ''_ICENTT-'',I1,''_ISORBT-'',I1,''_ITENSR-'',I1,
     *                  ''_IF-RAD-'',I1,''_IF-ORD-'',I1,''_IFDEEP-'',I1,
     *                  ''_IFPRON-'',I1,
     *                  ''_IFCORR-'',I1,''_IFRCVC-'',I1,''_IFRCAC-'',I1,
     *                  ''_IFVCAC-'',I1,''_IFRSVS-'',I1,''_IFRSAS-'',I1,
     *                  ''_IFVSAS-'',I1,
     *                  A15,''_'',A3,''.log'')')
     *
     *               NUCNAM(J1:J2),
     *               IFDENS,IFTENS,IF_PAI,ICENTT,ISORBT,ITENSR,
     *               IF_RAD,IF_INV,IFDEEP,IFPRON,IFCORR,IFRCVC,
     *               IFRCAC,IFVCAC,IFRSVS,IFRSAS,IFVSAS,TEXLAM,VERSIO
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
      SUBROUTINE HISTOG_FABRIC(LDCURV,SPEVEC,LDHIST,BINSIZ,XHISTO,
     *                                XHISTO_MINIMS,XHISTO_MAXIMS,
     *                                YHISTO_UNNORM,YHISTO_INTEGR,
     *                                              YHISTO_MAXIMS)
C
      INCLUDE 'MATDIM/NDSPEC.f'
      INCLUDE 'MATDIM/NDMONT.f'
      INCLUDE 'MATDIM/NDBINS.f'
C
      DIMENSION
     *          SPEVEC(1:NDMONT,1:NDSPEC)
      DIMENSION
     *          XAUXIL_HISTOG(1:NDBINS,1:NDSPEC),
     *          YAUXIL_HISTOG(1:NDBINS,1:NDSPEC)
      DIMENSION
     *          XHISTO(1:NDBINS,1:NDSPEC),
     *          XHISTO_MINIMS(1:NDSPEC),
     *          XHISTO_MAXIMS(1:NDSPEC)
      DIMENSION
     *          YHISTO_UNNORM(1:NDBINS,1:NDSPEC),
     *          YHISTO_INTEGR(1:NDBINS,1:NDSPEC),
     *          YHISTO_MAXIMS(1:NDBINS,1:NDSPEC)
      DIMENSION
     *          BINSIZ(1:NDSPEC)
      COMMON
     *       /MONTEC/ IFPSEU,IFPARA,LDMONT,LDBINS
      COMMON
     *       /PRFILE/ ISCREN,LOGFIL,LOGCHP,LOGCHN,LOGHIP,LOGHIN,
     *                LOGPRO,LOGNEU,LOGWRI,LSCREN,IRESUL,ICONVE,
     *                NOUTPT,LOGAUX,LOGENE,LOGRAD,LOGGRD,LOGMSH,
     *                                                   LOGBIS 
C
C=======================================================================
C
      IF (LOGWRI.GT.0) THEN
          WRITE(LOGFIL,'(''Entering to HISTOG_FABRIC'')')
      END IF
C
C=======================================================================
C
      DO K=1,NDBINS
         DO J=1,LDCURV
            YHISTO_UNNORM(K,J)=0.0
            YHISTO_INTEGR(K,J)=0.0
            YHISTO_MAXIMS(K,J)=0.0
            YAUXIL_HISTOG(K,J)=0.0
        END DO
      END DO
C
C=======================================================================
C
      DO J=1,LDCURV
C
         XMINIM=+1.0E+10
         XMAXIM=-1.0E+10
C
         DO I=1,LDMONT
C
            IF (XMINIM.GT.SPEVEC(I,J)) THEN
                XMINIM=SPEVEC(I,J)
            END IF
C
            IF (XMAXIM.LT.SPEVEC(I,J)) THEN
                XMAXIM=SPEVEC(I,J)
            END IF
C
         END DO
C_______________________________________________________________________
C
C                                        Defining the bins
         XHISTO_MINIMS(J)=XMINIM
         XHISTO_MAXIMS(J)=XMAXIM
C
         BINSIZ(J)=(XMAXIM-XMINIM)/LDBINS         
C
         DO K=1,LDBINS+1
            XAUXIL_HISTOG(K,J)=XMINIM+(K-1)*BINSIZ(J) !histogram X-axis
         END DO
C_______________________________________________________________________
C
C                                        Begining to count
         DO K=1,LDBINS+1   
            DO I=1,LDMONT
C         
               IF (SPEVEC(I,J).LT.XAUXIL_HISTOG(K+1,J).AND.
     *             SPEVEC(I,J).GE.XAUXIL_HISTOG( K ,J)) THEN
C     
                   YAUXIL_HISTOG(K,J)=YAUXIL_HISTOG(K,J)+1.0
C                   
               END IF
C         
            END DO     
         END DO
C
      END DO  !J=1,LDCURV
C
C=======================================================================
C
C     Beautifying the histograms: arraging x-axis and normalizations
C
      Y_MAXI=-1.0E+10
      DO J=1,LDCURV
C
         K=1
         XHISTO(K,J)=XHISTO_MINIMS(J)-BINSIZ(J)
         YHISTO_UNNORM(K,J)=0.0000
C
         K=K+1
         XHISTO(K,J)=XHISTO_MINIMS(J)
         YHISTO_UNNORM(K,J)=0.0000
C
         DO I=1,LDBINS
C
            K=K+1
            XHISTO(K,J)=XHISTO_MINIMS(J)+(I-1)*BINSIZ(J)
            YHISTO_UNNORM(K,J)=YAUXIL_HISTOG(I,J)
C
            K=K+1
            XHISTO(K,J)=XHISTO_MINIMS(J)+I*BINSIZ(J)
            YHISTO_UNNORM(K,J)=YAUXIL_HISTOG(I,J)
C
         END DO
C
         K=K+1
         XHISTO(K,J)=XHISTO_MINIMS(J)+LDBINS*BINSIZ(J)
         YHISTO_UNNORM(K,J)=0.0000
C
         K=K+1
         XHISTO(K,J)=XHISTO_MINIMS(J)+(LDBINS+1)*BINSIZ(J)
         YHISTO_UNNORM(K,J)=0.0000
C
         LDHIST=K
C
         IF (LDHIST.GT.NDBINS) THEN
C
             WRITE(LOGFIL,'(/,''Actual No. of bins, LDHIST='',I5,1x,
     *                        ''exceeds the limit NDBINS='',I5,/)')
     *                                     LDHIST,NDBINS
C
             WRITE(000000,'(/,''Actual No. of bins, LDHIST='',I5,1x,
     *                        ''exceeds the limit NDBINS='',I5,/)')
     *                                     LDHIST,NDBINS
C
             STOP 'Actual number of bins out of allowed limit'
C
         END IF 
C_______________________________________________________________________
C
C                                  Normalization 1: Integral ->> Area
C
         XINTEG=0.0
         DO K=1,LDHIST,2
            XINTEG=XINTEG+YHISTO_UNNORM(K,J)*BINSIZ(J)
         END DO
C
         DO K=1,LDHIST
            YHISTO_INTEGR(K,J)=YHISTO_UNNORM(K,J)/XINTEG
         END DO
C_______________________________________________________________________
C
C                                  Normalization 2: Histog. maxim
C
         DO K=1,LDHIST
            IF (Y_MAXI.LT.YHISTO_UNNORM(K,J)) THEN
                Y_MAXI=YHISTO_UNNORM(K,J)
            END IF
         END DO
C
      END DO
C
      DO J=1,LDCURV
         DO K=1,LDHIST
            YHISTO_MAXIMS(K,J)=YHISTO_UNNORM(K,J)/Y_MAXI
         END DO
      END DO
C
C=======================================================================
C
      IF (LOGWRI.GT.0) THEN
          WRITE(LOGFIL,'(''Exiting HISTOG_FABRIC'')')
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
      SUBROUTINE HISTOG2FABRIC(N_READ,NDCURV,LDCURV,LDHIST,
     *                         XHISTO_MAXIMS,XHISTO_MINIMS,
     *                         BINSIZ,XHISTO,YHISTO_UNNORM,
     *                         YHISTO_INTEGR,YHISTO_MAXIMS)
C
      INCLUDE 'MATDIM/NDBINS.f'
C
      DIMENSION
     *          VECTOR(1:NDCURV)
      DIMENSION
     *          XAUXIL_HISTOG(1:NDBINS,1:NDCURV),
     *          YAUXIL_HISTOG(1:NDBINS,1:NDCURV)
      DIMENSION
     *          XHISTO(1:NDBINS,1:NDCURV),
     *          XHISTO_MINIMS(1:NDCURV),
     *          XHISTO_MAXIMS(1:NDCURV)
      DIMENSION
     *          YHISTO_UNNORM(1:NDBINS,1:NDCURV),
     *          YHISTO_INTEGR(1:NDBINS,1:NDCURV),
     *          YHISTO_MAXIMS(1:NDBINS,1:NDCURV)
      DIMENSION
     *          BINSIZ(1:NDCURV)
      COMMON
     *       /MONTEC/ IFPSEU,IFPARA,LDMONT,LDBINS
      COMMON
     *       /PRFILE/ ISCREN,LOGFIL,LOGCHP,LOGCHN,LOGHIP,LOGHIN,
     *                LOGPRO,LOGNEU,LOGWRI,LSCREN,IRESUL,ICONVE,
     *                NOUTPT,LOGAUX,LOGENE,LOGRAD,LOGGRD,LOGMSH,
     *                                                   LOGBIS
C
C=======================================================================
C
C     This subroutine creates the histograms reading the data from
C     temporary files. As input we have to give the maxima and the
C     minima so that the subroutine can compute the bin sizes.
C
C=======================================================================
C
      IF (LOGWRI.GT.0) THEN
          WRITE(LOGFIL,'(''Entering to HISTOG2FABRIC'')')
      END IF
C
C=======================================================================
C
C     Initializing the necessary quantities
C
      DO K=1,NDBINS
         DO J=1,LDCURV
            YHISTO_UNNORM(K,J)=0.0
            YHISTO_INTEGR(K,J)=0.0
            YHISTO_MAXIMS(K,J)=0.0
            YAUXIL_HISTOG(K,J)=0.0
         END DO
      END DO
C
C=======================================================================
C
C     Creating the size of the bins
C
      DO ICURVE=1,LDCURV
C
         BINSIZ(ICURVE)=(XHISTO_MAXIMS(ICURVE)-XHISTO_MINIMS(ICURVE))
     *                 /LDBINS
C
         DO K_BINS=1,LDBINS+1
            XAUXIL_HISTOG(K_BINS,ICURVE)=XHISTO_MINIMS(ICURVE)
     *                                  +(K_BINS-1)*BINSIZ(ICURVE)
         END DO
C
      END DO
C
C=======================================================================
C
C     Reading the file and counting
C
      DO IMONTE=1,LDMONT
C
         READ(N_READ,*)JDUMMY,(VECTOR(I),I=1,LDCURV)
C
         DO ICURVE=1,LDCURV
            DO K_BINS=1,LDBINS
C
               IF (VECTOR(ICURVE).LT.XAUXIL_HISTOG(K_BINS+1,ICURVE)
     *        .AND.VECTOR(ICURVE).GE.XAUXIL_HISTOG(K_BINS,ICURVE)) THEN
C
                   YAUXIL_HISTOG(K_BINS,ICURVE)
     *            =YAUXIL_HISTOG(K_BINS,ICURVE)+1.0
C
               END IF
C
            END DO
         END DO
C
      END DO
C
C=======================================================================
C
C     Beautifying the histograms(I): arraging x-axis
C
      DO ICURVE=1,LDCURV
C
         K=1
         XHISTO(K,ICURVE)=XHISTO_MINIMS(ICURVE)-BINSIZ(ICURVE)
         YHISTO_UNNORM(K,ICURVE)=0.0000
C
         K=K+1
         XHISTO(K,ICURVE)=XHISTO_MINIMS(ICURVE)
         YHISTO_UNNORM(K,ICURVE)=0.0000
C
         DO I=1,LDBINS
C
            K=K+1
            XHISTO(K,ICURVE)=XHISTO_MINIMS(ICURVE)+(I-1)*BINSIZ(ICURVE)
            YHISTO_UNNORM(K,ICURVE)=YAUXIL_HISTOG(I,ICURVE)
C
            K=K+1
            XHISTO(K,ICURVE)=XHISTO_MINIMS(ICURVE)+I*BINSIZ(ICURVE)
            YHISTO_UNNORM(K,ICURVE)=YAUXIL_HISTOG(I,ICURVE)
C
         END DO
C
         K=K+1
         XHISTO(K,ICURVE)=XHISTO_MINIMS(ICURVE)+LDBINS*BINSIZ(ICURVE)
         YHISTO_UNNORM(K,ICURVE)=0.0000
C
         K=K+1
         XHISTO(K,ICURVE)=XHISTO_MINIMS(ICURVE)+(LDBINS+1)
     *                                         *BINSIZ(ICURVE)
         YHISTO_UNNORM(K,ICURVE)=0.0000
C
         LDHIST=K
C_______________________________________________________________________
C
         IF (LDHIST.GT.NDBINS) THEN
C
             WRITE(LOGFIL,'(/,''Actual No. of bins, LDHIST='',I5,1x,
     *                        ''exceeds the limit NDBINS='',I5,/)')
     *                                     LDHIST,NDBINS
C
             WRITE(000000,'(/,''Actual No. of bins, LDHIST='',I5,1x,
     *                        ''exceeds the limit NDBINS='',I5,/)')
     *                                     LDHIST,NDBINS
C
             STOP
     *      'HISTOG2FRABRIC: Actual number of bins out of allowed limit'
C
         END IF
C
      END DO
C
C=======================================================================
C
C     Beautifying the histograms (II): normalizations
C
C                                  Normalization 1: Integral ->> Area
C
      DO ICURVE=1,LDCURV
C
         XINTEG=0.0
         DO K_BINS=1,LDHIST,2
            XINTEG=XINTEG+YHISTO_UNNORM(K_BINS,ICURVE)*BINSIZ(ICURVE)
         END DO
C
         DO K_BINS=1,LDHIST
            YHISTO_INTEGR(K_BINS,ICURVE)=YHISTO_UNNORM(K_BINS,ICURVE)
     *                                  /XINTEG
         END DO
C
      END DO
C
C                                  Normalization 2: Histog. maxim
C
      Y_MAXI=-1.0E+10
      DO ICURVE=1,LDCURV
         DO K_BINS=1,LDHIST
            IF (Y_MAXI.LT.YHISTO_UNNORM(K_BINS,ICURVE)) THEN
                Y_MAXI=YHISTO_UNNORM(K_BINS,ICURVE)
            END IF
         END DO
      END DO
C
      DO ICURVE=1,LDCURV
         DO K_BINS=1,LDHIST
            YHISTO_MAXIMS(K_BINS,ICURVE)=YHISTO_UNNORM(K_BINS,ICURVE)
     *                                  /Y_MAXI
         END DO
      END DO
C
C=======================================================================
C
      IF (LOGWRI.GT.0) THEN
          WRITE(LOGFIL,'(''Exiting HISTOG2FABRIC'')')
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
      SUBROUTINE HISTOG_PRINTI(LDCURV,LDHIST,XHISTO,YHISTO,
     *                         XHISTO_MINIMS,XHISTO_MAXIMS,YEXPEC,
     *                         HISTIT,HX_TIT,HY_TIT,LABELS,INUCLI,
     *                                              ICOLOR_HISTOG)
C
      INCLUDE 'MATDIM/NDMESH.f'
      INCLUDE 'MATDIM/NDSPEC.f'
      INCLUDE 'MATDIM/NDPARS.f'
      INCLUDE 'MATDIM/NDNUCL.f'
      INCLUDE 'MATDIM/NDCOLO.f'
      INCLUDE 'MATDIM/NDTITL.f'
      INCLUDE 'MATDIM/NDBINS.f'
C
      PARAMETER
     *         (NDMES2=NDMESH*NDMESH)
      CHARACTER
     *          NUCSYM*6,INPSYM*6,TYPCHI*6,VERSIO*3
      CHARACTER
     *          FILNAM*256,LABELS*100,NUCNAM*06,TITLES*12,
     *          ISONAM*002,TEXLAM*020,TITPAR*13,FITNUC*02,
     *          TITLES_LATEXS*050,NUCNAM_LATEXS*010
      CHARACTER
     *          HISTIT*256,HX_TIT*256,HY_TIT*256
      CHARACTER
     *          FILNAM_IFITED*14,FILNAM_NUCFIT*30,
     *          FITNUC_AUXILI*02,PARNAM*8
      CHARACTER
     *          AUXTX1*30,AUXTX2*30,AUXTX3*30,AUXTX4*30,
     *          AUXTX5*30,AUXTX6*30
C
      DIMENSION
     *          XHISTO(1:NDBINS,1:NDSPEC),
     *          YHISTO(1:NDBINS,1:NDSPEC),
     *          LABELS(1:NDBINS,1:NDSPEC)
      DIMENSION
     *          XHISTO_MINIMS(1:NDSPEC),
     *          XHISTO_MAXIMS(1:NDSPEC)
      DIMENSION
     *          YEXPEC(1:NDSPEC),
     *          YHISTO_MAXIMS(1:NDSPEC),
     *          XPOSIT_YMAXIM(1:NDSPEC)
      DIMENSION
     *          ICOLOR_HISTOG(1:NDSPEC)
      DIMENSION
     *          FITNUC(1:NDNUCL),
     *          FITNUC_AUXILI(1:NDNUCL)
      DIMENSION
     *          PARPOT_PRINTI(1:NDPARS,1:NDMES2),
     *          RMSVAL(1:NDSPEC,1:NDNUCL)      
      COMMON
     *       /NUCLEU/ IZ_FIX,IN_FIX,ISOSPI 
      COMMON
     *       /DENSTR/ IFDENS
      COMMON
     *       /TENSOR/ IFTENS
      COMMON
     *       /TENOPT/ ICENTT,ISORBT,ITENSR 
      COMMON
     *       /DEPPAI/ IFDEEP,IFPRON
      COMMON
     *       /IFPAIR/ IF_PAI,IFGPAI,IFDPAI
      COMMON 
     *       /CHITAK/ IF_SPE,IF_RAD,IF_GAP,IF_FER,IF_DEN,IF_RHO,IF_INV,
     *                LDCHI2,TYPCHI(1:NDNUCL),ITKCHI(1:NDNUCL)
      COMMON
     *       /IFLAMB/ IFPAR1,IFPAR2,IFPAR3,IFPAR4,
     *                IFPAR5,IFPAR6,IFPAR7,IFPAR8,
     *                IFPA09,IFPA10,IFPA11,IFPA12
      COMMON
     *       /PARCOR/ IFCORR,IFRCVC,IFRCAC,IFVCAC,IFRSVS,IFRSAS,IFVSAS
      COMMON
     *       /POTPOT/ PARPOT(1:NDPARS)
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
     *       /SEARCH/ IFTAKE(1:NDPARS)
      COMMON
     *       /VERSIN/ VERSIO
      COMMON
     *       /MONTEC/ IFPSEU,IFPARA,LDMONT,LDBINS
      COMMON
     *       /GAUSSI/ PARPOT_XMEANS(1:NDPARS),
     *                PARPOT_SIGMAS(1:NDPARS)
      COMMON
     *       /PARAMT_TITLES/ TITLES(1:NDTITL),
     *                       TITLES_LATEXS(1:NDTITL),
     *                       NUCNAM_LATEXS(1:NDNUCL)

      COMMON
     *       /PRFILE/ ISCREN,LOGFIL,LOGCHP,LOGCHN,LOGHIP,LOGHIN,
     *                LOGPRO,LOGNEU,LOGWRI,LSCREN,IRESUL,ICONVE,
     *                NOUTPT,LOGAUX,LOGENE,LOGRAD,LOGGRD,LOGMSH,
     *                                                   LOGBIS
      DATA
     *       FITNUC(01) / 'Ox' /,
     *       FITNUC(02) / 'C0' /,
     *       FITNUC(03) / 'C8' /,
     *       FITNUC(04) / 'Ni' /,
     *       FITNUC(05) / 'Zr' /,
     *       FITNUC(06) / 'Sn' /,
     *       FITNUC(07) / 'Gd' /,
     *       FITNUC(08) / 'Pb' /
C
C=======================================================================
C
      IF (LOGWRI.GT.0) WRITE(LOGFIL,'(/,''Entering HISTOG_PRINTI'')')
C
      IF (LOGWRI.GT.0) 
     *    WRITE(LOGFIL,'(/,'' IZ= '',I3,'' IN= '',I3,/)') 
     *                        IZ_FIX,IN_FIX
C
C=======================================================================
C
      NRESUL=60
C      
      IF (ISOSPI.EQ.1) ISONAM='-P'
      IF (ISOSPI.EQ.0) ISONAM='-N'
C
C=======================================================================
C 
      KACTIV=0
      DO JNUCLI=1,LDNUCL
         IF (ITAKNU(JNUCLI).EQ.1) THEN
             KACTIV=KACTIV+1
             FITNUC_AUXILI(KACTIV)=FITNUC(JNUCLI)
         END IF
      END DO
      
      I1=3*KACTIV
      
      WRITE(FILNAM_NUCFIT,'(''_'',<NUCACT>(A2,''-''))')
     *                      (FITNUC_AUXILI(I),I=1,NUCACT)
C
C=======================================================================
C
      WRITE(NUCNAM,'(A6)') NUCSYM(INUCLI)
      
      J1=1
      J2=6
      IF (NUCNAM(1:1).EQ.' ') J1=2
      IF (NUCNAM(2:2).EQ.' ') J1=3
      IF (NUCNAM(3:3).EQ.' ') J1=4
C
C=======================================================================
C
      IF (IFDENS.EQ.0) THEN
          WRITE(FILNAM,'(''SPEMontC/IFDENS-0/SPE-MC_'',A,A2,
     *                  ''_IFDENS-'',I1,''_IFTENS-'',I1,''_IF-PAI-'',I1,
     *                  ''_IF-RAD-'',I1,''_IF-ORD-'',I1,''_IFDEEP-'',I1,
     *                  ''_IFPRON-'',I1,
     *                  ''_IFCORR-'',I1,''_IFRCVC-'',I1,''_IFRCAC-'',I1,
     *                  ''_IFVCAC-'',I1,''_IFRSVS-'',I1,''_IFRSAS-'',I1,
     *                  ''_IFVSAS-'',I1,
     *                  ''_'',A3,''.dat'')')
     *
     *               NUCNAM(J1:J2),ISONAM,IFDENS,IFTENS,IF_PAI,
     *               IF_RAD,IF_INV,IFDEEP,IFPRON,IFCORR,IFRCVC,
     *               IFRCAC,IFVCAC,IFRSVS,IFRSAS,IFVSAS,VERSIO
      END IF
C
C=======================================================================
C
      IF (IFDENS.EQ.1 .AND. IFTENS.EQ.0) THEN
C
          WRITE(TEXLAM,'(''_LAMBDA'',SP,I2,I2,I2,I2)')
     *                            IFPAR1,IFPAR2,IFPAR3,IFPAR4
C
          WRITE(FILNAM,'(''SPEMontC/IFDENS-1_IFTENS-0/SPE-MC_'',A,A2,
     *                  ''_IFDENS-'',I1,''_IFTENS-'',I1,''_IF-PAI-'',I1,
     *                  ''_IF-RAD-'',I1,''_IF-ORD-'',I1,''_IFDEEP-'',I1,
     *                  ''_IFPRON-'',I1,
     *                  ''_IFCORR-'',I1,''_IFRCVC-'',I1,''_IFRCAC-'',I1,
     *                  ''_IFVCAC-'',I1,''_IFRSVS-'',I1,''_IFRSAS-'',I1,
     *                  ''_IFVSAS-'',I1,
     *                  A15,''_'',A3,''.dat'')')
     *
     *               NUCNAM(J1:J2),ISONAM,
     *               IFDENS,IFTENS,IF_PAI,IF_RAD,IF_INV,IFDEEP,IFPRON,
     *               IFCORR,IFRCVC,IFRCAC,IFVCAC,IFRSVS,IFRSAS,IFVSAS,
     *                                                  TEXLAM,VERSIO
      END IF
C
C=======================================================================
C
      IF (IFDENS.EQ.1 .AND. IFTENS.EQ.1) THEN
C
          WRITE(TEXLAM,'(''_LAMBDA'',SP,I2,I2,I2,I2)')
     *                            IFPAR1,IFPAR2,IFPAR3,IFPAR4
C
          WRITE(FILNAM,'(''SPEMontC/IFDENS-1_IFTENS-0/SPE-MC_'',A,A2,
     *                  ''_IFDENS-'',I1,''_IFTENS-'',I1,''_IF-PAI-'',I1,
     *                  ''_ICENTT-'',I1,''_ISORBT-'',I1,''_ITENSR-'',I1,
     *                  ''_IF-RAD-'',I1,''_IF-ORD-'',I1,''_IFDEEP-'',I1,
     *                  ''_IFPRON-'',I1,
     *                  ''_IFCORR-'',I1,''_IFRCVC-'',I1,''_IFRCAC-'',I1,
     *                  ''_IFVCAC-'',I1,''_IFRSVS-'',I1,''_IFRSAS-'',I1,
     *                  ''_IFVSAS-'',I1,
     *                  A15,''_'',A3,''.dat'')')
     *
     *               NUCNAM(J1:J2),ISONAM,
     *               IFDENS,IFTENS,IF_PAI,ICENTT,ISORBT,ITENSR,
     *               IF_RAD,IF_INV,IFDEEP,IFPRON,IFCORR,IFRCVC,
     *               IFRCAC,IFVCAC,IFRSVS,IFRSAS,IFVSAS,TEXLAM,VERSIO
      END IF
C
C=======================================================================
C
      OPEN(NRESUL,FILE=FILNAM,STATUS='UNKNOWN')
C
C=======================================================================
C
      XMINIM=+1.0E+10
      XMAXIM=-1.0E+10
C
      DO ICURVE=1,LDCURV
         IF (XHISTO_MINIMS(ICURVE).LT.XMINIM) THEN
             XMINIM=XHISTO_MINIMS(ICURVE)
         END IF
         IF (XHISTO_MAXIMS(ICURVE).GT.XMAXIM) THEN
             XMAXIM=XHISTO_MAXIMS(ICURVE)
         END IF
      END DO
C
      YMINIM=+1.0E+10
      YMAXIM=-1.0E+10
C
      DO ICURVE=1,LDCURV
         DO IPOINT=1,LDHIST
            IF (YHISTO(IPOINT,ICURVE).LT.YMINIM) THEN
                YMINIM=YHISTO(IPOINT,ICURVE)
            END IF
            IF (YHISTO(IPOINT,ICURVE).GT.YMAXIM) THEN
                YMAXIM=YHISTO(IPOINT,ICURVE)
            END IF
         END DO
      END DO
C
C=======================================================================
C
      DO ICURVE=1,LDCURV
         YHISTO_MAXIMS(ICURVE)=-1.0E+10
         DO IPOINT=1,LDHIST
            IF (YHISTO(IPOINT,ICURVE).GT.YHISTO_MAXIMS(ICURVE)) THEN
                YHISTO_MAXIMS(ICURVE)=YHISTO(IPOINT,ICURVE)
                XPOSIT_YMAXIM(ICURVE)=XHISTO(IPOINT,ICURVE)
            END IF
         END DO
      END DO
C
C=======================================================================
C
      WRITE(NRESUL,'(/,''<<CODEVERSIO>> VERSIO'')')
C
      WRITE(NRESUL,'(16X,A3)') VERSIO
C
C=======================================================================
C
      WRITE(NRESUL,'(/,''<<MAIN_TITLE>>'')')
C      
      WRITE(NRESUL,'(15X,A,/)') HISTIT
C
C=======================================================================
C     
      WRITE(NRESUL,'(''<<NUCLIDINFO>>'',1X,''NUMB_Z'',2X,''NUMB_N'')')
C
      WRITE(NRESUL,'(15X,I3,6X,I3,/)') IZ_FIX,IN_FIX
C
      A_MASS=IZ_FIX+IN_FIX
C
C=======================================================================
C
      WRITE(NRESUL,'(''<<NoNUCACTIV>>'',1X,''NUCACT'')')
C      
      WRITE(NRESUL,'(15X,I3,/)') NUCACT
C
C=======================================================================
C
      WRITE(NRESUL,'(''<<EXPER_FILE>>'',1X,''IFDEEP'',2X
     *                                     ''IFPRON'')')
      
      WRITE(NRESUL,'(15X,2(I3,6X),/)') IFDEEP,IFPRON
C
C=======================================================================
C
      WRITE(NRESUL,'(''<<POTEN_INFO>>'',1X,''ISOSPI'',2X,''IFDENS'',
     *                                  2X,''IFTENS'',2X,''IF_PAI'')')
C      
      WRITE(NRESUL,'(15X,4(I3,5X),/)') ISOSPI,IFDENS,IFTENS,IF_PAI
C
C=======================================================================
C
      WRITE(NRESUL,'(''<<TENSR_INFO>>'',1X,''ICENTT'',2X
     *                                     ''ISORBT'',2X,''ITENSR'')')
C      
      WRITE(NRESUL,'(15X,3(I3,5X),/)') ICENTT,ISORBT,ITENSR
C
C=======================================================================
C
      WRITE(NRESUL,'(''<<TAKEN_CHI2>> IF_SPE  IF_RAD  IF_GAP  '',
     *               ''IF_FER  IF_DEN  IF_RHO  IF_INV'')')
      WRITE(NRESUL,'(15X,7(I3,5X),/)')IF_SPE,IF_RAD,IF_GAP,IF_FER,
     *                                       IF_DEN,IF_RHO,IF_INV
C
C=======================================================================
C      
      WRITE(NRESUL,'(''<<XAXIS_TEXT>>'',/,15x,A,/)') HX_TIT
C      
      WRITE(NRESUL,'(''<<YAXIS_TEXT>>'',/,15x,A,/)') HY_TIT
C
C=======================================================================
C
      WRITE(NRESUL,'(''<<SIDE_TEXTS>>'')')
C
      AUXTX1=' '
      AUXTX2=' '
      AUXTX3=' '
      AUXTX4=' '
      AUXTX5=' '
      AUXTX6=' '
C
      IF (IFDENS.EQ.1) AUXTX1='${\\rm DENS}\\,\\vert\\,$'
      IF (IFTENS.EQ.1) AUXTX2='${\\rm TENS}\\,\\vert\\,$'
      IF (IF_PAI.EQ.1) AUXTX3='${\\rm PAIR}\\,\\vert\\,$'
      IF (IF_RAD.EQ.1) AUXTX4='${\\rm RADI}\\,\\vert\\,$'
      IF (IFDEEP.EQ.1) AUXTX5='${\\rm DEEP}\\,\\vert\\,$'
      IF (IFPRON.EQ.1) AUXTX6='${\\rm PRON}\\,\\vert\\,$'
C
      WRITE(NRESUL,'(15X,''~~~Id: $\\rm LDMONT= '',I6,
     *                           ''\\,\\vert\\,$'',6A)')
     *            LDMONT,AUXTX1,AUXTX2,AUXTX3,AUXTX4,AUXTX5,AUXTX6
C
      IF (ISOSPI.EQ.1) THEN
          WRITE(NRESUL,'(15X,''~~~$\\sigma_{V_c}= '',f6.4,
     *                       ''\\,\\vert\\,$'',
     *                       '' $\\sigma_{r_c}= '',f6.4,
     *                       ''\\,\\vert\\,$'',
     *                       '' $\\sigma_{a_c}= '',f6.4,''$'')')
     *       PARPOT_SIGMAS(1),PARPOT_SIGMAS(2),PARPOT_SIGMAS(3)
          WRITE(NRESUL,'(15X,''~~~$\\sigma_{\\lambda_{so}}= '',f6.4,
     *                       ''\\,\\vert\\,$'',
     *                       '' $\\sigma_{r_{so}}= '',f6.4,
     *                       ''\\,\\vert\\,$'',
     *                       '' $\\sigma_{a_{so}}= '',f6.4,''$'')')
     *       PARPOT_SIGMAS(4),PARPOT_SIGMAS(5),PARPOT_SIGMAS(6)
      END IF
C
      IF (ISOSPI.EQ.0) THEN
          WRITE(NRESUL,'(15X,''~~~$\\sigma_{V_c}= '',f6.4,
     *                       ''\\,\\vert\\,$'',
     *                       '' $\\sigma_{r_c}= '',f6.4,
     *                       ''\\,\\vert\\,$'',
     *                       '' $\\sigma_{a_c}= '',f6.4,''$'')')
     *       PARPOT_SIGMAS(21),PARPOT_SIGMAS(22),PARPOT_SIGMAS(23)
          WRITE(NRESUL,'(15X,''~~~$\\sigma_{\\lambda_{so}}= '',f6.4,
     *                       ''\\,\\vert\\,$'',
     *                       '' $\\sigma_{r_{so}}= '',f6.4,
     *                       ''\\,\\vert\\,$'',
     *                       '' $\\sigma_{a_{so}}= '',f6.4,''$'')')
     *       PARPOT_SIGMAS(24),PARPOT_SIGMAS(25),PARPOT_SIGMAS(26)
      END IF
C
C=======================================================================
C     
      WRITE(NRESUL,'(/,15X,30(''=''),/)')
          
      WRITE(NRESUL,'(''<<X_AX_PARAM>> XMINIM_FIG  XMAXIM_FIG  '',
     *                   ''X_STEP_FIG'')')
          
      WRITE(NRESUL,'(15X,3(F10.2,2X),/)')XMINIM,XMAXIM,(XMAXIM-XMINIM)
     
      WRITE(NRESUL,'(''<<Y_AX_PARAM>> YMINIM_FIG  YMAXIM_FIG  '',
     *                   ''Y_STEP_FIG'')')
          
      WRITE(NRESUL,'(15X,3(F10.2,2X),/)')
     *                         YMINIM,YMAXIM,(YMAXIM-YMINIM)
     
      WRITE(NRESUL,'(15X,30(''=''),/)')
C
C=======================================================================
C
      I_COLM=0
      I_LEFT=0
      I_RIGH=0
C
      WRITE(NRESUL,'(''<<WHATLEGEND>>'',1X,''COLM_LABEL'',2X,
     *               ''LEFT_LABEL'',2X,''RIGH_LABEL'')')
      WRITE(NRESUL,'(15X,I5,7X,I5,7X,I5)')I_COLM,I_LEFT,I_RIGH
C
C=======================================================================
C      
      IFBULL=0
      IFHIST=1
      IFOVER=0
C      
      WRITE(NRESUL,'(/,''<<NO_OF_CURV>>'',1X,''NUMB_CURVE'',2X,
     *                 ''IF_BULLETS'',2X,
     *                 ''IF_HISTOGS'',2X,''IF_OVERFITING'')')
      WRITE(NRESUL,'(15X,I5,7X,I5,7X,I5,7X,I5,7X,I5)')LDCURV,
     *                                         IFBULL,IFHIST,IFOVER
C
C=======================================================================
C
      ISTYLE=0
      I_TYPE=3
      ITHICK=3
C      ICOLOR=70
C
      DO I_CURV=1,LDCURV
C
         WRITE(NRESUL,'(''<<>>'')')
         WRITE(NRESUL,'(''<<CURVE_'',I4.4,''>>'',1X,
     *                  ''THICK_LINE'',2X,''STYLE_LINE'',2X,
     *                  ''COLOR_LINE'',2X,''TYPE_POINT'')')I_CURV
         WRITE(NRESUL,'(15X,I5,7X,I5,7X,I5,7X,I5)')ITHICK,ISTYLE,
     *                                      ICOLOR_HISTOG(I_CURV),
     *                                                    I_TYPE
         WRITE(NRESUL,'(15X,''NUMB_POINT'',2X,''YXPECTED_V'',
     *                   2X,''XPOSI_YMAX'')')
         WRITE(NRESUL,'(15X,I10,2X,F10.4,2X,F10.4)')
     *                      LDHIST,YEXPEC(I_CURV),XPOSIT_YMAXIM(I_CURV)
C
         DO IPOINT=1,LDHIST
C
            WRITE(NRESUL,'(15X,F10.4,4X,E12.4,4X,A100)') 
     *                           XHISTO(IPOINT,I_CURV),
     *                           YHISTO(IPOINT,I_CURV),
     *                           LABELS(IPOINT,I_CURV)
C
         END DO
C
      END DO
C
C=======================================================================
C      
      WRITE(NRESUL,'(/,''<<GO_GETTHEM>>'',/)')
C
C=======================================================================
C
      RETURN
      END
C
C=======================================================================
C=======================================================================
C
      SUBROUTINE HISTOG2PRINTI(NDCURV,LDCURV,LDHIST,
     *                         XHISTO,YHISTO,LABELS,
     *                         YEXPEC,HISTIT,HX_TIT,
     *                         HY_TIT,ICOLOR_HISTOG,
     *                                       NRESUL)
C
      INCLUDE 'MATDIM/NDNUCL.f'
      INCLUDE 'MATDIM/NDPARS.f'
      INCLUDE 'MATDIM/NDCOLO.f'
      INCLUDE 'MATDIM/NDTITL.f'
      INCLUDE 'MATDIM/NDBINS.f'
C
      CHARACTER
     *          NUCSYM*6,INPSYM*6,TYPCHI*6,VERSIO*3
      CHARACTER
     *          FILNAM*256,LABELS*100,NUCNAM*05,TITLES*12,
     *          ISONAM*002,TEXLAM*020,TITPAR*13,FITNUC*02,
     *          TITLES_LATEXS*050,NUCNAM_LATEXS*010
      CHARACTER
     *          HISTIT*256,HX_TIT*256,HY_TIT*256
      CHARACTER
     *          FILNAM_IFITED*14,FILNAM_NUCFIT*30,
     *          FITNUC_AUXILI*02,PARNAM*8
      CHARACTER
     *          AUXTX1*30,AUXTX2*30,AUXTX3*30,AUXTX4*30,
     *          AUXTX5*30,AUXTX6*30
C
      DIMENSION
     *          XHISTO(1:NDBINS,1:NDCURV),
     *          YHISTO(1:NDBINS,1:NDCURV),
     *          LABELS(1:NDBINS,1:NDCURV)
      DIMENSION
     *          XHISTO_MINIMS(1:NDCURV),
     *          XHISTO_MAXIMS(1:NDCURV)
      DIMENSION
     *          YEXPEC(1:NDCURV),
     *          YHISTO_MAXIMS(1:NDCURV),
     *          XPOSIT_YMAXIM(1:NDCURV)
      DIMENSION
     *          ICOLOR_HISTOG(1:NDCURV)
      COMMON
     *       /NUCLEU/ IZ_FIX,IN_FIX,ISOSPI 
      COMMON
     *       /DENSTR/ IFDENS
      COMMON
     *       /TENSOR/ IFTENS
      COMMON
     *       /TENOPT/ ICENTT,ISORBT,ITENSR 
      COMMON
     *       /DEPPAI/ IFDEEP,IFPRON
      COMMON
     *       /IFPAIR/ IF_PAI,IFGPAI,IFDPAI
      COMMON 
     *       /CHITAK/ IF_SPE,IF_RAD,IF_GAP,IF_FER,IF_DEN,IF_RHO,IF_INV,
     *                LDCHI2,TYPCHI(1:NDNUCL),ITKCHI(1:NDNUCL)
      COMMON
     *       /IFLAMB/ IFPAR1,IFPAR2,IFPAR3,IFPAR4,
     *                IFPAR5,IFPAR6,IFPAR7,IFPAR8,
     *                IFPA09,IFPA10,IFPA11,IFPA12
      COMMON
     *       /PARCOR/ IFCORR,IFRCVC,IFRCAC,IFVCAC,IFRSVS,IFRSAS,IFVSAS
      COMMON
     *       /POTPOT/ PARPOT(1:NDPARS)
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
     *       /SEARCH/ IFTAKE(1:NDPARS)
      COMMON
     *       /VERSIN/ VERSIO
      COMMON
     *       /MONTEC/ IFPSEU,IFPARA,LDMONT,LDBINS
      COMMON
     *       /GAUSSI/ PARPOT_XMEANS(1:NDPARS),
     *                PARPOT_SIGMAS(1:NDPARS)
      COMMON
     *       /PARAMT_TITLES/ TITLES(1:NDTITL),
     *                       TITLES_LATEXS(1:NDTITL),
     *                       NUCNAM_LATEXS(1:NDNUCL)

      COMMON
     *       /PRFILE/ ISCREN,LOGFIL,LOGCHP,LOGCHN,LOGHIP,LOGHIN,
     *                LOGPRO,LOGNEU,LOGWRI,LSCREN,IRESUL,ICONVE,
     *                NOUTPT,LOGAUX,LOGENE,LOGRAD,LOGGRD,LOGMSH,
     *                                                   LOGBIS
C
C=======================================================================
C
      IF (LOGWRI.GT.0) WRITE(LOGFIL,'(/,''Entering HISTOG2PRINTI'')')
C
C=======================================================================
C
      XMINIM=+1.0E+10
      XMAXIM=-1.0E+10
C
      DO ICURVE=1,LDCURV
         DO IPOINT=1,LDHIST
             IF (XHISTO(IPOINT,ICURVE).LT.XMINIM) THEN
                 XMINIM=XHISTO(IPOINT,ICURVE)
             END IF
             IF (XHISTO(IPOINT,ICURVE).GT.XMAXIM) THEN
                 XMAXIM=XHISTO(IPOINT,ICURVE)
             END IF
          END DO
      END DO
C
      YMINIM=+1.0E+10
      YMAXIM=-1.0E+10
C
      DO ICURVE=1,LDCURV
         DO IPOINT=1,LDHIST
            IF (YHISTO(IPOINT,ICURVE).LT.YMINIM) THEN
                YMINIM=YHISTO(IPOINT,ICURVE)
            END IF
            IF (YHISTO(IPOINT,ICURVE).GT.YMAXIM) THEN
                YMAXIM=YHISTO(IPOINT,ICURVE)
            END IF
         END DO
      END DO
C
C=======================================================================
C
      DO ICURVE=1,LDCURV
         YHISTO_MAXIMS(ICURVE)=-1.0E+10
         DO IPOINT=1,LDHIST
            IF (YHISTO(IPOINT,ICURVE).GT.YHISTO_MAXIMS(ICURVE)) THEN
                YHISTO_MAXIMS(ICURVE)=YHISTO(IPOINT,ICURVE)
                XPOSIT_YMAXIM(ICURVE)=XHISTO(IPOINT,ICURVE)
            END IF
         END DO
      END DO
C
C=======================================================================
C
      WRITE(NRESUL,'(/,''<<CODEVERSIO>> VERSIO'')')
C
      WRITE(NRESUL,'(16X,A3)') VERSIO
C
C=======================================================================
C
      WRITE(NRESUL,'(/,''<<MAIN_TITLE>>'')')
C      
      WRITE(NRESUL,'(15X,A,/)') HISTIT
C
C=======================================================================
C     
      WRITE(NRESUL,'(''<<NUCLIDINFO>>'',1X,''NUMB_Z'',2X,''NUMB_N'')')
C
      WRITE(NRESUL,'(15X,I3,6X,I3,/)') IZ_FIX,IN_FIX
C
      A_MASS=IZ_FIX+IN_FIX
C
C=======================================================================
C
      WRITE(NRESUL,'(''<<NoNUCACTIV>>'',1X,''NUCACT'')')
C      
      WRITE(NRESUL,'(15X,I3,/)') NUCACT
C
C=======================================================================
C
      WRITE(NRESUL,'(''<<EXPER_FILE>>'',1X,''IFDEEP'',2X
     *                                     ''IFPRON'')')
      
      WRITE(NRESUL,'(15X,2(I3,6X),/)') IFDEEP,IFPRON
C
C=======================================================================
C
      WRITE(NRESUL,'(''<<POTEN_INFO>>'',1X,''ISOSPI'',2X,''IFDENS'',
     *                                  2X,''IFTENS'',2X,''IF_PAI'')')
C      
      WRITE(NRESUL,'(15X,4(I3,5X),/)') ISOSPI,IFDENS,IFTENS,IF_PAI
C
C=======================================================================
C
      WRITE(NRESUL,'(''<<TENSR_INFO>>'',1X,''ICENTT'',2X
     *                                     ''ISORBT'',2X,''ITENSR'')')
C      
      WRITE(NRESUL,'(15X,3(I3,5X),/)') ICENTT,ISORBT,ITENSR
C
C=======================================================================
C
      WRITE(NRESUL,'(''<<TAKEN_CHI2>> IF_SPE  IF_RAD  IF_GAP  '',
     *               ''IF_FER  IF_DEN  IF_RHO  IF_INV'')')
      WRITE(NRESUL,'(15X,7(I3,5X),/)')IF_SPE,IF_RAD,IF_GAP,IF_FER,
     *                                       IF_DEN,IF_RHO,IF_INV
C
C=======================================================================
C      
      WRITE(NRESUL,'(''<<XAXIS_TEXT>>'',/,15x,A,/)') HX_TIT
C      
      WRITE(NRESUL,'(''<<YAXIS_TEXT>>'',/,15x,A,/)') HY_TIT
C
C=======================================================================
C
      WRITE(NRESUL,'(''<<SIDE_TEXTS>>'')')
C
      AUXTX1=' '
      AUXTX2=' '
      AUXTX3=' '
      AUXTX4=' '
      AUXTX5=' '
      AUXTX6=' '
C
      IF (IFDENS.EQ.1) AUXTX1='${\\rm DENS}\\,\\vert\\,$'
      IF (IFTENS.EQ.1) AUXTX2='${\\rm TENS}\\,\\vert\\,$'
      IF (IF_PAI.EQ.1) AUXTX3='${\\rm PAIR}\\,\\vert\\,$'
      IF (IF_RAD.EQ.1) AUXTX4='${\\rm RADI}\\,\\vert\\,$'
      IF (IFDEEP.EQ.1) AUXTX5='${\\rm DEEP}\\,\\vert\\,$'
      IF (IFPRON.EQ.1) AUXTX6='${\\rm PRON}\\,\\vert\\,$'
C
      WRITE(NRESUL,'(15X,''~~~Id: $\\rm LDMONT= '',I6,
     *                           ''\\,\\vert\\,$'',6A)')
     *            LDMONT,AUXTX1,AUXTX2,AUXTX3,AUXTX4,AUXTX5,AUXTX6
C
      WRITE(NRESUL,'(15X,''~~~Fitted: '',$)')
      DO JNUCLI=1,LDNUCL
         IF (ITAKNU(JNUCLI).EQ.1) THEN
             WRITE(NRESUL,'(2X,A10,$)') NUCNAM_LATEXS(JNUCLI)
         END IF
      END DO
      WRITE(NRESUL,'()')
C
C=======================================================================
C     
      WRITE(NRESUL,'(/,15X,30(''=''),/)')
          
      WRITE(NRESUL,'(''<<X_AX_PARAM>> XMINIM_FIG  XMAXIM_FIG  '',
     *                   ''X_STEP_FIG'')')
          
      WRITE(NRESUL,'(15X,3(F10.2,2X),/)')XMINIM,XMAXIM,(XMAXIM-XMINIM)
     
      WRITE(NRESUL,'(''<<Y_AX_PARAM>> YMINIM_FIG  YMAXIM_FIG  '',
     *                   ''Y_STEP_FIG'')')
          
      WRITE(NRESUL,'(15X,3(F10.2,2X),/)')
     *                         YMINIM,YMAXIM,(YMAXIM-YMINIM)
     
      WRITE(NRESUL,'(15X,30(''=''),/)')
C
C=======================================================================
C
      I_COLM=0
      I_LEFT=0
      I_RIGH=0
C
      WRITE(NRESUL,'(''<<WHATLEGEND>>'',1X,''COLM_LABEL'',2X,
     *               ''LEFT_LABEL'',2X,''RIGH_LABEL'')')
      WRITE(NRESUL,'(15X,I5,7X,I5,7X,I5)')I_COLM,I_LEFT,I_RIGH
C
C=======================================================================
C      
      IFBULL=0
      IFHIST=1
      IFOVER=0
C      
      WRITE(NRESUL,'(/,''<<NO_OF_CURV>>'',1X,''NUMB_CURVE'',2X,
     *                 ''IF_BULLETS'',2X,
     *                 ''IF_HISTOGS'',2X,''IF_OVERFITING'')')
      WRITE(NRESUL,'(15X,I5,7X,I5,7X,I5,7X,I5,7X,I5)')LDCURV,
     *                                         IFBULL,IFHIST,IFOVER
C
C=======================================================================
C
      ISTYLE=0
      I_TYPE=3
      ITHICK=3
C      ICOLOR=70
C
      DO I_CURV=1,LDCURV
C
         WRITE(NRESUL,'(''<<>>'')')
         WRITE(NRESUL,'(''<<CURVE_'',I4.4,''>>'',1X,
     *                  ''THICK_LINE'',2X,''STYLE_LINE'',2X,
     *                  ''COLOR_LINE'',2X,''TYPE_POINT'')')I_CURV
         WRITE(NRESUL,'(15X,I5,7X,I5,7X,I5,7X,I5)')ITHICK,ISTYLE,
     *                                      ICOLOR_HISTOG(I_CURV),
     *                                                    I_TYPE
         WRITE(NRESUL,'(15X,''NUMB_POINT'',2X,''YXPECTED_V'',
     *                   2X,''XPOSI_YMAX'')')
         WRITE(NRESUL,'(15X,I10,2X,F10.4,2X,F10.4)')
     *                      LDHIST,YEXPEC(I_CURV),XPOSIT_YMAXIM(I_CURV)
C
         DO IPOINT=1,LDHIST
C
            WRITE(NRESUL,'(15X,F10.4,4X,E12.4,4X,A100)') 
     *                           XHISTO(IPOINT,I_CURV),
     *                           YHISTO(IPOINT,I_CURV)/YMAXIM,
     *                           LABELS(IPOINT,I_CURV)
C
         END DO

C
      END DO
C
C=======================================================================
C      
      WRITE(NRESUL,'(/,''<<GO_GETTHEM>>'',/)')
C
C=======================================================================
C
      RETURN
      END
C
C=======================================================================
C=======================================================================
C
      SUBROUTINE PLOTINFOPRINT(NRESUL,HISTIT,HX_TIT,HY_TIT,
     *                         LEVACT_PROTON,LEVACT_NEUTRS,
     *                         SIGACT_PROTON,SIGACT_NEUTRS,
     *                                LDCURV,YEXPEC,XUNITS)
C
      INCLUDE 'MATDIM/NDNUCL.f'
      INCLUDE 'MATDIM/NDPARS.f'
      INCLUDE 'MATDIM/NDCOLO.f'
      INCLUDE 'MATDIM/NDTITL.f'
      INCLUDE 'MATDIM/NDBINS.f'
C
      CHARACTER
     *          NUCSYM*6,INPSYM*6,TYPCHI*6,VERSIO*3
      CHARACTER
     *          FILNAM*256,LABELS*100,NUCNAM*05,TITLES*12,
     *          ISONAM*002,TEXLAM*020,TITPAR*13,FITNUC*02,
     *          TITLES_LATEXS*050,NUCNAM_LATEXS*010,XUNITS*40
      CHARACTER
     *          HISTIT*256,HX_TIT*256,HY_TIT*256
      CHARACTER
     *          FILNAM_IFITED*14,FILNAM_NUCFIT*30,
     *          FITNUC_AUXILI*02,PARNAM*8
      CHARACTER
     *          AUXTX1*30,AUXTX2*30,AUXTX3*30,AUXTX4*30,
     *          AUXTX5*30,AUXTX6*30
      COMMON
     *       /NUCLEU/ IZ_FIX,IN_FIX,ISOSPI
      COMMON
     *       /DENSTR/ IFDENS
      COMMON
     *       /TENSOR/ IFTENS
      COMMON
     *       /TENOPT/ ICENTT,ISORBT,ITENSR
      COMMON
     *       /DEPPAI/ IFDEEP,IFPRON
      COMMON
     *       /IFPAIR/ IF_PAI,IFGPAI,IFDPAI
      COMMON
     *       /CHITAK/ IF_SPE,IF_RAD,IF_GAP,IF_FER,IF_DEN,IF_RHO,IF_INV,
     *                LDCHI2,TYPCHI(1:NDNUCL),ITKCHI(1:NDNUCL)
      COMMON
     *       /IFLAMB/ IFPAR1,IFPAR2,IFPAR3,IFPAR4,
     *                IFPAR5,IFPAR6,IFPAR7,IFPAR8,
     *                IFPA09,IFPA10,IFPA11,IFPA12
      COMMON
     *       /PARCOR/ IFCORR,IFRCVC,IFRCAC,IFVCAC,IFRSVS,IFRSAS,IFVSAS
      COMMON
     *       /POTPOT/ PARPOT(1:NDPARS)
      COMMON
     *       /VLIMIT/ VMISTR(1:NDPARS),
     *                VMIMIN(1:NDPARS),
     *                VMIMAX(1:NDPARS),
     *                VMISTP(1:NDPARS)
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
     *       /SEARCH/ IFTAKE(1:NDPARS)
      COMMON
     *       /FACCOU/ COUFAC
      COMMON
     *       /VERSIN/ VERSIO
      COMMON
     *       /MONTEC/ IFPSEU,IFPARA,LDMONT,LDBINS
      COMMON
     *       /GAUSSI/ PARPOT_XMEANS(1:NDPARS),
     *                PARPOT_SIGMAS(1:NDPARS)
      COMMON
     *       /PARAMT_TITLES/ TITLES(1:NDTITL),
     *                       TITLES_LATEXS(1:NDTITL),
     *                       NUCNAM_LATEXS(1:NDNUCL)
      COMMON
     *       /PRFILE/ ISCREN,LOGFIL,LOGCHP,LOGCHN,LOGHIP,LOGHIN,
     *                LOGPRO,LOGNEU,LOGWRI,LSCREN,IRESUL,ICONVE,
     *                NOUTPT,LOGAUX,LOGENE,LOGRAD,LOGGRD,LOGMSH,
     *                                                   LOGBIS
C
C=======================================================================
C
      WRITE(NRESUL,'(/,''<<CODEVERSIO>> VERSIO'')')
C
      WRITE(NRESUL,'(16X,A3)') VERSIO
C
C=======================================================================
C
      WRITE(NRESUL,'(/,''<<MAIN_TITLE>>'')')
C
      WRITE(NRESUL,'(15X,A,/)') HISTIT
C
C=======================================================================
C
      WRITE(NRESUL,'(''<<NUCLIDINFO>>'',1X,''NUMB_Z'',2X,''NUMB_N'')')
C
      WRITE(NRESUL,'(15X,I3,6X,I3,/)') IZ_FIX,IN_FIX
C
      A_MASS=IZ_FIX+IN_FIX
C
C=======================================================================
C
      WRITE(NRESUL,'(''<<NoNUCACTIV>>'',1X,''NUCACT'')')
C
      WRITE(NRESUL,'(15X,I3,/)') NUCACT
C
C=======================================================================
C
      WRITE(NRESUL,'(''<<EXPER_FILE>>'',1X,''IFDEEP'',2X
     *                                     ''IFPRON'')')

      WRITE(NRESUL,'(15X,2(I3,6X),/)') IFDEEP,IFPRON
C
C=======================================================================
C
      WRITE(NRESUL,'(''<<POTEN_INFO>>'',1X,''ISOSPI'',2X,''IFDENS'',
     *                                  2X,''IFTENS'',2X,''IF_PAI'')')
C
      WRITE(NRESUL,'(15X,4(I3,5X),/)') ISOSPI,IFDENS,IFTENS,IF_PAI
C
C=======================================================================
C
      WRITE(NRESUL,'(''<<TENSR_INFO>>'',1X,''ICENTT'',2X
     *                                     ''ISORBT'',2X,''ITENSR'')')
C
      WRITE(NRESUL,'(15X,3(I3,5X),/)') ICENTT,ISORBT,ITENSR
C
C=======================================================================
C
      WRITE(NRESUL,'(''<<TAKEN_CHI2>> IF_SPE  IF_RAD  IF_GAP  '',
     *               ''IF_FER  IF_DEN  IF_RHO  IF_INV'')')
      WRITE(NRESUL,'(15X,7(I3,5X),/)')IF_SPE,IF_RAD,IF_GAP,IF_FER,
     *                                       IF_DEN,IF_RHO,IF_INV
C
C=======================================================================
C
      IF (ISOSPI.EQ.1) THEN
          WRITE(NRESUL,'(''<<COULOMBFAC>> COUFAC'')')
          WRITE(NRESUL,'(15X,F8.2,/)')COUFAC
      END IF
C
C=======================================================================
C
      WRITE(NRESUL,'(''<<XAXIS_TEXT>>'',/,15x,A,/)') HX_TIT
C
      WRITE(NRESUL,'(''<<YAXIS_TEXT>>'',/,15x,A,/)') HY_TIT
C
C=======================================================================
C
      WRITE(NRESUL,'(''<<SIDE_TEXTS>>'')')
C
      AUXTX1=' '
      AUXTX2=' '
      AUXTX3=' '
      AUXTX4=' '
      AUXTX5=' '
      AUXTX6=' '
C
      IF (IFDENS.EQ.1) AUXTX1='{\\rm DENS}\\,\\vert\\,'
      IF (IFTENS.EQ.1) AUXTX2='{\\rm TENS}\\,\\vert\\,'
      IF (IF_PAI.EQ.1) AUXTX3='{\\rm PAIR}\\,\\vert\\,'
      IF (IF_RAD.EQ.1) AUXTX4='{\\rm RADI}\\,\\vert\\,'
      IF (IFDEEP.EQ.1) AUXTX5='{\\rm DEEP}\\,\\vert\\,'
      IF (IFPRON.EQ.1) AUXTX6='{\\rm PRON}\\,\\vert\\,'
C
      IF (ISOSPI.EQ.1) THEN
C
          WRITE(NRESUL,'(15X,''~~~Id: $ N_{MC}= '',I6,
     *                       '',\\,N_p= '',i3,'',\\,\\sigma_p= '',f5.3,
     *                       '',\\,f_c= '',f6.2,
     *                       ''\\,\\vert\\,'',6A,''$'')')
     *            LDMONT,LEVACT_PROTON,SIGACT_PROTON,COUFAC,
     *            AUXTX1,AUXTX2,AUXTX3,AUXTX4,AUXTX5,AUXTX6
C
          WRITE(NRESUL,'(15X,''~~~Fitted: '',$)')
          DO JNUCLI=1,LDNUCL
             IF (ITAKNU(JNUCLI).EQ.1) THEN
                 WRITE(NRESUL,'(2X,A10,$)') NUCNAM_LATEXS(JNUCLI)
             END IF
          END DO
          WRITE(NRESUL,'()')
C
          IF (IFDENS.EQ.0) THEN
C
              WRITE(NRESUL,'(15X,''~~~$ V_p^c= '',f8.2,
     *                                           ''\\,{\\rm MeV},'',
     *                             ''\\,r_p^c= '',f8.2,
     *                                           ''\\,{\\rm fm},'',
     *                             ''\\,a_p^c= '',f8.2,
     *                                           ''\\,{\\rm fm}$'')')
     *                                  VMISTR(1),VMISTR(2),VMISTR(3)
              WRITE(NRESUL,'(15X,''~~~$ \\lambda_p^{so}= '',f8.2,
     *                           ''\\,{\\rm MeV\\,fm^2/\\hbar^2},'',
     *                           ''\\,r_p^{so}= '',f8.2,
     *                           ''\\,{\\rm fm},'',
     *                           ''\\,a_p^{so}= '',f8.2,
     *                           ''\\,{\\rm fm}$'')')
     *                                 VMISTR(4),VMISTR(5),VMISTR(6)
C
          END IF
C
          IF (IFDENS.EQ.1) THEN
C
              WRITE(NRESUL,'(15X,''~~~$ V_p^c= '',f8.2,
     *                                           ''\\,{\\rm MeV},'',
     *                             ''\\,r_p^c= '',f8.2,
     *                                           ''\\,{\\rm fm},'',
     *                             ''\\,a_p^c= '',f8.2,
     *                                           ''\\,{\\rm fm}$'')')
     *                                  VMISTR(1),VMISTR(2),VMISTR(3)
              WRITE(NRESUL,'(15X,''~~~$ \\lambda^{\\pi\\pi}= '',f10.2,
     *                           ''\\,{\\rm MeV\\,fm^5/\\hbar^2},'',
     *                           ''\\,\\lambda^{\\pi\\nu}= '',f10.2,
     *                           ''\\,{\\rm MeV\\,fm^5/\\hbar^2}$'')')
     *                                 VMISTR(39),VMISTR(40)
C
          END IF
C
      END IF
C
      IF (ISOSPI.EQ.0) THEN
C
          WRITE(NRESUL,'(15X,''~~~Id: $ N_{MC}= '',I6,
     *                       '',\\,N_n= '',i3,'',\\,\\sigma_n= '',f5.3,
     *                           ''\\,\\vert\\,'',6A,''$'')')
     *            LDMONT,LEVACT_NEUTRS,SIGACT_NEUTRS,
     *            AUXTX1,AUXTX2,AUXTX3,AUXTX4,AUXTX5,AUXTX6
C
          WRITE(NRESUL,'(15X,''~~~Fitted: '',$)')
          DO JNUCLI=1,LDNUCL
             IF (ITAKNU(JNUCLI).EQ.1) THEN
                 WRITE(NRESUL,'(2X,A10,$)') NUCNAM_LATEXS(JNUCLI)
             END IF
          END DO
          WRITE(NRESUL,'()')
C
          IF (IFDENS.EQ.0) THEN
C
              WRITE(NRESUL,'(15X,''~~~$ V_n^c= '',f8.2,
     *                                           ''\\,{\\rm MeV},'',
     *                             ''\\,r_n^c= '',f8.2,
     *                                           ''\\,{\\rm fm},'',
     *                             ''\\,a_n^c= '',f8.2,
     *                                           ''\\,{\\rm fm}$'')')
     *                                VMISTR(21),VMISTR(22),VMISTR(23)
              WRITE(NRESUL,'(15X,''~~~$ \\lambda_n^{so}= '',f8.2,
     *                           ''\\,{\\rm MeV\\,fm^2/\\hbar^2},'',
     *                           ''\\,r_n^{so}= '',f8.2,
     *                           ''\\,{\\rm fm},'',
     *                           ''\\,a_n^{so}= '',f8.2,
     *                           ''\\,{\\rm fm}$'')')
     *                                 VMISTR(24),VMISTR(25),VMISTR(26)
C
          END IF
C
          IF (IFDENS.EQ.1) THEN
C
              WRITE(NRESUL,'(15X,''~~~$ V_n^c= '',f8.2,
     *                                           ''\\,{\\rm MeV},'',
     *                             ''\\,r_n^c= '',f8.2,
     *                                           ''\\,{\\rm fm},'',
     *                             ''\\,a_n^c= '',f8.2,
     *                                           ''\\,{\\rm fm}$'')')
     *                                 VMISTR(21),VMISTR(22),VMISTR(23)
              WRITE(NRESUL,'(15X,''~~~$ \\lambda^{\\nu\\nu}= '',f10.2,
     *                           ''\\,{\\rm MeV\\,fm^5/\\hbar^2},'',
     *                           ''\\,\\lambda^{\\nu\\pi}= '',f10.2,
     *                           ''\\,{\\rm MeV\\,fm^5/\\hbar^2}$'')')
     *                                 VMISTR(42),VMISTR(41)
C
          END IF
C
      END IF
C
C=======================================================================
C
      WRITE(NRESUL,'(''<<VALUE_UNIT>>'')')
      WRITE(NRESUL,'(15X,A)') XUNITS
C
C=======================================================================
C
      IFBULL=0
      IFHIST=1
      IFOVER=0
C
      ITHICK=4
      ISTYLE=0
      ICOLOR=44
      I_TYPE=3
C
      XPOSIT_YMAXIM=99.9999
C
      WRITE(NRESUL,'(/,''<<NO_OF_CURV>>'',1X,''NUMB_CURVE'',2X,
     *                 ''IF_BULLETS'',2X,
     *                 ''IF_HISTOGS'',2X,''IF_OVERFITING'')')
      WRITE(NRESUL,'(15X,I5,7X,I5,7X,I5,7X,I5,7X,I5)')LDCURV,
     *                                         IFBULL,IFHIST,IFOVER
      WRITE(NRESUL,'(''<<>>'')')
      WRITE(NRESUL,'(''<<CURVE_0000>>'',1X,
     *                  ''THICK_LINE'',2X,''STYLE_LINE'',2X,
     *                  ''COLOR_LINE'',2X,''TYPE_POINT'')')
      WRITE(NRESUL,'(15X,I5,7X,I5,7X,I5,7X,I5)')ITHICK,ISTYLE,
     *                                          ICOLOR,I_TYPE
      WRITE(NRESUL,'(15X,''NUMB_POINT'',2X,''YXPECTED_V'',
     *                   2X,''XPOSI_YMAX'')')
      WRITE(NRESUL,'(15X,I10,2X,F10.4,2X,F10.4)')
     *                      LDMONT,YEXPEC,XPOSIT_YMAXIM
C
C=======================================================================
C
      RETURN
      END
C
C=======================================================================
C=======================================================================
C
      SUBROUTINE PLOTCURVPRINT(NRESUL,NDCURV,LDCURV,IMONTE,
     *                                LABELS,VECTOR_OFDATA)
C
      CHARACTER
     *          LABELS*100
      DIMENSION
     *          LABELS(1:NDCURV)
      DIMENSION
     *          VECTOR_OFDATA(1:NDCURV)
C
C=======================================================================
C
      IF (IMONTE.EQ.1) THEN
          WRITE(NRESUL,'(21X,<LDCURV>(12X,I6.6))')(I,I=1,LDCURV)
          WRITE(NRESUL,'(15X,''IMONTE'',$)')
          DO I=1,LDCURV
             LENGTH=0
             DO J=100,1,-1
                IF (LABELS(I)(J:J).NE.' ') THEN
                    LENGTH=J
                    GO TO 1
                END IF
             END DO
   1         CONTINUE
             IF (LENGTH.EQ.5) LENGTH=LENGTH+1
             WRITE(NRESUL,'(12X,A,$)')LABELS(I)(1:LENGTH)
          END DO
          WRITE(NRESUL,'()')
      END IF
C
      WRITE(NRESUL,'(15X,I6,<LDCURV>(3X,F15.8))')IMONTE,
     *                 (VECTOR_OFDATA(I),I=1,LDCURV)
C
C=======================================================================
C
      RETURN
      END
C
C=======================================================================
C=======================================================================
C
      SUBROUTINE PLOTEXTRPRINT(NRESUL,XMINIM,XMAXIM,YMINIM,YMAXIM)
C
C=======================================================================
C
      WRITE(NRESUL,'()')
C
      WRITE(NRESUL,'(''<<X_AX_PARAM>> XMINIM_FIG  XMAXIM_FIG  '',
     *                   ''X_STEP_FIG'')')
C
      WRITE(NRESUL,'(15X,3(F10.2,2X),/)')XMINIM,XMAXIM,(XMAXIM-XMINIM)
C
      WRITE(NRESUL,'(''<<Y_AX_PARAM>> YMINIM_FIG  YMAXIM_FIG  '',
     *                   ''Y_STEP_FIG'')')
C
      WRITE(NRESUL,'(15X,3(F10.2,2X),/)')
     *                         YMINIM,YMAXIM,(YMAXIM-YMINIM)
C
C=======================================================================
C
      WRITE(NRESUL,'(''<<GO_GETTHEM>>'',/)')
C
C=======================================================================
C
      RETURN
      END
C
C=======================================================================
C=======================================================================
C
      SUBROUTINE SELECT_VSODEN
C
      INCLUDE 'MATDIM/NDPARS.f'
C
      COMMON
     *       /DENSTR/ IFDENS
      COMMON
     *       /TENSOR/ IFTENS
      COMMON
     *       /TENOPT/ ICENTT,ISORBT,ITENSR
      COMMON
     *       /POTPOT/ PARPOT(1:NDPARS)
      COMMON
     *       /AUXPOT/ VKACEN_PROTON,RKACEN_PROTON,AKACEN_PROTON,
     *                VKASOR_PROTON,RKASOR_PROTON,AKASOR_PROTON,
     *                VKAEFM_PROTON,RKAEFM_PROTON,AKAEFM_PROTON,
     *                                            V0UNIT_PROTON,
     *                VKACEN_NEUTRS,RKACEN_NEUTRS,AKACEN_NEUTRS,
     *                VKASOR_NEUTRS,RKASOR_NEUTRS,AKASOR_NEUTRS,
     *                VKAEFM_NEUTRS,RKAEFM_NEUTRS,AKAEFM_NEUTRS,
     *                                            V0UNIT_NEUTRS,
     *                                                   RKACOU,
     *                              ALAMPP,ALAMPN,ALAMNP,ALAMNN,
     *                              TLAMPP,TLAMPN,TLAMNP,TLAMNN,
     *                              CLAMPP,CLAMPN,CLAMNP,CLAMNN
      COMMON
     *       /IFLAMB/ IFPAR1,IFPAR2,IFPAR3,IFPAR4,
     *                IFPAR5,IFPAR6,IFPAR7,IFPAR8,
     *                IFPA09,IFPA10,IFPA11,IFPA12
      COMMON
     *       /PRFILE/ ISCREN,LOGFIL,LOGCHP,LOGCHN,LOGHIP,LOGHIN,
     *                LOGPRO,LOGNEU,LOGWRI,LSCREN,IRESUL,ICONVE,
     *                NOUTPT,LOGAUX,LOGENE,LOGRAD,LOGGRD,LOGMSH,
     *                                                   LOGBIS  
C
C=======================================================================
C
C     This SUBROUTINE "corrects" the values of the SO-Potential
C     lambda's depending on the IFLAMB block from the NAME-LIST
C
C
C  spin-orbit and tensor concerned
C
C=======================================================================
C      
      IF (IFDENS.EQ.0) THEN
          IF (LOGWRI.GT.5) THEN
              WRITE(LOGFIL,'(15X,''Exiting  SELECT_VSODEN at '',
     *                           ''IFDENS=0'')')
          END IF
          RETURN
      END IF
C
C=======================================================================
C 
      IF (LOGWRI.GT.0) THEN
C
          WRITE(LOGFIL,'(18x,''Lambdas before modification:'',
     *               '' PARPOT(39)= '',F10.4,'' PARPOT(40)= '',F10.4,
     *               '' PARPOT(41)= '',F10.4,'' PARPOT(42)= '',F10.4)') 
     *                  PARPOT(39),PARPOT(40),PARPOT(41),PARPOT(42) 
          WRITE(LOGFIL,'(46X,'' ALAMPP    = '',F10.4,
     *                       '' ALAMPN    = '',F10.4,
     *                       '' ALAMNP    = '',F10.4,
     *                       '' ALAMNN    = '',F10.4,/)')
     *                  ALAMPP,ALAMPN,ALAMNP,ALAMNN  
      END IF
C
C=======================================================================
C     
      IF (IFPAR1.EQ.1.AND.IFPAR2.EQ.-1) THEN
          PARPOT(40)=PARPOT(39)
          ALAMPN=ALAMPP
      END IF
C      
      IF (IFPAR1.EQ.1.AND.IFPAR3.EQ.-1) THEN
          PARPOT(41)=PARPOT(39)
          ALAMNP=ALAMPP
      END IF
C      
      IF (IFPAR1.EQ.1.AND.IFPAR4.EQ.-1) THEN
          PARPOT(42)=PARPOT(39)
          ALAMNN=ALAMPP
      END IF
C_______________________________________________________________________
C     
      IF (IFPAR2.EQ.2.AND.IFPAR3.EQ.-2) THEN
          PARPOT(41)=PARPOT(40)
          ALAMNP=ALAMPN
      END IF
C
C=======================================================================
C
      IF (LOGWRI.GT.0) THEN
C
          WRITE(LOGFIL,'(18x,''Lambdas after modification:'',
     *               '' PARPOT(39)= '',F10.4,'' PARPOT(40)= '',F10.4,
     *               '' PARPOT(41)= '',F10.4,'' PARPOT(42)= '',F10.4)') 
     *                  PARPOT(39),PARPOT(40),PARPOT(41),PARPOT(42) 
          WRITE(LOGFIL,'(45X,'' ALAMPP    = '',F10.4,
     *                       '' ALAMPN    = '',F10.4,
     *                       '' ALAMNP    = '',F10.4,
     *                       '' ALAMNN    = '',F10.4)')
     *                  ALAMPP,ALAMPN,ALAMNP,ALAMNN 
      END IF
C
C=======================================================================
C      
      IF (IFTENS.EQ.0) THEN
          IF (LOGWRI.GT.5) THEN
              WRITE(LOGFIL,'(15X,''Exiting  SELECT_VSODEN at '',
     *                           ''IFTENS=0'')')
          END IF
          RETURN
      END IF
C
C=======================================================================
C 
      IF (LOGWRI.GT.0) THEN
C
          WRITE(LOGFIL,'(18x,''Lambdas before modification:'',
     *               '' PARPOT(43)= '',F10.4,'' PARPOT(44)= '',F10.4,
     *               '' PARPOT(45)= '',F10.4,'' PARPOT(46)= '',F10.4)') 
     *                  PARPOT(43),PARPOT(44),PARPOT(45),PARPOT(46) 
          WRITE(LOGFIL,'(46X,'' TLAMPP    = '',F10.4,
     *                       '' TLAMPN    = '',F10.4,
     *                       '' TLAMNP    = '',F10.4,
     *                       '' TLAMNN    = '',F10.4,/)')
     *                  TLAMPP,TLAMPN,TLAMNP,TLAMNN 
      END IF
C
C=======================================================================
C     
      IF (IFPAR5.EQ.1.AND.IFPAR6.EQ.-1) THEN
          PARPOT(44)=PARPOT(43)
          TLAMPN=TLAMPP
      END IF
C      
      IF (IFPAR5.EQ.1.AND.IFPAR7.EQ.-1) THEN
          PARPOT(45)=PARPOT(43)
          TLAMNP=TLAMPP
      END IF
C      
      IF (IFPAR5.EQ.1.AND.IFPAR8.EQ.-1) THEN
          PARPOT(46)=PARPOT(43)
          TLAMNN=TLAMPP
      END IF
C_______________________________________________________________________
C      
      IF (IFPAR6.EQ.2.AND.IFPAR7.EQ.-2) THEN
          PARPOT(45)=PARPOT(44)
          TLAMNP=TLAMPN
      END IF
C_______________________________________________________________________
C
      IF (IFPAR5.EQ.3.AND.IFPAR6.EQ.-3) THEN
          PARPOT(44)=-PARPOT(43)
          TLAMPN=-TLAMPP
      END IF
C      
      IF (IFPAR5.EQ.3.AND.IFPAR7.EQ.-3) THEN
          PARPOT(45)=-PARPOT(43)
          TLAMNP=-TLAMPP
      END IF
C      
      IF (IFPAR5.EQ.3.AND.IFPAR8.EQ.-3) THEN
          PARPOT(46)=PARPOT(43)
          TLAMNN=TLAMPP
      END IF
C
C=======================================================================
C
      IF (LOGWRI.GT.0) THEN
C
          WRITE(LOGFIL,'(18x,''Lambdas after modification:'',
     *               '' PARPOT(43)= '',F10.4,'' PARPOT(44)= '',F10.4,
     *               '' PARPOT(45)= '',F10.4,'' PARPOT(46)= '',F10.4)') 
     *                  PARPOT(43),PARPOT(44),PARPOT(45),PARPOT(46) 
          WRITE(LOGFIL,'(46X,'' TLAMPP    = '',F10.4,
     *                       '' TLAMPN    = '',F10.4,
     *                       '' TLAMNP    = '',F10.4,
     *                       '' TLAMNN    = '',F10.4,/)')
     *                  TLAMPP,TLAMPN,TLAMNP,TLAMNN
      END IF
C
C=======================================================================
C 
      IF (LOGWRI.GT.0) THEN
C
          WRITE(LOGFIL,'(18x,''Lambdas before modification:'',
     *               '' PARPOT(47)= '',F10.4,'' PARPOT(48)= '',F10.4,
     *               '' PARPOT(49)= '',F10.4,'' PARPOT(50)= '',F10.4)') 
     *                  PARPOT(47),PARPOT(48),PARPOT(49),PARPOT(50) 
          WRITE(LOGFIL,'(46X,'' CLAMPP    = '',F10.4,
     *                       '' CLAMPN    = '',F10.4,
     *                       '' CLAMNP    = '',F10.4,
     *                       '' CLAMNN    = '',F10.4,/)')
     *                  CLAMPP,CLAMPN,CLAMNP,CLAMNN
      END IF
C
C=======================================================================
C     
      IF (IFPA09.EQ.1.AND.IFPA10.EQ.-1) THEN
          PARPOT(48)=PARPOT(47)
          CLAMPN=CLAMPP
      END IF
C      
      IF (IFPA09.EQ.1.AND.IFPA11.EQ.-1) THEN
          PARPOT(49)=PARPOT(47)
          CLAMNP=CLAMPP
      END IF
C      
      IF (IFPA09.EQ.1.AND.IFPA12.EQ.-1) THEN
          PARPOT(50)=PARPOT(47)
          CLAMNN=CLAMPP
      END IF
C_______________________________________________________________________
C      
      IF (IFPA10.EQ.2.AND.IFPA11.EQ.-2) THEN
          PARPOT(49)=PARPOT(48)
          CLAMNP=CLAMPN
      END IF
C_______________________________________________________________________
C
      IF (IFPA09.EQ.3.AND.IFPA10.EQ.-3) THEN
          PARPOT(48)=-PARPOT(47)
          CLAMPN=-CLAMPP
      END IF
C      
      IF (IFPA09.EQ.3.AND.IFPA11.EQ.-3) THEN
          PARPOT(49)=-PARPOT(47)
          CLAMNP=-CLAMPP
      END IF
C      
      IF (IFPA09.EQ.3.AND.IFPA12.EQ.-3) THEN
          PARPOT(50)=PARPOT(47)
          CLAMNN=CLAMPP
      END IF
C
C=======================================================================
C
      IF (LOGWRI.GT.0) THEN
C
          WRITE(LOGFIL,'(18x,''Lambdas before modification:'',
     *               '' PARPOT(47)= '',F10.4,'' PARPOT(48)= '',F10.4,
     *               '' PARPOT(49)= '',F10.4,'' PARPOT(50)= '',F10.4)') 
     *                  PARPOT(47),PARPOT(48),PARPOT(49),PARPOT(50) 
          WRITE(LOGFIL,'(46X,'' CLAMPP    = '',F10.4,
     *                       '' CLAMPN    = '',F10.4,
     *                       '' CLAMNP    = '',F10.4,
     *                       '' CLAMNN    = '',F10.4,/)')
     *                  CLAMPP,CLAMPN,CLAMNP,CLAMNN 
      END IF
C
C=======================================================================
C
      IF (LOGWRI.GT.5) THEN
          WRITE(LOGFIL,'(18X,''Exiting SELECT_VSODEN, normal exit'')')
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
      SUBROUTINE SPHLAB(IZ_FIX,IN_FIX,EFERMI_PROTON,EFERMI_NEUTRS)
C      
      INCLUDE   'MATDIM/NDSPEC.f'
      INCLUDE   'MATDIM/NDPARS.f'
      INCLUDE   'MATDIM/NDRAUS.f'
      INCLUDE   'MATDIM/NDNUCL.f'
      INCLUDE   'MATDIM/NDMAIN.f'
C      
      PARAMETER
     *         (NDORBI=NDMAIN,NDJTOT=2*NDORBI+1) 
      CHARACTER
     *          LABTHE_PROTON*6,LABTHE_NEUTRS*6,LBSHEL*6
      CHARACTER
     *          LABPRO_REMOVE*6,LABNEU_REMOVE*6,
     *          REMOVE_PROTON*3,REMOVE_NEUTRS*3
C
      DIMENSION
     *          LBSHEL(0:NDMAIN,0:NDORBI,1:NDJTOT)
C
      COMMON
     *       /DIMSPH/ LDSING_PROTON,LEVTHE_PROTON,
     *                LEVEXP_PROTON(1:NDNUCL),
     *                LDSING_NEUTRS,LEVTHE_NEUTRS,
     *                LEVEXP_NEUTRS(1:NDNUCL)
      COMMON
     *       /PRFILE/ ISCREN,LOGFIL,LOGCHP,LOGCHN,LOGHIP,LOGHIN,
     *                LOGPRO,LOGNEU,LOGWRI,LSCREN,IRESUL,ICONVE,
     *                NOUTPT,LOGAUX,LOGENE,LOGRAD,LOGGRD,LOGMSH,
     *                                                   LOGBIS          
      COMMON
     *       /REMLAB/ REMOVE_PROTON,REMOVE_NEUTRS
      COMMON
     *       /OUTLAB/ LABPRO_REMOVE(1:NDRAUS),
     *                LABNEU_REMOVE(1:NDRAUS)
      COMMON
     *       /OUTIND/ INDPRO_LETOUT(1:NDRAUS),
     *                INDNEU_LETOUT(1:NDRAUS)
      COMMON
     *       /SPECTR/ ENETHE_PROTON(1:NDSPEC),
     *                ENETHE_NEUTRS(1:NDSPEC)
      COMMON
     *       /THELAB/ LABTHE_PROTON(1:NDSPEC),
     *                LABTHE_NEUTRS(1:NDSPEC)
      COMMON
     *       /DEGENE/ IDGLEV_PROTON(1:NDSPEC),
     *                ICUMUL_PROTON(0:NDSPEC),
     *                IDGLEV_NEUTRS(1:NDSPEC),
     *                ICUMUL_NEUTRS(0:NDSPEC)
      COMMON
     *       /QUANTN/ NWSSPH_PROTON(1:NDSPEC),
     *                LWSSPH_PROTON(1:NDSPEC),
     *                JWSSPH_PROTON(1:NDSPEC),
     *
     *                NWSSPH_NEUTRS(1:NDSPEC),
     *                LWSSPH_NEUTRS(1:NDSPEC),
     *                JWSSPH_NEUTRS(1:NDSPEC)
      COMMON
     *       /AUXPOT/ VKACEN_PROTON,RKACEN_PROTON,AKACEN_PROTON,
     *                VKASOR_PROTON,RKASOR_PROTON,AKASOR_PROTON,
     *                VKAEFM_PROTON,RKAEFM_PROTON,AKAEFM_PROTON,
     *                                            V0UNIT_PROTON,
     *                VKACEN_NEUTRS,RKACEN_NEUTRS,AKACEN_NEUTRS,
     *                VKASOR_NEUTRS,RKASOR_NEUTRS,AKASOR_NEUTRS,
     *                VKAEFM_NEUTRS,RKAEFM_NEUTRS,AKAEFM_NEUTRS,
     *                                            V0UNIT_NEUTRS,
     *                                                   RKACOU,
     *                              ALAMPP,ALAMPN,ALAMNP,ALAMNN,
     *                              TLAMPP,TLAMPN,TLAMNP,TLAMNN,
     *                              CLAMPP,CLAMPN,CLAMNP,CLAMNN
      COMMON
     *       /RANRAN/ IRMFST,I_RAND,IRANDO
      COMMON
     *       /POTPOT/ PARPOT(1:NDPARS) 
      COMMON
     *       /MINITE/ IDEFCN,ITECHI
      COMMON
     *       /ACTIVE/ IACTIV,IwMODE
      COMMON
     *       /ENEPRI/ ENEMAX
C_______________________________________________________________________
C
      DATA !      n  L  J
     *     LBSHEL(00,00,01) / '1s1/2 ' /   ! Shell no. 0
      DATA
     *     LBSHEL(01,01,03) / '1p3/2 ' /, 
     *     LBSHEL(01,01,01) / '1p1/2 ' /   ! Shell no. 1
      DATA
     *     LBSHEL(02,00,01) / '2s1/2 ' /, 
     *     LBSHEL(02,02,05) / '1d5/2 ' /, 
     *     LBSHEL(02,02,03) / '1d3/2 ' /   ! Shell no. 2
      DATA
     *     LBSHEL(03,01,03) / '2p3/2 ' /, 
     *     LBSHEL(03,01,01) / '2p1/2 ' /, 
     *     LBSHEL(03,03,07) / '1f7/2 ' /, 
     *     LBSHEL(03,03,05) / '1f5/2 ' /   ! Shell no. 3
      DATA
     *     LBSHEL(04,00,01) / '3s1/2 ' /, 
     *     LBSHEL(04,02,05) / '2d5/2 ' /, 
     *     LBSHEL(04,02,03) / '2d3/2 ' /, 
     *     LBSHEL(04,04,09) / '1g9/2 ' /, 
     *     LBSHEL(04,04,07) / '1g7/2 ' /   ! Shell no. 4
      DATA
     *     LBSHEL(05,01,03) / '3p3/2 ' /, 
     *     LBSHEL(05,01,01) / '3p1/2 ' /, 
     *     LBSHEL(05,03,07) / '2f7/2 ' /, 
     *     LBSHEL(05,03,05) / '2f5/2 ' /, 
     *     LBSHEL(05,05,11) / '1h11/2' /, 
     *     LBSHEL(05,05,09) / '1h9/2 ' /   ! Shell no. 5
      DATA
     *     LBSHEL(06,00,01) / '4s1/2 ' /, 
     *     LBSHEL(06,02,05) / '3d5/2 ' /, 
     *     LBSHEL(06,02,03) / '3d3/2 ' /, 
     *     LBSHEL(06,04,09) / '2g9/2 ' /, 
     *     LBSHEL(06,04,07) / '2g7/2 ' /, 
     *     LBSHEL(06,06,13) / '1i13/2' /, 
     *     LBSHEL(06,06,11) / '1i11/2' /   ! Shell no. 6
      DATA
     *     LBSHEL(07,01,03) / '4p3/2 ' /, 
     *     LBSHEL(07,01,01) / '4p1/2 ' /, 
     *     LBSHEL(07,03,07) / '3f7/2 ' /, 
     *     LBSHEL(07,03,05) / '3f5/2 ' /, 
     *     LBSHEL(07,05,11) / '2h11/2' /, 
     *     LBSHEL(07,05,09) / '2h9/2 ' /, 
     *     LBSHEL(07,07,15) / '1j15/2' /, 
     *     LBSHEL(07,07,13) / '1j13/2' /   ! Shell no. 7
      DATA
     *     LBSHEL(08,00,01) / '5s1/2 ' /, 
     *     LBSHEL(08,02,05) / '4d5/2 ' /, 
     *     LBSHEL(08,02,03) / '4d3/2 ' /, 
     *     LBSHEL(08,04,09) / '3g9/2 ' /, 
     *     LBSHEL(08,04,07) / '3g7/2 ' /, 
     *     LBSHEL(08,06,13) / '2i13/2' /, 
     *     LBSHEL(08,06,11) / '2i11/2' /, 
     *     LBSHEL(08,08,17) / '1k17/2' /, 
     *     LBSHEL(08,08,15) / '1k15/2' /   ! Shell no. 8
      DATA
     *     LBSHEL(09,01,03) / '5p3/2 ' /, 
     *     LBSHEL(09,01,01) / '5p1/2 ' /, 
     *     LBSHEL(09,03,07) / '4f7/2 ' /, 
     *     LBSHEL(09,03,05) / '4f5/2 ' /, 
     *     LBSHEL(09,05,11) / '3h11/2' /, 
     *     LBSHEL(09,05,09) / '3h9/2 ' /, 
     *     LBSHEL(09,07,15) / '2j15/2' /, 
     *     LBSHEL(09,07,13) / '2j13/2' /, 
     *     LBSHEL(09,09,19) / '1l19/2' /, 
     *     LBSHEL(09,09,17) / '1l17/2' /   ! Shell no. 9
      DATA
     *     LBSHEL(10,00,01) / '6s1/2 ' /, 
     *     LBSHEL(10,02,05) / '5d5/2 ' /, 
     *     LBSHEL(10,02,03) / '5d3/2 ' /, 
     *     LBSHEL(10,04,09) / '4g9/2 ' /, 
     *     LBSHEL(10,04,07) / '4g7/2 ' /, 
     *     LBSHEL(10,06,13) / '3i13/2' /, 
     *     LBSHEL(10,06,11) / '3i11/2' /, 
     *     LBSHEL(10,08,17) / '2k17/2' /, 
     *     LBSHEL(10,08,15) / '2k15/2' /, 
     *     LBSHEL(10,10,21) / '1m21/2' /, 
     *     LBSHEL(10,10,19) / '1m19/2' /   ! Shell no.10
      DATA
     *     LBSHEL(11,01,03) / '6p3/2 ' /, 
     *     LBSHEL(11,01,01) / '6p1/2 ' /, 
     *     LBSHEL(11,03,07) / '5f7/2 ' /, 
     *     LBSHEL(11,03,05) / '5f5/2 ' /, 
     *     LBSHEL(11,05,11) / '4h11/2' /, 
     *     LBSHEL(11,05,09) / '4h9/2 ' /, 
     *     LBSHEL(11,07,15) / '3j15/2' /, 
     *     LBSHEL(11,07,13) / '3j13/2' /, 
     *     LBSHEL(11,09,19) / '2l19/2' /, 
     *     LBSHEL(11,09,17) / '2l17/2' /, 
     *     LBSHEL(11,11,23) / '1n23/2' /, 
     *     LBSHEL(11,11,21) / '1n21/2' /   ! Shell no.11
      DATA
     *     LBSHEL(12,00,01) / '7s1/2 ' /, 
     *     LBSHEL(12,02,05) / '6d5/2 ' /, 
     *     LBSHEL(12,02,03) / '6d3/2 ' /, 
     *     LBSHEL(12,04,09) / '5g9/2 ' /, 
     *     LBSHEL(12,04,07) / '5g7/2 ' /, 
     *     LBSHEL(12,06,13) / '4i13/2' /, 
     *     LBSHEL(12,06,11) / '4i11/2' /, 
     *     LBSHEL(12,08,17) / '3k17/2' /, 
     *     LBSHEL(12,08,15) / '3k15/2' /, 
     *     LBSHEL(12,10,21) / '2m21/2' /, 
     *     LBSHEL(12,10,19) / '2m19/2' /, 
     *     LBSHEL(12,12,25) / '2o25/2' /, 
     *     LBSHEL(12,12,23) / '2o23/2' /   ! Shell no.12
      DATA
     *     LBSHEL(13,01,03) / '7p3/2 ' /, 
     *     LBSHEL(13,01,01) / '7p1/2 ' /, 
     *     LBSHEL(13,03,07) / '6f7/2 ' /, 
     *     LBSHEL(13,03,05) / '6f5/2 ' /, 
     *     LBSHEL(13,05,11) / '5h11/2' /, 
     *     LBSHEL(13,05,09) / '5h9/2 ' /, 
     *     LBSHEL(13,07,15) / '4j15/2' /, 
     *     LBSHEL(13,07,13) / '4j13/2' /, 
     *     LBSHEL(13,09,19) / '3l19/2' /, 
     *     LBSHEL(13,09,17) / '3l17/2' /, 
     *     LBSHEL(13,11,23) / '2n23/2' /, 
     *     LBSHEL(13,11,21) / '2n21/2' /, 
     *     LBSHEL(13,13,27) / '1r27/2' /, 
     *     LBSHEL(13,13,25) / '1r25/2' /   ! Shell no.13
      DATA
     *     LBSHEL(14,00,01) / '8s1/2 ' /, 
     *     LBSHEL(14,02,05) / '7d5/2 ' /, 
     *     LBSHEL(14,02,03) / '7d3/2 ' /, 
     *     LBSHEL(14,04,09) / '6g9/2 ' /, 
     *     LBSHEL(14,04,07) / '6g7/2 ' /, 
     *     LBSHEL(14,06,13) / '5i13/2' /, 
     *     LBSHEL(14,06,11) / '5i11/2' /, 
     *     LBSHEL(14,08,17) / '4k17/2' /, 
     *     LBSHEL(14,08,15) / '4k15/2' /, 
     *     LBSHEL(14,10,21) / '3m21/2' /, 
     *     LBSHEL(14,10,19) / '3m19/2' /, 
     *     LBSHEL(14,12,25) / '2o25/2' /, 
     *     LBSHEL(14,12,23) / '2o23/2' /,
     *     LBSHEL(14,14,29) / '1q29/2' /, 
     *     LBSHEL(14,14,27) / '1q27/2' /   ! Shell no.14
      DATA
     *     LBSHEL(15,01,03) / '8p3/2 ' /, 
     *     LBSHEL(15,01,01) / '8p1/2 ' /, 
     *     LBSHEL(15,03,07) / '7f7/2 ' /, 
     *     LBSHEL(15,03,05) / '7f5/2 ' /, 
     *     LBSHEL(15,05,11) / '6h11/2' /, 
     *     LBSHEL(15,05,09) / '6h9/2 ' /, 
     *     LBSHEL(15,07,15) / '5j15/2' /, 
     *     LBSHEL(15,07,13) / '5j13/2' /, 
     *     LBSHEL(15,09,19) / '4l19/2' /, 
     *     LBSHEL(15,09,17) / '4l17/2' /, 
     *     LBSHEL(15,11,23) / '3n23/2' /, 
     *     LBSHEL(15,11,21) / '3n21/2' /, 
     *     LBSHEL(15,13,27) / '2r27/2' /, 
     *     LBSHEL(15,13,25) / '2r25/2' /, 
     *     LBSHEL(15,15,31) / '1t31/2' /, 
     *     LBSHEL(15,15,29) / '1t29/2' /   ! Shell no.15
      DATA
     *     LBSHEL(16,00,01) / '9s1/2 ' /, 
     *     LBSHEL(16,02,05) / '8d5/2 ' /, 
     *     LBSHEL(16,02,03) / '8d3/2 ' /, 
     *     LBSHEL(16,04,09) / '7g9/2 ' /, 
     *     LBSHEL(16,04,07) / '7g7/2 ' /, 
     *     LBSHEL(16,06,13) / '6i13/2' /, 
     *     LBSHEL(16,06,11) / '6i11/2' /, 
     *     LBSHEL(16,08,17) / '5k17/2' /, 
     *     LBSHEL(16,08,15) / '5k15/2' /, 
     *     LBSHEL(16,10,21) / '4m21/2' /, 
     *     LBSHEL(16,10,19) / '4m19/2' /, 
     *     LBSHEL(16,12,25) / '3o25/2' /, 
     *     LBSHEL(16,12,23) / '3o23/2' /,
     *     LBSHEL(16,14,29) / '2q29/2' /, 
     *     LBSHEL(16,14,27) / '2q27/2' /,
     *     LBSHEL(16,16,33) / '1u33/2' /, 
     *     LBSHEL(16,16,31) / '1u31/2' /   ! Shell no.16
      DATA
     *     LBSHEL(17,01,03) / '9p3/2 ' /, 
     *     LBSHEL(17,01,01) / '9p1/2 ' /, 
     *     LBSHEL(17,03,07) / '8f7/2 ' /, 
     *     LBSHEL(17,03,05) / '8f5/2 ' /, 
     *     LBSHEL(17,05,11) / '7h11/2' /, 
     *     LBSHEL(17,05,09) / '7h9/2 ' /, 
     *     LBSHEL(17,07,15) / '6j15/2' /, 
     *     LBSHEL(17,07,13) / '6j13/2' /, 
     *     LBSHEL(17,09,19) / '5l19/2' /, 
     *     LBSHEL(17,09,17) / '5l17/2' /, 
     *     LBSHEL(17,11,23) / '4n23/2' /, 
     *     LBSHEL(17,11,21) / '4n21/2' /, 
     *     LBSHEL(17,13,27) / '3r27/2' /, 
     *     LBSHEL(17,13,25) / '3r25/2' /, 
     *     LBSHEL(17,15,31) / '2t31/2' /, 
     *     LBSHEL(17,15,29) / '2t29/2' /, 
     *     LBSHEL(17,17,35) / '1v35/2' /, 
     *     LBSHEL(17,17,33) / '1v33/2' /   ! Shell no.17
      DATA
     *     LBSHEL(18,00,01) / '10s1/2' /, 
     *     LBSHEL(18,02,05) / '9d5/2 ' /, 
     *     LBSHEL(18,02,03) / '9d3/2 ' /, 
     *     LBSHEL(18,04,09) / '8g9/2 ' /, 
     *     LBSHEL(18,04,07) / '8g7/2 ' /, 
     *     LBSHEL(18,06,13) / '7i13/2' /, 
     *     LBSHEL(18,06,11) / '7i11/2' /, 
     *     LBSHEL(18,08,17) / '6k17/2' /, 
     *     LBSHEL(18,08,15) / '6k15/2' /, 
     *     LBSHEL(18,10,21) / '5m21/2' /, 
     *     LBSHEL(18,10,19) / '5m19/2' /, 
     *     LBSHEL(18,12,25) / '4o25/2' /, 
     *     LBSHEL(18,12,23) / '4o23/2' /,
     *     LBSHEL(18,14,29) / '3q29/2' /, 
     *     LBSHEL(18,14,27) / '3q27/2' /,
     *     LBSHEL(18,16,33) / '2u33/2' /, 
     *     LBSHEL(18,16,31) / '2u31/2' /,
     *     LBSHEL(18,18,37) / '1w39/2' /, 
     *     LBSHEL(18,18,35) / '1w37/2' /   ! Shell no.18
      DATA
     *     LBSHEL(19,01,03) / '10p3/2' /, 
     *     LBSHEL(19,01,01) / '10p1/2' /, 
     *     LBSHEL(19,03,07) / '9f7/2 ' /, 
     *     LBSHEL(19,03,05) / '9f5/2 ' /, 
     *     LBSHEL(19,05,11) / '8h11/2' /, 
     *     LBSHEL(19,05,09) / '8h9/2 ' /, 
     *     LBSHEL(19,07,15) / '7j15/2' /, 
     *     LBSHEL(19,07,13) / '7j13/2' /, 
     *     LBSHEL(19,09,19) / '6l19/2' /, 
     *     LBSHEL(19,09,17) / '6l17/2' /, 
     *     LBSHEL(19,11,23) / '5n23/2' /, 
     *     LBSHEL(19,11,21) / '5n21/2' /, 
     *     LBSHEL(19,13,27) / '4r27/2' /, 
     *     LBSHEL(19,13,25) / '4r25/2' /, 
     *     LBSHEL(19,15,31) / '3t31/2' /, 
     *     LBSHEL(19,15,29) / '3t29/2' /, 
     *     LBSHEL(19,17,35) / '2v35/2' /, 
     *     LBSHEL(19,17,33) / '2v33/2' /, 
     *     LBSHEL(19,19,39) / '1x39/2' /, 
     *     LBSHEL(19,19,37) / '1x37/2' /   ! Shell no.19
      DATA
     *     LBSHEL(20,00,01) / '11s1/2' /, 
     *     LBSHEL(20,02,05) / '10d5/2' /, 
     *     LBSHEL(20,02,03) / '10d3/2' /, 
     *     LBSHEL(20,04,09) / '9g9/2 ' /, 
     *     LBSHEL(20,04,07) / '9g7/2 ' /, 
     *     LBSHEL(20,06,13) / '8i13/2' /, 
     *     LBSHEL(20,06,11) / '8i11/2' /, 
     *     LBSHEL(20,08,17) / '7k17/2' /, 
     *     LBSHEL(20,08,15) / '7k15/2' /, 
     *     LBSHEL(20,10,21) / '6m21/2' /, 
     *     LBSHEL(20,10,19) / '6m19/2' /, 
     *     LBSHEL(20,12,25) / '5o25/2' /, 
     *     LBSHEL(20,12,23) / '5o23/2' /,
     *     LBSHEL(20,14,29) / '4q29/2' /, 
     *     LBSHEL(20,14,27) / '4q27/2' /,
     *     LBSHEL(20,16,33) / '3u33/2' /, 
     *     LBSHEL(20,16,31) / '3u31/2' /,
     *     LBSHEL(20,18,37) / '2w39/2' /, 
     *     LBSHEL(20,18,35) / '2w37/2' /,
     *     LBSHEL(20,20,41) / '1y43/2' /, 
     *     LBSHEL(20,20,39) / '1y41/2' /   ! Shell no.20
C
C=======================================================================
C
C     This subroutine prepares the spherical labels e.g. for the
C     plotting system and for the chi^2  minimisation and prints
C     the corresponding table. 
C
C     Here we identify the  'spherical multiplets'   by a direct 
C     identification of the L and J quantum numbers. 
C
C     Called from WSSTAN and HAMMAT
C
C=======================================================================
C
      IF (IwMODE.EQ.1) RETURN ! No global table print when doing 
C                                                   the gradient
C=======================================================================
C
      IF (IwMODE.EQ.2 .AND. LOGWRI.GT.4) THEN
C
          WRITE(IRESUL,'(/,27x,''The very final Result Table'',/)')
C
          WRITE(LOGPRO,'(/,''IDEFCN='',I3,'' IRANDO='',I3)')
     *                       IDEFCN,         IRANDO
          WRITE(LOGNEU,'(/,''IDEFCN='',I3,'' IRANDO='',I3)')
     *                       IDEFCN,         IRANDO
C
          IF (REMOVE_PROTON.EQ.'YES') THEN
              WRITE(LOGPRO,'(/,''The following experimental proton ''
     *                         ''levels  were not included in the fit'',
     *                                                              /)')
C
              DO I=1,NDRAUS
                 IF (ISCREN.EQ.1.AND.INDPRO_LETOUT(I).EQ.1) THEN
                     WRITE(LOGPRO,'(''Removed proton '',a6,'' level'')')
     *                     LABPRO_REMOVE(I)                            
                 END IF
              END DO
          END IF
C
          IF (REMOVE_NEUTRS.EQ.'YES') THEN
              WRITE(LOGNEU,'(/,''The following experimental neutron '',
     *                         ''levels were not included in the fit:'',
     *                                                              /)')
              DO I=1,NDRAUS
                 IF (INDNEU_LETOUT(I).EQ.1) THEN
                     WRITE(LOGNEU,'(''Removed neutron '',a6,
     *                                          '' level'')')
     *                     LABNEU_REMOVE(I)                            
                 END IF
              END DO
          END IF
C
      END IF
C
C=======================================================================
C
      V0CENT_PROTON=PARPOT( 1)
      R0CENT_PROTON=PARPOT( 2)
      A0CENT_PROTON=PARPOT( 3)
      V0SORB_PROTON=PARPOT( 4)
      R0SORB_PROTON=PARPOT( 5)
      A0SORB_PROTON=PARPOT( 6)
      V0EFFM_PROTON=PARPOT( 7)
      R0EFFM_PROTON=PARPOT( 8)
      A0EFFM_PROTON=PARPOT( 9)
      R0COUL       =PARPOT(10)
     
      XK_V0C_PROTON=PARPOT(11)
      XK_R0C_PROTON=PARPOT(12)
      XK_A0C_PROTON=PARPOT(13)
      XK_LAM_PROTON=PARPOT(14)
      XK_RSO_PROTON=PARPOT(15)
      XK_ASO_PROTON=PARPOT(16)
      XK_LEF_PROTON=PARPOT(17)
      XK_REF_PROTON=PARPOT(18)
      XK_AEF_PROTON=PARPOT(19)
      XK_COU	   =PARPOT(20)
     
      V0CENT_NEUTRS=PARPOT(21)
      R0CENT_NEUTRS=PARPOT(22)
      A0CENT_NEUTRS=PARPOT(23)
      V0SORB_NEUTRS=PARPOT(24)
      R0SORB_NEUTRS=PARPOT(25)
      A0SORB_NEUTRS=PARPOT(26)
      V0EFFM_NEUTRS=PARPOT(27)
      R0EFFM_NEUTRS=PARPOT(28)
      A0EFFM_NEUTRS=PARPOT(29)
     
      XK_V0C_NEUTRS=PARPOT(30)
      XK_R0C_NEUTRS=PARPOT(31)
      XK_A0C_NEUTRS=PARPOT(32)
      XK_LAM_NEUTRS=PARPOT(33)
      XK_RSO_NEUTRS=PARPOT(34)
      XK_ASO_NEUTRS=PARPOT(35)
      XK_LEF_NEUTRS=PARPOT(36)
      XK_REF_NEUTRS=PARPOT(37)
      XK_AEF_NEUTRS=PARPOT(38)
     
      ALAMPP       =PARPOT(39)
      ALAMPN       =PARPOT(40)
      ALAMNP       =PARPOT(41)
      ALAMNN       =PARPOT(42)

      TLAMPP       =PARPOT(43)
      TLAMPN       =PARPOT(44)
      TLAMNP       =PARPOT(45)
      TLAMNN       =PARPOT(46)
      
      CLAMPP       =PARPOT(47)
      CLAMPN       =PARPOT(48)
      CLAMNP       =PARPOT(49)
      CLAMNN       =PARPOT(50)

      IF (LDSING_PROTON.GT.NDSPEC) THEN
          STOP 'LDSING_PROTON.GT.NDSPEC in SPHLAB'
      END IF
C
C=======================================================================
C
      N_PART=IZ_FIX
C
C=======================================================================
C
      ICUMUL_PROTON(0)=0
C
      DO JSTATE=1,LDSING_PROTON
C
         INDEXN=LWSSPH_PROTON(JSTATE)
     *         +NWSSPH_PROTON(JSTATE)
     *         +NWSSPH_PROTON(JSTATE)
         INDEXL=LWSSPH_PROTON(JSTATE)
         INDEXJ=JWSSPH_PROTON(JSTATE)
C
         LABTHE_PROTON(JSTATE)=LBSHEL(INDEXN,INDEXL,INDEXJ)
         IDGLEV_PROTON(JSTATE)=JWSSPH_PROTON(JSTATE)+1
C
         ICUMUL_PROTON(JSTATE)=ICUMUL_PROTON(JSTATE-1)
     *                        +IDGLEV_PROTON(JSTATE) 
C
         IF (ICUMUL_PROTON(JSTATE).EQ.N_PART) THEN
             EFERMI_PROTON=ENETHE_PROTON(JSTATE)+1.0e-8
         END IF
C
      END DO
C
C=======================================================================
C
      IF (LOGWRI.GT.4) THEN
C
          WRITE(IRESUL,'(/,80(''#''),/,''#'',78X,''#'',/,
     *               ''#   Central Potential:    V0CENT='',F8.4,
     *               ''   A0CENT='',F8.4,''  R0CENT='',F8.4,
     *               ''   #'',/,''#'',T62,''R0COUL='',F8.4,
     *               ''   #'')')
     *               V0CENT_PROTON,A0CENT_PROTON,R0CENT_PROTON,R0COUL
          WRITE(IRESUL,'(''#                         VKACEN='',F8.4,
     *               ''   AKACEN='',F8.4,''  RKACEN='',F8.4,
     *               ''   #'',/,''#'',T62,''RKACOU='',F8.4,
     *               ''   #'')')
     *               VKACEN_PROTON,AKACEN_PROTON,RKACEN_PROTON,
     *                                                  RKACOU
          WRITE(IRESUL,'(''#   S-Orbit Potential:    V0SORB='',F8.4,
     *               ''   A0SORB='',F8.4,''  R0SORB='',F8.4,
     *               ''   #'',78X,''#'')')
     *               V0SORB_PROTON,A0SORB_PROTON,R0SORB_PROTON
          WRITE(IRESUL,'(''#                         VKASOR='',F8.4,
     *               ''   AKASOR='',F8.4,''  RKASOR='',F8.4,
     *               ''   #'',/,''#'',78X,''#'',/,80(''#''))')
     *               VKASOR_PROTON,AKASOR_PROTON,RKASOR_PROTON

          WRITE(IRESUL,'(/,80(''#''),/,''#'',78X,''#'',/,
     *        ''#   Proton  Spherical Woods-Saxon Spectrum'',
     *        '' (LDSING='',I3,'')'',T63,
     *        ''Z='',I3,''   N='',I3,''    #'')') 
     *         LDSING_PROTON,IZ_FIX,IN_FIX
C
          WRITE(IRESUL,'(  ''#'',78X,''#'',/,80(''#''),             /,
     *            ''#'',78X,''#'',/,
     *            ''#   No)     Energy     State   Occup   '',
     *            ''    No)     Energy     State   Occup    #'',/,
     *            ''#'',78X,''#'')')
      END IF
C
C=======================================================================
C
      IF (IwMODE.EQ.2 .AND. LOGWRI.GT.4) THEN
C
          WRITE(LOGPRO,'(/,80(''#''),/,''#'',78X,''#'',/,
     *               ''#   Central Potential:     V0CENT='',F8.4,
     *               ''    A0CENT='',F6.4,''   R0CENT='',F6.4,
     *               ''    #'',/,''#'',T63,''R0COUL='',F6.4,
     *               ''    #'',/,
     *               ''#   S-Orbit Potential:     V0SORB='',F8.4,
     *               ''    A0SORB='',F6.4,''   R0SORB='',F6.4,
     *               ''    #'',/,''#'',78X,''#'',/,80(''#''))')
     *               V0CENT_PROTON,A0CENT_PROTON,R0CENT_PROTON,R0COUL,
     *               V0SORB_PROTON,A0SORB_PROTON,R0SORB_PROTON
C
          WRITE(LOGPRO,'(/,80(''#''),/,''#'',78X,''#'',/,
     *        ''#   Proton  Spherical Woods-Saxon Spectrum'',
     *        '' (LDSING='',I3,'')'',T63,
     *        ''Z='',I3,''   N='',I3,''    #'')') 
     *         LDSING_PROTON,IZ_FIX,IN_FIX
C
          WRITE(LOGPRO,'(  ''#'',78X,''#'',/,80(''#''),             /,
     *            ''#'',78X,''#'',/,
     *            ''#   No)     Energy     State   Occup   '',
     *            ''    No)     Energy     State   Occup    #'',/,
     *            ''#'',78X,''#'')')
      END IF
C
C=======================================================================
C
      DO ISTATE=1,LDSING_PROTON
         IF (ENETHE_PROTON(ISTATE).GT.ENEMAX) GO TO 1
      END DO
C
   1  CONTINUE
C   
      LDLINE=(ISTATE-1)/2
C
      DO LINE=1,LDLINE
C
         LINE_1=LINE
         LINE_2=LDLINE+LINE
C
         ICUM_1=ICUMUL_PROTON(LINE_1)
C
         ICUM_2=ICUMUL_PROTON(LINE_2)
C
         IF (LOGWRI.GT.4)
     *       WRITE(IRESUL,'(''#  '', I3,'')'',F12.5,4X,A6,I6,'' '',
     *                  ''      '',I3,'')'',F12.5,4X,A6,I6,'' '',
     *                                                 ''    #'')')
     *                LINE_1,ENETHE_PROTON(LINE_1),
     *                       LABTHE_PROTON(LINE_1),ICUM_1,
     *                LINE_2,ENETHE_PROTON(LINE_2),
     *                       LABTHE_PROTON(LINE_2),ICUM_2
C
         IF (IwMODE.EQ.2 .AND. LOGWRI.GT.4) THEN
C
             WRITE(LOGPRO,'(''#  '', I3,'')'',F12.5,4X,A6,I6,'' '',
     *                      ''      '',I3,'')'',F12.5,4X,A6,I6,'' '',
     *                                                   ''    #'')')
     *            LINE_1,ENETHE_PROTON(LINE_1),
     *                   LABTHE_PROTON(LINE_1),ICUM_1,
     *            LINE_2,ENETHE_PROTON(LINE_2),
     *                   LABTHE_PROTON(LINE_2),ICUM_2
         END IF
C
      END DO
C
C=======================================================================
C
      IF (LOGWRI.GT.4) WRITE(IRESUL,'(''#'',78X,''#'',/,80(''#''))')
C
      IF (IwMODE.EQ.2  .AND. LOGWRI.GT.4) THEN
          WRITE(LOGPRO,'(''#'',78X,''#'',/,80(''#''))')
      END IF
C
C=======================================================================
C=======================================================================
C     Repeating the same algorithm for the neutrons
C=======================================================================
C=======================================================================
C
      N_PART=IN_FIX
C
C=======================================================================
C
      ICUMUL_NEUTRS(0)=0
C
      DO JSTATE=1,LDSING_NEUTRS
C
         INDEXN=LWSSPH_NEUTRS(JSTATE)
     *         +NWSSPH_NEUTRS(JSTATE)
     *         +NWSSPH_NEUTRS(JSTATE)
         INDEXL=LWSSPH_NEUTRS(JSTATE)
         INDEXJ=JWSSPH_NEUTRS(JSTATE)
C
         LABTHE_NEUTRS(JSTATE)=LBSHEL(INDEXN,INDEXL,INDEXJ)
         IDGLEV_NEUTRS(JSTATE)=JWSSPH_NEUTRS(JSTATE)+1
C
         ICUMUL_NEUTRS(JSTATE)=ICUMUL_NEUTRS(JSTATE-1)
     *                        +IDGLEV_NEUTRS(JSTATE) 
C @@@
         IF (ICUMUL_NEUTRS(JSTATE).EQ.N_PART) THEN
             EFERMI_NEUTRS=ENETHE_NEUTRS(JSTATE)+1.0e-8
         END IF
C
      END DO
C
      IF (LOGWRI.GT.4) THEN
C
          WRITE(IRESUL,'(/,80(''#''),/,''#'',78X,''#'',/,
     *               ''#   Central Potential:    V0CENT='',F8.4,
     *               ''   A0CENT='',F8.4,''  R0CENT='',F8.4,
     *               ''   #'',78X,''#'')')
     *               V0CENT_NEUTRS,A0CENT_NEUTRS,R0CENT_NEUTRS
          WRITE(IRESUL,'(''#                         VKACEN='',F8.4,
     *               ''   AKACEN='',F8.4,''  RKACEN='',F8.4,
     *               ''   #'',78X,''#'')')
     *               VKACEN_NEUTRS,AKACEN_NEUTRS,RKACEN_NEUTRS
          WRITE(IRESUL,'(''#   S-Orbit Potential:    V0SORB='',F8.4,
     *               ''   A0SORB='',F8.4,''  R0SORB='',F8.4,
     *               ''   #'',78X,''#'')')
     *               V0SORB_NEUTRS,A0SORB_NEUTRS,R0SORB_NEUTRS
          WRITE(IRESUL,'(''#                         VKASOR='',F8.4,
     *               ''   AKASOR='',F8.4,''  RKASOR='',F8.4,
     *               ''   #'',/,''#'',78X,''#'',/,80(''#''))')
     *               VKASOR_NEUTRS,AKASOR_NEUTRS,RKASOR_NEUTRS
C
          WRITE(IRESUL,'(/,80(''#''),/,''#'',78X,''#'',/,
     *        ''#   Neutron Spherical Woods-Saxon Spectrum '',
     *        '' (LDSING='',I3,'')'',T63,
     *        ''Z='',I3,''   N='',I3,''    #'')') 
     *          LDSING_NEUTRS,IZ_FIX,IN_FIX
C
          WRITE(IRESUL,'(  ''#'',78X,''#'',/,80(''#''),             /,
     *            ''#'',78X,''#'',/,
     *            ''#   No)     Energy     State   Occup   '',
     *            ''    No)     Energy     State   Occup    #'',/,
     *            ''#'',78X,''#'')')
      END IF
C
      IF (IwMODE.EQ.2 .AND. LOGWRI.GT.4) THEN
C
          WRITE(LOGNEU,'(/,80(''#''),/,''#'',78X,''#'',/,
     *               ''#   Central Potential:     V0CENT='',F8.4,
     *               ''    A0CENT='',F6.4,''   R0CENT='',F6.4,
     *               ''    #'',/,''#'',78X,''#'',/,
     *               ''#   S-Orbit Potential:     V0SORB='',F8.4,
     *               ''    A0SORB='',F6.4,''   R0SORB='',F6.4,
     *               ''    #'',/,''#'',78X,''#'',/,80(''#''))')
     *               V0CENT_NEUTRS,A0CENT_NEUTRS,R0CENT_NEUTRS,
     *               V0SORB_NEUTRS,A0SORB_NEUTRS,R0SORB_NEUTRS
C
          WRITE(LOGNEU,'(/,80(''#''),/,''#'',78X,''#'',/,
     *               ''#   Neutron Spherical Woods-Saxon Spectrum '',
     *               '' (LDSING='',I3,'')'',T63,
     *               ''Z='',I3,''   N='',I3,''    #'')') 
     *               LDSING_NEUTRS,IZ_FIX,IN_FIX
C
          WRITE(LOGNEU,'(  ''#'',78X,''#'',/,80(''#''),             /,
     *            ''#'',78X,''#'',/,
     *            ''#   No)     Energy     State   Occup   '',
     *            ''    No)     Energy     State   Occup    #'',/,
     *            ''#'',78X,''#'')')
      END IF
C
C=======================================================================
C
      DO ISTATE=1,LDSING_NEUTRS
         IF (ENETHE_NEUTRS(ISTATE).GT.ENEMAX) GO TO 2
      END DO
C
   2  CONTINUE
C   
      LDLINE=(ISTATE-1)/2
C
      DO LINE=1,LDLINE
C
         LINE_1=LINE
         LINE_2=LDLINE+LINE
C
         ICUM_1=ICUMUL_NEUTRS(LINE_1)
C
         ICUM_2=ICUMUL_NEUTRS(LINE_2)
C
         IF (LOGWRI.GT.4)
     *       WRITE(IRESUL,'(''#  '', I3,'')'',F12.5,4X,A6,I6,'' '',
     *                  ''      '',I3,'')'',F12.5,4X,A6,I6,'' '',
     *                                                 ''    #'')')
     *         LINE_1,ENETHE_NEUTRS(LINE_1),
     *                LABTHE_NEUTRS(LINE_1),ICUM_1,
     *         LINE_2,ENETHE_NEUTRS(LINE_2),
     *                LABTHE_NEUTRS(LINE_2),ICUM_2
C
          IF (IwMODE.EQ.2 .AND. LOGWRI.GT.4) THEN
C
              WRITE(LOGNEU,'(''#  '', I3,'')'',F12.5,4X,A6,I6,'' '',
     *                       ''      '',I3,'')'',F12.5,4X,A6,I6,'' '',
     *                                                    ''    #'')')
     *              LINE_1,ENETHE_NEUTRS(LINE_1),
     *                     LABTHE_NEUTRS(LINE_1),ICUM_1,
     *              LINE_2,ENETHE_NEUTRS(LINE_2),
     *                     LABTHE_NEUTRS(LINE_2),ICUM_2
          END IF
C
      END DO
C
C=======================================================================
C
      IF (LOGWRI.GT.4) WRITE(IRESUL,'(''#'',78X,''#'',/,80(''#''))')
C
      IF (IwMODE.EQ.2 .AND. LOGWRI.GT.4) THEN
          WRITE(LOGNEU,'(''#'',78X,''#'',/,80(''#''))')
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
      SUBROUTINE READEG(IDEGEN,LBSPHE)
      CHARACTER
     *          LBSPHE*6
C
C=======================================================================
C     This subroutine defines  the level degeneracies
C     using for this purpose the spectroscopic labels
C=======================================================================
C
      IF (LBSPHE(4:4).EQ.'/') THEN
          READ(LBSPHE(3:3),'(I1)') JDEGEN
      END IF
C
      IF (LBSPHE(5:5).EQ.'/') THEN
          READ(LBSPHE(3:4),'(I2)') JDEGEN
      END IF
C
      IDEGEN=JDEGEN+1
C
C=======================================================================
C
      RETURN
      END
C
C=======================================================================
C=======================================================================
C
      SUBROUTINE PRTNG(INUCLI,ISOSPI)  
          
      INCLUDE   'MATDIM/NDSPEC.f'
      INCLUDE   'MATDIM/NITERA.f'
      INCLUDE   'MATDIM/NDNUCL.f'
      INCLUDE   'MATDIM/NDIM_N.f'
      INCLUDE   'MATDIM/NDLEXP.f'
      
      PARAMETER
     *         (NDIM_L=NDIM_N)
      CHARACTER
     *          LABTHE_PROTON*6,LABTHE_NEUTRS*6,LABTHE*6
      CHARACTER
     *          LABEXP_PROTON*6,LABEXP_NEUTRS*6,LABEXP*6
      DIMENSION
     *          NWSSPH(1:NDSPEC),
     *          LWSSPH(1:NDSPEC),
     *          JWSSPH(1:NDSPEC),
     *          ENETHE(1:NDSPEC),
     *          EXPEXP(1:NDSPEC),
     *          LABTHE(1:NDSPEC),
     *          LABEXP(1:NDSPEC),
     *          LABORD(0:NDIM_N,0:NDIM_L,0:2*NDIM_L+1)
      COMMON
     *       /DIMSPH/ LDSING_PROTON,LEVTHE_PROTON,
     *                LEVEXP_PROTON(1:NDNUCL),
     *                LDSING_NEUTRS,LEVTHE_NEUTRS,
     *                LEVEXP_NEUTRS(1:NDNUCL)
      COMMON
     *       /THELAB/ LABTHE_PROTON(1:NDSPEC),
     *                LABTHE_NEUTRS(1:NDSPEC)
      COMMON
     *       /SPECTR/ ENETHE_PROTON(1:NDSPEC),
     *                ENETHE_NEUTRS(1:NDSPEC)
      COMMON
     *       /QUANTN/ NWSSPH_PROTON(1:NDSPEC),
     *                LWSSPH_PROTON(1:NDSPEC),
     *                JWSSPH_PROTON(1:NDSPEC),
     *
     *                NWSSPH_NEUTRS(1:NDSPEC),
     *                LWSSPH_NEUTRS(1:NDSPEC),
     *                JWSSPH_NEUTRS(1:NDSPEC)
      COMMON
     *       /LEVORD/ LABORD_PROTON(1:NDNUCL,
     *                              0:NDIM_N,0:NDIM_L,0:2*NDIM_L+1),
     *                LABORD_NEUTRS(1:NDNUCL,
     *                              0:NDIM_N,0:NDIM_L,0:2*NDIM_L+1)
      COMMON
     *       /EX_ENE/ EXPEXP_PROTON(1:NDNUCL,1:NDLEXP),
     *                IDEGEX_PROTON(1:NDNUCL,1:NDLEXP),
     *
     *                EXPEXP_NEUTRS(1:NDNUCL,1:NDLEXP),
     *                IDEGEX_NEUTRS(1:NDNUCL,1:NDLEXP)
      COMMON
     *       /EX_LAB/ LABEXP_PROTON(1:NDNUCL,1:NDLEXP),
     *                LABEXP_NEUTRS(1:NDNUCL,1:NDLEXP)
C
C=======================================================================
C     
      
      OPEN(70,FILE='Eth_vs_Eexp.d',STATUS='UNKNOWN')
      
      WRITE(70,'(/,''LDSING'',7X,''n'',4X,''L'',4X,''J'',3X,
     *             ''LEVEL'',4X,''Y?'',9X,
     *             ''E_THEO'',8X,''E_EXPE'',/)')
C
C=======================================================================
C
      IF (ISOSPI.EQ.1) THEN
        
        LDSING=LDSING_PROTON
        LEVEXP=LEVEXP_PROTON(INUCLI)
        
        DO I=1,LDSING
            
            NWSSPH(I)=NWSSPH_PROTON(I)
            LWSSPH(I)=LWSSPH_PROTON(I)
            JWSSPH(I)=JWSSPH_PROTON(I)
            
            ENETHE(I)=ENETHE_PROTON(I)
            LABTHE(I)=LABTHE_PROTON(I)
            
            LABORD(NWSSPH(I),LWSSPH(I),JWSSPH(I))=
     *      LABORD_PROTON(INUCLI,NWSSPH(I),LWSSPH(I),JWSSPH(I))
            
        END DO
        
        DO I=1,LEVEXP
            
            EXPEXP(I)=EXPEXP_PROTON(INUCLI,I)
            LABEXP(I)=LABEXP_PROTON(INUCLI,I)
            
        END DO
     
      ELSEIF (ISOSPI.EQ.0) THEN
      
        LDSING=LDSING_NEUTRS
        LEVEXP=LEVEXP_NEUTRS(INUCLI)
        
        DO I=1,LDSING
            
            NWSSPH(I)=NWSSPH_NEUTRS(I)
            LWSSPH(I)=LWSSPH_NEUTRS(I)
            JWSSPH(I)=JWSSPH_NEUTRS(I)
            
            ENETHE(I)=ENETHE_NEUTRS(I)
            LABTHE(I)=LABTHE_NEUTRS(I)
            
            LABORD(NWSSPH(I),LWSSPH(I),JWSSPH(I))=
     *      LABORD_NEUTRS(INUCLI,NWSSPH(I),LWSSPH(I),JWSSPH(I))
            
        END DO
       
        DO I=1,LEVEXP
            
            EXPEXP(I)=EXPEXP_PROTON(INUCLI,I)
            LABEXP(I)=LABEXP_PROTON(INUCLI,I)
            
         END DO
       
      END IF
C
C=======================================================================
C      
      DO I=1,LDSING
        
        IF (LABORD(NWSSPH(I),LWSSPH(I),JWSSPH(I)).EQ.1) THEN 
        
         DO J=1,LEVEXP
C      
          IF (LABEXP(J).EQ.LABTHE(I)) THEN
      
           WRITE(70,'(I3,6X,3I5,3x,a6,3x,''Y'',3X,2F14.6)')
     *       I,NWSSPH(I),LWSSPH(I),JWSSPH(I),LABTHE(I),
     *                              ENETHE(I),EXPEXP(J)
          
           GO TO 3
          
          END IF
          
          END DO
           
           WRITE(70,'(I3,6X,3I5,3x,a6,3x,''Y'',3X,F14.6)')
     *       I,NWSSPH(I),LWSSPH(I),JWSSPH(I),LABTHE(I),
     *                                        ENETHE(I)
        
        ELSE
        
         DO J=1,LEVEXP
C      
          IF (LABEXP(J).EQ.LABTHE(I)) THEN
        
           WRITE(70,'(I3,6X,3I5,3x,a6,3x,''N'',3X,2F14.6)')
     *       I,NWSSPH(I),LWSSPH(I),JWSSPH(I),LABTHE(I),
     *                             ENETHE(I),EXPEXP(J)
           GO TO 3
          
          END IF
         
         END DO
         
         WRITE(70,'(I3,6X,3I5,3x,a6,3x,''N'',3X,F14.6)')
     *       I,NWSSPH(I),LWSSPH(I),JWSSPH(I),LABTHE(I),
     *                                        ENETHE(I)
          
        END IF
        
   3    CONTINUE     
      
      END DO
      
      CLOSE(70)
      
      RETURN
      END    
C     
C=======================================================================
C=======================================================================
C 

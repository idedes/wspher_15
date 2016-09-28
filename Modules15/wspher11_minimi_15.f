C FILE NAME = wspher11_minimi_15.f ! Keep this symbol:    $ident@string$
C
C=======================================================================
C=======================================================================
C         INTERFACE: USER-DEFINED PROBLEM VS LEVENBERG PACKAGE
C=======================================================================
C=======================================================================
C
      SUBROUTINE LMMINI(ARGPAR,I_SEED,LDRAND,NEWSED,NDLAST)
C
      INCLUDE   'MATDIM/NDSPEC.f'
      INCLUDE   'MATDIM/NDPARS.f'
      INCLUDE   'MATDIM/NDNUCL.f'
      INCLUDE   'MATDIM/NDLEXP.f'
      INCLUDE   'MATDIM/NDIM_M.f'
      INCLUDE   'MATDIM/MAXFEV.f'
C
      EXTERNAL 
     *          FUNMIN
C      
      PARAMETER 
     *         (LENGTH=256,INMODE=1,NPRINT=-1)
C
      CHARACTER
     *          LABPRO*6,LABNEU*6,LABPRO_PRINTG*006,LABNEU_PRINTG*6
      CHARACTER
     *          DIRNAM*256,GENAME*256,STRING*126
      CHARACTER
     *          LABTHE_PROTON*6,LABTHE_NEUTRS*6,
     *          LABEXP_PROTON*6,LABEXP_NEUTRS*6
      CHARACTER
     *          ACTION*1,CALLED*6,TYPCHI*6
      CHARACTER
     *          FILNAM*256,FILNAM_ENDING*13,INPSYM*6,
     *          LABTHE*006,LABTHE_PRINTG*06,TEXLAM*20,
     *          FILNA2*256,NUCSYM*6,FITNUC*2,FILNAM_NUCFIT*30,
     *          FITNUC_AUXILI*30,FILNA3*256,VERSIO*3
C    
      DIMENSION 
     *          ARGPAR(1:NDPARS),FUNVEC(1:NDIM_M),
     *          SCALFC(1:NDPARS),QTFARR(1:NDPARS)
      DIMENSION
     *          WORKA1(1:NDPARS),WORKA2(1:NDPARS),
     *          WORKA3(1:NDPARS),WORKA4(1:NDIM_M)
      DIMENSION 
     *          FUNJAC(1:NDIM_M,1:NDPARS)
      DIMENSION 
     *          IPIVOT(1:NDPARS)
      DIMENSION
     *          ENETHE(1:NDSPEC),LABTHE(1:NDSPEC)
      DIMENSION
     *          ENEPRO_PRINTG(1:NDSPEC),LABPRO_PRINTG(1:NDSPEC),
     *          ENENEU_PRINTG(1:NDSPEC),LABNEU_PRINTG(1:NDSPEC),
     *          ENETHE_PRINTG(1:NDSPEC),LABTHE_PRINTG(1:NDSPEC),
     *                                  PARPOT_PRINTG(1:NDPARS)
      DIMENSION
     *          ENEPRO(1:NDSPEC),LABPRO(1:NDSPEC),
     *          ENENEU(1:NDSPEC),LABNEU(1:NDSPEC)
      DIMENSION
     *          FITNUC(1:NDNUCL),
     *          FITNUC_AUXILI(1:NDNUCL)
C
      COMMON
     *       /DIMSPH/ LDSING_PROTON,LEVTHE_PROTON,
     *                              LEVEXP_PROTON(1:NDNUCL),
     *                LDSING_NEUTRS,LEVTHE_NEUTRS,
     *                              LEVEXP_NEUTRS(1:NDNUCL)
      COMMON
     *       /CHI_PN/ CHISQU_PROTON,CHISQU_NEUTRS
      COMMON
     *       /CHICHO/ VALCHI
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
     *       /TENOPT/ ICENTT,ISORBT,ITENSR
     *       /IFPAIR/ IF_PAI,IFGPAI,IFDPAI
      COMMON
     *       /CHIGRD/ CHISQU_GRDNRM,
     *                CHISQU_GRADIE(1:NDPARS)
      COMMON
     *       /STOPIN/ HIAUXI(1:MAXFEV)
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
     *       /PRFILE/ ISCREN,LOGFIL,LOGCHP,LOGCHN,LOGHIP,LOGHIN,
     *                LOGPRO,LOGNEU,LOGWRI,LSCREN,IRESUL,ICONVE,
     *                NOUTPT,LOGAUX,LOGENE,LOGRAD,LOGGRD,LOGMSH,
     *                                                   LOGBIS    
      COMMON
     *       /TOLERS/ TOLERF,TOLERX,TOLERG,FACTOR  
      COMMON  
     *       /VERSIN/ VERSIO
C
      DATA 
     *     NpUNIT / 80 /, 
     *     NnUNIT / 81 / 
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
C     Levenberg-Marquardt minimisation routine ??
C=======================================================================
C
C     ILEVEL - Number of functions
C @@@ ILEVEL should be the running value changing from case to case?
C @@@ NDPARS according to conventions the MAXIMUM dimension
C     NDPARS - Number of parameters (not greater than ILEVEL)
C
C            chi=sum_m(FUNJAC_m)^2 
C
C     Testing case setting ILEVEL=1
C     NDIM_M leading dimension (ILEVEL)
C
C=======================================================================
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
C  @@@ IRENE? CHECK THE PRIVILEGED ROLE OF 20
             IF (IPARAM.LE.20 .AND. IFDENS.EQ.0) THEN
                 IFPROT=1
             END IF
C
             IF (IPARAM.GT.20 .AND. IFDENS.EQ.0) THEN
                 IFNEUT=1
             END IF
C
             IF (IPARAM.GE.51 .AND. IFDENS.EQ.0) THEN
                 IFPROT=0
                 IFNEUT=0
                 IFBOTH=1
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
C
      IF (IFDENS.EQ.0 .AND. IFNEUT.EQ.1) FILNAM_ENDING='N'
C          
      IF (IFDENS.EQ.1 .OR. (IFDENS.EQ.0 .AND. IFBOTH.EQ.1)) THEN
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
          WRITE(FILNAM,'(''LogFiles/IFDENS-0/'',A,''_'',A1,
     *                   ''_IFDENS-'',I1,''_IFTENS-'',I1,
     *                   ''_IF-RAD-'',I1,''_IF-ORD-'',I1,
     *                   ''_IF-RHO-'',I1,
     *                   ''_IFDEEP-'',I1,''_IFPRON-'',I1,
     *                   ''_IFK-VC-'',I1,''_IFK-VS-'',I1,
     *                   ''_IFK-RC-'',I1,''_IFK-RS-'',I1,
     *                   ''_IFK-AC-'',I1,''_IFK-AS-'',I1,
     *                   ''_LDRAND-'',I5.5,''_'',A3,''.log'')') 
     *
     *                      FILNAM_NUCFIT(1:I1),FILNAM_ENDING,
     *                      IFDENS,IFTENS,IF_RAD,IF_INV,IF_RHO,
     *                      IFDEEP,IFPRON,IFK_VC,IFK_VS,IFK_RC,
     *                      IFK_RS,IFK_AC,IFK_AS,LDRAND,VERSIO
C          
          IF (NUCACT.GT.1 .AND. NUCACT.LT.LDNUCL) THEN
C              
              WRITE(FILNA2,'(''LogFiles/IFDENS-0/'',A,''_'',
     *                                      A1,''_Energies'',
     *                       ''_IFDENS-'',I1,''_IFTENS-'',I1,
     *                       ''_IF-RAD-'',I1,''_IF-ORD-'',I1,
     *                       ''_IF-RHO-'',I1,
     *                       ''_IFDEEP-'',I1,''_IFPRON-'',I1,
     *                       ''_IFK-VC-'',I1,''_IFK-VS-'',I1,
     *                       ''_IFK-RC-'',I1,''_IFK-RS-'',I1,
     *                       ''_IFK-AC-'',I1,''_IFK-AS-'',I1,
     *                       ''_LDRAND-'',I5.5,''_'',A3,''.log'')')
     * 
     *                          FILNAM_NUCFIT(1:I1),FILNAM_ENDING,
     *                          IFDENS,IFTENS,IF_RAD,IF_INV,IF_RHO,
     *                          IFDEEP,IFPRON,IFK_VC,IFK_VS,IFK_RC,
     *                          IFK_RS,IFK_AC,IFK_AS,LDRAND,VERSIO
             
              WRITE(FILNA3,'(''LogFiles/IFDENS-0/'',A,''_'',
     *                                         A1,''_Radii'',
     *                       ''_IFDENS-'',I1,''_IFTENS-'',I1,
     *                       ''_IF-RAD-'',I1,''_IF-ORD-'',I1,
     *                       ''_IF-RHO-'',I1,
     *                       ''_IFDEEP-'',I1,''_IFPRON-'',I1,
     *                       ''_IFK-VC-'',I1,''_IFK-VS-'',I1,
     *                       ''_IFK-RC-'',I1,''_IFK-RS-'',I1,
     *                       ''_IFK-AC-'',I1,''_IFK-AS-'',I1,
     *                       ''_LDRAND-'',I5.5,''_'',A3,''.log'')')
     * 
     *                          FILNAM_NUCFIT(1:I1),FILNAM_ENDING,
     *                          IFDENS,IFTENS,IF_RAD,IF_INV,IF_RHO,
     *                          IFDEEP,IFPRON,IFK_VC,IFK_VS,IFK_RC,
     *                          IFK_RS,IFK_AC,IFK_AS,LDRAND,VERSIO
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
          WRITE(FILNAM,'(''LogFiles/IFDENS-0/'',A,''_'',A1,
     *                   ''_IFDENS-'',I1,''_IFTENS-'',I1,
     *                   ''_IF-PAI-'',I1,
     *                   ''_IF-RAD-'',I1,''_IF-ORD-'',I1,
     *                   ''_IF-RHO-'',I1,
     *                   ''_IFDEEP-'',I1,''_IFPRON-'',I1,
     *                   ''_IFK-VC-'',I1,''_IFK-VS-'',I1,
     *                   ''_IFK-RC-'',I1,''_IFK-RS-'',I1,
     *                   ''_IFK-AC-'',I1,''_IFK-AS-'',I1,
     *                   ''_LDRAND-'',I5.5,''_'',A3,''.log'')') 
     *
     *                          FILNAM_NUCFIT(1:I1),FILNAM_ENDING,
     *                          IFDENS,IFTENS,IF_PAI,IF_RAD,IF_INV,
     *                          IF_RHO,IFDEEP,IFPRON,IFK_VC,IFK_VS,
     *                          IFK_RC,IFK_RS,IFK_AC,IFK_AS,LDRAND,
     *                                                      VERSIO
C          
          IF (NUCACT.GT.1 .AND. NUCACT.LT.LDNUCL) THEN
C              
              WRITE(FILNA2,'(''LogFiles/IFDENS-0/'',A,''_'',
     *                                      A1,''_Energies'',
     *                       ''_IFDENS-'',I1,''_IFTENS-'',I1,
     *                       ''_IF-PAI-'',I1,
     *                       ''_IF-RAD-'',I1,''_IF-ORD-'',I1,
     *                       ''_IF-RHO-'',I1,
     *                       ''_IFDEEP-'',I1,''_IFPRON-'',I1,
     *                       ''_IFK-VC-'',I1,''_IFK-VS-'',I1,
     *                       ''_IFK-RC-'',I1,''_IFK-RS-'',I1,
     *                       ''_IFK-AC-'',I1,''_IFK-AS-'',I1,
     *                       ''_LDRAND-'',I5.5,''_'',A3,''.log'')')
     * 
     *                          FILNAM_NUCFIT(1:I1),FILNAM_ENDING,
     *                          IFDENS,IFTENS,IF_PAI,IF_RAD,IF_INV,
     *                          IF_RHO,IFDEEP,IFPRON,IFK_VC,IFK_VS,
     *                          IFK_RC,IFK_RS,IFK_AC,IFK_AS,LDRAND,
     *                                                      VERSIO
             
              WRITE(FILNA3,'(''LogFiles/IFDENS-0/'',A,''_'',
     *                                         A1,''_Radii'',
     *                       ''_IFDENS-'',I1,''_IFTENS-'',I1,
     *                       ''_IF-PAI-'',I1,
     *                       ''_IF-RAD-'',I1,''_IF-ORD-'',I1,
     *                       ''_IF-RHO-'',I1,
     *                       ''_IFDEEP-'',I1,''_IFPRON-'',I1,
     *                       ''_IFK-VC-'',I1,''_IFK-VS-'',I1,
     *                       ''_IFK-RC-'',I1,''_IFK-RS-'',I1,
     *                       ''_IFK-AC-'',I1,''_IFK-AS-'',I1,
     *                       ''_LDRAND-'',I5.5,''_'',A3,''.log'')')
     * 
     *                          FILNAM_NUCFIT(1:I1),FILNAM_ENDING,
     *                          IFDENS,IFTENS,IF_PAI,IF_RAD,IF_INV,
     *                          IF_RHO,IFDEEP,IFPRON,IFK_VC,IFK_VS,
     *                          IFK_RC,IFK_RS,IFK_AC,IFK_AS,LDRAND,
     *                                                      VERSIO
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
          WRITE(FILNAM,'(''LogFiles/IFDENS-1_IFTENS-'',I1,''/'',
     *                                        A,''_'',A1,
     *                   ''_IFDENS-'',I1,''_IFTENS-'',I1,
     *                   ''_IF-PAI-'',I1,
     *                   ''_IF-RAD-'',I1,''_IF-ORD-'',I1,
     *                   ''_IFDEEP-'',I1,''_IFPRON-'',I1,
     *                   ''_IFK-VC-'',I1,''_IFK-VS-'',I1,
     *                   ''_IFK-RC-'',I1,''_IFK-RS-'',I1,
     *                   ''_IFK-AC-'',I1,''_IFK-AS-'',I1,
     *                   ''_LDRAND-'',I5.5,A15,''_'',A3,''.log'')') 
     *
     *                      IFTENS,FILNAM_NUCFIT(1:I1),FILNAM_ENDING,
     *                      IFDENS,IFTENS,IF_PAI,IF_RAD,IF_INV,
     *                      IFDEEP,IFPRON,IFK_VC,IFK_VS,
     *                      IFK_RC,IFK_RS,IFK_AC,IFK_AS,LDRAND,TEXLAM,
     *                                                         VERSIO
C______________________________________________________________________
C
          IF (IFTENS.EQ.1) THEN ! IFDENS=1
C
              WRITE(FILNAM,'(''LogFiles/IFDENS-1_IFTENS-'',I1,''/'',
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
     *                   ''_LDRAND-'',I5.5,A15,''_'',A3,''.log'')') 
     *
     *                      IFTENS,FILNAM_NUCFIT(1:I1),FILNAM_ENDING,
     *                      IFDENS,IFTENS,IF_PAI,ICENTT,ISORBT,ITENSR,
     *                      IF_RAD,IF_INV,IFDEEP,IFPRON,IFK_VC,IFK_VS,
     *                      IFK_RC,IFK_RS,IFK_AC,IFK_AS,LDRAND,TEXLAM,
     *                                                         VERSIO
          END IF
C______________________________________________________________________
C
          IF (NUCACT.GT.1 .AND. NUCACT.LT.LDNUCL) THEN ! IFDENS=1
C
              WRITE(FILNA2,'(''LogFiles/IFDENS-1_IFTENS-'',I1,''/'',
     *                              A,''_'',A1,''_Energies'',
     *                       ''_IFDENS-'',I1,''_IFTENS-'',I1,
     *                       ''_IF-PAI-'',I1,
     *                       ''_IF-RAD-'',I1,''_IF-ORD-'',I1,
     *                       ''_IFDEEP-'',I1,''_IFPRON-'',I1,
     *                       ''_IFK-VC-'',I1,''_IFK-VS-'',I1,
     *                       ''_IFK-RC-'',I1,''_IFK-RS-'',I1,
     *                       ''_IFK-AC-'',I1,''_IFK-AS-'',I1,
     *                       ''_LDRAND-'',I5.5,A15,''_'',A3,''.log'')') 
     *
     *                      IFTENS,FILNAM_NUCFIT(1:I1),FILNAM_ENDING,
     *                      IFDENS,IFTENS,IF_PAI,
     *                      IF_RAD,IF_INV,IFDEEP,IFPRON,IFK_VC,IFK_VS,
     *                      IFK_RC,IFK_RS,IFK_AC,IFK_AS,LDRAND,TEXLAM,
     *                                                         VERSIO
C______________________________________________________________________
C
              WRITE(FILNA3,'(''LogFiles/IFDENS-1_IFTENS-'',I1,''/'',
     *                                 A,''_'',A1,''_Radii'',
     *                       ''_IFDENS-'',I1,''_IFTENS-'',I1,
     *                       ''_IF-PAI-'',I1,
     *                       ''_IF-RAD-'',I1,''_IF-ORD-'',I1,
     *                       ''_IFDEEP-'',I1,''_IFPRON-'',I1,
     *                       ''_IFK-VC-'',I1,''_IFK-VS-'',I1,
     *                       ''_IFK-RC-'',I1,''_IFK-RS-'',I1,
     *                       ''_IFK-AC-'',I1,''_IFK-AS-'',I1,
     *                       ''_LDRAND-'',I5.5,A15,''_'',A3,''.log'')') 
     * 
     *                      IFTENS,FILNAM_NUCFIT(1:I1),FILNAM_ENDING,
     *                      IFDENS,IFTENS,IF_PAI,
     *                      IF_RAD,IF_INV,IFDEEP,IFPRON,IFK_VC,IFK_VS,
     *                      IFK_RC,IFK_RS,IFK_AC,IFK_AS,LDRAND,TEXLAM,
     *                                                         VERSIO
C______________________________________________________________________
C
              IF (IFTENS.EQ.1) THEN ! IFDENS=1
C
                  WRITE(FILNA2,'(''LogFiles/IFDENS-1_IFTENS-'',I1,''/'',
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
     *                           ''_LDRAND-'',I5.5,A15,''_'',A3,
     *                                                 ''.log'')') 
     *
     *                      IFTENS,FILNAM_NUCFIT(1:I1),FILNAM_ENDING,
     *                      IFDENS,IFTENS,IF_PAI,ICENTT,ISORBT,ITENSR,
     *                      IF_RAD,IF_INV,IFDEEP,IFPRON,IFK_VC,IFK_VS,
     *                      IFK_RC,IFK_RS,IFK_AC,IFK_AS,LDRAND,TEXLAM,
     *                                                         VERSIO
C
C______________________________________________________________________
C
                  WRITE(FILNA3,'(''LogFiles/IFDENS-1_IFTENS-'',I1,''/'',
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
     *                           ''_LDRAND-'',I5.5,A15,''_'',A3,
     *                                                 ''.log'')') 
     *
     *                      IFTENS,FILNAM_NUCFIT(1:I1),FILNAM_ENDING,
     *                      IFDENS,IFTENS,IF_PAI,ICENTT,ISORBT,ITENSR,
     *                      IF_RAD,IF_INV,IFDEEP,IFPRON,IFK_VC,IFK_VS,
     *                      IFK_RC,IFK_RS,IFK_AC,IFK_AS,LDRAND,TEXLAM,
     *                                                         VERSIO
C
              END IF
C          
              OPEN(LOGENE,FILE=FILNA2,STATUS='UNKNOWN',FORM='FORMATTED')
              OPEN(LOGRAD,FILE=FILNA3,STATUS='UNKNOWN',FORM='FORMATTED')
C
          END IF
C      
      END IF
C
C=======================================================================
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C=======================================================================
C @@@ IRENE WHAT IS THIS CONDITION EXACTLY? WHAT IS THE ROLE OF NUCACT?     
      IF (NUCACT.EQ.LDNUCL) THEN 
C
          IF (IFDENS.EQ.0) THEN
C
              WRITE(FILNAM,'(''LogFiles/IFDENS-0/All-Nuclei_'',A1,
     *                       ''_IFDENS-'',I1,''_IFTENS-'',I1,
     *                       ''_IF-PAI-'',I1,
     *                       ''_IF-RAD-'',I1,''_IF-ORD-'',I1,
     *                       ''_IFDEEP-'',I1,''_IFPRON-'',I1,
     *                       ''_IFK-VC-'',I1,''_IFK-VS-'',I1,
     *                       ''_IFK-RC-'',I1,''_IFK-RS-'',I1,
     *                       ''_IFK-AC-'',I1,''_IFK-AS-'',I1,
     *                       ''_LDRAND-'',I5.5,''_'',A3,''.log'')')
     *
     *                          FILNAM_ENDING,IFDENS,IFTENS,IF_PAI,
     *                          IF_RAD,IF_INV,IFDEEP,IFPRON,IFK_VC,
     *                          IFK_VS,IFK_RC,IFK_RS,IFK_AC,IFK_AS,
     *                                               LDRAND,VERSIO
C
              WRITE(FILNA2,'(''LogFiles/IFDENS-0/All-Nuclei_Energies_'',
     *                    A1,''_IFDENS-'',I1,''_IFTENS-'',I1,
     *                       ''_IF-PAI-'',I1,
     *                       ''_IF-RAD-'',I1,''_IF-ORD-'',I1,
     *                       ''_IFDEEP-'',I1,''_IFPRON-'',I1,
     *                       ''_IFK-VC-'',I1,''_IFK-VS-'',I1,
     *                       ''_IFK-RC-'',I1,''_IFK-RS-'',I1,
     *                       ''_IFK-AC-'',I1,''_IFK-AS-'',I1,
     *                       ''_LDRAND-'',I5.5,''_'',A3,''.log'')')
     *
     *                          FILNAM_ENDING,IFDENS,IFTENS,IF_PAI,
     *                          IF_RAD,IF_INV,IFDEEP,IFPRON,IFK_VC,
     *                          IFK_VS,IFK_RC,IFK_RS,IFK_AC,IFK_AS,
     *                                               LDRAND,VERSIO
C
              WRITE(FILNA3,'(''LogFiles/IFDENS-0/All-Nuclei_Radii_'',A1,
     *                       ''_IFDENS-'',I1,''_IFTENS-'',I1,
     *                       ''_IF-PAI-'',I1,
     *                       ''_IF-RAD-'',I1,''_IF-ORD-'',I1,
     *                       ''_IFDEEP-'',I1,''_IFPRON-'',I1,
     *                       ''_IFK-VC-'',I1,''_IFK-VS-'',I1,
     *                       ''_IFK-RC-'',I1,''_IFK-RS-'',I1,
     *                       ''_IFK-AC-'',I1,''_IFK-AS-'',I1,
     *                       ''_LDRAND-'',I5.5,''_'',A3,''.log'')')
     *
     *                          FILNAM_ENDING,IFDENS,IFTENS,IF_PAI,
     *                          IF_RAD,IF_INV,IFDEEP,IFPRON,IFK_VC,
     *                          IFK_VS,IFK_RC,IFK_RS,IFK_AC,IFK_AS,
     *                                               LDRAND,VERSIO
C              
          END IF
C
C=======================================================================
C          
          IF (IFDENS.EQ.1) THEN 
C              
              WRITE(FILNAM,'(''LogFiles/IFDENS-1_IFTENS-'',I1,
     *                       ''/All-Nuclei_'',A1,
     *                       ''_IFDENS-'',I1,''_IFTENS-'',I1,
     *                       ''_IF-PAI-'',I1,
     *                       ''_ICENTT-'',I1,''_ISORBT-'',I1,
     *                                       ''_ITENSR-'',I1,
     *                       ''_IF-RAD-'',I1,''_IF-ORD-'',I1,
     *                       ''_IFDEEP-'',I1,''_IFPRON-'',I1,
     *                       ''_IFK-VC-'',I1,''_IFK-VS-'',I1,
     *                       ''_IFK-RC-'',I1,''_IFK-RS-'',I1,
     *                       ''_IFK-AC-'',I1,''_IFK-AS-'',I1,
     *                       ''_LDRAND-'',I5.5,A15,''_'',A3,''.log'')')
     *
     *                          IFTENS,FILNAM_ENDING,IFDENS,IFTENS,
     *                          IF_PAI,ICENTT,ISORBT,ITENSR,IF_RAD,
     *                          IF_INV,IFDEEP,IFPRON,IFK_VC,IFK_VS,
     *                          IFK_RC,IFK_RS,IFK_AC,IFK_AS,LDRAND,
     *                                               TEXLAM,VERSIO
C
               WRITE(FILNA2,'(''LogFiles/IFDENS-1_IFTENS-'',I1,
     *                        ''/All-Nuclei_Energies_'',A1,
     *                        ''_IFDENS-'',I1,''_IFTENS-'',I1,
     *                        ''_IF-PAI-'',I1,
     *                        ''_ICENTT-'',I1,''_ISORBT-'',I1,
     *                                       ''_ITENSR-'',I1,
     *                        ''_IF-RAD-'',I1,''_IF-ORD-'',I1,
     *                        ''_IFDEEP-'',I1,''_IFPRON-'',I1,
     *                        ''_IFK-VC-'',I1,''_IFK-VS-'',I1,
     *                        ''_IFK-RC-'',I1,''_IFK-RS-'',I1,
     *                        ''_IFK-AC-'',I1,''_IFK-AS-'',I1,
     *                        ''_LDRAND-'',I5.5,A15,''_'',A3,''.log'')')
     *
     *                          IFTENS,FILNAM_ENDING,IFDENS,IFTENS,
     *                          IF_PAI,ICENTT,ISORBT,ITENSR,IF_RAD,
     *                          IF_INV,IFDEEP,IFPRON,IFK_VC,IFK_VS,
     *                          IFK_RC,IFK_RS,IFK_AC,IFK_AS,LDRAND,
     *                                               TEXLAM,VERSIO
C
               WRITE(FILNA3,'(''LogFiles/IFDENS-1_IFTENS-'',I1,
     *                        ''/All-Nuclei_Radii_'',A1,
     *                        ''_IFDENS-'',I1,''_IFTENS-'',I1,
     *                        ''_IF-PAI-'',I1,
     *                        ''_ICENTT-'',I1,''_ISORBT-'',I1,
     *                                       ''_ITENSR-'',I1,
     *                        ''_IF-RAD-'',I1,''_IF-ORD-'',I1,
     *                        ''_IFDEEP-'',I1,''_IFPRON-'',I1,
     *                        ''_IFK-VC-'',I1,''_IFK-VS-'',I1,
     *                        ''_IFK-RC-'',I1,''_IFK-RS-'',I1,
     *                        ''_IFK-AC-'',I1,''_IFK-AS-'',I1,
     *                        ''_LDRAND-'',I5.5,A15,''_'',A3,''.log'')')
     *
     *                          IFTENS,FILNAM_ENDING,IFDENS,IFTENS,
     *                          IF_PAI,ICENTT,ISORBT,ITENSR,IF_RAD,
     *                          IF_INV,IFDEEP,IFPRON,IFK_VC,IFK_VS,
     *                          IFK_RC,IFK_RS,IFK_AC,IFK_AS,LDRAND,
     *                                               TEXLAM,VERSIO
C
          END IF   
C
          OPEN(LOGENE,FILE=FILNA2,STATUS='UNKNOWN',FORM='FORMATTED')
          OPEN(LOGRAD,FILE=FILNA3,STATUS='UNKNOWN',FORM='FORMATTED')
C                    
      END IF
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
      NUMPAR=IACTIV
      INFOER=111
C
C=======================================================================
C    
      CHISQU_MINIML=1.0E+10
      
      IMASIV_PRODUC=0
      
      DO IRANDO=1,LDRAND
          
         IFUNCT_EVALUS=0
C     
         IDEFCN=1
C
         DO I=1,MAXFEV
             HIAUXI(I)=0.0
         END DO
C         
         IF (ISCREN.EQ.1) THEN
C
             WRITE(LSCREN,'(/,58X,''New run: NEWSED='',i14,3x,
     *                            ''No.'',i3,''/'',i3,/,101(''=''),/)') 
     *                                       NEWSED,IRANDO,LDRAND
         END IF
C         
         IF (LOGWRI.GT.0) THEN
C
             WRITE(LOGAUX,'(/,58X,''New run: NEWSED='',i14,3x,
     *                            ''No.'',i3,''/'',i3,/,101(''=''),/)') 
     *                                       NEWSED,IRANDO,LDRAND
C             
             WRITE(LOGENE,'(/,58X,''New run: NEWSED='',i14,3x,
     *                            ''No.'',i3,''/'',i3,/,101(''=''),/)') 
     *                                       NEWSED,IRANDO,LDRAND
C             
             WRITE(LOGRAD,'(/,58X,''New run: NEWSED='',i14,3x,
     *                            ''No.'',i3,''/'',i3,/,101(''=''),/)') 
     *                                       NEWSED,IRANDO,LDRAND
         END IF
C
         WRITE(NOUTPT,'(/,80(''#''),/,''#'',T80,''#'',/,
     *                  ''#'',5X,''IRANDO= '',I4,T80,''#'',/,
     *                  ''#'',T80,''#'',/,80(''#''))') IRANDO
C      
C=======================================================================
C
C        Generating random restarts
C
         IF (LOGWRI.GT.0) THEN
             WRITE(LOGFIL,'(//,9X,''Entering NEWRAN with IRANDO='',I3,
     *                         1X,''LDRAND='',I6,'' from LMMINI'')')
     *                              IRANDO,LDRAND
         END IF
C         
         CALL NEWRAN(NEWSED,IZ_FIX,IN_FIX,ISOSPI,ARGPAR,DIRNAM,
     *                                    GENAME,IFIRST,I_SEED)
C      
C=======================================================================
C        Calling the principal minimisation routine based on the
C        Levenberg-Marquardt algorithm (Argonne MINPACK project)  
C=======================================================================
C
         NDFUNC=NDIM_M    ! NDIM_M
         LDFUNC=LEVNUM    ! Actual number of functions, here: levels
C IRENE? NDPARS=NDPARS=48 ! ARGPAR(1:NDPARS)
         LDPARS=IACTIV
C 
C        MAXFEV=900       ! MAXFEV
C @@@ ??
         IFMODE=INMODE
C        NPRINT=-1        ! NPRINT
C @@@ WHAT DOES THIS MEAN???
         INFRUN=INFOER    ! In any case output variable
         NFCALL=0         ! NFUNEV    ! In any case output variable
         NJCALL=0         ! NJACEV    ! In any case output variable
C
C        I_PERM=IPIVOT    ! In any case output vector ?????????????
C                     
C        FJACOB(1:NDFUNC,1:NDPARS) <<=>> FUNJAC(1:NDIM_M,1:NDPARS)
C                 FARGUN(1:NDFUNC) <<=>> FUNVEC(1:NDIM_M)
C                 DIAGSC(1:NDPARS) <<=>> SCALFC(1:NDPARS)       
C                 I_PERM(1:NDPARS) <<=>> IPIVOT(1:NDPARS)       
C
C        Calling the Levenberg-Marquardt LAPACK minimisation routine 
C
         IF (LOGWRI.GT.0) THEN
             WRITE(LOGFIL,'(/,6X,''Entering LEVMAR from LMMINI      '',
     *                        5(3x,''IRANDO='',I2))')
     *                               IRANDO,IRANDO,IRANDO,IRANDO,IRANDO
         END IF
C
         CALL LEVMAR(FUNMIN,NDFUNC,LDFUNC,NDPARS,LDPARS,ARGPAR,
     *               FUNVEC,FUNJAC,TOLERF,TOLERX,TOLERG,MAXFEV,
     *               SCALFC,IFMODE,FACTOR,NPRINT,INFRUN,NFCALL,
     *               NJCALL,IPIVOT,QTFARR,WORKA1,WORKA2,WORKA3,
     *                                           WORKA4,NDLAST)
C
         INFOER=INFRUN 
C
         KK=0
         DO KPARAM=1,NDPARS
            IF (IFTAKE(KPARAM).EQ.1) THEN
                KK=KK+1
                PARPOT(KPARAM)=ARGPAR(KK)
            END IF
         END DO     
C
C=======================================================================
C
         WRITE(LOGAUX,'(/,''Number of WS_RUN Evaluations: '',I5,/)')
     *                      IFUNCT_EVALUS
C
         WRITE(LOGENE,'()') ! AESTHETICAL REASONS
         WRITE(LOGRAD,'()') ! AESTHETICAL REASONS
C
C=======================================================================
C      
         CALL TELLIT(INFOER,TOLERF,TOLERX,TOLERG,MAXFEV,LDLAST,EPSLAS)
C
C=======================================================================
C
         IF (LOGWRI.GT.0) THEN
             WRITE(LOGFIL,'(9X,''Entering JMATRX from LMMINI'')')
         END IF
C
         CALL JMATRX(LDFUNC,ARGPAR)
C
C=======================================================================
C  
C        Between the different restarts, we look for the one that
C        gives the smallest \chi^2
C
         IF (IFDENS.EQ.0 .AND. IFPROT.EQ.1) CHISQU=CHISQU_PROTON
C        
         IF (IFDENS.EQ.0 .AND. IFNEUT.EQ.1) CHISQU=CHISQU_NEUTRS
C
         IF (IFDENS.EQ.1 .OR. (IFDENS.EQ.0 .AND. IFBOTH.EQ.1)) THEN
             CHISQU=VALCHI
         END IF
C_______________________________________________________________________
C         
         IF (CHISQU.LT.CHISQU_MINIML) THEN
C
             CHISQU_MINIML=CHISQU
C
             IRANDO_PRINTG=IRANDO
             IEVALU_PRINTG=IFUNCT_EVALUS
             CHIGRD_PRINTG=CHISQU_GRDNRM
C
             DO IPARAM=1,NDPARS
                PARPOT_PRINTG(IPARAM)=PARPOT(IPARAM)
             END DO
C
         END IF
C
C=======================================================================
C
         IF (LOGWRI.GT.0) THEN
             WRITE(LOGFIL,'(9X,''Entering RADINF from LMMINI'')')
         END IF
C
         CALL RADINF(IFPROT,IFNEUT)
C
C=======================================================================
C
         I_SEED=NEWSED
C
      END DO  !IRANDO
C
C=======================================================================
      
      IMASIV_PRODUC=1
C
C=======================================================================
C
      CALL FINISH_LMMINI(IRANDO_PRINTG,LDRAND,PARPOT_PRINTG,
     *                   IEVALU_PRINTG,CHISQU_MINIML,CHIGRD_PRINTG)
C      
C=======================================================================
C 
      REWIND LOGPRO    
C      
C=======================================================================
C
   10 CONTINUE
C @@@ WHAT IS THIS?   
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
C
   20 CONTINUE
C @@@ IDEM  
      READ (LOGNEU,'(A120)',END=21) STRING
      IF (LOGWRI.GT.4) WRITE(NnUNIT,'(A80)'        ) STRING
C      
      GO TO 20
C
   21 CONTINUE  ! Ready with the copy for the neutrons     
C
C=======================================================================
C
      IF (LOGWRI.GT.0) THEN
          WRITE(LOGFIL,'(/,18X,''Exiting  LMMINI'')')
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
      SUBROUTINE FINISH_LMMINI(IRANDO_PRINTG,LDRAND,PARPOT_PRINTG,
     *                                IEVALU_PRINTG,CHISQU_MINIML,
     *                                              CHIGRD_PRINTG)
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
    
      EXTERNAL
     *         DENSIT,DENTHE_INTEGR
C      
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
     *          TITLES_LATEXS*100,NUCNAM_LATEXS*010
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
     *       /DIMSPH/ LDSING_PROTON,LEVTHE_PROTON,
     *                LEVEXP_PROTON(1:NDNUCL),
     *                LDSING_NEUTRS,LEVTHE_NEUTRS,
     *                LEVEXP_NEUTRS(1:NDNUCL)
      COMMON
     *       /MODIFR/ RMSEXP_PROTON(1:NDNUCL),
     *                RMSEXP_NEUTRS(1:NDNUCL),
     *                RMSTHE_PROTON,RMSTHE_NEUTRS,
     *                RMODIF_PROTON,RMODIF_NEUTRS,
     *                              PUSHIN,DELTAR
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
     *       /SPECTR/ ENETHE_PROTON(1:NDSPEC),
     *                ENETHE_NEUTRS(1:NDSPEC)
      COMMON
     *       /THELAB/ LABTHE_PROTON(1:NDSPEC),
     *                LABTHE_NEUTRS(1:NDSPEC)
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
     *       /P_N_LM/ IFPROT,IFNEUT,IFBOTH
      COMMON
     *       /DENSTR/ IFDENS
      COMMON
     *       /IFPAIR/ IF_PAI,IFGPAI,IFDPAI
      COMMON
     *       /FITTED/ IFITED
      COMMON
     *       /NUCLEU/ IZ_FIX,IN_FIX,ISOSPI
      COMMON
     *       /SIMUSE/ AOSCIL,RADCUT,INUCLI
      COMMON
     *       /CUTOFF/ RADIUS_CUTOFF(1:NDNUCL)
      COMMON
     *       /ARGUMN/ IARG_R,IARG_Z
      COMMON
     *       /PARAMT_TITLES/ TITLES(1:NDTITL),
     *                       TITLES_LATEXS(1:NDTITL),
     *                       NUCNAM_LATEXS(1:NDNUCL)
      COMMON
     *       /RMSGLO_NUCLEU/ RMSGLO_TAKPRO,
     *                       RMSGLO_TAKNEU,
     *                       RMSVAL_TAKPRO(1:NDNUCL),
     *                       RMSVAL_TAKNEU(1:NDNUCL)
      COMMON
     *       /PRFILE/ ISCREN,LOGFIL,LOGCHP,LOGCHN,LOGHIP,LOGHIN,
     *                LOGPRO,LOGNEU,LOGWRI,LSCREN,IRESUL,ICONVE,
     *                NOUTPT,LOGAUX,LOGENE,LOGRAD,LOGGRD,LOGMSH,
     *                                                   LOGBIS   
C
C=======================================================================
C
C     This routine calculates the energy levels for ALL the nuclei
C     with the optimal parameters resulting from the minimisation
C
C     The results are divided between "Fitted" & "Predic" folders
C
C=======================================================================
C  
      IPRINT_LMMINI=1
      
      DO IPARAM=1,NDPARS
         PARPOT(IPARAM)=PARPOT_PRINTG(IPARAM)
      END DO
      
      IF (NUCACT.LT.LDNUCL) RMSGLO_PRINTG=9.9999
          
      RMSGLO_PROTON=0.
      RMSGLO_NEUTRS=0.
C
      RMSGLO_PROAUX=0.
      RMSGLO_NEUAUX=0.
C
      RMSGLO_TAKPRO=0.
      RMSGLO_TAKNEU=0.
C
      IDEG_P=0
      IDEG_N=0
C
      IDEGEN_TAKPRO=0
      IDEGEN_TAKNEU=0
C
C=======================================================================
C      
      DO INUCLI=1,LDNUCL
          
         IZ_FIX=NUMB_Z(INUCLI)
         IN_FIX=NUMB_N(INUCLI)
C
         IF (LOGWRI.GT.4) THEN
             WRITE(LOGFIL,'(12X,''Entering WS_RUN from FINISH_LMMINI '',
     *                          ''INUCLI='',I2)') INUCLI
         END IF
                                                 I_MODE=0
         CALL WS_RUN(INUCLI,CHISQU_CHOICE,IABORT,I_MODE,I_FLAG,
     *                             CHISQU_AUXIL1,CHISQU_AUXIL2)
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
C
         IF (ITAKNU(INUCLI).EQ.1) THEN
C
             RMSGLO_TAKPRO=RMSGLO_TAKPRO+CHIDEG_PROTON
             RMSGLO_TAKNEU=RMSGLO_TAKNEU+CHIDEG_NEUTRS
C
             RMSVAL_TAKPRO(INUCLI)=RMSVAL_PRONUC(INUCLI)
             RMSVAL_TAKNEU(INUCLI)=RMSVAL_NEUNUC(INUCLI)
C
             DO IEXPER=1,LEVEXP_PROTON(INUCLI)
                IDEGEN_TAKPRO=IDEGEN_TAKPRO+IDEGEX_PROTON(INUCLI,IEXPER)
             END DO
C
             DO IEXPER=1,LEVEXP_NEUTRS(INUCLI)
                IDEGEN_TAKNEU=IDEGEN_TAKNEU+IDEGEX_NEUTRS(INUCLI,IEXPER)
             END DO
C
         END IF
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
      RMSGLO_TAKPRO=SQRT(RMSGLO_TAKPRO/IDEGEN_TAKPRO)
      RMSGLO_TAKNEU=SQRT(RMSGLO_TAKNEU/IDEGEN_TAKNEU)
C
C=======================================================================
C
C     Re-naming to have the right value for the lambdas
C     
      DO IPARAM=1,NDPARS
         PARPOT_PRINTG(IPARAM)=PARPOT(IPARAM)
      END DO
C
C=======================================================================
C     
      DO INUCLI=1,LDNUCL
          
         IZ_FIX=NUMB_Z(INUCLI)
         IN_FIX=NUMB_N(INUCLI)

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
                                  IFITED=0
         IF (ITAKNU(INUCLI).EQ.1) IFITED=1
             
         IF (IFDENS.EQ.0 .AND. IFPROT.EQ.1) THEN
         
             ISOSPI=1
             
             ERRMAX_PRINTG=ERRMAX_PRONUC(INUCLI)
C             EABSAV_PRINTG=EABSAV_PRONUC(INUCLI)
             RMSVAL_PRINTG=RMSVAL_PRONUC(INUCLI)
             
             RMSEXP_PRINTG=RMSEXP_PROTON(INUCLI)
             RMSTHE_PRINTG=RMSTHE_PRONUC(INUCLI)
             LEVTHE_PRINTG=LEVTHE_PRONUC(INUCLI)
             
             DO ITHEOR=1,LEVTHE_PRINTG
                ENETHE_PRINTG(ITHEOR)=ENETHE_PRONUC(ITHEOR,INUCLI)
                LABTHE_PRINTG(ITHEOR)=LABTHE_PRONUC(ITHEOR,INUCLI)
             END DO
             
             CALL WRITIN_ENELEV(ISOSPI,IRANDO_PRINTG,LDRAND,
     *                                        IEVALU_PRINTG,
     *                          PARPOT_PRINTG,CHISQU_MINIML,
     *                          CHIGRD_PRINTG,ERRMAX_PRINTG,
     *                          RMSVAL_PRINTG,RMSGLO_PROAUX,
     *                          RMSEXP_PRINTG,RMSTHE_PRINTG,
     *                          LEVTHE_PRINTG,ENETHE_PRINTG,
     *                          LABTHE_PRINTG,INUCLI)
         END IF
C_______________________________________________________________________
C         
         IF (IFDENS.EQ.0 .AND. IFNEUT.EQ.1) THEN
         
             ISOSPI=0
             
             ERRMAX_PRINTG=ERRMAX_NEUNUC(INUCLI)
C             EABSAV_PRINTG=EABSAV_NEUNUC(INUCLI)
             RMSVAL_PRINTG=RMSVAL_NEUNUC(INUCLI)
             
             RMSEXP_PRINTG=RMSEXP_NEUTRS(INUCLI)
             RMSTHE_PRINTG=RMSTHE_NEUNUC(INUCLI)
             LEVTHE_PRINTG=LEVTHE_NEUNUC(INUCLI)
             
             DO ITHEOR=1,LEVTHE_PRINTG
                ENETHE_PRINTG(ITHEOR)=ENETHE_NEUNUC(ITHEOR,INUCLI)
                LABTHE_PRINTG(ITHEOR)=LABTHE_NEUNUC(ITHEOR,INUCLI)
             END DO
             
             CALL WRITIN_ENELEV(ISOSPI,IRANDO_PRINTG,LDRAND,
     *                                        IEVALU_PRINTG,
     *                          PARPOT_PRINTG,CHISQU_MINIML,
     *                          CHIGRD_PRINTG,ERRMAX_PRINTG,
     *                          RMSVAL_PRINTG,RMSGLO_NEUAUX,
     *                          RMSEXP_PRINTG,RMSTHE_PRINTG,
     *                          LEVTHE_PRINTG,ENETHE_PRINTG,
     *                          LABTHE_PRINTG,INUCLI)
         END IF
C_______________________________________________________________________
C         
         IF (IFDENS.EQ.1 .OR. (IFDENS.EQ.0 .AND. IFBOTH.EQ.1)) THEN
         
             ISOSPI=1
             
             ERRMAX_PRINTG=ERRMAX_PRONUC(INUCLI)
C             EABSAV_PRINTG=EABSAV_PRONUC(INUCLI)
             RMSVAL_PRINTG=RMSVAL_PRONUC(INUCLI)
             
             RMSEXP_PRINTG=RMSEXP_PROTON(INUCLI)
             RMSTHE_PRINTG=RMSTHE_PRONUC(INUCLI)
             LEVTHE_PRINTG=LEVTHE_PRONUC(INUCLI)
             
             DO ITHEOR=1,LEVTHE_PRINTG
                ENETHE_PRINTG(ITHEOR)=ENETHE_PRONUC(ITHEOR,INUCLI)
                LABTHE_PRINTG(ITHEOR)=LABTHE_PRONUC(ITHEOR,INUCLI)
             END DO
             
             CALL WRITIN_ENELEV(ISOSPI,IRANDO_PRINTG,LDRAND,
     *                                        IEVALU_PRINTG,
     *                          PARPOT_PRINTG,CHISQU_MINIML,
     *                          CHIGRD_PRINTG,ERRMAX_PRINTG,
     *                          RMSVAL_PRINTG,RMSGLO_PROAUX,
     *                          RMSEXP_PRINTG,RMSTHE_PRINTG,
     *                          LEVTHE_PRINTG,ENETHE_PRINTG,
     *                          LABTHE_PRINTG,INUCLI)
C_______________________________________________________________________
C     
             ISOSPI=0
             
             ERRMAX_PRINTG=ERRMAX_NEUNUC(INUCLI)
C             EABSAV_PRINTG=EABSAV_NEUNUC(INUCLI)
             RMSVAL_PRINTG=RMSVAL_NEUNUC(INUCLI)
             
             RMSEXP_PRINTG=RMSEXP_NEUTRS(INUCLI)
             RMSTHE_PRINTG=RMSTHE_NEUNUC(INUCLI)
             LEVTHE_PRINTG=LEVTHE_NEUNUC(INUCLI)
             
             DO ITHEOR=1,LEVTHE_PRINTG
                ENETHE_PRINTG(ITHEOR)=ENETHE_NEUNUC(ITHEOR,INUCLI)
                LABTHE_PRINTG(ITHEOR)=LABTHE_NEUNUC(ITHEOR,INUCLI)
             END DO
             
             CALL WRITIN_ENELEV(ISOSPI,IRANDO_PRINTG,LDRAND,
     *                                        IEVALU_PRINTG,
     *                          PARPOT_PRINTG,CHISQU_MINIML,
     *                          CHIGRD_PRINTG,ERRMAX_PRINTG,
     *                          RMSVAL_PRINTG,RMSGLO_NEUAUX,
     *                          RMSEXP_PRINTG,RMSTHE_PRINTG,
     *                          LEVTHE_PRINTG,ENETHE_PRINTG,
     *                          LABTHE_PRINTG,INUCLI)
         END IF
C
C=======================================================================
C     
C        Printing the results in the standard output
C
         IF (LOGWRI.GT.0) THEN 
             WRITE(NOUTPT,'()')
       
             WRITE(NOUTPT,'(80(''#''),/,''#'',78X,''#'')')
             WRITE(NOUTPT,'(''#'',20X,''IZ_FIX='',I3,20X,''IN_FIX='',I3,
     *                                   18X,''#'',/,''#'',78X,''#'')')
     *                         IZ_FIX,IN_FIX
             WRITE(NOUTPT,'(80(''#''),/,''#'',78X,''#'')')
             WRITE(NOUTPT,'(''#'',2X,''NEUTRONS: Theoretical vs. '',
     *                      ''Experimental Energies'',29x,''#'')')
             WRITE(NOUTPT,'(''#'',78X,''#'',/,80(''#''))')
             WRITE(NOUTPT,'(''#'',78X,''#'')')
       
             WRITE(NOUTPT,'(''#'',3X,''No)'',7x,''Labels'',7x,
     *          ''EneThe'',14X,''No)'',7x,''Labels'',7x,''EneExp'',3X,
     *                                                    ''#'')')
       
             DO ITHEOR=1,LEVTHE_NEUNUC(INUCLI)
           
                DO IEXPER=1,LEVEXP_NEUTRS(INUCLI)
                   IF (LABEXP_NEUTRS(INUCLI,IEXPER).EQ.
     *                 LABTHE_NEUNUC(ITHEOR,INUCLI))THEN
                       WRITE(NOUTPT,'(''#'',i5,'')'',7x,a6,f20.13,
     *                          4X,i5,'')'',7x,a6,f15.4,1X,''#'')')
     *                                  ITHEOR,
     *                                  LABTHE_NEUNUC(ITHEOR,INUCLI),
     *                                  ENETHE_NEUNUC(ITHEOR,INUCLI),
     *                                  IEXPER,
     *                                  LABEXP_NEUTRS(INUCLI,IEXPER),
     *                                  EXPEXP_NEUTRS(INUCLI,IEXPER)
                       GO TO 3
                   END IF
                END DO
           
                WRITE(NOUTPT,'(''#'',i5,'')'',7x,a6,f20.13,39X,''#'')')    
     *                         ITHEOR,LABTHE_NEUNUC(ITHEOR,INUCLI),
     *                                ENETHE_NEUNUC(ITHEOR,INUCLI)
  3             CONTINUE
             END DO
       
             WRITE(NOUTPT,'(''#'',78X,''#'',/,80(''#''),/,
     *                                 ''#'',78X,''#'')')
       
             WRITE(NOUTPT,'(''#'',2X,''PROTONS: Theoretical vs. '',
     *                           ''Experimental Energies'',30x,''#'')')
             WRITE(NOUTPT,'(''#'',78X,''#'',/,80(''#''))')
             WRITE(NOUTPT,'(''#'',78X,''#'')')
       
             WRITE(NOUTPT,'(''#'',3X,''No)'',7x,''Labels'',7x,
     *           ''EneThe'',14X,''No)'',7x,''Labels'',7x,''EneExp'',3X,
     *                                                    ''#'')')
       
              DO ITHEOR=1,LEVTHE_PROTON
           
                 DO IEXPER=1,LEVEXP_PROTON(INUCLI)
                    IF (LABEXP_PROTON(INUCLI,IEXPER).EQ.
     *                  LABTHE_PRONUC(ITHEOR,INUCLI)) THEN
                        WRITE(NOUTPT,'(''#'',i5,'')'',7x,a6,f20.13,
     *                          4X,i5,'')'',7x,a6,f15.4,1X,''#'')')
     *                                  ITHEOR,
     *                                  LABTHE_PRONUC(ITHEOR,INUCLI),
     *                                  ENETHE_PRONUC(ITHEOR,INUCLI),
     *                                  IEXPER,
     *                                  LABEXP_PROTON(INUCLI,IEXPER),
     *                                  EXPEXP_PROTON(INUCLI,IEXPER)
                        GO TO 4
                    END IF
                 END DO
           
                 WRITE(NOUTPT,'(''#'',i5,'')'',7x,a6,f20.13,39X,''#'')')
     *                           ITHEOR,LABTHE_PRONUC(ITHEOR,INUCLI),
     *                                  ENETHE_PRONUC(ITHEOR,INUCLI)
  4             CONTINUE
             END DO
       
             WRITE(NOUTPT,'(''#'',78X,''#'',/,80(''#''))')
C
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
         DO IPOINT=1,ND_MAX(I_CURV)
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
          WRITE(LOGFIL,'(12X,''Exiting FINISH_LMMINI'')')
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
      SUBROUTINE MESHIN_MAPING(LDRAND)
C
      INCLUDE  'MATDIM/NDNUCL.f'
      INCLUDE  'MATDIM/NDPARS.f'
      INCLUDE  'MATDIM/NDTITL.f'
      INCLUDE  'MATDIM/NDIM_M.f'
      INCLUDE  'MATDIM/NDMESH.f'
      INCLUDE  'MATDIM/NDSPEC.f'
      INCLUDE  'MATDIM/NDLEXP.f'
      INCLUDE  'MATDIM/NDLAST.f'
C
      PARAMETER
     *         (NDMES2=NDMESH*NDMESH)
C
      CHARACTER
     *          INPSYM*6,TYPCHI*6,TAKCHI*3,NUCSYM*6
      CHARACTER
     *          TITPAR*13,TITMSH*8,TITLES*12,TITAUX*7,TAKPAR*3,
     *          TITLES_LATEXS*050,NUCNAM_LATEXS*010
C
      DIMENSION
     *          ARGMNT(1:NDPARS)
      DIMENSION
     *          INDEXI(1:NDIM_M),MCOUNT(1:NDIM_M),
     *          STPMSH(1:NDIM_M),MAXMSH(1:NDPARS)
      DIMENSION
     *          PARMSH(1:NDIM_M,1:NDMESH),
     *          AXSMAP(1:NDMESH,1:NDMESH)
      DIMENSION
     *          TITPAR(1:NDPARS),TITMSH(1:NDPARS),
     *          TAKCHI(1:NDNUCL),TITAUX(1:NDIM_M),
     *                           TAKPAR(1:NDIM_M)
      DIMENSION
     *          CHITOT_OFMESH(1:NDMESH,1:NDMESH),
     *          CHIPRO_OFMESH(1:NDMESH,1:NDMESH),
     *          CHINEU_OFMESH(1:NDMESH,1:NDMESH)
      DIMENSION
     *          PARPOT_OFMESH(1:NDMESH,1:NDMESH,1:NDPARS),
     *          PARPOT_PRINTI(1:NDPARS,1:NDMES2)
      DIMENSION
     *          RMSMAP_PRONUC(1:NDMESH,1:NDMESH,1:NDNUCL),
     *          RMSMAP_NEUNUC(1:NDMESH,1:NDMESH,1:NDNUCL)
      DIMENSION
     *          RMSMAP_PROGLO(1:NDMESH,1:NDMESH,1:NDNUCL),
     *          RMSMAP_NEUGLO(1:NDMESH,1:NDMESH,1:NDNUCL)
C
      COMMON
     *       /POTPOT/ PARPOT(1:NDPARS)
      COMMON
     *       /ARGINI/ XINITS(1:NDPARS)
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
     *       /NUCLEU/ IZ_FIX,IN_FIX,ISOSPI  
      COMMON
     *       /DENSTR/ IFDENS
     *       /TENSOR/ IFTENS  
      COMMON
     *       /TENOPT/ ICENTT,ISORBT,ITENSR
     *       /TESFIT/ IFTEST,IFFITS,IFMESH,ISPE_P,IFMOCA
      COMMON
     *       /FEVALS/ IFUNCT_EVALUS
     *       /P_N_LM/ IFPROT,IFNEUT,IFBOTH
     *       /MINITE/ IDEFCN,ITECHI
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
      IDEFCN=1
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
C=======================================================================     
C         Which parameters stay 'fixed' --> p_i and p_j --> PARMSH(N,M)
C=======================================================================
C          
      ICOUNT_OFMESH=0
C
      DO IPARAM=1,NDPARS
C
         IF (I_MESH(IPARAM).EQ.1) THEN
C
             IF (ICOUNT_OFPARS.EQ.0) THEN
                 IF (IPARAM.LE.20 .AND. IFDENS.EQ.0) IFPROT=1
                 IF (IPARAM.GT.20 .AND. IFDENS.EQ.0) IFNEUT=1
                 IF (IPARAM.GE.51 .AND. IFDENS.EQ.0) THEN
                     IFPROT=0
                     IFNEUT=0
                     IFBOTH=1
                 END IF
             END IF
C
             ICOUNT_OFMESH=ICOUNT_OFMESH+1
C
             STPMSH(IPARAM)=(XMAX_I(IPARAM)-XMIN_I(IPARAM))
     *                     /(MAXPAR(IPARAM)-1)
C
             MAXMSH(ICOUNT_OFMESH)=MAXPAR(IPARAM)
C
             INDEXI(ICOUNT_OFMESH)=IPARAM
C
             WRITE(TITMSH(ICOUNT_OFMESH),'(A8)')TITLES(IPARAM)(3:10)
C
             DO K=1,MAXPAR(IPARAM)
                AXSMAP(ICOUNT_OFMESH,K)=XMIN_I(IPARAM)
     *                                 +(K-1)*STPMSH(IPARAM)
                IF (MAXPAR(IPARAM).EQ.1) THEN
                    AXSMAP(ICOUNT_OFMESH,K)=XMIN_I(IPARAM)
                END IF
             END DO
C
         END IF
C
      END DO
C
      IF (ICOUNT_OFMESH.NE.2) THEN
C
          WRITE(LSCREN,'(''ALARM from MESHIN_MAPING: '',
     *                   ''ICOUNT_OFMESH='',I2,''when it has to be '',
     *                   ''equal to 2!! - Change the input file '')')
     *                    ICOUNT_OFMESH
          STOP 'Stop - Change the input file in MESH option'
C
      END IF
C
      ICMESH=ICOUNT_OFMESH
C
C=======================================================================     
C         
      NTOTAL_OFMESH=MAXMSH(1)*MAXMSH(2)
      NCOUNT_OFMESH=0
C
      IF (LOGWRI.GT.0) THEN
          WRITE(LOGFIL,'(''Entering OPENIN_MESHPM from Main'')')
      END IF
CID      CALL OPENIN_MESHPM(LDRAND,TITMSH,ICOUNT_OFPARS,NTOTAL_OFMESH,
CID     *                                                      MAXMSH)
C
C=======================================================================     
C
      RMSPRO_MAXIMA=-1E+10
      RMSPRO_MINIMI=+1E+10
C
      RMSNEU_MAXIMA=-1E+10
      RMSNEU_MINIMI=+1E+10
C
C=======================================================================     
C         
      IPARA1=INDEXI(1)
      IPARA2=INDEXI(2)
C
      DO INDEX1=1,MAXMSH(1)
C
         MCOUNT(1)=INDEX1
C
         XINITS(IPARA1)=XMIN_I(IPARA1)+(INDEX1-1)*STPMSH(IPARA1)
C
         IF (MAXMSH(1).EQ.1) XINITS(IPARA1)=XMIN_I(IPARA1)
C
         DO INDEX2=1,MAXMSH(2)
C
            MCOUNT(2)=INDEX2
C
            XINITS(IPARA2)=XMIN_I(IPARA2)+(INDEX2-1)*STPMSH(IPARA2)
C
            IF (MAXMSH(2).EQ.1) XINITS(IPARA2)=XMIN_I(IPARA2)
C
            NCOUNT_OFMESH=NCOUNT_OFMESH+1
C
C=======================================================================
C 
C           Calling the minimisation subroutine. 
C           We are inside the INDEX1 and INDEX2 do-loops.
C
            IF (ICOUNT_OFPARS.GT.0) THEN  ! Tabulating and Minimisation
C
                CALL LMMESH(NEWSED,ARGMNT,I_SEED,LDRAND,LEVNUM,
     *                      NTOTAL_OFMESH,NCOUNT_OFMESH,TITMSH,
     *                      ICMESH,ICOUNT_OFPARS,CHITOT_MINAUX,
     *                             CHISQU_PROAUX,CHISQU_NEUAUX,
     *                                    PARPOT_PRINTI,NDLAST)
C
                DO KPARAM=1,NDPARS
                   PARPOT(KPARAM)=PARPOT_PRINTI(KPARAM,NCOUNT_OFMESH)
                END DO
C_______________________________________________________________________
C
C              Once we have the vector of parameters, we enter to 
C              WS_RUN to calculate the RMSEnergy values that will
C                                                       be Mapped
C
               IDEG_P=0
               IDEG_N=0
               RMSGLO_TAKPRO=0.
               RMSGLO_TAKNEU=0.
C
               DO INUCLI=1,LDNUCL
                  IF (ITAKNU(INUCLI).EQ.1) THEN 
                      IZ_FIX=NUMB_Z(INUCLI)
                      IN_FIX=NUMB_N(INUCLI)
C
                      IF (LOGWRI.GT.4) THEN
                          WRITE(LOGFIL,'(12X,''Entering WS_RUN from '',
     *                                       ''MESHIN_MAPING '',
     *                                       ''INUCLI='',I2)') INUCLI
                      END IF
C
                      I_MODE=0
                      I_FLAG=1
C
                      CALL WS_RUN(INUCLI,CHISQU_CHOICE,IABORT,I_MODE,
     *                            I_FLAG,CHISQU_PROTON,CHISQU_NEUTRS)
C 
                      RMSMAP_PRONUC(INDEX1,INDEX2,INUCLI)
     *               =
     *                SQRT(CHIWEI_PROTON)
C
                      RMSMAP_NEUNUC(INDEX1,INDEX2,INUCLI)
     *               =
     *                SQRT(CHIWEI_NEUTRS)
C
C_______________________________________________________________________
C
C          Calculating the TOTAL degenerancy of the active nuclei
C
                      DO IEXPER=1,LEVEXP_PROTON(INUCLI)
                         IDEG_P=IDEG_P+IDEGEX_PROTON(INUCLI,IEXPER)
                      END DO
C
                      DO IEXPER=1,LEVEXP_NEUTRS(INUCLI)
                         IDEG_N=IDEG_N+IDEGEX_NEUTRS(INUCLI,IEXPER)
                      END DO
C
                      RMSGLO_TAKPRO=RMSGLO_TAKPRO+CHIDEG_PROTON
                      RMSGLO_TAKNEU=RMSGLO_TAKNEU+CHIDEG_NEUTRS
C
C=======================================================================
C
C                     Searching the RMSMAP maximum/minimum values
C
                      IF (RMSPRO_MAXIMA.LT.
     *                    RMSMAP_PRONUC(INDEX1,INDEX2,INUCLI)) THEN
                          RMSPRO_MAXIMA
     *                   =
     *                    RMSMAP_PRONUC(INDEX1,INDEX2,INUCLI)
                      END IF
C
                      IF (RMSPRO_MINIMI.GT.
     *                    RMSMAP_PRONUC(INDEX1,INDEX2,INUCLI)) THEN
                          RMSPRO_MINIMI
     *                   =
     *                    RMSMAP_PRONUC(INDEX1,INDEX2,INUCLI)
                      END IF
C
                      IF (RMSNEU_MAXIMA.LT.
     *                    RMSMAP_NEUNUC(INDEX1,INDEX2,INUCLI)) THEN
                          RMSNEU_MAXIMA
     *                   =
     *                    RMSMAP_NEUNUC(INDEX1,INDEX2,INUCLI)
                      END IF
C
                      IF (RMSNEU_MINIMI.GT.
     *                    RMSMAP_NEUNUC(INDEX1,INDEX2,INUCLI)) THEN
                          RMSNEU_MINIMI
     *                   =
     *                    RMSMAP_NEUNUC(INDEX1,INDEX2,INUCLI)
                      END IF
C
                  END IF
C
               END DO 
C
               RMSMAP_PROGLO(INDEX1,INDEX2,LDNUCL+1)
     *        =
     *         SQRT(RMSGLO_TAKPRO/IDEG_P)
C
               RMSMAP_NEUGLO(INDEX1,INDEX2,LDNUCL+1)
     *        =
     *         SQRT(RMSGLO_TAKNEU/IDEG_N)
C            
               DO IPARAM=1,NDPARS
                  PARPOT_OFMESH(INDEX1,INDEX2,IPARAM)
     *           =
     *            PARPOT_PRINTI(IPARAM,NCOUNT_OFMESH)
C
               END DO
C
            END IF
C
C=======================================================================
C
            IF (ICOUNT_OFPARS.EQ.0) THEN ! Only Tabulating,
C                                          NO Minimisation

                CALL TBMESH(IPARA1,IPARA2,INDEX1,INDEX2,LEVNUM,
     *                             CHITOT_MINAUX,CHISQU_PROAUX,
     *                             CHISQU_NEUAUX,PARPOT_PRINTI,
     *                             RMSMAP_PRONUC,RMSMAP_NEUNUC,
     *                             RMSMAP_PROGLO,RMSMAP_NEUGLO,
     *                             RMSPRO_MAXIMA,RMSPRO_MINIMI,
     *                             RMSNEU_MAXIMA,RMSNEU_MINIMI)
C
                DO IPARAM=1,NDPARS
                   PARPOT_OFMESH(INDEX1,INDEX2,IPARAM)
     *            =
     *             PARPOT(IPARAM)
                END DO
C                
            END IF
C
C=======================================================================
C            
         END DO !INDEX2
C        
      END DO !INDEX1
C
C=======================================================================
C         Writing MAP results in a file  
C=======================================================================
C
C     First, each nucleus individually
C
      DO INUCLI=1,LDNUCL
         IF(ITAKNU(INUCLI).EQ.1) THEN 
C
            WRITE(LOGFIL,'(''Entering WRITIN_CHIMAP with INUCLI= '',
     *                        I3,''from MESHIN_MAPING'')') INUCLI
C
            IZ_FIX=NUMB_Z(INUCLI)
            IN_FIX=NUMB_N(INUCLI)
C
            ISOSPI=1
            CALL WRITIN_CHIMAP(NDMESH,LDRAND,ICOUNT_OFPARS,
     *                         TITPAR,TITMSH,RMSPRO_MINIMI,
     *                         RMSPRO_MAXIMA,RMSNEU_MINIMI,
     *                         RMSNEU_MAXIMA,MAXMSH,AXSMAP,
     *                         RMSMAP_PRONUC,RMSMAP_NEUNUC,
     *                         PARPOT_OFMESH,INDEXI,INUCLI)
C
            ISOSPI=0
            CALL WRITIN_CHIMAP(NDMESH,LDRAND,ICOUNT_OFPARS,
     *                         TITPAR,TITMSH,RMSPRO_MINIMI,
     *                         RMSPRO_MAXIMA,RMSNEU_MINIMI,
     *                         RMSNEU_MAXIMA,MAXMSH,AXSMAP,
     *                         RMSMAP_PRONUC,RMSMAP_NEUNUC,
     *                         PARPOT_OFMESH,INDEXI,INUCLI)
         END IF
      END DO
C
C     Then, the global results
C
      INUCLI=LDNUCL+1
      IZ_FIX=000
      IN_FIX=000
C
      ISOSPI=11
      CALL WRITIN_CHIMAP(NDMESH,LDRAND,ICOUNT_OFPARS,
     *                   TITPAR,TITMSH,RMSPRO_MINIMI,
     *                   RMSPRO_MAXIMA,RMSNEU_MINIMI,
     *                   RMSNEU_MAXIMA,MAXMSH,AXSMAP,
     *                   RMSMAP_PROGLO,RMSMAP_NEUGLO,
     *                   PARPOT_OFMESH,INDEXI,INUCLI)
C
      ISOSPI=10
      CALL WRITIN_CHIMAP(NDMESH,LDRAND,ICOUNT_OFPARS,
     *                   TITPAR,TITMSH,RMSPRO_MINIMI,
     *                   RMSPRO_MAXIMA,RMSNEU_MINIMI,
     *                   RMSNEU_MAXIMA,MAXMSH,AXSMAP,
     *                   RMSMAP_PROGLO,RMSMAP_NEUGLO,
     *                   PARPOT_OFMESH,INDEXI,INUCLI)
C
C=======================================================================
C
      RETURN
      END     
C
C=======================================================================
C=======================================================================
C
      SUBROUTINE LMMESH(NEWSED,ARGPAR,I_SEED,LDRAND,LEVNUM,
     *                  NTOTAL_OFMESH,NCOUNT_OFMESH,TITMSH,
     *                  ICMESH,ICOUNT_OFPARS,CHITOT_MINAUX,
     *                         CHISQU_PROAUX,CHISQU_NEUAUX,
     *                                PARPOT_PRINTI,NDLAST)
C
      INCLUDE   'MATDIM/NDSPEC.f'
      INCLUDE   'MATDIM/NDPARS.f'
      INCLUDE   'MATDIM/NDNUCL.f'
      INCLUDE   'MATDIM/NDLEXP.f'
      INCLUDE   'MATDIM/NDIM_M.f'
      INCLUDE   'MATDIM/NDMESH.f'
      INCLUDE   'MATDIM/MAXFEV.f'
      INCLUDE   'MATDIM/NDTITL.f'
C      
      PARAMETER 
     *         (NDMES2=NDMESH*NDMESH,INMODE=1,NPRINT=-1)
C
      EXTERNAL 
     *          FUNMIN
C
      CHARACTER
     *          DIRNAM*256,GENAME*256,STRING*126
      CHARACTER
     *          LABTHE_PROTON*6,LABTHE_NEUTRS*6
      CHARACTER
     *          LABEXP_PROTON*6,LABEXP_NEUTRS*6
      CHARACTER
     *          ACTION*1,CALLED*6,TYPCHI*6
      CHARACTER
     *          INPSYM*6,TEXTRA*256,NUCSYM*6
      CHARACTER
     *          TITMSH*8,FILNAM_ENDING*1,FILNAM*256,TITLES*12
      CHARACTER
     *          FILNA2*256,FITNUC*2,FILNAM_NUCFIT*30,
     *          FITNUC_AUXILI*30,FILNA3*256,TITLES_NOTFIT*256,
     *          FILNAM_NOTTIT*256,TEXLAM*20,VERSIO*3,
     *          TITLES_LATEXS*050,NUCNAM_LATEXS*010
C    
      DIMENSION 
     *          ARGPAR(1:NDPARS),FUNVEC(1:NDIM_M),
     *          SCALFC(1:NDPARS),QTFARR(1:NDPARS)
      DIMENSION
     *          WORKA1(1:NDPARS),WORKA2(1:NDPARS),
     *          WORKA3(1:NDPARS),WORKA4(1:NDIM_M)
      DIMENSION 
     *          FUNJAC(1:NDIM_M,1:NDPARS)
      DIMENSION 
     *          IPIVOT(1:NDPARS)
      DIMENSION
     *          TITMSH(1:NDPARS),
     *          IAUXIL(1:NDPARS)
      DIMENSION
     *          PARPOT_PRINTI(1:NDPARS,1:NDMES2)
      DIMENSION
     *          FITNUC(1:NDNUCL),
     *          FITNUC_AUXILI(1:NDNUCL),
     *          TITLES_NOTFIT(1:NDPARS)
C
      COMMON
     *       /DIMSPH/ LDSING_PROTON,LEVTHE_PROTON,
     *                LEVEXP_PROTON(1:NDNUCL),
     *                LDSING_NEUTRS,LEVTHE_NEUTRS,
     *                LEVEXP_NEUTRS(1:NDNUCL)
      COMMON
     *       /CHICHO/ VALCHI
      COMMON
     *       /CHI_PN/ CHISQU_PROTON,CHISQU_NEUTRS
      COMMON
     *       /CHIGRD/ CHISQU_GRDNRM,
     *                CHISQU_GRADIE(1:NDPARS)
      COMMON 
     *       /CHITAK/ IF_SPE,IF_RAD,IF_GAP,IF_FER,IF_DEN,IF_RHO,IF_INV,
     *                LDCHI2,TYPCHI(1:NDNUCL),ITKCHI(1:NDNUCL)
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
     *       /MESHIN/ I_MESH(1:NDPARS),
     *                XMIN_I(1:NDPARS),
     *                XMAX_I(1:NDPARS),
     *                MAXPAR(1:NDPARS)
      COMMON
     *       /MSHPRI/ INFOER_PRINTI,
     *                IDEFCN_PRINTI,
     *                CHISQU_PRINTI(1:NDMES2),
     *                CHIGRD_NRMPRI(1:NDMES2),
     *                CHIGRD_COMPRI(1:NDPARS,1:NDMES2)     
      COMMON
     *       /POTPOT/ PARPOT(1:NDPARS)
      COMMON
     *       /MINITE/ IDEFCN,ITECHI
      COMMON
     *       /P_N_LM/ IFPROT,IFNEUT,IFBOTH
      COMMON
     *       /DENSTR/ IFDENS
      COMMON
     *       /TENSOR/ IFTENS
      COMMON
     *       /TENOPT/ ICENTT,ISORBT,ITENSR
      COMMON
     *       /IFPAIR/ IF_PAI,IFGPAI,IFDPAI
     *       /DEPPAI/ IFDEEP,IFPRON
      COMMON
     *       /IFKAPP/ IFK_VC,IFK_VS,IFK_RC,IFK_RS,IFK_AC,IFK_AS
      COMMON
     *       /IFLAMB/ IFPAR1,IFPAR2,IFPAR3,IFPAR4,
     *                IFPAR5,IFPAR6,IFPAR7,IFPAR8,
     *                IFPA09,IFPA10,IFPA11,IFPA12
      COMMON
     *       /FEVALS/ IFUNCT_EVALUS
      COMMON
     *       /STOPAR/ EPSLAS,LDLAST
      COMMON
     *       /STOPIN/ HIAUXI(1:MAXFEV)
      COMMON
     *       /PARAMT_TITLES/ TITLES(1:NDTITL),
     *                       TITLES_LATEXS(1:NDTITL),
     *                       NUCNAM_LATEXS(1:NDNUCL)
      COMMON
     *       /PRFILE/ ISCREN,LOGFIL,LOGCHP,LOGCHN,LOGHIP,LOGHIN,
     *                LOGPRO,LOGNEU,LOGWRI,LSCREN,IRESUL,ICONVE,
     *                NOUTPT,LOGAUX,LOGENE,LOGRAD,LOGGRD,LOGMSH,
     *                                                   LOGBIS  
      COMMON
     *       /TOLERS/ TOLERF,TOLERX,TOLERG,FACTOR   
      COMMON  
     *       /VERSIN/ VERSIO   
C
      DATA 
     *     NpUNIT / 80 / 
      DATA 
     *     NnUNIT / 81 /
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
C     Levenberg-Marquardt minimisation routine ??
C=======================================================================
C
C     ILEVEL - Number of functions
C     NDPARS - Number of parameters (not greater than ILEVEL)
C
C            chi=sum_m(FUNJAC_m)^2 
C
C     Testing case setting ILEVEL=1
C     NDIM_M leading dimension (ILEVEL)
C
C=======================================================================
C 
      CALL CPUTIM('LMMESH',1)
C
C=======================================================================
C      
      IF (LOGWRI.GT.0) THEN
          WRITE(LOGFIL,'(/,''Entering LMMESH'')')
      END IF
C
C=======================================================================
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
C  @@@ IRENE? CHECK THE PRIVILEGED ROLE OF 20
             IF (IPARAM.LE.20 .AND. IFDENS.EQ.0) THEN
                 IFPROT=1
             END IF
C
             IF (IPARAM.GT.20 .AND. IFDENS.EQ.0) THEN
                 IFNEUT=1
             END IF
C
             IF (IPARAM.GE.51 .AND. IFDENS.EQ.0) THEN
                 IFPROT=0
                 IFNEUT=0
                 IFBOTH=1
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
          WRITE(0,'(/,''ALERT in LMMESH!! Density is NOT activated'',1X,
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
              WRITE(LOGFIL,'(9X,''Minimisation for  protons'')')
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
      INDEXP=0
C      
      IF (IFDENS.EQ.0.AND.IFPROT.EQ.1) THEN
          DO IPARAM=1,6
              INDEXP=INDEXP+1
              IAUXIL(INDEXP)=IPARAM
          END DO
      END IF
C      
      IF (IFDENS.EQ.0.AND.IFNEUT.EQ.1) THEN
          DO IPARAM=21,26
              INDEXP=INDEXP+1
              IAUXIL(INDEXP)=IPARAM
          END DO
      END IF
C
      IF (IFDENS.EQ.0.AND.IFBOTH.EQ.1) THEN
C
          IF (IFK_VC.EQ.1) THEN
              INDEXP=INDEXP+1
              IAUXIL(INDEXP)=51   ! V_o central
              INDEXP=INDEXP+1
              IAUXIL(INDEXP)=52   ! kappa V_o central
          ELSE
              INDEXP=INDEXP+1 
              IAUXIL(INDEXP)=1    ! V_o central protons
              INDEXP=INDEXP+1
              IAUXIL(INDEXP)=21   ! V_o central neutrons
          END IF
C
          IF (IFK_RC.EQ.1) THEN
              INDEXP=INDEXP+1
              IAUXIL(INDEXP)=53   ! r_o central
              INDEXP=INDEXP+1
              IAUXIL(INDEXP)=54   ! kappa r_o central
          ELSE
              INDEXP=INDEXP+1
              IAUXIL(INDEXP)=2    ! r_o central protons
              INDEXP=INDEXP+1
              IAUXIL(INDEXP)=22   ! r_o central neutron
          END IF
C
          IF (IFK_AC.EQ.1) THEN
              INDEXP=INDEXP+1
              IAUXIL(INDEXP)=55   ! a_o central
              INDEXP=INDEXP+1
              IAUXIL(INDEXP)=56   ! kappa a_o central
          ELSE
              INDEXP=INDEXP+1
              IAUXIL(INDEXP)=3    ! a_o central protons
              INDEXP=INDEXP+1
              IAUXIL(INDEXP)=23   ! a_o central neutrons
          END IF
C
          IF (IFK_VS.EQ.1) THEN
              INDEXP=INDEXP+1
              IAUXIL(INDEXP)=57   ! V_o pure WS-SO
              INDEXP=INDEXP+1
              IAUXIL(INDEXP)=58   ! kappa V_o pure WS-SO
          ELSE
              INDEXP=INDEXP+1
              IAUXIL(INDEXP)=4    ! V_o pure WS-SO protons
              INDEXP=INDEXP+1
              IAUXIL(INDEXP)=24   ! V_o pure WS-SO neutrons
          END IF
C
          IF (IFK_RS.EQ.1) THEN
              INDEXP=INDEXP+1
              IAUXIL(INDEXP)=59   ! r_o pure WS-SO
              INDEXP=INDEXP+1
              IAUXIL(INDEXP)=60   ! kappa r_o pure WS-SO
          ELSE
              INDEXP=INDEXP+1
              IAUXIL(INDEXP)=5    ! r_o pure WS-SO protons
              INDEXP=INDEXP+1
              IAUXIL(INDEXP)=25   ! r_o pure WS-SO neutrons
          END IF
C
          IF (IFK_AS.EQ.1) THEN
              INDEXP=INDEXP+1
              IAUXIL(INDEXP)=61   ! a_o pure WS-SO
              INDEXP=INDEXP+1
              IAUXIL(INDEXP)=62   ! kappa a_o pure WS-SO
          ELSE
              INDEXP=INDEXP+1
              IAUXIL(INDEXP)=6    ! a_o pure WS-SO protons
              INDEXP=INDEXP+1
              IAUXIL(INDEXP)=26   ! a_o pure WS-SO neutrons
          END IF
C
      END IF
C_______________________________________________________________________
C      
      IF (IFDENS.EQ.1) THEN
C
          IF (IFK_VC.EQ.1) THEN
              INDEXP=INDEXP+1
              IAUXIL(INDEXP)=51   ! V_o central
              INDEXP=INDEXP+1
              IAUXIL(INDEXP)=52   ! kappa V_o central
          ELSE
              INDEXP=INDEXP+1 
              IAUXIL(INDEXP)=1    ! V_o central protons
              INDEXP=INDEXP+1
              IAUXIL(INDEXP)=21   ! V_o central neutrons
          END IF
C
          IF (IFK_RC.EQ.1) THEN
              INDEXP=INDEXP+1
              IAUXIL(INDEXP)=53   ! r_o central
              INDEXP=INDEXP+1
              IAUXIL(INDEXP)=54   ! kappa r_o central
          ELSE
              INDEXP=INDEXP+1
              IAUXIL(INDEXP)=2    ! r_o central protons
              INDEXP=INDEXP+1
              IAUXIL(INDEXP)=22   ! r_o central neutron
          END IF
C
          IF (IFK_AC.EQ.1) THEN
              INDEXP=INDEXP+1
              IAUXIL(INDEXP)=55   ! a_o central
              INDEXP=INDEXP+1
              IAUXIL(INDEXP)=56   ! kappa a_o central
          ELSE
              INDEXP=INDEXP+1
              IAUXIL(INDEXP)=3    ! a_o central protons
              INDEXP=INDEXP+1
              IAUXIL(INDEXP)=23   ! a_o central neutrons
          END IF
C
          DO IPARAM=39,42         !lambdas Density-SO
             INDEXP=INDEXP+1
             IAUXIL(INDEXP)=IPARAM
          END DO
C
          IF (IFTENS.EQ.1.AND.ISORBT.EQ.1) THEN
              DO IPARAM=43,46     !lambdas Tensor-SO
                 INDEXP=INDEXP+1
                 IAUXIL(INDEXP)=IPARAM
              END DO
          END IF
C
          IF (IFTENS.EQ.1.AND.ICENTT.EQ.1) THEN
              DO IPARAM=47,50     !lambdas Tensor-central
                 INDEXP=INDEXP+1
                 IAUXIL(INDEXP)=IPARAM
              END DO
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
          FILNAM_ENDING='B'
C              
          WRITE(TEXLAM,'(''_LAMBDA'',SP,I2,I2,I2,I2)')
     *                               IFPAR1,IFPAR2,IFPAR3,IFPAR4
C          
      END IF
C_______________________________________________________________________
C 
      KACTIV=0
      DO JNUCLI=1,LDNUCL
         IF (ITAKNU(JNUCLI).EQ.1) THEN
             KACTIV=KACTIV+1
             FITNUC_AUXILI(KACTIV)=FITNUC(JNUCLI)
         END IF
      END DO
C      
      I1=3*KACTIV-1
C      
      WRITE(FILNAM_NUCFIT,'(<NUCACT>(A2,''-''))')
     *                      (FITNUC_AUXILI(I),I=1,NUCACT)
C
C=======================================================================
C      
      MACTIV=0
      ICENTR_NOTFIT=0
      ISORBI_NOTFIT=0
      IDENSI_NOTFIT=0
C
      DO I=1,INDEXP
         IPARAM=IAUXIL(I)
         IF (IFTAKE(IPARAM).EQ.0 .AND. I_MESH(IPARAM).EQ.0) THEN
C
             MACTIV=MACTIV+1
             TITLES_NOTFIT(MACTIV)=TITLES(IPARAM)(3:10)
C
             IF (TITLES_NOTFIT(MACTIV)(3:6).EQ.'CENT') THEN
                 ICENTR_NOTFIT=ICENTR_NOTFIT+1
             END IF
             IF (TITLES_NOTFIT(MACTIV)(3:6).EQ.'SORB') THEN
                 ISORBI_NOTFIT=ISORBI_NOTFIT+1
             END IF
             IF (TITLES_NOTFIT(MACTIV)(1:5).EQ.'XLAMB') THEN
                 IDENSI_NOTFIT=IDENSI_NOTFIT+1
             END IF
C
         END IF
      END DO
C
      I2=9*MACTIV-1
C
      WRITE(FILNAM_NOTTIT,'(<MACTIV>(A8,''-''))')
     *                      (TITLES_NOTFIT(I),I=1,MACTIV)
C
      IF (IFDENS.EQ.1 .AND. ICENTR_NOTFIT.GT.0) THEN
          WRITE(FILNAM_NOTTIT,'(I1,''CentrPot-const'')')ICENTR_NOTFIT
          I2=15
          IF (IDENSI_NOTFIT.GT.0) THEN
              WRITE(FILNAM_NOTTIT,'(I1,''CentrPot-const_'',
     *                              I1,''SODenPot-const'')')
     *                                      ICENTR_NOTFIT,IDENSI_NOTFIT
              I2=31
          END IF
      END IF
C
C=======================================================================
C
      WRITE(0,'(<ICMESH>(''_'',A8))')(TITMSH(I),I=1,ICMESH)
C
      IF (IFDENS.EQ.0) THEN
C      
          IF (MACTIV.EQ.0) 
     *        WRITE(FILNAM,'(''LogFMesh/IFDENS-0/'',A,''_'',A1,
     *                   <ICMESH>(''_'',A8),
     *                   ''_IFDENS-'',I1,''_IFTENS-'',I1,
     *                   ''_IF-PAI-'',I1,
     *                   ''_IF-RAD-'',I1,''_IF-ORD-'',I1,
     *                   ''_IF-RHO-'',I1,
     *                   ''_IFDEEP-'',I1,''_IFPRON-'',I1,
     *                   ''_IFK-VC-'',I1,''_IFK-VS-'',I1,
     *                   ''_IFK-RC-'',I1,''_IFK-RS-'',I1,
     *                   ''_IFK-AC-'',I1,''_IFK-AS-'',I1,
     *                   ''_LDRAND-'',I5.5,''_'',A3,''.log'')') 
     *
     *                   FILNAM_NUCFIT(1:I1),FILNAM_ENDING,
     *                   (TITMSH(I),I=1,ICMESH),
     *                   IFDENS,IFTENS,IF_PAI,IF_RAD,IF_INV,
     *                   IF_RHO,IFDEEP,IFPRON,IFK_VC,IFK_VS,
     *                   IFK_RC,IFK_RS,IFK_AC,IFK_AS,LDRAND,
     *                                               VERSIO
C
          IF (MACTIV.GT.0) 
     *        WRITE(FILNAM,'(''LogFMesh/IFDENS-0/'',A,''_'',A1,
     *                   <ICMESH>(''_'',A8),
     *                   ''_IFDENS-'',I1,''_IFTENS-'',I1,
     *                   ''_IF-PAI-'',I1,
     *                   ''_IF-RAD-'',I1,''_IF-ORD-'',I1,
     *                   ''_IF-RHO-'',I1,
     *                   ''_IFDEEP-'',I1,''_IFPRON-'',I1,
     *                   ''_IFK-VC-'',I1,''_IFK-VS-'',I1,
     *                   ''_IFK-RC-'',I1,''_IFK-RS-'',I1,
     *                   ''_IFK-AC-'',I1,''_IFK-AS-'',I1,
     *                   ''_LDRAND-'',I5.5,''_'',A,''_'',A3,''.log'')') 
     *
     *                   FILNAM_NUCFIT(1:I1),FILNAM_ENDING,
     *                   (TITMSH(I),I=1,ICMESH),
     *                   IFDENS,IFTENS,IF_PAI,IF_RAD,IF_INV,
     *                   IF_RHO,IFDEEP,IFPRON,IFK_VC,IFK_VS,
     *                   IFK_RC,IFK_RS,IFK_AC,IFK_AS,LDRAND,
     *                   FILNAM_NOTTIT(1:I2),VERSIO
      END IF
C
C=======================================================================
C     
      IF (IFDENS.EQ.1) THEN
C      
          IF (MACTIV.EQ.0) 
     *        WRITE(FILNAM,'(''LogFMesh/IFDENS-1_IFTENS-'',I1,''/'',
     *                      A,''_'',A1,<ICMESH>(''_'',A8),
     *                   ''_IFDENS-'',I1,''_IFTENS-'',I1,
     *                   ''_IF-PAI-'',I1,
     *                   ''_IF-RAD-'',I1,''_IF-ORD-'',I1,
     *                   ''_IF-RHO-'',I1,
     *                   ''_IFDEEP-'',I1,''_IFPRON-'',I1,
     *                   ''_IFK-VC-'',I1,''_IFK-VS-'',I1,
     *                   ''_IFK-RC-'',I1,''_IFK-RS-'',I1,
     *                   ''_IFK-AC-'',I1,''_IFK-AS-'',I1,
     *                   ''_LDRAND-'',I5.5,A15,''_'',A3,''.log'')') 
     *
     *                   IFTENS,FILNAM_NUCFIT(1:I1),FILNAM_ENDING,
     *                   (TITMSH(I),I=1,ICMESH),
     *                   IFDENS,IFTENS,IF_PAI,IF_RAD,IF_INV,
     *                   IF_RHO,IFDEEP,IFPRON,IFK_VC,IFK_VS,
     *                   IFK_RC,IFK_RS,IFK_AC,IFK_AS,LDRAND,TEXLAM,
     *                                                      VERSIO
C
          IF (MACTIV.GT.0) 
     *        WRITE(FILNAM,'(''LogFMesh/IFDENS-1_IFTENS-'',I1,''/'',
     *                      A,''_'',A1,<ICMESH>(''_'',A8),
     *                   ''_IFDENS-'',I1,''_IFTENS-'',I1,
     *                   ''_IF-PAI-'',I1,
     *                   ''_IF-RAD-'',I1,''_IF-ORD-'',I1,
     *                   ''_IF-RHO-'',I1,
     *                   ''_IFDEEP-'',I1,''_IFPRON-'',I1,
     *                   ''_IFK-VC-'',I1,''_IFK-VS-'',I1,
     *                   ''_IFK-RC-'',I1,''_IFK-RS-'',I1,
     *                   ''_IFK-AC-'',I1,''_IFK-AS-'',I1,
     *                   ''_LDRAND-'',I5.5,''_'',A,A15,''_'',A3,
     *                                                 ''.log'')') 
     *
     *                   IFTENS,FILNAM_NUCFIT(1:I1),FILNAM_ENDING,
     *                   (TITMSH(I),I=1,ICMESH),
     *                   IFDENS,IFTENS,IF_PAI,IF_RAD,IF_INV,
     *                   IF_RHO,IFDEEP,IFPRON,IFK_VC,IFK_VS,
     *                   IFK_RC,IFK_RS,IFK_AC,IFK_AS,LDRAND,
     *                   FILNAM_NOTTIT(1:I2),TEXLAM,VERSIO
C
          IF (IFTENS.EQ.1) THEN
C      
              IF (MACTIV.EQ.0) 
     *            WRITE(FILNAM,'(''LogFMesh/IFDENS-1_IFTENS-'',I1,''/'',
     *                      A,''_'',A1,<ICMESH>(''_'',A8),
     *                   ''_IFDENS-'',I1,''_IFTENS-'',I1,
     *                   ''_IF-PAI-'',I1,
     *                   ''_ICENTT-'',I1,''_ISORBT-'',I1,
     *                                   ''_ITENSR-'',I1,
     *                   ''_IF-RAD-'',I1,''_IF-ORD-'',I1,
     *                   ''_IF-RHO-'',I1,
     *                   ''_IFDEEP-'',I1,''_IFPRON-'',I1,
     *                   ''_IFK-VC-'',I1,''_IFK-VS-'',I1,
     *                   ''_IFK-RC-'',I1,''_IFK-RS-'',I1,
     *                   ''_IFK-AC-'',I1,''_IFK-AS-'',I1,
     *                   ''_LDRAND-'',I5.5,A15,''_'',A3,''.log'')') 
     *
     *                   IFTENS,FILNAM_NUCFIT(1:I1),FILNAM_ENDING,
     *                   (TITMSH(I),I=1,ICMESH),
     *                   IFDENS,IFTENS,IF_PAI,ICENTT,ISORBT,ITENSR,
     *                   IF_RAD,IF_INV,IF_RHO,IFDEEP,IFPRON,IFK_VC,
     *                   IFK_VS,IFK_RC,IFK_RS,IFK_AC,IFK_AS,LDRAND,
     *                                               TEXLAM,VERSIO
C
              IF (MACTIV.GT.0) 
     *            WRITE(FILNAM,'(''LogFMesh/IFDENS-1_IFTENS-'',I1,''/'',
     *                      A,''_'',A1,<ICMESH>(''_'',A8),
     *                   ''_IFDENS-'',I1,''_IFTENS-'',I1,
     *                   ''_IF-PAI-'',I1,
     *                   ''_ICENTT-'',I1,''_ISORBT-'',I1,
     *                                   ''_ITENSR-'',I1,
     *                   ''_IF-RAD-'',I1,''_IF-ORD-'',I1,
     *                   ''_IF-RHO-'',I1,
     *                   ''_IFDEEP-'',I1,''_IFPRON-'',I1,
     *                   ''_IFK-VC-'',I1,''_IFK-VS-'',I1,
     *                   ''_IFK-RC-'',I1,''_IFK-RS-'',I1,
     *                   ''_IFK-AC-'',I1,''_IFK-AS-'',I1,
     *                   ''_LDRAND-'',I5.5,''_'',A,A15,''_'',A3,
     *                                                 ''.log'')') 
     *
     *                   IFTENS,FILNAM_NUCFIT(1:I1),FILNAM_ENDING,
     *                   (TITMSH(I),I=1,ICMESH),
     *                   IFDENS,IFTENS,IF_PAI,ICENTT,ISORBT,ITENSR,
     *                   IF_RAD,IF_INV,IF_RHO,IFDEEP,IFPRON,IFK_VC,
     *                   IFK_VS,IFK_RC,IFK_RS,IFK_AC,IFK_AS,LDRAND,
     *                   FILNAM_NOTTIT(1:I2),TEXLAM,VERSIO
C
          END IF
C
      END IF       
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
           WRITE(0,'(/,''LEVNUM= '',I4,''and NDIM_M= '',I4,
     *               ''--> LEVNUM should be .LT. NDIM_M!!'',/)')
           STOP 'STOP in LMMINI: LEVNUM.GT.NDIM_M!!'
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
      NUMPAR=IACTIV
      INFOER=111
C
C=======================================================================
C    
CID      NEWSED=0

      IF (ISCREN.EQ.1) THEN
          WRITE(LSCREN,'(/,101(''=''),/,101(''=''),/,
     *                     ''Point No.'',i4,''/'',i4,/,101(''=''))')
     *                                   NCOUNT_OFMESH,NTOTAL_OFMESH
      END IF
      
CID      IF (LOGWRI.GT.0) THEN
          WRITE(LOGAUX,'(/,101(''=''),/,101(''=''),/,
     *                     ''Point No.'',i4,''/'',i4,/,101(''=''))')
     *                                   NCOUNT_OFMESH,NTOTAL_OFMESH
CID      END IF
C      
C=======================================================================
C
      ITRMAX=30
CID      LDRAND=ITRMAX
      TESTGR=3.00
C
      CHITOT_MINAUX=1.0E+10
C
      IRANDO=0
      
CID      DO IRANDO=1,LDRAND
C
  100 CONTINUE
C
      IRANDO=IRANDO+1
C     
      IFUNCT_EVALUS=0
         
      IDEFCN=1
C      
C=======================================================================
C
      DO I=1,MAXFEV
         HIAUXI(I)=0.0
      END DO
C      
C=======================================================================
C
      IF (ISCREN.EQ.1) THEN
          WRITE(LSCREN,'(/,58X,''New run: NEWSED='',i14,3x,
     *                         ''No.'',i4,''/'',i4,/,101(''=''),/)')
     *                              NEWSED,IRANDO,ITRMAX
      END IF
C
      IF (LOGWRI.GT.0) THEN
C
          WRITE(LOGFIL,'(//,9X,''Entering NEWRAN loop with IRANDO='',
     *                      I3,1X,''LDRAND='',I3,'' from LMMESH'')')
     *                              IRANDO,ITRMAX
      END IF
C
      WRITE(LOGAUX,'(/,58X,''New run: NEWSED='',i14,3x,
     *                     ''No.'',i4,''/'',i4,/,101(''=''),/)')
     *                              NEWSED,IRANDO,ITRMAX
C      
C=======================================================================
C                  
      CALL NEWRAN(NEWSED,IZ_FIX,IN_FIX,ISOSPI,ARGPAR,DIRNAM,
     *                                 GENAME,IFIRST,I_SEED)
C      
C=======================================================================
C        Calling the principal minimisation routine based on the
C        Levenberg-Marquardt algorithm (Argonne MINPACK project)  
C=======================================================================
C
      NDFUNC=NDIM_M    ! NDIM_M
      LDFUNC=LEVNUM    ! Actual number of functions, here: levels
C     NDPARS=NDPARS=48 ! ARGPAR(1:NDPARS)
      LDPARS=IACTIV
C 
C     MAXFEV=900       ! MAXFEV
      IFMODE=INMODE
C     NPRINT=-1        ! NPRINT
      INFRUN=INFOER    ! In any case output variable
      NFCALL=0         !NFUNEV    ! In any case output variable
      NJCALL=0         !NJACEV    ! In any case output variable
C
C     I_PERM=IPIVOT    ! In any case output vector
C
C                        
C     FJACOB(1:NDFUNC,1:NDPARS) <<=>> FUNJAC(1:NDIM_M,1:NDPARS)
C     FARGUN(1:NDFUNC)          <<=>> FUNVEC(1:NDIM_M)
C     DIAGSC(1:NDPARS)          <<=>> SCALFC(1:NDPARS)
C     I_PERM(1:NDPARS)          <<=>> IPIVOT(1:NDPARS)
C                        
      CALL LEVMAR(FUNMIN,NDFUNC,LDFUNC,NDPARS,LDPARS,ARGPAR,
     *            FUNVEC,FUNJAC,TOLERF,TOLERX,TOLERG,MAXFEV,
     *            SCALFC,IFMODE,FACTOR,NPRINT,INFRUN,NFCALL,
     *            NJCALL,IPIVOT,QTFARR,WORKA1,WORKA2,WORKA3,
     *                                        WORKA4,NDLAST)
C
      INFOER=INFRUN
C
      JACTIV=0
      DO IPARAM=1,NDPARS
         IF (IFTAKE(IPARAM).EQ.1) THEN
             JACTIV=JACTIV+1
             PARPOT(IPARAM)=ARGPAR(JACTIV)
         END IF
      END DO
C
C=======================================================================
C         
      WRITE(LOGAUX,'(/,''Number of WS_RUN Evaluations: '',I5,/)')
     *                      IFUNCT_EVALUS
C
C=======================================================================
C      
      CALL TELLIT(INFOER,TOLERF,TOLERX,TOLERG,MAXFEV,LDLAST,EPSLAS)
C
C=======================================================================
C
      IF (LOGWRI.GT.0) THEN
          WRITE(LOGFIL,'(9X,''Entering JMATRX from LMMESH'')')
      END IF
C
cid   CALL JMATRX(LDFUNC,ARGPAR)
C
C=======================================================================
C
      IF (CHITOT_MINAUX.GT.VALCHI) THEN
C
          CHITOT_MINAUX=VALCHI
          CHISQU_PROAUX=CHISQU_PROTON
          CHISQU_NEUAUX=CHISQU_NEUTRS
C
          IF (IFDENS.EQ.0 .AND. IFPROT.EQ.1) THEN
              CHISQU_PRINTI(NCOUNT_OFMESH)=CHISQU_PROTON
          END IF
C
          IF (IFDENS.EQ.0 .AND. IFNEUT.EQ.1) THEN
              CHISQU_PRINTI(NCOUNT_OFMESH)=CHISQU_NEUTRS
          END IF
C
          IF (IFDENS.EQ.1 .OR. (IFDENS.EQ.0 .AND. IFBOTH.EQ.1)) THEN
              CHISQU_PRINTI(NCOUNT_OFMESH)=VALCHI
          END IF
C
          CHIGRD_NRMPRI(NCOUNT_OFMESH)=CHISQU_GRDNRM
C
          DO IPARAM=1,NDPARS
C
             PARPOT_PRINTI(IPARAM,NCOUNT_OFMESH)=PARPOT(IPARAM)
C
             CHIGRD_COMPRI(IPARAM,NCOUNT_OFMESH)
     *      =CHISQU_GRADIE(IPARAM)
C
          END DO
C
          IDEFCN_PRINTI=IDEFCN
C
          INFOER_PRINTI=INFOER
C
      END IF
C
C=======================================================================
C
      IF (IRANDO.LT.ITRMAX .AND. CHISQU_GRDNRM.GT.(3*TESTGR)) THEN
          I_SEED=NEWSED
          GO TO 100
      END IF
C
CID      END DO !IRANDO
C      
C=======================================================================
C      
      IF (LOGWRI.GT.4) WRITE(LOGFIL,'(''Entering WRITIN_MESHPM '',
     *                                                ''from LMMESH'')')
      CALL WRITIN_MESHPM(NCOUNT_OFMESH,PARPOT_PRINTI)
C      
C=======================================================================
C 
      REWIND LOGPRO    
C
   10 CONTINUE
C   
      READ (LOGPRO,'(A120)',END=11) STRING
      IF (LOGWRI.GT.4) WRITE(NpUNIT,'(A80)'       ) STRING
C      
      GO TO 10
   11 CONTINUE  ! Ready with the copy for the protons      
C      
C=======================================================================
C 
      REWIND LOGNEU    
C
   20 CONTINUE
C   
      READ (LOGNEU,'(A120)',END=21) STRING
      IF (LOGWRI.GT.4) WRITE(NnUNIT,'(A80)'       ) STRING
C      
      GO TO 20
   21 CONTINUE  ! Ready with the copy for the neutrons     
C
C=======================================================================
C
      IF (LOGWRI.GT.0) THEN
          WRITE(LOGFIL,'(/,18X,''Exiting  LMMESH'')')
      END IF   
C
C=======================================================================
C 
      CALL CPUTIM('LMMESH',0)
C
C=======================================================================
C      
      RETURN
      END
C
C=======================================================================
C=======================================================================
C
      SUBROUTINE TBMESH(IPARA1,IPARA2,INDEX1,INDEX2,LEVNUM,
     *                         CHITOT_MINAUX,CHISQU_PROAUX,
     *                         CHISQU_NEUAUX,PARPOT_PRINTI,
     *                         RMSMAP_PRONUC,RMSMAP_NEUNUC,
     *                         RMSMAP_PROGLO,RMSMAP_NEUGLO,
     *                         RMSPRO_MAXIMA,RMSPRO_MINIMI,
     *                         RMSNEU_MAXIMA,RMSNEU_MINIMI)
          
      INCLUDE   'MATDIM/NDNUCL.f'
      INCLUDE   'MATDIM/NDIM_P.f'
      INCLUDE   'MATDIM/NDPARS.f'
      INCLUDE   'MATDIM/NDMESH.f'
      INCLUDE   'MATDIM/NDLEXP.f'
      
      PARAMETER
     *          (NDMES2=NDMESH*NDMESH)
      CHARACTER  
     *          INPSYM*6,NUCSYM*6,TYPCHI*6
      DIMENSION
     *          PARPOT_PRINTI(1:NDPARS,1:NDMES2)
      DIMENSION
     *          RMSMAP_PRONUC(1:NDMESH,1:NDMESH,1:NDNUCL),
     *          RMSMAP_NEUNUC(1:NDMESH,1:NDMESH,1:NDNUCL)
      DIMENSION
     *          RMSMAP_PROGLO(1:NDMESH,1:NDMESH,1:NDNUCL),
     *          RMSMAP_NEUGLO(1:NDMESH,1:NDMESH,1:NDNUCL)
      COMMON
     *       /TESFIT/ IFTEST,IFFITS,IFMESH
      COMMON
     *       /DENSTR/ IFDENS
      COMMON
     *       /TENSOR/ IFTENS
      COMMON
     *       /TENOPT/ ICENTT,ISORBT,ITENSR
      COMMON
     *       /P_N_LM/ IFPROT,IFNEUT,IFBOTH
      COMMON
     *       /POTPOT/ PARPOT(1:NDIM_P)
      COMMON
     *       /ARGINI/ XINITS(1:NDIM_P)
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
      COMMON
     *       /DIMSPH/ LDSING_PROTON,LEVTHE_PROTON,
     *                LEVEXP_PROTON(1:NDNUCL),
     *                LDSING_NEUTRS,LEVTHE_NEUTRS,
     *                LEVEXP_NEUTRS(1:NDNUCL)
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
     *       /CHITAK/ IF_SPE,IF_RAD,IF_GAP,IF_FER,IF_DEN,IF_RHO,IF_INV,
     *                LDCHI2,TYPCHI(1:NDNUCL),ITKCHI(1:NDNUCL)
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
     *       /PRFILE/ ISCREN,LOGFIL,LOGCHP,LOGCHN,LOGHIP,LOGHIN,
     *                LOGPRO,LOGNEU,LOGWRI,LSCREN,IRESUL,ICONVE,
     *                NOUTPT,LOGAUX,LOGENE,LOGRAD,LOGGRD,LOGMSH,
     *                                                   LOGBIS  
C
C=======================================================================
C @@@
C     This SUBROUTINE what exactly???
C
C=======================================================================
C
      LEVNUM=0
       
      IF (IFDENS.EQ.1) THEN
       
          DO INUCLI=1,LDNUCL
        
             IF (ITAKNU(INUCLI).EQ.1) THEN
            
        
         IF (IF_SPE.EQ.1) LEVNUM=LEVNUM+LEVEXP_PROTON(INUCLI)
     *                                 +LEVEXP_NEUTRS(INUCLI)
     
         IF (IF_RAD.EQ.1) LEVNUM=LEVNUM+2 !Two radii for each nucleus (P+N)
C                                                              
         IF (IF_GAP.EQ.1) LEVNUM=LEVNUM+2 !Two Gap Energies for each nucleus
                                          !                            (P+N)
C                                                     
         IF (IF_FER.EQ.1) LEVNUM=LEVNUM+2 !Two Fermi Energies for each nucleus
                                          !                              (P+N)
        
         IF (IF_DEN.EQ.1) LEVNUM=LEVNUM+4 !Four densities for each nucleus: 
                                          !below and above the shell closure
                                          !                            (P+N)
C        IF (IF_RHO.EQ.1)
C        IF (IF_INV.EQ.1)
            
         END IF
        
        END DO
       
       END IF
C_______________________________________________________________________
C       
       IF ((IFDENS.EQ.0).AND.(IFPROT.EQ.1)) THEN 
           
        DO INUCLI=1,LDNUCL
        
         IF (ITAKNU(INUCLI).EQ.1) THEN
         
          IF (IF_SPE.EQ.1) LEVNUM=LEVNUM+LEVEXP_PROTON(INUCLI)
     
          IF (IF_RAD.EQ.1) LEVNUM=LEVNUM+1 !Proton radius
C                                                              
          IF (IF_GAP.EQ.1) LEVNUM=LEVNUM+1 !Proton Gap energy
C                                                     
          IF (IF_FER.EQ.1) LEVNUM=LEVNUM+1 !Proton Fermi energy
        
          IF (IF_DEN.EQ.1) LEVNUM=LEVNUM+2 !Proton up and down density energies

C         IF (IF_RHO.EQ.1)
C         IF (IF_INV.EQ.1)
            
         END IF
        
        END DO
       
       END IF
C_______________________________________________________________________
C       
       IF ((IFDENS.EQ.0).AND.(IFNEUT.EQ.1)) THEN 
       
        DO INUCLI=1,LDNUCL
        
         IF (ITAKNU(INUCLI).EQ.1) THEN
         
          IF (IF_SPE.EQ.1) LEVNUM=LEVNUM+LEVEXP_NEUTRS(INUCLI)
     
          IF (IF_RAD.EQ.1) LEVNUM=LEVNUM+1 !Neutron radius
C                                                              
          IF (IF_GAP.EQ.1) LEVNUM=LEVNUM+1 !Neutron Gap energy
C                                                     
          IF (IF_FER.EQ.1) LEVNUM=LEVNUM+1 !Neutron Fermi energy
        
          IF (IF_DEN.EQ.1) LEVNUM=LEVNUM+2 !Neutron up and down density energies

C         IF (IF_RHO.EQ.1)
C         IF (IF_INV.EQ.1)
            
         END IF
         
        END DO
       
       END IF       
C_______________________________________________________________________
C                
       DO IPARAM=1,NDPARS
          IF (IFMESH.EQ.1) THEN
              IF (I_MESH(IPARAM).EQ.1) THEN
                  PARPOT(IPARAM)=XINITS(IPARAM)
              END IF
          END IF
       END DO
                
       WRITE(0,'(''XINITS(1)= '',F9.3,'' ;XINITS(2)= '',F9.3)')
     *             XINITS(IPARA1),XINITS(IPARA2)
C                
       CHISQU_PROTON=0.
       CHISQU_NEUTRS=0.
       CHISQU_TOTALS=0.
C
       IDEG_P=0
       IDEG_N=0
       RMSGLO_TAKPRO=0.
       RMSGLO_TAKNEU=0.
C
       DO INUCLI=1,LDNUCL
C                    
          IF (ITAKNU(INUCLI).EQ.1) THEN
C           
              IZ_FIX=NUMB_Z(INUCLI)
              IN_FIX=NUMB_N(INUCLI)
C                    
              IDEFCN=0
              I_MODE=0
              I_FLAG=1
C_______________________________________________________________________
C
              IF (LOGWRI.GT.4) THEN
                  WRITE(LOGFIL,'(12X,''Entering WS_RUN from TBMESH '',
     *                               ''INUCLI='',I1,1X,
     *                               ''IZ_FIX='',I3,'' IN_FIX='',I3)')
     *                                 INUCLI,         IZ_FIX,IN_FIX
              END IF
C                     
              CALL WS_RUN(INUCLI,CHISQU_CHOICE,IABORT,I_MODE,
     *                    I_FLAG,CHISQU_AUXIL1,CHISQU_AUXIL2)
C     
              CHISQU_PROTON=CHISQU_PROTON+CHISQU_AUXIL1
              CHISQU_NEUTRS=CHISQU_NEUTRS+CHISQU_AUXIL2
              CHISQU_TOTALS=CHISQU_TOTALS+CHISQU_PROTON
     *                                   +CHISQU_NEUTRS
C
C             Storing the rms values of single nucleus
C
              RMSMAP_PRONUC(INDEX1,INDEX2,INUCLI)
     *       =
     *        SQRT(CHIWEI_PROTON)
C
              RMSMAP_NEUNUC(INDEX1,INDEX2,INUCLI)
     *       =
     *        SQRT(CHIWEI_NEUTRS)
C
C             Total degenerancy of the active nuclei (for the global)
C
              DO IEXPER=1,LEVEXP_PROTON(INUCLI)
                 IDEG_P=IDEG_P+IDEGEX_PROTON(INUCLI,IEXPER)
              END DO
C
              DO IEXPER=1,LEVEXP_NEUTRS(INUCLI)
                 IDEG_N=IDEG_N+IDEGEX_NEUTRS(INUCLI,IEXPER)
              END DO
C
              RMSGLO_TAKPRO=RMSGLO_TAKPRO+CHIDEG_PROTON
              RMSGLO_TAKNEU=RMSGLO_TAKNEU+CHIDEG_NEUTRS
C
C             Searching the RMSMAP maximum/minimum values
C
              IF (RMSPRO_MAXIMA.LT.
     *            RMSMAP_PRONUC(INDEX1,INDEX2,INUCLI)) THEN
                  RMSPRO_MAXIMA=RMSMAP_PRONUC(INDEX1,INDEX2,INUCLI)
              END IF
C
              IF (RMSPRO_MINIMI.GT.
     *            RMSMAP_PRONUC(INDEX1,INDEX2,INUCLI)) THEN
                  RMSPRO_MINIMI=RMSMAP_PRONUC(INDEX1,INDEX2,INUCLI)
              END IF
C
              IF (RMSNEU_MAXIMA.LT.
     *            RMSMAP_NEUNUC(INDEX1,INDEX2,INUCLI)) THEN
                  RMSNEU_MAXIMA=RMSMAP_NEUNUC(INDEX1,INDEX2,INUCLI)
              END IF
C
              IF (RMSNEU_MINIMI.GT.
     *            RMSMAP_NEUNUC(INDEX1,INDEX2,INUCLI)) THEN
                  RMSNEU_MINIMI=RMSMAP_NEUNUC(INDEX1,INDEX2,INUCLI)
              END IF
C
          END IF
C
       END DO
C
       RMSMAP_PROGLO(INDEX1,INDEX2,LDNUCL+1)=SQRT(RMSGLO_TAKPRO/IDEG_P)
       RMSMAP_NEUGLO(INDEX1,INDEX2,LDNUCL+1)=SQRT(RMSGLO_TAKNEU/IDEG_N)
C
       CHITOT_MINAUX=CHISQU_TOTALS
       CHISQU_PROAUX=CHISQU_PROTON
       CHISQU_NEUAUX=CHISQU_NEUTRS
C
C=======================================================================
C           
      RETURN
      END 
C
C=======================================================================
C=======================================================================
C
      SUBROUTINE LMMINI_MONTEC(IMONTE,LDMONT,IACTIV,ARGPAR,
     *                         LEVNUM,NDLAST,NEWSED,I_SEED)
C
      INCLUDE   'MATDIM/NDSPEC.f'
      INCLUDE   'MATDIM/NDPARS.f'
      INCLUDE   'MATDIM/NDNUCL.f'
      INCLUDE   'MATDIM/NDLEXP.f'
      INCLUDE   'MATDIM/NDIM_M.f'
      INCLUDE   'MATDIM/MAXFEV.f'
C
      EXTERNAL 
     *          FUNMIN
C      
      PARAMETER 
     *         (LENGTH=256,INMODE=1,NPRINT=-1)
C
      CHARACTER
     *          LABPRO*6,LABNEU*6,LABPRO_PRINTG*006,LABNEU_PRINTG*6
      CHARACTER
     *          DIRNAM*256,GENAME*256,STRING*126
      CHARACTER
     *          LABTHE_PROTON*6,LABTHE_NEUTRS*6,
     *          LABEXP_PROTON*6,LABEXP_NEUTRS*6
      CHARACTER
     *          ACTION*1,CALLED*6,TYPCHI*6
      CHARACTER
     *          FILNAM*256,FILNAM_ENDING*13,INPSYM*6,
     *          LABTHE*006,LABTHE_PRINTG*06,TEXLAM*20,
     *          FILNA2*256,NUCSYM*6,FITNUC*2,FILNAM_NUCFIT*30,
     *          FITNUC_AUXILI*30,FILNA3*256,VERSIO*3
C    
      DIMENSION 
     *          ARGPAR(1:NDPARS),FUNVEC(1:NDIM_M),
     *          SCALFC(1:NDPARS),QTFARR(1:NDPARS)
      DIMENSION
     *          WORKA1(1:NDPARS),WORKA2(1:NDPARS),
     *          WORKA3(1:NDPARS),WORKA4(1:NDIM_M)
      DIMENSION 
     *          FUNJAC(1:NDIM_M,1:NDPARS)
      DIMENSION 
     *          IPIVOT(1:NDPARS)
      DIMENSION
     *          ENETHE(1:NDSPEC),LABTHE(1:NDSPEC)
      DIMENSION
     *          ENEPRO_PRINTG(1:NDSPEC),LABPRO_PRINTG(1:NDSPEC),
     *          ENENEU_PRINTG(1:NDSPEC),LABNEU_PRINTG(1:NDSPEC),
     *          ENETHE_PRINTG(1:NDSPEC),LABTHE_PRINTG(1:NDSPEC),
     *                                  PARPOT_PRINTG(1:NDPARS)
      DIMENSION
     *          ENEPRO(1:NDSPEC),LABPRO(1:NDSPEC),
     *          ENENEU(1:NDSPEC),LABNEU(1:NDSPEC)
      DIMENSION
     *          FITNUC(1:NDNUCL),
     *          FITNUC_AUXILI(1:NDNUCL)
C
      COMMON
     *       /DIMSPH/ LDSING_PROTON,LEVTHE_PROTON,
     *                              LEVEXP_PROTON(1:NDNUCL),
     *                LDSING_NEUTRS,LEVTHE_NEUTRS,
     *                              LEVEXP_NEUTRS(1:NDNUCL)
      COMMON
     *       /CHI_PN/ CHISQU_PROTON,CHISQU_NEUTRS
      COMMON
     *       /CHICHO/ VALCHI
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
     *       /TENOPT/ ICENTT,ISORBT,ITENSR
     *       /IFPAIR/ IF_PAI,IFGPAI,IFDPAI
      COMMON
     *       /CHIGRD/ CHISQU_GRDNRM,
     *                CHISQU_GRADIE(1:NDPARS)
      COMMON
     *       /STOPIN/ HIAUXI(1:MAXFEV)
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
     *       /PRFILE/ ISCREN,LOGFIL,LOGCHP,LOGCHN,LOGHIP,LOGHIN,
     *                LOGPRO,LOGNEU,LOGWRI,LSCREN,IRESUL,ICONVE,
     *                NOUTPT,LOGAUX,LOGENE,LOGRAD,LOGGRD,LOGMSH,
     *                                                   LOGBIS    
      COMMON
     *       /TOLERS/ TOLERF,TOLERX,TOLERG,FACTOR  
      COMMON  
     *       /VERSIN/ VERSIO
C
C=======================================================================
C     This subroutine minimizes the \chi^2 over an external do-loop
C     od IMONTE=1,LDMONT
C=======================================================================
C
      NUMPAR=IACTIV
      INFOER=111
C
      IRANDO=IMONTE
C
C=======================================================================
C    
      IFUNCT_EVALUS=0
C     
      IDEFCN=1
C
      DO I=1,MAXFEV
         HIAUXI(I)=0.0
      END DO
C         
      IF (ISCREN.EQ.1) THEN
C
          WRITE(LSCREN,'(/,58X,''New run: NEWSED='',i14,3x,
     *                         ''No.'',i6,''/'',i6,/,101(''=''),/)') 
     *                                     NEWSED,IMONTE,LDMONT
      END IF
C
      WRITE(LOGAUX,'(/,58X,''New run: NEWSED='',i14,3x,
     *                     ''No.'',i6,''/'',i6,/,101(''=''),/)') 
     *                                     NEWSED,IMONTE,LDMONT
C         
      IF (LOGWRI.GT.0) THEN
C             
          WRITE(LOGENE,'(/,58X,''New run: NEWSED='',i14,3x,
     *                            ''No.'',i6,''/'',i6,/,101(''=''),/)') 
     *                                     NEWSED,IMONTE,LDMONT
C             
          WRITE(LOGRAD,'(/,58X,''New run: NEWSED='',i14,3x,
     *                            ''No.'',i6,''/'',i6,/,101(''=''),/)') 
     *                                     NEWSED,IMONTE,LDMONT
      END IF
C
      WRITE(NOUTPT,'(/,80(''#''),/,''#'',T80,''#'',/,
     *                  ''#'',5X,''IMONTE= '',I6,T80,''#'',/,
     *                  ''#'',T80,''#'',/,80(''#''))') IMONTE
C      
C=======================================================================
C
C     Generating random restarts
C
      IF (LOGWRI.GT.0) THEN
          WRITE(LOGFIL,'(//,9X,''Entering NEWRAN with IRANDO='',I3,
     *                      1X,''LDRAND='',I6,'' from LMMINI'')')
     *                              IMONTE,LDMONT
      END IF
C         
      CALL NEWRAN(NEWSED,IZ_FIX,IN_FIX,ISOSPI,ARGPAR,DIRNAM,
     *                                 GENAME,IFIRST,I_SEED)
C      
C=======================================================================
C     Calling the principal minimisation routine based on the
C     Levenberg-Marquardt algorithm (Argonne MINPACK project)  
C=======================================================================
C
      NDFUNC=NDIM_M    ! NDIM_M
      LDFUNC=LEVNUM    ! Actual number of functions, here: levels
C IRENE? NDPARS=NDPARS=48 ! ARGPAR(1:NDPARS)
      LDPARS=IACTIV
C 
C        MAXFEV=900       ! MAXFEV
C @@@ ??
      IFMODE=INMODE
C     NPRINT=-1        ! NPRINT
C @@@ WHAT DOES THIS MEAN???
      INFRUN=INFOER    ! In any case output variable
      NFCALL=0         ! NFUNEV    ! In any case output variable
      NJCALL=0         ! NJACEV    ! In any case output variable
C
C     I_PERM=IPIVOT    ! In any case output vector ?????????????
C                     
C     FJACOB(1:NDFUNC,1:NDPARS) <<=>> FUNJAC(1:NDIM_M,1:NDPARS)
C     FARGUN(1:NDFUNC) <<=>> FUNVEC(1:NDIM_M)
C     DIAGSC(1:NDPARS) <<=>> SCALFC(1:NDPARS)       
C     I_PERM(1:NDPARS) <<=>> IPIVOT(1:NDPARS)       
C
C     Calling the Levenberg-Marquardt LAPACK minimisation routine 
C
      IF (LOGWRI.GT.0) THEN
          WRITE(LOGFIL,'(/,6X,''Entering LEVMAR from LMMINI      '',
     *                   5(3x,''IMONTE='',I2))')
     *                          IMONTE,IMONTE,IMONTE,IMONTE,IMONTE
      END IF
C
      CALL LEVMAR(FUNMIN,NDFUNC,LDFUNC,NDPARS,LDPARS,ARGPAR,
     *            FUNVEC,FUNJAC,TOLERF,TOLERX,TOLERG,MAXFEV,
     *            SCALFC,IFMODE,FACTOR,NPRINT,INFRUN,NFCALL,
     *            NJCALL,IPIVOT,QTFARR,WORKA1,WORKA2,WORKA3,
     *                                        WORKA4,NDLAST)
C
      INFOER=INFRUN 
C
      KK=0
      DO KPARAM=1,NDPARS
         IF (IFTAKE(KPARAM).EQ.1) THEN
             KK=KK+1
             PARPOT(KPARAM)=ARGPAR(KK)
         END IF
      END DO     
C
C=======================================================================
C
      WRITE(LOGAUX,'(/,''Number of WS_RUN Evaluations: '',I5,/)')
     *                      IFUNCT_EVALUS
C
C=======================================================================
C      
      CALL TELLIT(INFOER,TOLERF,TOLERX,TOLERG,MAXFEV,LDLAST,EPSLAS)
C
C=======================================================================
C
      CALL RADINF(IFPROT,IFNEUT)
C
C=======================================================================
C
      IF (LOGWRI.GT.0) THEN
          WRITE(LOGFIL,'(9X,''Entering JMATRX from LMMINI'')')
      END IF
C
      CALL JMATRX(LDFUNC,ARGPAR)
C
C=======================================================================
C           
      RETURN
      END 
C
C=======================================================================
C=======================================================================
C 
      SUBROUTINE NEWRAN(NEWSED,IZ_FIX,IN_FIX,ISOSPI,ARGPAR,DIRNAM,
     *                                       GENAME,IFIRST,I_SEED)
      INCLUDE  'MATDIM/NDPARS.f'
      INCLUDE  'MATDIM/NDRAND.f'
C
      CHARACTER
     *          FILKEE*256,WHATEX*6
      CHARACTER
     *          DIRNAM*256,GENAME*256,STRAUX*256
C
      DIMENSION
     *          ARGPAR(1:NDPARS)
C
      COMMON
     *       /OLDVAL/ OVCENT,OXSORB,OXEFFM,ORCENT,ORSORB,OREFFM,
     *                                     OACENT,OASORB,OAEFFM,
     *                                                   O0COUL
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
     *       /PRFILE/ ISCREN,LOGFIL,LOGCHP,LOGCHN,LOGHIP,LOGHIN,
     *                LOGPRO,LOGNEU,LOGWRI,LSCREN,IRESUL,ICONVE,
     *                NOUTPT,LOGAUX,LOGENE,LOGRAD,LOGGRD,LOGMSH,
     *                                                   LOGBIS
      COMMON
     *       /CHIMIN/ HIWMIN,HICMIN,ERRMIN,SHEMIN,SHWMIN
      COMMON
     *       /CNTROL/ PARAMT(1:NDPARS),DPARAM(1:NDPARS)     
      COMMON
     *       /RANRAN/ IRMFST,I_RAND,IRANDO
      COMMON
     *       /SEARCH/ IFTAKE(1:NDPARS)
      COMMON
     *       /RANDSR/ ARANDO(1:NDPARS)
      COMMON
     *       /VLIMIT/ VMISTR(1:NDPARS),
     *                VMIMIN(1:NDPARS),
     *                VMIMAX(1:NDPARS),
     *                VMISTP(1:NDPARS)
      COMMON
     *       /VSTART/ STARTV(1:NDPARS)
      COMMON
     *       /SIZMIN/ LDSEAR
      COMMON
     *       /PRINAL/ IMTALK
      COMMON                
     *       /INPSEL/ WHATEX
C
C=======================================================================
C     This subroutine defines a complete set of random parameters;
C     It should be called at the beginning of one minimisation run 
C=======================================================================
C
      DIRNAM='Results'
C
C=======================================================================
C     Initialising the chi^2 functions for minimum calculations
C=======================================================================
C
      HIWMIN=1.0E+10
      HICMIN=1.0E+10
      ERRMIN=1.0E+10
      SHEMIN=1.0E+10
      SHWMIN=1.0E+10
C
C=======================================================================
C=======================================================================
C     Defining the filename that will contain the results for a
C     given initial (eventually) random set of Dirac parameters
C=======================================================================
C=======================================================================
C
      GENAME='ResFile' ! Write unit is called IRESUL  
C
      CALL STRLNG(DIRNAM,LENGTH,STRAUX)
C
      WRITE(FILKEE(1:LENGTH),'(A)') DIRNAM(1:LENGTH)
C
      LENGTH=LENGTH+1
C
C=======================================================================
C
      IF (ISOSPI.EQ.0) THEN
C
          WRITE(FILKEE(LENGTH:LENGTH+15),'(''/WS_Z'',
     *          I3.3,''_N'',I3.3,''_n_'')') IZ_FIX,IN_FIX
      END IF
C
      IF (ISOSPI.EQ.1) THEN
C
          WRITE(FILKEE(LENGTH:LENGTH+15),'(''/WS_Z'',
     *          I3.3,''_N'',I3.3,''_p_'')') IZ_FIX,IN_FIX
      END IF
C
C=======================================================================
C
      LENGTH=LENGTH+16
C
      IF (I_RAND.EQ.1) THEN
C
          WRITE(FILKEE(LENGTH:LENGTH+9),'(''rand'',I4.4,''_'')') IRANDO
C
          LENGTH=LENGTH+9
C
          CALL STRLNG(GENAME,LNGGEN,STRAUX)
C
          WRITE(FILKEE(LENGTH:LENGTH+LNGGEN),'(A)') GENAME(1:LNGGEN)
          WRITE(FILKEE(LENGTH+LNGGEN:LENGTH+LNGGEN+4),'(''.out'')')
C
      ELSE
C
          CALL STRLNG(GENAME,LNGGEN,STRAUX)
C
          WRITE(FILKEE(LENGTH:LENGTH+LNGGEN),'(A)') GENAME(1:LNGGEN)
          WRITE(FILKEE(LENGTH+LNGGEN:LENGTH+LNGGEN+4),'(''.out'')')
C
      END IF
C
      LENGTH=LENGTH+LNGGEN+4
C
C=======================================================================
C
      IF (LOGWRI.GT.0) THEN
C
          WRITE(LOGFIL,'(18X,''Opening UNIT=IRESUL='',I2,1X,
     *                       ''File:'',A)') IRESUL,FILKEE(1:LENGTH)
      END IF 
C
CID      OPEN (IRESUL,FILE=FILKEE(1:LENGTH),STATUS='UNKNOWN',
CID     *                                   FORM='FORMATTED')
C
C=======================================================================
C=======================================================================
C     Defining the filename that will contain the logfile of the
C     minimisation series  (these are useful at the begin stage)
C=======================================================================
C=======================================================================
C
      GENAME='ConvRep' ! Write unit is called LOGFIL   
C
      CALL STRLNG(DIRNAM,LENGTH,STRAUX)
C
      WRITE(FILKEE(1:LENGTH),'(A)') DIRNAM(1:LENGTH)
C
      LENGTH=LENGTH+1
C
C=======================================================================
C
      IF (ISOSPI.EQ.0) THEN
C
          WRITE(FILKEE(LENGTH:LENGTH+15),'(''/WS_Z'',
     *          I3.3,''_N'',I3.3,''_n_'')') IZ_FIX,IN_FIX
      END IF
C
      IF (ISOSPI.EQ.1) THEN
C
          WRITE(FILKEE(LENGTH:LENGTH+15),'(''/WS_Z'',
     *          I3.3,''_N'',I3.3,''_p_'')') IZ_FIX,IN_FIX
      END IF
C
C=======================================================================
C
      LENGTH=LENGTH+16
C
      IF (I_RAND.EQ.1) THEN
C
          WRITE(FILKEE(LENGTH:LENGTH+9),'(''rand'',I4.4,''_'')') IRANDO
C
          LENGTH=LENGTH+9
C
          CALL STRLNG(GENAME,LNGGEN,STRAUX)
C
          WRITE(FILKEE(LENGTH:LENGTH+LNGGEN),'(A)') GENAME(1:LNGGEN)
          WRITE(FILKEE(LENGTH+LNGGEN:LENGTH+LNGGEN+4),'(''.dat'')')
C
      ELSE
C
          CALL STRLNG(GENAME,LNGGEN,STRAUX)
C
          WRITE(FILKEE(LENGTH:LENGTH+LNGGEN),'(A)') GENAME(1:LNGGEN)
          WRITE(FILKEE(LENGTH+LNGGEN:LENGTH+LNGGEN+4),'(''.dat'')')
C
      END IF
C
      LENGTH=LENGTH+LNGGEN+4
C
C=======================================================================
C
      IF (LOGWRI.GT.0) THEN
C
          WRITE(LOGFIL,'(18X,''Opening UNIT=ICONVE='',I2,1X,
     *                       ''File:'',A)') ICONVE,FILKEE(1:LENGTH)
      END IF 
C
CID      OPEN (ICONVE,FILE=FILKEE(1:LENGTH),STATUS='UNKNOWN',
CID     *                                   FORM='FORMATTED')
C
      IF (ISCREN.GT.0) 
     *    WRITE(LSCREN,'(''Opening UNIT='',i2,'' File: '',A)')
     *                 ICONVE,            FILKEE(1:LENGTH)
C
      IF (LOGWRI.GT.0)
     *    WRITE(ICONVE,'(''Opening File: '',A)') FILKEE(1:LENGTH)
C
C=======================================================================
C
      CALL WARNIN(LOGFIL,WHATEX) ! Warning in case of the use
C                                  of 'not experimental' info
C
C=======================================================================
C
      IF (I_RAND.NE.1)  THEN
          IF (LOGWRI.GT.0) THEN
              WRITE(LOGFIL,'(''Exiting  NEWRAN - no action required'')')
          END IF 
          RETURN
      END IF 
C
C=======================================================================
C
      CALL RANDIN(I_SEED)        ! Initialize the seed
C
      CALL ZUFALL(NDPARS,ARANDO) ! Generates NDPARS random numbers 
C
C=======================================================================
C     First, determine the range of variation of each of the
C     parameters - then, select at RANDOM the starting point
C     variable. Each variable must be selected independently
C     of the other.                            =============
C=======================================================================
C
C     ... beginning with the central potential well, INDEX=1
C
C     The extremum (input) values of the parameters are kept
C     in the COMMON /EXTREM/
C
C     The random numbers in the range [0,1] are contained in
C     the COMMON /RANDSR/
C
      PARAMT(1)=V0CMIN+ARANDO(1)*(V0CMAX-V0CMIN)
C
      IF (IMTALK.EQ.1) THEN
C
          WRITE(IRESUL,'(''ARANDO(1) ='',F20.16,'' PARAMT(1) ='',F20.14,
     *                                         '' (OVCENT)'')')
     *                     ARANDO(1),              PARAMT(1)
      END IF
C
      IF (IFTAKE(1).EQ.1) OVCENT=PARAMT(1)
C_______________________________________________________________________
C
      PARAMT(11)=XKCMIN+ARANDO(11)*(XKCMAX-XKCMIN)
C
      IF (IMTALK.EQ.1) THEN
C
          WRITE(IRESUL,'(''ARANDO(11)='',F20.16,'' PARAMT(11)='',F20.14,
     *                                         '' (XK_V0C)'')')
     *                     ARANDO(11),             PARAMT(11)
      END IF
C @@@ CHECK CORRECTNESS
      IF (IFTAKE(11).EQ.1) OK_V0C=PARAMT(1)
C
C=======================================================================
C     Selecting STRENGTHS. The index values 4 and 7
C=======================================================================
C
      PARAMT(4)=XL_MIN+ARANDO(4)*(XL_MAX-XL_MIN)
C
      IF (IMTALK.EQ.1) THEN
C
          WRITE(IRESUL,'(''ARANDO(4) ='',F20.16,'' PARAMT(4) ='',F20.14,
     *                                         '' (OXSORB)'')')
     *                     ARANDO(4),              PARAMT(4)
      END IF
C
      IF (IFTAKE(4).EQ.1) OXSORB=PARAMT(4)
C_______________________________________________________________________
C
      PARAMT(14)=XKSMIN+ARANDO(14)*(XKSMAX-XKSMIN)
C
      IF (IMTALK.EQ.1) THEN
C
          WRITE(IRESUL,'(''ARANDO(14)='',F20.16,'' PARAMT(14)='',F20.14,
     *                                         '' (XK_LAM)'')')
     *                     ARANDO(14),              PARAMT(14)
      END IF
C  @@@ CHECK CORRECTNESS 
CID: OK!
      IF (IFTAKE(14).EQ.1) OK_XLM=PARAMT(14)
C
C=======================================================================
C
      PARAMT(7)=XEFMIN+ARANDO(7)*(XEFMAX-XEFMIN)
C
      IF (IRMFST.EQ.1) THEN
          PARAMT(7)=PARAMT(4)
      END IF
C
      IF (IMTALK.EQ.1) THEN
C
          WRITE(IRESUL,'(''ARANDO(7) ='',F20.16,'' PARAMT(7) ='',F20.14,
     *                                         '' (OXEFFM)'')')
     *                     ARANDO(7),              PARAMT(7)
      END IF
C      
      IF (IFTAKE(7).EQ.1) OXEFFM=PARAMT(7)
C_______________________________________________________________________
C
      PARAMT(17)=XKEMIN+ARANDO(17)*(XKEMAX-XKEMIN)
C
      IF (IRMFST.EQ.1) THEN
          PARAMT(17)=PARAMT(14)
      END IF
C
      IF (IMTALK.EQ.1) THEN
C
          WRITE(IRESUL,'(''ARANDO(17)='',F20.16,'' PARAMT(17)='',F20.14,
     *                                         '' (V0EFFM)'')')
     *                     ARANDO(17),              PARAMT(17)
      END IF
C      
      IF (IFTAKE(17).EQ.1) OEFMAS=PARAMT(17)
C
C=======================================================================
C     Selecting the RADII. The index values 2, 5 and 8
C=======================================================================
C
      PARAMT(2)=R0CMIN+ARANDO(2)*(R0CMAX-R0CMIN)
C
      IF (IMTALK.EQ.1) THEN
C
          WRITE(IRESUL,'(''ARANDO(2) ='',F20.16,'' PARAMT(2) ='',F20.14,
     *                                         '' (ORCENT)'')')
     *                     ARANDO(2),              PARAMT(2)
      END IF
C
      IF (IFTAKE(2).EQ.1) ORCENT=PARAMT(2)
C_______________________________________________________________________
C
      PARAMT(12)=XRCMIN+ARANDO(12)*(XRCMAX-XRCMIN)
C
      IF (IMTALK.EQ.1) THEN
C
          WRITE(IRESUL,'(''ARANDO(12)='',F20.16,'' PARAMT(12)='',F20.14,
     *                                         '' (XK_R0C)'')')
     *                     ARANDO(12),             PARAMT(12)
      END IF
C
      IF (IFTAKE(12).EQ.1) OK_R0C=PARAMT(12)
C
C=======================================================================
C
      PARAMT(5)=R0SMIN+ARANDO(5)*(R0SMAX-R0SMIN)
C
      IF (IMTALK.EQ.1) THEN
C
          WRITE(IRESUL,'(''ARANDO(5) ='',F20.16,'' PARAMT(5) ='',F20.14,
     *                                         '' (ORSORB)'')')
     *                     ARANDO(5),              PARAMT(5)
      END IF
C
      IF (IFTAKE(5).EQ.1) ORSORB=PARAMT(5)
C_______________________________________________________________________
C
      PARAMT(15)=XRSMIN+ARANDO(15)*(XRSMAX-XRSMIN)
C
      IF (IMTALK.EQ.1) THEN
C
          WRITE(IRESUL,'(''ARANDO(15)='',F20.16,'' PARAMT(15)='',F20.14,
     *                                         '' (XK_RSO)'')')
     *                     ARANDO(15),             PARAMT(15)
      END IF
C
      IF (IFTAKE(15).EQ.1) OK_RSO=PARAMT(15)
C
C=======================================================================
C
      PARAMT(8)=REFMIN+ARANDO(8)*(REFMAX-REFMIN)
C
      IF (IRMFST.EQ.1) THEN
          PARAMT(8)=PARAMT(5)
      END IF
C
      IF (IMTALK.EQ.1) THEN
C
          WRITE(IRESUL,'(''ARANDO(8) ='',F20.16,'' PARAMT(8) ='',F20.14,
     *                                         '' (OREFFM)'')')
     *                     ARANDO(8),              PARAMT(8)
      END IF
C
      IF (IFTAKE(8).EQ.1) OREFFM=PARAMT(8)
C_______________________________________________________________________
C
      PARAMT(18)=XREMIN+ARANDO(18)*(XREMAX-XREMIN)
C
      IF (IRMFST.EQ.1) THEN
          PARAMT(18)=PARAMT(15)
      END IF
C
      IF (IMTALK.EQ.1) THEN
C
          WRITE(IRESUL,'(''ARANDO(18)='',F20.16,'' PARAMT(18)='',F20.14,
     *                                         '' (XK_REF)'')')
     *                     ARANDO(18),             PARAMT(18)
      END IF
C
      IF (IFTAKE(18).EQ.1) OK_REF=PARAMT(18)
C
C=======================================================================
C     Selecting DIFFUSENESS; The index values 3, 6 and 9
C=======================================================================
C
      PARAMT(3)=A0CMIN+ARANDO(3)*(A0CMAX-A0CMIN)
C
      IF (IMTALK.EQ.1) THEN
C
          WRITE(IRESUL,'(''ARANDO(3) ='',F20.16,'' PARAMT(3) ='',F20.14,
     *                                         '' (OACENT)'')')
     *                     ARANDO(3),              PARAMT(3)
      END IF
C
      IF (IFTAKE(3).EQ.1) OACENT=PARAMT(3)
C_______________________________________________________________________
C
      PARAMT(13)=XACMIN+ARANDO(13)*(XACMAX-XACMIN)
C 
      IF (IMTALK.EQ.1) THEN
C
          WRITE(IRESUL,'(''ARANDO(13)='',F20.16,'' PARAMT(13)='',F20.14,
     *                                         '' (XK_A0C)'')')
     *                     ARANDO(13),             PARAMT(13)
      END IF
C
      IF (IFTAKE(13).EQ.1) OK_A0C=PARAMT(13)
C
C=======================================================================
C
      PARAMT(6)=A0SMIN+ARANDO(6)*(A0SMAX-A0SMIN)
C
      IF (IMTALK.EQ.1) THEN
C
          WRITE(IRESUL,'(''ARANDO(6) ='',F20.16,'' PARAMT(6) ='',F20.14,
     *                                         '' (OASORB)'')')
     *                     ARANDO(6),              PARAMT(6)
      END IF
C
      IF (IFTAKE(6).EQ.1) OASORB=PARAMT(6)
C_______________________________________________________________________
C
      PARAMT(16)=XASMIN+ARANDO(16)*(XASMAX-XASMIN)
C
      IF (IMTALK.EQ.1) THEN
C
          WRITE(IRESUL,'(''ARANDO(16)='',F20.16,'' PARAMT(16)='',F20.14,
     *                                         '' (XK_ASO)'')')
     *                     ARANDO(16),             PARAMT(16)
      END IF
C
      IF (IFTAKE(16).EQ.1) OK_ASO=PARAMT(16)
C
C=======================================================================
C
      PARAMT(9)=AEFMIN+ARANDO(9)*(AEFMAX-AEFMIN)
C
      IF (IRMFST.EQ.1) THEN
          PARAMT(9)=PARAMT(6)
      END IF
C
      IF (IMTALK.EQ.1) THEN
C
          WRITE(IRESUL,'(''ARANDO(9) ='',F20.16,'' PARAMT(9) ='',F20.14,
     *                                         '' (OAEFFM)'')')
     *                     ARANDO(9),              PARAMT(9)
      END IF
C
      IF (IFTAKE(9).EQ.1) OAEFFM=PARAMT(9)
C_______________________________________________________________________
C
      PARAMT(19)=XAEMIN+ARANDO(19)*(XAEMAX-XAEMIN)
C
      IF (IRMFST.EQ.1) THEN
          PARAMT(19)=PARAMT(16)
      END IF
C
      IF (IMTALK.EQ.1) THEN
C
          WRITE(IRESUL,'(''ARANDO(19)='',F20.16,'' PARAMT(19)='',F20.14,
     *                                         '' (XK_LEF)'')')
     *                     ARANDO(19),             PARAMT(19)
      END IF
C
      IF (IFTAKE(19).EQ.1) OK_XEF=PARAMT(19)
C
C=======================================================================
C         Selecting Coulomb radius; The index value 10
C=======================================================================
C
      PARAMT(10)=CouMIN+ARANDO(10)*(CouMAX-CouMIN)
C
      IF (IMTALK.EQ.1) THEN
C
          WRITE(IRESUL,'(''ARANDO(10)='',F20.16,'' PARAMT(10)='',F20.14,
     *                                         '' (O0COUL)'')')
     *                     ARANDO(10),             PARAMT(10)
      END IF
C
      IF (IFTAKE(10).EQ.1) O0COUL=PARAMT(10)
C_______________________________________________________________________
C
      PARAMT(20)=XCoMIN+ARANDO(20)*(XCoMAX-XCoMIN)
C
      IF (IMTALK.EQ.1) THEN
C
          WRITE(IRESUL,'(''ARANDO(20)='',F20.16,'' PARAMT(20)='',F20.14,
     *                                         '' (XK_COU)'')')
     *                     ARANDO(20),             PARAMT(20)
      END IF
C
      IF (IFTAKE(20).EQ.1) OK_COU=PARAMT(20)
C
C=======================================================================
C     Indices 21-37 are for neutrons...
C=======================================================================
C
      PARAMT(21)=V0CMIN+ARANDO(21)*(V0CMAX-V0CMIN)
C
      IF (IMTALK.EQ.1) THEN
C
          WRITE(IRESUL,'(''ARANDO(21)='',F20.16,'' PARAMT(21)='',F20.14,
     *                                         '' (OVCENT)'')')
     *                     ARANDO(21),              PARAMT(21)
      END IF
C @@@ WHY COMMENTED?
C      IF (IFTAKE(21).EQ.1) OVCENT=PARAMT(21)
C_______________________________________________________________________
C
      PARAMT(30)=XKCMIN+ARANDO(30)*(XKCMAX-XKCMIN)
C
      IF (IMTALK.EQ.1) THEN
C
          WRITE(IRESUL,'(''ARANDO(30)='',F20.16,'' PARAMT(30)='',F20.14,
     *                                         '' (XK_V0C)'')')
     *                     ARANDO(30),             PARAMT(30)
      END IF
C @@@ WHY COMMENTED???
      IF (IFTAKE(30).EQ.1) OK_V0C=PARAMT(30)
C
C=======================================================================
C     
C=======================================================================
C
      PARAMT(24)=XL_MIN+ARANDO(24)*(XL_MAX-XL_MIN)
C
      IF (IMTALK.EQ.1) THEN
C
          WRITE(IRESUL,'(''ARANDO(24)='',F20.16,'' PARAMT(24)='',F20.14,
     *                                          '' (OXSORB)'')')
     *                     ARANDO(24),              PARAMT(24)
      END IF
C
      IF (IFTAKE(24).EQ.1) OXSORB=PARAMT(24)
C_______________________________________________________________________
C
      PARAMT(33)=XKSMIN+ARANDO(33)*(XKSMAX-XKSMIN)
C
      IF (IMTALK.EQ.1) THEN
C
          WRITE(IRESUL,'(''ARANDO(33)='',F20.16,'' PARAMT(33)='',F20.14,
     *                                          '' (XK_LAM)'')')
     *                     ARANDO(33),              PARAMT(33)
      END IF
C @@@ VERIFY CORRECTNESS
CID: OK!
      IF (IFTAKE(33).EQ.1) OK_XLM=PARAMT(33)
C
C=======================================================================
C
      PARAMT(27)=XEFMIN+ARANDO(27)*(XEFMAX-XEFMIN)
C
      IF (IRMFST.EQ.1) THEN
          PARAMT(27)=PARAMT(24)
      END IF
C
      IF (IMTALK.EQ.1) THEN
C
          WRITE(IRESUL,'(''ARANDO(27)='',F20.16,'' PARAMT(27)='',F20.14,
     *                                         '' (OXEFFM)'')')
     *                     ARANDO(27),              PARAMT(27)
      END IF
C      
      IF (IFTAKE(27).EQ.1) OXEFFM=PARAMT(27)
C_______________________________________________________________________
C
      PARAMT(36)=XKEMIN+ARANDO(36)*(XKEMAX-XKEMIN)
C
      IF (IRMFST.EQ.1) THEN
          PARAMT(36)=PARAMT(33)
      END IF
C
      IF (IMTALK.EQ.1) THEN
C
          WRITE(IRESUL,'(''ARANDO(36)='',F20.16,'' PARAMT(36)='',F20.14,
     *                                         '' (V0EFFM)'')')
     *                     ARANDO(36),              PARAMT(36)
      END IF
C      
      IF (IFTAKE(36).EQ.1) OEFMAS=PARAMT(36)
C
C=======================================================================
C     Selecting the RADII. The index values 22, 5 and 8
C=======================================================================
C
      PARAMT(22)=R0CMIN+ARANDO(22)*(R0CMAX-R0CMIN)
C
      IF (IMTALK.EQ.1) THEN
C
          WRITE(IRESUL,'(''ARANDO(22)='',F20.16,'' PARAMT(22)='',F20.14,
     *                                         '' (ORCENT)'')')
     *                     ARANDO(22),              PARAMT(22)
      END IF
C
      IF (IFTAKE(22).EQ.1) ORCENT=PARAMT(22)
C_______________________________________________________________________
C
      PARAMT(31)=XRCMIN+ARANDO(31)*(XRCMAX-XRCMIN)
C
      IF (IMTALK.EQ.1) THEN
C
          WRITE(IRESUL,'(''ARANDO(31)='',F20.16,'' PARAMT(31)='',F20.14,
     *                                         '' (XK_R0C)'')')
     *                     ARANDO(31),             PARAMT(31)
      END IF
C
      IF (IFTAKE(31).EQ.1) OK_R0C=PARAMT(31)
C
C=======================================================================
C
      PARAMT(25)=R0SMIN+ARANDO(25)*(R0SMAX-R0SMIN)
C
      IF (IMTALK.EQ.1) THEN
C
          WRITE(IRESUL,'(''ARANDO(25)='',F20.16,'' PARAMT(25)='',F20.14,
     *                                         '' (ORSORB)'')')
     *                     ARANDO(25),              PARAMT(25)
      END IF
C
      IF (IFTAKE(25).EQ.1) ORSORB=PARAMT(25)
C_______________________________________________________________________
C
      PARAMT(34)=XRSMIN+ARANDO(34)*(XRSMAX-XRSMIN)
C
      IF (IMTALK.EQ.1) THEN
C
          WRITE(IRESUL,'(''ARANDO(34)='',F20.16,'' PARAMT(34)='',F20.14,
     *                                         '' (XK_RSO)'')')
     *                     ARANDO(34),             PARAMT(34)
      END IF
C
      IF (IFTAKE(34).EQ.1) OK_RSO=PARAMT(34)
C
C=======================================================================
C
      PARAMT(28)=REFMIN+ARANDO(28)*(REFMAX-REFMIN)
C
      IF (IRMFST.EQ.1) THEN
          PARAMT(28)=PARAMT(25)
      END IF
C
      IF (IMTALK.EQ.1) THEN
C
          WRITE(IRESUL,'(''ARANDO(28)='',F20.16,'' PARAMT(28)='',F20.14,
     *                                         '' (OREFFM)'')')
     *                     ARANDO(28),              PARAMT(28)
      END IF
C
      IF (IFTAKE(28).EQ.1) OREFFM=PARAMT(28)
C_______________________________________________________________________
C
      PARAMT(37)=XREMIN+ARANDO(37)*(XREMAX-XREMIN)
C
      IF (IRMFST.EQ.1) THEN
          PARAMT(37)=PARAMT(34)
      END IF
C
      IF (IMTALK.EQ.1) THEN
C
          WRITE(IRESUL,'(''ARANDO(37)='',F20.16,'' PARAMT(37)='',F20.14,
     *                                         '' (XK_REF)'')')
     *                     ARANDO(37),             PARAMT(37)
      END IF
C
      IF (IFTAKE(37).EQ.1) OK_REF=PARAMT(37)
C
C=======================================================================
C     
C=======================================================================
C
      PARAMT(23)=A0CMIN+ARANDO(23)*(A0CMAX-A0CMIN)
C
      IF (IMTALK.EQ.1) THEN
C
          WRITE(IRESUL,'(''ARANDO(23)='',F20.16,'' PARAMT(23)='',F20.14,
     *                                         '' (OACENT)'')')
     *                     ARANDO(23),              PARAMT(23)
      END IF
C
      IF (IFTAKE(23).EQ.1) OACENT=PARAMT(23)
C_______________________________________________________________________
C
      PARAMT(32)=XACMIN+ARANDO(32)*(XACMAX-XACMIN)
C 
      IF (IMTALK.EQ.1) THEN
C
          WRITE(IRESUL,'(''ARANDO(32)='',F20.16,'' PARAMT(32)='',F20.14,
     *                                         '' (XK_A0C)'')')
     *                     ARANDO(32),             PARAMT(32)
      END IF
C
      IF (IFTAKE(32).EQ.1) OK_A0C=PARAMT(32)
C
C=======================================================================
C
      PARAMT(26)=A0SMIN+ARANDO(26)*(A0SMAX-A0SMIN)
C
      IF (IMTALK.EQ.1) THEN
C
          WRITE(IRESUL,'(''ARANDO(26)='',F20.16,'' PARAMT(26)='',F20.14,
     *                                         '' (OASORB)'')')
     *                     ARANDO(26),              PARAMT(26)
      END IF
C
      IF (IFTAKE(26).EQ.1) OASORB=PARAMT(26)
C_______________________________________________________________________
C
      PARAMT(35)=XASMIN+ARANDO(35)*(XASMAX-XASMIN)
C
      IF (IMTALK.EQ.1) THEN
C
          WRITE(IRESUL,'(''ARANDO(35)='',F20.16,'' PARAMT(35)='',F20.14,
     *                                         '' (XK_ASO)'')')
     *                     ARANDO(35),             PARAMT(35)
      END IF
C
      IF (IFTAKE(35).EQ.1) OK_ASO=PARAMT(35)
C
C=======================================================================
C
      PARAMT(29)=AEFMIN+ARANDO(29)*(AEFMAX-AEFMIN)
C
      IF (IRMFST.EQ.1) THEN
          PARAMT(29)=PARAMT(26)
      END IF
C
      IF (IMTALK.EQ.1) THEN
C
          WRITE(IRESUL,'(''ARANDO(29)='',F20.16,'' PARAMT(29)='',F20.14,
     *                                         '' (OAEFFM)'')')
     *                     ARANDO(29),              PARAMT(29)
      END IF
C
      IF (IFTAKE(29).EQ.1) OAEFFM=PARAMT(29)
C_______________________________________________________________________
C
      PARAMT(38)=XAEMIN+ARANDO(38)*(XAEMAX-XAEMIN)
C
      IF (IRMFST.EQ.1) THEN
          PARAMT(38)=PARAMT(35)
      END IF
C
      IF (IMTALK.EQ.1) THEN
C
          WRITE(IRESUL,'(''ARANDO(38)='',F20.16,'' PARAMT(38)='',F20.14,
     *                                         '' (XK_LEF)'')')
     *                     ARANDO(38),             PARAMT(38)
      END IF
C
      IF (IFTAKE(38).EQ.1) OK_XEF=PARAMT(38)
C
C=======================================================================
C         Selecting the lambda values for density-dependent SO
C ###     and alpha_so and beta_so for the spin-current part
C=======================================================================
C
      PARAMT(39)=XPPMIN+ARANDO(39)*(XPPMAX-XPPMIN)
C
      IF (IMTALK.EQ.1) THEN
C
          WRITE(IRESUL,'(''ARANDO(39)='',F20.16,'' PARAMT(39)='',F20.14,
     *                                         '' (OLDXPP)'')')
     *                     ARANDO(39),             PARAMT(39)
      END IF
C
      IF (IFTAKE(39).EQ.1) OLDXPP=PARAMT(39)
C_______________________________________________________________________
C
      PARAMT(40)=XPNMIN+ARANDO(40)*(XPNMAX-XPNMIN)
C
      IF (IMTALK.EQ.1) THEN
C
          WRITE(IRESUL,'(''ARANDO(40)='',F20.16,'' PARAMT(40)='',F20.14,
     *                                         '' (OLDXPN)'')')
     *                     ARANDO(40),             PARAMT(40)
      END IF
C
      IF (IFTAKE(40).EQ.1) OLDXPN=PARAMT(40)
C
C=======================================================================
C
      PARAMT(41)=XNPMIN+ARANDO(41)*(XNPMAX-XNPMIN)
C
      IF (IMTALK.EQ.1) THEN
C
          WRITE(IRESUL,'(''ARANDO(41)='',F20.16,'' PARAMT(41)='',F20.14,
     *                                         '' (OLDXNP)'')')
     *                     ARANDO(41),             PARAMT(41)
      END IF
C
      IF (IFTAKE(41).EQ.1) OLDXNP=PARAMT(41)
C_______________________________________________________________________
C
      PARAMT(42)=XNNMIN+ARANDO(42)*(XNNMAX-XNNMIN)
C
      IF (IMTALK.EQ.1) THEN
C
          WRITE(IRESUL,'(''ARANDO(42)='',F20.16,'' PARAMT(42)='',F20.14,
     *                                         '' (OLDXNN)'')')
     *                     ARANDO(42),             PARAMT(42)
      END IF
C
      IF (IFTAKE(42).EQ.1) OLDXNN=PARAMT(42)
C
C=======================================================================
C     Selecting the lambda values for tensor part
C=======================================================================
C
      PARAMT(43)=YPPMIN+ARANDO(43)*(YPPMAX-YPPMIN)
C
      IF (IMTALK.EQ.1) THEN
C
          WRITE(IRESUL,'(''ARANDO(43)='',F20.16,'' PARAMT(43)='',F20.14,
     *                                         '' (OLDYPP)'')')
     *                     ARANDO(43),             PARAMT(43)
      END IF
C
      IF (IFTAKE(43).EQ.1) OLDYPP=PARAMT(43)
C_______________________________________________________________________
C
      PARAMT(44)=YPNMIN+ARANDO(44)*(YPNMAX-YPNMIN)
C
      IF (IMTALK.EQ.1) THEN
C
          WRITE(IRESUL,'(''ARANDO(44)='',F20.16,'' PARAMT(44)='',F20.14,
     *                                         '' (OLDYPN)'')')
     *                     ARANDO(44),             PARAMT(44)
      END IF
C
      IF (IFTAKE(44).EQ.1) OLDYPN=PARAMT(44)
C
C=======================================================================
C
      PARAMT(45)=YNPMIN+ARANDO(45)*(YNPMAX-YNPMIN)
C
      IF (IMTALK.EQ.1) THEN
C
          WRITE(IRESUL,'(''ARANDO(45)='',F20.16,'' PARAMT(45)='',F20.14,
     *                                         '' (OLDYNP)'')')
     *                     ARANDO(45),             PARAMT(45)
      END IF
C
      IF (IFTAKE(45).EQ.1) OLDYNP=PARAMT(45)
C_______________________________________________________________________
C
      PARAMT(46)=YNNMIN+ARANDO(46)*(YNNMAX-YNNMIN)
C
      IF (IMTALK.EQ.1) THEN
C
          WRITE(IRESUL,'(''ARANDO(46)='',F20.16,'' PARAMT(46)='',F20.14,
     *                                         '' (OLDYNN)'')')
     *                     ARANDO(46),             PARAMT(46)
      END IF
C
      IF (IFTAKE(46).EQ.1) OLDYNN=PARAMT(46)
C
C=======================================================================
C     Selecting the lambda values for tensor part
C=======================================================================
C
      PARAMT(47)=CPPMIN+ARANDO(47)*(CPPMAX-CPPMIN)
C
      IF (IMTALK.EQ.1) THEN
C
          WRITE(IRESUL,'(''ARANDO(47)='',F20.16,'' PARAMT(47)='',F20.14,
     *                                         '' (OLDCPP)'')')
     *                     ARANDO(47),             PARAMT(47)
      END IF
C
      IF (IFTAKE(47).EQ.1) OLDCPP=PARAMT(47)
C_______________________________________________________________________
C
      PARAMT(48)=CPNMIN+ARANDO(48)*(CPNMAX-CPNMIN)
C
      IF (IMTALK.EQ.1) THEN
C
          WRITE(IRESUL,'(''ARANDO(48)='',F20.16,'' PARAMT(48)='',F20.14,
     *                                         '' (OLDCPN)'')')
     *                     ARANDO(48),             PARAMT(48)
      END IF
C
      IF (IFTAKE(48).EQ.1) OLDCPN=PARAMT(48)
C
C=======================================================================
C
      PARAMT(49)=CNPMIN+ARANDO(49)*(CNPMAX-CNPMIN)
C
      IF (IMTALK.EQ.1) THEN
C
          WRITE(IRESUL,'(''ARANDO(49)='',F20.16,'' PARAMT(49)='',F20.14,
     *                                         '' (OLDCNP)'')')
     *                     ARANDO(49),             PARAMT(49)
      END IF
C
      IF (IFTAKE(49).EQ.1) OLDCNP=PARAMT(49)
C_______________________________________________________________________
C
      PARAMT(50)=CNNMIN+ARANDO(50)*(CNNMAX-CNNMIN)
C
      IF (IMTALK.EQ.1) THEN
C
          WRITE(IRESUL,'(''ARANDO(50)='',F20.16,'' PARAMT(50)='',F20.14,
     *                                         '' (OLDCNN)'')')
     *                     ARANDO(50),             PARAMT(50)
      END IF
C
      IF (IFTAKE(50).EQ.1) OLDCNN=PARAMT(50)
C
C=======================================================================
C=======================================================================
C
C     Real kappa parametrization
C
      PARAMT(51)=V0CMIN_KAPPAR+ARANDO(51)*(V0CMAX_KAPPAR-V0CMIN_KAPPAR)
C
      IF (IMTALK.EQ.1) THEN
C
          WRITE(IRESUL,'(''ARANDO(51)='',F20.16,'' PARAMT(51)='',F20.14,
     *                                         '' (OLDCNN)'')')
     *                     ARANDO(51),             PARAMT(51)
      END IF
C
      IF (IFTAKE(51).EQ.1) OLDV0C_KAPPAR=PARAMT(51)
C_______________________________________________________________________
C
      PARAMT(52)=XKCMIN_KAPPAR+ARANDO(52)*(XKCMAX_KAPPAR-XKCMIN_KAPPAR)
C
      IF (IMTALK.EQ.1) THEN
C
          WRITE(IRESUL,'(''ARANDO(52)='',F20.16,'' PARAMT(52)='',F20.14,
     *                                         '' (OLDCNN)'')')
     *                     ARANDO(52),             PARAMT(52)
      END IF
C
      IF (IFTAKE(52).EQ.1) OLDKVC_KAPPAR=PARAMT(52)
C_______________________________________________________________________
C
      PARAMT(53)=R0CMIN_KAPPAR+ARANDO(53)*(R0CMAX_KAPPAR-R0CMIN_KAPPAR)
C
      IF (IMTALK.EQ.1) THEN
C
          WRITE(IRESUL,'(''ARANDO(53)='',F20.16,'' PARAMT(53)='',F20.14,
     *                                         '' (OLDCNN)'')')
     *                     ARANDO(53),             PARAMT(53)
      END IF
C
      IF (IFTAKE(53).EQ.1) OLDR0C_KAPPAR=PARAMT(53)
C_______________________________________________________________________
C
      PARAMT(54)=XRCMIN_KAPPAR+ARANDO(54)*(XRCMAX_KAPPAR-XRCMIN_KAPPAR)
C
      IF (IMTALK.EQ.1) THEN
C
          WRITE(IRESUL,'(''ARANDO(54)='',F20.16,'' PARAMT(54)='',F20.14,
     *                                         '' (OLDCNN)'')')
     *                     ARANDO(54),             PARAMT(54)
      END IF
C
      IF (IFTAKE(54).EQ.1) OLDXRC_KAPPAR=PARAMT(54)
C_______________________________________________________________________
C
      PARAMT(55)=A0CMIN_KAPPAR+ARANDO(55)*(A0CMAX_KAPPAR-A0CMIN_KAPPAR)
C
      IF (IMTALK.EQ.1) THEN
C
          WRITE(IRESUL,'(''ARANDO(55)='',F20.16,'' PARAMT(55)='',F20.14,
     *                                         '' (OLDCNN)'')')
     *                     ARANDO(55),             PARAMT(55)
      END IF
C
      IF (IFTAKE(55).EQ.1) OLDA0C_KAPPAR=PARAMT(55)
C_______________________________________________________________________
C
      PARAMT(56)=XACMIN_KAPPAR+ARANDO(56)*(XACMAX_KAPPAR-XACMIN_KAPPAR)
C
      IF (IMTALK.EQ.1) THEN
C
          WRITE(IRESUL,'(''ARANDO(56)='',F20.16,'' PARAMT(56)='',F20.14,
     *                                         '' (OLDCNN)'')')
     *                     ARANDO(56),             PARAMT(56)
      END IF
C
      IF (IFTAKE(56).EQ.1) OLDXAC_KAPPAR=PARAMT(56)
C_______________________________________________________________________
C
      PARAMT(57)=XL_MIN_KAPPAR+ARANDO(57)*(XL_MAX_KAPPAR-XL_MIN_KAPPAR)
C
      IF (IMTALK.EQ.1) THEN
C
          WRITE(IRESUL,'(''ARANDO(57)='',F20.16,'' PARAMT(57)='',F20.14,
     *                                         '' (OLDCNN)'')')
     *                     ARANDO(57),             PARAMT(57)
      END IF
C
      IF (IFTAKE(57).EQ.1) OLD_XL_KAPPAR=PARAMT(57)
C_______________________________________________________________________
C
      PARAMT(58)=XKSMIN_KAPPAR+ARANDO(58)*(XKSMAX_KAPPAR-XKSMIN_KAPPAR)
C
      IF (IMTALK.EQ.1) THEN
C
          WRITE(IRESUL,'(''ARANDO(58)='',F20.16,'' PARAMT(58)='',F20.14,
     *                                         '' (OLDCNN)'')')
     *                     ARANDO(58),             PARAMT(58)
      END IF
C
      IF (IFTAKE(58).EQ.1) OLDXKS_KAPPAR=PARAMT(58)
C_______________________________________________________________________
C
      PARAMT(59)=R0SMIN_KAPPAR+ARANDO(59)*(R0SMAX_KAPPAR-R0SMIN_KAPPAR)
C
      IF (IMTALK.EQ.1) THEN
C
          WRITE(IRESUL,'(''ARANDO(59)='',F20.16,'' PARAMT(59)='',F20.14,
     *                                         '' (OLDCNN)'')')
     *                     ARANDO(59),             PARAMT(59)
      END IF
C
      IF (IFTAKE(59).EQ.1) OLDR0S_KAPPAR=PARAMT(59)
C_______________________________________________________________________
C
      PARAMT(60)=XRSMIN_KAPPAR+ARANDO(60)*(XRSMAX_KAPPAR-XRSMIN_KAPPAR)
C
      IF (IMTALK.EQ.1) THEN
C
          WRITE(IRESUL,'(''ARANDO(60)='',F20.16,'' PARAMT(60)='',F20.14,
     *                                         '' (OLDCNN)'')')
     *                     ARANDO(60),             PARAMT(60)
      END IF
C
      IF (IFTAKE(60).EQ.1) OLDXRS_KAPPAR=PARAMT(60)
C_______________________________________________________________________
C
      PARAMT(61)=A0SMIN_KAPPAR+ARANDO(61)*(A0SMAX_KAPPAR-A0SMIN_KAPPAR)
C
      IF (IMTALK.EQ.1) THEN
C
          WRITE(IRESUL,'(''ARANDO(61)='',F20.16,'' PARAMT(61)='',F20.14,
     *                                         '' (OLDCNN)'')')
     *                     ARANDO(61),             PARAMT(61)
      END IF
C
      IF (IFTAKE(61).EQ.1) OLDA0S_KAPPAR=PARAMT(61)
C_______________________________________________________________________
C
      PARAMT(62)=XASMIN_KAPPAR+ARANDO(62)*(XASMAX_KAPPAR-XASMIN_KAPPAR)
C
      IF (IMTALK.EQ.1) THEN
C
          WRITE(IRESUL,'(''ARANDO(62)='',F20.16,'' PARAMT(62)='',F20.14,
     *                                         '' (OLDCNN)'')')
     *                     ARANDO(62),             PARAMT(62)
      END IF
C
      IF (IFTAKE(62).EQ.1) OLDXAS_KAPPAR=PARAMT(62)
C
C=======================================================================
C
      IF (LOGWRI.GT.0) THEN
          IF (ISOSPI.EQ.0) THEN
C      
              WRITE(NOUTPT,'(/,80(''#''),/,''#'',T80,''#'',/,
     *            ''#                     Starting Parameter '',
     *            ''Values (Neutrons)                     #'',/,
     *            ''#'',T80,''#'',/,
     *            ''#    Central Potential         Spin-Orbit '',
     *            ''Potential        Effective Mass      #'',/,
     *            ''#'',T80,''#'',/,
     *            ''#    Vo      Ro      a       Lambda    Ro      '',
     *            ''a       Lambda     Ro     a     #'',/,
     *            ''# '',F7.3,2F8.4,3X,F7.3,2F8.4,3X,F7.3,2F8.4,''  #'',
     *            /,''#'',T80,''#'',/,80(''#''))')
     *
     *            PARAMT(1),PARAMT(2),PARAMT(3),
     *            PARAMT(4),PARAMT(5),PARAMT(6),
     *            PARAMT(7),PARAMT(8),PARAMT(9) 
C
          ELSE
C
              WRITE(NOUTPT,'(/,80(''#''),/,''#'',T80,''#'',/,
     *            ''#                     Starting Parameter '',
     *            ''Values (Protons)                      #'',/,
     *            ''#'',T80,''#'',/,
     *            ''#    Central Potential         Spin-Orbit '',
     *            ''Potential        Effective Mass      #'',/,
     *            ''#'',T80,''#'',/,
     *            ''#    Vo      Ro      a       Lambda    Ro      '',
     *            ''a       Lambda     Ro     a     #'',/,
     *            ''# '',F7.3,2F8.4,3X,F7.3,2F8.4,3X,F7.3,2F8.4,''  #'',
     *            /,''#'',T80,''#'',/,80(''#''))')
     *
     *            PARAMT(1),PARAMT(2),PARAMT(3),
     *            PARAMT(4),PARAMT(5),PARAMT(6),
     *            PARAMT(7),PARAMT(8),PARAMT(9) 
C
          END IF
C
      END IF      
C
C=======================================================================
C
      DO I=1,NDPARS
         STARTV(I)=PARAMT(I)
      END DO
C
C=======================================================================
C     Filtering the parameters which will be varied by the minimisation
C                                                           subroutine
C=======================================================================
C
      LDSEAR=0
C
C=======================================================================
C     Beginning with the central potential depth...
C=======================================================================
C
      DO I=1,NDPARS
         IF (IFTAKE(I).EQ.1) THEN
             LDSEAR=LDSEAR+1
             ARGPAR(LDSEAR)=PARAMT(I)
         END IF
      END DO
C  
C     In principle, the parameters of the spin-orbit term
C     and those related to the effective mass term should
C     be equal. On the phenomenological basis,  it may be 
C     of interest to let varying the three effective mass
C     parameters independently. For that reason we select
C     the control IRMFST (reduced mass fixed)  equal to 0
C
C     IF (IRMFST.EQ.0) THEN
C          
C         IF (LAMEFF.EQ.'YES') THEN
C             LDSEAR=LDSEAR+1             
C             ARGPAR(LDSEAR)=PARAMT(3)
C         END IF
C
C     END IF
C
C     IF (IRMFST.EQ.0) THEN
C          
C         See comments just above              
C          
C         IF (RADEFF.EQ.'YES') THEN
C             LDSEAR=LDSEAR+1             
C             ARGPAR(LDSEAR)=PARAMT(6)
C         END IF
C
C     END IF
C
C     IF (IRMFST.EQ.0) THEN
C          
C         See comments just above              
C          
C         IF (DIFEFF.EQ.'YES') THEN
C             LDSEAR=LDSEAR+1             
C             ARGPAR(LDSEAR)=PARAMT(9)
C         END IF
C
C     END IF
C
C=======================================================================
C
      LDSEAR=0
C
      DO I=1,NDPARS
C
         IF ((IFTAKE(I).EQ.1).AND.(IFIRST.NE.1)) THEN
C
             PARAMT(I)=VMIMIN(I)+ARANDO(I)*(VMIMAX(I)-VMIMIN(I))
             LDSEAR=LDSEAR+1
             ARGPAR(LDSEAR)=PARAMT(I)
C
         ELSE
             PARAMT(I)=VMISTR(I)
         END IF
C
      END DO
C
C=======================================================================
C     Reinitializing the random number generator to a new value
C=======================================================================
C
      PIOVE2=2.0*ATAN(1.0)
C
      NEWSED=NINT(10000.*ABS(COS(PIOVE2*ARANDO(1))))
C
C=======================================================================
C
      IF (LOGWRI.GT.0) THEN
          WRITE(LOGFIL,'(9X,''Exiting  NEWRAN'')')
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
      SUBROUTINE FUNMIN(NDFUNC,LDFUNC,NDPARS,LDPARS,ARGPAR,FUNVEC,
     *                                              FUNJAC,I_FLAG)
      INCLUDE   'MATDIM/NDSPEC.f'
      INCLUDE   'MATDIM/NDNUCL.f'
      INCLUDE   'MATDIM/NDLEXP.f'
      INCLUDE   'MATDIM/MAXFEV.f'
      INCLUDE   'MATDIM/NDIM_M.f'
      INCLUDE   'MATDIM/NDIM_P.f'
      INCLUDE   'MATDIM/ND_RHO.f'
      INCLUDE   'MATDIM/N_NOYX.f'
      INCLUDE   'MATDIM/NDTITL.f'
C_______________________________________________________________________
C      
      CHARACTER
     *          LABTHE_PROTON*6,LABTHE_NEUTRS*6
      CHARACTER
     *          LABEXP_PROTON*6,LABEXP_NEUTRS*6
      CHARACTER
     *          TITLES_EXPORT*12
      CHARACTER
     *          ACTION*1,TYPCHI*6
      CHARACTER
     *          BLANCK*1
      CHARACTER
     *          INPSYM*6,NUCSYM*6,SYMBNU*5
C_______________________________________________________________________
C      
      DIMENSION
     *          FUNCT0(1:NDIM_M),FUNCT1(1:NDIM_M),FUNEXP(1:NDIM_M),
     *          LEVELS(1:NDIM_M),LEXPLS(1:NDIM_M),DERIVF(1:NDIM_M)
      DIMENSION 
     *          ARGPAR(1:NDPARS),FUNVEC(1:NDFUNC)
      DIMENSION
     *          WEIGHT_FUNVEC(1:NDIM_M)
      DIMENSION 
     *          FUNJAC(1:NDFUNC,1:NDPARS)
      DIMENSION
     *          DEGENE(1:NDIM_M)
      DIMENSION
     *          WEISQR_FUNVEC(1:NDIM_M)
      DIMENSION
     *          CHISQU_AUXILP(1:NDPARS,1:NDNUCL),
     *          CHISQU_AUXILN(1:NDPARS,1:NDNUCL),
     *          CHISQU_GRADNT(1:NDPARS)
      DIMENSION
     *          FUNAUX(1:NDIM_M),
     *          IVECTR(1:NDIM_M)
      DIMENSION
     *          CHIWEI_NUCLEU(1:NDNUCL),
     *          CHIWEI_NUCPRO(1:NDNUCL),
     *          CHIWEI_NUCNEU(1:NDNUCL)
      DIMENSION        
     *          CHIRAD_NUCLEU(1:NDNUCL),
     *          CHIRAD_NUCPRO(1:NDNUCL),
     *          CHIRAD_NUCNEU(1:NDNUCL)
      DIMENSION
     *          DEGSUM_NUCPRO(1:NDNUCL),
     *          DEGSUM_NUCNEU(1:NDNUCL)
      DIMENSION
     *          IACTIV_VECTOR(1:NDPARS)
C_______________________________________________________________________
C      
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
     *       /NOWEIG/ CHIENE_PROTON(1:NDNUCL),CHIENE_NEUTRS(1:NDNUCL),
     *                CHIENE_TOTALS(1:NDNUCL),
     *                CHIRAD_PROTON(1:NDNUCL),CHIRAD_NEUTRS(1:NDNUCL),
     *                CHIRAD_TOTALS(1:NDNUCL)
      COMMON
     *       /NUCWEI/ WEINUC_PROTON(1:NDNUCL),
     *                WEINUC_NEUTRS(1:NDNUCL),
     *                WEIRAD_PROTON(1:NDNUCL),
     *                WEIRAD_NEUTRS(1:NDNUCL)
      COMMON
     *       /TITEXP/ TITLES_EXPORT(1:NDTITL)
      COMMON
     *       /WEIGHT/ WEIGHT_CORREL,WEIGHT_RADIUS,WEIGHT_INVERT,
     *                WEIGHT_EFERMI,WEIGHT_ENEGAP,WEIWEI,
     *                WEIGHT_ERRABS,WEIGHT_EABSAV,WEIGHT_ERRMAX,
     *                WEIGHT_DENSUP,WEIGHT_DENSDW,WEIGHT_RHODEN
      COMMON
     *       /DIMSPH/ LDSING_PROTON,LEVTHE_PROTON,
     *                LEVEXP_PROTON(1:NDNUCL),
     *                LDSING_NEUTRS,LEVTHE_NEUTRS,
     *                LEVEXP_NEUTRS(1:NDNUCL)
      COMMON
     *       /MODIFR/ RMSEXP_PROTON(1:NDNUCL),
     *                RMSEXP_NEUTRS(1:NDNUCL),
     *                RMSTHE_PROTON,RMSTHE_NEUTRS,
     *                RMODIF_PROTON,RMODIF_NEUTRS,
     *                              PUSHIN,DELTAR
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
     *       /SEARCH/ IFTAKE(1:NDIM_P)
      COMMON
     *       /CNTROL/ PARAMT(1:NDIM_P),
     *                DPARAM(1:NDIM_P)  
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
     *       /MESHIN/ I_MESH(1:NDIM_P),
     *                XMIN_I(1:NDIM_P),
     *                XMAX_I(1:NDIM_P),
     *                MAXPAR(1:NDIM_P)
      COMMON
     *       /TESFIT/ IFTEST,IFFITS,IFMESH,ISPE_P,IFMOCA
      COMMON 
     *       /CHITAK/ IF_SPE,IF_RAD,IF_GAP,IF_FER,IF_DEN,IF_RHO,IF_INV,
     *                LDCHI2,TYPCHI(1:NDNUCL),ITKCHI(1:NDNUCL)   
      COMMON
     *       /ACTACT/ ACTION(1:NDIM_P)
      COMMON
     *       /POTPOT/ PARPOT(1:NDIM_P)
      COMMON
     *       /ARGINI/ XINITS(1:NDIM_P)
      COMMON
     *       /P_N_LM/ IFPROT,IFNEUT,IFBOTH
      COMMON
     *       /MINITE/ IDEFCN,ITECHI
      COMMON
     *       /CHI_PN/ CHISQU_PROTON,CHISQU_NEUTRS
      COMMON
     *       /CHICHO/ VALCHI
      COMMON
     *       /CHIGRD/ CHISQU_GRDNRM,
     *                CHISQU_GRADIE(1:NDIM_P)
      COMMON
     *       /DENSTR/ IFDENS
      COMMON
     *       /TENSOR/ IFTENS
      COMMON
     *       /TENOPT/ ICENTT,ISORBT,ITENSR
      COMMON
     *       /IFKAPP/ IFK_VC,IFK_VS,IFK_RC,IFK_RS,IFK_AC,IFK_AS
      COMMON
     *       /STOPIN/ HIAUXI(1:MAXFEV)
      COMMON
     *       /THEDEN/ N_NUCL,I_NUCL(1:N_NOYX),
     *                RHOTHE(1:ND_RHO,1:N_NOYX)
      COMMON
     *       /EXPDEN/ LABELN(1:N_NOYX,1:2),
     *                R_MESH(1:ND_RHO,1:N_NOYX),
     *                RHOEXP(1:ND_RHO,1:N_NOYX)
     *
     *       /EXPSYM/ SYMBNU(1:N_NOYX)
      COMMON
     *       /CNTCNT/ ICOUNT_DENGRD,ICOUNT_FUNMIN,ICOUNT_HAMMAT,
     *                ICOUNT_EXPTHE,ICOUNT_CHOICE,ICOUNT_DIAMAT
      COMMON
     *       /PRFILE/ ISCREN,LOGFIL,LOGCHP,LOGCHN,LOGHIP,LOGHIN,
     *                LOGPRO,LOGNEU,LOGWRI,LSCREN,IRESUL,ICONVE,
     *                NOUTPT,LOGAUX,LOGENE,LOGRAD,LOGGRD,LOGMSH,
     *                                                   LOGBIS            
C
C=======================================================================
C
C     This subroutine constructs function FUNVEC, as the difference
C     between the experimental and theoretical data sets the latter
C     composed of, e.g., single particle energies, radii, densities
C     etc. FUNVEC is then used by Levenberg - Marquard minimisation
C     subroutine LEVMAR.
C 
C     This subroutine also constructs the Jacobian matrix FUNJAC by
C     using the method of finite differences.
C
C     If I_FLAG = 1, calculating the functions at ARGMNT and return
C     the function-vector in FUNVEC. Here we do not alter Jacobian.
C
C     If I_FLAG = 2, calculating the Jacobian at ARGMNT, and return
C     the corresponding matrix in FUNJAC.  Here not altering FUNVEC.
C
C     If IFMESH = 1, we prepare a map by minimising the chi^2, over 
C                                       a square-set of mesh points.
C
C=======================================================================
C
C     How we have defined the \chi^2 with the corresponding weights
C     and normalization factors:
C
C     w_k(e) - energy weight for nucleus k
C     w_k(r) - radius weight for nucleus k
C
C     \chi^2_k = w_k(e) \sum_i [ e_i(th) - e_i(exp) ]^2 * (2j_i+1)
C              + w_k(r) [ r_k(th)-r_k(exp) ]^2
C
C     N = \sum_k [ w_k(e) \sum_i (2j_i+1) + w_k(r) ] 
C
C     The total \chi^2 is:
C                          \chi^2 = [\sum_k (\chi^2_k)] / N
C
C=======================================================================
C
      CALL CPUTIM('FUNMIN',1)
C
      ICOUNT_FUNMIN=ICOUNT_FUNMIN+1
C
C=======================================================================
C
      IF (NDIM_P.NE.NDPARS) THEN
          WRITE(LSCREN,'(/,''Alarm in FUNMIN: NDIM_P= '',I3,1X,
     *                     ''and NDPARS= '',I3,''are not equal!!'',/)')
     *                           NDIM_P,NDPARS
          STOP 'STOP in FUNMIN: NDIM_P.NE.NDPARS'
      END IF
C
C=======================================================================
C
      IF (IFPROT.EQ.1) ISOSPI_MINIMI=1
      IF (IFNEUT.EQ.1) ISOSPI_MINIMI=0
C
C=======================================================================
C
      IF (LOGWRI.GT.4) THEN
          WRITE(LOGFIL,'(12X,''Entering VERIFY from FUNMIN, '',
     *                       ''IDEFCN='',I3)')
     *                         IDEFCN           
      END IF
C  @@@ IRENE - WHAT IS ABORTED AND WHY? MEANING OF THIS VARIABLE?     
      CALL VERIFY(ARGPAR,IABORT)
C
      IF (LOGWRI.GT.4) THEN
          WRITE(LOGFIL,'(12X,''After VERIFY in FUNMIN, IABORT= '',I1)') 
     *                                                 IABORT
      END IF
C
C=======================================================================
C
C     Here we collect the active (varying) parameters which are stored
C     in PARPOT(1), PARPOT(2), … PARPOT(MAX) whereas ARGPAR correspond
C     to a fixed storage in which active and non-active parameters are
C                                                                mixed
      IACTIV=0 ! IACTIV will count, how many active parameters we have
C
      IF (LOGWRI.GT.0) THEN
          WRITE(LOGFIL,'(12X,''Parameters in FUNMIN'',
     *                       ''(IDEFCN,IPARAM,IACTIV,'',
     *                       ''ARGPAR(IACTIV),PARPOT(IACTIV))'')')
      END IF
C
      DO IPARAM=1,NDPARS  
C                  
         ACTION(IPARAM)=' '
C
         IF (IFTAKE(IPARAM).EQ.1) THEN
C
             IACTIV=IACTIV+1
             PARPOT(IPARAM)=ARGPAR(IACTIV)
C
             IACTIV_VECTOR(IPARAM)=IACTIV
C
             IF (LOGWRI.GT.0) THEN
                 WRITE(LOGFIL,'(T13,''PARAM ='',A,T31,3I7,2F15.4)') 
     *                          TITLES_EXPORT(IPARAM),
     *                          IDEFCN,IPARAM,IACTIV,ARGPAR(IACTIV),
     *                                               PARPOT(IPARAM)
             END IF
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
      IF (LOGWRI.GT.0) WRITE(LOGFIL,'()')
C
C=======================================================================
C
C     Below: Setting to zero the future total chi^2 functions which
C            will include the individual, physicist-defined weight- 
C                                                           factors
C
      DO INUCLI=1,LDNUCL
C          
         CHIWEI_NUCLEU(INUCLI)=0.0
         CHIWEI_NUCPRO(INUCLI)=0.0
         CHIWEI_NUCNEU(INUCLI)=0.0
C                                 Similar for the radius chi^2 only
         CHIRAD_NUCLEU(INUCLI)=0.0
         CHIRAD_NUCPRO(INUCLI)=0.0
         CHIRAD_NUCNEU(INUCLI)=0.0
C       
      END DO
C
C=======================================================================
C
      ILEVEL=0
      LEVPRO=0
      LEVNEU=0
C
      INITI1=0
      INITI2=0
C
      IAUXIL=0 ! Counting the number of differences (exp-th) we have
C                                 (used in the Jacobian construction)
C
C     The meaning of the variables below:
C
C     CHISQU_PROTON and CHISQU_NEUTRS =>> total chi^2 summed on all
C                                         the components and nuclei
C
      CHISQU_PROTON=0.0 ! As above, for protons
      CHISQU_NEUTRS=0.0 ! As above, for neutrons
      CHISQU_TOTALS=0.0 ! As above, protons and neutrons together
C
      
C
      CHIRAD_WEIPRO=0.0
      CHIRAD_WEINEU=0.0
      CHIRHO_WEIPRO=0.0
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
         DEGSUM_NUCPRO(INUCLI)=0.0
         DEGSUM_NUCNEU(INUCLI)=0.0
      END DO
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
                 DEGSUM_NUCPRO(INUCLI)=DEGSUM_PROTON
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
                 DEGSUM_NUCNEU(INUCLI)=DEGSUM_NEUTRS
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
C     Starting the do-loop over the active nuclei
C
      CHIENE_WEIPRO=0.0 ! As above but including the s.p. energies only
      CHIENE_WEINEU=0.0 ! As above but including the s.p. energies only         
C
      DO INUCLI=1,LDNUCL
C     
         IF (ITAKNU(INUCLI).EQ.1) THEN
C
             IZ_FIX=NUMB_Z(INUCLI)
             IN_FIX=NUMB_N(INUCLI)
C
             IF (LOGWRI.GT.4) THEN
                 WRITE(LOGFIL,'(12X,''Entering WS_RUN from FUNMIN '',
     *                              ''[1] with IZ_FIX='',I3, 
     *                                      '' IN_FIX='',I3,
     *                                      '' INUCLI='',I1)')
     *                                         IZ_FIX,IN_FIX,INUCLI
c                WRITE(LOGFIL,'()')
             END IF
C
                                                     I_MODE=0
             CALL WS_RUN(INUCLI,CHISQU_CHOICE,IABORT,I_MODE,I_FLAG,
     *                                 CHISQU_AUXIL1,CHISQU_AUXIL2)
C_______________________________________________________________________
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
C            The following quantities are calculated SEPARATELY for
C            each nucleus; the \chi^2 contributions contain the weights
C
             CHIWEI_NUCPRO(INUCLI)=CHIWEI_PROTON*WEINUC_PROTON(INUCLI)
             CHIWEI_NUCNEU(INUCLI)=CHIWEI_NEUTRS*WEINUC_NEUTRS(INUCLI)
             CHIWEI_NUCLEU(INUCLI)=CHIWEI_NUCPRO(INUCLI)
     *                            +CHIWEI_NUCNEU(INUCLI)
C        
             CHIRAD_NUCPRO(INUCLI)=RADDIF_PROTON*WEIRAD_PROTON(INUCLI)
             CHIRAD_NUCNEU(INUCLI)=RADDIF_NEUTRS*WEIRAD_NEUTRS(INUCLI)
             CHIRAD_NUCLEU(INUCLI)=CHIRAD_NUCPRO(INUCLI)
     *                            +CHIRAD_NUCNEU(INUCLI)
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
C_______________________________________________________________________
C
C            The following quantities are calculated SEPARATELY for
C            each nucleus; the \chi^2 contributions do not contain 
C            any weights
C            These quantities are then pinted in LOGENE and LOGRAD
C        
             CHIENE_PROTON(INUCLI)=SQRT(CHIDEG_PROTON
     *                            /     DEGSUM_NUCPRO(INUCLI))
             CHIENE_NEUTRS(INUCLI)=SQRT(CHIDEG_NEUTRS
     *                            /     DEGSUM_NUCNEU(INUCLI))
             CHIENE_TOTALS(INUCLI)=CHIENE_PROTON(INUCLI)
     *                            +CHIENE_NEUTRS(INUCLI)
C
             CHIRAD_PROTON(INUCLI)=RMSTHE_PROTON-RMSEXP_PROTON(INUCLI)
             CHIRAD_NEUTRS(INUCLI)=RMSTHE_NEUTRS-RMSEXP_NEUTRS(INUCLI)
             CHIRAD_TOTALS(INUCLI)=RMSTHE_PROTON-RMSEXP_PROTON(INUCLI)
     *                            +RMSTHE_NEUTRS-RMSEXP_NEUTRS(INUCLI)
C
C=======================================================================                 
C=======================================================================                 
C
C            Option 01: Density-dependent spin-orbit
C
             IF (IFDENS.EQ.1 .OR. (IFDENS.EQ.0 .AND. IFBOTH.EQ.1)) THEN
C
C                Item No. 01: The single particle energies
C         
                 IF (IF_SPE.EQ.1) THEN ! The single-particle energy
C                                                      contribution
C                    PROTONS
C
                     DO IEXPER=1,LEVEXP_PROTON(INUCLI)
                        DO ITHEOR=1,LEVTHE_PROTON
C
                           IF (LABEXP_PROTON(INUCLI,IEXPER).EQ.
     *                         LABTHE_PROTON(ITHEOR)) THEN
C
                               ILEVEL=ILEVEL+1
                               LEVPRO=LEVPRO+1
                               LEVPRO=ILEVEL
C
                               IF (IEXPER.EQ.1) INITI1=ILEVEL
C
                               FUNCT0(ILEVEL)=ENETHE_PROTON(ITHEOR)
                               FUNEXP(ILEVEL)=EXPEXP_PROTON(INUCLI,
     *                                                      IEXPER)
                               LEVELS(ILEVEL)=ITHEOR
                               LEXPLS(ILEVEL)=IEXPER
C
                               DEGENE(ILEVEL)
     *                        =
     *                         REAL(IDEGEX_PROTON(INUCLI,IEXPER))
C
                               WEISQR_FUNVEC(ILEVEL)
     *                        =
     *                         SQRT(DEGENE(ILEVEL)*WEINUC_PROTON(INUCLI)
     *                                            /SUMWEI_PROTON)           
C
                           END IF
C
                        END DO ! ITHEOR
                     END DO ! IEXPER
C_______________________________________________________________________
C
C                    NEUTRONS
C
                     DO IEXPER=1,LEVEXP_NEUTRS(INUCLI)
                        DO ITHEOR=1,LEVTHE_NEUTRS
                   
                           IF (LABEXP_NEUTRS(INUCLI,IEXPER).EQ.
     *                         LABTHE_NEUTRS(ITHEOR))      THEN
C
                               ILEVEL=ILEVEL+1
                               LEVNEU=LEVNEU+1
                               LEVNEU=ILEVEL
C
                               IF (IEXPER.EQ.1) INITI2=ILEVEL
C
                               FUNCT0(ILEVEL)=ENETHE_NEUTRS(ITHEOR)
                               FUNEXP(ILEVEL)=EXPEXP_NEUTRS(INUCLI,
     *                                                      IEXPER)
                               LEVELS(ILEVEL)=ITHEOR
                               LEXPLS(ILEVEL)=IEXPER
C
                               DEGENE(ILEVEL)
     *                        =
     *                         REAL(IDEGEX_NEUTRS(INUCLI,IEXPER))
C
                               WEISQR_FUNVEC(ILEVEL)
     *                        =
     *                         SQRT(DEGENE(ILEVEL)*WEINUC_NEUTRS(INUCLI)
     *                                            /SUMWEI_NEUTRS)
                           END IF
C
                        END DO
                     END DO
C         
                 END IF ! IF_SPE=1
C         
C            IF (IFDENS.EQ.1) THEN  ==>>>  ==>>>
C 
C=======================================================================
C
C                Item No. 02: The proton and neutron r.m.s. radii
C          
                 IF (IF_RAD.EQ.1) THEN
C
                     ILEVEL=ILEVEL+1
C                    INITI3=ILEVEL
C
                     FUNCT0(ILEVEL)=RMSTHE_PROTON
                     FUNEXP(ILEVEL)=RMSEXP_PROTON(INUCLI)
C
                     WEISQR_FUNVEC(ILEVEL)=SQRT(WEIRAD_PROTON(INUCLI)
     *                                    *     WEINUC_PROTON(INUCLI)
     *                                    /     SUMWEI_PROTON) 
C
                     ILEVEL=ILEVEL+1
C
                     FUNCT0(ILEVEL)=RMSTHE_NEUTRS
                     FUNEXP(ILEVEL)=RMSEXP_NEUTRS(INUCLI)
C
                     WEISQR_FUNVEC(ILEVEL)=SQRT(WEIRAD_NEUTRS(INUCLI)
     *                                    *     WEINUC_NEUTRS(INUCLI)
     *                                    /     SUMWEI_NEUTRS) 
C
                 END IF ! IF_RAD=1
C
C=======================================================================
C         
C            IF (IFDENS.EQ.1) THEN  ==>>>  ==>>>
C
C                Item No. 03: The proton and neutron main gaps
C          
                 IF (IF_GAP.EQ.1) THEN
C
                     ILEVEL=ILEVEL+1
C                    INITI4=ILEVEL
C
                     FUNCT0(ILEVEL)=GAPTHE_PROTON
                     FUNEXP(ILEVEL)=GAPEXP_PROTON(INUCLI)
C
                     WEISQR_FUNVEC(ILEVEL)=SQRT(WEIGHT_ENEGAP
     *                                    /     SUMWEI_PROTON) 
C
                     ILEVEL=ILEVEL+1
C
                     FUNCT0(ILEVEL)=GAPTHE_NEUTRS
                     FUNEXP(ILEVEL)=GAPEXP_NEUTRS(INUCLI)
C
                     WEISQR_FUNVEC(ILEVEL)=SQRT(WEIGHT_ENEGAP
     *                                    /     SUMWEI_NEUTRS) 
C
                 END IF ! IF_GAP=1
C         
C======================================================================= 
C         
C            IF (IFDENS.EQ.1) THEN ==>>> ==>>>
C
C                Item No. 04: The proton and neutron "Fermi energies"
C         
                 IF (IF_FER.EQ.1) THEN
             
                     ILEVEL=ILEVEL+1
C                    INITI5=ILEVEL
C
                     FUNCT0(ILEVEL)=FERTHE_PROTON
                     FUNEXP(ILEVEL)=FERMEX_PROTON(INUCLI)
C
                     WEISQR_FUNVEC(ILEVEL)=SQRT(WEIGHT_EFERMI
     *                                    /     SUMWEI_PROTON) 
C
                     ILEVEL=ILEVEL+1
C
                     FUNCT0(ILEVEL)=FERTHE_NEUTRS
                     FUNEXP(ILEVEL)=FERMEX_NEUTRS(INUCLI)
C
                     WEISQR_FUNVEC(ILEVEL)=SQRT(WEIGHT_EFERMI
     *                                    /     SUMWEI_NEUTRS) 
C             
                 END IF ! IF_FER=1
C         
C=======================================================================
C         
C            IF (IFDENS.EQ.1) THEN ==>>> ==>>>
C
C                Item No. 05: The proton and neutron densities defined
C                                        separately below-, and above
                 IF (IF_DEN.EQ.1) THEN !                 the main gaps
             
                     ILEVEL=ILEVEL+1
C                    INITI6=ILEVEL
C     PROTONS
C
                     FUNCT0(ILEVEL)=DENLOW_PROTON
                     FUNEXP(ILEVEL)=DENSDW_PROTON(INUCLI)
C
                     WEISQR_FUNVEC(ILEVEL)=SQRT(WEIGHT_DENSDW
     *                                    /     SUMWEI_PROTON) 
C_______________________________________________________________________
C             
                     ILEVEL=ILEVEL+1
C
                     FUNCT0(ILEVEL)=DENUPP_PROTON
                     FUNEXP(ILEVEL)=DENSUP_PROTON(INUCLI)
C
                     WEISQR_FUNVEC(ILEVEL)=SQRT(WEIGHT_DENSUP
     *                                    /     SUMWEI_PROTON) 
C_______________________________________________________________________
C
C     NEUTRONS
C
                     ILEVEL=ILEVEL+1
C
                     FUNCT0(ILEVEL)=DENLOW_NEUTRS
                     FUNEXP(ILEVEL)=DENSDW_NEUTRS(INUCLI)
C
                     WEISQR_FUNVEC(ILEVEL)=SQRT(WEIGHT_DENSDW
     *                                    /     SUMWEI_NEUTRS) 
C_______________________________________________________________________
C 
                     ILEVEL=ILEVEL+1
C
                     FUNCT0(ILEVEL)=DENUPP_NEUTRS
                     FUNEXP(ILEVEL)=DENSUP_NEUTRS(INUCLI)
C
                     WEISQR_FUNVEC(ILEVEL)=SQRT(WEIGHT_DENSUP
     *                                    /     SUMWEI_NEUTRS) 
C                  
                 END IF ! IF_DEN=1
C
C=======================================================================
C         
C            IF (IFDENS.EQ.1) THEN ==>>> ==>>>
C
C                Item No. 06: The proton charge densities as functions
C                                               of the radial variable
                 IF (IF_RHO.EQ.1) THEN
*C            
*                     CALL RHOTHE_DENSIT(INUCLI,IZ_FIX,IN_FIX)
*C
*                     DO K_NOYX=1,N_NUCL
*C
*                        I_NOYX=I_NUCL(K_NOYX)
*C
*                        DO ID_RHO=1,ND_RHO
*C
*                           ILEVEL=ILEVEL+1
                   
*                           FUNCT0(ILEVEL)=RHOTHE(ID_RHO,I_NOYX)
                   
*                           FUNEXP(ILEVEL)=RHOEXP(ID_RHO,I_NOYX)
                   
*                           WEISQR_FUNVEC(ILEVEL)=SQRT(WEIGHT_RHODEN)
*C
*                        END DO
*                     END DO
C         
                 END IF ! IF_RHO=1
C
C======================================================================= 
C         
C            IF (IFDENS.EQ.1) THEN ==>>> ==>>>
C        
C                Item No. 07: The proton level order as criterion
C                                        Not tested, Neutrons missing
C                IF (IF_INV.EQ.1) THEN
C             
C                    LEVTOT=LEVPRO+LEVNEU
C             
C                    DO I=1,LEVPRO
C               
C                       IF (I.LT.LEVPRO) THEN
C               
C                           ILWEXP=LEXPLS(I)       !Experimental
C                           IUPEXP=LEXPLS(I+1)
C                           DIFEXP=EXPEXP_PROTON(INUCLI,IUPEXP)
C     *                           -EXPEXP_PROTON(ILWEXP)
C
C                           ILWTHE=LEVELS(I)       !Theoretical
C                           IUPTHE=LEVELS(I+1)               
C                           DIFTHE=ENETHE_PROTON(INUCLI,IUPTHE)
C     *                           -ENETHE_PROTON(ILWTHE)
C             
C                           IF (((DIFEXP.LT.0.).AND.(DIFTHE.GT.0.)).OR.
C     *                         ((DIFEXP.GT.0.).AND.(DIFTHE.LT.0.))) 
C                                                                 THEN  
C                               … incorrect order of the levels
C             
C                               ILEVEL=ILEVEL+1
C                  
C                               AUXIL1=ENETHE_PROTON(ILWTHE)*DEGENE(I)
C                               AUXIL2=ENETHE_PROTON(IUPTHE)*DEGENE(I+1)
C                 
C                               FUNVEC(ILEVEL)
C     *                        =
C     *                        (AUXIL2-AUXIL1)*WEIGHT_INVERT
C                  
C                           END IF
C
C                       END IF ! LEVPRO
C             
C                    END DO
C         
C                END IF ! IF_INV=1
C         
C=======================================================================
C=======================================================================                 
C         
C            IF (IFDENS.EQ.1) THEN ==>>> ==>>>
C  
             END IF ! IFDENS=1
C         
C=======================================================================
C=======================================================================
C
C            Option 02: Traditional Woods-Saxon, here for protons
C    
             IF ((IFDENS.EQ.0).AND.(IFPROT.EQ.1)) THEN
C         
C=======================================================================
C        
C                Item No. 01: The single particle energies
C
                 IF (IF_SPE.EQ.1) THEN 
C
                     DO IEXPER=1,LEVEXP_PROTON(INUCLI)
                        DO ITHEOR=1,LEVTHE_PROTON
C
                           IF (LABEXP_PROTON(INUCLI,IEXPER).EQ.
     *                         LABTHE_PROTON(ITHEOR))     THEN
C
                               ILEVEL=ILEVEL+1
                               LEVPRO=ILEVEL
C
                               IF (IEXPER.EQ.1) INITI1=ILEVEL
C
                               FUNCT0(ILEVEL)=ENETHE_PROTON(ITHEOR)
                               FUNEXP(ILEVEL)=EXPEXP_PROTON(INUCLI,
     *                                                      IEXPER)
                               LEVELS(ILEVEL)=ITHEOR
                               LEXPLS(ILEVEL)=IEXPER
C
                               DEGENE(ILEVEL)
     *                        =REAL(IDEGEX_PROTON(INUCLI,IEXPER))
C
                               WEIGHT_FUNVEC(ILEVEL)
     *                        =DEGENE(ILEVEL)*WEIWEI
C
                               WEISQR_FUNVEC(ILEVEL)
     *                        =
     *                         SQRT(DEGENE(ILEVEL)*WEINUC_PROTON(INUCLI)
     *                                            /SUMWEI_PROTON)
C
                           END IF
C
                        END DO ! ITHEOR
                     END DO ! IEXPER
C
                 END IF ! IF_SPE=1
C
C=======================================================================
C        
C            IF ((IFDENS.EQ.0).AND.(IFPROT.EQ.1)) THEN ==>>> ==>>>
C
C                Item No. 02: The proton r.m.s. radii
C 
                 IF (IF_RAD.EQ.1) THEN
C             
                     ILEVEL=ILEVEL+1
C                    INITI3=ILEVEL
             
                     FUNCT0(ILEVEL)=RMSTHE_PROTON
                     FUNEXP(ILEVEL)=RMSEXP_PROTON(INUCLI)
C             
                     WEISQR_FUNVEC(ILEVEL)
     *              =SQRT(WEIRAD_PROTON(INUCLI)
     *              *     WEINUC_PROTON(INUCLI)/SUMWEI_PROTON)
C             
                 END IF ! IF_RAD=1
C
C=======================================================================
C        
C            IF ((IFDENS.EQ.0).AND.(IFPROT.EQ.1)) THEN ==>>> ==>>>
C
C                Item No. 03: The proton main gap
C          
                 IF (IF_GAP.EQ.1) THEN
C
                     ILEVEL=ILEVEL+1
C                    INITI4=ILEVEL
             
                     FUNCT0(ILEVEL)=GAPTHE_PROTON
                     FUNEXP(ILEVEL)=GAPEXP_PROTON(INUCLI)
C
                     WEISQR_FUNVEC(ILEVEL)=SQRT(WEIGHT_ENEGAP
     *                                    /     SUMWEI_PROTON)
C
                 END IF ! IF_GAP=1
C         
C=======================================================================
C        
C            IF ((IFDENS.EQ.0).AND.(IFPROT.EQ.1)) THEN ==>>> ==>>>
C
C                Item No. 04: The proton "Fermi energy"
C         
                 IF (IF_FER.EQ.1) THEN
C
                     ILEVEL=ILEVEL+1
C                    INITI5=ILEVEL
C
                     FUNCT0(ILEVEL)=FERTHE_PROTON
                     FUNEXP(ILEVEL)=FERMEX_PROTON(INUCLI)
C
                     WEISQR_FUNVEC(ILEVEL)=SQRT(WEIGHT_EFERMI
     *                                    /     SUMWEI_PROTON)
C
                 END IF ! IF_FER=1
C         
C=======================================================================
C        
C            IF ((IFDENS.EQ.0).AND.(IFPROT.EQ.1)) THEN ==>>> ==>>>
C
C                Item No. 05: The proton level densities defined
C                                        separately below-, and above
C          
                 IF (IF_DEN.EQ.1) THEN
C
                     ILEVEL=ILEVEL+1
C                    INITI6=ILEVEL
C
                     FUNCT0(ILEVEL)=DENLOW_PROTON
                     FUNEXP(ILEVEL)=DENSDW_PROTON(INUCLI)
C
                     WEISQR_FUNVEC(ILEVEL)=SQRT(WEIGHT_DENSDW
     *                                    /     SUMWEI_PROTON)
C_______________________________________________________________________
C             
                     ILEVEL=ILEVEL+1
C
                     FUNCT0(ILEVEL)=DENUPP_PROTON
                     FUNEXP(ILEVEL)=DENSUP_PROTON(INUCLI)
C
                     WEISQR_FUNVEC(ILEVEL)=SQRT(WEIGHT_DENSUP
     *                                    /     SUMWEI_PROTON)
C
                 END IF ! IF_DEN=1
C
C=======================================================================
C        
C            IF ((IFDENS.EQ.0).AND.(IFPROT.EQ.1)) THEN ==>>> ==>>>
C
C                Item No. 06: The proton charge densities as functions
C                                               of the radial variable
C         
                 IF (IF_RHO.EQ.1) THEN
*C            
*                     CALL RHOTHE_DENSIT(INUCLI,IZ_FIX,IN_FIX)
*C
*                     DO K_NOYX=1,N_NUCL
*C
*                        I_NOYX=I_NUCL(K_NOYX)
*C
*                        DO ID_RHO=1,ND_RHO
*C
*                           ILEVEL=ILEVEL+1
*C
*                           FUNCT0(ILEVEL)=RHOTHE(ID_RHO,I_NOYX)
*C
*                           FUNEXP(ILEVEL)=RHOEXP(ID_RHO,I_NOYX)
*C
*                           WEISQR_FUNVEC(ILEVEL)=SQRT(WEIGHT_RHODEN)
*C
*                        END DO
*                     END DO
C
                 END IF ! IF_RHO=1
C
C======================================================================= 
C        
C            IF ((IFDENS.EQ.0).AND.(IFPROT.EQ.1)) THEN ==>>> ==>>>
C        
C                Item No. 07: The proton level order as criterion
C                                        Not tested
C        
C                IF (IF_INV.EQ.1) THEN
C             
C                    LEVTOT=LEVPRO+LEVNEU
C             
C                    DO I=1,LEVPRO
C               
C                       IF (I.LT.LEVPRO) THEN
C               
C                           ILWEXP=LEXPLS(I)       !Experimental
C                           IUPEXP=LEXPLS(I+1)
C                           DIFEXP=EXPEXP_PROTON(INUCLI,IUPEXP)
C     *                           -EXPEXP_PROTON(ILWEXP)
C
C                           ILWTHE=LEVELS(I)       !Theoretical
C                           IUPTHE=LEVELS(I+1)               
C                           DIFTHE=ENETHE_PROTON(INUCLI,IUPTHE)
C     *                           -ENETHE_PROTON(ILWTHE)
C             
C                           IF (((DIFEXP.LT.0.).AND.(DIFTHE.GT.0.)).OR.
C     *                         ((DIFEXP.GT.0.).AND.(DIFTHE.LT.0.)))
C                                                                  THEN  
C                               Incorrect order of the levels
C             
C                               ILEVEL=ILEVEL+1
C                  
C                               AUXIL1=ENETHE_PROTON(ILWTHE)*DEGENE(I)
C                               AUXIL2=ENETHE_PROTON(IUPTHE)*DEGENE(I+1)
C                 
C                               FUNVEC(ILEVEL)
C     *                        =
C     *                        (AUXIL2-AUXIL1)*WEIGHT_INVERT
C                  
C                           END IF
C                  
C                       END IF
C             
C                    END DO
C         
C                END IF ! IF_INV=1
C         
C=======================================================================
C=======================================================================
C        
C            IF ((IFDENS.EQ.0).AND.(IFPROT.EQ.1)) THEN ==>>> ==>>>
C
             END IF ! IFDENS=0 & IFPROT=1
C         
C=======================================================================
C
C            Option 03: Traditional Woods-Saxon, here for neutrons
C    
             IF ((IFDENS.EQ.0).AND.(IFNEUT.EQ.1)) THEN
C         
C=======================================================================
C
C                Item No. 01: The single particle energies, neutrons
C
                 IF (IF_SPE.EQ.1) THEN 
C
                     DO IEXPER=1,LEVEXP_NEUTRS(INUCLI)
                        DO ITHEOR=1,LEVTHE_NEUTRS
C
                           IF (LABEXP_NEUTRS(INUCLI,IEXPER).EQ.
     *                         LABTHE_NEUTRS(ITHEOR))     THEN
C
                               ILEVEL=ILEVEL+1
C                              LEVNEU=LEVNEU+1
                               LEVNEU=ILEVEL
C
                               IF (IEXPER.EQ.1) INITI1=ILEVEL
C
                               FUNCT0(ILEVEL)=ENETHE_NEUTRS(ITHEOR)
                               FUNEXP(ILEVEL)
     *                        =EXPEXP_NEUTRS(INUCLI,IEXPER)
C
                               LEVELS(ILEVEL)=ITHEOR
                               LEXPLS(ILEVEL)=IEXPER
C
                               DEGENE(ILEVEL)
     *                        =REAL(IDEGEX_NEUTRS(INUCLI,IEXPER))
C
                               WEISQR_FUNVEC(ILEVEL)
     *                        =SQRT(DEGENE(ILEVEL)
     *                        *WEINUC_NEUTRS(INUCLI)/SUMWEI_NEUTRS)
                 
                           END IF
C
                        END DO
                     END DO
C
                 END IF
C
C=======================================================================
C
C                Item No. 02: The neutron r.m.s. radii
C          
                 IF (IF_RAD.EQ.1) THEN
C
                     ILEVEL=ILEVEL+1
C                    INITI3=ILEVEL
C
                     FUNCT0(ILEVEL)=RMSTHE_NEUTRS
                     FUNEXP(ILEVEL)=RMSEXP_NEUTRS(INUCLI)
C
                     WEISQR_FUNVEC(ILEVEL)
     *              =SQRT(WEIRAD_NEUTRS(INUCLI)
     *              *     WEINUC_NEUTRS(INUCLI)/SUMWEI_NEUTRS)
C
                 END IF ! IF_RAD=1
C
C=======================================================================
C
C                Item No. 03: The neutron main gaps
C          
                 IF (IF_GAP.EQ.1) THEN
C
                     ILEVEL=ILEVEL+1
C                    INITI4=ILEVEL
C
                     FUNCT0(ILEVEL)=GAPTHE_NEUTRS
                     FUNEXP(ILEVEL)=GAPEXP_NEUTRS(INUCLI)
C
                     WEISQR_FUNVEC(ILEVEL)=SQRT(WEIGHT_ENEGAP
     *                                         /SUMWEI_NEUTRS)
C
                 END IF !IF_GAP=1
C         
C======================================================================= 
C
C                Item No. 04: The neutron "Fermi energies"
C         
                 IF (IF_FER.EQ.1) THEN
C
                     ILEVEL=ILEVEL+1
C                    INITI5=ILEVEL

                     FUNCT0(ILEVEL)=FERTHE_NEUTRS
                     FUNEXP(ILEVEL)=FERMEX_NEUTRS(INUCLI)
C
                     WEISQR_FUNVEC(ILEVEL)=SQRT(WEIGHT_EFERMI
     *                                         /SUMWEI_NEUTRS)
C
                 END IF ! IF_FER=1
C         
C=======================================================================
C
C                Item No. 05: The neutron densities defined
C                                        separately below-, and above
C          
                 IF (IF_DEN.EQ.1) THEN
C
                     ILEVEL=ILEVEL+1
C                    INITI6=ILEVEL
C
                     FUNCT0(ILEVEL)=DENLOW_NEUTRS
                     FUNEXP(ILEVEL)=DENSDW_NEUTRS(INUCLI)
C
                     WEISQR_FUNVEC(ILEVEL)=SQRT(WEIGHT_DENSDW
     *                                         /SUMWEI_NEUTRS)
C_______________________________________________________________________
C             
                     ILEVEL=ILEVEL+1
C
                     FUNCT0(ILEVEL)=DENUPP_NEUTRS
                     FUNEXP(ILEVEL)=DENSUP_NEUTRS(INUCLI)
C
                     WEISQR_FUNVEC(ILEVEL)=SQRT(WEIGHT_DENSUP
     *                                         /SUMWEI_NEUTRS)
C
                 END IF ! IF_DEN=1
C
C=======================================================================
C
C                Item No. 06: The neutron charge densities as functions
C                             of the radial variable not available from
C                                                            experiment         
C                IF (IF_RHO.EQ.1) THEN
C         
C                END IF !IF_RHO=1
C
C======================================================================= 
C        
C                Item No. 07: The neutron level order as criterion
C                                        Not tested, Neutrons missing
C        
C                IF (IF_INV.EQ.1) THEN
C             
C                    LEVTOT=LEVPRO+LEVNEU
C             
C                    DO I=1,LEVPRO
C               
C                       IF (I.LT.LEVPRO) THEN
C               
C                           ILWEXP=LEXPLS(I)       !Experimental
C                           IUPEXP=LEXPLS(I+1)
C                           DIFEXP=EXPEXP_PROTON(INUCLI,IUPEXP)
C     *                           -EXPEXP_PROTON(ILWEXP)
C
C                           ILWTHE=LEVELS(I)       !Theoretical
C                           IUPTHE=LEVELS(I+1)               
C                           DIFTHE=ENETHE_PROTON(INUCLI,IUPTHE)
C     *                           -ENETHE_PROTON(ILWTHE)
C             
C                           IF (((DIFEXP.LT.0.) .AND. 
C     *                          (DIFTHE.GT.0.))    .OR.
C     *                         ((DIFEXP.GT.0.) .AND.
C     *                          (DIFTHE.LT.0.)))   THEN 
C
C                               Incorrect order of levels
C             
C                               ILEVEL=ILEVEL+1
C                  
C                               AUXIL1
C     *                        =ENETHE_PROTON(ILWTHE)*DEGENE(I)
C                               AUXIL2
C     *                        =ENETHE_PROTON(IUPTHE)*DEGENE(I+1)
C                 
C                               FUNVEC(ILEVEL)
C     *                        =
C     *                        (AUXIL2-AUXIL1)*WEIGHT_INVERT
C                  
C                           END IF
C             
C                       END IF
C             
C                    END DO
C         
C                END IF ! IF_INV=1
C         
C=======================================================================
C=======================================================================
C        
C            IF ((IFDENS.EQ.0).AND.(IFNEUT.EQ.1)) THEN ==>>> ==>>>
C
             END IF ! IFDENS=0 & IFNEUT=1
C
             IF (LOGWRI.GT.0)
     *
     *       WRITE(LOGFIL,'(12X,''In FUNMIN - Finished I_FLAG=1 '',
     *                          ''option '',
     *                          ''single-function run, INUCLI= '',I2)')
     *                                                 INUCLI
C         
C=======================================================================
C        
             IPARTL=ILEVEL ! Partial dimension of FUNCT0, 
C                                                 FUNEXP and FUNVEC
C
C=======================================================================       
C============================ J A C O B I A N ========================== 
C============================ J A C O B I A N ========================== 
C============================ J A C O B I A N ========================== 
C============================ J A C O B I A N ========================== 
C=======================================================================
C         
C     DO INUCLI=1,LDNUCL  ==>>>  ==>>>
C     
C        IF (ITAKNU(INUCLI).EQ.1) THEN  ==>>>  ==>>>
C
             IF (I_FLAG.EQ.2) THEN ! I_FLAG=2 =>> Jacobian
C
C=======================================================================
C          
                 DO IPARAM=1,NDPARS
                    DERIVF(IPARAM)=0.0
                    DO INDEXX=INITI1,IPARTL
                       FUNJAC(INDEXX,IPARAM)=0.
                    END DO
                 END DO
C
C=======================================================================
C        
                 IACTIV=0  !Running over the number of parameters
C
                 DO IPARAM=1,NDPARS
         
                    IF (IFTAKE(IPARAM).EQ.1) THEN
C
                        PARPOT(IPARAM)=PARPOT(IPARAM)+DPARAM(IPARAM)
C
                        IF (LOGWRI.GT.4) THEN
C
                            WRITE(LOGFIL,'(12X,''Entering WS_RUN '',
     *                                         ''from FUNMIN [2] '',
     *                                         ''p+dp='',f10.5,1X,
     *                                         ''p='',f10.5,6X,
     *                                         ''IPARAM='',I2,1X,
     *                                         ''PARAM='',A)')
     *                                   PARPOT(IPARAM),PARPOT(IPARAM)
     *                                                 -DPARAM(IPARAM),
     *                                                         IPARAM,
     *                                           TITLES_EXPORT(IPARAM)
c                           WRITE(LOGFIL,'()')
                        END IF
C
C
C                       Below - I_MODE=1 means: Calculate the derivatives
                                                                I_MODE=1
                        CALL WS_RUN(INUCLI,CHISQU_CHOICE,IABORT,I_MODE,
     *                              I_FLAG,CHISQU_AUXIL3,CHISQU_AUXIL4)
C
                        IACTIV=IACTIV+1
C_______________________________________________________________________
C               
                        CHISQU_AUXILP(IPARAM,INUCLI)=CHISQU_AUXIL3 ! P
                        CHISQU_AUXILN(IPARAM,INUCLI)=CHISQU_AUXIL4 ! N
C_______________________________________________________________________
C    
                        ACTION(IPARAM)='*'
C
C=======================================================================
C
C                       Option 01: Density-dependent spin-orbit
C                      
                        IF (IFDENS.EQ.1 .OR.
     *                     (IFDENS.EQ.0 .AND. IFBOTH.EQ.1)) THEN 
C
C=======================================================================
C
C                           Item No. 01: The single particle energies
C           
                            IF (IF_SPE.EQ.1) THEN
C
                                IAUXIL=LEVNEU
C	       
                                DO INDEXX=INITI1,LEVPRO
C	       
                                   IWHICH=LEVELS(INDEXX)
C  
                                   FUNCT1(INDEXX)=ENETHE_PROTON(IWHICH)
C
                                   DERIVF(IACTIV)
     *                            =
     *                            (FUNCT1(INDEXX)-FUNCT0(INDEXX))
     *                            /DPARAM(IPARAM)
C
                                   FUNJAC(INDEXX,IACTIV)=DERIVF(IACTIV)
C
                                END DO
C_______________________________________________________________________
C
                                DO INDEXX=INITI2,LEVNEU
C
                                   IWHICH=LEVELS(INDEXX)
C
                                   FUNCT1(INDEXX)=ENETHE_NEUTRS(IWHICH)
C
                                   DERIVF(IACTIV)
     *                            =
     *                            (FUNCT1(INDEXX)-FUNCT0(INDEXX))                      
     *                            /DPARAM(IPARAM)
C
                                   FUNJAC(INDEXX,IACTIV)=DERIVF(IACTIV)
C                   
                                END DO
C
                            END IF ! IF_SPE=1
C
C=======================================================================
C
C                           Item No. 02: The proton and neutron
C                                                  r.m.s. radii
                            IF (IF_RAD.EQ.1) THEN 
C
                                IAUXIL=IAUXIL+1
C
                                FUNCT1(IAUXIL)=RMSTHE_PROTON
C
                                FUNJAC(IAUXIL,IACTIV)
     *                         =
     *                         (FUNCT1(IAUXIL)-FUNCT0(IAUXIL))
     *                         /DPARAM(IPARAM)
C_______________________________________________________________________
C            
                                IAUXIL=IAUXIL+1
C
                                FUNCT1(IAUXIL)=RMSTHE_NEUTRS
C
                                FUNJAC(IAUXIL,IACTIV)
     *                         =
     *                         (FUNCT1(IAUXIL)-FUNCT0(IAUXIL))
     *                         /DPARAM(IPARAM)  
C
                            END IF ! IF_RAD=1  
C
C=======================================================================
C
C                           Item No. 03: The proton and neutron
C                                                     main gaps
                            IF (IF_GAP.EQ.1) THEN 
C
                                IAUXIL=IAUXIL+1
C
                                FUNCT1(IAUXIL)=GAPTHE_PROTON
C
                                FUNJAC(IAUXIL,IACTIV)
     *                         =
     *                         (FUNCT1(IAUXIL)-FUNCT0(IAUXIL))
     *                         /DPARAM(IPARAM)
C_______________________________________________________________________
C           
                                IAUXIL=IAUXIL+1
C
                                FUNCT1(IAUXIL)=GAPTHE_NEUTRS
C
                                FUNJAC(IAUXIL,IACTIV)
     *                         =
     *                         (FUNCT1(IAUXIL)-FUNCT0(IAUXIL))
     *                         /DPARAM(IPARAM)
C
                            END IF ! IF_GAP=1   
C
C=======================================================================
C
C                           Item No. 04: The proton and neutron 
C                                              "Fermi energies"
C          
                            IF (IF_FER.EQ.1) THEN 
C
                                IAUXIL=IAUXIL+1
C
                                FUNCT1(IAUXIL)=FERTHE_PROTON
C
                                FUNJAC(IAUXIL,IACTIV)
     *                         =
     *                         (FUNCT1(IAUXIL)-FUNCT0(IAUXIL))
     *                         /DPARAM(IPARAM)
C_______________________________________________________________________
C          
                                IAUXIL=IAUXIL+1
C
                                FUNCT1(IAUXIL)=FERTHE_NEUTRS
C
                                FUNJAC(IAUXIL,IACTIV)
     *                         =
     *                         (FUNCT1(IAUXIL)-FUNCT0(IAUXIL))
     *                         /DPARAM(IPARAM)
C
                            END IF  ! IF_FER=1 
C
C=======================================================================
C
C                           Item No. 05: The proton and neutron 
C                                        densities defined separately
C                                        below-, and above the main gaps
C           
                            IF (IF_DEN.EQ.1) THEN 
C
                                IAUXIL=IAUXIL+1
C
                                FUNCT1(IAUXIL)=DENLOW_PROTON
C
                                FUNJAC(IAUXIL,IACTIV)
     *                         =
     *                         (FUNCT1(IAUXIL)-FUNCT0(IAUXIL))
     *                         /DPARAM(IPARAM)
C_______________________________________________________________________
C          
                                IAUXIL=IAUXIL+1
C
                                FUNCT1(IAUXIL)=DENUPP_PROTON
C
                                FUNJAC(IAUXIL,IACTIV)
     *                         =
     *                         (FUNCT1(IAUXIL)-FUNCT0(IAUXIL))
     *                         /DPARAM(IPARAM)               
C_______________________________________________________________________
C          
                                IAUXIL=IAUXIL+1
C
                                FUNCT1(IAUXIL)=DENLOW_NEUTRS
C
                                FUNJAC(IAUXIL,IACTIV)
     *                         =
     *                         (FUNCT1(IAUXIL)-FUNCT0(IAUXIL))
     *                         /DPARAM(IPARAM)
C_______________________________________________________________________
C
                                IAUXIL=IAUXIL+1
C
                                FUNCT1(IAUXIL)=DENUPP_NEUTRS
C
                                FUNJAC(IAUXIL,IACTIV)
     *                         =
     *                         (FUNCT1(IAUXIL)-FUNCT0(IAUXIL))
     *                         /DPARAM(IPARAM)
C
                            END IF  ! IF_DEN=1 
C
C=======================================================================
C 
                        END IF ! IFDENS=1
C
C=======================================================================
C
C                       Option 02: Traditional Woods-Saxon, 
C                                                    here for protons
C
                        IF ((IFDENS.EQ.0).AND.(IFPROT.EQ.1)) THEN
C
C=======================================================================
C
C                           Item No. 01: The single particle energies
C                     
                            IF (IF_SPE.EQ.1) THEN
C
                                IAUXIL=LEVPRO
C	       
                                DO INDEXX=INITI1,LEVPRO
C	       
                                   IWHICH=LEVELS(INDEXX)
C  
                                   FUNCT1(INDEXX)=ENETHE_PROTON(IWHICH)
C
                                   DERIVF(IACTIV)
     *                            =
     *                            (FUNCT1(INDEXX)-FUNCT0(INDEXX))
     *                            /DPARAM(IPARAM)
C
                                   FUNJAC(INDEXX,IACTIV)=DERIVF(IACTIV)
C
                                END DO
C
                            END IF ! IF_SPE=1
C
C=======================================================================
C
C                           Item No. 02: The proton r.m.s. radii
C           
                            IF (IF_RAD.EQ.1) THEN 
C
                                IAUXIL=IAUXIL+1
C
                                FUNCT1(IAUXIL)=RMSTHE_PROTON
C
                                FUNJAC(IAUXIL,IACTIV)
     *                         =
     *                         (FUNCT1(IAUXIL)-FUNCT0(IAUXIL))
     *                         /DPARAM(IPARAM)               
C
                            END IF ! IF_RAD=1  
C
C=======================================================================
C
C                           Item No. 03: The proton main gaps
C           
                            IF (IF_GAP.EQ.1) THEN 
C
                                IAUXIL=IAUXIL+1
C
                                FUNCT1(IAUXIL)=GAPTHE_PROTON
C
                                FUNJAC(IAUXIL,IACTIV)
     *                         =
     *                         (FUNCT1(IAUXIL)-FUNCT0(IAUXIL))
     *                         /DPARAM(IPARAM)
C
                            END IF ! IF_GAP=1   
C
C=======================================================================
C
C                           Item No. 04: The proton "Fermi energies"
C          
                            IF (IF_FER.EQ.1) THEN 
C
                                IAUXIL=IAUXIL+1
C
                                FUNCT1(IAUXIL)=FERTHE_PROTON
C
                                FUNJAC(IAUXIL,IACTIV)
     *                         =
     *                         (FUNCT1(IAUXIL)-FUNCT0(IAUXIL))
     *                         /DPARAM(IPARAM)
C
                            END IF ! IF_FER=1 
C 
C=======================================================================
C
C                           Item No. 05: The proton densities defined
C                                        separately below-, and above
C           
                            IF (IF_DEN.EQ.1) THEN 
C
                                IAUXIL=IAUXIL+1
C
                                FUNCT1(IAUXIL)=DENLOW_PROTON
C
                                FUNJAC(IAUXIL,IACTIV)
     *                         =
     *                         (FUNCT1(IAUXIL)-FUNCT0(IAUXIL))
     *                         /DPARAM(IPARAM)
C_______________________________________________________________________
C         
                                IAUXIL=IAUXIL+1
C
                                FUNCT1(IAUXIL)=DENUPP_PROTON
C
                                FUNJAC(IAUXIL,IACTIV)
     *                         =
     *                         (FUNCT1(IAUXIL)-FUNCT0(IAUXIL))
     *                         /DPARAM(IPARAM)              
C
                            END IF ! IF_DEN=1 
C
C=======================================================================
C
C                           Item No. 06: The proton charge densities 
C                                        as functions of the radial 
C                                                           variable
                            IF (IF_RHO.EQ.1) THEN
C            
*                                CALL RHOTHE_DENSIT(INUCLI,IZ_FIX,IN_FIX)
*C
*                                DO K_NOYX=1,N_NUCL
*C
*                                   I_NOYX=I_NUCL(K_NOYX)
*C
*                                   DO ID_RHO=1,ND_RHO
*C
*                                      IAUXIL=IAUXIL+1
                   
*                                      FUNCT1(IAUXIL)=RHOTHE(ID_RHO,
*     *                                                      I_NOYX)
*C
*                                      FUNJAC(IAUXIL,IACTIV)
*     *                               =
*     *                               (FUNCT1(IAUXIL)-FUNCT0(IAUXIL))
*     *                               /DPARAM(IPARAM)   
*C                   
*                                   END DO
*                                END DO
C
                            END IF
C
C=======================================================================
C
                        END IF ! IFDENS=0 & IFPROT=1
C
C=======================================================================
C=======================================================================
C=======================================================================
C
C                       Option 03: Traditional Woods-Saxon, for neutrons
C 
                        IF ((IFDENS.EQ.0) .AND. (IFNEUT.EQ.1)) THEN   
C
C                           Item No. 01: The single particle energies
C                                                   
                            IF (IF_SPE.EQ.1) THEN
C
                                IAUXIL=LEVNEU
C	       
                                DO INDEXX=INITI1,LEVNEU
C	       
                                   IWHICH=LEVELS(INDEXX)
C  
                                   FUNCT1(INDEXX)=ENETHE_NEUTRS(IWHICH)
C
                                   DERIVF(IACTIV)
     *                            =
     *                            (FUNCT1(INDEXX)-FUNCT0(INDEXX))
     *                            /DPARAM(IPARAM)
C
                                   FUNJAC(INDEXX,IACTIV)=DERIVF(IACTIV)
C
                                END DO
C
                            END IF ! IF_SPE=1
C
C=======================================================================
C
C                           Item No. 02: The neutron r.m.s. radii
C           
                            IF (IF_RAD.EQ.1) THEN 
C
                                IAUXIL=IAUXIL+1
C
                                FUNCT1(IAUXIL)=RMSTHE_NEUTRS
C
                                FUNJAC(IAUXIL,IACTIV)
     *                         =
     *                         (FUNCT1(IAUXIL)-FUNCT0(IAUXIL))
     *                         /DPARAM(IPARAM)             
C
                            END IF  ! IF_RAD=1  
C
C=======================================================================
C
C                           Item No. 03: The neutron main gaps
C           
                            IF (IF_GAP.EQ.1) THEN 
C
                                IAUXIL=IAUXIL+1
C
                                FUNCT1(IAUXIL)=GAPTHE_NEUTRS
C
                                FUNJAC(IAUXIL,IACTIV)
     *                         =
     *                         (FUNCT1(IAUXIL)-FUNCT0(IAUXIL))
     *                         /DPARAM(IPARAM)
C
                            END IF  ! IF_GAP=1   
C
C======================================================================= 
C
C                           Item No. 04: The neutron "Fermi energies"
C          
                            IF (IF_FER.EQ.1) THEN 
C
                                IAUXIL=IAUXIL+1
C
                                FUNCT1(IAUXIL)=FERTHE_NEUTRS
C
                                FUNJAC(IAUXIL,IACTIV)
     *                         =
     *                         (FUNCT1(IAUXIL)-FUNCT0(IAUXIL))
     *                         /DPARAM(IPARAM)
C
                            END IF  ! IF_FER=1 
C
C=======================================================================
C
C                           Item No. 05: The neutron densities defined
C                                         separately below-, and above
C           
                            IF (IF_DEN.EQ.1) THEN 
C
                                IAUXIL=IAUXIL+1
C 
                                FUNCT1(IAUXIL)=DENLOW_NEUTRS
C
                                FUNJAC(IAUXIL,IACTIV)
     *                         =
     *                         (FUNCT1(IAUXIL)-FUNCT0(IAUXIL))
     *                         /DPARAM(IPARAM)
C_______________________________________________________________________
C          
                                IAUXIL=IAUXIL+1
C 
                                FUNCT1(IAUXIL)=DENUPP_NEUTRS
C
                                FUNJAC(IAUXIL,IACTIV)
     *                         =
     *                         (FUNCT1(IAUXIL)-FUNCT0(IAUXIL))
     *                         /DPARAM(IPARAM)              
C
                            END IF  ! IF_DEN=1 
C
C=======================================================================
C
                        END IF ! IFDENS=0 & IFNEUT=1
C
C=======================================================================
C
                        ACTION(IPARAM)=' '
C
                        PARPOT(IPARAM)=PARPOT(IPARAM)-DPARAM(IPARAM)
C        
C=======================================================================
C            
                    END IF ! (IFTAKE(IPARAM).EQ.1)
C 
                 END DO ! IPARAM=1,NDPARS
C        
C=======================================================================
C         
                 IF (IPARTL.NE.IAUXIL) THEN
C
                     WRITE(NOUTPT,'(''Something went wrong in FUNMIN:'',
     *                              ''IPARTL='',I3,1X,
     *                              ''is different from IAUXIL='',I3,
     *                              '' -> They should be equal'')') 
     *                                IPARTL,IAUXIL
C
                     STOP 'STOP: Something went wrong in FUNMIN (I)'
C
                 END IF
C
C========================================================================         
C       
                 IF (LOGWRI.GT.0)
     *
     *           WRITE(LOGFIL,'(/,12x,''In FUNMIN - Finished I_FLAG=2'',
     *                                '' OPTION, '',
     *                                ''Jacobian run, INUCLI= '',I2)') 
     *                                                INUCLI
C
             END IF ! I_FLAG=2 ==> End of the Jacobian Storing
C
C========================================================================         
C
         END IF ! For active nuclei ITAKNU(INUCLI).EQ.1
      
      END DO ! Over all nuclei INUCLI
C
C========================================================================
C========================================================================
C========================================================================
C
      ITOTAL=ILEVEL ! Total dimension of FUNVEC and FUNJAC
C                     Should be equal to LDFUNC
C
      IF (LOGWRI.GT.0) THEN
          WRITE(LOGFIL,'(12X,''Out of the do-loop over INUCLI '',
     *                       ''in FUNMIN with number of active '',
     *                       ''parameters= '',I3)') ITOTAL
      END IF
C
C========================================================================
C       
      IF (ITOTAL.NE.LDFUNC) THEN
C
          WRITE(NOUTPT,'(''Something went wrong in FUNMIN: ITOTAL='',I3,
     *                   '' is different from LDFUNC='',I3,
     *                   '' -> They should be equal!!'')') ITOTAL,LDFUNC
C
          STOP 'STOP: Something went wrong in FUNMIN (II)'
C
      END IF
C________________________________________________________________________
C       
      IF (I_FLAG.EQ.2 .AND. IAUXIL.NE.LDFUNC) THEN
C
          WRITE(NOUTPT,'(''Something went wrong in FUNMIN: IAUXIL='',I3,
     *                   '' is different from LDFUNC='',I3,
     *                   '' -> They should be equal!!'')') IAUXIL,LDFUNC
C
          STOP 'STOP: Something went wrong in FUNMIN (III)'
C
      END IF
C
C=======================================================================
C      
      IF (I_FLAG.EQ.1) THEN
C
C         Storing FUNVEC(1:LDFUNC), which enters to LEVMAR routine
C           
          DO INDEXZ=1,LDFUNC
C
             FUNVEC(INDEXZ)=(FUNCT0(INDEXZ)-FUNEXP(INDEXZ))
     *                     * WEISQR_FUNVEC(INDEXZ)
          END DO
C
C=======================================================================
C         Printing actual parameters and corresponding Chi^2 (AS ORIG.)
C=======================================================================
C        
          VALCHI=CHISQU_TOTALS  ! chi^2=chi^(p)
C
          IF (IFDENS.EQ.0 .AND. IFPROT.EQ.1) 
     *                          HIAUXI(IDEFCN)=CHISQU_PROTON
          IF (IFDENS.EQ.0 .AND. IFNEUT.EQ.1) 
     *                          HIAUXI(IDEFCN)=CHISQU_NEUTRS
C
          IF (IFDENS.EQ.1 .OR. (IFDENS.EQ.0 .AND. IFBOTH.EQ.1)) THEN
              HIAUXI(IDEFCN)=VALCHI
          END IF
C
                                    I_MODE=0
          CALL INPRIN(IDEFCN,IACTIV,I_MODE,CHISQU_PROTON,
     *                                     CHISQU_NEUTRS,
     *                       CHIENE_WEIPRO,CHIENE_WEINEU,
     *                       CHIRAD_WEIPRO,CHIRAD_WEINEU,
     *                                     CHIRHO_WEIPRO)
          IF (NUCACT.GT.1) THEN
                                               I_MODE=0
              CALL INPRIN_GLOFIT(IDEFCN,IACTIV,I_MODE,CHIENE_PROTON,
     *                                  CHIENE_NEUTRS,CHIRAD_PROTON,
     *                                                CHIRAD_NEUTRS)
          END IF
C_______________________________________________________________________ 
C        
          IDEFCN=IDEFCN+1
C_______________________________________________________________________ 
C 
      END IF ! I_FLAG=1
C
C=======================================================================
C
      IF (I_FLAG.EQ.2) THEN 
C
C=======================================================================
C         Introducing weights to Jacobian FUNJAC(ifunc,ipars)
C=======================================================================
C           
          DO INDEXP=1,IACTIV
             DO INDEXF=1,LDFUNC
                FUNJAC(INDEXF,INDEXP)=FUNJAC(INDEXF,INDEXP)
     *                               *WEISQR_FUNVEC(INDEXF)
             END DO
          END DO
C
C=======================================================================
C
           INDEXP_AUXILI=0
      
           IF (IFDENS.EQ.0 .AND. IFPROT.EQ.1) THEN
               DO IPARAM=1,6
                  INDEXP_AUXILI=INDEXP_AUXILI+1
                  IVECTR(INDEXP_AUXILI)=IPARAM
               END DO
           END IF
C
           IF (IFDENS.EQ.0 .AND. IFNEUT.EQ.1) THEN
               DO IPARAM=21,26
                  INDEXP_AUXILI=INDEXP_AUXILI+1
                  IVECTR(INDEXP_AUXILI)=IPARAM
               END DO
           END IF
C
           IF (IFDENS.EQ.0.AND.IFBOTH.EQ.1) THEN
C
               IF (IFK_VC.EQ.1) THEN
                   INDEXP_AUXILI=INDEXP_AUXILI+1
                   IVECTR(INDEXP_AUXILI)=51   ! V_o central
                   INDEXP_AUXILI=INDEXP_AUXILI+1
                   IVECTR(INDEXP_AUXILI)=52   ! kappa V_o central
               ELSE
                   INDEXP_AUXILI=INDEXP_AUXILI+1 
                   IVECTR(INDEXP_AUXILI)=1    ! V_o central protons
                   INDEXP_AUXILI=INDEXP_AUXILI+1
                   IVECTR(INDEXP_AUXILI)=21   ! V_o central neutrons
               END IF
C
               IF (IFK_RC.EQ.1) THEN
                   INDEXP_AUXILI=INDEXP_AUXILI+1
                   IVECTR(INDEXP_AUXILI)=53   ! r_o central
                   INDEXP_AUXILI=INDEXP_AUXILI+1
                   IVECTR(INDEXP_AUXILI)=54   ! kappa r_o central
               ELSE
                   INDEXP_AUXILI=INDEXP_AUXILI+1
                   IVECTR(INDEXP_AUXILI)=2    ! r_o central protons
                   INDEXP_AUXILI=INDEXP_AUXILI+1
                   IVECTR(INDEXP_AUXILI)=22   ! r_o central neutron
               END IF
C
               IF (IFK_AC.EQ.1) THEN
                   INDEXP_AUXILI=INDEXP_AUXILI+1
                   IVECTR(INDEXP_AUXILI)=55   ! a_o central
                   INDEXP_AUXILI=INDEXP_AUXILI+1
                   IVECTR(INDEXP_AUXILI)=56   ! kappa a_o central
               ELSE
                   INDEXP_AUXILI=INDEXP_AUXILI+1
                   IVECTR(INDEXP_AUXILI)=3    ! a_o central protons
                   INDEXP_AUXILI=INDEXP_AUXILI+1
                   IVECTR(INDEXP_AUXILI)=23   ! a_o central neutrons
               END IF
C
               IF (IFK_VS.EQ.1) THEN
                   INDEXP_AUXILI=INDEXP_AUXILI+1
                   IVECTR(INDEXP_AUXILI)=57   ! V_o pure WS-SO
                   INDEXP_AUXILI=INDEXP_AUXILI+1
                   IVECTR(INDEXP_AUXILI)=58   ! kappa V_o pure WS-SO
               ELSE
                   INDEXP_AUXILI=INDEXP_AUXILI+1
                   IVECTR(INDEXP_AUXILI)=4    ! V_o pure WS-SO protons
                   INDEXP_AUXILI=INDEXP_AUXILI+1
                   IVECTR(INDEXP_AUXILI)=24   ! V_o pure WS-SO neutrons
               END IF
C
               IF (IFK_RS.EQ.1) THEN
                   INDEXP_AUXILI=INDEXP_AUXILI+1
                   IVECTR(INDEXP_AUXILI)=59   ! r_o pure WS-SO
                   INDEXP_AUXILI=INDEXP_AUXILI+1
                   IVECTR(INDEXP_AUXILI)=60   ! kappa r_o pure WS-SO
               ELSE
                   INDEXP_AUXILI=INDEXP_AUXILI+1
                   IVECTR(INDEXP_AUXILI)=5    ! r_o pure WS-SO protons
                   INDEXP_AUXILI=INDEXP_AUXILI+1
                   IVECTR(INDEXP_AUXILI)=25   ! r_o pure WS-SO neutrons
               END IF
C
               IF (IFK_AS.EQ.1) THEN
                   INDEXP_AUXILI=INDEXP_AUXILI+1
                   IVECTR(INDEXP_AUXILI)=61   ! a_o pure WS-SO
                   INDEXP_AUXILI=INDEXP_AUXILI+1
                   IVECTR(INDEXP_AUXILI)=62   ! kappa a_o pure WS-SO
               ELSE
                   INDEXP_AUXILI=INDEXP_AUXILI+1
                   IVECTR(INDEXP_AUXILI)=6    ! a_o pure WS-SO protons
                   INDEXP_AUXILI=INDEXP_AUXILI+1
                   IVECTR(INDEXP_AUXILI)=26   ! a_o pure WS-SO neutrons
               END IF
C
           END IF
C
C=======================================================================
C      
           IF (IFDENS.EQ.1) THEN
C
               IF (IFK_VC.EQ.1) THEN
                   INDEXP_AUXILI=INDEXP_AUXILI+1
                   IVECTR(INDEXP_AUXILI)=51   ! V_o central
                   INDEXP_AUXILI=INDEXP_AUXILI+1
                   IVECTR(INDEXP_AUXILI)=52   ! kappa V_o central
               ELSE
                   INDEXP_AUXILI=INDEXP_AUXILI+1 
                   IVECTR(INDEXP_AUXILI)=1    ! V_o central protons
                   INDEXP_AUXILI=INDEXP_AUXILI+1
                   IVECTR(INDEXP_AUXILI)=21   ! V_o central neutrons
               END IF
C
               IF (IFK_RC.EQ.1) THEN
                   INDEXP_AUXILI=INDEXP_AUXILI+1
                   IVECTR(INDEXP_AUXILI)=53   ! r_o central
                   INDEXP_AUXILI=INDEXP_AUXILI+1
                   IVECTR(INDEXP_AUXILI)=54   ! kappa r_o central
               ELSE
                   INDEXP_AUXILI=INDEXP_AUXILI+1
                   IVECTR(INDEXP_AUXILI)=2    ! r_o central protons
                   INDEXP_AUXILI=INDEXP_AUXILI+1
                   IVECTR(INDEXP_AUXILI)=22   ! r_o central neutron
               END IF
C
               IF (IFK_AC.EQ.1) THEN
                   INDEXP_AUXILI=INDEXP_AUXILI+1
                   IVECTR(INDEXP_AUXILI)=55   ! a_o central
                   INDEXP_AUXILI=INDEXP_AUXILI+1
                   IVECTR(INDEXP_AUXILI)=56   ! kappa a_o central
               ELSE
                   INDEXP_AUXILI=INDEXP_AUXILI+1
                   IVECTR(INDEXP_AUXILI)=3    ! a_o central protons
                   INDEXP_AUXILI=INDEXP_AUXILI+1
                   IVECTR(INDEXP_AUXILI)=23   ! a_o central neutrons
               END IF
C
               DO IPARAM=39,42
                  INDEXP_AUXILI=INDEXP_AUXILI+1
                  IVECTR(INDEXP_AUXILI)=IPARAM
               END DO
C          
               IF (IFTENS.EQ.1 .AND. ISORBT.EQ.1) THEN
                   DO IPARAM=43,46
                      INDEXP_AUXILI=INDEXP_AUXILI+1
                      IVECTR(INDEXP_AUXILI)=IPARAM
                   END DO
               END IF
C          
               IF (IFTENS.EQ.1 .AND. ICENTT.EQ.1) THEN
                   DO IPARAM=47,50
                      INDEXP_AUXILI=INDEXP_AUXILI+1
                      IVECTR(INDEXP_AUXILI)=IPARAM
                   END DO
               END IF
C          
           END IF
C
C=======================================================================
C          Calculating the chi^2 gradient
C=======================================================================
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
     *                         * WEISQR_FUNVEC(INDEXF)
C
                 CHISQU_GRADNT(INDEXP)=CHISQU_GRADNT(INDEXP)
     *                                +FUNAUX(INDEXF)
     *                                *FUNJAC(INDEXF,INDEXP)
C
             END DO
C         
             CHISQU_GRADNT(INDEXP)=2*CHISQU_GRADNT(INDEXP)
             CHISQU_GRDNRM=CHISQU_GRDNRM+(CHISQU_GRADNT(INDEXP))**2
C         
          END DO
C________________________________________________________________________
C        
          CHISQU_GRDNRM=SQRT(CHISQU_GRDNRM)
C        
CID          IF (LOGWRI.GT.0) THEN
C            
              JACTIV=0
C                
              IF (IFDENS.EQ.0 .AND. IFBOTH.EQ.0) THEN
                  WRITE(LOGAUX,'(3x,''G:'',148X,''  '',$)')
              END IF
C
              IF (IFDENS.EQ.0 .AND. IFBOTH.EQ.1) THEN
                  WRITE(LOGAUX,'(3x,''G:'',242X,''  '',$)')
              END IF
C
              IF (IFDENS.EQ.1 .AND. IFTENS.EQ.0) THEN 
                  WRITE(LOGAUX,'(3x,''G:'',220X,''  '',$)')
              END IF
C
              IF (IFDENS.EQ.1 .AND. IFTENS.EQ.1) THEN 
                  WRITE(LOGAUX,'(3x,''G:'',264X,''  '',$)')
              END IF
C                
              DO I=1,INDEXP_AUXILI
C
                 IPARAM=IVECTR(I)
C
                 IF (IFTAKE(IPARAM).EQ.1) THEN
                     JACTIV=IACTIV_VECTOR(IPARAM)
                     WRITE(LOGAUX,'(E11.3,$)') CHISQU_GRADNT(JACTIV)
                     CHISQU_GRADIE(IPARAM)=CHISQU_GRADNT(JACTIV)
                 ELSE
                     CHISQU_CNSTNT=0.0
                     WRITE(LOGAUX,'(E11.3,$)') CHISQU_CNSTNT
                     CHISQU_GRADIE(IPARAM)=CHISQU_CNSTNT
                 END IF
C
              END DO
C                
              WRITE(LOGAUX,'(F18.4)') CHISQU_GRDNRM
C            
CID          END IF
C        
          IF (ISCREN.NE.1) GO TO 1
C        
          WRITE(0,'(3X,''g:'',2X,<IACTIV>(E10.3,1X),F15.4)')
     *                        (CHISQU_GRADNT(I),I=1,IACTIV),
     *                         CHISQU_GRDNRM
     
   1      CONTINUE
C       
      END IF
C
C=======================================================================
C
c     IF (LOGWRI.GT.0) THEN
C
c         WRITE(LOGFIL,'(12X,''ACTIVE Parameters in FUNMIN'',
c    *                       ''(IDEFCN,IPARAM,IACTIV,'',
c    *                       ''ARGPAR(IACTIV),PARPOT(IACTIV))'')')
C       
c         KK=0
C       
c         DO IPARAM=1,NDPARS  
C               
c            IF (IFTAKE(IPARAM).EQ.1) THEN
c                KK=KK+1           
c                WRITE(LOGFIL,'(T20,''PARAM ='',A,
c    *                          T38,3I7,2F15.4)') 
c    *                          TITLES_EXPORT(IPARAM),IDEFCN,IPARAM,KK,
c    *                                        ARGPAR(KK),PARPOT(IPARAM)
c            END IF
C                  
c         END DO
C
c     END IF
C
C=======================================================================
C
      IF (LOGWRI.GT.0) THEN
          WRITE(LOGFIL,'(/,9X,''Exiting  FUNMIN'')')
      END IF 
C
C=======================================================================
C        
      CALL CPUTIM('FUNMIN',0)
C
C=======================================================================
C                          
      RETURN
      END 
C
C=======================================================================
C=======================================================================
C
      SUBROUTINE VERIFY(ARGMNT,IABORT)
C      
      INCLUDE   'MATDIM/NDPARS.f'
C
      CHARACTER
     *          ACTION*1,CALLED*6     
C
      DIMENSION
     *          ARGMNT(1:NDPARS)
C
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
     *       /EXTKAP/ V0CMIN_KAPPAR,V0CMAX_KAPPAR,XKCMIN_KAPPAR,
     *                XKCMAX_KAPPAR,A0CMIN_KAPPAR,A0CMAX_KAPPAR,
     *                XACMIN_KAPPAR,XACMAX_KAPPAR,R0CMIN_KAPPAR,
     *                R0CMAX_KAPPAR,XRCMIN_KAPPAR,XRCMAX_KAPPAR,
     *                XL_MIN_KAPPAR,XL_MAX_KAPPAR,XKSMIN_KAPPAR,
     *                XKSMAX_KAPPAR,A0SMIN_KAPPAR,A0SMAX_KAPPAR,
     *                XASMIN_KAPPAR,XASMAX_KAPPAR,R0SMIN_KAPPAR,
     *                R0SMAX_KAPPAR,XRSMIN_KAPPAR,XRSMAX_KAPPAR
      COMMON
     *       /VLIMIT/ VMISTR(1:NDPARS),
     *                VMIMIN(1:NDPARS),
     *                VMIMAX(1:NDPARS),
     *                VMISTP(1:NDPARS)
      COMMON
     *       /MEMORS/ MEMORY(1:NDPARS)
      COMMON
     *       /SEARCH/ IFTAKE(1:NDPARS)
      COMMON
     *       /ACTACT/ ACTION(1:NDPARS)
      COMMON
     *       /MINITE/ IDEFCN,ITECHI
      COMMON
     *       /PRFILE/ ISCREN,LOGFIL,LOGCHP,LOGCHN,LOGHIP,LOGHIN,
     *                LOGPRO,LOGNEU,LOGWRI,LSCREN,IRESUL,ICONVE,
     *                NOUTPT,LOGAUX,LOGENE,LOGRAD,LOGGRD,LOGMSH,
     *                                                   LOGBIS            
C
      DATA
     *        MAXMEM  / 150 /
C
C=======================================================================
C
C     This subroutine checks that all the values of the Hamiltonian 
C     parameters are not totally stupid (such as a negative radius, 
C     or a positive central well ... )
C  @@@ ???
C     LDSEAR: The actual number of the parameters searched
C


C @@@ DISCUSS THE USE OF THE PARAMETERS AND DELIMITERS




C=======================================================================
C
      IABORT=0
C
C=======================================================================
C 
      IF (IDEFCN.LE.1) THEN
C 
          IF (LOGWRI.GT.5) THEN
              WRITE(LOGFIL,
     *            '(27X,''IDEFCN.LE.1 => setting MEMORY to zero'')')
          END IF
C 
          DO I=1,NDPARS
             MEMORY(I)=0
          END DO
C 
      END IF     
C @@@    IS THE CHECK-UP BELOW TOTALLY USELESS?   ??? 
C=======================================================================
C     Below, we define an auxiliary vector that will allow
C     to keep all the parameters within given ranges
C=======================================================================
C
C     IACTIV=0
C
C     DO I=1,NDPARS/2
C
C        IF (IFTAKE(I).EQ.1) THEN
C
C            IACTIV=IACTIV+1
C
C            IF (ARGMNT(IACTIV).LT.0.0) THEN 
C 
C                IF (LOGWRI.GT.0) THEN
C                    WRITE(0,'(''Subroutine VERIFY: IACTIV='',I2,
C    *                         ''   OLD ARGMNT='',F12.6,
C    *                         ''   NEW ARGMNT='',F12.6)')
C    *               IACTIV,ARGMNT(IACTIV),0.5*STARTV(I) 
C                END IF
C             
C                ARGMNT(IACTIV)=0.5*STARTV(I)
C
C                ACTION(I)='#'
C
C                MEMORY(I)=MEMORY(I)+1                 
C                 
C            END IF 
C        END IF
C         
C     END DO
C
C=======================================================================
C
      IACTIV=0
C
      DO I=1,NDPARS
C
C=======================================================================
C        I=1 =>> Central potential depth
C=======================================================================
C
         IF ((I.EQ.1).AND.(IFTAKE(I).EQ.1)) THEN
C
             IACTIV=IACTIV+1
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).GT.V0CMAX) THEN
C
                 FACTOR=ABS(V0CMAX/ARGMNT(IACTIV))
C
                 IF (LOGWRI.GT.0) THEN                
C                 
                     WRITE(LOGFIL,'(18X,''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT='',E12.6,'' > '',F12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' Central Potential (P)'')')
     *                               I,IACTIV,ARGMNT(IACTIV),V0CMAX,
     *                                                FACTOR*V0CMAX
                 END IF
C
                 ARGMNT(IACTIV)=FACTOR*V0CMAX
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).LT.V0CMIN) THEN
C
                 FACTOR=ABS(ARGMNT(IACTIV))/(ABS(ARGMNT(IACTIV))+10)
C                 
                 IF (LOGWRI.GT.0) THEN
C                 
                     WRITE(LOGFIL,'(18X,''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT='',E12.6,'' < '',F12.6,
     *               '' NEWARG ='',F12.6,'' Central potential (P)'')')
     *                  I,IACTIV,ARGMNT(IACTIV),V0CMIN,
     *                                          FACTOR*V0CMIN
                 END IF
C
                 ARGMNT(IACTIV)=FACTOR*V0CMIN
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (MEMORY(I).GT.MAXMEM) THEN
C
                 WRITE(NOUTPT,'(/,80(''!''),/,''!'',78X,''!'',/,
     *                 ''!(P)Central potential depth out of range '',
     *                 ''more than MAXMEM='',I2,'' times. ABORT   !'',/,
     *                 ''!'',78X,''!'',/,80(''!''))')     MAXMEM
C
                 IF (LOGWRI.GT.0)
     *
     *           WRITE(LOGFIL,'(/,80(''!''),/,''!'',78X,''!'',/,
     *                 ''!(P)Central potential depth out of range '',
     *                 ''more than MAXMEM='',I2,'' times. ABORT   !'',/,
     *                 ''!'',78X,''!'',/,80(''!''))')     MAXMEM
C
                 IABORT=1
C
                 RETURN
C
             END IF
C
         END IF
C
C=======================================================================
C        I=4 =>> Spin-orbit strength coefficient Lambda
C=======================================================================
C
         IF (I.EQ.4.AND.IFTAKE(I).EQ.1) THEN
C
             IACTIV=IACTIV+1
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).GT.XL_MAX*UPFACT) THEN
C
                 FACTOR=XL_MAX/ARGMNT(IACTIV)             
C             
                 IF (LOGWRI.GT.0) THEN
C                 
                     WRITE(LOGFIL,'(18X,''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT='',E12.6,'' >  '',F12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' (Lambda s-o (P))'')')
     *                            I,IACTIV,ARGMNT(IACTIV),XL_MAX*UPFACT,
     *                                                    FACTOR*XL_MAX
                 END IF
C
                 ARGMNT(IACTIV)=FACTOR*XL_MAX
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).LT.XL_MIN*DWFACT) THEN
C
                 FACTOR=ABS(ARGMNT(IACTIV))/(ABS(ARGMNT(IACTIV))+5)
C                 
                 IF (LOGWRI.GT.0) THEN
C                 
                     WRITE(LOGFIL,'(18X,''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT ='',E12.6,'' < '',F12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' (Lambda s-o (P))'')')
     *                            I,IACTIV,ARGMNT(IACTIV),XL_MIN*DWFACT,
     *                                                    XL_MIN/FACTOR
                 END IF
C
                 ARGMNT(IACTIV)=XL_MIN/FACTOR
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (MEMORY(I).GT.MAXMEM) THEN
C
                 WRITE(NOUTPT,'(/,80(''!''),/,''!'',78X,''!'',/,
     *                 ''! (P) Lambda S-O parameter was out of range '',
     *                 ''more than MAXMEM='',I2,'' times. ABORT   !'',/,
     *                 ''!'',78X,''!'',/,80(''!''))')     MAXMEM
C
                 IF (LOGWRI.GT.0)
     *
     *           WRITE(LOGFIL,'(/,80(''!''),/,''!'',78X,''!'',/,
     *                 ''! (P) Lambda S-O parameter was out of range '',
     *                 ''more than MAXMEM='',I2,'' times. ABORT   !'',/,
     *                 ''!'',78X,''!'',/,80(''!''))')     MAXMEM
C
                 IABORT=1
C
                 RETURN
C
             END IF
C
         END IF
C
C=======================================================================
C        I=7 =>> Effective-mass strength coefficient Lambda
C=======================================================================
C @@@ DO WE HAVE EFFECTIVE MASS AND IN WHICH FORM?
         IF (I.EQ.7.AND.IFTAKE(I).EQ.1) THEN
C
             IACTIV=IACTIV+1
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).GT.XEFMAX*UPFACT) THEN
C
                 FACTOR=XEFMAX/ARGMNT(IACTIV)             
C
                 IF (LOGWRI.GT.0) THEN
C
                     WRITE(LOGFIL,'(18X,''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT ='',E12.6,'' > '',F12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' (Lambda eff-mass (P))'')')
     *                            I,IACTIV,ARGMNT(IACTIV),XEFMAX*UPFACT,
     *                                                    XEFMAX*FACTOR
                 END IF
C
                 ARGMNT(IACTIV)=XEFMAX*FACTOR
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).LT.XEFMIN*DWFACT) THEN
C
                 FACTOR=ABS(ARGMNT(IACTIV))/(ABS(ARGMNT(IACTIV))+5)
C
                 IF (LOGWRI.GT.0) THEN
C
                     WRITE(LOGFIL,'(18X,''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT='',E12.6,
     *                            ''<'',E12.6,'', NEWARG ='',F12.6,
     *                            '' (Lambda eff-mass (P))'')')
     *                            IACTIV,ARGMNT(IACTIV),XEFMIN*DWFACT,
     *                                                  XEFMIN/FACTOR
                 END IF
C
                 ARGMNT(IACTIV)=XEFMIN/FACTOR
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (MEMORY(I).GT.MAXMEM) THEN
C
                 WRITE(NOUTPT,'(/,80(''!''),/,''!'',78X,''!'',/,
     *                 ''!(P)Lambda EffMass paramtr was out of range '',
     *                 ''more than MAXMEM='',I2,'' times. ABORT   !'',/,
     *                 ''!'',78X,''!'',/,80(''!''))')     MAXMEM
C
                 IF (LOGWRI.GT.0)
     *
     *           WRITE(LOGFIL,'(/,80(''!''),/,''!'',78X,''!'',/,
     *                 ''!(P)Lambda EffMass paramtr was out of range '',
     *                 ''more than MAXMEM='',I2,'' times. ABORT   !'',/,
     *                 ''!'',78X,''!'',/,80(''!''))')     MAXMEM
C
                 IABORT=1
C
                 RETURN
C
             END IF
C
C=======================================================================
C
         END IF
C
C=======================================================================
C        I=2 =>> Central-potential radius
C=======================================================================
C
         IF ((I.EQ.2).AND.(IFTAKE(I).EQ.1)) THEN
C
             IACTIV=IACTIV+1
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).GT.R0CMAX*UPFACT) THEN
C
                 FACTOR=R0CMAX/ARGMNT(IACTIV)             
C
                 IF (LOGWRI.GT.0) THEN
C
                     WRITE(LOGFIL,'(18X,''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT ='',E12.6,
     *                            '' > '',E12.6, '' NEWARG ='',F12.6,
     *                            '' Central radius (P)'')')
     *                            I,IACTIV,ARGMNT(IACTIV),R0CMAX*UPFACT,
     *                                                    R0CMAX*FACTOR
                 END IF
C
                 ARGMNT(IACTIV)=R0CMAX*FACTOR
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).LT.R0CMIN*DWFACT) THEN
C
                 FACTOR=ABS(ARGMNT(IACTIV))/(ABS(ARGMNT(IACTIV))+0.1)                 
C
                 IF (LOGWRI.GT.0) THEN
C
                     WRITE(LOGFIL,'(18X,''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT ='',E12.6,
     *                            '' < '',E12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' Central radius (P)'')')
     *                            I,IACTIV,ARGMNT(IACTIV),R0CMIN*DWFACT,
     *                                                    R0CMIN/FACTOR
                 END IF
C
                 ARGMNT(IACTIV)=R0CMIN/FACTOR
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (MEMORY(I).GT.MAXMEM) THEN
C
                 WRITE(NOUTPT,'(/,80(''!''),/,''!'',78X,''!'',/,
     *                 ''! (P)Central potential r_0 was out of range '',
     *                 ''more than MAXMEM='',I2,'' times. ABORT   !'',/,
     *                 ''!'',78X,''!'',/,80(''!''))')     MAXMEM
C
                 IF (LOGWRI.GT.0)
     *
     *           WRITE(LOGFIL,'(/,80(''!''),/,''!'',78X,''!'',/,
     *                 ''! (P)Central potential r_0 was out of range '',
     *                 ''more than MAXMEM='',I2,'' times. ABORT   !'',/,
     *                 ''!'',78X,''!'',/,80(''!''))')     MAXMEM
C
                 IABORT=1
C
                 RETURN
C
             END IF
C
         END IF
C
C=======================================================================
C        I=5 =>> Spin-orbit radius 
C=======================================================================
C
         IF ((I.EQ.5).AND.(IFTAKE(I).EQ.1)) THEN
C
             IACTIV=IACTIV+1
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).GT.R0SMAX*UPFACT) THEN
C
                 FACTOR=R0SMAX/ARGMNT(IACTIV)             
C
                 IF (LOGWRI.GT.0) THEN
C
                     WRITE(LOGFIL,'(18X,''I='',I2,'' IACTIV ='',I3,
     *                              '' ARGMNT ='',E12.6,'' > '',F12.6,
     *                              '' NEWARG ='',F12.6,
     *                              '' S-O radius (P)'')')
     *                     I,IACTIV,ARGMNT(IACTIV),R0SMAX*UPFACT,
     *                                             R0SMAX*FACTOR
                 END IF
C
                 ARGMNT(IACTIV)=R0SMAX*FACTOR
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).LT.R0SMIN*DWFACT) THEN
C
                 FACTOR=ABS(ARGMNT(IACTIV))/(ABS(ARGMNT(IACTIV))+0.1)                 
C
                 IF (LOGWRI.GT.0) THEN
C
                     WRITE(LOGFIL,'(18X,''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT ='',E12.6,'' < '',F12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' S-O radius (P)'')')
     *                     I,IACTIV,ARGMNT(IACTIV),R0SMIN*DWFACT,
     *                                             R0SMIN/FACTOR
                 END IF
C
                 ARGMNT(IACTIV)=R0SMIN/FACTOR
C
                 MEMORY(I)=MEMORY(I)+1                 
C
                 ACTION(I)='#'
C
             END IF
C_______________________________________________________________________
C
             IF (MEMORY(I).GT.MAXMEM) THEN
C
                 WRITE(NOUTPT,'(/,80(''!''),/,''!'',78X,''!'',/,
     *                 ''! (P)Spin-Orbit radius r_0 was out of range '',
     *                 ''more than MAXMEM='',I2,'' times. ABORT   !'',/,
     *                 ''!'',78X,''!'',/,80(''!''))')     MAXMEM
C
                 IABORT=1
C
                 RETURN
C
             END IF
C
         END IF
C
C=======================================================================
C        I=8 =>> Effective mass radius
C=======================================================================
C
         IF ((I.EQ.8).AND.(IFTAKE(I).EQ.1)) THEN
C
             IACTIV=IACTIV+1
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).GT.REFMAX*UPFACT) THEN
C
                 FACTOR=REFMAX/ARGMNT(IACTIV)             
C
                 IF (LOGWRI.GT.0) THEN
C
                     WRITE(LOGFIL,'(''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT ='',E12.6,'' > '',F12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' Eff mass radius (P)'')')
     *                            I,IACTIV,ARGMNT(IACTIV),REFMAX*UPFACT,
     *                                                    REFMAX*FACTOR
                 END IF
C
                 ARGMNT(IACTIV)=REFMAX*FACTOR
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).LT.REFMIN*DWFACT) THEN
C
                 FACTOR=ABS(ARGMNT(IACTIV))/(ABS(ARGMNT(IACTIV))+0.1)                 
C
                 IF (LOGWRI.GT.0) THEN
C
                     WRITE(LOGFIL,'(''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT ='',E12.6,'' < '',F12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' Eff mass radius(P)'')')
     *                            I,IACTIV,ARGMNT(IACTIV),REFMIN*DWFACT,
     *                                                    REFMIN/FACTOR
                 END IF
C
                 ARGMNT(IACTIV)=REFMIN/FACTOR
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (MEMORY(I).GT.MAXMEM) THEN
C
                 WRITE(NOUTPT,'(/,80(''!''),/,''!'',78X,''!'',/,
     *                 ''! (P)Effective Mass radius was out of range '',
     *                 ''more than MAXMEM='',I2,'' times. ABORT   !'',/,
     *                 ''!'',78X,''!'',/,80(''!''))')     MAXMEM
C
                 IABORT=1
C
                 RETURN
C
             END IF
C
         END IF
C
C        Diffuseness parameters ...
C
         IF (I.EQ.3.AND.(IFTAKE(I).EQ.1)) THEN
C
             IACTIV=IACTIV+1
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).GT.A0CMAX*UPFACT) THEN
C
                 FACTOR=A0CMAX/ARGMNT(IACTIV)             
C
                 IF (LOGWRI.GT.0) THEN
C
                     WRITE(LOGFIL,'(''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT ='',E12.6,'' > '',F12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' Centr. diffuseness (P)'')')
     *                            I,IACTIV,ARGMNT(IACTIV),A0CMAX*UPFACT,
     *                                                    A0CMAX*FACTOR
                 END IF
C
                 ARGMNT(IACTIV)=A0CMAX*FACTOR
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).LT.A0CMIN*DWFACT) THEN
C
                 FACTOR=ABS(ARGMNT(IACTIV))/(ABS(ARGMNT(IACTIV))+0.1)                 
C
                 IF (LOGWRI.GT.0) THEN
C
                     WRITE(LOGFIL,'(''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT ='',E12.6,'' < '',F12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' Centr. diffuseness (P)'')')
     *                            I,IACTIV,ARGMNT(IACTIV),A0CMIN*DWFACT,
     *                                                    A0CMIN/FACTOR
                 END IF
C
                 ARGMNT(IACTIV)=A0CMIN/FACTOR
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (MEMORY(I).GT.MAXMEM) THEN
C
                 WRITE(NOUTPT,'(/,80(''!''),/,''!'',78X,''!'',/,
     *                 ''! (P)Central potential a_0 was out of range '',
     *                 ''more than MAXMEM='',I2,'' times. ABORT   !'',/,
     *                 ''!'',78X,''!'',/,80(''!''))')     MAXMEM
C
                 IABORT=1
C
                 RETURN
C
             END IF
C
         END IF
C
C=======================================================================
C
         IF ((I.EQ.6).AND.(IFTAKE(I).EQ.1)) THEN
C
             IACTIV=IACTIV+1
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).GT.A0SMAX*UPFACT) THEN
C
                 FACTOR=A0SMAX/ARGMNT(IACTIV)             
C
                 IF (LOGWRI.GT.0) THEN
C
                     WRITE(LOGFIL,'(''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT ='',E12.6,'' > '',F12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' S-O diffuseness (P)'')')
     *                            I,IACTIV,ARGMNT(IACTIV),A0SMAX*UPFACT,
     *                                                    A0SMAX*FACTOR
                 END IF
C
                 ARGMNT(IACTIV)=1.1*VMISTR(I)
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).LT.A0SMIN*DWFACT) THEN
C
                 FACTOR=ABS(ARGMNT(IACTIV))/(ABS(ARGMNT(IACTIV))+0.1)                 
C
                 IF (LOGWRI.GT.0) THEN
C
                     WRITE(LOGFIL,'(''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT ='',E12.6,'' < '',F12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' S-O diffuseness (P)'')')
     *                            I,IACTIV,ARGMNT(IACTIV),A0SMIN*DWFACT,
     *                                                    A0SMIN/FACTOR
                 END IF
C
                 ARGMNT(IACTIV)=A0SMIN/FACTOR
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (MEMORY(I).GT.MAXMEM) THEN
C
                 WRITE(NOUTPT,'(/,80(''!''),/,''!'',78X,''!'',/,
     *                 ''! (P)Spin-Orbit diffss a_0 was out of range '',
     *                 ''more than MAXMEM='',I2,'' times. ABORT   !'',/,
     *                 ''!'',78X,''!'',/,80(''!''))')     MAXMEM
C
                 IABORT=1
C
                 RETURN
C
             END IF
C
         END IF
C
C=======================================================================
C
         IF ((I.EQ.9).AND.(IFTAKE(I).EQ.1)) THEN
C
             IACTIV=IACTIV+1
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).GT.AEFMAX*UPFACT) THEN
C
                 FACTOR=AEFMAX/ARGMNT(IACTIV)             
C
                 IF (LOGWRI.GT.0) THEN
C
                     WRITE(LOGFIL,'(''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT ='',E12.6,'' > '',F12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' Meff diffuseness (P)'')')
     *                            I,IACTIV,ARGMNT(IACTIV),AEFMAX*UPFACT,
     *                                                    AEFMAX*FACTOR
                 END IF
C
                 ARGMNT(IACTIV)=AEFMAX*FACTOR
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).LT.AEFMIN*DWFACT) THEN
C
                 FACTOR=ABS(ARGMNT(IACTIV))/(ABS(ARGMNT(IACTIV))+0.1)                 
C
                 IF (LOGWRI.GT.0) THEN
C
                     WRITE(LOGFIL,'(''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT ='',E12.6,'' < '',F12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' Meff diffuseness (P)'')')
     *                            I,IACTIV,ARGMNT(IACTIV),AEFMIN*DWFACT,
     *                                                    AEFMIN/FACTOR
                 END IF
C
                 ARGMNT(IACTIV)=AEFMIN/FACTOR
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (MEMORY(I).GT.MAXMEM) THEN
C
                 WRITE(NOUTPT,'(/,80(''!''),/,''!'',78X,''!'',/,
     *                 ''!  (P)Eff-Mass diffuseness was out of range '',
     *                 ''more than MAXMEM='',I2,'' times. ABORT   !'',/,
     *                 ''!'',78X,''!'',/,80(''!''))')     MAXMEM
C
                 IABORT=1
C
                 RETURN
C
             END IF
C
         END IF
C
C=======================================================================
C
         IF ((I.EQ.10).AND.(IFTAKE(I).EQ.1)) THEN
C
             IACTIV=IACTIV+1
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).GT.CouMAX*UPFACT) THEN
C
                 FACTOR=CouMAX/ARGMNT(IACTIV)             
C
                 IF (LOGWRI.GT.0) THEN
C
                     WRITE(LOGFIL,'(''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT ='',E12.6,'' > '',F12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' Coulomb radius'')')
     *                            I,IACTIV,ARGMNT(IACTIV),CouMAX*UPFACT,
     *                                                    CouMAX*FACTOR
                 END IF
C
                 ARGMNT(IACTIV)=CouMAX*FACTOR
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).LT.CouMIN*DWFACT) THEN
C
                 FACTOR=ABS(ARGMNT(IACTIV))/(ABS(ARGMNT(IACTIV))+0.1)                 
C
                 IF (LOGWRI.GT.0) THEN
C
                     WRITE(LOGFIL,'(''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT ='',E12.6,'' < '',F12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' Coulomb radius'')')
     *                            I,IACTIV,ARGMNT(IACTIV),CouMIN*DWFACT,
     *                                                    CouMIN/FACTOR
                 END IF
C
                 ARGMNT(IACTIV)=CouMIN/FACTOR
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (MEMORY(I).GT.MAXMEM) THEN
C
                 WRITE(NOUTPT,'(/,80(''!''),/,''!'',78X,''!'',/,
     *                 ''!   Coulomb potential r_0  was out of range '',
     *                 ''more than MAXMEM='',I2,'' times. ABORT   !'',/,
     *                 ''!'',78X,''!'',/,80(''!''))')     MAXMEM
C
                 IABORT=1
C
                 RETURN
C
             END IF
C
         END IF
C
C=======================================================================
C        I=21 =>> Central potential depth
C=======================================================================
C
         IF ((I.EQ.21).AND.(IFTAKE(I).EQ.1)) THEN
C
             IACTIV=IACTIV+1
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).GT.V0CMAX) THEN
C
                 FACTOR=ABS(V0CMAX/ARGMNT(IACTIV))
C
                 IF (LOGWRI.GT.0) THEN                
C                 
                     WRITE(LOGFIL,'(''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT='',E12.6,'' > '',F12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' Central Potential (N)'')')
     *                            I,IACTIV,ARGMNT(IACTIV),V0CMAX,
     *                                             FACTOR*V0CMAX
                 END IF
C
                 ARGMNT(IACTIV)=FACTOR*V0CMAX
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).LT.V0CMIN) THEN
C
                 FACTOR=ABS(ARGMNT(IACTIV))/(ABS(ARGMNT(IACTIV))+10)
C                 
                 IF (LOGWRI.GT.0) THEN
C                 
                     WRITE(LOGFIL,'(''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT='',E12.6,'' <  '',F12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' Central potential (N)'')')
     *                            I,IACTIV,ARGMNT(IACTIV),V0CMIN,
     *                                             FACTOR*V0CMIN
                 END IF
C
                 ARGMNT(IACTIV)=FACTOR*V0CMIN
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (MEMORY(I).GT.MAXMEM) THEN
C
                 WRITE(NOUTPT,'(/,80(''!''),/,''!'',78X,''!'',/,
     *                 ''!(N)Cental potential depth was out of range '',
     *                 ''more than MAXMEM='',I2,'' times. ABORT   !'',/,
     *                 ''!'',78X,''!'',/,80(''!''))')     MAXMEM
C
                 IABORT=1
C
                 RETURN
C
             END IF
C
         END IF
C
C=======================================================================
C        I=24 =>> Spin-orbit     strength coefficient Lambda
C=======================================================================
C
         IF (I.EQ.24.AND.IFTAKE(I).EQ.1) THEN
C
             IACTIV=IACTIV+1
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).GT.XL_MAX*UPFACT) THEN
C
                 FACTOR=XL_MAX/ARGMNT(IACTIV)             
C             
                 IF (LOGWRI.GT.0) THEN
C                 
                     WRITE(LOGFIL,'(''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT='',E12.6,'' >  '',F12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' (Lambda s-o (N))'')')
     *                            I,IACTIV,ARGMNT(IACTIV),XL_MAX*UPFACT,
     *                                                    FACTOR*XL_MAX
                 END IF
C
                 ARGMNT(IACTIV)=FACTOR*XL_MAX
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).LT.XL_MIN*DWFACT) THEN
C
                 FACTOR=ABS(ARGMNT(IACTIV))/(ABS(ARGMNT(IACTIV))+5)                 
C                 
                 IF (LOGWRI.GT.0) THEN
C                 
                     WRITE(LOGFIL,'(''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT ='',E12.6,'' < '',F12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' (Lambda s-o (N))'')')
     *                            I,IACTIV,ARGMNT(IACTIV),XL_MIN*DWFACT,
     *                                                    XL_MIN/FACTOR
                 END IF
C
                 ARGMNT(IACTIV)=XL_MIN/FACTOR
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (MEMORY(I).GT.MAXMEM) THEN
C
                 WRITE(NOUTPT,'(/,80(''!''),/,''!'',78X,''!'',/,
     *                 ''! (N) Lambda S-O parameter was out of range '',
     *                 ''more than MAXMEM='',I2,'' times. ABORT   !'',/,
     *                 ''!'',78X,''!'',/,80(''!''))')     MAXMEM
C
                 IABORT=1
C
                 RETURN
C
             END IF
C
         END IF
C
C=======================================================================
C        I=27 =>> Effective-mass strength coefficient Lambda
C=======================================================================
C
         IF (I.EQ.27.AND.IFTAKE(I).EQ.1) THEN
C
             IACTIV=IACTIV+1
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).GT.XEFMAX*UPFACT) THEN
C
                 FACTOR=XEFMAX/ARGMNT(IACTIV)             
C
                 IF (LOGWRI.GT.0) THEN
                     WRITE(LOGFIL,'(''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT ='',E12.6,'' > '',F12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' (Lambda eff-mass (N))'')')
     *                            I,IACTIV,ARGMNT(IACTIV),XEFMAX*UPFACT,
     *                                                    XEFMAX*FACTOR
                 END IF
C
                 ARGMNT(IACTIV)=XEFMAX*FACTOR
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).LT.XEFMIN*DWFACT) THEN
C
                 FACTOR=ABS(ARGMNT(IACTIV))/(ABS(ARGMNT(IACTIV))+5)                 
C
                 IF (LOGWRI.GT.0) THEN
C
                     WRITE(LOGFIL,'(''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT='',E12.6,
     *                            ''<'',E12.6,'', NEWARG ='',F12.6,
     *                            '' (Lambda eff-mass (N))'')')
     *                            IACTIV,ARGMNT(IACTIV),XEFMIN*DWFACT,
     *                                                  XEFMIN/FACTOR
                 END IF
C
                 ARGMNT(IACTIV)=XEFMIN/FACTOR
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (MEMORY(I).GT.MAXMEM) THEN
C
                 WRITE(NOUTPT,'(/,80(''!''),/,''!'',78X,''!'',/,
     *                 ''!(N)Lambda EffMass paramtr was out of range '',
     *                 ''more than MAXMEM='',I2,'' times. ABORT   !'',/,
     *                 ''!'',78X,''!'',/,80(''!''))')     MAXMEM
C
                 IABORT=1
C
                 RETURN
C
             END IF
C
C=======================================================================
C
         END IF
C
C=======================================================================
C        I=22 =>> Central-potential radius
C=======================================================================
C
         IF ((I.EQ.22).AND.(IFTAKE(I).EQ.1)) THEN
C
             IACTIV=IACTIV+1
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).GT.R0CMAX*UPFACT) THEN
C
                 FACTOR=R0CMAX/ARGMNT(IACTIV)             
C
                 IF (LOGWRI.GT.0) THEN
C
                     WRITE(LOGFIL,'(''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT ='',E12.6,
     *                            '' > '',E12.6, '' NEWARG ='',F12.6,
     *                            '' Central radius (N)'')')
     *                            I,IACTIV,ARGMNT(IACTIV),R0CMAX*UPFACT,
     *                                                    R0CMAX*FACTOR
                 END IF
C
                 ARGMNT(IACTIV)=R0CMAX*FACTOR
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).LT.R0CMIN*DWFACT) THEN
C
                 FACTOR=ABS(ARGMNT(IACTIV))/(ABS(ARGMNT(IACTIV))+0.1)                 
C
                 IF (LOGWRI.GT.0) THEN
C
                     WRITE(LOGFIL,'(''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT ='',E12.6,
     *                            '' < '',E12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' Central radius (N)'')')
     *                            I,IACTIV,ARGMNT(IACTIV),R0CMIN*DWFACT,
     *                                                    R0CMIN/FACTOR
                 END IF
C
                 ARGMNT(IACTIV)=R0CMIN/FACTOR
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (MEMORY(I).GT.MAXMEM) THEN
C
                 WRITE(NOUTPT,'(/,80(''!''),/,''!'',78X,''!'',/,
     *                 ''! (N)Central potential r_0 was out of range '',
     *                 ''more than MAXMEM='',I2,'' times. ABORT   !'',/,
     *                 ''!'',78X,''!'',/,80(''!''))')     MAXMEM
C
                 IABORT=1
C
                 RETURN
C
             END IF
C
         END IF
C
C=======================================================================
C        I=25 =>> Spin-orbit radius 
C=======================================================================
C
         IF ((I.EQ.25).AND.(IFTAKE(I).EQ.1)) THEN
C
             IACTIV=IACTIV+1
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).GT.R0SMAX*UPFACT) THEN
C
                 FACTOR=R0SMAX/ARGMNT(IACTIV)             
C
                 IF (LOGWRI.GT.0) THEN
C
                     WRITE(LOGFIL,'(''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT ='',E12.6,'' > '',F12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' S-O radius (N)'')')
     *                            I,IACTIV,ARGMNT(IACTIV),R0SMAX*UPFACT,
     *                                                    R0SMAX*FACTOR
                 END IF
C
                 ARGMNT(IACTIV)=R0SMAX*FACTOR
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).LT.R0SMIN*DWFACT) THEN
C
                 FACTOR=ABS(ARGMNT(IACTIV))/(ABS(ARGMNT(IACTIV))+0.1)                 
C
                 IF (LOGWRI.GT.0) THEN
C
                     WRITE(LOGFIL,'(''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT ='',E12.6,'' < '',F12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' S-O radius (N)'')')
     *                            I,IACTIV,ARGMNT(IACTIV),R0SMIN*DWFACT,
     *                                                    R0SMIN/FACTOR
                 END IF
C
                 ARGMNT(IACTIV)=R0SMIN/FACTOR
C
                 MEMORY(I)=MEMORY(I)+1                 
C
                 ACTION(I)='#'
C
             END IF
C_______________________________________________________________________
C
             IF (MEMORY(I).GT.MAXMEM) THEN
C
                 WRITE(NOUTPT,'(/,80(''!''),/,''!'',78X,''!'',/,
     *                 ''! (N)Spin-Orbit radius r_0 was out of range '',
     *                 ''more than MAXMEM='',I2,'' times. ABORT   !'',/,
     *                 ''!'',78X,''!'',/,80(''!''))')     MAXMEM
C
                 IABORT=1
C
                 RETURN
C
             END IF
C
         END IF
C
C=======================================================================
C        I=28 =>> Effective mass radius
C=======================================================================
C
         IF ((I.EQ.28).AND.(IFTAKE(I).EQ.1)) THEN
C
             IACTIV=IACTIV+1
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).GT.REFMAX*UPFACT) THEN
C
                 FACTOR=REFMAX/ARGMNT(IACTIV)             
C
                 IF (LOGWRI.GT.0) THEN
C
                     WRITE(LOGFIL,'(''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT ='',E12.6,'' > '',F12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' Eff mass radius (N)'')')
     *                            I,IACTIV,ARGMNT(IACTIV),REFMAX*UPFACT,
     *                                                    REFMAX*FACTOR
                 END IF
C
                 ARGMNT(IACTIV)=REFMAX*FACTOR
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).LT.REFMIN*DWFACT) THEN
C
                 FACTOR=ABS(ARGMNT(IACTIV))/(ABS(ARGMNT(IACTIV))+0.1)                 
C
                 IF (LOGWRI.GT.0) THEN
C
                     WRITE(LOGFIL,'(''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT ='',E12.6,'' < '',F12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' Eff mass radius (N)'')')
     *                            I,IACTIV,ARGMNT(IACTIV),REFMIN*DWFACT,
     *                                                    REFMIN/FACTOR
                 END IF
C
                 ARGMNT(IACTIV)=REFMIN/FACTOR
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (MEMORY(I).GT.MAXMEM) THEN
C
                 WRITE(NOUTPT,'(/,80(''!''),/,''!'',78X,''!'',/,
     *                 ''! (N)Effective Mass radius was out of range '',
     *                 ''more than MAXMEM='',I2,'' times. ABORT   !'',/,
     *                 ''!'',78X,''!'',/,80(''!''))')     MAXMEM
C
                 IABORT=1
C
                 RETURN
C
             END IF
C
         END IF
C
C        Diffuseness parameters ...
C
         IF (I.EQ.23.AND.(IFTAKE(I).EQ.1)) THEN
C
             IACTIV=IACTIV+1
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).GT.A0CMAX*UPFACT) THEN
C
                 FACTOR=A0CMAX/ARGMNT(IACTIV)             
C
                 IF (LOGWRI.GT.0) THEN
C
                     WRITE(LOGFIL,'(''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT ='',E12.6,'' > '',F12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' Centr. diffuseness (N)'')')
     *                            I,IACTIV,ARGMNT(IACTIV),A0CMAX*UPFACT,
     *                                                    A0CMAX*FACTOR
                 END IF
C
                 ARGMNT(IACTIV)=A0CMAX*FACTOR
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).LT.A0CMIN*DWFACT) THEN
C
                 FACTOR=ABS(ARGMNT(IACTIV))/(ABS(ARGMNT(IACTIV))+0.1)                 
C
                 IF (LOGWRI.GT.0) THEN
C
                     WRITE(LOGFIL,'(''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT ='',E12.6,'' < '',F12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' Centr. diffuseness (N)'')')
     *                            I,IACTIV,ARGMNT(IACTIV),A0CMIN*DWFACT,
     *                                                    A0CMIN/FACTOR
                 END IF
C
                 ARGMNT(IACTIV)=A0CMIN/FACTOR
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (MEMORY(I).GT.MAXMEM) THEN
C
                 WRITE(NOUTPT,'(/,80(''!''),/,''!'',78X,''!'',/,
     *                 ''! (N)Central potential a_0 was out of range '',
     *                 ''more than MAXMEM='',I2,'' times. ABORT   !'',/,
     *                 ''!'',78X,''!'',/,80(''!''))')     MAXMEM
C
                 IABORT=1
C
                 RETURN
C
             END IF
C
         END IF
C
C=======================================================================
C
         IF ((I.EQ.26).AND.(IFTAKE(I).EQ.1)) THEN
C
             IACTIV=IACTIV+1
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).GT.A0SMAX*UPFACT) THEN
C
                 FACTOR=A0SMAX/ARGMNT(IACTIV)             
C
                 IF (LOGWRI.GT.0) THEN
C
                     WRITE(LOGFIL,'(''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT ='',E12.6,'' > '',F12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' S-O diffuseness (N)'')')
     *                            I,IACTIV,ARGMNT(IACTIV),A0SMAX*UPFACT,
     *                                                    A0SMAX*FACTOR
                 END IF
C
                 ARGMNT(IACTIV)=1.1*VMISTR(I)
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).LT.A0SMIN*DWFACT) THEN
C
                 FACTOR=ABS(ARGMNT(IACTIV))/(ABS(ARGMNT(IACTIV))+0.1)                 
C
                 IF (LOGWRI.GT.0) THEN
C
                     WRITE(LOGFIL,'(''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT ='',E12.6,'' < '',F12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' S-O diffuseness (N)'')')
     *                            I,IACTIV,ARGMNT(IACTIV),A0SMIN*DWFACT,
     *                                                    A0SMIN/FACTOR
                 END IF
C
                 ARGMNT(IACTIV)=A0SMIN/FACTOR
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (MEMORY(I).GT.MAXMEM) THEN
C
                 WRITE(NOUTPT,'(/,80(''!''),/,''!'',78X,''!'',/,
     *                 ''! (N)Spin-Orbit diffss a_0 was out of range '',
     *                 ''more than MAXMEM='',I2,'' times. ABORT   !'',/,
     *                 ''!'',78X,''!'',/,80(''!''))')     MAXMEM
C
                 IABORT=1
C
                 RETURN
C
             END IF
C
         END IF
C
C=======================================================================
C
         IF ((I.EQ.29).AND.(IFTAKE(I).EQ.1)) THEN
C
             IACTIV=IACTIV+1
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).GT.AEFMAX*UPFACT) THEN
C
                 FACTOR=AEFMAX/ARGMNT(IACTIV)             
C
                 IF (LOGWRI.GT.0) THEN
C
                     WRITE(LOGFIL,'(''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT ='',E12.6,'' > '',F12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' Meff diffuseness (N)'')')
     *                            I,IACTIV,ARGMNT(IACTIV),AEFMAX*UPFACT,
     *                                                    AEFMAX*FACTOR
                 END IF
C
                 ARGMNT(IACTIV)=AEFMAX*FACTOR
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).LT.AEFMIN*DWFACT) THEN
C
                 FACTOR=ABS(ARGMNT(IACTIV))/(ABS(ARGMNT(IACTIV))+0.1)                 
C
                 IF (LOGWRI.GT.0) THEN
C
                     WRITE(LOGFIL,'(''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT ='',E12.6,'' < '',F12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' Meff diffuseness (N)'')')
     *                            I,IACTIV,ARGMNT(IACTIV),AEFMIN*DWFACT,
     *                                                    AEFMIN/FACTOR
                 END IF
C
                 ARGMNT(IACTIV)=AEFMIN/FACTOR
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (MEMORY(I).GT.MAXMEM) THEN
C
                 WRITE(NOUTPT,'(/,80(''!''),/,''!'',78X,''!'',/,
     *                 ''!  (N)Eff-Mass diffuseness was out of range '',
     *                 ''more than MAXMEM='',I2,'' times. ABORT   !'',/,
     *                 ''!'',78X,''!'',/,80(''!''))')     MAXMEM
C
                 IABORT=1
C
                 RETURN
C
             END IF
C
         END IF
C
C=======================================================================
C        I=39 =>> Spin-orbit with density strength coefficient LambdaPP
C=======================================================================
C
         IF (I.EQ.39.AND.IFTAKE(I).EQ.1) THEN
C
             IACTIV=IACTIV+1
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).GT.XPPMAX*UPFACT) THEN
C
                 FACTOR=XPPMAX/ARGMNT(IACTIV)             
C             
                 IF (LOGWRI.GT.0) THEN
C                 
                     WRITE(LOGFIL,'(''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT='',E12.6,'' >  '',F12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' (Lambda PP)'')')
     *                            I,IACTIV,ARGMNT(IACTIV),XPPMAX*UPFACT,
     *                                                    FACTOR*XPPMAX
                 END IF
C
                 ARGMNT(IACTIV)=FACTOR*XPPMAX
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).LT.XPPMIN*DWFACT) THEN
C
                 FACTOR=ABS(ARGMNT(IACTIV))/(ABS(ARGMNT(IACTIV))+5)                 
C                 
                 IF (LOGWRI.GT.0) THEN
C                 
                     WRITE(LOGFIL,'(''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT ='',E12.6,'' < '',F12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' (Lambda PP)'')')
     *                            I,IACTIV,ARGMNT(IACTIV),XPPMIN*DWFACT,
     *                                                    XPPMIN/FACTOR
                 END IF
C
                 IF (ARGMNT(IACTIV).LT.0.0) THEN
                     ARGMNT(IACTIV)=XPPMIN*FACTOR
                 ELSE
                     ARGMNT(IACTIV)=XPPMIN/FACTOR
                 END IF
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (MEMORY(I).GT.MAXMEM) THEN
C
                 WRITE(NOUTPT,'(/,80(''!''),/,''!'',78X,''!'',/,
     *                 ''!     Lambda PP  parameter was out of range '',
     *                 ''more than MAXMEM='',I2,'' times. ABORT   !'',/,
     *                 ''!'',78X,''!'',/,80(''!''))')     MAXMEM
C
                 IABORT=1
C
                 RETURN
C
             END IF
C
         END IF
C
C=======================================================================
C        I=40 =>> Spin-orbit with density strength coefficient LambdaPN
C=======================================================================
C
         IF (I.EQ.40.AND.IFTAKE(I).EQ.1) THEN
C
             IACTIV=IACTIV+1
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).GT.XPNMAX*UPFACT) THEN
C
                 FACTOR=XPNMAX/ARGMNT(IACTIV)             
C             
                 IF (LOGWRI.GT.0) THEN
C                 
                     WRITE(LOGFIL,'(''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT='',E12.6,'' >  '',F12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' (Lambda PN)'')')
     *                            I,IACTIV,ARGMNT(IACTIV),XPNMAX*UPFACT,
     *                                                    FACTOR*XPNMAX
                 END IF
C
                 ARGMNT(IACTIV)=FACTOR*XPNMAX
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).LT.XPNMIN*DWFACT) THEN
C
                 FACTOR=ABS(ARGMNT(IACTIV))/(ABS(ARGMNT(IACTIV))+5)                 
C                 
                 IF (LOGWRI.GT.0) THEN
C                 
                     WRITE(LOGFIL,'(''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT ='',E12.6,'' < '',F12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' (Lambda PN)'')')
     *                            I,IACTIV,ARGMNT(IACTIV),XPNMIN*DWFACT,
     *                                                    XPNMIN/FACTOR
                 END IF
C
                 IF (ARGMNT(IACTIV).LT.0.0) THEN
                     ARGMNT(IACTIV)=XPNMIN*FACTOR
                 ELSE
                     ARGMNT(IACTIV)=XPNMIN/FACTOR
                 END IF
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (MEMORY(I).GT.MAXMEM) THEN
C
                 WRITE(NOUTPT,'(/,80(''!''),/,''!'',78X,''!'',/,
     *                 ''!     Lambda PN  parameter was out of range '',
     *                 ''more than MAXMEM='',I2,'' times. ABORT   !'',/,
     *                 ''!'',78X,''!'',/,80(''!''))')     MAXMEM
C
                 IABORT=1
C
                 RETURN
C
             END IF
C
         END IF
C
C=======================================================================
C        I=41 =>> Spin-orbit with density strength coefficient LambdaNP
C=======================================================================
C
         IF (I.EQ.41.AND.IFTAKE(I).EQ.1) THEN
C
             IACTIV=IACTIV+1
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).GT.XNPMAX*UPFACT) THEN
C
                 FACTOR=XNPMAX/ARGMNT(IACTIV)             
C             
                 IF (LOGWRI.GT.0) THEN
C                 
                     WRITE(LOGFIL,'(''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT='',E12.6,'' >  '',F12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' (Lambda NP)'')')
     *                            I,IACTIV,ARGMNT(IACTIV),XNPMAX*UPFACT,
     *                                                    FACTOR*XNPMAX
                 END IF
C
                 ARGMNT(IACTIV)=FACTOR*XNPMAX
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).LT.XNPMIN*DWFACT) THEN
C
                 FACTOR=ABS(ARGMNT(IACTIV))/(ABS(ARGMNT(IACTIV))+5)                 
C                 
                 IF (LOGWRI.GT.0) THEN
C                 
                     WRITE(LOGFIL,'(''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT ='',E12.6,'' < '',F12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' (Lambda NP)'')')
     *                            I,IACTIV,ARGMNT(IACTIV),XNPMIN*DWFACT,
     *                                                    XNPMIN/FACTOR
                 END IF
C
                 IF (ARGMNT(IACTIV).LT.0.0) THEN
                     ARGMNT(IACTIV)=XNPMIN*FACTOR
                 ELSE
                     ARGMNT(IACTIV)=XNPMIN/FACTOR
                 END IF
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (MEMORY(I).GT.MAXMEM) THEN
C
                 WRITE(NOUTPT,'(/,80(''!''),/,''!'',78X,''!'',/,
     *                 ''!     Lambda NP  parameter was out of range '',
     *                 ''more than MAXMEM='',I2,'' times. ABORT   !'',/,
     *                 ''!'',78X,''!'',/,80(''!''))')     MAXMEM
C
                 IABORT=1
C
                 RETURN
C
             END IF
C
         END IF
C
C=======================================================================
C        I=42 =>> Spin-orbit with density strength coefficient LambdaNN
C=======================================================================
C
         IF (I.EQ.42.AND.IFTAKE(I).EQ.1) THEN
C
             IACTIV=IACTIV+1
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).GT.XNNMAX*UPFACT) THEN
C
                 FACTOR=XNNMAX/ARGMNT(IACTIV)             
C             
                 IF (LOGWRI.GT.0) THEN
C                 
                     WRITE(LOGFIL,'(''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT='',E12.6,'' >  '',F12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' (Lambda NN)'')')
     *                            I,IACTIV,ARGMNT(IACTIV),XNNMAX*UPFACT,
     *                                                    FACTOR*XNNMAX
                 END IF
C
                 ARGMNT(IACTIV)=FACTOR*XNNMAX
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).LT.XNNMIN*DWFACT) THEN
C
                 FACTOR=ABS(ARGMNT(IACTIV))/(ABS(ARGMNT(IACTIV))+5)                 
C                 
                 IF (LOGWRI.GT.0) THEN
C                 
                     WRITE(LOGFIL,'(''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT ='',E12.6,'' < '',F12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' (Lambda NN)'')')
     *                            I,IACTIV,ARGMNT(IACTIV),XNNMIN*DWFACT,
     *                                                    XNNMIN/FACTOR
                 END IF
C
                 IF (ARGMNT(IACTIV).LT.0.0) THEN
                     ARGMNT(IACTIV)=XNNMIN*FACTOR
                 ELSE
                     ARGMNT(IACTIV)=XNNMIN/FACTOR
                 END IF
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (MEMORY(I).GT.MAXMEM) THEN
C
                 WRITE(NOUTPT,'(/,80(''!''),/,''!'',78X,''!'',/,
     *                 ''!     Lambda NN  parameter was out of range '',
     *                 ''more than MAXMEM='',I2,'' times. ABORT   !'',/,
     *                 ''!'',78X,''!'',/,80(''!''))')     MAXMEM
C
                 IABORT=1
C
                 RETURN
C
             END IF
C
         END IF
C
C=======================================================================
C        I=43 =>> Tensor part strength coefficient LambdaPP
C=======================================================================
C
         IF (I.EQ.43.AND.IFTAKE(I).EQ.1) THEN
C
             IACTIV=IACTIV+1
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).GT.YPPMAX*UPFACT) THEN
C
                 FACTOR=YPPMAX/ARGMNT(IACTIV)             
C             
                 IF (LOGWRI.GT.0) THEN
C                 
                     WRITE(LOGFIL,'(''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT='',E12.6,'' >  '',F12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' (Ylambd PP)'')')
     *                            I,IACTIV,ARGMNT(IACTIV),YPPMAX*UPFACT,
     *                                                    FACTOR*YPPMAX
                 END IF
C
                 ARGMNT(IACTIV)=FACTOR*YPPMAX
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).LT.YPPMIN*DWFACT) THEN
C
                 FACTOR=ABS(ARGMNT(IACTIV))/(ABS(ARGMNT(IACTIV))+5)                 
C                 
                 IF (LOGWRI.GT.0) THEN
C                 
                     WRITE(LOGFIL,'(''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT ='',E12.6,'' < '',F12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' (Ylambd PP)'')')
     *                            I,IACTIV,ARGMNT(IACTIV),YPPMIN*DWFACT,
     *                                                    YPPMIN/FACTOR
                 END IF
C
                 IF (ARGMNT(IACTIV).LT.0.0) THEN
                     ARGMNT(IACTIV)=YPPMIN*FACTOR
                 ELSE
                     ARGMNT(IACTIV)=YPPMIN/FACTOR
                 END IF
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (MEMORY(I).GT.MAXMEM) THEN
C
                 WRITE(NOUTPT,'(/,80(''!''),/,''!'',78X,''!'',/,
     *                 ''!     Ylambd PP  parameter was out of range '',
     *                 ''more than MAXMEM='',I2,'' times. ABORT   !'',/,
     *                 ''!'',78X,''!'',/,80(''!''))')     MAXMEM
C
                 IABORT=1
C
                 RETURN
C
             END IF
C
         END IF
C
C=======================================================================
C        I=44 =>> Tensor part strength coefficient LambdaPN
C=======================================================================
C
         IF (I.EQ.44.AND.IFTAKE(I).EQ.1) THEN
C
             IACTIV=IACTIV+1
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).GT.YPNMAX*UPFACT) THEN
C
                 FACTOR=YPNMAX/ARGMNT(IACTIV)             
C             
                 IF (LOGWRI.GT.0) THEN
C                 
                     WRITE(LOGFIL,'(''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT='',E12.6,'' >  '',F12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' (Ylambd PN)'')')
     *                            I,IACTIV,ARGMNT(IACTIV),YPNMAX*UPFACT,
     *                                                    FACTOR*YPNMAX
                 END IF
C
                 ARGMNT(IACTIV)=FACTOR*YPNMAX
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).LT.YPNMIN*DWFACT) THEN
C
                 FACTOR=ABS(ARGMNT(IACTIV))/(ABS(ARGMNT(IACTIV))+5)                 
C                 
                 IF (LOGWRI.GT.0) THEN
C                 
                     WRITE(LOGFIL,'(''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT ='',E12.6,'' < '',F12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' (Lambda PN)'')')
     *                            I,IACTIV,ARGMNT(IACTIV),YPNMIN*DWFACT,
     *                                                    YPNMIN/FACTOR
                 END IF
C
                 IF (ARGMNT(IACTIV).LT.0.0) THEN
                     ARGMNT(IACTIV)=YPNMIN*FACTOR
                 ELSE
                     ARGMNT(IACTIV)=YPNMIN/FACTOR
                 END IF
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (MEMORY(I).GT.MAXMEM) THEN
C
                 WRITE(NOUTPT,'(/,80(''!''),/,''!'',78X,''!'',/,
     *                 ''!     Ylambd PN  parameter was out of range '',
     *                 ''more than MAXMEM='',I2,'' times. ABORT   !'',/,
     *                 ''!'',78X,''!'',/,80(''!''))')     MAXMEM
C
                 IABORT=1
C
                 RETURN
C
             END IF
C
         END IF
C
C=======================================================================
C        I=45 =>> Tensor part strength coefficient LambdaNP
C=======================================================================
C
         IF (I.EQ.45.AND.IFTAKE(I).EQ.1) THEN
C
             IACTIV=IACTIV+1
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).GT.YNPMAX*UPFACT) THEN
C
                 FACTOR=YNPMAX/ARGMNT(IACTIV)             
C             
                 IF (LOGWRI.GT.0) THEN
C                 
                     WRITE(LOGFIL,'(''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT='',E12.6,'' >  '',F12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' (Ylambd NP)'')')
     *                            I,IACTIV,ARGMNT(IACTIV),YNPMAX*UPFACT,
     *                                                    FACTOR*YNPMAX
                 END IF
C
                 ARGMNT(IACTIV)=FACTOR*YNPMAX
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).LT.YNPMIN*DWFACT) THEN
C
                 FACTOR=ABS(ARGMNT(IACTIV))/(ABS(ARGMNT(IACTIV))+5)                 
C                 
                 IF (LOGWRI.GT.0) THEN
C                 
                     WRITE(LOGFIL,'(''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT ='',E12.6,'' < '',F12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' (Ylambd NP)'')')
     *                            I,IACTIV,ARGMNT(IACTIV),YNPMIN*DWFACT,
     *                                                    YNPMIN/FACTOR
                 END IF
C
                 IF (ARGMNT(IACTIV).LT.0.0) THEN
                     ARGMNT(IACTIV)=YNPMIN*FACTOR
                 ELSE
                     ARGMNT(IACTIV)=YNPMIN/FACTOR
                 END IF
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (MEMORY(I).GT.MAXMEM) THEN
C
                 WRITE(NOUTPT,'(/,80(''!''),/,''!'',78X,''!'',/,
     *                 ''!     Ylambd NP  parameter was out of range '',
     *                 ''more than MAXMEM='',I2,'' times. ABORT   !'',/,
     *                 ''!'',78X,''!'',/,80(''!''))')     MAXMEM
C
                 IABORT=1
C
                 RETURN
C
             END IF
C
         END IF
C
C=======================================================================
C        I=46 =>> Tensor part strength coefficient LambdaNN
C=======================================================================
C
         IF (I.EQ.46.AND.IFTAKE(I).EQ.1) THEN
C
             IACTIV=IACTIV+1
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).GT.YNNMAX*UPFACT) THEN
C
                 FACTOR=YNNMAX/ARGMNT(IACTIV)             
C             
                 IF (LOGWRI.GT.0) THEN
C                 
                     WRITE(LOGFIL,'(''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT='',E12.6,'' >  '',F12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' (Ylambd NN)'')')
     *                            I,IACTIV,ARGMNT(IACTIV),YNNMAX*UPFACT,
     *                                                    FACTOR*YNNMAX
                 END IF
C
                 ARGMNT(IACTIV)=FACTOR*YNNMAX
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).LT.YNNMIN*DWFACT) THEN
C
                 FACTOR=ABS(ARGMNT(IACTIV))/(ABS(ARGMNT(IACTIV))+5)                 
C                 
                 IF (LOGWRI.GT.0) THEN
C                 
                     WRITE(LOGFIL,'(''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT ='',E12.6,'' < '',F12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' (Ylambd NN)'')')
     *                            I,IACTIV,ARGMNT(IACTIV),YNNMIN*DWFACT,
     *                                                    YNNMIN/FACTOR
                 END IF
C
                 IF (ARGMNT(IACTIV).LT.0.0) THEN
                     ARGMNT(IACTIV)=YNNMIN*FACTOR
                 ELSE
                     ARGMNT(IACTIV)=YNNMIN/FACTOR
                 END IF
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (MEMORY(I).GT.MAXMEM) THEN
C
                 WRITE(NOUTPT,'(/,80(''!''),/,''!'',78X,''!'',/,
     *                 ''!     Ylambd NN  parameter was out of range '',
     *                 ''more than MAXMEM='',I2,'' times. ABORT   !'',/,
     *                 ''!'',78X,''!'',/,80(''!''))')     MAXMEM
C
                 IABORT=1
C
                 RETURN
C
             END IF
C
         END IF
C
C=======================================================================
C        I=47 =>> Central Tensor part strength coefficient LambdaPP
C=======================================================================
C
         IF (I.EQ.47.AND.IFTAKE(I).EQ.1) THEN
C
             IACTIV=IACTIV+1
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).GT.CPPMAX*UPFACT) THEN
C
                 FACTOR=CPPMAX/ARGMNT(IACTIV)             
C             
                 IF (LOGWRI.GT.0) THEN
C                 
                     WRITE(LOGFIL,'(''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT='',E12.6,'' >  '',F12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' (Clambd PP)'')')
     *                            I,IACTIV,ARGMNT(IACTIV),CPPMAX*UPFACT,
     *                                                    FACTOR*CPPMAX
                 END IF
C
                 ARGMNT(IACTIV)=FACTOR*CPPMAX
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).LT.CPPMIN*DWFACT) THEN
C
                 FACTOR=ABS(ARGMNT(IACTIV))/(ABS(ARGMNT(IACTIV))+5)                 
C                 
                 IF (LOGWRI.GT.0) THEN
C                 
                     WRITE(LOGFIL,'(''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT ='',E12.6,'' < '',F12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' (Clambd PP)'')')
     *                            I,IACTIV,ARGMNT(IACTIV),CPPMIN*DWFACT,
     *                                                    CPPMIN/FACTOR
                 END IF
C
                 IF (ARGMNT(IACTIV).LT.0.0) THEN
                     ARGMNT(IACTIV)=CPPMIN*FACTOR
                 ELSE
                     ARGMNT(IACTIV)=CPPMIN/FACTOR
                 END IF
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (MEMORY(I).GT.MAXMEM) THEN
C
                 WRITE(NOUTPT,'(/,80(''!''),/,''!'',78X,''!'',/,
     *                 ''!     Clambd PP  parameter was out of range '',
     *                 ''more than MAXMEM='',I2,'' times. ABORT   !'',/,
     *                 ''!'',78X,''!'',/,80(''!''))')     MAXMEM
C
                 IABORT=1
C
                 RETURN
C
             END IF
C
         END IF
C
C=======================================================================
C        I=48 =>> Central Tensor part strength coefficient LambdaPN
C=======================================================================
C
         IF (I.EQ.48.AND.IFTAKE(I).EQ.1) THEN
C
             IACTIV=IACTIV+1
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).GT.CPNMAX*UPFACT) THEN
C
                 FACTOR=CPNMAX/ARGMNT(IACTIV)             
C             
                 IF (LOGWRI.GT.0) THEN
C                 
                     WRITE(LOGFIL,'(''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT='',E12.6,'' >  '',F12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' (Clambd PN)'')')
     *                            I,IACTIV,ARGMNT(IACTIV),CPNMAX*UPFACT,
     *                                                    FACTOR*CPNMAX
                 END IF
C
                 ARGMNT(IACTIV)=FACTOR*CPNMAX
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).LT.CPNMIN*DWFACT) THEN
C
                 FACTOR=ABS(ARGMNT(IACTIV))/(ABS(ARGMNT(IACTIV))+5)                 
C                 
                 IF (LOGWRI.GT.0) THEN
C                 
                     WRITE(LOGFIL,'(''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT ='',E12.6,'' < '',F12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' (Clambd PN)'')')
     *                            I,IACTIV,ARGMNT(IACTIV),CPNMIN*DWFACT,
     *                                                    CPNMIN/FACTOR
                 END IF
C
                 IF (ARGMNT(IACTIV).LT.0.0) THEN
                     ARGMNT(IACTIV)=CPNMIN*FACTOR
                 ELSE
                     ARGMNT(IACTIV)=CPNMIN/FACTOR
                 END IF
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (MEMORY(I).GT.MAXMEM) THEN
C
                 WRITE(NOUTPT,'(/,80(''!''),/,''!'',78X,''!'',/,
     *                 ''!     Clambd PN  parameter was out of range '',
     *                 ''more than MAXMEM='',I2,'' times. ABORT   !'',/,
     *                 ''!'',78X,''!'',/,80(''!''))')     MAXMEM
C
                 IABORT=1
C
                 RETURN
C
             END IF
C
         END IF
C
C=======================================================================
C        I=49 =>> Tensor part strength coefficient LambdaNP
C=======================================================================
C
         IF (I.EQ.49.AND.IFTAKE(I).EQ.1) THEN
C
             IACTIV=IACTIV+1
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).GT.CNPMAX*UPFACT) THEN
C
                 FACTOR=CNPMAX/ARGMNT(IACTIV)             
C             
                 IF (LOGWRI.GT.0) THEN
C                 
                     WRITE(LOGFIL,'(''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT='',E12.6,'' >  '',F12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' (Clambd NP)'')')
     *                            I,IACTIV,ARGMNT(IACTIV),CNPMAX*UPFACT,
     *                                                    FACTOR*CNPMAX
                 END IF
C
                 ARGMNT(IACTIV)=FACTOR*CNPMAX
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).LT.CNPMIN*DWFACT) THEN
C
                 FACTOR=ABS(ARGMNT(IACTIV))/(ABS(ARGMNT(IACTIV))+5)                 
C                 
                 IF (LOGWRI.GT.0) THEN
C                 
                     WRITE(LOGFIL,'(''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT ='',E12.6,'' < '',F12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' (Clambd NP)'')')
     *                            I,IACTIV,ARGMNT(IACTIV),CNPMIN*DWFACT,
     *                                                    CNPMIN/FACTOR
                 END IF
C
                 IF (ARGMNT(IACTIV).LT.0.0) THEN
                     ARGMNT(IACTIV)=CNPMIN*FACTOR
                 ELSE
                     ARGMNT(IACTIV)=CNPMIN/FACTOR
                 END IF
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (MEMORY(I).GT.MAXMEM) THEN
C
                 WRITE(NOUTPT,'(/,80(''!''),/,''!'',78X,''!'',/,
     *                 ''!     Clambd NP  parameter was out of range '',
     *                 ''more than MAXMEM='',I2,'' times. ABORT   !'',/,
     *                 ''!'',78X,''!'',/,80(''!''))')     MAXMEM
C
                 IABORT=1
C
                 RETURN
C
             END IF
C
         END IF
C
C=======================================================================
C        I=50 =>> Tensor part strength coefficient LambdaNN
C=======================================================================
C
         IF (I.EQ.50.AND.IFTAKE(I).EQ.1) THEN
C
             IACTIV=IACTIV+1
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).GT.CNNMAX*UPFACT) THEN
C
                 FACTOR=CNNMAX/ARGMNT(IACTIV)             
C             
                 IF (LOGWRI.GT.0) THEN
C                 
                     WRITE(LOGFIL,'(''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT='',E12.6,'' >  '',F12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' (Clambd NN)'')')
     *                            I,IACTIV,ARGMNT(IACTIV),CNNMAX*UPFACT,
     *                                                    FACTOR*CNNMAX
                 END IF
C
                 ARGMNT(IACTIV)=FACTOR*CNNMAX
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).LT.CNNMIN*DWFACT) THEN
C
                 FACTOR=ABS(ARGMNT(IACTIV))/(ABS(ARGMNT(IACTIV))+5)                 
C                 
                 IF (LOGWRI.GT.0) THEN
C                 
                     WRITE(LOGFIL,'(''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT ='',E12.6,'' < '',F12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' (Clambd NN)'')')
     *                            I,IACTIV,ARGMNT(IACTIV),CNNMIN*DWFACT,
     *                                                    CNNMIN/FACTOR
                 END IF
C
                 IF (ARGMNT(IACTIV).LT.0.0) THEN
                     ARGMNT(IACTIV)=CNNMIN*FACTOR
                 ELSE
                     ARGMNT(IACTIV)=CNNMIN/FACTOR
                 END IF
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (MEMORY(I).GT.MAXMEM) THEN
C
                 WRITE(NOUTPT,'(/,80(''!''),/,''!'',78X,''!'',/,
     *                 ''!     Clambd NN  parameter was out of range '',
     *                 ''more than MAXMEM='',I2,'' times. ABORT   !'',/,
     *                 ''!'',78X,''!'',/,80(''!''))')     MAXMEM
C
                 IABORT=1
C
                 RETURN
C
             END IF
C
         END IF
C
C=======================================================================
C        I=51 =>> V_o central (kappa parametrization)
C=======================================================================
C
         IF (I.EQ.51.AND.IFTAKE(I).EQ.1) THEN
C
             IACTIV=IACTIV+1
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).GT.V0CMAX_KAPPAR*UPFACT) THEN
C
                 FACTOR=V0CMAX_KAPPAR/ARGMNT(IACTIV)             
C             
                 IF (LOGWRI.GT.0) THEN
C                 
                     WRITE(LOGFIL,'(''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT='',E12.6,'' >  '',F12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' (V0CENT_KAPPAR)'')')
     *                            I,IACTIV,ARGMNT(IACTIV),
     *                              V0CMAX_KAPPAR*UPFACT,
     *                              FACTOR*V0CMAX_KAPPAR
                 END IF
C
                 ARGMNT(IACTIV)=FACTOR*V0CMAX_KAPPAR
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).LT.V0CMIN_KAPPAR*DWFACT) THEN
C
                 FACTOR=ABS(ARGMNT(IACTIV))/(ABS(ARGMNT(IACTIV))+5)                 
C                 
                 IF (LOGWRI.GT.0) THEN
C                 
                     WRITE(LOGFIL,'(''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT ='',E12.6,'' < '',F12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' (Clambd NN)'')')
     *                            I,IACTIV,ARGMNT(IACTIV),
     *                              V0CMIN_KAPPAR*DWFACT,
     *                              V0CMIN_KAPPAR/FACTOR
                 END IF
C
                 IF (ARGMNT(IACTIV).LT.0.0) THEN
                     ARGMNT(IACTIV)=V0CMIN_KAPPAR*FACTOR
                 ELSE
                     ARGMNT(IACTIV)=V0CMIN_KAPPAR/FACTOR
                 END IF
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (MEMORY(I).GT.MAXMEM) THEN
C
                 WRITE(NOUTPT,'(/,80(''!''),/,''!'',78X,''!'',/,
     *                 ''! V0CENT_KAPPAR  parameter was out of range '',
     *                 ''more than MAXMEM='',I2,'' times. ABORT   !'',/,
     *                 ''!'',78X,''!'',/,80(''!''))')     MAXMEM
C
                 IABORT=1
C
                 RETURN
C
             END IF
C
         END IF
C
C=======================================================================
C        I=52 =>> kappa V_o central (kappa parametrization)
C=======================================================================
C
         IF (I.EQ.52.AND.IFTAKE(I).EQ.1) THEN
C
             IACTIV=IACTIV+1
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).GT.XKCMAX_KAPPAR*UPFACT) THEN
C
                 FACTOR=XKCMAX_KAPPAR/ARGMNT(IACTIV)             
C             
                 IF (LOGWRI.GT.0) THEN
C                 
                     WRITE(LOGFIL,'(''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT='',E12.6,'' >  '',F12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' (V0CENT_KAPPAR)'')')
     *                            I,IACTIV,ARGMNT(IACTIV),
     *                              XKCMAX_KAPPAR*UPFACT,
     *                              FACTOR*XKCMAX_KAPPAR
                 END IF
C
                 ARGMNT(IACTIV)=FACTOR*XKCMAX_KAPPAR
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).LT.XKCMIN_KAPPAR*DWFACT) THEN
C
                 FACTOR=ABS(ARGMNT(IACTIV))/(ABS(ARGMNT(IACTIV))+5)                 
C                 
                 IF (LOGWRI.GT.0) THEN
C                 
                     WRITE(LOGFIL,'(''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT ='',E12.6,'' < '',F12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' (XK_V0C_KAPPAR)'')')
     *                            I,IACTIV,ARGMNT(IACTIV),
     *                              XKCMIN_KAPPAR*DWFACT,
     *                              XKCMIN_KAPPAR/FACTOR
                 END IF
C
                 IF (ARGMNT(IACTIV).LT.0.0) THEN
                     ARGMNT(IACTIV)=XKCMIN_KAPPAR*FACTOR
                 ELSE
                     ARGMNT(IACTIV)=XKCMIN_KAPPAR/FACTOR
                 END IF
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (MEMORY(I).GT.MAXMEM) THEN
C
                 WRITE(IRESUL,'(/,80(''!''),/,''!'',78X,''!'',/,
     *                 ''! XK_V0C_KAPPAR  parameter was out of range '',
     *                 ''more than MAXMEM='',I2,'' times. ABORT   !'',/,
     *                 ''!'',78X,''!'',/,80(''!''))')     MAXMEM
C
                 IABORT=1
C
                 RETURN
C
             END IF
C
         END IF
C
C=======================================================================
C        I=53 =>> r_o central (kappa parametrization)
C=======================================================================
C
         IF (I.EQ.53.AND.IFTAKE(I).EQ.1) THEN
C
             IACTIV=IACTIV+1
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).GT.R0CMAX_KAPPAR*UPFACT) THEN
C
                 FACTOR=R0CMAX_KAPPAR/ARGMNT(IACTIV)             
C             
                 IF (LOGWRI.GT.0) THEN
C                 
                     WRITE(LOGFIL,'(''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT='',E12.6,'' >  '',F12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' (R0CENT_KAPPAR)'')')
     *                            I,IACTIV,ARGMNT(IACTIV),
     *                              R0CMAX_KAPPAR*UPFACT,
     *                              FACTOR*R0CMAX_KAPPAR
                 END IF
C
                 ARGMNT(IACTIV)=FACTOR*R0CMAX_KAPPAR
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).LT.R0CMIN_KAPPAR*DWFACT) THEN
C
                 FACTOR=ABS(ARGMNT(IACTIV))/(ABS(ARGMNT(IACTIV))+5)                 
C                 
                 IF (LOGWRI.GT.0) THEN
C                 
                     WRITE(LOGFIL,'(''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT ='',E12.6,'' < '',F12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' (R0CENT_KAPPAR)'')')
     *                            I,IACTIV,ARGMNT(IACTIV),
     *                              R0CMIN_KAPPAR*DWFACT,
     *                              R0CMIN_KAPPAR/FACTOR
                 END IF
C
                 IF (ARGMNT(IACTIV).LT.0.0) THEN
                     ARGMNT(IACTIV)=R0CMIN_KAPPAR*FACTOR
                 ELSE
                     ARGMNT(IACTIV)=R0CMIN_KAPPAR/FACTOR
                 END IF
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (MEMORY(I).GT.MAXMEM) THEN
C
                 WRITE(NOUTPT,'(/,80(''!''),/,''!'',78X,''!'',/,
     *                 ''! R0CENT_KAPPAR  parameter was out of range '',
     *                 ''more than MAXMEM='',I2,'' times. ABORT   !'',/,
     *                 ''!'',78X,''!'',/,80(''!''))')     MAXMEM
C
                 IABORT=1
C
                 RETURN
C
             END IF
C
         END IF
C
C=======================================================================
C        I=54 =>> kappa r_o central (kappa parametrization)
C=======================================================================
C
         IF (I.EQ.54.AND.IFTAKE(I).EQ.1) THEN
C
             IACTIV=IACTIV+1
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).GT.XRCMAX_KAPPAR*UPFACT) THEN
C
                 FACTOR=XRCMAX_KAPPAR/ARGMNT(IACTIV)             
C             
                 IF (LOGWRI.GT.0) THEN
C                 
                     WRITE(LOGFIL,'(''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT='',E12.6,'' >  '',F12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' (XK_R0C_KAPPAR)'')')
     *                            I,IACTIV,ARGMNT(IACTIV),
     *                              XRCMAX_KAPPAR*UPFACT,
     *                              FACTOR*XRCMAX_KAPPAR
                 END IF
C
                 ARGMNT(IACTIV)=FACTOR*XRCMAX_KAPPAR
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).LT.XRCMIN_KAPPAR*DWFACT) THEN
C
                 FACTOR=ABS(ARGMNT(IACTIV))/(ABS(ARGMNT(IACTIV))+5)                 
C                 
                 IF (LOGWRI.GT.0) THEN
C                 
                     WRITE(LOGFIL,'(''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT ='',E12.6,'' < '',F12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' (XK_R0C_KAPPAR)'')')
     *                            I,IACTIV,ARGMNT(IACTIV),
     *                              XRCMIN_KAPPAR*DWFACT,
     *                              XRCMIN_KAPPAR/FACTOR
                 END IF
C
                 IF (ARGMNT(IACTIV).LT.0.0) THEN
                     ARGMNT(IACTIV)=XRCMIN_KAPPAR*FACTOR
                 ELSE
                     ARGMNT(IACTIV)=XRCMIN_KAPPAR/FACTOR
                 END IF
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (MEMORY(I).GT.MAXMEM) THEN
C
                 WRITE(NOUTPT,'(/,80(''!''),/,''!'',78X,''!'',/,
     *                 ''! XK_R0C_KAPPAR  parameter was out of range '',
     *                 ''more than MAXMEM='',I2,'' times. ABORT   !'',/,
     *                 ''!'',78X,''!'',/,80(''!''))')     MAXMEM
C
                 IABORT=1
C
                 RETURN
C
             END IF
C
         END IF
C
C=======================================================================
C        I=55 =>> a_o central (kappa parametrization)
C=======================================================================
C
         IF (I.EQ.55.AND.IFTAKE(I).EQ.1) THEN
C
             IACTIV=IACTIV+1
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).GT.A0CMAX_KAPPAR*UPFACT) THEN
C
                 FACTOR=A0CMAX_KAPPAR/ARGMNT(IACTIV)             
C             
                 IF (LOGWRI.GT.0) THEN
C                 
                     WRITE(LOGFIL,'(''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT='',E12.6,'' >  '',F12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' (A0CENT_KAPPAR)'')')
     *                            I,IACTIV,ARGMNT(IACTIV),
     *                              A0CMAX_KAPPAR*UPFACT,
     *                              FACTOR*A0CMAX_KAPPAR
                 END IF
C
                 ARGMNT(IACTIV)=FACTOR*A0CMAX_KAPPAR
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).LT.A0CMIN_KAPPAR*DWFACT) THEN
C
                 FACTOR=ABS(ARGMNT(IACTIV))/(ABS(ARGMNT(IACTIV))+5)                 
C                 
                 IF (LOGWRI.GT.0) THEN
C                 
                     WRITE(LOGFIL,'(''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT ='',E12.6,'' < '',F12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' (A0CENT_KAPPAR)'')')
     *                            I,IACTIV,ARGMNT(IACTIV),
     *                              A0CMIN_KAPPAR*DWFACT,
     *                              A0CMIN_KAPPAR/FACTOR
                 END IF
C
                 IF (ARGMNT(IACTIV).LT.0.0) THEN
                     ARGMNT(IACTIV)=A0CMIN_KAPPAR*FACTOR
                 ELSE
                     ARGMNT(IACTIV)=A0CMIN_KAPPAR/FACTOR
                 END IF
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (MEMORY(I).GT.MAXMEM) THEN
C
                 WRITE(NOUTPT,'(/,80(''!''),/,''!'',78X,''!'',/,
     *                 ''! A0CENT_KAPPAR  parameter was out of range '',
     *                 ''more than MAXMEM='',I2,'' times. ABORT   !'',/,
     *                 ''!'',78X,''!'',/,80(''!''))')     MAXMEM
C
                 IABORT=1
C
                 RETURN
C
             END IF
C
         END IF
C
C=======================================================================
C        I=56 =>> kappa a_o central (kappa parametrization)
C=======================================================================
C
         IF (I.EQ.56.AND.IFTAKE(I).EQ.1) THEN
C
             IACTIV=IACTIV+1
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).GT.XACMAX_KAPPAR*UPFACT) THEN
C
                 FACTOR=XACMAX_KAPPAR/ARGMNT(IACTIV)             
C             
                 IF (LOGWRI.GT.0) THEN
C                 
                     WRITE(LOGFIL,'(''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT='',E12.6,'' >  '',F12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' (XK_A0C_KAPPAR)'')')
     *                            I,IACTIV,ARGMNT(IACTIV),
     *                              XACMAX_KAPPAR*UPFACT,
     *                              FACTOR*XACMAX_KAPPAR
                 END IF
C
                 ARGMNT(IACTIV)=FACTOR*XACMAX_KAPPAR
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).LT.XACMIN_KAPPAR*DWFACT) THEN
C
                 FACTOR=ABS(ARGMNT(IACTIV))/(ABS(ARGMNT(IACTIV))+5)                 
C                 
                 IF (LOGWRI.GT.0) THEN
C                 
                     WRITE(LOGFIL,'(''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT ='',E12.6,'' < '',F12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' (XK_A0C_KAPPAR)'')')
     *                            I,IACTIV,ARGMNT(IACTIV),
     *                              XACMIN_KAPPAR*DWFACT,
     *                              XACMIN_KAPPAR/FACTOR
                 END IF
C
                 IF (ARGMNT(IACTIV).LT.0.0) THEN
                     ARGMNT(IACTIV)=XACMIN_KAPPAR*FACTOR
                 ELSE
                     ARGMNT(IACTIV)=XACMIN_KAPPAR/FACTOR
                 END IF
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (MEMORY(I).GT.MAXMEM) THEN
C
                 WRITE(NOUTPT,'(/,80(''!''),/,''!'',78X,''!'',/,
     *                 ''! XK_A0C_KAPPAR  parameter was out of range '',
     *                 ''more than MAXMEM='',I2,'' times. ABORT   !'',/,
     *                 ''!'',78X,''!'',/,80(''!''))')     MAXMEM
C
                 IABORT=1
C
                 RETURN
C
             END IF
C
         END IF
C
C=======================================================================
C        I=57 =>> lambda_so pure WS (kappa parametrization)
C=======================================================================
C
         IF (I.EQ.57.AND.IFTAKE(I).EQ.1) THEN
C
             IACTIV=IACTIV+1
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).GT.XL_MAX_KAPPAR*UPFACT) THEN
C
                 FACTOR=XL_MAX_KAPPAR/ARGMNT(IACTIV)             
C             
                 IF (LOGWRI.GT.0) THEN
C                 
                     WRITE(LOGFIL,'(''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT='',E12.6,'' >  '',F12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' (V0SORB_KAPPAR)'')')
     *                            I,IACTIV,ARGMNT(IACTIV),
     *                              XL_MAX_KAPPAR*UPFACT,
     *                              FACTOR*XL_MAX_KAPPAR
                 END IF
C
                 ARGMNT(IACTIV)=FACTOR*XL_MAX_KAPPAR
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).LT.XL_MIN_KAPPAR*DWFACT) THEN
C
                 FACTOR=ABS(ARGMNT(IACTIV))/(ABS(ARGMNT(IACTIV))+5)                 
C                 
                 IF (LOGWRI.GT.0) THEN
C                 
                     WRITE(LOGFIL,'(''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT ='',E12.6,'' < '',F12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' (V0SORB_KAPPAR)'')')
     *                            I,IACTIV,ARGMNT(IACTIV),
     *                              XL_MIN_KAPPAR*DWFACT,
     *                              XL_MIN_KAPPAR/FACTOR
                 END IF
C
                 IF (ARGMNT(IACTIV).LT.0.0) THEN
                     ARGMNT(IACTIV)=XL_MIN_KAPPAR*FACTOR
                 ELSE
                     ARGMNT(IACTIV)=XL_MIN_KAPPAR/FACTOR
                 END IF
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (MEMORY(I).GT.MAXMEM) THEN
C
                 WRITE(NOUTPT,'(/,80(''!''),/,''!'',78X,''!'',/,
     *                 ''! V0SORB_KAPPAR  parameter was out of range '',
     *                 ''more than MAXMEM='',I2,'' times. ABORT   !'',/,
     *                 ''!'',78X,''!'',/,80(''!''))')     MAXMEM
C
                 IABORT=1
C
                 RETURN
C
             END IF
C
         END IF
C
C=======================================================================
C        I=58 =>> kappa lambda_so pure WS (kappa parametrization)
C=======================================================================
C
         IF (I.EQ.58.AND.IFTAKE(I).EQ.1) THEN
C
             IACTIV=IACTIV+1
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).GT.XKSMAX_KAPPAR*UPFACT) THEN
C
                 FACTOR=XKSMAX_KAPPAR/ARGMNT(IACTIV)             
C             
                 IF (LOGWRI.GT.0) THEN
C                 
                     WRITE(LOGFIL,'(''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT='',E12.6,'' >  '',F12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' (XK_V0S_KAPPAR)'')')
     *                            I,IACTIV,ARGMNT(IACTIV),
     *                              XKSMAX_KAPPAR*UPFACT,
     *                              FACTOR*XKSMAX_KAPPAR
                 END IF
C
                 ARGMNT(IACTIV)=FACTOR*XKSMAX_KAPPAR
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).LT.XKSMIN_KAPPAR*DWFACT) THEN
C
                 FACTOR=ABS(ARGMNT(IACTIV))/(ABS(ARGMNT(IACTIV))+5)                 
C                 
                 IF (LOGWRI.GT.0) THEN
C                 
                     WRITE(LOGFIL,'(''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT ='',E12.6,'' < '',F12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' (XK_V0S_KAPPAR)'')')
     *                            I,IACTIV,ARGMNT(IACTIV),
     *                              XKSMIN_KAPPAR*DWFACT,
     *                              XKSMIN_KAPPAR/FACTOR
                 END IF
C
                 IF (ARGMNT(IACTIV).LT.0.0) THEN
                     ARGMNT(IACTIV)=XKSMIN_KAPPAR*FACTOR
                 ELSE
                     ARGMNT(IACTIV)=XKSMIN_KAPPAR/FACTOR
                 END IF
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (MEMORY(I).GT.MAXMEM) THEN
C
                 WRITE(NOUTPT,'(/,80(''!''),/,''!'',78X,''!'',/,
     *                 ''! XK_V0S_KAPPAR  parameter was out of range '',
     *                 ''more than MAXMEM='',I2,'' times. ABORT   !'',/,
     *                 ''!'',78X,''!'',/,80(''!''))')     MAXMEM
C
                 IABORT=1
C
                 RETURN
C
             END IF
C
         END IF
C
C=======================================================================
C        I=59 =>> r_so pure WS (kappa parametrization)
C=======================================================================
C
         IF (I.EQ.59.AND.IFTAKE(I).EQ.1) THEN
C
             IACTIV=IACTIV+1
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).GT.R0SMAX_KAPPAR*UPFACT) THEN
C
                 FACTOR=R0SMAX_KAPPAR/ARGMNT(IACTIV)             
C             
                 IF (LOGWRI.GT.0) THEN
C                 
                     WRITE(LOGFIL,'(''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT='',E12.6,'' >  '',F12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' (R0SORB_KAPPAR)'')')
     *                            I,IACTIV,ARGMNT(IACTIV),
     *                              R0SMAX_KAPPAR*UPFACT,
     *                              FACTOR*R0SMAX_KAPPAR
                 END IF
C
                 ARGMNT(IACTIV)=FACTOR*R0SMAX_KAPPAR
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).LT.R0SMIN_KAPPAR*DWFACT) THEN
C
                 FACTOR=ABS(ARGMNT(IACTIV))/(ABS(ARGMNT(IACTIV))+5)                 
C                 
                 IF (LOGWRI.GT.0) THEN
C                 
                     WRITE(LOGFIL,'(''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT ='',E12.6,'' < '',F12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' (R0SORB_KAPPAR)'')')
     *                            I,IACTIV,ARGMNT(IACTIV),
     *                              R0SMIN_KAPPAR*DWFACT,
     *                              R0SMIN_KAPPAR/FACTOR
                 END IF
C
                 IF (ARGMNT(IACTIV).LT.0.0) THEN
                     ARGMNT(IACTIV)=R0SMIN_KAPPAR*FACTOR
                 ELSE
                     ARGMNT(IACTIV)=R0SMIN_KAPPAR/FACTOR
                 END IF
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (MEMORY(I).GT.MAXMEM) THEN
C
                 WRITE(NOUTPT,'(/,80(''!''),/,''!'',78X,''!'',/,
     *                 ''! R0SORB_KAPPAR  parameter was out of range '',
     *                 ''more than MAXMEM='',I2,'' times. ABORT   !'',/,
     *                 ''!'',78X,''!'',/,80(''!''))')     MAXMEM
C
                 IABORT=1
C
                 RETURN
C
             END IF
C
         END IF
C
C=======================================================================
C        I=60 =>> kappa r_so pure WS (kappa parametrization)
C=======================================================================
C
         IF (I.EQ.60.AND.IFTAKE(I).EQ.1) THEN
C
             IACTIV=IACTIV+1
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).GT.XRSMAX_KAPPAR*UPFACT) THEN
C
                 FACTOR=XRSMAX_KAPPAR/ARGMNT(IACTIV)             
C             
                 IF (LOGWRI.GT.0) THEN
C                 
                     WRITE(LOGFIL,'(''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT='',E12.6,'' >  '',F12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' (XK_R0S_KAPPAR)'')')
     *                            I,IACTIV,ARGMNT(IACTIV),
     *                              XRSMAX_KAPPAR*UPFACT,
     *                              FACTOR*XRSMAX_KAPPAR
                 END IF
C
                 ARGMNT(IACTIV)=FACTOR*XRSMAX_KAPPAR
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).LT.XRSMIN_KAPPAR*DWFACT) THEN
C
                 FACTOR=ABS(ARGMNT(IACTIV))/(ABS(ARGMNT(IACTIV))+5)                 
C                 
                 IF (LOGWRI.GT.0) THEN
C                 
                     WRITE(LOGFIL,'(''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT ='',E12.6,'' < '',F12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' (XK_R0S_KAPPAR)'')')
     *                            I,IACTIV,ARGMNT(IACTIV),
     *                              XRSMIN_KAPPAR*DWFACT,
     *                              XRSMIN_KAPPAR/FACTOR
                 END IF
C
                 IF (ARGMNT(IACTIV).LT.0.0) THEN
                     ARGMNT(IACTIV)=XRSMIN_KAPPAR*FACTOR
                 ELSE
                     ARGMNT(IACTIV)=XRSMIN_KAPPAR/FACTOR
                 END IF
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (MEMORY(I).GT.MAXMEM) THEN
C
                 WRITE(NOUTPT,'(/,80(''!''),/,''!'',78X,''!'',/,
     *                 ''! XK_R0S_KAPPAR  parameter was out of range '',
     *                 ''more than MAXMEM='',I2,'' times. ABORT   !'',/,
     *                 ''!'',78X,''!'',/,80(''!''))')     MAXMEM
C
                 IABORT=1
C
                 RETURN
C
             END IF
C
         END IF
C
C=======================================================================
C        I=61 =>> a_so pure WS (kappa parametrization)
C=======================================================================
C
         IF (I.EQ.61.AND.IFTAKE(I).EQ.1) THEN
C
             IACTIV=IACTIV+1
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).GT.A0SMAX_KAPPAR*UPFACT) THEN
C
                 FACTOR=A0SMAX_KAPPAR/ARGMNT(IACTIV)             
C             
                 IF (LOGWRI.GT.0) THEN
C                 
                     WRITE(LOGFIL,'(''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT='',E12.6,'' >  '',F12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' (A0SORB_KAPPAR)'')')
     *                            I,IACTIV,ARGMNT(IACTIV),
     *                              A0SMAX_KAPPAR*UPFACT,
     *                              FACTOR*A0SMAX_KAPPAR
                 END IF
C
                 ARGMNT(IACTIV)=FACTOR*A0SMAX_KAPPAR
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).LT.A0SMIN_KAPPAR*DWFACT) THEN
C
                 FACTOR=ABS(ARGMNT(IACTIV))/(ABS(ARGMNT(IACTIV))+5)                 
C                 
                 IF (LOGWRI.GT.0) THEN
C                 
                     WRITE(LOGFIL,'(''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT ='',E12.6,'' < '',F12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' (A0SORB_KAPPAR)'')')
     *                            I,IACTIV,ARGMNT(IACTIV),
     *                              A0SMIN_KAPPAR*DWFACT,
     *                              A0SMIN_KAPPAR/FACTOR
                 END IF
C
                 IF (ARGMNT(IACTIV).LT.0.0) THEN
                     ARGMNT(IACTIV)=A0SMIN_KAPPAR*FACTOR
                 ELSE
                     ARGMNT(IACTIV)=A0SMIN_KAPPAR/FACTOR
                 END IF
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (MEMORY(I).GT.MAXMEM) THEN
C
                 WRITE(NOUTPT,'(/,80(''!''),/,''!'',78X,''!'',/,
     *                 ''! A0SORB_KAPPAR  parameter was out of range '',
     *                 ''more than MAXMEM='',I2,'' times. ABORT   !'',/,
     *                 ''!'',78X,''!'',/,80(''!''))')     MAXMEM
C
                 IABORT=1
C
                 RETURN
C
             END IF
C
         END IF
C
C=======================================================================
C        I=62 =>> kappa a_so pure WS (kappa parametrization)
C=======================================================================
C
         IF (I.EQ.62.AND.IFTAKE(I).EQ.1) THEN
C
             IACTIV=IACTIV+1
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).GT.XASMAX_KAPPAR*UPFACT) THEN
C
                 FACTOR=XASMAX_KAPPAR/ARGMNT(IACTIV)             
C             
                 IF (LOGWRI.GT.0) THEN
C                 
                     WRITE(LOGFIL,'(''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT='',E12.6,'' >  '',F12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' (XK_A0S_KAPPAR)'')')
     *                            I,IACTIV,ARGMNT(IACTIV),
     *                              XASMAX_KAPPAR*UPFACT,
     *                              FACTOR*XASMAX_KAPPAR
                 END IF
C
                 ARGMNT(IACTIV)=FACTOR*XASMAX_KAPPAR
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (ARGMNT(IACTIV).LT.XASMIN_KAPPAR*DWFACT) THEN
C
                 FACTOR=ABS(ARGMNT(IACTIV))/(ABS(ARGMNT(IACTIV))+5)
C                 
                 IF (LOGWRI.GT.0) THEN
C                 
                     WRITE(LOGFIL,'(''I='',I2,'' IACTIV ='',I3,
     *                            '' ARGMNT ='',E12.6,'' < '',F12.6,
     *                            '' NEWARG ='',F12.6,
     *                            '' (XK_A0S_KAPPAR)'')')
     *                            I,IACTIV,ARGMNT(IACTIV),
     *                              XASMIN_KAPPAR*DWFACT,
     *                              XASMIN_KAPPAR/FACTOR
                 END IF
C
                 IF (ARGMNT(IACTIV).LT.0.0) THEN
                     ARGMNT(IACTIV)=XASMIN_KAPPAR*FACTOR
                 ELSE
                     ARGMNT(IACTIV)=XASMIN_KAPPAR/FACTOR
                 END IF
C
                 ACTION(I)='#'
C
                 MEMORY(I)=MEMORY(I)+1                 
C
             END IF
C_______________________________________________________________________
C
             IF (MEMORY(I).GT.MAXMEM) THEN
C
                 WRITE(NOUTPT,'(/,80(''!''),/,''!'',78X,''!'',/,
     *                 ''! XK_A0S_KAPPAR  parameter was out of range '',
     *                 ''more than MAXMEM='',I2,'' times. ABORT   !'',/,
     *                 ''!'',78X,''!'',/,80(''!''))')     MAXMEM
C
                 IABORT=1
C
                 RETURN
C
             END IF
C
         END IF
C
C=======================================================================
C
      END DO
C
C=======================================================================
C
      IF (LOGWRI.GT.4) THEN
          WRITE(LOGFIL,'(12X,''Exiting  VERIFY'')')
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

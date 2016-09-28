C FILE NAME = wspher06_potnts_15.f ! Keep this symbol: $ident@string$
C
C=======================================================================
C=======================================================================
C                MATRIX ELEMENTS, POTENTIALS, HAMILTONIANS
C=======================================================================
C=======================================================================
C
      SUBROUTINE INTROD(INUCLI)
      
      INCLUDE   'MATDIM/NDPARS.f'
      INCLUDE   'MATDIM/NDNUCL.f'
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
     *       /HBAR_V/ HBAR_C
      COMMON
     *       /DATA_V/ AOSCIL_PROTON(1:NDNUCL),
     *                AOSCIL_NEUTRS(1:NDNUCL),
     *                HOMEGA(1:NDNUCL)      
      COMMON
     *       /NUCLEU/ IZ_FIX,IN_FIX,ISOSPI 
     *       /DENSTR/ IFDENS
     *       /TENSOR/ IFTENS
      COMMON
     *       /TESFIT/ IFTEST,IFFITS,IFMESH,ISPE_P,IFMOCA
      COMMON
     *       /CTSHIF/ V0SHIF_PROTON(1:NDNUCL),
     *                V0SHIF_NEUTRS(1:NDNUCL)
      COMMON
     *       /IFKAPP/ IFK_VC,IFK_VS,IFK_RC,IFK_RS,IFK_AC,IFK_AS
     *       /MASSIV/ IMASIV_PRODUCC
      COMMON
     *       /PARCOR/ IFCORR,IFRCVC,IFRCAC,IFVCAC,IFRSVS,IFRSAS,IFVSAS
      COMMON
     *       /FACCOU/ COUFAC
      COMMON
     *       /PRFILE/ ISCREN,LOGFIL,LOGCHP,LOGCHN,LOGHIP,LOGHIN,
     *                LOGPRO,LOGNEU,LOGWRI,LSCREN,IRESUL,ICONVE,
     *                NOUTPT,LOGAUX,LOGENE,LOGRAD,LOGGRD,LOGMSH,
     *                                                   LOGBIS 
C
      DATA
     *     XMASSP /938.27231 /,! Units: MeV/c**2; ERROR: (+/-)0.00028
     *     XMASSN /939.56563 / ! Units: MeV/c**2; ERROR: (+/-)0.00028
C
C=======================================================================
C     Introducing all the constant constants
C=======================================================================
C
      IF (LOGWRI.GT.4) THEN
          WRITE(LOGFIL,'(/,9X,''Entering INTROD'')')
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
C
C     For one single run we optionally apply 
C     the shift to the central potential
C
      IF (IFTEST.EQ.1) THEN
          V0CENT_PROTON=V0CENT_PROTON+V0SHIF_PROTON(INUCLI)
          V0CENT_NEUTRS=V0CENT_NEUTRS+V0SHIF_NEUTRS(INUCLI)
      END IF
C      
C=======================================================================
C=======================================================================
C     Spin-orbit potential with density dependence
C=======================================================================
C=======================================================================
C     
      ALAMPP=PARPOT(39)
      ALAMPN=PARPOT(40)
      ALAMNP=PARPOT(41)
      ALAMNN=PARPOT(42)
C      
C=======================================================================
C=======================================================================
C     Spin-orbit TENSOR Potential
C=======================================================================
C=======================================================================
C  
      TLAMPP=PARPOT(43)
      TLAMPN=PARPOT(44)
      TLAMNP=PARPOT(45)
      TLAMNN=PARPOT(46)
C      
C=======================================================================
C=======================================================================
C     Central TENSOR Potential
C=======================================================================
C=======================================================================
C      
      CLAMPP=PARPOT(47)
      CLAMPN=PARPOT(48)
      CLAMNP=PARPOT(49)
      CLAMNN=PARPOT(50)
C      
C=======================================================================
C=======================================================================
C     KAPpa PARametrization = KAPPAR
C=======================================================================
C=======================================================================
C      
      V0CENT_KAPPAR=PARPOT(51)
      XK_V0C_KAPPAR=PARPOT(52)
      R0CENT_KAPPAR=PARPOT(53)
      XK_R0C_KAPPAR=PARPOT(54)
      A0CENT_KAPPAR=PARPOT(55)
      XK_A0C_KAPPAR=PARPOT(56)
C
      V0SORB_KAPPAR=PARPOT(57)
      XK_V0S_KAPPAR=PARPOT(58)
      R0SORB_KAPPAR=PARPOT(59)
      XK_R0S_KAPPAR=PARPOT(60)
      A0SORB_KAPPAR=PARPOT(61)
      XK_A0S_KAPPAR=PARPOT(62)
C
C=======================================================================
C=======================================================================
C 
      A_MASS=IZ_FIX+IN_FIX
C
C=======================================================================
C=======================================================================
C         Proton part
C=======================================================================
C=======================================================================
C
      ISOSPI=1
C
      SIGISO=(IN_FIX-IZ_FIX)/A_MASS
C
C-----------------------------------------------------------------------      
C
      IF (IFK_VC.EQ.1) THEN
          VKACEN_PROTON=V0CENT_KAPPAR*(1+SIGISO*XK_V0C_KAPPAR)
      ELSE
          VKACEN_PROTON=V0CENT_PROTON
          V0CENT_KAPPAR=0.0
          XK_V0C_KAPPAR=0.0
      END IF
C
      IF (IFK_RC.EQ.1) THEN
          RKACEN_PROTON=
     *    R0CENT_KAPPAR*(1+SIGISO*XK_R0C_KAPPAR)*A_MASS**(1./3.)
      ELSE
          RKACEN_PROTON=R0CENT_PROTON*A_MASS**(1./3.)
          R0CENT_KAPPAR=0.0
          XK_R0C_KAPPAR=0.0
      END IF 
C     
      IF (IFK_AC.EQ.1) THEN
          AKACEN_PROTON=A0CENT_KAPPAR*(1+SIGISO*XK_A0C_KAPPAR)
      ELSE
          AKACEN_PROTON=A0CENT_PROTON
          A0CENT_KAPPAR=0.0
          XK_A0C_KAPPAR=0.0
      END IF
C     
      R0COUL=R0CENT_PROTON
      RKACOU=COUFAC*R0COUL*(1+SIGISO*XK_COU)*A_MASS**(1./3.)
C
C-----------------------------------------------------------------------      
C
      IF (IFK_VS.EQ.1) THEN
          VKASOR_PROTON=V0SORB_KAPPAR*(1+SIGISO*XK_V0S_KAPPAR)
      ELSE
          VKASOR_PROTON=V0SORB_PROTON
          V0SORB_KAPPAR=0.0
          XK_V0S_KAPPAR=0.0
      END IF
C
      IF (IFK_RS.EQ.1) THEN
          RKASOR_PROTON=
     *    R0SORB_KAPPAR*(1+SIGISO*XK_R0S_KAPPAR)*A_MASS**(1./3.)
      ELSE
          RKASOR_PROTON=R0SORB_PROTON*A_MASS**(1./3.)
          R0SORB_KAPPAR=0.0
          XK_R0S_KAPPAR=0.0
      END IF 
C     
      IF (IFK_AS.EQ.1) THEN
          AKASOR_PROTON=A0SORB_KAPPAR*(1+SIGISO*XK_A0S_KAPPAR)
      ELSE
          AKASOR_PROTON=A0SORB_PROTON
          A0SORB_KAPPAR=0.0
          XK_A0S_KAPPAR=0.0
      END IF
C
C-----------------------------------------------------------------------      
C
      RKAEFM_PROTON=
     *R0EFFM_PROTON*(1-SIGISO*XK_REF_PROTON)*A_MASS**(1./3.)
C     
      AKAEFM_PROTON=A0EFFM_PROTON*(1-SIGISO*XK_AEF_PROTON)
      VKAEFM_PROTON=V0EFFM_PROTON*(1-SIGISO*XK_LEF_PROTON) 
C      
C=======================================================================
C
      IF (IFCORR.EQ.1 .AND. IFRCVC.EQ.1) THEN
C
          VKACEN_PROTON=V0CENT_PROTON
C
          R0CENT_PROTON=R0CV0C(VKACEN_PROTON)
          RKACEN_PROTON=R0CENT_PROTON*A_MASS**(1./3.)
C
          AKACEN_PROTON=A0CENT_PROTON
              
          R0COUL=R0CENT_PROTON
          RKACOU=COUFAC*R0COUL*(1+SIGISO*XK_COU)*A_MASS**(1./3.)
C
      END IF
C      
C=======================================================================
C=======================================================================
C         Neutron part
C=======================================================================
C=======================================================================
C
      ISOSPI=0
C
      SIGISO=(IN_FIX-IZ_FIX)/A_MASS
C
C-----------------------------------------------------------------------      
C
      IF (IFK_VC.EQ.1) THEN
          VKACEN_NEUTRS=V0CENT_KAPPAR*(1-SIGISO*XK_V0C_KAPPAR)
      ELSE
          VKACEN_NEUTRS=V0CENT_NEUTRS
          V0CENT_KAPPAR=0.0
          XK_V0C_KAPPAR=0.0
      END IF
C
      IF (IFK_RC.EQ.1) THEN
          RKACEN_NEUTRS=
     *    R0CENT_KAPPAR*(1-SIGISO*XK_R0C_KAPPAR)*A_MASS**(1./3.)
      ELSE
          RKACEN_NEUTRS=R0CENT_NEUTRS*A_MASS**(1./3.)
          R0CENT_KAPPAR=0.0
          XK_R0C_KAPPAR=0.0
      END IF
C
      IF (IFK_AC.EQ.1) THEN
          AKACEN_NEUTRS=A0CENT_KAPPAR*(1-SIGISO*XK_A0C_KAPPAR)
      ELSE
          AKACEN_NEUTRS=A0CENT_NEUTRS
          A0CENT_KAPPAR=0.0
          XK_A0C_KAPPAR=0.0
      END IF
C
C-----------------------------------------------------------------------      
C
      IF (IFK_VS.EQ.1) THEN
          VKASOR_NEUTRS=V0SORB_KAPPAR*(1-SIGISO*XK_V0S_KAPPAR)
      ELSE
          VKASOR_NEUTRS=V0SORB_NEUTRS
          V0SORB_KAPPAR=0.0
          XK_V0S_KAPPAR=0.0
      END IF
C
      IF (IFK_RS.EQ.1) THEN
          RKASOR_NEUTRS=
     *    R0SORB_KAPPAR*(1-SIGISO*XK_R0S_KAPPAR)*A_MASS**(1./3.)
      ELSE
          RKASOR_NEUTRS=R0SORB_NEUTRS*A_MASS**(1./3.)
          R0SORB_KAPPAR=0.0
          XK_R0S_KAPPAR=0.0
      END IF
C
      IF (IFK_AS.EQ.1) THEN
          AKASOR_NEUTRS=A0SORB_KAPPAR*(1-SIGISO*XK_A0S_KAPPAR)
      ELSE
          AKASOR_NEUTRS=A0SORB_NEUTRS
          A0SORB_KAPPAR=0.0
          XK_A0S_KAPPAR=0.0
      END IF
C
C-----------------------------------------------------------------------      
C
      RKAEFM_NEUTRS=
     *R0EFFM_NEUTRS*(1-SIGISO*XK_REF_NEUTRS)*A_MASS**(1./3.)
C     
      AKAEFM_NEUTRS=A0EFFM_NEUTRS*(1-SIGISO*XK_AEF_NEUTRS)
      VKAEFM_NEUTRS=V0EFFM_NEUTRS*(1-SIGISO*XK_LEF_NEUTRS)
C
C=======================================================================
C
      IF (IFCORR.EQ.1 .AND. IFRCVC.EQ.1) THEN
C
          VKACEN_NEUTRS=V0CENT_NEUTRS
C
          R0CENT_NEUTRS=R0CV0C(VKACEN_NEUTRS)
          RKACEN_NEUTRS=R0CENT_NEUTRS*A_MASS**(1./3.)
C
          AKACEN_NEUTRS=A0CENT_NEUTRS
C
      END IF
C
C=======================================================================
C
      IF (LOGWRI.GT.0) THEN
C
          WRITE(NOUTPT,'(/,80(''#''),/,''#'',T80,''#'',/,
     *           ''#'',27X,''P A R A M E T E R S '',T80,''#'',/,
     *           ''#'',T80,''#'',/,
     *           ''#'',20X,''IZ_FIX='',I3,10X,''IN_FIX='',I3,
     *                                           T80,''#'',/,
     *           ''#'',T80,''#'',/,80(''#''))')IZ_FIX,IN_FIX
      
          WRITE(NOUTPT,'(/,80(''#''),/,''#'',''  Protons'',T80,''#'',/,
     *            ''#  VKACEN  VKASOR    '',
     *            ''RKACEN  RKASOR    '',
     *            ''AKACEN  AKASOR    RKACOU'',T80,''#'',/,
     *            ''#'',2F8.3,2X,2F8.4,2X,2F8.4,2X,F8.4,T80,''#'',/,
     *            ''#'',T80,''#'',/,80(''#''))')
     *
     *            VKACEN_PROTON,VKASOR_PROTON, 
     *            RKACEN_PROTON,RKASOR_PROTON,
     *            AKACEN_PROTON,AKASOR_PROTON,RKACOU
C
          WRITE(NOUTPT,'(/,80(''#''),/,''#'',''  Neutrons'',T80,''#'',/,
     *            ''#  VKACEN  VKASOR    '',
     *            ''RKACEN  RKASOR    '',
     *            ''AKACEN  AKASOR'',T80,''#'',/,
     *            ''#'',2F8.3,2X,2F8.4,2X,2F8.4,T80,''#'',/,
     *            ''#'',T80,''#'',/,80(''#''))')
     *
     *            VKACEN_NEUTRS,VKASOR_NEUTRS,
     *            RKACEN_NEUTRS,RKASOR_NEUTRS,
     *            AKACEN_NEUTRS,AKASOR_NEUTRS
C
          WRITE(NOUTPT,'(/,80(''#''),/,''#'',''  Kappa '',
     *            ''Parametrisation'',T80,''#'',/,
     *            ''#  V0CENT  XK_V0C    '',
     *            ''R0CENT  XK_R0C    '',
     *            ''A0CENT  XK_R0C'',T80,''#'',/,
     *            ''#'',2F8.3,2X,2F8.4,2X,2F8.4,T80,''#'',/,
     *            ''#'',T80,''#'',/,
     *            ''#  V0SORB  XK_LAM    '',
     *            ''R0SORB  XK_RSO    '',
     *            ''A0SORB  XK_ASO'',T80,''#'',/,
     *            ''#'',2F8.3,2X,2F8.4,2X,2F8.4,T80,''#'',/,
     *            ''#'',T80,''#'',/,''#'',T64,''R0COUL  XK_COU  #'',/,
     *            ''#'',T64,F6.4,F8.3,''  #'',/,
     *            ''#'',T80,''#'',/,80(''#''))')
     *
     *            V0CENT_KAPPAR,XK_V0C_KAPPAR,
     *            R0CENT_KAPPAR,XK_R0C_KAPPAR,
     *            A0CENT_KAPPAR,XK_A0C_KAPPAR,
     *            V0SORB_KAPPAR,XK_V0S_KAPPAR,
     *            R0SORB_KAPPAR,XK_R0S_KAPPAR,
     *            A0SORB_KAPPAR,XK_A0S_KAPPAR,R0COUL,XK_COU
C
          WRITE(NOUTPT,'(/,80(''#''),/,''#'',''  No-Kappa Protons'',
     *                                                    T80,''#'',/,
     *            ''#  V0CENT  V0SORB    '',
     *            ''R0CENT  R0SORB    '',
     *            ''A0CENT  A0SORB'',T80,''#'',/,
     *            ''#'',2F8.3,2X,2F8.4,2X,2F8.4,T80,''#'',/,
     *            ''#'',T80,''#'',/,''#'',T64,''R0COUL'',T80,''#'',/,
     *            ''#'',T64,F6.4,T80,''#'',/,
     *            ''#'',T80,''#'',/,80(''#''))')
     *
     *            V0CENT_PROTON,V0SORB_PROTON,
     *            R0CENT_PROTON,R0SORB_PROTON,
     *            A0CENT_PROTON,A0SORB_PROTON,R0COUL
C
          WRITE(NOUTPT,'(/,80(''#''),/,''#'',''  No-Kappa Neutron'',
     *                                                    T80,''#'',/,
     *            ''#  V0CENT  V0SORB    '',
     *            ''R0CENT  R0SORB    '',
     *            ''A0CENT  A0SORB'',T80,''#'',/,
     *            ''#'',2F8.3,2X,2F8.4,2X,2F8.4,T80,''#'',/,
     *            ''#'',T80,''#'',/,80(''#''))')
     *
     *            V0CENT_NEUTRS,V0SORB_NEUTRS,
     *            R0CENT_NEUTRS,R0SORB_NEUTRS,
     *            A0CENT_NEUTRS,A0SORB_NEUTRS
C
          IF (IFDENS.EQ.1)
     *        WRITE(NOUTPT,'(/,80(''#''),/,''#'',
     *            ''  Lambdas'',T80,''#'',/,
     *            ''#'',9X,''ALAMPP    ALAMPN    ALAMNP    ALAMNN'',
     *                                                     T80,''#'',/,
     *            ''#'',7X,4(F8.4,2X),T80,''#'',/,
     *            ''#'',T80,''#'',/,
     *            ''#'',9X,''TLAMPP    TLAMPN    TLAMNP    TLAMNN'',
     *                                                     T80,''#'',/,
     *            ''#'',7X,4(F8.4,2X),T80,''#'',/,
     *            ''#'',T80,''#'',/,
     *            ''#'',9X,''CLAMPP    CLAMPN    CLAMNP    CLAMNN'',
     *                                                     T80,''#'',/,
     *            ''#'',7X,4(F8.4,2X),T80,''#'',/,
     *            ''#'',T80,''#'',/,
     *            ''#'',T80,''#'',/,80(''#''))')
     *
     *            ALAMPP,ALAMPN,ALAMNP,ALAMNN,
     *            TLAMPP,TLAMPN,TLAMNP,TLAMNN,
     *            CLAMPP,CLAMPN,CLAMNP,CLAMNN
C
      END IF
C
C=======================================================================
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
c      IF (LOGWRI.GT.0) THEN
c          WRITE(LOGFIL,'(/,9X,''Exiting INTROD'')')
c      END IF
C
C=======================================================================
C
      RETURN
      END
C
C=======================================================================
C=======================================================================
C
      FUNCTION R0CV0C(VARGUM)
C
      COMMON
     *       /NUCLEU/ IZ_FIX,IN_FIX,ISOSPI 
C
C=======================================================================
C     Parabolic parametric correlation between r_o^c and V_o^c
C=======================================================================
C
      IF (IZ_FIX.EQ.82 .AND. IN_FIX.EQ.126) THEN
          IF (ISOSPI.EQ.1) THEN
              ACOEFF=0.00023262
              BCOEFF=0.04561273
              CCOEFF=3.16142711
          ELSE
              ACOEFF=0.00020081
              BCOEFF=0.03822447
              CCOEFF=2.63782567
          END IF
      END IF
C
      R0CV0C=ACOEFF*VARGUM**2 + BCOEFF*VARGUM + CCOEFF
C
C=======================================================================
C
      END
C
C=======================================================================
C=======================================================================
C
      SUBROUTINE INTROD_KAPPAS(INUCLI)
      
      INCLUDE   'MATDIM/NDPARS.f'
      INCLUDE   'MATDIM/NDNUCL.f'
C
      CHARACTER
     *          NUCSYM*6,INPSYM*6
C
      COMMON
     *       /HBAR_V/ HBAR_C
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
     *       /KAPPAS/ V0CENK,XKAPPA_V0CENK,R0CENK,XKAPPA_R0CENK,
     *                                     A0CENK,XKAPPA_A0CENK,
     *                V0SORK,XKAPPA_V0SORK,R0SORK,XKAPPA_R0SORK,
     *                                     A0SORK,XKAPPA_A0SORK
      COMMON
     *       /NUCLEU/ IZ_FIX,IN_FIX,ISOSPI 
      COMMON
     *       /TESFIT/ IFTEST,IFFITS,IFMESH,ISPE_P,IFMOCA
      COMMON
     *       /CTSHIF/ V0SHIF_PROTON(1:NDNUCL),
     *                V0SHIF_NEUTRS(1:NDNUCL)
      COMMON
     *       /IFKAPP/ IFK_VC,IFK_VS,IFK_RC,IFK_RS,IFK_AC,IFK_AS
C
      COMMON
     *       /PRFILE/ ISCREN,LOGFIL,LOGCHP,LOGCHN,LOGHIP,LOGHIN,
     *                LOGPRO,LOGNEU,LOGWRI,LSCREN,IRESUL,ICONVE,
     *                NOUTPT,LOGAUX,LOGENE,LOGRAD,LOGGRD,LOGMSH,
     *                                                   LOGBIS 
C
C=======================================================================
C
      WRITE(LOGFIL,'(''Entering INTROD_KAPPAS with IZ_FIX= '',I3,
     *               '' and IN_FIX= '',I3)') IZ_FIX,IN_FIX
C
C=======================================================================
C
C     Initializing things
C
      V0CENK=0.
      XKAPPA_V0CENK=0.
C
      R0CENK=0.
      XKAPPA_R0CENK=0.
C
      A0CENK=0.
      XKAPPA_A0CENK=0.
C_______________________________________________________________________
C
      V0SORK=0.
      XKAPPA_V0SORK=0.
C
      R0SORK=0.
      XKAPPA_R0SORK=0.
C
      A0SORK=0.
      XKAPPA_A0SORK=0.
C
C=======================================================================
C
      IZ_FIX=NUMB_Z(INUCLI)
      IN_FIX=NUMB_N(INUCLI)
      A_FIXX=IZ_FIX+IN_FIX
      FRACTN=(IN_FIX-IZ_FIX)/A_FIXX
C
C=======================================================================
C
C     Identifying the nucleus over which we have minimized 
C     in order to calculate V_o and kappa from it.
C
      DO JNUCLI=1,LDNUCL
         IF (ITAKNU(JNUCLI).EQ.1) THEN
             IZ_AUX=NUMB_Z(JNUCLI)
             IN_AUX=NUMB_N(JNUCLI)
             A_AUXX=IZ_AUX+IN_AUX
         END IF
      END DO
C
C=======================================================================
C
      WRITE(LOGFIL,'(''IZ_AUX= '',I3,'' and IN_AUX= '',I3)') 
     *                 IZ_AUX,              IN_AUX
C
C=======================================================================
C=======================================================================
C
C     Introducing the "kappa" variant to the CENTRAL potential
C     parameters
C
C     As they come from the minimization
C
      V_PROT=PARPOT(01)
      V_NEUT=PARPOT(21)
C
      R_PROT=PARPOT(02)
      R_NEUT=PARPOT(22)
C
      A_PROT=PARPOT(03)
      A_NEUT=PARPOT(23)
C
C=======================================================================
C
      XKAPPA_V0CENK= (V_PROT-V_NEUT)*A_AUXX
     *             /((V_PROT+V_NEUT)*(IN_AUX-IZ_AUX))
C
      V0CENK=0.500*(V_PROT+V_NEUT)
C_______________________________________________________________________
C
C
      VKACEN_PROTON=V0CENK*(1+XKAPPA_V0CENK*FRACTN)
C
      VKACEN_NEUTRS=V0CENK*(1-XKAPPA_V0CENK*FRACTN)
C
      IF (IZ_AUX.EQ.IN_AUX) THEN
          XKAPPA_V0CENK=0.0
          VKACEN_PROTON=V_PROT
          VKACEN_NEUTRS=V_NEUT
      END IF
C
C=======================================================================
C
      IF (IFK_RC.EQ.1) THEN
C
          IF (R_PROT.GT.R_NEUT) THEN
C
              XKAPPA_R0CENK= (R_PROT-R_NEUT)*(IN_AUX+IZ_AUX)
     *                     /((R_PROT+R_NEUT)*(IN_AUX-IZ_AUX))
C
              R0CENK=0.500*(R_NEUT+R_PROT)
C_______________________________________________________________________
C
              RKACEN_PROTON=R0CENK*(1+XKAPPA_R0CENK*FRACTN)
     *                            *A_FIXX**(1./3.)
C
              RKACEN_NEUTRS=R0CENK*(1-XKAPPA_R0CENK*FRACTN)
     *                            *A_FIXX**(1./3.)
C
          ELSE
C
              XKAPPA_R0CENK= (R_NEUT-R_PROT)*(IN_AUX+IZ_AUX)
     *                     /((R_PROT+R_NEUT)*(IN_AUX-IZ_AUX))
C
              R0CENK=0.500*(R_NEUT+R_PROT)
C_______________________________________________________________________
C
              RKACEN_PROTON=R0CENK*(1-XKAPPA_R0CENK*FRACTN)
     *                            *A_FIXX**(1./3.)
C
              RKACEN_NEUTRS=R0CENK*(1+XKAPPA_R0CENK*FRACTN)
     *                            *A_FIXX**(1./3.)
C
          END IF
C
          IF (IZ_AUX.EQ.IN_AUX) THEN
              XKAPPA_R0CENK=0.0
              RKACEN_PROTON=R_PROT
              RKACEN_NEUTRS=R_NEUT
          END IF
C
      END IF
C
C=======================================================================
C
      IF (IFK_AC.EQ.1) THEN
C
          IF (A_PROT.GT.A_NEUT) THEN
C
              XKAPPA_A0CENK=(A_PROT-A_NEUT)*(IN_AUX+IZ_AUX)
     *                    /((A_PROT+A_NEUT)*(IN_AUX-IZ_AUX))
C
              A0CENK=0.500*(A_NEUT+A_PROT)
C_______________________________________________________________________
C
              AKACEN_PROTON=A0CENK*(1+XKAPPA_A0CENK*FRACTN)
C
              AKACEN_NEUTRS=A0CENK*(1-XKAPPA_A0CENK*FRACTN)
C
          ELSE
C
              XKAPPA_A0CENK= (A_NEUT-A_PROT)*(IN_AUX+IZ_AUX)
     *                     /((A_PROT+A_NEUT)*(IN_AUX-IZ_AUX))
C
              A0CENK=0.500*(A_NEUT+A_PROT)
C_______________________________________________________________________
C
              AKACEN_PROTON=A0CENK*(1-XKAPPA_A0CENK*FRACTN)
C
              AKACEN_NEUTRS=A0CENK*(1+XKAPPA_A0CENK*FRACTN)
C
          END IF
C
          IF (IZ_AUX.EQ.IN_AUX) THEN
              XKAPPA_A0CENK=0.0
              AKACEN_PROTON=A_PROT
              AKACEN_NEUTRS=A_NEUT
          END IF
C
      END IF
C
C=======================================================================
C
      WRITE(NOUTPT,'(/,80(''#''),/,''#'',T80,''#'',/,
     *           ''#'',27X,''PARAMETERS WITH KAPPAS'',T80,''#'',/,
     *           ''#'',T80,''#'',/,
     *           ''#'',25X,''IZ_FIX= '',I3,4X,''IN_FIX= '',I3,4X,/,
     *           ''#'',T80,''#'',/,80(''#''))') IZ_FIX,IN_FIX
C_______________________________________________________________________
C      
      WRITE(NOUTPT,'(/,80(''#''),/,''#'',''  Protons - Central'',
     *             T80,''#'',/,''#  '',17X,''V0CENT  XKAP_V  VKACEN'',
     *                        3X,''R0CENT  XKAP_R  RKACEN'',T80,''#'',/,
     *               ''#  '',15x,3f8.3,1x,3F8.3,T80,''#'',/,
     *             ''#'',T80,''#'',/,
     *             ''#'',44X,''A0CENT  XKAP_A  AKACEN'',T80,''#'',/,
     *             ''#'',42x,3F8.3,T80,''#'',/,
     *             ''#'',T80,''#'',/,80(''#''))')
     *
     *            V0CENK,XKAPPA_V0CENK,VKACEN_PROTON,
     *            R0CENK,XKAPPA_R0CENK,RKACEN_PROTON,
     *            A0CENK,XKAPPA_A0CENK,AKACEN_PROTON
C
      WRITE(NOUTPT,'(/,80(''#''),/,''#'',''  Neutron - Central'',
     *             T80,''#'',/,''#  '',17X,''V0CENT  XKAP_V  VKACEN'',
     *                        3X,''R0CENT  XKAP_R  RKACEN'',T80,''#'',/,
     *               ''#  '',15x,3f8.3,1x,3F8.3,T80,''#'',/,
     *             ''#'',T80,''#'',/,
     *             ''#'',44X,''A0CENT  XKAP_A  AKACEN'',T80,''#'',/,
     *             ''#'',42x,3F8.3,T80,''#'',/,
     *             ''#'',T80,''#'',/,80(''#''))')
     *
     *            V0CENK,XKAPPA_V0CENK,VKACEN_NEUTRS,
     *            R0CENK,XKAPPA_R0CENK,RKACEN_NEUTRS,
     *            A0CENK,XKAPPA_A0CENK,AKACEN_NEUTRS
C
C=======================================================================
C=======================================================================
C
C     Introducing the "kappa" variant to the SPIN-ORBIT potential
C     parameters (Pure Woods-Saxon)
C
C     As they come from the minimization
C
      V_PROT=PARPOT(04)
      V_NEUT=PARPOT(24)
C
      R_PROT=PARPOT(05)
      R_NEUT=PARPOT(25)
C
      A_PROT=PARPOT(06)
      A_NEUT=PARPOT(26)
C
C=======================================================================
C
      XKAPPA_V0SORK= (V_PROT-V_NEUT)*(IN_AUX+IZ_AUX)
     *             /((V_PROT+V_NEUT)*(IN_AUX-IZ_AUX))
C
      V0SORK=0.500*(V_NEUT+V_PROT)
C_______________________________________________________________________
C
C     Protons
C
      VKASOR_PROTON=V0SORK*(1+XKAPPA_V0SORK*FRACTN)
C
C     Neutrons
C
      VKASOR_NEUTRS=V0SORK*(1-XKAPPA_V0SORK*FRACTN)
C
      IF (IZ_AUX.EQ.IN_AUX) THEN
          XKAPPA_V0SORK=0.0
          VKASOR_PROTON=V_PROT
          VKASOR_NEUTRS=V_NEUT
      END IF
C
C=======================================================================
C
      IF (IFK_RS.EQ.1) THEN
C
          IF (R_PROT.GT.R_NEUT) THEN
C
              XKAPPA_R0SORK= (R_PROT-R_NEUT)*(IN_AUX+IZ_AUX)
     *                     /((R_PROT+R_NEUT)*(IN_AUX-IZ_AUX))
C
              R0SORK=0.500*(R_NEUT+R_PROT)
C_______________________________________________________________________
C
              RKASOR_PROTON=R0SORK*(1+XKAPPA_R0SORK*FRACTN)
     *                            *A_FIXX**(1./3.)
C
              RKASOR_NEUTRS=R0SORK*(1-XKAPPA_R0SORK*FRACTN)
     *                            *A_FIXX**(1./3.)
C
          ELSE
C
              XKAPPA_R0SORK= (R_NEUT-R_PROT)*(IN_AUX+IZ_AUX)
     *                     /((R_PROT+R_NEUT)*(IN_AUX-IZ_AUX))
C
              R0SORK=0.500*(R_NEUT+R_PROT)
C_______________________________________________________________________
C
              RKASOR_PROTON=R0SORK*(1-XKAPPA_R0SORK*FRACTN)
     *                            *A_FIXX**(1./3.)
C
              RKASOR_NEUTRS=R0SORK*(1+XKAPPA_R0SORK*FRACTN)
     *                            *A_FIXX**(1./3.)
C
          END IF
C
          IF (IZ_AUX.EQ.IN_AUX) THEN
              XKAPPA_R0SORK=0.0
              RKASOR_PROTON=R_PROT
              RKASOR_NEUTRS=R_NEUT
          END IF
C
      END IF
C
C=======================================================================
C
      IF (IFK_AS.EQ.1) THEN
C
          IF (A_PROT.GT.A_NEUT) THEN
C
              XKAPPA_A0SORK= (A_PROT-A_NEUT)*(IN_AUX+IZ_AUX)
     *                     /((A_PROT+A_NEUT)*(IN_AUX-IZ_AUX))
C
              A0SORK=0.500*(A_NEUT+A_PROT)
C_______________________________________________________________________
C
              AKASOR_PROTON=A0SORK*(1+XKAPPA_A0SORK*FRACTN)
C
              AKASOR_NEUTRS=A0SORK*(1-XKAPPA_A0SORK*FRACTN)
C
          ELSE
C
              XKAPPA_A0SORK= (A_NEUT-A_PROT)*(IN_AUX+IZ_AUX)
     *                     /((A_PROT+A_NEUT)*(IN_AUX-IZ_AUX))
C
              A0SORK=0.500*(A_NEUT+A_PROT)
C_______________________________________________________________________
C
              AKASOR_PROTON=A0SORK*(1-XKAPPA_A0SORK*FRACTN)
C
              AKASOR_NEUTRS=A0SORK*(1+XKAPPA_A0SORK*FRACTN)
C
          END IF
C
          IF (IZ_AUX.EQ.IN_AUX) THEN
              XKAPPA_A0SORK=0.0
              AKASOR_PROTON=A_PROT
              AKASOR_NEUTRS=A_NEUT
          END IF
C
      END IF
C
C=======================================================================
C      
      WRITE(NOUTPT,'(/,80(''#''),/,''#'',''  Protons - Central'',
     *             T80,''#'',/,''#  '',17X,''V0SORB  XKAP_V  VKASOR'',
     *                        3X,''R0SORB  XKAP_R  RKASOR'',T80,''#'',/,
     *               ''#  '',15x,3f8.3,1x,3F8.3,T80,''#'',/,
     *             ''#'',T80,''#'',/,
     *             ''#'',44X,''A0SORB  XKAP_A  AKASOR'',T80,''#'',/,
     *             ''#'',42x,3F8.3,T80,''#'',/,
     *             ''#'',T80,''#'',/,80(''#''))')
     *
     *            V0SORK,XKAPPA_V0SORK,VKASOR_PROTON,
     *            R0SORK,XKAPPA_R0SORK,RKASOR_PROTON,
     *            A0SORK,XKAPPA_A0SORK,AKASOR_PROTON
C
      WRITE(NOUTPT,'(/,80(''#''),/,''#'',''  Neutron - Central'',
     *             T80,''#'',/,''#  '',17X,''V0SORB  XKAP_V  VKASOR'',
     *                        3X,''R0SORB  XKAP_R  RKASOR'',T80,''#'',/,
     *               ''#  '',15x,3f8.3,1x,3F8.3,T80,''#'',/,
     *             ''#'',T80,''#'',/,
     *             ''#'',44X,''A0SORB  XKAP_A  AKASOR'',T80,''#'',/,
     *             ''#'',42x,3F8.3,T80,''#'',/,
     *             ''#'',T80,''#'',/,80(''#''))')
     *
     *            V0SORK,XKAPPA_V0SORK,VKASOR_NEUTRS,
     *            R0SORK,XKAPPA_R0SORK,RKASOR_NEUTRS,
     *            A0SORK,XKAPPA_A0SORK,AKASOR_NEUTRS
C
C=======================================================================
C
      RETURN
      END
C
C=======================================================================
C=======================================================================
C
      FUNCTION VCENTR(INUCLI,ZARGUM,IPOINT,LQNUMB)
C      
      INCLUDE    'MATDIM/NDNUCL.f'
C
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
     *        /HBAR_V/ HBAR_C
      COMMON
     *       /DATA_V/ AOSCIL_PROTON(1:NDNUCL),
     *                AOSCIL_NEUTRS(1:NDNUCL),
     *                HOMEGA(1:NDNUCL)      
      COMMON
     *       /NUCLEU/ IZ_FIX,IN_FIX,ISOSPI 
      COMMON
     *       /VPROEF/ V0EFFC
      COMMON
     *       /TENSOR/ IFTENS
      COMMON
     *       /TENOPT/ ICENTT,ISORBT,ITENSR 
C
      DATA
     *        ALFCON / 137.03602 /     
C
C=======================================================================
C
C      IF (LOGWRI.GT.0) THEN
C          WRITE(LOGFIL,'(/,9X,''Entering VCENTR'')')
C      END IF
C
C=======================================================================
C      
      IF (ISOSPI.EQ.1) THEN
          AOSCIL=AOSCIL_PROTON(INUCLI)
	      VKACEN=VKACEN_PROTON
	      RKACEN=RKACEN_PROTON
	      AKACEN=AKACEN_PROTON
      END IF
C
      IF (ISOSPI.EQ.0) THEN
          AOSCIL=AOSCIL_NEUTRS(INUCLI)
	      VKACEN=VKACEN_NEUTRS
	      RKACEN=RKACEN_NEUTRS
	      AKACEN=AKACEN_NEUTRS
      END IF
C
C=======================================================================
C
      RADIUS=AOSCIL*SQRT(ZARGUM) ! Radius in Fermi
C
C=======================================================================
C
      VCENWS=VKACEN/(1.0+EXP((RADIUS-RKACEN)/AKACEN))
C
C=======================================================================
C
C     Analytical exression for the Coulomb potential of the uniformly
C                                                      charged sphere
      V_COUL=0.0
C
      IF (ISOSPI.EQ.1) THEN
C            
          IF (RADIUS.LE.RKACOU) THEN
C            
              V_COUL=HBAR_C*(IZ_FIX-1)/ALFCON/RKACOU
     *              *(1.500-0.5000000*(RADIUS/RKACOU)**2)   
          ELSE
              V_COUL=HBAR_C*(IZ_FIX-1)/(ALFCON*RADIUS)
C            
          END IF
C          
          V0EFFC=VKACEN+HBAR_C*(IZ_FIX-1)/ALFCON/RKACOU*1.500
C
      END IF
C
C=======================================================================
C @@@
C     The tensor contribution to the central potential, cf. Dx Eq.y
C
      VCTENS=0.0
C      
      IF (IFTENS.EQ.1 .AND. ICENTT.EQ.1) THEN
C
          SPIN_J_PROTON=SPIN_J(INUCLI,1,ZARGUM,IPOINT,LQNUMB)
          SPIN_J_NEUTRS=SPIN_J(INUCLI,0,ZARGUM,IPOINT,LQNUMB)
C          
          IF (ISOSPI.EQ.1) THEN
C              
              VCTENS=CLAMPP*SPIN_J_PROTON+CLAMPN*SPIN_J_NEUTRS
C              
          END IF
C          
          IF (ISOSPI.EQ.0) THEN
C              
              VCTENS=CLAMNP*SPIN_J_PROTON+CLAMNN*SPIN_J_NEUTRS
C              
          END IF
C          
      END IF
C
C=======================================================================    
C
      VCENTR=VCENWS+V_COUL+VCTENS
C
C=======================================================================
C
C      IF (LOGWRI.GT.0) THEN
C          WRITE(LOGFIL,'(/,9X,''Exiting VCENTR'')')
C      END IF
C
C=======================================================================
C
      RETURN
      END      
C
C=======================================================================
C=======================================================================
C
      FUNCTION V_SORB(INUCLI,ZARGUM)
C
      INCLUDE    'MATDIM/NDNUCL.f'
C
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
     *       /DATA_V/ AOSCIL_PROTON(1:NDNUCL),
     *                AOSCIL_NEUTRS(1:NDNUCL),
     *                HOMEGA(1:NDNUCL)      
      COMMON
     *       /NUCLEU/ IZ_FIX,IN_FIX,ISOSPI  
C
C=======================================================================
C     Spherical spin-orbit Woods-Saxon form-factor;
C     the argument denoted ZARGUM is dimensionless
C=======================================================================
C
C      IF (LOGWRI.GT.0) THEN
C          WRITE(LOGFIL,'(/,9X,''Exiting V_SORB'')')
C      END IF
C
C=======================================================================
C
      IF (ISOSPI.EQ.1) THEN
C
          AOSCIL=AOSCIL_PROTON(INUCLI)
	  VKASOR=VKASOR_PROTON
          AKASOR=AKASOR_PROTON
          RKASOR=RKASOR_PROTON
          V0UNIT=V0UNIT_PROTON
C
      END IF
C
      IF (ISOSPI.EQ.0) THEN
C
          AOSCIL=AOSCIL_NEUTRS(INUCLI)
          VKASOR=VKASOR_NEUTRS
          AKASOR=AKASOR_NEUTRS
          RKASOR=RKASOR_NEUTRS
          V0UNIT=V0UNIT_NEUTRS
C
      END IF
C
C=======================================================================
C      
      RADIUS=AOSCIL*SQRT(ZARGUM) ! Radius in Fermis
C      WRITE(0,*)'RADIUS',RADIUS,AOSCIL,ZARGUM,INUCLI
C
C=======================================================================
C      
      V_SORB=-VKASOR/AKASOR*EXP((RADIUS-RKASOR)/AKASOR)
     *              /(1.000+EXP((RADIUS-RKASOR)/AKASOR))**2
     *              /RADIUS
C
C=======================================================================
C
C      IF (LOGWRI.GT.0) THEN
C          WRITE(LOGFIL,'(/,9X,''Exiting V_SORB'')')
C      END IF
C
C=======================================================================
C
      RETURN
      END         
C      
C=======================================================================
C=======================================================================
C
      FUNCTION V_SORB_DENSIT(INUCLI,ZARGUM,IPOINT,LQNUMB)
          
      INCLUDE 'MATDIM/NDNUCL.f'
C      
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
     *       /NUCLEU/ IZ_FIX,IN_FIX,ISOSPI  
      COMMON
     *       /TENSOR/ IFTENS
      COMMON
     *       /TENOPT/ ICENTT,ISORBT,ITENSR
      COMMON
     *       /DENSIF/ ITRNUM
      COMMON
     *       /DATA_V/ AOSCIL_PROTON(1:NDNUCL),
     *                AOSCIL_NEUTRS(1:NDNUCL),
     *                HOMEGA(1:NDNUCL)
C      
C=======================================================================
C     This FUNCTION computes the density dependent spin-orbit potential 
C     taking into account the spin-J current (optional: IFTENS=1)
C                                ===>> Eq.(5.1.39) of MFV_14D.pdf
C=======================================================================
C      
C      IF (LOGWRI.GT.0) THEN
C          WRITE(0,'(/,9X,''Entering V_SORB_DENSIT'')')
C      END IF
C
C=======================================================================
C
      CALL CPUTIM('V_SORB',1)
C
C=======================================================================
C     
      IF (ISOSPI.EQ.1) THEN ! Protons
C          
          VSOAUX=DENSIT_GRADNT(INUCLI,1,ZARGUM,IPOINT,LQNUMB)*ALAMPP
     *          +DENSIT_GRADNT(INUCLI,0,ZARGUM,IPOINT,LQNUMB)*ALAMPN
C     
          IF (IFTENS.EQ.1 .AND. ISORBT.EQ.1) THEN
C              
              VSOAUX=VSOAUX
     *              +SPINJR(INUCLI,1,ZARGUM,IPOINT,LQNUMB)*TLAMPP
     *              +SPINJR(INUCLI,0,ZARGUM,IPOINT,LQNUMB)*TLAMPN
C     
          END IF
C      
      END IF
C      
      IF (ISOSPI.EQ.0) THEN ! Neutrons
C       
          VSOAUX=DENSIT_GRADNT(INUCLI,1,ZARGUM,IPOINT,LQNUMB)*ALAMNP
     *          +DENSIT_GRADNT(INUCLI,0,ZARGUM,IPOINT,LQNUMB)*ALAMNN
     
          IF (IFTENS.EQ.1 .AND. ISORBT.EQ.1) THEN
C              
              VSOAUX=VSOAUX
     *              +SPINJR(INUCLI,1,ZARGUM,IPOINT,LQNUMB)*TLAMNP
     *              +SPINJR(INUCLI,0,ZARGUM,IPOINT,LQNUMB)*TLAMNN
C     
          END IF
     
      END IF
C      
      V_SORB_DENSIT=VSOAUX
C
C=======================================================================
C      
C      IF (LOGWRI.GT.0) THEN
C          WRITE(0,'(/,9X,''Exiting V_SORB_DENSIT'')')
C      END IF
C
C=======================================================================
C 
      CALL CPUTIM('V_SORB',0)
C
C=======================================================================
C      
      RETURN
      END	     
C      
C=======================================================================
C=======================================================================
C
      FUNCTION DENSIT(XARGUM)
C      
      INCLUDE  'MATDIM/NDSPEC.f'
      INCLUDE  'MATDIM/NDIM_N.f'
      INCLUDE  'MATDIM/NDGAUS.f'
      INCLUDE  'MATDIM/NDNUCL.f'
C      
      PARAMETER
     *         (NDIM_L=NDIM_N,NDBASE=NDIM_N+1) 
      REAL*16
     *          QARGUM,QFLAGN,QDFLAG,QNORNL
      DIMENSION
     *          LABORD(0:NDIM_N,0:NDIM_L,0:2*NDIM_L+1)
      DIMENSION
     *          CMATUP(0:NDIM_L,1:NDBASE,1:NDBASE),
     *          CMATDN(0:NDIM_L,1:NDBASE,1:NDBASE)
      DIMENSION
     *          ENERUP(0:NDIM_L,1:NDBASE),
     *          ENERDN(0:NDIM_L,1:NDBASE)
      DIMENSION
     *          VCOEUP(0:NDIM_L,1:NDBASE),
     *          VCOEDN(0:NDIM_L,1:NDBASE)
      DIMENSION
     *          NWSSPH(1:NDSPEC),
     *          LWSSPH(1:NDSPEC),
     *          JWSSPH(1:NDSPEC)
      DIMENSION
     *          DENSIT_VECTOR(1:NDGAUS,0:NDIM_L),
     *          DENAUX_TEST_L(0:NDIM_L)
      DIMENSION
     *          QFLAGN(0:NDIM_N,0:NDIM_L),
     *          QDFLAG(0:NDIM_N,0:NDIM_L),
     *          XFLAGN(0:NDIM_N,0:NDIM_L),
     *          XDFLAG(0:NDIM_N,0:NDIM_L)
C
      COMMON
     *       /DATA_V/ AOSCIL_PROTON(1:NDNUCL),
     *                AOSCIL_NEUTRS(1:NDNUCL),
     *                HOMEGA(1:NDNUCL)   
      COMMON
     *       /GAUSHE/ NGAUSS,NSHELL_PROTON,NSHELL_NEUTRS
      COMMON
     *       /NMLCNT/ XNORNL(0:NDIM_N,0:NDIM_L)
      COMMON
     *       /POLDEN/ 
     *        XPOLYN_DENSIT(1:NDGAUS,0:NDIM_N,0:NDIM_L,0:NDIM_L), 
     *        XDPOLY_DENSIT(1:NDGAUS,0:NDIM_N,0:NDIM_L,0:NDIM_L)
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
     *       /DIMSPH/ LDSING_PROTON,LEVTHE_PROTON,
     *                LEVEXP_PROTON(1:NDNUCL),
     *                LDSING_NEUTRS,LEVTHE_NEUTRS,
     *                LEVEXP_NEUTRS(1:NDNUCL)
      COMMON
     *       /QUANTN/ NWSSPH_PROTON(1:NDSPEC),
     *      	      LWSSPH_PROTON(1:NDSPEC),
     *      	      JWSSPH_PROTON(1:NDSPEC),
     *
     *      	      NWSSPH_NEUTRS(1:NDSPEC),
     *      	      LWSSPH_NEUTRS(1:NDSPEC),
     *      	      JWSSPH_NEUTRS(1:NDSPEC)
       COMMON
     *       /LEVORD/ LABORD_PROTON(1:NDNUCL,
     *                              0:NDIM_N,0:NDIM_L,0:2*NDIM_L+1),
     *                LABORD_NEUTRS(1:NDNUCL,
     *                              0:NDIM_N,0:NDIM_L,0:2*NDIM_L+1)
      COMMON
     *       /NUCLEU/ IZ_FIX,IN_FIX,ISOSPI 
      COMMON
     *       /ARGUMN/ IARG_R,IARG_Z
      COMMON
     *       /IFPAIR/ IF_PAI,IFGPAI,IFDPAI
      COMMON
     *       /BCSENE/ ENEMAX_BCSPAI(1:NDNUCL,1:2)
      COMMON
     *       /BCSDEL/ DELT_2(1:NDNUCL,1:2)
      COMMON
     *       /BCSLAM/ XLAMBD_PROTON(1:NDNUCL),
     *                XLAMBD_NEUTRS(1:NDNUCL)
      COMMON
     *       /SIMUSE/ AOSCIL,RADCUT,INUCLI
C
C=======================================================================   
C  
C       This FUNCTION calculates the TOTAL DENSITY FUNCTION D1. Eq.(4.7.14)
C          
C          It is prepared to distinct between neutrons and protons
C
C=======================================================================     
C
C             It is mainly used for plotting purposes
C
C=======================================================================     
C
C      WRITE(0,'(''Entering DENSIT with INUCLI= '',I3,''ISOSPI= '',I3)')
C     *                                 INUCLI,ISOSPI
C
C=======================================================================     
C
C                        P R O T O N S 
C
C=======================================================================    
C      
      IF (ISOSPI.EQ.1) THEN
C        
          N_PART=IZ_FIX
          NSHELL=NSHELL_PROTON
          LDSING=LDSING_PROTON
C        
          AOSCIL=AOSCIL_PROTON(INUCLI)
          
          IF (IF_PAI.EQ.1) THEN
              XLAMBD=XLAMBD_PROTON(INUCLI)
              DELTA2=DELT_2(INUCLI,1)
              ENEMAX=ENEMAX_BCSPAI(INUCLI,1)
          END IF
C        
          DO INDEXX=1,LDSING
C          
             LWSSPH(INDEXX)=LWSSPH_PROTON(INDEXX)
             NWSSPH(INDEXX)=NWSSPH_PROTON(INDEXX)
             JWSSPH(INDEXX)=JWSSPH_PROTON(INDEXX)
C          
          END DO
C        
          DO INDEXX=1,LDSING
C
             JQNAUX=JWSSPH(INDEXX)
             LQNAUX=LWSSPH(INDEXX)
	     NNBAUX=NWSSPH(INDEXX)
             NBBASE=(NSHELL-LQNAUX)/2  
C             
             LABORD(NNBAUX,LQNAUX,JQNAUX)
     *	    =LABORD_PROTON(INUCLI,NNBAUX,LQNAUX,JQNAUX)
C     
             ENERUP(LQNAUX,NNBAUX+1)
     *      =ENERUP_AUXPRO(LQNAUX,NNBAUX+1)          
C
             ENERDN(LQNAUX,NNBAUX+1)
     *      =ENERDN_AUXPRO(LQNAUX,NNBAUX+1)
C
             DO N1NUMB=0,NBBASE
C		  
                CMATUP(LQNAUX,N1NUMB+1,NNBAUX+1)
     *         =CMATUP_AUXPRO(LQNAUX,N1NUMB+1,NNBAUX+1)
C		  
                CMATDN(LQNAUX,N1NUMB+1,NNBAUX+1)
     *         =CMATDN_AUXPRO(LQNAUX,N1NUMB+1,NNBAUX+1)
C
             END DO
         END DO                  
C
      END IF
C
C=======================================================================     
C
C                        N E U T R O N S 
C
C=======================================================================          
C      
      IF (ISOSPI.EQ.0) THEN
C        
          N_PART=IN_FIX
          NSHELL=NSHELL_NEUTRS
          LDSING=LDSING_NEUTRS
C        
          AOSCIL=AOSCIL_NEUTRS(INUCLI)
          
          IF (IF_PAI.EQ.1) THEN
              XLAMBD=XLAMBD_NEUTRS(INUCLI)
              DELTA2=DELT_2(INUCLI,2)
              ENEMAX=ENEMAX_BCSPAI(INUCLI,2)
          END IF
C        
          DO INDEXX=1,LDSING
C          
             LWSSPH(INDEXX)=LWSSPH_NEUTRS(INDEXX)
             NWSSPH(INDEXX)=NWSSPH_NEUTRS(INDEXX)
             JWSSPH(INDEXX)=JWSSPH_NEUTRS(INDEXX)
C          
          END DO
C        
          DO INDEXX=1,LDSING
C
             JQNAUX=JWSSPH(INDEXX)
             LQNAUX=LWSSPH(INDEXX)
	     NNBAUX=NWSSPH(INDEXX)
             NBBASE=(NSHELL-LQNAUX)/2  
C             
             LABORD(NNBAUX,LQNAUX,JQNAUX)
     *	    =LABORD_NEUTRS(INUCLI,NNBAUX,LQNAUX,JQNAUX)
     
             ENERUP(LQNAUX,NNBAUX+1)
     *      =ENERUP_AUXNEU(LQNAUX,NNBAUX+1)          
C
             ENERDN(LQNAUX,NNBAUX+1)
     *      =ENERDN_AUXNEU(LQNAUX,NNBAUX+1)
C
             DO N1NUMB=0,NBBASE
C		  
                CMATUP(LQNAUX,N1NUMB+1,NNBAUX+1)
     *         =CMATUP_AUXNEU(LQNAUX,N1NUMB+1,NNBAUX+1)
C		  
                CMATDN(LQNAUX,N1NUMB+1,NNBAUX+1)
     *         =CMATDN_AUXNEU(LQNAUX,N1NUMB+1,NNBAUX+1)
     
C
             END DO
         END DO                  
C
      END IF
C
C=======================================================================     
C
C       Calculating the polynomials in the desired ZARGUM point
C
C=======================================================================
C
C     IARG_R and IARG_Z are used in order to choose if we want to 
C     integrate over small r or over z.
C
C     IARG_R=1 - the integration will be done over r
C
C     IARG_Z=1 - the integration will be done over z
C
C=======================================================================
C
      IF (IARG_Z.EQ.1) ZARGUM=XARGUM
      IF (IARG_R.EQ.1) ZARGUM=(XARGUM/AOSCIL)**2
C
      KIND_L=1
      NORDER=NSHELL/2
      LORDER=NSHELL
C      
      QARGUM=QEXT(ZARGUM)    !passing the argument to real*16
C      
      CALL LAGUER(QARGUM,QFLAGN,QDFLAG,NORDER,LORDER,KIND_L,
     *                                        NDIM_L,NDIM_N)
C     
      DO NAUXIL=0,NORDER
        DO LAUXIL=0,LORDER
            XFLAGN(NAUXIL,LAUXIL)=REAL(QFLAGN(NAUXIL,LAUXIL))
            XDFLAG(NAUXIL,LAUXIL)=REAL(QDFLAG(NAUXIL,LAUXIL))
        END DO
      END DO
C
C=======================================================================     
C
C                   Calculating the DENSITY Eq.(4.7.14)
C
C=======================================================================        
C      
      PINUMB=4.0*ATAN(1.0)
      FACTOR=1.0/(PINUMB*(AOSCIL**3)) ! multiplied by 2/a 
C                                             from N_{n1,l}*N_{n2,l}
C
C  Here we calculate J_UP block
C      
      DEN_UP=0.0
C
      TESTV2=0.0
C      
      DO LQPRIM=NSHELL,0,-1
C          
         JNUMUP=2*LQPRIM+1
C
         DO NWSNUM=1,(NSHELL-LQPRIM)/2+1
        
            ENERGY=ENERUP(LQPRIM,NWSNUM)
           
            IF (IF_PAI.EQ.1.AND.ENERGY.GT.ENEMAX) GO TO 1
           
            IF (IF_PAI.EQ.0
     *          .AND.LABORD(NWSNUM-1,LQPRIM,JNUMUP).EQ.0) GO TO 1
     
            DNL_UP=0.0
C            
            DO N1NUMB=0,(NSHELL-LQPRIM)/2
C                
               DO N2NUMB=0,(NSHELL-LQPRIM)/2
C                  
                  CFACT1=(LQPRIM+1)*CMATUP(LQPRIM,N1NUMB+1,NWSNUM)
     *                             *CMATUP(LQPRIM,N2NUMB+1,NWSNUM)
C    
                  POLPOL=XNORNL(N1NUMB,LQPRIM)*XFLAGN(N1NUMB,LQPRIM)
     *                  *XNORNL(N2NUMB,LQPRIM)*XFLAGN(N2NUMB,LQPRIM)
C                
                  DNL_UP=DNL_UP+ZARGUM**LQPRIM*EXP(-ZARGUM)
     *                         *POLPOL*CFACT1
C                  
               END DO !N2NUMB
C              
            END DO !N1NUMB
C
                IF (IF_PAI.EQ.1) THEN
                    V2_PUP=V2PAIR(ENERGY,XLAMBD,DELTA2)
                    DEN_UP=DEN_UP+DNL_UP*V2_PUP
                    TESTV2=TESTV2+V2_PUP*(JNUMUP+1)
                ELSE
                    DEN_UP=DEN_UP+DNL_UP
                END IF
C
   1         CONTINUE
C            
         END DO !NWSNUM
C        
      END DO !LQPRIM
C
C  Here we calculate J_DOWN block
C
      DEN_DW=0.0

      DO LQPRIM=NSHELL,0,-1
C
         JNUMDW=2*LQPRIM-1     
C          
         DO NWSNUM=1,(NSHELL-LQPRIM)/2+1
            
            ENERGY=ENERDN(LQPRIM,NWSNUM)
           
            IF (JNUMDW.GE.0) THEN
                
                IF (IF_PAI.EQ.1.AND.ENERGY.GT.ENEMAX) GO TO 2
           
                IF (IF_PAI.EQ.0
     *             .AND.LABORD(NWSNUM-1,LQPRIM,JNUMDW).EQ.0) GO TO 2
            
                DNL_DW=0.0
C            
                DO N1NUMB=0,(NSHELL-LQPRIM)/2
C                
                   DO N2NUMB=0,(NSHELL-LQPRIM)/2
C
                      CFACT2=(LQPRIM  )*CMATDN(LQPRIM,N1NUMB+1,NWSNUM)
     *                                 *CMATDN(LQPRIM,N2NUMB+1,NWSNUM)
C    
                      POLPOL=XNORNL(N1NUMB,LQPRIM)*XFLAGN(N1NUMB,LQPRIM)
     *                      *XNORNL(N2NUMB,LQPRIM)*XFLAGN(N2NUMB,LQPRIM)
C                
                      DNL_DW=DNL_DW+ZARGUM**LQPRIM*EXP(-ZARGUM)
     *                             *POLPOL*CFACT2
C                  
                   END DO !N2NUMB
C              
                END DO !N1NUMB
               
                IF (IF_PAI.EQ.1) THEN
                    V2_PDW=V2PAIR(ENERGY,XLAMBD,DELTA2)
                    DEN_DW=DEN_DW+DNL_DW*V2_PDW
                    TESTV2=TESTV2+V2_PDW*(JNUMDW+1)
                ELSE
                    DEN_DW=DEN_DW+DNL_DW
                END IF             
                                                
            END IF  !JNUMDW.GE.0
           
   2        CONTINUE            
      
         END DO !NWSNUM
C        
      END DO !LQPRIM
C      
      DENSIT=(DEN_UP+DEN_DW)*FACTOR
C
C=======================================================================
C     
      RETURN
      END
C      
C=======================================================================
C=======================================================================
C
      FUNCTION DENSIT_LAGUER(INUCLI,ISOSPI_AUXILI,ZARGUM,IPOINT,LQNUMB)
C      
      INCLUDE  'MATDIM/NDSPEC.f'
      INCLUDE  'MATDIM/NDIM_N.f'
      INCLUDE  'MATDIM/NDGAUS.f'
      INCLUDE  'MATDIM/NDNUCL.f'
C      
      PARAMETER
     *         (NDIM_L=NDIM_N,NDBASE=NDIM_N+1) 
      DIMENSION
     *          LABORD(0:NDIM_N,0:NDIM_L,0:2*NDIM_L+1)
      DIMENSION
     *          CMATUP(0:NDIM_L,1:NDBASE,1:NDBASE),
     *          CMATDN(0:NDIM_L,1:NDBASE,1:NDBASE)
      DIMENSION
     *          ENERUP(0:NDIM_L,1:NDBASE),
     *          ENERDN(0:NDIM_L,1:NDBASE)
      DIMENSION
     *          VCOEUP(0:NDIM_L,1:NDBASE),
     *          VCOEDN(0:NDIM_L,1:NDBASE)
      DIMENSION
     *          NWSSPH(1:NDSPEC),
     *          LWSSPH(1:NDSPEC),
     *          JWSSPH(1:NDSPEC)
      DIMENSION
     *          DENSIT_VECTOR(1:NDGAUS,0:NDIM_L),
     *          DENAUX_TEST_L(0:NDIM_L)
C
      COMMON
     *       /DATA_V/ AOSCIL_PROTON(1:NDNUCL),
     *                AOSCIL_NEUTRS(1:NDNUCL),
     *                HOMEGA(1:NDNUCL)   
      COMMON
     *       /GAUSHE/ NGAUSS,NSHELL_PROTON,NSHELL_NEUTRS
      COMMON
     *       /POLPOL/ XPOLYN(1:NDGAUS,0:NDIM_N,0:NDIM_L), 
     *                XDPOLY(1:NDGAUS,0:NDIM_N,0:NDIM_L)
      COMMON
     *       /POLDEN/ 
     *        XPOLYN_DENSIT(1:NDGAUS,0:NDIM_N,0:NDIM_L,0:NDIM_L), 
     *        XDPOLY_DENSIT(1:NDGAUS,0:NDIM_N,0:NDIM_L,0:NDIM_L)
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
     *       /DIMSPH/ LDSING_PROTON,LEVTHE_PROTON,
     *                LEVEXP_PROTON(1:NDNUCL),
     *                LDSING_NEUTRS,LEVTHE_NEUTRS,
     *                LEVEXP_NEUTRS(1:NDNUCL)
      COMMON
     *       /QUANTN/ NWSSPH_PROTON(1:NDSPEC),
     *      	      LWSSPH_PROTON(1:NDSPEC),
     *      	      JWSSPH_PROTON(1:NDSPEC),
     *
     *      	      NWSSPH_NEUTRS(1:NDSPEC),
     *      	      LWSSPH_NEUTRS(1:NDSPEC),
     *      	      JWSSPH_NEUTRS(1:NDSPEC)
       COMMON
     *       /LEVORD/ LABORD_PROTON(1:NDNUCL,
     *                              0:NDIM_N,0:NDIM_L,0:2*NDIM_L+1),
     *                LABORD_NEUTRS(1:NDNUCL,
     *                              0:NDIM_N,0:NDIM_L,0:2*NDIM_L+1)
      COMMON
     *       /NUCLEU/ IZ_FIX,IN_FIX,ISOSPI 
      COMMON
     *       /IFPAIR/ IF_PAI,IFGPAI,IFDPAI
      COMMON
     *       /BCSENE/ ENEMAX_BCSPAI(1:NDNUCL,1:2)
      COMMON
     *       /BCSDEL/ DELT_2(1:NDNUCL,1:2)
      COMMON
     *       /BCSLAM/ XLAMBD_PROTON(1:NDNUCL),
     *                XLAMBD_NEUTRS(1:NDNUCL)
C
C=======================================================================   
C  
C       This FUNCTION calculates the TOTAL DENSITY FUNCTION Eq.(4.7.14)
C       at LAGUERRE NODES
C          
C          It is prepared to distinct between netrons and protons 
C
C=======================================================================     
C
C                        P R O T O N S 
C
C=======================================================================    
C     
      IF (ISOSPI_AUXILI.EQ.1) THEN
C        
          N_PART=IZ_FIX
          NSHELL=NSHELL_PROTON
          LDSING=LDSING_PROTON
C        
          AOSCIL=AOSCIL_PROTON(INUCLI)
          
          IF (IF_PAI.EQ.1) THEN
              XLAMBD=XLAMBD_PROTON(INUCLI)
              DELTA2=DELT_2(INUCLI,1)
              ENEMAX=ENEMAX_BCSPAI(INUCLI,1)
          END IF
C        
          DO INDEXX=1,LDSING
C          
             LWSSPH(INDEXX)=LWSSPH_PROTON(INDEXX)
             NWSSPH(INDEXX)=NWSSPH_PROTON(INDEXX)
             JWSSPH(INDEXX)=JWSSPH_PROTON(INDEXX)
C          
          END DO
C        
          DO INDEXX=1,LDSING
C
             JQNAUX=JWSSPH(INDEXX)
             LQNAUX=LWSSPH(INDEXX)
	     NNBAUX=NWSSPH(INDEXX)
             NBBASE=(NSHELL-LQNAUX)/2  
C             
             LABORD(NNBAUX,LQNAUX,JQNAUX)
     *	    =LABORD_PROTON(INUCLI,NNBAUX,LQNAUX,JQNAUX)
C     
             ENERUP(LQNAUX,NNBAUX+1)
     *      =ENERUP_AUXPRO(LQNAUX,NNBAUX+1)          
C
             ENERDN(LQNAUX,NNBAUX+1)
     *      =ENERDN_AUXPRO(LQNAUX,NNBAUX+1)
C
             DO N1NUMB=0,NBBASE
C		  
                CMATUP(LQNAUX,N1NUMB+1,NNBAUX+1)
     *         =CMATUP_AUXPRO(LQNAUX,N1NUMB+1,NNBAUX+1)
C		  
                CMATDN(LQNAUX,N1NUMB+1,NNBAUX+1)
     *         =CMATDN_AUXPRO(LQNAUX,N1NUMB+1,NNBAUX+1)
C
             END DO
          END DO                  
C
      END IF
C
C=======================================================================     
C
C                        N E U T R O N S 
C
C=======================================================================          
C      
      IF (ISOSPI_AUXILI.EQ.0) THEN
        
          N_PART=IN_FIX
          NSHELL=NSHELL_NEUTRS
          LDSING=LDSING_NEUTRS
        
          AOSCIL=AOSCIL_NEUTRS(INUCLI)
          
          IF (IF_PAI.EQ.1) THEN
              XLAMBD=XLAMBD_NEUTRS(INUCLI)
              DELTA2=DELT_2(INUCLI,2)
              ENEMAX=ENEMAX_BCSPAI(INUCLI,2)
          END IF
        
          DO INDEXX=1,LDSING
          
             LWSSPH(INDEXX)=LWSSPH_NEUTRS(INDEXX)
             NWSSPH(INDEXX)=NWSSPH_NEUTRS(INDEXX)
             JWSSPH(INDEXX)=JWSSPH_NEUTRS(INDEXX)
          
          END DO
        
          DO INDEXX=1,LDSING
C
             JQNAUX=JWSSPH(INDEXX)
             LQNAUX=LWSSPH(INDEXX)
	     NNBAUX=NWSSPH(INDEXX)
             NBBASE=(NSHELL-LQNAUX)/2  
             
             LABORD(NNBAUX,LQNAUX,JQNAUX)
     *	    =LABORD_NEUTRS(INUCLI,NNBAUX,LQNAUX,JQNAUX)
     
             ENERUP(LQNAUX,NNBAUX+1)
     *      =ENERUP_AUXNEU(LQNAUX,NNBAUX+1)          
C
             ENERDN(LQNAUX,NNBAUX+1)
     *      =ENERDN_AUXNEU(LQNAUX,NNBAUX+1)
C
             DO N1NUMB=0,NBBASE
C		  
                CMATUP(LQNAUX,N1NUMB+1,NNBAUX+1)
     *         =CMATUP_AUXNEU(LQNAUX,N1NUMB+1,NNBAUX+1)
C		  
                CMATDN(LQNAUX,N1NUMB+1,NNBAUX+1)
     *         =CMATDN_AUXNEU(LQNAUX,N1NUMB+1,NNBAUX+1)
C
             END DO
          END DO                  
C
      END IF
C
C=======================================================================     
C
C                   Calculating the DENSITY Eq.(4.7.14)
C
C=======================================================================        
C      
      PINUMB=4.0*ATAN(1.0)
      FACTOR=1.0/(PINUMB*(AOSCIL**3)) ! multiplied by 2/a 
C                                             from N_{n1,l}*N_{n2,l}
C
C  Here we calculate J_UP block
C      
      DEN_UP=0.0
C
C      DO LQPRIM=NSHELL,0,-1
C          
        JNUMUP=2*LQNUMB+1
        
         DO NWSNUM=1,(NSHELL-LQNUMB)/2+1
        
            ENERGY=ENERUP(LQNUMB,NWSNUM)
           
            IF (IF_PAI.EQ.1.AND.ENERGY.GT.ENEMAX) GO TO 1
           
            IF (IF_PAI.EQ.0
     *          .AND.LABORD(NWSNUM-1,LQNUMB,JNUMUP).EQ.0) GO TO 1
     
                 DNL_UP=0.0
            
                 DO N1NUMB=0,(NSHELL-LQNUMB)/2
                
                    DO N2NUMB=0,(NSHELL-LQNUMB)/2
                  
                       CFACT1=(LQNUMB+1)*CMATUP(LQNUMB,N1NUMB+1,NWSNUM)
     *                                  *CMATUP(LQNUMB,N2NUMB+1,NWSNUM)
               
                      POLPOL=XPOLYN(IPOINT,N1NUMB,LQNUMB)
     *                      *XPOLYN(IPOINT,N2NUMB,LQNUMB)
    
C              POLPOL=XPOLYN_DENSIT(IPOINT,N1NUMB,LQPRIM,LQNUMB) 
C     *              *XPOLYN_DENSIT(IPOINT,N2NUMB,LQPRIM,LQNUMB)
                
                      DNL_UP=DNL_UP+POLPOL*CFACT1
CID     *                             *ZARGUM**LQNUMB*EXP(-ZARGUM)
                  
                   END DO !N2NUMB
              
                END DO !N1NUMB
          
                IF (IF_PAI.EQ.1) THEN
                    V2_PUP=V2PAIR(ENERGY,XLAMBD,DELTA2)
                    DEN_UP=DEN_UP+DNL_UP*V2_PUP
                ELSE
                    DEN_UP=DEN_UP+DNL_UP
                END IF
              
   1         CONTINUE
            
        END DO !NWSNUM
        
C      END DO !LQPRIM
C
C  Here we calculate J_DOWN block
C
      DEN_DW=0.0
C      DO LQPRIM=NSHELL,0,-1

        JNUMDW=2*LQNUMB-1          
                              
        DO NWSNUM=1,(NSHELL-LQNUMB)/2+1
            
           ENERGY=ENERDN(LQNUMB,NWSNUM)
           
           IF (JNUMDW.GE.0) THEN
                
               IF (IF_PAI.EQ.1.AND.ENERGY.GT.ENEMAX) GO TO 2
           
               IF (IF_PAI.EQ.0
     *             .AND.LABORD(NWSNUM-1,LQNUMB,JNUMDW).EQ.0) GO TO 2
            
                    DNL_DW=0.0
            
                    DO N1NUMB=0,(NSHELL-LQNUMB)/2
                
                       DO N2NUMB=0,(NSHELL-LQNUMB)/2
    
                          CFACT2=(LQNUMB)
     *                          *CMATDN(LQNUMB,N1NUMB+1,NWSNUM)
     *                          *CMATDN(LQNUMB,N2NUMB+1,NWSNUM)
               
                          POLPOL=XPOLYN(IPOINT,N1NUMB,LQNUMB)
     *                          *XPOLYN(IPOINT,N2NUMB,LQNUMB)
    
C              POLPOL=XPOLYN_DENSIT(IPOINT,N1NUMB,LQPRIM,LQNUMB) 
C     *              *XPOLYN_DENSIT(IPOINT,N2NUMB,LQPRIM,LQNUMB)
                
                          DNL_DW=DNL_DW+POLPOL*CFACT2
CID     *                          *ZARGUM**LQNUMB*EXP(-ZARGUM)
                  
                  END DO !N2NUMB
              
               END DO !N1NUMB
               
               IF (IF_PAI.EQ.1) THEN
                   V2_PDW=V2PAIR(ENERGY,XLAMBD,DELTA2)
                   DEN_DW=DEN_DW+DNL_DW*V2_PDW
               ELSE
                   DEN_DW=DEN_DW+DNL_DW
               END IF             
                                                
           END IF  !JNUMDW.GE.0
           
   2       CONTINUE            
      
        END DO !NWSNUM
        
C      END DO !LQPRIM
      
      DENSIT_LAGUER=(DEN_UP+DEN_DW)!*FACTOR
C
C=======================================================================
C     
      RETURN
      END
C      
C=======================================================================
C=======================================================================
C
      FUNCTION DENSIT_STATES(INUCLI,ISOSPI_AUXILI,ZARGUM,IPOINT,
     *                                     LQNUMB,NWSNUM,JQNUMB)
C      
      INCLUDE  'MATDIM/NDSPEC.f'
      INCLUDE  'MATDIM/NDIM_N.f'
      INCLUDE  'MATDIM/NDGAUS.f'
      INCLUDE  'MATDIM/NDNUCL.f'
C      
      PARAMETER
     *         (NDIM_L=NDIM_N,NDBASE=NDIM_N+1) 
      DIMENSION
     *          LABORD(0:NDIM_N,0:NDIM_L,0:2*NDIM_L+1)
      DIMENSION
     *          CMATUP(0:NDIM_L,1:NDBASE,1:NDBASE),
     *          CMATDN(0:NDIM_L,1:NDBASE,1:NDBASE)
      DIMENSION
     *          ENERUP(0:NDIM_L,1:NDBASE),
     *          ENERDN(0:NDIM_L,1:NDBASE)
      DIMENSION
     *          VCOEUP(0:NDIM_L,1:NDBASE),
     *          VCOEDN(0:NDIM_L,1:NDBASE)
      DIMENSION
     *          NWSSPH(1:NDSPEC),
     *          LWSSPH(1:NDSPEC),
     *          JWSSPH(1:NDSPEC)
      DIMENSION
     *          DENSIT_VECTOR(1:NDGAUS,0:NDIM_L),
     *          DENAUX_TEST_L(0:NDIM_L)
C
      COMMON
     *       /DATA_V/ AOSCIL_PROTON(1:NDNUCL),
     *                AOSCIL_NEUTRS(1:NDNUCL),
     *                HOMEGA(1:NDNUCL)   
      COMMON
     *       /GAUSHE/ NGAUSS,NSHELL_PROTON,NSHELL_NEUTRS
      COMMON
     *       /POLPOL/ XPOLYN(1:NDGAUS,0:NDIM_N,0:NDIM_L), 
     *                XDPOLY(1:NDGAUS,0:NDIM_N,0:NDIM_L)
      COMMON
     *       /POLDEN/ 
     *        XPOLYN_DENSIT(1:NDGAUS,0:NDIM_N,0:NDIM_L,0:NDIM_L), 
     *        XDPOLY_DENSIT(1:NDGAUS,0:NDIM_N,0:NDIM_L,0:NDIM_L)
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
     *       /V2_COE/ V2_PUP_VECTOR(0:NDIM_L,1:NDBASE),
     *                V2_PDW_VECTOR(0:NDIM_L,1:NDBASE),
     *                DUP_LN(0:NDIM_L,1:NDBASE),
     *                DDW_LN(0:NDIM_L,1:NDBASE),
     *                ENERGY,XLAMBD,DELTA2
      COMMON
     *       /V2COEF/ VCOEUP_PROTON(0:NDIM_L,1:NDBASE),
     *                VCOEDN_PROTON(0:NDIM_L,1:NDBASE),
     *
     *                VCOEUP_NEUTRS(0:NDIM_L,1:NDBASE),
     *                VCOEDN_NEUTRS(0:NDIM_L,1:NDBASE)
      COMMON
     *       /DIMSPH/ LDSING_PROTON,LEVTHE_PROTON,
     *                LEVEXP_PROTON(1:NDNUCL),
     *                LDSING_NEUTRS,LEVTHE_NEUTRS,
     *                LEVEXP_NEUTRS(1:NDNUCL)
      COMMON
     *       /QUANTN/ NWSSPH_PROTON(1:NDSPEC),
     *      	      LWSSPH_PROTON(1:NDSPEC),
     *      	      JWSSPH_PROTON(1:NDSPEC),
     *
     *      	      NWSSPH_NEUTRS(1:NDSPEC),
     *      	      LWSSPH_NEUTRS(1:NDSPEC),
     *      	      JWSSPH_NEUTRS(1:NDSPEC)
       COMMON
     *       /LEVORD/ LABORD_PROTON(1:NDNUCL,
     *                              0:NDIM_N,0:NDIM_L,0:2*NDIM_L+1),
     *                LABORD_NEUTRS(1:NDNUCL,
     *                              0:NDIM_N,0:NDIM_L,0:2*NDIM_L+1)
      COMMON
     *       /NUCLEU/ IZ_FIX,IN_FIX,ISOSPI 
      COMMON
     *       /IFPAIR/ IF_PAI,IFGPAI,IFDPAI
      COMMON
     *       /BCSENE/ ENEMAX_BCSPAI(1:NDNUCL,1:2)
      COMMON
     *       /BCSDEL/ DELT_2(1:NDNUCL,1:2)
      COMMON
     *       /BCSLAM/ XLAMBD_PROTON(1:NDNUCL),
     *                XLAMBD_NEUTRS(1:NDNUCL)
C
C=======================================================================   
C  
C       This FUNCTION calculates the PARTIAL DENSITY FUNCTION
C        Eq.(4.7.14) at the LAGUERRE NODES.
C          
C          It is prepared to distinct between neutrons and protons 
C
C=======================================================================     
C
C                        P R O T O N S 
C
C=======================================================================    
C     
      IF (ISOSPI_AUXILI.EQ.1) THEN
C        
          N_PART=IZ_FIX
          NSHELL=NSHELL_PROTON
          LDSING=LDSING_PROTON
C        
          AOSCIL=AOSCIL_PROTON(INUCLI)
          
          IF (IF_PAI.EQ.1) THEN
              XLAMBD=XLAMBD_PROTON(INUCLI)
              DELTA2=DELT_2(INUCLI,1)
              ENEMAX=ENEMAX_BCSPAI(INUCLI,1)
          END IF
C        
          DO INDEXX=1,LDSING
C          
             LWSSPH(INDEXX)=LWSSPH_PROTON(INDEXX)
             NWSSPH(INDEXX)=NWSSPH_PROTON(INDEXX)
             JWSSPH(INDEXX)=JWSSPH_PROTON(INDEXX)
C          
          END DO
C        
          DO INDEXX=1,LDSING
C
             JQNAUX=JWSSPH(INDEXX)
             LQNAUX=LWSSPH(INDEXX)
	     NNBAUX=NWSSPH(INDEXX)
             NBBASE=(NSHELL-LQNAUX)/2  
C             
             LABORD(NNBAUX,LQNAUX,JQNAUX)
     *	    =LABORD_PROTON(INUCLI,NNBAUX,LQNAUX,JQNAUX)
C     
             ENERUP(LQNAUX,NNBAUX+1)
     *      =ENERUP_AUXPRO(LQNAUX,NNBAUX+1)          
C
             ENERDN(LQNAUX,NNBAUX+1)
     *      =ENERDN_AUXPRO(LQNAUX,NNBAUX+1)
C
             DO N1NUMB=0,NBBASE
C		  
                CMATUP(LQNAUX,N1NUMB+1,NNBAUX+1)
     *         =CMATUP_AUXPRO(LQNAUX,N1NUMB+1,NNBAUX+1)
C		  
                CMATDN(LQNAUX,N1NUMB+1,NNBAUX+1)
     *         =CMATDN_AUXPRO(LQNAUX,N1NUMB+1,NNBAUX+1)
C
             END DO
          END DO                  
C
      END IF
C
C=======================================================================     
C
C                        N E U T R O N S 
C
C=======================================================================          
C      
      IF (ISOSPI_AUXILI.EQ.0) THEN
        
          N_PART=IN_FIX
          NSHELL=NSHELL_NEUTRS
          LDSING=LDSING_NEUTRS
        
          AOSCIL=AOSCIL_NEUTRS(INUCLI)
          
          IF (IF_PAI.EQ.1) THEN
              XLAMBD=XLAMBD_NEUTRS(INUCLI)
              DELTA2=DELT_2(INUCLI,2)
              ENEMAX=ENEMAX_BCSPAI(INUCLI,2)
          END IF
        
          DO INDEXX=1,LDSING
          
             LWSSPH(INDEXX)=LWSSPH_NEUTRS(INDEXX)
             NWSSPH(INDEXX)=NWSSPH_NEUTRS(INDEXX)
             JWSSPH(INDEXX)=JWSSPH_NEUTRS(INDEXX)
          
          END DO
        
          DO INDEXX=1,LDSING
C
             JQNAUX=JWSSPH(INDEXX)
             LQNAUX=LWSSPH(INDEXX)
	     NNBAUX=NWSSPH(INDEXX)
             NBBASE=(NSHELL-LQNAUX)/2  
             
             LABORD(NNBAUX,LQNAUX,JQNAUX)
     *	    =LABORD_NEUTRS(INUCLI,NNBAUX,LQNAUX,JQNAUX)
     
             ENERUP(LQNAUX,NNBAUX+1)
     *      =ENERUP_AUXNEU(LQNAUX,NNBAUX+1)          
C
             ENERDN(LQNAUX,NNBAUX+1)
     *      =ENERDN_AUXNEU(LQNAUX,NNBAUX+1)
C
             DO N1NUMB=0,NBBASE
C		  
                CMATUP(LQNAUX,N1NUMB+1,NNBAUX+1)
     *         =CMATUP_AUXNEU(LQNAUX,N1NUMB+1,NNBAUX+1)
C		  
                CMATDN(LQNAUX,N1NUMB+1,NNBAUX+1)
     *         =CMATDN_AUXNEU(LQNAUX,N1NUMB+1,NNBAUX+1)
C
             END DO
          END DO                  
C
      END IF
C
C=======================================================================     
C
C                   Calculating the DENSITY Eq.(4.7.14)
C
C=======================================================================        
C      
      PINUMB=4.0*ATAN(1.0)
      FACTOR=1.0/(PINUMB*(AOSCIL**3)) ! multiplied by 2/a 
C                                             from N_{n1,l}*N_{n2,l}

      DENAUX=0.0
      
      V2_PUP_VECTOR(LQNUMB,NWSNUM)=0.0
      V2_PDW_VECTOR(LQNUMB,NWSNUM)=0.0

C
C  Here we calculate J_UP block
C      
      DEN_UP=0.0
      
      V2_SUM=0.0
C
C      DO LQPRIM=NSHELL,0,-1
C          
        JNUMUP=2*LQNUMB+1
        
C         DO NWSNUM=1,(NSHELL-LQNUMB)/2+1
        
        IF (JQNUMB.EQ.JNUMUP) THEN
            
            ENERGY=ENERUP(LQNUMB,NWSNUM)
           
            IF (IF_PAI.EQ.1.AND.ENERGY.GT.ENEMAX) GO TO 1
           
            IF (IF_PAI.EQ.0
     *          .AND.LABORD(NWSNUM-1,LQNUMB,JNUMUP).EQ.0) GO TO 1
     
                 DNL_UP=0.0
            
                 DO N1NUMB=0,(NSHELL-LQNUMB)/2
                
                    DO N2NUMB=0,(NSHELL-LQNUMB)/2
                  
                       CFACT1=(LQNUMB+1)*CMATUP(LQNUMB,N1NUMB+1,NWSNUM)
     *                                  *CMATUP(LQNUMB,N2NUMB+1,NWSNUM)
               
                      POLPOL=XPOLYN(IPOINT,N1NUMB,LQNUMB)
     *                      *XPOLYN(IPOINT,N2NUMB,LQNUMB)
    
C              POLPOL=XPOLYN_DENSIT(IPOINT,N1NUMB,LQPRIM,LQNUMB) 
C     *              *XPOLYN_DENSIT(IPOINT,N2NUMB,LQPRIM,LQNUMB)
                
                      DNL_UP=DNL_UP+POLPOL*CFACT1
CID     *                             *ZARGUM**LQNUMB*EXP(-ZARGUM)
                  
                   END DO !N2NUMB
              
                END DO !N1NUMB
          
                IF (IF_PAI.EQ.1) THEN
                    
                    V2_PUP=V2PAIR(ENERGY,XLAMBD,DELTA2)
                    V2_PUP_VECTOR(LQNUMB,NWSNUM)=V2_PUP
                    
                    DEN_UP=DEN_UP+DNL_UP*V2_PUP
                    DUP_LN(LQNUMB,NWSNUM)=DNL_UP
                    
                    V2_SUM=V2_SUM+V2_PUP
                    
                ELSE
                    DEN_UP=DEN_UP+DNL_UP
                END IF
              
   1         CONTINUE
   
        DENAUX=DEN_UP
   
        END IF 
            
C        END DO !NWSNUM
        
C      END DO !LQPRIM
C      
C
C  Here we calculate J_DOWN block
C
      DEN_DW=0.0
C      DO LQPRIM=NSHELL,0,-1

        JNUMDW=2*LQNUMB-1          
                              
C        DO NWSNUM=1,(NSHELL-LQNUMB)/2+1
            
           ENERGY=ENERDN(LQNUMB,NWSNUM)
           
           IF (JNUMDW.EQ.JQNUMB.AND.JNUMDW.GE.0) THEN
                
               IF (IF_PAI.EQ.1.AND.ENERGY.GT.ENEMAX) GO TO 2
           
               IF (IF_PAI.EQ.0
     *             .AND.LABORD(NWSNUM-1,LQNUMB,JNUMDW).EQ.0) GO TO 2
            
                    DNL_DW=0.0
            
                    DO N1NUMB=0,(NSHELL-LQNUMB)/2
                
                       DO N2NUMB=0,(NSHELL-LQNUMB)/2
    
                          CFACT2=(LQNUMB)
     *                          *CMATDN(LQNUMB,N1NUMB+1,NWSNUM)
     *                          *CMATDN(LQNUMB,N2NUMB+1,NWSNUM)
               
                          POLPOL=XPOLYN(IPOINT,N1NUMB,LQNUMB)
     *                          *XPOLYN(IPOINT,N2NUMB,LQNUMB)
    
C              POLPOL=XPOLYN_DENSIT(IPOINT,N1NUMB,LQPRIM,LQNUMB) 
C     *              *XPOLYN_DENSIT(IPOINT,N2NUMB,LQPRIM,LQNUMB)
                
                          DNL_DW=DNL_DW+POLPOL*CFACT2
CID     *                          *ZARGUM**LQNUMB*EXP(-ZARGUM)
                  
                  END DO !N2NUMB
              
               END DO !N1NUMB
               
               IF (IF_PAI.EQ.1) THEN
                   
                   V2_PDW=V2PAIR(ENERGY,XLAMBD,DELTA2)
                   V2_PDW_VECTOR(LQNUMB,NWSNUM)=V2_PDW
                   
                   DEN_DW=DEN_DW+DNL_DW*V2_PDW
                   DDW_LN(LQNUMB,NWSNUM)=DNL_DW
                   
                   V2_SUM=V2_SUM+V2_PDW
                   
               ELSE
                   DEN_DW=DEN_DW+DNL_DW
               END IF   
               
               DENAUX=DEN_DW
                                                
           END IF  !JNUMDW.GE.0
           
   2       CONTINUE            
      
C        END DO !NWSNUM
        
C      END DO !LQPRIM
      
      DENSIT_STATES=DENAUX!*FACTOR
C
C=======================================================================
C     
      RETURN
      END
C      
C=======================================================================
C=======================================================================
C
      FUNCTION DENSIT_GRADNT(INUCLI,ISOSPI_AUXILI,ZARGUM,IPOINT,LQNUMB)
C          
      INCLUDE  'MATDIM/NDSPEC.f'
      INCLUDE  'MATDIM/NDIM_N.f'
      INCLUDE  'MATDIM/NDGAUS.f'
      INCLUDE  'MATDIM/NDNUCL.f'
C      
      PARAMETER
     *         (NDIM_L=NDIM_N,NDBASE=NDIM_N+1) 
C      
      DIMENSION
     *          LABORD(0:NDIM_N,0:NDIM_L,0:2*NDIM_L+1)
      DIMENSION
     *          CMATUP(0:NDIM_L,1:NDBASE,1:NDBASE),
     *          CMATDN(0:NDIM_L,1:NDBASE,1:NDBASE)
      DIMENSION
     *          ENERUP(0:NDIM_L,1:NDBASE),
     *          ENERDN(0:NDIM_L,1:NDBASE)
      DIMENSION
     *          VCOEUP(0:NDIM_L,1:NDBASE),
     *          VCOEDN(0:NDIM_L,1:NDBASE)
      DIMENSION
     *          NWSSPH(1:NDSPEC),
     *          LWSSPH(1:NDSPEC),
     *          JWSSPH(1:NDSPEC)
      DIMENSION
     *          DENDER(1:NDGAUS,0:NDIM_N,0:NDIM_N,0:NDIM_L),
     *          AUXILI(1:NDGAUS)
C
      COMMON
     *       /POLDEN/ 
     *        XPOLYN_DENSIT(1:NDGAUS,0:NDIM_N,0:NDIM_L,0:NDIM_L), 
     *        XDPOLY_DENSIT(1:NDGAUS,0:NDIM_N,0:NDIM_L,0:NDIM_L)
      COMMON
     *       /DATA_V/ AOSCIL_PROTON(1:NDNUCL),
     *                AOSCIL_NEUTRS(1:NDNUCL),
     *                HOMEGA(1:NDNUCL)    
      COMMON
     *       /GAUSHE/ NGAUSS,NSHELL_PROTON,NSHELL_NEUTRS
      COMMON
     *       /NMLCNT/ XNORNL(0:NDIM_N,0:NDIM_L)
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
     *       /DIMSPH/ LDSING_PROTON,LEVTHE_PROTON,
     *                LEVEXP_PROTON(1:NDNUCL),
     *                LDSING_NEUTRS,LEVTHE_NEUTRS,
     *                LEVEXP_NEUTRS(1:NDNUCL)
      COMMON
     *       /QUANTN/ NWSSPH_PROTON(1:NDSPEC),
     *      	      LWSSPH_PROTON(1:NDSPEC),
     *      	      JWSSPH_PROTON(1:NDSPEC),
     *
     *      	      NWSSPH_NEUTRS(1:NDSPEC),
     *      	      LWSSPH_NEUTRS(1:NDSPEC),
     *      	      JWSSPH_NEUTRS(1:NDSPEC)
       COMMON
     *       /LEVORD/ LABORD_PROTON(1:NDNUCL,
     *                              0:NDIM_N,0:NDIM_L,0:2*NDIM_L+1),
     *                LABORD_NEUTRS(1:NDNUCL,
     *                              0:NDIM_N,0:NDIM_L,0:2*NDIM_L+1)
      COMMON
     *       /DENSIF/ ITRNUM
      COMMON
     *       /NUCLEU/ IZ_FIX,IN_FIX,ISOSPI      
      COMMON
     *       /IFPAIR/ IF_PAI,IFGPAI,IFDPAI
      COMMON
     *       /BCSENE/ ENEMAX_BCSPAI(1:NDNUCL,1:2)
      COMMON
     *       /BCSDEL/ DELT_2(1:NDNUCL,1:2)
      COMMON
     *       /BCSLAM/ XLAMBD_PROTON(1:NDNUCL),
     *                XLAMBD_NEUTRS(1:NDNUCL)
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
       CALL CPUTIM('DENGRD',1)
C       
       ICOUNT_DENGRD=ICOUNT_DENGRD+1
C
C=======================================================================   
C  
CID        This FUNCTION calculates the GRADIENT OF THE DENSITY 
C               Eq.(4.7.21)=(1/r)*(drho/dr) from MFV12.pdf
C
C=======================================================================     
C
c      IF (LOGWRI.GT.0) THEN
c          WRITE(LOGFIL,'(/,9X,''Entering DENSIT_GRADNT'',/,
c     *                   9X,/,9X,''ISOSPI_AUXILI='',I1)')ISOSPI_AUXILI
c      END IF
C
C=======================================================================     
C
C                        P R O T O N S 
C
C=======================================================================  
C    
      IF (ISOSPI_AUXILI.EQ.1) THEN
C        
          N_PART=IZ_FIX
          NSHELL=NSHELL_PROTON
          LDSING=LDSING_PROTON
C        
          AOSCIL=AOSCIL_PROTON(INUCLI)
          
          IF (IF_PAI.EQ.1) THEN
              XLAMBD=XLAMBD_PROTON(INUCLI)
              DELTA2=DELT_2(INUCLI,1)
              ENEMAX=ENEMAX_BCSPAI(INUCLI,1)
          END IF
C        
          DO INDEXX=1,LDSING
C          
             LWSSPH(INDEXX)=LWSSPH_PROTON(INDEXX)
             NWSSPH(INDEXX)=NWSSPH_PROTON(INDEXX)
             JWSSPH(INDEXX)=JWSSPH_PROTON(INDEXX)
C          
          END DO
          
          DO INDEXX=1,LDSING
C
             JQNAUX=JWSSPH(INDEXX)
             LQNAUX=LWSSPH(INDEXX)
	         NNBAUX=NWSSPH(INDEXX)
             NBBASE=(NSHELL-LQNAUX)/2  
C             
             LABORD(NNBAUX,LQNAUX,JQNAUX)
     *	    =LABORD_PROTON(INUCLI,NNBAUX,LQNAUX,JQNAUX)
C     
             ENERUP(LQNAUX,NNBAUX+1)
     *      =ENERUP_AUXPRO(LQNAUX,NNBAUX+1)          

             ENERDN(LQNAUX,NNBAUX+1)
     *      =ENERDN_AUXPRO(LQNAUX,NNBAUX+1)
C
             DO N1NUMB=0,NBBASE
C		  
                CMATUP(LQNAUX,N1NUMB+1,NNBAUX+1)
     *         =CMATUP_AUXPRO(LQNAUX,N1NUMB+1,NNBAUX+1)
C		  
                CMATDN(LQNAUX,N1NUMB+1,NNBAUX+1)
     *         =CMATDN_AUXPRO(LQNAUX,N1NUMB+1,NNBAUX+1)
C
             END DO
         END DO                  
C
      END IF
C
C=======================================================================     
C
C                        N E U T R O N S 
C
C=======================================================================          
C      
      IF (ISOSPI_AUXILI.EQ.0) THEN
        
        N_PART=IN_FIX
        NSHELL=NSHELL_NEUTRS
        LDSING=LDSING_NEUTRS
        
        AOSCIL=AOSCIL_NEUTRS(INUCLI)
        
        IF (IF_PAI.EQ.1) THEN
            XLAMBD=XLAMBD_NEUTRS(INUCLI)
            DELTA2=DELT_2(INUCLI,2)
            ENEMAX=ENEMAX_BCSPAI(INUCLI,2)
        END IF
        
        DO INDEXX=1,LDSING
          
          LWSSPH(INDEXX)=LWSSPH_NEUTRS(INDEXX)
          NWSSPH(INDEXX)=NWSSPH_NEUTRS(INDEXX)
          JWSSPH(INDEXX)=JWSSPH_NEUTRS(INDEXX)
          
        END DO
C
        DO INDEXX=1,LDSING
C
             JQNAUX=JWSSPH(INDEXX)
             LQNAUX=LWSSPH(INDEXX)
             NNBAUX=NWSSPH(INDEXX)
             NBBASE=(NSHELL-LQNAUX)/2  
             
             LABORD(NNBAUX,LQNAUX,JQNAUX)
     *	    =LABORD_NEUTRS(INUCLI,NNBAUX,LQNAUX,JQNAUX)
C
             ENERUP(LQNAUX,NNBAUX+1)
     *      =ENERUP_AUXNEU(LQNAUX,NNBAUX+1)          
C
             ENERDN(LQNAUX,NNBAUX+1)
     *      =ENERDN_AUXNEU(LQNAUX,NNBAUX+1)
C     
             DO N1NUMB=0,NBBASE
C		  
                CMATUP(LQNAUX,N1NUMB+1,NNBAUX+1)
     *         =CMATUP_AUXNEU(LQNAUX,N1NUMB+1,NNBAUX+1)
C		  
                CMATDN(LQNAUX,N1NUMB+1,NNBAUX+1)
     *         =CMATDN_AUXNEU(LQNAUX,N1NUMB+1,NNBAUX+1)
     
C
             END DO
        END DO                  
C
      END IF
C
C=======================================================================     
C
C                 Calculating the DENSITY GRADIENT Eq.(4.7.21)
C
C=======================================================================        
C   
      PINUMB=4.*ATAN(1.0)
      FACTOR=2./(PINUMB*AOSCIL**5)
C
C  Here we calculate J_UP block
C      
      DEN_UP=0.0
      TESTV2=0.0
      
      DO LQPRIM=NSHELL,0,-1 
         
         JNUMUP=2*LQPRIM+1  
          
         DO NWSNUM=1,(NSHELL-LQPRIM)/2+1
            
            ENERGY=ENERUP(LQPRIM,NWSNUM)
            
            IF (IF_PAI.EQ.1.AND.ENERGY.GT.ENEMAX) GO TO 1
           
            IF (IF_PAI.EQ.0
     *          .AND.LABORD(NWSNUM-1,LQPRIM,JNUMUP).EQ.0) GO TO 1
                
                DNL_UP=0.0
                
                DO N1NUMB=0,(NSHELL-LQPRIM)/2
                
                   DO N2NUMB=0,(NSHELL-LQPRIM)/2
                                         
                      CFACT1=(LQPRIM+1)*CMATUP(LQPRIM,N1NUMB+1,NWSNUM)
     *                                 *CMATUP(LQPRIM,N2NUMB+1,NWSNUM)
              
                      POLPOL=XPOLYN_DENSIT(IPOINT,N1NUMB,LQPRIM,LQNUMB) 
     *                      *XPOLYN_DENSIT(IPOINT,N2NUMB,LQPRIM,LQNUMB)
     
                      POLDER=XPOLYN_DENSIT(IPOINT,N1NUMB,LQPRIM,LQNUMB) 
     *                      *XDPOLY_DENSIT(IPOINT,N2NUMB,LQPRIM,LQNUMB)
     
                      DERPOL=XDPOLY_DENSIT(IPOINT,N1NUMB,LQPRIM,LQNUMB) 
     *                      *XPOLYN_DENSIT(IPOINT,N2NUMB,LQPRIM,LQNUMB)
              
                      ZLTLAG=(LQPRIM-ZARGUM)*POLPOL
     *                      +ZARGUM*(POLDER+DERPOL)
                
                      DNL_UP=DNL_UP+EXP(-ZARGUM)*ZARGUM**(LQPRIM-1)
     *                             *CFACT1*ZLTLAG
     
                  END DO !N2NUMB
              
               END DO !N1NUMB
               
               IF (IF_PAI.EQ.1) THEN
                   
                   V2_PUP=V2PAIR(ENERGY,XLAMBD,DELTA2)
                   DEN_UP=DEN_UP+DNL_UP*V2_PUP
                   
                   TESTV2=TESTV2+V2_PUP*(JNUMUP+1)   
                   
               ELSE
                   DEN_UP=DEN_UP+DNL_UP
               END IF
              
   1         CONTINUE
           
         END DO !NWSNUM
      
      END DO !LQPRIM
C
C  Here we calculate J_DOWN block
C
      DEN_DW=0.0
      
      DO LQPRIM=NSHELL,0,-1 
          
         JNUMDW=2*LQPRIM-1
          
         DO NWSNUM=1,(NSHELL-LQPRIM)/2+1
             
            ENERGY=ENERDN(LQPRIM,NWSNUM)
           
            IF (JNUMDW.GE.0) THEN
                
                IF (IF_PAI.EQ.1.AND.ENERGY.GT.ENEMAX) GO TO 2
           
                IF (IF_PAI.EQ.0
     *              .AND.LABORD(NWSNUM-1,LQPRIM,JNUMDW).EQ.0) GO TO 2
                    
                    DNL_DW=0.0
                    
                    DO N1NUMB=0,(NSHELL-LQPRIM)/2
                
                       DO N2NUMB=0,(NSHELL-LQPRIM)/2
    
                          CFACT2=(LQPRIM)
     *                          *CMATDN(LQPRIM,N1NUMB+1,NWSNUM)
     *                          *CMATDN(LQPRIM,N2NUMB+1,NWSNUM)
              
                          POLPOL=XPOLYN_DENSIT(IPOINT,N1NUMB,LQPRIM,
     *                                                       LQNUMB) 
     *                          *XPOLYN_DENSIT(IPOINT,N2NUMB,LQPRIM,
     *                                                       LQNUMB) 
     
                          POLDER=XPOLYN_DENSIT(IPOINT,N1NUMB,LQPRIM,
     *                                                       LQNUMB)  
     *                          *XDPOLY_DENSIT(IPOINT,N2NUMB,LQPRIM,
     *                                                       LQNUMB) 
     
                          DERPOL=XDPOLY_DENSIT(IPOINT,N1NUMB,LQPRIM,
     *                                                       LQNUMB) 
     *                          *XPOLYN_DENSIT(IPOINT,N2NUMB,LQPRIM,
     *                                                       LQNUMB) 
              
                          ZLTLAG=(LQPRIM-ZARGUM)*POLPOL
     *                          +ZARGUM*(POLDER+DERPOL)
     
                          DNL_DW=DNL_DW+EXP(-ZARGUM)*ZARGUM**(LQPRIM-1)
     *                                 *CFACT2*ZLTLAG
                
                       END DO !N2NUMB
              
                    END DO !N1NUMB
                    
                    IF (IF_PAI.EQ.1) THEN
                        
                        V2_PDW=V2PAIR(ENERGY,XLAMBD,DELTA2)
                        DEN_DW=DEN_DW+DNL_DW*V2_PDW
                        
                        TESTV2=TESTV2+V2_PDW*(JNUMDW+1) 
                        
                   ELSE
                        DEN_DW=DEN_DW+DNL_DW
                   END IF

            END IF
            
   2        CONTINUE 
                
         END DO !NWSNUM
      
      END DO !LQPRIM
C
C=======================================================================     
C      
      DENSIT_GRADNT=(DEN_UP+DEN_DW)*FACTOR
C
C=======================================================================     
C
      IF (IF_PAI.EQ.1) THEN
          IF (ABS(TESTV2-N_PART).GT.1.E-4) THEN
              WRITE(LOGFIL,'(''ITRNUM= '',I2)')ITRNUM
              WRITE(LOGFIL,'(''ISOSPI_AUXILI='',I2)')ISOSPI_AUXILI
              WRITE(LOGFIL,'(''In DENSIT_GRADNT: TESTV2='',F20.13)')
     *                                           TESTV2
              WRITE(LOGFIL,'(18X,''N_PART='',I6)')N_PART
     
              IF (ISCREN.GE.5) THEN
                  WRITE(LSCREN,'(''ITRNUM= '',I2)')ITRNUM
                  WRITE(LSCREN,'(''ISOSPI_AUXILI='',I2)')ISOSPI_AUXILI
                  WRITE(LSCREN,'(''In DENSIT_GRADNT: TESTV2='',F20.13)')
     *                                          TESTV2
                  WRITE(LSCREN,'(18X,''N_PART='',I6)')N_PART
              END IF
     
C              STOP 'Stop in DENSIT_GRADNT: TESTV2 not correct!'
          END IF
      END IF
C
C=======================================================================     
C
c      IF (LOGWRI.GT.0) THEN
c          WRITE(LOGFIL,'(/,9X,''Exiting DENSIT_GRADNT'')')
c      END IF
C
C=======================================================================
C
      CALL CPUTIM('DENGRD',0)
C
C=======================================================================
C      
      RETURN
      END

C      
C=======================================================================
C=======================================================================
C
      FUNCTION DENGRD_PLOTIN(INUCLI,ISOSPI_AUXILI,ZARGUM)
C          
      INCLUDE  'MATDIM/NDSPEC.f'
      INCLUDE  'MATDIM/NDIM_N.f'
      INCLUDE  'MATDIM/NDGAUS.f'
      INCLUDE  'MATDIM/NDNUCL.f'
C      
      PARAMETER
     *         (NDIM_L=NDIM_N,NDBASE=NDIM_N+1) 
C      
      REAL*16
     *          QARGUM,QFLAGN,QDFLAG,QNORNL
      DIMENSION
     *          LABORD(0:NDIM_N,0:NDIM_L,0:2*NDIM_L+1)
      DIMENSION
     *          CMATUP(0:NDIM_L,1:NDBASE,1:NDBASE),
     *          CMATDN(0:NDIM_L,1:NDBASE,1:NDBASE)
      DIMENSION
     *          ENERUP(0:NDIM_L,1:NDBASE),
     *          ENERDN(0:NDIM_L,1:NDBASE)
      DIMENSION
     *          VCOEUP(0:NDIM_L,1:NDBASE),
     *          VCOEDN(0:NDIM_L,1:NDBASE)
      DIMENSION
     *          NWSSPH(1:NDSPEC),
     *          LWSSPH(1:NDSPEC),
     *          JWSSPH(1:NDSPEC)
      DIMENSION
     *          DENDER(1:NDGAUS,0:NDIM_N,0:NDIM_N,0:NDIM_L),
     *          AUXILI(1:NDGAUS)
      DIMENSION
     *          QFLAGN(0:NDIM_N,0:NDIM_L),
     *          QDFLAG(0:NDIM_N,0:NDIM_L),
     *          XFLAGN(0:NDIM_N,0:NDIM_L),
     *          XDFLAG(0:NDIM_N,0:NDIM_L)
C
      COMMON
     *       /POLDEN/ 
     *        XPOLYN_DENSIT(1:NDGAUS,0:NDIM_N,0:NDIM_L,0:NDIM_L), 
     *        XDPOLY_DENSIT(1:NDGAUS,0:NDIM_N,0:NDIM_L,0:NDIM_L)
      COMMON
     *       /DATA_V/ AOSCIL_PROTON(1:NDNUCL),
     *                AOSCIL_NEUTRS(1:NDNUCL),
     *                HOMEGA(1:NDNUCL)    
      COMMON
     *       /GAUSHE/ NGAUSS,NSHELL_PROTON,NSHELL_NEUTRS
      COMMON
     *       /NMLCNT/ XNORNL(0:NDIM_N,0:NDIM_L)
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
     *       /DIMSPH/ LDSING_PROTON,LEVTHE_PROTON,
     *                LEVEXP_PROTON(1:NDNUCL),
     *                LDSING_NEUTRS,LEVTHE_NEUTRS,
     *                LEVEXP_NEUTRS(1:NDNUCL)
      COMMON
     *       /QUANTN/ NWSSPH_PROTON(1:NDSPEC),
     *      	      LWSSPH_PROTON(1:NDSPEC),
     *      	      JWSSPH_PROTON(1:NDSPEC),
     *
     *      	      NWSSPH_NEUTRS(1:NDSPEC),
     *      	      LWSSPH_NEUTRS(1:NDSPEC),
     *      	      JWSSPH_NEUTRS(1:NDSPEC)
       COMMON
     *       /LEVORD/ LABORD_PROTON(1:NDNUCL,
     *                              0:NDIM_N,0:NDIM_L,0:2*NDIM_L+1),
     *                LABORD_NEUTRS(1:NDNUCL,
     *                              0:NDIM_N,0:NDIM_L,0:2*NDIM_L+1)
      COMMON
     *       /DENSIF/ ITRNUM
      COMMON
     *       /NUCLEU/ IZ_FIX,IN_FIX,ISOSPI      
      COMMON
     *       /IFPAIR/ IF_PAI,IFGPAI,IFDPAI
      COMMON
     *       /BCSENE/ ENEMAX_BCSPAI(1:NDNUCL,1:2)
      COMMON
     *       /BCSDEL/ DELT_2(1:NDNUCL,1:2)
      COMMON
     *       /BCSLAM/ XLAMBD_PROTON(1:NDNUCL),
     *                XLAMBD_NEUTRS(1:NDNUCL)
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
       CALL CPUTIM('DENGRD',1)
C       
       ICOUNT_DENGRD=ICOUNT_DENGRD+1
C
C=======================================================================   
C  
CID        This FUNCTION calculates the GRADIENT OF THE DENSITY 
C               Eq.(4.7.21)=(1/r)*(drho/dr) from MFV12.pdf
C
C=======================================================================     
C
c      IF (LOGWRI.GT.0) THEN
c          WRITE(LOGFIL,'(/,9X,''Entering DENSIT_GRADNT'',/,
c     *                   9X,/,9X,''ISOSPI_AUXILI='',I1)')ISOSPI_AUXILI
c      END IF
C
C=======================================================================     
C
C                        P R O T O N S 
C
C=======================================================================  
C    
      IF (ISOSPI_AUXILI.EQ.1) THEN
C        
          N_PART=IZ_FIX
          NSHELL=NSHELL_PROTON
          LDSING=LDSING_PROTON
C        
          AOSCIL=AOSCIL_PROTON(INUCLI)
          
          IF (IF_PAI.EQ.1) THEN
              XLAMBD=XLAMBD_PROTON(INUCLI)
              DELTA2=DELT_2(INUCLI,1)
              ENEMAX=ENEMAX_BCSPAI(INUCLI,1)
          END IF
C        
          DO INDEXX=1,LDSING
C          
             LWSSPH(INDEXX)=LWSSPH_PROTON(INDEXX)
             NWSSPH(INDEXX)=NWSSPH_PROTON(INDEXX)
             JWSSPH(INDEXX)=JWSSPH_PROTON(INDEXX)
C          
          END DO
          
          DO INDEXX=1,LDSING
C
             JQNAUX=JWSSPH(INDEXX)
             LQNAUX=LWSSPH(INDEXX)
	     NNBAUX=NWSSPH(INDEXX)
             NBBASE=(NSHELL-LQNAUX)/2  
C             
             LABORD(NNBAUX,LQNAUX,JQNAUX)
     *	    =LABORD_PROTON(INUCLI,NNBAUX,LQNAUX,JQNAUX)
C     
             ENERUP(LQNAUX,NNBAUX+1)
     *      =ENERUP_AUXPRO(LQNAUX,NNBAUX+1)          

             ENERDN(LQNAUX,NNBAUX+1)
     *      =ENERDN_AUXPRO(LQNAUX,NNBAUX+1)
C
             DO N1NUMB=0,NBBASE
C		  
                CMATUP(LQNAUX,N1NUMB+1,NNBAUX+1)
     *         =CMATUP_AUXPRO(LQNAUX,N1NUMB+1,NNBAUX+1)
C		  
                CMATDN(LQNAUX,N1NUMB+1,NNBAUX+1)
     *         =CMATDN_AUXPRO(LQNAUX,N1NUMB+1,NNBAUX+1)
C
             END DO
         END DO                  
C
      END IF
C
C=======================================================================     
C
C                        N E U T R O N S 
C
C=======================================================================          
C      
      IF (ISOSPI_AUXILI.EQ.0) THEN
        
        N_PART=IN_FIX
        NSHELL=NSHELL_NEUTRS
        LDSING=LDSING_NEUTRS
        
        AOSCIL=AOSCIL_NEUTRS(INUCLI)
        
        IF (IF_PAI.EQ.1) THEN
            XLAMBD=XLAMBD_NEUTRS(INUCLI)
            DELTA2=DELT_2(INUCLI,2)
            ENEMAX=ENEMAX_BCSPAI(INUCLI,2)
        END IF
        
        DO INDEXX=1,LDSING
          
          LWSSPH(INDEXX)=LWSSPH_NEUTRS(INDEXX)
          NWSSPH(INDEXX)=NWSSPH_NEUTRS(INDEXX)
          JWSSPH(INDEXX)=JWSSPH_NEUTRS(INDEXX)
          
        END DO
C
        DO INDEXX=1,LDSING
C
             JQNAUX=JWSSPH(INDEXX)
             LQNAUX=LWSSPH(INDEXX)
	     NNBAUX=NWSSPH(INDEXX)
             NBBASE=(NSHELL-LQNAUX)/2  
             
             LABORD(NNBAUX,LQNAUX,JQNAUX)
     *	    =LABORD_NEUTRS(INUCLI,NNBAUX,LQNAUX,JQNAUX)
C
             ENERUP(LQNAUX,NNBAUX+1)
     *      =ENERUP_AUXNEU(LQNAUX,NNBAUX+1)          
C
             ENERDN(LQNAUX,NNBAUX+1)
     *      =ENERDN_AUXNEU(LQNAUX,NNBAUX+1)
C     
             DO N1NUMB=0,NBBASE
C		  
                CMATUP(LQNAUX,N1NUMB+1,NNBAUX+1)
     *         =CMATUP_AUXNEU(LQNAUX,N1NUMB+1,NNBAUX+1)
C		  
                CMATDN(LQNAUX,N1NUMB+1,NNBAUX+1)
     *         =CMATDN_AUXNEU(LQNAUX,N1NUMB+1,NNBAUX+1)
     
C
             END DO
        END DO                  
C
      END IF
C
C=======================================================================     
C
C       Calculating the polynomials in the desired ZARGUM point
C
C=======================================================================          
C      
      KIND_L=1
      NORDER=NSHELL/2
      LORDER=NSHELL
      
      QARGUM=QEXT(ZARGUM)    !passing the argument to quadruple precision
      
      CALL LAGUER(QARGUM,QFLAGN,QDFLAG,NORDER,LORDER,KIND_L,
     *                                        NDIM_L,NDIM_N)
     
      DO NAUXIL=0,NORDER
        DO LAUXIL=0,LORDER
            XFLAGN(NAUXIL,LAUXIL)=REAL(QFLAGN(NAUXIL,LAUXIL))
            XDFLAG(NAUXIL,LAUXIL)=REAL(QDFLAG(NAUXIL,LAUXIL))
        END DO
      END DO
C
C=======================================================================     
C
C                 Calculating the DENSITY GRADIENT Eq.(4.7.21)
C
C=======================================================================        
C   
      PINUMB=4.*ATAN(1.0)
      FACTOR=2./(PINUMB*AOSCIL**5)
C
C  Here we calculate J_UP block
C      
      DEN_UP=0.0
      TESTV2=0.0
      
      DO LQPRIM=NSHELL,0,-1 
         
         JNUMUP=2*LQPRIM+1  
          
         DO NWSNUM=1,(NSHELL-LQPRIM)/2+1
            
            ENERGY=ENERUP(LQPRIM,NWSNUM)
            
            IF (IF_PAI.EQ.1.AND.ENERGY.GT.ENEMAX) GO TO 1
           
            IF (IF_PAI.EQ.0
     *          .AND.LABORD(NWSNUM-1,LQPRIM,JNUMUP).EQ.0) GO TO 1
                
                DNL_UP=0.0
                
                DO N1NUMB=0,(NSHELL-LQPRIM)/2
                
                   DO N2NUMB=0,(NSHELL-LQPRIM)/2
                                         
                      CFACT1=(LQPRIM+1)*CMATUP(LQPRIM,N1NUMB+1,NWSNUM)
     *                                 *CMATUP(LQPRIM,N2NUMB+1,NWSNUM)
              
                      POLPOL=XNORNL(N1NUMB,LQPRIM)*XFLAGN(N1NUMB,LQPRIM)
     *                      *XNORNL(N2NUMB,LQPRIM)*XFLAGN(N2NUMB,LQPRIM)
     
                      POLDER=XNORNL(N1NUMB,LQPRIM)*XFLAGN(N1NUMB,LQPRIM)
     *                      *XNORNL(N2NUMB,LQPRIM)*XDFLAG(N2NUMB,LQPRIM)
     
                      DERPOL=XNORNL(N1NUMB,LQPRIM)*XDFLAG(N1NUMB,LQPRIM)
     *                      *XNORNL(N2NUMB,LQPRIM)*XFLAGN(N2NUMB,LQPRIM)
              
                      ZLTLAG=(LQPRIM-ZARGUM)*POLPOL
     *                      +ZARGUM*(POLDER+DERPOL)
                
                      DNL_UP=DNL_UP+EXP(-ZARGUM)*ZARGUM**(LQPRIM-1)
     *                             *CFACT1*ZLTLAG
     
                  END DO !N2NUMB
              
               END DO !N1NUMB
               
               IF (IF_PAI.EQ.1) THEN
                   
                   V2_PUP=V2PAIR(ENERGY,XLAMBD,DELTA2)
                   DEN_UP=DEN_UP+DNL_UP*V2_PUP
                   
                   TESTV2=TESTV2+V2_PUP*(JNUMUP+1)   
                   
               ELSE
                   DEN_UP=DEN_UP+DNL_UP
               END IF
              
   1         CONTINUE
           
         END DO !NWSNUM
      
      END DO !LQPRIM
C
C  Here we calculate J_DOWN block
C
      DEN_DW=0.0
      
      DO LQPRIM=NSHELL,0,-1 
          
         JNUMDW=2*LQPRIM-1
          
         DO NWSNUM=1,(NSHELL-LQPRIM)/2+1
             
            ENERGY=ENERDN(LQPRIM,NWSNUM)
           
            IF (JNUMDW.GE.0) THEN
                
                IF (IF_PAI.EQ.1.AND.ENERGY.GT.ENEMAX) GO TO 2
           
                IF (IF_PAI.EQ.0
     *              .AND.LABORD(NWSNUM-1,LQPRIM,JNUMDW).EQ.0) GO TO 2
                    
                    DNL_DW=0.0
                    
                    DO N1NUMB=0,(NSHELL-LQPRIM)/2
                
                       DO N2NUMB=0,(NSHELL-LQPRIM)/2
    
                          CFACT2=(LQPRIM)
     *                          *CMATDN(LQPRIM,N1NUMB+1,NWSNUM)
     *                          *CMATDN(LQPRIM,N2NUMB+1,NWSNUM)
              
                          POLPOL=
     *                    XNORNL(N1NUMB,LQPRIM)*XFLAGN(N1NUMB,LQPRIM)
     *                   *XNORNL(N2NUMB,LQPRIM)*XFLAGN(N2NUMB,LQPRIM)
     
                          POLDER=
     *                    XNORNL(N1NUMB,LQPRIM)*XFLAGN(N1NUMB,LQPRIM)
     *                   *XNORNL(N2NUMB,LQPRIM)*XDFLAG(N2NUMB,LQPRIM)
     
                          DERPOL=
     *                    XNORNL(N1NUMB,LQPRIM)*XDFLAG(N1NUMB,LQPRIM)
     *                   *XNORNL(N2NUMB,LQPRIM)*XFLAGN(N2NUMB,LQPRIM)
              
                          ZLTLAG=(LQPRIM-ZARGUM)*POLPOL
     *                          +ZARGUM*(POLDER+DERPOL)
     
                          DNL_DW=DNL_DW+EXP(-ZARGUM)*ZARGUM**(LQPRIM-1)
     *                                 *CFACT2*ZLTLAG
                
                       END DO !N2NUMB
              
                    END DO !N1NUMB
                    
                    IF (IF_PAI.EQ.1) THEN
                        
                        V2_PDW=V2PAIR(ENERGY,XLAMBD,DELTA2)
                        DEN_DW=DEN_DW+DNL_DW*V2_PDW
                        
                        TESTV2=TESTV2+V2_PDW*(JNUMDW+1) 
                        
                   ELSE
                        DEN_DW=DEN_DW+DNL_DW
                   END IF

            END IF
            
   2        CONTINUE 
                
         END DO !NWSNUM
      
      END DO !LQPRIM
C
C=======================================================================     
C      
      DENGRD_PLOTIN=(DEN_UP+DEN_DW)*FACTOR
C
C=======================================================================     
C
      IF (IF_PAI.EQ.1) THEN
          IF (ABS(TESTV2-N_PART).GT.1.E-4) THEN
              WRITE(LOGFIL,'(''ITRNUM= '',I2)')ITRNUM
              WRITE(LOGFIL,'(''ISOSPI_AUXILI='',I2)')ISOSPI_AUXILI
              WRITE(LOGFIL,'(''In DENSIT_GRADNT: TESTV2='',F20.13)')
     *                                           TESTV2
              WRITE(LOGFIL,'(18X,''N_PART='',I6)')N_PART
     
              IF (ISCREN.GE.5) THEN
                  WRITE(LSCREN,'(''ITRNUM= '',I2)')ITRNUM
                  WRITE(LSCREN,'(''ISOSPI_AUXILI='',I2)')ISOSPI_AUXILI
                  WRITE(LSCREN,'(''In DENSIT_GRADNT: TESTV2='',F20.13)')
     *                                          TESTV2
                  WRITE(LSCREN,'(18X,''N_PART='',I6)')N_PART
              END IF
     
C              STOP 'Stop in DENSIT_GRADNT: TESTV2 not correct!'
          END IF
      END IF
C
C=======================================================================     
C
c      IF (LOGWRI.GT.0) THEN
c          WRITE(LOGFIL,'(/,9X,''Exiting DENSIT_GRADNT'')')
c      END IF
C
C=======================================================================
C
      CALL CPUTIM('DENGRD',0)
C
C=======================================================================
C      
      RETURN
      END
C      
C=======================================================================
C=======================================================================
C
      FUNCTION SPINJR(INUCLI,ISOSPI_AUXILI,ZARGUM,IPOINT,LQNUMB)
C          
      INCLUDE  'MATDIM/NDSPEC.f'
      INCLUDE  'MATDIM/NDIM_N.f'
      INCLUDE  'MATDIM/NDGAUS.f'
      INCLUDE  'MATDIM/NDNUCL.f'
C      
      PARAMETER
     *         (NDIM_L=NDIM_N,NDBASE=NDIM_N+1) 
C      
      DIMENSION
     *          LABORD(0:NDIM_N,0:NDIM_L,0:2*NDIM_L+1)
      DIMENSION
     *          CMATUP(0:NDIM_L,1:NDBASE,1:NDBASE),
     *          CMATDN(0:NDIM_L,1:NDBASE,1:NDBASE)
      DIMENSION
     *          ENERUP(0:NDIM_L,1:NDBASE),
     *          ENERDN(0:NDIM_L,1:NDBASE)
      DIMENSION
     *          VCOEUP(0:NDIM_L,1:NDBASE),
     *          VCOEDN(0:NDIM_L,1:NDBASE)
      DIMENSION
     *          NWSSPH(1:NDSPEC),
     *          LWSSPH(1:NDSPEC),
     *          JWSSPH(1:NDSPEC)
      DIMENSION
     *          DENDER(1:NDGAUS,0:NDIM_N,0:NDIM_N,0:NDIM_L),
     *          AUXILI(1:NDGAUS)
C
      COMMON
     *       /POLDEN/ 
     *        XPOLYN_DENSIT(1:NDGAUS,0:NDIM_N,0:NDIM_L,0:NDIM_L), 
     *        XDPOLY_DENSIT(1:NDGAUS,0:NDIM_N,0:NDIM_L,0:NDIM_L)
      COMMON
     *       /DATA_V/ AOSCIL_PROTON(1:NDNUCL),
     *                AOSCIL_NEUTRS(1:NDNUCL),
     *                HOMEGA(1:NDNUCL)    
      COMMON
     *       /GAUSHE/ NGAUSS,NSHELL_PROTON,NSHELL_NEUTRS
      COMMON
     *       /NMLCNT/ XNORNL(0:NDIM_N,0:NDIM_L)
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
     *       /DIMSPH/ LDSING_PROTON,LEVTHE_PROTON,
     *                LEVEXP_PROTON(1:NDNUCL),
     *                LDSING_NEUTRS,LEVTHE_NEUTRS,
     *                LEVEXP_NEUTRS(1:NDNUCL)
      COMMON
     *       /QUANTN/ NWSSPH_PROTON(1:NDSPEC),
     *      	      LWSSPH_PROTON(1:NDSPEC),
     *      	      JWSSPH_PROTON(1:NDSPEC),
     *
     *      	      NWSSPH_NEUTRS(1:NDSPEC),
     *      	      LWSSPH_NEUTRS(1:NDSPEC),
     *      	      JWSSPH_NEUTRS(1:NDSPEC)
       COMMON
     *       /LEVORD/ LABORD_PROTON(1:NDNUCL,
     *                              0:NDIM_N,0:NDIM_L,0:2*NDIM_L+1),
     *                LABORD_NEUTRS(1:NDNUCL,
     *                              0:NDIM_N,0:NDIM_L,0:2*NDIM_L+1)
      COMMON
     *       /NUCLEU/ IZ_FIX,IN_FIX,ISOSPI 
      COMMON
     *       /IFPAIR/ IF_PAI,IFGPAI,IFDPAI
      COMMON
     *       /BCSENE/ ENEMAX_BCSPAI(1:NDNUCL,1:2)
      COMMON
     *       /BCSDEL/ DELT_2(1:NDNUCL,1:2)
      COMMON
     *       /BCSLAM/ XLAMBD_PROTON(1:NDNUCL),
     *                XLAMBD_NEUTRS(1:NDNUCL)
      COMMON
     *       /DENSIF/ ITRNUM
      COMMON
     *       /PRFILE/ ISCREN,LOGFIL,LOGCHP,LOGCHN,LOGHIP,LOGHIN,
     *                LOGPRO,LOGNEU,LOGWRI,LSCREN,IRESUL,ICONVE,
     *                NOUTPT,LOGAUX,LOGENE,LOGRAD,LOGGRD,LOGMSH,
     *                                                   LOGBIS  
C
C=======================================================================   
C
       CALL CPUTIM('SPINJR',1)
C
C=======================================================================   
C  
CID  This FUNCTION calculates the SPIN CURRENT J(r) over r: (1/r)*J(r)
C
C=======================================================================     
C
c      IF (LOGWRI.GT.0) THEN
c          WRITE(LOGFIL,'(/,9X,''Entering SPINJR'',/,
c     *                   9X,/,9X,''ISOSPI_AUXILI='',I1)')ISOSPI_AUXILI
c      END IF
C
C=======================================================================     
C
C                        P R O T O N S 
C
C=======================================================================    
C
c      WRITE(0,*)'ENTERING IN THE DENSITY FUNCTION'      
                  
      IF (ISOSPI_AUXILI.EQ.1) THEN
        
          N_PART=IZ_FIX
          NSHELL=NSHELL_PROTON
          LDSING=LDSING_PROTON
        
          AOSCIL=AOSCIL_PROTON(INUCLI)
      
          IF (IF_PAI.EQ.1) THEN
              XLAMBD=XLAMBD_PROTON(INUCLI)
              DELTA2=DELT_2(INUCLI,1)
              ENEMAX=ENEMAX_BCSPAI(INUCLI,1)
          END IF
        
          DO INDEXX=1,LDSING
          
             LWSSPH(INDEXX)=LWSSPH_PROTON(INDEXX)
             NWSSPH(INDEXX)=NWSSPH_PROTON(INDEXX)
             JWSSPH(INDEXX)=JWSSPH_PROTON(INDEXX)
          
          END DO
        
          DO INDEXX=1,LDSING
C
             JQNAUX=JWSSPH(INDEXX)
             LQNAUX=LWSSPH(INDEXX)
	     NNBAUX=NWSSPH(INDEXX)
             NBBASE=(NSHELL-LQNAUX)/2  
             
             LABORD(NNBAUX,LQNAUX,JQNAUX)
     *	    =LABORD_PROTON(INUCLI,NNBAUX,LQNAUX,JQNAUX)
C     
             ENERUP(LQNAUX,NNBAUX+1)
     *      =ENERUP_AUXPRO(LQNAUX,NNBAUX+1)          

             ENERDN(LQNAUX,NNBAUX+1)
     *      =ENERDN_AUXPRO(LQNAUX,NNBAUX+1)
C
             DO N1NUMB=0,NBBASE
C		  
                CMATUP(LQNAUX,N1NUMB+1,NNBAUX+1)
     *         =CMATUP_AUXPRO(LQNAUX,N1NUMB+1,NNBAUX+1)
C		  
                CMATDN(LQNAUX,N1NUMB+1,NNBAUX+1)
     *         =CMATDN_AUXPRO(LQNAUX,N1NUMB+1,NNBAUX+1)
C
             END DO
         END DO                  
C
      END IF
C
C=======================================================================     
C
C                        N E U T R O N S 
C
C=======================================================================          
C      
      IF (ISOSPI_AUXILI.EQ.0) THEN
        
          N_PART=IN_FIX
          NSHELL=NSHELL_NEUTRS
          LDSING=LDSING_NEUTRS
        
          AOSCIL=AOSCIL_NEUTRS(INUCLI)
C
          IF (IF_PAI.EQ.1) THEN
              XLAMBD=XLAMBD_NEUTRS(INUCLI)
              DELTA2=DELT_2(INUCLI,2)
              ENEMAX=ENEMAX_BCSPAI(INUCLI,2)
          END IF
        
          DO INDEXX=1,LDSING
          
             LWSSPH(INDEXX)=LWSSPH_NEUTRS(INDEXX)
             NWSSPH(INDEXX)=NWSSPH_NEUTRS(INDEXX)
             JWSSPH(INDEXX)=JWSSPH_NEUTRS(INDEXX)
          
          END DO
        
          DO INDEXX=1,LDSING
C
             JQNAUX=JWSSPH(INDEXX)
             LQNAUX=LWSSPH(INDEXX)
	     NNBAUX=NWSSPH(INDEXX)
             NBBASE=(NSHELL-LQNAUX)/2  
             
             LABORD(NNBAUX,LQNAUX,JQNAUX)
     *	    =LABORD_NEUTRS(INUCLI,NNBAUX,LQNAUX,JQNAUX)
C     
             ENERUP(LQNAUX,NNBAUX+1)
     *      =ENERUP_AUXNEU(LQNAUX,NNBAUX+1)          
C
             ENERDN(LQNAUX,NNBAUX+1)
     *      =ENERDN_AUXNEU(LQNAUX,NNBAUX+1)
C     
             DO N1NUMB=0,NBBASE
C		  
                CMATUP(LQNAUX,N1NUMB+1,NNBAUX+1)
     *         =CMATUP_AUXNEU(LQNAUX,N1NUMB+1,NNBAUX+1)
C		  
                CMATDN(LQNAUX,N1NUMB+1,NNBAUX+1)
     *         =CMATDN_AUXNEU(LQNAUX,N1NUMB+1,NNBAUX+1)
     
C
             END DO
         END DO                  
C
      END IF
C
C=======================================================================     
C
C                     Calculating J(r)/r
C
C=======================================================================        
C   
      PINUMB=4.*ATAN(1.0)
      FACTOR=1./(PINUMB*AOSCIL**5)
C
C  Here we calculate J_UP block
C       
      TESTV2=0.0
      XJR_UP=0.0
      
      DO LQPRIM=NSHELL,0,-1 
         
         JNUMUP=2*LQPRIM+1
          
         DO NWSNUM=1,(NSHELL-LQPRIM)/2+1
             
            ENERGY=ENERUP(LQPRIM,NWSNUM)
            
            IF (IF_PAI.EQ.1.AND.ENERGY.GT.ENEMAX) GO TO 1
           
            IF (IF_PAI.EQ.0
     *          .AND.LABORD(NWSNUM-1,LQPRIM,JNUMUP).EQ.0) GO TO 1
            
            XNL_UP=0.0
            
            DO N1NUMB=0,(NSHELL-LQPRIM)/2
                
               DO N2NUMB=0,(NSHELL-LQPRIM)/2
                  
                  CFACT1=LQPRIM*(LQPRIM+1)
     *                  *CMATUP(LQPRIM,N1NUMB+1,NWSNUM)
     *                  *CMATUP(LQPRIM,N2NUMB+1,NWSNUM)
              
                  POLPOL=XPOLYN_DENSIT(IPOINT,N1NUMB,LQPRIM,LQNUMB) 
     *                  *XPOLYN_DENSIT(IPOINT,N2NUMB,LQPRIM,LQNUMB)
                
                  XNL_UP=XNL_UP+EXP(-ZARGUM)*ZARGUM**(LQPRIM-1)
     *                         *CFACT1*POLPOL
                  
               END DO !N2NUMB
              
            END DO !N1NUMB
            
            IF (IF_PAI.EQ.1) THEN
                V2_PUP=V2PAIR(ENERGY,XLAMBD,DELTA2)
                XJR_UP=XJR_UP+XNL_UP*V2_PUP  
                TESTV2=TESTV2+V2_PUP*(JNUMUP+1)   
            ELSE
                XJR_UP=XJR_UP+XNL_UP
            END IF
            
   1        CONTINUE
            
         END DO !NWSNUM
      
      END DO !LQPRIM
C
C  Here we calculate J_DOWN block
C      
      XJR_DW=0.0
      
      DO LQPRIM=NSHELL,0,-1 
         
         JNUMDW=2*LQPRIM-1
          
         DO NWSNUM=1,(NSHELL-LQPRIM)/2+1
             
            ENERGY=ENERDN(LQPRIM,NWSNUM)
            
            IF (JNUMDW.GE.0) THEN
                
                IF (IF_PAI.EQ.1.AND.ENERGY.GT.ENEMAX) GO TO 2
           
                IF (IF_PAI.EQ.0
     *              .AND.LABORD(NWSNUM-1,LQPRIM,JNUMDW).EQ.0) GO TO 2
                    
                XNL_DW=0.0
            
                DO N1NUMB=0,(NSHELL-LQPRIM)/2
                
                   DO N2NUMB=0,(NSHELL-LQPRIM)/2
    
                      CFACT2=(-LQPRIM)*(LQPRIM+1)
     *                      *CMATDN(LQPRIM,N1NUMB+1,NWSNUM)
     *                      *CMATDN(LQPRIM,N2NUMB+1,NWSNUM)
              
                      POLPOL=XPOLYN_DENSIT(IPOINT,N1NUMB,LQPRIM,LQNUMB) 
     *                      *XPOLYN_DENSIT(IPOINT,N2NUMB,LQPRIM,LQNUMB)
                
                      XNL_DW=XNL_DW+EXP(-ZARGUM)*ZARGUM**(LQPRIM-1)
     *                             *CFACT2*POLPOL
                
                   END DO !N2NUMB
              
                END DO !N1NUMB
                
                IF (IF_PAI.EQ.1) THEN
                    V2_PDW=V2PAIR(ENERGY,XLAMBD,DELTA2)
                    XJR_DW=XJR_DW+XNL_DW*V2_PDW
                    TESTV2=TESTV2+V2_PDW*(JNUMDW+1) 
                ELSE
                    XJR_DW=XJR_DW+XNL_DW
                END IF

            END IF
            
   2        CONTINUE 
            
        END DO !NWSNUM
      
      END DO !LQPRIM
      
      SPINJR=(XJR_UP+XJR_DW)*FACTOR
C
C=======================================================================     
C
      IF (IF_PAI.EQ.1) THEN
          IF (ABS(TESTV2-N_PART).GT.1.E-4) THEN
              WRITE(LOGFIL,'(''ITRNUM= '',I2)')ITRNUM
              WRITE(LOGFIL,'(''ISOSPI_AUXILI='',I2)')ISOSPI_AUXILI
              WRITE(LOGFIL,'(''In SPINJR: TESTV2='',F20.13)')TESTV2
              WRITE(LOGFIL,'(11X,''N_PART='',I6)')N_PART
     
              IF (ISCREN.GE.5) THEN
                  WRITE(0,'(''ITRNUM= '',I2)')ITRNUM
                  WRITE(0,'(''ISOSPI_AUXILI='',I2)')ISOSPI_AUXILI
                  WRITE(0,'(''In SPINJR: TESTV2='',F20.13)')TESTV2
                  WRITE(0,'(11X,''N_PART='',I6)')N_PART
              END IF
     
C              STOP 'Stop in DENSIT_GRADNT: TESTV2 not correct!'
          END IF
      END IF
C
C=======================================================================     
C
c      IF (LOGWRI.GT.0) THEN
c          WRITE(LOGFIL,'(/,9X,''Exiting SPINJR'')')
c      END IF
C
C=======================================================================
C
      CALL CPUTIM('SPINJR',0)
C
C=======================================================================
C      
      RETURN
      END   
C      
C=======================================================================
C=======================================================================
C
      FUNCTION SPIN_J(INUCLI,ISOSPI_AUXILI,ZARGUM,IPOINT,LQNUMB)
C          
      INCLUDE  'MATDIM/NDSPEC.f'
      INCLUDE  'MATDIM/NDIM_N.f'
      INCLUDE  'MATDIM/NDGAUS.f'
      INCLUDE  'MATDIM/NDNUCL.f'
C      
      PARAMETER
     *         (NDIM_L=NDIM_N,NDBASE=NDIM_N+1) 
C      
      DIMENSION
     *          LABORD(0:NDIM_N,0:NDIM_L,0:2*NDIM_L+1)
      DIMENSION
     *          CMATUP(0:NDIM_L,1:NDBASE,1:NDBASE),
     *          CMATDN(0:NDIM_L,1:NDBASE,1:NDBASE)
      DIMENSION
     *          ENERUP(0:NDIM_L,1:NDBASE),
     *          ENERDN(0:NDIM_L,1:NDBASE)
      DIMENSION
     *          VCOEUP(0:NDIM_L,1:NDBASE),
     *          VCOEDN(0:NDIM_L,1:NDBASE)
      DIMENSION
     *          NWSSPH(1:NDSPEC),
     *          LWSSPH(1:NDSPEC),
     *          JWSSPH(1:NDSPEC)
      DIMENSION
     *          DENDER(1:NDGAUS,0:NDIM_N,0:NDIM_N,0:NDIM_L),
     *          AUXILI(1:NDGAUS)
C
      COMMON
     *       /POLDEN/ 
     *        XPOLYN_DENSIT(1:NDGAUS,0:NDIM_N,0:NDIM_L,0:NDIM_L), 
     *        XDPOLY_DENSIT(1:NDGAUS,0:NDIM_N,0:NDIM_L,0:NDIM_L)
      COMMON
     *       /DATA_V/ AOSCIL_PROTON(1:NDNUCL),
     *                AOSCIL_NEUTRS(1:NDNUCL),
     *                HOMEGA(1:NDNUCL)    
      COMMON
     *       /GAUSHE/ NGAUSS,NSHELL_PROTON,NSHELL_NEUTRS
      COMMON
     *       /NMLCNT/ XNORNL(0:NDIM_N,0:NDIM_L)
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
     *       /DIMSPH/ LDSING_PROTON,LEVTHE_PROTON,
     *                LEVEXP_PROTON(1:NDNUCL),
     *                LDSING_NEUTRS,LEVTHE_NEUTRS,
     *                LEVEXP_NEUTRS(1:NDNUCL)
      COMMON
     *       /QUANTN/ NWSSPH_PROTON(1:NDSPEC),
     *      	      LWSSPH_PROTON(1:NDSPEC),
     *      	      JWSSPH_PROTON(1:NDSPEC),
     *
     *      	      NWSSPH_NEUTRS(1:NDSPEC),
     *      	      LWSSPH_NEUTRS(1:NDSPEC),
     *      	      JWSSPH_NEUTRS(1:NDSPEC)
       COMMON
     *       /LEVORD/ LABORD_PROTON(1:NDNUCL,
     *                              0:NDIM_N,0:NDIM_L,0:2*NDIM_L+1),
     *                LABORD_NEUTRS(1:NDNUCL,
     *                              0:NDIM_N,0:NDIM_L,0:2*NDIM_L+1)
      COMMON
     *       /NUCLEU/ IZ_FIX,IN_FIX,ISOSPI 
      COMMON
     *       /IFPAIR/ IF_PAI,IFGPAI,IFDPAI
      COMMON
     *       /BCSENE/ ENEMAX_BCSPAI(1:NDNUCL,1:2)
      COMMON
     *       /BCSDEL/ DELT_2(1:NDNUCL,1:2)
      COMMON
     *       /BCSLAM/ XLAMBD_PROTON(1:NDNUCL),
     *                XLAMBD_NEUTRS(1:NDNUCL)
      COMMON
     *       /DENSIF/ ITRNUM
      COMMON
     *       /PRFILE/ ISCREN,LOGFIL,LOGCHP,LOGCHN,LOGHIP,LOGHIN,
     *                LOGPRO,LOGNEU,LOGWRI,LSCREN,IRESUL,ICONVE,
     *                NOUTPT,LOGAUX,LOGENE,LOGRAD,LOGGRD,LOGMSH,
     *                                                   LOGBIS  
C
C=======================================================================   
C
       CALL CPUTIM('SPIN_J',1)
C
C=======================================================================   
C  
CID        This FUNCTION calculates the SPIN CURRENT J(r) 
C
C=======================================================================     
C
c      IF (LOGWRI.GT.0) THEN
c          WRITE(LOGFIL,'(/,9X,''Entering SPIN_J'',/,
c     *                   9X,/,9X,''ISOSPI_AUXILI='',I1)')ISOSPI_AUXILI
c      END IF
C
C=======================================================================     
C
C                        P R O T O N S 
C
C=======================================================================    
C
c      WRITE(0,*)'ENTERING IN THE DENSITY FUNCTION'      
                  
      IF (ISOSPI_AUXILI.EQ.1) THEN
        
        N_PART=IZ_FIX
        NSHELL=NSHELL_PROTON
        LDSING=LDSING_PROTON
        
        AOSCIL=AOSCIL_PROTON(INUCLI)
C
        IF (IF_PAI.EQ.1) THEN
              XLAMBD=XLAMBD_PROTON(INUCLI)
              DELTA2=DELT_2(INUCLI,1)
              ENEMAX=ENEMAX_BCSPAI(INUCLI,1)
          END IF
        
        DO INDEXX=1,LDSING
          
          LWSSPH(INDEXX)=LWSSPH_PROTON(INDEXX)
          NWSSPH(INDEXX)=NWSSPH_PROTON(INDEXX)
          JWSSPH(INDEXX)=JWSSPH_PROTON(INDEXX)
          
        END DO
        
        DO INDEXX=1,LDSING
C
           JQNAUX=JWSSPH(INDEXX)
           LQNAUX=LWSSPH(INDEXX)
	   NNBAUX=NWSSPH(INDEXX)
           NBBASE=(NSHELL-LQNAUX)/2  
             
           LABORD(NNBAUX,LQNAUX,JQNAUX)
     *	  =LABORD_PROTON(INUCLI,NNBAUX,LQNAUX,JQNAUX)
             
           ENERUP(LQNAUX,NNBAUX+1)
     *    =ENERUP_AUXPRO(LQNAUX,NNBAUX+1)          

           ENERDN(LQNAUX,NNBAUX+1)
     *    =ENERDN_AUXPRO(LQNAUX,NNBAUX+1)
C
           DO N1NUMB=0,NBBASE
C		  
              CMATUP(LQNAUX,N1NUMB+1,NNBAUX+1)
     *       =CMATUP_AUXPRO(LQNAUX,N1NUMB+1,NNBAUX+1)
C		  
              CMATDN(LQNAUX,N1NUMB+1,NNBAUX+1)
     *       =CMATDN_AUXPRO(LQNAUX,N1NUMB+1,NNBAUX+1)
C
           END DO
        END DO                  
C
      END IF
C
C=======================================================================     
C
C                        N E U T R O N S 
C
C=======================================================================          
C      
      IF (ISOSPI_AUXILI.EQ.0) THEN
        
        N_PART=IN_FIX
        NSHELL=NSHELL_NEUTRS
        LDSING=LDSING_NEUTRS
        
        AOSCIL=AOSCIL_NEUTRS(INUCLI)
C     
        IF (IF_PAI.EQ.1) THEN
              XLAMBD=XLAMBD_NEUTRS(INUCLI)
              DELTA2=DELT_2(INUCLI,2)
              ENEMAX=ENEMAX_BCSPAI(INUCLI,2)
          END IF
        
        DO INDEXX=1,LDSING
          
          LWSSPH(INDEXX)=LWSSPH_NEUTRS(INDEXX)
          NWSSPH(INDEXX)=NWSSPH_NEUTRS(INDEXX)
          JWSSPH(INDEXX)=JWSSPH_NEUTRS(INDEXX)
          
        END DO
        
        DO INDEXX=1,LDSING
C
           JQNAUX=JWSSPH(INDEXX)
           LQNAUX=LWSSPH(INDEXX)
	   NNBAUX=NWSSPH(INDEXX)
           NBBASE=(NSHELL-LQNAUX)/2  
             
           LABORD(NNBAUX,LQNAUX,JQNAUX)
     *	  =LABORD_NEUTRS(INUCLI,NNBAUX,LQNAUX,JQNAUX)
C     
           ENERUP(LQNAUX,NNBAUX+1)
     *    =ENERUP_AUXNEU(LQNAUX,NNBAUX+1)          
C
           ENERDN(LQNAUX,NNBAUX+1)
     *    =ENERDN_AUXNEU(LQNAUX,NNBAUX+1)
C     
           DO N1NUMB=0,NBBASE
C		  
              CMATUP(LQNAUX,N1NUMB+1,NNBAUX+1)
     *       =CMATUP_AUXNEU(LQNAUX,N1NUMB+1,NNBAUX+1)
C		  
              CMATDN(LQNAUX,N1NUMB+1,NNBAUX+1)
     *       =CMATDN_AUXNEU(LQNAUX,N1NUMB+1,NNBAUX+1)
C
             END DO
        END DO                  
C
      END IF
C
C=======================================================================     
C
C                          Calculating J(r)
C
C=======================================================================        
C   
      PINUMB=4.*ATAN(1.0)
      FACTOR=1./(PINUMB*AOSCIL**4)      
C
C  Here we calculate J_UP block
C       
      TESTV2=0.0
      XJR_UP=0.0
C      
      DO LQPRIM=NSHELL,0,-1 
C         
         JNUMUP=2*LQPRIM+1
C          
         DO NWSNUM=1,(NSHELL-LQPRIM)/2+1
C             
            ENERGY=ENERUP(LQPRIM,NWSNUM)
C            
            IF (IF_PAI.EQ.1.AND.ENERGY.GT.ENEMAX) GO TO 1
C           
            IF (IF_PAI.EQ.0
     *          .AND.LABORD(NWSNUM-1,LQPRIM,JNUMUP).EQ.0) GO TO 1
C            
            XNL_UP=0.0
C            
            DO N1NUMB=0,(NSHELL-LQPRIM)/2
C                
               DO N2NUMB=0,(NSHELL-LQPRIM)/2
C                  
                  CFACT1=LQPRIM*(LQPRIM+1)
     *                  *CMATUP(LQPRIM,N1NUMB+1,NWSNUM)
     *                  *CMATUP(LQPRIM,N2NUMB+1,NWSNUM)
C              
                  POLPOL=XPOLYN_DENSIT(IPOINT,N1NUMB,LQPRIM,LQNUMB) 
     *                  *XPOLYN_DENSIT(IPOINT,N2NUMB,LQPRIM,LQNUMB)
C                
                  XNL_UP=XNL_UP+EXP(-ZARGUM)*ZARGUM**(LQPRIM-0.5)
     *                         *CFACT1*POLPOL
C                  
               END DO !N2NUMB
C              
            END DO !N1NUMB
C            
            IF (IF_PAI.EQ.1) THEN
                V2_PUP=V2PAIR(ENERGY,XLAMBD,DELTA2)
                XJR_UP=XJR_UP+XNL_UP*V2_PUP  
                TESTV2=TESTV2+V2_PUP*(JNUMUP+1)   
            ELSE
                XJR_UP=XJR_UP+XNL_UP
            END IF
C            
   1        CONTINUE
C            
         END DO !NWSNUM
C      
      END DO !LQPRIM
C
C  Here we calculate J_DOWN block
C      
      XJR_DW=0.0
C      
      DO LQPRIM=NSHELL,0,-1 
C         
         JNUMDW=2*LQPRIM-1
C
         IF (JNUMDW.GE.0) THEN
C          
             DO NWSNUM=1,(NSHELL-LQPRIM)/2+1
C             
                ENERGY=ENERDN(LQPRIM,NWSNUM)
C                
                IF (IF_PAI.EQ.1.AND.ENERGY.GT.ENEMAX) GO TO 2
C           
                IF (IF_PAI.EQ.0
     *              .AND.LABORD(NWSNUM-1,LQPRIM,JNUMDW).EQ.0) GO TO 2
C                    
                XNL_DW=0.0
C            
                DO N1NUMB=0,(NSHELL-LQPRIM)/2
C                
                   DO N2NUMB=0,(NSHELL-LQPRIM)/2
C    
                      CFACT2=(-LQPRIM)*(LQPRIM+1)
     *                      *CMATDN(LQPRIM,N1NUMB+1,NWSNUM)
     *                      *CMATDN(LQPRIM,N2NUMB+1,NWSNUM)
C              
                      POLPOL=XPOLYN_DENSIT(IPOINT,N1NUMB,LQPRIM,LQNUMB) 
     *                      *XPOLYN_DENSIT(IPOINT,N2NUMB,LQPRIM,LQNUMB)
C                
                      XNL_DW=XNL_DW+EXP(-ZARGUM)*ZARGUM**(LQPRIM-0.5)
     *                             *CFACT2*POLPOL
C                
                   END DO !N2NUMB
C              
                END DO !N1NUMB
C                
                IF (IF_PAI.EQ.1) THEN
                    V2_PDW=V2PAIR(ENERGY,XLAMBD,DELTA2)
                    XJR_DW=XJR_DW+XNL_DW*V2_PDW
                    TESTV2=TESTV2+V2_PDW*(JNUMDW+1) 
                ELSE
                    XJR_DW=XJR_DW+XNL_DW
                END IF
C            
   2            CONTINUE 
C            
             END DO !NWSNUM
C
         END IF
C      
      END DO !LQPRIM
C      
      SPIN_J=(XJR_UP+XJR_DW)*FACTOR
C
C=======================================================================     
C
      IF (IF_PAI.EQ.1) THEN
          IF (ABS(TESTV2-N_PART).GT.1.E-4) THEN
              WRITE(LOGFIL,'(''ITRNUM= '',I2)')ITRNUM
              WRITE(LOGFIL,'(''ISOSPI_AUXILI='',I2)')ISOSPI_AUXILI
              WRITE(LOGFIL,'(''In SPIN_J: TESTV2='',F20.13)')TESTV2
              WRITE(LOGFIL,'(11X,''N_PART='',I6)')N_PART
     
              IF (ISCREN.GE.5) THEN
                  WRITE(0,'(''ITRNUM= '',I2)')ITRNUM
                  WRITE(0,'(''ISOSPI_AUXILI='',I2)')ISOSPI_AUXILI
                  WRITE(0,'(''In SPIN_J: TESTV2='',F20.13)')TESTV2
                  WRITE(0,'(11X,''N_PART='',I6)')N_PART
              END IF
     
C              STOP 'Stop in DENSIT_GRADNT: TESTV2 not correct!'
          END IF
      END IF
C
C=======================================================================     
C
c      IF (LOGWRI.GT.0) THEN
c          WRITE(LOGFIL,'(/,9X,''Exiting SPIN_J'')')
c      END IF
C
C=======================================================================
C
      CALL CPUTIM('SPIN_J',0)
C
C=======================================================================
C      
      RETURN
      END   
C      
C=======================================================================
C=======================================================================
C
      FUNCTION SPIN_J_PLOTIN(INUCLI,ISOSPI_AUXILI,ZARGUM)
C          
      INCLUDE  'MATDIM/NDSPEC.f'
      INCLUDE  'MATDIM/NDIM_N.f'
      INCLUDE  'MATDIM/NDGAUS.f'
      INCLUDE  'MATDIM/NDNUCL.f'
C      
      PARAMETER
     *         (NDIM_L=NDIM_N,NDBASE=NDIM_N+1) 
C      
      REAL*16
     *          QARGUM,QFLAGN,QDFLAG,QNORNL
C
      DIMENSION
     *          LABORD(0:NDIM_N,0:NDIM_L,0:2*NDIM_L+1)
      DIMENSION
     *          CMATUP(0:NDIM_L,1:NDBASE,1:NDBASE),
     *          CMATDN(0:NDIM_L,1:NDBASE,1:NDBASE)
      DIMENSION
     *          ENERUP(0:NDIM_L,1:NDBASE),
     *          ENERDN(0:NDIM_L,1:NDBASE)
      DIMENSION
     *          VCOEUP(0:NDIM_L,1:NDBASE),
     *          VCOEDN(0:NDIM_L,1:NDBASE)
      DIMENSION
     *          NWSSPH(1:NDSPEC),
     *          LWSSPH(1:NDSPEC),
     *          JWSSPH(1:NDSPEC)
      DIMENSION
     *          DENDER(1:NDGAUS,0:NDIM_N,0:NDIM_N,0:NDIM_L),
     *          AUXILI(1:NDGAUS)
      DIMENSION
     *          QFLAGN(0:NDIM_N,0:NDIM_L),
     *          QDFLAG(0:NDIM_N,0:NDIM_L),
     *          XFLAGN(0:NDIM_N,0:NDIM_L),
     *          XDFLAG(0:NDIM_N,0:NDIM_L)
C
      COMMON
     *       /POLDEN/ 
     *        XPOLYN_DENSIT(1:NDGAUS,0:NDIM_N,0:NDIM_L,0:NDIM_L), 
     *        XDPOLY_DENSIT(1:NDGAUS,0:NDIM_N,0:NDIM_L,0:NDIM_L)
      COMMON
     *       /DATA_V/ AOSCIL_PROTON(1:NDNUCL),
     *                AOSCIL_NEUTRS(1:NDNUCL),
     *                HOMEGA(1:NDNUCL)    
      COMMON
     *       /GAUSHE/ NGAUSS,NSHELL_PROTON,NSHELL_NEUTRS
      COMMON
     *       /NMLCNT/ XNORNL(0:NDIM_N,0:NDIM_L)
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
     *       /DIMSPH/ LDSING_PROTON,LEVTHE_PROTON,
     *                LEVEXP_PROTON(1:NDNUCL),
     *                LDSING_NEUTRS,LEVTHE_NEUTRS,
     *                LEVEXP_NEUTRS(1:NDNUCL)
      COMMON
     *       /QUANTN/ NWSSPH_PROTON(1:NDSPEC),
     *      	      LWSSPH_PROTON(1:NDSPEC),
     *      	      JWSSPH_PROTON(1:NDSPEC),
     *
     *      	      NWSSPH_NEUTRS(1:NDSPEC),
     *      	      LWSSPH_NEUTRS(1:NDSPEC),
     *      	      JWSSPH_NEUTRS(1:NDSPEC)
       COMMON
     *       /LEVORD/ LABORD_PROTON(1:NDNUCL,
     *                              0:NDIM_N,0:NDIM_L,0:2*NDIM_L+1),
     *                LABORD_NEUTRS(1:NDNUCL,
     *                              0:NDIM_N,0:NDIM_L,0:2*NDIM_L+1)
      COMMON
     *       /NUCLEU/ IZ_FIX,IN_FIX,ISOSPI 
      COMMON
     *       /IFPAIR/ IF_PAI,IFGPAI,IFDPAI
      COMMON
     *       /BCSENE/ ENEMAX_BCSPAI(1:NDNUCL,1:2)
      COMMON
     *       /BCSDEL/ DELT_2(1:NDNUCL,1:2)
      COMMON
     *       /BCSLAM/ XLAMBD_PROTON(1:NDNUCL),
     *                XLAMBD_NEUTRS(1:NDNUCL)
      COMMON
     *       /DENSIF/ ITRNUM
      COMMON
     *       /PRFILE/ ISCREN,LOGFIL,LOGCHP,LOGCHN,LOGHIP,LOGHIN,
     *                LOGPRO,LOGNEU,LOGWRI,LSCREN,IRESUL,ICONVE,
     *                NOUTPT,LOGAUX,LOGENE,LOGRAD,LOGGRD,LOGMSH,
     *                                                   LOGBIS  
C
C=======================================================================   
C
       CALL CPUTIM('SPIN_J',1)
C
C=======================================================================   
C  
CID        This FUNCTION calculates the SPIN CURRENT J(r) 
C
C=======================================================================     
C
c      IF (LOGWRI.GT.0) THEN
c          WRITE(LOGFIL,'(/,9X,''Entering SPIN_J'',/,
c     *                   9X,/,9X,''ISOSPI_AUXILI='',I1)')ISOSPI_AUXILI
c      END IF
C
C=======================================================================     
C
C                        P R O T O N S 
C
C=======================================================================    
C
c      WRITE(0,*)'ENTERING IN THE DENSITY FUNCTION'      
                  
      IF (ISOSPI_AUXILI.EQ.1) THEN
        
        N_PART=IZ_FIX
        NSHELL=NSHELL_PROTON
        LDSING=LDSING_PROTON
        
        AOSCIL=AOSCIL_PROTON(INUCLI)
C
        IF (IF_PAI.EQ.1) THEN
              XLAMBD=XLAMBD_PROTON(INUCLI)
              DELTA2=DELT_2(INUCLI,1)
              ENEMAX=ENEMAX_BCSPAI(INUCLI,1)
          END IF
        
        DO INDEXX=1,LDSING
          
          LWSSPH(INDEXX)=LWSSPH_PROTON(INDEXX)
          NWSSPH(INDEXX)=NWSSPH_PROTON(INDEXX)
          JWSSPH(INDEXX)=JWSSPH_PROTON(INDEXX)
          
        END DO
        
        DO INDEXX=1,LDSING
C
           JQNAUX=JWSSPH(INDEXX)
           LQNAUX=LWSSPH(INDEXX)
	   NNBAUX=NWSSPH(INDEXX)
           NBBASE=(NSHELL-LQNAUX)/2  
             
           LABORD(NNBAUX,LQNAUX,JQNAUX)
     *	  =LABORD_PROTON(INUCLI,NNBAUX,LQNAUX,JQNAUX)
             
           ENERUP(LQNAUX,NNBAUX+1)
     *    =ENERUP_AUXPRO(LQNAUX,NNBAUX+1)          

           ENERDN(LQNAUX,NNBAUX+1)
     *    =ENERDN_AUXPRO(LQNAUX,NNBAUX+1)
C
           DO N1NUMB=0,NBBASE
C		  
              CMATUP(LQNAUX,N1NUMB+1,NNBAUX+1)
     *       =CMATUP_AUXPRO(LQNAUX,N1NUMB+1,NNBAUX+1)
C		  
              CMATDN(LQNAUX,N1NUMB+1,NNBAUX+1)
     *       =CMATDN_AUXPRO(LQNAUX,N1NUMB+1,NNBAUX+1)
C
           END DO
        END DO                  
C
      END IF
C
C=======================================================================     
C
C                        N E U T R O N S 
C
C=======================================================================          
C      
      IF (ISOSPI_AUXILI.EQ.0) THEN
        
        N_PART=IN_FIX
        NSHELL=NSHELL_NEUTRS
        LDSING=LDSING_NEUTRS
        
        AOSCIL=AOSCIL_NEUTRS(INUCLI)
C     
        IF (IF_PAI.EQ.1) THEN
              XLAMBD=XLAMBD_NEUTRS(INUCLI)
              DELTA2=DELT_2(INUCLI,2)
              ENEMAX=ENEMAX_BCSPAI(INUCLI,2)
          END IF
        
        DO INDEXX=1,LDSING
          
          LWSSPH(INDEXX)=LWSSPH_NEUTRS(INDEXX)
          NWSSPH(INDEXX)=NWSSPH_NEUTRS(INDEXX)
          JWSSPH(INDEXX)=JWSSPH_NEUTRS(INDEXX)
          
        END DO
        
        DO INDEXX=1,LDSING
C
           JQNAUX=JWSSPH(INDEXX)
           LQNAUX=LWSSPH(INDEXX)
	   NNBAUX=NWSSPH(INDEXX)
           NBBASE=(NSHELL-LQNAUX)/2  
             
           LABORD(NNBAUX,LQNAUX,JQNAUX)
     *	  =LABORD_NEUTRS(INUCLI,NNBAUX,LQNAUX,JQNAUX)
C     
           ENERUP(LQNAUX,NNBAUX+1)
     *    =ENERUP_AUXNEU(LQNAUX,NNBAUX+1)          
C
           ENERDN(LQNAUX,NNBAUX+1)
     *    =ENERDN_AUXNEU(LQNAUX,NNBAUX+1)
C     
           DO N1NUMB=0,NBBASE
C		  
              CMATUP(LQNAUX,N1NUMB+1,NNBAUX+1)
     *       =CMATUP_AUXNEU(LQNAUX,N1NUMB+1,NNBAUX+1)
C		  
              CMATDN(LQNAUX,N1NUMB+1,NNBAUX+1)
     *       =CMATDN_AUXNEU(LQNAUX,N1NUMB+1,NNBAUX+1)
C
             END DO
        END DO                  
C
      END IF

C
C=======================================================================     
C
C       Calculating the polynomials in the desired ZARGUM point
C
C=======================================================================          
C      
      KIND_L=1
      NORDER=NSHELL/2
      LORDER=NSHELL
      
      QARGUM=QEXT(ZARGUM)    !passing the argument to real*16
      
      CALL LAGUER(QARGUM,QFLAGN,QDFLAG,NORDER,LORDER,KIND_L,
     *                                        NDIM_L,NDIM_N)
     
      DO NAUXIL=0,NORDER
        DO LAUXIL=0,LORDER
            XFLAGN(NAUXIL,LAUXIL)=REAL(QFLAGN(NAUXIL,LAUXIL))
            XDFLAG(NAUXIL,LAUXIL)=REAL(QDFLAG(NAUXIL,LAUXIL))
        END DO
      END DO
C
C=======================================================================     
C
C                          Calculating J(r)
C
C=======================================================================        
C   
      PINUMB=4.*ATAN(1.0)
      FACTOR=1./(PINUMB*AOSCIL**4)      
C
C  Here we calculate J_UP block
C       
      TESTV2=0.0
      XJR_UP=0.0
C      
      DO LQPRIM=NSHELL,0,-1 
C         
         JNUMUP=2*LQPRIM+1
C          
         DO NWSNUM=1,(NSHELL-LQPRIM)/2+1
C             
            ENERGY=ENERUP(LQPRIM,NWSNUM)
C            
            IF (IF_PAI.EQ.1.AND.ENERGY.GT.ENEMAX) GO TO 1
C           
            IF (IF_PAI.EQ.0
     *          .AND.LABORD(NWSNUM-1,LQPRIM,JNUMUP).EQ.0) GO TO 1
C            
            XNL_UP=0.0
C            
            DO N1NUMB=0,(NSHELL-LQPRIM)/2
C                
               DO N2NUMB=0,(NSHELL-LQPRIM)/2
C                  
                  CFACT1=LQPRIM*(LQPRIM+1)
     *                  *CMATUP(LQPRIM,N1NUMB+1,NWSNUM)
     *                  *CMATUP(LQPRIM,N2NUMB+1,NWSNUM)
C              
                  POLPOL=XNORNL(N1NUMB,LQPRIM)*XFLAGN(N1NUMB,LQPRIM)
     *                  *XNORNL(N2NUMB,LQPRIM)*XFLAGN(N2NUMB,LQPRIM)
C                
                  XNL_UP=XNL_UP+EXP(-ZARGUM)*ZARGUM**(LQPRIM-0.5)
     *                         *CFACT1*POLPOL
C                  
               END DO !N2NUMB
C              
            END DO !N1NUMB
C            
            IF (IF_PAI.EQ.1) THEN
                V2_PUP=V2PAIR(ENERGY,XLAMBD,DELTA2)
                XJR_UP=XJR_UP+XNL_UP*V2_PUP  
                TESTV2=TESTV2+V2_PUP*(JNUMUP+1)   
            ELSE
                XJR_UP=XJR_UP+XNL_UP
            END IF
C            
   1        CONTINUE
C            
         END DO !NWSNUM
C      
      END DO !LQPRIM
C
C  Here we calculate J_DOWN block
C      
      XJR_DW=0.0
C      
      DO LQPRIM=NSHELL,0,-1 
C         
         JNUMDW=2*LQPRIM-1
C
         IF (JNUMDW.GE.0) THEN
C          
             DO NWSNUM=1,(NSHELL-LQPRIM)/2+1
C             
                ENERGY=ENERDN(LQPRIM,NWSNUM)
C                
                IF (IF_PAI.EQ.1.AND.ENERGY.GT.ENEMAX) GO TO 2
C           
                IF (IF_PAI.EQ.0
     *              .AND.LABORD(NWSNUM-1,LQPRIM,JNUMDW).EQ.0) GO TO 2
C                    
                XNL_DW=0.0
C            
                DO N1NUMB=0,(NSHELL-LQPRIM)/2
C                
                   DO N2NUMB=0,(NSHELL-LQPRIM)/2
C    
                      CFACT2=(-LQPRIM)*(LQPRIM+1)
     *                      *CMATDN(LQPRIM,N1NUMB+1,NWSNUM)
     *                      *CMATDN(LQPRIM,N2NUMB+1,NWSNUM)
C              
                      POLPOL=XNORNL(N1NUMB,LQPRIM)*XFLAGN(N1NUMB,LQPRIM)
     *                      *XNORNL(N2NUMB,LQPRIM)*XFLAGN(N2NUMB,LQPRIM)
C                
                      XNL_DW=XNL_DW+EXP(-ZARGUM)*ZARGUM**(LQPRIM-0.5)
     *                             *CFACT2*POLPOL
C                
                   END DO !N2NUMB
C              
                END DO !N1NUMB
C                
                IF (IF_PAI.EQ.1) THEN
                    V2_PDW=V2PAIR(ENERGY,XLAMBD,DELTA2)
                    XJR_DW=XJR_DW+XNL_DW*V2_PDW
                    TESTV2=TESTV2+V2_PDW*(JNUMDW+1) 
                ELSE
                    XJR_DW=XJR_DW+XNL_DW
                END IF
C            
   2            CONTINUE 
C            
             END DO !NWSNUM
C
         END IF
C      
      END DO !LQPRIM
C      
      SPIN_J_PLOTIN=(XJR_UP+XJR_DW)*FACTOR
C
C=======================================================================     
C
      IF (IF_PAI.EQ.1) THEN
          IF (ABS(TESTV2-N_PART).GT.1.E-4) THEN
              WRITE(LOGFIL,'(''ITRNUM= '',I2)')ITRNUM
              WRITE(LOGFIL,'(''ISOSPI_AUXILI='',I2)')ISOSPI_AUXILI
              WRITE(LOGFIL,'(''In SPIN_J: TESTV2='',F20.13)')TESTV2
              WRITE(LOGFIL,'(11X,''N_PART='',I6)')N_PART
     
              IF (ISCREN.GE.5) THEN
                  WRITE(0,'(''ITRNUM= '',I2)')ITRNUM
                  WRITE(0,'(''ISOSPI_AUXILI='',I2)')ISOSPI_AUXILI
                  WRITE(0,'(''In SPIN_J: TESTV2='',F20.13)')TESTV2
                  WRITE(0,'(11X,''N_PART='',I6)')N_PART
              END IF
     
C              STOP 'Stop in DENSIT_GRADNT: TESTV2 not correct!'
          END IF
      END IF
C
C=======================================================================     
C
c      IF (LOGWRI.GT.0) THEN
c          WRITE(LOGFIL,'(/,9X,''Exiting SPIN_J'')')
c      END IF
C
C=======================================================================
C
      CALL CPUTIM('SPIN_J',0)
C
C=======================================================================
C      
      RETURN
      END   
C      
C=======================================================================
C=======================================================================
C
      FUNCTION DENSIT_PLOTIN(INUCLI,ISOSPI_AUXILI,ZARGUM)
C      
      INCLUDE  'MATDIM/NDSPEC.f'
      INCLUDE  'MATDIM/NDIM_N.f'
      INCLUDE  'MATDIM/NDGAUS.f'
      INCLUDE  'MATDIM/NDNUCL.f'
C      
      PARAMETER
     *         (NDIM_L=NDIM_N,NDBASE=NDIM_N+1) 
      REAL*16
     *          QARGUM,QFLAGN,QDFLAG,QNORNL
      DIMENSION
     *          LABORD(0:NDIM_N,0:NDIM_L,0:2*NDIM_L+1)
      DIMENSION
     *          CMATUP(0:NDIM_L,1:NDBASE,1:NDBASE),
     *          CMATDN(0:NDIM_L,1:NDBASE,1:NDBASE)
      DIMENSION
     *          ENERUP(0:NDIM_L,1:NDBASE),
     *          ENERDN(0:NDIM_L,1:NDBASE)
      DIMENSION
     *          VCOEUP(0:NDIM_L,1:NDBASE),
     *          VCOEDN(0:NDIM_L,1:NDBASE)
      DIMENSION
     *          NWSSPH(1:NDSPEC),
     *          LWSSPH(1:NDSPEC),
     *          JWSSPH(1:NDSPEC)
      DIMENSION
     *          DENSIT_VECTOR(1:NDGAUS,0:NDIM_L),
     *          DENAUX_TEST_L(0:NDIM_L)
      DIMENSION
     *          QFLAGN(0:NDIM_N,0:NDIM_L),
     *          QDFLAG(0:NDIM_N,0:NDIM_L),
     *          XFLAGN(0:NDIM_N,0:NDIM_L),
     *          XDFLAG(0:NDIM_N,0:NDIM_L)
C
      COMMON
     *       /DATA_V/ AOSCIL_PROTON(1:NDNUCL),
     *                AOSCIL_NEUTRS(1:NDNUCL),
     *                HOMEGA(1:NDNUCL)   
      COMMON
     *       /GAUSHE/ NGAUSS,NSHELL_PROTON,NSHELL_NEUTRS
      COMMON
     *       /NMLCNT/ XNORNL(0:NDIM_N,0:NDIM_L)
      COMMON
     *       /POLDEN/ 
     *        XPOLYN_DENSIT(1:NDGAUS,0:NDIM_N,0:NDIM_L,0:NDIM_L), 
     *        XDPOLY_DENSIT(1:NDGAUS,0:NDIM_N,0:NDIM_L,0:NDIM_L)
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
     *       /DIMSPH/ LDSING_PROTON,LEVTHE_PROTON,
     *                LEVEXP_PROTON(1:NDNUCL),
     *                LDSING_NEUTRS,LEVTHE_NEUTRS,
     *                LEVEXP_NEUTRS(1:NDNUCL)
      COMMON
     *       /QUANTN/ NWSSPH_PROTON(1:NDSPEC),
     *      	      LWSSPH_PROTON(1:NDSPEC),
     *      	      JWSSPH_PROTON(1:NDSPEC),
     *
     *      	      NWSSPH_NEUTRS(1:NDSPEC),
     *      	      LWSSPH_NEUTRS(1:NDSPEC),
     *      	      JWSSPH_NEUTRS(1:NDSPEC)
       COMMON
     *       /LEVORD/ LABORD_PROTON(1:NDNUCL,
     *                              0:NDIM_N,0:NDIM_L,0:2*NDIM_L+1),
     *                LABORD_NEUTRS(1:NDNUCL,
     *                              0:NDIM_N,0:NDIM_L,0:2*NDIM_L+1)
      COMMON
     *       /NUCLEU/ IZ_FIX,IN_FIX,ISOSPI 
      COMMON
     *       /IFPAIR/ IF_PAI,IFGPAI,IFDPAI
      COMMON
     *       /BCSENE/ ENEMAX_BCSPAI(1:NDNUCL,1:2)
      COMMON
     *       /BCSDEL/ DELT_2(1:NDNUCL,1:2)
      COMMON
     *       /BCSLAM/ XLAMBD_PROTON(1:NDNUCL),
     *                XLAMBD_NEUTRS(1:NDNUCL)
C
C=======================================================================   
C  
C  ###  This FUNCTION calculates the TOTAL DENSITY FUNCTION Eq.(4.7.14)
C          
C  ###     It is prepared to distinct between netrons and protons and
C          the way of integrating:  Gauss-Laguerre  or  Gauss-Hermite
C
C=======================================================================     
C
C                        P R O T O N S 
C
C=======================================================================    
C      
      IF (ISOSPI_AUXILI.EQ.1) THEN
        
          N_PART=IZ_FIX
          NSHELL=NSHELL_PROTON
          LDSING=LDSING_PROTON
        
          AOSCIL=AOSCIL_PROTON(INUCLI)
          
          IF (IF_PAI.EQ.1) THEN
              XLAMBD=XLAMBD_PROTON(INUCLI)
              DELTA2=DELT_2(INUCLI,1)
              ENEMAX=ENEMAX_BCSPAI(INUCLI,1)
          END IF
        
          DO INDEXX=1,LDSING
          
             LWSSPH(INDEXX)=LWSSPH_PROTON(INDEXX)
             NWSSPH(INDEXX)=NWSSPH_PROTON(INDEXX)
             JWSSPH(INDEXX)=JWSSPH_PROTON(INDEXX)
          
          END DO
        
          DO INDEXX=1,LDSING
C
             JQNAUX=JWSSPH(INDEXX)
             LQNAUX=LWSSPH(INDEXX)
	     NNBAUX=NWSSPH(INDEXX)
             NBBASE=(NSHELL-LQNAUX)/2  
             
             LABORD(NNBAUX,LQNAUX,JQNAUX)
     *	    =LABORD_PROTON(INUCLI,NNBAUX,LQNAUX,JQNAUX)
     
             ENERUP(LQNAUX,NNBAUX+1)
     *      =ENERUP_AUXPRO(LQNAUX,NNBAUX+1)          

             ENERDN(LQNAUX,NNBAUX+1)
     *      =ENERDN_AUXPRO(LQNAUX,NNBAUX+1)
C
             DO N1NUMB=0,NBBASE
C		  
                CMATUP(LQNAUX,N1NUMB+1,NNBAUX+1)
     *         =CMATUP_AUXPRO(LQNAUX,N1NUMB+1,NNBAUX+1)
C		  
                CMATDN(LQNAUX,N1NUMB+1,NNBAUX+1)
     *         =CMATDN_AUXPRO(LQNAUX,N1NUMB+1,NNBAUX+1)
C
             END DO
         END DO                  
C
      END IF
C
C=======================================================================     
C
C                        N E U T R O N S 
C
C=======================================================================          
C      
      IF (ISOSPI_AUXILI.EQ.0) THEN
        
          N_PART=IN_FIX
          NSHELL=NSHELL_NEUTRS
          LDSING=LDSING_NEUTRS
        
          AOSCIL=AOSCIL_NEUTRS(INUCLI)
          
          IF (IF_PAI.EQ.1) THEN
              XLAMBD=XLAMBD_NEUTRS(INUCLI)
              DELTA2=DELT_2(INUCLI,2)
              ENEMAX=ENEMAX_BCSPAI(INUCLI,2)
          END IF
        
          DO INDEXX=1,LDSING
          
             LWSSPH(INDEXX)=LWSSPH_NEUTRS(INDEXX)
             NWSSPH(INDEXX)=NWSSPH_NEUTRS(INDEXX)
             JWSSPH(INDEXX)=JWSSPH_NEUTRS(INDEXX)
          
          END DO
        
          DO INDEXX=1,LDSING
C
             JQNAUX=JWSSPH(INDEXX)
             LQNAUX=LWSSPH(INDEXX)
	     NNBAUX=NWSSPH(INDEXX)
             NBBASE=(NSHELL-LQNAUX)/2  
             
             LABORD(NNBAUX,LQNAUX,JQNAUX)
     *	    =LABORD_NEUTRS(INUCLI,NNBAUX,LQNAUX,JQNAUX)
     
             ENERUP(LQNAUX,NNBAUX+1)
     *      =ENERUP_AUXNEU(LQNAUX,NNBAUX+1)          
C
             ENERDN(LQNAUX,NNBAUX+1)
     *      =ENERDN_AUXNEU(LQNAUX,NNBAUX+1)
C     
             DO N1NUMB=0,NBBASE
C		  
                CMATUP(LQNAUX,N1NUMB+1,NNBAUX+1)
     *         =CMATUP_AUXNEU(LQNAUX,N1NUMB+1,NNBAUX+1)
C		  
                CMATDN(LQNAUX,N1NUMB+1,NNBAUX+1)
     *         =CMATDN_AUXNEU(LQNAUX,N1NUMB+1,NNBAUX+1)
     
C
             END DO
         END DO                  
C
      END IF
C
C=======================================================================     
C
C       Calculating the polynomials in the desired ZARGUM point
C
C=======================================================================          
C      
      KIND_L=1
      NORDER=NSHELL/2
      LORDER=NSHELL
      
      QARGUM=QEXT(ZARGUM)    !passing the argument to real*16
      
      CALL LAGUER(QARGUM,QFLAGN,QDFLAG,NORDER,LORDER,KIND_L,
     *                                        NDIM_L,NDIM_N)
     
      DO NAUXIL=0,NORDER
        DO LAUXIL=0,LORDER
            XFLAGN(NAUXIL,LAUXIL)=REAL(QFLAGN(NAUXIL,LAUXIL))
            XDFLAG(NAUXIL,LAUXIL)=REAL(QDFLAG(NAUXIL,LAUXIL))
        END DO
      END DO
C
C=======================================================================     
C
C                   Calculating the DENSITY Eq.(4.7.14)
C
C=======================================================================        
C      
      PINUMB=4.0*ATAN(1.0)
      FACTOR=1.0/(PINUMB*(AOSCIL**3)) ! multiplied by 2/a 
C                                             from N_{n1,l}*N_{n2,l}
C
C  Here we calculate J_UP block
C
      DEN_UP=0.0
      
      DO LQPRIM=NSHELL,0,-1
         
         JNUMUP=2*LQPRIM+1  
          
         DO NWSNUM=1,(NSHELL-LQPRIM)/2+1
             
            ENERGY=ENERUP(LQPRIM,NWSNUM)
            
            IF (IF_PAI.EQ.1.AND.ENERGY.GT.ENEMAX) GO TO 1
           
            IF (IF_PAI.EQ.0
     *          .AND.LABORD(NWSNUM-1,LQPRIM,JNUMUP).EQ.0) GO TO 1
                
                DNL_UP=0.0
            
            DO N1NUMB=0,(NSHELL-LQPRIM)/2
                
               DO N2NUMB=0,(NSHELL-LQPRIM)/2
                  
                  CFACT1=(LQPRIM+1)*CMATUP(LQPRIM,N1NUMB+1,NWSNUM)
     *                             *CMATUP(LQPRIM,N2NUMB+1,NWSNUM)
    
                  POLPOL=XNORNL(N1NUMB,LQPRIM)*XFLAGN(N1NUMB,LQPRIM)
     *                  *XNORNL(N2NUMB,LQPRIM)*XFLAGN(N2NUMB,LQPRIM)
                
                  DNL_UP=DNL_UP+ZARGUM**LQPRIM*EXP(-ZARGUM)
     *                         *POLPOL*CFACT1
                  
               END DO !N2NUMB
              
            END DO !N1NUMB
            
            IF (IF_PAI.EQ.1) THEN
                V2_PUP=V2PAIR(ENERGY,XLAMBD,DELTA2)
                DEN_UP=DEN_UP+DNL_UP*V2_PUP  
            ELSE
                DEN_UP=DEN_UP+DNL_UP
            END IF
              
   1        CONTINUE
            
         END DO !NWSNUM
        
      END DO !LQPRIM
C
C  Here we calculate J_DOWN block
C
      DEN_DW=0.0
      
      DO LQPRIM=NSHELL,0,-1 
          
         JNUMDW=2*LQPRIM-1
          
         DO NWSNUM=1,(NSHELL-LQPRIM)/2+1
             
            ENERGY=ENERDN(LQPRIM,NWSNUM)
            
            IF (JNUMDW.GE.0) THEN
                
                IF (IF_PAI.EQ.1.AND.ENERGY.GT.ENEMAX) GO TO 2
           
                IF (IF_PAI.EQ.0
     *              .AND.LABORD(NWSNUM-1,LQPRIM,JNUMDW).EQ.0) GO TO 2
                    
                    DNL_DW=0.0
            
                    DO N1NUMB=0,(NSHELL-LQPRIM)/2
                
                       DO N2NUMB=0,(NSHELL-LQPRIM)/2
    
                          CFACT2=(LQPRIM)
     *                          *CMATDN(LQPRIM,N1NUMB+1,NWSNUM)
     *                          *CMATDN(LQPRIM,N2NUMB+1,NWSNUM)
    
                          POLPOL
     *                   =
     *                    XNORNL(N1NUMB,LQPRIM)*XFLAGN(N1NUMB,LQPRIM)
     *                   *XNORNL(N2NUMB,LQPRIM)*XFLAGN(N2NUMB,LQPRIM)
                
                          DNL_DW=DNL_DW+ZARGUM**LQPRIM*EXP(-ZARGUM)
     *                                 *POLPOL*CFACT2
                  
                       END DO !N2NUMB
              
                    END DO !N1NUMB
                    
                    IF (IF_PAI.EQ.1) THEN
                        V2_PDW=V2PAIR(ENERGY,XLAMBD,DELTA2)
                        DEN_DW=DEN_DW+DNL_DW*V2_PDW
                   ELSE
                        DEN_DW=DEN_DW+DNL_DW
                   END IF

            END IF
            
   2        CONTINUE 
            
         END DO !NWSNUM
        
      END DO !LQPRIM
C
C=======================================================================
C      
      DENSIT_PLOTIN=(DEN_UP+DEN_DW)*FACTOR
C
C=======================================================================
C     
      RETURN
      END
C      
C=======================================================================
C=======================================================================
C
      SUBROUTINE RHONLJ_SUMMED(INUCLI,ISOSPI,ZARGUM,IPOINT,
     *                                N_CURV,Y_CURV,L_CURV)
C    
      INCLUDE 'MATDIM/NDMESH.f'
      INCLUDE 'MATDIM/NDSPEC.f'
      INCLUDE 'MATDIM/NDMAIN.f'
      INCLUDE 'MATDIM/NDNUCL.f'
      INCLUDE 'MATDIM/NDIM_N.f'
      
      PARAMETER
     *         (NDORBI=NDMAIN,NDJTOT=2*NDORBI+1,NDIM_L=NDIM_N,
     *                                          NDBASE=NDIM_N+1) 
      CHARACTER
     *          L_CURV*100,LBSHEL*6,JOPTIO*2,AUXLAB*100
      CHARACTER 
     *          LABTEX*11,LABEXP*6
      DIMENSION
     *          LBSHEL(0:NDMAIN,0:NDORBI,1:NDJTOT)
      DIMENSION
     *          LABORD(0:NDIM_N,0:NDIM_L,0:2*NDIM_L+1)
      DIMENSION
     *          NWSSPH(1:NDSPEC),
     *          LWSSPH(1:NDSPEC),
     *          JWSSPH(1:NDSPEC)
      DIMENSION
     *          ENERUP(0:NDIM_L,1:NDBASE),
     *          ENERDN(0:NDIM_L,1:NDBASE)
      DIMENSION
     *          Y_CURV(1:NDMESH,1:NDSPEC),
     *          L_CURV(1:NDMESH,1:NDSPEC)
      COMMON
     *       /GAUSHE/ NGAUSS,NSHELL_PROTON,NSHELL_NEUTRS
      COMMON
     *       /EUDAUX/ ENERUP_AUXPRO(0:NDIM_L,1:NDBASE),
     *                ENERDN_AUXPRO(0:NDIM_L,1:NDBASE),
     *
     *                ENERUP_AUXNEU(0:NDIM_L,1:NDBASE),
     *                ENERDN_AUXNEU(0:NDIM_L,1:NDBASE)
      COMMON
     *       /DIMSPH/ LDSING_PROTON,LEVTHE_PROTON,
     *                LEVEXP_PROTON(1:NDNUCL),
     *                LDSING_NEUTRS,LEVTHE_NEUTRS,
     *                LEVEXP_NEUTRS(1:NDNUCL)
      COMMON
     *       /QUANTN/ NWSSPH_PROTON(1:NDSPEC),
     *      	      LWSSPH_PROTON(1:NDSPEC),
     *      	      JWSSPH_PROTON(1:NDSPEC),
     *
     *      	      NWSSPH_NEUTRS(1:NDSPEC),
     *      	      LWSSPH_NEUTRS(1:NDSPEC),
     *      	      JWSSPH_NEUTRS(1:NDSPEC)
      COMMON
     *       /LEVORD/ LABORD_PROTON(1:NDNUCL,
     *                              0:NDIM_N,0:NDIM_L,0:2*NDIM_L+1),
     *                LABORD_NEUTRS(1:NDNUCL,
     *                              0:NDIM_N,0:NDIM_L,0:2*NDIM_L+1)
      COMMON
     *       /IFPAIR/ IF_PAI,IFGPAI,IFDPAI
      COMMON
     *       /BCSENE/ ENEMAX_BCSPAI(1:NDNUCL,1:2)
      COMMON
     *       /BCSDEL/ DELT_2(1:NDNUCL,1:2)
      COMMON
     *       /BCSLAM/ XLAMBD_PROTON(1:NDNUCL),
     *                XLAMBD_NEUTRS(1:NDNUCL)
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
C     This subroutine sums "state by state" the single 
C     state densities from the DENSIT_STATES subroutine.
C
C     In this way we will be able to plot the "partial" densities in 
C     such a way that we will have as many curves as occupied states
C     and the first curve is the first state curve, and the last one
C     is the total sum (i.e. total density) of the state densities.
C=======================================================================     
C
C                        P R O T O N S 
C
C=======================================================================    
C      
      IF (ISOSPI.EQ.1) THEN
C        
        NSHELL=NSHELL_PROTON
        LDSING=LDSING_PROTON
        
        IF (IF_PAI.EQ.1) THEN
            XLAMBD=XLAMBD_PROTON(INUCLI)
            DELTA2=DELT_2(INUCLI,1)
            ENEMAX=ENEMAX_BCSPAI(INUCLI,1)
        END IF
C        
        DO INDEXX=1,LDSING
C          
          LWSSPH(INDEXX)=LWSSPH_PROTON(INDEXX)
          NWSSPH(INDEXX)=NWSSPH_PROTON(INDEXX)
          JWSSPH(INDEXX)=JWSSPH_PROTON(INDEXX)
          
        END DO
C        
        DO INDEXX=1,LDSING
C
           JQNAUX=JWSSPH(INDEXX)
           LQNAUX=LWSSPH(INDEXX)
	   NNBAUX=NWSSPH(INDEXX)
           NBBASE=(NSHELL-LQNAUX)/2  
C             
           LABORD(NNBAUX,LQNAUX,JQNAUX)
     *	  =
     *     LABORD_PROTON(INUCLI,NNBAUX,LQNAUX,JQNAUX)
     
           ENERUP(LQNAUX,NNBAUX+1)
     *    =ENERUP_AUXPRO(LQNAUX,NNBAUX+1)          
C
           ENERDN(LQNAUX,NNBAUX+1)
     *    =ENERDN_AUXPRO(LQNAUX,NNBAUX+1)
C
        END DO                  
C
      END IF
C
C=======================================================================     
C
C                        N E U T R O N S 
C
C=======================================================================          
C      
      IF (ISOSPI.EQ.0) THEN
C        
          NSHELL=NSHELL_NEUTRS
          LDSING=LDSING_NEUTRS
        
          IF (IF_PAI.EQ.1) THEN
            XLAMBD=XLAMBD_NEUTRS(INUCLI)
            DELTA2=DELT_2(INUCLI,2)
            ENEMAX=ENEMAX_BCSPAI(INUCLI,2)
          END IF
C        
          DO INDEXX=1,LDSING
C          
             LWSSPH(INDEXX)=LWSSPH_NEUTRS(INDEXX)
             NWSSPH(INDEXX)=NWSSPH_NEUTRS(INDEXX)
             JWSSPH(INDEXX)=JWSSPH_NEUTRS(INDEXX)
C          
          END DO
C        
        DO INDEXX=1,LDSING
C
           JQNAUX=JWSSPH(INDEXX)
           LQNAUX=LWSSPH(INDEXX)
	   NNBAUX=NWSSPH(INDEXX)
           NBBASE=(NSHELL-LQNAUX)/2  
C             
           LABORD(NNBAUX,LQNAUX,JQNAUX)
     *	  =
     *     LABORD_NEUTRS(INUCLI,NNBAUX,LQNAUX,JQNAUX)
     
           ENERUP(LQNAUX,NNBAUX+1)
     *    =ENERUP_AUXNEU(LQNAUX,NNBAUX+1)          
C
           ENERDN(LQNAUX,NNBAUX+1)
     *    =ENERDN_AUXNEU(LQNAUX,NNBAUX+1)
C     
        END DO                  
C
      END IF
C
C=======================================================================          
C
      RHOSUM=0.0
C      
      I_CURV=0
C      
      DO INDEXX=1,LDSING
C
         JSMALL=0
C            
         LQNUMB=LWSSPH(INDEXX)
         NNUMBR=NWSSPH(INDEXX)
         JQNUMB=JWSSPH(INDEXX)
            
         NWSNUM=NNUMBR+1
            
         J1AUXI=2*LQNUMB+1
            
         IF (J1AUXI.EQ.JQNUMB) THEN
             JOPTIO='UP'
             J2AUXI=JQNUMB
             ENERGY=ENERUP(LQNUMB,NNUMBR+1)
         END IF
            
         J1AUXI=2*LQNUMB-1
            
         IF (J1AUXI.EQ.JQNUMB) THEN
             JOPTIO='DW'
             J2AUXI=JQNUMB
             ENERGY=ENERDN(LQNUMB,NNUMBR+1)
         END IF
C            
         IF (IF_PAI.EQ.1.AND.ENERGY.GT.ENEMAX) GO TO 1
         
         IF (IF_PAI.EQ.0.AND.LABORD(NWSNUM-1,LQNUMB,JQNUMB).EQ.0)
     *                                                         GO TO 1
C            
         CALL RHONLJ_STATES(INUCLI,ISOSPI,LQNUMB,NWSNUM,
     *                             JOPTIO,ZARGUM,RHONLJ)
         I_CURV=I_CURV+1
                
         NSMALL=NNUMBR
         LSMALL=LQNUMB
         JSMALL=JQNUMB
         N_BIGG=2*NSMALL+LSMALL
C            
         IF(IF_PAI.EQ.1) THEN
            V2FACT=V2PAIR(ENERGY,XLAMBD,DELTA2)
            RHOSUM=RHOSUM+RHONLJ!*V2FACT
         ELSE
            RHOSUM=RHOSUM+RHONLJ
         END IF
C            
         Y_CURV(IPOINT,I_CURV)=RHOSUM
C                
         LABEXP=LBSHEL(N_BIGG,LSMALL,JSMALL)
         CALL LATEXS_LABELS(LABEXP,LABTEX)
C                
         WRITE(AUXLAB,'(''\\boldmath'',a)') LABTEX
C     
         L_CURV(IPOINT,I_CURV)=AUXLAB
C
   1     CONTINUE
C                     
      END DO
C      
      N_CURV=I_CURV
C
C=======================================================================
C      
      RETURN
      END
C      
C=======================================================================
C=======================================================================
C
      SUBROUTINE RHONLJ_STATES(INUCLI,ISOSPI,LQNUMB,NWSNUM,
     *                                JOPTIO,ZARGUM,RHONLJ)
     
C
      INCLUDE 'MATDIM/NDGAUS.f'
      INCLUDE 'MATDIM/NDSPEC.f'
      INCLUDE 'MATDIM/NDIM_N.f'
      INCLUDE 'MATDIM/NDNUCL.f'
      
      PARAMETER
     *         (NDIM_L=NDIM_N,NDBASE=NDIM_N+1) 
      REAL*16
     *          QARGUM,QFLAGN,QDFLAG,QNORNL
      CHARACTER
     *          JOPTIO*2
      DIMENSION
     *          LABORD(0:NDIM_N,0:NDIM_L,0:2*NDIM_L+1)
      DIMENSION
     *          QFLAGN(0:NDIM_N,0:NDIM_L),
     *          QDFLAG(0:NDIM_N,0:NDIM_L),
     *          XFLAGN(0:NDIM_N,0:NDIM_L),
     *          XDFLAG(0:NDIM_N,0:NDIM_L)
      DIMENSION
     *          CMATUP(0:NDIM_L,1:NDBASE,1:NDBASE),
     *          CMATDN(0:NDIM_L,1:NDBASE,1:NDBASE)
      DIMENSION
     *          VCOEUP(0:NDIM_L,1:NDBASE),
     *          VCOEDN(0:NDIM_L,1:NDBASE)
      DIMENSION
     *          ENERUP(0:NDIM_L,1:NDBASE),
     *          ENERDN(0:NDIM_L,1:NDBASE)
      DIMENSION
     *          NWSSPH(1:NDSPEC),
     *          LWSSPH(1:NDSPEC),
     *          JWSSPH(1:NDSPEC)
      DIMENSION
     *          DENJUP(1:NDGAUS,0:NDIM_L),
     *          DENJDW(1:NDGAUS,0:NDIM_L)
      COMMON
     *       /GAUSHE/ NGAUSS,NSHELL_PROTON,NSHELL_NEUTRS
      COMMON
     *       /NMLCNT/ XNORNL(0:NDIM_N,0:NDIM_L)
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
     *       /DIMSPH/ LDSING_PROTON,LEVTHE_PROTON,
     *                LEVEXP_PROTON(1:NDNUCL),
     *                LDSING_NEUTRS,LEVTHE_NEUTRS,
     *                LEVEXP_NEUTRS(1:NDNUCL)
      COMMON
     *       /QUANTN/ NWSSPH_PROTON(1:NDSPEC),
     *      	      LWSSPH_PROTON(1:NDSPEC),
     *      	      JWSSPH_PROTON(1:NDSPEC),
     *
     *      	      NWSSPH_NEUTRS(1:NDSPEC),
     *      	      LWSSPH_NEUTRS(1:NDSPEC),
     *      	      JWSSPH_NEUTRS(1:NDSPEC)
      COMMON
     *       /DATA_V/ AOSCIL_PROTON(1:NDNUCL),
     *                AOSCIL_NEUTRS(1:NDNUCL),
     *                HOMEGA(1:NDNUCL)  
      COMMON
     *       /LEVORD/ LABORD_PROTON(1:NDNUCL,
     *                              0:NDIM_N,0:NDIM_L,0:2*NDIM_L+1),
     *                LABORD_NEUTRS(1:NDNUCL,
     *                              0:NDIM_N,0:NDIM_L,0:2*NDIM_L+1)
C
C=======================================================================     
C   This subroutine calculates the density \rho(r) of every
C   occupied state  s e p a r a t e l y (for fixed NLJ numbers) and
c   the sum runs only for the indices n_1 and n_2
C=======================================================================     
C
C                        P R O T O N S 
C
C=======================================================================    
C      
      IF (ISOSPI.EQ.1) THEN
C        
        NSHELL=NSHELL_PROTON
        LDSING=LDSING_PROTON
C        
        AOSCIL=AOSCIL_PROTON(INUCLI)
C        
        DO INDEXX=1,LDSING
C          
          LWSSPH(INDEXX)=LWSSPH_PROTON(INDEXX)
          NWSSPH(INDEXX)=NWSSPH_PROTON(INDEXX)
          JWSSPH(INDEXX)=JWSSPH_PROTON(INDEXX)
          
        END DO
C        
        DO INDEXX=1,LDSING
C
           JQNAUX=JWSSPH(INDEXX)
           LQNAUX=LWSSPH(INDEXX)
	   NNBAUX=NWSSPH(INDEXX)
           NBBASE=(NSHELL-LQNAUX)/2  
C             
           LABORD(NNBAUX,LQNAUX,JQNAUX)
     *	  =LABORD_PROTON(INUCLI,NNBAUX,LQNAUX,JQNAUX)
C     
           ENERUP(LQNAUX,NNBAUX+1)
     *    =ENERUP_AUXPRO(LQNAUX,NNBAUX+1)          
C
           ENERDN(LQNAUX,NNBAUX+1)
     *    =ENERDN_AUXPRO(LQNAUX,NNBAUX+1)
C
           DO N1NUMB=0,NBBASE
C		  
              CMATUP(LQNAUX,N1NUMB+1,NNBAUX+1)
     *       =CMATUP_AUXPRO(LQNAUX,N1NUMB+1,NNBAUX+1)
C		  
              CMATDN(LQNAUX,N1NUMB+1,NNBAUX+1)
     *       =CMATDN_AUXPRO(LQNAUX,N1NUMB+1,NNBAUX+1)
C
           END DO
        END DO                  
C
      END IF
C
C=======================================================================     
C
C                        N E U T R O N S 
C
C=======================================================================          
C      
      IF (ISOSPI.EQ.0) THEN
C        
          NSHELL=NSHELL_NEUTRS
          LDSING=LDSING_NEUTRS
C        
          AOSCIL=AOSCIL_NEUTRS(INUCLI)
C        
          DO INDEXX=1,LDSING
C          
             LWSSPH(INDEXX)=LWSSPH_NEUTRS(INDEXX)
             NWSSPH(INDEXX)=NWSSPH_NEUTRS(INDEXX)
             JWSSPH(INDEXX)=JWSSPH_NEUTRS(INDEXX)
C          
          END DO
C        
          DO INDEXX=1,LDSING
C
             JQNAUX=JWSSPH(INDEXX)
             LQNAUX=LWSSPH(INDEXX)
	     NNBAUX=NWSSPH(INDEXX)
             NBBASE=(NSHELL-LQNAUX)/2  
C             
             LABORD(NNBAUX,LQNAUX,JQNAUX)
     *	    =LABORD_NEUTRS(INUCLI,NNBAUX,LQNAUX,JQNAUX)
C     
             ENERUP(LQNAUX,NNBAUX+1)
     *      =ENERUP_AUXNEU(LQNAUX,NNBAUX+1)          
C
             ENERDN(LQNAUX,NNBAUX+1)
     *      =ENERDN_AUXNEU(LQNAUX,NNBAUX+1)
C
             DO N1NUMB=0,NBBASE
C		  
                CMATUP(LQNAUX,N1NUMB+1,NNBAUX+1)
     *         =CMATUP_AUXNEU(LQNAUX,N1NUMB+1,NNBAUX+1)
C		  
                CMATDN(LQNAUX,N1NUMB+1,NNBAUX+1)
     *         =CMATDN_AUXNEU(LQNAUX,N1NUMB+1,NNBAUX+1)
C
             END DO
          END DO                  
C
      END IF
C
C=======================================================================     
C
C       Calculating the polynomials in the desired ZARGUM point
C
C=======================================================================          
C      
      KIND_L=1
      NORDER=NSHELL/2
      LORDER=NSHELL
      
      QARGUM=QEXT(ZARGUM)    !passing the argument to real*16
      
      CALL LAGUER(QARGUM,QFLAGN,QDFLAG,NORDER,LORDER,KIND_L,
     *                                        NDIM_L,NDIM_N)
     
      DO NAUXIL=0,NORDER
        DO LAUXIL=0,LORDER
            XFLAGN(NAUXIL,LAUXIL)=REAL(QFLAGN(NAUXIL,LAUXIL))
            XDFLAG(NAUXIL,LAUXIL)=REAL(QDFLAG(NAUXIL,LAUXIL))
        END DO
      END DO
C
C=======================================================================     
C
      PINUMB=4.*ATAN(1.0)
      FACTOR=1./(PINUMB*AOSCIL**3)
      
      DENAUX=0.0
            
      DO N1NUMB=0,(NSHELL-LQNUMB)/2
                
         DO N2NUMB=0,(NSHELL-LQNUMB)/2
                
            CFACT1=(LQNUMB+1)*CMATUP(LQNUMB,N1NUMB+1,NWSNUM)
     *                       *CMATUP(LQNUMB,N2NUMB+1,NWSNUM)
    
            CFACT2=(LQNUMB  )*CMATDN(LQNUMB,N1NUMB+1,NWSNUM)
     *                       *CMATDN(LQNUMB,N2NUMB+1,NWSNUM)
     
            POLPRD=ZARGUM**LQNUMB*EXP(-ZARGUM)
     *            *XNORNL(N1NUMB,LQNUMB)*XFLAGN(N1NUMB,LQNUMB)
     *            *XNORNL(N2NUMB,LQNUMB)*XFLAGN(N2NUMB,LQNUMB)
              
            IF (JOPTIO.EQ.'UP') THEN
               
                DENAUX=DENAUX+CFACT1*POLPRD
               
            END IF
               
            IF (JOPTIO.EQ.'DW') THEN
               
                DENAUX=DENAUX+CFACT2*POLPRD

            END IF
                
         END DO
          
      END DO
      
      RHONLJ=DENAUX*FACTOR
C
C=======================================================================
C      
      RETURN
      END
C      
C=======================================================================
C=======================================================================
C      
      FUNCTION V2PAIR(ENERGY,XLAMBD,DELTA2)
      
      COMMON
     *       /PRFILE/ ISCREN,LOGFIL,LOGCHP,LOGCHN,LOGHIP,LOGHIN,
     *                LOGPRO,LOGNEU,LOGWRI,LSCREN,IRESUL,ICONVE,
     *                NOUTPT,LOGAUX,LOGENE,LOGRAD,LOGGRD,LOGMSH,
     *                                                   LOGBIS  
C
C=======================================================================
C 
C      WRITE(LOGFIL,'(''in FUNCTION V2PAIR: XLAMBD= '',F20.13,
C     *               '' DELTA2= '',F20.13,'' ENERGY= '',F20.13)')
C     *                  XLAMBD,DELTA2,ENERGY
C
C=======================================================================
C         
      V2PAIR=0.5*(1-((ENERGY-XLAMBD)/SQRT((ENERGY-XLAMBD)**2+DELTA2)))
C
C=======================================================================
C      
      RETURN
      END
C      
C=======================================================================
C=======================================================================
C
      SUBROUTINE CONSIS_VERIFS(IABORT)
C
      INCLUDE   'MATDIM/NDNUCL.f'
      INCLUDE   'MATDIM/NDSPEC.f'
      INCLUDE   'MATDIM/NITERA.f'
C
      COMMON
     *       /DIMSPH/ LDSING_PROTON,LEVTHE_PROTON,
     *                              LEVEXP_PROTON(1:NDNUCL),
     *                LDSING_NEUTRS,LEVTHE_NEUTRS,
     *                              LEVEXP_NEUTRS(1:NDNUCL)
      COMMON
     *       /ENRGYY/ ENERGY_PROTON(1:NDSPEC,0:NITERA), 
     *                ENERGY_NEUTRS(1:NDSPEC,0:NITERA)
      COMMON
     *       /QUANTN/ NWSSPH_PROTON(1:NDSPEC),
     *                LWSSPH_PROTON(1:NDSPEC),
     *                JWSSPH_PROTON(1:NDSPEC),
     *
     *                NWSSPH_NEUTRS(1:NDSPEC),
     *                LWSSPH_NEUTRS(1:NDSPEC),
     *                JWSSPH_NEUTRS(1:NDSPEC)
      COMMON
     *       /TESFIT/ IFTEST,IFFITS,IFMESH,ISPE_P,IFMOCA
      COMMON
     *       /DENSIF/ ITRNUM
     *       /WHENST/ DIFFVA
      COMMON
     *       /PRFILE/ ISCREN,LOGFIL,LOGCHP,LOGCHN,LOGHIP,LOGHIN,
     *                LOGPRO,LOGNEU,LOGWRI,LSCREN,IRESUL,ICONVE,
     *                NOUTPT,LOGAUX,LOGENE,LOGRAD,LOGGRD,LOGMSH,
     *                                                   LOGBIS   
C      
C=======================================================================
C     This subroutine verifies the self-consistency condition in terms
C     of single-particle eigenvalues for the density-dependent Hamilt.     
C=======================================================================
C
      IABORT=0
C                        Starting with protons
      DIFFVA=0.0
C      
      IF (IFTEST.EQ.1) THEN
      
          IF (LOGWRI.GT.1) THEN
      
              WRITE(LOGFIL,'(/,''Testing the Self-Consistency'')')
      
              WRITE(LOGFIL,'(/,''ITRNUM    ENERGY(ITRUM-1)     '',
     *                                         ''ENERGY(ITRNUM)'',/)')
      
              WRITE(LOGFIL,'(''Protons'')')
       
          END IF
      
      END IF
C      
C=======================================================================
C
      DO I=1,LDSING_PROTON
C
         IF (ITRNUM.EQ.0) THEN
             DIFFVA=DIFFVA+ENERGY_PROTON(I,ITRNUM)**2
             GO TO 1
         END IF
C	    
         DO J=1,LDSING_PROTON
C	    
            IF (NWSSPH_PROTON(I).EQ.NWSSPH_PROTON(J).AND.
     *	        LWSSPH_PROTON(I).EQ.LWSSPH_PROTON(J).AND.
     *	        JWSSPH_PROTON(I).EQ.JWSSPH_PROTON(J)) THEN
C     
                DIFFVA=DIFFVA+(ENERGY_PROTON(I,ITRNUM)
     *                        -ENERGY_PROTON(J,ITRNUM-1))**2
                
                IF (IFTEST.EQ.1) THEN
                
                IF (LOGWRI.GT.1)  WRITE(LOGFIL,'(I5.3,3ES20.10)')ITRNUM,
     *                                       ENERGY_PROTON(I,ITRNUM-1),
     *                                       ENERGY_PROTON(J,ITRNUM),
     *                                       ENERGY_PROTON(I,ITRNUM-1)-
     *                                       ENERGY_PROTON(J,ITRNUM)
                END IF
C     
	        GO TO 1
C     
            END IF
C	    
         END DO
C	    
   1     CONTINUE
C	    
      END DO	         
C      
C=======================================================================      
C                      ... continuing with neutrons
C=======================================================================      
C     
  
      IF (IFTEST.EQ.1) THEN
                
      IF (LOGWRI.GT.1) WRITE(21,'(''Neutrons'')')
      
      END IF
      
      DO I=1,LDSING_NEUTRS
C
         IF (ITRNUM.EQ.0) THEN
	     DIFFVA=DIFFVA+ENERGY_NEUTRS(I,ITRNUM)**2
	     GO TO 2
	 END IF
C	    
         DO J=1,LDSING_NEUTRS
C	    
            IF (NWSSPH_NEUTRS(I).EQ.NWSSPH_NEUTRS(J).AND.
     *	        LWSSPH_NEUTRS(I).EQ.LWSSPH_NEUTRS(J).AND.
     *	        JWSSPH_NEUTRS(I).EQ.JWSSPH_NEUTRS(J)) THEN
C     
                DIFFVA=DIFFVA+(ENERGY_NEUTRS(I,ITRNUM)
     *                        -ENERGY_NEUTRS(J,ITRNUM-1))**2
                
                IF (IFTEST.EQ.1) THEN
                
                IF (LOGWRI.GT.1) WRITE(21,'(I5.3,2ES20.10)')ITRNUM,
     *                                   ENERGY_NEUTRS(I,ITRNUM-1),
     *                                     ENERGY_NEUTRS(J,ITRNUM)
     
                END IF
C     
	        GO TO 2
C     
            END IF
C	    
         END DO
C	    
   2     CONTINUE
C	    
      END DO	 
C           
C=======================================================================
C      
      IF (ITRNUM.EQ.NITERA) THEN
          WRITE(0,'(''ITRNUM>'',I3,'' No convergence in HAMMAT'')')
     *                NITERA
	  IABORT=1
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
      FUNCTION CHOICE(ISOSPI,IDEFCN,INUCLI)
    
      INCLUDE   'MATDIM/NDNUCL.f'
      INCLUDE   'MATDIM/NDSPEC.f'
      INCLUDE   'MATDIM/NDLEXP.f'
      INCLUDE   'MATDIM/NDIM_M.f'
      INCLUDE   'MATDIM/ND_RHO.f'
      INCLUDE   'MATDIM/N_NOYX.f'
      
      CHARACTER
     *          TYPCHI*6,SYMBNU*5
      CHARACTER
     *          LABTHE_PROTON*6,LABTHE_NEUTRS*6,
     *          LABEXP_PROTON*6,LABEXP_NEUTRS*6
      CHARACTER
     *          LABEXP*6,LABTHE*6
C      DIMENSION
C     *          MEXPER(1:NDIM_M),MTHEOR(1:NDIM_M)
      DIMENSION
     *          ENETHE(1:NDSPEC),LABTHE(1:NDSPEC)
      DIMENSION
     *          EXPEXP(1:NDLEXP),IDEGEX(1:NDLEXP)
      DIMENSION
     *          IDGLEV(1:NDSPEC)
      DIMENSION
     *          LABEXP(1:NDLEXP)
      DIMENSION
     *          DEGENE(1:NDIM_M)
      DIMENSION
     *          GAPEXP(1:NDNUCL),FERMEX(1:NDNUCL),
     *          DENSUP(1:NDNUCL),DENSDW(1:NDNUCL)
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
     *       /MODIFR/ RMSEXP_PROTON(1:NDNUCL),
     *                RMSEXP_NEUTRS(1:NDNUCL),
     *                RMSTHE_PROTON,RMSTHE_NEUTRS,
     *                RMODIF_PROTON,RMODIF_NEUTRS,
     *                              PUSHIN,DELTAR
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
     *       /CHIVAL/ CHISQU_CORREL,DIFSQU_RADIUS,CHISQU_INVERT,
     *                DIFSQU_EFERMI,DIFSQU_ENEGAP,CHIWEI_ENERGY,
     *                CHIWEI_ENEDEG,ERRABS_WEIDEG,EABSAV,ERRMAX,
     *                       DIFSQU_DENSUP,DIFSQU_DENSDW,CHIRHO
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
     *       /DECOMP/ VEICOR,VEIRAD,VEIINV,VEIFER,VEIGAP,VEIWEI,
     *                       VEIDIF,VEIABS,VEIMAX,VEIDUP,VEIDDW
      COMMON 
     *       /CHITAK/ IF_SPE,IF_RAD,IF_GAP,IF_FER,IF_DEN,IF_RHO,IF_INV,
     *                LDCHI2,TYPCHI(1:NDNUCL),ITKCHI(1:NDNUCL) 
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
     *       /CNTCNT/ ICOUNT_DENGRD,ICOUNT_FUNMIN,ICOUNT_HAMMAT,
     *                ICOUNT_EXPTHE,ICOUNT_CHOICE,ICOUNT_DIAMAT
C
C=======================================================================
C
C     This SUBROUTINE collects all the contributions to the \chi^2.
C
C     It is called from EXPTHE and it is mainly used for printing.
C
C=======================================================================
C
      ICOUNT_CHOICE=ICOUNT_CHOICE+1

      IF (ISOSPI.EQ.1) THEN
C
          LEVTHE=LDSING_PROTON ! Number of levels below user defined 
C                                                       limit SPHMAX
          LEVEXP=LEVEXP_PROTON(INUCLI)
          RMSEXP=RMSEXP_PROTON(INUCLI)
          RMSTHE=RMSTHE_PROTON
C
          DO ITHEOR=1,LEVTHE
             LABTHE(ITHEOR)=LABTHE_PROTON(ITHEOR)
             ENETHE(ITHEOR)=ENETHE_PROTON(ITHEOR)
             IDGLEV(ITHEOR)=IDGLEV_PROTON(ITHEOR)
          END DO
C          
          DO IEXPER=1,LEVEXP
             LABEXP(IEXPER)=LABEXP_PROTON(INUCLI,IEXPER)
             EXPEXP(IEXPER)=EXPEXP_PROTON(INUCLI,IEXPER)
             IDEGEX(IEXPER)=IDEGEX_PROTON(INUCLI,IEXPER)
          END DO
          
          WEINUC=WEINUC_PROTON(INUCLI)
          WEIRAD=WEIRAD_PROTON(INUCLI)
C          
      END IF 
C
      IF (ISOSPI.EQ.0) THEN
C
          LEVTHE=LDSING_NEUTRS  
          LEVEXP=LEVEXP_NEUTRS(INUCLI)
C
          RMSEXP=RMSEXP_NEUTRS(INUCLI)
          RMSTHE=RMSTHE_NEUTRS
C          
          DO ITHEOR=1,LEVTHE
             LABTHE(ITHEOR)=LABTHE_NEUTRS(ITHEOR)
             ENETHE(ITHEOR)=ENETHE_NEUTRS(ITHEOR)
             IDGLEV(ITHEOR)=IDGLEV_NEUTRS(ITHEOR)
          END DO
C          
          DO IEXPER=1,LEVEXP
             LABEXP(IEXPER)=LABEXP_NEUTRS(INUCLI,IEXPER)
             EXPEXP(IEXPER)=EXPEXP_NEUTRS(INUCLI,IEXPER)
             IDEGEX(IEXPER)=IDEGEX_NEUTRS(INUCLI,IEXPER)
          END DO
C         
          WEINUC=WEINUC_NEUTRS(INUCLI)
          WEIRAD=WEIRAD_NEUTRS(INUCLI)
  
      END IF
C
C=======================================================================
C
      VEICOR=0.0
      VEIRAD=0.0
      VEIINV=0.0
      VEIFER=0.0
      VEIGAP=0.0
      VEIWEI=0.0
      VEIDIF=0.0
      VEIABS=0.0
      VEIMAX=0.0
      VEIDUP=0.0
      VEIDDW=0.0
C      
      NO_ACT=1
C
C=======================================================================
C     Initialising all the chi^2 values 
C=======================================================================
C      
      CHI2_1=0.  ! Single Particle Energies
      CHI2_2=0.  ! RMS
      CHI2_3=0.  ! Gap Energy
      CHI2_4=0.  ! Fermi Energy
      CHI2_5=0.  ! Densities below and above the shell closure
      CHI2_6=0.  ! Density rho(r)
      CHI2_7=0.  ! Level inversions      
C
C=======================================================================
C     Starting with CHI2_1, related to SPE
C=======================================================================
C        
       IF (IF_SPE.EQ.1) THEN
C       
C       WRITE(LOGFIL,'(''Calculating chi2_1...'')')
C           
         CHI2_1=CHI2_1+CHIWEI_ENEDEG*WEINUC
C         
         VEIWEI=1
C         
         NO_ACT=0
C       
       END IF !IF_SPE=1
C       
C       WRITE(LOGFIL,'(''Calculating chi2_1...OK!'')')
C
C=======================================================================
C               Continuating with CHI2_2, realted to RMS
C=======================================================================
C      
      IF (IF_RAD.EQ.1) THEN 
C       
       CHI2_2=CHI2_2+DIFSQU_RADIUS*WEIRAD*WEINUC
C       
       VEIRAD=1
C       
       NO_ACT=0
C      
      END IF ! IF_RAD=1
C
C=======================================================================
C                  CHI2_3, realted to the Gap Energy
C=======================================================================
C      
      IF (IF_GAP.EQ.1) THEN 
C       
       CHI2_3=CHI2_3+DIFSQU_ENEGAP*WEIGHT_ENEGAP
C       
       VEIGAP=1
C       
       NO_ACT=0
C      
      END IF ! IF_GAP=1  
C
C=======================================================================
C                  CHI2_4, realted to the Fermi Energy
C                   (first protons, then neutrons)
C=======================================================================
C      
      IF (IF_FER.EQ.1) THEN 
C       
       CHI2_4=CHI2_4+DIFSQU_EFERMI*WEIGHT_EFERMI
C       
       VEIFER=1
C       
       NO_ACT=0
C      
      END IF ! IF_FER=1  
C
C=======================================================================
C   CHI2_5, realted to the Densities below and above the shell closure
C                   (first protons, then neutrons)
C=======================================================================
C      
      IF (IF_DEN.EQ.1) THEN 
C       
          CHI2_5=CHI2_5+DIFSQU_DENSUP*WEIGHT_DENSDW
     *                 +DIFSQU_DENSDW*WEIGHT_DENSUP
C     
          VEIDUP=1
          VEIDDW=1
C       
          NO_ACT=0
C      
      END IF ! IF_DEN=1   
C
C=======================================================================
C              CHI2_6, realted to the Densities rho(r)
C                   (first protons, then neutrons)
C=======================================================================
C      
      IF (IF_RHO.EQ.1) THEN 
C         
C          CHI2_6=CHI2_6+CHIRHO*WEIGHT_RHODEN
C
C          VEIRHO=1
C          NO_ACT=0
C      
      END IF ! IF_RHO=1
C
C=======================================================================
C              CHI2_7, realted to the Level Inversions
C                   (first protons, then neutrons)
C=======================================================================
C      
      IF (IF_INV.EQ.1) THEN 
       
*       DO I=1,LEVTOT
               
*        IF (I.LT.LEVTOT) THEN
               
*             IEXPE1=MEXPER(I)       !Experimental
*             IEXPE2=MEXPER(I+1)
*             DIFEXP=EXPEXP(IEXPE2)-EXPEXP(IEXPE1)
             
*             ITHEO1=MTHEOR(I)       !Theoretical
*             ITHEO2=MTHEOR(I+1)               
*             DIFTHE=ENETHE(ITHEO2)-ENETHE(ITHEO1)
             
*             IF (((DIFEXP.LT.0.).AND.(DIFTHE.GT.0.)).OR.
*     *           ((DIFEXP.GT.0.).AND.(DIFTHE.LT.0.))) THEN  !bad order of the levels
             
*                  AUXIL1=ENETHE(ILWTHE)*DEGENE(I)
*                  AUXIL2=ENETHE(IUPTHE)*DEGENE(I+1)
                  
*                  AUXILI=AUXIL2-AUXIL1
                 
*C                  CHI2_7=CHI2_7+AUXILI*AUXILI*WEIGHT_INVERT
                  
*             END IF
         
*        END IF
             
*       END DO
       
*       CHI2_7=CHI2_7+CHISQU_INVERT
       
*       VEIINV=1
       
*       NO_ACT=0
      
      END IF ! IF_INV=1 
C
C=======================================================================
C           Finally we sum all the contributions to the total chi^2
C=======================================================================
C
      CHOICE=CHI2_1+CHI2_2+CHI2_3+CHI2_4+CHI2_5+CHI2_6+CHI2_7
C
C=======================================================================
C      
      RETURN
      END
C
C=======================================================================
C=======================================================================
C
      SUBROUTINE COUNTI_RHOEXP(INUCLI,NRHOEX)
      
      INCLUDE  'MATDIM/NDNUCL.f'
      INCLUDE  'MATDIM/ND_RHO.f'
      INCLUDE  'MATDIM/N_NOYX.f'
    
      CHARACTER 
     *          NUCSYM*6,SYMBNU*5
      COMMON
     *       /NUCINF/ LDNUCL,
     *                NUCSYM(1:NDNUCL),
     *                NUMB_Z(1:NDNUCL),
     *                NUMB_N(1:NDNUCL)
      COMMON
     *       /EXPDEN/ LABELN(1:N_NOYX,1:2),
     *                R_MESH(1:ND_RHO,1:N_NOYX),
     *                RHOEXP(1:ND_RHO,1:N_NOYX)
     *
     *       /EXPSYM/ SYMBNU(1:N_NOYX)
C
C=======================================================================
C      
      IZ_FIX=NUMB_Z(INUCLI)
      IN_FIX=NUMB_N(INUCLI)
C      
      DO I_NOYX=1,N_NOYX
C          
         IZ_AUX=LABELN(I_NOYX,1)
         IN_AUX=LABELN(I_NOYX,2)
C         
         IF (IZ_FIX.EQ.IZ_AUX .AND. IN_FIX.EQ.IN_AUX) THEN
C             
             DO ID_RHO=1,ND_RHO
                NRHOEX=NRHOEX+1
             END DO
C             
         END IF
C         
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
      SUBROUTINE EXPTHE(IZ_FIX,IN_FIX,DENSUP,DENSDW,FERMEX,GAPEXP,
     *                  IDEFCN,ISOSPI,IPARAM,CHISQU_CHOICE,I_FLAG,
     *                                                     INUCLI)
C
      INCLUDE  'MATDIM/NDSPEC.f'
      INCLUDE  'MATDIM/NDLEXP.f'
      INCLUDE  'MATDIM/NDPARS.f'
      INCLUDE  'MATDIM/NDITEH.f'
      INCLUDE  'MATDIM/NDNUCL.f'
      INCLUDE  'MATDIM/NDISOS.f'
      INCLUDE  'MATDIM/ND_RHO.f'
      INCLUDE  'MATDIM/N_NOYX.f'
C
      CHARACTER
     *          LABEXP*6,LABTHE*6,LWINDW*6,CHIDEF*6
      CHARACTER
     *          LABTHE_PROTON*6,LABTHE_NEUTRS*6
      CHARACTER
     *          LABEXP_PROTON*6,LABEXP_NEUTRS*6
      CHARACTER
     *          NUCLID*10,TYPCHI*6,SYMBNU*5
C
      DIMENSION
     *          CORELE(1:NDLEXP,1:NDLEXP),CORELT(1:NDLEXP,1:NDLEXP)
      DIMENSION
     *          CHNORM(1:NDLEXP),CHWEIG(1:NDLEXP),INDTRU(1:NDLEXP)
      DIMENSION
     *          EWINDW(1:NDLEXP),LWINDW(1:NDLEXP)
      DIMENSION
     *          ENETHE(1:NDSPEC),LABTHE(1:NDSPEC)
      DIMENSION
     *          EXPEXP(1:NDLEXP),IDEGEX(1:NDLEXP)
      DIMENSION
     *          IDGLEV(1:NDSPEC)
      DIMENSION
     *          LABEXP(1:NDLEXP)
      DIMENSION
     *          ICOLEC(1:NDSPEC)
      DIMENSION
     *          NUCLID(0:NDISOS)
      DIMENSION
     *          GAPEXP(1:NDNUCL),FERMEX(1:NDNUCL),
     *          DENSUP(1:NDNUCL),DENSDW(1:NDNUCL)
C
C=======================================================================
C     Follows the full information about  theoretical spherical
C                                                   l e v e l s
C=======================================================================
C
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
C
C=======================================================================
C     Follows the full information about experimental spherical
C                                                   l e v e l s
C=======================================================================
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
     *       /RMSVAL/ HINORM_PROTON,HINORM_NEUTRS
      COMMON
     *       /CHIVAL/ CHISQU_CORREL,DIFSQU_RADIUS,CHISQU_INVERT,
     *                DIFSQU_EFERMI,DIFSQU_ENEGAP,CHIWEI_ENERGY,
     *                CHIWEI_ENEDEG,ERRABS_WEIDEG,EABSAV,ERRMAX,
     *                       DIFSQU_DENSUP,DIFSQU_DENSDW,CHIRHO
      COMMON
     *       /MODIFR/ RMSEXP_PROTON(1:NDNUCL),
     *                RMSEXP_NEUTRS(1:NDNUCL),
     *                RMSTHE_PROTON,RMSTHE_NEUTRS,
     *                RMODIF_PROTON,RMODIF_NEUTRS,
     *                              PUSHIN,DELTAR
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
C                     <<<=== Here:  Storing the minimum of chi2 
C                                   &other discrepancy measures
      COMMON
     *       /CHIMIN/ HIWMIN,HICMIN,ERRMIN,SHEMIN,SHWMIN
C     
C     Starting ('old') and obtained ('new') potential parameters     
C     
      COMMON
     *       /CRITR1/ SCHIWE,SDIFWE,SCHICO,SABSAV,SERMAX,
     *                SDENUP,SDENSU,SDENLO,SDENSD,CHITOT,
     *                                            RADDIF
      COMMON
     *       /CRITR2/ ICHIWE,IDIFWE,ICHICO,IABSAV,IERMAX,
     *                IDENUP,IDENSU,IDENLO,IDENSD,IHITOT
       COMMON
     *       /THEORF/ FERTHE_PROTON,GAPTHE_PROTON,
     *                DENUPP_PROTON,DENLOW_PROTON,
     *
     *                FERTHE_NEUTRS,GAPTHE_NEUTRS,
     *                DENUPP_NEUTRS,DENLOW_NEUTRS 
      COMMON
     *       /INVINV/ INVRSN_PROTON(0:NDITEH),
     *                INVRSN_NEUTRS(0:NDITEH)
C
C=======================================================================
C     Follows the full information about theoretical and experimental
C                                                   d e n s i t i e s
C=======================================================================
C      
      COMMON
     *       /EXPDEN/ LABELN(1:N_NOYX,1:2),
     *                R_MESH(1:ND_RHO,1:N_NOYX),
     *                RHOEXP(1:ND_RHO,1:N_NOYX)
     *
     *       /EXPSYM/ SYMBNU(1:N_NOYX)
      COMMON
     *       /THEDEN/ N_NUCL,I_NUCL(1:N_NOYX),
     *                RHOTHE(1:ND_RHO,1:N_NOYX)
C
C=======================================================================
C
      COMMON
     *       /PRFILE/ ISCREN,LOGFIL,LOGCHP,LOGCHN,LOGHIP,LOGHIN,
     *                LOGPRO,LOGNEU,LOGWRI,LSCREN,IRESUL,ICONVE,
     *                NOUTPT,LOGAUX,LOGENE,LOGRAD,LOGGRD,LOGMSH,
     *                                                   LOGBIS  
      COMMON
     *       /DIMSPH/ LDSING_PROTON,LEVTHE_PROTON,
     *                LEVEXP_PROTON(1:NDNUCL),
     *                LDSING_NEUTRS,LEVTHE_NEUTRS,
     *                LEVEXP_NEUTRS(1:NDNUCL)
      COMMON
     *       /OCCREC/ N_ACTU_PROTON,N_CORR_PROTON,
     *                N_ACTU_NEUTRS,N_CORR_NEUTRS
      COMMON
     *       /POTPOT/ PARPOT(1:NDPARS)
      COMMON
     *       /WHERBE/ IRUNMI,NEWOLD
      COMMON
     *       /MAXIER/ EMAXER
      COMMON
     *       /CNTCNT/ ICOUNT_DENGRD,ICOUNT_FUNMIN,ICOUNT_HAMMAT,
     *                ICOUNT_EXPTHE,ICOUNT_CHOICE,ICOUNT_DIAMAT
C
      DATA
     *      NUCLID(0) / '[Neutrons]' /,
     *      NUCLID(1) / '[Protons] ' /
C
      DATA
     *      HIWLIM / 0.7500 /
      DATA
     *      VALLIM / 1.5000 /
      DATA
     *      EPSLIM / 1.0e-9 /
C
C=======================================================================
C
C     This subroutine compares  the theoretical  and experimental
C     spectra for spherical nuclei and calculates the mean-square
C     deviations
C
C=======================================================================
C
C     Meaning of some parameters:
C
C     HISQUA - CHI-SQUARE as a measure of the difference between the
C              experimental and theoretical sets of energy levels
C
C     ERRMAX - maximum absolute value of the error (theory-experim.)
C
C     EABSAV - arithmetical average of the absolute values of errors
C
C     FERMEX - experimental estimate of the Fermi level
C     FERMTH - theoretical  position of the Fermi level
C
C     GAPEXP - experimental estimate  of the shell closure
C     GAPTHE - theoretical  result  for  the shell closure
C
C     DENSUP - experimental estimate of the level  density A B O V E
C                                             the main shell closure
C     DENSDW - experimental estimate of the level  density B E L O W
C                                             the main shell closure
C
C     DNTHUP - theoretical result for the level density    A B O V E
C                                             the main shell closure
C     DNTHDW - theoretical result for the level  density   B E L O W
C                                             the main shell closure
C
C     N_THUP - the number of theoretical levels  A B O V E  the main
C              closure which appear  too high compared to experiment
C
C     N_THDW - the number of theoretical levels  B E L O W  the main
C              closure which appear  too high compared to experiment
C
C=======================================================================
C
      ICOUNT_EXPTHE=ICOUNT_EXPTHE+1
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
C     
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
C     
      V0CENT_NEUTRS=PARPOT(21)
      R0CENT_NEUTRS=PARPOT(22)
      A0CENT_NEUTRS=PARPOT(23)
      V0SORB_NEUTRS=PARPOT(24)
      R0SORB_NEUTRS=PARPOT(25)
      A0SORB_NEUTRS=PARPOT(26)
      V0EFFM_NEUTRS=PARPOT(27)
      R0EFFM_NEUTRS=PARPOT(28)
      A0EFFM_NEUTRS=PARPOT(29)
C     
      XK_V0C_NEUTRS=PARPOT(30)
      XK_R0C_NEUTRS=PARPOT(31)
      XK_A0C_NEUTRS=PARPOT(32)
      XK_LAM_NEUTRS=PARPOT(33)
      XK_RSO_NEUTRS=PARPOT(34)
      XK_ASO_NEUTRS=PARPOT(35)
      XK_LEF_NEUTRS=PARPOT(36)
      XK_REF_NEUTRS=PARPOT(37)
      XK_AEF_NEUTRS=PARPOT(38)
C     
      ALAMPP       =PARPOT(39)
      ALAMPN       =PARPOT(40)
      ALAMNP       =PARPOT(41)
      ALAMNN       =PARPOT(42)
C
      TLAMPP       =PARPOT(43)
      TLAMPN       =PARPOT(44)
      TLAMNP       =PARPOT(45)
      TLAMNN       =PARPOT(46)
C      
      CLAMPP       =PARPOT(47)
      CLAMPN       =PARPOT(48)
      CLAMNP       =PARPOT(49)
      CLAMNN       =PARPOT(50)
C
C=======================================================================
C
      IF (IDEFCN.EQ.1) THEN
C
C=======================================================================
C         Initialising the chi^2  type minimum control variables
C=======================================================================
C
          HIWMIN=1.0E+10
          HICMIN=1.0E+10
          ERRMIN=1.0E+10
          SHEMIN=1.0E+10
          SHWMIN=1.0E+10
C          
          CHITOT=1.0E+10
C
      END IF
C
C=======================================================================
C
      IF (ISOSPI.EQ.1) THEN
C
          LEVTHE=LDSING_PROTON ! Number of levels below user defined 
C                                                       limit SPHMAX
          LEVEXP=LEVEXP_PROTON(INUCLI)
          RMSEXP=RMSEXP_PROTON(INUCLI)
          RMSTHE=RMSTHE_PROTON
C
          DO ITHEOR=1,LEVTHE
             LABTHE(ITHEOR)=LABTHE_PROTON(ITHEOR)
             ENETHE(ITHEOR)=ENETHE_PROTON(ITHEOR)
             IDGLEV(ITHEOR)=IDGLEV_PROTON(ITHEOR)
          END DO
C          
          DO IEXPER=1,LEVEXP
             LABEXP(IEXPER)=LABEXP_PROTON(INUCLI,IEXPER)
             EXPEXP(IEXPER)=EXPEXP_PROTON(INUCLI,IEXPER)
             IDEGEX(IEXPER)=IDEGEX_PROTON(INUCLI,IEXPER)
          END DO
C          
      END IF 
C
C=======================================================================
C
      IF (ISOSPI.EQ.0) THEN
C
          LEVTHE=LDSING_NEUTRS ! Number of levels below user defined 
C                                                       limit SPHMAX  
          LEVEXP=LEVEXP_NEUTRS(INUCLI)
C
          RMSEXP=RMSEXP_NEUTRS(INUCLI)
          RMSTHE=RMSTHE_NEUTRS
C          
          DO ITHEOR=1,LEVTHE
             LABTHE(ITHEOR)=LABTHE_NEUTRS(ITHEOR)
             ENETHE(ITHEOR)=ENETHE_NEUTRS(ITHEOR)
             IDGLEV(ITHEOR)=IDGLEV_NEUTRS(ITHEOR)
          END DO
C          
          DO IEXPER=1,LEVEXP
             LABEXP(IEXPER)=LABEXP_NEUTRS(INUCLI,IEXPER)
             EXPEXP(IEXPER)=EXPEXP_NEUTRS(INUCLI,IEXPER)
             IDEGEX(IEXPER)=IDEGEX_NEUTRS(INUCLI,IEXPER)
          END DO
C          
      END IF
C
C=======================================================================
C     Calculating the lowest (EMINTH) and the highest (EMAXTH) 
C     theoretical energies and comparing with all the existing 
C     experimental data
C=======================================================================
C
      EMINTH=+1.0E+10
      EMAXTH=-1.0E+10
C
      KTHEOR=0
C
      DO IEXPER=1,LEVEXP
         DO ITHEOR=1,LEVTHE
            ENERGY=ENETHE(ITHEOR)
            IF (LABEXP(IEXPER).EQ.LABTHE(ITHEOR)) THEN
                KTHEOR=KTHEOR+1
                IF (EMINTH.GE.ENERGY) EMINTH=ENERGY
                IF (EMAXTH.LE.ENERGY) EMAXTH=ENERGY
            END IF
         END DO
      END DO
C
      IF (KTHEOR.NE.LEVEXP) THEN
C 
          WRITE(LOGFIL,'(/,''The program does not find all theoretical''
     *                '' levels - according to experimental data !'')')
          WRITE(LOGFIL,'(/,''KTHEOR='',I3,''; LEVEXP='',I3)')
     *                       KTHEOR,          LEVEXP
C
          WRITE(LSCREN,'(/,''The program does not find all theoretical''
     *                '' levels - according to experimental data !'')')
          WRITE(LSCREN,'(/,''KTHEOR='',I3,''; LEVEXP='',I3)')
     *                       KTHEOR,          LEVEXP
          STOP 'STOP in EXPTHE: Missing correspondence with experiment'
C
      END IF
C
      EMINTH=EMINTH-EPSLIM
      EMAXTH=EMAXTH+EPSLIM
C
      DO ITHEOR=1,LEVTHE
         IF (ENETHE(ITHEOR).GT.EMAXTH) THEN
             GO TO 1
         END IF
      END DO
C
   1  CONTINUE
C
      LEVTHE=ITHEOR-1   
C
C=======================================================================
C
C     Attention: The instructions below are meant for the fully
C                closed shells and sub-shells (our experimental
C                data bank contains only those and consequently
C                this should be sufficient for us (so far))
C
C=======================================================================
C
C     Calculating the theoretical gap and Fermi level
C
      IOCCUP=IDGLEV(1)
                       NPARTI=IZ_FIX
      IF (ISOSPI.EQ.0) NPARTI=IN_FIX
C
      DO ITHEOR=2,LEVTHE
C
         IOCCUP=IOCCUP+IDGLEV(ITHEOR)
C
         IF (IOCCUP.GT.NPARTI) THEN
C
             GAPTHE=ENETHE(ITHEOR)-ENETHE(ITHEOR-1)
             FERTHE=0.5*(ENETHE(ITHEOR)+ENETHE(ITHEOR-1))
C
             GO TO 2
C
         END IF
C
      END DO
C
C=======================================================================
C
   2  CONTINUE
C
      IOCCUP=IOCCUP-IDGLEV(ITHEOR)
C
      N_ACTU=NPARTI
      N_CORR=NPARTI
C
      IF (IOCCUP.NE.NPARTI) THEN
C
          N_ACTU=IOCCUP
C
          IF (LOGWRI.GT.4) THEN
C
              DO I=1,20          
              WRITE(LOGFIL,'(''Occupation number ='',i3,'' in EXPTHE '',
     *                       ''should be equal to NPARTI='',i3)') 
     *                                     IOCCUP,NPARTI
              END DO
C
C             STOP 'Occupation number inconsistent, STOP from EXPTHE'
C
          END IF
C
      ELSE
C
          IF (LOGWRI.GT.5) THEN
C
              WRITE(LOGFIL,'(18X,''Occupation number ='',i3,1x,
     *                           ''should be equal to NPARTI='',i3)') 
     *                                         IOCCUP,NPARTI
          END IF
C
      END IF
C
C=======================================================================
C
      IPRIN1=0
      IPRIN2=0
      IPRIN3=0
C
C=======================================================================
C
      DO IEXPER=1,LEVEXP
         INDTRU(IEXPER)=0
      END DO
C
      ITHUPP=0
      ITHLOW=0
      IDEGER=0
C
      ELOMIN=+10000.
      EUPMIN=+10000.
C
      ELOMAX=-10000.
      EUPMAX=-10000.
C
      VALMAX=-1.0E+10
C
      DIFSUM=0.0
      CHISUM=0.0
C
      CHIWEI_ENERGY=0.0
      CHIWEI_ENEDEG=0.0
      ERRABS_WEIDEG=0.0
C
      EABSAV=0.0
C
      IHISUM=0
C
      DO ITHEOR=1,LEVTHE
C
         ETHEOR=ENETHE(ITHEOR)
C
         DO IEXPER=1,LEVEXP
C
            IF (LABEXP(IEXPER).EQ.LABTHE(ITHEOR)) THEN
C
                INDTRU(IEXPER)=1
C
                IHISUM=IHISUM+1
                EWINDW(IHISUM)=ETHEOR                
                LWINDW(IHISUM)=LABTHE(ITHEOR)             
C
                IDEGER=IDEGER+IDEGEX(IEXPER)
C
                ENEDIF=ETHEOR-EXPEXP(IEXPER)
                ESQUAR=ENEDIF*ENEDIF
C
                ERRABS=ABS(ENEDIF)
C
                EABSAV=EABSAV+ERRABS
C
                CHISUM=CHISUM+ESQUAR
                CHIWEI_ENERGY=CHIWEI_ENERGY+ESQUAR*IDEGEX(IEXPER)
C
                DIFSUM=DIFSUM+ERRABS
                ERRABS_WEIDEG=ERRABS_WEIDEG+ERRABS*IDEGEX(IEXPER)
C
                IF (ERRABS.GT.VALMAX) VALMAX=ERRABS
C
                CHNORM(IEXPER)=SQRT(CHISUM/IHISUM)
                CHWEIG(IEXPER)=SQRT(CHIWEI_ENERGY/IDEGER)
C
C=======================================================================
C               Below, we calculate average level densities
C               above the Fermi level (particle states) and
C               below the Fermi level (hole states)
C=======================================================================
C
                IF (ETHEOR.GT.FERTHE) THEN
C                                                      Particles:
                    ITHUPP=ITHUPP+IDGLEV(ITHEOR)
C
                    IF (EUPMIN.GT.ETHEOR) EUPMIN=ETHEOR
                    IF (EUPMAX.LT.ETHEOR) EUPMAX=ETHEOR
C
                ELSE
C                                                      Holes:
                    ITHLOW=ITHLOW+IDGLEV(ITHEOR)
C
                    IF (ELOMIN.GT.ETHEOR) ELOMIN=ETHEOR
                    IF (ELOMAX.LT.ETHEOR) ELOMAX=ETHEOR
C
                END IF
C
C=======================================================================
C            
            END IF
C
         END DO
C
      END DO
C
      INWIND=IHISUM 
C
      HWNORM=SQRT(CHIWEI_ENERGY/IDEGER)
      HINORM=SQRT(CHISUM/IHISUM)
C
      EABSAV=EABSAV/LEVEXP
      ERRMAX=VALMAX
      
      EMAXER=ERRMAX ! changed for LMMINI
C
C=======================================================================
C     Calculating and storing the inversions: theory vs. experiment
C=======================================================================
C
      ITHBIS=0
      CHISQU_INVERT=0.0
C
      DO ITHEOR=1,INWIND
C
         ETHEOR=EWINDW(ITHEOR)
C
         DO IEXPER=1,LEVEXP
C
            IF (LABEXP(IEXPER).EQ.LWINDW(ITHEOR)) THEN
C
                ITHBIS=ITHBIS+1
C
                ICOLEC(ITHBIS)=IEXPER
C
                IF (IEXPER.NE.ITHEOR) THEN
C
                    CHISQU_INVERT=CHISQU_INVERT
     *                           +
     *                           (ETHEOR-EXPEXP(IEXPER))**2
     *                           *IDEGEX(IEXPER)
                END IF
C
            END IF	    
C
         END DO
C
      END DO
C
      IF (ITHBIS.NE.LEVEXP) THEN
C
          WRITE(LSCREN,'(/,''Alarm in EXPTHE with INUCLI= '',I2,
     *                     '' and ISOSPI= '',I1,
     *                     '' : ITHBIS= '',I3,'' and LEVEXP= '',I3,
     *                     '' ==> They should be equal!'',/)')
     *                        INUCLI,ISOSPI,ITHBIS,LEVEXP
C
          STOP 'STOP in EXPTHE: ITHBIS.NE.LEVEXP'  
C  
      END IF
C
C=======================================================================
C     Estimating the number of inversions theory vs. experiment
C=======================================================================
C
      CALL INVERT(ICOLEC,NDSPEC,ITHBIS,NOFINV)
C
      IF (ISOSPI.EQ.1) THEN
          INVRSN_PROTON(IDEFCN)=NOFINV
      END IF
C
      IF (ISOSPI.EQ.0) THEN
          INVRSN_NEUTRS(IDEFCN)=NOFINV
      END IF
C
C=======================================================================
C     Defining the 'normalised' final expression of the 'chi2' 
C=======================================================================
C
      CHIWEI_ENEDEG=CHIWEI_ENERGY
      CHIWEI_ENERGY=CHIWEI_ENERGY/IDEGER
      CHINOR=SQRT(CHISUM/IHISUM)
C
C=======================================================================
C     Calculating energy densities below and above the Fermi level
C=======================================================================
C
      ERRABS_WEIDEG=ERRABS_WEIDEG/IDEGER
C
      DENLOW=ITHLOW/(ELOMAX-ELOMIN)
      DENUPP=ITHUPP/(EUPMAX-EUPMIN)
C
C=======================================================================
C=======================================================================
C     Calculating the correlation matrices for experiment and theory
C=======================================================================
C=======================================================================
C
      DO IEXPEL=1,LEVEXP
         DO IEXPER=1,LEVEXP
            CORELE(IEXPEL,IEXPER)=(EXPEXP(IEXPEL)-EXPEXP(IEXPER))
         END DO
      END DO
C
      DO ITHEOL=1,INWIND
         DO ITHEOR=1,INWIND
            CORELT(ITHEOL,ITHEOR)=(EWINDW(ITHEOL)-EWINDW(ITHEOR))
         END DO
      END DO
C
C=======================================================================
C     Below taking into account level degeneracies when calculating
C     correlation matrices and the corresponding chi^2 contribution 
C=======================================================================
C
      CHISQU_CORREL=0.0
      XNORMA=0.0
C
      DO IEXPEL=1,LEVEXP
         DO IEXPER=1,LEVEXP
C
            DEGDEG=IDEGEX(IEXPEL)*IDEGEX(IEXPER)
            DEGDEG=SQRT(DEGDEG)
            DIFFER=CORELT(IEXPEL,IEXPER)-CORELE(IEXPEL,IEXPER)
            DIFFER=DIFFER*DIFFER
C
            CHISQU_CORREL=CHISQU_CORREL+DEGDEG*DIFFER
            XNORMA=XNORMA+DEGDEG
C
         END DO
      END DO
C
      CHISQU_CORREL=SQRT(0.5*CHISQU_CORREL/XNORMA)
C
C=======================================================================
C     Quantities below are used in the minimisation procedure
C=======================================================================
C
      DIFSQU_EFERMI=(FERMEX(INUCLI)-FERTHE)**2
      DIFSQU_ENEGAP=(GAPEXP(INUCLI)-GAPTHE)**2
      DIFSQU_DENSUP=(DENUPP-DENSUP(INUCLI))**2
      DIFSQU_DENSLW=(DENLOW-DENSDW(INUCLI))**2
C
      DIFSQU_RADIUS=(RMSEXP-RMSTHE)**2
C
C=======================================================================
C     Calculating the proposed shift of the central radius r0
C=======================================================================
C
      CALL SHIFTR(INUCLI,IZ_FIX,IN_FIX,ISOSPI)
C
C=======================================================================
C     Calculating the \rho(r) \chi^2
C=======================================================================
*C      
*      CHIRHO=0.0
*C          
*      DO I_AUXI=1,N_NUCL
*C             
*         I_NOYX=I_NUCL(I_AUXI)
*C             
*         DO ID_RHO=1,ND_RHO
                
*            CHIRHO=CHIRHO+(RHOTHE(ID_RHO,I_NOYX)
*     *                   - RHOEXP(ID_RHO,I_NOYX))**2
*         END DO
*C          
*      END DO
C
C=======================================================================
C     Calculating the actual (weighted) chi^2 function value
C=======================================================================
C
      HISQUA=CHOICE(ISOSPI,IDEFCN,INUCLI)
      CHISQU_CHOICE=HISQUA
C
C=======================================================================
C
      IF (HICMIN.GT.CHISQU_CORREL) THEN
          HICMIN=CHISQU_CORREL
          IPRIN3=1
      END IF
C
      IF (HIWMIN.GT.CHIWEI_ENERGY) THEN
          HIWMIN=CHIWEI_ENERGY
          IPRIN1=1
      END IF
C
      IF (CHIWEI_ENERGY.LE.HIWLIM) IPRIN1=1
C
      IF (VALMAX.LE.VALLIM) IPRIN2=1
C
      IF (ERRMIN.GT.VALMAX) THEN
          ERRMIN=VALMAX
          IPRIN2=1
      END IF
C
C=======================================================================
C
      IF (CHITOT.GT.HISQUA) THEN
          CHITOT=HISQUA
          IHITOT=IDEFCN
          RADDIF=SQRT(DIFSQU_RADIUS)
          IPRIN3=1
      END IF
C
      IPRIN3=1
C
C=======================================================================
C     If the results on chi^2 and maximum error were improved, 
C     we print them ...
C=======================================================================
C
      IF (I_FLAG.EQ.1) THEN
C
          IF ((IPRIN1.EQ.1).OR.(IPRIN2.EQ.1).OR.(IPRIN3.EQ.1)) THEN
C
              CALL CMPRSN(INUCLI,IZ_FIX,IN_FIX,ISOSPI,NDISOS,NUCLID,
     *                    IRUNMI,IDEFCN,EMINTH,EMAXTH,R0COUL,CHNORM,
     *                    FERMEX,FERTHE,GAPEXP,GAPTHE,CHIWEI_ENERGY,
     *                    ERRABS_WEIDEG,CHISQU_CORREL,EABSAV,ERRMAX,
     *                    DENUPP,DENSUP,DENLOW,DENSDW,CHWEIG,
     *                    CHISQU_INVERT,CHISQU_CHOICE,INDTRU,IPARAM)
C
          END IF
C
      END IF
C
C=======================================================================
C     Collecting the minimum values of the criterion functions
C=======================================================================
C
C     ... first the total chi^2 i.e. the one defined by CHOICE
C
C     IF (CHITOT.GT.HISQUA) THEN
C         CHITOT=HISQUA
C         IHITOT=IDEFCN
C     END IF
C
C         ... then the rest ...
C
      IF (SCHIWE.GT.CHIWEI_ENERGY) THEN
          SCHIWE=CHIWEI_ENERGY
          ICHIWE=IDEFCN
      END IF
C
      IF (SDIFWE.GT.ERRABS_WEIDEG) THEN
          SDIFWE=ERRABS_WEIDEG
          IDIFWE=IDEFCN
      END IF
C      
      IF (SCHICO.GT.CHISQU_CORREL) THEN
          SCHICO=CHISQU_CORREL
          ICHICO=IDEFCN
      END IF
C      
      IF (SABSAV.GT.EABSAV) THEN
          SABSAV=EABSAV
          IABSAV=IDEFCN
      END IF
      
      IF (SERMAX.GT.ERRMAX) THEN
          SERMAX=ERRMAX
          IERMAX=IDEFCN
      END IF
      
      IF (SDENUP.GT.DENUPP) THEN
          SDENUP=DENUPP
          IDENUP=IDEFCN
      END IF
      
      IF (SDENSU.GT.DENSUP(INUCLI)) THEN
          SDENSU=DENSUP(INUCLI)
          IDENSU=IDEFCN
      END IF
      
      IF (SDENLO.GT.DENLOW) THEN
          SDENLO=DENLOW
          IDENLO=IDEFCN
      END IF
      
      IF (SDENSD.GT.DENSDW(INUCLI)) THEN
          SDENSD=DENSDW(INUCLI)
          IDENSD=IDEFCN
      END IF
C
C=======================================================================
C     Storing the results on disk
C=======================================================================
C
      IF (ISOSPI.EQ.1) THEN
C      
          CHICOR_PROTON=CHISQU_CORREL
          RADDIF_PROTON=DIFSQU_RADIUS
          CHIINV_PROTON=CHISQU_INVERT
          FERDIF_PROTON=DIFSQU_EFERMI
          GAPDIF_PROTON=DIFSQU_ENEGAP
          CHIWEI_PROTON=CHIWEI_ENERGY
          CHIDEG_PROTON=CHIWEI_ENEDEG
          EABSWD_PROTON=ERRABS_WEIDEG
          EABSAV_PROTON=EABSAV
          ERRMAX_PROTON=ERRMAX
          DIFFUP_PROTON=DIFSQU_DENSUP
          DIFFDW_PROTON=DIFSQU_DENSDW
          N_ACTU_PROTON=N_ACTU
          N_CORR_PROTON=N_CORR
C
          GAPTHE_PROTON=GAPTHE
          FERTHE_PROTON=FERTHE
          DENLOW_PROTON=DENLOW
          DENUPP_PROTON=DENUPP
C
          HINORM_PROTON=HINORM
C
          CHIRHO_PROTON=CHIRHO
C
          IF (LOGWRI.GT.0) THEN
              WRITE(NOUTPT,'(80(''#''),/,''#'',T80,''#'',/,
     *                   ''#  Proton Energy RMS= '',f15.6,T80,''#'',/,
     *                   ''#'',T80,''#'',/,80(''#''))')
     *                              SQRT(CHIWEI_PROTON)
          END IF
C
      END IF          
C
C=======================================================================
C
      IF (ISOSPI.EQ.0) THEN
C      
          CHICOR_NEUTRS=CHISQU_CORREL
          RADDIF_NEUTRS=DIFSQU_RADIUS
          CHIINV_NEUTRS=CHISQU_INVERT
          FERDIF_NEUTRS=DIFSQU_EFERMI
          GAPDIF_NEUTRS=DIFSQU_ENEGAP
          CHIWEI_NEUTRS=CHIWEI_ENERGY
          CHIDEG_NEUTRS=CHIWEI_ENEDEG
          EABSWD_NEUTRS=ERRABS_WEIDEG
          EABSAV_NEUTRS=EABSAV
          ERRMAX_NEUTRS=ERRMAX
          DIFFUP_NEUTRS=DIFSQU_DENSUP
          DIFFDW_NEUTRS=DIFSQU_DENSDW
          N_ACTU_NEUTRS=N_ACTU
          N_CORR_NEUTRS=N_CORR
C
          GAPTHE_NEUTRS=GAPTHE
          FERTHE_NEUTRS=FERTHE
          DENLOW_NEUTRS=DENLOW
          DENUPP_NEUTRS=DENUPP
C
          HINORM_NEUTRS=HINORM
C
          IF (LOGWRI.GT.0) THEN
              WRITE(NOUTPT,'(80(''#''),/,''#'',T80,''#'',/,
     *                   ''#  Neutron Energy RMS= '',f15.6,T80,''#'',/,
     *                   ''#'',T80,''#'',/,80(''#''))')
     *                              SQRT(CHIWEI_NEUTRS)
          END IF
C
      END IF 
C
C=======================================================================
C
                
C
C=======================================================================
C
      IMODUL=IDEFCN-1
C
      IF ((IDEFCN.EQ.1).OR.(MOD(IMODUL,20).EQ.0)) THEN
C
          IF (ISOSPI.EQ.1) THEN
C
              IF (LOGWRI.GT.0) THEN               
C
                  WRITE(LOGCHP,'(/,''From EXPTHE'',/)')
C
                  WRITE(LOGCHP,'(/,80(''#''),/,''#'',78X,''#'',/,
     *                ''#  Test of the quality of the spherical '',
     *                ''s.p. spectrum  PROTONS Z ='',I3,
     *                ''  N ='',I3,''  #'',/,''#'',78X,''#'',/,
     *                ''#  No:  HiWei  HiCor  HiNor  ErrAv  ErrMx'',
     *                ''  GapT   D_UpT  D_DwT   R_th    F_th  #'',/,
     *                ''#                                 Ex -->>'',
     *                F6.2,F8.3,F7.3,F7.3,F8.2,''  #'')')
     *
     *                IZ_FIX,IN_FIX,GAPEXP(INUCLI),DENSUP(INUCLI),
     *                DENSDW(INUCLI),RMSEXP,FERMEX(INUCLI)
              END IF
C
          END IF
C
          IF (ISOSPI.EQ.0) THEN
C
              IF (LOGWRI.GT.0) THEN
C
                  WRITE(LOGCHN,'(/,''From EXPTHE'',/)')
C                        
                  WRITE(LOGCHN,'(/,80(''#''),/,''#'',78X,''#'',/,
     *                ''#  Test of the quality of the spherical '',
     *                ''s.p. spectrum NEUTRONS Z ='',I3,
     *                ''  N ='',I3,''  #'',/,''#'',78X,''#'',/,
     *                ''#  No:  HiWei  HiCor  HiNor  ErrAv  ErrMx'',
     *                ''  GapT   D_UpT  D_DwT   R_th    F_th  #'',/,
     *                ''#                                 Ex -->>'',
     *                F6.2,F8.3,F7.3,F7.3,F8.2,''  #'')')
     *
     *                IZ_FIX,IN_FIX,GAPEXP(INUCLI),DENSUP(INUCLI),
     *                DENSDW(INUCLI),RMSEXP,FERMEX(INUCLI)
              END IF
C              
          END IF
C
      END IF
C
C=======================================================================
C
      IF (LOGWRI.GT.0) THEN
C
          IF (ISOSPI.EQ.1) THEN
C                              
              WRITE(LOGCHP,'(''#'',I5,5F7.3,F6.2,F8.3,2F7.3,F8.2,
     *                                                 ''  #'')')
     *
     *           IDEFCN,CHIWEI_ENERGY,CHISQU_CORREL,CHINOR,EABSAV,
     *                  ERRMAX,GAPTHE,DENUPP,DENLOW,RMSTHE,FERTHE
          END IF
C
          IF (ISOSPI.EQ.0) THEN
C                              
              WRITE(LOGCHN,'(''#'',I5,5F7.3,F6.2,F8.3,2F7.3,F8.2,
     *                                                 ''  #'')')
     *
     *           IDEFCN,CHIWEI_ENERGY,CHISQU_CORREL,CHINOR,EABSAV,
     *                  ERRMAX,GAPTHE,DENUPP,DENLOW,RMSTHE,FERTHE
          END IF
C
      END IF
C
C=======================================================================
C
      IF (LOGWRI.GT.5) THEN
          WRITE(LOGFIL,'(15X,''Exiting  EXPTHE'')')
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
      SUBROUTINE SHIFTR(INUCLI,IZ_FIX,IN_FIX,ISOSPI)
C          
      INCLUDE  'MATDIM/NDNUCL.f'
C
      COMMON
     *       /MODIFR/ RMSEXP_PROTON(1:NDNUCL),
     *                RMSEXP_NEUTRS(1:NDNUCL),
     *                RMSTHE_PROTON,RMSTHE_NEUTRS,
     *                RMODIF_PROTON,RMODIF_NEUTRS,
     *                              PUSHIN,DELTAR
      COMMON
     *       /DERIVD/ IDERID
      COMMON
     *       /PRINAL/ IMTALK
      COMMON
     *       /PRFILE/ ISCREN,LOGFIL,LOGCHP,LOGCHN,LOGHIP,LOGHIN,
     *                LOGPRO,LOGNEU,LOGWRI,LSCREN,IRESUL,ICONVE,
     *                NOUTPT,LOGAUX,LOGENE,LOGRAD,LOGGRD,LOGMSH,
     *                                                   LOGBIS 
C
C=======================================================================      
C     This function performs dynamical shift of the central
C     radius constant when minimising the ch2^2
C=======================================================================      
C
C     PUSHIN - A small proportion of the actual discrepancy
C              between Rtheo and Rexp;
C              expected to be generally between 0.1 and 0.5
C
C     RMODIF - Accumulated modification of the radius; each
C              time the new iteration has been performed it
C              will be increased  as long as the difference
C                         DELTAR=RMSTHE-RMSEXP
C              is positive - otherwise it will be decreased
C
C=======================================================================      
C
      IF (ISOSPI.EQ.1) THEN
         RMSTHE=RMSTHE_PROTON
	 RMSEXP=RMSEXP_PROTON(INUCLI)
      END IF
C      
      IF (ISOSPI.EQ.0) THEN
         RMSTHE=RMSTHE_NEUTRS
	 RMSEXP=RMSEXP_NEUTRS(INUCLI)
      END IF
C
      IF (IDERID.EQ.1) THEN
          RMODIF=0
      ELSE
C
          A_MASS=REAL(IZ_FIX+IN_FIX)
          FACTOR=A_MASS**(1.0/3.0)
C
C         Below we calculate an estimate of the central-radius par.
C                                                  modification
          DELTAR=RMSTHE-RMSEXP
          DELTAR=DELTAR/FACTOR
C
          RMODIF=PUSHIN*DELTAR
C
          IF (ISOSPI.EQ.1) THEN
	  
	     RMODIF_PROTON=RMODIF
	  
             IF (IMTALK.EQ.1 .AND. LOGWRI.GT.0) THEN
                 WRITE(NOUTPT,'(''In SHIFTR (P): RMSEXP='',f7.4,
     *		               '' RMSTHE='',f10.7,'' DELTAR = '',F14.7,
     *                             '' RMODIF = '',F14.7)') 
     *                     RMSEXP,RMSTHE,DELTAR,RMODIF
             END IF
	  END IF
C
          IF (ISOSPI.EQ.0) THEN
	  
	     RMODIF_NEUTRS=RMODIF
	  
             IF (IMTALK.EQ.1 .AND. LOGWRI.GT.0) THEN
                 WRITE(NOUTPT,'(''In SHIFTR (N): RMSEXP='',f7.4,
     *		             '' RMSTHE='',f10.7,'' DELTAR = '',F14.7,
     *                             '' RMODIF = '',F14.7)') 
     *                     RMSEXP,RMSTHE,DELTAR,RMODIF
             END IF
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
      SUBROUTINE INVERT(INVERS,NDSPEC,LEVACT,NOFINV)
      DIMENSION
     *          INVERS(1:NDSPEC)
C
      COMMON
     *       /PRINAL/ IMTALK
      COMMON
     *       /PRFILE/ ISCREN,LOGFIL,LOGCHP,LOGCHN,LOGHIP,LOGHIN,
     *                LOGPRO,LOGNEU,LOGWRI,LSCREN,IRESUL,ICONVE,
     *                NOUTPT,LOGAUX,LOGENE,LOGRAD,LOGGRD,LOGMSH,
     *                                                   LOGBIS
C
C=======================================================================
C     This routine  e s t i m a t e s   the number of inversions
C     of the theoretical levels with respect to the experimental 
C     order
C=======================================================================
C
      IF (LEVACT.LT.1.OR.LEVACT.GT.NDSPEC) THEN
C
          IF (IMTALK.EQ.1) THEN
              WRITE(NOUTPT,'(/,''LEVACT='',I5,'' is out of range, '',
     *                    ''NDSPEC='',I5)') LEVACT,NDSPEC
          END IF
C      
          LEVACT=NDSPEC
C
      END IF          
C
C=======================================================================
C
      NOFINV=0
C
      IF (LEVACT.EQ.1) RETURN
C
      IF (LEVACT.EQ.2.AND.INVERS(1).GT.INVERS(2)) NOFINV=1      
C
C=======================================================================
C
      IF (LEVACT.GT.2) THEN
C
          NOFINV=0      
C      
          DO I=1,LEVACT-1
             I_DIFF=INVERS(I+1)-INVERS(I)
             IF (I_DIFF.NE.1) NOFINV=NOFINV+1
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
C=======================================================================
C
      SUBROUTINE CMPRSN(INUCLI,IZ_FIX,IN_FIX,ISOSPI,NDISOS,NUCLID,
     *                  IRUNMI,IDEFCN,EMINTH,EMAXTH,RCOULO,CHNORM,
     *                  FERMEX,FERTHE,GAPEXP,GAPTHE,CHIWEI,DIFWEI,
     *                  CHICOR,EABSAV,ERRMAX,DENUPP,DENSUP,DENLOW,
     *                  DENSDW,CHWEIG,CHIINV,CHISQU_CHOICE,INDTRU,
     *                                                     IPARAM)
C
      INCLUDE  'MATDIM/NDSPEC.f'
      INCLUDE  'MATDIM/NDLEXP.f'
      INCLUDE  'MATDIM/NDITEH.f'
      INCLUDE  'MATDIM/NDPARS.f'
      INCLUDE  'MATDIM/NDNUCL.f'
C
      CHARACTER
     *          LABEXP*6,LABTHE*6,CHIDEF*6,NUCLID*10
      CHARACTER
     *          LABEXP_PROTON*6,LABEXP_NEUTRS*6,
     *          LABTHE_PROTON*6,LABTHE_NEUTRS*6
C
      DIMENSION
     *          CHNORM(1:NDLEXP),CHWEIG(1:NDLEXP),INDTRU(1:NDLEXP)
      DIMENSION
     *          ENETHE(1:NDSPEC),
     *          LABTHE(1:NDSPEC)
      DIMENSION
     *          EXPEXP(1:NDLEXP),
     *          LABEXP(1:NDLEXP)
      DIMENSION
     *          NUCLID(0:NDISOS)
      DIMENSION
     *          GAPEXP(1:NDNUCL),FERMEX(1:NDNUCL),
     *          DENSUP(1:NDNUCL),DENSDW(1:NDNUCL)
C     
      COMMON
     *       /MODIFR/ RMSEXP_PROTON(1:NDNUCL),
     *                RMSEXP_NEUTRS(1:NDNUCL),
     *                RMSTHE_PROTON,RMSTHE_NEUTRS,
     *                RMODIF_PROTON,RMODIF_NEUTRS,
     *                              PUSHIN,DELTAR
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
     *       /PRFILE/ ISCREN,LOGFIL,LOGCHP,LOGCHN,LOGHIP,LOGHIN,
     *                LOGPRO,LOGNEU,LOGWRI,LSCREN,IRESUL,ICONVE,
     *                NOUTPT,LOGAUX,LOGENE,LOGRAD,LOGGRD,LOGMSH,
     *                                                   LOGBIS                                                  
      COMMON
     *       /POTPOT/ PARPOT(1:NDPARS)
      COMMON
     *       /DIMSPH/ LDSING_PROTON,LEVTHE_PROTON,
     *                LEVEXP_PROTON(1:NDNUCL),
     *                LDSING_NEUTRS,LEVTHE_NEUTRS,
     *                LEVEXP_NEUTRS(1:NDNUCL) 
      COMMON
     *       /SPECTR/ ENETHE_PROTON(1:NDSPEC),
     *                ENETHE_NEUTRS(1:NDSPEC)
      COMMON
     *       /THELAB/ LABTHE_PROTON(1:NDSPEC),
     *                LABTHE_NEUTRS(1:NDSPEC)
      COMMON
     *       /INVINV/ INVRSN_PROTON(0:NDITEH),
     *                INVRSN_NEUTRS(0:NDITEH)
      COMMON
     *       /ACTIVE/ IACTIV,IwMODE
      COMMON
     *       /WHICHI/ CHIDEF
C
C=======================================================================
C     Printing routine only - but what does it do ???
C=======================================================================
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
C
      IF (ISOSPI.EQ.1) THEN
          LEVTHE_PROTON=LDSING_PROTON
C
          LEVTHE=LDSING_PROTON ! Number of levels below user defined 
C                                  limit SPHMAX
          LEVEXP=LEVEXP_PROTON(INUCLI)
c	  
	  V0CENT=V0CENT_PROTON
	  R0CENT=R0CENT_PROTON
	  A0CENT=A0CENT_PROTON
	  V0SORB=V0SORB_PROTON
	  R0SORB=R0SORB_PROTON
	  A0SORB=A0SORB_PROTON
	  V0EFFM=V0EFFM_PROTON
	  R0EFFM=R0EFFM_PROTON
	  A0EFFM=A0EFFM_PROTON
	  XK_V0C=XK_V0C_PROTON
	  XK_R0C=XK_R0C_PROTON
	  XK_A0C=XK_A0C_PROTON
     	  XK_LAM=XK_LAM_PROTON
          XK_RSO=XK_RSO_PROTON
	  XK_ASO=XK_ASO_PROTON
	  XK_LEF=XK_LEF_PROTON
	  XK_REF=XK_REF_PROTON
	  XK_AEF=XK_AEF_PROTON
C
          RMSEXP=RMSEXP_PROTON(INUCLI)
	  RMSTHE=RMSTHE_PROTON
	  RMODIF=RMODIF_PROTON
C
          DO ITHEOR=1,LEVTHE
             LABTHE(ITHEOR)=LABTHE_PROTON(ITHEOR)
             ENETHE(ITHEOR)=ENETHE_PROTON(ITHEOR)
          END DO
C          
          DO IEXPER=1,LEVEXP
             LABEXP(IEXPER)=LABEXP_PROTON(INUCLI,IEXPER)
             EXPEXP(IEXPER)=EXPEXP_PROTON(INUCLI,IEXPER)
          END DO
C          
      END IF 
C
      IF (ISOSPI.EQ.0) THEN
          LEVTHE_NEUTRS=LDSING_NEUTRS
C
          LEVTHE=LDSING_NEUTRS  
          LEVEXP=LEVEXP_NEUTRS(INUCLI)
C          
	  V0CENT=V0CENT_NEUTRS
	  R0CENT=R0CENT_NEUTRS
	  A0CENT=A0CENT_NEUTRS
	  V0SORB=V0SORB_NEUTRS
	  R0SORB=R0SORB_NEUTRS
	  A0SORB=A0SORB_NEUTRS
	  V0EFFM=V0EFFM_NEUTRS
	  R0EFFM=R0EFFM_NEUTRS
	  A0EFFM=A0EFFM_NEUTRS
	  XK_V0C=XK_V0C_NEUTRS
	  XK_R0C=XK_R0C_NEUTRS
	  XK_A0C=XK_A0C_NEUTRS
     	  XK_LAM=XK_LAM_NEUTRS
          XK_RSO=XK_RSO_NEUTRS
	  XK_ASO=XK_ASO_NEUTRS
	  XK_LEF=XK_LEF_NEUTRS
	  XK_REF=XK_REF_NEUTRS
	  XK_AEF=XK_AEF_NEUTRS
C
	  RMSEXP=RMSEXP_NEUTRS(INUCLI)
	  RMSTHE=RMSTHE_NEUTRS
	  RMODIF=RMODIF_NEUTRS
C
          DO ITHEOR=1,LEVTHE
             LABTHE(ITHEOR)=LABTHE_NEUTRS(ITHEOR)
             ENETHE(ITHEOR)=ENETHE_NEUTRS(ITHEOR)
          END DO
C          
          DO IEXPER=1,LEVEXP
             LABEXP(IEXPER)=LABEXP_NEUTRS(INUCLI,IEXPER)
             EXPEXP(IEXPER)=EXPEXP_NEUTRS(INUCLI,IEXPER)
          END DO
C          
      END IF
C
      CALL PRINTI(IDEFCN)
C
      IF (ISOSPI.EQ.0 .AND. LOGWRI.GT.4) THEN
C
          WRITE(IRESUL,'(80(''#''),/,''#'',78X,''#'',/,
     *                    ''#  Comparison between the calculated '',
     *                    ''and experimental results:  Z='',I2,
     *                    ''   N='',I3,''   #'',/,
     *                    ''#  '',A10,''                            '',
     *                    ''                RST='',I2,
     *                    '' IDEFCN='',I5,''   #'',/,
     *                    ''#'',78X,''#'',/,
     *                    80(''#''),/,''#'',78X,''#'',/,
     *                    ''#   No)   E_calc    E_expe    Lab_th    '',
     *                    ''Lab_exp    Eth-Eex   Chisqr   ChisqW   #'',
     *                    /,''#'',78X,''#'')')
     *
     *                                    IZ_FIX,IN_FIX,NUCLID(ISOSPI),
     *                                                  IRUNMI,IDEFCN
C
          IF (IwMODE.EQ.2 .AND. LOGWRI.GT.4)
     *
     *    WRITE(LOGNEU,'(/,80(''#''),/,''#'',78X,''#'',/,
     *                    ''#  Comparison between the calculated '',
     *                    ''and experimental results:  Z='',I2,
     *                    ''   N='',I3,''   #'',/,
     *                    ''#  '',A10,''                            '',
     *                    ''                RST='',I2,
     *                    '' IDEFCN='',I5,''   #'',/,
     *                    ''#'',78X,''#'',/,
     *                    80(''#''),/,''#'',78X,''#'',/,
     *                    ''#   No)   E_calc    E_expe    Lab_th    '',
     *                    ''Lab_exp    Eth-Eex   Chisqr   ChisqW   #'',
     *                    /,''#'',78X,''#'')')
     *
     *                                    IZ_FIX,IN_FIX,NUCLID(ISOSPI),
     *                                                  IRUNMI,IDEFCN
C
      END IF
C
      IF (ISOSPI.EQ.1 .AND. LOGWRI.GT.4) THEN
C
          WRITE(IRESUL,'(80(''#''),/,''#'',78X,''#'',/,
     *                    ''#  Comparison between the calculated '',
     *                    ''and experimental results:  Z='',I2,
     *                    ''   N='',I3,''   #'',/,
     *                    ''#  '',A10,'' {Coulomb Effective Radius ='',
     *                    F6.3,''*A^(1/3)}'','' RST='',I2,
     *                    '' IDEFCN='',I5,''   #'',/,''#'',78X,''#'',/,
     *                    80(''#''),/,''#'',78X,''#'',/,
     *                    ''#   No)   E_calc    E_expe    Lab_th    '',
     *                    ''Lab_exp    Eth-Eex   Chisqr   ChisqW   #'',
     *                    /,''#'',78X,''#'')')
     *
     *                                    IZ_FIX,IN_FIX,NUCLID(ISOSPI),
     *                                           RCOULO,IRUNMI,IDEFCN
C
          IF (IwMODE.EQ.2 .AND. LOGWRI.GT.4)
     *
     *    WRITE(LOGPRO,'(/,80(''#''),/,''#'',78X,''#'',/,
     *                    ''#  Comparison between the calculated '',
     *                    ''and experimental results:  Z='',I2,
     *                    ''   N='',I3,''   #'',/,
     *                    ''#  '',A10,'' {Coulomb Effective Radius ='',
     *                    F6.3,''*A^(1/3)}'','' RST='',I2,
     *                    '' IDEFCN='',I5,''   #'',/,''#'',78X,''#'',/,
     *                    80(''#''),/,''#'',78X,''#'',/,
     *                    ''#   No)   E_calc    E_expe    Lab_th    '',
     *                    ''Lab_exp    Eth-Eex   Chisqr   ChisqW   #'',
     *                    /,''#'',78X,''#'')')
     *
     *                                    IZ_FIX,IN_FIX,NUCLID(ISOSPI),
     *                                           RCOULO,IRUNMI,IDEFCN
C
      END IF
C
C=======================================================================
C
      DO ITHEOR=1,LEVTHE
C
         ETHEOR=ENETHE(ITHEOR)
C
         DO IEXPER=1,LEVEXP
C
            IF ((LABEXP(IEXPER).EQ.LABTHE(ITHEOR)).AND.
     *         ((ETHEOR.GE.EMINTH).AND.(ETHEOR.LE.EMAXTH))) THEN
C
                 IF (LOGWRI.GT.4)
     *
     *           WRITE(IRESUL,'(''#   '',I2,'') '',F8.3,2X,F8.3,4X,
     *                             A6,5X,A6,F11.4,2F9.4,''   #'')')
     *
     *                 IEXPER,ETHEOR,EXPEXP(IEXPER),LABTHE(ITHEOR),
     *                               LABEXP(IEXPER),
     *                        ETHEOR-EXPEXP(IEXPER),CHNORM(IEXPER),
     *                                              CHWEIG(IEXPER)
C
                 IF (ISOSPI.EQ.1 .AND. IwMODE.EQ.2 .AND. LOGWRI.GT.4)
     *                 
     *           WRITE(LOGPRO,'(''#   '',I2,'') '',F8.3,2X,F8.3,4X,
     *                             A6,5X,A6,F11.4,2F9.4,''   #'')')
     *
     *                 IEXPER,ETHEOR,EXPEXP(IEXPER),LABTHE(ITHEOR),
     *                               LABEXP(IEXPER),
     *                        ETHEOR-EXPEXP(IEXPER),CHNORM(IEXPER),
     *                                              CHWEIG(IEXPER)
C
                 IF (ISOSPI.EQ.0 .AND. IwMODE.EQ.2 .AND. LOGWRI.GT.4)
     *                 
     *           WRITE(LOGNEU,'(''#   '',I2,'') '',F8.3,2X,F8.3,4X,
     *                             A6,5X,A6,F11.4,2F9.4,''   #'')')
     *
     *                 IEXPER,ETHEOR,EXPEXP(IEXPER),LABTHE(ITHEOR),
     *                               LABEXP(IEXPER),
     *                        ETHEOR-EXPEXP(IEXPER),CHNORM(IEXPER),
     *                                              CHWEIG(IEXPER)
C
            END IF
C
         END DO
C
      END DO
C
C=======================================================================
C
      DO IEXPER=1,LEVEXP
         IF (INDTRU(IEXPER).EQ.0 .AND. LOGWRI.GT.0) THEN
             WRITE(IRESUL,'(80(''!''))')
             WRITE(IRESUL,'(''!!!  Warning : experimental level '',
     *                      ''number '',I3.3,'' has no theoretical '',
     *                      ''counterpart   !!!'')') IEXPER
             WRITE(IRESUL,'(80(''!''))')
         END IF
      END DO
C
C=======================================================================
C     Printing a short table if weighted chi^2 or Max. error or
C     else the correlation matrix test when improved
C=======================================================================
C
      IF (ISOSPI.EQ.0 .AND. LOGWRI.GT.4) THEN
C
          WRITE(IRESUL,'(/,80(''#''),/,''#'',78X,''#'',/,
     *                 ''#   Numerical measures of discrepancies '',
     *                 '' theory/experiment    Z ='',I3,
     *                 ''   N='',I3,''   #'',/,''#   '',A10,
     *                 ''                                        '',
     *                 ''RST='',I2,
     *                 ''  IDEFCN ='',I6,''   #'',/,
     *                 ''#   FermEx='',F7.3,
     *                 ''  FermTh='',F7.3,''   GapEx='',F5.2,
     *                 '' GapTh='',F5.2,'' Re='',F4.2,'' Rt='',F4.2,
     *                 ''   #'',/,
     *                 ''#                                    '',
     *                 ''INVER='',I3,''               Ct='',F7.4,
     *                 ''        #'',/,
     *                 ''#   ChiWei  DifWei  ChiCor  EAbsAv   ErrMx'',
     *                 ''   DenUpT  DenUpE   DenDwT  DenDwE   #'')')
     *
     *                                   IZ_FIX,IN_FIX,NUCLID(ISOSPI),
     *                                                 IRUNMI,IDEFCN,
     *                                   FERMEX(INUCLI),FERTHE,
     *                                   GAPEXP(INUCLI),GAPTHE,
     *                                   RMSEXP,RMSTHE,
     *                                   INVRSN_NEUTRS(IDEFCN),RMODIF
C
          IF (IwMODE.EQ.2 .AND. LOGWRI.GT.4)
     *          
     *    WRITE(LOGNEU,'(/,80(''#''),/,''#'',78X,''#'',/,
     *                 ''#   Numerical measures of discrepancies '',
     *                 '' theory/experiment    Z ='',I3,
     *                 ''   N='',I3,''   #'',/,''#   '',A10,
     *                 ''                                        '',
     *                 ''RST='',I2,
     *                 ''  IDEFCN ='',I6,''   #'',/,
     *                 ''#   FermEx='',F7.3,
     *                 ''  FermTh='',F7.3,''   GapEx='',F5.2,
     *                 '' GapTh='',F5.2,'' Re='',F4.2,'' Rt='',F4.2,
     *                 ''   #'',/,
     *                 ''#                                    '',
     *                 ''INVER='',I3,''               Ct='',F7.4,
     *                 ''        #'',/,
     *                 ''#   ChiWei  DifWei  ChiCor  EAbsAv   ErrMx'',
     *                 ''   DenUpT  DenUpE   DenDwT  DenDwE   #'')')
     *
     *                                   IZ_FIX,IN_FIX,NUCLID(ISOSPI),
     *                                                 IRUNMI,IDEFCN,
     *                                   FERMEX(INUCLI),FERTHE,
     *                                   GAPEXP(INUCLI),GAPTHE,
     *                                   RMSEXP,RMSTHE,
     *                                   INVRSN_NEUTRS(IDEFCN),RMODIF
C
      END IF
C
C=======================================================================
C
      IF (ISOSPI.EQ.1 .AND. LOGWRI.GT.4) THEN
C
          WRITE(IRESUL,'(/,80(''#''),/,''#'',78X,''#'',/,
     *                 ''#   Numerical measures of discrepancies '',
     *                 '' theory/experiment    Z ='',I3,
     *                 ''   N='',I3,''   #'',/,
     *                 ''#   '',A10,''  {Coulomb Eff. Radius ='',
     *                 F6.3,''*A^(1/3)} RST='',I2,''  IDEFCN ='',
     *                 I6,''   #'',/,
     *                 ''#   FermEx='',F7.3,
     *                 ''  FermTh='',F7.3,''   GapEx='',F5.2,
     *                 '' GapTh='',F5.2,'' Re='',F4.2,'' Rt='',F4.2,
     *                 ''   #'',/,
     *                 ''#                                    '',
     *                 ''INVER='',I3,''               Ct='',F7.4,
     *                 ''        #'',/,
     *                 ''#   ChiWei  DifWei  ChiCor  EAbsAv   ErrMx'',
     *                 ''   DenUpT  DenUpE   DenDwT  DenDwE   #'')')
     *
     *                                   IZ_FIX,IN_FIX,NUCLID(ISOSPI),
     *                                          RCOULO,IRUNMI,IDEFCN,
     *                                   FERMEX(INUCLI),FERTHE,
     *                                   GAPEXP(INUCLI),GAPTHE,
     *                                   RMSEXP,RMSTHE,
     *                                   INVRSN_PROTON(IDEFCN),RMODIF
C
          IF (IwMODE.EQ.2 .AND. LOGWRI.GT.4)
     *          
     *    WRITE(LOGPRO,'(/,80(''#''),/,''#'',78X,''#'',/,
     *                 ''#   Numerical measures of discrepancies '',
     *                 '' theory/experiment    Z ='',I3,
     *                 ''   N='',I3,''   #'',/,
     *                 ''#   '',A10,''  {Coulomb Eff. Radius ='',
     *                 F6.3,''*A^(1/3)} RST='',I2,''  IDEFCN ='',
     *                 I6,''   #'',/,
     *                 ''#   FermEx='',F7.3,
     *                 ''  FermTh='',F7.3,''   GapEx='',F5.2,
     *                 '' GapTh='',F5.2,'' Re='',F4.2,'' Rt='',F4.2,
     *                 ''   #'',/,
     *                 ''#                                    '',
     *                 ''INVER='',I3,''               Ct='',F7.4,
     *                 ''        #'',/,
     *                 ''#   ChiWei  DifWei  ChiCor  EAbsAv   ErrMx'',
     *                 ''   DenUpT  DenUpE   DenDwT  DenDwE   #'')')
     *
     *                                   IZ_FIX,IN_FIX,NUCLID(ISOSPI),
     *                                          RCOULO,IRUNMI,IDEFCN,
     *                                   FERMEX(INUCLI),FERTHE,
     *                                   GAPEXP(INUCLI),GAPTHE,
     *                                   RMSEXP,RMSTHE,
     *                                   INVRSN_PROTON(IDEFCN),RMODIF
C
      END IF
C
C=======================================================================
C
      IF (LOGWRI.GT.4) 
     *    WRITE(IRESUL,'(''# '',4F8.4,F8.3,F9.4,F8.4,F9.4,F8.4,
     *      ''   #'')')CHIWEI,DIFWEI,CHICOR,EABSAV,ERRMAX,DENUPP,DENSUP,
     *                                                    DENLOW,DENSDW
C
      IF (LOGWRI.GT.4) 
     *    WRITE(IRESUL,'(''#'',T37,A6,'':'',T46,''CHISQU_CHOICE='',
     *             F12.4,T80,''#'',/, 80(''#''))')CHIDEF,CHISQU_CHOICE
C
      IF (ISOSPI.EQ.0 .AND. IwMODE.EQ.2 .AND. LOGWRI.GT.4) THEN
C      
          WRITE(LOGNEU,'(''# '',4F8.4,F8.3,F9.4,F8.4,F9.4,F8.4,
     *                                                       ''   #'')')
     *          CHIWEI,DIFWEI,CHICOR,EABSAV,ERRMAX,DENUPP,DENSUP,
     *                                             DENLOW,DENSDW
C
          WRITE(LOGNEU,'(''#'',T37,A6,'':'',T46,''CHISQU_CHOICE='',
     *                   F12.4,T80,''#'',/, 80(''#''))')
     *                                            CHIDEF,CHISQU_CHOICE
      END IF
C
      IF (ISOSPI.EQ.1 .AND. IwMODE.EQ.2 .AND. LOGWRI.GT.4) THEN
C      
          WRITE(LOGPRO,'(''# '',4F8.4,F8.3,F9.4,F8.4,F9.4,F8.4,
     *                                                       ''   #'')')
     *          CHIWEI,DIFWEI,CHICOR,EABSAV,ERRMAX,DENUPP,DENSUP,
     *                                             DENLOW,DENSDW
C
          WRITE(LOGPRO,'(''#'',T37,A6,'':'',T46,''*CHISQU_CHOICE='',
     *                  F12.4,T80,''#'',/, 80(''#''))')
     *                                            CHIDEF,CHISQU_CHOICE
C
      END IF
C
C=======================================================================    
C
      CALL WEIPRI(IRESUL)
C
C=======================================================================
C
      RETURN
      END
C
C=======================================================================
C=======================================================================
C
      SUBROUTINE WARNIN(N_UNIT,WHATEX)
      CHARACTER
     *          WHATEX*6
C
C=======================================================================
C     This is a WARNING PRINT routine only
C=======================================================================
C
      IF (WHATEX.EQ.'EXTRUE'.or.WHATEX.EQ.'EXPNEW') RETURN
C
      WRITE(N_UNIT,'(''  N   N   OOO  TTTTT    EEEEE  X   X  PPPP '',
     *               ''  EEEEE  RRRR   III  M   M  !!!!!!'',/,
     *               ''  NN  N  O   O   T      E       X X   P   P'',
     *               ''  E      R   R   I   MM MM   !!!! '',/,
     *               ''  N N N  O   O   T      EEE      X    PPPP '',
     *               ''  EEE    RRRR    I   M M M    !!  '',/,
     *               ''  N  NN  O   O   T      E       X X   P    '',
     *               ''  E      R  R    I   M   M        '',/,
     *               ''  N   N   OOO    T      EEEEE  X   X  P    '',
     *               ''  EEEEE  R   R  III  M   M    !!  '',//,
     *               80(''=''),/,80(''=''),/)')            
C
C=======================================================================
C
      RETURN
      END      
C
C=======================================================================
C=======================================================================
C
      SUBROUTINE WEIPRI(N_UNIT)
C
      COMMON
     *       /WEIGHT/ WEIGHT_CORREL,WEIGHT_RADIUS,WEIGHT_INVERT,
     *                WEIGHT_EFERMI,WEIGHT_ENEGAP,WEIWEI,
     *                WEIGHT_ERRABS,WEIGHT_EABSAV,WEIGHT_ERRMAX,
     *                WEIGHT_DENSUP,WEIGHT_DENSDW,WEIGHT_RHODEN
      COMMON
     *       /CHIVAL/ CHISQU_CORREL,DIFSQU_RADIUS,CHISQU_INVERT,
     *                DIFSQU_EFERMI,DIFSQU_ENEGAP,CHIWEI_ENERGY,
     *                CHIWEI_ENEDEG,ERRABS_WEIDEG,EABSAV,ERRMAX,
     *                       DIFSQU_DENSUP,DIFSQU_DENSDW,CHIRHO
      COMMON
     *       /DECOMP/ VEICOR,VEIRAD,VEIINV,VEIFER,VEIGAP,VEIWEI,
     *                       VEIDIF,VEIABS,VEIMAX,VEIDUP,VEIDDW
      COMMON
     *         /PRFILE/ ISCREN,LOGFIL,LOGCHP,LOGCHN,LOGHIP,LOGHIN,
     *                  LOGPRO,LOGNEU,LOGWRI,LSCREN,IRESUL,ICONVE,
     *                  NOUTPT,LOGAUX,LOGENE,LOGRAD,LOGGRD,LOGMSH,
     *                                                     LOGBIS  
C
C=======================================================================
C     This subroutine ???
C=======================================================================
C
      IF (LOGWRI.GT.4) THEN
          WRITE(N_UNIT,'(/,80(''#''),/,''#'',78X,''#'',/,
     *                 ''# WEICOR WEIRAD WEIINV WEIFER WEIGAP WEIWEI'',
     *                 '' WEIDIF WEIABS WEIMAX WEIDUP WEIDDW #'',/,
     *                 ''#'',11(1X,F6.2),'' #'')')
     *
     *                 WEIGHT_CORREL,WEIGHT_RADIUS,WEIGHT_INVERT,
     *                 WEIGHT_EFERMI,WEIGHT_ENEGAP,WEIWEI,
     *                 WEIGHT_ERRABS,WEIGHT_EABSAV,WEIGHT_ERRMAX,
     *                 WEIGHT_DENSUP,WEIGHT_DENSDW
C
          WRITE(N_UNIT,'(  ''# CHICOR DIFRAD CHIINV FERDIF GAPDIF'',
     *                 '' CHIWEI DIFWEI EABSAV ERRMAX DIFFUP DIFFDW #'',
     *                 /,''#'',1X,2F6.3,1X,9(F7.1),'' #'')')
     *
     *                 CHISQU_CORREL,DIFSQU_RADIUS,CHISQU_INVERT,
     *                 DIFSQU_EFERMI,DIFSQU_ENEGAP,CHIWEI_ENERGY,
     *                 ERRABS_WEIDEG,EABSAV,ERRMAX,DIFSQU_DENSUP,
     *                                             DIFSQU_DENSDW 
C
      END IF
C
      FEICOR=VEICOR*WEIGHT_CORREL*CHISQU_CORREL
      FEIRAD=VEIRAD*WEIGHT_RADIUS*DIFSQU_RADIUS
      FEIINV=VEIINV*WEIGHT_INVERT*CHISQU_INVERT
      FEIFER=VEIFER*WEIGHT_EFERMI*DIFSQU_EFERMI
      FEIGAP=VEIGAP*WEIGHT_ENEGAP*DIFSQU_ENEGAP
      FEIWEI=VEIWEI*WEIWEI*CHIWEI_ENERGY
      FEIDIF=VEIDIF*WEIGHT_ERRABS*ERRABS_WEIDEG
      FEIABS=VEIABS*WEIGHT_EABSAV*EABSAV
      FEIMAX=VEIMAX*WEIGHT_ERRMAX*ERRMAX
      FEIDUP=VEIDUP*WEIGHT_DENSUP*DIFSQU_DENSUP
      FEIDDW=VEIDDW*WEIGHT_DENSDW*DIFSQU_DENSDW
C
      IF (LOGWRI.GT.4) THEN
          WRITE(N_UNIT,'( ''# FEICOR FEIRAD FEIINV FEIFER FEIGAP'',
     *                 '' FEIWEI FEIDIF FEIABS FEIMAX FEIDUP FEIDDW #'',
     *                 /,''#'',11(F7.2),'' #'',/,
     *                 ''#'',78X,''#'',/,80(''#''))')
     *
     *                    FEICOR,FEIRAD,FEIINV,FEIFER,FEIGAP,FEIWEI,
     *           	         FEIDIF,FEIABS,FEIMAX,FEIDUP,FEIDDW  
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
      SUBROUTINE STRLNG(STRING,LENGTH,STRAUX)
C      
      INCLUDE  'MATDIM/NDLNGT.f'
      INCLUDE  'MATDIM/NDX256.f'
C      
      CHARACTER
     *          STRING*256,BLANCK*1
      CHARACTER
     *          STRAUX*256
      DATA
     *     BLANCK / ' ' /
C
C=======================================================================
C     This subroutine gives the length of a string "STRING"
C     counting all non-blank characters
C=======================================================================
C
      LENGTH=0
C
      DO I=1,NDLNGT
         IF (STRING(I:I).NE.BLANCK) THEN
             LENGTH=LENGTH+1
         END IF
      END DO
C
      DO I=1,NDX256
         STRAUX(I:I)='.'
      END DO
C
      DO I=1,MIN(LENGTH,NDX256)
         STRAUX(I:I)=STRING(I:I)
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

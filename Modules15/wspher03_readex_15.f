C FILE NAME = wspher03_readex_15.f ! Keep this symbol:    $ident@string$
C=======================================================================
C=======================================================================
C=======================================================================
C      
C=======================================================================
C=======================================================================
C                READING EXPERIMENTAL SINGLE PARTICLE LEVELS
C=======================================================================
C=======================================================================
C
      SUBROUTINE READEX_LEVELS(IFEXPE,IDEFCN)
C
      INCLUDE  'MATDIM/NDLEXP.f'
      INCLUDE  'MATDIM/NDNUEX.f'
C
      CHARACTER
     *          FILEXP*256,SYMBNU_AUXILI*5,KEYWOR*5
      CHARACTER
     *                     SYMBEX*6,SYMBPL*10
      DIMENSION
     *          SYMBNU_AUXILI(1:NDNUEX)
C
      COMMON
     *       /DEPPAI/ IFDEEP,IFPRON
      COMMON
     *       /EXPZ_N/ IZEXPE(1:NDNUEX),INEXPE(1:NDNUEX),
     *                                 NOFEXP(1:NDNUEX)
      COMMON
     *       /EXPLMB/ BINDLV(1:NDNUEX),BINDOR(1:NDNUEX),
     *                                 DELPAI(1:NDNUEX)
      COMMON
     *       /EXPENE/ ENEXPE(1:NDNUEX,1:NDLEXP),
     *                E_OROS(1:NDNUEX,1:NDLEXP),
     *                ESHELM(1:NDNUEX,1:NDLEXP),
     *                ENEWEX(1:NDNUEX,1:NDLEXP)
      COMMON
     *       /EXPIND/ LNEXPE(1:NDNUEX,1:NDLEXP),
     *                L_OROS(1:NDNUEX,1:NDLEXP),
     *                LSHELM(1:NDNUEX,1:NDLEXP)
      COMMON
     *       /EXPSLV/ SYMBEX(1:NDNUEX,1:NDLEXP),
     *                SYMBPL(1:NDNUEX,1:NDLEXP)
      COMMON
     *       /EXPNUC/ NUCEXP
C
      COMMON
     *       /PRFILE/ ISCREN,LOGFIL,LOGCHP,LOGCHN,LOGHIP,LOGHIN,
     *                LOGPRO,LOGNEU,LOGWRI,LSCREN,IRESUL,ICONVE,
     *                NOUTPT,LOGAUX,LOGENE,LOGRAD,LOGGRD,LOGMSH,
     *                                                   LOGBIS
C_______________________________________________________________________
C
      DATA
     *     N_UNIT / 10 /
C
C=======================================================================
C     Subroutine reading the experimental data in NAMELIST style
C=======================================================================
C
      IF (IDEFCN.GT.1) THEN
          IF (LOGWRI.GT.4) THEN
              WRITE(LOGFIL,'(12X,''Exiting  READEX_LEVELS, no action''
     *                                                             )')
          END IF
          RETURN
      END IF
C
C     IDEFCN above =>> iteration count in minimising routines
C
C=======================================================================
C
      IF (IFDEEP.EQ.0 .AND. IFPRON.EQ.0) THEN
      
          FILEXP='ws_exper/exp_lev_sph.d'
C
          OPEN(N_UNIT,FILE=FILEXP,STATUS='UNKNOWN',FORM='FORMATTED')
C
          IF (LOGWRI.GT.4) THEN
              WRITE(LOGFIL,'(12X,''Opening Unit='',I2,1X,
     *                           ''File=ws_exper/exp_lev_sph.d'')')
     *                                                  N_UNIT
          END IF    
      
      END IF
C
C=======================================================================
C
      IF (IFDEEP.EQ.1 .AND. IFPRON.EQ.0) THEN
      
          FILEXP='ws_exper/exp_lev_sph_ifdeep.d'
C
          OPEN(N_UNIT,FILE=FILEXP,STATUS='UNKNOWN',FORM='FORMATTED')
C
          IF (LOGWRI.GT.4) THEN
              WRITE(LOGFIL,'(12X,''Opening Unit='',I2,1X,
     *                         ''File=ws_exper/exp_lev_sph_ifdeep.d''
     *                                                            )') 
     *                                                      N_UNIT
          END IF    
      
      END IF
C
C=======================================================================
C
      IF (IFDEEP.EQ.0 .AND. IFPRON.EQ.1) THEN
      
          FILEXP='ws_exper/exp_lev_sph_ifpron.d'
C
          OPEN(N_UNIT,FILE=FILEXP,STATUS='UNKNOWN',FORM='FORMATTED')
C
          IF (LOGWRI.GT.4) THEN
              WRITE(LOGFIL,'(12X,''Opening Unit='',I2,1X,
     *                        ''File=ws_exper/exp_lev_sph_ifpron.d'')') 
     *                                                      N_UNIT
          END IF    
      
      END IF
C
C=======================================================================
C      
      IF (IFDEEP.EQ.1 .AND. IFPRON.EQ.1) THEN
      
          FILEXP='ws_exper/exp_lev_sph_ifdeep_ifpron.d'
C
          OPEN(N_UNIT,FILE=FILEXP,STATUS='UNKNOWN',FORM='FORMATTED')
C
          IF (LOGWRI.GT.4) THEN
              WRITE(LOGFIL,'(12X,''Opening Unit='',I2,1X,
     *                ''File=ws_exper/exp_lev_sph_ifdeep_ifpron.d'')') 
     *                                                      N_UNIT
          END IF    
      
      END IF
C
C=======================================================================
C
      DO NO_EXP=1,NDNUEX
         DO IEXPER=1,NDLEXP
            SYMBEX(NO_EXP,IEXPER)='      '
            SYMBPL(NO_EXP,IEXPER)='yyyyyy'
         END DO
      END DO
C
C=======================================================================
C
      NO_EXP=0
C
   1  CONTINUE
C
C=======================================================================
C
      IF (NO_EXP.GT.NDNUEX) THEN
C
          WRITE(LSCREN,'(//,''NO_EXP='',I3,'' exceeds NDNUEX='',I2)')
     *                        NO_EXP,                 NDNUEX
          STOP 'NO_EXP.GT.NDNUEX in READEX_LEVELS'
C
      END IF
C
C=======================================================================
C
      KEYWOR='     '
C
      READ (N_UNIT,'(A5)',END=2) KEYWOR
C
C=======================================================================
C
      IF (KEYWOR.EQ.'ENDGO') THEN
          NUCEXP=NO_EXP
          CLOSE(N_UNIT)
          IF (LOGWRI.GT.4) THEN
              WRITE(LOGFIL,'(12X,''Closing UNIT='',I2)') N_UNIT
              WRITE(LOGFIL,'( 9X,''Exiting  READEX_LEVELS'')')
          END IF
          RETURN
      END IF
C
C=======================================================================
C
      IF (KEYWOR(2:3).NE.'  ') THEN
C
          NO_EXP=NO_EXP+1
C
          SYMBNU_AUXILI(NO_EXP)=KEYWOR
C
          READ(N_UNIT,*) IZEXPE(NO_EXP),INEXPE(NO_EXP),
     *                   BINDOR(NO_EXP),BINDLV(NO_EXP),
     *                                  DELPAI(NO_EXP)
C_______________________________________________________________________
C
          WRITE(NOUTPT,'(/,''Z='',I2,'' N='',I3,/,T10,''Bexp='',f7.3,
     *                     ''   DPair='',F6.4,/)') 
     *
     *                       IZEXPE(NO_EXP),INEXPE(NO_EXP),
     *                                      BINDLV(NO_EXP),
     *                                      DELPAI(NO_EXP)
C_______________________________________________________________________
C
          IF (LOGWRI.GT.5)
     *          
     *    WRITE(LOGFIL,'(/,9X,''Z='',I2,'' N='',I3,/,T10,''Bexp='',f7.3,
     *                                           ''   DPair='',F6.4)') 
     *
     *                          IZEXPE(NO_EXP),INEXPE(NO_EXP),
     *                                         BINDLV(NO_EXP),
     *                                         DELPAI(NO_EXP)
C_______________________________________________________________________
C
          WRITE(NOUTPT,'(T16,''Oros    Eexp   E_shm   E_new'')')
C
          IF (LOGWRI.GT.5)
     *          
     *    WRITE(LOGFIL,'(T16,''Oros    Eexp   E_shm   E_new'')')
C_______________________________________________________________________
C
          DO IEXPER=1,NDLEXP
C
             READ(N_UNIT,'(20X,A6,43X,A10)') SYMBEX(NO_EXP,IEXPER),
     *                                       SYMBPL(NO_EXP,IEXPER)
C
             IF (SYMBEX(NO_EXP,IEXPER).EQ.'xxxxxx') THEN
                 NOFEXP(NO_EXP)=IEXPER-1
                 GO TO 1
             ELSE
                 READ(N_UNIT,*) E_OROS(NO_EXP,IEXPER),
     *                          L_OROS(NO_EXP,IEXPER),
     *                          ENEXPE(NO_EXP,IEXPER),
     *                          LNEXPE(NO_EXP,IEXPER),
     *                          ESHELM(NO_EXP,IEXPER),
     *                          ENEWEX(NO_EXP,IEXPER)
             END IF
C
             WRITE(NOUTPT,'(''#   '',i1,'')'',4f8.2,2x,a6)')
     *                                                 IEXPER,
     *             E_OROS(NO_EXP,IEXPER),ENEXPE(NO_EXP,IEXPER),
     *             ESHELM(NO_EXP,IEXPER),ENEWEX(NO_EXP,IEXPER),
     *                                   SYMBEX(NO_EXP,IEXPER)
C
          IF (LOGWRI.GT.5)
     *          
     *       WRITE(LOGFIL,'(T10,i1,'')'',4f8.2,2x,a6)')
     *                                                 IEXPER,
     *             E_OROS(NO_EXP,IEXPER),ENEXPE(NO_EXP,IEXPER),
     *             ESHELM(NO_EXP,IEXPER),ENEWEX(NO_EXP,IEXPER),
     *                                   SYMBEX(NO_EXP,IEXPER)
C
          END DO
C
      END IF
C
C=======================================================================
C
      GO TO 1
C
C=======================================================================
C
   2  CONTINUE
C
      WRITE(LSCREN,'(//,''Reading experimental data - '',
     *             ''END OF FILE encountered'',//,
     *             ''This file should end with  ENDGO  statement'')')
C
C=======================================================================
C
      STOP 'ENDGO statement missing READEX_LEVELS'
      END      
C
C=======================================================================
C=======================================================================
C
      SUBROUTINE PREPAR_EXPLEV(IZ_FIX,IN_FIX,FERMEX,GAPEXP,DENSUP,
     *                                                     DENSDW,
     *                         IDEFCN,WHATEX,ISOSPI,IFEXPE,INUCLI)
C      
      INCLUDE   'MATDIM/NDSPEC.f'
      INCLUDE   'MATDIM/NDLEXP.f'
      INCLUDE   'MATDIM/NDNUEX.f'
      INCLUDE   'MATDIM/NDRAUS.f'
      INCLUDE   'MATDIM/NDNUCL.f'
C_______________________________________________________________________
C
      CHARACTER
     *          SYMBPL*10,SYMBEX*6,LABEXP*6,LABAUX*6,WHATEX*06
      CHARACTER
     *          LABEXP_PROTON*6,LABEXP_NEUTRS*6
      CHARACTER
     *          LABPRO_REMOVE*6,LABNEU_REMOVE*6,
     *          REMOVE_PROTON*3,REMOVE_NEUTRS*3
C_______________________________________________________________________
C
      DIMENSION
     *          ENEAUX(1:NDSPEC),INDAUX(1:NDSPEC),
     *          LABAUX(1:NDSPEC),IDEAUX(1:NDSPEC)
      DIMENSION
     *          EXPEXP(1:NDLEXP),IDEGEX(1:NDLEXP)
      DIMENSION
     *          LABEXP(1:NDLEXP)
C_______________________________________________________________________
C
      COMMON
     *       /REMLAB/ REMOVE_PROTON,REMOVE_NEUTRS
      COMMON
     *       /OUTLAB/ LABPRO_REMOVE(1:NDRAUS),
     *                LABNEU_REMOVE(1:NDRAUS)
      COMMON
     *       /OUTIND/ INDPRO_LETOUT(1:NDRAUS),
     *                INDNEU_LETOUT(1:NDRAUS)
C
C=======================================================================
C     Follows the full information about experimental spherical
C                                                   l e v e l s
C=======================================================================
C
      COMMON
     *       /EXPZ_N/ IZEXPE(1:NDNUEX),INEXPE(1:NDNUEX),
     *                                 NOFEXP(1:NDNUEX)
      COMMON
     *       /EXPLMB/ BINDLV(1:NDNUEX),BINDOR(1:NDNUEX),
     *                                 DELPAI(1:NDNUEX)
      COMMON
     *       /EXPENE/ ENEXPE(1:NDNUEX,1:NDLEXP),
     *                E_OROS(1:NDNUEX,1:NDLEXP),
     *                ESHELM(1:NDNUEX,1:NDLEXP),
     *                ENEWEX(1:NDNUEX,1:NDLEXP)
      COMMON
     *       /EXPSLV/ SYMBEX(1:NDNUEX,1:NDLEXP),
     *                SYMBPL(1:NDNUEX,1:NDLEXP)
      COMMON
     *       /EXPNUC/ NUCEXP
C                               <<<=== HERE: Final experimental
C                                      energies and level SYMBS
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
     *       /DIMSPH/ LDSING_PROTON,LEVTHE_PROTON,
     *                LEVEXP_PROTON(1:NDNUCL),
     *                LDSING_NEUTRS,LEVTHE_NEUTRS,
     *                LEVEXP_NEUTRS(1:NDNUCL)
      COMMON
     *       /RANRAN/ IRMFST,I_RAND,IRANDO
      COMMON
     *       /WHERBE/ IRUNMI,NEWOLD 
      COMMON
     *       /PRFILE/ ISCREN,LOGFIL,LOGCHP,LOGCHN,LOGHIP,LOGHIN,
     *                LOGPRO,LOGNEU,LOGWRI,LSCREN,IRESUL,ICONVE,
     *                NOUTPT,LOGAUX,LOGENE,LOGRAD,LOGGRD,LOGMSH,
     *                                                   LOGBIS
C
C=======================================================================
C
C     This subroutine provides experimental spectra for magic spherical
C     nuclei (based on the data-bank information) read by READEX_LEVELS
C
C     More precisely, binding energies are transformed into s.p. levels
C
C=======================================================================
C
C     Meaning of some parameters:
C
C     FERMEX - experimental estimate  of  the  Fermi level
C
C     GAPEXP - experimental estimate of the main shell gap
C
C     DENSUP - experimental estimate of the level  density
C                      a b o v e    the main shell closure
C     DENSDW - experimental estimate of the level  density
C                      b e l o w    the main shell closure
C
C=======================================================================
C 
C         IDEFCN=1 ->> for the frist iteration in MINIMI
C 
      IF (IDEFCN.GT.1) THEN
          IF (LOGWRI.GT.4) THEN
              WRITE(LOGFIL,'(9X,''Exiting  PREPAR_EXPLEV, no action'')')
          END IF
          RETURN
      END IF
C
C=======================================================================
C
      NEXIST=1
C
C=======================================================================
C
      IF (MOD(IZ_FIX,2).NE.0) THEN
C
          WRITE(LSCREN,'(/,80(''|''),/,''|'',78X,''|'',/,
     *                   ''|   Comparison with experiment for odd Z='',
     *                     I3,
     *                   '' impossible with the existing data  |'',/,
     *                   ''|'',78X,''|'',/,80(''|''),/)') IZ_FIX
C
          IF (LOGWRI.GT.5)
     *
     *    WRITE(LOGFIL,'(/,80(''|''),/,''|'',78X,''|'',/,
     *                   ''|   Comparison with experiment for odd Z='',
     *                     I3,
     *                   '' impossible with the existing data  |'',/,
     *                   ''|'',78X,''|'',/,80(''|''),/)') IZ_FIX
C
          STOP 'Comparison with experiment impossible (Z-case)'
C
      END IF
C
C=======================================================================
C
      IF (MOD(IN_FIX,2).NE.0) THEN
C
          WRITE(LSCREN,'(/,80(''|''),/,''|'',78X,''|'',/,
     *                   ''|   Comparison with experiment for odd N='',
     *                     I3,
     *                   '' impossible with the existing data  |'',/,
     *                   ''|'',78X,''|'',/,80(''|''),/)') IN_FIX
C
          IF (LOGWRI.GT.5)
     *
     *    WRITE(LOGFIL,'(/,80(''|''),/,''|'',78X,''|'',/,
     *                   ''|   Comparison with experiment for odd N='',
     *                     I3,
     *                   '' impossible with the existing data  |'',/,
     *                   ''|'',78X,''|'',/,80(''|''),/)') IN_FIX
C
          STOP 'Comparison with experiment impossible (N-case)'
C
      END IF
C                
C=======================================================================
C
      DO LEVELS=1,NDLEXP
         EXPEXP(LEVELS)=999.999
      END DO
C
      LEVEXP=0
      LEVEXD=0
C
      DO NUCLEU=1,NUCEXP
C
C=======================================================================
C        Extracting the energies of the PROTON hole-states
C=======================================================================
C
         IF (IZEXPE(NUCLEU).EQ.IZ_FIX-1 .AND. ISOSPI.EQ.1
     *                                  .AND.
     *       INEXPE(NUCLEU).EQ.IN_FIX)               THEN
C
             DO LEVELS=1,NOFEXP(NUCLEU)
C
C=======================================================================
C               True experimental levels
C=======================================================================
C
                IF (WHATEX.EQ.'EXTRUE'.AND.
     *                         ENEXPE(NUCLEU,LEVELS).GE.0.000) THEN
C
                    LEVEXP=LEVEXP+1
                    EXPEXP(LEVEXP)=BINDLV(NUCLEU)
     *                            -ENEXPE(NUCLEU,LEVELS)
C @@@ IRENE, What are EXPLOW??? LEAVE A COMMENT
                    IF (LEVEXP.EQ.1) EXPLOW=EXPEXP(LEVEXP)
C
                END IF
C
C=======================================================================
C               Experiment plus corrections from A. Oros
C=======================================================================
C
                IF (WHATEX.EQ.'EXOROS'.AND.
     *                         E_OROS(NUCLEU,LEVELS).GE.0.000) THEN
C
                    LEVEXP=LEVEXP+1
                    EXPEXP(LEVEXP)=BINDOR(NUCLEU)
     *                            -E_OROS(NUCLEU,LEVELS)
C
                    IF (LEVEXP.EQ.1) EXPLOW=EXPEXP(LEVEXP)
C
                END IF
C
C=======================================================================
C               Using the shell model instead of experiment
C=======================================================================
C
                IF (WHATEX.EQ.'SHELLM'.AND.
     *                         ESHELM(NUCLEU,LEVELS).GE.0.000) THEN
C
                    LEVEXP=LEVEXP+1
                    EXPEXP(LEVEXP)=BINDLV(NUCLEU)
     *                            -ESHELM(NUCLEU,LEVELS)
C
                    IF (LEVEXP.EQ.1) EXPLOW=EXPEXP(LEVEXP)
C
                END IF
C  @@@ IRENE - PLEASE CONFIRM: PORQUET?
C=======================================================================
C               NEW experimental levels MG PORQUET ??
C=======================================================================
C
                IF (WHATEX.EQ.'EXPNEW'.AND.
     *                         ENEWEX(NUCLEU,LEVELS).NE.-.99)  THEN
C
                    LEVEXP=LEVEXP+1
                    EXPEXP(LEVEXP)=ENEWEX(NUCLEU,LEVELS)
C
                    IF (LEVEXP.EQ.1) EXPLOW=EXPEXP(LEVEXP)

                END IF
C
C=======================================================================
C               Defining the spectroscopic labels
C=======================================================================
C
                IF (LEVEXD.NE.LEVEXP) THEN
C
                    LEVEXD=LEVEXP
                    LABEXP(LEVEXP)=SYMBEX(NUCLEU,LEVELS)
C
                END IF
C
C=======================================================================
C
             END DO
C
             NEXIST=0
C
         END IF
C
C=======================================================================
C
         LVPRIM=0
C
C=======================================================================
C        Extracting the energies of the PROTON particle states
C=======================================================================
C
         IF (IZEXPE(NUCLEU).EQ.IZ_FIX+1.AND.ISOSPI.EQ.1
     *                                 .AND.
     *       INEXPE(NUCLEU).EQ.IN_FIX) THEN
C
             DO LEVELS=1,NOFEXP(NUCLEU)
C
C=======================================================================
C               True experimental levels
C=======================================================================
C
                IF (WHATEX.EQ.'EXTRUE'.AND.
     *                         ENEXPE(NUCLEU,LEVELS).GE.0.000) THEN
C
                    LEVEXP=LEVEXP+1
                    EXPEXP(LEVEXP)=BINDLV(NUCLEU)
     *                            +ENEXPE(NUCLEU,LEVELS)
                    LVPRIM=LVPRIM+1
C
                    IF (LVPRIM.EQ.1) EXPUPP=EXPEXP(LEVEXP)
C
                END IF
C
C=======================================================================
C               Experiment plus corrections from A. Oros
C=======================================================================
C
C @@@ IRENE, WHAT ARE THESE "CORRECTIONS" SUPPOSED OT BE ???
C
                IF (WHATEX.EQ.'EXOROS'.AND.
     *                         E_OROS(NUCLEU,LEVELS).GE.0.000) THEN
C
                    LEVEXP=LEVEXP+1
                    EXPEXP(LEVEXP)=BINDOR(NUCLEU)
     *                            +E_OROS(NUCLEU,LEVELS)
                    LVPRIM=LVPRIM+1
C
                    IF (LVPRIM.EQ.1) EXPUPP=EXPEXP(LEVEXP)
C
                END IF
C
C=======================================================================
C               Using the shell model instead of experiment
C=======================================================================
C
                IF (WHATEX.EQ.'SHELLM'.AND.
     *                         ESHELM(NUCLEU,LEVELS).GE.0.000) THEN
C
                    LEVEXP=LEVEXP+1
                    EXPEXP(LEVEXP)=BINDLV(NUCLEU)
     *                            +ESHELM(NUCLEU,LEVELS)
                    LVPRIM=LVPRIM+1
C
                    IF (LVPRIM.EQ.1) EXPUPP=EXPEXP(LEVEXP)
C
                END IF
C
C=======================================================================
C               NEW experimental levels
C=======================================================================
C
C @@@ IRENE - WHAT ARE "NEW" LEVELS ???
C
                IF (WHATEX.EQ.'EXPNEW'.AND.
     *                         ENEWEX(NUCLEU,LEVELS).NE.-.99) THEN
C
                    LEVEXP=LEVEXP+1
                    EXPEXP(LEVEXP)=ENEWEX(NUCLEU,LEVELS)
C
                    LVPRIM=LVPRIM+1
C
                    IF (LVPRIM.EQ.1) EXPUPP=EXPEXP(LEVEXP)

                END IF
C
C=======================================================================
C
                IF (LEVEXD.NE.LEVEXP) THEN
C
                    LEVEXD=LEVEXP
                    LABEXP(LEVEXP)=SYMBEX(NUCLEU,LEVELS)
C
                END IF
C
C=======================================================================
C
             END DO
C
             NEXIST=0
C
         END IF
C
C=======================================================================
C        Extracting the energies of the NEUTRON hole states
C=======================================================================
C
         IF (INEXPE(NUCLEU).EQ.IN_FIX-1 .AND. ISOSPI.EQ.0
     *                                  .AND.
     *       IZEXPE(NUCLEU).EQ.IZ_FIX)               THEN
C
             DO LEVELS=1,NOFEXP(NUCLEU)
C
C=======================================================================
C               Here true experimental data
C=======================================================================
C
                IF (WHATEX.EQ.'EXTRUE'.AND.
     *                         ENEXPE(NUCLEU,LEVELS).GE.0.0) THEN
C
                    LEVEXP=LEVEXP+1
                    EXPEXP(LEVEXP)=BINDLV(NUCLEU)
     *                            -ENEXPE(NUCLEU,LEVELS)
C
                    IF (LEVEXP.EQ.1) EXPLOW=EXPEXP(LEVEXP)
C
                END IF
C
C=======================================================================
C               Here corrections from A. Oros
C=======================================================================
C
                IF (WHATEX.EQ.'EXOROS'.AND.
     *                         E_OROS(NUCLEU,LEVELS).GE.0.0) THEN
C
                    LEVEXP=LEVEXP+1
                    EXPEXP(LEVEXP)=BINDOR(NUCLEU)
     *                            -E_OROS(NUCLEU,LEVELS)
C
                    IF (LEVEXP.EQ.1) EXPLOW=EXPEXP(LEVEXP)
C
                END IF
C
C=======================================================================
C               Here shell model instead of experiment
C=======================================================================
C
                IF (WHATEX.EQ.'SHELLM'.AND.
     *                         ESHELM(NUCLEU,LEVELS).GE.0.0) THEN
C
                    LEVEXP=LEVEXP+1
                    EXPEXP(LEVEXP)=BINDLV(NUCLEU)
     *                            -ESHELM(NUCLEU,LEVELS)
C
                    IF (LEVEXP.EQ.1) EXPLOW=EXPEXP(LEVEXP)
C
                END IF
C
C=======================================================================
C               NEW experimental levels
C=======================================================================
C
                IF (WHATEX.EQ.'EXPNEW'.AND.
     *                         ENEWEX(NUCLEU,LEVELS).NE.-.99)THEN
C
                    LEVEXP=LEVEXP+1
                    EXPEXP(LEVEXP)=ENEWEX(NUCLEU,LEVELS)
C
                   IF (LEVEXP.EQ.1) EXPLOW=EXPEXP(LEVEXP)

                END IF
C
C=======================================================================
C
                IF (LEVEXD.NE.LEVEXP) THEN
C
                    LEVEXD=LEVEXP
                    LABEXP(LEVEXP)=SYMBEX(NUCLEU,LEVELS)
C
                END IF
C
C=======================================================================
C
             END DO
C
             NEXIST=0
C
         END IF
C
C=======================================================================
C
         LVPRIM=0
C
C=======================================================================
C        Extracting the energies of the NEUTRON particle states
C=======================================================================
C
         IF (INEXPE(NUCLEU).EQ.IN_FIX+1.AND.ISOSPI.EQ.0
     *                                 .AND.
     *       IZEXPE(NUCLEU).EQ.IZ_FIX)             THEN
C
             DO LEVELS=1,NOFEXP(NUCLEU)
C
C=======================================================================
C               Here true experimental data
C=======================================================================
C
                IF (WHATEX.EQ.'EXTRUE'.AND.
     *                         ENEXPE(NUCLEU,LEVELS).GE.0.0) THEN
C
                    LEVEXP=LEVEXP+1
                    EXPEXP(LEVEXP)=BINDLV(NUCLEU)
     *                            +ENEXPE(NUCLEU,LEVELS)
C
                    LVPRIM=LVPRIM+1
C
                    IF (LVPRIM.EQ.1) EXPUPP=EXPEXP(LEVEXP)
C
                END IF
C
C=======================================================================
C               Here corrections from A. Oros
C=======================================================================
C
                IF (WHATEX.EQ.'EXOROS'.AND.
     *                         E_OROS(NUCLEU,LEVELS).GE.0.0) THEN
C
                    LEVEXP=LEVEXP+1
                    EXPEXP(LEVEXP)=BINDOR(NUCLEU)
     *                            +E_OROS(NUCLEU,LEVELS)
C
                    LVPRIM=LVPRIM+1
C
                    IF (LVPRIM.EQ.1) EXPUPP=EXPEXP(LEVEXP)
C
                END IF
C
C=======================================================================
C               Here shell model instead of experiment
C=======================================================================
C
                IF (WHATEX.EQ.'SHELLM'.AND.
     *                         ESHELM(NUCLEU,LEVELS).GE.0.0) THEN
C
                    LEVEXP=LEVEXP+1
                    EXPEXP(LEVEXP)=BINDLV(NUCLEU)
     *                            +ESHELM(NUCLEU,LEVELS)
C
                    LVPRIM=LVPRIM+1
C
                    IF (LVPRIM.EQ.1) EXPUPP=EXPEXP(LEVEXP)
C
                END IF
C
C=======================================================================
C               NEW experimental levels
C=======================================================================
C
                IF (WHATEX.EQ.'EXPNEW'.AND.
     *                         ENEWEX(NUCLEU,LEVELS).NE.-.99) THEN
C
                    LEVEXP=LEVEXP+1
                    EXPEXP(LEVEXP)=ENEWEX(NUCLEU,LEVELS)
C
                    LVPRIM=LVPRIM+1
C
                    IF (LVPRIM.EQ.1) EXPUPP=EXPEXP(LEVEXP)

                END IF
C
C=======================================================================
C
                IF (LEVEXD.NE.LEVEXP) THEN
C
                    LEVEXD=LEVEXP
                    LABEXP(LEVEXP)=SYMBEX(NUCLEU,LEVELS)
C
                END IF
C
C=======================================================================
C
             END DO
C
             NEXIST=0
C
         END IF
C
C=======================================================================
C
      END DO
C
C=======================================================================
C     Ready with the determination of the experimental spectrum
C=======================================================================
C
      FERMEX=0.5*(EXPUPP+EXPLOW)
      GAPEXP=     EXPUPP-EXPLOW
C
      IEXPUP=0
      IEXPDW=0
C
      UPPMAX=-1000.0
      UPPMIN=+1000.0
      DWNMIN=+1000.0
      DWNMAX=-1000.0
C
      DO IEXPER=1,LEVEXP
C
         IF (EXPEXP(IEXPER).GT.FERMEX) THEN
C
             IF (UPPMAX.LE.EXPEXP(IEXPER)) UPPMAX=EXPEXP(IEXPER)
             IF (UPPMIN.GT.EXPEXP(IEXPER)) UPPMIN=EXPEXP(IEXPER)
C
C            Definig level degeneracies out of spectrscopic labels
C
             CALL READEG(IDEGEX(IEXPER),LABEXP(IEXPER))
C
             IEXPUP=IEXPUP+IDEGEX(IEXPER)
C
         END IF
C
         IF (EXPEXP(IEXPER).LT.FERMEX) THEN
C
             IF (DWNMIN.GE.EXPEXP(IEXPER)) DWNMIN=EXPEXP(IEXPER)
             IF (DWNMAX.LT.EXPEXP(IEXPER)) DWNMAX=EXPEXP(IEXPER)
C
C            Definig level degeneracies out of spectrscopic labels
C
             CALL READEG(IDEGEX(IEXPER),LABEXP(IEXPER))
C
             IEXPDW=IEXPDW+IDEGEX(IEXPER)
C
         END IF
C
      END DO
C
      RANGUP=UPPMAX-UPPMIN
      RANGDW=DWNMAX-DWNMIN
C
                       DENSUP=0.0
      IF (IEXPUP.NE.0) DENSUP=REAL(IEXPUP)/RANGUP
C
                       DENSDW=0.0
      IF (IEXPDW.NE.0) DENSDW=REAL(IEXPDW)/RANGDW
C
C=======================================================================
C=======================================================================
C     Ready with the reconstruction of the experimental spectra
C=======================================================================
C=======================================================================
C
      IF (NEXIST.EQ.1) THEN
C
          WRITE(LSCREN,'(//,80(''|''),/,''|'',78X,''|'',/,
     *                   ''|   WARNING: The experimental data for '',
     *                   ''the nucleus Z='',I3,'' N='',I3,
     *                   '' <=> NONEXISTING  |'',/,''|'',78X,''|'',/,
     *                                                   80(''|''),/)')
     *
     *                                                   IZ_FIX,IN_FIX
C
          IF (LOGWRI.GT.0)
     *    WRITE(LOGFIL,'(//,80(''|''),/,''|'',78X,''|'',/,
     *                   ''|   WARNING: The experimental data for '',
     *                   ''the nucleus Z='',I3,'' N='',I3,
     *                   '' <=> NONEXISTING  |'',/,''|'',78X,''|'',/,
     *                                                   80(''|''),/)')
     *
     *                                                   IZ_FIX,IN_FIX
C
          STOP 'No experimental data to fit to for this nucleus'
C
      END IF
C
C=======================================================================
C
      NORMAL=0
C
      DO IEXPER=1,LEVEXP
         NORMAL=NORMAL+IDEGEX(IEXPER)
      END DO
C
C=======================================================================
C     Ordering the experimental levels before writing them on disc
C=======================================================================
C
      DO LEVELS=1,LEVEXP
         ENEAUX(LEVELS)=EXPEXP(LEVELS)
         INDAUX(LEVELS)=LEVELS
      END DO
C
      CALL ORDHEA(ENEAUX,INDAUX,LEVEXP,NDSPEC)
C
C=======================================================================
C
      IF (ISOSPI.EQ.0) THEN
C
          IF (LOGWRI.GT.0) WRITE(LOGCHN,'(/,''From PREPAR_EXPLEV'',/)')
C
          IF (LOGWRI.GT.0) 
     * 
     *    WRITE(LOGCHN,'(/,80(''#''),/,''#'',78X,''#'',/,
     *           ''#  Experimental Levels for the Spherical '',
     *           ''Nucleus    Z ='',I3,''  N ='',I3,
     *           '' [Neutrons]  #'',
     *           /,''#                                         '',
     *           ''          Version ===>>>>   '',A6,''   #'',/,
     *           ''#'',78X,''#'',/,80(''#''),/,''#'',78X,''#'',/,
     *           ''#     Level              Energy        '',
     *           ''    State               Degen.          #'',/,
     *           ''#'',78X,''#'')')
     *
     *                                           IZ_FIX,IN_FIX,WHATEX
C
          IF (LOGWRI.GT.5)
     *
     *        WRITE(LOGFIL,'(9X,
     *                       ''Experimental Levels for the Spherical '',
     *                       ''Nucleus    Z ='',I3,''  N ='',I3,
     *                       '' [Neutrons]'',/,T59,
     *                       ''Version ===>>>>   '',A6,/,9X,
     *                       ''Level   Energy    State    Degen.'')')
     *
     *                                           IZ_FIX,IN_FIX,WHATEX
      ELSE
C
          IF (LOGWRI.GT.0) WRITE(LOGCHP,'(/,''From PREPAR_EXPLEV'',/)')
C
          IF (LOGWRI.GT.0)
     * 
     *    WRITE(LOGCHP,'(/,80(''#''),/,''#'',78X,''#'',/,
     *          ''#  Experimental Levels for the Spherical '',
     *          ''Nucleus    Z ='',I3,''  N ='',I3,
     *          ''  [Protons]  #'',
     *          /,''#                                         '',
     *          ''          Version ===>>>>   '',A6,''   #'',/,
     *          ''#'',78X,''#'',/,80(''#''),/,''#'',78X,''#'',/,
     *          ''#     Level              Energy        '',
     *          ''    State               Degen.          #'',/,
     *          ''#'',78X,''#'')')
     *
     *                                          IZ_FIX,IN_FIX,WHATEX
C
          IF (LOGWRI.GT.5)
     *
     *        WRITE(LOGFIL,'(9X,
     *                       ''Experimental Levels for the Spherical '',
     *                       ''Nucleus    Z ='',I3,''  N ='',I3,
     *                       ''  [Protons]'',/,T59,
     *                       ''Version ===>>>>   '',A6,/,9X,
     *                       ''Level   Energy    State    Degen.'')')
     *
     *                                          IZ_FIX,IN_FIX,WHATEX
C
      END IF
C
C=======================================================================
C
      DO IEXPER=1,LEVEXP
         LABAUX(IEXPER)=LABEXP(INDAUX(IEXPER))
         IDEAUX(IEXPER)=IDEGEX(INDAUX(IEXPER))
      END DO
C
      DO IEXPER=1,LEVEXP
         LABEXP(IEXPER)=LABAUX(IEXPER)
         IDEGEX(IEXPER)=IDEAUX(IEXPER)
         EXPEXP(IEXPER)=ENEAUX(IEXPER)
      END DO
C
C=======================================================================
C     This is the normal run with full experimental information
C=======================================================================
C
      IF (ISOSPI.EQ.1) THEN
C
          LEVEXP_PROTON(INUCLI)=LEVEXP
C	  
          DO IEXPER=1,LEVEXP
             LABEXP_PROTON(INUCLI,IEXPER)=LABEXP(IEXPER)
             IDEGEX_PROTON(INUCLI,IEXPER)=IDEGEX(IEXPER)
             EXPEXP_PROTON(INUCLI,IEXPER)=EXPEXP(IEXPER)
          END DO
C
      END IF
C
      IF (ISOSPI.EQ.0) THEN
C
          LEVEXP_NEUTRS(INUCLI)=LEVEXP
C	  
          DO IEXPER=1,LEVEXP
             LABEXP_NEUTRS(INUCLI,IEXPER)=LABEXP(IEXPER)
             IDEGEX_NEUTRS(INUCLI,IEXPER)=IDEGEX(IEXPER)
             EXPEXP_NEUTRS(INUCLI,IEXPER)=EXPEXP(IEXPER)
          END DO
C
      END IF
C
C=======================================================================
C     and this is a user-requested massacre of the real experiment
C     It consists in removing certain experimental levels from the
C     fit in order to study their influence on the parametrisation
C=======================================================================
C @@@ IRENE!  Is it in agreement with WHATNOTAKE?
      IF (REMOVE_PROTON.EQ.'YES') THEN
C
          IF (ISOSPI.EQ.1) THEN
C
              IXREAL=0
C	  
              DO IEXPER=1,LEVEXP
                 DO I=1,NDRAUS
                    IF (LABEXP(IEXPER).EQ.LABPRO_REMOVE(I)) THEN
                        IF (INDPRO_LETOUT(I).EQ.0) THEN
                            IXREAL=IXREAL+1 
                            LABEXP_PROTON(INUCLI,IXREAL)=LABEXP(IEXPER)
                            IDEGEX_PROTON(INUCLI,IXREAL)=IDEGEX(IEXPER)
                            EXPEXP_PROTON(INUCLI,IXREAL)=EXPEXP(IEXPER)
                        END IF
                    END IF
                 END DO
              END DO
C
              LEVEXP_PROTON(INUCLI)=IXREAL
C
          END IF
C
      END IF
C_______________________________________________________________________
C
      IF (REMOVE_NEUTRS.EQ.'YES') THEN
C
          IF (ISOSPI.EQ.0) THEN
C
              IXREAL=0
C	  
              DO IEXPER=1,LEVEXP
                 DO I=1,NDRAUS
                    IF (LABEXP(IEXPER).EQ.LABNEU_REMOVE(I)) THEN
                        IF (INDNEU_LETOUT(I).EQ.0) THEN
                            IXREAL=IXREAL+1 
                            LABEXP_NEUTRS(INUCLI,IXREAL)=LABEXP(IEXPER)
                            IDEGEX_NEUTRS(INUCLI,IXREAL)=IDEGEX(IEXPER)
                            EXPEXP_NEUTRS(INUCLI,IXREAL)=EXPEXP(IEXPER)
                        END IF
                    END IF
                 END DO
              END DO
C
              LEVEXP_NEUTRS(INUCLI)=IXREAL
C
          END IF
C
      END IF
C
C=======================================================================
C     Final printout Table
C=======================================================================
C
      DO IEXPER=1,LEVEXP
C        
         IF (LOGWRI.GT.0) THEN
             IF (ISOSPI.EQ.1) THEN
C
                 WRITE(LOGCHP,'(''#'',I8,14X,F8.4,12X,A6,15X,I2,
     *                                         ''             #'')')
     *                                       IEXPER,EXPEXP(IEXPER),
     *                                              LABEXP(IEXPER),
     *                                              IDEGEX(IEXPER)
             END IF
C
             IF (ISOSPI.EQ.0) THEN
C
                 WRITE(LOGCHN,'(''#'',I8,14X,F8.4,12X,A6,15X,I2,
     *                                         ''             #'')')
     *                                       IEXPER,EXPEXP(IEXPER),
     *                                              LABEXP(IEXPER),
     *                                              IDEGEX(IEXPER)
             END IF 
         END IF
C
         IF (LOGWRI.GT.5) THEN
C
             WRITE(LOGFIL,'(I12,4X,F7.3,4X,A6,4X,I2)')
     *                                       IEXPER,EXPEXP(IEXPER),
     *                                              LABEXP(IEXPER),
     *                                              IDEGEX(IEXPER)
         END IF
C
      END DO ! End IEXPER
C
C=======================================================================
C
      IF (LOGWRI.GT.0) THEN
          IF (ISOSPI.EQ.1) THEN
              WRITE(LOGCHP,'(''#'',78X,''#'',/,80(''#''),/)')
          ELSE
              WRITE(LOGCHN,'(''#'',78X,''#'',/,80(''#''),/)')
          END IF
      END IF
C
C=======================================================================
C     The final level choice and possibly modified printout Table
C=======================================================================
C
      IF (REMOVE_PROTON.EQ.'YES') THEN
C
          IF (ISOSPI.EQ.1 .AND. LOGWRI.GT.0) THEN
C
              WRITE(LOGCHP,'(/,''Instead of the above: '',
     *                         ''The Actual Proton Table used'',/)')        
C
              WRITE(LOGCHP,'(/,80(''#''),/,''#'',78X,''#'',/,
     *              ''#  Experimental Levels for the Spherical '',
     *              ''Nucleus    Z ='',I3,''  N ='',I3,
     *              '' [Protons]   #'',
     *              /,''#                                         '',
     *              ''          Version ===>>>>   '',A6,''   #'',/,
     *              ''#'',78X,''#'',/,80(''#''),/,''#'',78X,''#'',/,
     *              ''#     Level              Energy        '',
     *              ''    State               Degen.          #'',/,
     *              ''#'',78X,''#'')')
     *                                       IZ_FIX,IN_FIX,WHATEX
C
              DO IEXPER=1,LEVEXP_PROTON(INUCLI)
C  
                 WRITE(LOGCHP,'(''#'',I8,14X,F8.4,12X,A6,15X,I2,
     *                                         ''             #'')')
     *                 IEXPER,EXPEXP_PROTON(INUCLI,IEXPER),
     *                        LABEXP_PROTON(INUCLI,IEXPER),
     *                        IDEGEX_PROTON(INUCLI,IEXPER)
C
              END DO
C
              WRITE(LOGCHP,'(''#'',78X,''#'',/,80(''#''),/)')
C
          END IF          
C
      END IF          
C_______________________________________________________________________
C
      IF (REMOVE_NEUTRS.EQ.'YES') THEN
C
          IF (ISOSPI.EQ.0 .AND. LOGWRI.GT.0) THEN
C
              WRITE(LOGCHN,'(/,''Instead of the above: '',
     *                         ''The Actual Neutron Table used'',/)')        
C
              WRITE(LOGCHN,'(/,80(''#''),/,''#'',78X,''#'',/,
     *              ''#  Experimental Levels for the Spherical '',
     *              ''Nucleus    Z ='',I3,''  N ='',I3,
     *              '' [Neutrons]  #'',
     *              /,''#                                         '',
     *              ''          Version ===>>>>   '',A6,''   #'',/,
     *              ''#'',78X,''#'',/,80(''#''),/,''#'',78X,''#'',/,
     *              ''#     Level              Energy        '',
     *              ''    State               Degen.          #'',/,
     *              ''#'',78X,''#'')')
     *
     *                                       IZ_FIX,IN_FIX,WHATEX
C        
              DO IEXPER=1,LEVEXP_NEUTRS(INUCLI)
C      
                 WRITE(LOGCHN,'(''#'',I8,14X,F8.4,12X,A6,15X,I2,
     *                                         ''             #'')')
     *                 IEXPER,EXPEXP_NEUTRS(INUCLI,IEXPER),
     *                        LABEXP_NEUTRS(INUCLI,IEXPER),
     *                        IDEGEX_NEUTRS(INUCLI,IEXPER)
C
              END DO
C         
              WRITE(LOGCHN,'(''#'',78X,''#'',/,80(''#''),/)')
C
          END IF         
C
      END IF         
C
C=======================================================================
C
      IF (LOGWRI.GT.4) THEN
          WRITE(LOGFIL,'(9X,''Exiting  PREPAR_EXPLEV'')')
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
      SUBROUTINE RDPROT_RMSRAD(NDRADZ,RADPRO,CHGRMS,RERROR,IZ_FIX,
     *                                                     IN_FIX,
     *                                              IFEXPE,IPRINT)
C
      CHARACTER
     *          FILNAM*256
C
      DIMENSION
     *          IZNUCL(1:NDRADZ),INNUCL(1:NDRADZ),I_MASS(1:NDRADZ)
      DIMENSION
     *          RMSCHG(1:NDRADZ),CHGERR(1:NDRADZ),RMSPRO(1:NDRADZ)
      DIMENSION
     *          DELRAD(1:NDRADZ),DELERR(1:NDRADZ),
     *          AUXRMS(1:NDRADZ),AUXERR(1:NDRADZ)
      DIMENSION
     *          SPINGS(1:NDRADZ)
C
      COMMON
     *       /RANRAN/ IRMFST,I_RAND,IRANDO
      COMMON
     *       /WHERBE/ IRUNMI,NEWOLD 
      COMMON
     *       /PRFILE/ ISCREN,LOGFIL,LOGCHP,LOGCHN,LOGHIP,LOGHIN,
     *                LOGPRO,LOGNEU,LOGWRI,LSCREN,IRESUL,ICONVE,
     *                NOUTPT,LOGAUX,LOGENE,LOGRAD,LOGGRD,LOGMSH,
     *                                                   LOGBIS
C
      DATA
     *     CHPROT / +0.743 /,
     *     CHNEUT / -0.119 /
      DATA
     *     N_READ / 80 /
C
C=======================================================================
C
C     Subroutine reading the experimental proton root-mean-square
C     radii from the data file ws_exper/exp_rms_radii one nucleus
C                                                       at a time
C=======================================================================
C
C     CHPROT and CHNEUT are the mean square radii of the charge
C     distributions in a proton and a neutron, respectively.
C     (cf. Eq.(3) J. Dobaczewski et al., Z.Phys A354, 27-35 (1996))
C
C=======================================================================
C
      RADPRO=-1.0
      CHGRMS=-1.0
      RERROR= 0.0
C
      IFOUND=0
C
C=======================================================================
C
      WRITE(NOUTPT,'(//3(80(1H*),/),                             /,
     *         9(1H*),''   READING  THE  EXPERIMENTAL  RESULTS '',
     *                 '' FROM  RDPROT_RMSRAD   ''
     *          9(1H*),//,                              3(80(1H*),/))')
C
C=======================================================================
C
      FILNAM='ws_exper/exp_rms_radii.d'
C
      OPEN(N_READ,FILE=FILNAM,STATUS='UNKNOWN',FORM='FORMATTED')
C
      IF (LOGWRI.GT.4) THEN
          WRITE(LOGFIL,'(12X,''Opening the file '',
     *                       ''ws_exper/exp_rms_radii for protons'')')
      END IF
C
      READ(N_READ,'(A80)') STRING
C
      NCOUNT=1
C
C=======================================================================
C
   1  CONTINUE
C
C=======================================================================
C
      READ (N_READ,*,END=2) IZNUCL(NCOUNT),INNUCL(NCOUNT),
     *                                     I_MASS(NCOUNT),
     *                                     SPINGS(NCOUNT),
     *                      DELRAD(NCOUNT),DELERR(NCOUNT),
     *                      AUXRMS(NCOUNT),AUXERR(NCOUNT)
C
C=======================================================================
C
      RMSCHG(NCOUNT)=AUXRMS(NCOUNT)+DELRAD(NCOUNT)
C
      IF (AUXRMS(NCOUNT).LT.0.0001) RMSCHG(NCOUNT)=-1.0
C
      CHGERR(NCOUNT)=AUXERR(NCOUNT)+DELERR(NCOUNT)
C
C=======================================================================
C
      IF (RMSCHG(NCOUNT).GT.0.0001) THEN
C
          RMSPRO(NCOUNT)=SQRT(RMSCHG(NCOUNT)**2
     *                  -CHPROT
     *                  -CHNEUT*REAL(INNUCL(NCOUNT))
     *                         /REAL(IZNUCL(NCOUNT)))
      ELSE
          RMSPRO(NCOUNT)=-1.0
      END IF
C
C=======================================================================
C
      IF (IZNUCL(NCOUNT).EQ.0) GO TO 2
C
C=======================================================================
C
      NCOUNT=NCOUNT+1
C
      GO TO 1
C
   2  CONTINUE
C
C=======================================================================
C     If the experimental information about our nucleus exists in
C     the data-bank, this is perfect
C=======================================================================
C
      DO I=1,NCOUNT
C
         IF ((IZNUCL(I).EQ.IZ_FIX).AND.(INNUCL(I).EQ.IN_FIX)) THEN
C
             IFOUND=1
C
             RADPRO=RMSPRO(I)
             CHGRMS=RMSCHG(I)
             RERROR=CHGERR(I)
             IFOUNX=IFOUND
C
         END IF
      END DO
C 
C=======================================================================
C     If the experimental information does not exist we calculate
C     the RMS radii using the Pomorski & Pomorski formula. 
C
C     ID has verified (23/03/2016) that for the PROTONS, 
C     the Pomorski & Pomorski formula fits better to the 
C     experimental results that the one from Dobaczewski 
C                                                article.
C=======================================================================
C
      IF (IFOUND.EQ.0) THEN
C
          A_MASS=REAL(IN_FIX+IZ_FIX)
          GEOMFA=SQRT(3.0000/5.0000)
C
          RADPOM=1.240*GEOMFA*(1.0+1.646/A_MASS
     *                            -0.191/A_MASS*REAL(IN_FIX-IZ_FIX))
     *                       *A_MASS**(1.0/3.0)
C
          RADPRO=SQRT(RADPOM**2-0.743
     *                         +0.119*REAL(IN_FIX)/REAL(IZ_FIX))     
          CHGRMS=RADPOM
C
C         Using the typical formula of RMSD with the ones known 
C         experimentally and the theory, we get RERROR=0.04
C
          RERROR=0.04   
C
          IFOUNX=IFOUND
C
          IFOUND=1
C      
      END IF
C
C=======================================================================
C
      IF (IFOUND.EQ.1) THEN
C
          IF (LOGWRI.GT.4) THEN
C
              WRITE(LOGFIL,'(12X,''Experimental proton radius '',
     *                           ''information: RADPRO='',F5.3,1X,
     *                           ''CHGRMS='',F5.3,1X,
     *                           ''RERROR='',F5.3,2X,''IFOUND='',I2)')
     *                             RADPRO,CHGRMS,RERROR,IFOUND
C
          END IF
          
          WRITE(NOUTPT,'(''IZ_FIX= '',I3,3X,''IN_FIX= '',I3)') IZ_FIX,
     *                                                         IN_FIX
          WRITE(NOUTPT,'(''RADPRO='',F5.3,'' CHGRMS='',F5.3,1X,
     *                   ''RERROR='',F5.3)') RADPRO,CHGRMS,RERROR
C
      ELSE
C
          IF (LOGWRI.GT.4) THEN
C
              WRITE(LOGFIL,'(9X,''No experimental information about '',
     *                          ''proton radii for the nucleus Z='',I3,
     *                                                     ''  N='',I3
     *                                                              )')
     *
     *                                                    IZ_FIX,IN_FIX
C
          END IF
C
      END IF
C
C=======================================================================
C
      CLOSE(N_READ)
C
      IF (LOGWRI.GT.4) THEN
          WRITE(LOGFIL,'(9X,''Exiting  RDPROT_RMSRAD'')')
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
      SUBROUTINE PROVID_RMS_NR(NDRADN,RADNEU,ERRNEU,IZ_FIX,IN_FIX,
     *                                       IDEFCN,IFEXPE,IPRINT)
C
      DIMENSION
     *          IZNUCL(1:NDRADN),INNUCL(1:NDRADN)
      DIMENSION
     *          ERRORN(1:NDRADN),RMSNEU(1:NDRADN)
C
      COMMON
     *       /RANRAN/ IRMFST,I_RAND,IRANDO
      COMMON
     *       /WHERBE/ IRUNMI,NEWOLD 
      COMMON
     *       /PRFILE/ ISCREN,LOGFIL,LOGCHP,LOGCHN,LOGHIP,LOGHIN,
     *                LOGPRO,LOGNEU,LOGWRI,LSCREN,IRESUL,ICONVE,
     *                NOUTPT,LOGAUX,LOGENE,LOGRAD,LOGGRD,LOGMSH,
     *                                                   LOGBIS 
C
C=======================================================================
C     Subroutine providing the experimental neutron r.m.s. radii
C=======================================================================
C
      IF (IDEFCN.GT.1) THEN
          RETURN
      END IF
C
C=======================================================================
C 
C     Experimental neutron rms radii: Table 1, J. Dobaczewski et al.,
C                                              Z.Phys. A354, 27 (1996) 
C
      RADNEU=-1.0
      ERRNEU= 0.0
C
      IFOUND=0
C
      NCOUNT=1
      IZNUCL(NCOUNT)=20
      INNUCL(NCOUNT)=20
      RMSNEU(NCOUNT)=3.49
      ERRORN(NCOUNT)=0.05
C
      NCOUNT=2
      IZNUCL(NCOUNT)=20
      INNUCL(NCOUNT)=28
      RMSNEU(NCOUNT)=3.63
      ERRORN(NCOUNT)=0.05
C
      NCOUNT=3
      IZNUCL(NCOUNT)=28
      INNUCL(NCOUNT)=30
      RMSNEU(NCOUNT)=3.70
      ERRORN(NCOUNT)=0.05
C
      NCOUNT=4
      IZNUCL(NCOUNT)=28
      INNUCL(NCOUNT)=36
      RMSNEU(NCOUNT)=3.91
      ERRORN(NCOUNT)=0.05
C
      NCOUNT=5
      IZNUCL(NCOUNT)=40
      INNUCL(NCOUNT)=50
      RMSNEU(NCOUNT)=4.29
      ERRORN(NCOUNT)=0.07
C
      NCOUNT=6
      IZNUCL(NCOUNT)=50
      INNUCL(NCOUNT)=66
      RMSNEU(NCOUNT)=4.69
      ERRORN(NCOUNT)=0.05
C
      NCOUNT=7
      IZNUCL(NCOUNT)=50
      INNUCL(NCOUNT)=74
      RMSNEU(NCOUNT)=4.85
      ERRORN(NCOUNT)=0.05
C
      NCOUNT=8
      IZNUCL(NCOUNT)=82
      INNUCL(NCOUNT)=126
      RMSNEU(NCOUNT)=5.59
      ERRORN(NCOUNT)=0.04
C
C=======================================================================
C     This in the case where the experimental information on the current
C                                                         nucleus exists
C=======================================================================
C
      DO I=1,NCOUNT
C
         IF (IZNUCL(I).EQ.IZ_FIX.AND.INNUCL(I).EQ.IN_FIX) THEN
C
             RADNEU=RMSNEU(I)
             ERRNEU=ERRORN(I)
C
             IFOUND=1
C
             IF (LOGWRI.GT.4) THEN
C
                 WRITE(LOGFIL,'(12X,
     *                         ''The corresponding neutron radius'',
     *                          3X,''RMSNEU='',F5.3,''  ERRNEU='',F5.3
     *                                                              )')
     *                               RADNEU,            ERRNEU
             END IF
C
         END IF
C
      END DO
C
C=======================================================================
C     In the opposite case i.e. the experimental information does not 
C     exist, we estimate the neutron radius using Pomorski & Pomorski
C     formula
C=======================================================================
C
      IF (IFOUND.EQ.0) THEN
C
C=======================================================================
C
C         It has been verified,  that the experimentally known  neutron
C         data points agree slightly better with results from  
C
C                              'PHENOM_RMSRAD',
C
C         rather than  with the formula(s) by Pomorski & Pomorski, from 
C         Nucl.Phys. A635, 1998, 484. We keep it however as a comment
C
CID:
C         The authors have verified it...: citing textually:
C
C         Formula (6) improves the expression used in [30] by adding
C         terms proportional to A^{-2}
C
C         Where formula (6) is the one from Dobaczewski et al., 
C         and [30] refers to Pomorska and Pomorski article.
CID.
C
C=======================================================================
C
C         Neutron RMS Radii coming from Nucl.Phys. A635, 1998, 484
C
C         A_MASS=REAL(IN_FIX+IZ_FIX)
C         GEOMFA=SQRT(3.0/5.0)
C
C         RADNEU=1.176*GEOMFA*(1.0+0.250*REAL(IN_FIX-IZ_FIX)/A_MASS
C     *                           +2.806/A_MASS)*A_MASS**(1.0/3.0)
C
C=======================================================================
C 
          IF (LOGWRI.GT.4) THEN
              WRITE(LOGFIL,'(12X,''Entering PHENOM_RMSRAD from '',
     *                           ''PROVID_RMS_NR - case IFOUND=0'')')
          END IF
C
          CALL PHENOM_RMSRAD(IZ_FIX,IN_FIX,RADNEU,RADUMY)
C
C         Using the typical formula of RMSD with the ones known 
C         experimentally and the theory, we get RERROR=0.07
C 
          ERRNEU=0.07
C
          IFOUND=1
C
      END IF
C
C=======================================================================
C
      IF (IFOUND.EQ.0) THEN
C
          IF (LOGWRI.GT.4) THEN
              WRITE(LOGFIL,'(9X,''No neutron  information  found '',
     *                          ''case IFOUND=0'')')
          END IF
C
      END IF
C
C=======================================================================
C
      IF (LOGWRI.GT.4) THEN
          WRITE(LOGFIL,'(9X,''Exiting  PROVID_RMS_NR'')')
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
      SUBROUTINE PHENOM_RMSRAD(IZ_FIX,IN_FIX,RADNEU,RADPRO)
C
      COMMON
     *       /PRFILE/ ISCREN,LOGFIL,LOGCHP,LOGCHN,LOGHIP,LOGHIN,
     *                LOGPRO,LOGNEU,LOGWRI,LSCREN,IRESUL,ICONVE,
     *                NOUTPT,LOGAUX,LOGENE,LOGRAD,LOGGRD,LOGMSH,
     *                                                   LOGBIS
C
      DATA
     *     R0PARN, R0PARP / 1.1760,+1.2140 /
      DATA
     *     ALFA1N, ALFA1P / 0.1341,-0.1233 /
      DATA
     *     ALFA2N, ALFA2P / 4.8280,-3.4840 /
      DATA
     *     XAPA1N, XAPA1P / 3.2640,+2.6390 /
      DATA
     *     XAPA2N, XAPA2P / -.7121,+0.2543 /
C
C=======================================================================
C
C     This subroutine defines simple fit expressions, based on rather 
C     extensive HF Skyrme calculations, with the SkP parameterisation
C     from J. Dobaczewski et al., Z.Phys. A354, 27 (1996) 
C
C @@@ IRENE, please verify on some occasion whether this is true???????????
C
C     It has been verified, that the experimentally known data points
C     agree slightly better with this,  rather than with the formulas 
C     by Pomorska & Pomorski from Nucl. Phys. A635, 1998, 484.
C
C=======================================================================
C
      A_MASS=REAL(IN_FIX+IZ_FIX)
      DIFFER=REAL(IN_FIX-IZ_FIX)
C
      EXPRES=1.0+XAPA1N/A_MASS
     *          +ALFA1N*DIFFER/A_MASS
     *          +XAPA2N/A_MASS**2
     *          +ALFA2N*DIFFER/A_MASS**2
C
      FACTOR=SQRT(3./5.)*R0PARN*A_MASS**(1./3.)
C      
      RADNEU=FACTOR*EXPRES     
C
C=======================================================================
C
      EXPRES=1.0+XAPA1P/A_MASS
     *          +ALFA1P*DIFFER/A_MASS
     *          +XAPA2P/A_MASS**2
     *          +ALFA2P*DIFFER/A_MASS**2
C
      FACTOR=SQRT(3./5.)*R0PARP*A_MASS**(1./3.)
C      
      RADPRO=FACTOR*EXPRES     
C
C=======================================================================
C
      IF (LOGWRI.GT.4) THEN
          WRITE(LOGFIL,'(16X,''From PHENOM_RMSRAD: Fit to HF with SkP'',
     *                 /,16X,''(since no real experimental data found)''
     *                                                               )')
          WRITE(LOGFIL,'(16X,''RADPRO='',F6.3,'' RADNEU='',F6.3)')
     *                         RADPRO,           RADNEU
      END IF
C
      IF (LOGWRI.GT.4) THEN
          WRITE(LOGFIL,'(12X,''Exiting  PHENOM_RMSRAD'')')
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
      SUBROUTINE PREPAR_DENEXP
C      
      INCLUDE  'MATDIM/N_SUMS.f'
      INCLUDE  'MATDIM/ND_RHO.f'
      INCLUDE  'MATDIM/N_NOYX.f'
      INCLUDE  'MATDIM/NDGAUS.f'
      INCLUDE  'MATDIM/NDNUCL.f'
      INCLUDE  'MATDIM/NDIM_N.f'
C
      EXTERNAL
     *         DENEXP_INTEGR,RMSEXP_INTEGR
C      
      CHARACTER 
     *          SYMBNU*5,NUCSYM*6  
C
      DIMENSION
     *          COEFFI(1:N_SUMS,1:N_NOYX),RADIUS(1:N_NOYX)
C
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
     *
     *       /EXPSYM/ SYMBNU(1:N_NOYX)
      COMMON
     *       /BESSEL/ COEFFI_BESSEL(1:N_SUMS,1:NDNUCL)
      COMMON
     *       /CUTOFF/ RADIUS_CUTOFF(1:NDNUCL) 
      COMMON
     *       /ARGUMN/ IARG_R,IARG_Z
      COMMON
     *       /SIMUSE/ AOSCIL,RADCUT,INUCLI
      COMMON
     *       /NUCLEU/ IX_FIX,IN_FIX,ISOSPI
      COMMON
     *       /MODIFR/ RMSEXP_PROTON(1:NDNUCL),
     *                RMSEXP_NEUTRS(1:NDNUCL),
     *                RMSTHE_PROTON,RMSTHE_NEUTRS,
     *                RMODIF_PROTON,RMODIF_NEUTRS,
     *                              PUSHIN,DELTAR
C
      COMMON
     *       /PRFILE/ ISCREN,LOGFIL,LOGCHP,LOGCHN,LOGHIP,LOGHIN,
     *                LOGPRO,LOGNEU,LOGWRI,LSCREN,IRESUL,ICONVE,
     *                NOUTPT,LOGAUX,LOGENE,LOGRAD,LOGGRD,LOGMSH,
     *                                                   LOGBIS
C
C=======================================================================
C
C     This SUBROUTINE reads the experimental data for charge-densities
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
      PINUMB=4.0*ATAN(1.0)
C
C=======================================================================
C 
      IF (LOGWRI.GT.4) THEN
          WRITE(LOGFIL,'(/,12X,''Entering READEX_DENSIT from '',
     *                         ''PREPAR_DENEXP, protons only'')')
      END IF
C
      CALL READEX_DENSIT(N_SUMS,N_NOYX,COEFFI,SYMBNU,LABELN,RADIUS)
C
C=======================================================================
C
      ISOSPI=1
C
C=======================================================================
C
C     Rearranging the vectors COEFFI and RADIUS to the new vectors
C     COEFFI_BESSEL and RADIUS_CUTOFF, respectively, in order to
C     have them referenced to INUCLI and not to I_NOYX.
C
      DO INUCLI=1,LDNUCL
         RADIUS_CUTOFF(INUCLI)=0.0
         DO I_SUMS=1,N_SUMS
            COEFFI_BESSEL(I_SUMS,INUCLI)=0.0
         END DO
      END DO
C
      DO INUCLI=1,LDNUCL
         DO I_NOYX=1,N_NOYX
C
            IF (LABELN(I_NOYX,1).EQ.NUMB_Z(INUCLI) .AND.
     *          LABELN(I_NOYX,2).EQ.NUMB_N(INUCLI)) THEN
C
                RADIUS_CUTOFF(INUCLI)=RADIUS(I_NOYX)
C
                DO I_SUMS=1,N_SUMS
                   COEFFI_BESSEL(I_SUMS,INUCLI)=COEFFI(I_SUMS,I_NOYX)
                END DO
C
            END IF
         END DO
      END DO
C
C=======================================================================
C
C     Checking if we obtain the particle number with the experimental
C     charge density. We also compare the rms radius that we obtain
C     using this density with RMSEXP_PROTON (read from RDPROT_RMSRAD)
C
C     The integration is done with the Simpson Method.
C
      IARG_R=1
      IARG_Z=0
C
      DO INUCLI=1,LDNUCL
C                                      INUCLI,AOSCIL,RADCUT pass via 
         AOSCIL=AOSCIL_PROTON(INUCLI)                       ! COMMON
         RADCUT=RADIUS_CUTOFF(INUCLI)
C
         DO I_NOYX=1,N_NOYX
C
            IF (LABELN(I_NOYX,1).EQ.NUMB_Z(INUCLI) .AND.
     *          LABELN(I_NOYX,2).EQ.NUMB_N(INUCLI)) THEN
C
                IZ_FIX=NUMB_Z(INUCLI)
                A=0.00001
                B=12.0
                N=40
C_______________________________________________________________________
C
                IF (LOGWRI.GT.4) THEN
                    WRITE(LOGFIL,'(/,12X,''Entering SIMPSN_INTEGR '',
     *                             ''from PREPAR_DENEXP to test '',
     *                             ''DENEXP'')')
                END IF
C
                CALL SIMPSN_INTEGR(A,B,N,DENEXP_INTEGR,SIMDEN,ERROR1)
C
                DIFFER=ABS(SIMDEN-IZ_FIX)
C
                IF (LOGWRI.GT.4) THEN
C
                    WRITE(LOGFIL,'(/,12X,''Test summary: DENSITY'')')
C
                    WRITE(LOGFIL,'(15X,''Integral of DENEXP for '',
     *                 ''INUCLI= '',I1,'' with N='',i3,1X,
     *                 ''from A=0 to B=12 is:'',2X,f20.13)') INUCLI,N,
     *                                                       SIMDEN
C
                    WRITE(LOGFIL,'(15X,''Nominal integration error '',
     *                          ''df/dx=1, in SIMPSN_INTEGR:      '',
     *                          ''ERRORS= '',F20.13)') ERROR1
C
                    IF (DIFFER.GT.10.0*ERROR1) THEN
C
                        WRITE(LOGFIL,'(15X,''WARNING: DIFFER='',F20.13,
     *                                 1X,''is > 10*ERROR1= '',F20.13)')
     *                                                  DIFFER,10*ERROR1
                    END IF
C
                END IF
C_______________________________________________________________________
C
                IF (LOGWRI.GT.4) THEN
                    WRITE(LOGFIL,'(/,12X,''Entering SIMPSN_INTEGR '',
     *                                   ''from PREPAR_DENEXP for '',
     *                                   ''RMSEXP'')')
                END IF
C
                CALL SIMPSN_INTEGR(A,B,N,RMSEXP_INTEGR,SIMRMS,ERROR2)
C
                ERRORT=(SIMRMS/SIMDEN)*(ERROR1/SIMDEN+ERROR2/SIMRMS)
C
                SIMRMS=SIMRMS/SIMDEN 
C
                DIFFER=(SQRT(SIMRMS)-RMSEXP_PROTON(INUCLI))
C
                IF (LOGWRI.GT.4) THEN
C
                    WRITE(LOGFIL,'(/,12X,''Test summary: RADIUS'')')
C
                    WRITE(LOGFIL,'(12X,''R.M.S. from RMSEXP for '',
     *               ''INUCLI= '',I3,'' from A=0 to B=12 is:'',2X,
     *               f20.13,/,T57,''R.M.S.-exp='',f22.13)')INUCLI,
     *                                               SQRT(SIMRMS),
     *                                      RMSEXP_PROTON(INUCLI)
C
                    WRITE(LOGFIL,'(12X,''The difference between '',
     *                               ''expected and obtained is:'',
     *                                               9X,f20.13)') DIFFER
C
                    WRITE(LOGFIL,'(9X,''Nominal integration error '',
     *                         ''df/dx=1 in SIMPSN_INTEGR: ERRORS= '',
     *                                               F20.13)') ERRORT
C
                END IF
C_______________________________________________________________________
C
            END IF
C
         END DO
C
      END DO
C 
C=======================================================================
C 
      IF (LOGWRI.GT.4) THEN
          WRITE(LOGFIL,'(/,''Exiting  PREPAR_DENEXP'')')
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
      SUBROUTINE READEX_DENSIT(N_SUMS,N_NOYX,COEFFI,SYMBNU,LABELN,
     *                                                     RADIUS)
C 
      PARAMETER 
     *         (N_UNIT=60)
C       
      CHARACTER 
     *          FILENA*80,KEYWOR*5,ENDNOY*5,SYMBNU*5
C
      DIMENSION
     *          RADIUS(1:N_NOYX),COEFFI(1:N_SUMS,1:N_NOYX)
      DIMENSION
     *          LABELN(1:N_NOYX,1:2)
      DIMENSION
     *          SYMBNU(1:N_NOYX)
C
      COMMON
     *         /PRFILE/ ISCREN,LOGFIL,LOGCHP,LOGCHN,LOGHIP,LOGHIN,
     *                  LOGPRO,LOGNEU,LOGWRI,LSCREN,IRESUL,ICONVE,
     *                  NOUTPT,LOGAUX,LOGENE,LOGRAD,LOGGRD,LOGMSH,
     *                                                     LOGBIS                                                 
C
C=======================================================================
C
C     This subroutine  reads the coefficients of the Bessel expansion 
C     expressing the nuclear charge density. Data in "data_densities"
C
C=======================================================================
C        
      FILENA='dens_exper/data_density.d'
C       
      OPEN(N_UNIT,FILE=FILENA,FORM='FORMATTED',IOSTAT=IRC) 
C     
      IF (IRC.NE.0) STOP 'Error in opening the file - READEX_DENSIT'
C
C=======================================================================
C       
      I_NOYX=0
      RADAUX=0.0
C
C=======================================================================
C      
   1  CONTINUE
                                                    KEYWOR='      '
      READ(N_UNIT,'(A5,2X,I3,2X,I3,4X,F7.4)',END=2) KEYWOR,
     *                                IZ_FIX,IN_FIX,RADAUX
C 
      IF (KEYWOR.EQ.'ENDGO') THEN
          IF (LOGWRI.GT.4) THEN
              WRITE(LOGFIL,'(12X,''Exiting  READEX_DENSIT'')')
          END IF
          RETURN
      END IF      
C      
      IF (KEYWOR(3:4).NE.'  ') THEN
          
          I_NOYX=I_NOYX+1
          LABELN(I_NOYX,1)=IZ_FIX
          LABELN(I_NOYX,2)=IN_FIX         
          SYMBNU(I_NOYX)=KEYWOR
          RADIUS(I_NOYX)=RADAUX
C          
          DO I_SERI=1,N_SUMS
             READ(N_UNIT,'(9X,E13.5)') COEFFI(I_SERI,I_NOYX)
          END DO
         
          READ(N_UNIT,*) ENDNOY
         
          IF (ENDNOY.EQ.'XXXXXX') THEN
              GO TO 1
          END IF
C     
      END IF
C      
      GO TO 1
C      
   2  CONTINUE
C
      STOP 'ENDGO STATEMENT MISSING'
C
C=======================================================================
C      
      RETURN
      END
C
C=======================================================================
C=======================================================================
C
      SUBROUTINE WRITIT_DENPLO(ND_RHO,N_NOYX,R_MESH,RHOEXP,SYMBNU,
     *                                                     LABELN)
C
      PARAMETER 
     *         (N_UNIT=70)
C           
      CHARACTER 
     *          FILENA*256,FILTAB*256,SYMBNU*5
C
      DIMENSION 
     *          RHOEXP(1:ND_RHO,1:N_NOYX),R_MESH(1:ND_RHO,1:N_NOYX)
      DIMENSION
     *          LABELN(1:N_NOYX,1:2)
      DIMENSION 
     *          SYMBNU(1:N_NOYX),
     *          FILTAB(1:N_NOYX)
C
      COMMON
     *         /PRFILE/ ISCREN,LOGFIL,LOGCHP,LOGCHN,LOGHIP,LOGHIN,
     *                  LOGPRO,LOGNEU,LOGWRI,LSCREN,IRESUL,ICONVE,
     *                  NOUTPT,LOGAUX,LOGENE,LOGRAD,LOGGRD,LOGMSH,
     *                                                     LOGBIS
C
C======================================================================
C
C     This subroutine writes the nuclear charge density in the files
C     "density_zXXX_nXXX.dat" in the format ready for plotting
C
C======================================================================
C       
      FILTAB(1)='dens_exper/density_z008_n008.dat'
      FILTAB(2)='dens_exper/density_z020_n020.dat'
      FILTAB(3)='dens_exper/density_z020_n028.dat'
      FILTAB(4)='dens_exper/density_z040_n050.dat'
      FILTAB(5)='dens_exper/density_z082_n126.dat'
C
C======================================================================
C     
      DO I_NOYX=1,N_NOYX 
C   
         FILENA=FILTAB(I_NOYX)
C
         IF (LOGWRI.GT.0) THEN
             WRITE(LOGFIL,'(15X,''Writing densities to: '',A)') FILENA
         END IF
C          
         OPEN(N_UNIT,FILE=FILENA,FORM='FORMATTED',IOSTAT=IRC) 
C     
         IF (IRC.NE.0) STOP 'Error in opening the file - WRITIT_DENPLO'
           
         WRITE(N_UNIT,'(A5,2X,I3,2X,I3)') SYMBNU(I_NOYX),
     *                   LABELN(I_NOYX,1),LABELN(I_NOYX,2)       
C          
         DO ID_RHO=1,ND_RHO
            WRITE(N_UNIT,'(10X,F7.4,2X,F7.4)') R_MESH(ID_RHO,I_NOYX),
     *                                         RHOEXP(ID_RHO,I_NOYX)
         END DO
C
         WRITE(N_UNIT,'(A6)') 'XXXXXX'
C          
         REWIND(N_UNIT)
C
         CLOSE (N_UNIT)
C        
      END DO
C
      IF (LOGWRI.GT.0) THEN
          WRITE(LOGFIL,'(12X,''Exiting  WRITIT_DENPLO'')')
      END IF
C
C======================================================================
C      
      RETURN
      END        
C
C=======================================================================
C=======================================================================
C
      FUNCTION DENEXP(XARGUM)
C      
      INCLUDE  'MATDIM/NDNUCL.f'
      INCLUDE  'MATDIM/N_SUMS.f'
C
      COMMON
     *       /BESSEL/ COEFFI_BESSEL(1:N_SUMS,1:NDNUCL)
      COMMON
     *       /ARGUMN/ IARG_R,IARG_Z
      COMMON
     *       /SIMUSE/ AOSCIL,RADCUT,INUCLI
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
      PINUMB=4*ATAN(1.0)
C
C=======================================================================
C
      IF (IARG_R.EQ.1) THEN
          RARGUM=XARGUM
      END IF
      IF (IARG_Z.EQ.1) THEN
          RARGUM=AOSCIL*SQRT(XARGUM)
      END IF
C
C=======================================================================
C
      IF (RARGUM.LE.RADCUT) THEN
          AUXILI=0.0
          DO I_SUMS=1,N_SUMS
             AUXARG=I_SUMS*PINUMB*RARGUM/RADCUT
             AUXILI=AUXILI+COEFFI_BESSEL(I_SUMS,INUCLI)
     *                    *SIN(AUXARG)/AUXARG
          END DO
      END IF
C
      IF (RARGUM.GT.RADCUT) THEN
          AUXILI=0.0
      END IF
C
C=======================================================================
C
      DENEXP=AUXILI
C
C=======================================================================
C
      RETURN
      END
C      
C=======================================================================
C=======================================================================
C
      SUBROUTINE READIN_PSEUDO
C
      INCLUDE 'MATDIM/NDNUCL.f'
      INCLUDE 'MATDIM/NDLEXP.f'
C
      CHARACTER
     *          KEYWOR*10,STRING*10,
     *          LABPRO_PSEUDO*6,LABNEU_PSEUDO*6
C
      COMMON
     *       /DENSTR/ IFDENS
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
     *         /PRFILE/ ISCREN,LOGFIL,LOGCHP,LOGCHN,LOGHIP,LOGHIN,
     *                  LOGPRO,LOGNEU,LOGWRI,LSCREN,IRESUL,ICONVE,
     *                  NOUTPT,LOGAUX,LOGENE,LOGRAD,LOGGRD,LOGMSH,
     *                                                     LOGBIS
C     
C=======================================================================
C
      N_READ_PSEUDO=8
C
      IF (IFDENS.EQ.0)
     *    OPEN(UNIT=N_READ_PSEUDO,FILE='ws15_pseudo-exper_IFDENS-0.d')
      IF (IFDENS.EQ.1)
     *    OPEN(UNIT=N_READ_PSEUDO,FILE='ws15_pseudo-exper_IFDENS-1.d')
C
      WRITE(LOGFIL,'(/,''Opening ws15_pseudo-exper.d for reading it'',
     *                                                            /)')
C
   1  CONTINUE     
C     
C=======================================================================
C     
      READ (N_READ_PSEUDO,'(A10)',END=2)  KEYWOR  
C
C=======================================================================
C
      IF (KEYWOR.EQ.'HOWMANYNUC') THEN
C
          READ(N_READ_PSEUDO,*) LDNUCL
          WRITE(LOGFIL,'(A10,I4)') KEYWOR,LDNUCL
C
          DO I=1,LDNUCL
C
             READ(N_READ_PSEUDO,*)STRING                  ! NUCLEUS-XX
             WRITE(LOGFIL,'(A,/)') STRING
C
             READ(N_READ_PSEUDO,*)IZ_PSE(I),IN_PSE(I)
             READ(N_READ_PSEUDO,*)STRING                  ! ISOSPI  NUMLEV
             READ(N_READ_PSEUDO,*)ISOSPI,LEVNEU_PSEUDO(I)
             READ(N_READ_PSEUDO,*)STRING                  ! PSEUDO...
C
             DO J=1,LEVNEU_PSEUDO(I)
                READ(N_READ_PSEUDO,'(11X,F17.13,2X,I6,2X,A6,2X,F6.4,
     *                                4X,I1)') ENENEU_PSEUDO(J,I),
     *                                         IDENEU_PSEUDO(J,I),
     *                                         LABNEU_PSEUDO(J,I),
     *                                         SIGNEU_PSEUDO(J,I),
     *                                         LEVTAK_NEUTRS(J,I)
C
                WRITE(LOGFIL,'(11X,F17.13,2X,I6,2X,A6,2X,F6.4,4X,I1)')
     *                               ENENEU_PSEUDO(J,I),
     *                               IDENEU_PSEUDO(J,I),
     *                               LABNEU_PSEUDO(J,I),
     *                               SIGNEU_PSEUDO(J,I),
     *                               LEVTAK_NEUTRS(J,I)
             END DO
C
             READ(N_READ_PSEUDO,*)STRING                  ! ISOSPI  NUMLEV
             READ(N_READ_PSEUDO,*)ISOSPI,LEVPRO_PSEUDO(I)
             READ(N_READ_PSEUDO,*)STRING                  ! PSEUDO...
C
             WRITE(LOGFIL,'()')
C
             DO J=1,LEVPRO_PSEUDO(I)
                READ(N_READ_PSEUDO,'(11X,F17.13,2X,I6,2X,A6,2X,F6.4,
     *                                4X,I1)')
     *                               ENEPRO_PSEUDO(J,I),
     *                               IDEPRO_PSEUDO(J,I),
     *                               LABPRO_PSEUDO(J,I),
     *                               SIGPRO_PSEUDO(J,I),
     *                               LEVTAK_PROTON(J,I)
C
                WRITE(LOGFIL,'(11X,F17.13,2X,I6,2X,A6,2X,F6.4,4X,I1)')
     *                               ENEPRO_PSEUDO(J,I),
     *                               IDEPRO_PSEUDO(J,I),
     *                               LABPRO_PSEUDO(J,I),
     *                               SIGPRO_PSEUDO(J,I),
     *                               LEVTAK_PROTON(J,I)
             END DO
C
             READ(N_READ_PSEUDO,*)STRING                  ! <<>>
C
          END DO
C
      END IF
C
C=======================================================================
C
      IF (KEYWOR.EQ.'EXECUTE-GO') THEN
C
C=======================================================================
C
          IF (LOGWRI.GT.0) THEN
              WRITE(LOGFIL,'(/,''Exiting  READIN_PSEUDO'')')
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
      WRITE(NOUTPT,'(''No "EXECUTE-GO" in the input data stream '',
     *               ''(in READIN_PSEUDO)'')')
C   
      STOP
     *
     *   'Wrong input structure, eof not allowed here, STOP from PSEUDO'
C
      END
C      
C=======================================================================
C=======================================================================
C
      SUBROUTINE CREATG_TITLES(NDTITL,NDNUCL,TITLES,TITLES_LATEXS,
     *                                PARPOT_UNITSS,NUCNAM_LATEXS)
C
      INCLUDE
     *         'MATDIM/NDPARS.f'
      CHARACTER
     *          TITLES*012,TITLES_LATEXS*050,NUCNAM_LATEXS*010,
     *                                       PARPOT_UNITSS*040
      DIMENSION
     *          TITLES(1:NDTITL),
     *          TITLES_LATEXS(1:NDTITL),
     *          NUCNAM_LATEXS(1:NDNUCL),
     *          PARPOT_UNITSS(1:NDPARS)
C
C=======================================================================
C
C     This subroutines defines the name of the titles once forever.
C     It is called from the main part of the program. In this way,
C     we avoid having useless time the large column of the defini-
C     tion of the titles each time we want them. 
C
C     It also defines the parameter and nuclei names under latex 
C     format (plotting-output purposes)
C
C=======================================================================
C
      TITLES(01)='  V0CENT_p  '  
      TITLES(02)='  R0CENT_p  '   
      TITLES(03)='  A0CENT_p  '  
      TITLES(04)='  V0SORB_p  '
      TITLES(05)='  R0SORB_p  '
      TITLES(06)='  A0SORB_p  '         
      TITLES(07)='  V0EFFM_p  '        
      TITLES(08)='  R0EFFM_p  '        
      TITLES(09)='  A0EFFM_p  '        
      TITLES(10)='  R0COUL_p  '         
C
      TITLES(11)='  XK_V0C_p  '	  
      TITLES(12)='  XK_R0C_p  '	  
      TITLES(13)='  XK_A0C_p  '	  
      TITLES(14)='  XK_LAM_p  '	  
      TITLES(15)='  XK_RSO_p  '	  
      TITLES(16)='  XK_ASO_p  '	  
      TITLES(17)='  XK_LEF_p  '	 
      TITLES(18)='  XK_REF_p  '	 
      TITLES(19)='  XK_AEF_p  '	 
      TITLES(20)='  XK_COU_p  '
C
      TITLES(21)='  V0CENT_n  '	  
      TITLES(22)='  R0CENT_n  '	  
      TITLES(23)='  A0CENT_n  '	  
      TITLES(24)='  V0SORB_n  '	  
      TITLES(25)='  R0SORB_n  '	  
      TITLES(26)='  A0SORB_n  '	  
      TITLES(27)='  V0EFFM_n  '	 
      TITLES(28)='  R0EFFM_n  '	 
      TITLES(29)='  A0EFFM_n  '
C
      TITLES(30)='  XK_V0C_n  '	  
      TITLES(31)='  XK_R0C_n  '	  
      TITLES(32)='  XK_A0C_n  '	  
      TITLES(33)='  XK_LAM_n  '	  
      TITLES(34)='  XK_RSO_n  '	  
      TITLES(35)='  XK_ASO_n  '	  
      TITLES(36)='  XK_LEF_n  '	 
      TITLES(37)='  XK_REF_n  '	 
      TITLES(38)='  XK_AEF_n  '
C
      TITLES(39)='  XLAMB_pp  '	  
      TITLES(40)='  XLAMB_pn  '	  
      TITLES(41)='  XLAMB_np  '	  
      TITLES(42)='  XLAMB_nn  '
C
      TITLES(43)='  YLAMB_PP  '	  
      TITLES(44)='  YLAMB_PN  '	  
      TITLES(45)='  YLAMB_NP  '	  
      TITLES(46)='  YLAMB_NN  '
C
      TITLES(47)='  CLAMB_PP  '	  
      TITLES(48)='  CLAMB_PN  '	  
      TITLES(49)='  CLAMB_NP  '	  
      TITLES(50)='  CLAMB_NN  '  
C
      TITLES(51)='  V0CENTRL  '
      TITLES(52)='  KAPA_V0C  '
      TITLES(53)='  R0CENTRL  ' 
      TITLES(54)='  KAPA_R0C  '
      TITLES(55)='  A0CENTRL  ' 
      TITLES(56)='  KAPA_A0C  '
C
      TITLES(57)='  LAMBDASO  '
      TITLES(58)='  KAPA_LSO  '
      TITLES(59)='  R0_SORBT  '
      TITLES(60)='  KAPA_RSO  '
      TITLES(61)='  A0_SORBT  '
      TITLES(62)='  KAPA_A0S  '
C_______________________________________________________________________
C
      TITLES(63)='  GAPEXP_P  '	  
      TITLES(64)='  GAPTHE_P  '	  
      TITLES(65)='  RMSEXP_P  '	  
      TITLES(66)='  RMSTHE_P  '	  
C
      TITLES(67)='  EABSAV_P  '	  
      TITLES(68)='  ERRMAX_P  '	  
      TITLES(69)='  DENUPP_P  '	  
      TITLES(70)='  DENSUP_P  '	  
      TITLES(71)='  DENLOW_P  '	  
      TITLES(72)='  DENSDW_P  '	  
C
      TITLES(73)='  CHISQU_P  '	  
      TITLES(74)='  INVERS_P  '
C
      TITLES(75)='  GAPEXP_N  '	  
      TITLES(76)='  GAPTHE_N  '	  
      TITLES(77)='  RMSEXP_N  '	  
      TITLES(78)='  RMSTHE_N  '	  
C
      TITLES(79)='  EABSAV_N  '	  
      TITLES(80)='  ERRMAX_N  '	  
      TITLES(81)='  DENUPP_N  '	  
      TITLES(82)='  DENSUP_N  '	  
      TITLES(83)='  DENLOW_N  '	  
      TITLES(84)='  DENSDW_N  '	  
C
      TITLES(85)='  CHISQU_N  '	  
      TITLES(86)='  INVERS_N  '	  
C
      TITLES(87)='  CHI__SUM  '
C
C=======================================================================
C
C     Parameter names with LaTeX
C
      TITLES_LATEXS(01)='$V_{p}^c$'  
      TITLES_LATEXS(02)='$r_{p}^c$'   
      TITLES_LATEXS(03)='$a_{p}^c$'  
      TITLES_LATEXS(04)='$\\lambda_{p}^{so}$'  
      TITLES_LATEXS(05)='$r_{p}^{so}$'
      TITLES_LATEXS(06)='$a_{p}^{so}$'         
      TITLES_LATEXS(07)='$V_{p}^{eff}$'        
      TITLES_LATEXS(08)='$r_{p}^{eff}$'           
      TITLES_LATEXS(09)='$a_{p}^{eff}$'           
      TITLES_LATEXS(10)='$r_{Coul}$'         
C
      TITLES_LATEXS(11)='$\\kappa_{V_{p}^c}$'  
      TITLES_LATEXS(12)='$\\kappa_{r_{p}^c}$'  
      TITLES_LATEXS(13)='$\\kappa_{a_{p}^c}$'  
      TITLES_LATEXS(14)='$\\kappa_{\\lambda_{p}^{so}}$'  
      TITLES_LATEXS(15)='$\\kappa_{r_{p}^{so}}$'
      TITLES_LATEXS(16)='$\\kappa_{a_{p}^{so}}$'
      TITLES_LATEXS(17)='$\\kappa_{V_{p}^{eff}}$'  
      TITLES_LATEXS(18)='$\\kappa_{r_{p}^{eff}}$'   
      TITLES_LATEXS(19)='$\\kappa_{a_{p}^{eff}}$'  
      TITLES_LATEXS(20)='$\\kappa_{r_{Coul}}$'  
C
      TITLES_LATEXS(21)='$V_{n}^c$'  
      TITLES_LATEXS(22)='$r_{n}^c$'   
      TITLES_LATEXS(23)='$a_{n}^c$'  
      TITLES_LATEXS(24)='$\\lambda_{n}^{so}$'  
      TITLES_LATEXS(25)='$r_{n}^{so}$'
      TITLES_LATEXS(26)='$a_{n}^{so}$'         
      TITLES_LATEXS(27)='$V_{n}^{eff}$'        
      TITLES_LATEXS(28)='$r_{n}^{eff}$'           
      TITLES_LATEXS(29)='$a_{n}^{eff}$'      
C
      TITLES_LATEXS(30)='$\\kappa_{V_{n}^c}$'  
      TITLES_LATEXS(31)='$\\kappa_{r_{n}^c}$'  
      TITLES_LATEXS(32)='$\\kappa_{a_{n}^c}$'  
      TITLES_LATEXS(33)='$\\kappa_{\\lambda_{n}^{so}}$'  
      TITLES_LATEXS(34)='$\\kappa_{r_{n}^{so}}$'
      TITLES_LATEXS(35)='$\\kappa_{a_{n}^{so}}$'
      TITLES_LATEXS(36)='$\\kappa_{V_{n}^{eff}}$'  
      TITLES_LATEXS(37)='$\\kappa_{r_{n}^{eff}}$'   
      TITLES_LATEXS(38)='$\\kappa_{a_{n}^{eff}}$' 
C
      TITLES_LATEXS(39)='$\\lambda_{\\pi\\pi}$'	  
      TITLES_LATEXS(40)='$\\lambda_{\\pi\\nu}$'	  
      TITLES_LATEXS(41)='$\\lambda_{\\nu\\pi}$'	  
      TITLES_LATEXS(42)='$\\lambda_{\\nu\\nu}$'
C
      TITLES_LATEXS(43)='$\\lambda_{\\pi\\pi}^{so}$'	  
      TITLES_LATEXS(44)='$\\lambda_{\\pi\\nu}^{so}$'	  
      TITLES_LATEXS(45)='$\\lambda_{\\nu\\pi}^{so}$'	  
      TITLES_LATEXS(46)='$\\lambda_{\\nu\\nu}^{so}$'
C
      TITLES_LATEXS(47)='$\\lambda_{\\pi\\pi}^{c}$'	  
      TITLES_LATEXS(48)='$\\lambda_{\\pi\\nu}^{c}$'	  
      TITLES_LATEXS(49)='$\\lambda_{\\nu\\pi}^{c}$'	  
      TITLES_LATEXS(50)='$\\lambda_{\\nu\\nu}^{c}$'
C
      TITLES_LATEXS(51)='$V_0^c$'
      TITLES_LATEXS(52)='$\\kappa_{V_0^c}$'
      TITLES_LATEXS(53)='$r_0^c$' 
      TITLES_LATEXS(54)='$\\kappa_{r_0^c}$'
      TITLES_LATEXS(55)='$a_0^c$' 
      TITLES_LATEXS(56)='$\\kappa_{a_0^c}$'
C
      TITLES_LATEXS(57)='$\\lambda_0^{so}$'
      TITLES_LATEXS(58)='$\\kappa_{\\lambda_0^{so}}$'
      TITLES_LATEXS(59)='$r_0^{so}$'
      TITLES_LATEXS(60)='$\\kappa_{r_0^{so}}$'
      TITLES_LATEXS(61)='$a_0^{so}$' 
      TITLES_LATEXS(62)='$\\kappa_{a_0^{so}}$'
C
C=======================================================================
C
C     Parameter units
C
      PARPOT_UNITSS(01)='MeV'
      PARPOT_UNITSS(02)='fm'
      PARPOT_UNITSS(03)='fm'
      PARPOT_UNITSS(04)='MeV\\,fm$^2/\\hbar^2$'
      PARPOT_UNITSS(05)='fm'
      PARPOT_UNITSS(06)='fm'
      PARPOT_UNITSS(07)='XX'
      PARPOT_UNITSS(08)='XX'
      PARPOT_UNITSS(09)='XX'
      PARPOT_UNITSS(10)='fm'
C
      DO I=11,20
         PARPOT_UNITSS(I)='XX'
      END DO
C
      PARPOT_UNITSS(21)='MeV'
      PARPOT_UNITSS(22)='fm'
      PARPOT_UNITSS(23)='fm'
      PARPOT_UNITSS(24)='MeV\\,fm$^2/\\hbar^2$'
      PARPOT_UNITSS(25)='fm'
      PARPOT_UNITSS(26)='fm'
      PARPOT_UNITSS(27)='XX'
      PARPOT_UNITSS(28)='XX'
      PARPOT_UNITSS(29)='XX'
C
      DO I=30,38
         PARPOT_UNITSS(I)='XX'
      END DO
C
      PARPOT_UNITSS(39)='MeV\\,fm$^5/\\hbar^2$'
      PARPOT_UNITSS(40)='MeV\\,fm$^5/\\hbar^2$'
      PARPOT_UNITSS(41)='MeV\\,fm$^5/\\hbar^2$'
      PARPOT_UNITSS(42)='MeV\\,fm$^5/\\hbar^2$'
C
      DO I=43,62
         PARPOT_UNITSS(I)='XX'
      END DO
C
C=======================================================================
C
C     Nuclei names
C
      NUCNAM_LATEXS(01)='$^{16}$O'
      NUCNAM_LATEXS(02)='$^{40}$Ca'
      NUCNAM_LATEXS(03)='$^{48}$Ca'
      NUCNAM_LATEXS(04)='$^{56}$Ni'
      NUCNAM_LATEXS(05)='$^{90}$Zr'
      NUCNAM_LATEXS(06)='$^{132}$Sn'
      NUCNAM_LATEXS(07)='$^{146}$Gd'
      NUCNAM_LATEXS(08)='$^{208}$Pb'
C
C=======================================================================
C
      RETURN
      END
C      
C=======================================================================
C=======================================================================
C
      SUBROUTINE LATEXS_LABELS(LABEXP,LABTEX)
C
      CHARACTER
     *          LABEXP*6,LABTEX*11
C
C=======================================================================
C
C     This subroutine transforms the energy labels written in the form
C     like 1h11/2 to latex form, i.e.: $1h_{11/2}$
C
C     LABEXP - input label 1h11/2
C
C     LABTEX - output label $1h_{11/2}$
C
C=======================================================================
C   
      IF (LABEXP(4:4).EQ.'/') THEN
C
          WRITE(LABTEX,'(''$'',A2,''_{'',A3,''}$'')')
     *                                    LABEXP(1:2),
     *                                    LABEXP(3:5)
      ELSE
C
          WRITE(LABTEX,'(''$'',A2,''_{'',A4,''}$'')')
     *                                    LABEXP(1:2),
     *                                    LABEXP(3:6)
      END IF
        
      DO K=1,11
         IF (LABTEX(K:K).EQ.'/') LABTEX(K:K)='@'
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
      SUBROUTINE LATEXS_TITLES(AUXTIT,TITTEX)
C          
      CHARACTER
     *          AUXTIT*10,TITTEX*14
C
C=======================================================================
C
C     This subroutine transforms the parameter titles written according
C     to the style e.g. V0CENT_P to the latex form i.e.: $V_{0,cent}^p$
C
C     AUXTIT - input title V0CENT_P
C
C     TITTEX - output title $V_{0,cent}^p$
C
C=======================================================================
C   
      IF (AUXTIT(4:7).EQ.'CENT') THEN
C
          IF (AUXTIT(8:9).EQ.'_p') THEN
              WRITE(TITTEX,'(''$'',A1,''_{0,cent}^p$'')')AUXTIT(2:2)
          END IF
          IF (AUXTIT(8:9).EQ.'_n') THEN
              WRITE(TITTEX,'(''$'',A1,''_{0,cent}^n$'')')AUXTIT(2:2)
          END IF
C
      END IF
C_______________________________________________________________________
C        
      IF (AUXTIT(4:7).EQ.'SORB') THEN
C
          IF (AUXTIT(8:9).EQ.'_p') THEN
             WRITE(TITTEX,'(''$'',A1,''_{0,so}^p$'')')AUXTIT(2:2)
          END IF
          IF (AUXTIT(8:9).EQ.'_n') THEN
            WRITE(TITTEX,'(''$'',A1,''_{0,so}^n$'')')AUXTIT(2:2)
          END IF
C
      END IF
C_______________________________________________________________________
C        
      IF (AUXTIT(4:7).EQ.'COUL') THEN
          WRITE(TITTEX,'(''$'',A1,''_{0,coul}$'')')AUXTIT(2:2)
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

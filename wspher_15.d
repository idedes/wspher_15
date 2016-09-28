
                  INPUT DATA SET FOR THE WSPHER_09.f PROGRAM

====================================================================================================                                      

IFHOWTORUN  ISIMPL IFTEST IFFITS IFMESH ISPE_P IFMOCA IF_DIR IFDENS IFTENS ISCREN LOGWRI
               0      0      1      0     0      0      0      1      0       1      1

IFTAKECHI2  IF_SPE IF_RMS IF_GAP IF_FER IF_DEN IF_RHO IF_INV
              1      1      0      0      0      0      0

EXPER_FILE  IFDEEP IFPRON
               0      1

HOWMANYNUC  LDNUCL
              8

ITAKENUCLS  TAKE01 TAKE02 TAKE03 TAKE04 TAKE05 TAKE06 TAKE07 TAKE08 TAKE09 TAKE10 TAKE11 TAKE12 TAKE13 TAKE14
               0      0      0      0      0      0      0      1      0      0      0      0      0      0

NUCLEIDATA  NUCL01 NUCL02 NUCL03 NUCL04 NUCL05 NUCL06 NUCL07 NUCL08 NUCL09 NUCL10 NUCL11 NUCL12 NUCL13 NUCL14
               16O   40Ca   48Ca   56Ni   90Zr  132Sn  146Gd  208Pb  278Fl  300Fl  342Fl 290Ubh 312Ubh 354Ubh
               8     20     20     28     40      50     64     82     114    114    114    126    126    126
               8     20     28     28     50      82     82    126     164    184    228    164    184    228

IF-PAIRING  IFGPAI IFDPAI
               0      0

BISECTSTEP  BISTEP
            +5.00

DELTAFACTO  DEL_Ox DEL_C0 DEL_C8 DEL_Ni DEL_Zr DEL_Sn DEL_Gd DEL_Pb DEL_09 DEL_10 DEL_11 DEL_12 DEL_13 DEL_14 (1-protons, 2-neutrons)
              1.00   1.00   1.00   1.00   1.00   1.00   1.00  15.00   1.00   1.00   1.00   1.00   1.00   1.00  protons
              1.00   1.00   1.00   1.00   1.00   1.00   1.00  15.00   1.00   1.00   1.00   1.00   1.00   1.00  neutrons

TENSOPTION  ICENTT ISORBT ITENSR
               0      0      0  
          
PARAMSBASE  NSHELL_PROTON NSHELL_NEUTRS NGAUSS HOMEG0 
                  18            18        40   41.000 
              
EXPOPTIONS  WHATEX (EXTRUE, EXOROS, SHELLM, EXPNEW)             
            EXPNEW
              
WHATPRINTS  ENEMAX             
            +5.000

SVDWHATCUT  SVDCUT             
            5.0E-1

GNUPLPRINT  PRIMAX (maximum chi^2)
             150.0

IFLAMBVSOD  IFPAR1  IFPAR2  IFPAR3  IFPAR4
              +1      -1      -1      -1

IFLAMBVSOT  IFPAR5  IFPAR6  IFPAR7  IFPAR8
              +3      -3      -3      -3

IFLAMBVCNT  IFPA09  IFPA10  IFPA11  IFPA12
              +3      -3      -3      -3
               
UNITLAMBDA  UNITLA
             100.0

IFKAPPAPAR  IFK_VC IFK_VS IFK_RC IFK_RS IFK_AC IFK_AS
               1      0      0      0      0      0

IFPARAMCOR  IFCORR
               0

PARAMCORRE  IFRCVC IFRCAC IFVCAC IFRSVS IFRSAS IFVSAS
               0      0      0      0      0      0
             
            
            Single run or starting values (default below: Old Universal)
            ============================================================

CENTRALPOT  V0CENT  R0CENT  A0CENT  XK_V0C  XK_R0C  XK_A0C  R0COUL  XK_COU  COUFAC     
           -58.238711841    1.269887550    0.750000000   0.00    0.00    0.00   1.2400   0.00     0.90
           -43.364600579    1.287135308    0.750000000   0.00    0.00    0.00

SPIN-ORBIT  V0SORB  R0SORB  A0SORB  XK_LAM  XK_RSO  XK_ASO 
            28.695689746    1.225143723    0.700000000   0.00    0.00    0.00
            32.527698991    0.875100158    0.700000000   0.00    0.00    0.00

KAPPA-CENT  V0CENT  XK_V0C  R0CENT  XK_R0C  A0CENT  XK_A0C 
            -50.801656       0.692044  1.2736  0.700   0.700     0.1

KAPPA-SORB  V0SORB  XK_LAM  R0SORB  XK_RSO  A0SORB  XK_ASO
            30.611694      -0.295883  1.2736  0.700   0.700     0.1

EFFECTMASS  V0EFFM  R0EFFM  A0EFFM  XK_LEF  XK_REF  XK_AEF 
            0.0000  0.0000  0.0000   0.00    0.00    0.00       
            0.0000  0.0000  0.0000   0.00    0.00    0.00   

DENSTSORBI  XLAMPP  XLAMPN  XLAMNP  XLAMNN
            159.2206       159.2206       159.2206       159.2206

TENSOSORBI  YLAMPP  YLAMPN  YLAMNP  YLAMNN
            -62.41  +54.91  +54.91  -62.41

TENSOCENTR  CLAMPP  CLAMPN  CLAMNP  CLAMNN
            42.7083      -42.7083      -42.7083       42.7083

SHIFTCENTR  SHIFOx  SHIFC0  SHIFC8  SHIFNi  SHIFZr  SHIFSn  SHIFGd  SHIFPb  SHIFPb  SHIFPb  SHIFPb  SHIFPb  SHIFPb  SHIFPb
            +0.000  +0.000  +0.000  +0.000  +0.000  +0.000  +0.000  +0.000  +0.000  +0.000  +0.000  +0.000  +0.000  +0.000
            +0.000  +0.000  +0.000  +0.000  +0.000  +0.000  +0.000  +0.000  +0.000  +0.000  +0.000  +0.000  +0.000  +0.000

    
            Control parameters for the minimisation subroutines
            ===================================================
            
            Role of STPMAX:   DIRNWT(I)=DIRNWT(I)*STPMAX/DIRNOR
                              STPMAX=STEPMX*MAX(Length(X),LDIMEN)
 
WHICHMINIM  IRMFST  I_RAND  LDRAND  I_SEED  NEWSED  IMTALK  IRUNMI
              1       1       20     7893     977      1       1

MINIMEVALS  ITECHI
              900 

TOLERANCES  DWFACT  UPFACT  PUSHIN
            0.9000  1.1000  0.0000

STOPTOLERS  TOLERF  TOLERX  TOLERG  FACTOR
            1.E-08  1.E-08  1.E-03  100.00

STOPINCOND  LDLAST  EPSLAS
              10    1.E-04
            
            Choice and definition of the chi-square subroutines
            ===================================================
            
MINWEIGHTS  WEICOR  WEIRAD  WEIINV  WEIFER  WEIGAP  WEIWEI  WEIDIF  WEIABS  WEIMAX  WEIDUP  WEIDDW  WEIRHO
             0.000    0.00   0.000   0.000   0.000   0.000   0.000  00.000   0.000   0.000   0.000   01.000   

CHOICEOFHI  CHIDEF Possible values: HICORR, HIWEIG, DIFWEI, EABSAV, ERRMAX
            HIWEIG
            
CHIWEI_PRO  WEI_Ox  WEI_Ca  WEI_Ca  WEI_Ni  WEI_Zr  WEI_Sn  WEI_Gd  WEI_Pb
            13.00    5.20   4.3333  3.7143  2.3111  1.5758  1.4247    1.00

CHIWEI_NEU  WEI_Ox  WEI_Ca  WEI_Ca  WEI_Ni  WEI_Zr  WEI_Sn  WEI_Gd  WEI_Pb
            13.00    5.20   4.3333  3.7143  2.3111  1.5758  1.4247    1.00

RADWEI_PRO  WEI_Ox  WEI_Ca  WEI_Ca  WEI_Ni  WEI_Zr  WEI_Sn  WEI_Gd  WEI_Pb
            50.000  50.000  50.000  0.0000  0.0000  50.000  0.0000  500.00

RADWEI_NEU  WEI_Ox  WEI_Ca  WEI_Ca  WEI_Ni  WEI_Zr  WEI_Sn  WEI_Gd  WEI_Pb
            50.000  50.000  50.000  0.0000  0.0000  50.000  0.0000  500.00
            
            RMS values when fitting only with 16O, 40-48Ca, 132Sn and 208Pb with CHIWEI=208/A
            (without radius)
            ===================================================

RMSVAL_PRO  Pro_Ox  Pro_Ca  Pro_Ca  Pro_Ni  Pro_Zr  Pro_Sn  Pro_Gd  Pro_Pb
            0.5896  0.8072  0.8609  1.0000  1.0000  0.4359  1.0000  0.4637

RMSVAL_NEU  Neu_Ox  Neu_Ca  Neu_Ca  Neu_Ni  Neu_Zr  Neu_Sn  Neu_Gd  Neu_Pb
            0.3583  0.4460  0.8150  1.0000  1.0000  0.8120  1.0000  0.6390

RMSGLO_REF  RMSPRO  RMSNEU
            0.6090  0.6762
            
            RMS values when fitting Individually each nucleus (without radius)
            ===================================================

RMSIND_PRO  Pro_Ox  Pro_Ca  Pro_Ca  Pro_Ni  Pro_Zr  Pro_Sn  Pro_Gd  Pro_Pb
            0.0998  0.0368  0.0866  0.1102  0.2428  0.1005  0.1995  0.0831

RMSIND_NEU  Neu_Ox  Neu_Ca  Neu_Ca  Neu_Ni  Neu_Zr  Neu_Sn  Neu_Gd  Neu_Pb
            0.1043  0.1821  0.2403  0.1303  0.5463  0.2848  0.0280  0.1423
            
            Minimisation Section (Central / Spin-Orbit / Effective Mass)
            ============================================================

CENTLIMITS  V0CMIN  V0CMAX  A0CMIN  A0CMAX  R0CMIN  R0CMAX  CouMIN  CouMAX        
            -80.00  -30.00  0.6000  0.8000  1.0000  1.3000  1.0000  1.3000

SORBLIMITS  XL_MIN  XL_MAX  A0SMIN  A0SMAX  R0SMIN  R0SMAX        
            20.000  40.000  0.6000  0.8000  0.8000  1.3000

EFFMLIMITS  XEFMIN  XEFMAX  AEFMIN  AEFMAX  REFMIN  REFMAX        
            20.000  50.000  0.2000  1.4000  0.3500  1.4500

DENSLIMITS  XPPMIN  XPPMAX  XPNMIN  XPNMAX  XNPMIN  XNPMAX  XNNMIN  XNNMAX  
            -200.0   200.0  -200.0   200.0  -200.0   200.0  -200.0   200.0 

TENSLIMITS  YPPMIN  YPPMAX  YPNMIN  YPNMAX  YNPMIN  YNPMAX  YNNMIN  YNNMAX      
            -200.0   200.0  -200.0   200.0  -200.0   200.0  -200.0   200.0 

TNCNLIMITS  CPPMIN  CPPMAX  CPNMIN  CPNMAX  CNPMIN  CNPMAX  CNNMIN  CNNMAX      
            -200.0   200.0  -200.0   200.0  -200.0   200.0  -200.0   200.0 
            
                    Steps for the derivatives [potentials]
                    ======================================

DERPARSTPS  DL_V0C  DL_A0C  DL_R0C  DL_XLA  DL_A0S  DL_R0S  DL_XRM  DL_ARM  DL_RRM  DL_COU 
             0.001  0.0001  0.0001   0.001  0.0001  0.0001    0.25    0.02    0.05    0.01  

             0.02    0.02    0.02    0.02    0.02    0.02    0.02    0.02     0.02    0.02    0.02    0.02      
            DL_XPP  DL_XPN  DL_XNP  DL_XNN  DL_YPP  DL_YPN  DL_YNP  DL_YNN  DL_CPP  DL_CPN  DL_CNP  DL_CNN

            Minimisation Section (Central / Spin-Orbit / Effective Mass)
                         [Isospin Dependence - kappas]
                         =============================

CENTKAPPAS  XKCMIN  XKCMAX  XACMIN  XACMAX  XRCMIN  XRCMAX  XCoMIN  XCoMAX
            -1.000  +1.000  -1.000  +1.000  -1.000  +1.000  -1.000  +1.000        

SORBKAPPAS  XKSMIN  XKSMAX  XASMIN  XASMAX  XRSMIN  XRSMAX
            -1.000  +1.000  -1.000  +1.000  -1.000  +1.000         

EFFMKAPPAS  XKEMIN  XKEMAX  XAEMIN  XAEMAX  XREMIN  XREMAX        
            -1.000  +1.000  -1.000  +1.000  -1.000  +1.000         
            
                  Steps for the derivatives [isospin/kappas] for protons
                  ======================================================

DERKAPSTPS  DLKV0C  DLKA0C  DLKR0C  DLKXLA  DLKA0S  DLKR0S  DLKXRM  DLKARM  DLKRRM  DLKCOU
            0.0100  0.0100  0.0100  0.0100  0.0100  0.0100  0.0100  0.0100  0.0100  0.0100

            Minimisation Section (Central / Spin-Orbit)
                         [Isospin Dependence - kappas (new version)]
                         =============================

KAPA-C-LIM  V0CMIN  V0CMAX  XKCMIN  XKCMAX  A0CMIN  A0CMAX  XACMIN  XACMAX  R0CMIN  R0CMAX  XRCMIN  XRCMAX
            -70.00  -30.00  -1.000  +1.000  0.4000  1.1000  -1.000  +1.000  0.9000  1.5000  -1.000  +1.000

KAPASO-LIM  XL_MIN  XL_MAX  XKSMIN  XKSMAX  A0SMIN  A0SMAX  XASMIN  XASMAX  R0SMIN  R0SMAX  XRSMIN  XRSMAX        
            20.000  40.000  -1.000  +1.000  0.4000  1.1000  -1.000  +1.000  0.5000  1.5000  -1.000  +1.000

STEPSKAPPA  DL_V0C  DLKV0C  DL_A0C  DLKA0C  DL_R0C  DLKR0C  DL_VSO  DLKXLA  DL_A0S  DLKA0S  DL_R0S  DLKR0S
             0.001  0.0100  0.0001  0.0100  0.0001  0.0100   0.001  0.0100  0.0001  0.0100  0.0001  0.0100
                    

                  Proton parameters
		  =================
		  
IFTAKEPARP  IF_V0C  IF_A0C  IF_R0C  IF_XLA  IF_A0S  IF_R0S  IF_XRM  IF_ARM  IF_RRM  IF_COU
               0       1       1       0       0       0       0       0       0       0

IFTAKEKAPP  IFKV0C  IFKA0C  IFKR0C  IFKXLA  IFKA0S  IFKR0S  IFKXRM  IFKARM  IFKRRM  IFKCOU
               0       0       0       0       0       0       0       0       0       0  
                    
                  Neutron parameters
                  ==================
		  
IFTAKEPARN  IF_V0C  IF_A0C  IF_R0C  IF_XLA  IF_A0S  IF_R0S  IF_XRM  IF_ARM  IF_RRM
               0       1       1       0       0       0       0       0       0

IFTAKEKAPN  IFKV0C  IFKA0C  IFKR0C  IFKXLA  IFKA0S  IFKR0S  IFKXRM  IFKARM  IFKRRM
               0       0       0       0       0       0       0       0       0  
            
                  Spin-orbit and tensor lambda (density dependent)
                  ================================================
		  
IFTAKEDENS  IF_XPP  IF_XPN  IF_XNP  IF_XNN
              1       0       0       0

IFTAKETENS  IF_YPP  IF_YPN  IF_YNP  IF_YNN  
              0       0       0       0 
              
IFTAKECNTN  IF_CPP  IF_CPN  IF_CNP  IF_CNN
              0       0       0       0  
                    
                  Kappa parameters
                  ================

IFTAKEKAPC  IFTV0C  IFTKVC  IFTR0C  IFTKRC  IFTA0C  IFTKAC
              1       1       0       0       0       0

IFTAKEKAPS  IFTV0S  IFTKVS  IFTR0S  IFTKRS  IFTA0S  IFTKAS
              0       0       0       0       0       0

            ==================================================
            ==================================================
            
    
            Control parameters for the MESH option subroutines (NOPXXX - Number Of Points for XXX)
            ==================================================
            
CENTMESHLM  V0CMIN  V0CMAX  NOPV0C  A0CMIN  A0CMAX  NOPA0C  R0CMIN  R0CMAX  NOPR0C  CouMIN  CouMAX  NOPCou
            -90.00  -10.00    81    0.1000  1.5000    29    0.8000  2.0000   241    1.0000  1.3000    20

SORBMESHLM  XL_MIN  XL_MAX  NOP_XL  A0SMIN  A0SMAX  NOPA0S  R0SMIN  R0SMAX  NOPR0S      
            05.000  70.000   261    0.1000  1.7000   321    0.4000  2.0000   321

EFFMMESHLM  XEFMIN  XEFMAX  NOPXEF  AEFMIN  AEFMAX  NOPAEF  REFMIN  REFMAX  NOPREF      
            20.000  50.000    20    0.2000  1.4000    20    0.3500  1.4500    20  

DENSMESHLM  XPPMIN  XPPMAX  NOPXPP  XPNMIN  XPNMAX  NOPXPN  XNPMIN  XNPMAX  NOPXNP  XNNMIN  XNNMAX  NOPXNN  
            -5.000  +15.00    11    -5.000  +15.00    11    -5.000  +15.00    11    -5.000  +15.00    11    

TENSMESHLM  YPPMIN  YPPMAX  NOPYPP  YPNMIN  YPNMAX  NOPYPN  YNPMIN  YNPMAX  NOPYNP  YNNMIN  YNNMAX  NOPYNN
            -5.000  +15.00    21    -5.000  +15.00    21    -5.000  +15.00    21    -5.000  +15.00    21     

CNTNMESHLM  CPPMIN  CPPMAX  NOPCPP  CPNMIN  CPNMAX  NOPCPN  CNPMIN  CNPMAX  NOPCNP  CNNMIN  CNNMAX  NOPCNN
            -3.000  23.500    20    -3.000  23.500    20    -3.000  23.500    20    -3.000  23.500    20    

CENTKAPMSH  XKCMIN  XKCMAX  NOPXKC  XACMIN  XACMAX  NOPXKC  XRCMIN  XRCMAX  NOPXRC  XCoMIN  XCoMAX  NOPXCo
            -1.000  +1.000    20    -1.000  +1.000    20    -1.000  +1.000    20    -1.000  +1.000    20

SORBKAPMSH  XKSMIN  XKSMAX  NOPXKS  XASMIN  XASMAX  NOPXAS  XRSMIN  XRSMAX  NOPXRS
            -1.000  +1.000    20    -1.000  +1.000    20    -1.000  +1.000    20           

EFFMKAPMSH  XKEMIN  XKEMAX  NOPXKE  XAEMIN  XAEMAX  NOPXAE  XREMIN  XREMAX  NOPXAE      
            -1.000  +1.000    20    -1.000  +1.000    20    -1.000  +1.000    20  

                  Kappa parametrization (central + spin-orbit new version)
		          =================

CNTVRAMESH  V0CMIN  V0CMAX  NOPV0C  A0CMIN  A0CMAX  NOPA0C  R0CMIN  R0CMAX  NOPR0C 
            -110.0  -10.00    21    0.1000  1.5000    29    0.8000  2.0000    25  

CNTKAPMESH  XKCMIN  XKCMAX  NOPXKC  XACMIN  XACMAX  NOPXKC  XRCMIN  XRCMAX  NOPXRC
            -1.000  +1.000    20    -1.000  +1.000    20    -1.000  +1.000    20 

SORLRAMESH  XL_MIN  XL_MAX  NOP_XL  A0SMIN  A0SMAX  NOPA0S  R0SMIN  R0SMAX  NOPR0S      
            05.000  70.000    27    0.1000  1.5000    29    0.4000  2.0000    33

SORKAPMESH  XKSMIN  XKSMAX  NOPXKS  XASMIN  XASMAX  NOPXAS  XRSMIN  XRSMAX  NOPXRS
            -1.000  +1.000    20    -1.000  +1.000    20    -1.000  +1.000    20   
                  
                  Proton parameters
		          =================
		  
IFMESHPARP  IM_V0C  IM_A0C  IM_R0C  IM_XLA  IM_A0S  IM_R0S  IM_XRM  IM_ARM  IM_RRM  IM_COU
              0       0       0       0       0       0       0       0       0       0

IFMESHKAPP  IMKV0C  IMKA0C  IMKR0C  IMKXLA  IMKA0S  IMKR0S  IMKXRM  IMKARM  IMKRRM  IMKCOU
              0       0       0       0       0       0       0       0       0       0  
                    
                  Neutron parameters
		          ==================
		  
IFMESHPARN  IM_V0C  IM_A0C  IM_R0C  IM_XLA  IM_A0S  IM_R0S  IM_XRM  IM_ARM  IM_RRM
              0       0       0       0       0       0       0       0       0

IFMESHKAPN  IMKV0C  IMKA0C  IMKR0C  IMKXLA  IMKA0S  IMKR0S  IMKXRM  IMKARM  IMKRRM
              0       0       0       0       0       0       0       0       0   
            
                  Spin-orbit and tensor lambda (density dependent)
		          ================================================
		  
IFMESHDENS  IM_XPP  IM_XPN  IM_XNP  IM_XNN  
              0       0       0       0      

IFMESHTENS  IM_YPP  IM_YPN  IM_YNP  IM_YNN  
              0       0       0       0 
              
IFMESHCNTN  IM_CPP  IM_CPN  IM_CNP  IM_CNN  
              0       0       0       0   

                  Kappa parametrization (central + spin-orbit new version)
		          =================

IFMESHCNTK  IM_V0C  IM_KVC  IM_R0C  IM_KRC  IM_A0C  IM_KAC
              0       0       0       0       0       0

IFMESHSORK  IM_VSO  IM_KVS  IM_RSO  IM_KRS  IM_ASO  IM_KAS
              0       0       0       0       0       0  

            ==================================================
            ==================================================
            
    
                  Monte-Carlo Section (restarts, bins, mean and sigma values)
		          =================

MONTECARLO  IFPSEU  IFPARA  LDMONT  LDBINS
               0       0    050000    40

                  Mean values for the gaussian restarts
		          =================

CENTPROT_M  V0CENT  R0CENT  A0CENT  XK_V0C  XK_R0C  XK_A0C  R0COUL  XK_COU
            -62.21  1.2136  0.6337  0.0000  0.0000  0.0000  0.0000  0.0000

CENTNEUT_M  V0CENT  R0CENT  A0CENT  XK_V0C  XK_R0C  XK_A0C
            -39.93  1.3607  0.6784  0.0000  0.0000  0.0000

SORBPROT_M  V0SORB  R0SORB  A0SORB  XK_LAM  XK_RSO  XK_ASO 
            24.007  1.0007  0.7985  0.0000  0.0000  0.0000 

SORBNEUT_M  V0SORB  R0SORB  A0SORB  XK_LAM  XK_RSO  XK_ASO   
            26.843  1.2362  0.3722  0.0000  0.0000  0.0000

EFFEPROT_M  V0EFFM  R0EFFM  A0EFFM  XK_LEF  XK_REF  XK_AEF 
            0.0000  0.0000  0.0000  0.0000  0.0000  0.0000

EFFENEUT_M  V0EFFM  R0EFFM  A0EFFM  XK_LEF  XK_REF  XK_AEF  
            0.0000  0.0000  0.0000  0.0000  0.0000  0.0000

SORBDENS_M  ALAMPP  ALAMPN  ALAMNP  ALAMNN
            0.0000  0.0000  0.0000  0.0000

SORBTENS_M  YLAMPP  YLAMPN  YLAMNP  YLAMNN
            0.0000  0.0000  0.0000  0.0000

CENTTENS_M  CLAMPP  CLAMPN  CLAMNP  CLAMNN
            0.0000  0.0000  0.0000  0.0000

KAPACENT_M  V0CENT  XK_V0C  R0CENT  XK_R0C  A0CENT  XK_A0C 
            0.0000  0.0000  0.0000  0.0000  0.0000  0.0000

KAPASORB_M  V0SORB  XK_V0S  R0SORB  XK_R0S  A0SORB  XK_A0S 
            0.0000  0.0000  0.0000  0.0000  0.0000  0.0000

                  Sigma values for the gaussian restarts
		          =================

CENTPROT_S  V0CENT  R0CENT  A0CENT  XK_V0C  XK_R0C  XK_A0C  R0COUL  XK_COU
            0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000

CENTNEUT_S  V0CENT  R0CENT  A0CENT  XK_V0C  XK_R0C  XK_A0C
            0.0000  0.0000  0.0000  0.0000  0.0000  0.0000

SORBPROT_S  V0SORB  R0SORB  A0SORB  XK_LAM  XK_RSO  XK_ASO 
            0.0000  1.0000  0.0000  0.0000  0.0000  0.0000

SORBNEUT_S  V0SORB  R0SORB  A0SORB  XK_LAM  XK_RSO  XK_ASO 
            0.0000  1.0000  0.0000  0.0000  0.0000  0.0000

EFFEPROT_S  V0EFFM  R0EFFM  A0EFFM  XK_LEF  XK_REF  XK_AEF 
            0.0000  0.0000  0.0000  0.0000  0.0000  0.0000

EFFENEUT_S  V0EFFM  R0EFFM  A0EFFM  XK_LEF  XK_REF  XK_AEF  
            0.0000  0.0000  0.0000  0.0000  0.0000  0.0000

SORBDENS_S  ALAMPP  ALAMPN  ALAMNP  ALAMNN
            0.0000  0.0000  0.0000  0.0000

SORBTENS_S  YLAMPP  YLAMPN  YLAMNP  YLAMNN
            0.0000  0.0000  0.0000  0.0000

CENTTENS_S  CLAMPP  CLAMPN  CLAMNP  CLAMNN
            0.0000  0.0000  0.0000  0.0000

KAPACENT_S  V0CENT  XK_V0C  R0CENT  XK_R0C  A0CENT  XK_A0C 
            0.0000  0.0000  0.0000  0.0000  0.0000  0.0000

KAPASORB_S  V0SORB  XK_V0S  R0SORB  XK_R0S  A0SORB  XK_A0S 
            0.0000  0.0000  0.0000  0.0000  0.0000  0.0000

            
            A T T E N T I O N: THIS IS A "THROW AWAY TABLE"
            CONVENTION BELOW: 1=letout 0=keep it in the run  
            Attention Achtung Attention POSITION  SENSITIVE
WHATNOTAKE  N  E  U  T  R  O  N         P  R  O  T  O  N  S
            0  1s1/2   Nshell=0         0  1s1/2   Nshell=0
            0  1p3/2   Nshell=1         0  1p3/2   Nshell=1
            0  1p1/2   Nshell=1         0  1p1/2   Nshell=1
            0  1d5/2   Nshell=2         0  1d5/2   Nshell=2
            0  2s1/2   Nshell=2         0  2s1/2   Nshell=2
            0  1d3/2   Nshell=2         0  1d3/2   Nshell=2
            0  1f7/2   Intruder         0  1f7/2   Intruder
            0  2p3/2   Nshell=3         0  2p3/2   Nshell=3
            0  1f5/2   Nshell=3         0  1f5/2   Nshell=3
            0  2p1/2   Nshell=3         0  2p1/2   Nshell=3
            0  1g9/2   Intruder         0  1g9/2   Intruder
            0  1g7/2   Nshell=4         0  1g7/2   Nshell=4
            0  2d5/2   Nshell=4         0  2d5/2   Nshell=4
            0  2d3/2   Nshell=4         0  2d3/2   Nshell=4
            0  3s1/2   Nshell=4         0  3s1/2   Nshell=4
            0  1h11/2  Intruder         0  1h11/2  Intruder
            0  1h9/2   Nshell=5         0  1h9/2   Nshell=5
            0  2f7/2   Nshell=5         0  2f7/2   Nshell=5
            0  2f5/2   Nshell=5         0  2f5/2   Nshell=5
            0  3p3/2   Nshell=5         0  3p3/2   Nshell=5
            0  3p1/2   Nshell=5         0  3p1/2   Nshell=5
            0  1i13/2  Intruder         0  1i13/2  Intruder 
            0  2g9/2   Nshell=6         0  2g9/2   Nshell=6
            0  3d5/2   Nshell=6         0  3d5/2   Nshell=6
            0  1i11/2  Nshell=6         0  1i11/2  Nshell=6
            0  2g7/2   Nshell=6         0  2g7/2   Nshell=6
            0  4s1/2   Nshell=6         0  4s1/2   Nshell=6
            0  3d3/2   Nshell=6         0  3d3/2   Nshell=6
            0  1j15/2  Intruder         0  1j15/2  Intruder
                         

IDENTIFIER  EXTENT
            no-1s1.1p3.1p1              <<== comment for neutrons    
            in_fact_dummy!              <<== comment for protons  



IFTAKLEVEL  Which levels we consider of each nucleus
            16O
            P1p3/2  P1p1/2  P1d5/2  P2s1/2  P1d3/2
              1       1       1       1       1
            N1p3/2  N1p1/2  N1d5/2  N2s1/2  N1d3/2
              1       1       1       1       1     
            40Ca
            P1d3/2  P2s1/2  P1d5/2  P1f7/2  P2p3/2  P2p1/2  P1f5/2
              1       1       1       1       1       1       1
            N1d3/2  N2s1/2  N1d5/2  N1p1/2  N1f7/2  N2p3/2  N2p1/2  N1f5/2  N1g9/2
              1       1       1       1       1       1       1        1      1
            48Ca
            P2s1/2  P1d3/2  P1d5/2  P1p1/2  P1f7/2  P2p3/2  P2p1/2  P1f5/2  P1g9/2
              1       1       1       1       1       1       1        1      1
            N1d5/2  N2s1/2  N1d3/2  N1f7/2  N2p3/2  N2p1/2  N1f5/2  N1g9/2  N2d5/2
              1       1       1       1       1       1       1        1      1
            56Ni
            P1f7/2  P1d3/2  P2s1/2  P2p3/2  P1f5/2  P2p1/2  P1g9/2
              1       1       1       1       1       1       1   
            N1f7/2  N2s1/2  N1d3/2  N1d5/2  N2p3/2  N1f5/2  N2p1/2  N1g9/2  N2d5/2
              1       1       1       1       1       1       1        1      1
            90Zr
            P2p1/2  P2p3/2  P1f5/2  P1f7/2  P1d3/2  P1g9/2  P2d5/2  P3s1/2  P2d3/2  P1g7/2  P1h11/2
              1       1       1       1       1       1       1        1      1       1       1
            N1g9/2  N2p1/2  N2p3/2  N1f5/2  N1f7/2  N2d5/2  N3s1/2  N2d3/2  N1g7/2  N1h11/2 N2f7/2
              1       1       1       1       1       1       1        1      1       1       1
            132Sn
            P1g9/2  P2p1/2  P2p3/2  P1f5/2  P1f7/2  P1g7/2  P2d5/2  P3s1/2  P1h11/2 P2d3/2  P2f7/2
              1       1       1       1       1       1       1        1      1       1       1
            N2d3/2  N1h11/2 N3s1/2  N2d5/2  N1g7/2  N1g9/2  N2f7/2  N3p3/2  N1h9/2  N3p1/2  N2f5/2  N1i13/2 N3p1/2  N2g9/2
              1       1       1       1       1       1       1        1      1       1       1       1        1       1
            146Gd
            P2d5/2  P1g7/2  P3s1/2  P1h11/2 P2d3/2
              1       1       1       1       1   
            N3s1/2  N2d3/2  N1h11/2 N2f7/2  N1i13/2 N3p3/2  N1h9/2
              1       1       1       1       1       1       1   
            208Pb
            P3s1/2  P2d3/2  P1h11/2 P2d5/2  P1g7/2  P1g9/2  P1h9/2  P2f7/2  P1i13/2 P3p3/2 P2f5/2  P3p1/2  P2g9/2
              1       1       1       1       1       1       1        1      1       1       1       1        1  
            N3p1/2  N2f5/2  N3p3/2  N1i13/2 N2f7/2  N1h9/2  N1h11/2 N2g9/2  N1i11/2 N3d5/2 N1j15/2 N4s1/2  N2g7/2  N3d3/2
              1       1       1       1       1       1       1        1      1       1       1       1        1     1

LEV_LENGTH  XMIN_T  XMAX_T  XMIN_E  XMAX_E
              2.5     3.5    4.3      5.3

ATTRIBUTES  STYLE1  TYPEPN
               0       3
            COLOR1  COLOR2  COLOR3  COLOR4  COLOR5  COLOR6  COLOR7  COLOR8  COLOR9  COLO10  COLO11  COLO12  COLO13  COLO14  COLO15  COLO16  COLO17  COLO18  COLO19  COLO20  COLO21  COLO22  COLO23  COLO24  COLO25  COLO26  COLO27  COLO28  COLO29  COLO30  COLO31  COLO32  COLO33  COLO34  COLO35  COLO36  COLO37  COLO38  COLO39  COLO40
              44      47      53      56      59      62      65      68      45      51      54      58      59      60      61      62      63      64      65      66      67      68       44      45      46      47      51      52      53      54      55      56      57      58      59      60      61      62      63      64
            THICK1  THICK2  THICK3  THICK4  THICK5  THICK6  THICK7  THICK8  THICK9  THIC10  THIC11  THIC12  THIC13  THIC14  THIC15  THIC16  THIC17  THIC18  THIC19  THIC20  THIC21  THIC22  THIC23  THIC24  THIC25  THIC26  THIC27  THIC28  THIC29  THIC30  THIC31  THIC32  THIC33  THIC34  THIC35  THIC36  THIC37  THIC38  THIC39  THIC40
               5       5        5      5       5       5       5       5       5       5       5       5       5       5       5       5       5       5       5       5       5       5       5       5       5       5       5       5       5       5       5       5       5       5       5       5       5       5       5       5
                 
EXECUTE!!!

==========================================================================================
          
                                         

C FILE NAME = wspher12_minpac_15.f ! Keep this symbol:    $ident@string$
C
C=======================================================================
C=======================================================================
C             LEVENBERG-MARQUARDT MINIMISATION ROUTINE PACKAGE   
C=======================================================================
C=======================================================================
C
      SUBROUTINE LEVMAR(FUNMIN,NDFUNC,LDFUNC,NDPARS,LDPARS,XARGUM,
     *                  FARGUM,FJACOB,TOLERF,TOLERX,TOLERG,MAXFEV,
     *                  DIAGSC,IFMODE,FACTOR,NPRINT,INFRUN,NFCALL,
     *                  NJCALL,I_PERM,QTRANF,WORKA1,WORKA2,WORKA3,
     *                                              WORKA4,NDLAST)
C
      INCLUDE   'MATDIM/MAXFEV_AUXILI.f'
      INCLUDE   'MATDIM/NDPARS_STORAG.f'
      INCLUDE   'MATDIM/NDRAND.f'
C      
      DIMENSION 
     *          XARGUM(1:NDPARS),FARGUM(1:NDFUNC),
     *          DIAGSC(1:NDPARS),QTRANF(1:NDPARS) 
      DIMENSION 
     *          WORKA1(1:NDPARS),WORKA2(1:NDPARS),
     *          WORKA3(1:NDPARS),WORKA4(1:NDFUNC)
      DIMENSION
     *          FJACOB(1:NDFUNC,1:NDPARS) 
      DIMENSION 
     *          I_PERM(1:NDPARS)
C
      COMMON
     *       /STOPAR/ EPSLAS,LDLAST
CID      COMMON
C     *      /COLECT_GRADMI/ GRADMI_ITERAT(1:MAXFEV_AUXILI,1:NDRAND)
C     *      /COLECT_PARSMI/ PARSMI_ITERAT(1:MAXFEV_AUXILI,1:NDRAND,
C     *                                    1:NDPARS_STORAG)
      COMMON
     *       /PRFILE/ ISCREN,LOGFIL,LOGCHP,LOGCHN,LOGHIP,LOGHIN,
     *                LOGPRO,LOGNEU,LOGWRI,LSCREN,IRESUL,ICONVE,
     *                NOUTPT,LOGAUX,LOGENE,LOGRAD,LOGGRD,LOGMSH,
     *                                                   LOGBIS  
      COMMON
     *       /RANRAN/ IRMFST,I_RAND,IRANDO
      COMMON 
     *       /LAST_G/ GRAD_N                                                 
C
      DATA 
     *       R_UNIT, CONST1, CONST5, CONS25, CONS75, CN0001, R_ZERO
     *     / 1.0000, 1.0E-1, 5.0E-1, 2.5E-1, 7.5E-1, 1.0E-4, 0.0000/
C
C=======================================================================
C
C     The purpose of LEVMAR is to minimise the sum of the squares of
C     LDFUNC nonlinear functions of LDPARS variables  by employing a
C     modification  of the standard Levenberg-Marquardt method.  The 
C     user must provide a subroutine, which calculates the functions 
C     and the Jacobian.
C
C     The subroutine statement is
C
C     SUBROUTINE LEVMAR(FUNMIN,MDFUNC,LDFUNC,NDPARS,LDPARS,XARGUM,
C                       FARGUM,FJACOB,TOLERF,TOLERX,TOLERG,MAXFEV,
C                       DIAGSC,IFMODE,FACTOR,NPRINT,INFRUN,NFCALL,
C                       NJCALL,I_PERM,QTRANF,WORKA1,WORKA2,WORKA3,
C                                                          WORKA4)
C
C     where
C
C     FUNMIN - the name of the user-supplied subroutine calculating 
C              the functions and the Jacobian.  It must be declared 
C              in an external statement in the user calling program, 
C              and should be written as follows.
C
C              SUBROUTINE FUNMIN(NDFUNC,LDFUNC,NDPARS,LDPARS,XARGUM,
C                                              FARGUM,FJACOB,I_FLAG)
C
C              DIMENSION 
C                        XARGUM(1:NDPARS),FARGUM(1:NDFUNC),
C                        FJACOB(1:NDFUNC,1:NDPARS)
C              ---------
C              IF I_FLAG = 1 calculate the functions at XARGUM, and
C                          return this vector in FARGUM. 
C                          Do not alter FJACOB.
C              IF I_FLAG = 2 calculate  the corresponding  Jacobian 
C                          at XARGUM,  return this matrix in FJACOB. 
C                          Do not alter FARGUM.
C              ---------
C              RETURN
C              END
C
C     The value of I_FLAG should not be changed by FUNMIN unless the
C     user wants to terminate execution of LEVMAR.  In this case set 
C     I_FLAG to a negative integer.
C
C     LDFUNC - a positive integer input variable;  set to the number
C              of functions.
C
C     LDPARS - a positive integer input variable;  set to the number
C              of variables.  LDPARS must not exceed NDPARS, and the 
C              latter must not exceed LDFUNC.
C
C     XARGUM - a vector of length NDPARS.  On input, it must contain
C              an initial estimate of the solution vector. On output 
C              XARGUM will contain the final result for the solution 
C                                                             vector.
C
C     FARGUM - an output array of length NDFUNC,  which contains the 
C              functions evaluated at the output XARGUM.
C
C     FJACOB - an output NDFUNC by NDPARS array. The upper LDPARS by 
C              LDPARS submatrix of the square matrix FJACOB contains 
C              an upper-triangular matrix  RUPTRI;  it has  diagonal 
C              elements of non-increasing magnitude such that
C
C                        T          T                T
C                       P *(Jacobian *Jacobian)*P = R *R,
C
C              where P is a permutation matrix, and Jacobian  is the 
C              final calculated Jacobian.  Column `J' of P is column 
C              I_PERM(J), (below), of the identity matrix. The lower 
C              trapezoidal part of FJACOB contains information gene-
C              rated during the computation of RUPTRI.
C
C     NDFUNC - the maximum number of the functions, that can be used
C              for the minimisation algorithm
C
C     LDFUNC - a positive integer, that must not exceed NDFUNC,  the 
C              actual number of functions, and the leading dimension 
C              of the array FJACOB
C
C     TOLERF - a nonnegative input variable. Termination occurs when 
C              both the actual, and predicted relative reductions in 
C              the sum of squares are at most TOLERF. Therefore, the
C              latter measures the relative error desired in the sum 
C              of squares.
C
C     TOLERX - a nonnegative input variable. Termination occurs when 
C              the relative error  between  two consecutive iterates 
C              is at most TOLERX. Therefore, the latter measures the
C              relative error desired in the approximate solution.
C
C     TOLERG - a nonnegative input variable. Termination occurs when 
C              the cosine of the angle between FARGUM and any column 
C              of the Jacobian is at most  TOLERG  in absolute value. 
C              Therefore, TOLERG  measures the orthogonality desired 
C              between the function vector and the resulting columns 
C              of the Jacobian.
C
C     MAXFEV - a positive integer input variable. Termination occurs 
C              when the number of calls to FUNMIN  with I_FLAG=1 has
C              reached MAXFEV.
C
C     DIAGSC - a vector of length NDPARS. If IFMODE = 1  (see below), 
C              DIAGSC is internally set.  If IFMODE = 2, DIAGSC must 
C              contain positive entries that serve as multiplicative 
C              scale factors for the variables
C
C     IFMODE - an integer input variable. If IFMODE = 1,   variables 
C              will be scaled internally. If IFMODE = 2, the scaling 
C              is specified by the input DIAGSC. Other values of pa-
C              rameter IFMODE are equivalent to IFMODE = 1.
C
C     FACTOR - a positive input variable that determines the initial 
C              step bound.  This bound is defined as  the product of
C              FACTOR and the Euclidean norm of DIAGSC*XARGUM if non-
C              zero, or else to FACTOR itself.  In most cases FACTOR 
C              should lie in the interval  [0.1, 100.0]. This  upper 
C              value (100.0) is a generally recommended.
C
C     NPRINT - an integer input variable that is used to control the
C              printing of iterates if it is positive.  In this case,
C              FUNMIN is called with I_FLAG = 0  at the beginning of 
C              the first iteration and every NPRINT iterations there-
C              after - and immediately prior to return,  with XARGUM, 
C              FARGUM, and FJACOB available for printing. FARGUM and
C              FJACOB should not be altered. If  NPRINT  is negative, 
C              no special calls of FUNMIN with I_FLAG = 0 are made.
C
C     INFRUN - an integer output variable.  If the user  has stopped 
C              execution,  INFRUN is set  to the (negative) value of 
C              I_FLAG. see description of FUNMIN.  Otherwise, INFRUN 
C              is set as follows:
C
C              INFRUN = 0  Improper input parameters
C
C              INFRUN = 1  Both the actual  and  predicted  relative  
C                          reductions  in the sum  of squares are at 
C                          most TOLERF
C
C              INFRUN = 2  Relative error  between  two  consecutive 
C                          iterates is at most TOLERX
C
C              INFRUN = 3  Conditions for INFRUN=1 and INFRUN=2 both 
C                                                               hold
C
C              INFRUN = 4  The cosine  of the angle  between  FARGUM 
C                          and any column of the Jacobian is at most 
C                          TOLERG in absolute value
C
C              INFRUN = 5  Number of calls to FUNMIN with I_FLAG = 1 
C                          has reached MAXFEV
C
C              INFRUN = 6  TOLERF is too small. No further reduction 
C                          in the sum of squares is possible
C
C              INFRUN = 7  TOLERX is too small.  In such  a case  no 
C                          further  improvement  in  the approximate 
C                          solution XARGUM is possible
C
C              INFRUN = 8  TOLERG is too small. FARGUM is orthogonal 
C                          to the columns of the Jacobian to machine 
C                          precision
C
C     The following ones are introduced by JD and ID:
C
C              INFRUN = 9  The difference between  the  last  LDLAST
C                          chi^2 is less than EPSLAS
C
C              INFRUN = 10 The norm of the real gradien is less than
C                          TOLERG
C
C     NFCALL - an integer output variable set to the number of calls 
C              to FUNMIN with I_FLAG = 1
C
C     NJCALL - an integer output variable set to the number of calls 
C              to FUNMIN with I_FLAG = 2
C
C     I_PERM - an integer output array of length NDPARS.  It defines 
C              a permutation matrix P such that Jacobian * P = Q * R,
C              where  "Jacobian"  stands for  the finally calculated 
C              Jacobian, Q is orthogonal (not stored), and R denoted
C              RUPTRI in the code is upper triangular, with diagonal 
C              elements of nonincreasing magnitude. Column J of P is 
C              column I_PERM(J) of the identity matrix.
C
C     QTRANF - an output array of length NDPARS,  which contains the 
C              first LDPARS elements  of the vector  obtained as the
C              result of multiplication (Q Transpose)*FARGUM.
C
C     WORKA1,  WORKA2, and WORKA3 are auxiliary work arrays of length 
C              NDPARS
C
C     WORKA4 - a work array of length NDFUNC.
C
C     Subprograms called:
C
C              User-supplied ...... FUNMIN
C
C              MINPACK-supplied ... DPMPAR,EUNORM,LM_PAR,QRFACT
C
C              FORTRAN-supplied ... ABS,MAX,MIN,SQRT,MOD
C
C     Argonne National Laboratory. MINPACK PROJECT. March 1980.
C     Burton S. Garbow, Kenneth E. Hillstrom and Jorge J. More
C
C=======================================================================
C
      CALL CPUTIM('LEVMAR',1)
C
C=======================================================================
C
      IF (NDPARS.NE.NDPARS_STORAG) THEN
          WRITE(NOUTPT,'(''Please make NDPARS='',I3,'' coincide with '',
     *                   ''NDPARS_STORAG='',I3)') NDPARS,
     *                     NDPARS_STORAG
          STOP 'NDPARS must be equal to NDPARS_STORAG, STOP from LEVMAR'
      END IF
C
C=======================================================================
C
      IF (MAXFEV.NE.MAXFEV_AUXILI) THEN
          WRITE(0,'(/,''Alarm in LEVMAR: MAXFEV= '',I4,
     *              '' and MAXFEV_AUXILI= '',I4,'' are not equal!'',/)')
     *                 MAXFEV,MAXFEV_AUXILI
          STOP 'STOP in LEVMAR: MAXFEV.NE.MAXFEV_AUXILI'
      END IF
C
C=======================================================================
C
CID      DO I=1,MAXFEV
CID         DO J=1,NDRAND
CID            GRADMI_ITERAT(I,J)=999999.999
CID            DO K=1, NDPARS
CID               PARSMI_ITERAT(I,J,K)=99999.999
CID            END DO
CID         END DO
CID      END DO
C
C=======================================================================
C
      EPSMCH=DPMPAR(1) !  EPSMCH is the machine precision
C
      I_STEP_GRAD_F=0
C
      I_ZERO=0
      INFRUN=0
      I_FLAG=0
C
      NFCALL=0 ! Number of function evaluations
      NJCALL=0 ! Number of Jacobian evaluations
C
C     Check the input parameters for inconsistency/errors
C
      IF (LOGWRI.GT.0) THEN
C
          WRITE(LOGFIL,'(/,09X,''Inside LEVMAR: Verification of '',
     *                         ''the input'',
     *                   /,09X,''We should have '',6x,
     *                             ''   Actually we have'')')
C_______________________________________________________________________
C
          WRITE(LOGFIL,'(09X,''LDPARS .GT. I_ZERO:'',9x,I3,
     *                                     '' and'',I3,'' ?'')')
     *                         LDPARS,     I_ZERO
C
          WRITE(LOGFIL,'(09X,''LDFUNC .GE. LDPARS:'',9x,I3,
     *                                     '' and'',I3,'' ?'')')
     *                         LDFUNC,     LDPARS
C
          WRITE(LOGFIL,'(09X,''MAXFEV .GT. I_ZERO:'',9x,I3,
     *                                   '' and'',I3,'' ?'',/)')
     *                         MAXFEV,     I_ZERO
C_______________________________________________________________________
C
          WRITE(LOGFIL,'(09X,''TOLERF .GE. R_ZERO: '',3x,E8.2,
     *                                     '' and '',E8.2,'' ?'')')
     *                         TOLERF,     R_ZERO
C
          WRITE(LOGFIL,'(09X,''TOLERX .GE. R_ZERO: '',3x,E8.2,
     *                                     '' and '',E8.2,'' ?'')')
     *                         TOLERX,     R_ZERO
C
          WRITE(LOGFIL,'(09X,''TOLERG .GE. R_ZERO: '',3x,E8.2,
     *                                     '' and '',E8.2,'' ?'')')
     *                         TOLERG,     R_ZERO
C
          WRITE(LOGFIL,'(09X,''FACTOR .GT. R_ZERO: '',3x,E8.2,
     *                                     '' and '',E8.2,'' ?'')')
     *                         FACTOR,     R_ZERO
      END IF
C
C @@@ WHAT IS THIS FORM OF DANGEROUS PROGRAMMING?
C @@@ IRENE - PLEASE INTRODUCE STOP STOP STOP IN THE CASE OF COMMENTED IF
C
C     IF (LDFUNC.LT.LDPARS) THEN
C
C         WRITE(NOUTPT,'(''Nonsense request: Number of parameters '',
C    *                   ''exceeds the number of data points'',/)')
C
C         STOP 'Number of parameters > number of data points, in LEVMAR'
C
C     END IF
C
CID:  We have 'silenced' on of the conditions in the followinf IF
C     in order to be able to fit all the cases without the radius
C     We know that this will probably mean the appearance of over
C     fitting in our results. And this is what we want to show and
C     we will show how to fix it.
C
      IF (LDPARS.LE.I_ZERO .OR. !LDFUNC.LT.LDPARS.OR.
     *    TOLERF.LT.R_ZERO .OR. TOLERX.LT.R_ZERO.OR.
     *    TOLERG.LT.R_ZERO .OR. MAXFEV.LE.I_ZERO.OR.
     *                          FACTOR.LE.R_ZERO) THEN
C_______________________________________________________________________
C
          IF (LOGWRI.GT.0) THEN
              WRITE(LOGFIL,'(9X,''One of the above conditions was '',
     *                          ''true =>> terminating and exiting'')')
          END IF
C
          GO TO 300
C
      END IF
C
C=======================================================================
C
      IF (IFMODE.EQ.2) THEN
          DO J=1,LDPARS
             IF (DIAGSC(J).LE.R_ZERO) GO TO 300
          END DO
      END IF
C
C     We evaluate the function at the starting point and calculate 
C                                                        its norm.
      IF (LOGWRI.GT.0) THEN
          WRITE(LOGFIL,'(/,9X,''Entering FUNMIN from LEVMAR [1] '',
     *                        '' starting point'')')
      END IF 
C 
      I_FLAG=1
C
      CALL FUNMIN(NDFUNC,LDFUNC,NDPARS,LDPARS,XARGUM,FARGUM,
     *                                        FJACOB,I_FLAG)
C
C=======================================================================
C
      NFCALL=1
C
      IF (I_FLAG.LT.0) GO TO 300
C
      F_NORM=EUNORM(NDFUNC,LDFUNC,FARGUM)
C
C     Initialise Levenberg-Marquardt parameter and iteration counter
C
      PARAUX=R_ZERO
      ITERAT=1
C
C     Beginning of the outer loop <<<== jumping here from below ...
C
   30 CONTINUE
C
      I_FLAG=2 !  Calculate the Jacobian matrix
C
      IF (LOGWRI.GT.0) THEN
          WRITE(LOGFIL,'(/,9X,''Entering FUNMIN from LEVMAR [2] '',
     *                        ''to calculate the Jacobian'')')
      END IF 
C
      CALL FUNMIN(NDFUNC,LDFUNC,NDPARS,LDPARS,XARGUM,FARGUM,
     *                                        FJACOB,I_FLAG)
C
C=======================================================================
C
C     Calculating the 'real' \chi^2 gradient norm, namely:
C
C     \chi^2 = \sum_j [ f_j(x) - y_j ] ^ 2
C
C     d\chi^2/dx_i = 2 * \sum_j ( f_j(x) - y_j ) * df_j(x_i)/dx_i
C
C     norm = sqrt[ \sum_i (d\chi^2/dx_i) ^ 2 ]
C
C     From FUNMIN we have:
C
C     ( f_j(x) - y_j ) = FARGUM(J_FUNC)
C
C     df_j(x_i)/dx_i = FJACOB(J_FUNC,I_PARS)
C
      GRAD_N=0.0
C
      DO I_PARS=1,LDPARS
C
         GRAD_V=0.0
C
         DO I_FUNC=1,LDFUNC
            GRAD_V=GRAD_V+FARGUM(I_FUNC)*FJACOB(I_FUNC,I_PARS)
         END DO
C
         GRAD_V=2*GRAD_V
         GRAD_N=GRAD_N+GRAD_V**2
C
      END DO
C
      GRAD_N=SQRT(GRAD_N)
C
C=======================================================================
C
      IF (LOGWRI.GT.0) THEN
C
          I_STEP_GRAD_F=I_STEP_GRAD_F+1
C
CID          GRADMI_ITERAT(I_STEP_GRAD_F,IRANDO)=GRAD_N
C
          WRITE(LOGFIL,'()')
C
          DO I=1,3
             WRITE(LOGFIL,'(9X,39(''#''),''   Evaluation number='',i3,
     *                                         1X,''for IRANDO='',I2)')
     *                                        I_STEP_GRAD_F,IRANDO
          END DO
C
CID          WRITE(LOGFIL,'(/,9X,''In LEVMAR: GRAD_N= '',f20.13,/,20x,
CID     *                        ''TOLERG= '',f20.13,/)') GRAD_N,TOLERG
          DO I=1,3
             WRITE(LOGFIL,'(9X,39(''#''),''   Evaluation number='',i3,
     *                                         1X,''for IRANDO='',I2)')
     *                                        I_STEP_GRAD_F,IRANDO
          END DO
C
      END IF
C
C=======================================================================
C
C     Checking if the gradient is less than our tolerance
C
      IF (GRAD_N.LT.TOLERG) INFRUN=10
C
      IF (INFRUN.NE.0) GO TO 300
C
C=======================================================================
C
      NJCALL=NJCALL+1
C
      IF (I_FLAG.LT.0) GO TO 300
C
C=======================================================================
C
C     If requested, call FUNMIN to enable printing of iterates
C
      IF (NPRINT.GT.0) THEN
C
          I_FLAG=0
C
          IF (MOD(ITERAT-1,NPRINT).EQ.0) THEN
C
              IF (LOGWRI.GT.0) THEN
                  WRITE(LOGFIL,'(9X,''Entering FUNMIN from LEVMAR [3]'',
     *                          '' to enable printing iterates, '',
     *                          ''as they say...'')')
              END IF 
C
              CALL FUNMIN(NDFUNC,LDFUNC,NDPARS,LDPARS,XARGUM,FARGUM,
     *                                                FJACOB,I_FLAG)
          END IF
C      
          IF (I_FLAG.LT.0) GO TO 300
C
      END IF
C
C=======================================================================
C
C     Compute the QR factorisation of the Jacobian
C
      CALL QRFACT(NDFUNC,LDFUNC,NDPARS,LDPARS,FJACOB,.TRUE.,
     *                          I_PERM,WORKA1,WORKA2,WORKA3)
C
C=======================================================================
C
C     On the first iteration and if IFMODE=1, scale according
C     to the norms of the columns of the initial Jacobian.
C
      IF (ITERAT.EQ.1) THEN
C
          IF (IFMODE.NE.2) THEN
C
              DO J=1,LDPARS
                                          DIAGSC(J)=WORKA2(J)
                 IF (WORKA2(J).EQ.R_ZERO) DIAGSC(J)=R_UNIT
              END DO
C
          END IF
C
C         On the first iteration, calculate the norm of the scaled 
C         XARGUM and initialise the step bound DELTAB
C
          DO J=1,LDPARS
             WORKA3(J)=DIAGSC(J)*XARGUM(J)
          END DO
C 
          XNORM=EUNORM(NDPARS,LDPARS,WORKA3)
          DELTAB=FACTOR*XNORM
C 
          IF (DELTAB.EQ.R_ZERO) DELTAB=FACTOR
C 
      END IF
C
C=======================================================================
C
C     Form (Q TRANSPOSE)*FARGUM and store the first LDPARS components 
C                                                           in QTRANF 
      DO I=1,LDFUNC
         WORKA4(I)=FARGUM(I)
      END DO
C
      DO J=1,LDPARS
C
         IF (FJACOB(J,J).NE.R_ZERO) THEN
C
             SUMAUX=R_ZERO
C
             DO I=J,LDFUNC
                SUMAUX=SUMAUX+FJACOB(I,J)*WORKA4(I)
             END DO
C
             TEMPOR=-SUMAUX/FJACOB(J,J)
C
             DO I=J,LDFUNC
                WORKA4(I)=WORKA4(I)+FJACOB(I,J)*TEMPOR
             END DO
C
         END IF
C
         FJACOB(J,J)=WORKA1(J)
         QTRANF(J)=WORKA4(J)
C
      END DO
C
C=======================================================================
C
C     Compute the norm of the scaled gradient
C
      G_NORM=R_ZERO
C
      IF (F_NORM.NE.R_ZERO) THEN
C
          DO J=1,LDPARS
C
             L=I_PERM(J)
C
             IF (WORKA2(L).NE.R_ZERO) THEN
C
                 SUMAUX=R_ZERO
                 DO I=1,J
                    SUMAUX=SUMAUX+FJACOB(I,J)*(QTRANF(I)/F_NORM)
                 END DO
C
                 G_NORM=MAX(G_NORM,ABS(SUMAUX/WORKA2(L)))
C
             END IF
C
          END DO
C
      END IF

      IF (LOGWRI.GT.0) THEN
          WRITE(LOGFIL,'(/,9X,''In LEVMAR: G_NORM= '',f20.13,/,20x,
     *                        ''TOLERG= '',f20.13,/)') G_NORM,TOLERG
      END IF
C
C=======================================================================
C
C     Test for convergence of the gradient norm
C
      IF (G_NORM.LE.TOLERG) INFRUN=4
C
      IF (INFRUN.NE.0) GO TO 300
C
C
C=======================================================================
C
C     Rescale if necessary
C
      IF (IFMODE.NE.2) THEN
          DO J=1,LDPARS
             DIAGSC(J) =MAX(DIAGSC(J),WORKA2(J))
          END DO
      END IF
C
C=======================================================================
C
C     Beginning of the inner loop
C
  200 CONTINUE
C
C=======================================================================
C
C     Determine the Levenberg-Marquardt parameter
C
      CALL LM_PAR(NDFUNC,LDFUNC,NDPARS,LDPARS,FJACOB,I_PERM,
     *            DIAGSC,QTRANF,DELTAB,PARAUX,WORKA1,WORKA2,
     *                                        WORKA3,WORKA4)
C
C=======================================================================
C
C     Store the direction P and XARGUM+P. Calculate the norm of P.
C
      DO J=1,LDPARS
         WORKA1(J)=-WORKA1(J)
         WORKA2(J)=XARGUM(J)+WORKA1(J)
         WORKA3(J)=DIAGSC(J)*WORKA1(J)
      END DO
C
      P_NORM=EUNORM(NDPARS,LDPARS,WORKA3)
C
C=======================================================================
C
C     On the first iteration, adjust the initial step bound
C
      IF (ITERAT.EQ.1) DELTAB=MIN(DELTAB,P_NORM)
C
      I_FLAG=1 !  Evaluate the function at XARGUM + P
C                              and calculate its norm
C
      IF (LOGWRI.GT.0) THEN
          WRITE(LOGFIL,'(9X,''Entering FUNMIN from LEVMAR [4] '',
     *                      ''to evaluate the functions at p+dp'')')
      END IF 
C
      CALL FUNMIN(NDFUNC,LDFUNC,NDPARS,LDPARS,WORKA2,WORKA4,
     *                                        FJACOB,I_FLAG)
C
      NFCALL=NFCALL+1
C
      IF (I_FLAG.LT.0) GO TO 300
C
      FNORM1=EUNORM(NDFUNC,LDFUNC,WORKA4)
C
C=======================================================================
C
C     Compute the scaled actual reduction
C
      ACTRED=-R_UNIT
C
      IF (CONST1*FNORM1.LT.F_NORM) THEN
          ACTRED=R_UNIT-(FNORM1/F_NORM)**2
      END IF
C
C=======================================================================
C
C     Compute the scaled predicted reduction and the scaled 
C                                    directional derivative
      DO J=1,LDPARS
         WORKA3(J)=R_ZERO
         L=I_PERM(J)
         TEMPOR=WORKA1(L)
         DO I=1,J
            WORKA3(I)=WORKA3(I)+FJACOB(I,J)*TEMPOR
         END DO
      END DO
C
      TEMP_1=EUNORM(NDPARS,LDPARS,WORKA3)/F_NORM
      TEMP_2=(SQRT(PARAUX)*P_NORM)/F_NORM
C
      PRERED=TEMP_1**2+TEMP_2**2/CONST5
      DIRDER=-(TEMP_1**2+TEMP_2**2)
C
C=======================================================================
C
C     Compute the ratio of the actual to the predicted reduction
C
                            RATIO=R_ZERO
      IF (PRERED.NE.R_ZERO) RATIO=ACTRED/PRERED
C
C=======================================================================
C
C     Update the step bound
C
      IF (RATIO.LE.CONS25) THEN
C
          IF (ACTRED.GE.R_ZERO) TEMPOR=CONST5
C
          IF (ACTRED.LT.R_ZERO) THEN
              TEMPOR=CONST5*DIRDER/(DIRDER+CONST5*ACTRED)
          END IF
C
          IF (CONST1*FNORM1.GE.F_NORM .OR. TEMPOR.LT.CONST1) THEN
              TEMPOR=CONST1
          END IF
C
          DELTAB=TEMPOR*MIN(DELTAB,P_NORM/CONST1)
          PARAUX=PARAUX/TEMPOR
C
          GO TO 260
C
      END IF
C
      IF (.NOT.(PARAUX.NE.R_ZERO.AND.RATIO.LT.CONS75)) THEN
          DELTAB=P_NORM/CONST5
          PARAUX=CONST5*PARAUX
      END IF
C
  260 CONTINUE
C
C=======================================================================
C
C     Test for successful iteration
C
      IF (RATIO.GE.CN0001) THEN
C
C         Successful iteration. update XARGUM, FARGUM
C                                     and their norms
          DO J=1,LDPARS
             XARGUM(J)=WORKA2(J)
             WORKA2(J)=DIAGSC(J)*XARGUM(J)
          END DO
C
          DO I=1,LDFUNC
             FARGUM(I)=WORKA4(I)
          END DO
C
          XNORM=EUNORM(NDPARS,LDPARS,WORKA2)
          F_NORM=FNORM1
          ITERAT=ITERAT+1
C
      END IF
C
C=======================================================================
C
C     Tests for convergence
C
      IF (ABS(ACTRED).LE.TOLERF .AND. PRERED.LE.TOLERF
     *                          .AND. CONST5*RATIO.LE.R_UNIT) INFRUN=1
C
      IF (DELTAB.LE.TOLERX*XNORM) INFRUN=2
C
      IF (ABS(ACTRED).LE.TOLERF .AND. PRERED.LE.TOLERF .AND.
     *              CONST5*RATIO.LE.R_UNIT .AND. INFRUN.EQ.2) INFRUN=3
C
C=======================================================================
C
      IF (LOGWRI.GT.0) THEN
          WRITE(LOGFIL,'(/,9X,''Entering VERIFY_STOPIN from LEVMAR'')')
      END IF
C
      CALL VERIFY_STOPIN(NDLAST,I_STOP)
      
      IF (I_STOP.EQ.LDLAST-1) INFRUN=9
C
      IF (INFRUN.NE.0) GO TO 300
C
C=======================================================================
C
C     Tests for termination and stringent tolerances
C
      IF (NFCALL.GE.MAXFEV) INFRUN=5
C 
      IF (ABS(ACTRED).LE.EPSMCH .AND. PRERED.LE.EPSMCH.AND.
     *                          CONST5*RATIO.LE.R_UNIT) INFRUN=6
C
      IF (DELTAB.LE.EPSMCH*XNORM) INFRUN=7
C
      IF (G_NORM.LE.EPSMCH) INFRUN=8
C
      IF (INFRUN.NE.0) GO TO 300
C
C         End of the inner loop. Repeat if iteration unsuccessful
C                                                   (go up to 200)
          IF (RATIO.LT.CN0001) GO TO 200
C
C         End of the outer loop
C
      GO TO 30
C
  300 CONTINUE
C
C=======================================================================
C
C     Termination, either normal or user imposed.
C
      IF (I_FLAG.LT.0) THEN 
          INFRUN=I_FLAG
      END IF
C
      I_FLAG=0
C
      IF (NPRINT.GT.0) THEN
C
          IF (LOGWRI.GT.0) THEN
              WRITE(LOGFIL,'(9X,''Entering FUNMIN from LEVMAR [5] '',
     *                           '' for final print reasons'')')
          END IF 
C
          CALL FUNMIN(NDFUNC,LDFUNC,NDPARS,LDPARS,XARGUM,FARGUM,
     *                                            FJACOB,I_FLAG)
      END IF
C
C=======================================================================
C
      IF (LOGWRI.GT.0) THEN
          WRITE(LOGFIL,'(/,9X,''Exiting  LEVMAR'')')
      END IF
C
C=======================================================================
C      
      CALL CPUTIM('LEVMAR',0)
C
C=======================================================================
C
      RETURN
      END     
C
C=======================================================================
C=======================================================================
C
      SUBROUTINE VERIFY_STOPIN(NDLAST,ICOUNT)
C      
      INCLUDE  'MATDIM/MAXFEV.f'
C
      DIMENSION
     *          HILAST(1:NDLAST),
     *          DIFLAS(1:NDLAST)
C_______________________________________________________________________
C      
      COMMON
     *       /STOPAR/ EPSLAS,LDLAST
      COMMON
     *       /STOPIN/ HIAUXI(1:MAXFEV)
      COMMON
     *       /MINITE/ IDEFCN,ITECHI
      COMMON
     *       /PRFILE/ ISCREN,LOGFIL,LOGCHP,LOGCHN,LOGHIP,LOGHIN,
     *                LOGPRO,LOGNEU,LOGWRI,LSCREN,IRESUL,ICONVE,
     *                NOUTPT,LOGAUX,LOGENE,LOGRAD,LOGGRD,LOGMSH,
     *                                                   LOGBIS 
C
C=======================================================================
C @@@ IRENE, PLEASE ADD COMMENTS - WHAT IS THIS ROUTINE SUPPOSED TO DO ???     
C @@@ IRENE, WHAT IS EPSLAS??? HOW MUCH? WHY?
C=======================================================================
C      
      I_LAST=IDEFCN-1
C
      IF (I_LAST.LT.LDLAST) THEN
          IF (LOGWRI.GT.0) THEN
              WRITE(LOGFIL,'(12X,''I_LAST='',I3,'' < LDLAST='',I3)')
     *                             I_LAST,           LDLAST
              WRITE(LOGFIL,'(09X,''Exiting  VERIFY_STOPIN with no '',
     *                           ''action'')')
          END IF
          RETURN
      END IF
C
C=======================================================================
C          
      CALL CPUTIM('VERSTP',1)
C
C=======================================================================
C      
      IF (NDLAST.LT.LDLAST) THEN
          WRITE(0,'(/,''Alarm in VERIFY_STOPIN!! NDLAST= '',I2,
     *              '' and LDLAST= '',I2,/)')NDLAST,LDLAST
          STOP 'STOP in VERIFY_SOTPIN: NDLAST.LT.LDLAST'
      END IF
C
C=======================================================================
C      
      DO I=1,LDLAST
          HILAST(I)=0.0
      END DO
C
      HILAST(1)=HIAUXI(I_LAST)
C
      DO I=2,LDLAST
         J=I-1
         HILAST(I)=HIAUXI(I_LAST-J)
      END DO
C
      IF (LOGWRI.GT.0) THEN
          WRITE(LOGFIL,'(12X,''Differences: '',$)')
      END IF
C
      DO I=2,LDLAST
          DIFLAS(I-1)=ABS(HILAST(I-1)-HILAST(I))
          IF (LOGWRI.GT.0) THEN
              IF (I.EQ.2) WRITE(LOGFIL,'(    I3,F8.4)')I-1,DIFLAS(I-1)
              IF (I.GT.2) WRITE(LOGFIL,'(15X,I1,F8.4)')I-1,DIFLAS(I-1)
          END IF
      END DO
C
      ICOUNT=0
      DO I=1,LDLAST-1
          IF (DIFLAS(I).LT.EPSLAS) THEN
              ICOUNT=ICOUNT+1
          END IF
      END DO
C
C=======================================================================
C
      IF (LOGWRI.GT.0) THEN
          WRITE(LOGFIL,'(9X,''Exiting  VERIFY_STOPIN with ICOUNT='',
     *                                              I2)') ICOUNT 
      END IF
C
C=======================================================================
C                          
      CALL CPUTIM('VERSTP',0)
C
C=======================================================================
C      
      RETURN
      END
C
C=======================================================================
C=======================================================================
C 
      SUBROUTINE TELLIT(INFOER,TOLERF,TOLERX,TOLERG,MAXFEV,
     *                                       LDLAST,EPSLAS) 
C
      COMMON 
     *       /LAST_G/ GRAD_N                                   
      COMMON
     *       /PRFILE/ ISCREN,LOGFIL,LOGCHP,LOGCHN,LOGHIP,LOGHIN,
     *                LOGPRO,LOGNEU,LOGWRI,LSCREN,IRESUL,ICONVE,
     *                NOUTPT,LOGAUX,LOGENE,LOGRAD,LOGGRD,LOGMSH,
     *                                                   LOGBIS  
C
C=======================================================================
C     Routine printing the interpretation of the minimisation result
C=======================================================================
C
C     Printing the information on screen
C
      IF (ISCREN.GT.0) THEN
C      
          WRITE(LSCREN,'()')
C
          IF (INFOER.EQ.0) 
     *    WRITE(LSCREN,'(''INFOER=0 ==>> Minimisation failed -'',
     *                   '' improper input parameters'')')
C
          IF (INFOER.EQ.1) 
     *    WRITE(LSCREN,'(''INFOER=1 ==>> Minimisation finished - both'',
     *                 '' actual and predicted relative reductions'',
     *                 '' in the sum of squares are at most ftol'',
     *                 '' (FTOL='',E12.5,'')'',
     *                 '' The gradient norm is GRAD_N= '',f10.4,
     *                 '' TOLERG= '',es12.5)') TOLERF,GRAD_N,TOLERG
C
          IF (INFOER.EQ.2) 
     *    WRITE(LSCREN,'(''INFOER=2 ==>> Minimisation finished -'',
     *    '' relative error between two consecutive iterates is at'',
     *    '' most xtol (XTOL='',ES12.5,'')'',
     *    '' The gradient norm is GRAD_N= '',f10.4,
     *    '' TOLERG= '',es12.5)') TOLERX,GRAD_N,TOLERG
C
          IF (INFOER.EQ.3) 
     *    WRITE(LSCREN,'(''INFOER=3 ==>> Minimisation finished - both'',
     *    '' actual and predicted relative reductions'',
     *    '' in the sum of squares are at most ftol AND relative'',
     *    '' error between two consecutive iterates is at most xtol'',
     *    '' (FTOL='',ES12.5,'')'','' (XTOL='',ES12.5,'')'',
     *    '' The gradient norm is GRAD_N= '',f10.4,
     *    '' TOLERG= '',es12.5)')TOLERF,TOLERX,GRAD_N,TOLERG
C
          IF (INFOER.EQ.4) 
     *    WRITE(LSCREN,'(''INFOER=4 ==>> Minimisation finished - the'',
     *    '' cosine of the angle between fvec and any column of the'',
     *    '' jacobian is at most gtol in absolute value'',
     *    '' (GTOL='',ES12.5,'')'',
     *    '' The gradient norm is GRAD_N= '',f10.4,
     *    '' TOLERG= '',es12.5)')TOLERG,GRAD_N,TOLERG
C
          IF (INFOER.EQ.5) 
     *    WRITE(LSCREN,'(''INFOER=5 ==>> Minimisation finished -'',
     *    '' number of calls to fcn with iflag = 1 has reached maxfev'',
     *    '' (MAXFEV='',I6,'')'',
     *    '' The gradient norm is GRAD_N= '',f10.4,
     *    '' TOLERG= '',es12.5)')MAXFEV,GRAD_N,TOLERG
C
          IF (INFOER.EQ.6) 
     *    WRITE(LSCREN,'(''INFOER=6 ==>> Minimisation finished - ftol'',
     *    '' is too small. no further reduction in the sum of'',
     *    '' squares is possible'','' (FTOL='',ES12.5,'')'',
     *    '' The gradient norm is GRAD_N= '',f10.4,
     *    '' TOLERG= '',es12.5)')TOLERF,GRAD_N,TOLERG
C
          IF (INFOER.EQ.7) 
     *    WRITE(LSCREN,'(''INFOER=7 ==>> Minimisation finished - xtol'',
     *    '' is too small. No further improvement in the approximate'',
     *    '' solution x is possible'','' (XTOL='',ES12.5,'')'',
     *    '' The gradient norm is GRAD_N= '',f10.4,
     *    '' TOLERG= '',es12.5)')TOLERX,GRAD_N,TOLERG
C
          IF (INFOER.EQ.8) 
     *    WRITE(LSCREN,'(''INFOER=8 ==>> Minimisation finished - gtol'',
     *    '' is too small. fvec is orthogonal to the columns of the'',
     *    '' jacobian to machine precision'',
     *    '' (GTOL='',ES12.5,'')'',
     *    '' The gradient norm is GRAD_N= '',f10.4,
     *    '' TOLERG= '',es12.5)')TOLERG,GRAD_N,TOLERG
C
          IF (INFOER.EQ.9) 
     *    WRITE(LSCREN,'(''INFOER=9 ==>> Minimisation finished - the'',
     *           '' difference between the last '',I2,'' chi^2 is'',
     *           '' less than EPSLAS= '',ES8.2,
     *           '' The gradient norm is GRAD_N= '',f10.4,
     *           '' TOLERG= '',es12.5)')
     *              LDLAST,EPSLAS,GRAD_N,TOLERG
C
          IF (INFOER.EQ.10) 
     *    WRITE(LSCREN,'(''INFOER=10 ==>> Minimisation finished - the'',
     *           '' \chi^2 gradient is less than TOLERG= '',ES12.5,
     *           '' The gradient norm is GRAD_N= '',f10.4)')
     *             TOLERG,GRAD_N
      
          WRITE(LSCREN,'()')
C      
      END IF
C
C=======================================================================
C
C     Printing the information in the logfile with the minimisation
C                                                             steps
CID      IF (LOGWRI.GT.0) THEN
C      
          WRITE(LOGAUX,'()')
C
          IF (INFOER.EQ.0) 
     *    WRITE(LOGAUX,'(''INFOER=0 ==>> Minimisation failed -'',
     *                   '' improper input parameters'')')
C
          IF (INFOER.EQ.1) 
     *    WRITE(LOGAUX,'(''INFOER=1 ==>> Minimisation finished - both'',
     *                 '' actual and predicted relative reductions'',
     *                 '' in the sum of squares are at most ftol'',
     *                 '' (FTOL='',E12.5,'')'',
     *                 '' The gradient norm is GRAD_N= '',f10.4,
     *                 '' TOLERG= '',es12.5)') TOLERF,GRAD_N,TOLERG
C
          IF (INFOER.EQ.2) 
     *    WRITE(LOGAUX,'(''INFOER=2 ==>> Minimisation finished -'',
     *    '' relative error between two consecutive iterates is at'',
     *    '' most xtol (XTOL='',ES12.5,'')'',
     *    '' The gradient norm is GRAD_N= '',f10.4,
     *    '' TOLERG= '',es12.5)') TOLERX,GRAD_N,TOLERG
C
          IF (INFOER.EQ.3) 
     *    WRITE(LOGAUX,'(''INFOER=3 ==>> Minimisation finished - both'',
     *    '' actual and predicted relative reductions'',
     *    '' in the sum of squares are at most ftol AND relative'',
     *    '' error between two consecutive iterates is at most xtol'',
     *    '' (FTOL='',ES12.5,'')'','' (XTOL='',ES12.5,'')'',
     *    '' The gradient norm is GRAD_N= '',f20.13,
     *    '' TOLERG= '',es12.5)')TOLERF,TOLERX,GRAD_N,TOLERG
C
          IF (INFOER.EQ.4) 
     *    WRITE(LOGAUX,'(''INFOER=4 ==>> Minimisation finished - the'',
     *    '' cosine of the angle between fvec and any column of the'',
     *    '' jacobian is at most gtol in absolute value'',
     *    '' (GTOL='',ES12.5,'')'',
     *    '' The gradient norm is GRAD_N= '',f20.13,
     *    '' TOLERG= '',es12.5)')TOLERG,GRAD_N,TOLERG
C
          IF (INFOER.EQ.5) 
     *    WRITE(LOGAUX,'(''INFOER=5 ==>> Minimisation finished -'',
     *    '' number of calls to fcn with iflag = 1 has reached maxfev'',
     *    '' (MAXFEV='',I6,'')'',
     *    '' The gradient norm is GRAD_N= '',f20.13,
     *    '' TOLERG= '',es12.5)')MAXFEV,GRAD_N,TOLERG
C
          IF (INFOER.EQ.6) 
     *    WRITE(LOGAUX,'(''INFOER=6 ==>> Minimisation finished - ftol'',
     *    '' is too small. no further reduction in the sum of'',
     *    '' squares is possible'','' (FTOL='',ES12.5,'')'',
     *    '' The gradient norm is GRAD_N= '',f20.13,
     *    '' TOLERG= '',es12.5)')TOLERF,GRAD_N,TOLERG
C
          IF (INFOER.EQ.7) 
     *    WRITE(LOGAUX,'(''INFOER=7 ==>> Minimisation finished - xtol'',
     *    '' is too small. No further improvement in the approximate'',
     *    '' solution x is possible'','' (XTOL='',ES12.5,'')'',
     *    '' The gradient norm is GRAD_N= '',f20.13,
     *    '' TOLERG= '',es12.5)')TOLERX,GRAD_N,TOLERG
C
          IF (INFOER.EQ.8) 
     *    WRITE(LOGAUX,'(''INFOER=8 ==>> Minimisation finished - gtol'',
     *    '' is too small. fvec is orthogonal to the columns of the'',
     *    '' jacobian to machine precision'',
     *    '' (GTOL='',ES12.5,'')'',
     *    '' The gradient norm is GRAD_N= '',f20.13,
     *    '' TOLERG= '',es12.5)')TOLERG,GRAD_N,TOLERG
C
          IF (INFOER.EQ.9) 
     *    WRITE(LOGAUX,'(''INFOER=9 ==>> Minimisation finished - the'',
     *           '' difference between the last '',I2,'' chi^2 is'',
     *           '' less than EPSLAS= '',ES8.2,
     *           '' The gradient norm is GRAD_N= '',f20.13,
     *           '' TOLERG= '',es12.5)')
     *              LDLAST,EPSLAS,GRAD_N,TOLERG
C
          IF (INFOER.EQ.10) 
     *    WRITE(LOGAUX,'(''INFOER=10 ==>> Minimisation finished - the'',
     *           '' \chi^2 gradient is less than TOLERG= '',ES12.5,
     *           '' The gradient norm is GRAD_N= '',f20.13)')
     *             TOLERG,GRAD_N
      
          WRITE(LOGAUX,'()')
C      
CID      END IF
C
C=======================================================================
C
C     Printing the information in the main logfile
C
      IF (LOGWRI.GT.0) THEN
C      
          WRITE(LOGFIL,'()')
C
          IF (INFOER.EQ.0) 
     *    WRITE(LOGFIL,'(''INFOER=0 ==>> Minimisation failed -'',
     *                   '' improper input parameters'')')
C
          IF (INFOER.EQ.1) 
     *    WRITE(LOGFIL,'(''INFOER=1 ==>> Minimisation finished - both'',
     *                 '' actual and predicted relative reductions'',
     *                 '' in the sum of squares are at most ftol'',
     *                 '' (FTOL='',E12.5,'')'',
     *                 '' The gradient norm is GRAD_N= '',f20.13,
     *                 '' TOLERG= '',es12.5)') TOLERF,GRAD_N,TOLERG
C
          IF (INFOER.EQ.2) 
     *    WRITE(LOGFIL,'(''INFOER=2 ==>> Minimisation finished -'',
     *    '' relative error between two consecutive iterates is at'',
     *    '' most xtol (XTOL='',ES12.5,'')'',
     *    '' The gradient norm is GRAD_N= '',f20.13,
     *    '' TOLERG= '',es12.5)') TOLERX,GRAD_N,TOLERG
C
          IF (INFOER.EQ.3) 
     *    WRITE(LOGFIL,'(''INFOER=3 ==>> Minimisation finished - both'',
     *    '' actual and predicted relative reductions'',
     *    '' in the sum of squares are at most ftol AND relative'',
     *    '' error between two consecutive iterates is at most xtol'',
     *    '' (FTOL='',ES12.5,'')'','' (XTOL='',ES12.5,'')'',
     *    '' The gradient norm is GRAD_N= '',f20.13,
     *    '' TOLERG= '',es12.5)')TOLERF,TOLERX,GRAD_N,TOLERG
C
          IF (INFOER.EQ.4) 
     *    WRITE(LOGFIL,'(''INFOER=4 ==>> Minimisation finished - the'',
     *    '' cosine of the angle between fvec and any column of the'',
     *    '' jacobian is at most gtol in absolute value'',
     *    '' (GTOL='',ES12.5,'')'',
     *    '' The gradient norm is GRAD_N= '',f20.13,
     *    '' TOLERG= '',es12.5)')TOLERG,GRAD_N,TOLERG
C
          IF (INFOER.EQ.5) 
     *    WRITE(LOGFIL,'(''INFOER=5 ==>> Minimisation finished -'',
     *    '' number of calls to fcn with iflag = 1 has reached maxfev'',
     *    '' (MAXFEV='',I6,'')'',
     *    '' The gradient norm is GRAD_N= '',f20.13,
     *    '' TOLERG= '',es12.5)')MAXFEV,GRAD_N,TOLERG
C
          IF (INFOER.EQ.6) 
     *    WRITE(LOGFIL,'(''INFOER=6 ==>> Minimisation finished - ftol'',
     *    '' is too small. no further reduction in the sum of'',
     *    '' squares is possible'','' (FTOL='',ES12.5,'')'',
     *    '' The gradient norm is GRAD_N= '',f20.13,
     *    '' TOLERG= '',es12.5)')TOLERF,GRAD_N,TOLERG
C
          IF (INFOER.EQ.7) 
     *    WRITE(LOGFIL,'(''INFOER=7 ==>> Minimisation finished - xtol'',
     *    '' is too small. No further improvement in the approximate'',
     *    '' solution x is possible'','' (XTOL='',ES12.5,'')'',
     *    '' The gradient norm is GRAD_N= '',f20.13,
     *    '' TOLERG= '',es12.5)')TOLERX,GRAD_N,TOLERG
C
          IF (INFOER.EQ.8) 
     *    WRITE(LOGFIL,'(''INFOER=8 ==>> Minimisation finished - gtol'',
     *    '' is too small. fvec is orthogonal to the columns of the'',
     *    '' jacobian to machine precision'',
     *    '' (GTOL='',ES12.5,'')'',
     *    '' The gradient norm is GRAD_N= '',f20.13,
     *    '' TOLERG= '',es12.5)')TOLERG,GRAD_N,TOLERG
C
          IF (INFOER.EQ.9) 
     *    WRITE(LOGFIL,'(''INFOER=9 ==>> Minimisation finished - the'',
     *           '' difference between the last '',I2,'' chi^2 is'',
     *           '' less than EPSLAS= '',ES8.2,
     *           '' The gradient norm is GRAD_N= '',f20.13,
     *           '' TOLERG= '',es12.5)')
     *              LDLAST,EPSLAS,GRAD_N,TOLERG
C
          IF (INFOER.EQ.10) 
     *    WRITE(LOGFIL,'(''INFOER=10 ==>> Minimisation finished - the'',
     *           '' \chi^2 gradient is less than TOLERG= '',ES12.5,
     *           '' The gradient norm is GRAD_N= '',f20.13)')
     *             TOLERG,GRAD_N
      
          WRITE(LOGFIL,'()')
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
      FUNCTION EUNORM(NDPARS,LDPARS,XARGUM)
C
      DIMENSION 
     *          XARGUM(1:NDPARS)
C
      DATA 
     *       R_UNIT,R_ZERO, RDWARF  , RGIANT 
     *     / 1.0000,0.0000,3.834E-20,1.304E19 /
C
C=======================================================================
C
C     Given an (n=LDPARS)-vector XARGUM -> this function calculates 
C     the Euclidean norm of XARGUM.
C
C     We introduce the definitions of small, intermediate and large 
C     components depend on two constants, RDWARF and RGIANT.  There
C     are restrictions on these constants, such that RDWARF**2  not
C     underflow and RGIANT**2 not overflow. The values of constants 
C     given here are suitable for every known computer.
C
C     The Euclidean norm is then computed  by collecting the sum of
C     squares in three different sums.  The sums of squares for the
C     small and large components  are scaled, so that no  overflows
C     occur.  Non-destructive underflows are permitted.  Underflows
C     and overflows do not occur in the computation of the unscaled
C     sum of squares for the intermediate components.
C
C     The function statement is
C
C                  FUNCTION EUNORM(NDPARS,LDPARS,XARGUM)
C
C     where
C
C     LDPARS - a positive integer input variable, the actual number
C              of parameters (dimension of the problem)
C
C     XARGUM - an input array of length  LDPARS  that is not larger 
C              than NDPARS
C
C=======================================================================
C
C     SUBPROGRAMS CALLED
C
C                 FORTRAN-SUPPLIED ... ABS, SQRT
C
C     Argonne National Laboratory. MINPACK PROJECT, March 1980
C     Burton S. Garbow, Kenneth E. Hillstrom and Jorge J. More
C
C=======================================================================
C 
      CALL CPUTIM('EUNORM',1)
C
C=======================================================================
C
      S1TERM=R_ZERO
      S2TERM=R_ZERO
      S3TERM=R_ZERO
C
      X1_MAX=R_ZERO
      X3_MAX=R_ZERO
C
      FLOATN=LDPARS
      AGIANT=RGIANT/FLOATN
C
      DO I=1,LDPARS
C
         XARABS=ABS(XARGUM(I))
C
         IF (.NOT.(XARABS.GT.RDWARF.AND.XARABS.LT.AGIANT)) THEN
C
             IF (XARABS.GT.RDWARF) THEN
C
C                Sum for large components
C
                 IF (XARABS.GT.X1_MAX) THEN
C
                     S1TERM=R_UNIT+S1TERM*(X1_MAX/XARABS)**2
                     X1_MAX=XARABS
C
                     GO TO 20
C
                 END IF
C
                 S1TERM=S1TERM+(XARABS/X1_MAX)**2
C
   20            CONTINUE
                 GO TO 60
             END IF
C
C            Sum for small components
C
             IF (XARABS.GT.X3_MAX) THEN
                 S3TERM=R_UNIT+S3TERM*(X3_MAX/XARABS)**2
                 X3_MAX=XARABS
                 GO TO 50
             END IF
C
             IF (XARABS.NE.R_ZERO) S3TERM=S3TERM+(XARABS/X3_MAX)**2
C
   50        CONTINUE
   60        CONTINUE
C
             GO TO 80
C
         END IF
C
C        Sum for intermediate components
C
         S2TERM=S2TERM+XARABS**2
C
   80    CONTINUE
C
      END DO
C
C     Calculation of norm
C
      IF (S1TERM.NE.R_ZERO) THEN
          EUNORM=X1_MAX*SQRT(S1TERM+(S2TERM/X1_MAX)/X1_MAX)
          GO TO 130
      END IF
C
      IF (S2TERM.NE.R_ZERO) THEN
C
          IF (S2TERM.GE.X3_MAX) THEN
              EUNORM=SQRT(S2TERM*(R_UNIT+(X3_MAX/S2TERM)
     *                          *(X3_MAX*S3TERM)))
          END IF
C
          IF (S2TERM.LT.X3_MAX) THEN
              EUNORM=SQRT(X3_MAX*((S2TERM/X3_MAX)+(X3_MAX*S3TERM)))
          END IF
C
          GO TO 120
C
      END IF
C
      EUNORM=X3_MAX*SQRT(S3TERM)
C
  120 CONTINUE
  130 CONTINUE
C
C=======================================================================
C 
      CALL CPUTIM('EUNORM',0)
C
C=======================================================================
C
      RETURN
      END
C
C=======================================================================
C=======================================================================
C         
      SUBROUTINE LM_PAR(NDFUNC,LDFUNC,NDPARS,LDPARS,RUPTRI,I_PERM,
     *                  DIAGSC,QTRANB,DELTAB,PARLEV,XARGUM,S_DIAG,
     *                                              WORKA1,WORKA2)
C 
      DIMENSION 
     *          DIAGSC(1:NDPARS),QTRANB(1:NDPARS),
     *          XARGUM(1:NDPARS),S_DIAG(1:NDPARS),
     *          WORKA1(1:NDPARS),WORKA2(1:NDPARS)
      DIMENSION 
     *          RUPTRI(1:NDFUNC,1:NDPARS)
      DIMENSION 
     *          I_PERM(1:NDPARS)
C
      DATA 
     *      CONST1,P001,R_ZERO 
     *     /1.0D-1,1.0D-3,0.0D0/
C
C=======================================================================
C
C     Given an m by NDPARS matrix A, an NDPARS by NDPARS nonsingular 
C     diagonal matrix D, an m-vector B, and a positive number DELTAB,
C     The problem is to determine  a value for the parameter  PARLEV 
C     such that if XARGUM solves the system
C
C               D*XARGUM = D   and   sqrt(PARLEV)*D*XARGUM = 0,
C
C     in the least squares sense and DXNORM is the Euclidean norm of 
C     D*XARGUM, then either PARLEV is zero and
C
C                     (DXNORM-DELTAB).LE.0.1*DELTAB,
C
C     or PARLEV is positive and
C
C                   abs(DXNORM-DELTAB).LE.0.1*DELTAB 
C
C     This subroutine completes the solution of the problem if it is 
C     provided with the necessary information from  QR factorization, 
C     with column pivoting, of A,  that is, if A*P = Q*RUPTRI, where 
C     P is a permutation matrix, Q has orthogonal columns and RUPTRI 
C     is an upper triangular matrix  with diagonal elements that are 
C     of nonincreasing magnitude, then LM_PAR expects the full upper 
C     triangle of  RUPTRI,  the permutation matrix P,  and the first 
C     LDPARS components of (Q Transpose) * B. On output, LM_PAR also 
C     provides an upper triangular matrix s such that
C
C                     t   t                   t
C                    P *(A *A + PARLEV*D*D)*P = S *S .
C
C     S is employed within LM_PAR ->> it may be of separate interest
C
C     Only a few iterations  are generally needed for convergence of 
C     the algorithm.  If the limit of 10 iterations is reached, then 
C     the output  PARLEV will contain the best value obtained so far.
C
C     The subroutine statement is
C
C     SUBROUTINE LM_PAR(NDFUNC,LDFUNC,NDPARS,LDPARS,RUPTRI,I_PERM,
C                       DIAGSC,QTRANB,DELTAB,PARLEV,XARGUM,S_DIAG,
C                                                   WORKA1,WORKA2)
C
C     where
C
C     LDPARS - a positive integer input variable set to the order of 
C                                                             RUPTRI
C
C     RUPTRI - a square  NDPARS by NDPARS  array.  On input the full 
C              upper triangle  must contain  the full upper triangle 
C              of the matrix RUPTRI;  then on output, the full upper 
C              triangle is unaltered,  and the strict lower triangle 
C              contains the strict  upper triangle  (transposed)  of 
C              the upper triangular matrix S.
C
C     NDFUNC - a positive,  integer input variable  not smaller than 
C              NDPARS, which specifies  the leading dimension of the 
C              array RUPTRI.
C
C     I_PERM - an integer input array of length NDPARS which defines 
C              the permutation matrix p such that A * P = Q * RUPTRI. 
C              Column J of P is column I_PERM(j) of identity matrix.
C
C     DIAGSC - an input array of length NDPARS which must contain the
C              diagonal elements of the matrix D
C
C     QTRANB - an input array of length NDPARS which must contain the 
C              first LDPARS elements of the vector (Q Transpose) * B.
C
C     DELTAB - a positive input variable, it specifies an upper bound 
C              on the Euclidean norm of D*XARGUM
C
C     PARLEV - a nonnegative variable. containing on input an initial 
C              estimate of the Levenberg-Marquardt parameter.  On the 
C              output PARLEV contains the final estimate.
C
C     XARGUM - an output array of length  LDPARS,  which contains the 
C              least squares solution of the system 
C
C                  A*XARGUM = B   and   sqrt(PARLEV)*d*XARGUM = 0,
C
C              for the output PARLEV.
C
C     S_DIAG - an output array of length LDPARS.  It contains all the
C              diagonal elements of the upper triangular matrix S.
C
C              WORKA1 and WORKA2 are work vectors of length NDPARS.
C
C=======================================================================
C
C     SUBPROGRAMS CALLED
C
C                 MINPACK-SUPPLIED ... DPMPAR, EUNORM, QRSOLV
C
C                 FORTRAN-SUPPLIED ... ABS, MAX, MIN, SQRT
C
C     Argonne National Laboratory. MINPACK PROJECT, March 1980
C     Burton S. Garbow, Kenneth E. Hillstrom and Jorge J. More
C
C=======================================================================
C 
      CALL CPUTIM('LM_PAR',1)
C
C=======================================================================
C
C     DWARF is the smallest positive magnitude.
C
      DWARF=DPMPAR(2)
C
C     Compute and store in XARGUM the Gauss-Newton direction. If the
C     Jacobian is rank-deficient, obtain a least squares solution.
C
      N_SING=LDPARS
C
      DO J=1,LDPARS
         WORKA1(J)=QTRANB(J)
         IF (RUPTRI(J,J).EQ.R_ZERO.AND.N_SING.EQ.LDPARS) N_SING=J-1
         IF (N_SING.LT.LDPARS) WORKA1(J)=R_ZERO
      END DO
C
      IF (N_SING.GE.1) THEN
C
          DO K=1,N_SING
C
             J=N_SING-K+1
             WORKA1(J)=WORKA1(J)/RUPTRI(J,J)
             TEMPOR=WORKA1(J)
             JMIN_1=J-1
C
             IF (JMIN_1.GE.1) THEN
                 DO I=1,JMIN_1
                    WORKA1(I)=WORKA1(I)-RUPTRI(I,J)*TEMPOR
                 END DO
             END IF
C
          END DO
C
      END IF
C
      DO J=1,LDPARS
         L=I_PERM(J)
         XARGUM(L)=WORKA1(J)
      END DO
C
C=======================================================================
C
C     Initialize the iteration counter. Evaluate the function at
C     the origin,  and test  for acceptance  of the gauss-newton 
C     direction.
C
      ITERAT=0
C
      DO J=1,LDPARS
         WORKA2(J)=DIAGSC(J)*XARGUM(J)
      END DO
C
      DXNORM=EUNORM(NDPARS,LDPARS,WORKA2)
      FP_AUX=DXNORM-DELTAB
C
      IF (FP_AUX.LE.CONST1*DELTAB) GO TO 220
C
C=======================================================================
C
C     If the Jacobian is not rank deficient, the Newton step provides 
C     a lower bound, PARLOW, for the zero of the function.  Otherwise 
C     set this bound to zero.
C
      PARLOW=R_ZERO
C
      IF (N_SING.GE.LDPARS) THEN
C
          DO J=1,LDPARS
             L=I_PERM(J)
             WORKA1(J)=DIAGSC(L)*(WORKA2(L)/DXNORM)
          END DO
C
          DO J=1,LDPARS
C
             SUMAUX=R_ZERO
             JMIN_1=J-1
C
             IF (JMIN_1.GE.1) THEN
                 DO I=1,JMIN_1
                    SUMAUX=SUMAUX+RUPTRI(I,J)*WORKA1(I)
                 END DO
             END IF
C
             WORKA1(J)=(WORKA1(J)-SUMAUX)/RUPTRI(J,J)
C
          END DO
C
          TEMPOR=EUNORM(NDPARS,LDPARS,WORKA1)
          PARLOW=((FP_AUX/DELTAB)/TEMPOR)/TEMPOR
C
      END IF
C
C=======================================================================
C
C     Calculate an upper bound, PARUPP, for the zero of the function
C
      DO J=1,LDPARS
         SUMAUX=R_ZERO
         DO I=1,J
            SUMAUX=SUMAUX+RUPTRI(I,J)*QTRANB(I)
         END DO
         L=I_PERM(J)
         WORKA1(J)=SUMAUX/DIAGSC(L)
      END DO
C
      G_NORM=EUNORM(NDPARS,LDPARS,WORKA1)
      PARUPP=G_NORM/DELTAB
C
      IF (PARUPP.EQ.R_ZERO) PARUPP=DWARF/MIN(DELTAB,CONST1)
C
C=======================================================================
C
C     If the input PARLEV lies outside of the interval (PARLOW,PARUPP)
C     set PARLEV to the closer endpoint
C
      PARLEV=MAX(PARLEV,PARLOW)
      PARLEV=MIN(PARLEV,PARUPP)
C
      IF (PARLEV.EQ.R_ZERO) PARLEV=G_NORM/DXNORM
C
C     Beginning of an iteration
C
  150 CONTINUE
C
      ITERAT=ITERAT+1
C
C=======================================================================
C
C     Evaluate the function at the current value of PARLEV
C
      IF (PARLEV.EQ.R_ZERO) PARLEV=MAX(DWARF,P001*PARUPP)
C
      TEMPOR=SQRT(PARLEV)
C
      DO J=1,LDPARS
         WORKA1(J)=TEMPOR*DIAGSC(J)
      END DO
C  
      CALL QRSOLV(NDFUNC,LDFUNC,NDPARS,LDPARS,RUPTRI,I_PERM,
     *                   WORKA1,QTRANB,XARGUM,S_DIAG,WORKA2)
C  
      DO J=1,LDPARS
         WORKA2(J)=DIAGSC(J)*XARGUM(J)
      END DO
C  
      DXNORM=EUNORM(NDPARS,LDPARS,WORKA2)
      TEMPOR=FP_AUX
      FP_AUX=DXNORM-DELTAB
C
C=======================================================================
C
C     If the function is small enough,  accept the current  value
C     of PARLEV. also test for the exceptional cases where PARLOW
C     is zero or the number of iterations has reached 10
C
      IF (ABS(FP_AUX).LE.CONST1*DELTAB.OR.PARLOW.EQ.R_ZERO.AND.
     *    FP_AUX.LE.TEMPOR.AND.TEMPOR.LT.R_ZERO.OR.ITERAT.EQ.10) 
     *                                                GO TO 220
C
C=======================================================================
C
C     Compute the Newton correction
C
      DO J=1,LDPARS
         L=I_PERM(J)
         WORKA1(J)=DIAGSC(L)*(WORKA2(L)/DXNORM)
      END DO
C
      DO J=1,LDPARS
C
         WORKA1(J)=WORKA1(J)/S_DIAG(J)
         TEMPOR=WORKA1(J)
         JP1=J+1
C
         IF (LDPARS.GE.JP1) THEN
C
             DO I=JP1,LDPARS
                WORKA1(I)=WORKA1(I)-RUPTRI(I,J)*TEMPOR
             END DO
C
          END IF
C
      END DO
C
      TEMPOR=EUNORM(NDPARS,LDPARS,WORKA1)
      PARC=((FP_AUX/DELTAB)/TEMPOR)/TEMPOR
C
C=======================================================================
C
C     Depending on the sign of the function, update PARLOW or PARUPP
C
      IF (FP_AUX.GT.R_ZERO) PARLOW=MAX(PARLOW,PARLEV)
      IF (FP_AUX.LT.R_ZERO) PARUPP=MIN(PARUPP,PARLEV)
C
C=======================================================================
C
C     Compute an improved estimate for PARLEV.
C
      PARLEV=MAX(PARLOW,PARLEV+PARC)
C
C     End of an iteration.
C
      GO TO 150
C
  220 CONTINUE
C
      IF (ITERAT.EQ.0) PARLEV=R_ZERO
C
C=======================================================================
C 
      CALL CPUTIM('LM_PAR',0)
C
C=======================================================================      
C      
      RETURN
      END 
C
C=======================================================================      
C=======================================================================      
C
      SUBROUTINE QRFACT(NDFUNC,LDFUNC,NDPARS,LDPARS,AMATRX,LPIVOT,
     *                                I_PERM,R_DIAG,ACNORM,WORK_A)
C
      LOGICAL 
     *          LPIVOT
      DIMENSION 
     *          R_DIAG(1:NDPARS),ACNORM(1:NDPARS),WORK_A(1:NDPARS)
      DIMENSION 
     *          AMATRX(1:NDFUNC,1:NDPARS)
      DIMENSION 
     *          I_PERM(1:NDPARS)
C
      DATA 
     *     R_UNIT,P05,R_ZERO 
     *    /1.0000,5.0e-2,0.000/
C
C=======================================================================
C
C     This subroutine uses  Householder transformations  with column
C     pivoting to compute a QR factorisation of the LDFUNC by LDPARS 
C     matrix AMATRX. that is, QRFACT determines an orthogonal matrix 
C     Q, a permutation matrix P,  together with an upper trapezoidal 
C     matrix RUPTRI, with diagonal elements  of nonincreasing magni-
C     tude, such that A*P = Q*RUPTRI. The Householder transformation 
C     for column K, K = 1,2, ... ,MIN(LDFUNC,LDPARS), is of the form
C
C                                          t
C                          I - (1/U(K))*U*U
C
C     where U has zeros in the first K-1 positions. the form of this 
C     transformation, and the method of pivoting  has first appeared 
C     in the corresponding LINPACK subroutine.
C
C     The subroutine statement is
C
C       SUBROUTINE QRFACT(NDFUNC,LDFUNC,NDPARS,LDPARS,AMATRX,LPIVOT,
C                                       I_PERM,R_DIAG,ACNORM,WORK_A)
C
C     where
C
C     LDFUNC - a positive integer input variable,  set to the number
C              of rows of AMATRX
C
C     LDPARS - a positive integer input variable,  set to the number
C              of columns of AMATRX
C
C     AMATRX - an LDFUNC by LDPARS array.  On input  AMATRX contains 
C              the matrix,  for which the QR factorization  is to be 
C              computed. On output the strict upper trapezoidal part 
C              of AMATRX  contains the strict upper trapezoidal part 
C              of RUPTRI, and the lower  trapezoidal part  of AMATRX 
C              contains  a factored form of Q  (non-trivial elements 
C              of the U vectors described above).
C
C     LDFUNC - a positive integer input variable;  it  specifies the 
C              leading dimension of the array AMATRX not larger than 
C                                                             NDFUNC
C
C     LPIVOT - a logical input variable. If LPIVOT is set true, then 
C              column pivoting is enforced.  If LPIVOT is set false,
C              then no column pivoting is done.
C
C     I_PERM - an integer output array of length NDPARS;  it defines 
C              the permutation matrix P such that AMATRX*P=Q*RUPTRI.
C              Column  J of P  is column  I_PERM(J)  of the identity 
C              matrix. If LPIVOT is false, I_PERM is not referenced.
C
C     LPIVOT - a positive integer input variable. If LPIVOT is false
C              then LPIVOT may be as small as 1.  If LPIVOT is true, 
C              then LPIVOT must be at least LDPARS.
C
C     R_DIAG - an output array of length LDPARS,  which contains the
C              diagonal elements of RUPTRI.
C
C     ACNORM - an output array of length LDPARS,  which contains the
C              norms  of  the corresponding  columns  of  the  input 
C              matrix AMATRX. If this information is not needed then 
C              ACNORM can coincide with R_DIAG.
C
C     WORK_A - a work  array  of length  LDPARS. If LPIVOT is false, 
C              then WORK_A can coincide with R_DIAG.
C
C     SUBPROGRAMS CALLED
C
C                 MINPACK-SUPPLIED ... DPMPAR, EUNORM
C
C                 FORTRAN-SUPPLIED ... MAX, SQRT, MIN0
C
C=======================================================================
C
C     Argonne National Laboratory. MINPACK PROJECT, March 1980
C     Burton S. Garbow, Kenneth E. Hillstrom and Jorge J. More
C
C=======================================================================
C
      CALL CPUTIM('QRFACT',1)
C
C=======================================================================
C
      EPSMCH=DPMPAR(1) !  EPSMCH is the machine precision
C
C     Compute the initial column norms and initialize several arrays
C
      DO J=1,LDPARS
C
         ACNORM(J)=EUNORM(NDFUNC,LDFUNC,AMATRX(1,J))
         R_DIAG(J)=ACNORM(J)
         WORK_A(J)=R_DIAG(J)
C
         IF (LPIVOT) I_PERM(J)=J
C
      END DO
C
C     Reduce AMATRX to RUPTRI with Householder transformations
C
      MINMN=MIN0(LDFUNC,LDPARS)
C
      DO J=1,MINMN
C
         IF (LPIVOT) THEN
C
C            Bring the column of largest norm into the LPIVOT position
C
             KMAX=J
C
             DO K=J,LDPARS
                IF (R_DIAG(K).GT.R_DIAG(KMAX)) KMAX=K
             END DO
C
             IF (KMAX.NE.J) THEN
C
                 DO I=1,LDFUNC
                    TEMPOR=AMATRX(I,J)
                    AMATRX(I,J)=AMATRX(I,KMAX)
                    AMATRX(I,KMAX)=TEMPOR
                 END DO
C
                 R_DIAG(KMAX)=R_DIAG(J)
                 WORK_A(KMAX)=WORK_A(J)
                 K=I_PERM(J)
                 I_PERM(J)=I_PERM(KMAX)
                 I_PERM(KMAX)=K
C
             END IF
C
         END IF
C
C        Compute the Householder transformation to reduce the
C        J-th column of AMATRX to a multiple of the J-th unit 
C                                                      vector
C
         AJNORM=EUNORM(NDFUNC,LDFUNC-J+1,AMATRX(J,J))
C
         IF (AJNORM.NE.R_ZERO) THEN
C
             IF (AMATRX(J,J).LT.R_ZERO) AJNORM=-AJNORM
C
             DO I=J,LDFUNC
                AMATRX(I,J)=AMATRX(I,J)/AJNORM
             END DO
C
             AMATRX(J,J)=AMATRX(J,J)+R_UNIT
C
C            Apply the transformation to the remaining columns
C                                         and update the norms
             JP1=J+1
C
             IF (LDPARS.GE.JP1) THEN
C
                 DO K=JP1,LDPARS
C
                    SUMAUX=R_ZERO
C
                    DO I=J,LDFUNC
                       SUMAUX=SUMAUX+AMATRX(I,J)*AMATRX(I,K)
                    END DO
C
                    TEMPOR=SUMAUX/AMATRX(J,J)
C
                    DO I=J,LDFUNC
                       AMATRX(I,K)=AMATRX(I,K)-TEMPOR*AMATRX(I,J)
                    END DO
C
                    IF (.NOT.LPIVOT.OR.R_DIAG(K).EQ.R_ZERO) GO TO 80
C
                    TEMPOR=AMATRX(J,K)/R_DIAG(K)
                    R_DIAG(K)
     *             =
     *              R_DIAG(K)*SQRT(MAX(R_ZERO,R_UNIT-TEMPOR**2))
C
                    IF (P05*(R_DIAG(K)/WORK_A(K))**2.GT.EPSMCH) GO TO 80
C
                    R_DIAG(K)=EUNORM(NDFUNC,LDFUNC-J,AMATRX(JP1,K))
                    WORK_A(K)=R_DIAG(K)
C
   80               CONTINUE
C
                 END DO
C
             END IF
C
         END IF
C
         R_DIAG(J)=-AJNORM
C
      END DO
C
C=======================================================================
C
      CALL CPUTIM('QRFACT',0)
      
      RETURN
      END 
C
C=======================================================================
C=======================================================================
C  
      SUBROUTINE QRSOLV(NDFUNC,LDFUNC,NDPARS,LDPARS,RUPTRI,I_PERM,
     *                         DIAGSC,QTRANB,XARGUM,S_DIAG,WAUXIL)
C
      DIMENSION 
     *          DIAGSC(1:NDPARS),QTRANB(1:NDPARS),
     *          XARGUM(1:NDPARS),S_DIAG(1:NDPARS),
     *                           WAUXIL(1:NDPARS)
      DIMENSION 
     *          RUPTRI(1:NDFUNC,1:NDPARS)
      DIMENSION 
     *          I_PERM(1:NDPARS)
C
      DATA 
     *          CONST5,CONS25,R_ZERO 
     *        / 5.0e-1,2.5e-1,0.0000/
C
C=======================================================================
C
C     Given an LDFUNC by LDPARS matrix, say A, an LDPARS by LDPARS 
C     diagonal matrix D,  and an (m=LDFUNC)-vector B;  the problem 
C     is to determine an XARGUM which solves the system
C
C                      A*XARGUM=B,   D*XARGUM=0,
C
C     in the least squares sense.
C
C     This subroutine completes the solution of the problem, if it 
C     is provided with the necessary information from the QR facto-
C     rization, with column pivoting, of A, that is, if
C
C                          A*P = Q*RUPTRI, 
C
C     where P  is a permutation matrix,  Q  has orthogonal columns, 
C     and RUPTRI  is an upper-triangular  matrix  with all diagonal
C     elements of nonincreasing magnitude,  then QRSOLV expects the 
C     full upper triangle of RUPTRI, the permutation matrix P, and 
C     the first LDPARS components of (Q Transpose)*B. The system
C
C                      A*XARGUM = B, D*XARGUM = 0, 
C
C     is then equivalent to
C
C                               t       t
C                 RUPTRI * Z = Q *B ,  P * D * P * Z = 0 ,
C
C     where XARGUM = P*Z.  If this system does  not have full rank,
C     then a least squares solution  is obtained.  On output QRSOLV
C     also provides an upper triangular matrix S such that
C
C                          t   t               t
C                         P *(A *A + D*D)*P = S *S .
C 
C     S is computed within QRSOLV (and may be of separate interest)
C
C     The subroutine statement is
C
C     SUBROUTINE QRSOLV(NDFUNC,LDFUNC,NDPARS,LDPARS,RUPTRI,I_PERM,
C                              DIAGSC,QTRANB,XARGUM,S_DIAG,WAUXIL)
C
C     where
C
C     LDPARS - a positive, integer input variable  set to the order 
C                                                         of RUPTRI
C
C     RUPTRI - an LDPARS by LDPARS array.  On input, the full upper
C              triangle must contain the full upper triangle of the 
C              matrix RUPTRI.  On output the full upper triangle is 
C              unaltered and the strict lower triangle contains the 
C              strict upper triangle transposed of the upper trian-
C              gular matrix S.
C
C     NDFUNC - defines the number of functions not less than NDPARS;
C              at the same time, it specifies the leading dimension 
C              of the array RUPTRI
C
C     I_PERM - an integer input array of length LDPARS;  it defines 
C              the permutation matrix P such that 
C
C                               A*P = Q*RUPTRI
C
C              Column  J of P  is column I_PERM(j)  of the identity 
C                                                            matrix
C
C     DIAGSC - an input array of length LDPARS,  which must contain 
C              the diagonal elements of the matrix D
C
C     QTRANB - an input array of length LDPARS,  which must contain 
C              the first LDPARS elements of vector (Q Transpose)*B.
C
C     XARGUM - an output array of length LDPARS, which contains the 
C              least squares solution of the system 
C
C                          A*XARGUM = B, D*XARGUM = 0
C
C     S_DIAG - an output array of length LDPARS, which contains the
C              diagonal elements of the upper triangular matrix S
C
C     WAUXIL - a work array of length NDPARS.
C
C     SUBPROGRAMS CALLED
C
C                 FORTRAN-SUPPLIED ... ABS,SQRT
C
C=======================================================================
C
C     Argonne National Laboratory. MINPACK PROJECT, March 1980
C     Burton S. Garbow, Kenneth E. Hillstrom and Jorge J. More
C
C=======================================================================
C
C     Copy RUPTRI and (Q Transpose)*B to preserve input and initialize 
C     S. In particular, save the diagonal elements of RUPTRI in XARGUM
C
      DO J=1,LDPARS
         DO I=J,LDPARS
            RUPTRI(I,J)=RUPTRI(J,I)
         END DO
         XARGUM(J)=RUPTRI(J,J)
         WAUXIL(J)=QTRANB(J)
      END DO
C
C     Eliminate the diagonal matrix d using a givens rotation.
C
      DO J=1,LDPARS
C
C        Prepare the row of D to be eliminated, locating the
C        diagonal element using P from the QR factorization.
C
         L=I_PERM(J)
C
         IF (DIAGSC(L).EQ.R_ZERO) GO TO 90
C
         DO K=J,LDPARS
            S_DIAG(K)=R_ZERO
         END DO
C
         S_DIAG(J)=DIAGSC(L)
C
C        The transformations  to eliminate  the row  of D
C        modify only a single element of  (Q Transpose)*B
C        beyond the first LDPARS, which is initially zero.
C
         QTBPJ=R_ZERO
C
         DO K=J,LDPARS
C
C           Determine a Givens rotation which eliminates the
C           appropriate element in the current row of D
C
            IF (S_DIAG(K).NE.R_ZERO) THEN
C
                IF (ABS(RUPTRI(K,K)).LT.ABS(S_DIAG(K))) THEN
                    COTAN=RUPTRI(K,K)/S_DIAG(K)
                    SIN=CONST5/SQRT(CONS25+CONS25*COTAN**2)
                    COS=SIN*COTAN
                    GO TO 50
                END IF
C
                TAN=S_DIAG(K)/RUPTRI(K,K)
                COS=CONST5/SQRT(CONS25+CONS25*TAN**2)
                SIN=COS*TAN
C
   50           CONTINUE
C
C               Compute the modified diagonal element of RUPTRI 
C               and the modified element of ((Q Transpose) * B)
C
                RUPTRI(K,K)=COS*RUPTRI(K,K)+SIN*S_DIAG(K)
C
                TEMPOR=COS*WAUXIL(K)+SIN*QTBPJ
                QTBPJ=-SIN*WAUXIL(K)+COS*QTBPJ
                WAUXIL(K)=TEMPOR
C
C               Accumulate the tranformation in the row of s
C
                KP1=K+1
C
                IF (LDPARS.GE.KP1) THEN
                    DO I=KP1,LDPARS
                       TEMPOR=COS*RUPTRI(I,K)+SIN*S_DIAG(I)
                       S_DIAG(I)=-SIN*RUPTRI(I,K)+COS*S_DIAG(I)
                       RUPTRI(I,K)=TEMPOR
                    END DO
                END IF
            END IF
C
         END DO
   90    CONTINUE
C
C        Store the diagonal element of s and restore
C        the corresponding diagonal element of RUPTRI.
C
         S_DIAG(j)=RUPTRI(j,j)
         RUPTRI(j,j)=XARGUM(j)
C
      END DO
C
C     Solve the triangular system for Z. If the system issingular, 
C     then obtain a least squares solution.
C
      N_SING=LDPARS
C
      DO J=1,LDPARS
         IF (S_DIAG(J).EQ.R_ZERO.AND.N_SING.EQ.LDPARS) N_SING=J-1
         IF (N_SING.LT.LDPARS) WAUXIL(J)=R_ZERO
      END DO
C
      IF (N_SING.GE.1) THEN
C
          DO K=1,N_SING
C
             J=N_SING-K+1
             SUMAUX=R_ZERO
             JP1=J+1
C
             IF (N_SING.GE.JP1) THEN
                 DO I=JP1,N_SING
                    SUMAUX=SUMAUX+RUPTRI(I,J)*WAUXIL(I)
                 END DO
             END IF
C
             WAUXIL(J)=(WAUXIL(J)-SUMAUX)/S_DIAG(J)
C
          END DO
C
      END IF
C
C     Permute the components of z back to components of XARGUM
C
      DO J=1,LDPARS
         L=I_PERM(J)
         XARGUM(L)=WAUXIL(J)
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
      SUBROUTINE CHKDER(NDFUNC,LDFUNC,NDPARS,LDPARS,XARGUM,FARGUM,
     *                         FJACOB,XNEIGH,F_VECP,IFMODE,ERRGRA)
C
      DIMENSION XARGUM(1:NDPARS),FARGUM(1:NDFUNC),
     *          XNEIGH(1:NDPARS),F_VECP(1:NDFUNC),
     *                           ERRGRA(1:NDFUNC)
      DIMENSION 
     *          FJACOB(1:NDFUNC,1:NDPARS)
C
      DATA 
     *      FACTOR,R_UNIT,R_ZERO 
     *     /1.00e2,1.0000,0.0000/
C
C=======================================================================
C
C     This subroutine checks whether the gradients of LDFUNC nonlinear 
C     functions in LDPARS variables, evaluated at a point XARGUM,  are
C     consistent  with the functions  themselves.  The user  must call 
C     CHKDER twice, First with IFMODE = 1 and then with IFMODE = 2.
C
C     IFMODE = 1; On input XARGUM must contain the point of evaluation;
C                 on output, XNEIGH is set to a neighbouring point.
C
C     IFMODE = 2. on input, FARGUM  must contain the functions and the
C                 rows of FJACOB must contain the gradients of the res-
C                 pective functions, each evaluated at XARGUM;  F_VECP
C                 must contain the functions evaluated at XNEIGH. Next,
C                 on output,  ERRGRA contains measures  of correctness 
C                 of the respective gradients.
C
C     The subroutine does not perform reliably  if the cancellation or
C     rounding errors cause a severe loss of significance  in terms of
C     evaluation  of functionS.  Therefore,  none of the components of 
C     XARGUM  should be unusually small  (in particular, zero)  or any
C     other value which may cause loss of significance.
C
C     The subroutine statement is
C
C     SUBROUTINE CHKDER(NDFUNC,LDFUNC,NDPARS,LDPARS,XARGUM,FARGUM,
C                              FJACOB,XNEIGH,F_VECP,IFMODE,ERRGRA)
C
C     where
C
C     LDFUNC - a positive integer input variable, set to the number of 
C                                                            functions
C
C     LDPARS - a positive integer input variable, set to the number of
C                                                            variables
C
C     XARGUM - an input array of maximum length NDPARS
C
C     FARGUM - a vector  of length LDFUNC.  On input, when  IFMODE = 2,
C              FARGUM must contain the functions evaluated at XARGUM
C
C     FJACOB - an LDFUNC by LDPARS array.  On input when IFMODE=2, the
C              the rows of  FJACOB  must contain the gradients  of the 
C              respective functions evaluated at XARGUM.
C
C     NDFUNC - a positive integer input parameter not less than LDFUNC
C              which specifies the maximum  dimension of  array FJACOB
C
C     XNEIGH - a vector of length LDPARS.  On output, when IFMODE = 1,
C              XNEIGH is set to a neighboring point of XARGUM.
C
C     F_VECP - an array of length LDFUNC.  On input, when  IFMODE = 2,
C              F_VECP must contain the functions evaluated at XNEIGH.
C
C     IFMODE - an integer input variable,  set to 1  on the first call
C              and 2 on the second.  Other values of IFMODE  are equi-
C              valent to IFMODE = 1.
C
C     ERRGRA - a vector of length LDFUNC.  On output, when IFMODE = 2,
C              ERRGRA  contains measures of correctness of the respec-
C              tive gradients.  If there is  no severe  loss of signi-
C              ficance,  then if ERRGRA(i) is 1.0 the i-th gradient is 
C              correct, while if ERRGRA(i) is 0.0 the i-th gradient is 
C              incorrect. For values of ERRGRA between 0.0 and 1.0 the 
C              categorization is less certain.  In general, a value of 
C              ERRGRA(i)  greater  than 0.5  indicates  that  the i-th 
C              gradient is probably correct while a value of ERRGRA(i) 
C              less than 0.5 indicates that the i-th gradient is proba-
C              bly incorrect.
C
C     SUBPROGRAMS CALLED
C
C                 MINPACK SUPPLIED ... DPMPAR
C
C                 FORTRAN SUPPLIED ... ABS, LOG10, SQRT
C
C=======================================================================
C
C     Argonne National Laboratory. MINPACK PROJECT, March 1980
C     Burton S. Garbow, Kenneth E. Hillstrom and Jorge J. More
C
C=======================================================================
C
C     epsmch is the machine precision.
C
      EPSMCH=DPMPAR(1) !  EPSMCH is the machine precision
C
      SQREPS=SQRT(EPSMCH)
C
      IF (IFMODE.NE.2) THEN
C
C         IFMODE=1:
C
         DO J=1,LDPARS
            TEMPOR=SQREPS*ABS(XARGUM(J))
            IF (TEMPOR.EQ.R_ZERO) TEMPOR=SQREPS
            XNEIGH(J)=XARGUM(J)+TEMPOR
         END DO
C
         RETURN
C
      END IF
C
C     IFMODE=2:
C
      EPSF=FACTOR*EPSMCH
      EPSLOG=LOG10(EPS)
C
      DO I=1,LDFUNC
         ERRGRA(I)=R_ZERO
      END DO
C
      DO J=1,LDPARS
C
         TEMPOR=ABS(XARGUM(J))
C
         IF (TEMPOR.EQ.R_ZERO) TEMPOR=R_UNIT
C
         DO I=1,LDFUNC
            ERRGRA(I)=ERRGRA(I)+TEMPOR*FJACOB(I,J)
         END DO
C
      END DO
C
      DO I=1,LDFUNC
C
         TEMPOR=R_UNIT
C
         IF (FARGUM(I).NE.R_ZERO.AND.F_VECP(I).NE.R_ZERO.AND.
     *       ABS(F_VECP(I)-FARGUM(I)).GE.EPSF*ABS(FARGUM(I))) THEN
C
             TEMPOR=SQREPS*ABS((F_VECP(I)-FARGUM(I))
     *             /SQREPS-ERRGRA(I))
     *             /(ABS(FARGUM(I))+ABS(F_VECP(I)))
C
         END IF
C
         ERRGRA(I)=R_UNIT
C
         IF (TEMPOR.GT.EPSMCH.AND.TEMPOR.LT.SQREPS) THEN
             ERRGRA(I)=(LOG10(TEMPOR)-EPSLOG)/EPSLOG
         END IF
     
         IF (TEMPOR.GE.SQREPS) ERRGRA(I)=R_ZERO
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
      FUNCTION DPMPAR(I)
C
C     This function provides double precision machine parameters when 
C     the appropriate set of data statements is activated by removing 
C     the c from column 1, and all other data statements are rendered 
C     inactive.  Most of the parameter values  were obtained from the 
C     corresponding Bell Laboratories Port Library function.
C
C     The function statement is
C
C                            FUNCTION DPMPAR(I)
C
C     where
C
C     I - an integer input variable set to 1, 2, or 3,  which selects 
C         the desired machine parameter.  If the machine has t base b 
C         digits, and its smallest and largest exponents are EMIN and 
C         EMAX, respectively, then these parameters are
C
C         DPMPAR(1) = b**(1 - t), the machine precision,
C
C         DPMPAR(2) = b**(emin - 1), the smallest magnitude, and:
C
C         DPMPAR(3) = b**emax*(1 - b**(-t)), the largest magnitude.
C
C     Argonne National Laboratory. MINPACK Project. November 1996.
C     Burton S. Garbow, Kenneth E. Hillstrom and Jorge J. More'
C

      DIMENSION 
     *            MCHEPS(4),MINMAG(4),MAXMAG(4),DMACH(3)
C                 In fact:     double precision DMACH(3)
      EQUIVALENCE 
     *            (DMACH(1),MCHEPS(1))
      EQUIVALENCE 
     *            (DMACH(2),MINMAG(1))
      EQUIVALENCE 
     *            (DMACH(3),MAXMAG(1))
C
C     Machine constants for the IBM 360/370 series,
C     the Amdahl 470/V6, the ICL 2900, the Itel AS/6,
C     the Xerox Sigma 5/7/9 and the Sel systems 85/86.
C
C     data mcheps(1),mcheps(2) / z34100000, z00000000 /
C     data minmag(1),minmag(2) / z00100000, z00000000 /
C     data maxmag(1),maxmag(2) / z7fffffff, zffffffff /
C
C     Machine constants for the Honeywell 600/6000 series.
C
C     data mcheps(1),mcheps(2) / o606400000000, o000000000000 /
C     data minmag(1),minmag(2) / o402400000000, o000000000000 /
C     data maxmag(1),maxmag(2) / o376777777777, o777777777777 /
C
C     Machine constants for the CDC 6000/7000 series.
C
C     data mcheps(1) / 15614000000000000000b /
C     data mcheps(2) / 15010000000000000000b /
C
C     data minmag(1) / 00604000000000000000b /
C     data minmag(2) / 00000000000000000000b /
C
C     data maxmag(1) / 37767777777777777777b /
C     data maxmag(2) / 37167777777777777777b /
C
C     Machine constants for the PDP-10 (KA processor).
C
C     data mcheps(1),mcheps(2) / "114400000000, "000000000000 /
C     data minmag(1),minmag(2) / "033400000000, "000000000000 /
C     data maxmag(1),maxmag(2) / "377777777777, "344777777777 /
C
C     Machine constants for the PDP-10 (KI processor).
C
C     data mcheps(1),mcheps(2) / "104400000000, "000000000000 /
C     data minmag(1),minmag(2) / "000400000000, "000000000000 /
C     data maxmag(1),maxmag(2) / "377777777777, "377777777777 /
C
C     Machine constants for the PDP-11. 
C
C     data mcheps(1),mcheps(2) /   9472,      0 /
C     data mcheps(3),mcheps(4) /      0,      0 /
C
C     data minmag(1),minmag(2) /    128,      0 /
C     data minmag(3),minmag(4) /      0,      0 /
C
C     data maxmag(1),maxmag(2) /  32767,     -1 /
C     data maxmag(3),maxmag(4) /     -1,     -1 /
C
C     Machine constants for the Burroughs 6700/7700 systems.
C
C     data mcheps(1) / o1451000000000000 /
C     data mcheps(2) / o0000000000000000 /
C
C     data minmag(1) / o1771000000000000 /
C     data minmag(2) / o7770000000000000 /
C
C     data maxmag(1) / o0777777777777777 /
C     data maxmag(2) / o7777777777777777 /
C
C     Machine constants for the Burroughs 5700 system.
C
C     data mcheps(1) / o1451000000000000 /
C     data mcheps(2) / o0000000000000000 /
C
C     data minmag(1) / o1771000000000000 /
C     data minmag(2) / o0000000000000000 /
C
C     data maxmag(1) / o0777777777777777 /
C     data maxmag(2) / o0007777777777777 /
C
C     Machine constants for the Burroughs 1700 system.
C
C     data mcheps(1) / zcc6800000 /
C     data mcheps(2) / z000000000 /
C
C     data minmag(1) / zc00800000 /
C     data minmag(2) / z000000000 /
C
C     data maxmag(1) / zdffffffff /
C     data maxmag(2) / zfffffffff /
C
C     Machine constants for the Univac 1100 series.
C
C     data mcheps(1),mcheps(2) / o170640000000, o000000000000 /
C     data minmag(1),minmag(2) / o000040000000, o000000000000 /
C     data maxmag(1),maxmag(2) / o377777777777, o777777777777 /
C
C     Machine constants for the Data General Eclipse S/200.
C
C     Note - it may be appropriate to include the following card -
C     static dmach(3)
C
C     data minmag/20k,3*0/,maxmag/77777k,3*177777k/
C     data mcheps/32020k,3*0/
C
C     Machine constants for the Harris 220.
C
C     data mcheps(1),mcheps(2) / '20000000, '00000334 /
C     data minmag(1),minmag(2) / '20000000, '00000201 /
C     data maxmag(1),maxmag(2) / '37777777, '37777577 /
C
C     Machine constants for the Cray-1.
C
C     data mcheps(1) / 0376424000000000000000b /
C     data mcheps(2) / 0000000000000000000000b /
C
C     data minmag(1) / 0200034000000000000000b /
C     data minmag(2) / 0000000000000000000000b /
C
C     data maxmag(1) / 0577777777777777777777b /
C     data maxmag(2) / 0000007777777777777776b /
C
C     Machine constants for the Prime 400.
C
C     data mcheps(1),mcheps(2) / :10000000000, :00000000123 /
C     data minmag(1),minmag(2) / :10000000000, :00000100000 /
C     data maxmag(1),maxmag(2) / :17777777777, :37777677776 /
C
C     Machine constants for the VAX-11.
C
C     data mcheps(1),mcheps(2) /   9472,  0 /
C     data minmag(1),minmag(2) /    128,  0 /
C     data maxmag(1),maxmag(2) / -32769, -1 /
C
C     Machine constants for IEEE machines.
C
      DATA DMACH(1) /2.22044604926D-016/
      DATA DMACH(2) /2.22507385852D-308/
      DATA DMACH(3) /1.79769313485D+308/
C
C=======================================================================      
C
      DPMPAR=DMACH(I)
C
C=======================================================================      
C
      RETURN
      END 
C
C=======================================================================
C=======================================================================
C
C           WRITINGS FOR GOOD-LOOKING OUTPUTS AND XFIG PLOTTING
C
C=======================================================================
C=======================================================================
C
      SUBROUTINE PRINTI(IDEFCN)
      
      INCLUDE  'MATDIM/NDPARS.f'

      CHARACTER
     *          NUCLEO*8
      DIMENSION
     *          NUCLEO(0:1)
C
      COMMON
     *       /PRFILE/ ISCREN,LOGFIL,LOGCHP,LOGCHN,LOGHIP,LOGHIN,
     *                LOGPRO,LOGNEU,LOGWRI,LSCREN,IRESUL,ICONVE,
     *                NOUTPT,LOGAUX,LOGENE,LOGRAD,LOGGRD,LOGMSH,
     *                                                   LOGBIS                                                  
      COMMON
     *       /NUCLEU/ IZ_FIX,IN_FIX,ISOSPI
C
      COMMON
     *       /POTPOT/ PARPOT(1:NDPARS)
C
      DATA      
     *       NUCLEO(0) / 'Neutrons' /         
      DATA
     *       NUCLEO(1) / ' Protons' /         
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
      IF (ISOSPI.EQ.1 .AND. LOGWRI.EQ.2) THEN
C      
         WRITE(IRESUL,'(/,80(''#''),//,80(''#''),/,''#'',t80,''#'',/,
     *             ''#   It'',i4,t69,a8,t80,''#'',/,
     *             ''#   V0CENT R0CENT A0CENT'',
     *             ''  V0SORB R0SORB A0SORB'',
     *             ''  V0EFFM R0EFFM A0EFFM  R0COUL   #'')') IDEFCN,
     *                                                NUCLEO(ISOSPI)
C
         WRITE(IRESUL,'(''#'',F9.3,2F7.4,2(F8.3,2F7.4),1x,f7.4,t80,
     *                  ''#'')')
     *               V0CENT_PROTON,R0CENT_PROTON,A0CENT_PROTON,
     *               V0SORB_PROTON,R0SORB_PROTON,A0SORB_PROTON,
     *               V0EFFM_PROTON,R0EFFM_PROTON,A0EFFM_PROTON,R0COUL
C
         WRITE(IRESUL,'(''#   XK_V0C XK_R0C XK_A0C'',
     *                  ''  XK_LAM XK_RSO XK_ASO'',
     *                  ''  XK_LEF XK_REF XK_AEF  XK_COU   #'')')
C
         WRITE(IRESUL,'(''#'',F9.3,2F7.4,2(F8.3,2F7.3),1x,f7.3,t80,
     *                  ''#'')')          
     *               XK_V0C_PROTON,XK_R0C_PROTON,XK_A0C_PROTON,
     *               XK_LAM_PROTON,XK_RSO_PROTON,XK_ASO_PROTON,
     *               XK_LEF_PROTON,XK_REF_PROTON,XK_AEF_PROTON,XK_COU
C
         WRITE(IRESUL,'(''#'',T80,''#'')')     
C
      END IF
      
      IF (ISOSPI.EQ.0 .AND. LOGWRI.EQ.2) THEN
C      
         WRITE(IRESUL,'(/,80(''#''),//,80(''#''),/,''#'',t80,''#'',/,
     *             ''#   It'',i4,t69,a8,t80,''#'',/,
     *             ''#   V0CENT R0CENT A0CENT'',
     *             ''  V0SORB R0SORB A0SORB'',
     *             ''  V0EFFM R0EFFM A0EFFM  R0COUL   #'')') IDEFCN,
     *                                                NUCLEO(ISOSPI)
C
         WRITE(IRESUL,'(''#'',F9.3,2F7.4,2(F8.3,2F7.4),1x,f7.4,t80,
     *                  ''#'')')
     *               V0CENT_NEUTRS,R0CENT_NEUTRS,A0CENT_NEUTRS,
     *               V0SORB_NEUTRS,R0SORB_NEUTRS,A0SORB_NEUTRS,
     *               V0EFFM_NEUTRS,R0EFFM_NEUTRS,A0EFFM_NEUTRS,R0COUL
C
         WRITE(IRESUL,'(''#   XK_V0C XK_R0C XK_A0C'',
     *             ''  XK_LAM XK_RSO XK_ASO'',
     *             ''  XK_LEF XK_REF XK_AEF  XK_COU   #'')')
C
         WRITE(IRESUL,'(''#'',F9.3,2F7.4,2(F8.3,2F7.3),1x,f7.3,t80,
     *                  ''#'')')          
     *               XK_V0C_NEUTRS,XK_R0C_NEUTRS,XK_A0C_NEUTRS,
     *               XK_LAM_NEUTRS,XK_RSO_NEUTRS,XK_ASO_NEUTRS,
     *               XK_LEF_NEUTRS,XK_REF_NEUTRS,XK_AEF_NEUTRS,XK_COU
C
         WRITE(IRESUL,'(''#'',T80,''#'')')
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
      SUBROUTINE INPRIN(IDEFCN,IACTIV,I_MODE,
     *                         CHISQU_PROTON,CHISQU_NEUTRS,
     *                         CHIENE_WEIPRO,CHIENE_WEINEU,
     *                         CHIRAD_WEIPRO,CHIRAD_WEINEU,
     *                                       CHIRHO_WEIPRO)
      
      INCLUDE  'MATDIM/NDPARS.f'
      INCLUDE  'MATDIM/NDTITL.f'
      INCLUDE  'MATDIM/NDNUCL.f'
      
      CHARACTER
     *          FROMXX*1,ACTION*1,TITLES*12,STRING_PROTON*8,
     *                                      STRING_NEUTRS*8
      CHARACTER
     *          TITLES_LATEXS*100,NUCNAM_LATEXS*010
      DIMENSION
     *          FROMXX(0:2)
      DIMENSION
     *          IAUXIL(1:NDPARS)
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
     *       /PRFILE/ ISCREN,LOGFIL,LOGCHP,LOGCHN,LOGHIP,LOGHIN,
     *                LOGPRO,LOGNEU,LOGWRI,LSCREN,IRESUL,ICONVE,
     *                NOUTPT,LOGAUX,LOGENE,LOGRAD,LOGGRD,LOGMSH,
     *                                                   LOGBIS                                                  
      COMMON
     *       /OCCREC/ N_ACTU_PROTON,N_CORR_PROTON,
     *                N_ACTU_NEUTRS,N_CORR_NEUTRS
      COMMON
     *       /POTPOT/ PARPOT(1:NDPARS)
      COMMON
     *       /ACTACT/ ACTION(1:NDPARS)
      COMMON
     *       /SEARCH/ IFTAKE(1:NDPARS)
CID      COMMON
CID     *       /GSTORE/ DERPRI(1:NDPARS) 
      COMMON
     *       /CALFRO/ I_FROM
      COMMON
     *       /CHICHO/ VALCHI
      COMMON
     *       /DENSTR/ IFDENS
      COMMON
     *       /TENSOR/ IFTENS
      COMMON
     *       /TENOPT/ ICENTT,ISORBT,ITENSR
      COMMON
     *       /IFKAPP/ IFK_VC,IFK_VS,IFK_RC,IFK_RS,IFK_AC,IFK_AS
      COMMON
     *       /PARAMT_TITLES/ TITLES(1:NDTITL),
     *                       TITLES_LATEXS(1:NDTITL),
     *                       NUCNAM_LATEXS(1:NDNUCL)
C
      DATA
     *       FROMXX(00) / ':' /         
      DATA
     *       FROMXX(01) / 'f' /         
      DATA
     *       FROMXX(02) / 'g' / 
C
C=======================================================================
C     This subroutine makes a dynamical construction of the output
C     line on the user screen showing the advancement of the mini-
C     misation in function of the processing iterations 
C=======================================================================
C
      LOGWRI_AUXILI=1
      IF (LOGWRI.EQ.0) THEN
          LOGWRI_AUXILI=LOGWRI
          LOGWRI=1 !Irene says: I want this kind of files always...
      END IF
C
C=======================================================================
C
      IFPROT=0
      IFNEUT=0
      IFBOTH=0
C      
      DO IPARAM=1,NDPARS
         IF (IFTAKE(IPARAM).EQ.1) THEN
             IF (IPARAM.LE.20.AND.IFDENS.EQ.0) THEN
                 IFPROT=1
             END IF
             IF (IPARAM.GT.20.AND.IFDENS.EQ.0) THEN
                 IFNEUT=1
             END IF
             IF (IPARAM.GE.51) THEN
                 IFPROT=0
                 IFNEUT=0
                 IFBOTH=1
             END IF
         END IF
      END DO   
C
      IF ((IFDENS.EQ.0).AND.(IFNEUT.EQ.1).AND.(IFPROT.EQ.1)) THEN
C  
           IFPROT=0
           IFNEUT=0
           IFBOTH=1
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
C     Titles of columns are printed every  20  active output lines
C
      IF ((MOD(IDEFCN,20).EQ.1.AND.I_MODE.EQ.0).OR.I_MODE.EQ.2) THEN
C
CID          WRITE(ICONVE,'(/,'' Iter '',$)')
C
CID          DO I=1,NDPARS
C             IF (IFTAKE(I).EQ.1) THEN
C                 WRITE(ICONVE,'(A10,$)') TITLES(I)
C             END IF
C          END DO
C
C          WRITE(ICONVE,'(''   ChiSquar'')')

C_______________________________________________________________________
C
          IF (LOGWRI.GT.0) THEN
              
              WRITE(LOGAUX,'(/,'' Iter   '',$)')
C
              IF (I_MODE.EQ.2) THEN
C
                  IF (IFDENS.EQ.0 .AND. IFBOTH.EQ.0) THEN
                      WRITE(LOGAUX,'(1X,''r.m.s.-final'',6x,
     *                                  ''ChiSquar-final'',6X,
     *                                  ''Chi2-EnerFinal'',6X,
     *                                  ''Chi2-RadiFinal'',6X,
     *                                  '' Chi2-RhoFinal'',''    '',$)')
                  END IF
                  IF (IFDENS.EQ.1 .OR. 
     *               (IFDENS.EQ.0 .AND. IFBOTH.EQ.1)) THEN 
                      WRITE(LOGAUX,'(5X,''r.m.s.--{P,N}-final'',9x,
     *                                  ''ChiSquar{P,N}-final'',9X,
     *                                  ''ChiSqrtE{P,N}-final'',9X,
     *                                  ''ChiSqrtR{P,N}-final'',5x,
     *                                                   '' '',$)')
                  END IF
C
              ELSE
C
                  IF (IFDENS.EQ.0 .AND. IFPROT.EQ.1) THEN
                      WRITE(LOGAUX,'(1X,''r.m.s.-value'',2x,
     *                                  ''ChiSquar'',6X,
     *                                  ''ChiSqrt-Ener'',2X,
     *                                  ''ChiSqrt-Radi'',2X,
     *                                  ''ChiSqrt-Rho'',2X,''   '',$)')
                  END IF
C
                  IF (IFDENS.EQ.0 .AND. IFNEUT.EQ.1) THEN
                      WRITE(LOGAUX,'(1X,''r.m.s.-value'',2x,
     *                              ''ChiSquar'',6X,
     *                              ''ChiSqrt-Ener'',2X,
     *                              ''ChiSqrt-Radi'',2X,''   '',$)')
                  END IF
C
                  IF (IFDENS.EQ.1 .OR. 
     *               (IFDENS.EQ.0 .AND. IFBOTH.EQ.1)) THEN 
                      WRITE(LOGAUX,'(7X,''r.m.s.-{P,N}'',2x,
     *                              ''ChiSquar{P,N}'',17X,
     *                              ''ChiSqrtE{P,N}'',17X,
     *                              ''ChiSqrtR{P,N}'',5X,''  '',$)')
                  END IF
C
              END IF
C
              DO I=1,INDEXP
                 IPARAM=IAUXIL(I)
                 WRITE(LOGAUX,'(A11,$)') TITLES(IPARAM)
              END DO
C
              IF (I_MODE.NE.2) THEN
              
                  IF (IFDENS.EQ.0 .AND. IFBOTH.EQ.0) THEN
                      WRITE(LOGAUX,'(''        '',$)')
                  END IF
C
                  IF (IFDENS.EQ.1 .OR. 
     *               (IFDENS.EQ.0 .AND. IFBOTH.EQ.1)) THEN  
                      WRITE(LOGAUX,'(18X,''    '',$)')
                  END IF
C
                  DO I=1,INDEXP
                     IPARAM=IAUXIL(I)
                     WRITE(LOGAUX,'(2X,''D'',A8,$)') 
     *                                            TITLES(IPARAM)(3:10)
                  END DO
C
                  WRITE(LOGAUX,'(''     Gradient_Norm'')')
C
              END IF
C
              IF (I_MODE.EQ.2) WRITE(LOGAUX,'()')
C
          END IF
C_______________________________________________________________________
C
          IF (ISCREN.NE.1) GO TO 1
C          
          WRITE(LSCREN,'(/,'' Iter   '',$)')
C
          DO I=1,NDPARS
             IF (IFTAKE(I).EQ.1) THEN
                 WRITE(LSCREN,'(A11,$)') TITLES(I)
             END IF
          END DO
C
          IF (I_MODE.EQ.2) THEN
              WRITE(LSCREN,'(''   ChiSquar{P,N}-final'')')
          ELSE
              WRITE(LSCREN,'(''      ChiSquar{P,N}'',/)')
          END IF
C
   1      CONTINUE          
C
      END IF
C
C=======================================================================
C
      IF (N_ACTU_PROTON.NE.N_CORR_PROTON) THEN
          WRITE(STRING_PROTON,'(I4,'':'',I3)') 
     *          N_ACTU_PROTON,N_CORR_PROTON
      ELSE
          WRITE(STRING_PROTON,'(''        '')') 
      END IF
C
      IF (N_ACTU_NEUTRS.NE.N_CORR_NEUTRS) THEN
          WRITE(STRING_NEUTRS,'(I4,'':'',I3)') 
     *          N_ACTU_NEUTRS,N_CORR_NEUTRS
      ELSE
          WRITE(STRING_NEUTRS,'(''        '')') 
      END IF
C
C=======================================================================
C
CID      WRITE(ICONVE,'(i4,A1,$)') IDEFCN,FROMXX(I_FROM)
C
C      DO I=1,NDPARS
C         IF (IFTAKE(I).EQ.1) THEN
C             WRITE(ICONVE,'(F9.4,A1,$)') PARPOT(I),ACTION(I)
C         END IF
C      END DO
C
C      WRITE(ICONVE,'(F16.4,A8,A8)') VALCHI,STRING_PROTON,STRING_NEUTRS 
C
C=======================================================================
C                            Constructing dynamically the output lines
      IF (LOGWRI.GT.0) THEN
          
        IF (I_MODE.EQ.0) THEN
            
          WRITE(LOGAUX,'(I4,A1,''  '',$)') IDEFCN,FROMXX(I_FROM)
            
          IF (IFDENS.EQ.0 .AND. IFPROT.EQ.1) THEN
            WRITE(LOGAUX,'(5E14.6,''    '',$)') SQRT(CHIWEI_PROTON),
     *                                          CHISQU_PROTON,
     *                                          CHIENE_WEIPRO,
     *                                          CHIRAD_WEIPRO,
     *                                          CHIRHO_WEIPRO
            GO TO 3
          END IF
C          
          IF (IFDENS.EQ.0 .AND. IFNEUT.EQ.1) THEN
            WRITE(LOGAUX,'(4E14.6,''    '',$)')  SQRT(CHIWEI_NEUTRS),
     *                                           CHISQU_NEUTRS,
     *                                           CHIENE_WEINEU,
     *                                           CHIRAD_WEINEU
            GO TO 3
          END IF
C         
          WRITE(LOGAUX,'(8E14.6,''    '',$)') SQRT(CHIWEI_PROTON),
     *                                        SQRT(CHIWEI_NEUTRS),
     *                                        CHISQU_PROTON,
     *                                        CHISQU_NEUTRS,
     *                                        CHIENE_WEIPRO,
     *                                        CHIENE_WEINEU,
     *                                        CHIRAD_WEIPRO,
     *                                        CHIRAD_WEINEU
          
   3      CONTINUE     
        
        END IF
        
        IF (I_MODE.EQ.2) THEN
            
          WRITE(LOGAUX,'(I4,A1,''  '',$)') IDEFCN,FROMXX(I_FROM)
            
          IF (IFDENS.EQ.0.AND.IFPROT.EQ.1) THEN
              WRITE(LOGAUX,'(2E16.6,3E20.6,''    '',$)')
     *                                             SQRT(CHIWEI_PROTON),
     *                                                   CHISQU_PROTON,
     *                                                   CHIENE_WEIPRO,
     *                                                   CHIRAD_WEIPRO,
     *                                                   CHIRHO_WEIPRO
            GO TO 4
          END IF
C          
          IF (IFDENS.EQ.0.AND.IFNEUT.EQ.1) THEN
            WRITE(LOGAUX,'(2E16.6,2E20.6,''    '',$)')
     *                                           SQRT(CHIWEI_NEUTRS),
     *                                                 CHISQU_NEUTRS,
     *                                                 CHIENE_WEINEU,
     *                                                 CHIRAD_WEINEU
            GO TO 4
          END IF
C         
          WRITE(LOGAUX,'(8E14.6,''    '',$)') SQRT(CHIWEI_PROTON),
     *                                        SQRT(CHIWEI_NEUTRS),
     *                                             CHISQU_PROTON,
     *                                             CHISQU_NEUTRS,
     *                                             CHIENE_WEIPRO,
     *                                             CHIENE_WEINEU,
     *                                             CHIRAD_WEIPRO,
     *                                             CHIRAD_WEINEU
          
   4      CONTINUE     
        
        END IF
C_______________________________________________________________________  
C      
        IF (I_MODE.NE.1) THEN
C
            JSTORE=0
          
            WRITE(LOGAUX,'(7X,$)')
C
             DO I=1,INDEXP
                IPARAM=IAUXIL(I)
                WRITE(LOGAUX,'(F10.5,A1,$)') PARPOT(IPARAM),
     *                                       ACTION(IPARAM)
             END DO
C
        END IF
C_______________________________________________________________________  
C
        IF (I_MODE.NE.1) THEN
            
            IF (IFDENS.EQ.0 .AND. IFPROT.EQ.1) THEN
                WRITE(LOGAUX,'(A8,A8)') STRING_PROTON,STRING_NEUTRS
                GO TO 5                    
            END IF
C          
            IF (IFDENS.EQ.0.AND.IFNEUT.EQ.1) THEN
                WRITE(LOGAUX,'(A8,A8)') STRING_PROTON,STRING_NEUTRS
                GO TO 5
            END IF
C         
            WRITE(LOGAUX,'(A8,A8)') STRING_PROTON,STRING_NEUTRS
          
   5        CONTINUE     
        
        END IF
      
      END IF
C
C=======================================================================
C
      IF (ISCREN.NE.1) THEN
          GO TO 2
      END IF
C
C=======================================================================
C                            Constructing dynamically the output lines
      IF (I_MODE.NE.1) THEN
C
          WRITE(0,'(I4,A1,''  '',$)') IDEFCN,FROMXX(I_FROM)
C
          JSTORE=0
C
          DO I=1,NDPARS
             IF (IFTAKE(I).EQ.1) THEN
                 WRITE(0,'(F10.4,A1,$)') PARPOT(I),ACTION(I)
             END IF
          END DO
C
      END IF
C_______________________________________________________________________
C
C     IF (IACTIV.NE.0.AND.I_MODE.EQ.1.AND.I_FROM.EQ.2) THEN
C
C         WRITE(0,'(F16.4,F10.4,A8,A8)') 
C    *         VALCHI,DERPRI(IACTIV),STRING_PROTON,STRING_NEUTRS
c      ELSE
c          WRITE(0,'(F16.4,A8,A8)') VALCHI,STRING_PROTON,STRING_NEUTRS
C     END IF  
C
      IF (I_MODE.NE.1) THEN
C         WRITE(0,'(F10.4,A8,A8)') VALCHI,STRING_PROTON,STRING_NEUTRS
          WRITE(0,'(2F10.4,A8,A8)') CHISQU_PROTON,CHISQU_NEUTRS,
     *                              STRING_PROTON,STRING_NEUTRS
      END IF
C
   2  CONTINUE
C
C=======================================================================
C
      IF (LOGWRI_AUXILI.EQ.0) THEN
          LOGWRI=0
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
      SUBROUTINE INPRIN_GLOFIT(IDEFCN,IACTIV,I_MODE,CHIENE_PROTON,
     *                                CHIENE_NEUTRS,CHIRAD_PROTON,
     *                                              CHIRAD_NEUTRS)
C
      INCLUDE  'MATDIM/NDPARS.f'
      INCLUDE  'MATDIM/NDNUCL.f'
      INCLUDE  'MATDIM/NDTITL.f'
C
      CHARACTER
     *          FROMXX*1,ACTION*1,TITLES*12,STRING_PROTON*8,
     *                                      STRING_NEUTRS*8
      CHARACTER
     *          CHITIT*12,NUCSYM*6,INPSYM*6,RADTIT*12,
     *          TITLES_LATEXS*050,NUCNAM_LATEXS*010
C
      DIMENSION
     *          CHITIT(1:010),RADTIT(1:010)
      DIMENSION
     *          FROMXX(0:2)
      DIMENSION
     *          IAUXIL(1:NDPARS)
      DIMENSION
     *          CHIENE_PROTON(1:NDNUCL),
     *          CHIENE_NEUTRS(1:NDNUCL)
      DIMENSION
     *          CHIRAD_PROTON(1:NDNUCL),
     *          CHIRAD_NEUTRS(1:NDNUCL)
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
     *       /PRFILE/ ISCREN,LOGFIL,LOGCHP,LOGCHN,LOGHIP,LOGHIN,
     *                LOGPRO,LOGNEU,LOGWRI,LSCREN,IRESUL,ICONVE,
     *                NOUTPT,LOGAUX,LOGENE,LOGRAD,LOGGRD,LOGMSH,
     *                                                   LOGBIS 
      COMMON
     *       /POTPOT/ PARPOT(1:NDPARS)
      COMMON
     *       /ACTACT/ ACTION(1:NDPARS)
      COMMON
     *       /SEARCH/ IFTAKE(1:NDPARS)
      COMMON
     *       /CALFRO/ I_FROM
      COMMON
     *       /DENSTR/ IFDENS
      COMMON
     *       /TENSOR/ IFTENS
      COMMON
     *       /TENOPT/ ICENTT,ISORBT,ITENSR
      COMMON
     *       /IFKAPP/ IFK_VC,IFK_VS,IFK_RC,IFK_RS,IFK_AC,IFK_AS
      COMMON
     *       /PARAMT_TITLES/ TITLES(1:NDTITL),
     *                       TITLES_LATEXS(1:NDTITL),
     *                       NUCNAM_LATEXS(1:NDNUCL)
C
      DATA
     *       FROMXX(00) / ':' /         
      DATA
     *       FROMXX(01) / 'f' /         
      DATA
     *       FROMXX(02) / 'g' / 
      DATA
     *       CHITIT(01) / '  RMSENE_Ox ' /,
     *       CHITIT(02) / '  RMSENE_Ca ' /,
     *       CHITIT(03) / '  RMSENE_Ca ' /,
     *       CHITIT(04) / '  RMSENE_Ni ' /,
     *       CHITIT(05) / '  RMSENE_Zr ' /,
     *       CHITIT(06) / '  RMSENE_Sn ' /,
     *       CHITIT(07) / '  RMSENE_Gd ' /,
     *       CHITIT(08) / '  RMSENE_Pb ' /
      DATA
     *       RADTIT(01) / '  DIFRAD_Ox ' /,
     *       RADTIT(02) / '  DIFRAD_Ca ' /,
     *       RADTIT(03) / '  DIFRAD_Ca ' /,
     *       RADTIT(04) / '  DIFRAD_Ni ' /,
     *       RADTIT(05) / '  DIFRAD_Zr ' /,
     *       RADTIT(06) / '  DIFRAD_Sn ' /,
     *       RADTIT(07) / '  DIFRAD_Gd ' /,
     *       RADTIT(08) / '  DIFRAD_Pb ' /
C
C=======================================================================
C     This subroutine makes a dynamical construction of the output
C     line on the user screen showing the advancement of the mini-
C     misation in function of the processing iterations 
C=======================================================================
C
      LOGWRI_AUXILI=1
      IF (LOGWRI.EQ.0) THEN
          LOGWRI_AUXILI=LOGWRI
          LOGWRI=1 !Irene says: I want this kind of files always...
      END IF
C
C=======================================================================
C
      IFPROT=0
      IFNEUT=0
      IFBOTH=0
C      
      DO IPARAM=1,NDPARS
         IF (IFTAKE(IPARAM).EQ.1) THEN
             IF (IPARAM.LE.20 .AND. IFDENS.EQ.0) THEN
                 IFPROT=1
             END IF
             IF (IPARAM.GT.20 .AND. IFDENS.EQ.0) THEN
                 IFNEUT=1
             END IF
             IF (IPARAM.GE.51 .AND. IFDENS.EQ.0) THEN
                 IFPROT=0
                 IFNEUT=0
                 IFBOTH=1
             END IF
         END IF
      END DO  
C
      IF ((IFDENS.EQ.0).AND.(IFNEUT.EQ.1).AND.(IFPROT.EQ.1)) THEN
           IFPROT=0
           IFNEUT=0
           IFBOTH=1
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
C     Titles of columns are printed every  20  active output lines
C
      IF ((MOD(IDEFCN,20).EQ.1.AND.I_MODE.EQ.0).OR.I_MODE.EQ.2) THEN
C
          IF (LOGWRI.GT.0) THEN
              
              WRITE(LOGENE,'(/,/,'' Iter   '',$)')
              WRITE(LOGRAD,'(/,/,'' Iter   '',$)')
C
              DO I=1,INDEXP
                 IPARAM=IAUXIL(I)
                 WRITE(LOGENE,'(A14,$)') TITLES(IPARAM)
                 WRITE(LOGRAD,'(A14,$)') TITLES(IPARAM)
              END DO
            
              DO I=1,LDNUCL
                 IF (ITAKNU(I).EQ.1) THEN
                     IF (IFDENS.EQ.0 .AND. IFBOTH.EQ.0) THEN
                         WRITE(LOGENE,'(A14,$)') CHITIT(I)
                         WRITE(LOGRAD,'(A14,$)') RADTIT(I)
                     END IF
                     IF ((IFDENS.EQ.1) .OR.  
     *                   (IFDENS.EQ.0 .AND. IFBOTH.EQ.1))THEN
                         WRITE(LOGENE,'(8X,A9,''{P,N}'',8X,''  '',$)') 
     *                                            CHITIT(I)(3:11)
                         WRITE(LOGRAD,'(8X,A9,''{P,N}'',8X,''  '',$)') 
     *                                            RADTIT(I)(3:11)
                     END IF
                 END IF
              END DO
C
              IF (I_MODE.NE.2) THEN
C
                  IF (IFDENS.EQ.0 .AND. IFBOTH.EQ.0) THEN
                      WRITE(LOGENE,'(''    '',$)')
                      WRITE(LOGRAD,'(''    '',$)')
                  END IF
C
                  IF ((IFDENS.EQ.1) .OR. 
     *                (IFDENS.EQ.0 .AND. IFBOTH.EQ.1))THEN
                    WRITE(LOGENE,'(18X,''    '',$)')
                    WRITE(LOGRAD,'(18X,''    '',$)')
                  END IF
C
              END IF
C
          END IF          
C
      END IF 
C
      WRITE(LOGENE,'()')
      WRITE(LOGRAD,'()')
C
C=======================================================================
C                            Constructing dynamically the output lines
      IF (LOGWRI.GT.0) THEN
C_______________________________________________________________________  
C          
       IF (I_MODE.NE.1) THEN
           
           WRITE(LOGENE,'(I4,A1,''  '',$)') IDEFCN,FROMXX(I_FROM)
           WRITE(LOGRAD,'(I4,A1,''  '',$)') IDEFCN,FROMXX(I_FROM)
C 
           JSTORE=0
          
C          WRITE(LOGENE,'(7X,$)')
C
           DO I=1,INDEXP
              IPARAM=IAUXIL(I)
              WRITE(LOGENE,'(F13.5,A1,$)') PARPOT(IPARAM),
     *                                     ACTION(IPARAM)
              WRITE(LOGRAD,'(F13.5,A1,$)') PARPOT(IPARAM),
     *                                     ACTION(IPARAM)
           END DO
C
        END IF
C_______________________________________________________________________  
C          
        IF (I_MODE.NE.1) THEN
C
            DO INUCLI=1,LDNUCL
C
               IF (ITAKNU(INUCLI).EQ.1) THEN
C
                   IF (IFDENS.EQ.0 .AND. IFPROT.EQ.1) THEN
                       WRITE(LOGENE,'(E14.4,$)') CHIENE_PROTON(INUCLI)
                       WRITE(LOGRAD,'(E14.4,$)') CHIRAD_PROTON(INUCLI)
                   END IF
C          
                   IF (IFDENS.EQ.0 .AND. IFNEUT.EQ.1) THEN
                       WRITE(LOGENE,'(E14.4,$)') CHIENE_NEUTRS(INUCLI)
                       WRITE(LOGRAD,'(E14.4,$)') CHIRAD_NEUTRS(INUCLI)
                   END IF
C
                   IF ((IFDENS.EQ.1)
     *            .OR. (IFDENS.EQ.0 .AND. IFBOTH.EQ.1)) THEN
C         
                       WRITE(LOGENE,'(2E14.4,''    '',$)') 
     *                   CHIENE_PROTON(INUCLI),CHIENE_NEUTRS(INUCLI)
C
                      WRITE(LOGRAD,'(2E14.4,''    '',$)') 
     *                   CHIRAD_PROTON(INUCLI),CHIRAD_NEUTRS(INUCLI)
C
                   END IF
C
               END IF
C
            END DO
C
        END IF
C_______________________________________________________________________  
C      
      END IF
C
C=======================================================================
C
      IF (LOGWRI_AUXILI.EQ.0) THEN
          LOGWRI=0
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
      SUBROUTINE RADINF(IFPROT,IFNEUT)
C   
      INCLUDE  'MATDIM/NDNUCL.f'
C
      CHARACTER
     *          NUCSYM*6,INPSYM*6
C
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
C      COMMON
C     *       /P_N_LM/ IFPROT,IFNEUT,IFBOTH
      COMMON
     *       /MODIFR/ RMSEXP_PROTON(1:NDNUCL),
     *                RMSEXP_NEUTRS(1:NDNUCL),
     *                RMSTHE_PROTON,RMSTHE_NEUTRS,
     *                RMODIF_PROTON,RMODIF_NEUTRS,
     *                              PUSHIN,DELTAR
      COMMON
     *       /NUCWEI/ WEINUC_PROTON(1:NDNUCL),
     *                WEINUC_NEUTRS(1:NDNUCL),
     *                WEIRAD_PROTON(1:NDNUCL),
     *                WEIRAD_NEUTRS(1:NDNUCL)
      COMMON
     *       /DENSTR/ IFDENS
      COMMON
     *       /PRFILE/ ISCREN,LOGFIL,LOGCHP,LOGCHN,LOGHIP,LOGHIN,
     *                LOGPRO,LOGNEU,LOGWRI,LSCREN,IRESUL,ICONVE,
     *                NOUTPT,LOGAUX,LOGENE,LOGRAD,LOGGRD,LOGMSH,
     *                                                   LOGBIS    
C
C=======================================================================
C      
      DO INUCLI=1,LDNUCL
          
        IF (ITAKNU(INUCLI).EQ.1) THEN
      
        IF (IFDENS.EQ.0 .AND. IFPROT.EQ.1) THEN
          
          WRITE(LOGAUX,'(/,A6,'' --> Radii Info - Protons:'',
     *                           '' Rex= '',f6.4,'' fm ;  '',
     *                           '' Rth= '',f6.4,'' fm ;  '',
     *                           '' Rwt= '',f6.2,'' ;  '',
     *                           '' Dif= '',e10.4,'' fm'',/)')
     *                                          NUCSYM(INUCLI),
     *                                   RMSEXP_PROTON(INUCLI),
     *                                           RMSTHE_PROTON,
     *                                   WEIRAD_PROTON(INUCLI),
     *                     RMSEXP_PROTON(INUCLI)-RMSTHE_PROTON
        END IF
      
        IF (IFDENS.EQ.0.AND.IFNEUT.EQ.1) THEN
C
          WRITE(LOGAUX,'(/,A6,'' --> Radii Info - Neutrons:'',
     *                            '' Rex= '',f6.4,'' fm ;  '',
     *                            '' Rth= '',f6.4,'' fm ;  '',
     *                            '' Rwt= '',f6.2,'' ;  '',
     *                            '' Dif= '',e10.4,'' fm'',/)')
     *                                          NUCSYM(INUCLI),
     *                                   RMSEXP_NEUTRS(INUCLI),
     *                                           RMSTHE_NEUTRS,
     *                                   WEIRAD_NEUTRS(INUCLI),
     *                     RMSEXP_NEUTRS(INUCLI)-RMSTHE_NEUTRS
        END IF
      
        IF (IFDENS.EQ.1) THEN
C
          WRITE(LOGAUX,'(/,A6,'' --> Radii Info - Protons:'',
     *                           '' Rex= '',f6.4,'' fm ;  '',
     *                           '' Rth= '',f6.4,'' fm ;  '',
     *                           '' Rwt= '',f6.2,'' ;  '',
     *                           '' Dif= '',e10.4,'' fm'',/)')
     *                                          NUCSYM(INUCLI),
     *                                   RMSEXP_PROTON(INUCLI),
     *                                           RMSTHE_PROTON,
     *                                   WEIRAD_PROTON(INUCLI),
     *                     RMSEXP_PROTON(INUCLI)-RMSTHE_PROTON
C
          WRITE(LOGAUX,'(/,A6,'' --> Radii Info - Neutrons:'',
     *                            '' Rex= '',f6.4,'' fm ;  '',
     *                            '' Rth= '',f6.4,'' fm ;  '',
     *                            '' Rwt= '',f6.2,'' ;  '',
     *                            '' Dif= '',e10.4,'' fm'',/)')
     *                                          NUCSYM(INUCLI),
     *                                   RMSEXP_NEUTRS(INUCLI),
     *                                           RMSTHE_NEUTRS,
     *                                   WEIRAD_NEUTRS(INUCLI),
     *                     RMSEXP_NEUTRS(INUCLI)-RMSTHE_NEUTRS
          
        END IF
        
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
      SUBROUTINE WRITIN_ENELEV(ISOSPI,IRANDO_PRINTG,LDRAND,
     *                         IEVALU_PRINTG,PARPOT_PRINTG,
     *                         CHISQU_PRINTG,CHIGRD_PRINTG,
     *                         ERRMAX_PRINTG,RMSVAL_PRINTG,
     *                         RMSGLO_PRINTG,RMSEXP_PRINTG,
     *                         RMSTHE_PRINTG,LEVTHE_PRINTG,
     *                         ENETHE_PRINTG,LABTHE_PRINTG,
     *                                              INUCLI)
c          
      INCLUDE  'MATDIM/NDNUCL.f'
      INCLUDE  'MATDIM/NDPARS.f'
      INCLUDE  'MATDIM/NDLEXP.f'
      INCLUDE  'MATDIM/NDSPEC.f'
      INCLUDE  'MATDIM/NDIM_M.f'
      INCLUDE  'MATDIM/NDCOLO.f'
      INCLUDE  'MATDIM/NDTITL.f'
c      
      CHARACTER
     *          FILNAM*256,TITPAR*13,AUXLAB*6,LABTEX*11,TITLES*12,
     *          FILNAM_ENDING*7,NUCNAM*256,TEXLAM*20
      CHARACTER
     *          FILNAM_IFITED*14,FITNUC*2,FILNAM_NUCFIT*30,
     *          FITNUC_AUXILI*2,TITLES_LATEXS*050,
     *          NUCNAM_LATEXS*010
      CHARACTER
     *          INPSYM*6,TYPCHI*6,NUCSYM*6,VERSIO*3
      CHARACTER
     *          LABEXP*6,LABTHE*6,LABTHE_PRINTG*6
      CHARACTER
     *          LABTHE_PROTON*6,LABTHE_NEUTRS*6
      CHARACTER
     *          LABEXP_PROTON*6,LABEXP_NEUTRS*6
      DIMENSION
     *          ENETHE(1:NDSPEC),LABTHE(1:NDSPEC)
      DIMENSION
     *          EXPEXP(1:NDLEXP),IDEGEX(1:NDLEXP)
      DIMENSION
     *          LABEXP(1:NDLEXP)
      DIMENSION
     *          ENETHE_PRINTG(1:NDSPEC),LABTHE_PRINTG(1:NDSPEC),
     *          PARPOT_PRINTG(1:NDPARS)
      DIMENSION
     *          ARGPAR(1:NDPARS)
      DIMENSION
     *          TITPAR(1:NDPARS),
     *          PARPOT_AUXILI(1:NDPARS)
      DIMENSION
     *          FITNUC(1:NDNUCL),
     *          FITNUC_AUXILI(1:NDNUCL)
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
     *       /ER_RMS/ RMSERR_PROTON(1:NDNUCL),
     *                RMSERR_NEUTRS(1:NDNUCL)
      COMMON 
     *       /RMSVAL/ HINORM_PROTON,HINORM_NEUTRS
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
      COMMON
     *       /DEPPAI/ IFDEEP,IFPRON
      COMMON
     *       /CHI_PN/ CHISQU_PROTON,CHISQU_NEUTRS
      COMMON
     *       /IFKAPP/ IFK_VC,IFK_VS,IFK_RC,IFK_RS,IFK_AC,IFK_AS
      COMMON
     *       /XLENGT/ XMIN_T,XMAX_T,XMIN_E,XMAX_E
      COMMON
     *       /ATTRIB/ ISTYLE,I_TYPE,ICOLOR(1:NDCOLO),ITHICK(1:NDCOLO)
      COMMON 
     *       /CHITAK/ IF_SPE,IF_RAD,IF_GAP,IF_FER,IF_DEN,IF_RHO,IF_INV,
     *                LDCHI2,TYPCHI(1:NDNUCL),ITKCHI(1:NDNUCL)
      COMMON
     *       /IFLAMB/ IFPAR1,IFPAR2,IFPAR3,IFPAR4,
     *                IFPAR5,IFPAR6,IFPAR7,IFPAR8,
     *                IFPA09,IFPA10,IFPA11,IFPA12
      COMMON
     *       /PRFILE/ ISCREN,LOGFIL,LOGCHP,LOGCHN,LOGHIP,LOGHIN,
     *                LOGPRO,LOGNEU,LOGWRI,LSCREN,IRESUL,ICONVE,
     *                NOUTPT,LOGAUX,LOGENE,LOGRAD,LOGGRD,LOGMSH,
     *                                                   LOGBIS    
      COMMON
     *       /FITTED/ IFITED
      COMMON
     *       /VPROEF/ V0EFFC
      COMMON
     *       /LAMUNI/ UNITLA
      COMMON
     *       /FACCOU/ COUFAC
      COMMON
     *       /VERSIN/ VERSIO
      COMMON
     *       /CTSHIF/ V0SHIF_PROTON(1:NDNUCL),
     *                V0SHIF_NEUTRS(1:NDNUCL)
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
     *       /PARAMT_TITLES/ TITLES(1:NDTITL),
     *                       TITLES_LATEXS(1:NDTITL),
     *                       NUCNAM_LATEXS(1:NDNUCL)
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
      IF (LOGWRI.GT.0) THEN
          WRITE(LOGFIL,'(/,''Entering WRITING_ENELEV'')')
      
          WRITE(LOGFIL,'(/,''IFITED='',I2,'' IZ= '',I3,'' IN= '',I3,
     *                 '' ISOSPI= '',I1,/)') IFITED,NUMB_Z(INUCLI),
     *                                      NUMB_N(INUCLI),ISOSPI
      END IF
C
C=======================================================================
C
      IACTIV=0
C
C=======================================================================
C      
      IF (IFDENS.EQ.0 .AND. IFPROT.EQ.1) THEN
          
          IF (LOGWRI.GT.0) 
     *    WRITE(LOGFIL,'(''Assigning thing for pure WS and protons'',$
     *                                                             )')
          
          ISOSPI=1
          FILNAM_ENDING='-P'
C_______________________________________________________________________
C          
          LEVTHE=LDSING_PROTON 
          LEVEXP=LEVEXP_PROTON(INUCLI)
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
C_______________________________________________________________________
C          
          IPARAM_AUXILI=0
          
          DO IPARAM=1,6
            IPARAM_AUXILI=IPARAM_AUXILI+1
            TITPAR(IPARAM_AUXILI)=TITLES(IPARAM)(2:8)//'PROTON'
            PARPOT_AUXILI(IPARAM_AUXILI)=PARPOT_PRINTG(IPARAM)
          END DO
          
          NINITI=1
          NFINAL=IPARAM_AUXILI
          ITOTAL=IPARAM_AUXILI
          
          NINCNT=1
          NFICNT=NFINAL/2
          
          NIN_SO=NFICNT+1
          NFI_SO=NFINAL
          
          CHISQU=CHISQU_PROTON
          CHICOR=CHICOR_PROTON
          RADDIF=RADDIF_PROTON
          CHIINV=CHIINV_PROTON
          FERDIF=FERDIF_PROTON
          GAPDIF=GAPDIF_PROTON
          CHIWEI=CHIWEI_PROTON
          EABSAV=EABSAV_PROTON
          ERRMAX=ERRMAX_PROTON
          DIFFUP=DIFFUP_PROTON
          DIFFDW=DIFFDW_PROTON
          
          HINORM=HINORM_PROTON
          
          RMSERR=RMSERR_PROTON(INUCLI)
          
          IF (LOGWRI.GT.0) WRITE(LOGFIL,'('' ... OK!'')')
          
      END IF
C
C=======================================================================
C      
      IF (IFDENS.EQ.0 .AND. IFNEUT.EQ.1) THEN

          IF (LOGWRI.GT.0) 
     *    WRITE(LOGFIL,'(''Assigning thing for pure WS and neutrons'',$
     *                                                             )')
          
          ISOSPI=0
          FILNAM_ENDING='-N'
C_______________________________________________________________________
C          
          LEVTHE=LDSING_NEUTRS 
          LEVEXP=LEVEXP_NEUTRS(INUCLI)
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
C_______________________________________________________________________
C          
          IPARAM_AUXILI=0
          
          DO IPARAM=21,26
            IPARAM_AUXILI=IPARAM_AUXILI+1
            TITPAR(IPARAM_AUXILI)=TITLES(IPARAM)(2:8)//'NEUTRS'
            PARPOT_AUXILI(IPARAM_AUXILI)=PARPOT_PRINTG(IPARAM)
          END DO
          
          NINITI=1
          NFINAL=IPARAM_AUXILI
          ITOTAL=IPARAM_AUXILI
          
          NINCNT=1
          NFICNT=NFINAL/2
          
          NIN_SO=NFICNT+1
          NFI_SO=NFINAL
          
          CHISQU=CHISQU_NEUTRS
          CHICOR=CHICOR_NEUTRS
          RADDIF=RADDIF_NEUTRS
          CHIINV=CHIINV_NEUTRS
          FERDIF=FERDIF_NEUTRS
          GAPDIF=GAPDIF_NEUTRS
          CHIWEI=CHIWEI_NEUTRS
          EABSAV=EABSAV_NEUTRS
          ERRMAX=ERRMAX_NEUTRS
          DIFFUP=DIFFUP_NEUTRS
          DIFFDW=DIFFDW_NEUTRS
          
          HINORM=HINORM_NEUTRS
          
          RMSERR=RMSERR_NEUTRS(INUCLI)
          
          IF (LOGWRI.GT.0) WRITE(LOGFIL,'('' ... OK!'')')
          
      END IF
C
C=======================================================================
C
      IF (IFDENS.EQ.0 .AND. IFBOTH.EQ.1) THEN
C
          IF (ISOSPI.EQ.1) THEN
C
              IF (LOGWRI.GT.0) 
     *        WRITE(LOGFIL,'(''Assigning thing for pure WS '',
     *                                           ''and protons'',$)')
          
              FILNAM_ENDING='-P'
C_______________________________________________________________________
C          
              LEVTHE=LDSING_PROTON 
              LEVEXP=LEVEXP_PROTON(INUCLI)
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
C_______________________________________________________________________
C          
              IPARAM_AUXILI=0
          
              DO IPARAM=1,6
                 IPARAM_AUXILI=IPARAM_AUXILI+1
                 TITPAR(IPARAM_AUXILI)=TITLES(IPARAM)(2:8)//'PROTON'
                 PARPOT_AUXILI(IPARAM_AUXILI)=PARPOT_PRINTG(IPARAM)
              END DO
          
              NINITI=1
              NFINAL=IPARAM_AUXILI
              ITOTAL=IPARAM_AUXILI
          
              NINCNT=1
              NFICNT=NFINAL/2
          
              NIN_SO=NFICNT+1
              NFI_SO=NFINAL
          
              CHISQU=CHISQU_PROTON
              CHICOR=CHICOR_PROTON
              RADDIF=RADDIF_PROTON
              CHIINV=CHIINV_PROTON
              FERDIF=FERDIF_PROTON
              GAPDIF=GAPDIF_PROTON
              CHIWEI=CHIWEI_PROTON
              EABSAV=EABSAV_PROTON
              ERRMAX=ERRMAX_PROTON
              DIFFUP=DIFFUP_PROTON
              DIFFDW=DIFFDW_PROTON
          
              HINORM=HINORM_PROTON
          
              RMSERR=RMSERR_PROTON(INUCLI)
          
              IF (LOGWRI.GT.0) WRITE(LOGFIL,'('' ... OK!'')')
C
          END IF
C_______________________________________________________________________
C 
          IF (ISOSPI.EQ.0) THEN
              IF (LOGWRI.GT.0) 
     *            WRITE(LOGFIL,'(''Assigning thing for pure WS '',
     *                                           ''and neutrons'',$)')
              ISOSPI=0
              FILNAM_ENDING='-N'
C_______________________________________________________________________
C          
              LEVTHE=LDSING_NEUTRS 
              LEVEXP=LEVEXP_NEUTRS(INUCLI)
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
C_______________________________________________________________________
C          
              IPARAM_AUXILI=0
          
              DO IPARAM=21,26
                 IPARAM_AUXILI=IPARAM_AUXILI+1
                 TITPAR(IPARAM_AUXILI)=TITLES(IPARAM)(2:8)//'NEUTRS'
                 PARPOT_AUXILI(IPARAM_AUXILI)=PARPOT_PRINTG(IPARAM)
              END DO
          
              NINITI=1
              NFINAL=IPARAM_AUXILI
              ITOTAL=IPARAM_AUXILI
          
              NINCNT=1
              NFICNT=NFINAL/2
          
              NIN_SO=NFICNT+1
              NFI_SO=NFINAL
          
              CHISQU=CHISQU_NEUTRS
              CHICOR=CHICOR_NEUTRS
              RADDIF=RADDIF_NEUTRS
              CHIINV=CHIINV_NEUTRS
              FERDIF=FERDIF_NEUTRS
              GAPDIF=GAPDIF_NEUTRS
              CHIWEI=CHIWEI_NEUTRS
              EABSAV=EABSAV_NEUTRS
              ERRMAX=ERRMAX_NEUTRS
              DIFFUP=DIFFUP_NEUTRS
              DIFFDW=DIFFDW_NEUTRS
          
              HINORM=HINORM_NEUTRS
          
              RMSERR=RMSERR_NEUTRS(INUCLI)
          
              IF (LOGWRI.GT.0) WRITE(LOGFIL,'('' ... OK!'')')
          END IF
C
      END IF
C
C=======================================================================
C      
      IF (IFDENS.EQ.1) THEN
          
          IF (ISOSPI.EQ.1) THEN
              
              FILNAM_ENDING='-P'          
              
              LEVEXP=LEVEXP_PROTON(INUCLI)
C          
              DO IEXPER=1,LEVEXP
                 LABEXP(IEXPER)=LABEXP_PROTON(INUCLI,IEXPER)
                 EXPEXP(IEXPER)=EXPEXP_PROTON(INUCLI,IEXPER)
              END DO
              
              RMSERR=RMSERR_PROTON(INUCLI)
          
          END IF
C_______________________________________________________________________
C          
          IF (ISOSPI.EQ.0) THEN
              
              FILNAM_ENDING='-N'          
              
              LEVEXP=LEVEXP_NEUTRS(INUCLI)
C          
              DO IEXPER=1,LEVEXP
                 LABEXP(IEXPER)=LABEXP_NEUTRS(INUCLI,IEXPER)
                 EXPEXP(IEXPER)=EXPEXP_NEUTRS(INUCLI,IEXPER)
              END DO
              
              RMSERR=RMSERR_NEUTRS(INUCLI)
          
          END IF
C_______________________________________________________________________
C          
          IPARAM_AUXILI=0
C          
          DO IPARAM=1,6
             IPARAM_AUXILI=IPARAM_AUXILI+1
             TITPAR(IPARAM_AUXILI)=TITLES(IPARAM)(2:8)//'PROTON'
             PARPOT_AUXILI(IPARAM_AUXILI)=PARPOT_PRINTG(IPARAM)
             IF (IPARAM.GE.4) PARPOT_AUXILI(IPARAM_AUXILI)=9999.
          END DO 
C          
          NINCNT_PROTON=1
          NFICNT_PROTON=IPARAM_AUXILI/2
          
          NIN_WS_PROTON=IPARAM_AUXILI/2+1
          NFI_WS_PROTON=IPARAM_AUXILI
C          
          DO IPARAM=21,26
             IPARAM_AUXILI=IPARAM_AUXILI+1
             TITPAR(IPARAM_AUXILI)=TITLES(IPARAM)(2:8)//'NEUTRS'
             PARPOT_AUXILI(IPARAM_AUXILI)=PARPOT_PRINTG(IPARAM)
             IF (IPARAM.GE.24) PARPOT_AUXILI(IPARAM_AUXILI)=9999.
          END DO 
C          
          NINCNT_NEUTRS=NFI_WS_PROTON+1
          NFICNT_NEUTRS=NINCNT_NEUTRS+2
          
          NIN_WS_NEUTRS=NFICNT_NEUTRS+1
          NFI_WS_NEUTRS=IPARAM_AUXILI
C          
          DO IPARAM=39,42
             IPARAM_AUXILI=IPARAM_AUXILI+1
             TITPAR(IPARAM_AUXILI)=TITLES(IPARAM)(2:9)
             PARPOT_AUXILI(IPARAM_AUXILI)=PARPOT_PRINTG(IPARAM)/UNITLA
          END DO
C          
          NIN_SO_LAMBDA=NFI_WS_NEUTRS+1
          NFI_SO_LAMBDA=IPARAM_AUXILI
C          
          IF (IFTENS.EQ.1.AND.ISORBT.EQ.1) THEN
              DO IPARAM=43,46
                 IPARAM_AUXILI=IPARAM_AUXILI+1
                 TITPAR(IPARAM_AUXILI)=TITLES(IPARAM)(2:9)
                 PARPOT_AUXILI(IPARAM_AUXILI)=PARPOT_PRINTG(IPARAM)
     *                                       /UNITLA
              END DO
              NIN_SO_TENLAM=IPARAM_AUXILI-3
              NFI_SO_TENLAM=IPARAM_AUXILI
          END IF
C          
          IF (IFTENS.EQ.1.AND.ICENTT.EQ.1) THEN
              DO IPARAM=47,50
                 IPARAM_AUXILI=IPARAM_AUXILI+1
                 TITPAR(IPARAM_AUXILI)=TITLES(IPARAM)(2:9)
                 PARPOT_AUXILI(IPARAM_AUXILI)=PARPOT_PRINTG(IPARAM)
     *                                       /UNITLA
              END DO
              NINCNT_TENLAM=IPARAM_AUXILI-3
              NFICNT_TENLAM=IPARAM_AUXILI
          END IF
C          
          NINITI=1
          NFINAL=IPARAM_AUXILI
          ITOTAL=IPARAM_AUXILI
C          
          WRITE(TEXLAM,'(''_LAMBDA'',SP,I2,I2,I2,I2)')
     *                            IFPAR1,IFPAR2,IFPAR3,IFPAR4
C          
      END IF
C
C=======================================================================
C
      YMINIM=1.0E+10
      YMAXIM=-1.0E+10
      
      DO IEXPER=1,LEVEXP
        DO ITHEOR=1,LEVTHE_PRINTG
                
          IF (LABEXP(IEXPER).EQ.LABTHE_PRINTG(ITHEOR)) THEN

              IF (ENETHE_PRINTG(ITHEOR).LT.YMINIM) THEN
              
                  YMINIM=ENETHE_PRINTG(ITHEOR)
              
              END IF
          
              IF (ENETHE_PRINTG(ITHEOR).GT.YMAXIM) THEN
               
                  YMAXIM=ENETHE_PRINTG(ITHEOR)
              
              END IF
              
              IF (EXPEXP(IEXPER).LT.YMINIM) THEN
              
                  YMINIM=EXPEXP(IEXPER)
              
              END IF
          
              IF (EXPEXP(IEXPER).GT.YMAXIM) THEN
              
                  YMAXIM=EXPEXP(IEXPER)
              
              END IF
              
          END IF
                
        END DO
      END DO
C
C=======================================================================
C      
      NRESUL=55+INUCLI
C_______________________________________________________________________
C  
      IF (IFITED.EQ.3) FILNAM_IFITED='WSUniv/WSUniv-' !from WSuniv 
      IF (IFITED.EQ.2) FILNAM_IFITED='OneRun/OneRun-' !from single run
      IF (IFITED.EQ.1) FILNAM_IFITED='Fitted/Fitted-'
      IF (IFITED.EQ.0) FILNAM_IFITED='Predic/Predic-'
C_______________________________________________________________________
C
      I2=1
      DO II=1,6
         IF (NUCSYM(INUCLI)(II:II).EQ.' ') THEN
             I2=II+1
         END IF
      END DO
C_______________________________________________________________________
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
C_______________________________________________________________________
C      
      IF (IFDENS.EQ.0) THEN
C      
          WRITE(FILNAM,'(''EnergyLevels/IFDENS-0/'',A14,A,A2,A,
     *                   ''_IFDENS-'',I1,''_IFTENS-'',I1,
     *                   ''_IF-PAI-'',I1,
     *                   ''_IF-RAD-'',I1,''_IF-ORD-'',I1,
     *                   ''_IF-RHO-'',I1,
     *                   ''_IFDEEP-'',I1,''_IFPRON-'',I1,
     *                   ''_IFK-VC-'',I1,''_IFK-VS-'',I1,
     *                   ''_IFK-RC-'',I1,''_IFK-RS-'',I1,
     *                   ''_IFK-AC-'',I1,''_IFK-AS-'',I1,
     *                   ''_LDRAND-'',I5.5,''_'',A3,''.dat'')') 
     *                      FILNAM_IFITED,NUCSYM(INUCLI)(I2:6),
     *                      FILNAM_ENDING,FILNAM_NUCFIT(1:I1),
     *                      IFDENS,IFTENS,IF_PAI,IF_RAD,IF_INV,
     *                      IF_RHO,IFDEEP,IFPRON,IFK_VC,IFK_VS,
     *                      IFK_RC,IFK_RS,IFK_AC,IFK_AS,LDRAND,
     *                                                  VERSIO
C      
      END IF
C_______________________________________________________________________
C      
      IF (IFDENS.EQ.1) THEN
C      
          WRITE(FILNAM,'(''EnergyLevels/IFDENS-1_IFTENS-'',I1,''/'',
     *                                       A14,A,A2,A,
     *                   ''_IFDENS-'',I1,''_IFTENS-'',I1,
     *                   ''_IF-PAI-'',I1,
     *                   ''_IF-RAD-'',I1,''_IF-ORD-'',I1,
     *                   ''_IFDEEP-'',I1,''_IFPRON-'',I1,
     *                   ''_IFK-VC-'',I1,''_IFK-VS-'',I1,
     *                   ''_IFK-RC-'',I1,''_IFK-RS-'',I1,
     *                   ''_IFK-AC-'',I1,''_IFK-AS-'',I1,
     *                   ''_LDRAND-'',I5.5,A15,''_'',A3,''.dat'')') 
     *                              IFTENS,
     *                              FILNAM_IFITED,NUCSYM(INUCLI)(I2:6),
     *                              FILNAM_ENDING,FILNAM_NUCFIT(1:I1),
     *                              IFDENS,IFTENS,IF_PAI,IF_RAD,IF_INV,
     *                              IFDEEP,IFPRON,IFK_VC,IFK_VS,
     *                              IFK_RC,IFK_RS,IFK_AC,IFK_AS,
     *                              LDRAND,TEXLAM,VERSIO
C
          IF (IFTENS.EQ.1) THEN
C              
              WRITE(FILNAM,'(''EnergyLevels/IFDENS-1_IFTENS-'',I1,''/'',
     *                                       A14,A,A2,A,
     *                   ''_IFDENS-'',I1,''_IFTENS-'',I1,
     *                   ''_IF-PAI-'',I1,
     *                   ''_ICENTT-'',I1,''_ISORBT-'',I1,
     *                                   ''_ITENSR-'',I1,
     *                   ''_IF-RAD-'',I1,''_IF-ORD-'',I1,
     *                   ''_IFDEEP-'',I1,''_IFPRON-'',I1,
     *                   ''_IFK-VC-'',I1,''_IFK-VS-'',I1,
     *                   ''_IFK-RC-'',I1,''_IFK-RS-'',I1,
     *                   ''_IFK-AC-'',I1,''_IFK-AS-'',I1,
     *                   ''_LDRAND-'',I5.5,A15,''_'',A3,''.dat'')') 
     *                              IFTENS,
     *                              FILNAM_IFITED,NUCSYM(INUCLI)(I2:6),
     *                              FILNAM_ENDING,FILNAM_NUCFIT(1:I1),
     *                              IFDENS,IFTENS,IF_PAI,
     *                              ICENTT,ISORBT,ITENSR,
     *                              IF_RAD,IF_INV,
     *                              IFDEEP,IFPRON,IFK_VC,IFK_VS,
     *                              IFK_RC,IFK_RS,IFK_AC,IFK_AS,
     *                              LDRAND,TEXLAM,VERSIO
C
          END IF
C      
      END IF
C
C=======================================================================
C        
      OPEN(NRESUL,FILE=FILNAM,STATUS='UNKNOWN',FORM='FORMATTED')
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
          
      WRITE(NRESUL,'(15X,''Comparison Theory-Experiment'',/)')
C
C=======================================================================
C     
      WRITE(NRESUL,'(''<<NUCLIDINFO>>'',1X,''NUMB_Z'',3X,''NUMB_N'')')
C
      WRITE(NRESUL,'(15X,I3,6X,I3,/)')NUMB_Z(INUCLI),NUMB_N(INUCLI)
C_______________________________________________________________________
C
      A_MASS=NUMB_Z(INUCLI)+NUMB_N(INUCLI)
C
C=======================================================================
C
      WRITE(NRESUL,'(''<<NoNUCACTIV>>'',1X,''NUCACT'')')
C      
      WRITE(NRESUL,'(15X,I3,/)') NUCACT
C
C=======================================================================
C
      WRITE(NRESUL,'(''<<NUCLEI_FIT>> '',$)')
      
      JACTIV=0
      DO JNUCLI=1,LDNUCL
         IF (ITAKNU(JNUCLI).EQ.1) THEN
             JACTIV=JACTIV+1
             WRITE(NRESUL,'(''Z_No'',i2.2,3X,''N_No'',i2.2,''   '',$)')
     *                        JACTIV,JACTIV
         END IF
      END DO
      WRITE(NRESUL,'()')
      
      WRITE(NRESUL,'(14X,''  '',$)')
      DO JNUCLI=1,LDNUCL
         IF (ITAKNU(JNUCLI).EQ.1) THEN
             WRITE(NRESUL,'(I4,5X,I4,2X,''   '',$)')NUMB_Z(JNUCLI),
     *                                              NUMB_N(JNUCLI)
         END IF
      END DO
      WRITE(NRESUL,'()')
      WRITE(NRESUL,'()')
C
C=======================================================================
C
      WRITE(NRESUL,'(''<<EXPER_FILE>>'',1X,''IFDEEP'',3X
     *                                     ''IFPRON'')')
      
      WRITE(NRESUL,'(15X,2(I3,6X),/)') IFDEEP,IFPRON
C
C=======================================================================
C
      WRITE(NRESUL,'(''<<POTEN_INFO>>'',1X,''ISOSPI'',3X
     *                                     ''IFDENS'',3X,''IFTENS'',
     *                                                3X,''IF_PAI'')')
      
      WRITE(NRESUL,'(15X,4(I3,6X),/)') ISOSPI,IFDENS,IFTENS,IF_PAI
C
C=======================================================================
C
      WRITE(NRESUL,'(''<<TENSR_INFO>>'',1X,''ICENTT'',3X
     *                                     ''ISORBT'',3X,''ITENSR'')')
      
      WRITE(NRESUL,'(15X,3(I3,6X),/)') ICENTT,ISORBT,ITENSR
C
C=======================================================================
C
      WRITE(NRESUL,'(''<<TAKEN_CHI2>> IF_SPE   IF_RAD   IF_GAP   '',
     *               ''IF_FER   IF_DEN   IF_RHO   IF_INV'')')
      WRITE(NRESUL,'(15X,7(I3,6X),/)')IF_SPE,IF_RAD,IF_GAP,IF_FER,
     *                                       IF_DEN,IF_RHO,IF_INV
C
C=======================================================================
C
      WRITE(NRESUL,'(''<<IFKAPPAPAR>> IFK_VC   IFK_VS   IFK_RC   '',
     *                              ''IFK_RS   IFK_AC   IFK_AS'')')
      WRITE(NRESUL,'(15X,6(I3,6X),/)') IFK_VC,IFK_VS,IFK_RC,IFK_RS,
     *                                               IFK_AC,IFK_AS
C
C=======================================================================
C
      WRITE(NRESUL,'(''<<RANDOMNESS>> LDRAND   IRANDO'')')
C
      WRITE(NRESUL,'(15X,I6,3X,I6,/)') LDRAND,IRANDO_PRINTG
C
C=======================================================================
C          
      IF (IFDENS.EQ.0 .AND. IFPROT.EQ.1) THEN
C
          WRITE(NRESUL,'(''<<CENTRL_PRO>>'',1X,
     *                       ''VCENTR_PROTON  RCENTR_PROTON  '',
     *                       ''ACENTR_PROTON'')')
          WRITE(NRESUL,'(15X,3(F13.9,2X),/)') 
     *                   VKACEN_PROTON,RKACEN_PROTON/A_MASS**(1./3.),
     *                   AKACEN_PROTON
     
          WRITE(NRESUL,'(''<<SOR-WS_PRO>>'',1X,
     *                       ''VSORBT_PROTON  RSORBT_PROTON  '',
     *                       ''ASORBT_PROTON'')')
          WRITE(NRESUL,'(15X,3(F13.9,2X),/)') 
     *                   VKASOR_PROTON,RKASOR_PROTON/A_MASS**(1./3.),
     *                   AKASOR_PROTON
C
      END IF
C_______________________________________________________________________
C      
      IF (IFDENS.EQ.0 .AND. IFNEUT.EQ.1) THEN
C
          WRITE(NRESUL,'(''<<CENTRL_NEU>>'',1X,
     *                       ''VCENTR_NEUTRS  RCENTR_NEUTRS  '',
     *                       ''ACENTR_NEUTRS'')')
          WRITE(NRESUL,'(15X,3(F13.9,2X),/)') 
     *                   VKACEN_NEUTRS,RKACEN_NEUTRS/A_MASS**(1./3.),
     *                   AKACEN_NEUTRS
     
          WRITE(NRESUL,'(''<<SOR-WS_NEU>>'',1X,
     *                       ''VSORBT_NEUTRS  RSORBT_NEUTRS  '',
     *                       ''ASORBT_NEUTRS'')')
          WRITE(NRESUL,'(15X,3(F13.9,2X),/)') 
     *                   VKASOR_NEUTRS,RKASOR_NEUTRS/A_MASS**(1./3.),
     *                   AKASOR_NEUTRS
C
      END IF
C_______________________________________________________________________
C
      IF (IFDENS.EQ.0 .AND. IFBOTH.EQ.1) THEN
C
          IF (ISOSPI.EQ.1) THEN

              WRITE(NRESUL,'(''<<CENTRL_PRO>>'',1X,
     *                       ''VCENTR_PROTON  RCENTR_PROTON  '',
     *                       ''ACENTR_PROTON'')')
              WRITE(NRESUL,'(15X,3(F13.9,2X),/)') 
     *                   VKACEN_PROTON,RKACEN_PROTON/A_MASS**(1./3.),
     *                   AKACEN_PROTON
     
              WRITE(NRESUL,'(''<<SOR-WS_PRO>>'',1X,
     *                       ''VSORBT_PROTON  RSORBT_PROTON  '',
     *                       ''ASORBT_PROTON'')')
              WRITE(NRESUL,'(15X,3(F13.9,2X),/)') 
     *                   VKASOR_PROTON,RKASOR_PROTON/A_MASS**(1./3.),
     *                   AKASOR_PROTON
          END IF
C
          IF (ISOSPI.EQ.0) THEN

              WRITE(NRESUL,'(''<<CENTRL_NEU>>'',1X,
     *                       ''VCENTR_NEUTRS  RCENTR_NEUTRS  '',
     *                       ''ACENTR_NEUTRS'')')
              WRITE(NRESUL,'(15X,3(F13.9,2X),/)') 
     *                   VKACEN_NEUTRS,RKACEN_NEUTRS/A_MASS**(1./3.),
     *                   AKACEN_NEUTRS
     
              WRITE(NRESUL,'(''<<SOR-WS_NEU>>'',1X,
     *                       ''VSORBT_NEUTRS  RSORBT_NEUTRS  '',
     *                       ''ASORBT_NEUTRS'')')
              WRITE(NRESUL,'(15X,3(F13.9,2X),/)') 
     *                   VKASOR_NEUTRS,RKASOR_NEUTRS/A_MASS**(1./3.),
     *                   AKASOR_NEUTRS
          END IF
C
      END IF
C_______________________________________________________________________
C     
      IF (IFDENS.EQ.1) THEN
C
          WRITE(NRESUL,'(''<<CENTRL_PRO>>'',1X,
     *                       ''VCENTR_PROTON  RCENTR_PROTON  '',
     *                       ''ACENTR_PROTON'')')
          WRITE(NRESUL,'(15X,3(F13.9,2X),/)') 
     *                 VKACEN_PROTON,RKACEN_PROTON/A_MASS**(1./3.),
     *                 AKACEN_PROTON
C
          WRITE(NRESUL,'(''<<CENTRL_NEU>>'',1X,
     *                       ''VCENTR_NEUTRS  RCENTR_NEUTRS  '',
     *                       ''ACENTR_NEUTRS'')')
          WRITE(NRESUL,'(15X,3(F13.9,2X),/)') 
     *                 VKACEN_NEUTRS,RKACEN_NEUTRS/A_MASS**(1./3.),
     *                 AKACEN_NEUTRS
C          
          WRITE(NRESUL,'(''<<SOR_DENSIT>>'',1X,4(A13,2X))')
     *                 (TITPAR(I),I=NIN_SO_LAMBDA,NFI_SO_LAMBDA)
          WRITE(NRESUL,'(15X,4(F13.6,2X),/)') 
     *                 (PARPOT_AUXILI(I),I=NIN_SO_LAMBDA,NFI_SO_LAMBDA)
C          
          IF (IFTENS.EQ.1 .AND. ISORBT.EQ.1) THEN
              WRITE(NRESUL,'(''<<SOR_TENSOR>>'',1X,4(A13,2X))')
     *                 (TITPAR(I),I=NIN_SO_TENLAM,NFI_SO_TENLAM)
              WRITE(NRESUL,'(15X,4(F13.6,2X),/)') 
     *                 (PARPOT_AUXILI(I),I=NIN_SO_TENLAM,NFI_SO_TENLAM)
          END IF
C          
          IF (IFTENS.EQ.1 .AND. ICENTT.EQ.1) THEN
              WRITE(NRESUL,'(''<<CNT_TENSOR>>'',1X,4(A13,2X))')
     *                 (TITPAR(I),I=NINCNT_TENLAM,NFICNT_TENLAM)
              WRITE(NRESUL,'(15X,4(F13.6,2X),/)') 
     *                 (PARPOT_AUXILI(I),I=NINCNT_TENLAM,NFICNT_TENLAM)
          END IF
C
      END IF
C
C=======================================================================
C
      WRITE(NRESUL,'(''<<KAPPA_CENT>> V0CENT  KAP_V0  '',
     *              ''R0CENT  KAP_R0  A0CENT  KAP_A0'')')
C
      WRITE(NRESUL,'(15X,6(F13.6,2X),/)')V0CENT_KAPPAR,XK_V0C_KAPPAR,
     *                                  R0CENT_KAPPAR,XK_R0C_KAPPAR,
     *                                  A0CENT_KAPPAR,XK_A0C_KAPPAR
C
C=======================================================================
C
      WRITE(NRESUL,'(''<<KAPPA_SORB>> V0SORB  KAP_V0  '',
     *              ''R0SORB  KAP_R0  A0SORB  KAP_A0'')')
C
      WRITE(NRESUL,'(15X,6(F13.6,2X),/)')V0SORB_KAPPAR,XK_V0S_KAPPAR,
     *                                  R0SORB_KAPPAR,XK_R0S_KAPPAR,
     *                                  A0SORB_KAPPAR,XK_A0S_KAPPAR
C
C=======================================================================
C 
      IF (ISOSPI.EQ.1) THEN
C
          WRITE(NRESUL,'(''<<V0P_EFFECT>> V0EFFC'')')
          WRITE(NRESUL,'(15X,F8.4,/)') V0EFFC
C
          WRITE(NRESUL,'(''<<COULOMBFAC>> COUFAC'')')
          WRITE(NRESUL,'(15X,F8.4,/)') COUFAC
C
      END IF
C
C=======================================================================
C
      WRITE(NRESUL,'(''<<UNITLAMBDA>>'',1X,''UNITLA'')')
      WRITE(NRESUL,'(15X,F6.2,/)') UNITLA
C
C=======================================================================
C
      WRITE(NRESUL,'(''<<WEIGHT_NUC>>'',1X,<LDNUCL>(''NUC_'',I2.2,3X))')
     *                                               (I,I=1,LDNUCL)
      WRITE(NRESUL,'(15X,<LDNUCL>(F6.2,3X))') 
     *                                   (WEINUC_PROTON(I),I=1,LDNUCL)
      WRITE(NRESUL,'(15X,<LDNUCL>(F6.2,3X),/)') 
     *                                   (WEINUC_NEUTRS(I),I=1,LDNUCL)
C
C=======================================================================
C
      WRITE(NRESUL,'(''<<WEIGHT_RAD>>'',1X,<LDNUCL>(''NUC_'',I2.2,3X))')
     *                                               (I,I=1,LDNUCL)
      WRITE(NRESUL,'(15X,<LDNUCL>(F6.0,3X))') 
     *                                   (WEIRAD_PROTON(I),I=1,LDNUCL)
      WRITE(NRESUL,'(15X,<LDNUCL>(F6.0,3X),/)') 
     *                                   (WEIRAD_NEUTRS(I),I=1,LDNUCL)
C
C=======================================================================
C
      WRITE(NRESUL,'(''<<RADII_INFO>> RMSEXP  RMSTHE'',
     *                            ''  EXPERR  DIFFER'')')
      
      WRITE(NRESUL,'(15X,F6.4,2X,F6.4,2X,F6.4,2X,E11.4/)')
     *                                             RMSEXP_PRINTG,
     *                                             RMSTHE_PRINTG,
     *                                                    RMSERR,
     *                               RMSTHE_PRINTG-RMSEXP_PRINTG
C
C=======================================================================
C 
      WRITE(NRESUL,'(''<<ERROR_INFO>> MX_ABS  RMSVAL  RMSGLO'')')
      
      WRITE(NRESUL,'(15X,F6.3,2X,F6.3,2X,F6.3,/)')ERRMAX_PRINTG,
     *                                            RMSVAL_PRINTG,
     *                                            RMSGLO_PRINTG
C
C=======================================================================
C
      WRITE(NRESUL,'(''<<XAXIS_TEXT>>'',/,16X,
     *                   ''Spherical Woods-Saxon Hamiltonian'',/)')
          
      WRITE(NRESUL,'(''<<YAXIS_TEXT>>'',/,16X,
     *                   ''Single Particle Energies'',/)')
          
      WRITE(NRESUL,'(''<<SIDE_TEXTS>>'')')
          
      WRITE(NRESUL,'(16X,''Id.:~$N_{\\rm r-rest}='',I4,'',\\,'',
     *                   ''N_{\\rm f-eval}='',I5,'',\\,\\chi^2='',F10.4,
     *                       '',\\,\\nabla\\chi^2='',F10.4,''$'',/)')
     *                         LDRAND,IEVALU_PRINTG,
     *                                CHISQU_PRINTG,CHIGRD_PRINTG
      WRITE(NRESUL,'(16X,''~~~~~${\\rm IFDENS}= '',I1,'',\\,'',
     *                         ''{\\rm IF-RAD}= '',I1,/)')IFDENS,IF_RAD
C
C=======================================================================
C     
      WRITE(NRESUL,'(16X,30(''=''),/)')
          
      WRITE(NRESUL,'(''<<X_AX_PARAM>>  XMINIM_FIG  XMAXIM_FIG  '',
     *                   ''X_STEP_FIG'')')
          
      WRITE(NRESUL,'(16X,3(F10.2,2X),/)')
     *                         XMIN_T,XMAX_E,(XMAX_E-XMIN_E)
     
      WRITE(NRESUL,'(''<<Y_AX_PARAM>>  YMINIM_FIG  YMAXIM_FIG  '',
     *                   ''Y_STEP_FIG'')')
          
      WRITE(NRESUL,'(16X,3(F10.2,2X),/)')
     *                         YMINIM,YMAXIM,(YMAXIM-YMINIM)
     
      WRITE(NRESUL,'(16X,30(''=''),/)')
C
C=======================================================================
C
      NCURVE=LEVEXP*2
      NTESTG=0
      NUMB_POINT=2
          
      WRITE(NRESUL,'(''<<NUM_LEVELS>>'',2X,''LEVNUM'',2X,
     *                   ''NUMB_POINT'')')
          
      WRITE(NRESUL,'(16X,I6,2X,I6/)')LEVEXP,NUMB_POINT
          
      WRITE(NRESUL,'(''<<NO_OF_CURV>>  NUMB_CURVE  IF_BULLETS  '',
     *                   ''IF_HISTOGS  IF_OVERFITING'')')
          
      WRITE(NRESUL,'(16X,4(I6,6X))')NCURVE,NTESTG,NTESTG,NTESTG
C
C=======================================================================
C         
      ICOUNT_CURVES=0
      I_COLO=0
          
      DO IEXPER=1,LEVEXP
        DO ITHEOR=1,LEVTHE_PRINTG
C
          IF (LABEXP(IEXPER).EQ.LABTHE_PRINTG(ITHEOR)) THEN
     
              ICOUNT_CURVES=ICOUNT_CURVES+1
              
              AUXLAB=LABEXP(IEXPER)
                 
              CALL LATEXS_LABELS(AUXLAB,LABTEX)
     
              WRITE(NRESUL,'(''<<>>'')')
          
              WRITE(NRESUL,'(''<<CURVE_'',I4.4,''>>'',2X,
     *                       ''THICK_LINE  STYLE_LINE  '',
     *                       ''COLOR_LINE  TYPE_POINT'')')ICOUNT_CURVES
          
              WRITE(NRESUL,'(16X,4(I6,6X))')ITHICK(ICOUNT_CURVES),
     *                                      ISTYLE,I_COLO,I_TYPE
              
              WRITE(NRESUL,'(15X,F7.3,3X,F7.3,4X,A11)')
     *                  XMIN_T,ENETHE_PRINTG(ITHEOR),LABTEX
              WRITE(NRESUL,'(15X,F7.3,3X,F7.3,4X,A11)')
     *                  XMAX_T,ENETHE_PRINTG(ITHEOR),LABTEX
             
          END IF
                
        END DO
      END DO
         
      DO IEXPER=1,LEVEXP
        DO ITHEOR=1,LEVTHE_PRINTG
                
          IF (LABEXP(IEXPER).EQ.LABTHE_PRINTG(ITHEOR)) THEN
     
              ICOUNT_CURVES=ICOUNT_CURVES+1
              
              AUXLAB=LABEXP(IEXPER)
                 
              CALL LATEXS_LABELS(AUXLAB,LABTEX)
     
              WRITE(NRESUL,'(''<<>>'')')
          
              WRITE(NRESUL,'(''<<CURVE_'',I4.4,''>>'',2X,
     *                       ''THICK_LINE  STYLE_LINE  '',
     *                       ''COLOR_LINE  TYPE_POINT'')')ICOUNT_CURVES
          
              WRITE(NRESUL,'(16X,4(I6,6X))')ITHICK(ICOUNT_CURVES),
     *                                      ISTYLE,I_COLO,I_TYPE
              
              WRITE(NRESUL,'(15X,F7.3,3X,F7.3,4X,A11)')
     *                  XMIN_E,EXPEXP(IEXPER),LABTEX
              WRITE(NRESUL,'(15X,F7.3,3X,F7.3,4X,A11)')
     *                  XMAX_E,EXPEXP(IEXPER),LABTEX
             
          END IF
                
        END DO
      END DO
C
C=======================================================================
C          
      WRITE(NRESUL,'(''<<GO_GETTHEM>>'')')
C
C=======================================================================
C
      IF (LOGWRI.GT.0) WRITE(LOGFIL,'(''Exiting WRITING_ENELEV'',/)')
C
C=======================================================================
C      
      RETURN
      END 
C
C=======================================================================
C=======================================================================
C
      SUBROUTINE DENSIT_OUTPUT(CURNAM,N_CURV,ND_MAX,X_CURV,Y_CURV,
     *                         L_CURV,ISOSPI,IZ_NUM,IN_NUM,MAIN_T,
     *                         XTITLE,YTITLE,ISPLIT,INUCLI,LDRAND)
      
      INCLUDE 'MATDIM/NDMESH.f'
      INCLUDE 'MATDIM/NDSPEC.f'
      INCLUDE 'MATDIM/NDPARS.f'
      INCLUDE 'MATDIM/NDNUCL.f'
      INCLUDE 'MATDIM/NDCOLO.f'
      INCLUDE 'MATDIM/NDTITL.f'
      
      CHARACTER
     *          NUCSYM*6,INPSYM*6,TYPCHI*6,VERSIO*3
      CHARACTER
     *          FILNAM*256,L_CURV*100,NUCNAM*6,TITLES*12,FITNUC*2,
     *          ISONAM*002,AUXTIT*12,TITTEX*30,TEXLAM*20,TITPAR*13,
     *          CURNAM*020,MAIN_T*256,XTITLE*256,YTITLE*256
      CHARACTER
     *          FILNAM_IFITED*14,FILNAM_NUCFIT*30,
     *          FITNUC_AUXILI*2,TITLES_LATEXS*050,
     *                          NUCNAM_LATEXS*010
      DIMENSION
     *          ND_MAX(1:NDSPEC),
     *          X_CURV(1:NDMESH,1:NDSPEC),
     *          Y_CURV(1:NDMESH,1:NDSPEC),
     *          L_CURV(1:NDMESH,1:NDSPEC)
      DIMENSION 
     *          TITPAR(1:NDPARS),
     *          PARPOT_AUXILI(1:NDPARS)
      DIMENSION
     *          FITNUC(1:NDNUCL),
     *          FITNUC_AUXILI(1:NDNUCL)
      COMMON
     *       /DENSTR/ IFDENS
      COMMON
     *       /TENSOR/ IFTENS
      COMMON
     *       /TENOPT/ ICENTT,ISORBT,ITENSR 
      COMMON
     *       /P_N_LM/ IFPROT,IFNEUT,IFBOTH
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
     *       /DELTAF/ DFACTO(1:NDNUCL,1:2),
     *       /BCSDEL/ DELT_2(1:NDNUCL,1:2)
      COMMON
     *       /BCSLAM/ XLAMBD_PROTON(1:NDNUCL),
     *                XLAMBD_NEUTRS(1:NDNUCL)
      COMMON
     *       /PRFILE/ ISCREN,LOGFIL,LOGCHP,LOGCHN,LOGHIP,LOGHIN,
     *                LOGPRO,LOGNEU,LOGWRI,LSCREN,IRESUL,ICONVE,
     *                NOUTPT,LOGAUX,LOGENE,LOGRAD,LOGGRD,LOGMSH,
     *                                                   LOGBIS 
      COMMON
     *       /ATTRIB/ ISTYLE,I_TYPE,ICOLOR(1:NDCOLO),ITHICK(1:NDCOLO)
      COMMON
     *       /VPROEF/ V0EFFC
      COMMON
     *       /LAMUNI/ UNITLA
      COMMON
     *       /FITTED/ IFITED
      COMMON
     *       /CTSHIF/ V0SHIF_PROTON(1:NDNUCL),
     *                V0SHIF_NEUTRS(1:NDNUCL)
      COMMON
     *       /VERSIN/ VERSIO
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
      IF (LOGWRI.GT.0) WRITE(LOGFIL,'(/,''Entering DENSIT_OUTPUT'')')
      
      IF (LOGWRI.GT.0) 
     *    WRITE(LOGFIL,'(/,''IFITED='',I2,'' IZ= '',I3,'' IN= '',I3,/)') 
     *                  IFITED,NUMB_Z(INUCLI),NUMB_N(INUCLI)
C
C=======================================================================
C
      YMINIM=+1.0E+10
      YMAXIM=-1.0E+10
      
      DO I_CURV=1,N_CURV
         DO IPOINT=1,ND_MAX(I_CURV)
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
      
      DO I_CURV=1,N_CURV
         DO IPOINT=1,ND_MAX(I_CURV)
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
      IF (ISOSPI.EQ.2) ISONAM='-B'
C
C=======================================================================
C          
      IF (IFITED.EQ.2) FILNAM_IFITED='OneRun/OneRun-'  !from single run
      IF (IFITED.EQ.1) FILNAM_IFITED='Fitted/Fitted-'
      IF (IFITED.EQ.0) FILNAM_IFITED='Predic/Predic-'
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
C      
          WRITE(FILNAM,'(''DensitCurves/IFDENS-0/'',A14,A5,A,A2,A,
     *                   ''_IFDENS-'',I1,''_IFTENS-'',I1,
     *                   ''_IF-PAI-'',I1,
     *                   ''_IF-RAD-'',I1,''_IF-ORD-'',I1,
     *                   ''_IFDEEP-'',I1,''_IFPRON-'',I1,
     *                   ''_LDRAND-'',I5.5,''_'',A3,''.dat'')') 
     *                                    FILNAM_IFITED,CURNAM,
     *                                    NUCNAM(J1:J2),ISONAM,
     *                                    FILNAM_NUCFIT(01:I1),
     *                                    IFDENS,IFTENS,IF_PAI,
     *                                    IF_RAD,IF_INV,
     *                                    IFDEEP,IFPRON,LDRAND,
     *                                                  VERSIO
C     
      END IF
C
C=======================================================================
C      
      IF (IFDENS.EQ.1.AND.IFTENS.EQ.0) THEN 
C
          WRITE(TEXLAM,'(''_LAMBDA'',SP,I2,I2,I2,I2)')
     *                            IFPAR1,IFPAR2,IFPAR3,IFPAR4
C      
          WRITE(FILNAM,'(''DensitCurves/IFDENS-1_IFTENS-'',I1,''/'',
     *                                           A14,A5,A,A2,A,    
     *                         ''_IFDENS-'',I1,''_IFTENS-'',I1,
     *                         ''_IF-PAI-'',I1
     *                         ''_IF-RAD-'',I1,''_IF-ORD-'',I1,
     *                         ''_IFDEEP-'',I1,''_IFPRON-'',I1,
     *                         ''_LDRAND-'',I5.5,A15,''_'',A3,
     *                                               ''.dat'')') 
     *                      IFTENS,FILNAM_IFITED,CURNAM,NUCNAM(J1:J2),
     *                      ISONAM,FILNAM_NUCFIT(01:I1),
     *                      IFDENS,IFTENS,IF_PAI,IF_RAD,IF_INV,
     *                      IFDEEP,IFPRON,LDRAND,TEXLAM,VERSIO
C     
      END IF
C
C=======================================================================
C      
      IF (IFDENS.EQ.1.AND.IFTENS.EQ.1) THEN 
C          
          WRITE(TEXLAM,'(''_LAMBDA'',SP,I2,I2,I2,I2)')
     *                            IFPAR1,IFPAR2,IFPAR3,IFPAR4
C      
          WRITE(FILNAM,'(''DensitCurves/IFDENS-1_IFTENS-'',I1,''/'',
     *                                           A14,A5,A,A2,A,    
     *                         ''_IFDENS-'',I1,''_IFTENS-'',I1,
     *                         ''_IF-PAI-'',I1,
     *                         ''_ICENTT-'',I1,''_ISORBT-'',I1,
     *                                         ''_ITENSR-'',I1,
     *                         ''_IF-RAD-'',I1,''_IF-ORD-'',I1,
     *                         ''_IFDEEP-'',I1,''_IFPRON-'',I1,
     *                         ''_LDRAND-'',I5.5,A15,''_'',A3,
     *                                               ''.dat'')') 
     *                      IFTENS,FILNAM_IFITED,CURNAM,NUCNAM(J1:J2),
     *                      ISONAM,FILNAM_NUCFIT(01:I1),
     *                                    IFDENS,IFTENS,IF_PAI,
     *                                    ICENTT,ISORBT,ITENSR,
     *                             IF_RAD,IF_INV,IFDEEP,IFPRON,
     *                                           LDRAND,TEXLAM,
     *                                                  VERSIO
C     
      END IF
C
C=======================================================================
C
      OPEN(NRESUL,FILE=FILNAM,STATUS='UNKNOWN')
C
C=======================================================================
C      
      IF (IFDENS.EQ.0 .AND. IFPROT.EQ.1) THEN
C                    
          IPARAM_AUXILI=0
          
          DO IPARAM=1,6
             IPARAM_AUXILI=IPARAM_AUXILI+1
             TITPAR(IPARAM_AUXILI)=TITLES(IPARAM)(2:8)//'PROTON'
             PARPOT_AUXILI(IPARAM_AUXILI)=PARPOT(IPARAM)
          END DO
          
          NINITI=1
          NFINAL=IPARAM_AUXILI
          ITOTAL=IPARAM_AUXILI

          NINCNT=1
          NFICNT=NFINAL/2
          
          NIN_SO=NFICNT+1
          NFI_SO=NFINAL
          
      END IF
C
C=======================================================================
C      
      IF (IFDENS.EQ.0 .AND. IFNEUT.EQ.1) THEN
C                    
          IPARAM_AUXILI=0
          
          DO IPARAM=21,26
             IPARAM_AUXILI=IPARAM_AUXILI+1
             TITPAR(IPARAM_AUXILI)=TITLES(IPARAM)(2:8)//'NEUTRS'
             PARPOT_AUXILI(IPARAM_AUXILI)=PARPOT(IPARAM)
          END DO
          
          NINITI=1
          NFINAL=IPARAM_AUXILI
          ITOTAL=IPARAM_AUXILI
          
          NINCNT=1
          NFICNT=NFINAL/2
          
          NIN_SO=NFICNT+1
          NFI_SO=NFINAL
          
      END IF
C
C=======================================================================
C
      IF (IFDENS.EQ.0 .AND. ISOSPI.EQ.2) THEN
C
         IPARAM_AUXILI=0
C          
          DO IPARAM=1,6
             IPARAM_AUXILI=IPARAM_AUXILI+1
             TITPAR(IPARAM_AUXILI)=TITLES(IPARAM)(2:8)//'PROTON'
             PARPOT_AUXILI(IPARAM_AUXILI)=PARPOT(IPARAM)
          END DO 
C
          NINCNT_PROTON=1
          NFICNT_PROTON=IPARAM_AUXILI/2
          
          NIN_WS_PROTON=IPARAM_AUXILI/2+1
          NFI_WS_PROTON=IPARAM_AUXILI
C          
          DO IPARAM=21,26
             IPARAM_AUXILI=IPARAM_AUXILI+1
             TITPAR(IPARAM_AUXILI)=TITLES(IPARAM)(2:8)//'NEUTRS'
             PARPOT_AUXILI(IPARAM_AUXILI)=PARPOT(IPARAM)
          END DO 
C          
          NINCNT_NEUTRS=NFI_WS_PROTON+1
          NFICNT_NEUTRS=NINCNT_NEUTRS+2
          
          NIN_WS_NEUTRS=NFICNT_NEUTRS+1
          NFI_WS_NEUTRS=IPARAM_AUXILI
C          
          IFNEUT=0
          IFPROT=0
C
      END IF
C
C=======================================================================
C
      IF (IFDENS.EQ.1) THEN   
          
          IPARAM_AUXILI=0
C          
          DO IPARAM=1,6
             IPARAM_AUXILI=IPARAM_AUXILI+1
             TITPAR(IPARAM_AUXILI)=TITLES(IPARAM)(2:8)//'PROTON'
             PARPOT_AUXILI(IPARAM_AUXILI)=PARPOT(IPARAM)
             IF (IPARAM.GE.4) PARPOT_AUXILI(IPARAM_AUXILI)=9999.
          END DO 
C
          NINCNT_PROTON=1
          NFICNT_PROTON=IPARAM_AUXILI/2
          
          NIN_WS_PROTON=IPARAM_AUXILI/2+1
          NFI_WS_PROTON=IPARAM_AUXILI
C          
          DO IPARAM=21,26
             IPARAM_AUXILI=IPARAM_AUXILI+1
             TITPAR(IPARAM_AUXILI)=TITLES(IPARAM)(2:8)//'NEUTRS'
             PARPOT_AUXILI(IPARAM_AUXILI)=PARPOT(IPARAM)
             IF (IPARAM.GE.24) PARPOT_AUXILI(IPARAM_AUXILI)=9999.
          END DO 
C          
          NINCNT_NEUTRS=NFI_WS_PROTON+1
          NFICNT_NEUTRS=NINCNT_NEUTRS+2
          
          NIN_WS_NEUTRS=NFICNT_NEUTRS+1
          NFI_WS_NEUTRS=IPARAM_AUXILI
C         
          DO IPARAM=39,42
             IPARAM_AUXILI=IPARAM_AUXILI+1
             TITPAR(IPARAM_AUXILI)=TITLES(IPARAM)(2:9)
             PARPOT_AUXILI(IPARAM_AUXILI)=PARPOT(IPARAM)/UNITLA
          END DO
C          
          NIN_SO_LAMBDA=NFI_WS_NEUTRS+1
          NFI_SO_LAMBDA=IPARAM_AUXILI
C          
          IF (IFTENS.EQ.1.AND.ISORBT.EQ.1) THEN
              DO IPARAM=43,46
                 IPARAM_AUXILI=IPARAM_AUXILI+1
                 TITPAR(IPARAM_AUXILI)=TITLES(IPARAM)(2:9)
                 PARPOT_AUXILI(IPARAM_AUXILI)=PARPOT(IPARAM)/UNITLA
              END DO
              NIN_SO_TENLAM=IPARAM_AUXILI-3
              NFI_SO_TENLAM=IPARAM_AUXILI
          END IF
C          
          IF (IFTENS.EQ.1.AND.ICENTT.EQ.1) THEN
              DO IPARAM=47,50
                 IPARAM_AUXILI=IPARAM_AUXILI+1
                 TITPAR(IPARAM_AUXILI)=TITLES(IPARAM)(2:9)
                 PARPOT_AUXILI(IPARAM_AUXILI)=PARPOT(IPARAM)/UNITLA
              END DO
              NINCNT_TENLAM=IPARAM_AUXILI-3
              NFICNT_TENLAM=IPARAM_AUXILI
          END IF
C          
          NINITI=1
          NFINAL=IPARAM_AUXILI
          ITOTAL=IPARAM_AUXILI
C          
          WRITE(TEXLAM,'(''_LAMBDA'',SP,I2,I2,I2,I2)')
     *                            IFPAR1,IFPAR2,IFPAR3,IFPAR4
C          
      END IF   
C
C=======================================================================
C      
      DO I_CURV=1,N_CURV
         DO IPOINT=1,ND_MAX(I_CURV)
            DO K=1,100
               IF (L_CURV(IPOINT,I_CURV)(K:K).EQ.'/') THEN
                   L_CURV(IPOINT,I_CURV)(K:K)='@'
               END IF
            END DO
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
      WRITE(NRESUL,'(15X,A,/)') MAIN_T
C_______________________________________________________________________
C     
      WRITE(NRESUL,'(''<<NUCLIDINFO>>'',1X,''NUMB_Z'',2X,''NUMB_N'')')
C
      WRITE(NRESUL,'(15X,I3,6X,I3,/)') IZ_NUM,IN_NUM
C
      A_MASS=IZ_NUM+IN_NUM
C_______________________________________________________________________
C
      WRITE(NRESUL,'(''<<NoNUCACTIV>>'',1X,''NUCACT'')')
C      
      WRITE(NRESUL,'(15X,I3,/)') NUCACT
C_______________________________________________________________________
C
      WRITE(NRESUL,'(''<<EXPER_FILE>>'',1X,''IFDEEP'',2X
     *                                     ''IFPRON'')')
      
      WRITE(NRESUL,'(15X,2(I3,6X),/)') IFDEEP,IFPRON
C_______________________________________________________________________
C
      WRITE(NRESUL,'(''<<POTEN_INFO>>'',1X,''ISOSPI'',2X,''IFDENS'',
     *                                  2X,''IFTENS'',2X,''IF_PAI'')')
C      
      WRITE(NRESUL,'(15X,4(I3,5X),/)') ISOSPI,IFDENS,IFTENS,IF_PAI
C_______________________________________________________________________
C
      WRITE(NRESUL,'(''<<TENSR_INFO>>'',1X,''ICENTT'',2X
     *                                     ''ISORBT'',2X,''ITENSR'')')
C      
      WRITE(NRESUL,'(15X,3(I3,5X),/)') ICENTT,ISORBT,ITENSR
C_______________________________________________________________________
C
      WRITE(NRESUL,'(''<<TAKEN_CHI2>> IF_SPE  IF_RAD  IF_GAP  '',
     *               ''IF_FER  IF_DEN  IF_RHO  IF_INV'')')
      WRITE(NRESUL,'(15X,7(I3,5X),/)')IF_SPE,IF_RAD,IF_GAP,IF_FER,
     *                                       IF_DEN,IF_RHO,IF_INV
C_______________________________________________________________________
C          
      IF (IFDENS.EQ.0 .AND. IFPROT.EQ.1) THEN
          
          WRITE(NRESUL,'(''<<CENTRL_PRO>>'',1X,3(A13,2X))')
     *                          (TITPAR(I),I=NINCNT,NFICNT)
          WRITE(NRESUL,'(15X,3(F13.4,2X),/)') 
     *                   VKACEN_PROTON,RKACEN_PROTON/A_MASS**(1./3.),
     *                   AKACEN_PROTON
     
          WRITE(NRESUL,'(''<<SOR-WS_PRO>>'',1X,3(A13,2X))')
     *                          (TITPAR(I),I=NIN_SO,NFI_SO)
          WRITE(NRESUL,'(15X,3(F13.4,2X),/)') 
     *                   VKASOR_PROTON,RKASOR_PROTON/A_MASS**(1./3.),
     *                   AKASOR_PROTON
          
      END IF
C_______________________________________________________________________
C      
      IF (IFDENS.EQ.0 .AND. IFNEUT.EQ.1) THEN
          
          WRITE(NRESUL,'(''<<CENTRL_NEU>>'',1X,3(A13,4X))')
     *                          (TITPAR(I),I=NINCNT,NFICNT)
          WRITE(NRESUL,'(15X,3(F13.4,2X),/)') 
     *                   VKACEN_NEUTRS,RKACEN_NEUTRS/A_MASS**(1./3.),
     *                   AKACEN_NEUTRS
     
          WRITE(NRESUL,'(''<<SOR-WS_NEU>>'',1X,3(A13,4X))')
     *                          (TITPAR(I),I=NIN_SO,NFI_SO)
          WRITE(NRESUL,'(15X,3(F13.4,2X),/)') 
     *                   VKASOR_NEUTRS,RKASOR_NEUTRS/A_MASS**(1./3.),
     *                   AKASOR_NEUTRS
          
      END IF
C_______________________________________________________________________
C
      IF (IFDENS.EQ.0 .AND. ISOSPI.EQ.2) THEN
C
          WRITE(NRESUL,'(''<<CENTRL_PRO>>'',1X,3(A13,2X))')
     *                 (TITPAR(I),I=NINCNT_PROTON,NFICNT_PROTON)
          WRITE(NRESUL,'(15X,3(F13.4,2X),/)') 
     *                 VKACEN_PROTON,RKACEN_PROTON/A_MASS**(1./3.),
     *                 AKACEN_PROTON
C
          WRITE(NRESUL,'(''<<CENTRL_NEU>>'',1X,3(A13,2X))')
     *                 (TITPAR(I),I=NINCNT_NEUTRS,NFICNT_NEUTRS)
          WRITE(NRESUL,'(15X,3(F13.4,2X),/)') 
     *                   VKACEN_NEUTRS,RKACEN_NEUTRS/A_MASS**(1./3.),
     *                   AKACEN_NEUTRS
C
          WRITE(NRESUL,'(''<<SOR-WS_PRO>>'',1X,3(A13,2X))')
     *                  (TITPAR(I),I=NIN_WS_PROTON,NFI_WS_PROTON)
          WRITE(NRESUL,'(15X,3(F13.4,2X),/)')
     *                   VKASOR_PROTON,RKASOR_PROTON/A_MASS**(1./3.),
     *                   AKASOR_PROTON
C
          WRITE(NRESUL,'(''<<SOR-WS_NEU>>'',1X,3(A13,2X))')
     *                 (TITPAR(I),I=NIN_WS_NEUTRS,NFI_WS_NEUTRS)
          WRITE(NRESUL,'(15X,3(F13.4,2X),/)') 
     *                   VKASOR_NEUTRS,RKASOR_NEUTRS/A_MASS**(1./3.),
     *                   AKASOR_NEUTRS
C
      END IF
C_______________________________________________________________________
C      
      IF (IFDENS.EQ.1) THEN
C
          WRITE(NRESUL,'(''<<CENTRL_PRO>>'',1X,3(A13,2X))')
     *                 (TITPAR(I),I=NINCNT_PROTON,NFICNT_PROTON)
          WRITE(NRESUL,'(15X,3(F13.4,2X),/)') 
     *                 VKACEN_PROTON,RKACEN_PROTON/A_MASS**(1./3.),
     *                 AKACEN_PROTON
C
          WRITE(NRESUL,'(''<<CENTRL_NEU>>'',1X,3(A13,2X))')
     *                 (TITPAR(I),I=NINCNT_NEUTRS,NFICNT_NEUTRS)
          WRITE(NRESUL,'(15X,3(F13.4,2X),/)') 
     *                 VKASOR_PROTON,RKASOR_PROTON/A_MASS**(1./3.),
     *                 AKASOR_PROTON
C
          WRITE(NRESUL,'(''<<SOR-WS_PRO>>'',1X,3(A13,2X))')
     *                  (TITPAR(I),I=NIN_WS_PROTON,NFI_WS_PROTON)
          WRITE(NRESUL,'(15X,3(F13.4,2X),/)') 
     *                 VKACEN_NEUTRS,RKACEN_NEUTRS/A_MASS**(1./3.),
     *                 AKACEN_NEUTRS
C
          WRITE(NRESUL,'(''<<SOR-WS_NEU>>'',1X,3(A13,4X))')
     *                 (TITPAR(I),I=NIN_WS_NEUTRS,NFI_WS_NEUTRS)
          WRITE(NRESUL,'(15X,3(F13.4,2X),/)') 
     *                 VKASOR_NEUTRS,RKASOR_NEUTRS/A_MASS**(1./3.),
     *                 AKASOR_NEUTRS
          
          WRITE(NRESUL,'(''<<SOR_DENSIT>>'',1X,4(A13,4X))')
     *                 (TITPAR(I),I=NIN_SO_LAMBDA,NFI_SO_LAMBDA)
          WRITE(NRESUL,'(15X,4(F13.4,2X),/)') 
     *                 (PARPOT_AUXILI(I),I=NIN_SO_LAMBDA,NFI_SO_LAMBDA)
          
          IF (IFTENS.EQ.1 .AND. ISORBT.EQ.1) THEN
              WRITE(NRESUL,'(''<<SOR_TENSOR>>'',1X,4(A13,4X))')
     *                 (TITPAR(I),I=NIN_SO_TENLAM,NFI_SO_TENLAM)
              WRITE(NRESUL,'(15X,4(F13.4,2X),/)') 
     *                 (PARPOT_AUXILI(I),I=NIN_SO_TENLAM,NFI_SO_TENLAM)
          END IF
          
          IF (IFTENS.EQ.1 .AND. ICENTT.EQ.1) THEN
              WRITE(NRESUL,'(''<<CNT_TENSOR>>'',1X,4(A13,4X))')
     *                 (TITPAR(I),I=NINCNT_TENLAM,NFICNT_TENLAM)
              WRITE(NRESUL,'(15X,4(F13.4,2X),/)') 
     *                 (PARPOT_AUXILI(I),I=NINCNT_TENLAM,NFICNT_TENLAM)
          END IF
C
      END IF
C
C=======================================================================
C 
      IF (ISOSPI.EQ.1) THEN
          WRITE(NRESUL,'(''<<V0P_EFFECT>>  V0EFFC'')')
          WRITE(NRESUL,'(16X,F8.4,/)') V0EFFC
      END IF
C
C=======================================================================
C
      IF (IFITED.EQ.2) THEN
          WRITE(NRESUL,'(''<<SHIFTCENTR>>'',1X,''V0SHIF_PROTON'',
     *                                      2X,''V0SHIF_NEUTRS'')')
          WRITE(NRESUL,'(15X,F13.4,2X,F13.4,/)')V0SHIF_PROTON(INUCLI),
     *                                          V0SHIF_NEUTRS(INUCLI)
      END IF
C
C=======================================================================
C
      WRITE(NRESUL,'(''<<UNITLAMBDA>>'',2X,''UNITLA'')')
      WRITE(NRESUL,'(16X,F6.2,/)') UNITLA
C
C=======================================================================
C 
      IF (IF_PAI.EQ.1) THEN
C
          WRITE(NRESUL,'(''<<DELTAFACTO>>'',$)')
          DO JNUCLI=1,LDNUCL
             IF (ITAKNU(JNUCLI).EQ.1) THEN
                 WRITE(NRESUL,'(3X,A2,''   '',$)') FITNUC(JNUCLI)
             END IF
          END DO
          WRITE(NRESUL,'()')
          WRITE(NRESUL,'(14X,'' '',$)')
          DO JNUCLI=1,LDNUCL
             IF (ITAKNU(JNUCLI).EQ.1) THEN
                 WRITE(NRESUL,'(F6.4,''  '',$)') DFACTO(JNUCLI,1)
             END IF
          END DO
          WRITE(NRESUL,'()')
          WRITE(NRESUL,'(14X,'' '',$)')
          DO JNUCLI=1,LDNUCL
             IF (ITAKNU(JNUCLI).EQ.1) THEN
                 WRITE(NRESUL,'(F6.4,''  '',$)') DFACTO(JNUCLI,2)
             END IF
          END DO
C_______________________________________________________________________
C
          WRITE(NRESUL,'()')
          WRITE(NRESUL,'(/,''<<BCS_DELTA2>>'',$)')
          DO JNUCLI=1,LDNUCL
             IF (ITAKNU(JNUCLI).EQ.1) THEN
                 WRITE(NRESUL,'(3X,A2,''   '',$)') FITNUC(JNUCLI)
             END IF
          END DO
          WRITE(NRESUL,'()')
          WRITE(NRESUL,'(14X,'' '',$)')
          DO JNUCLI=1,LDNUCL
             IF (ITAKNU(JNUCLI).EQ.1) THEN
                 WRITE(NRESUL,'(F6.4,''  '',$)') DELT_2(JNUCLI,1)
             END IF
          END DO
          WRITE(NRESUL,'()')
          WRITE(NRESUL,'(14X,'' '',$)')
          DO JNUCLI=1,LDNUCL
             IF (ITAKNU(JNUCLI).EQ.1) THEN
                 WRITE(NRESUL,'(F6.4,''  '',$)') DELT_2(JNUCLI,2)
             END IF
          END DO
C_______________________________________________________________________
C
          WRITE(NRESUL,'()')
          WRITE(NRESUL,'(/,''<<BCS_XLAMBD>>'',$)')
          DO JNUCLI=1,LDNUCL
             IF (ITAKNU(JNUCLI).EQ.1) THEN
                 WRITE(NRESUL,'(3X,A2,''   '',$)') FITNUC(JNUCLI)
             END IF
          END DO
          WRITE(NRESUL,'()')
          WRITE(NRESUL,'(10X,'' '',$)')
          DO JNUCLI=1,LDNUCL
             IF (ITAKNU(JNUCLI).EQ.1) THEN
                 WRITE(NRESUL,'(F10.4,''  '',$)') XLAMBD_PROTON(JNUCLI)
             END IF
          END DO
          WRITE(NRESUL,'()')
          WRITE(NRESUL,'(10X,'' '',$)')
          DO JNUCLI=1,LDNUCL
             IF (ITAKNU(JNUCLI).EQ.1) THEN
                 WRITE(NRESUL,'(F10.4,''  '',$)') XLAMBD_NEUTRS(JNUCLI)
             END IF
          END DO
          WRITE(NRESUL,'()')
          WRITE(NRESUL,'()')
c          
      END IF
C
C=======================================================================
C      
      WRITE(NRESUL,'(''<<XAXIS_TEXT>>'',/,15x,A,/)') XTITLE
C      
      WRITE(NRESUL,'(''<<YAXIS_TEXT>>'',/,15x,A,/)') YTITLE
C      
      WRITE(NRESUL,'(''<<SIDE_TEXTS>>'')')
C      
      WRITE(NRESUL,'(15X,''~~~'',$)')
C
      IF (IFDENS.EQ.1) WRITE(NRESUL,'(''${\\rm DENS}\\,\\vert\\,$'',$)')
      IF (IFTENS.EQ.1) WRITE(NRESUL,'(''${\\rm TENS}\\,\\vert\\,$'',$)')
      IF (IF_PAI.EQ.1) WRITE(NRESUL,'(''${\\rm PAIR}\\,\\vert\\,$'',$)')
      IF (IF_RAD.EQ.1) WRITE(NRESUL,'(''${\\rm RADI}\\,\\vert\\,$'',$)')
      IF (IFPRON.EQ.1) WRITE(NRESUL,'(''${\\rm PRON}\\,\\vert\\,$'',$)')
      IF (IFDEEP.EQ.1) WRITE(NRESUL,'(''${\\rm DEEP}\\,\\vert\\,$'',$)')
C
      WRITE(NRESUL,'(''$\\Delta V_{0,p}^{cent}='',f9.2,''\\,\\vert\\,'',
     *                ''\\Delta V_{0,n}^{cent}='',f0.2,''$'')')
     *                      V0SHIF_PROTON(INUCLI),V0SHIF_NEUTRS(INUCLI)
C     
      IF (IFTENS.EQ.1) THEN
          WRITE(NRESUL,'(15X,''~~~~~~~~'',$)')
          IF (ISORBT.EQ.1) 
     *        WRITE(NRESUL,'(''${\\rm TSOR}\\,\\vert\\,$'',$)')
          IF (ICENTT.EQ.1) 
     *        WRITE(NRESUL,'(''${\\rm TCNT}\\,\\vert\\,$'',$)')
          IF (ITENSR.EQ.1) 
     *        WRITE(NRESUL,'(''${\\rm TTOT}\\,\\vert\\,$'',$)')
          WRITE(NRESUL,'()')
      END IF
C_______________________________________________________________________
C
      IF (IF_PAI.EQ.1) THEN
          WRITE(NRESUL,'(15X,''~~~${\\rm DFACTO_p}='',F6.2,'',\\,'',
     *                   ''\\Delta^2_p='',f6.2,''~{\\rm MeV^2},\\,'',
     *                   ''\\lambda_p= '',f6.2,''~{\\rm MeV}$'')')
     *                       DFACTO(INUCLI,1),DELT_2(INUCLI,1),
     *                                 XLAMBD_PROTON(INUCLI)
C
          WRITE(NRESUL,'(15X,''~~~${\\rm DFACTO_n}='',F6.2,'',\\,'',
     *                   ''\\Delta^2_n='',f6.2,''~{\\rm MeV^2},\\,'',
     *                   ''\\lambda_n= '',f6.2,''~{\\rm MeV}$'')')
     *                       DFACTO(INUCLI,2),DELT_2(INUCLI,2),
     *                                 XLAMBD_NEUTRS(INUCLI)
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
      IF (ISPLIT.EQ.1) THEN
          I_COLM=1
          I_LEFT=0
          I_RIGH=0
      END IF
C
      IF (ISPLIT.EQ.0) THEN
          I_COLM=0
          I_LEFT=1
          I_RIGH=0
      END IF
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
C      
      WRITE(NRESUL,'(/,''<<NO_OF_CURV>>'',1X,''NUMB_CURVE'',2X,
     *                 ''IF_SPLITIN'',2X,''IF_BULLETS'',2X,
     *                 ''IF_HISTOGS'',2X,''IF_OVERFITING'')')
      WRITE(NRESUL,'(15X,I5,7X,I5,7X,I5,7X,I5,7X,I5)')N_CURV,ISPLIT,
     *                                         IFBULL,IFHIST,IFOVER
C
C=======================================================================
C
      ISTYLE=0
      J_TYPE=I_TYPE
      
      DO I_CURV=1,N_CURV
          
         ICOLOR_AUXIL=ICOLOR(I_CURV)
         
         IF (ISPLIT.EQ.1) THEN
         
             ISTYLE=2
         
             IF (MOD(I_CURV+1,3).EQ.0) ISTYLE=1
             IF (MOD(I_CURV,3).EQ.0) ISTYLE=0
         
             IF (I_CURV.EQ.N_CURV) THEN
                 ISTYLE=0
                 ICOLOR(I_CURV)=0
             END IF
         
         END IF
         
         IF (ISPLIT.EQ.0) THEN
             
             ICOLOR(I_CURV)=ICOLOR(I_CURV+2*(I_CURV-1)) !To choose the
C                                                        dark color of
C                                                        of the triplet            
         END IF
         
         WRITE(NRESUL,'(''<<>>'')')
         WRITE(NRESUL,'(''<<CURVE_'',I4.4,''>>'',1X,
     *                  ''THICK_LINE'',2X,''STYLE_LINE'',2X,
     *                  ''COLOR_LINE'',2X,''TYPE_POINT'')')I_CURV
         WRITE(NRESUL,'(15X,I5,7X,I5,7X,I5,7X,I5)')ITHICK(I_CURV),
     *                                             ISTYLE,
     *                                             ICOLOR(I_CURV),
     *                                             J_TYPE
         WRITE(NRESUL,'(15X,''NUMB_POINT'',/,15X,I5)')ND_MAX(I_CURV)
         
         DO IPOINT=1,ND_MAX(I_CURV)
            
            WRITE(NRESUL,'(15X,F10.4,4X,E12.4,4X,A100)') 
     *                           X_CURV(IPOINT,I_CURV),
     *                           Y_CURV(IPOINT,I_CURV),
     *                           L_CURV(IPOINT,I_CURV)
            
         END DO
         
         J_TYPE=J_TYPE+2
C         
         ICOLOR(I_CURV)=ICOLOR_AUXIL
         
      END DO
C
C=======================================================================
C      
      WRITE(NRESUL,'(/,''<<GO_GETTHEM>>'',/)')
C
C=======================================================================
C
      IF (LOGWRI.GT.0) WRITE(LOGFIL,'(''Exiting DENSIT_OUTPUT'',/)')
C
C=======================================================================
C      
      RETURN
      END
C
C=======================================================================
C=======================================================================
C
      SUBROUTINE WRITIN_CHIMAP(NDMESH,LDRAND,ICOUNT_OFPARS,
     *                         TITPAR,TITMSH,RMSPRO_MINIMI,
     *                         RMSPRO_MAXIMA,RMSNEU_MINIMI,
     *                         RMSNEU_MAXIMA,MAXMSH,AXSMAP,
     *                         RMSMAP_PRONUC,RMSMAP_NEUNUC,
     *                         PARPOT_OFMESH,INDEXI,INUCLI)
C     
      INCLUDE  'MATDIM/NDNUCL.f'
      INCLUDE  'MATDIM/NDPARS.f'
      INCLUDE  'MATDIM/NDTITL.f'
      INCLUDE  'MATDIM/NDIM_M.f'
C      
      CHARACTER
     *          FILNAM*300,FILNAM_ENDING*3,INPSYM*6,TEXTRA*256
      CHARACTER
     *          TITLES*12,TITPAR*13,TITMSH*8,TITAUX*7
      CHARACTER
     *          FILNA2*256,FITNUC*2,FILNAM_NUCFIT*30,
     *          FITNUC_AUXILI*30,FILNA3*256,TITLES_NOTFIT*256,
     *          FILNAM_NOTTIT*256,TEXLAM*20
      CHARACTER
     *          TITLES_LATEXS*050,NUCNAM_TITLES*010
      CHARACTER
     *          TYPCHI*6,TAKCHI*3
      CHARACTER 
     *          NUCSYM*6,NUCSYM_AUXILI*6
      DIMENSION 
     *          TITPAR(1:NDPARS),
     *          TITMSH(1:NDPARS),
     *          TITAUX(1:NDPARS)
      DIMENSION
     *          TAKPAR(1:NDPARS),
     *          TAKCHI(1:NDNUCL)
      DIMENSION
     *          MAXMSH(1:NDPARS),
     *          AXSMAP(1:NDMESH,1:NDMESH)
      DIMENSION
     *          CHITOT_OFMESH(1:NDMESH,1:NDMESH),
     *          CHIPRO_OFMESH(1:NDMESH,1:NDMESH),
     *          CHINEU_OFMESH(1:NDMESH,1:NDMESH)
      DIMENSION
     *          PARPOT_OFMESH(1:NDMESH,1:NDMESH,1:NDPARS),
     *          PARPOT_AUXILI(1:NDPARS)
      DIMENSION
     *          IAUXIL(1:NDPARS),
     *          INDEXI(1:NDIM_M)
      DIMENSION
     *          FITNUC(1:NDNUCL),
     *          FITNUC_AUXILI(1:NDNUCL),
     *          TITLES_NOTFIT(1:NDPARS)
      DIMENSION
     *          RMSMAP_PRONUC(1:NDMESH,1:NDMESH,1:NDNUCL),
     *          RMSMAP_NEUNUC(1:NDMESH,1:NDMESH,1:NDNUCL)
C
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
     *       /MESHIN/ I_MESH(1:NDPARS),
     *                XMIN_I(1:NDPARS),
     *                XMAX_I(1:NDPARS),
     *                MAXPAR(1:NDPARS)
      COMMON
     *       /VLIMIT/ VMISTR(1:NDPARS),
     *                VMIMIN(1:NDPARS),
     *                VMIMAX(1:NDPARS),
     *                VMISTP(1:NDPARS)
      COMMON
     *       /NUCLEU/ IZ_FIX,IN_FIX,ISOSPI 
      COMMON
     *       /DENSTR/ IFDENS
      COMMON
     *       /TENSOR/ IFTENS
      COMMON
     *       /TENOPT/ ICENTT,ISORBT,ITENSR 
      COMMON
     *       /P_N_LM/ IFPROT,IFNEUT,IFBOTH
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
     *       /CHITAK/ IF_SPE,IF_RAD,IF_GAP,IF_FER,IF_DEN,IF_RHO,IF_INV,
     *                LDCHI2,TYPCHI(1:NDNUCL),ITKCHI(1:NDNUCL)
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
C                Writing MAP results in a file  
C=======================================================================
C
      IF (LOGWRI.GT.0) WRITE(LOGFIL,'(''Entering WRITING_CHIMAP'')')
C
C=======================================================================
C
      INDEXP=0
C      
      IF (IFDENS.EQ.0 .AND. IFPROT.EQ.1) THEN
          DO IPARAM=1,6
              INDEXP=INDEXP+1
              IAUXIL(INDEXP)=IPARAM
          END DO
      END IF
C      
      IF (IFDENS.EQ.0 .AND. IFNEUT.EQ.1) THEN
          DO IPARAM=21,26
              INDEXP=INDEXP+1
              IAUXIL(INDEXP)=IPARAM
          END DO
      END IF
C
      IF (IFDENS.EQ.0 .AND. IFBOTH.EQ.1) THEN
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
          IF (IFTENS.EQ.1 .AND. ISORBT.EQ.1) THEN
              DO IPARAM=43,46     !lambdas Tensor-SO
                 INDEXP=INDEXP+1
                 IAUXIL(INDEXP)=IPARAM
              END DO
          END IF
C
          IF (IFTENS.EQ.1 .AND. ICENTT.EQ.1) THEN
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
      N_UNIT=30
C
C=======================================================================
C      
      IF (ISOSPI.EQ.1 .OR. ISOSPI.EQ.11) FILNAM_ENDING='P'
      IF (ISOSPI.EQ.0 .OR. ISOSPI.EQ.10) FILNAM_ENDING='N'
C
      IF (IFDENS.EQ.1) THEN
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
      IF (ISOSPI.LE.1) THEN 
          NUCSYM_AUXILI=NUCSYM(INUCLI)
      ELSE
          NUCSYM_AUXILI='Global'
      END IF
C
      I3=6
      I4=1
      IF (NUCSYM_AUXILI(1:1).EQ.' ') I4=2
      IF (NUCSYM_AUXILI(2:2).EQ.' ') I4=3
      IF (NUCSYM_AUXILI(3:3).EQ.' ') I4=4
C
C=======================================================================
C     
      IF (IFDENS.EQ.0) THEN
C      
          IF (MACTIV.EQ.0)
     *        WRITE(FILNAM,'(''Chi2Maps/IFDENS-0/Chi2Map_'',
     *                        A,''_'',A1,''_'',A,''_'',A8,''_'',A8,
     *                      ''_IFDENS-'',I1,''_IFTENS-'',I1,
     *                      ''_IF-PAI-'',I1,
     *                      ''_IF-RAD-'',I1,''_IF-ORD-'',I1,
     *                      ''_IF-RHO-'',I1,
     *                      ''_IFDEEP-'',I1,''_IFPRON-'',I1,
     *                      ''_IFK-VC-'',I1,''_IFK-VS-'',I1,
     *                      ''_IFK-RC-'',I1,''_IFK-RS-'',I1,
     *                      ''_IFK-AC-'',I1,''_IFK-AS-'',I1,
     *                      ''_LDRAND-'',I5.5,''_'',A3,''.dat'')') 
     *
     *                      NUCSYM_AUXILI(I4:I3),
     *                      FILNAM_ENDING,FILNAM_NUCFIT(1:I1),
     *                      TITMSH(1),TITMSH(2),
     *                      IFDENS,IFTENS,IF_PAI,IF_RAD,IF_INV,
     *                      IF_RHO,IFDEEP,IFPRON,IFK_VC,IFK_VS,
     *                      IFK_RC,IFK_RS,IFK_AC,IFK_AS,LDRAND,
     *                                                  VERSIO
C
          IF (MACTIV.GT.0) 
     *        WRITE(FILNAM,'(''Chi2Maps/IFDENS-0/Chi2Map_'',
     *                     A,''_'',A1,''_'',A,''_'',A8,''_'',A8,
     *                   ''_IFDENS-'',I1,''_IFTENS-'',I1,
     *                   ''_IF-PAI-'',I1,
     *                   ''_IF-RAD-'',I1,''_IF-ORD-'',I1,
     *                   ''_IF-RHO-'',I1,
     *                   ''_IFDEEP-'',I1,''_IFPRON-'',I1,
     *                   ''_IFK-VC-'',I1,''_IFK-VS-'',I1,
     *                   ''_IFK-RC-'',I1,''_IFK-RS-'',I1,
     *                   ''_IFK-AC-'',I1,''_IFK-AS-'',I1,
     *                   ''_LDRAND-'',I5.5,''_'',A,''_'',A3,''.dat'')') 
     *
     *                   NUCSYM_AUXILI(I4:I3),
     *                   FILNAM_ENDING,FILNAM_NUCFIT(1:I1),
     *                   TITMSH(1),TITMSH(2),
     *                   IFDENS,IFTENS,IF_PAI,IF_RAD,IF_INV,
     *                   IF_RHO,IFDEEP,IFPRON,IFK_VC,IFK_VS,
     *                   IFK_RC,IFK_RS,IFK_AC,IFK_AS,LDRAND,
     *                   FILNAM_NOTTIT(1:I2),VERSIO
C
      END IF
C
C=======================================================================
C     
      IF (IFDENS.EQ.1) THEN
C      
          IF (MACTIV.EQ.0) 
     *        WRITE(FILNAM,'(''Chi2Maps/IFDENS-1_IFTENS-'',I1,''/'',
     *                   ''Chi2Map_'',A,''_'',A1,''_'',A,
     *                   ''_'',A8,''_'',A8,
     *                   ''_IFDENS-'',I1,''_IFTENS-'',I1,
     *                   ''_IF-PAI-'',I1,
     *                   ''_IF-RAD-'',I1,''_IF-ORD-'',I1,
     *                   ''_IF-RHO-'',I1,
     *                   ''_IFDEEP-'',I1,''_IFPRON-'',I1,
     *                   ''_IFK-VC-'',I1,''_IFK-VS-'',I1,
     *                   ''_IFK-RC-'',I1,''_IFK-RS-'',I1,
     *                   ''_IFK-AC-'',I1,''_IFK-AS-'',I1,
     *                   ''_LDRAND-'',I5.5,A15,''_'',A3,''.dat'')') 
     *
     *                   IFTENS,NUCSYM_AUXILI(I4:I3),
     *                   FILNAM_ENDING,FILNAM_NUCFIT(1:I1),
     *                   TITMSH(1),TITMSH(2),
     *                   IFDENS,IFTENS,IF_PAI,IF_RAD,IF_INV,
     *                   IF_RHO,IFDEEP,IFPRON,IFK_VC,IFK_VS,
     *                   IFK_RC,IFK_RS,IFK_AC,IFK_AS,LDRAND,
     *                                        TEXLAM,VERSIO
C
          IF (MACTIV.GT.0) 
     *        WRITE(FILNAM,'(''Chi2Maps/IFDENS-1_IFTENS-'',I1,''/'',
     *                   ''Chi2Map_'',A,''_'',A1,''_'',A,
     *                   ''_'',A8,''_'',A8,
     *                   ''_IFDENS-'',I1,''_IFTENS-'',I1,
     *                   ''_IF-PAI-'',I1,
     *                   ''_IF-RAD-'',I1,''_IF-ORD-'',I1,
     *                   ''_IF-RHO-'',I1,
     *                   ''_IFDEEP-'',I1,''_IFPRON-'',I1,
     *                   ''_IFK-VC-'',I1,''_IFK-VS-'',I1,
     *                   ''_IFK-RC-'',I1,''_IFK-RS-'',I1,
     *                   ''_IFK-AC-'',I1,''_IFK-AS-'',I1,
     *                   ''_LDRAND-'',I5.5,''_'',A,A15,''_'',A3,
     *                                                  ''.dat'')') 
     *
     *                   IFTENS,NUCSYM_AUXILI(I4:I3),
     *                   FILNAM_ENDING,FILNAM_NUCFIT(1:I1),
     *                   TITMSH(1),TITMSH(2),
     *                   IFDENS,IFTENS,IF_PAI,IF_RAD,IF_INV,
     *                   IF_RHO,IFDEEP,IFPRON,IFK_VC,IFK_VS,
     *                   IFK_RC,IFK_RS,IFK_AC,IFK_AS,LDRAND,
     *                   FILNAM_NOTTIT(1:I2),TEXLAM,VERSIO
C
          IF (IFTENS.EQ.1) THEN
C      
              IF (MACTIV.EQ.0) 
     *            WRITE(FILNAM,'(''Chi2Maps/IFDENS-1_IFTENS-'',I1,''/'',
     *                   ''Chi2Map_'',A,''_'',A1,''_'',A,
     *                   ''_'',A8,''_'',A8,
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
     *                   ''_LDRAND-'',I5.5,A15,''_'',A3,''.dat'')') 
     *
     *                   IFTENS,NUCSYM_AUXILI(I4:I3),
     *                   FILNAM_ENDING,FILNAM_NUCFIT(1:I1),
     *                   TITMSH(1),TITMSH(2),
     *                   IFDENS,IFTENS,IF_PAI,ICENTT,ISORBT,ITENSR,
     *                   IF_RAD,IF_INV,IF_RHO,IFDEEP,IFPRON,IFK_VC,
     *                   IFK_VS,IFK_RC,IFK_RS,IFK_AC,IFK_AS,LDRAND,
     *                                               TEXLAM,VERSIO
C
              IF (MACTIV.GT.0) 
     *            WRITE(FILNAM,'(''Chi2Maps/IFDENS-1_IFTENS-'',I1,''/'',
     *                   ''Chi2Map_'',A,''_'',A1,''_'',A,
     *                   ''_'',A8,''_'',A8,
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
     *                                                 ''.dat'')') 
     *
     *                   IFTENS,NUCSYM_AUXILI(I4:I3),
     *                   FILNAM_ENDING,FILNAM_NUCFIT(1:I1),
     *                   TITMSH(1),TITMSH(2),
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
      OPEN (N_UNIT,FILE=FILNAM,STATUS='UNKNOWN',FORM='FORMATTED')
C
C=======================================================================
C      
      WRITE(N_UNIT,'(/,''<<Strasbourg Nuclear Inverse Problem Study '',
     *                 ''- Standard Output - WS-SPH>>'',/)')
C
C=======================================================================
C
      WRITE(N_UNIT,'(''<<VersionOfCode>>'',3X,''VERSIO'')')
C
      WRITE(N_UNIT,'(21X,A3,/)') VERSIO
C
C=======================================================================
C     
      WRITE(N_UNIT,'(''<<Nuclid_Inform>>'',3X,''NUMB_Z'',
     *                                     3X,''NUMB_N'')')
C
      WRITE(N_UNIT,'(19X,I3,6X,I3,/)')IZ_FIX,IN_FIX
C
C=======================================================================
C
      WRITE(N_UNIT,'(''<<HowManyNuclei>>'',3X,''NUCACT'')')   
      
      WRITE(N_UNIT,'(19X,I4,/)')NUCACT
C
C=======================================================================
C
      WRITE(N_UNIT,'(''<<Nuclei_Inform>>'',3x,
     *               <NUCACT>(''IZ_No'',i1,3x,
     *               ''IN_No'',i1,5x))')(I,I=1,NUCACT),(J,J=1,NUCACT)
      WRITE(N_UNIT,'(20X,<NUCACT>(I4,5X,I4,7X),/)')
     *                                 (INPUTZ(I),INPUTN(I),I=1,NUCACT)
C
C=======================================================================
C
      WRITE(N_UNIT,'(''<<Experimt_File>>'',3X,''IFDEEP'',2X
     *                                     ''IFPRON'')')
      
      WRITE(N_UNIT,'(19X,2(I6,3X),/)') IFDEEP,IFPRON
C
C=======================================================================
C
      ISOSPI_PRINTI=ISOSPI
      IF (ISOSPI.GT.1) ISOSPI_PRINTI=ISOSPI-10
C
      WRITE(N_UNIT,'(''<<PotentialInfo>>'',3x,
     *               ''ISOSPI   IFDENS   IFTENS   IF_PAI'')')
      
      WRITE(N_UNIT,'(19X,I6,3X,I6,3X,I6,3X,I6,3X,/)')ISOSPI_PRINTI,
     *                                        IFDENS,IFTENS,IF_PAI
C
C=======================================================================
C
      WRITE(N_UNIT,'(''<<TensorialInfo>>'',3X,''ICENTT'',3X
     *                                     ''ISORBT'',3X,''ITENSR'')')
      
      WRITE(N_UNIT,'(19X,3(I6,3X),/)') ICENTT,ISORBT,ITENSR
C
C=======================================================================
C
      WRITE(N_UNIT,'(''<<KappasParamet>>   IFK_VC   IFK_VS   '',
     *                      ''IFK_RC   IFK_RS   IFK_AC   IFK_AS'')')
      WRITE(N_UNIT,'(19X,6(I6,3X),/)') IFK_VC,IFK_VS,IFK_RC,IFK_RS,
     *                                               IFK_AC,IFK_AS
C
C=======================================================================
C      
      WRITE(N_UNIT,'(''<<MinimizedChi2>>'',3X,
     *                     <LDCHI2>(A6,3X))')(TYPCHI(I),I=1,LDCHI2)
     
      WRITE(N_UNIT,'(19X,<LDCHI2>(I6,3X),/)')(ITKCHI(I),I=1,LDCHI2)
C
C=======================================================================
C
      WRITE(N_UNIT,'(''<<HowManyParams>>'',3X,''ICOUNT'',
     *                                     3X,''LDPARS'')')
     
      WRITE(N_UNIT,'(19X,I6,3X,I6,/)')ICOUNT_OFPARS,INDEXP
C
C=======================================================================
C
      DO I=1,INDEXP
         IPARAM=IAUXIL(I)
         IF (IFTAKE(IPARAM).EQ.0 .AND. I_MESH(IPARAM).EQ.0) THEN
             IFTAKE(IPARAM)=2
         END IF
      END DO
C_______________________________________________________________________
C
      WRITE(N_UNIT,'(''<<MinimizedOver>>'',$)')
      
      DO I=1,INDEXP
          IPARAM=IAUXIL(I)
          WRITE(N_UNIT,'(A11,$)')TITLES(IPARAM)(3:10)
      END DO
      
      WRITE(N_UNIT,'()')
      
      WRITE(N_UNIT,'(16X,'' '',$)')
      
      DO I=1,INDEXP
          IPARAM=IAUXIL(I)
          WRITE(N_UNIT,'(I11,$)') IFTAKE(IPARAM)
      END DO
      
      WRITE(N_UNIT,'()')
      WRITE(N_UNIT,'()')
C
C=======================================================================
C
      IF (MACTIV.NE.0) THEN  !which means that there are fixed parameters
C
          WRITE(N_UNIT,'(''<<Fixed_PValues>>'',$)')
C
          DO I=1,INDEXP
             IPARAM=IAUXIL(I)
             IF (IFTAKE(IPARAM).EQ.2 .AND. I_MESH(IPARAM).EQ.0) THEN
                 WRITE(N_UNIT,'(A11,$)')TITLES(IPARAM)(3:10)
             END IF
          END DO
C
          WRITE(N_UNIT,'()')
C
          WRITE(N_UNIT,'(16X,'' '',$)')
C
          DO I=1,INDEXP
             IPARAM=IAUXIL(I)
             IF (IFTAKE(IPARAM).EQ.2 .AND. I_MESH(IPARAM).EQ.0) THEN
                 WRITE(N_UNIT,'(F11.4,$)')VMISTR(IPARAM)
             END IF
          END DO
C
          WRITE(N_UNIT,'()')
          WRITE(N_UNIT,'()')
C
      END IF
C_______________________________________________________________________
C
      DO I=1,INDEXP
         IPARAM=IAUXIL(I)
         IF (IFTAKE(IPARAM).EQ.2) THEN
             IFTAKE(IPARAM)=0
         END IF
      END DO
C
C=======================================================================
C
      WRITE(N_UNIT,'(''<<X--Axis--Mult>>'',3X,''ParamI'')')
      
      WRITE(N_UNIT,'(20X,A,3X,/)')TITMSH(1)  
C
C=======================================================================
C
      WRITE(N_UNIT,'(''<<Y--Axis--Mult>>'',3X,''ParamJ'')')
      
      WRITE(N_UNIT,'(20X,A,3X,/)')TITMSH(2)  
C
C=======================================================================
C
      WRITE(N_UNIT,'(''<<ChiSqrt_Extrm>>'',3X,''PROMIN'',4X,''PROMAX'',
     *                                   4X,''NEUMIN'',4X,''NEUMAX'')')
      
      WRITE(N_UNIT,'(20X,4(ES8.2,2X),/)')RMSPRO_MINIMI,RMSPRO_MAXIMA,
     *                                   RMSNEU_MINIMI,RMSNEU_MAXIMA
C
C=======================================================================
C      
      WRITE(N_UNIT,'(''<<Plot--X--axis>>'',3X,''NofPnts'')')
      
      WRITE(N_UNIT,'(21X,I3)')MAXMSH(1)
      
      WRITE(N_UNIT,'(20X,''Values'')')
      
      WRITE(N_UNIT,'(18X,<MAXMSH(1)>F11.4,/)')
     *              (AXSMAP(1,K1),K1=1,MAXMSH(1))
C
C=======================================================================
C     
      WRITE(N_UNIT,'(''<<Plot--Y--axis>>'',3X,''NofPnts'')')
      
      WRITE(N_UNIT,'(21X,I3)')MAXMSH(2)
      
      WRITE(N_UNIT,'(20X,''Values'')')
      
      WRITE(N_UNIT,'(18X,<MAXMSH(2)>F11.4,/)')
     *              (AXSMAP(2,K2),K2=1,MAXMSH(2))
C
C=======================================================================
C
C     Writing the results obtained for the PROTON RMS
C
      IF (ISOSPI.EQ.1 .OR. ISOSPI.EQ.11) THEN
C
          WRITE(N_UNIT,'(''<<RMSDevi_Protn>>'')')
      
          WRITE(N_UNIT,'(10X,<MAXMSH(1)>F12.4)')
     *              (AXSMAP(1,K1),K1=1,MAXMSH(1))
      
          DO K2=1,MAXMSH(2)
        
             WRITE(N_UNIT,'(F10.4,$)')AXSMAP(2,K2)
        
             DO K1=1,MAXMSH(1)
         
                WRITE(N_UNIT,'(ES12.4,$)')RMSMAP_PRONUC(K1,K2,INUCLI)
        
             END DO
        
             WRITE(N_UNIT,'()')
       
          END DO
      
          WRITE(N_UNIT,'()')
C
      END IF
C
C=======================================================================
C
C     Writing the results obtained for the NEUTRON RMS
C
      IF (ISOSPI.EQ.0 .OR. ISOSPI.EQ.10) THEN    
C 
          WRITE(N_UNIT,'(''<<RMSDevi_Neutr>>'')')
      
          WRITE(N_UNIT,'(10X,<MAXMSH(1)>F12.4)')
     *              (AXSMAP(1,K1),K1=1,MAXMSH(1))
      
          DO K2=1,MAXMSH(2)
        
             WRITE(N_UNIT,'(F10.4,$)')AXSMAP(2,K2)
        
             DO K1=1,MAXMSH(1)
         
                WRITE(N_UNIT,'(ES12.4,$)')RMSMAP_NEUNUC(K1,K2,INUCLI)
        
             END DO
        
             WRITE(N_UNIT,'()')
       
          END DO
      
          WRITE(N_UNIT,'()') 
C
      END IF
C
C=======================================================================
C
C     Writing the results obtained for the parameters with ITAKE.EQ.1
C     
      DO IPARAM=1,NDPARS
         IF (IFTAKE(IPARAM).EQ.1) THEN
C
             WRITE(N_UNIT,'(''<<Parm_'',A8,''>>'')')TITLES(IPARAM)(3:10)
C
             WRITE(N_UNIT,'(10X,<MAXMSH(1)>F12.4)')
     *                          (AXSMAP(1,K1),K1=1,MAXMSH(1))
C
             DO K2=1,MAXMSH(2)
C
                WRITE(N_UNIT,'(F10.4,$)')AXSMAP(2,K2)
C
                DO K1=1,MAXMSH(1)
C
                   WRITE(N_UNIT,'(ES12.4,$)')PARPOT_OFMESH(K1,K2,IPARAM)
C
                END DO
C
                WRITE(N_UNIT,'()')
C
             END DO
C
             WRITE(N_UNIT,'()')
C
         END IF
      END DO
C           
C=======================================================================
C 
      WRITE(N_UNIT,'(''<<EXECUTESET_GO>>'')')      
C
C=======================================================================
C
      IF (LOGWRI.GT.0) WRITE(LOGFIL,'(''Exiting WRITING_CHIMAP'')')
C
C=======================================================================
C      
      RETURN
      END
C
C=======================================================================
C=======================================================================
C
      SUBROUTINE OPENIN_MESHPM(LDRAND,TITMSH,ICOUNT_OFPARS,
     *                                NTOTAL_OFMESH,MAXMSH)
C
      INCLUDE   'MATDIM/NDPARS.f'
      INCLUDE   'MATDIM/NDNUCL.f'
      INCLUDE   'MATDIM/NDMESH.f'
      INCLUDE   'MATDIM/NDTITL.f'
C      
      PARAMETER
     *          (NDMES2=NDMESH*NDMESH)
      CHARACTER
     *          TITLES*12,TITAUX*12,TITMSH*8,TEXTRA*256
      CHARACTER
     *          FILNAM*256,FILNAM_ENDING*1,INPSYM*6,NUCSYM*6,TYPCHI*6
      CHARACTER
     *          FILNAM_NUCFIT*30,FITNUC*2,
     *          FITNUC_AUXILI*30,TITLES_NOTFIT*256,FILNAM_NOTTIT*256,
     *          TEXLAM*20,VERSIO*3
      CHARACTER
     *          TITLES_LATEXS*050,NUCNAM_LATEXS*010
      DIMENSION
     *          TITAUX(1:NDPARS),
     *          TITMSH(1:NDPARS)
      DIMENSION
     *          IAUXIL(1:NDPARS)
      DIMENSION
     *          MAXMSH(1:NDPARS)
      DIMENSION
     *          FITNUC(1:NDNUCL),
     *          FITNUC_AUXILI(1:NDNUCL),
     *          TITLES_NOTFIT(1:NDPARS)
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
     *       /DENSTR/ IFDENS
      COMMON
     *       /TENSOR/ IFTENS
      COMMON
     *       /TENOPT/ ICENTT,ISORBT,ITENSR
      COMMON
     *       /IFPAIR/ IF_PAI,IFGPAI,IFDPAI
     *       /DEPPAI/ IFDEEP,IFPRON
      COMMON
     *       /SEARCH/ IFTAKE(1:NDPARS)
      COMMON
     *       /MESHIN/ I_MESH(1:NDPARS),
     *                XMIN_I(1:NDPARS),
     *                XMAX_I(1:NDPARS),
     *                MAXPAR(1:NDPARS)
      COMMON
     *       /MSHPRI/ INOFER_PRINTI,
     *                IDEFCN_PRINTI,
     *                CHISQU_PRINTI(1:NDMES2),
     *                CHIGRD_NRMPRI(1:NDMES2),
     *                CHIGRD_COMPRI(1:NDPARS,1:NDMES2)
      COMMON
     *       /P_N_LM/ IFPROT,IFNEUT,IFBOTH
      COMMON
     *       /IFKAPP/ IFK_VC,IFK_VS,IFK_RC,IFK_RS,IFK_AC,IFK_AS
      COMMON
     *       /CHITAK/ IF_SPE,IF_RAD,IF_GAP,IF_FER,IF_DEN,IF_RHO,IF_INV,
     *                LDCHI2,TYPCHI(1:NDNUCL),ITKCHI(1:NDNUCL)
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
     *       /VERSIN/ VERSIO
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
      IF (LOGWRI.GT.0) WRITE(LOGFIL,'(''Entering OPENIN_MESHP'')')
C
C=======================================================================
C
      INDEXP=0
C      
      IF (IFDENS.EQ.0 .AND. IFPROT.EQ.1) THEN
          ISOSPI=1
          DO IPARAM=1,6
              INDEXP=INDEXP+1
              IAUXIL(INDEXP)=IPARAM
          END DO
      END IF
C      
      IF (IFDENS.EQ.0 .AND. IFNEUT.EQ.1) THEN
          ISOSPI=0
          DO IPARAM=21,26
              INDEXP=INDEXP+1
              IAUXIL(INDEXP)=IPARAM
          END DO
      END IF
C
      IF (IFDENS.EQ.0 .AND. IFBOTH.EQ.1) THEN
C
          ISOSPI=2
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
          ISOSPI=2
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
      LOGMSH=39
C
C=======================================================================
C      
      IF (IFDENS.EQ.0.AND.IFPROT.EQ.1) FILNAM_ENDING='P'
      IF (IFDENS.EQ.0.AND.IFNEUT.EQ.1) FILNAM_ENDING='N'
      IF (IFDENS.EQ.0.AND.IFBOTH.EQ.1) FILNAM_ENDING='B'
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
      DO I=1,INDEXP
         IPARAM=IAUXIL(I)
         IF (IFTAKE(IPARAM).EQ.0.AND.I_MESH(IPARAM).EQ.0) THEN
             MACTIV=MACTIV+1
             TITLES_NOTFIT(MACTIV)=TITLES(IPARAM)(3:10)
         END IF
      END DO
C
      I2=9*MACTIV-1
C
      WRITE(FILNAM_NOTTIT,'(<MACTIV>(A8,''-''))')
     *                      (TITLES_NOTFIT(I),I=1,MACTIV)
C
C=======================================================================
C     
      IF (IFDENS.EQ.0) THEN
C      
          IF (MACTIV.EQ.0) 
     *        WRITE(FILNAM,'(''MeshTables/IFDENS-0/'',A,''_'',A1,
     *                   ''_'',A8,''_'',A8,
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
     *                   TITMSH(1),TITMSH(2),
     *                   IFDENS,IFTENS,IF_PAI,IF_RAD,IF_INV,
     *                   IF_RHO,IFDEEP,IFPRON,IFK_VC,IFK_VS,
     *                   IFK_RC,IFK_RS,IFK_AC,IFK_AS,LDRAND,
     *                                               VERSIO
C
          IF (MACTIV.GT.0) 
     *        WRITE(FILNAM,'(''MeshTables/IFDENS-0/'',A,''_'',A1,
     *                   ''_'',A8,''_'',A8,
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
     *                   TITMSH(1),TITMSH(2),
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
     *        WRITE(FILNAM,'(''MeshTables/IFDENS-1_IFTENS-'',I1,''/'',
     *                      A,''_'',A1,''_'',A8,''_'',A8,
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
     *                   TITMSH(1),TITMSH(2),
     *                   IFDENS,IFTENS,IF_PAI,IF_RAD,IF_INV,
     *                   IF_RHO,IFDEEP,IFPRON,IFK_VC,IFK_VS,
     *                   IFK_RC,IFK_RS,IFK_AC,IFK_AS,LDRAND,TEXLAM,
     *                                                      VERSIO
C
          IF (MACTIV.GT.0) 
     *        WRITE(FILNAM,'(''MeshTables/IFDENS-1_IFTENS-'',I1,
     *                   ''/'',A,''_'',A1,''_'',A8,''_'',A8,
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
     *                   TITMSH(1),TITMSH(2),
     *                   IFDENS,IFTENS,IF_PAI,IF_RAD,IF_INV,
     *                   IF_RHO,IFDEEP,IFPRON,IFK_VC,IFK_VS,
     *                   IFK_RC,IFK_RS,IFK_AC,IFK_AS,LDRAND,
     *                   FILNAM_NOTTIT(1:I2),TEXLAM,VERSIO
C
          IF (IFTENS.EQ.1) THEN
C      
              IF (MACTIV.EQ.0) 
     *            WRITE(FILNAM,'(''MeshTables/IFDENS-1_IFTENS-'',I1,
     *                   ''/'',A,''_'',A1,''_'',A8,''_'',A8,
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
     *                   TITMSH(1),TITMSH(2),
     *                   IFDENS,IFTENS,IF_PAI,ICENTT,ISORBT,ITENSR,
     *                   IF_RAD,IF_INV,IF_RHO,IFDEEP,IFPRON,IFK_VC,
     *                   IFK_VS,IFK_RC,IFK_RS,IFK_AC,IFK_AS,LDRAND,
     *                                               TEXLAM,VERSIO
C
              IF (MACTIV.GT.0) 
     *            WRITE(FILNAM,'(''MeshTables/IFDENS-1_IFTENS-'',I1,
     *                   ''/'',A,''_'',A1,''_'',A8,''_'',A8,
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
     *                   TITMSH(1),TITMSH(2),
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
      OPEN (LOGMSH,FILE=FILNAM,STATUS='UNKNOWN',FORM='FORMATTED')
C
C=======================================================================
C
      WRITE(LOGMSH,'(/,''<<HowManyNuclei>>'',3x,''NUCACT'')')
      WRITE(LOGMSH,'(20X,I4)')NUCACT
C_______________________________________________________________________
C
      WRITE(LOGMSH,'(/,''<<Nuclei_Inform>>'',3x,
     *               <NUCACT>(''IZ_No'',i1,3x,
     *               ''IN_No'',i1,5x))')(I,I=1,NUCACT),(J,J=1,NUCACT)
      WRITE(LOGMSH,'(20X,<NUCACT>(I4,5X,I4,7X))')
     *                                 (INPUTZ(I),INPUTN(I),I=1,NUCACT)
C_______________________________________________________________________
C
      WRITE(LOGMSH,'(/,''<<PotentialInfo>>'',3x,''ISOSPI'',
     *              3X,''IFDENS'',3x,''IFTENS'')')
      WRITE(LOGMSH,'(20X,I4,5X,I4,5X,I4)')ISOSPI,IFDENS,IFTENS
C_______________________________________________________________________
C
      WRITE(LOGMSH,'(/,''<<Chi2_Informat>>'',3x,''IF_SPE   IF_RAD   '',
     *               ''IF_GAP   IF_FER   IF_DEN   IF_RHO   IF_INV'')')
      WRITE(LOGMSH,'(20X,7(I4,5X))')IF_SPE,IF_RAD,IF_GAP,IF_FER,IF_DEN,
     *                                                   IF_RHO,IF_INV
C_______________________________________________________________________
C
      WRITE(LOGMSH,'(/,''<<NoOf_Restarts>>'',3x,''LDRAND'')')
      WRITE(LOGMSH,'(20X,I4)') LDRAND
C_______________________________________________________________________
C
      WRITE(LOGMSH,'(/,''<<Dim_MeshPoint>>'',3X,''MESH01'',
     *                                       3X,''MESH02'')')
      WRITE(LOGMSH,'(20X,I4,5X,I4)')MAXMSH(1),MAXMSH(2)
C_______________________________________________________________________
C
      INDAUX=0
C_______________________________________________________________________
C
      WRITE(LOGMSH,'(/,''<<Tabuld_Params>>'',$)')
C      
      DO I=1,INDEXP
         IPARAM=IAUXIL(I)
         IF (I_MESH(IPARAM).EQ.1) THEN
             WRITE(LOGMSH,'(A11,$)')TITLES(IPARAM)(3:10)
             INDAUX=INDAUX+1
             TITAUX(INDAUX)=TITLES(IPARAM)(3:10)
         END IF
      END DO
C      
      WRITE(LOGMSH,'()')
      WRITE(LOGMSH,'()')
C_______________________________________________________________________
C
      WRITE(LOGMSH,'(''<<Constn_Params>>'',$)')
C      
      DO I=1,INDEXP
         IPARAM=IAUXIL(I)
         IF (IFTAKE(IPARAM).EQ.0 .AND. I_MESH(IPARAM).EQ.0) THEN
             WRITE(LOGMSH,'(A11,$)')TITLES(IPARAM)(3:10)
             INDAUX=INDAUX+1
             TITAUX(INDAUX)=TITLES(IPARAM)(3:10)
         END IF
      END DO
C      
      WRITE(LOGMSH,'()')
      WRITE(LOGMSH,'()')
C_______________________________________________________________________
C
      WRITE(LOGMSH,'(''<<Fitted_Params>>'',$)')
C      
      DO I=1,INDEXP
         IPARAM=IAUXIL(I)
         IF (IFTAKE(IPARAM).EQ.1) THEN
             WRITE(LOGMSH,'(A11,$)')TITLES(IPARAM)(3:10)
             INDAUX=INDAUX+1
             TITAUX(INDAUX)=TITLES(IPARAM)(3:10)
         END IF
      END DO
C      
      WRITE(LOGMSH,'()')
      WRITE(LOGMSH,'()')
C_______________________________________________________________________
C
      WRITE(LOGMSH,'(''<<MinimSolution>>'',$)')
C      
      DO I=1,INDAUX
         WRITE(LOGMSH,'(3X,A8,''   '',$)') TITAUX(I)
      END DO
C      
      WRITE(LOGMSH,'(3X,''CHISQURT'',6X,''ITERATIN'',6X,''STOPPING'',
     *               6X,''GRADNORM   '',$)')
C      
      DO I=1,INDAUX
         WRITE(LOGMSH,'(3X,A8,''   '',$)') TITAUX(I)
      END DO
C      
      WRITE(LOGMSH,'()')
C
C=======================================================================
C
      IF (LOGWRI.GT.0) WRITE(LOGFIL,'(''Exiting OPENIN_MESHP'')')
C
C=======================================================================
C      
      RETURN
      END
C
C=======================================================================
C=======================================================================
C
      SUBROUTINE WRITIN_MESHPM(NCOUNT_OFMESH,PARPOT_PRINTI)
C
      INCLUDE   'MATDIM/NDPARS.f'
      INCLUDE   'MATDIM/NDNUCL.f'
      INCLUDE   'MATDIM/NDMESH.f'
      INCLUDE   'MATDIM/NDTITL.f'
C      
      PARAMETER
     *          (NDMES2=NDMESH*NDMESH)
      CHARACTER
     *          TITLES*12,TITAUX*12,TITMSH*8,TEXTRA*256
      CHARACTER
     *          FILNAM*256,FILNAM_ENDING*3,INPSYM*6,NUCSYM*6,TYPCHI*6,
     *          TITLES_LATEXS*050,NUCNAM_LATEXS*010
      DIMENSION
     *          TITAUX(1:NDPARS),
     *          TITMSH(1:NDPARS)
      DIMENSION
     *          IAUXIL(1:NDPARS),
     *          PARAUX(1:NDPARS,1:NDMES2),
     *          PARPOT_PRINTI(1:NDPARS,1:NDMES2),
     *          CHIGRD_COMAUX(1:NDPARS,1:NDMES2)
      DIMENSION
     *          MAXMSH(1:NDPARS)
      COMMON
     *       /DENSTR/ IFDENS
      COMMON
     *       /TENSOR/ IFTENS
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
     *       /P_N_LM/ IFPROT,IFNEUT,IFBOTH
      COMMON
     *       /IFKAPP/ IFK_VC,IFK_VS,IFK_RC,IFK_RS,IFK_AC,IFK_AS
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
      IF (LOGWRI.GT.0) THEN
          WRITE(LOGFIL,'(''Entering WRITIN_MESHPM with IFDENS= '',I1,
     *                   ''; IFPROT= '',I1,''; IFNEUT= '',I1,
     *                   ''; IFBOTH= '',I1)')IFDENS,IFPROT,
     *                                       IFNEUT,IFBOTH
      END IF
C
C=======================================================================
C
      INDEXP=0
C      
      IF (IFDENS.EQ.0 .AND. IFPROT.EQ.1) THEN
          DO IPARAM=1,6
              INDEXP=INDEXP+1
              IAUXIL(INDEXP)=IPARAM
          END DO
      END IF
C      
      IF (IFDENS.EQ.0 .AND. IFNEUT.EQ.1) THEN
          DO IPARAM=21,26
              INDEXP=INDEXP+1
              IAUXIL(INDEXP)=IPARAM
          END DO
      END IF
C
      IF (IFDENS.EQ.0 .AND. IFBOTH.EQ.1) THEN
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
          IF (IFTENS.EQ.1 .AND. ISORBT.EQ.1) THEN
              DO IPARAM=43,46     !lambdas Tensor-SO
                 INDEXP=INDEXP+1
                 IAUXIL(INDEXP)=IPARAM
              END DO
          END IF
C
          IF (IFTENS.EQ.1 .AND. ICENTT.EQ.1) THEN
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
      WRITE(LOGMSH,'(10X,''    '',$)')
      
      DO I=1,INDEXP
         IPARAM=IAUXIL(I)
         IF (I_MESH(IPARAM).EQ.1) THEN
             WRITE(LOGMSH,'(E14.4,$)') 
     *             PARPOT_PRINTI(IPARAM,NCOUNT_OFMESH)
         END IF
      END DO
      
      DO I=1,INDEXP
         IPARAM=IAUXIL(I)
         IF (I_MESH(IPARAM).EQ.0 .AND. IFTAKE(IPARAM).EQ.0) THEN
             WRITE(LOGMSH,'(E14.4,$)')  
     *             PARPOT_PRINTI(IPARAM,NCOUNT_OFMESH)
         END IF
      END DO
      
      DO I=1,INDEXP
         IPARAM=IAUXIL(I)
         IF (IFTAKE(IPARAM).EQ.1) THEN
             WRITE(LOGMSH,'(E14.4,$)')  
     *             PARPOT_PRINTI(IPARAM,NCOUNT_OFMESH)
         END IF
      END DO
      
      WRITE(LOGMSH,'(E14.4,6X,I4,10X,I4,4X,E14.4,$)') 
     *                                    CHISQU_PRINTI(NCOUNT_OFMESH),
     *                                    IDEFCN_PRINTI,INFOER_PRINTI,
     *                                    CHIGRD_NRMPRI(NCOUNT_OFMESH)
      
      DO I=1,INDEXP
         IPARAM=IAUXIL(I)
         WRITE(LOGMSH,'(E14.4,$)') CHIGRD_COMPRI(IPARAM,NCOUNT_OFMESH)
      END DO
      
      WRITE(LOGMSH,'()')
C
C=======================================================================
C
      IF (LOGWRI.GT.0) WRITE(LOGFIL,'(''Exiting WRITIN_MESHPM'')')
C
C=======================================================================
C      
      RETURN
      END
C      
C=======================================================================
C=======================================================================
C
      FUNCTION DENEXP_INTEGR(XARGUM)
C
      COMMON
     *       /ARGUMN/ IARG_R,IARG_Z
      COMMON
     *       /SIMUSE/ AOSCIL,RADCUT,INUCLI
      COMMON
     *       /PRFILE/ ISCREN,LOGFIL,LOGCHP,LOGCHN,LOGHIP,LOGHIN,
     *                LOGPRO,LOGNEU,LOGWRI,LSCREN,IRESUL,ICONVE,
     *                NOUTPT,LOGAUX,LOGENE,LOGRAD,LOGGRD,LOGMSH,
     *                                                   LOGBIS 
C
C=======================================================================
C     Expression for the function of the volume integral of the 
C     experimental density. Then, the integral is done via SIMPSON
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
          DENEXP_INTEGR=4*PINUMB*DENEXP(XARGUM)*XARGUM**2
      END IF
C
      IF (IARG_Z.EQ.1) THEN
          DENEXP_INTEGR=2*PINUMB*AOSCIL**3*DENEXP(XARGUM)*SQRT(XARGUM)
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
      FUNCTION DENTHE_INTEGR(XARGUM)
C
      COMMON
     *       /ARGUMN/ IARG_R,IARG_Z
      COMMON
     *       /SIMUSE/ AOSCIL,RADCUT,INUCLI
      COMMON
     *       /PRFILE/ ISCREN,LOGFIL,LOGCHP,LOGCHN,LOGHIP,LOGHIN,
     *                LOGPRO,LOGNEU,LOGWRI,LSCREN,IRESUL,ICONVE,
     *                NOUTPT,LOGAUX,LOGENE,LOGRAD,LOGGRD,LOGMSH,
     *                                                   LOGBIS 
C
C=======================================================================
C     Expression for the function of the volume integral of the 
C     theoretical density. Then, the integral is done via SIMPSON
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
          DENTHE_INTEGR=4*PINUMB*DENSIT(XARGUM)*XARGUM**2
      END IF
C
      IF (IARG_Z.EQ.1) THEN
          DENTHE_INTEGR=2*PINUMB*AOSCIL**3*DENSIT(XARGUM)*SQRT(XARGUM)
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
      FUNCTION DENTHE_DENEXP(XARGUM)
C
C=======================================================================
C     Expression for the function of the integral of the difference
C     experimental - theory density. The integral is done via SIMPSON
C=======================================================================
C
      DENTHE_DENEXP=DENSIT(XARGUM)-DENEXP(XARGUM)
C
C=======================================================================
C
      RETURN
      END
C      
C=======================================================================
C=======================================================================
C
      FUNCTION RMSEXP_INTEGR(XARGUM)
C
      COMMON
     *       /ARGUMN/ IARG_R,IARG_Z
      COMMON
     *       /SIMUSE/ AOSCIL,RADCUT,INUCLI
      COMMON
     *       /PRFILE/ ISCREN,LOGFIL,LOGCHP,LOGCHN,LOGHIP,LOGHIN,
     *                LOGPRO,LOGNEU,LOGWRI,LSCREN,IRESUL,ICONVE,
     *                NOUTPT,LOGAUX,LOGENE,LOGRAD,LOGGRD,LOGMSH,
     *                                                   LOGBIS 
C
C=======================================================================
C     Expression for the function of the integral of the volume 
C     integral of the rms radius. The integral is done via SIMPSON
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
          RMSEXP_INTEGR=4*PINUMB*DENEXP(XARGUM)*XARGUM**4
      END IF
C
      IF (IARG_Z.EQ.1) THEN
          RMSEXP_INTEGR=2*PINUMB*AOSCIL**3*DENEXP(XARGUM)
     *                 *SQRT(XARGUM)*XARGUM
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
      SUBROUTINE SIMPSN_INTEGR(A,B,N,F,SIMINT,ERRORS)
C
      EXTERNAL
     *         F 
C   
      COMMON
     *       /PRFILE/ ISCREN,LOGFIL,LOGCHP,LOGCHN,LOGHIP,LOGHIN,
     *                LOGPRO,LOGNEU,LOGWRI,LSCREN,IRESUL,ICONVE,
     *                NOUTPT,LOGAUX,LOGENE,LOGRAD,LOGGRD,LOGMSH,
     *                                                   LOGBIS 
C
C=======================================================================
C
C     Integration using the SIMPSON METHOD. The user have to give
C     the integration limits A and B,  the number of points N and 
C     the function to integrate F. For clarification the integral
C     is:
C                     SIMINT = \int F(X) dX
C
C     ERRORS - The nominal error of the Simpson method without f'
C
C=======================================================================
C
      STEP_H=(B-A)/N

      SIMINT=0.0

      DO I=0,N
         XARGUM=A+I*STEP_H
         FARGUM=F(XARGUM)
         IF ((I.EQ.0).OR.(I.EQ.N)) SIMINT=SIMINT+FARGUM
         IF ((I.GT.0).AND.(I.LT.N).AND.(MOD(I,2).EQ.0)) THEN
             SIMINT=SIMINT+2*FARGUM
         END IF
         IF ((I.GT.0).AND.(I.LT.N).AND.(MOD(I,2).NE.0)) THEN
             SIMINT=SIMINT+4*FARGUM
         END IF
      END DO
C
C=======================================================================
C
      SIMINT=SIMINT*STEP_H/3
C
      ERRORS=(B-A)*STEP_H**4/180.0
C
C=======================================================================
C
      IF (LOGWRI.GT.4) THEN
          WRITE(LOGFIL,'(12X,''Exiting  SIMPSN_INTEGR'')')
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
C                      TIMING SUBROUTINES
C
C=======================================================================
C=======================================================================
C
      SUBROUTINE TIMPRI
C
      PARAMETER
     *          (NDSUBR=30)
C
      PARAMETER
     *          (NDCOLU=2)
C
      CHARACTER
     *          OUITEM(NDCOLU)*35
C
      CHARACTER
     *          NAMALL*6
C
      COMMON
     *       /TIMPRI_TIMEXE/ TIMALL(1:NDSUBR),TIMACT(1:NDSUBR)
     *       /TIMEXE_KONOFF/ KONOFF(1:NDSUBR)
     *       /TIMEXE_NAMALL/ NAMALL(1:NDSUBR)
     *       /CPUDAT_TIMNOW/ NUSUBR
      COMMON
     *       /PRFILE/ ISCREN,LOGFIL,LOGCHP,LOGCHN,LOGHIP,LOGHIN,
     *                LOGPRO,LOGNEU,LOGWRI,LSCREN,IRESUL,ICONVE,
     *                NOUTPT,LOGAUX,LOGENE,LOGRAD,LOGGRD,LOGMSH,
     *                                                   LOGBIS    
C
C=======================================================================
C     This subroutine prints execution times for predefined subroutines
C=======================================================================
C
C      CALL CPUTIM('WSUNIV',0)
C
      IF (LOGWRI.GT.0) THEN
C
          WRITE(LOGFIL,'(/,''Running TIMPRI'',/,''NUSUBR,NDCOLU'',/,
     *                                2(3X,I4))') NUSUBR,NDCOLU
      END IF
C
C=======================================================================
C
      S_PERD=86400.0
      S_PERH=3600.0
      S_PMIN=60.0
      NSPERD=86400
      NSPERH=3600
      NSPMIN=60
C
      IF (NUSUBR.GT.0) THEN
C
          IF (LOGWRI.GT.0)
     *        WRITE(LOGFIL,'(//,80(''#''),/,         ''#'',78X,''#'',/,
     *                   ''#    CPU and system times for the routines'',
     *                   '' and their down-calling trees T1     #'',/,
     *                    ''#'',78X,''#'',/,               80(''#''),/,
     *                                               ''#'',78X,''#'')')
C_______________________________________________________________________
C
          DO ITSUBR=1,NUSUBR,NDCOLU
C
             DO NUITEM=1,NDCOLU
                 CALL BLANKS_STRING(OUITEM(NUITEM),LENGTH)
             END DO
C
             MXITEM=MIN(NUSUBR-ITSUBR+1,NDCOLU)
C
             DO NUITEM=1,MXITEM
C
                MUSUBR=ITSUBR+NUITEM-1
C
                IF (KONOFF(MUSUBR).EQ.1) THEN
C
                    IF (LOGWRI.GT.0)
     *                  WRITE(LOGFIL,'(/,1X,19(1H/),
     *                    '' In TIMPRI time is still on for '',
     *                                       A6,3X,19(1H/),/)')
     *                                           NAMALL(MUSUBR)
C
                    STOP 'In TIMPRI time is still on <<= '
C
                END IF
C
                N_DAYS=INT(TIMALL(MUSUBR)/S_PERD)
                T_DAYS=FLOAT(N_DAYS*NSPERD)
                N_HOUR=INT((TIMALL(MUSUBR)-T_DAYS)/S_PERH)
                T_HOUR=FLOAT(N_HOUR*NSPERH+N_DAYS*NSPERD)
                N_MINS=INT((TIMALL(MUSUBR)-T_HOUR)/S_PMIN)
                T_MINS=FLOAT(N_MINS*NSPMIN+N_HOUR*NSPERH+N_DAYS*NSPERD)
                SECOND=TIMALL(MUSUBR)-T_MINS
C
                WRITE(OUITEM(NUITEM),'(I2,'' d '',I2,'' h '',I2,
     *                                 '' min '',F6.3,'' s => '',A6)')
     *                N_DAYS,N_HOUR,N_MINS,SECOND,NAMALL(MUSUBR)
C
             END DO
C
             IF (LOGWRI.GT.0)
     *           WRITE(LOGFIL,'(''#'',2X,<NDCOLU>(A35,2X),2X,''#'',/,
     *                                             ''#'',78X,''#'')')
     *                               (OUITEM(NUITEM),NUITEM=1,NDCOLU)
C
          END DO
C
          IF (LOGWRI.GT.0) WRITE(LOGFIL,'(80(''#''),/)')
C_______________________________________________________________________
C
          IF (LOGWRI.GT.0)
     *        WRITE(LOGFIL,'(//,80(''#''),/,         ''#'',78X,''#'',/,
     *                   ''#    CPU and system times for the routines'',
     *                   '' and their down-calling trees T2     #'',/,
     *                    ''#'',78X,''#'',/,               80(''#''),/,
     *                                               ''#'',78X,''#'')')
C
          DO ITSUBR=1,NUSUBR,NDCOLU
C
             DO NUITEM=1,NDCOLU
                 CALL BLANKS_STRING(OUITEM(NUITEM),LENGTH)
             END DO
C
             MXITEM=MIN(NUSUBR-ITSUBR+1,NDCOLU)
C
             DO NUITEM=1,MXITEM
C
                MUSUBR=ITSUBR+NUITEM-1
C
                IF (KONOFF(MUSUBR).EQ.1) THEN
C
                    IF (LOGWRI.GT.0)
     *                  WRITE(LOGFIL,'(/,1X,19(1H/),
     *                    '' In TIMPRI time is still on for '',
     *                                       A6,3X,19(1H/),/)')
     *                                           NAMALL(MUSUBR)
C
                    STOP 'In TIMPRI time is still on <<= '
C
                END IF
C
                WRITE(OUITEM(NUITEM),'(F23.3,'' s => '',A6)')
     *                                 TIMALL(MUSUBR),NAMALL(MUSUBR)
C
             END DO
C
             IF (LOGWRI.GT.0)
     *           WRITE(LOGFIL,'(''#'',2X,<NDCOLU>(A35,2X),2X,''#'',/,
     *                                         ''#'',78X,''#'')')
     *                               (OUITEM(NUITEM),NUITEM=1,NDCOLU)
C
          END DO
C
          IF (LOGWRI.GT.0) WRITE(LOGFIL,'(80(''#''),/)')
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
      BLOCK DATA CPUDAT
C
      COMMON
     *       /CPUDAT_TIMNOW/ NUSUBR
C
      DATA
     *       NUSUBR / 0 /
C
      END
C
C=======================================================================
C=======================================================================
C
      SUBROUTINE CPUTIM(NAMSUB,IONOFF)
C
      PARAMETER
     *         (NDSUBR=30)
C
      CHARACTER
     *          NAMALL*6
C
      CHARACTER
     *          NAMSUB*6
C
      COMMON
     *       /TIMPRI_TIMEXE/ TIMALL(NDSUBR),TIMACT(NDSUBR)
     *       /TIMEXE_KONOFF/ KONOFF(NDSUBR)
     *       /TIMEXE_NAMALL/ NAMALL(NDSUBR)
C
      COMMON
     *       /CPUDAT_TIMNOW/ NUSUBR
C
      COMMON
     *       /PRFILE/ ISCREN,LOGFIL,LOGCHP,LOGCHN,LOGHIP,LOGHIN,
     *                LOGPRO,LOGNEU,LOGWRI,LSCREN,IRESUL,ICONVE,
     *                NOUTPT,LOGAUX,LOGENE,LOGRAD,LOGGRD,LOGMSH,
     *                                                   LOGBIS   
C
C=======================================================================
C     This subroutine collects execution times for different
C     subroutines
C=======================================================================
C
      I_SUBR=0
C
C=======================================================================
C     Each subroutine is identified by its name. Consecutive
C     numbers are attributed to subroutines  in the sequence
C     of consecutive calls to CPUTIM
C=======================================================================
C
      IF (NUSUBR.GT.0) THEN
C
          DO MUSUBR=1,NUSUBR
             IF (NAMALL(MUSUBR).EQ.NAMSUB) I_SUBR=MUSUBR
          END DO
C
      ELSE
C
          DO MUSUBR=1,NDSUBR
             KONOFF(MUSUBR)=0
             TIMACT(MUSUBR)=0.0
          END DO
C
      END IF
C
      IF (I_SUBR.EQ.0) THEN
C
          IF (NUSUBR.GE.NDSUBR) STOP 'Too large NUSUBR in CPUTIM <<= '
C
          NUSUBR=NUSUBR+1
          I_SUBR=NUSUBR
          NAMALL(I_SUBR)=NAMSUB
          TIMALL(I_SUBR)=0.0
C
      END IF
C
C=======================================================================
C
C     Here a call is made to CPU_TIME which returns the cpu time
C
C=======================================================================
C
      CALL CPU_TIME(SECNDS) 
C
C=======================================================================
C
      IF (IONOFF.NE.0.AND.IONOFF.NE.1)
     *    STOP 'Wrong IONOFF in CPUTIM <<= '
C
C     Here the time is set
C
      IF (IONOFF.EQ.1) THEN
C
          IF (KONOFF(I_SUBR).EQ.1) THEN
              IF (LOGWRI.GT.0)
     *            WRITE(LOGFIL,'(''Subroutine '',A6)') NAMSUB
              IF (ISCREN.GT.0)
     *            WRITE(LSCREN,'(''Subroutine '',A6)') NAMSUB
              STOP 'CPU_TIME is already on in CPUTIM <<= '
          END IF
C
          TIMACT(I_SUBR)=SECNDS
          KONOFF(I_SUBR)=1
C
      ELSE
C
C         Here the time is collected
C
          IF (KONOFF(I_SUBR).EQ.0) THEN
              IF (LOGWRI.GT.0)
     *            WRITE(LOGFIL,'(''Subroutine '',A6)') NAMSUB
              IF (ISCREN.GT.0)
     *            WRITE(LSCREN,'(''Subroutine '',A6)') NAMSUB
              STOP 'CPU_TIME is already off in CPUTIM <<= '
          END IF
C
          TIMALL(I_SUBR)=TIMALL(I_SUBR)+SECNDS-TIMACT(I_SUBR)
          KONOFF(I_SUBR)=0
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
      SUBROUTINE BLANKS_STRING(STRING,LENGTH)
C
      CHARACTER
     *           STRING*(*)
C
C=======================================================================
C
      LENGTH=LEN(STRING)
C
      DO I=1,LENGTH
         STRING(I:I)=' '
      END DO
C
C=======================================================================
C
      RETURN
      END 
C
C=======================================================================
C
C
C=======================================================================
C=======================================================================
C

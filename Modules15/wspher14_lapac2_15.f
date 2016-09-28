C FILE NAME = wspher14_lapac2_15.f ! Keep this symbol:    $ident@string$
C
C=======================================================================
C=======================================================================
C                SINGULAR-VALUE DECOMPOSITION PACKAGE
C=======================================================================
C=======================================================================
C
      SUBROUTINE DBDSQR( UPLO, N, NCVT, NRU, NCC, D, E, VT, LDVT, U,
     $                   LDU, C, LDC, WORK, INFO )
*
*  -- LAPACK routine (version 3.1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     January 2007
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDC, LDU, LDVT, N, NCC, NCVT, NRU
*     ..
*     .. Array Arguments ..
CID      DOUBLE PRECISION   C( LDC, * ), D( * ), E( * ), U( LDU, * ),
CID     $                   VT( LDVT, * ), WORK( * )
      DIMENSION   C( LDC, * ), D( * ), E( * ), U( LDU, * ),
     $                   VT( LDVT, * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DBDSQR computes the singular values and, optionally, the right and/or
*  left singular vectors from the singular value decomposition (SVD) of
*  a real N-by-N (upper or lower) bidiagonal matrix B using the implicit
*  zero-shift QR algorithm.  The SVD of B has the form
* 
*     B = Q * S * P**T
* 
*  where S is the diagonal matrix of singular values, Q is an orthogonal
*  matrix of left singular vectors, and P is an orthogonal matrix of
*  right singular vectors.  If left singular vectors are requested, this
*  subroutine actually returns U*Q instead of Q, and, if right singular
*  vectors are requested, this subroutine returns P**T*VT instead of
*  P**T, for given real input matrices U and VT.  When U and VT are the
*  orthogonal matrices that reduce a general matrix A to bidiagonal
*  form:  A = U*B*VT, as computed by DGEBRD, then
*
*     A = (U*Q) * S * (P**T*VT)
*
*  is the SVD of A.  Optionally, the subroutine may also compute Q**T*C
*  for a given real input matrix C.
*
*  See "Computing  Small Singular Values of Bidiagonal Matrices With
*  Guaranteed High Relative Accuracy," by J. Demmel and W. Kahan,
*  LAPACK Working Note #3 (or SIAM J. Sci. Statist. Comput. vol. 11,
*  no. 5, pp. 873-912, Sept 1990) and
*  "Accurate singular values and differential qd algorithms," by
*  B. Parlett and V. Fernando, Technical Report CPAM-554, Mathematics
*  Department, University of California at Berkeley, July 1992
*  for a detailed description of the algorithm.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  B is upper bidiagonal;
*          = 'L':  B is lower bidiagonal.
*
*  N       (input) INTEGER
*          The order of the matrix B.  N >= 0.
*
*  NCVT    (input) INTEGER
*          The number of columns of the matrix VT. NCVT >= 0.
*
*  NRU     (input) INTEGER
*          The number of rows of the matrix U. NRU >= 0.
*
*  NCC     (input) INTEGER
*          The number of columns of the matrix C. NCC >= 0.
*
*  D       (input/output) DOUBLE PRECISION array, dimension (N)
*          On entry, the n diagonal elements of the bidiagonal matrix B.
*          On exit, if INFO=0, the singular values of B in decreasing
*          order.
*
*  E       (input/output) DOUBLE PRECISION array, dimension (N-1)
*          On entry, the N-1 offdiagonal elements of the bidiagonal
*          matrix B. 
*          On exit, if INFO = 0, E is destroyed; if INFO > 0, D and E
*          will contain the diagonal and superdiagonal elements of a
*          bidiagonal matrix orthogonally equivalent to the one given
*          as input.
*
*  VT      (input/output) DOUBLE PRECISION array, dimension (LDVT, NCVT)
*          On entry, an N-by-NCVT matrix VT.
*          On exit, VT is overwritten by P**T * VT.
*          Not referenced if NCVT = 0.
*
*  LDVT    (input) INTEGER
*          The leading dimension of the array VT.
*          LDVT >= max(1,N) if NCVT > 0; LDVT >= 1 if NCVT = 0.
*
*  U       (input/output) DOUBLE PRECISION array, dimension (LDU, N)
*          On entry, an NRU-by-N matrix U.
*          On exit, U is overwritten by U * Q.
*          Not referenced if NRU = 0.
*
*  LDU     (input) INTEGER
*          The leading dimension of the array U.  LDU >= max(1,NRU).
*
*  C       (input/output) DOUBLE PRECISION array, dimension (LDC, NCC)
*          On entry, an N-by-NCC matrix C.
*          On exit, C is overwritten by Q**T * C.
*          Not referenced if NCC = 0.
*
*  LDC     (input) INTEGER
*          The leading dimension of the array C.
*          LDC >= max(1,N) if NCC > 0; LDC >=1 if NCC = 0.
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (2*N)
*          if NCVT = NRU = NCC = 0, (max(1, 4*N)) otherwise
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  If INFO = -i, the i-th argument had an illegal value
*          > 0:  the algorithm did not converge; D and E contain the
*                elements of a bidiagonal matrix which is orthogonally
*                similar to the input matrix B;  if INFO = i, i
*                elements of E have not converged to zero.
*
*  Internal Parameters
*  ===================
*
*  TOLMUL  DOUBLE PRECISION, default = max(10,min(100,EPS**(-1/8)))
*          TOLMUL controls the convergence criterion of the QR loop.
*          If it is positive, TOLMUL*EPS is the desired relative
*             precision in the computed singular values.
*          If it is negative, abs(TOLMUL*EPS*sigma_max) is the
*             desired absolute accuracy in the computed singular
*             values (corresponds to relative accuracy
*             abs(TOLMUL*EPS) in the largest singular value.
*          abs(TOLMUL) should be between 1 and 1/EPS, and preferably
*             between 10 (for fast convergence) and .1/EPS
*             (for there to be some accuracy in the results).
*          Default is to lose at either one eighth or 2 of the
*             available decimal digits in each computed singular value
*             (whichever is smaller).
*
*  MAXITR  INTEGER, default = 6
*          MAXITR controls the maximum number of passes of the
*          algorithm through its inner loop. The algorithms stops
*          (and so fails to converge) if the number of passes
*          through the inner loop exceeds MAXITR*N**2.
*
*  =====================================================================
*
*     .. Parameters ..
CID      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
CID      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
CID      DOUBLE PRECISION   NEGONE
      REAL                 NEGONE
      PARAMETER          ( NEGONE = -1.0D0 )
CID      DOUBLE PRECISION   HNDRTH
      PARAMETER          ( HNDRTH = 0.01D0 )
CID      DOUBLE PRECISION   TEN
      PARAMETER          ( TEN = 10.0D0 )
CID      DOUBLE PRECISION   HNDRD
      PARAMETER          ( HNDRD = 100.0D0 )
CID      DOUBLE PRECISION   MEIGTH
      REAL                 MEIGTH
      PARAMETER          ( MEIGTH = -0.125D0 )
      INTEGER            MAXITR
      PARAMETER          ( MAXITR = 6 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LOWER, ROTATE
      REAL               MU
      INTEGER            I, IDIR, ISUB, ITER, J, LL, LLL, M, MAXIT, NM1,
     $                   NM12, NM13, OLDLL, OLDM
CID      DOUBLE PRECISION   ABSE, ABSS, COSL, COSR, CS, EPS, F, G, H, MU,
CID     $                   OLDCS, OLDSN, R, SHIFT, SIGMN, SIGMX, SINL,
CID     $                   SINR, SLL, SMAX, SMIN, SMINL, SMINOA,
CID     $                   SN, THRESH, TOL, TOLMUL, UNFL
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
CID      DOUBLE PRECISION   DLAMCH
      EXTERNAL           LSAME, DLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLARTG, DLAS2, DLASQ1, DLASR, DLASV2, DROT,
     $                   DSCAL, DSWAP, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, MAX, MIN, SIGN, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      LOWER = LSAME( UPLO, 'L' )
      IF( .NOT.LSAME( UPLO, 'U' ) .AND. .NOT.LOWER ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NCVT.LT.0 ) THEN
         INFO = -3
      ELSE IF( NRU.LT.0 ) THEN
         INFO = -4
      ELSE IF( NCC.LT.0 ) THEN
         INFO = -5
      ELSE IF( ( NCVT.EQ.0 .AND. LDVT.LT.1 ) .OR.
     $         ( NCVT.GT.0 .AND. LDVT.LT.MAX( 1, N ) ) ) THEN
         INFO = -9
      ELSE IF( LDU.LT.MAX( 1, NRU ) ) THEN
         INFO = -11
      ELSE IF( ( NCC.EQ.0 .AND. LDC.LT.1 ) .OR.
     $         ( NCC.GT.0 .AND. LDC.LT.MAX( 1, N ) ) ) THEN
         INFO = -13
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DBDSQR', -INFO )
         RETURN
      END IF
      IF( N.EQ.0 )
     $   RETURN
      IF( N.EQ.1 )
     $   GO TO 160
*
*     ROTATE is true if any singular vectors desired, false otherwise
*
      ROTATE = ( NCVT.GT.0 ) .OR. ( NRU.GT.0 ) .OR. ( NCC.GT.0 )
*
*     If no singular vectors desired, use qd algorithm
*
      IF( .NOT.ROTATE ) THEN
         CALL DLASQ1( N, D, E, WORK, INFO )
         RETURN
      END IF
*
      NM1 = N - 1
      NM12 = NM1 + NM1
      NM13 = NM12 + NM1
      IDIR = 0
*
*     Get machine constants
*
      EPS = DLAMCH( 'Epsilon' )
      UNFL = DLAMCH( 'Safe minimum' )
*
*     If matrix lower bidiagonal, rotate to be upper bidiagonal
*     by applying Givens rotations on the left
*
      IF( LOWER ) THEN
         DO 10 I = 1, N - 1
            CALL DLARTG( D( I ), E( I ), CS, SN, R )
            D( I ) = R
            E( I ) = SN*D( I+1 )
            D( I+1 ) = CS*D( I+1 )
            WORK( I ) = CS
            WORK( NM1+I ) = SN
   10    CONTINUE
*
*        Update singular vectors if desired
*
         IF( NRU.GT.0 )
     $      CALL DLASR( 'R', 'V', 'F', NRU, N, WORK( 1 ), WORK( N ), U,
     $                  LDU )
         IF( NCC.GT.0 )
     $      CALL DLASR( 'L', 'V', 'F', N, NCC, WORK( 1 ), WORK( N ), C,
     $                  LDC )
      END IF
*
*     Compute singular values to relative accuracy TOL
*     (By setting TOL to be negative, algorithm will compute
*     singular values to absolute accuracy ABS(TOL)*norm(input matrix))
*
      TOLMUL = MAX( TEN, MIN( HNDRD, EPS**MEIGTH ) )
      TOL = TOLMUL*EPS
*
*     Compute approximate maximum, minimum singular values
*
      SMAX = ZERO
      DO 20 I = 1, N
         SMAX = MAX( SMAX, ABS( D( I ) ) )
   20 CONTINUE
      DO 30 I = 1, N - 1
         SMAX = MAX( SMAX, ABS( E( I ) ) )
   30 CONTINUE
      SMINL = ZERO
      IF( TOL.GE.ZERO ) THEN
*
*        Relative accuracy desired
*
         SMINOA = ABS( D( 1 ) )
         IF( SMINOA.EQ.ZERO )
     $      GO TO 50
         MU = SMINOA
         DO 40 I = 2, N
            MU = ABS( D( I ) )*( MU / ( MU+ABS( E( I-1 ) ) ) )
            SMINOA = MIN( SMINOA, MU )
            IF( SMINOA.EQ.ZERO )
     $         GO TO 50
   40    CONTINUE
   50    CONTINUE
         SMINOA = SMINOA / SQRT( DBLE( N ) )
         THRESH = MAX( TOL*SMINOA, MAXITR*N*N*UNFL )
      ELSE
*
*        Absolute accuracy desired
*
         THRESH = MAX( ABS( TOL )*SMAX, MAXITR*N*N*UNFL )
      END IF
*
*     Prepare for main iteration loop for the singular values
*     (MAXIT is the maximum number of passes through the inner
*     loop permitted before nonconvergence signalled.)
*
      MAXIT = MAXITR*N*N
      ITER = 0
      OLDLL = -1
      OLDM = -1
*
*     M points to last element of unconverged part of matrix
*
      M = N
*
*     Begin main iteration loop
*
   60 CONTINUE
*
*     Check for convergence or exceeding iteration count
*
      IF( M.LE.1 )
     $   GO TO 160
      IF( ITER.GT.MAXIT )
     $   GO TO 200
*
*     Find diagonal block of matrix to work on
*
      IF( TOL.LT.ZERO .AND. ABS( D( M ) ).LE.THRESH )
     $   D( M ) = ZERO
      SMAX = ABS( D( M ) )
      SMIN = SMAX
      DO 70 LLL = 1, M - 1
         LL = M - LLL
         ABSS = ABS( D( LL ) )
         ABSE = ABS( E( LL ) )
         IF( TOL.LT.ZERO .AND. ABSS.LE.THRESH )
     $      D( LL ) = ZERO
         IF( ABSE.LE.THRESH )
     $      GO TO 80
         SMIN = MIN( SMIN, ABSS )
         SMAX = MAX( SMAX, ABSS, ABSE )
   70 CONTINUE
      LL = 0
      GO TO 90
   80 CONTINUE
      E( LL ) = ZERO
*
*     Matrix splits since E(LL) = 0
*
      IF( LL.EQ.M-1 ) THEN
*
*        Convergence of bottom singular value, return to top of loop
*
         M = M - 1
         GO TO 60
      END IF
   90 CONTINUE
      LL = LL + 1
*
*     E(LL) through E(M-1) are nonzero, E(LL-1) is zero
*
      IF( LL.EQ.M-1 ) THEN
*
*        2 by 2 block, handle separately
*
         CALL DLASV2( D( M-1 ), E( M-1 ), D( M ), SIGMN, SIGMX, SINR,
     $                COSR, SINL, COSL )
         D( M-1 ) = SIGMX
         E( M-1 ) = ZERO
         D( M ) = SIGMN
*
*        Compute singular vectors, if desired
*
         IF( NCVT.GT.0 )
     $      CALL DROT( NCVT, VT( M-1, 1 ), LDVT, VT( M, 1 ), LDVT, COSR,
     $                 SINR )
         IF( NRU.GT.0 )
     $      CALL DROT( NRU, U( 1, M-1 ), 1, U( 1, M ), 1, COSL, SINL )
         IF( NCC.GT.0 )
     $      CALL DROT( NCC, C( M-1, 1 ), LDC, C( M, 1 ), LDC, COSL,
     $                 SINL )
         M = M - 2
         GO TO 60
      END IF
*
*     If working on new submatrix, choose shift direction
*     (from larger end diagonal element towards smaller)
*
      IF( LL.GT.OLDM .OR. M.LT.OLDLL ) THEN
         IF( ABS( D( LL ) ).GE.ABS( D( M ) ) ) THEN
*
*           Chase bulge from top (big end) to bottom (small end)
*
            IDIR = 1
         ELSE
*
*           Chase bulge from bottom (big end) to top (small end)
*
            IDIR = 2
         END IF
      END IF
*
*     Apply convergence tests
*
      IF( IDIR.EQ.1 ) THEN
*
*        Run convergence test in forward direction
*        First apply standard test to bottom of matrix
*
         IF( ABS( E( M-1 ) ).LE.ABS( TOL )*ABS( D( M ) ) .OR.
     $       ( TOL.LT.ZERO .AND. ABS( E( M-1 ) ).LE.THRESH ) ) THEN
            E( M-1 ) = ZERO
            GO TO 60
         END IF
*
         IF( TOL.GE.ZERO ) THEN
*
*           If relative accuracy desired,
*           apply convergence criterion forward
*
            MU = ABS( D( LL ) )
            SMINL = MU
            DO 100 LLL = LL, M - 1
               IF( ABS( E( LLL ) ).LE.TOL*MU ) THEN
                  E( LLL ) = ZERO
                  GO TO 60
               END IF
               MU = ABS( D( LLL+1 ) )*( MU / ( MU+ABS( E( LLL ) ) ) )
               SMINL = MIN( SMINL, MU )
  100       CONTINUE
         END IF
*
      ELSE
*
*        Run convergence test in backward direction
*        First apply standard test to top of matrix
*
         IF( ABS( E( LL ) ).LE.ABS( TOL )*ABS( D( LL ) ) .OR.
     $       ( TOL.LT.ZERO .AND. ABS( E( LL ) ).LE.THRESH ) ) THEN
            E( LL ) = ZERO
            GO TO 60
         END IF
*
         IF( TOL.GE.ZERO ) THEN
*
*           If relative accuracy desired,
*           apply convergence criterion backward
*
            MU = ABS( D( M ) )
            SMINL = MU
            DO 110 LLL = M - 1, LL, -1
               IF( ABS( E( LLL ) ).LE.TOL*MU ) THEN
                  E( LLL ) = ZERO
                  GO TO 60
               END IF
               MU = ABS( D( LLL ) )*( MU / ( MU+ABS( E( LLL ) ) ) )
               SMINL = MIN( SMINL, MU )
  110       CONTINUE
         END IF
      END IF
      OLDLL = LL
      OLDM = M
*
*     Compute shift.  First, test if shifting would ruin relative
*     accuracy, and if so set the shift to zero.
*
      IF( TOL.GE.ZERO .AND. N*TOL*( SMINL / SMAX ).LE.
     $    MAX( EPS, HNDRTH*TOL ) ) THEN
*
*        Use a zero shift to avoid loss of relative accuracy
*
         SHIFT = ZERO
      ELSE
*
*        Compute the shift from 2-by-2 block at end of matrix
*
         IF( IDIR.EQ.1 ) THEN
            SLL = ABS( D( LL ) )
            CALL DLAS2( D( M-1 ), E( M-1 ), D( M ), SHIFT, R )
         ELSE
            SLL = ABS( D( M ) )
            CALL DLAS2( D( LL ), E( LL ), D( LL+1 ), SHIFT, R )
         END IF
*
*        Test if shift negligible, and if so set to zero
*
         IF( SLL.GT.ZERO ) THEN
            IF( ( SHIFT / SLL )**2.LT.EPS )
     $         SHIFT = ZERO
         END IF
      END IF
*
*     Increment iteration count
*
      ITER = ITER + M - LL
*
*     If SHIFT = 0, do simplified QR iteration
*
      IF( SHIFT.EQ.ZERO ) THEN
         IF( IDIR.EQ.1 ) THEN
*
*           Chase bulge from top to bottom
*           Save cosines and sines for later singular vector updates
*
            CS = ONE
            OLDCS = ONE
            DO 120 I = LL, M - 1
               CALL DLARTG( D( I )*CS, E( I ), CS, SN, R )
               IF( I.GT.LL )
     $            E( I-1 ) = OLDSN*R
               CALL DLARTG( OLDCS*R, D( I+1 )*SN, OLDCS, OLDSN, D( I ) )
               WORK( I-LL+1 ) = CS
               WORK( I-LL+1+NM1 ) = SN
               WORK( I-LL+1+NM12 ) = OLDCS
               WORK( I-LL+1+NM13 ) = OLDSN
  120       CONTINUE
            H = D( M )*CS
            D( M ) = H*OLDCS
            E( M-1 ) = H*OLDSN
*
*           Update singular vectors
*
            IF( NCVT.GT.0 )
     $         CALL DLASR( 'L', 'V', 'F', M-LL+1, NCVT, WORK( 1 ),
     $                     WORK( N ), VT( LL, 1 ), LDVT )
            IF( NRU.GT.0 )
     $         CALL DLASR( 'R', 'V', 'F', NRU, M-LL+1, WORK( NM12+1 ),
     $                     WORK( NM13+1 ), U( 1, LL ), LDU )
            IF( NCC.GT.0 )
     $         CALL DLASR( 'L', 'V', 'F', M-LL+1, NCC, WORK( NM12+1 ),
     $                     WORK( NM13+1 ), C( LL, 1 ), LDC )
*
*           Test convergence
*
            IF( ABS( E( M-1 ) ).LE.THRESH )
     $         E( M-1 ) = ZERO
*
         ELSE
*
*           Chase bulge from bottom to top
*           Save cosines and sines for later singular vector updates
*
            CS = ONE
            OLDCS = ONE
            DO 130 I = M, LL + 1, -1
               CALL DLARTG( D( I )*CS, E( I-1 ), CS, SN, R )
               IF( I.LT.M )
     $            E( I ) = OLDSN*R
               CALL DLARTG( OLDCS*R, D( I-1 )*SN, OLDCS, OLDSN, D( I ) )
               WORK( I-LL ) = CS
               WORK( I-LL+NM1 ) = -SN
               WORK( I-LL+NM12 ) = OLDCS
               WORK( I-LL+NM13 ) = -OLDSN
  130       CONTINUE
            H = D( LL )*CS
            D( LL ) = H*OLDCS
            E( LL ) = H*OLDSN
*
*           Update singular vectors
*
            IF( NCVT.GT.0 )
     $         CALL DLASR( 'L', 'V', 'B', M-LL+1, NCVT, WORK( NM12+1 ),
     $                     WORK( NM13+1 ), VT( LL, 1 ), LDVT )
            IF( NRU.GT.0 )
     $         CALL DLASR( 'R', 'V', 'B', NRU, M-LL+1, WORK( 1 ),
     $                     WORK( N ), U( 1, LL ), LDU )
            IF( NCC.GT.0 )
     $         CALL DLASR( 'L', 'V', 'B', M-LL+1, NCC, WORK( 1 ),
     $                     WORK( N ), C( LL, 1 ), LDC )
*
*           Test convergence
*
            IF( ABS( E( LL ) ).LE.THRESH )
     $         E( LL ) = ZERO
         END IF
      ELSE
*
*        Use nonzero shift
*
         IF( IDIR.EQ.1 ) THEN
*
*           Chase bulge from top to bottom
*           Save cosines and sines for later singular vector updates
*
            F = ( ABS( D( LL ) )-SHIFT )*
     $          ( SIGN( ONE, D( LL ) )+SHIFT / D( LL ) )
            G = E( LL )
            DO 140 I = LL, M - 1
               CALL DLARTG( F, G, COSR, SINR, R )
               IF( I.GT.LL )
     $            E( I-1 ) = R
               F = COSR*D( I ) + SINR*E( I )
               E( I ) = COSR*E( I ) - SINR*D( I )
               G = SINR*D( I+1 )
               D( I+1 ) = COSR*D( I+1 )
               CALL DLARTG( F, G, COSL, SINL, R )
               D( I ) = R
               F = COSL*E( I ) + SINL*D( I+1 )
               D( I+1 ) = COSL*D( I+1 ) - SINL*E( I )
               IF( I.LT.M-1 ) THEN
                  G = SINL*E( I+1 )
                  E( I+1 ) = COSL*E( I+1 )
               END IF
               WORK( I-LL+1 ) = COSR
               WORK( I-LL+1+NM1 ) = SINR
               WORK( I-LL+1+NM12 ) = COSL
               WORK( I-LL+1+NM13 ) = SINL
  140       CONTINUE
            E( M-1 ) = F
*
*           Update singular vectors
*
            IF( NCVT.GT.0 )
     $         CALL DLASR( 'L', 'V', 'F', M-LL+1, NCVT, WORK( 1 ),
     $                     WORK( N ), VT( LL, 1 ), LDVT )
            IF( NRU.GT.0 )
     $         CALL DLASR( 'R', 'V', 'F', NRU, M-LL+1, WORK( NM12+1 ),
     $                     WORK( NM13+1 ), U( 1, LL ), LDU )
            IF( NCC.GT.0 )
     $         CALL DLASR( 'L', 'V', 'F', M-LL+1, NCC, WORK( NM12+1 ),
     $                     WORK( NM13+1 ), C( LL, 1 ), LDC )
*
*           Test convergence
*
            IF( ABS( E( M-1 ) ).LE.THRESH )
     $         E( M-1 ) = ZERO
*
         ELSE
*
*           Chase bulge from bottom to top
*           Save cosines and sines for later singular vector updates
*
            F = ( ABS( D( M ) )-SHIFT )*( SIGN( ONE, D( M ) )+SHIFT /
     $          D( M ) )
            G = E( M-1 )
            DO 150 I = M, LL + 1, -1
               CALL DLARTG( F, G, COSR, SINR, R )
               IF( I.LT.M )
     $            E( I ) = R
               F = COSR*D( I ) + SINR*E( I-1 )
               E( I-1 ) = COSR*E( I-1 ) - SINR*D( I )
               G = SINR*D( I-1 )
               D( I-1 ) = COSR*D( I-1 )
               CALL DLARTG( F, G, COSL, SINL, R )
               D( I ) = R
               F = COSL*E( I-1 ) + SINL*D( I-1 )
               D( I-1 ) = COSL*D( I-1 ) - SINL*E( I-1 )
               IF( I.GT.LL+1 ) THEN
                  G = SINL*E( I-2 )
                  E( I-2 ) = COSL*E( I-2 )
               END IF
               WORK( I-LL ) = COSR
               WORK( I-LL+NM1 ) = -SINR
               WORK( I-LL+NM12 ) = COSL
               WORK( I-LL+NM13 ) = -SINL
  150       CONTINUE
            E( LL ) = F
*
*           Test convergence
*
            IF( ABS( E( LL ) ).LE.THRESH )
     $         E( LL ) = ZERO
*
*           Update singular vectors if desired
*
            IF( NCVT.GT.0 )
     $         CALL DLASR( 'L', 'V', 'B', M-LL+1, NCVT, WORK( NM12+1 ),
     $                     WORK( NM13+1 ), VT( LL, 1 ), LDVT )
            IF( NRU.GT.0 )
     $         CALL DLASR( 'R', 'V', 'B', NRU, M-LL+1, WORK( 1 ),
     $                     WORK( N ), U( 1, LL ), LDU )
            IF( NCC.GT.0 )
     $         CALL DLASR( 'L', 'V', 'B', M-LL+1, NCC, WORK( 1 ),
     $                     WORK( N ), C( LL, 1 ), LDC )
         END IF
      END IF
*
*     QR iteration finished, go back and check convergence
*
      GO TO 60
*
*     All singular values converged, so make them positive
*
  160 CONTINUE
      DO 170 I = 1, N
         IF( D( I ).LT.ZERO ) THEN
            D( I ) = -D( I )
*
*           Change sign of singular vectors, if desired
*
            IF( NCVT.GT.0 )
     $         CALL DSCAL( NCVT, NEGONE, VT( I, 1 ), LDVT )
         END IF
  170 CONTINUE
*
*     Sort the singular values into decreasing order (insertion sort on
*     singular values, but only one transposition per singular vector)
*
      DO 190 I = 1, N - 1
*
*        Scan for smallest D(I)
*
         ISUB = 1
         SMIN = D( 1 )
         DO 180 J = 2, N + 1 - I
            IF( D( J ).LE.SMIN ) THEN
               ISUB = J
               SMIN = D( J )
            END IF
  180    CONTINUE
         IF( ISUB.NE.N+1-I ) THEN
*
*           Swap singular values and vectors
*
            D( ISUB ) = D( N+1-I )
            D( N+1-I ) = SMIN
            IF( NCVT.GT.0 )
     $         CALL DSWAP( NCVT, VT( ISUB, 1 ), LDVT, VT( N+1-I, 1 ),
     $                     LDVT )
            IF( NRU.GT.0 )
     $         CALL DSWAP( NRU, U( 1, ISUB ), 1, U( 1, N+1-I ), 1 )
            IF( NCC.GT.0 )
     $         CALL DSWAP( NCC, C( ISUB, 1 ), LDC, C( N+1-I, 1 ), LDC )
         END IF
  190 CONTINUE
      GO TO 220
*
*     Maximum number of iterations exceeded, failure to converge
*
  200 CONTINUE
      INFO = 0
      DO 210 I = 1, N - 1
         IF( E( I ).NE.ZERO )
     $      INFO = INFO + 1
  210 CONTINUE
  220 CONTINUE
      RETURN
*
*     End of DBDSQR
*
      END
      SUBROUTINE DGEBD2( M, N, A, LDA, D, E, TAUQ, TAUP, WORK, INFO )
*
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
*     ..
*     .. Array Arguments ..
CID      DOUBLE PRECISION   A( LDA, * ), D( * ), E( * ), TAUP( * ),
CID     $                   TAUQ( * ), WORK( * )
      DIMENSION          A( LDA, * ), D( * ), E( * ), TAUP( * ),
     $                   TAUQ( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DGEBD2 reduces a real general m by n matrix A to upper or lower
*  bidiagonal form B by an orthogonal transformation: Q' * A * P = B.
*
*  If m >= n, B is upper bidiagonal; if m < n, B is lower bidiagonal.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows in the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns in the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the m by n general matrix to be reduced.
*          On exit,
*          if m >= n, the diagonal and the first superdiagonal are
*            overwritten with the upper bidiagonal matrix B; the
*            elements below the diagonal, with the array TAUQ, represent
*            the orthogonal matrix Q as a product of elementary
*            reflectors, and the elements above the first superdiagonal,
*            with the array TAUP, represent the orthogonal matrix P as
*            a product of elementary reflectors;
*          if m < n, the diagonal and the first subdiagonal are
*            overwritten with the lower bidiagonal matrix B; the
*            elements below the first subdiagonal, with the array TAUQ,
*            represent the orthogonal matrix Q as a product of
*            elementary reflectors, and the elements above the diagonal,
*            with the array TAUP, represent the orthogonal matrix P as
*            a product of elementary reflectors.
*          See Further Details.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  D       (output) DOUBLE PRECISION array, dimension (min(M,N))
*          The diagonal elements of the bidiagonal matrix B:
*          D(i) = A(i,i).
*
*  E       (output) DOUBLE PRECISION array, dimension (min(M,N)-1)
*          The off-diagonal elements of the bidiagonal matrix B:
*          if m >= n, E(i) = A(i,i+1) for i = 1,2,...,n-1;
*          if m < n, E(i) = A(i+1,i) for i = 1,2,...,m-1.
*
*  TAUQ    (output) DOUBLE PRECISION array dimension (min(M,N))
*          The scalar factors of the elementary reflectors which
*          represent the orthogonal matrix Q. See Further Details.
*
*  TAUP    (output) DOUBLE PRECISION array, dimension (min(M,N))
*          The scalar factors of the elementary reflectors which
*          represent the orthogonal matrix P. See Further Details.
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (max(M,N))
*
*  INFO    (output) INTEGER
*          = 0: successful exit.
*          < 0: if INFO = -i, the i-th argument had an illegal value.
*
*  Further Details
*  ===============
*
*  The matrices Q and P are represented as products of elementary
*  reflectors:
*
*  If m >= n,
*
*     Q = H(1) H(2) . . . H(n)  and  P = G(1) G(2) . . . G(n-1)
*
*  Each H(i) and G(i) has the form:
*
*     H(i) = I - tauq * v * v'  and G(i) = I - taup * u * u'
*
*  where tauq and taup are real scalars, and v and u are real vectors;
*  v(1:i-1) = 0, v(i) = 1, and v(i+1:m) is stored on exit in A(i+1:m,i);
*  u(1:i) = 0, u(i+1) = 1, and u(i+2:n) is stored on exit in A(i,i+2:n);
*  tauq is stored in TAUQ(i) and taup in TAUP(i).
*
*  If m < n,
*
*     Q = H(1) H(2) . . . H(m-1)  and  P = G(1) G(2) . . . G(m)
*
*  Each H(i) and G(i) has the form:
*
*     H(i) = I - tauq * v * v'  and G(i) = I - taup * u * u'
*
*  where tauq and taup are real scalars, and v and u are real vectors;
*  v(1:i) = 0, v(i+1) = 1, and v(i+2:m) is stored on exit in A(i+2:m,i);
*  u(1:i-1) = 0, u(i) = 1, and u(i+1:n) is stored on exit in A(i,i+1:n);
*  tauq is stored in TAUQ(i) and taup in TAUP(i).
*
*  The contents of A on exit are illustrated by the following examples:
*
*  m = 6 and n = 5 (m > n):          m = 5 and n = 6 (m < n):
*
*    (  d   e   u1  u1  u1 )           (  d   u1  u1  u1  u1  u1 )
*    (  v1  d   e   u2  u2 )           (  e   d   u2  u2  u2  u2 )
*    (  v1  v2  d   e   u3 )           (  v1  e   d   u3  u3  u3 )
*    (  v1  v2  v3  d   e  )           (  v1  v2  e   d   u4  u4 )
*    (  v1  v2  v3  v4  d  )           (  v1  v2  v3  e   d   u5 )
*    (  v1  v2  v3  v4  v5 )
*
*  where d and e denote diagonal and off-diagonal elements of B, vi
*  denotes an element of the vector defining H(i), and ui an element of
*  the vector defining G(i).
*
*  =====================================================================
*
*     .. Parameters ..
CID      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLARF, DLARFG, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.LT.0 ) THEN
         CALL XERBLA( 'DGEBD2', -INFO )
         RETURN
      END IF
*
      IF( M.GE.N ) THEN
*
*        Reduce to upper bidiagonal form
*
         DO 10 I = 1, N
*
*           Generate elementary reflector H(i) to annihilate A(i+1:m,i)
*
            CALL DLARFG( M-I+1, A( I, I ), A( MIN( I+1, M ), I ), 1,
     $                   TAUQ( I ) )
            D( I ) = A( I, I )
            A( I, I ) = ONE
*
*           Apply H(i) to A(i:m,i+1:n) from the left
*
            IF( I.LT.N )
     $         CALL DLARF( 'Left', M-I+1, N-I, A( I, I ), 1, TAUQ( I ),
     $                     A( I, I+1 ), LDA, WORK )
            A( I, I ) = D( I )
*
            IF( I.LT.N ) THEN
*
*              Generate elementary reflector G(i) to annihilate
*              A(i,i+2:n)
*
               CALL DLARFG( N-I, A( I, I+1 ), A( I, MIN( I+2, N ) ),
     $                      LDA, TAUP( I ) )
               E( I ) = A( I, I+1 )
               A( I, I+1 ) = ONE
*
*              Apply G(i) to A(i+1:m,i+1:n) from the right
*
               CALL DLARF( 'Right', M-I, N-I, A( I, I+1 ), LDA,
     $                     TAUP( I ), A( I+1, I+1 ), LDA, WORK )
               A( I, I+1 ) = E( I )
            ELSE
               TAUP( I ) = ZERO
            END IF
   10    CONTINUE
      ELSE
*
*        Reduce to lower bidiagonal form
*
         DO 20 I = 1, M
*
*           Generate elementary reflector G(i) to annihilate A(i,i+1:n)
*
            CALL DLARFG( N-I+1, A( I, I ), A( I, MIN( I+1, N ) ), LDA,
     $                   TAUP( I ) )
            D( I ) = A( I, I )
            A( I, I ) = ONE
*
*           Apply G(i) to A(i+1:m,i:n) from the right
*
            IF( I.LT.M )
     $         CALL DLARF( 'Right', M-I, N-I+1, A( I, I ), LDA,
     $                     TAUP( I ), A( I+1, I ), LDA, WORK )
            A( I, I ) = D( I )
*
            IF( I.LT.M ) THEN
*
*              Generate elementary reflector H(i) to annihilate
*              A(i+2:m,i)
*
               CALL DLARFG( M-I, A( I+1, I ), A( MIN( I+2, M ), I ), 1,
     $                      TAUQ( I ) )
               E( I ) = A( I+1, I )
               A( I+1, I ) = ONE
*
*              Apply H(i) to A(i+1:m,i+1:n) from the left
*
               CALL DLARF( 'Left', M-I, N-I, A( I+1, I ), 1, TAUQ( I ),
     $                     A( I+1, I+1 ), LDA, WORK )
               A( I+1, I ) = E( I )
            ELSE
               TAUQ( I ) = ZERO
            END IF
   20    CONTINUE
      END IF
      RETURN
*
*     End of DGEBD2
*
      END
      SUBROUTINE DGEBRD( M, N, A, LDA, D, E, TAUQ, TAUP, WORK, LWORK,
     $                   INFO )
*
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LWORK, M, N
*     ..
*     .. Array Arguments ..
CID      DOUBLE PRECISION   A( LDA, * ), D( * ), E( * ), TAUP( * ),
CID     $                   TAUQ( * ), WORK( * )
      DIMENSION          A( LDA, * ), D( * ), E( * ), TAUP( * ),
     $                   TAUQ( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DGEBRD reduces a general real M-by-N matrix A to upper or lower
*  bidiagonal form B by an orthogonal transformation: Q**T * A * P = B.
*
*  If m >= n, B is upper bidiagonal; if m < n, B is lower bidiagonal.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows in the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns in the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the M-by-N general matrix to be reduced.
*          On exit,
*          if m >= n, the diagonal and the first superdiagonal are
*            overwritten with the upper bidiagonal matrix B; the
*            elements below the diagonal, with the array TAUQ, represent
*            the orthogonal matrix Q as a product of elementary
*            reflectors, and the elements above the first superdiagonal,
*            with the array TAUP, represent the orthogonal matrix P as
*            a product of elementary reflectors;
*          if m < n, the diagonal and the first subdiagonal are
*            overwritten with the lower bidiagonal matrix B; the
*            elements below the first subdiagonal, with the array TAUQ,
*            represent the orthogonal matrix Q as a product of
*            elementary reflectors, and the elements above the diagonal,
*            with the array TAUP, represent the orthogonal matrix P as
*            a product of elementary reflectors.
*          See Further Details.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  D       (output) DOUBLE PRECISION array, dimension (min(M,N))
*          The diagonal elements of the bidiagonal matrix B:
*          D(i) = A(i,i).
*
*  E       (output) DOUBLE PRECISION array, dimension (min(M,N)-1)
*          The off-diagonal elements of the bidiagonal matrix B:
*          if m >= n, E(i) = A(i,i+1) for i = 1,2,...,n-1;
*          if m < n, E(i) = A(i+1,i) for i = 1,2,...,m-1.
*
*  TAUQ    (output) DOUBLE PRECISION array dimension (min(M,N))
*          The scalar factors of the elementary reflectors which
*          represent the orthogonal matrix Q. See Further Details.
*
*  TAUP    (output) DOUBLE PRECISION array, dimension (min(M,N))
*          The scalar factors of the elementary reflectors which
*          represent the orthogonal matrix P. See Further Details.
*
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The length of the array WORK.  LWORK >= max(1,M,N).
*          For optimum performance LWORK >= (M+N)*NB, where NB
*          is the optimal blocksize.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*
*  Further Details
*  ===============
*
*  The matrices Q and P are represented as products of elementary
*  reflectors:
*
*  If m >= n,
*
*     Q = H(1) H(2) . . . H(n)  and  P = G(1) G(2) . . . G(n-1)
*
*  Each H(i) and G(i) has the form:
*
*     H(i) = I - tauq * v * v'  and G(i) = I - taup * u * u'
*
*  where tauq and taup are real scalars, and v and u are real vectors;
*  v(1:i-1) = 0, v(i) = 1, and v(i+1:m) is stored on exit in A(i+1:m,i);
*  u(1:i) = 0, u(i+1) = 1, and u(i+2:n) is stored on exit in A(i,i+2:n);
*  tauq is stored in TAUQ(i) and taup in TAUP(i).
*
*  If m < n,
*
*     Q = H(1) H(2) . . . H(m-1)  and  P = G(1) G(2) . . . G(m)
*
*  Each H(i) and G(i) has the form:
*
*     H(i) = I - tauq * v * v'  and G(i) = I - taup * u * u'
*
*  where tauq and taup are real scalars, and v and u are real vectors;
*  v(1:i) = 0, v(i+1) = 1, and v(i+2:m) is stored on exit in A(i+2:m,i);
*  u(1:i-1) = 0, u(i) = 1, and u(i+1:n) is stored on exit in A(i,i+1:n);
*  tauq is stored in TAUQ(i) and taup in TAUP(i).
*
*  The contents of A on exit are illustrated by the following examples:
*
*  m = 6 and n = 5 (m > n):          m = 5 and n = 6 (m < n):
*
*    (  d   e   u1  u1  u1 )           (  d   u1  u1  u1  u1  u1 )
*    (  v1  d   e   u2  u2 )           (  e   d   u2  u2  u2  u2 )
*    (  v1  v2  d   e   u3 )           (  v1  e   d   u3  u3  u3 )
*    (  v1  v2  v3  d   e  )           (  v1  v2  e   d   u4  u4 )
*    (  v1  v2  v3  v4  d  )           (  v1  v2  v3  e   d   u5 )
*    (  v1  v2  v3  v4  v5 )
*
*  where d and e denote diagonal and off-diagonal elements of B, vi
*  denotes an element of the vector defining H(i), and ui an element of
*  the vector defining G(i).
*
*  =====================================================================
*
*     .. Parameters ..
CID      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IINFO, J, LDWRKX, LDWRKY, LWKOPT, MINMN, NB,
     $                   NBMIN, NX
CID      DOUBLE PRECISION   WS
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEBD2, DGEMM, DLABRD, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MAX, MIN
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      INFO = 0
      NB = MAX( 1, ILAENV( 1, 'DGEBRD', ' ', M, N, -1, -1 ) )
      LWKOPT = ( M+N )*NB
      WORK( 1 ) = DBLE( LWKOPT )
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      ELSE IF( LWORK.LT.MAX( 1, M, N ) .AND. .NOT.LQUERY ) THEN
         INFO = -10
      END IF
      IF( INFO.LT.0 ) THEN
         CALL XERBLA( 'DGEBRD', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      MINMN = MIN( M, N )
      IF( MINMN.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
      WS = MAX( M, N )
      LDWRKX = M
      LDWRKY = N
*
      IF( NB.GT.1 .AND. NB.LT.MINMN ) THEN
*
*        Set the crossover point NX.
*
         NX = MAX( NB, ILAENV( 3, 'DGEBRD', ' ', M, N, -1, -1 ) )
*
*        Determine when to switch from blocked to unblocked code.
*
         IF( NX.LT.MINMN ) THEN
            WS = ( M+N )*NB
            IF( LWORK.LT.WS ) THEN
*
*              Not enough work space for the optimal NB, consider using
*              a smaller block size.
*
               NBMIN = ILAENV( 2, 'DGEBRD', ' ', M, N, -1, -1 )
               IF( LWORK.GE.( M+N )*NBMIN ) THEN
                  NB = LWORK / ( M+N )
               ELSE
                  NB = 1
                  NX = MINMN
               END IF
            END IF
         END IF
      ELSE
         NX = MINMN
      END IF
*
      DO 30 I = 1, MINMN - NX, NB
*
*        Reduce rows and columns i:i+nb-1 to bidiagonal form and return
*        the matrices X and Y which are needed to update the unreduced
*        part of the matrix
*
         CALL DLABRD( M-I+1, N-I+1, NB, A( I, I ), LDA, D( I ), E( I ),
     $                TAUQ( I ), TAUP( I ), WORK, LDWRKX,
     $                WORK( LDWRKX*NB+1 ), LDWRKY )
*
*        Update the trailing submatrix A(i+nb:m,i+nb:n), using an update
*        of the form  A := A - V*Y' - X*U'
*
         CALL DGEMM( 'No transpose', 'Transpose', M-I-NB+1, N-I-NB+1,
     $               NB, -ONE, A( I+NB, I ), LDA,
     $               WORK( LDWRKX*NB+NB+1 ), LDWRKY, ONE,
     $               A( I+NB, I+NB ), LDA )
         CALL DGEMM( 'No transpose', 'No transpose', M-I-NB+1, N-I-NB+1,
     $               NB, -ONE, WORK( NB+1 ), LDWRKX, A( I, I+NB ), LDA,
     $               ONE, A( I+NB, I+NB ), LDA )
*
*        Copy diagonal and off-diagonal elements of B back into A
*
         IF( M.GE.N ) THEN
            DO 10 J = I, I + NB - 1
               A( J, J ) = D( J )
               A( J, J+1 ) = E( J )
   10       CONTINUE
         ELSE
            DO 20 J = I, I + NB - 1
               A( J, J ) = D( J )
               A( J+1, J ) = E( J )
   20       CONTINUE
         END IF
   30 CONTINUE
*
*     Use unblocked code to reduce the remainder of the matrix
*
      CALL DGEBD2( M-I+1, N-I+1, A( I, I ), LDA, D( I ), E( I ),
     $             TAUQ( I ), TAUP( I ), WORK, IINFO )
      WORK( 1 ) = WS
      RETURN
*
*     End of DGEBRD
*
      END
      SUBROUTINE DGELQ2( M, N, A, LDA, TAU, WORK, INFO )
*
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
*     ..
*     .. Array Arguments ..
CID      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
      DIMENSION A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DGELQ2 computes an LQ factorization of a real m by n matrix A:
*  A = L * Q.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the m by n matrix A.
*          On exit, the elements on and below the diagonal of the array
*          contain the m by min(m,n) lower trapezoidal matrix L (L is
*          lower triangular if m <= n); the elements above the diagonal,
*          with the array TAU, represent the orthogonal matrix Q as a
*          product of elementary reflectors (see Further Details).
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  TAU     (output) DOUBLE PRECISION array, dimension (min(M,N))
*          The scalar factors of the elementary reflectors (see Further
*          Details).
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (M)
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*
*  Further Details
*  ===============
*
*  The matrix Q is represented as a product of elementary reflectors
*
*     Q = H(k) . . . H(2) H(1), where k = min(m,n).
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v'
*
*  where tau is a real scalar, and v is a real vector with
*  v(1:i-1) = 0 and v(i) = 1; v(i+1:n) is stored on exit in A(i,i+1:n),
*  and tau in TAU(i).
*
*  =====================================================================
*
*     .. Parameters ..
CID      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, K
CID      DOUBLE PRECISION   AII
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLARF, DLARFG, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGELQ2', -INFO )
         RETURN
      END IF
*
      K = MIN( M, N )
*
      DO 10 I = 1, K
*
*        Generate elementary reflector H(i) to annihilate A(i,i+1:n)
*
         CALL DLARFG( N-I+1, A( I, I ), A( I, MIN( I+1, N ) ), LDA,
     $                TAU( I ) )
         IF( I.LT.M ) THEN
*
*           Apply H(i) to A(i+1:m,i:n) from the right
*
            AII = A( I, I )
            A( I, I ) = ONE
            CALL DLARF( 'Right', M-I, N-I+1, A( I, I ), LDA, TAU( I ),
     $                  A( I+1, I ), LDA, WORK )
            A( I, I ) = AII
         END IF
   10 CONTINUE
      RETURN
*
*     End of DGELQ2
*
      END
      SUBROUTINE DGELQF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
*
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LWORK, M, N
*     ..
*     .. Array Arguments ..
CID      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
      DIMENSION A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DGELQF computes an LQ factorization of a real M-by-N matrix A:
*  A = L * Q.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the M-by-N matrix A.
*          On exit, the elements on and below the diagonal of the array
*          contain the m-by-min(m,n) lower trapezoidal matrix L (L is
*          lower triangular if m <= n); the elements above the diagonal,
*          with the array TAU, represent the orthogonal matrix Q as a
*          product of elementary reflectors (see Further Details).
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  TAU     (output) DOUBLE PRECISION array, dimension (min(M,N))
*          The scalar factors of the elementary reflectors (see Further
*          Details).
*
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.  LWORK >= max(1,M).
*          For optimum performance LWORK >= M*NB, where NB is the
*          optimal blocksize.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  Further Details
*  ===============
*
*  The matrix Q is represented as a product of elementary reflectors
*
*     Q = H(k) . . . H(2) H(1), where k = min(m,n).
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v'
*
*  where tau is a real scalar, and v is a real vector with
*  v(1:i-1) = 0 and v(i) = 1; v(i+1:n) is stored on exit in A(i,i+1:n),
*  and tau in TAU(i).
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IB, IINFO, IWS, K, LDWORK, LWKOPT, NB,
     $                   NBMIN, NX
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGELQ2, DLARFB, DLARFT, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      NB = ILAENV( 1, 'DGELQF', ' ', M, N, -1, -1 )
      LWKOPT = M*NB
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      ELSE IF( LWORK.LT.MAX( 1, M ) .AND. .NOT.LQUERY ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGELQF', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      K = MIN( M, N )
      IF( K.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
      NBMIN = 2
      NX = 0
      IWS = M
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
*
*        Determine when to cross over from blocked to unblocked code.
*
         NX = MAX( 0, ILAENV( 3, 'DGELQF', ' ', M, N, -1, -1 ) )
         IF( NX.LT.K ) THEN
*
*           Determine if workspace is large enough for blocked code.
*
            LDWORK = M
            IWS = LDWORK*NB
            IF( LWORK.LT.IWS ) THEN
*
*              Not enough workspace to use optimal NB:  reduce NB and
*              determine the minimum value of NB.
*
               NB = LWORK / LDWORK
               NBMIN = MAX( 2, ILAENV( 2, 'DGELQF', ' ', M, N, -1,
     $                 -1 ) )
            END IF
         END IF
      END IF
*
      IF( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) THEN
*
*        Use blocked code initially
*
         DO 10 I = 1, K - NX, NB
            IB = MIN( K-I+1, NB )
*
*           Compute the LQ factorization of the current block
*           A(i:i+ib-1,i:n)
*
            CALL DGELQ2( IB, N-I+1, A( I, I ), LDA, TAU( I ), WORK,
     $                   IINFO )
            IF( I+IB.LE.M ) THEN
*
*              Form the triangular factor of the block reflector
*              H = H(i) H(i+1) . . . H(i+ib-1)
*
               CALL DLARFT( 'Forward', 'Rowwise', N-I+1, IB, A( I, I ),
     $                      LDA, TAU( I ), WORK, LDWORK )
*
*              Apply H to A(i+ib:m,i:n) from the right
*
               CALL DLARFB( 'Right', 'No transpose', 'Forward',
     $                      'Rowwise', M-I-IB+1, N-I+1, IB, A( I, I ),
     $                      LDA, WORK, LDWORK, A( I+IB, I ), LDA,
     $                      WORK( IB+1 ), LDWORK )
            END IF
   10    CONTINUE
      ELSE
         I = 1
      END IF
*
*     Use unblocked code to factor the last or only block.
*
      IF( I.LE.K )
     $   CALL DGELQ2( M-I+1, N-I+1, A( I, I ), LDA, TAU( I ), WORK,
     $                IINFO )
*
      WORK( 1 ) = IWS
      RETURN
*
*     End of DGELQF
*
      END
      SUBROUTINE DGEQR2( M, N, A, LDA, TAU, WORK, INFO )
*
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
*     ..
*     .. Array Arguments ..
CID      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
      DIMENSION   A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DGEQR2 computes a QR factorization of a real m by n matrix A:
*  A = Q * R.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the m by n matrix A.
*          On exit, the elements on and above the diagonal of the array
*          contain the min(m,n) by n upper trapezoidal matrix R (R is
*          upper triangular if m >= n); the elements below the diagonal,
*          with the array TAU, represent the orthogonal matrix Q as a
*          product of elementary reflectors (see Further Details).
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  TAU     (output) DOUBLE PRECISION array, dimension (min(M,N))
*          The scalar factors of the elementary reflectors (see Further
*          Details).
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (N)
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*
*  Further Details
*  ===============
*
*  The matrix Q is represented as a product of elementary reflectors
*
*     Q = H(1) H(2) . . . H(k), where k = min(m,n).
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v'
*
*  where tau is a real scalar, and v is a real vector with
*  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
*  and tau in TAU(i).
*
*  =====================================================================
*
*     .. Parameters ..
CID      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, K
CID      DOUBLE PRECISION   AII
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLARF, DLARFG, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGEQR2', -INFO )
         RETURN
      END IF
*
      K = MIN( M, N )
*
      DO 10 I = 1, K
*
*        Generate elementary reflector H(i) to annihilate A(i+1:m,i)
*
         CALL DLARFG( M-I+1, A( I, I ), A( MIN( I+1, M ), I ), 1,
     $                TAU( I ) )
         IF( I.LT.N ) THEN
*
*           Apply H(i) to A(i:m,i+1:n) from the left
*
            AII = A( I, I )
            A( I, I ) = ONE
            CALL DLARF( 'Left', M-I+1, N-I, A( I, I ), 1, TAU( I ),
     $                  A( I, I+1 ), LDA, WORK )
            A( I, I ) = AII
         END IF
   10 CONTINUE
      RETURN
*
*     End of DGEQR2
*
      END
      SUBROUTINE DGEQRF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
*
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LWORK, M, N
*     ..
*     .. Array Arguments ..
CID      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
      DIMENSION    A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DGEQRF computes a QR factorization of a real M-by-N matrix A:
*  A = Q * R.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the M-by-N matrix A.
*          On exit, the elements on and above the diagonal of the array
*          contain the min(M,N)-by-N upper trapezoidal matrix R (R is
*          upper triangular if m >= n); the elements below the diagonal,
*          with the array TAU, represent the orthogonal matrix Q as a
*          product of min(m,n) elementary reflectors (see Further
*          Details).
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  TAU     (output) DOUBLE PRECISION array, dimension (min(M,N))
*          The scalar factors of the elementary reflectors (see Further
*          Details).
*
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.  LWORK >= max(1,N).
*          For optimum performance LWORK >= N*NB, where NB is
*          the optimal blocksize.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  Further Details
*  ===============
*
*  The matrix Q is represented as a product of elementary reflectors
*
*     Q = H(1) H(2) . . . H(k), where k = min(m,n).
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v'
*
*  where tau is a real scalar, and v is a real vector with
*  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
*  and tau in TAU(i).
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IB, IINFO, IWS, K, LDWORK, LWKOPT, NB,
     $                   NBMIN, NX
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEQR2, DLARFB, DLARFT, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      NB = ILAENV( 1, 'DGEQRF', ' ', M, N, -1, -1 )
      LWKOPT = N*NB
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      ELSE IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGEQRF', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      K = MIN( M, N )
      IF( K.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
      NBMIN = 2
      NX = 0
      IWS = N
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
*
*        Determine when to cross over from blocked to unblocked code.
*
         NX = MAX( 0, ILAENV( 3, 'DGEQRF', ' ', M, N, -1, -1 ) )
         IF( NX.LT.K ) THEN
*
*           Determine if workspace is large enough for blocked code.
*
            LDWORK = N
            IWS = LDWORK*NB
            IF( LWORK.LT.IWS ) THEN
*
*              Not enough workspace to use optimal NB:  reduce NB and
*              determine the minimum value of NB.
*
               NB = LWORK / LDWORK
               NBMIN = MAX( 2, ILAENV( 2, 'DGEQRF', ' ', M, N, -1,
     $                 -1 ) )
            END IF
         END IF
      END IF
*
      IF( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) THEN
*
*        Use blocked code initially
*
         DO 10 I = 1, K - NX, NB
            IB = MIN( K-I+1, NB )
*
*           Compute the QR factorization of the current block
*           A(i:m,i:i+ib-1)
*
            CALL DGEQR2( M-I+1, IB, A( I, I ), LDA, TAU( I ), WORK,
     $                   IINFO )
            IF( I+IB.LE.N ) THEN
*
*              Form the triangular factor of the block reflector
*              H = H(i) H(i+1) . . . H(i+ib-1)
*
               CALL DLARFT( 'Forward', 'Columnwise', M-I+1, IB,
     $                      A( I, I ), LDA, TAU( I ), WORK, LDWORK )
*
*              Apply H' to A(i:m,i+ib:n) from the left
*
               CALL DLARFB( 'Left', 'Transpose', 'Forward',
     $                      'Columnwise', M-I+1, N-I-IB+1, IB,
     $                      A( I, I ), LDA, WORK, LDWORK, A( I, I+IB ),
     $                      LDA, WORK( IB+1 ), LDWORK )
            END IF
   10    CONTINUE
      ELSE
         I = 1
      END IF
*
*     Use unblocked code to factor the last or only block.
*
      IF( I.LE.K )
     $   CALL DGEQR2( M-I+1, N-I+1, A( I, I ), LDA, TAU( I ), WORK,
     $                IINFO )
*
      WORK( 1 ) = IWS
      RETURN
*
*     End of DGEQRF
*
      END
      SUBROUTINE DLABRD( M, N, NB, A, LDA, D, E, TAUQ, TAUP, X, LDX, Y,
     $                   LDY )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            LDA, LDX, LDY, M, N, NB
*     ..
*     .. Array Arguments ..
CID      DOUBLE PRECISION   A( LDA, * ), D( * ), E( * ), TAUP( * ),
CID     $                   TAUQ( * ), X( LDX, * ), Y( LDY, * )
      DIMENSION          A( LDA, * ), D( * ), E( * ), TAUP( * ),
     $                   TAUQ( * ), X( LDX, * ), Y( LDY, * )
*     ..
*
*  Purpose
*  =======
*
*  DLABRD reduces the first NB rows and columns of a real general
*  m by n matrix A to upper or lower bidiagonal form by an orthogonal
*  transformation Q' * A * P, and returns the matrices X and Y which
*  are needed to apply the transformation to the unreduced part of A.
*
*  If m >= n, A is reduced to upper bidiagonal form; if m < n, to lower
*  bidiagonal form.
*
*  This is an auxiliary routine called by DGEBRD
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows in the matrix A.
*
*  N       (input) INTEGER
*          The number of columns in the matrix A.
*
*  NB      (input) INTEGER
*          The number of leading rows and columns of A to be reduced.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the m by n general matrix to be reduced.
*          On exit, the first NB rows and columns of the matrix are
*          overwritten; the rest of the array is unchanged.
*          If m >= n, elements on and below the diagonal in the first NB
*            columns, with the array TAUQ, represent the orthogonal
*            matrix Q as a product of elementary reflectors; and
*            elements above the diagonal in the first NB rows, with the
*            array TAUP, represent the orthogonal matrix P as a product
*            of elementary reflectors.
*          If m < n, elements below the diagonal in the first NB
*            columns, with the array TAUQ, represent the orthogonal
*            matrix Q as a product of elementary reflectors, and
*            elements on and above the diagonal in the first NB rows,
*            with the array TAUP, represent the orthogonal matrix P as
*            a product of elementary reflectors.
*          See Further Details.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  D       (output) DOUBLE PRECISION array, dimension (NB)
*          The diagonal elements of the first NB rows and columns of
*          the reduced matrix.  D(i) = A(i,i).
*
*  E       (output) DOUBLE PRECISION array, dimension (NB)
*          The off-diagonal elements of the first NB rows and columns of
*          the reduced matrix.
*
*  TAUQ    (output) DOUBLE PRECISION array dimension (NB)
*          The scalar factors of the elementary reflectors which
*          represent the orthogonal matrix Q. See Further Details.
*
*  TAUP    (output) DOUBLE PRECISION array, dimension (NB)
*          The scalar factors of the elementary reflectors which
*          represent the orthogonal matrix P. See Further Details.
*
*  X       (output) DOUBLE PRECISION array, dimension (LDX,NB)
*          The m-by-nb matrix X required to update the unreduced part
*          of A.
*
*  LDX     (input) INTEGER
*          The leading dimension of the array X. LDX >= M.
*
*  Y       (output) DOUBLE PRECISION array, dimension (LDY,NB)
*          The n-by-nb matrix Y required to update the unreduced part
*          of A.
*
*  LDY     (input) INTEGER
*          The leading dimension of the array Y. LDY >= N.
*
*  Further Details
*  ===============
*
*  The matrices Q and P are represented as products of elementary
*  reflectors:
*
*     Q = H(1) H(2) . . . H(nb)  and  P = G(1) G(2) . . . G(nb)
*
*  Each H(i) and G(i) has the form:
*
*     H(i) = I - tauq * v * v'  and G(i) = I - taup * u * u'
*
*  where tauq and taup are real scalars, and v and u are real vectors.
*
*  If m >= n, v(1:i-1) = 0, v(i) = 1, and v(i:m) is stored on exit in
*  A(i:m,i); u(1:i) = 0, u(i+1) = 1, and u(i+1:n) is stored on exit in
*  A(i,i+1:n); tauq is stored in TAUQ(i) and taup in TAUP(i).
*
*  If m < n, v(1:i) = 0, v(i+1) = 1, and v(i+1:m) is stored on exit in
*  A(i+2:m,i); u(1:i-1) = 0, u(i) = 1, and u(i:n) is stored on exit in
*  A(i,i+1:n); tauq is stored in TAUQ(i) and taup in TAUP(i).
*
*  The elements of the vectors v and u together form the m-by-nb matrix
*  V and the nb-by-n matrix U' which are needed, with X and Y, to apply
*  the transformation to the unreduced part of the matrix, using a block
*  update of the form:  A := A - V*Y' - X*U'.
*
*  The contents of A on exit are illustrated by the following examples
*  with nb = 2:
*
*  m = 6 and n = 5 (m > n):          m = 5 and n = 6 (m < n):
*
*    (  1   1   u1  u1  u1 )           (  1   u1  u1  u1  u1  u1 )
*    (  v1  1   1   u2  u2 )           (  1   1   u2  u2  u2  u2 )
*    (  v1  v2  a   a   a  )           (  v1  1   a   a   a   a  )
*    (  v1  v2  a   a   a  )           (  v1  v2  a   a   a   a  )
*    (  v1  v2  a   a   a  )           (  v1  v2  a   a   a   a  )
*    (  v1  v2  a   a   a  )
*
*  where a denotes an element of the original matrix which is unchanged,
*  vi denotes an element of the vector defining H(i), and ui an element
*  of the vector defining G(i).
*
*  =====================================================================
*
*     .. Parameters ..
CID      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMV, DLARFG, DSCAL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( M.LE.0 .OR. N.LE.0 )
     $   RETURN
*
      IF( M.GE.N ) THEN
*
*        Reduce to upper bidiagonal form
*
         DO 10 I = 1, NB
*
*           Update A(i:m,i)
*
            CALL DGEMV( 'No transpose', M-I+1, I-1, -ONE, A( I, 1 ),
     $                  LDA, Y( I, 1 ), LDY, ONE, A( I, I ), 1 )
            CALL DGEMV( 'No transpose', M-I+1, I-1, -ONE, X( I, 1 ),
     $                  LDX, A( 1, I ), 1, ONE, A( I, I ), 1 )
*
*           Generate reflection Q(i) to annihilate A(i+1:m,i)
*
            CALL DLARFG( M-I+1, A( I, I ), A( MIN( I+1, M ), I ), 1,
     $                   TAUQ( I ) )
            D( I ) = A( I, I )
            IF( I.LT.N ) THEN
               A( I, I ) = ONE
*
*              Compute Y(i+1:n,i)
*
               CALL DGEMV( 'Transpose', M-I+1, N-I, ONE, A( I, I+1 ),
     $                     LDA, A( I, I ), 1, ZERO, Y( I+1, I ), 1 )
               CALL DGEMV( 'Transpose', M-I+1, I-1, ONE, A( I, 1 ), LDA,
     $                     A( I, I ), 1, ZERO, Y( 1, I ), 1 )
               CALL DGEMV( 'No transpose', N-I, I-1, -ONE, Y( I+1, 1 ),
     $                     LDY, Y( 1, I ), 1, ONE, Y( I+1, I ), 1 )
               CALL DGEMV( 'Transpose', M-I+1, I-1, ONE, X( I, 1 ), LDX,
     $                     A( I, I ), 1, ZERO, Y( 1, I ), 1 )
               CALL DGEMV( 'Transpose', I-1, N-I, -ONE, A( 1, I+1 ),
     $                     LDA, Y( 1, I ), 1, ONE, Y( I+1, I ), 1 )
               CALL DSCAL( N-I, TAUQ( I ), Y( I+1, I ), 1 )
*
*              Update A(i,i+1:n)
*
               CALL DGEMV( 'No transpose', N-I, I, -ONE, Y( I+1, 1 ),
     $                     LDY, A( I, 1 ), LDA, ONE, A( I, I+1 ), LDA )
               CALL DGEMV( 'Transpose', I-1, N-I, -ONE, A( 1, I+1 ),
     $                     LDA, X( I, 1 ), LDX, ONE, A( I, I+1 ), LDA )
*
*              Generate reflection P(i) to annihilate A(i,i+2:n)
*
               CALL DLARFG( N-I, A( I, I+1 ), A( I, MIN( I+2, N ) ),
     $                      LDA, TAUP( I ) )
               E( I ) = A( I, I+1 )
               A( I, I+1 ) = ONE
*
*              Compute X(i+1:m,i)
*
               CALL DGEMV( 'No transpose', M-I, N-I, ONE, A( I+1, I+1 ),
     $                     LDA, A( I, I+1 ), LDA, ZERO, X( I+1, I ), 1 )
               CALL DGEMV( 'Transpose', N-I, I, ONE, Y( I+1, 1 ), LDY,
     $                     A( I, I+1 ), LDA, ZERO, X( 1, I ), 1 )
               CALL DGEMV( 'No transpose', M-I, I, -ONE, A( I+1, 1 ),
     $                     LDA, X( 1, I ), 1, ONE, X( I+1, I ), 1 )
               CALL DGEMV( 'No transpose', I-1, N-I, ONE, A( 1, I+1 ),
     $                     LDA, A( I, I+1 ), LDA, ZERO, X( 1, I ), 1 )
               CALL DGEMV( 'No transpose', M-I, I-1, -ONE, X( I+1, 1 ),
     $                     LDX, X( 1, I ), 1, ONE, X( I+1, I ), 1 )
               CALL DSCAL( M-I, TAUP( I ), X( I+1, I ), 1 )
            END IF
   10    CONTINUE
      ELSE
*
*        Reduce to lower bidiagonal form
*
         DO 20 I = 1, NB
*
*           Update A(i,i:n)
*
            CALL DGEMV( 'No transpose', N-I+1, I-1, -ONE, Y( I, 1 ),
     $                  LDY, A( I, 1 ), LDA, ONE, A( I, I ), LDA )
            CALL DGEMV( 'Transpose', I-1, N-I+1, -ONE, A( 1, I ), LDA,
     $                  X( I, 1 ), LDX, ONE, A( I, I ), LDA )
*
*           Generate reflection P(i) to annihilate A(i,i+1:n)
*
            CALL DLARFG( N-I+1, A( I, I ), A( I, MIN( I+1, N ) ), LDA,
     $                   TAUP( I ) )
            D( I ) = A( I, I )
            IF( I.LT.M ) THEN
               A( I, I ) = ONE
*
*              Compute X(i+1:m,i)
*
               CALL DGEMV( 'No transpose', M-I, N-I+1, ONE, A( I+1, I ),
     $                     LDA, A( I, I ), LDA, ZERO, X( I+1, I ), 1 )
               CALL DGEMV( 'Transpose', N-I+1, I-1, ONE, Y( I, 1 ), LDY,
     $                     A( I, I ), LDA, ZERO, X( 1, I ), 1 )
               CALL DGEMV( 'No transpose', M-I, I-1, -ONE, A( I+1, 1 ),
     $                     LDA, X( 1, I ), 1, ONE, X( I+1, I ), 1 )
               CALL DGEMV( 'No transpose', I-1, N-I+1, ONE, A( 1, I ),
     $                     LDA, A( I, I ), LDA, ZERO, X( 1, I ), 1 )
               CALL DGEMV( 'No transpose', M-I, I-1, -ONE, X( I+1, 1 ),
     $                     LDX, X( 1, I ), 1, ONE, X( I+1, I ), 1 )
               CALL DSCAL( M-I, TAUP( I ), X( I+1, I ), 1 )
*
*              Update A(i+1:m,i)
*
               CALL DGEMV( 'No transpose', M-I, I-1, -ONE, A( I+1, 1 ),
     $                     LDA, Y( I, 1 ), LDY, ONE, A( I+1, I ), 1 )
               CALL DGEMV( 'No transpose', M-I, I, -ONE, X( I+1, 1 ),
     $                     LDX, A( 1, I ), 1, ONE, A( I+1, I ), 1 )
*
*              Generate reflection Q(i) to annihilate A(i+2:m,i)
*
               CALL DLARFG( M-I, A( I+1, I ), A( MIN( I+2, M ), I ), 1,
     $                      TAUQ( I ) )
               E( I ) = A( I+1, I )
               A( I+1, I ) = ONE
*
*              Compute Y(i+1:n,i)
*
               CALL DGEMV( 'Transpose', M-I, N-I, ONE, A( I+1, I+1 ),
     $                     LDA, A( I+1, I ), 1, ZERO, Y( I+1, I ), 1 )
               CALL DGEMV( 'Transpose', M-I, I-1, ONE, A( I+1, 1 ), LDA,
     $                     A( I+1, I ), 1, ZERO, Y( 1, I ), 1 )
               CALL DGEMV( 'No transpose', N-I, I-1, -ONE, Y( I+1, 1 ),
     $                     LDY, Y( 1, I ), 1, ONE, Y( I+1, I ), 1 )
               CALL DGEMV( 'Transpose', M-I, I, ONE, X( I+1, 1 ), LDX,
     $                     A( I+1, I ), 1, ZERO, Y( 1, I ), 1 )
               CALL DGEMV( 'Transpose', I, N-I, -ONE, A( 1, I+1 ), LDA,
     $                     Y( 1, I ), 1, ONE, Y( I+1, I ), 1 )
               CALL DSCAL( N-I, TAUQ( I ), Y( I+1, I ), 1 )
            END IF
   20    CONTINUE
      END IF
      RETURN
*
*     End of DLABRD
*
      END
      SUBROUTINE DLACPY( UPLO, M, N, A, LDA, B, LDB )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, LDB, M, N
*     ..
*     .. Array Arguments ..
CID      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
      DIMENSION A( LDA, * ), B( LDB, * )
*     ..
*
*  Purpose
*  =======
*
*  DLACPY copies all or part of a two-dimensional matrix A to another
*  matrix B.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          Specifies the part of the matrix A to be copied to B.
*          = 'U':      Upper triangular part
*          = 'L':      Lower triangular part
*          Otherwise:  All of the matrix A
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
*          The m by n matrix A.  If UPLO = 'U', only the upper triangle
*          or trapezoid is accessed; if UPLO = 'L', only the lower
*          triangle or trapezoid is accessed.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  B       (output) DOUBLE PRECISION array, dimension (LDB,N)
*          On exit, B = A in the locations specified by UPLO.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,M).
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, J
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN
*     ..
*     .. Executable Statements ..
*
      IF( LSAME( UPLO, 'U' ) ) THEN
         DO 20 J = 1, N
            DO 10 I = 1, MIN( J, M )
               B( I, J ) = A( I, J )
   10       CONTINUE
   20    CONTINUE
      ELSE IF( LSAME( UPLO, 'L' ) ) THEN
         DO 40 J = 1, N
            DO 30 I = J, M
               B( I, J ) = A( I, J )
   30       CONTINUE
   40    CONTINUE
      ELSE
         DO 60 J = 1, N
            DO 50 I = 1, M
               B( I, J ) = A( I, J )
   50       CONTINUE
   60    CONTINUE
      END IF
      RETURN
*
*     End of DLACPY
*
      END
C      
CID      DOUBLE PRECISION FUNCTION DLANGE( NORM, M, N, A, LDA, WORK )
      FUNCTION DLANGE( NORM, M, N, A, LDA, WORK )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          NORM
      INTEGER            LDA, M, N
*     ..
*     .. Array Arguments ..
CID      DOUBLE PRECISION   A( LDA, * ), WORK( * )
      DIMENSION    A( LDA, * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DLANGE  returns the value of the one norm,  or the Frobenius norm, or
*  the  infinity norm,  or the  element of  largest absolute value  of a
*  real matrix A.
*
*  Description
*  ===========
*
*  DLANGE returns the value
*
*     DLANGE = ( max(abs(A(i,j))), NORM = 'M' or 'm'
*              (
*              ( norm1(A),         NORM = '1', 'O' or 'o'
*              (
*              ( normI(A),         NORM = 'I' or 'i'
*              (
*              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
*
*  where  norm1  denotes the  one norm of a matrix (maximum column sum),
*  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
*  normF  denotes the  Frobenius norm of a matrix (square root of sum of
*  squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.
*
*  Arguments
*  =========
*
*  NORM    (input) CHARACTER*1
*          Specifies the value to be returned in DLANGE as described
*          above.
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.  When M = 0,
*          DLANGE is set to zero.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.  When N = 0,
*          DLANGE is set to zero.
*
*  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
*          The m by n matrix A.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(M,1).
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (MAX(1,LWORK)),
*          where LWORK >= M when NORM = 'I'; otherwise, WORK is not
*          referenced.
*
* =====================================================================
*
*     .. Parameters ..
CID      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J
CID      DOUBLE PRECISION   SCALE, SUM, VALUE
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLASSQ
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
*
      IF( MIN( M, N ).EQ.0 ) THEN
         VALUE = ZERO
      ELSE IF( LSAME( NORM, 'M' ) ) THEN
*
*        Find max(abs(A(i,j))).
*
         VALUE = ZERO
         DO 20 J = 1, N
            DO 10 I = 1, M
               VALUE = MAX( VALUE, ABS( A( I, J ) ) )
   10       CONTINUE
   20    CONTINUE
      ELSE IF( ( LSAME( NORM, 'O' ) ) .OR. ( NORM.EQ.'1' ) ) THEN
*
*        Find norm1(A).
*
         VALUE = ZERO
         DO 40 J = 1, N
            SUM = ZERO
            DO 30 I = 1, M
               SUM = SUM + ABS( A( I, J ) )
   30       CONTINUE
            VALUE = MAX( VALUE, SUM )
   40    CONTINUE
      ELSE IF( LSAME( NORM, 'I' ) ) THEN
*
*        Find normI(A).
*
         DO 50 I = 1, M
            WORK( I ) = ZERO
   50    CONTINUE
         DO 70 J = 1, N
            DO 60 I = 1, M
               WORK( I ) = WORK( I ) + ABS( A( I, J ) )
   60       CONTINUE
   70    CONTINUE
         VALUE = ZERO
         DO 80 I = 1, M
            VALUE = MAX( VALUE, WORK( I ) )
   80    CONTINUE
      ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN
*
*        Find normF(A).
*
         SCALE = ZERO
         SUM = ONE
         DO 90 J = 1, N
            CALL DLASSQ( M, A( 1, J ), 1, SCALE, SUM )
   90    CONTINUE
         VALUE = SCALE*SQRT( SUM )
      END IF
*
      DLANGE = VALUE
      RETURN
*
*     End of DLANGE
*
      END
CID      DOUBLE PRECISION FUNCTION DLAPY2( X, Y )
      FUNCTION DLAPY2( X, Y )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
CID      DOUBLE PRECISION   X, Y
*     ..
*
*  Purpose
*  =======
*
*  DLAPY2 returns sqrt(x**2+y**2), taking care not to cause unnecessary
*  overflow.
*
*  Arguments
*  =========
*
*  X       (input) DOUBLE PRECISION
*  Y       (input) DOUBLE PRECISION
*          X and Y specify the values x and y.
*
*  =====================================================================
*
*     .. Parameters ..
CID      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
CID      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
CID      DOUBLE PRECISION   W, XABS, YABS, Z
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
*
      XABS = ABS( X )
      YABS = ABS( Y )
      W = MAX( XABS, YABS )
      Z = MIN( XABS, YABS )
      IF( Z.EQ.ZERO ) THEN
         DLAPY2 = W
      ELSE
         DLAPY2 = W*SQRT( ONE+( Z / W )**2 )
      END IF
      RETURN
*
*     End of DLAPY2
*
      END
      SUBROUTINE DLARF( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          SIDE
      INTEGER            INCV, LDC, M, N
CID      DOUBLE PRECISION   TAU
*     ..
*     .. Array Arguments ..
CID      DOUBLE PRECISION   C( LDC, * ), V( * ), WORK( * )
      DIMENSION    C( LDC, * ), V( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DLARF applies a real elementary reflector H to a real m by n matrix
*  C, from either the left or the right. H is represented in the form
*
*        H = I - tau * v * v'
*
*  where tau is a real scalar and v is a real vector.
*
*  If tau = 0, then H is taken to be the unit matrix.
*
*  Arguments
*  =========
*
*  SIDE    (input) CHARACTER*1
*          = 'L': form  H * C
*          = 'R': form  C * H
*
*  M       (input) INTEGER
*          The number of rows of the matrix C.
*
*  N       (input) INTEGER
*          The number of columns of the matrix C.
*
*  V       (input) DOUBLE PRECISION array, dimension
*                     (1 + (M-1)*abs(INCV)) if SIDE = 'L'
*                  or (1 + (N-1)*abs(INCV)) if SIDE = 'R'
*          The vector v in the representation of H. V is not used if
*          TAU = 0.
*
*  INCV    (input) INTEGER
*          The increment between elements of v. INCV <> 0.
*
*  TAU     (input) DOUBLE PRECISION
*          The value tau in the representation of H.
*
*  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
*          On entry, the m by n matrix C.
*          On exit, C is overwritten by the matrix H * C if SIDE = 'L',
*          or C * H if SIDE = 'R'.
*
*  LDC     (input) INTEGER
*          The leading dimension of the array C. LDC >= max(1,M).
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension
*                         (N) if SIDE = 'L'
*                      or (M) if SIDE = 'R'
*
*  =====================================================================
*
*     .. Parameters ..
CID      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMV, DGER
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Executable Statements ..
*
      IF( LSAME( SIDE, 'L' ) ) THEN
*
*        Form  H * C
*
         IF( TAU.NE.ZERO ) THEN
*
*           w := C' * v
*
            CALL DGEMV( 'Transpose', M, N, ONE, C, LDC, V, INCV, ZERO,
     $                  WORK, 1 )
*
*           C := C - v * w'
*
            CALL DGER( M, N, -TAU, V, INCV, WORK, 1, C, LDC )
         END IF
      ELSE
*
*        Form  C * H
*
         IF( TAU.NE.ZERO ) THEN
*
*           w := C * v
*
            CALL DGEMV( 'No transpose', M, N, ONE, C, LDC, V, INCV,
     $                  ZERO, WORK, 1 )
*
*           C := C - w * v'
*
            CALL DGER( M, N, -TAU, WORK, 1, V, INCV, C, LDC )
         END IF
      END IF
      RETURN
*
*     End of DLARF
*
      END
      SUBROUTINE DLARFB( SIDE, TRANS, DIRECT, STOREV, M, N, K, V, LDV,
     $                   T, LDT, C, LDC, WORK, LDWORK )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          DIRECT, SIDE, STOREV, TRANS
      INTEGER            K, LDC, LDT, LDV, LDWORK, M, N
*     ..
*     .. Array Arguments ..
CID      DOUBLE PRECISION   C( LDC, * ), T( LDT, * ), V( LDV, * ),
CID     $                   WORK( LDWORK, * )
      DIMENSION    C( LDC, * ), T( LDT, * ), V( LDV, * ),
     $                   WORK( LDWORK, * )
*     ..
*
*  Purpose
*  =======
*
*  DLARFB applies a real block reflector H or its transpose H' to a
*  real m by n matrix C, from either the left or the right.
*
*  Arguments
*  =========
*
*  SIDE    (input) CHARACTER*1
*          = 'L': apply H or H' from the Left
*          = 'R': apply H or H' from the Right
*
*  TRANS   (input) CHARACTER*1
*          = 'N': apply H (No transpose)
*          = 'T': apply H' (Transpose)
*
*  DIRECT  (input) CHARACTER*1
*          Indicates how H is formed from a product of elementary
*          reflectors
*          = 'F': H = H(1) H(2) . . . H(k) (Forward)
*          = 'B': H = H(k) . . . H(2) H(1) (Backward)
*
*  STOREV  (input) CHARACTER*1
*          Indicates how the vectors which define the elementary
*          reflectors are stored:
*          = 'C': Columnwise
*          = 'R': Rowwise
*
*  M       (input) INTEGER
*          The number of rows of the matrix C.
*
*  N       (input) INTEGER
*          The number of columns of the matrix C.
*
*  K       (input) INTEGER
*          The order of the matrix T (= the number of elementary
*          reflectors whose product defines the block reflector).
*
*  V       (input) DOUBLE PRECISION array, dimension
*                                (LDV,K) if STOREV = 'C'
*                                (LDV,M) if STOREV = 'R' and SIDE = 'L'
*                                (LDV,N) if STOREV = 'R' and SIDE = 'R'
*          The matrix V. See further details.
*
*  LDV     (input) INTEGER
*          The leading dimension of the array V.
*          If STOREV = 'C' and SIDE = 'L', LDV >= max(1,M);
*          if STOREV = 'C' and SIDE = 'R', LDV >= max(1,N);
*          if STOREV = 'R', LDV >= K.
*
*  T       (input) DOUBLE PRECISION array, dimension (LDT,K)
*          The triangular k by k matrix T in the representation of the
*          block reflector.
*
*  LDT     (input) INTEGER
*          The leading dimension of the array T. LDT >= K.
*
*  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
*          On entry, the m by n matrix C.
*          On exit, C is overwritten by H*C or H'*C or C*H or C*H'.
*
*  LDC     (input) INTEGER
*          The leading dimension of the array C. LDA >= max(1,M).
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (LDWORK,K)
*
*  LDWORK  (input) INTEGER
*          The leading dimension of the array WORK.
*          If SIDE = 'L', LDWORK >= max(1,N);
*          if SIDE = 'R', LDWORK >= max(1,M).
*
*  =====================================================================
*
*     .. Parameters ..
CID      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      CHARACTER          TRANST
      INTEGER            I, J
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DCOPY, DGEMM, DTRMM
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( M.LE.0 .OR. N.LE.0 )
     $   RETURN
*
      IF( LSAME( TRANS, 'N' ) ) THEN
         TRANST = 'T'
      ELSE
         TRANST = 'N'
      END IF
*
      IF( LSAME( STOREV, 'C' ) ) THEN
*
         IF( LSAME( DIRECT, 'F' ) ) THEN
*
*           Let  V =  ( V1 )    (first K rows)
*                     ( V2 )
*           where  V1  is unit lower triangular.
*
            IF( LSAME( SIDE, 'L' ) ) THEN
*
*              Form  H * C  or  H' * C  where  C = ( C1 )
*                                                  ( C2 )
*
*              W := C' * V  =  (C1'*V1 + C2'*V2)  (stored in WORK)
*
*              W := C1'
*
               DO 10 J = 1, K
                  CALL DCOPY( N, C( J, 1 ), LDC, WORK( 1, J ), 1 )
   10          CONTINUE
*
*              W := W * V1
*
               CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Unit', N,
     $                     K, ONE, V, LDV, WORK, LDWORK )
               IF( M.GT.K ) THEN
*
*                 W := W + C2'*V2
*
                  CALL DGEMM( 'Transpose', 'No transpose', N, K, M-K,
     $                        ONE, C( K+1, 1 ), LDC, V( K+1, 1 ), LDV,
     $                        ONE, WORK, LDWORK )
               END IF
*
*              W := W * T'  or  W * T
*
               CALL DTRMM( 'Right', 'Upper', TRANST, 'Non-unit', N, K,
     $                     ONE, T, LDT, WORK, LDWORK )
*
*              C := C - V * W'
*
               IF( M.GT.K ) THEN
*
*                 C2 := C2 - V2 * W'
*
                  CALL DGEMM( 'No transpose', 'Transpose', M-K, N, K,
     $                        -ONE, V( K+1, 1 ), LDV, WORK, LDWORK, ONE,
     $                        C( K+1, 1 ), LDC )
               END IF
*
*              W := W * V1'
*
               CALL DTRMM( 'Right', 'Lower', 'Transpose', 'Unit', N, K,
     $                     ONE, V, LDV, WORK, LDWORK )
*
*              C1 := C1 - W'
*
               DO 30 J = 1, K
                  DO 20 I = 1, N
                     C( J, I ) = C( J, I ) - WORK( I, J )
   20             CONTINUE
   30          CONTINUE
*
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
*
*              Form  C * H  or  C * H'  where  C = ( C1  C2 )
*
*              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
*
*              W := C1
*
               DO 40 J = 1, K
                  CALL DCOPY( M, C( 1, J ), 1, WORK( 1, J ), 1 )
   40          CONTINUE
*
*              W := W * V1
*
               CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Unit', M,
     $                     K, ONE, V, LDV, WORK, LDWORK )
               IF( N.GT.K ) THEN
*
*                 W := W + C2 * V2
*
                  CALL DGEMM( 'No transpose', 'No transpose', M, K, N-K,
     $                        ONE, C( 1, K+1 ), LDC, V( K+1, 1 ), LDV,
     $                        ONE, WORK, LDWORK )
               END IF
*
*              W := W * T  or  W * T'
*
               CALL DTRMM( 'Right', 'Upper', TRANS, 'Non-unit', M, K,
     $                     ONE, T, LDT, WORK, LDWORK )
*
*              C := C - W * V'
*
               IF( N.GT.K ) THEN
*
*                 C2 := C2 - W * V2'
*
                  CALL DGEMM( 'No transpose', 'Transpose', M, N-K, K,
     $                        -ONE, WORK, LDWORK, V( K+1, 1 ), LDV, ONE,
     $                        C( 1, K+1 ), LDC )
               END IF
*
*              W := W * V1'
*
               CALL DTRMM( 'Right', 'Lower', 'Transpose', 'Unit', M, K,
     $                     ONE, V, LDV, WORK, LDWORK )
*
*              C1 := C1 - W
*
               DO 60 J = 1, K
                  DO 50 I = 1, M
                     C( I, J ) = C( I, J ) - WORK( I, J )
   50             CONTINUE
   60          CONTINUE
            END IF
*
         ELSE
*
*           Let  V =  ( V1 )
*                     ( V2 )    (last K rows)
*           where  V2  is unit upper triangular.
*
            IF( LSAME( SIDE, 'L' ) ) THEN
*
*              Form  H * C  or  H' * C  where  C = ( C1 )
*                                                  ( C2 )
*
*              W := C' * V  =  (C1'*V1 + C2'*V2)  (stored in WORK)
*
*              W := C2'
*
               DO 70 J = 1, K
                  CALL DCOPY( N, C( M-K+J, 1 ), LDC, WORK( 1, J ), 1 )
   70          CONTINUE
*
*              W := W * V2
*
               CALL DTRMM( 'Right', 'Upper', 'No transpose', 'Unit', N,
     $                     K, ONE, V( M-K+1, 1 ), LDV, WORK, LDWORK )
               IF( M.GT.K ) THEN
*
*                 W := W + C1'*V1
*
                  CALL DGEMM( 'Transpose', 'No transpose', N, K, M-K,
     $                        ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
               END IF
*
*              W := W * T'  or  W * T
*
               CALL DTRMM( 'Right', 'Lower', TRANST, 'Non-unit', N, K,
     $                     ONE, T, LDT, WORK, LDWORK )
*
*              C := C - V * W'
*
               IF( M.GT.K ) THEN
*
*                 C1 := C1 - V1 * W'
*
                  CALL DGEMM( 'No transpose', 'Transpose', M-K, N, K,
     $                        -ONE, V, LDV, WORK, LDWORK, ONE, C, LDC )
               END IF
*
*              W := W * V2'
*
               CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Unit', N, K,
     $                     ONE, V( M-K+1, 1 ), LDV, WORK, LDWORK )
*
*              C2 := C2 - W'
*
               DO 90 J = 1, K
                  DO 80 I = 1, N
                     C( M-K+J, I ) = C( M-K+J, I ) - WORK( I, J )
   80             CONTINUE
   90          CONTINUE
*
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
*
*              Form  C * H  or  C * H'  where  C = ( C1  C2 )
*
*              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
*
*              W := C2
*
               DO 100 J = 1, K
                  CALL DCOPY( M, C( 1, N-K+J ), 1, WORK( 1, J ), 1 )
  100          CONTINUE
*
*              W := W * V2
*
               CALL DTRMM( 'Right', 'Upper', 'No transpose', 'Unit', M,
     $                     K, ONE, V( N-K+1, 1 ), LDV, WORK, LDWORK )
               IF( N.GT.K ) THEN
*
*                 W := W + C1 * V1
*
                  CALL DGEMM( 'No transpose', 'No transpose', M, K, N-K,
     $                        ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
               END IF
*
*              W := W * T  or  W * T'
*
               CALL DTRMM( 'Right', 'Lower', TRANS, 'Non-unit', M, K,
     $                     ONE, T, LDT, WORK, LDWORK )
*
*              C := C - W * V'
*
               IF( N.GT.K ) THEN
*
*                 C1 := C1 - W * V1'
*
                  CALL DGEMM( 'No transpose', 'Transpose', M, N-K, K,
     $                        -ONE, WORK, LDWORK, V, LDV, ONE, C, LDC )
               END IF
*
*              W := W * V2'
*
               CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Unit', M, K,
     $                     ONE, V( N-K+1, 1 ), LDV, WORK, LDWORK )
*
*              C2 := C2 - W
*
               DO 120 J = 1, K
                  DO 110 I = 1, M
                     C( I, N-K+J ) = C( I, N-K+J ) - WORK( I, J )
  110             CONTINUE
  120          CONTINUE
            END IF
         END IF
*
      ELSE IF( LSAME( STOREV, 'R' ) ) THEN
*
         IF( LSAME( DIRECT, 'F' ) ) THEN
*
*           Let  V =  ( V1  V2 )    (V1: first K columns)
*           where  V1  is unit upper triangular.
*
            IF( LSAME( SIDE, 'L' ) ) THEN
*
*              Form  H * C  or  H' * C  where  C = ( C1 )
*                                                  ( C2 )
*
*              W := C' * V'  =  (C1'*V1' + C2'*V2') (stored in WORK)
*
*              W := C1'
*
               DO 130 J = 1, K
                  CALL DCOPY( N, C( J, 1 ), LDC, WORK( 1, J ), 1 )
  130          CONTINUE
*
*              W := W * V1'
*
               CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Unit', N, K,
     $                     ONE, V, LDV, WORK, LDWORK )
               IF( M.GT.K ) THEN
*
*                 W := W + C2'*V2'
*
                  CALL DGEMM( 'Transpose', 'Transpose', N, K, M-K, ONE,
     $                        C( K+1, 1 ), LDC, V( 1, K+1 ), LDV, ONE,
     $                        WORK, LDWORK )
               END IF
*
*              W := W * T'  or  W * T
*
               CALL DTRMM( 'Right', 'Upper', TRANST, 'Non-unit', N, K,
     $                     ONE, T, LDT, WORK, LDWORK )
*
*              C := C - V' * W'
*
               IF( M.GT.K ) THEN
*
*                 C2 := C2 - V2' * W'
*
                  CALL DGEMM( 'Transpose', 'Transpose', M-K, N, K, -ONE,
     $                        V( 1, K+1 ), LDV, WORK, LDWORK, ONE,
     $                        C( K+1, 1 ), LDC )
               END IF
*
*              W := W * V1
*
               CALL DTRMM( 'Right', 'Upper', 'No transpose', 'Unit', N,
     $                     K, ONE, V, LDV, WORK, LDWORK )
*
*              C1 := C1 - W'
*
               DO 150 J = 1, K
                  DO 140 I = 1, N
                     C( J, I ) = C( J, I ) - WORK( I, J )
  140             CONTINUE
  150          CONTINUE
*
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
*
*              Form  C * H  or  C * H'  where  C = ( C1  C2 )
*
*              W := C * V'  =  (C1*V1' + C2*V2')  (stored in WORK)
*
*              W := C1
*
               DO 160 J = 1, K
                  CALL DCOPY( M, C( 1, J ), 1, WORK( 1, J ), 1 )
  160          CONTINUE
*
*              W := W * V1'
*
               CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Unit', M, K,
     $                     ONE, V, LDV, WORK, LDWORK )
               IF( N.GT.K ) THEN
*
*                 W := W + C2 * V2'
*
                  CALL DGEMM( 'No transpose', 'Transpose', M, K, N-K,
     $                        ONE, C( 1, K+1 ), LDC, V( 1, K+1 ), LDV,
     $                        ONE, WORK, LDWORK )
               END IF
*
*              W := W * T  or  W * T'
*
               CALL DTRMM( 'Right', 'Upper', TRANS, 'Non-unit', M, K,
     $                     ONE, T, LDT, WORK, LDWORK )
*
*              C := C - W * V
*
               IF( N.GT.K ) THEN
*
*                 C2 := C2 - W * V2
*
                  CALL DGEMM( 'No transpose', 'No transpose', M, N-K, K,
     $                        -ONE, WORK, LDWORK, V( 1, K+1 ), LDV, ONE,
     $                        C( 1, K+1 ), LDC )
               END IF
*
*              W := W * V1
*
               CALL DTRMM( 'Right', 'Upper', 'No transpose', 'Unit', M,
     $                     K, ONE, V, LDV, WORK, LDWORK )
*
*              C1 := C1 - W
*
               DO 180 J = 1, K
                  DO 170 I = 1, M
                     C( I, J ) = C( I, J ) - WORK( I, J )
  170             CONTINUE
  180          CONTINUE
*
            END IF
*
         ELSE
*
*           Let  V =  ( V1  V2 )    (V2: last K columns)
*           where  V2  is unit lower triangular.
*
            IF( LSAME( SIDE, 'L' ) ) THEN
*
*              Form  H * C  or  H' * C  where  C = ( C1 )
*                                                  ( C2 )
*
*              W := C' * V'  =  (C1'*V1' + C2'*V2') (stored in WORK)
*
*              W := C2'
*
               DO 190 J = 1, K
                  CALL DCOPY( N, C( M-K+J, 1 ), LDC, WORK( 1, J ), 1 )
  190          CONTINUE
*
*              W := W * V2'
*
               CALL DTRMM( 'Right', 'Lower', 'Transpose', 'Unit', N, K,
     $                     ONE, V( 1, M-K+1 ), LDV, WORK, LDWORK )
               IF( M.GT.K ) THEN
*
*                 W := W + C1'*V1'
*
                  CALL DGEMM( 'Transpose', 'Transpose', N, K, M-K, ONE,
     $                        C, LDC, V, LDV, ONE, WORK, LDWORK )
               END IF
*
*              W := W * T'  or  W * T
*
               CALL DTRMM( 'Right', 'Lower', TRANST, 'Non-unit', N, K,
     $                     ONE, T, LDT, WORK, LDWORK )
*
*              C := C - V' * W'
*
               IF( M.GT.K ) THEN
*
*                 C1 := C1 - V1' * W'
*
                  CALL DGEMM( 'Transpose', 'Transpose', M-K, N, K, -ONE,
     $                        V, LDV, WORK, LDWORK, ONE, C, LDC )
               END IF
*
*              W := W * V2
*
               CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Unit', N,
     $                     K, ONE, V( 1, M-K+1 ), LDV, WORK, LDWORK )
*
*              C2 := C2 - W'
*
               DO 210 J = 1, K
                  DO 200 I = 1, N
                     C( M-K+J, I ) = C( M-K+J, I ) - WORK( I, J )
  200             CONTINUE
  210          CONTINUE
*
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
*
*              Form  C * H  or  C * H'  where  C = ( C1  C2 )
*
*              W := C * V'  =  (C1*V1' + C2*V2')  (stored in WORK)
*
*              W := C2
*
               DO 220 J = 1, K
                  CALL DCOPY( M, C( 1, N-K+J ), 1, WORK( 1, J ), 1 )
  220          CONTINUE
*
*              W := W * V2'
*
               CALL DTRMM( 'Right', 'Lower', 'Transpose', 'Unit', M, K,
     $                     ONE, V( 1, N-K+1 ), LDV, WORK, LDWORK )
               IF( N.GT.K ) THEN
*
*                 W := W + C1 * V1'
*
                  CALL DGEMM( 'No transpose', 'Transpose', M, K, N-K,
     $                        ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
               END IF
*
*              W := W * T  or  W * T'
*
               CALL DTRMM( 'Right', 'Lower', TRANS, 'Non-unit', M, K,
     $                     ONE, T, LDT, WORK, LDWORK )
*
*              C := C - W * V
*
               IF( N.GT.K ) THEN
*
*                 C1 := C1 - W * V1
*
                  CALL DGEMM( 'No transpose', 'No transpose', M, N-K, K,
     $                        -ONE, WORK, LDWORK, V, LDV, ONE, C, LDC )
               END IF
*
*              W := W * V2
*
               CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Unit', M,
     $                     K, ONE, V( 1, N-K+1 ), LDV, WORK, LDWORK )
*
*              C1 := C1 - W
*
               DO 240 J = 1, K
                  DO 230 I = 1, M
                     C( I, N-K+J ) = C( I, N-K+J ) - WORK( I, J )
  230             CONTINUE
  240          CONTINUE
*
            END IF
*
         END IF
      END IF
*
      RETURN
*
*     End of DLARFB
*
      END
      SUBROUTINE DLARFG( N, ALPHA, X, INCX, TAU )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INCX, N
CID      DOUBLE PRECISION   ALPHA, TAU
*     ..
*     .. Array Arguments ..
CID      DOUBLE PRECISION   X( * )
      DIMENSION X( * )
*     ..
*
*  Purpose
*  =======
*
*  DLARFG generates a real elementary reflector H of order n, such
*  that
*
*        H * ( alpha ) = ( beta ),   H' * H = I.
*            (   x   )   (   0  )
*
*  where alpha and beta are scalars, and x is an (n-1)-element real
*  vector. H is represented in the form
*
*        H = I - tau * ( 1 ) * ( 1 v' ) ,
*                      ( v )
*
*  where tau is a real scalar and v is a real (n-1)-element
*  vector.
*
*  If the elements of x are all zero, then tau = 0 and H is taken to be
*  the unit matrix.
*
*  Otherwise  1 <= tau <= 2.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the elementary reflector.
*
*  ALPHA   (input/output) DOUBLE PRECISION
*          On entry, the value alpha.
*          On exit, it is overwritten with the value beta.
*
*  X       (input/output) DOUBLE PRECISION array, dimension
*                         (1+(N-2)*abs(INCX))
*          On entry, the vector x.
*          On exit, it is overwritten with the vector v.
*
*  INCX    (input) INTEGER
*          The increment between elements of X. INCX > 0.
*
*  TAU     (output) DOUBLE PRECISION
*          The value tau.
*
*  =====================================================================
*
*     .. Parameters ..
CID      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            J, KNT
CID      DOUBLE PRECISION   BETA, RSAFMN, SAFMIN, XNORM
*     ..
*     .. External Functions ..
CID      DOUBLE PRECISION   DLAMCH, DLAPY2, DNRM2
      EXTERNAL           DLAMCH, DLAPY2, DNRM2
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, SIGN
*     ..
*     .. External Subroutines ..
      EXTERNAL           DSCAL
*     ..
*     .. Executable Statements ..
*
      IF( N.LE.1 ) THEN
         TAU = ZERO
         RETURN
      END IF
*
      XNORM = DNRM2( N-1, X, INCX )
*
      IF( XNORM.EQ.ZERO ) THEN
*
*        H  =  I
*
         TAU = ZERO
      ELSE
*
*        general case
*
         BETA = -SIGN( DLAPY2( ALPHA, XNORM ), ALPHA )
         SAFMIN = DLAMCH( 'S' ) / DLAMCH( 'E' )
         IF( ABS( BETA ).LT.SAFMIN ) THEN
*
*           XNORM, BETA may be inaccurate; scale X and recompute them
*
            RSAFMN = ONE / SAFMIN
            KNT = 0
   10       CONTINUE
            KNT = KNT + 1
            CALL DSCAL( N-1, RSAFMN, X, INCX )
            BETA = BETA*RSAFMN
            ALPHA = ALPHA*RSAFMN
            IF( ABS( BETA ).LT.SAFMIN )
     $         GO TO 10
*
*           New BETA is at most 1, at least SAFMIN
*
            XNORM = DNRM2( N-1, X, INCX )
            BETA = -SIGN( DLAPY2( ALPHA, XNORM ), ALPHA )
            TAU = ( BETA-ALPHA ) / BETA
            CALL DSCAL( N-1, ONE / ( ALPHA-BETA ), X, INCX )
*
*           If ALPHA is subnormal, it may lose relative accuracy
*
            ALPHA = BETA
            DO 20 J = 1, KNT
               ALPHA = ALPHA*SAFMIN
   20       CONTINUE
         ELSE
            TAU = ( BETA-ALPHA ) / BETA
            CALL DSCAL( N-1, ONE / ( ALPHA-BETA ), X, INCX )
            ALPHA = BETA
         END IF
      END IF
*
      RETURN
*
*     End of DLARFG
*
      END
      SUBROUTINE DLARFT( DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          DIRECT, STOREV
      INTEGER            K, LDT, LDV, N
*     ..
*     .. Array Arguments ..
CID      DOUBLE PRECISION   T( LDT, * ), TAU( * ), V( LDV, * )
      DIMENSION    T( LDT, * ), TAU( * ), V( LDV, * )
*     ..
*
*  Purpose
*  =======
*
*  DLARFT forms the triangular factor T of a real block reflector H
*  of order n, which is defined as a product of k elementary reflectors.
*
*  If DIRECT = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular;
*
*  If DIRECT = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular.
*
*  If STOREV = 'C', the vector which defines the elementary reflector
*  H(i) is stored in the i-th column of the array V, and
*
*     H  =  I - V * T * V'
*
*  If STOREV = 'R', the vector which defines the elementary reflector
*  H(i) is stored in the i-th row of the array V, and
*
*     H  =  I - V' * T * V
*
*  Arguments
*  =========
*
*  DIRECT  (input) CHARACTER*1
*          Specifies the order in which the elementary reflectors are
*          multiplied to form the block reflector:
*          = 'F': H = H(1) H(2) . . . H(k) (Forward)
*          = 'B': H = H(k) . . . H(2) H(1) (Backward)
*
*  STOREV  (input) CHARACTER*1
*          Specifies how the vectors which define the elementary
*          reflectors are stored (see also Further Details):
*          = 'C': columnwise
*          = 'R': rowwise
*
*  N       (input) INTEGER
*          The order of the block reflector H. N >= 0.
*
*  K       (input) INTEGER
*          The order of the triangular factor T (= the number of
*          elementary reflectors). K >= 1.
*
*  V       (input/output) DOUBLE PRECISION array, dimension
*                               (LDV,K) if STOREV = 'C'
*                               (LDV,N) if STOREV = 'R'
*          The matrix V. See further details.
*
*  LDV     (input) INTEGER
*          The leading dimension of the array V.
*          If STOREV = 'C', LDV >= max(1,N); if STOREV = 'R', LDV >= K.
*
*  TAU     (input) DOUBLE PRECISION array, dimension (K)
*          TAU(i) must contain the scalar factor of the elementary
*          reflector H(i).
*
*  T       (output) DOUBLE PRECISION array, dimension (LDT,K)
*          The k by k triangular factor T of the block reflector.
*          If DIRECT = 'F', T is upper triangular; if DIRECT = 'B', T is
*          lower triangular. The rest of the array is not used.
*
*  LDT     (input) INTEGER
*          The leading dimension of the array T. LDT >= K.
*
*  Further Details
*  ===============
*
*  The shape of the matrix V and the storage of the vectors which define
*  the H(i) is best illustrated by the following example with n = 5 and
*  k = 3. The elements equal to 1 are not stored; the corresponding
*  array elements are modified but restored on exit. The rest of the
*  array is not used.
*
*  DIRECT = 'F' and STOREV = 'C':         DIRECT = 'F' and STOREV = 'R':
*
*               V = (  1       )                 V = (  1 v1 v1 v1 v1 )
*                   ( v1  1    )                     (     1 v2 v2 v2 )
*                   ( v1 v2  1 )                     (        1 v3 v3 )
*                   ( v1 v2 v3 )
*                   ( v1 v2 v3 )
*
*  DIRECT = 'B' and STOREV = 'C':         DIRECT = 'B' and STOREV = 'R':
*
*               V = ( v1 v2 v3 )                 V = ( v1 v1  1       )
*                   ( v1 v2 v3 )                     ( v2 v2 v2  1    )
*                   (  1 v2 v3 )                     ( v3 v3 v3 v3  1 )
*                   (     1 v3 )
*                   (        1 )
*
*  =====================================================================
*
*     .. Parameters ..
CID      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J
CID      DOUBLE PRECISION   VII
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMV, DTRMV
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      IF( LSAME( DIRECT, 'F' ) ) THEN
         DO 20 I = 1, K
            IF( TAU( I ).EQ.ZERO ) THEN
*
*              H(i)  =  I
*
               DO 10 J = 1, I
                  T( J, I ) = ZERO
   10          CONTINUE
            ELSE
*
*              general case
*
               VII = V( I, I )
               V( I, I ) = ONE
               IF( LSAME( STOREV, 'C' ) ) THEN
*
*                 T(1:i-1,i) := - tau(i) * V(i:n,1:i-1)' * V(i:n,i)
*
                  CALL DGEMV( 'Transpose', N-I+1, I-1, -TAU( I ),
     $                        V( I, 1 ), LDV, V( I, I ), 1, ZERO,
     $                        T( 1, I ), 1 )
               ELSE
*
*                 T(1:i-1,i) := - tau(i) * V(1:i-1,i:n) * V(i,i:n)'
*
                  CALL DGEMV( 'No transpose', I-1, N-I+1, -TAU( I ),
     $                        V( 1, I ), LDV, V( I, I ), LDV, ZERO,
     $                        T( 1, I ), 1 )
               END IF
               V( I, I ) = VII
*
*              T(1:i-1,i) := T(1:i-1,1:i-1) * T(1:i-1,i)
*
               CALL DTRMV( 'Upper', 'No transpose', 'Non-unit', I-1, T,
     $                     LDT, T( 1, I ), 1 )
               T( I, I ) = TAU( I )
            END IF
   20    CONTINUE
      ELSE
         DO 40 I = K, 1, -1
            IF( TAU( I ).EQ.ZERO ) THEN
*
*              H(i)  =  I
*
               DO 30 J = I, K
                  T( J, I ) = ZERO
   30          CONTINUE
            ELSE
*
*              general case
*
               IF( I.LT.K ) THEN
                  IF( LSAME( STOREV, 'C' ) ) THEN
                     VII = V( N-K+I, I )
                     V( N-K+I, I ) = ONE
*
*                    T(i+1:k,i) :=
*                            - tau(i) * V(1:n-k+i,i+1:k)' * V(1:n-k+i,i)
*
                     CALL DGEMV( 'Transpose', N-K+I, K-I, -TAU( I ),
     $                           V( 1, I+1 ), LDV, V( 1, I ), 1, ZERO,
     $                           T( I+1, I ), 1 )
                     V( N-K+I, I ) = VII
                  ELSE
                     VII = V( I, N-K+I )
                     V( I, N-K+I ) = ONE
*
*                    T(i+1:k,i) :=
*                            - tau(i) * V(i+1:k,1:n-k+i) * V(i,1:n-k+i)'
*
                     CALL DGEMV( 'No transpose', K-I, N-K+I, -TAU( I ),
     $                           V( I+1, 1 ), LDV, V( I, 1 ), LDV, ZERO,
     $                           T( I+1, I ), 1 )
                     V( I, N-K+I ) = VII
                  END IF
*
*                 T(i+1:k,i) := T(i+1:k,i+1:k) * T(i+1:k,i)
*
                  CALL DTRMV( 'Lower', 'No transpose', 'Non-unit', K-I,
     $                        T( I+1, I+1 ), LDT, T( I+1, I ), 1 )
               END IF
               T( I, I ) = TAU( I )
            END IF
   40    CONTINUE
      END IF
      RETURN
*
*     End of DLARFT
*
      END
      SUBROUTINE DLARTG( F, G, CS, SN, R )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
CID      DOUBLE PRECISION   CS, F, G, R, SN
*     ..
*
*  Purpose
*  =======
*
*  DLARTG generate a plane rotation so that
*
*     [  CS  SN  ]  .  [ F ]  =  [ R ]   where CS**2 + SN**2 = 1.
*     [ -SN  CS  ]     [ G ]     [ 0 ]
*
*  This is a slower, more accurate version of the BLAS1 routine DROTG,
*  with the following other differences:
*     F and G are unchanged on return.
*     If G=0, then CS=1 and SN=0.
*     If F=0 and (G .ne. 0), then CS=0 and SN=1 without doing any
*        floating point operations (saves work in DBDSQR when
*        there are zeros on the diagonal).
*
*  If F exceeds G in magnitude, CS will be positive.
*
*  Arguments
*  =========
*
*  F       (input) DOUBLE PRECISION
*          The first component of vector to be rotated.
*
*  G       (input) DOUBLE PRECISION
*          The second component of vector to be rotated.
*
*  CS      (output) DOUBLE PRECISION
*          The cosine of the rotation.
*
*  SN      (output) DOUBLE PRECISION
*          The sine of the rotation.
*
*  R       (output) DOUBLE PRECISION
*          The nonzero component of the rotated vector.
*
*  This version has a few statements commented out for thread safety
*  (machine parameters are computed on each entry). 10 feb 03, SJH.
*
*  =====================================================================
*
*     .. Parameters ..
CID      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
CID      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
CID      DOUBLE PRECISION   TWO
      PARAMETER          ( TWO = 2.0D0 )
*     ..
*     .. Local Scalars ..
*     LOGICAL            FIRST
      INTEGER            COUNT, I
CID      DOUBLE PRECISION   EPS, F1, G1, SAFMIN, SAFMN2, SAFMX2, SCALE
*     ..
*     .. External Functions ..
CID      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, INT, LOG, MAX, SQRT
*     ..
*     .. Save statement ..
*     SAVE               FIRST, SAFMX2, SAFMIN, SAFMN2
*     ..
*     .. Data statements ..
*     DATA               FIRST / .TRUE. /
*     ..
*     .. Executable Statements ..
*
*     IF( FIRST ) THEN
         SAFMIN = DLAMCH( 'S' )
         EPS = DLAMCH( 'E' )
         SAFMN2 = DLAMCH( 'B' )**INT( LOG( SAFMIN / EPS ) /
     $            LOG( DLAMCH( 'B' ) ) / TWO )
         SAFMX2 = ONE / SAFMN2
*        FIRST = .FALSE.
*     END IF
      IF( G.EQ.ZERO ) THEN
         CS = ONE
         SN = ZERO
         R = F
      ELSE IF( F.EQ.ZERO ) THEN
         CS = ZERO
         SN = ONE
         R = G
      ELSE
         F1 = F
         G1 = G
         SCALE = MAX( ABS( F1 ), ABS( G1 ) )
         IF( SCALE.GE.SAFMX2 ) THEN
            COUNT = 0
   10       CONTINUE
            COUNT = COUNT + 1
            F1 = F1*SAFMN2
            G1 = G1*SAFMN2
            SCALE = MAX( ABS( F1 ), ABS( G1 ) )
            IF( SCALE.GE.SAFMX2 )
     $         GO TO 10
            R = SQRT( F1**2+G1**2 )
            CS = F1 / R
            SN = G1 / R
            DO 20 I = 1, COUNT
               R = R*SAFMX2
   20       CONTINUE
         ELSE IF( SCALE.LE.SAFMN2 ) THEN
            COUNT = 0
   30       CONTINUE
            COUNT = COUNT + 1
            F1 = F1*SAFMX2
            G1 = G1*SAFMX2
            SCALE = MAX( ABS( F1 ), ABS( G1 ) )
            IF( SCALE.LE.SAFMN2 )
     $         GO TO 30
            R = SQRT( F1**2+G1**2 )
            CS = F1 / R
            SN = G1 / R
            DO 40 I = 1, COUNT
               R = R*SAFMN2
   40       CONTINUE
         ELSE
            R = SQRT( F1**2+G1**2 )
            CS = F1 / R
            SN = G1 / R
         END IF
         IF( ABS( F ).GT.ABS( G ) .AND. CS.LT.ZERO ) THEN
            CS = -CS
            SN = -SN
            R = -R
         END IF
      END IF
      RETURN
*
*     End of DLARTG
*
      END
      SUBROUTINE DLAS2( F, G, H, SSMIN, SSMAX )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
CID      DOUBLE PRECISION   F, G, H, SSMAX, SSMIN
*     ..
*
*  Purpose
*  =======
*
*  DLAS2  computes the singular values of the 2-by-2 matrix
*     [  F   G  ]
*     [  0   H  ].
*  On return, SSMIN is the smaller singular value and SSMAX is the
*  larger singular value.
*
*  Arguments
*  =========
*
*  F       (input) DOUBLE PRECISION
*          The (1,1) element of the 2-by-2 matrix.
*
*  G       (input) DOUBLE PRECISION
*          The (1,2) element of the 2-by-2 matrix.
*
*  H       (input) DOUBLE PRECISION
*          The (2,2) element of the 2-by-2 matrix.
*
*  SSMIN   (output) DOUBLE PRECISION
*          The smaller singular value.
*
*  SSMAX   (output) DOUBLE PRECISION
*          The larger singular value.
*
*  Further Details
*  ===============
*
*  Barring over/underflow, all output quantities are correct to within
*  a few units in the last place (ulps), even in the absence of a guard
*  digit in addition/subtraction.
*
*  In IEEE arithmetic, the code works correctly if one matrix element is
*  infinite.
*
*  Overflow will not occur unless the largest singular value itself
*  overflows, or is within a few ulps of overflow. (On machines with
*  partial overflow, like the Cray, overflow may occur if the largest
*  singular value is within a factor of 2 of overflow.)
*
*  Underflow is harmless if underflow is gradual. Otherwise, results
*  may correspond to a matrix modified by perturbations of size near
*  the underflow threshold.
*
*  ====================================================================
*
*     .. Parameters ..
CID      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
CID      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
CID      DOUBLE PRECISION   TWO
      PARAMETER          ( TWO = 2.0D0 )
*     ..
*     .. Local Scalars ..
CID      DOUBLE PRECISION   AS, AT, AU, C, FA, FHMN, FHMX, GA, HA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
*
      FA = ABS( F )
      GA = ABS( G )
      HA = ABS( H )
      FHMN = MIN( FA, HA )
      FHMX = MAX( FA, HA )
      IF( FHMN.EQ.ZERO ) THEN
         SSMIN = ZERO
         IF( FHMX.EQ.ZERO ) THEN
            SSMAX = GA
         ELSE
            SSMAX = MAX( FHMX, GA )*SQRT( ONE+
     $              ( MIN( FHMX, GA ) / MAX( FHMX, GA ) )**2 )
         END IF
      ELSE
         IF( GA.LT.FHMX ) THEN
            AS = ONE + FHMN / FHMX
            AT = ( FHMX-FHMN ) / FHMX
            AU = ( GA / FHMX )**2
            C = TWO / ( SQRT( AS*AS+AU )+SQRT( AT*AT+AU ) )
            SSMIN = FHMN*C
            SSMAX = FHMX / C
         ELSE
            AU = FHMX / GA
            IF( AU.EQ.ZERO ) THEN
*
*              Avoid possible harmful underflow if exponent range
*              asymmetric (true SSMIN may not underflow even if
*              AU underflows)
*
               SSMIN = ( FHMN*FHMX ) / GA
               SSMAX = GA
            ELSE
               AS = ONE + FHMN / FHMX
               AT = ( FHMX-FHMN ) / FHMX
               C = ONE / ( SQRT( ONE+( AS*AU )**2 )+
     $             SQRT( ONE+( AT*AU )**2 ) )
               SSMIN = ( FHMN*C )*AU
               SSMIN = SSMIN + SSMIN
               SSMAX = GA / ( C+C )
            END IF
         END IF
      END IF
      RETURN
*
*     End of DLAS2
*
      END
      SUBROUTINE DLASCL( TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          TYPE
      INTEGER            INFO, KL, KU, LDA, M, N
CID      DOUBLE PRECISION   CFROM, CTO
*     ..
*     .. Array Arguments ..
CID      DOUBLE PRECISION   A( LDA, * )
      DIMENSION A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  DLASCL multiplies the M by N real matrix A by the real scalar
*  CTO/CFROM.  This is done without over/underflow as long as the final
*  result CTO*A(I,J)/CFROM does not over/underflow. TYPE specifies that
*  A may be full, upper triangular, lower triangular, upper Hessenberg,
*  or banded.
*
*  Arguments
*  =========
*
*  TYPE    (input) CHARACTER*1
*          TYPE indices the storage type of the input matrix.
*          = 'G':  A is a full matrix.
*          = 'L':  A is a lower triangular matrix.
*          = 'U':  A is an upper triangular matrix.
*          = 'H':  A is an upper Hessenberg matrix.
*          = 'B':  A is a symmetric band matrix with lower bandwidth KL
*                  and upper bandwidth KU and with the only the lower
*                  half stored.
*          = 'Q':  A is a symmetric band matrix with lower bandwidth KL
*                  and upper bandwidth KU and with the only the upper
*                  half stored.
*          = 'Z':  A is a band matrix with lower bandwidth KL and upper
*                  bandwidth KU.
*
*  KL      (input) INTEGER
*          The lower bandwidth of A.  Referenced only if TYPE = 'B',
*          'Q' or 'Z'.
*
*  KU      (input) INTEGER
*          The upper bandwidth of A.  Referenced only if TYPE = 'B',
*          'Q' or 'Z'.
*
*  CFROM   (input) DOUBLE PRECISION
*  CTO     (input) DOUBLE PRECISION
*          The matrix A is multiplied by CTO/CFROM. A(I,J) is computed
*          without over/underflow if the final result CTO*A(I,J)/CFROM
*          can be represented without over/underflow.  CFROM must be
*          nonzero.
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          The matrix to be multiplied by CTO/CFROM.  See TYPE for the
*          storage type.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  INFO    (output) INTEGER
*          0  - successful exit
*          <0 - if INFO = -i, the i-th argument had an illegal value.
*
*  =====================================================================
*
*     .. Parameters ..
CID      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            DONE
      INTEGER            I, ITYPE, J, K1, K2, K3, K4
CID      DOUBLE PRECISION   BIGNUM, CFROM1, CFROMC, CTO1, CTOC, MUL, SMLNUM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
CID      DOUBLE PRECISION   DLAMCH
      EXTERNAL           LSAME, DLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
*
      IF( LSAME( TYPE, 'G' ) ) THEN
         ITYPE = 0
      ELSE IF( LSAME( TYPE, 'L' ) ) THEN
         ITYPE = 1
      ELSE IF( LSAME( TYPE, 'U' ) ) THEN
         ITYPE = 2
      ELSE IF( LSAME( TYPE, 'H' ) ) THEN
         ITYPE = 3
      ELSE IF( LSAME( TYPE, 'B' ) ) THEN
         ITYPE = 4
      ELSE IF( LSAME( TYPE, 'Q' ) ) THEN
         ITYPE = 5
      ELSE IF( LSAME( TYPE, 'Z' ) ) THEN
         ITYPE = 6
      ELSE
         ITYPE = -1
      END IF
*
      IF( ITYPE.EQ.-1 ) THEN
         INFO = -1
      ELSE IF( CFROM.EQ.ZERO ) THEN
         INFO = -4
      ELSE IF( M.LT.0 ) THEN
         INFO = -6
      ELSE IF( N.LT.0 .OR. ( ITYPE.EQ.4 .AND. N.NE.M ) .OR.
     $         ( ITYPE.EQ.5 .AND. N.NE.M ) ) THEN
         INFO = -7
      ELSE IF( ITYPE.LE.3 .AND. LDA.LT.MAX( 1, M ) ) THEN
         INFO = -9
      ELSE IF( ITYPE.GE.4 ) THEN
         IF( KL.LT.0 .OR. KL.GT.MAX( M-1, 0 ) ) THEN
            INFO = -2
         ELSE IF( KU.LT.0 .OR. KU.GT.MAX( N-1, 0 ) .OR.
     $            ( ( ITYPE.EQ.4 .OR. ITYPE.EQ.5 ) .AND. KL.NE.KU ) )
     $             THEN
            INFO = -3
         ELSE IF( ( ITYPE.EQ.4 .AND. LDA.LT.KL+1 ) .OR.
     $            ( ITYPE.EQ.5 .AND. LDA.LT.KU+1 ) .OR.
     $            ( ITYPE.EQ.6 .AND. LDA.LT.2*KL+KU+1 ) ) THEN
            INFO = -9
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLASCL', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. M.EQ.0 )
     $   RETURN
*
*     Get machine parameters
*
      SMLNUM = DLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
*
      CFROMC = CFROM
      CTOC = CTO
*
   10 CONTINUE
      CFROM1 = CFROMC*SMLNUM
      CTO1 = CTOC / BIGNUM
      IF( ABS( CFROM1 ).GT.ABS( CTOC ) .AND. CTOC.NE.ZERO ) THEN
         MUL = SMLNUM
         DONE = .FALSE.
         CFROMC = CFROM1
      ELSE IF( ABS( CTO1 ).GT.ABS( CFROMC ) ) THEN
         MUL = BIGNUM
         DONE = .FALSE.
         CTOC = CTO1
      ELSE
         MUL = CTOC / CFROMC
         DONE = .TRUE.
      END IF
*
      IF( ITYPE.EQ.0 ) THEN
*
*        Full matrix
*
         DO 30 J = 1, N
            DO 20 I = 1, M
               A( I, J ) = A( I, J )*MUL
   20       CONTINUE
   30    CONTINUE
*
      ELSE IF( ITYPE.EQ.1 ) THEN
*
*        Lower triangular matrix
*
         DO 50 J = 1, N
            DO 40 I = J, M
               A( I, J ) = A( I, J )*MUL
   40       CONTINUE
   50    CONTINUE
*
      ELSE IF( ITYPE.EQ.2 ) THEN
*
*        Upper triangular matrix
*
         DO 70 J = 1, N
            DO 60 I = 1, MIN( J, M )
               A( I, J ) = A( I, J )*MUL
   60       CONTINUE
   70    CONTINUE
*
      ELSE IF( ITYPE.EQ.3 ) THEN
*
*        Upper Hessenberg matrix
*
         DO 90 J = 1, N
            DO 80 I = 1, MIN( J+1, M )
               A( I, J ) = A( I, J )*MUL
   80       CONTINUE
   90    CONTINUE
*
      ELSE IF( ITYPE.EQ.4 ) THEN
*
*        Lower half of a symmetric band matrix
*
         K3 = KL + 1
         K4 = N + 1
         DO 110 J = 1, N
            DO 100 I = 1, MIN( K3, K4-J )
               A( I, J ) = A( I, J )*MUL
  100       CONTINUE
  110    CONTINUE
*
      ELSE IF( ITYPE.EQ.5 ) THEN
*
*        Upper half of a symmetric band matrix
*
         K1 = KU + 2
         K3 = KU + 1
         DO 130 J = 1, N
            DO 120 I = MAX( K1-J, 1 ), K3
               A( I, J ) = A( I, J )*MUL
  120       CONTINUE
  130    CONTINUE
*
      ELSE IF( ITYPE.EQ.6 ) THEN
*
*        Band matrix
*
         K1 = KL + KU + 2
         K2 = KL + 1
         K3 = 2*KL + KU + 1
         K4 = KL + KU + 1 + M
         DO 150 J = 1, N
            DO 140 I = MAX( K1-J, K2 ), MIN( K3, K4-J )
               A( I, J ) = A( I, J )*MUL
  140       CONTINUE
  150    CONTINUE
*
      END IF
*
      IF( .NOT.DONE )
     $   GO TO 10
*
      RETURN
*
*     End of DLASCL
*
      END
      SUBROUTINE DLASET( UPLO, M, N, ALPHA, BETA, A, LDA )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, M, N
CID      DOUBLE PRECISION   ALPHA, BETA
*     ..
*     .. Array Arguments ..
CID      DOUBLE PRECISION   A( LDA, * )
      DIMENSION A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  DLASET initializes an m-by-n matrix A to BETA on the diagonal and
*  ALPHA on the offdiagonals.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          Specifies the part of the matrix A to be set.
*          = 'U':      Upper triangular part is set; the strictly lower
*                      triangular part of A is not changed.
*          = 'L':      Lower triangular part is set; the strictly upper
*                      triangular part of A is not changed.
*          Otherwise:  All of the matrix A is set.
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  ALPHA   (input) DOUBLE PRECISION
*          The constant to which the offdiagonal elements are to be set.
*
*  BETA    (input) DOUBLE PRECISION
*          The constant to which the diagonal elements are to be set.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On exit, the leading m-by-n submatrix of A is set as follows:
*
*          if UPLO = 'U', A(i,j) = ALPHA, 1<=i<=j-1, 1<=j<=n,
*          if UPLO = 'L', A(i,j) = ALPHA, j+1<=i<=m, 1<=j<=n,
*          otherwise,     A(i,j) = ALPHA, 1<=i<=m, 1<=j<=n, i.ne.j,
*
*          and, for all UPLO, A(i,i) = BETA, 1<=i<=min(m,n).
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
* =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, J
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN
*     ..
*     .. Executable Statements ..
*
      IF( LSAME( UPLO, 'U' ) ) THEN
*
*        Set the strictly upper triangular or trapezoidal part of the
*        array to ALPHA.
*
         DO 20 J = 2, N
            DO 10 I = 1, MIN( J-1, M )
               A( I, J ) = ALPHA
   10       CONTINUE
   20    CONTINUE
*
      ELSE IF( LSAME( UPLO, 'L' ) ) THEN
*
*        Set the strictly lower triangular or trapezoidal part of the
*        array to ALPHA.
*
         DO 40 J = 1, MIN( M, N )
            DO 30 I = J + 1, M
               A( I, J ) = ALPHA
   30       CONTINUE
   40    CONTINUE
*
      ELSE
*
*        Set the leading m-by-n submatrix to ALPHA.
*
         DO 60 J = 1, N
            DO 50 I = 1, M
               A( I, J ) = ALPHA
   50       CONTINUE
   60    CONTINUE
      END IF
*
*     Set the first min(M,N) diagonal elements to BETA.
*
      DO 70 I = 1, MIN( M, N )
         A( I, I ) = BETA
   70 CONTINUE
*
      RETURN
*
*     End of DLASET
*
      END
      SUBROUTINE DLASQ1( N, D, E, WORK, INFO )
*
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INFO, N
*     ..
*     .. Array Arguments ..
CID      DOUBLE PRECISION   D( * ), E( * ), WORK( * )
      DIMENSION D( * ), E( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DLASQ1 computes the singular values of a real N-by-N bidiagonal
*  matrix with diagonal D and off-diagonal E. The singular values
*  are computed to high relative accuracy, in the absence of
*  denormalization, underflow and overflow. The algorithm was first
*  presented in
*
*  "Accurate singular values and differential qd algorithms" by K. V.
*  Fernando and B. N. Parlett, Numer. Math., Vol-67, No. 2, pp. 191-230,
*  1994,
*
*  and the present implementation is described in "An implementation of
*  the dqds Algorithm (Positive Case)", LAPACK Working Note.
*
*  Arguments
*  =========
*
*  N     (input) INTEGER
*        The number of rows and columns in the matrix. N >= 0.
*
*  D     (input/output) DOUBLE PRECISION array, dimension (N)
*        On entry, D contains the diagonal elements of the
*        bidiagonal matrix whose SVD is desired. On normal exit,
*        D contains the singular values in decreasing order.
*
*  E     (input/output) DOUBLE PRECISION array, dimension (N)
*        On entry, elements E(1:N-1) contain the off-diagonal elements
*        of the bidiagonal matrix whose SVD is desired.
*        On exit, E is overwritten.
*
*  WORK  (workspace) DOUBLE PRECISION array, dimension (4*N)
*
*  INFO  (output) INTEGER
*        = 0: successful exit
*        < 0: if INFO = -i, the i-th argument had an illegal value
*        > 0: the algorithm failed
*             = 1, a split was marked by a positive value in E
*             = 2, current block of Z not diagonalized after 30*N
*                  iterations (in inner while loop)
*             = 3, termination criterion of outer while loop not met 
*                  (program created more than N unreduced blocks)
*
*  =====================================================================
*
*     .. Parameters ..
CID      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IINFO
CID      DOUBLE PRECISION   EPS, SCALE, SAFMIN, SIGMN, SIGMX
*     ..
*     .. External Subroutines ..
      EXTERNAL           DCOPY, DLAS2, DLASCL, DLASQ2, DLASRT, XERBLA
*     ..
*     .. External Functions ..
CID      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SQRT
*     ..
*     .. Executable Statements ..
*
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -2
         CALL XERBLA( 'DLASQ1', -INFO )
         RETURN
      ELSE IF( N.EQ.0 ) THEN
         RETURN
      ELSE IF( N.EQ.1 ) THEN
         D( 1 ) = ABS( D( 1 ) )
         RETURN
      ELSE IF( N.EQ.2 ) THEN
         CALL DLAS2( D( 1 ), E( 1 ), D( 2 ), SIGMN, SIGMX )
         D( 1 ) = SIGMX
         D( 2 ) = SIGMN
         RETURN
      END IF
*
*     Estimate the largest singular value.
*
      SIGMX = ZERO
      DO 10 I = 1, N - 1
         D( I ) = ABS( D( I ) )
         SIGMX = MAX( SIGMX, ABS( E( I ) ) )
   10 CONTINUE
      D( N ) = ABS( D( N ) )
*
*     Early return if SIGMX is zero (matrix is already diagonal).
*
      IF( SIGMX.EQ.ZERO ) THEN
         CALL DLASRT( 'D', N, D, IINFO )
         RETURN
      END IF
*
      DO 20 I = 1, N
         SIGMX = MAX( SIGMX, D( I ) )
   20 CONTINUE
*
*     Copy D and E into WORK (in the Z format) and scale (squaring the
*     input data makes scaling by a power of the radix pointless).
*
      EPS = DLAMCH( 'Precision' )
      SAFMIN = DLAMCH( 'Safe minimum' )
      SCALE = SQRT( EPS / SAFMIN )
      CALL DCOPY( N, D, 1, WORK( 1 ), 2 )
      CALL DCOPY( N-1, E, 1, WORK( 2 ), 2 )
      CALL DLASCL( 'G', 0, 0, SIGMX, SCALE, 2*N-1, 1, WORK, 2*N-1,
     $             IINFO )
*         
*     Compute the q's and e's.
*
      DO 30 I = 1, 2*N - 1
         WORK( I ) = WORK( I )**2
   30 CONTINUE
      WORK( 2*N ) = ZERO
*
      CALL DLASQ2( N, WORK, INFO )
*
      IF( INFO.EQ.0 ) THEN
         DO 40 I = 1, N
            D( I ) = SQRT( WORK( I ) )
   40    CONTINUE
         CALL DLASCL( 'G', 0, 0, SCALE, SIGMX, N, 1, D, N, IINFO )
      END IF
*
      RETURN
*
*     End of DLASQ1
*
      END
      SUBROUTINE DLASQ2( N, Z, INFO )
*
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     Modified to call DLAZQ3 in place of DLASQ3, 13 Feb 03, SJH.
*
*     .. Scalar Arguments ..
      INTEGER            INFO, N
*     ..
*     .. Array Arguments ..
CID      DOUBLE PRECISION   Z( * )
      DIMENSION  Z( * )
*     ..
*
*  Purpose
*  =======
*
*  DLASQ2 computes all the eigenvalues of the symmetric positive 
*  definite tridiagonal matrix associated with the qd array Z to high
*  relative accuracy are computed to high relative accuracy, in the
*  absence of denormalization, underflow and overflow.
*
*  To see the relation of Z to the tridiagonal matrix, let L be a
*  unit lower bidiagonal matrix with subdiagonals Z(2,4,6,,..) and
*  let U be an upper bidiagonal matrix with 1's above and diagonal
*  Z(1,3,5,,..). The tridiagonal is L*U or, if you prefer, the
*  symmetric tridiagonal to which it is similar.
*
*  Note : DLASQ2 defines a logical variable, IEEE, which is true
*  on machines which follow ieee-754 floating-point standard in their
*  handling of infinities and NaNs, and false otherwise. This variable
*  is passed to DLAZQ3.
*
*  Arguments
*  =========
*
*  N     (input) INTEGER
*        The number of rows and columns in the matrix. N >= 0.
*
*  Z     (workspace) DOUBLE PRECISION array, dimension ( 4*N )
*        On entry Z holds the qd array. On exit, entries 1 to N hold
*        the eigenvalues in decreasing order, Z( 2*N+1 ) holds the
*        trace, and Z( 2*N+2 ) holds the sum of the eigenvalues. If
*        N > 2, then Z( 2*N+3 ) holds the iteration count, Z( 2*N+4 )
*        holds NDIVS/NIN^2, and Z( 2*N+5 ) holds the percentage of
*        shifts that failed.
*
*  INFO  (output) INTEGER
*        = 0: successful exit
*        < 0: if the i-th argument is a scalar and had an illegal
*             value, then INFO = -i, if the i-th argument is an
*             array and the j-entry had an illegal value, then
*             INFO = -(i*100+j)
*        > 0: the algorithm failed
*              = 1, a split was marked by a positive value in E
*              = 2, current block of Z not diagonalized after 30*N
*                   iterations (in inner while loop)
*              = 3, termination criterion of outer while loop not met 
*                   (program created more than N unreduced blocks)
*
*  Further Details
*  ===============
*  Local Variables: I0:N0 defines a current unreduced segment of Z.
*  The shifts are accumulated in SIGMA. Iteration count is in ITER.
*  Ping-pong is controlled by PP (alternates between 0 and 1).
*
*  =====================================================================
*
*     .. Parameters ..
CID      DOUBLE PRECISION   CBIAS
      PARAMETER          ( CBIAS = 1.50D0 )
CID      DOUBLE PRECISION   ZERO, HALF, ONE, TWO, FOUR, HUNDRD
      PARAMETER          ( ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0,
     $                     TWO = 2.0D0, FOUR = 4.0D0, HUNDRD = 100.0D0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            IEEE
      INTEGER            I0, I4, IINFO, IPN4, ITER, IWHILA, IWHILB, K, 
     $                   N0, NBIG, NDIV, NFAIL, PP, SPLT, TTYPE
CID      DOUBLE PRECISION   D, DESIG, DMIN, DMIN1, DMIN2, DN, DN1, DN2, E,
CID     $                   EMAX, EMIN, EPS, OLDEMN, QMAX, QMIN, S, SAFMIN,
CID     $                   SIGMA, T, TAU, TEMP, TOL, TOL2, TRACE, ZMAX
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLAZQ3, DLASRT, XERBLA
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
CID      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH, ILAENV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
*      
*     Test the input arguments.
*     (in case DLASQ2 is not called by DLASQ1)
*
      INFO = 0
      EPS = DLAMCH( 'Precision' )
      SAFMIN = DLAMCH( 'Safe minimum' )
      TOL = EPS*HUNDRD
      TOL2 = TOL**2
*
      IF( N.LT.0 ) THEN
         INFO = -1
         CALL XERBLA( 'DLASQ2', 1 )
         RETURN
      ELSE IF( N.EQ.0 ) THEN
         RETURN
      ELSE IF( N.EQ.1 ) THEN
*
*        1-by-1 case.
*
         IF( Z( 1 ).LT.ZERO ) THEN
            INFO = -201
            CALL XERBLA( 'DLASQ2', 2 )
         END IF
         RETURN
      ELSE IF( N.EQ.2 ) THEN
*
*        2-by-2 case.
*
         IF( Z( 2 ).LT.ZERO .OR. Z( 3 ).LT.ZERO ) THEN
            INFO = -2
            CALL XERBLA( 'DLASQ2', 2 )
            RETURN
         ELSE IF( Z( 3 ).GT.Z( 1 ) ) THEN
            D = Z( 3 )
            Z( 3 ) = Z( 1 )
            Z( 1 ) = D
         END IF
         Z( 5 ) = Z( 1 ) + Z( 2 ) + Z( 3 )
         IF( Z( 2 ).GT.Z( 3 )*TOL2 ) THEN
            T = HALF*( ( Z( 1 )-Z( 3 ) )+Z( 2 ) ) 
            S = Z( 3 )*( Z( 2 ) / T )
            IF( S.LE.T ) THEN
               S = Z( 3 )*( Z( 2 ) / ( T*( ONE+SQRT( ONE+S / T ) ) ) )
            ELSE
               S = Z( 3 )*( Z( 2 ) / ( T+SQRT( T )*SQRT( T+S ) ) )
            END IF
            T = Z( 1 ) + ( S+Z( 2 ) )
            Z( 3 ) = Z( 3 )*( Z( 1 ) / T )
            Z( 1 ) = T
         END IF
         Z( 2 ) = Z( 3 )
         Z( 6 ) = Z( 2 ) + Z( 1 )
         RETURN
      END IF
*
*     Check for negative data and compute sums of q's and e's.
*
      Z( 2*N ) = ZERO
      EMIN = Z( 2 )
      QMAX = ZERO
      ZMAX = ZERO
      D = ZERO
      E = ZERO
*
      DO 10 K = 1, 2*( N-1 ), 2
         IF( Z( K ).LT.ZERO ) THEN
            INFO = -( 200+K )
            CALL XERBLA( 'DLASQ2', 2 )
            RETURN
         ELSE IF( Z( K+1 ).LT.ZERO ) THEN
            INFO = -( 200+K+1 )
            CALL XERBLA( 'DLASQ2', 2 )
            RETURN
         END IF
         D = D + Z( K )
         E = E + Z( K+1 )
         QMAX = MAX( QMAX, Z( K ) )
         EMIN = MIN( EMIN, Z( K+1 ) )
         ZMAX = MAX( QMAX, ZMAX, Z( K+1 ) )
   10 CONTINUE
      IF( Z( 2*N-1 ).LT.ZERO ) THEN
         INFO = -( 200+2*N-1 )
         CALL XERBLA( 'DLASQ2', 2 )
         RETURN
      END IF
      D = D + Z( 2*N-1 )
      QMAX = MAX( QMAX, Z( 2*N-1 ) )
      ZMAX = MAX( QMAX, ZMAX )
*
*     Check for diagonality.
*
      IF( E.EQ.ZERO ) THEN
         DO 20 K = 2, N
            Z( K ) = Z( 2*K-1 )
   20    CONTINUE
         CALL DLASRT( 'D', N, Z, IINFO )
         Z( 2*N-1 ) = D
         RETURN
      END IF
*
      TRACE = D + E
*
*     Check for zero data.
*
      IF( TRACE.EQ.ZERO ) THEN
         Z( 2*N-1 ) = ZERO
         RETURN
      END IF
*         
*     Check whether the machine is IEEE conformable.
*         
      IEEE = ILAENV( 10, 'DLASQ2', 'N', 1, 2, 3, 4 ).EQ.1 .AND.
     $       ILAENV( 11, 'DLASQ2', 'N', 1, 2, 3, 4 ).EQ.1      
*         
*     Rearrange data for locality: Z=(q1,qq1,e1,ee1,q2,qq2,e2,ee2,...).
*
      DO 30 K = 2*N, 2, -2
         Z( 2*K ) = ZERO 
         Z( 2*K-1 ) = Z( K ) 
         Z( 2*K-2 ) = ZERO 
         Z( 2*K-3 ) = Z( K-1 ) 
   30 CONTINUE
*
      I0 = 1
      N0 = N
*
*     Reverse the qd-array, if warranted.
*
      IF( CBIAS*Z( 4*I0-3 ).LT.Z( 4*N0-3 ) ) THEN
         IPN4 = 4*( I0+N0 )
         DO 40 I4 = 4*I0, 2*( I0+N0-1 ), 4
            TEMP = Z( I4-3 )
            Z( I4-3 ) = Z( IPN4-I4-3 )
            Z( IPN4-I4-3 ) = TEMP
            TEMP = Z( I4-1 )
            Z( I4-1 ) = Z( IPN4-I4-5 )
            Z( IPN4-I4-5 ) = TEMP
   40    CONTINUE
      END IF
*
*     Initial split checking via dqd and Li's test.
*
      PP = 0
*
      DO 80 K = 1, 2
*
         D = Z( 4*N0+PP-3 )
         DO 50 I4 = 4*( N0-1 ) + PP, 4*I0 + PP, -4
            IF( Z( I4-1 ).LE.TOL2*D ) THEN
               Z( I4-1 ) = -ZERO
               D = Z( I4-3 )
            ELSE
               D = Z( I4-3 )*( D / ( D+Z( I4-1 ) ) )
            END IF
   50    CONTINUE
*
*        dqd maps Z to ZZ plus Li's test.
*
         EMIN = Z( 4*I0+PP+1 )
         D = Z( 4*I0+PP-3 )
         DO 60 I4 = 4*I0 + PP, 4*( N0-1 ) + PP, 4
            Z( I4-2*PP-2 ) = D + Z( I4-1 )
            IF( Z( I4-1 ).LE.TOL2*D ) THEN
               Z( I4-1 ) = -ZERO
               Z( I4-2*PP-2 ) = D
               Z( I4-2*PP ) = ZERO
               D = Z( I4+1 )
            ELSE IF( SAFMIN*Z( I4+1 ).LT.Z( I4-2*PP-2 ) .AND.
     $               SAFMIN*Z( I4-2*PP-2 ).LT.Z( I4+1 ) ) THEN
               TEMP = Z( I4+1 ) / Z( I4-2*PP-2 )
               Z( I4-2*PP ) = Z( I4-1 )*TEMP
               D = D*TEMP
            ELSE
               Z( I4-2*PP ) = Z( I4+1 )*( Z( I4-1 ) / Z( I4-2*PP-2 ) )
               D = Z( I4+1 )*( D / Z( I4-2*PP-2 ) )
            END IF
            EMIN = MIN( EMIN, Z( I4-2*PP ) )
   60    CONTINUE 
         Z( 4*N0-PP-2 ) = D
*
*        Now find qmax.
*
         QMAX = Z( 4*I0-PP-2 )
         DO 70 I4 = 4*I0 - PP + 2, 4*N0 - PP - 2, 4
            QMAX = MAX( QMAX, Z( I4 ) )
   70    CONTINUE
*
*        Prepare for the next iteration on K.
*
         PP = 1 - PP
   80 CONTINUE
*
*     Initialise variables to pass to DLAZQ3
*
      TTYPE = 0
      DMIN1 = ZERO
      DMIN2 = ZERO
      DN    = ZERO
      DN1   = ZERO
      DN2   = ZERO
      TAU   = ZERO
*
      ITER = 2
      NFAIL = 0
      NDIV = 2*( N0-I0 )
*
      DO 140 IWHILA = 1, N + 1
         IF( N0.LT.1 ) 
     $      GO TO 150
*
*        While array unfinished do 
*
*        E(N0) holds the value of SIGMA when submatrix in I0:N0
*        splits from the rest of the array, but is negated.
*      
         DESIG = ZERO
         IF( N0.EQ.N ) THEN
            SIGMA = ZERO
         ELSE
            SIGMA = -Z( 4*N0-1 )
         END IF
         IF( SIGMA.LT.ZERO ) THEN
            INFO = 1
            RETURN
         END IF
*
*        Find last unreduced submatrix's top index I0, find QMAX and
*        EMIN. Find Gershgorin-type bound if Q's much greater than E's.
*
         EMAX = ZERO 
         IF( N0.GT.I0 ) THEN
            EMIN = ABS( Z( 4*N0-5 ) )
         ELSE
            EMIN = ZERO
         END IF
         QMIN = Z( 4*N0-3 )
         QMAX = QMIN
         DO 90 I4 = 4*N0, 8, -4
            IF( Z( I4-5 ).LE.ZERO )
     $         GO TO 100
            IF( QMIN.GE.FOUR*EMAX ) THEN
               QMIN = MIN( QMIN, Z( I4-3 ) )
               EMAX = MAX( EMAX, Z( I4-5 ) )
            END IF
            QMAX = MAX( QMAX, Z( I4-7 )+Z( I4-5 ) )
            EMIN = MIN( EMIN, Z( I4-5 ) )
   90    CONTINUE
         I4 = 4 
*
  100    CONTINUE
         I0 = I4 / 4
*
*        Store EMIN for passing to DLAZQ3.
*
         Z( 4*N0-1 ) = EMIN
*
*        Put -(initial shift) into DMIN.
*
         DMIN = -MAX( ZERO, QMIN-TWO*SQRT( QMIN )*SQRT( EMAX ) )
*
*        Now I0:N0 is unreduced. PP = 0 for ping, PP = 1 for pong.
*
         PP = 0 
*
         NBIG = 30*( N0-I0+1 )
         DO 120 IWHILB = 1, NBIG
            IF( I0.GT.N0 ) 
     $         GO TO 130
*
*           While submatrix unfinished take a good dqds step.
*
            CALL DLAZQ3( I0, N0, Z, PP, DMIN, SIGMA, DESIG, QMAX, NFAIL,
     $                   ITER, NDIV, IEEE, TTYPE, DMIN1, DMIN2, DN, DN1,
     $                   DN2, TAU )
*
            PP = 1 - PP
*
*           When EMIN is very small check for splits.
*
            IF( PP.EQ.0 .AND. N0-I0.GE.3 ) THEN
               IF( Z( 4*N0 ).LE.TOL2*QMAX .OR.
     $             Z( 4*N0-1 ).LE.TOL2*SIGMA ) THEN
                  SPLT = I0 - 1
                  QMAX = Z( 4*I0-3 )
                  EMIN = Z( 4*I0-1 )
                  OLDEMN = Z( 4*I0 )
                  DO 110 I4 = 4*I0, 4*( N0-3 ), 4
                     IF( Z( I4 ).LE.TOL2*Z( I4-3 ) .OR.
     $                   Z( I4-1 ).LE.TOL2*SIGMA ) THEN
                        Z( I4-1 ) = -SIGMA
                        SPLT = I4 / 4
                        QMAX = ZERO
                        EMIN = Z( I4+3 )
                        OLDEMN = Z( I4+4 )
                     ELSE
                        QMAX = MAX( QMAX, Z( I4+1 ) )
                        EMIN = MIN( EMIN, Z( I4-1 ) )
                        OLDEMN = MIN( OLDEMN, Z( I4 ) )
                     END IF
  110             CONTINUE
                  Z( 4*N0-1 ) = EMIN
                  Z( 4*N0 ) = OLDEMN
                  I0 = SPLT + 1
               END IF
            END IF
*
  120    CONTINUE
*
         INFO = 2
         RETURN
*
*        end IWHILB
*
  130    CONTINUE
*
  140 CONTINUE
*
      INFO = 3
      RETURN
*
*     end IWHILA   
*
  150 CONTINUE
*      
*     Move q's to the front.
*      
      DO 160 K = 2, N
         Z( K ) = Z( 4*K-3 )
  160 CONTINUE
*      
*     Sort and compute sum of eigenvalues.
*
      CALL DLASRT( 'D', N, Z, IINFO )
*
      E = ZERO
      DO 170 K = N, 1, -1
         E = E + Z( K )
  170 CONTINUE
*
*     Store trace, sum(eigenvalues) and information on performance.
*
      Z( 2*N+1 ) = TRACE 
      Z( 2*N+2 ) = E
      Z( 2*N+3 ) = DBLE( ITER )
      Z( 2*N+4 ) = DBLE( NDIV ) / DBLE( N**2 )
      Z( 2*N+5 ) = HUNDRD*NFAIL / DBLE( ITER )
      RETURN
*
*     End of DLASQ2
*
      END
      SUBROUTINE DLASQ5( I0, N0, Z, PP, TAU, DMIN, DMIN1, DMIN2, DN,
     $                   DNM1, DNM2, IEEE )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      LOGICAL            IEEE
      INTEGER            I0, N0, PP
CID      DOUBLE PRECISION   DMIN, DMIN1, DMIN2, DN, DNM1, DNM2, TAU
*     ..
*     .. Array Arguments ..
CID      DOUBLE PRECISION   Z( * )
      DIMENSION    Z( * )
*     ..
*
*  Purpose
*  =======
*
*  DLASQ5 computes one dqds transform in ping-pong form, one
*  version for IEEE machines another for non IEEE machines.
*
*  Arguments
*  =========
*
*  I0    (input) INTEGER
*        First index.
*
*  N0    (input) INTEGER
*        Last index.
*
*  Z     (input) DOUBLE PRECISION array, dimension ( 4*N )
*        Z holds the qd array. EMIN is stored in Z(4*N0) to avoid
*        an extra argument.
*
*  PP    (input) INTEGER
*        PP=0 for ping, PP=1 for pong.
*
*  TAU   (input) DOUBLE PRECISION
*        This is the shift.
*
*  DMIN  (output) DOUBLE PRECISION
*        Minimum value of d.
*
*  DMIN1 (output) DOUBLE PRECISION
*        Minimum value of d, excluding D( N0 ).
*
*  DMIN2 (output) DOUBLE PRECISION
*        Minimum value of d, excluding D( N0 ) and D( N0-1 ).
*
*  DN    (output) DOUBLE PRECISION
*        d(N0), the last value of d.
*
*  DNM1  (output) DOUBLE PRECISION
*        d(N0-1).
*
*  DNM2  (output) DOUBLE PRECISION
*        d(N0-2).
*
*  IEEE  (input) LOGICAL
*        Flag for IEEE or non IEEE arithmetic.
*
*  =====================================================================
*
*     .. Parameter ..
CID      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
*     ..
*     .. Local Scalars ..
      INTEGER            J4, J4P2
CID      DOUBLE PRECISION   D, EMIN, TEMP
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN
*     ..
*     .. Executable Statements ..
*
      IF( ( N0-I0-1 ).LE.0 )
     $   RETURN
*
      J4 = 4*I0 + PP - 3
      EMIN = Z( J4+4 )
      D = Z( J4 ) - TAU
      DMIN = D
      DMIN1 = -Z( J4 )
*
      IF( IEEE ) THEN
*
*        Code for IEEE arithmetic.
*
         IF( PP.EQ.0 ) THEN
            DO 10 J4 = 4*I0, 4*( N0-3 ), 4
               Z( J4-2 ) = D + Z( J4-1 )
               TEMP = Z( J4+1 ) / Z( J4-2 )
               D = D*TEMP - TAU
               DMIN = MIN( DMIN, D )
               Z( J4 ) = Z( J4-1 )*TEMP
               EMIN = MIN( Z( J4 ), EMIN )
   10       CONTINUE
         ELSE
            DO 20 J4 = 4*I0, 4*( N0-3 ), 4
               Z( J4-3 ) = D + Z( J4 )
               TEMP = Z( J4+2 ) / Z( J4-3 )
               D = D*TEMP - TAU
               DMIN = MIN( DMIN, D )
               Z( J4-1 ) = Z( J4 )*TEMP
               EMIN = MIN( Z( J4-1 ), EMIN )
   20       CONTINUE
         END IF
*
*        Unroll last two steps.
*
         DNM2 = D
         DMIN2 = DMIN
         J4 = 4*( N0-2 ) - PP
         J4P2 = J4 + 2*PP - 1
         Z( J4-2 ) = DNM2 + Z( J4P2 )
         Z( J4 ) = Z( J4P2+2 )*( Z( J4P2 ) / Z( J4-2 ) )
         DNM1 = Z( J4P2+2 )*( DNM2 / Z( J4-2 ) ) - TAU
         DMIN = MIN( DMIN, DNM1 )
*
         DMIN1 = DMIN
         J4 = J4 + 4
         J4P2 = J4 + 2*PP - 1
         Z( J4-2 ) = DNM1 + Z( J4P2 )
         Z( J4 ) = Z( J4P2+2 )*( Z( J4P2 ) / Z( J4-2 ) )
         DN = Z( J4P2+2 )*( DNM1 / Z( J4-2 ) ) - TAU
         DMIN = MIN( DMIN, DN )
*
      ELSE
*
*        Code for non IEEE arithmetic.
*
         IF( PP.EQ.0 ) THEN
            DO 30 J4 = 4*I0, 4*( N0-3 ), 4
               Z( J4-2 ) = D + Z( J4-1 )
               IF( D.LT.ZERO ) THEN
                  RETURN
               ELSE
                  Z( J4 ) = Z( J4+1 )*( Z( J4-1 ) / Z( J4-2 ) )
                  D = Z( J4+1 )*( D / Z( J4-2 ) ) - TAU
               END IF
               DMIN = MIN( DMIN, D )
               EMIN = MIN( EMIN, Z( J4 ) )
   30       CONTINUE
         ELSE
            DO 40 J4 = 4*I0, 4*( N0-3 ), 4
               Z( J4-3 ) = D + Z( J4 )
               IF( D.LT.ZERO ) THEN
                  RETURN
               ELSE
                  Z( J4-1 ) = Z( J4+2 )*( Z( J4 ) / Z( J4-3 ) )
                  D = Z( J4+2 )*( D / Z( J4-3 ) ) - TAU
               END IF
               DMIN = MIN( DMIN, D )
               EMIN = MIN( EMIN, Z( J4-1 ) )
   40       CONTINUE
         END IF
*
*        Unroll last two steps.
*
         DNM2 = D
         DMIN2 = DMIN
         J4 = 4*( N0-2 ) - PP
         J4P2 = J4 + 2*PP - 1
         Z( J4-2 ) = DNM2 + Z( J4P2 )
         IF( DNM2.LT.ZERO ) THEN
            RETURN
         ELSE
            Z( J4 ) = Z( J4P2+2 )*( Z( J4P2 ) / Z( J4-2 ) )
            DNM1 = Z( J4P2+2 )*( DNM2 / Z( J4-2 ) ) - TAU
         END IF
         DMIN = MIN( DMIN, DNM1 )
*
         DMIN1 = DMIN
         J4 = J4 + 4
         J4P2 = J4 + 2*PP - 1
         Z( J4-2 ) = DNM1 + Z( J4P2 )
         IF( DNM1.LT.ZERO ) THEN
            RETURN
         ELSE
            Z( J4 ) = Z( J4P2+2 )*( Z( J4P2 ) / Z( J4-2 ) )
            DN = Z( J4P2+2 )*( DNM1 / Z( J4-2 ) ) - TAU
         END IF
         DMIN = MIN( DMIN, DN )
*
      END IF
*
      Z( J4+2 ) = DN
      Z( 4*N0-PP ) = EMIN
      RETURN
*
*     End of DLASQ5
*
      END
      SUBROUTINE DLASQ6( I0, N0, Z, PP, DMIN, DMIN1, DMIN2, DN,
     $                   DNM1, DNM2 )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            I0, N0, PP
CID      DOUBLE PRECISION   DMIN, DMIN1, DMIN2, DN, DNM1, DNM2
*     ..
*     .. Array Arguments ..
CID      DOUBLE PRECISION   Z( * )
      DIMENSION    Z( * )
*     ..
*
*  Purpose
*  =======
*
*  DLASQ6 computes one dqd (shift equal to zero) transform in
*  ping-pong form, with protection against underflow and overflow.
*
*  Arguments
*  =========
*
*  I0    (input) INTEGER
*        First index.
*
*  N0    (input) INTEGER
*        Last index.
*
*  Z     (input) DOUBLE PRECISION array, dimension ( 4*N )
*        Z holds the qd array. EMIN is stored in Z(4*N0) to avoid
*        an extra argument.
*
*  PP    (input) INTEGER
*        PP=0 for ping, PP=1 for pong.
*
*  DMIN  (output) DOUBLE PRECISION
*        Minimum value of d.
*
*  DMIN1 (output) DOUBLE PRECISION
*        Minimum value of d, excluding D( N0 ).
*
*  DMIN2 (output) DOUBLE PRECISION
*        Minimum value of d, excluding D( N0 ) and D( N0-1 ).
*
*  DN    (output) DOUBLE PRECISION
*        d(N0), the last value of d.
*
*  DNM1  (output) DOUBLE PRECISION
*        d(N0-1).
*
*  DNM2  (output) DOUBLE PRECISION
*        d(N0-2).
*
*  =====================================================================
*
*     .. Parameter ..
CID      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
*     ..
*     .. Local Scalars ..
      INTEGER            J4, J4P2
CID      DOUBLE PRECISION   D, EMIN, SAFMIN, TEMP
*     ..
*     .. External Function ..
CID      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN
*     ..
*     .. Executable Statements ..
*
      IF( ( N0-I0-1 ).LE.0 )
     $   RETURN
*
      SAFMIN = DLAMCH( 'Safe minimum' )
      J4 = 4*I0 + PP - 3
      EMIN = Z( J4+4 ) 
      D = Z( J4 )
      DMIN = D
*
      IF( PP.EQ.0 ) THEN
         DO 10 J4 = 4*I0, 4*( N0-3 ), 4
            Z( J4-2 ) = D + Z( J4-1 ) 
            IF( Z( J4-2 ).EQ.ZERO ) THEN
               Z( J4 ) = ZERO
               D = Z( J4+1 )
               DMIN = D
               EMIN = ZERO
            ELSE IF( SAFMIN*Z( J4+1 ).LT.Z( J4-2 ) .AND.
     $               SAFMIN*Z( J4-2 ).LT.Z( J4+1 ) ) THEN
               TEMP = Z( J4+1 ) / Z( J4-2 )
               Z( J4 ) = Z( J4-1 )*TEMP
               D = D*TEMP
            ELSE 
               Z( J4 ) = Z( J4+1 )*( Z( J4-1 ) / Z( J4-2 ) )
               D = Z( J4+1 )*( D / Z( J4-2 ) )
            END IF
            DMIN = MIN( DMIN, D )
            EMIN = MIN( EMIN, Z( J4 ) )
   10    CONTINUE
      ELSE
         DO 20 J4 = 4*I0, 4*( N0-3 ), 4
            Z( J4-3 ) = D + Z( J4 ) 
            IF( Z( J4-3 ).EQ.ZERO ) THEN
               Z( J4-1 ) = ZERO
               D = Z( J4+2 )
               DMIN = D
               EMIN = ZERO
            ELSE IF( SAFMIN*Z( J4+2 ).LT.Z( J4-3 ) .AND.
     $               SAFMIN*Z( J4-3 ).LT.Z( J4+2 ) ) THEN
               TEMP = Z( J4+2 ) / Z( J4-3 )
               Z( J4-1 ) = Z( J4 )*TEMP
               D = D*TEMP
            ELSE 
               Z( J4-1 ) = Z( J4+2 )*( Z( J4 ) / Z( J4-3 ) )
               D = Z( J4+2 )*( D / Z( J4-3 ) )
            END IF
            DMIN = MIN( DMIN, D )
            EMIN = MIN( EMIN, Z( J4-1 ) )
   20    CONTINUE
      END IF
*
*     Unroll last two steps. 
*
      DNM2 = D
      DMIN2 = DMIN
      J4 = 4*( N0-2 ) - PP
      J4P2 = J4 + 2*PP - 1
      Z( J4-2 ) = DNM2 + Z( J4P2 )
      IF( Z( J4-2 ).EQ.ZERO ) THEN
         Z( J4 ) = ZERO
         DNM1 = Z( J4P2+2 )
         DMIN = DNM1
         EMIN = ZERO
      ELSE IF( SAFMIN*Z( J4P2+2 ).LT.Z( J4-2 ) .AND.
     $         SAFMIN*Z( J4-2 ).LT.Z( J4P2+2 ) ) THEN
         TEMP = Z( J4P2+2 ) / Z( J4-2 )
         Z( J4 ) = Z( J4P2 )*TEMP
         DNM1 = DNM2*TEMP
      ELSE
         Z( J4 ) = Z( J4P2+2 )*( Z( J4P2 ) / Z( J4-2 ) )
         DNM1 = Z( J4P2+2 )*( DNM2 / Z( J4-2 ) )
      END IF
      DMIN = MIN( DMIN, DNM1 )
*
      DMIN1 = DMIN
      J4 = J4 + 4
      J4P2 = J4 + 2*PP - 1
      Z( J4-2 ) = DNM1 + Z( J4P2 )
      IF( Z( J4-2 ).EQ.ZERO ) THEN
         Z( J4 ) = ZERO
         DN = Z( J4P2+2 )
         DMIN = DN
         EMIN = ZERO
      ELSE IF( SAFMIN*Z( J4P2+2 ).LT.Z( J4-2 ) .AND.
     $         SAFMIN*Z( J4-2 ).LT.Z( J4P2+2 ) ) THEN
         TEMP = Z( J4P2+2 ) / Z( J4-2 )
         Z( J4 ) = Z( J4P2 )*TEMP
         DN = DNM1*TEMP
      ELSE
         Z( J4 ) = Z( J4P2+2 )*( Z( J4P2 ) / Z( J4-2 ) )
         DN = Z( J4P2+2 )*( DNM1 / Z( J4-2 ) )
      END IF
      DMIN = MIN( DMIN, DN )
*
      Z( J4+2 ) = DN
      Z( 4*N0-PP ) = EMIN
      RETURN
*
*     End of DLASQ6
*
      END
      SUBROUTINE DLASR( SIDE, PIVOT, DIRECT, M, N, C, S, A, LDA )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          DIRECT, PIVOT, SIDE
      INTEGER            LDA, M, N
*     ..
*     .. Array Arguments ..
CID      DOUBLE PRECISION   A( LDA, * ), C( * ), S( * )
      DIMENSION    A( LDA, * ), C( * ), S( * )
*     ..
*
*  Purpose
*  =======
*
*  DLASR applies a sequence of plane rotations to a real matrix A,
*  from either the left or the right.
*  
*  When SIDE = 'L', the transformation takes the form
*  
*     A := P*A
*  
*  and when SIDE = 'R', the transformation takes the form
*  
*     A := A*P**T
*  
*  where P is an orthogonal matrix consisting of a sequence of z plane
*  rotations, with z = M when SIDE = 'L' and z = N when SIDE = 'R',
*  and P**T is the transpose of P.
*  
*  When DIRECT = 'F' (Forward sequence), then
*  
*     P = P(z-1) * ... * P(2) * P(1)
*  
*  and when DIRECT = 'B' (Backward sequence), then
*  
*     P = P(1) * P(2) * ... * P(z-1)
*  
*  where P(k) is a plane rotation matrix defined by the 2-by-2 rotation
*  
*     R(k) = (  c(k)  s(k) )
*          = ( -s(k)  c(k) ).
*  
*  When PIVOT = 'V' (Variable pivot), the rotation is performed
*  for the plane (k,k+1), i.e., P(k) has the form
*  
*     P(k) = (  1                                            )
*            (       ...                                     )
*            (              1                                )
*            (                   c(k)  s(k)                  )
*            (                  -s(k)  c(k)                  )
*            (                                1              )
*            (                                     ...       )
*            (                                            1  )
*  
*  where R(k) appears as a rank-2 modification to the identity matrix in
*  rows and columns k and k+1.
*  
*  When PIVOT = 'T' (Top pivot), the rotation is performed for the
*  plane (1,k+1), so P(k) has the form
*  
*     P(k) = (  c(k)                    s(k)                 )
*            (         1                                     )
*            (              ...                              )
*            (                     1                         )
*            ( -s(k)                    c(k)                 )
*            (                                 1             )
*            (                                      ...      )
*            (                                             1 )
*  
*  where R(k) appears in rows and columns 1 and k+1.
*  
*  Similarly, when PIVOT = 'B' (Bottom pivot), the rotation is
*  performed for the plane (k,z), giving P(k) the form
*  
*     P(k) = ( 1                                             )
*            (      ...                                      )
*            (             1                                 )
*            (                  c(k)                    s(k) )
*            (                         1                     )
*            (                              ...              )
*            (                                     1         )
*            (                 -s(k)                    c(k) )
*  
*  where R(k) appears in rows and columns k and z.  The rotations are
*  performed without ever forming P(k) explicitly.
*
*  Arguments
*  =========
*
*  SIDE    (input) CHARACTER*1
*          Specifies whether the plane rotation matrix P is applied to
*          A on the left or the right.
*          = 'L':  Left, compute A := P*A
*          = 'R':  Right, compute A:= A*P**T
*
*  PIVOT   (input) CHARACTER*1
*          Specifies the plane for which P(k) is a plane rotation
*          matrix.
*          = 'V':  Variable pivot, the plane (k,k+1)
*          = 'T':  Top pivot, the plane (1,k+1)
*          = 'B':  Bottom pivot, the plane (k,z)
*
*  DIRECT  (input) CHARACTER*1
*          Specifies whether P is a forward or backward sequence of
*          plane rotations.
*          = 'F':  Forward, P = P(z-1)*...*P(2)*P(1)
*          = 'B':  Backward, P = P(1)*P(2)*...*P(z-1)
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  If m <= 1, an immediate
*          return is effected.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  If n <= 1, an
*          immediate return is effected.
*
*  C       (input) DOUBLE PRECISION array, dimension
*                  (M-1) if SIDE = 'L'
*                  (N-1) if SIDE = 'R'
*          The cosines c(k) of the plane rotations.
*
*  S       (input) DOUBLE PRECISION array, dimension
*                  (M-1) if SIDE = 'L'
*                  (N-1) if SIDE = 'R'
*          The sines s(k) of the plane rotations.  The 2-by-2 plane
*          rotation part of the matrix P(k), R(k), has the form
*          R(k) = (  c(k)  s(k) )
*                 ( -s(k)  c(k) ).
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          The M-by-N matrix A.  On exit, A is overwritten by P*A if
*          SIDE = 'R' or by A*P**T if SIDE = 'L'.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  =====================================================================
*
*     .. Parameters ..
CID      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, INFO, J
CID      DOUBLE PRECISION   CTEMP, STEMP, TEMP
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      INFO = 0
      IF( .NOT.( LSAME( SIDE, 'L' ) .OR. LSAME( SIDE, 'R' ) ) ) THEN
         INFO = 1
      ELSE IF( .NOT.( LSAME( PIVOT, 'V' ) .OR. LSAME( PIVOT,
     $         'T' ) .OR. LSAME( PIVOT, 'B' ) ) ) THEN
         INFO = 2
      ELSE IF( .NOT.( LSAME( DIRECT, 'F' ) .OR. LSAME( DIRECT, 'B' ) ) )
     $          THEN
         INFO = 3
      ELSE IF( M.LT.0 ) THEN
         INFO = 4
      ELSE IF( N.LT.0 ) THEN
         INFO = 5
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLASR ', INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( ( M.EQ.0 ) .OR. ( N.EQ.0 ) )
     $   RETURN
      IF( LSAME( SIDE, 'L' ) ) THEN
*
*        Form  P * A
*
         IF( LSAME( PIVOT, 'V' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 20 J = 1, M - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 10 I = 1, N
                        TEMP = A( J+1, I )
                        A( J+1, I ) = CTEMP*TEMP - STEMP*A( J, I )
                        A( J, I ) = STEMP*TEMP + CTEMP*A( J, I )
   10                CONTINUE
                  END IF
   20          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 40 J = M - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 30 I = 1, N
                        TEMP = A( J+1, I )
                        A( J+1, I ) = CTEMP*TEMP - STEMP*A( J, I )
                        A( J, I ) = STEMP*TEMP + CTEMP*A( J, I )
   30                CONTINUE
                  END IF
   40          CONTINUE
            END IF
         ELSE IF( LSAME( PIVOT, 'T' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 60 J = 2, M
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 50 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = CTEMP*TEMP - STEMP*A( 1, I )
                        A( 1, I ) = STEMP*TEMP + CTEMP*A( 1, I )
   50                CONTINUE
                  END IF
   60          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 80 J = M, 2, -1
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 70 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = CTEMP*TEMP - STEMP*A( 1, I )
                        A( 1, I ) = STEMP*TEMP + CTEMP*A( 1, I )
   70                CONTINUE
                  END IF
   80          CONTINUE
            END IF
         ELSE IF( LSAME( PIVOT, 'B' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 100 J = 1, M - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 90 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = STEMP*A( M, I ) + CTEMP*TEMP
                        A( M, I ) = CTEMP*A( M, I ) - STEMP*TEMP
   90                CONTINUE
                  END IF
  100          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 120 J = M - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 110 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = STEMP*A( M, I ) + CTEMP*TEMP
                        A( M, I ) = CTEMP*A( M, I ) - STEMP*TEMP
  110                CONTINUE
                  END IF
  120          CONTINUE
            END IF
         END IF
      ELSE IF( LSAME( SIDE, 'R' ) ) THEN
*
*        Form A * P'
*
         IF( LSAME( PIVOT, 'V' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 140 J = 1, N - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 130 I = 1, M
                        TEMP = A( I, J+1 )
                        A( I, J+1 ) = CTEMP*TEMP - STEMP*A( I, J )
                        A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
  130                CONTINUE
                  END IF
  140          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 160 J = N - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 150 I = 1, M
                        TEMP = A( I, J+1 )
                        A( I, J+1 ) = CTEMP*TEMP - STEMP*A( I, J )
                        A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
  150                CONTINUE
                  END IF
  160          CONTINUE
            END IF
         ELSE IF( LSAME( PIVOT, 'T' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 180 J = 2, N
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 170 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = CTEMP*TEMP - STEMP*A( I, 1 )
                        A( I, 1 ) = STEMP*TEMP + CTEMP*A( I, 1 )
  170                CONTINUE
                  END IF
  180          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 200 J = N, 2, -1
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 190 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = CTEMP*TEMP - STEMP*A( I, 1 )
                        A( I, 1 ) = STEMP*TEMP + CTEMP*A( I, 1 )
  190                CONTINUE
                  END IF
  200          CONTINUE
            END IF
         ELSE IF( LSAME( PIVOT, 'B' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 220 J = 1, N - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 210 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = STEMP*A( I, N ) + CTEMP*TEMP
                        A( I, N ) = CTEMP*A( I, N ) - STEMP*TEMP
  210                CONTINUE
                  END IF
  220          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 240 J = N - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 230 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = STEMP*A( I, N ) + CTEMP*TEMP
                        A( I, N ) = CTEMP*A( I, N ) - STEMP*TEMP
  230                CONTINUE
                  END IF
  240          CONTINUE
            END IF
         END IF
      END IF
*
      RETURN
*
*     End of DLASR
*
      END
      SUBROUTINE DLASRT( ID, N, D, INFO )
*
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          ID
      INTEGER            INFO, N
*     ..
*     .. Array Arguments ..
CID      DOUBLE PRECISION   D( * )
      DIMENSION D( * )
*     ..
*
*  Purpose
*  =======
*
*  Sort the numbers in D in increasing order (if ID = 'I') or
*  in decreasing order (if ID = 'D' ).
*
*  Use Quick Sort, reverting to Insertion sort on arrays of
*  size <= 20. Dimension of STACK limits N to about 2**32.
*
*  Arguments
*  =========
*
*  ID      (input) CHARACTER*1
*          = 'I': sort D in increasing order;
*          = 'D': sort D in decreasing order.
*
*  N       (input) INTEGER
*          The length of the array D.
*
*  D       (input/output) DOUBLE PRECISION array, dimension (N)
*          On entry, the array to be sorted.
*          On exit, D has been sorted into increasing order
*          (D(1) <= ... <= D(N) ) or into decreasing order
*          (D(1) >= ... >= D(N) ), depending on ID.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            SELECT
      PARAMETER          ( SELECT = 20 )
*     ..
*     .. Local Scalars ..
      INTEGER            DIR, ENDD, I, J, START, STKPNT
CID      DOUBLE PRECISION   D1, D2, D3, DMNMX, TMP
*     ..
*     .. Local Arrays ..
      INTEGER            STACK( 2, 32 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     ..
*     .. Executable Statements ..
*
*     Test the input paramters.
*
      INFO = 0
      DIR = -1
      IF( LSAME( ID, 'D' ) ) THEN
         DIR = 0
      ELSE IF( LSAME( ID, 'I' ) ) THEN
         DIR = 1
      END IF
      IF( DIR.EQ.-1 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLASRT', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.LE.1 )
     $   RETURN
*
      STKPNT = 1
      STACK( 1, 1 ) = 1
      STACK( 2, 1 ) = N
   10 CONTINUE
      START = STACK( 1, STKPNT )
      ENDD = STACK( 2, STKPNT )
      STKPNT = STKPNT - 1
      IF( ENDD-START.LE.SELECT .AND. ENDD-START.GT.0 ) THEN
*
*        Do Insertion sort on D( START:ENDD )
*
         IF( DIR.EQ.0 ) THEN
*
*           Sort into decreasing order
*
            DO 30 I = START + 1, ENDD
               DO 20 J = I, START + 1, -1
                  IF( D( J ).GT.D( J-1 ) ) THEN
                     DMNMX = D( J )
                     D( J ) = D( J-1 )
                     D( J-1 ) = DMNMX
                  ELSE
                     GO TO 30
                  END IF
   20          CONTINUE
   30       CONTINUE
*
         ELSE
*
*           Sort into increasing order
*
            DO 50 I = START + 1, ENDD
               DO 40 J = I, START + 1, -1
                  IF( D( J ).LT.D( J-1 ) ) THEN
                     DMNMX = D( J )
                     D( J ) = D( J-1 )
                     D( J-1 ) = DMNMX
                  ELSE
                     GO TO 50
                  END IF
   40          CONTINUE
   50       CONTINUE
*
         END IF
*
      ELSE IF( ENDD-START.GT.SELECT ) THEN
*
*        Partition D( START:ENDD ) and stack parts, largest one first
*
*        Choose partition entry as median of 3
*
         D1 = D( START )
         D2 = D( ENDD )
         I = ( START+ENDD ) / 2
         D3 = D( I )
         IF( D1.LT.D2 ) THEN
            IF( D3.LT.D1 ) THEN
               DMNMX = D1
            ELSE IF( D3.LT.D2 ) THEN
               DMNMX = D3
            ELSE
               DMNMX = D2
            END IF
         ELSE
            IF( D3.LT.D2 ) THEN
               DMNMX = D2
            ELSE IF( D3.LT.D1 ) THEN
               DMNMX = D3
            ELSE
               DMNMX = D1
            END IF
         END IF
*
         IF( DIR.EQ.0 ) THEN
*
*           Sort into decreasing order
*
            I = START - 1
            J = ENDD + 1
   60       CONTINUE
   70       CONTINUE
            J = J - 1
            IF( D( J ).LT.DMNMX )
     $         GO TO 70
   80       CONTINUE
            I = I + 1
            IF( D( I ).GT.DMNMX )
     $         GO TO 80
            IF( I.LT.J ) THEN
               TMP = D( I )
               D( I ) = D( J )
               D( J ) = TMP
               GO TO 60
            END IF
            IF( J-START.GT.ENDD-J-1 ) THEN
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
            ELSE
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
            END IF
         ELSE
*
*           Sort into increasing order
*
            I = START - 1
            J = ENDD + 1
   90       CONTINUE
  100       CONTINUE
            J = J - 1
            IF( D( J ).GT.DMNMX )
     $         GO TO 100
  110       CONTINUE
            I = I + 1
            IF( D( I ).LT.DMNMX )
     $         GO TO 110
            IF( I.LT.J ) THEN
               TMP = D( I )
               D( I ) = D( J )
               D( J ) = TMP
               GO TO 90
            END IF
            IF( J-START.GT.ENDD-J-1 ) THEN
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
            ELSE
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
            END IF
         END IF
      END IF
      IF( STKPNT.GT.0 )
     $   GO TO 10
      RETURN
*
*     End of DLASRT
*
      END
      SUBROUTINE DLASSQ( N, X, INCX, SCALE, SUMSQ )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INCX, N
CID      DOUBLE PRECISION   SCALE, SUMSQ
*     ..
*     .. Array Arguments ..
CID      DOUBLE PRECISION   X( * )
      DIMENSION   X( * )
*     ..
*
*  Purpose
*  =======
*
*  DLASSQ  returns the values  scl  and  smsq  such that
*
*     ( scl**2 )*smsq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,
*
*  where  x( i ) = X( 1 + ( i - 1 )*INCX ). The value of  sumsq  is
*  assumed to be non-negative and  scl  returns the value
*
*     scl = max( scale, abs( x( i ) ) ).
*
*  scale and sumsq must be supplied in SCALE and SUMSQ and
*  scl and smsq are overwritten on SCALE and SUMSQ respectively.
*
*  The routine makes only one pass through the vector x.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The number of elements to be used from the vector X.
*
*  X       (input) DOUBLE PRECISION array, dimension (N)
*          The vector for which a scaled sum of squares is computed.
*             x( i )  = X( 1 + ( i - 1 )*INCX ), 1 <= i <= n.
*
*  INCX    (input) INTEGER
*          The increment between successive values of the vector X.
*          INCX > 0.
*
*  SCALE   (input/output) DOUBLE PRECISION
*          On entry, the value  scale  in the equation above.
*          On exit, SCALE is overwritten with  scl , the scaling factor
*          for the sum of squares.
*
*  SUMSQ   (input/output) DOUBLE PRECISION
*          On entry, the value  sumsq  in the equation above.
*          On exit, SUMSQ is overwritten with  smsq , the basic sum of
*          squares from which  scl  has been factored out.
*
* =====================================================================
*
*     .. Parameters ..
CID      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            IX
CID      DOUBLE PRECISION   ABSXI
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS
*     ..
*     .. Executable Statements ..
*
      IF( N.GT.0 ) THEN
         DO 10 IX = 1, 1 + ( N-1 )*INCX, INCX
            IF( X( IX ).NE.ZERO ) THEN
               ABSXI = ABS( X( IX ) )
               IF( SCALE.LT.ABSXI ) THEN
                  SUMSQ = 1 + SUMSQ*( SCALE / ABSXI )**2
                  SCALE = ABSXI
               ELSE
                  SUMSQ = SUMSQ + ( ABSXI / SCALE )**2
               END IF
            END IF
   10    CONTINUE
      END IF
      RETURN
*
*     End of DLASSQ
*
      END
      SUBROUTINE DLASV2( F, G, H, SSMIN, SSMAX, SNR, CSR, SNL, CSL )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
CID      DOUBLE PRECISION   CSL, CSR, F, G, H, SNL, SNR, SSMAX, SSMIN
*     ..
*
*  Purpose
*  =======
*
*  DLASV2 computes the singular value decomposition of a 2-by-2
*  triangular matrix
*     [  F   G  ]
*     [  0   H  ].
*  On return, abs(SSMAX) is the larger singular value, abs(SSMIN) is the
*  smaller singular value, and (CSL,SNL) and (CSR,SNR) are the left and
*  right singular vectors for abs(SSMAX), giving the decomposition
*
*     [ CSL  SNL ] [  F   G  ] [ CSR -SNR ]  =  [ SSMAX   0   ]
*     [-SNL  CSL ] [  0   H  ] [ SNR  CSR ]     [  0    SSMIN ].
*
*  Arguments
*  =========
*
*  F       (input) DOUBLE PRECISION
*          The (1,1) element of the 2-by-2 matrix.
*
*  G       (input) DOUBLE PRECISION
*          The (1,2) element of the 2-by-2 matrix.
*
*  H       (input) DOUBLE PRECISION
*          The (2,2) element of the 2-by-2 matrix.
*
*  SSMIN   (output) DOUBLE PRECISION
*          abs(SSMIN) is the smaller singular value.
*
*  SSMAX   (output) DOUBLE PRECISION
*          abs(SSMAX) is the larger singular value.
*
*  SNL     (output) DOUBLE PRECISION
*  CSL     (output) DOUBLE PRECISION
*          The vector (CSL, SNL) is a unit left singular vector for the
*          singular value abs(SSMAX).
*
*  SNR     (output) DOUBLE PRECISION
*  CSR     (output) DOUBLE PRECISION
*          The vector (CSR, SNR) is a unit right singular vector for the
*          singular value abs(SSMAX).
*
*  Further Details
*  ===============
*
*  Any input parameter may be aliased with any output parameter.
*
*  Barring over/underflow and assuming a guard digit in subtraction, all
*  output quantities are correct to within a few units in the last
*  place (ulps).
*
*  In IEEE arithmetic, the code works correctly if one matrix element is
*  infinite.
*
*  Overflow will not occur unless the largest singular value itself
*  overflows or is within a few ulps of overflow. (On machines with
*  partial overflow, like the Cray, overflow may occur if the largest
*  singular value is within a factor of 2 of overflow.)
*
*  Underflow is harmless if underflow is gradual. Otherwise, results
*  may correspond to a matrix modified by perturbations of size near
*  the underflow threshold.
*
* =====================================================================
*
*     .. Parameters ..
CID      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
CID      DOUBLE PRECISION   HALF
      PARAMETER          ( HALF = 0.5D0 )
CID      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
CID      DOUBLE PRECISION   TWO
      PARAMETER          ( TWO = 2.0D0 )
CID      DOUBLE PRECISION   FOUR
      PARAMETER          ( FOUR = 4.0D0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            GASMAL, SWAP
      INTEGER            PMAX
      REAL               L,M,MM
CID      DOUBLE PRECISION   A, CLT, CRT, D, FA, FT, GA, GT, HA, HT, L, M,
CID     $                   MM, R, S, SLT, SRT, T, TEMP, TSIGN, TT
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, SIGN, SQRT
*     ..
*     .. External Functions ..
CID      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
*     ..
*     .. Executable Statements ..
*
      FT = F
      FA = ABS( FT )
      HT = H
      HA = ABS( H )
*
*     PMAX points to the maximum absolute element of matrix
*       PMAX = 1 if F largest in absolute values
*       PMAX = 2 if G largest in absolute values
*       PMAX = 3 if H largest in absolute values
*
      PMAX = 1
      SWAP = ( HA.GT.FA )
      IF( SWAP ) THEN
         PMAX = 3
         TEMP = FT
         FT = HT
         HT = TEMP
         TEMP = FA
         FA = HA
         HA = TEMP
*
*        Now FA .ge. HA
*
      END IF
      GT = G
      GA = ABS( GT )
      IF( GA.EQ.ZERO ) THEN
*
*        Diagonal matrix
*
         SSMIN = HA
         SSMAX = FA
         CLT = ONE
         CRT = ONE
         SLT = ZERO
         SRT = ZERO
      ELSE
         GASMAL = .TRUE.
         IF( GA.GT.FA ) THEN
            PMAX = 2
            IF( ( FA / GA ).LT.DLAMCH( 'EPS' ) ) THEN
*
*              Case of very large GA
*
               GASMAL = .FALSE.
               SSMAX = GA
               IF( HA.GT.ONE ) THEN
                  SSMIN = FA / ( GA / HA )
               ELSE
                  SSMIN = ( FA / GA )*HA
               END IF
               CLT = ONE
               SLT = HT / GT
               SRT = ONE
               CRT = FT / GT
            END IF
         END IF
         IF( GASMAL ) THEN
*
*           Normal case
*
            D = FA - HA
            IF( D.EQ.FA ) THEN
*
*              Copes with infinite F or H
*
               L = ONE
            ELSE
               L = D / FA
            END IF
*
*           Note that 0 .le. L .le. 1
*
            M = GT / FT
*
*           Note that abs(M) .le. 1/macheps
*
            T = TWO - L
*
*           Note that T .ge. 1
*
            MM = M*M
            TT = T*T
            S = SQRT( TT+MM )
*
*           Note that 1 .le. S .le. 1 + 1/macheps
*
            IF( L.EQ.ZERO ) THEN
               R = ABS( M )
            ELSE
               R = SQRT( L*L+MM )
            END IF
*
*           Note that 0 .le. R .le. 1 + 1/macheps
*
            A = HALF*( S+R )
*
*           Note that 1 .le. A .le. 1 + abs(M)
*
            SSMIN = HA / A
            SSMAX = FA*A
            IF( MM.EQ.ZERO ) THEN
*
*              Note that M is very tiny
*
               IF( L.EQ.ZERO ) THEN
                  T = SIGN( TWO, FT )*SIGN( ONE, GT )
               ELSE
                  T = GT / SIGN( D, FT ) + M / T
               END IF
            ELSE
               T = ( M / ( S+T )+M / ( R+L ) )*( ONE+A )
            END IF
            L = SQRT( T*T+FOUR )
            CRT = TWO / L
            SRT = T / L
            CLT = ( CRT+SRT*M ) / A
            SLT = ( HT / FT )*SRT / A
         END IF
      END IF
      IF( SWAP ) THEN
         CSL = SRT
         SNL = CRT
         CSR = SLT
         SNR = CLT
      ELSE
         CSL = CLT
         SNL = SLT
         CSR = CRT
         SNR = SRT
      END IF
*
*     Correct signs of SSMAX and SSMIN
*
      IF( PMAX.EQ.1 )
     $   TSIGN = SIGN( ONE, CSR )*SIGN( ONE, CSL )*SIGN( ONE, F )
      IF( PMAX.EQ.2 )
     $   TSIGN = SIGN( ONE, SNR )*SIGN( ONE, CSL )*SIGN( ONE, G )
      IF( PMAX.EQ.3 )
     $   TSIGN = SIGN( ONE, SNR )*SIGN( ONE, SNL )*SIGN( ONE, H )
      SSMAX = SIGN( SSMAX, TSIGN )
      SSMIN = SIGN( SSMIN, TSIGN*SIGN( ONE, F )*SIGN( ONE, H ) )
      RETURN
*
*     End of DLASV2
*
      END
      SUBROUTINE DLAZQ3( I0, N0, Z, PP, DMIN, SIGMA, DESIG, QMAX, NFAIL,
     $                   ITER, NDIV, IEEE, TTYPE, DMIN1, DMIN2, DN, DN1,
     $                   DN2, TAU )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      LOGICAL            IEEE
      INTEGER            I0, ITER, N0, NDIV, NFAIL, PP, TTYPE
CID      DOUBLE PRECISION   DESIG, DMIN, DMIN1, DMIN2, DN, DN1, DN2, QMAX,
CID     $                   SIGMA, TAU
*     ..
*     .. Array Arguments ..
CID      DOUBLE PRECISION   Z( * )
      DIMENSION  Z( * )
*     ..
*
*  Purpose
*  =======
*
*  DLAZQ3 checks for deflation, computes a shift (TAU) and calls dqds.
*  In case of failure it changes shifts, and tries again until output
*  is positive.
*
*  Arguments
*  =========
*
*  I0     (input) INTEGER
*         First index.
*
*  N0     (input) INTEGER
*         Last index.
*
*  Z      (input) DOUBLE PRECISION array, dimension ( 4*N )
*         Z holds the qd array.
*
*  PP     (input) INTEGER
*         PP=0 for ping, PP=1 for pong.
*
*  DMIN   (output) DOUBLE PRECISION
*         Minimum value of d.
*
*  SIGMA  (output) DOUBLE PRECISION
*         Sum of shifts used in current segment.
*
*  DESIG  (input/output) DOUBLE PRECISION
*         Lower order part of SIGMA
*
*  QMAX   (input) DOUBLE PRECISION
*         Maximum value of q.
*
*  NFAIL  (output) INTEGER
*         Number of times shift was too big.
*
*  ITER   (output) INTEGER
*         Number of iterations.
*
*  NDIV   (output) INTEGER
*         Number of divisions.
*
*  IEEE   (input) LOGICAL
*         Flag for IEEE or non IEEE arithmetic (passed to DLASQ5).
*
*  TTYPE  (input/output) INTEGER
*         Shift type.  TTYPE is passed as an argument in order to save
*         its value between calls to DLAZQ3
*
*  DMIN1  (input/output) REAL
*  DMIN2  (input/output) REAL
*  DN     (input/output) REAL
*  DN1    (input/output) REAL
*  DN2    (input/output) REAL
*  TAU    (input/output) REAL
*         These are passed as arguments in order to save their values
*         between calls to DLAZQ3
*
*  This is a thread safe version of DLASQ3, which passes TTYPE, DMIN1,
*  DMIN2, DN, DN1. DN2 and TAU through the argument list in place of
*  declaring them in a SAVE statment.
*
*  =====================================================================
*
*     .. Parameters ..
CID      DOUBLE PRECISION   CBIAS
      PARAMETER          ( CBIAS = 1.50D0 )
CID      DOUBLE PRECISION   ZERO, QURTR, HALF, ONE, TWO, HUNDRD
      PARAMETER          ( ZERO = 0.0D0, QURTR = 0.250D0, HALF = 0.5D0,
     $                     ONE = 1.0D0, TWO = 2.0D0, HUNDRD = 100.0D0 )
*     ..
*     .. Local Scalars ..
      INTEGER            IPN4, J4, N0IN, NN
CID      DOUBLE PRECISION   EPS, G, S, SAFMIN, T, TEMP, TOL, TOL2
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLASQ5, DLASQ6, DLAZQ4
*     ..
*     .. External Function ..
CID      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MIN, SQRT
*     ..
*     .. Executable Statements ..
*
      N0IN   = N0
      EPS    = DLAMCH( 'Precision' )
      SAFMIN = DLAMCH( 'Safe minimum' )
      TOL    = EPS*HUNDRD
      TOL2   = TOL**2
      G      = ZERO
*
*     Check for deflation.
*
   10 CONTINUE
*
      IF( N0.LT.I0 )
     $   RETURN
      IF( N0.EQ.I0 )
     $   GO TO 20
      NN = 4*N0 + PP
      IF( N0.EQ.( I0+1 ) )
     $   GO TO 40
*
*     Check whether E(N0-1) is negligible, 1 eigenvalue.
*
      IF( Z( NN-5 ).GT.TOL2*( SIGMA+Z( NN-3 ) ) .AND.
     $    Z( NN-2*PP-4 ).GT.TOL2*Z( NN-7 ) )
     $   GO TO 30
*
   20 CONTINUE
*
      Z( 4*N0-3 ) = Z( 4*N0+PP-3 ) + SIGMA
      N0 = N0 - 1
      GO TO 10
*
*     Check  whether E(N0-2) is negligible, 2 eigenvalues.
*
   30 CONTINUE
*
      IF( Z( NN-9 ).GT.TOL2*SIGMA .AND.
     $    Z( NN-2*PP-8 ).GT.TOL2*Z( NN-11 ) )
     $   GO TO 50
*
   40 CONTINUE
*
      IF( Z( NN-3 ).GT.Z( NN-7 ) ) THEN
         S = Z( NN-3 )
         Z( NN-3 ) = Z( NN-7 )
         Z( NN-7 ) = S
      END IF
      IF( Z( NN-5 ).GT.Z( NN-3 )*TOL2 ) THEN
         T = HALF*( ( Z( NN-7 )-Z( NN-3 ) )+Z( NN-5 ) )
         S = Z( NN-3 )*( Z( NN-5 ) / T )
         IF( S.LE.T ) THEN
            S = Z( NN-3 )*( Z( NN-5 ) /
     $          ( T*( ONE+SQRT( ONE+S / T ) ) ) )
         ELSE
            S = Z( NN-3 )*( Z( NN-5 ) / ( T+SQRT( T )*SQRT( T+S ) ) )
         END IF
         T = Z( NN-7 ) + ( S+Z( NN-5 ) )
         Z( NN-3 ) = Z( NN-3 )*( Z( NN-7 ) / T )
         Z( NN-7 ) = T
      END IF
      Z( 4*N0-7 ) = Z( NN-7 ) + SIGMA
      Z( 4*N0-3 ) = Z( NN-3 ) + SIGMA
      N0 = N0 - 2
      GO TO 10
*
   50 CONTINUE
*
*     Reverse the qd-array, if warranted.
*
      IF( DMIN.LE.ZERO .OR. N0.LT.N0IN ) THEN
         IF( CBIAS*Z( 4*I0+PP-3 ).LT.Z( 4*N0+PP-3 ) ) THEN
            IPN4 = 4*( I0+N0 )
            DO 60 J4 = 4*I0, 2*( I0+N0-1 ), 4
               TEMP = Z( J4-3 )
               Z( J4-3 ) = Z( IPN4-J4-3 )
               Z( IPN4-J4-3 ) = TEMP
               TEMP = Z( J4-2 )
               Z( J4-2 ) = Z( IPN4-J4-2 )
               Z( IPN4-J4-2 ) = TEMP
               TEMP = Z( J4-1 )
               Z( J4-1 ) = Z( IPN4-J4-5 )
               Z( IPN4-J4-5 ) = TEMP
               TEMP = Z( J4 )
               Z( J4 ) = Z( IPN4-J4-4 )
               Z( IPN4-J4-4 ) = TEMP
   60       CONTINUE
            IF( N0-I0.LE.4 ) THEN
               Z( 4*N0+PP-1 ) = Z( 4*I0+PP-1 )
               Z( 4*N0-PP ) = Z( 4*I0-PP )
            END IF
            DMIN2 = MIN( DMIN2, Z( 4*N0+PP-1 ) )
            Z( 4*N0+PP-1 ) = MIN( Z( 4*N0+PP-1 ), Z( 4*I0+PP-1 ),
     $                            Z( 4*I0+PP+3 ) )
            Z( 4*N0-PP ) = MIN( Z( 4*N0-PP ), Z( 4*I0-PP ),
     $                          Z( 4*I0-PP+4 ) )
            QMAX = MAX( QMAX, Z( 4*I0+PP-3 ), Z( 4*I0+PP+1 ) )
            DMIN = -ZERO
         END IF
      END IF
*
      IF( DMIN.LT.ZERO .OR. SAFMIN*QMAX.LT.MIN( Z( 4*N0+PP-1 ),
     $    Z( 4*N0+PP-9 ), DMIN2+Z( 4*N0-PP ) ) ) THEN
*
*        Choose a shift.
*
         CALL DLAZQ4( I0, N0, Z, PP, N0IN, DMIN, DMIN1, DMIN2, DN, DN1,
     $                DN2, TAU, TTYPE, G )
*
*        Call dqds until DMIN > 0.
*
   80    CONTINUE
*
         CALL DLASQ5( I0, N0, Z, PP, TAU, DMIN, DMIN1, DMIN2, DN,
     $                DN1, DN2, IEEE )
*
         NDIV = NDIV + ( N0-I0+2 )
         ITER = ITER + 1
*
*        Check status.
*
         IF( DMIN.GE.ZERO .AND. DMIN1.GT.ZERO ) THEN
*
*           Success.
*
            GO TO 100
*
         ELSE IF( DMIN.LT.ZERO .AND. DMIN1.GT.ZERO .AND.
     $            Z( 4*( N0-1 )-PP ).LT.TOL*( SIGMA+DN1 ) .AND.
     $            ABS( DN ).LT.TOL*SIGMA ) THEN
*
*           Convergence hidden by negative DN.
*
            Z( 4*( N0-1 )-PP+2 ) = ZERO
            DMIN = ZERO
            GO TO 100
         ELSE IF( DMIN.LT.ZERO ) THEN
*
*           TAU too big. Select new TAU and try again.
*
            NFAIL = NFAIL + 1
            IF( TTYPE.LT.-22 ) THEN
*
*              Failed twice. Play it safe.
*
               TAU = ZERO
            ELSE IF( DMIN1.GT.ZERO ) THEN
*
*              Late failure. Gives excellent shift.
*
               TAU = ( TAU+DMIN )*( ONE-TWO*EPS )
               TTYPE = TTYPE - 11
            ELSE
*
*              Early failure. Divide by 4.
*
               TAU = QURTR*TAU
               TTYPE = TTYPE - 12
            END IF
            GO TO 80
         ELSE IF( DMIN.NE.DMIN ) THEN
*
*           NaN.
*
            TAU = ZERO
            GO TO 80
         ELSE
*
*           Possible underflow. Play it safe.
*
            GO TO 90
         END IF
      END IF
*
*     Risk of underflow.
*
   90 CONTINUE
      CALL DLASQ6( I0, N0, Z, PP, DMIN, DMIN1, DMIN2, DN, DN1, DN2 )
      NDIV = NDIV + ( N0-I0+2 )
      ITER = ITER + 1
      TAU = ZERO
*
  100 CONTINUE
      IF( TAU.LT.SIGMA ) THEN
         DESIG = DESIG + TAU
         T = SIGMA + DESIG
         DESIG = DESIG - ( T-SIGMA )
      ELSE
         T = SIGMA + TAU
         DESIG = SIGMA - ( T-TAU ) + DESIG
      END IF
      SIGMA = T
*
      RETURN
*
*     End of DLAZQ3
*
      END
      SUBROUTINE DLAZQ4( I0, N0, Z, PP, N0IN, DMIN, DMIN1, DMIN2, DN,
     $                   DN1, DN2, TAU, TTYPE, G )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            I0, N0, N0IN, PP, TTYPE
CID      DOUBLE PRECISION   DMIN, DMIN1, DMIN2, DN, DN1, DN2, G, TAU
*     ..
*     .. Array Arguments ..
CID      DOUBLE PRECISION   Z( * )
      DIMENSION Z ( * )
*     ..
*
*  Purpose
*  =======
*
*  DLAZQ4 computes an approximation TAU to the smallest eigenvalue 
*  using values of d from the previous transform.
*
*  I0    (input) INTEGER
*        First index.
*
*  N0    (input) INTEGER
*        Last index.
*
*  Z     (input) DOUBLE PRECISION array, dimension ( 4*N )
*        Z holds the qd array.
*
*  PP    (input) INTEGER
*        PP=0 for ping, PP=1 for pong.
*
*  N0IN  (input) INTEGER
*        The value of N0 at start of EIGTEST.
*
*  DMIN  (input) DOUBLE PRECISION
*        Minimum value of d.
*
*  DMIN1 (input) DOUBLE PRECISION
*        Minimum value of d, excluding D( N0 ).
*
*  DMIN2 (input) DOUBLE PRECISION
*        Minimum value of d, excluding D( N0 ) and D( N0-1 ).
*
*  DN    (input) DOUBLE PRECISION
*        d(N)
*
*  DN1   (input) DOUBLE PRECISION
*        d(N-1)
*
*  DN2   (input) DOUBLE PRECISION
*        d(N-2)
*
*  TAU   (output) DOUBLE PRECISION
*        This is the shift.
*
*  TTYPE (output) INTEGER
*        Shift type.
*
*  G     (input/output) DOUBLE PRECISION
*        G is passed as an argument in order to save its value between
*        calls to DLAZQ4
*
*  Further Details
*  ===============
*  CNST1 = 9/16
*
*  This is a thread safe version of DLASQ4, which passes G through the
*  argument list in place of declaring G in a SAVE statment.
*
*  =====================================================================
*
*     .. Parameters ..
CID      DOUBLE PRECISION   CNST1, CNST2, CNST3
      PARAMETER          ( CNST1 = 0.5630D0, CNST2 = 1.010D0,
     $                   CNST3 = 1.050D0 )
CID      DOUBLE PRECISION   QURTR, THIRD, HALF, ZERO, ONE, TWO, HUNDRD
      PARAMETER          ( QURTR = 0.250D0, THIRD = 0.3330D0,
     $                   HALF = 0.50D0, ZERO = 0.0D0, ONE = 1.0D0,
     $                   TWO = 2.0D0, HUNDRD = 100.0D0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I4, NN, NP
CID      DOUBLE PRECISION   A2, B1, B2, GAM, GAP1, GAP2, S
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
*
*     A negative DMIN forces the shift to take that absolute value
*     TTYPE records the type of shift.
*
      IF( DMIN.LE.ZERO ) THEN
         TAU = -DMIN
         TTYPE = -1
         RETURN
      END IF
*       
      NN = 4*N0 + PP
      IF( N0IN.EQ.N0 ) THEN
*
*        No eigenvalues deflated.
*
         IF( DMIN.EQ.DN .OR. DMIN.EQ.DN1 ) THEN
*
            B1 = SQRT( Z( NN-3 ) )*SQRT( Z( NN-5 ) )
            B2 = SQRT( Z( NN-7 ) )*SQRT( Z( NN-9 ) )
            A2 = Z( NN-7 ) + Z( NN-5 )
*
*           Cases 2 and 3.
*
            IF( DMIN.EQ.DN .AND. DMIN1.EQ.DN1 ) THEN
               GAP2 = DMIN2 - A2 - DMIN2*QURTR
               IF( GAP2.GT.ZERO .AND. GAP2.GT.B2 ) THEN
                  GAP1 = A2 - DN - ( B2 / GAP2 )*B2
               ELSE
                  GAP1 = A2 - DN - ( B1+B2 )
               END IF
               IF( GAP1.GT.ZERO .AND. GAP1.GT.B1 ) THEN
                  S = MAX( DN-( B1 / GAP1 )*B1, HALF*DMIN )
                  TTYPE = -2
               ELSE
                  S = ZERO
                  IF( DN.GT.B1 )
     $               S = DN - B1
                  IF( A2.GT.( B1+B2 ) )
     $               S = MIN( S, A2-( B1+B2 ) )
                  S = MAX( S, THIRD*DMIN )
                  TTYPE = -3
               END IF
            ELSE
*
*              Case 4.
*
               TTYPE = -4
               S = QURTR*DMIN
               IF( DMIN.EQ.DN ) THEN
                  GAM = DN
                  A2 = ZERO
                  IF( Z( NN-5 ) .GT. Z( NN-7 ) )
     $               RETURN
                  B2 = Z( NN-5 ) / Z( NN-7 )
                  NP = NN - 9
               ELSE
                  NP = NN - 2*PP
                  B2 = Z( NP-2 )
                  GAM = DN1
                  IF( Z( NP-4 ) .GT. Z( NP-2 ) )
     $               RETURN
                  A2 = Z( NP-4 ) / Z( NP-2 )
                  IF( Z( NN-9 ) .GT. Z( NN-11 ) )
     $               RETURN
                  B2 = Z( NN-9 ) / Z( NN-11 )
                  NP = NN - 13
               END IF
*
*              Approximate contribution to norm squared from I < NN-1.
*
               A2 = A2 + B2
               DO 10 I4 = NP, 4*I0 - 1 + PP, -4
                  IF( B2.EQ.ZERO )
     $               GO TO 20
                  B1 = B2
                  IF( Z( I4 ) .GT. Z( I4-2 ) )
     $               RETURN
                  B2 = B2*( Z( I4 ) / Z( I4-2 ) )
                  A2 = A2 + B2
                  IF( HUNDRD*MAX( B2, B1 ).LT.A2 .OR. CNST1.LT.A2 ) 
     $               GO TO 20
   10          CONTINUE
   20          CONTINUE
               A2 = CNST3*A2
*
*              Rayleigh quotient residual bound.
*
               IF( A2.LT.CNST1 )
     $            S = GAM*( ONE-SQRT( A2 ) ) / ( ONE+A2 )
            END IF
         ELSE IF( DMIN.EQ.DN2 ) THEN
*
*           Case 5.
*
            TTYPE = -5
            S = QURTR*DMIN
*
*           Compute contribution to norm squared from I > NN-2.
*
            NP = NN - 2*PP
            B1 = Z( NP-2 )
            B2 = Z( NP-6 )
            GAM = DN2
            IF( Z( NP-8 ).GT.B2 .OR. Z( NP-4 ).GT.B1 )
     $         RETURN
            A2 = ( Z( NP-8 ) / B2 )*( ONE+Z( NP-4 ) / B1 )
*
*           Approximate contribution to norm squared from I < NN-2.
*
            IF( N0-I0.GT.2 ) THEN
               B2 = Z( NN-13 ) / Z( NN-15 )
               A2 = A2 + B2
               DO 30 I4 = NN - 17, 4*I0 - 1 + PP, -4
                  IF( B2.EQ.ZERO )
     $               GO TO 40
                  B1 = B2
                  IF( Z( I4 ) .GT. Z( I4-2 ) )
     $               RETURN
                  B2 = B2*( Z( I4 ) / Z( I4-2 ) )
                  A2 = A2 + B2
                  IF( HUNDRD*MAX( B2, B1 ).LT.A2 .OR. CNST1.LT.A2 ) 
     $               GO TO 40
   30          CONTINUE
   40          CONTINUE
               A2 = CNST3*A2
            END IF
*
            IF( A2.LT.CNST1 )
     $         S = GAM*( ONE-SQRT( A2 ) ) / ( ONE+A2 )
         ELSE
*
*           Case 6, no information to guide us.
*
            IF( TTYPE.EQ.-6 ) THEN
               G = G + THIRD*( ONE-G )
            ELSE IF( TTYPE.EQ.-18 ) THEN
               G = QURTR*THIRD
            ELSE
               G = QURTR
            END IF
            S = G*DMIN
            TTYPE = -6
         END IF
*
      ELSE IF( N0IN.EQ.( N0+1 ) ) THEN
*
*        One eigenvalue just deflated. Use DMIN1, DN1 for DMIN and DN.
*
         IF( DMIN1.EQ.DN1 .AND. DMIN2.EQ.DN2 ) THEN 
*
*           Cases 7 and 8.
*
            TTYPE = -7
            S = THIRD*DMIN1
            IF( Z( NN-5 ).GT.Z( NN-7 ) )
     $         RETURN
            B1 = Z( NN-5 ) / Z( NN-7 )
            B2 = B1
            IF( B2.EQ.ZERO )
     $         GO TO 60
            DO 50 I4 = 4*N0 - 9 + PP, 4*I0 - 1 + PP, -4
               A2 = B1
               IF( Z( I4 ).GT.Z( I4-2 ) )
     $            RETURN
               B1 = B1*( Z( I4 ) / Z( I4-2 ) )
               B2 = B2 + B1
               IF( HUNDRD*MAX( B1, A2 ).LT.B2 ) 
     $            GO TO 60
   50       CONTINUE
   60       CONTINUE
            B2 = SQRT( CNST3*B2 )
            A2 = DMIN1 / ( ONE+B2**2 )
            GAP2 = HALF*DMIN2 - A2
            IF( GAP2.GT.ZERO .AND. GAP2.GT.B2*A2 ) THEN
               S = MAX( S, A2*( ONE-CNST2*A2*( B2 / GAP2 )*B2 ) )
            ELSE 
               S = MAX( S, A2*( ONE-CNST2*B2 ) )
               TTYPE = -8
            END IF
         ELSE
*
*           Case 9.
*
            S = QURTR*DMIN1
            IF( DMIN1.EQ.DN1 )
     $         S = HALF*DMIN1
            TTYPE = -9
         END IF
*
      ELSE IF( N0IN.EQ.( N0+2 ) ) THEN
*
*        Two eigenvalues deflated. Use DMIN2, DN2 for DMIN and DN.
*
*        Cases 10 and 11.
*
         IF( DMIN2.EQ.DN2 .AND. TWO*Z( NN-5 ).LT.Z( NN-7 ) ) THEN 
            TTYPE = -10
            S = THIRD*DMIN2
            IF( Z( NN-5 ).GT.Z( NN-7 ) )
     $         RETURN
            B1 = Z( NN-5 ) / Z( NN-7 )
            B2 = B1
            IF( B2.EQ.ZERO )
     $         GO TO 80
            DO 70 I4 = 4*N0 - 9 + PP, 4*I0 - 1 + PP, -4
               IF( Z( I4 ).GT.Z( I4-2 ) )
     $            RETURN
               B1 = B1*( Z( I4 ) / Z( I4-2 ) )
               B2 = B2 + B1
               IF( HUNDRD*B1.LT.B2 )
     $            GO TO 80
   70       CONTINUE
   80       CONTINUE
            B2 = SQRT( CNST3*B2 )
            A2 = DMIN2 / ( ONE+B2**2 )
            GAP2 = Z( NN-7 ) + Z( NN-9 ) -
     $             SQRT( Z( NN-11 ) )*SQRT( Z( NN-9 ) ) - A2
            IF( GAP2.GT.ZERO .AND. GAP2.GT.B2*A2 ) THEN
               S = MAX( S, A2*( ONE-CNST2*A2*( B2 / GAP2 )*B2 ) )
            ELSE 
               S = MAX( S, A2*( ONE-CNST2*B2 ) )
            END IF
         ELSE
            S = QURTR*DMIN2
            TTYPE = -11
         END IF
      ELSE IF( N0IN.GT.( N0+2 ) ) THEN
*
*        Case 12, more than two eigenvalues deflated. No information.
*
         S = ZERO 
         TTYPE = -12
      END IF
*
      TAU = S
      RETURN
*
*     End of DLAZQ4
*
      END
      SUBROUTINE DORG2R( M, N, K, A, LDA, TAU, WORK, INFO )
*
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, M, N
*     ..
*     .. Array Arguments ..
CID      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
      DIMENSION A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DORG2R generates an m by n real matrix Q with orthonormal columns,
*  which is defined as the first n columns of a product of k elementary
*  reflectors of order m
*
*        Q  =  H(1) H(2) . . . H(k)
*
*  as returned by DGEQRF.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix Q. M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix Q. M >= N >= 0.
*
*  K       (input) INTEGER
*          The number of elementary reflectors whose product defines the
*          matrix Q. N >= K >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the i-th column must contain the vector which
*          defines the elementary reflector H(i), for i = 1,2,...,k, as
*          returned by DGEQRF in the first k columns of its array
*          argument A.
*          On exit, the m-by-n matrix Q.
*
*  LDA     (input) INTEGER
*          The first dimension of the array A. LDA >= max(1,M).
*
*  TAU     (input) DOUBLE PRECISION array, dimension (K)
*          TAU(i) must contain the scalar factor of the elementary
*          reflector H(i), as returned by DGEQRF.
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (N)
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument has an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
CID      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J, L
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLARF, DSCAL, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORG2R', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.LE.0 )
     $   RETURN
*
*     Initialise columns k+1:n to columns of the unit matrix
*
      DO 20 J = K + 1, N
         DO 10 L = 1, M
            A( L, J ) = ZERO
   10    CONTINUE
         A( J, J ) = ONE
   20 CONTINUE
*
      DO 40 I = K, 1, -1
*
*        Apply H(i) to A(i:m,i:n) from the left
*
         IF( I.LT.N ) THEN
            A( I, I ) = ONE
            CALL DLARF( 'Left', M-I+1, N-I, A( I, I ), 1, TAU( I ),
     $                  A( I, I+1 ), LDA, WORK )
         END IF
         IF( I.LT.M )
     $      CALL DSCAL( M-I, -TAU( I ), A( I+1, I ), 1 )
         A( I, I ) = ONE - TAU( I )
*
*        Set A(1:i-1,i) to zero
*
         DO 30 L = 1, I - 1
            A( L, I ) = ZERO
   30    CONTINUE
   40 CONTINUE
      RETURN
*
*     End of DORG2R
*
      END
      SUBROUTINE DORGBR( VECT, M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
*
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          VECT
      INTEGER            INFO, K, LDA, LWORK, M, N
*     ..
*     .. Array Arguments ..
CID      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
      DIMENSION   A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DORGBR generates one of the real orthogonal matrices Q or P**T
*  determined by DGEBRD when reducing a real matrix A to bidiagonal
*  form: A = Q * B * P**T.  Q and P**T are defined as products of
*  elementary reflectors H(i) or G(i) respectively.
*
*  If VECT = 'Q', A is assumed to have been an M-by-K matrix, and Q
*  is of order M:
*  if m >= k, Q = H(1) H(2) . . . H(k) and DORGBR returns the first n
*  columns of Q, where m >= n >= k;
*  if m < k, Q = H(1) H(2) . . . H(m-1) and DORGBR returns Q as an
*  M-by-M matrix.
*
*  If VECT = 'P', A is assumed to have been a K-by-N matrix, and P**T
*  is of order N:
*  if k < n, P**T = G(k) . . . G(2) G(1) and DORGBR returns the first m
*  rows of P**T, where n >= m >= k;
*  if k >= n, P**T = G(n-1) . . . G(2) G(1) and DORGBR returns P**T as
*  an N-by-N matrix.
*
*  Arguments
*  =========
*
*  VECT    (input) CHARACTER*1
*          Specifies whether the matrix Q or the matrix P**T is
*          required, as defined in the transformation applied by DGEBRD:
*          = 'Q':  generate Q;
*          = 'P':  generate P**T.
*
*  M       (input) INTEGER
*          The number of rows of the matrix Q or P**T to be returned.
*          M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix Q or P**T to be returned.
*          N >= 0.
*          If VECT = 'Q', M >= N >= min(M,K);
*          if VECT = 'P', N >= M >= min(N,K).
*
*  K       (input) INTEGER
*          If VECT = 'Q', the number of columns in the original M-by-K
*          matrix reduced by DGEBRD.
*          If VECT = 'P', the number of rows in the original K-by-N
*          matrix reduced by DGEBRD.
*          K >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the vectors which define the elementary reflectors,
*          as returned by DGEBRD.
*          On exit, the M-by-N matrix Q or P**T.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A. LDA >= max(1,M).
*
*  TAU     (input) DOUBLE PRECISION array, dimension
*                                (min(M,K)) if VECT = 'Q'
*                                (min(N,K)) if VECT = 'P'
*          TAU(i) must contain the scalar factor of the elementary
*          reflector H(i) or G(i), which determines Q or P**T, as
*          returned by DGEBRD in its array argument TAUQ or TAUP.
*
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK. LWORK >= max(1,min(M,N)).
*          For optimum performance LWORK >= min(M,N)*NB, where NB
*          is the optimal blocksize.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
CID      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY, WANTQ
      INTEGER            I, IINFO, J, LWKOPT, MN, NB
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL           DORGLQ, DORGQR, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      WANTQ = LSAME( VECT, 'Q' )
      MN = MIN( M, N )
      LQUERY = ( LWORK.EQ.-1 )
      IF( .NOT.WANTQ .AND. .NOT.LSAME( VECT, 'P' ) ) THEN
         INFO = -1
      ELSE IF( M.LT.0 ) THEN
         INFO = -2
      ELSE IF( N.LT.0 .OR. ( WANTQ .AND. ( N.GT.M .OR. N.LT.MIN( M,
     $         K ) ) ) .OR. ( .NOT.WANTQ .AND. ( M.GT.N .OR. M.LT.
     $         MIN( N, K ) ) ) ) THEN
         INFO = -3
      ELSE IF( K.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -6
      ELSE IF( LWORK.LT.MAX( 1, MN ) .AND. .NOT.LQUERY ) THEN
         INFO = -9
      END IF
*
      IF( INFO.EQ.0 ) THEN
         IF( WANTQ ) THEN
            NB = ILAENV( 1, 'DORGQR', ' ', M, N, K, -1 )
         ELSE
            NB = ILAENV( 1, 'DORGLQ', ' ', M, N, K, -1 )
         END IF
         LWKOPT = MAX( 1, MN )*NB
         WORK( 1 ) = LWKOPT
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORGBR', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
      IF( WANTQ ) THEN
*
*        Form Q, determined by a call to DGEBRD to reduce an m-by-k
*        matrix
*
         IF( M.GE.K ) THEN
*
*           If m >= k, assume m >= n >= k
*
            CALL DORGQR( M, N, K, A, LDA, TAU, WORK, LWORK, IINFO )
*
         ELSE
*
*           If m < k, assume m = n
*
*           Shift the vectors which define the elementary reflectors one
*           column to the right, and set the first row and column of Q
*           to those of the unit matrix
*
            DO 20 J = M, 2, -1
               A( 1, J ) = ZERO
               DO 10 I = J + 1, M
                  A( I, J ) = A( I, J-1 )
   10          CONTINUE
   20       CONTINUE
            A( 1, 1 ) = ONE
            DO 30 I = 2, M
               A( I, 1 ) = ZERO
   30       CONTINUE
            IF( M.GT.1 ) THEN
*
*              Form Q(2:m,2:m)
*
               CALL DORGQR( M-1, M-1, M-1, A( 2, 2 ), LDA, TAU, WORK,
     $                      LWORK, IINFO )
            END IF
         END IF
      ELSE
*
*        Form P', determined by a call to DGEBRD to reduce a k-by-n
*        matrix
*
         IF( K.LT.N ) THEN
*
*           If k < n, assume k <= m <= n
*
            CALL DORGLQ( M, N, K, A, LDA, TAU, WORK, LWORK, IINFO )
*
         ELSE
*
*           If k >= n, assume m = n
*
*           Shift the vectors which define the elementary reflectors one
*           row downward, and set the first row and column of P' to
*           those of the unit matrix
*
            A( 1, 1 ) = ONE
            DO 40 I = 2, N
               A( I, 1 ) = ZERO
   40       CONTINUE
            DO 60 J = 2, N
               DO 50 I = J - 1, 2, -1
                  A( I, J ) = A( I-1, J )
   50          CONTINUE
               A( 1, J ) = ZERO
   60       CONTINUE
            IF( N.GT.1 ) THEN
*
*              Form P'(2:n,2:n)
*
               CALL DORGLQ( N-1, N-1, N-1, A( 2, 2 ), LDA, TAU, WORK,
     $                      LWORK, IINFO )
            END IF
         END IF
      END IF
      WORK( 1 ) = LWKOPT
      RETURN
*
*     End of DORGBR
*
      END
      SUBROUTINE DORGL2( M, N, K, A, LDA, TAU, WORK, INFO )
*
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, M, N
*     ..
*     .. Array Arguments ..
CID      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
      DIMENSION   A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DORGL2 generates an m by n real matrix Q with orthonormal rows,
*  which is defined as the first m rows of a product of k elementary
*  reflectors of order n
*
*        Q  =  H(k) . . . H(2) H(1)
*
*  as returned by DGELQF.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix Q. M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix Q. N >= M.
*
*  K       (input) INTEGER
*          The number of elementary reflectors whose product defines the
*          matrix Q. M >= K >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the i-th row must contain the vector which defines
*          the elementary reflector H(i), for i = 1,2,...,k, as returned
*          by DGELQF in the first k rows of its array argument A.
*          On exit, the m-by-n matrix Q.
*
*  LDA     (input) INTEGER
*          The first dimension of the array A. LDA >= max(1,M).
*
*  TAU     (input) DOUBLE PRECISION array, dimension (K)
*          TAU(i) must contain the scalar factor of the elementary
*          reflector H(i), as returned by DGELQF.
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (M)
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument has an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
CID      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J, L
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLARF, DSCAL, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.M ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORGL2', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.LE.0 )
     $   RETURN
*
      IF( K.LT.M ) THEN
*
*        Initialise rows k+1:m to rows of the unit matrix
*
         DO 20 J = 1, N
            DO 10 L = K + 1, M
               A( L, J ) = ZERO
   10       CONTINUE
            IF( J.GT.K .AND. J.LE.M )
     $         A( J, J ) = ONE
   20    CONTINUE
      END IF
*
      DO 40 I = K, 1, -1
*
*        Apply H(i) to A(i:m,i:n) from the right
*
         IF( I.LT.N ) THEN
            IF( I.LT.M ) THEN
               A( I, I ) = ONE
               CALL DLARF( 'Right', M-I, N-I+1, A( I, I ), LDA,
     $                     TAU( I ), A( I+1, I ), LDA, WORK )
            END IF
            CALL DSCAL( N-I, -TAU( I ), A( I, I+1 ), LDA )
         END IF
         A( I, I ) = ONE - TAU( I )
*
*        Set A(i,1:i-1) to zero
*
         DO 30 L = 1, I - 1
            A( I, L ) = ZERO
   30    CONTINUE
   40 CONTINUE
      RETURN
*
*     End of DORGL2
*
      END
      SUBROUTINE DORGLQ( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
*
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, LWORK, M, N
*     ..
*     .. Array Arguments ..
CID      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
      DIMENSION A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DORGLQ generates an M-by-N real matrix Q with orthonormal rows,
*  which is defined as the first M rows of a product of K elementary
*  reflectors of order N
*
*        Q  =  H(k) . . . H(2) H(1)
*
*  as returned by DGELQF.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix Q. M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix Q. N >= M.
*
*  K       (input) INTEGER
*          The number of elementary reflectors whose product defines the
*          matrix Q. M >= K >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the i-th row must contain the vector which defines
*          the elementary reflector H(i), for i = 1,2,...,k, as returned
*          by DGELQF in the first k rows of its array argument A.
*          On exit, the M-by-N matrix Q.
*
*  LDA     (input) INTEGER
*          The first dimension of the array A. LDA >= max(1,M).
*
*  TAU     (input) DOUBLE PRECISION array, dimension (K)
*          TAU(i) must contain the scalar factor of the elementary
*          reflector H(i), as returned by DGELQF.
*
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK. LWORK >= max(1,M).
*          For optimum performance LWORK >= M*NB, where NB is
*          the optimal blocksize.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument has an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
CID      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IB, IINFO, IWS, J, KI, KK, L, LDWORK,
     $                   LWKOPT, NB, NBMIN, NX
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLARFB, DLARFT, DORGL2, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      NB = ILAENV( 1, 'DORGLQ', ' ', M, N, K, -1 )
      LWKOPT = MAX( 1, M )*NB
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.M ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      ELSE IF( LWORK.LT.MAX( 1, M ) .AND. .NOT.LQUERY ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORGLQ', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.LE.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
      NBMIN = 2
      NX = 0
      IWS = M
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
*
*        Determine when to cross over from blocked to unblocked code.
*
         NX = MAX( 0, ILAENV( 3, 'DORGLQ', ' ', M, N, K, -1 ) )
         IF( NX.LT.K ) THEN
*
*           Determine if workspace is large enough for blocked code.
*
            LDWORK = M
            IWS = LDWORK*NB
            IF( LWORK.LT.IWS ) THEN
*
*              Not enough workspace to use optimal NB:  reduce NB and
*              determine the minimum value of NB.
*
               NB = LWORK / LDWORK
               NBMIN = MAX( 2, ILAENV( 2, 'DORGLQ', ' ', M, N, K, -1 ) )
            END IF
         END IF
      END IF
*
      IF( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) THEN
*
*        Use blocked code after the last block.
*        The first kk rows are handled by the block method.
*
         KI = ( ( K-NX-1 ) / NB )*NB
         KK = MIN( K, KI+NB )
*
*        Set A(kk+1:m,1:kk) to zero.
*
         DO 20 J = 1, KK
            DO 10 I = KK + 1, M
               A( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
      ELSE
         KK = 0
      END IF
*
*     Use unblocked code for the last or only block.
*
      IF( KK.LT.M )
     $   CALL DORGL2( M-KK, N-KK, K-KK, A( KK+1, KK+1 ), LDA,
     $                TAU( KK+1 ), WORK, IINFO )
*
      IF( KK.GT.0 ) THEN
*
*        Use blocked code
*
         DO 50 I = KI + 1, 1, -NB
            IB = MIN( NB, K-I+1 )
            IF( I+IB.LE.M ) THEN
*
*              Form the triangular factor of the block reflector
*              H = H(i) H(i+1) . . . H(i+ib-1)
*
               CALL DLARFT( 'Forward', 'Rowwise', N-I+1, IB, A( I, I ),
     $                      LDA, TAU( I ), WORK, LDWORK )
*
*              Apply H' to A(i+ib:m,i:n) from the right
*
               CALL DLARFB( 'Right', 'Transpose', 'Forward', 'Rowwise',
     $                      M-I-IB+1, N-I+1, IB, A( I, I ), LDA, WORK,
     $                      LDWORK, A( I+IB, I ), LDA, WORK( IB+1 ),
     $                      LDWORK )
            END IF
*
*           Apply H' to columns i:n of current block
*
            CALL DORGL2( IB, N-I+1, IB, A( I, I ), LDA, TAU( I ), WORK,
     $                   IINFO )
*
*           Set columns 1:i-1 of current block to zero
*
            DO 40 J = 1, I - 1
               DO 30 L = I, I + IB - 1
                  A( L, J ) = ZERO
   30          CONTINUE
   40       CONTINUE
   50    CONTINUE
      END IF
*
      WORK( 1 ) = IWS
      RETURN
*
*     End of DORGLQ
*
      END
      SUBROUTINE DORGQR( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
*
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, LWORK, M, N
*     ..
*     .. Array Arguments ..
CID      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
      DIMENSION   A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DORGQR generates an M-by-N real matrix Q with orthonormal columns,
*  which is defined as the first N columns of a product of K elementary
*  reflectors of order M
*
*        Q  =  H(1) H(2) . . . H(k)
*
*  as returned by DGEQRF.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix Q. M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix Q. M >= N >= 0.
*
*  K       (input) INTEGER
*          The number of elementary reflectors whose product defines the
*          matrix Q. N >= K >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the i-th column must contain the vector which
*          defines the elementary reflector H(i), for i = 1,2,...,k, as
*          returned by DGEQRF in the first k columns of its array
*          argument A.
*          On exit, the M-by-N matrix Q.
*
*  LDA     (input) INTEGER
*          The first dimension of the array A. LDA >= max(1,M).
*
*  TAU     (input) DOUBLE PRECISION array, dimension (K)
*          TAU(i) must contain the scalar factor of the elementary
*          reflector H(i), as returned by DGEQRF.
*
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK. LWORK >= max(1,N).
*          For optimum performance LWORK >= N*NB, where NB is the
*          optimal blocksize.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument has an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
CID      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IB, IINFO, IWS, J, KI, KK, L, LDWORK,
     $                   LWKOPT, NB, NBMIN, NX
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLARFB, DLARFT, DORG2R, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      NB = ILAENV( 1, 'DORGQR', ' ', M, N, K, -1 )
      LWKOPT = MAX( 1, N )*NB
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      ELSE IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORGQR', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.LE.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
      NBMIN = 2
      NX = 0
      IWS = N
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
*
*        Determine when to cross over from blocked to unblocked code.
*
         NX = MAX( 0, ILAENV( 3, 'DORGQR', ' ', M, N, K, -1 ) )
         IF( NX.LT.K ) THEN
*
*           Determine if workspace is large enough for blocked code.
*
            LDWORK = N
            IWS = LDWORK*NB
            IF( LWORK.LT.IWS ) THEN
*
*              Not enough workspace to use optimal NB:  reduce NB and
*              determine the minimum value of NB.
*
               NB = LWORK / LDWORK
               NBMIN = MAX( 2, ILAENV( 2, 'DORGQR', ' ', M, N, K, -1 ) )
            END IF
         END IF
      END IF
*
      IF( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) THEN
*
*        Use blocked code after the last block.
*        The first kk columns are handled by the block method.
*
         KI = ( ( K-NX-1 ) / NB )*NB
         KK = MIN( K, KI+NB )
*
*        Set A(1:kk,kk+1:n) to zero.
*
         DO 20 J = KK + 1, N
            DO 10 I = 1, KK
               A( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
      ELSE
         KK = 0
      END IF
*
*     Use unblocked code for the last or only block.
*
      IF( KK.LT.N )
     $   CALL DORG2R( M-KK, N-KK, K-KK, A( KK+1, KK+1 ), LDA,
     $                TAU( KK+1 ), WORK, IINFO )
*
      IF( KK.GT.0 ) THEN
*
*        Use blocked code
*
         DO 50 I = KI + 1, 1, -NB
            IB = MIN( NB, K-I+1 )
            IF( I+IB.LE.N ) THEN
*
*              Form the triangular factor of the block reflector
*              H = H(i) H(i+1) . . . H(i+ib-1)
*
               CALL DLARFT( 'Forward', 'Columnwise', M-I+1, IB,
     $                      A( I, I ), LDA, TAU( I ), WORK, LDWORK )
*
*              Apply H to A(i:m,i+ib:n) from the left
*
               CALL DLARFB( 'Left', 'No transpose', 'Forward',
     $                      'Columnwise', M-I+1, N-I-IB+1, IB,
     $                      A( I, I ), LDA, WORK, LDWORK, A( I, I+IB ),
     $                      LDA, WORK( IB+1 ), LDWORK )
            END IF
*
*           Apply H to rows i:m of current block
*
            CALL DORG2R( M-I+1, IB, IB, A( I, I ), LDA, TAU( I ), WORK,
     $                   IINFO )
*
*           Set rows 1:i-1 of current block to zero
*
            DO 40 J = I, I + IB - 1
               DO 30 L = 1, I - 1
                  A( L, J ) = ZERO
   30          CONTINUE
   40       CONTINUE
   50    CONTINUE
      END IF
*
      WORK( 1 ) = IWS
      RETURN
*
*     End of DORGQR
*
      END
      SUBROUTINE DORM2R( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,
     $                   WORK, INFO )
*
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER            INFO, K, LDA, LDC, M, N
*     ..
*     .. Array Arguments ..
CID      DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
      DIMENSION  A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DORM2R overwrites the general real m by n matrix C with
*
*        Q * C  if SIDE = 'L' and TRANS = 'N', or
*
*        Q'* C  if SIDE = 'L' and TRANS = 'T', or
*
*        C * Q  if SIDE = 'R' and TRANS = 'N', or
*
*        C * Q' if SIDE = 'R' and TRANS = 'T',
*
*  where Q is a real orthogonal matrix defined as the product of k
*  elementary reflectors
*
*        Q = H(1) H(2) . . . H(k)
*
*  as returned by DGEQRF. Q is of order m if SIDE = 'L' and of order n
*  if SIDE = 'R'.
*
*  Arguments
*  =========
*
*  SIDE    (input) CHARACTER*1
*          = 'L': apply Q or Q' from the Left
*          = 'R': apply Q or Q' from the Right
*
*  TRANS   (input) CHARACTER*1
*          = 'N': apply Q  (No transpose)
*          = 'T': apply Q' (Transpose)
*
*  M       (input) INTEGER
*          The number of rows of the matrix C. M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix C. N >= 0.
*
*  K       (input) INTEGER
*          The number of elementary reflectors whose product defines
*          the matrix Q.
*          If SIDE = 'L', M >= K >= 0;
*          if SIDE = 'R', N >= K >= 0.
*
*  A       (input) DOUBLE PRECISION array, dimension (LDA,K)
*          The i-th column must contain the vector which defines the
*          elementary reflector H(i), for i = 1,2,...,k, as returned by
*          DGEQRF in the first k columns of its array argument A.
*          A is modified by the routine but restored on exit.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.
*          If SIDE = 'L', LDA >= max(1,M);
*          if SIDE = 'R', LDA >= max(1,N).
*
*  TAU     (input) DOUBLE PRECISION array, dimension (K)
*          TAU(i) must contain the scalar factor of the elementary
*          reflector H(i), as returned by DGEQRF.
*
*  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
*          On entry, the m by n matrix C.
*          On exit, C is overwritten by Q*C or Q'*C or C*Q' or C*Q.
*
*  LDC     (input) INTEGER
*          The leading dimension of the array C. LDC >= max(1,M).
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension
*                                   (N) if SIDE = 'L',
*                                   (M) if SIDE = 'R'
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
CID      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LEFT, NOTRAN
      INTEGER            I, I1, I2, I3, IC, JC, MI, NI, NQ
CID      DOUBLE PRECISION   AII
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLARF, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      LEFT = LSAME( SIDE, 'L' )
      NOTRAN = LSAME( TRANS, 'N' )
*
*     NQ is the order of Q
*
      IF( LEFT ) THEN
         NQ = M
      ELSE
         NQ = N
      END IF
      IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( K.LT.0 .OR. K.GT.NQ ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, NQ ) ) THEN
         INFO = -7
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORM2R', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 )
     $   RETURN
*
      IF( ( LEFT .AND. .NOT.NOTRAN ) .OR. ( .NOT.LEFT .AND. NOTRAN ) )
     $     THEN
         I1 = 1
         I2 = K
         I3 = 1
      ELSE
         I1 = K
         I2 = 1
         I3 = -1
      END IF
*
      IF( LEFT ) THEN
         NI = N
         JC = 1
      ELSE
         MI = M
         IC = 1
      END IF
*
      DO 10 I = I1, I2, I3
         IF( LEFT ) THEN
*
*           H(i) is applied to C(i:m,1:n)
*
            MI = M - I + 1
            IC = I
         ELSE
*
*           H(i) is applied to C(1:m,i:n)
*
            NI = N - I + 1
            JC = I
         END IF
*
*        Apply H(i)
*
         AII = A( I, I )
         A( I, I ) = ONE
         CALL DLARF( SIDE, MI, NI, A( I, I ), 1, TAU( I ), C( IC, JC ),
     $               LDC, WORK )
         A( I, I ) = AII
   10 CONTINUE
      RETURN
*
*     End of DORM2R
*
      END
      SUBROUTINE DORMBR( VECT, SIDE, TRANS, M, N, K, A, LDA, TAU, C,
     $                   LDC, WORK, LWORK, INFO )
*
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS, VECT
      INTEGER            INFO, K, LDA, LDC, LWORK, M, N
*     ..
*     .. Array Arguments ..
CID      DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
      DIMENSION  A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  If VECT = 'Q', DORMBR overwrites the general real M-by-N matrix C
*  with
*                  SIDE = 'L'     SIDE = 'R'
*  TRANS = 'N':      Q * C          C * Q
*  TRANS = 'T':      Q**T * C       C * Q**T
*
*  If VECT = 'P', DORMBR overwrites the general real M-by-N matrix C
*  with
*                  SIDE = 'L'     SIDE = 'R'
*  TRANS = 'N':      P * C          C * P
*  TRANS = 'T':      P**T * C       C * P**T
*
*  Here Q and P**T are the orthogonal matrices determined by DGEBRD when
*  reducing a real matrix A to bidiagonal form: A = Q * B * P**T. Q and
*  P**T are defined as products of elementary reflectors H(i) and G(i)
*  respectively.
*
*  Let nq = m if SIDE = 'L' and nq = n if SIDE = 'R'. Thus nq is the
*  order of the orthogonal matrix Q or P**T that is applied.
*
*  If VECT = 'Q', A is assumed to have been an NQ-by-K matrix:
*  if nq >= k, Q = H(1) H(2) . . . H(k);
*  if nq < k, Q = H(1) H(2) . . . H(nq-1).
*
*  If VECT = 'P', A is assumed to have been a K-by-NQ matrix:
*  if k < nq, P = G(1) G(2) . . . G(k);
*  if k >= nq, P = G(1) G(2) . . . G(nq-1).
*
*  Arguments
*  =========
*
*  VECT    (input) CHARACTER*1
*          = 'Q': apply Q or Q**T;
*          = 'P': apply P or P**T.
*
*  SIDE    (input) CHARACTER*1
*          = 'L': apply Q, Q**T, P or P**T from the Left;
*          = 'R': apply Q, Q**T, P or P**T from the Right.
*
*  TRANS   (input) CHARACTER*1
*          = 'N':  No transpose, apply Q  or P;
*          = 'T':  Transpose, apply Q**T or P**T.
*
*  M       (input) INTEGER
*          The number of rows of the matrix C. M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix C. N >= 0.
*
*  K       (input) INTEGER
*          If VECT = 'Q', the number of columns in the original
*          matrix reduced by DGEBRD.
*          If VECT = 'P', the number of rows in the original
*          matrix reduced by DGEBRD.
*          K >= 0.
*
*  A       (input) DOUBLE PRECISION array, dimension
*                                (LDA,min(nq,K)) if VECT = 'Q'
*                                (LDA,nq)        if VECT = 'P'
*          The vectors which define the elementary reflectors H(i) and
*          G(i), whose products determine the matrices Q and P, as
*          returned by DGEBRD.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.
*          If VECT = 'Q', LDA >= max(1,nq);
*          if VECT = 'P', LDA >= max(1,min(nq,K)).
*
*  TAU     (input) DOUBLE PRECISION array, dimension (min(nq,K))
*          TAU(i) must contain the scalar factor of the elementary
*          reflector H(i) or G(i) which determines Q or P, as returned
*          by DGEBRD in the array argument TAUQ or TAUP.
*
*  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
*          On entry, the M-by-N matrix C.
*          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q
*          or P*C or P**T*C or C*P or C*P**T.
*
*  LDC     (input) INTEGER
*          The leading dimension of the array C. LDC >= max(1,M).
*
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.
*          If SIDE = 'L', LWORK >= max(1,N);
*          if SIDE = 'R', LWORK >= max(1,M).
*          For optimum performance LWORK >= N*NB if SIDE = 'L', and
*          LWORK >= M*NB if SIDE = 'R', where NB is the optimal
*          blocksize.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            APPLYQ, LEFT, LQUERY, NOTRAN
      CHARACTER          TRANST
      INTEGER            I1, I2, IINFO, LWKOPT, MI, NB, NI, NQ, NW
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL           DORMLQ, DORMQR, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      APPLYQ = LSAME( VECT, 'Q' )
      LEFT = LSAME( SIDE, 'L' )
      NOTRAN = LSAME( TRANS, 'N' )
      LQUERY = ( LWORK.EQ.-1 )
*
*     NQ is the order of Q or P and NW is the minimum dimension of WORK
*
      IF( LEFT ) THEN
         NQ = M
         NW = N
      ELSE
         NQ = N
         NW = M
      END IF
      IF( .NOT.APPLYQ .AND. .NOT.LSAME( VECT, 'P' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
         INFO = -2
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) ) THEN
         INFO = -3
      ELSE IF( M.LT.0 ) THEN
         INFO = -4
      ELSE IF( N.LT.0 ) THEN
         INFO = -5
      ELSE IF( K.LT.0 ) THEN
         INFO = -6
      ELSE IF( ( APPLYQ .AND. LDA.LT.MAX( 1, NQ ) ) .OR.
     $         ( .NOT.APPLYQ .AND. LDA.LT.MAX( 1, MIN( NQ, K ) ) ) )
     $          THEN
         INFO = -8
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -11
      ELSE IF( LWORK.LT.MAX( 1, NW ) .AND. .NOT.LQUERY ) THEN
         INFO = -13
      END IF
*
      IF( INFO.EQ.0 ) THEN
         IF( APPLYQ ) THEN
            IF( LEFT ) THEN
               NB = ILAENV( 1, 'DORMQR', SIDE // TRANS, M-1, N, M-1,
     $              -1 )
            ELSE
               NB = ILAENV( 1, 'DORMQR', SIDE // TRANS, M, N-1, N-1,
     $              -1 )
            END IF
         ELSE
            IF( LEFT ) THEN
               NB = ILAENV( 1, 'DORMLQ', SIDE // TRANS, M-1, N, M-1,
     $              -1 )
            ELSE
               NB = ILAENV( 1, 'DORMLQ', SIDE // TRANS, M, N-1, N-1,
     $              -1 )
            END IF
         END IF
         LWKOPT = MAX( 1, NW )*NB
         WORK( 1 ) = LWKOPT
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORMBR', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      WORK( 1 ) = 1
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
*
      IF( APPLYQ ) THEN
*
*        Apply Q
*
         IF( NQ.GE.K ) THEN
*
*           Q was determined by a call to DGEBRD with nq >= k
*
            CALL DORMQR( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,
     $                   WORK, LWORK, IINFO )
         ELSE IF( NQ.GT.1 ) THEN
*
*           Q was determined by a call to DGEBRD with nq < k
*
            IF( LEFT ) THEN
               MI = M - 1
               NI = N
               I1 = 2
               I2 = 1
            ELSE
               MI = M
               NI = N - 1
               I1 = 1
               I2 = 2
            END IF
            CALL DORMQR( SIDE, TRANS, MI, NI, NQ-1, A( 2, 1 ), LDA, TAU,
     $                   C( I1, I2 ), LDC, WORK, LWORK, IINFO )
         END IF
      ELSE
*
*        Apply P
*
         IF( NOTRAN ) THEN
            TRANST = 'T'
         ELSE
            TRANST = 'N'
         END IF
         IF( NQ.GT.K ) THEN
*
*           P was determined by a call to DGEBRD with nq > k
*
            CALL DORMLQ( SIDE, TRANST, M, N, K, A, LDA, TAU, C, LDC,
     $                   WORK, LWORK, IINFO )
         ELSE IF( NQ.GT.1 ) THEN
*
*           P was determined by a call to DGEBRD with nq <= k
*
            IF( LEFT ) THEN
               MI = M - 1
               NI = N
               I1 = 2
               I2 = 1
            ELSE
               MI = M
               NI = N - 1
               I1 = 1
               I2 = 2
            END IF
            CALL DORMLQ( SIDE, TRANST, MI, NI, NQ-1, A( 1, 2 ), LDA,
     $                   TAU, C( I1, I2 ), LDC, WORK, LWORK, IINFO )
         END IF
      END IF
      WORK( 1 ) = LWKOPT
      RETURN
*
*     End of DORMBR
*
      END
      SUBROUTINE DORML2( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,
     $                   WORK, INFO )
*
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER            INFO, K, LDA, LDC, M, N
*     ..
*     .. Array Arguments ..
CID      DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
      DIMENSION  A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DORML2 overwrites the general real m by n matrix C with
*
*        Q * C  if SIDE = 'L' and TRANS = 'N', or
*
*        Q'* C  if SIDE = 'L' and TRANS = 'T', or
*
*        C * Q  if SIDE = 'R' and TRANS = 'N', or
*
*        C * Q' if SIDE = 'R' and TRANS = 'T',
*
*  where Q is a real orthogonal matrix defined as the product of k
*  elementary reflectors
*
*        Q = H(k) . . . H(2) H(1)
*
*  as returned by DGELQF. Q is of order m if SIDE = 'L' and of order n
*  if SIDE = 'R'.
*
*  Arguments
*  =========
*
*  SIDE    (input) CHARACTER*1
*          = 'L': apply Q or Q' from the Left
*          = 'R': apply Q or Q' from the Right
*
*  TRANS   (input) CHARACTER*1
*          = 'N': apply Q  (No transpose)
*          = 'T': apply Q' (Transpose)
*
*  M       (input) INTEGER
*          The number of rows of the matrix C. M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix C. N >= 0.
*
*  K       (input) INTEGER
*          The number of elementary reflectors whose product defines
*          the matrix Q.
*          If SIDE = 'L', M >= K >= 0;
*          if SIDE = 'R', N >= K >= 0.
*
*  A       (input) DOUBLE PRECISION array, dimension
*                               (LDA,M) if SIDE = 'L',
*                               (LDA,N) if SIDE = 'R'
*          The i-th row must contain the vector which defines the
*          elementary reflector H(i), for i = 1,2,...,k, as returned by
*          DGELQF in the first k rows of its array argument A.
*          A is modified by the routine but restored on exit.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A. LDA >= max(1,K).
*
*  TAU     (input) DOUBLE PRECISION array, dimension (K)
*          TAU(i) must contain the scalar factor of the elementary
*          reflector H(i), as returned by DGELQF.
*
*  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
*          On entry, the m by n matrix C.
*          On exit, C is overwritten by Q*C or Q'*C or C*Q' or C*Q.
*
*  LDC     (input) INTEGER
*          The leading dimension of the array C. LDC >= max(1,M).
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension
*                                   (N) if SIDE = 'L',
*                                   (M) if SIDE = 'R'
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
CID      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LEFT, NOTRAN
      INTEGER            I, I1, I2, I3, IC, JC, MI, NI, NQ
CID      DOUBLE PRECISION   AII
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLARF, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      LEFT = LSAME( SIDE, 'L' )
      NOTRAN = LSAME( TRANS, 'N' )
*
*     NQ is the order of Q
*
      IF( LEFT ) THEN
         NQ = M
      ELSE
         NQ = N
      END IF
      IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( K.LT.0 .OR. K.GT.NQ ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, K ) ) THEN
         INFO = -7
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORML2', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 )
     $   RETURN
*
      IF( ( LEFT .AND. NOTRAN ) .OR. ( .NOT.LEFT .AND. .NOT.NOTRAN ) )
     $     THEN
         I1 = 1
         I2 = K
         I3 = 1
      ELSE
         I1 = K
         I2 = 1
         I3 = -1
      END IF
*
      IF( LEFT ) THEN
         NI = N
         JC = 1
      ELSE
         MI = M
         IC = 1
      END IF
*
      DO 10 I = I1, I2, I3
         IF( LEFT ) THEN
*
*           H(i) is applied to C(i:m,1:n)
*
            MI = M - I + 1
            IC = I
         ELSE
*
*           H(i) is applied to C(1:m,i:n)
*
            NI = N - I + 1
            JC = I
         END IF
*
*        Apply H(i)
*
         AII = A( I, I )
         A( I, I ) = ONE
         CALL DLARF( SIDE, MI, NI, A( I, I ), LDA, TAU( I ),
     $               C( IC, JC ), LDC, WORK )
         A( I, I ) = AII
   10 CONTINUE
      RETURN
*
*     End of DORML2
*
      END
      SUBROUTINE DORMLQ( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,
     $                   WORK, LWORK, INFO )
*
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER            INFO, K, LDA, LDC, LWORK, M, N
*     ..
*     .. Array Arguments ..
CID      DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
      DIMENSION A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DORMLQ overwrites the general real M-by-N matrix C with
*
*                  SIDE = 'L'     SIDE = 'R'
*  TRANS = 'N':      Q * C          C * Q
*  TRANS = 'T':      Q**T * C       C * Q**T
*
*  where Q is a real orthogonal matrix defined as the product of k
*  elementary reflectors
*
*        Q = H(k) . . . H(2) H(1)
*
*  as returned by DGELQF. Q is of order M if SIDE = 'L' and of order N
*  if SIDE = 'R'.
*
*  Arguments
*  =========
*
*  SIDE    (input) CHARACTER*1
*          = 'L': apply Q or Q**T from the Left;
*          = 'R': apply Q or Q**T from the Right.
*
*  TRANS   (input) CHARACTER*1
*          = 'N':  No transpose, apply Q;
*          = 'T':  Transpose, apply Q**T.
*
*  M       (input) INTEGER
*          The number of rows of the matrix C. M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix C. N >= 0.
*
*  K       (input) INTEGER
*          The number of elementary reflectors whose product defines
*          the matrix Q.
*          If SIDE = 'L', M >= K >= 0;
*          if SIDE = 'R', N >= K >= 0.
*
*  A       (input) DOUBLE PRECISION array, dimension
*                               (LDA,M) if SIDE = 'L',
*                               (LDA,N) if SIDE = 'R'
*          The i-th row must contain the vector which defines the
*          elementary reflector H(i), for i = 1,2,...,k, as returned by
*          DGELQF in the first k rows of its array argument A.
*          A is modified by the routine but restored on exit.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A. LDA >= max(1,K).
*
*  TAU     (input) DOUBLE PRECISION array, dimension (K)
*          TAU(i) must contain the scalar factor of the elementary
*          reflector H(i), as returned by DGELQF.
*
*  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
*          On entry, the M-by-N matrix C.
*          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.
*
*  LDC     (input) INTEGER
*          The leading dimension of the array C. LDC >= max(1,M).
*
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.
*          If SIDE = 'L', LWORK >= max(1,N);
*          if SIDE = 'R', LWORK >= max(1,M).
*          For optimum performance LWORK >= N*NB if SIDE = 'L', and
*          LWORK >= M*NB if SIDE = 'R', where NB is the optimal
*          blocksize.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            NBMAX, LDT
      PARAMETER          ( NBMAX = 64, LDT = NBMAX+1 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LEFT, LQUERY, NOTRAN
      CHARACTER          TRANST
      INTEGER            I, I1, I2, I3, IB, IC, IINFO, IWS, JC, LDWORK,
     $                   LWKOPT, MI, NB, NBMIN, NI, NQ, NW
*     ..
*     .. Local Arrays ..
CID      DOUBLE PRECISION   T( LDT, NBMAX )
      DIMENSION   T( LDT, NBMAX )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLARFB, DLARFT, DORML2, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      LEFT = LSAME( SIDE, 'L' )
      NOTRAN = LSAME( TRANS, 'N' )
      LQUERY = ( LWORK.EQ.-1 )
*
*     NQ is the order of Q and NW is the minimum dimension of WORK
*
      IF( LEFT ) THEN
         NQ = M
         NW = N
      ELSE
         NQ = N
         NW = M
      END IF
      IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( K.LT.0 .OR. K.GT.NQ ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, K ) ) THEN
         INFO = -7
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -10
      ELSE IF( LWORK.LT.MAX( 1, NW ) .AND. .NOT.LQUERY ) THEN
         INFO = -12
      END IF
*
      IF( INFO.EQ.0 ) THEN
*
*        Determine the block size.  NB may be at most NBMAX, where NBMAX
*        is used to define the local array T.
*
         NB = MIN( NBMAX, ILAENV( 1, 'DORMLQ', SIDE // TRANS, M, N, K,
     $        -1 ) )
         LWKOPT = MAX( 1, NW )*NB
         WORK( 1 ) = LWKOPT
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORMLQ', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
      NBMIN = 2
      LDWORK = NW
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
         IWS = NW*NB
         IF( LWORK.LT.IWS ) THEN
            NB = LWORK / LDWORK
            NBMIN = MAX( 2, ILAENV( 2, 'DORMLQ', SIDE // TRANS, M, N, K,
     $              -1 ) )
         END IF
      ELSE
         IWS = NW
      END IF
*
      IF( NB.LT.NBMIN .OR. NB.GE.K ) THEN
*
*        Use unblocked code
*
         CALL DORML2( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK,
     $                IINFO )
      ELSE
*
*        Use blocked code
*
         IF( ( LEFT .AND. NOTRAN ) .OR.
     $       ( .NOT.LEFT .AND. .NOT.NOTRAN ) ) THEN
            I1 = 1
            I2 = K
            I3 = NB
         ELSE
            I1 = ( ( K-1 ) / NB )*NB + 1
            I2 = 1
            I3 = -NB
         END IF
*
         IF( LEFT ) THEN
            NI = N
            JC = 1
         ELSE
            MI = M
            IC = 1
         END IF
*
         IF( NOTRAN ) THEN
            TRANST = 'T'
         ELSE
            TRANST = 'N'
         END IF
*
         DO 10 I = I1, I2, I3
            IB = MIN( NB, K-I+1 )
*
*           Form the triangular factor of the block reflector
*           H = H(i) H(i+1) . . . H(i+ib-1)
*
            CALL DLARFT( 'Forward', 'Rowwise', NQ-I+1, IB, A( I, I ),
     $                   LDA, TAU( I ), T, LDT )
            IF( LEFT ) THEN
*
*              H or H' is applied to C(i:m,1:n)
*
               MI = M - I + 1
               IC = I
            ELSE
*
*              H or H' is applied to C(1:m,i:n)
*
               NI = N - I + 1
               JC = I
            END IF
*
*           Apply H or H'
*
            CALL DLARFB( SIDE, TRANST, 'Forward', 'Rowwise', MI, NI, IB,
     $                   A( I, I ), LDA, T, LDT, C( IC, JC ), LDC, WORK,
     $                   LDWORK )
   10    CONTINUE
      END IF
      WORK( 1 ) = LWKOPT
      RETURN
*
*     End of DORMLQ
*
      END
      SUBROUTINE DORMQR( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,
     $                   WORK, LWORK, INFO )
*
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER            INFO, K, LDA, LDC, LWORK, M, N
*     ..
*     .. Array Arguments ..
CID      DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
      DIMENSION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DORMQR overwrites the general real M-by-N matrix C with
*
*                  SIDE = 'L'     SIDE = 'R'
*  TRANS = 'N':      Q * C          C * Q
*  TRANS = 'T':      Q**T * C       C * Q**T
*
*  where Q is a real orthogonal matrix defined as the product of k
*  elementary reflectors
*
*        Q = H(1) H(2) . . . H(k)
*
*  as returned by DGEQRF. Q is of order M if SIDE = 'L' and of order N
*  if SIDE = 'R'.
*
*  Arguments
*  =========
*
*  SIDE    (input) CHARACTER*1
*          = 'L': apply Q or Q**T from the Left;
*          = 'R': apply Q or Q**T from the Right.
*
*  TRANS   (input) CHARACTER*1
*          = 'N':  No transpose, apply Q;
*          = 'T':  Transpose, apply Q**T.
*
*  M       (input) INTEGER
*          The number of rows of the matrix C. M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix C. N >= 0.
*
*  K       (input) INTEGER
*          The number of elementary reflectors whose product defines
*          the matrix Q.
*          If SIDE = 'L', M >= K >= 0;
*          if SIDE = 'R', N >= K >= 0.
*
*  A       (input) DOUBLE PRECISION array, dimension (LDA,K)
*          The i-th column must contain the vector which defines the
*          elementary reflector H(i), for i = 1,2,...,k, as returned by
*          DGEQRF in the first k columns of its array argument A.
*          A is modified by the routine but restored on exit.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.
*          If SIDE = 'L', LDA >= max(1,M);
*          if SIDE = 'R', LDA >= max(1,N).
*
*  TAU     (input) DOUBLE PRECISION array, dimension (K)
*          TAU(i) must contain the scalar factor of the elementary
*          reflector H(i), as returned by DGEQRF.
*
*  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
*          On entry, the M-by-N matrix C.
*          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.
*
*  LDC     (input) INTEGER
*          The leading dimension of the array C. LDC >= max(1,M).
*
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.
*          If SIDE = 'L', LWORK >= max(1,N);
*          if SIDE = 'R', LWORK >= max(1,M).
*          For optimum performance LWORK >= N*NB if SIDE = 'L', and
*          LWORK >= M*NB if SIDE = 'R', where NB is the optimal
*          blocksize.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            NBMAX, LDT
      PARAMETER          ( NBMAX = 64, LDT = NBMAX+1 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LEFT, LQUERY, NOTRAN
      INTEGER            I, I1, I2, I3, IB, IC, IINFO, IWS, JC, LDWORK,
     $                   LWKOPT, MI, NB, NBMIN, NI, NQ, NW
*     ..
*     .. Local Arrays ..
CID      DOUBLE PRECISION   T( LDT, NBMAX )
      DIMENSION   T( LDT, NBMAX )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLARFB, DLARFT, DORM2R, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      LEFT = LSAME( SIDE, 'L' )
      NOTRAN = LSAME( TRANS, 'N' )
      LQUERY = ( LWORK.EQ.-1 )
*
*     NQ is the order of Q and NW is the minimum dimension of WORK
*
      IF( LEFT ) THEN
         NQ = M
         NW = N
      ELSE
         NQ = N
         NW = M
      END IF
      IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( K.LT.0 .OR. K.GT.NQ ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, NQ ) ) THEN
         INFO = -7
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -10
      ELSE IF( LWORK.LT.MAX( 1, NW ) .AND. .NOT.LQUERY ) THEN
         INFO = -12
      END IF
*
      IF( INFO.EQ.0 ) THEN
*
*        Determine the block size.  NB may be at most NBMAX, where NBMAX
*        is used to define the local array T.
*
         NB = MIN( NBMAX, ILAENV( 1, 'DORMQR', SIDE // TRANS, M, N, K,
     $        -1 ) )
         LWKOPT = MAX( 1, NW )*NB
         WORK( 1 ) = LWKOPT
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORMQR', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
      NBMIN = 2
      LDWORK = NW
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
         IWS = NW*NB
         IF( LWORK.LT.IWS ) THEN
            NB = LWORK / LDWORK
            NBMIN = MAX( 2, ILAENV( 2, 'DORMQR', SIDE // TRANS, M, N, K,
     $              -1 ) )
         END IF
      ELSE
         IWS = NW
      END IF
*
      IF( NB.LT.NBMIN .OR. NB.GE.K ) THEN
*
*        Use unblocked code
*
         CALL DORM2R( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK,
     $                IINFO )
      ELSE
*
*        Use blocked code
*
         IF( ( LEFT .AND. .NOT.NOTRAN ) .OR.
     $       ( .NOT.LEFT .AND. NOTRAN ) ) THEN
            I1 = 1
            I2 = K
            I3 = NB
         ELSE
            I1 = ( ( K-1 ) / NB )*NB + 1
            I2 = 1
            I3 = -NB
         END IF
*
         IF( LEFT ) THEN
            NI = N
            JC = 1
         ELSE
            MI = M
            IC = 1
         END IF
*
         DO 10 I = I1, I2, I3
            IB = MIN( NB, K-I+1 )
*
*           Form the triangular factor of the block reflector
*           H = H(i) H(i+1) . . . H(i+ib-1)
*
            CALL DLARFT( 'Forward', 'Columnwise', NQ-I+1, IB, A( I, I ),
     $                   LDA, TAU( I ), T, LDT )
            IF( LEFT ) THEN
*
*              H or H' is applied to C(i:m,1:n)
*
               MI = M - I + 1
               IC = I
            ELSE
*
*              H or H' is applied to C(1:m,i:n)
*
               NI = N - I + 1
               JC = I
            END IF
*
*           Apply H or H'
*
            CALL DLARFB( SIDE, TRANS, 'Forward', 'Columnwise', MI, NI,
     $                   IB, A( I, I ), LDA, T, LDT, C( IC, JC ), LDC,
     $                   WORK, LDWORK )
   10    CONTINUE
      END IF
      WORK( 1 ) = LWKOPT
      RETURN
*
*     End of DORMQR
*
      END

      SUBROUTINE DCOPY(N,DX,INCX,DY,INCY)
*     .. Scalar Arguments ..
      INTEGER INCX,INCY,N
*     ..
*     .. Array Arguments ..
CID      DOUBLE PRECISION DX(*),DY(*)
      DIMENSION DX(*),DY(*)
*     ..
*
*  Purpose
*  =======
*
*     copies a vector, x, to a vector, y.
*     uses unrolled loops for increments equal to one.
*     jack dongarra, linpack, 3/11/78.
*     modified 12/3/93, array(1) declarations changed to array(*)
*
*
*     .. Local Scalars ..
      INTEGER I,IX,IY,M,MP1
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MOD
*     ..
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20
*
*        code for unequal increments or equal increments
*          not equal to 1
*
      IX = 1
      IY = 1
      IF (INCX.LT.0) IX = (-N+1)*INCX + 1
      IF (INCY.LT.0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
          DY(IY) = DX(IX)
          IX = IX + INCX
          IY = IY + INCY
   10 CONTINUE
      RETURN
*
*        code for both increments equal to 1
*
*
*        clean-up loop
*
   20 M = MOD(N,7)
      IF (M.EQ.0) GO TO 40
      DO 30 I = 1,M
          DY(I) = DX(I)
   30 CONTINUE
      IF (N.LT.7) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,7
          DY(I) = DX(I)
          DY(I+1) = DX(I+1)
          DY(I+2) = DX(I+2)
          DY(I+3) = DX(I+3)
          DY(I+4) = DX(I+4)
          DY(I+5) = DX(I+5)
          DY(I+6) = DX(I+6)
   50 CONTINUE
      RETURN
      END
      
      
      SUBROUTINE DROT(N,DX,INCX,DY,INCY,C,S)
*     .. Scalar Arguments ..
CID      DOUBLE PRECISION C,S
      INTEGER INCX,INCY,N
*     ..
*     .. Array Arguments ..
CID      DOUBLE PRECISION DX(*),DY(*)
      DIMENSION  DX(*),DY(*)
*     ..
*
*  Purpose
*  =======
*
*     applies a plane rotation.
*     jack dongarra, linpack, 3/11/78.
*     modified 12/3/93, array(1) declarations changed to array(*)
*
*
*     .. Local Scalars ..
CID      DOUBLE PRECISION DTEMP
      INTEGER I,IX,IY
*     ..
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20
*
*       code for unequal increments or equal increments not equal
*         to 1
*
      IX = 1
      IY = 1
      IF (INCX.LT.0) IX = (-N+1)*INCX + 1
      IF (INCY.LT.0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
          DTEMP = C*DX(IX) + S*DY(IY)
          DY(IY) = C*DY(IY) - S*DX(IX)
          DX(IX) = DTEMP
          IX = IX + INCX
          IY = IY + INCY
   10 CONTINUE
      RETURN
*
*       code for both increments equal to 1
*
   20 DO 30 I = 1,N
          DTEMP = C*DX(I) + S*DY(I)
          DY(I) = C*DY(I) - S*DX(I)
          DX(I) = DTEMP
   30 CONTINUE
      RETURN
      END



CID      DOUBLE PRECISION FUNCTION DNRM2(N,X,INCX)
      FUNCTION DNRM2(N,X,INCX)
*     .. Scalar Arguments ..
      INTEGER INCX,N
*     ..
*     .. Array Arguments ..
CID      DOUBLE PRECISION X(*)
      DIMENSION X(*)
*     ..
*
*  Purpose
*  =======
*
*  DNRM2 returns the euclidean norm of a vector via the function
*  name, so that
*
*     DNRM2 := sqrt( x'*x )
*
*
*  -- This version written on 25-October-1982.
*     Modified on 14-October-1993 to inline the call to DLASSQ.
*     Sven Hammarling, Nag Ltd.
*
*
*     .. Parameters ..
CID      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
*     ..
*     .. Local Scalars ..
CID      DOUBLE PRECISION ABSXI,NORM,SCALE,SSQ
      REAL NORM
      INTEGER IX
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC ABS,SQRT
*     ..
      IF (N.LT.1 .OR. INCX.LT.1) THEN
          NORM = ZERO
      ELSE IF (N.EQ.1) THEN
          NORM = ABS(X(1))
      ELSE
          SCALE = ZERO
          SSQ = ONE
*        The following loop is equivalent to this call to the LAPACK
*        auxiliary routine:
*        CALL DLASSQ( N, X, INCX, SCALE, SSQ )
*
          DO 10 IX = 1,1 + (N-1)*INCX,INCX
              IF (X(IX).NE.ZERO) THEN
                  ABSXI = ABS(X(IX))
                  IF (SCALE.LT.ABSXI) THEN
                      SSQ = ONE + SSQ* (SCALE/ABSXI)**2
                      SCALE = ABSXI
                  ELSE
                      SSQ = SSQ + (ABSXI/SCALE)**2
                  END IF
              END IF
   10     CONTINUE
          NORM = SCALE*SQRT(SSQ)
      END IF
*
      DNRM2 = NORM
      RETURN
*
*     End of DNRM2.
*
      END
C
C=======================================================================
C=======================================================================
C                    RANDOM NUMBER GENERATION PACKAGE
C=======================================================================
C=======================================================================
C
C     Note for human users: how to use these subroutines
C     ====================
C
C     First we initialize the "seed" by calling RANDIN(I_SEED), where
C     I_SEED can be set to any value 
C
C     Next generate your random numbers (uniformly distributed within
C     0 and 1) by calling ZUFALL(N,A);  N is the actual number of the
C     random numbers needed and A is the vector that contains them 
C
C     If you want to restore the seed to its original state, you must
C     call ZUFALLSV,   b e f o r e  the last call of ZUFALL  (see the 
C     example below)
C
C     Typical use
C     ===========
C
C     PARAMETER
C    *         (N=500)
C     DIMENSION
C    *          A(1:N),B(1:N),SVBLOK(608)
C
C
C     CALL RANDIN(0)          ! initialize the seed
C 
C     CALL ZUFALLSV(SVBLOK)   ! save the original seed for future
C                             ! restart
C
C     CALL ZUFALL(500,A)      ! generates the random numbers in the
C                             ! vector A
C
C      ....
C
C     CALL ZUFALLRS(SVBLOK)   ! restores the original seed as it was
C                             ! before
C
C     CALL ZUFALL(500,B)      ! B should be identical to A
C
C=======================================================================
C=======================================================================
C
C                 README FOR ZUFALL RANDOM NUMBER PACKAGE
C                 ------ --- ------ ------ ------ -------
C
C     This package contains a portable random number generator set for: 
C
C          a. uniform (u in [0,1)), 
C          b. normal (<g> = 0, <g^2> = 1), and
C          c. Poisson distributions. 
C
C     The basic module, the uniform generator, uses a lagged Fibonacci 
C     series generator:
C     
C                   t    = u(n-273) + u(n-607)
C                   u(n) = t - float(int(t))
C     
C     where each number generated, u(k), is floating point.  Since the 
C     numbers are floating point,  the left-end boundary  of the range 
C     contains zero.  This package  is nearly portable, except for the 
C     following:
C     
C         (1) It is written in lower case, 
C         (2) the test package contains a timer (second), which is not 
C             portable, and 
C         (3) there are cycle times (in seconds) in data statements 
C             for NEC SX-3, Fujitsu VP2200, and Cray Y-MP. 
C
C     Select  your  favorite  and comment-out the others.  Replacement 
C     functions for 'second' are included - comment out the others. 
C     Otherwise the package is portable and returns the same set of 
C     floating point numbers up to word precision on any machine. 
C     There are compiler directives ($cdir for Cray, Cvdir for SX-3, 
C     and VOCL for Fujitsu VP2200) which should be otherwise ignored.
C 
C     To compile this beast, note that all floating point numbers
C     are declared 'double precision'. On Cray X-MP, Y-MP, and C-90
C     machines, use the cft77 (cf77) option -dp to run this in 64
C     bit mode (not 128 bit double).
C 
C     External documentation, 
C
C    "Lagged Fibonacci Random Number Generators for the NEC SX-3", 
C     is to be published in the International Journal of High Speed 
C     Computing (1994). Otherwise, ask the author: 
C 
C          W. P. Petersen 
C          IPS, RZ F-5
C          ETHZ
C          CH 8092, Zurich
C          Switzerland
C 
C          e-mail:  wpp@ips.ethz.ch.
C 
C======================================================================= 
C 
C     THE PACKAGE CONTAINS THE FOLLOWING ROUTINES:
C 
C======================================================================= 
C    
C          UNIFORM generator routines:
C          =======
C
C          subroutine RANDIN(seed)
C          integer seed
C    
C                  initializes common block containing seeds. if seed=0,
C                  the default value is 1802.
C    
C          subroutine zufall(n,u)
C          integer n
C          double precision u(n)
C    
C                  returns set of n uniforms u(1), ..., u(n).
C    
C          subroutine zufallsv(zusave)
C          double precision zusave(608)
C    
C                  saves buffer and pointer in zusave, for later 
C                  restarts
C    
C          subroutine zufallrs(zusave)
C          double precision zusave(608)
C    
C                  restores seed buffer and pointer from zusave
C    
C=======================================================================
C    
C     NORMAL generator routines:
C     ======
C
C          subroutine normalen(n,g)
C          integer n
C          double precision g(n)
C    
C             returns set of n normals g(1), ..., g(n) such that
C             mean <g> = 0, and variance <g**2> = 1.
C    
C          subroutine normalsv(normsv)
C          double precision normsv(1634)
C    
C             saves zufall seed buffer and pointer in normsv
C             buffer/pointer for normalen restart also in normsv
C    
C          subroutine normalrs(normsv)
C          double precision normsv(1634)
C    
C             restores zufall seed buffer/pointer and 
C             buffer/pointer for normalen restart from normsv
C    
C=======================================================================
C    
C    POISSON generator routine:
C    =======
C
C          subroutine fische(n,mu,q)
C          integer n,q(n)
C          double precision mu
C    
C              returns set of n integers q, with poisson
C              distribution, density p(q,mu) = exp(-mu) mu**q/q!
C      
C          USE zufallsv and zufallrs for stop/restart sequence
c
c          this program tests three random number generators:
c------------------------------------------------------------
c          IMPORTANT!!     be sure to compile with -dp option on
c          the Cray cft77 compiler. This disables 128-bit double
c          precision, converting it to 64-bit single.  Likewise,
c          comment out the 'second' function. On other machines,
c          uncomment lines for appropriate version of 'second'.
c------------------------------------------------------------
c      zufallt = a uniform distribution of r.v.'s test
c
c
c=====================================================================
c=====================================================================
c
      subroutine zufall(n,a)
      
C     implicit none

c
c portable lagged Fibonacci series uniform random number
c generator with "lags" -273 und -607:
c
c       t    = u(i-273)+buff(i-607)  (floating pt.)
c       u(i) = t - float(int(t))
c
c W.P. Petersen, IPS, ETH Zuerich, 19 Mar. 92
c
C     double precision a(*)
      DIMENSION a(*)
      
C     double precision buff(607)
      DIMENSION buff(607)
      
C     double precision t

      integer i,k,ptr,VL,k273,k607
      
      integer buffsz,nn,n,left,q,qq
      
      integer aptr,aptr0,bptr
c
      common /klotz0/buff,ptr
      data buffsz/607/
c
      aptr = 0
      nn   = n
c
1     continue
c
      if(nn .le. 0) return
c
c factor nn = q*607 + r
c
      q    = (nn-1)/607
      left = buffsz - ptr
c
      if(q .le. 1) then
c
c only one or fewer full segments
c
         if(nn .lt. left) then
            do 2 i=1,nn
               a(i+aptr) = buff(ptr+i)
2           continue
            ptr  = ptr + nn
            return
         else
            do 3 i=1,left
               a(i+aptr) = buff(ptr+i)
3           continue
            ptr  = 0
            aptr = aptr + left
            nn   = nn - left
c  buff -> buff case
            VL   = 273
            k273 = 334
            k607 = 0
            do 4 k=1,3
cdir$ ivdep
*vdir nodep
*VOCL LOOP, TEMP(t), NOVREC(buff)
               do 5 i=1,VL
                  t            = buff(k273+i) + buff(k607+i)
                  buff(k607+i) = t - float(int(t))
5              continue
               k607 = k607 + VL
               k273 = k273 + VL
               VL   = 167
               if(k.eq.1) k273 = 0
4           continue
c
            goto 1
         endif
      else
c
c more than 1 full segment
c 
          do 6 i=1,left
             a(i+aptr) = buff(ptr+i)
6         continue
          nn   = nn - left
          ptr  = 0
          aptr = aptr+left
c 
c buff -> a(aptr0)
c 
          VL   = 273
          k273 = 334
          k607 = 0
          do 7 k=1,3
             if(k.eq.1)then
*VOCL LOOP, TEMP(t)
                do 8 i=1,VL
                   t         = buff(k273+i) + buff(k607+i)
                   a(aptr+i) = t - float(int(t))
8               continue
                k273 = aptr
                k607 = k607 + VL
                aptr = aptr + VL
                VL   = 167
             else
cdir$ ivdep
*vdir nodep
*VOCL LOOP, TEMP(t)
                do 9 i=1,VL
                   t         = a(k273+i) + buff(k607+i)
                   a(aptr+i) = t - float(int(t))
9               continue
                k607 = k607 + VL
                k273 = k273 + VL
                aptr = aptr + VL
             endif
7         continue
          nn = nn - 607
c
c a(aptr-607) -> a(aptr) for last of the q-1 segments
c
          aptr0 = aptr - 607
          VL    = 607
c
*vdir novector
          do 10 qq=1,q-2
             k273 = 334 + aptr0
cdir$ ivdep
*vdir nodep
*VOCL LOOP, TEMP(t), NOVREC(a)
             do 11 i=1,VL
                t         = a(k273+i) + a(aptr0+i)
                a(aptr+i) = t - float(int(t))
11           continue
             nn    = nn - 607
             aptr  = aptr + VL
             aptr0 = aptr0 + VL
10        continue
c
c a(aptr0) -> buff, last segment before residual
c
          VL   = 273
          k273 = 334 + aptr0
          k607 = aptr0
          bptr = 0
          do 12 k=1,3
             if(k.eq.1) then
*VOCL LOOP, TEMP(t)
                do 13 i=1,VL
                   t            = a(k273+i) + a(k607+i)
                   buff(bptr+i) = t - float(int(t))
13              continue
                k273 = 0
                k607 = k607 + VL
                bptr = bptr + VL
                VL   = 167
             else
cdir$ ivdep
*vdir nodep
*VOCL LOOP, TEMP(t), NOVREC(buff)
                do 14 i=1,VL
                   t            = buff(k273+i) + a(k607+i)
                   buff(bptr+i) = t - float(int(t))
14              continue
                k607 = k607 + VL
                k273 = k273 + VL
                bptr = bptr + VL
             endif
12        continue
          goto 1
      endif
      end
c
c=====================================================================
c=====================================================================
c
      subroutine RANDIN(seed)
      
C     implicit none
c
c  generates initial seed buffer by linear congruential
c  method. Taken from Marsaglia, FSU report FSU-SCRI-87-50
c  variable seed should be 0 < seed <31328
c
      integer seed
      integer ptr
      
C     double precision s,t
      
C     double precision buff(607)
      DIMENSION buff(607)
      
      integer ij,kl,i,ii,j,jj,k,l,m
      common /klotz0/buff,ptr
      data ij/1802/,kl/9373/
c
c=====================================================================
c
      if(seed.ne.0) ij = seed
c
      i = mod(ij/177,177) + 2
      j = mod(ij,177) + 2
      k = mod(kl/169,178) + 1
      l = mod(kl,169)
      do 1 ii=1,607
         s = 0.0
         t = 0.5
         do 2 jj=1,24
            m = mod(mod(i*j,179)*k,179)
            i = j
            j = k
            k = m
            l = mod(53*l+1,169)
            if(mod(l*m,64).ge.32) s = s+t
            t = .5*t
2        continue
         buff(ii) = s
1     continue
c
c=====================================================================
c
      return
      end
c
c=====================================================================
c=====================================================================
c
      subroutine zufallsv(svblk)
      
C     implicit none
c
c  saves common blocks klotz0, containing seeds and 
c  pointer to position in seed block. IMPORTANT: svblk must be
c  dimensioned at least 608 in driver. The entire contents
c  of klotz0 (pointer in buff, and buff) must be saved.
c
C     double precision buff(607)
      DIMENSION buff(607)
      
      integer ptr,i
      
C     double precision svblk(*)
      DIMENSION svblk(*)
c
      common /klotz0/ buff,ptr
c
c=====================================================================
c
      svblk(1) = ptr
      do 1 i=1,607
         svblk(i+1) = buff(i)
1     continue
c
c=====================================================================
c
      return
      end
c
c=====================================================================
c=====================================================================
c
      subroutine zufallrs(svblk)
      
C     implicit none
c
c  restores common block klotz0, containing seeds and pointer
c  to position in seed block. IMPORTANT: svblk must be
c  dimensioned at least 608 in driver. The entire contents
c  of klotz0 must be restored.
c
C     double precision buff(607)
      DIMENSION buff(607)
      
      integer i,ptr
      
C     double precision svblk(*)
      DIMENSION svblk(*)
c
      common /klotz0/ buff,ptr
c
      ptr = svblk(1)
      do 1 i=1,607
         buff(i) = svblk(i+1)
1     continue
c
      return
      end
c
c=====================================================================
c=====================================================================
c
*******************************************************************
********	FILE: randgen.f				***********
********	AUTHORS: Richard Chandler		***********
********		 (richard@stats.ucl.ac.uk)	***********
********		 Paul Northrop 			***********
********		 (northrop@stats.ox.ac.uk)	***********
********	LAST MODIFIED: 26/8/03			***********
********	See file randgen.txt for details	***********
*******************************************************************
*
      BLOCK DATA ZBQLBD01
*
*       Initializes seed array etc. for random number generator.
*       The values below have themselves been generated using the
*       NAG generator.
*
      COMMON /ZBQL0001/ ZBQLIX,B,C
      DOUBLE PRECISION ZBQLIX(43),B,C
      INTEGER I
      DATA (ZBQLIX(I),I=1,43) /8.001441D7,5.5321801D8,
     +1.69570999D8,2.88589940D8,2.91581871D8,1.03842493D8,
     +7.9952507D7,3.81202335D8,3.11575334D8,4.02878631D8,
     +2.49757109D8,1.15192595D8,2.10629619D8,3.99952890D8,
     +4.12280521D8,1.33873288D8,7.1345525D7,2.23467704D8,
     +2.82934796D8,9.9756750D7,1.68564303D8,2.86817366D8,
     +1.14310713D8,3.47045253D8,9.3762426D7 ,1.09670477D8,
     +3.20029657D8,3.26369301D8,9.441177D6,3.53244738D8,
     +2.44771580D8,1.59804337D8,2.07319904D8,3.37342907D8,
     +3.75423178D8,7.0893571D7 ,4.26059785D8,3.95854390D8,
     +2.0081010D7,5.9250059D7,1.62176640D8,3.20429173D8,
     +2.63576576D8/
      DATA B / 4.294967291D9 /
      DATA C / 0.0D0 /
      END
******************************************************************
******************************************************************
******************************************************************
      SUBROUTINE ZBQLINI(SEED)
******************************************************************
*       To initialize the random number generator - either
*       repeatably or nonrepeatably. Need double precision
*       variables because integer storage can't handle the
*       numbers involved
******************************************************************
*	ARGUMENTS
*	=========
*	SEED	(integer, input). User-input number which generates
*		elements of the array ZBQLIX, which is subsequently used 
*		in the random number generation algorithm. If SEED=0,
*		the array is seeded using the system clock if the 
*		FORTRAN implementation allows it.
******************************************************************
*	PARAMETERS
*	==========
*	LFLNO	(integer). Number of lowest file handle to try when
*		opening a temporary file to copy the system clock into.
*		Default is 80 to keep out of the way of any existing
*		open files (although the program keeps searching till
*		it finds an available handle). If this causes problems,
*               (which will only happen if handles 80 through 99 are 
*               already in use), decrease the default value.
******************************************************************
      INTEGER LFLNO
      PARAMETER (LFLNO=80)
******************************************************************
*	VARIABLES
*	=========
*	SEED	See above
*	ZBQLIX	Seed array for the random number generator. Defined
*		in ZBQLBD01
*	B,C	Used in congruential initialisation of ZBQLIX
*	SS,MM,}	System clock secs, mins, hours and days
*	HH,DD }
*	FILNO	File handle used for temporary file
*	INIT	Indicates whether generator has already been initialised
*
      INTEGER SEED,SS,MM,HH,DD,FILNO,I
      INTEGER INIT
      DOUBLE PRECISION ZBQLIX(43),B,C
      DOUBLE PRECISION TMPVAR1,DSS,DMM,DHH,DDD

      COMMON /ZBQL0001/ ZBQLIX,B,C
      COMMON
     *       /JDINFO/ NO_DAY,N_HOUR,MINUTE,NSECND
      SAVE INIT

*
*	Ensure we don't call this more than once in a program
*
      IF (INIT.GE.1) THEN
       IF(INIT.EQ.1) THEN
        WRITE(*,1)
        INIT = 2
       ENDIF
       RETURN
      ELSE
       INIT = 1
      ENDIF
*
*       If SEED = 0, cat the contents of the clock into a file
*       and transform to obtain ZQBLIX(1), then use a congr.
*       algorithm to set remaining elements. Otherwise take
*       specified value of SEED.
*
*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*>>>>>>>	NB FOR SYSTEMS WHICH DO NOT SUPPORT THE  >>>>>>>
*>>>>>>>	(NON-STANDARD) 'CALL SYSTEM' COMMAND,    >>>>>>>
*>>>>>>>	THIS WILL NOT WORK, AND THE FIRST CLAUSE >>>>>>>
*>>>>>>>	OF THE FOLLOWING IF BLOCK SHOULD BE	 >>>>>>>
*>>>>>>>	COMMENTED OUT.				 >>>>>>>
*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      IF (SEED.EQ.0) THEN
*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*>>>>>>>	COMMENT OUT FROM HERE IF YOU DON'T HAVE  >>>>>>>
*>>>>>>>	'CALL SYSTEM' CAPABILITY ...		 >>>>>>>
*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
       CALL SYSTEM(' date +%S%M%H%j > zbql1234.tmp')
*
*       Try all file numbers for LFLNO to 999 
*
       FILNO = LFLNO
 10    OPEN(FILNO,FILE='zbql1234.tmp',ERR=11)
       GOTO 12
 11    FILNO = FILNO + 1
       IF (FILNO.GT.999) THEN
        WRITE(*,2)
        RETURN
       ENDIF
       GOTO 10
 12    READ(FILNO,'(3(I2),I3)') SS,MM,HH,DD
C
       NO_DAY=DD
       N_HOUR=HH
       MINUTE=MM
       NSECND=SS
C       
       CLOSE(FILNO)
       CALL SYSTEM('rm zbql1234.tmp')
       DSS = DINT((DBLE(SS)/6.0D1) * B)
       DMM = DINT((DBLE(MM)/6.0D1) * B)
       DHH = DINT((DBLE(HH)/2.4D1) * B)
       DDD = DINT((DBLE(DD)/3.65D2) * B)
       TMPVAR1 = DMOD(DSS+DMM+DHH+DDD,B)
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
*<<<<<<<<	... TO HERE (END OF COMMENTING OUT FOR 	  <<<<<<<
*<<<<<<<<	USERS WITHOUT 'CALL SYSTEM' CAPABILITY	  <<<<<<<
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      ELSE
       TMPVAR1 = DMOD(DBLE(SEED),B)
      ENDIF
      ZBQLIX(1) = TMPVAR1
      DO 100 I = 2,43
       TMPVAR1 = ZBQLIX(I-1)*3.0269D4
       TMPVAR1 = DMOD(TMPVAR1,B)       
       ZBQLIX(I) = TMPVAR1
 100  CONTINUE

 1    FORMAT(//5X,'****WARNING**** You have called routine ZBQLINI ',
     +'more than',/5X,'once. I''m ignoring any subsequent calls.',//)
 2    FORMAT(//5X,'**** ERROR **** In routine ZBQLINI, I couldn''t',
     +' find an',/5X,
     +'available file number. To rectify the problem, decrease the ',
     +'value of',/5X,
     +'the parameter LFLNO at the start of this routine (in file ',
     +'randgen.f)',/5X,
     +'and recompile. Any number less than 100 should work.')
      END
******************************************************************
      FUNCTION ZBQLU01(DUMMY)
*
*       Returns a uniform random number between 0 & 1, using
*       a Marsaglia-Zaman type subtract-with-borrow generator.
*       Uses double precision, rather than integer, arithmetic 
*       throughout because MZ's integer constants overflow
*       32-bit integer storage (which goes from -2^31 to 2^31).
*       Ideally, we would explicitly truncate all integer 
*       quantities at each stage to ensure that the double
*       precision representations do not accumulate approximation
*       error; however, on some machines the use of DNINT to
*       accomplish this is *seriously* slow (run-time increased
*       by a factor of about 3). This double precision version 
*       has been tested against an integer implementation that
*       uses long integers (non-standard and, again, slow) -
*       the output was identical up to the 16th decimal place
*       after 10^10 calls, so we're probably OK ...
*
      DOUBLE PRECISION ZBQLU01,DUMMY,B,C,ZBQLIX(43),X,B2,BINV
      INTEGER CURPOS,ID22,ID43

      COMMON /ZBQL0001/ ZBQLIX,B,C
      SAVE /ZBQL0001/
      SAVE CURPOS,ID22,ID43
      DATA CURPOS,ID22,ID43 /1,22,43/

      B2 = B
      BINV = 1.0D0/B
 5    X = ZBQLIX(ID22) - ZBQLIX(ID43) - C
      IF (X.LT.0.0D0) THEN
       X = X + B
       C = 1.0D0
      ELSE
       C = 0.0D0
      ENDIF
      ZBQLIX(ID43) = X
*
*     Update array pointers. Do explicit check for bounds of each to
*     avoid expense of modular arithmetic. If one of them is 0 the others
*     won't be
*
      CURPOS = CURPOS - 1
      ID22 = ID22 - 1
      ID43 = ID43 - 1
      IF (CURPOS.EQ.0) THEN
       CURPOS=43
      ELSEIF (ID22.EQ.0) THEN
       ID22 = 43
      ELSEIF (ID43.EQ.0) THEN
       ID43 = 43
      ENDIF
*
*     The integer arithmetic there can yield X=0, which can cause 
*     problems in subsequent routines (e.g. ZBQLEXP). The problem
*     is simply that X is discrete whereas U is supposed to 
*     be continuous - hence if X is 0, go back and generate another
*     X and return X/B^2 (etc.), which will be uniform on (0,1/B). 
*
      IF (X.LT.BINV) THEN
       B2 = B2*B
       GOTO 5
      ENDIF

      ZBQLU01 = X/B2

      END
******************************************************************
******************************************************************
******************************************************************
      FUNCTION ZBQLNOR(MU,SIGMA)
*
*       Returns a random number Normally distributed with mean
*       MU and standard deviation |SIGMA|, using the Box-Muller
*       algorithm
*
      DOUBLE PRECISION THETA,R,ZBQLNOR,ZBQLU01,PI,MU,SIGMA
      DOUBLE PRECISION SPARE
      INTEGER STATUS
      SAVE STATUS,SPARE,PI
      DATA STATUS /-1/
     
      IF (STATUS.EQ.-1) PI = 4.0D0*DATAN(1.0D0)

      IF (STATUS.LE.0) THEN
       THETA = 2.0D0*PI*ZBQLU01(0.0D0)
       R = DSQRT( -2.0D0*DLOG(ZBQLU01(0.0D0)) )
       ZBQLNOR = (R*DCOS(THETA))
       SPARE = (R*DSIN(THETA))
       STATUS = 1
      ELSE
       ZBQLNOR = SPARE
       STATUS = 0
      ENDIF
      
      ZBQLNOR = MU + (SIGMA*ZBQLNOR)

      END
C
C=======================================================================
C=======================================================================
C

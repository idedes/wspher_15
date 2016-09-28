C FILE NAME = wspher08_diamat_15.f ! Keep this symbol:    $ident@string$
C       
C=======================================================================
C=======================================================================
C                SYMMETRIC MATRIX DIAGONALISATION PACKAGE     
C=======================================================================
C=======================================================================
C
      SUBROUTINE DIAMAT(HAMILT,SPENER,AUXVEC,NDBASE,LDBASE,EIGVEC)
C
      LOGICAL
     *          ISFULL,EIGVEC
C
      DIMENSION
     *          HAMILT(1:NDBASE,1:NDBASE)
      DIMENSION
     *          SPENER(1:NDBASE),
     *          AUXVEC(1:NDBASE)
      COMMON
     *       /CNTCNT/ ICOUNT_DENGRD,ICOUNT_FUNMIN,ICOUNT_HAMMAT,
     *                ICOUNT_EXPTHE,ICOUNT_CHOICE,ICOUNT_DIAMAT
C
C=======================================================================
C
C     This routine diagonalises a real symmetric LDBASExLDBASE matrix 
C
C     EIGVEC => To be set .FALSE.  if  wave functions  are not needed
C     AUXVEC => Auxiliary vector,  its starting values are irrelevant
C     SPENER => At output the eigen-values, on the entry - irrelevant
C
C     HAMILT => The matrix to diagonalize -> must be symmetric though
C
C     NDBASE => The maximum dimension of matrices to be diagonalised
C     LDBASE => The actual  dimension  of the basis  for  our problem
C
C=======================================================================
C
      CALL CPUTIM('DIAMAT',1)
      ICOUNT_DIAMAT=ICOUNT_DIAMAT+1
C      
      ISFULL=.TRUE.
C
      CALL DIAGON(HAMILT,LDBASE,NDBASE,SPENER,AUXVEC,EIGVEC,ISFULL)
C      
      CALL CPUTIM('DIAMAT',0)
C
C=======================================================================
C
      RETURN
      END
C
C=======================================================================
C=======================================================================
C
      SUBROUTINE DIAGON(A,N,NP,D,E,EIGVEC,ISFULL)
C
C     This routine diagonalises a REAL and SYMMETRIC matrix; based
C     on "Numerical Recipies"
C
C_______________________________________________________________________
C_______________________________________________________________________
C
C                 TO BE COMPILED IN AUTODOUBLE MODE !!!!
C                 TO BE COMPILED IN AUTODOUBLE MODE !!!!
C                 TO BE COMPILED IN AUTODOUBLE MODE !!!!
C_______________________________________________________________________
C_______________________________________________________________________
C
C
C        Parameters of the routine:
C
C           A  - NxN matrix, declared as NPxNP in the caller
C           N  - integer: actual dimension of A
C           NP - integer: physical dimension of A declared in the caller
C           D  - vector of dimension at least N
C           E  - vector of dimension at least N
C           EIGVEC  - logical
C           ISFULL  - logical
C
C
C                The routine  can perform  two tasks:
C                ------------------------------------
C
C    1.) ISFULL = .TRUE. (i.e., matrix A is full,  n o t  tridiagonal)
C
C        Finds  eigenvalues  and (optionally)  eigenvectors of a real,
C        symmetric, NxN matrix  declared  as NPxNP  matrix in calling
C        program.  D on output will  contain eigenvalues in ascending
C        order. On input the contents od D and E is irrelevant.
C        Then:
C
C           If EIGVEC is .FALSE., eigenvectors will not be calculated
C           (but matrix A will be changed anyway).
C
C           If EIGVEC is  .TRUE., the normalised  eigenvector corres-
C           ponding to eigenvalue D(K) will sit in A(I,K), I = 1, N.
C
C           Vector E provides working space; it is irrelevant on both
C           input and output.
C
C        In this case the routine will call "TRES2"
C
C    2.) ISFULL = .FALSE. (i.e., matrix A is already tridiagonal)
C
C        Finds  eigenvalues  and (optionally)  eigenvectors of a real,
C        symmetric, tridiagonal, NxN matrix declared as NPxNP  matrix
C        in calling program. D on output will  contain eigenvalues in
C        ascending order.  On input D  contains  diagonal elements of
C        the matrix  to be diagonalised,  while vector E(I), I = 2, N
C        contains its subdiagonal elements - E(1) is irrelevant.
C        Values in matrix A are also irrelevant on input in this case.
C        On output D and A will have the same meaning as in case 1.
C        The vector E will be changed, but its contents is irrelevant
C        on output.
C
C_______________________________________________________________________
C
C        THE ROUTINE USES "PYTHAG" AND "TRES2" (INCLUDED BELOW)
C_______________________________________________________________________
C
C        BASED ON  "NUMERICAL  RECIPES"
C                                                  T.R. Werner, May 1991
C                                           Last modified: November 1993
C_______________________________________________________________________
C
      PARAMETER 
     *         (NDITER=100,ONEONE=1.0)
C
      DIMENSION  
     *          A(NP,N),D(N),E(N)
C
      LOGICAL    
     *          NONLY,EIGVEC,ISFULL
C
C=======================================================================
C
C      CALL CPUTIM('DIAGON',1)
C
C=======================================================================
C
      IF (N.GT.NP.OR.N.LT.1) STOP 'Incorrrect matrix size in DIAGON'
C
C=======================================================================
C
      IF (N.EQ.1) THEN
          D(1)=A(1,1)
          A(1,1)=1
          RETURN
      END IF
C
C=======================================================================
C
      IF (EIGVEC) THEN
C
          NONLY=.TRUE.
C
          IF (.NOT.ISFULL) THEN
              DO J=1,N
                 DO I=1,N
                    A(I,J)=0
                 END DO
                 A(J,J)=1
              END DO
          ELSE
              CALL TRES2(A,N,NP,D,E,NONLY)
          END IF
C
      ELSE
                                           NONLY=.FALSE.
         IF (ISFULL) CALL TRES2(A,N,NP,D,E,NONLY)
C
      END IF
C
      DO I=2,N
         E(I-1)=E(I)
      END DO
C
      E(N)=0
C
      DO 4 L=1,N
         ITER=0
    1    DO M=L,N-1
            DD=ABS(D(M))+ABS(D(M+1))
            IF (ABS(E(M))+DD.EQ.DD) GO TO 2
         END DO
         M=N
C
    2    IF (M.NE.L) THEN
C
            IF (ITER.EQ.NDITER) THEN
               PRINT '(''No convergence in DIAGON. ITER='',I3)',
     *                                                    ITER
               STOP 'No convergence in DIAGON'
            END IF
C
            ITER=ITER+1
            G=(D(L+1)-D(L))/(2*E(L))
            R=PYTHAG(G,ONEONE)
            G=D(M)-D(L)+E(L)/(G+SIGN(R,G))
            S=1
            C=1
            P=0
C
            DO 3 I=M-1,L,-1
               F=S*E(I)
               B=C*E(I)
               R=PYTHAG(F,G)
               E(I+1)=R
               IF (R.EQ.0) THEN
                  D(I+1)=D(I+1)-P
                  E(M)=0
                  GO TO 1
               END IF
               S=F/R
               C=G/R
               G=D(I+1)-P
               R=(D(I)-G)*S+2*C*B
               P=S*R
               D(I+1)=G+P
               G=C*R-B
               IF (NONLY) THEN
                  DO  K=1,N
                     F=A(K,I+1)
                     A(K,I+1)=S*A(K,I)+C*F
                     A(K,I)=C*A(K,I)-S*F
                  END DO
               END IF
    3       END DO
C
            D(L)=D(L)-P
            E(L)=G
            E(M)=0
            GO TO 1
         ENDIF
    4 END DO
C
      DO 5 I=1,N-1
         K=I
         P=D(I)
         DO J=I+1,N
            IF (D(J).LT.P) THEN
               K=J
               P=D(J)
            END IF
         END DO
         IF (K.NE.I) THEN
            D(K)=D(I)
            D(I)=P
            DO J=1,N
               P=A(J,I)
               A(J,I)=A(J,K)
               A(J,K)=P
            END DO
         END IF
    5 END DO
C
C=======================================================================
C
C      CALL CPUTIM('DIAGON',0)
C
C=======================================================================
C      
      RETURN
      END
C
C=======================================================================
C
      SUBROUTINE TRES2(A,N,NP,D,E,NONLY)
C
      DIMENSION 
     *          A(NP,N),D(N),E(N)
      LOGICAL   
     *          NONLY
C
      IF (N.GT.NP.OR.N.LT.1) THEN
          STOP 'Incorrect matrix dimension in TRES2'
      ELSE IF (N.EQ.1) THEN
          D(1)=A(1,1)
          E(1)=0
          A(1,1)=1
          RETURN
      END IF
C
      DO 1 I=N,2,-1
         L=I-1
         H=0
         SCALE=0
         IF (L.GT.1) THEN
            DO K=1,L
               SCALE=SCALE+ABS(A(I,K))
            END DO
            IF (SCALE.EQ.0) THEN
               E(I)=A(I,L)
            ELSE
               DO K=1,L
                  A(I,K)=A(I,K)/SCALE
                  H=H+A(I,K)**2
               END DO
               F=A(I,L)
               G=-SIGN(SQRT(H),F)
               E(I)=SCALE*G
               H=H-F*G
               A(I,L)=F-G
               F=0
               DO J=1,L
                  IF (NONLY) A(J,I)=A(I,J)/H
                  G=0
                  DO K=1,J
                     G=G+A(J,K)*A(I,K)
                  END DO
                  IF (L.GT.J) THEN
                     DO K=J+1,L
                        G=G+A(K,J)*A(I,K)
                     END DO
                  ENDIF
                  E(J)=G/H
                  F=F+E(J)*A(I,J)
               END DO
               HH=F/(H+H)
               DO J=1,L
                  F=A(I,J)
                  G=E(J)-HH*F
                  E(J)=G
                  DO K=1,J
                     A(J,K)=A(J,K)-F*E(K)-G*A(I,K)
                  END DO
               END DO
            ENDIF
         ELSE
            E(I)=A(I,L)
         ENDIF
         D(I)=H
    1 END DO
C
      E(1)=0
C
      IF (.NOT.NONLY) THEN
         DO I=1,N
            D(I)=A(I,I)
         END DO
         RETURN
      ELSE
C
         D(1)=0
C
         DO 2 I=1,N
            L=I-1
            IF (D(I).NE.0) THEN
               DO J=1,L
                  G=0
                  DO K=1,L
                     G=G+A(I,K)*A(K,J)
                  END DO
                  DO K=1,L
                     A(K,J)=A(K,J)-G*A(K,I)
                  END DO
               END DO
            ENDIF
            D(I)=A(I,I)
            A(I,I)=1
            IF (L.GE.1) THEN
               DO J=1,L
                  A(I,J)=0
                  A(J,I)=0
               END DO
            END IF
    2    END DO
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
      FUNCTION PYTHAG(A,B)
C
C     THE ROUTINE CALCULATES SQRT(A**2+B**2) ACCURATELY
C
      ABSA=ABS(A)
      ABSB=ABS(B)
C
      IF (ABSA.GT.ABSB) THEN
          PYTHAG=ABSA*SQRT(1+(ABSB/ABSA)**2)
      ELSE
         IF (ABSB.EQ.0) THEN
             PYTHAG=0
         ELSE
             PYTHAG=ABSB*SQRT(1+(ABSA/ABSB)**2)
         END IF
      END IF
C
      END
C      
C=======================================================================
C=======================================================================
C

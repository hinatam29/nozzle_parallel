c ======================================================================
c
c         Solving Elliptic-Type Partial Diffrential Equation
c
c                with 7 ponits Finite Difference Method
c
c                    to obtain Electrical Potential
c
c    Incomplete LU decomposition Bi Conjugate Gradient Stable method
c
c     programer :  H. Kobayashi
c
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     (i-1,j,k)   = ijk - nt*nz
c     (i,j-1,k)   = ijk - nz
c     (i,j,k-1)   = ijk - 1
c     (i,j,k)     = ijk
c     (i,j,k+1)   = ijk + 1
c     (i,j+1,k)   = ijk + nz
c     (i+1,j,k)   = ijk + nt*nz
c ----------------------------------------------------------------------
      subroutine pfield_bicg
      include 'subcom.inc'

      parameter( n3 = (nx-1)*(ny-1) )

      parameter( M3=ny-1, M6=1 )
      parameter( N1=n3+M3, N2=-M3 )
      dimension F1(n3), F2(n3), F3(n3)
      dimension A(N2:N1,9)
      dimension B(n3), R0(n3)
      dimension DD(N2:N1), PP(N2:n3+M3), RR(N2:n3+M3)
      dimension S(N2:n3+M3), T(N2:n3+M3), V(N2:n3+M3)
      dimension XX(N2:n3+M3)

      dimension convr(nx,ny), convz(nx,ny)

      ITR = 10000
      EPS = 1.0E-40
      SGM = 0.0

      pfix = 1.013e5

c
c: --- initialize ---
c
      do 100 i=1, n3 
         F1(i)  = 0.0
         F2(i)  = 0.0
         F3(i)  = 0.0
         B(i) = 0.0
  100 continue


c#COLLAPSE
      do 110 i=1, 9
         do 110 j=N2, N1
            A(j,i)= 0.0
  110 continue

      do 115 i=1, nx
      do 115 j=1, ny
        p0(i,j) = 1.013e5
115   continue


c
c: ----- set coefficient Fk and Pe -----
c

      do 118 i=1, nx-2
      do 118 k=1, ny-2
        ijk = jj(i,k)

        convr(i,k) = ux0(i,k)*(ux_sx(i+1,k) - ux_sx(i,k))
     &              / (x(i+1,k) - x(i,k))
     &             + uy0(i,k)*(ux_sy(i,k+1)-ux_sy(i,k))
     &              / (y(i,k+1) - y(i,k))

        convz(i,k) = ux0(i,k)*(uy_sx(i+1,k) - uy_sx(i,k))
     &              / (x(i+1,k) - x(i,k))
     &             + uy0(i,k)*(uy_sy(i,k+1) - uy_sy(i,k))
     &              / (y(i,k+1) - y(i,k))
 118  continue

      do 120 i=2, nx-2
      do 120 k=2, ny-2
        ijk = jj(i,k)

        xx1 = 0.5*(x(i,k) + x(i-1,k))
        xx2 = 0.5*(x(i+1,k) + x(i,k))
        xx3 = 0.5*(x(i+2,k) + x(i+1,k))
        yy1 = 0.5*(y(i,k) + y(i,k-1))
        yy2 = 0.5*(y(i,k+1) + y(i,k))
        yy3 = 0.5*(y(i,k+2) + y(i,k+1))

        divu = (-ux_sx(i,k)+ux_sx(i+1,k))/(x(i+1,k)-x(i,k))
     &         +   ux0(i,k)/xx2
     &         +  (-uy_sy(i,k)+uy_sy(i,k+1))/(y(i,k+1)-y(i,k))

c       B(ijk) = -( (convr(i+1,k) - convr(i-1,k))/(xx3 - xx1)
c    &           +  convr(i,k)/xx2
c    &            + (convz(i,k+1) - convz(i,k-1))/(yy3 - yy1) )
c    &           +  divu/dt
        B(ijk) = divu/dt
c    &           + (viscr(i+1,k) - viscr(i-1,k))/(xx3 - xx1)
c    &           + viscr(i,k)/xx2
c    &           + (viscz(i,k+1) - viscz(i,k-1))/(yy3 - yy1)
     &           - p0(i,k)/trho(i,k)/tsv(i,k)/tsv(i,k)/dt/dt
c       B(ijk) = - p0(i,k)/trho(i,k)/tsv(i,k)/tsv(i,k)/dt/dt
 120  continue


      do 200 i=1, nx-2
      do 200 k=2, ny-2
         ijk = jj(i,k)

         trhom1 = 0.5*(trho(i,k) + trho(i,k-1))
         trhom2 = 0.5*(trho(i,k+1) + trho(i,k))

         dzs = 1./(y(i,k+1)-y(i,k))
         rhoi = 0.5*( 1./trhom1 + 1./trhom2 )
         F3(ijk) = rhoi*dzs
  200 continue

      do 240 i=2, nx-2
      do 240 k=1, ny-2
         ijk = jj(i,k)

         trhom1 = 0.5*(trho(i,k) + trho(i-1,k))
         trhom2 = 0.5*(trho(i+1,k) + trho(i,k))

         drs = 1./(x(i+1,k)-x(i,k))
c        rhoi = 0.5*( 1./trhom1
c    &              + 1./trhom2 )
         rhoi = 0.5*( x(i,k)/trhom1
     &              + x(i+1,k)/trhom2 )
         F1(ijk) = drs*rhoi
  240 continue

c
c: ----- set band matrix (A) and left-hand vector (B) -----
c
c: ***** set matrix A (inner point) ************************************
c
      do 300 k=2, ny-2
      do 300 i=2, nx-2
         ijk  = jj(i,k)

         xx2 = 0.5*(x(i+1,k) + x(i,k))

         rrkps = 1./xx2
c        rrkps = 1.
         drs = 0.5*(x(i+1,k) - x(i-1,k))
         drs = 1./drs
         dzs = 0.5*(y(i,k+1) - y(i,k-1))
         dzs = 1./dzs

         pcrr = -1./trho(i,k)/tsv(i,k)/tsv(i,k)/dt/dt

         A(ijk,2) = rrkps*drs*F1(ijk-(ny-1))
         A(ijk,4) = dzs*F3(ijk-1)
         A(ijk,5) = rrkps*drs*(-F1(ijk)-F1(ijk-(ny-1)))
     &              +dzs*(-F3(ijk)-F3(ijk-1))
     &              +pcrr
         A(ijk,6) = dzs*F3(ijk)
         A(ijk,8) = rrkps*drs*F1(ijk)
  300 continue


c: ***** Wall sides boundary *******************************************

c
c: --- bottom boundary ---
c
c     do 400 i=1, ihol
      do 400 i=1, nx-1
         ijk = jj(i,1)
         A(ijk,5)  = 1.0
         A(ijk,6)  = -1.0
         B(ijk)    = 0.

         A(ijk,2) = 0.
         A(ijk,4) = 0.
         A(ijk,8) = 0.
  400 continue

c     do 410 i=ihol+1,nx-1
c        ijk = jj(i,1)
c        A(ijk,5) = 1.0
c        A(ijk,6) = -1.0
c        B(ijk) = 0.

c        A(ijk,2) = 0.
c        A(ijk,4) = 0.
c        A(ijk,8) = 0.
c 410 continue

c
c: --- top boundary ---
c
      do 450 i=1, iinj
c        k = ny-1
c        ijk = jj(i,k)
         ijk = jj(i,ny-1)

         A(ijk,5)  = 1.0
c        B(ijk)    = -uy0(i,k)*(uy_sy(i,k)-uy_sy(i,k-1))*trho(i,k)
c        B(ijk)    = 1.013e5
         B(ijk)    = 0
         A(ijk,4) = -1.

c        A(ijk,4) = 0.
         A(ijk,2) = 0.
         A(ijk,6) = 0.
         A(ijk,8) = 0.

  450 continue

c     do 500 i=iinj+1, inoz
      do 500 i=iinj+1, nx-1
         ijk = jj(i,ny-1)
         A(ijk,5)  = 1.0
         B(ijk)    = 0.
c        B(ijk)    = 1.013e5
         A(ijk,4)  =-1.

         A(ijk,2) = 0.
c        A(ijk,4) = 0.
         A(ijk,6) = 0.
         A(ijk,8) = 0.
  500 continue

c     do 550 i=inoz, nx-1
c        ijk = jj(i,ny-1)
c        A(ijk,5) = 1.
c        B(ijk) = 1.013e5
c        B(ijk) = 0.
c        A(ijk,4) = -1.
   
c        A(ijk,2) = 0.
c        A(ijk,4) = 0.
c        A(ijk,6) = 0.
c        A(ijk,8) = 0.
c 550 continue

c: ***** axis  boundary ***********************************************
      do 600 k=1, ny-1
         ijk = jj(1,k)
         A(ijk,5)  = 1.0
         A(ijk,8)  = -1.0
         B(ijk)    = 0.

         A(ijk,2) = 0.
         A(ijk,4) = 0.
         A(ijk,6) = 0.
  600 continue


c: ***** out boundary **********************************************
c     do 700 k=1, jhol
c        ijk = jj(nx-1,k)
c        A(ijk,5) = 1.0
c        A(ijk,2) = -1.0
c        B(ijk) = 0.

c        A(ijk,4) = 0.
c        A(ijk,6) = 0.
c        A(ijk,8) = 0.
c 700 continue

      do 730 k=1, ny-1
         ijk = jj(nx-1,k)
         A(ijk,5) = 1.0
         A(ijk,2) = -1.0
         B(ijk)   = 0.

         A(ijk,4) = 0.
         A(ijk,6) = 0.
         A(ijk,8) = 0.
  730 continue

c     do 770 k=jnoz+1, ny-1
c        ijk = jj(nx-1,k)
c        A(ijk,5) = 1.0
c        A(ijk,2) =-1.
c        B(ijk)   = 0.

c        A(ijk,4) = 0.
c        A(ijk,6) = 0.
c        A(ijk,8) = 0.
c 770 continue
c
c: --- nozzle boundary ---
c
      do 800 k= jnoz+1, ny-1
         ijk = jj(iinj,k)
         A(ijk,5)  = 1.0
         A(ijk,2)  =-1.0
         B(ijk)    = 0.

         A(ijk,4) = 0.
         A(ijk,6) = 0.
         A(ijk,8) = 0.
 800   continue

c     do 810 k=jnoz+1, ny-1
c        ijk = jj(inoz,k)
c        A(ijk,5) = 1.0
c        A(ijk,8) = -1.0
c        B(ijk) = 0.

c        A(ijk,2) = 0.
c        A(ijk,4) = 0.
c        A(ijk,6) = 0.
c810   continue

c     do 850 i= iinj+1, nx-1
      do 850 i= iinj, nx-1
         ijk = jj(i,jnoz)
         A(ijk,5)  = 1.0
         A(ijk,4)  =-1.0
         B(ijk)    = 0.

         A(ijk,2) = 0.
         A(ijk,6) = 0.
         A(ijk,8) = 0.
 850   continue

c     do 870 k= 1, jhol
c        ijk = jj(ihol,k)
c        A(ijk,5)  = 1.0
c        A(ijk,2)  =-1.0
c        B(ijk)    = 0.

c        A(ijk,4) = 0.
c        A(ijk,6) = 0.
c        A(ijk,8) = 0.
c870   continue

c     do 890 i= ihol+1, nx-1
c        ijk = jj(i,jhol)
c        ijk = jj(i,jhol-1)
c        A(ijk,5)  = 1.0
c        A(ijk,6)  =-1.0
c        B(ijk)    = 0.

c        A(ijk,2) = 0.
c        A(ijk,4) = 0.
c        A(ijk,8) = 0.
c890   continue


c: ----- assume initial solution -----
c
      do 900 i=1, nx-1
      do 900 k=1, ny-1
         ijk = jj(i,k)
         XX(ijk) = p0(i,k)
  900 continue


990   format(13f12.3)

c
c: ----- solve band matrix -----
c
      call DBICGST(A, N3, N1, N2, M3, M6, 
     &             B, R0, DD, PP, RR,
     &             S, T, V, XX, EPS, ITR, IER, SGM)

      do 1000 i=1, nx-1
      do 1000 k=1, ny-1
         ijk = jj(i,k)
         p0(i,k) = XX(ijk)
 1000 continue
      
      return
      end


      SUBROUTINE DBICGST(A,N3,N1,N2,M3,M6,
     &                   B,R0,D,P,R,S,T,V,XX,
     &                   EPS,ITR,IER,SGM)
************************************************************************
*                                                                      *
*   INCOMPLETE LU DECOMPOSITION BI CONJUGATED GRADIENT STABLE METHOD   *
*        FOR FINITE DIFFERENCE METHOD   (9 POINT DISCRETIZATION)       *
*                                                                      *
*                                                                      *
*     --- input ---                                                    *
*                                                                      *
*     A   .. 2-DIM. ARRAY CONTAINING THE MATRIX                        *
*     N   .. ROW SIZE OF THE MATRIX (A)                                *
*     N1  .. UPPER ROW SIZE OF THE ARRAY (A)                           *
*     N2  .. LOWER ROW SIZE OF THE ARRAY (A)                           *
*     M1  .. NUMBER OF MESH POINTS FROM (I,J,K) TO (I,J-1,K+1)         *
*     M2  .. NUMBER OF MESH POINTS FROM (I,J,K) TO (I,J,K+1)           *
*     M3  .. NUMBER OF MESH POINTS FROM (I,J,K) TO (I,J+1,K+1)         *
*     M4  .. NUMBER OF MESH POINTS FROM (I,1,K) TO (I,NT,K+1)(PERIODIC)*
*     M5  .. NUMBER OF MESH POINTS FROM (I,J,K) TO (I,J+1,K)           *
*     M6  .. NUMBER OF MESH POINTS FROM (I,J,K) TO (I+1,J,K)           *
*     B   .. 1-DIM. ARRAY CONTAINING THE RIGHT HAND SIDE VECTOR        *
*     X   .. 1-DIM. ARRAY CONTAINING INITIAL SOLUTION VECTOR           *
*     D   .. 1-DIM. WORKING ARRAY                                      *
*     P   .. 1-DIM. WORKING ARRAY                                      *
*     R   .. 1-DIM. WORKING ARRAY                                      *
*     R0  .. 1-DIM. WORKING ARRAY                                      *
*     S   .. 1-DIM. WORKING ARRAY                                      *
*     T   .. 1-DIM. WORKING ARRAY                                      *
*     V   .. 1-DIM. WORKING ARRAY                                      *
*     EPS .. TOLERANCE FOR CONVERGENCE                                 *
*     ITR .. MAXIMUM NUMBER OF ITERATION                               *
*     SGM .. METHOD CODE  (SGM=0 : ILU     SGM<>0 : MILU)              *
*                                                                      *
*     --- output ---                                                   *
*                                                                      *
*     X   .. 1-DIM. ARRAY CONTAINING THE SOLUTION VECTOR               *
*     EPS .. RELATIVE ERROR AT RETURN                                  *
*     ITR .. NUMBER OF ITERATION AT RETURN                             *
*     IER .. ERROR CODE                                                *
*                                                                      *
************************************************************************
      IMPLICIT REAL(A-H,O-Z)
      DIMENSION A(N2:N1,9), B(N3), R0(N3),
     &          D(-M3:N3+M3), P(N2:N3+M3), R(N2:N3+M3), S(N2:N3+M3),
     &          T(N2:N3+M3), V(N2:N3+M3), XX(N2:N3+M3)
C
C  PARAMETER CHECK
C
c     IF( (N1 .LT. N+M6) .OR. (M3 .LE. 1)  .OR.
c    &    (M1 .GE. M2)   .OR. (M2 .GE. M5) .OR.
c    &    (M5 .GE. M3)   .OR. (M6 .GT. N)  .OR.
c    &    (N2 .GT. -M6)  .OR. (SGM .LT. 0.0) ) THEN
c        WRITE(*,*) '(SUBR)BICGST: INVALID ARGUMENT. ',
c    &              N,N1,N2,M1,M2,M3,M4,M5,SGM
c        IER = 2
c        RETURN
c     ENDIF
C
      TH = 1.0E0
      IF( SGM .GT. 0.0 .AND. SGM .LT. 1.0 ) THEN
         TH = SGM
         SGM = 1.0E0
      ENDIF
C
C  INITIALIZATION
C
C#ASSERT DEPENDENCY(IGNORE)

      DO 10 I=1-M3, 0
         D(I) = 0.0E0
         P(I) = 0.0E0
         R(I) = 0.0E0
         S(I) = 0.0E0
         T(I) = 0.0E0
         V(I) = 0.0E0
         XX(I) = 0.0E0
         D(I+N3+M3) = 0.0E0
         P(I+N3+M3) = 0.0E0
         R(I+N3+M3) = 0.0E0
         S(I+N3+M3) = 0.0E0
         T(I+N3+M3) = 0.0E0
         XX(I+N3+M3) = 0.0E0
         V(I+N3+M3) = 0.0E0

   10 CONTINUE

      DO 20 I=1, N3
         D(I) = 0.0E0
         P(I) = 0.0E0
         R(I) = 0.0E0
         S(I) = 0.0E0
         T(I) = 0.0E0
         V(I) = 0.0E0
   20 CONTINUE

C
C  MODIFIED INCOMPLETE CHOLESKY DECOMPOSITION
C
      IF( SGM .NE. 0.0 ) THEN
         DO 30 I=1, N3
            SS = SGM*A(I,5)
     &           - A(I,4)*D(I-M6)
     &           * (A(I-M6,6) +A(I-M6,8) )
     &           - A(I,2)*D(I-M3)
     &           * (A(I-M3,6) +A(I-M3,8) )
            D(I) = 1.0E0 / SS

   30    CONTINUE
C
C  INCOMPLETE CHOLESKY DECOMPOSITION
C
      ELSE
         DO 40 I=1, N3
            SS = A(I,5) - A(I,4)*A(I-M6,6)*D(I-M6)
     &                  - A(I,2)*A(I-M3,8)*D(I-M3)
            D(I) = 1.0E0 / SS

   40    CONTINUE
      ENDIF

      DO 50 I=1, N3
         V(I) = A(I,2)*XX(I-M3) 
     &        + A(I,4)*XX(I-M6) 
     &        + A(I,5)*XX(I)
     &        + A(I,6)*XX(I+M6)  
     &        + A(I,8)*XX(I+M3) 
   50 CONTINUE
      DO 60 I=1, N3
         R(I) = B(I) - V(I)
   60 CONTINUE
C
C  INCOMPLETE LU DECOMPOSITION
C
      DO 70 I=1, N3
         R(I) = D(I)*( R(I)
     &        - A(I,4)*R(I-M6) - A(I,2)*R(I-M3) )
   70 CONTINUE
      DO 80 I=N3, 1, -1
         R(I) = R(I)
     &     - D(I)*( A(I,6)*R(I+M6) + A(I,8)*R(I+M3) )
   80 CONTINUE
      RHO = 0.0E0
      BNORM = 0.0E0
      DO 90 I=1, N3
         R0(I) = R(I)
         P(I) = R0(I)
         RHO = RHO + R0(I)*R0(I)
         BNORM = BNORM + B(I)*B(I)
   90 CONTINUE
C
C  ITERATION PHASE
C
      DO 300 K=1, ITR
         DO 100 I=1, N3
            V(I) = 
     &           + A(I,2)*P(I-M3) + A(I,4)*P(I-M6)
     &           + A(I,5)*P(I)
     &           + A(I,6)*P(I+M6) + A(I,8)*P(I+M3)
  100    CONTINUE
         DO 110 I=1, N3
            V(I) = D(I)*( V(I)
     &           - A(I,4)*V(I-M6) - A(I,2)*V(I-M3) )
  110    CONTINUE
         DO 120 I=N3, 1, -1
            V(I) = V(I)
     &      - D(I)*( A(I,6)*V(I+M6) + A(I,8)*V(I+M3) )
  120    CONTINUE
C
         RV = 0.0E0
         DO 130 I=1, N3
            RV = RV + V(I)*R0(I)
  130    CONTINUE
         ALPHA = RHO / RV
C
         DO 140 I=1, N3
            S(I) = R(I) - ALPHA*V(I)
  140    CONTINUE
C
         DO 150 I=1, N3
            T(I) = 
     &           + A(I,2)*S(I-M3)  +  A(I,4)*S(I-M6)
     &           + A(I,5)*S(I)
     &           + A(I,6)*S(I+M6)  +  A(I,8)*S(I+M3)
  150    CONTINUE
         DO 160 I=1, N3
            T(I) = D(I)*( T(I)
     &           - A(I,4)*T(I-M6)
     &           - A(I,2)*T(I-M3) )
  160    CONTINUE
         DO 170 I=N3, 1, -1
            T(I) = T(I)
     &        - D(I)*( A(I,6)*T(I+M6)+ A(I,8)*T(I+M3) )
  170    CONTINUE
C
         ST = 0.0E0
         TT = 0.0E0
         DO 180 I=1, N3
            ST = ST + S(I)*T(I)
  180    CONTINUE
         DO 190 I=1, N3
            TT = TT + T(I)*T(I)
  190    CONTINUE
         OMEGA = ST / TT
C
         RNORM = 0.0E0
         RHON = 0.0E0
         Y  = 0.0E0
         X1 = 0.0E0
         X2 = 0.0E0
         DO 200 I=1, N3
            Y = XX(I)
            XX(I) = XX(I) + ALPHA*P(I) + OMEGA*S(I)
            R(I) = S(I) - OMEGA*T(I)
            RHON = RHON + R0(I)*R(I)
            RNORM = RNORM + R(I)*R(I)
            X1 = X1 + Y*Y
            X2 = X2 + (XX(I) - Y)**2
  200    CONTINUE
C
C
C  CONVERGENCE CHECK
C
         IF( X1 .NE. 0.0 ) THEN
            RES = SQRT(X2 / X1)
            IF( RES .LE. EPS .OR. RNORM .LE. EPS*BNORM ) THEN
               ITR = K
               IER = 0
               EPS = RES
               IF( TH .NE. 1.0E0 ) SGM = TH
               RETURN
            ENDIF
         ENDIF
C
         BETA = (RHON/RHO)*(ALPHA/OMEGA)
         RHO  = RHON
         DO 210 I=1, N3
            P(I) = R(I) + BETA*(P(I) - OMEGA*V(I))
  210    CONTINUE
  300 CONTINUE
C
C  NOT CONVERGENCE
C
      IER = 1
      WRITE(*,*) '(SUBR)BICGST: NO CONVERGENCE'
       stop
      IF( TH .NE. 1.0E0 ) SGM = TH

      RETURN
      END


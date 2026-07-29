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
      subroutine potential
      include 'subcom.inc'

      parameter( n3 = (nx-1)*(ny-1))

      parameter( M3=ny-1, M6=1 )
      parameter( N1=n3+M3, N2=-M3 )
      dimension F1(n3), F2(n3), F3(n3)
      dimension A(N2:N1,9)
      dimension B(n3), R0(n3)
      dimension DD(N2:N1), PP(N2:n3+M3), RR(N2:n3+M3)
      dimension S(N2:n3+M3), T(N2:n3+M3), V(N2:n3+M3)
      dimension XX(N2:n3+M3)

      dimension teps(nx,ny)


      ITR = 10000
c     EPS = 1.0E-40   ! 到達不能な値。毎ステップ上限反復を回し切り極端に遅かった
      EPS = 1.0E-6
      SGM = 0.0
 
      eps0  = 8.8542e-12
      epsil = 10.
c     epsil = 9.5
      epsig = 1.
c     epsil = 9.5
c     epsig = 1.


c     sigl = 1.4
      sigl = 0.277
c     sigl = 0.365
c     sigg = 2.9e-15
c     sigg = 0.0025
      sigg = 2.5e-10
c     sigg = 1.4


c     vanode = 0.7e3
c     vanode = 1.
c     vanode = -6.e3
c     vanode = 7.e3
      vanode = 3.e3
c     vanode = 0.01e3
c     vanode = 1.e3

c
c: --- initialize ---
c
      do 100 i=1, n3
         F1(i)  = 0.0
         F2(i)  = 0.0
         F3(i)  = 0.0
         B(i)   = 0.0
         phi(i) = 0.0
  100 continue

c     do 101 i=1, nx
c     do 101 j=1, ny
c       if( vf0(i,j) .gt. 1.0 ) vf0(i,j) = 1.0
c       if( vf0(i,j) .lt. 0. ) vf0(i,j) = 0.

c       if( vf_sx(i,j) .gt. 1.0 ) vf_sx(i,j) = 1.0
c       if( vf_sx(i,j) .lt. 0. ) vf_sx(i,j) = 0.

c       if( vf_sy(i,j) .gt. 1.0 ) vf_sy(i,j) = 1.0
c       if( vf_sy(i,j) .lt. 0. ) vf_sy(i,j) = 0.
c101   continue

      do 102 i=2, nx-2
      do 102 k=2, ny-2
c     sdx  = x(i+1,k) - x(i,k)
c     sdy  = y(i,k+1) - y(i,k)

c     alph = sqrt(sdx**2. + sdy**2.)
c     xhev = 0.5*max(-1.,
c    &           min(1.,dphi(i,k)/alph+1./pi*sin(pi*dphi(i,k)/alph)))

c     epsi(i,k) = 0.5*(epsil+epsig)+(epsil-epsig)*xhev
c     sig(i,k)  = 0.5*(sigl+sigg)+(sigl-sigg)*xhev

c     if( vf0(i,k) .ge. 0.5 ) then
c     alph = sqrt(sdx**2. + sdy**2.)
       
c     if( vf0(i,k) .ge. 0.5 .and. dphi(i,k) .ge. 0.  ) then
c     if( vf0(i,k) .ge. 0.5 ) then
c     if( dphi(i,k) .ge. 0. ) then
c     if( k .ge. jnoz .and. i .le. iinj ) then
c        epsi(i,k) = epsil
c        sig(i,k) = sigl
c     else
c        epsi(i,k) = epsig
c        sig(i,k) = sigg
c     end if

102   continue

      do 103 i=1, nx-1
      epsi(i,ny-1) = epsi(i,ny-2)
      sig(i,ny-1) = sig(i,ny-2)
103   continue

      do 104 k=1, ny-1
      epsi(nx-1,k) = epsi(nx-2,k)
      sig(nx-1,k) = sig(nx-2,k)
104   continue

c     do 105 i=1, nx
c     do 105 k=1, ny
      do 105 i=1, nx-1
      do 105 k=1, ny-1


c       if( vf0(i,k) .ge. 0.9 ) then
c          sig(i,k) = 1.4
c       else
c          sig(i,k) = 2.9e-15
c          sig(i,k) = 7.88e-15
c          sig(i,k) = 0.025
c       end if

c       if( i .le. inoz .and.  k .ge. jnoz ) then
c          epsi(i,k) = 2.
c          sig(i,k) = 6.e7
c       end if

c       if( i .ge. ihol .and. k .ge. jhol2 .and. k .le. jhol1 ) then
c          epsi(i,k) = 2.
c          sig(i,k) = 6.e7
c       end if

  105 continue

c#COLLAPSE
      do 110 i=1, 9
         do 110 j=N2, N1
            A(j,i)= 0.0
  110 continue


c
c: ----- set coefficient Fk and Pe -----
c

      do 120 i=1, nx-1
      do 120 k=1, ny-1 
        ijk = jj(i,k)

        divrh = (rhoc0(i+1,k)-rhoc0(i,k))/(x(i+1,k)-x(i,k))
     &        + rhoc0(i,k)/x(i,k)
     &        + (rhoc0(i,k+1)-rhoc0(i,k))/(y(i,k+1)-y(i,k))

c       B(ijk) = - divrh / eps0
c       B(ijk) = 0.
c       B(ijk) = - rhoc0(i,k)/eps0
c       B(ijk) = - rhoc0(i,k)*1.e2 / eps0
 120   continue


c     do 200 i=1, nx-2
      do 200 i=1, nx-1
      do 200 k=1, ny-2
         ijk = jj(i,k)
c        yy2 = 0.5*(y(i,k+1) + y(i,k))
c        yy3 = 0.5*(y(i,k+2) + y(i,k+1))

         dzs = 1./(y(i,k+1)-y(i,k))
c        diel = 0.5*( epsi(i,k) + epsi(i,k+1) )
c        diel = 0.5*( x(i,k) + x(i,k+1) )
         diel = epsig*eps0

c        diel = 1.
         F3(ijk) = diel*dzs
  200 continue

      do 240 i=1, nx-2
c     do 240 k=1, ny-2
      do 240 k=1, ny-1
         ijk = jj(i,k)
c        xx2 = 0.5*(x(i+1,k) + x(i,k))
c        xx3 = 0.5*(x(i+2,k) + x(i+1,k))

c        drs = 1./(xx3 - xx2)
         drs = 1./(x(i+1,k)-x(i,k))
c        diel = 0.5*( xx2*sig(i,k)
c    &              + xx3*sig(i+1,k) )
c        diel = 0.5*( x(i,k)*epsi(i,k)
c    &              + x(i+1,k)*epsi(i+1,k) )
c        diel = 0.5*( x(i,k)
c    &              + x(i+1,k) )
c        円筒(2次元軸対称): diel = 0.5*(x(i,k)+x(i+1,k))*epsig*eps0
c        平面2次元(デカルト xy): r重みなし
         diel = epsig*eps0
c        diel = 1.
         F1(ijk) = drs*diel
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

c        rrkps = 1./xx2   ! 円筒(2次元軸対称): 1/r 係数
         rrkps = 1.       ! 平面2次元(デカルト xy)
c        drs = 0.5*(x(i+1,k) - x(i-1,k))
         drs = sdx1(i,k)
         drs = 1./drs
c        dzs = 0.5*(y(i,k+1) - y(i,k-1))
         dzs = sdy1(i,k)
         dzs = 1./dzs

c        xx1 = 0.5*(x(i,k) + x(i-1,k))
c        xx2 = 0.5*(x(i+1,k) + x(i,k))
c        xx3 = 0.5*(x(i+2,k) + x(i+1,k))
c        yy1 = 0.5*(y(i,k) + y(i,k-1))
c        yy2 = 0.5*(y(i,k+1) + y(i,k))
c        yy3 = 0.5*(y(i,k+2) + y(i,k+1))

c        rrkps = 1. / xx2
c        drs = 0.5*(xx3 - xx1)
c        drs = 1./drs
c        dzs = 0.5*(yy3 - yy1)
c        dzs = 1./dzs

         A(ijk,2) = rrkps*drs*F1(ijk-(ny-1))
         A(ijk,4) = dzs*F3(ijk-1)
         A(ijk,5) = rrkps*drs*(-F1(ijk)-F1(ijk-(ny-1)))
     &              +dzs*(-F3(ijk)-F3(ijk-1))
         A(ijk,6) = dzs*F3(ijk)
         A(ijk,8) = rrkps*drs*F1(ijk)
  300 continue



c: ***** Wall sides boundary *******************************************

c
c: --- bottom boundary ---
c

      do 400 i=1, nx-1
c     do 400 i=1, iinj-1
         ijk = jj(i,1)
         A(ijk,5)  = 1.0
         A(ijk,6)  =-1.0
         B(ijk)    = 0.

         A(ijk,2) = 0.
c        A(ijk,6) = 0.
         A(ijk,4) = 0.
         A(ijk,8) = 0.
  400 continue

c
c: --- top boundary ---
c

c     do 500 i=1, inoz-1
c     do 500 i=1, iinj-1
      do 500 i=1, nx-1
         ijk = jj(i,ny-1)
         A(ijk,5) = 1.0
         A(ijk,4) =-1.0
         B(ijk)   = 0.

         A(ijk,2) = 0.
         A(ijk,6) = 0.
         A(ijk,8) = 0.
c        A(ijk,4) = 0.
  500 continue

c     do 550 i=inoz, nx-1
c        ijk = jj(i,ny-1)
c        A(ijk,5) = 1.0
c        A(ijk,4) = -1.0
c        B(ijk)   = 0.

c        A(ijk,2) = 0.
c        A(ijk,6) = 0.
c        A(ijk,8) = 0.
c 550 continue

c     --- axis center ---
c     do 600 k=1, jnoz-1
      do 600 k=1, ny-1
         ijk = jj(1,k)
         A(ijk,5)  = 1.0
         A(ijk,8)  =-1.0
         B(ijk)    = 0.

         A(ijk,2) = 0.
         A(ijk,4) = 0.
         A(ijk,6) = 0.
  600 continue

c     do 650 k= jnoz, ny-1
c        ijk = jj(1,k)
c        A(ijk,5) = 1.0
c        A(ijk,8) = -1.0
c        B(ijk)   = 0.

c        A(ijk,2) = 0.
c        A(ijk,4) = 0.
c        A(ijk,6) = 0.
c 650 continue

c     --- out side of computational area --- 
c     do 700 k=1, jhol2-1
c     do 700 k=1, jnoz-1
      do 700 k=1, ny-1
         ijk = jj(nx-1,k)
         A(ijk,5) = 1.0
         A(ijk,2) = -1.0
         B(ijk)   = 0.

c        A(ijk,2) = 0.
         A(ijk,4) = 0.
         A(ijk,6) = 0.
         A(ijk,8) = 0.
  700 continue

c     do 705 k=jhol2, jhol1-1
c        ijk = jj(i,k)
c        A(ijk,5) = 1.0
c        A(ijk,2) = -1.0
c        B(ijk)   = 0.

c        A(ijk,4) = 0.
c        A(ijk,6) = 0.
c        A(ijk,8) = 0.
c 705 continue

c     do 710 k=jhol1, ny-1
c        ijk = jj(nx-1,k)
c        A(ijk,5) = 1.0
c        A(ijk,2) = -1.0
c        B(ijk)   = 0.

c        A(ijk,2) = 0.
c        A(ijk,4) = 0.
c        A(ijk,6) = 0.
c        A(ijk,8) = 0.
c 710 continue

c     ------ nozzle ------
c     --- bottom nozzle ---
      do 720 i=1, inoz-1
         k = jnoz
         ijk = jj(i,k)
         A(ijk,5) = 1.0
         A(ijk,6) = -1.0
         B(ijk)   = 0.

         A(ijk,2) = 0.
         A(ijk,4) = 0.
         A(ijk,8) = 0.
  720 continue     

c     --- outside nozzle ---
      do 730 k=jnoz, ny-1
         i = inoz
         ijk = jj(i,k)
         A(ijk,5) = 1.0
         A(ijk,2) = -1.0

         A(ijk,8) = 0.
         A(ijk,4) = 0.
         A(ijk,6) = 0.
  730 continue 

c     --- inside nozzle ---
      do 740 i=1, inoz
      do 740 k=jnoz, ny-1
         ijk = jj(i,k)
         A(ijk,5) = 1.0
c        B(ijk) = vanode   ! 旧: ノズル +7kV
         B(ijk) = 0.       ! ノズルは接地(0V)

         A(ijk,2) = 0.
         A(ijk,4) = 0.
         A(ijk,6) = 0.
         A(ijk,8) = 0.
  740 continue

c     ------ nozzle 2 (mirror image, right boundary) ------
c     --- bottom nozzle 2 ---
      do 721 i=inoz2+1, nx-1
         k = jnoz
         ijk = jj(i,k)
         A(ijk,5) = 1.0
         A(ijk,6) = -1.0
         B(ijk)   = 0.

         A(ijk,2) = 0.
         A(ijk,4) = 0.
         A(ijk,8) = 0.
  721 continue

c     --- inside face of nozzle 2 ---
      do 731 k=jnoz, ny-1
         i = inoz2
         ijk = jj(i,k)
         A(ijk,5) = 1.0
         A(ijk,8) = -1.0

         A(ijk,2) = 0.
         A(ijk,4) = 0.
         A(ijk,6) = 0.
  731 continue

c     --- inside nozzle 2 ---
      do 741 i=inoz2, nx-1
      do 741 k=jnoz, ny-1
         ijk = jj(i,k)
         A(ijk,5) = 1.0
c        B(ijk) = vanode   ! 旧: ノズル +7kV
         B(ijk) = 0.       ! ノズルは接地(0V)

         A(ijk,2) = 0.
         A(ijk,4) = 0.
         A(ijk,6) = 0.
         A(ijk,8) = 0.
  741 continue

c     ------ counter electrode (holes mirrored on both sides) ------
c     --- bottom  ---
      do 750 i=ihol, ihol2
         k = jhol2
         ijk = jj(i,k)
         A(ijk,5) = 1.0
         A(ijk,4) = -1.0
c        A(ijk,6) = -1.0
         B(ijk)   = 0.

         A(ijk,2) = 0.
c        A(ijk,4) = 0.
         A(ijk,6) = 0.
         A(ijk,8) = 0.
  750 continue

c     --- side ---
      do 760 k=jhol2, jhol1
         i = ihol
         ijk = jj(i,k)
         A(ijk,5) = 1.0
         A(ijk,8) = -1.0
         B(ijk)   = 0.

         A(ijk,4) = 0.
         A(ijk,6) = 0.
         A(ijk,2) = 0.
  760 continue

c     --- side (right hole, mirror of do760) ---
      do 761 k=jhol2, jhol1
         i = ihol2
         ijk = jj(i,k)
         A(ijk,5) = 1.0
         A(ijk,2) = -1.0
         B(ijk)   = 0.

         A(ijk,4) = 0.
         A(ijk,6) = 0.
         A(ijk,8) = 0.
  761 continue

c     --- top ---
      do 770 i=ihol, ihol2
         k = jhol1
         ijk = jj(i,k)
         A(ijk,5) = 1.0
         A(ijk,6) = -1.0
c        A(ijk,4) = -1.0
         B(ijk)   = 0.

         A(ijk,2) = 0.
         A(ijk,4) = 0.
c        A(ijk,6) = 0.
         A(ijk,8) = 0.
  770 continue

c     --- inside ---
      do 780 i=ihol, ihol2
      do 780 k=jhol2, jhol1
         ijk = jj(i,k)
         A(ijk,5) = 1.0
c        B(ijk) = 0.        ! 旧: 対向電極 接地
         B(ijk) = -vanode   ! 対向電極 -7kV

         A(ijk,2) = 0.
         A(ijk,4) = 0.
         A(ijk,6) = 0.
         A(ijk,8) = 0.
  780 continue


c     iedg = iinj+6
c     do 750 k=jnoz, ny-1
c     if( k .le. jnoz+5 ) then
c        istrt = iedg
c        iedg = iedg - 1
c     else
c        istrt = iinj
c     end if
c     do 750 i=iinj, nx-1

c     do 750 i=istrt, nx-1
c        ijk = jj(i,k)
c        A(ijk,5) = 1.0
c        B(ijk)   = vanode

c        A(ijk,2) = 0.
c        A(ijk,4) = 0.
c        A(ijk,6) = 0.
c        A(ijk,8) = 0.
c 750 continue

c     do 760 k=jnoz, ny-1
c        i = iinj
c        ijk = jj(i,k)
c        A(ijk,5) = 1.0
c        B(ijk)   = 0.

c        A(ijk,8) = -1.0
c        A(ijk,2) = 0.
c        A(ijk,4) = 0.
c        A(ijk,6) = 0.
c 760 continue

c     do 770 i=iinj, nx-1
c        k = jnoz
c        ijk = jj(i,k)
c        A(ijk,5) = 1.0
c        B(ijk)   = 0.

c        A(ijk,6) = -1.0
c        A(ijk,2) = 0.
c        A(ijk,4) = 0.
c        A(ijk,8) = 0.
c 770 continue

c     do 800 i=2, iinj-1
c     do 800 i=2, iinj-9
c     do 800 k=2, jnoz+2
c       ijk = jj(i,k)
c       alp = 1.e-5
c       alp = sdy(i,k)
c       if( abs(dphi(i,k)) .le. alp .and. vf0(i,k) .ge. 0.5 ) then
c          teps(i,k)  =  (1.-vf0(i,k))*epsig+vf0(i,k)*epsil
c          sdyn = y(i,k+1)-y(i,k)           
c          sdys = y(i,k)-y(i,k-1)           
c          sdxe = x(i+1,k)-x(i,k)           
c          sdxw = x(i,k)-x(i-1,k)           
 
 
c          gval = sqrt(gradx(i,k)*gradx(i,k)+grady(i,k)*grady(i,k))  
c          xnr = gradx(i,k)/gval
c          xnz = grady(i,k)/gval

c          A(ijk,5) = teps(i,k)/sdyn+teps(i,k)/sdys
c          A(ijk,6) = -teps(i,k)/sdyn
c          A(ijk,4) = -teps(i,k)/sdys
c          A(ijk,2) = 0.
c          A(ijk,8) = 0.

c          B(ijk) = rhoc0(i,k)/eps0

c          A(ijk,5) = A(ijk,5) + teps(i,k)/sdxw+teps(i,k)/sdxe
c          A(ijk,2) = -teps(i,k)/sdxw
c          A(ijk,8) = -teps(i,k)/sdxe
c          A(ijk,6) = 0.
c          A(ijk,4) = 0.

c          write(*,*) A(ijk,5), A(ijk,2), A(ijk,8)

c          B(ijk) = rhoc0(i,k)/eps0
c          B(ijk) = 0.

c       end if
c800   continue


c
c: ----- assume initial solution -----
c
      do 900 i=1, nx-1
      do 900 k=1, ny-1
         ijk = jj(i,k)
         XX(ijk) = phi(ijk)
  900 continue



990   format(13f12.3)

c
c: ----- solve band matrix -----
c
      call DBICGST(A, N3, N1, N2, M3, M6, 
     &             B, R0, DD, PP, RR,
     &             S, T, V, XX, EPS, ITR, IER, SGM)
c

      do 1000 i=1, nx-1
      do 1000 k=1, ny-1
         ijk = jj(i,k)
         phi(ijk) = XX(ijk)
 1000 continue

c     ----- 左右対称化 -----
c     配置は左右対称(同一ノズル・同電圧)なので真の解も左右対称。
c     反復ソルバが残す左右非対称の数値誤差を、鏡像 i<->nx-i の
c     平均で除去する。（左右非対称な設定にする場合はこのループを外す）
      do 1100 k=1, ny-1
      do 1100 i=1, nx-1
         im = nx - i
         if( im .gt. i )then
            ijk1 = jj(i,k)
            ijk2 = jj(im,k)
            ph = 0.5*( phi(ijk1) + phi(ijk2) )
            phi(ijk1) = ph
            phi(ijk2) = ph
         end if
 1100 continue

      return
      end




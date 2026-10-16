      subroutine initial
      include 'subcom.inc'

      real vp(100000)

c     --- particle initial position and velocity ---

c     dp1 = 1.e-6
c     dp2 = 1.e-6

c     ip = 10
      pi = 4.*atan(1.)


      do 10 k=1,ipmx

      xp(k) = 0.
c     xp(k) = 1.e-3
c     xp(k) = 2.e-3
c     yp(k) = 98.5e-3
      yp(k) = 98.2e-3

c     uxp(k) = 2.
      uxp(k) = 0.
c     uyp(k) = 0.
c     uyp(k) = -9.
      uyp(k) = -2.4
c     uyp(k) = -2.5

 10   continue

c      write(*,*) xp(1), yp(1), uxp(1), uyp(1)

      eps0 = 8.8542e-12
      epsig = 1.
 
      do 20 i=1, nx-1
      do 20 j=1, ny-1
           epsi(i,j) = epsig*eps0           
 20   continue

      do 30 k=1, ipmx

      pout(k) = 0.

 30   continue

      do 40 k=1,ipmx

      cl(k) = 0.
      clt(k) = 0.

 40   continue

      do 50 k=1,ipmx

      alp = 0.28
      pi = 4.*atan(1.)
      rhol = 1027.

      vp(k) = rhol*pi*dp(k)*dp(k)*dp(k)/6.

      cltmx(k) = alp*vp(k)*44.01/170.21
c     cltmx(k) = cltmx(k)*1.e7
 
c     cltmx(k) = 1.e-1 
c     write(*,*) cltmx(1)

      cltmx(k) = cltmx(k)/1.977

c     write(*,*) cltmx(1)


 50   continue

c     write(*,*) cltmx(1)

      do 60 i=1, nx-1
      do 60 j=1, ny-1

      cg(i,j) = 0.445
      cgn(i,j) = 0.445

c     cg(i,j) = 1.
c     cgn(i,j) = 1.

 60   continue


      vtotal = 0.   

      do 70 k=1, ipmx

      pi = 4.*atan(1.)
c     rhol = 1240.

      vp(k) = pi*dp(k)*dp(k)*dp(k)/6.
c     vp(k) = rhol*pi*dp(k)*dp(k)*dp(k)/6.


c     --- Q = 0.1mL/h = 2.78e-11m3/s , time = 1.e-3s --- 
c     --- Q = 0.01mL/h = 2.78e-12m3/s , time = 1.e-3s --- 
c     --- Q = 0.001mL/h = 2.78e-13m3/s , time = 1.e-4s --- 
c     --- Q = 0.0001mL/h = 2.78e-14m3/s , time = 1.e-4s --- 
c     qtotal = 2.8e-14
c     qtotal = 2.78e-15
c     qtotal = 2.78e-16
c     qtotal = 2.78e-17
c     qtotal = 5.56e-17
c     qtotal = 8.34e-17

c     qtotal = 5.56e-17
      qtotal = 1.11e-16


      vtotal = vtotal + vp(k)
c     write(*,*) vtotal

      if( vtotal .ge. qtotal )then
      ipmx2 = k
      write(*,*) ipmx2, qtotal, vtotal
      go to 80
      end if

 70   continue

 80   continue

c     ipmx = ipmx2
c     write(*,*) ipmx

c    eps0 = 8.8542e-12
c    epsig = 1.
 
c    epsif = ( 3.*epsig*(epsil-epsig)*eps0 )/( 2*epsig+epsil )

c    q1 = 8*pi*pi*abs(epsif)*fs*dp(1)*dp(1)*dp(1)
c    q1 = sqrt( q1 )  


c     do 20 i=1, inoz
c     do 20 j=jnoz, ny-1
c          epsi(i,j) = eps0
c20   continue

c     do 30 i=ihol, nx-1
c     do 30 j=jhol2, jhol1
c          epsi(i,j) = eps0
c30   continue

c     do 50 i=0, nx+1
c     do 50 j=0, ny+1
c          pp_sx(i,j) = 1.013e5
c          pp_sy(i,j) = 1.013e5
c          pp0(i,j) = 1.013e5
c          p0(i,j) = 1.013e5
c50    continue

c     do 70 i=0, nx+1
c     do 70 j=0, ny+1
c          rhoc0_sx(i,j) = 0.
c          rhoc0_sy(i,j) = 0.
c          rhoc0(i,j) = 0.
c70   continue

      do 100 i=1, nx
      do 100 j=1, ny
           ux_sx(i,j) = 0.
           uy_sx(i,j) = 0.
 100    continue
           

      do 200 i=1, nx
      do 200 j=1, ny
           ux_sy(i,j) = 0.
           uy_sy(i,j) = 0.
 200    continue

c     do 300 i=0, nx+1
c     do 300 j=0, ny+1
c       xx = 0.5*(x(i,j)+x(i+1,j))
c       yy = 0.5*(y(i,j)+y(i,j+1))
c       if( xx .le. 0.00025 .and. xx .gt. 0.
c    &      .and.  yy .le. 0.0025 .and. yy .gt. 0.0018) then


c       if( i .le. iinj+1 .and.  j .le. jnoz ) then
c       if(  j .ge. jnoz-5 ) then
c       if( i .le. iinj .and.  j .ge. jnoz  .or.
c    &      x(i,j)**2+(y(i,j)-0.00121)**2 .le. (0.000121)**2  ) then

c       aa = 0.1e-4
c       rad = 0.5*(x(iinj,jnoz)*x(iinj,jnoz)+aa*aa)/aa
c       yy1 = y(iinj,jnoz)+(rad-aa)
c       dis2 = (x(i,j)-x(1,jnoz))**2.+(y(i,j)-yy1)**2.
c       if( dis2 .le. rad*rad ) then
c       if( i .le. iinj+1 .and.  j .ge. ny-2 ) then
c       if( i .le. iinj+1 .and.  j .ge. jnoz ) then
c          vf0(i,j) = 1.0
c          trho(i,j) = rhol
c          tmu(i,j) = xmul
c          tsig(i,j) = sigl
c       else
c          vf0(i,j) = 0.
c          trho(i,j) = rhog
c          tmu(i,j) = xmug
c          tsig(i,j) = sigg
c       end if


c       ux0(i,j)  = 0.25*( ux_sx(i,j)+ux_sx(i+1,j)
c    &                  +  ux_sy(i,j)+ux_sy(i,j+1) )
c       uy0(i,j)  = 0.25*( uy_sx(i,j)+uy_sx(i+1,j)
c    &                  +  uy_sy(i,j)+uy_sy(i,j+1) )

c       ux0(i,j) = 0.
c       uy0(i,j) = 0.

c       end if
c300    continue


      return
      end


      subroutine intraction
      include 'subcom.inc'

      real q(100000), vp(100000), vq(100000)
      real rx(100000), ry(100000), rr1(100000), rr3(100000)
      real ftx(100000), fty(100000)
      real rx1(100000), ry1(100000), rr2(100000), rr6(100000)
      real ftx1(100000), fty1(100000)
      real ftxt1(100000), ftyt1(100000)
c     real ftxt(100000), ftyt(100000)

      eps0 = 8.8542e-12
      epsil = 10.
c     epsil = 9.5
      epsig = 1.
      
      fs = 42.9e-3
c     fs = 30.e-3
      pi = 4.*atan(1.)

      epsif = ( 3.*epsig*(epsil-epsig)*eps0 ) / ( 2.*epsig+epsil )


      do 2 k=1, ip      

      if( pout(k) .ne. 0. )then

c     if( pout(k) .eq. 1. )then
      q(k) = 0.

      else 
c     q(k) = 8.*pi*pi*abs(epsif)*fs*dp(k)*dp(k)*dp(k)
c     q(k) = sqrt( q(k) )

      q(k) = 8.*eps0*fs*dp(k)*dp(k)*dp(k)
      q(k) = pi*sqrt( q(k) )

      end if

 2    continue

      do 7 k=1, ip

      ftxt(k) = 0.
      ftyt(k) = 0.

      ftx(k) = 0.
      fty(k) = 0.

      ftxt1(k) = 0.
      ftyt1(k) = 0.

      ftx1(k) = 0.
      fty1(k) = 0.


 7    continue

c     --- relative position ---

      do 10 k=1, ip
      do 20 l=1, ip

c     if( l .ne. k )then
      if( l .eq. k )then

      ftx(l) = 0.
      fty(l) = 0.

      else

      rx(l) = xp(k) - xp(l)
      ry(l) = yp(k) - yp(l)
      
      rr1(l) = sqrt( rx(l)*rx(l) + ry(l)*ry(l) )

          if( rr1(l) .le. dp(k)+dp(l) )then
          rr1(l) = dp(k) + dp(l)
          else
          rr1(l) = rr1(l)
          end if

c     --- intraction ---

      rr3(l) = rr1(l)*rr1(l)*rr1(l)


      ftx(l) = q(l) / rr3(l) 
      ftx(l) = ftx(l)*( xp(k) - xp(l) )


      fty(l) = q(l) / rr3(l)
      fty(l) = fty(l)*( yp(k) - yp(l) )


      ftxt(k) = ftxt(k) + ftx(l)
      ftyt(k) = ftyt(k) + fty(l)

      end if

 20   continue
 10   continue



      do 22 k=1, ip
      do 24 l=1, ip

      if( l .eq. k )then

      ftx1(l) = 0.
      fty1(l) = 0.

      else

      rx1(l) = xp(k) + xp(l)
      ry1(l) = yp(k) - yp(l)

      rr2(l) = sqrt( rx1(l)*rx1(l) + ry1(l)*ry1(l) )

          if( rr2(l) .le. dp(k)+dp(l) )then
          rr2(l) = dp(k) + dp(l)
          else
          rr2(l) = rr2(l)
          end if

      rr6(l) = rr2(l)*rr2(l)*rr2(l)

      ftx1(l) = q(l) / rr6(l)
      ftx1(l) = ftx1(l)*( xp(k) + xp(l) )

      fty1(l) = q(l) / rr6(l)
      fty1(l) = fty1(l)*( yp(k) - yp(l) )

      ftxt1(k) = ftxt1(k) + ftx1(l)
      ftyt1(k) = ftyt1(k) + fty1(l)

      end if

 24   continue
 22   continue      


      do 26 k=1, ip

      ftxt(k) = ftxt(k) + ftxt1(k)
      ftyt(k) = ftyt(k) + ftyt1(k)

 26   continue
          

      do 30 k=1, ip

      fs = 42.9e-3
c     fs = 30.e-3
      epsig = epsig
      epsil = epsil
      rhol = 1027.
      pi = 4.*atan(1.)

      epsif = ( 3.*epsig*(epsil-epsig)*eps0 ) / ( 2.*epsig+epsil )
      vp(k) = pi*dp(k)*dp(k)*dp(k)/6.
c     vp(k) = rhol*pi*dp(k)*dp(k)*dp(k)/6.

c     vq(k) = 48.*pi*fs*abs(epsif)/vp(k)
c     vq(k) = sqrt( vq(k) )
c     vq(k) = vq(k) / rhol
c     write(*,*) vq(k)

      vq(k) = 8.*eps0*fs*dp(k)*dp(k)*dp(k)
      vq(k) = pi*sqrt( vq(k) )
      vq(k) = vq(k) / vp(k)
      vq(k) = vq(k) / rhol

 30   continue

      do 40 k=1, ip

c     if( pout(k) .ne. 0. )then
c     if( pout(k) .eq. 1. )then

c     uxp(k) = 0.
c     uyp(k) = 0.
c     go to 45

c     else if( pout(k) .eq. 2. )then

c     uxp(k) = 0.
c     uyp(k) = 0.
c     go to 45

c     else

c     uxp(k) = uxp(k) + dt*vq(k)*ftxt(k)/4./pi/eps0
c     uyp(k) = uyp(k) + dt*vq(k)*ftyt(k)/4./pi/eps0

c     end if

c45   continue     
 40   continue


c     xlen = 1.e-1
c     ylen = 1.e-1

      do 50 k=1, ip

c     if( pout(k) .ne. 1. )then
c     if( pout(k) .eq. 1. )then

c     xp(k) = xp(k)
c     yp(k) = yp(k)

c     else if( pout(k) .eq. 2. )then

c     xp(k) = xlen
c     yp(k) = 0.

c     else if( xp(k) + dt*uxp(k) .lt. 0. )then

c     xp(k) = 0.
c     xp(k) = xp(k)
c     yp(k) = yp(k) + dt*uyp(k)

c     else if( xp(k) .lt. 0. )then

c     xp(k) = 0.
c     yp(k) = yp(k) + dt*uyp(k)

c     else

c     xp(k) = xp(k) + dt*uxp(k)
c     yp(k) = yp(k) + dt*uyp(k)

c     end if

 50   continue


c     call pcount


      return 

      end




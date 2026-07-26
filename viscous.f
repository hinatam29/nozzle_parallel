
      subroutine viscous
      include 'subcom.inc'

       real kB
       real kn(100000), Cc(100000)
       real Rep(100000), Cdf(100000)
       real Cd(100000)
       real sp(100000)
c      real vp(100000)
c      real fvisx(100000), fvisy(100000)

       kB = 1.3806e-23
       T = 297.
       p1 = 1.013e5
       pi = 4.*atan(1.0)

c      call pcount

c      --- kn --- 
       do 10 k=1,ip

       kB = 1.3806e-23
       T = 297.
       p1 = 1.013e5
       pi = 4.*atan(1.0)

       kn(k) = sqrt(8.)*kB*T 
c      kn(k) = kn(k) / ( pi*p1*dp(k)*dp(k) )
       kn(k) = kn(k) / ( pi*p1*dp(k)*dp(k)*dp(k) )

       Cc(k) = 1. / ( 1.+kn(k)*(1.165+0.483*exp(-0.997/kn(k))) )

 10    continue


c      ---particle reynolds number---

       do 20 k=1,ip

       Rep(k) = rhog*dp(k)
     &        * sqrt( uxp(k)*uxp(k)+uyp(k)*uyp(k) ) / xmug


       if( Rep(k) .lt. 0.2 )then

       Cdf(k) = 24./Rep(k)

       else if( Rep(k) .ge. 0.2 .or. Rep(k) .lt. 2. )then

       Cdf(k) = 24./Rep(k)
     &     *( 1.+3.*Rep(k)/16. )

       else if( Rep(k) .gt. 2. .or. Rep(k) .lt. 21. )then

       Cdf(k) = 24./Rep(k)
     &     *( 1.+0.11*(Rep(k)**(0.81)) ) 

       else if( Rep(k) .gt. 21. .or. Rep(k) .lt. 200. )then

       Cdf(k) = 24./Rep(k)
     &     *( 1.+0.189*(Rep(k)**(0.62)) ) 

       end if 


 20    continue


c      ---resistance coefficient---

       do 30 k=1,ip

       Cd(k) = Cdf(k)*Cc(k)

 30    continue

c      write(*,*) Cd(1), Cdf(1), Cc(1)

c      ---viscous drag---

       do 40 k=1, ip

       rhog = 1.25
       sp(k) = pi*dp(k)*dp(k)/4.
c      sp1 = pi*dp1*dp1/4.
c      sp2 = pi*dp2*dp2/4.

c      write(*,*) sp(k)

       fvisx(k) = 0.5*Cd(k)*sp(k)*rhog
     &          * sqrt( uxp(k)*uxp(k)+uyp(k)*uyp(k) )
     &          * uxp(k)
       fvisy(k) = 0.5*Cd(k)*sp(k)*rhog
     &          * sqrt( uxp(k)*uxp(k)+uyp(k)*uyp(k) )
     &          * uyp(k)

 40    continue

c      write(*,*) fvisy(1)

c      ---next particle velocity---
       do 50 k=1, ip

c      if( pout(k) .ne. 0. )then
c      if( pout(k) .eq. 1.  .or.  pout(k) .eq. 2. )then

c      uxp(k) = 0.
c      uyp(k) = 0.
c      go to 55

c      else

c      rhol = 1027.
c      pi = 4.*atan( 1. )
c      vp(k) = rhol*pi*dp(k)*dp(k)*dp(k)/6. 

c      uxp(k) = uxp(k) - dt/vp(k)*fvisx(k)
c      uyp(k) = uyp(k) - dt/vp(k)*fvisy(k)

c      end if

c55    continue
 50    continue

c      write(*,*) fvisy(1)*dt/vp(1)


c      ---next particle position---
c      xlen = 1.e-1
c      ylen = 1.e-1

       do 60 k=1, ip

c      if( pout(k) .ne. 0. )then
c      if( pout(k) .eq. 1. )then

c      xp(k) = xp(k)
c      yp(k) = yp(k)

c      else if( pout(k) .eq. 2. )then

c      xp(k) = xlen
c      yp(k) = 0.

c      else if( xp(k) + dt*uxp(k) .lt. 0. )then

c      xp(k) = xp(k)
c      yp(k) = yp(k) + dt*uyp(k)

c      else if( xp(k) .lt. 0. )then
      
c      xp(k) = xlen
c      yp(k) = 0.

c      else if( uyp(k) .gt. 0. )then

c      xp(k) = xp(k) - dt*uxp(k)
c      yp(k) = yp(k) - dt*uxp(k)

c      else 

c      xp(k) = xp(k) + dt*uxp(k)
c      yp(k) = yp(k) + dt*uyp(k)

c      end if

 60    continue


c      call pcount

c      write(*,*) xp(1), yp(1), uxp(1), uyp(1)

      return
      end
       

                     


      subroutine eforce
      include 'subcom.inc'

      real vp(100000), vq(100000)
      real inex(100000), jnex(100000)
      real l1(100000), l2(100000), l3(100000), l4(100000)
      real Ernex(100000), Eznex(100000)
c     real fer(100000), fez(100000)

      eps0 = 8.8542e-12
      epsil = 10.
c     epsil = 9.5
      epsig = 1.



c     --- volume density ---

        do 10 k=1, ip

        fs = 42.9e-3
c       fs = 30.e-3
        epsig = epsig
        epsil = epsil
        pi = 4.*atan(1.)
        vp(k) = pi*dp(k)*dp(k)*dp(k)/6.
c       vp(k) = rhol*pi*dp(k)*dp(k)*dp(k)/6.

        epsif = ( 3.*epsig*(epsil-epsig)*eps0 ) / ( 2.*epsig+epsil )

c       vq(k) = 48.*pi*fs*abs(epsif)/vp(k)
c       vq(k) = sqrt( vq(k) )
c       vq(k) = vq(k) / rhol

        rhol = 1027.

        vq(k) = 8.*eps0*fs*dp(k)*dp(k)*dp(k)
        vq(k) = pi*sqrt( vq(k) )
        vq(k) = vq(k) / vp(k)
        vq(k) = vq(k) / rhol

c       write(*,*) vq(k)
 10     continue

c      write(*,*) vq(1)

c     --- particle position search ---

c      write(*,*) xp(k),yp(k)

       do 11 k=1, ip

c      fallback cell for particles at/beyond the right or top edge
c      (e.g. nozzle 2 particles on the mirror plane x=xlen)
       inex(k) = nx-2
       jnex(k) = ny-2

       do 12 i=1, nx-1

       if( x(i,1) .gt. xp(k) )then
         inex(k) = i-1
         go to 13
       end if

 12    continue
 13    continue

       do 14 j=1, ny-1

       if( y(1,j) .gt. yp(k) )then
         jnex(k) = j-1
         go to 15
       end if

 14    continue
 15    continue

 11    continue

c      write(*,*) inex(1), jnex(1)

c     --- electrostatic force interpolation ---
         
       do 20 k=1, ip

       l1(k) = sqrt( (xp(k)-x(inex(k),jnex(k)))
     &              *(xp(k)-x(inex(k),jnex(k)))
     &              +(yp(k)-y(inex(k),jnex(k)))
     &              *(yp(k)-y(inex(k),jnex(k))) )

       l2(k) = sqrt( (x(inex(k)+1,jnex(k))-xp(k))
     &              *(x(inex(k)+1,jnex(k))-xp(k))
     &              +(yp(k)-y(inex(k)+1,jnex(k)))
     &              *(yp(k)-y(inex(k)+1,jnex(k))) )

       l3(k) = sqrt( (x(inex(k)+1,jnex(k)+1)-xp(k))
     &              *(x(inex(k)+1,jnex(k)+1)-xp(k))
     &              +(y(inex(k)+1,jnex(k)+1)-yp(k))
     &              *(y(inex(k)+1,jnex(k)+1)-yp(k)) )

       l4(k) = sqrt( (xp(k)-x(inex(k),jnex(k)+1))
     &              *(xp(k)-x(inex(k),jnex(k)+1))
     &              +(y(inex(k),jnex(k)+1)-yp(k))
     &              *(y(inex(k),jnex(k)+1)-yp(k)) )

c      write(*,*) l1(1), l2(1), l3(1), l4(1)

       l1(k) = 1. / l1(k)
       l2(k) = 1. / l2(k)
       l3(k) = 1. / l3(k)
       l4(k) = 1. / l4(k)

c      write(*,*) inex(k),jnex(k)
c      write(*,*) l1(1), l2(1), l3(1), l4(1)

       ijk = jj(inex(k),jnex(k))

c      write(*,*) ijk,inex(k),jnex(k)


          Ernex(k) = ( Er(ijk)*l1(k)+Er(ijk+(ny-1))*l2(k)
     &             +   Er(ijk+ny)*l3(k)+Er(ijk+1)*l4(k) )
     &             / ( l1(k)+l2(k)+l3(k)+l4(k) )

          Eznex(k) = ( Ez(ijk)*l1(k)+Ez(ijk+(ny-1))*l2(k)
     &             +   Ez(ijk+ny)*l3(k)+Ez(ijk+1)*l4(k) )
     &             / ( l1(k)+l2(k)+l3(k)+l4(k) ) 

c      write(*,*) Ernex(1),Eznex(1),Er(ijk)

 20    continue


       do 30 k=1,ip

          fer(k) = vq(k)*Ernex(k)
          fez(k) = vq(k)*Eznex(k)

 30    continue

c      write(*,*) fer(1),fez(1)

c     --- next particle velocity ---
       do 40 k=1, ip 
  
c      if( pout(k) .ne. 0. )then
c      if( pout(k) .eq. 1. )then

c        uxp(k) = 0.
c        uyp(k) = 0.
c        go to 45

c      else

c        uxp(k) = uxp(k) + dt*fer(k)
c        uyp(k) = uyp(k) + dt*fez(k)    

c      end if

c45    continue
 40    continue

c      write(*,*) uxp(1),uyp(1)

c      do 45 k=1, ip

c      if( xp(k) .lt. 0. )then

c      uxp(k) = 0.
c      uyp(k) = uyp(k)

c      end if

c45    continue


c     --- next particle position ---
c      write(*,*) xp(1), yp(1)

c      xlen = 1.e-1
c      ylen = 1.e-1

       do 50 k=1, ip
 
c      if( pout(k) .ne. 0. )then
c      if( pout(k) .eq. 1. )then
 
c         xp(k) = xp(k)
c         yp(k) = yp(k)

c      else if( pout(k) .eq. 2. )then

c         xp(k) = xlen
c         yp(k) = 0.

c      else if( xp(k) + dt*uxp(k) .lt. 0. )then

c         xp(k) = 0.
c         xp(k) = xp(k)
c         yp(k) = yp(k) + dt*uyp(k)

c      else if( xp(k) .lt. 0. )then

c         xp(k) = 0.
c         yp(k) = yp(k) + dt*uyp(k)

c      else

c         xp(k) = xp(k) + dt*uxp(k)
c         yp(k) = yp(k) + dt*uyp(k)

c      end if

 50    continue


c      call pcount


      return

      end


       subroutine diffusion_gas
       include 'subcom.inc'

       real sp(100000), vp(100000)
       real inex(100000), jnex(100000)
       real l1(100000), l2(100000), l3(100000), l4(100000)
c      real cltmx(1000)


c      kk = 9.e-1
c      kk = 1.6

       rhol = 1027.
c      rhol = 1240.
       pi = 4.*atan(1.)


c      ------ particle position search ------

       do 10 k=1, ip

       do 11 i=1, nx-1

       if( x(i,1) .gt. xp(k) )then
         inex(k) = i-1
         go to 12
       end if

 11    continue
 12    continue

       do 13 j=1, ny-1

       if( y(1,j) .gt. yp(k) )then
         jnex(k) = j-1
         go to 14
       end if

 13    continue
 14    continue


 10    continue


c      ------ closest lattice search ------ 

       do 20 k=1, ip

       l1(k) = sqrt( (xp(k)-x(inex(k),jnex(k)))
     &              *(xp(k)-x(inex(k),jnex(k)))
     &              +(yp(k)-y(inex(k),jnex(k)))
     &              *(yp(k)-y(inex(k),jnex(k))) )

       l2(k) = sqrt( (x(inex(k)+1,jnex(k))-xp(k))
     &              *(x(inex(k)+1,jnex(k))-xp(k))
     &              +(yp(k)-y(inex(k)+1,jnex(k)))
     &              *(yp(k)-y(inex(k)+1,jnex(k))) )

       l3(K) = sqrt( (x(inex(k)+1,jnex(k)+1)-xp(k))
     &              *(x(inex(k)+1,jnex(k)+1)-xp(k))
     &              +(y(inex(k)+1,jnex(k)+1)-yp(k))
     &              *(y(inex(k)+1,jnex(k)+1)-yp(k)) )

       l4(k) = sqrt( (xp(k)-x(inex(k),jnex(k)+1))
     &              *(xp(k)-x(inex(k),jnex(k)+1))
     &              +(y(inex(k),jnex(k)+1)-yp(k))
     &              *(y(inex(k),jnex(k)+1)-yp(k)) )


c      ------ concentration interpolation ------

       sp(k) = pi*dp(k)*dp(k)
c      sp(k) = pi*dp(k)*dp(k)/4.
       vp(k) = pi*dp(k)*dp(k)*dp(k)/6.

       if( clt(k) .gt. cltmx(k) )then
c        write(*,*) k, clt(k), cltmx(k)
         go to 24
       end if

c      abco2 = 1.e-2
c      abco2 = 1.e-12
c      abco2 = 1.e-20
c      abco2 = 1.e-21
       abco2 = 1.e-25

       if( l1(k) .le. l2(k)  .and.  l1(k) .le. l4(k) )then
         cg(inex(k),jnex(k)) = cg(inex(k),jnex(k))
     &                       - dt*abco2*cg(inex(k),jnex(k))*sp(k)/vp(k)
         go to 22
       end if

       if( l2(k) .le. l3(k)  .and.  l2(k) .le. l3(k) )then
         cg(inex(k)+1,jnex(k)) = cg(inex(k)+1,jnex(k))
     &                      - dt*abco2*cg(inex(k)+1,jnex(k))*sp(k)/vp(k)
         go to 22
       end if

       if( l3(k) .lt. l4(k)  .and.  l3(k) .le. l2(k) )then
         cg(inex(k)+1,jnex(k)+1) = cg(inex(k)+1,jnex(k)+1)
     &                    - dt*abco2*cg(inex(k)+1,jnex(k)+1)*sp(k)/vp(k)
         go to 22
       end if

       if( l4(k) .lt. l3(k)  .and.  l4(k) .lt. l1(k) )then
         cg(inex(k),jnex(k)+1) = cg(inex(k),jnex(k)+1)
     &                      - dt*abco2*cg(inex(k),jnex(k)+1)*sp(k)/vp(k)
         go to 22
       end if


 22    continue

 24    continue

 20    continue

   
       call bound1


       return 


       end 



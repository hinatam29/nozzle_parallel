
       subroutine pcount
       include 'subcom.inc'

       real vp(100000),vq(100000)
       real fe(100000), fint(100000)

       save nbad
       data nbad /0/


       xhol = 5.e-3 
c      yhol1 = 96.e-3
       yhol1 = 93.e-3
c      yhol2 = 95.e-3
       yhol2 = 92.e-3
c      xlen = 1.e-1   旧: 10cm。geo.f と同じ値にすること
       xlen = 4.e-2
c      xlen = 50.e-3
       ylen = 1.e-1


       do 5 k=1, ip

       if( pout(k) .ne. 0. )then
       
       go to 6
       
       else

       rhol = 1027.
       pi = 4.*atan(1.)

       vp(k) = rhol*pi*dp(k)*dp(k)*dp(k)/6.

       eps0 = 8.8542e-12
       fs = 42.9e-3

       vq(k) = 8.*eps0*fs*dp(k)*dp(k)*dp(k)
       vq(k) = pi*sqrt( vq(k) )
       vq(k) = vq(k)/vp(k)
c      vq(k) = vq(k)/rhol


       gv = 9.8

       uxp(k) = uxp(k) - dt*fvisx(k)/vp(k)
     &        + dt*fer(k) 
     &        + dt*vq(k)*ftxt(k)/4./pi/eps0

       uyp(k) = uyp(k) - dt*fvisy(k)/vp(k)
     &        + dt*gv
     &        + dt*fez(k)
     &        + dt*vq(k)*ftyt(k)/4./pi/eps0

c      --- 発散粒子の検出・隔離（1粒子の発散が全粒子へ伝播するのを防ぐ）---
c      粘性/電場/クーロンのどの力が原因かを最大20回だけ出力する
       vmag = sqrt( uxp(k)*uxp(k) + uyp(k)*uyp(k) )
       if( vmag.ne.vmag .or. vmag.gt.1.e4 )then
          if( nbad.lt.20 )then
          write(*,*) 'BADP k=',k,' istp=',istp,
     &      ' x=',xp(k),' y=',yp(k),
     &      ' vis=',-dt*fvisx(k)/vp(k),
     &      ' fer=',dt*fer(k),
     &      ' coul=',dt*vq(k)*ftxt(k)/4./pi/eps0
          end if
          nbad = nbad + 1
          uxp(k) = 0.
          uyp(k) = 0.
          xp(k)  = 0.
          yp(k)  = 0.
          pout(k) = 9.
       end if


c      fe(k) = sqrt( fer(k)*fer(k) + fez(k)*fez(k) )
c      fint(k) = vq(k)*sqrt(ftxt(k)*ftxt(k)+ftyt(k)*ftyt(k))/4./pi/eps0

c           if( fint(k) .eq. 0. )then
c            fcom(k) = 0.
c           else
c            fcom(k) = fe(k) / fint(k)
c           end if

       end if

 6     continue
 5     continue

c      write(*,*) xp(1),yp(1)
c      write(*,*) fvisx(1),fvisy(1)
c      write(*,*) fer(1),fez(1)
c      write(*,*) ftxt(1),ftyt(1)

       do 7 k=1, ip

c      --- reflect at left axis (x=0) ---
       if( xp(k) + dt*uxp(k) .lt. 0. )then
       uxp(k) = -uxp(k)
       end if

c      --- reflect at right mirror plane (x=xlen, nozzle 2) ---
       if( xp(k) + dt*uxp(k) .gt. xlen )then
       uxp(k) = -uxp(k)
       end if

 7     continue


       do 8 k=1, ip

       if( pout(k) .ne. 0. )then

       go to 9

       else

c           if( xp(k) + dt*uxp(k) .lt. 0. )then
       
c           xp(k) = xp(k)
c           yp(k) = yp(k) + dt*uyp(k)

c           else 

            xp(k) = xp(k) + dt*uxp(k)
            yp(k) = yp(k) + dt*uyp(k)

c           end if

       end if

 9     continue
 8     continue       


       do 10 k=1, ip


c      ------ counter electrode ------

c      if( xp(k) .ge. xhol  .and.  yp(k) .ge. yhol2
c    &     .and.  yp(k) .le. yhol1 )then

       if( i .ge. ihol  .and.  j .ge. jhol2  .and.  j .le. jhol1 )then

       pout(k) = 1.
c      write(*,*) k

c      ------ out of area ------

       else if( xp(k) .ge. 0.  .and.  yp(k) .ge. ylen )then  

       pout(k) = 2.
       xp(k) = xp(k)
       yp(k) = ylen

       else if( yp(k) .ge. 0.  .and.  xp(k) .ge. xlen )then

c      --- right side is a mirror plane (nozzle 2): keep in domain ---
       pout(k) = 0.
       xp(k) = xlen
       yp(k) = yp(k)

       else if( xp(k) .ge. 0.  .and.  yp(k) .le. 0. )then

       pout(k) = 4.
       xp(k) = xp(k)
       yp(k) = 0.

c      ------ axial ------

c      else if( xp(k) .le. 1.e-7 ) then
c      else if( xp(k) .lt. 0. ) then

c      pout(k) = 1.

       else if( xp(k) .lt. 0. )then

       pout(k) = 0. 
       xp(k) = 0.
       yp(k) = yp(k)

       else
    
       pout(k) = 0.

       end if

 10    continue



c      ------ out of area ------

       do 20 k=1,ip

c      if( xp(k) .lt. 0. )then
c      write(*,*) k, dp(k), xp(k), yp(k)
c      end if

c      if( yp(k) .le. 0. )then
c      pout(k) = 1.
c      end if

c      if( xp(k) .ge. xlen )then
c      pout(k) = 1.
c      end if

 20    continue


       return 

       end




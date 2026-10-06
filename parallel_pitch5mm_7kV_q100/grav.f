
      subroutine grav
      include 'subcom.inc'

c     do 100 i=2, nx-1
c     do 110 j=2, ny-2
c       uy_sx(i,j) = uy_sx(i,j)-dt*gv
c110   continue
c100   continue

c     do 200 i=2, nx-2
c     do 210 j=2, ny-1
c       uy_sy(i,j) = uy_sy(i,j)-dt*gv
c210   continue
c200   continue

c     do 300 i=2, nx-2
c     do 300 j=2, ny-2
c       uy0(i,j) = uy0(i,j)-dt*gv
c       uy0p(i,j) = uy0(i,j)-dt*gv
c       uy0p(i,j) = uy0(i,j)-dt*gv
c       uyp = uyp-dt*gv
c300   continue

      do 10 k=1,ip

        if( pout(k) .ne. 0. )then
c       if( pout(k) .eq. 1. )then

        uxp(k) = 0.
        uyp(k) = 0.
c       yp(k) = yp(k)

        else

        uxp(k) = uxp(k)
        uyp(k) = uyp(k) - dt*gv
c       yp(k) = yp(k) + dt*uyp(k)

        end if

 10   continue


      xlen = 1.e-1
      ylen = 1.e-1

      do 20 k=1, ip

        if( pout(k) .ne. 0. )then
c       if( pout(k) .eq. 1. )then
        
        xp(k) = xp(k)
        yp(k) = yp(k) 

c       else if( pout(k) .eq. 2. )then

c       xp(k) = xlen
c       yp(k) = 0.

c       else if( xp(k) .lt. 0. )then

c       xp(k) = 0.
c       yp(k) = yp(k)

        else

        xp(k) = xp(k)
        yp(k) = yp(k) + dt*uyp(k)

        end if

 20   continue


c     call pcount

c       write(*,*) xp(1),yp(1),uxp(1),uyp(1)

c       uyp2 = uyp2-dt*gv
c       yp2 = yp2+dt*uyp2


c       write(*,*) yp1, yp2

c     call result(uy_sx,uy_sy,uy0,uypx,uypy,uy0p)

c     call bound1

      return
      end

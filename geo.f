
      subroutine geo
      include 'subcom.inc'

c     xlen = 0.165e-3 
      xlen = 1.e-1 
c     xlen = 50.e-3 
c     xlen = 20.e-3 
c     ylen = 1.43e-3
c     ylen = 0.53e-3
      ylen = 1.e-1
      dx = xlen/real(nx)
      dy = ylen/real(ny)

      xnoz = 0.4e-3
c     xnoz = 1.e-3
      ynoz = 1.e-3
c     ynoz = 1.5e-3
      xhol = 5.e-3
c     yhol1 = 96.e-3
      yhol1 = 93.e-3
c     yhol2 = 95.e-3
      yhol2 = 92.e-3

      do 50 j=1, ny-1
      do 50 i=1, nx-1
        ijk = (i-1)*(ny-1)+j
        jj(i,j) = ijk
50    continue

c     do 55 i=1, nx
c     do 55 j=1, ny
c       jj(i,j) = ijk
c       write(*,*) i,j, ijk
c55    continue

c     open(10,file='./grid79.dat',status='old')

      do 100 j=1, ny
      do 100 i=1, nx
c     do 100 j=1, ny
c       read(10,*) x(i,j), y(i,j)
c       x(i,j) = 0.5*x(i,j)
c       y(i,j) = 0.37*y(i,j)
        x(i,j) = 1.e-15 + dx*(i-1)
        y(i,j) = 0. + dy*(j-1)
 100   continue

      do 120 j=1, ny
        x(0,j) = x(1,j) - (x(2,j)-x(1,j))
        y(0,j) = y(1,j) - (y(2,j)-y(1,j))
        x(nx+1,j) = x(nx,j) + (x(nx,j)-x(nx-1,j))
        y(nx+1,j) = y(nx,j) + (y(nx,j)-y(nx-1,j))
120   continue

      do 130 i=1, nx
        x(i,0) = x(i,1) - (x(i,2)-x(i,1))
        y(i,0) = y(i,1) - (y(i,2)-y(i,1))
        x(i,ny+1) = x(i,ny) + (x(i,ny)-x(i,ny-1))
        y(i,ny+1) = y(i,ny) + (y(i,ny)-y(i,ny-1))
130   continue

      do 300 i=0, nx
      do 300 j=0, ny+1
        sdx(i,j) = x(i+1,j)-x(i,j) 
c       sdx(i,j) = x(i,j)-x(i-1,j) 
c       sdx(i,j) = x(1,j)-x(0,j) 
c       write(*,*) x(i+1,j)-x(i,j)-sdx(i,j)
 300   continue

      do 305 i=1, nx
      do 305 j=1, ny-1
        sdx1(i,j) = 0.5*(x(i+1,j)+x(i-1,j))
 305   continue

      do 400 i=0, nx+1
      do 400 j=0, ny
        sdy(i,j) = y(i,j+1)-y(i,j)
c       sdy(i,j) = y(i,2)-y(i,1) 
c       write(*,*) y(i,j+1)-y(i,j), sdy(i,j)
400   continue

      do 405 i=1, nx-1
      do 405 j=1, ny
       sdy1(i,j) = 0.5*(y(i,j+1)+y(i,j-1))
405   continue

      do 500 j=0, ny+1
       sdx(nx+1,j) = x(nx+1,j)-x(nx,j)
500   continue

      do 600 i=0, nx+1
       sdy(i,ny+1) = y(i,ny+1)-y(i,ny)
600   continue

      do 700 i=1, nx
        if( x(i,1) .lt. xnoz ) then
            inoz = i-1
        end if
 700   continue

      do 800 j=1, ny
        if( y(1,j) .lt. ylen-ynoz) then
            jnoz = j
         end if
 800   continue

      do 900 i=1, nx
        if( x(i,1) .lt. xhol )then
            ihol = i-1
        end if
 900   continue

      do 1000 j=1, ny
        if( y(nx,j) .lt. yhol1 )then
            jhol1 = j-1
        end if
 1000  continue

      do 1100 j=1, ny
        if( y(nx,j) .lt. yhol2 )then
            jhol2 = j
        end if
 1100  continue


c     --- 2nd nozzle (mirror image at right boundary, x=xlen) ---
c     right nozzle inner edge : mirror of left nozzle (i=1..inoz)
c     about the right boundary i=nx-1  ->  i=inoz2..nx-1
      inoz2 = nx - inoz
c     right hole of counter electrode : mirror of left hole (i<ihol)
      ihol2 = nx - ihol

c     write(*,*) inoz,jnoz,ihol,jhol1,jhol2

      write(*,*) 'nozzle1: ', x(inoz,ny)
      write(*,*) 'nozzle2 (mirror): ', x(inoz2,ny)
      write(*,*) 'counter electrode: ', x(ihol,jhol1)
      write(*,*) 'electrode right hole: ', x(ihol2,jhol1)

      return
      end


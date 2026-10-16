
      subroutine bound1
      include 'subcom.inc'


c     --- right side ---

      do 100  j=1, ny

      cg(nx,j) = -cg(nx-1,j)

c     ux0(nx-1,j) = ux0(nx-2,j)
c     ux_sx(nx-1,j) = ux_sx(nx-2,j)
c     ux_sy(nx-1,j) = ux_sy(nx-2,j)

 100  continue      

c     --- left side ---

      do 200 j=1, ny-1

c     cg(1,j) = cg(2,j)

c     ux0(1,j) = ux0(2,j)
c     ux_sx(1,j) = ux_sx(2,j)
c     ux_sy(1,j) = ux_sy(2,j)

 200  continue

c     --- top ---

      do 300 i=1, nx-1

      cg(i,ny) = -cg(i,ny-1)

c     ux0(i,ny-1) = ux0(i,ny-2)
c     ux_sx(i,ny-1) = ux_sx(i,ny-2)
c     ux_sy(i,ny-1) = ux_sy(i,ny-2)

 300  continue

c     --- bottom ---

      do 400 i=1, nx-1

      cg(i,0) = -cg(i,1)

c     ux0(i,1) = ux0(i,2)
c     ux_sx(i,1) = ux_sx(i,2)
c     ux_sy(i,1) = ux_sy(i,2)

 400  continue

c     --- side nozzle ---

      do 500 j=jnoz, ny-1

      cg(inoz,j) = cg(inoz+1,j)

 500  continue

c     --- bottom nozzle ---

      do 600 i=1, inoz-1

      cg(i,jnoz) = cg(i,jnoz-1)

 600  continue

c     --- inside nozzle ---

      do 650 i=1, inoz-1
      do 650 j=jnoz+1, ny-1

c     cg(i,j) = 0.

 650  continue

c     --- top counter nozzle ---

      do 700 i=ihol+1, nx-1

      cg(i,jhol1) = cg(i,jhol1+1)

 700  continue

c     --- bottom counter nozzle ---

      do 800 i=ihol+1, nx-1

      cg(i,jhol2) = cg(i,jhol2-1)

 800  continue

c     --- side counter nozzle ---

      do 900 j=jhol2, jhol1

      cg(ihol,j) = cg(ihol-1,j)

 900  continue

c     --- inside counter nozzle --- 

      do 950 i=ihol+1, nx-1
      do 950 j=jhol2+1, jhol1-1

      cg(i,j) = 0.

 950  continue


      end




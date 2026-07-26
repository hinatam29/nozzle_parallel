
      subroutine dfunc
      include 'subcom.inc'

      dimension isrfx(nx*ny), isrfy(nx*ny)
      dimension dphi0(nx,ny), dval(nx,ny), val(nx,ny)
      dimension xsrf(nx*ny), ysrf(nx*ny)
      dimension dphis(nx,ny)

      k = 1

      isrf = 11

      do 10 i=2, nx-2
      do 10 j=2, ny-2
        ie = i+1
        iw = i
        jn = j+1
        js = j


c       if( vf0(iw,j) .le. 0.9 .and. vf0(ie,j) .ge. 0.9
c    &  .or. vf0(iw,j) .ge. 0.9 .and. vf0(ie,j) .le. 0.9 ) then
        if( vf0(iw,j) .le. 0.5 .and. vf0(ie,j) .ge. 0.5
     &  .or. vf0(iw,j) .ge. 0.5 .and. vf0(ie,j) .le. 0.5 ) then
             xx = 0.5*(x(i,j)+x(i+1,j))
             yy = 0.5*(y(i,j)+y(i,j+1))
             xsrf(k) = sdx(i,j)/(vf0(ie,j)-vf0(iw,j))*(0.5-vf0(iw,j))+xx
             ysrf(k) = yy
           k = k+1
        end if


        ii = i
        jp = j+1

c       if( vf0(i,jn) .le. 0.9 .and. vf0(i,js) .ge. 0.9
c    &  .or. vf0(i,jn) .ge. 0.9 .and. vf0(i,js) .le. 0.9 ) then
        if( vf0(i,jn) .le. 0.5 .and. vf0(i,js) .ge. 0.5
     &  .or. vf0(i,jn) .ge. 0.5 .and. vf0(i,js) .le. 0.5 ) then
             xx = 0.5*(x(i,j)+x(i+1,j))
             yy = 0.5*(y(i,j)+y(i,j+1))
             xsrf(k) = xx
             ysrf(k) = sdy(i,j)/(vf0(i,jn)-vf0(i,js))*(0.5-vf0(i,js))+yy
           k = k+1
        end if
 10    continue


c:    --------------------------------------
c:    calculating distance function phi
c:    --------------------------------------

      do 20 i=2, nx-2
      do 20 j=2, ny-2
        if( vf0(i,j) .gt. 0.5 ) then
c       if( vf0(i,j) .gt. 0.9 ) then
           sgn = 1.
        else if( vf0(i,j) .le. 0.5 ) then
c       else if( vf0(i,j) .le. 0.9 ) then
           sgn = -1.
        end if

        dphimin = 10000.   

        do 30 ik=1, k-1
        xx = 0.5*(x(i,j)+x(i+1,j))
        yy = 0.5*(y(i,j)+y(i,j+1))
        dis = sqrt((xx-xsrf(ik))**2.+(yy-ysrf(ik))**2.) 

        if( dis .le. dphimin ) dphimin = dis
30      continue

        dphi(i,j) = sgn*dphimin
20    continue

      do 33 j=1, ny
        dphi(1,j) = dphi(2,j)
        dphi(nx-1,j) = dphi(nx-2,j)
33    continue

      do 36 i=1, nx
        dphi(i,1) = dphi(i,2)
        dphi(i,ny-1) = dphi(i,ny-2)
36    continue

      go to 60
     
      do 42 i=2, nx-2
      do 42 j=2, ny-2

        eps0 = sdx(i,j)*sdx(i,j) + sdy(i,j)*sdy(i,j)
        dval(i,j) = dphi(i,j)/sqrt(dphi(i,j)*dphi(i,j)+eps0)
        dphi0(i,j) = dphi(i,j)
42    continue


45    continue
      sumd = 0.
      inum = 1 

      do 50 i=2, nx-2
      do 50 j=2, ny-2
        xx1 = 0.5*(x(i,j) + x(i-1,j))
        xx2 = 0.5*(x(i+1,j) + x(i,j))
        xx3 = 0.5*(x(i+2,j) + x(i+1,j))
        yy1 = 0.5*(y(i,j) + y(i,j-1))
        yy2 = 0.5*(y(i,j+1) + y(i,j))
        yy3 = 0.5*(y(i,j+2) + y(i,j+1))

        if( dphi0(i,j) .gt. 0. ) then
        ga = (dphi(i,j)-dphi(i-1,j))/(xx2 - xx1)
        gb = -(dphi(i+1,j)-dphi(i,j))/(xx3 - xx2)
        gc = (dphi(i,j)-dphi(i,j-1))/(yy2 - yy1)
        gd = -(dphi(i,j+1)-dphi(i,j))/(yy3 - yy2)
           val(i,j) = sqrt( max(ga*ga,gb*gb)+max(gc*gc,gd*gd) )-1. 
        else 
           val(i,j) = 0.
        end if

50    continue

c     sdx  = x(i+1,j) - x(i,j)
      dtt = 0.1*sdx(i,j)

      do 55 i=2, nx-1 
      do 55 j=2, ny-1
        sdx  = x(i+1,j) - x(i,j)

        dphip = dphi(i,j)
        dphi(i,j) = dphi(i,j) - dtt*dval(i,j)*val(i,j)
        if( abs(dphi(i,j)) .le. sdx(i,j) ) then
        sumd = sumd + abs(dphi(i,j)-dphip)
        inum = inum + 1
        end if
55    continue

      do 57 ik=1, k-1
        ii = isrfx(ik)
        jp = isrfy(ik)
        dphi(ii,jp) = 0.
57    continue


c     sdx  = x(i+1,j) - x(i,j)
      if( sumd/(inum-1) .gt. sdx(i,j)*sdx(i,j)*dtt ) go to 45

 
60    continue


5000  return
      end


      subroutine diffusion
      include 'subcom.inc'

      do 100 i=2, nx-2
      do 100 j=2, ny-2

c     if( i .ge. 2 .and. i .le. iwall .and. j .ge. jnoz ) go to 100

        if( trho(i,j) .le. 0. ) write(*,*) '0!!! ', trho(i,j), i, j
        xmux0 = 0.5*(tmu(i,j) + tmu(i-1,j))
        xmux1 = 0.5*(tmu(i,j) + tmu(i+1,j))
        xmuy0 = 0.5*(tmu(i,j) + tmu(i,j-1))
        xmuy1 = 0.5*(tmu(i,j) + tmu(i,j+1))

        sdx = x(i+1,j) - x(i,j)
        sdy = y(i,j+1) - y(i,j)
 
        xx1 = 0.5*(x(i,j) + x(i-1,j))
        xx2 = 0.5*(x(i+1,j) + x(i,j))
        xx3 = 0.5*(x(i+2,j) + x(i+1,j))
        yy1 = 0.5*(y(i,j) + y(i,j-1))
        yy2 = 0.5*(y(i,j+1) + y(i,j))
        yy3 = 0.5*(y(i,j+2) + y(i,j+1))

        uxx_x0 = (ux0(i,j)-ux0(i-1,j))/(xx2-xx1) 
        uxx_x1 = (ux0(i+1,j)-ux0(i,j))/(xx3-xx2)

        uxy_y0 = (ux0(i,j)-ux0(i,j-1))/(yy2-yy1)
        uxy_y1 = (ux0(i,j+1)-ux0(i,j))/(yy3-yy2) 


        uyx_x0 = (uy0(i,j)-uy0(i-1,j))/(xx2-xx1) 
        uyx_x1 = (uy0(i+1,j)-uy0(i,j))/(xx3-xx2) 

        uyy_y0 = (uy0(i,j)-uy0(i,j-1))/(yy2-yy1) 
        uyy_y1 = (uy0(i,j+1)-uy0(i,j))/(yy3-yy2) 
       
       
        t11_x0 = xmux0*uxx_x0
        t12_x0 = xmux0*uyx_x0
        t12_y0 = xmuy0*uxy_y0
        t22_y0 = xmuy0*uyy_y0 

        t11_x1 = xmux1*uxx_x1 
        t12_x1 = xmux1*uyx_x1
        t12_y1 = xmuy1*uxy_y1
        t22_y1 = xmuy1*uyy_y1 
       
        viscr = 1./trho(i,j)*( (t11_x1-t11_x0)/sdx(i,j)
     &                            + (t12_y1-t12_y0)/sdy(i,j) )
     &             + 1./trho(i,j)*1.0/xx2
     &             * ( vf0(i,j)*xmul + (1.-vf0(i,j))*xmug )
     &             * ( (ux0(i+1,j)-ux0(i-1,j))/(xx3 - xx1) 
     &               - ux0(i,j)/xx2 )

        viscz = 1./trho(i,j)*( (t12_x1-t12_x0)/sdx(i,j)
     &                            + (t22_y1-t22_y0)/sdy(i,j) )
     &             + 1./trho(i,j)*1.0/xx2
     &             * ( vf0(i,j)*xmul + (1.-vf0(i,j))*xmug )
     &             * (uy0(i+1,j)-uy0(i-1,j))/(xx3 - xx1)

        ux0p(i,j) = ux0(i,j) + dt * viscr
        uy0p(i,j) = uy0(i,j) + dt * viscz

100   continue


      call TECSIAx(ux_sx,ux0,uxpx,ux0p)
      call TECSIAx(uy_sx,uy0,uypx,uy0p)
      call TECSIAy(ux_sy,ux0,uxpy,ux0p)
      call TECSIAy(uy_sy,uy0,uypy,uy0p)


      call result(ux_sx,ux_sy,ux0,uxpx,uxpy,ux0p)
      call result(uy_sx,uy_sy,uy0,uypx,uypy,uy0p)

      call bound

      return
      end
        
      

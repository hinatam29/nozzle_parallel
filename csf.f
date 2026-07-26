
      subroutine csf
      include 'subcom.inc'


      dimension grad1x(nx,ny), grad1y(nx,ny)
c     dimension gradx(nx,ny), grady(nx,ny)
      dimension gradx2(nx,ny), grady2(nx,ny), gradxy(nx,ny)
      dimension vcnx3(nx,ny)
      dimension vcny3(nx,ny)
      dimension isrfx(nx*ny), isrfy(nx*ny)
      dimension dphi0(nx,ny), dval(nx,ny), val(nx,ny)
      dimension xsrf(nx*ny), ysrf(nx*ny)

      real nomfct

      sgm_sf = 23.61e-3

      do 100 i=1, nx-3
      do 100 j=1, ny-3
         xx2 = 0.5*(x(i+1,j) + x(i,j))
         xx3 = 0.5*(x(i+2,j) + x(i+1,j))
         yy2 = 0.5*(y(i,j+1) + y(i,j))
         yy3 = 0.5*(y(i,j+2) + y(i,j+1))

         xx3 = x(i+1,j)
         xx2 = x(i,j)
         yy3 = y(i,j+1)
         yy2 = y(i,j)


         grad1x(i,j) = 0.5*(dphi(i+1,j)-dphi(i,j)
     &                     +dphi(i+1,j+1)-dphi(i,j+1))/(xx3-xx2)
         grad1y(i,j) = 0.5*(dphi(i,j+1)-dphi(i,j)
     &                     +dphi(i+1,j+1)-dphi(i+1,j))/(yy3-yy2)

100   continue 

      do 105 j=1, ny-3
         xx2 = 0.5*(x(nx-1,j) + x(nx-2,j))
         xx3 = 0.5*(x(nx,j) + x(nx-1,j))
         yy2 = 0.5*(y(nx,j+1) + y(nx,j))
         yy3 = 0.5*(y(nx,j+2) + y(nx,j+1))

         xx3 = x(nx-1,j)
         xx2 = x(nx-2,j)
         yy3 = y(nx,j+2)
         yy2 = y(nx,j+1)
 

         grad1x(nx-1,j) = 0.5*(dphi(nx-1,j)-dphi(nx-2,j)
     &                     +dphi(nx-1,j+1)-dphi(nx-2,j+1))/(xx3-xx2)
         grad1y(nx-1,j) = 0.5*(dphi(nx-2,j+1)-dphi(nx-2,j)
     &                     +dphi(nx-1,j+1)-dphi(nx-1,j))/(yy3-yy2)

105   continue 

      do 110 i=1, nx-3
         xx2 = 0.5*(x(i+1,ny) + x(i,ny))
         xx3 = 0.5*(x(i+2,ny) + x(i+1,ny))
         yy2 = 0.5*(y(i,ny-1) + y(i,ny-2))
         yy3 = 0.5*(y(i,ny) + y(i,ny-1))

         xx3 = x(i+1,ny)
         xx2 = x(i+2,ny)
         yy3 = y(i,ny-1)
         yy2 = y(i,ny-2)

         grad1x(i,ny-1) = 0.5*(dphi(i+1,ny-2)-dphi(i,ny-2)
     &                     +dphi(i+1,ny-1)-dphi(i,ny-1))/(xx3-xx2)
         grad1y(i,ny-1) = 0.5*(dphi(i,ny-1)-dphi(i,ny-2)
     &                     +dphi(i+1,ny-1)-dphi(i+1,ny-2))/(yy3-yy2)
110   continue 

      do 120 i=2, nx-2
      do 120 j=2, ny-2
         gradx(i,j) = 0.25*(grad1x(i,j)+grad1x(i-1,j)
     &                     +grad1x(i,j-1)+grad1x(i-1,j-1))
         grady(i,j) = 0.25*(grad1y(i,j)+grad1y(i-1,j)
     &                     +grad1y(i,j-1)+grad1y(i-1,j-1))
120   continue

      do 130 i=2, nx-2
      do 130 j=2, ny-2
         xx1 = 0.5*(x(i,j) + x(i-1,j))
         xx3 = 0.5*(x(i+2,j) + x(i+1,j))
         yy1 = 0.5*(y(i,j) + y(i,j-1))
         yy3 = 0.5*(y(i,j+2) + y(i,j+1))

         xx3 = x(i+1,j)
         xx1 = x(i-1,j)
         yy3 = y(i,j+1)
         yy1 = y(i,j-1)

         gradxy(i,j) = (dphi(i+1,j+1)+dphi(i-1,j-1)
     &                 -dphi(i+1,j-1)-dphi(i-1,j+1))
     &                 /(xx3 - xx1)/(yy3 - yy1)
130   continue

      do 140 i=2, nx-2
      do 140 j=2, ny-2
c        sdx  = x(i+1,j) - x(i,j)

         xx1 = 0.5*(x(i,j) + x(i-1,j))
         xx2 = 0.5*(x(i+1,j) + x(i,j))
         xx3 = 0.5*(x(i+2,j) + x(i+1,j))

         xx1 = x(i-1,j)
         xx2 = x(i,j)
         xx3 = x(i+1,j)

         dphi_x  = ( dphi(i+1,j) - dphi(i,j) )/(xx3-xx2)
         dphi_x1 = ( dphi(i,j) - dphi(i-1,j) )/(xx2-xx1)

         gradx2(i,j) = ( dphi_x - dphi_x1 )/sdx(i,j)
140   continue

      do 150 i=2, nx-2
      do 150 j=2, ny-2
c        sdy  = y(i,j+1) - y(i,j)

         yy1 = 0.5*(y(i,j) + y(i,j-1))
         yy2 = 0.5*(y(i,j+1) + y(i,j))
         yy3 = 0.5*(y(i,j+2) + y(i,j+1))

         yy1 = y(i,j-1)
         yy2 = y(i,j)
         yy3 = y(i,j+1)


         dphi_y  = ( dphi(i,j+1) - dphi(i,j) )/(yy3-yy2)
         dphi_y1 = ( dphi(i,j) - dphi(i,j-1) )/(yy2-yy1)

         grady2(i,j) = ( dphi_y - dphi_y1 )/sdy(i,j)
150   continue


c     cntr = 0.
c     fct0 = 0.25/eps0*(1.+cos(pi*cntr/eps0))
c    &       * (sign(1.,cntr/eps0+1.)-sign(1.,cntr/eps0-1.))

      do 200 i=2, nx-2
      do 200 j=2, ny-2
c     sdx  = x(i+1,j) - x(i,j)
c     sdy  = y(i,j+1) - y(i,j)

      alp = 1.e-5

c     if( i .ge. iinj .and. j .ge. jnoz ) then
c        fsvx(i,j) = 0.
c        fsvy(i,j) = 0.
c     else

c     if( vf0(i,j) .lt. 0.5 ) then
      if( abs(dphi(i,j)) .gt. 1.e-5 ) then
c     if( abs(dphi(i,j)) .gt. 1.e-5 .and. j .ge. jnoz-1 ) then
c     if( vf0(i,j) .eq. 0. ) then
         fsvx(i,j) = 0.
         fsvy(i,j) = 0.
      else

      if( gradx(i,j) .eq. 0. ) then
         fsvx(i,j) = 0.
      else if( grady(i,j) .eq. 0. ) then
         fsvy(i,j) = 0.
      else



         crv1 = gradx2(i,j)*grady(i,j)*grady(i,j)
     &       - 2.*gradx(i,j)*grady(i,j)*gradxy(i,j)
     &       + grady2(i,j)*gradx(i,j)*gradx(i,j) 

c        crv2 = (gradx(i,j)*gradx(i,j)+grady(i,j)*grady(i,j))**1.5
         crv2 = (1.+gradx(i,j)*gradx(i,j)+grady(i,j)*grady(i,j))**1.5
        
c        crv1 = (1.-2.*gradx(i,j)*gradx(i,j)+grady(i,j)*grady(i,j))
c    &          *gradx2(i,j)
c    &          -6.*gradx(i,j)*grady(i,j)*gradxy(i,j)
c    &         +(1.+gradx(i,j)*gradx(i,j)-2.*grady(i,j)*grady(i,j))
c    &          *grady2(i,j) 

c        crv2 = (1.+gradx(i,j)*gradx(i,j)+grady(i,j)*grady(i,j))**2.5
        
c        crv(i,j) = -crv1/crv2
         crv(i,j) = crv1/crv2


c        fct = 0.5*(1./alp + dphi(i,j)/alp*cos(pi*dphi(i,j)/alp))
         fct = 0.5*(1./alp + 1./alp*cos(pi*dphi(i,j)/alp))
         rhave = 0.5*(rhog+rhol)
         gval = sqrt(gradx(i,j)*gradx(i,j)+grady(i,j)*grady(i,j))
c        fsvx(i,j) = trho(i,j)/rhave*sgm_sf*crv*fct
         fsvx(i,j) = -trho(i,j)/rhave*sgm_sf*crv(i,j)*fct
c        fsvx(i,j) = -trho(i,j)/rhave*sgm_sf*crv(i,j)
     &                 *gradx(i,j)/gval
c        fsvy(i,j) = trho(i,j)/rhave*sgm_sf*crv*fct
         fsvy(i,j) = -trho(i,j)/rhave*sgm_sf*crv(i,j)*fct
c        fsvy(i,j) = -trho(i,j)/rhave*sgm_sf*crv(i,j)
     &                 *grady(i,j)/gval
         ux0p(i,j) = ux0(i,j)+dt/trho(i,j)*fsvx(i,j)
         uy0p(i,j) = uy0(i,j)+dt/trho(i,j)*fsvy(i,j)



         end if
      end if
c     end if
c     end if

200   continue
 
c      go to 5000

      call TECSIAx(ux_sx,ux0,uxpx,ux0p)
      call TECSIAx(uy_sx,uy0,uypx,uy0p)
      call TECSIAy(ux_sy,ux0,uxpy,ux0p)
      call TECSIAy(uy_sy,uy0,uypy,uy0p)

      call result(ux_sx,ux_sy,ux0,uxpx,uxpy,ux0p)
      call result(uy_sx,uy_sy,uy0,uypx,uypy,uy0p)

c     call boundu(ux_sx,ux_sy,ux0)
c     call boundv(uy_sx,uy_sy,uy0)
      call bound


5000  return
      end

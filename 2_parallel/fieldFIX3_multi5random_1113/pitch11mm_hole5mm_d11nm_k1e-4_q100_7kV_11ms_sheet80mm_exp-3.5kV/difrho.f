
      subroutine difrho
      include 'subcom.inc'

      Di = 4.0e-10

      do 100 i=2, nx-2
      do 100 j=2, ny-2
        rdifr = (rhoc0(i+1,j)-2.*rhoc0(i,j)
     &                       +rhoc0(i-1,j))/sdx(i,j)/sdx(i,j)
     &        + (rhoc0(i+1,j)-rhoc0(i-1,j))/sdx(i,j)/x(i,j)
        rdifz = (rhoc0(i,j+1)-2.*rhoc0(i,j)
     &                       +rhoc0(i,j-1))/sdy(i,j)/sdy(i,j)
  
        rhoc0(i,j) = rhoc0(i,j) + dt*Di*(rdifr+rdifz)
100   continue

c     call TECSIAx(rhoc0_sx,rhoc0,rhocpx,rhoc0p)
c     call TECSIAy(rhoc0_sy,rhoc0,rhocpy,rhoc0p)

c     call result(rhoc0_sx,rhoc0_sy,rhoc0,rhocpx,rhocpy,rhoc0p)

c     call bound      

      return
      end 


        

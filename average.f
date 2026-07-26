
      subroutine average
      include 'subcom.inc'

      dimension dphis(nx,ny)

      C1 = 1./(1.+1./sqrt(2.))
      C2 = C1/sqrt(2.)
      iind = 1.

      do 5 i=1, nx 
      do 5 j=1, ny
        dphis(i,j) = 0.
5     continue

10    continue

      
      do 20 i=2, nx-2
      do 20 j=2, ny-2

      eps0 = 10*sdx(i,j)

c     if( abs(dphi(i,j)) .le. eps0 ) then
      dphis(i,j) = 0.5*dphi(i,j)
     &            + 0.5/(1.+4.*C1+4.*C2)*(dphi(i,j)
     &            + C1*(dphi(i-1,j)+dphi(i+1,j)+dphi(i,j-1)+dphi(i,j+1))
     &            + C2*(dphi(i-1,j-1)+dphi(i-1,j+1)
     &                  +dphi(i+1,j-1)+dphi(i+1,j+1)))
c     end if

20    continue

      do 30 i=2, nx-2
      do 30 j=2, ny-2

      dphi(i,j) = dphis(i,j)

30    continue 

      iind = iind + 1.

      if ( iind .le. 10) go to 10

c     call TECSIAx(vf_sx,vf0,vfpx,vf0p)
c     call TECSIAy(vf_sy,vf0,vfpy,vf0p)
c     call result(vf_sx,vf_sy,vf0,vfpx,vfpy,vf0p)

c     call bound


      return
      end 

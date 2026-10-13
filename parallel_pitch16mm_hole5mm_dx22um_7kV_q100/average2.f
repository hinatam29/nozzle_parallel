
      subroutine average2
      include 'subcom.inc'

      dimension vfs(nx,ny)

      C1 = 1./(1.+1./sqrt(2.))
      C2 = C1/sqrt(2.)
      iind = 1.

10    continue
      
      do 20 i=2, nx-2
      do 20 j=2, ny-2


c     if( abs(vf0(i,j)) .le. eps0 ) then
      vfs(i,j) = 0.5*vf0(i,j)
     &            + 0.5/(1.+4.*C1+4.*C2)*(vf0(i,j)
     &            + C1*(vf0(i-1,j)+vf0(i+1,j)+vf0(i,j-1)+vf0(i,j+1))
     &            + C2*(vf0(i-1,j-1)+vf0(i-1,j+1)
     &                  +vf0(i+1,j-1)+vf0(i+1,j+1)))
c     end if

20    continue

      do 30 i=2, nx-2
      do 30 j=2, ny-2

      vf0(i,j) = vfs(i,j)

30    continue 

      iind = iind + 1.

      if ( iind .le. 6) go to 10

      call TECSIAx(vf_sx,vf0,vfpx,vf0p)
      call TECSIAy(vf_sy,vf0,vfpy,vf0p)
      call result(vf_sx,vf_sy,vf0,vfpx,vfpy,vf0p)

      call bound

      return
      end 

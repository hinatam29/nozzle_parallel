
      subroutine fixvf0
      include 'subcom.inc'

      do 100 i=1,nx-1
      do 100 j=1,ny-1

      if( dphi(i,j) .gt. 0. ) then
          vf0(i,j) = 1.0
          trho(i,j) = rhol
          tmu(i,j) = xmul
      else if( dphi(i,j) .lt. 0. ) then
          vf0(i,j) = 0.
          trho(i,j) = rhog
          tmu(i,j) = xmug
      end if

  100  continue

      
      do 200 i=1,nx-1
      do 200 j=1,ny-1

      vf_sx(i,j) = 0.5*(vf0(i,j)+vf0(i+1,j))
      vf_sy(i,j) = 0.5*(vf0(i,j)+vf0(i,j+1))

 200  continue

      call bound

      return
      end


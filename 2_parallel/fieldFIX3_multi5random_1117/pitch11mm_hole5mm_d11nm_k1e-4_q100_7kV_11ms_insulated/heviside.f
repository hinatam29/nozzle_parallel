
      subroutine heviside
      include 'subcom.inc'


      do 100 i=1, nx
      do 100 j=1, ny
        alp = 1.e-5
c       alp = 1.5*sdy(i,j)
        if( dphi(i,j) .lt. -alp ) then
           vf0(i,j) = 0.
        else if( abs(dphi(i,j)) .le. alp ) then
           vf0(i,j) = 0.5*(1.+dphi(i,j)/alp+1./pi*sin(pi*dphi(i,j)/alp))
        else 
           vf0(i,j) = 1.
        end if
100   continue

      do 200 i=1, nx-1
      do 200 j=1, ny-1
        vf_sx(i,j) = 0.5*(vf0(i,j)+vf0(i+1,j))
        vf_sy(i,j) = 0.5*(vf0(i,j)+vf0(i,j+1))
200   continue

      return
      end  


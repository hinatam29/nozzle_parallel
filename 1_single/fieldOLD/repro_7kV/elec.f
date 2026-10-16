
      subroutine elec
      include 'subcom.inc'

c     sigl = 0.365
c     sigg = 2.5e-10

      dt2 = 1.e-10

       do 200 i=2, iinj-1
c      do 200 j=2, jnoz-5
       do 200 j=2, jnoz+2

       ijk = jj(i,j)
       alp = 1.e-5
       if( abs(dphi(i,j)) .le. alp .and. vf0(i,j) .ge. 0.5) then
c      if( vf0(i,j) .ge. 0.9 ) then
       gval = sqrt(gradx(i,j)*gradx(i,j)+grady(i,j)*grady(i,j))
       xnr = gradx(i,j)/gval
       xnz = grady(i,j)/gval
       Erl = Er(ijk+1)
       Ezl = Ez(ijk+1)
       efld = Erl*Erl+Ezl*Ezl
       efld = sqrt(efld)


       divqv = (rhoc0(i,j)*ux0(i,j)-rhoc0(i-1,j)*ux0(i-1,j))/sdx(i,j)
     &       + rhoc0(i,j)*ux0(i,j)/x(i,j)
     &       + (rhoc0(i,j+1)*uy0(i,j+1)-rhoc0(i,j)*uy0(i,j))/sdy(i,j)
       
c      rhoc0(i,j) = rhoc0(i,j)-dt*divqv
       rhoc0(i,j) = rhoc0(i,j)-dt2*tsig(i,j)*Ezl
       rhoc0(i,j) = rhoc0(i,j)+dt2*tsig(i,j)*Erl

c      if( abs(rhoc0(i,j)) .ge. 5.e-2 ) then
c         rhoc0(i,j) = 5.e-2*rhoc0(i,j)/abs(rhoc0(i,j))
c      end if

c      write(*,*) rhoc0(i,j)
c      rmax = 4.e-6 

       end if

 200  continue

c     call bound

      return
      end




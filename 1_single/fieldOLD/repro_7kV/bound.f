      subroutine bound
      include 'subcom.inc'

******************* u ********************

****** axitial ******
      do 100 j=1, ny-1
      ux_sx(2,j) = 0.
      ux_sy(1,j) = -ux_sy(2,j)
      ux0(1,j) = -ux0(2,j)

      vf0(1,j) = vf0(2,j)

      vf_sx(2,j) = vf_sx(3,j)
c     vf_sy(2,j) = vf_sy(3,j)
      vf_sy(1,j) = vf_sy(2,j)

      rhoc0_sx(2,j) = rhoc0_sx(3,j)
      rhoc0_sy(1,j) = rhoc0_sy(2,j)
      rhoc0(1,j) = rhoc0(2,j)

 100   continue

****** top ******
      do 110 i=1, iinj
      ux_sx(i,ny-1) = -ux_sx(i,ny-2) 
      ux_sy(i,ny-1) = 0. 
      ux0(i,ny-1) = -ux0(i,ny-2) 

      vf0(i,ny-1) = vf0(i,ny-2)
      vf_sx(i,ny-1) = vf_sx(i,ny-2)
      vf_sy(i,ny-1) = vf_sy(i,ny-2)

      vf0(i,ny-1) = 1.
      vf0(i,ny-2) = 1.

      rhoc0_sx(i,ny-1) = rhoc0_sx(i,ny-2)
      rhoc0_sy(i,ny-1) = rhoc0_sy(i,ny-2)
      rhoc0(i,ny-1) = rhoc0(i,ny-2)

 110   continue

****** top2 ******
c     do 115 i=inoz+1, nx-1
c     ux_sx(i,ny-1) = ux_sx(i,ny-2)
c     ux_sy(i,ny-1) = ux_sy(i,ny-2)
c     ux0(i,ny-1) = ux0(i,ny-2)
c     ux_sx(i,ny-1) = -ux_sx(i,ny-2)
c     ux_sy(i,ny-1) = 0.
c     ux0(i,ny-1) = -ux0(i,ny-2)

c     vf0(i,ny-1) = vf0(i,ny-2)
c     vf_sx(i,ny-1) = vf_sx(i,ny-2)
c     vf_sy(i,ny-1) = vf_sy(i,ny-2)

c115  continue

****** bottom ******
c     do 120 i=1, ihol
      do 120 i=1, nx-1
      ux_sx(i,1) = ux_sx(i,2) 
      ux_sy(i,2) = ux_sy(i,3)
      ux0(i,1) = ux0(i,2)

      vf0(i,2) = vf0(i,3)
      vf0(i,1) = vf0(i,2)
c     vf0(i,0) = vf0(i,3)

      vf_sx(i,2) = vf_sx(i,3)
      vf_sx(i,1) = vf_sx(i,2)
      vf_sy(i,2) = vf_sy(i,3)

      rhoc0_sx(i,1) = rhoc0_sx(i,2)
      rhoc0_sy(i,2) = rhoc0_sy(i,3)
      rhoc0(i,1) = rhoc0(i,2)

 120   continue

****** right side ******
c     do 130 j=jhol+1, jnoz
      do 130 j=1, ny-1
c     ux_sx(nx-1,j) = ux_sx(nx-2,j) 
c     ux_sy(nx-1,j) = ux_sy(nx-2,j)
c     ux0(nx-1,j) = ux0(nx-2,j)

      ux_sx(nx-1,j) = 0.
      ux_sy(nx-1,j) = -ux_sy(nx-2,j)
      ux0(nx-1,j) = -ux0(nx-2,j)

      vf0(nx-1,j) = vf0(nx-2,j)
c     vf0(nx,j) = vf0(nx-3,j)

      vf_sx(nx-1,j) = vf_sx(nx-2,j)
      vf_sy(nx-1,j) = vf_sy(nx-2,j)

      rhoc0_sx(nx-1,j) = rhoc0_sx(nx-2,j)
      rhoc0_sy(nx-1,j) = rhoc0_sy(nx-2,j)
      rhoc0(nx-1,j) = rhoc0(nx-2,j)

 130   continue

****** inside nozzle  ******
      do 140 j=jnoz+1, ny-1
      i = iinj

      ux_sx(i,j) = 0.
      ux_sy(i,j) = -ux_sy(i-1,j)
      ux0(i,j) = -ux0(i-1,j)

c     vf0(i,j) = vf0(i-1,j)
c     vf_sx(i,j) = vf_sx(i-1,j)
c     vf_sy(i,j) = vf_sy(i-1,j)
      vf0(i,j) = 1.0
      vf_sx(i,j) = 1.0
      vf_sy(i,j) = 1.0

      rhoc0(i,j) = rhoc0(i-1,j)
      rhoc0_sx(i,j) = rhoc0_sx(i-1,j)
      rhoc0_sy(i,j) = rhoc0_sy(i-1,j)

 140   continue

****** outside nozzle ******
c     do 150 j=jnoz+1, ny-1
c     i=inoz
c     ux_sx(i,j) = 0.
c     ux_sy(i,j) = -ux_sy(i-1,j)
c     ux0(i,j) = -ux0(i-1,j) 

c     vf0(i,j) = vf0(i-1,j)
c     vf0(nx,j) = vf0(nx-1,j)

c     vf_sx(i,j) = vf_sx(i-1,j)
c     vf_sy(i,j) = vf_sy(i-1,j)

c150   continue

****** bottom nozzle ******
c     do 160 i=iinj+1, inoz
      do 160 i=iinj+1, nx-1
      j = jnoz

      ux_sx(i,j) = -ux_sx(i,j-1)
      ux_sy(i,j) = 0.
      ux0(i,j) = -ux0(i,j-1)

      vf0(i,j) = vf0(i,j-1)
      vf_sx(i,j) = vf_sx(i,j-1)
      vf_sy(i,j) = vf_sy(i,j-1)

      rhoc0(i,j) = rhoc0(i,j-1)
      rhoc0_sx(i,j) = rhoc0_sx(i,j-1)
      rhoc0_sy(i,j) = rhoc0_sy(i,j-1) 

 160   continue

****** top nozzle ******
c     do 170 i=iinj+1, inoz
      do 170 i=iinj+1, nx-1
      ux_sx(i,ny-1) = -ux_sx(i,ny-2)
      ux_sy(i,ny-1) = 0.
      ux0(i,ny-1) = -ux0(i,ny-2)

      vf0(i,ny-1) = vf0(i,ny-2)
      vf_sx(i,ny-1) = vf_sx(i,ny-2)
      vf_sy(i,ny-1) = vf_sy(i,ny-2)

      rhoc0(i,ny-1) = rhoc0(i,ny-2)
      rhoc0_sx(i,ny-1) = rhoc0_sx(i,ny-2)
      rhoc0_sy(i,ny-1) = rhoc0_sy(i,ny-2) 

 170   continue

c***** inside plate ******
c     do 175 j=1, jhol
c     i=ihol
c     ux_sx(i,j) = 0.
c     ux_sy(i,j) = -ux_sy(i-1,j)
c     ux0(i,j) = -ux0(i-1,j)

c     vf0(i,j) = vf0(i-1,j)
c     vf_sx(i,j) = vf_sx(i-1,j)
c     vf_sy(i,j) = vf_sy(i-1,j)

c175   continue

****** outside plate ******
c     do 180 j=1, jhol
c     ux_sx(nx-1,j) = 0.
c     ux_sy(nx-1,j) = -ux_sy(nx-2,j)
c     ux0(nx-1,j) = -ux0(nx-2,j)

c     vf0(nx-1,j) = vf0(nx-2,j)
c     vf_sx(nx-1,j) = vf_sx(nx-2,j)
c     vf_sy(nx-1,j) = vf_sy(nx-2,j)

c180   continue

****** bottom plate ******
c     do 185 i=ihol+1, nx-1
c     ux_sx(i,1) = -ux_sx(i,2)
c     ux_sy(i,2) = 0.
c     ux0(i,1) = -ux0(i,2)

c     vf0(i,1) = vf0(i,2)
c     vf_sx(i,1) = vf_sx(i,2)
c     vf_sy(i,1) = vf_sy(i,2)

c185   continue

****** top plate ******
c     do 190 i=ihol+1, nx-1
c     j=jhol
c     ux_sx(i,j) = -ux_sx(i,j-1)
c     ux_sy(i,j) = 0.
c     ux0(i,j) = -ux0(i,j-1)

c     vf0(i,j) = vf0(i,j-1)
c     vf_sx(i,j) = vf_sx(i,j-1)
c     vf_sy(i,j) =vf_sy(i,j-1)

c190   continue

****** nozzle ******
      do 195 j=jnoz+1, ny-1
c     do 195 i=iinj+1, inoz
      do 195 i=iinj+1, nx-1
      ux_sx(i,j) = 0.
      ux_sy(i,j) = 0.
      ux0(i,j) = 0.

      vf0(i,j) = 0.
      vf_sx(i,j) = 0.
      vf_sy(i,j) = 0.

      rhoc0(i,j) = 0.
      rhoc0_sx(i,j) = 0.
      rhoc0_sy(i,j) = 0.

      p0(i,j) = 1.013e5 
 195   continue

****** plate ******
c      do 197 j=1, jhol
c      do 197 i=ihol+1, nx-1
c      ux_sx(i,j) = 0.
c      ux_sy(i,j) = 0.
c      ux0(i,j) = 0.

c      vf0(i,j) = 0.
c      vf_sx(i,j) = 0.
c      vf_sy(i,j) = 0.

c197    continue



******************* v *******************

c     uin = -10.e-2
c     uin = -15.e-2
      uin = -20.e-2
      umax = 2.*uin

****** axitial ******
      do 200 j=1, ny-1
      uy_sx(2,j) = uy_sx(3,j) 
      uy_sy(1,j) = uy_sy(2,j)
      uy0(1,j) = uy0(2,j)
 200   continue

****** top ******
      do 210 i=1,iinj
c     uin = umax*(1.-x(i,ny)*x(i,ny)/x(iinj,ny)/x(iinj,ny))
      uy_sx(i,ny-1) = 2.*uin-uy_sx(i,ny-2)
      uy_sy(i,ny-1) = uin
      uy0(i,ny-1) = 2.*uin-uy0(i,ny-2)
c     uy_sx(i,ny-1) = uy_sx(i,ny-2)
c     uy_sy(i,ny-1) = uy_sy(i,ny-2)
c     uy0(i,ny-1) = uy0(i,ny-2)

 210   continue

****** top2 ******
c     do 215 i=inoz+1, nx-1
c     uy_sx(i,ny-1) = uy_sx(i,ny-2)
c     uy_sy(i,ny-1) = uy_sy(i,ny-2)
c     uy0(i,ny-1) = uy0(i,ny-2)
c     uy_sx(i,ny-1) = -uy_sx(i,ny-2)
c     uy_sy(i,ny-1) = 0.
c     uy0(i,ny-1) = -uy0(i,ny-2)

c215  continue

****** bottom ******
c     do 220 i=1, ihol
      do 220 i=1, ny-1
      uy_sx(i,1) = uy_sx(i,2)
      uy_sy(i,2) = uy_sy(i,3) 
      uy0(i,1) = uy0(i,2)

 220   continue 

****** right side ******
c     do 230 j=jhol+1, jnoz
      do 230 j=1, ny-1
c     uy_sx(nx-1,j) = uy_sx(nx-2,j)
c     uy_sy(nx-1,j) = uy_sy(nx-2,j)
c     uy0(nx-1,j) = uy0(nx-2,j)
      uy_sx(nx-1,j) = 0.
      uy_sy(nx-1,j) = -uy_sy(nx-2,j)
      uy0(nx-1,j) = -uy0(nx-2,j)

 230   continue

****** inside nozzle ******
      do 240 j=jnoz+1, ny-1
      i = iinj

      uy_sx(i,j) = 0.
      uy_sy(i,j) = -uy_sy(i-1,j)
      uy0(i,j) = -uy0(i-1,j)

 240   continue

****** outside nozzle ******
c     do 250 j=jnoz+1, ny-1
c     i=inoz
c     uy_sx(i,j) = 0.
c     uy_sy(i,j) = -uy_sy(i-1,j)
c     uy0(i,j) = -uy0(i-1,j)

c250   continue

****** bottom nozzle ******
c     do 260 i=iinj+1, inoz
      do 260 i=iinj+1, nx-1
      j = jnoz
   
      uy_sx(i,j) = -uy_sx(i,j-1)
      uy_sy(i,j) = 0.
      uy0(i,j) = -uy0(i,j-1)

 260   continue

****** top nozzle ******
c     do 270 i=iinj+1, inoz
      do 270 i=iinj+1, nx-1
      uy_sx(i,ny-1) = -uy_sx(i,ny-2)
      uy_sy(i,ny-1) = 0.
      uy0(i,ny-1) = -uy0(i,ny-2)

 270   continue

****** inside plate ******
c     do 280 j=1, jhol
c     i=ihol
c     uy_sx(i,j) = 0.
c     uy_sy(i,j) = -uy_sy(i-1,j)
c     uy0(i,j) = -uy0(i-1,j)

c280   continue

****** outiside plate ******
c     do 290 j=1, jhol
c     uy_sx(nx-1,j) = 0.
c     uy_sy(nx-1,j) = -uy_sy(nx-2,j)
c     uy0(nx-1,j) = -uy0(nx-2,j)

c290   continue

****** bottom plate ******
c     do 300 i=ihol+1, nx-1
c     uy_sx(1,j) = -uy_sx(2,j)
c     uy_sy(2,j) = 0.
c     uy0(1,j) = -uy0(2,j)

c300   continue

****** top plate ******
c     do 310 i=ihol+1, nx-1
c     j=jhol
c     uy_sx(i,j) = -uy_sx(i,j-1)
c     uy_sy(i,j) = 0.
c     uy0(i,j) = -uy0(i,j-1)

c310   continue

c***** nozzle ******
      do 320 j=jnoz+1, ny-1
c     do 320 i=iinj+1, inoz
      do 320 i=iinj+1, nx-1
 
      uy_sx(i,j) = 0.
      uy_sy(i,j) = 0.
      uy0(i,j) = 0.

 320   continue

c     do 330 i=1, nx-1
c     do 330 j=1, ny-1
      do 330 i=2, iinj-1
      do 330 j=1, jnoz+1


      alp = sqrt(sdx1(i,j)*sdx1(i,j)+sdy1(i,j)*sdy1(i,j))

      alp = 5.e-6
      if( abs(dphi(i,j)) .gt. alp ) then
         rhoc0(i,j) = 0.
         rhoc0_sx(i,j) = 0.
         rhoc0_sy(i,j) = 0.
      else
         gval = sqrt(gradx(i,j)*gradx(i,j)+grady(i,j)*grady(i,j))
         xnr = gradx(i,j)/gval
         xnz = grady(i,j)/gval
         if( gval .eq. 0. .or. xnz .eq. 0. ) go to 330

c        fct1 = 1./sdx(i,j)+1./x(i,j)+1./sdy(i,j)*xnr/xnz
c        RHS = 1./sdx(i,j)*ux0(i-1,j)-1./sdy(i,j)*uy0(i,j+1)
c        ux0(i,j) = RHS/fct1
c        uy0(i,j) = -xnr/xnz*ux0(i,j)

         fct1 = 1. + (xnr/xnz)**2.
         RHS = ux0(i,j+1)*ux0(i,j+1)+uy0(i,j+1)*uy0(i,j+1)
c        ux0(i,j) = -sqrt(RHS/fct1)
c        uy0(i,j) = -xnr/xnz*ux0(i,j)

c        write(*,*) istp, xnr, xnz
c        ux_sx(i,j) = 0.5*(ux0(i,j)+ux0(i-1,j))
c        ux_sy(i,j) = 0.5*(ux0(i,j)+ux0(i,j+1))
c        uy_sx(i,j) = 0.5*(uy0(i,j)+uy0(i-1,j))
c        uy_sy(i,j) = 0.5*(uy0(i,j)+uy0(i,j+1))
      end if
c     if( i .gt.iinj-2 ) then
c        rhoc0(i,j) = 0.
c        rhoc0_sx(i,j) = 0.
c        rhoc0_sy(i,j) = 0.
c     end if
330   continue

c     do 330 j=jnoz-1, ny-1
c     do 320 i=iinj+1, inoz
c     do 330 i=1, nx-1
 
c     rhoc0_sx(i,j) = 0.
c     rhoc0_sy(i,j) = 0.
c     rhoc0(i,j) = 0.
c330   continue


****** plate ******
c     do 330 j=1, jhol
c     do 330 i=ihol+1,nx-1
c     uy_sx(i,j) = 0.
c     uy_sy(i,j) = 0.
c     uy0(i,j) = 0.

c330   continue


      end

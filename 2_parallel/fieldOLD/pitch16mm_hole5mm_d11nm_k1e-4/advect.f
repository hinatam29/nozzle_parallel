      subroutine advect
      include 'subcom.inc'

c:    x-direction
c:    ------------------------------------
c     call cipcsl3x(trho_sx,trho0,trhpx,trh0p)
      call cipcsl3x(ux_sx,ux0,uxpx,ux0p)
      call cipcsl3x(uy_sx,uy0,uypx,uy0p)
c     call cipcsl3x(pp_sx,pp0,pppx,pp0p)
      call cipcsl3x(vf_sx,vf0,vfpx,vf0p)
c     call cipcsl3x(dp_sx,dphi,dppx,dp0p)
c     call cipcsl3x(rhoc0_sx,rhoc0,rhocpx,rhoc0p)

      call crrctx(uxpx)
      call crrctx(uypx)
      call crrctx(vfpx)
c     call crrctx(rhocpx)

c     call TECSIAy(trho_sy,trho0,trhpy,trh0p)
      call TECSIAy(ux_sy,ux0,uxpy,ux0p)
      call TECSIAy(uy_sy,uy0,uypy,uy0p)
c     call TECSIAy(pp_sy,pp0,pppy,pp0p)
      call TECSIAy(vf_sy,vf0,vfpy,vf0p)
c     call TECSIAy(dp_sy,dphi,dppy,dp0p)
c     call TECSIAy(rhoc0_sy,rhoc0,rhocpy,rhoc0p)

c     call result(trho_sx,trho_sy,trho0,trhpx,trhpy,trh0p)
      call result(ux_sx,ux_sy,ux0,uxpx,uxpy,ux0p)
      call result(uy_sx,uy_sy,uy0,uypx,uypy,uy0p)
c     call result(pp_sx,pp_sy,pp0,pppx,pppy,pp0p)
      call result(vf_sx,vf_sy,vf0,vfpx,vfpy,vf0p)
c     call result(dp_sx,dp_sy,dphi,dppx,dppy,dp0p)
c     call result(rhoc0_sx,rhoc0_sy,rhoc0,rhocpx,rhocpy,rhoc0p)

      call bound

c:    y-direction
c:    ------------------------------------
c     call cipcsl3y(trho_sy,trho0,trhpy,trh0p)
      call cipcsl3y(ux_sy,ux0,uxpy,ux0p)
      call cipcsl3y(uy_sy,uy0,uypy,uy0p)
c     call cipcsl3y(pp_sy,pp0,pppy,pp0p)
      call cipcsl3y(vf_sy,vf0,vfpy,vf0p)
c     call cipcsl3y(dp_sy,dphi,dppy,dp0p)
c     call cipcsl3y(rhoc0_sy,rhoc0,rhocpy,rhoc0p)

      call crrcty(uxpy)
      call crrcty(uypy)
      call crrcty(vfpy)
c     call crrcty(rhocpy)

c     call TECSIAx(trho_sx,trho0,trhpx,trh0p)
      call TECSIAx(ux_sx,ux0,uxpx,ux0p)
      call TECSIAx(uy_sx,uy0,uypx,uy0p)
c     call TECSIAx(pp_sx,pp0,pppx,pp0p)
      call TECSIAx(vf_sx,vf0,vfpx,vf0p)
c     call TECSIAx(dp_sx,dphi,dppx,dp0p)
c     call TECSIAx(rhoc0_sx,rhoc0,rhocpx,rhoc0p)

c     call result(trho_sx,trho_sy,trho0,trhpx,trhpy,trh0p)
      call result(ux_sx,ux_sy,ux0,uxpx,uxpy,ux0p)
      call result(uy_sx,uy_sy,uy0,uypx,uypy,uy0p)
c     call result(pp_sx,pp_sy,pp0,pppx,pppy,pp0p)
      call result(vf_sx,vf_sy,vf0,vfpx,vfpy,vf0p)
c     call result(dp_sx,dp_sy,dphi,dppx,dppy,dp0p)
c     call result(rhoc0_sx,rhoc0_sy,rhoc0,rhocpx,rhocpy,rhoc0p)

      call bound


      return
      end
     

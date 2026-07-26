
      subroutine proj
      include 'subcom.inc'

      call proj_siax(ux_sx,uxpx,p0)
      call TECVIAx(ux_sx,ux0,uxpx,ux0p)
      call TECSIAy(ux_sy,ux0,uxpy,ux0p)

      call result(ux_sx,ux_sy,ux0,uxpx,uxpy,ux0p)
      call bound


      call proj_siay(uy_sy,uypy,p0)
      call TECVIAy(uy_sy,uy0,uypy,uy0p)
      call TECSIAx(uy_sx,uy0,uypx,uy0p)

      call result(uy_sx,uy_sy,uy0,uypx,uypy,uy0p)
      call bound


      return
      end
     

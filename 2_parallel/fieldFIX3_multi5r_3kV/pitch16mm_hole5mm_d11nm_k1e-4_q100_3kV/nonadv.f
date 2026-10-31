
      subroutine nonadv
      include 'subcom.inc'

      call crrctx(trhpx)
      call TECSIAy(trho_sy,trho0,trhpy,trh0p)
      call result(trho_sx,trho_sy,trho0,trhpx,trhpy,trh0p)
      call bound

      call crrcty(rhopy)
      call TECSIAx(trho_sy,trho0,trhpy,trh0p)
      call result(trho_sx,trho_sy,trho0,trhpx,trhpy,trh0p)
      call bound

      return
      end


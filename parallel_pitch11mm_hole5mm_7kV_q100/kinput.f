
      subroutine kinput
      include 'subcom.inc'

      do 100 i=1, nx-1
      do 100 j=1, ny-1
        ijk = jj(i,j)

        read(300,*) istp, time
        read(300,*) p0(i,j),vf0(i,j),ux0(i,j),uy0(i,j)
        read(300,*) ux_sx(i,j),uy_sx(i,j),ux_sy(i,j),uy_sy(i,j)
        read(300,*) vf_sx(i,j),vf_sy(i,j)
c       read(300,*) phi(ijk), Er(ik), Ez(ijk), rhoc(ijk)
        read(300,*) phi(ijk), Er(ijk), Ez(ijk), rhoc0(i,k)
c       read(300,*) fer(ijk), fez(ijk)
        read(300,*) dphi(i,j), fsvx(i,j), fsvy(i,j)
        read(300,*) trho(i,j),tmu(i,j)
        read(300,*) epsi(i,j), sig(i,j), rhoc0(i,j)
100   continue

      return
      end

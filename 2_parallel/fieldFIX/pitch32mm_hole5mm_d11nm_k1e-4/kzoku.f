
      subroutine kzoku
      include 'subcom.inc'

      do 100 i=1, nx-1
      do 100 j=1, ny-1
        ijk = jj(i,j)

        write(200,*) istp, time
        write(200,*) p0(i,j),vf0(i,j),ux0(i,j),uy0(i,j)
        write(200,*) ux_sx(i,j),uy_sx(i,j),ux_sy(i,j),uy_sy(i,j)
        write(200,*) vf_sx(i,j),vf_sy(i,j)
c       write(200,*) phi(ijk), Er(ijk), Ez(ijk), rhoc(ijk)
        write(200,*) phi(ijk), Er(ijk), Ez(ijk), rhoc0(i,k)
c       write(200,*) fer(ijk), fez(ijk)
        write(200,*) dphi(i,j), fsvx(i,j),fsvy(i,j)
        write(200,*) trho(i,j), tmu(i,j)
        write(200,*) epsi(i,j), sig(i,j), rhoc0(i,j)
100   continue

      rewind 200

      return
      end

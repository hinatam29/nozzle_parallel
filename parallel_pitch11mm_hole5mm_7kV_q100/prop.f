
      subroutine prop
      include 'subcom.inc'

      do 5 i=1, nx
      do 5 j=1, ny
        if( vf0(i,j) .gt. 1.0 ) then
          vf0(i,j) = 1.0
        else
          if( vf0(i,j) .lt. 0. ) then
            vf0(i,j) = 0.
          end if
        end if

c       if( vf_sx(i,j) .gt. 1.0 ) then
c         vf_sx(i,j) = 1.0
c       else
c         if( vf_sx(i,j) .lt. 0. ) then
c           vf_sx(i,j) = 0.
c         end if
c       end if

c       if( vf_sy(i,j) .gt. 1.0 ) then
c         vf_sy(i,j) = 1.0
c       else
c         if( vf_sy(i,j) .lt. 0. ) then
c           vf_sy(i,j) = 0.
c         end if
c       end if
5     continue


      do 140 i=1, nx-1
      trho(i,1)  = (1-vf0(i,1))*rhog  + vf0(i,1)*rhol
      trho(i,ny-1) = (1-vf0(i,ny-1))*rhog + vf0(i,ny-1)*rhol
      tmu(i,1)  = (1-vf0(i,1))*xmug  + vf0(i,1)*xmul
      tmu(i,ny-1) = (1-vf0(i,ny-1))*xmug + vf0(i,ny-1)*xmul

      svg = sqrt(1.402*p0(i,1)/rhog)
      tsv(i,1)  = (1-vf0(i,1))*svg  + vf0(i,1)*svl
      tsv(i,ny-1) = (1-vf0(i,ny-1))*svg + vf0(i,ny-1)*svl

      tsig (i,1) = (1-vf0(i,1))*sigg + vf0(i,1)*sigl
      tsig (i,ny-1) = (1-vf0(i,ny-1))*sigg + vf0(i,ny-1)*sigl
140   continue

      do 150 j=1, ny-1
      trho(1,j)  = (1.-vf0(1,j))*rhog  + vf0(1,j)*rhol
      trho(nx-1,j) = (1.-vf0(nx-1,j))*rhog + vf0(nx-1,j)*rhol
      tmu(1,j)  = (1.-vf0(1,j))*xmug  + vf0(1,j)*xmul
      tmu(nx-1,j) = (1.-vf0(nx-1,j))*xmug + vf0(nx-1,j)*xmul

c     svg = 1.402*p0(1,j)/rhog
c     write(*,*) j, p0(1,j), uy0(1,j)
c     write(*,*) j, p0(1,j), svg 
c     write(*,*) x(i,j), y(i,j), p0(i,j)

      svg = sqrt(1.402*p0(1,j)/rhog)
      tsv(1,j)  = (1-vf0(1,j))*svg  + vf0(1,j)*svl
      tsv(nx-1,j) = (1-vf0(nx-1,j))*svg + vf0(nx-1,j)*svl

      tsig(1,j) = (1-vf0(1,j))*sigg + vf0(1,j)*sigl
      tsig(nx-1,j) = (1-vf0(nx-1,j))*sigg + vf0(nx-1,j)*sigl
150   continue

      do 200 i=2, nx-2
      do 200 j=2, ny-2

      eps0 = 2.0*sqrt(sdx(i,j)**2. + sdy(i,j)**2.)
c     eps0 = 0.5*sdx(i,j)
      xhev = 0.5*max(-1.,
     &           min(1.,dphi(i,j)/eps0+1./pi*sin(pi*dphi(i,j)/eps0)))

      svg = sqrt(1.402*p0(i,j)/rhog)
c     trho(i,j) = 0.5*(rhol+rhog)+(rhol-rhog)*xhev
c     tmu(i,j)  = 0.5*(xmul+xmug)+(xmul-xmug)*xhev
c     tsv(i,j)  = 0.5*(svl+svg)+(svl-svg)*xhev
c     write(*,*) xhev, trho(i,j)

c     write(*,*) i, j, p0(i,j)
      trho(i,j)  = (1.-vf0(i,j))*rhog  + vf0(i,j)*rhol
      tmu(i,j)  = (1.-vf0(i,j))*xmug  + vf0(i,j)*xmul
      tsv(i,j)  = (1.-vf0(i,j))*svg  + vf0(i,j)*svl
      tsig(i,j) = (1-vf0(i,j))*sigg + vf0(i,j)*sigl
200   continue
 
5000  return
      end

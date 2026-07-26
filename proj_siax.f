
      subroutine proj_siax(q_s,qp,q0)
      include 'subcom.inc'

      dimension q_s(0:nn,0:nn), qp(0:nn,0:nn), q0(0:nn,0:nn)

      do 100 i=2,nx-1
      do 100 j=2,ny-2
c       rx0 = vf0(i-1,j)*rhol+(1.-vf0(i-1,j))*rhog
c       rx1 = vf0(i,j)*rhol+(1.-vf0(i,j))*rhog
c       trho = 0.5*(rx0+rx1)

        trhom = 0.5*(trho(i-1,j)+trho(i,j))

c       qp(i,j) = q_s(i,j)-dt*(-q0(i-1,j)+q0(i,j))/dx(i,j)
        qp(i,j) = q_s(i,j)-dt*(-q0(i-1,j)+q0(i,j))/trhom/sdx(i,j)
100   continue

      return
      end
 


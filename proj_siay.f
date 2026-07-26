
      subroutine proj_siay(q_s,qp,q0)
      include 'subcom.inc'

      dimension q_s(0:nn,0:nn), qp(0:nn,0:nn), q0(0:nn,0:nn)

      do 100 i=2,nx-2
      do 100 j=2,ny-1
        ry0 = vf0(i,j-1)*rhol+(1.-vf0(i,j-1))*rhog
        ry1 = vf0(i,j)*rhol+(1.-vf0(i,j))*rhog
        trho = 0.5*(ry0+ry1)

        trhom = 0.5*(trho(i,j-1)+trho(i,j))

c       qp(i,j) = q_s(i,j)-dt*(-q0(i,j-1)+q0(i,j))/dy(i,j)
        qp(i,j) = q_s(i,j)-dt*(-q0(i,j-1)+q0(i,j))/trhom/sdy(i,j)
100   continue

      return
      end
 


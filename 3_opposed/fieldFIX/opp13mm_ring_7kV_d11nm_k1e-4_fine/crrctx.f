
      subroutine crrctx(qp)
      include 'subcom.inc'

      dimension qp(0:nn,0:nn)

      do 100 j=2, ny-2
      do 110 i=2, nx-1

        xx = 0.5*(x(i+1,j)+x(i,j))

        ux=(ux0(i,j)-ux0(i-1,j))/sdx1(i,j)
     &     +ux_sx(i,j)/x(i,j)

        vy=0.5*((uy0(i-1,j+1)+uy0(i,j+1))
     &         -(uy0(i-1,j-1)+uy0(i,j-1)))
     &        / (2.*sdy1(i,j))
 
        qp(i,j) = qp(i,j)-dt*qp(i,j)*(ux+vy)
c       qp(i,j) = qp(i,j)-dt*qp(i,j)*(ux)

110   continue
100   continue

      return
      end



      subroutine crrcty(qp)
      include 'subcom.inc'

      dimension qp(0:nn,0:nn)


      do 100 i=2, nx-2
      do 110 j=2, ny-1

        xx = 0.5*(x(i+1,j)+x(i,j))

        ux=0.5*((ux0(i+1,j-1)+ux0(i+1,j))
     &         -(ux0(i-1,j-1)+ux0(i-1,j)))
     &        / (2.*sdx1(i,j))
     &     +ux_sx(i,j)/x(i,j)
 
        vy=(uy0(i,j)-uy0(i,j-1))/sdy1(i,j)

        qp(i,j) = qp(i,j)-dt*qp(i,j)*(ux+vy)
c       qp(i,j) = qp(i,j)-dt*qp(i,j)*(vy)
110   continue
100   continue

      return
      end





      subroutine TECVIAx(q_s,q0,qp,q0p)
      include 'subcom.inc'

      dimension q_s(0:nn,0:nn), q0(0:nn,0:nn)
      dimension qp(0:nn,0:nn), q0p(0:nn,0:nn) 

      do 100 i=2, nx-2
      do 100 j=2, ny-2
        q0p(i,j) = q0(i,j)
     &           +0.5*(qp(i,j)-q_s(i,j)+qp(i+1,j)-q_s(i+1,j))
100   continue

      return
      end 
    


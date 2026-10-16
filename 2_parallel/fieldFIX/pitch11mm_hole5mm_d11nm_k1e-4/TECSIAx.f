
      subroutine TECSIAx(q_s,q0,qp,q0p)
      include 'subcom.inc'

      dimension q_s(0:nn,0:nn), q0(0:nn,0:nn) 
      dimension qp(0:nn,0:nn), q0p(0:nn,0:nn) 

      do 100 i=2, nx-1
      do 100 j=2, ny-2
        qp(i,j) = TEC*(q_s(i,j)
     &           +0.5*(q0p(i-1,j)-q0(i-1,j)+q0p(i,j)-q0(i,j)))
     &           +(1.-TEC)*0.5*(q0p(i-1,j)+q0p(i,j))
100   continue

      return
      end 
    


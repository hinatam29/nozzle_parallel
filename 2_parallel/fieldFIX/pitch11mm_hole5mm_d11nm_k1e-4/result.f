
      subroutine result(q_sx,q_sy,q0,qpx,qpy,q0p)
      include 'subcom.inc'

      dimension q_sx(0:nn,0:nn), q_sy(0:nn,0:nn), q0(0:nn,0:nn) 
      dimension qpx(0:nn,0:nn), qpy(0:nn,0:nn), q0p(0:nn,0:nn)

  
      do 100 i=2, nx-1
      do 100 j=2, ny-2
        q_sx(i,j) = qpx(i,j)
100   continue

      do 200 i=2, nx-2
      do 200 j=2, ny-1
        q_sy(i,j) = qpy(i,j)
200   continue 
      
      do 300 i=2, nx-2
      do 300 j=2, ny-2
        q0(i,j) = q0p(i,j)
300   continue 

      return
      end

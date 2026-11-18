
       subroutine pserch
       include 'subcom.inc'


       real l1,l2,l3,l4

c      alp = 0. 
c      alp = 0.15e-3 
c      bet = 0.1e-3

c      --- initialize ---

       do 5 i=1, nx-1
       do 5 j=1, ny-1
       ijk = jj(i,j)

c      Ernex(ijk) = Er(ijk)
c      Eznex(ijk) = Ez(ijk)*0.1

 5     continue

c      --- grid shift ---

       do 10 i=1, nx-1
       do 10 j=1, ny-1

c      xnex(i,j) = x(i,j) + alp
c      ynex(i,j) = y(i,j) + bet

c10    continue

c      --- particle position search---

c      do 100 k=1, nx-1
c      do 100 l=1, ny-1

c      if( x(k,1) .gt. xnex(i,j) )then
c        inex = k-1
c        go to 110
c      end if 
c100   continue

c110   continue

c      do 120 l=1, ny-1

c      if( y(1,l) .gt. ynex(i,j) )then
c        jnex = l-1
c        go to 130
c      end if
c120   continue

c130   continue

c      write(*,*) inex, jnex

c      --- electrostatic force interpolation ---

c      write(*,*) inex, jnex, xnex-x(inex,jnex), x(2,1)-x(1,1)
c      write(*,*) sqrt((xnex-x(inex,jnex))**2+(ynex-y(inex,jnex))**2)

c      l1 = sqrt((xnex(i,j)-x(inex,jnex))*(xnex(i,j)-x(inex,jnex)) 
c    &          +(ynex(i,j)-y(inex,jnex))*(ynex(i,j)-y(inex,jnex)))

c      l2 = sqrt((x(inex+1,jnex)-xnex(i,j))*(x(inex+1,jnex)-xnex(i,j))
c    &          +(ynex(i,j)-y(inex+1,jnex))*(ynex(i,j)-y(inex+1,jnex)))

c      l3 = sqrt( (x(inex+1,jnex+1)-xnex(i,j))
c    &           *(x(inex+1,jnex+1)-xnex(i,j))
c    &           +(y(inex+1,jnex+1)-ynex(i,j))
c    &           *(y(inex+1,jnex+1)-ynex(i,j)) )

c      l4 = sqrt((xnex(i,j)-x(inex,jnex+1))*(xnex(i,j)-x(inex,jnex+1))
c    &          +(y(inex,jnex+1)-ynex(i,j))*(y(inex,jnex+1)-ynex(i,j)))

c      write(*,*) l1,l2,l3,l4

c      l1 = 1./l1
c      l2 = 1./l2
c      l3 = 1./l3
c      l4 = 1./l4

c      write(*,*) l1,l2,l3,l4

c      ijk = jj(i,j)
c      ijkd = jj(inex,jnex)

c      write(*,*) ijk, ijkd

c      phinex(ijk) = ( phi(ijkd)*l1+phi(ijkd+(ny-1))*l2
c    &            + phi(ijkd+ny)*l3+phi(ijkd+1)*l4 )
c    &            / ( l1+l2+l3+l4 )       


c      Ernex(ijk) = ( Er(ijkd)*l1+Er(ijkd+(ny-1))*l2
c    &              + Er(ijkd+ny)*l3+Er(ijkd+1)*l4 )
c    &              / ( l1+l2+l3+l4 )

c      write(*,*) Er(ijkd), Ernex(ijk)

c      write(*,*) Eznex(ijk)

c      Eznex(ijk) = ( Ez(ijkd)*l1+Ez(ijkd+(ny-1))*l2
c    &              + Ez(ijkd+ny)*l3+Ez(ijkd+1)*l4 )
c    &              / ( l1+l2+l3+l4 )

c100   continue
       
c      ijk = jj(2,j)
c      write (*,*) y(2,j), ynex(2,j), Ez(ijk), Eznex(ijk)

 10    continue

c      do 300 j=1, ny-1
c      ijk = jj(4,j)

c      write(*,*) y(4,j), ynex(4,j), Ez(ijk), Eznex(ijk)
c    &            0.5*(Ez(ijk+1)-Ez(ijkd))

c300   continue

       return
       end

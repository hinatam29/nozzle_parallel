
       subroutine psearch
       include 'subcom.inc'



c      --- particle position search ---

c      do 10 k=1, ip


c      do 100 i=1, nx-1

c      if( x(i,1) .gt. xp(k) )then
c        inex(k) = i-1
c        go to 110
c      end if

c100   continue
c110   continue

c      do 120 j=1, ny-1

c      if( y(1,j) .gt. yp(k) )then
c        jnex(k) = j-1
c        go to 130
c      end if

c120   continue
c130   continue


c10    continue

c      do 140 i=1, nx-1

c      if( x(i,1) .gt. xp2 )then
c        inex2 = i-1
c        go to 150
c      end if

c140   continue

c150   continue


c      do 160 j=1, ny-1

c      if( y(1,j) .gt. yp2 )then
c        jnex2 = j-1
c        go to 170
c      end if

c160   continue

c170   continue


c      write(*,*) inex2, jnex2

       return 

       end

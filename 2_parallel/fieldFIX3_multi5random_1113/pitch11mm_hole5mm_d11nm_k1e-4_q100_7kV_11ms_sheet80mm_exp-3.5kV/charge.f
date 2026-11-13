
      subroutine charge
      include 'subcom.inc'

      eps0 = 8.8542e-12

c     do 100 i=1, nx-1
c     do 100 k=1, ny-1
      do 100 i=2, nx-2
      do 100 k=2, ny-2
        ijk = jj(i,k) 

c       sdx = x(i+1,k) - x(i,k)
c       sdy = y(i,k+1) - y(i,k)

c       xx2 = 0.5*(x(i+1,k) + x(i,k))
        xx2 = 0.5*(x(i,k) + x(i-1,k))
        xx3 = 0.5*(x(i+2,k) + x(i+1,k))
c       yy2 = 0.5*(y(i,k+1) + y(i,k))
        yy2 = 0.5*(y(i,k) + y(i,k-1))
        yy3 = 0.5*(y(i,k+2) + y(i,k+1))


       alp = 1.e-6
c      if( i .ge. iinj .and. k .ge. jnoz ) then
c         rhoc(ijk) = 0.
c      else

c      if( vf0(i,k) .lt. 0.5 ) then
c      if( vf0(i,k) .eq. 0. ) then
c         rhoc(ijk) = 0.
c      else
       if( abs(dphi(i,k)) .gt. alp ) then
             rhoc0(i,k) = 0.
          else

        rhoc0(i,k) = ( epsi(i+1,k)*eps0*Er(ijk+(ny-1))
     &              - epsi(i-1,k)*eps0*Er(ijk-(ny-1)) )
     &              /( xx3 - xx2 )
     &            + epsi(i,k)*eps0*Er(ijk)/xx2
     &            + ( epsi(i,k+1)*eps0*Ez(ijk+1)
     &              - epsi(i,k-1)*eps0*Ez(ijk-1) )
     &              /( yy3 - yy2 )
c       rhoc0(i,k) = ( epsi(i+1,k)*eps0*Er(ijk+(ny-1))
c    &              - epsi(i-1,k)*eps0*Er(ijk-(ny-1)) )
c    &            + epsi(i,k)*eps0*Er(ijk)/xx2
c       rhoc0(i,k) = ( epsi(i,k+1)*eps0*Ez(ijk+1)
c    &              - epsi(i,k-1)*eps0*Ez(ijk-1) )
       end if

       rhoc0(i,k) = -rhoc0(i,k)

       if( rhoc0(i,k) .le. 0. ) rhoc0(i,k) = 0.        
c      if( vf0(i,k) .le. 0.5)then
c         rhoc(ijk) = 0.
c      end if  
c      write(*,*) rhoc0(i,k)
 
c      if( i .ge. iinj .and. k .ge. jnoz ) then
c         rhoc(ijk) = 0.
c      else

c       grdpx1 = x(i+1,k)*0.5*(epsi(i+1,k)+epsi(i,k))
c    &                  *(phi(ijk+ny)-phi(ijk))/(x(i+1,k)-x(i,k))
c       grdpx2 = x(i-1,k)*0.5*(epsi(i-1,k)+epsi(i,k))
c    &                  *(phi(ijk)-phi(ijk-ny))/(x(i,k)-x(i-1,k))
c       grdpx1 = 0.5*(epsi(i+1,k)+epsi(i,k))
c    &                  *(phi(ijk+ny)-phi(ijk))/(x(i+1,k)-x(i,k))
c       grdpx2 = 0.5*(epsi(i-1,k)+epsi(i,k))
c    &                  *(phi(ijk)-phi(ijk-ny))/(x(i,k)-x(i-1,k))
c       xx = 0.5*(x(i+1,k)-x(i-1,k))

c       grdpy1 = 0.5*(epsi(i,k+1)+epsi(i,k))
c    &                  *(phi(ijk+1)-phi(ijk))/(y(i,k+1)-y(i,k))
c       grdpy2 = 0.5*(epsi(i,k-1)+epsi(i,k))
c    &                  *(phi(ijk)-phi(ijk-1))/(y(i,k)-y(i,k-1))
c       yy = 0.5*(y(i,k+1)-y(i,k-1))

c       rhoc(ijk) = 1./x(i,k)*(grdpx1-grdpx2)/xx 
c       rhoc(ijk) = (grdpx1-grdpx2)/xx 
c    &            + (grdpy1-grdpy2)/yy 
c    &            + epsi(i,k)/x(i,k)
c    &             *(phi(ijk+ny)-phi(ijk-ny))
c    &             *(phi(ijk+ny)-phi(ijk))
c    &              /(x(i+1,k)-x(i-1,k))
c    &              /(x(i+1,k)-x(i,k))
c       rhoc(ijk) = -eps0*rhoc(ijk) 

c         end if
c      end if
 100  continue

c     do 150 k=jnoz, ny-1
c      ijk = jj(iinj-1,k)

c      if( vf0(iinj-1,k) .lt. 0.5 ) then
c         rhoc(ijk) = 0.
c      else
c       rhoc(ijk) = ( epsi(iinj-1,k)*eps0*Er(ijk)
c    &              - epsi(iinj-2,k)*eps0*Er(ijk-ny) )
c    &              /( x(iinj-1,k) - x(iinj-2,k) )
c    &            + epsi(iinj-1,k)*eps0*Er(ijk)/x(iinj-1,k)
c    &            + ( epsi(iinj-1,k+1)*eps0*Ez(ijk+1)
c    &              - epsi(iinj-1,k)*eps0*Ez(ijk) )
c    &              /( y(iinj-1,k+1) - y(iinj-1,k) )

c         end if
c150  continue

c     do 200 k=2, ny-1
c       ijk = jj(1,k)
c       rhoc(ijk) = rhoc(ijk+ny)
c200  continue

      return
      end

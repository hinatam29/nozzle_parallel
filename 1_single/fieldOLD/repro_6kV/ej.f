
      subroutine ej
c ======================================================================
c
c     purpose   :  calculate Hall(Faraday) electric field
c                  and  Hall(Faraday) current density
c
c     CPU       :  CRAY,  MIPS,  PC-9801
c     date      :  July  25, 1996
c     programer :  H. Kobayashi
c
c ----------------------------------------------------------------------
      include 'subcom.inc'

      dimension Ez1(n), Er1(n)

c     echg = 1.602e-19
 
c
c: ----- calculate Ez of inner point -----
c
      alp = 1.e-5
      do 100 i=1, nx-1
      do 100 k=2, ny-2
         ijk = jj(i,k)

         yy1 = 0.5*(y(i,k) + y(i,k-1))
         yy3 = 0.5*(y(i,k+2) + y(i,k+1))

         dzs = 1./( yy3 - yy1 )
         Ez(ijk) = - ( phi(ijk+1) - phi(ijk-1) )*dzs

         if( abs(dphi(i,k)) .le. alp ) then
         yy1 = 0.5*(y(i,k) + y(i,k-1))
         yy3 = 0.5*(y(i,k+1) + y(i,k))
         dzs = 1./( yy3 - yy1 )
         Ez(ijk) = -( phi(ijk)-phi(ijk-1) )*dzs
         end if

c       if( Ez(ijk) .gt. 0. ) write(*,*) phi(ijk+1), phi(ijk-1)

  100 continue



c
c
c: ----- calculate Ez of inlet and outlet boundary -----
c
      do 110 i=1, nx-1
         ijk = jj(i,1)

         yy1 = 0.5*(y(i,2) + y(i,1))
         yy2 = 0.5*(y(i,3) + y(i,2))

         dzs = 1./( yy2 - yy1 )
         Ez(ijk) = - ( phi(ijk+1) - phi(ijk) )*dzs

         ijk = jj(i,ny-1)

         yy1 = 0.5*(y(i,ny-1) + y(i,ny-2))
         yy2 = 0.5*(y(i,ny) + y(i,ny-1))

         dzs = 1./( yy2 - yy1 )
         Ez(ijk) = - ( phi(ijk) - phi(ijk-1) )*dzs


  110 continue

c
c: ----- calculate Er of inner point -----
c
      do 200 i=2, nx-2
      do 200 k=1, ny-1
         ijk = jj(i,k)

         xx1 = 0.5*(x(i,k) + x(i-1,k))
         xx3 = 0.5*(x(i+2,k) + x(i+1,k))
 
         Er(ijk) = - ( phi(ijk+(ny-1)) - phi(ijk-(ny-1)) )
     &             / ( xx3 - xx1 )
  200 continue

c     do 205 k=jnoz, ny-1
c        ijk = jj(iinj,k)

c        xx1 = 0.5*(x(iinj,k) + x(iinj-1,k))
c        xx2 = 0.5*(x(iinj+1,k) + x(iinj,k))

c        Er(ijk) = - ( phi(ijk) - phi(ijk-(ny-1)) )
c    &             / ( xx2 - xx1 )

c 205 continue

c     do 208 k=2, jhol
c        ijk = jj(ihol,k)

c        xx1 = 0.5*(x(ihol,k) + x(ihol-1,k))
c        xx2 = 0.5*(x(ihol+1,k) + x(ihol,k))

c        Er(ijk) = - ( phi(ijk) - phi(ijk-(ny-1)) )
c    &             / ( xx2 - xx1 )

c 208 continue

c
c: ----- calculate Er of inlet and outlet boundary -----
c
      do 210 k=1, ny-1
         ijk = jj(1,k)
         xx1 = 0.5*(x(2,k) + x(1,k))
         xx2 = 0.5*(x(3,k) + x(2,k))

         Er(ijk) = -( phi(ijk+(ny-1)) - phi(ijk) )/( xx2 - xx1 )

         ijk = jj(nx-1,k)
         xx1 = 0.5*(x(nx-2,k) + x(nx-3,k))
         xx2 = 0.5*(x(nx-1,k) + x(nx-2,k))

         Er(ijk) = - ( phi(ijk) - phi(ijk-(ny-1)) )
     &              /( xx2 - xx1 )

  210 continue


c------------------ current ---------------------

c      xjz = 0.
c      crntz = 0.

c      do 300 k=icedge+1,nz-iaedge
c      do 310 i=1,nr-1
c        ijk = jj(i,k)
c        dnp = dnO2p(ijk)+dnN2p(ijk)+dnO4p(ijk)+dnN4p(ijk)
c     &      + dnH2Op(ijk)+dnH3Op(ijk)+dnCH4p(ijk)+dnCH3p(ijk)
c        dnng = dne(ijk)
c        dnng2 = dnO2ng(ijk)

c        pcrtz(ijk) =  echg*dnp*pmbl(ijk)*Ez(ijk)
c        ecrtz(ijk) =  echg*dnng*embl(ijk)*Ez(ijk)
c     &             +  echg*dnng2*pmbl(ijk)*Ez(ijk)

c        pcrtr(ijk) =  echg*dnp*pmbl(ijk)*Er(ijk)
c        ecrtr(ijk) =  echg*dnng*embl(ijk)*Er(ijk)
c     &             +  echg*dnng2*pmbl(ijk)*Er(ijk)

c        xjz = xjz+(pcrtz(ijk)+ecrtz(ijk))
c     &          *2.*pi*Rgrid(i,k)*(Rgrid(i+1,k)-Rgrid(i,k))
c310   continue
c      crntz = crntz + xjz
c      xjz = 0.
c  300 continue

c      crntz = crntz/(nz-iaedge-icedge)


      return
      end 


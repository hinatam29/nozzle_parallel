
      subroutine advVOF
      include 'subcom.inc'

c     dimension vfl(0:nn,0:nn), vf0p(0:nn,0:nn)
      dimension vfl(0:nn,0:nn)

c:    x-direction
c:    -----------------------

      do 100 j=2, ny-2
      do 100 i=2, nx-1
        uu = ux_sx(i,j)

        jc=j

        if( uu .ge. 0.) then
           ic=i-1
        else 
           ic=i
        end if

c       xnx = -vf0(ic-1,jc)+vf0(ic+1,jc)
        xnx = -dphi(ic-1,jc)+dphi(ic+1,jc)

        xi = uu*dt

        vf01 = vf0(ic,jc)
        vf02 = vf0(ic-1,jc)
        vf03 = vf0(ic+1,jc)

c       sdx1 = sdx(i,j)

        vfl(i,j) = THINC(vf01,vf02,vf03,xnx,xi,sdx)
c       vfl(i,j) = vf0(i,j)
100   continue

      do 150 j=2, ny-2
      do 150 i=2, nx-2
        vf0p(i,j) = vf0(i,j)-(vfl(i+1,j)-vfl(i,j))/sdx(i,j)
     &            + dt*vf0(i,j)*(vf_sx(i+1,j)-vf_sx(i,j))/sdx(i,j)
150   continue

      call bound

      do 170 i=2, nx-2
      do 170 j=2, ny-2
        vf0(i,j) = vf0p(i,j)
170   continue
      


c:    y-direction
c:    -----------------------

      do 200 i=2, nx-2
      do 200 j=2, ny-1
        vv = uy_sy(i,j)
       
        ic=i

        if( vv .ge. 0.) then
           jc=j-1
        else 
           jc=j
        end if

c       xny = -vf0(ic,jc-1)+vf0(ic,jc+1)
        xny = -dphi(ic,jc-1)+dphi(ic,jc+1)

        xi = vv*dt

        vf01 = vf0(ic,jc)
        vf02 = vf0(ic,jc-1)
        vf03 = vf0(ic,jc+1)

c       sdy1 = sdy(i,j)
        vfl(i,j) = THINC(vf01,vf02,vf03,xny,xi,sdy)
c       vfl(i,j) = vf0(i,j)

c       write(*,*) istp, vfl(i,j)
200   continue

      do 250 j=2, ny-2
      do 250 i=2, nx-2
        vf0p(i,j) = vf0(i,j)-(vfl(i,j+1)-vfl(i,j))/sdy(i,j)
     &            + dt*vf0(i,j)*(vf_sy(i,j+1)-vf_sy(i,j))/sdy(i,j)
250   continue


      do 270 i=2, nx-2
      do 270 j=2, ny-2
        vf0(i,j) = vf0p(i,j)
270   continue

      return
      end


      real function THINC(vf01i,vf02i,vf03i,xn1,xi1,dx1)
      include 'subcom.inc'

c     eps = 1.0e-10 
      eps = 1.0e-3
      sl = 20.0
      bet = 3.5

      if( abs(vf01i-1.) .le. eps .or. abs(vf01i) .le. eps ) then
         THINC = vf01i*xi1
      else if( abs(vf01i-vf02i) .le. eps .or. 
     &         abs(vf01i-vf03i) .le. eps) then
         THINC = vf01i*xi1
c     else if( vf02i .eq. 0. .or. vf03i .eq. 0. .or. 
c    &       (vf01i-vf02i)*(vf01i-vf03i) .ge. 0. )  then
      else if( (vf01i-vf02i)*(vf01i-vf03i) .ge. 0. )  then
         THINC = vf01i*xi1
      else
         if( xn1 .ge. 1.0 ) then
            sgm = 1.0
            alp = min(2.0,vf03i/(1.+eps))
         else
            sgm = -1.0
            alp = min(2.0,vf02i/(1.+eps))
         end if


         if( xi1 .ge. 0. ) then
            zz = 1.0
         else
            zz = 0.0
         end if

         a0 = max(sgm*2.*bet/alp*(vf01i-0.5*alp), -sl)
         a0 = min(a0,sl)

         xa = exp(a0)
         xb = exp(bet)
         xc = exp(-bet)

         xm = 0.5/bet*log(abs(xa-xb)/abs(xc-xa))

         xd = cosh(bet*(zz-xm))
         xe = cosh(bet*(zz-xi1/dx1-xm))

         THINC = 0.5*(sgm*alp/bet*log(abs(xd/xe))*dx1+alp*xi1)
      end if
     
      return
      end

  
          

      

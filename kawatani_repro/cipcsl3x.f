
      subroutine cipcsl3x(q_s,q0,qp,q0p)
      include 'subcom.inc'

      dimension q_s(0:nn,0:nn), q0(0:nn,0:nn) 
      dimension qp(0:nn,0:nn), q0p(0:nn,0:nn)
      dimension fl(0:nn,0:nn)
c     dimension fl(nn)

      do 100 j=2, ny-2
      do 200 i=2, nx-1
        uu = ux_sx(i,j)
c       uu = ux0(i,j)
c       uu = 0.

        if( uu .ge. 0.) then
        ip = i-1
        iv = i-1
        dd = -sdx(i,j)
        else
        ip = i+1
        iv = i
        dd = sdx(i,j)
        end if

        xi = -uu*dt

        gr1 = q_s(iv+1,j)-q_s(iv,j)
        gr2 = q0(iv,j)-q0(iv-1,j)
        gr3 = q0(iv+1,j)-q0(iv,j)

        if( abs(gr1) .le. TVB*sdx(i,j)*sdx(i,j) ) then
        grad = gr1
        else

        if( gr1*gr2 .le. 0. .or. gr1*gr3 .le. 0. ) then
        grad = 0.
        else
        grad = min(gr1,gr2,gr3)
        end if

        end if

        grad = grad*dd/sdx(i,j)


        c3 = 4.*(-q_s(i,j)+q_s(ip,j)-grad)/(dd*dd*dd) 
        c2 = 3.*(3.*q_s(i,j)-q_s(ip,j)-2.*q0(iv,j)+2.*grad)/(dd*dd) 
        c1 = 2.*(-3.*q_s(i,j)+3.*q0(iv,j)-grad)/dd
        c0 = q_s(i,j)

        qp(i,j) = c3*xi*xi*xi
     &          + c2*xi*xi  
     &          + c1*xi  
     &          + c0 
        fl(i,j) = -(c3/4.*xi*xi*xi*xi 
c       fl(i) = -(c3/4.*xi*xi*xi*xi 
     &          + c2/3.*xi*xi*xi 
     &          + c1/2.*xi*xi 
     &          + c0*xi) 

200   continue
100   continue

      do 250 j=2, ny-2
      fl(2,j) = ux_sx(2,j)*q_s(2,j)*dt
      fl(nx-1,j) = ux_sx(nx-1,j)*q_s(nx-1,j)*dt
250   continue


      do 300 j=2, ny-2
      do 300 i=2, nx-2
        q0p(i,j) = q0(i,j)-(fl(i+1,j)-fl(i,j))/sdx(i,j)
c       q0p(i,j) = q0(i,j)-(fl(i+1)-fl(i))/dx(i,j)
300   continue

      return
      end      

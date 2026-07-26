
      subroutine cipcsl3y(q_s,q0,qp,q0p)
      include 'subcom.inc'

      dimension q_s(0:nn,0:nn), q0(0:nn,0:nn)
      dimension qp(0:nn,0:nn), q0p(0:nn,0:nn)
      dimension fl(0:nn,0:nn)

      do 100 i=2, nx-2
      do 200 j=2, ny-1
        vv = uy_sy(i,j)

        if( vv .ge. 0.) then
        jp = j-1
        jv = j-1
c       dd = -dy(i,j)
        dd = -sdy(i,j)
        else
        jp = j+1
        jv = j
c       dd = dy(i,j)
        dd = sdy(i,j)
        end if

        xi = -vv*dt

        gr1 = q_s(i,jv+1)-q_s(i,jv)
        gr2 = q0(i,jv)-q0(i,jv-1)
        gr3 = q0(i,jv+1)-q0(i,jv)

c       if( abs(gr1) .le. TVB*dy(i,j)*dy(i,j) ) then
        if( abs(gr1) .le. TVB*sdy(i,j)*sdy(i,j) ) then
        grad = gr1
        else

        if( gr1*gr2 .le. 0. .or. gr1*gr3 .le. 0. ) then
        grad = 0.
        else
        grad = min(gr1,gr2,gr3)
        end if

        end if

c       grad = grad*dd/dy(i,j)
        grad = grad*dd/sdy(i,j)


        c3 = 4.*(-q_s(i,j)+q_s(i,jp)-grad)/(dd*dd*dd) 
        c2 = 3.*(3.*q_s(i,j)-q_s(i,jp)-2.*q0(i,jv)+2.*grad)/(dd*dd) 
        c1 = 2.*(-3.*q_s(i,j)+3.*q0(i,jv)-grad)/dd
        c0 = q_s(i,j)

        qp(i,j) = c3*xi*xi*xi
     &          + c2*xi*xi  
     &          + c1*xi  
     &          + c0 
        fl(i,j) = -(c3/4.*xi*xi*xi*xi 
     &          + c2/3.*xi*xi*xi 
     &          + c1/2.*xi*xi 
     &          + c0*xi) 

200   continue
100   continue

      do 250 i=2, nx-2
      fl(i,2) = uy_sy(i,2)*q_s(i,2)*dt
      fl(i,ny-1) = uy_sy(i,ny-1)*q_s(i,ny-1)*dt
250   continue

      do 300 i=2, nx-2
      do 300 j=2, ny-2
c       q0p(i,j) = q0(i,j)-(fl(i,j+1)-fl(i,j))/dy(i,j)
        q0p(i,j) = q0(i,j)-(fl(i,j+1)-fl(i,j))/sdy(i,j)
300   continue

      return
      end      

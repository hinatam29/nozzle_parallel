
      subroutine pfield_sor
      include 'subcom.inc'

      dimension prhs(0:nn,0:nn)

      omg = 1.5
      itr = 1. 

      itpmax = 5000

      pgerr = 0.
      pgt = 0. 

      do 100 i=2, nx-2
      do 100 j=2, ny-2
      prhs(i,j) =((-ux_sx(i,j)+ux_sx(i+1,j))/sdx(i,j)
     &        +   (-uy_sy(i,j)+uy_sy(i,j+1))/sdy(i,j))/dt
100   continue

150   continue

c     do 200 itrp=1, itpmax
        perr = 0.0
        do 300 i=2, nx-2
        do 300 j=2, ny-2
        pb = p0(i,j)

        rv = vf0(i,j)*rhol+(1.-vf0(i,j))*rhog
        rx0 = vf0(i-1,j)*rhol+(1.-vf0(i-1,j))*rhog
        rx1 = vf0(i+1,j)*rhol+(1.-vf0(i+1,j))*rhog
        ry0 = vf0(i,j-1)*rhol+(1.-vf0(i,j-1))*rhog
        ry1 = vf0(i,j+1)*rhol+(1.-vf0(i,j+1))*rhog

        r_sx0 = 0.5*(rv+rx0)
        r_sx1 = 0.5*(rv+rx1)
        r_sy0 = 0.5*(rv+ry0)
        r_sy1 = 0.5*(rv+ry1)

        dpm = -((1./r_sx0+1./r_sx1)/sdx(i,j)/sdx(i,j)
     &         +(1./r_sy0+1./r_sy1)/sdy(i,j)/sdy(i,j))
        dpxx = -(p0(i-1,j)/r_sx0+p0(i+1,j)/r_sx1)/sdx(i,j)/sdx(i,j)
        dpyy = -(p0(i,j-1)/r_sy0+p0(i,j+1)/r_sy1)/sdy(i,j)/sdy(i,j)

        p0(i,j) = omg*(dpxx+dpyy+prhs(i,j))/dpm
     &          + (1.-omg)*pb

        if(abs(p0(i,j)-pb) .ge. perr) then
           perr = abs(p0(i,j)-pb)
        end if
300     continue

        do 350 j=2,ny-2
          p0(1,j)=p0(2,j) 
          p0(nx-1,j)=p0(nx-2,j) 
          p0(nx,j)=p0(nx-1,j) 
350     continue

        do 370 i=2,nx-2
          p0(i,1)=p0(i,2) 
          p0(i,ny-1)=p0(i,ny-2) 
          p0(i,ny)=p0(i,ny-1) 
370     continue

c       write(*,*) perr

        if( perr .le. 1.0e-3 ) then
c       write(*,*) perr
           go to 1000
        end if
 200   continue

      go to 150

1000  return
      end 
                  


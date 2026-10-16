
       subroutine diffusion1
       include 'subcom.inc'

       double precision cg0, cg2, cg4, cg6, cg8
       double precision cgn0
       double precision cgnn(0:nn,0:nn)
c      double precision cvolume(0:nn,0:nn)
       double precision crate(0:nn,0:nn)
       double precision vcell(0:nn,0:nn)

       
       do 10 i=1, nx-1
       do 10 j=1, ny-1

       cgnn(i,j) = 0.
       
 10    continue


       do 100 i=2, nx-2
       do 100 j=2, ny-2

       if( i .ge. ihol  .and.  j .ge. jhol2  .and.  j .le. jhol1 )then
       go to 110
       end if

       sdx(i,j) = x(i+1,j) - x(i,j)
       sdy(i,j) = y(i,j+1) - y(i,j)


c      dd = 1.e-4

       cg0 = cg(i,j)
       cg2 = cg(i-1,j)
       cg4 = cg(i,j-1)
       cg6 = cg(i,j+1)
       cg8 = cg(i+1,j)


c      平面2次元(デカルト xy): 円筒の第1次微分項を除去
       cgn0 = cg0
     &      + dt*1.6e-5*( (cg2-2.*cg0+cg8)/sdx(i,j)/sdx(i,j)
     &                   +(cg4-2.*cg0+cg6)/sdy(i,j)/sdy(i,j) )

       cgnn(i,j) = cgn0

c      cgn(i,j) = cg(i,j)
c    &          + dt*1.6e-5*( (cg(i-1,j)-2.*cg(i,j)+cg(i+1,j))
c    &                    /sdx(i,j)/sdx(i,j)
c    &                   +1./sdx(i,j)
c    &                    *(cg(i+1,j)-cg(i,j))/sdx(i,j)
c    &                   +(cg(i,j-1)-2.*cg(i,j)+cg(i,j+1))
c    &                    /sdy(i,j)/sdy(i,j) )

c      write(*,*) cg(i-1,j)-2.*cg(i,j)+cg(i+1,j)
c      write(*,*) cgn(i,j)

 110   continue

 100   continue



       do 200 i=2, nx-2
       do 200 j=2, ny-2

       cg(i,j) = cgnn(i,j)

 200   continue

       xhol = 5.e-3
       yhol1 = 93.e-3
       yhol2 = 92.e-3

       do 300 i=2, nx-2
       do 300 j=2, ny-2

c      if( x(i,j) .gt. xhol  .and.  y(i,j) .gt. yhol2
c    &     .and.  y(i,j) .lt. yhol1 )then

       if( i .gt. ihol .and. j .gt. jhol2 .and. j .lt. jhol1 )then

c      go to 310
       cgg(i,j) = 0.445
c      cgg(i,j) = 1.

       else 

       cgg(i,j) = 0.445 - cg(i,j)
c      cgg(i,j) = 1. - cg(i,j)

       end if

 310   continue
 300   continue


       do 350 i=2, nx-2
       do 350 j=2, ny-2

       crate(i,j) = 0.
       cvolume(i,j) = 0.
       vcell(i,j) = 0.

 350   continue       


       do 400 i=2, nx-2
       do 400 j=2, ny-2

       pi = 4.*atan( 1. )

c      if( x(i,j) .ge. xhol  .and.  y(i,j) .ge. yhol2
c    &     .and.  y(i,j) .le. yhol1 )then

       if( i .gt. ihol .and. j .gt. jhol2 .and. j .lt. jhol1 )then

       go to 410
       
       else 

       crate(i,j) = cgg(i,j) / 0.445
c      crate(i,j) = cgg(i,j) / 1.
c      crate(i,j) = crate(i,j)*100

          if( crate(i,j) .ne. 0.   .and. crate(i,j) .ne. 1. )then
c         write(*,*) i,j,crate(i,j)
c         end if

c         円筒(2次元軸対称): 2πr*sdx*sdy
c         平面2次元(デカルト xy): セル面積 sdx*sdy
          vcell(i,j) = sdx(i,j)*sdy(i,j)
       
          cvolume(i,j) = 0.445 * crate(i,j) * vcell(i,j)
c         cvolume(i,j) = crate(i,j) * vcell(i,j)

          else 
       
          cvolume(i,j) = 0.

          end if
       end if

 410   continue      
 400   continue

       ctotal = 0.

       do 500 i=2, nx-2
       do 500 j=2, ny-2

       ctotal = ctotal + cvolume(i,j)

 500   continue       


c      write(*,*)  'volume integral concentration', ctotal       


       call bound1

       return

       end


  


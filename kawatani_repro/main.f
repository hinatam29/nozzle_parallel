
      program main
      include 'subcom.inc'

      gam = 1.402

c     dt = 5.0e-11
c     dt = 1.0e-10
c     dt = 2.0e-9
c     dt = 1.0e-8
c     dt = 1.0e-6
c     dt = 5.0e-7
      dt = 2.0e-7
      ftime = 11.e-3
c     ftime = 1.e2
      otime = 2.5e-2
      dtime = 1.e-2

      etime = 1.e-6

      pi = 4.0*atan(1.0)

      Re = 1000.

      gv = 9.8
c     gv = 0.

c     fss = 1. / ( 70.e3 )
c     fss = 70.e3
c     fss = fss/dt
c     fss = 3.
      fss = 7.
c     fss = 17.
c     fss = 14.
c     fss = 47.

c     write(*,*) fss

c     ipmx = 500.
c     ipmx = 2000.
      ipmx = 100000.
      ip = 1.

c     rhol = 797.88 
      rhol = 1027.
      rhog = 1.25

c     xmul = 0.8544e-3
c     xmul = 50e-3
      xmul = 10e-3
c     xmul = 85.44e-3
c     xmul = 128e-3
      xmug = 1.862e-5

c     svl = sqrt(1.18e8/789.)
      svl = 1500.
      svg = sqrt(1.402*287.*300.)

      sigl = 0.277
c     sigl = 100.e-6
      sigg = 2.5e-10 

c     vxmax = 0.
c     vymax = 0.


      TVB = 0. 
c     TEC = 1. 
      TEC = 0.95 
     

      iend = int(ftime/dt)
c     iend = 212505
      iout = int(otime/dt)
c     idata = int(dtime/dt)
c     idata = 1
c     idata = 10
c     idata = 100
      idata = 1000
      ikzoku = 500

      ielec = int(etime/dt)

      time = 0.
      istp = 1

      open(100,file='./flow.dat',status='unknown',blank='null')
      open(200,file='./kzoku.dat',status='unknown',blank='null')
      open(300,file='./spray.dat',status='unknown',blank='null')

      write(300,*) ' title = " SPRAY " '
      write(300,*) ' variables= "xp","yp","dp","clt","pout" '

      write(100,*) ' title = " FLOW FIELD  " '
c     write(100,*) ' variables= "x","y","vf0(-)","ux0(m/s)","uy0(m/s)",
c    & "p0(N/m2)","dphi(m)","fsvx","fsvy","phi","Er","Ez","rhoc0",
c    & "fer","fez"'
c     write(100,*) ' variables= "x","y","phi","Er","Ez","Etotal","cg"'
c     write(100,*) ' variables= "xp","yp","fcom"'
c     write(100,*) ' variables= "xp","yp","dp"'
      write(100,*) ' variables= "x","y","phi","cvolume","cgg"'
c     write(100,*) ' variables= "xp","yp","phi","cgg"'

c     iend = 10
c     iout = 1

      call geo

      call rancal
      call initial

c     open(300,file='./cont.dat',status='old')
c     call kinput

c     call bound1

 100   continue

c     call dfunc
c     call prop

c100   continue
c     call advect
c     call advVOF
c     call dfunc
c     call heviside
c     end if
c     call pfield_bicg
c     call proj

c     call dfunc
c     if( istp .eq. 1 .or. mod(real(istp),1000.) .eq. 0. ) then
c     call heviside
c     end if
c     call prop

c     call psearch
c     call diffusion
c     call difrho
      call viscous
c     call grav
c     call csf

c     irep = 1
110   continue
c     if( istp .eq. 1 .or. mod(real(istp),20.) .eq. 0. ) then
c     if( irep .le. 10 ) then
      call potential
      call ej
c     call charge
c     end if
c     call elec
c     irep = irep + 1
c     go to 110
c     end if

c     write(*,*) ip
      call eforce

      if( ip .eq. 1. ) then
      go to 130
      end if

      call intraction

 130  continue

      call pcount

      call diffusion_liq
      call diffusion_gas
      call diffusion1
      


      if( istp .eq. 1 .or. mod(real(istp),real(ikzoku)) .eq. 0. ) then
         call kzoku
      end if

      write(*,*) istp,time,ip,ctotal
      if( istp .eq. 1 .or. mod(real(istp),real(idata)) .eq. 0. ) then
c     write(*,*) istp, time  
c     if( istp .eq. iend ) then
c     if( time .eq. ftime ) then
c     write(100,*) 'zone t="n= ',istp,'", i=',nx-3,',j=',ny-3,',f=point'
      write(100,*) 'zone t="n= ',istp,'", i=',nx-1,',j=',ny-1,',f=point'
c     write(100,*) 'zone t="n= ',istp,'", f=point'
         do 150 j=1, ny-1
         do 160 i=1, nx-1
         xx = 0.5*(x(i,j)+x(i+1,j))
         yy = 0.5*(y(i,j)+y(i,j+1))
         ij = jj(i,j)         

         Etotal(ij) = sqrt( (Er(ij)*Er(ij)+Ez(ij)*Ez(ij)) )

c        write(100,1000) xx, yy, vf0(i,j), ux0(i,j), uy0(i,j), p0(i,j),
c    &                   dphi(i,j),fsvx(i,j),fsvy(i,j),phi(ij),
c    &                   dphi(i,j),gradx(i,j),grady(i,j),phi(ij),
c    &                   Er(ij),Ez(ij),rhoc0(i,j),fer(ij),fez(ij) 
c        write(100,*) xx, yy, vf0(i,j)
c        write(100,*) xx, yy, p0(i,j)
         write(100,1000) xx, yy, phi(ij), cvolume(i,j), cgg(i,j)
c        write(100,1000) xx, yy, phi(ij),Er(ij),Ez(ij),
c    &                   Etotal(ij),cg(i,j)
c        write(100,*) xx, yy, phi(ij),Etotal(ij),cg(i,j)
c    &                   Ernex(ij), Eznex(ij) 
 160      continue
c        write(100,*)
 150      continue

c        --- spray (particle) output for visualization ---
         write(300,*) 'zone t="n= ',istp,'", i=',ip,',f=point'
         do 170 k=1, ip
         write(300,1000) xp(k), yp(k), dp(k), clt(k), pout(k)
 170      continue

c        --- 各出力ごとにディスクへ確実に書き出す ---
c        （途中で止めても完成フレームがファイルに残る）
         call flush(100)
         call flush(300)

1000  format(15E20.8)
      end if

c     vxmax = 0.
c     vymax = 0.

c     ux = 0.
c     uy = 0.
c     cfl = 0.3
c     if( mod(real(istp),100.) .eq. 0. ) then

c     do 200  i=1, nx
c     do 200  j=1, ny

c     ux = sqrt(ux_sx(i,j)*ux_sx(i,j)+
c    &          ux_sy(i,j)*ux_sy(i,j))
       
c     uy = sqrt(uy_sx(i,j)*uy_sx(i,j)+
c    &          uy_sy(i,j)*uy_sy(i,j))
      
c     ux = ux_sx(i,j)
c     uy = uy_sy(i,j)
      
c     uxsdx = ux/sdx(i,j)
c     uysdy = uy/sdy(i,j) 

c     if( vxmax .lt. uxsdx) then
c     vxmax = uxsdx
c     end if
      
c     if( vymax .lt. uysdy) then
c     vymax = uysdy
c     end if
c00   continue

c     dt = min( cfl/vxmax, cfl/vymax)
c     if( dt .gt. 1.e-6) dt = 1.e-6
        
c     end if

c     write(*,*) real(istp), fss

      if( ip .eq. ipmx )then
         go to 200
      end if
   
      if( mod(real(istp),real(fss)) .eq. 0. )then
         ip = ip + 1
      end if

 200  continue

      if( istp .le. iend ) then
      time = time + dt
      istp = istp + 1
      go to 100
      end if

 
      stop
      end      

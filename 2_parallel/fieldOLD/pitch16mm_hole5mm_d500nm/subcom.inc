      implicit real(a-b,d-h,o-z)
c     implicit real(a-h,o-z)
      implicit double precision(c)

c     parameter( nx = 90, ny = 70,  nn = 200 ) 
c     parameter( nx = 30, ny = 120,  nn = 400, n= nx*ny ) 
      parameter( nx = 500, ny = 500,  nn = 550, n= nx*ny, m=100000 ) 
c     parameter( nx = 250, ny = 500,  nn = 550, n= nx*ny, m=100000 ) 
c     parameter( nx = 50, ny = 360,  nn = 390, n= nx*ny ) 
c     parameter( nx = 50, ny = 160,  nn = 390, n= nx*ny ) 
c     parameter( nx = 15, ny = 15,  nn = 200 ) 
c     parameter( nx = 153, ny = 153,  nn = 200 ) 
      parameter( nx1 = nx+1, ny1 = ny+1 ) 

      common /arry0/ time,pi,Re,istp,iinj,jnoz,inoz,ihol,jhol1,jhol2,
     &               inoz2,ihol2
      common /arry1/ dt,gam,sdx(0:nn,0:nn),sdy(0:nn,0:nn),
     &               sdx1(0:nn,0:nn),sdy1(0:nn,0:nn),
c    &               iinj(0:nn),inoz(0:nn),
     &               x(0:nn,0:nn),TVB,TEC,
     &               y(0:nn,0:nn),jj(nn,nn)
      common /arry2/ rho_sx(0:nn,0:nn),ux_sx(0:nn,0:nn),
     &               uy_sx(0:nn,0:nn),
     &               trho_sx(0:nn,0:nn),
     &               pp_sx(0:nn,0:nn),
     &               rhoc0_sx(0:nn,0:nn),
     &               pg_sx(0:nn,0:nn),ene_sx(0:nn,0:nn),
     &               vf_sx(0:nn,0:nn),
     &               dp_sx(0:nn,0:nn),
     &               ene_sy(0:nn,0:nn),
     &               rho_sy(0:nn,0:nn),ux_sy(0:nn,0:nn),
     &               uy_sy(0:nn,0:nn),
     &               trho_sy(0:nn,0:nn),
     &               pp_sy(0:nn,0:nn),
     &               rhoc0_sy(0:nn,0:nn),
     &               pg_sy(0:nn,0:nn),
     &               vf_sy(0:nn,0:nn),
     &               dp_sy(0:nn,0:nn),
     &               xmdt_sx(0:nn,0:nn),ymdt_sx(0:nn,0:nn), 
     &               xmdt_sy(0:nn,0:nn),ymdt_sy(0:nn,0:nn), 
     &               upv(0:nn,0:nn),vpv(0:nn,0:nn),
     &               trho(0:nn,0:nn),tmu(0:nn,0:nn),tsv(0:nn,0:nn),
     &               tsig(0:nn,0:nn)
      common /arry3/ rho0(0:nn,0:nn),xmdt0(0:nn,0:nn),ymdt0(0:nn,0:nn),
     &               ene0(0:nn,0:nn),
     &               p0(0:nn,0:nn),rhopx(0:nn,0:nn),rhopy(0:nn,0:nn),
     &               rho0p(0:nn,0:nn),
     &               ux0(0:nn,0:nn),uy0(0:nn,0:nn),
     &               trho0(0:nn,0:nn),pp0(0:nn,0:nn),
     &               uxpx(0:nn,0:nn),uxpy(0:nn,0:nn),ux0p(0:nn,0:nn),
     &               uypx(0:nn,0:nn),uypy(0:nn,0:nn),uy0p(0:nn,0:nn),
     &               vfpx(0:nn,0:nn),vfpy(0:nn,0:nn),vf0p(0:nn,0:nn),
     &               trhpx(0:nn,0:nn),trhpy(0:nn,0:nn),trh0p(0:nn,0:nn),
     &               pppx(0:nn,0:nn),pppy(0:nn,0:nn),pp0p(0:nn,0:nn),
     &               rhoc0(0:nn,0:nn),
     &               rhocpx(0:nn,0:nn),rhocpy(0:nn,0:nn),
     &               rhoc0p(0:nn,0:nn),
c    &               dppx(0:nn,0:nn),dppy(0:nn,0:nn),
c    &               dp0p(0:nn,0:nn),
     &               xmdtpx(0:nn,0:nn),xmdtpy(0:nn,0:nn),
     &               xmdt0p(0:nn,0:nn),
     &               ymdtpx(0:nn,0:nn),ymdtpy(0:nn,0:nn),
     &               ymdt0p(0:nn,0:nn),
     &               enepx(0:nn,0:nn),enepy(0:nn,0:nn),ene0p(0:nn,0:nn)
      common /arry4/ sv(0:nn,0:nn)
      common /arry5/ gv,rhol,rhog,xmul,xmug,svl,svg,uin,sigg,sigl
      common /arry6/ vf0(0:nn,0:nn),dphi(0:nn,0:nn),
     &               fsvx(nx,ny),fsvy(nx,ny), 
     &               gradx(nx,ny),grady(nx,ny)
      common /arry7/ phi(n),epsi(0:nn,0:nn),sig(0:nn,0:nn),
c    &               rhoc(n),
     &               ele(0:nn,0:nn),
     &               Er(n),Ez(n),
c    &               fer(n),fez(n),
     &               crv(nx,ny)
c    &               alph(0:nn,0:nn)
      common /arry8/ ip,ipmx,dp(m),xp(m),yp(m),uxp(m),uyp(m),pout(m),
     &               Etotal(n),fcom(m)
c    &               cg(0:nn,0:nn),cgn(0:nn,0:nn)
c    &               xnex(0:nn,0:nn),ynex(0:nn,0:nn),
c    &               Ernex(n),Eznex(n)
      common /arry9/ cl(m),clt(m),cltmx(m),cvolume(0:nn,0:nn),ctotal,
     &               cg(0:nn,0:nn),cgn(0:nn,0:nn),cgg(0:nn,0:nn),
     &               fvisx(m),fvisy(m),fer(m),fez(m),ftxt(m),ftyt(m)

  

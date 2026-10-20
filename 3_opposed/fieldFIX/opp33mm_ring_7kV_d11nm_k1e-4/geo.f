
      subroutine geo
      include 'subcom.inc'

c     xlen = 0.165e-3 
      xlen = 1.e-1 
c     xlen = 50.e-3 
c     xlen = 20.e-3 
c     ylen = 1.43e-3
c     ylen = 0.53e-3
      ylen = 1.e-1
      dx = xlen/real(nx)
      dy = ylen/real(ny)

c     ================= 形状の設定（ここだけ変えれば形が変わる）=================
c     nopp  : ノズルの本数（1: 上から1本 / 2: 上下から向かい合わせ）
c     xnoz  : ノズルの半径方向の幅（川谷さんと同じ 0.4 mm）
c     xhol  : リング電極の内半径（実験: 内径 10 mm）
c     xring : リング電極の外半径（実験: 外径 12 mm）
c     yrc   : リング電極の厚さ方向の中心（厚さ 1 mm、格子点にそろえる）
c     ytip1 : 上のノズルの先端の高さ
c     ytip2 : 下のノズルの先端の高さ（nopp = 2 のとき）
      nopp  = 2
      xnoz  = 0.4e-3
      xhol  = 5.0e-3
      xring = 6.0e-3
      yrc   = 50.0e-3
      ytip1 = 66.5e-3
      ytip2 = 33.5e-3
c     ==========================================================================
c     xnoz = 1.e-3
c     ynoz = 1.e-3   ! 旧: 上端から 1 mm（今は ytip1 で決める）
      ynoz = ylen - ytip1
c     ynoz = 1.5e-3
c     xhol = 5.e-3   ! 上の設定で決める
c     yhol1 = 96.e-3
      yhol1 = 93.e-3
c     yhol2 = 95.e-3
      yhol2 = 92.e-3

      do 50 j=1, ny-1
      do 50 i=1, nx-1
        ijk = (i-1)*(ny-1)+j
        jj(i,j) = ijk
50    continue

c     do 55 i=1, nx
c     do 55 j=1, ny
c       jj(i,j) = ijk
c       write(*,*) i,j, ijk
c55    continue

c     open(10,file='./grid79.dat',status='old')

      do 100 j=1, ny
      do 100 i=1, nx
c     do 100 j=1, ny
c       read(10,*) x(i,j), y(i,j)
c       x(i,j) = 0.5*x(i,j)
c       y(i,j) = 0.37*y(i,j)
        x(i,j) = 1.e-15 + dx*(i-1)
        y(i,j) = 0. + dy*(j-1)
 100   continue

      do 120 j=1, ny
        x(0,j) = x(1,j) - (x(2,j)-x(1,j))
        y(0,j) = y(1,j) - (y(2,j)-y(1,j))
        x(nx+1,j) = x(nx,j) + (x(nx,j)-x(nx-1,j))
        y(nx+1,j) = y(nx,j) + (y(nx,j)-y(nx-1,j))
120   continue

      do 130 i=1, nx
        x(i,0) = x(i,1) - (x(i,2)-x(i,1))
        y(i,0) = y(i,1) - (y(i,2)-y(i,1))
        x(i,ny+1) = x(i,ny) + (x(i,ny)-x(i,ny-1))
        y(i,ny+1) = y(i,ny) + (y(i,ny)-y(i,ny-1))
130   continue

      do 300 i=0, nx
      do 300 j=0, ny+1
        sdx(i,j) = x(i+1,j)-x(i,j) 
c       sdx(i,j) = x(i,j)-x(i-1,j) 
c       sdx(i,j) = x(1,j)-x(0,j) 
c       write(*,*) x(i+1,j)-x(i,j)-sdx(i,j)
 300   continue

      do 305 i=1, nx
      do 305 j=1, ny-1
        sdx1(i,j) = 0.5*(x(i+1,j)+x(i-1,j))
 305   continue

      do 400 i=0, nx+1
      do 400 j=0, ny
        sdy(i,j) = y(i,j+1)-y(i,j)
c       sdy(i,j) = y(i,2)-y(i,1) 
c       write(*,*) y(i,j+1)-y(i,j), sdy(i,j)
400   continue

      do 405 i=1, nx-1
      do 405 j=1, ny
       sdy1(i,j) = 0.5*(y(i,j+1)+y(i,j-1))
405   continue

      do 500 j=0, ny+1
       sdx(nx+1,j) = x(nx+1,j)-x(nx,j)
500   continue

      do 600 i=0, nx+1
       sdy(i,ny+1) = y(i,ny+1)-y(i,ny)
600   continue

      do 700 i=1, nx
        if( x(i,1) .lt. xnoz ) then
            inoz = i-1
        end if
 700   continue

      do 800 j=1, ny
        if( y(1,j) .lt. ylen-ynoz) then
            jnoz = j
         end if
 800   continue

c     --- リング電極の格子番号（四捨五入。旧: 穴の縁が 0.4 mm 内側にずれていた）---
      ihol  = nint( xhol/dx ) + 1
      iring = nint( xring/dx ) + 1
      jhol2 = nint( (yrc-0.4e-3)/dy ) + 1
      jhol1 = nint( (yrc+0.4e-3)/dy ) + 1
c     --- 下のノズル（nopp = 2）: 先端より上で最初の格子点まで ---
      jnozb = 0
      do 950 j=1, ny
        if( y(1,j) .le. ytip2 ) jnozb = j + 1
 950  continue
      xnoz0 = xnoz
      xhol0 = xhol
      xring0 = xring
      yrc0 = yrc



c     write(*,*) inoz,jnoz,ihol,jhol1,jhol2 

      write(*,*) 'nozzle: ', x(inoz,ny)
      write(*,*) 'counter electrode: ', x(ihol,jhol1)
      write(*,*) 'nopp=',nopp,' nozzle1 from y=',y(1,jnoz),
     &           ' ring r=',x(ihol,1),'-',x(iring,1),
     &           ' z=',y(1,jhol2),'-',y(1,jhol1)
      if( nopp .eq. 2 ) write(*,*) 'nozzle2 up to y=',y(1,jnozb)
c     可視化（visualize.py）用に形状を書き出す
      open(77,file='./geom.txt',status='unknown')
      write(77,*) nopp, xlen, ylen, xnoz, xhol, xring,
     &            yrc-0.5e-3, yrc+0.5e-3, ytip1, ytip2
      close(77)

      return
      end


       subroutine diffusion_gas
       include 'subcom.inc'
c     =================================================================
c     液滴が吸った二酸化炭素を気体から減らす（v3 版、2026-10-22）
c     旧版からの変更:
c       - 減らす格子点の選び方を、液滴の側（diffusion_liq.f）と同じにした
c         （旧版は2つ目の条件の書き間違いで、液滴と違う格子点を選ぶことがあった）
c       - 鏡面の列（i=1, i=nx-1）の格子は幅が半分なので、体積も半分にする
c         （diffusion1.f と合わせて、気体の減り方＝液滴の吸収量になる）
c     =================================================================
       real sp(100000), vp(100000)
       real inex(100000), jnex(100000)
       real l1(100000), l2(100000), l3(100000), l4(100000)
       integer ii, jjj

       pi = 4.*atan(1.)
       abco2 = 1.0e-4

       do 10 k=1, ip
       inex(k) = nx-2
       jnex(k) = ny-2
       do 11 i=1, nx-1
       if( x(i,1) .gt. xp(k) )then
         inex(k) = i-1
         go to 12
       end if
 11    continue
 12    continue
       do 13 j=1, ny-1
       if( y(1,j) .gt. yp(k) )then
         jnex(k) = j-1
         go to 14
       end if
 13    continue
 14    continue
       if( inex(k) .lt. 1. ) inex(k) = 1.
       if( inex(k) .gt. nx-2 ) inex(k) = nx-2
       if( jnex(k) .lt. 1. ) jnex(k) = 1.
       if( jnex(k) .gt. ny-2 ) jnex(k) = ny-2
 10    continue

       do 20 k=1, ip
       i = inex(k)
       j = jnex(k)
       l1(k) = sqrt( (xp(k)-x(i,j))**2 + (yp(k)-y(i,j))**2 )
       l2(k) = sqrt( (x(i+1,j)-xp(k))**2 + (yp(k)-y(i+1,j))**2 )
       l3(k) = sqrt( (x(i+1,j+1)-xp(k))**2 + (y(i+1,j+1)-yp(k))**2 )
       l4(k) = sqrt( (xp(k)-x(i,j+1))**2 + (y(i,j+1)-yp(k))**2 )
       sp(k) = pi*dp(k)*dp(k)
       vp(k) = pi*dp(k)*dp(k)*dp(k)/6.
       if( clt(k) .gt. cltmx(k) ) go to 20

c      diffusion_liq.f と同じ順番・同じ条件で格子点を選ぶ
       if( l1(k) .le. l2(k)  .and.  l1(k) .le. l4(k) )then
         ii = i
         jjj = j
       else if( l2(k) .le. l1(k)  .and.  l2(k) .le. l3(k) )then
         ii = i+1
         jjj = j
       else if( l3(k) .lt. l4(k)  .and.  l3(k) .lt. l2(k) )then
         ii = i+1
         jjj = j+1
       else if( l4(k) .lt. l3(k)  .and.  l4(k) .lt. l1(k) )then
         ii = i
         jjj = j+1
       else
         go to 20
       end if

       vol = sdx(ii,jjj)*sdy(ii,jjj)*1.0
       if( ii .eq. 1 .or. ii .eq. nx-1 ) vol = 0.5*vol
       cg(ii,jjj) = cg(ii,jjj) - dt*abco2*cg(ii,jjj)*sp(k)/vol
 20    continue

       call bound1

       return
       end

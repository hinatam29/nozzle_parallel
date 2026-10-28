       subroutine diffusion1
       include 'subcom.inc'
c     =================================================================
c     気体（二酸化炭素）の拡散（v3 版、2026-10-22）
c     旧版からの変更:
c       - 鏡面の列（左のノズルの軸 i=1 と右のノズルの軸 i=nx-1）も
c         拡散を計算する。鏡の条件（外側の値＝内側の隣の値）を使う。
c         旧版はこの2列が計算されず、液滴が吸うと濃度が下がったまま
c         周りの気体を吸い込む「穴」になり、気体の減り方が液滴の吸収量の
c         約 1.5 倍になっていた
c       - ctotal（気体側で減った量）にこの2列も数える（格子の幅は半分）
c       - 格子のループを OpenMP で並列化
c     =================================================================
       double precision cg0, cg2, cg4, cg6, cg8
       double precision cgn0
       double precision cgnn(0:nn,0:nn)
       double precision crate(0:nn,0:nn)
       double precision vcell(0:nn,0:nn)

       do 10 i=1, nx-1
       do 10 j=1, ny-1
       cgnn(i,j) = 0.
 10    continue

c$omp parallel do private(i,j,cg0,cg2,cg4,cg6,cg8,cgn0)
       do 100 i=1, nx-1
       do 105 j=2, ny-2
       if( i .ge. ihol  .and.  i .le. ihol2  .and.
     &     j .ge. jhol2  .and.  j .le. jhol1 ) go to 105
       sdx(i,j) = x(i+1,j) - x(i,j)
       sdy(i,j) = y(i,j+1) - y(i,j)
       cg0 = cg(i,j)
       if( i .eq. 1 )then
         cg2 = cg(2,j)
       else
         cg2 = cg(i-1,j)
       end if
       if( i .eq. nx-1 )then
         cg8 = cg(nx-2,j)
       else
         cg8 = cg(i+1,j)
       end if
       cg4 = cg(i,j-1)
       cg6 = cg(i,j+1)
       cgn0 = cg0
     &      + dt*1.6e-5*( (cg2-2.*cg0+cg8)/sdx(i,j)/sdx(i,j)
     &                   +(cg4-2.*cg0+cg6)/sdy(i,j)/sdy(i,j) )
       cgnn(i,j) = cgn0
 105   continue
 100   continue
c$omp end parallel do

       do 200 i=1, nx-1
       do 200 j=2, ny-2
       cg(i,j) = cgnn(i,j)
 200   continue

       do 300 i=1, nx-1
       do 300 j=2, ny-2
       if( i .gt. ihol .and. i .lt. ihol2 .and.
     &     j .gt. jhol2 .and. j .lt. jhol1 )then
         cgg(i,j) = 0.445
       else
         cgg(i,j) = 0.445 - cg(i,j)
       end if
 300   continue

       do 350 i=1, nx-1
       do 350 j=2, ny-2
       crate(i,j) = 0.
       cvolume(i,j) = 0.
       vcell(i,j) = 0.
 350   continue

       do 400 i=1, nx-1
       do 400 j=2, ny-2
       if( i .gt. ihol .and. i .lt. ihol2 .and.
     &     j .gt. jhol2 .and. j .lt. jhol1 ) go to 400
       crate(i,j) = cgg(i,j) / 0.445
       if( crate(i,j) .ne. 0.   .and. crate(i,j) .ne. 1. )then
c         平面2次元(デカルト xy): セル面積。鏡面の列は幅が半分
          vcell(i,j) = sdx(i,j)*sdy(i,j)
          if( i .eq. 1 .or. i .eq. nx-1 ) vcell(i,j) = 0.5*vcell(i,j)
          cvolume(i,j) = 0.445 * crate(i,j) * vcell(i,j)
       else
          cvolume(i,j) = 0.
       end if
 400   continue

       ctotal = 0.
       do 500 i=1, nx-1
       do 500 j=2, ny-2
       ctotal = ctotal + cvolume(i,j)
 500   continue

       call bound1

       return
       end

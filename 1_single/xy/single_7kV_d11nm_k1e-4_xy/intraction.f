      subroutine intraction
      include 'subcom.inc'
c     =================================================================
c     液滴どうしのクーロン力（高速化版、2026-10-22、patch_fast.py）
c     物理は旧版と同じ:
c       (1) 直接の力
c       (2) 左の鏡面 x=0 の鏡像の力（自分自身の鏡像は含めない＝旧版と同じ）
c       (3) 右の鏡面 x=xlen の鏡像の力（自分自身の鏡像は含めない）
c     変えたところ（結果は旧版と一致することを確認済み）:
c       - 二重ループを1つにまとめた
c       - 飛行中の液滴（pout=0）だけを集めて計算する
c         （付着などした液滴は電荷 0 として扱われ、力も使われないため）
c       - 内側のループをベクトル化（omp simd）した
c     元のファイルは intraction.f.orig
c     =================================================================
      real xa(100000), ya(100000), da(100000), qa(100000)
      integer ia(100000)
      real fxs, fys, fx1, fy1, fx2, fy2
      real rx, ry, r, r3, dd, w, xk, yk, dk

      eps0 = 8.8542e-12
      fs = 42.9e-3
      pi = 4.*atan(1.)
      xlen = 1.e-1      ! xy の Single: 領域の横幅 100 mm

c     --- 飛行中の液滴を集める ---
      na = 0
      do 2 k=1, ip
         ftxt(k) = 0.
         ftyt(k) = 0.
         if( pout(k) .eq. 0. )then
            na = na + 1
            ia(na) = k
            xa(na) = xp(k)
            ya(na) = yp(k)
            da(na) = dp(k)
            qa(na) = pi*sqrt( 8.*eps0*fs*dp(k)*dp(k)*dp(k) )
         end if
 2    continue

c$omp parallel do schedule(static)
c$omp& private(kk,ll,xk,yk,dk,fxs,fys,fx1,fy1,fx2,fy2,
c$omp&         rx,ry,r,r3,dd,w)
      do 10 kk=1, na
         xk = xa(kk)
         yk = ya(kk)
         dk = da(kk)
         fxs = 0.
         fys = 0.
         fx1 = 0.
         fy1 = 0.
         fx2 = 0.
         fy2 = 0.
c$omp simd reduction(+:fxs,fys,fx1,fy1,fx2,fy2)
c$omp& private(rx,ry,r,r3,dd,w)
         do 20 ll=1, na
            dd = dk + da(ll)
            ry = yk - ya(ll)
c           自分自身（ll=kk）は重み 0
            w = merge( 0., 1., ll .eq. kk )
c           (1) 直接
            rx = xk - xa(ll)
            r  = sqrt( max( rx*rx + ry*ry, dd*dd ) )
            r3 = w*qa(ll)/(r*r*r)
            fxs = fxs + r3*rx
            fys = fys + r3*ry
c           (2) 左の鏡面 x=0 の鏡像
            rx = xk + xa(ll)
            r  = sqrt( max( rx*rx + ry*ry, dd*dd ) )
            r3 = w*qa(ll)/(r*r*r)
            fx1 = fx1 + r3*rx
            fy1 = fy1 + r3*ry
c           (3) 右の鏡面 x=xlen の鏡像
            rx = xk + xa(ll) - 2.*xlen
            r  = sqrt( max( rx*rx + ry*ry, dd*dd ) )
            r3 = w*qa(ll)/(r*r*r)
            fx2 = fx2 + r3*rx
            fy2 = fy2 + r3*ry
 20      continue
         ftxt(ia(kk)) = fxs + fx1 + fx2
         ftyt(ia(kk)) = fys + fy1 + fy2
 10   continue
c$omp end parallel do

      return
      end

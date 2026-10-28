      subroutine intraction
      include 'subcom.inc'
c     =================================================================
c     液滴どうしのクーロン力（v3 版、2026-10-22、patch_v3.py）
c     物理は旧版と同じ:
c       (1) 直接の力
c       (2) 左の鏡面 x=0 の鏡像の力
c       (3) 右の鏡面 x=xlen の鏡像の力
c     v3 の変更: 横に並べて出した液滴（軸から 0.025, 0.05 mm）は、自分自身の
c       鏡像の力も入れる。鏡像は「軸の反対側に並べて出した同じノズルの液滴」
c       にあたり、本当は反発を受ける（旧版は含めていなかった）。
c       軸の上で出した液滴（番号 k で (k-1) mod 3 = 0）は1個だけの液滴なので、
c       自分の鏡像は入れない（入れると、軸からわずかにずれたときに自分の
c       鏡像から非常に強く押し出されてしまう）
c     高速化（patch_fast.py と同じ）:
c       - 二重ループを1つにまとめた
c       - 飛行中の液滴（pout=0）だけを集めて計算する
c         （付着などした液滴は電荷 0 として扱われ、力も使われないため）
c       - 内側のループをベクトル化（omp simd）した
c     元のファイルは intraction.f.orig
c     （patch_v3）
c     =================================================================
      real xa(100000), ya(100000), da(100000), qa(100000), sa(100000)
      integer ia(100000)
      real fxs, fys, fx1, fy1, fx2, fy2
      real rx, ry, r, r3, dd, w, wi, xk, yk, dk, ws

      eps0 = 8.8542e-12
      fs = 42.9e-3
      pi = 4.*atan(1.)
      xlen = 3.2e-2     ! ノズル間隔 32mm

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
c           --- v3r: 帯電量を限界帯電量の 0.8 倍にする ---
            qa(na) = 0.8 * qa(na)
c           自分の鏡像の重み：軸の上で出た液滴は 0、横に出た液滴は 1
            if( mod(k-1,3) .eq. 0 )then
               sa(na) = 0.
            else
               sa(na) = 1.
            end if
         end if
 2    continue

c$omp parallel do schedule(static)
c$omp& private(kk,ll,xk,yk,dk,fxs,fys,fx1,fy1,fx2,fy2,
c$omp&         rx,ry,r,r3,dd,w,wi,ws)
      do 10 kk=1, na
         xk = xa(kk)
         yk = ya(kk)
         dk = da(kk)
         ws = sa(kk)
         fxs = 0.
         fys = 0.
         fx1 = 0.
         fy1 = 0.
         fx2 = 0.
         fy2 = 0.
c$omp simd reduction(+:fxs,fys,fx1,fy1,fx2,fy2)
c$omp& private(rx,ry,r,r3,dd,w,wi)
         do 20 ll=1, na
            dd = dk + da(ll)
            ry = yk - ya(ll)
c           自分自身（ll=kk）は重み 0
            w = merge( 0., 1., ll .eq. kk )
c           鏡像の重み（自分自身の鏡像は ws）
            wi = merge( ws, 1., ll .eq. kk )
c           (1) 直接
            rx = xk - xa(ll)
            r  = sqrt( max( rx*rx + ry*ry, dd*dd ) )
            r3 = w*qa(ll)/(r*r*r)
            fxs = fxs + r3*rx
            fys = fys + r3*ry
c           (2) 左の鏡面 x=0 の鏡像
            rx = xk + xa(ll)
            r  = sqrt( max( rx*rx + ry*ry, dd*dd ) )
            r3 = wi*qa(ll)/(r*r*r)
            fx1 = fx1 + r3*rx
            fy1 = fy1 + r3*ry
c           (3) 右の鏡面 x=xlen の鏡像
            rx = xk + xa(ll) - 2.*xlen
            r  = sqrt( max( rx*rx + ry*ry, dd*dd ) )
            r3 = wi*qa(ll)/(r*r*r)
            fx2 = fx2 + r3*rx
            fy2 = fy2 + r3*ry
 20      continue
         ftxt(ia(kk)) = fxs + fx1 + fx2
         ftyt(ia(kk)) = fys + fy1 + fy2
 10   continue
c$omp end parallel do

      return
      end

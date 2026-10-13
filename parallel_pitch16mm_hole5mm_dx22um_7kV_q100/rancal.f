
      subroutine rancal
      include 'subcom.inc'

      do 10 k = 1, ipmx

c     --- uniform random number ---

      x1 = rand()
      x2 = rand()

c     --- standard normal distribution ---


      sig1 = 1.
c     sig1 = 1.2
c     myu1 = 1.
c     myu1 = 4.
      myu1 = 5.
c     myu1 = 9.
c     myu1 = 7.
c     myu1 = 6.

      pi = 4.*atan(1.)
      xnorm1 = sqrt( -2.*log(x1) )*cos( 2.*pi*x2 )
      xnorm2 = sqrt( -2.*log(x1) )*sin( 2.*pi*x2 )

c     xnorm1 = exp( xnorm1 )

      xnorm1 = xnorm1*sig1 + myu1
      xnorm1 = xnorm1*1.e-7
c     xnorm1 = xnorm1*1.e-6


      dp(k) = xnorm1

c     dp(1) = dp(1)*0.1
c     dp(2) = 1.e-6

c     if( dp(k) .gt. 1.e-4 )then
c     dp(k) = dp(k)*1.e-1
c     end if

c     if( dp(k) .lt. 9.e-6 )then
c     dp(k) = dp(k)*10.
c     end if

      if( dp(k) .le. 0. )then
      dp(k) = -dp(k)
      end if

c     write(*,*) k, dp(k)
 10   continue


c     dp(1) = 4.e-6
c     dp(1) = 9.e-7
      dp(1) = 5.e-7
c     dp(1) = 7.e-7
c     dp(1) = 5.e-6
c     dp(1) = 6.e-6
c     write(*,*) dp(1)

      return
      end




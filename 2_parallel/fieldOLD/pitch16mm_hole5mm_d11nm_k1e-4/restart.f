c ======================================================================
c  途中から再開するための保存・読み込み（restart.dat）
c   - wrestart : 今の状態を restart.dat に保存する
c                （一時ファイルに書いてから名前を変えるので、保存中に
c                  止まっても前回の restart.dat は壊れない）
c   - rrestart : restart.dat から状態を読み込む
c   - prand    : 液滴を出すタイミング用の乱数（状態を保存できる自作版）
c ======================================================================
      subroutine wrestart(iseed, ns1, ns2)
      include 'subcom.inc'
      integer*8 iseed
      integer ns1, ns2

      open(400, file='./restart.tmp', form='unformatted',
     &     status='replace')
      write(400) nx, ny, ipmx
      write(400) istp, time, ip, iseed, ns1, ns2
      write(400) (dp(k),k=1,ipmx), (xp(k),k=1,ipmx),
     &           (yp(k),k=1,ipmx)
      write(400) (uxp(k),k=1,ipmx), (uyp(k),k=1,ipmx),
     &           (pout(k),k=1,ipmx)
      write(400) (cl(k),k=1,ipmx), (clt(k),k=1,ipmx),
     &           (cltmx(k),k=1,ipmx)
      write(400) ((cg(i,j),i=0,nn),j=0,nn)
      write(400) ((cgn(i,j),i=0,nn),j=0,nn)
      close(400)
      call rename('./restart.tmp', './restart.dat')

      return
      end


      subroutine rrestart(iseed, ns1, ns2)
      include 'subcom.inc'
      integer*8 iseed
      integer ns1, ns2
      integer nx0, ny0, ipmx0

      open(400, file='./restart.dat', form='unformatted',
     &     status='old')
      read(400) nx0, ny0, ipmx0
      if( nx0 .ne. nx .or. ny0 .ne. ny .or. ipmx0 .ne. ipmx )then
        write(*,*) 'ERROR: restart.dat is from a different setup'
        write(*,*) '(grid or particle number differs).'
        write(*,*) 'Remove it and run again:  rm restart.dat'
        stop
      end if
      read(400) istp, time, ip, iseed, ns1, ns2
      read(400) (dp(k),k=1,ipmx), (xp(k),k=1,ipmx),
     &          (yp(k),k=1,ipmx)
      read(400) (uxp(k),k=1,ipmx), (uyp(k),k=1,ipmx),
     &          (pout(k),k=1,ipmx)
      read(400) (cl(k),k=1,ipmx), (clt(k),k=1,ipmx),
     &          (cltmx(k),k=1,ipmx)
      read(400) ((cg(i,j),i=0,nn),j=0,nn)
      read(400) ((cgn(i,j),i=0,nn),j=0,nn)
      close(400)

      return
      end


c     --- Park-Miller の線形合同法（0以上1未満の乱数）---
      function prand(iseed)
      integer*8 iseed
      iseed = mod( 16807_8*iseed, 2147483647_8 )
      prand = real(iseed) / 2147483647.
      return
      end

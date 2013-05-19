      program RESAMPLE
c
c  A driver program to test the new resampler
c
      implicit real*8 (a-h,o-z)
      parameter (nmax = 10000, lwork = 100*nmax+100)
      dimension x0(nmax), y0(nmax), rnlx(nmax), rnly(nmax), 
     *          rkappa(nmax), dsdth(nmax)
      dimension xy(2,nmax), der1(2,nmax), der2(2,nmax), 
     *          der3(2,nmax), fis(nmax), work(lwork)
      complex*16 zk, z1, z2
c
         ai = 2.d0
         bi = 1.d0
         zk = dcmplx(0.d0,0.d0)
         thetax = 0.d0  
         nd = 1024   
         call PRINI (6,13)
         open (unit = iplot, file = 'res.m')
         call GEOMETRY (nd,1,ai,bi,thetax,zk,h,x0,y0,iplot,rnlx,
     *                  rnly,rkappa,dsdth)
         NMIN = 60
         NDERS = 2
         n = 5096
         n0 = nd+1
         call RSRESA (ier,x0,y0,n0,nmin,nders,xy,der1,der2,
     *                der3,n,fis,h,acc,err,work,ltot)
         do i = 1,n
            x0(i) = xy(1,i)
            y0(i) = xy(2,i)
            xdot = der1(1,i)
            ydot = der1(2,i)
            xddot = der2(1,i)
            yddot = der2(2,i)
            rkappa(i) = (-xdot*yddot+ydot*xddot)
         end do
         call RSPLOT(x0,rkappa,n,1,iplot)
         if (ier.ne.0) then
            call PRINF ('  IER = *',ier,1)
ccc            stop
         end if
         call PRIN2 ('  X1 = *',xy(1,1),1)
         call PRIN2 ('  XN = *',xy(1,n),1)
         call PRIN2 ('  Y1 = *',xy(2,1),1)
         call PRIN2 ('  YN = *',xy(2,n),1)
         call PRIN2 ('  H = *',h,1)
         call PRIN2 ('  ACC = *',acc,1)
         call PRIN2 ('  ERR = *',err,1)
c
      stop
      end
c
c
c---------------
      subroutine GEOMETRY (nd,k,ai,bi,thetax,zk,h,x,y,iplot,rnlx,
     *                     rnly,rkappa,dsdth)
c---------------
c  Calculates the geometry based on ai and bi (the axes of the ellipse)
c  and zk (the centre of the ellipse)
c
      implicit real*8 (a-h,o-z)
      dimension x(nd+1), y(nd+1), rnlx(nd), rnly(nd), rkappa(nd),
     *          dsdth(nd)
      complex*16 zk, zf(2000), work(10000)
c
         eye = dcmplx(0.d0,1.d0)
         pi = 4.d0*datan(1.d0)
c
c     NMIN, NDERS, IFERR, N1 are parameters of the resampling
c           subroutine RERES3.
c
c
         istart = 0
         do kbod = 1,k
            cosx = dcos(thetax)
            sinx = dsin(thetax)
c
            h = 2.d0*pi/nd
            DO I=1,nd
               th = (i-1.d0)*h
               cx = dreal(zk)
               cy = dimag(zk)
c------------------
c  Star-shaped objects
c     rotated coordinates
c
                  eps = bi
                  NCYC = 4
                  snn = dsin(NCYC*th)
                  csn = dcos(NCYC*th)
                  radius = ai + eps*csn
                  rdot = -NCYC*eps*snn
                  rddot = -NCYC**2*eps*csn
c
                  xp = radius*dcos(th)
                  yp = radius*dsin(th)
                  cs = dcos(th)
                  sn = dsin(th)
                  xpdot = rdot*cs - radius*sn
                  ypdot = rdot*sn + radius*cs
                  xpddot = rddot*cs - 2.d0*rdot*sn - radius*cs
                  ypddot = rddot*sn + 2.d0*rdot*cs - radius*sn
c
                  xdot = xpdot*cosx - ypdot*sinx
                  ydot = xpdot*sinx + ypdot*cosx
                  xddot = xpddot*cosx - ypddot*sinx
                  yddot = xpddot*sinx + ypddot*cosx
c------------------
c 
c  ELLIPSES               
c
c  rotated coordinates
c
               xp = ai*dcos(th)
               yp = bi*dsin(th)
c
c  fixed coordinates
c
               xdot = -ai*dsin(th)*cosx - bi*dcos(th)*sinx
               ydot = -ai*dsin(th)*sinx + bi*dcos(th)*cosx
               xddot = -ai*dcos(th)*cosx + bi*dsin(th)*sinx
               yddot = -ai*dcos(th)*sinx - bi*dsin(th)*cosx
c
c          Make sure normal points towards centre (i.e. outside of
c          the domain
c
c                  
               x(i) = cx + xp*cosx - yp*sinx
               y(i) = cy + xp*sinx + yp*cosx
                rnorm = dsqrt(xdot**2+ydot**2)
              rnlx(i) = -xdot/rnorm
               rnly(i) = ydot/rnorm
               dsdth(i) = rnorm
               rkappa(i) = (-xdot*yddot+ydot*xddot)/(rnorm**3)
               zf(i) = rkappa(i)
c
            ENDDO
            x(nd+1) = x(1)
            y(nd+1) = y(1)
            istart = istart+nd
         end do
ccc         call RSPLOT(x,rkappa,nd,1,iplot)
         call DCFFTI (nd,work)
         call DCFFTF (nd,zf,work)
         do i = 1,nd
            zf(i) = zf(i)/nd
         end do
ccc         call PRIN2 ('  FFT OF KAPPA = *',zf,2*nd)
c
      return
      end
            
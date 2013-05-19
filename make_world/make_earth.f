      program MAKE_EARTH
c
c  Reads in coastline data and uses resampler to generate a smooth curve
      implicit double precision (a-h,o-z)
      parameter (nmax = 15000, lwork = 100*nmax+100)
      dimension xp(nmax), yp(nmax), x0(nmax), y0(nmax), rnlx(nmax),  
     *          rnly(nmax), rkappa(nmax), dsdth(nmax)
      dimension xy(2,nmax), der1(2,nmax), der2(2,nmax), 
     *          der3(2,nmax), fis(nmax), work(lwork), s(nmax)
      complex*16 zk, z1, z2
      dimension x(nmax), y(nmax) 
      complex*16 dzdt(nmax), z(nmax)
c
         call PRINI (6,13)
         call READ_IN_DATA (np, xp, yp)
            call prin2 (' xp = *', xp, np)
            call prin2 (' yp = *', yp, np)
            call PRINF (' np = *', np, 1)
c
         call POLYGON (np, n0, xp, yp, x0, y0)
            call PRINF (' n0 = *', n0, 1)
ccc            call PRIN2 (' x0 = *', x0, n0)
ccc            call PRIN2 (' y0 = *', y0, n0)
            open (unit = 25, file = 'polygon.m')
            call RSPLOT (x0, y0, n0, 1, 25)
            close (25)
ccc         call SQUARE (n0, x0, y0)
ccc            open (unit = 25, file = 'polygon.m')
ccc            call RSPLOT (x0, y0, n0, 1, 25)
ccc            close (25)
c
         NMIN = 128
         NDERS = 2
         n = 4*n0
         call PRINF (' n = *', n, 1)
         call RSRESA (ier,x0,y0,n0,nmin,nders,xy,der1,der2,
     *                der3,n,fis,h,acc,err,work,ltot)
         call PRINF (' ier = *', ier, 1)
         call PRIN2 (' acc = *', acc, 1)
         call PRIN2 (' err = *', err, 1)
         call prin2 (' h = *', h, 1)
         do i = 1,n
            s(i) = (i-1.d0)*h
            x0(i) = xy(1,i)
            y0(i) = xy(2,i)
            xdot = der1(1,i)
            ydot = der1(2,i)
            xddot = der2(1,i)
            yddot = der2(2,i)
            rkappa(i) = (-xdot*yddot+ydot*xddot)
         end do
         call PRIN2 (' ltot = *', ltot, 1)
         call PRIN2 (' x0(1) = *', x0(1), 1)
         call PRIN2 (' y0(1) = *', y0(1), 1)
         call PRIN2 (' x0(end) = *', x0(n), 1)
         call PRIN2 (' y0(end) = *', y0(n), 1)
c
c dump out data
         open (unit=31, file = 'coast.m')
         open (unit=32, file = 'kappa.m')
         call RSPLOT (s, rkappa, n, 1,32)
         call RSPLOT (x0, y0, n, 1, 31)
         close (31)
         close (32)
c
      stop
      end
c
c
c---------------
      subroutine READ_IN_DATA (np, x0, y0)
c---------------
c  read in data for x0, y0 - this contains the vertices of a polygon
      implicit real*8 (a-h,o-z)
      dimension x0(*), y0(*)
c
         open (unit=21, file = 'xy_poly_antarctica.dat', status = 'old')
c
         read (21,*) np
         do i = 1, np
            read(21,*) x, y
            x0(i) = x
            y0(i) = y
         end do
         x0(np+1) = x0(1)
         y0(np+1) = y0(1)
         np = np+1
         close (21)
c
      return
      end
c
c
c---------------
      subroutine POLYGON (np, n0, xp, yp, x0, y0)
c---------------
c  Constructs a polygon out of vertices stored in xp, yp
c  fills in each line segment with ns points
      implicit real*8 (a-h,o-z)
      dimension xp(np), yp(np), x0(*), y0(*)
c
         ns = 20
         do i = 1, np-1
            ind = (ns+1)*(i-1)+1
            x0(ind) = xp(i)
            y0(ind) = yp(i)
            dx = xp(i+1)-xp(i)
            dy = yp(i+1)-yp(i)
            do j = 1, ns
               x0(ind+j) = xp(i) + dx*j/(ns+1)
               y0(ind+j) = yp(i) + dy*j/(ns+1)
            end do
         end do
         n0 = (np-1)*(ns+1)+1
         x0(n0) = x0(1)
         y0(n0) = y0(1)
c
      return
      end
c
c
c---------------
      subroutine SQUARE (nd, x, y)
c---------------
c  construct a resampled square
c
      implicit double precision (a-h,o-z)
      dimension x(*), y(*)
c
c Construct right side
         N_side = 256
         ds = 2.d0/N_side
         do i = 1, N_side
            x(i) = 1.d0
            y(i) = -1.d0 + ds*(i-1.d0)
         end do
c
c Construct top
         N_start = N_side 
         call prinf ('N_start = *', N_start, 1)
         do i = 1, N_side
            x(N_start+i) = 1.d0 - ds*(i-1.d0)
            y(N_start+i) = 1.d0
         end do
c
c Construct left
         N_start = N_start + N_side
         call prinf ('N_start = *', N_start, 1)
         do i = 1, N_side
            x(N_start+i) = -1.d0
            y(N_start+i) = 1.d0 - ds*(i-1.d0)
         end do
c
c Construct bottom
         N_start = N_start + N_side
         call prinf ('N_start = *', N_start, 1)
         do i = 1, N_side
            x(N_start+i) = -1.d0 + ds*(i-1.d0)
            y(N_start+i) = -1.d0
         end do
c
c close
         nd = N_start + N_side
         x(nd+1) = x(1)
         y(nd+1) = y(1)
         nd = nd+1
ccc         call prinf ('nd = *', nd, 1)
ccc         call prin2 (' x = *', x, nd)
ccc         call prin2 (' y = *', y, nd)
c
      return
      end

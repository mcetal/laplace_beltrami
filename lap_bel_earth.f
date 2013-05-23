      PROGRAM LAP_BEL_EARTH
c     ------------------------------------------------------------------
c
c
c  Solves integral equation for Laplace-Beltrami on the Sphere by 
c  using a stereographic projection and the FMM in the complex plane.
c
c  This version uses resampler to define curves parameterized wrt 
c  arc length
c
c
c     ------------------------------------------------------------------
      implicit real*8 (a-h, o-z)
      parameter (kmax = 5, npmax = 75000, nmax = kmax*npmax)
c
c Geometry of holes
      dimension xs(nmax), ys(nmax), zs(nmax)
      dimension diag(nmax)
      complex*16 dzeta(nmax), zeta(nmax)
c
c Hole dimensions
      dimension nd(kmax), n0(kmax), np(kmax)
      dimension ds(kmax), arcl(kmax)
c
c For resampler
      parameter (ireswork = 9*npmax+100)
      dimension xp(nmax), yp(nmax), x0(nmax), y0(nmax)
      dimension xy(2*npmax), der1(2*npmax), der2(2*npmax), 
     *          der3(2*npmax), fis(npmax), reswork(ireswork)
c
c  Grid variables
      parameter (nth_max = 1000, nphi_max = 1000, 
     1          ng_max = nth_max*nphi_max)
      dimension igrid(ng_max), 
     1          u_gr(ng_max), x_gr(ng_max), y_gr(ng_max), z_gr(ng_max),
     2          xzeta_gr(ng_max), yzeta_gr(ng_max), uex_gr(ng_max)
      complex*16 zeta_gr(ng_max)
c
c target points are used to check accuracy
      dimension xz_tar(ng_max), yz_tar(ng_max), u_tar(ng_max)
      complex*16 zeta_tar(ng_max)    
c
c System
      dimension aKu(nmax), density(nmax), A_k(kmax)
c
c Exact potential
      dimension phi_ex(nmax), phi_1e(nmax), phi_2e(nmax), phi_3e(nmax),
     1          phi_1(nmax), phi_2(nmax), phi_3(nmax)
      complex*16 zfield_ex(nmax)
c
c Point Vorticies
      parameter (nvortmax = 100)
      dimension vort_k(nvortmax), x1_vort(nvortmax), x2_vort(nvortmax), 
     1          x3_vort(nvortmax)
      complex*16 zk_vort(nvortmax), zvel(nvortmax)
c
c  Matrix equation variables for GMRES
c  MAXL is the maximum nubmer of GMRES iterations performed
c       before restarting.
c  LRWORK is the dimension of a real workspace needed by DGMRES.
c  LIWORK is the dimension of an integer workspace needed by DGMRES.
c  GMWORK and IWORK are work arrays used by DGMRES
c
      parameter (maxl = 50,liwork=30,  
     1           lrwork=10+(nmax+kmax)*(maxl+6)+maxl*(maxl+3))
      dimension gmwork(lrwork), igwork(liwork)
      dimension rhs(nmax+kmax), soln(nmax+kmax)
c
c  Preconditioner arrays
      dimension IPVTBF(kmax)
      DIMENSION SCHUR(kmax*kmax),WB(kmax)

c
c Fast Multipole Arrays
      parameter (nsp = 20*nmax + 20*ng_max)
      dimension x_zeta(nmax+ng_max), y_zeta(nmax+ng_max)
      complex*16 qa(nmax+ng_max), cfield(nmax+ng_max)
      dimension poten(nmax+ng_max), wksp(nsp)
c
c Other arrays
      dimension alpha(nmax), w(nmax), u(nmax)
      complex*16 zeta_k(kmax), zf(nmax), wsave(20*nmax)
      REAL*4 TIMEP(2), ETIME
c
c common blocks
      common /geometry/ x_zeta, y_zeta, zeta, dzeta, zeta_k
      common /inteqn/ diag, ds, arcl
      common /sys_size/ k, nd, nbk
      common /fasblk2/ schur,wb,ipvtbf
c
c Open output files
ccc         open (unit = 11, file = 'movie/blob.m')
ccc         open (unit = 21, file = 'density.m')
ccc         open (unit = 22, file = 'vort_loc.m')
ccc         open (unit = 31, file = 'igrid.dat')
ccc         open (unit = 32, file = 'xgrid.dat')
ccc         open (unit = 33, file = 'ygrid.dat')
ccc         open (unit = 34, file = 'zgrid.dat')
ccc         open (unit = 35, file = 'xzeta_grid.dat')
ccc         open (unit = 36, file = 'yzeta_grid.dat')
ccc         open (unit = 37, file = 'movie/u_movie.m')
ccc         open (unit = 38, file = 'uex_grid.dat')
ccc         open (unit = 41, file = 'movie/vort_path.m')
ccc         open (unit = 43, file = 'ugrid.dat')
ccc         open (unit = 51, file = 'movie/isl_grid.dat')
ccc         open (unit = 52, file = 'movie/targets.m')
c
c Read hole geometry data
         k = 4
         call PRINI (6,13)
         call READ_IN_DATA (k, np, xp, yp)
         call POLYGON (k, np, n0, xp, yp, x0, y0)
         call READ_IN_VORTICES (nvort, vort_k, zk_vort, x1_vort, 
     1                          x2_vort, x3_vort, gamma_tot)
ccc         call ELLIPSE (k, n0, x0, y0)
c
c Get stereo graphic projection - use resampler in stereographic plane
         nvort = 0
         call STEREO (k, n0, npmax, x0, y0, xy, der1, der2, der3, fis, 
     1                reswork, nd, nbk, zeta, dzeta, x_zeta, y_zeta, 
     2                diag, zeta_k, ds, arcl, xs, ys, zs)
        stop 
ccc         call GEO_TEST (k, n0, npmax, x0, y0, xy, der1, der2, der3, fis, 
ccc     1                reswork, nd, nbk, zeta, dzeta, x_zeta, y_zeta, 
ccc     2                diag, zeta_k, ds, arcl, xs, ys, zs) 
         call TARGET_POINTS (ntar, xz_tar, yz_tar, zeta_tar)
         nx = 500
         ny = 500
         call STEREO_GRID (k, nd, nbk, nx, ny, ds, x_zeta, y_zeta, 
     1                     zeta, dzeta, ntarg, xzeta_gr, yzeta_gr,
     *                     zeta_gr, igrid, qa, cfield, poten, nsp, wksp) 
         call SURFACE_GRID (nx, ny, xzeta_gr, yzeta_gr, x_gr, y_gr,
     *                      z_gr) 
ccc         call DCFFTI (nd(1), wsave)
ccc         call RSCPLOT (zk_vort, nvort, 1, 41)
            call DUMP (nx, ny, u_gr, igrid, 0, 31)
            call DUMP (nx, ny, x_gr, igrid, 1, 32)
            call DUMP (nx, ny, y_gr, igrid, 1, 33)
            call DUMP (nx, ny, z_gr, igrid, 1, 34)
            call DUMP (nx, ny, xzeta_gr, igrid, 1, 35)
            call DUMP (nx, ny, yzeta_gr, igrid, 1, 36)
ccc         stop
ccc         call PRIn2 (' diag_stereo = *', diag_stereo, nbk)
c
c Time loop for vortex path
         tbeg = etime(timep)
         dt = 0.0001d0
         ntime =  1
         do it = 1, ntime
            time = it*dt  
            call PRIN2 (' time = *', time, 1)       
c
c Construct the RHS and solve
         call GETRHS (k, nd, nbk, zeta_k, zeta, rhs,
     1                nvort, vort_k, zk_vort, gamma_tot)
ccc         call PRIN2 (' rhs = *', rhs, nbk)
c
c Construct system matrices to be dumped into matlab
         write (6,*) 'here am i'
         call SOLVE (nd, k, kmax, nbk, rhs, soln, density, A_k, gmwork, 
     1               lrwork, igwork, liwork, maxl, schur, wb,
     2               ipvtbf, zeta, zeta_k)
cccc
cccc Overresolve the data on the boundary (for plotting purposes only)
ccc         call RESAMPLE (nd,k,nbk,nsamp,h,z,dz,xat,yat,u,xsamp,
ccc     *                  ysamp,zwork,ncwork)
c
c Construct solution on surface grid
ccc         call SOL_GRID (nd, k, nbk, nth, nphi, density, A_k, zeta_k,   
ccc     1                  zeta, dzeta, igrid, zeta_gr, u_gr)
         nplot = mod(it,100)
         call PRINF (' nplot = *', nplot, 1)
ccc         if (mod(it,100).eq.0) then
         call SOL_GRID_FMM (nd, k, nbk, nx, ny, ds, density, A_k,    
     1                      zeta_k, zeta, dzeta, igrid, zeta_gr,  
     2                      u_gr, x_zeta, y_zeta, qa, cfield,   
     3                      poten, nsp, wksp, nvort, vort_k, zk_vort)
         call DUMP (nx, ny, u_gr, igrid, 1, 43)
         call SOL_TAR_FMM (nd, k, nbk, ntar, ds, density, A_k, zeta_k,   
     1                     x_zeta, y_zeta, zeta, dzeta, zeta_tar, u_tar,
     2                     xz_tar, yz_tar, qa, cfield, poten, nsp, 
     3                     wksp, nvort, vort_k, zk_vort)
c
c for a vortex in presence of cap with radius r0, check solution
ccc         call SOL_VORT_CHECK (nd, k, nbk, nth, nphi, nvort, zeta_gr,  
ccc     1                        igrid, u_gr, uex_gr, zk_vort, r0)
ccc            call DUMP_MOVIE_ALL (nth, nphi, time, u_gr, it, 37)
ccc            call DUMP_MOVIE_VORT (nth, nphi, time, zk_vort(1), u_gr, 
ccc     1                            it, 37)
ccc          end if
ccc            call DUMP (nth, nphi, uex_gr, igrid, 1, 38)
c
c Calculate velocity at a point
ccc         call CALC_VEL (k, nd, nbk, nvort, density, A_k, zeta, dzeta, 
ccc     1                  zeta_k, vort_k, zk_vort, zvel, zf, wsave)
ccc         do ivort = 1, nvort
ccc            zk_vort(ivort) = zk_vort(ivort) + dt*zvel(ivort)
ccc         end do
         call PRIn2 (' zk_vort = *', zk_vort, 2*nvort)
         call RSCPLOT (zk_vort, nvort, 1, 41)
         end do
         tend = etime(timep)
         call PRIN2 (' TOTAL CPU TIME = *', tend-tbeg, 1)
c
c calculate 
c
c Check error in solution
ccc         call CHECK_ERROR_GRID (k, nx, ny, zeta_k, zeta_gr, igrid, u_gr)
         call CHECK_ERROR_TAR (nd, k, nbk, ntar, zeta_k, zeta_tar,  
     1                         u_tar)
c
         close (31)
         close (32)
         close (33)
         close (34)
         close (35)
         close (36)
         close (37)
c
      stop
      end
c
c
c---------------
      subroutine READ_IN_DATA (k, np, x0, y0)
c---------------
c  read in data for x0, y0 - this contains the vertices of a polygon
      implicit real*8 (a-h,o-z)
      dimension x0(*), y0(*), np(k)
c
c The world is flipped upside down so antarctica contains the north pole
         open (unit=21, file = 'make_world/xy_poly_antarctica.dat', 
     1         status = 'old')
ccc         open (unit=21, file = 'myxy.dat', status = 'old')
c
         read (21,*) npoints
         np(1) = npoints
         do i = 1, np(1)
            read(21,*) x, y
            x0(i) = x
            y0(i) = y
         end do
         x0(np(1)+1) = x0(1)
         y0(np(1)+1) = y0(1)
         np(1) = np(1)+1
         close (21)
c
         open (unit = 25, file = 'original_poly.m')
         call RSPLOT (x0, y0, np(1), 1, 25)
c
c North and South America
c
         open (unit=21, file = 'make_world/xy_poly_northamerica.dat', 
     1         status = 'old')
ccc         open (unit=21, file = 'myxy.dat', status = 'old')
c
         read (21,*) npoints
         np(2) = npoints
         do i = 1, np(2)
            read(21,*) x, y
            x0(i+np(1)) = x
            y0(i+np(1)) = y
         end do
         x0(np(2)+1+np(1)) = x0(1+np(1))
         y0(np(2)+1+np(1)) = y0(1+np(1))
         np(2) = np(2)+1
         close (21)
c
         call RSPLOT (x0(np(1)+1), y0(np(1)+1), np(2), 1, 25)
c
c Australia
c
         open (unit=21, file = 'make_world/xy_poly_australia.dat', 
     1         status = 'old')
ccc         open (unit=21, file = 'myxy.dat', status = 'old')
c
         read (21,*) npoints
         np(3) = npoints
         istart = np(1) + np(2)
         do i = 1, np(3)
            read(21,*) x, y
            x0(istart+i) = x
            y0(istart+i) = y
         end do
         x0(istart+np(3)+1) = x0(istart+1)
         y0(istart+np(3)+1) = y0(istart+1)
         np(3) = np(3)+1
         close (21)
c
         call RSPLOT (x0(istart+1), y0(istart+1), np(3), 1, 25)
c
c Europe and Africa
c
         open (unit=21, file = 'make_world/xy_poly_europe_africa.dat', 
     1         status = 'old')
ccc         open (unit=21, file = 'myxy.dat', status = 'old')
c
         read (21,*) npoints
         np(4) = npoints
         istart = np(1) + np(2) + np(3)
         do i = 1, np(4)
            read(21,*) x, y
            x0(istart+i) = x
            y0(istart+i) = y-0.05d0
         end do
         x0(istart+np(4)+1) = x0(istart+1)
         y0(istart+np(4)+1) = y0(istart+1)
         np(4) = np(4)+1
         close (21)
c
         call RSPLOT (x0(istart+1), y0(istart+1), np(4), 1, 25)
         close (25)
      return
      end
c
c
c---------------
      subroutine READ_IN_VORTICES (nvort, vort_k, zk_vort, x1_vort, 
     1                             x2_vort, x3_vort, gamma_tot)
c---------------
c  Read in vortices and assign random strengths to them
      implicit real*8 (a-h,o-z)
      dimension vort_k(*), x1_vort(*), x2_vort(*), x3_vort(*)
      complex*16 zk_vort(*)
c
         pi = 4.d0*datan(1.d0)
c 
         open (unit = 21, file = 'vort_loc.dat')
         read (21,*) nvort
         do ivort = 1, nvort
            read (21,*) xz, yz
            zk_vort(ivort) = dcmplx(xz,yz)
            call STEREO_TO_PHYS (zk_vort(ivort), x1_vort(ivort),
     1                           x2_vort(ivort), x3_vort(ivort))
         end do
c
c assign random vortex strengths on (-2pi,2pi)
         call zufalli (0)
         call zufall (nvort,vort_k)
         gamma_tot = 0.d0
         do ivort = 1, nvort
            vort_k(ivort) = -2.d0*pi + 4*pi*vort_k(ivort)
            gamma_tot = gamma_tot + vort_k(ivort)
         end do
         close (21)
c
c save to *.m file for plotting
         open (unit=22, file = 'vort_loc.m')
         call RSC_STAR_PLOT(zk_vort,nvort,22)
         close (22)
c
         call prin2 (' gamma_tot = *', gamma_tot, 1)
         call prin2 (' vort_k = *', vort_k, nvort)
         call prin2 (' zk_vort = *', zk_vort, 2*nvort)
c
      return
      end
c
c
c---------------
      subroutine POLYGON (k, np, n0, xp, yp, x0, y0)
c---------------
c  Constructs a polygon out of vertices stored in xp, yp
c  fills in each line segment with ns points
      implicit real*8 (a-h,o-z)
      dimension xp(*), yp(*), x0(*), y0(*), np(k), n0(k)
c
c  ns is the number of points per line segment
         ns = 30
         open (unit = 15, file = 'polygon.m')
         call PRINF (' np = *', np, k)
c
c  loop over all the bodies
         istartp = 0
         istart0 = 0
         do kbod = 1, k
            do i = 1, np(kbod)-1
               ind = (ns+1)*(i-1)+1
               x0(istart0+ind) = xp(istartp+i)
               y0(istart0+ind) = yp(istartp+i)
               dx = xp(istartp+i+1)-xp(istartp+i)
               dy = yp(istartp+i+1)-yp(istartp+i)
               do j = 1, ns
                  x0(istart0+ind+j) = xp(istartp+i) + dx*j/(ns+1)
                  y0(istart0+ind+j) = yp(istartp+i) + dy*j/(ns+1)
               end do
            end do
            istartp = istartp + np(kbod)
            n0(kbod) = (np(kbod)-1)*(ns+1)+1
            x0(istart0+n0(kbod)) = x0(istart0+1)
            y0(istart0+n0(kbod)) = y0(istart0+1)
            call RSPLOT (x0(istart0+1), y0(istart0+1), n0(kbod), 1, 25)
            istart0 = istart0 + n0(kbod)
         end do
         call PRINF (' n0 = *', n0, k)
c
         close (25)
c
      return
      end
c
c
c---------------
      subroutine ELLIPSE (k, n0, x0, y0)
c---------------
c  Constructs a polygon out of vertices stored in xp, yp
c  fills in each line segment with ns points
      implicit real*8 (a-h,o-z)
      dimension x0(*), y0(*), n0(k)
c 
         pi = 4.d0*datan(1.d0)
         open (unit = 25, file = 'polygon.m')
c
         n0(1) = 128
         n0(2) = 128
         dth = 2.d0*pi/n0(1)
         do i = 1, n0(1)
            theta = dth*(i-1.d0)
            x0(i) = 2.d0*dcos(theta)
            y0(i) = 2.d0*dsin(theta)
         end do
         x0(n0(1)+1) = x0(1)
         y0(n0(1)+1) = y0(1)
         n0(1) = n0(1)+1
         call RSPLOT (x0, y0, n0(1), 1, 25)
c
         dth = 2.d0*pi/n0(2)
         do i = 1, n0(2)
            theta = dth*(i-1.d0)
            x0(n0(1)+i) = 0.25d0*dcos(theta)
            y0(n0(1)+i) = 0.25d0*dsin(theta)
         end do
         x0(n0(1)+n0(2)+1) = x0(n0(1)+1)
         y0(n0(1)+n0(2)+1) = y0(n0(1)+1)
         n0(2) = n0(2)+1
c
         call RSPLOT (x0(n0(1)+1), y0(n0(1)+1), n0(2), 1, 25)
         close (25)
c
      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine READ_DATA (k, nd, nbk, nth, nphi, ak, bk, cx, cy, cz, 
     1                      th_k, phi_k, nvort, x1_vort, x2_vort, 
     2                      x3_vort, vort_k, gamma_tot, r0)
c---------------
      implicit real*8 (a-h,o-z)
      dimension ak(*), bk(*), cx(*), cy(*), cz(*), th_k(*), phi_k(*)
      dimension x1_vort(*), x2_vort(*), x3_vort(*), vort_k(*)
c
         open (unit = 12, file = 'input.data')
c
         read (12,*) k, nd, nvort
         nbk = k*nd
         call PRINF (' nbk = *', nbk, 1)
         read(12,*) nth, nphi
         do kbod = 1, k
            read(12,*) ak(kbod), bk(kbod), th_k(kbod), phi_k(kbod)
            call SPH2CART (th_k(kbod),phi_k(kbod), 1.d0, cx(kbod),
     1                     cy(kbod), cz(kbod))
         end do
         call PRINF (' nvort = *', nvort, 1)
         gamma_tot = 0.d0
         do ivort = 1, nvort
            read(12,*) theta, phi, vort_k(ivort)
            call SPH2CART (theta, phi, 1.d0, x1_vort(ivort),
     1                     x2_vort(ivort), x3_vort(ivort))
            gamma_tot = gamma_tot + vort_k(ivort)
         end do
         call PRINI (6,13) 
         call PRIN2 (' ak = *', ak, k)
         call PRIN2 (' bk = *', bk, k)
         call PRIN2 (' cx = *', cx, k)
         call PRIN2 (' cy = *', cy, k)
         call PRIN2 (' cz = *', cz, k)
         call PRINF (' nth = *', nth, 1)
         call PRINF (' nphi = *', nphi, 1)
         call PRIN2 ('    vort_k = *', vort_k, nvort)
         call PRIN2 ('    x1_vort = *', x1_vort, nvort)
         call PRIN2 ('    x2_vort = *', x2_vort, nvort)
         call PRIN2 ('    x3_vort = *', x3_vort, nvort)
         r0 = ak(1)/(1.d0-dsqrt(1.d0-(ak(1))**2))
         call PRIN2 (' r0 = *', r0, 1)
c
         close(12)
c
      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine SPH2CART (theta, phi, r, x, y, z)
c---------------
      implicit real*8 (a-h,o-z)
c
         x = r*dcos(phi)*dcos(theta)
         y = r*dcos(phi)*dsin(theta)
         z = r*dsin(phi)
c
      return
      end      
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine CROSS (u, v, w)
c---------------
      implicit real*8 (a-h,o-z)
      dimension u(3), v(3), w(3)
c
         w(1) = u(2)*v(3) - v(2)*u(3)
         w(2) = v(1)*u(3) - u(1)*v(3)
         w(3) = u(1)*v(2) - v(1)*u(2)
c
      return
      end      
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine TARGET_POINTS (ntar, xz_tar, yz_tar, zeta_tar)
c---------------
c  Define target points in complex plane to check error
      implicit real*8 (a-h,o-z)
      dimension xz_tar(*), yz_tar(*)
      complex*16 zeta_tar(*), eye
c
         eye = dcmplx(0.d0,1.d0)
         pi = 4.d0*datan(1.d0)
c
         ntar = 20
         dth = 2.d0*pi/ntar
         do i = 1, ntar
            theta = dth*(i-1.d0)
            zeta_tar(i) = dcmplx(-14.d0,1.5d0)+ 1.d0*cdexp(eye*theta)
ccc            zeta_tar(i) = 1.d0*cdexp(eye*theta)
            xz_tar(i) = dreal(zeta_tar(i))
            yz_tar(i) = dimag(zeta_tar(i))
         end do
c
         open (unit = 25, file = 'stereo_targets.m')
         call RSCPLOT (zeta_tar, ntar, 1, 25)
         close (25)
                  
c
      return
      end      
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine DOT (u, v, u_dot_v)
c---------------
      implicit real*8 (a-h,o-z)
      dimension u(3), v(3)
c
         u_dot_v = 0.d0
         do i = 1, 3
            u_dot_v = u_dot_v + u(i)*v(i)
         end do
c
      return
      end      
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine R_FUNC (alpha, A, B, th, phi, x, y, z)
c---------------
      implicit real*8 (a-h,o-z)
      dimension zaxis(3), xaxis(3), yaxis(3)
c
         pi = 4.d0*datan(1.d0)
c
         call SPH2CART (th, phi, 1.d0, zaxis(1), zaxis(2), zaxis(3))
         call SPH2CART (th, phi-0.5d0*pi, 1.d0, xaxis(1), xaxis(2), 
     1                  xaxis(3))
         call CROSS (zaxis, xaxis, yaxis)
         xp = A*dcos(alpha)
         yp = B*dsin(alpha)
         zp = dsqrt(1.d0 - xp**2 - yp**2)
         x = xp*xaxis(1) + yp*yaxis(1) + zp*zaxis(1)
         y = xp*xaxis(2) + yp*yaxis(2) + zp*zaxis(2)
         z = xp*xaxis(3) + yp*yaxis(3) + zp*zaxis(3)
c
      return
      end      
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine DR_FUNC (alpha, A, B, th, phi, dx, dy, dz)
c---------------
      implicit real*8 (a-h,o-z)
      dimension zaxis(3), xaxis(3), yaxis(3)
c
         pi = 4.d0*datan(1.d0)
c
         call SPH2CART (th, phi, 1.d0, zaxis(1), zaxis(2), zaxis(3))
         call SPH2CART (th, phi-0.5d0*pi, 1.d0, xaxis(1), xaxis(2), 
     1                  xaxis(3))
         call CROSS (zaxis, xaxis, yaxis)
         xp = A*dcos(alpha)
          dxp = -A*dsin(alpha)
         yp = B*dsin(alpha)
          dyp = B*dcos(alpha)
         zp = dsqrt(1.d0 - xp**2 - yp**2)
          dzp = (-xp*dxp-yp*dyp)/zp
         dx = dxp*xaxis(1) + dyp*yaxis(1) + dzp*zaxis(1)
         dy = dxp*xaxis(2) + dyp*yaxis(2) + dzp*zaxis(2)
         dz = dxp*xaxis(3) + dyp*yaxis(3) + dzp*zaxis(3)
c
      return
      end      
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine D2R_FUNC (alpha, A, B, th, phi, d2x, d2y, d2z)
c---------------
      implicit real*8 (a-h,o-z)
      dimension zaxis(3), xaxis(3), yaxis(3)
c
         pi = 4.d0*datan(1.d0)
c
         call SPH2CART (th, phi, 1.d0, zaxis(1), zaxis(2), zaxis(3))
         call SPH2CART (th, phi-0.5d0*pi, 1.d0, xaxis(1), xaxis(2), 
     1                  xaxis(3))
         call CROSS (zaxis, xaxis, yaxis)
         xp = A*dcos(alpha)
          dxp = -A*dsin(alpha)
          d2xp = -A*dcos(alpha)
         yp = B*dsin(alpha)
          dyp = B*dcos(alpha)
          d2yp = -B*dsin(alpha)
         zp = dsqrt(1.d0 - xp**2 - yp**2)
          dzp = (-xp*dxp-yp*dyp)/zp
          d2zp = (-dxp**2 - xp*d2xp - dyp**2 - yp*d2yp - dzp**2)/zp
         d2x = d2xp*xaxis(1) + d2yp*yaxis(1) + d2zp*zaxis(1)
         d2y = d2xp*xaxis(2) + d2yp*yaxis(2) + d2zp*zaxis(2)
         d2z = d2xp*xaxis(3) + d2yp*yaxis(3) + d2zp*zaxis(3)
c
      return
      end      
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine MAKE_GEO (k, nd, nbk, ak, bk, th_k, phi_k, xs, ys, zs,
     1                     dx, dy, dz, d2x, d2y, d2z, xn, yn, zn, dsda, 
     2                     diag)
c---------------
      implicit real*8 (a-h,o-z)
      dimension ak(k), bk(k), th_k(k), phi_k(k), r(3), t(3), pn(3),
     1          vn(3), diag(nbk), d2x(nbk), d2y(nbk), d2z(nbk)
      dimension xs(nbk), ys(nbk), zs(nbk), dx(nbk), dy(nbk), dz(nbk),
     1          xn(nbk), yn(nbk), zn(nbk), dsda(nbk)
c
         pi = 4.d0*datan(1.d0)
c
         dalph = 2.d0*pi/nd
         istart = 0
         do kbod = 1, k
            do i = 1, nd
               alpha = dalph*(i-1.d0)
               call R_FUNC (alpha, ak(kbod), bk(kbod), th_k(kbod), 
     1                      phi_k(kbod), xs(istart+i), ys(istart+i), 
     2                      zs(istart+i))
               call DR_FUNC (alpha, ak(kbod), bk(kbod), th_k(kbod), 
     1                      phi_k(kbod), dx(istart+i), dy(istart+i), 
     2                      dz(istart+i))
c
c Construct normal to surface of sphere
               r(1) = xs(istart+i)
               r(2) = ys(istart+i)
               r(3) = zs(istart+i)
c
c Construct tangent to curve lying in plane of sphere
               t(1) = dx(istart+i)
               t(2) = dy(istart+i)
               t(3) = dz(istart+i)
               dsda(istart+i) = dsqrt((t(1))**2 + (t(2))**2 + (t(3))**2)
               do j = 1, 3
                  t(j) = t(j)/dsda(istart+i)
               end do
c
c Construct normal to curve lying in plane of sphere
               call CROSS (t,r,vn)
               xn(istart+i) = vn(1)
               yn(istart+i) = vn(2)
               zn(istart+i) = vn(3)
c
c Construct the diagonal entry
               call D2R_FUNC (alpha, ak(kbod), bk(kbod), th_k(kbod), 
     1                        phi_k(kbod), pn(1), pn(2), pn(3))
               d2x(istart+i) = pn(1)
               d2y(istart+i) = pn(2)
               d2z(istart+i) = pn(3)
               pn(1) = pn(1)/(dsda(istart+i))**2
               pn(2) = pn(2)/(dsda(istart+i))**2
               pn(3) = pn(3)/(dsda(istart+i))**2
               call CROSS(pn,r,vn)
               call DOT (t,vn,t_dot_vn)
               diag(istart+i) = -t_dot_vn*dsda(istart+i)/(4.d0*pi)
            end do
            istart = istart + nd
         end do
         is = 1
         do kbod = 1, k
            call RS_3D_PLOT (xs(is),ys(is),zs(is), nd, 1, 42)
            is = is + nd
         end do
c
      return
      end      
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine IN_OR_OUT (x, y, z, eps, ak, bk, th_k, phi_k, itest)
c---------------
      implicit real*8 (a-h,o-z)
      dimension p(3), x_ax(3), y_ax(3), z_ax(3)
c
         pi = 4.d0*datan(1.d0)
c
         p(1) = x
         p(2) = y
         p(3) = z
         call SPH2CART (th_k, phi_k, 1.d0, z_ax(1), z_ax(2), z_ax(3))
         call SPH2CART (th_k, phi_k-0.5d0*pi, 1.d0, x_ax(1), x_ax(2), 
     1                  x_ax(3))
         call CROSS (z_ax, x_ax, y_ax)
         call DOT (p, x_ax, x1)
         call DOT (p, y_ax, y1)
         call DOT (p, z_ax, z1)
         rad = dsqrt((x1/ak)**2 + (y1/bk)**2)
         if ((rad<1.d0+eps).and.(z1>0.d0)) then 
            itest = 1
           else
            itest = 0
         end if
c
      return
      end      
c
c*************************************************
c
      subroutine STEREO_GRID (k, nd, nbk, nx, ny, ds, x_zeta, y_zeta, 
     1                        zeta, dzeta, ntarg, xzeta_gr, yzeta_gr,
     2                        zeta_gr, igrid, qa, cfield, poten, nsp, 
     3                        wksp) 
c---------------
c  Construct grid points in stereographic plane
c
      implicit real*8 (a-h,o-z)
      parameter (npmax = 100000, ntmax = 100000)
      dimension xzeta_gr(nx,ny), yzeta_gr(nx,ny), igrid(nx,ny)
      dimension nd(k), ds(k), x_zeta(*), y_zeta(*)
      complex*16 zeta(nbk), dzeta(nbk)
      complex*16 eye
      complex*16 qa(*), cfield(*), ztar, zeta_gr(nx,ny)
      dimension poten(*), wksp(nsp)
      integer*4 iout(2), inform(10), ierr(10)
      REAL*4 TIMEP(2), ETIME
c
         pi = 4.0d0*datan(1.d0)
         eye = dcmplx(0.d0,1.d0)
c
C  DIMENSIONS OF EMBEDDING RECTANGLE
         call MINMAX (nbk, x_zeta, xl, xr)
         call MINMAX (nbk, y_zeta, yb, yt)
cccc
cccc  bounding box for the continents
ccc         xl = -3.d0
ccc         xr = 3.d0
ccc         yb = -2.d0
ccc         yt = 3.4d0
c
         call prin2 ('xmin = *', xl, 1)
         call prin2 ('xmax = *', xr, 1)
         call prin2 ('ymin = *', yb, 1)
         call prin2 ('ymax = *', yt, 1)
         hx=(xr-xl)/DFLOAT(nx-1)
         hy=(yt-yb)/DFLOAT(ny-1)
         call prin2 (' hx = *', hx, 1)
         call prin2 (' hy = *', hy, 1)
c
c  Zero igrid values
c
         do i = 1,nx
            do j = 1,ny
               igrid(i,j) = 0
            end do
         end do
c
c  Determine which grid points are in the domain and which are outside
c     igrid(i,j) = 1 if (x,y) is outside all holes, ie INSIDE the
c                    domain
c     igrid(i,j) = 0 if (x,y) is inside a hole, ie OUTSIDE the 
c                    domain
c  First sort through each boundary point and remove all grid points that
c  are within 5*DS of that point.
c  Set 
c     igrid(i,j) = -1 if grid points close to a boundary
c
         istart = 0
         do nbod = 1, k
            do i = 1,nd(nbod)
               eps = ds(nbod)
               xmin = dreal(zeta(istart+i)) - eps
               xmax = dreal(zeta(istart+i)) + eps
               ymin = dimag(zeta(istart+i)) - eps
               ymax = dimag(zeta(istart+i)) + eps
               ileft = IDINT((xmin - xl)/hx + 1.d0)
                  if (ileft.lt.1) ileft = 1
               iright = DNINT((xmax - xl)/hx + 1.d0)
                  if (iright.gt.nx) iright = nx
               jbot = IDINT((ymin - yb)/hy + 1.d0)
                  if (jbot.lt.1) jbot = 1
               jtop = DNINT((ymax - yb)/hy + 1.d0)
                  if (jtop.gt.ny) jtop = ny
               do ix = ileft,iright
                  do iy = jbot,jtop
                     igrid(ix,iy) = -1
                  end do
               end do
            end do
            istart = istart + nd(nbod)
         end do
c
c  Pack unmarked points into xtarg and ytarg
c
         itarg = 1
         do i = 1,nx
            x = xl + (i-1)*hx
            do j = 1,ny
               y = yb + (j-1)*hy
               xzeta_gr(i,j) = x
               yzeta_gr(i,j) = y
               zeta_gr(i,j) = dcmplx(x,y)
               if (igrid(i,j).ge.0) then
                  x_zeta(nbk+itarg) = x
                  y_zeta(nbk+itarg) = y
                  qa(nbk+itarg) = 0.d0
                  itarg = itarg + 1
               end if
            end do
         end do
         ntarg = itarg - 1
         call prinf ('ntarg after test 1 = *', ntarg, 1)
c
c  Now evaluate cauchy integral with density = z using CADAP
c  (only PHIP is of interest here)
c
         istart = 0
         do kbod = 1, k
            do i = 1, nd(kbod)
               qa(istart+i) = ds(kbod)*1.d0
     1                           *dzeta(istart+i)/(2.d0*pi)
            end do
            istart = istart + nd(kbod)
         end do
         tbeg = etime(timep)
         iout(1) = 0
         iout(2) = 13
         iflag7 = 3
         napb = 50
         ninire = 2
         mex = 300
         eps7 = 1.d-14
         tol = 1.d-14
         nnn = nbk+ntarg
ccc         nnn = nbk
ccc         nnn = nbk+100
         call DAPIF2 (iout, iflag7, nnn, napb, ninire, mex, ierr, 
     &                inform, tol, eps7, x_zeta, y_zeta, qa, poten,  
     &                cfield, wksp, nsp, CLOSE)
         call PRINI (6, 13)
ccc         call PRIN2 (' cfield = *', cfield, 2*nnn)
         call PRINF (' Number of Levels = *', inform(3), 1)
         if (ierr(1).ne.0) then
            write (6,*) '  ERROR IN DAPIF2, IERR = ', (ierr(ii),ii=1,6)
            write(6,*) '  INFORM = ', (inform(ii),ii=1,6)
            stop
         end if
c
c  Test to see of points are in the domain, using the value of 
c     imag(cfield)
c     PHIP(z) = -1  if  z is in a hole (outside the domain)
c     PHIP(z) = 1   if z is in the bounded domain (k0 = 0), 
c                      but not in a hole
c     PHIP(z) = 0   if z is outside the domain if k0 = 0, inside the
c                      domain if k0 = 1
c  Let TEST = PHIP + k0, then if 
c     TEST = 1,        the point is in the domain and igrid = 0
c     TEST = 0 or -1,  the point is outside the domain and igrid = -1
c
         itarg = 1
         eps = .1
         do i = 1,nx
            x = xl + (i-1)*hx
            do j = 1,ny
               y = yb + (j-1)*hy
               if (igrid(i,j).ge.0) then
                  test = dimag(cfield(nbk+itarg))
                  igrid(i,j) = DNINT(test) 
ccc                  call prin2 (' test = *', test, 1)
ccc                  call prinf (' igrid = *', igrid(i,j), 1)
ccc                  call prin2 (' x = *', x, 1)
ccc                  call prin2 (' y = *', y, 1)
                  itarg = itarg + 1
               end if
            end do
         end do
c
      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine SURFACE_GRID (nx, ny, xzeta_gr, yzeta_gr, x_gr, y_gr,
     *                         z_gr)
c---------------
      implicit real*8 (a-h,o-z)
      dimension xzeta_gr(nx,ny), yzeta_gr(nx,ny), x_gr(nx,ny),
     1          y_gr(nx,ny), z_gr(nx,ny)
      complex*16 eye, zeta
c
         pi = 4.d0*datan(1.d0)
         eye = dcmplx(0.d0,1.d0)
c
         do i = 1, nx
            do j = 1, ny
               zeta = dcmplx(xzeta_gr(i,j),yzeta_gr(i,j))
               call STEREO_TO_PHYS (zeta, x_gr(i,j), y_gr(i,j), 
     1                              z_gr(i,j))
               y_gr(i,j) = -y_gr(i,j)
               z_gr(i,j) = -z_gr(i,j)
            end do
         end do
c
      return
      end      
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine STEREO (k, n0, npmax, x0, y0, xy, der1, der2, der3,  
     1                fis, reswork, nd, nbk, zeta, dzeta, x_zeta,  
     2                y_zeta, diag, zeta_k, ds, arcl, xs, ys, zs)
c---------------
c  takes curve data in x0(n0), y0(n0) and resamples to construct 
c  zeta, dzeta, x_zeta, y_zeta, diag and zeta_k
c  also calculates system size nd, k, nbk
c  curves on sphere given in (xs,ys,zs)
c  xy, der1, der2, der3, fix, reswork used by resampler
c
      implicit real*8 (a-h,o-z)
      dimension n0(k), nd(k)
      dimension xs(nbk), ys(nbk), zs(nbk), dx(nbk), dy(nbk), dz(nbk),
     1          x_zeta(nbk), y_zeta(nbk), d2x(nbk), d2y(nbk), d2z(nbk),
     2          diag(nbk), ds(k), arcl(k)
      complex*16 zeta(nbk), dzeta(nbk), eye, d2zeta, zextra, zeta_k(k)
c
c resampler arrays
      dimension x0(*), y0(*), xy(2,npmax), der1(2,npmax), 
     1          der2(2,npmax), der3(2,npmax), fis(*), reswork(*)
c
         pi = 4.d0*datan(1.d0)
         eye = dcmplx(0.d0,1.d0)
c
c  parameters for resampler
         NMIN = 512
ccc         NMIN = 2
         NDERS = 2
c
c  open file for zeta_k
         open (unit = 26, file = 'make_world/zeta_k.dat', 
     1         status = 'old')
c
c  open file for plotting resampled curve
         open (unit = 25, file = 'geo_stereo.m')
         open (unit = 27, file = 'geo_3d.m')
c
c  loop over each body
         istart = 0
         istart0 = 0
         nbk = 0
         do kbod = 1, k
            call PRINF ('In resampler, kbod = *', kbod, 1)
            nd(kbod) = 4*n0(kbod)
            nbk = nbk+nd(kbod)
            call RSRESA (ier, x0(istart0+1), y0(istart0+1), n0(kbod),  
     1                   nmin, nders, xy, der1, der2, der3, nd(kbod),   
     2                   fis, ds(kbod), acc, err, reswork, ltot)
            if (ier.ne.0) then
               call PRINF ('ERROR IN RESAMPLER, IER = *', ier, 1)
               stop
            end if
            call PRIN2 ('   acc = *', acc, 1)
            call PRIN2 ('   err = *', err, 1)
            call prin2 ('   ds = *', ds(kbod), 1)         
            do i = 1,nd(kbod)
               x_zeta(istart+i) = xy(1,i)
               y_zeta(istart+i) = xy(2,i)
               zeta(istart+i) = 
     1                   dcmplx(x_zeta(istart+i),y_zeta(istart+i))
               if (kbod.eq.1) then
                  xdot = -der1(1,i)
                  ydot = -der1(2,i)
                 else
                  xdot = der1(1,i)
                  ydot = der1(2,i)
               end if
               dzeta(istart+i) = dcmplx(xdot,ydot)
               xddot = der2(1,i)
               yddot = der2(2,i)
               rkappa = (xdot*yddot-ydot*xddot)
ccc               call PRIN2 (' rkappa = *', rkappa, 1)
               zextra = dzeta(istart+i)*dconjg(zeta(istart+i))
     1                      /(1.d0+cdabs(zeta(istart+i))**2)
               zextra = zextra/(2.d0*pi)
               diag(istart+i) = 0.25d0*rkappa/pi - dimag(zextra)
               call STEREO_TO_PHYS (zeta(istart+i), xs(istart+i), 
     1                              ys(istart+i), zs(istart+i))
               ys(istart+i) = -ys(istart+i)
               zs(istart+i) = -zs(istart+i)
            end do
ccc            call prin2 (' dzeta = *', dzeta(istart+1), 2*nd(kbod))
ccc            call prin2 (' diag = *', diag(istart+1), nd(kbod))
            arcl(kbod) = nd(kbod)*ds(kbod)
            call PRIN2 ('   arc length = *', arcl(kbod), 1)
            read (26, *) xzetak, yzetak
            zeta_k(kbod) = dcmplx(xzetak,yzetak)
            call RSPLOT (x_zeta(istart+1), y_zeta(istart+1), nd(kbod), 
     1                   1, 25)
            call RS_3D_PLOT (xs(istart+1), ys(istart+1), zs(istart+1),  
     1                       nd(kbod), 1, 27)
            istart = istart + nd(kbod)
            istart0 = istart0 + n0(kbod)
         end do
         call PRIN2 (' zeta_k = *', zeta_k, 2*k)
         call PRINF (' nd = *', nd, k)
         call PRINF (' nbk = *', nbk, 1)
         close (25)
         close (26)
         close (27)
c
      return
      end  
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine GEO_TEST (k, n0, npmax, x0, y0, xy, der1, der2, der3,  
     1                fis, reswork, nd, nbk, zeta, dzeta, x_zeta,  
     2                y_zeta, diag, zeta_k, ds, arcl, xs, ys, zs)
c---------------
c  trying to debug
c
      implicit real*8 (a-h,o-z)
      dimension n0(k), nd(k)
      dimension xs(nbk), ys(nbk), zs(nbk), dx(nbk), dy(nbk), dz(nbk),
     1          x_zeta(nbk), y_zeta(nbk), d2x(nbk), d2y(nbk), d2z(nbk),
     2          diag(nbk), ds(k), arcl(k)
      complex*16 zeta(nbk), dzeta(nbk), eye, d2zeta, zextra, zeta_k(k)
c
c resampler arrays
      dimension x0(*), y0(*), xy(2,npmax), der1(2,npmax), 
     1          der2(2,npmax), der3(2,npmax), fis(*), reswork(*)
c
         pi = 4.d0*datan(1.d0)
         eye = dcmplx(0.d0,1.d0)
c
c  parameters for resampler
ccc         NMIN = 128
         NMIN = 2
         NDERS = 2
c
c  open file for zeta_k
         open (unit = 26, file = 'zeta_k.m', status = 'old')
c
c  open file for plotting resampled curve
         open (unit=25, file = 'coast.m')
c
c
         nd(1) = 256
         rad1 = 2.d0
         nd(2) = 256
         rad2 = 0.25d0
         nbk = 0
         do kbod = 1, k
            nbk = nbk+nd(kbod)
         end do
         call prinf (' nbk = *', nbk, 1)
         call prinf (' nd = *', nd, 2)
         dth = 2.d0*pi/nd(1)
         ds(1) = rad1*dth
         arcl(1) = 2.d0*pi*rad1
         do i = 1, nd(1)
            theta = dth*(i-1.d0)
            x_zeta(i) = rad1*dcos(theta)
            y_zeta(i) = rad1*dsin(theta)
            zeta(i) = dcmplx(x_zeta(i),y_zeta(i))
            dzeta(i) = -dcmplx(-dsin(theta),dcos(theta))
            rkappa = -1.d0/rad1
            zextra = dzeta(i)*dconjg(zeta(i))
     1                      /(1.d0+cdabs(zeta(i))**2)
            zextra = zextra/(2.d0*pi)
            diag(i) = 0.25d0*rkappa/pi - dimag(zextra)
            call STEREO_TO_PHYS (zeta(i), xs(i), 
     1                              ys(i), zs(i))
         end do
         call RSPLOT (x_zeta, y_zeta, nd(1), 1, 25)
         read(26,*) xx, yy
         zeta_k(1) = dcmplx(xx,yy)
         call PRINf ('kbod = *', 1, 1)
         call prin2 (' zeta = *', zeta, 2*nd(1))
         call prin2 (' x_zeta = *', x_zeta, nd(1))
         call prin2 (' y_zeta = *', y_zeta, nd(1))
         call prin2 (' diag = *', diag, nd(1))
         call prin2 (' dzeta = *', dzeta, 2*nd(1))
c
         istart = nd(1)
         dth = 2.d0*pi/nd(2)
         ds(2) = rad2*dth
         arcl(2) = 2.d0*pi*rad2
         do i = 1, nd(2)
            theta = dth*(i-1.d0)
            x_zeta(istart+i) = rad2*dcos(theta)
            y_zeta(istart+i) = rad2*dsin(theta)
            zeta(istart+i) = dcmplx(x_zeta(istart+i),y_zeta(istart+i))
            dzeta(istart+i) = dcmplx(-dsin(theta),dcos(theta))
            rkappa = 1.d0/rad2
            zextra = dzeta(istart+i)*dconjg(zeta(istart+i))
     1                      /(1.d0+cdabs(zeta(istart+i))**2)
            zextra = zextra/(2.d0*pi)
            diag(istart+i) = 0.25d0*rkappa/pi - dimag(zextra)
            call STEREO_TO_PHYS (zeta(istart+i), xs(istart+i), 
     1                              ys(istart+i), zs(istart+i))
         end do
         read(26,*) xx, yy
         zeta_k(2) = dcmplx(xx,yy)
         call RSPLOT (x_zeta(istart+1), y_zeta(istart+1), nd(2), 1, 25)
         close (25)
         close (26)
         call PRINf ('kbod = *', 2, 1)
         call prin2 (' zeta = *', zeta(istart+1), 2*nd(2))
         call prin2 (' x_zeta = *', x_zeta(istart+1), nd(2))
         call prin2 (' y_zeta = *', y_zeta(istart+1), nd(2))
         call prin2 (' diag = *', diag(istart+1), nd(2))
         call prin2 (' dzeta = *', dzeta(istart+1), 2*nd(2))
         call prin2 (' ds = *', ds, k)
         call prin2 (' zeta_k = *', zeta_k, 2*k)
c
      return
      end  
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine MINMAX (nd, x, xmin, xmax)
c---------------
c  find min and max of array x(nd)
      implicit real*8 (a-h,o-z)
      dimension x(nd)
c
         xmin = 1.d10
         xmax = -1.d10
         do i = 1, nd
            xmin = min(xmin, x(i))
            xmax = max(xmax, x(i))
         end do
c     
      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine STEREO_TO_PHYS (zeta, x, y, z)
c---------------
c mapping zeta --> (x,y,z) on sphere
c
      implicit real*8 (a-h,o-z)
      complex*16 zeta, eye
c
         eye = dcmplx(0.d0,1.d0)
c
         den = 1.d0 + dreal(zeta*dconjg(zeta))
         x = dreal(zeta+dconjg(zeta))/den
         y = dreal((zeta-dconjg(zeta))/(eye*den))
         z = dreal((-1.d0+zeta*dconjg(zeta))/den)
c
      return
      end     
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine PHYS_TO_STEREO (x, y, z, zeta)
c---------------
c mapping (x,y,z) on sphere --> zeta on complex plane
c
      implicit real*8 (a-h,o-z)
      complex*16 zeta, eye
c
         eye = dcmplx(0.d0,1.d0)
c
         zeta = (x+eye*y)/(1.d0-z)
c
      return
      end     
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine GEO_INIT (nphi, nth, n, theta, phi, x1, x2, x3, z, q)
c---------------
      implicit real*8 (a-h,o-z)
      dimension theta(*), phi(*), x1(*), x2(*), x3(*)
      complex*16 q(*), z(*), eye
c
         pi = 4.d0*datan(1.d0)
         eye = dcmplx(0.d0,1.d0)
c         
         nth = 1
         nphi = 50
         n = nth*nphi
         dth = 0.5d0*pi/nth
         dphi = 2.d0*pi/nphi
c         
         istart = 0
         do i = 1, nth
            th = 0.25d0*pi + (i-1)*dth
            do j = 1, nphi
               theta(istart+j) = th
               ph = (j-1)*dphi
               phi(istart+j) = ph
               x1(istart+j) = dcos(ph)*dsin(th)
               x2(istart+j) = dsin(ph)*dsin(th)
               x3(istart+j) = dcos(th)
               z(istart+j) = (x1(istart+j) + eye*x2(istart+j))/
     1                         (1 - x3(istart+j))
               q(istart+j) = 1.d0
            end do
            istart = istart+nphi
         end do
c
         call PRINI (6,13)
c
      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine GETRHS (k, nd, nbk, zeta_k, zeta, rhs,
     1                   nvort, vort_k, zk_vort, gamma_tot)
c---------------
      implicit real*8 (a-h,o-z)
      dimension rhs(nbk), xs(nbk), ys(nbk), nd(k), 
     1          zs(nbk), vort_k(nvort)
      complex*16 zeta_k(k), eye, zeta(nbk), zk_vort(nvort), zdis
c
         eye = dcmplx(0.d0,1.d0)
c
         istart = 0
         call PRIN2 (' in rhs, zeta_k = *', zeta_k, 2*k)
         do kbod = 1, k
            do j = 1, nd(kbod)
               psi = 0.d0
               do mbod = 1, k
                  psi = psi + 1.d0/(zeta(istart+j)-zeta_k(mbod)) 
     1                   + dconjg(1.d0/(zeta(istart+j)-zeta_k(mbod)))
               end do
ccc               do ivort = 1, nvort
ccc                  zdis = zeta(istart+j) - zk_vort(ivort)
ccc                  az1 = (cdabs(zeta(istart+j)))**2
ccc                  az2 = (cdabs(zk_vort(ivort)))**2
ccc                  arg = dreal(zdis*conjg(zdis)/((1.d0+az1)*(1.d0+az2)))
ccc                  psi = psi + vort_k(ivort)*dlog(arg)
ccc               end do
               rhs(istart+j) = psi 
            end do
            istart = istart+nd(kbod)
         end do
         rhs(nbk+1) = gamma_tot
         do kbod = 2, k
            rhs(nbk+kbod) = 0.d0
         end do
         call PRIN2 (' rhs at end = *', rhs(nbk+1), k)
ccc         rhs(nbk+k) = gamma_tot
ccc         rhs(nbk+k) = -gamma_tot
c
c Dump it out
         open (unit = 24, file = 'rhs.dat')
         do i = 1, nbk+k
               write(24,'(e20.13,$)')(rhs(i))
               write (24,'(a)')  ''
         end do
         close (24) 
c
      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine EXACT_POT (n, nsp, x, y, z, qa, phi)
c---------------
      implicit real*8 (a-h,o-z)
      dimension x(n), y(n), phi(n)
      complex*16 qa(n), z(n), eye, zdis, z1, z2, zarg
c
      REAL*4 TIMEP(2), ETIME
c
         tbeg = etime(timep)
         do i = 1, n
            phi(i) = 0.d0
            x(i) = dreal(z(i))
            y(i) = dimag(z(i))
            z1 = 1 + (cdabs(z(i)))**2
            do j = 1, i-1
               z2 = 1 + (cdabs(z(j)))**2
               zdis = z(j)-z(i)
               zarg = zdis*dconjg(zdis)/(z1*z2)
               phi(i) = phi(i) + real(qa(j)*dlog(cdabs(zarg)))
            end do
            do j = i+1, n
               z2 = 1 + (cdabs(z(j)))**2
               zdis = z(j)-z(i)
               zarg = zdis*dconjg(zdis)/(z1*z2)
               phi(i) = phi(i) + real(qa(j)*dlog(cdabs(zarg)))
            end do
         end do
c
         tend = etime(timep)
         call PRIN2 ('TIME FOR DIRECT POTENTIAL = *',tend-tbeg,1)
c
      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine EXACT_FIELD (n, nsp, x, y, z, qa, zfield)
c---------------
      implicit real*8 (a-h,o-z)
      dimension x(n), y(n)
      complex*16 qa(n), z(n), zfield(n), eye, zdis, z1, z2, zarg
c
      REAL*4 TIMEP(2), ETIME
c
         do i = 1, n
            zfield(i) = 0.d0
            x(i) = dreal(z(i))
            y(i) = dimag(z(i))
            do j = 1, i-1
               zfield(i) = zfield(i) - dconjg(qa(j)/(z(j) - z(i)))
            end do
            do j = i+1, n
               zfield(i) = zfield(i) - dconjg(qa(j)/(z(j) - z(i)))
            end do
         end do
         call PRIN2 (' zfield direct = *', zfield, 2*n)
c
         tend = etime(timep)
c
      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine GRADIENT_EX (n, x1, x2, x3, z, qa, phi_1, phi_2, 
     1                        phi_3, phi_1e, phi_2e, phi_3e)
c---------------
      implicit real*8 (a-h,o-z)
      dimension x1(n), x2(n), x3(n), phi_1(n), phi_2(n), phi_3(n), 
     1          phi_1e(n), phi_2e(n), phi_3e(n)
      complex*16 z(n), qa(n), eye, zdis, z1, z2, zarg, dphi_dz
c
         do i = 1, n
            phi_1(i) = 0.d0
            phi_2(i) = 0.d0
            phi_3(i) = 0.d0
            do j = 1, i-1
               r = (x1(i)-x1(j))**2 + (x2(i)-x2(j))**2 +
     1             (x3(i)-x3(j))**2
               phi_1(i) = phi_1(i) + qa(j)*(x1(j)-x1(i))/r
               phi_2(i) = phi_2(i) + qa(j)*(x2(j)-x2(i))/r
               phi_3(i) = phi_3(i) + qa(j)*(x3(j)-x3(i))/r
            end do
            do j = i+1, n
               r = (x1(i)-x1(j))**2 + (x2(i)-x2(j))**2 +
     1             (x3(i)-x3(j))**2
               phi_1(i) = phi_1(i) + qa(j)*(x1(j)-x1(i))/r
               phi_2(i) = phi_2(i) + qa(j)*(x2(j)-x2(i))/r
               phi_3(i) = phi_3(i) + qa(j)*(x3(j)-x3(i))/r
            end do
         end do
         call PRIN2 (' dphi dx1 1 = *', phi_1, n)
c
         do i = 1, n
            phi_1e(i) = 0.d0
            phi_2e(i) = 0.d0
            phi_3e(i) = 0.d0
            do j = 1, i-1
               dz_dx1 = 1.d0/(1.d0-x3(j))
               dphi_dz = 1.d0/(z(j)-z(i)) 
     1                    - dconjg(z(j))/(1.d0+(cdabs(z(j)))**2)
               phi_1e(i) = phi_1e(i) 
     1                      + 2.d0*dreal(qa(j)*dphi_dz)*dz_dx1
            end do
            do j = i+1, n
               dz_dx1 = 1.d0/(1.d0-x3(j))
               dphi_dz = 1.d0/(z(j)-z(i)) 
     1                    - dconjg(z(j))/(1.d0+(cdabs(z(j)))**2)
               phi_1e(i) = phi_1e(i) 
     1                      + 2.d0*dreal(qa(j)*dphi_dz)*dz_dx1
            end do
         end do
         call PRIN2 (' dphi dx1 2 = *', phi_1e, n)
ccc         call PRIN2 (' dphi dx2 = *', phi_2, n)
ccc         call PRIN2 (' dphi dx3 = *', phi_3, n)
c
         tend = etime(timep)
c
      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine GRAD_BUILD (n, x1, x2, x3, z, cfield, phi_1, phi_2, 
     1                       phi_3)
c---------------
      implicit real*8 (a-h,o-z)
      dimension x1(n), x2(n), x3(n), phi_1(n), phi_2(n), phi_3(n)
      complex*16 z(n), cfield(n), eye, zdis, z1, z2, zarg
c
         eye = dcmplx(0.d0,1.d0)
c
         do i = 1, n
            phi_1(i) = 2.d0*dreal(cfield(i))/(1.d0-x3(i))
            phi_2(i) = 2.d0*dreal(eye*cfield(i))/(1.d0-x3(i))
            phi_3(i) = -2.d0*dreal(z(i)*cfield(i))/(1.d0-x3(i))
         end do
         call PRIN2 (' dphi dx1 = *', phi_1, n)
         call PRIN2 (' dphi dx2 = *', phi_2, n)
         call PRIN2 (' dphi dx3 = *', phi_3, n)
c
         tend = etime(timep)
c
      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine  FASMVP (k, nd, nbk, nsp, ds, x_zeta, y_zeta, zeta,    
     1                    dzeta, zeta_k, diag, A_k, u, w, qa,   
     2                    cfield, poten, wksp)
c---------------
c
      implicit real*8 (a-h, o-z)
      integer*4 iout(2), inform(10), ierr(10)
      complex*16 zeta(nbk), dzeta(nbk), qa(nbk), cfield(nbk), 
     1           zQsum, zQ2sum, eye, zeta_k(k), zdis
      dimension u(*), w(*), poten(nbk), wksp(nsp), x_zeta(nbk), 
     1          y_zeta(nbk), diag(nbk), dsda(nbk), A_k(k), ds(k),
     2          nd(k)
      REAL*4 TIMEP(2), ETIME
c
         pi = 4.D0*DATAN(1.D0)
         eye = DCMPLX(0.D0,1.D0)
c
c Extract off A_k
         do kbod = 1, k
            A_k(kbod) = u(nbk+kbod)
         end do
c
         zQsum = 0.d0
         istart = 0
         do kbod = 1, k
            do i = 1, nd(kbod)
               qa(istart+i) = ds(kbod)*u(istart+i)*dzeta(istart+i)
     1                         /(2.d0*pi)
               zQ2sum = ds(kbod)*u(istart+i)*dzeta(istart+i)*
     1                   dconjg(zeta(istart+i))
     2                  /(1.d0+cdabs(zeta(istart+i))**2)
               zQsum = zQsum - zQ2sum/(2.d0*pi)
            end do
            istart = istart+nd(kbod)
         end do
ccc         call PRIn2 (' zQsum = *', zQsum, 2)
c
         tbeg = etime(timep)
         iout(1) = 0
         iout(2) = 13
         iflag7 = 3
         napb = 20
         ninire = 2
         mex = 300
         eps7 = 1.d-14
         tol = 1.d-14
         nnn = nbk
         call DAPIF2 (iout, iflag7, nnn, napb, ninire, mex, ierr, 
     &                inform, tol, eps7, x_zeta, y_zeta, qa, poten,  
     &                cfield, wksp, nsp, CLOSE)
         call PRINI (6, 13)
ccc         call PRIN2 (' qa = *', qa, 2*nnn)
ccc         call PRIN2 (' cfielf = *', cfield, 2*nnn)
         if (ierr(1).ne.0) then
            write (6,*) '  ERROR IN DAPIF2, IERR = ', (ierr(ii),ii=1,6)
            write(6,*) '  INFORM = ', (inform(ii),ii=1,6)
            stop
         end if
ccc         call PRINF (' Number of Levels used = *', inform(3), 1)
ccc         call PRIN2 (' cfield = *', cfield, 2*nbk)
c
c Fix up field
         istart = 0
         do kbod = 1, k
            do i = 1, nd(kbod)
               zQ2sum = ds(kbod)*u(istart+i)*dzeta(istart+i)
     1                     *dconjg(zeta(istart+i))
     1                      /(1.d0+cdabs(zeta(istart+i))**2)
               zQ2sum = - zQ2sum/(2.d0*pi)
               cfield(istart+i) = cfield(istart+i) - zQsum + zQ2sum
               w(istart+i) = 0.5d0*u(istart+i) 
     1                       + dimag(cfield(istart+i)) 
     2                       - ds(kbod)*diag(istart+i)*u(istart+i)
c
c Add on log singularities
               do mbod = 1, k
                  zdis = zeta(istart+i) - zeta_k(mbod)
                  rad = 2.d0*(cdabs(zdis))**2
     1                      /((1+(cdabs(zeta(istart+i)))**2)
     2                     *((1+(cdabs(zeta_k(mbod)))**2)))
                  w(istart+i) = w(istart+i) + A_k(mbod)*0.5d0*dlog(rad)
               end do
            end do
            istart = istart+nd(kbod)
         end do
c
c
c Constraints for multiple log but with same-valued streamlines
         w(nbk+1) = 0.d0
         do kbod = 1, k
            w(nbk+1) = w(nbk+1) + A_k(kbod)
         end do
         istart = nd(1)
         do kbod = 2, k
            w(nbk+kbod) = 0.d0
            do i = istart+1, istart+nd(kbod)
               w(nbk + kbod) = w(nbk + kbod) + u(i)
            end do
            istart = istart+nd(kbod)
         end do 
         tend = etime(timep)
ccc         call PRIN2 (' poten = *', poten, n)
ccc         call PRIN2 (' TIME FOR FMM  = *',tend-tbeg,1)
ccc         call PRIN2 (' cfield = *', cfield, 2*n)
c
      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine  BUILD_MAT_STEREO (k, nd, nbk, x_zeta, y_zeta,     
     1                              zeta, dzeta, zeta_k, diag, amat)
c---------------
c
      implicit real*8 (a-h, o-z)
      integer*4 iout(2), inform(10), ierr(10)
      complex*16 zeta(nbk), dzeta(nbk),  
     1           zQsum, zQ2sum, eye, zeta_k(k), zdis, zsrc
      dimension x_zeta(nbk), 
     1          y_zeta(nbk), diag(nbk), amat(nbk+k,nbk+k)
      REAL*4 TIMEP(2), ETIME
c
         pi = 4.D0*DATAN(1.D0)
         eye = DCMPLX(0.D0,1.D0)
         dalph = 2.d0*pi/nd
c
ccc         call PRIn2 (' zQsum = *', zQsum, 2)
c
c Fix up field
         do i = 1, nbk
            do j = 1, nbk
               if (i.ne.j) then
                  zdis = zeta(i)-zeta(j)
                  zsrc = dzeta(j)/zdis
                  zsrc = dalph*zsrc/(2.d0*pi*eye)
                  zQsum = dzeta(j)*dconjg(zeta(j))
     1                   /(1.d0+cdabs(zeta(j))**2)
                  zQsum = dalph*zQsum/(2.d0*pi*eye)
                  amat(i,j) = dreal(zsrc+zQsum)
                else
                  amat(i,i) = 0.5d0-dalph*diag(i) 
               end if
            end do
c
c Add on log singularities
            do kbod = 1, k
ccc            do kbod = 1, 1
               zdis = zeta(i) - zeta_k(kbod)
               rad = 2.d0*(cdabs(zdis))**2/((1+(cdabs(zeta(i)))**2)
     1                  *((1+(cdabs(zeta_k(kbod)))**2)))
               amat(i,nbk+kbod)= 0.5d0*dlog(rad)
            end do
         end do
c
c
c Constraints for multiple log but with same-valued streamlines
         do kbod = 1, k
            amat(nbk+1,nbk+kbod) = 1.d0 
         end do
         do kbod = 2, k
            do i = (kbod-1)*nd+1, kbod*nd
ccc               amat(nbk+kbod,i) = cdabs(dzeta(i))
               amat(nbk+kbod,i) = 1.d0
            end do
         end do 
c
c Dump it out
         open (unit = 24, file = 'amat_stereo.dat')
         do i = 1, nbk+k
            do j = 1, nbk+k
               write(24,'(e20.13,$)')(amat(i,j))
               write (24,'(a)')  ''
            end do
         end do
         close (24) 
c
      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine  FASMVP_TEST (k, nd, nbk, nsp, x_zeta, y_zeta, zeta,    
     1                    dzeta, zeta_k, diag, dsda, A_k, u, w, qa,   
     2                    cfield, poten, wksp)
c---------------
c
      implicit real*8 (a-h, o-z)
      integer*4 iout(2), inform(10), ierr(10)
      complex*16 zeta(nbk), dzeta(nbk), qa(nbk), cfield(nbk), 
     1           zQsum, zQ2sum, eye, zeta_k(k), zdis
      dimension u(*), w(*), poten(nbk), wksp(nsp), x_zeta(nbk), 
     1          y_zeta(nbk), diag(nbk), dsda(nbk), A_k(k)
      REAL*4 TIMEP(2), ETIME
c
         pi = 4.D0*DATAN(1.D0)
         eye = DCMPLX(0.D0,1.D0)
         dalph = 2.d0*pi/nd
c
c Extract off A_k
         do kbod = 1, k
            A_k(kbod) = u(nbk+kbod)
         end do
c
         zQsum = 0.d0
         do i = 1, nbk
            qa(i) = dalph*u(i)*dzeta(i)/(2.d0*pi)
            zQ2sum = dalph*u(i)*dzeta(i)*dconjg(zeta(i))
     1                   /(1.d0+cdabs(zeta(i))**2)
            zQsum = zQsum - zQ2sum/(2.d0*pi)
         end do
ccc         call PRIn2 (' zQsum = *', zQsum, 2)
c
         tbeg = etime(timep)
         iout(1) = 0
         iout(2) = 13
         iflag7 = 3
         napb = 20
         ninire = 2
         mex = 300
         eps7 = 1.d-14
         tol = 1.d-14
         nnn = nbk
         call DAPIF2 (iout, iflag7, nnn, napb, ninire, mex, ierr, 
     &                inform, tol, eps7, x_zeta, y_zeta, qa, poten,  
     &                cfield, wksp, nsp, CLOSE)
         call PRINI (6, 13)
ccc         call PRIN2 (' qa = *', qa, 2*nnn)
ccc         call PRIN2 (' cfielf = *', cfield, 2*nnn)
         if (ierr(1).ne.0) then
            write (6,*) '  ERROR IN DAPIF2, IERR = ', (ierr(ii),ii=1,6)
            write(6,*) '  INFORM = ', (inform(ii),ii=1,6)
            stop
         end if
ccc         call PRINF (' Number of Levels used = *', inform(3), 1)
ccc         call PRIN2 (' cfield = *', cfield, 2*nbk)
c
c Fix up field
         do i = 1, nbk
            zQ2sum = dalph*u(i)*dzeta(i)*dconjg(zeta(i))
     1                   /(1.d0+cdabs(zeta(i))**2)
            zQ2sum = - zQ2sum/(2.d0*pi)
            cfield(i) = cfield(i) - zQsum + zQ2sum
            w(i) = 0.5d0*u(i) 
c
c Add on log singularities
            do kbod = 1, k
ccc            do kbod = 1, 1
               zdis = zeta(i) - zeta_k(kbod)
ccc               rad = 4.d0*(cdabs(zdis))**2/((1+(cdabs(zeta(i)))**2)
ccc     1                  *((1+(cdabs(zeta_k(kbod)))**2)))
ccc               w(i) = w(i) + A_k(kbod)*0.5d0*dlog(rad)
               rad = 2.d0*(cdabs(zdis))**2/((1+(cdabs(zeta(i)))**2)
     1                  *((1+(cdabs(zeta_k(kbod)))**2)))
               w(i) = w(i) + A_k(kbod)*0.5d0*dlog(rad)
            end do
         end do
c
c
c Constraints for multiple log but with same-valued streamlines
         w(nbk+1) = 0.d0
         do kbod = 1, k
            w(nbk+1) = w(nbk+1) + A_k(kbod)
         end do
         do kbod = 2, k
            w(nbk+kbod) = 0.d0
            do i = (kbod-1)*nd+1, kbod*nd
               w(nbk + kbod) = w(nbk + kbod) + u(i)
            end do
         end do 
         tend = etime(timep)
ccc         call PRIN2 (' poten = *', poten, n)
ccc         call PRIN2 (' TIME FOR FMM  = *',tend-tbeg,1)
ccc         call PRIN2 (' cfield = *', cfield, 2*n)
c
      return
      end
c
c
c---------------
      subroutine SOLVE (nd, k, kmax, nbk, rhs, soln, u, A_k, rwork, 
     1                  lrwork, iwork, liwork, maxl, schur, bnew,
     2                  ipvtbf, z, zk)
c---------------
c
      implicit real*8 (a-h,o-z)
      external MATVEC_LAPL, MSOLVE
c
c  System
      dimension nd(k), soln(*), rhs(*), u(nbk), A_k(k)
c
c  DGMRES work arrays
      dimension rwork(lrwork),iwork(liwork)
c
c  Pivots
      dimension schur(kmax*kmax),bnew(kmax),ipvtbf(kmax)
c
      complex*16 z(nbk), zk(k)
c
c  Timings
c
      real*4 timep(2), etime
c
         pi = 4.d0*datan(1.d0)
c
c  solve linear system using GMRES.
c
c
c     parameters for DGMRES
         itol = 0
         tol = 1.0d-10
         isym = 0
         iwork(1) = maxl
         do i=2,liwork
            iwork(i) = 0
         enddo
c
c  Preconditioner flag
c
ccc         iwork(4) = -1
         iwork(4) = 0
c
c  Restart flag
c  
         iwork(5) = 5      
c
c  factor preconditioner
ccc         t0 = etime(timep)
ccc         call SCHUR_FACTOR (z,ND,nbk,K,zk,SCHUR,bnew,IPVTBF)
ccc         t1 = etime(timep)
         call prin2 (' time in preconditioner = *', t1-t0, 1)
c
c     provide initial guess soln
         norder = nbk + k
ccc         norder = nbk
         do i=1,norder
            soln(i) = rhs(i)
         enddo
c
         t0 = etime(timep)
         call DGMRES (norder, rhs, soln, nelt, ia, ja, a, isym,
     1               MATVEC_LAPL, MSOLVE, itol, tol, itmax, iter, err,  
     1               ierr, 6, sb, sx, rwork, lrwork, iwork, 
     1               liwork, rw, iw)
         call Prin2 (' after laplace solve, err = *', err, 1)
         call PrinF (' after laplace solve, ierr = *', ierr, 1)
         call PRINI (6,13)
         call PRINF ('  # GMRES ITERATIONS = *',iter,1)
         if (ierr.gt.2) then
            call PRINF ('  SOMETHING WRONG IN GMRES, IERR = *',ierr,1)
            call PRINF ('  iwork = *',iwork,10)
            stop
           elseif (ierr.ge.0) then
            t1 = etime(timep)
            tsec = t1 - t0
            call PRIn2 (' time in solve = *', tsec, 1)
c
c  unpack soln into U
            istart = 0
            do i = 1, nbk
               u(i) = soln(i)
            end do  
ccc            call PRIn2 (' density = *', u, nbk)
            istart = 0
            do kbod = 1, k
ccc               call PRINF (' kbod = *', kbod, 1)
ccc               call PRIn2 ('     u = *', u((kbod-1)*nd(kbod-1)+1), 
ccc     1                                   nd(kbod))
               A_k(kbod) = soln(nbk+kbod)
ccc               utest = 0.d0
ccc               do i = 1, nd(kbod)
ccc                  utest = utest + u(istart+i)
ccc               end do
ccc               call prin2 (' this should be zero = *', utest, 1)
ccc               istart = istart + nd(1)
            end do 
            call PRIN2 (' A_k = *', A_k, k)
         end if
c
c Dump it out
         open (unit = 24, file = 'solution.dat')
         do i = 1, nbk+k
               write(24,'(e20.13,$)')(soln(i))
               write (24,'(a)')  ''
         end do
         close (24) 
c
      return
      end
c
c*************************************************
c
      subroutine RESAMPLE (ns,nsk,k0,kk,nsamp,h,z,dz,xat,yat,u,x,y,
     *                     wksp,nsp)
c
c  Overresolve the data on the boundary
c
      implicit real*8 (a-h,o-z)
      dimension h(k0:kk),ns(k0:kk),xat(*),yat(*),x(*),y(*)
      complex*16 z(*),dz(*),u(*),wksp(nsp)
c
         pi = 4.d0*datan(1.d0)
c
c  do z first
c
         istart = 0
         istart2 = 0
         do nbod = k0,kk
            nd2 = nsamp*ns(nbod)
            ndm1 = ns(nbod)-1
            nd2m1 = nd2-1
            call FTRPIN (wksp,nsp,ip2,ipt,ndm1,nd2m1)
            call FINTER (xat(istart+1),x(istart2+1),ndm1,nd2m1,wksp,
     *                   nsp,ip2,ipt)
            call FINTER (yat(istart+1),y(istart2+1),ndm1,nd2m1,wksp,
     *                   nsp,ip2,ipt)
            do i = 1,nsamp*ns(nbod)
               z(istart2+i) = dcmplx(x(istart2+i),y(istart2+i))
            end do
            istart = istart + ns(nbod)
            istart2 = istart2 + nsamp*ns(nbod)
         end do
c
c  now do dz
c
         do i = 1,nsk
            xat(i) = dreal(dz(i))
            yat(i) = dimag(dz(i))
         end do
         istart = 0
         istart2 = 0
         do nbod = k0,kk
            nd2 = nsamp*ns(nbod)
            ndm1 = ns(nbod)-1
            nd2m1 = nd2-1
            call FTRPIN (wksp,nsp,ip2,ipt,ndm1,nd2m1)
            call FINTER (xat(istart+1),x(istart2+1),ndm1,nd2m1,wksp,
     *                   nsp,ip2,ipt)
            call FINTER (yat(istart+1),y(istart2+1),ndm1,nd2m1,wksp,
     *                   nsp,ip2,ipt)
            do i = 1,nsamp*ns(nbod)
               dz(istart2+i) = dcmplx(x(istart2+i),y(istart2+i))
            end do
            istart = istart + ns(nbod)
            istart2 = istart2 + nsamp*ns(nbod)
         end do
c
c  now do u
c
         do i = 1,nsk
            xat(i) = dreal(u(i))
            yat(i) = dimag(u(i))
         end do
         istart = 0
         istart2 = 0
         do nbod = k0,kk
            nd2 = nsamp*ns(nbod)
            ndm1 = ns(nbod)-1
            nd2m1 = nd2-1
            call FTRPIN (wksp,nsp,ip2,ipt,ndm1,nd2m1)
            call FINTER (xat(istart+1),x(istart2+1),ndm1,nd2m1,wksp,
     *                   nsp,ip2,ipt)
            call FINTER (yat(istart+1),y(istart2+1),ndm1,nd2m1,wksp,
     *                   nsp,ip2,ipt)
            do i = 1,nsamp*ns(nbod)
               u(istart2+i) = dcmplx(x(istart2+i),y(istart2+i))
            end do
            istart = istart + ns(nbod)
            istart2 = istart2 + nsamp*ns(nbod)
         end do
c
c  Update points and stuff
c
         nsk = nsamp*nsk
         do nbod = k0,kk
            ns(nbod) = nsamp*ns(nbod)
            h(nbod) = 2.d0*pi/ns(nbod)
         end do
         do i = 1,nsk
            xat(i) = dreal(z(i))
            yat(i) = dimag(z(i))
         end do
c
      return
      end
c
c
c---------------
      subroutine SOL_GRID (nd, k, nbk, nth, nphi, u, A_k, zeta_k, zeta,  
     1                     dzeta, igrid, zeta_gr, u_gr)
c---------------
c
      implicit real*8 (a-h,o-z)
      dimension u(nbk), igrid(nth,nphi), u_gr(nth,nphi), A_k(k)
      complex*16 zeta(nbk), dzeta(nbk), zeta_gr(nth,nphi), zkern, 
     1           zeta_k(k), zdis, eye
      REAL*4 TIMEP(2), ETIME
c
         pi = 4.d0*datan(1.d0)
         eye = dcmplx(0.d0, 1.d0)
         dalph = 2.d0*pi/nd
c
         tbeg = etime(timep)
         do i = 1, nth
            do j = 1, nphi
               u_gr(i,j) = 0.d0
               if (igrid(i,j).ne.0) then 
               do ip = 1, nbk
                  zdis = zeta(ip) - zeta_gr(i,j)
                  zkern = dzeta(ip)/zdis - dconjg(zeta(ip))*dzeta(ip)
     1                         /(1.d0+cdabs(zeta(ip))**2)
                  u_gr(i,j) = u_gr(i,j) 
     1                   - dalph*u(ip)*dimag(zkern)/(2.d0*pi)
               end do
               do kbod = 1, k
                  zdis = zeta_gr(i,j) - zeta_k(kbod)
                  rad = 
     1              2.d0*(cdabs(zdis))**2/((1+(cdabs(zeta_gr(i,j)))**2)
     2                  *((1+(cdabs(zeta_k(kbod)))**2)))
                  u_gr(i,j) = u_gr(i,j) + A_k(kbod)*0.5d0*dlog(rad)
               end do
               end if
            end do
         end do
         tend = etime(timep)
ccc         call PRIN2 (' poten = *', poten, n)
         call PRIN2 (' TIME FOR GRID = *',tend-tbeg,1)
c
      return
      end
c
c
c---------------
      subroutine SOL_GRID_FMM (nd, k, nbk, nx, ny, ds, u, A_k, zeta_k,   
     1                         zeta, dzeta, igrid, zeta_gr, 
     2                         u_gr, x_zeta, y_zeta, qa, cfield, poten,  
     3                         nsp, wksp, nvort, vort_k, zk_vort)
c---------------
c
      implicit real*8 (a-h,o-z)
      dimension u(nbk), igrid(nx,ny), u_gr(nx,ny), A_k(k),
     1          vort_k(nvort)
      dimension nd(k), ds(k)
      complex*16 zeta(nbk), dzeta(nbk), zkern, zeta_gr(nx,ny),
     1           zeta_k(k), zdis, eye, qa(*), cfield(*), zQsum, 
     2           zQ2sum, zk_vort(nvort)
      dimension x_zeta(*), y_zeta(*), poten(*), wksp(nsp)
      integer*4 iout(2), inform(10), ierr(10)
      REAL*4 TIMEP(2), ETIME
c
         pi = 4.d0*datan(1.d0)
         eye = dcmplx(0.d0, 1.d0)
c
c pack zeta and zeta_gr into x_zeta and y_zeta
         istart = 0
         zQsum = 0.d0
         do kbod = 1, k
            do i = 1, nd(kbod)
               x_zeta(istart+i) = dreal(zeta(istart+i))
               y_zeta(istart+i) = dimag(zeta(istart+i))
               qa(istart+i) = ds(kbod)*u(istart+i)
     1                            *dzeta(istart+i)/(2.d0*pi)
               zQ2sum = ds(kbod)*u(istart+i)*dzeta(istart+i)
     1                    *dconjg(zeta(istart+i))
     2                    /(1.d0+cdabs(zeta(istart+i))**2)
               zQsum = zQsum - zQ2sum/(2.d0*pi)
            end do
            istart = istart + nd(kbod)
         end do
         ij = nbk
         do i = 1, nx
            do j = 1, ny
               if (igrid(i,j).eq.1) then
                  ij = ij + 1
                  x_zeta(ij) = dreal(zeta_gr(i,j))
                  y_zeta(ij) = dimag(zeta_gr(i,j))
                  qa(ij) = 0.d0
               end if
            end do
         end do
         call PRINF (' in SOL_GRID_FMM, ntar = *', ij-nbk, 1)
ccc         call PRIn2 (' zQsum = *', zQsum, 2)
c
         tbeg = etime(timep)
         iout(1) = 0
         iout(2) = 13
         iflag7 = 3
         napb = 50
         ninire = 2
         mex = 300
         eps7 = 1.d-14
         tol = 1.d-14
         nnn = ij
ccc         nnn = nbk
ccc         nnn = nbk+100
         call DAPIF2 (iout, iflag7, nnn, napb, ninire, mex, ierr, 
     &                inform, tol, eps7, x_zeta, y_zeta, qa, poten,  
     &                cfield, wksp, nsp, CLOSE)
         call PRINI (6, 13)
ccc         call PRIN2 (' cfield = *', cfield, 2*nnn)
         call PRINF (' Number of Levels = *', inform(3), 1)
         if (ierr(1).ne.0) then
            write (6,*) '  ERROR IN DAPIF2, IERR = ', (ierr(ii),ii=1,6)
            write(6,*) '  INFORM = ', (inform(ii),ii=1,6)
            stop
         end if
ccc         call PRINF (' Number of Levels used = *', inform(3), 1)
ccc         call PRIN2 (' cfield = *', cfield, 2*nbk)
c
         call PRIN2 (' a_k in sol_GRID_FMM = *', A_k, k)
c Fix up field
         ij = nbk
         umax = -1.d10
         umin = 1.d10
         do i = 1, nx
            do j = 1, ny
               if (igrid(i,j).eq.1) then     
                  ij = ij + 1           
                  u_gr(i,j) = dimag(cfield(ij) - zQsum)
                  do kbod = 1, k
                     zdis = zeta_gr(i,j) - zeta_k(kbod)
                     rad = 2.d0*
     1                  (cdabs(zdis))**2/((1+(cdabs(zeta_gr(i,j)))**2)
     2                  *((1+(cdabs(zeta_k(kbod)))**2)))
                     u_gr(i,j) = u_gr(i,j) + A_k(kbod)*0.5d0*dlog(rad)
                  end do
                  umax = max(umax,u_gr(i,j))
                  umin = min(umin,u_gr(i,j))
                 else
                  u_gr(i,j) = -1000.d0
               end if
            end do
         end do
         call PRIN2 (' Max solution = *', umax, 1)
         call PRIN2 (' Min solution = *', umin, 1)
c
c Add in vortex singularities
ccc         call PRIN2 (' zk_vort = *', zk_vort, 2*nvort)
         call PRIN2 (' vort_k = *', vort_k, nvort)
         do i = 1, nx
            do j = 1, ny
cccccc               u_gr(i,j) = 0.d0
               if (igrid(i,j).eq.1) then   
                  psi = 0.d0  
                  do ivort = 1, nvort
                     zdis = zeta_gr(i,j) - zk_vort(ivort)
                     az1 = (cdabs(zeta_gr(i,j)))**2
                     az2 = (cdabs(zk_vort(ivort)))**2
                     arg = dreal(zdis*conjg(zdis)
     1                        /((1.d0+az1)*(1.d0+az2)))
                     psi = psi + vort_k(ivort)*dlog(arg)
                  end do
cccccc                  call PRIN2 (' psi = *', psi, 1)
                  u_gr(i,j) = u_gr(i,j) - psi
               end if
            end do
         end do
c
         tend = etime(timep)
ccc         call PRIN2 (' poten = *', poten, n)
         call PRIN2 (' TIME FOR FMM  ON GRID = *',tend-tbeg,1)
ccc         call PRIN2 (' cfield = *', cfield, 2*n)
c
      return
      end
c
c
c---------------
      subroutine SOL_TAR_FMM (nd, k, nbk, ntar, ds, u, A_k, zeta_k,   
     1                        x_zeta, y_zeta, zeta, dzeta, zeta_tar, 
     2                        u_tar, xz_tar, yz_tar, qa, cfield, poten,  
     3                        nsp, wksp, nvort, vort_k, zk_vort)
c---------------
c
      implicit real*8 (a-h,o-z)
      dimension u(nbk), u_tar(ntar), A_k(k), vort_k(nvort), x_zeta(*),
     1          y_zeta(*), nd(k), ds(k)
      complex*16 zeta(nbk), dzeta(nbk), zeta_tar(ntar), zkern, 
     1           zeta_k(k), zdis, eye, qa(*), cfield(*), zQsum, 
     2           zQ2sum, zk_vort(nvort), ztar
      dimension xz_tar(ntar), yz_tar(ntar), poten(*), wksp(nsp)
      integer*4 iout(2), inform(10), ierr(10)
      REAL*4 TIMEP(2), ETIME
c
         pi = 4.d0*datan(1.d0)
         eye = dcmplx(0.d0, 1.d0)
         call PRIN2 (' in sol_tar_fmm, ds = *', ds, k)
         call PRIN2 (' A_k = *', A_k, k)
         call prin2 (' zeta_k = *', zeta_k, 2*k)
c
c pack zeta and zeta_gr into x_zeta and y_zeta
         istart = 0
         do kbod = 1, k
            do i = 1, nd(kbod)
               qa(istart+i) = ds(kbod)*u(istart+i)
     1                           *dzeta(istart+i)/(2.d0*pi)
            end do
            istart = istart + nd(kbod)
         end do
         ij = nbk
         do i = 1, ntar
            ij = ij + 1
            x_zeta(ij) = xz_tar(i)
            y_zeta(ij) = yz_tar(i)
            qa(ij) = 0.d0
         end do
         call PRINF (' ij = *', ij, 1)
c
         zQsum = 0.d0
         istart = 0
         do kbod = 1, k
            do i = 1, nd(kbod)
               zQ2sum = ds(kbod)*u(istart+i)*dzeta(istart+i)
     1                    *dconjg(zeta(istart+i))
     2                      /(1.d0+cdabs(zeta(istart+i))**2)
               zQsum = zQsum - zQ2sum/(2.d0*pi)
            end do
            istart = istart + nd(kbod)
         end do
ccc         call prin2 (' qa in sol_tar_fmm = *', qa, 2*nbk)
ccc         call PRIn2 (' zQsum = *', zQsum, 2)
c
         tbeg = etime(timep)
         iout(1) = 0
         iout(2) = 13
         iflag7 = 3
         napb = 50
         ninire = 2
         mex = 300
         eps7 = 1.d-14
         tol = 1.d-14
         nnn = ij
ccc         nnn = nbk
ccc         nnn = nbk+100
         call DAPIF2 (iout, iflag7, nnn, napb, ninire, mex, ierr, 
     &                inform, tol, eps7, x_zeta, y_zeta, qa, poten,  
     &                cfield, wksp, nsp, CLOSE)
         call PRINI (6, 13)
         call PRIN2 (' cfield = *', cfield(nbk+1), 2*ntar)
         call PRINF (' Number of Levels = *', inform(3), 1)
         if (ierr(1).ne.0) then
            write (6,*) '  ERROR IN DAPIF2, IERR = ', (ierr(ii),ii=1,6)
            write(6,*) '  INFORM = ', (inform(ii),ii=1,6)
            stop
         end if
ccc         call PRINF (' Number of Levels used = *', inform(3), 1)
ccc         call PRIN2 (' cfield = *', cfield, 2*nbk)
c
ccc         call PRIN2 (' a_k in sol_GRID_FMM = *', A_k, k)
c Fix up field
         ij = nbk
         do i = 1, ntar
            ij = ij + 1           
            u_tar(i) = dimag(cfield(ij) - zQsum)
            ztar = dcmplx(xz_tar(i),yz_tar(i))
            do kbod = 1, k
               zdis = ztar - zeta_k(kbod)
                     rad = 
     1              (cdabs(zdis))**2/((1+(cdabs(ztar))**2)
     2                  *((1+(cdabs(zeta_k(kbod)))**2)))
                     u_tar(i) = u_tar(i) + A_k(kbod)*dlog(rad)
            end do
         end do
c
c Add in vortex singularities
ccc         call PRIN2 (' zk_vort = *', zk_vort, 2*nvort)
ccc         call PRIN2 (' vort_k = *', vort_k, nvort)
         do i = 1, ntar
            psi = 0.d0  
            ztar = dcmplx(xz_tar(i),yz_tar(i))
            do ivort = 1, nvort
               zdis = ztar - zk_vort(ivort)
               az1 = (cdabs(ztar))**2
               az2 = (cdabs(zk_vort(ivort)))**2
               arg = dreal(zdis*conjg(zdis)
     1                        /((1.d0+az1)*(1.d0+az2)))
               psi = psi + vort_k(ivort)*dlog(arg)
               u_tar(i) = u_tar(i) - psi
            end do
         end do
         call prin2 (' u_tar = *', u_tar, ntar)
c
         tend = etime(timep)
ccc         call PRIN2 (' poten = *', poten, n)
         call PRIN2 (' TIME FOR FMM  ON GRID = *',tend-tbeg,1)
ccc         call PRIN2 (' cfield = *', cfield, 2*n)
c
      return
      end
c
c
c---------------
      subroutine SOL_VORT_CHECK (nd, k, nbk, nth, nphi, nvort, zeta_gr,  
     1                           igrid, u_gr, uex_gr, zk_vort, r0)
c---------------
c
      implicit real*8 (a-h,o-z)
      dimension u_gr(nth,nphi), uex_gr(nth,nphi), igrid(nth,nphi)
      complex*16 zeta_gr(nth,nphi), zk_vort(nvort), zet, z, zvort, eye,
     1           zeta_vort, zdis
c
         pi = 4.d0*datan(1.d0)
         eye = dcmplx(0.d0, 1.d0)
c
         err = 0.d0
         call PRIN2 (' zk_vort = *', zk_vort, 2*nvort)
         do i = 1, 1
            do j = 1, 1
ccc               if (igrid(i,j).eq.1) then 
                  z = zeta_gr(i,j)
                  zeta_vort = eye*(r0-zk_vort(1))/(r0+zk_vort(1))
                  zet = eye*(r0-z)/(r0+z)
                  arg = cdabs((zet-zeta_vort)/(zet-dconjg(zeta_vort)))
                  uex_gr(i,j) = -dlog(arg)
                  err = max(err,dabs(uex_gr(i,j)-u_gr(i,j)))
ccc                  call PRIN2 (' uex = *', uex_gr(i,j), 1)
ccc                  call PRIn2 (' u = *', u_gr(i,j), 1)
ccc               end if
            end do
         end do
         call PRIN2 (' error in vorticity cap solution = *', err, 1)
c
      return
      end
c
c
c---------------
      subroutine CALC_VEL (k, nd, nbk, nvort, u, A_k, zeta, dzeta, 
     1                     zeta_k, vort_k, zk_vort, zvel, zf, wsave)
c---------------
c
      implicit real*8 (a-h,o-z)
      dimension A_k(k), u(nbk), vort_k(nvort)
      complex*16 zeta(nbk), dzeta(nbk), zeta_k(k), 
     1           zk_vort(k), zf(*), zvel(nvort), wsave(*), eye
c
         pi = 4.d0*datan(1.d0)
         eye = dcmplx(0.d0, 1.d0)
         dalph = 2.d0*pi/nd
c
         istart = 0
         do ivort = 1, nvort
            zvel(ivort) = 0.d0
            vfactor = 1.d0 + (cdabs(zk_vort(ivort)))**2
            do kbod = 1, k
               do i = 1, nd
                  zf(i) = u(istart+i)*dzeta(istart+i)/
     1                     (zk_vort(ivort)-zeta(istart+i))**2
ccc                  zf(i) = u(istart+i)
               end do
               call DCFFTF (nd, zf, wsave)
               zvel(ivort) = zvel(ivort) + 0.5d0*eye*zf(1)/nd
               istart = istart + nd
            end do
ccc               do mbod = 1, k
               do mbod = 1, 1
ccc                  zvel(ivort) = zvel(ivort) + 0.5d0*A_k(mbod)*
ccc     1                (1.d0/(zk_vort(ivort) - zeta_k(mbod)) -
ccc     2                 dconjg(zk_vort(ivort))/vfactor)
                  zvel(ivort) = zvel(ivort) + A_k(mbod)*
     1                (1.d0/(zk_vort(ivort) - zeta_k(mbod)) -
     2                 dconjg(zk_vort(ivort))/vfactor)
               end do
               zvel(ivort) = dconjg(-eye*0.5*zvel(ivort)*vfactor**2)
          end do
         call PRIN2 (' zvel = *', zvel, 2*nvort)
c
      return
      end
c
c
c---------------
      subroutine CHECK_ERROR_GRID (k, nx, ny, zeta_k, zeta_gr, igrid, 
     1                             u_gr)
c---------------
c
      implicit real*8 (a-h,o-z)
      dimension igrid(nx,ny), u_gr(nx,ny)
      complex*16 zeta_gr(nx,ny), zeta_k(k), zdis, eye
c
         err = 0.d0
         umax = 0.d0
         do i = 1, nx
            do j = 1, ny
               if (igrid(i,j).eq.1) then
               u_ex = 0.d0
               do mbod = 1, k
                  u_ex = u_ex + 1.d0/(zeta_gr(i,j)-zeta_k(mbod)) 
     1                   + dconjg(1.d0/(zeta_gr(i,j)-zeta_k(mbod)))
               end do
               err = max(err,dabs(u_ex-u_gr(i,j)))
               umax = max(umax, dabs(u_ex))
               end if
            end do
         end do
         call PRIN2 ('abs error on grid = *', err, 1)
         call PRIN2 ('rel error on grid = *', err/umax, 1)
c
c
      return
      end
c
c
c---------------
      subroutine CHECK_ERROR_TAR (nd, k, nbk, ntar, zeta_k, zeta_tar,  
     1                            u_tar)
c---------------
c
      implicit real*8 (a-h,o-z)
      dimension u_tar(ntar), nd(k)
      complex*16 zeta_tar(ntar), zkern, zeta_k(k), zdis, eye
c
         pi = 4.d0*datan(1.d0)
         eye = dcmplx(0.d0, 1.d0)
c
ccc         call PRIn2 (' zeta_k in check_ERROR = *', zeta_k, 2*k)
ccc         call PRIN2 (' zeta_tar = *', zeta_tar, 2)
c         call PRINF (' In CHECK_ERROR_TAR, NTAR = *', ntar, 1)
         umax = 0.d0
         err = 0.d0
         do i = 1, ntar
            u_ex = 0.d0
            do mbod = 1, k
               u_ex = u_ex + dreal(1.d0/(zeta_tar(i)-zeta_k(mbod)) 
     1                + dconjg(1.d0/(zeta_tar(i)-zeta_k(mbod))))
            end do
            umax = max(umax, dabs(u_ex))
ccc               u_ex = u_ex + 66.d0 
            err = max(err,dabs(u_ex-u_tar(i)))
ccc            call PRIN2 ('### u_ex = *', u_ex, 1)
ccc            call PRIN2 ('    u_tar  = *', u_tar(i), 1)
         end do
         call PRIN2 (' abs error in targets = *', err, 1)
         call PRIN2 (' rel error in targets = *', err/umax, 1)
c
c
      return
      end
c
C********************************************************************
      SUBROUTINE MSOLVE(N, R, Z, NELT, IA, JA, A, ISYM, RWORK, IWORK)
C********************************************************************
c
c     Another routine required by DGMRES. It allows the use of a
c     preconditioner.
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c
      DIMENSION R(N), Z(N)
c
      parameter (kmax = 5, npmax = 75000, nmax = kmax*npmax)
      parameter (nth_max = 1000, nphi_max = 1000, 
     1          ng_max = nth_max*nphi_max)
      parameter (nsp = 20*nmax + 20*ng_max)
c      
      common /geometry/ x_zeta, y_zeta, zeta, dzeta, dsda
      common /inteqn/ diag, cx, cy, cz
      common /sys_size/ k, nd, nbk
      common /fasblk2/ schur,wb,ipvtbf
      dimension schur(kmax*kmax),wb(kmax),ipvtbf(kmax)
      dimension x_zeta(nmax+ng_max), y_zeta(nmax+ng_max), dsda(nmax),  
     1          diag(nmax), cx(kmax), cy(kmax), cz(kmax)
      complex*16 zeta(nmax), dzeta(nmax), zeta_k(kmax), eye
c
         eye = dcmplx(0.d0,1.d0)
         do kbod = 1, k
            zeta_k(kbod) = (cx(kbod) + eye*cy(kbod))/(1.d0-cz(kbod))
         end do
c
         job = 0
         call SCHUR_APPLY (zeta,ND,nbk,K,zeta_k,r,z,JOB,
     1                     SCHUR,k,wb,IPVTBF)
c
      RETURN
      END
c
c*********************************************
c
      SUBROUTINE SCHUR_APPLY (z,ND,nbk,K,zk,U,W,JOB,
     1                        SCHUR,LDS,BNEW,IPVTBF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER *4 IPVTBF(k)
      DIMENSION U(*),W(*)
      DIMENSION SCHUR(k,k),BNEW(k)
      complex*16 zk(k), z(nbk), zdis
c
c     For exterior problem
C
C
C     This subroutine applies the PRECONDITIONER matrix 
C     associated with the modified boundary integral approach
C     to the exterior Dirichlet problem in a fast manner.
C     In the new integral equation formulation derived from
C     Mikhlin, the discrete system takes the form  Mx = b where
C     
C           M = |Mp | C |   Mp is N x N (N = total # bdry points)
C               |---|---|   C is N x k, R is k x N and L is k x k.
C               | R | L |
C
C     Mp is of the form I + Q where Q is an integral operator 
C     with continuous kernel. A good preconditioner, therefore,
C     is
C     
C           B = | I | C |   I is N x N (N = total # bdry points)
C               |---|---|   C is N x k, R is k x N and L is k x k.
C               | R | L |
C
C     The vector U is ordered as 1) dipole densities 2) coeffs Ak.
C
C ********************    ON INPUT: ****************************
C
C     X = x-coordinates of boundary discretization points
C     Y = y-coordinates of boundary discretization points
C     ND(i) = number of points in discretization of ith body
C     NDMAX = max number of points in discretization of each body
C     K = number of bodies
C
C     CX,CY = coordinates of source point inside Kth body
C     H = mesh width on Kth body
C     U = given vector to which matrix is to be applied
C     JOB  INTEGER
C          = 0 to solve B*W = U 
C          = NONZERO to solve TRANSPOSE(B)*W= U.
C
C ********************    ON OUTPUT: ****************************
C
C     if (JOB .EQ. 0) then W = B^(-1)*U
C     if (JOB .NEQ. 0) then W = B^(-T)*U
C-----------------------------------------------------------------
C
C     begin by doing Gausian elimination on block R. 
C     Each of the first (K-1) rows of Block R is composed of 
C     o single sequence of 1's, of length ND(K). (see paper).
C
C     form rhs corresponding to Schur complement.
C
c
       ISTART = nd
       bnew(1) = u(nbk+1)
       DO KBOD = 2,K
         SUM1 = 0.0d0
         DO i = 1,ND
             SUM1 = SUM1 + U(ISTART+i)
         end do
         BNEW(KBOD) = U(nbk+kbod)-2.d0*SUM1
         ISTART = ISTART+ND
      end do
c
C
C     solve for coefficients Ak.
C
      CALL DGESL(SCHUR,k,k,IPVTBF,BNEW,JOB)
      DO 1700 KBOD = 1,K
         W(nbk+kbod) = BNEW(kbod)
1700  CONTINUE
C
C     solve for remaining variables.
C
         ISTART = 0
         DO 2500 NBOD = 1,K
         DO 2400 I = 1,ND
            SUM1 = 0.0d0
	    DO 2350 KBOD = 1,K
               dis = cdabs(z(istart+i) - zk(kbod))
	       arg_log = 2.d0*dis**2/((1+cdabs(z(istart+i))**2)
     1                          *(1+cdabs(zk(kbod))**2))
               SUM1 = SUM1 + 0.5d0*dlog(arg_log)*bnew(kbod)
2350        CONTINUE
            W(ISTART+i) = 2.d0*U(ISTART+i) - 2.d0*SUM1
2400     CONTINUE
         ISTART = ISTART+ND
2500     CONTINUE
C
ccc      call prin2(' w is *',w,norder)
      RETURN
      END
c
c*********************************************
c
      SUBROUTINE SCHUR_FACTOR (z,ND,NMAX,K,zk,SCHUR,WB,IPVTBF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER *4 IPVTBF(k)
      DIMENSION SCHUR(k,k),WB(k)
      complex*16 zk(k), z(nmax), zdis
C
c
c     For exterior problem
C
C     This subroutine factors the Schur complement
C     of the PRECONDITIONER matrix 
C     associated with the modified boundary integral approach.
C     In the new integral equation formulation derived from
C     Mikhlin, the discrete system takes the form  Mx = b where
C     
C           M = |Mp | C |   Mp is N x N (N = total # bdry points)
C               |---|---|   C is N x k, R is k x N and L is k x k.
C               | R | L |
C
C     Mp is of the form I + Q where Q is an integral operator 
C     with continuous kernel. A good preconditioner, therefore,
C     is
C     
C           B = | I | C |   I is N x N (N = total # bdry points)
C               |---|---|   C is N x k, R is k x N and L is k x k.
C               | R | L |
C
C     The vector U is ordered as 1) dipole densities 2) coeffs Ak.
C
C ********************    ON INPUT: ****************************
C
C     X = x-coordinates of boundary discretization points
C     Y = y-coordinates of boundary discretization points
C     ND(i) = number of points in discretization of ith body
C     NDMAX = max number of points in discretization of each body
C     K = number of bodies
C
C     CX,CY = coordinates of source point inside Kth body
C     H = mesh width on Kth body
C     ZWB -s a work array of dimension LDS
C
C ********************    ON OUTPUT: ****************************
C
C     SCHUR contains the LU factors of the Schur complement
C     IPVTBF contains pivoting information from LINPACK.
C
C-----------------------------------------------------------------
C
C     begin by doing Gausian elimination on block R. 
C     Each of the first (K-1) rows of Block R is composed of 
C     o single sequence of 1's, of length ND(K). (see paper).
C
C     I.e. form the Schur complement and corresponding rhs.
C
      write (6,*) '** PRECONDITIONER  **'
c
ccc      NORDER = ND*K + K
      istart = nd
      DO IBOD = 2,K
      DO KBOD = 1,K
         SUM1 = 0.0d0
         DO i = 1,ND
            dis = cdabs(z(istart+i) - zk(kbod))
	    arg_log = 2.d0*dis**2/((1+cdabs(z(istart+i))**2)
     1                          *(1+cdabs(zk(kbod))**2))
            SUM1 = SUM1 + 0.5d0*dlog(arg_log)
         end do
         SCHUR(IBOD,KBOD) = -2.d0*SUM1
      end do
      istart = istart+nd
      end do
C
C     Next construct last row of Schur complement.
C 
      DO KBOD = 1,K
	 SCHUR(1,KBOD) = 1.0d0
      end do
C 
      CALL DGECO(SCHUR,K,K,IPVTBF,RCOND,WB)
      call prin2(' rcond = *',rcond,1)
      DO KBOD = 1,K
	 call prinf(' column *',KBOD,1)
	 call prin2(' SCHUR = *',schur(1,KBOD),K)
      END DO
      call PRINf (' ipvtbf = *', ipvtbf, K)
      call PRIN2 (' wb = *', wb, k)
C
ccc      call prin2(' w is *',w,norder)
      RETURN
      END
c
c
c---------------
      subroutine MATVEC_LAPL (N, XX, YY, NELT, IA, JA, A, ISYM)
c---------------
c
c  Required by DGMRES with this precise calling sequence.
c  We ignore most of the parameters except N, XX and YY
c  and call the fast multipole routine FASMVP to do the actual
c  work.
c
      implicit double precision (a-h,o-z)
      dimension xx(n), yy(n)
      parameter (kmax = 5, npmax = 75000, nmax = kmax*npmax)
      parameter (nth_max = 1000, nphi_max = 1000, 
     1          ng_max = nth_max*nphi_max)
      parameter (nsp = 20*nmax + 20*ng_max)
c      
      common /geometry/ x_zeta, y_zeta, zeta, dzeta, zeta_k
      common /inteqn/ diag, ds, arcl
      common /sys_size/ k, nd, nbk
      common /fasblk2/ schur,wb,ipvtbf
c
      dimension nd(kmax), ds(kmax), arcl(kmax)
c
      dimension x_zeta(nmax+ng_max), y_zeta(nmax+ng_max),   
     1          diag(nmax)
      complex*16 zeta(nmax), dzeta(nmax), zeta_k(kmax)
c
c Fast Multipole Arrays
c
      complex*16 qa(nmax), cfield(nmax)
      dimension poten(nmax), wksp(nsp)
c
c  local work arrays
c
      real*4 timep(2), etime
      dimension A_k(kmax)
      complex*16 eye
c
         eye = dcmplx(0.d0,1.d0)
c
         t0 = etime(timep)
c
         call PRINI (6,13)
ccc         call prin2 ('zetak in fasmvp = *', zeta_k, 2*k)
c
         call FASMVP (k, nd, nbk, nsp, ds, x_zeta, y_zeta, zeta, dzeta,   
     1                zeta_k, diag, A_k, xx, yy, qa, cfield, poten, 
     2                wksp)
         t1 = etime(timep)
         tsec = t1 - t0
ccc         WRITE(13,*) 'TIME IN SECONDS FOR MATVEC = ',TSEC
ccc         WRITE(6,*) 'TIME IN SECONDS FOR MATVEC = ',TSEC
c
      RETURN
      END
c
c*********************
c
      subroutine DUMP (nx,ny,ugrid,igrid,ireal,if)
      implicit real*8 (a-h,o-z)
      dimension ugrid(nx,ny), igrid(nx,ny)
c
         DO i = 1,NX
            do j = 1, NY
               if (ireal.eq.1) then 
                  write(if,'(e20.13,$)')(ugrid(I,J))
                  write (if,'(a)')  ''
                 else
                  write(if,'(i4,$)') (igrid(i,j))
                  write (if,'(a)')  ''
               end if
            end do
         ENDDO
c
      return
      end
c
c*********************
c
      subroutine DUMP_MOVIE (nx, ny, time, ugrid, iframe, if)
      implicit real*8 (a-h,o-z)
      dimension ugrid(nx,ny), igrid(nx,ny)
c
         write (if,*) 'NX = ',ny,';'
         write (if,*) 'NY = ',nx+1,';'
         write (if,*) 'a = zeros(NX,NY);'
         WRITE(if,*) 'sol = ['
         DO i = 1,NX
            do j = 1, NY
               write(if,'(e20.13,$)')(ugrid(I,J))
               write (if,'(a)')  ''
            end do
         ENDDO
c
c periodic
         do j = 1, NY
               write(if,'(e20.13,$)')(ugrid(1,J))
               write (if,'(a)')  ''
         end do
         write (if,*) '];'
         write (if,*) 'a(:) = sol(:);'
         write (if,1200) time
         write (if,*) 'figure(1); clf;'
         write (if,*) '   surf(xgrid,ygrid,zgrid,a)'
         write (if,*) '   view([64,-4])'
         write (if,*) '   shading interp'
         write (if,*) '   colormap(jet2)'
         write (if,*) '   lighting phong'
         write (if,*) '   material dull'
         write (if,*) '   camlight(''headlight'')'
         write (if,*) '   caxis([-3 5])'
         write (if,*) '   hold on'
c         write (if,*) '   geo_3d'
         write (if,*) '   axis([-1 1 -1 1 -1 1])'
         write (if,*) '   axis equal'
         write (if,*) '   axis off'
         write (if,*) '   title([ ''t = '', num2str(time)]);'
         write (if,*) 
     1          "fname = sprintf('pngfiles/frame%.4d',", iframe,");"
         write (if,*) "print('-dpng', '-r150', fname);"
c
1200  FORMAT ('time = ',F6.2)
c
c
      return
      end
c
c*********************
c
      subroutine DUMP_MOVIE_ALL (nx, ny, time, ugrid, iframe, if)
      implicit real*8 (a-h,o-z)
      dimension ugrid(nx,ny), igrid(nx,ny)
c
         write (if,*) 'NX = ',ny,';'
         write (if,*) 'NY = ',nx+1,';'
         write (if,*) 'a = zeros(NX,NY);'
         WRITE(if,*) 'sol = ['
         DO i = 1,NX
            do j = 1, NY
               write(if,'(e20.13,$)')(ugrid(I,J))
               write (if,'(a)')  ''
            end do
         ENDDO
c
c periodic
         do j = 1, NY
               write(if,'(e20.13,$)')(ugrid(1,J))
               write (if,'(a)')  ''
         end do
         write (if,*) '];'
         write (if,*) 'a(:) = sol(:);'
         write (if,1200) time
         write (if,*) 'figure(1); clf;'
         write (if,*) 'subplot(2,2,1)'
         write (if,*) '   surf(xgrid,ygrid,zgrid,a)'
         write (if,*) '   view([64,-4])'
         write (if,*) '   shading interp'
         write (if,*) '   colormap(jet2)'
         write (if,*) '   lighting phong'
         write (if,*) '   material dull'
         write (if,*) '   camlight(''headlight'')'
         write (if,*) '   caxis([-3 5])'
         write (if,*) '   hold on'
c         write (if,*) '   geo_3d'
         write (if,*) '   axis([-1 1 -1 1 -1 1])'
         write (if,*) '   axis equal'
         write (if,*) '   axis off'
         write (if,*) '   title([ ''t = '', num2str(time)]);'
         write (if,*) 'subplot(2,2,2)'
         write (if,*) '   surf(xgrid,ygrid,zgrid,a)'
         write (if,*) '   view([-116,-4])'
         write (if,*) '   shading interp'
         write (if,*) '   colormap(jet2)'
         write (if,*) '   lighting phong'
         write (if,*) '   material dull'
         write (if,*) '   camlight(''headlight'')'
         write (if,*) '   caxis([-3 5])'
         write (if,*) '   hold on'
c         write (if,*) '   geo_3d'
         write (if,*) '   axis([-1 1 -1 1 -1 1])'
         write (if,*) '   axis equal'
         write (if,*) '   axis off'
c         write (if,*) '   title([ ''t = '', num2str(time)]);'
         write (if,*) 'subplot(2,2,3)'
         write (if,*) '   v = [0:.1:6.2];'
         write (if,*) '   contour(xzeta_grid,yzeta_grid,a,v)'
         write (if,*) '   hold on'
         write (if,*) '   blob'
         write (if,*) '   axis([-2.5 2.5 -2.5 2.5])'
         write (if,*) '   axis equal'
         write (if,*) '   axis off'
c         write (if,*) '   title([ ''t = '', num2str(time)]);'
         write (if,*) 
     1          "fname = sprintf('pngfiles/frame%.4d',", iframe,");"
         write (if,*) "print('-dpng', '-r150', fname);"
c
1200  FORMAT ('time = ',F6.2)
c
c
      return
      end
c
c*********************
c
      subroutine DUMP_MOVIE_VORT (nx, ny, time, xi, ugrid, iframe, if)
      implicit real*8 (a-h,o-z)
      dimension ugrid(nx,ny), igrid(nx,ny)
      complex*16 xi, eye
c
         eye = dcmplx(0.d0,1.d0)
c
c compute view point
         factor = 1.d0+(cdabs(xi))**2
         x1 = (xi+dconjg(xi))/factor
         x2 = (xi-dconjg(xi))/(eye*factor)
         x3 = -(1.d0-(cdabs(xi))**2)/factor
c
         write (if,*) 'NX = ',ny,';'
         write (if,*) 'NY = ',nx+1,';'
         write (if,*) 'a = zeros(NX,NY);'
         WRITE(if,*) 'sol = ['
         DO i = 1,NX
            do j = 1, NY
               write(if,'(e20.13,$)')(ugrid(I,J))
               write (if,'(a)')  ''
            end do
         ENDDO
c
c periodic
         do j = 1, NY
               write(if,'(e20.13,$)')(ugrid(1,J))
               write (if,'(a)')  ''
         end do
         write (if,*) '];'
         write (if,*) 'a(:) = sol(:);'
         write (if,1200) time
         write (if,*) 'figure(1); clf;'
         write (if,*) 'subplot(1,2,1)'
         write (if,*) '   surf(xgrid,ygrid,zgrid,a)'
         write (if,*) '   shading interp'
         write (if,*) '   colormap(jet2)'
         write (if,*) '   caxis([-3 5])'
         write (if,*) '   hold on'
c         write (if,*) '   geo_3d'
         write (if,*) '   axis([-1 1 -1 1 -1 1])'
         write (if,*) '   axis equal'
         write (if,*) '   view([',x1,',',x2,',',x3,'])'
         write (if,*) '   lighting phong'
         write (if,*) '   material dull'
         write (if,*) '   camlight(''headlight'')'
         write (if,*) '   axis off'
         write (if,*) '   title([ ''t = '', num2str(time)]);'
         write (if,*) 'subplot(1,2,2)'
         write (if,*) '   v = [0:.1:6.2];'
         write (if,*) '   contour(xzeta_grid,yzeta_grid,a,v)'
         write (if,*) '   hold on'
         write (if,*) '   blob'
         write (if,*) '   axis([-2.5 2.5 -2.5 2.5])'
         write (if,*) '   axis equal'
         write (if,*) '   axis off'
c         write (if,*) '   title([ ''t = '', num2str(time)]);'
         write (if,*) 
     1          "fname = sprintf('pngfiles/frame%.4d',", iframe,");"
         write (if,*) "print('-dpng', '-r150', fname);"
c
1200  FORMAT ('time = ',F6.2)
c
c
      return
      end
c
C**********************************************************************
C
	SUBROUTINE CLOSE(EPS,I1,I2,Q1,Q2,X1,Y1,X2,Y2,
     *                   FI1,FI2,POT1,POT2)
C
C   *** DESCRIPTION :
C
C       Close approch subroutine provided by the user.
C   Called when two particles are closer to each other than EPS.
C
C   *** INPUT PARAMETERS :
C
C   EPS     =  close approach distance
C   I1      =  number of first  particle
C   I2      =  number of second particle
C   Q1      =  charge of first  particle
C   Q2      =  charge of second particle
C   X1      =  x-component of position of first  particle
C   Y1      =  y-component of position of first  particle
C   X2      =  x-component of position of second particle
C   Y2      =  y-component of position of second particle
C
C   *** OUTPUT PARAMETERS :
C
C   FI1     =  field contribution from first particle on second
C   FI2     =  field contribution from second particle on first
C   POT1    =  potential contribution from first particle on second
C   POT2    =  potential contribution from second particle on first
C
C
C**********************************************************************
C
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	INTEGER II(2)
	DOUBLE PRECISION POS1(2),POS2(2)
	DOUBLE COMPLEX FI1,FI2,OUT2,Q1,Q2
	DATA DONE/1.0/,DHALF/0.5/
	II(1) = I1
	II(2) = I2
	POS1(1) = X1
	POS1(2) = Y1
	POS2(1) = X2
	POS2(2) = Y2
	FI1 = 0
	FI2 = 0
CCC     POT1 = 0
CCC     POT2 = 0
	POT1=DLOG( (X2-X1)**2 + (Y2-Y1)**2 ) /2
	POT2=POT1
ccc	   CALL PRINF('CLOSE PARTICLES ARE:*',II,2)
ccc	   CALL PRIN2(' COORDINATES OF FIRST  PART=*',POS1,2)
ccc	   CALL PRIN2(' COORDINATES OF SECOND PART=*',POS2,2)
	RETURN
	END

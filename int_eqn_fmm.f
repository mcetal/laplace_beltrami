      PROGRAM INT_EQN_FMM
c     ------------------------------------------------------------------
c
c
c  Solves integral equation for Laplace-Beltrami on the Sphere by 
c  using a stereographic projection and the FMM in the complex plane.
c
c  this solves the Dirichlet bvp in multiply connected domains
c  (each hole has a log source)
c
c  can build a grid of holes on surface of sphere
c
c
c     ------------------------------------------------------------------
      implicit real*8 (a-h, o-z)
      parameter (kmax = 500, npmax = 2048, nmax = kmax*npmax)
c
c Geometry of holes
      dimension ak(kmax), bk(kmax), th_k(kmax), phi_k(kmax), cx(kmax),
     1          cy(kmax), cz(kmax)
      dimension xs(nmax), ys(nmax), zs(nmax), dx(nmax), dy(nmax),
     1          dz(nmax), xn(nmax), yn(nmax), zn(nmax)
      dimension dsda(nmax), diag(nmax)
      dimension diag_stereo(nmax), d2x(nmax), d2y(nmax), d2z(nmax)
      complex*16 dzeta(nmax), zeta(nmax)
c
c  Grid variables
      parameter (nth_max = 1000, nphi_max = 1000, 
     1          ng_max = nth_max*nphi_max)
      dimension igrid(ng_max), th_gr(ng_max), phi_gr(ng_max),
     1          u_gr(ng_max), x_gr(ng_max), y_gr(ng_max), z_gr(ng_max),
     2          xzeta_gr(ng_max), yzeta_gr(ng_max), uex_gr(ng_max),
     3          alph_gr(ng_max)
      complex*16 zeta_gr(ng_max)
c
c target points are used to check accuracy
      dimension xtar(ng_max), ytar(ng_max), ztar(ng_max), 
     1          xz_tar(ng_max), yz_tar(ng_max), u_tar(ng_max)
      complex*16 zeta_tar(ng_max)    
cccc
cccc system matrices
cccc uncomment this only for smaller problems; memory problems otherwise
ccc      dimension asphere((nmax+kmax)**2), astereo((nmax+kmax)**2)      
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
      common /geometry/ x_zeta, y_zeta, zeta, dzeta, dsda
      common /inteqn/ diag_stereo, cx, cy, cz, zeta_k
      common /sys_size/ k, nd, nbk
      common /fasblk2/ schur,wb,ipvtbf
      common /sphere_int/ xs, ys, zs, xn, yn, zn, diag
c
c Open output files
c
c Initial Hole Geometry is given by reading in data
c   if iflag = 1, read in from input.data
c   if iflag = 2, read in from input_gridhole.data 
c   if iflag = 3, construct grid of holes
         iflag = 2
         if ((iflag.eq.1).or.(iflag.eq.2)) then  
            call READ_DATA (k, nd, nbk, nth, nphi, ak, bk, cx, cy, cz, 
     1                      th_k, phi_k, nvort, x1_vort, x2_vort, 
     2                      x3_vort, vort_k, gamma_tot, r0, zeta_k, 
     3                      iflag)
           else
c
c initial hole geometry is given by distributing random holes about 
c sphere
            call INIT_HOLE_GRID (k, nd, nbk, nth, nphi, ak, bk, cx, cy,  
     1                      cz, th_k, phi_k, nvort, x1_vort, x2_vort, 
     2                      x3_vort, vort_k, gamma_tot, r0, zeta_k)
         end if
         call DCFFTI (nd, wsave)
c
c Construct hole geometry and grid on surface of sphere
         call MAKE_GEO (k, nd, nbk, ak, bk, th_k, phi_k, xs, ys, zs,
     1                  dx, dy, dz, d2x, d2y, d2z, xn, yn, zn, dsda, 
     2                  diag)
c
c Get stereo graphic projection
         call STEREO (k, nd, nbk, xs, ys, zs, dx, dy, dz, d2x, d2y, d2z,
     1                zeta, dzeta, x_zeta, y_zeta, diag_stereo, 
     2                nvort, x1_vort, x2_vort, x3_vort, zk_vort) 
       stop
c
c Construct grid on surface of sphere
ccc         call SURFACE_GRID (k, nd, nbk, nth, nphi, ak, bk, cx, cy, cz, 
ccc     1                      th_k, phi_k, th_gr, phi_gr, x_gr, y_gr,
ccc     2                      z_gr, zeta_gr, xzeta_gr, yzeta_gr, igrid,
ccc     3                      alph_gr, xtar, ytar, ztar, ntar)
         call TARGET_POINTS (k, nd, nbk, ak, bk, cx, cy, cz, 
     1                      th_k, phi_k, xtar, ytar, ztar, ntar,
     2                      xz_tar, yz_tar, zeta_tar)
c
c Construct the RHS and solve
         call GETRHS (k, nd, nbk, cx, cy, cz, zeta_k, zeta, rhs,
     1                nvort, vort_k, zk_vort, gamma_tot)
ccc         call PRIN2 (' rhs = *', rhs, nbk)
         call SOLVE (nd, k, kmax, nbk, rhs, soln, density, A_k, gmwork, 
     1               lrwork, igwork, liwork, dsda, maxl, schur, wb,
     2               ipvtbf, zeta, zeta_k)
cccc
cccc Construct system matrices to be dumped into matlab
ccc         call BUILD_MAT_SPHERE (k, nd, nbk, xs, ys, zs, xn, yn, zn,
ccc     1                          dsda, diag, cx, cy, cz, asphere)
ccc         call BUILD_MAT_STEREO (k, nd, nbk, x_zeta, y_zeta, zeta,     
ccc     1                          dzeta, zeta_k, diag_stereo, astereo)
cccc
cccc Overresolve the data on the boundary (for plotting purposes only)
ccc         call RESAMPLE (nd,k,nbk,nsamp,h,z,dz,xat,yat,u,xsamp,
ccc     *                  ysamp,zwork,ncwork)
c
c Construct solution on surface grid
ccc         call SOL_GRID (nd, k, nbk, nth, nphi, density, A_k, zeta_k,   
ccc     1                  zeta, dzeta, igrid, zeta_gr, u_gr)
ccc         call SOL_GRID_FMM (nd, k, nbk, nth, nphi, density, A_k, zeta_k,   
ccc     1                      zeta, dzeta, igrid, zeta_gr, u_gr,
ccc     2                      x_zeta, y_zeta, qa, cfield, poten, nsp, 
ccc     3                      wksp, nvort, vort_k, zk_vort)
         call SOL_TAR_FMM (nd, k, nbk, ntar, density, A_k, zeta_k,   
     1                     x_zeta, y_zeta, zeta, dzeta, zeta_tar, u_tar,
     2                     xz_tar, yz_tar, qa, cfield, poten, nsp, 
     3                     wksp, nvort, vort_k, zk_vort)
c
c calculate 
c
c Check error in solution
ccc         call CHECK_ERROR (nd, k, nbk, nth, nphi, zeta_k, igrid, 
ccc     1                     zeta_gr, u_gr)
         call CHECK_ERROR_TAR (nd, k, nbk, ntar, zeta_k, zeta_tar,  
     1                         u_tar)
c
      stop
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine READ_DATA (k, nd, nbk, nth, nphi, ak, bk, cx, cy, cz, 
     1                      th_k, phi_k, nvort, x1_vort, x2_vort, 
     2                      x3_vort, vort_k, gamma_tot, r0, zeta_k,
     3                      iflag)
c---------------
      implicit real*8 (a-h,o-z)
      dimension ak(*), bk(*), cx(*), cy(*), cz(*), th_k(*), phi_k(*)
      dimension x1_vort(*), x2_vort(*), x3_vort(*), vort_k(*)
      complex*16 zeta_k(*), eye
c
         eye = dcmplx(0.d0,1.d0)
c
         if (iflag.eq.1) then
            open (unit = 12, file = 'input.data')
           else
            open (unit = 12, file = 'input_gridhole.data')
         end if
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
c
c a fudge for pole
ccc         call SPH2CART (0.d0, 1.4d0, 1.d0, cx(1), cy(1), cz(1))
         eps = 1.d-6
         do kbod = 1, k
            check = dabs(cz(kbod)-1.d0)
c
c be careful if one of the hole centres is at the north pole
c if it is, nudge it a little
            if (check.lt.eps) then
               cz_s = 0.9999d0
               cx_s = dsqrt(0.5d0*(1-cz_s**2))
               cy_s = cx_s
               zeta_k(kbod) = (cx_s + eye*cy_s)/(1.d0-cz_s)
               write (6,*) 'Fudging hole centre a little'
               call prin2 (' old centre, cx = *', cx(kbod), 1)
               call prin2 (' old centre, cy = *', cy(kbod), 1)
               call prin2 (' old centre, cz = *', cz(kbod), 1)
               call prin2 (' new centre, cx = *', cx_s, 1)
               call prin2 (' new centre, cy = *', cy_s, 1)
               call prin2 (' new centre, cz = *', cz_s, 1)
            else
               zeta_k(kbod) = (cx(kbod) + eye*cy(kbod))/(1.d0-cz(kbod))
            end if
         end do
ccc         call PRIN2 (' zeta_k = *', zeta_k, 2*k)
         close(12)
c
      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine INIT_HOLE_GRID (k, nd, nbk, nth, nphi, ak, bk, cx, cy,  
     1                      cz, th_k, phi_k, nvort, x1_vort, x2_vort, 
     2                      x3_vort, vort_k, gamma_tot, r0, zeta_k)
c---------------
      implicit real*8 (a-h,o-z)
      dimension ak(*), bk(*), cx(*), cy(*), cz(*), th_k(*), phi_k(*)
      dimension x1_vort(*), x2_vort(*), x3_vort(*), vort_k(*)
      dimension arand(1000), brand(1000)
      complex*16 zeta_k(*), eye
c
         pi = 4.d0*datan(1.d0)
         eye = dcmplx(0.d0,1.d0)
c
         open (unit = 12, file = 'input.data')
         call PRINI (6,13) 
c
c read in hole data from file to get nd and vortices (we will ignore 
c k)
         read (12,*) k, nd, nvort
         read(12,*) nth, nphi
         do kbod = 1, k
            read(12,*) ak(kbod), bk(kbod), th_k(kbod), phi_k(kbod)
         end do
         call PRINF (' nvort = *', nvort, 1)
         gamma_tot = 0.d0
         do ivort = 1, nvort
            read(12,*) theta, phi, vort_k(ivort)
            call SPH2CART (theta, phi, 1.d0, x1_vort(ivort),
     1                     x2_vort(ivort), x3_vort(ivort))
            gamma_tot = gamma_tot + vort_k(ivort)
         end do
c
         close(12)
c
         call PRINF (' nth = *', nth, 1)
         call PRINF (' nphi = *', nphi, 1)
         call PRIN2 ('    vort_k = *', vort_k, nvort)
         call PRIN2 ('    x1_vort = *', x1_vort, nvort)
         call PRIN2 ('    x2_vort = *', x2_vort, nvort)
         call PRIN2 ('    x3_vort = *', x3_vort, nvort)
c
         np = 17
         nt = 27
ccc         np = 5
ccc         nt = 5
         k = (np-2)*nt + 2
         call PRINF (' k = *', k, 1)
         nbk = k*nd
         call PRINF (' nd = *', nd, 1)
         call PRINF (' NBK = *', nbk, 1)
         dphi = pi/(np-1)
         dth = 2*pi/nt
c
c Generate random numbers
         seed = 0
         call zufalli(seed)
         call zufallt(k,arand)
         call zufallt(k,brand)
         call PRIN2 (' arand = *', arand, k)
         call PRIN2 (' brand = *', brand, k)
c
c hole at north pole
         rmax = 0.5d0*dphi
ccc         rmax = dphi
         ak(1) = (0.01d0 + 0.98d0*arand(1))*rmax
         bk(1) = -(0.01d0+0.98d0*brand(1))*rmax
ccc         ak(1) = rmax
ccc         bk(1) = -rmax
         phi_k(1) = 0.5d0*pi
         th_k(1) = 0.d0
         call SPH2CART (th_k(1), phi_k(1), 1.d0, cx(1), cy(1), cz(1))
c
c grid holes between the poles
         kbod = 2
         do ip = 2, np-1
            call prinf('ip = *', ip, 1) 
            phi = 0.5d0*pi-dphi*(ip-1.d0)
            call prin2 (' phi = *', phi, 1)
            call prin2 (' dcos(phi) = *', dcos(phi), 1)
            rmax = 0.5d0*min(dphi,dabs(dcos(phi))*dth)
ccc            rmax = min(dphi,dabs(dcos(phi))*dth)
            call PRIN2 ('  rmax = *', rmax, 1)
            do it = 1, nt
               theta = dth*(it-1.d0)
               phi_k(kbod) = phi
               th_k(kbod) = theta
               ak(kbod) = (0.01d0 + 0.98d0*arand(kbod))*rmax
               bk(kbod) = -(0.01d0 + 0.98d0*brand(kbod))*rmax
ccc               ak(kbod) = rmax
ccc               bk(kbod) = -rmax
               call SPH2CART (th_k(kbod), phi_k(kbod), 1.d0, cx(kbod),
     1                        cy(kbod), cz(kbod))
               kbod = kbod+1
            end do
         end do
         call prinf (' kbod should be k = *', kbod, 1)
c
c hole at South pole
         rmax = 0.5d0*dphi
         ak(k) = (0.01d0 + 0.98d0*arand(k))*rmax
         bk(k) = -(0.01d0 + 0.98d0*brand(k))*rmax
ccc         ak(k) = rmax
ccc         bk(k) = -rmax
         phi_k(k) = -0.5d0*pi
         th_k(k) = 0.d0
         call SPH2CART (th_k(k), phi_k(k), 1.d0, cx(k), cy(k), cz(k))
         call prin2 (' ak = *', ak, k)
         call prin2 (' bk = *', bk, k)
c
c a fudge for pole
ccc         call SPH2CART (0.d0, 1.4d0, 1.d0, cx(1), cy(1), cz(1))
         eps = 1.d-6
         do kbod = 1, k
            check = dabs(cz(kbod)-1.d0)
c
c be careful if one of the hole centres is at the north pole
c if it is, nudge it a little
            if (check.lt.eps) then
               cz_s = 0.9999d0
               cx_s = dsqrt(0.5d0*(1-cz_s**2))
               cy_s = cx_s
               zeta_k(kbod) = (cx_s + eye*cy_s)/(1.d0-cz_s)
               write (6,*) 'Fudging hole centre a little'
               call prin2 (' old centre, cx = *', cx(kbod), 1)
               call prin2 (' old centre, cy = *', cy(kbod), 1)
               call prin2 (' old centre, cz = *', cz(kbod), 1)
               call prin2 (' new centre, cx = *', cx_s, 1)
               call prin2 (' new centre, cy = *', cy_s, 1)
               call prin2 (' new centre, cz = *', cz_s, 1)
               call prin2 (' new zeta_k = *', zeta_k(kbod), 2)
            else
               zeta_k(kbod) = (cx(kbod) + eye*cy(kbod))/(1.d0-cz(kbod))
            end if
         end do
ccc         call PRIN2 (' zeta_k = *', zeta_k, 2*k)
c
c write out everything to be read in the next time
         open (unit=12, file = 'input_gridhole.data')
         write (12,*) k, nd, nvort
         write (12,*) nth, nphi
         do kbod = 1, k
            write (12,'(4(e24.15,$))') ak(kbod), bk(kbod), th_k(kbod), 
     1                 phi_k(kbod)
            write (12,'(a)')  ''
         end do
         do ivort = 1, nvort
            read(12,'(4(e24.15,$))') theta, phi, vort_k(ivort)
            write (12,'(a)')  ''
         end do
         close (12)
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
         open (unit = 42, file = 'geo_3d.m')
         do kbod = 1, k
            call RS_3D_PLOT (xs(is),ys(is),zs(is), nd, 1, 42)
            is = is + nd
         end do
         close (42)
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
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine SURFACE_GRID (k, nd, nbk, nth, nphi, ak, bk, cx, cy,  
     1                         cz, th_k, phi_k, th_gr, phi_gr, x_gr, 
     2                         y_gr, z_gr, zeta_gr, xzeta_gr, yzeta_gr,
     3                         igrid, alph_gr, xtar, ytar, ztar, ntar)
c---------------
      implicit real*8 (a-h,o-z)
      dimension ak(k), bk(k), cx(k), cy(k), cz(k), th_k(k), phi_k(k)
      dimension igrid(nth,nphi), th_gr(nth,nphi), phi_gr(nth,nphi),
     1          x_gr(nth,nphi), y_gr(nth,nphi), z_gr(nth,nphi),
     2          xzeta_gr(nth,nphi), yzeta_gr(nth,nphi), 
     3          alph_gr(nth,nphi)
      dimension xtar(*), ytar(*), ztar(*)
      complex*16 zeta_gr(nth,nphi), eye
c
         pi = 4.d0*datan(1.d0)
         eye = dcmplx(0.d0,1.d0)
c
c Calculate epsilon
         radmax = 0.d0
         do kbod = 1, k
            radmax = max(radmax, dabs(ak(kbod)))
            radmax = max(radmax, dabs(bk(kbod)))
         end do
         call PRIN2 (' maximum radius = *', radmax, 1)
         eps = 2*2.d0*pi*radmax/nd
ccc         eps = 0.05d0
         call PRIN2 (' Epsilon = *', eps, 1) 
c
         dth = 2.d0*pi/nth
         dphi = pi/nphi
         ntar = 0
         do i = 1, nth
            theta = (i-1)*dth
            do j = 1, nphi
               phi = (j-0.5d0)*dphi - 0.5d0*pi
               th_gr(i,j) = theta
               phi_gr(i,j) = phi
               call SPH2CART (th_gr(i,j), phi_gr(i,j), 1.d0, 
     1                        x_gr(i,j), y_gr(i,j), z_gr(i,j))
               zeta_gr(i,j) = (x_gr(i,j) + eye*y_gr(i,j))/
     1                             (1.d0 - z_gr(i,j))
               xzeta_gr(i,j) = dreal(zeta_gr(i,j))
               yzeta_gr(i,j) = dimag(zeta_gr(i,j))
               in_out = 0
               in_out2 = 0
               do kbod = 1, k
                  call IN_OR_OUT (x_gr(i,j), y_gr(i,j), z_gr(i,j), eps,  
     1                            ak(kbod), bk(kbod), th_k(kbod), 
     2                            phi_k(kbod), itest)
                  call IN_OR_OUT (x_gr(i,j), y_gr(i,j), z_gr(i,j), 0.d0,  
     1                            ak(kbod), bk(kbod), th_k(kbod), 
     2                            phi_k(kbod), itest2)
                  in_out = in_out + itest
                  in_out2 = in_out2 + itest2
               end do
               if (in_out>0) then
                  igrid(i,j) = 0
                 else
                  igrid(i,j) = 1
               end if
               if (in_out2>0) then
                  alph_gr(i,j) = 0.9
                 else
                  alph_gr(i,j) = 0.1
               end if
            end do
         end do 
         open (unit = 31, file = 'igrid.dat')
         open (unit = 32, file = 'xgrid.dat')
         open (unit = 33, file = 'ygrid.dat')
         open (unit = 34, file = 'zgrid.dat')
         open (unit = 35, file = 'xzeta_grid.dat')
         open (unit = 36, file = 'yzeta_grid.dat')
         call DUMP (nth, nphi, x_gr, igrid, 0, 31)
         call DUMP (nth, nphi, x_gr, igrid, 1, 32)
         call DUMP (nth, nphi, y_gr, igrid, 1, 33)
         call DUMP (nth, nphi, z_gr, igrid, 1, 34)
         call DUMP (nth, nphi, xzeta_gr, igrid, 1, 35)
         call DUMP (nth, nphi, yzeta_gr, igrid, 1, 36)
         close (31)
         close (32)
         close (33)
         close (34)
         close (35)
         close (36)
c
      return
      end      
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine TARGET_POINTS (k, nd, nbk, ak, bk, cx, cy, cz, 
     1                          th_k, phi_k, xtar, ytar, ztar, ntar,
     2                          xz_tar, yz_tar, zeta_tar)
c---------------
      implicit real*8 (a-h,o-z)
      dimension ak(k), bk(k), cx(k), cy(k), cz(k), th_k(k), phi_k(k)
      dimension xtar(*), ytar(*), ztar(*), xz_tar(*), yz_tar(*)
      complex*16 zeta, eye, zeta_tar(*)
c
         pi = 4.d0*datan(1.d0)
         eye = dcmplx(0.d0,1.d0)
c
c Calculate epsilon
         radmax = 0.d0
         do kbod = 1, k
            radmax = max(radmax, dabs(ak(kbod)))
            radmax = max(radmax, dabs(bk(kbod)))
         end do
         call PRIN2 (' maximum radius = *', radmax, 1)
ccc         eps = 20*2.d0*pi*radmax/nd
         eps = 1.0d0
         call PRIN2 (' Epsilon = *', eps, 1) 
c
         nth = 10
         nphi = 15
         dth = 2.d0*pi/nth
         dphi = pi/nphi
         ntar = 0
         do i = 1, nth
            theta = (i-1)*dth
            do j = 1, nphi
               phi = (j-0.5d0)*dphi - 0.5d0*pi
               call SPH2CART (theta, phi, 1.d0, x, y, z)
               zeta = (x + eye*y)/(1.d0 - z)
               xzeta = dreal(zeta)
               yzeta = dimag(zeta)
               in_out = 0
               do kbod = 1, k
                  call IN_OR_OUT (x, y, z, eps,  
     1                            ak(kbod), bk(kbod), th_k(kbod), 
     2                            phi_k(kbod), itest)
                  in_out = in_out + itest
               end do
               if (in_out.eq.0) then
                  ntar = ntar+1
                  xtar(ntar) = x
                  ytar(ntar) = y
                  ztar(ntar) = z
                  zeta_tar(ntar) = zeta
                  xz_tar(ntar) = xzeta
                  yz_tar(ntar) = yzeta
               end if
            end do
         end do 
c
c blah
ccc         ntar = 20
ccc         dth = 2.d0*pi/ntar
ccc         do i = 1, ntar
ccc            theta = dth*(i-1)
ccc            zeta_tar(i) = 0.5d0*cdexp(eye*theta)
ccc            xz_tar(i) = dreal(zeta_tar(i))
ccc            yz_tar(i) = dimag(zeta_tar(i))
ccc         end do
c
         open (unit = 52, file = 'targets_3d.m')
         open (unit = 53, file = 'stereo_targets.m')
         call RSCPLOT (zeta_tar, ntar, 1, 53)
         call RS_3D_PLOT (xtar, ytar, ztar, ntar, 1, 52)
         close (32)
         close (33)
c
      return
      end      
c
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine STEREO (k, nd, nbk, xs, ys, zs, dx, dy, dz, d2x, d2y,  
     1                   d2z, zeta, dzeta, x_zeta,y_zeta, diag, 
     2                   nvort, x1_vort, x2_vort, x3_vort, zk_vort)
c---------------
      implicit real*8 (a-h,o-z)
      dimension xs(nbk), ys(nbk), zs(nbk), dx(nbk), dy(nbk), dz(nbk),
     1          x_zeta(nbk), y_zeta(nbk), d2x(nbk), d2y(nbk), d2z(nbk),
     2          diag(nbk), x1_vort(nvort), x2_vort(nvort), 
     3          x3_vort(nvort)
      complex*16 zeta(nbk), dzeta(nbk), eye, d2zeta, zextra, 
     1           zk_vort(nvort)
c
         pi = 4.d0*datan(1.d0)
         eye = dcmplx(0.d0,1.d0)
c
         do i = 1, nbk
            zeta(i) = (xs(i) + eye*ys(i))/(1.d0-zs(i))
            x_zeta(i) = dreal(zeta(i))
            y_zeta(i) = dimag(zeta(i))
            dzeta(i) = (dx(i) + eye*dy(i) + zeta(i)*dz(i))/(1.d0-zs(i))
            xdot = dreal(dzeta(i))
            ydot = dimag(dzeta(i))
            ds = cdabs(dzeta(i))
            d2zeta = (d2x(i) + eye*d2y(i) + 2.d0*dzeta(i)*dz(i) + 
     1                zeta(i)*d2z(i))/(1.d0-zs(i))
            xddot = dreal(d2zeta)
            yddot = dimag(d2zeta)
            rkappa = (xdot*yddot-ydot*xddot)/ds**3
ccc            call prin2 (' rkappa = *', rkappa, 1)
            zextra = dzeta(i)*dconjg(zeta(i))/(1.d0+cdabs(zeta(i))**2)
            zextra = zextra/(2.d0*pi)
            diag(i) = 0.25d0*rkappa*ds/pi - dimag(zextra)
ccc            diag(i) = 0.25d0*rkappa*ds/pi 
         end do
ccc         call prin2 (' zeta 1 = *', zeta, 2*nd)
ccc         call prin2 (' zeta 2 = *', zeta(1+nd), 2*nd)
ccc         call prin2 (' dzeta 1 = *', dzeta, 2*nd)
ccc         call prin2 (' dzeta 2 = *', dzeta(1+nd), 2*nd)
ccc         call prin2 (' diag = *', diag, nbk)
         open (unit = 22, file = 'stereo_vortices.m')
         open (unit = 23, file = 'vortices_3d.m')
         do ivort = 1, nvort
            zk_vort(ivort) = (x1_vort(ivort) + eye*x2_vort(ivort))/
     1                       (1.d0 - x3_vort(ivort))
         end do
         call RSCPLOT (zk_vort, nvort, 1, 22)
         call RS_3D_PLOT (x1_vort,x2_vort,x3_vort,nvort, 1, 23)
         close (22)
         close (23)
         open (unit = 11, file = 'stereo_geo.m')
         do kbod = 1, k
            call RSCPLOT (zeta((kbod-1)*nd+1), nd, 1, 11)
         end do
         close (11)
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
      subroutine GETRHS (k, nd, nbk, cx, cy, cz, zeta_k, zeta, rhs,
     1                   nvort, vort_k, zk_vort, gamma_tot)
c---------------
      implicit real*8 (a-h,o-z)
      dimension cx(k), cy(k), cz(k), rhs(nbk), xs(nbk), ys(nbk), 
     1          zs(nbk), vort_k(nvort)
      complex*16 zeta_k(k), eye, zeta(nbk), zk_vort(nvort), zdis
c
         eye = dcmplx(0.d0,1.d0)
         call PRIN2 (' zeta_k = *', zeta_k, 2*k)
c
         istart = 0
         do kbod = 1, k
            do j = 1, nd
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
            istart = istart+nd
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
      subroutine MATVEC_SPHERE (k, nd, nbk, xs, ys, zs, xn, yn, zn,
     1                          dsda, diag, cx, cy, cz, A_k, u, w)
c---------------
      implicit real*8 (a-h,o-z)
      dimension u(nbk), xs(nbk), ys(nbk), zs(nbk), xn(nbk), yn(nbk),
     1          zn(nbk), diag(nbk), w(nbk), dsda(nbk), cx(k), cy(k),
     2          cz(k), A_k(k)
      dimension gradu(3), vn(3), R0(3), R(3), delR(3)
c
         pi = 4.d0*datan(1.d0)
         dalph = 2.d0*pi/nd
c
c
c Extract off A_k
         do kbod = 1, k
            A_k(kbod) = u(nbk+kbod)
         end do
c
         do i = 1, nbk
            w(i) = 0.d0
            R0(1) = xs(i)
            R0(2) = ys(i)
            R0(3) = zs(i)
            do j = 1, nbk
               if (i.ne.j) then  
                  vn(1) = xn(j)
                  vn(2) = yn(j)
                  vn(3) = zn(j)
                  R(1) = xs(j)
                  R(2) = ys(j)
                  R(3) = zs(j)
                  do ik = 1, 3
                     delR(ik) = R(ik) - R0(ik)
                  end do
                  call DOT (delR, delR, rad2)
                  call DOT (delR, vn, du_dn)
                   du_dn = du_dn/rad2
                   du = dalph*u(j)*du_dn*dsda(j)/(2.d0*pi)
ccc               aKu(i) = aKu(i) - dalph*qa(j)*du_dn*dsda(j)/(2.d0*pi)
                  w(i) = w(i) +
     1                       dalph*u(j)*du_dn*dsda(j)/(2.d0*pi)
                else
                  w(i) = w(i) + 0.5d0*u(i) - dalph*u(i)*diag(i)
               end if
               do kbod = 1, k
                  R(1) = cx(kbod)
                  R(2) = cy(kbod)
                  R(3) = cz(kbod)
                  do ik = 1, 3
                     delR(ik) = R(ik) - R0(ik)
                  end do
                  call DOT (delR, delR, rad2)
                  w(i) = w(i) + A_k(kbod)*dlog(dsqrt(rad2/2.d0))
               end do
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
               w(nbk + kbod) = w(nbk + kbod) + u(i)*dsda(i)
            end do
         end do 
c
      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine BUILD_MAT_SPHERE (k, nd, nbk, xs, ys, zs, xn, yn, zn,
     1                             dsda, diag, cx, cy, cz, amat)
c---------------
      implicit real*8 (a-h,o-z)
      dimension xs(nbk), ys(nbk), zs(nbk), xn(nbk), yn(nbk),
     1          zn(nbk), diag(nbk), dsda(nbk), cx(k), cy(k),
     2          cz(k), amat(nbk+k,nbk+k)
      dimension gradu(3), vn(3), R0(3), R(3), delR(3)
c
         pi = 4.d0*datan(1.d0)
         dalph = 2.d0*pi/nd
c
         do i = 1, nbk
            R0(1) = xs(i)
            R0(2) = ys(i)
            R0(3) = zs(i)
            do j = 1, nbk
               if (i.ne.j) then  
                  vn(1) = xn(j)
                  vn(2) = yn(j)
                  vn(3) = zn(j)
                  R(1) = xs(j)
                  R(2) = ys(j)
                  R(3) = zs(j)
                  do ik = 1, 3
                     delR(ik) = R(ik) - R0(ik)
                  end do
                  call DOT (delR, delR, rad2)
                  call DOT (delR, vn, du_dn)
                   du_dn = du_dn/rad2
ccc               aKu(i) = aKu(i) - dalph*qa(j)*du_dn*dsda(j)/(2.d0*pi)
                  amat(i,j) = dalph*du_dn*dsda(j)/(2.d0*pi)
                else
                  amat(i,i) =  0.5d0 - dalph*diag(i)
               end if
               do kbod = 1, k
                  R(1) = cx(kbod)
                  R(2) = cy(kbod)
                  R(3) = cz(kbod)
                  do ik = 1, 3
                     delR(ik) = R(ik) - R0(ik)
                  end do
                  call DOT (delR, delR, rad2)
                  amat(i,nbk+kbod) = dlog(dsqrt(rad2/2.d0))
               end do
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
               amat(nbk+kbod,i) = dsda(i)
            end do
         end do 
c
c Dump it out
         open (unit = 24, file = 'amat_sphere.dat')
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
      subroutine MATVEC_STEREO (k, nd, nbk, qa, zeta, zs, dx, dy, dz,
     1                          dsda, diag, aKu)
c---------------
      implicit real*8 (a-h,o-z)
      dimension qa(nbk), dx(nbk), dy(nbk), dz(nbk), diag(nbk), 
     1          aKu(nbk), dsda(nbk), zs(nbk)
      dimension gradu(3), vn(3), R0(3), R(3), delR(3)
      complex*16 zeta(nbk), dzeta, zdis, eye, zkern
c
         pi = 4.d0*datan(1.d0)
         eye = dcmplx(0.d0,1.d0)
         dalph = 2.d0*pi/nd
c
         do i = 1, nbk
            aKu(i) = 0.d0
            do j = 1, i-1
               zdis = zeta(j) - zeta(i)
               dzeta = (dx(j) + eye*dy(j) + zeta(j)*dz(j))/(1.d0-zs(j))
ccc               zkern = dzeta/zdis 
               zkern = dzeta/zdis - dconjg(zeta(j))*dzeta
     1                   /(1.d0+cdabs(zeta(j))**2)
               aKu(i) = aKu(i) - dalph*qa(j)*dimag(zkern)/(2.d0*pi)
            end do
            do j = i+1, nbk
               zdis = zeta(j) - zeta(i)
               dzeta = (dx(j) + eye*dy(j) + zeta(j)*dz(j))/(1.d0-zs(j))
ccc               zkern = dzeta/zdis 
               zkern = dzeta/zdis - dconjg(zeta(j))*dzeta
     1                   /(1.d0+cdabs(zeta(j))**2)
               aKu(i) = aKu(i) - dalph*qa(j)*dimag(zkern)/(2.d0*pi)
            end do
         end do
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
      subroutine  FASMVP (k, nd, nbk, nsp, x_zeta, y_zeta, zeta,    
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
         napb = 40
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
            w(i) = 0.5d0*u(i) + dimag(cfield(i)) - dalph*diag(i)*u(i)
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
     1                  lrwork, iwork, liwork, dsda, maxl, schur, bnew,
     2                  ipvtbf, z, zk)
c---------------
c
      implicit real*8 (a-h,o-z)
      external MATVEC_LAPL, MSOLVE
c
c  System
      dimension soln(*), rhs(*), u(nbk), A_k(k), dsda(nbk)
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
         tol = 1.0d-11
         isym = 0
         iwork(1) = maxl
         do i=2,liwork
            iwork(i) = 0
         enddo
c
c  Preconditioner flag
c
         iwork(4) = -1
ccc         iwork(4) = 0
c
c  Restart flag
c  
         iwork(5) = 5      
c
c  factor preconditioner
         t0 = etime(timep)
         call SCHUR_FACTOR (z,ND,nbk,K,zk,SCHUR,bnew,IPVTBF)
         t1 = etime(timep)
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
            do kbod = 1, k
ccc               call PRINF (' kbod = *', kbod, 1)
ccc               call PRIn2 ('     u = *', u((kbod-1)*nd+1), nd)
               A_k(kbod) = soln(nbk+kbod)
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
      subroutine SOL_GRID_FMM (nd, k, nbk, nth, nphi, u, A_k, zeta_k,   
     1                         zeta, dzeta, igrid, zeta_gr, u_gr,
     2                         x_zeta, y_zeta, qa, cfield, poten, nsp, 
     3                         wksp, nvort, vort_k, zk_vort)
c---------------
c
      implicit real*8 (a-h,o-z)
      dimension u(nbk), igrid(nth,nphi), u_gr(nth,nphi), A_k(k),
     1          vort_k(nvort)
      complex*16 zeta(nbk), dzeta(nbk), zeta_gr(nth,nphi), zkern, 
     1           zeta_k(k), zdis, eye, qa(*), cfield(*), zQsum, 
     2           zQ2sum, zk_vort(nvort)
      dimension x_zeta(*), y_zeta(*), poten(*), wksp(nsp)
      integer*4 iout(2), inform(10), ierr(10)
      REAL*4 TIMEP(2), ETIME
c
         pi = 4.d0*datan(1.d0)
         eye = dcmplx(0.d0, 1.d0)
         dalph = 2.d0*pi/nd
c
c pack zeta and zeta_gr into x_zeta and y_zeta
         do i = 1, nbk
            x_zeta(i) = dreal(zeta(i))
            y_zeta(i) = dimag(zeta(i))
            qa(i) = dalph*u(i)*dzeta(i)/(2.d0*pi)
         end do
         ij = nbk
         do i = 1, nth
            do j = 1, nphi
               if (igrid(i,j).ne.0) then
                  ij = ij + 1
                  x_zeta(ij) = dreal(zeta_gr(i,j))
                  y_zeta(ij) = dimag(zeta_gr(i,j))
                  qa(ij) = 0.d0
               end if
            end do
         end do
         call PRINF (' ij = *', ij, 1)
c
         zQsum = 0.d0
         do i = 1, nbk
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
ccc         call PRIN2 (' a_k in sol_GRID_FMM = *', A_k, k)
c Fix up field
         ij = nbk
         umax = -1.d10
         umin = 1.d10
         do i = 1, nth
            do j = 1, nphi
               u_gr(i,j) = -10.d0
               if (igrid(i,j).ne.0) then     
                  ij = ij + 1           
                  u_gr(i,j) = dimag(cfield(ij) - zQsum)
                  do kbod = 1, k
ccc                  do kbod = 1, 1
                     zdis = zeta_gr(i,j) - zeta_k(kbod)
                     rad = 2.d0*
     1              (cdabs(zdis))**2/((1+(cdabs(zeta_gr(i,j)))**2)
     2                  *((1+(cdabs(zeta_k(kbod)))**2)))
                     u_gr(i,j) = u_gr(i,j) + A_k(kbod)*0.5d0*dlog(rad)
                  end do
                  umax = max(umax,u_gr(i,j))
                  umin = min(umin,u_gr(i,j))
                 else
                  u_gr(i,j) = -100000.d0
               end if
            end do
         end do
         call PRIN2 (' Max solution = *', umax, 1)
         call PRIN2 (' Min solution = *', umin, 1)
c
c Add in vortex singularities
ccc         call PRIN2 (' zk_vort = *', zk_vort, 2*nvort)
         call PRIN2 (' vort_k = *', vort_k, nvort)
         do i = 1, nth
            do j = 1, nphi
cccccc               u_gr(i,j) = 0.d0
               if (igrid(i,j).ne.0) then   
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
         open (unit=43, file = 'ugrid.dat')
         call DUMP (nth, nphi, u_gr, igrid, 1, 43)
         close (43)
c
      return
      end
c
c
c---------------
      subroutine SOL_TAR_FMM (nd, k, nbk, ntar, u, A_k, zeta_k,   
     1                        x_zeta, y_zeta, zeta, dzeta, zeta_tar, 
     2                        u_tar, xz_tar, yz_tar, qa, cfield, poten,  
     3                        nsp, wksp, nvort, vort_k, zk_vort)
c---------------
c
      implicit real*8 (a-h,o-z)
      dimension u(nbk), u_tar(ntar), A_k(k), vort_k(nvort), x_zeta(*),
     1          y_zeta(*)
      complex*16 zeta(nbk), dzeta(nbk), zeta_tar(ntar), zkern, 
     1           zeta_k(k), zdis, eye, qa(*), cfield(*), zQsum, 
     2           zQ2sum, zk_vort(nvort), ztar
      dimension xz_tar(ntar), yz_tar(ntar), poten(*), wksp(nsp)
      integer*4 iout(2), inform(10), ierr(10)
      REAL*4 TIMEP(2), ETIME
c
         pi = 4.d0*datan(1.d0)
         eye = dcmplx(0.d0, 1.d0)
         dalph = 2.d0*pi/nd
c
c pack zeta and zeta_gr into x_zeta and y_zeta
         do i = 1, nbk
            x_zeta(i) = dreal(zeta(i))
            y_zeta(i) = dimag(zeta(i))
            qa(i) = dalph*u(i)*dzeta(i)/(2.d0*pi)
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
         do i = 1, nbk
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
ccc         call prin2 (' u_tar = *', u_tar, ntar)
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
      subroutine CHECK_ERROR (nd, k, nbk, nth, nphi, zeta_k, igrid, 
     1                        zeta_gr, u_gr)
c---------------
c
      implicit real*8 (a-h,o-z)
      dimension igrid(nth,nphi), u_gr(nth,nphi)
      complex*16 zeta_gr(nth,nphi), zkern, zeta_k(k), zdis, eye
c
         pi = 4.d0*datan(1.d0)
         eye = dcmplx(0.d0, 1.d0)
         dalph = 2.d0*pi/nd
c
ccc         call PRIn2 (' zeta_k in check_ERROR = *', zeta_k, 2*k)
ccc         call PRIN2 (' zeta_gr = *', zeta_gr(1,1), 2)
         err = 0.d0
         do i = 1, nth
            do j = 1, nphi
               if (igrid(i,j).eq.1) then
               u_ex = 0.d0
               do mbod = 1, k
                  u_ex = u_ex + 1.d0/(zeta_gr(i,j)-zeta_k(mbod)) 
     1                   + dconjg(1.d0/(zeta_gr(i,j)-zeta_k(mbod)))
               end do
ccc               u_ex = u_ex + 66.d0 
               err = max(err,dabs(u_ex-u_gr(i,j)))
               call PRIN2 ('### u_ex = *', u_ex, 1)
               call PRIN2 ('    u_gr  = *', u_gr(i,j), 1)
               end if
            end do
         end do
         call PRIN2 (' max error in solution on grid = *', err, 1)
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
      dimension u_tar(ntar)
      complex*16 zeta_tar(ntar), zkern, zeta_k(k), zdis, eye
c
         pi = 4.d0*datan(1.d0)
         eye = dcmplx(0.d0, 1.d0)
         dalph = 2.d0*pi/nd
c
ccc         call PRIn2 (' zeta_k in check_ERROR = *', zeta_k, 2*k)
ccc         call PRIN2 (' zeta_tar = *', zeta_tar, 2)
ccc         zeta_k(1) = dcmplx(-10.d0,-10.d0)
ccc         call prin2 (' zeta_k in check error = *', zeta_k, 2)
ccc         call prin2 (' zeta_tar = *', zeta_tar, 2*ntar)
         call PRINF (' In CHECK_ERROR_TAR, NTAR = *', ntar, 1)
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
         call PRIN2 (' max abs error in solution = *', err, 1)
         call PRIN2 (' max rel error in solution = *', err/umax, 1)
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
      parameter (kmax = 500, npmax = 2048, nmax = kmax*npmax)
      parameter (nth_max = 1000, nphi_max = 1000, 
     1          ng_max = nth_max*nphi_max)
      parameter (nsp = 20*nmax + 20*ng_max)
c      
      common /geometry/ x_zeta, y_zeta, zeta, dzeta, dsda
      common /inteqn/ diag, cx, cy, cz, zeta_k
      common /sys_size/ k, nd, nbk
      common /fasblk2/ schur,wb,ipvtbf
      dimension schur(kmax*kmax),wb(kmax),ipvtbf(kmax)
      dimension x_zeta(nmax+ng_max), y_zeta(nmax+ng_max), dsda(nmax),  
     1          diag(nmax), cx(kmax), cy(kmax), cz(kmax)
      complex*16 zeta(nmax), dzeta(nmax), zeta_k(kmax), eye
c
         eye = dcmplx(0.d0,1.d0)
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
       call prin2 (' zk = *', zk, 2*k)
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
ccc      DO KBOD = 1,K
ccc	 call prinf(' column *',KBOD,1)
ccc	 call prin2(' SCHUR = *',schur(1,KBOD),K)
ccc      END DO
ccc      call PRINf (' ipvtbf = *', ipvtbf, K)
ccc      call PRIN2 (' wb = *', wb, k)
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
      parameter (kmax = 500, npmax = 2048, nmax = kmax*npmax)
      parameter (nth_max = 1000, nphi_max = 1000, 
     1          ng_max = nth_max*nphi_max)
      parameter (nsp = 20*nmax + 20*ng_max)
c      
      common /geometry/ x_zeta, y_zeta, zeta, dzeta, dsda
      common /inteqn/ diag_stereo, cx, cy, cz, zeta_k
      common /sys_size/ k, nd, nbk
      common /sphere_int/ xs, ys, zs, xn, yn, zn, diag
c
      dimension x_zeta(nmax+ng_max), y_zeta(nmax+ng_max), dsda(nmax),  
     1          diag_stereo(nmax), cx(kmax), cy(kmax), cz(kmax)
      dimension xs(nmax), ys(nmax), zs(nmax), xn(nmax), yn(nmax), 
     1          zn(nmax), diag(nmax)
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
c
ccc         call MATVEC_SPHERE (k, nd, nbk, xs, ys, zs, xn, yn, zn,
ccc     1                       dsda, diag, cx, cy, cz, A_k, xx, yy)
         call FASMVP (k, nd, nbk, nsp, x_zeta, y_zeta, zeta, dzeta,   
     1                zeta_k, diag_stereo, dsda, A_k, xx, yy, qa,   
     2                cfield, poten, wksp)
ccc         call FASMVP_TEST (k, nd, nbk, nsp, x_zeta, y_zeta, zeta,    
ccc     1                dzeta, zeta_k, diag, dsda, A_k, xx, yy, qa,   
ccc     2                cfield, poten, wksp)
ccc         call FASMVP (nd, nbk, k, x, y, dz, dsdth, rkappa, 
ccc     *                zk, xx, yy, h, za, pk, stress, dens)
ccc         call MSOLVE(N, yy, xx, NELT, IA, JA, A, ISYM, RWORK, IWORK)
ccc         call PRIn2 ('  xx = *', xx, n)
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
c periodic
            do j = 1, NY
               if (ireal.eq.1) then 
                  write(if,'(e20.13,$)')(ugrid(1,J))
                  write (if,'(a)')  ''
                 else
                  write(if,'(i4,$)') (igrid(1,j))
                  write (if,'(a)')  ''
               end if
            end do
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

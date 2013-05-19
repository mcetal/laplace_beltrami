      PROGRAM MAIN
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c
c  NDMAX is the maximum number of points allowed on each boundary
c  KMAX is the maximum number of holes plus outer boundary.
c  NMAX is the maximum number of unknowns
c  NCMAX is the maximum dimension of the unknown density vector.
c          It is twice the size of the maximum number of points, since
c          the density is a complex quantity.
c
c     IMPORTANT NOTE: By convention, we store complex vectors as a
c                     sequence of real,imaginary parts
c
      parameter (nmax = 5000)
c
c  maxl is the maximum nubmer of gmres iterations performed
c     before restarting.
c  lrwork is the dimension of a real workspace needed by dgmres.
c  liwork is the dimension of an integer workspace needed by dgmres.
c  rwork and iwork are work arrays used by dgmres
c
      parameter (maxl = 200)
      parameter (lrwork=2+nmax*(maxl+6)+maxl*(maxl+3), liwork=20)
      dimension rwork(lrwork),iwork(liwork)
c
c  X,Y are used to store coordinates of hole and outer boundaries
c  DZDT is dz/dt on boundary, where curve is equispaced with respect
c    to parameter t.
c  RKAPPA is the curvature
c
      DIMENSION X(nmax),Y(nmax),x2(2*nmax), y2(2*nmax), rkappa(nmax)
      complex*16 dzdt(nmax), z(nmax)
c
c  ZK is the hole centres
c  ZK0 is a point inside of the domain 
c  AI,BI describe the axes of the elliptical holes
c
      complex*16 zk, zk0, z0(nmax), C(nmax)
      dimension R_xx(nmax), R_x(nmax)
c
c  RHS stores the right-hand side of the integral equation,
c     namely W_x + i W_y. 
c  SOLN is the solution to the integral equation
c  UCMPLX is the complex form of the solution:
c     UCMPLX(I) = (SOLN(2*I-1),SOLN(2*i))
c
      DIMENSION RHS(2*nmax), SOLN(2*nmax), stream_ex(nmax)
      COMPLEX *16 UCMPLX(NMAX), u2(2*nmax), psi(nmax), chidl(nmax),
     1            chisl(nmax), phi(nmax), up(nmax), phi_ex(nmax), 
     2            psi_ex(nmax), phip(nmax), phip_ex(nmax)
      complex*16 eye, zshear, zb, zc, zc_exact
c
c  arrays needed by resampler
      dimension xy(2,nmax), der1(2,nmax), der2(2,nmax), 
     *          der3(2,nmax), fis(nmax), res_work(101*nmax)
c    
c  ZB stores the coefficients b_k in the Sherman-Lauricella equation.
c  ZC stores the coefficient of the C_k log (z-z_k) terms
c  zw,zn,zx,zu are used in calculating b_k and C_k
c
      COMPLEX *16 zw,zn,zx,zu,zval,boundary
c
c  workspace to check residual after solving the system
c
      dimension wmat(2*nmax)
      complex*16 dens(2*nmax),ddens(2*nmax), gradw(nmax)
c
c  FFTs
      dimension amode(nmax), fmode(nmax)
      complex*16 zf(nmax), zfi(nmax), wsave(20*nmax)
c
c   work arrays 
c
      parameter (nrwork = 2000000, ncwork = 2000000)
      complex*16 zwork(ncwork), work(nrwork)
c
c   double and single layer potentials
c
      dimension rho(nmax), f(nmax), dsdth(nmax), A_int(nmax), Gr(nmax)
      complex*16 qa(nmax), cfield(nmax), phicmplx(nmax)
c
c   Goursat functions
c
      dimension rechi_dl(nmax), rechi(nmax), rechi_calc(nmax), 
     1          stream(nmax), rechi_ex(nmax), err_chi(nmax)
c
c  other stuff
c
      dimension x0(nmax), c1(nmax), rtilde(nmax)
c
      PI = 4.0D0*DATAN(1.0D0)
      eye = dcmplx(0.d0,1.d0)
c 
c  Open plotting files
         call PRINI (6,13)
         open (unit = 10, file = 'boundary.m')
c
c  Initialize input, output files, get size of problem and 
c  choice of boundary conditions (ITEST)
c
         call SQUARE (nmax, nd, eps, bi, h, x, y, z, dzdt, dsdth, 
     1                  rkappa, xy, der1, der2, der3, fis, acc, err, 
     2                  res_work, ltot, aparam_length)
c
      STOP
      END
c
c
c---------------
      subroutine SQUARE (nmax, nd, eps, bi, h, x, y, z, dzdt, dsdth, 
     1                          rkappa, xy, der1, der2, der3, fis, acc,  
     2                              err, work, ltot, aparam_length)
c---------------
c  construct a resampled square
c
      implicit real*8 (a-h,o-z)
      dimension x(nmax), y(nmax), dsdth(nmax), rkappa(nmax), 
     1          xy(2,nmax), der1(2,nmax), der2(2,nmax), der3(2,nmax),
     2          fis(nmax), work(*)
      complex*16 dzdt(nmax), z(nmax)
c
         eye = dcmplx(0.d0,1.d0)
         pi = 4.d0*datan(1.d0)
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
         call prinf ('nd = *', nd, 1)
ccc         call RSPLOT (x, y, nd, 1, 10)
ccc        stop
cccc
cccc  test with an ellipse
ccc         nd = 512
ccc         ai = eps
ccc         do i = 1, nd
ccc            th = 2.d0*pi*(i-1)/nd
ccc            cs = dcos(th-thetax)
ccc            sn = dsin(th-thetax)
ccc            a2 = (eps)**2
ccc            b2 = (bi)**2
ccc            rnorm = dsqrt(b2*cs**2 + a2*sn**2)
ccc            radius = ai*bi/rnorm
ccc            rdot = -ai*bi*cs*sn*(-b2+a2)/rnorm**3
ccc            rddot =  -ai*bi*(2.d0*a2*b2*cs**2*sn**2
ccc     *                     +a2**2*cs**4 + a2*b2-a2**2+b2**2*cs**4
ccc     *                     -2.d0*b2**2*cs**2)/rnorm**5
ccc            x(i) = radius*dcos(th)
ccc            y(i) = radius*dsin(th)
ccc         end do
ccc         x(nd+1) = x(1)
ccc         y(nd+1) = y(1)
ccc         nd = nd+1
c
         NMIN = 64
         NDERS = 2
         n = 4096
c
         call RSRESA (ier,x,y,nd,nmin,nders,xy,der1,der2,
     *                der3,n,fis,h,acc,err,work,ltot)
         do i = 1,n
            x(i) = xy(1,i)
            y(i) = xy(2,i)
            z(i) = dcmplx(x(i),y(i))
            xdot = der1(1,i)
            ydot = der1(2,i)
            xddot = der2(1,i)
            yddot = der2(2,i)
            dsdth(i) = dsqrt(xdot**2 + ydot**2)
            dzdt(i) = dcmplx(xdot, ydot)
            rkappa(i) = -(-xdot*yddot+ydot*xddot)
         end do
         if (ier.ne.0) then
            call PRINF ('  IER = *',ier,1)
            stop
         end if
         call prin2 (' acc = *', acc, 1)
         call prin2 (' err = *', err, 1)
         nd = n
         call RSPLOT (x, y, nd, 1, 10)
ccc         stop
      return
      end

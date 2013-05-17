      PROGRAM STEREO_DRIVER
c     ------------------------------------------------------------------
c
c
c  Calculates FMM on stereographic projection of points on the sphere
c
c
c     ------------------------------------------------------------------
      implicit real*8 (a-h, o-z)
      parameter (npmax = 1000, nthmax = 1000, nmax = npmax*nthmax)
      parameter (nsp = 20*nmax)
c
c Geometry
      dimension theta(nmax), phi(nmax), x1(nmax), x2(nmax), x3(nmax),
     1          x(nmax), y(nmax)
      complex*16 z(nmax)
c
c Exact potential
      dimension phi_ex(nmax), phi_1e(nmax), phi_2e(nmax), phi_3e(nmax),
     1          phi_1(nmax), phi_2(nmax), phi_3(nmax)
      complex*16 zfield_ex(nmax)
c
c Fast Multipole Arrays
c
      complex*16 qa(nmax), cfield(nmax)
      dimension poten(nmax), wksp(nsp)
c
c Open output files
c
         open (unit = 11, file = 'blob.m')
c
         call GEO_INIT (nphi, nth, n, theta, phi, x1, x2, x3, z, qa)
         call RSCPLOT (z, n, 1, 11)
c
c Calculate exact potential
         call EXACT_POT (n, nsp, x, y, z, qa, phi_ex)
         call EXACT_FIELD (n, nsp, x, y, z, qa, zfield_ex)
         call GRADIENT_EX (n, x1, x2, x3, z, qa, phi_1, phi_2, phi_3,  
     1                     phi_1e, phi_2e, phi_3e) 
c
c Try FMM
         call FMM_TEST (n, nsp, x, y, z, x1, x2, x3, qa, cfield, 
     1                  poten, wksp)
c
c Construct gradient field
         call GRAD_BUILD (n, x1, x2, x3, z, cfield, phi_1, phi_2, 
     1                    phi_3) 

c
c check error
         errmax = 0.d0
         do i = 1, n
            err = cdabs(cfield(i)-zfield_ex(i))
            errmax = max(errmax, err)
         end do
         call PRIN2 (' error = *', errmax, 1)
c
      stop
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
      subroutine FMM_TEST (n, nsp, x, y, z, x1, x2, x3, qa, cfield, 
     1                     poten, wksp)
c---------------
c
      implicit real*8 (a-h, o-z)
      integer*4 iout(2), inform(10), ierr(10)
      complex*16 z(n), qa(n), cfield(n), zQsum, zQ2sum
      dimension xp(n), yp(n), poten(n), wksp(nsp), x1(n), x2(n), x3(n)
      REAL*4 TIMEP(2), ETIME
c
         pi = 4.D0*DATAN(1.D0)
         eye = DCMPLX(0.D0,1.D0)
c
         do i = 1, n
            qa(i) = qa(i)/(1-x3(i))
         end do
c
         tbeg = etime(timep)
         iout(1) = 6
         iout(2) = 13
         iflag7 = 2
         napb = 10
         ninire = 2
         mex = 300
         eps7 = 1.d-14
         tol = 1.d-14
         nnn = n
         call DAPIF2 (iout, iflag7, nnn, napb, ninire, mex, ierr, 
     &                inform, tol, eps7, x, y, qa, poten, cfield, 
     &                wksp, nsp, CLOSE)
         if (ierr(1).ne.0) then
            write (6,*) '  ERROR IN DAPIF2, IERR = ', (ierr(ii),ii=1,6)
            write(6,*) '  INFORM = ', (inform(ii),ii=1,6)
            stop
         end if
         call PRINF (' Number of Levels used = *', inform(3), 1)
         call PRIN2 (' cfield = *', cfield, 2*n)
c
c Fix up log
         zQsum = 0.d0
         zQ2sum = 0.d0
         do i = 1, n
            zQsum = zQsum + qa(i)
            zQ2sum = zQ2sum + qa(i)*dlog(1.d0+(cdabs(z(i)))**2)
         end do
         call PRIN2 (' zQsum = *', zQsum, 2)
         call PRIN2 (' zQ2sum = *', zQ2sum, 2)
         do i = 1, n
            poten(i) = 2.d0*poten(i)  
     1                  - cdabs(zQsum)*dlog(1.d0+(cdabs(z(i)))**2)
     2                  - dreal(zQ2sum)
            poten(i) = poten(i) 
     1                    + 2.d0*qa(i)*dlog(1.d0+(cdabs(z(i)))**2)
         end do
c
c Fix up field
         zQsum = 0.d0
         do i = 1, n
            zQsum = zQsum + qa(i)*dconjg(z(i))/(1.d0+(cdabs(z(i)))**2)
         end do
         do i = 1, n
            cfield(i) = -dconjg(cfield(i)) - zQsum 
     1                    + qa(i)*dconjg(z(i))/(1.d0+(cdabs(z(i)))**2)
            cfield(i) = 2.d0*cfield(i)
         end do
         call PRIN2 (' cfield = *', cfield, 2*n)
c
         tend = etime(timep)
ccc         call PRIN2 (' poten = *', poten, n)
         call PRIN2 (' TIME FOR FMM  = *',tend-tbeg,1)
ccc         call PRIN2 (' cfield = *', cfield, 2*n)
c
      return
      end
c
c*********************
c
      subroutine DUMP (nx,ny,vel,vort,press,igrid)
      implicit real*8 (a-h,o-z)
      complex*16 vel(nx,ny)
      dimension vort(nx,ny),press(nx,ny),igrid(nx,ny)
c
         open (unit = 44,file = 'u.data')
         open (unit = 45,file = 'v.data')
         open (unit = 46,file = 'vort.data')
         open (unit = 47,file = 'press.data')
ccc         open (unit = 48,file = 'grid.data')
         DO i = 1,NX
            write(44,'(e20.13,$)')(dreal(vel(I,J)),J=1,NY)
            write (44,'(a)')  ''
            write(45,'(e20.13,$)')(dimag(vel(I,J)),J=1,NY)
            write (45,'(a)')  ''
            write(46,'(e20.13,$)')(vort(I,J),J=1,NY)
            write (46,'(a)')  ''
            write(47,'(e20.13,$)')(press(I,J),J=1,NY)
            write (47,'(a)')  ''
ccc            write (48,'(i4,$)') (igrid(i,j),J=1,NY)
ccc            write (48,'(a)')  ''
         ENDDO
         close (44)
         close (45)
         close (46)
         close (47)
ccc         close (48)
         write (6,1001) vmin,vmax
c
 1001    format (' RANGE OF VORTICITY:  ',e12.6,2x,e12.6)
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

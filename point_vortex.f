      PROGRAM POINT_VORTEX
c     ------------------------------------------------------------------
c
c     point_vortex.f
c
c
c     Created on Wed Jan  9 09:13:52 2008
c     Copyright (c) 2008 MyCompany. All rights reserved.
c
c     Define an initial vorticity distribution at points inside 
c     an ellipse
c
c
c     ------------------------------------------------------------------
      implicit real*8 (a-h, o-z)
      parameter (npmax = 100000)
      parameter (nsp = 20*npmax)
c
c Geometry
c
      dimension xp(npmax), yp(npmax)
c
c Vorticity and velocity
c
      dimension vort(npmax), up(npmax), vp(npmax), up_dir(npmax), 
     &          vp_dir(npmax)
c
c Fast Multipole Arrays
c
      complex*16 qa(npmax), cfield(npmax)
      dimension poten(npmax), wksp(nsp)
c
c Open output files
c
         open (unit = 11, file = 'blob.m')
c
         call READINI (nx, ny, nt, nplot, Tfinal, delta, ai, bi)
         call PRINI (6, 13)
         call CREATE_ELLIPSE (nx, ny, np, dA, ai, bi, xp, yp)
         call PRINF (' NP = *', np, 1)
         call VORTICITY_DIST (np, xp, yp, vort)
         call RSPLOT (xp, yp, np, 1, 11)
c
c Time stepping loop
c
         dt = Tfinal/nt
         call PRIN2 ('   DT = *', dt, 1)
         do itime = 1, nt
            time = itime*dt
           call CALC_VEL_DIR(np, dA, xp, yp, vort, up, vp)
ccc            call CALC_VEL_FMP(np, nsp, dA, xp, yp, vort, qa, cfield,  
ccc     &                        poten, wksp, up, vp)
            call EULER (np, dt, xp, yp, up, vp)
            if (mod(itime,nplot).eq.0) then 
               call PRIN2 (' TIME = *', time, 1)
               call RSPLOT (xp, yp, np, 1, 11)
            end if
         end do
c         
      stop
      end     
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine READINI (nx, ny, nt, nplot, Tfinal, delta, ai, bi)
c---------------
c
      implicit real*8 (a-h,o-z)
c
         open (2, file = 'input_data')
         read (2,*) nx, ny
         read (2,*) ai, bi
         read (2,*) delta
         read (2,*) nt, nplot
         read (2,*) Tfinal
         close (2)
c
      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine CREATE_ELLIPSE (nx, ny, np, dA, ai, bi, xp, yp)
c---------------
c
      implicit real*8 (a-h,o-z)
      dimension xp(*), yp(*)
c
          dx = 2.d0*ai/nx
          dy = 2.d0*bi/ny
          ip = 0
          do i = 1, nx
             xgrid = -ai + i*dx
             do j = 1, ny
                ygrid = -bi + j*dy
                ellipse_check = (xgrid/ai)**2 + (ygrid/bi)**2
                if (ellipse_check.lt.1.d0) then
                   ip = ip + 1
                   xp(ip) = xgrid
                   yp(ip) = ygrid
                endif
             end do
          end do
          np = ip 
          dA = dx*dy 
c
      return 
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine VORTICITY_DIST(np, xp, yp, vort)
c---------------
c
      implicit real*8 (a-h,o-z)
      dimension xp(np), yp(np), vort(np)
c
         do i = 1, np
            rad = dsqrt((xp(i))**2 + (yp(i))**2)
            vort(i) = GAUSS(rad) 
ccc            vort(i) = 1.d0 
         end do
c
      return 
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      real*8 function GAUSS(r)
c---------------
c
      implicit real*8 (a-h,o-z)
c
          omega_0 = 1.d0
          delta = 0.5d0
          GAUSS = omega_0*dexp(-(r/delta)**2) 
c
      return 
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine CALC_VEL_DIR(np, dA, xp, yp, vort, up, vp)
c---------------
c
      implicit real*8 (a-h,o-z)
      dimension xp(np), yp(np), vort(np), up(np), vp(np), aKern(2)
c
         do i = 1, np
            up(i) = 0.d0
            vp(i) = 0.d0
            do j = 1, np
               if (j.ne.i) then
                  del_x = xp(i) - xp(j)
                  del_y = yp(i) - yp(j)
                  call FIELD (del_x,del_y,aKern)
                  up(i) = up(i) + dA*vort(j)*aKern(1)
                  vp(i) = vp(i) + dA*vort(j)*aKern(2)
               end if
            end do
         end do
ccc         call PRIN2 ('up direct = *', up, np)
ccc         call PRIN2 ('vp direct = *', vp, np)
c
      return 
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine FIELD (x,y,aK)
c---------------
c
      implicit real*8 (a-h,o-z)
      dimension aK(2)
c
         pi = 4.d0*datan(1.d0)
c
          r = dsqrt(x**2 + y**2)
          aK(1) = -y/(2.d0*pi*r**2)
          aK(2) = x/(2.d0*pi*r**2)
c
      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine EULER (np, dt, xp, yp, up, vp)
c---------------
c
      implicit real*8 (a-h,o-z)
      dimension xp(np), yp(np), up(np), vp(np)
c
         do i = 1, np
            xp(i) = xp(i) + dt*up(i)
            yp(i) = yp(i) + dt*vp(i)
         end do
c
      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine FMP_TEST (np, xp, yp, vort)
c---------------
c
      implicit real*8 (a-h, o-z)
      integer*4 iout(2), inform(10), ierr(10)
      complex*16 qa(1000), cfield(1000), cfield_dir(1000), eye, zi, zj, 
     &           zfield
      dimension xp(np), yp(np), vort(np), poten(1000), wksp(10000)
c
         pi = 4.D0*DATAN(1.D0)
         eye = DCMPLX(0.D0,1.D0)
c
         do i = 1, np
            qa(i) = 1.d0
         end do
c
         iout(1) = 6
         iout(2) = 13
         iflag7 = 2
         napb = 30
         ninire = 2
         mex = 300
         eps7 = 100
         tol = 1.d-14
         nnn = np
         nsp = 10000
         call DAPIF2 (iout, iflag7, nnn, napb, ninire, mex, ierr, 
     &                inform, tol, eps7, xp, yp, qa, poten, cfield, 
     &                wksp, nsp, CLOSE)
         if (ierr(1).ne.0) then
            write (6,*) '  ERROR IN DAPIF2, IERR = *', (ierr(ii),ii=1,6)
            write(6,*) '  INFORM = *', (inform(ii),ii=1,6)
            stop
         end if
         call PRIN2 (' CFIELD = *', cfield, 2*np)
c
         do i = 1, np
            cfield_dir(i) = 0.d0
            zi = dcmplx(xp(i),yp(i))
            do j = 1, i-1
               zj = dcmplx(xp(j), yp(j))
               zfield = (zj-zi)/(cdabs(zj-zi))**2
               cfield_dir(i) = cfield_dir(i) - zfield
            end do
            do j = i+1, np
               zj = dcmplx(xp(j), yp(j))
               zfield = (zj-zi)/(cdabs(zj-zi))**2
               cfield_dir(i) = cfield_dir(i) - zfield
            end do
         end do
         call PRIN2 (' CFIELD DIRECT = *', cfield_dir, 2*np)
c
      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine CALC_VEL_FMP (np, nsp, dA, xp, yp, vort, qa, cfield,  
     &                         poten, wksp, up, vp)
c---------------
c
      implicit real*8 (a-h,o-z)
      integer*4 iout(2), inform(10), ierr(10)
      complex*16 qa(np), cfield(np), eye, zi, zj
      dimension xp(np), yp(np), vort(np), poten(np), wksp(nsp), up(np), 
     &          vp(np)
      real*4 TIMEP(2), ETIME
      external CLOSE
c
      pi = 4.d0*datan(1.d0)
c
         do i = 1, np
            qa(i) = dA*vort(i)/(2.d0*pi)
         end do
c
         iout(1) = 6
         iout(2) = 13
         iflag7 = 2
         napb = 20
         ninire = 2
         mex = 300
         eps7 = 1.d-14
         tol = 1.d-7
         nnn = np
         T0 = ETIME(TIMEP)
         call DAPIF2 (iout, iflag7, nnn, napb, ninire, mex, ierr, 
     &                inform, tol, eps7, xp, yp, qa, poten, cfield, 
     &                wksp, nsp, CLOSE)
         T1 = ETIME(TIMEP)
         call PRINI (6, 13)
         call PRIN2 (' TIME IN DAPIF = *', T1 - T0, 1)
         if (ierr(1).ne.0) then
            write (6,*) '  ERROR IN DAPIF2, IERR = ', (ierr(ii),ii=1,6)
            write(6,*) '  INFORM = ', (inform(ii),ii=1,6)
            stop
         end if
ccc         call PRIN2 (' CFIELD = *', cfield, 2*np)
         do i = 1, np
            up(i) = -dimag(cfield(i))
            vp(i) = dreal(cfield(i))
         end do
ccc         call PRIN2 ('up = *', up, np)
ccc         call PRIN2 ('vp = *', vp, np)
c
      return 
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
        SUBROUTINE RSPLOT(X,Y,N,NBOD,IW)
        IMPLICIT REAL *8 (A-H,O-Z)
c---------------
c
        DIMENSION X(1),Y(1)
C
C        THIS SUBROUTINE PLOTS THE CURVE SPECIFIED BY ITS NODES
C        X,Y (CONCEPTUALLY, THE CURVE CONNECTS THE POINT (X(1),Y(1))
C        WITH THE POINT (X(N),Y(N)). 
C
C  INPUT PARAMETERS:
C   X,Y - THE COORDINATES OF THE CURVE TO PLOT
C   N - THE NUMBER OF ELEMENTS IN ARAYS X,Y.
C   IW - THE FORTRAN UNIT NUMBER ON WHICH THE OUTPUT DATA SET IS WRITTEN
C
C  OUTPUT PARAMETERS: NONE
C
C
      WRITE(IW,*) 'x = ['
      WRITE(IW,1400) (X(I),Y(I),I=1,N)
ccc      write (iw,1400) x(1),y(1)
      WRITE(IW,*) '];'
1300  FORMAT(D15.6)
1400  FORMAT(2D15.6)
      write(iw,*) ' xx = x(:,1);'
      write(iw,*) ' yy = x(:,2);'
ccc      write(iw, *) ' plot(xx,yy)'
      write(iw, *) " plot(xx,yy,'o')"
      write(iw,*) 'axis equal'
      write(iw, *) 'pause'
ccc      IF (NBOD.eq.1) write(iw,*) ' hold on'
      RETURN
      end
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

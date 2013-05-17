close all; clear;
%
% entire world
NX = 500; NY = 500;
load xzeta_grid_world.dat
 a = zeros(NX,NY);
 a(:) = xzeta_grid_world(:);
 xzeta_grid_world = a;
load yzeta_grid_world.dat
 a = zeros(NX,NY);
 a(:) = yzeta_grid_world(:);
 yzeta_grid_world = a;
load xgrid_world.dat
 a = zeros(NX,NY);
 a(:) = xgrid_world(:);
 xgrid_world = a;
load ygrid_world.dat
 a = zeros(NX,NY);
 a(:) = ygrid_world(:);
 ygrid_world = a;
load zgrid_world.dat
 a = zeros(NX,NY);
 a(:) = zgrid_world(:);
 zgrid_world = a;
load ugrid_world.dat
 a = zeros(NX,NY);
 a(:) = ugrid_world(:);
 ugrid_world = a;
%
% close up of continents
load xzeta_grid.dat
 a = zeros(NX,NY);
 a(:) = xzeta_grid(:);
 xzeta_grid = a;
load yzeta_grid.dat
 a = zeros(NX,NY);
 a(:) = yzeta_grid(:);
 yzeta_grid = a;
load xgrid.dat
 a = zeros(NX,NY);
 a(:) = xgrid(:);
 xgrid = a;
load ygrid.dat
 a = zeros(NX,NY);
 a(:) = ygrid(:);
 ygrid = a;
load zgrid.dat
 a = zeros(NX,NY);
 a(:) = zgrid(:);
 zgrid = a;
load ugrid.dat
 a = zeros(NX,NY);
 a(:) = ugrid(:);
 ugrid = a;
%
%
umin = -6.
umax = 6.
figure(1)
%subplot(1,2,1)
   hold on
   %surf(xgrid_world,ygrid_world,zgrid_world,ugrid_world)
   surf(xgrid,ygrid,zgrid,ugrid)
      shading interp
   jet2 = colormap(jet);
   jet2 = [0. 0. 0.;jet2];
   colormap(jet2)
   hold on
   geo_3d
   view (-180,25)
   lighting phong
   camlight('headlight')
   camlight('right')
   material dull
   %light
   caxis([umin umax])
   axis equal
   axis off
figure(2)
%subplot(1,2,2)
   hold on
   %surf(xgrid_world,ygrid_world,zgrid_world,ugrid_world)
   surf(xgrid,ygrid,zgrid,ugrid)
      shading interp
   jet2 = colormap(jet);
   jet2 = [0. 0. 0.;jet2];
   colormap(jet2)
   hold on
   geo_3d
   view (71,29)
   lighting phong
   camlight('headlight')
   camlight('right')
   material dull
   %light
   caxis([umin umax])
   axis equal
   axis off
figure(3)
vc = [umin:.1:umax];
contour (xzeta_grid,yzeta_grid,ugrid,vc)
hold on
coast
%axis ([-2 2 -2 2])
close all; clear;
%
% Define stereographic projection
   stereo = inline('(x+1i*y)./(1-z)','x','y','z');

load xy_africa.dat;
[x,y,z]=continent(xy_africa); z = -z;
%
% plot the world
figure(1)
   sphere
   colormap([0.5 0.5 0.5])
   shading flat
   alpha(0.5)
   hold on
   plot3(0,0,1,'r*')
   plot3(x,y,z,'k','LineWidth',2)
load xy_mpoly_antarctica.dat;
[x,y,z]=continent(xy_mpoly_antarctica); z = -z;
   plot3(x,y,z,'k','LineWidth',2)
load xy_mpoly_australia.dat;
[x,y,z]=continent(xy_mpoly_australia);z = -z;
   plot3(x,y,z,'k','LineWidth',2)
%load xy_baffin.dat;
%[x,y,z]=continent(xy_baffin);
   %plot3(x,y,z,'k','LineWidth',2)
load xy_mpoly_europe_africa.dat;
[x,y,z]=continent(xy_mpoly_europe_africa);z = -z;
   plot3(x,y,z,'k','LineWidth',2)
load xy_mpoly_northamerica.dat;
[x,y,z]=continent(xy_mpoly_northamerica);z = -z;
   plot3(x,y,z,'k','LineWidth',2)
%
%
figure(2)
plot(xy_mpoly_antarctica(:,1),xy_mpoly_antarctica(:,2),'k')
hold on
plot(xy_mpoly_australia(:,1),xy_mpoly_australia(:,2),'k')
plot(xy_mpoly_europe_africa(:,1),xy_mpoly_europe_africa(:,2),'k')
plot(xy_mpoly_northamerica(:,1),xy_mpoly_northamerica(:,2),'k')
   

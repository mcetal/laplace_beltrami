   close all; clear;
%
% Antarctica
%   worldmap('antarctica')
   antarctica = shaperead('landareas', 'UseGeoCoords', true,...
  'Selector',{@(name) strcmp(name,'Antarctica'), 'Name'});
   theta = 90 - antarctica.Lat;
   theta = theta*2*pi/360;
   theta = theta(1:end-1);
   phi = antarctica.Lon;
   phi = phi*2*pi/360;
   phi = phi(1:end-1);
   x_ant = cos(phi).*sin(theta);
   y_ant = sin(phi).*sin(theta);
   z_ant = cos(theta);
   n = size(z_ant);
   disp(['This map has ',num2str(n(2)),' points'])
%
% North America
   NorthAmerica = shaperead('landareas', 'UseGeoCoords', true,...
  'Selector',{@(name) strcmp(name,'NorthAmerica'), 'Name'});
   theta = 90 - NorthAmerica.Lat;
   theta = theta*2*pi/360;
   theta = theta(1:end-1);
   phi = NorthAmerica.Lon;
   phi = phi*2*pi/360;
   phi = phi(1:end-1);
   x_NA = cos(phi).*sin(theta);
   y_NA = sin(phi).*sin(theta);
   z_NA = cos(theta);
   n = size(z_NA);
   disp(['This map has ',num2str(n(2)),' points'])

%
% plot the world
   sphere
   colormap([0.5 0.5 0.5])
   shading flat
   alpha(0.5)
   hold on
   plot3(x_ant,y_ant,z_ant,'k','LineWidth',2)
   plot3(x_NA,y_NA,z_NA,'k','LineWidth',2)
%
% Define stereographic projection
   stereo = inline('(x+1i*y)./(1-z)','x','y','z');
%
% Calculate and plot stereographic projection
   zeta = stereo(x_ant,y_ant,z_ant);
   figure(2)
   plot(real(zeta),imag(zeta),'k')
   r_zeta = real(zeta'); i_zeta = imag(zeta');
   hold on
   plot(r_zeta(1),i_zeta(1),'r*')
   plot(r_zeta(end),i_zeta(end),'b*')
   save('r_zeta.dat','r_zeta','-ascii','-double')
   save('i_zeta.dat','i_zeta','-ascii','-double')
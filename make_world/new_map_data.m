   close all; clear;
   load coast
%
% map data in lat and long
   theta = 90 - lat;
   theta = theta*2*pi/360;
   phi = long;
   phi = phi*2*pi/360;
   x = cos(phi).*sin(theta);
   y = sin(phi).*sin(theta);
   z = cos(theta);
   n = size(z);
   disp(['This map has ',num2str(n(1)),' points'])
%
% plot the world
   sphere
   colormap([0.5 0.5 0.5])
   shading flat
   alpha(0.5)
   hold on
   plot3(x, y, z,'k','LineWidth',2)
%
% Define stereographic projection
   stereo = inline('(x+1i*y)./(1-z)','x','y','z');
%
% Calculate and plot stereographic projection
   zeta = stereo(x,y,z);
   figure(2)
   plot(real(zeta),imag(zeta),'k')
   r_zeta = real(zeta'); i_zeta = imag(zeta');
   hold on
   plot(r_zeta(1),i_zeta(1),'r*')
   plot(r_zeta(end),i_zeta(end),'b*')

%
% Loop through all the islands
   test = isnan(theta);
   ic = 1
   for i=1:n(1)
       i
       test(i)
       ic
       if (test(i)==0)
           xi(ic) = x(i);
           yi(ic) = y(i);
           zi(ic) = z(i);
           ic = ic+1;
       else
           if (ic>50) 
           figure(3)
           subplot(1,2,1)
           sphere
              colormap([0.5 0.5 0.5])
              shading flat
              alpha(0.5)
              hold on
              plot3(x,y,z,'k')
              plot3(xi(1:ic-1), yi(1:ic-1), zi(1:ic-1),'r','LineWidth',2)
              zeta = stereo(xi,yi,zi);
              r_zeta = real(zeta'); i_zeta=imag(zeta');
           subplot(1,2,2)
              plot(r_zeta(1:ic-1),i_zeta(1:ic-1))
              xy = zeros(ic-1,2);
              xy(:,1) = r_zeta(1:ic-1); xy(:,2) = i_zeta(1:ic-1);
           pause
           close all
              save('xy.dat','xy','-ascii')
           end
           ic = 1
        end
   end
       

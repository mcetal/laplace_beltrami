function [x,y,z]=continent(xy)
   zeta = xy(:,1)+1i*xy(:,2);
   x = (zeta+conj(zeta))./(1+abs(zeta).^2);
   y = (zeta-conj(zeta))./(1i*(1+abs(zeta).^2));
   z = (-1+abs(zeta).^2)./(1+abs(zeta).^2);

%% Mimimum size of element
function [lmin_wtharm,lmin_h_x,lmin_h_y] = minimum_l_fem(E,rho,nu,fmax,el_side_x,el_side_y,m,n,theta,phi,c_air)
% E, rho, nu are properties of the material
% nn number of harmonics

wmax = 2*pi*fmax;

%pressure wave
% aux = real(E)*(1-nu)/((1+nu)*(1-2*nu));
% cp = sqrt(aux/rho);

%shear wave
aux2 = real(E)/(2*(1+nu));
cs = sqrt(aux2/rho);
  
k_max = wmax/cs;
lambda = 2*pi/k_max;
 
lmin_wtharm = lambda/6;

%Structural wavenumbers with harmonics in x direction
k=wmax/c_air;

kx = k*sin(theta)*cos(phi); 
k_harm_x = max(kx)+2*pi*m(end);
lambda_h_x = 2*pi/k_harm_x; 
lmin_h_x = lambda_h_x/6;


if el_side_x  > lmin_h_x
   disp('FE discretization in x direction is not appropriate for these quantitify of harmonics') 
end
 
%Structural wavenumbers with harmonics in y direction
ky = k*sin(theta)*sin(phi); 
k_harm_y = max(ky)+2*pi*n(end);
lambda_h_y = 2*pi/k_harm_y; 
lmin_h_y = lambda_h_y/6;


if el_side_y  > lmin_h_y
   disp('FE discretization in y direction is not appropriate for these quantitify of harmonics') 
end

end
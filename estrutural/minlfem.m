%% Mimimum size of element
function lmin = minlfem(E,h,rho, nu,fmax)
wmax = 2*pi*fmax;
D = (E*h^3)/(12*(1-nu^2));
 k = (((wmax^2)*rho*h)/D)^0.25;
 lambda = (2*pi)/k;
 lmin = lambda/6;
 
end
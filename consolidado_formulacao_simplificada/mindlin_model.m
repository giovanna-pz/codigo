% Model 2 - Mindlin Ferreira 2008
I =h^3/12;

% kapa=0.8601; % cccc / cccf case
% kapa=0.822; % scsc case
 %kapa=5/6; % ssss case
kapa=1; % teste
 % bending part
C_bending=I*E/(1-nu^2)*...
[1 nu 0;nu 1 0;0 0 (1-nu)/2];
% shear part
C_shear=kapa*h*E/2/(1+nu)*eye(2);

[KG,MG]= mindlin_stiffness_mass_matrices(GDof,nelem,elem,nnode,nodes,h,rho,I,C_shear,C_bending);
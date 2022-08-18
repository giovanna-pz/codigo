function [omega,N,modes_full]=dispersion_2D_rectangular(Lx,Ly,nc,K_new,M_new,int_dofs,act_dofs,dof,fatx,faty)
%Lx,Ly: cell size
%int_dofs:vector of internal dofs
%act_dofs:vector of active dofs


% %% Dispersion relation - considering a square lattice
% Lx = Ly  
 
% N points of discretization - border of the first brillouin zone
n = 30; % n intervals per line;
N = 4*n;

%allowed frequencies - solution of the eigenproblem

t= nc + length(int_dofs); % dofs of qc + internal dofs

omega = zeros(t,N);

%Verificar como vou guardar os modos
% modes = zeros(t,t,N); 
%Para cada ponto eu tenho uma matriz t,t -> cada omega tem um autovetor
GDof= length([int_dofs,act_dofs]);
 modes_full = zeros(GDof,t,N); 
%each plane (x,y) are the wavemodes for one wavenumber N

%line T to X
for i = 1:n+1
    ky =0;
    kx = (i-1)*((pi/Lx)/n);
    [omega_aux,V_full] = redKM(dof,kx,ky,Lx,Ly, K_new, M_new, int_dofs, act_dofs,fatx,faty,t) ;
    modes_full(:,:,i)= V_full;
    omega(:,i) = omega_aux;
end

%line X to M
for j= n+1:2*n

 kx =pi/Lx;
 ky = (j-n)*((pi/Ly)/n);
 
 [omega_aux,V_full] = redKM(dof,kx,ky,Lx,Ly, K_new, M_new, int_dofs, act_dofs,fatx,faty,t) ;
 modes_full(:,:,j)= V_full;
 omega(:,j) = omega_aux;
 
end

%line M to K
for k= (2*n+1):(3*n)

 kx =((3*n)-k)*((pi/Lx)/n);
 ky = pi/Ly;
 
[omega_aux,V_full] = redKM(dof,kx,ky,Lx,Ly, K_new, M_new, int_dofs, act_dofs,fatx,faty,t);
 modes_full(:,:,k)= V_full;
 omega(:,k) = omega_aux;
 
end


%line K to T
for kk= (3*n+1):(4*n)

 ky =((4*n)-kk)*((pi/Ly)/n);
 kx = 0;
 
[omega_aux,V_full] = redKM(dof,kx,ky,Lx,Ly, K_new, M_new, int_dofs, act_dofs,fatx,faty,t);
 modes_full(:,:,kk)= V_full;
 omega(:,kk) = omega_aux;
 
end


end
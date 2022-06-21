 function [omega,V_full] = redKM(dof,kx,ky,Lx,Ly, K, M, int_dofs, act_dofs,fat,t) 
%function [omega,V] = redKM_3(dof,kx,ky,Lx,Ly, K, M, int_dofs, act_dofs) 
% function to calculate the reduced matrcies for a given problem with
% periodic condition
% matrices lambda_R and lambda_L are determined case by case

    mu_x = kx*Lx;
    mu_y = ky*Ly;
    
    lambda_x = exp(-1i*mu_x);
    lambda_y = exp(-1i*mu_y);
    
    
    %The Lambda_R matrix in the article contains typography errors
%Matriz Giovanna
Lbd_R = [         eye(dof,dof), zeros(dof,fat*dof), zeros(dof,fat*dof);...
            lambda_x*eye(dof,dof), zeros(dof,fat*dof), zeros(dof,fat*dof);...
            lambda_y*eye(dof,dof), zeros(dof,fat*dof), zeros(dof,fat*dof);...
            lambda_x*lambda_y*eye(dof,dof), zeros(dof,fat*dof), zeros(dof,fat*dof);
            zeros(fat*dof,dof), eye(fat*dof,fat*dof),zeros(fat*dof,fat*dof);
            zeros(fat*dof,dof), zeros(fat*dof,fat*dof), eye(fat*dof,fat*dof);
            zeros(fat*dof,dof),lambda_x*eye(fat*dof,fat*dof), zeros(fat*dof,fat*dof);
            zeros(fat*dof,dof), zeros(fat*dof,fat*dof),lambda_y*eye(fat*dof,fat*dof)];
         

    Lbd_L = Lbd_R'; %conjugate transpose, lambda_r is complex
      
    %Separating the matrices components (if necessary)
    Kaa = K(act_dofs, act_dofs);
    Kai = K(act_dofs, int_dofs);
    Kia = K(int_dofs,act_dofs);
    Kii = K(int_dofs, int_dofs);
    
    Maa = M(act_dofs, act_dofs);
    Mai = M(act_dofs, int_dofs);
    Mia = M(int_dofs,act_dofs);
    Mii = M(int_dofs, int_dofs); 
  
  
  %% Eigenproblem formulation
    G1 = [Lbd_L*Kaa*Lbd_R, Lbd_L*Kai;
          Kia*Lbd_R, Kii];
    G2 = [Lbd_L*Maa*Lbd_R, Lbd_L*Mai;
         Mia*Lbd_R, Mii];
    
    
    [V1,D] =eig(G1,G2);
    [omega_2,index] = sort(diag(D));
    omega = sqrt(omega_2);
    
    V = V1(:,index); %V= t x t
  
%     % Recovering the global eigenvector (for plots)
%     V_full = GDof+1,t
    na = length(act_dofs);
    ni = length(int_dofs);
    nc = 2*fat*dof + dof; %nc is determined case by case dofs of qc
    A = [Lbd_R, zeros(na,ni);
    zeros(ni,nc), eye(ni,ni)];
   
   for kk = 1:t
   V_full(:,kk) = A*V(:,kk);
        
   end
end
%Not condensed version of the WFE problem
%Partitioning the dynamic stiffness matrix - equation 14         
D_part = [D_til_bb,D_til_bi;
       D_til_ib,D_til_ii];
   
Lambda_R_p = [Lambda_R, zeros(na,ni);
    zeros(ni,nc), eye(ni,ni)];

Lambda_L_p = [Lambda_L, zeros(nc,ni);
    zeros(ni, na),eye(ni,ni)];

A_p = Lambda_L_p*D_part*Lambda_R_p;

%Partitioning the external force vector
e_p = [Fext_b;Fext_i];

B_p = Lambda_L_p*e_p;


% Solving Equation 42
Sol = A_p\B_p;

%q = [qc,qi]
% the order here is (q1,q2,q3,q4,ql,qb,qr,qt)
q_bound = Lambda_R*Sol(1:nc);
q_internal = Sol(nc+1:end);
%%%%%%% end of the not condensed WFE problem

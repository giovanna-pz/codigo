% Dynamic condensed version of the WFE problem

%Dynamic condensation - Equation 34
%Fluid loading effects added
D_r_til = D_til_bb-D_til_bi*(D_til_ii\D_til_ib);
F_r_til = Fext_b - (D_til_bi/D_til_ii)*Fext_i;

%Equation 43
D_t = Lambda_L*D_r_til*Lambda_R;
e_t = Lambda_L*F_r_til;

%Solving - Equation 42

Sol = D_t\e_t;

%q = [qc,qi]
% the order here is (q1,q2,q3,q4,ql,qb,qr,qt)
%Equation 38
q_bound = Lambda_R*Sol;

%Recovering internal nodes
%Equation 

q_internal = (D_til_ii\Fext_i) -  (D_til_ii\D_til_ib)*q_bound;


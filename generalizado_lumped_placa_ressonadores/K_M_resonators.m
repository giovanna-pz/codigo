function [K_new, M_new, numberRes] = K_M_resonators(KG,MG,node_res_first,dof,GDof,kr,mr)

%This function adds the degrees of freedom of the resonators and calculates
%the mass and stiffness matrices for an strucuture
%with attached resonators

% Input
%KG and MG: Mass and Stiffness matrices for the substrate material/strcuture
%nelem: number of finite elements in the structure
%node_res_first: first row with resonators in the bottom part of the
%structure
%Ngrid: number of nodes in one side of the strcuture
%dof: number of dofs per node
%GDof: number of dofs in the strucuture

numberRes = length(kr);
%lr = length(node_res_first);

node_res = node_res_first;

% for g=1:lr-1
%     
%     l_aux = Ngrid*2*ones(1,lr)+node_res(g,:);
%     node_res =[node_res; l_aux];
%     
% end

node_res = sort(reshape(node_res,[1,numberRes]));
z_dof_res = (node_res-ones(1,numberRes))*dof+ones(1,numberRes);
dof_res = GDof+1:GDof+numberRes;


% increasing matrix size
K_new= zeros(GDof+numberRes);
M_new= zeros(GDof+numberRes);

K_new(1:GDof, 1: GDof) = KG;
M_new(1:GDof, 1: GDof) = MG;

% inserting data from spring-mass system
for rr=1:numberRes
    K_new(z_dof_res(rr),z_dof_res(rr)) = K_new(z_dof_res(rr),z_dof_res(rr))+kr(rr);

    K_new(z_dof_res(rr),dof_res(rr)) = K_new(z_dof_res(rr),dof_res(rr)) - kr(rr);

    K_new(dof_res(rr),z_dof_res(rr)) = K_new(dof_res(rr),z_dof_res(rr))-kr(rr);

    K_new(dof_res(rr),dof_res(rr)) = K_new(dof_res(rr),dof_res(rr))+kr(rr);
    
    M_new(dof_res(rr),dof_res(rr)) = M_new(dof_res(rr),dof_res(rr))+mr(rr);

end


end
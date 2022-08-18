% Partitioning of the problem - used in WFE model
% Notation
%%3--- qt ----4
%ql    qi     qr
%%1--- qb ----2

%Boundary nodes - change if I change the number of elements
q1 = 1;
q2 = NgridX;
q4 = nnode;
q3 = nnode-NgridX+1;
ql= NgridX+1:NgridX:q3-NgridX;
qb = 2:1:NgridX-1;
qr = q2+NgridX:NgridX:q4-NgridX;
qt = q3+1:1:q4-1;

q_all = 1:1:nnode; % Vector with all nodes
qB = [q1,q2,q3,q4,ql,qb,qr,qt]; % boundary nodes

% internal nodes
qi = setdiff(q_all,qB); 

% Describing internal and boundary dofs

 dofs_b = findNodeDofs(qB,dof); %boundary dofs
 dofs_i = sort(findNodeDofs(qi,dof)); %internal dofs - increasing order
 
 dofs_res = GDof+1:1:GDof+numberRes; % dofs of resonators
 
 % updating internal dofs with resonators
 dofs_i = [dofs_i, dofs_res];
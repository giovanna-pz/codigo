function [index] = findNodeDofs(nodes,ndof)

%% [index] = findNodeDofs(nodes,ndof)
%
% Computes global dofs associated to a list of nodes
%
% index - vector of global dofs associated to nodes
% nodes - nodes numbers whose dof must be determined
% ndof - nuber of dofs per node

nnodes = length(nodes);
index = zeros(1,nnodes*ndof);

for i=1:nnodes % nodes loop
    
    ni = nodes(i);
    index ( (i-1)*ndof+1 : i*ndof ) = (ni-1)*ndof+1 : ni*ndof ;

end
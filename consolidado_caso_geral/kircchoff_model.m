%kircchoff_model

% %Model 1 - Kirchoff plate - 2022
materials = [ E nu rho ];
%inci = [ node1, node 2, node3, node4, type of material] 
inci = [elem, ones(nelem,1)];
% nodes_with_z = [xcoordinate, ycoordinate, z coordinate]
nodes_with_z = [nodes, zeros(nnode,1)];
% Mass ans Stiffness matrices of the plate using Kirchoff model
[KG,MG] = fePlateMatrixAssembly(nodes_with_z,inci,h,materials,'kirchhoff');
function [stiff,mass]= mindlin_stiffness_mass_matrices(GDof,nelem,elem,nnode,nodes,h,rho,I,C_shear,C_bending) 

[stiff]=...
formStiffnessMatrixMindlinQ4(GDof,nelem,...
elem,nnode,nodes,C_shear,...
C_bending,h,I);

[mass]=...
formMassMatrixMindlinQ4(GDof,nelem,...
elem,nnode,nodes,h,rho,I);


end 

%................................................................
function [mass]=...
formMassMatrixMindlinQ4(GDof,numberElements,...
elementNodes,numberNodes,nodeCoordinates,thickness,rho,I)
% computation of mass matrix
% for Mindlin plate element
% mass : mass matrix
mass=zeros(GDof);
% Gauss quadrature for bending part
[gaussWeights,gaussLocations]=gaussQuadrature('complete');
% cycle for element
for e=1:numberElements
% indice : nodal condofectivities for each element
indice=elementNodes(e,:);
ndof=length(indice);
% cycle for Gauss point
for q=1:size(gaussWeights,1)
GaussPoint=gaussLocations(q,:);
xi=GaussPoint(1);
eta=GaussPoint(2);
% shape functions and derivatives
[shapeFunction,naturalDerivatives]=shapeFunctionQ4(xi,eta);
% Jacobian matrix, inverse of Jacobian,
% derivatives w.r.t. x,y
[Jacob,~,XYderivatives]=...
Jacobian(nodeCoordinates(indice,:),naturalDerivatives);
% [B] matrix bending
B_b=zeros(3,3*ndof);
B_b(1,ndof+1:2*ndof) = XYderivatives(:,1)';
B_b(2,2*ndof+1:3*ndof)= XYderivatives(:,2)';
B_b(3,ndof+1:2*ndof) = XYderivatives(:,2)';
B_b(3,2*ndof+1:3*ndof)= XYderivatives(:,1)';
% mass matrix
mass(indice,indice)=mass(indice,indice)+...
shapeFunction*shapeFunction'*thickness*...
rho*gaussWeights(q)*det(Jacob);
mass(indice+numberNodes,indice+numberNodes)=...
mass(indice+numberNodes,indice+numberNodes)+...
shapeFunction*shapeFunction'*I*...
rho*gaussWeights(q)*det(Jacob);
mass(indice+2*numberNodes,indice+2*numberNodes)=...
mass(indice+2*numberNodes,indice+2*numberNodes)+...
shapeFunction*shapeFunction'*I*...
rho*gaussWeights(q)*det(Jacob);
end % Gauss point
end % element

end

%................................................................
function [K]=...
formStiffnessMatrixMindlinQ4(GDof,numberElements,...
elementNodes,numberNodes,nodeCoordinates,C_shear,...
C_bending,thickness,I)
% computation of stiffness matrix
% for Mindlin plate element
% K : stiffness matrix
K=zeros(GDof);
% Gauss quadrature for bending part
[gaussWeights,gaussLocations]=gaussQuadrature('complete');
% cycle for element
% cycle for element
for e=1:numberElements
% indice : nodal condofectivities for each element
% elementDof: element degrees of freedom
indice=elementNodes(e,:);
elementDof=[indice indice+numberNodes indice+2*numberNodes];
ndof=length(indice);
% cycle for Gauss point
for q=1:size(gaussWeights,1)
GaussPoint=gaussLocations(q,:);
xi=GaussPoint(1);
eta=GaussPoint(2);
% shape functions and derivatives
[shapeFunction,naturalDerivatives]=shapeFunctionQ4(xi,eta);
% Jacobian matrix, inverse of Jacobian,
% derivatives w.r.t. x,y
[Jacob,~,XYderivatives]=...
Jacobian(nodeCoordinates(indice,:),naturalDerivatives);
% [B] matrix bending
B_b=zeros(3,3*ndof);
B_b(1,ndof+1:2*ndof) = XYderivatives(:,1)';
B_b(2,2*ndof+1:3*ndof)= XYderivatives(:,2)';
B_b(3,ndof+1:2*ndof) = XYderivatives(:,2)';
B_b(3,2*ndof+1:3*ndof)= XYderivatives(:,1)';
% stiffness matrix bending
K(elementDof,elementDof)=K(elementDof,elementDof)+ ...
B_b'*C_bending*B_b*gaussWeights(q)*det(Jacob);
end % Gauss point
end % element
% shear stiffness matrix
% Gauss quadrature for shear part
[gaussWeights,gaussLocations]=gaussQuadrature('reduced');
% cycle for element
% cycle for element
for e=1:numberElements
% indice : nodal condofectivities for each element
% elementDof: element degrees of freedom
indice=elementNodes(e,:);
elementDof=[indice indice+numberNodes indice+2*numberNodes];
ndof=length(indice);
% cycle for Gauss point
for q=1:size(gaussWeights,1)
GaussPoint=gaussLocations(q,:);
xi=GaussPoint(1);
eta=GaussPoint(2);
% shape functions and derivatives
[shapeFunction,naturalDerivatives]=shapeFunctionQ4(xi,eta);
% Jacobian matrix, inverse of Jacobian,
% derivatives w.r.t. x,y
[Jacob,~,XYderivatives]=...
Jacobian(nodeCoordinates(indice,:),naturalDerivatives);
% [B] matrix shear
B_s=zeros(2,3*ndof);
B_s(1,1:ndof) = XYderivatives(:,1)';
B_s(2,1:ndof) = XYderivatives(:,2)';
B_s(1,ndof+1:2*ndof) = shapeFunction;
B_s(2,2*ndof+1:3*ndof)= shapeFunction;
% stiffness matrix shear
K(elementDof,elementDof)=K(elementDof,elementDof)+...
B_s'*C_shear *B_s*gaussWeights(q)*det(Jacob);
end % gauss point
end % element

end

% .............................................................
function [JacobianMatrix,invJacobian,XYDerivatives]=...
Jacobian(nodeCoordinates,naturalDerivatives)
% JacobianMatrix : Jacobian matrix
% invJacobian : inverse of Jacobian Matrix
% XYDerivatives : derivatives w.r.t. x and y
% naturalDerivatives : derivatives w.r.t. xi and eta
% nodeCoordinates : nodal coordinates at element level
JacobianMatrix=nodeCoordinates'*naturalDerivatives;
invJacobian=inv(JacobianMatrix);
XYDerivatives=naturalDerivatives*invJacobian;
end % end function Jacobian

% .............................................................
function [shape,naturalDerivatives]=shapeFunctionQ4(xi,eta)
% shape function and derivatives for Q4 elements
% shape : Shape functions
% naturalDerivatives: derivatives w.r.t. xi and eta
% xi, eta: natural coordinates (-1 ... +1)
shape=1/4*[ (1-xi)*(1-eta);(1+xi)*(1-eta);
(1+xi)*(1+eta);(1-xi)*(1+eta)];
naturalDerivatives=...
1/4*[-(1-eta), -(1-xi);1-eta, -(1+xi);
1+eta, 1+xi;-(1+eta), 1-xi];
end % end function shapeFunctionQ4


function [weights,locations]=gaussQuadrature(option)
% Gauss quadrature for Q4 elements
% option ‘complete’ (2x2)
% option ‘reduced’ (1x1)
% locations: Gauss point locations
% weights: Gauss point weights
switch option
case 'complete'
locations=...
[ -0.577350269189626 -0.577350269189626;
0.577350269189626 -0.577350269189626;
0.577350269189626 0.577350269189626;
-0.577350269189626 0.577350269189626];
weights=[ 1;1;1;1];
case 'reduced'
locations=[0 0];
weights=[4];
end
end % end function gaussQuadrature

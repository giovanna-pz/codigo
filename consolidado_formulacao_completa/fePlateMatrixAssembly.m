function [K,M] = fePlateMatrixAssembly(nodes,elems,h,materials,plateType)
%
%% [K,M] = fePlateMatrixAssembly(nodes,elems,h,materials,plateType)
%
% assembles stiffness and mass matrix
%
% K - stiffness matrix
% M - mass matrix
% nodes - nodal coordinate values
% elems - nodal connectivity
% h - plate thickness
% materials - material properties [E nu rho]
% plateType - 'mindlin', 'kirchhoff', 'dkt'

% obs.: Nc - funcao de forma para coordenadas (e Jacobiano)
%       Nw - funcao de forma para deslocamentos

% number of dofs per node
% Kirchoff (thin): w, dw/dx, dw/dy
% Mindlin (thick): w, psi_x, psi_y

% choice of integration points
nGL = 2;
[points,weights] = feGaussQuadrature2(nGL,'quad');

% dimensions
nDofsNode = 3;  % number of dofs per node
nNodesElem = 4; % number of nodes per element

nNodes = size(nodes,1); % total number of nodes in the system
nElems = size(elems,1); % number of elements

nDofsSys = nNodes*nDofsNode;      % total number of system dofs

%% material properties

nMaterials = size(materials,1); % number of material properties
Db = cell(1,nMaterials);         % preallocate bending constitutive matrices
Ds = cell(1,nMaterials);         % preallocate shear constitutive matrices

for iMaterial = 1:nMaterials
    
    E  = materials(iMaterial,1);
    nu = materials(iMaterial,2);
    Db{iMaterial} = E/(1-nu^2)*[1 nu 0; nu 1 0; 0 0 (1-nu)/2 ];
    Ds{iMaterial} = E/(2*(1+nu))*eye(2);
    
end

%%

K = zeros(nDofsSys,nDofsSys);
M = zeros(nDofsSys,nDofsSys);

for iElem=1:nElems

    elemNodes = elems(iElem,1:nNodesElem); % element nodes
    
    xCoords = nodes(elemNodes,1);  % x-coordinates of nodes
    yCoords = nodes(elemNodes,2);  % y-coordinates of nodes
    
    iMaterial = elems(iElem,end);
    
    thisDb = Db{iMaterial};
    thisDs = Ds{iMaterial};
    
    rho = materials(iMaterial,3);
    
    switch plateType
        
        case 'kirchhoff'
            
            [kElem,mElem] = kirchhoffQ4PlateMatrices(thisDb,rho,points,weights,...
                nNodesElem,nDofsNode,xCoords,yCoords,h);
            
        case 'mindlin'
            
            [kElem,mElem] = mindlinPlateMatrices(thisDb,thisDs,rho,points,weights,...
                nNodesElem,nDofsNode,xCoords,yCoords,h);
            
        case 'kirchhoffPRZ'
            
            lx = max(xCoords) - min(xCoords);
            ly = max(yCoords) - min(yCoords);
            
            E   = materials(iMaterial,1);
            nu  = materials(iMaterial,2);
            [kElem,mElem]=kirchhoffPRZPlateMatrices(lx,ly,E,nu,rho,h);
            
    end
            
    index = findNodeDofs(elemNodes,nDofsNode);
    
    K(index,index) = K(index,index) + kElem;
    M(index,index) = M(index,index) + mElem;

end

% force symmetry
K = triu(K) + triu(K,1).';
M = triu(M) + triu(M,1).';

end

%% testing plate element

function [kElem,mElem]=kirchhoffPRZPlateMatrices(lx,ly,E,nu,rho,h)
%       This function produces rectangular thin plate classical theory finite element 
%       stiffness and mass matrices following the expressions given by
%       Przemieniecki in the book entitled Introduction do Matriz Structural Analysis 
%       noncompatible displacements were used
% lx - element length in x
% ly - element length in y
% E - Young's modulus
% nu - Poisson's ratio
% rho  mass density (volume)
% h - plate thickness

kElem=zeros(12,12);

CT1=E*h^3/(12*(1-nu^2)*lx*ly);
C1=ly/lx;    % beta
C2=C1*C1;  % beta^2
C3=1/C2;   % beta^-2
C4=14-4*nu;
C5=1+4*nu;
C6=1-nu;

kElem(1,1)=4*(C2+C3)+(1/5)*C4;
kElem(2,1)=(2*C3+(1/5)*C5)*ly;
kElem(3,1)=-(2*C2+(1/5)*C5)*lx;
kElem(4,1)=2*(C2-2*C3)-(1/5)*C4;
kElem(5,1)=(2*C3+(1/5)*C6)*ly;
kElem(6,1)=(-C2+(1/5)*C5)*lx;
kElem(7,1)=-2*(C2+C3)+(1/5)*C4;
kElem(8,1)=(C3-(1/5)*C6)*ly;
kElem(9,1)=(-C2+(1/5)*C6)*lx;
kElem(10,1)=-2*(2*C2-C3)-(1/5)*C4;
kElem(11,1)=(C3-(1/5)*C5)*ly;
kElem(12,1)=-(2*C2+(1/5)*C6)*lx;

kElem(2,2)=((4/3)*C3+(4/15)*C6)*ly*ly;
kElem(3,2)=-nu*lx*ly;
kElem(4,2)=-kElem(5,1);
kElem(5,2)=((2/3)*C3-(1/15)*C6)*ly*ly;
kElem(7,2)=-kElem(8,1);
kElem(8,2)=((1/3)*C3+(1/15)*C6)*ly*ly;
kElem(10,2)=kElem(11,1);
kElem(11,2)=((2/3)*C3-(4/15)*C6)*ly*ly;

kElem(3,3)=((4/3)*C2+(4/15)*C6)*lx*lx;
kElem(4,3)=kElem(6,1);
kElem(6,3)=((2/3)*C2-(4/15)*C6)*lx*lx;
kElem(7,3)=-kElem(9,1);
kElem(9,3)=((1/3)*C2+(1/15)*C6)*lx*lx;
kElem(10,3)=-kElem(12,1);
kElem(12,3)=((2/3)*C2-(1/15)*C6)*lx*lx;

kElem(4,4)=kElem(1,1);
kElem(5,4)=-kElem(2,1);
kElem(6,4)=kElem(3,1);
kElem(7,4)=kElem(10,1);
kElem(8,4)=-kElem(10,2);
kElem(9,4)=kElem(12,1);
kElem(10,4)=kElem(7,1);
kElem(11,4)=kElem(7,2);
kElem(12,4)=-kElem(7,3);

kElem(5,5)=kElem(2,2);
kElem(6,5)=-kElem(3,2);
kElem(7,5)=-kElem(11,1);
kElem(8,5)=kElem(11,2);
kElem(10,5)=kElem(8,1);
kElem(11,5)=kElem(8,2);

kElem(6,6)=kElem(3,3);
kElem(7,6)=-kElem(12,1);
kElem(9,6)=kElem(12,3);
kElem(10,6)=-kElem(9,1);
kElem(12,6)=kElem(9,3);

kElem(7,7)=kElem(1,1);
kElem(8,7)=-kElem(2,1);
kElem(9,7)=-kElem(3,1);
kElem(10,7)=kElem(4,1);
kElem(11,7)=kElem(4,2);
kElem(12,7)=-kElem(4,3);

kElem(8,8)=kElem(2,2);
kElem(9,8)=kElem(3,2);
kElem(10,8)=-kElem(4,2);
kElem(11,8)=kElem(5,2);

kElem(9,9)=kElem(3,3);
kElem(10,9)=-kElem(4,3);
kElem(12,9)=kElem(6,3);

kElem(10,10)=kElem(1,1);
kElem(11,10)=kElem(2,1);
kElem(12,10)=-kElem(3,1);

kElem(11,11)=kElem(2,2);
kElem(12,11)=-kElem(3,2);

kElem(12,12)=kElem(3,3);

% for I=1:12                   
%    for J=1:12
%         if J>I
%            kElem(I,J)=kElem(J,I);
%         end
%    end
% end

kElem = kElem + tril(kElem,-1).';

kElem=CT1*kElem;

% Assembles the mass matrix

mElem=zeros(12,12);

CT1=(rho*lx*ly*h)/176400.;

mElem(1,1)=24178.;
mElem(2,1)=3227*ly;
mElem(3,1)=-3227*lx;
mElem(4,1)=8582.;
mElem(5,1)=-1918.*ly;
mElem(6,1)=-1393.*lx;
mElem(7,1)=2758.;
mElem(8,1)=-812.*ly;
mElem(9,1)=812.*lx;
mElem(10,1)=8582.;
mElem(11,1)=1393*ly;
mElem(12,1)=1918*lx;

mElem(2,2)=560.*ly*ly;
mElem(3,2)=-441.*lx*ly;
mElem(4,2)=1918.*ly;
mElem(5,2)=-420.*ly*ly;
mElem(6,2)=-294.*lx*ly;
mElem(7,2)=812.*ly;
mElem(8,2)=-210.*ly*ly;
mElem(9,2)=196.*lx*ly;
mElem(10,2)=1393.*ly;
mElem(11,2)=280.*ly*ly;
mElem(12,2)=294.*lx*ly;

mElem(3,3)=560.*lx*lx;
mElem(4,3)=-1393.*lx;
mElem(5,3)=294.*lx*ly;
mElem(6,3)=280.*lx*lx;
mElem(7,3)=-812.*lx;
mElem(8,3)=196.*lx*ly;
mElem(9,3)=-210.*lx*lx;
mElem(10,3)=-1918.*lx;
mElem(11,3)=-294.*lx*ly;
mElem(12,3)=-420.*lx*lx;

mElem(4,4)=24178;
mElem(5,4)=-3227*ly;
mElem(6,4)=-3227*lx;
mElem(7,4)=8582;
mElem(8,4)=-1393*ly;
mElem(9,4)=1918*lx;
mElem(10,4)=2758;
mElem(11,4)=812*ly;
mElem(12,4)=812*lx;

mElem(5,5)=560*ly*ly;
mElem(6,5)=441*lx*ly;
mElem(7,5)=-1393*ly;
mElem(8,5)=280*ly*ly;
mElem(9,5)=-294*lx*ly;
mElem(10,5)=-812*ly;
mElem(11,5)=-210*ly*ly;
mElem(12,5)=-196*lx*ly;

mElem(6,6)=560*lx*lx;
mElem(7,6)=-1918*lx;
mElem(8,6)=294*lx*ly;
mElem(9,6)=-420*lx*lx;
mElem(10,6)=-812*lx;
mElem(11,6)=-196*lx*ly;
mElem(12,6)=-210*lx*lx;

mElem(7,7)=24178;
mElem(8,7)=-3227*ly;
mElem(9,7)=3227*lx;
mElem(10,7)=8582;
mElem(11,7)=1918*ly;
mElem(12,7)=1393*lx;

mElem(8,8)=560*ly*ly;
mElem(9,8)=-441*lx*ly;
mElem(10,8)=-1918*ly;
mElem(11,8)=-420*ly*ly;
mElem(12,8)=-294*lx*ly;

mElem(9,9)=560*lx*lx;
mElem(10,9)=1393*lx;
mElem(11,9)=294*lx*ly;
mElem(12,9)=280*lx*lx;

mElem(10,10)=24178;
mElem(11,10)=3227*ly;
mElem(12,10)=3227*lx;

mElem(11,11)=560*ly*ly;
mElem(12,11)=441*lx*ly;

mElem(12,12)=560*lx*lx;

% for I=1:12            
%     for J=1:12
%         if J>I
%             mElem(I,J)=mElem(J,I);
%         end
%     end
% end

mElem = mElem + tril(mElem,-1).';

mElem=CT1*mElem;
end

%% element bending and mass calculation for Mindlin plate

function [kElem,mElem] = mindlinPlateMatrices(Db,Ds,rho,points,weights,...
    nNodesElem,nDofsNode,xCoords,yCoords,h)

% [kElem,mElem] = kirchhoffQ4PlateMatrices(D,rho,points,weights,...
%    nNodesElem,nDofsNode,xCoords,yCoords,h)
%
% Computes plate element stiffness and mass matrix considering generic node
% coordinates in the xy plane. Does not assume a rectangular shape.
%
% kElem - element stiffness matrix
% mElem - element consistent mass matrix
% Db - plate bending stiffness
% Ds - plate shear stiffness
% rho - mass density
% points - numerical integration points
% weights - numerial integration weights
% nNodesElem - number of nodes per element
% nDofsNode - number of dofs per node
% xCoords - x-coordinates of nodes
% yCoords - y-coordinates of nodes
% h - plate thickness

nDofsElem = nNodesElem*nDofsNode; % dofs per element

kElem = zeros(nDofsElem,nDofsElem);
mElem = zeros(nDofsElem,nDofsElem);

Mdiag = [ h 0 0;
         0 h^3/12 0;
         0 0 h^3/12 ]; % for consistent mass matrix
kappa = 5/6; % shear correction factor

for iPoint = 1:size(points,1)
    
    xi  = points(iPoint,1);
    eta = points(iPoint,2);
    wxi_weta = weights(iPoint);

    % interpolation and shape
    [N,dN_dxi,dN_deta] = feIsoQ4(xi,eta); % coordinates
    Np = feShape(nNodesElem,nDofsNode,N);
    
    % Jacobian
    J = feJacob2(dN_dxi,dN_deta,xCoords,yCoords);
    detJ = det(J);
    
    % kinematic
    [dN_dx,dN_dy] = feDerivative2(dN_dxi,dN_deta,J);
    [Bb,Bs] = feKinematicMindlin(nNodesElem,N,dN_dx,dN_dy);
    
    % bending element matrix
    kElem = kElem + (h^3/12)*(Bb.'*Db*Bb)*wxi_weta*detJ + ...
        kappa*h*(Bs.'*Ds*Bs)*wxi_weta*detJ;
    
    % consistent mass element matrix
    mElem = mElem + rho*(Np.'*Mdiag*Np)*wxi_weta*detJ;
            
end

end

%% element bending and mass calculation for 4-node quadrangular Kirchhoff plate

function [kElem,mElem] = kirchhoffQ4PlateMatrices(D,rho,points,weights,...
    nNodesElem,nDofsNode,xCoords,yCoords,h)

% [kElem,mElem] = kirchhoffQ4PlateMatrices(D,rho,points,weights,...
%    nNodesElem,nDofsNode,xCoords,yCoords,h)
%
% Computes plate element stiffness and mass matrix considering generic node
% coordinates in the xy plane. Does not assume a rectangular shape.
%
% kElem - element stiffness matrix
% mElem - element consistent mass matrix
% D - plate bending stiffness
% rho - mass density
% points - numerical integration points
% weights - numerial integration weights
% nNodesElem - number of nodes per element
% nDofsNode - number of dofs per node
% xCoords - x-coordinates of nodes
% yCoords - y-coordinates of nodes
% h - plate thickness

nDofsElem = nNodesElem*nDofsNode; % dofs per element

kElem = zeros(nDofsElem,nDofsElem);
mElem = zeros(nDofsElem,nDofsElem);

A = [ evalKirchhoffA(xCoords(1),yCoords(1));
      evalKirchhoffA(xCoords(2),yCoords(2));
      evalKirchhoffA(xCoords(3),yCoords(3));
      evalKirchhoffA(xCoords(4),yCoords(4)) ];

nDofsElem = nNodesElem*nDofsNode;
invA = A\eye(nDofsElem);

for iPoint = 1:size(points,1)
    
    % integration points
    xi  = points(iPoint,1);
    eta = points(iPoint,2);
    wxi_weta = weights(iPoint);
    
    % Jacobian using bilinear
    [Nc,dNc_dxi,dNc_deta] = feIsoQ4(xi,eta);
    J = feJacob2(dNc_dxi,dNc_deta,xCoords,yCoords);
    detJ = det(J);

    % kinematic
    x = Nc*xCoords;
    y = Nc*yCoords;
    B = feKinematicKirchhoff(x,y,invA);
    
    % bending element matrix
    kElem = kElem + (h^3/12)*(B.'*D*B)*wxi_weta*detJ;
    
    % consistent mass element matrix
    N = evalKirchhoffA(x,y);
    Nw = N(1,:)*invA; % w shape function
    mElem = mElem + rho*h*(Nw.'*Nw)*wxi_weta*detJ; % w.^2
    
end

end

%% auxiliary functions

function [N] = feShape(nNodesElem,nDofsNode,Ns)

%% [N] = feShape(nNodesElem,nDofsNode,Ns)
%
% Determine the matrix containing shape functions
%
% ex: [ N1 0 N2 0 N3 0 N4 0 ]
%     [ 0 N1 0 N2 0 N3 0 N4 ]
%
% nNodesElem - number of nodes per element
% nDofsNode - number of dofs per node
% Ns - shape function (for each node)
% dofs: uz psi_x psi_y, ux uy uz...

N = zeros(nDofsNode,nDofsNode*nNodesElem);

for i=1:nDofsNode
    N(i,i:nDofsNode:end) = Ns;
end

end

function [Bb,Bs] = feKinematicMindlin(nNodesElem,N,dN_dx,dN_dy)

%% [Bb,Bs] = feKinematicMindlin(nnel,dNd_x,dN_dy)
%
% Determine the kinematic matrix expression for bending and shear matrices
% of the Mindlin plate
%
% nNodesElem - number of nodes per element
% dN_dx - derivatives of shape functions with respect to x
% dN_dy - derivatives of shape functions with respect to y

Bb = zeros(3,3*nNodesElem);
Bs = zeros(2,3*nNodesElem);

nDofsNode = 3;

Bb(1,2:nDofsNode:end) = dN_dx;
Bb(2,3:nDofsNode:end) = dN_dy;
Bb(3,2:nDofsNode:end) = dN_dy;
Bb(3,3:nDofsNode:end) = dN_dx;

Bs(1,1:nDofsNode:end) = dN_dx;
Bs(1,2:nDofsNode:end) = -N;
Bs(2,1:nDofsNode:end) = dN_dy;
Bs(2,3:nDofsNode:end) = -N;

end

function [B] = feKinematicKirchhoff(x,y,invA)

%% [B] = feKinematicKirchhoff(x,y,invA)
%
% Kirchhoff plate element for use with inverse of matrix evaluated using
% generalized displacements.
% 
% {k} = {d2w/dx2 d2w/dy2 2 d2w/dxdy } = [B]{ w1 w1,x w1,y ... }
  
B = [0, 0, 0, 2, 0, 0, 6*x, 2*y, 0, 0, 6*x*y, 0;
     0, 0, 0, 0, 0, 2, 0, 0, 2*x, 6*y, 0, 6*x*y;
     0, 0, 0, 0, 2, 0, 0, 4*x, 4*y, 0, 6*x^2, 6*y^2]*invA;

end

function [A3] = evalKirchhoffA(x,y)

%% [A3] = evalKirchhoffA(x,y)
%
% Evaluate matrix containing polynomial with generalized degrees-of-freedom
% to determine kinematic matrix.

A = [ 1 x y x^2 x*y y^2 x^3 x^2*y x*y^2 y^3 x^3*y x*y^3 ];
dA_dx = [0, 1, 0, 2*x, y, 0, 3*x^2, 2*x*y, y^2, 0, 3*x^2*y, y^3];
dA_dy = [0, 0, 1, 0, x, 2*y, 0, x^2, 2*x*y, 3*y^2, x^3, 3*x*y^2];

A3 = [ A; dA_dx; dA_dy ];

end

function [J] = feJacob2(dN_dxi,dN_deta,xCoords,yCoords)

%% [J] = feJacob2(nnel,dN_dqsi,dhds,xcoord,ycoord)
%
% Determine de Jacobian for two-dimensional mapping
%
% J - Jacobian for two-dimension
% dN_dxi - derivative of shape functions wrt natural coordinate xi
% dN_deta - derivative of shape functions wrt natural coordinate eta
% xCoords - x axis coordinate values of nodes
% yCoords - y axis coordinate values of nodes

J = [ dN_dxi; dN_deta ]*[ xCoords yCoords ];

end

function [N,dN_dxi,dN_deta] = feIsoQ4(xi,eta)

%% [N,dN_dxi,dN_deta] = feisoq4(xi,eta)
%
% Compute isoparametric four-node quadrilateral shape functions and
% their derivatives at the selected integration point in terms of
% natural coordinates (Kwon), for independent displacement fields
%
% N - shape functions for four-node element
% dN_dxi - derivatives of the shape function wrt xi
% dN_deta - derivatives of the shape function wrt eta
% xi - xi coordinate value of the selected point
% eta - eta coordinate value of the selected point

% coordinates: 1 (-1,-1), 2 (1,-1), 3 (1,1), 4 (-1,1)

N = 1/4*[ (1-xi)*(1-eta) (1+xi)*(1-eta) (1+xi)*(1+eta) (1-xi)*(1+eta) ];
dN_dxi = 1/4*[ -(1-eta) (1-eta) (1+eta) -(1+eta) ];
dN_deta = 1/4*[ -(1-xi) -(1+xi) (1+xi) (1-xi) ];

end

function [points,weights] = feGaussQuadrature2(nGL,type)

%% [points,weights] = feGaussQuad2(nGL,type)
%
% Determine integration points and weighting coefficients of
% Gauss-Legendre quadrature for two-dimensional integration.
%
% points - matrix containing integration points (xi, eta)
% weights - vector containing weight coefficients (wxi*weta)
% nGL - number of integration points in one axis
% type - 'quad' or 'tria'

switch type
    
    case 'quad' % use basic 1D points

        [point1,weight1] = feGaussQuadrature1(nGL);
        
        [point_xi,point_eta] = meshgrid(point1);
        [weight_xi,weight_eta] = meshgrid(weight1);
        
        points = [ point_xi(:) point_eta(:) ];
        weights = weight_xi(:).*weight_eta(:);
        
    case 'tria'
        
        switch nGL

            case 1
                points = [ 1/3 1/3 ];
                weights = 1.0;

            case 2 % quadratic
                points = [ 2/3 1/6;
                           1/6 1/6;
                           1/6 2/3 ];
                weights = [ 1/3; 1/3; 1/3 ];

            case 3
                points = [ 1/2 0;
                           0 1/2;
                           1/2 1/2 ];
                weights = [ 1/3; 1/3; 1/3 ];

            case 4 % cubic
                points = [ 1/3 1/3;
                           3/5 1/5;
                           1/5 1/5;
                           1/5 3/5 ];
                weights = [ -27/48; 25/48; 25/48; 25/48 ];

        end
        
end

end

function [point1,weight1] = feGaussQuadrature1(ngl)

%% [point1,weight1] = feGaussQuadrature1(ngl)
%
% Determine integration points and weighting coefficients of
% Gauss-Legendre quadrature for one-dimensional integration (Kwon)
%
% ngl - number of integration points
% point1 - vector containing integration points
% weight1 - vector containing weight coefficients

switch ngl
    
    case 1
        
        point1 = 0;
        weight1 = 2;
        
    case 2
        
        a = sqrt(1/3);
        point1 = [-a a];
        weight1 = [1 1];
        
    case 3
        
        a = sqrt(3/5);
        point1 = [-a 0 a];
        weight1 = [5/9 8/9 5/9];
        
    case 4
        
        a = sqrt(3/7+2/7*sqrt(6/5));
        b = sqrt(3/7-2/7*sqrt(6/5));
        
        c = (18-sqrt(30))/36;
        d = (18+sqrt(30))/36;
        
        point1 = [-a -b b a];
        weight1 = [c d d c];
        
    case 5
        
        a = 1/3*sqrt(5+2*sqrt(10/7));
        b = 1/3*sqrt(5-2*sqrt(10/7));
        
        c = (322-13*sqrt(70))/900;
        d = (322+13*sqrt(70))/900;
        
        point1 = [ -a -b 0 b a ];
        weight1 = [ c d 128/225 d c ];
        
    otherwise
        
        point1 = [];
        weight1 = [];
        disp('no defined Gauss-Legendre integration rule.')
        
end

point1 = point1.';
weight1 = weight1.';

end

function [nodesDofs] = findNodeDofs(nodeList,nDofsNode)

% Computes global dofs associated to a list of nodes
%
% nodesDofs - vector of global DOFs associated with nodes
% nodeList - list of node numbers whose DOFs must be determined
% nDofsNode - nuber of dofs per node

%%

nNodes = length(nodeList);
nodesDofs = zeros(1,nNodes*nDofsNode);

for iNode=1:nNodes % nodes loop
    
    node = nodeList(iNode);
    nodesDofs ( (iNode-1)*nDofsNode+1 : iNode*nDofsNode ) = (node-1)*nDofsNode+1 : node*nDofsNode ;

end

end

function [dN_dx,dN_dy] = feDerivative2(dN_dxi,dN_deta,J)

%% [dN_dx,dN_dy] = feDerivative2(dN_dxi,dN_deta,J)
%
% Determine derivative of 2D isoparametric shape functions with respect
% to physical coordinate system.
%
% dN_dx - derivative of shape function wrt physical coordinate x
% dN_dy - derivative of shape function wrt physical coordinate y
% dN_dxi - derivative of shape function wrt natural coordinate xi
% dN_deta - derivative of shape function wrt natural coordinate eta
% J - Jacobian 2D matrix

%(2 x nNodesEleme) = (2x2) x (2 x nNodesEleme)
dN_dxy = J\[ dN_dxi; dN_deta ];

dN_dx = dN_dxy(1,:); % row of x-derivatives
dN_dy = dN_dxy(2,:); % row of y-derivatives

end
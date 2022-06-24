function [K,M] = feKirchhoffPlateMatrix(nodes,elem,ndof,E,nu,rho,h)
%
%% [K,M] = feKirchhoffPlateMatrix(nodes,elem,ndof,E,nu,t,rho)
%
% assembles stiffness and mass matrix
%
% K - stiffness matrix
% M - mass matrix
% nodes - nodal coordinate values
% elem - nodal connectivity
% ndof - numbers of dof per node
% E - Young's modulus
% nu - Poisson ratio
% rho - specific mass
% h - thickness

% number of dofs per node
% Kirchoff (thin): w, dw/dx, dw/dy

% choice of integration points
nglb = 2; % complete (bending)
nglm = 2; % complete (mass)

nglxb = nglb; nglyb = nglb;
nglxm = nglm; nglym = nglm;

[pointb,weightb] = feGaussQuad2(nglxb,nglyb);
[pointm,weightm] = feGaussQuad2(nglxm,nglym);

% dimensions
nnode = size(nodes,1); % total number of nodes in the system
nel = size(elem,1);    % number of elements
nnel = size(elem,2);   % number of nodes per element
edof = nnel*ndof;      % dofs per element
sdof = nnode*ndof;     % total system dofs

K = zeros(sdof,sdof);
M = zeros(sdof,sdof);

%inserting matrices with zeros in sparse form
% K = sparse(sdof,sdof);
% M = sparse(sdof,sdof);
% 
% K = spalloc(sdof,sdof,1);
% M = spalloc(sdof,sdof,1);


D = E/(1-nu^2)*[1 nu 0;
                nu 1 0;
                0 0 (1-nu)/2 ];

for iel=1:nel

    nd = elem(iel,1:nnel);   % element nodes   
    xcoord = nodes(nd,1);    % x-coordinates of nodes
    ycoord = nodes(nd,2);    % y-coordinates of nodes

	k = plateBending(D,nglxb,nglyb,pointb,weightb,nnel,edof,xcoord,ycoord,h);     % numerical integration for bending   
	m = plateMass(nglxm,nglym,pointm,weightm,rho,nnel,ndof,edof,xcoord,ycoord,h); % numerical integration for mass
    
    index = findNodeDofs(nd,ndof);

    K(index,index) = K(index,index) + k;
    M(index,index) = M(index,index) + m;

end

end

%% bending calculation

function [k] = plateBending(D,nglxb,nglyb,pointb,weightb,nnel,edof,xcoord,ycoord,h)

k = zeros(edof,edof);

for intx = 1:nglxb
    
    xi = pointb(intx,1);
    wxi = weightb(intx,1);

    for inty = 1:nglyb

        eta = pointb(inty,2);
        weta = weightb(inty,2);
        
        [~,dNc_dxi,dNc_deta] = feIsoQ4(xi,eta);
        
        J = feJacob2(nnel,dNc_dxi,dNc_deta,xcoord,ycoord);
        detJ = det(J);
        
        B = feKinematicBend_Kirchhoff_Parametric(xi,eta,J);               
        
        % bending element matrix
        k = k + B.'*D*B*wxi*weta*detJ;
            
    end
end

k = h^3/12*k;

end

%% mass calculation 

function [m] = plateMass(nglxm,nglym,pointm,weightm,rho,nnel,ndof,edof,xcoord,ycoord,h)

m = zeros(edof,edof);

Mdiag = [ h 0 0 ; 0 eps 0 ; 0 0 eps ];

for intx = 1:nglxm
    
    xi = pointm(intx,1);
    wxi = weightm(intx,1);
    
    for inty = 1:nglym
        
        eta = pointm(inty,2);
        weta = weightm(inty,2);
        
        [Nc,dNc_dxi,dNc_deta] = feIsoQ4(xi,eta);
                        
        J = feJacob2(nnel,dNc_dxi,dNc_deta,xcoord,ycoord);        
        detJ = det(J);

        Np = feShape(nnel,ndof,Nc);
                           
        % consistent mass element matrix                    
        m = m + Np.'*Mdiag*Np*wxi*weta*detJ;
                    
    end
end

m = rho*m;
    
end

%% 

function [N,dN_dxi,dN_deta] = feIsoQ4(xi,eta)

%% [N,dN_dxi,dN_deta] = feisoq4(xi,eta)
%
% Compute isoparametric four-node quadrilateral shape functions and
% their derivatives at the selected integration point in terms of
% natural coordinates (Kwon)
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

function [Bb] = feKinematicBend_Kirchhoff_Parametric(xi,eta,J)

%% [Bb] = feKinematicBend_Kirchhoff_Parametric(xi,eta,J)
%
% Determine the kinematic matrix expression relating bending curvatures
% to rotations and displacements for Kirchhoff plates
%
% xi - xi coordinate
% eta - eta coordinate
% J - jacobian matrix

invJ = J\eye(2);
Bb = transfdiff2w(invJ)*diff2N(xi,eta)*transfd(J);

end

function d2w = transfdiff2w(invJ)

% transformation between second derivatives in xi,eta and x,y spaces

T11 = invJ(1,1);
T12 = invJ(1,2);
T21 = invJ(2,1);
T22 = invJ(2,2);

d2w = [ T11^2 T12^2 2*T11*T12;
       T21^2 T22^2 2*T21*T22;
       2*T11*T21 2*T12*T22 2*(T11*T22+T12*T21) ];

end

function d2N = diff2N(xi,eta)

% calculates second derivative of shape functions wrt xi, eta

d2N = ...
[        (3*xi)/4 - (3*eta*xi)/4, eta/4 + (3*xi)/4 - (3*eta*xi)/4 - 1/4,                                     0,        (3*eta*xi)/4 - (3*xi)/4, (3*xi)/4 - eta/4 - (3*eta*xi)/4 + 1/4,                                     0,      - (3*xi)/4 - (3*eta*xi)/4, eta/4 + (3*xi)/4 + (3*eta*xi)/4 + 1/4,                                     0,        (3*xi)/4 + (3*eta*xi)/4, (3*xi)/4 - eta/4 + (3*eta*xi)/4 - 1/4,                                     0;
       (3*eta)/4 - (3*eta*xi)/4,                                     0, (3*eta)/4 + xi/4 - (3*eta*xi)/4 - 1/4,       (3*eta)/4 + (3*eta*xi)/4,                                     0, (3*eta)/4 - xi/4 + (3*eta*xi)/4 - 1/4,     - (3*eta)/4 - (3*eta*xi)/4,                                     0, (3*eta)/4 + xi/4 + (3*eta*xi)/4 + 1/4,       (3*eta*xi)/4 - (3*eta)/4,                                     0, (3*eta)/4 - xi/4 - (3*eta*xi)/4 + 1/4;
 1/2 - (3*xi^2)/8 - (3*eta^2)/8,               xi/4 - (3*xi^2)/8 + 1/8,             eta/4 - (3*eta^2)/8 + 1/8, (3*eta^2)/8 + (3*xi^2)/8 - 1/2,               1/8 - (3*xi^2)/8 - xi/4,             (3*eta^2)/8 - eta/4 - 1/8, 1/2 - (3*xi^2)/8 - (3*eta^2)/8,               (3*xi^2)/8 + xi/4 - 1/8,             (3*eta^2)/8 + eta/4 - 1/8, (3*eta^2)/8 + (3*xi^2)/8 - 1/2,               (3*xi^2)/8 - xi/4 - 1/8,             1/8 - (3*eta^2)/8 - eta/4];

end

function transfd = transfd(J)

% tansformation between nodal coordinates in xi,eta and x,y space

transfd = eye(12);

transfd(2:3,2:3) = J;
transfd(5:6,5:6) = J;
transfd(8:9,8:9) = J;
transfd(11:12,11:12) = J;

end

function [J] = feJacob2(nnel,dN_dxi,dN_deta,xcoord,ycoord)

%% [J] = feJacob2(nnel,dN_dqsi,dhds,xcoord,ycoord)
%
% Determine de Jacobian for two-dimensional mapping
%
% J - Jacobian for two-dimension
% nnel - number of nodes per element
% dN_dxi - derivative of shape functions wrt natural coordinate xi
% dN_deta - derivative of shape functions wrt natural coordinate eta
% xcoord - x axis coordinate values of nodes
% ycoord - y axis coordinate values of nodes

J = zeros(2,2);

for i=1:nnel
    J = J + [ dN_dxi(i); dN_deta(i) ]*[xcoord(i) ycoord(i)];
end

%J = [ xcoord.'; ycoord.' ] * [ dN_dxi.' dN_deta.' ];

end

function [N] = feShape(nnel,ndof,Ns)

%% [N] = feShape(nnel,ndof,Ns)
%
% Determine the matrix containing shape functions
%
% ex: [ N1 0 N2 0 N3 0 N4 0 ]
%     [ 0 N1 0 N2 0 N3 0 N4 ]
%
% nnel - number of nodes per element
% ndof - number of dofs per node
% Ns - shape function (for each node)
% dofs: uz psi_x psi_y, ux uy uz...

if( length(Ns) == nnel )

    N = zeros(ndof,ndof*nnel);
    I = eye(ndof);

    for i=1:nnel

        il = (i-1)*ndof+1;
        iu = (i-1)*ndof+ndof;

        N(1:ndof,il:iu)=Ns(i)*I;

    end
    
else
    
	disp('Number of shape functions is not correct.')
    N = [];
    return

end

end

function [point1,weight1] = feGaussQuad1(ngl)

%% [point1,weight1] = feGaussQuad1(ngl)
%
% Determine integration points and weighting coefficients of
% Gauss-Legendre quadrature for one-dimensional integration (Kwon)
%
% ngl - number of integration points
% point1 - vector containing integration points
% weight1 - vector containing weight coefficients

switch ngl
    
    case 1
        
        point1 = [0];
        weight1 = [2];
        
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
        disp('sem regra definida.')
        
end

point1 = point1';
weight1 = weight1';

end

function [point2,weight2] = feGaussQuad2(ngl_xi,ngl_eta)

%% [point2,weight2] = feGaussQuad2(nglx,ngly)
%
% Determine integration points and weighting coefficients of
% Gauss-Legendre quadrature for two-dimensional integration (Kwon)
%
% ngl_xi - number of integration points in the xi-axis
% ngl_eta - number of integration points in the eta-axis
% point2 - vectors containing integration points
% weight2 - vectors containing weight coefficients

ngl = max(ngl_xi,ngl_eta);

point2 = zeros(ngl,2);
weight2 = zeros(ngl,2);

[point_xi,weight_xi] = feGaussQuad1(ngl_xi);
[point_eta,weight_eta] = feGaussQuad1(ngl_eta);

% [ xi eta ]

for int_xi=1:ngl_xi
    point2(int_xi,1) = point_xi(int_xi);
    weight2(int_xi,1) = weight_xi(int_xi);
end

for int_eta=1:ngl_eta
    point2(int_eta,2) = point_eta(int_eta);
    weight2(int_eta,2) = weight_eta(int_eta);
end

end
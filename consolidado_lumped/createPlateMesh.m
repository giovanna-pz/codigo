function [nodes,elems] = createPlateMesh(varargin)

%% [nodes,elems] = createPlateMesh(Lx,[Ly],nElemsX,[nElemsY])
%
% Creates a plate mesh of size Lx x Ly, with prescribed number of elements
% for each direction.
% 
% nodes - node coordinates
% elems - nodal connectivity
% Lx - x-direction length
% Ly - y-direction length (optional)
% nElemsX - number of elements in the x-direction
% nElemsY - number of elements in the y-direction (optional)

%% separate between square and rectangular cases

switch nargin
    
    case 2 %square
        
        Lx = varargin{1};
        Ly = Lx;
        
        nElemsX = varargin{2};
        nElemsY = nElemsX;
        
    case 4 %rectangular
        
        Lx = varargin{1};
        Ly = varargin{2};
        
        nElemsX = varargin{3};
        nElemsY = varargin{4};
        
end

%% node coordinates

gridX = -Lx/2:Lx/nElemsX:Lx/2;
gridY = -Ly/2:Ly/nElemsY:Ly/2;

nGridX = length(gridX);
nGridY = length(gridY);

nodes = zeros(nGridX*nGridY,3);

for iy=1:nGridY
    y = gridY(iy);
    for ix=1:nGridX
        x = gridX(ix);
        nodes( (iy-1)*nGridX + ix, : ) = [ x y 0 ];
    end
end

%% nodal connectivty

elems = zeros((nGridX-1)*(nGridY-1),5); % [ node1 node2 node3 node4 material ]

for iy=1:nGridY-1
    for ix=1:nGridX-1
        
        % [ 4 3 ]
        % [ 1 2 ]

        n1 = (iy-1)*nGridX + ix;
        n2 = (iy-1)*nGridX + ix+1;
        n3 = iy*nGridX + ix+1;
        n4 = iy*nGridX + ix;
        
        iElem = (iy-1)*(nGridX-1) + ix;
        
        elems(iElem,:) = [ n1 n2 n3 n4 1 ];
       
        
    end
end

end
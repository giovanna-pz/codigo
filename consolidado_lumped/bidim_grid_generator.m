function [nodes,elem,nnode,nelem,Ngrid]= bidim_grid_generator(L,elem_side)
% This function generates a equally space mesh in a 2D square structure,
% assuming that the size of elements in x and y direction are the same

% Author: Vinicius Dal Poggeto

%Input:
% L: Size of strcuture
% elem_size: size of one finite element

%Output
%nodes: matrix containing the nodes coordinates;
%elem: m,atrix of the incidence of the nodes in the elements;
% nnode: number of total nodes;
%nelem: number of total elements.


%% grid

Nel_side = L/elem_side;

%Nx = Nel_side + 1;

Lgrid = -L/2:L/Nel_side:L/2;
Ngrid = length(Lgrid);
nodes = zeros(Ngrid^2,2);

for iry=1:Ngrid
    y = Lgrid(iry);
    for irx=1:Ngrid
        x = Lgrid(irx);
        nodes( (iry-1)*Ngrid + irx, : ) = [ x y ];
    end
end

%% nodal connectivity

elem = zeros((Ngrid-1)^2,4);

iel=1;

for iy=1:Ngrid-1
    for ix=1:Ngrid-1
        
        % [ 4 3 ]
        % [ 1 2 ]
        
        n1 = (iy-1)*Ngrid + ix;
        n2 = (iy-1)*Ngrid + ix+1;
        n3 = iy*Ngrid + ix+1;
        n4 = iy*Ngrid + ix;
        
        elem(iel,:) = [ n1 n2 n3 n4 ]; % material 1
       
        iel = iel+1;
        
    end
    
end

nnode = length(nodes(:,1));
nelem = length(elem(:,1));

end
clear
% close all
clc

%% This program calculates the displacements and the STL of a metamaterial plate 

% Infinite plate with periodical resonators using the WFE method and the assumption that in a periodical
%material space harmonics are created.

%Case of 4 resonator in the unit cell

%Reference: Y.Yang et al., Journal of Sound and Vibration, 2019, Vibroacoustic analysis of periodic structures using a wave and
% Reference: Fahy and Gardonio, Sound and structural vibration
%finite element\ method

%The formulation used in this program is an adaptation of Y.Yang et al.
%using the concepts discussed in the book of Fahy and Gardonio

%Example from section 9.1 of the paper

% The plate formulation used here is the Kirchoff plate using a Q4 element
% or  Mindlin plate using a Q4 element
% bilinear element

format long

tic

%% Part 1 - Parameters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% BEGGINING OF INPUT DATA %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Paramaters changed for analysis

% FREQUENCY VECTOR DATA

%same frequency vector as reference
fmin = 50; %[Hz]
fmax = 5000; %[Hz]

freq = logspace(log10(fmin),log10(fmax),1e4);
omega = 2*pi*freq; %[rad/s]
nfreq = length(freq);

% Number of elements
%minimum 4 elements
%other quantitites, multiple of four
nel_x =8;
nel_y =8;

% Other fixed parameters

subprogram_inputdata_onlykx

% subprogram_inputdata_onlyky

% subprogram_inputdata_kxandky

%%%%%%%%%%%%%%%%%%%%%%%%
%%% END OF INPUT DATA %%
%%%%%%%%%%%%%%%%%%%%%%%%

% end Part1

%% Part 2 - Finite Element Model

% %% Testing the meshing size - optional
% % lmin_wtharm: structural wavenumbers without harmonics
% % lmin_h:structural wavenumbers with harmonics
%[lmin_wtharm,lmin_h_x,lmin_h_y] = minimum_l_fem(E,rho,nu,fmax,el_side_x,el_side_y,m,n,theta,phi,c0);

% Generating FE mesh
[nodes,elem,nnode,nelem,Ngrid]= bidim_grid_generator(Lx,el_side_x);
GDof =nnode*dof; % total number of dofs in strucuture 

% Plot mesh
plotMesh(nodes,elem)

% Mass and stiffness matrices - substructure plate
% Choose for desired model
kircchoff_model
%mindlin_model

% Updating K and M with resonators

%center node, where the resonator is attached
%[~,node_res_first]=min(sqrt(nodes(:,1).^2)+nodes(:,2).^2);


%Nodes at 1/4 of the unit cell
[node_res1, ~] = find(nodes(:,1)==-Lx/4 & nodes(:,2)==-Lx/4);
[node_res2, ~] = find(nodes(:,1)==Lx/4 & nodes(:,2)==-Lx/4);
[node_res3, ~] = find(nodes(:,1)==-Lx/4 & nodes(:,2)==Lx/4);
[node_res4, ~] = find(nodes(:,1)==Lx/4 & nodes(:,2)==Lx/4);


node_res_first =[node_res1,node_res2,node_res3,node_res4];

[K_new, M_new, numberRes] = K_M_resonators(KG,MG,node_res_first,dof,GDof,kr,mr);

%Subprogram for partitioning boundary and internal dofs
subprogram_partitioningdofs
% end Part 2

%% Part 3 - Natural Frequencies and modes of Vibration - OPTIONAL
% Calculating natural frequencies and modes of vibration - structural only 
%  calculate_naturalfreqs_modes 
% end Part 3

%% Part 4 - DISPERSION ANALYSIS %% OPTIONAL
% Calculate dispersion diagram and plot wave modes  
calc_disp_plot_modes(Lx,Ly,K_new,M_new,dofs_i,dofs_b,dof,Ngrid,el_side_x,el_side_y,nnode,numberRes)
% end Part 4




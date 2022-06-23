%% INPUT DATA
clear all;
close all;
clc
%% Geometry and materials

% lx: element size in x direction
% ly: element size in y direction
% h: plate thickness
% E: modulus of elasticity 
% nu: Poisson
% eta: loss factor

%Data for a square plate
L= 0.3; %[m]
nel_lat=20; %number of element in one lateral of the plate
lx = L/nel_lat;
ly=lx;
h=0.001; 

% Material parameters (for aluminum plate)

E=7e10;
rho=2700;
nu=0.3;
%eta=0.005;

%% Resonator data
A = 4*lx*ly; % area of a unitary cell
ms = h*rho*A; % [kg] 

%Resonator 1 - node 5
fr = 200; %[Hz]
mr = 0.5*ms; % [kg]

%Nodes where the resonators will be placed in the first bottom line
node_res = [23,25,27,29,31,33,35,37,39,41]; 

%% Boundary conditions
%prescribedNodes = [] % No B.C.
%% Frequencies
fmin = 10;
df = 0.5;
fmax = 1000;

%% Nodes of interest

node_force = 190; %node where the force is applied.
node_obs = 336; % response node

%% Force amplitude
F0 =1;

%% Saving variables
filename = 'input_data.mat';
save(filename)
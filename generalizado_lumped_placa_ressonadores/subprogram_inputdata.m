%DATA SUBSTRATE
E = 7e10; % [Pa]
rho = 2700; %  density[kg/m^3]
h = 1e-3; % thickness [m]
nu = 0.3; % poisson
eta = 0.01; % loss factor

E = E*(1+1j*eta); % structural damping added

% DATA FOR FLUID LAYERS
rho0 = 1.21; % air density[kg/m^3]
c0 = 340; % speed of sound  [m/s]

%INCIDENT WAVE
P_inc=1; % Amplitude of incidence wave

% FINITE ELEMENT MODEL OF THE PERIODIC CELL
% type of finite element: plate with 4 nodes and 3 dofs per node
dof=3;
nnos_el=4;

% Length of one finite element
el_side_x = Lx/nel_x;
el_side_y = Ly/nel_y;

% Area of a single finite element
A_element =el_side_x*el_side_y;
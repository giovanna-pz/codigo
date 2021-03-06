%DATA SUBSTRATE
E = 7e10; % [Pa]
rho = 2700; %  density[kg/m^3]
h = 1e-3; % thickness [m]
nu = 0.3;
%Damp 0
% eta = 0;
% %Damp1
% eta = 0.001;
% %Damp2
 eta = 0.01;

E = E*(1+1j*eta); % structural damping added


% DATA FOR FLUID LAYERS
rho0 = 1.21; % air density[kg/m^3]
c0 = 340; % speed of sound  [m/s]

% DATA FOR INCIDENCE  WAVE 
theta = deg2rad(50);
%phi=0, so tem contribuicao em kx
phi = deg2rad(0);
%phi=0;
P_inc=1; % Ampltude of incidence wave

% DATA FOR RESONATORS - 4 resonators case

%fr1 = 250 Hz
kr1 = 1.665495742683829e+04;
mr1= 0.027/4;

%fr2 =300 Hz
kr2 = 9.593e4; %[N/m]
mr2 = 0.027/4; %[kg]

%fr3 =350 Hz
kr3 = 3.264371655660304e+04;
mr3 = 0.027/4;

%fr4 =400 Hz
kr4 = 4.263669101270603e+04;
mr4 = 0.027/4;

kr = [kr1,kr2,kr3,kr4];
mr = [mr1,mr2,mr3,mr4];

% FINITE ELEMENT MODEL OF THE PERIODIC CELL
% type of finite element: plate with 4 nodes and 3 dofs per node
dof=3;
nnos_el=4;

%Length of structure
Lx = 0.1; %[m]
Ly = 0.1; %[m]


% Length of one finite element
el_side_x = Lx/nel_x;
el_side_y = Ly/nel_y;

% SPACE HARMONICS

%index of harmonics
m= -6:1:6;
n=-6:1:6;
%
[MM,NN]=meshgrid(m,n);

nm = length(m);
nn = length(n);

%NUMERICAL INTEGRATION
% Gauss points for numerical integration
% 2 points integration
csi_aux = [-1/sqrt(3),1/sqrt(3)];
wcsi_aux = [1 , 1];
eta_aux = [-1/sqrt(3),1/sqrt(3)];
weta_aux = [1, 1];


% 3 points integration
% Best results so far!
% csi_aux =  [-0.77459667 0           0.77459667];
% wcsi_aux = [ 0.55555556 0.88888889  0.55555556];
% eta_aux = csi_aux;
% weta_aux = wcsi_aux;

% % 4 points integration
% csi_aux =  [-0.86113631 -0.33998104 0.33998104 0.86113631 ];
% wcsi_aux = [ 0.3478548   0.6521452  0.6521452  0.3478548 ];
% eta_aux = csi_aux;
% weta_aux = wcsi_aux;

% % 5 points integration
% csi_aux =  [-0.9061798459 -0.5384693101  0          0.5384693101  0.9061798459  ];
% wcsi_aux = [ 0.2369269     0.4786287     0.5688889  0.4786287     0.2369269     ];
% eta_aux = csi_aux;
% weta_aux = wcsi_aux;


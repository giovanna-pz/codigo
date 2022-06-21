
[~,node_interest]=min(sqrt(nodes(:,1).^2)+nodes(:,2).^2);
node_force = node_interest; %node where the force is applied.
node_obs = node_interest; % node where the response is measured

% Force vector 
F = zeros(GDof,1);
input = (node_force-1)*dof+1;
F0 = 1000; %[N]
F(input)= -F0; % force in z direction
% Initializating variables for FRF 
u_z = zeros(1,nfreq);  % Dispalcement at the observed node
u_f = zeros(GDof,nfreq); %Displacements/rotations for all frequencies
output = (node_obs-1)*dof+1;

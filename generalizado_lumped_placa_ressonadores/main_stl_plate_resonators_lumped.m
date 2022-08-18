clear
close all
clc

%% This program calculates the displacements and the STL of a plate with
%attached resonators using the 2D WFE method
% The plate formulation used here is the Kirchoff plate using a Q4 element
% or bilinear element

format long

tic

%% Part 1 - Parameters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% BEGGINING OF INPUT DATA %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Paramaters changed for analysis

% FREQUENCY VECTOR DATA
fmin = 50; %[Hz]
fmax = 6000; %[Hz]

freq=logspace(log10(fmin),log10(fmax),1e3);
omega = 2*pi*freq; %[rad/s]
nfreq = length(freq);

% Number of finite elements in x and y directions
nel_x = 8;
nel_y =4;

% nel_x = 16;
% nel_y =16;

%Number of harmonics
N = 5;
m = -N:1:N;
n = -N:1:N;
[MM,NN]=meshgrid(m,n);

nm = length(m);
nn = length(n);

%Length of substructure or periodical cell
Lx = 0.06; %[m]
Ly = 0.03; %[m]

% Lx =0.1;
% Ly = 0.1; %[m]


% Other fixed parameters
subprogram_inputdata

%Oblique incidence wave - angles of incidence

%only_kx
theta = 50;
phi =0;

%only_ky
% phi = 50;
% theta =90;
        
%kx_and_ky
% phi = 50;
% theta =30;

theta =deg2rad(theta);
phi =deg2rad(phi);

%Resonator(s) data
A_cell = Lx*Ly; % area of the periodic cell [m^2]
m_cell = h*rho*A_cell; % [kg] 

% 1 resonator case
%fr =300; %natural frequencie of resonator [Hz]
%mr = 0.2*m_cell;
%rr=1;

% 2 resonators case
fr =[200,250]; %natural frequencie(s) of resonator(s) [Hz]
mr = [0.2*m_cell,0.25*m_cell]; % masse(s) of resonator(s)
rr=2;

%%%%%%%%%%%%%%%%%%%%%%%%
%%% END OF INPUT DATA %%
%%%%%%%%%%%%%%%%%%%%%%%%

% end Part1

%% Part 2 - Finite Element Model

% %% Testing the meshing size - optional
% % lmin_wtharm: structural wavenumbers without harmonics
% % lmin_h:structural wavenumbers with harmonics
[lmin_wtharm,lmin_h_x,lmin_h_y] = minimum_l_fem(E,rho,nu,h,fmax,el_side_x,el_side_y,m,n,theta,phi,c0);

%Mesh generator for square or rectangular periodic cells
[nodes,elem, NgridX, NgridY,nnode,nelem] = createPlateMesh(Lx,Ly,nel_x,nel_y);

GDof =nnode*dof; % total number of dofs in strucuture 

% Plot mesh
% plotMesh(nodes,elem)

% Mass and stiffness matrices - substructure plate
% Choose for desired model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kircchoff_model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Attached resonator(s) - Mass - spring system

kr = 4*(pi^2)*(fr.^2).*mr; % stiffness of resonator [N/m]

% Updating K and M with resonators


switch rr
    case 1

    % 1 resonator case
    %center of the square unit cell
    [~,node_res_first]=min(sqrt(nodes(:,1).^2)+nodes(:,2).^2);
    case 2 
    % 2 resonators case
    [node_res1, ~] = find(nodes(:,1)==-Lx/4 & nodes(:,2)==0);
    [node_res2, ~] = find(nodes(:,1)==Lx/4 & nodes(:,2)==0);
    node_res_first =[node_res1,node_res2];

end

[K_new, M_new, numberRes] = K_M_resonators(KG,MG,node_res_first,dof,GDof,kr,mr);


%Subprogram for partitioning boundary and internal dofs
subprogram_partitioningdofs

% end Part 2

%% Part 3 - Natural Frequencies and modes of Vibration - OPTIONAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating natural frequencies and modes of vibration - structural only 
%  calculate_naturalfreqs_modes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end Part 3

%% Part 4 - DISPERSION ANALYSIS %% OPTIONAL
% Calculate dispersion diagram and plot wave modes  
 calc_disp_plot_modes(Lx,Ly,K_new,M_new,dofs_i,dofs_b,dof,NgridX,NgridY,el_side_x,el_side_y,theta,fmax,numberRes)
% end Part 4

%% Part 5 - FRF calculation, material in vaccum
% subprogram_frf_for_validation_optional
% end Part 5

% Definitions for calculation of the wavenumbers
% Equation 8

% index for summation of harmonic components
m_index = nm;
n_index= nn;
        
% kx, ky calculation
kx_f =@(k,theta,phi) k*sin(theta)*cos(phi);
ky_f =@(k,theta,phi) k*sin(theta)*sin(phi);
        
% kx_aux and ky_aux - Matrices
kx_aux_f = @(mu_x,m,Lx,m_index,n_index)(ones(m_index,n_index)*mu_x+2*MM*pi)/Lx;
ky_aux_f = @(mu_y,n,Ly,m_index,n_index)(ones(m_index,n_index)*mu_y+2*NN*pi)/Ly;
k_matricial_f= @(k,m_index,n_index) k*ones(m_index,n_index);
        
% Wavenumber - Paper - Equation  
f_kz = @(k_matricial,kx_aux,ky_aux) -(1j)*sqrt(kx_aux.^2+ky_aux.^2-k_matricial.^2);
    
% Initializating variables for STL and displacements
tau_total = zeros(nfreq,1);
disp = zeros(nfreq,1);

%Saving the amplitude of the harmonics to analyze their values
W1_mn_store= zeros(m_index,n_index,nfreq);


%% FREQUENCY LOOP - FRF(vacuum), DISPLACEMENT (fluid loading included) AND STL ANALYSIS %%%%%%%%
% Begin frequency loop
parfor i=1:nfreq

    A1 = freq(i);   
    formatSpec = ' \n Frequency =%4.2f Hz \n';
    fprintf(formatSpec,A1)    

    w=omega(i);

    %Dynamic stiffness matrix of the problem - plate with resonators
    D = K_new - w^2*M_new;

    % % Part 5 - FRF calculation, material in vaccum - OPTIONAL
    % % bare plate
    % D_struct = KG - w^2*MG;
    % u = D_struct\F;
    % u_f(:,i)=u;
    % u_z(i) = u(output);

    %end Part 5

    %% Part 6 - Wavenumber calculation in air and at the material

    %Wavenumbers of the incident wave (air)
    k= w/c0;

    kx = kx_f(k,theta,phi);
    ky = ky_f(k,theta,phi);

    %Propagation constants in x and y directions in the material
    mu_x = kx*Lx;
    mu_y = ky*Ly;

    %Wavenumber in z direction (air) - depending on harmonics - Matricial form

    kx_aux=kx_aux_f(mu_x,m,Lx,m_index,n_index);
    ky_aux=ky_aux_f(mu_y,n,Ly,m_index,n_index);

    %Calling wavenumber function to calculate kz1mn
    k_matricial= k_matricial_f(k,m_index,n_index);
    kz1mn = f_kz(k_matricial,kx_aux,ky_aux); 
    kz2mn = kz1mn; %same fluid

    %Dynamic Stiffness of the fluids
    rho1=rho0;
    rho2=rho0;
    %Equation 10
    Df1mn = (1j*rho1*w^2)./(kz1mn);
    Df2mn = (1j*rho2*w^2)./(kz2mn); 

    % end Part 6

    %% Part 7 - External force calculation (Acoustic force)

    % initializating force/moments vector related to blocked pressure wave
    Fext2=zeros(GDof+numberRes,1);
    
    % end Part 7

    %% Part 8 - Fluid loading calculation - Reaction of the fluid to the plate
    %displacement - This calculation depends on the spatial harmonics created
    %on the plate

    %Summation for harmonic components to calculate Df - Equation 40    
    [Df,V_store_H,Fext] = force_fluid_lumped2(GDof,numberRes,m_index,n_index,nnode,nelem,...
    elem,nodes,dof, Df1mn,kx_aux,ky_aux,kx,ky,A_element);
    
    %Adding rotational degrees of freedom
    Fext2(1:3:GDof)= -2*P_inc*Fext;
    
    %Adding fluid effects to dynamic stiffness matrix
    %Until this point the index of D is in the original order - matrices K, M
    %Equation 32
    D_til = D+(1/(Lx*Ly))*Df;

    %end Part 8

    %% Part 9 - WFE Model - Applying periodical conditions and the equilibrium
    % of forces at the boundary of the periodic cell

    %Separating internal and boundary dofs
    %Change order of boundary dofs for dofs in this order (q1,q2,q3,q4,ql,qb,qr,qt)

    %Partitioning Dynamic Stiffness Matrix
    D_til_bb= D_til(dofs_b,dofs_b);
    D_til_bi= D_til(dofs_b,dofs_i);
    D_til_ib= D_til(dofs_i,dofs_b);
    D_til_ii= D_til(dofs_i,dofs_i);

    % Partitioning external force (Acoustic)
    Fext_b = Fext2(dofs_b);
    Fext_i = Fext2(dofs_i);


    % Wave propagation constants - equation 36
    lambda_x = exp(-1j*mu_x);
    lambda_y = exp(-1j*mu_y);

    % fat - number of nodes at the lateral of the cell, excluding the corners
    fatx = NgridX-2; 
    faty = NgridY-2;
    
    % Generalized matrix for both square and rectangular unit cells 
    Lambda_R = [         eye(dof,dof), zeros(dof,faty*dof), zeros(dof,fatx*dof);...
            lambda_x*eye(dof,dof), zeros(dof,faty*dof), zeros(dof,fatx*dof);...
            lambda_y*eye(dof,dof), zeros(dof,faty*dof), zeros(dof,fatx*dof);...
            lambda_x*lambda_y*eye(dof,dof), zeros(dof,faty*dof), zeros(dof,fatx*dof);
            zeros(faty*dof,dof), eye(faty*dof,faty*dof),zeros(faty*dof,fatx*dof);
            zeros(fatx*dof,dof), zeros(fatx*dof,faty*dof), eye(fatx*dof,fatx*dof);
            zeros(faty*dof,dof),lambda_x*eye(faty*dof,faty*dof), zeros(faty*dof,fatx*dof);
            zeros(fatx*dof,dof), zeros(fatx*dof,faty*dof),lambda_y*eye(fatx*dof,fatx*dof)];

    % Equation 41
    Lambda_L = Lambda_R';        

    % na - number of dofs at the boundary
    na = length(dofs_b);
    % ni - number of internal dofs
    ni = length(dofs_i);
    % nc- numbers of dofs in qc - qc = [q1,ql,qb] - equation 38 
    nc = dof + fatx*dof +faty*dof;

    % Swicth between Not condensed problem and condensed problem
    % wfeproblemnotcondensed
%     wfeproblemdynamiccondensed
    D_r_til = D_til_bb-D_til_bi*(D_til_ii\D_til_ib);
    F_r_til = Fext_b - (D_til_bi/D_til_ii)*Fext_i;

    %Equation 43
    D_t = Lambda_L*D_r_til*Lambda_R;
    e_t = Lambda_L*F_r_til;

    %Solving - Equation 42

    Sol = D_t\e_t;

    %q = [qc,qi]
    % the order here is (q1,q2,q3,q4,ql,qb,qr,qt)
    %Equation 38
    q_bound = Lambda_R*Sol;

    %Recovering internal nodes
    %Equation 
    q_internal = (D_til_ii\Fext_i) -  (D_til_ii\D_til_ib)*q_bound;

    %Changing to the order of the original FE problem
    q_total = zeros(GDof+numberRes,1);
    q_total(dofs_b) = q_bound;
    q_total(dofs_i) = q_internal;

    % end Part 9


    %% Part 10 - Spatial harmonics amplitude calculation
    
    % Calculation of the amplitudes of the harmonics
    % Equation 28 and 29
    W1_mn= zeros(m_index,n_index);
    for ii =1:m_index
        for jj = 1:n_index
            v_aux = V_store_H(:,ii,jj);
            v_aux1 = zeros(GDof+numberRes,1);
            v_aux1(1:dof:GDof) = v_aux;
            %Equation 29
            % without the complex conjugate because of the reformulation
            %using shape functions in equation 28
            W1_mn(ii,jj) = (1/(Lx*Ly))*v_aux1.'*q_total;
        end
    end

    W2_mn = W1_mn; %same degrees of freedom on the two sides of the plate
    W1_mn_store(:,:,i)=W1_mn;

    % Amplitude of the sound pressure harmonics
    % Equation 22 - deduction
    P2_mn =  -Df2mn.*W2_mn./exp(-1j*kz2mn*h);


    %% Part 11 -Displacement calculation - Summation of harmonics components

    %Positions where displacement is evaluated 
    %center of the plate
    x_disp=0;
    y_disp=0;

    %Equation 4
    aux_disp = exp(-1j*(kx_aux*x_disp+ky_aux*y_disp));
    disps_mn = W1_mn.* aux_disp;

    %summation
    aux = sum(sum(disps_mn,2));
    disp(i) = aux;

    % end Part 11

    %Part 12 - sound transmission coefficient
    %Equation 45
    tau2_aux=(1/rho2)*kz2mn;
    tau1_aux =(1/rho1)*kz1mn;

    tau_mn_num= real(tau2_aux.*abs(P2_mn.^2));
    tau_mn_den = tau1_aux.*abs((P_inc*ones(m_index,n_index)).^2);

    tau_mn = tau_mn_num./tau_mn_den;
    %summation
    tau = sum(sum(tau_mn),2);     

    tau_total(i)=tau;

end

%% HARMONIC AMPLITUDE GRAPHS
figure
W1_mn_interest = W1_mn_store(:,:,57);
max_W1_mn=max(max(W1_mn_interest));
heatmap(m,n,log10(abs(W1_mn_interest./max_W1_mn)))

%% STL AND PLATE DISPLACEMENT (WITH FLUID LOADING) GRAPHS 

%Part 13 - STL calculation
 STL = -10*log10(tau_total);

save('STL','STL')

figure
loglog(freq,abs(disp))
title('Displacement at the center of the plate with resonator(s)')
xlabel('Frequency [Hz]')
ylabel('w [m]')
legend('Displacement - WFE')
%axis([fmin,fmax,1*10-6,1*10-4])

grid on

% saveas(gcf,'displacement.jpg')


figure
semilogx(freq,STL)
title('Sound transmission loss for a plate with resonator(s)')
legend('STL - WFE')
xlabel('Frequency [Hz]')
ylabel('STL [dB]')
axis([fmin,fmax,0,70])
grid on

% saveas(gcf,'STL.jpg')

toc
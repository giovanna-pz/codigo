clear
% close all
clc

%% This program calculates the displacements and the STL of a metamaterial plate 

% Infinite plate with periodical resonators using the WFE method and the assumption that in a periodical
%material space harmonics are created.
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
%plotMesh(nodes,elem)

% Mass and stiffness matrices - substructure plate
% Choose for desired model
kircchoff_model
%mindlin_model

% Updating K and M with resonators
[~,node_res_first]=min(sqrt(nodes(:,1).^2)+nodes(:,2).^2);
%center node, where the resonator is attached
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
% calc_disp_plot_modes(Lx,Ly,K_new,M_new,dofs_i,dofs_b,dof,Ngrid,el_side_x,el_side_y,nnode)
% end Part 4

%%% FREQUENCY LOOP - FRF(vacuum), DISPLACEMENT (fluid loading included) AND STL ANALYSIS %%%%%%%%

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

%teste 2 - numero de onda
 %f_kz = @(k_matricial,kx_aux,ky_aux) -sqrt(k_matricial.^2-kx_aux.^2-ky_aux.^2).*(k_matricial.^2-kx_aux.^2-ky_aux.^2 > 0)+...
  %    (+1j)*sqrt(k_matricial.^2-kx_aux.^2-ky_aux.^2).*(k_matricial.^2-kx_aux.^2-ky_aux.^2 < 0);

%teste 3 - numero de onda
%  f_kz = @(k_matricial,kx_aux,ky_aux) sqrt(k_matricial.^2-kx_aux.^2-ky_aux.^2).*(k_matricial.^2-kx_aux.^2-ky_aux.^2 >= 0)+...
%   (-1j)*sqrt(kx_aux.^2-ky_aux.^2-k_matricial.^2).*(k_matricial.^2-kx_aux.^2-ky_aux.^2 < 0);  
  
% Initializating variables for STL and displacements
tau_total = zeros(nfreq,1);
disp = zeros(nfreq,1);

% Saving wavenumbers vector and matrices to analyze their behavior in plots
kx_p = zeros(1,nfreq);
ky_p = zeros(1,nfreq);
kx_aux_p = zeros(nfreq,m_index,n_index);
ky_aux_p = zeros(nfreq,m_index,n_index);
kzmn_p = zeros(nfreq,m_index,n_index);


% Saving the amplitude of the harmonics to analyze their values
W1_mn_store= zeros(m_index,n_index,nfreq);

% Begin frequency loop
for i=1:nfreq
    
    % Printing the current frequency step
    A1 = freq(i);   
    formatSpec = ' \n Frequency =%4.2f Hz \n';
    fprintf(formatSpec,A1)    

    w=omega(i);

    % Dynamic stiffness matrix of the problem - plate with resonators
    D = K_new - w^2*M_new;

    %% Part 5 - FRF calculation, material in vaccum - OPTIONAL
    % % bare plate
    % D_struct = KG - w^2*MG;
    % u = D_struct\F;
    % u_f(:,i)=u;
    % u_z(i) = u(output);

    % end Part 5

    %% Part 6 - Wavenumber calculation in air and at the material
    % Wavenumbers of the incident wave (air)
    k= w/c0;

    kx = kx_f(k,theta,phi);
    ky = ky_f(k,theta,phi);

    % Storing wavenumbers to plot
    kx_p(i) = kx;
    ky_p(i) = ky;

    % Propagation constants in x and y directions in the material
    mu_x = kx*Lx;
    mu_y = ky*Ly;

    % Wavenumber in z direction (air) - depending on harmonics - Matricial form
    kx_aux=kx_aux_f(mu_x,m,Lx,m_index,n_index);
    ky_aux=ky_aux_f(mu_y,n,Ly,m_index,n_index);

    % Storing wavenumbers to plot
    kx_aux_p(i,:,:) = kx_aux;
    ky_aux_p(i,:,:) = ky_aux;

    % Calling wavenumber function to calculate kz1mn
    k_matricial= k_matricial_f(k,m_index,n_index);
    kz1mn = f_kz(k_matricial,kx_aux,ky_aux); 
    kz2mn = kz1mn; %same fluid

    % Storing wavenumber to plot
    kzmn_p(i,:,:) = kz1mn;

    % Dynamic Stiffness of the fluids
    rho1=rho0;
    rho2=rho0;
    
    %Equation 
    Df1mn = (1j*rho1*w^2)./(kz1mn);
    Df2mn = (1j*rho2*w^2)./(kz2mn); 

    % end Part 6

    %% Part 7 - External force calculation (Acoustic force)

    % Initializating force/moments vector
    Fext2=zeros(GDof+numberRes,1);

    % Vector of external forces  calculated only for
    % displacement degrees of freedom
    % Equation 18/19, but using shape functions
    alpha = force_external(nnode,nelem,elem,nodes,csi_aux,eta_aux,wcsi_aux,weta_aux,kx,ky);

    % Adding rotational degrees of freedom
    Fext2(1:3:GDof)= alpha;
    
    % end Part 7

    %% Part 8 - Fluid loading calculation - Reaction of the fluid to the plate
    % displacement - This calculation depends on the spatial harmonics created
    % on the plate

    % Summation for harmonic components to calculate Df - Equation 33

    % Based on equations 25,26 and 28 using shape functions
    [Df,V_store_H,s_pinc] = force_fluid(GDof, numberRes,m_index,n_index,nnode,nelem,elem,nodes,...
    csi_aux,eta_aux, wcsi_aux, weta_aux, dof, Df1mn,kx_aux,ky_aux,kz2mn);

    % Adding fluid effects to dynamic stiffness matrix

    % Until this point the index of D is in the original order - matrices K, M
    % Equation 32
    D_til = D+(1/(Lx*Ly))*Df;

    
    %Appling correction to FEXt 2 - Equation 3 deduction Giovanna
    M=n_index*m_index;
    Fext2=(-P_inc-(P_inc*k*cos(theta)/M)*s_pinc)*Fext2;
    
    % end Part 8

    %% Part 9 - WFE Model - Applying periodical conditions and the equilibrium
    % of forces at the boundary of the periodic cell

    % Separating internal and boundary dofs
    % Change order of boundary dofs for dofs in this order (q1,q2,q3,q4,ql,qb,qr,qt)

    % Partitioning Dynamic Stiffness Matrix
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
    fat = Ngrid-2; 

    % The Lambda_R matrix in the article contains typography errors 
    % Equation 38
    %Matriz Giovanna 
    Lambda_R = [         eye(dof,dof), zeros(dof,fat*dof), zeros(dof,fat*dof);...
            lambda_x*eye(dof,dof), zeros(dof,fat*dof), zeros(dof,fat*dof);...
            lambda_y*eye(dof,dof), zeros(dof,fat*dof), zeros(dof,fat*dof);...
            lambda_x*lambda_y*eye(dof,dof), zeros(dof,fat*dof), zeros(dof,fat*dof);
            zeros(fat*dof,dof), eye(fat*dof,fat*dof),zeros(fat*dof,fat*dof);
            zeros(fat*dof,dof), zeros(fat*dof,fat*dof), eye(fat*dof,fat*dof);
            zeros(fat*dof,dof),lambda_x*eye(fat*dof,fat*dof), zeros(fat*dof,fat*dof);
            zeros(fat*dof,dof), zeros(fat*dof,fat*dof),lambda_y*eye(fat*dof,fat*dof)];

    % Equation 41
    Lambda_L = Lambda_R';        


    % na - number of dofs at the boundary
    na = length(dofs_b);
    
    % ni - number of internal dofs
    ni = length(dofs_i);
    
    % nc- numbers of dofs in qc - qc = [q1,ql,qb] - equation 38 
    nc = dof + 2*fat*dof;

    % Swicth between Not condensed problem and condensed problem
    %wfeproblemnotcondensed    % Not condensed
    wfeproblemdynamiccondensed % Condesed
   
    % Changing to the order of the original FE problem
    q_total = zeros(GDof+numberRes,1);
    q_total(dofs_b) = q_bound;
    q_total(dofs_i) = q_internal;

    % end Part 9


    %% Part 10 - Spatial harmonics amplitude calculation

    % Calculation of the amplitudes of the harmonics
    %Equation 28 and 29
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
    % Equation 22
    P2_mn =  -Df2mn.*W2_mn.*(1./exp(-1j*kz2mn*h));

    %end Part 10


    %% Part 11 - Displacement calculation - Summation of harmonics components

    % Positions where displacement is evaluated 
    % center of the plate
    x_disp=0;
    y_disp=0;

    % Equation 4
    aux_disp = exp(-1j*(kx_aux*x_disp+ky_aux*y_disp));
    disps_mn = W1_mn.* aux_disp;

    % Summation
    disp(i) = sum(sum(disps_mn,2));

    % end Part 11

    %% Part 12 - sound transmission coefficient
    % Equation 45

    tau2_aux=(1/rho2)*kz2mn;
    tau1_aux =(1/rho1)*kz1mn;

    tau_mn_num= real(tau2_aux.*abs(P2_mn.^2));
    tau_mn_den = tau1_aux.*abs((P_inc*ones(m_index,n_index)).^2);
    % tau_mn_num= real(abs(P2_mn.^2));
    % tau_mn_den = abs((P_inc*ones(m_index,n_index)).^2);

    tau_mn = tau_mn_num./tau_mn_den;
    
    % Summation
    tau = sum(sum(tau_mn),2);     


    tau_total(i)=tau;

end

%%  Plotting section 

% %% Calculating FRFs 
% load freq_vector_comsol
% load Recpt_comsol
% 
% %Receptance
% alfa = u_z/F0;
% semilogx(f,20*log10(abs(alfa)),'k','LineWidth',2)
% hold on
% plot(freq_vector_comsol,Recpt_comsol,'r--','LineWidth',2)
% title('Receptance at the center of the bare plate')
% xlabel('Frequency [Hz]')
% ylabel('Receptance [dB]')
% legend('FEM - Matlab','Reference - Comsol')
% grid on

%% Amplitude of harmonics

% figure
% W1_mn_interest = W1_mn_store(:,:,20);
% stem3(m,n,real(W1_mn_interest),'filled')
% hold on
% stem3(m,n,imag(W1_mn_interest),'--*r')
% title('Absolute values of harmonics for a given frequency')

figure
W1_mn_interest = W1_mn_store(:,:,57);
max_W1_mn = W1_mn_interest(5,5);
heatmap(m,n,log10(abs(real(W1_mn_interest./max_W1_mn))))
title('Normalized values of harmonics for a given frequency')


%% Wavenumbers graphs

% figure
% plot(f,kx_p,'LineWidth',2)
% hold on
% plot(f,ky_p,'LineWidth',2)
% title('Wavenumbers of the excitation waves')
% legend('k_x','k_y')
% xlabel('Frequency [Hz]')
% ylabel('k [2\pi/m]')
% 
% 
% figure
% for pp =1:m_index
%     for pp2 = 1:n_index
%         hold on       
%         plot(f,kx_aux_p(:,pp,pp2),'LineWidth',2)
%         hold on
%     end
% end
% title(' Spatial harmonics of kx at the plate')
% xlabel('Frequency [Hz]')
% ylabel('k_x_aux [2\pi/m]')
% 
% figure
% for pp3 =1:m_index
%     for pp4 = 1:n_index
%         hold on       
%         plot(f,ky_aux_p(:,pp3,pp4),'LineWidth',2)
%         hold on
%     end
% end
% title(' Spatial harmonics of ky at the plate')
% xlabel('Frequency [Hz]')
% ylabel('k_y_aux [2\pi/m]')
% 
% 
% figure
% % for pp5 =1:m_index
% %     for pp6 = 1:n_index
% % 
% % hold on       
% % plot(f,abs(real(kzmn_p(:,pp5,pp6))),'LineWidth',2)
% % hold on
% % plot(f,-abs(imag(kzmn_p(:,pp5,pp6))),'LineWidth',2)
% %     end
% % end
% 
% plot(f,abs(real(kzmn_p(:,5,1))),'LineWidth',2)
% hold on
% plot(f,-abs(imag(kzmn_p(:,5,1))),'LineWidth',2)
% hold on
% plot(f,abs(real(kzmn_p(:,6,1))),'LineWidth',2)
% hold on
% plot(f,-abs(imag(kzmn_p(:,6,1))),'LineWidth',2)
% hold on
% plot(f,abs(real(kzmn_p(:,7,1))),'LineWidth',2)
% hold on
% plot(f,-abs(imag(kzmn_p(:,7,1))),'LineWidth',2)
% hold on
% plot(f,abs(real(kzmn_p(:,8,1))),'LineWidth',2)
% hold on
% plot(f,-abs(imag(kzmn_p(:,8,1))),'LineWidth',2)
% hold on
% plot(f,abs(real(kzmn_p(:,9,1))),'LineWidth',2)
% hold on
% plot(f,-abs(imag(kzmn_p(:,9,1))),'LineWidth',2)
% 
% legend('Real m=0,n=0','Imag m=0,n=0','Real m=1,n=0','Imag m=1,n=0',...
%     'Real m=2,n=0','Imag m=2,n=0','Real m=3,n=0','Imag m=3,n=0', 'Real m=4,n=0','Imag m=4,n=0')
% 
% title('Postive Harmonics of kzmn at the fluid')
% xlabel('Frequency [Hz]')
% ylabel('kz_mn [2\pi/m]')  
% 
% figure
% 
% plot(f,abs(real(kzmn_p(:,1,1))),'LineWidth',2)
% hold on
% plot(f,-abs(imag(kzmn_p(:,1,1))),'LineWidth',2)
% hold on
% plot(f,abs(real(kzmn_p(:,2,1))),'LineWidth',2)
% hold on
% plot(f,-abs(imag(kzmn_p(:,2,1))),'LineWidth',2)
% hold on
% plot(f,abs(real(kzmn_p(:,3,1))),'LineWidth',2)
% hold on
% plot(f,-abs(imag(kzmn_p(:,3,1))),'LineWidth',2)
% hold on
% plot(f,abs(real(kzmn_p(:,4,1))),'LineWidth',2)
% hold on
% plot(f,-abs(imag(kzmn_p(:,4,1))),'LineWidth',2)
% 
% 
% 
% legend('Real m=-4,n=0','Imag m=-4,n=0','Real m=-3,n=0','Imag m=-3,n=0',...
%     'Real m=-2,n=0','Imag m=-2,n=0','Real m=-1,n=0','Imag m=-1,n=0')
% 
% title('Negative Harmonics of kzmn at the fluid')
% xlabel('Frequency [Hz]')
% ylabel('kz_mn [2\pi/m]')  



%% STL and pltae displacement(with fluid loading) graphs

%Part 13 - STL calculation
 STL = -10*log10(tau_total);

 save('STL','STL')

%To be compared with figure 6b of the paper
figure
load f 
load uz_damp
loglog(freq,abs(disp),f,abs(uz_damp))
title('Displacement at the center of the metamaterial plate')
xlabel('Frequency [Hz]')
ylabel('w [m]')
legend('Displacement - WFE','Displacement - Reference - \phi=0, \theta = 50 ')
%axis([fmin,fmax,1*10-6,1*10-4])

grid on

 saveas(gcf,'displacement.jpg')

%To be compared with figure 10 of the paper
figure
semilogx(freq,STL)
hold on

load STL_damp
semilogx(f,STL_damp)
title('Sound Transmission Loss for a Metamaterial')
legend('STL - WFE','Reference - \phi=0, \theta = 50 ')
xlabel('Frequency [Hz]')
ylabel('STL [dB]')
axis([fmin,fmax,0,70])
grid on


% saveas(gcf,'STL.jpg')
toc
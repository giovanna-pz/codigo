function calc_disp_plot_modes(Lx,Ly,K_new,M_new,dofs_i,dofs_b,dof,NgridX,NgridY,el_side_x,el_side_y,theta,fmax,numberRes)

% %Plate in vacuum - no fluid loading
% % No structural damping

%alterar aqui pra ficar generico

fatx = NgridX-2; %number of nodes at qb and qt
faty = NgridY-2; %number of nodes at qb and qt
%qc = [q1,ql,qb]
nc = dof + fatx*dof +faty*dof; %number of dofs at vector qc
na = length(dofs_b);

if Lx == Ly
    [omega_disp,N, modes_full]=dispersion_2D_square(Lx,Ly,nc,K_new,M_new,dofs_i,dofs_b,dof,fatx,faty);
    n= N/3;
    % Graphs
    freq_disp =real(omega_disp)/(2*pi);

    %Complete graph
    %graph Frequency vs. k 
    figure
    plot(1:N,freq_disp, 'r.-' )
    grid on
    axis([1 N 0 fmax])

    set(gca,'XTick',[1 n 2*n 3*n-1],'XTickLabel',{'T','X','M','T'})
    ylabel('Frequency [Hz]')
    grid on
    % saveas(gcf,'dispersion.jpg')
else
    [omega_disp,N, modes_full]=dispersion_2D_rectangular(Lx,Ly,nc,K_new,M_new,dofs_i,dofs_b,dof,fatx,faty);
    n=N/4;
        
    freq_disp =real(omega_disp)/(2*pi);

    %Complete graph
    %graph Frequency vs. k 
    figure
    plot(1:N,freq_disp, 'r.-' )
    grid on
    axis([1 N 0 fmax])

    set(gca,'XTick',[1 n 2*n 3*n 4*n-1],'XTickLabel',{'T','X','M','K','T'})
    ylabel('Frequency [Hz]')
    grid on
    % saveas(gcf,'dispersion.jpg')
end




% Graph only in x direction
%graph k vs. frequency 

%Converting points to wavenumber
k_plot = ((1:n)-ones(1,n))*(1/(n-1));

figure
plot(freq_disp(:,[1:n]),k_plot, 'r.-' )
hold on
f_air = 1:1:fmax;
c_air=340;
k_air = (2*pi*f_air/(c_air))*(Lx/pi)*sin(theta);
plot(f_air,(k_air),'k')
grid on
axis([0 fmax 0 k_plot(n)])
ylabel('Wavenumber [1/m]')
xlabel('Frequency [Hz]')
grid on
%saveas(gcf,'dispersion.jpg')

% WAVE MODES GRAPHS

% F = 330 HZ -> N=8, FREQ_DISP(:,8)= 330 na posicao 2, logo, o modo
% equivalente e modes_full(:,2,8);

% F = 725 HZ -> N=8, FREQ_DISP(:,17)= 729 na posicao 3, logo, o modo
% equivalente e modes_full(:,3,17);

x = -Lx/2:el_side_x:Lx/2;
y = -Ly/2:el_side_y:Ly/2;

[XX,YY]=meshgrid(x,y);

mode = modes_full(:,3,10);

% modo - order: q = [qc,qi]
%back to the order of the FE Model
mode_bound = mode(1:na);
mode_internal = mode(na+1:end);


mode_fem(dofs_b) = mode_bound ;
mode_fem(dofs_i) = mode_internal;

%removing translational dofs from resonator
trans_mode = mode_fem(1:3:end-numberRes);


modo1 = reshape(trans_mode,[NgridX,NgridY]);


figure
surf(XX,YY,(real(modo1)).')

end
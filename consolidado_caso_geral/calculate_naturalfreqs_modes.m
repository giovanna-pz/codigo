% Determing natural frequencies
[V,D1]= eig(real(KG),MG); 
%d = eigs(A,10,'smallestabs')
[omega2_nat,index] = sort(diag(D1));
%Modos
V_new = V(:,index);
omega_nat = sqrt(omega2_nat);
freqs = omega_nat/(2*pi);
freqs_int=freqs(1:17);

% Comsol results for natural frequencies
load freqs_comsol.mat

%Compare results

nat_freqs = [freqs_int,freqs_comsol]
%% Plotting the modes
x = -Lx/2:el_side_x:Lx/2;
y = -Ly/2:el_side_y:Ly/2;

[XX,YY]=meshgrid(x,y);
desl = V_new(1:3:end,:);
modo = 4;
modo1 = reshape(desl(:,modo),[sqrt(nnode),sqrt(nnode)]);
figure
surf(XX,YY,modo1)
function [Df,V_store_H,Fext] = force_fluid_lumped2(GDof,numberRes,m_index,n_index,nnode,nelem,...
    elem,nodes,dof, Df1mn,kx_aux,ky_aux,kx,ky,A_element)

%kx,ky
%V_store_H it is a 3D Matrix: 
% lines: nodes
% columns: harmonics in x (m_index)
% Layers: harmonics in y (n_index)
V_store_H =zeros(nnode,m_index,n_index);

%Df - fluid loading to be added at Dynamic stiffness matrix
Df=zeros(GDof+numberRes,GDof+numberRes); %start summation

% Fluid forces calculated only for displacement degrees of freedom        
Ffluid = zeros(m_index,n_index,nnode);
Ffluid_H = zeros(m_index,n_index,nnode);
Fext=zeros(nnode,1);
          
for uu = 1:nelem
    %Identifying the nodes of the element 
    index = elem(uu,:);

    %Identifying coordinates of the nodes
    xcoord2 = nodes(index,1); %column vector
    ycoord2 = nodes(index,2);

    
    xmean= mean(xcoord2);
    ymean = mean(ycoord2);
    
    % Computing the lumped fluid force at each element
    f_mean= 0.25*A_element*exp(-1j*((kx_aux*xmean)+(ky_aux*ymean)));
    f_mean_H = 0.25*A_element*exp(+1j*((kx_aux*xmean)+(ky_aux*ymean)));

    Ffluid(:,:,index(1)) = Ffluid(:,:, index(1)) + f_mean;
    Ffluid(:,:,index(2)) = Ffluid(:,:, index(2)) + f_mean;
    Ffluid(:,:,index(3)) = Ffluid(:,:, index(3)) + f_mean;
    Ffluid(:,:,index(4)) = Ffluid(:,:, index(4)) + f_mean;

    Ffluid_H(:,:,index(1)) = Ffluid_H(:,:, index(1)) + f_mean_H;
    Ffluid_H(:,:,index(2)) = Ffluid_H(:,:, index(2)) + f_mean_H;
    Ffluid_H(:,:,index(3)) = Ffluid_H(:,:, index(3)) + f_mean_H;
    Ffluid_H(:,:,index(4)) = Ffluid_H(:,:, index(4)) + f_mean_H;
    
    % Computing the external fluid force at each element
    f_mean = A_element*exp(-1j*(kx*xmean+ky*ymean));
    fe = [f_mean/4; f_mean/4;f_mean/4;f_mean/4];
    Fext(index) = Fext(index)+fe;
end
            
for kk = 1:m_index
    for ll =1:n_index
        v1_mn=zeros(GDof+numberRes,1);
        v1_mn_H = zeros(GDof+numberRes,1);
        
        %adding rotation degrees of freedom
        v1_mn(1:dof:GDof) = Ffluid(kk,ll,:);
        v1_mn_H(1:dof:GDof) = Ffluid_H(kk,ll,:);
         
          
        %Storing this vector - To be used in future calculation
        V_store_H(:,kk,ll) = Ffluid_H(kk,ll,:);
          
        %Equation 33
        %(v1_mn_H.')- transpose without the complex conjugate
        Df = Df + 2*Df1mn(kk,ll)*(v1_mn*v1_mn_H.'); %Equation 33
          
    end
end
 
end
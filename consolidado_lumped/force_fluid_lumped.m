function [Df,V_store_H] = force_fluid_lumped(GDof, numberRes,m_index,n_index,nnode,nelem,elem,nodes,...
     dof, Df1mn,Df2mn,kx_aux,ky_aux,A_element)

%V_store_H: 3D Matrix: lines: nodes
%columns: harmonics in x (m_index)
%Layers: harmonics in y (n_index)
V_store_H =zeros(nnode,m_index,n_index);

%Df - fluid loading to be added at Dynamic stiffness matrix
Df=zeros(GDof+numberRes,GDof+numberRes); %start summation

 for kk = 1:m_index
  for ll =1:n_index
          
          kx_aux_aux = kx_aux(kk,ll);
          ky_aux_aux = ky_aux(kk,ll);
      
          v1_mn=zeros(GDof+numberRes,1);
          v1_mn_H = zeros(GDof+numberRes,1);
          %Fluid forces calculated only for displacement degrees of freedom
          Ffluid = zeros(nnode,1);
          Ffluid_H = zeros(nnode,1);
          
          for uu = 1:nelem
              
              %Identifying the nodes of the element 
              index2 = elem(uu,:);
              %%Identifying coordinates of the nodes
              xcoord2 = nodes(index2,1); %column vector
              ycoord2 = nodes(index2,2);

              xmean= mean(xcoord2);
              ymean = mean(ycoord2);
              
              f_mean= A_element*exp(-1j*((kx_aux_aux*xmean)+(ky_aux_aux*ymean)));
              f_mean_H = A_element*exp(+1j*((kx_aux_aux*xmean)+(ky_aux_aux*ymean)));
              
              f_fluid = [f_mean/4; f_mean/4;f_mean/4;f_mean/4];
              f_fluid_H = [f_mean_H/4; f_mean_H/4;f_mean_H/4;f_mean_H/4];

              Ffluid(index2) = Ffluid(index2)+f_fluid;
              Ffluid_H(index2) = Ffluid_H(index2)+f_fluid_H;
          end
          
          %adding rotation degrees of freedom
          v1_mn(1:dof:GDof) = Ffluid;
          v1_mn_H(1:dof:GDof) = Ffluid_H;
          
          %expression for the displacement at the resonator mass
%           v1_mn(GDof+numberRes)=v1_mn(node_res_first*dof-2)/(1-(w/sqrt(kr/mr))); 
%           v1_mn_H(GDof+numberRes)=v1_mn_H(node_res_first*dof-2)/(1-(w/sqrt(kr/mr)));
% nao usar, gera uma descontinuidade na curva
          
          %Storing this vector - To be used in future calculation
          V_store_H(:,kk,ll) = Ffluid_H;
          
          %Equation 33
          %(v1_mn_H.')- transpose without the complex conjugate
          %Equation 42 deduction Giovanna
          Df = Df + Df1mn(kk,ll)*(v1_mn*v1_mn_H.')+Df2mn(kk,ll)*(v1_mn*v1_mn_H.') ; %Equation 33%Df = Df + 2*Df1mn(kk,ll)*(v1_mn*v1_mn_H.'); %Equation 33
  end
 end
 
end
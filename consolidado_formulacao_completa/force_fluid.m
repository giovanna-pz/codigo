function [Df,V_store_H,s_pinc] = force_fluid(GDof, numberRes,m_index,n_index,nnode,nelem,elem,nodes,...
    csi_aux,eta_aux,wcsi_aux, weta_aux, dof, Df1mn,kx_aux,ky_aux,kz1mn)

%V_store_H: 3D Matrix: 
% lines: nodes
% columns: harmonics in x (m_index)
% Layers: harmonics in y (n_index)
V_store_H =zeros(nnode,m_index,n_index);

%Df - fluid loading to be added at Dynamic stiffness matrix
Df=zeros(GDof+numberRes,GDof+numberRes); %start summation

%Summation of term added to incidence force

s_pinc=0;

 for kk = 1:m_index
  for ll =1:n_index
          
      
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

              
              %numerical integration
              f_fluid=zeros(4,1);
              f_fluid_H=zeros(4,1);
              
              for cc =1:length(csi_aux)
                  for dd =1:length(eta_aux)
                      %Equation 26
                      csi2 = csi_aux(cc);
                      eta2 = eta_aux(dd);
                      [N2,detJ2]= Quad(csi2,eta2,xcoord2,ycoord2);
                      x2 = N2.'*xcoord2;
                      y2 = N2.'*ycoord2;
                      
                      aux_N2 = wcsi_aux(cc)*weta_aux(dd)*N2;
                      
                      f_fluid = aux_N2*exp(-1j*((kx_aux(kk,ll)*x2)+(ky_aux(kk,ll)*y2)))*detJ2+f_fluid;
                      f_fluid_H = aux_N2*exp(+1j*((kx_aux(kk,ll)*x2)+(ky_aux(kk,ll)*y2)))*...
                          detJ2+f_fluid_H;                  
                  
                  end
              end 

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
          Df = Df + 2*Df1mn(kk,ll)*(v1_mn*v1_mn_H.'); %Equation 33
          
          %Equation 31 - deduction Giovanna
          
          s_pinc = s_pinc + (1/kz1mn(kk,ll));
         
      
  
  
  end
 end
 
end
function Fext = force_external(nnode,nelem,elem,nodes,csi_aux,eta_aux,wcsi_aux,weta_aux,kx,ky)

Fext=zeros(nnode,1);

for jj=1:nelem
    
    %Identifying the nodes of the element 
    index = elem(jj,:);
    %%Identifying coordinates of the nodes
    xcoord = nodes(index,1); %column vector
    ycoord = nodes(index,2);
    
    fe=zeros(4,1);
    
    %Numerical integration for external forces 
    for aa=1:length(csi_aux)
        for bb=1:length(eta_aux)
    
            csi = csi_aux(aa);
            eta = eta_aux(bb);
            wcsi = wcsi_aux(aa);
            weta = weta_aux(bb);
            
            [N,detJ] = Quad(csi,eta,xcoord,ycoord);
            x = N.'*xcoord;
            y = N.'*ycoord;
            
         
            fe = wcsi*weta*N*exp(-1j*(kx*x+ky*y))*detJ + fe;
    
        end
    end
    
    
    Fext(index) = Fext(index)+fe;
    
end
 
end
function Fext = force_external_lumped(nnode,nelem,elem,nodes,P_inc,kx,ky, A_element)

Fext=zeros(nnode,1);

for jj=1:nelem
    
    %Identifying the nodes of the element 
    index = elem(jj,:);
    %%Identifying coordinates of the nodes
    xcoord = nodes(index,1); %column vector
    ycoord = nodes(index,2);
 
    xmean= mean(xcoord);
    ymean = mean(ycoord);
    
    
    f_mean = A_element*P_inc*exp(-1j*(kx*xmean+ky*ymean));
    
    fe = [f_mean/4; f_mean/4;f_mean/4;f_mean/4];
    
    Fext(index) = Fext(index)+fe;
    
end
 
end
                      % 
%-------------------------------------------------------------------------%
% Plot malha: Plota a malha de elementos finitos a partir dos nos e       %
%             conectividade                                               %
%-------------------------------------------------------------------------%
function plotmalha(coord,inci,numbering)
%se numbering é
% 1: Plota a malha em preto com os nós e os números de elementos
% 0: Plota a malha em preto somente
nel=length(inci(:,1));  %numero de elementos
nnos=length(coord(:,1));%numero de nós
nnel=size(inci(:,2:end),2);      %Numero de elementos por elemento
dimension=size(coord(:,2:end),2); %Dimensao da malha
%Initializacao das matrices
X = zeros(nnel,nel) ;
Y = zeros(nnel,nel) ;
if dimension==2            %For 2D plots
    for iel=1:nel
        for i=1:nnel
            X(i,iel)=coord(inci(iel,i+1),2);
            Y(i,iel)=coord(inci(iel,i+1),3);
        end
    end
    %Plotting the FEM mesh, display Node numbers and Element numbers
        plot(X,Y,'k') 
        title('Finite Element Mesh') ;
        fill(X,Y,'w') % Preenche os elementos de branco para ficar mas visivel
        axis off 
    %Numeracao
    if numbering==1
    text(coord(:,2),coord(:,3),int2str(coord(:,1)),'fontsize',16,'color','k');
       for i = 1:nel
            text(sum(X(:,i))/4,sum(Y(:,i))/4,int2str(i),'fontsize',16,'color','r');
       end
    end
  end
    




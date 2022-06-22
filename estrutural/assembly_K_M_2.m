function [KG,MG]= assembly_K_M_2(Ke, Me, dof,nel,nnos_el,GDof,inci)
KG = zeros(GDof);
MG = zeros(GDof);

for i = 1: nel 
    nodes = inci(i,2:nnos_el+1);
    
    elementDof = findNodeDofs(nodes,dof);
    
    
    KG(elementDof,elementDof) = KG(elementDof,elementDof) + Ke;
    MG(elementDof,elementDof) = MG(elementDof,elementDof) + Me;

end

end
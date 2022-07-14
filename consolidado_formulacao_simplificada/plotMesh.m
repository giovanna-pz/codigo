%Plot mesh - Verification
function plotMesh(nodes,elem)

figure, hold on

% element
for iel=1:size(elem,1)
    nodes_el = elem(iel,1:4);
    coord_nodes = nodes(nodes_el,:);
    fill(coord_nodes(:,1),coord_nodes(:,2), [1 1 1]);
    text(mean(coord_nodes(:,1)),mean(coord_nodes(:,2)),num2str(iel),'Color',[1 0 0])
end

% nodes
for in=1:size(nodes,1)
    scatter(nodes(in,1),nodes(in,2),'k')
    text(nodes(in,1),nodes(in,2),num2str(in))
end


end
function boundary_node_map = find_boundary (node, element);
%Copyright: Lampros Kourtis, 2007, Stanford University
%email    : kourtis@stanford.edu
node_count=zeros(length(node),1);

for t=1:length(element)
    for g=1:4
        node_count(element(t,g))=node_count(element(t,g))+1;
    end
end

boundary_node_map = node_count<3;
function [n, boundary_node_map] = mesh_normals (node, nodes_th, canal_flag)
%Copyright: Lampros Kourtis, 2007, Stanford University
%email    : kourtis@stanford.edu
n_nodes=length(node);
t=zeros(n_nodes,2);
boundary_node_map = zeros(n_nodes,1);

%find tangent of nodes by differentiating the contour
t(1:nodes_th,:)=diff([node(1:nodes_th,:);node(1,:)]);
t(1:nodes_th,:)=t(1:nodes_th,:)./repmat(sqrt(t(1:nodes_th,1).^2 + t(1:nodes_th,2).^2),1,2);
boundary_node_map(1:nodes_th)=1;

%do same for inner canal (if existent)
if canal_flag==1   
 t(n_nodes-nodes_th+1:n_nodes,:)=diff([node(n_nodes-nodes_th+1:n_nodes,:);node(n_nodes-nodes_th+1,:)]);   
 t(n_nodes-nodes_th+1:n_nodes,:)=-t(n_nodes-nodes_th+1:n_nodes,:)./repmat(sqrt(t(n_nodes-nodes_th+1:n_nodes,1).^2 + t(n_nodes-nodes_th+1:n_nodes,2).^2),1,2);   
 boundary_node_map(n_nodes-nodes_th+1:n_nodes)=1;
end

%invert for outwards
t=-t;

%Find normals from tangent vectors
n=[t(:,2),-t(:,1)];

%figure(7);quiver (node(:,1), node(:,2), n(:,1),n(:,2)); % plot normals
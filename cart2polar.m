function [node_polar] = cart2polar (node_cart, center)
%Copyright: Lampros Kourtis, 2007, Stanford University
%email    : kourtis@stanford.edu
node_cart=node_cart-repmat(center,length(node_cart),1);    %move to center
[node_polar(:,1) node_polar(:,2)]=cart2pol(node_cart(:,1),node_cart(:,2));   %convert to polar

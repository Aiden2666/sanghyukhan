function [node_cart] = polar2cart (node_polar, center);
%Copyright: Lampros Kourtis, 2007, Stanford University
%email    : kourtis@stanford.edu

%% Offset data to have a center at 0,0 and then use pol2cart
%center=fliplr(center);


% Assign cartesian values for nodes
    [node_cart(:,1) node_cart(:,2)] = pol2cart(node_polar(:,1), node_polar(:,2));
    node_cart(:,1) = node_cart(:,1)+center(1);
    node_cart(:,2) = node_cart(:,2)+center(2);



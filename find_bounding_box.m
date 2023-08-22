function [x_dimension, y_dimension ] = find_bounding_box (node_ext_cart, resolution);
%Copyright: Lampros Kourtis, 2007, Stanford University
%email    : kourtis@stanford.edu
x_dimension = resolution(1)*abs(max(node_ext_cart(:,1)) - min(node_ext_cart(:,1)));
y_dimension = resolution(2)*abs(max(node_ext_cart(:,2)) - min(node_ext_cart(:,2)));


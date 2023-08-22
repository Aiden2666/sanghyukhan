function [nodal_value] = elemental2nodal (elemental_value, node, element);
%Copyright: Lampros Kourtis, 2007, Stanford University
%email    : kourtis@stanford.edu
num_elements = length(element);
num_nodes    = length(node);

mag          = zeros(num_nodes,1); 
mag_count    = zeros(num_nodes,1);


     for t=1:num_elements                                    % convert elemental values to nodal  
         mag(element(t,:))=mag(element(t,:))+elemental_value(t);   % values by averaging   
         mag_count(element(t,:))=mag_count(element(t,:))+1;  % element values over each node
     end
   
nodal_value=mag./mag_count;
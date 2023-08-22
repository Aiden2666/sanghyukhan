function [YM, v, G, rho, area, centroid, center_mesh, center_mesh_E ] = assign_material (node, node_real, element, A_YM, A_G, A_v, A_d, res, nodes_th);

%% This routine creates for each element a uniformly seeded sampling map in natural (element) coordinates
%  Then using the element shape functions (isoparametric), it maps these
%  points on the image to obtain a series of values of pixels that belong
%  in the space of the element. 
%  Then, the average of these values is selected to be assigned for this
%  particular element. The routine also calculates the centroid of each
%  element and the geometric center of the mesh
%Copyright: Lampros Kourtis, 2007, Stanford University
%email    : kourtis@stanford.edu
%% 

%number of sampling points along parametric edge
num_sample_pts=5;
num_elements = length(element);

aa=0;
xx=zeros(25*length(element),1);
yy=zeros(25*length(element),1);

%Define seed points isoparametric element map (natural coordinates)
[ksi, eta] = meshgrid(-0.9:0.3:0.9,-0.9:0.3:0.9);

%Initialize variables
YM = zeros(num_elements,1);
G  = zeros(num_elements,1);
v  = zeros(num_elements,1);
rho= zeros(num_elements,1);

centroid = zeros(num_elements,2);
area = zeros(num_elements,1);

%% Repeat loop for every element to define centroid, area, and elastic properties of each
%  Use isoparametric map to find correspondance between natural (element)
%  coordinates and actual (real) coordinates of seed points ksi & eta.

for t=1:num_elements
    
%calculate the centroid of each element (bone coordinates)    
centroid(t,1) = mean(node_real(element(t,:),1)); 
centroid(t,2) = mean(node_real(element(t,:),2));

%Define Shape functions and evaluate them at the grid points
N1 = (1/4) * (1-ksi) .* (1-eta);
N2 = (1/4) * (1+ksi) .* (1-eta);
N3 = (1/4) * (1+ksi) .* (1+eta);
N4 = (1/4) * (1-ksi) .* (1+eta);

%Get the isoparametric map. x,y are the mapped variables in actual image coordinates, See Hughes T.J. 1984
x = N1*node(element(t,1),1) + N2*node(element(t,2),1) + N3*node(element(t,3),1) + N4*node(element(t,4),1);
y = N1*node(element(t,1),2) + N2*node(element(t,2),2) + N3*node(element(t,3),2) + N4*node(element(t,4),2);

x=reshape(round(x),length(x)^2,1);
y=reshape(round(y),length(y)^2,1);

%Assign mean values to the elastic properties
%First, Check to see if this is a boundary element 

    YM(t)   = ( trimmean (diag(A_YM(y , x )),20));    %Elastic Modulus     % Note that y and x are inversed 
    G(t)    = ( trimmean (diag(A_G (y , x )),20));    %Shear Modulus       % because of matlab's inverse
    v(t)    = ( mean (diag(A_v (y , x ))));    %Poisson's Ratio     % definition of x,y for plots and images
    rho(t)  = ( trimmean (diag(A_d (y , x )),20));    %Density (rho)
 
    %in case partial voluming correction is off, the external boundary
    %elements are to be treated differently than the internal cause
    %assigning the mean value usually underestimates the assigned values

     if t<=nodes_th      %case where it is a boundary element
       %Elastic Modulus  
     pts     = diag(A_YM(y , x )); mean_ = mean(pts); upper_sample_index = find(pts>mean_);
     if isempty(upper_sample_index)==0
        YM(t)   = mean(pts(upper_sample_index));       
     end
       %Shear Modulus
     pts     = diag(A_G(y , x )); mean_ = mean(pts); upper_sample_index = find(pts>mean_);  
     if isempty(upper_sample_index)==0
        G(t)    = mean(pts(upper_sample_index));;     
     end
      %Density (rho)
     pts     = diag(A_d(y , x )); mean_ = mean(pts); upper_sample_index = find(pts>mean_);  
     if isempty(upper_sample_index)==0
        rho(t)    = mean(pts(upper_sample_index));; 
     end        
      %Poisson's Ratio
     v(t)    = ( max (diag(A_v (y , x ))));    %Poisson's Ratio   
     
%     YM(t)   = ( max (diag(A_YM(y , x ))));    %Elastic Modulus 
%     G(t)    = ( max (diag(A_G (y , x ))));    %Shear Modulus     
%     rho(t)  = ( max (diag(A_d (y , x ))));    %Density (rho) 

     end

    
%Calculate the area of each element by summing the areas of the two
%triangles that from it.
area(t) = 0.5*abs( det ( [  node_real(element(t,1),1) node_real(element(t,1),2) 1;...
                        node_real(element(t,2),1) node_real(element(t,2),2) 1;...
                        node_real(element(t,3),1) node_real(element(t,3),2) 1] )) +...
          0.5*abs( det ( [  node_real(element(t,3),1) node_real(element(t,3),2) 1;...
                        node_real(element(t,4),1) node_real(element(t,4),2) 1;...
                        node_real(element(t,1),1) node_real(element(t,1),2) 1] ));

%Make sure that there are no elements with zero density
 if rho(t)<0.05;rho(t)=0.05;end
 if YM(t)<0.07   ;YM(t) =0.07  ;end
%                     
                    
end

center_mesh(1,1)=sum(centroid(:,1).*area)/sum(area);
center_mesh(1,2)=sum(centroid(:,2).*area)/sum(area);

center_mesh_E(1,1)=sum(centroid(:,1).*area.*YM)/sum(area.*YM);
center_mesh_E(1,2)=sum(centroid(:,2).*area.*YM)/sum(area.*YM);

%Extra Correction for partial voluming. The first-to-external ring of elements gets
%assigned to the external ones (in case they were of lower value)
 
     ext_elements = 1:nodes_th;
     first_ring_elements = nodes_th+1:2*nodes_th;
     YM(1:nodes_th) = max([YM(ext_elements) YM(first_ring_elements)]');
     G(1:nodes_th)  = max([G(ext_elements) G(first_ring_elements)]');
     rho(1:nodes_th)  = max([rho(ext_elements) rho(first_ring_elements)]');

    
 

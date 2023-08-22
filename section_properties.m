function [I, J, I_principal, V,   I_E, J_E, I_E_principal, V_E,   I_G, J_G, I_G_principal, V_G] = section_properties (YM, G, rho, centroid, area, center_mesh, center_mesh_E)
%Copyright: Lampros Kourtis, 2007, Stanford University
%email    : kourtis@stanford.edu                                                                                                                     


Ix=0;Iy=0;Ix_J=0;Iy_J=0;Ixy=0;
Ix_E=0;Iy_E=0;Ixy_E=0;
Ix_G=0;Iy_G=0;Ixy_G=0;

%Sectional props calculation: ATTENTION:Need to add the element's own Ix. For example
%for a rectangular : bh^3/12. Split up elements in triangles and add them.
for t=1:length(rho)
    Ix  = Ix   + ((centroid(t,2)-center_mesh(2))^2)*area(t);
    Iy  = Iy   + ((centroid(t,1)-center_mesh(1))^2)*area(t);
    Ixy = Ixy  + ((centroid(t,2)-center_mesh(2))*(centroid(t,1)-center_mesh(1)))*area(t);
    
    %Elastic modulus weighted 
    Ix_E  = Ix_E   + ((centroid(t,2)-center_mesh_E(2))^2)*area(t)*YM(t);
    Iy_E  = Iy_E   + ((centroid(t,1)-center_mesh_E(1))^2)*area(t)*YM(t);
    Ixy_E = Ixy_E  + ((centroid(t,2)-center_mesh_E(2))*(centroid(t,1)-center_mesh_E(1)))*area(t)*YM(t);
    
    %(G=mu) weighted J 
    Ix_G  = Ix_G   + ((centroid(t,2)-center_mesh_E(2))^2)*area(t)*G(t);
    Iy_G  = Iy_G   + ((centroid(t,1)-center_mesh_E(1))^2)*area(t)*G(t);
    Ixy_G = Ixy_G  + ((centroid(t,2)-center_mesh_E(2))*(centroid(t,1)-center_mesh_E(1)))*area(t)*G(t);
    
end

% Traditional homogenous section properties (not accounting for inhomogeneity)
I   =   [Ix Ixy;Ixy Iy];
J   =   Ix   +  Iy;
[V,I12]=eig(I);I12 =   diag(I12);
I_principal = [max(I12) 0; min(I12) 0];
if I12(1)<I12(2);   V = [V(:,2) V(:,1)];         end % sort axes

% Elastic modulus weighted section properties
I_E =  [Ix_E Ixy_E;Ixy_E Iy_E]; 
%I_E =  sum(area) * I_E / sum(area.*YM)
J_E =  ( Ix_E + Iy_E ) ;
[V_E,I12]=eig(I_E);I12 =   diag(I12);
I_E_principal = [max(I12) 0; min(I12) 0];
if I12(1)<I12(2);   V_E = [V_E(:,2) V_E(:,1)];   end % sort axes

% Shear modulus weighted section properties
I_G =  [Ix_G Ixy_G;Ixy_G Iy_G]; 
%I_G =  sum(area) * I_G / sum(area.*G);
J_G =  (Ix_G + Iy_G) ;
[V_G,I12]=eig(I_G);I12 =   diag(I12);
I_G_principal = [max(I12) 0; min(I12) 0];
if I12(1)<I12(2);   V_G = [V_G(:,2) V_G(:,1)];   end % sort axes

I
J
I_principal
V

I_E
I_E_principal
V_E


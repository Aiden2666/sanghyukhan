function [bending_stress] = calculate_bending (M_x, M_y, Q_x, Q_y, z, P, I, node, area,centroid, YM, center)
%Copyright: Lampros Kourtis, 2007, Stanford University
%email    : kourtis@stanford.edu     

% calculate the bending for a given moment of inertia tensor I
% I can be either the geometric (I) or the modulus weighted tensor (I_E)
% from main.m

num_nodes=length(node);
num_elements = length(area);
bending_stress = zeros (num_elements,1);

Ix  = I(1,1);
Iy  = I(2,2);
Ixy = I(2,1);

x = centroid(:,1);
y = centroid(:,2);

A = sum(area.*YM);

if or(M_x~=0,M_y~=0)
    bending_stress = - ( (M_y*Ix + M_x*Ixy) / (Ix*Iy-Ixy^2) ) .* YM .*x + ( (M_x*Iy + M_y*Ixy) / (Ix*Iy-Ixy^2) ) .* YM .*y ;
end

if or(Q_x~=0,Q_y~=0)
    bending_stress = bending_stress + ( ( (Q_x*Ix + Q_y*Ixy) / (Ix*Iy-Ixy^2) ) .* YM .*x + ( (Q_y*Iy + Q_x*Ixy) / (Ix*Iy-Ixy^2) ) .* YM .*y ) * z ;
end

axial_stress   = ( P / A ) * YM;

bending_stress = bending_stress + axial_stress;



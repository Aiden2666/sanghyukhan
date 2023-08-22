function [center] = find_center (img)
%Copyright: Lampros Kourtis, 2007, Stanford University
%email    : kourtis@stanford.edu
[size_x,size_y]=size(img);

%Set the mid point as the initial center for the estimation of the actual center
initial_x_center = size_x/2;
initial_y_center = size_y/2;

% %Create blob
% bin_image = im2bw(img,0);
% %Fill holes in blobs
% bin_image = imfill(bin_image,'holes');
% %Select bone blob from all others
% bin_image = bwselect(bin_image,initial_x_center, initial_y_center);

%% Weighted center calculation 
% To make calculations faster, use matrix operations instead of exhaustive
% pixel to pixel calculations

% Create gradient masks for the 2 directions
kernel=linspace(1,size_x, size_x);
y_kernel=kernel'*ones(1,size_y);

kernel=linspace(1,size_y, size_y);
x_kernel=(kernel'*ones(size_x,1)')';


% Multiply by the image to find gradient based image in x and y directions
x_image = x_kernel.*img;
y_image = y_kernel.*img;

% Simply find the mean of the gradient image which is the center
center_x = sum(sum(x_image))/sum(sum(img));
center_y = sum(sum(y_image))/sum(sum(img));

center=[center_x, center_y];
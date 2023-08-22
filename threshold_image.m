function [A, A_bin] = threshold_image (A,segmentation_threshold)
%Copyright: Lampros Kourtis, 2007, Stanford University
%email    : kourtis@stanford.edu
%% Threshold image based on argument value

% A is the original HU image
% A_d is the density image
% A_YM is the Young's Modulus image
% A_v is the Poisson's ratio image

[size_x,size_y]=size(A);

% threshold the image. im2bw needs a 0 to 1 value, therefore the conversion
% first create a binary image bin_A
A_bin = im2bw(A/max(max(A)),segmentation_threshold/max(max(A)));

% now apply the thresholding to the original image. A remains in HU
A = A_bin.*A;

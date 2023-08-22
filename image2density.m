function [A_d] = image2density (A, bin_A)
%Copyright: Lampros Kourtis, 2007, Stanford University
%email    : kourtis@stanford.edu

%% image to density map

% A is the original HU image
% A_d is the density image
[size_x,size_y]=size(A);

% convert image values to apparent wet density (gr/cm3)
% density is given using a y=ax+b relation  for which:

%% Input
A=A/1000;
a=2.2092;
b=-0.5852;

%% Convert

A_d=bin_A.*abs(A*a + b*ones(size_x,size_y)); %in gr/cm3


% A_d now has bone density values in gr/cm3 (no water & soft tissue present)

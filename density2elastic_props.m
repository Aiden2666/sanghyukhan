function [A_YM, A_G, A_v] = density2elastic_props (A_d)
%Copyright: Lampros Kourtis, 2007, Stanford University
%email    : kourtis@stanford.edu

%% Assign material properties YM, v based on density based on Keller 1994

%%
% A_d is the density image
[size_x,size_y]=size(A_d);


%% Need some input, the two exponents and the first linear term
   exp1=3.46;
   lin_term1=1990;

%% For each pixel, define properties

     A_YM = lin_term1 * A_d.^exp1 ;
     
% Here the shear modulus is given as a percentage of the Elastic modulus as
% defined in Reilly & Burstein 1976? for the long axis of the bone
     A_G = 0.19 * A_YM;

     A_v = 0.20 * (A_d./A_d);
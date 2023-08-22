function [node_ext_polar, node_ext_cart, node_int_polar, node_int_cart] = detect_boundary ( center, nodes_th, canal_flag, ext_blob, int_blob)
%Copyright: Lampros Kourtis, 2007, Stanford University
%email    : kourtis@stanford.edu

%% External bone contour detection

%Get binary perimeter (boundary) map
boundary = bwperim (ext_blob);

%Get location of boundary pixels
[y,x,v] = find (boundary);

%Translate so that x,y have an origin at the center of the image (in order
%to convert to polar coordinates
 x=x-center(1);
 y=y-center(2);
 
%Convert to polar coordinates
[th,r]=cart2pol(x,y);

%Sort points in polar th direction
[th,order]=sort(th);
r=r(order); %reorder r coords accordingly
%Now r,th are ordered around the clock

%[th,r] = reducem(th,r,0.5) %Activate in case of large (unwanted) number of
                                %contour points. It downsamples the contour

%Perform smoothing spline operation
pp=csaps(th,r,0.999,[],repmat(30,length(r),1));

    %Define (re)sampling points for spline evaluation 
 pt_step = 2*pi / nodes_th;
 theta = -pi:pt_step:pi;

    %Evaluate Spline
radius=fnval(pp,theta);
 %radius=radius-0.5;      %Apply correction for binary image pixel roundoff

[x,y]=pol2cart(theta,radius); %the pi/2-  is there cause matlab is treating x,y in an inverse way
 x=x+center(1);      %Apply correction for center offset (again because 
 y=y+center(2);     %of binary image pixel roundoff
  
 node_ext_polar = [theta',  radius'];  %Assign values to 
 node_ext_cart  = [x'    ,  y'     ];            %output arguments

 node_ext_polar(end,:) = node_ext_polar(1,:);   %Circle points by substituting 
 node_ext_cart (end,:) = node_ext_cart (1,:);   %the last point with the first

%% Intramedulary canal boundary detection
if canal_flag==1  %case where there is an intramedulary canal
    
%Get binary perimeter (boundary) map
boundary_int = bwperim (int_blob);

%Get location of boundary pixels
[y_int,x_int,v] = find (boundary_int);
%Translate so that x2,y2 have an origin at the center of the image (in order
%to convert to polar coordinates
 x_int=x_int-center(1);
 y_int=y_int-center(2);
%Convert to polar coordinates
[th_int,r_int]=cart2pol(x_int,y_int);

%Sort points in polar th direction
[th_int,order]=sort(th_int);
r_int=r_int(order);                         %reorder r coords accordingly
%Now r2,th2 are ordered around the clock

%[th_int,r_int] = reducem(th_int,r_int,0.5) %Activate in case of large (unwanted) number of
                                            %contour points. It downsamples the contour

%Perform smoothing spline operation
pp_int=csaps(th_int,r_int,0.9,[],repmat(30,length(r_int),1));

    %Define (re)sampling points for spline evaluation 
 pt_step = 2*pi / nodes_th;
 theta = -pi:pt_step:pi;

    %Evaluate Spline
radius_int=fnval(pp_int,theta);
 radius_int=radius_int+0.5;                 %Apply correction for binary image pixel roundoff

 
[x_int,y_int]=pol2cart(theta,radius_int); %the pi/2-  is there cause matlab is treating x,y in an inverse way
 x_int=x_int+center(1);                   %Apply correction for center offset (again because 
 y_int=y_int+center(2);                 %of binary image pixel roundoff
        
 node_int_polar = [theta' ,radius_int'];    %Assign values to 
 node_int_cart  = [x_int' ,y_int'     ];    %output arguments
 
 node_int_polar(end,:) = node_int_polar(1,:);   %Circle points by substituting 
 node_int_cart (end,:) = node_int_cart (1,:);   %the last point with the first
 
else                                        %case where there is an intramedulary canal
    
 node_int_polar = [];                       %Assign values to 
 node_int_cart  = [];                       %output arguments
end

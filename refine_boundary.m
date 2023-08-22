function [node] = refine_boundary (node, A, center_A, z, line_handle)
%Copyright: Lampros Kourtis, 2007, Stanford University
%email    : kourtis@stanford.edu

%node   is given and output in polar coordinates
%z      is the zone width in pixels. The existing boundary will serve as the midline.
%A      is the image
%center_A = [cx, cy] is the center of image A

num_pts = length(node);
sampling_pts = z;

%% Start radial search by acquiring the profile of the image along the search
%zone width for the internal and the external boundary

%initialize optimization parameters
ops = optimset('Display','off','MaxIter',200,'MaxSQPiter',10,'TolFun',5e-2);
ops.LargeScale = 'off';
A_ = []; B_ = [];Aeq = []; Beq = [];LB = [ 1000 0.05  1]; UB = [max(max(A)) 4 z-1];

h = waitbar(0,'Boundary Refinement');

for t=1:num_pts
    waitbar(t/num_pts,h)
    z_start = polar2cart ([node(t,1), node(t,2)+z/2], center_A);   %starting (outer) point of zone search
    z_end   = polar2cart ([node(t,1), node(t,2)-z/2], center_A);   %ending   (inner) point of zone search    
    % Take the profile of the image intensity along the path
    [cx,cy,zone_profile] = improfile(A, [ z_start(1)  z_end(1) ], [ z_start(2)  z_end(2) ],sampling_pts,'bicubic');
    % find the radial distance of each of these points
    path=[cx,cy]; path_dist=(sum((path-repmat(z_start,length(path),1))'.^2).^0.5)';
    % fit a gaussian distribution error function (erf) to the profile data :
    % 0.5*A*(1+erf((x-mu)/(sigma*sqrt(2))));
    if t==1
        x0 = [ zone_profile(end) 1.6 sampling_pts/2 ];  % first time run intial guesses [A   sigma  mu]
    else
        x0 = [ x(1) x(2) x(3) ];  % now apply previously found coefs to help optimization      
    end

    [x, fval, exitflag] = fmincon('fun', x0, A_, B_, Aeq, Beq, LB, UB, [], ops, zone_profile, path_dist);
    node(t,2)=node(t,2)+z/2-x(3);
  
node_ext_cart = polar2cart(node,center_A);
set(line_handle,'XData',node_ext_cart(1:end-1,1),'YData',node_ext_cart(1:end-1,2));
drawnow
end

    close(h)
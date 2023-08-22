function [r, center,error] = Calculate_Curvature (x,y);
%Copyright: Lampros Kourtis, 2007, Stanford University
%email    : kourtis@stanford.edu
% flatten arrays into n x 1 vectors.
  x = x(:);
  y = y(:);
 
 n=length(x);
  
% compute centroid.
  xbar = mean(x);
  ybar = mean(y);
  
% shift centroid to origin.
  xrel = x - xbar;
  yrel = y - ybar;
  
% scaling.
  if 1
    sc = mean(sqrt(xrel.^2 + yrel.^2));
    xrel = xrel / sc;
    yrel = yrel / sc;
  else
    sc = 1.0;
  end
  
% make n x 4 'design matrix' D.
  D = [xrel.^2+yrel.^2   -2*xrel   -2*yrel   ones(n, 1)];
  
% find approximate null-vector
  [U, S, V] = svd(D);
  nv = V(:, end);
  
  if nv(1)
    % dividing nv by nv(1) gives [1; *; *; *]
    a = nv(2) / nv(1);
    b = nv(3) / nv(1);
    r2 = a^2 + b^2 - nv(4) / nv(1);
    if r2 >= 0
      % OK.
      r = sqrt(r2);
    else
      % no good.
      a = nan;
      b = nan;
      r = nan;
    end
  else
    % no good. probably caused by all points lying on a line.
    a = nan;
    b = nan;
    r = nan;
  end
  
  % adjust for shifted centroid, and scaling.
  a = sc * a + xbar;
  b = sc * b + ybar;
  r = sc * r;

  center=[a,b];
  error = abs((r- sqrt(   sum   (  (center-[x(round(median(1:n))), y(round(median(1:n)))] ).^2)    )   )/r); %eucledian distance
  
end


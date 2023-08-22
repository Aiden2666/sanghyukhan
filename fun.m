function [d] = fun(coef, y, x)



myest = 0.5*coef(1)*(1+erf((x'-coef(3))/(coef(2)*sqrt(2))));

d = sqrt( sum( (myest - y').^2 ) );

% A       = coef(1);
% sigma   = coef(2);
% mu      = coef(3);
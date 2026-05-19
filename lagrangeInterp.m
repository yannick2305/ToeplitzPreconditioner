function v = lagrangeInterp(u,x,y)
%
% Given distinct points x0, x1, ..., xn and values
% y0, y1, ..., yn, and u,
% return value v of the polynomial interpolant of degree <= n
% to the data (x_i,y_i), i=0,1,...,n, at u.
%
% Inputs:
% x  : [x_0, x_1, ..., x_n], interpolation points
% y  : [y_0, y_1, ..., y_n], data values
% u  : point where interpolant should be evaluated.
%
% Outputs:
% v  : value of interpolant.
%
% Author: Michael Floater
% e-mail: michaelf@math.uio.no
% Created: May 2026, using Matlab R2024b
% Copyright 2026

n = length(x);
lag = zeros(1,n);
for i=1:n
  lag(i) = 1.0;
  for j=1:n
    if j ~= i
      lag(i) = lag(i) * (u-x(j))/(x(i)-x(j));
    end
  end
end

v = lag * y';
 

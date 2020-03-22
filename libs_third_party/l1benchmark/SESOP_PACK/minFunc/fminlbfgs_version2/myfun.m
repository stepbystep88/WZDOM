function [f,g] = myfun(x)
% This function MYFUN is used for demonstrating the optimizer
% on a quadratic function.
%
% [f,g] = myfun(x)
%
%  x : input parameters of equation length 5
%  f : equation value
%  g : equation gradient
%
% example,
%  [f,g] = myfun([1 1 1 1 1]);
%
f = 3*x(1)^2 + 2*x(1)*x(2) + x(2)^2 + x(3)^2 + x(4)^2 + x(5)^4;
if nargout > 1
   g(1) = 6*x(1)+2*x(2);
   g(2) = 2*x(1)+2*x(2);
   g(3) = 2*x(3);
   g(4) = 2*x(4);
   g(5) = 4*x(5)^3;
end

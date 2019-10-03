function [f, g] = bsShiftedRastrigin(x, isGradient, o)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------
% Shifted Rastrigin Function 
% Programmed by Bin She (bin.stepbystep@gmail.com)
% Programming dates: May 2019   
% -------------------------------------------------------------------------
% For function details and reference information, see:
% see http://benchmarkfcns.xyz/benchmarkfcns/rastriginfcn.html
%
% Global minimum: min(f) = 0 at x = o'
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Input
%
% x: it must be a column vector representing the variable to be estimated
%
% isGradient: a logical variabel. The function will calculate gradient
% information of this objective function if it is equal to true, otherwise,
% the subroutine will not compute the gradient information. 
% -------------------------------------------------------------------------
% Output
%
% f: objective function value.
% 
% g: gradient information. g = [] if isGradient = 0.
% -------------------------------------------------------------------------

    if nargin == 1 || isempty(isGradient)
        isGradient = true;
    end
    
    
    [n, m] = size(x);
    o = repmat(o, 1, m);
    A = 10;
    z = x - o;
    
    pi2z = 2 * pi * z;
    
    z2 = z .^ 2;
    f = A * n + sum( z2 - 10 * cos(pi2z), 1);
    
    if isGradient
        g = 2*z + 20 * pi * sin(pi2z);
    else
        g = [];
    end
end
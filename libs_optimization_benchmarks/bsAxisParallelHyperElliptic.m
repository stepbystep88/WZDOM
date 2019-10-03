function [f, g] = bsAxisParallelHyperElliptic(x, isGradient)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------
% Axis parallel hyper-ellipsoid function
% Programmed by Bin She (bin.stepbystep@gmail.com)
% Programming dates: June 2019   
% -------------------------------------------------------------------------
% For function details and reference information, see:
% see http://www.geatbx.com/ver_3_5/fcnfun1a.html
%
% Global minimum: min(f) = 0 at x = (0, 0, ...)'; -600<=x(i)<= 600.
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
    
    seq = (1 : n)';
    seq = repmat(seq, 1, m);
    
    f = sum(seq .* x .^ 2, 1);
    
    if isGradient
        
        g = 2 * seq .* x;
    else
        g = [];
    end
end
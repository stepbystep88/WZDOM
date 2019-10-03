function [f, g] = bsShiftedSchwefel(z, isGradient, o)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------
% Schwefel Function
% Programmed by Bin She (bin.stepbystep@gmail.com)
% Programming dates: June 2019   
% -------------------------------------------------------------------------
% For function details and reference information, see:
% see http://benchmarkfcns.xyz/benchmarkfcns/schwefelfcn.html
%
% Global minimum: min(f) = 0 at x = (0, 0, ...)'; -500<=x(i)<= 500.
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

    
    [n, m] = size(z);
    o = repmat(o, 1, m);
    
    x = z - o;
    
    tmp = sqrt(abs(x));
    f = 418.9829 * n - sum(x.*sin(tmp), 1);
    
    if isGradient
        g = -sin(tmp) - (x.*cos(tmp).*sign(x))./(2*tmp); 
    else
        g = [];
    end
end
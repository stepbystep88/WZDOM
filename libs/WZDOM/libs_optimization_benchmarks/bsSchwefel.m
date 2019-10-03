function [f, g] = bsSchwefel(x, isGradient)
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

    
    n = size(x, 1);
    
    tmp = sqrt(abs(x));
    
    a = 4.209687462275036e+002;
    b = 418.9829 - a*sin(sqrt(a));
    
    f = (418.9829 - b) * n - sum(x.*sin(tmp), 1); 
    
    if isGradient
        g = -sin(tmp) - (x.*cos(tmp).*sign(x))./(2*tmp); 
    else
        g = [];
    end
end
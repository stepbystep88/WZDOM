function [f, g] = bsDropWave(z, isGradient)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------
% bsDropWave Function 
% Programmed by Bin She (bin.stepbystep@gmail.com)
% Programming dates: Sep 2019   
% -------------------------------------------------------------------------
% For function details and reference information, see:
% see https://www.sfu.ca/~ssurjano/Code/dropm.html
%
% Global minimum: min(f) = 0 at x = (0, 0, ...)'
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


    n = size(z, 1);
    assert(n == 2, 'Cross-Leg Table function is only defined on a 2D space.')
    x = z(1, :);
    y = z(2, :);
    
    frac1 = 1 + cos(12*sqrt(x.^2+y.^2));
    frac2 = 0.5*(x.^2+y.^2) + 2;

    f = -frac1./frac2;
    
    if isGradient
        
        tmp1 = cos(12*(x.^2 + y.^2).^(1/2));
        tmp2 = (x.^2/2 + y.^2/2 + 2).^2;
        tmp3 = sin(12*(x.^2 + y.^2).^(1/2));
        tmp4 = ((x.^2 + y.^2).^(1/2).*(x.^2/2 + y.^2/2 + 2));
        
        g1 = (x.*(tmp1 + 1))./tmp2 + (12*x.*tmp3)./tmp4;
        g2 = (y.*(tmp1 + 1))./tmp2 + (12*y.*tmp3)./tmp4;
        g = [g1; g2];
    else
        g = [];
    end
end
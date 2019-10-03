function [f, g] = bsCrossLegTable(x, isGradient)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------
% Cross-Leg Table function
% Programmed by Bin She (bin.stepbystep@gmail.com)
% Programming dates: May 2019   
% -------------------------------------------------------------------------
% For function details and reference information, see:
% https://al-roomi.org/benchmarks/unconstrained/2-dimensions/45-cross-leg-table-function
%
% Global minimum: min(f) = -1 at x = (0, 0)'
%
% Booth's function is only defined on a 2D space.
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
    assert(n == 2, 'Cross-Leg Table function is only defined on a 2D space.')
    X = x(1, :);
    Y = x(2, :);
    sinsin = sin(X) .* sin(Y);
    sigma2 = sqrt(X.^2 + Y.^2);
    sigma1 = exp( abs(sigma2/pi - 100));
    
    f = -(abs(sinsin .* sigma1) + 1) .^ -0.1;
    
    if isGradient
        tmp1 = (sigma1 .* abs(sinsin) + 1).^1.1;
        tmp2 = sign(sigma2./pi-100).*abs(sinsin)./(pi.*sigma2);
        
        g1 = 0.1*sigma1.*(sign(sinsin).*cos(X).*sin(Y) + (X.*tmp2))./tmp1;
        g2 = 0.1*sigma1.*(sign(sinsin).*cos(Y).*sin(X) + (Y.*tmp2))./tmp1;
        g = [g1; g2];

 
    else
        g = [];
    end
end
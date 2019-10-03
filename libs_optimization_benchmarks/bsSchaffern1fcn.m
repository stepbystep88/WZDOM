function [f, g] = bsSchaffern1fcn(x, isGradient)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------
% Schaffer's function
% Programmed by Bin She (bin.stepbystep@gmail.com)
% Programming dates: June 2019   
% -------------------------------------------------------------------------
% For function details and reference information, see:
% http://benchmarkfcns.xyz/benchmarkfcns/schaffern1fcn.html
%
% Global minimum: min(f) = 0 at x = (0, 0)'
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
    assert(n == 2, 'Booth''s function is only defined on a 2D space.')
    X = x(1, :);
    Y = x(2, :);
    
    tmp0 = (X .^ 2 + Y .^ 2) .^ 2;
    
    numeratorcomp = (sin(tmp0) .^ 2) - 0.5; 
    tmp1 = 1 + 0.001 * (X .^2 + Y .^2);
    denominatorcomp = tmp1 .^2 ;
    f = 0.5 + numeratorcomp ./ denominatorcomp;
    tmp2 = cos(tmp0) .* sin(tmp0) .* (X.^2 + Y.^2);
    
    if isGradient
        g1 = (8*X.*tmp2)./tmp1.^2 - (X.*numeratorcomp)./(250*tmp1.^3);
        g2 = (8*Y.*tmp2)./tmp1.^2 - (Y.*numeratorcomp)./(250*tmp1.^3);
        g = [g1; g2];
    else
        g = [];
    end
end
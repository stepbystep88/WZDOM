function [f, g] = bsPeak(z, isGradient)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------
% MATLAB PEAK Function 
% Programmed by Bin She (bin.stepbystep@gmail.com)
% Programming dates: Sep 2019   
% -------------------------------------------------------------------------
% For function details and reference information, see:
% see http://benchmarkfcns.xyz/benchmarkfcns/ackleyfcn
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
    
    f = 3*(1-x).^2.*exp(-(x.^2) - (y+1).^2) ... 
   - 10*(x/5 - x.^3 - y.^5).*exp(-x.^2-y.^2) ... 
   - 1/3*exp(-(x+1).^2 - y.^2);
    
    if isGradient
        
        tmp1 = exp(- (x + 1).^2 - y.^2);
        tmp2 = exp(- (y + 1).^2 - x.^2);
        tmp3 = exp(- x.^2 - y.^2);
        
        g1 = (tmp1.*(2*x + 2))/3 + 3*tmp2.*(2*x - 2) + tmp3.*(30*x.^2 - 2) - 6*x.*tmp2.*(x - 1).^2 - 2*x.*tmp3.*(10*x.^3 - 2*x + 10*y.^5);
        g2 = (2*y.*tmp1)/3 + 50*y.^4.*tmp3 - 3*tmp2.*(2*y + 2).*(x - 1).^2 - 2*y.*tmp3.*(10*x.^3 - 2*x + 10*y.^5);
        g = [g1; g2];
    else
        g = [];
    end
end
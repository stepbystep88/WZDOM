function [f, g] = bsAckley(x, isGradient)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------
% Ackley Function 
% Programmed by Bin She (bin.stepbystep@gmail.com)
% Programming dates: May 2019   
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


    n = size(x, 1);
    
    sum1 = sum(x .^ 2, 1);
    sum2 = sum(cos(2 * pi * x), 1);
    sqrtTerm = sqrt(1/n*sum1);
    expSqrtTerm = exp(-0.2*sqrtTerm);
    expCosTerm = exp(1/n*sum2);
    f = -20*expSqrtTerm - expCosTerm + exp(1) + 20;
    
    if isGradient
        
        g1 = (4/n./sqrtTerm) .* (expSqrtTerm  .* x);
        g2 = pi * (expCosTerm .* sin(2 * pi * x));
        g = g1 + g2;
    else
        g = [];
    end
end
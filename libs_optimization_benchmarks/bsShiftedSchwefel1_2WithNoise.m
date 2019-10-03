function [f, g] = bsShiftedSchwefel1_2WithNoise(x, isGradient, data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------
% Shifted Schwefel’s Problem 1.2 with Noise in Fitness
% Programmed by Bin She (bin.stepbystep@gmail.com)
% Programming dates: June 2019   
% -------------------------------------------------------------------------
% For function details and reference information, see:
% see http://benchmarkfcns.xyz/benchmarkfcns/schwefel222fcn.html
%
% Global minimum: min(f) = 0 at x = o; -600<=x(i)<= 600.
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
    sigma = repmat(data(1), 1, m);
    o = repmat(data(2:end), 1, m);
    
    z = x - o;
    
    sumzid = zeros(n, m);
    for i = 1 : n
        sumzid(i, :) = sum(z(1:i, :), 1);
    end
    
    
    f = sum(sumzid .^ 2, 1) .* (1 + 0.4*abs(sigma));
    
    
    
    if isGradient
        g = zeros(n, m);
        for i = 1 : n
            g(i, :) = 2 * sum(sumzid(i:n, :), 1);
        end
        g = (1 + 0.4*abs(sigma)) .* g; 
    else
        g = [];
    end
end
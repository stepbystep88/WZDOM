function [f, g] = bsStochasticRosenbrock(x, isGradient, sigma)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------
% StochasticRosenbrock Function 
% Programmed by Bin She (bin.stepbystep@gmail.com)
% Programming dates: May 2019   
% -------------------------------------------------------------------------
% For function details and reference information, see:
% see https://www.sfu.ca/~ssurjano/rosen.html
%
% Global minimum: min(f) = 0 at x = (1, 1, ...)'
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
    
    x2 = x.^2;
    
    xi = x(1:n-1, :);
    xip1 = x(2:n, :);
    xi_2 = x2(1:n-1, :);
    
%     sigma = repmat(n, m);
    
    sigma = repmat(sigma, 1, m);
%     sigma = rand(n, m);
    
    f = sum ( 100 * (sigma(1:n-1, :) .* (xip1 - xi_2).^2) +  (xi - 1).^2, 1);
%     f = f';
    
    if isGradient
        g1 = 2*(x(1, :) - 1) - 400*sigma(1, :).*x(1, :).*(x(2, :) - x(1, :).^2);
        gm = 200*sigma(n-1, :).*(x(n, :) - x(n-1, :).^2);

        xk = x(2:n-1, :);
        xkm1_2 = x2(1:n-2, :);
        xkp1 = x(3:n, :);
        xk_2 = x2(2:n-1, :);

        gk = 200.*sigma(2:n-1, :).*(xk-xkm1_2) + 2*xk - 400.*sigma(2:n-1, :).*(xk.*(xkp1 - xk_2)) - 2;

        g = [g1; gk; gm];
    else
        g = [];
    end
    
end
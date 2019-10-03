function [f, g] = bsCFSingle(x, isGradient, o)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------
% Griewank Function 
% Programmed by Bin She (bin.stepbystep@gmail.com)
% Programming dates: May 2019   
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
    f = zeros(1, m);
    g = zeros(size(x));
    
    for i = 1 : 3
        switch i
            case 1
                [tf, gf] = bsRosenbrock(x+1, isGradient);
            case 2
                [tf, gf] = bsSphere(x, isGradient);
            case 3
                [tf, gf] = bsSchwefel2_22(x, isGradient);
%             case 2
%                 [tf, gf] = bsStochasticRosenbrock(x, isGradient, o(3:3+n-1));
%             case 3
%                 [tf, gf] = bsStochasticRosenbrock(x, isGradient, o(3+n:end));    
        end
        
        f = f + o(i) * tf;
        if isGradient
            g = g + o(i) * gf;
        end
    end
end
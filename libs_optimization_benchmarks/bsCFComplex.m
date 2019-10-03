function [f, g] = bsCFComplex(x, isGradient, o)
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

    n = size(x, 1);
    
    f = zeros(1, size(x, 2));
    g = zeros(size(x));
    
    for i = 1 : 4
        switch i
            case 1
                [tf, gf] = bsSphere(x, isGradient);
            case 2
                [tf, gf] = bsXinSheYang(x, isGradient);
            case 3
                [tf, gf] = bsAckley(x, isGradient);
            case 4
                [tf, gf] = bsRastrigin(x, isGradient);
            case 5
                [tf, gf] = Griewank(x, isGradient);
%             case 6
%                 [tf, gf] = bsShiftedRastrigin(x, isGradient, o(7+n:end));   
        end
        
        f = f + o(i) * tf;
        if isGradient
            g = g + o(i) * gf;
        end
    end
end
function [f, g, data] = bsReg1DZeroOrderTK(x, data, isInitial)
%% Regularization function, return the value and gradient of function $|x|_2^2$
% Bin She, bin.stepbystep@gmail.com, October, 2019
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
%
% x             is a column vector; refers to the parameter to be estimated.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    f = x' * x;
    
    g = 2 * x;
end
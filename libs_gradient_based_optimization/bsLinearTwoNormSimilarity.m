function [f, g] = bsLinearTwoNormSimilarity(x, data)
%% calculate the value of function |Ax - cB|_2^2 and its gradient as well 
% only compare the similarity of two items
% Bin She, bin.stepbystep@gmail.com, February, 2019
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
%
% x                         is a column vector; refers to the parameter to be estimated.
%
% data.A               is a matrix; refers to the kernel/deblur/sampling matrix.
%
% data.B               is a column vector; refers to the observed data.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT
% f                     	is a scalar; refers to |Ax - B|_2^2.
%
% g                         is a column vector; refers to the gradient of
% |Ax - B|_2^2 with respect to x.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    d1 = data.A * x;
    d2 = data.B;
    d1d1 = d1' * d1;
    
    c = (d1'*d2) / d1d1;
    
    z = c*d1 - d2;

    f = sum(z.^2, 1);
    
    g1 = 2*c*z;
    g2 = 2*(d1'*z)/d1d1*(-2*c*d1 + d2);
%     g = 2*(c - (d1'*z)/d1d1)*(data.A' * z);
    g = data.A' * (g1 + g2);
end
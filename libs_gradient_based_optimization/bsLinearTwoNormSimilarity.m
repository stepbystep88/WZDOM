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
    c = 1 / norm(d1);
    
    % method 1
    z = c*d1 - d2;
    f = sum(z.^2, 1);    
    g = 2*data.A'*(c*z - c^3*d1*(d1'*z));
    
    
end
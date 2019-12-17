function [f, g] = bsLinearCorrelation(m, data)
%% calculate the correlation of Ax, B and its gradient as well 
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
    
    x = data.A * m;
    y = data.B;
    
    mxy = mean(x .* y);
    mx = mean(x);
    my = mean(y);
    mx2 = mean(x.^2);
    my2 = mean(y.^2);
    
    f = (mxy - mx.*my)/(sqrt(mx2-mx^2)*sqrt(my2-my^2));
    
    g = 2*c*(data.A' * z) ;
end
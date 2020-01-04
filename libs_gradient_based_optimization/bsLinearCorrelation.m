function [f, g] = bsLinearCorrelation(x, data)
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
    
%     x = data.A * m;
%     y = data.B;
    
%     mxy = mean(x .* y);
%     mx = mean(x);
%     my = mean(y);
%     mx2 = mean(x.^2);
%     my2 = mean(y.^2);
%     
%     f = (mxy - mx.*my)/(sqrt(mx2-mx^2)*sqrt(my2-my^2));
    
    
    f = fcn(x, data);
    fcn_ = @(x)fcn(x, data);
    g = mb_numDiff(fcn_, x, 0.000001);
%     g1 = 2*c*z;
%     g2 = 2*(d1'*z)/d1d1*(-2*c*d1 + d2);
% %     g = 2*(c - (d1'*z)/d1d1)*(data.A' * z);
%     g = data.A' * (g1 + g2);
    
%     g = 2*c*(data.A' * z) ;
end

function f = fcn(x, data)
    d1 = data.A * x;
    d2 = data.B;
    
    c = bsComputeGain(d2, d1);
%     c = d1'*d2/(d1'*d1);
%     c = norm(d2)/norm(d1);
    z = c*d1 - d2;

    f = sum(z.^2, 1);
end
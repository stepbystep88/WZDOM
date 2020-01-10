function [f, g] = bsXinSheYang(x, isGradient, c)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------
% Xin-She Yang Function Function
% Programmed by Bin She (bin.stepbystep@gmail.com)
% Programming dates: June 2019   
% -------------------------------------------------------------------------
% For function details and reference information, see:
% see http://benchmarkfcns.xyz/benchmarkfcns/xinsheyangn1fcn.html
%
% Global minimum: min(f) = 0 at x = (0, 0, ...)'; -600<=x(i)<= 600.
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
    
    absx = abs(x);
    seq = (1 : n)';
    seq = repmat(seq, 1, m);
    
    if length(c) == 1
        o = rand(n, m) * c;
    else
        o = repmat(c, 1, m);
    end
    
    tmp = o .* absx .^ (seq - 1);
    
    f = sum(tmp .* absx, 1);
    
    if isGradient
        g = sign(x) .* tmp .* seq;
    else
        g = [];
    end
end
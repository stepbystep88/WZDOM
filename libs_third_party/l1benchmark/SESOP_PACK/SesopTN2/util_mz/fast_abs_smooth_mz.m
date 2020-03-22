function [f,df,d2f]=fast_abs_smooth_mz(x,eps)
% Smoothed approximation of absolute value function and its 1st and 2nd derivatives.  
% Appoximation becomes more sharp,  when eps --> 0
%
%               f=eps*(abs(x/eps) + 1./(abs(x/eps)+1) - 1);
%
%       In other words,
%                h(u) = abs(u) + 1./(abs(u)+1) - 1;  
%                f(x) = eps*h(x/eps)
%
% Call:
%   [f,df,d2f]=fast_abs_smooth_mz(x,eps)
%
% Input:
%   x      -  vector or matrix argument 
%   eps  -   smoothing parameter
%
% Output:
%   f       - absolute value approximation 
%   df    - derivative of f
%   d2f    - second derivative of f
%

% Michael Zibulevsky 10.07.2008  05.06.2008

tau       = (1/eps)*x;
abs_tau   = abs(tau);
tmp=1./(abs_tau+1);

f       = eps*(abs_tau + tmp - 1);

if nargout>1,
	tmp1=tmp.*tmp;
    df = tau.*(abs_tau + 2).*tmp1;
end
if nargout>2,
    d2f = (2/eps).*tmp1.*tmp;
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a=5;
da=0.001;
f=[-a:da:a]';
n=length(f);

[phi,dphi,d2phi]= fast_abs_smooth_mz(f, 0.1);
figure;plot(f,[phi,dphi,d2phi]);grid
figure;plot(f(1:n-1),[(phi(2:n)-phi(1:n-1))/da dphi(1:n-1)]);grid
figure;plot(f(1:n-1),[(dphi(2:n)-dphi(1:n-1))/da d2phi(1:n-1)]);grid
figure;plot(f(1:n-1),[d2phi(1:n-1)]);grid
figure;plot(f(1:n-1),[dphi(1:n-1)]);grid
figure;plot(f(1:n-1),(phi(2:n)-phi(1:n-1))/da);grid
figure;plot(f(1:n-1),(dphi(2:n)-dphi(1:n-1))/da);grid
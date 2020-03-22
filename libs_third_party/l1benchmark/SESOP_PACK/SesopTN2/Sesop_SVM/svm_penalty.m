function [fsum,g,HZ]=svm_penalty(u,Z,par)
%SVM penalty for violation of separating strip: mu*sum(|x|+x), 
%    where x=u+1;
%          mu=0.5*par.c/length(u)
%          f(x)=|.| is smoothed with parameter par.eps_smooth_const_linear:
%                          h(t) = abs(t) + 1./(abs(t)+1) - 1;
%                          f(x) = eps*h(x/eps)
%
% Call: [f,g,HZ]=svm_penalty(u,Z,par)
%
% Input:
%    u - argument
%    Z - matrix to be multiplied by the Hessian (if needed)
%
% Output: function value, gradient and Hessian*Z
%
% if par.flagXnew==0, then old persistent values of inner variables are reused 

% Michael Zibulevsky, 05.08.2008        
%
% Copyright (c) 2008. All rights reserved. Free for academic use. No warranty 

persistent flagF flagG flagH fsum_old g_old  tau abs_tau tmp tmp1 mud2f 

if par.flagXnew
	flagF=1;flagG=1;flagH=1;
end

n=length(u);
mu=0.5*par.c;
eps=par.eps_smooth_const_linear;

if flagF,
	x=u+1;  % shift of the penalty
	tau       = (1/eps)*x;
	abs_tau   = abs(tau);
	tmp=1./(abs_tau+1);
	f = eps*(abs_tau + tmp - 1);
	fsum=mu*(sum(f+x));
    fsum_old=fsum;
	flagF=0;
else
	fsum=fsum_old;
end


if nargout>1,    % Gradient
	if flagG
		tmp1=tmp.*tmp;
		df = tau.*(abs_tau + 2).*tmp1;
		g= mu*(df + 1);
		g_old=g;
		flagG=0;
	else
		g=g_old;
	end
end

if nargout>2,  	 % Hessian-matrix product
	if  flagH
		mud2f = mu*(2/eps).*tmp1.*tmp;
		flagH=0;
	end
	%HZ=(mud2f*ones(1,size(Z,2))).*Z;
	HZ = muldm(mud2f(:), Z);
end




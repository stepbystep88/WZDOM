function [fsum,g,HZ]=svm_penalty_huber(u,Z,par)
%SVM quadratic-linear penalty for violation of separating strip: mu*sum_f(x), 
%    where x=u+1;
%          mu=0.5*par.c/length(u)
%          h(x)= 0,       x<=0; 
%                x^2,   0<x<1/2; 
%                x-1/4,   x>=1/2
%
%          f(x) = eps*h(x/eps)
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

persistent flagF flagG flagH fsum_old g_old  tau ind2 ind3 abs_tau tmp tmp1 mud2f diag_mud2f_sparse

if par.flagXnew
	flagF=1;flagG=1;flagH=1;
end

n=length(u);
mu=par.c;
eps=par.eps_smooth_const_linear;

if flagF,
	%x=u+1;  % shift of the penalty
	tau       = (1/eps)*(u+1);

	tau_less_0=(tau<=0);
	tau_gt_1 =(tau>=1);
	
	%ind1=find(tau_less_0);
	ind2=find(~(tau_less_0 |tau_gt_1));
	ind3=find(tau_gt_1);
	fsum= (mu*eps)*(0.5*sum(tau(ind2).^2) + sum(tau(ind3) - 0.5));

	fsum_old=fsum;
	flagF=0;
else
	fsum=fsum_old;
end


if nargout>1,    % Gradient
	if flagG
		g=zeros(n,1);
		g(ind2)= mu*tau(ind2);
		g(ind3)=mu;
		
		g_old=g;
		flagG=0;
	else
		g=g_old;
	end
end

if nargout>2,  	 % Hessian-matrix product
	if  flagH
		
% 		mud2f=zeros(n,1);
% 		mud2f=spalloc(n,1,n);
% 		mud2f(ind2) = 2*mu/eps;
% 		mud2f_sparse=spalloc(n,1,n);
% 		mud2f_sparse(ind2) = 2*mu/eps;
%       diag_mud2f_sparse=spdiags(mud2f_sparse,0,n,n);

		diag_mud2f_sparse=sparse(ind2,ind2,mu/eps,n,n);
		
		flagH=0;
	end
	%%HZ=(mud2f*ones(1,size(Z,2))).*Z;
	%HZ = muldm(mud2f(:), Z);
	HZ=diag_mud2f_sparse*Z;
end




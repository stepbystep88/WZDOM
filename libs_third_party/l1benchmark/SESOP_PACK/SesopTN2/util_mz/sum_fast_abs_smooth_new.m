function  [fsum,g,HZ,diagH]=sum_fast_abs_smooth_new(x,Z,par);
% Penalty for  vector x:  mu*sum(|x|), 
%    |.| is smoothed with parameter par.eps_smooth_abs:
%                       h(u) = abs(u) + 1./(abs(u)+1) - 1;
%                       f(x) = eps*h(x/eps)
%
% Call: [fsum,g,HZ,diagH]=sum_fast_abs_smooth_new(x,Z,par)
%
% Input:
%    x - separating vector w=[w;b] in SVM
%    Z - matrix to be multiplied by the Hessian (if needed)
%
% Output: function value, gradient, Hessian*Z, and diagonal of Hessian
%
% if par.flagXnew==0, then old persistent values of inner variables are reused 

% Michael Zibulevsky, 05.08.2008; 07.09.2008      
%
% Copyright (c) 2008. All rights reserved. Free for academic use. No warranty 

persistent flagF flagG flagH fsum_old g_old  tau abs_tau tmp tmp1 mud2f 

if par.flagXnew
	flagF=1;flagG=1;flagH=1;
end

mu=par.weight_abs_penalty;
eps=par.eps_smooth_abs;

if flagF,
	tau       = (1/eps)*x;
	abs_tau   = abs(tau);
	tmp=1./(abs_tau+1);
	f = eps*(abs_tau + tmp - 1);
	fsum=mu*(sum(f));
    fsum_old=fsum;
	flagF=0;
else
	fsum=fsum_old;
end


if nargout>1,    % Gradient
	if flagG
		tmp1=tmp.*tmp;
		df = tau.*(abs_tau + 2).*tmp1;
		g= mu*df;
		g_old=g;
		flagG=0;
	else
		g=g_old;
	end
end

if nargout>2   	   % Hessian-matrix product
   if  flagH
      mud2f = mu*(2/eps).*tmp1.*tmp;
      flagH=0;
   end
   if  numel(Z)>0,
      %HZ=(mud2f*ones(1,size(Z,2))).*Z;
      HZ = muldm(mud2f(:), Z);
   else
      HZ=[];
   end
end

if nargout>3,       % diagonal of Hessian
	diagH=mud2f;
end















% % [f,g,HZ,d2phi]=diag_penalty(x,Z,par);
% % Input:
% % x - vector argument
% % Z - matrix to be multiplied by the Hessian (if needed)
% % par.weight_abs_penalty
% % par.eps_smooth_abs
% 
% mu=1;
% n=length(x(:));
% x=x(1:n-1);
% 
% 
% if nargout == 1   % Objective function
%   [phi]=fast_abs_smooth_mz(x,par.eps_smooth_abs);
%   f=mu*sum(phi(:));
% end
% 
% if nargout==2,   % Gradient
%   [phi,dphi]=fast_abs_smooth_mz(x,par.eps_smooth_abs);
%   f= mu*sum(phi(:));
%   g= [mu*dphi;0];
% end
% 
% if nargout>2,  % Objective function, Gradient,
%                % Hess-vector product; Hess diagonal
%   [phi,dphi,d2phi]=fast_abs_smooth_mz(x,par.eps_smooth_abs);
%   f=mu*sum(phi(:));
%  g= [mu*dphi;0];
%   d2phi=[mu*d2phi(:);0];
%   if   isempty(Z), HZ=[];
%   else HZ = muldm(d2phi,Z);
%  end
% end


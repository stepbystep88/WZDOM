function [f,g,HZ,diagH]=diag_quadr_penalty(x,Z,par)
% f=0.5*par.quadrpenpar*sumsqr(x(:));
%
%Call: [f,g,HZ,diagH]=diag_quadr_penalty(x,Z,par)

% Michael Zibulevsky 04.08.2008

mu=par.quadrpenpar;
f=0.5*mu*sumsqr(x(:));
if nargout >1,
	g=mu*x;
end

if nargout >2,
	HZ=mu*Z;
end

if nargout >3,
	diagH=mu*ones(length(x(:)),1);
end


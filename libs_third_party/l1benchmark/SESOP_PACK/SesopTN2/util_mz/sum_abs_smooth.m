function  [f,g,HW,d2phi]=sum_abs_smooth(x,W,par);
%  Smoothed L1-norm penalty:  mu*sum(abs_smooth(x));  
%   calls fast_abs_smooth_mz or abs_smoothed_eps
%
%Call:
%    [f,g,HW,d2phi]=sum_abs_smooth(x,W,par);
%
% Input:
%   x - vector argument
%  W - matrix to be multiplied by the Hessian (if needed)
%   par.weight_abs_penalty  -  defines value of mu
%   par.eps_smooth_abs      -  smoothing parameter of absolute value approximation
%   par.flagXnew:  1 - x has been changed from the previous call of the function
%                          0 - x the same as in previous call; persistent precomputed variables are valid for use
%
% Output:
%   f - function value
%   g - gradient
%   HW - Hessian-matrix product
%   d2phi - diagonal of the Hessian
%
%Calls other functions:
%    abs_smoothed_eps 

% Michael Zibulevsky 10.07.2008

global GlobalTimeFunc_x;

cputime0=cputime;

persistent f_old g_old  d2phi_old   flagFgh

if par.flagXnew,  
	flagFgh=[ 1 1 1];   % flagFgh indicates need to compute new values f,g,d2phi (when 1) 
	%                               or use stored persistent values from previous iterations (when 0)
end


mu=par.weight_abs_penalty;




if nargout == 1   % Objective function
	if flagFgh(1)
		%[phi]=par.abs_approx(x,par.eps_smooth_abs);
		%[phi]=fast_abs_smooth_mz(x,par.eps_smooth_abs);
		[phi]=abs_smoothed_eps(x,par.eps_smooth_abs);
		f=mu*sum(phi(:));
		f_old=f;
		flagFgh(1)=0;
	else
		f=f_old;
	end
end

if nargout==2,   % Gradient
	if any(flagFgh(1:2))
		%[phi,dphi]=par.abs_approx(x,par.eps_smooth_abs);
		%[phi,dphi]=fast_abs_smooth_mz(x,par.eps_smooth_abs);
		[phi,dphi]=abs_smoothed_eps(x,par.eps_smooth_abs);
		f= mu*sum(phi(:));
		g= mu*dphi;
		f_old=f;g_old=g;
		flagFgh(1:2)=0;
	else
		f=f_old;
		g=g_old;
	end
end

if nargout>2,  % Objective function, Gradient, Hess-vector product; Hess diagonal
	if any(flagFgh)
		%[phi,dphi,d2phi]=par.abs_approx(x,par.eps_smooth_abs);
		%[phi,dphi, d2phi]=fast_abs_smooth_mz(x,par.eps_smooth_abs);
		[phi,dphi, d2phi]=abs_smoothed_eps(x,par.eps_smooth_abs);
		f=mu*sum(phi(:));
		g=mu*dphi;
		d2phi=mu*d2phi(:);
		 f_old=f;g_old=g;d2phi_old=d2phi;
		 flagFgh(1:3)=0;
	else
		 f=f_old;g=g_old; d2phi=d2phi_old;
	end
	if   isempty(W), HW=[];
	else HW = muldm(d2phi,W);
	end
end

GlobalTimeFunc_x=GlobalTimeFunc_x+cputime-cputime0;


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function B=muldm(d,A)  % Compute efficiently   diag(d)*A
%B=muldm(d,A)
%
% B= diag(d)*A
d=d(:);
[M,N]=size(A);
n=length(d);
if n ~= M, error('Wrong matrix sizes');end

B=A;
for j=1:N
B(:,j)=d.*B(:,j);
end

end



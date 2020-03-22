function [f,g,HZ]= ls_fgh(u,Z,par)
% Least-squares objective function:  f=1/2*||u-par.y||^2
%
%Call:
%    [f,g,HZ]=ls_fgh(u,Z,par):   
%
% Input:
%     u - least-square argument:   f= 0.5*||u-par.y||^2
%     Z - matrix to be multiplied by the Hessian (if needed)
%     par.flagXnew:  1 - x has been changed from the previous call of the function
%                    0 - x the same as in previous call; persistent precomputed variables are valid for use
%
% Output:
%   f - function value
%   g - gradient
%   HZ - Hessian-matrix product


persistent f_old  tmp;

if par.flagXnew
	tmp=u - par.y;
	f= 0.5*  (tmp(:)' *tmp(:));
	f_old=f;
else
	f=f_old;   %disp('flagFgh=0')
end

if nargout>1 
  g=tmp;
end

if nargout>2, % multiplication Z by the Hessian
  HZ=Z;
end


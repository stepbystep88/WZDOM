function x = CGmz(func,x0,b,par,varargin)
% Solve equation Hx = b with Conjugate Gradients
%  Hx is computed by  user function func(x,par,varargin{:});

global Niter

% begin initialization
epsilon = 1e-9;
epsilonsq=epsilon^2;

beta=0;
p=0;
x=x0;
Hx=feval(func,x,par,varargin{:});
r= Hx - b;
sqnormr=norm(r(:))^2;
sqnormr0=sqnormr;

for Niter=1:par.max_iter,
    p=-r+beta*p;
    Hp=feval(func,p,par,varargin{:});
    alpha = sqnormr/(p(:)' * Hp(:));
    x=x+alpha*p;
    if mod(Niter,5)==0,
                % Call user-defined function for progress report 
                if isfield(par, 'report'), par.report(x,Niter,par,varargin{:}); end  
    end
    r=r+alpha*Hp;
    sqnormr_new=norm(r(:))^2;
    beta=sqnormr_new/sqnormr;
    sqnormr=sqnormr_new;
    if sqnormr/sqnormr0<epsilonsq, break;end
end

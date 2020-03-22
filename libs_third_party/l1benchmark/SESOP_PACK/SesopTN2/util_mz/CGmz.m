function [x,report] = CGmz(multH,x0,b,options,par,varargin)
% Solve equation Hx = b with Conjugate Gradients, 
%  where H - symmetric, positive semidefinite matrix
%  Hx is computed by  user function multH(x,par,varargin{:});

global Niter

epsilon = 1e-10; % Required accuracy in norm of residual r=Hx-b
epsilonsq=epsilon^2;

if isfield(options,'cg_figure_handle') && ~isempty(options.cg_figure_handle)
	cg_figure_handle=options.cg_figure_handle;
elseif options.period_show_progress>0
	cg_figure_handle=figure('name', 'Linear CG grad. norms');
end


beta=0;
p=0;
x=x0;
Hx= multH(x,par,varargin{:});
r= Hx - b;
normr=norm(r(:));
sqnormr=normr^2;
sqnormr0=sqnormr;

report.gradnorms=[];
for Niter=1:options.max_iter_cg,

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %          Progress report
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  report.Niter=Niter;
  report.gradnorms=[report.gradnorms normr];
  if mod(Niter,options.period_show_progress)==0,
	figure(cg_figure_handle);
	semilogy([0:length(report.gradnorms)-1],report.gradnorms);grid
    if isfield(options, 'report_func'),
      options.report_func(x,report,par,varargin{:});
    end
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %   CG iteration
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  p=-r+beta*p;
  Hp=multH(p,par,varargin{:});
  alpha = sqnormr/(p(:)' * Hp(:));        
  %alpha = -(r(:)'*p(:))/(p(:)' * Hp(:)); % This line is from my lecture, but Gill-Murray-Wright (previous line) seems to be a bit better
  x=x+alpha*p;
  r=r+alpha*Hp;
  normr=norm(r(:));
  sqnormr_new=normr^2;
  beta=sqnormr_new/sqnormr;
  sqnormr=sqnormr_new;
  if sqnormr/sqnormr0<epsilonsq, break;end

end

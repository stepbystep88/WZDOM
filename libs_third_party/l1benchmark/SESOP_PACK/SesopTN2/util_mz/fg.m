function [f,g]=fg(x,par)
% Value and gradient of f=par.f_u(Ax)+ par.f_x(x); For use by fminunc 

%global GlobalNiter GlobalGradNorms GlobalFuncValues GlobalTimes

u=par.multA(x,par);

if nargout==1,
	f_u=par.func_u(u,[],par);   % Function value
	f2=par.func_x(x,[],par);
	f=f_u+f2;
else
	[f_u,g_u]=par.func_u(u,[],par);
	[f2,g2_x]=par.func_x(x,[],par);
	f=f_u+f2;
	
	g1_x=par.multAt(g_u,par);        % Gradient
	g = g1_x+g2_x;
   

%    GlobalGradNorms=[GlobalGradNorms, norm(g)];
%    GlobalFuncValues=[GlobalFuncValues,f];
%    GlobalTimes=[GlobalTimes, cputime];
% 
%    report.gradnorms=GlobalGradNorms;    % array of gradient norms with iteration
%    report.func_values=GlobalFuncValues; % array of objective function values with iteration;
%    
%    if mod(GlobalNiter,par.period_show_progress)==0, % || GlobalNiter==1,
%       report.Niter = GlobalNiter;
%       par.report_func(x,report,par);
%    end
%    GlobalNiter=GlobalNiter+1;
end

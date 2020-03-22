function [f,g]=report_global(x,f,g,par)
% For use by fminunc 

global GlobalNiter GlobalGradNorms GlobalFuncValues GlobalTimes


   GlobalGradNorms=[GlobalGradNorms, norm(g)];
   GlobalFuncValues=[GlobalFuncValues,f];
   GlobalTimes=[GlobalTimes, cputime];

   report.gradnorms=GlobalGradNorms;    % array of gradient norms with iteration
   report.func_values=GlobalFuncValues; % array of objective function values with iteration;
   
   if mod(GlobalNiter,par.period_show_progress)==0, % || GlobalNiter==1,
      report.Niter = GlobalNiter;
      par.report_func(x,report,par);
   end
   GlobalNiter=GlobalNiter+1;
end

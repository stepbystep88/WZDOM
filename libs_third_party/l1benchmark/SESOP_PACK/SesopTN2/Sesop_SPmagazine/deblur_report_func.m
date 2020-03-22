function deblur_report_func(c,report,par)
%plot current optimization results 

% Input:  c                   - current synthesis coefficients
%         report.gradnorms    - array of gradient norms with iteration
%         report.func_values  - array of objective function values with iteration;
%         report.Niter        - current sesop_iter;   ;

global GlobalXsnr GlobalY_RNR GlobalNNiter GlobalSNRtime
GlobalNNiter=[GlobalNNiter report.Niter];

x=par.SynthesisFunc(c,par);
X=reshape(x,par.imagesize);

%Xsnr=10*log10(sumsqr(par.x00)/sumsqr(x-par.x00));
Xsnr= -10*log10(sumsqr( par.x00/norm(par.x00) - x/norm(x)));

y= par.projector(x,par);
Y=reshape(y,par.proj_size);
Y_RNR=10*log10(sumsqr(par.y-y)/(numel(y)*par.sigma_noise^2));  % Residual-to-noise ratio


GlobalXsnr=[GlobalXsnr Xsnr];
GlobalY_RNR=[GlobalY_RNR Y_RNR]; % Residual-to-noise ratio
GlobalSNRtime=[GlobalSNRtime cputime];

figure(par.BPblur_report_figure_handle);
subplot(221);semilogy(report.gradnorms);title('Norm of Gradient');grid
subplot(222);fbest=min(report.func_values);semilogy(report.func_values-fbest);grid
title(sprintf('Function value minus best (%g)',fbest));
subplot(223);plot(GlobalNNiter,GlobalXsnr,GlobalNNiter,GlobalY_RNR);%legend('Xsnr','Yrnr');
title('SNR (blue) of recovered image and YRNR');grid
subplot(224);imagesc(X);title('Recovered Image');colorbar
drawnow
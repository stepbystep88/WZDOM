function c=analysis_blur_adj(y,par)
%Adjoint to blur_synthesis
kern=flipud(fliplr(par.blur_kernel));
x=convn(y,kern,'same');
c=par.AnalysisFunc(x,par);


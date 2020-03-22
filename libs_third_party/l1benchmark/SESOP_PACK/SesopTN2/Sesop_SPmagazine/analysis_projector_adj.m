function c=analysis_projector_adj(y,par)
%Adjoint to projector_synthesis

global GlobalTimeMultAt

cputime0=cputime;

x=par.projector_adj(y,par);
c=par.AnalysisFunc(x,par);


GlobalTimeMultAt=GlobalTimeMultAt+cputime-cputime0;


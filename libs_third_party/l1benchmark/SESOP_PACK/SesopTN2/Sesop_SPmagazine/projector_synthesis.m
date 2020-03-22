function y=projector_synthesis(c,par)
%Apply synthesis operator (e.g. wavelets) and then projector (e.g. Blur, Radon)
global GlobalTimeMultA

cputime0=cputime;

x=par.SynthesisFunc(c,par);
y=par.projector(x,par);

GlobalTimeMultA=GlobalTimeMultA+cputime-cputime0;


function y=blur_synthesis(c,par)
%Combine blur and [wavelet] synthesis

x=par.SynthesisFunc(c,par);
y=convn(x,par.blur_kernel,'same');


function z=TI2Dsynthesismz(x,par)

cqf    =par.wavelet_cqf;
levels =par.wavelet_levels;  

n=par.imagesize(1);

coefs=reshape(x,n,length(x(:))/n);
z = mirdwt_TI2D(coefs, cqf, levels);
z=z(:);


function [coefs,coefs_size]=TI2Danalysismz(x,par)

cqf    =par.wavelet_cqf;
levels =par.wavelet_levels; 



X=reshape(x,par.imagesize);
coefs = mrdwt_TI2D(X, cqf, levels);
coefs_size =size(coefs);
coefs =coefs(:);


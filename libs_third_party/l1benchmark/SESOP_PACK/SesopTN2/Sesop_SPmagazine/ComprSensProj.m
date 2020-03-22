function y=ComprSensProj(x,par)
% Compute sub-sampled 2d FFT (indeces in par.ComprSens.ind)

X=reshape(x,par.imagesize);
[m,n]=size(X);
FX=(1/sqrt(m*n))*fftn(X);
y1=FX(par.ComprSens.ind);
y=[real(y1(:));imag(y1(:))];
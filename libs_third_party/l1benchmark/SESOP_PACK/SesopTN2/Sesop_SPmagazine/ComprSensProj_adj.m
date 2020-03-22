function x=ComprSensProj_adj(y,par)

k=numel(y)/2;
Y=zeros(par.imagesize);
Y(par.ComprSens.ind)=y(1:k)+1i*y(k+1:end);

[m,n]=size(Y);
X=sqrt(m*n)*ifftn(Y);

x=real(X(:));
%x=[real(X(:));imag(X(:))];
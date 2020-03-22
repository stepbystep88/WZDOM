function x=blurring_adj(y,par)
% Adjoint convolution
kern=flipud(fliplr(par.blur_kernel));
x=convn(reshape(y,par.imagesize) ,kern,'same');
x=x(:);


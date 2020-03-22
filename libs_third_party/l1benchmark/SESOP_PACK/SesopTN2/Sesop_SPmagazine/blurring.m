function y=blurring(x,par)
% Convolve signal/image with blur kernel 

y=convn(reshape(x,par.imagesize),par.blur_kernel,'same');
y=y(:);

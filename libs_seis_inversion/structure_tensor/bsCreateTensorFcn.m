function [T, rotate] = bsCreateTensorFcn(f)
% https://github.com/gpeyre/numerical-tours/blob/bf78b5774b9d5a1e6c8089746c89661004896c3d/matlab/pde_3_diffusion_tensor.ipynb

    [n1, n2] = size(f);
    
    
    t = @(n)interp1([1 round(n/2) n], [0 round(n/2) 0], 1:n, 'linear');
    [X2,X1] = meshgrid(t(n2), t(n1));
%     t = [0:n1/2 -n1/2+1:-1];
%     [X2,X1] = meshgrid(t,t);

    normalize = @(h)h/sum(h(:));
    
    cconv = @(f,h)real(ifft2(fft2(f).*repmat(fft2(h),[1 1 size(f,3)])));
    
    h = @(sigma)normalize( exp( -(X1.^2+X2.^2)/(2*sigma^2) ) );
    
%     blur = real(ifft2(fft2(f).*repmat(fft2(h),[1 1 size(f,3)])));
    
    blur = @(f,sigma)cconv(f,h(sigma));

    options.order = 2;
    nabla = @(f)grad(f,options);
    
    tensorize = @(u)cat(3, u(:,:,1).^2, u(:,:,2).^2, u(:,:,1).*u(:,:,2));
    
    rotate = @(T)cat(3, T(:,:,2), T(:,:,1), -T(:,:,3));
    
    T = @(f,sigma)blur( tensorize( nabla(f) ), sigma);
end
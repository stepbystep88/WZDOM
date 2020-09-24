function blur = bsCreateBlurFcn(f)
    [n1, n2, n3] = size(f);
    
    t = @(n)interp1([1 round(n/2) n], [0 round(n/2) 0], 1:n, 'linear');
    
    normalize = @(h)h/sum(h(:));
    
    if n3 ~= 1
        cconv = @(f,h)real(ifftn(fftn(f).* fftn(h)));
        [X2, X1, X3] = meshgrid(t(n2), t(n1), t(n3));
        h = @(sigma)normalize( exp( -(X1.^2+X2.^2+X3.^2)/(2*sigma^2) ) );
    else
        cconv = @(f,h)real(ifft2(fft2(f).*repmat(fft2(h),[1 1 size(f,3)])));
        [X2,X1] = meshgrid(t(n2), t(n1));
        h = @(sigma)normalize( exp( -(X1.^2+X2.^2)/(2*sigma^2) ) );
	end
    
    blur = @(f,sigma)cconv(f,h(sigma));
    
end
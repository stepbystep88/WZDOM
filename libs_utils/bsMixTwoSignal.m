function z = bsMixTwoSignal(x, y, fs1, fs2, dt)
    
    n = length(x);
    fs = 1/dt/n*(0:n/2);
    
    fftx = fft(x, n);
%     fftx = fftshift(fftx);
    
    ffty = fft(y, n) * 0.7;
%     ffty = fftshift(ffty);
    
    fz = zeros(size(fftx));
    fz(1) = fftx(1);
    
    [~, i1] = min(abs(fs - fs1));
    [~, i2] = min(abs(fs - fs2));
    coef = linspace(0, 1, i2-i1+1)';
    fz(2:i1-1) = fftx(2:i1-1);
    fz(i1:i2) = fftx(i1:i2).*(1-coef) + ffty(i1:i2).*coef;
    fz(i2+1:n/2+1) = ffty(i2+1:n/2+1);
    fz(n/2+2:end) = conj(flipud(fz(2:n/2)));
    z = ifft(fz, n);
end
function freq = bsGetMainFreq(postSeisData, dt)
    N = 512; 
    j = 0 : N - 1; 

    wellNum = size(postSeisData, 2);
    bestft = zeros(wellNum, 1);
    for i = 1 : wellNum
        Fs = 1000.0 / dt; 
        fr = fft(postSeisData(:, i), N);
        mag = sqrt(real(fr).^2 + imag(fr).^2);
        f = j * Fs / N;  
        bestft(i, 1) = min( f(  mag(:, 1) == max(mag)  ) );
    end

    freq = mean( bestft );
    fprintf('The main frequency is :%.2f\n', freq);
        
end
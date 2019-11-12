function [ W ] = bsWaveletMatrix(sampNum, wavelet, waveletFreq, dt)

    if isempty(wavelet) && nargin == 4
        % generate ricker wavelet
        wavelet = bsGenRickerWavelet(waveletFreq, dt);
    end
    
    [~, ix] = max(abs(wavelet));

    Wmatrix = convmtx(wavelet, sampNum);
    W = Wmatrix(ix : ix+sampNum-1, :);
    
end

function [wavelet] = bsGenRickerWavelet(freq, dt)
    wave = s_create_wavelet({'type','ricker'}, {'frequencies', freq}, {'step', dt}, {'wlength', 120});   
    wavelet = wave.traces;                                      
end


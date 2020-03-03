function wavelet = bsGenWavelet(type, freq, dt, wavelet)
    switch type
        case 'ricker'
            wave = s_create_wavelet({'type','ricker'}, {'frequencies', freq}, {'step', dt}, {'wlength', 80});
        case 'zero-phase'
            wave = s_create_wavelet({'type','zero-phase'}, ...
                {'frequencies', freq-20,freq-10,freq+10,freq+20}, ...
                {'step', dt}, ...
                {'wlength', 100});
        case 'min-phase'
            wave = s_create_wavelet({'type','min-phase'}, ...
                {'frequencies', freq-10,freq-5,freq+5,freq+10}, ...
                {'step', dt}, ...
                {'wlength', 100});
        case 'input'
            if isempty(wavelet)
                disp(wavelet);
                error('When type is input, the input struct must contain wavelet field.');
            else
                wave.traces = wavelet;
            end
        otherwise
            validatestring(type, {'rikcer', 'zero-phase', 'min-phase', 'input'});
    end
    wavelet = wave.traces;          
end
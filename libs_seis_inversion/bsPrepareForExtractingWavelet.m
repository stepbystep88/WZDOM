function [GInvParam, postSeisData, horizon] = bsPrepareForExtractingWavelet(GInvParam, usedTimeLine, wellLogs, type)

    sampNum = GInvParam.upNum + GInvParam.downNum; 
    dt = GInvParam.dt;
    
    wells = cell2mat(wellLogs);
    inIds = [wells.inline];
    crossIds = [wells.crossline];
    
    % horizon of the traces at well location
    horizon = bsCalcHorizonTime(usedTimeLine, inIds, crossIds);
    startTime = horizon - GInvParam.upNum * dt;
    % read seismic data
    postSeisData = bsGetPostSeisData(GInvParam, inIds, crossIds, startTime, sampNum-1);
    
    freq = bsGetMainFreq(postSeisData, dt);
    GInvParam.waveletFreq = freq;
    
    % create wavelet
    switch type
        case 'ricker'
            wave = s_create_wavelet({'type','ricker'}, {'frequencies', freq}, {'step', GInvParam.dt}, {'wlength', 80});
        case 'zero-phase'
            wave = s_create_wavelet({'type','zero-phase'}, ...
                {'frequencies', freq-10,freq-5,freq+5,freq+10}, ...
                {'step', GInvParam.dt}, ...
                {'wlength', 120});
        case 'min-phase'
            wave = s_create_wavelet({'type','min-phase'}, ...
                {'frequencies', freq-10,freq-5,freq+5,freq+10}, ...
                {'step', GInvParam.dt}, ...
                {'wlength', 120});
        case 'input'
            if isempty(GInvParam.wavelet)
                disp(GInvParam);
                error('When type is input, the input struct must contain wavelet field.');
            else
                wave.traces = GInvParam.wavelet;
            end
        otherwise
            validatestring(type, {'rikcer', 'zero-phase', 'min-phase', 'input'});
    end
    wavelet = wave.traces;          
    GInvParam.wavelet = wavelet;
end

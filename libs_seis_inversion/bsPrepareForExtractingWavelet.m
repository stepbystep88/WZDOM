function [GInvParam, postSeisData, horizon] = bsPrepareForExtractingWavelet(GInvParam, usedTimeLine, wellLogs, type)

    sampNum = GInvParam.upNum + GInvParam.downNum; 
    dt = GInvParam.dt;
    
    wells = cell2mat(wellLogs);
    inIds = [wells.inline];
    crossIds = [wells.crossline];
    
    % horizon of the traces at well location
    horizon = bsGetHorizonTime(usedTimeLine, inIds, crossIds);
    startTime = horizon - GInvParam.upNum * dt;
    % read seismic data
    postSeisData = bsGetPostSeisData(GInvParam, inIds, crossIds, startTime, sampNum-1);
    
    freq = bsGetMainFreq(postSeisData, dt);
    GInvParam.waveletFreq = freq;
    
    % create wavelet
    GInvParam.wavelet = bsGenWavelet(type, freq, dt, GInvParam.wavelet);

end

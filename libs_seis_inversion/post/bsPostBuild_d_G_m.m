function [d, G, m] = bsPostBuild_d_G_m(GPostInvParam, inline, crossline, startTime, initLog, model)

    sampNum = GPostInvParam.upNum + GPostInvParam.downNum; 
    
    [postSeisData, GPostInvParam.postSeisData.segyInfo] = bsReadTracesByIds(...
        GPostInvParam.postSeisData.fileName, ...
        GPostInvParam.postSeisData.segyInfo, ...
        inline, ...
        crossline, ...
        startTime, ...
        sampNum,...
        GPostInvParam.dt);
    
    d = postSeisData(1 : end-1);
    if isfield(GPostInvParam, 'seismicFiltCoef') && ~isempty(GPostInvParam.seismicFiltCoef)
        d = bsButtLowPassFilter(d, GPostInvParam.seismicFiltCoef);
    end
    
    if isempty(model)
        G = bsPostGenGMatrix(GPostInvParam.wavelet, sampNum);
    else
        G = model.orginal_G;
    end
    
    m = log(initLog);
end
function [d, G, m] = bsPostBuild_d_G_m(GPostInvParam, inline, crossline, startTime, initLog, model)

    sampNum = GPostInvParam.upNum + GPostInvParam.downNum; 
    
    switch lower(GPostInvParam.postSeisData.mode)
        case 'segy'
            [postSeisData, GPostInvParam.postSeisData.segyInfo] = bsReadTracesByIds(...
            GPostInvParam.postSeisData.fileName, ...
            GPostInvParam.postSeisData.segyInfo, ...
            inline, ...
            crossline, ...
            startTime, ...
            sampNum,...
            GPostInvParam.dt);
        case 'function' % get initial model by calling a function
            
            if isempty(GPostInvParam.postSeisData.fcn)
                error('When GInvParam.postSeisData.mode is function, the GInvParam.postSeisData.fcn could not be empty!\n');
            end
            postSeisData = GPostInvParam.postSeisData.fcn(inline, crossline, startTime);
    end
    
    d = postSeisData(1 : sampNum-1);
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
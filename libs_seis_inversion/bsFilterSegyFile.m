function bsFilterSegyFile(GInvParam, timeLine, inputFileName, outFileName, GSegyInfo, varargin)

    p = inputParser;
    
    [rangeInline, rangeCrossline] = bsGetWorkAreaRange(GSegyInfo, inputFileName);
    
    addParameter(p, 'rangeInline', rangeInline);
    addParameter(p, 'rangeCrossline', rangeCrossline);
    addParameter(p, 'expandNum', 30);
    addParameter(p, 'filtCoef', 0.1);
    
    p.parse(varargin{:});  
    options = p.Results;
    
    GInvParam.upNum = GInvParam.upNum + options.expandNum;
    GInvParam.downNum = GInvParam.downNum + options.expandNum;
    
    rangeInline = options.rangeInline;
    rangeCrossline = options.rangeCrossline;
    
    [inIds, crossIds] = bsGetCDPsByRange(rangeInline, rangeCrossline);
    horizonTime = bsGetHorizonTime(timeLine{GInvParam.usedTimeLineId}, inIds, crossIds, 1);
    startTime = horizonTime - GInvParam.upNum*GInvParam.dt;
    sampNum = GInvParam.upNum + GInvParam.downNum;

    data = bsReadTracesByIds(inputFileName, GSegyInfo, inIds, crossIds, startTime, sampNum, GInvParam.dt);
    
    filteredData = bsFilterProfileData(data, options.filtCoef);
    
    res.inIds = inIds;
    res.crossIds = crossIds;
    res.horizon = horizonTime;
    res.upNum = GInvParam.upNum;
    res.dt = GInvParam.dt;
    bsWriteInvResultIntoSegyFile(res, filteredData, inputFileName, GSegyInfo, outFileName);
end
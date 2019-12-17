function wellLogs = bsUpdateWelllogsByHorizon(GInvParam, timeLine, wellLogs)
    for i = 1 : length(wellLogs)
        
        wellInfo = wellLogs{i};
        
        [~, ~, horizon] = bsGetWellBaseInfo(timeLine{GInvParam.usedTimeLineId}, ...
        wellInfo.inline, wellInfo.crossline, 1, 2, 1, 2, 3);
    
        [~, wellLogs{i}.wellLog] ...
            = bsExtractWellDataByHorizon(wellInfo.wellLog, horizon, 1, ...
            GInvParam.indexInWellData.time, GInvParam.upNum, GInvParam.downNum, 1, GInvParam.dt);
    end
end
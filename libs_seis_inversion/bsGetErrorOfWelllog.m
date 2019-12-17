function [wellLogs, dataIndex, type] = bsGetErrorOfWelllog(GInvParam, timeLine, wellLogs, varargin)
    nWell = length(wellLogs);
    
    p = inputParser;
    addParameter(p, 'expandNum', 0);
    
    p.parse(varargin{:});  
    options = p.Results;
    
    wells = cell2mat(wellLogs);
    wellInIds = [wells.inline];
    wellCrossIds = [wells.crossline];
    
    usedTimeLine = timeLine{GInvParam.usedTimeLineId};
    wellHorizonTimes = bsGetHorizonTime(usedTimeLine, ...
        wellInIds, wellCrossIds);
    % re-organize the data, and add the error data at the last column;
    for i = 1 : nWell
        
        wellInfo = wellLogs{i};
        
        wellInfo = bsGetSynTrace(GInvParam, ...
            wellInfo.inline, ...
            wellInfo.crossline, ...
            wellHorizonTimes(i), ...
            wellInfo, ...
            options.expandNum);
        
        wellLogs{i} = wellInfo;
    end
    
    dataIndex = size(wellInfo.wellLog, 2);
    type = {'error'};
end

function wellInfo = bsGetSynTrace(GInvParam, inline, crossline, horizonTime, wellInfo, expandNum)
%     GInvParam.isNormal = false;
    
    indexInWellData = GInvParam.indexInWellData;
    GInvParam.upNum = GInvParam.upNum + expandNum;
    GInvParam.downNum = GInvParam.downNum + expandNum;
    
    wellData = bsExtractWellDataByHorizon(...
            wellInfo.wellLog, ...
            horizonTime, ...
            1:size(wellInfo.wellLog, 2), ...
            GInvParam.indexInWellData.time, ...
            GInvParam.upNum, ...
            GInvParam.downNum, ...
            1);
        
    startTime = horizonTime - GInvParam.upNum * GInvParam.dt;
    
    switch lower(GInvParam.flag)
    case {'prestack', 'pre-stack'}
        trueLog = wellData(:, ...
            [   indexInWellData.depth, ...
                        indexInWellData.vp, ...
                        indexInWellData.vs, ...
                        indexInWellData.rho]);

        [d, G, m, ~] ...
            = bsPreBuild_d_G_m(GInvParam, inline, crossline, startTime, trueLog);
    
        synData = G * m;
        realData = d;
        
    case {'poststack', 'post-stack'}
        trueLog = wellData(:, ...
            [indexInWellData.ip]);
        
        [d, G, m] = bsPostBuild_d_G_m(GInvParam, inline, crossline, startTime, trueLog, []);
        synData = G * m;
        realData = d;
    end
    
    errorData = realData - synData;
    
    wellInfo.wellLog = [wellData(1:end-1, :), realData, synData, errorData];
end
function [wellLogs, dataIndex, type] = bsGetErrorOfWelllog(GInvParam, timeLine, wellLogs, varargin)
    nWell = length(wellLogs);
    
    p = inputParser;
    addParameter(p, 'expandNum', 0);
    addParameter(p, 'filtCoef', 1);
    
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
            options.expandNum, ...
            options.filtCoef);
        
        wellLogs{i} = wellInfo;
    end
    
    dataIndex = size(wellInfo.wellLog, 2);
    type = {'error'};
end

function wellInfo = bsGetSynTrace(GInvParam, inline, crossline, horizonTime, wellInfo, expandNum, filtCoef)
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
            filtCoef);
        
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
        
%         maxAbsD = norm(d);
%         realData = realData / maxAbsD;
%         synData = synData / maxAbsD;
    end
    
    errorData = realData - synData;
    
%     tmp = [realData, synData, errorData];
%     figure; plot(realData, 'k'); hold on; plot(errorData, 'r'); plot(synData, 'b');
    
    wellInfo.wellLog = [wellData(1:end-1, :), realData, synData, errorData];
end
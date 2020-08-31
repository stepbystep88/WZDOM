function [basicInfo] = bsInitBasicInfoForShowingProfile(GShowProfileParam, GInvParam, wellLogs, timeLine, profile)

    % get horizon information
    horizons = bsGetHorizons(GInvParam, timeLine, profile);
    
    [wellPos, wellIndex, wellNames] = bsFindWellLocation(wellLogs, profile.inIds, profile.crossIds);
    
    if length(unique(profile.inIds)) < length(unique(profile.crossIds))
        traceIds = profile.crossIds;
    else
        traceIds = profile.inIds;
    end
    
    % smooth the horizon
    horizon = bsButtLowPassFilter(profile.horizon, 0.1);
    
    if isempty(wellPos)
        left = 1;
        right = length(profile.inIds);
    else
        left = min(wellPos) - GShowProfileParam.showLeftTrNumByWells;
        right = max(wellPos) + GShowProfileParam.showRightTrNumByWells;

        if left < 1
            left = 1;
        end

        if right > length(profile.inIds)
            right = length(profile.inIds);
        end
    end
        
    basicInfo.horizon = horizon(left:right);
    basicInfo.horizons = horizons(:, left:right);
    basicInfo.inIds = profile.inIds(:, left:right);
    basicInfo.crossIds = profile.crossIds(:, left:right);
    basicInfo.traceIds = traceIds(left:right);
    basicInfo.wellPos = wellPos - left + 1;
    basicInfo.left = left;
    basicInfo.right = right;
    basicInfo.wellIndex = wellIndex;
    basicInfo.wellNames = wellNames;
    basicInfo.upNum = GInvParam.upNum;
    basicInfo.downNum = GInvParam.downNum;
    
    % scale the horiozon and get the time information
    scaleFactor = GShowProfileParam.scaleFactor;
    traceNum = length(basicInfo.traceIds);
    sampNum = GInvParam.upNum + GInvParam.downNum;
    
    dt = GInvParam.dt;
    newDt = dt / scaleFactor;
    
    time0 = basicInfo.horizon - GInvParam.upNum * dt;
    minTime = min(time0) - GShowProfileParam.edgeOffsetNum * dt;
    maxTime = max(time0) + sampNum * dt + GShowProfileParam.edgeOffsetNum * dt;
    
    sequence = linspace(0, dt*(sampNum-1), sampNum);
    timeGrid = repmat(sequence', 1, traceNum)...
        + repmat(time0, sampNum, 1);
    
    nHorizon = size(horizons, 1);
    if GShowProfileParam.isScaleHorizon
        newTraceIds = linspace(basicInfo.traceIds(1), basicInfo.traceIds(end), traceNum * scaleFactor);
        basicInfo.newHorizon = interp1(basicInfo.traceIds, basicInfo.horizon, newTraceIds, 'spline');
        basicInfo.newTraceIds = newTraceIds;
        
        newHorizons = zeros(nHorizon, length(newTraceIds));
        for i = 1 : nHorizon
            newHorizons(i, :) = interp1(basicInfo.traceIds, basicInfo.horizons(i, :), newTraceIds, 'spline');
        end

        basicInfo.newHorizons = newHorizons;
    else
        basicInfo.newHorizon = basicInfo.horizon;
        basicInfo.newTraceIds = basicInfo.traceIds;
        basicInfo.newHorizons = horizons;
    end
    
    
    basicInfo.dt = dt;
    basicInfo.minTime = minTime;
    basicInfo.newTimeSeq = minTime : newDt : maxTime;
    basicInfo.newTime0 = basicInfo.newHorizon - GInvParam.upNum * dt;
    basicInfo.timeSeq = minTime : dt : maxTime;
    
    basicInfo.newDt = newDt;
    basicInfo.scaleFactor = scaleFactor;
    basicInfo.timeGrid = timeGrid;
    
end

function horizons = bsGetHorizons(GInvParam, timeLine, profile)
    inIds = profile.inIds;
    crossIds = profile.crossIds;
    nHorizon = length(timeLine);
    
    horizons = zeros(nHorizon, length(inIds));

    for i = 1 : nHorizon
        horizons(i, :) = bsGetHorizonTime(timeLine{i}, inIds, crossIds, GInvParam.isParallel, GInvParam.numWorkers);
    end
end

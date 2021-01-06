function outWells =bsPreProcessWelllogs(GInvParam, timeLine, wellLogs, wellInfo)

    wellb = wellInfo.wellLog;
    outWells = wellLogs;
    
    switch GInvParam.flag
        case {'prestack', 'pre-stack'}
            dataIndex = [GInvParam.indexInWellData.vp, GInvParam.indexInWellData.vs, GInvParam.indexInWellData.rho];
        case {'poststack', 'post-stack'}
            dataIndex = GInvParam.indexInWellData.ip;
    end
    
    boffset = round((bsGetHorizonTime(timeLine{GInvParam.usedTimeLineId}, wellInfo.inline, wellInfo.crossline) - wellb(1, GInvParam.indexInWellData.time))/GInvParam.dt);
    
    for i =  1 : length(wellLogs)
        horizon = bsGetHorizonTime(timeLine{GInvParam.usedTimeLineId}, wellLogs{i}.inline, wellLogs{i}.crossline);
        outWells{i}.wellLog = bsJoint(wellLogs{i}.wellLog, wellb, horizon, boffset, GInvParam.indexInWellData.time, dataIndex, GInvParam.dt);
    end
            
end

function wella = bsJoint(wella, wellb, horizon, boffset, timeIndex, dataIndex, dt)
%     [~, index] = min(abs(wellb(:, timeIndex) - wella(1, timeIndex)));
    aoffset = round((horizon - wella(1, timeIndex)) / dt);
    index = boffset - aoffset + 1;
    
    x = wella(:, dataIndex);
    na = size(wella, 1);
    y = wellb(index:index+na-1, dataIndex);
    fs = 5/(500/dt);
    
    for i = 1 : length(dataIndex)
        if mod(na, 2) == 0
            wella(:, dataIndex(i)) = bsMixTwoSignal(y(:, i), x(:, i), fs, fs, dt);
        else
            wella(1:end-1, dataIndex(i)) = bsMixTwoSignal(y(1:end-1, i), x(1:end-1, i), fs, fs, dt);
            wella(end, dataIndex(i)) = wella(end-1, dataIndex(i));
        end
    end
    
end
function subWellData = bsExtractWellDataByHorizon(wellLog, horizon, dataIndex, timeIndex, upNum, downNum, filtCoef)
%% Extract target data of length (upNum+sampNum) from welllog data by given horizon information
% Programmed by: Bin She (Email: bin.stepbystep@gmail.com)
% Programming dates: Nov 2019
% -------------------------------------------------------------------------
    data = wellLog(:, dataIndex);
    [nData, nAtt] = size(data);
    
    for i = 1 : nAtt
        data(:, i) = bsButtLowPassFilter(data(:, i), filtCoef);
    end
    
    time = wellLog(:, timeIndex);
    dist = horizon - time;
    [~, index] = min(abs(dist));
    iPos = index - upNum;
    sampNum = upNum + downNum;
    
    subWellData = zeros(sampNum, nAtt);
%                 t0 = horizon(i) - GPostInvParam.upNum * GPostInvParam.dt;
%                 iPos = round((wellLogs{j}.t0 - t0) / GPostInvParam.dt);

    if iPos < 0
        % start time of origianl welllog data is below the
        % start time of the profile to be shown
        sPos = 1;
        lsPos = abs(iPos) + 1;
    else
        sPos = iPos + 1;
        lsPos = 1;
    end

    if iPos + sampNum > nData
        % end time of origianl welllog data is above the
        % end time of the profile to be shown
        ePos = nData;
        lePos = nData - iPos;
    else
        ePos = iPos+sampNum;
        lePos = sampNum;
    end

    subWellData(lsPos : lePos, :) = data(sPos : ePos, :);
    subWellData(1:lsPos, :) = repmat(data(sPos, :), lsPos, 1);
    subWellData(lePos:end, :) = repmat(data(ePos, :), sampNum-lePos+1, 1);
   
end
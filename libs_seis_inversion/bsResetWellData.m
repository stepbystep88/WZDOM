function reCreateWelllogs = bsResetWellData(GInvParam, timeLine, wellLogs, invResult)
    inIds = invResult.inIds;
    crossIds = invResult.crossIds;
    
    [wellPos, wellIndex, wellNames] = bsFindWellLocation(wellLogs, inIds, crossIds);
    reCreateWelllogs = wellLogs;

    for i = 1 : length(wellIndex)
        iWell = wellIndex(i);
        wellInfo = reCreateWelllogs{iWell};

        [horizonTimes] = bsGetHorizonTime(timeLine{GInvParam.usedTimeLineId}, wellInfo.inline, wellInfo.crossline);
        dataIndex = [GInvParam.indexInWellData.vp, GInvParam.indexInWellData.vs, GInvParam.indexInWellData.rho];

        startTime = horizonTimes - GInvParam.upNum * GInvParam.dt;

        timeData = wellInfo.wellLog(:, GInvParam.indexInWellData.time);

        for j = 1 : 3
            x = invResult.data{j}(:, wellPos(i));
            y = bsGetWellData(GInvParam, {wellInfo}, horizonTimes, dataIndex(j), 1);
            newData = bsMixTwoSignal(x, y, 10, 50, GInvParam.dt/1000);
    %         newData = x;

            [~, ps] = min(abs(startTime - timeData));
            pe = ps + GInvParam.upNum + GInvParam.downNum - 1;

            if pe > size(timeData, 1)
                pe = size(timeData, 1);
            end

            wellInfo.wellLog(ps:pe, dataIndex(j)) = newData(1:pe-ps+1, :);


        end

        reCreateWelllogs{iWell} = wellInfo;
    end
end
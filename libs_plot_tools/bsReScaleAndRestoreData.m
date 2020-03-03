function [newProfileData] = bsReScaleAndRestoreData(basicInfo, profileData, isParallel)

    
    [sampNum, ~] = size(profileData);
    minTime = basicInfo.minTime;
    newTime0 = basicInfo.newTime0;
    newT = basicInfo.newTimeSeq;
    newDt = basicInfo.newDt;
    newSampNum = round(basicInfo.scaleFactor*sampNum);
    newTrNum = length(basicInfo.newTraceIds);

    [X, Y] = meshgrid(basicInfo.traceIds, 1:sampNum);
    [Xq,Yq] = meshgrid(basicInfo.newTraceIds, linspace(1, sampNum, newSampNum));
    
    if length(basicInfo.newTraceIds) == length(basicInfo.traceIds)
        Z = profileData;
    else
        Z = interp2(X, Y, profileData, Xq, Yq,'cubic');
    end
    
    newProfileData = inf(length(newT), newTrNum);
    for i = 1 : newTrNum
        s = round((newTime0(i) - minTime) / newDt);
        newProfileData(s:s+newSampNum-1, i) = Z(:, i);
    end
    
end
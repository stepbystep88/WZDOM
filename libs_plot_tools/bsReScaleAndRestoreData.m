function [newProfileData] = bsReScaleAndRestoreData(basicInfo, profileData, isParallel)

    time = basicInfo.timeGrid;
    newT = basicInfo.newTimeSeq;
    newSampNum = length(newT);
    
    [X, Y] = meshgrid(basicInfo.traceIds, newT);
    Z = inf(size(X));
    [Xq,Yq] = meshgrid(basicInfo.newTraceIds, newT);
    
    
    
    [sampNum, traceNum] = size(profileData);
    
    if isParallel
        parfor j = 1 : traceNum
            for i = 1 : newSampNum
                ti = newT(i);

                if ti >= time(1, j) && ti<= time(sampNum, j)
                    Z(i, j) = bsCalVal(ti, time(:, j), profileData(:, j));
                end
            end
        end
    else
        for j = 1 : traceNum
            for i = 1 : newSampNum
                ti = newT(i);

                if ti >= time(1, j) && ti<= time(sampNum, j)
                    Z(i, j) = bsCalVal(ti, time(:, j), profileData(:, j));
                end
            end
        end
    end
    
    if length(basicInfo.newTraceIds) == length(basicInfo.traceIds)
        newProfileData = Z;
    else
        newProfileData = interp2(X, Y, Z, Xq, Yq,'cubic');
    end
end
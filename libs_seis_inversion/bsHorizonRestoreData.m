function [newProfileData, minTime] = bsHorizonRestoreData(profileData, horizon, upNum, dt, isZeroMinTime)

    if ~exist('isZeroMinTime', 'var')
        isZeroMinTime = 0;
    end
    
    [sampNum, trNum] = size(profileData);
    
    % fill data based on horizon
    time0 = horizon - upNum * dt;
    if isZeroMinTime 
        minTime = 0;
    else
        minTime = min(time0) - 5 * dt;
    end
    maxTime = max(time0) + sampNum * dt + 5 * dt;
    newSampNum = round((maxTime - minTime)/dt);
    
    poses = round((time0 - minTime)/dt);
    newProfileData = zeros(newSampNum, trNum);
    newProfileData(:) = nan;
    
    for i = 1 : trNum
        newProfileData(poses(i):poses(i)+sampNum-1, i) = profileData(:, i);
    end
    
end
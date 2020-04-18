function [newProfileData, minTime] = bsHorizonRestoreData(profileData, horizon, upNum, dt, minTime, isShowProgress)

    if ~exist('isShowProgress', 'var')
        isShowProgress = 0;
    end
    
    [sampNum, trNum] = size(profileData);
    
    % fill data based on horizon
    time0 = horizon - upNum * dt;
    if ~exist('minTime', 'var')
        minTime = min(time0) - 5 * dt;
    end
    maxTime = max(time0) + sampNum * dt + 5 * dt;
    newSampNum = round((maxTime - minTime)/dt);
    
    poses = round((time0 - minTime)/dt);
    newProfileData = zeros(newSampNum, trNum);
    newProfileData(:) = nan;
    
    for i = 1 : trNum
        if isShowProgress && mod(i, 10000) == 0
            % print information
            fprintf('Fill data by horion information: %d%%...\n', round(i/trNum*100));
        end
        
        try
            newProfileData(poses(i)+1:poses(i)+sampNum, i) = profileData(:, i);
        catch
            i;
        end
    end
    
end
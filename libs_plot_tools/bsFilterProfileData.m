function profileData = bsFilterProfileData(profileData, showFiltCoef)
    [sampNum, ~] = size(profileData);
    
    % filter data along with horizon
    if showFiltCoef > 0 && showFiltCoef < 1
        try
            for i = 1 : sampNum
                profileData(i, :) = bsButtLowPassFilter(profileData(i, :), showFiltCoef);
            end
        catch
        end
    end
end
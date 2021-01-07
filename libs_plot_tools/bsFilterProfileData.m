function profileData = bsFilterProfileData(profileData, showFiltCoef, isHorizonFilter)

    if nargin < 3
        isHorizonFilter = 1;
    end
    
    
    
    if ~isHorizonFilter
        profileData = profileData';
    end
    
    [sampNum, ~] = size(profileData);
    
    % filter data along with horizon
    if showFiltCoef > 0 && showFiltCoef < 1
        [b, a] = butter(10, showFiltCoef, 'low');
        for i = 1 : sampNum
            if mod(i, 10000) == 0
                fprintf('Filtering data progress information: %d/%d...\n', i, sampNum);
            end  
            try
                profileData(i, :) = filtfilt(b, a, profileData(i, :));
%                     profileData(i, :) = bsButtLowPassFilter(profileData(i, :), showFiltCoef);
            catch
                warning('Filtering failed.');
            end
        end
    end
    
    if ~isHorizonFilter
        profileData = profileData';
    end
    
end
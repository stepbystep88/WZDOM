function [horizonTimes] = bsCalcHorizonTime(usedTimeLine, inIds, crossIds, isPrintInfo)
%% get the horizon of given traces

    if nargin < 4
        isPrintInfo = 0;
    end
    
    nTrace = length(inIds);
    horizonTimes = zeros(1, nTrace);
    
    for i = 1 : nTrace
        [~, ~, horizonTimes(i)] = bsCalcWellBaseInfo(usedTimeLine, ...
            inIds(i), crossIds(i), 1, 2, 1, 2, 3);
        
        percent = i/nTrace*100;
        if isPrintInfo
            if mod(i, 100) == 0
                fprintf('Calculating horizon information %.2f%%\n', percent);
            end
        end
        
    end
end
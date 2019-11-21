function [horizonTimes] = bsCalcHorizonTimeParallel(usedTimeLine, inIds, crossIds, numWorkers)
%% get the horizon of given traces
    
    nTrace = length(inIds);
    horizonTimes = zeros(1, nTrace);
    
    pbm = bsInitParforProgress(numWorkers, nTrace, 'Calculating horizon information');
            
    
    parfor i = 1 : nTrace
        try
            [~, ~, horizonTimes(i)] = bsCalcWellBaseInfo(usedTimeLine, ...
                inIds(i), crossIds(i), 1, 2, 1, 2, 3);
        catch
            error('%d trace is failed!!!', i);
        end

        bsIncParforProgress(pbm, i, 100);
        
    end
    
    
    
end
function [horizonTimes] = bsCalcHorizonTime(usedTimeLine, inIds, crossIds, isParallel, numWorkers)
%% get the horizon of given traces

    if ~exist('isParallel', 'var')
        isParallel = 0;
    end
    
    nTrace = length(inIds);
    horizonTimes = zeros(1, nTrace);
    
    if isParallel
        
        if ~exist('numWorkers', 'var')
            numWorkers = bsGetMaxNumWorkers();
        end
        
        pbm = bsInitParforProgress(numWorkers, nTrace, 'Calculating horizon information', [], 0);
            
        parfor i = 1 : nTrace
            try
                [~, ~, horizonTimes(i)] = bsCalcWellBaseInfo(usedTimeLine, ...
                    inIds(i), crossIds(i), 1, 2, 1, 2, 3);
            catch
                error('%d trace is failed!!!', i);
            end

            bsIncParforProgress(pbm, i, 10000);

        end
    else
        for i = 1 : nTrace
            try
                [~, ~, horizonTimes(i)] = bsCalcWellBaseInfo(usedTimeLine, ...
                    inIds(i), crossIds(i), 1, 2, 1, 2, 3);
            catch
                error('%d trace is failed!!!', i);
            end

        end
    end
    
    
    
end
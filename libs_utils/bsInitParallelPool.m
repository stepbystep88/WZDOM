function numWorkers = bsInitParallelPool(numWorkers)
    
    p = gcp('nocreate');
    
    if isempty(numWorkers)
        if ~isempty(p)
            
            numWorkers = p.NumWorkers;
        end
    end
    
    if isempty(p)
        if isempty(numWorkers)
            parpool('local');
        else
            parpool('local', numWorkers);
        end
    elseif p.NumWorkers < numWorkers
        delete p;
        parpool('local', numWorkers);
    end
    
    p = gcp('nocreate');
    numWorkers = p.NumWorkers;
end
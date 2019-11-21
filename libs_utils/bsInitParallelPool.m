function numWorkers = bsInitParallelPool(numWorkers)
    
    p = gcp('nocreate');
    
    if isempty(p) 
        parpool('local', numWorkers);
    elseif p.NumWorkers < numWorkers
        delete p;
        parpool('local', numWorkers);
    end
    
    p = gcp('nocreate');
    numWorkers = p.NumWorkers;
end
function pbm = bsInitParforProgress(numWorkers, nLoop, title, basePath, isPrintBySavingFile)

    bsInitParallelPool(numWorkers);
    
    if ~exist('title', 'var')
        pbm.title = 'Progress information';
    else
        pbm.title = title;
    end
    
    if ~exist('isPrintBySavingFile', 'var')
        pbm.isPrintBySavingFile = 0;
    else
        pbm.isPrintBySavingFile = isPrintBySavingFile;
    end
    
    pbm.name = [basePath, sprintf('/parfor_progress_%d.txt', round(rand()*10000))];
    p = gcp('nocreate');
    pbm.numWorkers = p.NumWorkers;
    pbm.nLoop = nLoop;
    pbm.fid = 0;
    
    
    if isPrintBySavingFile
        fid = fopen(pbm.name, 'w');
        fclose(fid);
    end
end
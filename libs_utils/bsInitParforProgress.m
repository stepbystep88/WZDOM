function pbm = bsInitParforProgress(numWorkers, nLoop, title)

    bsInitParallelPool(numWorkers);
    
    if ~exist('title', 'var')
        pbm.title = 'Progress information';
    else
        pbm.title = title;
    end
    
    pbm.name = './parfor_progress.txt';
    p = gcp('nocreate');
    pbm.nWorkers = p.NumWorkers;
    pbm.nLoop = nLoop;
    pbm.fid = 0;
    
%     parfor_progress = 0;
%     assignin('base', 'parfor_progress', parfor_progress);
    fid = fopen(pbm.name, 'w');
    fclose(fid);
end
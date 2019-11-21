function pbm = bsInitParforProgress(nLoop, title)

    if ~exist('title', 'var')
        pbm.title = 'Progress information';
    else
        pbm.title = title;
    end
    
    pbm.name = './parfor_progress.txt';
    p = gcp('nocreate');
    pbm.nWorkers = p.NumWorkers;
    pbm.nLoop = nLoop;
    pbm.fids = zeros(1, p.NumWorkers);
    
    fid = fopen(pbm.name, 'w');
    fclose(fid);
    
end
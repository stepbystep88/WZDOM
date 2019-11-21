function bsPrintParforProgress(pbm)

%     s = dir(pbm.name);
%     the_size = s.bytes;
    id = get(getCurrentTask(), 'ID');
    
    parfor_progress = evalin('base','parfor_progress');
    nFinished = parfor_progress * pbm.numWorkers;
    fprintf('[%s]: %.2f%% of %d loops (told by worker [%d])...\n', ...
        pbm.title, nFinished/pbm.nLoop*100, pbm.nLoop, id);
    
end
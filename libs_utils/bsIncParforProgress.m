function bsIncParforProgress(pbm, iTrace, interval)
    
    
    fid = fopen(pbm.name, 'a');
    fwrite(fid, uint8(0), 'uint8');
    fclose(fid);
    
%     parfor_progress = evalin('base','parfor_progress');
%     parfor_progress = parfor_progress + 1;
%     assignin('base', 'parfor_progress', parfor_progress);
    
    if mod(iTrace, interval) == 0
        id = get(getCurrentTask(), 'ID');

%         nFinished = parfor_progress * pbm.numWorkers;
        s = dir(pbm.name);
        nFinished = s.bytes;

        fprintf('[%s]: %.2f%% of %d loops (told by worker [%d])...\n', ...
            pbm.title, nFinished/pbm.nLoop*100, pbm.nLoop, id);
    end
    
end
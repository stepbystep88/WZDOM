function bsIncParforProgress(pbm, iTrace, interval)
    
    
    if mod(iTrace, interval) == 0

        if pbm.isPrintBySavingFile
            fid = fopen(pbm.name, 'a');
            fwrite(fid, uint8(0), 'uint8');
            fclose(fid);
        
            id = get(getCurrentTask(), 'ID');

            s = dir(pbm.name);
            nFinished = s.bytes * interval;

            fprintf('[%s]: %.2f%% of %d loops (told by worker [%d])...\n', ...
                pbm.title, nFinished/pbm.nLoop*100, pbm.nLoop, id);
        else
            fprintf('[%s]: dealing with the %d-th loop of %d loops ...\n', ...
                pbm.title, iTrace, pbm.nLoop);
        end
    end
    
end
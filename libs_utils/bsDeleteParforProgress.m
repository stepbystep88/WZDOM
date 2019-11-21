function bsDeleteParforProgress(pbm)
    
    for i = 1 : length(pbm.fids)
        if pbm.fids(i) > 0
            try
                fclose(pbm.fids(i));
            catch
            end
        end
    end
end
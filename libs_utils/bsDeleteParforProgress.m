function bsDeleteParforProgress(pbm)
    
    try
        if pbm.isPrintBySavingFile
            delete(pbm.name);
        end
    catch
    end
    
end
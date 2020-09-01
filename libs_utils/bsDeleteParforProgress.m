function bsDeleteParforProgress(pbm)
    warning('off');
    try
        if pbm.isPrintBySavingFile
            delete(pbm.name);
        end
    catch
    end
    warning('on');
    
end
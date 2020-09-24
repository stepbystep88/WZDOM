function bsDeleteParforProgress(pbm)
    warning('off');
    try
%         if pbm.isPrintBySavingFile
%             delete(pbm.name);
%         end
        delete(pbm.name);
    catch
    end
    warning('on');
    
end
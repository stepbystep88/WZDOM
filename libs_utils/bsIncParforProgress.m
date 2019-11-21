function pbm = bsIncParforProgress(pbm)
    
%     id = get(getCurrentTask(), 'ID');
    
    if pbm.fid == 0
        pbm.fid = fopen(pbm.name, 'a');
    end
    
    fwrite(fid, uint8(0), 'uint8');
    
end
function bsIncParforProgress(pbm)
    
    id = get(getCurrentTask(), 'ID');
    if pbm.fids(id) == 0
        pbm.fids(id) = fopen(pbm.name, 'a');
        fwrite(pbm.fids(id), uint8(0), 'uint8');
    end
end
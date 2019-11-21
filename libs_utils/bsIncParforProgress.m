function pbm = bsIncParforProgress(pbm)
    
%     id = get(getCurrentTask(), 'ID');
    fid = fopen(pbm.name, 'a');
    fwrite(fid, uint8(0), 'uint8');
    fclose(fid);
end
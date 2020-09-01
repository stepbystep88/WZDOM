function pbm = bsResetParforProgress(pbm, title)
    bsDeleteParforProgress(pbm);
    pbm.title = title;
    
    try
        fid = fopen(pbm.name, 'w');
        fclose(fid);
    catch
    end
end
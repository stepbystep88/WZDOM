function bsPrintParforProgress(pbm)

    s = dir(pbm.name);
    the_size = s.bytes;

    fprintf('\n[%s]: %.2f%% of %d loops..\n', pbm.title, the_size/pbm.nLoop*100, pbm.nLoop);
    
end
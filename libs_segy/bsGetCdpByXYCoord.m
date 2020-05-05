function [inIds, crossIds] = bsGetCdpByXYCoord(Xs, Ys, infos)
  
    nTrace = length(Xs);
    inIds = zeros(1, nTrace);
    crossIds = zeros(1, nTrace);
    
    for i = 1 : nTrace
        dist = (Xs(i) - infos(:, 3)).^2 + (Ys(i) - infos(:, 4)).^2;
        [~, index] = min(dist);
        
        inIds(i) = infos(index, 1);
        crossIds(i) = infos(index, 2);
    end
    
end
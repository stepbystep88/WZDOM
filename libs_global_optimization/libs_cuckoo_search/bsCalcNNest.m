function nNest = bsCalcNNest(minN, maxN, maxNFE, NFE)
    % this funciton is set for population reduction strategy
    nNest = round((minN - maxN)/maxNFE*NFE + maxN);
end
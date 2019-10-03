function [nests, gradient, fitness] = bsSortNests(nests, gradient, fitness, nNewNest)

    [nDim, nNest] = size(nests);
    % the first dimension is the order information
    nDim = nDim - 1;
    
    piFi = [nests', gradient', fitness];
    % Sort by Fi in assending order
    piFiS = sortrows(piFi, nDim*2+2);

    % Decrease number of nests, we only need lots of nests initially to get
    % a good sampling of the objective function
    nests = piFiS(1:nNewNest, 1:nDim+1)';
    gradient = piFiS(1:nNewNest, end-nDim:end-1)';
    fitness = piFiS(1:nNewNest, end);

        
end


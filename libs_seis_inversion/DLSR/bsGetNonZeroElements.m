function [gamma_vals, gamma_locs, gamma_K] = bsGetNonZeroElements(gammas, sparsity, Idx)
% 将单道多个小块的非零元素全部获取到并列成一个vector

    ncell = size(gammas, 2);
    gamma_vals = zeros(sparsity*ncell, 1);
    gamma_locs = zeros(sparsity*ncell, 1) - 1;
    gamma_K = zeros(sparsity*ncell, 1);
    
    sPos = 1;
    for i = 1 : ncell
        
        index = find(gammas(:, i) ~= 0);
        ePos = sPos + length(index) - 1;
        
        gamma_vals(sPos:ePos) = gammas(index, i);
        gamma_locs(sPos:ePos) = index;
        
        if(nargin >= 3)
            gamma_K(i) = Idx(index);
        end
        
        sPos = sPos + sparsity;
    end
end
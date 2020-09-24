function [nlm_ps] = bsGetNonLocalSimilarPatches(self, data, sizeAtom, K)
    [sampNum, nTrace] = size(data);
    ncell = sampNum - sizeAtom + 1;
    
    D = inf(ncell * nTrace, 3);
    
    
    v_ones = ones(1, nTrace);
    v_lin = 1 : nTrace;
    nlm_ps = zeros(K, ncell);
    
    for k = 1 : ncell
        kpatch = repmat(self(k : k+sizeAtom-1), 1, nTrace);
        sp = 1;
        for j = 1 : ncell
            cmp_patches = data(j:j+sizeAtom-1, :);
            
            residual = kpatch - cmp_patches;
            
            ep = sp + nTrace - 1;
            D(sp:ep, 1) = sum(residual.^2, 1);
            D(sp:ep, 3) = v_ones * j;
            D(sp:ep, 2) = v_lin;
            sp = ep + 1;
        end
        
        [B, ~] = bsMinK(D, K, 1);
        nlm_ps(:, k) = (B(:, 2) - 1) * ncell + B(:, 3);
    end
   
    
%     weight = B(:, 1).^e;
%     weight = weight / sum(weight);
%     index = B(:, 2);
    
end
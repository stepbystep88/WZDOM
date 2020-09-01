function [weight, index] = bsGetNonLocalSimilarPatches(self, data, sizeAtom, K)
    [sampNum, nTrace] = size(data);
    ncell = sampNum - sizeAtom + 1;
    
    D = inf(ncell * nTrace, 2);
    iter = 0;

    for k = 1 : ncell
        kpatch = repmat(self(k : k+sizeAtom-1), 1, nTrace);
        for j = 1 : ncell
            cmp_patches = data(j:j+sizeAtom-1, :);
            
            residual = kpatch - cmp_patches;
            
            D(k, :, 1) = sum(residual.^2, 1);
        end
    end
   
    bsSparseTransInput2Patches
        for c2 = i2-N : searchStride(2) : i2+N
            if c2<=0 || c2>n2 || (c2==i2 && c1==i1)
                continue;
            end

            iter = iter + 1;

            j = (c1-1)*n2 + c2;

            
            D(iter, 1) = norm(difference/normCoef);
            D(iter, 2) = j;
        end
    end

    [B, ~] = bsMinK(D, K, 1);
    weight = B(:, 1).^e;
    weight = weight / sum(weight);
    index = B(:, 2);
end
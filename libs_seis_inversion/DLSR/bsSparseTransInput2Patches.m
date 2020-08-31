function [all_patches] = bsSparseTransInput2Patches(GSParam, input, inline, crossline)
    
    [~, nBlock] = size(input{1});
    nAtt = length(input);
    
    sizeAtom = GSParam.sizeAtom;
    ncell = GSParam.ncell;
    
    all_patches = zeros(GSParam.nSpecialFeat*nBlock+sizeAtom*nBlock*nAtt, ncell);
    patches = zeros(sizeAtom, ncell);
    

    sPos = 1;
    for k = 1 : nBlock
        ePos = sPos + GSParam.nSpecialFeat - 1;
        % add location or/and time information
        if GSParam.trainDICParam.isAddLocInfo && GSParam.trainDICParam.isAddTimeInfo
            all_patches(sPos : ePos, :) = [ones(1, ncell) * inline(k); ones(1, ncell) * crossline(k); 1 : ncell];
        elseif GSParam.trainDICParam.isAddLocInfo
            all_patches(sPos : ePos, :) = [ones(1, ncell) * inline(k); ones(1, ncell) * crossline(k);];
        elseif GSParam.trainDICParam.isAddTimeInfo
            all_patches(sPos : ePos, :) = 1 : ncell;
        end
        
        sPos = ePos + 1;
    
        for i = 1 : nAtt
            for j = 1 : ncell
                js = GSParam.index(j);
                patches(:, j) = input{i}(js : js+sizeAtom-1, k);
            end
            
            ePos = sPos + sizeAtom - 1;
            all_patches(sPos:ePos, :) = patches;
            sPos = sPos + sizeAtom;
        end
    end
end
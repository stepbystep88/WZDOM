function [all_patches] = bsSparseTransInput2PatchesByNLM(GSParam, input, nlm_ps, inline, crossline)
    
    nBlock = size(nlm_ps, 1);
    nAtt = length(input);
    nTrace = size(input{1}, 2);
    
    sizeAtom = GSParam.sizeAtom;
    ncell = GSParam.ncell;
    
    all_patches = zeros(GSParam.nSpecialFeat*nBlock+sizeAtom*nBlock*nAtt, ncell);
    
    for j = 1 : ncell
        sPos = 1;
        
        for k = 1 : nBlock

            iTrace = floor((nlm_ps(k, j) + ncell - 1) / ncell);
            iPatch = nlm_ps(k, j) - (iTrace - 1) * ncell;

            ePos = sPos + GSParam.nSpecialFeat - 1;
            % add location or/and time information
            if GSParam.trainDICParam.isAddLocInfo && GSParam.trainDICParam.isAddTimeInfo
                all_patches(sPos : ePos, :) = [ones(1, ncell) * inline(iTrace); ones(1, ncell) * crossline(iTrace); 1 : ncell];
            elseif GSParam.trainDICParam.isAddLocInfo
                all_patches(sPos : ePos, :) = [ones(1, ncell) * inline(iTrace); ones(1, ncell) * crossline(iTrace);];
            elseif GSParam.trainDICParam.isAddTimeInfo
                all_patches(sPos : ePos, :) = 1 : ncell;
            end

            sPos = ePos + 1;
            ePos = sPos + sizeAtom - 1;
            
            for i = 1 : nAtt
                all_patches(sPos:ePos, j) = input{i}(iPatch : iPatch+sizeAtom-1, iTrace);
            end
            
            sPos = sPos + sizeAtom;
        end
    end
end
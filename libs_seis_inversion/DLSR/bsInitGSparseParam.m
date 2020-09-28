function [GSparseParam] = bsInitGSparseParam(GSparseParam, sampNum, nBlock, lowIndex, highIndex)
    GSparseParam.sizeAtom = GSparseParam.trainDICParam.sizeAtom;
    GSparseParam.nAtom = GSparseParam.trainDICParam.nAtom;
    
    GSparseParam.nrepeat = GSparseParam.sizeAtom - GSparseParam.stride;
    
    index = 1 : GSparseParam.stride : sampNum - GSparseParam.sizeAtom + 1;
    if(index(end) ~= sampNum - GSparseParam.sizeAtom + 1)
        index = [index, sampNum - GSparseParam.sizeAtom + 1];
    end
    
    GSparseParam.nSpecialFeat = GSparseParam.trainDICParam.isAddLocInfo * 2 + GSparseParam.trainDICParam.isAddTimeInfo;
    GSparseParam.index = index;
    GSparseParam.ncell = length(index);
    
    
    
    GSparseParam.ODIC = repmat(GSparseParam.DIC, nBlock, 1);
    if nBlock > 1
        for i = 1 : size(GSparseParam.ODIC, 2)
            tmp = norm(GSparseParam.ODIC(:, i));
            GSparseParam.ODIC(:, i) = GSparseParam.ODIC(:, i) / tmp;
            GSparseParam.DIC(:, i) = GSparseParam.DIC(:, i) / tmp;
        end
    end
    
    GSparseParam.omp_G = GSparseParam.ODIC' * GSparseParam.ODIC;
    
    rangeCoef = GSparseParam.rangeCoef;
    ncell = GSparseParam.ncell;
    
    if ~isempty(rangeCoef)
        GSparseParam.min_values = repmat(rangeCoef(:, 1), nBlock, ncell);
        GSparseParam.max_values = repmat(rangeCoef(:, 2), nBlock, ncell);

%         GSparseParam.min_values_one = repmat(rangeCoef(:, 1), 1, ncell);
%         GSparseParam.max_values_one = repmat(rangeCoef(:, 2), 1, ncell);
    else
        GSparseParam.min_values = [];
        GSparseParam.max_values = [];

    end
    
    if exist('highIndex', 'var')
        
        sPos = GSparseParam.nSpecialFeat + (highIndex - 1)*GSparseParam.sizeAtom + 1;
        ePos = sPos + GSparseParam.sizeAtom - 1;
        
        highs = sPos : ePos;
        
        if exist('lowIndex', 'var') && ~isempty(lowIndex)
            sPos = GSparseParam.nSpecialFeat + (lowIndex - 1)*GSparseParam.sizeAtom + 1;
            ePos = sPos + GSparseParam.sizeAtom - 1;

            lows = sPos : ePos;
        else
            lows = setdiff(1:size(GSparseParam.DIC, 1), highs);
        end
        
        
        if strcmpi(GSparseParam.trainDICParam.feature_reduction, 'off')
            D1 = GSparseParam.DIC(lows, :);
            D2 = GSparseParam.DIC(highs, :);
            
        else
            nfeat = size(GSparseParam.output.B, 2);
            D1 = GSparseParam.DIC(1:nfeat, :);
            D2 = GSparseParam.DIC(nfeat+1:end, :);
            
            B = cell(nBlock, nBlock);
            for i = 1 : nBlock
                for j = 1 : nBlock
                    if i == j
                        B{i, i} = GSparseParam.output.B;
                    else
                        B{i, j} = zeros(size(GSparseParam.output.B));
                    end
                end
            end
            
            GSparseParam.output.B = sparse(cell2mat(B));
        end

        D1 = repmat(D1, nBlock, 1);

        [D1, D2] = bsNormalDIC(D1, D2);
        
        GSparseParam.lowDIC = D1;
        GSparseParam.highDIC = D2;
        GSparseParam.omp_low_G = D1' * D1;
        
        if ~isempty(rangeCoef)
            GSparseParam.low_min_values = repmat(rangeCoef(lows, 1), nBlock, ncell);
            GSparseParam.low_max_values = repmat(rangeCoef(lows, 2), nBlock, ncell);

            GSparseParam.high_min_values = repmat(rangeCoef(highs, 1), 1, ncell);
            GSparseParam.high_max_values = repmat(rangeCoef(highs, 2), 1, ncell);
        else
            GSparseParam.low_min_values = [];
            GSparseParam.low_max_values = [];
            GSparseParam.high_min_values = [];
            GSparseParam.high_max_values = [];
        end
    end
    
    
end

function [D1, D2, C] = bsNormalDIC(D1, D2)
    C = zeros(size(D1, 2), 1);
    
    for i = 1 : size(D1, 2)
        normCoef = norm(D1(:, i));
        D1(:, i) = D1(:, i) / normCoef;
        D2(:, i) = D2(:, i) / normCoef;
        C(i) = normCoef;
    end
end
function [DIC, rangeCoef] = bsTrain1DSparseJointDIC(datas, GTrainDICParam)
%% To train a joint dictionary by using K-SVD algorithm
%
% Programmed by: Bin She (Email: bin.stepbystep@gmail.com)
% Programming dates: Dec 2019
% -------------------------------------------------------------------------

    nWell = length(datas);
    allPatchs = [];
    
    nAtt = size(datas{1}, 2);
    rangeCoef = zeros(nAtt, 2);
    rangeCoef(:, 1) = inf;
    rangeCoef(:, 2) = -inf;
    
    % get the maximum and minimum value of different attributes
    for j = 1 : nWell
        data = datas{j};
        
        for i = 1 : nAtt
            if( GTrainDICParam.filtCoef < 1)
                data(:, i) = bsButtLowPassFilter(data(:, i), GTrainDICParam.filtCoef);
            end
            rangeCoef(i, 1) = min(rangeCoef(i, 1), min(data(:, i)));
            rangeCoef(i, 2) = max(rangeCoef(i, 2), max(data(:, i)));
        end
        
        datas{j} = data;
    end
    
    for i = 1 : nAtt
        patchs = [];

        for j = 1 : nWell
            data = datas{j};
            num = size(data, 1);
            index = 1 : GTrainDICParam.stride : num - GTrainDICParam.sizeAtom + 1;
         
            for k = index
                subData = data(k : k+GTrainDICParam.sizeAtom-1, i);
                
                subData = (subData - rangeCoef(i, 1)) / (rangeCoef(i, 2) - rangeCoef(i, 1));

                patchs = [patchs, subData];
            end
        end
        
        allPatchs = [allPatchs; patchs];
        
    end
    
    params.data = allPatchs;
    params.Tdata = GTrainDICParam.sparsity;
    params.dictsize = GTrainDICParam.nAtom;
    params.iternum = GTrainDICParam.iterNum;
    params.memusage = 'high';

    % using the third-party toolbox to train the dictionary
    DIC = ksvd(params, GTrainDICParam.isShowIterInfo);
end

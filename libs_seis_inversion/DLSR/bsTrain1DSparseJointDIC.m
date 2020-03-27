function [DIC, rangeCoef, output] = bsTrain1DSparseJointDIC(datas, GTrainDICParam, xs, ys)
%% To train a joint dictionary by using K-SVD algorithm
%
% Programmed by: Bin She (Email: bin.stepbystep@gmail.com)
% Programming dates: Dec 2019
% -------------------------------------------------------------------------
    
    nWell = length(datas);
    allPatchs = [];
    
    nAtt = size(datas{1}, 2);
    
    % 预处理
    for j = 1 : nWell
        data = datas{j};
        
        for i = 1 : nAtt
            if( GTrainDICParam.filtCoef < 1)
                data(:, i) = bsButtLowPassFilter(data(:, i), GTrainDICParam.filtCoef);
            end
        end
        
        datas{j} = data;
    end
    
    validatestring(string(GTrainDICParam.normalizationMode), {'whole_data_max_min', 'feat_max_min', 'off'});
    
    % 小块拼接
    AP = cell(1, nAtt);
    nPatch = 0;
    for j = 1 : nWell
        nPatch = nPatch + length(1 : GTrainDICParam.stride : size(datas{j}, 1) - GTrainDICParam.sizeAtom + 1);
    end
    
    for i = 1 : nAtt
        n = GTrainDICParam.sizeAtom;
        
        if i == 1 
            n = n + GTrainDICParam.isAddLocInfo*2 + GTrainDICParam.isAddTimeInfo;
        end
        
        AP{i} = zeros(n, nPatch);
        counter = 1;
        
        for j = 1 : nWell
            data = datas{j};
            num = size(data, 1);
            index = 1 : GTrainDICParam.stride : num - GTrainDICParam.sizeAtom + 1;
         
            for k = index
                patch = data(k : k+GTrainDICParam.sizeAtom-1, i);
                
                if i == 1
                    if GTrainDICParam.isAddLocInfo && GTrainDICParam.isAddTimeInfo
                        subData = [patch; xs(j); ys(j); k];
                    elseif GTrainDICParam.isAddLocInfo
                        subData = [patch; xs(j); ys(j)];
                    elseif GTrainDICParam.isAddTimeInfo
                        subData = [patch; k];
                    else
                        subData = patch;
                    end
                else
                    subData = patch;
                end
                
                
                AP{i}(:, counter) = subData;
                counter = counter + 1;
            end
        end
    end
    
    % 归一化处理
    switch GTrainDICParam.normalizationMode
        case 'whole_data_max_min'
            rangeCoef = zeros(nAtt, 2);
            for i = 1 : nAtt
                rangeCoef(i, 1) = min(AP{i}(:));
                rangeCoef(i, 2) = max(AP{i}(:));

                AP{i} = (AP{i} - rangeCoef(i, 1)) / (rangeCoef(i, 2) - rangeCoef(i, 1));
            end
        case 'feat_max_min'
            rangeCoef = cell(1, nAtt);
            for i = 1 : nAtt
                rangeCoef{i} = zeros(size(AP{i}, 1), 2);
                
                rangeCoef{i}(:, 1) = min(AP{i}, [], 2);
                rangeCoef{i}(:, 2) = max(AP{i}, [], 2);
                
                nPatch = size(AP{i}, 2);
                
                min_val = repmat(rangeCoef{i}(:, 1), 1, nPatch);
                max_val = repmat(rangeCoef{i}(:, 2), 1, nPatch);
                
                AP{i} = (AP{i} - min_val) ./ (max_val - min_val);
            end
            
        case 'off'
            rangeCoef = [];
            
    end
    
   
    if GTrainDICParam.feature_reduction
        R = AP{1}*AP{1}'; 
        [B, SS] = eig(R); 
        Permute = fliplr(eye(size(R, 1))); 
        SS = Permute * SS * Permute; % so that the eigenvalues are sorted descending
        B = B * Permute; 
        energy = cumsum(diag(SS))/sum(diag(SS)); 
        % figure(1); clf; plot(energy)
        pos=find(energy>0.999, 1);
        B = B(:, 1:pos);
        AP{1} = B' * AP{1};
        output.B = B;
    else
        output = [];
    end
    
    %整合所有的patches到一起
    allPatchs = [];
    for i = 1 : nAtt
        allPatchs = [allPatchs; AP{i}];
    end
    
    params.data = allPatchs;
    params.Tdata = GTrainDICParam.sparsity;
    params.dictsize = GTrainDICParam.nAtom;
    params.iternum = GTrainDICParam.iterNum;
    params.memusage = 'high';

    % using the third-party toolbox to train the dictionary
    DIC = ksvd(params, GTrainDICParam.isShowIterInfo);
end

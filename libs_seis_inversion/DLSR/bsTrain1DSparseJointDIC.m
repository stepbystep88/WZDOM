function [DIC, rangeCoef, output] = bsTrain1DSparseJointDIC(datas, GTrainDICParam, xs, ys)
%% To train a joint dictionary by using K-SVD algorithm
%
% Programmed by: Bin She (Email: bin.stepbystep@gmail.com)
% Programming dates: Dec 2019
% -------------------------------------------------------------------------
    
    nWell = length(datas);
    allPatchs = [];
    
    nAtt = size(datas{1}, 2);
    
    % Ԥ����
    for j = 1 : nWell
        data = datas{j};
        
        for i = 1 : nAtt
            if( GTrainDICParam.filtCoef < 1)
                data(:, i) = bsButtLowPassFilter(data(:, i), GTrainDICParam.filtCoef);
            end
        end
        
        datas{j} = data;
    end
    
    validatestring(string(GTrainDICParam.normalizationMode), {'whole_data_max_min', 'whole_data_mean_sigma', 'self_sigma', ...
        'feat_max_min', 'feat_mean_sigma', 'self_mean_sigma', 'off', 'none'});
    
    % С��ƴ��
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
                        subData = [xs(j); ys(j); k; patch];
                    elseif GTrainDICParam.isAddLocInfo
                        subData = [xs(j); ys(j); patch];
                    elseif GTrainDICParam.isAddTimeInfo
                        subData = [k; patch];
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
    
    switch GTrainDICParam.normalizationMode
        case {'off', 'none'}
            
            rangeCoef = [];
        case 'feat_max_min'
            rangeCoef = [];
            
            for i = 1 : length(AP)
                min_value = min(AP{i}, [], 2);
                max_value = max(AP{i}, [], 2);
                
                rangeCoef = [rangeCoef; [min_value, max_value]];
                
                min_value = repmat(min_value, 1, size(AP{i}, 2));
                max_value = repmat(max_value, 1, size(AP{i}, 2));
                
                AP{i} = (AP{i} - min_value) ./ (max_value - min_value);
            end
        case 'feat_mean_sigma'
            rangeCoef = [];
            
            for i = 1 : length(AP)
                mean_value = mean(AP{i}, 2);
                sigma_value = var(AP{i}, 0, 2);
                
                rangeCoef = [rangeCoef; [mean_value, sigma_value]];
                
                mean_value = repmat(mean_value, 1, size(AP{i}, 2));
                sigma_value = repmat(sigma_value, 1, size(AP{i}, 2));
                
                AP{i} = (AP{i} - mean_value) ./ sigma_value;
            end
        case 'self_mean_sigma'
            rangeCoef = [];
            for i = 1 : length(AP)
                mean_value = mean(AP{i}, 1);
                sigma_value = var(AP{i}, 0, 1);
                
%                 rangeCoef = [rangeCoef; [mean_value, sigma_value]];
                
                mean_value = repmat(mean_value, size(AP{i}, 1), 1);
                sigma_value = repmat(sigma_value, size(AP{i}, 1), 1);
                
                AP{i} = (AP{i} - mean_value) ./ sigma_value;
            end
            
        case 'whole_data_max_min'
            rangeCoef = [];
            
            for i = 1 : length(AP)
                nFeat = size(AP{i}, 1);


%                 min_value = prctile(reshape(AP{i}(:),[], 1), 10);
%                 max_value = prctile(reshape(AP{i}(:),[], 1), 90);
                min_value = min(AP{i}(:));
                max_value = max(AP{i}(:));
                rangeCoef = [rangeCoef; [repmat(min_value, nFeat, 1), repmat(max_value, nFeat, 1)]];
                
                AP{i} = (AP{i} - min_value) ./ (max_value - min_value);
                
            end
        case 'whole_data_mean_sigma'
            rangeCoef = [];
            
            for i = 1 : length(AP)
                nFeat = size(AP{i}, 1);

                mean_value = mean(AP{i}(:));
                var_value = var(AP{i}(:));
                rangeCoef = [rangeCoef; [repmat(mean_value, nFeat, 1), repmat(var_value, nFeat, 1)]];
                
                AP{i} = (AP{i} - mean_value) ./ var_value;
                
            end
    end
    
    
    switch GTrainDICParam.feature_reduction
        case 'high_resolution'
            [AP{1}, output.B] = bsReduction(AP{1}, GTrainDICParam.n_reduced_feature);
            %�������е�patches��һ��
            allPatchs = [];
            for i = 1 : nAtt
                allPatchs = [allPatchs; AP{i}];
            end
    
        case 'all'
            %�������е�patches��һ��
            allPatchs = [];
            for i = 1 : nAtt
                allPatchs = [allPatchs; AP{i}];
            end

            [allPatchs, output.B] = bsReduction(allPatchs, GTrainDICParam.n_reduced_feature);
        
        case 'off'
            %�������е�patches��һ��
            allPatchs = [];
            for i = 1 : nAtt
                allPatchs = [allPatchs; AP{i}];
            end
            output = [];
    end
    
    
    
    params.data = allPatchs;
    params.Tdata = GTrainDICParam.sparsity;
    if GTrainDICParam.nAtom <= 1 && GTrainDICParam.nAtom > 0
        params.dictsize = round(size(allPatchs, 2)* GTrainDICParam.nAtom);
    else
        params.dictsize = GTrainDICParam.nAtom;
    end
    params.iternum = GTrainDICParam.iterNum;
    params.memusage = 'high';

    % using the third-party toolbox to train the dictionary
    DIC = ksvd(params, GTrainDICParam.isShowIterInfo);
    
   
    [DIC, Idx] = bsClusterDIC(DIC, 5);
%     DIC = newDIC;
end




        
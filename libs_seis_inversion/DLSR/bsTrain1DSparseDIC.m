function [DIC, rangeCoef, output] = bsTrain1DSparseDIC(datas, GTrainDICParam, xs, ys)
%% To train a dictionary by using K-SVD algorithm
%
% Programmed by: Bin She (Email: bin.stepbystep@gmail.com)
% Programming dates: Nov 2019
% -------------------------------------------------------------------------
%
% Input 
% datas                     training samples, cell data type
% GTrainDICParam.sizeAtom               	the length of an atom
% GTrainDICParam.nAtom                     the number of atoms
% GTrainDICParam.iterNum                   the number of iteration of the training process
% GTrainDICParam.stride                    the step size to cut the data into pieces
% 
%
% Output
% dic                       the trained dictionary
% initDic                   the 
% -------------------------------------------------------------------------

    nCols = length(datas);
    
    %sizeAtom, nAtom, stride, iterNum, filtCoef, sparsity, isShowIterInfo
    
    % patches with the same length which are used to train the spare
    % dictionary
    patches = [];
    rangeCoef = [];
    output = [];
    
    % generate patches
    for i = 1 : nCols
        data = datas{i};
        num = size(data, 1);
        
        if( length(GTrainDICParam.filtCoef) == 1)
            % using low-pass filter
            if(GTrainDICParam.filtCoef < 1)
                data = bsButtLowPassFilter(data, GTrainDICParam.filtCoef);
            end
        else
            % using band-pass filter
            [b, a] = butter(10, GTrainDICParam.filtCoef, 'stop');
            data = filtfilt(b, a, data);
        end
        
        
    
        % gather patches toghter
        for j = 1 : GTrainDICParam.stride : num - GTrainDICParam.sizeAtom + 1
            
            patch = data(j:j+GTrainDICParam.sizeAtom-1);
        
            if GTrainDICParam.isAddLocInfo && GTrainDICParam.isAddTimeInfo
                subData = [xs(i); ys(i); j; patch];
            elseif GTrainDICParam.isAddLocInfo
                subData = [xs(i); ys(i); patch];
            elseif GTrainDICParam.isAddTimeInfo
                subData = [j; patch];
            else
                subData = patch;
            end
        
            patches = [patches, subData];
        end
    end
   
    switch GTrainDICParam.normalizationMode
        case 'off'
            rangeCoef = [];
            
        case 'feat_max_min'
            min_value = min(patches, [], 2);
            max_value = max(patches, [], 2);

            rangeCoef = [min_value, max_value];

            min_value = repmat(min_value, 1, size(patches, 2));
            max_value = repmat(max_value, 1, size(patches, 2));

            patches = (patches - min_value) ./ (max_value - min_value);
    end
    
    if ~strcmp(GTrainDICParam.feature_reduction, 'off')
        [new_patches, output.B] = bsReduction(patches, GTrainDICParam.n_reduced_feature);
    else
        new_patches = patches;
    end
    
    
    params.data = new_patches;
    params.Tdata = GTrainDICParam.sparsity;
%     params.dictsize = GTrainDICParam.nAtom;
    if GTrainDICParam.nAtom <= 1 && GTrainDICParam.nAtom > 0
        params.dictsize = round(size(new_patches, 2)* GTrainDICParam.nAtom);
    else
        params.dictsize = GTrainDICParam.nAtom;
    end
    params.iternum = GTrainDICParam.iterNum;
    params.memusage = 'high';

    % using the third-party toolbox to train the dictionary
    DIC = ksvd(params, GTrainDICParam.isShowIterInfo);
    
end



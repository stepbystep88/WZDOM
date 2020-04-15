function DIC = bsTrain1DSparseDIC(datas, GTrainDICParam, xs, ys)
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
            patches = [patches, data(j : j+GTrainDICParam.sizeAtom-1)];
        end
    end
   
    
    params.data = patches;
    params.Tdata = GTrainDICParam.sparsity;
    params.dictsize = GTrainDICParam.nAtom;
    params.iternum = GTrainDICParam.iterNum;
    params.memusage = 'high';

    % using the third-party toolbox to train the dictionary
    DIC = ksvd(params, GTrainDICParam.isShowIterInfo);
    
end



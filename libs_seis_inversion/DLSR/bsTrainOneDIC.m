function [DIC, trainIndex] = bsTrainOneDIC(trueModel, trainIndex, GTrainDICParam)
%% This function is used to learn one dictionary from given model
% Programmed by: Bin She (Email: bin.stepbystep@gmail.com)
% Programming dates: Nov 2019
% -------------------------------------------------------------------------

    % create a folder to save the trained dictionary
    if ~exist(GTrainDICParam.dicSavePath,'dir')
       mkdir(GTrainDICParam.dicSavePath);
    end
    
    % get the number of train wells
    if iscell(trueModel)
        trueModel = trueModel(:, trainIndex);
        trainNum = length(trueModel);
    else
        trueModel = trueModel{trainIndex};
        trainNum = size(trueModel, 2);
    end
    
    %% set a file name automatically based on the trianing parameter and stop the training
    % progress if the file already exits
    fileName = sprintf('%s/DIC_trainNum_%d_sizeAtom_%d_nAtom_%d_filtCoef_%.2f.mat', ...
        GTrainDICParam.dicSavePath, ...
        trainNum, ...
        GTrainDICParam.sizeAtom, ...
        GTrainDICParam.nAtom, ...
        GTrainDICParam.filtCoef);
    
    if exist(fileName, 'file')
       load(fileName);
       return;
    end
    
    %% cut the welllog data into pieces and train a dictionary
    GTrainDICParam.iterNum = 5;        

    trainLogs = cell(1, trainNum);
    for i = 1 : trainNum
        if iscell(trueModel)
            trainLogs{i} = trueModel{i}.wellLog(:, 1);
        else
            trainLogs{i} = trueModel(:, i);
        end
    end

    
    % learn dictionary 
    [DIC] = bsLearn1DSparseDIC(trainLogs, GTrainDICParam);
    
    % save the results
    save(fileName, 'DIC', 'trainIndex');
        
end
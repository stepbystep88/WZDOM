function [DIC, trainIndex, rangeCoef] = bsTrainDics(GTrainDICParam, trueModel, trainIndex, attIndex, isReTrain)
%% This function is used to learn one dictionary from given model
% Programmed by: Bin She (Email: bin.stepbystep@gmail.com)
% Programming dates: Nov 2019
% -------------------------------------------------------------------------
    if ~exist('isReTrain', 'var')
        isReTrain = 0;
    end
    
    % create a folder to save the trained dictionary
    if ~exist(GTrainDICParam.dicSavePath, 'dir')
       mkdir(GTrainDICParam.dicSavePath);
    end
    
    % get the number of train wells
    trueModel = trueModel(:, trainIndex);
    trainNum = length(trainIndex);
    
    %% set a file name automatically based on the trianing parameter and stop the training
    % progress if the file already exits
    if isempty(GTrainDICParam.title)
        fileName = sprintf('%s/DIC_%s_trainNum_%d_sizeAtom_%d_nAtom_%d_filtCoef_%.2f.mat', ...
            GTrainDICParam.dicSavePath, ...
            GTrainDICParam.flag, ...
            trainNum, ...
            GTrainDICParam.sizeAtom, ...
            GTrainDICParam.nAtom, ...
            GTrainDICParam.filtCoef);
    else
        fileName = sprintf('%s/DIC_%s_%s__trainNum_%d_sizeAtom_%d_nAtom_%d_filtCoef_%.2f.mat', ...
            GTrainDICParam.dicSavePath, ...
            GTrainDICParam.flag, ...
            GTrainDICParam.title, ...
            trainNum, ...
            GTrainDICParam.sizeAtom, ...
            GTrainDICParam.nAtom, ...
            GTrainDICParam.filtCoef);
    end
    
    if exist(fileName, 'file') && isReTrain==0
       res = load(fileName);
       DIC = res.DIC;
       trainIndex = res.trainIndex;
       rangeCoef = res.rangeCoef;
       
       return;
    end
    
    %% cut the welllog data into pieces and train a dictionary
%     GTrainDICParam.iterNum = 5; 
    nAtt = length(attIndex);
    rangeCoef = [];
    
    switch lower(GTrainDICParam.flag)
        case {'one', 'ssr'}
            DIC = cell(1, nAtt);
            for k = 1 : nAtt
                trainLogs = cell(1, trainNum);
                for i = 1 : trainNum
                    trainLogs{i} = trueModel{i}.wellLog(:, attIndex(k));
                end

                % learn dictionary 
                [DIC{k}] = bsTrain1DSparseDIC(trainLogs, GTrainDICParam);
            end
            
            if strcmpi(GTrainDICParam.flag, 'one')
                % only one dictionary
                DIC = DIC{1};
            end
            
        case 'csr'
            trainLogs = cell(1, trainNum);
            for i = 1 : trainNum
                trainLogs{i} = trueModel{i}.wellLog(:, attIndex);
            end
            [DIC, rangeCoef] = bsTrain1DSparseJointDIC(trainLogs, GTrainDICParam);
            
    end
    
    % save the results
    save(fileName, 'DIC', 'trainIndex', 'rangeCoef');
        
end
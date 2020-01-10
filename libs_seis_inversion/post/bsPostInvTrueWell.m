function [invVals, outputs, model] = bsPostInvTrueWell(GPostInvParam, wellInfo, timeLine, methods)
%% inverse welllog data
% Programmed by: Bin She (Email: bin.stepbystep@gmail.com)
% Programming dates: Nov 2019
% -------------------------------------------------------------------------
%
% Input 
% GPostInvParam     all information of the inverse task
% wellInfo      information and data of target well
% timeLine      horizon information
% methods       the methods to solve the inverse task
% 
% Output
% invVals       inverted results
% outputs       outputs of the iteration process including some intermediate
% results
% -------------------------------------------------------------------------

%     sampNum = GPostInvParam.upNum + GPostInvParam.downNum;
    [~, ~, horizonTime] = bsGetWellBaseInfo(timeLine{GPostInvParam.usedTimeLineId}, ...
        wellInfo.inline, wellInfo.crossline, 1, 2, 1, 2, 3);
    
    % obtain target welllog data based on horizon information
    trueLog = bsExtractWellDataByHorizon(...
                    wellInfo.wellLog, ...
                    horizonTime, ...
                    GPostInvParam.indexInWellData.ip, ...
                    GPostInvParam.indexInWellData.time, ...
                    GPostInvParam.upNum, ...
                    GPostInvParam.downNum, ...
                    1);
                
%     welllog = wellInfo.wellLog;
%     dist = horizonTime - welllog(:, GPostInvParam.indexOfTimeInWellData);
%     [~, index] = min(abs(dist));
%     s = index - GPostInvParam.upNum;
%     trueLog = welllog(s : s+sampNum-1, 1);
    
    % create model data
    model = bsPostPrepareModel(GPostInvParam, wellInfo.inline, wellInfo.crossline, horizonTime, trueLog, []);
    
    %% perform invsersion process
    nMethod = size(methods, 1);
    invVals = cell(1, nMethod);
    outputs = cell(1, nMethod);
    
    for i = 1 : nMethod
        
        method = methods{i};
        
        fprintf('Solving the trace of inline=%d and crossline=%d by using method %s...\n', ...
            wellInfo.inline, wellInfo.crossline, ...
            method.name);
        
%         method.options.optimalX = model.trueX;
        
        tic;
        [xOut, fval, exitFlag, output] = bsPostInv1DTrace(...
            model.d, model.G, model.initX, model.Lb, model.Ub, method);       
 

        invVals{i}.Ip = exp(xOut);
        invVals{i}.model = model;
        invVals{i}.name = method.name;
        output.timeCost = toc;
        output.MRRMSE = bsCalcRRSE(trueLog, model.initLog, invVals{i}.Ip);
        outputs{i} = output;
        
        fprintf('[RRMSE=%.3f, fval=%.3f, timeCost=%.3f, exitFlag=%d]\n', output.MRRMSE, fval, output.timeCost, exitFlag);
    end
    
end
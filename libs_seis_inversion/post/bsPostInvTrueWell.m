function [invVals, model, outputs] = bsPostInvTrueWell(GPostInvParam, wellInfo, timeLine, methods)
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
% model         the data including d, G, m, m_0 of the inverse task
% outputs       outputs of the iteration process including some intermediate
% results
% -------------------------------------------------------------------------

    sampNum = GPostInvParam.upNum + GPostInvParam.downNum;
    [~, ~, horizonTime] = bsCalcWellBaseInfo(timeLine{GPostInvParam.usedTimeLineId}, ...
        wellInfo.inline, wellInfo.crossline, 1, 2, 1, 2, 3);
    
    % obtain true welllog dataq
    welllog = wellInfo.wellLog;
    dist = horizonTime - welllog(:, 2);
    [~, index] = min(abs(dist));
    s = index - GPostInvParam.upNum;
    trueLog = welllog(s : s+sampNum-1, 1);
    
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
        
        tic;
        [xOut, fval, exitFlag, output] = bsPostInv1DTrace(method.flag, ...  % method flag
            model.d, model.G, model.initX, model.Lb, model.Ub, ...      % data
            method.regParam, method.parampkgs, method.options);       
 

        invVals{i} = exp(xOut);
        output.timeCost = toc;
        output.MRRMSE = sqrt(mse(invVals{i} - trueLog));
        outputs{i} = output;
        
        fprintf('[RRMSE=%.2e, fval=%.3e, timeCost=%.2f, exitFlag=%d]\n', output.MRRMSE, fval, output.timeCost, exitFlag);
    end
    
end
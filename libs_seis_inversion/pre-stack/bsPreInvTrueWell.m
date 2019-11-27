function [invVals, model, outputs] = bsPreInvTrueWell(GPreInvParam, wellInfo, timeLine, methods)
%% prestack inversion of welllog data
% Programmed by: Bin She (Email: bin.stepbystep@gmail.com)
% Programming dates: Nov 2019
% -------------------------------------------------------------------------
%
% Input 
% GPreInvParam     all information of the inverse task
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

    sampNum = GPreInvParam.upNum + GPreInvParam.downNum;
    [~, ~, horizonTime] = bsCalcWellBaseInfo(timeLine{GPreInvParam.usedTimeLineId}, ...
        wellInfo.inline, wellInfo.crossline, 1, 2, 1, 2, 3);
    
    % obtain target welllog data
    % obtain target welllog data based on horizon information
    indexInWellData = GPreInvParam.indexInWellData;
    
    trueLog = bsExtractWellDataByHorizon(...
                    wellInfo.wellLog, ...
                    horizonTime, ...
                    [   indexInWellData.depth, ...
                        indexInWellData.vp, ...
                        indexInWellData.vs, ...
                        indexInWellData.rho], ...
                    indexInWellData.time, ...
                    GPreInvParam.upNum, ...
                    GPreInvParam.downNum, ...
                    1);
    
    % create model data
    model = bsPrePrepareModel(GPreInvParam, wellInfo.inline, wellInfo.crossline, horizonTime, trueLog, []);
    
    %% perform invsersion process
    nMethod = size(methods, 1);
    invVals = cell(1, nMethod);
    outputs = cell(1, nMethod);
    
    for i = 1 : nMethod
        
        method = methods{i};
        
        fprintf('Solving the trace of inline=%d and crossline=%d by using method %s...\n', ...
            wellInfo.inline, wellInfo.crossline, ...
            method.name);
        
        method.mode = GPreInvParam.mode;
        method.lsdCoef = model.lsdCoef;
        
        tic;
        [xOut, fval, exitFlag, output] = bsPreInv1DTrace(...
            model.d, model.G, model.initX, model.Lb, model.Ub, method);       
 
        [res.vp, res.vs, res.rho] = bsPreRecoverElasticParam(xOut, GPreInvParam.mode, model.lsdCoef);
        res.t0 = model.t0;
        res.name = method.name;
        res.model = model;
        res.inline = wellInfo.inline;
        res.crossline = wellInfo.crossline;
        
        invVals{i} = res;
        
        output.timeCost = toc;
        output.RMSEVp = sqrt(mse(res.vp - trueLog(:, 2)));
        output.RMSEVs = sqrt(mse(res.vs - trueLog(:, 3)));
        output.RMSERho = sqrt(mse(res.rho - trueLog(:, 4)));
        outputs{i} = output;
        
        fprintf('[RMSEVp=%.2e, RMSEVs=%.2e, RMSERho=%.2e, fval=%.3e, timeCost=%.2f, exitFlag=%d]\n', ...
            output.RMSEVp, output.RMSEVs, output.RMSERho, fval, output.timeCost, exitFlag);
    end
    
end
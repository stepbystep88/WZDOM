function bestLambda = bsBestParameterByBisection(trueModel, logLeft, logRight, nIter, inputObjFcnPkgs, xInit, Lb, Ub, options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate the crossplot of objective function value and regularization term for different lambdas
%
% Programmed by: Bin She (Email: bin.stepbystep@gmail.com)
% Programming dates: May 2019
% -------------------------------------------------------------------------
% INPUT
%
% trueModel: true model, must be given
%
% logLeft: logarithm of the minimum regularization pramameter
%
% logRight: logarithm of the maximum regularization pramameter
% 
% nIter: the number of iteration for bisection 
%
% inputObjFcnPkgs       is a function handel, or a cell made up by a seris of
% functions, the data for calculating each function, and the weight
% coefficients. It could be like the following three forms:
% 1. function (only one function, no need for extra data)
% 2. {function, struct} (only one function, need extra data)
% 3. {function1, struct1, weight1; function2 struct2; weight2; ...} 
% (mutiple functions)
% 
% xInit : the initial guess of the parameter x.
% 
% Lb: lower boundary
%
% Ub: upper boundary
%
% options     is a struct which contains the data and infomation of the way
% of calculating the objective function. 
%
% -------------------------------------------------------------------------
% OUTPUT
%
% bestLambda: the best lambda choosen by L-curve criteria
%
% -------------------------------------------------------------------------

    mseOld = zeros(1, 3);
    
    [mseOld(1), xOut] = bsGetMse(logLeft, inputObjFcnPkgs, xInit, Lb, Ub, options, trueModel);
    [mseOld(3), xOut] = bsGetMse(logRight, inputObjFcnPkgs, xInit, Lb, Ub, options, trueModel);
    
    while nIter > 0
        
        logMid = 0.5 * (logLeft + logRight);
        
        [mseOld(2), xOut] = bsGetMse(logMid, inputObjFcnPkgs, xInit, Lb, Ub, options, trueModel);
        
        [minMse, minIndex] = min(mseOld);
        
        switch minIndex
            case 1
                bestLambda = logLeft;
                logRight = logMid;
                mseOld(3) = mseOld(2);
            case 2
                bestLambda = logMid;
                break;
            case 3
                bestLambda = logRight;
                logLeft = logMid;
                mseOld(1) = mseOld(2);
        end
        nIter = nIter - 1;
    end
    
    for i = 1 : nIter
        logLeftNew = 0.5 * (logLeft + logMid);
        [mseLeft, xOut] = bsGetMse(logLeftNew, inputObjFcnPkgs, xInit, Lb, Ub, options, trueModel);
        
        if minMse > mseLeft
            logRight = logMid;
            logMid = logLeftNew;
            minMse = mseLeft;
            bestLambda = logLeftNew;
            continue;
        end
        
        logRightNew = 0.5 * (logRight + logMid);
        [mseRight, xOut] = bsGetMse(logRightNew, inputObjFcnPkgs, xInit, Lb, Ub, options, trueModel);
        
        if minMse > mseRight
            logLeft = logMid;
            logMid = logRightNew;
            minMse = mseRight;
            bestLambda = logRightNew;
            continue;
        end
        
        logLeft = logLeftNew;
        logRight = logRightNew;
        
    end
    
    bestLambda = exp(bestLambda);
    
end

function [mseModel, xOut] = bsGetMse(logLambda, inputObjFcnPkgs, xInit, Lb, Ub, options, trueModel)
    
    if ~isempty(options.plotFcn)
        fprintf('In bisection method. Runing the optimization process with regularization parameter=%d\n', exp(logLambda));
    end
        
    % update the lambda
    inputObjFcnPkgs{2, 3} = exp(logLambda);
        
    % call bsGBSolverByOptions to solve the problem at current lambda
    xOut = bsGBSolveByOptions(inputObjFcnPkgs, xInit, Lb, Ub, options);
%     mseModel(1) = sqrt(mse(xOut - trueModel));
    mseModel(1) = sum(abs(xOut - trueModel)) + 2*sum(abs(gradient(xOut)));
end
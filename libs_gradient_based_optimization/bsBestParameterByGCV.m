function [bestLambda, curveData] = bsBestParameterByGCV(lambdas, inputObjFcnPkgs, xInit, Lb, Ub, options, isShowFigure)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate the best regularization parameter by GCV
%
% Programmed by: Bin She (Email: bin.stepbystep@gmail.com)
% Programming dates: October 2019
% -------------------------------------------------------------------------
% INPUT
%
% lambdas       is a array; refers to a seris of regularization parameters.
% 
% inputObjFcnPkgs       is a function handel, or a cell made up by a seris of
% functions, the data for calculating each function, and the weight
% coefficients. It could be like the following three forms:
% 1. function (only one function, no need for extra data)
% 2. {function, struct} (only one function, need extra data)
% 3. {function1, struct1, weight1; function2 struct2; weight2; ...} 
% (mutiple functions)
%
% xInit         is a colmun vector; refers to the initial guess of the
% parameter x.
% 
% Lb: lower boundary
%
% Ub: upper boundary
%
% options     is a struct which contains the data and infomation of the way
% of calculating the objective function. 
% 
% isShowFigure is a logical; When it is true, the L-Curve will be
% displayed.
%
% -------------------------------------------------------------------------
% OUTPUT
%
% bestLambda: the best lambda choosen by GCV criteri
% 
% curveData is a array of size n*2, the first column corresponds to the
% data residual error, the second column is the regularization term.
%
% -------------------------------------------------------------------------

    nLambdas = length(lambdas);
    curveData = zeros(nLambdas, 3);
    mainData = inputObjFcnPkgs{1, 2};
    
    s = svd(mainData.A);
    m = size(mainData.A, 1);
    
    for i = 1 : nLambdas
        
        if ~isempty(options.plotFcn)
            fprintf('In linear search method. Runing the optimization process with regularization parameter=%d\n', lambdas(i));
        end
        
        % update the lambda
        inputObjFcnPkgs{2, 3} = lambdas(i);
        
        % call bsGBSolverByOptions to solve the problem at current lambda
        [xOut, fval, exitFlag, output] = bsGBSolveByOptions(inputObjFcnPkgs, xInit, Lb, Ub, options);
        
        % update the inputObjFckPkgs, this is due to we initialize some
        % data of the inputObjFckPkgs through GBSolver
        inputObjFcnPkgs = output.inputObjFcnPkgs;
        
        % save the residual error and regularization term.
        curveData(i, 1) = inputObjFcnPkgs{1, 1}(xOut, mainData);
        curveData(i, 2) = inputObjFcnPkgs{2, 1}(xOut, inputObjFcnPkgs{2, 2});
        
        residual = (mainData.A * xOut - mainData.B);
        Nume = residual' * residual;
        
        tmp1 = sum(s.^2 ./ (s.^2 + lambdas(i)));
        Deno = (m - tmp1).^2;
        
        curveData(i, 3) = Nume / Deno;
    end
    
    
    
    
    [~, index] = min(curveData(:, 3));
    bestLambda = lambdas(index);
    
    
    if isShowFigure == 1
        
        plotSet = bsGetDefaultPlotSet();
        
        figure;
        plot(log(lambdas), curveData(:, 1), 'k-*', 'linewidth', plotSet.linewidth);
        xlabel('\lambda (log scale)');
        ylabel('Objective function');
        bsPlotSetDefault(bsGetDefaultPlotSet());
        
        figure;
        plot(log(lambdas), curveData(:, 3), 'k-*', 'linewidth', plotSet.linewidth);
        xlabel('\lambda (log scale)');
        ylabel('GCV(\lambda)');
        
        bsPlotSetDefault(bsGetDefaultPlotSet());
        
    end
    
end
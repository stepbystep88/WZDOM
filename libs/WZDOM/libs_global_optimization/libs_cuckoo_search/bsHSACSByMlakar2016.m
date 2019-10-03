function [x,fval,exitFlag,output] = bsHSACSByMlakar2016(objFunc, Lb, Ub, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% This code modifies the Modified Cuckoo Search (MCS) 
% algorithm implemented by Walton
%
% Programmed by Bin She (bin.stepbystep@gmail.com)
% Programmed dates: May 2019
%
% -------------------------------------------------------------------------
% The work is based on the following paper:
%
% 1) %S.Walton, O.Hassan, K.Morgan and M.R.Brown "Modified cuckoo search: A
% new gradient free optimisation algorithm" Chaos, Solitons & Fractals Vol
% 44 Issue 9, Sept 2011 pp. 710-718 DOI:10.1016/j.chaos.2011.06.004
% -------------------------------------------------------------------------

    % parse some basic parameters for this process
    p = inputParser;
%     rng(125789);

    p = bsAddBasicalParameters(p, length(Lb));
    
%     addParameter(p, 'nNest', 25 );      % The number of nests
    addParameter(p, 'minNests', 10 );      
    addParameter(p, 'maxNests', 100 );        
    addParameter(p, 'pa', 0.1 ); 
    addParameter(p, 'pb', 0.8 ); 
    addParameter(p, 'pe', 0.1 ); 
    addParameter(p, 'pd', 0.9 ); 
    addParameter(p, 'cr', 0.9 ); 
    addParameter(p, 'alpha1', 0.9 ); 
    addParameter(p, 'lambda', 1.5);
    addParameter(p, 'pbest', 4 ); 
    addParameter(p, 'Fi', 0.5);
    addParameter(p, 'nHistory', 20);
    
    % set default function of generation initial population as @bsGenerateInitialPopulationByLHS
    addParameter(p, 'initialPopulationFcn', @bsGenerateInitialPopulationByRandom );
    
    addParameter(p, 'interval', 5 );    % for each interval iterations, perform gradient descent algorithm once
    
    addParameter(p, 'innerMaxIter', 100 ); % the max number iterations for gradient descent algorithm
    addParameter(p, 'optionsForGB', [] );
    
    p.parse(varargin{:});  
    params = p.Results;
    
    nNest = params.maxNests;
    
    
    % call initial function to generate the initial population
    nests = params.initialPopulationFcn(Lb, Ub, nNest);
    
    % Get the current bestNest
    fitness = inf * ones(nNest, 1);
    nDim = length(Lb);
    gradient = zeros(nDim, nNest);
    
%     gradientObjFunc = @(x)(objFunc(x, 0));
    gradientObjFunc = @(x)(objFunc(x, 1));
%     [nests, fitness, gradient] = bsGetFitnessWithGradient(gradientObjFunc, nests, nests, fitness, gradient);
    [globalMinFVal, globalBestNest, nests, fitness, gradient, K] = bsGetBestNestWithGradient(gradientObjFunc, nests, nests, fitness, gradient);
%     xInit = globalBestNest;
    
    nfev = nNest;   %count the number of evaluations of objective functions
    fs = [];
    
    paramsConfigure = {
        params.pa, 0.05, 0.2, 'pa';
        params.pb, 0.5, 1.0, 'pb';
        params.pe, 0.05, 0.2, 'pe';
        params.pd, 0, 1, 'pd';
        params.cr, 0, 1, 'cr';
        params.alpha1, 0, 1, 'alpha1';
        params.lambda, 1.2, 1.8, 'lambda';
        params.Fi, 0.05, 1, 'Fi';
    };

%     nParameter = size(paramsConfigure, 1);
    algParams = bsCreateHistoryStruct(paramsConfigure, params.nHistory);

    %% Starting iterations
    for iter = 1 : params.maxIter
        
        % calculate the new parameters based on the history information
        algParams = bsCalNewValueOfAlgParamsSpecail(algParams);
        
        % sort the current nests in order of fitness
        nNest = bsCalcNNest(params.minNests, params.maxNests, params.maxFunctionEvaluations, nfev);
        [nests, gradient, fitness] = bsSortNests(nests, gradient, fitness, nNest);
        
        if iter == 1
            % get the current best nest
            globalMinFVal = fitness(1);
            globalBestNest = nests(:, 1);
            xInit = globalBestNest;
            fs = [0; 0; globalMinFVal];
        end
        
        
        % Generate new solutions (but keep the current bestNest)
        newNest = bsGetCuckoosByTwoStrategies(nests, algParams, params.pbest, Lb, Ub);
        [nests, newFitness, gradient] = bsGetFitnessWithGradient(gradientObjFunc, nests, newNest, fitness, gradient);
        algParams = bsCalcScoreOfParametersSpecail(algParams, fitness, newFitness, 1, params.nHistory);
        fitness = newFitness;

        % Update the counter
        nfev = nfev + nNest; 
       
        % Discovery and randomization
        newNest = bsEmptyNests(nests, Lb, Ub, algParams.pa.value) ;
        
        [nests, newFitness, gradient] = bsGetFitnessWithGradient(gradientObjFunc, nests, newNest, fitness, gradient);
        algParams = bsCalcScoreOfParametersSpecail(algParams, fitness, newFitness, 2, params.nHistory);
        fitness = newFitness;
        
        % Update the counter
        nfev = nfev + nNest; 
        [~, index] = min(fitness);
        bestFVal = fitness(index);
        bestNest = nests(:, index);
        
        
        % Find the bestNest objective so far  
        if bestFVal < globalMinFVal
            globalMinFVal = bestFVal;
            globalBestNest = bestNest;
        end
        
        % Update the counter again
        if params.isSaveMiddleRes
            fs = [fs, [iter; nfev; globalMinFVal]];
        end
        
        data.fNew = globalMinFVal;
        data.nfev = nfev;
        data.iter = iter;
        
        exitFlag = bsCheckStopCriteria(data, params);
        [data] = bsPlotMidResults(xInit, data, params, Lb, Ub, exitFlag > 0);
        
        if exitFlag > 0
            break;
        end
        
%         data.fOld = globalMinFVal;
%         data.xOld = globalBestNest;
        
    end %% End of iterations

    x = globalBestNest;
    fval = globalMinFVal;
    output.funcCount = nfev;
    output.iterations = iter;
    output.midResults = fs;  
    output.frames = data.frames;
end

function algParams = bsCalNewValueOfAlgParamsSpecail(algParams)
    
    
    names = fieldnames(algParams);
    nField = length(names);
    
    for i = 1 : nField
        t = getfield(algParams, names{i});
        
        nTrain = t.nTrain;
        
        if nTrain <= 1
            continue;
        end
          
%         weight = t.store(1:nTrain, 2) ./ sum(t.store(1:nTrain, 2));
%         S = t.store(1:nTrain, 1);
%         t.value = sum(weight .* S .^ 2) / sum(weight .* S);
    
        t.value = t.min + randi(floor((t.max-t.min)*10000))/10000;
        if rand() < 0.4
            algParams = setfield(algParams, names{i}, t);
        end
            
%         algParams = setfield(algParams, names{i}, t);
    end
end

function algParams = bsCalcScoreOfParametersSpecail(algParams, fitness, newFitNess, flag, nHistory)

    delatFit = fitness - newFitNess;
    increasedFit = sum(delatFit(delatFit>0));
%     increasedFit = sum(delatFit>0);

    if increasedFit > 0
        switch flag
            case 1
                algParams.alpha1 = bsUpdateParameter(algParams.alpha1, increasedFit, nHistory);
                algParams.Fi = bsUpdateParameter(algParams.Fi, increasedFit, nHistory);
%                 algParams.pe = bsUpdateParameter(algParams.pe, increasedFit, nHistory);
%                 algParams.pb = bsUpdateParameter(algParams.pb, increasedFit, nHistory);
                algParams.cr = bsUpdateParameter(algParams.cr, increasedFit, nHistory);
                algParams.lambda = bsUpdateParameter(algParams.lambda, increasedFit, nHistory);
                
%                 algParams.alpha.history(iHistory, 2) = increasedFit;
%                 algParams.Fi.history(iHistory, 2) = increasedFit;
%                 algParams.pe.history(iHistory, 2) = increasedFit;
%                 algParams.pb.history(iHistory, 2) = increasedFit;
%                 algParams.cr.history(iHistory, 2) = increasedFit;
%                 algParams.lambda.history(iHistory, 2) = increasedFit;
            case 2
%                 algParams.alpha.pa(iHistory, 2) = increasedFit;
%                 algParams.pa = bsUpdateParameter(algParams.pa, increasedFit, nHistory);
%                 algParams.pd = bsUpdateParameter(algParams.pd, increasedFit, nHistory);
        end
    end
    
end

function t = bsUpdateParameter(t, increasedFit, nHistory)
    t.index = t.index + 1;
    t.nTrain = t.nTrain + 1;
    
    if (t.index > nHistory)
        t.index = 1;
    end
    
    if (t.nTrain > nHistory)
        t.nTrain = nHistory;
    end
    
    t.store(t.index, 2) = increasedFit;
    t.store(t.index, 1) = t.value;
end

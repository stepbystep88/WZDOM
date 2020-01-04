function [x,fval,exitFlag,output] = bsHSACSByShe2019(objFunc, Lb, Ub, varargin)
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
    addParameter(p, 'nHistory', 20);
    
    % set default function of generation initial population as @bsGenerateInitialPopulationByLHS
    addParameter(p, 'initialPopulationFcn', @bsGenerateInitialPopulationByRandom );
    
    p.parse(varargin{:});  
    params = p.Results;
    
    nNest = params.maxNests;
    
    
    % call initial function to generate the initial population
    nests = params.initialPopulationFcn(Lb, Ub, nNest);
    
    % Get the current bestNest
    fitness = inf * ones(nNest, 1);
    nDim = length(Lb);
    
    gradientObjFunc = @(x)(objFunc(x, 0));
%     gradientObjFunc = @(x)(objFunc(x, 1));
%     [nests, fitness, gradient] = bsGe(gradientObjFunc, nests, nests, fitness, gradient);
    [~, ~, nests, fitness, ~] = bsGetBestNest(gradientObjFunc, nests, nests, fitness);
%     [globalMinFVal, globalBestNest, nests, fitness, gradient, K] = bsGetBestNestWithGradient(gradientObjFunc, nests, nests, fitness, gradient);
%     xInit = globalBestNest;
    
    nfev = nNest;   %count the number of evaluations of objective functions
    fs = [];

    paramsConfigure = {
        
        0.8, 0.66, 1, 'pb', 1, 0.1;
        0.9, 0.0, 1, 'pe', 1, 0.1;
        1.5, 1.1, 1.9, 'lambda', 1, 0.1;
        0.9, 0.05, 2, 'alpha1', 1, 0.3;
        0.9, 0.05, 2, 'alpha3', 1, 0.3;
        0.9, 0, 1, 'cr', 1, 0.1;
        0.1, 0, 1, 'pf', 1, 0.1;
        0.9, 0.05, 2, 'alpha4', 1, 0.3;
        
        0.9, 0, 1, 'pa', 2, 0.1;
        0.8, 0, 1, 'pg', 2, 0.1;
        0.9, 0.05, 2, 'alpha2', 2, 0.3;
    };


    pbIndex = 1;
    peIndex = 2;
    lambdaIndex = 3;
    alpha1Index = 4; 
    alpha3Index = 5;
    
    crIndex = 6;
    pfIndex = 7;
%     

    paIndex = 1;
    alpha2Index = 3;

    nHistory = params.nHistory;
%     nParameter = size(paramsConfigure, 1);
    algParams = bsCreateHistoryStructNew(paramsConfigure, nHistory, nNest);
    nField = size(paramsConfigure, 1);
    paramHistory = zeros(nField, params.maxIter);
    maxFEV = params.maxFunctionEvaluations;
%     pbest = params.pbest;

    nNoUpDate = 0;

    %% Starting iterations
    for iter = 1 : params.maxIter
        

        nNest = bsCalcNNest(params.minNests, params.maxNests, maxFEV, nfev);

        % calculate the new parameters based on the history information
        
%         paramHistory = bsSaveParamInformation(paramHistory, algParams, iter, params.maxIter);
% 
        [algParams, nNoUpDate] = bsCalNewValueOfAlgParamsNew(algParams, nNest, nNoUpDate, 0, iter);

        st = 1;
        for i = 1 : 2
            t = algParams{i};
            
            paramHistory(st:st+t.nparam-1, iter) = mean(t.values, 2);
            st = st + t.nparam;
        end
        
        
        
        [nests, fitness] = bsSortNestsWithoutGradient(nests, fitness, nNest);
        
        if iter == 1
            % get the current best nest
            globalMinFVal = fitness(1);
            globalBestNest = nests(:, 1);
            xInit = globalBestNest;
            fs = [0; 0; globalMinFVal];
        end
        
        
        % Generate new solutions (but keep the current bestNest)
%         newNest = bsGetCuckoosByTwoStrategiesNew(nests, algParams, params.pbest, Lb, Ub);
        d1 = algParams{1}.values;
        
%         pbest, alphas, c1s, c2s, lambda, pbs, pes, crs

        c1s = randRepmat(d1(alpha3Index, :), 0, 1, 0.1, nDim); 
        c2s = randRepmat(d1(alpha3Index, :), 0, 1, 0.1, nDim);
        alpha1s = randRepmat(d1(alpha1Index, :), 0, 1, 0.1, nDim);
        newNest = bsGetCuckoosByTwoStrategiesNew(nests, Lb, Ub, ...
            alpha1s, ...
            c1s, ...
            c2s, ...
            d1(lambdaIndex, :), ...
            d1(pbIndex, :), ...
            d1(peIndex, :), ...
            d1(crIndex, :), ...
            d1(pfIndex, :)...
        );
        [~, ~, nests, newFitness, ~] = bsGetBestNest(gradientObjFunc, nests, newNest, fitness);
        algParams = bsCalcScoreOfParametersNew(algParams, fitness, newFitness, 1, nHistory);
        fitness = newFitness;

        % Update the counter
        nfev = nfev + nNest; 
       
        d2 = algParams{2}.values;

        a2pha1s = randRepmat(d2(alpha2Index, :), 0, 1, 0.3, nDim);
        newNest = bsEmptyNestsTwoStrategiesNew(nests, Lb, Ub, ...
            d2(paIndex, :), ...
            a2pha1s);
        [bestFVal, bestNest, nests, newFitness, ~] = bsGetBestNest(gradientObjFunc, nests, newNest, fitness);
        algParams = bsCalcScoreOfParametersNew(algParams, fitness, newFitness, 2, nHistory);
        fitness = newFitness;
        
        % Update the counter
        nfev = nfev + nNest; 

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

    paramHistory(:, iter+1:end) = [];
    x = globalBestNest;
    fval = globalMinFVal;
    output.funcCount = nfev;
    output.iterations = iter;
    output.midResults = fs;  
    output.frames = data.frames;
    output.paramHistory = paramHistory;
end

function [nests, fitness] = bsSortNestsWithoutGradient(nests, fitness, nNewNest)

    [nDim, nNest] = size(nests);

    piFi = [nests', fitness];
    % Sort by Fi in assending order
    piFiS = sortrows(piFi, nDim+1);

    % Decrease number of nests, we only need lots of nests initially to get
    % a good sampling of the objective function
    nests = piFiS(1:nNewNest, 1:nDim)';
    fitness = piFiS(1:nNewNest, end);
        
end


function params = randRepmat(params, minVal, maxVal, coef, nDim)   
    nNest = length(params);
    params = repmat(params, nDim, 1);

    params = params - coef + 2*coef*rand(nDim, nNest);

    params(params>maxVal) = maxVal;
    params(params<minVal) = minVal;
end

function [parameters] = bsCreateHistoryStructNew(paramsConfigure, nHistory, nNest)
%     nParameter = size(paramsConfigure, 1);
    
    t = cell(1, 2);
    flags = cell2mat(paramsConfigure(:, 5));
    seqs = 1 : size(paramsConfigure, 1);
    for i = 1 : 2
        t{i}.index = 0;
        t{i}.nTrain = 0;
        index = flags == i;
        n = sum(index);
        t{i}.nparam = n;
        t{i}.names = paramsConfigure(index, 4);
        t{i}.store = zeros(nHistory, n+1);
        
        t{i}.values = zeros(n, nNest);
        t{i}.max = zeros(n, 1);
        t{i}.min = zeros(n, 1);
        t{i}.coef = zeros(n, 1);
    
        seq = seqs(index);

        for j = 1 : n
            k = seq(j);
            
            minVal = paramsConfigure{k, 2};
            maxVal = paramsConfigure{k, 3};
            coef = paramsConfigure{k, 6};

            t{i}.min(j) = minVal;
            t{i}.max(j) = maxVal;
            t{i}.coef(j) = coef;

%             t{i}.values(j, :) = bsGenValsByMinMax(minVal, maxVal, 1, nNest);
%             t{i}.values(j, :) = paramsConfigure{k, 1} - coef + coef*2*rand(1, nNest);
        end

%         t{i}.values = bsGenerateInitialPopulationByLHS(t{i}.min, t{i}.max, nNest);
        t{i}.values = bsGenValsByMinMax(t{i}.min, t{i}.max, nNest, t{i}.coef);
    end
    
    parameters = t;

end

function valSeq = bsGenValsByMinMax(minVal, maxVal, n, coef)
    m = length(minVal);

    value = minVal + (maxVal-minVal).*rand(m, 1);
    value = repmat(value, 1, n);
    coef = repmat(coef, 1, n);
    valSeq = value - coef + 2*coef.*rand(m, n);

%     minVals = repmat(minVal, 1, n);
%     maxVals = repmat(maxVal, 1, n);
%     valSeq = minVals + (maxVals - minVals) .* rand(m, n);

    valSeq = bsSimpleBounds(valSeq, minVal, maxVal);
end

function algParams = bsCalcScoreOfParametersNew(algParams, fitness, newFitNess, flag, nHistory)

    deltaFit = fitness - newFitNess;
    deltaFit(deltaFit<0) = 0;

%     increasedFit = sum(deltaFit(deltaFit>0)) / length(fitness);
%     increasedFit = sum(delatFit>0);
    increasedFit = deltaFit;
    
    algParams{flag} = bsUpdateParameterNew(algParams{flag}, increasedFit, nHistory);
    
end

function t = bsUpdateParameterNew(t, increasedFit, nHistory)

%     I = increasedFit<=0;
%     increasedFit(I) = [];
%     I = increasedFit>=0;

    nNest = length(increasedFit);
    
    spos = t.index + 1;
    t.index = t.index + nNest;
    t.nTrain = t.nTrain + nNest;
    
%     values = mean(t.values, 2);
% %     values = values';
    values = t.values';

    if (t.index <= nHistory)
        t.store(spos:t.index, end) = increasedFit;
        t.store(spos:t.index, 1:t.nparam) = values;
    else
        n1 = nHistory - spos + 1;
        n2 = nNest - n1;
%         value = t.values';

        t.store(spos:nHistory, end) = increasedFit(1:n1);
        t.store(spos:nHistory, 1:t.nparam) = values(1:n1, :);
        t.store(1:n2, end) = increasedFit(n1+1:end);
        t.store(1:n2, 1:t.nparam) = values(n1+1:end, :);
        
        t.index = n2;
    end
    
    if (t.nTrain > nHistory)
        t.nTrain = nHistory;
    end
    
end

function [algParams, nNoUpDate] = bsCalNewValueOfAlgParamsNew(algParams, nNest, nNoUpDate, isRest, iter)
    
    isNoUpdate = false;

    for i = 1 : length(algParams)
        t = algParams{i};
        

        if isRest
            t.nTrain = 0;
            t.index = 0;
        end

        nTrain = t.nTrain;
        [nHistory, nParam] = size(t.store);
        nParam = nParam - 1;

        if nTrain >= size(t.store, 1)
            
            
            sortData = sortrows(t.store(1:nTrain, :), -(nParam+1));
            sortData = sortData(1:5, :);

            sumWeight = sum(sortData(:, nParam+1));

            if sumWeight > 0
%                 weight = sortData(:, nParam+1) ./ sumWeight; 
                S = sortData(:, 1:nParam);

%                 value = S' * weight;
                
%                 [~, index] = max(t.store(:, 2));
    %             value = t.store(index, 1);

%                 trueVal = normrnd(value, 0.1, 1, nNest);
                value = mean(S, 1)';
                value = repmat(value, 1, nNest);
                coef = repmat(t.coef, 1, nNest);
                trueVal = value - coef + 2*coef.*rand(nParam, nNest);

                trueVal = bsSimpleBounds(trueVal, t.min, t.max);

                t.values = trueVal;
            else
                
%                 value = mod(iter, 100)/100;
%                 value = value + randi([-0.1, 0.1]*1e6, 1, nNest)/1e6 + 0.05;
%                 t.value = bsSimpleBounds(value, t.min, t.max);
%                 t.values = t.min + randi(floor(1e6*(t.max-t.min)), 1, nNest)/1e6;
%                 t.values = bsGenerateInitialPopulationByLHS(t.min, t.max, nNest);
                t.values = bsGenValsByMinMax(t.min, t.max, nNest, t.coef);
                isNoUpdate = true;
            end
            
        else
%             t.value = t.min + randi(floor(1e6*(t.max-t.min)), 1, nNest)/1e6;
%             t.values = bsGenerateInitialPopulationByLHS(t.min, t.max, nNest);
            t.values = bsGenValsByMinMax(t.min, t.max, nNest, t.coef);
        end
        
        
        algParams{i} = t;

        
    end

    
        
    if isNoUpdate
        nNoUpDate = nNoUpDate + 1;
    else
        nNoUpDate = 0;
    end
        

    
end

function newNests = bsGetCuckoosByTwoStrategiesNew(nests, Lb, Ub, alphas, c1s, c2s, lambda, pbs, pes, crs, pfs)

    % Levy flights
    [nDim, nNest] = size(nests);

    %% Levy flights by Mantegna'nest algorithm
%     steps = bsLevy(size(nests, 1), nNest, beta);
    crs = repmat(crs, nDim, 1);
%     c1s = repmat(c1s, nDim, 1);
    
    newNests = nests;
    
%     pbestNest = nests(:, 1);
    pbest = floor(pes * nNest);
    pbest(pbest < 1) = 1;

    methods = rand(1, nNest) < pbs;
    steps = bsLevy(nDim, sum(methods), lambda(methods));
    it = 1;

    for i = 1 : nNest
        nest = nests(:, i);

%         choose method
        if methods(i)
    
            pbestNest = nests(:, randi(pbest(i)));

            if rand() < pfs(i)
                stepsize = alphas(i) .* steps(:, it) .* (nest - pbestNest);
                it = it + 1;
            else
                stepsize = alphas(:, i) .* steps(:, it);
                it = it + 1;
            end
            
        else
            % randomly choose one pbest from nests (nests has been ordered by fitness)
            pbestNest = nests(:, randi(pbest(i)));
            
            r1 = randi(nNest);
            r2 = randi(nNest);
            
            stepsize = c1s(:, i) .* (pbestNest - nest + nests(:, r1) - nests(:, r2));
        end
        
        
        newNests(:, i) = nest + stepsize;
    end
    
    % boundary 
    newNests = bsBetterBounds(newNests, Lb, Ub, nests);
%     newNests = bsSimpleBounds(newNests, Lb, Ub);
    
%     random copy the original nests
    K1 = rand(size(nests)) <= crs;
%     K2 = repmat((1:nDim)', 1, nNest) == randi(nDim, nDim, nNest);
%     K = K1 | K2;
    K = K1;
%     
    newNests = newNests .* K + nests .* (~K);
    
end

function newNests = bsEmptyNestsTwoStrategiesNew(nests, Lb, Ub, pa, alpha2)

    % A fraction of worse nests are discovered with a probability pa
    [nDim, nNest] = size(nests);
    
    pa = repmat(pa, nDim, 1);
%     alpha2 = repmat(alpha2, nDim, 1);
%     alpha2 = alpha2 - 0.1 + 0.2*rand(nDim, nNest);
%     alpha2 = rand(nDim, nNest);

    K = rand(size(nests)) <= pa;

    stepsize = alpha2 .* (nests(:, randperm(nNest)) - nests(:, randperm(nNest)));
    newNests = nests + stepsize .* K;
    % call bsSimpleBounds to project the new nests into the reasonable range
%     newNests = bsBetterBounds(newNests, Lb, Ub, nests);
    newNests = bsSimpleBounds(newNests, Lb, Ub);

end
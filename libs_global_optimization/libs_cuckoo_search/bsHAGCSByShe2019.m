function [x,fval,exitFlag,output] = bsHAGCSByShe2019(objFunc, Lb, Ub, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Programmed by Bin She (bin.stepbystep@gmail.com)
% Programmed dates: May 2019
%
% -------------------------------------------------------------------------

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                           prepare data and parameters                              %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % parse some basic parameters for this process
    p = inputParser;
    % set some basic parameters to control the iteration process, for
    % example, to control the stop criteria incluing the maximum number of
    % function evaluations, and the maximum number of iterations
    p = bsAddBasicalParameters(p, length(Lb));
    
    % the mimimum number of nests exists in the population
    addParameter(p, 'minNests', 5 );       
    % the maximum number of nests exists in the population
    addParameter(p, 'maxNests', 100 );      
    % the number of history iterations saved for generating new control parameters
    addParameter(p, 'nHistory', 500);        
    
    % set default function of generation initial population as @bsGenerateInitialPopulationByLHS
    addParameter(p, 'initialPopulationFcn', @bsGenerateInitialPopulationByRandom );
    % for every [interval] iterations, perform gradient descent algorithm once
    addParameter(p, 'interval', 50 );    
    
    % the max number iterations for gradient-based algorithm to perform local search 
    addParameter(p, 'innerMaxIter', 100 ); 
    % the options parameters for the gradient-based algorithm, we use
    % GBSolver to solve the optimization problems
    addParameter(p, 'GBOptions', [] );
    
    % whether to save the detail update information of all population. The
    % value is set to yes when we need to display an animation of
    % optimization process in general.
    addParameter(p, 'isSaveDetailUpdates', false);
    
    p.parse(varargin{:});  
    params = p.Results;
    
    % the current number of nests is set to the maximum number of nest, it
    % will decrease to params.minNests as iteration goes on
    nNest = params.maxNests;
    
    % function handle of objective function, 0 means we don't need to
    % calculate the gradient information, see bsRosenbrock.m for example
    noGradientFunc = objFunc;
    % call initial function to generate the initial population
    if nargin(params.initialPopulationFcn) == 3
        nests = params.initialPopulationFcn(Lb, Ub, nNest);
    else
        % in this case, the objective function will be used to create the initial population 
        nests = params.initialPopulationFcn(Lb, Ub, nNest, noGradientFunc);
    end
    
    % create the optiions struct for the GBSolver, see bsGBSolver.m
    if isempty(params.GBOptions)
        GBOptions = bsCreateGBOptions(length(Lb), ...
            'optimalFunctionTolerance', params.optimalFunctionTolerance, ...
            'maxIter', params.innerMaxIter, ...
            'isSaveMiddleRes', params.isSaveMiddleRes, ... % whether to save the middle results during the iterations
            'display', 'off', ...
            'optimalF', params.optimalF);                  % optimalF is the groundtruth of the global minimum function value
    else
        [~, GBOptions] = bsGetFieldsWithDefaults(params.GBOptions, {...
            'optimalFunctionTolerance', params.optimalFunctionTolerance; ...
            'maxIter', params.innerMaxIter; ...
            'isSaveMiddleRes', params.isSaveMiddleRes; ... % whether to save the middle results during the iterations
            'display', 'off'; ...
            'optimalF', params.optimalF ...
        });
    end
    
    %% ---------------------------------------------------------------------------------------------------------------%%
    % initialzie a population
    % ----------------------------------------------------------------------------------------------------------------%%
    nDim = length(Lb);                  % the dimension of the problem to be solved
    % fitness records the score of each nest
    % in this case, the smaller its value is, the better the corresponding solution is.    
    fitness = inf * ones(nNest, 1);   
    gradient = zeros(nDim, nNest);      % used to save the gradient information of each nest
    
    % objective function, this function will calcuate the gradient information
    gradientObjFunc = objFunc;
    % calculate the gradient information, and fitness of the initilized population. 
    [nests, fitness, gradient] = bsGetFitnessWithGradient(gradientObjFunc, nests, nests, fitness, gradient);
    nfev = nNest;   % count the number of evaluations of objective functions

    %% ---------------------------------------------------------------------------------------------------------------%%
    % configure the contorl parameters, see the Table 1 of our paper
    % ----------------------------------------------------------------------------------------------------------------%%
    paramsConfigure = {
        % each row represents min, max, name, group ID, variation
        % here, if the gropu ID = 1, the parameters will be used only in
        % the first for loop, see Algorithm 2 of our papaer. Otherwise, the
        % parameters will be only used in the second for loop.
        
        0.0, 1, 'pe', 1, 0.1;           % Elitist parameter determining the greediness of the elitist selection
        0.66, 1, 'pb', 1, 0.1;          % Balance probability controlling different exploration strategies
        0, 1, 'pc', 1, 0.1;             % Balance probability controlling different exploration strategies
        1.1, 1.9, 'lambda', 1, 0.1;     % Stability factor used in the Levy flight
        0.05, 2, 'alpha2', 1, 0.3;      % Scaling factor
        0.05, 2, 'alpha3', 1, 0.3;      % Scaling factor
        0.05, 2, 'alpha4', 1, 0.3;      % Scaling factor
        0, 1, 'cr', 1, 0.1;             % Crossover rate
        
            
        0, 1, 'pa', 2, 0.1;        % Switching parameter for the replacement of worse solutions
        0, 1, 'pg', 2, 0.1;        % Launching probability of theGBlocalrandomwalkstrategy
        0.05, 2, 'alpha1', 2, 0.3; % Scaling factor
    };

    % the indices of corresponding control parameters
    % group I
    peIndex = 1;
    pbIndex = 2;
    pcIndex = 3;
    lambdaIndex = 4;
    alpha2Index = 5; 
    alpha3Index = 6;
    alpha4Index = 7;
    crIndex = 8;
    
    % group II
    paIndex = 1;
    pgIndex = 2;
    alpha1Index = 3;

    nHistory = params.nHistory;
    % a struct saving all detail information of the control parameters
    algParams = bsCreateHistoryStructNew(paramsConfigure, nHistory, nNest);
    nField = size(paramsConfigure, 1);              % the number of control parameters
    % used to record the mean values of all parameters at each iterations,
    % it will be outputed to find the change trend of each parameter
    paramHistory = zeros(nField, params.maxIter);   
    
    % when cost every intervalFEV function evalutions, the gradient-based
    % local search will be performed
    maxFEV = params.maxFunctionEvaluations;
    intervalFEV = maxFEV / params.interval;
    stopFEV = intervalFEV;
    % track the convergence path of each nest
    detailUpdates = cell(1, nNest);
    % set the id of each nest
    nests_ids = 1 : nNest;
    
    mid_results = [];
    
    %% ---------------------------------------------------------------------------------------------------------------%%
    % Starting iterations, see Algorithm 2
    % ----------------------------------------------------------------------------------------------------------------%%
    for iter = 1 : params.maxIter
        
        % calculate the current number of nests
        nNest = bsCalcCurNNest(params.minNests, params.maxNests, maxFEV, nfev);

        % calculate the new parameters based on the history information
        [algParams] = bsCalNewValueOfAlgParamsNew(algParams, nNest, iter);

        % record the statistic history information
        if params.isSaveMiddleRes
            st = 1;
            for i = 1 : 2
                t = algParams{i};
                paramHistory(st:st+t.nparam-1, iter) = mean(t.values, 2);
                st = st + t.nparam;
            end
        end
        
        
        
        % sort the nests with ID information of each nest, so that we can
        % trach how each nest (has a fixed ID) updates
        nests_with_ids = [nests_ids; nests];
        [nests_with_ids, gradient, fitness] = bsSortNestsWithId(nests_with_ids, gradient, fitness, nNest);
        nests = nests_with_ids(2:end, :);       % the sorted nests follows the ascending order of fitness
        nests_ids = nests_with_ids(1, 1:nNest);
        
        % save the middle iteration results
        if iter == 1
            % get the current best nest
            globalMinFVal = fitness(1);
            globalBestNest = nests(:, 1);
            xInit = globalBestNest;
        end
        
        if params.isSaveMiddleRes
            
            if iter == 1
                mid_results = [0; 0; globalMinFVal];
%             else
%                 mid_results = [mid_results, [iter; nfev; globalMinFVal]];
            end
        end
        
        if params.isSaveDetailUpdates
            for inest = 1 : nNest
                detailUpdates{nests_ids(inest)} = [detailUpdates{nests_ids(inest)}, [0; fitness(inest); nests(:, inest)]];
            end
        end
            
        % Generate new solutions (but keep the current bestNest)
        d1 = algParams{1}.values;
        alpha2s = randRepmat(d1(alpha2Index, :), 0, 1, 0.1, nDim);
        alpha3s = randRepmat(d1(alpha3Index, :), 0, 1, 0.1, nDim); 
        alpha4s = randRepmat(d1(alpha4Index, :), 0, 1, 0.1, nDim);
        newNest = bsGetCuckoosByThreeStrategiesNew(nests, Lb, Ub, ...
            alpha2s, ...
            alpha3s, ...
            alpha4s, ...
            d1(lambdaIndex, :), ...
            d1(pbIndex, :), ...
            d1(peIndex, :), ...
            d1(crIndex, :), ...
            d1(pcIndex, :)...
        );
        % calculate the gradient information, and fitness of the new
        % population, the old individuals will not be replaced if their new
        % individuals are worse than the old ones
        [nests, newFitness, gradient] = bsGetFitnessWithGradient(gradientObjFunc, nests, newNest, fitness, gradient);
        % record the feedback of the control parameters of the first group
        algParams = bsCalcScoreOfParametersNew(algParams, fitness, newFitness, 1, nHistory);
        fitness = newFitness;

        % Update the counter
        nfev = nfev + nNest; 
       
        d2 = algParams{2}.values;
        a1pha1s = randRepmat(d2(alpha1Index, :), 0, 1, 0.3, nDim);
        newNest = bsEmptyNestsTwoStrategiesNew(nests, gradient, Lb, Ub, ...
            d2(pgIndex, :), ...
            d2(paIndex, :), ...
            a1pha1s);
        % calculate the gradient information, and fitness of the new
        % population, the old individuals will not be replaced if their new
        % individuals are worse than the old ones
        [nests, newFitness, gradient] = bsGetFitnessWithGradient(gradientObjFunc, nests, newNest, fitness, gradient);
        % record the feedback of the control parameters of the second group
        algParams = bsCalcScoreOfParametersNew(algParams, fitness, newFitness, 2, nHistory);
        fitness = newFitness;
        
        % Update the counter
        nfev = nfev + nNest; 
        
        [~, index] = min(fitness);
        bestFVal = fitness(index);
        bestNest = nests(:, index);
        

        % record the global best solution obtained so far  
        if bestFVal < globalMinFVal
            globalMinFVal = bestFVal;
            globalBestNest = bestNest;
        end
        
        if ((nfev>stopFEV) || iter == 1)
            % call the gradient-based algorithm to find the nearest local
            % minima of the pbest solution
            
            
            nests_with_ids = [nests_ids; nests];
            [nests_with_ids, gradient, fitness] = bsSortNestsWithId(nests_with_ids, gradient, fitness, nNest);
            nests = nests_with_ids(2:end, :);
            nests_ids = nests_with_ids(1, 1:nNest);
        
            % randomly choose a top 3 solution, and then search its local
            % minima by GBSlover
            index = randi([1, 3]);
           
            [nests(:, index), fitness(index), ~, output] = bsGBSolveByOptions(gradientObjFunc, nests(:, index), Lb, Ub, GBOptions);
            gradient(:, index) = output.gradient;   % update gradient information
            nfev = nfev + output.funcCount;         % count the number of function evaluation called in the gradient-based local search

            % update the global information
            if fitness(index) < globalMinFVal
                globalMinFVal = fitness(index);
                globalBestNest = nests(:, index);
            end
            
            % the local search will be activated at next stopFEV iteration
            % (after intervalFEV iterations)
            stopFEV = stopFEV + intervalFEV;
        end

        if params.isSaveMiddleRes
            mid_results = [mid_results, [iter; nfev; globalMinFVal]];
        end
%         
        
        data.fNew = globalMinFVal;
        data.nfev = nfev;
        data.iter = iter;
        
        [exitFlag, output.message] = bsCheckStopCriteria(data, params);
        [data] = bsPlotMidResults(xInit, data, params, Lb, Ub, exitFlag > 0);
        
        if exitFlag > 0
            
            if params.isSaveDetailUpdates
                for inest = 1 : nNest
                    detailUpdates{nests_ids(inest)} = [detailUpdates{nests_ids(inest)}, [0; fitness(inest); nests(:, inest)]];
                end
            end
            
            break;
        end
        
    end %% End of iterations
    
    
    % save the search results
    paramHistory(:, iter+1:end) = [];
    x = globalBestNest;
    fval = globalMinFVal;
    output.funcCount = nfev;
    output.iterations = iter;
    output.midResults = mid_results;  
    output.frames = data.frames;
    output.paramHistory = paramHistory;
    output.detailUpdates = detailUpdates;
end

function params = randRepmat(params, minVal, maxVal, variation, nDim)   
    
    % parameters promotion strategy, in other words, for different nests,
    % different dimension, we use different value of control parameters.
    nNest = length(params);
    params = repmat(params, nDim, 1);

    params = params - variation + 2*variation*rand(nDim, nNest);

    params(params>maxVal) = maxVal;
    params(params<minVal) = minVal;
end

function [parameters] = bsCreateHistoryStructNew(paramsConfigure, nHistory, nNest)
    % create a structure data saving the detail information of control
    % parameters
    t = cell(1, 2);
    % the name of each parameter
    flags = cell2mat(paramsConfigure(:, 4));
    % the index of the parameters
    seqs = 1 : size(paramsConfigure, 1);
    
    % there is two group
    for i = 1 : 2
        t{i}.index = 0;     % index refers to the location where to save the current iteration in *.store
        t{i}.nTrain = 0;    % the number of training iterations that have been saved
        t{i}.iter = 0;      % the number of iterations
        index = flags == i; % get the logical index indicating whether the parameters belongs to the current group
        
        n = sum(index);
        t{i}.nparam = n;    % n is the number of control parameters in this group
        t{i}.names = paramsConfigure(index, 3); 
        % note that the second dimension of store is n+1, not n
        % this is because the last column of sotre will save the increased
        % fitness in each iteration
        t{i}.store = zeros(nHistory, n+1);  
        
        % save the range and variation of each control parameter
        t{i}.values = zeros(n, nNest);
        t{i}.max = zeros(n, 1);
        t{i}.min = zeros(n, 1);
        t{i}.variation = zeros(n, 1);
    
        seq = seqs(index);

        for j = 1 : n
            k = seq(j);
            
            minVal = paramsConfigure{k, 1};
            maxVal = paramsConfigure{k, 2};
            variation = paramsConfigure{k, 5};

            t{i}.min(j) = minVal;
            t{i}.max(j) = maxVal;
            t{i}.variation(j) = variation;
        end

        % initialize the values of all control parameters
        t{i}.values = bsGenValsByMinMax(t{i}.min, t{i}.max, nNest, t{i}.variation, t{i}.iter);
    end
    
    % return the parameters
    parameters = t;

end

function valSeq = bsGenValsByMinMax(minVal, maxVal, n, variation, iter)
    % randomly generate the values of a control parameter given its range and variation 
    m = length(minVal);

    value = minVal + (maxVal-minVal).*rand(m, 1);
    value = repmat(value, 1, n);
    variation = repmat(variation, 1, n);
    valSeq = value - variation + 2*variation.*rand(m, n);
    
    % limit the value sequences in the range of [min, max]
    valSeq = bsBetterBounds(valSeq, minVal, maxVal, []);
    
end

function algParams = bsCalcScoreOfParametersNew(algParams, fitness, newFitNess, flag, nHistory)

    % find that whose fitness is updated to a smaller value, which means better
    % solution achieved (the smaller fitness is, the better the solution is).
    deltaFit = fitness - newFitNess;
    deltaFit(deltaFit<0) = 0; % only consider fitness > newFitNess
   
    algParams{flag} = bsUpdateParameterNew(algParams{flag}, deltaFit, nHistory);
    
end

function t = bsUpdateParameterNew(t, deltaFit, nHistory)
    % update the history information saved in store based on the change of
    % fitness, 

    nNest = length(deltaFit);
    
    spos = t.index + 1;
    t.index = t.index + nNest;
    t.nTrain = t.nTrain + nNest;
    t.iter = t.iter + 1;
%     values = mean(t.values, 2);
% %     values = values';
    values = t.values';

    if (t.index <= nHistory)
        % if the number of history information hasn't been exceeds the size
        % of history store
        t.store(spos:t.index, end) = deltaFit;
        t.store(spos:t.index, 1:t.nparam) = values;
    else
        % when exceeds the history store, we cover the oldest history
        % information with the newest information
        n1 = nHistory - spos + 1;
        n2 = nNest - n1;

        t.store(spos:nHistory, end) = deltaFit(1:n1);
        t.store(spos:nHistory, 1:t.nparam) = values(1:n1, :);
        t.store(1:n2, end) = deltaFit(n1+1:end);
        t.store(1:n2, 1:t.nparam) = values(n1+1:end, :);
        
        t.index = n2;
    end
    
    if (t.nTrain > nHistory)
        t.nTrain = nHistory;
    end
    
end

function [algParams] = bsCalNewValueOfAlgParamsNew(algParams, nNest, iter)
    
    % there are two groups in algParams
    for i = 1 : length(algParams)
        % get the information of the i-th group
        t = algParams{i};
        
        nTrain = t.nTrain;  
        [~, nParam] = size(t.store);
        nParam = nParam - 1;

        if nTrain >= size(t.store, 1)
            % if the number of training iterations has been greater than
            % the size of history store
            
            % sort store with descending order of the last dimension (the increased fitness)
            sortData = sortrows(t.store(1:nTrain, :), -(nParam+1));
            % we only use the first five history records for obtaining the
            % parameters of next generation
            sortData = sortData(1:5, :);

            % sum the total increased fitness of the pass 5 best iterations
            sumWeight = sum(sortData(:, nParam+1));

            if sumWeight > 0
                % if the total increased fitness > 0, the value of
                % parameter in next genration is generated following the
                % equation of p_i = p_mean - p_var + 2 * p_var * rand(0,1)
                S = sortData(:, 1:nParam);
                value = mean(S, 1)';
                value = repmat(value, 1, nNest);
                variation = repmat(t.variation, 1, nNest);
                trueVal = value - variation + 2*variation.*rand(nParam, nNest);
                trueVal = bsBetterBounds(trueVal, t.min, t.max, []);
                t.values = trueVal;
            else
                % if there is no increased fitness during the pass
                % iterations, we just generate the parameters randomly
                % within the given range
                t.values = bsGenValsByMinMax(t.min, t.max, nNest, t.variation, t.iter);
            end
            
        else
            t.values = bsGenValsByMinMax(t.min, t.max, nNest, t.variation, t.iter);
        end
        
        
        algParams{i} = t;
    end
    
end

function newNests = bsGetCuckoosByThreeStrategiesNew(nests, Lb, Ub, alphas2, alphas3, alphas4, lambda, pbs, pes, crs, pcs)

    % Levy flights
    [nDim, nNest] = size(nests);

    %% Levy flights by Mantegna'nest algorithm
    crs = repmat(crs, nDim, 1);
    
    newNests = nests;
    
%     pbestNest = nests(:, 1);
    pbest = floor(pes * nNest);
    pbest(pbest < 1) = 1;

    % two types of the flighs, one is levy-flight-based, the other one is
    % DE-based
    methods = rand(1, nNest) <= pbs;
    
    % randomly generate levy flight for all levy-flight-based strategies
    steps = bsLevy(nDim, sum(methods), lambda(methods));
    it = 1;     % indicate the location of levy-flight-based walks

    for i = 1 : nNest
        nest = nests(:, i);

%         choose method
        if methods(i)
            % the following two strategies are levy-flight-based walk
            
            % randomly choose one pbest from nests (nests has been ordered by fitness)
            pbestNest = nests(:, randi(pbest(i)));

            if rand() <= pcs(i)
                % stratety 1
                stepsize = alphas3(:, i) .* steps(:, it) .* (nest - pbestNest);
                it = it + 1;
            else
                % stratety 2
                stepsize = alphas2(:, i) .* steps(:, it);
                it = it + 1;
            end
            
        else
            % this strategy is DE-based walk
            
            % stratety 3
            % randomly choose one pbest from nests (nests has been ordered by fitness)
            pbestNest = nests(:, randi(pbest(i)));
            
            r1 = randi(nNest);
            r2 = randi(nNest);
            
            stepsize = alphas4(:, i) .* (pbestNest - nest + nests(:, r1) - nests(:, r2));
        end
        
        % obtain the i-th new nests, it may come from arbitrary one of
        % the above three strategies
        newNests(:, i) = nest + stepsize;
    end
    
    % call bsBetterBounds to project the new nests into the given range
    newNests = bsBetterBounds(newNests, Lb, Ub, nests);
    
	% crossover ratio matrix
    K = rand(size(nests)) <= crs;

    % generate new nests with crossover ratio matrix K
    newNests = newNests .* K + nests .* (~K);
    
end

function newNests = bsEmptyNestsTwoStrategiesNew(nests, gradient, Lb, Ub, pg, pa, alpha2)

    % A fraction of worse nests are discovered with a probability pa
    [nDim, nNest] = size(nests);
    
    pa = repmat(pa, nDim, 1);

    % replacement parabability matrix    
    K = rand(size(nests)) <= pa;
    
    % exploration strategy choosing matrix
    U = repmat((rand(1, nNest) <= pg), nDim, 1);

    stepsize = alpha2 .* (nests(:, randperm(nNest)) - nests(:, randperm(nNest)));
    
    % new popluation generated without gradient information
    newNests2 = nests + stepsize .* K;
    % new population generated with gradient information
    newNests1 = nests + abs(stepsize) .* sign(-gradient) .* K;
    
    % combine two different new population together with exploration strategy choosing matrix U 
    newNests = newNests1 .* U + (~U) .* newNests2;
    
    % call bsBetterBounds to project the new nests into the given range
    newNests = bsBetterBounds(newNests, Lb, Ub, nests);

end

function [nests, gradient, fitness] = bsSortNestsWithId(nests, gradient, fitness, nNewNest)

    [nDim, nNest] = size(nests);
    % the first dimension is the ID information, we use the ID to trach how
    % each nest changes
    nDim = nDim - 1;
    
    piFi = [nests', gradient', fitness]; % the column looks like (nDim+1, nDim, 1)
    % Sort piFi with assending order of fitness (at the nDim*2+2-th column)
    piFiS = sortrows(piFi, nDim*2+2);

    % get the sorted nests with id information, and their corresponding
    % gradient and fitness.
    nests = piFiS(1:nNewNest, 1:nDim+1)';
    gradient = piFiS(1:nNewNest, end-nDim:end-1)';
    fitness = piFiS(1:nNewNest, end);

        
end

function nNest = bsCalcCurNNest(minN, maxN, maxNFE, NFE)
    % this funciton is set for population reduction strategy
    nNest = round((minN - maxN)/maxNFE*NFE + maxN);
end



function [x,fval,exitFlag,output] = bsHGBCSByShe2019(objFunc, Lb, Ub, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% This code implements the hybrid gradient based Cuckoo search (HGBCS) 
% algorithm by Bin She
%
% Programmed by Bin She (bin.stepbystep@gmail.com)
% Programmed dates: May 2019
%
% -------------------------------------------------------------------------
% The work is based on the following paper:
%
% 1) Fateen S E K, Bonilla-Petriciolet A. Gradient-based cuckoo search for 
% global optimization[J]. Mathematical Problems in Engineering, 2014, 2014.
%
% 2) X.-S. Yang, S. Deb, Cuckoo search via Levy flights,
% in: Proc. of World Congress on Nature & Biologically Inspired
% Computing (NaBIC 2009), December 2009, India,
% IEEE Publications, USA,  pp. 210-214 (2009).
% http://arxiv.org/PS_cache/arxiv/pdf/1003/1003.1594v1.pdf 
% -------------------------------------------------------------------------

    % parse some basic parameters for this process
    p = inputParser;
%     rng(125789);

    p = bsAddBasicalParameters(p, length(Lb));
    
    addParameter(p, 'nNest', 25 );      % The number of nests
    addParameter(p, 'pa', 0.5 );        % Discovery rate of alien eggs/solutions
    addParameter(p, 'beta', 3/2 );      % beta, for levy flight
    addParameter(p, 'alpha', 0.1 );     % for generating local solutions
    addParameter(p, 'interval', 5 );    % for each interval iterations, perform gradient descent algorithm once
    addParameter(p, 'innerMaxIter', 100 ); % the max number iterations for gradient descent algorithm
    
    % set default function of generation initial population as @bsGenerateInitialPopulationByLHS
    addParameter(p, 'initialPopulationFcn', @bsGenerateInitialPopulationByLHS );
    addParameter(p, 'optionsForGB', [] );
    % whether to save the detail update information of all population. The
    % value is set to yes when we need to display an animation of
    % optimization process in general.
    addParameter(p, 'isSaveDetailUpdates', false);
    
    p.parse(varargin{:});  
    params = p.Results;
    
    nNest = params.nNest;
    
    % call initial function to generate the initial population
    nests = params.initialPopulationFcn(Lb, Ub, nNest);
    
    % Get the current bestNest
    fitness = inf * ones(nNest, 1);
    nDim = length(Lb);
    gradient = zeros(nDim, nNest);
    
%     gradientObjFunc = @(x)(objFunc(x, 0));
    gradientObjFunc = objFunc;
    [globalMinFVal, globalBestNest, nests, fitness, gradient] = bsGetBestNestWithGradient(gradientObjFunc, nests, nests, fitness, gradient);
    
    nfev = nNest;   %count the number of evaluations of objective functions
    fs = [];
    xInit = globalBestNest;
    
    GBOptions = bsCreateGBOptions(length(Lb), ...
        'optimalFunctionTolerance', params.optimalFunctionTolerance, ...
        'maxIter', params.innerMaxIter, ...
        'isSaveMiddleRes', params.isSaveMiddleRes, ...
        'optimalF', params.optimalF);
    
    maxFEV = params.maxFunctionEvaluations;
    intervalFEV = maxFEV / params.interval;
    stopFEV = intervalFEV;
    fs = [0; 0; globalMinFVal];
    % track the convergence path of each nest
    detailUpdates = cell(1, nNest);
    if params.isSaveDetailUpdates
        for inest = 1 : nNest
            detailUpdates{inest} = [detailUpdates{inest}, [0; fitness(inest); nests(:, inest)]];
        end
    end
    
    %% Starting iterations
    for iter = 1 : params.maxIter
        
        % Generate new solutions (but keep the current bestNest)
        newNests = bsGetCuckoos(params.beta, nests, globalBestNest, Lb, Ub, params.alpha);   
        
        % get the bestNest solution among the new solutions
        [bestFVal, bestNest, nests, fitness, gradient] = bsGetBestNestWithGradient(gradientObjFunc, nests, newNests, fitness, gradient);
        
        % Update the counter
        nfev = nfev + nNest; 
          
        % Discovery and randomization
        newNests = bsEmptyNestsWithGradient(nests, gradient, Lb, Ub, params.pa) ;

        % Evaluate this set of solutions
        [bestFVal, bestNest, nests, fitness, gradient, K] = bsGetBestNestWithGradient(gradientObjFunc, nests, newNests, fitness, gradient);
        
        % Update the counter
        nfev = nfev + nNest; 
            
         if nfev > stopFEV
            [bestNest, bestFVal, ~, output] = bsGBSolveByOptions(gradientObjFunc, bestNest, Lb, Ub, GBOptions);
            % K is the index of the bestNest in the population nests
            nests(:, K) = bestNest;
            fitness(K) = bestFVal;
            gradient(:,K) = output.gradient;
            
            % update the number of function evaluations taken
            nfev = nfev + output.funcCount;
            
            stopFEV = stopFEV + intervalFEV;
        end

        
        % Find the bestNest objective so far  
        if bestFVal < globalMinFVal
            globalMinFVal = bestFVal;
            globalBestNest = bestNest;
        end
        
        % Update the counter again
        if params.isSaveMiddleRes
            fs = [fs, [iter; nfev; globalMinFVal]];
        end
        
        if params.isSaveDetailUpdates
            for inest = 1 : nNest
                detailUpdates{inest} = [detailUpdates{inest}, [0; fitness(inest); nests(:, inest)]];
            end
        end
        
        data.fNew = globalMinFVal;
        data.nfev = nfev;
        data.iter = iter;
        
        exitFlag = bsCheckStopCriteria(data, params);
        [data] = bsPlotMidResults(xInit, data, params, Lb, Ub, exitFlag > 0);
        
        if exitFlag > 0
            break;
        end
        
        
    end %% End of iterations

    x = globalBestNest;
    fval = globalMinFVal;
    output.funcCount = nfev;
    output.iterations = iter;
    output.midResults = fs;  
    output.frames = data.frames;
    output.detailUpdates = detailUpdates;
end

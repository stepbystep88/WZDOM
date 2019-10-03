function [x,fval,exitFlag,output] = bsGBCSByFateen2014(objFunc, Lb, Ub, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% This code implements the gradient based Cuckoo search (GBCS) algorithm by Fateen
% Programmed by Bin She (bin.stepbystep@gmail.com)
% Programmed dates: May 2019
%
% -------------------------------------------------------------------------
% see the follow paper for details:
% 1) Fateen S E K, Bonilla-Petriciolet A. Gradient-based cuckoo search for 
% global optimization[J]. Mathematical Problems in Engineering, 2014, 2014.
% -------------------------------------------------------------------------

    % parse some basic parameters for this process
    p = inputParser;
%     rng(125789);

    p = bsAddBasicalParameters(p, length(Lb));
    
    addParameter(p, 'nNest', 25 );     % The number of nests
    addParameter(p, 'pa', 0.5 );       % Discovery rate of alien eggs/solutions
    addParameter(p, 'beta', 3/2 );
    addParameter(p, 'alpha', 0.1 );
    % set default function of generation initial population as @bsGenerateInitialPopulationByLHS
    addParameter(p, 'initialPopulationFcn', @bsGenerateInitialPopulationByLHS );
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
    gradientObjFunc = @(x)(objFunc(x, 1));
    
    [globalMinFVal, globalBestNest, nests, fitness, gradient] = bsGetBestNestWithGradient(gradientObjFunc, nests, nests, fitness, gradient);
    xInit = globalBestNest;
    
    nfev = nNest;   %count the number of evaluations of objective functions
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
        
        % Update the counter again
        nfev = nfev + nNest;
        
        
        % Find the bestNest objective so far  
        if bestFVal < globalMinFVal
            globalMinFVal = bestFVal;
            globalBestNest = bestNest;
        end
        
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
        
%         data.fOld = globalMinFVal;
%         data.xOld = globalBestNest;
        
    end %% End of iterations

    x = globalBestNest;
    fval = globalMinFVal;
    output.funcCount = nfev;
    output.iterations = iter;
    output.midResults = fs;    
    output.frames = data.frames;   
    output.detailUpdates = detailUpdates;
end





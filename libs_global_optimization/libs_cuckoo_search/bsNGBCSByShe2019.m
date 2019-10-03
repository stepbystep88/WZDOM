function [x,fval,exitFlag,output] = bsNGBCSByShe2019(objFunc, Lb, Ub, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% This code implements the novel gradient based Cuckoo search (NGBCS) 
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
    addParameter(p, 'gbStepsize', 0.01 ); % stepsize for gradient descent algorithm
    
    % set default function of generation initial population as @bsGenerateInitialPopulationByLHS
    addParameter(p, 'initialPopulationFcn', @bsGenerateInitialPopulationByLHS );
    addParameter(p, 'optionsForGB', [] );
    
    p.parse(varargin{:});  
    params = p.Results;
    
    nNest = params.nNest;
    
    % call initial function to generate the initial population
    nests = params.initialPopulationFcn(Lb, Ub, nNest);
    
    % Get the current bestNest
    fitness = 10^10 * ones(nNest, 1);
    noGradientObjFunc = @(x)(objFunc(x, 0));
    gradientObjFunc = @(x)(objFunc(x, 1));
    [globalMinFVal, globalBestNest, nests, fitness] = bsGetBestNest(noGradientObjFunc, nests, nests, fitness);
    
    nfev = nNest;   %count the number of evaluations of objective functions
    xInit = globalBestNest;
    
%     data.fOld = globalMinFVal;
%     data.xOld = globalBestNest;
   
    GBOptions = bsCreateGBOptions(length(Lb), ...
        'optimalFunctionTolerance', params.optimalFunctionTolerance, ...
        'maxIter', params.innerMaxIter, ...
        'isSaveMiddleRes', params.isSaveMiddleRes, ...
        'optimalF', params.optimalF);
    
    %% Starting iterations
    numIter = 0; % this iter inclues the nubmer of iterations performed in the gradient based optimization process.
%     midResults.x = [];
    fs = [];
    
    for iter = 1 : params.maxIter
        
        numIter = numIter + 1;
        
        % Generate new solutions (but keep the current bestNest)
        newNests = bsGetCuckoos(params.beta, nests, globalBestNest, Lb, Ub, params.alpha);   
        
        % get the bestNest solution among the new solutions
        [bestFVal, bestNest, nests, fitness] = bsGetBestNest(noGradientObjFunc, nests, newNests, fitness);
        
        % Update the counter
        nfev = nfev + nNest; 
          
        % Discovery and randomization
        newNests = bsEmptyNests(nests, Lb, Ub, params.pa) ;

        % Evaluate this set of solutions
        [bestFVal, bestNest, nests, fitness, K] = bsGetBestNest(noGradientObjFunc, nests, newNests, fitness);
        
        
        % Update the counter again
        nfev = nfev + nNest;
        
        % Update the counter again
        if params.isSaveMiddleRes
            fs = [fs, [nfev; min(bestFVal, globalMinFVal)]];
        end
        
        
        if mod(iter, params.interval) == 0 && bestFVal < globalMinFVal
            
            % it has to be different from the current globalBestNest
            [bestNest, bestFVal, ~, output] = bsGBSolveByOptions(gradientObjFunc, bestNest, Lb, Ub, GBOptions);
            
%             K = randi(nNest);
            nests(:, K) = bestNest;
            fitness(K) = bestFVal;
            
            nfev = nfev + output.funcCount;
            numIter = numIter + output.iterations;
            
            if params.isSaveMiddleRes
                fs = [fs, [nfev; bestFVal]];
            end
            
        end
            
        
        % Find the bestNest objective so far  
        if bestFVal < globalMinFVal
            globalMinFVal = bestFVal;
            globalBestNest = bestNest;
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
    output.iterations = numIter;
    output.midResults = fs;    
    output.frames = data.frames;
end





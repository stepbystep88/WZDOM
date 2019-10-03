function [x,fval,exitFlag,output] = bsCSByYangAndDeb2009(objFunc, Lb, Ub, varargin)
%% This code just re-organizes the code implemented by Xin-She Yang and Suash Deb
% Organized by Bin She (bin.stepbystep@gmail.com)
% Organizing dates: May 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -----------------------------------------------------------------
% Cuckoo Search (CS) algorithm by Xin-She Yang and Suash Deb      %
% Programmed by Xin-She Yang at Cambridge University              %
% Programming dates: Nov 2008 to June 2009                        %
% Last revised: Dec  2009   (simplified version for demo only)    %
% -----------------------------------------------------------------
% Papers -- Citation Details:
% 1) X.-S. Yang, S. Deb, Cuckoo search via Levy flights,
% in: Proc. of World Congress on Nature & Biologically Inspired
% Computing (NaBIC 2009), December 2009, India,
% IEEE Publications, USA,  pp. 210-214 (2009).
% http://arxiv.org/PS_cache/arxiv/pdf/1003/1003.1594v1.pdf 
% 2) X.-S. Yang, S. Deb, Engineering optimization by cuckoo search,
% Int. J. Mathematical Modelling and Numerical Optimisation, 
% Vol. 1, No. 4, 330-343 (2010). 
% http://arxiv.org/PS_cache/arxiv/pdf/1005/1005.2908v2.pdf
% ----------------------------------------------------------------%
% This demo program only implements a standard version of         %
% Cuckoo Search (CS), as the Levy flights and generation of       %
% new solutions may use slightly different methods.               %
% The pseudo code was given sequentially (select a cuckoo etc),   %
% but the implementation here uses Matlab'nest vector capability,    %
% which results in neater/better codes and shorter running time.  % 
% This implementation is different and more efficient than the    %
% the demo code provided in the book by 
%    "Yang X. S., Nature-Inspired Metaheuristic Algoirthms,       % 
%     2nd Edition, Luniver Press, (2010).                 "       %
% --------------------------------------------------------------- %
% =============================================================== %
% Notes:                                                          %
% Different implementations may lead to slightly different        %
% behavour and/or results, but there is nothing wrong with it,    %
% as this is the nature of random walks and all metaheuristics.   %
% -----------------------------------------------------------------

    % parse some basic parameters for this process
    p = inputParser;
%     rng(125789);

    p = bsAddBasicalParameters(p, length(Lb));
    
    addParameter(p, 'nNest', 25 );     % The number of nests
    addParameter(p, 'pa', 0.5 );       % Discovery rate of alien eggs/solutions
    addParameter(p, 'beta', 3/2 );
    addParameter(p, 'alpha', 3/2 );
    % set default function of generation initial population as @bsGenerateInitialPopulationByLHS
    addParameter(p, 'initialPopulationFcn', @bsGenerateInitialPopulationByRandom );
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
    [globalMinFVal, globalBestNest, nests, fitness] = bsGetBestNest(objFunc, nests, nests, fitness);
    xInit = globalBestNest;
    
    nfev = nNest;   %count the number of function evaluations
    
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
        [bestFVal, bestNest, nests, fitness] = bsGetBestNest(objFunc, nests, newNests, fitness);
        
        % Update the counter
        nfev = nfev + nNest; 
          
        % Discovery and randomization
        newNests = bsEmptyNests(nests, Lb, Ub, params.pa) ;

        % Evaluate this set of solutions
        [bestFVal, bestNest, nests, fitness] = bsGetBestNest(objFunc, nests, newNests, fitness);
        
        
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
        
        
    end %% End of iterations

    x = globalBestNest;
    fval = globalMinFVal;
    output.funcCount = nfev;
    output.iterations = iter;
    output.midResults = fs;
    output.frames = data.frames;
    output.detailUpdates = detailUpdates;
end




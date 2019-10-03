function [x,fval,exitFlag,output] = bsMCSByWalton2011(objFunc, Lb, Ub, varargin)
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
    
    addParameter(p, 'nNest', 25 );      % The number of nests
    addParameter(p, 'pamax', 1 );        % Discovery rate of alien eggs/solutions
    addParameter(p, 'pamin', 0.25 );        % Discovery rate of alien eggs/solutions
    addParameter(p, 'beta', 3/2 );      % beta, for levy flight
    addParameter(p, 'alphamin', 0.02 );         % Step size factor, increase this to decrease step size
    addParameter(p, 'alphamax', 0.5 );         % Step size factor, increase this to decrease step size
%     addParameter(p, 'pwr', 0.5 );       % Power that step size is reduced by each generation
    addParameter(p, 'nesD', 1 );        % Number of eggs deleted at each generation
    addParameter(p, 'minNests', 10 );
    % set default function of generation initial population as @bsGenerateInitialPopulationByLHS
    addParameter(p, 'initialPopulationFcn', @bsGenerateInitialPopulationByRandom );
    
    p.parse(varargin{:});  
    params = p.Results;
    
    nNest = params.nNest;
    
    % call initial function to generate the initial population
    nests = params.initialPopulationFcn(Lb, Ub, nNest);
    
    % Get the current bestNest
    fitness = inf * ones(nNest, 1);
    nDim = length(Lb);

    %1 - Calculate fitness for initial nests
    [globalMinFVal, globalBestNest, nests, fitness] = bsGetBestNest(objFunc, nests, nests, fitness);
    
    nfev = nNest;   %count the number of evaluations of objective functions
    fs = [];
    
    xInit = globalBestNest;
%     A = (Ub - Lb) ./ params.A;
%     pa = params.pa;               
%     ptop = 1 - pa;   
    beta = params.beta;

    %% Starting iterations
    for iter = 1 : params.maxIter
        
        % a) sort the current nests in order of fitness
    
        % First put vectors into matrix form
        piFi = [nests', fitness];
        % Sort by Fi in assending order
        piFiS = sortrows(piFi, nDim+1);
    
        % Decrease number of nests, we only need lots of nests initially to get
        % a good sampling of the objective function
        nNest = max(params.minNests, nNest - params.nesD);
        nests = piFiS(1:nNest, 1:nDim)';
        fitness = piFiS(1:nNest, nDim+1);
%         nTop = max(3, round(nNest*ptop));     %%%%%%%%%%%%%%%
    
%         a = A ./ (iter ^ (params.pwr));
%         
%         a = 0.05;
        
        
        if iter > 500
            pa = params.pamin;
%             a = params.pamin;
        else
            pa = params.pamax - iter/500 * (params.pamax - params.pamin);
%             c = 1/100 * log(params.alphamin / params.alphamax);
%             a = params.alphamax * exp(c * iter);
        end
        
        % 2. Loop over each Cuckoo which has been discarded
        % Generate new solutions (but keep the current bestNest)
        newNests = bsGetCuckoosByChaos(beta, nests, globalBestNest, Lb, Ub);   
        
        % get the bestNest solution among the new solutions
        [bestFVal, bestNest, nests, fitness] = bsGetBestNest(objFunc, nests, newNests, fitness);
        
        % Update the counter
        nfev = nfev + nNest; 
        
        % Discovery and randomization
        newNests = bsEmptyNests(nests, Lb, Ub, pa) ;
        
        % Evaluate this set of solutions
        [bestFVal, bestNest, nests, fitness] = bsGetBestNest(objFunc, nests, newNests, fitness);
        
        % Find the bestNest objective so far  
        if bestFVal < globalMinFVal
            globalMinFVal = bestFVal;
            globalBestNest = bestNest;
        end
        
        % Update the counter again
        nfev = nfev + nNest;
        
        % Update the counter again
        if params.isSaveMiddleRes
            fs = [fs, [nfev; globalMinFVal]];
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



function [outImpedances, regParam] = bsSeisInv1DMultiTraces(regFlag, ds, G, initImpedans, Lb, Ub, regParam, parampkgs, options, isParallel)
    
    [nSampNum, nTrace] = size(initImpedans);
    
    %% deal with different cases of input arguments
    if ~exist('Lb', 'var')
        Lb = [];
    else 
        if size(Lb, 2) == 1
            Lb = repmat(Lb, 1, nTrace);
        end
        
        Lb = log(Lb);
    end
    
    if ~exist('Ub', 'var')
        Ub = [];
    else
        if size(Ub, 2) == 1
            Ub = repmat(Ub, 1, nTrace);
        end
        
        Ub = log(Ub);
    end
    
    if ~exist('regParam', 'var')
        regParam = [];
    end
    
    if ~exist('parampkgs', 'var') || isempty(parampkgs)
        parampkgs = [];
    end
    
    if ~exist('options', 'var') || isempty(options)
        options = bsCreateSeisInv1DOptions(nSampNum);
    end
    
    if ~exist('isParallel', 'var') || isempty(isParallel)
        isParallel = true;
    end
    
    if nTrace ~= size(initImpedans, 2)
        error('The number of traces of obseved seismic data ds is not consistent with that of initial impedances.');
    end
    
    options.GBOptions.display = 'off';
    options.GBOptions.plotFcn = [];
    options.GBOptions.isSaveMiddleRes = false;
    
    % inversion process is performed in logrithm domain
    initXs = log(initImpedans);
    outXs = zeros(nSampNum, nTrace);
    
    % if the regularization parameter is not assigned, we search the best
    % parameter at first
    % select the regularization parameter from random 10 traces
    nSelectRegTrace = 10;
    randIndex = randi(nTrace, 1, nSelectRegTrace);
    regParams = [];
    invertedTraces = zeros(1, nTrace);
    
    fprintf('Selecting regularization parameter from random 10 traces...\n');
    
    for i = 1 : nSelectRegTrace
        iTrace = randIndex(i);
        
        if isempty(Lb) 
            [outXs(:, iTrace), fval, exitFlag, output] = bsSeisInv1DTrace(regFlag, ds(:, iTrace), G, initXs(:, iTrace), ...
                [], [], regParam, parampkgs, options);
        else
            [outXs(:, iTrace), fval, exitFlag, output] = bsSeisInv1DTrace(regFlag, ds(:, iTrace), G, initXs(:, iTrace), ...
                Lb(:, iTrace), Ub(:, iTrace), regParam, parampkgs, options);
        end
        
        regParams = [regParams, output.regParam];
        % even if the regParam is not empty, we still inverse one trace at
        % first is because we need to update the parampkgs (initial some
        % difference matrix, so that the initial processing will not be
        % performed in the inversion of the other traces).
        parampkgs = output.parampkgs;
        invertedTraces(iTrace) = 1;
    end
  
    if isa(regParams(1), 'struct')
        fields = fieldnames(regParams(1));
        avgRegParam = regParams(1);
        
        T = struct2table(regParams);  % turn it into a table
        for i = 1 : length(fields)
            avgRegParam = setfield(avgRegParam, fields{i}, sum(T{:, i}) / nSelectRegTrace);
        end
        regParam = avgRegParam;
    else
        regParam = mean(regParams);
    end
    
    
    if isParallel
        parfor iTrace = 1 : nTrace
            
            if invertedTraces(iTrace) 
                continue;
            end
            
            outXs(:, iTrace) = bsSeisInv1DTrace(regFlag, ds(:, iTrace), G, initXs(:, iTrace), Lb, Ub, regParam, parampkgs, options);
            fprintf('Runing %dth trace using method %s...\n', iTrace, regFlag);
        end
    else
        
        
        for iTrace = 1 : nTrace
            
            if invertedTraces(iTrace) 
                continue;
            end
            
            outXs(:, iTrace) = bsSeisInv1DTrace(regFlag, ds(:, iTrace), G, initXs(:, iTrace), Lb, Ub, regParam, parampkgs, options);
            invertedTraces(iTrace) = 1;
            nFinishedTrace = sum(invertedTraces);
            
            fprintf('Runing %dth trace using method %s, in progress (%d/%d)...\n', iTrace, regFlag, nFinishedTrace, nTrace);
        end
    end
    
    outImpedances = exp(outXs);
    
end
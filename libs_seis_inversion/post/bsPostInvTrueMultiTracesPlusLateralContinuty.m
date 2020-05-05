function [invResults] = bsPostInvTrueMultiTracesPlusLateralContinuty(GInvParam, inIds, crossIds, timeLine, methods)
%% inverse multiple traces
% Programmed by: Bin She (Email: bin.stepbystep@gmail.com)
% Programming dates: Nov 2019
% -------------------------------------------------------------------------
%
% Input 
% GInvParam     all information of the inverse task
% inIds      	inline ids to be inverted
% crossIds      crossline ids to be inverted
% timeLine      horizon information
% methods       the methods to solve the inverse task, each method is a
% struct depicting the details of the methods. For example, method.load
% indicates whether load the result from mat or segy directly.
% 
% Output
% invVals       inverted results, a cell includ
% -------------------------------------------------------------------------
    
    assert(length(inIds) == length(crossIds), 'The length of inline ids and crossline ids must be the same.');
    nMethod = size(methods, 1);
    traceNum = length(inIds);
    sampNum = GInvParam.upNum + GInvParam.downNum;
    % save the inverted results
    rangeIn = [min(inIds), max(inIds)];
    rangeCross = [min(crossIds), max(crossIds)];
    
    % horion of the whole volume
    usedTimeLine = timeLine{GInvParam.usedTimeLineId};
    
    % create folder to save the intermediate results
    try
        warning('off');
        mkdir([GInvParam.modelSavePath,'/models/']);
        mkdir([GInvParam.modelSavePath,'/mat_results/']);
        mkdir([GInvParam.modelSavePath,'/sgy_results/']);
        warning('on');
    catch
    end
    
    invResults = cell(1, nMethod);
    % horizon of given traces
    horizonTimes = bsGetHorizonTime(usedTimeLine, inIds, crossIds, GInvParam.isParallel, GInvParam.numWorkers);
    
    startTimes = horizonTimes - GInvParam.dt * GInvParam.upNum;
    
%     startTimes = bsButtLowPassFilter(startTimes, 0.1);
    
    for i = 1 : nMethod
        method = methods{i};
        methodName = method.name;
        matFileName = bsGetFileName('mat');
        segyFileName = bsGetFileName('segy');
        
        res.source = [];
        
        if isfield(method, 'load')
            loadInfo = method.load;
            switch loadInfo.mode
                % load results directly
                case 'mat'
                    try
                        % from mat file
                        if isfield(loadInfo, 'fileName') && ~isempty(loadInfo.fileName)
                            load(GInvParam.load.fileName);
                        else
                            load(matFileName);
                        end

                        res.source = 'mat';
                    catch
                        warning('load mat file failed.');
                    end
                    
                case 'segy'
                    % from sgy file
                    try
                        if ~isfield(loadInfo, 'fileName') || isempty(loadInfo.fileName)
                            loadInfo.fileName = segyFileName;
                        end

                        [data, loadInfo.segyInfo] = bsReadTracesByIds(...
                            loadInfo.fileName, ...
                            loadInfo.segyInfo, ...
                            inIds, crossIds, startTimes, sampNum, GInvParam.dt);

                        res.source = 'segy';
                    catch
                        warning('load segy file failed.');
                    end
                    
                case 'assign'
                    data = loadInfo.data;
                    res.source = 'assign';
            end
            
        end
        
        if isempty(res.source)
            % obtain results by computing
            data = bsCallInvFcn();
            res.source = 'compution';
        end
        
        res.data = data;            % inverted results
        res.inIds = inIds;          % inverted inline ids 
        res.crossIds = crossIds;    % inverted crossline ids
        res.horizon = horizonTimes; % the horizon of inverted traces
        res.name = method.name;     % the name of this method
        res.dt = GInvParam.dt;
        res.upNum = GInvParam.upNum;
        res.downNum = GInvParam.downNum;
        
        if isfield(method, 'type')  % the type of inverted volume, IP, Seismic, VP, VS, Rho
            res.type = method.type;
        else
            res.type = 'Ip';
        end
        
        if isfield(method, 'showFiltCoef')
            res.showFiltCoef = method.showFiltCoef;
        else
            res.showFiltCoef = 0;
        end
        
        invResults{i} = res;
        
        
        % save sgy file
        if isfield(method, 'isSaveSegy') && method.isSaveSegy && ~strcmp(res.source, 'segy')
            fprintf('Writing segy file:%s ....\n', segyFileName);
            bsWriteInvResultIntoSegyFile( ...
                res, data, ...
                GInvParam.postSeisData.fileName, ...
                GInvParam.postSeisData.segyInfo, ...
                segyFileName);
            fprintf('Write segy file:%s successfully!\n', segyFileName);
        end
        
        % save mat file
        if isfield(method, 'isSaveMat') && method.isSaveMat && ~strcmp(res.source, 'mat')
            fprintf('Writing mat file:%s...\n', matFileName);
            try
                save(matFileName, 'data', 'horizonTimes', 'inIds', 'crossIds', 'GInvParam', 'method');
            catch
                save(matFileName, 'data', 'horizonTimes', 'inIds', 'crossIds', 'GInvParam', 'method', '-v7.3');
            end
            fprintf('Write mat file:%s successfully!\n', matFileName);
        end
    end
    
    function fileName = bsGetFileName(type)
        switch type
            case 'mat'
                fileName = sprintf('%s/mat_results/Ip_%s_inline_[%d_%d]_crossline_[%d_%d].mat', ...
                    GInvParam.modelSavePath, methodName, rangeIn(1), rangeIn(2), rangeCross(1), rangeCross(2));
            case 'segy'
                fileName = sprintf('%s/sgy_results/Ip_%s_inline_[%d_%d]_crossline_[%d_%d].sgy', ...
                    GInvParam.modelSavePath, methodName, rangeIn(1), rangeIn(2), rangeCross(1), rangeCross(2));
        end
        
    end

    function data = bsCallInvFcn()
        % tackle the inverse task
%         data = zeros(sampNum, traceNum);
        GSParam = bsInitDLSRPkgs(method.parampkgs, sampNum);
        options = method.options;
        seisOption = GInvParam.seisInvOptions;
        GBOptions = seisOption.GBOptions;
        mainFunc = seisOption.mainFunc;
        initRegParam = options.initRegParam;
        models = cell(1, traceNum);
        xs = zeros(sampNum, traceNum);
        avg_xs = zeros(sampNum, traceNum);
        gamma = method.regParam.gamma;
        lambda = method.regParam.lambda;
        GBOptions.maxIter = options.innerIter;
        neiboors = cell(1, traceNum);
        
        if isfield(GSParam, 'nMultipleTrace')
            KTrace = GSParam.nMultipleTrace;
        else
            KTrace = 1;
        end
        
        [rangeInline, rangeCrossline] = bsGetWorkAreaRangeByParam(GInvParam);
        nRangeInline = rangeInline(2) - rangeInline(1) + 1;
        nRangeCrossline = rangeCrossline(2) - rangeCrossline(1) + 1;

        if nRangeInline > nRangeCrossline
            nTracePerLine = nRangeInline;
        else
            nTracePerLine = nRangeCrossline;
        end
    
        pbm = bsInitParforProgress(GInvParam.numWorkers, ...
            traceNum, ...
            sprintf('Post inversion progress information by method %s', method.name), ...
            GInvParam.modelSavePath, ...
            GInvParam.isPrintBySavingFile);
            
            
        % 第一步获取所欲的模型
        preModel = bsPostPrepareModel(GInvParam, inIds(1), crossIds(1), horizonTimes(1), [], []);
            
        
        models{1} = preModel;
        xs(:, 1) = models{1}.initX;
        
        for iTrace = 2 : traceNum
            models{iTrace} = bsPostPrepareModel(GInvParam, inIds(iTrace), crossIds(iTrace), horizonTimes(iTrace), [], preModel);
            xs(:, iTrace) = models{iTrace}.initX;
            bsIncParforProgress(pbm, iTrace, 100);
        end
        
        parfor iTrace = 1 : traceNum
            % 找当前当的所有邻近道
            neiboors{iTrace} = bsFindNearestKTrace(iTrace, inIds, crossIds, KTrace, nTracePerLine);
        end
        
        % 第二步：大循环
        for iter = 1 : options.maxIter
            fprintf("第%d次迭代：常规反演\n", iter);
            parfor iTrace = 1 : traceNum
                xs(:, iTrace) = invNormalOneTrace(xs(:, iTrace), models{iTrace}, mainFunc, lambda, initRegParam, GBOptions);
            end
            
            Ips = exp(xs);
            
            fprintf("第%d次迭代：稀疏重构\n", iter);
            for iTrace = 1 : traceNum
                
                avg_xs(:, iTrace) = sparseRebuildOneTrace(GSParam, Ips(:, neiboors{iTrace}));
            end
            
            fprintf("第%d次迭代：反演结果合并\n", iter);
            parfor iTrace = 1 : traceNum
                xNew = Ips(:, iTrace) * (1 - gamma) + avg_xs(:, iTrace) * gamma;
                xs(:, iTrace) = log(xNew);
            end
            
        end
        
        bsDeleteParforProgress(pbm);
        
        data = exp(xs);

    end
    
end

function xOut = invNormalOneTrace(xInit, model, mainFunc, lambda, initRegParam, GBOptions)
    inputObjFcnPkgs = {
        mainFunc,    struct('A', model.G, 'B', model.d),   1; 
        @bsReg1DTKInitModel,    struct('xInit', []), lambda;
        @bsReg1DTKInitModel,    struct('xInit', model.initX), initRegParam;
    };
    [xOut, ~, ~, ~] = bsGBSolveByOptions(inputObjFcnPkgs, xInit, [], [], GBOptions);
end

function avgLog = sparseRebuildOneTrace(GSParam, input)
    % sparse reconstruction

    nBlock = size(input, 2);
    sizeAtom = GSParam.sizeAtom;
    ncell = GSParam.ncell;
    
    all_patches = zeros(sizeAtom*nBlock, ncell);
    patches = zeros(sizeAtom, ncell);
    
    for k = 1 : nBlock
        sPos = (k-1)*sizeAtom + 1;
        ePos = sPos + sizeAtom - 1;
        
        for j = 1 : ncell
            js = GSParam.index(j);
            patches(:, j) = input(js : js+sizeAtom-1, k);
        end
        
        all_patches(sPos:ePos, :) = patches;
        
    end
   
    gammas = omp(GSParam.ODIC'*all_patches, GSParam.omp_G, GSParam.sparsity);
    new_patches = GSParam.DIC *  gammas;
    
    avgLog = bsAvgPatches(new_patches, GSParam.index, size(input, 1));
    
end


function GSParam = bsInitDLSRPkgs(GSParam, sampNum)
    
    if isfield(GSParam, 'omp_G')
        return;
    end

    validatestring(string(GSParam.reconstructType), {'equation', 'simpleAvg'});

    
    trainDICParam = GSParam.trainDICParam;
    
    sizeAtom = trainDICParam.sizeAtom;
    nAtom= trainDICParam.nAtom;
    GSParam.sizeAtom = sizeAtom;
    GSParam.nAtom = nAtom;
    GSParam.nrepeat = sizeAtom - GSParam.stride;
    
    index = 1 : GSParam.stride : sampNum - sizeAtom + 1;
    if(index(end) ~= sampNum - sizeAtom + 1)
        index = [index, sampNum - sizeAtom + 1];
    end
    
    GSParam.index = index;
    GSParam.ncell = length(index);
    
    if isfield(GSParam, 'nMultipleTrace')
        nBlock = GSParam.nMultipleTrace;
    else
        nBlock = 1;
    end
    
    GSParam.ODIC = repmat(GSParam.DIC, nBlock, 1);
    
    for i = 1 : size(GSParam.ODIC, 2)
        tmp = norm(GSParam.ODIC(:, i));
        GSParam.ODIC(:, i) = GSParam.ODIC(:, i) / tmp;
        GSParam.DIC(:, i) = GSParam.DIC(:, i) / tmp;
    end
    
    GSParam.omp_G = GSParam.ODIC' * GSParam.ODIC;

end


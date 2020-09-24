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
            [data, estimatedError] = bsCallInvFcn();
            res.source = 'compution';
        end
        
        if exist('estimatedError', 'var')
            res.error = estimatedError;
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
                save(matFileName, 'data', 'res', 'horizonTimes', 'inIds', 'crossIds', 'method');
            catch
                save(matFileName, 'data', 'res', 'horizonTimes', 'inIds', 'crossIds', 'method', '-v7.3');
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

    function [data, ys] = bsCallInvFcn()
        % tackle the inverse task
%         data = zeros(sampNum, traceNum);
        
%         models = cell(1, traceNum);
        xs = zeros(sampNum, traceNum);
        ds = zeros(sampNum-1, traceNum);
        scaleFactors = zeros(1, traceNum);
        
        neiboors = cell(1, traceNum);
        
        if isfield(method.parampkgs, 'nMultipleTrace')
            KTrace = method.parampkgs.nMultipleTrace;
        else
            KTrace = 1;
        end
        
        try
            [rangeInline, rangeCrossline] = bsGetWorkAreaRangeByParam(GInvParam);
            nRangeInline = rangeInline(2) - rangeInline(1) + 1;
            nRangeCrossline = rangeCrossline(2) - rangeCrossline(1) + 1;

            if nRangeInline > nRangeCrossline
                nTracePerLine = nRangeInline;
            else
                nTracePerLine = nRangeCrossline;
            end
        catch
            nTracePerLine = max(crossIds) - min(crossIds) + 1;
        end
    
        pbm = bsInitParforProgress(GInvParam.numWorkers, ...
            traceNum, ...
            '', ...
            GInvParam.modelSavePath, ...
            GInvParam.isPrintBySavingFile);
            
            
        % 第一步获取所欲的模型
        preModel = bsPostPrepareModel(GInvParam, inIds(1), crossIds(1), horizonTimes(1), [], []);
            
        
%         models{1} = preModel;
        xs(:, 1) = preModel.initX;
        ds(:, 1) = preModel.d;
        scaleFactors(1) = preModel.scaleFactor;
        
        pbm.title = sprintf('Prepare model... %s', method.name);
        
        parfor iTrace = 2 : traceNum
            model = bsPostPrepareModel(GInvParam, inIds(iTrace), crossIds(iTrace), horizonTimes(iTrace), [], preModel);
%             xs(:, iTrace) = models{iTrace}.initX;
            xs(:, iTrace) = model.initX;
            ds(:, iTrace) = model.d;
            scaleFactors(iTrace) = model.scaleFactor;
            
            bsIncParforProgress(pbm, iTrace, 100);
        end
        
%         ds = ds .* repmat(scaleFactors, sampNum-1, 1);
%         GST_options = bsCreateGSTParam(2, 'sigma', 4, 'iterNum', 20);
%         ds = bsSmoothByGST2D(ds, [], GST_options);
%         ds = ds ./ repmat(scaleFactors, sampNum-1, 1);
        
        parfor iTrace = 1 : traceNum
            % 找当前当的所有邻近道
            neiboors{iTrace} = bsFindNearestKTrace(iTrace, inIds, crossIds, KTrace, nTracePerLine);
        end
        
        switch method.flag
            case 'DLSR'
                [data, ys] = bsPostInvMultiTracesByDLSR(GInvParam, neiboors, ds, preModel.orginal_G, xs, scaleFactors, inIds, crossIds, method);
            case {'DLSR-DEM', 'DLSR-EM'}
                [data, ys] = bsPostInvMultiTracesByDLSR_EM(GInvParam, neiboors, ds, preModel.orginal_G, xs, scaleFactors, inIds, crossIds, method);
        end
    end
    
end






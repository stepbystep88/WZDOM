function [invResults] = bsPostInvTrueMultiTraces(GInvParam, inIds, crossIds, timeLine, methods)
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
        
        % create folder to save the intermediate results
        try
            if ~isfield(method, 'load') || strcmpi(method.load.mode, 'off') || (~exist(method.load.fileName, 'file') && ~exist(matFileName, 'file'))

                warning('off');
                mkdir([GInvParam.modelSavePath, '/', methodName, '/mat_results/']);
                mkdir([GInvParam.modelSavePath, '/', methodName, '/sgy_results/']);
                warning('on');
            end
        catch
        end
    
        res.source = [];
        
        if isfield(method, 'load')
            loadInfo = method.load;
            switch loadInfo.mode
                % load results directly
                case 'mat'
                    try
                        % from mat file
                        if isfield(loadInfo, 'fileName') && ~isempty(loadInfo.fileName)
                            load(loadInfo.fileName);
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
        
        if isfield(res, 'pickle')
            pickle = res.pickle;
        else
            pickle = [];
        end
        
        if isempty(res.source)
            % obtain results by computing
%             data = bsCallInvFcn();
            
            % obtain results by computing
            if (strcmpi(method.load.mode, 'off') && any(strcmpi({'DLSR-GST', 'DLSR-EM', 'DLSR-DEM'}, method.flag))) || ...
                    ( isfield(method, 'parampkgs') ...
                    && isfield(method.parampkgs, 'nMultipleTrace') ...
                    && method.parampkgs.nMultipleTrace > 1) 
                
                if ~isfield(method.parampkgs, 'nMultipleTrace')
                    method.parampkgs.nMultipleTrace = 1;
                end
                
                KTrace = method.parampkgs.nMultipleTrace;
                [data, pickle] = bsCallInvFcnMultiTraces(KTrace);
    
            else
                data = bsCallInvFcnTraceByTrace();
            end
            
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
        res.pickle = pickle;
        
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
        if isfield(method, 'isSaveSegy') && method.isSaveSegy && ~strcmp(res.source, 'segy') && ~strcmp(res.source, 'mat')
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
                fileName = sprintf('%s/%s/mat_results/Ip_%s_inline_[%d_%d]_crossline_[%d_%d].mat', ...
                    GInvParam.modelSavePath, methodName, methodName, rangeIn(1), rangeIn(2), rangeCross(1), rangeCross(2));
            case 'segy'
                fileName = sprintf('%s/%s/sgy_results/Ip_%s_inline_[%d_%d]_crossline_[%d_%d].sgy', ...
                    GInvParam.modelSavePath, methodName, methodName, rangeIn(1), rangeIn(2), rangeCross(1), rangeCross(2));
        end
        
    end

    function data = bsCallInvFcnTraceByTrace()
        % tackle the inverse task
        data = zeros(sampNum, traceNum);
        
        
        % obtain a preModel avoid calculating matrix G again and again.
        % see line 20 of function bsPostPrepareModel for details
        [data(:, 1), preModel, output] = bsPostInvOneTrace(GInvParam, horizonTimes(1), method, inIds(1), crossIds(1), [], 0, 1);
        if ~isempty(output)
            method.parampkgs = output.parampkgs;
        end
        
        if GInvParam.isParallel
            
            pbm = bsInitParforProgress(GInvParam.numWorkers, ...
                traceNum, ...
                sprintf('Post inversion progress information by method %s', method.name), ...
                GInvParam.modelSavePath, ...
                GInvParam.isPrintBySavingFile);
            
            % parallel computing
            parfor iTrace = 2 : traceNum
                
                data(:, iTrace) = bsPostInvOneTrace(GInvParam, horizonTimes(iTrace), method, inIds(iTrace), crossIds(iTrace), preModel, 0, iTrace);
                
                bsIncParforProgress(pbm, iTrace, 100);
            end
            

        else
            % non-parallel computing 
            for iTrace = 2 : traceNum
                data(:, iTrace) = bsPostInvOneTrace(GInvParam, horizonTimes(iTrace), method, inIds(iTrace), crossIds(iTrace), preModel, 1, iTrace);
            end
        end
    end
    
    function [data, ys] = bsCallInvFcnMultiTraces(KTrace)
        % tackle the inverse task
        
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
            nTracePerLine = max(max(crossIds) - min(crossIds), max(inIds) - min(inIds)) + 1;
        end
    
        pbm = bsInitParforProgress(GInvParam.numWorkers, ...
            traceNum, ...
            '', ...
            GInvParam.modelSavePath, ...
            GInvParam.isPrintBySavingFile);
        
        try
            [~,name,~] = fileparts(GInvParam.initModel.ip.fileName);
        catch
            name = '';
        end
        
        tmp_file_name = sprintf('%s/mid_nTrace_%d_inIds_[%d, %d]_crossIds_[%d %d]_filtCoef_%.2f_init_model_name_%s.mat', ...
            GInvParam.modelSavePath, length(inIds), min(inIds), max(inIds), ...
            min(crossIds), max(crossIds), GInvParam.initModel.filtCoef, name);
        
        % 第一步获取所欲的模型
        preModel = bsPostPrepareModel(GInvParam, inIds(1), crossIds(1), horizonTimes(1), [], []);
            
        if exist(tmp_file_name, 'file')
            load(tmp_file_name, 'xs', 'ds', 'scaleFactors');
        else
            xs = zeros(sampNum, traceNum);
            ds = zeros(sampNum-1, traceNum);
            scaleFactors = zeros(1, traceNum);

            neiboors = cell(1, traceNum);
        
    %         models{1} = preModel;
            xs(:, 1) = preModel.initX;
            ds(:, 1) = preModel.d;
            scaleFactors(1) = preModel.scaleFactor;

            pbm = bsResetParforProgress(pbm, sprintf('Preparing model... %s', method.name));
            parfor iTrace = 2 : traceNum
                model = bsPostPrepareModel(GInvParam, inIds(iTrace), crossIds(iTrace), horizonTimes(iTrace), [], preModel);
    %             xs(:, iTrace) = models{iTrace}.initX;
                xs(:, iTrace) = model.initX;
                ds(:, iTrace) = model.d;
                scaleFactors(iTrace) = model.scaleFactor;

                bsIncParforProgress(pbm, iTrace, 1000);
            end
            try
                save(tmp_file_name, 'xs', 'ds', 'scaleFactors');
            catch
                save(tmp_file_name, 'xs', 'ds', 'scaleFactors', '-v7.3');
            end
            
        end
        
        
        pbm = bsResetParforProgress(pbm, sprintf('Finding neiboors... %s', method.name));
        parfor iTrace = 1 : traceNum
            % 找当前当的所有邻近道
            neiboors{iTrace} = bsFindNearestKTrace(iTrace, inIds, crossIds, KTrace, nTracePerLine);
            bsIncParforProgress(pbm, iTrace, 1000);
        end    
        bsDeleteParforProgress(pbm);    
        
        
        switch method.flag
            case 'DLSR'
                [data, ys] = bsPostInvMultiTracesByDLSR(GInvParam, neiboors, ds, preModel.orginal_G, xs, scaleFactors, inIds, crossIds, method);
            case {'DLSR-EM', 'DLSR-DEM'}
                [data, ys] = bsPostInvMultiTracesByDLSR_EM(GInvParam, neiboors, ds, preModel.orginal_G, xs, scaleFactors, inIds, crossIds, method);
            case 'DLSR-GST'
                [data, ys] = bsPostInvMultiTracesByDLSR_GST(GInvParam, neiboors, ds, preModel.orginal_G, xs, scaleFactors, inIds, crossIds, method);
        end
        
        
    end
    
end


function fileName = bsGetModelFileName(modelSavePath, inId, crossId)
    % get the file name of model (to save the )
    fileName = sprintf('%s/models/model_inline_%d_crossline_%d.mat', ...
        modelSavePath, inId, crossId);

end

function [idata, model, output] = bsPostInvOneTrace(GInvParam, horizonTime, method, inId, crossId, preModel, isprint, iTrace)

    if isprint
        fprintf('Solving the trace of inline=%d and crossline=%d by using method %s...\n', ...
            inId, crossId, method.name);
    end
    
%     try
        % create model data
        if GInvParam.isReadMode
            % in read mode, model is loaded from local file
            modelFileName = bsGetModelFileName(GInvParam.modelSavePath, inId, crossId);
            parLoad(modelFileName);
        else
            % otherwise, create the model by computing. Note that we input
            % argment preModel is a pre-calculated model, by doing this, we
            % avoid wasting time on calculating the common data of different
            % traces such as forward matrix G.
            model = bsPostPrepareModel(GInvParam, inId, crossId, horizonTime, [], preModel);
            if GInvParam.isSaveMode
                % in save mode, mode should be saved as local file
                modelFileName = bsGetModelFileName(GInvParam.modelSavePath, inId, crossId);
                parSave(modelFileName, model);
            end
        end

        method.options.inline = inId;
        method.options.crossline = crossId;
        method.options.scaleFactor = model.scaleFactor;
        
        if isa(method.regParam, 'function_handle')
            method.regParam = method.regParam(iTrace);
        end
        
        [xOut, ~, ~, output] = bsPostInv1DTrace(model.d, model.G, model.initX, model.Lb, model.Ub, method);                       

        idata = exp(xOut);
        
%     catch err
%         fprintf(getReport(err));
%         
%         sampNum = GInvParam.upNum + GInvParam.downNum;
%         idata = zeros(sampNum, 1);
%         model = [];
%         output = [];
%     end
    
end


function [invResults] = bsPreInvTrueMultiTraces(GInvParam, inIds, crossIds, timeLine, methods)
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
    horizonTimes = bsGetHorizonTime(usedTimeLine, inIds, crossIds, ...
            GInvParam.isParallel, GInvParam.numWorkers);
        
    if ~isempty(GInvParam.smooth_horizon_fcn) && ~isempty(GInvParam.smooth_horizon_fcn)
        horizonTimes = GInvParam.smooth_horizon_fcn(horizonTimes);
    end
    
    startTimes = horizonTimes - GInvParam.dt * GInvParam.upNum;
    
    
    for i = 1 : nMethod
        method = methods{i};
        methodName = method.name;
        matFileName = bsGetFileName('mat');
        
        % create folder to save the intermediate results
        try
            if ~isfield(method, 'load') || strcmpi(method.load.mode, 'off') || (ischar(method.load.fileName) && ~exist(method.load.fileName, 'file') && ~exist(matFileName, 'file'))
                warning('off');
                mkdir([GInvParam.modelSavePath, methodName, '/mat_results/']);
                mkdir([GInvParam.modelSavePath, methodName, '/sgy_results/']);
                warning('on');
            end
        catch
            warning('off');
            mkdir([GInvParam.modelSavePath, methodName, '/mat_results/']);
            mkdir([GInvParam.modelSavePath, methodName, '/sgy_results/']);
            warning('on');
        end
        
        % try loading the results from mat or segy file
        res = loadResults();
        
        
        if isempty(res.source)
            % obtain results by computing
            if ( isfield(method, 'flag') && any(strcmpi({'CSR-GST', 'SSR-GST'}, method.flag)) ) || ...
                    (isfield(method, 'parampkgs') ...
                    && isfield(method.parampkgs, 'nMultipleTrace') ...
                    && method.parampkgs.nMultipleTrace > 1 ...
                    && any(strcmpi({'CSR', 'CSR-EM'}, method.flag)))
                
                KTrace = method.parampkgs.nMultipleTrace;
                [vp, vs, rho] = bsCallInvFcnMultiTraces(KTrace);
            else
                [vp, vs, rho] = bsCallInvFcnTraceByTrace();
            end
        
            
            res.source = 'computation';
            res.type = {'vp', 'vs', 'rho'};
            res.data = {vp, vs, rho};
        end
        
        if isfield(method, 'showFiltCoef')
            res.showFiltCoef = method.showFiltCoef;
        else
            res.showFiltCoef = 0;
        end
        
        res.inIds = inIds;          % inverted inline ids 
        res.crossIds = crossIds;    % inverted crossline ids
        res.horizon = horizonTimes; % the horizon of inverted traces
        res.name = method.name;     % the name of this method
        res.dt = GInvParam.dt;
        res.upNum = GInvParam.upNum;
        res.downNum = GInvParam.downNum;
        
        invResults{i} = res;
        
        % save the results into mat or segy file
        saveResults();
        
    end
    
    function res = loadResults()
        res.source = [];
        
        if isfield(method, 'load')
            loadInfo = method.load;
            switch loadInfo.mode
                % load results directly
                case 'assign'
                    res.source = 'assign';
                    res.type = loadInfo.type;
                    res.data = loadInfo.data;
                case 'mat'
                    try
                        % from mat file
                        if isfield(loadInfo, 'fileName') && ~isempty(loadInfo.fileName)
                            load_mat = load(GInvParam.load.fileName);
                        else
                            load_mat = load(matFileName);
                        end

                        res.source = 'mat';
                        res.type = {'vp', 'vs', 'rho'};
                        vp = load_mat.vp;
                        vs = load_mat.vs;
                        rho = load_mat.rho;
                        res.data = {vp, vs, rho};
                        
                    catch
                        warning('load mat file failed.');
                    end
                    
                case 'segy'
                    % from sgy file
                    if ~iscell(loadInfo.fileName)
                        loadInfo.fileName = {loadInfo.fileName};
                        nFile = 1;
                    else
                        nFile = length(loadInfo.fileName);
                    end
                    data = cell(1, nFile);
                    for iFile = 1 : nFile
                        if length(loadInfo.segyInfo) == 1
                            % share the common segy basic info
                            iSegyInfo = loadInfo.segyInfo;
                        else
                            % use different segy info for different files
                            iSegyInfo = loadInfo.segyInfo{iFile};
                        end
                        
                        if isfield(loadInfo, 'prestack') && loadInfo.prestack == 1
                            data{iFile} = bsStackPreSeisData(loadInfo.fileName{iFile}, ...
                                iSegyInfo, ...
                                inIds, crossIds, startTimes, sampNum, GInvParam.dt);
                        else
                            [data{iFile}, ~] = bsReadTracesByIds(...
                                loadInfo.fileName{iFile}, ...
                                iSegyInfo, ...
                                inIds, crossIds, startTimes, sampNum, GInvParam.dt);
                        end
                    end
                    
                    res.source = 'segy';
                    res.type = method.type;
                    if nFile == 1
                        res.data = data{1};
                    else
                        res.data = data;
                    end
                    
                    
            end
            
        end
    end

    function saveResults()
        % save mat file
        if isfield(method, 'isSaveMat') && method.isSaveMat && strcmp(res.source, 'computation')
            fprintf('Writing mat file:%s...\n', matFileName);
            tmethod.name = method.name;
            tmethod.flag = method.flag;
            
            try
                save(matFileName, 'vp', 'vs', 'rho', 'GInvParam', 'inIds', 'crossIds', 'horizonTimes', 'tmethod');
            catch
                save(matFileName, 'vp', 'vs', 'rho', 'GInvParam', 'inIds', 'crossIds', 'horizonTimes', 'tmethod', '-v7.3');
            end
            fprintf('Write mat file:%s successfully!\n', matFileName);
        end
        
        % save sgy file
        if isfield(method, 'isSaveSegy') && method.isSaveSegy 
            for k = 1 : 3
                switch k
                    case 1
                        segyFileName = bsGetFileName('segy', 'vp');
                        kdata = vp;
                    case 2
                        segyFileName = bsGetFileName('segy', 'vs');
                        kdata = vs;
                    case 3
                        segyFileName = bsGetFileName('segy', 'rho');
                        kdata = rho;
                end
                fprintf('Writing segy file:%s ....\n', segyFileName);
                bsWriteInvResultIntoSegyFile(res, kdata, ...
                    GInvParam.preSeisData.fileName, ...
                    GInvParam.preSeisData.segyInfo, ...
                    segyFileName, 0);
                fprintf('Write segy file:%s successfully!\n', segyFileName);
            end
            res.data = [];
        end
        
        
    end

    function fileName = bsGetFileName(type, attName)
        switch type
            case 'mat'
                fileName = sprintf('%s/%s/mat_results/vpsrho_%s_inline_[%d_%d]_crossline_[%d_%d].mat', ...
                    GInvParam.modelSavePath, methodName, methodName, rangeIn(1), rangeIn(2), rangeCross(1), rangeCross(2));
            case 'segy'
                fileName = sprintf('%s/%s/sgy_results/%s_%s_inline_[%d_%d]_crossline_[%d_%d].sgy', ...
                    GInvParam.modelSavePath, methodName, attName, methodName, rangeIn(1), rangeIn(2), rangeCross(1), rangeCross(2));
        end
        
    end

    function [vp, vs, rho] = bsCallInvFcnTraceByTrace()
        % tackle the inverse task
        vp = zeros(sampNum, traceNum);
        vs = zeros(sampNum, traceNum);
        rho = zeros(sampNum, traceNum);
        
        % obtain a preModel avoid calculating matrix G again and again.
        % see line 20 of function bsPrePrepareModel for details
        [vp(:, 1), vs(:, 1), rho(:, 1), preModel, output] = bsPreInvOneTrace(GInvParam, horizonTimes(1), method, inIds(1), crossIds(1), [], 0);
        if ~isempty(output)
            method.parampkgs = output.parampkgs;
        end
        
        if GInvParam.isParallel
            
            pbm = bsInitParforProgress(GInvParam.numWorkers, ...
                traceNum, ...
                sprintf('Pre inversion progress information by method %s', method.name), ...
                GInvParam.modelSavePath, ...
                GInvParam.isPrintBySavingFile);
            
            % parallel computing
            parfor iTrace = 2 : traceNum
%                     gpuDevice(1);
                    [vp(:, iTrace), vs(:, iTrace), rho(:, iTrace)] ...
                        = bsPreInvOneTrace(GInvParam, horizonTimes(iTrace), method, ...
                            inIds(iTrace), crossIds(iTrace), preModel, 0);

                    bsIncParforProgress(pbm, iTrace, 101);
            end
            
            bsDeleteParforProgress(pbm);
        else
            % non-parallel computing 
            for iTrace = 2 : traceNum
                [vp(:, iTrace), vs(:, iTrace), rho(:, iTrace)] ...
                    = bsPreInvOneTrace(GInvParam, horizonTimes(iTrace), method, ...
                        inIds(iTrace), crossIds(iTrace), preModel, 1);
            end
        end
    end
    
    function [vp, vs, rho] = bsCallInvFcnMultiTraces(KTrace)
        % tackle the inverse task
        neiboors = cell(1, traceNum);
        
        nTracePerLine = max(max(crossIds) - min(crossIds), max(inIds) - min(inIds)) + 1;
    
        pbm = bsInitParforProgress(GInvParam.numWorkers, ...
            traceNum, ...
            '', ...
            GInvParam.modelSavePath, ...
            GInvParam.isPrintBySavingFile);
            
        tmp_file_name = sprintf('%s/mid_nTrace_%d_inIds_[%d, %d]_crossIds_[%d %d]_filtCoef_%.2f_angleNum_%d.mat', ...
            GInvParam.modelSavePath, length(inIds), min(inIds), max(inIds), ...
            min(crossIds), max(crossIds), GInvParam.initModel.filtCoef, GInvParam.angleTrNum);
        
        if exist(tmp_file_name, 'file')
            try
                load(tmp_file_name, 'xs', 'ds', 'lsdCoefs', 'scaleFactors', 'G');
            catch
            end
        end
        
        if ~exist('ds', 'var')
            xs = zeros(sampNum*3, traceNum);
            ds = zeros( (sampNum-1)*GInvParam.angleTrNum, traceNum);
            scaleFactors = zeros(1, traceNum);
%             Gs = cell(1, traceNum);
            lsdCoefs = cell(1, traceNum);
            
            % 第一步获取所需的模型
            preModel = bsPrePrepareModel(GInvParam, inIds(1), crossIds(1), horizonTimes(1), [], []);

            xs(:, 1) = preModel.initX;
            ds(:, 1) = preModel.d;
    %         Gs{1} = preModel.G;
            G = preModel.orginal_G;
            scaleFactors(1) = preModel.scaleFactor;
            lsdCoefs{1} = preModel.lsdCoef;

            pbm.title = sprintf('Prepare model... %s', method.name);


            parfor iTrace = 2 : traceNum
                model = bsPrePrepareModel(GInvParam, inIds(iTrace), crossIds(iTrace), horizonTimes(iTrace), [], preModel);
    %             xs(:, iTrace) = models{iTrace}.initX;
                xs(:, iTrace) = model.initX;
                ds(:, iTrace) = model.d;
    %             Gs{iTrace} = model.G;
                lsdCoefs{iTrace} = model.lsdCoef;
                scaleFactors(iTrace) = model.scaleFactor;

                bsIncParforProgress(pbm, iTrace, 100);
            end
            try
                save(tmp_file_name, 'xs', 'ds', 'lsdCoefs', 'scaleFactors', 'G', '-v7.3');
            catch
                save(tmp_file_name, 'xs', 'ds', 'lsdCoefs', 'scaleFactors', 'G', '-v7.3');
            end
            
        end
        
        parfor iTrace = 1 : traceNum
            % 找当前当的所有邻近道
            neiboors{iTrace} = bsFindNearestKTrace(iTrace, inIds, crossIds, KTrace, nTracePerLine);
        end
        
        switch method.flag
            case 'CSR'
                
                [vp, vs, rho] = bsPreInvMultiTracesByCSR(GInvParam, neiboors, ds, G, xs, scaleFactors, lsdCoefs, inIds, crossIds, method);
            case 'CSR-GST'
                horizonTime = bsGetHorizonTime(timeLine{GInvParam.usedTimeLineId}, inIds, crossIds, 1);
                startTime = horizonTime - GInvParam.upNum * GInvParam.dt;
                postSeisData = bsGetPostSeisData(GInvParam, inIds, crossIds, startTime, sampNum);
                shiftedData = bsPhase90Shift(postSeisData);
    
                [vp, vs, rho] = bsPreInvMultiTracesByCSR_GST(GInvParam, neiboors, ds, G, xs, scaleFactors, lsdCoefs, shiftedData, inIds, crossIds, method);
%             case 'DLSR-EM'
%                 [data, ys] = bsPostInvMultiTracesByDLSR_EM(GInvParam, neiboors, ds, preModel.orginal_G, xs, scaleFactors, inIds, crossIds, method);
        end
    end
end


function fileName = bsGetModelFileName(modelSavePath, inId, crossId)
    % get the file name of model (to save the )
    fileName = sprintf('%s/models/model_inline_%d_crossline_%d.mat', ...
        modelSavePath, inId, crossId);

end

function [vp, vs, rho, model, output] = bsPreInvOneTrace(GInvParam, horizonTime, method, inId, crossId, preModel, isprint)

    if isprint
        fprintf('Solving the trace of inline=%d and crossline=%d by using method %s...\n', ...
            inId, crossId, method.name);
    end
    
    try
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
            model = bsPrePrepareModel(GInvParam, inId, crossId, horizonTime, [], preModel);
            if GInvParam.isSaveMode
                % in save mode, model should be saved as local file
                modelFileName = bsGetModelFileName(GInvParam.modelSavePath, inId, crossId);
                parSave(modelFileName, model);
            end
        end

        method.mode = GInvParam.mode;
        method.lsdCoef = model.lsdCoef;
        method.options.inline = inId;
        method.options.crossline = crossId;
        [xOut, ~, ~, output] = bsPreInv1DTrace(model.d, model.G, model.initX, model.Lb, model.Ub, method);                       

        [vp, vs, rho] = bsPreRecoverElasticParam(xOut, GInvParam.mode, model.lsdCoef);
    catch err
        fprintf(getReport(err));
        fprintf('The trace of inId=%d, crossId=%d has errors.\n', inId, crossId);
        sampNum = GInvParam.upNum + GInvParam.downNum;
        vp = zeros(sampNum, 1);
        vs = zeros(sampNum, 1);
        rho = zeros(sampNum, 1);
        model = [];
        output = [];
    end
end

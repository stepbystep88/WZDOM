function [invResults] = bsPreInvTrueMultiTraces(GPreInvParam, inIds, crossIds, timeLine, methods)
%% inverse multiple traces
% Programmed by: Bin She (Email: bin.stepbystep@gmail.com)
% Programming dates: Nov 2019
% -------------------------------------------------------------------------
%
% Input 
% GPreInvParam     all information of the inverse task
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
    sampNum = GPreInvParam.upNum + GPreInvParam.downNum;
    % save the inverted results
    rangeIn = [min(inIds), max(inIds)];
    rangeCross = [min(crossIds), max(crossIds)];
    
    % horion of the whole volume
    usedTimeLine = timeLine{GPreInvParam.usedTimeLineId};
    
    % create folder to save the intermediate results
    try
        mkdir([GPreInvParam.modelSavePath,'/models/']);
        mkdir([GPreInvParam.modelSavePath,'/mat_results/']);
        mkdir([GPreInvParam.modelSavePath,'/sgy_results/']);
    catch
    end
    
    invResults = cell(1, nMethod);
    
    
    % horizon of given traces
    horizonTimes = bsGetHorizonTime(usedTimeLine, inIds, crossIds, ...
            GPreInvParam.isParallel, GPreInvParam.numWorkers);
        
    startTimes = horizonTimes - GPreInvParam.dt * GPreInvParam.upNum;
    
    for i = 1 : nMethod
        method = methods{i};
        methodName = method.name;
        matFileName = bsGetFileName('mat');
        
        % try loading the results from mat or segy file
        res = loadResults();
        
        if isempty(res.source)
            % obtain results by computing
            [vp, vs, rho] = bsCallInvFcn();
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
        res.dt = GPreInvParam.dt;
        res.upNum = GPreInvParam.upNum;
        res.downNum = GPreInvParam.downNum;
        
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
                case 'mat'
                    try
                        % from mat file
                        if isfield(loadInfo, 'fileName') && ~isempty(loadInfo.fileName)
                            load_mat = load(GPreInvParam.load.fileName);
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
                            iSegyInfo = loadInfo.segyInfo(iFile);
                        end
                        
                        if isfield(loadInfo, 'prestack') && loadInfo.prestack == 1
                            data{iFile} = bsStackPreSeisData(loadInfo.fileName{iFile}, ...
                                iSegyInfo, ...
                                inIds, crossIds, startTimes, sampNum, GPreInvParam.dt);
                        else
                            [data{iFile}, ~] = bsReadTracesByIds(...
                                loadInfo.fileName{iFile}, ...
                                iSegyInfo, ...
                                inIds, crossIds, startTimes, sampNum, GPreInvParam.dt);
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
                    GPreInvParam.preSeisData.fileName, ...
                    GPreInvParam.preSeisData.segyInfo, ...
                    segyFileName);
                fprintf('Write segy file:%s successfully!\n', segyFileName);
            end
            res.data = [];
        end
        
        % save mat file
        if isfield(method, 'isSaveMat') && method.isSaveMat && strcmp(res.source, 'computation')
            fprintf('Writing mat file:%s...\n', matFileName);
            try
                save(matFileName, 'vp', 'vs', 'rho', 'GPreInvParam', 'inIds', 'crossIds', 'horizonTimes', 'method');
            catch
                save(matFileName, 'vp', 'vs', 'rho', 'GPreInvParam', 'inIds', 'crossIds', 'horizonTimes', 'method', '-v7.3');
            end
            fprintf('Write mat file:%s successfully!\n', matFileName);
        end
    end

    function fileName = bsGetFileName(type, attName)
        switch type
            case 'mat'
                fileName = sprintf('%s/mat_results/vpsrho_%s_inline_[%d_%d]_crossline_[%d_%d].mat', ...
                    GPreInvParam.modelSavePath, methodName, rangeIn(1), rangeIn(2), rangeCross(1), rangeCross(2));
            case 'segy'
                fileName = sprintf('%s/sgy_results/%s_%s_inline_[%d_%d]_crossline_[%d_%d].sgy', ...
                    GPreInvParam.modelSavePath, attName, methodName, rangeIn(1), rangeIn(2), rangeCross(1), rangeCross(2));
        end
        
    end

    function [vp, vs, rho] = bsCallInvFcn()
        % tackle the inverse task
        vp = zeros(sampNum, traceNum);
        vs = zeros(sampNum, traceNum);
        rho = zeros(sampNum, traceNum);
        
        % obtain a preModel avoid calculating matrix G again and again.
        % see line 20 of function bsPrePrepareModel for details
        [vp(:, 1), vs(:, 1), rho(:, 1), preModel, output] = bsPreInvOneTrace(GPreInvParam, horizonTimes(1), method, inIds(1), crossIds(1), [], 0);
        if ~isempty(output)
            method.parampkgs = output.parampkgs;
        end
        
        if GPreInvParam.isParallel
            
            pbm = bsInitParforProgress(GPreInvParam.numWorkers, ...
                traceNum, ...
                sprintf('Pre inversion progress information by method %s', method.name), ...
                GPreInvParam.modelSavePath, ...
                GPreInvParam.isPrintBySavingFile);
            
            % parallel computing
            parfor iTrace = 2 : traceNum
                    [vp(:, iTrace), vs(:, iTrace), rho(:, iTrace)] ...
                        = bsPreInvOneTrace(GPreInvParam, horizonTimes(iTrace), method, ...
                            inIds(iTrace), crossIds(iTrace), preModel, 0);

                    bsIncParforProgress(pbm, iTrace, 101);
            end
            

        else
            % non-parallel computing 
            for iTrace = 2 : traceNum
                [vp(:, iTrace), vs(:, iTrace), rho(:, iTrace)] ...
                    = bsPreInvOneTrace(GPreInvParam, horizonTimes(iTrace), method, ...
                        inIds(iTrace), crossIds(iTrace), preModel, 1);
            end
        end
    end
    
end


function fileName = bsGetModelFileName(modelSavePath, inId, crossId)
    % get the file name of model (to save the )
    fileName = sprintf('%s/models/model_inline_%d_crossline_%d.mat', ...
        modelSavePath, inId, crossId);

end

function [vp, vs, rho, model, output] = bsPreInvOneTrace(GPreInvParam, horizonTime, method, inId, crossId, preModel, isprint)

    if isprint
        fprintf('Solving the trace of inline=%d and crossline=%d by using method %s...\n', ...
            inId, crossId, method.name);
    end
    
    try
        % create model data
        if GPreInvParam.isReadMode
            % in read mode, model is loaded from local file
            modelFileName = bsGetModelFileName(GPreInvParam.modelSavePath, inId, crossId);
            parLoad(modelFileName);
        else
            % otherwise, create the model by computing. Note that we input
            % argment preModel is a pre-calculated model, by doing this, we
            % avoid wasting time on calculating the common data of different
            % traces such as forward matrix G.
            model = bsPrePrepareModel(GPreInvParam, inId, crossId, horizonTime, [], preModel);
            if GPreInvParam.isSaveMode
                % in save mode, model should be saved as local file
                modelFileName = bsGetModelFileName(GPreInvParam.modelSavePath, inId, crossId);
                parSave(modelFileName, model);
            end
        end

        method.mode = GPreInvParam.mode;
        method.lsdCoef = model.lsdCoef;
        [xOut, ~, ~, output] = bsPreInv1DTrace(model.d, model.G, model.initX, model.Lb, model.Ub, method);                       

        [vp, vs, rho] = bsPreRecoverElasticParam(xOut, GPreInvParam.mode, model.lsdCoef);
    catch err
        fprintf(getReport(err));
        
        sampNum = GPreInvParam.upNum + GPreInvParam.downNum;
        vp = zeros(sampNum, 1);
        vs = zeros(sampNum, 1);
        rho = zeros(sampNum, 1);
        model = [];
        output = [];
    end
end

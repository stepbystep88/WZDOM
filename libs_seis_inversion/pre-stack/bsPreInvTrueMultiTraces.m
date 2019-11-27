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
    if GPreInvParam.isParallel
        horizonTimes = bsCalcHorizonTimeParallel(usedTimeLine, inIds, crossIds, GPreInvParam.numWorkers);
    else
        horizonTimes = bsCalcHorizonTime(usedTimeLine, inIds, crossIds);
    end
    
    startTimes = horizonTimes - GPreInvParam.dt * GPreInvParam.upNum;
    
    for i = 1 : nMethod
        method = methods{i};
        methodName = method.name;
        matFileName = bsGetFileName('mat');
        
        res.source = [];
        
        if isfield(method, 'load')
            loadInfo = method.load;
            switch loadInfo.mode
                % load results directly
                case 'mat'
                    try
                        % from mat file
                        if isfield(loadInfo, 'fileName') && ~isempty(loadInfo.fileName)
                            load(GPreInvParam.load.fileName);
                        else
                            load(matFileName);
                        end

                        res.source = 'mat';
                    catch
                        warning('load mat file failed.');
                    end
                    
                case 'segy'
                    % from sgy file
                    [vp, loadInfo.segyInfo, ~] = bsReadTracesByIds(...
                        loadInfo.vp.fileName, ...
                        loadInfo.vp.segyInfo, ...
                        inIds, crossIds, startTimes, sampNum, GPreInvParam.dt);
                    
                    [vs, loadInfo.segyInfo, ~] = bsReadTracesByIds(...
                        loadInfo.vs.fileName, ...
                        loadInfo.vs.segyInfo, ...
                        inIds, crossIds, startTimes, sampNum, GPreInvParam.dt);
                    
                    [rho, loadInfo.segyInfo, ~] = bsReadTracesByIds(...
                        loadInfo.rho.fileName, ...
                        loadInfo.rho.segyInfo, ...
                        inIds, crossIds, startTimes, sampNum, GPreInvParam.dt);
                    
                    res.source = 'segy';
            end
            
        end
        
        if isempty(res.source)
            % obtain results by computing
            [vp, vs, rho] = bsCallInvFcn();
            res.source = 'computation';
        end
        
        res.data = {vp, vs, rho};
        res.inIds = inIds;          % inverted inline ids 
        res.crossIds = crossIds;    % inverted crossline ids
        res.horizon = horizonTimes; % the horizon of inverted traces
        res.name = method.name;     % the name of this method
        res.dt = GPreInvParam.dt;
        res.upNum = GPreInvParam.upNum;
        res.downNum = GPreInvParam.downNum;
        
        if isfield(method, 'type')  % the type of inverted volume, IP, Seismic, VP, VS, Rho
            res.type = method.type;
        else
            res.type = {'vp', 'vs', 'rho'};
        end
        
        if isfield(method, 'showFiltCoef')
            res.showFiltCoef = method.showFiltCoef;
        else
            res.showFiltCoef = 0;
        end
        
        invResults{i} = res;
        
        
        % save sgy file
        if isfield(method, 'isSaveSegy') && method.isSaveSegy && strcmp(res.source, 'computation')
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
                    GPreInvParam.preSeisData.segyFileName, ...
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
                save(matFileName, 'vp', 'vs', 'rho', 'horizonTimes', 'inIds', 'crossIds', 'GPreInvParam', 'method');
            catch
                save(matFileName, 'vp', 'vs', 'rho', 'horizonTimes', 'inIds', 'crossIds', 'GPreInvParam', 'method', '-v7.3');
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
        method.parampkgs = output.parampkgs;
        
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
                
                bsIncParforProgress(pbm, iTrace, 50);
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
end

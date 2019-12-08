function model = bsPostPrepareModel(GPostInvParam, inline, crossline, horizonTime, trueLog, model)
%% create model package which involves d, G, m, etc.
% Programmed by: Bin She (Email: bin.stepbystep@gmail.com)
% Programming dates: Nov 2019
% -------------------------------------------------------------------------
    
    
    
    sampNum = GPostInvParam.upNum + GPostInvParam.downNum; 
    startTime = horizonTime - GPostInvParam.dt * GPostInvParam.upNum;
    
    [postSeisData, GPostInvParam.postSeisData.segyInfo] = bsReadTracesByIds(...
        GPostInvParam.postSeisData.fileName, ...
        GPostInvParam.postSeisData.segyInfo, ...
        inline, ...
        crossline, ...
        startTime, ...
        sampNum,...
        GPostInvParam.dt);
    
    if isempty(model)
        model.orginal_G = bsPostGenGMatrix(GPostInvParam.wavelet, sampNum);
        model.G = model.orginal_G;
    end
    
    % start time of the inverted time interval
    model.t0 = round(startTime / GPostInvParam.dt) * GPostInvParam.dt;
    model.inId = inline;
    model.crossId = crossline;
    model.d = postSeisData(1 : end-1);
    model.origianl_d = model.d;
    
    if isfield(GPostInvParam, 'seismicFiltCoef') && ~isempty(GPostInvParam.seismicFiltCoef)
        model.d = bsButtLowPassFilter(model.d, GPostInvParam.seismicFiltCoef);
    end
    
    switch lower(GPostInvParam.initModel.mode)
        % the source of initial model
        
        case 'segy' % get initial model from segy file
            % start location of the inverted time interval
            [initLog, GPostInvParam.initModel.segyInfo] = bsReadTracesByIds(...
                GPostInvParam.initModel.ip.fileName, ...
                GPostInvParam.initModel.ip.segyInfo, ...
                inline, ...
                crossline, ...
                startTime, ...
                sampNum, ...
                GPostInvParam.dt);
            initLog = bsButtLowPassFilter(initLog, GPostInvParam.initModel.filtCoef);
            
        case 'filter_from_true_log' % get initial model by filtering the true model
            if isempty(trueLog)
                error('When initModel.mode is filter_from_true_log, true welllog data must be inputed.');
            end
            initLog = bsButtLowPassFilter(trueLog, GPostInvParam.initModel.filtCoef);
            
        case 'function' % get initial model by calling a function
            
            if isempty(GPostInvParam.initModel.fcn)
                error('When GPostInvParam.initModel.mode is function, the GPreInvParam.initModel.fcn could not be empty!\n');
            end
            initLog = GPostInvParam.initModel.fcn(GPostInvParam, inline, crossline, startTime);
            
        otherwise
            validatestring(GPostInvParam.initModel.mode, ['segy', 'filter_from_true_log', 'function']);
    end
    
    if exist('trueLog', 'var') && ~isempty(trueLog)
        model.trueLog = trueLog;
        model.trueX = log(model.trueLog);
    end
    
    model.dTrue = model.d;
    model.initLog = initLog;
    model.initX = log(model.initLog);
       
    % set boundary information
    switch lower(GPostInvParam.bound.mode)
        case 'off'
            model.Lb = [];
            model.Ub = [];
        case 'fixed'
            if length(GPostInvParam.bound.Lb) == 1
                model.Lb = ones(sampNum, 1) * log(GPostInvParam.bound.Lb);
                model.Ub = ones(sampNum, 1) * log(GPostInvParam.bound.Ub);
            else
                model.Lb = log(GPostInvParam.bound.Lb);
                model.Ub = log(GPostInvParam.bound.Ub);
            end
        case 'based_on_init'
            model.Lb = log(model.initLog - GPostInvParam.bound.offset_init);
            model.Ub = log(model.initLog + GPostInvParam.bound.offset_init);
        otherwise
            validatestring(GPostInvParam.bound.mode, ['off', 'fixed', 'based_on_init']);
    end
    
    % normalize
    if GPostInvParam.isNormal
        model.maxAbsD = norm(model.d);
        model.d = model.d / model.maxAbsD;
        model.G = model.orginal_G / model.maxAbsD;    % we have to use the original G to normalize
    end
end
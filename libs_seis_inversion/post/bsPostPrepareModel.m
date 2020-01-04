function model = bsPostPrepareModel(GPostInvParam, inline, crossline, horizonTime, trueLog, model)
%% create model package which involves d, G, m, etc.
% Programmed by: Bin She (Email: bin.stepbystep@gmail.com)
% Programming dates: Nov 2019
% -------------------------------------------------------------------------
    
    
    
    sampNum = GPostInvParam.upNum + GPostInvParam.downNum; 
    startTime = horizonTime - GPostInvParam.dt * GPostInvParam.upNum;
    
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
                error('When GInvParam.initModel.mode is function, the GInvParam.initModel.fcn could not be empty!\n');
            end
            initLog = GPostInvParam.initModel.fcn(inline, crossline, startTime);
            
        otherwise
            validatestring(GPostInvParam.initModel.mode, ['segy', 'filter_from_true_log', 'function']);
    end
    
    [model.d, model.orginal_G, model.initX] = bsPostBuild_d_G_m(GPostInvParam, inline, crossline, startTime, initLog, model);
    
    % start time of the inverted time interval
    model.t0 = round(startTime / GPostInvParam.dt) * GPostInvParam.dt;
    model.inId = inline;
    model.crossId = crossline;
    model.dTrue = model.d;
    model.initLog = initLog;
       
    if exist('trueLog', 'var') && ~isempty(trueLog)
        model.trueLog = trueLog;
        model.trueX = log(model.trueLog);
        
        if GPostInvParam.errorModel.isUse
            model.d = model.orginal_G * model.trueX;
        end
    end
    
    if GPostInvParam.errorModel.isUse && isempty(trueLog)
        errorData = bsReadTracesByIds(...
                GPostInvParam.errorModel.fileName, ...
                GPostInvParam.errorModel.segyInfo, ...
                inline, ...
                crossline, ...
                startTime, ...
                sampNum-1, ...
                GPostInvParam.dt);
        model.d = model.d - errorData;
    end
    
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
            model.Lb = log(model.initLog - GPostInvParam.bound.offset_init.ip);
            model.Ub = log(model.initLog + GPostInvParam.bound.offset_init.ip);
        otherwise
            validatestring(GPostInvParam.bound.mode, ['off', 'fixed', 'based_on_init']);
    end
    
    % normalize
    if GPostInvParam.isNormal
        model.maxAbsD = norm(model.d);
        model.d = model.d / model.maxAbsD;
        model.G = model.orginal_G / model.maxAbsD;    % we have to use the original G to normalize
    else
        model.G = model.orginal_G;
    end
end
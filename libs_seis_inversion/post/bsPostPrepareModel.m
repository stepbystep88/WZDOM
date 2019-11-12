function model = bsPostPrepareModel(GPostInvParam, inline, crossline, horizonTime, trueLog, model)
%% create model package which involves d, G, m, etc.
% Programmed by: Bin She (Email: bin.stepbystep@gmail.com)
% Programming dates: Nov 2019
% -------------------------------------------------------------------------
    
    
    
    sampNum = GPostInvParam.upNum + GPostInvParam.downNum; 
    pos = bsCalcT0Pos(GPostInvParam, GPostInvParam.postSeisData.segyInfo, horizonTime);
    
    [postSeisData, GPostInvParam.postSeisData.segyInfo] = bsReadTracesByIds(...
        GPostInvParam.postSeisData.segyFileName, ...
        GPostInvParam.postSeisData.segyInfo, ...
        inline, ...
        crossline, ...
        pos, ...
        sampNum);
    
    if isempty(model)
        D = bsGen1DDiffOperator(sampNum, 1, 1);
        W = bsWaveletMatrix(sampNum-1, GPostInvParam.wavelet);
        model.orgig_G = 0.5 * W * D;
        model.G = model.orgig_G;
    end
    
    % start time of the inverted time interval
    model.t0 = round(horizonTime / GPostInvParam.dt) * GPostInvParam.dt - GPostInvParam.upNum * GPostInvParam.dt;
    model.inId = inline;
    model.crossId = crossline;
    model.d = postSeisData(1 : end-1);
    
    if isfield(GPostInvParam, 'seismicFiltCoef') && ~isempty(GPostInvParam.seismicFiltCoef)
        [b, a] = butter(10, GPostInvParam.seismicFiltCoef);
        model.d = filtfilt(b, a, model.d);
    end
    
    switch GPostInvParam.initModel.mode
        % the source of initial model
        
        case 'segy' % get initial model from segy file
            % start location of the inverted time interval
            initPos = bsCalcT0Pos(GPostInvParam, GPostInvParam.initModel.segyInfo, horizonTime);
            [initLog, GPostInvParam.initModel.segyInfo] = bsReadTracesByIds(...
                GPostInvParam.initModel.segyFileName, ...
                GPostInvParam.initModel.segyInfo, ...
                inline, ...
                crossline, ...
                initPos, ...
                sampNum);
            initLog = bsButtLowPassFilter(initLog, GPostInvParam.initModel.filtCoef);
        case 'filter_from_true_log' % get initial model by filtering the true model
            initLog = bsButtLowPassFilter(trueLog, GPostInvParam.initModel.filtCoef);
        case 'directly' % get initial model from an initial data directly
            initLog = GPostInvParam.initModel.initLog;
    end
    
    if exist('trueLog', 'var') && ~isempty(trueLog)
        model.trueLog = trueLog;
        model.trueX = log(model.trueLog);
    end
    
    model.dTrue = model.d;
    model.initLog = initLog;
    model.initX = log(model.initLog);
    model.pos = pos;
       
    % set boundary information
    switch GPostInvParam.bound.mode
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
    end
    
    % normalize
    if GPostInvParam.isNormal
        model.maxAbsD = max(abs(model.d));
        model.d = model.d / model.maxAbsD;
        model.G = model.orgig_G / model.maxAbsD;    % we have to use the original G to normalize
    end
end
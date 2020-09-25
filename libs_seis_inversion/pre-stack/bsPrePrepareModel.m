function model = bsPrePrepareModel(GPreInvParam, inline, crossline, horizonTime, trueLog, model)
%% create model package for prestack inversion which involves d, G, m, etc.
% Programmed by: Bin She (Email: bin.stepbystep@gmail.com)
% Programming dates: Nov 2019
% -------------------------------------------------------------------------
    
    
    
    sampNum = GPreInvParam.upNum + GPreInvParam.downNum; 
    startTime = horizonTime - GPreInvParam.dt * GPreInvParam.upNum;
    
    % load model
    initModel = GPreInvParam.initModel;
    switch lower(initModel.mode)
        % the source of initial model
        % return initLog: 2D matrix, first-fourth columns are depth, vp,
        % vs, rho, respectively.
        case 'segy' % get initial model from segy file
            % start location of the inverted time interval
            initLog = bsReadMultiSegyFiles(...
                [initModel.vp, initModel.vs, initModel.rho], ...
                inline, crossline, startTime, sampNum, GPreInvParam.dt);
            
            if isfield(initModel, 'format') && strcmp(initModel.format, 'ipis')
                initLog(:, 1) = initLog(:, 1) ./  initLog(:, 3);
                initLog(:, 2) = initLog(:, 2) ./  initLog(:, 3);
            end
            
            depth = bsGetDepth(initLog(:, 1), GPreInvParam.dt);
            initLog = [depth, initLog];
            initLog = bsFiltWelllog(initLog, initModel.filtCoef);
            
        case 'filter_from_true_log' % get initial model by filtering the true model
            initLog = bsFiltWelllog(trueLog, initModel.filtCoef);
            
        case 'function' % get initial model by calling a function
            
            if isempty(initModel.fcn)
                error('When GPreInvParam.initModel.mode is function, the GPreInvParam.initModel.fcn could not be empty!\n');
            end
            initLog = initModel.fcn(GPreInvParam, inline, crossline, startTime);
            
        otherwise
            validatestring(GPreInvParam.initModel.mode, {'segy', 'filter_from_true_log', 'function'});
    end
    
    
% -------------------------------------------------------------------------
    if ~isempty(model) && isfield(model, 'orginal_G')
        G = model.orginal_G;
    else
        G = [];
    end
    
    % build d, G, m
    [model.d, model.G, model.initX, model.lsdCoef, model.angleData] ...
        = bsPreBuild_d_G_m(GPreInvParam, inline, crossline, startTime, initLog, G);
    GPreInvParam.lsdCoef = model.lsdCoef;
    
    % check data
    if isnan(model.initX)
        error('There is nan data in model.initX');
    end
    
    if isnan(model.G)
        error('There is nan data in model.G');
    end
    
    if isnan(model.d)
        error('There is nan data in model.d');
    end
% -------------------------------------------------------------------------            
    
    % start time of the inverted time interval
    model.t0 = round(startTime / GPreInvParam.dt) * GPreInvParam.dt;
    model.inline = inline;
    model.crossline = crossline;
    model.initLog = initLog;
    model.original_d = model.d;
    
    if exist('trueLog', 'var') && ~isempty(trueLog)
        model.trueLog = trueLog;
        
        [model.trueX, ~] = bsPreBuildModelParam(trueLog, GPreInvParam.mode, model.lsdCoef);
    end
    
    % set boundary information
    [model.Lb, model.Ub] = bsGetBound(GPreInvParam, initLog);
    
    if GPreInvParam.isScale && ~isempty(GPreInvParam.bestPostSeisData)
        postData = bsGetPostSeisData(GPreInvParam, inline, crossline, startTime, sampNum-1);
        c = bsComputeGain(postData, GPreInvParam.bestPostSeisData);
        model.G = model.G * c;
    end
        
    % normalize
    if GPreInvParam.isNormal
        model.scaleFactor = norm(model.d);
        model.orginal_G = model.G;
        model.d = model.d / model.scaleFactor;
        model.G = model.G / model.scaleFactor;    % we have to use the original G to normalize
    end
end



function [Lb, Ub] = bsGetBound(GPreInvParam, initLog)
    bound = GPreInvParam.bound;
    sampNum = size(initLog, 1);
    
    switch lower(bound.mode)
        case 'off'
            Lb = [];
            Ub = [];
        case 'fixed'
            if length(bound.vp.Lb) == 1
                OneVector = ones(sampNum, 1);
                lvp = OneVector * bound.vp.Lb;
                uvp = OneVector * bound.vp.Ub;
                lvs = OneVector * bound.vs.Lb;
                uvs = OneVector * bound.vs.Ub;
                lrho = OneVector * bound.rho.Lb;
                urho = OneVector * bound.rho.Ub;
                
                Lb = bsPreBuildModelParam(...
                        [initLog(:, 1), lvp, lvs, lrho], ...
                        GPreInvParam.mode, ...
                        GPreInvParam.lsdCoef);
                Ub = bsPreBuildModelParam(...
                        [initLog(:, 1), uvp, uvs, urho], ...
                        GPreInvParam.mode, ...
                        GPreInvParam.lsdCoef);
            else
                Lb = bsPreBuildModelParam(...
                        [initLog(:, 1), bound.vp.Lb, bound.vs.Lb, bound.rho.Lb], ...
                        GPreInvParam.mode, ...
                        GPreInvParam.lsdCoef);
                Ub = bsPreBuildModelParam(...
                        [initLog(:, 1), bound.vp.Ub, bound.vs.Ub, bound.rho.Ub], ...
                        GPreInvParam.mode, ...
                        GPreInvParam.lsdCoef);
            end
        case 'based_on_init'
            offset = repmat([0, bound.offset_init.vp, bound.offset_init.vs, bound.offset_init.rho], sampNum, 1);
            
            Lb = bsPreBuildModelParam(...
                    initLog - offset, ...
                    GPreInvParam.mode, ...
                    GPreInvParam.lsdCoef);
            Ub = bsPreBuildModelParam(...
                    initLog + offset, ...
                    GPreInvParam.mode, ...
                    GPreInvParam.lsdCoef);
        otherwise
            validatestring(GPreInvParam.bound.mode, ['off', 'fixed', 'based_on_init']);
    end
end
